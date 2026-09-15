"""
2026-09-02 结果偏低根因：质量平衡（库存）取证脚本（独立复核，不采信历史结论）。

阶段1（不跑2D，秒级）：对每口井核算
  - V_cem_design : 施工泵注中领/中/尾浆的设计总体积
  - V_ann        : 评价域物理环空体积（solver 同口径）
  - 库存比       = V_cem_design / V_ann
  - V_inj        : 按鞋口 provider 逐秒积分 [0,stop]，真正跨过鞋口进入环空的水泥体积
                   （含 erf 过渡带稀释，因此 V_inj <= V_cem_design）
  - 入库完成度   = V_inj / V_cem_design（检验 F4 入口稀释 / 停止时刻是否把水泥截在管内）

阶段2（代表井消融）：同一井在不同开关下跑 2D，比较最终域内水泥体积 V_dom、η_E、窄边、
  壁面冻结占比、前缘位置，量化每条“质量损耗/窜槽通道”的贡献。

只读取/计算，不改任何求解器代码；结果落 results/_质量平衡取证_2026-09-02/。
"""
from __future__ import annotations

import json
import time as _time
from pathlib import Path

import numpy as np

from cemdisp.data.fluid_spec import FluidRole
from cemdisp.data.loaders import (
    load_hu101_tailpipe, load_hu102_tailpipe, load_hu103_tailpipe,
    load_hu1_tailpipe, load_hu2_tailpipe, load_ht1_001_tailpipe,
    load_ht1_003_tailpipe, load_ht1_004_tailpipe,
)
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver, _trapez2d
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT = PROJECT_ROOT / "results" / "_质量平衡取证_2026-09-02"
OUT.mkdir(parents=True, exist_ok=True)

WELLS = {
    "hu101": load_hu101_tailpipe,
    "hu102": load_hu102_tailpipe,
    "hu103": load_hu103_tailpipe,
    "hu1": load_hu1_tailpipe,
    "hu2": load_hu2_tailpipe,
    "ht1_001": load_ht1_001_tailpipe,
    "ht1_003": load_ht1_003_tailpipe,
    "ht1_004": load_ht1_004_tailpipe,
}
CEMENT_ROLES = {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}


def build_case(loader):
    """复刻 runner：加载→套管1D→耦合 provider→停止时刻。"""
    well_spec, fluids, schedule, _ = loader()
    role_by_name = {f.name: f.role for f in fluids}
    casing_solver = CasingFlowSolver(enable_gravity=True)
    casing_result = casing_solver.run(well_spec, fluids, schedule)
    provider = build_coupled_annulus_inlet_provider(
        casing_result, casing_solver, fluids, split_cement_phases=True)
    stop = float(casing_result.cement_end_time_s)
    # 设计水泥体积（按泵注步骤流体角色汇总）
    v_cem_design = sum(s.volume_m3 for s in schedule.steps
                       if role_by_name.get(s.fluid_name) in CEMENT_ROLES)
    return well_spec, fluids, schedule, provider, stop, v_cem_design


def integrate_injection(provider, stop, dt=1.0):
    """逐秒积分 [0,stop] 跨过鞋口的各相体积。"""
    vols = {"cement": 0.0, "spacer": 0.0, "flusher": 0.0, "mud": 0.0}
    n = int(np.ceil(stop / dt))
    for k in range(n + 1):
        t = min(k * dt, stop)
        st = provider(t)
        q = st.flow_rate_m3_s
        frac = dict(st.phase_fractions)
        vols["cement"] += q * (frac.get("lead", 0.0) + frac.get("tail", 0.0)
                               + frac.get("cement", 0.0)) * dt
        vols["spacer"] += q * frac.get("spacer", 0.0) * dt
        vols["flusher"] += q * frac.get("flusher", 0.0) * dt
        vols["mud"] += q * frac.get("mud", 0.0) * dt
    return vols


def stage1():
    rows = {}
    for name, loader in WELLS.items():
        well_spec, fluids, schedule, provider, stop, v_cem_design = build_case(loader)
        solver_tmp = AnnulusD2DGASolver()
        v_ann = solver_tmp._physical_annular_volume(well_spec)
        inj = integrate_injection(provider, stop)
        rows[name] = {
            "停止时刻_s": round(stop, 1),
            "设计水泥量_m3": round(v_cem_design, 3),
            "物理环空体积_m3": round(v_ann, 3),
            "库存比_设计水泥除以环空": round(v_cem_design / v_ann, 4),
            "实际入环空水泥_m3": round(inj["cement"], 3),
            "入库完成度_入环空除以设计": round(inj["cement"] / max(v_cem_design, 1e-9), 4),
            "入环空前置液_m3": round(inj["spacer"], 3),
            "域顶m": well_spec.top_md_m, "域底m": well_spec.bottom_md_m,
        }
        print(name, json.dumps(rows[name], ensure_ascii=False))
    (OUT / "阶段1_库存核算.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    return rows


def run_variant(name, loader, nz, **solver_kw):
    """跑单个消融配置，返回质量平衡与最终场诊断。"""
    well_spec, fluids, schedule, provider, stop, v_cem_design = build_case(loader)
    v_ann = AnnulusD2DGASolver()._physical_annular_volume(well_spec)
    inj = integrate_injection(provider, stop)
    t0 = _time.time()
    solver = AnnulusD2DGASolver(total_t=stop, nz=nz, ny=40, **solver_kw)
    res = solver.run(well_spec, fluids, provider, schedule=schedule)
    geom = res.geom
    cement = res.cement_field
    half_v = _trapez2d(geom["b"], geom)
    v_dom = 2.0 * _trapez2d(geom["b"] * cement, geom)   # 域内最终水泥体积
    s_max = float(geom["s"][-1])
    fin = res.metrics.iloc[-1]
    out = {
        "配置": name,
        "nz": nz,
        "eta_E": round(float(fin["effective_efficiency"]), 4),
        "域内水泥_m3": round(v_dom, 3),
        "入环空水泥_m3": round(inj["cement"], 3),
        "守恒率_域内除以入环空": round(v_dom / max(inj["cement"], 1e-9), 4),
        "损耗水泥_m3": round(inj["cement"] - v_dom, 3),
        "窄边etaN": round(float(res.summary["eta_narrow"]), 4),
        "宽边均水泥浓度": round(float(np.mean(cement[0])), 4),
        "中线均水泥浓度": round(float(np.mean(cement[cement.shape[0] // 2])), 4),
        "壁面冻结占比": round(float(fin["mean_wall_mud"]), 4),
        "宽边前缘m": round(float(fin["front_wide_m"]), 1),
        "窄边前缘m": round(float(fin["front_narrow_m"]), 1),
        "域长m": round(s_max, 1),
        "宽边到顶": bool(fin["front_wide_m"] >= s_max - 1.0),
        "耗时_s": round(_time.time() - t0, 1),
        "物理环空体积_m3": round(v_ann, 3),
        "库存比": round(v_cem_design / v_ann, 3),
    }
    print(f"  [{name}] " + json.dumps(out, ensure_ascii=False))
    return out


def stage2(well_names, nz):
    all_res = {}
    # 消融配置：逐项关闭可疑通道
    # ⚠️ 2026-09-15 Task 13 死参清理：e_clip_max 传参已删（Task 10 弃用，e_clip 硬截断
    # 移除后传值不生效），dispersion_* 同理（Task 7 弃用）——"近同心e0.05"与
    # "生产口径e0.90"两变体因此与 BASE基线 等效（仅保留键名供历史 JSON 结构对比）。
    # c_min 形参已于 2026-09-07（B2）从 solver 删除："无壁面冻结"改用活参数
    # enable_yield_gate=False 表达同一意图（wall 恒零，见 solver __init__ docstring）。
    variants = {
        "BASE基线": dict(),
        "无壁面冻结": dict(enable_yield_gate=False),
        "近同心e0.05": dict(),
        "无弥散": dict(),
        "无I3浮力通量": dict(enable_d2dga_i3_flux=False),
        "无D2DGA放大": dict(enable_d2dga=False, enable_true_buoyancy=False),
        "理想活塞_同心无壁无弥散无I3无放大": dict(
            enable_yield_gate=False,
            enable_d2dga_i3_flux=False,
            enable_d2dga=False, enable_true_buoyancy=False),
        "生产口径e0.90": dict(),
    }
    for wname in well_names:
        loader = WELLS[wname]
        print(f"\n===== {wname} (nz={nz}) =====")
        all_res[wname] = {}
        for vname, kw in variants.items():
            try:
                all_res[wname][vname] = run_variant(vname, loader, nz, **kw)
            except Exception as exc:  # 单个配置失败不中断
                all_res[wname][vname] = {"配置": vname, "error": f"{type(exc).__name__}: {exc}"}
                print(f"  [{vname}] ERROR {exc}")
    (OUT / f"阶段2_通道消融_nz{nz}.json").write_text(
        json.dumps(all_res, ensure_ascii=False, indent=2), encoding="utf-8")
    return all_res


if __name__ == "__main__":
    import sys
    mode = sys.argv[1] if len(sys.argv) > 1 else "stage1"
    if mode == "stage1":
        stage1()
    elif mode == "stage2":
        nz = int(sys.argv[3]) if len(sys.argv) > 3 else 120
        stage2(sys.argv[2].split(","), nz)
