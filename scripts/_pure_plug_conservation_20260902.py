"""
2026-09-02 决定性守恒实验：合成“纯水泥单相塞流”入口，剥离多相/过量/偏心干扰。

构造恒定排量、恒为 tail=1 的入口（无 spacer/flusher/mud 从入口进入），
在给定目标填充率 target = 注入体积/环空体积 下停止，检验：
  守恒率 = 域内水泥体积 / 注入体积  —— 恒定速度场半拉格朗日平流理应≈1。

分别扫描：
  - target 0.7 / 1.0 / 1.3（1.3 时多余水泥本应流出顶部，但顶部为反射边界）
  - 同心 e0.05 vs 偏心 e0.55
  - 无弥散 vs 默认弥散
  - 是否存在前置液占位（对照：先注入一段 spacer 再注水泥）
若“纯水泥单相、target<1、同心、无弥散”仍显著<1，则核心平流/注入本身不守恒。
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from cemdisp.data.loaders import load_hu102_tailpipe
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver, _trapez2d
from cemdisp.models2d.boundary_bridge import AnnulusInletState

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT = PROJECT_ROOT / "results" / "_质量平衡取证_2026-09-02"
OUT.mkdir(parents=True, exist_ok=True)


def make_pure_cement_provider(q_m3s: float):
    """恒定排量、纯 tail 水泥单相入口。"""
    def provider(t):
        return AnnulusInletState(time_s=t, flow_rate_m3_s=q_m3s,
                                 stage_name="纯水泥", phase_fractions=(("tail", 1.0),))
    return provider


def run_one(tag, target, e_clip, dispersion, nz, q_m3s=0.008, spacer_first_m3=0.0):
    well_spec, fluids, _, _ = load_hu102_tailpipe()
    v_ann = AnnulusD2DGASolver()._physical_annular_volume(well_spec)
    # 若先注前置液：前 t_sp 注 spacer，之后注水泥；总模拟到水泥注入 target*V_ann
    t_sp = spacer_first_m3 / q_m3s if spacer_first_m3 > 0 else 0.0
    v_inj_cement = target * v_ann
    t_cement = v_inj_cement / q_m3s
    total_t = t_sp + t_cement

    def provider(t):
        if t < t_sp:
            return AnnulusInletState(t, q_m3s, "前置液", (("spacer", 1.0),))
        return AnnulusInletState(t, q_m3s, "纯水泥", (("tail", 1.0),))

    kw = dict(e_clip_max=e_clip, c_min=0.0, enable_d2dga_i3_flux=False,
              enable_d2dga=False, enable_true_buoyancy=False,
              enable_cfl_adaptive=False, dt=2.0)
    if not dispersion:
        kw.update(dispersion_axial=0.0, dispersion_azimuthal=0.0)
    solver = AnnulusD2DGASolver(total_t=total_t, nz=nz, ny=40, **kw)
    res = solver.run(well_spec, fluids, provider)
    g = res.geom
    cement = res.cement_field
    v_dom = 2.0 * _trapez2d(g["b"] * cement, g)
    spacer = res.spacer_field
    v_sp = 2.0 * _trapez2d(g["b"] * np.clip(spacer, 0, 1), g)
    fin = res.metrics.iloc[-1]
    out = {
        "实验": tag, "target填充率": target, "e上限": e_clip, "弥散": dispersion,
        "前置液_m3": spacer_first_m3,
        "环空体积": round(v_ann, 3), "注入水泥": round(v_inj_cement, 3),
        "域内水泥": round(v_dom, 3), "域内前置液": round(v_sp, 3),
        "水泥守恒率": round(v_dom / v_inj_cement, 4),
        "eta_E": round(float(fin["effective_efficiency"]), 4),
        "宽边前缘": round(float(fin["front_wide_m"]), 1),
        "窄边前缘": round(float(fin["front_narrow_m"]), 1),
        "域长": round(float(g["s"][-1]), 1),
    }
    print(json.dumps(out, ensure_ascii=False))
    return out


if __name__ == "__main__":
    nz = 120
    results = []
    # 1) 纯水泥单相，target<1，同心，无弥散 —— 必须≈1 才算核心守恒
    results.append(run_one("纯水泥_0.7_同心_无弥散", 0.7, 0.05, False, nz))
    # 2) 同上，target=1.0
    results.append(run_one("纯水泥_1.0_同心_无弥散", 1.0, 0.05, False, nz))
    # 3) target=1.3（过量，顶部反射）
    results.append(run_one("纯水泥_1.3_同心_无弥散", 1.3, 0.05, False, nz))
    # 4) 纯水泥 0.7 偏心无弥散
    results.append(run_one("纯水泥_0.7_偏心0.55_无弥散", 0.7, 0.55, False, nz))
    # 5) 纯水泥 0.7 同心 开弥散
    results.append(run_one("纯水泥_0.7_同心_有弥散", 0.7, 0.05, True, nz))
    # 6) 先注 0.5V_ann 前置液，再注 0.7 水泥，同心无弥散（界面归一化对照）
    well_spec, _, _, _ = load_hu102_tailpipe()
    v_ann = AnnulusD2DGASolver()._physical_annular_volume(well_spec)
    results.append(run_one("前置液+水泥_0.7_同心_无弥散", 0.7, 0.05, False, nz,
                           spacer_first_m3=0.5 * v_ann))
    (OUT / "阶段3_纯水泥守恒实验.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
