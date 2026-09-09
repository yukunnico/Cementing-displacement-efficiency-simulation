# -*- coding: utf-8 -*-
"""
探针脚本（只读，不改任何求解器代码）：量化 1D 套管段混浆过渡带宽度现状。

目的（2026-09-08 套管内浓度剖面真解调研）：
1. 对 3 口代表井（hu101 / hu102 / ht1_003）按生产口径（CasingFlowSolver(enable_gravity=True)，
   其余默认：alpha=0.25、混浆增强开、has_plug=False）复算每个流体界面的：
   - Taylor-Aris 未截断 D（用 alpha=1e30 的求解器读 min 内部的渐近值）
   - 对流尺度上限 d_cap = 0.25·U·R 及截断比（taylor/d_cap）
   - 混浆增强因子（At>0 且 Re>100 判据）
   - 生效 σ_t（含 0.5·t_travel / dt 双向 clamp）
   - 过渡带等效管内长度 L_mix = 2σ_t·U 与体积 2σ_t·Q（m³ / bbl）
2. 汇总每井"全部过渡带体积 / 环空评价域体积"百分比（评价域体积与质量平衡取证同口径
   _physical_annular_volume）。
3. 对照生产弥散时间线中实际生成的多相过渡事件数（交叉验证）。

只读：不改 cemdisp 包内任何文件；结果打印 + 落 JSON 到本目录。
"""
from __future__ import annotations

import json
import math
from pathlib import Path

from cemdisp.data.fluid_spec import FluidRole
from cemdisp.data.loaders import (
    load_hu101_tailpipe,
    load_hu102_tailpipe,
    load_ht1_003_tailpipe,
)
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
from cemdisp.transport1d import CasingFlowSolver
from cemdisp.transport1d.shoe_timeline import ShoeEventKind

OUT_DIR = Path(__file__).resolve().parent

WELLS = {
    "hu101": load_hu101_tailpipe,
    "hu102": load_hu102_tailpipe,
    "ht1_003": load_ht1_003_tailpipe,
}

BBL_M3 = 0.1589873  # 1 bbl = 0.1589873 m³


def probe_well(name: str, loader) -> dict:
    well_spec, fluids, schedule, _ = loader()
    fluid_by_name = {f.name: f for f in fluids}

    # 生产口径（与 8 井 runner 一致）
    solver = CasingFlowSolver(enable_gravity=True)
    res = solver.run(well_spec, fluids, schedule)

    # 关弥散重跑：shoe_timeline 此时 = _apply_dispersion_to_timeline 的原始输入事件序列
    solver_sharp = CasingFlowSolver(enable_gravity=True, enable_axial_dispersion=False)
    res_sharp = solver_sharp.run(well_spec, fluids, schedule)
    sharp_events = res_sharp.shoe_timeline.events

    # 无上限口径：alpha 巨大 → _compute_dispersion_coefficient 返回 min(taylor, 1e30·U·R) = taylor
    solver_uncapped = CasingFlowSolver(enable_gravity=True, dispersion_alpha=1.0e30)

    # 与 casing_flow.py:577 完全一致的半径口径
    R = (well_spec.liner_id_mm or 100.0) / 2000.0
    dt = solver.dt
    alpha = solver.dispersion_alpha

    interfaces: list[dict] = []
    for i, event in enumerate(sharp_events):
        if event.kind != ShoeEventKind.FRONT_ARRIVAL or not event.phase_fractions:
            continue
        fluid_name = event.phase_fractions[0][0]
        fluid = fluid_by_name.get(fluid_name)
        if fluid is None:
            continue
        # 与 casing_flow.py:600-609 一致的"向前找异名流体"逻辑
        prev_fluid = ""
        for j in range(i - 1, -1, -1):
            if sharp_events[j].phase_fractions and sharp_events[j].phase_fractions[0][0] != fluid_name:
                prev_fluid = sharp_events[j].phase_fractions[0][0]
                break
        if not prev_fluid:
            continue
        Q = event.flow_rate_m3_s
        U = Q / (math.pi * R ** 2)
        if U < 1e-9:
            continue
        d_cap = alpha * U * R
        taylor_uncapped = solver_uncapped._compute_dispersion_coefficient(R, fluid, U)
        d_eff_capped = solver._compute_dispersion_coefficient(R, fluid, U)
        prev_spec = fluid_by_name[prev_fluid]
        inst = solver._interface_instability_factor(fluid, prev_spec, R, U)
        d_eff = d_eff_capped * inst
        t_travel = well_spec.shoe_md_m / U
        raw_sigma = math.sqrt(2.0 * d_eff * t_travel) / U
        sigma_hi = min(raw_sigma, 0.5 * t_travel)
        sigma_t = max(sigma_hi, dt)
        is_capped = taylor_uncapped > d_cap * 1.001  # 截断是否生效
        band_vol = 2.0 * sigma_t * Q  # 过渡窗内 Q 近似恒定
        interfaces.append({
            "界面": f"{prev_fluid}→{fluid_name}",
            "到达时刻_s": round(event.time_s, 1),
            "Q_m3_min": round(Q * 60.0, 3),
            "U_m_s": round(U, 3),
            "taylor_未截断_D_m2_s": float(f"{taylor_uncapped:.3e}"),
            "d_cap_m2_s": float(f"{d_cap:.3e}"),
            "截断比_taylor除以cap": round(taylor_uncapped / d_cap, 1) if d_cap > 0 else None,
            "上限生效": is_capped,
            "混浆增强因子": round(inst, 2),
            "D_eff_生效": float(f"{d_eff:.3e}"),
            "σ_t_原始_s": round(raw_sigma, 2),
            "σ_t_生效_s": round(sigma_t, 2),
            "σ_被dt抬高": sigma_hi < dt,
            "带长L_mix_m": round(2.0 * sigma_t * U, 1),
            "带体积_m3": round(band_vol, 3),
            "带体积_bbl": round(band_vol / BBL_M3, 2),
        })

    # 环空评价域体积（与 scripts/_mass_balance_diag_20260902.py 同口径）
    v_ann = AnnulusD2DGASolver()._physical_annular_volume(well_spec)
    # 管内容积（时间轴口径，双内径感知）
    legacy_area = CasingFlowSolver._pipe_cross_section_area(well_spec)
    pipe_vol = CasingFlowSolver._timeline_pipe_volume(well_spec, well_spec.shoe_md_m * legacy_area)
    # 设计水泥体积
    v_cem_design = sum(
        s.volume_m3 for s in schedule.steps
        if fluid_by_name.get(s.fluid_name) is not None
        and fluid_by_name[s.fluid_name].role in {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}
    )
    total_band_vol = sum(itf["带体积_m3"] for itf in interfaces)
    # 生产弥散时间线中多相过渡事件数（交叉验证）
    n_multiphase = sum(
        1 for ev in res.shoe_timeline.events if len(ev.phase_fractions) > 1
    )
    return {
        "井名": name,
        "鞋深_m": well_spec.shoe_md_m,
        "管内容积_m3": round(pipe_vol, 2),
        "环空评价域体积_m3": round(v_ann, 2),
        "设计水泥量_m3": round(v_cem_design, 2),
        "stop_cement_end_s": round(float(res.cement_end_time_s), 1),
        "界面数_有过渡带": len(interfaces),
        "弥散时间线多相事件数": n_multiphase,
        "全井过渡带体积合计_m3": round(total_band_vol, 3),
        "全井过渡带体积合计_bbl": round(total_band_vol / BBL_M3, 2),
        "过渡带体积占环空评价域百分比": round(100.0 * total_band_vol / v_ann, 3) if v_ann > 0 else None,
        "过渡带体积占设计水泥百分比": round(100.0 * total_band_vol / v_cem_design, 3) if v_cem_design > 0 else None,
        "各界面明细": interfaces,
    }


def main() -> None:
    out = {name: probe_well(name, loader) for name, loader in WELLS.items()}
    print(json.dumps(out, ensure_ascii=False, indent=2))
    (OUT_DIR / "probe_mixing_width_20260908_results.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    print("\n[done] JSON -> probe_mixing_width_20260908_results.json")


if __name__ == "__main__":
    main()
