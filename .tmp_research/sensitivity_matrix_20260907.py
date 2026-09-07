"""敏感性分析实验（2026-09-07）：现场可控杠杆 → η_E/η_N 响应矩阵。

设计原则：只扫【现场可操作/可设计的参数】，回答"模型敏感性分析能否指导现场"。
杠杆与现场对应：
- 排量 Q（±20%/±40%）：泵注程序排量（现场最直接可控）
- 隔离液密度（±0.10 g/cc）：隔离液配方设计
- 水泥浆 n（±0.1）：水泥浆配方设计（分散剂剂量）
- 泥浆 PV（±30%）：钻井液流变调整（固井前调整）
- 居中度 standoff（±0.1）：扶正器方案（hu103 已有 09-06 数据，本脚本补 hu101/ht1_003）

井选择：hu103（低偏心、最快）+ hu101（强偏心、判别力最强）全杠杆；ht1_003 补 standoff。
每变体 nz=250 ny=40 生产口径（三项修复默认开）；排量/流体变化传导到 1D 时序（1D 重跑）。
输出：results/敏感性分析_2026-09-07/敏感性矩阵.csv
"""
from __future__ import annotations

import dataclasses
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib  # noqa: E402

matplotlib.use("Agg")

from cemdisp.data.fluid_spec import FluidRole, FluidSpec  # noqa: E402
from cemdisp.data.well_spec import DepthValuePoint  # noqa: E402
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver  # noqa: E402
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider  # noqa: E402
from cemdisp.transport1d import CasingFlowSolver  # noqa: E402
from scripts._mass_balance_diag_20260902 import WELLS, build_case  # noqa: E402

OUT = PROJECT_ROOT / "results" / "敏感性分析_2026-09-07"
OUT.mkdir(parents=True, exist_ok=True)
NZ = 250


def _scale_rate(schedule, factor: float):
    """返回排量缩放后的 schedule 副本（体积不变、时长按比例变化）。"""
    steps = [dataclasses.replace(s, rate_m3_min=s.rate_m3_min * factor) if s.rate_m3_min > 0 else s
             for s in schedule.steps]
    return dataclasses.replace(schedule, steps=tuple(steps))


def _scale_density(fluids, role, delta_kg_m3: float):
    """返回指定角色流体密度平移后的 fluids 副本。"""
    return tuple(
        FluidSpec(f.name, f.role, f.density_kg_m3 + delta_kg_m3, f.rheology_model,
                  plastic_viscosity_pa_s=f.plastic_viscosity_pa_s,
                  yield_stress_pa=f.yield_stress_pa,
                  power_law_n=f.power_law_n, consistency_k=f.consistency_k)
        if f.role == role else f for f in fluids
    )


def _scale_power_law_n(fluids, roles: set, delta: float):
    """幂律 n 平移（只对指定角色的幂律流体，n 平移后 clip 到 [0.3, 1.2]）。"""
    out = []
    for f in fluids:
        if f.role in roles and f.rheology_model.value == "power_law" and f.power_law_n is not None:
            n_new = float(np.clip(f.power_law_n + delta, 0.3, 1.2))
            out.append(FluidSpec(f.name, f.role, f.density_kg_m3, f.rheology_model,
                                 plastic_viscosity_pa_s=f.plastic_viscosity_pa_s,
                                 yield_stress_pa=f.yield_stress_pa,
                                 power_law_n=n_new, consistency_k=f.consistency_k))
        else:
            out.append(f)
    return tuple(out)


def _scale_bingham_pv(fluids, roles: set, factor: float):
    """宾汉 PV 缩放（指定角色）。"""
    out = []
    for f in fluids:
        if f.role in roles and f.rheology_model.value == "bingham" and f.plastic_viscosity_pa_s:
            out.append(FluidSpec(f.name, f.role, f.density_kg_m3, f.rheology_model,
                                 plastic_viscosity_pa_s=f.plastic_viscosity_pa_s * factor,
                                 yield_stress_pa=f.yield_stress_pa,
                                 power_law_n=f.power_law_n, consistency_k=f.consistency_k))
        else:
            out.append(f)
    return tuple(out)


def run_variant(tag: str, well_name: str, *, schedule=None, fluids=None, well_override=None,
                solver_kw=None) -> dict:
    """跑单变体：1D 重跑（schedule/fluid 变化会影响 1D 时序）→ 2D 求解。"""
    t0 = time.time()
    well_spec, fluids0, schedule0, _provider, _stop, _v_design = build_case(WELLS[well_name])
    sched = schedule if schedule is not None else schedule0
    flu = fluids if fluids is not None else fluids0
    spec = well_override if well_override is not None else well_spec

    casing_solver = CasingFlowSolver(enable_gravity=True)
    casing_result = casing_solver.run(spec, flu, sched)
    provider = build_coupled_annulus_inlet_provider(
        casing_result, casing_solver, flu, split_cement_phases=True)
    stop = float(casing_result.cement_end_time_s)

    solver = AnnulusD2DGASolver(total_t=stop, nz=NZ, ny=40, **(solver_kw or {}))
    res = solver.run(spec, flu, provider, schedule=sched)
    fin = res.metrics.iloc[-1]
    g = res.geom
    s_max = float(g["s"][-1])
    row = {
        "井": well_name, "变体": tag,
        "η_E": round(float(fin["effective_efficiency"]), 4),
        "η_N": round(float(res.summary.get("eta_narrow", np.nan)), 4),
        "stop_s": round(stop, 1),
        "窄边前缘m": round(float(fin["front_narrow_m"]), 1),
        "窄边到位率": round(float(fin["front_narrow_m"]) / s_max, 4),
        "窜槽指数": round(float(fin["channeling_index"]), 4),
        "混浆指数": round(float(fin["mixing_index"]), 4),
        "耗时s": round(time.time() - t0, 1),
    }
    print(json.dumps(row, ensure_ascii=False), flush=True)
    return row


def sweep(well_name: str, base_row: dict | None = None) -> list[dict]:
    """单井全杠杆扫描。base_row 传入时跳过重跑基线。"""
    rows = [] if base_row is None else [base_row]
    well_spec, fluids, schedule, _provider, _stop, _v_design = build_case(WELLS[well_name])

    # --- 排量 ±20%/±40% ---
    for factor, tag in [(0.6, "排量-40%"), (0.8, "排量-20%"), (1.2, "排量+20%"), (1.4, "排量+40%")]:
        rows.append(run_variant(tag, well_name, schedule=_scale_rate(schedule, factor)))

    # --- 隔离液密度 ±0.10 g/cc ---
    for delta, tag in [(-100.0, "隔离液密度-0.10"), (+100.0, "隔离液密度+0.10")]:
        rows.append(run_variant(tag, well_name, fluids=_scale_density(fluids, FluidRole.SPACER, delta)))

    # --- 水泥浆 n ±0.1 ---
    cement_roles = {FluidRole.LEAD, FluidRole.TAIL, FluidRole.INTERMEDIATE}
    for delta, tag in [(-0.1, "水泥n-0.1"), (+0.1, "水泥n+0.1")]:
        rows.append(run_variant(tag, well_name, fluids=_scale_power_law_n(fluids, cement_roles, delta)))

    # --- 泥浆 PV ±30% ---
    for factor, tag in [(0.7, "泥浆PV-30%"), (1.3, "泥浆PV+30%")]:
        rows.append(run_variant(tag, well_name, fluids=_scale_bingham_pv(fluids, {FluidRole.MUD}, factor)))

    return rows


def standoff_sweep(well_name: str) -> list[dict]:
    """standoff ±0.1 扫描。"""
    rows = []
    well_spec, fluids, schedule, _provider, _stop, _v_design = build_case(WELLS[well_name])
    for delta, tag in [(-0.1, "standoff-0.1"), (+0.1, "standoff+0.1")]:
        pts = tuple(DepthValuePoint(depth_md_m=p.depth_md_m,
                                    value=float(np.clip(p.value + delta, 0.2, 1.0)))
                    for p in well_spec.standoff_profile)
        rows.append(run_variant(tag, well_name, well_override=dataclasses.replace(well_spec, standoff_profile=pts)))
    return rows


def main():
    all_rows: list[dict] = []
    for well in ("hu103", "hu101"):
        print(f"\n########## {well} ##########", flush=True)
        try:
            base = run_variant("基线", well)
            all_rows.extend(sweep(well, base_row=base))
        except Exception as exc:
            traceback.print_exc()
            all_rows.append({"井": well, "变体": "ERROR", "错误": f"{type(exc).__name__}: {exc}"})
        pd.DataFrame(all_rows).to_csv(OUT / "敏感性矩阵.csv", index=False, encoding="utf-8-sig")
    # ht1_003 补 standoff（hu103±0.1 已有 09-06 数据）
    try:
        all_rows.extend(standoff_sweep("ht1_003"))
    except Exception:
        traceback.print_exc()
    df = pd.DataFrame(all_rows)
    df.to_csv(OUT / "敏感性矩阵.csv", index=False, encoding="utf-8-sig")
    (OUT / "README.md").write_text(
        "# 敏感性分析矩阵（2026-09-07）\n\n"
        "- 生产口径 nz250 ny40（三项修复默认开）；每变体 1D 重跑（排量/流体变化传导到时序与停止时刻）\n"
        "- 杠杆：排量±20%/±40%、隔离液密度±0.10g/cc、水泥浆n±0.1、泥浆PV±30%、standoff±0.1\n"
        "- 井：hu103（低偏心）+hu101（强偏心）全杠杆；ht1_003 补 standoff；hu103 standoff±0.1 见 results/_增量取证_2026-09-06/\n"
        "- 用途：论文敏感性章节 + '能否指导现场'裁定\n",
        encoding="utf-8")
    print("\n[done] 敏感性矩阵已落盘", OUT)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
