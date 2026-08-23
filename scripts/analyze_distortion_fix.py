"""环空失真修正前后对比分析（M0-M4+I3）。

对呼101三剖面（假设/扶正器间/扶正器处）跑两组：
  - baseline: 全部新开关关（逐位复现 2832a08 行为）
  - corrected: M1弥散dt归一 + M3屈服门槛 + M2流态修正 + I3局部化 + M4 e=0.90
对比 η_E(全井/评价窗)、η_N(窄四分位)、mixing、channeling、失稳、wall占比。

一次性结果分析脚本（results-phase），非生产代码。
"""
from __future__ import annotations

import csv
import time
from dataclasses import replace
from pathlib import Path
from typing import cast

import numpy as np

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.loaders import load_hu101_tailpipe
from cemdisp.data.pumping_schedule import PumpingSchedule
from cemdisp.data.well_spec import DepthValuePoint, WellSpec
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowResult, CasingFlowSolver

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PROJECT_ROOT / "results" / "失真修正前后对比"
TOP_MD, BOTTOM_MD = 5400.0, 7868.0
NZ = 80  # 快速对比（设计网格为140/250；80用于快速看方向）

ASSUMED = [(5400, .45), (6100, .38), (6796, .44), (7200, .48), (7600, .42), (7868, .46)]
BETWEEN = [(5400, .78), (5700, .72), (6000, .65), (6300, .58), (6600, .52),
           (6900, .48), (7200, .42), (7500, .32), (7700, .25), (7868, .22)]
ATCENT = [(5400, .88), (5700, .85), (6000, .80), (6300, .76), (6600, .72),
          (6900, .70), (7200, .68), (7500, .65), (7700, .62), (7868, .60)]


def _total_t(schedule: PumpingSchedule) -> float:
    return sum(0. if s.rate_m3_min <= 0 else s.volume_m3 / s.rate_m3_min * 60. for s in schedule.steps)


def _stop_t(cr: CasingFlowResult, fluids: tuple[FluidSpec, ...]) -> float:
    roles = {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}
    by = {f.name: f for f in fluids}
    fronts = sorted(cr.fronts, key=lambda f: f.time_s)
    last = next((f.time_s for f in fronts
                 if by.get(f.fluid_name) is not None and by[f.fluid_name].role in roles), None)
    if last is not None:
        for f in fronts:
            fl = by.get(f.fluid_name)
            if fl is None or fl.role in roles:
                continue
            if f.time_s >= last - 1e-9:
                return float(f.time_s)
    return float(cr.cement_end_time_s)


def _well(base: WellSpec, profile) -> WellSpec:
    return replace(base, standoff_profile=tuple(DepthValuePoint(md, v) for md, v in profile))


def run_case(label, well, fluids, schedule, *, corrected):
    print(f"\n{'='*64}\n  {label}  corrected={corrected}\n{'='*64}")
    t0 = time.perf_counter()
    cr = CasingFlowSolver(enable_gravity=True).run(well, fluids, schedule)
    inlet = build_coupled_annulus_inlet_provider(cr, CasingFlowSolver(enable_gravity=True),
                                                 fluids, split_cement_phases=True)
    total_t = min(_total_t(schedule) + 1200.0, _stop_t(cr, fluids) + 600.0)
    kw = dict(total_t=total_t, nz=NZ)
    if corrected:
        kw.update(dispersion_dt_scale=1.0, enable_yield_gate=True,
                  enable_regime_split=True, enable_local_i3=True, e_clip_max=0.90)
    solver = AnnulusD2DGASolver(**kw)
    res = solver.run(well, fluids, inlet)
    dt = time.perf_counter() - t0

    fr = cast(dict, res.summary["最终结果"])
    we = res.summary.get("评价窗效率", {})
    lt = res.summary.get("低尾指标", {})
    m = {
        "case": label, "corrected": corrected,
        "eta_E": float(fr["全井段最终有效顶替效率"]),
        "eta_N": float(fr.get("窄四分位效率", float("nan"))),
        "mixing": float(fr["最终混浆指数"]),
        "channeling": float(fr["最终窜槽指数"]),
        "instability": float(fr["最终失稳指数"]),
        "instab_lin": float(fr.get("最终失稳指数_线性", float("nan"))),
        "cement_occ": float(fr["最终水泥浆占据率"]),
        "wall_frac": float(np.mean(res.wall_field)) if res.wall_field is not None else float("nan"),
        "standoff_lt05_frac": float(lt.get("standoff低于0.5段占比", float("nan"))),
        "narrow_lowtail_frac": float(lt.get("窄边效率低于0.05域占比", float("nan"))),
        "elapsed_s": round(dt, 1),
    }
    for wname, wd in we.items():
        m[f"win::{wname}::eta_E"] = wd["eta_E"]
        m[f"win::{wname}::eta_N"] = wd["eta_N"]
    print(f"  η_E={m['eta_E']:.4f}  η_N={m['eta_N']:.4f}  mix={m['mixing']:.4f}  "
          f"chan={m['channeling']:.4f}  wall%={m['wall_frac']*100:.1f}  ({dt:.0f}s)")
    return m


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    base_well, fluids, schedule, _ = load_hu101_tailpipe()
    profiles = [("assumed", ASSUMED), ("between_cent", BETWEEN), ("at_cent", ATCENT)]
    rows = []
    for name, prof in profiles:
        well = _well(base_well, prof)
        rows.append(run_case(f"{name}_BASELINE", well, fluids, schedule, corrected=False))
        rows.append(run_case(f"{name}_CORRECTED", well, fluids, schedule, corrected=True))

    keys = sorted(set().union(*[r.keys() for r in rows]))
    csv_path = OUTPUT_DIR / "前后对比.csv"
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)

    print("\n" + "=" * 70 + "\n  关键指标：BASELINE → CORRECTED\n" + "=" * 70)
    for name, _ in profiles:
        b = next(r for r in rows if r["case"] == f"{name}_BASELINE")
        c = next(r for r in rows if r["case"] == f"{name}_CORRECTED")
        print(f"  {name:<16} η_E {b['eta_E']:.3f}→{c['eta_E']:.3f}   "
              f"η_N {b['eta_N']:.3f}→{c['eta_N']:.3f}   "
              f"mix {b['mixing']:.3f}→{c['mixing']:.3f}   "
              f"chan {b['channeling']:.3f}→{c['channeling']:.3f}")
    eN = {n: next(r for r in rows if r["case"] == f"{n}_CORRECTED")["eta_N"] for n, _ in profiles}
    mono = eN['at_cent'] >= eN['assumed'] >= eN['between_cent']
    print(f"\n  η_N 单调性（corrected，居中度越高η_N越高）: {mono}  "
          f"[between={eN['between_cent']:.4f} assumed={eN['assumed']:.4f} at={eN['at_cent']:.4f}]")
    print(f"  CSV: {csv_path}")


if __name__ == "__main__":
    main()
