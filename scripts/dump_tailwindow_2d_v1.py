"""dump_tailwindow_2d_v1: 验证 RR +600s 尾窗对 eta_E 的实际影响（2D 小网格实跑对比）。

对 hu101 / ht1_001 两代表井，固定 adopted（corrected）口径 + nz=120/ny=24 + CFL 自适应，
分别以 tt = stop（runner 冻结口径）与 tt = stop+600（RR 尾窗口径）跑 2D 环空求解器，
比较 eta_E / cement_occ / 混浆指数。回答：+600s 尾窗里泵入的替浆液对 η_E 是抬高、压低还是中性。

用法：PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/dump_tailwindow_2d_v1.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import cemdisp.data.loaders as L
from cemdisp.data.fluid_spec import FluidRole
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver

NZ, NY = 120, 24

CORRECTED_KW = dict(
    dispersion_dt_scale=1.0,
    enable_yield_gate=True,
    enable_regime_split=True,
    enable_local_i3=True,
    e_clip_max=0.90,
)
roles = {FluidRole.LEAD, FluidRole.TAIL, FluidRole.INTERMEDIATE}


def runner_stop(cr, fluids) -> float:
    by = {f.name: f for f in fluids}
    fs = sorted(cr.fronts, key=lambda f: f.time_s)
    last_cem = None
    for f in fs:
        fl = by.get(f.fluid_name)
        if fl is not None and fl.role in roles:
            last_cem = f.time_s
    if last_cem is not None:
        for f in fs:
            fl = by.get(f.fluid_name)
            if fl is None or fl.role not in roles:
                continue
            if f.time_s >= last_cem - 1e-9:
                return float(f.time_s)
    return float(cr.cement_end_time_s)


def main() -> None:
    wells = [("hu101", L.load_hu101_tailpipe), ("ht1_001", L.load_ht1_001_tailpipe)]
    for label, loader in wells:
        well, fluids, schedule, _ = loader()
        cr = CasingFlowSolver(enable_gravity=True).run(well, fluids, schedule)
        by = {f.name: f for f in fluids}
        fs = sorted(cr.fronts, key=lambda f: f.time_s)
        last_cem = None
        for f in fs:
            fl = by.get(f.fluid_name)
            if fl is not None and fl.role in roles:
                last_cem = f.time_s
        stop = float(cr.cement_end_time_s)
        tt_rr = min(
            sum(0.0 if s.rate_m3_min <= 0 else s.volume_m3 / s.rate_m3_min * 60.0 for s in schedule.steps) + 1200.0,
            stop + 600.0,
        )
        print(f"{label}: stop={stop:.0f}s RR tt={tt_rr:.0f}s", flush=True)
        for tag, tt in (("stop", stop), ("stop+600", stop + 600.0)):
            solver = AnnulusD2DGASolver(total_t=tt, nz=NZ, ny=NY, enable_cfl_adaptive=True, **CORRECTED_KW)
            inlet = build_coupled_annulus_inlet_provider(
                cr, CasingFlowSolver(enable_gravity=True), fluids, split_cement_phases=True)
            t0 = time.perf_counter()
            res = solver.run(well, fluids, inlet)
            fr = res.summary["最终结果"]
            print(f"  {tag}: tt={tt:.0f}s eta_E={float(fr['全井段最终有效顶替效率']):.4f} "
                  f"cement_occ={float(fr['最终水泥浆占据率']):.4f} mixing={float(fr['最终混浆指数']):.4f} "
                  f"({time.perf_counter() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
