"""修复验证：hu1/hu103 baseline（无 kwargs）nz=60 复跑，确认死锁解除。"""
from __future__ import annotations

import numpy as np

import cemdisp.data.loaders as L
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver

WELLS = [
    ("hu1", L.load_hu1_tailpipe),
    ("hu103", L.load_hu103_tailpipe),
]

for name, loader in WELLS:
    well, fluids, schedule, _ = loader()
    cs = CasingFlowSolver(enable_gravity=True)
    cr = cs.run(well, fluids, schedule)
    stop = float(cr.cement_end_time_s or cr.pumping_end_time_s)
    prov = build_coupled_annulus_inlet_provider(cr, cs, fluids, split_cement_phases=True)
    solver = AnnulusD2DGASolver(total_t=stop, nz=60, enable_cfl_adaptive=True)
    res = solver.run(well, fluids, prov)
    final = res.summary["最终结果"]
    wall = res.wall_field
    col_frozen = int(((wall > 0.5).all(axis=0)).sum())
    col_mean = res.cement_field.mean(axis=0)
    nz_cols = int((col_mean > 1e-3).sum())
    print(
        f"[{name}] stop={stop:.0f}s eta_E={final['全井段最终有效顶替效率']:.6f} "
        f"occ={final['最终水泥浆占据率']:.6f} wall_mean={float(wall.mean()):.4f} "
        f"全1列={col_frozen}/{wall.shape[1]} cement非零列={nz_cols}/{len(col_mean)}"
    )
