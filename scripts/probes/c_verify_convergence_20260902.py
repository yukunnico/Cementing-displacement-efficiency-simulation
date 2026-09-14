"""C 修复验收 2-4：corrected 逐位回归 + baseline nz 收敛（hu1/hu103 × nz=60/150/250）。"""
from __future__ import annotations

import cemdisp.data.loaders as L
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver

CORRECTED_KW = dict(
    dispersion_dt_scale=1.0,
    enable_yield_gate=True,
    enable_regime_split=True,
    enable_local_i3=True,
    e_clip_max=0.90,
)

REF = {"hu1": 0.574314858968, "hu103": 0.314363677407}

for name, loader in (("hu1", L.load_hu1_tailpipe), ("hu103", L.load_hu103_tailpipe)):
    well, fluids, schedule, _ = loader()
    cs = CasingFlowSolver(enable_gravity=True)
    cr = cs.run(well, fluids, schedule)
    stop = float(cr.cement_end_time_s or cr.pumping_end_time_s)
    prov = build_coupled_annulus_inlet_provider(cr, cs, fluids, split_cement_phases=True)

    solver = AnnulusD2DGASolver(total_t=stop, nz=60, enable_cfl_adaptive=True, **CORRECTED_KW)
    res = solver.run(well, fluids, prov)
    eta_c = float(res.summary["最终结果"]["全井段最终有效顶替效率"])
    match = abs(eta_c - REF[name]) < 5e-13
    print(f"[{name}] corrected nz=60 eta_E={eta_c:.12f} ref={REF[name]:.12f} match={match}")

    etas = {}
    for nz in (60, 150, 250):
        solver = AnnulusD2DGASolver(total_t=stop, nz=nz, enable_cfl_adaptive=True)
        res = solver.run(well, fluids, prov)
        eta = float(res.summary["最终结果"]["全井段最终有效顶替效率"])
        wall = res.wall_field
        col_frozen = int(((wall > 0.5).all(axis=0)).sum())
        etas[nz] = eta
        print(f"[{name}] baseline nz={nz} eta_E={eta:.6f} wall_mean={float(wall.mean()):.4f} 全1列={col_frozen}/{nz}")

    spread = (max(etas.values()) - min(etas.values())) * 100.0
    print(f"[{name}] nz极差={spread:.2f}pp (nz60={etas[60]:.4f} nz150={etas[150]:.4f} nz250={etas[250]:.4f})")
