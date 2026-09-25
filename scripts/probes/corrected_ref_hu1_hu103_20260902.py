"""corrected_full 配置参照值（hu1/hu103 nz=60）：用于验收"corrected 口径不得被修复改坏"。"""
from __future__ import annotations

import cemdisp.data.loaders as L
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver

CORRECTED_KW = dict(
    # 2026-09-26：`dispersion_dt_scale`（M1）形参已随 Task 7 删除自创弥散而移除，
    # 传了运行期即 TypeError ⇒ 键已删，corrected 口径只剩存活开关。
    enable_yield_gate=True,
    enable_regime_split=True,
    enable_local_i3=True,
    e_clip_max=0.90,
)

for name, loader in (("hu1", L.load_hu1_tailpipe), ("hu103", L.load_hu103_tailpipe)):
    well, fluids, schedule, _ = loader()
    cs = CasingFlowSolver(enable_gravity=True)
    cr = cs.run(well, fluids, schedule)
    stop = float(cr.cement_end_time_s or cr.pumping_end_time_s)
    prov = build_coupled_annulus_inlet_provider(cr, cs, fluids, split_cement_phases=True)
    solver = AnnulusD2DGASolver(total_t=stop, nz=60, enable_cfl_adaptive=True, **CORRECTED_KW)
    res = solver.run(well, fluids, prov)
    final = res.summary["最终结果"]
    print(
        f"[{name}] corrected_full nz=60 eta_E={final['全井段最终有效顶替效率']:.12f} "
        f"occ={final['最终水泥浆占据率']:.12f}"
    )
