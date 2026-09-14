"""hu103 eta≈0 二分定位（2026-09-02）。

runner 基线配置重跑后 hu1 eta_E≈3e-6 / hu103≈0.0013，同一 loader、同一 tt=stop、
RR corrected 口径正常。本脚本 nz=60 下对 hu103 逐一叠加 CORRECTED_KW 成员做二分，
定位"救活"开关，并对 baseline 附加 wall/cement 场诊断验证机制假设。
"""

from __future__ import annotations

import numpy as np

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


def _fmt_row(label, final, wall, cement):
    col_frozen = int(((wall > 0.5).all(axis=0)).sum())
    col_mean = cement.mean(axis=0)
    nz_cols = int((col_mean > 1e-3).sum())
    eta = float(final["全井段最终有效顶替效率"])
    occ = float(final["最终水泥浆占据率"])
    return (
        label,
        eta,
        occ,
        float(wall.mean()),
        col_frozen,
        wall.shape[1],
        nz_cols,
        int(len(col_mean)),
        float(col_mean[0]),
    )


def main():
    well, fluids, schedule, _ = L.load_hu103_tailpipe()
    cs = CasingFlowSolver(enable_gravity=True)
    cr = cs.run(well, fluids, schedule)
    stop = float(cr.cement_end_time_s or cr.pumping_end_time_s)
    print(f"hu103 stop={stop:.0f}s")
    prov = build_coupled_annulus_inlet_provider(cr, cs, fluids, split_cement_phases=True)

    configs = [("baseline", {})]
    for k, v in CORRECTED_KW.items():
        configs.append((f"+{k}={v}", {k: v}))
    configs.append(("corrected_full", dict(CORRECTED_KW)))

    print("label | eta_E | occ | wall_mean | 全1列 | nz | cement非零列 | col0")
    for label, kw in configs:
        solver = AnnulusD2DGASolver(total_t=stop, nz=60, enable_cfl_adaptive=True, **kw)
        res = solver.run(well, fluids, prov)
        row = _fmt_row(label, res.summary["最终结果"], res.wall_field, res.cement_field)
        print(row)
        import sys
        sys.stdout.flush()

    # baseline 诊断：死锁时间序列
    solver = AnnulusD2DGASolver(total_t=stop, nz=60, enable_cfl_adaptive=True)
    res = solver.run(well, fluids, prov)
    m = res.metrics
    mc = m["mean_cement"].to_numpy()
    mw = m["mean_wall_mud"].to_numpy()
    t = m["time_s"].to_numpy()
    idx = np.linspace(0, len(mc) - 1, 12).astype(int)
    print("--- baseline mean_cement/mean_wall 时间序列 ---")
    for i in idx:
        print(f"t={t[i]:.0f}s mean_cement={mc[i]:.5f} mean_wall={mw[i]:.4f}")
    growth = np.nonzero(np.diff(mc) > 1e-8)[0]
    if len(growth):
        print(f"mean_cement 最后增长 t={t[growth[-1] + 1]:.0f}s / {stop:.0f}s")


if __name__ == "__main__":
    main()
