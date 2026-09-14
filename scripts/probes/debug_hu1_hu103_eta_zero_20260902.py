"""hu1/hu103 η≈0 回归定位（2026-09-02）。

现象：runner 基线配置重跑后 hu1 eta_E≈3e-6、hu103≈0.0013，而同一 loader/同一 tt
在 RR corrected 口径下正常（hu1 0.558 / hu103 0.307）。
本脚本在 nz=60 快跑下做判别：
1) 打印 1D 时间轴（fronts、cement_end、pumping_end）；
2) 打印环空入口 provider 在整个窗口的相分数采样（判断水泥是否被喂进环空）；
3) baseline vs corrected 两配置各跑一次 2D，定位差异是否在 CORRECTED_KW。
schedule 参数已排除（仅喂 Tier0 诊断，annulus_d2dga.py:943-945）。
"""

from __future__ import annotations

import cemdisp.data.loaders as L
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver

WELLS = {
    "hu1": L.load_hu1_tailpipe,
    "hu103": L.load_hu103_tailpipe,
}

CORRECTED_KW = dict(
    dispersion_dt_scale=1.0,
    enable_yield_gate=True,
    enable_regime_split=True,
    enable_local_i3=True,
    e_clip_max=0.90,
)


def probe(name: str, loader) -> None:
    well, fluids, schedule, _ = loader()
    cs = CasingFlowSolver(enable_gravity=True)
    cr = cs.run(well, fluids, schedule)
    stop = float(cr.cement_end_time_s or cr.pumping_end_time_s)
    print(f"\n=== {name}: pumping_end={cr.pumping_end_time_s:.0f}s cement_end={cr.cement_end_time_s} -> stop={stop:.0f}s")
    for fr in sorted(cr.fronts, key=lambda x: x.time_s):
        print(f"  front {fr.fluid_name}: t={fr.time_s:.1f}s")
    prov = build_coupled_annulus_inlet_provider(cr, cs, fluids, split_cement_phases=True)
    print("  入口相分数采样：")
    for frac in (0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.97):
        ts = frac * stop
        st = prov(ts)
        parts = []
        for attr in ("flow_rate_m3_s", "stage_name", "phase_fractions"):
            parts.append(f"{attr}={getattr(st, attr, '?')}")
        print(f"    t={ts:8.0f}s  " + "  ".join(parts))
    for label, kw in (("baseline", {}), ("corrected", CORRECTED_KW)):
        solver = AnnulusD2DGASolver(total_t=stop, nz=60, enable_cfl_adaptive=True, **kw)
        res = solver.run(well, fluids, prov)
        final = res.summary["最终结果"]
        print(
            f"  [{label}] nz=60 eta_E={final['全井段最终有效顶替效率']:.6f} "
            f"occ={final['最终水泥浆占据率']:.6f} "
            f"mixing={final['最终混浆指数']:.4f} channeling={final['最终窜槽指数']:.4f}"
        )
    # 二分：逐个 corrected kwarg 加到 baseline 上，定位"救活"hu103 的开关
    print("  --- 二分 corrected kwargs ---")
    for key, val in CORRECTED_KW.items():
        solver = AnnulusD2DGASolver(total_t=stop, nz=60, enable_cfl_adaptive=True, **{key: val})
        res = solver.run(well, fluids, prov)
        final = res.summary["最终结果"]
        print(
            f"  [{key}={val}] nz=60 eta_E={final['全井段最终有效顶替效率']:.6f}"
        )


def main() -> None:
    for name, loader in WELLS.items():
        probe(name, loader)


if __name__ == "__main__":
    main()
