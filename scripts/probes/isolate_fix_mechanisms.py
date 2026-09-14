"""隔离各修正机制，找出 corrected 配置 wall=100%/η_E=0.03 的元凶。"""
from __future__ import annotations
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
NZ = 80
ASSUMED = [(5400,.45),(6100,.38),(6796,.44),(7200,.48),(7600,.42),(7868,.46)]

def _total_t(s): return sum(0. if x.rate_m3_min<=0 else x.volume_m3/x.rate_m3_min*60. for x in s.steps)
def _stop_t(cr, fluids):
    roles={FluidRole.LEAD,FluidRole.INTERMEDIATE,FluidRole.TAIL}; by={f.name:f for f in fluids}
    fs=sorted(cr.fronts,key=lambda f:f.time_s)
    last=next((f.time_s for f in fs if by.get(f.fluid_name) and by[f.fluid_name].role in roles),None)
    if last is not None:
        for f in fs:
            fl=by.get(f.fluid_name)
            if fl is None or fl.role in roles: continue
            if f.time_s>=last-1e-9: return float(f.time_s)
    return float(cr.cement_end_time_s)

def run(label, well, fluids, schedule, **kw):
    t0=time.perf_counter()
    cr=CasingFlowSolver(enable_gravity=True).run(well,fluids,schedule)
    inlet=build_coupled_annulus_inlet_provider(cr,CasingFlowSolver(enable_gravity=True),fluids,split_cement_phases=True)
    tt=min(_total_t(schedule)+1200.,_stop_t(cr,fluids)+600.)
    res=AnnulusD2DGASolver(total_t=tt,nz=NZ,**kw).run(well,fluids,inlet)
    fr=cast(dict,res.summary["最终结果"])
    wf=float(np.mean(res.wall_field)) if res.wall_field is not None else float("nan")
    print(f"  {label:<22} η_E={fr['全井段最终有效顶替效率']:.4f} η_N={fr.get('窄四分位效率',float('nan')):.4f} "
          f"mix={fr['最终混浆指数']:.4f} chan={fr['最终窜槽指数']:.4f} wall%={wf*100:5.1f}  ({time.perf_counter()-t0:.0f}s)")

def main():
    base_well,fluids,schedule,_=load_hu101_tailpipe()
    well=replace(base_well,standoff_profile=tuple(DepthValuePoint(m,v) for m,v in ASSUMED))
    print(f"呼101 assumed 剖面, nz={NZ}\n"+"="*90)
    configs=[
        ("BASELINE(全关)", {}),
        ("M1 only", dict(dispersion_dt_scale=1.0)),
        ("M3 only", dict(enable_yield_gate=True)),
        ("M2 only", dict(enable_regime_split=True)),
        ("I3 only", dict(enable_local_i3=True)),
        ("M4 only(e=.90)", dict(e_clip_max=0.90)),
        ("M1+M3", dict(dispersion_dt_scale=1.0, enable_yield_gate=True)),
        ("M3+M4", dict(enable_yield_gate=True, e_clip_max=0.90)),
        ("M3+M2+M4", dict(enable_yield_gate=True, enable_regime_split=True, e_clip_max=0.90)),
        ("ALL corrected", dict(dispersion_dt_scale=1.0, enable_yield_gate=True,
                               enable_regime_split=True, enable_local_i3=True, e_clip_max=0.90)),
    ]
    for label,kw in configs:
        run(label,well,fluids,schedule,**kw)

if __name__=="__main__":
    main()
