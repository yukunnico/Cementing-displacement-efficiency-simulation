"""对比 hu103 重力修正开关 / 弥散开关 下尾浆跨鞋口体积，定位尾浆被截原因。"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np
from cemdisp.data.loaders import load_hu103_tailpipe
from cemdisp.transport1d import CasingFlowSolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider

def integrate(grav, axial_disp=True):
    ws, fluids, sched, _ = load_hu103_tailpipe()
    cs = CasingFlowSolver(enable_gravity=grav, enable_axial_dispersion=axial_disp)
    res = cs.run(ws, fluids, sched)
    prov = build_coupled_annulus_inlet_provider(res, cs, fluids, split_cement_phases=True)
    tend = res.pumping_end_time_s; dt = 2.0; k = 0; acc = {}
    while k*dt <= tend:
        st = prov(k*dt)
        for key, v in st.phase_fractions:
            acc[key] = acc.get(key, 0) + st.flow_rate_m3_s*v*dt
        k += 1
    print(f"gravity={grav} axialDisp={axial_disp} cement_end={res.cement_end_time_s:.0f} "
          f"pump_end={tend:.0f} lead={acc.get('lead',0):.1f} tail={acc.get('tail',0):.1f} "
          f"spacer={acc.get('spacer',0):.1f} mud={acc.get('mud',0):.1f} 合计={sum(acc.values()):.1f}")

for g in (True, False):
    for d in (True, False):
        integrate(g, d)
