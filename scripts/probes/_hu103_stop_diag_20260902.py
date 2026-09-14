"""诊断 hu103（及对照井）套管1D停止时刻：为何 cement_end_time 时只有部分水泥跨过鞋口。"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cemdisp.data.fluid_spec import FluidRole
from cemdisp.data.loaders import load_hu103_tailpipe, load_hu102_tailpipe
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver

CEMENT = {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}

def diag(name, loader):
    ws, fluids, sched, _ = loader()
    role = {f.name: f.role for f in fluids}
    cs = CasingFlowSolver(enable_gravity=True)
    cr = cs.run(ws, fluids, sched)
    prov = build_coupled_annulus_inlet_provider(cr, cs, fluids, split_cement_phases=True)
    v_cem = sum(s.volume_m3 for s in sched.steps if role.get(s.fluid_name) in CEMENT)
    stop = cr.cement_end_time_s
    pend = cr.pumping_end_time_s
    # 逐秒积分跨过鞋口水泥，找达到 99% 设计水泥量的时刻
    dt = 1.0
    cum = 0.0; t_full = None
    n = int(np.ceil(pend, )) + 1
    for k in range(n+1):
        t = k*dt
        st = prov(t); fr = dict(st.phase_fractions)
        cum += st.flow_rate_m3_s*(fr.get("lead",0)+fr.get("tail",0)+fr.get("cement",0))*dt
        if t_full is None and cum >= 0.99*v_cem:
            t_full = t
    # stop 时刻已入环空比例
    cum_stop = 0.0
    for k in range(int(stop)+1):
        st = prov(k*dt); fr=dict(st.phase_fractions)
        cum_stop += st.flow_rate_m3_s*(fr.get("lead",0)+fr.get("tail",0)+fr.get("cement",0))*dt
    print(f"\n=== {name} ===")
    print(f"设计水泥量={v_cem:.2f} m3  当前cement_end={stop:.1f}s  泵注结束={pend:.1f}s")
    print(f"cement_end 时已入环空水泥={cum_stop:.2f} m3 = {cum_stop/v_cem:.1%}")
    print(f"99%水泥跨过鞋口的物理时刻 t_full={t_full}s  (cement_end 提前 {t_full-stop:.1f}s, 比值 {stop/t_full:.3f})")
    print("各水泥泵注步：")
    for s in sched.steps:
        if role.get(s.fluid_name) in CEMENT:
            print(f"  {s.fluid_name:>10} 起{s.start_time_s:8.1f} 止{s.end_time_s:8.1f} 体积{s.volume_m3:6.2f} 排量{s.rate_m3_min:.3f}")

diag("hu103", load_hu103_tailpipe)
diag("hu102", load_hu102_tailpipe)
