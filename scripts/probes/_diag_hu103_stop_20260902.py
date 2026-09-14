"""诊断 hu103（及对照 hu102）1D 停止时刻：为何 cement_end_time 时只有部分水泥跨过鞋口。"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from cemdisp.data.fluid_spec import FluidRole
from cemdisp.data.loaders import load_hu103_tailpipe, load_hu102_tailpipe
from cemdisp.transport1d import CasingFlowSolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider

def diag(name, loader):
    ws, fluids, sched, _ = loader()
    role = {f.name: f.role for f in fluids}
    solver = CasingFlowSolver(enable_gravity=True)
    res = solver.run(ws, fluids, sched)
    prov = build_coupled_annulus_inlet_provider(res, solver, fluids, split_cement_phases=True)
    v_cem = sum(s.volume_m3 for s in sched.steps if role.get(s.fluid_name) in
                {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL})
    print(f"\n===== {name} =====")
    print(f"泵注结束 pumping_end = {res.pumping_end_time_s:.1f} s ; cement_end = {res.cement_end_time_s:.1f} s")
    print("泵注步骤（流体[角色] 体积 累计 起止时间）:")
    cum = 0.0
    for s in sched.steps:
        cum += s.volume_m3
        r = role.get(s.fluid_name)
        tag = "  <== 水泥" if r in {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL} else ""
        print(f"  {s.fluid_name:<12}[{str(r).split('.')[-1]:<9}] V={s.volume_m3:7.2f} 累计={cum:7.2f} "
              f"rate={s.rate_m3_min if s.rate_m3_min is not None else float('nan'):6.3f}{tag}")
    # 逐秒积分跨鞋口水泥，找达到 99% 设计水泥量的时刻
    dt = 2.0
    tmax = res.pumping_end_time_s
    c = 0.0; t_star = None; c_at_end = 0.0
    k = 0
    while k*dt <= tmax:
        t = k*dt
        st = prov(t); fr = dict(st.phase_fractions)
        dcem = st.flow_rate_m3_s*(fr.get("lead",0)+fr.get("tail",0)+fr.get("cement",0))*dt
        c += dcem
        if t <= res.cement_end_time_s: c_at_end += dcem
        if t_star is None and c >= 0.99*v_cem: t_star = t
        k += 1
    print(f"设计水泥量={v_cem:.2f} m3; cement_end时累计跨鞋口={c_at_end:.2f} ({c_at_end/v_cem:.1%}); "
          f"达到99%设计量的真实时刻 t*={t_star}")
    print(f"cement_end 相对 t* 提前 {None if t_star is None else round(res.cement_end_time_s-t_star,1)} s")

diag("hu103", load_hu103_tailpipe)
diag("hu102", load_hu102_tailpipe)
