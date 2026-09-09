"""Task1 预检：hu103 界面 t_contact vs t_travel 比值（裁定2 风险验证）"""
import math
from cemdisp.data.loaders import load_hu103_tailpipe
from cemdisp.transport1d.casing_flow import CasingFlowSolver

well, fluids, schedule, _ = load_hu103_tailpipe()
solver = CasingFlowSolver(enable_gravity=True)
area = solver._pipe_cross_section_area(well)
pipe_vol = solver._timeline_pipe_volume(well, well.shoe_md_m * area)
cs_vol = well.shoe_md_m * area
steps = solver._build_scheduled_steps(schedule)
steps_cut, _ = solver._displacement_sequence_cutoff(steps)
initial_fluid = solver._initial_fluid_name(fluids, schedule)
r = (well.liner_id_mm or 100.0) / 2000.0

def inject_time(cum):
    for s in steps_cut:
        if cum <= s.cumulative_volume_end_m3 + 1e-12:
            if s.step.rate_m3_min <= 0.0:
                return s.end_time_s
            return s.start_time_s + max(cum - s.cumulative_volume_start_m3, 0.0) / s.step.rate_m3_min * 60.0
    return None

print(f"pipe_vol(timeline)={pipe_vol:.4f}  cs_vol={cs_vol:.4f}  ratio={pipe_vol/cs_vol:.6f}")
print(f"travel-equiv vol (shoe*pi*r^2, r=lower ID)={well.shoe_md_m*math.pi*r*r:.4f}  ratio={pipe_vol/(well.shoe_md_m*math.pi*r*r):.6f}")
print()
print(f"{'step':<14}{'t_arr':>10}{'t_inj':>9}{'t_contact':>11}{'Q':>8}{'t_travel':>10}{'ratio':>8}{'sig_rel%':>9}")
for i, s in enumerate(steps_cut):
    t = solver._front_arrival_time(s, steps_cut, pipe_vol)
    if t is None:
        continue
    if solver.enable_gravity:
        t = solver._gravity_corrected_arrival_time(
            t, s.step.fluid_name, solver._displaced_fluid_name(steps_cut, i, initial_fluid), fluids, well)
    t_inj = inject_time(s.cumulative_volume_start_m3)
    # Q at arrival: active step rate
    q = 0.0
    for s2 in steps_cut:
        if s2.start_time_s <= t < s2.end_time_s:
            q = s2.step.rate_m3_min / 60.0
            break
    U = q / (math.pi * r * r) if q > 0 else float('nan')
    t_travel = well.shoe_md_m / U
    tc = t - t_inj
    print(f"{s.step.fluid_name:<14}{t:>10.1f}{t_inj:>9.1f}{tc:>11.1f}{q:>8.4f}{t_travel:>10.1f}{tc/t_travel:>8.4f}{(math.sqrt(tc/t_travel)-1)*100:>8.1f}%")
