"""变体验证：若 t_travel 改用时间线体积/去重力，2a/2b/2c 是否可通过"""
import math
from cemdisp.data.loaders import load_hu103_tailpipe
from cemdisp.transport1d.casing_flow import CasingFlowSolver

well, fluids, schedule, _ = load_hu103_tailpipe()
solver = CasingFlowSolver(enable_gravity=True)
area = solver._pipe_cross_section_area(well)
pipe_vol = solver._timeline_pipe_volume(well, well.shoe_md_m * area)
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

t_travel_tl = pipe_vol / 0.03  # 变体：t_travel = 时间线体积/Q（Q=0.03 恒定段）
print(f"variant t_travel (timeline vol/Q) = {t_travel_tl:.1f}")
print(f"{'step':<12}{'t_arr_grav':>11}{'t_arr_pist':>11}{'t_inj':>9}{'tc_grav':>9}{'tc_pist':>9}{'2a_grav':>9}{'2a_pist':>9}{'2c%_grav':>10}{'2c%_pist':>10}")
for i, s in enumerate(steps_cut):
    t_pist = solver._front_arrival_time(s, steps_cut, pipe_vol)
    if t_pist is None:
        continue
    t_grav = solver._gravity_corrected_arrival_time(
        t_pist, s.step.fluid_name, solver._displaced_fluid_name(steps_cut, i, initial_fluid), fluids, well)
    t_inj = inject_time(s.cumulative_volume_start_m3)
    tc_g = t_grav - t_inj
    tc_p = t_pist - t_inj
    ok_g = "PASS" if tc_g <= t_travel_tl + 1e-9 else "FAIL"
    ok_p = "PASS" if tc_p <= t_travel_tl + 1e-9 else "FAIL"
    d_g = (math.sqrt(tc_g / t_travel_tl) - 1) * 100
    d_p = (math.sqrt(tc_p / t_travel_tl) - 1) * 100
    print(f"{s.step.fluid_name:<12}{t_grav:>11.1f}{t_pist:>11.1f}{t_inj:>9.1f}{tc_g:>9.1f}{tc_p:>9.1f}{ok_g:>9}{ok_p:>9}{d_g:>9.2f}%{d_p:>9.2f}%")
