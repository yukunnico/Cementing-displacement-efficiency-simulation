import sys; sys.path.insert(0, 'tests')
from test_casing_mixing_contact_time import _synthetic_well, _synthetic_fluids
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep, PumpingStageEvent
from cemdisp.transport1d.casing_flow import CasingFlowSolver
from cemdisp.transport1d.shoe_timeline import ShoeEventKind
import math
shoe, area = 100.0, 0.01
well = _synthetic_well(shoe, area)
fluids = _synthetic_fluids(spacer_rho=None)
sched = PumpingSchedule(steps=(
  PumpingScheduleStep(step_name='注入隔离液', fluid_name='隔离液', volume_m3=0.4, rate_m3_min=1.0, start_time_s=0.0, end_time_s=24.0, event_tag=PumpingStageEvent.INJECT_SPACER),
  PumpingScheduleStep(step_name='注入尾浆', fluid_name='尾浆', volume_m3=1.0, rate_m3_min=1.0, start_time_s=24.0, end_time_s=84.0, event_tag=PumpingStageEvent.INJECT_CEMENT),
  PumpingScheduleStep(step_name='替浆', fluid_name='顶替液', volume_m3=2.0, rate_m3_min=1.0, start_time_s=84.0, end_time_s=204.0, event_tag=PumpingStageEvent.INJECT_DISPLACEMENT)))
s_on = CasingFlowSolver(enable_gravity=False, mixing_contact_time=True)
# 复刻 _build_shoe_timeline 内部，查看事件生成的 t_arrival 事件点
well_spec = well
pipe_area = s_on._pipe_cross_section_area(well_spec)
legacy_pv = well_spec.shoe_md_m * pipe_area
pv = s_on._timeline_pipe_volume(well_spec, legacy_pv)
steps_full = s_on._build_scheduled_steps(sched)
steps, _ = s_on._displacement_sequence_cutoff(steps_full)
initial_fluid = s_on._initial_fluid_name(fluids, sched)
for i, sch in enumerate(steps):
    ft = s_on._front_arrival_time(sch, steps, pv)
    print(f'step {i} {sch.step.fluid_name}: cum_start={sch.cumulative_volume_start_m3} front_arr={ft}')
# 关键：隔离液步（index 0）前缘到达 60s；尾浆步（index 1）前缘到达 84s
# 但事件循环里 FRONT_ARRIVAL 事件 phase=(尾浆,1.0) 只在 84s 生成一个
# 弥散处理时按流体消耗：尾浆列表 [84]，cursor 0 → t_arrival_matched=84
# t_inject = _inject_start_time(steps[1]) = 24 → t_contact=60 = t_travel ✓ σ应相等
# 而实测 σ_on == σ_off —— 恒速井上这是对的！问题在 subtest 断言 σ_on==σ_off 全过，
# 失败的是 '尾浆<-顶替液' σ_off=0.0 的 band —— 那是 REAR_EXIT 同刻 fragment
