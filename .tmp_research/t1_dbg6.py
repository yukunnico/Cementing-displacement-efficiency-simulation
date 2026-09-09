import sys; sys.path.insert(0, 'tests')
from test_casing_mixing_contact_time import _synthetic_well, _synthetic_fluids
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep, PumpingStageEvent
from cemdisp.transport1d.casing_flow import CasingFlowSolver
well = _synthetic_well(100.0, 0.01)
fluids = _synthetic_fluids(spacer_rho=None)
sched = PumpingSchedule(steps=(
  PumpingScheduleStep(step_name='尾浆一段', fluid_name='尾浆', volume_m3=0.5, rate_m3_min=1.0, start_time_s=0.0, end_time_s=30.0, event_tag=PumpingStageEvent.INJECT_CEMENT),
  PumpingScheduleStep(step_name='注入隔离液', fluid_name='隔离液', volume_m3=0.5, rate_m3_min=1.0, start_time_s=30.0, end_time_s=60.0, event_tag=PumpingStageEvent.INJECT_SPACER),
  PumpingScheduleStep(step_name='尾浆二段', fluid_name='尾浆', volume_m3=1.0, rate_m3_min=1.0, start_time_s=60.0, end_time_s=120.0, event_tag=PumpingStageEvent.INJECT_CEMENT),
  PumpingScheduleStep(step_name='替浆', fluid_name='顶替液', volume_m3=2.0, rate_m3_min=1.0, start_time_s=120.0, end_time_s=240.0, event_tag=PumpingStageEvent.INJECT_DISPLACEMENT)))
s = CasingFlowSolver(enable_gravity=False, mixing_contact_time=True)
steps, _ = s._displacement_sequence_cutoff(s._build_scheduled_steps(sched))
for i, st in enumerate(steps):
    print(i, st.step.fluid_name, 'cum:', st.cumulative_volume_start_m3, '->', st.cumulative_volume_end_m3, 'front:', s._front_arrival_time(st, steps, 1.0))
print('t_inject seg2:', s._inject_start_time(steps[2], steps))
