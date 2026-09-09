import sys; sys.path.insert(0, 'tests')
from test_casing_mixing_contact_time import _synthetic_well, _synthetic_fluids, _dispersion_bands
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep, PumpingStageEvent
from cemdisp.transport1d.casing_flow import CasingFlowSolver
shoe, area = 100.0, 0.01
well = _synthetic_well(shoe, area)
fluids = _synthetic_fluids(spacer_rho=None)
sched = PumpingSchedule(steps=(
  PumpingScheduleStep(step_name='注入隔离液', fluid_name='隔离液', volume_m3=0.4, rate_m3_min=1.0, start_time_s=0.0, end_time_s=24.0, event_tag=PumpingStageEvent.INJECT_SPACER),
  PumpingScheduleStep(step_name='注入尾浆', fluid_name='尾浆', volume_m3=1.0, rate_m3_min=1.0, start_time_s=24.0, end_time_s=84.0, event_tag=PumpingStageEvent.INJECT_CEMENT),
  PumpingScheduleStep(step_name='替浆', fluid_name='顶替液', volume_m3=2.0, rate_m3_min=1.0, start_time_s=84.0, end_time_s=204.0, event_tag=PumpingStageEvent.INJECT_DISPLACEMENT)))
for flag in (False, True):
    s = CasingFlowSolver(enable_gravity=False, mixing_contact_time=flag)
    r = s.run(well, fluids, sched)
    print(f'--- flag={flag} events={len(r.shoe_timeline.events)}')
    for e in r.shoe_timeline.events:
        pf = e.phase_fractions
        print(f'  {e.time_s:10.4f} {e.kind.value:13s} {pf}')
