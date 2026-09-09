import sys; sys.path.insert(0, 'tests')
from test_casing_mixing_contact_time import _synthetic_well, _synthetic_fluids
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep, PumpingStageEvent
from cemdisp.transport1d.casing_flow import CasingFlowSolver
shoe, area = 100.0, 0.01
well = _synthetic_well(shoe, area)
fluids = _synthetic_fluids(spacer_rho=None)
sched = PumpingSchedule(steps=(
  PumpingScheduleStep(step_name='注入隔离液', fluid_name='隔离液', volume_m3=0.4, rate_m3_min=1.0, start_time_s=0.0, end_time_s=24.0, event_tag=PumpingStageEvent.INJECT_SPACER),
  PumpingScheduleStep(step_name='注入尾浆', fluid_name='尾浆', volume_m3=1.0, rate_m3_min=1.0, start_time_s=24.0, end_time_s=84.0, event_tag=PumpingStageEvent.INJECT_CEMENT),
  PumpingScheduleStep(step_name='替浆', fluid_name='顶替液', volume_m3=2.0, rate_m3_min=1.0, start_time_s=84.0, end_time_s=204.0, event_tag=PumpingStageEvent.INJECT_DISPLACEMENT)))
s = CasingFlowSolver(enable_gravity=False, mixing_contact_time=True)
steps_full = s._build_scheduled_steps(sched)
steps, _ = s._displacement_sequence_cutoff(steps_full)
tail = steps[1]
print('tail cum start/end:', tail.cumulative_volume_start_m3, tail.cumulative_volume_end_m3)
print('t_inject tail:', s._inject_start_time(tail, steps))
pv = s._timeline_pipe_volume(well, well.shoe_md_m * 0.01)
print('pipe_vol:', pv, 't_arrival tail:', s._front_arrival_time(tail, steps, pv))
r_liner = well.liner_id_mm/2000.0
U = (1.0/60.0)/(3.141592653589793*r_liner**2)
print('U:', U, 't_travel:', shoe/U)
import math
D = s._compute_dispersion_coefficient(r_liner, fluids[2], U)
inst = s._interface_instability_factor(fluids[2], fluids[3], r_liner, U)
print('D:', D, 'instability:', inst, 'D*inst:', D*inst)
D_eff = D*inst
sig_off = math.sqrt(2*D_eff*(shoe/U))/U
sig_on = s._contact_time_integrated_sigma(60.0, s._inject_start_time(tail, steps), shoe/U, D_eff, U, s.dt)
print('sig_off raw:', sig_off, '-> clamp:', max(min(sig_off, 0.5*(shoe/U)), s.dt))
print('sig_on  raw:', sig_on, '-> clamp:', max(min(sig_on, 0.5*(shoe/U)), s.dt))
