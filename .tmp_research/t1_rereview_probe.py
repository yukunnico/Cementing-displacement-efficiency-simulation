# rereview probe: verify F1 fix (dt=0.05 restores discriminating power)
# method: sigma = sqrt(2*D_eff*t)/U is dt-independent; only max(sigma, dt) floor
# depends on dt. If sigma(dt=0.05) == sigma(dt=0.001) exactly, raw sigma > 0.05
# (no clamp). Also reproduce original F1 finding at dt=2.0.
import sys, math
sys.path.insert(0, r"D:\users\desktop\research\控压固井项目\cement model")
sys.path.insert(0, r"D:\users\desktop\research\控压固井项目\cement model\tests")

from test_casing_mixing_contact_time import _synthetic_well, _synthetic_fluids, _dispersion_bands
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep, PumpingStageEvent
from cemdisp.transport1d.casing_flow import CasingFlowSolver

well = _synthetic_well(shoe_md_m=100.0, area_m2=0.01)
fluids = _synthetic_fluids(spacer_rho=None)
schedule = PumpingSchedule(steps=(
    PumpingScheduleStep(step_name="注入隔离液", fluid_name="隔离液", volume_m3=0.4,
                        rate_m3_min=1.0, start_time_s=0.0, end_time_s=24.0,
                        event_tag=PumpingStageEvent.INJECT_SPACER),
    PumpingScheduleStep(step_name="注入尾浆", fluid_name="尾浆", volume_m3=1.0,
                        rate_m3_min=1.0, start_time_s=24.0, end_time_s=84.0,
                        event_tag=PumpingStageEvent.INJECT_CEMENT),
    PumpingScheduleStep(step_name="替浆", fluid_name="顶替液", volume_m3=2.0,
                        rate_m3_min=1.0, start_time_s=84.0, end_time_s=204.0,
                        event_tag=PumpingStageEvent.INJECT_DISPLACEMENT),
))

def bands_at(dt):
    s = CasingFlowSolver(enable_gravity=False, mixing_contact_time=False, dt=dt)
    ev = s.run(well, fluids, schedule).shoe_timeline.events
    return _dispersion_bands(ev)

b20 = bands_at(2.0)
b005 = bands_at(0.05)
b0001 = bands_at(0.001)

print("=== per-interface sigma vs dt ===")
for i, (x20, x005, x0001) in enumerate(zip(b20, b005, b0001)):
    tag = f"{x005[0]}<-{x005[1]}"
    print(f"[{i}] {tag}: sigma(dt=2.0)={x20[2]!r}  sigma(dt=0.05)={x005[2]!r}  sigma(dt=0.001)={x0001[2]!r}")
    print(f"    clamped@2.0: {x20[2] == 2.0}   stable 0.05->0.001: {x005[2] == x0001[2]}   raw>0.05: {x005[2] > 0.05}")

# report claims raw sigma first interface = 1.0077s
print(f"\nfirst-band raw sigma reported 1.0077 -> actual {b005[0][2]}")

# band centers dt-invariance (event skeleton sanity across dt)
centers_ok = all(a[3] == b[3] == c[3] for a, b, c in zip(b20, b005, b0001))
print(f"band centers identical across dt=2.0/0.05/0.001: {centers_ok}")
print(f"n bands: {len(b20)}/{len(b005)}/{len(b0001)}")
