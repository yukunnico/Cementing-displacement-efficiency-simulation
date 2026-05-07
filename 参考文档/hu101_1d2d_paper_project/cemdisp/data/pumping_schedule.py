
from __future__ import annotations
from dataclasses import dataclass
from typing import Tuple

@dataclass(frozen=True)
class PumpingScheduleStep:
    step_name: str
    fluid_name: str
    volume_m3: float
    rate_m3_min: float
    remarks: str = ""

@dataclass(frozen=True)
class PumpingSchedule:
    steps: Tuple[PumpingScheduleStep, ...]

    def total_volume_m3(self) -> float:
        return sum(s.volume_m3 for s in self.steps)

    def total_time_s(self) -> float:
        return sum(s.volume_m3 / s.rate_m3_min * 60.0 for s in self.steps)

    def cumulative_boundaries(self):
        boundaries=[]
        v=0.0
        t=0.0
        for step in self.steps:
            duration=step.volume_m3/step.rate_m3_min*60.0
            boundaries.append((v, v+step.volume_m3, t, t+duration, step))
            v+=step.volume_m3
            t+=duration
        return boundaries

    def time_at_cumulative_volume(self, volume_m3: float) -> float:
        vtarget=max(0.0, min(float(volume_m3), self.total_volume_m3()))
        for v0,v1,t0,t1,step in self.cumulative_boundaries():
            if vtarget <= v1 + 1e-12:
                return t0 + (vtarget - v0) / step.rate_m3_min * 60.0
        return self.total_time_s()

    def state_at_cumulative_volume(self, volume_m3: float):
        vtarget=float(volume_m3)
        if vtarget < 0.0:
            return None
        for v0,v1,t0,t1,step in self.cumulative_boundaries():
            if v0 <= vtarget < v1 or (abs(vtarget-v1)<1e-9 and v1>=self.total_volume_m3()):
                return step
        return self.steps[-1]

    def cumulative_volume_at_time(self, time_s: float) -> float:
        t=float(time_s)
        if t <= 0: return 0.0
        for v0,v1,t0,t1,step in self.cumulative_boundaries():
            if t <= t1 + 1e-12:
                return v0 + (t-t0)*step.rate_m3_min/60.0
        return self.total_volume_m3()
