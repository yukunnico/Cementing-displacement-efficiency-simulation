
from __future__ import annotations
from dataclasses import dataclass
from typing import Tuple
from cemdisp.data.fluid_spec import FluidSpec, FluidRole
from cemdisp.data.pumping_schedule import PumpingSchedule
from cemdisp.data.well_spec import WellSpec

@dataclass(frozen=True)
class InterfaceFront:
    fluid_name: str
    role: str
    start_volume_m3: float
    end_volume_m3: float
    shoe_start_time_s: float
    shoe_end_time_s: float

@dataclass(frozen=True)
class PipeExitState:
    time_s: float
    fluid_name: str
    role: str
    density_kg_m3: float
    viscosity_pa_s: float
    flow_rate_m3_s: float
    cement_fraction: float
    lead_fraction: float
    tail_fraction: float
    spacer_fraction: float

@dataclass(frozen=True)
class CasingFlowResult:
    pipe_volume_m3: float
    fronts: Tuple[InterfaceFront, ...]
    cement_arrival_time_s: float
    cement_end_time_s: float
    displacement_arrival_time_s: float

class CasingFlowSolver:
    """套管段一维体积前沿追踪。

    采用不可压缩活塞流假设，地面施工程序整体滞后 pipe_volume 后在鞋口输出。
    该层只负责把地面泵注时序转换为鞋口环空入口边界，不引入环空额外工程修正。
    """
    def run(self, well: WellSpec, fluids: Tuple[FluidSpec, ...], schedule: PumpingSchedule) -> CasingFlowResult:
        by_name={f.name:f for f in fluids}
        fronts=[]
        cement_arrival=None
        cement_end=None
        displacement_arrival=None
        for v0,v1,t0,t1,step in schedule.cumulative_boundaries():
            fluid=by_name[step.fluid_name]
            ts=schedule.time_at_cumulative_volume(v0 + well.casing_lag_volume_m3)
            te=schedule.time_at_cumulative_volume(v1 + well.casing_lag_volume_m3)
            fronts.append(InterfaceFront(fluid.name, fluid.role.value, v0, v1, ts, te))
            if fluid.role in {FluidRole.LEAD, FluidRole.TAIL} and cement_arrival is None:
                cement_arrival=ts
            if fluid.role in {FluidRole.LEAD, FluidRole.TAIL}:
                cement_end=te
            if fluid.role == FluidRole.DISPLACEMENT and displacement_arrival is None:
                displacement_arrival=ts
        return CasingFlowResult(
            pipe_volume_m3=well.casing_lag_volume_m3,
            fronts=tuple(fronts),
            cement_arrival_time_s=float(cement_arrival or 0.0),
            cement_end_time_s=float(cement_end or schedule.total_time_s()),
            displacement_arrival_time_s=float(displacement_arrival or schedule.total_time_s()),
        )

    def exit_state(self, time_s: float, well: WellSpec, fluids: Tuple[FluidSpec, ...], schedule: PumpingSchedule) -> PipeExitState:
        by_name={f.name:f for f in fluids}
        pumped=schedule.cumulative_volume_at_time(time_s)
        outlet_volume=pumped - well.casing_lag_volume_m3
        if outlet_volume < 0.0:
            fluid=by_name["钻井液"]
            q=0.0
            step_name="套管滞后期"
        else:
            step=schedule.state_at_cumulative_volume(outlet_volume)
            fluid=by_name[step.fluid_name]
            q=step.rate_m3_min/60.0
            step_name=step.step_name
        lead=1.0 if fluid.role == FluidRole.LEAD else 0.0
        tail=1.0 if fluid.role == FluidRole.TAIL else 0.0
        spacer=1.0 if fluid.role in {FluidRole.WASH, FluidRole.SPACER} else 0.0
        cement=lead+tail
        return PipeExitState(
            time_s=float(time_s), fluid_name=fluid.name, role=fluid.role.value,
            density_kg_m3=fluid.density_kg_m3, viscosity_pa_s=fluid.apparent_viscosity(100.0),
            flow_rate_m3_s=q, cement_fraction=cement, lead_fraction=lead, tail_fraction=tail, spacer_fraction=spacer,
        )
