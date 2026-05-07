
from __future__ import annotations
from dataclasses import dataclass
from typing import Tuple, Callable
from cemdisp.transport1d.casing_flow import CasingFlowSolver
from cemdisp.data.well_spec import WellSpec
from cemdisp.data.fluid_spec import FluidSpec
from cemdisp.data.pumping_schedule import PumpingSchedule

@dataclass(frozen=True)
class AnnulusInletState:
    time_s: float
    fluid_name: str
    role: str
    density_kg_m3: float
    viscosity_pa_s: float
    flow_rate_m3_s: float
    lead_fraction: float
    tail_fraction: float
    spacer_fraction: float


def pipe_exit_to_annulus_inlet(exit_state) -> AnnulusInletState:
    return AnnulusInletState(
        time_s=exit_state.time_s,
        fluid_name=exit_state.fluid_name,
        role=exit_state.role,
        density_kg_m3=exit_state.density_kg_m3,
        viscosity_pa_s=exit_state.viscosity_pa_s,
        flow_rate_m3_s=exit_state.flow_rate_m3_s,
        lead_fraction=exit_state.lead_fraction,
        tail_fraction=exit_state.tail_fraction,
        spacer_fraction=exit_state.spacer_fraction,
    )


def build_coupled_annulus_inlet_provider(
    solver: CasingFlowSolver,
    well: WellSpec,
    fluids: Tuple[FluidSpec, ...],
    schedule: PumpingSchedule,
) -> Callable[[float], AnnulusInletState]:
    def provider(time_s: float) -> AnnulusInletState:
        return pipe_exit_to_annulus_inlet(solver.exit_state(time_s, well, fluids, schedule))
    return provider
