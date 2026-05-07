
from __future__ import annotations
from dataclasses import dataclass
from enum import Enum

class FluidRole(str, Enum):
    MUD = "mud"
    WASH = "wash"
    SPACER = "spacer"
    LEAD = "lead"
    TAIL = "tail"
    DISPLACEMENT = "displacement"
    OTHER = "other"

class RheologyModel(str, Enum):
    NEWTONIAN = "newtonian"
    BINGHAM = "bingham"
    POWER_LAW = "power_law"

@dataclass(frozen=True)
class FluidSpec:
    name: str
    role: FluidRole
    density_kg_m3: float
    rheology_model: RheologyModel
    plastic_viscosity_pa_s: float = 0.0
    yield_stress_pa: float = 0.0
    power_law_n: float = 1.0
    consistency_k: float = 0.001

    def apparent_viscosity(self, shear_rate_s: float = 100.0) -> float:
        gamma = max(float(shear_rate_s), 1.0e-6)
        if self.rheology_model == RheologyModel.NEWTONIAN:
            return max(self.plastic_viscosity_pa_s, 1e-5)
        if self.rheology_model == RheologyModel.BINGHAM:
            return min(max(self.plastic_viscosity_pa_s + self.yield_stress_pa / gamma, 1e-5), 5.0)
        if self.rheology_model == RheologyModel.POWER_LAW:
            return min(max(self.consistency_k * gamma ** (self.power_law_n - 1.0), 1e-5), 5.0)
        raise ValueError(f"Unsupported rheology model: {self.rheology_model}")
