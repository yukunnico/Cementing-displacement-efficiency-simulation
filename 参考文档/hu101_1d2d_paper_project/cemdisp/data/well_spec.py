
from __future__ import annotations
from dataclasses import dataclass
from typing import Tuple

@dataclass(frozen=True)
class DepthValuePoint:
    depth_md_m: float
    value: float

@dataclass(frozen=True)
class EvaluationWindow:
    name: str
    top_md_m: float
    bottom_md_m: float
    window_type: str

@dataclass(frozen=True)
class WellSpec:
    well_name: str
    top_md_m: float
    bottom_md_m: float
    shoe_md_m: float
    hanger_md_m: float
    liner_od_mm: float
    liner_id_mm: float
    casing_lag_volume_m3: float
    upper_liner_od_mm: float
    lower_liner_od_mm: float
    upper_lower_transition_md_m: float
    hole_diameter_profile: Tuple[DepthValuePoint, ...]
    inclination_profile: Tuple[DepthValuePoint, ...]
    standoff_profile: Tuple[DepthValuePoint, ...]
    evaluation_windows: Tuple[EvaluationWindow, ...]
    notes: str = ""
