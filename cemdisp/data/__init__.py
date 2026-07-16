"""
cemdisp.data - 标准输入数据结构模块

本模块定义了固井顶替模型所需的所有标准输入数据结构，包括：
- WellSpec: 单口井的井筒几何参数和评价井段
- FluidSpec: 单种流体的物性参数（密度、流变模型参数）
- PumpingSchedule: 地面泵注施工程序（注入步骤时序）
- ValidationData: 与单井相关的现场校验资料路径（CBL报告、施工记录等）

这些数据结构采用 frozen dataclass 形式，确保数据不可变性和合法性校验。
所有数据均不硬编码单井具体数值，便于多井复用。
"""

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.provenance import (
    FluidProvenance,
    SectionProvenance,
    WellProvenance,
    WELL_PROVENANCE,
    build_injected_fluid_provenance_summary,
    format_injected_fluid_provenance_markdown,
)
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep
from cemdisp.data.validation_data import ValidationData
from cemdisp.data.well_spec import DepthValuePoint, EvaluationWindow, WellSpec

__all__ = [
    "DepthValuePoint",
    "EvaluationWindow",
    "FluidRole",
    "FluidProvenance",
    "FluidSpec",
    "PumpingSchedule",
    "PumpingScheduleStep",
    "RheologyModel",
    "SectionProvenance",
    "ValidationData",
    "WellProvenance",
    "WellSpec",
    "WELL_PROVENANCE",
    "build_injected_fluid_provenance_summary",
    "format_injected_fluid_provenance_markdown",
]
