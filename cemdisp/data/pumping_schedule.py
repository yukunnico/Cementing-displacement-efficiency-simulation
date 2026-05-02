"""
施工程序与地面泵注时序数据结构

本模块定义了地面泵注施工程序的数据结构，用于表达现场分段注入、排量切换、停泵/续泵等操作。

主要类：
- PumpingScheduleStep: 地面施工程序中的单个注入步骤
    * step_name: 步骤名称（如"注入尾管水泥浆"、"替浆推进"）
    * fluid_name: 流体名称（需与FluidSpec.name对应）
    * volume_m3: 注入体积（立方米）
    * rate_m3_min: 泵注排量（立方米/分钟）
    * start_time_s: 可选的步骤开始时间（秒），不提供则按前序步骤结束时间自动计算
    * end_time_s: 可选的步骤结束时间（秒）
    * remarks: 备注信息

- PumpingSchedule: 整套地面施工程序，由多个PumpingScheduleStep组成

约束条件：
- 至少需要一个施工步骤
- 所有步骤要么全部提供显式时间，要么全部由排量自动计算
- 步骤的开始时间必须非递减
- 若有实际注入体积，排量必须大于零
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Tuple
import math


def _require_finite_non_negative(name: str, value: float) -> None:
    """校验数值是否非负且有限的内部辅助函数。"""
    if not math.isfinite(value) or value < 0.0:
        raise ValueError(f"{name}必须为非负有限数值")


def _require_finite_positive(name: str, value: float) -> None:
    """校验数值是否大于零且有限的内部辅助函数。"""
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name}必须为大于0的有限数值")


@dataclass(frozen=True)
class PumpingScheduleStep:
    """地面施工程序中的单个注入步骤。"""

    step_name: str
    fluid_name: str
    volume_m3: float
    rate_m3_min: float
    start_time_s: Optional[float] = None
    end_time_s: Optional[float] = None
    remarks: str = ""

    def __post_init__(self) -> None:
        if not self.step_name.strip():
            raise ValueError("step_name不能为空")
        if not self.fluid_name.strip():
            raise ValueError("fluid_name不能为空")
        _require_finite_non_negative("volume_m3", self.volume_m3)
        _require_finite_non_negative("rate_m3_min", self.rate_m3_min)
        if self.start_time_s is not None:
            _require_finite_non_negative("start_time_s", self.start_time_s)
        if self.end_time_s is not None:
            _require_finite_non_negative("end_time_s", self.end_time_s)
        if self.start_time_s is not None and self.end_time_s is not None:
            if self.end_time_s <= self.start_time_s:
                raise ValueError("end_time_s 必须大于 start_time_s")
        if self.volume_m3 > 0.0 and self.rate_m3_min <= 0.0:
            raise ValueError("存在实际注入体积时，rate_m3_min 必须大于0")


@dataclass(frozen=True)
class PumpingSchedule:
    """整套地面施工程序。"""

    steps: Tuple[PumpingScheduleStep, ...] = field(default_factory=tuple)
    notes: Tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if not self.steps:
            raise ValueError("PumpingSchedule 至少需要一个步骤")
        explicit_start_times = [step.start_time_s for step in self.steps if step.start_time_s is not None]
        if explicit_start_times and len(explicit_start_times) != len(self.steps):
            raise ValueError("若使用显式时间，所有步骤都必须提供 start_time_s")
        if len(explicit_start_times) == len(self.steps):
            ordered_starts = [step.start_time_s for step in self.steps]
            assert ordered_starts[0] is not None
            for previous, current in zip(ordered_starts, ordered_starts[1:]):
                assert previous is not None and current is not None
                if current < previous:
                    raise ValueError("施工步骤的 start_time_s 必须非递减")

    @property
    def total_volume_m3(self) -> float:
        return sum(step.volume_m3 for step in self.steps)
