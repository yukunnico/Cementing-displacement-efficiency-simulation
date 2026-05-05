"""
流体物性与流变数据结构

本模块定义了固井施工中各种流体的物性参数和流变模型。

主要类：
- RheologyModel: 支持的流变模型枚举
    * NEWTONIAN: 牛顿流体（粘度恒定）
    * BINGHAM: Bingham塑性流体（屈服应力+塑性粘度）
    * POWER_LAW: 幂律流体（流性指数+稠度系数）
    * HERSCHEL_BULKLEY: Herschel-Bulkley流体（屈服应力+幂律参数）

- FluidRole: 施工程序中流体角色的枚举
    * MUD: 钻井液（被顶替的原地流体）
    * WASH: 清洗液
    * SPACER: 隔离液
    * LEAD: 领浆（前置水泥浆）
    * TAIL: 尾浆（水泥浆主体）
    * DISPLACEMENT: 顶替液
    * OTHER: 其他流体

- FluidSpec: 单种流体的完整物性规格

流体参数约束：
- 密度必须为正数且有限
- 不同流变模型有不同的必填参数要求
- 屈服应力可以为0（表示无屈服特性）
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Optional
import math


class RheologyModel(str, Enum):
    """支持的流变模型枚举。"""

    NEWTONIAN = "newtonian"
    BINGHAM = "bingham"
    POWER_LAW = "power_law"
    HERSCHEL_BULKLEY = "herschel_bulkley"


class FluidRole(str, Enum):
    """施工程序中的流体角色。"""

    MUD = "mud"
    WASH = "wash"
    SPACER = "spacer"
    LEAD = "lead"
    INTERMEDIATE = "intermediate"
    TAIL = "tail"
    DISPLACEMENT = "displacement"
    OTHER = "other"


def _require_finite_positive(name: str, value: float) -> None:
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name}必须为大于0的有限数值")


@dataclass(frozen=True)
class FluidSpec:
    """单种流体的标准物性输入。"""

    name: str
    role: FluidRole
    density_kg_m3: float
    rheology_model: RheologyModel = RheologyModel.NEWTONIAN
    plastic_viscosity_pa_s: Optional[float] = None
    yield_stress_pa: Optional[float] = None
    power_law_n: Optional[float] = None
    consistency_k: Optional[float] = None

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("流体名称不能为空")
        _require_finite_positive("density_kg_m3", self.density_kg_m3)

        if self.rheology_model == RheologyModel.NEWTONIAN:
            if self.plastic_viscosity_pa_s is None:
                raise ValueError("牛顿流体必须提供 plastic_viscosity_pa_s")
            _require_finite_positive("plastic_viscosity_pa_s", self.plastic_viscosity_pa_s)

        elif self.rheology_model == RheologyModel.BINGHAM:
            if self.plastic_viscosity_pa_s is None or self.yield_stress_pa is None:
                raise ValueError("Bingham 流体必须提供 plastic_viscosity_pa_s 和 yield_stress_pa")
            _require_finite_positive("plastic_viscosity_pa_s", self.plastic_viscosity_pa_s)
            if not math.isfinite(self.yield_stress_pa) or self.yield_stress_pa < 0.0:
                raise ValueError("yield_stress_pa 必须为非负有限数值")

        elif self.rheology_model == RheologyModel.POWER_LAW:
            if self.power_law_n is None or self.consistency_k is None:
                raise ValueError("幂律流体必须提供 power_law_n 和 consistency_k")
            _require_finite_positive("power_law_n", self.power_law_n)
            _require_finite_positive("consistency_k", self.consistency_k)

        elif self.rheology_model == RheologyModel.HERSCHEL_BULKLEY:
            if self.power_law_n is None or self.consistency_k is None or self.yield_stress_pa is None:
                raise ValueError("Herschel-Bulkley 流体必须提供 power_law_n、consistency_k 和 yield_stress_pa")
            _require_finite_positive("power_law_n", self.power_law_n)
            _require_finite_positive("consistency_k", self.consistency_k)
            if not math.isfinite(self.yield_stress_pa) or self.yield_stress_pa < 0.0:
                raise ValueError("yield_stress_pa 必须为非负有限数值")
