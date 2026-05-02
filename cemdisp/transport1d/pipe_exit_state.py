"""
鞋口处出流状态数据结构

本模块定义了套管鞋口处实时出流状态的数据结构。

主要字段：
- time_s: 当前时刻（秒）
- flow_rate_m3_s: 当前排量（立方米/秒）
- stage_name: 当前施工阶段名称
- phase_fractions: 各相流体的体积分数元组

用途：
- 作为套管内1D输运层向环空2D层传递的边界信息
- 记录某个时刻从套管鞋口流入环空的流体类型和排量
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Tuple


@dataclass(frozen=True)
class PipeExitState:
    """鞋口处实时出流状态。

    记录任意时刻从套管鞋口流出的流体状态，
    供环空二维模型作为入口边界条件使用。
    """

    time_s: float
    flow_rate_m3_s: float
    stage_name: str
    phase_fractions: Tuple[Tuple[str, float], ...] = field(default_factory=tuple)
