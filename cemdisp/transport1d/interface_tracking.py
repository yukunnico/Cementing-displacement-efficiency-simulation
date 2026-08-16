"""
套管内流体前缘追踪数据结构

本模块定义套管内流体前缘的位置快照数据结构。

主要类：
- InterfaceFront: 单个流体前缘的位置快照
    * fluid_name: 流体名称
    * distance_m: 从地面算起的前缘深度（米）
    * time_s: 到达该位置的时间（秒）
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class InterfaceFront:
    """单个流体前缘在套管内的位置快照。

    记录某个时刻某流体前缘的深度和到达时间。
    """

    fluid_name: str
    distance_m: float
    time_s: float
