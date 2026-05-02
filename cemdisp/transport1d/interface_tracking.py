"""
套管内流体前缘追踪数据结构与轻量追踪器

本模块定义了套管内流体前缘追踪的核心数据结构和辅助工具类。

主要类：
- InterfaceFront: 单个流体前缘的位置快照
    * fluid_name: 流体名称
    * distance_m: 从地面算起的前缘深度（米）
    * time_s: 到达该位置的时间（秒）

- InterfaceTracker: 时间步长形式的前缘推进追踪器
    * advance_front(): 按排量推进指定流体的前缘
    * fronts_snapshot(): 获取所有前缘的当前状态快照
    * fluid_at_shoe: 当前鞋口处正在流出的流体名

设计说明：
- 主求解器(CasingFlowSolver)当前采用解析计算
- InterfaceTracker保留用于可视化或非恒定排量细化场景
- 位移增量 = 排量(m³/s) × 时间步长(s) / 管内截面积(m²)
"""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class InterfaceFront:
    """单个流体前缘在套管内的位置快照。

    记录某个时刻某流体前缘的深度和到达时间。
    """

    fluid_name: str
    distance_m: float
    time_s: float


class InterfaceTracker:
    """追踪多个流体前缘在套管内的推进。

    该类用于时间步长形式的前缘推进；主求解器当前采用解析计算，
    但保留这个轻量追踪器，便于后续接入可视化或非恒定排量细化。
    """

    def __init__(self, shoe_depth_m: float, pipe_area_m2: float) -> None:
        """初始化前缘追踪器。

        Args:
            shoe_depth_m: 鞋口深度（米）
            pipe_area_m2: 管内截面积（平方米）
        """
        if not math.isfinite(shoe_depth_m) or shoe_depth_m <= 0.0:
            raise ValueError("shoe_depth_m 必须为大于0的有限数值")
        if not math.isfinite(pipe_area_m2) or pipe_area_m2 <= 0.0:
            raise ValueError("pipe_area_m2 必须为大于0的有限数值")
        self.shoe_depth_m: float = shoe_depth_m
        self.pipe_area_m2: float = pipe_area_m2
        self._front_distance_by_fluid: dict[str, float] = {}
        self._arrival_time_by_fluid: dict[str, float] = {}
        self._current_time_s: float = 0.0
        self._fluid_at_shoe: str = ""

    def advance_front(self, fluid_name: str, rate_m3_s: float, dt: float) -> None:
        """推进指定流体的前缘。

        Args:
            fluid_name: 需要推进的流体名称。
            rate_m3_s: 当前排量，单位 m³/s。
            dt: 时间步长，单位 s。
        """

        if not fluid_name.strip():
            raise ValueError("fluid_name 不能为空")
        if not math.isfinite(rate_m3_s) or rate_m3_s < 0.0:
            raise ValueError("rate_m3_s 必须为非负有限数值")
        if not math.isfinite(dt) or dt < 0.0:
            raise ValueError("dt 必须为非负有限数值")

        # 体积推进法：位移增量 = 排量 × 时间 / 管内截面积。
        previous_distance = self._front_distance_by_fluid.get(fluid_name, 0.0)
        distance_increment = rate_m3_s * dt / self.pipe_area_m2
        new_distance = min(self.shoe_depth_m, previous_distance + distance_increment)
        self._current_time_s += dt
        self._front_distance_by_fluid[fluid_name] = new_distance

        # 首次到达鞋口时记录到达时间，并把鞋口流体切换为最近到达流体。
        if new_distance >= self.shoe_depth_m and fluid_name not in self._arrival_time_by_fluid:
            self._arrival_time_by_fluid[fluid_name] = self._current_time_s
            self._fluid_at_shoe = fluid_name

    def fronts_snapshot(self, time_s: float) -> tuple[InterfaceFront, ...]:
        """返回当前所有前缘的状态快照。"""

        if not math.isfinite(time_s) or time_s < 0.0:
            raise ValueError("time_s 必须为非负有限数值")
        return tuple(
            InterfaceFront(fluid_name=fluid_name, distance_m=distance_m, time_s=time_s)
            for fluid_name, distance_m in sorted(self._front_distance_by_fluid.items())
        )

    @property
    def fluid_at_shoe(self) -> str:
        """返回当前鞋口处正在流出的流体名。"""

        return self._fluid_at_shoe
