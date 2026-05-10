"""
鞋口出流时间轴抽象

本模块定义鞋口事件的时序抽象，用于为环空二维模型提供鞋口出流边界条件。
鞋口时间轴（ShoeTimeline）按时间顺序记录一系列鞋口事件（ ShoeEvent），
并支持在任意时刻查询对应的鞋口出流状态（PipeExitState）。

核心概念：
- ShoeEventKind：鞋口事件类型（FRONT_ARRIVAL / REAR_EXIT / RATE_SWITCH / SHUTDOWN / RESTART / END）
- ShoeEvent：单个鞋口事件快照，包含时间、类型、排量、阶段名、液相分数
- ShoeTimeline：按时间轴管理事件，提供 .at(time_s) 查询

用途：
- 作为套管内1D输运层向环空2D层传递的鞋口出流边界
- 支持施工阶段切换、排量变化、启停泵等事件的建模
- 为后续 CasingFlowSolver 提供统一的鞋口出流状态查询接口
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from cemdisp.transport1d.pipe_exit_state import PipeExitState


class ShoeEventKind(Enum):
    """鞋口事件类型枚举。

    用于标识 ShoeEvent 的具体事件类别。
    """

    FRONT_ARRIVAL = "FRONT_ARRIVAL"   # 前缘到达鞋口
    REAR_EXIT = "REAR_EXIT"           # 尾部离开鞋口
    RATE_SWITCH = "RATE_SWITCH"       # 排量切换
    SHUTDOWN = "SHUTDOWN"             # 停泵
    RESTART = "RESTART"               # 重新启动
    END = "END"                       # 时间轴结束


@dataclass(frozen=True)
class ShoeEvent:
    """鞋口处单个事件快照。

    记录某个时刻鞋口处发生的单一事件及其当时的出流状态。
    设计为不可变对象，方便在时间轴中安全传递和比较。

    属性：
    - time_s: 事件发生时间（秒）
    - kind: 事件类型（ShoeEventKind）
    - flow_rate_m3_s: 事件发生时的排量（立方米/秒）
    - stage_name: 当前施工阶段名称
    - phase_fractions: 各相流体的体积分数元组，每项为 (流体名, 分数)
    """

    time_s: float
    kind: ShoeEventKind
    flow_rate_m3_s: float
    stage_name: str
    phase_fractions: tuple[tuple[str, float], ...] = ()  # 默认空元组


class ShoeTimeline:
    """鞋口事件时间轴。

    按时间顺序存储一系列 ShoeEvent，支持在任意时刻查询对应的鞋口出流状态。
    时间轴内部按 time_s 升序排列，.at() 查询时通过线性扫描定位最近事件。

    设计约束：
    - 事件时间必须单调递增，不允许乱序插入
    - 空时间轴在任意时刻返回零排量、空阶段名的默认 PipeExitState

    示例：
        events = [
            ShoeEvent(time_s=0.0, kind=ShoeEventKind.FRONT_ARRIVAL,
                      flow_rate_m3_s=0.0, stage_name="初始"),
            ShoeEvent(time_s=60.0, kind=ShoeEventKind.FRONT_ARRIVAL,
                      flow_rate_m3_s=0.02, stage_name="注入领浆",
                      phase_fractions=(("领浆", 1.0),)),
            ShoeEvent(time_s=180.0, kind=ShoeEventKind.REAR_EXIT,
                      flow_rate_m3_s=0.0, stage_name="停泵"),
        ]
        tl = ShoeTimeline(events=events)
        state = tl.at(time_s=100.0)  # 返回 time_s=60.0 事件的状态
    """

    def __init__(self, events: list[ShoeEvent]) -> None:
        """初始化鞋口事件时间轴。

        Args:
            events: 按 time_s 升序排列的 ShoeEvent 列表。
                   空列表表示空时间轴。

        Raises:
            ValueError: 如果事件列表中任意两个事件的时间不满足递增顺序。
        """
        # 将列表转换为不可变元组，保证安全性
        if events:
            # 验证时间单调递增
            for i in range(1, len(events)):
                if events[i].time_s < events[i - 1].time_s:
                    msg = f"事件时间必须单调递增，但 events[{i-1}].time_s={events[i-1].time_s} > events[{i}].time_s={events[i].time_s}"
                    raise ValueError(msg)
        self._events: tuple[ShoeEvent, ...] = tuple(events)

    @property
    def events(self) -> tuple[ShoeEvent, ...]:
        """返回时间轴中所有事件的只读视图。"""
        return self._events

    def at(self, time_s: float) -> PipeExitState:
        """查询任意时刻的鞋口出流状态。

        查找时间轴中最近过去（time_s 之前，或恰好等于 time_s）的事件，
        并返回该事件对应的 PipeExitState。

        两种边界情况：
        - 若时间轴为空：返回零排量、空阶段名的默认 PipeExitState（time_s 强制为 0.0）。
        - 若查询时刻早于所有事件：返回首个事件的状态（最近过去）。

        算法：线性扫描，时间复杂度 O(n)。事件数量有限时足够高效。

        Args:
            time_s: 查询时刻（秒），允许任意数值。

        Returns:
            PipeExitState：最近过去事件的出流状态。
        """
        if not self._events:
            # 空时间轴：返回默认零状态
            return PipeExitState(
                time_s=0.0,
                flow_rate_m3_s=0.0,
                stage_name="",
                phase_fractions=(),
            )

        # 查找最近过去事件
        candidate: ShoeEvent | None = None
        for event in self._events:
            if event.time_s <= time_s:
                candidate = event
            else:
                # 事件按时间升序排列，一旦遇到第一个大于 time_s 的事件，
                # 之后的都不会更接近，直接 break
                break

        if candidate is None:
            # 查询时刻早于所有事件：返回第一个事件（最近过去）
            candidate = self._events[0]

        return PipeExitState(
            time_s=candidate.time_s,
            flow_rate_m3_s=candidate.flow_rate_m3_s,
            stage_name=candidate.stage_name,
            phase_fractions=candidate.phase_fractions,
        )