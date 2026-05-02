"""
套管内一维输运层

本模块实现套管内流体一维输运计算，用于追踪从地面注入的流体在套管内的下行过程。

核心功能：
- 流体前缘追踪（InterfaceTracker）：按体积推进法追踪各流体前缘位置
- 前沿计算结果（CasingFlowResult）：存储各流体前缘的最终位置和到达时间
- 鞋口出流状态（PipeExitState）：给定时刻鞋口处流出的流体类型和排量

设计原则：
- 套管内1D输运层提供更真实的鞋口出流边界，不污染环空2D核心
- 环空求解器可独立运行，无需依赖套管内1D层
- 1D→2D边界桥接接口清晰，便于后续扩展

主要类：
- InterfaceFront: 单个流体前缘在套管内的位置（深度+时间）
- InterfaceTracker: 时间步长形式的前缘推进追踪器
- CasingFlowResult: 套管内1D输运的解析计算结果
- CasingFlowSolver: 套管内1D前沿追踪求解器（支持解析计算）
- PipeExitState: 鞋口处实时出流状态
"""

from cemdisp.transport1d.casing_flow import CasingFlowResult, CasingFlowSolver
from cemdisp.transport1d.interface_tracking import InterfaceFront, InterfaceTracker
from cemdisp.transport1d.pipe_exit_state import PipeExitState

__all__ = [
    "CasingFlowResult",
    "CasingFlowSolver",
    "InterfaceFront",
    "InterfaceTracker",
    "PipeExitState",
]
