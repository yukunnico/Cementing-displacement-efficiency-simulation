"""
鞋口出流到环空入口边界的桥接结构

本模块实现从套管内1D输运层到环空2D层的边界桥接功能。

核心功能：
1. AnnulusInletState：环空入口边界状态数据结构
    - time_s: 当前时间（秒）
    - flow_rate_m3_s: 排量（立方米/秒）
    - stage_name: 施工阶段名称
    - phase_fractions: 各相体积分数（如水泥100%或钻井液100%）

2. pipe_exit_to_annulus_inlet()：直接将鞋口出流状态映射为环空入口
    - 用于简单的直接边界传递

3. build_coupled_annulus_inlet_provider()：构建1D-2D耦合的边界提供器
    - 将套管鞋口出流状态转换为环空入口边界
    - 自动将流体名称映射为环空三相模型（cement/spacer/mud）
    - 保持1D层可插拔，不改变环空2D求解器接口

耦合边界口径：
- 水泥前缘尚未到达鞋口前，环空入口看到的是被推出套管的初始钻井液；
- 水泥前缘到达鞋口后，环空入口才切换为水泥相；
- 因此效率偏低时应优先检查环空二维求解总时长是否足够，
  而不是提前把地面注入的水泥直接施加到环空入口。

设计原则：
- 环空2D求解器不感知1D层实现细节
- 流体名称到相分数的映射通过FluidRole枚举自动完成
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.transport1d.casing_flow import CasingFlowResult, CasingFlowSolver
from cemdisp.transport1d.pipe_exit_state import PipeExitState


@dataclass(frozen=True)
class AnnulusInletState:
    """环空入口边界状态。

    描述环空顶部入口处任意时刻的流体状态，
    用于环空二维求解器的边界条件输入。
    """

    time_s: float
    flow_rate_m3_s: float
    stage_name: str
    phase_fractions: tuple[tuple[str, float], ...] = field(default_factory=tuple)


def pipe_exit_to_annulus_inlet(pipe_exit_state: PipeExitState) -> AnnulusInletState:
    """把鞋口出流状态直接映射为环空入口边界状态。"""

    return AnnulusInletState(
        time_s=pipe_exit_state.time_s,
        flow_rate_m3_s=pipe_exit_state.flow_rate_m3_s,
        stage_name=pipe_exit_state.stage_name,
        phase_fractions=pipe_exit_state.phase_fractions,
    )


def build_coupled_annulus_inlet_provider(
    casing_result: CasingFlowResult,
    casing_solver: CasingFlowSolver,
    fluids: tuple[FluidSpec, ...],
    *,
    split_cement_phases: bool = False,
) -> Callable[[float], AnnulusInletState]:
    """构建1D-2D耦合的环空入口边界提供器。

    该桥接函数只负责把套管鞋口出流状态转换为环空入口边界，
    不改变环空二维求解器接口，从而保持 1D 层可插拔。

    注意：这里采用“鞋口实际出流”作为环空入口边界。水泥从地面进入套管
    到到达鞋口之间存在管内容积延迟，因此延迟期间环空入口仍为初始钻井液；
    水泥前缘到达鞋口后才映射为 cement 相。若需要让水泥在环空内有足够
    顶替时间，应延长环空二维求解器 total_t（例如 12000 s），而不是修改
    本桥接逻辑或相分数映射。
    """

    role_by_name: dict[str, FluidRole] = {fluid.name: fluid.role for fluid in fluids}

    def _phase_fractions_for_fluid(fluid_name: str) -> tuple[tuple[str, float], ...]:
        """根据流体角色映射为环空二维模型入口相分数。"""

        role = role_by_name.get(fluid_name, FluidRole.MUD)
        # 默认保持历史三相口径：领浆/尾浆归为 cement，相容于 Hu102 既有逻辑。
        # 当 split_cement_phases=True 时，领浆/尾浆分开送入 2D 层，供 Hu101 更真实地跟踪
        # lead → tail → mud invasion 的现场顺序。
        if split_cement_phases and role == FluidRole.LEAD:
            return (("lead", 1.0),)
        if split_cement_phases and role == FluidRole.TAIL:
            return (("tail", 1.0),)
        if role in {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}:
            return (("cement", 1.0),)
        if role in {FluidRole.WASH, FluidRole.SPACER}:
            return (("spacer", 1.0),)
        return (("mud", 1.0),)

    def _provider(time_s: float) -> AnnulusInletState:
        # 先查询套管内 1D 鞋口出流，再把流体名称映射为环空三相分数。
        pipe_exit_state = casing_solver.pipe_exit_state_at(casing_result, time_s)
        fluid_name = pipe_exit_state.phase_fractions[0][0] if pipe_exit_state.phase_fractions else ""
        mapped_pipe_exit = PipeExitState(
            time_s=pipe_exit_state.time_s,
            flow_rate_m3_s=pipe_exit_state.flow_rate_m3_s,
            stage_name=pipe_exit_state.stage_name,
            phase_fractions=_phase_fractions_for_fluid(fluid_name),
        )
        return pipe_exit_to_annulus_inlet(mapped_pipe_exit)

    return _provider
