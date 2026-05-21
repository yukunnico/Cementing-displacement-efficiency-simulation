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
from typing import overload

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.provenance import WellProvenance
from cemdisp.transport1d.casing_flow import CasingFlowResult, CasingFlowSolver
from cemdisp.transport1d.pipe_exit_state import PipeExitState
from cemdisp.transport1d.shoe_timeline import ShoeEventKind, ShoeTimeline


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


def _phase_fractions_for_fluid(
    fluid_name: str,
    fluids: tuple[FluidSpec, ...],
    *,
    split_cement_phases: bool = False,
) -> tuple[tuple[str, float], ...]:
    """根据流体角色映射为环空二维模型入口相分数。

    默认保持历史三相口径：领浆/中间浆/尾浆统一归为 cement。
    当 split_cement_phases=True 时，把领浆与中间浆并入 lead，相当于"前导水泥相"；
    尾浆单独进入 tail，相容于现有 2D 层的 lead/tail 两段水泥跟踪能力。
    """

    role_by_name: dict[str, FluidRole] = {fluid.name: fluid.role for fluid in fluids}
    role = role_by_name.get(fluid_name, FluidRole.MUD)
    if split_cement_phases and role in {FluidRole.LEAD, FluidRole.INTERMEDIATE}:
        return (("lead", 1.0),)
    if split_cement_phases and role == FluidRole.TAIL:
        return (("tail", 1.0),)
    if role in {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}:
        return (("cement", 1.0),)
    if role in {FluidRole.WASH, FluidRole.SPACER}:
        return (("spacer", 1.0),)
    return (("mud", 1.0),)


def _phase_fractions_from_state(
    pipe_exit_state: PipeExitState,
    fluids: tuple[FluidSpec, ...],
    *,
    split_cement_phases: bool = False,
) -> tuple[tuple[str, float], ...]:
    """将鞋口出流状态映射为环空相分数，支持多相共存。

    当管内轴向弥散开启时，PipeExitState 可能包含多个相的共存分数。
    本函数将每个流体名称映射为环空相后，按体积分数加权合并。
    """
    mapped: dict[str, float] = {}
    for fluid_name, frac in pipe_exit_state.phase_fractions:
        sub_fractions = _phase_fractions_for_fluid(fluid_name, fluids, split_cement_phases=split_cement_phases)
        for phase_name, phase_frac in sub_fractions:
            mapped[phase_name] = mapped.get(phase_name, 0.0) + frac * phase_frac
    return tuple(sorted(mapped.items()))


@overload
def build_coupled_annulus_inlet_provider(
    shoe_timeline: ShoeTimeline,
    provenance: WellProvenance,
    fluids: tuple[FluidSpec, ...],
    *,
    split_cement_phases: bool = False,
) -> Callable[[float], AnnulusInletState]: ...


@overload
def build_coupled_annulus_inlet_provider(
    casing_result: CasingFlowResult,
    casing_solver: CasingFlowSolver,
    fluids: tuple[FluidSpec, ...],
    *,
    split_cement_phases: bool = False,
) -> Callable[[float], AnnulusInletState]: ...


def build_coupled_annulus_inlet_provider(
    arg1,
    arg2,
    fluids,
    *,
    split_cement_phases: bool = False,
):
    """构建1D-2D耦合的环空入口边界提供器。

    支持两种调用方式：
    1. 新方式：build_coupled_annulus_inlet_provider(timeline, provenance, fluids)
    2. 旧方式：build_coupled_annulus_inlet_provider(casing_result, casing_solver, fluids)

    该桥接函数把套管鞋口出流状态转换为环空入口边界，
    不改变环空二维求解器接口，从而保持 1D 层可插拔。
    """

    if isinstance(arg1, ShoeTimeline):
        shoe_timeline = arg1
        provenance = arg2

        def _provider(time_s: float) -> AnnulusInletState:
            pipe_exit_state = shoe_timeline.at(time_s)
            # 使用多相映射支持管内轴向弥散后的多相共存状态
            mapped_fractions = _phase_fractions_from_state(
                pipe_exit_state, fluids, split_cement_phases=split_cement_phases
            )
            return AnnulusInletState(
                time_s=pipe_exit_state.time_s,
                flow_rate_m3_s=pipe_exit_state.flow_rate_m3_s,
                stage_name=pipe_exit_state.stage_name,
                phase_fractions=mapped_fractions,
            )

        return _provider

    # 旧方式：arg1=casing_result, arg2=casing_solver
    casing_result = arg1
    casing_solver = arg2
    role_by_name: dict[str, FluidRole] = {fluid.name: fluid.role for fluid in fluids}

    def _legacy_phase_fractions_for_fluid(fluid_name: str) -> tuple[tuple[str, float], ...]:
        role = role_by_name.get(fluid_name, FluidRole.MUD)
        if split_cement_phases and role in {FluidRole.LEAD, FluidRole.INTERMEDIATE}:
            return (("lead", 1.0),)
        if split_cement_phases and role == FluidRole.TAIL:
            return (("tail", 1.0),)
        if role in {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}:
            return (("cement", 1.0),)
        if role in {FluidRole.WASH, FluidRole.SPACER}:
            return (("spacer", 1.0),)
        return (("mud", 1.0),)

    def _legacy_provider(time_s: float) -> AnnulusInletState:
        pipe_exit_state = casing_solver.pipe_exit_state_at(casing_result, time_s)
        # 使用多相映射支持管内轴向弥散后的多相共存状态
        mapped_fractions = _phase_fractions_from_state(pipe_exit_state, fluids, split_cement_phases=split_cement_phases)
        mapped_pipe_exit = PipeExitState(
            time_s=pipe_exit_state.time_s,
            flow_rate_m3_s=pipe_exit_state.flow_rate_m3_s,
            stage_name=pipe_exit_state.stage_name,
            phase_fractions=mapped_fractions,
        )
        return pipe_exit_to_annulus_inlet(mapped_pipe_exit)

    return _legacy_provider


def build_sync_card(
    well_name: str,
    shoe_timeline: ShoeTimeline,
    provenance: WellProvenance,
) -> dict[str, object]:
    """构造单井同步画像卡。

    汇总鞋口时间轴事件统计与 provenance 中的代理提醒，
    用于快速判断本井边界同步的主要不确定性来源。
    """

    events = shoe_timeline.events
    sync_note = ""
    if provenance.sync.status != "field":
        sync_note = provenance.sync.note

    return {
        "井名": well_name,
        "鞋口同步口径": {
            "事件数": len(events),
            "首事件时间_s": events[0].time_s if events else None,
            "末事件时间_s": events[-1].time_s if events else None,
        },
        "代理提醒": sync_note,
    }