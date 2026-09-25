# -*- coding: utf-8 -*-
"""T5：停泵沉降的 τ_c 必须用本井内径（R = liner_id/2），不得用硬编码 0.05 m。

背景：`_settled_exit_fluid_name_enhanced` 内以无参形式调用
`_effective_pipe_radius_m()` ⇒ 半径恒取硬编码兜底 0.05 m，而八井真实内径
对应 0.0539–0.0556 m ⇒ 临界屈服应力 τ_c = Δρ·g·R 偏低 8–10%，屈服抑制被
系统性高估。修复：把本井 well_spec 传进沉降路径。

零位移语义：该沉降分支只在停泵（flow≈0）+ 密度倒置 + 有屈服应力时进入，
生产 2D 时间窗内不可达（见 task-4 报告"不可达判据"），故本修复不改默认路径
任何数值——由 `test_default_path_bitwise_anchor.py` + `test_anchor_integrity.py`
两锚守护。
"""

import inspect

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.pumping_schedule import (
    PumpingSchedule,
    PumpingScheduleStep,
    PumpingStageEvent,
)
from cemdisp.data.well_spec import WellSpec
from cemdisp.transport1d.casing_flow import CasingFlowSolver


class _WellStub:
    """只带 liner_id_mm 的最小 stub（`_effective_pipe_radius_m` 的唯一消费字段）。"""

    liner_id_mm = 111.16


def test_effective_radius_uses_well_spec():
    """半径助手：有 well_spec 用本井内径；无 well_spec 退回 0.05 m（向后兼容）。"""

    s = CasingFlowSolver(enable_gravity=True)
    assert s._effective_pipe_radius_m(_WellStub()) == 111.16 / 2000.0
    assert s._effective_pipe_radius_m(None) == 0.05


def test_settled_path_passes_well_spec():
    """源码守卫：沉降路径不得再以无参形式取半径。"""

    src = inspect.getsource(CasingFlowSolver._settled_exit_fluid_name_enhanced)
    assert "_effective_pipe_radius_m()" not in src, (
        "停泵沉降仍以无参调用取半径（恒 0.05 m，T5 缺陷）")


def _fluids() -> tuple[FluidSpec, ...]:
    """泥浆（初始管内流体）→ 隔离液 → 尾浆（重且带屈服应力）→ 顶替液。

    停泵时刻鞋口相 = 尾浆，必须是"重于初始泥浆且有屈服应力"的流体
    （delta_rho > 0 且 yield_stress_pa > 0 才会进入屈服抑制计算）。
    顶替液步只为把尾浆前缘推进到泵注段内（否则前缘到达时刻退化为泵注终点）。
    """

    return (
        FluidSpec(name="泥浆", role=FluidRole.MUD, density_kg_m3=1200.0, plastic_viscosity_pa_s=0.02),
        FluidSpec(name="隔离液", role=FluidRole.SPACER, density_kg_m3=1100.0, plastic_viscosity_pa_s=0.01),
        FluidSpec(
            name="尾浆",
            role=FluidRole.TAIL,
            density_kg_m3=1900.0,
            plastic_viscosity_pa_s=0.08,
            yield_stress_pa=50.0,
        ),
        FluidSpec(name="顶替液", role=FluidRole.DISPLACEMENT, density_kg_m3=1100.0, plastic_viscosity_pa_s=0.02),
    )


def _well() -> WellSpec:
    """鞋深 100 m、管内容积 1 m³ 的单径测试井（内径 ≈ 112.8 mm）。"""

    return WellSpec(
        well_name="测试井",
        top_md_m=1.0,
        bottom_md_m=120.0,
        shoe_md_m=100.0,
        liner_id_mm=112.84,
    )


def _schedule() -> PumpingSchedule:
    """隔离液 0.5 → 尾浆 1.5 → 顶替液 1.0 m³ → 停泵候凝 60 s。

    管内容积 1 m³：前缘按"该步出发体积 + 管容"到达，故停泵（180–240 s）时
    鞋口相 = 尾浆（重、带屈服应力），满足 delta_rho>0 且 shutdown_duration>0。
    """

    def _step(name, fluid, volume, rate, t0, t1, tag):
        return PumpingScheduleStep(
            step_name=name,
            fluid_name=fluid,
            volume_m3=volume,
            rate_m3_min=rate,
            start_time_s=t0,
            end_time_s=t1,
            event_tag=tag,
        )

    return PumpingSchedule(
        steps=(
            _step("注入隔离液", "隔离液", 0.5, 1.0, 0.0, 30.0, PumpingStageEvent.INJECT_SPACER),
            _step("注入尾浆", "尾浆", 1.5, 1.0, 30.0, 120.0, PumpingStageEvent.INJECT_CEMENT),
            _step("注入顶替液", "顶替液", 1.0, 1.0, 120.0, 180.0, PumpingStageEvent.INJECT_DISPLACEMENT),
            _step("停泵候凝", "顶替液", 0.0, 0.0, 180.0, 240.0, PumpingStageEvent.SHUTDOWN),
        )
    )


class _RadiusSpy(CasingFlowSolver):
    """记录 `_effective_pipe_radius_m` 每次收到的 well_spec 实参。"""

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.radius_well_specs: list[object] = []

    def _effective_pipe_radius_m(self, well_spec=None) -> float:
        self.radius_well_specs.append(well_spec)
        return super()._effective_pipe_radius_m(well_spec)


def test_pipe_exit_state_at_hands_own_well_spec_to_settled_path():
    """行为守卫：run() 后的 result，其停泵沉降查询必须收到本井 well_spec。

    停泵窗（180–240 s）内 flow=0、鞋口相为尾浆（ρ1900 > 初始泥浆 ρ1200）、
    尾浆带屈服应力 ⇒ 进入沉降分支的屈服抑制段。修复前 `pipe_exit_state_at`
    无法提供 well_spec，半径助手收到 None；修复后应收到 run() 时那口井。
    """

    well = _well()
    solver = _RadiusSpy(enable_gravity=True, enable_axial_dispersion=False)
    result = solver.run(well, _fluids(), _schedule())

    solver.radius_well_specs.clear()
    state = solver.pipe_exit_state_at(result, 210.0)
    assert state.flow_rate_m3_s < 1.0e-9, "210 s 应为停泵（零排量）时刻"
    assert solver.radius_well_specs, "pipe_exit_state_at 未进入停泵沉降分支"
    assert solver.radius_well_specs[-1] is well, (
        "停泵沉降路径未收到 run() 的本井 well_spec（T5）")
    assert solver._effective_pipe_radius_m(well) == 112.84 / 2000.0
