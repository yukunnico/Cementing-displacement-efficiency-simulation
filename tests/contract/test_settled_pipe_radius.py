# -*- coding: utf-8 -*-
"""T5：停泵沉降的 τ_c 必须用本井内径（R = liner_id/2），不得用硬编码 0.05 m。

背景：`_settled_exit_fluid_name_enhanced` 内曾以无参形式调用
`_effective_pipe_radius_m()` ⇒ 半径恒取硬编码兜底 0.05 m。τ_c = Δρ·g·R
对半径线性，故 R 偏低即按同比例压低 τ_c。八井实测有**两个口径**，量级不同
（均为"0.05 m 相对真实半径偏低"，分母取真实半径）：

  · 相对 `liner_id_mm/2`（八井 0.05397–0.05558 m）：偏低 7.4–10.0%；
  · 相对模型自己的截面积口径 R_A = sqrt(pipe_cross_section_m2 / π)
    （八井 0.0581–0.0644 m；`pipe_cross_section_m2` 八井全部来自
    `shoe_lag_volume_m3 / shoe_md_m`）：偏低 13.9–22.4%
    （等价地，0.05 m 比 R_A 小 16–29%）。

两个口径彼此不一致属**已登记的独立口径裁定**（需用户先定），本任务不"修正"它，
只保证把**本井 well_spec** 送到停泵沉降路径，使 R 不再是无参兜底。

屈服抑制的高估**并非"系统性"**：`yield_ratio` 被 `min(τ_y/τ_c, 1.0)` 截顶、
`total_suppression` 被 `max(…, 0.99)` 吸收，故只有**未饱和**（τ_y/τ_c < 1 且
总抑制未触顶）的组合才出现位移；饱和组合下修复前后逐位相同。

零位移语义：该沉降分支只在停泵（flow≈0）+ 密度倒置 + 有屈服应力时进入，
生产 2D 时间窗 [0, cement_end_time_s] 内不可达（静态判据与行号见 task-4 报告
"修复轮 1"节）。⚠️ 逐位锚两文件（`test_default_path_bitwise_anchor.py` /
`test_anchor_integrity.py`）虽然全绿，但**对本分支零判别力**：锚把 inlet provider
接在一个**从未 run()** 的 `CasingFlowSolver` 上，其 id 表全空 ⇒ `delta_rho ≡ 0.0`、
`well_spec ≡ None`，即便四道闸门全开，本修复在锚里也是 no-op。锚绿是**结构性先验**，
不含本修复信息，只作"未把别处改坏"的回归护栏。
"""

import pytest

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


class _TauCriticalProbe(float):
    """浮点子类：在 `(Δρ·g) * R` 连乘处捕获**字面** τ_c（不靠测试端重算）。

    源码是 `tau_critical = delta_rho * self.g_constant * pipe_radius_m`，
    左结合先算出 `Δρ·g`（普通 float），再乘本 helper 返回的 R。因 float 子类的
    **反射方法优先**，`_TauCriticalProbe.__rmul__` 会收到 `Δρ·g` 并被调用，
    从而拿到求解器真正算出的 τ_c。
    """

    captured: float | None = None

    def __mul__(self, other: float) -> float:
        result = float(self) * other
        _TauCriticalProbe.captured = result
        return result

    __rmul__ = __mul__


class _RadiusSpy(CasingFlowSolver):
    """记录 `_effective_pipe_radius_m` 每次收到的实参与实际返回的半径。"""

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.radius_well_specs: list[object] = []
        self.radius_values: list[float] = []

    def _effective_pipe_radius_m(self, well_spec=None) -> float:
        self.radius_well_specs.append(well_spec)
        radius = super()._effective_pipe_radius_m(well_spec)
        self.radius_values.append(radius)
        return _TauCriticalProbe(radius)


def test_pipe_exit_state_at_hands_own_well_spec_to_settled_path():
    """行为守卫：停泵沉降必须收到本井 well_spec，且 τ_c 由本井半径算出。

    停泵窗（180–240 s）内 flow=0、鞋口相为尾浆（ρ1900 > 初始泥浆 ρ1200，
    Δρ=700）、尾浆带屈服应力 ⇒ 进入沉降分支的屈服抑制段。修复前
    `pipe_exit_state_at` 无法提供 well_spec，半径恒取兜底 0.05 m；修复后
    应收到 run() 那口井，τ_c = Δρ·g·R 随之抬到本井口径。

    本测试同时承担原先"源码字符串检查"的职责（`inspect.getsource` 断言已删）：
    `self._effective_pipe_radius_m(None)` / `( )` / `getattr(...)()` 这类
    绕过写法的语义都等价于缺陷本身（半径恒 0.05），会在下面的实参断言与
    τ_c 数值断言上同时变红。
    """

    well = _well()
    solver = _RadiusSpy(enable_gravity=True, enable_axial_dispersion=False)
    result = solver.run(well, _fluids(), _schedule())

    solver.radius_well_specs.clear()
    solver.radius_values.clear()
    _TauCriticalProbe.captured = None
    state = solver.pipe_exit_state_at(result, 210.0)

    assert state.flow_rate_m3_s < 1.0e-9, "210 s 应为停泵（零排量）时刻"
    assert solver.radius_well_specs, "pipe_exit_state_at 未进入停泵沉降分支"

    # ① 收到的实参：非 None，且就是本井（不是别井、不是空壳）
    received = solver.radius_well_specs[-1]
    assert received is well, "停泵沉降路径未收到 run() 的本井 well_spec（T5）"
    assert received.liner_id_mm == 112.84

    # ② 真正喂进 τ_c 的半径 = 本井 liner_id/2，而非硬编码 0.05
    assert solver.radius_values[-1] == 112.84 / 2000.0, (
        f"停泵沉降取的半径是 {solver.radius_values[-1]!r} m，不是本井内径半径")

    # ③ 字面 τ_c 确实变了（spy 在连乘处捕获，非测试端重算）
    delta_rho = 1900.0 - 1200.0
    tau_c_fixed = delta_rho * solver.g_constant * (112.84 / 2000.0)
    tau_c_defect = delta_rho * solver.g_constant * 0.05
    assert _TauCriticalProbe.captured is not None, (
        "未捕获到 τ_c = Δρ·g·R 的连乘（沉降屈服段未走到半径消费处）")
    assert _TauCriticalProbe.captured == pytest.approx(tau_c_fixed, rel=1e-12), (
        f"τ_c 实测 {_TauCriticalProbe.captured} Pa ≠ 本井口径 {tau_c_fixed} Pa")
    assert _TauCriticalProbe.captured > tau_c_defect, (
        f"τ_c 未相对缺陷口径（0.05 m）抬高：实测 {_TauCriticalProbe.captured} Pa，"
        f"缺陷口径 {tau_c_defect} Pa")
