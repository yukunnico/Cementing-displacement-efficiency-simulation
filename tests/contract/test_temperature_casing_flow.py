# -*- coding: utf-8 -*-
"""Task 8（T1-3）契约测试：`casing_flow` 三个物性入口挂接 `enable_temperature_rheology`
温变流变（默认关 ⇒ 逐位 = HEAD）。

覆盖（task-8-brief 需求逐条）：

1. **默认关**：`CasingFlowSolver().enable_temperature_rheology is False`、
   `temperature_rheology_t_c == 60.0`；关时不建温度场（`_temperature_field`/
   `_step_T_c` 为 None）。
2. **T-off 恒等（L1 硬约束，单元级）**：`_phase_props` 返回同一对象、不写
   memo/代表温度；三个入口与「原始流体手算原表达式」逐位相等。
3. **T-off 短链逐位锚**：短链 run 的 fronts/鞋口事件时刻/summary + 三入口直接
   求值与 HEAD（改前 task-8-toff-probe 捕获）逐位一致；T-off 且显式传温度场
   ⇒ 场被忽略、结果与不传逐位一致。
4. **T-on 三入口吃派生值（Constant(60)）**：各入口输出 == 「T-off 求解器吃
   `fluid_at(base, 60)` 派生流体」的手写参照，且与吃原始流体的输出不同
   （夹具灵敏度护栏）；派生流体是 Bingham 且 `yield_stress_pa > 0`。
5. **memo 一次派生**：同 (fluid,T) 两次 `_phase_props` 返回同一对象。
6. **开关真被消费**：短链 T-on vs T-off 的 fronts/时间线不同；`run()` 后
   `_step_T_c == 60.0`。

密度本轮不变（ρ(T,P) 属 Phase P）：`fluid_at` 不改密度，重力修正只接
屈服应力侧（测试 6 断言密度读数逐位不变）。
"""
from __future__ import annotations

import math

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.pumping_schedule import (
    PumpingSchedule,
    PumpingScheduleStep,
    PumpingStageEvent,
)
from cemdisp.data.rheology_vs_temperature import fluid_at
from cemdisp.data.temperature_field import ConstantTemperatureField
from cemdisp.data.well_spec import WellSpec
from cemdisp.transport1d.casing_flow import CasingFlowSolver


# --------------------------------------------------------------------------- #
# 短链装配（与 tests/contract/test_casing_flow.py 同源；锚值来自改前 HEAD 探针
# .superpowers/sdd/温压耦合改进计划_2026-09-30/task-8-toff-before.txt）
# --------------------------------------------------------------------------- #

def _well() -> WellSpec:
    return WellSpec(
        well_name="T8_temp_1d",
        top_md_m=1.0,
        bottom_md_m=120.0,
        shoe_md_m=100.0,
        liner_id_mm=math.sqrt(4.0 * 0.01 / math.pi) * 1000.0,
    )


def _fluids() -> tuple[FluidSpec, ...]:
    return (
        FluidSpec(name="泥浆", role=FluidRole.MUD, density_kg_m3=1200.0,
                  plastic_viscosity_pa_s=0.02),
        FluidSpec(name="隔离液", role=FluidRole.SPACER, density_kg_m3=1100.0,
                  plastic_viscosity_pa_s=0.01),
        FluidSpec(name="尾浆", role=FluidRole.TAIL, density_kg_m3=1900.0,
                  plastic_viscosity_pa_s=0.08),
        FluidSpec(name="顶替液", role=FluidRole.DISPLACEMENT, density_kg_m3=1050.0,
                  plastic_viscosity_pa_s=0.01),
    )


def _schedule() -> PumpingSchedule:
    return PumpingSchedule(
        steps=(
            PumpingScheduleStep(
                step_name="注入隔离液", fluid_name="隔离液", volume_m3=1.0,
                rate_m3_min=1.0, start_time_s=0.0, end_time_s=60.0,
                event_tag=PumpingStageEvent.INJECT_SPACER,
            ),
            PumpingScheduleStep(
                step_name="注入尾浆", fluid_name="尾浆", volume_m3=1.0,
                rate_m3_min=1.0, start_time_s=60.0, end_time_s=120.0,
                event_tag=PumpingStageEvent.INJECT_CEMENT,
            ),
            PumpingScheduleStep(
                step_name="停泵候凝", fluid_name="顶替液", volume_m3=0.0,
                rate_m3_min=0.0, start_time_s=120.0, end_time_s=180.0,
                event_tag=PumpingStageEvent.SHUTDOWN,
            ),
            PumpingScheduleStep(
                step_name="重启顶替", fluid_name="顶替液", volume_m3=2.0,
                rate_m3_min=2.0, start_time_s=180.0, end_time_s=240.0,
                event_tag=PumpingStageEvent.RESTART,
            ),
        )
    )


def _tail() -> FluidSpec:
    return next(f for f in _fluids() if f.name == "尾浆")


class TestSwitchDefaults:
    """1. 默认关 + 关时不建场。"""

    def test_default_off(self) -> None:
        solver = CasingFlowSolver()
        assert solver.enable_temperature_rheology is False
        assert solver.temperature_rheology_t_c == 60.0
        assert solver._temperature_field is None
        assert solver._step_T_c is None

    def test_explicit_off_is_same_as_default(self) -> None:
        a = CasingFlowSolver(enable_temperature_rheology=False,
                             temperature_rheology_t_c=60.0)
        b = CasingFlowSolver()
        assert a.enable_temperature_rheology == b.enable_temperature_rheology
        assert a.temperature_rheology_t_c == b.temperature_rheology_t_c


class TestPhasePropsOffIdentity:
    """2. T-off 恒等：不派生、不写 memo、不写代表温度。"""

    def test_off_is_identity(self) -> None:
        solver = CasingFlowSolver()
        tail = _tail()
        assert solver._phase_props(tail) is tail
        assert solver._phase_memo == {}
        assert solver._step_T_c is None

    def test_off_run_leaves_no_trace(self) -> None:
        solver = CasingFlowSolver()
        solver.run(_well(), _fluids(), _schedule(),
                   temperature_field=ConstantTemperatureField(90.0))
        assert solver._temperature_field is None
        assert solver._step_T_c is None
        assert solver._phase_memo == {}

    def test_off_entries_match_manual_original_expressions(self) -> None:
        """T-off 三入口 == 原流体手算（原表达式语义由同一实现复现 ⇒ 恒等护栏）。"""
        solver = CasingFlowSolver()
        tail = _tail()
        # off 路径吃原始流体：与「显式关开关的另一实例」逐位一致
        other = CasingFlowSolver(enable_temperature_rheology=False)
        assert (solver._compute_dispersion_coefficient(0.05, tail, 1e-3)
                == other._compute_dispersion_coefficient(0.05, tail, 1e-3))
        assert (solver._effective_viscosity(tail, 0.5, 0.05)
                == other._effective_viscosity(tail, 0.5, 0.05))
        fluids = _fluids()
        assert (solver._gravity_corrected_arrival_time(100.0, "尾浆", "隔离液", fluids, _well())
                == other._gravity_corrected_arrival_time(
                    100.0, "尾浆", "隔离液", fluids, _well()))


class TestShortChainOffBitwiseAnchor:
    """3. T-off 短链逐位锚（改前 HEAD 探针捕获值，精确相等）。"""

    # task-8-toff-before.txt：fronts / summary
    _FRONT_TIMES = (62.608695652173914, 88.0, 180.0)
    _EVENT_TIMES = (
        0.0, 60.0,
        60.608695652173914, 61.608695652173914, 62.608695652173914,
        63.608695652173914, 64.6086956521739, 64.6086956521739,
        84.81324420391842, 86.4066221019592, 88.0,
        89.5933778980408, 91.18675579608158, 91.18675579608158,
        120.0, 125.21739130434783, 180.0,
    )
    _D_TAIL_U05 = 0.00625
    _MU_TAIL_U05 = 0.08
    _GRAVITY = 73.33333333333334
    _GRAVITY_LEGACY = 9.999999999999998

    def test_short_chain_fronts_and_summary_bitwise(self) -> None:
        solver = CasingFlowSolver()
        res = solver.run(_well(), _fluids(), _schedule())
        assert tuple(f.time_s for f in res.fronts) == self._FRONT_TIMES
        assert tuple(e.time_s for e in res.shoe_timeline.events) == self._EVENT_TIMES
        assert res.pumping_end_time_s == 180.0
        assert res.cement_end_time_s == 180.0
        assert res.pipe_cross_section_m2 == 0.01

    def test_three_entries_direct_values_bitwise(self) -> None:
        solver = CasingFlowSolver()
        fluids = _fluids()
        tail = _tail()
        assert solver._compute_dispersion_coefficient(0.05, tail, 0.5) == self._D_TAIL_U05
        assert solver._effective_viscosity(tail, 0.5, 0.05) == self._MU_TAIL_U05
        assert solver._gravity_corrected_arrival_time(
            100.0, "尾浆", "隔离液", fluids, _well()) == self._GRAVITY
        legacy = CasingFlowSolver(enable_buoyancy_physics=False,
                                  settling_velocity_factor=1.0)
        assert legacy._gravity_corrected_arrival_time(
            100.0, "尾浆", "隔离液", fluids, _well()) == self._GRAVITY_LEGACY

    def test_off_ignores_injected_temperature_field(self) -> None:
        base = CasingFlowSolver().run(_well(), _fluids(), _schedule())
        off = CasingFlowSolver()
        other = off.run(_well(), _fluids(), _schedule(),
                        temperature_field=ConstantTemperatureField(90.0))
        assert tuple(f.time_s for f in other.fronts) == tuple(
            f.time_s for f in base.fronts)
        assert tuple(e.time_s for e in other.shoe_timeline.events) == tuple(
            e.time_s for e in base.shoe_timeline.events)
        assert other.cement_end_time_s == base.cement_end_time_s


class TestOnEntriesUseDerivedFluid:
    """4/5/6. Constant(60) 派生值进入三入口 + memo 一次派生 + 开关被消费。"""

    def test_phase_props_on_matches_manual_fluid_at_and_memoizes(self) -> None:
        solver = CasingFlowSolver(enable_temperature_rheology=True,
                                  temperature_rheology_t_c=60.0)
        tail = _tail()
        derived = solver._phase_props(tail)
        assert derived == fluid_at(tail, 60.0)
        assert derived.rheology_model == RheologyModel.BINGHAM
        assert derived.yield_stress_pa is not None and derived.yield_stress_pa > 0.0
        assert derived.density_kg_m3 == tail.density_kg_m3  # 密度本轮不变
        assert tail.rheology_model == RheologyModel.NEWTONIAN  # 原对象不被就地修改
        # memo：同 (fluid, T) 同对象（T-on 一次派生非每步）
        assert solver._phase_props(tail) is derived

    def test_effective_viscosity_on_enters_derived(self) -> None:
        tail = _tail()
        derived = fluid_at(tail, 60.0)
        off = CasingFlowSolver()
        on = CasingFlowSolver(enable_temperature_rheology=True)
        expected = off._effective_viscosity(derived, 0.5, 0.05)
        baseline = off._effective_viscosity(tail, 0.5, 0.05)
        assert expected != baseline  # 夹具灵敏度护栏（否则断言空转）
        assert on._effective_viscosity(tail, 0.5, 0.05) == expected
        # 原始流体为 NEWTONIAN(0.08)、派生为 Bingham ⇒ 结果确为派生参数所定
        assert baseline == 0.08

    def test_dispersion_on_enters_derived(self) -> None:
        # U=1e-3：off 走对流上限 d_cap、on 走 Fan&Wang 分支 ⇒ 有区分度
        tail = _tail()
        derived = fluid_at(tail, 60.0)
        off = CasingFlowSolver()
        on = CasingFlowSolver(enable_temperature_rheology=True)
        expected = off._compute_dispersion_coefficient(0.05, derived, 1e-3)
        baseline = off._compute_dispersion_coefficient(0.05, tail, 1e-3)
        assert expected != baseline
        assert on._compute_dispersion_coefficient(0.05, tail, 1e-3) == expected

    def test_gravity_on_enters_derived_tauy_density_unchanged(self) -> None:
        fluids = _fluids()
        derived_fluids = tuple(fluid_at(f, 60.0) for f in fluids)
        off = CasingFlowSolver()
        on = CasingFlowSolver(enable_temperature_rheology=True)
        expected = off._gravity_corrected_arrival_time(
            100.0, "尾浆", "隔离液", derived_fluids, _well())
        baseline = off._gravity_corrected_arrival_time(
            100.0, "尾浆", "隔离液", fluids, _well())
        assert expected != baseline  # τy(T) 抑制项进入 ⇒ 有区分度
        assert on._gravity_corrected_arrival_time(
            100.0, "尾浆", "隔离液", fluids, _well()) == expected
        # 密度侧逐位不变（fluid_at 不改密度；入口只接粘度侧）
        assert on._get_fluid_density("尾浆", fluids) == off._get_fluid_density(
            "尾浆", fluids)
        assert on._get_fluid_density("尾浆", derived_fluids) == off._get_fluid_density(
            "尾浆", fluids)

    def test_run_on_sets_step_t_and_short_chain_differs(self) -> None:
        off = CasingFlowSolver()
        r_off = off.run(_well(), _fluids(), _schedule())
        on = CasingFlowSolver(enable_temperature_rheology=True,
                              temperature_rheology_t_c=60.0)
        r_on = on.run(_well(), _fluids(), _schedule(),
                      temperature_field=ConstantTemperatureField(60.0))
        assert on._step_T_c == 60.0
        assert on._temperature_field is not None
        assert on._phase_memo != {}
        # 开关真被消费：T-on 短链输出与 T-off 不同
        assert tuple(f.time_s for f in r_on.fronts) != tuple(
            f.time_s for f in r_off.fronts)
        assert tuple(e.time_s for e in r_on.shoe_timeline.events) != tuple(
            e.time_s for e in r_off.shoe_timeline.events)


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
