# -*- coding: utf-8 -*-
"""Phase 4a+4c 契约测试：1D casing 逐深温度求值点 + 1D/2D 同源化（2026-10-07）。

权威 = docs/superpowers/specs/2026-10-07-phase4-depthwise-temperature-design.md
§1 项1 / §2「4a+4c 验收」；行号与波及口径 = 同日期测绘调研（项 1-4 / 12-13）。

三段式关2 + 入口口径 + B8 + 4c 同源 + 回退：

1. **关2① 不传场（T-off）⇒ 逐位 = HEAD**：短链 fronts/events、三入口与 B8
   两路直接值 == 改前 HEAD 捕获锚（与 test_temperature_casing_flow 同源字面量）；
   `_phase_props(f, md, t)` 恒等、memo 空；`_dispersion_depthwise_records`
   恒空（诊断零痕迹）。
2. **关2② Constant(60) ⇒ 逐位**：md-mode 与标量回退路径同 memo 键 ⇒ 同一派生
   对象；"强制标量路径"（monkeypatch 抹掉 (md,t)）与 md-mode 的 run 输出逐位
   相等；fluid_at 吃到的 T 恒 60.0。
3. **关2③ Geothermal ⇒ 消费点温度随 md 变**：鞋口 vs 域顶差 = GEO_GRAD×Δmd
   （手算断言）；`_phase_props(f, md, t)` == `fluid_at(f, T0+g·md, P)` 精确。
4. **三入口 (md,t) 口径**：重力=（鞋深, 到达时刻）终点单点；弥散=体积链反推
   沿程中点（常排量短链 ⇒ md=鞋深/2 手算 + 闭合残差≈0 落盘）；有效黏度
   （has_plug=False 才可达）与弥散同点。**两口径不共用**：弥散点 md 严格 < 鞋深。
5. **B8 两路包装**：legacy 重力、停泵沉降 raw τy 消费点走 `_phase_props`
   （记录 (md,t)，且 T-on 输出 ≠ T-off ⇒ 包装真实生效；T-off 逐位锚不变）。
6. **4c 同源**：Geothermal 短链各鞋口消费点 T_1d − field.T(shoe, t) ≡ 0；
   主消费 md 集合不含域顶（该值仅存在于回退标量 `_step_T_c`）⇒ 1D 不再用
   「域顶单点」，21.7°C 断裂语义（测绘项12）消除。
7. **回退路径**：md_m 缺省 ⇒ `_step_T_c`/`_step_P_mpa` 标量语义保留（无 run、
   无压力场、Constant P 同值三种情形）。
"""
from __future__ import annotations

import math

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.pressure_field import ConstantPressureField
from cemdisp.data.pumping_schedule import (
    PumpingSchedule,
    PumpingScheduleStep,
    PumpingStageEvent,
)
from cemdisp.data.rheology_vs_temperature import fluid_at
from cemdisp.data.temperature_field import (
    GEO_GRAD_C_PER_M,
    GEO_T0_C,
    ConstantTemperatureField,
    GeothermalTemperatureField,
)
from cemdisp.data.well_spec import WellSpec
from cemdisp.transport1d import casing_flow as casing_module
from cemdisp.transport1d.casing_flow import CasingFlowSolver

# --------------------------------------------------------------------------- #
# 短链装配（与 test_temperature_casing_flow 同源：top=1 / shoe=100 / area=0.01）
# --------------------------------------------------------------------------- #

# HEAD 逐位锚（改前 HEAD 探针捕获，同 test_temperature_casing_flow 断言值）
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


def _well() -> WellSpec:
    return WellSpec(
        well_name="P4a_depthwise_1d",
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


def _on_solver(**kwargs) -> CasingFlowSolver:
    kwargs.setdefault("enable_temperature_rheology", True)
    kwargs.setdefault("temperature_rheology_t_c", 60.0)
    return CasingFlowSolver(**kwargs)


def _spy_method(obj, name: str, log: list) -> None:
    """包装实例方法：按调用实参记录（bound method ⇒ args 不含 self）。"""
    original = getattr(obj, name)

    def spy(*args, **kwargs):
        log.append(args)
        return original(*args, **kwargs)

    setattr(obj, name, spy)


def _spy_fluid_at(log: list):
    """包装模块级 `fluid_at`：记录 (fluid, T, P) 后委托原函数，返回原函数。"""
    original = casing_module.fluid_at

    def spy(fluid, temperature_c, pressure_mpa=None, **kwargs):
        log.append((fluid, temperature_c, pressure_mpa))
        return original(fluid, temperature_c, pressure_mpa, **kwargs)

    casing_module.fluid_at = spy
    return original


def _times(result):
    return (
        tuple(f.time_s for f in result.fronts),
        tuple(e.time_s for e in result.shoe_timeline.events),
        result.pumping_end_time_s,
        result.cement_end_time_s,
        result.pipe_cross_section_m2,
    )


# --------------------------------------------------------------------------- #
# 关2①：不传场 ⇒ 逐位 = HEAD（T-off 恒等，与 md_m 取值无关）
# --------------------------------------------------------------------------- #
class TestOffBitwiseHead:
    def test_short_chain_bitwise_head_anchor(self) -> None:
        solver = CasingFlowSolver()
        res = solver.run(_well(), _fluids(), _schedule())
        assert tuple(f.time_s for f in res.fronts) == _FRONT_TIMES
        assert tuple(e.time_s for e in res.shoe_timeline.events) == _EVENT_TIMES

    def test_entries_and_b8_bitwise_head_anchor(self) -> None:
        solver = CasingFlowSolver()
        fluids = _fluids()
        tail = _tail()
        # T-off 即使给出 (md,t) 也恒等（回退红线：md 参数不被消费）
        assert solver._phase_props(tail, 55.0, 3.0) is tail
        assert solver._phase_memo == {}
        assert solver._compute_dispersion_coefficient(0.05, tail, 0.5) == _D_TAIL_U05
        assert solver._compute_dispersion_coefficient(0.05, tail, 0.5, 50.0, 10.0) == _D_TAIL_U05
        assert solver._effective_viscosity(tail, 0.5, 0.05) == _MU_TAIL_U05
        assert solver._effective_viscosity(tail, 0.5, 0.05, 50.0, 10.0) == _MU_TAIL_U05
        assert solver._gravity_corrected_arrival_time(
            100.0, "尾浆", "隔离液", fluids, _well()) == _GRAVITY
        legacy = CasingFlowSolver(enable_buoyancy_physics=False,
                                  settling_velocity_factor=1.0)
        assert legacy._gravity_corrected_arrival_time(
            100.0, "尾浆", "隔离液", fluids, _well()) == _GRAVITY_LEGACY

    def test_off_leaves_no_depthwise_records(self) -> None:
        solver = CasingFlowSolver()
        solver.run(_well(), _fluids(), _schedule())
        assert solver._dispersion_depthwise_records == []
        assert solver._phase_memo == {}
        # T-off 显式传场同样零痕迹（场不被读）
        off2 = CasingFlowSolver()
        off2.run(_well(), _fluids(), _schedule(),
                 temperature_field=GeothermalTemperatureField())
        assert off2._temperature_field is None
        assert off2._step_T_c is None
        assert off2._dispersion_depthwise_records == []


# --------------------------------------------------------------------------- #
# 关2②：Constant(60) 场 ⇒ md-mode 与标量回退路径逐位
# --------------------------------------------------------------------------- #
class TestConstantFieldBitwise:
    def test_memo_key_equal_returns_same_derived_object(self) -> None:
        solver = _on_solver()
        tail = _tail()
        scalar = solver._phase_props(tail)                      # 回退：_step_T_c=60
        depthwise = solver._phase_props(tail, 7868.0, 1234.5)   # Constant ⇒ 同 60.0
        assert depthwise is scalar  # 同键 ⇒ 同对象（下游算术序列不变 ⇒ 逐位构造性证明）

    def test_constant_run_vs_no_field_run_bitwise(self) -> None:
        a = _on_solver()
        ra = a.run(_well(), _fluids(), _schedule())  # 无场 ⇒ 回退 Constant(60)
        b = _on_solver()
        rb = b.run(_well(), _fluids(), _schedule(),
                   temperature_field=ConstantTemperatureField(60.0))
        assert _times(ra) == _times(rb)

    def test_forced_scalar_path_bitwise_vs_depthwise_run(self) -> None:
        """抹掉 (md,t)（= 4a 前纯标量消费）与 md-mode run 输出逐位相等。"""
        r_ref = _on_solver().run(_well(), _fluids(), _schedule(),
                                 temperature_field=ConstantTemperatureField(60.0))
        forced = _on_solver()
        original = forced._phase_props
        forced._phase_props = lambda fluid, md_m=None, t_s=None: original(fluid, None, None)
        r_forced = forced.run(_well(), _fluids(), _schedule(),
                              temperature_field=ConstantTemperatureField(60.0))
        assert _times(r_forced) == _times(r_ref)

    def test_fluid_at_t_always_constant_value(self) -> None:
        solver = _on_solver()
        log = []
        original = _spy_fluid_at(log)
        try:
            solver.run(_well(), _fluids(), _schedule(),
                       temperature_field=ConstantTemperatureField(60.0))
        finally:
            casing_module.fluid_at = original
        assert log, "T-on 短链必有派生消费"
        assert all(t == 60.0 for _, t, _ in log)


# --------------------------------------------------------------------------- #
# 关2③ + 4c：Geothermal 静温场——消费点温度随 md 变、与场求值同源
# --------------------------------------------------------------------------- #
class TestGeothermalDepthwiseResponse:
    def test_consumption_temperature_varies_with_md_handcomputed(self) -> None:
        solver = _on_solver()
        solver.run(_well(), _fluids(), _schedule(),
                   temperature_field=GeothermalTemperatureField())
        tail = _tail()
        log = []
        original = _spy_fluid_at(log)
        try:
            solver._phase_memo.clear()
            solver._phase_props(tail, 1.0, 0.0)     # 域顶
            solver._phase_memo.clear()
            solver._phase_props(tail, 100.0, 0.0)   # 鞋口
        finally:
            casing_module.fluid_at = original
        assert len(log) == 2
        t_top, t_shoe = log[0][1], log[1][1]
        # 手算断言：鞋口 vs 域顶求值差 = GEO_GRAD×Δmd（Δmd=99）
        assert abs((t_shoe - t_top) - GEO_GRAD_C_PER_M * 99.0) < 1e-12
        assert abs(t_top - (GEO_T0_C + GEO_GRAD_C_PER_M * 1.0)) < 1e-12

    def test_phase_props_equals_manual_fluid_at_per_point(self) -> None:
        solver = _on_solver()
        field = GeothermalTemperatureField()
        solver._temperature_field = field
        tail = _tail()
        for md, t in ((50.0, 0.0), (100.0, 60.5), (78.6, 120.0)):
            solver._phase_memo.clear()
            got = solver._phase_props(tail, md, t)
            expected = fluid_at(tail, field.T(md, t), None)
            assert got == expected  # 求值点=场对象本身（4c 同源的构造性表述）

    def test_shoe_consumption_same_source_diff_zero_and_no_domaintop(self) -> None:
        """4c 操作化：鞋口消费点 T_1d − field.T(shoe,t) ≡ 0；主消费不含域顶。"""
        solver = _on_solver()
        field = GeothermalTemperatureField()
        calls = []
        _spy_method(solver, "_phase_props", calls)
        solver.run(_well(), _fluids(), _schedule(), temperature_field=field)
        assert calls, "Geothermal 短链应有物性求值"
        # ① 鞋口消费点（重力/B8 口径）：1D 实际吃到的温度 vs 同一场对象在
        # (shoe_md, t) 的求值 ⇒ 差恒 0（同源操作化定义）
        shoe_points = sorted({(100.0, t) for _, md, t in calls if md == 100.0})
        assert shoe_points, "应有鞋口单点消费（重力入口=现成鞋深+到达时刻）"
        tail = _tail()
        for md, t in shoe_points:
            t_q = 0.0 if t is None else t
            solver._phase_memo.clear()
            log = []
            original = _spy_fluid_at(log)
            try:
                solver._phase_props(tail, md, t_q)
            finally:
                casing_module.fluid_at = original
            assert log and log[0][1] - field.T(md, t_q) == 0.0
        # ② 主消费 md 集合不含域顶（1.0）⇒「域顶单点标量」不再是消费口径
        used_md = {md for _, md, _ in calls if md is not None}
        assert used_md and 1.0 not in used_md
        # ③ 回退标量仍=域顶 t=0（保留但仅回退；与鞋口值差>1°C 证明断裂语义已移位）
        assert solver._step_T_c == field.T(1.0, 0.0)
        assert abs(solver._step_T_c - field.T(100.0, 0.0)) > 1.0

    def test_closure_residual_recorded_and_near_zero(self) -> None:
        """关5：等截面短链体积链闭合残差 ≈0（仅浮点），T-on 落盘非空。"""
        solver = _on_solver()
        solver.run(_well(), _fluids(), _schedule(),
                   temperature_field=GeothermalTemperatureField())
        recs = solver._dispersion_depthwise_records
        assert recs, "T-on 弥散事件应有 (md,t) 反推记录"
        for t_evt, md, t_rep, closure in recs:
            assert abs(closure) < 1e-9            # 管容/等效截面 ≡ 鞋深
            assert 0.0 <= md <= 100.0
            assert 0.0 <= t_rep <= t_evt


# --------------------------------------------------------------------------- #
# 三入口 (md,t) 口径（monkeypatch 记录 `_phase_props` 实参）
# --------------------------------------------------------------------------- #
class TestEntryCalibers:
    def test_gravity_shoe_endpoint_single_point(self) -> None:
        solver = _on_solver()
        calls = []
        _spy_method(solver, "_phase_props", calls)
        solver._gravity_corrected_arrival_time(
            250.0, "尾浆", "隔离液", _fluids(), _well())
        assert calls and calls[-1][1:] == (100.0, 250.0)
        # 测绘风险条：重力口径=鞋口终点单点（md ≡ 鞋深、t ≡ 到达时刻）
        calls.clear()
        # well_spec 缺失 ⇒ md None ⇒ 标量回退（不得凭空造"前缘深度"）
        solver._gravity_corrected_arrival_time(250.0, "尾浆", "隔离液", _fluids(), None)
        assert calls[-1][1:] == (None, 250.0)

    def test_dispersion_pathwise_midpoint_volume_chain(self) -> None:
        """弥散=体积链反推沿程中点：常排量短链 ⇒ md 手算 = 鞋深/2、t=到达−30s。

        独立运动学（由夹具日程推得，不复用实现内部）：1 m³/min 恒排量、
        管容=1.0 m³（=鞋深 100 m × 0.01 m²）⇒ 界面沿程中点恒在半深 50 m、
        中点时刻 = 到达时刻 − 行程时长/2 = t_arr − 30 s（t_arr>60 时）。
        """
        solver = _on_solver()
        calls = []
        _spy_method(solver, "_phase_props", calls)
        solver.run(_well(), _fluids(), _schedule(),
                   temperature_field=GeothermalTemperatureField())
        recs = solver._dispersion_depthwise_records
        assert recs, "T-on 弥散入口应有沿程中点反推记录"
        for t_arr, md, t_rep, _closure in recs:
            if t_arr >= 60.0:
                expect_t = 0.5 * ((t_arr - 60.0) + t_arr)   # t_dep=t_arr−60（同排量回推）
                expect_md = 50.0                            # 沿程中点=半深
            else:
                expect_t = 0.5 * t_arr                      # v_dep 越界钳 0 ⇒ t_dep=0
                expect_md = (expect_t / 60.0) / 0.01
            assert abs(md - expect_md) < 1e-6, (t_arr, md, expect_md)
            assert abs(t_rep - expect_t) < 1e-6, (t_arr, t_rep, expect_t)
        # 弥散消费点存在且全为沿程点（0 < md < 鞋深）
        disp = [c for c in calls if c[1] is not None and 0.0 < c[1] < 100.0]
        assert disp, "弥散入口应有沿程中点消费"
        # 口径不共用：没有任何"沿程点"取鞋深终点值（那是重力入口口径）
        assert all(0.0 < c[1] < 100.0 for c in disp)

    def test_effective_viscosity_same_point_as_dispersion(self) -> None:
        """有效黏度（has_plug=False 才可达）与弥散吃同一 (md,t)——死路仍改齐。"""
        solver = _on_solver(has_plug=False)  # 生产 8 井全 True ⇒ 仅探针/非胶塞井可达
        visc_calls = []
        phase_calls = []
        _spy_method(solver, "_effective_viscosity", visc_calls)
        _spy_method(solver, "_phase_props", phase_calls)
        solver.run(_well(), _fluids(), _schedule(),
                   temperature_field=GeothermalTemperatureField())
        assert visc_calls, "has_plug=False 时混浆增强应触达 `_effective_viscosity`"
        for _fluid, _u, _r, md, t in visc_calls:
            assert md is not None and t is not None
            # 签名一致性：透传的 (md,t) 与弥散入口同点（测绘项3 口径同源）
            assert any(c[1] == md and c[2] == t for c in phase_calls)
        # 生产口径复核：has_plug=True ⇒ 短路在前，本入口不被调用
        plugged = _on_solver(has_plug=True)
        v2 = []
        _spy_method(plugged, "_effective_viscosity", v2)
        plugged.run(_well(), _fluids(), _schedule(),
                    temperature_field=GeothermalTemperatureField())
        assert v2 == []

    def test_dispersion_direct_args_passthrough(self) -> None:
        solver = _on_solver()
        calls = []
        _spy_method(solver, "_phase_props", calls)
        solver._compute_dispersion_coefficient(0.05, _tail(), 0.5, 37.0, 11.0)
        assert calls[-1][1:] == (37.0, 11.0)
        calls.clear()
        solver._compute_dispersion_coefficient(0.05, _tail(), 0.5)  # 缺省 ⇒ 回退
        assert calls[-1][1:] == (None, None)


# --------------------------------------------------------------------------- #
# B8 两路 raw τy 包 `_phase_props`（口径一致性修补，测绘项4）
# --------------------------------------------------------------------------- #
class TestB8Wrapping:
    def test_legacy_gravity_wrapped(self) -> None:
        fluids = _fluids()
        off = CasingFlowSolver(enable_buoyancy_physics=False,
                               settling_velocity_factor=1.0)
        # T-off 逐位锚（HEAD 值）不因包装改变
        assert off._gravity_corrected_arrival_time(
            100.0, "尾浆", "隔离液", fluids, _well()) == _GRAVITY_LEGACY
        # T-on（构造层 90°C）：raw τy 消费点经 `_phase_props`，(md,t)=(鞋深,到达时刻)
        on = _on_solver(enable_buoyancy_physics=False,
                        settling_velocity_factor=1.0,
                        temperature_rheology_t_c=90.0)
        calls = []
        _spy_method(on, "_phase_props", calls)
        out = on._gravity_corrected_arrival_time(250.0, "尾浆", "隔离液", fluids, _well())
        assert calls and calls[-1][1:] == (100.0, 250.0)
        # 包装真实生效：派生 τy>0 进入屈服抑制 ⇒ 与 T-off 同参输出不同
        assert out != off._gravity_corrected_arrival_time(
            250.0, "尾浆", "隔离液", fluids, _well())

    def test_settled_exit_wrapped(self) -> None:
        solver = _on_solver()
        res = solver.run(_well(), _fluids(), _schedule(),
                         temperature_field=GeothermalTemperatureField())
        calls = []
        _spy_method(solver, "_phase_props", calls)
        # 停泵段 t=150：flow=0 ⇒ 走增强沉降；当前出流流体=尾浆（延迟体积 1.0 m³）
        state = solver.pipe_exit_state_at(res, 150.0)
        assert calls, "停泵沉降 raw τy 消费点应经 `_phase_props`"
        assert any(c[1] == 100.0 and c[2] == 150.0 for c in calls)  # 鞋深+查询时刻
        assert state.time_s == 150.0
        # T-off 同调用零痕迹：包装点在 T-off 恒等（返回原对象、memo 不写）；
        # 且停泵出流名与 HEAD 语义一致（锚 = T-off 短链的顶替后出流）。
        off = CasingFlowSolver()
        res_off = off.run(_well(), _fluids(), _schedule())
        state_off = off.pipe_exit_state_at(res_off, 150.0)
        assert off._phase_memo == {}
        assert state_off.phase_fractions[0][0] == _settled_head_name(off, res_off)


def _settled_head_name(solver, res) -> str:
    """HEAD 语义手工参照：T-off 时 `_settled_exit_fluid_name_enhanced` 吃原始 τy。

    尾浆 τy=None ⇒ yield_suppression=0；gel_factor=0.95*(1−exp(−30/600))；
    v_eff=0.0015*700*(1−gel)；settled=cum(150)+v_eff*30*0.01 vs 管容 1.0。
    """
    gel = 0.95 * (1.0 - math.exp(-30.0 / 600.0))
    v_eff = 0.0015 * (1900.0 - 1200.0) * (1.0 - min(gel, 0.99))
    cum_150 = 2.0  # 隔离液1+尾浆1（停泵步 0 体积）
    settled = cum_150 + v_eff * 30.0 * res.pipe_cross_section_m2
    return "泥浆" if settled < 1.0 else "尾浆"


# --------------------------------------------------------------------------- #
# 回退路径：标量 `_step_T_c` / `_step_P_mpa` 语义保留
# --------------------------------------------------------------------------- #
class TestFallbackScalarPath:
    def test_no_run_uses_construction_scalar(self) -> None:
        solver = CasingFlowSolver(enable_temperature_rheology=True,
                                  temperature_rheology_t_c=80.0)
        tail = _tail()
        assert solver._phase_props(tail) == fluid_at(tail, 80.0)
        assert solver._step_T_c == 80.0

    def test_no_md_after_run_uses_step_scalar(self) -> None:
        solver = _on_solver()
        solver.run(_well(), _fluids(), _schedule(),
                   temperature_field=GeothermalTemperatureField())
        tail = _tail()
        # md 缺省 ⇒ `_step_T_c`（域顶 t=0，=改前口径）；P 无场 ⇒ None
        assert solver._phase_props(tail) == fluid_at(tail, solver._step_T_c, None)
        assert solver._phase_props(tail) is solver._phase_props(tail)  # memo

    def test_pressure_field_md_mode_and_fallback(self) -> None:
        pfield = ConstantPressureField(35.0)
        solver = _on_solver(pressure_field=pfield)
        field = GeothermalTemperatureField()
        solver._temperature_field = field
        tail = _tail()
        # md-mode：P 与 T 同模式现场求值（P(md,0)，静压与 t 无关）
        solver._phase_memo.clear()
        got = solver._phase_props(tail, 100.0, 9.0)
        assert got == fluid_at(tail, field.T(100.0, 9.0), 35.0)
        # md 缺省 ⇒ 标量 `_step_P_mpa`（未 run ⇒ None ⇒ 与关2 旧口径同键族）
        solver._phase_memo.clear()
        assert solver._phase_props(tail) == fluid_at(tail, solver._step_T_c, None)
        # 无压力场：md-mode 也只吃 T，P=None
        plain = _on_solver()
        plain._temperature_field = field
        assert plain._phase_props(tail, 100.0, 9.0) == fluid_at(tail, field.T(100.0, 9.0), None)


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
