"""
T1-10 混浆段增强分散 + T1-11 浮力滑移物理化 单元测试。

对应设计规格：docs/superpowers/specs/2026-08-07-1d-casing-mixing-buoyancy-design.md §7.1

覆盖项：
- T1-11：Δρ=0（等密度流体）→ 修正为零，到达时间不变
- T1-11：重驱轻（水泥驱泥浆）→ 到达时间缩短
- T1-11：轻驱重（泥浆驱水泥）→ 到达时间延长
- T1-11：τ_y→0 极限 → 回到无屈服抑制的浮力修正
- T1-11：enable_buoyancy_physics=False → 回退到旧经验乘子
- T1-10：At=0（等密度界面）→ 增强因子=1
- T1-10：重驱轻+高Re → 增强因子>1
- T1-10：has_plug=True → 增强因子=1
- T1-10：enable_mixing_enhancement=False → 不增强
- 辅助方法 _displaced_fluid_name 与被顶替流体接线
"""

import math
import unittest

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.pumping_schedule import (
    PumpingSchedule,
    PumpingScheduleStep,
    PumpingStageEvent,
)
from cemdisp.data.well_spec import DepthValuePoint, WellSpec
from cemdisp.transport1d.casing_flow import CasingFlowSolver
from cemdisp.transport1d.shoe_timeline import ShoeEventKind


def _fluids() -> tuple[FluidSpec, ...]:
    """测试用流体；名称与施工步骤保持一致。"""

    return (
        FluidSpec(name="泥浆", role=FluidRole.MUD, density_kg_m3=1200.0, plastic_viscosity_pa_s=0.02),
        FluidSpec(name="隔离液", role=FluidRole.SPACER, density_kg_m3=1100.0, plastic_viscosity_pa_s=0.01),
        FluidSpec(name="尾浆", role=FluidRole.TAIL, density_kg_m3=1900.0, plastic_viscosity_pa_s=0.08),
    )


def _simple_well() -> WellSpec:
    """鞋深 100 m、管内容积 1 m³（面积 0.01 m²）的垂直测试井，无井斜剖面。"""

    return WellSpec(
        well_name="测试井",
        top_md_m=1.0,
        bottom_md_m=120.0,
        shoe_md_m=100.0,
        liner_id_mm=math.sqrt(4.0 * 0.01 / math.pi) * 1000.0,
    )


def _inclined_well(inclination_deg: float) -> WellSpec:
    """带单一井斜角的测试井。"""

    return WellSpec(
        well_name="斜井测试井",
        top_md_m=1.0,
        bottom_md_m=120.0,
        shoe_md_m=100.0,
        liner_id_mm=math.sqrt(4.0 * 0.01 / math.pi) * 1000.0,
        inclination_profile=(DepthValuePoint(1.0, inclination_deg),),
    )


def _single_step_schedule() -> PumpingSchedule:
    """单步骤注水泥程序（水泥驱泥浆，重驱轻）。"""

    return PumpingSchedule(
        steps=(
            PumpingScheduleStep(
                step_name="注水泥",
                fluid_name="尾浆",
                volume_m3=20.0,
                rate_m3_min=2.0,
                event_tag=PumpingStageEvent.INJECT_CEMENT,
            ),
        )
    )


class TestBuoyancyPhysics(unittest.TestCase):
    """T1-11 浮力滑移物理化（Atwood 数公式）测试。"""

    def _solver(self, **kwargs) -> CasingFlowSolver:
        return CasingFlowSolver(**kwargs)

    def test_equal_density_no_correction(self) -> None:
        """Δρ=0（等密度流体）：At=0 → 修正为零，到达时间不变。"""

        solver = self._solver()
        # 等密度：泥浆(1200) 驱 等密度流体(1200)
        fluids = (
            FluidSpec(name="泥浆A", role=FluidRole.MUD, density_kg_m3=1200.0, plastic_viscosity_pa_s=0.02),
            FluidSpec(name="泥浆B", role=FluidRole.OTHER, density_kg_m3=1200.0, plastic_viscosity_pa_s=0.02),
        )
        result = solver._gravity_corrected_arrival_time(1000.0, "泥浆B", "泥浆A", fluids)
        self.assertAlmostEqual(result, 1000.0, places=9)

    def test_heavy_displaces_light_speeds_up(self) -> None:
        """重驱轻（水泥 1900 驱泥浆 1200）：到达时间缩短。"""

        solver = self._solver()
        fluids = _fluids()
        at = (1900.0 - 1200.0) / (1900.0 + 1200.0)  # = 700/3100
        result = solver._gravity_corrected_arrival_time(1000.0, "尾浆", "泥浆", fluids)
        self.assertAlmostEqual(result, 1000.0 * (1.0 - at), places=9)
        self.assertLess(result, 1000.0)

    def test_light_displaces_heavy_slows_down(self) -> None:
        """轻驱重（隔离液 1100 驱水泥 1900）：到达时间延长。"""

        solver = self._solver()
        fluids = _fluids()
        at = (1900.0 - 1100.0) / (1900.0 + 1100.0)  # = 800/3000
        result = solver._gravity_corrected_arrival_time(1000.0, "隔离液", "尾浆", fluids)
        self.assertAlmostEqual(result, 1000.0 * (1.0 + at), places=9)
        self.assertGreater(result, 1000.0)

    def test_inclination_projection_halves_correction(self) -> None:
        """井斜角 60° 投影：重力因子按 cos(60°)=0.5 折减。"""

        solver = self._solver()
        fluids = _fluids()
        well = _inclined_well(60.0)
        at = 700.0 / 3100.0
        result = solver._gravity_corrected_arrival_time(1000.0, "尾浆", "泥浆", fluids, well)
        self.assertAlmostEqual(result, 1000.0 * (1.0 - 0.5 * at), places=9)

    def test_yield_stress_zero_limit(self) -> None:
        """τ_y→0 极限：回到无屈服抑制的浮力修正。"""

        solver = self._solver()
        tail_newtonian = FluidSpec(name="尾浆", role=FluidRole.TAIL, density_kg_m3=1900.0, plastic_viscosity_pa_s=0.08)
        tail_bingham = FluidSpec(
            name="尾浆B", role=FluidRole.TAIL, density_kg_m3=1900.0,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=0.08, yield_stress_pa=0.0,
        )
        fluids = (FluidSpec(name="泥浆", role=FluidRole.MUD, density_kg_m3=1200.0, plastic_viscosity_pa_s=0.02), tail_newtonian)

        base = solver._gravity_corrected_arrival_time(1000.0, "尾浆", "泥浆", fluids)
        fluids_bingham = (fluids[0], tail_bingham)
        bingham_zero = solver._gravity_corrected_arrival_time(1000.0, "尾浆B", "泥浆", fluids_bingham)
        self.assertAlmostEqual(bingham_zero, base, places=9)

        # 正屈服应力应抑制浮力（修正幅度变小 → 到达时间更接近未修正值）
        tail_bingham_y = FluidSpec(
            name="尾浆BY", role=FluidRole.TAIL, density_kg_m3=1900.0,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=0.08, yield_stress_pa=50.0,
        )
        fluids_yield = (fluids[0], tail_bingham_y)
        with_yield = solver._gravity_corrected_arrival_time(1000.0, "尾浆BY", "泥浆", fluids_yield)
        self.assertGreater(with_yield, bingham_zero)
        self.assertLess(with_yield, 1000.0)

    def test_legacy_fallback_matches_old_formula(self) -> None:
        """enable_buoyancy_physics=False：回退到旧经验乘子公式。"""

        solver = self._solver(enable_buoyancy_physics=False, settling_velocity_factor=0.0015)
        fluids = _fluids()
        # 旧公式：gravity_factor = 0.0015 * (1900-1000)/1000 = 0.00135
        expected = 1000.0 * (1.0 - 0.0015 * 900.0 / 1000.0)
        result = solver._gravity_corrected_arrival_time(1000.0, "尾浆", "泥浆", fluids)
        self.assertAlmostEqual(result, expected, places=9)

        # 旧公式不依赖被顶替流体
        result_other = solver._gravity_corrected_arrival_time(1000.0, "尾浆", "隔离液", fluids)
        self.assertAlmostEqual(result_other, expected, places=9)

    def test_displaced_fluid_name_helper(self) -> None:
        """_displaced_fluid_name：首步回退初始泥浆，后续为上一步流体。"""

        fluids = _fluids()
        schedule = PumpingSchedule(
            steps=(
                PumpingScheduleStep(step_name="注隔离液", fluid_name="隔离液", volume_m3=1.0, rate_m3_min=1.0, event_tag=PumpingStageEvent.INJECT_SPACER),
                PumpingScheduleStep(step_name="注尾浆", fluid_name="尾浆", volume_m3=1.0, rate_m3_min=1.0, event_tag=PumpingStageEvent.INJECT_CEMENT),
            )
        )
        scheduled_steps = CasingFlowSolver()._build_scheduled_steps(schedule)
        self.assertEqual(CasingFlowSolver._displaced_fluid_name(scheduled_steps, 0, "泥浆"), "泥浆")
        self.assertEqual(CasingFlowSolver._displaced_fluid_name(scheduled_steps, 1, "泥浆"), "隔离液")


class TestMixingEnhancement(unittest.TestCase):
    """T1-10 混浆段增强分散测试。"""

    def _solver(self, **kwargs) -> CasingFlowSolver:
        return CasingFlowSolver(**kwargs)

    def _mud_cement(self) -> tuple[FluidSpec, FluidSpec]:
        mud = FluidSpec(name="泥浆", role=FluidRole.MUD, density_kg_m3=1200.0, plastic_viscosity_pa_s=0.02)
        cement = FluidSpec(name="尾浆", role=FluidRole.TAIL, density_kg_m3=1900.0, plastic_viscosity_pa_s=0.08)
        return mud, cement

    def test_equal_density_interface_factor_one(self) -> None:
        """At=0（等密度界面）：增强因子=1。"""

        solver = self._solver()
        a = FluidSpec(name="A", role=FluidRole.MUD, density_kg_m3=1200.0, plastic_viscosity_pa_s=0.02)
        b = FluidSpec(name="B", role=FluidRole.MUD, density_kg_m3=1200.0, plastic_viscosity_pa_s=0.02)
        factor = solver._interface_instability_factor(b, a, 0.05, 0.5)
        self.assertEqual(factor, 1.0)

    def test_heavy_displaces_light_high_re_enhances(self) -> None:
        """重驱轻 + 高 Re：增强因子>1（与规格公式一致）。"""

        solver = self._solver()
        mud, cement = self._mud_cement()
        R, U = 0.05, 0.5
        # 手算：mu_mean = sqrt(0.08*0.02)=0.04；Re = 1550*0.5*0.1/0.04 = 1937.5
        at = 700.0 / 3100.0
        Re = 1550.0 * U * 2.0 * R / math.sqrt(0.08 * 0.02)
        expected = min(1.0 + solver.mixing_enhancement_factor * at * math.sqrt(Re / 100.0), solver.max_mixing_enhancement)
        factor = solver._interface_instability_factor(cement, mud, R, U)
        self.assertAlmostEqual(factor, expected, places=9)
        self.assertGreater(factor, 1.0)
        self.assertLessEqual(factor, solver.max_mixing_enhancement)

    def test_has_plug_suppresses_enhancement(self) -> None:
        """has_plug=True：胶塞刮拭阻止混合 → 增强因子=1。"""

        solver = self._solver(has_plug=True)
        mud, cement = self._mud_cement()
        factor = solver._interface_instability_factor(cement, mud, 0.05, 0.5)
        self.assertEqual(factor, 1.0)

    def test_mixing_disabled_keeps_taylor_aris_only(self) -> None:
        """enable_mixing_enhancement=False：不注入增强因子（过渡带较窄）。"""

        fluids = _fluids()
        well = WellSpec(
            well_name="测试井",
            top_md_m=1.0,
            bottom_md_m=1100.0,
            shoe_md_m=1000.0,
            liner_id_mm=124.3,
        )

        def _span(solver: CasingFlowSolver) -> float:
            result = solver.run(well, fluids, _single_step_schedule())
            transition_times = [
                e.time_s for e in result.shoe_timeline.events
                if e.kind == ShoeEventKind.FRONT_ARRIVAL and len(e.phase_fractions) > 1
            ]
            return (max(transition_times) - min(transition_times)) if transition_times else 0.0

        off = _span(self._solver(enable_mixing_enhancement=False))
        on = _span(self._solver(enable_mixing_enhancement=True))
        self.assertGreater(on, off)  # 混浆增强 → 过渡带更宽

    def test_effective_viscosity_newtonian(self) -> None:
        """_effective_viscosity：牛顿流体直接返回塑性黏度。"""

        solver = self._solver()
        mud, _ = self._mud_cement()
        self.assertAlmostEqual(solver._effective_viscosity(mud, 0.5, 0.05), 0.02, places=9)
        self.assertAlmostEqual(solver._effective_viscosity(mud, 0.0, 0.05), 0.02, places=9)  # 零流速回退

    def test_effective_viscosity_bingham_shear_thinning(self) -> None:
        """_effective_viscosity：Bingham 表观黏度 = τ_y/γ + k_cons·γ^(n-1)。

        按 spec 4.4：Bingham 无 power-law 参数时 k_cons 兜底 0.01、n=1.0。
        """

        solver = self._solver()
        bingham = FluidSpec(
            name="B", role=FluidRole.TAIL, density_kg_m3=1900.0,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=0.05, yield_stress_pa=4.0,
        )
        R = 0.05
        for U in (0.1, 0.5, 1.0):
            shear_rate = 8.0 * U / (2.0 * R)
            expected = 4.0 / shear_rate + 0.01 * shear_rate ** 0.0
            self.assertAlmostEqual(solver._effective_viscosity(bingham, U, R), expected, places=9)

    def test_effective_viscosity_herschel_bulkley(self) -> None:
        """_effective_viscosity：HB 表观黏度 = τ_y/γ + k·γ^(n-1)。"""

        solver = self._solver()
        hb = FluidSpec(
            name="HB", role=FluidRole.TAIL, density_kg_m3=1800.0,
            rheology_model=RheologyModel.HERSCHEL_BULKLEY,
            yield_stress_pa=2.0, power_law_n=0.7, consistency_k=1.5,
        )
        U = 0.3
        R = 0.05
        shear_rate = 8.0 * U / (2.0 * R)
        expected = 2.0 / shear_rate + 1.5 * shear_rate ** (0.7 - 1.0)
        self.assertAlmostEqual(solver._effective_viscosity(hb, U, R), expected, places=9)


class TestRunIntegration(unittest.TestCase):
    """run() 调用点接线集成验证。"""

    def test_default_solver_exercises_atwood_correction_heavy_drives_light(self) -> None:
        """默认参数（buoyancy physics on）：重驱轻 → run() 前缘时间按 Atwood 缩短。"""

        fluids = _fluids()
        # 隔离液(1100) 驱 泥浆(1200)：轻驱重 → Atwood 减速修正
        at = (1200.0 - 1100.0) / (1200.0 + 1100.0)  # 100/2300
        schedule = PumpingSchedule(
            steps=(
                PumpingScheduleStep(step_name="注隔离液", fluid_name="隔离液", volume_m3=1.0, rate_m3_min=1.0, event_tag=PumpingStageEvent.INJECT_SPACER),
            )
        )
        result = CasingFlowSolver().run(_simple_well(), fluids, schedule)
        self.assertAlmostEqual(result.fronts[0].time_s, 60.0 * (1.0 + at), places=6)

    def test_legacy_solver_preserves_old_front_times(self) -> None:
        """enable_buoyancy_physics=False：恢复旧公式的前缘时间（向后兼容）。"""

        fluids = _fluids()
        # 隔离液: 旧公式 gravity_factor = 1.0 * (1100-1000)/1000 = 0.1（settling_velocity_factor=1.0）
        schedule = PumpingSchedule(
            steps=(
                PumpingScheduleStep(step_name="注隔离液", fluid_name="隔离液", volume_m3=1.0, rate_m3_min=1.0, event_tag=PumpingStageEvent.INJECT_SPACER),
            )
        )
        solver = CasingFlowSolver(enable_buoyancy_physics=False, settling_velocity_factor=1.0)
        result = solver.run(_simple_well(), fluids, schedule)
        self.assertAlmostEqual(result.fronts[0].time_s, 54.0, places=6)

    def test_default_parameters(self) -> None:
        """新增参数默认值核对（spec §6）。"""

        solver = CasingFlowSolver()
        self.assertTrue(solver.enable_buoyancy_physics)
        self.assertEqual(solver.buoyancy_correction_factor, 1.0)
        self.assertTrue(solver.enable_mixing_enhancement)
        self.assertEqual(solver.mixing_enhancement_factor, 5.0)
        self.assertEqual(solver.max_mixing_enhancement, 10.0)
        self.assertFalse(solver.has_plug)

    def test_parameter_validation(self) -> None:
        """非法参数应拒绝。"""

        with self.assertRaises(ValueError):
            CasingFlowSolver(buoyancy_correction_factor=-0.1)
        with self.assertRaises(ValueError):
            CasingFlowSolver(mixing_enhancement_factor=-1.0)
        with self.assertRaises(ValueError):
            CasingFlowSolver(max_mixing_enhancement=0.5)


if __name__ == "__main__":
    unittest.main()