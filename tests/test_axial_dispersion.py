"""
管内轴向弥散模型测试。

验证 Taylor-Aris 型轴向弥散在套管内 1D 求解器中的行为：
- 弥散系数公式对不同流变模型正确
- 弥散后时间线包含平滑过渡事件
- 边界桥接正确处理多相共存
"""

import unittest
import math

from cemdisp.data.fluid_spec import FluidSpec, FluidRole, RheologyModel
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep, PumpingStageEvent
from cemdisp.data.well_spec import WellSpec, DepthValuePoint, EvaluationWindow
from cemdisp.transport1d.casing_flow import CasingFlowSolver
from cemdisp.transport1d.shoe_timeline import ShoeEventKind
from cemdisp.transport1d.pipe_exit_state import PipeExitState
from cemdisp.models2d.boundary_bridge import _phase_fractions_from_state, _phase_fractions_for_fluid


class TestDispersionCoefficient(unittest.TestCase):
    """测试弥散系数计算。"""

    def _make_solver(self, alpha: float = 0.2) -> CasingFlowSolver:
        return CasingFlowSolver(
            enable_axial_dispersion=True,
            dispersion_alpha=alpha,
        )

    def test_newtonian_dispersion_formula(self) -> None:
        """牛顿流体 Taylor-Aris: D_eff = D_mol·Pe²/192 = U²R²/(48·D_mol)。"""
        solver = self._make_solver()
        fluid = FluidSpec(
            name="water",
            role=FluidRole.MUD,
            density_kg_m3=1000.0,
            rheology_model=RheologyModel.NEWTONIAN,
            plastic_viscosity_pa_s=0.001,
        )
        R = 0.05

        # 低速（Pe 极小）处于 Taylor-Aris 公式自洽区：U²R²/(48·D_mol)
        U_low = 1.0e-8
        D_eff = solver._compute_dispersion_coefficient(R, fluid, U_low)
        expected = (U_low ** 2) * (R ** 2) / (48.0 * 1.0e-9)
        self.assertAlmostEqual(D_eff, expected, places=15)

        # 固井流速（Pe~10^8）：Taylor 渐近值不可用，截断到对流尺度上限 α·U·R
        U_high = 1.0
        D_eff = solver._compute_dispersion_coefficient(R, fluid, U_high)
        self.assertAlmostEqual(D_eff, 0.2 * U_high * R, places=12)

    def test_power_law_shear_thinning_reduces_dispersion(self) -> None:
        """幂律 n<1（剪切稀化）→ 速度剖面变平 → 弥散减小（Batot 2016）。"""
        solver = self._make_solver()
        fluid_n05 = FluidSpec(
            name="pl05", role=FluidRole.MUD, density_kg_m3=1200.0,
            rheology_model=RheologyModel.POWER_LAW, power_law_n=0.5, consistency_k=0.5,
        )
        fluid_n10 = FluidSpec(
            name="pl10", role=FluidRole.MUD, density_kg_m3=1200.0,
            rheology_model=RheologyModel.POWER_LAW, power_law_n=1.0, consistency_k=0.5,
        )
        R = 0.05
        # 低速构造 Taylor-Aris 公式自洽区（固井流速下两分支同被对流上限截断，
        # 无法区分 Batot 流变因子）
        U = 1.0e-8
        D_eff_05 = solver._compute_dispersion_coefficient(R, fluid_n05, U)
        D_eff_10 = solver._compute_dispersion_coefficient(R, fluid_n10, U)
        # n=0.5 速度剖面更扁平 → 弥散更小
        self.assertLess(D_eff_05, D_eff_10)

    def test_zero_velocity_returns_zero(self) -> None:
        """停泵时 (U=0) 弥散系数为 0。"""
        solver = self._make_solver()
        fluid = FluidSpec(
            name="mud",
            role=FluidRole.MUD,
            density_kg_m3=1200.0,
            rheology_model=RheologyModel.NEWTONIAN,
            plastic_viscosity_pa_s=0.01,
        )
        D_eff = solver._compute_dispersion_coefficient(0.05, fluid, 0.0)
        self.assertEqual(D_eff, 0.0)

    def test_bingham_dispersion_fan_wang(self) -> None:
        """Bingham 流体 → Fan & Wang (1966) 公式，τ_y=0 时退化为牛顿。"""
        solver = self._make_solver()
        fluid = FluidSpec(
            name="mud",
            role=FluidRole.MUD,
            density_kg_m3=1200.0,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=0.05,
            yield_stress_pa=0.0,  # 零屈服 → 牛顿极限
        )
        R = 0.05
        # 低速区验证 τ_y=0 的牛顿退化（未被对流上限截断）
        U = 1.0e-8

        D_eff = solver._compute_dispersion_coefficient(R, fluid, U)
        expected = (U ** 2) * (R ** 2) / (48.0 * 1.0e-9)
        self.assertAlmostEqual(D_eff, expected, places=15)

    def test_bingham_yield_stress_reduces_dispersion(self) -> None:
        """Bingham 屈服抑制弥散：低速区 Fan & Wang 可测，固井流速被上限截断。"""
        solver = self._make_solver()
        def _make_bingham(tau_y: float) -> FluidSpec:
            return FluidSpec(
                name="bh", role=FluidRole.MUD, density_kg_m3=1200.0,
                rheology_model=RheologyModel.BINGHAM,
                plastic_viscosity_pa_s=0.05, yield_stress_pa=tau_y,
            )
        R = 0.05
        # 低速构造 Fan & Wang 公式自洽区：屈服应力存在时塞流增大、弥散减小
        # （τ_y=1 vs 0 相差约 4 个数量级；τ_y=1 vs 20 因 ξ₀ 饱和 0.999 不可区分）
        U = 1.0e-8
        D_newtonian = solver._compute_dispersion_coefficient(R, _make_bingham(0.0), U)
        D_yield = solver._compute_dispersion_coefficient(R, _make_bingham(1.0), U)
        self.assertLess(D_yield, D_newtonian)  # 屈服应力存在 -> 塞流增大 -> 弥散减小

        # 固井流速下屈服梯度不可测：均被对流尺度上限截断到 α·U·R
        D_cap_newtonian = solver._compute_dispersion_coefficient(R, _make_bingham(0.0), 1.0)
        D_cap_yield = solver._compute_dispersion_coefficient(R, _make_bingham(20.0), 1.0)
        self.assertAlmostEqual(D_cap_newtonian, D_cap_yield, places=12)

class TestDispersionTimeline(unittest.TestCase):
    """测试弥散后时间线生成。"""

    def _make_well_spec(self) -> WellSpec:
        return WellSpec(
            well_name="test",
            top_md_m=1.0,
            bottom_md_m=1000.0,
            shoe_md_m=1000.0,
            hole_diameter_profile=(DepthValuePoint(1.0, 216.0), DepthValuePoint(1000.0, 216.0)),
            inclination_profile=(DepthValuePoint(1.0, 0.0), DepthValuePoint(1000.0, 0.0)),
            standoff_profile=(DepthValuePoint(1.0, 0.8), DepthValuePoint(1000.0, 0.8)),
            liner_od_mm=139.7,
            liner_id_mm=124.3,
            evaluation_windows=(EvaluationWindow("target", 1.0, 1000.0),),
        )

    def _make_schedule(self) -> PumpingSchedule:
        # 注入 20 m³，速率 2 m³/min，确保水泥前缘能到达鞋口
        return PumpingSchedule(steps=(
            PumpingScheduleStep(
                step_name="注水泥",
                fluid_name="cement",
                volume_m3=20.0,
                rate_m3_min=2.0,
                event_tag=PumpingStageEvent.INJECT_CEMENT,
            ),
        ))

    def test_dispersion_disabled_returns_original(self) -> None:
        """enable_axial_dispersion=False 时时间线不变。"""
        solver = CasingFlowSolver(enable_axial_dispersion=False)
        well = self._make_well_spec()
        fluid = FluidSpec(
            name="cement",
            role=FluidRole.LEAD,
            density_kg_m3=1900.0,
            rheology_model=RheologyModel.NEWTONIAN,
            plastic_viscosity_pa_s=0.3,
        )
        schedule = self._make_schedule()

        result = solver.run(well, (fluid,), schedule)
        events = result.shoe_timeline.events

        # 原始时间线应只有 FRONT_ARRIVAL 等离散事件
        arrival_events = [e for e in events if e.kind == ShoeEventKind.FRONT_ARRIVAL]
        self.assertEqual(len(arrival_events), 1)
        self.assertEqual(arrival_events[0].phase_fractions, (("cement", 1.0),))

    def test_dispersion_enabled_creates_transition_events(self) -> None:
        """enable_axial_dispersion=True 时产生多相过渡事件。"""
        solver = CasingFlowSolver(enable_axial_dispersion=True)
        well = self._make_well_spec()
        cement = FluidSpec(
            name="cement",
            role=FluidRole.LEAD,
            density_kg_m3=1900.0,
            rheology_model=RheologyModel.NEWTONIAN,
            plastic_viscosity_pa_s=0.3,
        )
        mud = FluidSpec(
            name="mud",
            role=FluidRole.MUD,
            density_kg_m3=1200.0,
            rheology_model=RheologyModel.NEWTONIAN,
            plastic_viscosity_pa_s=0.01,
        )
        schedule = self._make_schedule()

        result = solver.run(well, (mud, cement), schedule)
        events = result.shoe_timeline.events

        # 弥散后应有多相共存事件（分数不全是 0 或 1）
        transition_events = [
            e for e in events
            if e.kind == ShoeEventKind.FRONT_ARRIVAL
            and len(e.phase_fractions) > 1
        ]
        self.assertGreater(len(transition_events), 0)

        # 验证存在中间浓度事件
        has_intermediate = any(
            any(0.0 < frac < 1.0 for _, frac in e.phase_fractions)
            for e in transition_events
        )
        self.assertTrue(has_intermediate, "应存在中间浓度的过渡事件")

    def test_dispersed_events_sorted_by_time(self) -> None:
        """弥散事件按时间升序排列。"""
        solver = CasingFlowSolver(enable_axial_dispersion=True)
        well = self._make_well_spec()
        cement = FluidSpec(
            name="cement",
            role=FluidRole.LEAD,
            density_kg_m3=1900.0,
            rheology_model=RheologyModel.NEWTONIAN,
            plastic_viscosity_pa_s=0.3,
        )
        mud = FluidSpec(
            name="mud",
            role=FluidRole.MUD,
            density_kg_m3=1200.0,
            rheology_model=RheologyModel.NEWTONIAN,
            plastic_viscosity_pa_s=0.01,
        )
        schedule = self._make_schedule()

        result = solver.run(well, (mud, cement), schedule)
        events = result.shoe_timeline.events
        times = [e.time_s for e in events]
        self.assertEqual(times, sorted(times))


class TestBoundaryBridgeMultiPhase(unittest.TestCase):
    """测试边界桥接多相映射。"""

    def test_phase_fractions_from_state_single_phase(self) -> None:
        """单相状态与原有 _phase_fractions_for_fluid 行为一致。"""
        fluids = (
            FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1200.0, rheology_model=RheologyModel.NEWTONIAN, plastic_viscosity_pa_s=0.01),
            FluidSpec(name="cement", role=FluidRole.LEAD, density_kg_m3=1900.0, rheology_model=RheologyModel.NEWTONIAN, plastic_viscosity_pa_s=0.3),
        )
        state = PipeExitState(
            time_s=0.0,
            flow_rate_m3_s=0.01,
            stage_name="test",
            phase_fractions=(("cement", 1.0),),
        )

        result_from_state = _phase_fractions_from_state(state, fluids)
        result_from_fluid = _phase_fractions_for_fluid("cement", fluids)
        self.assertEqual(result_from_state, result_from_fluid)

    def test_phase_fractions_from_state_multi_phase(self) -> None:
        """多相状态按体积分数加权合并。"""
        fluids = (
            FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1200.0, rheology_model=RheologyModel.NEWTONIAN, plastic_viscosity_pa_s=0.01),
            FluidSpec(name="cement", role=FluidRole.LEAD, density_kg_m3=1900.0, rheology_model=RheologyModel.NEWTONIAN, plastic_viscosity_pa_s=0.3),
        )
        # 50% cement + 50% mud
        state = PipeExitState(
            time_s=0.0,
            flow_rate_m3_s=0.01,
            stage_name="test",
            phase_fractions=(("cement", 0.5), ("mud", 0.5)),
        )

        result = _phase_fractions_from_state(state, fluids)
        result_dict = dict(result)
        self.assertAlmostEqual(result_dict.get("cement", 0.0), 0.5, places=6)
        self.assertAlmostEqual(result_dict.get("mud", 0.0), 0.5, places=6)

    def test_phase_fractions_from_state_split_cement(self) -> None:
        """split_cement_phases=True 时多相映射仍正确。"""
        fluids = (
            FluidSpec(name="lead", role=FluidRole.LEAD, density_kg_m3=1800.0, rheology_model=RheologyModel.NEWTONIAN, plastic_viscosity_pa_s=0.3),
            FluidSpec(name="tail", role=FluidRole.TAIL, density_kg_m3=1900.0, rheology_model=RheologyModel.NEWTONIAN, plastic_viscosity_pa_s=0.4),
        )
        state = PipeExitState(
            time_s=0.0,
            flow_rate_m3_s=0.01,
            stage_name="test",
            phase_fractions=(("lead", 0.7), ("tail", 0.3)),
        )

        result = _phase_fractions_from_state(state, fluids, split_cement_phases=True)
        result_dict = dict(result)
        self.assertAlmostEqual(result_dict.get("lead", 0.0), 0.7, places=6)
        self.assertAlmostEqual(result_dict.get("tail", 0.0), 0.3, places=6)


class TestDispersionAlphaParameter(unittest.TestCase):
    """测试弥散系数参数。"""

    def test_default_alpha_is_025(self) -> None:
        """默认 dispersion_alpha = 0.25。"""
        solver = CasingFlowSolver()
        self.assertEqual(solver.dispersion_alpha, 0.25)

    def test_custom_alpha(self) -> None:
        """自定义 dispersion_alpha 正确保存。"""
        solver = CasingFlowSolver(dispersion_alpha=0.5)
        self.assertEqual(solver.dispersion_alpha, 0.5)

    def test_negative_alpha_raises(self) -> None:
        """负的 dispersion_alpha 应抛出 ValueError。"""
        with self.assertRaises(ValueError):
            CasingFlowSolver(dispersion_alpha=-0.1)


if __name__ == "__main__":
    unittest.main()
