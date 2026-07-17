"""
regime_classifiers 单元测试（T0-3 浮力数阈值分类 + T0-6 停泵有限时间衰减）。

使用 SimpleNamespace mock result/fluids/schedule，不跑真实求解器。

公式来源：
- T0-3：Zhang & Frigaard 2023, §3.1.3（b 阈值 0/20/80）
- T0-6：Moyers-González et al. 2007, 式 3.35（冻结判据）/ 式 3.40（冻结时间）
"""

import math
import unittest
from types import SimpleNamespace

import numpy as np

from cemdisp.data.pumping_schedule import PumpingStageEvent
from cemdisp.diagnostics.regime_classifiers import (
    BuoyancyRegimeResult,
    ShutdownDecayResult,
    classify_buoyancy_regime,
    classify_buoyancy_regime_from_result,
    compute_shutdown_decay,
)


def _make_geom(e=0.3, gap_m=0.04, inc_deg=0.0, hole_mm=220.0, od_mm=139.7):
    """构造 mock geom（键与 annulus_d2dga._build_geom 一致）。"""
    return {
        "e": np.array([e, e]),
        "b": np.full((2, 2), gap_m),
        "inc_deg": np.array([inc_deg, inc_deg]),
        "hole_mm": np.array([hole_mm, hole_mm]),
        "od_mm": np.array([od_mm, od_mm]),
    }


def _make_result(geom=None, b_number=50.0):
    return SimpleNamespace(
        geom=geom if geom is not None else _make_geom(),
        summary={"buoyancy_number": b_number},
    )


def _make_fluid(name, density, tau_y=None, mu_p=0.05):
    return SimpleNamespace(
        name=name,
        density_kg_m3=density,
        yield_stress_pa=tau_y,
        plastic_viscosity_pa_s=mu_p,
    )


def _make_schedule(with_shutdown=True, rate_m3_min=1.2, string_tag=False):
    """构造 mock schedule；string_tag=True 时用裸字符串 event_tag 测试兼容性。"""
    cement_tag = "INJECT_CEMENT" if string_tag else PumpingStageEvent.INJECT_CEMENT
    shutdown_tag = "SHUTDOWN" if string_tag else PumpingStageEvent.SHUTDOWN
    steps = [
        SimpleNamespace(
            step_name="注水泥",
            fluid_name="cement",
            volume_m3=10.0,
            rate_m3_min=rate_m3_min,
            event_tag=cement_tag,
        )
    ]
    if with_shutdown:
        steps.append(
            SimpleNamespace(
                step_name="停泵",
                fluid_name="none",
                volume_m3=0.0,
                rate_m3_min=0.0,
                event_tag=shutdown_tag,
            )
        )
    return SimpleNamespace(steps=tuple(steps))


class TestBuoyancyRegime(unittest.TestCase):
    """T0-3 浮力数阈值分类（Zhang 2023 §3.1.3）。"""

    def test_forbidden(self) -> None:
        """b=-5 → forbidden（密度不稳定，禁止）。"""
        res = classify_buoyancy_regime(-5.0)
        self.assertEqual(res.regime, "forbidden")
        self.assertIn("密度", res.design_advice)
        self.assertIn("重浆", res.design_advice)

    def test_highly_dispersive(self) -> None:
        """b=10 → highly_dispersive。"""
        res = classify_buoyancy_regime(10.0)
        self.assertEqual(res.regime, "highly_dispersive")

    def test_steady_capable(self) -> None:
        """b=50 → steady_capable。"""
        res = classify_buoyancy_regime(50.0)
        self.assertEqual(res.regime, "steady_capable")

    def test_non_dispersive(self) -> None:
        """b=100 → non_dispersive。"""
        res = classify_buoyancy_regime(100.0)
        self.assertEqual(res.regime, "non_dispersive")

    def test_boundary_values(self) -> None:
        """边界归属：b=0→highly_dispersive，b=20→steady_capable，b=80→non_dispersive。"""
        self.assertEqual(classify_buoyancy_regime(0.0).regime, "highly_dispersive")
        self.assertEqual(classify_buoyancy_regime(19.999).regime, "highly_dispersive")
        self.assertEqual(classify_buoyancy_regime(20.0).regime, "steady_capable")
        self.assertEqual(classify_buoyancy_regime(79.999).regime, "steady_capable")
        self.assertEqual(classify_buoyancy_regime(80.0).regime, "non_dispersive")
        # 负值即使接近 0 也 forbidden
        self.assertEqual(classify_buoyancy_regime(-1.0e-12).regime, "forbidden")

    def test_nan_raises(self) -> None:
        """NaN 输入抛 ValueError。"""
        with self.assertRaises(ValueError):
            classify_buoyancy_regime(float("nan"))

    def test_from_result(self) -> None:
        """便捷函数从 result.summary['buoyancy_number'] 读取。"""
        res = classify_buoyancy_regime_from_result(_make_result(b_number=50.0))
        self.assertEqual(res.regime, "steady_capable")
        self.assertAlmostEqual(res.b_number, 50.0)

    def test_to_dict(self) -> None:
        res = classify_buoyancy_regime(50.0)
        d = res.to_dict()
        self.assertEqual(
            set(d.keys()), {"b_number", "regime", "design_advice"}
        )
        self.assertEqual(d["regime"], "steady_capable")

    def test_result_frozen(self) -> None:
        """结果 dataclass 不可变。"""
        res = classify_buoyancy_regime(50.0)
        with self.assertRaises(Exception):
            res.regime = "forbidden"  # type: ignore[misc]


class TestShutdownDecay(unittest.TestCase):
    """T0-6 停泵有限时间衰减（Moyers-González 2007 式 3.35/3.40）。"""

    def test_satisfied_vertical_well(self) -> None:
        """垂直井（β=0 → S_buoy=0）+ 高屈服应力 → 判据满足，t_s 有限且量级合理。"""
        fluids = [
            _make_fluid("mud", density=1400.0, tau_y=10.0, mu_p=0.04),
            _make_fluid("cement", density=1900.0, tau_y=15.0, mu_p=0.06),
        ]
        res = compute_shutdown_decay(
            _make_result(geom=_make_geom(e=0.3, inc_deg=0.0)),
            fluids,
            _make_schedule(rate_m3_min=1.2),
        )
        self.assertIsInstance(res, ShutdownDecayResult)
        self.assertTrue(res.condition_satisfied)
        self.assertTrue(math.isfinite(res.freeze_time_s))
        self.assertGreater(res.freeze_time_s, 0.0)
        # 量级检查：典型参数下 t_s 为秒级到分钟级，不应爆炸
        self.assertLess(res.freeze_time_s, 1.0e4)
        self.assertAlmostEqual(res.tau_y_min, 10.0)
        self.assertIn("式 3.35 满足", res.physical_interpretation)

    def test_satisfied_inclined_small_density_contrast(self) -> None:
        """斜井小密度差：S_buoy≈13.9 Pa < τ_Y,min/(1+e)≈15.4 Pa → 满足。"""
        # S_buoy = π·0.0899·20·9.81·sin30°/2 ≈ 13.86 Pa
        # τ_support = 20/1.3 ≈ 15.38 Pa
        fluids = [
            _make_fluid("mud", density=1580.0, tau_y=20.0, mu_p=0.05),
            _make_fluid("cement", density=1600.0, tau_y=25.0, mu_p=0.06),
        ]
        res = compute_shutdown_decay(
            _make_result(geom=_make_geom(e=0.3, inc_deg=30.0)),
            fluids,
            _make_schedule(rate_m3_min=1.2),
        )
        self.assertTrue(res.condition_satisfied)
        self.assertTrue(math.isfinite(res.freeze_time_s))
        self.assertGreater(res.freeze_time_s, 0.0)

    def test_not_satisfied_newtonian(self) -> None:
        """τ_y=0（牛顿流体）→ 不满足 + inf（无屈服则不停滞）。"""
        fluids = [
            _make_fluid("mud", density=1400.0, tau_y=None, mu_p=0.04),
            _make_fluid("cement", density=1900.0, tau_y=0.0, mu_p=0.06),
        ]
        res = compute_shutdown_decay(
            _make_result(geom=_make_geom(e=0.3, inc_deg=0.0)),
            fluids,
            _make_schedule(),
        )
        self.assertFalse(res.condition_satisfied)
        self.assertEqual(res.freeze_time_s, math.inf)
        self.assertAlmostEqual(res.tau_y_min, 0.0)
        self.assertIn("牛顿", res.physical_interpretation)

    def test_not_satisfied_large_density_contrast(self) -> None:
        """斜井大密度差：S_buoy≈346 Pa >> τ_Y,min/(1+e)≈3.8 Pa → 不满足 + inf。"""
        # S_buoy = π·0.0899·500·9.81·sin30°/2 ≈ 346 Pa
        fluids = [
            _make_fluid("mud", density=1100.0, tau_y=5.0, mu_p=0.04),
            _make_fluid("cement", density=1600.0, tau_y=8.0, mu_p=0.06),
        ]
        res = compute_shutdown_decay(
            _make_result(geom=_make_geom(e=0.3, inc_deg=30.0)),
            fluids,
            _make_schedule(),
        )
        self.assertFalse(res.condition_satisfied)
        self.assertEqual(res.freeze_time_s, math.inf)
        self.assertAlmostEqual(res.tau_y_min, 5.0)
        self.assertIn("式 3.35 不满足", res.physical_interpretation)

    def test_no_shutdown_event(self) -> None:
        """无 SHUTDOWN 事件 → 诊断不适用（False + inf + 说明）。"""
        fluids = [_make_fluid("mud", density=1400.0, tau_y=10.0)]
        res = compute_shutdown_decay(
            _make_result(),
            fluids,
            _make_schedule(with_shutdown=False),
        )
        self.assertFalse(res.condition_satisfied)
        self.assertEqual(res.freeze_time_s, math.inf)
        # tau_y_min 字段仍应正确返回（信息有用）
        self.assertAlmostEqual(res.tau_y_min, 10.0)
        self.assertIn("不适用", res.physical_interpretation)

    def test_string_event_tag_compatible(self) -> None:
        """裸字符串 event_tag='SHUTDOWN' 也能触发诊断（mock 兼容性）。"""
        fluids = [_make_fluid("mud", density=1400.0, tau_y=10.0)]
        res = compute_shutdown_decay(
            _make_result(geom=_make_geom(inc_deg=0.0)),
            fluids,
            _make_schedule(string_tag=True),
        )
        self.assertTrue(res.condition_satisfied)

    def test_zero_rate_before_shutdown(self) -> None:
        """停泵前排量为零 → w₀=0 → t_s=0（已静止，立即冻结）。"""
        fluids = [_make_fluid("mud", density=1400.0, tau_y=10.0)]
        res = compute_shutdown_decay(
            _make_result(geom=_make_geom(inc_deg=0.0)),
            fluids,
            _make_schedule(rate_m3_min=0.0),
        )
        self.assertTrue(res.condition_satisfied)
        self.assertEqual(res.freeze_time_s, 0.0)

    def test_tau_y_min_takes_minimum(self) -> None:
        """τ_Y,min 取各流体最小值；None（牛顿）按 0 计。"""
        fluids = [
            _make_fluid("a", density=1400.0, tau_y=8.0),
            _make_fluid("b", density=1500.0, tau_y=None),
            _make_fluid("c", density=1600.0, tau_y=15.0),
        ]
        res = compute_shutdown_decay(
            _make_result(), fluids, _make_schedule()
        )
        self.assertAlmostEqual(res.tau_y_min, 0.0)
        self.assertFalse(res.condition_satisfied)

    def test_freeze_time_formula_order_of_magnitude(self) -> None:
        """核对式 3.40 重建的数值：手算 t_s≈4.1 s（垂直井算例）。"""
        # e=0.3, d=0.02, μ=0.04, ρ_min=1400, τ_Y,min=10, β=0
        # w₀ = (1.2/60) / (π/4·(0.22²−0.1397²)) ≈ 0.882 m/s
        # α₁ = 0.04·1.3/(1400·0.02²) ≈ 0.0929 1/s
        # α₂ = 2·1.3·(10/1.3)/(1400·π²·0.02) ≈ 0.0724 m/s²
        # z₀ = 0.882·1.3/π ≈ 0.365 m/s
        # t_s = ln(1 + z₀·α₁/α₂)/α₁ ≈ 4.13 s
        fluids = [
            _make_fluid("mud", density=1400.0, tau_y=10.0, mu_p=0.04),
            _make_fluid("cement", density=1900.0, tau_y=15.0, mu_p=0.06),
        ]
        res = compute_shutdown_decay(
            _make_result(geom=_make_geom(e=0.3, gap_m=0.04, inc_deg=0.0)),
            fluids,
            _make_schedule(rate_m3_min=1.2),
        )
        self.assertTrue(res.condition_satisfied)
        self.assertAlmostEqual(res.freeze_time_s, 4.13, delta=0.5)

    def test_to_dict(self) -> None:
        fluids = [_make_fluid("mud", density=1400.0, tau_y=10.0)]
        res = compute_shutdown_decay(
            _make_result(geom=_make_geom(inc_deg=0.0)), fluids, _make_schedule()
        )
        d = res.to_dict()
        self.assertEqual(
            set(d.keys()),
            {"condition_satisfied", "freeze_time_s", "tau_y_min", "physical_interpretation"},
        )
        self.assertTrue(d["condition_satisfied"])

    def test_result_frozen(self) -> None:
        """结果 dataclass 不可变。"""
        fluids = [_make_fluid("mud", density=1400.0, tau_y=10.0)]
        res = compute_shutdown_decay(
            _make_result(geom=_make_geom(inc_deg=0.0)), fluids, _make_schedule()
        )
        with self.assertRaises(Exception):
            res.condition_satisfied = False  # type: ignore[misc]


if __name__ == "__main__":
    unittest.main()
