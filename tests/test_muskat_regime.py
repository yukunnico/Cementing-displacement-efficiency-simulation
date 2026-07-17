"""
muskat_regime 单元测试（T0-2 Muskat 三 regime 稳定性诊断）。

使用 SimpleNamespace mock result/fluids/well_spec，不跑真实求解器；
判据测试优先使用纯函数入口 classify_muskat_regime(m, b, e)。

公式来源：
- Bararpour & Frigaard 2025, JFM 1009, A15, 式 3.8-3.13（w_f/w_finger/Δw 与三 regime 判据）
- Pelipenko & Frigaard 2004, 式 3.33/3.34（宽/窄边补充判据）
- 牛顿等密度极限 m_c≈1.5（Bararpour 2025 图 15a）
"""

import json
import math
import unittest
from types import SimpleNamespace

import numpy as np

from cemdisp.diagnostics.muskat_regime import (
    MuskatRegimeResult,
    _buoyant_mobility_I2,
    _finger_velocity,
    _front_flux,
    _isotropic_flux_q0,
    _mean_mobility_I1,
    classify_muskat_regime,
    compute_muskat_regime,
)


def _make_geom(e=0.3, gap_m=0.04):
    """构造 mock geom（键与 annulus_d2dga._build_geom 一致）。"""
    return {
        "e": np.array([e, e]),
        "b": np.full((2, 2), gap_m),
    }


def _make_result(geom=None, b_number=100.0):
    return SimpleNamespace(
        geom=geom if geom is not None else _make_geom(),
        summary={"buoyancy_number": b_number},
    )


def _make_fluid(role, density, mu_p=0.05, consistency_k=None):
    return SimpleNamespace(
        role=role,
        density_kg_m3=density,
        plastic_viscosity_pa_s=mu_p,
        consistency_k=consistency_k,
    )


def _make_fluids(mu_mud=0.1, mu_cement=0.05, rho_mud=1400.0, rho_cement=1900.0):
    return [
        _make_fluid("mud", rho_mud, mu_p=mu_mud),
        _make_fluid("lead", rho_cement, mu_p=mu_cement),
    ]


class TestNewtonianIsoDensityLimit(unittest.TestCase):
    """牛顿等密度极限：m > m_c≈1.5 → unstable（Bararpour 2025 图 15a）。"""

    def test_unstable_above_critical_viscosity_ratio(self) -> None:
        """m=2.0 > m_c，b=0 → unstable，且 Δw 全区间为正（无 c_critical）。"""
        res = classify_muskat_regime(2.0, 0.0, 0.0, n_grid=501)
        self.assertEqual(res.regime, "unstable")
        self.assertGreater(res.delta_w_min, 0.0)
        self.assertTrue(math.isnan(res.c_critical))
        profile = np.asarray(res.delta_w_profile)
        self.assertEqual(profile.size, 501)
        self.assertTrue(np.all(profile > 0.0))

    def test_critical_ratio_bracketed(self) -> None:
        """临界黏度比 m_c≈1.5：m=1.4 时 Δw_min<0，m=1.6 时 Δw_min>0。"""
        below = classify_muskat_regime(1.4, 0.0, 0.0, n_grid=501, include_profiles=False)
        above = classify_muskat_regime(1.6, 0.0, 0.0, n_grid=501, include_profiles=False)
        self.assertLess(below.delta_w_min, 0.0)
        self.assertNotEqual(below.regime, "unstable")
        self.assertGreater(above.delta_w_min, 0.0)
        self.assertEqual(above.regime, "unstable")

    def test_equal_viscosity_not_unstable(self) -> None:
        """m=1.0（等黏等密度）→ partial_penetration，Δw_min<0，无激波。"""
        res = classify_muskat_regime(1.0, 0.0, 0.0, n_grid=501, include_profiles=False)
        self.assertEqual(res.regime, "partial_penetration")
        self.assertLess(res.delta_w_min, 0.0)
        self.assertFalse(res.shock_detected)


class TestBuoyancyStabilization(unittest.TestCase):
    """浮力数 b 的稳定化作用（b>0 密度稳定，重顶替轻）。"""

    def test_large_buoyancy_stable(self) -> None:
        """m=2.0（失稳黏度比）+ b=100 → stable，且 c* < c_critical。"""
        res = classify_muskat_regime(2.0, 100.0, 0.0, n_grid=501, include_profiles=False)
        self.assertEqual(res.regime, "stable")
        self.assertLess(res.delta_w_min, 0.0)
        self.assertTrue(math.isfinite(res.c_critical))
        self.assertLess(res.c_star, res.c_critical)

    def test_buoyancy_monotonicity(self) -> None:
        """浮力单调性：b=0/10/100 的 Δw_min 单调递减（越大越稳定）。"""
        dw0 = classify_muskat_regime(2.0, 0.0, 0.0, n_grid=501, include_profiles=False).delta_w_min
        dw10 = classify_muskat_regime(2.0, 10.0, 0.0, n_grid=501, include_profiles=False).delta_w_min
        dw100 = classify_muskat_regime(2.0, 100.0, 0.0, n_grid=501, include_profiles=False).delta_w_min
        self.assertGreater(dw0, dw10)
        self.assertGreater(dw10, dw100)

    def test_b10_less_stable_than_b100(self) -> None:
        """b=10 比 b=100 更不稳定（Δw_min 更大）。"""
        r10 = classify_muskat_regime(2.0, 10.0, 0.0, n_grid=501, include_profiles=False)
        r100 = classify_muskat_regime(2.0, 100.0, 0.0, n_grid=501, include_profiles=False)
        self.assertGreater(r10.delta_w_min, r100.delta_w_min)
        self.assertEqual(r100.regime, "stable")

    def test_moderate_buoyancy_partial(self) -> None:
        """m=2.0 + b=5 → partial_penetration（介于 unstable 与 stable 之间）。"""
        res = classify_muskat_regime(2.0, 5.0, 0.0, n_grid=501, include_profiles=False)
        self.assertEqual(res.regime, "partial_penetration")
        self.assertLess(res.delta_w_min, 0.0)
        self.assertTrue(math.isfinite(res.c_critical))

    def test_shock_detected_when_unstable(self) -> None:
        """失稳算例（m=2.0, b=0）通量存在凸区间 → shock_detected=True。"""
        res = classify_muskat_regime(2.0, 0.0, 0.0, n_grid=501, include_profiles=False)
        self.assertTrue(res.shock_detected)


class TestWideNarrowSideCriteria(unittest.TestCase):
    """宽/窄边补充判据（Pelipenko 2004 式 3.33/3.34 牛顿形式）。"""

    def test_high_m_both_sides_unstable(self) -> None:
        """m=2.0, b=0：宽边 m>1 失稳、窄边 1/m<1 失稳 → 均 True。"""
        res = classify_muskat_regime(2.0, 0.0, 0.0, n_grid=101, include_profiles=False)
        self.assertTrue(res.wide_side_unstable)
        self.assertTrue(res.narrow_side_unstable)

    def test_low_m_both_sides_stable(self) -> None:
        """m=0.5, b=0：宽边 m<1、窄边 1/m>1 → 均 False。"""
        res = classify_muskat_regime(0.5, 0.0, 0.0, n_grid=101, include_profiles=False)
        self.assertFalse(res.wide_side_unstable)
        self.assertFalse(res.narrow_side_unstable)

    def test_buoyancy_stabilizes_both_sides(self) -> None:
        """m=2.0 + b=100：浮力使宽/窄边判据均转为稳定。"""
        res = classify_muskat_regime(2.0, 100.0, 0.3, n_grid=101, include_profiles=False)
        self.assertFalse(res.wide_side_unstable)
        self.assertFalse(res.narrow_side_unstable)


class TestNumericalSafety(unittest.TestCase):
    """端点奇异保护与数值安全。"""

    def test_endpoints_finite_no_nan(self) -> None:
        """c̄→0/1 端点不崩溃、无 NaN；网格端点恰为 0 与 1。"""
        for m, b in [(0.5, 0.0), (1.0, 0.0), (2.0, 0.0), (2.0, 100.0), (3.0, 10.0)]:
            res = classify_muskat_regime(m, b, 0.0, n_grid=201)
            grid = np.asarray(res.c_bar_grid)
            profile = np.asarray(res.delta_w_profile)
            self.assertAlmostEqual(grid[0], 0.0)
            self.assertAlmostEqual(grid[-1], 1.0)
            self.assertTrue(np.all(np.isfinite(profile)), msg=f"m={m} b={b} 剖面含非有限值")
            self.assertTrue(math.isfinite(res.delta_w_min))
            self.assertTrue(math.isfinite(res.delta_w_max))
            self.assertTrue(math.isfinite(res.c_star))

    def test_closure_endpoint_analytic_values(self) -> None:
        """牛顿闭包端点解析值：q0(0)=0、q0(1)=1、I2(0)=I2(1)=0、I1 恒正、F 端点 0/1。"""
        c = np.array([0.0, 1.0])
        for m in [0.5, 1.0, 2.0]:
            i1 = _mean_mobility_I1(c, m)
            i2 = _buoyant_mobility_I2(c, m)
            q0 = _isotropic_flux_q0(c, m)
            flux = _front_flux(c, m, 10.0)
            self.assertTrue(np.all(i1 > 0.0))
            self.assertAlmostEqual(float(i2[0]), 0.0, places=12)
            self.assertAlmostEqual(float(i2[1]), 0.0, places=12)
            self.assertAlmostEqual(float(q0[0]), 0.0, places=12)
            self.assertAlmostEqual(float(q0[1]), 1.0, places=12)
            self.assertAlmostEqual(float(flux[0]), 0.0, places=12)
            self.assertAlmostEqual(float(flux[1]), 1.0, places=12)

    def test_finger_velocity_endpoint_analytic_values(self) -> None:
        """指进速度端点解析值：w_finger(1)=1；w_finger(0)=m−b·√m/3（式 3.12）。"""
        c = np.array([0.0, 1.0])
        wf = _finger_velocity(c, 2.0, 0.0)
        self.assertAlmostEqual(float(wf[0]), 2.0, places=10)
        self.assertAlmostEqual(float(wf[1]), 1.0, places=10)
        wf_b = _finger_velocity(c, 2.0, 10.0)
        self.assertAlmostEqual(float(wf_b[0]), 2.0 - 10.0 * math.sqrt(2.0) / 3.0, places=10)
        self.assertAlmostEqual(float(wf_b[1]), 1.0, places=10)

    def test_invalid_inputs_raise(self) -> None:
        """非法输入抛 ValueError：m≤0、m 非有限、b 非有限、e 非有限、网格过粗。"""
        with self.assertRaises(ValueError):
            classify_muskat_regime(0.0, 0.0)
        with self.assertRaises(ValueError):
            classify_muskat_regime(-1.0, 0.0)
        with self.assertRaises(ValueError):
            classify_muskat_regime(float("nan"), 0.0)
        with self.assertRaises(ValueError):
            classify_muskat_regime(1.0, float("nan"))
        with self.assertRaises(ValueError):
            classify_muskat_regime(1.0, 0.0, float("nan"))
        with self.assertRaises(ValueError):
            classify_muskat_regime(1.0, 0.0, 0.0, n_grid=10)


class TestResultContainer(unittest.TestCase):
    """MuskatRegimeResult 容器行为：to_dict / frozen / profiles 开关。"""

    def test_to_dict_keys_and_nan_conversion(self) -> None:
        """to_dict 键完整；unstable 时 c_critical 由 NaN 转 None。"""
        res = classify_muskat_regime(2.0, 0.0, 0.0, n_grid=201)
        d = res.to_dict()
        self.assertEqual(
            set(d.keys()),
            {
                "regime", "c_critical", "c_star",
                "wide_side_unstable", "narrow_side_unstable",
                "delta_w_min", "delta_w_max",
                "viscosity_ratio", "buoyancy_number", "eccentricity",
                "shock_detected", "c_bar_grid", "delta_w_profile", "notes",
            },
        )
        self.assertEqual(d["regime"], "unstable")
        self.assertIsNone(d["c_critical"])
        self.assertIsInstance(d["c_bar_grid"], list)
        self.assertIsInstance(d["notes"], list)
        self.assertEqual(len(d["c_bar_grid"]), 201)

    def test_to_dict_json_serializable(self) -> None:
        """to_dict 结果可 JSON 序列化（含中文 notes）。"""
        res = classify_muskat_regime(2.0, 100.0, 0.3, n_grid=201)
        s = json.dumps(res.to_dict())
        self.assertIn('"stable"', s)

    def test_result_frozen(self) -> None:
        """结果 dataclass 不可变。"""
        res = classify_muskat_regime(2.0, 0.0, 0.0, n_grid=101, include_profiles=False)
        self.assertIsInstance(res, MuskatRegimeResult)
        with self.assertRaises(Exception):
            res.regime = "stable"  # type: ignore[misc]

    def test_include_profiles_false(self) -> None:
        """include_profiles=False → 网格与剖面为空 tuple，其余字段正常。"""
        res = classify_muskat_regime(2.0, 0.0, 0.0, n_grid=101, include_profiles=False)
        self.assertEqual(res.c_bar_grid, ())
        self.assertEqual(res.delta_w_profile, ())
        self.assertEqual(res.regime, "unstable")
        d = res.to_dict()
        self.assertEqual(d["c_bar_grid"], [])
        self.assertEqual(d["delta_w_profile"], [])

    def test_notes_disclose_approximations(self) -> None:
        """notes 披露牛顿近似与 I3/b 标度约定。"""
        res = classify_muskat_regime(2.0, 0.0, 0.0, n_grid=101, include_profiles=False)
        joined = "".join(res.notes)
        self.assertIn("牛顿近似", joined)
        self.assertIn("I3", joined)


class TestComputeMuskatRegime(unittest.TestCase):
    """compute_muskat_regime：从 mock result/fluids/well_spec 提取输入（Layer 5 入口）。"""

    def test_basic_extraction(self) -> None:
        """m=μ_mud/μ_cement、b 取 summary、e 取 geom['e'] 均值；结果与纯函数一致。"""
        fluids = _make_fluids(mu_mud=0.1, mu_cement=0.05)  # m=2.0
        res = compute_muskat_regime(
            _make_result(geom=_make_geom(e=0.3), b_number=0.0),
            fluids,
            n_grid=501,
            include_profiles=False,
        )
        self.assertAlmostEqual(res.viscosity_ratio, 2.0)
        self.assertAlmostEqual(res.buoyancy_number, 0.0)
        self.assertAlmostEqual(res.eccentricity, 0.3)
        self.assertEqual(res.regime, "unstable")
        ref = classify_muskat_regime(2.0, 0.0, 0.3, n_grid=501, include_profiles=False)
        self.assertEqual(res.regime, ref.regime)
        self.assertAlmostEqual(res.delta_w_min, ref.delta_w_min)
        self.assertIn("输入提取", res.notes[-1])

    def test_role_enum_compatible(self) -> None:
        """role 为带 .value 的枚举式对象也能识别（与 FluidRole 兼容）。"""
        fluids = [
            SimpleNamespace(
                role=SimpleNamespace(value="mud"),
                density_kg_m3=1400.0,
                plastic_viscosity_pa_s=0.1,
                consistency_k=None,
            ),
            SimpleNamespace(
                role=SimpleNamespace(value="lead"),
                density_kg_m3=1900.0,
                plastic_viscosity_pa_s=0.05,
                consistency_k=None,
            ),
        ]
        res = compute_muskat_regime(
            _make_result(b_number=100.0), fluids, n_grid=201, include_profiles=False
        )
        self.assertAlmostEqual(res.viscosity_ratio, 2.0)
        self.assertEqual(res.regime, "stable")

    def test_tail_used_when_no_lead(self) -> None:
        """无 lead 时退用 tail 作为水泥浆。"""
        fluids = [
            _make_fluid("mud", 1400.0, mu_p=0.1),
            _make_fluid("tail", 1900.0, mu_p=0.05),
        ]
        res = compute_muskat_regime(
            _make_result(b_number=0.0), fluids, n_grid=201, include_profiles=False
        )
        self.assertAlmostEqual(res.viscosity_ratio, 2.0)

    def test_buoyancy_fallback_recompute(self) -> None:
        """summary 缺 buoyancy_number 时按 Zhang 2022 p.8 定义重算。

        b = (ρ_c−ρ_m)·g·(gap/2)²/(μ_mud·ŵ₀) = 500·9.81·0.02²/(0.05·1.0) = 39.24
        """
        result = SimpleNamespace(summary={}, geom=_make_geom(e=0.3, gap_m=0.04))
        fluids = _make_fluids(mu_mud=0.05, mu_cement=0.05, rho_mud=1400.0, rho_cement=1900.0)
        res = compute_muskat_regime(
            result, fluids, n_grid=201, mean_velocity_m_s=1.0, include_profiles=False
        )
        self.assertAlmostEqual(res.buoyancy_number, 39.24, delta=0.01)
        self.assertIn("重算", res.notes[-1])

    def test_buoyancy_missing_no_fallback_raises(self) -> None:
        """summary 无 b 且缺回退输入（mean_velocity_m_s）→ ValueError。"""
        result = SimpleNamespace(summary={}, geom=_make_geom())
        with self.assertRaises(ValueError):
            compute_muskat_regime(result, _make_fluids(), n_grid=201)

    def test_eccentricity_standoff_fallback(self) -> None:
        """geom 缺 'e' 时由 well_spec.standoff_profile 换算 e=1−mean(standoff)。"""
        result = SimpleNamespace(
            summary={"buoyancy_number": 100.0},
            geom={"b": np.full((2, 2), 0.04)},
        )
        well_spec = SimpleNamespace(
            standoff_profile=[SimpleNamespace(value=0.7), SimpleNamespace(value=0.9)]
        )
        res = compute_muskat_regime(
            result, _make_fluids(), well_spec, n_grid=201, include_profiles=False
        )
        self.assertAlmostEqual(res.eccentricity, 0.2)
        self.assertIn("standoff", res.notes[-1])

    def test_eccentricity_default_zero(self) -> None:
        """无偏心信息（geom 无 'e'、无 well_spec）→ e=0 同心近似。"""
        result = SimpleNamespace(
            summary={"buoyancy_number": 100.0},
            geom={"b": np.full((2, 2), 0.04)},
        )
        res = compute_muskat_regime(
            result, _make_fluids(), n_grid=201, include_profiles=False
        )
        self.assertAlmostEqual(res.eccentricity, 0.0)

    def test_consistency_k_fallback(self) -> None:
        """塑性黏度缺失时用一致性系数 consistency_k 作等效牛顿黏度。"""
        fluids = [
            _make_fluid("mud", 1400.0, mu_p=None, consistency_k=0.1),
            _make_fluid("lead", 1900.0, mu_p=None, consistency_k=0.05),
        ]
        res = compute_muskat_regime(
            _make_result(b_number=0.0), fluids, n_grid=201, include_profiles=False
        )
        self.assertAlmostEqual(res.viscosity_ratio, 2.0)

    def test_missing_mud_raises(self) -> None:
        """fluids 中无 role='mud' → ValueError。"""
        fluids = [_make_fluid("lead", 1900.0)]
        with self.assertRaises(ValueError):
            compute_muskat_regime(_make_result(), fluids, n_grid=201)

    def test_missing_cement_raises(self) -> None:
        """fluids 中无 lead/tail → ValueError。"""
        fluids = [_make_fluid("mud", 1400.0)]
        with self.assertRaises(ValueError):
            compute_muskat_regime(_make_result(), fluids, n_grid=201)

    def test_missing_viscosity_raises(self) -> None:
        """黏度字段均不可用 → ValueError。"""
        fluids = [
            _make_fluid("mud", 1400.0, mu_p=None, consistency_k=None),
            _make_fluid("lead", 1900.0, mu_p=0.05),
        ]
        with self.assertRaises(ValueError):
            compute_muskat_regime(_make_result(), fluids, n_grid=201)


if __name__ == "__main__":
    unittest.main()
