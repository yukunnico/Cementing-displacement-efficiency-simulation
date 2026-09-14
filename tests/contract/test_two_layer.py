"""two_layer 两层牛顿闭包契约测试（Z&F22 式 4.21a/(4.21b)/(4.26)/(4.28)）。

测试内所有解析复算均为**独立转写**（不调用被测函数自身）：
- I3 参考式采用代数恒等变形（括号展开 3+c̄(4m−3)、分母 1+c̄³(m−1)）；
- I1/I2 归一化式按原运算顺序转写（用于逐位迁移验证）；
- I1/I2 量纲式按 Z&F22 (4.21a)/(4.21b) 原文直接复算（换算因子 H³/H⁴ ÷ √(η1η2)）；
- q0 = c̄·(4.28) 独立转写（(4.25) 第一项除以基准通量）。
"""
import math

import numpy as np
import pytest

from cemdisp.models2d.two_layer import (
    buoyancy_flux_distribution_i3,
    isotropic_flux_q0,
    layer_thickness_fraction,
    mobility_i1,
    mobility_i2,
)

# 简报规定的验证网格：m × c̄ 全组合（15 组）
M_GRID = (0.2, 0.5, 1.0, 2.0, 5.0)
C_GRID = (0.25, 0.5, 0.75)


# ---------------------------------------------------------------------------
# 独立解析复算（不调用被测函数）
# ---------------------------------------------------------------------------
def _i3_reference(c_bar: float, m: float) -> float:
    """Z&F22 (4.26) 独立复算：I3 = c̄²(1−c̄)³[4m·c̄+3(1−c̄)] / {2m[m·c̄³+1−c̄³]}。

    转写时做了恒等变形：4m·c̄+3(1−c̄) = 3 + c̄(4m−3)；m·c̄³+1−c̄³ = 1 + c̄³(m−1)。
    """
    num = c_bar**2 * (1.0 - c_bar) ** 3 * (3.0 + c_bar * (4.0 * m - 3.0))
    den = 2.0 * m * (1.0 + c_bar**3 * (m - 1.0))
    return num / den


def _i1_normalized_reference(c: np.ndarray, m: float) -> np.ndarray:
    """归一化 I1（Bararpour25 式 2.24 / Z&F22 (4.21a) 的 H³/√(η1η2) 归一），按原运算顺序转写。"""
    sq_m = math.sqrt(m)
    return (sq_m * c**3 + (1.0 - c**3) / sq_m) / 3.0


def _i2_normalized_reference(c: np.ndarray, m: float) -> np.ndarray:
    """归一化 I2（Bararpour25 式 2.25 / Z&F22 (4.21b) 的 H⁴/√(η1η2) 归一），按原运算顺序转写。"""
    sq_m = math.sqrt(m)
    return (2.0 * sq_m * c**3 * (1.0 - c) + c * (1.0 - c) ** 2 * (1.0 + 2.0 * c) / sq_m) / 6.0


def _i1_dimensional_reference(c_bar: float, m: float, eta1: float, eta2: float, H: float) -> float:
    """Z&F22 (4.21a) 原文量纲式直接复算：I1 = y_i³/(3η2) + (H³−y_i³)/(3η1)，y_i = c̄·H。"""
    y_i = c_bar * H
    return y_i**3 / (3.0 * eta2) + (H**3 - y_i**3) / (3.0 * eta1)


def _i2_dimensional_reference(c_bar: float, m: float, eta1: float, eta2: float, H: float) -> float:
    """Z&F22 (4.21b) 原文量纲式直接复算：I2 = (H−y_i)y_i³/(3η2) + y_i(H−y_i)²(H+2y_i)/(6η1)。"""
    y_i = c_bar * H
    return (H - y_i) * y_i**3 / (3.0 * eta2) + y_i * (H - y_i) ** 2 * (H + 2.0 * y_i) / (6.0 * eta1)


def _q0_reference(c_bar: float, m: float) -> float:
    """q0 独立复算：q0 = c̄·[m·c̄²+1.5(1−c̄²)] / [m·c̄³+1−c̄³]（(4.25) 第一项 / 基准通量，(4.28)）。"""
    num = m * c_bar**2 + 1.5 * (1.0 - c_bar**2)
    den = m * c_bar**3 + (1.0 - c_bar**3)
    return c_bar * num / den


# ---------------------------------------------------------------------------
# I3（Z&F22 式 4.26）——15 组合独立复算到 1e-12
# ---------------------------------------------------------------------------
class TestI3AgainstIndependentRecomputation:
    @pytest.mark.parametrize("m", M_GRID)
    @pytest.mark.parametrize("c_bar", C_GRID)
    def test_i3_matches_reference_1e12(self, m, c_bar):
        got = buoyancy_flux_distribution_i3(c_bar, m)
        ref = _i3_reference(c_bar, m)
        assert abs(got - ref) <= 1e-12, f"I3({c_bar}, {m}) 偏差 {abs(got - ref):.3e} 超过 1e-12"

    @pytest.mark.parametrize("m", M_GRID)
    def test_i3_grid_values_finite_and_positive(self, m):
        vals = np.array([buoyancy_flux_distribution_i3(c, m) for c in C_GRID])
        assert np.all(np.isfinite(vals))
        assert np.all(vals > 0.0)  # 混合区内 I3 恒正

    def test_i3_grid_anchor_value_m1_c05(self):
        # 手算锚点（沿用既有契约测试口径）：I3(0.5,1)=0.0546875
        assert buoyancy_flux_distribution_i3(0.5, 1.0) == pytest.approx(0.0546875, rel=1e-12)

    @pytest.mark.parametrize("m", M_GRID)
    def test_i3_endpoints_zero(self, m):
        # 边界置零行为随迁移保留：c̄=0 / c̄=1 精确为 0
        assert buoyancy_flux_distribution_i3(0.0, m) == 0.0
        assert buoyancy_flux_distribution_i3(1.0, m) == 0.0


# ---------------------------------------------------------------------------
# I1（Z&F22 式 4.21a）——牛顿退化 + 量纲换算
# ---------------------------------------------------------------------------
class TestMobilityI1:
    def test_newton_degeneracy_i1_0_equals_i1_1(self):
        # m=1（牛顿流体对）时 I1 退化为与 c̄ 无关的常数：I1(0)=I1(1)=1/3
        i1_0 = mobility_i1(0.0, 1.0)
        i1_1 = mobility_i1(1.0, 1.0)
        assert i1_0 == i1_1  # 端点严格相等（逐位）
        assert i1_0 == pytest.approx(1.0 / 3.0, abs=1e-15)

    @pytest.mark.parametrize("c_bar", (0.0,) + C_GRID + (1.0,))
    def test_newton_degeneracy_constant_in_c(self, c_bar):
        # m=1 时全 c̄ 域内 I1 ≡ 1/3（归一化口径）
        assert mobility_i1(c_bar, 1.0) == pytest.approx(1.0 / 3.0, abs=1e-15)

    @pytest.mark.parametrize("m", M_GRID)
    def test_normalized_matches_independent_transcription_bitwise(self, m):
        # 默认形参（eta1=eta2=H=1）下与归一化闭式独立转写**逐位**一致（迁移保真）
        c = np.array([0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0])
        got = mobility_i1(c, m)
        expected = _i1_normalized_reference(c, m)
        np.testing.assert_array_equal(got, expected)

    @pytest.mark.parametrize("eta1,eta2,H", [(2.0, 0.5, 2.0), (0.5, 2.0, 1.0), (1.7, 0.9, 0.05)])
    def test_dimensional_matches_zf22_4_21a(self, eta1, eta2, H):
        # 传入一致形参（m = eta1/eta2）时，结果须与 (4.21a) 原文量纲式独立复算一致
        m = eta1 / eta2
        for c_bar in C_GRID:
            got = mobility_i1(c_bar, m, eta1=eta1, eta2=eta2, H=H)
            ref = _i1_dimensional_reference(c_bar, m, eta1, eta2, H)
            assert got == pytest.approx(ref, rel=1e-12, abs=1e-15)


# ---------------------------------------------------------------------------
# I2（Z&F22 式 4.21b）——端点 + 量纲换算
# ---------------------------------------------------------------------------
class TestMobilityI2:
    @pytest.mark.parametrize("m", M_GRID)
    def test_endpoints_zero(self, m):
        # I2(0)=I2(1)=0（浮力流度在端点消失），随迁移保留
        assert mobility_i2(0.0, m) == 0.0
        assert mobility_i2(1.0, m) == 0.0

    @pytest.mark.parametrize("m", M_GRID)
    def test_normalized_matches_independent_transcription_bitwise(self, m):
        c = np.array([0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0])
        got = mobility_i2(c, m)
        expected = _i2_normalized_reference(c, m)
        np.testing.assert_array_equal(got, expected)

    @pytest.mark.parametrize("eta1,eta2,H", [(2.0, 0.5, 2.0), (0.5, 2.0, 1.0), (1.7, 0.9, 0.05)])
    def test_dimensional_matches_zf22_4_21b(self, eta1, eta2, H):
        m = eta1 / eta2
        for c_bar in C_GRID:
            got = mobility_i2(c_bar, m, eta1=eta1, eta2=eta2, H=H)
            ref = _i2_dimensional_reference(c_bar, m, eta1, eta2, H)
            assert got == pytest.approx(ref, rel=1e-12, abs=1e-15)


# ---------------------------------------------------------------------------
# q0（Z&F22 式 4.28 / (4.25) 第一项）
# ---------------------------------------------------------------------------
class TestIsotropicFluxQ0:
    @pytest.mark.parametrize("m", M_GRID)
    def test_q0_endpoints_exact(self, m):
        # q0(0)=0（无水泥无通量）、q0(1)=1（全水泥全通量），须精确成立
        assert isotropic_flux_q0(0.0, m) == 0.0
        assert isotropic_flux_q0(1.0, m) == 1.0

    @pytest.mark.parametrize("m", M_GRID)
    @pytest.mark.parametrize("c_bar", C_GRID)
    def test_q0_matches_reference_1e12(self, m, c_bar):
        got = isotropic_flux_q0(c_bar, m)
        ref = _q0_reference(c_bar, m)
        assert abs(got - ref) <= 1e-12

    @pytest.mark.parametrize("c_bar", C_GRID)
    def test_q0_m1_closed_form(self, c_bar):
        # m=1 闭式：q0 = c̄(1.5−0.5c̄²) = (3c̄−c̄³)/2（与抛物型剖面中心带份额精确一致）
        expected = (3.0 * c_bar - c_bar**3) / 2.0
        assert isotropic_flux_q0(c_bar, 1.0) == pytest.approx(expected, rel=1e-12)

    def test_q0_interior_consistent_with_legacy_amplification(self):
        # 内部点（[0.01,0.99]）与旧 d2dga_flux_amplification × c̄ 一致（同式 4.28）
        from cemdisp.models2d.d2dga_flux import d2dga_flux_amplification

        c = np.linspace(0.01, 0.99, 33)
        for m in M_GRID:
            got = isotropic_flux_q0(c, m)
            legacy = np.asarray(c, dtype=float) * d2dga_flux_amplification(c, m)
            np.testing.assert_allclose(got, legacy, rtol=1e-12, atol=0.0)


# ---------------------------------------------------------------------------
# layer_thickness_fraction（c̄ = y_i/H 接缝）
# ---------------------------------------------------------------------------
class TestLayerThicknessFraction:
    @pytest.mark.parametrize("c_bar", (0.0, 0.25, 0.5, 0.75, 1.0))
    def test_identity_mapping_returns_float(self, c_bar):
        # 牛顿两层平界面：流体2体积分数 = 界面位置分数 c̄ = y_i/H（恒等映射的显式接缝）
        out = layer_thickness_fraction(c_bar)
        assert isinstance(out, float)
        assert out == c_bar

    def test_returns_python_float_for_numpy_scalar(self):
        out = layer_thickness_fraction(np.float64(0.4))
        assert isinstance(out, float)
        assert out == 0.4


# ---------------------------------------------------------------------------
# 迁移保真：d2dga_flux 薄层 re-export + 标量/数组契约
# ---------------------------------------------------------------------------
class TestReExportIntegrity:
    def test_legacy_names_are_aliases_of_two_layer(self):
        import cemdisp.models2d.d2dga_flux as legacy

        assert legacy.d2dga_dispersion_I1 is mobility_i1
        assert legacy.d2dga_dispersion_I2 is mobility_i2
        assert legacy.d2dga_dispersion_function_I3 is buoyancy_flux_distribution_i3

    def test_scalar_input_returns_float(self):
        assert isinstance(mobility_i1(0.5, 2.0), float)
        assert isinstance(mobility_i2(0.5, 2.0), float)
        assert isinstance(buoyancy_flux_distribution_i3(0.5, 2.0), float)
        assert isinstance(isotropic_flux_q0(0.5, 2.0), float)

    def test_array_input_preserves_shape(self):
        c = np.array([0.2, 0.5, 0.8])
        assert mobility_i1(c, 2.0).shape == c.shape
        assert mobility_i2(c, 2.0).shape == c.shape
        assert buoyancy_flux_distribution_i3(c, 2.0).shape == c.shape
        assert isotropic_flux_q0(c, 2.0).shape == c.shape


# ---------------------------------------------------------------------------
# q₀ 参数化不变量（Task 9，controller 裁定第 4 条，2026-09-15）
# 随机 m ∈ [1e-3, 1e4] × c̄ 网格断言 0 ≤ q₀ ≤ 1——q₀ 若接线进任何通量层，
# 这是它的入口护栏；本 Task 未接线 q₀（半拉格朗日速度平流 + I₃ 通量层既有
# 口径不变，见 task-9-report），护栏先行固化。
# ---------------------------------------------------------------------------
class TestIsotropicFluxQ0Invariants:
    def test_q0_bounded_unit_interval_random_grid(self):
        rng = np.random.default_rng(20260915)
        m_values = np.exp(rng.uniform(np.log(1.0e-3), np.log(1.0e4), size=64))
        c_grid = np.linspace(0.0, 1.0, 201)
        for m in m_values:
            q0 = np.asarray(isotropic_flux_q0(c_grid, float(m)), dtype=float)
            assert np.all(q0 >= 0.0), f"m={m}: q₀ 出现负值（min={q0.min()}）"
            assert np.all(q0 <= 1.0), f"m={m}: q₀ 超过 1（max={q0.max()}）"

    def test_q0_endpoints_exact_all_m(self):
        # q₀(0)=0、q₀(1)=1 对任意 m 精确成立（two_layer docstring 契约）
        for m in (1.0e-3, 0.2, 1.0, 5.0, 1.0e4):
            assert isotropic_flux_q0(0.0, m) == 0.0
            assert isotropic_flux_q0(1.0, m) == 1.0

    def test_q0_equals_c_bar_times_flux_amplification_on_interior(self):
        # 内点 q₀ = c̄·f（d2dga_flux_amplification 的未裁剪恒等关系）
        from cemdisp.models2d.d2dga_flux import d2dga_flux_amplification

        rng = np.random.default_rng(915)
        for m in (0.2, 1.0, 5.0):
            c = rng.uniform(0.02, 0.98, size=32)
            f = np.asarray(d2dga_flux_amplification(c, m), dtype=float)
            q0 = np.asarray(isotropic_flux_q0(c, m), dtype=float)
            assert np.allclose(q0, c * f, rtol=1e-12)
