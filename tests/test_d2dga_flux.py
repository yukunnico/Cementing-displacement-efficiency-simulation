"""D2DGA 通量纯函数测试。覆盖：式4.28放大因子（含数组m）+ 式4.26 I₃弥散函数。"""
import numpy as np
import pytest
from cemdisp.models2d.d2dga_flux import d2dga_flux_amplification, d2dga_dispersion_function_I3


class TestFluxAmplificationArrayM:
    def test_scalar_m_equals_1_at_c05_returns_1p375(self):
        # m=1, c=0.5: f(0.5,1) = [0.25 + 1.5*0.75] / [0.125 + 0.875] = 1.375/1.0 = 1.375
        # NOTE: brief's original "returns 1.5 constant" is a math error — the formula
        # gives 1.5 only at c=0; f(c,1) decreases monotonically from 1.5 to 1.0 as c goes 0->1.
        c = np.array([0.5])
        f = d2dga_flux_amplification(c, viscosity_ratio=1.0)
        assert np.allclose(f, 1.375, rtol=1e-4)

    def test_array_m_different_from_scalar(self):
        # 同一 c 场，不同 m 应给不同 f
        c = np.array([0.3, 0.5, 0.7])
        m_array = np.array([0.5, 1.0, 2.0])
        f_array = d2dga_flux_amplification(c, viscosity_ratio=m_array)
        f_scalar = d2dga_flux_amplification(c, viscosity_ratio=1.0)
        assert not np.allclose(f_array, f_scalar)
        # 小 m（顶替液更粘）应使放大因子更接近 1（弥散更弱）在高浓度区
        # m=0.5 在 c=0.9 处 f 应 < m=1 的 1.5
        f_m05 = d2dga_flux_amplification(np.array([0.9]), viscosity_ratio=0.5)
        assert f_m05[0] < 1.5

    def test_output_shape_matches_input(self):
        c = np.linspace(0.01, 0.99, 40)
        m = np.full(40, 0.8)
        f = d2dga_flux_amplification(c, viscosity_ratio=m)
        assert f.shape == c.shape


class TestDispersionFunctionI3:
    def test_zero_concentration_returns_zero(self):
        # c=0 -> I3 = 0 (分子含 c²)
        assert d2dga_dispersion_function_I3(0.0, m=1.0) == pytest.approx(0.0, abs=1e-9)

    def test_full_concentration_returns_zero(self):
        # c=1 -> I3 = 0 (分子含 (1-c)³)
        assert d2dga_dispersion_function_I3(1.0, m=1.0) == pytest.approx(0.0, abs=1e-9)

    def test_mid_concentration_positive(self):
        # c=0.5 时 I3 应为正且较大（峰值附近）
        i3 = d2dga_dispersion_function_I3(0.5, m=1.0)
        assert i3 > 0.0
        # 手算核对 m=1, c=0.5:
        # 分子=0.25*0.125*(4*0.5+3*0.5)=0.25*0.125*3.5=0.109375
        # 分母=2*1*(1*0.125+1-0.125)=2*1.0=2 -> I3=0.0546875
        # NOTE: brief's "0.04296875" came from a 2.75 vs 3.5 math slip.
        assert i3 == pytest.approx(0.0546875, rel=1e-4)

    def test_zero_density_implied_zero_flux_handled_separately(self):
        # I3 函数本身不含 Δρ；Δρ=0 时通量在调用方置零。这里只验 I3 数值正确。
        c = np.array([0.2, 0.5, 0.8])
        i3 = d2dga_dispersion_function_I3(c, m=1.0)
        assert i3.shape == c.shape
        assert np.all(i3 >= 0.0)  # I3 在 [0,1] 非负


class TestScalarCementArrayM:
    def test_scalar_cement_array_m_returns_array(self):
        # scalar cement + array m must NOT crash; returns array shaped like m
        c = 0.5
        m = np.array([0.5, 1.0, 2.0])
        f = d2dga_flux_amplification(c, viscosity_ratio=m)
        assert isinstance(f, np.ndarray)
        assert f.shape == m.shape

    def test_scalar_cement_scalar_m_returns_float(self):
        f = d2dga_flux_amplification(0.5, viscosity_ratio=1.0)
        assert isinstance(f, float)
        assert f == pytest.approx(1.375)
