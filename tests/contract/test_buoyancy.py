"""浮力口径统一模块契约测试。

守护 `cemdisp.models2d.buoyancy` 的四条规格（Task 3）：
1. 浮力数 b 的定义与符号（Z&F22 p.8）；
2. 密度倒置 → b<0；
3. 幂律/HB 表观黏度用 K·γ̇^(n−1)，**不得静默回退**到硬编码 0.05；
4. Froude 数平方的量级（Z&F22 (2.6)），不是硬编码 1.0。
"""

import pytest

from cemdisp.models2d.buoyancy import (buoyancy_number, froude_squared,
                                       fluid_apparent_viscosity)
from cemdisp.data.fluid_spec import FluidSpec, FluidRole, RheologyModel


def test_buoyancy_number_formula_and_sign():
    b = buoyancy_number(2100.0, 2020.0, 0.0257, 0.066, 0.563)
    assert b == pytest.approx(80 * 9.81 * 0.0257 ** 2 / (0.066 * 0.563), rel=1e-12)
    assert b > 0


def test_density_inversion_is_negative():
    assert buoyancy_number(1900.0, 1960.0, 0.045, 0.058, 0.764) < 0


def test_power_law_viscosity_not_silent_fallback():
    pl = FluidSpec(name="泥浆", role=FluidRole.MUD, density_kg_m3=1960.0,
                   rheology_model=RheologyModel.POWER_LAW,
                   power_law_n=0.72, consistency_k=0.381)
    v = fluid_apparent_viscosity(pl, shear_rate=10.0)
    assert v == pytest.approx(0.381 * 10.0 ** (0.72 - 1.0), rel=1e-12)
    assert v != 0.05


def test_froude_squared_matches_manual_scale():
    f2 = froude_squared(mu_displaced=0.058, w0_mps=0.764, half_gap_m=0.0458,
                        rho_displaced=1200.0, gap_scale_m=0.0917,
                        mean_radius_m=0.1071)
    assert 1e-4 < f2 < 1.0     # 实测量级 ~4.7e-3，绝不是 1.0
