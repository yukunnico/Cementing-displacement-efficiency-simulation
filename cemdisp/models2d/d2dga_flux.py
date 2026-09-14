"""D2DGA 通量闭包的兼容薄层（re-export）。

I₁/I₂/I₃ 两层牛顿闭包已迁移至 :mod:`cemdisp.models2d.two_layer`（Task 6，
迁移保逐位），本模块保留历史函数名作为**同一函数对象的别名**，避免破坏
既有 import（``annulus_d2dga.py`` / ``muskat_regime.py`` /
``zhang2022_benchmark.py`` / ``tests/contract/test_d2dga_flux.py``）。

仍在本模块实现的两个函数：
- ``d2dga_flux_amplification``：式 (4.28) 放大因子的裁剪包装（c̄ 裁到
  [0.01,0.99] 后取 f；未裁剪的精确端点口径见 ``two_layer.isotropic_flux_q0``）；
- ``d2dga_buoyancy_flux``：式 (4.25) 第二项浮力弥散通量的组装。
"""

from __future__ import annotations

import math
from typing import overload

import numpy as np
from numpy.typing import NDArray

from cemdisp.models2d.two_layer import (
    buoyancy_flux_distribution_i3,
    mobility_i1,
    mobility_i2,
)

Array = NDArray[np.float64]
FloatOrArray = float | Array

# 历史函数名 -> two_layer 新名（同一函数对象，数值行为逐位一致）
d2dga_dispersion_I1 = mobility_i1
d2dga_dispersion_I2 = mobility_i2
d2dga_dispersion_function_I3 = buoyancy_flux_distribution_i3


@overload
def d2dga_flux_amplification(
    cement_fraction: float,
    viscosity_ratio: FloatOrArray = 1.0,
    *,
    min_fraction: float = 0.01,
    max_fraction: float = 0.99,
    min_amplification: float = -np.inf,
    max_amplification: float = np.inf,
) -> float: ...


@overload
def d2dga_flux_amplification(
    cement_fraction: Array,
    viscosity_ratio: FloatOrArray = 1.0,
    *,
    min_fraction: float = 0.01,
    max_fraction: float = 0.99,
    min_amplification: float = -np.inf,
    max_amplification: float = np.inf,
) -> Array: ...


def d2dga_flux_amplification(
    cement_fraction: FloatOrArray,
    viscosity_ratio: FloatOrArray = 1.0,
    *,
    min_fraction: float = 0.01,
    max_fraction: float = 0.99,
    min_amplification: float = -np.inf,
    max_amplification: float = np.inf,
) -> FloatOrArray:
    """计算 D2DGA 水泥相通量放大因子（Zhang & Frigaard 2022, 式 4.28）。

    公式沿用当前求解器中的 Zhang & Frigaard (2022) 口径：
    ``f(c, m) = [m*c² + 1.5*(1-c²)] / [m*c³ + (1-c³)]``。
    这里仅做数值安全裁剪，不引入任何 CBL、泥饼、温度、凝胶或湍流等工程修正。
    未裁剪、端点精确的通量函数口径（q₀ = c̄·f）见
    :func:`cemdisp.models2d.two_layer.isotropic_flux_q0`。

    Args:
        cement_fraction: 水泥相局部体积分数，可以是标量或 NumPy 数组。
        viscosity_ratio: 被顶替液/顶替液的黏度比 ``m``，标量或与 ``cement_fraction``
            形状兼容的数组（用于 R1：随空间变化的 m 场）。
        min_fraction: 公式计算前的水泥体积分数下限，避免零浓度奇异行为。
        max_fraction: 公式计算前的水泥体积分数上限，避免完全充满时的数值尖点。
        min_amplification: 放大因子的下限裁剪值。
        max_amplification: 放大因子的上限裁剪值。

    Returns:
        与输入形状一致的通量放大因子；标量输入返回 ``float``。
    """

    # 将输入转为数组统一计算；copy=False 保持轻量，后续 clip 会生成安全结果。
    fraction = np.asarray(cement_fraction, dtype=float)
    m = np.asarray(viscosity_ratio, dtype=float)
    c_safe = np.clip(fraction, min_fraction, max_fraction)
    numerator = m * c_safe**2 + 1.5 * (1.0 - c_safe**2)
    denominator = m * c_safe**3 + (1.0 - c_safe**3)
    amplification = np.clip(numerator / denominator, min_amplification, max_amplification)

    if amplification.ndim == 0:
        return float(amplification)
    return amplification.astype(float, copy=False)


def d2dga_buoyancy_flux(
    c_bar: FloatOrArray,
    m: float,
    delta_rho: FloatOrArray,
    H: FloatOrArray,
    eta2: FloatOrArray,
    f_phi: FloatOrArray,
    f_xi: FloatOrArray,
) -> tuple[FloatOrArray, FloatOrArray]:
    """计算 D2DGA 浮力驱动弥散通量 q_buoy（Zhang & Frigaard 2022, 式 4.25 第二项）。

    q_buoy = (Δρ H³ / (6 η2)) · I3(ḉ, m) · [-f_xi, f_phi]

    delta_rho（密度差, kg/m³）与 eta2（顶替液黏度）均可为标量或数组；
    数组时与 H/c_bar/f_phi/f_xi 按广播规则对齐（局部化 I3）。

    返回 (q_phi, q_xi)。
    """
    i3 = d2dga_dispersion_function_I3(c_bar, m)
    coef = (delta_rho * np.asarray(H, dtype=float) ** 3) / (6.0 * np.maximum(eta2, 1.0e-9))
    q_phi = coef * i3 * np.asarray(f_phi, dtype=float)
    q_xi = -coef * i3 * np.asarray(f_xi, dtype=float)
    return q_phi, q_xi
