"""D2DGA 通量放大因子的独立计算函数。

本模块只放置与 D2DGA 通量修正公式直接相关的纯函数，避免把公式散落在
环空二维求解器内部，便于后续单独测试和复用。函数不读写外部状态，也不改变
输入数组。
"""

from __future__ import annotations

from typing import overload

import numpy as np
from numpy.typing import NDArray


Array = NDArray[np.float64]
FloatOrArray = float | Array


@overload
def d2dga_flux_amplification(
    cement_fraction: float,
    viscosity_ratio: FloatOrArray = 1.0,
    *,
    min_fraction: float = 0.01,
    max_fraction: float = 0.99,
    min_amplification: float = 0.5,
    max_amplification: float = 2.0,
) -> float: ...


@overload
def d2dga_flux_amplification(
    cement_fraction: Array,
    viscosity_ratio: FloatOrArray = 1.0,
    *,
    min_fraction: float = 0.01,
    max_fraction: float = 0.99,
    min_amplification: float = 0.5,
    max_amplification: float = 2.0,
) -> Array: ...


def d2dga_flux_amplification(
    cement_fraction: FloatOrArray,
    viscosity_ratio: FloatOrArray = 1.0,
    *,
    min_fraction: float = 0.01,
    max_fraction: float = 0.99,
    min_amplification: float = 0.5,
    max_amplification: float = 2.0,
) -> FloatOrArray:
    """计算 D2DGA 水泥相通量放大因子。

    公式沿用当前求解器中的 Zhang & Frigaard (2022) 口径：
    ``f(c, m) = [m*c² + 1.5*(1-c²)] / [m*c³ + (1-c³)]``。
    这里仅做数值安全裁剪，不引入任何 CBL、泥饼、温度、凝胶或湍流等工程修正。

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


def d2dga_dispersion_function_I3(
    c_bar: FloatOrArray,
    m: float = 1.0,
    *,
    min_fraction: float = 0.01,
    max_fraction: float = 0.99,
) -> FloatOrArray:
    """计算 D2DGA 浮力弥散函数 I3(ḉ, m)（Zhang & Frigaard 2022, 式 4.26）。

    公式：I3 = ḉ²(1-ḉ)³[4m·ḉ + 3(1-ḉ)] / {2m[m·ḉ³ + 1 - ḉ³]}

    性质：ḉ=0 或 ḉ=1 时 I3=0；ḉ≈0.5 附近达峰。用于 R2 浮力驱动弥散通量。

    Args:
        c_bar: 间隙平均水泥浓度（0~1），标量或数组。
        m: 黏度比 η_displaced/η_displacing。
        min_fraction: 计算前浓度下限，避免零浓度奇异。
        max_fraction: 计算前浓度上限，避免充满时奇异。
    """
    c = np.asarray(c_bar, dtype=float)
    c_safe = np.clip(c, min_fraction, max_fraction)
    c2 = c_safe ** 2
    c3 = c_safe ** 3
    one_minus_c = 1.0 - c_safe
    numerator = c2 * (one_minus_c ** 3) * (4.0 * m * c_safe + 3.0 * one_minus_c)
    denominator = 2.0 * m * (m * c3 + 1.0 - c3)
    i3 = numerator / denominator
    # 边界处置零（c=0 或 c=1 的精确值，clip 之外）
    i3 = np.where((c < min_fraction) | (c > max_fraction), 0.0, i3)
    if np.isscalar(c_bar):
        return float(i3)
    return i3.astype(float, copy=False)


def d2dga_buoyancy_flux(
    c_bar: FloatOrArray,
    m: float,
    delta_rho: float,
    H: FloatOrArray,
    eta2: float,
    f_phi: FloatOrArray,
    f_xi: FloatOrArray,
) -> tuple[FloatOrArray, FloatOrArray]:
    """计算 D2DGA 浮力驱动弥散通量 q_buoy（Zhang & Frigaard 2022, 式 4.25 第二项）。

    q_buoy = (Δρ H³ / (6 η2)) · I3(ḉ, m) · [-f_xi, f_phi]

    返回 (q_phi, q_xi)。
    """
    i3 = d2dga_dispersion_function_I3(c_bar, m)
    coef = (delta_rho * np.asarray(H, dtype=float) ** 3) / (6.0 * max(eta2, 1.0e-9))
    q_phi = coef * i3 * np.asarray(f_phi, dtype=float)
    q_xi = -coef * i3 * np.asarray(f_xi, dtype=float)
    return q_phi, q_xi
