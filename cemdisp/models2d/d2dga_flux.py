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
    viscosity_ratio: float = 1.0,
    *,
    min_fraction: float = 0.01,
    max_fraction: float = 0.99,
    min_amplification: float = 0.5,
    max_amplification: float = 2.0,
) -> float: ...


@overload
def d2dga_flux_amplification(
    cement_fraction: Array,
    viscosity_ratio: float = 1.0,
    *,
    min_fraction: float = 0.01,
    max_fraction: float = 0.99,
    min_amplification: float = 0.5,
    max_amplification: float = 2.0,
) -> Array: ...


def d2dga_flux_amplification(
    cement_fraction: FloatOrArray,
    viscosity_ratio: float = 1.0,
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
        viscosity_ratio: 被顶替液/顶替液的黏度比 ``m``。
        min_fraction: 公式计算前的水泥体积分数下限，避免零浓度奇异行为。
        max_fraction: 公式计算前的水泥体积分数上限，避免完全充满时的数值尖点。
        min_amplification: 放大因子的下限裁剪值。
        max_amplification: 放大因子的上限裁剪值。

    Returns:
        与输入形状一致的通量放大因子；标量输入返回 ``float``。
    """

    # 将输入转为数组统一计算；copy=False 保持轻量，后续 clip 会生成安全结果。
    fraction = np.asarray(cement_fraction, dtype=float)
    c_safe = np.clip(fraction, min_fraction, max_fraction)
    numerator = viscosity_ratio * c_safe**2 + 1.5 * (1.0 - c_safe**2)
    denominator = viscosity_ratio * c_safe**3 + (1.0 - c_safe**3)
    amplification = np.clip(numerator / denominator, min_amplification, max_amplification)

    if np.isscalar(cement_fraction):
        return float(amplification)
    return amplification.astype(float, copy=False)
