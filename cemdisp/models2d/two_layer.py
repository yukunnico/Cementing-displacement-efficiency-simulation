"""两层牛顿闭包（Z&F22 §4.2）——I₁/I₂/I₃/q₀ 的收敛模块。

本模块把 Zhang & Frigaard (2022) 两层间隙闭包的文献公式集中在一处
（I₁/I₂/I₃ 自 ``d2dga_flux.py`` 迁移而来，迁移保逐位；q₀ 与
``layer_thickness_fraction`` 为新增接缝），供环空二维求解器与后续
Task 8/9（两层闭包接线）统一消费。函数均为纯函数：不读写外部状态，
不改变输入数组。

几何与符号约定（依 Z&F22 (4.11)–(4.21) 的推导口径）
--------------------------------------------------
- 间隙坐标 y∈[0, H]，y=0 为零剪应力面（两层解的应力 τ(0)=0），y=H 为壁面
  （无滑移，(4.19) 的分部积分要求 u(H)=0）；
- 流体 2 = 顶替液（水泥，黏度 η₂）占据 [0, y_i] 中线带；
  流体 1 = 被顶替液（黏度 η₁）占据 (y_i, H] 壁面带；
- c̄ = y_i/H 为间隙平均水泥体积分数；m = η₁/η₂ = 被顶替/顶替黏度比。

归一化约定
----------
- 默认形参（``eta1=eta2=H=1``）下 ``mobility_i1``/``mobility_i2`` 返回文献
  归一化闭式（Bararpour 2025 式 2.24/2.25；与迁移前 ``d2dga_dispersion_I1/I2``
  逐位一致）；
- 传入一致形参（``m = eta1/eta2``）时按 Z&F22 (4.21a)/(4.21b) 量纲式复现：
  代数恒等关系 I₁ = (H³/√(η₁η₂))·n₁(c̄,m)、I₂ = (H⁴/√(η₁η₂))·n₂(c̄,m)，
  其中 n₁/n₂ 为归一化闭式。

硬约束
------
- ``buoyancy_flux_distribution_i3`` **必须**用 Z&F22 (4.26)；B&F25 (2.27)
  的印刷式（括号 ``3(1−c̄²)``）已核实不自洽，禁止替换。

References
----------
Zhang & Frigaard (2022), JFM 947 A32：式 (4.21a,b)、(4.25)、(4.26)、(4.28)。
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]
FloatOrArray = float | Array


def mobility_i1(
    c_bar: FloatOrArray,
    m: float,
    eta1: float = 1.0,
    eta2: float = 1.0,
    H: float = 1.0,
) -> FloatOrArray:
    """牛顿平均流动度 I₁(c̄, m)（Z&F22 式 4.21a；归一化口径 = Bararpour 2025 式 2.24）。

    归一化闭式（默认 ``eta1=eta2=H=1``）::

        n₁ = [√m·c̄³ + (1−c̄³)/√m] / 3

    恒正，最小值 ≥ min(√m, 1/√m)/3。传入一致形参（``m = eta1/eta2``）时，
    乘以量纲换算因子 ``H³/√(eta1·eta2)`` 后严格等于式 (4.21a) 原文::

        I₁ = y_i³/(3η₂) + (H³−y_i³)/(3η₁)

    Args:
        c_bar: 间隙平均水泥浓度 c̄ = y_i/H ∈ [0,1]，标量或数组。
        m: 黏度比 η₁/η₂ = 被顶替液/顶替液（水泥）。
        eta1: 流体 1（被顶替液，壁面带）黏度。
        eta2: 流体 2（顶替液，中线带）黏度。
        H: 间隙特征高度（Z&F22 两层解的 y∈[0,H] 口径）。

    Returns:
        默认形参下与迁移前 ``d2dga_dispersion_I1`` 逐位一致；
        标量输入返回 ``float``，数组输入返回 ``Array``。
    """
    c = np.asarray(c_bar, dtype=float)
    sq_m = math.sqrt(m)
    out = (sq_m * c**3 + (1.0 - c**3) / sq_m) / 3.0
    # 量纲换算（Z&F22 (4.21a) ↔ 归一化闭式的代数恒等因子）；默认形参下因子=1.0，逐位无扰。
    out = out * (H**3 / math.sqrt(eta1 * eta2))
    return float(out) if np.isscalar(c_bar) else out.astype(float, copy=False)


def mobility_i2(
    c_bar: FloatOrArray,
    m: float,
    eta1: float = 1.0,
    eta2: float = 1.0,
    H: float = 1.0,
) -> FloatOrArray:
    """牛顿浮力流动度 I₂(c̄, m)（Z&F22 式 4.21b；归一化口径 = Bararpour 2025 式 2.25）。

    归一化闭式（默认 ``eta1=eta2=H=1``）::

        n₂ = [2√m·c̄³(1−c̄) + c̄(1−c̄)²(1+2c̄)/√m] / 6

    I₂(0)=I₂(1)=0，在混合区为正。传入一致形参（``m = eta1/eta2``）时，
    乘以量纲换算因子 ``H⁴/√(eta1·eta2)`` 后严格等于式 (4.21b) 原文::

        I₂ = (H−y_i)·y_i³/(3η₂) + y_i·(H−y_i)²·(H+2y_i)/(6η₁)

    Args:
        c_bar: 间隙平均水泥浓度 c̄ = y_i/H ∈ [0,1]，标量或数组。
        m: 黏度比 η₁/η₂ = 被顶替液/顶替液（水泥）。
        eta1: 流体 1（被顶替液，壁面带）黏度。
        eta2: 流体 2（顶替液，中线带）黏度。
        H: 间隙特征高度（Z&F22 两层解的 y∈[0,H] 口径）。

    Returns:
        默认形参下与迁移前 ``d2dga_dispersion_I2`` 逐位一致；
        标量输入返回 ``float``，数组输入返回 ``Array``。
    """
    c = np.asarray(c_bar, dtype=float)
    sq_m = math.sqrt(m)
    out = (2.0 * sq_m * c**3 * (1.0 - c) + c * (1.0 - c) ** 2 * (1.0 + 2.0 * c) / sq_m) / 6.0
    # 量纲换算（Z&F22 (4.21b) ↔ 归一化闭式的代数恒等因子）；默认形参下因子=1.0，逐位无扰。
    out = out * (H**4 / math.sqrt(eta1 * eta2))
    return float(out) if np.isscalar(c_bar) else out.astype(float, copy=False)


def buoyancy_flux_distribution_i3(
    c_bar: FloatOrArray,
    m: float = 1.0,
    *,
    min_fraction: float = 0.01,
    max_fraction: float = 0.99,
) -> FloatOrArray:
    """浮力弥散分布函数 I₃(c̄, m)（**Z&F22 式 4.26**，硬约束口径）。

    公式：I₃ = c̄²(1−c̄)³[4m·c̄ + 3(1−c̄)] / {2m[m·c̄³ + 1 − c̄³]}

    性质：c̄=0 或 c̄=1 时 I₃=0；c̄≈0.5 附近达峰。用于 (4.25) 第二项的
    浮力驱动弥散通量（见 ``d2dga_flux.d2dga_buoyancy_flux``）。

    ⚠️ 不得替换为 Bararpour 2025 (2.27) 的印刷式（括号 ``3(1−c̄²)``）——
    该式与 Z&F22/Z&F23/解析解三者均不自洽（2026-09-14 源文献审查裁定）。

    Args:
        c_bar: 间隙平均水泥浓度（0~1），标量或数组。
        m: 黏度比 η₁/η₂ = 被顶替液/顶替液。
        min_fraction: 计算前浓度下限，避免零浓度奇异。
        max_fraction: 计算前浓度上限，避免充满时奇异。

    Returns:
        标量输入返回 ``float``，数组输入返回 ``Array``。
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


def isotropic_flux_q0(c_bar: FloatOrArray, m: float) -> FloatOrArray:
    """各向同性通量函数 q₀(c̄, m)（Z&F22 式 (4.25) 第一项 / 式 (4.28)）。

    (4.25) 第一项 = 基准通量 × c̄ × f(c̄,m)，其中 f 为 (4.28) 放大因子；
    以基准通量归一后的顶替液（水泥）通量份额即::

        q₀ = c̄ · [m·c̄² + 1.5(1−c̄²)] / [m·c̄³ + 1 − c̄³]

    性质：q₀(0)=0、q₀(1)=1 精确成立；m=1 时 q₀ = c̄(1.5−0.5c̄²) =
    (3c̄−c̄³)/2，与抛物型剖面的中线带通量份额精确一致。c̄→0 时
    q₀ ≈ 1.5·c̄（水泥带位于剖面最快的中线区）。

    与 ``d2dga_flux.d2dga_flux_amplification`` 的关系：后者是 (4.28) 的
    裁剪包装（c̄ 裁到 [0.01,0.99] 后取 f）；本函数不做 c̄ 裁剪，以换取
    端点精确值，内部点两者满足 q₀ = c̄·f。两层实现分离是因为
    ``d2dga_flux`` 依赖本模块做 re-export，不允许反向依赖。

    Args:
        c_bar: 间隙平均水泥浓度 c̄ ∈ [0,1]，标量或数组（不裁剪，越界值
            的输出无物理意义，由调用方保证定义域）。
        m: 黏度比 η₁/η₂ = 被顶替液/顶替液，须 > 0（否则分母在 c̄→1 处退化）。

    Returns:
        标量输入返回 ``float``，数组输入返回 ``Array``。
    """
    c = np.asarray(c_bar, dtype=float)
    c2 = c ** 2
    c3 = c ** 3
    numerator = m * c2 + 1.5 * (1.0 - c2)
    denominator = m * c3 + (1.0 - c3)
    out = c * (numerator / denominator)
    return float(out) if np.isscalar(c_bar) else out.astype(float, copy=False)


def layer_thickness_fraction(c_bar) -> float:
    """界面厚度分数 c̄ = y_i/H（牛顿两层平界面的显式接缝）。

    牛顿两层闭包中界面是平的：流体 2（顶替液）体积分数恰等于界面位置
    分数 c̄ = y_i/H，故本函数为恒等映射。单设此接缝是为了给后续 HB
    （Herschel–Bulkley）推广留替换点——届时"浓度 → 界面厚度分数"不再
    是恒等式（屈服核/塞流区需另行求解），Task 8/9 只依赖本接口即可。

    Args:
        c_bar: 间隙平均水泥浓度，标量（plan 签名约定返回 ``float``）。

    Returns:
        界面厚度分数 y_i/H（= c̄）。
    """
    return float(c_bar)
