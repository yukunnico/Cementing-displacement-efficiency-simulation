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
- ``mobility_i1``/``mobility_i2``/``isotropic_flux_q0`` 为牛顿闭式（阶段 B 逐位
  护栏的锚），HB 相关改动只做**纯增量**（新增 :func:`hb_groups`）。

非牛顿增量与适用域
------------------
:func:`hb_groups` 实现 B&F25 (2.29)—(2.34) 的 HB 无量纲群换算（供论文报告与制表
索引；真实闭包值由 ``hb_closure.HBClosure`` 在注入的局部 G 上现算）。该推广属
**扩展应用，未获外部验证**；Z&F22/23 的判据与结论仅严格适用于**竖直井 + 牛顿
流体**，不得引其作 HB 依据（见 ``docs/源模型口径与适用域声明.md``）。

References
----------
Zhang & Frigaard (2022), JFM 947 A32：式 (4.21a,b)、(4.25)、(4.26)、(4.28)。
Bararpour & Frigaard (2025), JFM 1022 A15：(2.29)—(2.34)（``hb_groups``）。
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


def hb_groups(
    kappa1: float,
    kappa2: float,
    n1: float,
    n2: float,
    tauY1: float,
    tauY2: float,
    gamma0: float,
) -> dict:
    """HB 无量纲群换算（Bararpour & Frigaard 2025 式 (2.29)—(2.34)）。

    下标约定与 :func:`mobility_i1` 一致：``1`` = 被顶替液（壁面带）、``2`` = 顶替液
    （水泥，中线带）；``κ̂_k/n_k/τ̂_{Y,k}`` 为量纲稠度系数/幂律指数/屈服应力。

    公式（逐式对照文献，``γ̇₀`` 为**调用方给定的参考剪切率**）::

        μ̂e = [κ̂₁γ̇₀^{n₁−1}·κ̂₂γ̇₀^{n₂−1}]^{1/2}      (2.29)
        τ̂₀ = μ̂eγ̇₀ + max{τ̂_Y,1, τ̂_Y,2}                (2.30)
        m  = κ̂₁γ̇₀^{n₁−n₂}/κ̂₂                            (2.32)
        B  = max{τ̂_Y,1, τ̂_Y,2}/(μ̂eγ̇₀)                    (2.33)
        κ_k = κ̂_kγ̇₀^{n_k}/τ̂₀ = m^{1/2}/(1+B)           (2.32)
        τ_{Y,k} = τ̂_{Y,k}/τ̂₀                              (2.32)

    ``μ̂e`` 是 γ̇₀ 处的有效黏度尺度（= ``√(η₁η₂)``，``η_k = κ̂_kγ̇₀^{n_k−1}`` 为表观
    黏度），故 ``m`` 就是 γ̇₀ 处的**表观黏度比**——与牛顿 ``η₁/η₂`` 同一意义的推广。
    用 ``max{τ̂_Y}``（而非逐相）是为了让缩放屈服应力落在 ``[0,1)``（``max τ_{Y,k}
    = B/(1+B)``），便于数值求值与制表。

    ⚠️ **``γ̇₀`` 的选择决定整组数**（文献口径 ``γ̇₀ = ŵ₀/d̂*``：速度尺度/长度尺度），
    且 ``B`` 里的 ``γ̇₀`` 与**局部应力尺度**绑定（``B`` 大 ⇒ 屈服项重要）。本函数
    只做代数换算、不选 ``γ̇₀``；D2DGA 的真实闭包值**不**由这些群决定，而由
    ``hb_closure.HBClosure`` 在**注入的局部 G** 上现算（``(c̄, B)`` 表只在固定 ``G``
    下自洽）。本组数的用途是**论文报告与制表索引**。
    ⚠️ 本推广属**扩展应用，未获外部验证**（B&F25 自述无外部验证）；Z&F22/23 的判据
    与结论仅严格适用于**竖直井 + 牛顿流体**，不得引其作 HB 依据。

    Args:
        kappa1, kappa2: 稠度系数 κ̂₁/κ̂₂（> 0）。
        n1, n2: 幂律指数 n₁/n₂（> 0）。
        tauY1, tauY2: 屈服应力 τ̂_{Y,1}/τ̂_{Y,2}（≥ 0）。
        gamma0: 参考剪切率 γ̇₀（> 0），口径见上（文献 ``ŵ₀/d̂*``）。

    Returns:
        dict：``mu_e``（μ̂e）、``tau_0``（τ̂₀）、``m``、``B``、``kappa_k``（2-元组，
        逐相缩放稠度系数）、``tau_Yk``（2-元组，逐相缩放屈服应力）。
        牛顿退化 ``n₁=n₂=1, τ̂_Y≡0`` ⇒ ``B=0``、``m=κ̂₁/κ̂₂``、
        ``kappa_k=(√m, 1/√m)``、``tau_Yk=(0, 0)``。

    Raises:
        ValueError: κ̂ ≤ 0、n ≤ 0、τ̂_Y < 0 或 γ̇₀ ≤ 0（含 NaN/非有限值）。
    """
    for name, val, positive in (("kappa1", kappa1, True), ("kappa2", kappa2, True),
                                ("n1", n1, True), ("n2", n2, True),
                                ("gamma0", gamma0, True)):
        v = float(val)
        if not (math.isfinite(v) and (v > 0.0 if positive else v >= 0.0)):
            raise ValueError(f"HB 群换算：{name} 须为正的有限值，得到 {val!r}")
    tau_max = 0.0
    for name, val in (("tauY1", tauY1), ("tauY2", tauY2)):
        v = float(val)
        if not (math.isfinite(v) and v >= 0.0):
            raise ValueError(f"HB 群换算：{name} 须为非负的有限值，得到 {val!r}")
        tau_max = max(tau_max, v)

    k1, k2, e1, e2, g = float(kappa1), float(kappa2), float(n1), float(n2), float(gamma0)
    # (2.29)：μ̂e 用表观黏度 η_k = κ̂_kγ̇₀^{n_k−1} 的几何平均（等价于文献乘积式）
    mu_e = math.sqrt((k1 * g ** (e1 - 1.0)) * (k2 * g ** (e2 - 1.0)))
    tau_0 = mu_e * g + tau_max                                   # (2.30)
    m = k1 * g ** (e1 - e2) / k2                                 # (2.32)
    B = tau_max / (mu_e * g)                                     # (2.33)
    return {
        "mu_e": mu_e,
        "tau_0": tau_0,
        "m": m,
        "B": B,
        "kappa_k": (k1 * g ** e1 / tau_0, k2 * g ** e2 / tau_0),  # (2.32) 左式
        "tau_Yk": (float(tauY1) / tau_0, float(tauY2) / tau_0),   # (2.32) 末式
    }
