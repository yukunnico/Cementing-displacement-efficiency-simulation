"""流变闭包提供者协议（Phase B-1）。

把 ``stream_function`` 求解器与具体闭包解耦：求解器只依赖 :class:`ClosureProvider`
协议。:class:`NewtonianClosure` 逐位等于 ``two_layer.mobility_i1/i2``
（Z&F22 (4.21a)/(4.21b)）。Phase A 的 HB 查表闭包（``HBClosure``）实现同一协议，
额外依赖局部压力梯度 G（由外层迭代注入）。

适用域提醒（务必随论文口径一致）
--------------------------------
源模型 Z&F22/23 的判据与结论**仅严格适用于竖直井 + 牛顿流体**；本协议的任何
幂律/HB 实现均属**扩展应用，未获外部验证**。见
``docs/源模型口径与适用域声明.md`` §1 声明 3。

References
----------
Zhang & Frigaard (2022), *JFM* **947**, A32：(4.21a)/(4.21b)。
Bararpour & Frigaard (2025), *JFM* **1022**, A15：(2.14)/(2.15) HB 推广。
"""

from __future__ import annotations

from typing import Protocol, runtime_checkable

import numpy as np
from numpy.typing import NDArray

from cemdisp.models2d.two_layer import mobility_i1, mobility_i2

Array = NDArray[np.float64]


@runtime_checkable
class ClosureProvider(Protocol):
    """两层闭包协议：返回间隙平均流动度 I₁ 与浮力流动度 I₂。

    约定（Z&F22 (4.21a)/(4.21b) 口径）：``c_bar`` 为 (ny,nz) 场或标量；
    ``m``/``eta1``/``eta2`` 为标量；``H`` 为 (ny,nz) 场或标量；返回与
    ``c_bar``/``H`` 广播后的 float 数组。
    """

    def mobility(self, c_bar, m: float, eta1: float, eta2: float, H) -> Array:
        """间隙平均流动度 I₁(c̄, m; η₁, η₂, H)（式 (4.21a)）。"""
        ...

    def buoyant_mobility(self, c_bar, m: float, eta1: float, eta2: float, H) -> Array:
        """浮力流动度 I₂(c̄, m; η₁, η₂, H)（式 (4.21b)）。"""
        ...


class NewtonianClosure:
    """牛顿两层闭包（Z&F22 (4.21a)/(4.21b)）——逐位等于 ``two_layer`` 原函数。

    本类不做任何数值加工，仅把 ``mobility_i1``/``mobility_i2`` 适配到
    :class:`ClosureProvider` 协议。它的存在保证「closure 缺省注入」与 HEAD
    逐位一致（B-1 的 L1 硬约束）。
    """

    def mobility(self, c_bar, m: float, eta1: float, eta2: float, H) -> Array:
        return np.asarray(mobility_i1(c_bar, m, eta1=eta1, eta2=eta2, H=H), dtype=float)

    def buoyant_mobility(self, c_bar, m: float, eta1: float, eta2: float, H) -> Array:
        return np.asarray(mobility_i2(c_bar, m, eta1=eta1, eta2=eta2, H=H), dtype=float)


class PowerLawGapClosure:
    """牛顿闭包 + 幂律间隙一阶修正（B-3，**近似**，非 B&F25 闭包口径）。

    I₁ → I₁ · (H/H̄)^{1/n − 1}（H̄ = 全场均值）。依据：单流体幂律槽流
    ū ∝ H^{1+1/n}G^{1/n} ⇒ I₁ = Hū/G ∝ H^{2+1/n}G^{1/n−1}；相对牛顿 H³ 的
    场修正即 (H/H̄)^{1/n−1}。G 的纯标量因子在椭圆解 flux 归一化中消去，
    故只有 H 依赖进入算子形状。n=1 时因子恒为 1（逐位退化为牛顿）。

    ⚠️ 适用域：一阶近似，不得引 B&F25 作为方法学依据（见设计规格 §1.5）。
    """

    def __init__(self, n: float, base: "ClosureProvider | None" = None) -> None:
        self._n = float(n)
        if not (self._n > 0.0):
            # n ≤ 0（含 0.0/NaN）会使 1/n 无定义或非物理——构造即拒绝。
            raise ValueError(f"幂律指数 n 必须为正（1/n 需有定义），得到 n={n!r}")
        self._base = base if base is not None else NewtonianClosure()

    def mobility(self, c_bar, m: float, eta1: float, eta2: float, H) -> Array:
        I1 = np.asarray(self._base.mobility(c_bar, m, eta1, eta2, H), dtype=float)
        if abs(self._n - 1.0) < 1e-12:
            return I1
        Harr = np.asarray(H, dtype=float)
        Hbar = float(np.mean(Harr))
        return I1 * (Harr / Hbar) ** (1.0 / self._n - 1.0)

    def buoyant_mobility(self, c_bar, m: float, eta1: float, eta2: float, H) -> Array:
        """浮力流动度 I₂ —— **故意不修正**。

        本类只覆盖 I₁ 的 H 依赖（I₁ ∝ H^{2+1/n}）；I₂ = H⁴ 的幂律标度未推导，
        故直接委托 base（牛顿闭包）。当前 solver 只调 :meth:`mobility`（I₂ 经
        ``χ`` 的 ``i2_field`` 另有牛顿闭式路径消费），该口子暂无影响；若未来
        I₂ 也进算子，须先补幂律标度推导。
        """
        return self._base.buoyant_mobility(c_bar, m, eta1, eta2, H)
