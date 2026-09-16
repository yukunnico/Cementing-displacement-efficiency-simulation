"""一维间隙 Herschel–Bulkley 弱解求解器（B&F25 附录 A；A-1a/A-1c）。

``solve_fixed_G`` 求解附录 A.2.1 的两流体层 1D 间隙问题（增广拉格朗日 / Uzawa）；
``solve_fixed_mean_velocity`` 求解附录 A.2.2 的**定均速**变体（给 ``ū*`` 反求 ``G``，
D2DGA 外迭代实际需要的入口）。两者都返回 (2.13)—(2.15) 定义的闭包量 ``I₁ / I₂``
与 ``q₀``。它们是 Phase A 的 ``HBClosure``（Task 4）的底层求解器。

坐标与层约定
------------
- 内部在归一化坐标 ỹ∈[0,1] 上求解（附录 A 的 tilde 系，H=1）：ỹ=0 为零剪应力面
  （对称面，τ̃(0)=0），ỹ=1 为壁面（无滑移 ũ(1)=0）；界面在 ỹ = c̄。
- 流体 2 = 顶替液（中线带 ỹ<c̄），用 ``n[1]/kappa[1]/tau_y[1]``；流体 1 =
  被顶替液（壁面带 ỹ>c̄），用下标 ``0``——与 ``two_layer`` 模块 docstring 一致。
- ``H`` 为物理半隙；由 (A4) 换到 tilde 系：``κ̃_k = κ_k/H^{n_k}``、``G̃ = H·G``、
  ``τ̃_{Y,k} = τ_{Y,k}``（应力与速度在 (A4) 下与量纲值同值）。

返回值口径（**量纲**，与 ``two_layer.mobility_i1/mobility_i2`` 同口径）
---------------------------------------------------------------------
- ``I1 = H²·Ĩ₁``、``I2 = H³·Ĩ₂``（Ĩ 为 tilde 系下 (2.14)/(2.15) 的积分）；
  ``H=1`` 时与 Z&F22 (4.21a)/(4.21b) 逐项一致（见 ``two_layer.mobility_i1``
  docstring 的量纲换算段）。
- ``q0`` 无量纲，按 (2.22) 归一：``q0(0)=0``、``q0(1)=1``；本实现按**通量比**定义
  （见下"q₀ 口径"段）。
- ``u``/``y``：同一套格边网格（都长 ``ny``，``y = linspace(0,1,ny)``）；``u`` 是
  tilde 系速度 ũ()，即量纲速度在归一化位置处的取值。
- **``u`` 的形状契约**：轴向/共线输入（``G`` 与 ``Gb`` 共线，含全零 ``Gb``）下 ``u``
  是长 ``ny`` 的 1 维数组（这是 A.2.1/A.2.2 的标量情形，与 A-1a 起的行为逐位一致）；
  **非共线**输入下问题是 2 维向量场，``u`` 是 ``(2, ny)`` 数组（第 0/1 分量 = 输入
  坐标系的 φ/ξ 分量）。
- **ū 的口径是 :attr:`GapSolution.u_bar`（``np.trapezoid(u, y)``），不是 ``u.mean()``**：
  格边网格含 ỹ=0 与 ỹ=1 两端点，``mean`` 带 O(h) 偏差（``ny=201`` 时实测 ~1.3e-3），
  ``trapezoid`` 为 O(h²)（~2.1e-5）。全项目只认 ``u_bar``。

q₀ 口径（通量比，**取代** B&F25 印刷式 (2.22)）
---------------------------------------------
``q₀ = |∫₀^{c̄} ũ dỹ| / |∫^1 ũ dỹ|``（逐分量求积后取模之比）——通量比就是 q₀ 的定义。
实现与精度见 :func:`_flux_ratio_q0`；``Gb=0`` 时退回 Z&F22 等密度闭式
``isotropic_flux_q0``（实测 rel≤2.4e-16），是本改动的回归锚。

⚠️ B&F25 **印刷式 (2.22)** 与本口径不一致（``m=3, c̄=0.5`` 牛顿：本口径/解析槽流
= 0.75，印刷式 = 0.9），且其分子用 ``y_i·∫ỹ/η₁`` 一类与 (2.14) 不同族的权重、对
``Gb≠0`` 结构上不足（``q₀`` 实为**各层厚度**的函数）⇒ 不采用印刷式。

实现要点（A.2.3 交错网格）
--------------------------
- ``ũ`` 在格边（``ny`` 点）、``λ/q/c̄`` 在格心（``ny−1`` 点）；格心处 (A18) 的
  右端为 ``dũ^{k+1}/dỹ = q^k − λ̃^k/r``（对称面在 ỹ=0 ⇒ 积分常数 C=0，见报告）；
  ``ũ`` 由壁面 ũ(1)=0 向回累积（中点求积，对线性 γ̇ 精确）。
- ``Ĩ₁/Ĩ₂/q₀`` 由 (2.14)/(2.15)/(2.22) 在**解析应力场** τ̃(ỹ)=λ₀(ỹ)（(A1)/(A2)/
  (A17)，逐层仿射）上 quadrature 求得：按屈服面切段后逐段 Gauss–Legendre，故
  界面/屈服面处无离散误差（牛顿锚可命中 rel=1e-10，与网格 ``ny`` 无关）。
- Uzawa 固定 ``ρ=r=1``（A.2.3，满足 (A31)）；``θ`` 由 (A20) 二分解出；判据为
  ``‖Δũ‖_p, ‖Δq‖_p, ‖λ̃‖_p < tol``，``p = 1+min(n₁,n₂)``；达 ``max_iter`` 仍未
  收敛**抛 RuntimeError**（不静默回退）。

向量版（R-T1-3，解除共线限制）
------------------------------
B&F25 的 ``G``（(2.11)）与 ``Gb``（(2.12)）本身是 **2 维向量**，附录 A 只写了共线
标量式。本实现的推广（逐层应力是 2 维**仿射**场 ``τ̃(ỹ)=c·ỹ+d``）：

- 中线带（流体 2，∈[0,c̄]）：``c=G̃−(1−c̄)Gb̃``、``d=0``（τ̃(0)=0 对称面）；
- 壁面带（流体 1，∈[c̄,1]）：``c=G̃+c̄Gb̃``、``d=(c₂−c₁)c̄``（界面剪应力连续）。
  共线时 ``d ∥ c``（可归约为有符号标量）；**非共线时 ``d`` 不平行于 ``c``**——
  这正是标量式跑不了的原因。
- 关键恒等式：``|τ̃(ỹ)|`` 仍是 ỹ 的**标量**函数，且 HB 本构各向同性（``η_eff`` 只
  依赖 ``|τ̃|``）⇒ 按屈服面/零应力面切段 + 逐段 Gauss–Legendre 的求积结构直接复用，
  只是切点改为 ``|c·ỹ+d|`` 的零点（一般化为二次方程）。
- Uzawa 侧：``η_eff`` 各向同性 ⇒ (A18) 的算子逐分量同一系数 ⇒ ``ũ``/``q`` 逐分量求
  解；两组分唯一的耦合处是 (A19)/(A29) 的 ``|m|``（向量模）。``_bisect_theta`` 因此
  按逐格**向量模**工作。

⚠️ 共线输入：（i）先做"旋转到共同方向"的**精确**归约（问题各向同性），（ii）归约后
第 2 分量恒为 0，故所有模/切点走与原标量实现逐式相同的分支 ⇒ ``I₁/I₂/q₀/iters/u``
与向量化改造前**逐位一致**（回归基线见 ``.tmp_research/task3_baseline/``）。若让
``hypot``/二次求根无条件接管，会有 ≤1 ulp 差异，破坏该回归锚。

适用范围（本骨架的边界，勿静默外推）
------------------------------------
全场未屈服（处处 ``|τ̃|≤τ̃_Y``）⇒ 间隙无流动、``I₁=0``、``q₀`` 无定义，抛
``ValueError``。定均速模式（A.2.2）的收敛判据用 ``‖Δλ̃‖_p`` 替代定 G 模式的
``‖λ̃‖_p``（后者在定均速下 → ``−G̃ỹ ≠ 0``，照搬会永不收敛；见该函数 docstring）。

References
----------
Bararpour & Frigaard (2025), JFM 1022 A15：附录 A.2.1（(A16)—(A21)）、A.2.2
（(A22)—(A30)）、A.2.3（交错网格、ρ=r=1、p 范数判据）；闭包 (2.13)—(2.15)。q₀
按通量比定义（Z&F22 (4.25)/(4.28) 等密度闭式在 ``Gb=0`` 下为回归锚）。
Zhang & Frigaard (2022), JFM 947 A32：(4.21a)/(4.21b)、(4.25)/(4.28)。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]

# Gauss–Legendre 节点（每个光滑段；屈服面/零应力点处已切段 ⇒ 每段解析光滑）
_GL_ORDER = 16
_GL_X, _GL_W = np.polynomial.legendre.leggauss(_GL_ORDER)
# (A20)/(A29) 内层二分的固定步数（区间宽度至少缩 2⁻⁶⁰）
_BISECT_STEPS = 60


@dataclass(frozen=True)
class GapSolution:
    """一维间隙解与闭包量（口径见模块 docstring）。

    Attributes:
        I1: 平均流动度（**量纲**，``H²·Ĩ₁``；= Z&F22 (4.21a) / ``mobility_i1`` 口径）。
        I2: 浮力流动度（**量纲**，``H³·Ĩ₂``；= Z&F22 (4.21b) / ``mobility_i2`` 口径）。
        q0: 各向同性通量函数（无量纲，**通量比**口径；``q0(0)=0``、``q0(1)=1``）。
        u: tilde 系速度剖面 ；共线/轴向输入下是长 ``ny`` 的 1 维数组，非共线输入下
            是 ``(2, ny)`` 数组（分量 = 输入坐标系）。``u(1)=0`` 恒成立。
        y: 归一化坐标 ỹ∈[0,1]（``linspace(0,1,ny)``）。
        G: 定 G 模式 = 输入的修改压力梯度；定均速模式 = **解出**的 G
            （归一化为 2-元组，可直接传回 ``solve_fixed_G``）。
        iters: Uzawa 实际迭代次数（含收敛判定的那一步）。
        converged: 恒为 ``True``——不收敛时抛 ``RuntimeError``，不会返回 ``False``。
    """

    I1: float
    I2: float
    q0: float
    u: Array
    y: Array
    G: tuple[float, float]
    iters: int
    converged: bool

    @property
    def u_bar(self):
        """间隙平均流速 ū = ``np.trapezoid(u, y)``（**全项目唯一口径**）。

        ⚠️ 不要用 ``u.mean()``：格边网格含 ỹ=0 与 =1 两端点，``mean`` 是 O(h) 偏差
        （``ny=201`` 牛顿槽流实测 1.3e-3），``trapezoid`` 为 O(h²)（同例 2.1e-5）。
        非共线输入下 ``u`` 是 ``(2, ny)`` ⇒ 返回 ``(2,)`` 的逐分量均值。
        """
        return np.trapezoid(self.u, self.y)


def solve_fixed_G(
    c_bar: float,
    n: tuple[float, float],
    kappa: tuple[float, float],
    tau_y: tuple[float, float],
    G,
    Gb=(0.0, 0.0),
    *,
    H: float = 1.0,
    ny: int = 201,
    r: float = 1.0,
    tol: float = 1e-10,
    max_iter: int = 2000,
) -> GapSolution:
    """B&F25 附录 A.2.1：固定 G/Gb 求间隙弱解 ũ，再算 Ĩ₁/Ĩ₂/q₀。

    Args:
        c_bar: 间隙平均顶替液（水泥）体积分数 c̄ = ỹ_i ∈ [0,1]（界面位置）。
        n: 幂律指数 ``(n₁, n₂)``——流体 1（壁面带，被顶替液）/ 流体 2（中线带）。
        kappa: 稠度系数 ``(κ₁, κ₂)``，与 n 同下标对应。
        tau_y: 屈服应力 ``(τ_{Y,1}, τ_{Y,2})``，与 n 同下标对应。
        G: 修改压力梯度（2 维向量或标量；标量视作第 1 分量非零）。
        Gb: 浮力向量（同上）；须与 ``G`` 共线（见模块 docstring 的适用域）。
        H: 物理半隙（(A4) 的 H，量纲化因子与 Z&F22 (4.21a,b) 同口径）。
        ny: 输出/格边网格点数（格心数 = ``ny−1``）。
        r: Uzawa 数值参数（(A16)/(A18)/(A19) 的 r；固定 1 满足 (A31)），
            非几何量（环空半径为 ``r_a``，与 ``H`` 无关）。
        tol: 收敛容差（三个 ``p`` 范数同时低于此值）。
        max_iter: Uzawa 最大迭代步数。

    Returns:
        :class:`GapSolution`。

    Raises:
        ValueError: 输入越界（c̄、n、κ、τ_Y、H、ny、r、max_iter）；``G``/``Gb`` 非共线或
            同时为零；全场未屈服（``I₁=0``，闭包无定义）。
        RuntimeError: ``max_iter`` 步内未达 ``tol``（不静默返回未收敛解）。
    """
    c_bar = float(c_bar)
    H = float(H)
    r = float(r)
    tol = float(tol)
    n = (float(n[0]), float(n[1]))
    kappa = (float(kappa[0]), float(kappa[1]))
    tau_y = (float(tau_y[0]), float(tau_y[1]))
    _validate_inputs(c_bar, n, kappa, tau_y, H, ny, r, max_iter)

    G_vec, Gb_vec = _as_2vec(G, "G"), _as_2vec(Gb, "Gb")
    # 方向归约：共线 ⇒ 旋转到共同方向（精确 + 逐位回归）；非共线 ⇒ 原样向量
    G_t, Gb_t = _reduce_directions(G_vec, Gb_vec)

    # (A4) 换到 tilde 系（H=1）
    kappa_t = (kappa[0] / H**n[0], kappa[1] / H**n[1])
    G_t, Gb_t = H * G_t, H * Gb_t

    # 闭包量：由解析应力场 quadrature（(2.14)/(2.15) + q₀ 通量比）
    I1_t, I2_t, q0 = _closure_integrals(c_bar, n, kappa_t, tau_y, G_t, Gb_t)

    # 定 G 弱解（(A16)—(A21)）
    y, u, iters = _uzawa_fixed_G(c_bar, n, kappa_t, tau_y, G_t, Gb_t, ny, r, tol, max_iter)
    if _is_axial(G_t, Gb_t):
        u = u[0]  # 轴向/共线 ⇒ 1 维剖面（A-1a 起的形状契约，见模块 docstring）

    return GapSolution(
        I1=H**2 * I1_t,
        I2=H**3 * I2_t,
        q0=q0,
        u=u,
        y=y,
        G=(float(G_vec[0]), float(G_vec[1])),
        iters=iters,
        converged=True,
    )


def solve_fixed_mean_velocity(
    c_bar,
    n,
    kappa,
    tau_y,
    u_bar,
    Gb=(0.0, 0.0),
    *,
    H: float = 1.0,
    ny: int = 201,
    r: float = 1.0,
    tol: float = 1e-10,
    max_iter: int = 2000,
) -> GapSolution:
    """B&F25 附录 A.2.2：给定平均流速 u_bar（与 Gb）反求 G 并解出剖面。

    D2DGA 外迭代（Task 5）实际需要的入口：外迭代给出 ū（体积守恒）与 Gb，需要
    反求对应的修改压力梯度 G 与剖面。形参名/顺序与 :func:`solve_fixed_G` 一致，
    仅把输入 ``G`` 换成 ``u_bar``。

    算法（(A22)—(A30)；符号约定沿用 A-1a 起的镜像约定，勿改）
    -------------------------------------------------------
    1. (A22)：``u_P,0 = (1-y^2)/(2r)``（单位压力梯度槽流的解析解，ū_P,0 = 1/(3r)）；
    2. (A23)：``lam0`` **只平衡浮力**（与定 G 模式的 (A17) 不同：那里还要平衡 G）。
       本实现的镜像式：``lam0 = -(1-cbar)*Gb*y`` / ``-cbar*Gb*(1-y)``；
    3. (A26)：``du_I/dy = q^k - lam_t^k/r``（积分常数沿用定 G 模式的 C=0），由壁面
       ``u_I(1)=0`` 回积；
    4. (A27)/(A28)：``alpha = (u_bar* - ubar_I)/ubar_P0``（**逐分量**），
       ``u = u_I + alpha*u_P,0``，``G^{k+1} = alpha`` —— G 由**平均流速约束**定出，
       不靠外层迭代；
    5. (A29)：``q^{k+1}`` 与 (A20) 同式，``m = lam + lam_t + r*(du^{k+1}/dy)``，注意用
       **更新后**的总速度梯度（定 G 模式里用 q^k 的等价形式在这里不成立：定均速下
       ``lam_t`` 不趋于零）；
    6. (A30)：``lam_t^{k+1} = lam_t^k + rho(du^{k+1}/dy - q^{k+1})``，rho=1。

    自洽性：``solve_fixed_mean_velocity(u_bar=u*)`` 解出的 ``G`` 代回
    :func:`solve_fixed_G` 必复现 ``u*``（``u_bar`` 口径）——两模式的不动点逐式相同。

    ⚠️ 收敛判据的**唯一偏离**：文献 (A.2.3) 监控 ``||du||_p, ||dq||_p, ||lam_t||_p``。
    定 G 模式下 ``lam_t -> 0`` ⇒ 第三项可用；定均速模式下 ``lam_t -> -G*y != 0``
    （lam0 只含浮力，G 的有效部分全在 lam_t 里）⇒ 照搬 ``||lam_t||_p < tol``
    **永不成立**，本实现改用 ``||d(lam_t)||_p``；前两项与 ``p = 1+min(n1,n2)`` 完全沿用。

    Args:
        u_bar: 目标平均流速 ū*（2 维向量或标量，标量视作第 1 分量）。约束在格边网格的
            ``np.trapezoid`` 口径上精确满足（与 :attr:`GapSolution.u_bar` 同一口径）；
            相应地把 (A22) 的 ū_P,0 取其同网格 trapezoid 值（解析 1/(3r) 的 O(h^2) 版本）。
        Gb: 浮力向量（同上）。
        其余同 :func:`solve_fixed_G`；``G`` 不再是输入而是解出量。

    Returns:
        :class:`GapSolution`：``G`` = **解出**的 G（量纲，= tilde 解 / H）；``u`` = 最终
        剖面（求解过程中就地得到，非事后反算）；``I1/I2/q0`` 在解出的 G 处按定 G 模式的
        同一闭包求积——``Gb=0`` 时恒有 ``ū = I1·G``（即"平均流动度"），供 B&F25 图 6 用。

    Raises:
        ValueError: 输入越界；``u_bar`` 形参非法；``u_bar`` 与 ``Gb`` 同时为零；
            全场未屈服（``I1=0``，闭包无定义）。
        RuntimeError: ``max_iter`` 步内未达 ``tol``，或数值发散（不静默回退）。
    """
    c_bar = float(c_bar)
    H = float(H)
    r = float(r)
    tol = float(tol)
    n = (float(n[0]), float(n[1]))
    kappa = (float(kappa[0]), float(kappa[1]))
    tau_y = (float(tau_y[0]), float(tau_y[1]))
    _validate_inputs(c_bar, n, kappa, tau_y, H, ny, r, max_iter)

    ubar_vec = _as_2vec(u_bar, "u_bar")
    Gb_vec = _as_2vec(Gb, "Gb")
    if float(np.linalg.norm(ubar_vec)) <= 0.0 and float(np.linalg.norm(Gb_vec)) <= 0.0:
        raise ValueError(
            "u_bar 与 Gb 不能同时为零：无驱动、无目标流量 ⇒ G 无定义（闭包无意义）"
        )

    # (A4) 换到 tilde 系（H=1）：速度与应力同值，唯 κ/压力梯度按 H 缩放
    kappa_t = (kappa[0] / H**n[0], kappa[1] / H**n[1])
    Gb_t = H * Gb_vec

    y, u, G_t, iters = _uzawa_fixed_mean_velocity(
        c_bar, n, kappa_t, tau_y, ubar_vec, Gb_t, ny, r, tol, max_iter
    )
    I1_t, I2_t, q0 = _closure_integrals(c_bar, n, kappa_t, tau_y, G_t, Gb_t)
    if _is_axial(G_t, Gb_t):
        u = u[0]  # 轴向/共线 ⇒ 1 维剖面（形状契约见模块 docstring）

    return GapSolution(
        I1=H**2 * I1_t,
        I2=H**3 * I2_t,
        q0=q0,
        u=u,
        y=y,
        G=(float(G_t[0] / H), float(G_t[1] / H)),
        iters=iters,
        converged=True,
    )


def _uzawa_fixed_mean_velocity(c_bar, n, kappa_t, tau_y, ubar_vec, Gb_t, ny, r, tol, max_iter):
    """(A22)—(A30) 定均速 Uzawa，返回 ``(y, u, G_t, iters)``（G_t 在 tilde 系）。

    与定 G 模式的差别（① 更新次序、② m 的梯度取更新后的总值、③ 收敛判据用 d(lam_t)）
    见 :func:`solve_fixed_mean_velocity` 的 docstring。
    """
    n_cell = ny - 1
    h = 1.0 / n_cell
    y_edge = np.linspace(0.0, 1.0, ny)
    y_cen = (np.arange(n_cell) + 0.5) * h

    in_fluid2 = y_cen < c_bar  # 中线带（顶替液）
    kappa_c = np.where(in_fluid2, kappa_t[1], kappa_t[0])
    n_c = np.where(in_fluid2, n[1], n[0])
    tauy_c = np.where(in_fluid2, tau_y[1], tau_y[0])

    # (A23) + 镜像符号：lam0 只平衡浮力（band2 斜率 -(1-cbar)Gb，band1 斜率 +cbar*Gb）
    lam0 = np.where(
        in_fluid2,
        -(1.0 - c_bar) * Gb_t[:, None] * y_cen,
        -c_bar * Gb_t[:, None] * (1.0 - y_cen),
    )

    p_norm = 1.0 + min(n[0], n[1])  # (A.2.3)：p = 1 + min(n1,n2)
    u_p0 = (1.0 - y_edge**2) / (2.0 * r)  # (A22) 解析解
    ubar_p0 = float(np.trapezoid(u_p0, y_edge))  # = 1/(3r) 的同一网格 trapezoid 版本
    u_p0_du = y_cen / r  # u_P,0 对 du/dy = -u' 的贡献（= +y/r）

    q = np.zeros((2, n_cell))
    lam_t = np.zeros((2, n_cell))  # lam_t^k = lam^k - lam0
    u = np.zeros((2, ny))
    alpha = np.zeros(2)
    with np.errstate(over="ignore", invalid="ignore"):
        for iters in range(1, max_iter + 1):
            # ① (A26)：du_I/dy = q^k - lam_t^k/r，由壁面 u_I(1)=0 回积
            du_I_dy = q - lam_t / r
            u_I = _integrate_from_wall(du_I_dy, h, ny)

            # ② (A27)/(A28)：alpha 逐分量；(A25)：u = u_I + alpha*u_P,0；G^{k+1} = alpha
            alpha = (ubar_vec - np.trapezoid(u_I, y_edge, axis=-1)) / ubar_p0
            u_new = u_I + alpha[:, None] * u_p0[None, :]
            du_dy = du_I_dy + alpha[:, None] * u_p0_du[None, :]  # 总梯度（= -u'）

            # ③ (A29)：q^{k+1}（局部极小化，|m|<=tau_Y ⇒ 0），m 用**更新后**的总梯度
            q_new = _local_q_min(lam0 + lam_t + r * du_dy, kappa_c, n_c, tauy_c, r)

            # ④ (A30)：lam_t^{k+1} = lam_t^k + rho(du^{k+1}/dy - q^{k+1})，rho=1
            lam_new = lam_t + (du_dy - q_new)

            res_u = _pnorm(_vec_abs(u_new - u), p_norm)
            res_q = _pnorm(_vec_abs(q_new - q), p_norm)
            res_lam = _pnorm(_vec_abs(lam_new - lam_t), p_norm)  # 定均速：用 d(lam_t)

            if not (np.isfinite(res_u) and np.isfinite(res_q) and np.isfinite(res_lam)):
                raise RuntimeError(
                    f"Uzawa（定均速）数值发散（溢出/NaN，第 {iters} 步）："
                    f"||du||_p={res_u:.3e}, ||dq||_p={res_q:.3e}, ||d(lam_t)||_p={res_lam:.3e}。"
                    f"r={r!r} 过小（lam_t 模态收缩因子 ~|1-1/r|，r<~0.5 不收敛）；"
                    "补救方向：取 r~1（r>1 反而更慢）并增大 max_iter。"
                )
            u, q, lam_t = u_new, q_new, lam_new
            if res_u < tol and res_q < tol and res_lam < tol:
                return y_edge, u, alpha, iters

    raise RuntimeError(
        f"Uzawa（定均速）在 max_iter={max_iter} 步内未收敛："
        f"||du||_p={res_u:.3e}, ||dq||_p={res_q:.3e}, ||d(lam_t)||_p={res_lam:.3e}, "
        f"tol={tol:.1e}（三分量须同时低于 tol；补救方向是取 r~1 并增大 max_iter）"
    )


# --------------------------------------------------------------------------- #
# 输入校验与向量归约
# --------------------------------------------------------------------------- #


def _validate_inputs(c_bar, n, kappa, tau_y, H, ny, r, max_iter) -> None:
    if not 0.0 <= c_bar <= 1.0:
        raise ValueError(f"c_bar 须在 [0,1] 内，得到 c_bar={c_bar!r}")
    if not H > 0.0:
        raise ValueError(f"H 须为正（物理半隙），得到 H={H!r}")
    if not r > 0.0:
        raise ValueError(f"r 须为正（Uzawa 数值参数），得到 r={r!r}")
    if not isinstance(ny, int) or isinstance(ny, bool) or ny < 3:
        raise ValueError(f"ny 须为 ≥3 的整数（格心数 = ny−1），得到 ny={ny!r}")
    if not isinstance(max_iter, int) or isinstance(max_iter, bool) or max_iter < 1:
        raise ValueError(f"max_iter 须为 ≥1 的整数，得到 max_iter={max_iter!r}")
    if any(not x > 0.0 for x in n):
        raise ValueError(f"幂律指数须为正（1/n 需有定义），得到 n={n!r}")
    if any(not k > 0.0 for k in kappa):
        raise ValueError(f"稠度系数须为正，得到 kappa={kappa!r}")
    if any(t < 0.0 for t in tau_y):
        raise ValueError(f"屈服应力须非负，得到 tau_y={tau_y!r}")


def _as_2vec(v, name: str) -> Array:
    """标量或 2 维向量 → 长度为 2 的 float 数组（标量视作第 1 分量）。

    ⚠️ 沿用本项目 ``two_layer``/``stream_function`` 的二维向量口径：标量输入
    按"沿第 1 分量"（φ 方向）解释——调用方若用标量，应确保这与几何方向一致。
    """
    arr = np.asarray(v, dtype=float).ravel()
    if arr.size == 1:
        return np.array([float(arr[0]), 0.0])
    if arr.size == 2:
        return arr
    raise ValueError(f"{name} 须为标量或 2 维向量，得到 shape={np.shape(v)}")


def _reduce_directions(G_vec: Array, Gb_vec: Array) -> tuple[Array, Array]:
    """把 ``(G, Gb)`` 归约到内部求解坐标系（返回两个长度 2 的数组）。

    共线（含 ``Gb=0``）时旋转到共同方向 ``e = Ĝ``（``G=0`` 时取 ``Ĝb``），返回
    ``(Gs, 0)`` 与 ``(Gbs, 0)``：``Gs=|G|``、``Gbs=Gb·e``（有号）。共线旋转对各向同性
    问题是**精确等价**的，且归约后第 2 分量恒 0 ⇒ 下游全部走标量分支，与向量化改造前
    **逐位一致**（R-T1-3 的回归锚，基线见 ``.tmp_research/task3_baseline/``）。
    非共线时原样返回 ``(G, Gb)``（无共同方向，下游走向量分支）。
    """
    G_mag = float(np.linalg.norm(G_vec))
    Gb_mag = float(np.linalg.norm(Gb_vec))
    if G_mag <= 0.0:
        if Gb_mag <= 0.0:
            raise ValueError(
                "G 与 Gb 不能同时为零：无缝驱动应力 ⇒ 流动度为 0，闭包无定义"
            )
        e = Gb_vec / Gb_mag
    else:
        e = G_vec / G_mag
    Gbs = float(np.dot(Gb_vec, e))
    perp = Gb_vec - Gbs * e
    if float(np.linalg.norm(perp)) > 1e-9 * max(1.0, Gb_mag):
        return G_vec, Gb_vec  # 非共线 ⇒ 向量分支（不旋转）
    return np.array([G_mag, 0.0]), np.array([Gbs, 0.0])


# --------------------------------------------------------------------------- #
# 定 G 模式：Uzawa（A.2.1）
# --------------------------------------------------------------------------- #


def _stress_vectors(c_bar, G_t: Array, Gb_t: Array):
    """逐层应力**向量**系数 ``(c2, d2, c1, d1)``（(A1)/(A2)/(A17)，**符号见下**）。

    ``τ̃₂(ỹ) = c₂· + d₂``（中线带，ỹ∈[0,c̄]，τ̃₂(0)=0）；
    ``τ̃₁(ỹ) = c₁·ỹ``（壁面带，ỹ∈[c̄,1]，界面剪应力连续；``d₁=(c₂−c₁)c̄``，非共线时 ``d₁`` 不平行于 ``c₁``）。

    ⚠️ 符号约定：``G`` 是驱动方向的修改压力梯度（(2.13) 口径）：
    ``Ī₁>0`` ⇒ 间隙平均流速与 G 同向。故本实现的 λ₀ 取附录 A (A17) 的**相反号**
    （即 ``λ₀ = −λ₀^{(A17)}``）：附录 A 的 (A1)/(A16) 印刷式与 (2.13) 反号，
    照抄会使 ``G>0`` 时 ``ũ<0``（与 (2.13) 及 Z&F22 口径矛盾）。闭包量
    ``Ĩ₁/Ĩ₂/q₀`` 只依赖 ``|τ̃|``，此选择不影响它们。
    """
    c2 = G_t - (1.0 - c_bar) * Gb_t
    c1 = G_t + c_bar * Gb_t
    return c2, np.zeros(2), c1, (c2 - c1) * c_bar


def _is_axial(c, d) -> bool:
    """应力仿射场是否轴对齐（第 2 分量恒 0 ⇒ 走标量分支 = 逐位回归锚）。"""
    return c[1] == 0.0 and d[1] == 0.0


def _tau_abs(c, d, y):
    """``|tau(ỹ)| = |c·ỹ + d|``（2 维仿射场的模；HB 本构各向同性 ⇒ 只依赖它）。

    ⚠️ 轴对齐时返回 ``|c₀·ỹ+d₀|``——与向量化前的标量实现逐位一致；若让 ``hypot``
    无条件接管会引入 ≤1 ulp 差异，破坏共线回归锚。
    """
    if _is_axial(c, d):
        return np.abs(c[0] * y + d[0])
    return np.hypot(c[0] * y + d[0], c[1] * y + d[1])


def _vec_abs(v):
    """分量数组逐格模 ``|v|``（轴对齐时 = ``|v₀|``，逐位回归锚）。"""
    if not np.any(v[1] != 0.0):
        return np.abs(v[0])
    return np.hypot(v[0], v[1])


def _local_q_min(m, kappa_c, n_c, tauy_c, r):
    """(A20)/(A29) 局部极小化：``q = θ·m``（``|m|≤τ_Y`` 取 0），逐格向量模。

    ``m`` 为 ``(2, n_cell)``，返回同形状 ``q``——两组分**唯一**的耦合处就是 ``|m|``。
    """
    q_new = np.zeros_like(m)
    abs_m = _vec_abs(m)
    yielded = abs_m > tauy_c
    if np.any(yielded):
        theta = _bisect_theta(
            abs_m[yielded], kappa_c[yielded], n_c[yielded], tauy_c[yielded], r
        )
        q_new[:, yielded] = theta * m[:, yielded]
    return q_new


def _integrate_from_wall(du_dy, h, ny):
    """由壁面 ``u(1)=0`` 向回累积：``u(y_i) = h·Σ_{j≥i} du_dy_j``（(A18)/(A26)）。

    输入 ``(2, ny-1)`` → 输出 ``(2, ny)``（第 1 维 = 分量）；中点求积，对线性 γ̇ 精确。
    """
    return np.concatenate(
        [h * np.cumsum(du_dy[:, ::-1], axis=1)[:, ::-1], np.zeros((2, 1))], axis=1
    )


def _uzawa_fixed_G(c_bar, n, kappa_t, tau_y, G_t, Gb_t, ny, r, tol, max_iter):
    """(A16)—(A21) 的 Uzawa 迭代（ρ=r=1），返回 ``(y, u, iters)``。

    收敛快慢由 ``r`` 决定：``q`` 模态收缩因子 ``r/(κ̃+r)``（r≈1 最快，r>1 更慢），
    ``λ̃`` 模态因子 ~``|1−1/r|``（r≲0.5 不收敛）⇒ 不收敛时的补救方向是
    **r≈1 + 增大 max_iter**（不是增大 r）。
    """
    n_cell = ny - 1
    h = 1.0 / n_cell
    y_edge = np.linspace(0.0, 1.0, ny)
    y_cen = (np.arange(n_cell) + 0.5) * h

    in_fluid2 = y_cen < c_bar  # 中线带（顶替液）
    kappa_c = np.where(in_fluid2, kappa_t[1], kappa_t[0])
    n_c = np.where(in_fluid2, n[1], n[0])
    tauy_c = np.where(in_fluid2, tau_y[1], tau_y[0])

    c2, d2, c1, d1 = _stress_vectors(c_bar, G_t, Gb_t)
    lam0 = np.where(in_fluid2, c2[:, None] * y_cen, c1[:, None] * y_cen + d1[:, None])

    p_norm = 1.0 + min(n[0], n[1])  # (A.2.3)：p = 1 + min(n₁,n₂)
    q = np.zeros((2, n_cell))
    lam_t = np.zeros((2, n_cell))  # λ̃^k（λ^k = λ₀ + λ̃^k）
    u = np.zeros((2, ny))
    with np.errstate(over="ignore", invalid="ignore"):
        for iters in range(1, max_iter + 1):
            # ① (A19)/(A20)：q^{k+1}（局部极小化；|m|≤τ_Y ⇒ q=0）
            m = lam0 + r * q
            q_new = _local_q_min(m, kappa_c, n_c, tauy_c, r)

            # ② (A18)：dũ^{k+1}/dỹ|_格心 = q^k − λ̃^k/r（对称面在 ỹ=0 ⇒ C=0），
            #    由壁面 (1)=0 向回累积（中点求积；与格心处的向后差分逐位一致）
            du_dy = q - lam_t / r
            u_new = _integrate_from_wall(du_dy, h, ny)

            # ③ (A21)：λ̃^{k+1} = λ̃^k + ρ(dũ^{k+1}/dỹ − q^{k+1})，ρ=1
            lam_new = lam_t + (du_dy - q_new)

            res_u = _pnorm(_vec_abs(u_new - u), p_norm)
            res_q = _pnorm(_vec_abs(q_new - q), p_norm)
            res_lam = _pnorm(_vec_abs(lam_new), p_norm)

            if not (np.isfinite(res_u) and np.isfinite(res_q) and np.isfinite(res_lam)):
                raise RuntimeError(
                    f"Uzawa（定 G）数值发散（溢出/NaN，第 {iters} 步）："
                    f"‖Δũ‖_p={res_u:.3e}, ‖Δq‖_p={res_q:.3e}, ‖λ̃‖_p={res_lam:.3e}。"
                    f"r={r!r} 过小（λ̃ 模态收缩因子 ~|1−1/r|，r≲0.5 不收敛）；"
                    "补救方向：取 r≈1（r>1 反而更慢）并增大 max_iter。"
                )
            u, q, lam_t = u_new, q_new, lam_new
            if res_u < tol and res_q < tol and res_lam < tol:
                return y_edge, u, iters

    raise RuntimeError(
        f"Uzawa（定 G）在 max_iter={max_iter} 步内未收敛："
        f"‖Δũ‖_p={res_u:.3e}, ‖Δq‖_p={res_q:.3e}, ‖λ̃‖_p={res_lam:.3e}, tol={tol:.1e}"
        f"（三分量须同时低于 tol。收敛快慢由 r 决定：q 模态收缩因子 r/(κ̃+r) ⇒ r≈1 最快，"
        "r>1 反而更慢，r≲0.5 不收敛 ⇒ 补救方向是取 r≈1 并增大 max_iter，不是增大 r）"
    )


def _bisect_theta(abs_m, kappa, n, tau_y, r):
    """(A20)：解 ``κ̃θⁿ|m|ⁿ + rθ|m| = |m| − τ_Y`` 的 θ>0 根（逐格向量化二分）。

    区间 ``[0, (|m|−τ_Y)/(r|m|)]``：左端 LHS=0<RHS，右端由 ``κ̃θⁿ|m|ⁿ ≥ 0``
    保证 LHS≥RHS；固定 ``_BISECT_STEPS`` 步（区间宽度缩 2⁻⁶⁰）。调用方保证
    ``|m| > τ_Y``（否则 q 直接取 0）。
    """
    lo = np.zeros_like(abs_m)
    hi = (abs_m - tau_y) / (r * abs_m)
    for _ in range(_BISECT_STEPS):
        mid = 0.5 * (lo + hi)
        f = kappa * mid**n * abs_m**n + r * mid * abs_m - (abs_m - tau_y)
        below = f < 0.0
        lo = np.where(below, mid, lo)
        hi = np.where(below, hi, mid)
    return 0.5 * (lo + hi)


def _pnorm(v, p: float) -> float:
    """离散 p 范数（格点平均口径）：``(mean|v_i|^p)^{1/p}``（A.2.3 的 ‖·‖_p）。"""
    arr = np.abs(np.asarray(v, dtype=float))
    if arr.size == 0:
        return 0.0
    return float(np.mean(arr**p) ** (1.0 / p))


def _weight_y2(yv):
    """ỹ² 权重（(2.14)/(2.22) 第一项的被积权重）。"""
    return yv * yv


# --------------------------------------------------------------------------- #
# 闭包量：(2.14)/(2.15)/(2.22) 在解析应力场上求积
# --------------------------------------------------------------------------- #


def _closure_integrals(c_bar, n, kappa_t, tau_y, G_t, Gb_t):
    """tilde 系下的 ``(Ĩ₁, Ĩ₂, q₀)``（(2.14)/(2.15)/(2.22)）。

    逐层积分被 ỹ=c̄ 与屈服面切开；被积函数 ``weight(ỹ)·(1/η̃)``，其中
    ``1/η̃ = γ̇(τ̃)/|τ̃|``（(A3) 有效黏度），未屈服处为 0。
    """
    c2, d2, c1, d1 = _stress_vectors(c_bar, G_t, Gb_t)
    # 中线带（流体 2）：τ̃ = a₂·ỹ，ỹ∈[0,c̄]
    int2_i1 = _gap_integral(0.0, c_bar, c2, d2, kappa_t[1], n[1], tau_y[1], _weight_y2)
    # 壁面带（流体 1）：τ̃ = (a₂−a₁)c̄ + a₁·ỹ，ỹ∈[c̄,1]
    int1_i1 = _gap_integral(c_bar, 1.0, c1, d1, kappa_t[0], n[0], tau_y[0], _weight_y2)
    int1_i2 = _gap_integral(
        c_bar, 1.0, c1, d1, kappa_t[0], n[0], tau_y[0], lambda yv: yv * (1.0 - yv)
    )

    I1_t = int2_i1 + int1_i1
    if not I1_t > 0.0:
        raise ValueError(
            "全场未屈服（|τ̃(ỹ)| ≤ τ̃_Y 处处成立）：间隙无流动 ⇒ Ĩ₁=0、q₀ 无定义"
            "（static plug 情形，须由调用方另行处理）"
        )
    I2_t = (1.0 - c_bar) * int2_i1 + c_bar * int1_i2
    q0 = _flux_ratio_q0(c_bar, c2, d2, c1, d1, kappa_t, n, tau_y)
    return I1_t, I2_t, q0


def _gap_integral(lo, hi, c, d, kappa, n, tau_y, weight) -> float:
    """``∫_lo^hi weight()·(1/η̃)(ỹ) dỹ``，应力仿射 ``τ̃ = p_ + q_·ỹ``。

    按屈服面（``|τ̃|=τ̃_Y``）与零应力点（``τ̃=0``）切段，逐段 Gauss–Legendre；
    每段被积函数解析光滑 ⇒ 界面/屈服面无离散误差。
    """
    total = 0.0
    for a, b in _yield_pieces(lo, hi, c, d, tau_y):
        yq = 0.5 * (b - a) * _GL_X + 0.5 * (a + b)
        inv_eta = _inv_eff_viscosity(_tau_abs(c, d, yq), kappa, n, tau_y)
        total += 0.5 * (b - a) * float(np.sum(_GL_W * weight(yq) * inv_eta))
    return total


def _yield_pieces(lo, hi, c, d, tau_y) -> list[tuple[float, float]]:
    """把 ``[lo,hi]`` 在屈服面/零应力点处切开（返回可能为空的段列表）。

    轴对齐时与向量化前的标量实现逐式相同（切点 ``(+-tau_Y-d0)/c0``、``-d0/c0``）；
    非轴向时解二次方程 ``|c*y+d|^2 = tau_Y^2`` 与 ``|c*y+d|^2 = 0``（零应力点，
    被积函数在此不可导）。``c=0`` 而 ``d!=0`` 的层退化为一元线性/无根。
    """
    if not hi > lo:
        return []
    if _is_axial(c, d):
        p_, q_ = d[0], c[0]
        cuts = [] if q_ == 0.0 else [
            (tau_y - p_) / q_, (-tau_y - p_) / q_, -p_ / q_
        ]
    else:
        A = float(c @ c)
        B = 2.0 * float(c @ d)
        cuts = _quad_roots(A, B, float(d @ d) - tau_y * tau_y)
        cuts += _quad_roots(A, B, float(d @ d))
    edges = [lo] + sorted(x for x in cuts if lo < x < hi) + [hi]
    return [(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]


def _quad_roots(A, B, C) -> list[float]:
    """``A*y^2 + B*y + C = 0`` 的实根（``A<=0`` 退化为线性；无实根返回空表）。"""
    if A <= 0.0:
        return [] if B == 0.0 else [-C / B]
    disc = B * B - 4.0 * A * C
    if disc < 0.0:
        return []
    r = float(np.sqrt(disc))
    return [(-B - r) / (2.0 * A), (-B + r) / (2.0 * A)]


def _gap_integral_vec(lo, hi, c, d, kappa, n, tau_y, weight) -> Array:
    """``int weight(y)*gamma_vec(y) dy``（2 维向量值，返回长度 2 数组）。

    ``gamma_vec = (1/eta_tilde)*tau_tilde``（(A3)；HB 本构各向同性 ⇒ ``1/eta_tilde``
    只依赖 ``|tau_tilde|``），未屈服处为 0。段落切法与 ``_gap_integral`` 相同。
    """
    total = np.zeros(2)
    for a, b in _yield_pieces(lo, hi, c, d, tau_y):
        yq = 0.5 * (b - a) * _GL_X + 0.5 * (a + b)
        inv_eta = _inv_eff_viscosity(_tau_abs(c, d, yq), kappa, n, tau_y)
        total = total + 0.5 * (b - a) * np.sum(
            _GL_W * weight(yq) * (inv_eta * (c[:, None] * yq + d[:, None])), axis=-1
        )
    return total


def _flux_ratio_q0(c_bar, c2, d2, c1, d1, kappa_t, n, tau_y) -> float:
    """``q0 = |int_0^cbar u dy| / |int_0^1 u dy|``（通量比 = 顶替液通量份额）。

    被积函数 ``u(y) = int_y^1 gamma dy'`` 的二重积分换序化为单重::

        int_0^1 u dy = int_0^1 y*gamma dy
        int_0^cbar u dy = int_0^cbar y*gamma dy + cbar*int_cbar^1 gamma dy

    两者都在**解析应力场**上按屈服面切段 + 逐段 Gauss-Legendre（各段光滑）⇒ 与
    ``ny`` 无关、端点精确：``q0(0)=0``（分子是空积分）、``q0(1)=1``（分子分母是同一
    浮点和）。``Gb=0`` 时与 Z&F22 等密度闭式 ``isotropic_flux_q0`` 一致（回归锚）；
    非共线时按题注口径取**逐分量求积后的模之比**。
    """
    w_y = lambda yv: yv  # noqa: E731
    w_one = lambda yv: np.ones_like(yv)  # noqa: E731

    den2 = _gap_integral_vec(0.0, c_bar, c2, d2, kappa_t[1], n[1], tau_y[1], w_y)
    den1 = _gap_integral_vec(c_bar, 1.0, c1, d1, kappa_t[0], n[0], tau_y[0], w_y)
    int1_one = _gap_integral_vec(c_bar, 1.0, c1, d1, kappa_t[0], n[0], tau_y[0], w_one)

    num_abs = float(_vec_abs(den2 + c_bar * int1_one))
    den_abs = float(_vec_abs(den2 + den1))
    return num_abs / den_abs


def _inv_eff_viscosity(tau_abs, kappa, n, tau_y) -> Array:
    """有效黏度倒数 ``1/η̃ = γ̇/|τ̃|``（(A3)）；未屈服（``|τ̃|≤τ̃_Y``）取 0。"""
    tau = np.asarray(tau_abs, dtype=float)
    out = np.zeros_like(tau)
    yielded = tau > tau_y
    if np.any(yielded):
        gamma_dot = ((tau[yielded] - tau_y) / kappa) ** (1.0 / n)
        out[yielded] = gamma_dot / tau[yielded]
    return out