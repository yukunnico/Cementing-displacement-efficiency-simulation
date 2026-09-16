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
  tilde 系速度 ũ(ỹ)，即量纲速度在归一化位置处的取值。
- **``u`` 的形状契约**：轴向/共线输入（``G`` 与 ``Gb`` 共线，含全零 ``Gb``）下 ``u``
  是长 ``ny`` 的 1 维数组（这是 A.2.1/A.2.2 的标量情形，与 A-1a 起的行为逐位一致）；
  **非共线**输入下问题是 2 维向量场，``u`` 是 ``(2, ny)`` 数组（第 0/1 分量 = 输入
  坐标系的 φ/ξ 分量）。
  ⚠️ **2 维时 ``u[-1]`` 不是壁面值（它是 ξ 分量的整条剖面）、``u.mean()`` 也不是 ū**（返回
  标量而不是逐分量）。取壁面/末点一律用 **``u[..., -1]``**（两种形状下语义相同；
  一维时还是标量，2 维时是 ``(2,)``）；逐分量通道一律用 ``u[..., i]``，避免 ``u[i]``。
- **ū 的口径是 :attr:`GapSolution.u_bar`（``np.trapezoid(u, y)``），不是 ``u.mean()``**：
  格边网格含 ỹ=0 与 ỹ=1 两端点，``mean`` 带 O(h) 偏差（``ny=201`` 时实测 ~1.3e-3），
  ``trapezoid`` 为 O(h²)（~2.1e-5）。全项目只认 ``u_bar``。

q₀ 口径（通量比，**取代** B&F25 印刷式 (2.22)）
---------------------------------------------
``q₀ = |∫₀^{c̄} ũ dỹ| / |∫^1 ũ dỹ|``（逐分量求积后取模之比）——通量比就是 q₀ 的定义。
实现与精度见 :func:`_flux_ratio_q0`；``Gb=0`` 时退回 Z&F22 等密度闭式
``isotropic_flux_q0``（实测 rel≤2.4e-16），是本改动的回归锚。

⚠️ B&F25 **印刷式 (2.22)** 的分子用 ``y_i·∫ỹ/η₁`` 一类与 (2.14) 不同族的权重，
**等密度（``Gb=0``）归一化**下尚可用，但在 ``τ̃ ≠ G̃·ỹ``（即 ``Gb≠0``）时结构上不足
（``q₀`` 实为**各层厚度**的函数）⇒ 不采用印刷式。字面代入印刷式（分子取
``y_i∫ỹ/η₁``）在 ``m=3, c̄=0.5`` 牛顿下确实得 **0.9**；而本实现/解析槽流/Z&F22 (4.25)
三者均为 **0.75**——两个数字都对，差别在定义而非算错。

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

- 中线带（流体 2，ỹ∈[0,c̄]）：``c=G̃−(1−c̄)Gb̃``、``d=0``（τ̃(0)=0 对称面）；
- 壁面带（流体 1，ỹ∈[c̄,1]）：``c=G̃+c̄Gb̃``、``d=(c₂−c₁)c̄``（界面剪应力连续）。
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

批次轴（Task 4.5）
------------------
``solve_fixed_G_batch`` / ``closure_integrals_batch`` 把上述标量内核整体 numpy 化到
``(n_points, n_cells)``：一次求解 ``n_points`` 个独立间隙问题（各自 ``c̄_i`` ⇒ 各自的
界面掩码与逐格 ``κ̃/n/τ̃_Y`` 数组）。标量路径**完全未改**（其 25 条契约测试逐条不变）。

两条路径的一致性契约是 **rel ≤ 1e-12**：轴向时逐位相同（``‖·‖``/``dot`` 的逐元素式
与标量路径的 ``np.linalg.norm``/``np.dot`` 在含零分量时逐位相等，实测 20000/20000），
非轴向点因 BLAS 末位差可差 ~1 ulp。``u`` 在全体轴向时 ``(n_points, ny)``、存在非共线
点时为 ``(n_points, 2, ny)``。退化格点（全场未屈服 / ``G=Gb=0``）在批量里进
``undefined`` 掩码而**不抛错**（逐点版仍抛 ``ValueError``）——窄边 static wall layer
可能成片出现，须由调用方（如 ``HBClosure`` 的 R-T1-6 地板）处置；不收敛/发散仍抛
``RuntimeError``（与逐点同纪律）。

⚠️ 实测（呼101 真实场 6647 唯一键、``ny=201``）：**闭包-only 0.34 s**（对逐点 83.7 s
⇒ 246×），**含 Uzawa 的全量批量 49.4 s**（仅 1.7×）——瓶颈是 ``_BISECT_STEPS=60``
的逐格二分在批次上的内存带宽（60 步 × ~8 次全数组读写 × Uzawa 迭代数）。凡只需
``I₁/I₂/q₀`` 的调用方（如 ``HBClosure``）走 ``closure_integrals_batch``；需要 ``u``
的调用方读 ``solve_fixed_G_batch``。详见 task-4-report.md §10。

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

        ⚠️ 不要用 ``u.mean()``：格边网格含 ỹ=0 与 ỹ=1 两端点，``mean`` 是 O(h) 偏差
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
        同一闭包求积——``Gb=0`` 时恒有 ``H·ū = I1·G``（即"平均流动度"；``H=1`` 时
        退化为 ``ū = I1·G``，B&F25 图 6 取 ``H=1`` 故写作 ``ū = I1G``），供 B&F25 图 6 用。
        ⚠️ 订正留痕（Task 4.5）：原句写"``ū = I1·G``"只在 ``H=1`` 成立；一般口径由
        ``I1 = H²Ĩ₁`` 与 ``G̃ = H·G`` 给出 ``ū = Ĩ₁G̃ = (I1/H)·G/H·H`` ⇒ ``H·ū = I1·G``
        （数值核验：``H=0.008`` 时 ``ū/(I1·G/H) = 0.99997``，偏离全部来自 ``u_bar`` 的
        trapezoid O(h²) 离散误差；Task 4 报告 §3）。

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

    ``τ̃₂(ỹ) = c₂·ỹ + d₂``（中线带，ỹ∈[0,c̄]，τ̃₂(0)=0）；
    ``τ̃₁(ỹ) = c₁·ỹ + d₁``（壁面带，ỹ∈[c̄,1]，界面剪应力连续；
    ``d₁=(c₂−c₁)c̄``，非共线时 ``d₁`` 不平行于 ``c₁``）。

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

    ⚠️ 轴对齐时返回 ``|c₀·ỹ+d₀|``——与向量化前的标量实现逐位一致。（注：``np.hypot(x, 0.0)``
    与 ``abs(x)`` **严格相等**，所以 ulp 风险**不在**模的计算；需要标量分支的其实是
    **切点求法**（见 :func:`_yield_pieces`）：标量支用一次式 ``(±τ_Y−d₀)/c₀``，
    向量支解二次方程（``_quad_roots``）——两者算术不同。）
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
            #    由壁面 ũ(1)=0 向回累积（中点求积；与格心处的向后差分逐位一致）
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
    """ỹ² 权重（(2.14) 的被积权重；q₀ 已改走通量比口径，不再用它）。"""
    return yv * yv


# --------------------------------------------------------------------------- #
# 闭包量：(2.14)/(2.15) + q₀ 通量比，在解析应力场上求积
# --------------------------------------------------------------------------- #


def _closure_integrals(c_bar, n, kappa_t, tau_y, G_t, Gb_t):
    """tilde 系下的 ``(Ĩ₁, Ĩ₂, q₀)``（(2.14)/(2.15) + q₀ 通量比）。

    逐层积分被 ỹ=c̄ 与屈服面切开；被积函数 ``weight(ỹ)·(1/η̃)``，其中
    ``1/η̃ = γ̇(τ̃)/|τ̃|``（(A3) 有效黏度），未屈服处为 0。
    """
    c2, d2, c1, d1 = _stress_vectors(c_bar, G_t, Gb_t)
    # 中线带（流体 2）：τ̃ = c₂·ỹ + d₂，ỹ∈[0,c̄]
    int2_i1 = _gap_integral(0.0, c_bar, c2, d2, kappa_t[1], n[1], tau_y[1], _weight_y2)
    # 壁面带（流体 1）：τ̃ = c₁·ỹ + d₁，ỹ∈[c̄,1]
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
    """``∫_lo^hi weight(ỹ)·(1/η̃)(ỹ) dỹ``，应力仿射 ``τ̃(ỹ) = c·ỹ + d``（2 维）。

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


# --------------------------------------------------------------------------- #
# Task 4.5：批次轴——n_points 个独立间隙问题一次求解
# --------------------------------------------------------------------------- #
#
# 动机：Phase A 的 ``HBClosure`` 逐格调用 :func:`solve_fixed_G`，呼101 真实场
# （(40,250)=10000 点）单次 ``mobility()`` 实测 82–116 s（Task 4 报告 §4.2）。
# 瓶颈不在"每格一次解"本身，而在**每次解内部的 numpy 调用次数**：``_BISECT_STEPS=60``
# 的逐格二分 × 每个 Uzawa 迭代 × 迭代数 ⇒ 单次解 ~600 次小数组 numpy 调用
# （天花板由调用开销而非计算量决定）。批次版把这些调用换成同一批调用、数组多一条
# 批次轴 ⇒ 调用开销摊薄（见 Task 4.5 报告 §10 的实测加速比）。
#
# 形状契约（与标量路径的 `_as_2vec` 口径一致，歧义处已钉死）
# --------------------------------------------------------
# - ``c_bar``/``H``：标量 或 ``(n_points,)``；
# - ``G``/``Gb``：标量（⇒ 第 1 分量常向量）、``(2,)`` **一律按公共常向量**解释、
#   ``(n_points,)``（⇒ 每点沿第 1 分量）、``(n_points,2)``（⇒ 逐点向量）。
#   ⚠️ ``n_points == 2`` 的批次若要逐点标量 G，请传 ``(2,2)``（``(2,)`` 归公共向量）。
# - ``n``/``kappa``/``tau_y``：``(2,)`` 公共元组 或 ``(n_points,2)`` 逐点。
# - 批次大小 N 由上述"带长度"的输入取最大值，彼此不一致即 ``ValueError``。
#
# 与标量路径的关系（**冻结标量路径**，只新增）
# ------------------------------------------
# 批量核函数与标量逐式对应，且尽量保持**相同的运算顺序**（同为分段 GL-16 求积、
# 同为 60 步二分、同为 ρ=r=1 的 Uzawa、逐点独立收敛即冻结）。两处**已知的末位差**：
# ``‖·‖``/``dot`` 用逐元素式（``sqrt(a²+b²)``/``a*c+b*d``）替代标量路径的
# ``np.linalg.norm``/``np.dot``（BLAS，末位可差 1 ulp）；**轴向时逐位相同**
# （实测 20000/20000）。故本模块的契约是**批量 vs 标量 rel ≤ 1e-12**，不是逐位
# （见 tests/contract/test_gap_solver_batch.py）。
#
# 退化格点（**不中断整批**）
# ------------------------
# 标量路径对"全场未屈服"与"G=Gb=0"逐点抛 ``ValueError``（无闭包）。批次里这两种
# 格点可能**合法地**成片出现（窄边 static wall layer）⇒ 批量 API **不抛错**，改为
# 逐点标记 ``undefined``（其 ``I1/I2/q0/u`` 为 NaN、``converged=False``、``iters=0``），
# 由调用方（如 ``HBClosure`` 的 R-T1-6 地板）自行处置。不收敛/发散**仍抛
# ``RuntimeError``**（与标量同纪律，不静默回退）。


@dataclass(frozen=True)
class ClosureIntegralsBatch:
    """批次闭包量（(2.14)/(2.15)/(2.22) 求积）——**不含** Uzawa 解 ``u``。

    闭包量由**解析应力场** ``τ̃(ỹ)=c·ỹ+d`` 求积得到，与 Uzawa 迭代解无关
    （见模块 docstring 的"实现要点"）。凡只需 ``I₁/I₂/q₀`` 的调用方（如
    ``HBClosure``）走本入口即可省掉迭代。

    Attributes:
        I1/I2/q0: ``(n_points,)``，口径同 :class:`GapSolution`；``undefined`` 处为 NaN。
        undefined: ``(n_points,)`` bool——闭包无定义（``I₁≤0``：全场未屈服，含
            ``G=Gb=0`` 的零驱动格点）。
    """

    I1: Array
    I2: Array
    q0: Array
    undefined: Array


@dataclass(frozen=True)
class GapSolutionBatch:
    """批次间隙解（逐项口径与 :class:`GapSolution` 相同，多一条批次轴）。

    Attributes:
        I1/I2: ``(n_points,)`` 量纲闭包量（= ``GapSolution.I1/I2`` 口径）。
        q0: ``(n_points,)`` 通量比口径 ``q₀``。
        u: ``(n_points, ny)``（**全体轴向**）或 ``(n_points, 2, ny)``（存在非共线点）。
            ``undefined`` 点的整行为 NaN。取末点一律用 ``u[..., -1]``。
        y: ``(ny,)`` 归一化坐标（全批次共用）。
        G: ``(n_points, 2)`` 输入的（量纲、原坐标系）修改压力梯度——同
            :class:`GapSolution` 的 ``G`` 口径。
        iters: ``(n_points,)`` int，各点自身收敛步（``undefined`` 点为 0）。
        converged: ``(n_points,)`` bool；``undefined`` 点为 False，其余恒 True
            （仍有未收敛点则整批抛 ``RuntimeError``，不返回 False）。
        undefined: ``(n_points,)`` bool——闭包无定义（全场未屈服 / ``G=Gb=0``）。

    ``H·ū = I₁·G``（等密度；(2.13)）逐点成立（``ū`` 见 :attr:`u_bar`；离散偏差为
    trapezoid 的 O(h²)，同 :class:`GapSolution`）。
    """

    I1: Array
    I2: Array
    q0: Array
    u: Array
    y: Array
    G: Array
    iters: Array
    converged: Array
    undefined: Array

    @property
    def u_bar(self):
        """间隙平均流速 ū = ``np.trapezoid(u, y, axis=-1)``（**全项目唯一口径**）。

        ``u`` 为 ``(n_points, ny)`` 时返回 ``(n_points,)``；为 ``(n_points, 2, ny)``
        时返回 ``(n_points, 2)``（逐分量）。⚠️ 不要用 ``u.mean(axis=-1)``。
        """
        return np.trapezoid(self.u, self.y, axis=-1)


def closure_integrals_batch(c_bar, n, kappa, tau_y, G, Gb=(0.0, 0.0), *, H=1.0):
    """批次闭包量 ``(I₁, I₂, q₀)``（(2.14)/(2.15)/(2.22)；**不做 Uzawa 迭代**）。

    Args:
        c_bar: ``(n_points,)`` 或标量——间隙平均顶替液体积分数（界面位置 ỹ_i）。
        n: ``(n_points,2)`` 或 ``(2,)``——幂律指数（1=壁面带/被顶替液、2=中线带/顶替液）。
        kappa: ``(n_points,2)`` 或 ``(2,)``——稠度系数。
        tau_y: ``(n_points,2)`` 或 ``(2,)``——屈服应力（≥0）。
        G: 修改压力梯度；``(n_points,2)``/``(n_points,)``/``(2,)``/标量（形状契约见模块 docstring）。
        Gb: 浮力向量，形状规则同 ``G``。
        H: ``(n_points,)`` 或标量——物理半隙。

    Returns:
        :class:`ClosureIntegralsBatch`。

    Raises:
        ValueError: 形状不一致/取值越界（``c̄∉[0,1]``、``H≤0``、``n≤0``、``κ≤0``、
            ``τ_Y<0``）。
    """
    ctx = _batch_inputs(c_bar, n, kappa, tau_y, G, Gb, H)
    I1_t, I2_t, q0 = _closure_integrals_batch(
        ctx["c"], ctx["n"], ctx["kappa_t"], ctx["tau_y"], ctx["Gt"], ctx["Gbt"], ctx["axial"]
    )
    undef = ~(I1_t > 0.0)                      # 标量路径的判据：`not I1_t > 0` ⇒ 抛错
    nan = np.full(I1_t.shape, np.nan)
    return ClosureIntegralsBatch(
        I1=np.where(undef, nan, ctx["H2"] * I1_t),
        I2=np.where(undef, nan, ctx["H3"] * I2_t),
        q0=np.where(undef, nan, q0),
        undefined=undef,
    )


def solve_fixed_G_batch(
    c_bar,
    n,
    kappa,
    tau_y,
    G,
    Gb=(0.0, 0.0),
    *,
    H=1.0,
    ny: int = 201,
    r: float = 1.0,
    tol: float = 1e-10,
    max_iter: int = 2000,
) -> GapSolutionBatch:
    """B&F25 附录 A.2.1 的**批次版**：一次求解 ``n_points`` 个独立间隙问题。

    逐点语义与 :func:`solve_fixed_G` **完全相同**（同离散、同 ``ρ=r=1``、同 60 步二分、
    同 ``p=1+min(n₁,n₂)`` 判据）；差别只在"沿批次轴一次算完"与退化格点处理：

    - 每点用**自己的** ``c̄_i`` ⇒ 各自的界面掩码 ``in_fluid2``、各自的逐格 ``κ̃/n/τ̃_Y``；
    - 逐点独立收敛即**冻结**（其 ``u/q/λ`` 不再更新）⇒ 各点返回的仍是"它自己收敛那一步"
      的解，``iters`` 逐点给出；
    - 退化格点（全场未屈服 / ``G=Gb=0``）标记 ``undefined`` 而不抛错（见模块 docstring）；
    - 任一**有效**格点发散或 ``max_iter`` 内不收敛 ⇒ ``RuntimeError``（不静默）。

    Args:
        c_bar: ``(n_points,)`` 或标量。
        n/kappa/tau_y: ``(n_points,2)`` 或 ``(2,)``（逐相，下标同 :func:`solve_fixed_G`）。
        G/Gb: 形状契约见模块 docstring。
        H: ``(n_points,)`` 或标量。
        ny/r/tol/max_iter: 与 :func:`solve_fixed_G` 同名同义、全批次共用。

    Returns:
        :class:`GapSolutionBatch`。

    Raises:
        ValueError: 形状不一致或取值越界。
        RuntimeError: 有有效格点发散或未在 ``max_iter`` 内收敛（消息含点数与最差残差）。
    """
    ctx = _batch_inputs(c_bar, n, kappa, tau_y, G, Gb, H, ny=ny, r=r, max_iter=max_iter)
    I1_t, I2_t, q0 = _closure_integrals_batch(
        ctx["c"], ctx["n"], ctx["kappa_t"], ctx["tau_y"], ctx["Gt"], ctx["Gbt"], ctx["axial"]
    )
    undef = ~(I1_t > 0.0)
    y, u, iters = _uzawa_fixed_G_batch(
        ctx["c"], ctx["n"], ctx["kappa_t"], ctx["tau_y"], ctx["Gt"], ctx["Gbt"],
        ctx["axial"], undef, ny, r, tol, max_iter,
    )
    if np.all(ctx["axial"]):                 # 全体轴向/共线 ⇒ (n_points, ny)（形状契约）
        u = u[:, 0, :]
    nan = np.full(I1_t.shape, np.nan)
    return GapSolutionBatch(
        I1=np.where(undef, nan, ctx["H2"] * I1_t),
        I2=np.where(undef, nan, ctx["H3"] * I2_t),
        q0=np.where(undef, nan, q0),
        u=u,
        y=y,
        G=ctx["G_raw"],
        iters=iters,
        converged=~undef,
        undefined=undef,
    )


# --------------------------------------------------------------------------- #
# 批次内部：输入规范化
# --------------------------------------------------------------------------- #


def _batch_as_scalar(v, n_pts, name):
    """``(n_points,)`` 或标量 → ``(n_points,)``。"""
    arr = np.asarray(v, dtype=float)
    if arr.ndim == 0:
        return np.full(n_pts, float(arr))
    if arr.ndim == 1 and arr.shape[0] == n_pts:
        return arr.astype(float, copy=False)
    raise ValueError(
        f"{name} 须为标量或 (n_points,)={n_pts} 的数组（批接口径，见模块 docstring），"
        f"得到 shape={arr.shape}"
    )


def _batch_as_pair(v, n_pts, name):
    """``(2,)`` 公共元组 或 ``(n_points,2)`` → ``(n_points,2)``。"""
    arr = np.asarray(v, dtype=float)
    if arr.shape == (2,):
        return np.broadcast_to(arr, (n_pts, 2)).copy()
    if arr.shape == (n_pts, 2):
        return arr
    raise ValueError(
        f"{name} 须为 (2,) 公共元组或 (n_points,2)={n_pts}×2 的数组（批接口径），"
        f"得到 shape={arr.shape}"
    )


def _batch_as_vector(v, n_pts, name):
    """标量 / ``(2,)`` / ``(n_points,)`` / ``(n_points,2)`` → ``(n_points,2)``。

    ⚠️ ``(2,)`` **一律**按公共常向量解释（与标量路径 ``_as_2vec`` 口径一致）；
    ``n_points == 2`` 的批次要逐点标量请传 ``(2,2)``。
    """
    arr = np.asarray(v, dtype=float)
    if arr.ndim == 0:
        out = np.zeros((n_pts, 2))
        out[:, 0] = float(arr)
        return out
    if arr.shape == (2,):
        return np.broadcast_to(arr, (n_pts, 2)).copy()
    if arr.ndim == 1 and arr.shape[0] == n_pts:
        out = np.zeros((n_pts, 2))
        out[:, 0] = arr
        return out
    if arr.shape == (n_pts, 2):
        return arr
    raise ValueError(
        f"{name} 须为标量、(2,) 公共向量、(n_points,)={n_pts} 或 (n_points,2)；"
        f"得到 shape={arr.shape}"
    )


def _batch_n_points(c_bar, n, kappa, tau_y, G, Gb, H) -> int:
    """由输入形状定出批次大小 N（带长度的输入须彼此一致）。"""
    cands = []
    for name, v, is_pair in (("c_bar", c_bar, False), ("H", H, False),
                             ("n", n, True), ("kappa", kappa, True),
                             ("tau_y", tau_y, True), ("G", G, True), ("Gb", Gb, True)):
        arr = np.asarray(v, dtype=float)
        if arr.ndim == 2:
            if arr.shape[1] != 2:
                raise ValueError(
                    f"{name} 的第二维须为 2（逐相/二维向量），得到 shape={arr.shape}"
                )
            cands.append((name, arr.shape[0]))
        elif arr.ndim == 1 and not (is_pair and arr.shape == (2,)):
            cands.append((name, arr.shape[0]))
        elif arr.ndim > 2:
            raise ValueError(f"{name} 的维度须 ≤ 2，得到 shape={arr.shape}")
    if not cands:
        return 1
    n_pts = cands[0][1]
    for name, n_i in cands:
        if n_i != n_pts:
            raise ValueError(
                f"批接口径：各输入的批次长度须一致，得到 {name}={n_i} 与 {cands[0][0]}={n_pts}"
            )
    return n_pts


def _batch_inputs(c_bar, n, kappa, tau_y, G, Gb, H, *, ny=None, r=None, max_iter=None) -> dict:
    """规范化批量输入 → 逐点数组 + 方向归约 + (A4) tilde 换算。

    Returns:
        dict：``c``(N,)、``n``(N,2)、``kappa_t``(N,2)、``tau_y``(N,2)、``Gt``/``Gbt``
        (N,2)（tilde 系）、``axial``(N,) bool、``G_raw``(N,2)（原坐标系输入的 G，供回显）、
        ``H2``/``H3``(N,)（``H²``/``H³`` 量纲因子）。
    """
    n_pts = _batch_n_points(c_bar, n, kappa, tau_y, G, Gb, H)
    if ny is not None and (not isinstance(ny, int) or isinstance(ny, bool) or ny < 3):
        raise ValueError(f"ny 须为 ≥3 的整数（格心数 = ny−1），得到 ny={ny!r}")
    c = _batch_as_scalar(c_bar, n_pts, "c_bar")
    Hv = _batch_as_scalar(H, n_pts, "H")
    n_arr = _batch_as_pair(n, n_pts, "n")
    kap = _batch_as_pair(kappa, n_pts, "kappa")
    tau = _batch_as_pair(tau_y, n_pts, "tau_y")
    G_raw = _batch_as_vector(G, n_pts, "G")
    Gb_raw = _batch_as_vector(Gb, n_pts, "Gb")
    if r is not None and not (float(r) > 0.0):
        raise ValueError(f"r 须为正（Uzawa 数值参数），得到 r={r!r}")
    if max_iter is not None and (not isinstance(max_iter, int)
                                 or isinstance(max_iter, bool) or max_iter < 1):
        raise ValueError(f"max_iter 须为 ≥1 的整数，得到 max_iter={max_iter!r}")
    # 逐点取值校验（与 _validate_inputs 同判据，逐点化）
    _check_range(c, 0.0, 1.0, "c_bar")
    _check_positive(Hv, "H")
    _check_positive(n_arr, "n")
    _check_positive(kap, "kappa")
    _check_nonneg(tau, "tau_y")

    # 方向归约（逐点；共线 ⇒ 旋转到共同方向，非共线 ⇒ 原样）——与 _reduce_directions 逐式同
    Gx, Gy = G_raw[:, 0], G_raw[:, 1]
    Gbx, Gby = Gb_raw[:, 0], Gb_raw[:, 1]
    G_mag = np.sqrt(Gx * Gx + Gy * Gy)      # 轴向时与 np.linalg.norm 逐位同（实测）
    Gb_mag = np.sqrt(Gbx * Gbx + Gby * Gby)
    safe_G = np.where(G_mag > 0.0, G_mag, 1.0)
    safe_Gb = np.where(Gb_mag > 0.0, Gb_mag, 1.0)
    ex = np.where(G_mag > 0.0, Gx / safe_G, np.where(Gb_mag > 0.0, Gbx / safe_Gb, 0.0))
    ey = np.where(G_mag > 0.0, Gy / safe_G, np.where(Gb_mag > 0.0, Gby / safe_Gb, 0.0))
    Gbs = Gbx * ex + Gby * ey
    perp_x = Gbx - Gbs * ex
    perp_y = Gby - Gbs * ey
    perp = np.sqrt(perp_x * perp_x + perp_y * perp_y)
    collinear = perp <= 1e-9 * np.maximum(1.0, Gb_mag)
    Gt = np.where(collinear[:, None], np.stack([G_mag, np.zeros(n_pts)], axis=1), G_raw)
    Gbt = np.where(collinear[:, None], np.stack([Gbs, np.zeros(n_pts)], axis=1), Gb_raw)

    # (A4) 换到 tilde 系：κ̃ = κ/H^n、G̃ = H·G（与标量路径同一乘法顺序）
    kappa_t = kap / Hv[:, None] ** n_arr
    Gt = Hv[:, None] * Gt
    Gbt = Hv[:, None] * Gbt
    axial = (Gt[:, 1] == 0.0) & (Gbt[:, 1] == 0.0)
    return {"c": c, "n": n_arr, "kappa_t": kappa_t, "tau_y": tau, "Gt": Gt, "Gbt": Gbt,
            "axial": axial, "G_raw": G_raw, "H2": Hv ** 2, "H3": Hv ** 3, "n_pts": n_pts}


def _check_range(v, lo, hi, name) -> None:
    arr = np.asarray(v, dtype=float)
    if not np.all(np.isfinite(arr)) or np.any(arr < lo) or np.any(arr > hi):
        raise ValueError(f"{name} 须在 [{lo},{hi}] 内的有限值（gap_solver 定义域）")


def _check_positive(v, name) -> None:
    arr = np.asarray(v, dtype=float)
    if not np.all(np.isfinite(arr)) or np.any(arr <= 0.0):
        raise ValueError(f"{name} 须为正的有限值")


def _check_nonneg(v, name) -> None:
    arr = np.asarray(v, dtype=float)
    if not np.all(np.isfinite(arr)) or np.any(arr < 0.0):
        raise ValueError(f"{name} 须为非负的有限值")


# --------------------------------------------------------------------------- #
# 批次内部：闭包求积（(2.14)/(2.15)/(2.22)，逐点解析应力场）
# --------------------------------------------------------------------------- #

# 每层切段上限：非轴向 4 个候选切点（屈服面 2 + 零应力面 2）⇒ 5 段
_BATCH_SEGMENTS = 5


def _stress_vectors_batch(c, Gt, Gbt):
    """逐点应力仿射系数 ``(c₂, d₂, c₁, d₁)``（(A1)/(A2)/(A17)；与 ``_stress_vectors`` 逐式同）。"""
    c2 = Gt - (1.0 - c)[:, None] * Gbt
    c1 = Gt + c[:, None] * Gbt
    return c2, np.zeros_like(c2), c1, (c2 - c1) * c[:, None]


def _quad_roots_batch(A, B, C):
    """``A y² + B y + C = 0`` 的逐点实根（``(N,2)``；无根/退化处填 ``+inf``）。

    与 ``_quad_roots`` 逐分支对应：``A ≤ 0`` ⇒ 退化为线性（``B=0`` 时无根）；
    ``disc < 0`` ⇒ 无根。无根填 ``+inf`` ⇒ 后续 clip 到区间右端 ⇒ **零宽段**，
    与标量路径"该切点不存在"等价（零宽段贡献恒 0）。
    """
    linear = A <= 0.0
    safe_B = np.where(B != 0.0, B, 1.0)
    lin_root = np.where(B != 0.0, -C / safe_B, np.inf)
    disc = B * B - 4.0 * A * C
    has = (~linear) & (disc >= 0.0)
    r = np.sqrt(np.maximum(disc, 0.0))
    safe_A = np.where(A > 0.0, A, 1.0)
    q1 = (-B - r) / (2.0 * safe_A)
    q2 = (-B + r) / (2.0 * safe_A)
    return np.stack([np.where(linear, lin_root, np.where(has, q1, np.inf)),
                     np.where(linear, np.inf, np.where(has, q2, np.inf))], axis=1)


def _layer_edges_batch(lo, hi, cv, dv, tau_y, axial):
    """逐层切段端点 ``(N, _BATCH_SEGMENTS+1)``（升序、含端点）。

    轴向点用一次式切点 ``(±τ_Y−d₀)/c₀``、``−d₀/c₀``（与标量 ``_yield_pieces`` 的标量支
    逐式相同）；非轴向点用二次求根（``|c·ỹ+d|²=τ_Y²`` 与 ``=0``）。两套切换点都算、
    按 ``axial`` 选一路（避免分支循环）。切点 clip 到 ``[lo,hi]`` ⇒ 越界/重复切点
    退化为**零宽段**（贡献 0），与标量路径的"筛掉区间外切点"等价。
    """
    n_pts = lo.shape[0]
    c0, d0 = cv[:, 0], dv[:, 0]
    with np.errstate(divide="ignore", invalid="ignore"):
        c_safe = np.where(c0 != 0.0, c0, 1.0)
        ax_pts = np.stack([(tau_y - d0) / c_safe, (-tau_y - d0) / c_safe, -d0 / c_safe,
                           np.full(n_pts, np.inf)], axis=1)
        # 标量支：c₀ == 0 ⇒ 无切点（`cuts = [] if q_ == 0.0`）
        no_cut = np.broadcast_to((c0 == 0.0)[:, None], (n_pts, 4))
        ax_pts = np.where(no_cut, np.inf, ax_pts)
    A = (cv * cv).sum(axis=1)
    B = 2.0 * (cv * dv).sum(axis=1)
    dv2 = (dv * dv).sum(axis=1)
    vec_pts = np.concatenate([_quad_roots_batch(A, B, dv2 - tau_y * tau_y),
                              _quad_roots_batch(A, B, dv2)], axis=1)
    cuts = np.where(axial[:, None], ax_pts, vec_pts)
    cuts = np.clip(cuts, lo[:, None], hi[:, None])
    edges = np.sort(np.concatenate([lo[:, None], cuts, hi[:, None]], axis=1), axis=1)
    return edges


def _tau_abs_batch(cv, dv, yq, axial):
    """``|τ̃(ỹ)|``（``(N,16)``）：轴向 = ``|c₀+d₀|``；非轴向 = ``hypot``。"""
    x = cv[:, 0:1] * yq + dv[:, 0:1]
    if np.all(axial):
        return np.abs(x)
    y = cv[:, 1:2] * yq + dv[:, 1:2]
    return np.where(axial[:, None], np.abs(x), np.hypot(x, y))


def _inv_eff_viscosity_batch(tau_abs, kappa, n_pow, tau_y):
    """逐点 ``1/η̃ = γ̇/|τ̃|``（(A3)）；未屈服取 0（与 ``_inv_eff_viscosity`` 同式）。"""
    yielded = tau_abs > tau_y
    base = np.maximum(tau_abs - tau_y, 0.0)               # 未屈服处取 0（避免负底开方）
    gamma_dot = (base / kappa) ** (1.0 / n_pow)
    return np.where(yielded, gamma_dot / np.where(yielded, tau_abs, 1.0), 0.0)


def _gap_integral_batch(lo, hi, cv, dv, kappa, n_pow, tau_y, axial, weight):
    """``∫_lo^hi weight()·(1/η̃) dỹ``（``(N,)``）——逐点区间 + 逐点应力仿射。

    与 ``_gap_integral`` 同结构：按屈服面/零应力面切段 + 逐段 GL-16；求和顺序为
    段序升序（零宽段贡献恒 0 ⇒ 与标量"段不存在"等价）。
    """
    edges = _layer_edges_batch(lo, hi, cv, dv, tau_y, axial)
    total = np.zeros(lo.shape[0])
    for s in range(_BATCH_SEGMENTS):
        a, b = edges[:, s], edges[:, s + 1]
        span = 0.5 * (b - a)
        yq = span[:, None] * _GL_X[None, :] + 0.5 * (a + b)[:, None]
        inv_eta = _inv_eff_viscosity_batch(_tau_abs_batch(cv, dv, yq, axial),
                                           kappa[:, None], n_pow[:, None], tau_y[:, None])
        seg = span * np.sum(_GL_W[None, :] * weight(yq) * inv_eta, axis=1)
        total = total + np.where(b > a, seg, 0.0)
    return total


def _gap_integral_vec_batch(lo, hi, cv, dv, kappa, n_pow, tau_y, axial, weight):
    """``∫ weight·(1/η̃)·τ̃ dỹ``（``(N,2)``）——``_gap_integral_vec`` 的批版本。"""
    edges = _layer_edges_batch(lo, hi, cv, dv, tau_y, axial)
    total = np.zeros_like(cv)
    for s in range(_BATCH_SEGMENTS):
        a, b = edges[:, s], edges[:, s + 1]
        span = 0.5 * (b - a)
        yq = span[:, None] * _GL_X[None, :] + 0.5 * (a + b)[:, None]
        tau_vec = cv[:, :, None] * yq[:, None, :] + dv[:, :, None]
        inv_eta = _inv_eff_viscosity_batch(_tau_abs_batch(cv, dv, yq, axial),
                                           kappa[:, None], n_pow[:, None], tau_y[:, None])
        seg = span[:, None] * np.sum(_GL_W[None, None, :] * weight(yq)[:, None, :]
                                     * (inv_eta[:, None, :] * tau_vec), axis=-1)
        total = total + np.where((b > a)[:, None], seg, 0.0)
    return total


def _flux_ratio_q0_batch(c, c2, d2, c1, d1, kappa_t, n_pow, tau_y, axial):
    """``q0 = |∫₀^{c̄} ũ dỹ| / |∫₀¹ ũ dỹ|``（逐点；``_flux_ratio_q0`` 的批版本）。"""
    w_y = lambda yv: yv                                       # noqa: E731
    w_one = lambda yv: np.ones_like(yv)                       # noqa: E731
    den2 = _gap_integral_vec_batch(np.zeros_like(c), c, c2, d2, kappa_t[:, 1],
                                   n_pow[:, 1], tau_y[:, 1], axial, w_y)
    den1 = _gap_integral_vec_batch(c, np.ones_like(c), c1, d1, kappa_t[:, 0],
                                   n_pow[:, 0], tau_y[:, 0], axial, w_y)
    int1_one = _gap_integral_vec_batch(c, np.ones_like(c), c1, d1, kappa_t[:, 0],
                                       n_pow[:, 0], tau_y[:, 0], axial, w_one)
    num = den2 + c[:, None] * int1_one
    den = den2 + den1
    num_abs = _vec_abs_batch(num, axial)
    den_abs = _vec_abs_batch(den, axial)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(den_abs > 0.0, num_abs / np.where(den_abs > 0.0, den_abs, 1.0), 0.0)


def _vec_abs_batch(v, axial):
    """``(N,2,...)`` 的逐点逐格模 ``|v|``：轴向点 = ``|v|``，非轴向 = ``hypot``。

    （``hypot(x,0) == |x|`` 严格成立 ⇒ 轴向点选哪一支都逐位相同；此处的轴向掩码只是
    与标量 ``_vec_abs`` 的结构对应。）
    """
    if np.all(axial):
        return np.abs(v[:, 0])
    mask = axial.reshape((axial.shape[0],) + (1,) * (v.ndim - 2))
    return np.where(mask, np.abs(v[:, 0]), np.hypot(v[:, 0], v[:, 1]))


def _closure_integrals_batch(c, n_arr, kappa_t, tau_y, Gt, Gbt, axial):
    """逐点 ``(Ĩ₁, Ĩ₂, q₀)``（tilde 系；含端点/屈服面切段求积）。"""
    c2, d2, c1, d1 = _stress_vectors_batch(c, Gt, Gbt)
    zero, one = np.zeros_like(c), np.ones_like(c)
    # 中线带（流体 2，ỹ∈[0,c̄]）与壁面带（流体 1，∈[c̄,1]）——下标 1/0 对应 κ/n/τ_Y
    int2_i1 = _gap_integral_batch(zero, c, c2, d2, kappa_t[:, 1], n_arr[:, 1],
                                  tau_y[:, 1], axial, _weight_y2)
    int1_i1 = _gap_integral_batch(c, one, c1, d1, kappa_t[:, 0], n_arr[:, 0],
                                  tau_y[:, 0], axial, _weight_y2)
    int1_i2 = _gap_integral_batch(c, one, c1, d1, kappa_t[:, 0], n_arr[:, 0],
                                  tau_y[:, 0], axial,
                                  lambda yv: yv * (1.0 - yv))
    I1_t = int2_i1 + int1_i1
    I2_t = (1.0 - c) * int2_i1 + c * int1_i2
    q0 = _flux_ratio_q0_batch(c, c2, d2, c1, d1, kappa_t, n_arr, tau_y, axial)
    return I1_t, I2_t, q0


# --------------------------------------------------------------------------- #
# 批次内部：Uzawa（A.2.1，逐点独立 + 收敛冻结）
# --------------------------------------------------------------------------- #


def _pnorm_batch(v, p):
    """逐点离散 p 范数 ``(mean|v_i|^p)^{1/p}``（``(N,)``；``_pnorm`` 的批版本）。

    对第 2 维以后**全部**求均值（与 ``_pnorm`` 对 ``(2,n_cell)``/``(n_cell,)`` 的
    ``np.mean`` 同口径）。⚠️ 标量路径先经 ``_vec_abs`` 把逐格 2-向量取模再求 p 范数
    （``_pnorm(_vec_abs(...))``）——调用方也必须先取模，否则判据口径不同（会改变
    收敛步数）。
    """
    arr = np.abs(np.asarray(v, dtype=float))
    axes = tuple(range(1, arr.ndim))
    expo = p.reshape((p.shape[0],) + (1,) * (arr.ndim - 1))
    return np.mean(arr ** expo, axis=axes) ** (1.0 / p)


def _integrate_from_wall_batch(du_dy, h, ny):
    """由壁面 ``u(1)=0`` 向回累积（``(N,2,ny-1) → (N,2,ny)``；同 ``_integrate_from_wall``）。"""
    n_pts = du_dy.shape[0]
    return np.concatenate(
        [h * np.cumsum(du_dy[..., ::-1], axis=-1)[..., ::-1], np.zeros((n_pts, 2, 1))],
        axis=-1,
    )


def _local_q_min_batch(m, kappa_c, n_c, tauy_c, r, axial):
    """(A20)/(A29) 逐点逐格局部极小化（只对屈服格做二分；未屈服取 0）。"""
    q_new = np.zeros_like(m)
    abs_m = _vec_abs_batch(m, axial)                       # (N, n_cell)：模是**逐格**量
    yielded = abs_m > tauy_c                               # (N, n_cell)
    if np.any(yielded):
        # 二分按**屈服格**压缩（与标量 `abs_m[yielded]` 同一压缩口径）；θ 逐格一个，
        # 写回时对该格两个分量同乘（= 标量 `q_new[:, yielded] = theta * m[:, yielded]`）。
        flat = yielded.reshape(-1)
        theta_cell = np.zeros_like(abs_m)
        theta_cell.reshape(-1)[flat] = _bisect_theta(
            abs_m.reshape(-1)[flat],
            np.broadcast_to(kappa_c, abs_m.shape).reshape(-1)[flat],
            np.broadcast_to(n_c, abs_m.shape).reshape(-1)[flat],
            np.broadcast_to(tauy_c, abs_m.shape).reshape(-1)[flat],
            r,
        )
        q_new = theta_cell[:, None, :] * m
    return q_new


def _uzawa_fixed_G_batch(c, n_arr, kappa_t, tau_y, Gt, Gbt, axial, undef, ny, r, tol, max_iter):
    """(A16)—(A21) 批量化 Uzawa（ρ=r=1）；逐点独立收敛即冻结。

    Returns:
        ``(y, u, iters)``：``u`` 为 ``(N,2,ny)``（调用方按需裁成 ``(N,ny)``）；
        ``undef`` 点为 NaN、``iters=0``。
    """
    n_pts = c.shape[0]
    n_cell = ny - 1
    h = 1.0 / n_cell
    y_edge = np.linspace(0.0, 1.0, ny)
    y_cen = (np.arange(n_cell) + 0.5) * h

    in_fluid2 = y_cen[None, :] < c[:, None]                      # (N, n_cell)
    kappa_c = np.where(in_fluid2, kappa_t[:, 1][:, None], kappa_t[:, 0][:, None])
    n_c = np.where(in_fluid2, n_arr[:, 1][:, None], n_arr[:, 0][:, None])
    tauy_c = np.where(in_fluid2, tau_y[:, 1][:, None], tau_y[:, 0][:, None])

    c2, d2, c1, d1 = _stress_vectors_batch(c, Gt, Gbt)
    lam0 = np.where(in_fluid2[:, None, :],
                    c2[:, :, None] * y_cen[None, None, :],
                    c1[:, :, None] * y_cen[None, None, :] + d1[:, :, None])
    p_norm = 1.0 + np.minimum(n_arr[:, 0], n_arr[:, 1])          # (A.2.3)

    q = np.zeros((n_pts, 2, n_cell))
    lam_t = np.zeros((n_pts, 2, n_cell))
    u = np.zeros((n_pts, 2, ny))
    active = ~undef
    iters = np.zeros(n_pts, dtype=int)
    res_u = res_q = res_lam = np.zeros(n_pts)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        for k in range(1, max_iter + 1):
            m_ = lam0 + r * q
            q_new = _local_q_min_batch(m_, kappa_c, n_c, tauy_c, r, axial)

            du_dy = q - lam_t / r
            u_new = _integrate_from_wall_batch(du_dy, h, ny)
            lam_new = lam_t + (du_dy - q_new)

            res_u = _pnorm_batch(_vec_abs_batch(u_new - u, axial), p_norm)
            res_q = _pnorm_batch(_vec_abs_batch(q_new - q, axial), p_norm)
            res_lam = _pnorm_batch(_vec_abs_batch(lam_new, axial), p_norm)

            finite = np.isfinite(res_u) & np.isfinite(res_q) & np.isfinite(res_lam)
            if np.any(active & ~finite):
                bad = int(np.count_nonzero(active & ~finite))
                raise RuntimeError(
                    f"Uzawa（批量定 G）数值发散：{bad}/{n_pts} 个有效格点出现非有限残差"
                    f"（第 {k} 步）：‖Δũ‖_p、‖Δq‖_p、‖λ̃‖_p 中出现 NaN/inf。"
                    "补救方向，与标量路径相同：取 r≈1（r≲0.5 不收敛）并增大 max_iter。"
                )
            # 活跃点推进；本步刚收敛的点保留"收敛那一步"的值（与标量 return 语义等价）
            step = active[:, None, None]
            u = np.where(step, u_new, u)
            q = np.where(step, q_new, q)
            lam_t = np.where(step, lam_new, lam_t)
            now_done = active & (res_u < tol) & (res_q < tol) & (res_lam < tol)
            iters[now_done] = k
            active = active & ~now_done
            if not np.any(active):
                break
    if np.any(active):
        bad = int(np.count_nonzero(active))
        raise RuntimeError(
            f"Uzawa（批量定 G）在 max_iter={max_iter} 步内未收敛：{bad}/{n_pts} 个有效格点"
            f"（最差 ‖Δũ‖_p={np.max(res_u[active]):.3e}, ‖Δq‖_p={np.max(res_q[active]):.3e}, "
            f"‖λ̃‖_p={np.max(res_lam[active]):.3e}），tol={tol:.1e}"
            "（三分量须同时低于 tol；补救方向是取 r≈1 并增大 max_iter）"
        )
    if np.any(undef):
        u = np.where(undef[:, None, None], np.nan, u)
    return y_edge, u, iters