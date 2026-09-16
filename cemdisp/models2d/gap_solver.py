"""一维间隙 Herschel–Bulkley 弱解求解器（B&F25 附录 A；A-1a 定 G 模式）。

``solve_fixed_G`` 求解附录 A.2.1 的两流体层 1D 间隙问题（增广拉格朗日 / Uzawa），
并返回 (2.13)—(2.15)、(2.22) 定义的闭包量 ``I₁ / I₂ / q₀``。它是 Phase A 的
``HBClosure``（Task 4）的底层求解器；定均速模式（A.2.2）另在 Task 3。

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
- ``q0`` 无量纲，按 (2.22) 归一：``q0(0)=0``、``q0(1)=1``。
- ``u``/``y``：同一套格边网格（都长 ``ny``，``y = linspace(0,1,ny)``）；``u`` 是
  tilde 系速度 ũ(ỹ)，即量纲速度在归一化位置处的取值。

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

适用范围（本骨架的边界，勿静默外推）
------------------------------------
附录 A 的 (A1)/(A2)/(A16)—(A21) 是**标量（共线）**形式。本实现把 ``G``/``Gb``
投影到共同方向后按标量求解；``G`` 与 ``Gb`` **非共线**时抛 ``ValueError``
（向量版、以及非共线时 q₀ 的严格推广，留待后续任务）。全场未屈服（处处
``|τ̃|≤τ̃_Y``）⇒ 间隙无流动、``I₁=0``、``q₀`` 无定义，本实现抛 ``ValueError``。

References
----------
Bararpour & Frigaard (2025), JFM 1022 A15：附录 A.2.1（(A16)—(A21)）、A.2.3
（交错网格、ρ=r=1、p 范数判据）；闭包 (2.13)—(2.15)、q₀ (2.22)。
Zhang & Frigaard (2022), JFM 947 A32：(4.21a)/(4.21b)。
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
    """定 G 模式的一维间隙解与闭包量（口径见模块 docstring）。

    Attributes:
        I1: 平均流动度（**量纲**，``H²·Ĩ₁``；= Z&F22 (4.21a) / ``mobility_i1`` 口径）。
        I2: 浮力流动度（**量纲**，``H³·Ĩ₂``；= Z&F22 (4.21b) / ``mobility_i2`` 口径）。
        q0: 各向同性通量函数（无量纲，(2.22)；``q0(0)=0``、``q0(1)=1``）。
        u: tilde 系速度剖面 ũ(ỹ)，格边网格（长 ``ny``），``u(1)=0``。
        y: 归一化坐标 ỹ∈[0,1]（``linspace(0,1,ny)``）。
        G: 输入的修改压力梯度（归一化为 2-元组，可直接回传本函数）。
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
    # 共线标量归约（附录 A 为标量式）：Gs = |G|、Gbs = Gb 在 Ĝ 上的有号投影
    Gs, Gbs = _collinear_scalars(G_vec, Gb_vec)

    # (A4) 换到 tilde 系（H=1）
    kappa_t = (kappa[0] / H**n[0], kappa[1] / H**n[1])
    G_t, Gb_t = H * Gs, H * Gbs

    # 闭包量：由解析应力场 quadrature（(2.14)/(2.15)/(2.22)）
    I1_t, I2_t, q0 = _gap_closures(c_bar, n, kappa_t, tau_y, G_t, Gb_t)

    # 定 G 弱解（(A16)—(A21)）
    y, u, iters = _uzawa_fixed_G(c_bar, n, kappa_t, tau_y, G_t, Gb_t, ny, r, tol, max_iter)

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


def _collinear_scalars(G_vec: Array, Gb_vec: Array) -> tuple[float, float]:
    """把 ``(G, Gb)`` 归约到附录 A 的共线标量 ``(Gs, Gbs)``。

    方向 ``e = Ĝ``（``G=0`` 时取 ``Ĝb``）；``Gs = |G|``、``Gbs = Gb·e``（有号）。
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
    perp = Gb_vec - float(np.dot(Gb_vec, e)) * e
    if float(np.linalg.norm(perp)) > 1e-9 * max(1.0, Gb_mag):
        raise ValueError(
            "本骨架只实现附录 A 的共线标量情形（(A1)/(A2)/(A16)—(A21) 均为标量式）："
            f"G={tuple(G_vec)} 与 Gb={tuple(Gb_vec)} 不共线 ⇒ 需向量版（后续任务）"
        )
    return G_mag, float(np.dot(Gb_vec, e))


# --------------------------------------------------------------------------- #
# 定 G 模式：Uzawa（A.2.1）
# --------------------------------------------------------------------------- #


def _stress_slopes(c_bar: float, G_t: float, Gb_t: float) -> tuple[float, float]:
    """逐层应力斜率（(A1)/(A2)/(A17)，**符号见下**）。

    ``τ̃₂(ỹ) = a₂·ỹ``（中线带，ỹ∈[0,c̄]，τ̃₂(0)=0）；
    ``τ̃₁(ỹ) = (a₂−a₁)·c̄ + a₁·ỹ``（壁面带，ỹ∈[c̄,1]，界面剪应力连续）。

    ⚠️ 符号约定：``G`` 是驱动方向的修改压力梯度（(2.13) 口径）：
    ``Ī₁>0`` ⇒ 间隙平均流速与 G 同向。故本实现的 λ₀ 取附录 A (A17) 的**相反号**
    （即 ``λ₀ = −λ₀^{(A17)}``）：附录 A 的 (A1)/(A16) 印刷式与 (2.13) 反号，
    照抄会使 ``G>0`` 时 ``ũ<0``（与 (2.13) 及 Z&F22 口径矛盾）。闭包量
    ``Ĩ₁/Ĩ₂/q₀`` 只依赖 ``|τ̃|``，此选择不影响它们。
    """
    a2 = G_t - (1.0 - c_bar) * Gb_t
    a1 = G_t + c_bar * Gb_t
    return a2, a1


def _uzawa_fixed_G(c_bar, n, kappa_t, tau_y, G_t, Gb_t, ny, r, tol, max_iter):
    """(A16)—(A21) 的 Uzawa 迭代（ρ=r=1），返回 ``(y, u, iters)``。"""
    n_cell = ny - 1
    h = 1.0 / n_cell
    y_edge = np.linspace(0.0, 1.0, ny)
    y_cen = (np.arange(n_cell) + 0.5) * h

    in_fluid2 = y_cen < c_bar  # 中线带（顶替液）
    kappa_c = np.where(in_fluid2, kappa_t[1], kappa_t[0])
    n_c = np.where(in_fluid2, n[1], n[0])
    tauy_c = np.where(in_fluid2, tau_y[1], tau_y[0])

    a2, a1 = _stress_slopes(c_bar, G_t, Gb_t)
    lam0 = np.where(in_fluid2, a2 * y_cen, (a2 - a1) * c_bar + a1 * y_cen)

    p_norm = 1.0 + min(n[0], n[1])  # (A.2.3)：p = 1 + min(n₁,n₂)
    q = np.zeros(n_cell)
    lam_t = np.zeros(n_cell)  # λ̃^k（λ^k = λ₀ + λ̃^k）
    u = np.zeros(ny)
    for iters in range(1, max_iter + 1):
        # ① (A19)/(A20)：q^{k+1}（局部极小化；|m|≤τ_Y ⇒ q=0）
        m = lam0 + r * q
        abs_m = np.abs(m)
        q_new = np.zeros(n_cell)
        yielded = abs_m > tauy_c
        if np.any(yielded):
            theta = _bisect_theta(
                abs_m[yielded], kappa_c[yielded], n_c[yielded], tauy_c[yielded], r
            )
            q_new[yielded] = theta * m[yielded]

        # ② (A18)：dũ^{k+1}/d|_格心 = q^k − λ̃^k/r（对称面在 ỹ=0 ⇒ C=0），
        #    由壁面 ũ(1)=0 向回累积（中点求积；与格心处的向后差分逐位一致）
        du_dy = q - lam_t / r
        u_new = np.concatenate([h * np.cumsum(du_dy[::-1])[::-1], np.zeros(1)])

        # ③ (A21)：λ̃^{k+1} = λ̃^k + ρ(dũ^{k+1}/dỹ − q^{k+1})，ρ=1
        lam_new = lam_t + (du_dy - q_new)

        res_u = _pnorm(u_new - u, p_norm)
        res_q = _pnorm(q_new - q, p_norm)
        res_lam = _pnorm(lam_new, p_norm)
        u, q, lam_t = u_new, q_new, lam_new
        if res_u < tol and res_q < tol and res_lam < tol:
            return y_edge, u, iters

    raise RuntimeError(
        f"Uzawa（定 G）在 max_iter={max_iter} 步内未收敛："
        f"‖Δũ‖_p={res_u:.3e}, ‖Δq‖_p={res_q:.3e}, ‖λ̃‖_p={res_lam:.3e}, tol={tol:.1e}"
        "（三分量须同时低于 tol；可增大 max_iter 或 r）"
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


def _gap_closures(c_bar, n, kappa_t, tau_y, G_t, Gb_t):
    """tilde 系下的 ``(Ĩ₁, Ĩ₂, q₀)``（(2.14)/(2.15)/(2.22)）。

    逐层积分被 ỹ=c̄ 与屈服面切开；被积函数 ``weight(ỹ)·(1/η̃)``，其中
    ``1/η̃ = γ̇(τ̃)/|τ̃|``（(A3) 有效黏度），未屈服处为 0。
    """
    a2, a1 = _stress_slopes(c_bar, G_t, Gb_t)
    # 中线带（流体 2）：τ̃ = a₂·ỹ，ỹ∈[0,c̄]
    int2_i1 = _gap_integral(0.0, c_bar, 0.0, a2, kappa_t[1], n[1], tau_y[1], _weight_y2)
    # 壁面带（流体 1）：τ̃ = (a₂−a₁)c̄ + a₁·ỹ，ỹ∈[c̄,1]
    p1 = (a2 - a1) * c_bar
    int1_i1 = _gap_integral(c_bar, 1.0, p1, a1, kappa_t[0], n[0], tau_y[0], _weight_y2)
    int1_i2 = _gap_integral(
        c_bar, 1.0, p1, a1, kappa_t[0], n[0], tau_y[0], lambda yv: yv * (1.0 - yv)
    )
    int1_q0 = _gap_integral(c_bar, 1.0, p1, a1, kappa_t[0], n[0], tau_y[0], lambda yv: yv)

    I1_t = int2_i1 + int1_i1
    if not I1_t > 0.0:
        raise ValueError(
            "全场未屈服（|τ̃(ỹ)| ≤ τ̃_Y 处处成立）：间隙无流动 ⇒ Ĩ₁=0、q₀ 无定义"
            "（static plug 情形，须由调用方另行处理）"
        )
    I2_t = (1.0 - c_bar) * int2_i1 + c_bar * int1_i2
    q0 = (int2_i1 + c_bar * int1_q0) / I1_t
    return I1_t, I2_t, q0


def _gap_integral(lo, hi, p_, q_, kappa, n, tau_y, weight) -> float:
    """``∫_lo^hi weight()·(1/η̃)(ỹ) dỹ``，应力仿射 ``τ̃ = p_ + q_·ỹ``。

    按屈服面（``|τ̃|=τ̃_Y``）与零应力点（``τ̃=0``）切段，逐段 Gauss–Legendre；
    每段被积函数解析光滑 ⇒ 界面/屈服面无离散误差。
    """
    total = 0.0
    for a, b in _yield_pieces(lo, hi, p_, q_, tau_y):
        yq = 0.5 * (b - a) * _GL_X + 0.5 * (a + b)
        inv_eta = _inv_eff_viscosity(np.abs(p_ + q_ * yq), kappa, n, tau_y)
        total += 0.5 * (b - a) * float(np.sum(_GL_W * weight(yq) * inv_eta))
    return total


def _yield_pieces(lo, hi, p_, q_, tau_y) -> list[tuple[float, float]]:
    """把 ``[lo,hi]`` 在屈服面/零应力点处切开（返回可能为空的段列表）。"""
    if not hi > lo:
        return []
    cuts = []
    if q_ != 0.0:
        for root in ((tau_y - p_) / q_, (-tau_y - p_) / q_, -p_ / q_):
            if lo < root < hi:
                cuts.append(root)
    edges = [lo] + sorted(cuts) + [hi]
    return [(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]


def _inv_eff_viscosity(tau_abs, kappa, n, tau_y) -> Array:
    """有效黏度倒数 ``1/η̃ = γ̇/|τ̃|``（(A3)）；未屈服（``|τ̃|≤τ̃_Y``）取 0。"""
    tau = np.asarray(tau_abs, dtype=float)
    out = np.zeros_like(tau)
    yielded = tau > tau_y
    if np.any(yielded):
        gamma_dot = ((tau[yielded] - tau_y) / kappa) ** (1.0 / n)
        out[yielded] = gamma_dot / tau[yielded]
    return out