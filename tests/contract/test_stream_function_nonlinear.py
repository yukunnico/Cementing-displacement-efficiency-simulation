"""Task 5（A-3a）契约测试：椭圆方程非线性外迭代 + 批次"由均速反求 G"工具。

覆盖（brief 的验收项 + controller 裁定的数学前提）：

1. **牛顿极限逐位**（硬验收）：``HBClosure`` 牛顿极限下 ``solve_stream_function_nonlinear``
   与线性 ``solve_stream_function`` 逐位一致，且**首轮即刻返回**（不依赖注入的 G）。
2. **非牛顿不动点自洽**：收敛后 ① 返回的 Ψ 与"当前闭包状态（注入 G）下的线性解"逐位一致；
   ② 逐格闭合关系 ``|ū| = Ī₁(G̃)·G̃`` 成立（残差 = 反求容差）；③ 与**独立路径**
   ``gap_solver.solve_fixed_G``（标量 Uzawa）交叉核对 ``u_bar``（口径链的独立验证）。
3. **失败路径**：``max_outer`` 用尽 ⇒ ``RuntimeError``（不静默回退）；零驱动初值
   （G≡0）⇒ 闭包 R-T1-6 地板 + 单次 ``RuntimeWarning``，且外迭代仍能恢复收敛。
4. **求根工具**（``gap_solver.solve_g_from_mean_velocity_batch``）：流的单调性前提、
   ∥ū∥=0 的退化格（``undefined`` 掩码）、反求 vs 标量路径、I₁ 与
   ``closure_integrals_batch`` 逐位同口径。

⚠️ 适用域：HB 闭包属**扩展应用，未获外部验证**（B&F25 自述无外部验证）；Z&F22/23
的判据与结论仅严格适用于竖直井 + 牛顿流体，不得引作 HB/斜井依据。
"""
from __future__ import annotations

import numpy as np
import pytest

from cemdisp.models2d import gap_solver
from cemdisp.models2d.hb_closure import HBClosure
from cemdisp.models2d.stream_function import (solve_stream_function,
                                              solve_stream_function_nonlinear,
                                              velocity_from_stream_function)

# 非牛顿参数组（探针与测试共用；τ_Y>0 ⇒ 有屈服门槛，反求非平凡）
HB_N = (0.7, 0.8)
HB_KAPPA = (1.4, 0.9)
HB_TAUY = (2.0, 0.5)
ETA1, ETA2 = 1.4, 0.9                      # 牛顿极限参照黏度（κ 的 n=1 对应量）
M_RATIO = ETA1 / ETA2


def _geom(ny=40, nz=5):
    """测试几何（与 test_stream_function.py 同族；H 为 Z&F22 半隙）。"""
    y = np.linspace(0, np.pi * 0.1071, ny)
    phi = y / y[-1]
    H = np.broadcast_to(0.0458 * (1 + 0.3 * np.cos(np.pi * phi))[:, None], (ny, nz)).copy()
    return {"y": y, "phi": phi, "H": H, "b": 2 * H, "s": np.linspace(0, 10.0, nz),
            "hole_mm": np.full((1, nz), 260.0), "od_mm": np.full((1, nz), 168.3)}


def _b_field(g, amp=5.0):
    """竖直井 + c̄ 的 φ-梯度（A32-22 页机制：源项由 b 的梯度驱动）。"""
    c = (0.5 + 0.2 * np.cos(np.pi * g["phi"]))[:, None] * np.ones_like(g["H"])
    chi = 1.0 + 2.0 * c - c ** 3
    b = np.zeros((2,) + g["H"].shape)
    b[0] = amp * chi
    return c, b


def _hb_closure(n=HB_N, kappa=HB_KAPPA, tau_y=HB_TAUY, *, m=M_RATIO, B=0.3, ny_gap=201):
    return HBClosure(n=n, kappa=kappa, tau_y=tau_y, m=m, B=B, ny_gap=ny_gap)


def _newtonian_limit_closure():
    """牛顿极限 HBClosure：``n≡1、τ_y≡0`` ⇒ 结构性短路到 NewtonianClosure；
    κ 在 n=1 时即黏度 ⇒ κ=(η₁,η₂)、m=η₁/η₂ 与线性参照同口径。"""
    return HBClosure(n=(1.0, 1.0), kappa=(ETA1, ETA2), tau_y=(0.0, 0.0),
                     m=M_RATIO, B=0.0)


# --------------------------------------------------------------------------- #
# 1. 牛顿极限逐位（R2 硬验收）
# --------------------------------------------------------------------------- #


def test_newtonian_limit_is_bitwise_equal_to_linear_solver():
    g = _geom()
    c, b = _b_field(g)
    psi_lin = solve_stream_function(g, c, ETA1, ETA2, M_RATIO, b)
    hb = _newtonian_limit_closure()
    # 未注入任何 G：牛顿极限短路不消费 G（首轮即刻返回的直接证据）
    assert hb._G is None
    psi_nl = solve_stream_function_nonlinear(g, c, hb, b, max_outer=1)
    assert np.array_equal(psi_lin, psi_nl)


def test_newtonian_limit_needs_only_one_outer_round():
    """首轮 ΔĪ₁≡0 ⇒ 立即返回（max_outer=1 不得抛错）。"""
    g = _geom()
    c, b = _b_field(g)
    hb = _newtonian_limit_closure()
    psi_nl = solve_stream_function_nonlinear(g, c, hb, b, omega=0.5, tol=1e-12, max_outer=1)
    assert np.all(np.isfinite(psi_nl))
    # 首轮确实走到了"反求 + 注入"（G 已被写入）——却仍返回首轮 Ψ（ΔĪ₁≡0 短路）；
    # 牛顿极限的结构性短路不消费缓存（cache_stats 恒为 0，证据：委托 NewtonianClosure）。
    assert hb._G is not None
    assert hb.cache_stats == {"hits": 0, "misses": 0}


# --------------------------------------------------------------------------- #
# 2. 非牛顿：不动点自洽 + 独立路径交叉核对
# --------------------------------------------------------------------------- #


def _converged_hb_run(max_outer=200, omega=0.5, G0=300.0):
    """默认 tol=1e-10 下实测 48 轮收敛 ⇒ 留 4× 余量（50 轮过紧）。"""
    g = _geom()
    c, b = _b_field(g)
    hb = _hb_closure()
    hb.set_pressure_gradient(G0)
    psi = solve_stream_function_nonlinear(g, c, hb, b, omega=omega, tol=1e-10,
                                          max_outer=max_outer)
    return g, c, b, hb, psi


def test_nonlinear_solution_is_fixed_point_of_current_closure_state():
    """收敛 ⇒ 返回的 Ψ 与"当前注入 G 下的线性解"逐位一致（不动点）。"""
    g, c, b, hb, psi = _converged_hb_run()
    psi_check = solve_stream_function(g, c, ETA1, ETA2, M_RATIO, b, closure=hb)
    assert np.array_equal(psi, psi_check)


def test_nonlinear_velocity_closes_against_injected_gradient():
    """逐格闭合关系 |ū| = Ī₁(G̃)·G̃（残差 = 反求容差量级）。"""
    g, c, b, hb, psi = _converged_hb_run()
    H = g["H"]
    w, v = velocity_from_stream_function(psi, g)
    u_mag = np.hypot(w, v)
    G_inj = np.asarray(hb._G, dtype=float)            # 白盒：读回注入的 G 场（唯一出口）
    I1 = np.asarray(hb.mobility(c, M_RATIO, ETA1, ETA2, H), dtype=float)  # 量纲 I₁
    pred = (I1 / H ** 2) * (H * G_inj)                # tilde: Ī₁·G̃ = (I₁/H²)·(H·G)
    rel = np.max(np.abs(u_mag - pred) / np.abs(pred))
    assert rel < 1e-6, f"闭合关系残差 {rel:.3e} 超容差"


def test_nonlinear_velocity_matches_scalar_uzawa_path():
    """独立路径核对：把逐格 G 送标量 ``solve_fixed_G``（Uzawa）⇒ u_bar ≈ 离散 ū。

    这是**口径链的独立验证**（不复用反求路径）：``solve_fixed_G`` 的 ``u_bar`` 是
    解析求积口径的间隙均值，与 (2.2) 有限差分的离散 ū 只差 Ψ 的离散误差。
    """
    g, c, b, hb, psi = _converged_hb_run()
    H = g["H"]
    w, v = velocity_from_stream_function(psi, g)
    u_mag = np.hypot(w, v)
    G_inj = np.asarray(hb._G, dtype=float)
    i, j = 20, 2                                     # 取一个非退化格点
    sol = gap_solver.solve_fixed_G(float(c[i, j]), HB_N, HB_KAPPA, HB_TAUY,
                                   float(G_inj[i, j]), H=float(H[i, j]), ny=201)
    rel = abs(float(sol.u_bar) - float(u_mag[i, j])) / abs(float(sol.u_bar))
    assert rel < 0.05, f"独立路径相对偏差 {rel:.3e}（口径链或离散误差异常）"


def test_nonlinear_differs_from_newtonian_solution():
    """判别性：非牛顿解的 I₁ 场与牛顿参照显著不同（测试不是恒真）。"""
    g, c, b, hb, psi = _converged_hb_run()
    psi_lin = solve_stream_function(g, c, ETA1, ETA2, M_RATIO, b)
    assert not np.allclose(psi, psi_lin, rtol=1e-6, atol=1e-9)


def test_under_relaxation_reduces_step():
    """ω 的作用：ω=1.0 与 ω=0.5 都收敛，但 I₁ 场路径不同（欠松弛改变步长）。"""
    _, _, _, hb_full, psi_full = _converged_hb_run(omega=1.0)
    _, _, _, hb_half, psi_half = _converged_hb_run(omega=0.5)
    assert np.allclose(psi_full, psi_half, rtol=1e-6, atol=1e-8)   # 同一不动点
    assert not np.array_equal(psi_full, psi_half)                  # 不同路径


# --------------------------------------------------------------------------- #
# 3. 失败路径（不得静默）
# --------------------------------------------------------------------------- #


def test_max_outer_exhausted_raises():
    """max_outer 用尽（tol=0 不可能满足）⇒ RuntimeError，不静默回退牛顿闭包。"""
    g = _geom()
    c, b = _b_field(g)
    hb = _hb_closure()
    hb.set_pressure_gradient(300.0)
    with pytest.raises(RuntimeError, match="max_outer"):
        solve_stream_function_nonlinear(g, c, hb, b, tol=0.0, max_outer=1)


def test_zero_drive_initial_gradient_floors_and_warns_then_recovers():
    """初值 G≡0 ⇒ 闭包 R-T1-6 地板 + 单次 RuntimeWarning；外迭代仍恢复收敛。"""
    g = _geom()
    c, b = _b_field(g)
    hb = _hb_closure()
    hb.set_pressure_gradient(0.0)
    with pytest.warns(RuntimeWarning, match="地板"):
        psi = solve_stream_function_nonlinear(g, c, hb, b, tol=1e-10, max_outer=50)
    assert np.all(np.isfinite(psi))
    # 地板告警即"首轮 G≡0 ⇒ 闭包无定义"的直接证据（HBClosure 每实例只告警一次）；
    # 收敛后注入的是反求出的非零 G（该处已无地板格 ⇒ n_static_wall_points 回到 0）
    assert np.max(np.abs(np.asarray(hb._G, dtype=float))) > 0.0


# --------------------------------------------------------------------------- #
# 4. 求根工具（gap_solver 新增公共入口）
# --------------------------------------------------------------------------- #


def test_inverse_utility_matches_scalar_uzawa_path():
    """反求的 G 送标量 ``solve_fixed_G`` 复现目标 ū（独立实现交叉核对）。"""
    H = 0.008
    for c_bar, n, kap, tau in ((0.45, HB_N, HB_KAPPA, HB_TAUY),
                               (0.20, (1.0, 1.0), (0.058, 0.171), (0.0, 0.0)),
                               (0.62, (0.55, 0.90), (3.1, 0.4), (6.0, 1.0))):
        u_target = 0.02
        inv = gap_solver.solve_g_from_mean_velocity_batch(c_bar, n, kap, tau,
                                                          u_target, H=H)
        assert not bool(inv.undefined[0])
        sol = gap_solver.solve_fixed_G(c_bar, n, kap, tau, float(inv.G[0]), H=H, ny=2001)
        assert abs(float(sol.u_bar) - u_target) / u_target < 1e-5


def test_inverse_utility_i1_is_same_caliber_as_closure_integrals():
    """I₁ 口径一致：与 ``closure_integrals_batch`` 在解出的 G 处**到 1 ulp** 相同。

    ⚠️ 这里不能断言逐位：反求输出的是**量纲 G = g̃/H**，闭包侧再做 (A4) 的
    ``g̃ = H·G`` 回乘——该往返在浮点上并非恒等（实测随机 x 上 ``H*(x/H) == x``
    仅 85.7%），故 1 ulp 级差异是口径**换算**的固有误差、非口径不一致
    （实测 field 上 max rel = 1.2e-16，36/37 点逐位相同）。
    """
    H = 0.008
    inv = gap_solver.solve_g_from_mean_velocity_batch(0.45, HB_N, HB_KAPPA, HB_TAUY,
                                                      np.array([1e-3, 1e-1, 1.0]), H=H)
    for k in range(3):
        cb = gap_solver.closure_integrals_batch(0.45, HB_N, HB_KAPPA, HB_TAUY,
                                                float(inv.G[k]), H=H)
        assert float(inv.I1[k]) == pytest.approx(float(cb.I1[0]), rel=1e-15)


def test_inverse_utility_i1_caliber_on_field_to_one_ulp():
    """场级口径核对（37 个不同 c̄）：反求 I₁ 与闭包求值 max rel ≤ 1e-15。"""
    H = 0.0458
    c = np.linspace(0.05, 0.95, 37)
    inv = gap_solver.solve_g_from_mean_velocity_batch(c, HB_N, HB_KAPPA, HB_TAUY,
                                                      np.linspace(0.5, 6.0, 37), H=H)
    cb = gap_solver.closure_integrals_batch(c, HB_N, HB_KAPPA, HB_TAUY, inv.G, H=H)
    assert np.max(np.abs(inv.I1 - cb.I1) / cb.I1) < 1e-15


def test_inverse_utility_zero_velocity_is_undefined():
    """ū≡0（零驱动格）⇒ ``undefined`` 掩码（无正根），G=0、I₁=NaN，不抛错。"""
    inv = gap_solver.solve_g_from_mean_velocity_batch(np.array([0.3, 0.5]),
                                                      HB_N, HB_KAPPA, HB_TAUY,
                                                      np.array([0.0, 1e-2]), H=0.008)
    assert bool(inv.undefined[0]) and not bool(inv.undefined[1])
    assert float(inv.G[0]) == 0.0 and np.isnan(inv.I1[0])


def test_inverse_utility_rejects_negative_mean_velocity():
    with pytest.raises(ValueError, match="非负"):
        gap_solver.solve_g_from_mean_velocity_batch(0.4, HB_N, HB_KAPPA, HB_TAUY,
                                                    -1e-3, H=0.008)


def test_flow_curve_scale_covariance():
    """标度不变性（模块/函数 docstring 的"标度口径"段所引命题）：

    ``κ' = κ·sⁿ/σ、τ_Y' = τ_Y/σ、ū' = ū/s ⇒ G' = G/σ、I₁' = (σ/s)·I₁``。
    （误用 ``κ' = κ/sⁿ`` 会破坏不变性；``σ = s`` 时该标度**可**由闭包参数
    ``(κ·ŵ^{n−1}, τ_Y/ŵ)`` 共同吸收——即"等价地"Option B，生产口径下
    rel ≈ 1e-14，见 ``test_velocity_scale_alignment_equivalence``。）
    """
    H = 0.008
    base = gap_solver.solve_g_from_mean_velocity_batch(0.45, HB_N, HB_KAPPA, HB_TAUY,
                                                       0.02, H=H)
    for s_scale, sigma in ((10.0, 10.0), (10.0, 1.0), (0.1, 1.0), (3.7, 2.499)):
        kap_s = tuple(k * s_scale ** n / sigma for k, n in zip(HB_KAPPA, HB_N))
        tau_s = tuple(t / sigma for t in HB_TAUY)
        inv = gap_solver.solve_g_from_mean_velocity_batch(
            0.45, HB_N, kap_s, tau_s, 0.02 / s_scale, H=H)
        assert float(inv.G[0]) == pytest.approx(float(base.G[0]) / sigma, rel=1e-10)
        assert float(inv.I1[0]) == pytest.approx(
            float(base.I1[0]) * sigma / s_scale, rel=1e-10)


def test_velocity_scale_alignment_equivalence():
    """接线"等价地"公式钉子（R-T5-1 REVISED；**生产口径** ū_phys = ŵ·ū_mod）：

    同一物理问题两种对齐给出**同一个** ``I₁``：① 物理参数 + 物理速度（推荐，
    velocity_scale=ŵ 乘入）；② 模块口径参数 ``(κ·ŵ^{n−1}, τ_Y/ŵ)`` + 模块速度
    （= 物理速度/ŵ）。模块口径速度比物理**大** 1/ŵ 倍（呼101 量级实测 377×，
    生产锚 ``annulus_d2dga.py:1433``：``w = w_unit·(q_half/π)``），故模块口径的
    应力与 ``τ_Y`` 同比大 1/ŵ 倍、``κ`` 大 ``ŵ^{n−1}`` 倍——无量纲群不变。
    （fix round 1 曾按反向口径误写 ``(κ·ŵ^{1−n}, τ_Y·ŵ)`` 并"证伪"生产约定，
    fix round 2 订正；旧公式作为反例钉在本测试尾部，防止再写反。）
    """
    H = 0.008
    w_hat = 2.653e-3                               # 呼101 量级 q_half/π（生产锚）
    u_phys = 0.02
    u_mod = u_phys / w_hat                         # 生产口径：模块速度 = 物理/ŵ
    base = gap_solver.solve_g_from_mean_velocity_batch(0.45, HB_N, HB_KAPPA, HB_TAUY,
                                                       u_phys, H=H)
    kap_mod = tuple(k * w_hat ** (n - 1.0) for k, n in zip(HB_KAPPA, HB_N))
    tau_mod = tuple(t / w_hat for t in HB_TAUY)
    mod = gap_solver.solve_g_from_mean_velocity_batch(0.45, HB_N, kap_mod, tau_mod,
                                                      u_mod, H=H)
    assert float(mod.I1[0]) == pytest.approx(float(base.I1[0]), rel=1e-10)
    assert float(mod.G[0]) == pytest.approx(float(base.G[0]) / w_hat, rel=1e-10)
    # 反例钉子：fix round 1 的翻转式 (κ·ŵ^{1−n}, τ_Y·ŵ) 在生产口径下必不成立
    # （探针实测 rel ≈ 40.6，O(1) 级错误）——防止"反向警示"再写反方向。
    kap_old = tuple(k * w_hat ** (1.0 - n) for k, n in zip(HB_KAPPA, HB_N))
    tau_old = tuple(t * w_hat for t in HB_TAUY)
    mod_old = gap_solver.solve_g_from_mean_velocity_batch(0.45, HB_N, kap_old, tau_old,
                                                          u_mod, H=H)
    rel_old = abs(float(mod_old.I1[0]) - float(base.I1[0])) / float(base.I1[0])
    assert rel_old > 0.1, f"旧公式 rel={rel_old:.3e} 意外地小，反例钉子失效"


def test_flow_curve_is_monotone_in_gradient():
    """数学前提：``F(g̃) = Ĩ₁(g̃)·g̃`` 在 g̃ 上严格递增（含 n>1 剪切增稠与有屈服）。

    这是"逐格标量求根"赖以成立的单调性前提（与 brief 的"先验证单调性再依赖"一致）。
    """
    H = 0.008
    g_t = np.geomspace(1e-8, 1e4, 600)
    for c_bar, n, kap, tau in ((0.45, HB_N, HB_KAPPA, HB_TAUY),
                               (0.35, (1.5, 1.3), (0.7, 1.1), (0.0, 0.0)),   # 剪切增稠
                               (0.80, (1.2, 0.6), (0.9, 2.2), (0.5, 3.0))):
        cb = gap_solver.closure_integrals_batch(np.full(600, c_bar), n, kap, tau,
                                                g_t / H, H=H)
        F = (cb.I1 / H ** 2) * g_t
        assert np.all(np.diff(F[~cb.undefined]) > 0.0)


# --------------------------------------------------------------------------- #
# 5. Gb 口径显式化 + 斜井/强浮力算例（R-T5-4 / R-T5-5，fix round 1）
# --------------------------------------------------------------------------- #


def test_nonzero_gb_argument_is_rejected_not_silently_ignored():
    """形参 ``Gb ≠ 0`` ⇒ ``NotImplementedError``（显式拒绝，两级入口均不静默忽略）。"""
    g = _geom()
    c, b = _b_field(g)
    hb = _hb_closure()
    hb.set_pressure_gradient(300.0)
    with pytest.raises(NotImplementedError, match="Gb"):
        solve_stream_function_nonlinear(g, c, hb, b, Gb=(0.0, 1.0))
    with pytest.raises(NotImplementedError, match="Gb"):
        gap_solver.solve_g_from_mean_velocity_batch(0.45, HB_N, HB_KAPPA, HB_TAUY,
                                                    0.02, Gb=(1.0, 0.0), H=0.008)


def test_nonfinite_gb_argument_is_rejected():
    """``Gb`` 含 NaN/Inf ⇒ ``ValueError``（非有限值不得进入口径判断）。"""
    g = _geom()
    c, b = _b_field(g)
    hb = _hb_closure()
    hb.set_pressure_gradient(300.0)
    with pytest.raises(ValueError, match="有限"):
        solve_stream_function_nonlinear(g, c, hb, b, Gb=(0.0, np.nan))
    with pytest.raises(ValueError, match="有限"):
        gap_solver.solve_g_from_mean_velocity_batch(0.45, HB_N, HB_KAPPA, HB_TAUY,
                                                    0.02, Gb=(np.inf, 0.0), H=0.008)


def test_closure_with_injected_nonzero_gb_is_rejected():
    """闭包已注入非零 ``Gb`` 而反演按 Gb=0 ⇒ 口径分裂，入口显式拒绝（防静默）。"""
    g = _geom()
    c, b = _b_field(g)
    hb = _hb_closure()
    hb.set_pressure_gradient(300.0, Gb=(0.5, 0.2))
    with pytest.raises(NotImplementedError, match="非零 Gb"):
        solve_stream_function_nonlinear(g, c, hb, b)


def _inclined_b_field(g, beta_deg=45.0, amp_phi=1.0, amp_xi=1.0e6):
    """斜井 b 场：φ-槽（轴向重力）+ ξ-槽（方位浮力，沿 ξ 取半波、两端严格为零）。

    ⚠️ **合成应力算例（非物理量级）**：``amp_xi=1e6`` 约为生产浮力量级
    （b ~ O(10–70)，0708 八井实测）的 1e4 倍——目的只有一个：让 ``π·|v̄| > |w̄|``
    的 **v 主导格**真实出现（本模块度量下 ξ-槽响应天然被 1/r_a 放大的 φ-槽与
    单位通量 BC 压制，生产量级永远达不到 v 主导，实测见 task-5-report.md
    §Important-3）。``amp_xi=1`` 起的**生产量级**收敛性由
    ``test_inclined_production_scale_converges`` 覆盖。

    ξ-槽沿 ξ 取 ``sin(πs/s_max)``（两端为零 ⇒ 散度源项不含边界贡献，线性叠加
    自洽）；``tan(β)`` 为井斜引起的重力方位分量比例（R29：ξ-槽配 sin πφ·sin β）。
    """
    c = (0.5 + 0.2 * np.cos(np.pi * g["phi"]))[:, None] * np.ones_like(g["H"])
    chi = 1.0 + 2.0 * c - c ** 3
    b = np.zeros((2,) + g["H"].shape)
    b[0] = amp_phi * chi
    b[1] = (amp_xi * chi * np.sin(np.pi * g["phi"])[:, None]
            * np.tan(np.deg2rad(beta_deg))
            * np.sin(np.pi * g["s"] / g["s"][-1])[None, :])
    return c, b


def _converged_inclined_run(amp_xi=1.0e6, beta_deg=45.0, tol=1e-10, max_outer=300):
    g = _geom(nz=9)
    c, b = _inclined_b_field(g, beta_deg=beta_deg, amp_phi=1.0, amp_xi=amp_xi)
    hb = _hb_closure()
    hb.set_pressure_gradient(300.0)
    psi = solve_stream_function_nonlinear(g, c, hb, b, Gb=(0.0, 0.0), omega=0.5,
                                          tol=tol, max_outer=max_outer)
    return g, c, b, hb, psi


def test_inclined_strong_buoyancy_converges_with_v_dominant_cells():
    """斜井（b_ξ≠0）强浮力应力算例收敛，且存在 v 分量为主的格点（π 加权被真实触发）。

    生产含义：本仓井斜 1.8–15° 内方位浮力经 (4.22) 的 b 向量承担、闭包反演口径
    Gb=0（R-T5-4 ①的显式契约）；该算例证明反演口径在 v-主导格仍自洽收敛。
    """
    g, c, b, hb, psi = _converged_inclined_run()
    assert np.all(np.isfinite(psi))
    w, v = velocity_from_stream_function(psi, g)
    u_w = np.hypot(w, np.pi * v)
    dominant = (np.pi * np.abs(v) > np.abs(w)) & (u_w > 0.01 * float(u_w.max()))
    assert int(dominant.sum()) > 0, "算例未产生 v 主导格点，π 加权未被触发"


def test_inclined_case_closure_relation_uses_weighted_u_caliber():
    """闭合关系必须用**加权口径** ``ū = hypot(w̄, π·v̄)``（R-T5-5）：

    收敛场上 v-主导格点处，加权口径闭合残差 ≤ 1e-6，而未加权口径偏差达 O(0.1)
    （纯 v 格理论上限 (π−1)/π ≈ 0.68，实测场混合格约 0.16）——两口径相差
    ~9 个量级，钉住"两分量差 π 倍被 hypot 混合"这一错误。
    """
    g, c, b, hb, psi = _converged_inclined_run()
    H = g["H"]
    w, v = velocity_from_stream_function(psi, g)
    G_inj = np.asarray(hb._G, dtype=float)
    I1 = np.asarray(hb.mobility(c, M_RATIO, ETA1, ETA2, H), dtype=float)
    pred = (I1 / H ** 2) * (H * G_inj)                 # Ī₁·G̃（量纲链）
    u_w = np.hypot(w, np.pi * v)
    u_p = np.hypot(w, v)                               # 未加权（错误口径）
    keep = (G_inj > 0) & (u_w > 0.05 * float(u_w.max()))
    rel_w = float(np.max(np.abs(u_w - pred)[keep] / np.abs(pred[keep])))
    rel_p = float(np.max(np.abs(u_p - pred)[keep] / np.abs(pred[keep])))
    assert rel_w < 1e-6, f"加权口径闭合残差 {rel_w:.3e} 超容差"
    assert rel_p > 0.05, f"未加权口径残差 {rel_p:.3e} 过小，判别力不足"


def test_inclined_production_scale_converges():
    """生产量级斜井算例（β=15°、amp_xi=1 ≈ 现场浮力量级）收敛且有限。

    与应力算例的分工：本条证明**生产参数**下外迭代行为正常（v 分量不主导，
    π 加权属正确性硬化）；应力算例负责触发 v 主导格。
    """
    g, c, b, hb, psi = _converged_inclined_run(amp_xi=1.0, beta_deg=15.0)
    assert np.all(np.isfinite(psi))
    w, v = velocity_from_stream_function(psi, g)
    assert float(np.max(np.abs(v) / np.maximum(np.abs(w), 1e-300))) < 1.0
