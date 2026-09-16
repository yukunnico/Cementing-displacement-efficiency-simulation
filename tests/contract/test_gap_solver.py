"""A-1a：一维间隙弱解（定 G 模式）。先钉解析极限锚，再实现。

锚点（全部有闭式解，可逐位核对）：

① 单流体牛顿 ⇒ 槽流 ``u(ỹ)=G(1−ỹ²)/(2η)``、``Ī₁=1/(3η)``；
② 两层牛顿 ⇒ ``I₁`` = Z&F22 (4.21a)（= ``two_layer.mobility_i1`` 输出口径）；
③ Bingham 单流体 ⇒ 塞流区 ``du/dỹ≡0``，屈服带与 Bingham 槽流解逐位一致；
④ 护栏（台账要求，brief 之外补）⇒ 两层牛顿恒等式在 **H≠1** 下仍成立，锁死
   ``I₁`` 的 H 量纲口径（brief 只测 H=1，漏乘/错乘 H 因子不会变红）。

⚠️ ③ 的塞流位置：brief 原文写 ``y_plug = 1 − τY/|G|``（即断言 ``y>0.7`` 处
``du/dỹ≡0``）——这与同一 brief 的 ①（``u=(1−ỹ²)G/(2η)`` ⇒ 对称面在 ỹ=0、
壁面在 ỹ=1）以及 B&F25 (A1)/(A2)（``τ̃₂(0)=0``、``ũ(1)=0``）矛盾：零剪应力面在
ỹ=0 ⇒ 塞流必紧邻 ỹ=0（本例 ``y_plug = τY/|G| = 0.3``），而 ``y∈(0.7,1]`` 恰是
壁面屈服带（``|du/dỹ| = (Gỹ−τY)/κ ∈ [4,7]``，见 ``task-1-report.md``）。本测试
按 (A1)/(A2) 与 ① 修正为 ``y_plug = τY/|G|``，并把 brief 标题里的
"matches_slot_solution" 落成逐位断言。
"""

import warnings

import numpy as np
import pytest

from cemdisp.models2d.gap_solver import solve_fixed_G


def test_newtonian_single_fluid_matches_slot_solution():
    """单流体牛顿（c̄=0）⇒ Ī₁=1/(3η)，槽流 u(y)=G(1−y²)/(2η)。"""
    sol = solve_fixed_G(c_bar=0.0, n=(1.0, 1.0), kappa=(2.0, 2.0),
                        tau_y=(0.0, 0.0), G=(1.0, 0.0))
    assert sol.converged
    assert sol.I1 == pytest.approx(1.0 / (3.0 * 2.0), rel=1e-8)
    y = sol.y
    assert np.allclose(sol.u, (1.0 - y**2) / (2.0 * 2.0), rtol=1e-6, atol=1e-8)


def test_newtonian_two_layer_matches_two_layer_closure():
    """牛顿两层 ⇒ Ī₁ 与 two_layer.mobility_i1 的归一化闭式一致。"""
    from cemdisp.models2d.two_layer import mobility_i1
    c, e1, e2 = 0.4, 3.0, 1.0
    m = e1 / e2
    sol = solve_fixed_G(c_bar=c, n=(1.0, 1.0), kappa=(e1, e2),
                        tau_y=(0.0, 0.0), G=(1.0, 0.0))
    expect = float(np.asarray(mobility_i1(c, m, eta1=e1, eta2=e2, H=1.0)))
    assert sol.I1 == pytest.approx(expect, rel=1e-7)


def test_bingham_plug_exists_and_matches_slot_solution():
    """Bingham 单流体：|τ|≤τY 区 du/dy=0（塞流），塞流起于 y=τY/|G|。

    塞流位置按 B&F25 (A1)（零剪应力面在 ỹ=0）修正，详见模块 docstring。
    屈服带内逐位对照 Bingham 槽流解：
    ``u(ỹ) = ∫_ỹ^1 (Gỹ′−τY)/κ dỹ′ = (G/2κ)(1−ỹ²) − (τY/κ)(1−ỹ)``。
    """
    G, tauY, kappa = 10.0, 3.0, 1.0
    sol = solve_fixed_G(c_bar=0.0, n=(1.0, 1.0), kappa=(kappa, kappa),
                        tau_y=(tauY, tauY), G=(G, 0.0))
    assert sol.converged
    y_plug = tauY / G
    du = np.gradient(sol.u, sol.y)
    mask_plug = sol.y < y_plug - 1e-6
    assert np.max(np.abs(du[mask_plug])) < 1e-6
    # 屈服带（壁面侧）：与 Bingham 槽流解逐位一致
    mask_yield = sol.y > y_plug + 1e-6
    u_exact = (G / (2.0 * kappa)) * (1.0 - sol.y**2) - (tauY / kappa) * (1.0 - sol.y)
    assert np.allclose(sol.u[mask_yield], u_exact[mask_yield], rtol=1e-6, atol=1e-8)


def test_newtonian_two_layer_identity_holds_at_H_not_one():
    """两层牛顿在任意 H 下等于 mobility_i1（Z&F22 (4.21a)）——锁死 H 口径。

    H≠1 时若实现漏乘/错乘 H 的量纲因子，H=1 的测试不会变红，这条会。
    c̄ 取 0.2/0.5/0.8 三档；m=1.0 时闭式 n₁≡1/3（与 c̄ 无关），
    I₁ = H³/(3η)·(Hū/G) 的关系式在此恒成立。
    """
    from cemdisp.models2d.two_layer import mobility_i1
    H = 7.0
    for c in (0.2, 0.5, 0.8):
        sol = solve_fixed_G(c_bar=c, n=(1.0, 1.0), kappa=(2.0, 2.0),
                            tau_y=(0.0, 0.0), G=(1.0, 0.0), H=H)
        expect = float(np.asarray(mobility_i1(c, 1.0, eta1=2.0, eta2=2.0, H=H)))
        assert sol.I1 == pytest.approx(expect, rel=1e-10)


def test_newtonian_two_layer_i2_and_q0_match_two_layer_closed_forms():
    """两层牛顿的 ``I₂``/``q₀`` 与 Z&F22/``two_layer`` 闭式一致——锁 (2.15)/(2.22) 权重。

    (2.15) 第二项被积函数是 ``ỹ(H−ỹ)/η₁``（含 ỹ 一次因子）；(2.22) 的 q₀ 须退回
    Z&F22 (4.25) 的牛顿闭式 ``isotropic_flux_q0``。brief 的三条锚不含这两个量，
    但它们是 Task 4 查表的直接输入（本条测试在实现期抓到过 I₂ 权重漏 ỹ 的缺陷）。
    """
    from cemdisp.models2d.two_layer import isotropic_flux_q0, mobility_i2
    e1, e2, m = 3.0, 1.0, 3.0
    for H in (1.0, 3.0):
        for c in (0.0, 0.1, 0.4, 0.7, 0.9, 1.0):
            sol = solve_fixed_G(c_bar=c, n=(1.0, 1.0), kappa=(e1, e2),
                                tau_y=(0.0, 0.0), G=(1.0, 0.0), H=H)
            expect_i2 = float(mobility_i2(c, m, eta1=e1, eta2=e2, H=H))
            assert sol.I2 == pytest.approx(expect_i2, rel=1e-10, abs=1e-14)
            assert sol.q0 == pytest.approx(float(isotropic_flux_q0(c, m)),
                                           rel=1e-10, abs=1e-14)


def test_failure_paths_are_loud():
    """不收敛抛 ``RuntimeError``；非共线 G/Gb、全场未屈服抛 ``ValueError``。

    设计规格 §3 A-1：“超时/不收敛**抛错**（不静默回退）”。
    """
    base = dict(c_bar=0.0, n=(0.5, 0.5), kappa=(1.0, 1.0), tau_y=(0.0, 0.0))
    with pytest.raises(RuntimeError, match="未收敛"):
        solve_fixed_G(G=(1.0, 0.0), max_iter=3, **base)
    with pytest.raises(ValueError, match="不共线"):
        solve_fixed_G(G=(1.0, 0.0), Gb=(0.0, 1.0), **base)
    with pytest.raises(ValueError, match="未屈服"):
        solve_fixed_G(c_bar=0.0, n=(1.0, 1.0), kappa=(1.0, 1.0),
                      tau_y=(100.0, 100.0), G=(1.0, 0.0))


def test_divergence_guard_raises_runtimeerror_and_leaks_no_warning():
    """``r≲0.5`` 数值发散：必须**先**抛 ``RuntimeError``，且不泄漏 numpy 运行时告警。

    发散机制：λ̃ 模态的收缩因子 ~``|1−1/r|``（r=0.2 ⇒ ~4 倍/步）。
    未加护栏时实测先是 `RuntimeWarning: overflow encountered in reduce` 逃逸
    ——在 `-W error::RuntimeWarning` 下会把本意的 RuntimeError 变成 RuntimeWarning，
    破坏"响亮但干净"的契约（下游 Task 2/3 会看到错误的异常类型）。
    """
    args = dict(c_bar=0.0, n=(1.0, 1.0), kappa=(1.0, 1.0), tau_y=(0.0, 0.0),
                G=(1.0, 0.0), r=0.2, max_iter=2000)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        with pytest.raises(RuntimeError, match="发散"):
            solve_fixed_G(**args)