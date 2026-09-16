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

A-1b（Task 2）追加：⑤ ``τ_Y>0`` 的闭包**值**断言（解析求积闭包是 Phase A 的全部
要点，之前无任何测试钉住它）；⑥ H 标度 ``I₁ ∝ H^{2+1/n}``；⑦ Uzawa 收敛性
（步数有限 + 幂律指数越大越快，B&F25 图 18 的定性趋势）；⑧ ``Gb≠0`` 的判别性
交叉检查（三条锚全部 ``Gb=0`` ⇒ 实现里的界面项与两处 Gb 拆分恒为零，有 bug 也看不见）。
"""

import math
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
    # 闭包值也钉住（τ_Y>0 ⇒ 解析求积闭包，正是 Phase A 的要点）：
    #   I₁ = ∫_0^1 ỹ²·(γ/|τ̃|) d = ∫_{τY/G}^1 (ỹ² − (τY/G)·ỹ) dỹ
    #      = [ỹ³/3 − (τY/G)ỹ²/2]_{τY/G}^{1}；G=10、τY=3 ⇒ 1/3 − 0.15 − 0.009 + 0.0135
    #      = 0.1878333…
    assert sol.I1 == pytest.approx(0.1878333333333333, rel=1e-10)


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


# --------------------------------------------------------------------------- #
# A-1b（Task 2）：H 标度、τ_Y>0 闭包值、Gb≠0 判别性、Uzawa 收敛性
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("n_exp,expect", [(1.0, 3.0), (0.7, 2 + 1/0.7), (0.4, 2 + 1/0.4)])
def test_single_fluid_power_law_I1_scales_as_H_pow_2_plus_1_over_n(n_exp, expect):
    """阶段 B 已数值核验：I₁ ∝ H^{2+1/n}（**固定量纲 G** 扫 H，直接看 ``sol.I1``）。

    单流体幂律用 (2.14) + 幂律流变有闭式 ``Ĩ₁ = κ̃^{−1/n}G̃^{1/n−1}/(2+1/n)``；
    代入 (A4)（``κ̃=κ/Hⁿ``、``G̃=H·G``）⇒ ``I₁ = H^{2+1/n}·κ^{−1/n}G^{1/n−1}/(2+1/n)``。

    ⚠️ brief 原文写 ``solve_fixed_G(..., G=(1.0/H, 0.0)).I1 * H**3``：它**同时**改了
    量纲 G（⇒ ``G̃=H·G≡1`` ⇒ 只测出 ``I₁ = H³/(2+1/n)``，斜率 ``4−1/n``）。实测
    n=1 偶然为 3.0（假绿）、n=0.7 为 2.5714（应 3.4286）、n=0.4 因 ``G̃/κ̃`` 过大
    在 ``max_iter=2000`` 内不收敛（RuntimeError）⇒ 三种失败。故按裁定改为本条。
    """
    Hs = np.array([0.006, 0.008, 0.010, 0.012, 0.014])
    sols = [solve_fixed_G(c_bar=0.0, n=(n_exp, n_exp), kappa=(1.0, 1.0),
                          tau_y=(0.0, 0.0), G=(1.0, 0.0), H=float(H)) for H in Hs]
    assert all(sol.converged for sol in sols)
    I1s = np.array([sol.I1 for sol in sols])
    slope = np.polyfit(np.log(Hs), np.log(I1s), 1)[0]
    assert slope == pytest.approx(expect, abs=1e-3)
    # 加固（brief 之外）：纯幂律斜率对"与 H 无关的常数错误"不敏感 ⇒ 再钉解析前因子
    # （κ=1、G=1 ⇒ I₁ = H^{2+1/n}/(2+1/n)；实测偏差 ≤2.3e-9，为分段 GL-16 求积误差）
    assert I1s == pytest.approx(Hs**(2.0 + 1.0 / n_exp) / (2.0 + 1.0 / n_exp), rel=1e-6)


def test_uzawa_converges_and_reports_iters():
    sol = solve_fixed_G(0.5, (0.7, 0.7), (1.0, 1.0), (0.0, 0.0), G=(5.0, 0.0))
    assert sol.converged and 0 < sol.iters <= 2000


def test_uzawa_iters_decrease_with_power_law_index():
    """幂律指数越大收敛越快（B&F25 §A.2.3 图 18 的定性趋势）。

    参数按图 18 的题注：「two identical fluids … isodense flows (b = 0) at
    ``κ₁ = κ₂ = 1``, ``c̄ = 0.5``」——两种相同流体（等 κ、等 n、τ_Y=0 的
    shear-thinning 族）、等密度（Gb=0）、H=1。实测步数 n=0.2→1.0：
    ``96 / 58 / 44 / 37 / 33``（牛顿最快，与论文「increasing the power-law index
    (n) improves the convergence as does reducing the yield stress」一致）。

    ⚠️ **定 G 口径的依赖**：同一组参数下该趋势随 ``G`` 变化（单位参数 G=1 时
    单调下降如上，G=0.5 时反转为 ``14/25/30/32/32``，G=2 时为 ``1440/150/66/43/34``）。
    图 18 的原始口径是**定均速**（A.2.2，``solve_fixed_mean_velocity``，Task 3），
    在那里比的是同一 ``w̄`` 下的 n 族。故本条只钉「G=1 单位参数下 n 越大越快」。
    """
    iters = {}
    for n_exp in (0.2, 0.4, 0.6, 0.8, 1.0):
        sol = solve_fixed_G(c_bar=0.5, n=(n_exp, n_exp), kappa=(1.0, 1.0),
                            tau_y=(0.0, 0.0), G=(1.0, 0.0), H=1.0)
        iters[n_exp] = sol.iters
    assert iters[1.0] <= iters[0.2]


def test_buoyancy_term_enters_closures_when_gb_nonzero():
    """``Gb`` 真的进入闭包（判别性交叉检查）。

    三条解析锚全部 ``Gb=0``——此时实现里 ``a₂≡a₁``、界面项 ``(a₂−a₁)c̄`` 与两处
    Gb 拆分恒等于零，Gb 路径有 bug 也看不见。本条取 ``τ_Y>0``（牛顿且 ``τ_Y=0``
    时 ``1/η̃≡1/κ`` 与应力场无关、对 Gb 完全盲），并要求 ``Gb`` 沿 ``G`` 共线。
    """
    args = dict(c_bar=0.5, n=(1.0, 1.0), kappa=(1.0, 1.0), tau_y=(0.1, 0.1),
                G=(1.0, 0.0), H=1.0)
    dense = solve_fixed_G(Gb=(0.0, 0.0), **args)
    buoy = solve_fixed_G(Gb=(0.3, 0.0), **args)
    assert buoy.converged
    # 判别性：I₁ 与等密度情形不同（且非浮点噪声量级 ⇒ 实测差 ~1.5%）
    assert not math.isclose(dense.I1, buoy.I1)
    assert abs(buoy.I1 / dense.I1 - 1.0) > 1e-3
    # 结构性：闭包仍为正、无滑移仍成立（Gb≠0 无独立文献闭式 ⇒ 不钉数字）
    assert dense.I1 > 0.0 and dense.I2 > 0.0
    assert buoy.I1 > 0.0 and buoy.I2 > 0.0
    assert buoy.u[-1] == 0.0
