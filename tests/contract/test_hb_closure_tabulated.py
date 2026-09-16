"""A-2：HB 无量纲群 (2.29)—(2.34) 与 ``HBClosure``（B&F25 闭包提供者）。

文件名沿用计划（``test_hb_closure_tabulated.py``）：A-2 原设计为"预计算
``(c̄, B)`` 二维表 + 插值"。控制器的任务 4 裁定**覆盖**了该设计——闭包值必须由
``gap_solver`` 用**注入的 G** 与真实应力场现算（``(c̄, B)`` 表只在固定 G 下自洽，
B 里的 ``γ̇₀`` 依赖局部应力尺度）；查表/缓存只作为**加速**手段保留。故本文件的
断言全部针对"真实 G 求值"口径，实测耗时与缓存命中率见 ``task-4-report.md``。

锚点分层（避免"用产出该值的同一函数再算一遍"的恒真断言）：

① 牛顿退化：``hb_groups`` 在 ``n=1, τ_Y=0`` 下 ``B=0``、``m=κ₁/κ₂``（brief 原条），
   且 ``κ_k``/``τ_{Y,k}`` 命中 (2.32)/(2.30) 的闭式；
② 群换算的**独立复算**（物理定义式：表观黏度 ``η_k=κ_kγ̇₀^{n_k−1}`` 与
   ``max τ_{Y,k}=B/(1+B)``），不复用被测函数的中间量；
③ 牛顿极限**逐位**（R2 结构性短路：monkeypatch 钉死"不查表、不进 gap_solver"）；
④ 非牛顿解析锚：Bingham 单流体 ``Q/Q_N = 1−1.5B+0.5B³``（B&F25 §5.1 同款）、
   纯幂律 ``I₁ = H^{2+1/n}G^{1/n−1}κ^{−1/n}/(2+1/n)``；
⑤ 与 ``gap_solver`` 的交叉核对：定 G 模式钉**接线口径**（逐位）、**定均速模式**
   钉数值正确性（不同算法、同一不动点，``ū = I₁G/H`` 的唯一自洽对）；
⑥ R-T1-6：全场未屈服 ⇒ 地板值 + 单次告警，``ValueError`` 不逃逸。

References
----------
Bararpour & Frigaard (2025), *JFM* **1022**, A15：(2.13)—(2.15)、(2.29)—(2.34)、
附录 A.2.1/A.2.2（``gap_solver`` 消费）。
"""

import math
import warnings

import numpy as np
import pytest

from cemdisp.models2d.gap_solver import solve_fixed_G, solve_fixed_mean_velocity
from cemdisp.models2d.hb_closure import (
    HBClosure,
    NewtonianClosure,
    _STATIC_WALL_MOBILITY_FLOOR,
)
from cemdisp.models2d.two_layer import hb_groups, mobility_i1

# 两层 HB 参数（量纲自洽：κ~O(1) Pa·sⁿ、τ_Y~O(1) Pa、H~8 mm、G~10² Pa/m）
HB_N = (0.7, 0.8)
HB_KAPPA = (1.4, 0.9)
HB_TAU_Y = (2.0, 0.5)
HB_H = 0.008
HB_G = 500.0


# --------------------------------------------------------------------------- #
# ① hb_groups 牛顿退化（brief 原条 + 闭式钉住）
# --------------------------------------------------------------------------- #
def test_hb_groups_reduce_to_newtonian():
    """``n₁=n₂=1, τ_Y≡0`` ⇒ ``B=0``、``m=κ₁/κ₂``（(2.29)—(2.34) 的牛顿退化）。"""
    g = hb_groups(kappa1=0.1, kappa2=0.05, n1=1.0, n2=1.0,
                  tauY1=0.0, tauY2=0.0, gamma0=100.0)
    assert g["B"] == pytest.approx(0.0)
    assert g["m"] == pytest.approx(0.1 / 0.05)
    # 牛顿下 μ̂e = √(κ₁κ₂)γ̇₀（(2.29) 的 n−1=0 ⇒ 与 γ̇₀ 无关），τ̂₀ = μ̂eγ̇₀ + 0
    mu_e = math.sqrt(0.1 * 0.05)
    assert g["mu_e"] == pytest.approx(mu_e)
    assert g["tau_0"] == pytest.approx(mu_e * 100.0)
    # (2.32) 右式在 B=0、n=1 下：κ_k = m^{∓1/2}
    assert np.asarray(g["kappa_k"]) == pytest.approx(
        (math.sqrt(2.0), 1.0 / math.sqrt(2.0))
    )
    assert np.asarray(g["tau_Yk"]) == pytest.approx((0.0, 0.0))


def test_hb_groups_match_physical_definitions():
    """非牛顿：群换算命中物理定义式（独立复算，不复用被测函数中间量）。

    - ``μ̂e = √(η₁η₂)``、``η_k = κ_kγ̇₀^{n_k−1}``（(2.29)：μ̂e 是 γ̇₀ 处的有效黏度尺度）；
    - ``m = η₁/η₂``（(2.32) 的物理内容 = γ̇₀ 处表观黏度比）；
    - ``κ_k = m^{∓1/2}/(1+B)``（(2.32) 右式，与左式 ``κ_kγ̇₀^{n_k}/τ̂₀`` 代数等价）；
    - ``max τ_{Y,k} = B/(1+B) < 1``（文献用 ``max`` 的目的：缩放屈服应力落在 [0,1)）。
    """
    k1, k2, n1, n2, ty1, ty2, gam = 1.4, 0.9, 0.7, 0.8, 2.0, 0.5, 10.0
    g = hb_groups(kappa1=k1, kappa2=k2, n1=n1, n2=n2, tauY1=ty1, tauY2=ty2, gamma0=gam)

    eta1_app = k1 * gam ** (n1 - 1.0)          # 表观黏度 η₁(γ̇)
    eta2_app = k2 * gam ** (n2 - 1.0)          # 表观黏度 η₂(γ₀)
    mu_e_exp = math.sqrt(eta1_app * eta2_app)  # (2.29)
    tau_0_exp = mu_e_exp * gam + max(ty1, ty2)  # (2.30)
    b_exp = max(ty1, ty2) / (mu_e_exp * gam)   # (2.33)
    m_exp = eta1_app / eta2_app                # (2.32) 物理内容

    assert g["mu_e"] == pytest.approx(mu_e_exp, rel=1e-14)
    assert g["tau_0"] == pytest.approx(tau_0_exp, rel=1e-14)
    assert g["B"] == pytest.approx(b_exp, rel=1e-14)
    assert g["m"] == pytest.approx(m_exp, rel=1e-14)
    assert np.asarray(g["kappa_k"]) == pytest.approx(
        (math.sqrt(m_exp) / (1.0 + b_exp), 1.0 / (math.sqrt(m_exp) * (1.0 + b_exp))),
        rel=1e-14,
    )
    assert np.asarray(g["tau_Yk"]) == pytest.approx(
        (ty1 / tau_0_exp, ty2 / tau_0_exp), rel=1e-14
    )
    assert max(g["tau_Yk"]) == pytest.approx(b_exp / (1.0 + b_exp), rel=1e-14)
    assert max(g["tau_Yk"]) < 1.0


def test_hb_groups_validates_inputs():
    """``γ̇₀ ≤ 0``/``κ ≤ 0``/``n ≤ 0``/``τ_Y < 0`` 一律构造期拒绝。"""
    base = dict(kappa1=1.0, kappa2=1.0, n1=0.8, n2=0.8, tauY1=0.1, tauY2=0.1, gamma0=10.0)
    hb_groups(**base)                                   # 合法基准不抛错
    for bad in ({"gamma0": 0.0}, {"gamma0": -1.0}, {"kappa1": 0.0},
                {"kappa2": -1.0}, {"n1": 0.0}, {"tauY1": -1e-9}):
        with pytest.raises(ValueError):
            hb_groups(**{**base, **bad})


# --------------------------------------------------------------------------- #
# ③ HBClosure：牛顿极限逐位 + R2 结构性短路
# --------------------------------------------------------------------------- #
def test_hb_closure_newtonian_limit():
    """``n=1, τ_Y≡0`` ⇒ 逐位等于 ``NewtonianClosure``（收紧自 brief 的 rtol=1e-6）。

    R2：**结构性短路**——直接委托 ``NewtonianClosure``，不查表、不进 ``gap_solver``
    （Task 5 与 L1 硬约束要求牛顿极限逐位 = HEAD，任何数值路径都做不到逐位）。
    """
    c = np.linspace(0.0, 1.0, 7)
    H = np.full(7, 0.012)
    hb = HBClosure(n=(1.0, 1.0), kappa=(1.0, 1.0), tau_y=(0.0, 0.0), m=1.0, B=0.0)
    nw = NewtonianClosure()
    assert hb.is_newtonian_limit
    assert np.array_equal(hb.mobility(c, 1.0, 1.0, 1.0, H), nw.mobility(c, 1.0, 1.0, 1.0, H))
    assert np.array_equal(hb.buoyant_mobility(c, 1.0, 1.0, 1.0, H),
                          nw.buoyant_mobility(c, 1.0, 1.0, 1.0, H))
    # 口径仍与 two_layer 逐位（HBClosure 不得引入任何算术）
    assert np.array_equal(hb.mobility(c, 1.0, 1.0, 1.0, H),
                          np.asarray(mobility_i1(c, 1.0, eta1=1.0, eta2=1.0, H=H)))


def test_hb_closure_newtonian_limit_never_calls_gap_solver(monkeypatch):
    """R2 的结构性证明：短路路径**不调用** ``gap_solver``（未注入 G 也照常工作）。

    牛顿极限的逐位要求不可能由"进 gap_solver 再算回来"满足（数值解与解析闭式
    差 ~1e-9 量级）；本条用 monkeypatch 把 ``gap_solver`` 的批量闭包入口打哑
    （Task 4.5 起非牛顿路径走 ``closure_integrals_batch``），短路若被改回数值路径即红。
    """
    import cemdisp.models2d.hb_closure as hbc

    def _boom(*args, **kwargs):  # pragma: no cover - 触发即失败
        raise AssertionError("牛顿极限不应进入 gap_solver 数值路径")

    hb = HBClosure(n=(1.0, 1.0), kappa=(0.058, 0.171), tau_y=(0.0, 0.0),
                   m=0.058 / 0.171, B=0.0)
    c = np.linspace(0.0, 1.0, 5)
    H = np.full(5, 0.01)
    expect = NewtonianClosure().mobility(c, 0.34, 0.058, 0.171, H)   # 打哑前先取期望

    monkeypatch.setattr(hbc, "closure_integrals_batch", _boom)

    # 未注入 G：短路路径不需要（也不会去找）压力梯度
    assert np.array_equal(hb.mobility(c, 0.34, 0.058, 0.171, H), expect)

    # 对照：非牛顿（τ_Y>0）必须走 gap_solver —— 打哑后应触发 AssertionError
    hb_hb = HBClosure(n=(1.0, 1.0), kappa=(2.0, 2.0), tau_y=(0.5, 0.5), m=1.0, B=0.2)
    hb_hb.set_pressure_gradient(1.0)
    assert not hb_hb.is_newtonian_limit
    with pytest.raises(AssertionError, match="gap_solver"):
        hb_hb.mobility(c, 1.0, 2.0, 2.0, H)


def test_hb_closure_requires_injected_gradient():
    """非牛顿：未注入 G ⇒ 显式 ``RuntimeError``（不静默回退牛顿闭包）。"""
    hb = HBClosure(n=HB_N, kappa=HB_KAPPA, tau_y=HB_TAU_Y, m=1.0, B=0.3)
    c = np.full((3, 2), 0.5)
    H = np.full((3, 2), HB_H)
    with pytest.raises(RuntimeError, match="set_pressure_gradient"):
        hb.mobility(c, 1.0, 1.4, 0.9, H)
    with pytest.raises(RuntimeError, match="set_pressure_gradient"):
        hb.buoyant_mobility(c, 1.0, 1.4, 0.9, H)
    hb.set_pressure_gradient(HB_G)
    assert np.all(np.asarray(hb.mobility(c, 1.0, 1.4, 0.9, H)) > 0.0)


def test_hb_closure_implements_protocol():
    """``HBClosure`` 满足阶段 B 冻结的 ``ClosureProvider`` 协议。"""
    from cemdisp.models2d.hb_closure import ClosureProvider

    assert isinstance(HBClosure(n=HB_N, kappa=HB_KAPPA, tau_y=HB_TAU_Y, m=1.0, B=0.3),
                      ClosureProvider)


# --------------------------------------------------------------------------- #
# ④ 非牛顿解析锚（独立于实现）
# --------------------------------------------------------------------------- #
def test_hb_closure_bingham_slot_analytic_anchor():
    """单流体 Bingham（``c̄=0``）⇒ ``I₁ = I₁_N·(1 − 1.5B + 0.5B³)``（槽流 Q/Q_N）。

    ``B = τ_Y/(G·H)`` 为壁面 Bingham 数（tilde 系壁面应力 ``G̃=G·H``）；该因子是
    Bingham 槽流的标准结果（与设计规格 §5.1 引用的 ``Q/Q_N`` 同式），**不调用**被测
    函数。``H=1`` 下取 ``I₁_N = 1/(3κ)``。
    """
    kappa, tau_y, G = 2.0, 0.4, 1.0
    hb = HBClosure(n=(1.0, 1.0), kappa=(kappa, kappa), tau_y=(tau_y, tau_y),
                   m=1.0, B=tau_y / (kappa * G))
    hb.set_pressure_gradient(G)
    I1 = float(hb.mobility(0.0, 1.0, kappa, kappa, 1.0))
    B_wall = tau_y / (G * 1.0)
    expect = (1.0 / (3.0 * kappa)) * (1.0 - 1.5 * B_wall + 0.5 * B_wall ** 3)
    assert I1 == pytest.approx(expect, rel=1e-9)


def test_hb_closure_power_law_H_scaling_and_prefactor():
    """单流体幂律（``τ_Y=0``）：``I₁ = H^{2+1/n}G^{1/n−1}κ^{−1/n}/(2+1/n)``。

    ``Ĩ₁ = ∫ỹ²/η̃ dỹ``、``1/η̃ = γ̇/|τ̃| = κ̃^{−1/n}G̃^{1/n−1}`` 代入 (A4)
    （``κ̃=κ/Hⁿ``、``G̃=H·G``）⇒ 上式（与 Task 2 的 gap_solver 同款解析锚，
    此处钉的是 HBClosure 对 H/κ 的**传递口径**）。
    """
    n_exp, kappa = 0.6, 2.0
    hb = HBClosure(n=(n_exp, n_exp), kappa=(kappa, kappa), tau_y=(0.0, 0.0),
                   m=1.0, B=0.0)
    hb.set_pressure_gradient(1.0)
    H1, H2 = 0.01, 0.02
    I1_a = float(hb.mobility(0.5, 1.0, kappa, kappa, H1))
    I1_b = float(hb.mobility(0.5, 1.0, kappa, kappa, H2))
    for H, got in ((H1, I1_a), (H2, I1_b)):
        expect = H ** (2.0 + 1.0 / n_exp) * kappa ** (-1.0 / n_exp) / (2.0 + 1.0 / n_exp)
        assert got == pytest.approx(expect, rel=1e-7)
    assert I1_b / I1_a == pytest.approx((H2 / H1) ** (2.0 + 1.0 / n_exp), rel=1e-12)


# --------------------------------------------------------------------------- #
# ⑤ 与 gap_solver 的交叉核对（接线口径 + 定均速模式数值核对）
# --------------------------------------------------------------------------- #
def test_hb_closure_matches_gap_solver_fixed_G_and_mean_velocity():
    """定 G 逐位（**接线口径**）+ 定均速模式（**数值正确性**，不同算法同一不动点）。

    - 定 G 逐位：``HBClosure`` 必须把注入的 ``G`` 与调用方 ``H`` 原样送进
      ``solve_fixed_G``（钉 G/H 的传递，不引入任何算术）。
    - 定均速核对：``ū = I₁·G/H``（(2.13)，``Hū = Ī₁G``）是闭包唯一的自洽对；
      ``solve_fixed_mean_velocity(ū)`` 经**另一条迭代路径**反求 G 并给出 I₁，
      两条路径必须在 Uzawa 容差内一致。这是本文件里唯一"不重跑同一函数"的数值核对。
    """
    args = dict(c_bar=0.45, n=HB_N, kappa=HB_KAPPA, tau_y=HB_TAU_Y)
    H = 1.0
    ref = solve_fixed_G(G=(3.755, 0.0), H=H, ny=201, **args)

    hb = HBClosure(n=HB_N, kappa=HB_KAPPA, tau_y=HB_TAU_Y, m=1.0, B=0.3)
    hb.set_pressure_gradient(ref.G)
    got = float(hb.mobility(0.45, 1.0, 1.4, 0.9, H))
    assert got == ref.I1                       # 接线口径：同参同值（逐位）
    assert float(hb.buoyant_mobility(0.45, 1.0, 1.4, 0.9, H)) == ref.I2

    mean = solve_fixed_mean_velocity(u_bar=ref.u_bar, H=H, ny=201, **args)
    assert mean.G[0] == pytest.approx(ref.G[0], rel=1e-8)   # G 往返
    assert got == pytest.approx(mean.I1, rel=1e-8)          # 跨模式 I₁ 一致


def test_hb_closure_gb_is_consumed_and_non_collinear_supported():
    """(2.13) 的 ``Gb`` 不是报告量：注入后进真实求解；非共线 ``G``/``Gb`` 亦可解。

    判别性（沿用 Task 2/3 的 Gb 检查风格）：``τ_Y>0`` 时 ``Gb`` 改变 I₁（牛顿且
    ``τ_Y=0`` 时 ``1/η̃≡1/κ`` 对 Gb 盲，测不出）。
    """
    c, H = 0.5, 1.0
    hb = HBClosure(n=(1.0, 1.0), kappa=(1.0, 1.0), tau_y=(0.1, 0.1), m=1.0, B=0.05)

    hb.set_pressure_gradient((1.0, 0.0))
    dense = float(hb.mobility(c, 1.0, 1.0, 1.0, H))
    hb.set_pressure_gradient((1.0, 0.0), Gb=(0.3, 0.0))
    buoy = float(hb.mobility(c, 1.0, 1.0, 1.0, H))
    assert not math.isclose(buoy, dense)
    assert abs(buoy / dense - 1.0) > 1e-3
    # 与直接调用一致（钉 Gb 也原样传递）
    ref = solve_fixed_G(c_bar=c, n=(1.0, 1.0), kappa=(1.0, 1.0), tau_y=(0.1, 0.1),
                        G=(1.0, 0.0), Gb=(0.3, 0.0), H=H, ny=201)
    assert buoy == ref.I1

    # 非共线（向量版，R-T1-3）：可解且与轴向/共线结果不同
    hb.set_pressure_gradient((500.0, 150.0), Gb=(100.0, 400.0))
    nc = np.asarray(hb.mobility(0.45, 1.0, 1.4, 0.9, HB_H))
    assert np.all(np.isfinite(nc)) and np.all(nc > 0.0)
    hb.set_pressure_gradient((500.0, 150.0))
    ax = float(hb.mobility(0.45, 1.0, 1.4, 0.9, HB_H))
    assert abs(float(nc.ravel()[0]) / ax - 1.0) > 1e-3


# --------------------------------------------------------------------------- #
# ⑥ R-T1-6：全场未屈服 ⇒ 地板 + 单次告警（ValueError 不逃逸）
# --------------------------------------------------------------------------- #
def test_hb_closure_static_wall_floor_and_warn_once():
    """全场未屈服（static wall layer）⇒ I₁ 地板值 + 单次 ``RuntimeWarning``。

    窄边强屈服格可能**合法地**全场未屈服（``gap_solver`` 抛 ``ValueError``）；该
    异常不得逃逸到求解器。地板口径复用 Phase B 先例
    ``stream_function._WALL_CONDUCTANCE_FLOOR`` 的同式思路：``I₁_eff =
    max(I₁_HB, 地板×I₁_牛顿)``（相对地板 ⇒ 与量纲/量级无关）；``I₂ ≡ 0``
    （``1/η̃≡0`` 处处 ⇒ ``Ĩ₂`` 的积分恒为 0，是精确值而非地板）。
    """
    c = np.linspace(0.0, 1.0, 5)
    H = np.full(5, 0.01)                      # G̃ = G·H = 1e−5 ≪ τ_Y = 5 ⇒ 处处未屈服
    # 地板**取值**也钉住（改动须显式）：与 Phase B 先例同量级（1e-6）
    assert _STATIC_WALL_MOBILITY_FLOOR == 1.0e-6
    hb = HBClosure(n=(1.0, 1.0), kappa=(2.0, 2.0), tau_y=(5.0, 5.0), m=1.0, B=0.9)
    hb.set_pressure_gradient(1.0e-3)
    with pytest.warns(RuntimeWarning, match="未屈服"):
        I1 = np.asarray(hb.mobility(c, 1.0, 2.0, 2.0, H))
    assert hb.n_static_wall_points == 5
    # 地板 = 地板系数 × 牛顿流动度（同一 c̄/m/η/H 口径）
    assert np.array_equal(
        I1, _STATIC_WALL_MOBILITY_FLOOR * np.asarray(mobility_i1(c, 1.0, eta1=2.0, eta2=2.0, H=H))
    )
    assert np.array_equal(
        np.asarray(hb.buoyant_mobility(c, 1.0, 2.0, 2.0, H)), np.zeros(5)
    )
    # 单次告警：同一实例再调不再告警（计数照常）
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        hb.mobility(c, 1.0, 2.0, 2.0, H)
        assert not [w for w in rec if issubclass(w.category, RuntimeWarning)]


def test_hb_closure_zero_gradient_point_is_floored_with_warning():
    """``G`` 与 ``Gb`` 同时为零的退化点同样地板化 + 告警（另有别于"未屈服"的成因）。"""
    hb = HBClosure(n=HB_N, kappa=HB_KAPPA, tau_y=HB_TAU_Y, m=1.0, B=0.3)
    hb.set_pressure_gradient(0.0)
    with pytest.warns(RuntimeWarning, match="不能同时为零"):
        I1 = np.asarray(hb.mobility(0.5, 1.0, 1.4, 0.9, HB_H))
    assert I1 == pytest.approx(
        _STATIC_WALL_MOBILITY_FLOOR * float(mobility_i1(0.5, 1.0, eta1=1.4, eta2=0.9, H=HB_H))
    )


def test_hb_closure_validates_inputs():
    """输入校验：``c̄ ∉ [0,1]``/``H ≤ 0``/``G`` 形状非法 ⇒ 显式 ``ValueError``（不地板化）。"""
    hb = HBClosure(n=HB_N, kappa=HB_KAPPA, tau_y=HB_TAU_Y, m=1.0, B=0.3)
    hb.set_pressure_gradient(HB_G)
    with pytest.raises(ValueError, match="c_bar"):
        hb.mobility(1.5, 1.0, 1.4, 0.9, HB_H)
    with pytest.raises(ValueError, match="H"):
        hb.mobility(0.5, 1.0, 1.4, 0.9, 0.0)
    for bad in ((1.0, 2.0, 3.0), np.zeros((3, 2, 5))):
        hb.set_pressure_gradient(bad)
        with pytest.raises(ValueError, match="G"):
            hb.mobility(np.full((3, 2), 0.5), 1.0, 1.4, 0.9, np.full((3, 2), HB_H))
    # 构造期校验
    for bad in ({"n": (0.0, 0.7)}, {"kappa": (1.0, 0.0)}, {"tau_y": (-1.0, 0.5)},
                {"m": 0.0}, {"B": -0.1}, {"ny_gap": 2}):
        with pytest.raises(ValueError):
            HBClosure(**{"n": HB_N, "kappa": HB_KAPPA, "tau_y": HB_TAU_Y, "m": 1.0,
                         "B": 0.3, **bad})


# --------------------------------------------------------------------------- #
# 形状契约与缓存
# --------------------------------------------------------------------------- #
def test_hb_closure_shape_contract():
    """返回与 ``c̄``/``H`` 广播后的数组；``(2,)+shape`` 的 G 场逐点非共线可用。"""
    ny, nz = 4, 3
    c = np.linspace(0.1, 0.9, ny)[:, None] * np.ones((1, nz))
    H = np.full((ny, nz), HB_H)
    hb = HBClosure(n=HB_N, kappa=HB_KAPPA, tau_y=HB_TAU_Y, m=1.0, B=0.3)
    # 轴向 (2,ny,nz) 向量场：第 2 分量恒零 ⇔ 标量口径
    axial = np.stack([np.full((ny, nz), HB_G), np.zeros((ny, nz))])
    hb.set_pressure_gradient(axial)
    I1 = hb.mobility(c, 1.0, 1.4, 0.9, H)
    assert isinstance(I1, np.ndarray) and I1.shape == (ny, nz)
    assert np.all(I1 > 0.0)
    # G 为标量时按第 1 分量（轴向）解释（``gap_solver._as_2vec`` 口径）⇒ 与 (G,0) 逐位同
    hb.set_pressure_gradient(HB_G)
    assert np.array_equal(hb.mobility(c, 1.0, 1.4, 0.9, H), I1)
    # 双向 G（第 2 分量非零）与轴向不同（逐点向量场确实进了求解）
    hb.set_pressure_gradient(np.stack([np.full((ny, nz), HB_G), np.full((ny, nz), 0.3 * HB_G)]))
    assert not np.array_equal(hb.mobility(c, 1.0, 1.4, 0.9, H), I1)
    # 标量 c̄/H（标量 G 注入）⇒ 0 维数组（与 NewtonianClosure 的返回类型一致）
    hb.set_pressure_gradient(HB_G)
    out = hb.mobility(0.5, 1.0, 1.4, 0.9, HB_H)
    assert isinstance(out, np.ndarray) and out.shape == ()
    # 场形状的 G 配标量 c̄ ⇒ 形状不匹配，显式 ValueError（不静默广播）
    hb.set_pressure_gradient(axial)
    with pytest.raises(ValueError, match="G"):
        hb.mobility(0.5, 1.0, 1.4, 0.9, HB_H)


def test_hb_closure_cache_reuses_solutions_and_is_keyed_by_gradient():
    """缓存：同一 (c̄,G,Gb,H) 只解一次；``buoyant_mobility`` 复用 ``mobility`` 的解。

    ⚠️ 键含 **G**（不只是 c̄）——``(c̄,B)`` 表只在固定 G 下自洽（裁定），故 G 不同
    必须重解；G 场均匀时缓存自动退化为"按唯一 c̄ 值缓存"。
    """
    c = np.array([0.2, 0.5, 0.2, 0.5])       # 两个唯一 c̄
    H = np.full(4, HB_H)
    hb = HBClosure(n=HB_N, kappa=HB_KAPPA, tau_y=HB_TAU_Y, m=1.0, B=0.3)
    hb.set_pressure_gradient(HB_G)
    hb.mobility(c, 1.0, 1.4, 0.9, H)
    s1 = hb.cache_stats
    assert s1 == {"hits": 2, "misses": 2}    # 4 点 → 2 解 + 2 命中
    hb.buoyant_mobility(c, 1.0, 1.4, 0.9, H)
    s2 = hb.cache_stats
    assert s2["hits"] == 6 and s2["misses"] == 2   # 复用：不再新增解

    # G 场改变 ⇒ 键改变 ⇒ 必须重解（不得拿旧 G 的闭包冒充）
    hb.set_pressure_gradient(2.0 * HB_G)
    hb.mobility(c, 1.0, 1.4, 0.9, H)
    assert hb.cache_stats["misses"] == 4