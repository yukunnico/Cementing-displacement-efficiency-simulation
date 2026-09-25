# -*- coding: utf-8 -*-
"""闭式 Ĩ₁ 快路径的契约测试（2026-09-18）。

守护三件事：
1. **判据正确**：`_i1_closed_applicable` 只在「壁面带 n₁≡1（宾汉）且中线带 τ_Y2≡0」
   时为真；H3 臂（注入水泥 τ_y）必须为假 ⇒ 自动回落求积（不静默外推）。
2. **闭式 = 求积**（同一物理量、两种算法）：在 H1 生产口径下相对差 < 1e-8。
   实测（10000 格真实呼101 场、逐点变 H、c̄∈[0,1]）GL-16 口径 max rel = 3.8e-9，
   该残差是**求积自身的误差**（被积函数 ∝ ỹ^{1/n₂+1} 在 ỹ=0 端点奇异）。
3. **退化/边界格**：g̃=0（零驱动/冻结）⇒ 0；c̄=0 / c̄=1 与求积一致。
"""
from __future__ import annotations

import numpy as np
import pytest

import cemdisp.models2d.gap_solver as gs


def _case(n1=1.0, n2=0.844, k1=0.058, k2=0.381, t1=9.2, t2=0.0, N=400):
    """构造一批格点：真实量级的 (κ, n, τ_Y, H) + 覆盖 [0,1] 的 c̄ + 覆盖量级的 g̃。"""
    rng = np.random.default_rng(20260918)
    H = 0.0069 + 0.038 * rng.random(N)                 # 半隙 6.9–45.6 mm（呼101 实测域）
    n_arr = np.broadcast_to(np.array([n1, n2]), (N, 2))
    kappa_t = np.broadcast_to(np.array([k1, k2]), (N, 2)) / H[:, None] ** n_arr
    tau = np.stack([np.full(N, t1), np.full(N, t2)], axis=1)
    c = np.linspace(0.0, 1.0, N)
    g = np.full(N, 766.2)
    return c, n_arr, kappa_t, tau, g


def _quad(c, n_arr, kappa_t, tau, g, monkeypatch):
    """强制回落求积（用于对照）。"""
    monkeypatch.setattr(gs, "_i1_closed_applicable", lambda *a, **k: False)
    return gs._i1_tilde_axial_batch(c, n_arr, kappa_t, tau, g)


def test_applicability_judge():
    _, n_arr, _, tau_h1, _ = _case(t2=0.0)
    assert gs._i1_closed_applicable(n_arr, tau_h1) is True
    _, n_arr, _, tau_h3, _ = _case(t2=0.97)          # H3 臂：水泥 τy 注入
    assert gs._i1_closed_applicable(n_arr, tau_h3) is False
    _, n_arr, _, tau_h1, _ = _case(n1=0.7, t2=0.0)   # 壁面带非宾汉
    assert gs._i1_closed_applicable(n_arr, tau_h1) is False


def test_closed_form_matches_quadrature(monkeypatch):
    c, n_arr, kappa_t, tau, g = _case(t2=0.0)        # H1 生产口径
    v_closed = gs._i1_tilde_axial_batch(c, n_arr, kappa_t, tau, g)
    v_quad = _quad(c, n_arr, kappa_t, tau, g, monkeypatch)
    rel = np.abs(v_closed - v_quad) / np.maximum(np.abs(v_quad), 1e-300)
    assert np.max(rel) < 1e-8, f"闭式 vs 求积 max rel = {np.max(rel):.3e}"


def test_closed_form_is_bitwise_identical_to_layer_decomposition():
    """模块分派路径 == 手工两段闭式之和（逐位）。"""
    c, n_arr, kappa_t, tau, g = _case(t2=0.0)
    v = gs._i1_tilde_axial_batch(c, n_arr, kappa_t, tau, g)
    i2, i1 = gs._i1_layer_pieces_closed(c, kappa_t, tau, n_arr, g)
    assert np.array_equal(v, i2 + i1)


def test_zero_drive_returns_zero(monkeypatch):
    """g̃ = 0（零驱动/冻结格）⇒ Ĩ₁ = 0，与求积路径一致。

    ⚠️ 契约域：``g̃ ≥ 0``。反求的支架下界恒为 0（``lo = 0``），生产路径不会产生
    负 ``g̃``；负值不在本函数契约内（求积路径取 ``|τ̃|``，对负 g̃ 的行为未被定义）。
    """
    c, n_arr, kappa_t, tau, g = _case(t2=0.0, N=50)
    g = g.copy(); g[:10] = 0.0
    v_closed = gs._i1_tilde_axial_batch(c, n_arr, kappa_t, tau, g)
    v_quad = _quad(c, n_arr, kappa_t, tau, g, monkeypatch)
    assert np.all(v_closed[:10] == 0.0)
    rel = np.abs(v_closed - v_quad) / np.maximum(np.abs(v_quad), 1e-300)
    assert np.max(rel) < 1e-8, f"max rel = {np.max(rel):.3e}"


def test_endpoints_c_zero_and_one(monkeypatch):
    """c̄ = 0（全泥浆）与 c̄ = 1（全水泥）两端点两路径一致。"""
    c, n_arr, kappa_t, tau, _ = _case(t2=0.0, N=2)
    c = np.array([0.0, 1.0])
    g = np.array([766.2, 766.2])
    v_closed = gs._i1_tilde_axial_batch(c, n_arr, kappa_t, tau, g)
    v_quad = _quad(c, n_arr, kappa_t, tau, g, monkeypatch)
    rel = np.abs(v_closed - v_quad) / np.maximum(np.abs(v_quad), 1e-300)
    # 容差取 1e-8：残差来源是求积（GL-16）而非闭式，实测该口径 max rel ≈ 3.8e-9。
    assert np.max(rel) < 1e-8, f"端点 rel = {np.max(rel):.3e}"


def test_h3_caliper_never_uses_closed_form(monkeypatch):
    """τ_Y2 ≠ 0 时必须走求积（闭式推导不成立），两路径结果一致且判据为假。"""
    c, n_arr, kappa_t, tau, g = _case(t2=0.97)
    assert gs._i1_closed_applicable(n_arr, tau) is False
    calls = {"n": 0}
    orig = gs._i1_layer_pieces_closed

    def spy(*a, **k):
        calls["n"] += 1
        return orig(*a, **k)

    monkeypatch.setattr(gs, "_i1_layer_pieces_closed", spy)
    _ = gs._i1_tilde_axial_batch(c, n_arr, kappa_t, tau, g)
    assert calls["n"] == 0, "H3 口径不应进入闭式路径"


# --------------------------------------------------------------------------- #
# I₂ 与 q₀ 的闭式（2026-09-18 追加）
# --------------------------------------------------------------------------- #
def _closure_case(t2=0.0, N=500, seed=20260918):
    """`closure_integrals_batch` 的输入三元组（真实量级；覆盖 c̄∈[0,1]）。"""
    rng = np.random.default_rng(seed)
    H = 0.0069 + 0.038 * rng.random(N)
    n_arr = np.broadcast_to(np.array([1.0, 0.844]), (N, 2))
    kap = np.broadcast_to(np.array([0.058, 0.381]), (N, 2))
    tau = np.stack([np.full(N, 9.2), np.full(N, t2)], axis=1)
    c = np.linspace(0.0, 1.0, N)
    G = np.full(N, 766.2)
    return c, n_arr, kap, tau, G, H


def test_i2_and_q0_match_quadrature(monkeypatch):
    """I₂ 与 q₀ 闭式 = 求积（H1 生产口径）。残差来源是求积，实测 max < 1e-8。"""
    c, n_arr, kap, tau, G, H = _closure_case(t2=0.0)
    a = gs.closure_integrals_batch(c, n_arr, kap, tau, G, H=H)
    monkeypatch.setattr(gs, "_i1_closed_applicable", lambda *x, **k: False)
    b = gs.closure_integrals_batch(c, n_arr, kap, tau, G, H=H)
    for name, u, v in (("I₁", a.I1, b.I1), ("I₂", a.I2, b.I2), ("q₀", a.q0, b.q0)):
        rel = np.abs(u - v) / np.maximum(np.abs(v), 1e-300)
        assert np.max(rel) < 1e-8, f"{name} max rel = {np.max(rel):.3e}"
    assert np.array_equal(a.undefined, b.undefined)


def test_q0_component_caliber():
    """q₀ 的分量口径回归锚：(A + c̄·C)/(A + B)，不得回退到含 τ̃ 的错误组合。

    ⚠️ 这是 2026-09-18 踩过的坑：`_gap_integral_vec_batch` 的被积核含 τ̃ 因子
    （故其 den2/den1/int1_one 各带一个 g̃），错用会把 q₀ 算错 ~10%。
    """
    c, n_arr, kap, tau, G, H = _closure_case(t2=0.0, N=4)   # c̄ = 0, 1/3, 2/3, 1
    c = np.array([0.0, 0.25, 0.5, 1.0])
    H = np.full(4, 0.02); G = np.full(4, 766.2)
    n_arr = np.broadcast_to(np.array([1.0, 0.844]), (4, 2))
    kap = np.broadcast_to(np.array([0.058, 0.381]), (4, 2))
    tau = np.stack([np.full(4, 9.2), np.zeros(4)], axis=1)
    got = gs.closure_integrals_batch(c, n_arr, kap, tau, G, H=H).q0
    # 独立复算（用模块的求积取 A/B/C，再按标定式组合）
    # ⚠️ closure_integrals_batch 收**量纲** G，内部换算 g̃ = H·G（见 _batch_inputs）；
    #    下面的手工复算必须用同一 g̃，否则会拿两个不同应力尺度对比。
    kt = kap / H[:, None] ** n_arr
    axial = np.ones(4, bool)
    Gt = np.stack([G * H, np.zeros(4)], 1)
    c2, d2, c1, d1 = gs._stress_vectors_batch(c, Gt, np.zeros((4, 2)))
    z, o = np.zeros(4), np.ones(4)
    A = gs._gap_integral_batch(z, c, c2, d2, kt[:, 1], n_arr[:, 1], tau[:, 1], axial, gs._weight_y2)
    B = gs._gap_integral_batch(c, o, c1, d1, kt[:, 0], n_arr[:, 0], tau[:, 0], axial, gs._weight_y2)
    C = gs._gap_integral_batch(c, o, c1, d1, kt[:, 0], n_arr[:, 0], tau[:, 0], axial, lambda y: y)
    expect = (A + c * C) / np.maximum(A + B, 1e-300)
    # 容差 1e-8：expect 由模块**求积**给出、got 由**闭式**给出，差即求积误差（实测 ~2e-9）
    assert np.allclose(got, expect, rtol=1e-8, atol=1e-10), f"q₀={got} 期望={expect}"
    assert got[0] == 0.0 and got[-1] == 1.0, "端点必须精确为 0 / 1"


def test_q0_endpoints_and_zero_drive(monkeypatch):
    """q₀ 边界：c̄=0 ⇒ 0；c̄=1 ⇒ 1；g̃=0 ⇒ 0——闭式与求积一致。"""
    c = np.array([0.0, 1.0, 0.5, 0.5])
    H = np.full(4, 0.02)
    n_arr = np.broadcast_to(np.array([1.0, 0.844]), (4, 2))
    kap = np.broadcast_to(np.array([0.058, 0.381]), (4, 2))
    tau = np.stack([np.full(4, 9.2), np.zeros(4)], axis=1)
    G = np.array([766.2, 766.2, 0.0, 1e-6])
    a = gs.closure_integrals_batch(c, n_arr, kap, tau, G, H=H)
    monkeypatch.setattr(gs, "_i1_closed_applicable", lambda *x, **k: False)
    b = gs.closure_integrals_batch(c, n_arr, kap, tau, G, H=H)
    assert a.q0[0] == 0.0 and a.q0[1] == 1.0
    # ⚠️ g̃=0（零驱动）时 Ĩ₁=0 ⇒ closure_integrals_batch 按设计把三个量都置 NaN
    #    （标量路径的 undefined 判据），**不是** 0——两条路径必须同为 NaN。
    assert np.isnan(a.q0[2]) and np.isnan(b.q0[2])
    assert np.array_equal(a.undefined, b.undefined)
    ok = ~a.undefined
    assert np.max(np.abs(a.q0[ok] - b.q0[ok])) < 1e-8


def test_h3_caliper_never_uses_i2_q0_closed(monkeypatch):
    """τ_Y2 ≠ 0 时 I₂/q₀ 同样必须走求积（闭式不成立）。"""
    c, n_arr, kap, tau, G, H = _closure_case(t2=0.97)
    assert gs._i1_closed_applicable(n_arr, tau) is False
    calls = {"n": 0}
    orig = gs._i2_q0_closed_extra

    def spy(*a, **k):
        calls["n"] += 1
        return orig(*a, **k)

    monkeypatch.setattr(gs, "_i2_q0_closed_extra", spy)
    _ = gs.closure_integrals_batch(c, n_arr, kap, tau, G, H=H)
    assert calls["n"] == 0, "H3 口径不应进入 I₂/q₀ 闭式路径"
