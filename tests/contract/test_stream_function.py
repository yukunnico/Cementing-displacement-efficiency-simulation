"""stream_function 流函数椭圆方程求解契约测试（Z&F22 (2.3)/(4.22)/(2.2)）。

前两条为 Task 8 简报**逐字给定**的验收测试（BC+单调性 / 每列通量与深度
无关）；第三条为 **controller 裁定 R29（2026-09-15）改写版**：原简报版
（均匀 c̄ + 均匀 b_field）在 (4.22) 字面分组下源项恒为零、不可判别，改用
c̄ 带 φ-梯度的界面形态 + 忠实 (2,ny,nz) 全向量路径激励 A32-22 页的竖直井
真机制（"the elliptic equation for the stream function is driven by
gradients in b... these gradients arise from φ-gradients of the
concentration-dependent expression"）。

其后为本实现补充的边界情形：

- b_field 非零时每列通量仍守恒（∫2r_aH·w̄ dφ = 1 逐列成立）；
- ξ 两端零梯度下，ξ-均匀数据必给 ξ-均匀解（Neumann 行自洽性）；
- 方程对 Ψ 线性（变系数线性 Poisson 的叠加原理——本 Task 的关键性质）；
- (2.2) 第二式 v̄ = −∂ξΨ/(2H) 的契约锁死（ξ 向变 Ψ）。

b_field 一律走 (4.22) 字面分组 (2,ny,nz) 全向量（R29：唯一忠实口径，
(ny,nz) 标量路径已删除）。
"""
import numpy as np
import pytest
from cemdisp.models2d.stream_function import (solve_stream_function,
                                              velocity_from_stream_function)

def _geom(ny=40, nz=5):
    y = np.linspace(0, np.pi * 0.1071, ny)
    phi = y / y[-1]
    H = np.broadcast_to(0.0458 * (1 + 0.3 * np.cos(np.pi * phi))[:, None], (ny, nz)).copy()
    return {"y": y, "phi": phi, "H": H, "b": 2 * H, "s": np.linspace(0, 10.0, nz),
            "hole_mm": np.full((1, nz), 260.0), "od_mm": np.full((1, nz), 168.3)}

def test_boundary_conditions_and_monotonicity():
    g = _geom()
    psi = solve_stream_function(g, np.zeros_like(g["H"]), np.ones_like(g["H"]),
                                np.ones_like(g["H"]), 1.0, np.zeros((2,) + g["H"].shape))
    assert psi[0, :] == pytest.approx(0.0, abs=1e-12)
    assert psi[-1, :] == pytest.approx(1.0, rel=1e-9)
    assert np.all(np.diff(psi[:, 0]) >= -1e-12)

def test_flux_is_depth_independent():
    """同心/无浮力时每列通量必须相同。"""
    g = _geom()
    c = np.full_like(g["H"], 0.5)
    psi = solve_stream_function(g, c, np.ones_like(c), np.ones_like(c), 1.0,
                                np.zeros((2,) + g["H"].shape))
    w, _ = velocity_from_stream_function(psi, g)
    flux = np.trapezoid(2.0 * g["H"] * w, x=g["phi"], axis=0)
    assert np.allclose(flux, flux[0], rtol=1e-6)

def test_buoyancy_shifts_flow_to_narrow_side():
    """竖直井浮力机制（Z&F22 A32-22 页）：椭圆方程由 b 的梯度驱动，均匀 c̄ 下
    源项为零；须给 c̄ 的 φ-梯度（宽边 c̄ 高、窄边低的界面形态），经
    χ = ρ − Δρ·I₂/(H·I₁) 产生 b_φ = χ·f_φ 的 φ-梯度，源项
    −(1/r_a)∂φ(χ·r_a·cosβ/F²) 才非零。等黏度（η₁=η₂、m=1）下
    I₂/(H·I₁) = c̄(1−c̄²)，取 ρ₁=1、Δρ_paper=−1（密集水泥）的解析形
    χ = 1 + (2c̄ − c̄³)——其中常数部分源项为零（自动验证其退化）。
    断言（判别性保住）：浮力开启后窄边（φ→1）流速份额上升（b=0 vs b≠0）。
    """
    g = _geom()
    # 宽边 c̄ 高、窄边低：c̄ = 0.5 + 0.2cos(πφ) ∈ [0.3, 0.7]（ξ-均匀）
    c = (0.5 + 0.2 * np.cos(np.pi * g["phi"]))[:, None] * np.ones_like(g["H"])
    z = np.zeros((2,) + g["H"].shape)
    chi = 1.0 + 2.0 * c - c**3          # χ 的解析形（随 φ 单调递减）
    b = z.copy()
    b[0] = 5.0 * chi                    # φ-槽：b_φ = χ·f_φ（cosβ 分量幅值折入 5.0）
    b[1] = 0.0                          # ξ-槽：竖直井 f_ξ = sin(πφ)sinβ = 0
    w0, _ = velocity_from_stream_function(
        solve_stream_function(g, c, np.ones_like(c), np.ones_like(c), 1.0, z), g)
    w1, _ = velocity_from_stream_function(
        solve_stream_function(g, c, np.ones_like(c), np.ones_like(c), 1.0, b), g)
    ny = g["H"].shape[0]
    narrow = slice(3 * ny // 4, ny)
    assert w1[narrow].mean() > w0[narrow].mean()


# ---------------------------------------------------------------------------
# 以下为本实现补充的边界情形（controller 裁定第 7 条授权 + R29 修复轮）
# ---------------------------------------------------------------------------
def test_flux_conserved_with_buoyancy():
    """b_field 非零且 c̄ 沿 ξ 变化（真实前缘形态）时，不变量
    ∫₀¹ 2r_a·H·w̄ dφ = Ψ(1)−Ψ(0) = 1 逐列成立（等价于 trapezoid(2H·w̄) = 1/r_a）。

    依据：内点中心差分 + φ 端一阶单侧 + 梯形积分构成差分-梯形恒等式
    trapezoid(∂φΨ) = Ψ(1)−Ψ(0)，与 b_field 无关——浮力只在列内重分配
    轴向通量，不改变每列总通量。
    """
    g = _geom()
    c = np.full_like(g["H"], 0.5)
    c = c * (0.8 + 0.4 * np.tanh(np.linspace(-1.5, 1.5, g["H"].shape[1]))[None, :])
    b = np.zeros((2,) + g["H"].shape)
    b[0] = 5.0 * np.sin(np.pi * g["phi"])[:, None]                    # b_φ（φ-梯度）
    b[1] = 0.3 * np.tanh(np.linspace(-1.0, 1.0, g["H"].shape[1]))[None, :]  # b_ξ（ξ-梯度）
    psi = solve_stream_function(g, c, np.ones_like(c), np.ones_like(c), 1.0, b)
    w, _ = velocity_from_stream_function(psi, g)
    flux = np.trapezoid(2.0 * g["H"] * w, x=g["phi"], axis=0)
    r_a = (260.0 + 168.3) / 4.0 / 1000.0
    assert np.allclose(flux, 1.0 / r_a, rtol=1e-9)

def test_xi_uniform_data_gives_xi_uniform_solution():
    """ξ 两端零梯度 + ξ-均匀数据 ⇒ 解 ξ-均匀（Neumann 行自洽、无非物理 ξ 漂移）。"""
    g = _geom()
    c = np.full_like(g["H"], 0.5)
    b = np.zeros((2,) + g["H"].shape)
    b[0] = 3.0 * np.sin(np.pi * g["phi"])[:, None]   # ξ-均匀（φ-梯度）浮力
    psi = solve_stream_function(g, c, np.ones_like(c), np.ones_like(c), 1.0, b)
    assert np.allclose(psi, psi[:, :1], atol=1e-10)

def test_linearity_in_b_field():
    """方程对 Ψ 线性：Ψ(b1+b2) = Ψ(b1)+Ψ(b2)−Ψ(0)。

    每步 I₁、I₂ 只依赖 c̄ 与 m（本时间步已知），矩阵与 b 无关、仅右端
    随 b 线性变化——这是"每步只需解一次变系数线性 Poisson"的核心性质。
    """
    g = _geom()
    c = np.full_like(g["H"], 0.5)
    ones = np.ones_like(c)
    z = np.zeros((2,) + g["H"].shape)
    b1 = z.copy()
    b1[0] = 2.0 * np.sin(np.pi * g["phi"])[:, None]
    b2 = z.copy()
    b2[0] = 7.0 * (0.5 + 0.2 * np.cos(np.pi * g["phi"]))[:, None]
    psi_0 = solve_stream_function(g, c, ones, ones, 1.0, z)
    psi_1 = solve_stream_function(g, c, ones, ones, 1.0, b1)
    psi_2 = solve_stream_function(g, c, ones, ones, 1.0, b2)
    psi_12 = solve_stream_function(g, c, ones, ones, 1.0, b1 + b2)
    assert np.allclose(psi_1 + psi_2 - psi_0, psi_12, atol=1e-10)

def test_vector_b_field_boundary_and_flux():
    """(2,ny,nz) 全向量路径（(4.22) 字面分组：φ-槽 b_φ、ξ-槽 b_ξ）
    同样满足 Dirichlet BC 与每列通量守恒。"""
    g = _geom()
    c = np.full_like(g["H"], 0.5)
    b = np.zeros((2,) + g["H"].shape)
    b[0] = 1.5 * np.sin(np.pi * g["phi"])[:, None]   # b_φ（φ-槽）
    b[1] = 0.3                                        # b_ξ（ξ-槽）
    psi = solve_stream_function(g, c, np.ones_like(c), np.ones_like(c), 1.0, b)
    assert psi[0, :] == pytest.approx(0.0, abs=1e-12)
    assert psi[-1, :] == pytest.approx(1.0, rel=1e-9)
    w, _ = velocity_from_stream_function(psi, g)
    flux = np.trapezoid(2.0 * g["H"] * w, x=g["phi"], axis=0)
    r_a = (260.0 + 168.3) / 4.0 / 1000.0
    assert np.allclose(flux, 1.0 / r_a, rtol=1e-9)

def test_v_is_negative_xi_derivative_over_2h():
    """(2.2) 第二式契约锁死：v̄ = −∂ξΨ/(2H)。

    构造 ξ 向变化的 Ψ：c̄ 沿 ξ 变化（前缘形态）且 **m=2**（m=1 时 I₁ 与 c̄
    无关，c̄ 的 ξ-变化对方程不可见——这正是需要黏度比的原因），使算子系数
    I₁ 沿 ξ 变化 ⇒ Ψ 沿 ξ 变化 ⇒ v̄ 非零。以独立转写的差分（内点中心、
    ξ 端一阶单侧）复算 −∂ξΨ/(2H) 对照。
    """
    g = _geom()
    c = np.full_like(g["H"], 0.5)
    c = c * (0.8 + 0.4 * np.tanh(np.linspace(-1.5, 1.5, g["H"].shape[1]))[None, :])
    b = np.zeros((2,) + g["H"].shape)
    b[0] = 5.0 * np.sin(np.pi * g["phi"])[:, None]
    psi = solve_stream_function(g, c, np.ones_like(c), np.ones_like(c), 2.0, b)
    w, v = velocity_from_stream_function(psi, g)
    assert w.shape == v.shape == psi.shape
    assert not np.allclose(psi, psi[:, :1], atol=1e-8)   # 前置：Ψ 确实沿 ξ 变化
    # 独立转写 −∂ξΨ/(2H)（内点中心差分，ξ 端一阶单侧）
    dxi = g["s"][1] - g["s"][0]
    dpsi = np.empty_like(psi)
    dpsi[:, 1:-1] = (psi[:, 2:] - psi[:, :-2]) / (2.0 * dxi)
    dpsi[:, 0] = (psi[:, 1] - psi[:, 0]) / dxi
    dpsi[:, -1] = (psi[:, -1] - psi[:, -2]) / dxi
    v_ref = -dpsi / (2.0 * g["H"])
    assert np.allclose(v, v_ref, rtol=1e-12)
    assert not np.allclose(v, 0.0)   # v̄ 真实非零（排除恒零假通过）
