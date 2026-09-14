"""stream_function 流函数椭圆方程求解契约测试（Z&F22 (2.3)/(4.22)/(2.2)）。

前三条为 Task 8 简报**逐字给定**的验收测试（BC+单调性 / 每列通量与深度
无关 / 浮力把流推向窄边）；其后为本实现补充的边界情形：

- b_field 非零时每列通量仍守恒（∫∂φΨ dφ = 1，controller 裁定第 7 条）；
- ξ 两端零梯度下，ξ-均匀数据必给 ξ-均匀解（Neumann 行自洽性）；
- 方程对 Ψ 线性（变系数线性 Poisson 的叠加原理——本 Task 的关键性质）；
- (2,ny,nz) 全向量路径（(4.22) 忠实分组）同样满足 BC 与每列通量守恒。
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
                                np.ones_like(g["H"]), 1.0, np.zeros_like(g["H"]))
    assert psi[0, :] == pytest.approx(0.0, abs=1e-12)
    assert psi[-1, :] == pytest.approx(1.0, rel=1e-9)
    assert np.all(np.diff(psi[:, 0]) >= -1e-12)

def test_flux_is_depth_independent():
    """同心/无浮力时每列通量必须相同。"""
    g = _geom()
    c = np.full_like(g["H"], 0.5)
    psi = solve_stream_function(g, c, np.ones_like(c), np.ones_like(c), 1.0, np.zeros_like(c))
    w, _ = velocity_from_stream_function(psi, g)
    flux = np.trapezoid(2.0 * g["H"] * w, x=g["phi"], axis=0)
    assert np.allclose(flux, flux[0], rtol=1e-6)

def test_buoyancy_shifts_flow_to_narrow_side():
    g = _geom()
    c = np.full_like(g["H"], 0.5)
    z = np.zeros_like(c)
    w0, _ = velocity_from_stream_function(
        solve_stream_function(g, c, np.ones_like(c), np.ones_like(c), 1.0, z), g)
    w1, _ = velocity_from_stream_function(
        solve_stream_function(g, c, np.ones_like(c), np.ones_like(c), 1.0, z + 5.0), g)
    ny = g["H"].shape[0]
    narrow = slice(3 * ny // 4, ny)
    assert w1[narrow].mean() > w0[narrow].mean()


# ---------------------------------------------------------------------------
# 以下为本实现补充的边界情形（brief Step 1 未覆盖、controller 裁定第 7 条授权）
# ---------------------------------------------------------------------------
def test_flux_conserved_with_buoyancy():
    """b_field 非零且 c̄ 沿 ξ 变化（真实前缘形态）时，每列通量仍恒等于 1/r_a。

    依据：内点中心差分 + φ 端一阶单侧 + 梯形积分构成差分-梯形恒等式
    trapezoid(∂φΨ) = Ψ(1)−Ψ(0) = 1，故 2H·w̄ = ∂φΨ/r_a 的每列积分恒为
    1/r_a，与 b_field 无关（浮力只在列内重分配通量）。
    """
    g = _geom()
    c = np.full_like(g["H"], 0.5)
    c = c * (0.8 + 0.4 * np.tanh(np.linspace(-1.5, 1.5, g["H"].shape[1]))[None, :])
    b = np.full_like(g["H"], 5.0)
    psi = solve_stream_function(g, c, np.ones_like(c), np.ones_like(c), 1.0, b)
    w, _ = velocity_from_stream_function(psi, g)
    flux = np.trapezoid(2.0 * g["H"] * w, x=g["phi"], axis=0)
    r_a = (260.0 + 168.3) / 4.0 / 1000.0
    assert np.allclose(flux, 1.0 / r_a, rtol=1e-9)

def test_xi_uniform_data_gives_xi_uniform_solution():
    """ξ 两端零梯度 + ξ-均匀数据 ⇒ 解 ξ-均匀（Neumann 行自洽、无非物理 ξ 漂移）。"""
    g = _geom()
    c = np.full_like(g["H"], 0.5)
    psi = solve_stream_function(g, c, np.ones_like(c), np.ones_like(c), 1.0,
                                np.full_like(g["H"], 3.0))
    assert np.allclose(psi, psi[:, :1], atol=1e-10)

def test_linearity_in_b_field():
    """方程对 Ψ 线性：Ψ(b1+b2) = Ψ(b1)+Ψ(b2)−Ψ(0)。

    每步 I₁、I₂ 只依赖 c̄ 与 m（本时间步已知），矩阵与 b 无关、仅右端
    随 b 线性变化——这是"每步只需解一次变系数线性 Poisson"的核心性质。
    """
    g = _geom()
    c = np.full_like(g["H"], 0.5)
    ones = np.ones_like(c)
    z = np.zeros_like(c)
    b1 = np.full_like(c, 2.0)
    b2 = np.full_like(c, 7.0)
    psi_0 = solve_stream_function(g, c, ones, ones, 1.0, z)
    psi_1 = solve_stream_function(g, c, ones, ones, 1.0, b1)
    psi_2 = solve_stream_function(g, c, ones, ones, 1.0, b2)
    psi_12 = solve_stream_function(g, c, ones, ones, 1.0, b1 + b2)
    assert np.allclose(psi_1 + psi_2 - psi_0, psi_12, atol=1e-10)

def test_vector_b_field_boundary_and_flux():
    """(2,ny,nz) 全向量路径（(4.22) 忠实分组：φ-槽 b_φ、ξ-槽 b_ξ）
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
