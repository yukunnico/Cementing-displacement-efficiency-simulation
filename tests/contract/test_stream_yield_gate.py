"""B-2：屈服门接入流函数算子（wall=None 逐位无扰；wall>0 冻结窄边）。"""
import numpy as np

from cemdisp.models2d.stream_function import solve_stream_function, _WALL_CONDUCTANCE_FLOOR


def _geom(ny=9, nz=4):
    return {"phi": np.linspace(0.0, 1.0, ny), "s": np.linspace(0.0, 30.0, nz),
            "H": np.full((ny, nz), 0.01), "hole_mm": np.full(nz, 215.9),
            "od_mm": np.full(nz, 168.3), "y": np.linspace(0.0, np.pi * 0.1, ny),
            "b": np.full((ny, nz), 0.02)}


def test_wall_none_is_bitwise_unchanged():
    g = _geom(); c = np.full((9, 4), 0.5); b = np.zeros((2, 9, 4))
    a = solve_stream_function(g, c, 0.058, 0.171, 0.34, b)
    z = solve_stream_function(g, c, 0.058, 0.171, 0.34, b, wall=np.zeros((9, 4)))
    assert np.array_equal(a, z)


def test_wall_freezes_narrow_side():
    """窄边（φ→1）wall=1 ⇒ 该侧轴向速度显著下降。"""
    from cemdisp.models2d.stream_function import velocity_from_stream_function
    g = _geom(); c = np.full((9, 4), 0.5); b = np.zeros((2, 9, 4))
    wall = np.zeros((9, 4)); wall[-3:, :] = 1.0
    psi0 = solve_stream_function(g, c, 0.058, 0.171, 0.34, b)
    psi1 = solve_stream_function(g, c, 0.058, 0.171, 0.34, b, wall=wall)
    w0, _ = velocity_from_stream_function(psi0, g)
    w1, _ = velocity_from_stream_function(psi1, g)
    assert np.mean(np.abs(w1[-3:])) < np.mean(np.abs(w0[-3:]))
    assert np.all(np.isfinite(psi1))
    # controller 裁定（2026-09-16）：断言须落在冻结区**内部**（行 7–8）并核验物理
    # ——冻结 ⇒ w→0（量级显著小，而非仅"比基线小"）。行 6 是冻结/活跃交界，
    # 其梯度属正常界面过渡（实测 ~0.36×基线），不参与断言。
    assert np.mean(np.abs(w1[-2:])) < 1.0e-3 * np.mean(np.abs(w0[-2:]))


def test_wall_redistributes_flux_without_changing_total():
    """守恒不变量（设计规格 B-2 验收「守恒仍成立」）：屈服门只把通量从冻结区
    重新分配到活跃区，不改变每列总通量 ∫2r_aH·w̄ dφ = 1/r_a（单位 BC）。"""
    from cemdisp.models2d.stream_function import velocity_from_stream_function
    g = _geom(); c = np.full((9, 4), 0.5); b = np.zeros((2, 9, 4))
    wall = np.zeros((9, 4)); wall[-3:, :] = 1.0
    psi0 = solve_stream_function(g, c, 0.058, 0.171, 0.34, b)
    psi1 = solve_stream_function(g, c, 0.058, 0.171, 0.34, b, wall=wall)
    w0, _ = velocity_from_stream_function(psi0, g)
    w1, _ = velocity_from_stream_function(psi1, g)
    flux0 = np.trapezoid(2.0 * g["H"] * w0, x=g["phi"], axis=0)
    flux1 = np.trapezoid(2.0 * g["H"] * w1, x=g["phi"], axis=0)
    assert np.allclose(flux1, flux0, rtol=1.0e-9)      # 总通量逐列守恒不变
    assert np.mean(w1[:4]) > np.mean(w0[:4])           # 活跃（宽）边流速上升


def test_conductance_floor_prevents_singular():
    g = _geom(); c = np.full((9, 4), 0.5); b = np.zeros((2, 9, 4))
    psi = solve_stream_function(g, c, 0.058, 0.171, 0.34, b, wall=np.ones((9, 4)))
    assert np.all(np.isfinite(psi))
