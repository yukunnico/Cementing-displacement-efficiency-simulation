"""B-1：闭包协议与牛顿实现（逐位等价 + 求解器注入无扰）。

对应设计规格 docs/superpowers/specs/2026-09-16-hb-closure-route2-design.md §3 阶段 B-1。
L1 硬约束：closure 缺省注入必须与 HEAD 的 two_layer.mobility_i1/i2 逐位一致。
"""
import numpy as np

from cemdisp.models2d.hb_closure import ClosureProvider, NewtonianClosure
from cemdisp.models2d.two_layer import mobility_i1, mobility_i2


def test_newtonian_closure_implements_protocol():
    assert isinstance(NewtonianClosure(), ClosureProvider)


def test_newtonian_closure_bitwise_equals_two_layer():
    c = np.linspace(0.0, 1.0, 7)
    H = np.full(7, 0.0123)
    n = NewtonianClosure()
    assert np.array_equal(
        n.mobility(c, 0.7, 0.058, 0.171, H),
        np.asarray(mobility_i1(c, 0.7, eta1=0.058, eta2=0.171, H=H)),
    )
    assert np.array_equal(
        n.buoyant_mobility(c, 0.7, 0.058, 0.171, H),
        np.asarray(mobility_i2(c, 0.7, eta1=0.058, eta2=0.171, H=H)),
    )


def _tiny_geom(ny=5, nz=4):
    phi = np.linspace(0.0, 1.0, ny)
    s = np.linspace(0.0, 30.0, nz)
    H = np.full((ny, nz), 0.01)
    return {
        "phi": phi,
        "s": s,
        "H": H,
        "hole_mm": np.full(nz, 215.9),
        "od_mm": np.full(nz, 168.3),
        "y": np.linspace(0.0, np.pi * 0.1, ny),
        "b": np.full((ny, nz), 0.02),
    }


def test_solve_stream_default_equals_explicit_newtonian():
    """closure=None 与 closure=NewtonianClosure() 逐位一致。"""
    from cemdisp.models2d.stream_function import solve_stream_function

    g = _tiny_geom()
    c = np.full((5, 4), 0.5)
    b = np.zeros((2, 5, 4))
    psi_a = solve_stream_function(g, c, 0.058, 0.171, 0.34, b)
    psi_b = solve_stream_function(g, c, 0.058, 0.171, 0.34, b, closure=NewtonianClosure())
    assert np.array_equal(psi_a, psi_b)
