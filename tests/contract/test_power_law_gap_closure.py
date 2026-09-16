"""B-3：幂律间隙一阶修正（n=1 逐位退化；n<1 增强偏心分流）。"""
import numpy as np

from cemdisp.models2d.hb_closure import PowerLawGapClosure, NewtonianClosure


def test_n_equals_one_is_bitwise_newtonian():
    c = np.linspace(0.0, 1.0, 7); H = np.linspace(0.008, 0.014, 7)
    a = PowerLawGapClosure(1.0).mobility(c, 0.7, 0.058, 0.171, H)
    b = NewtonianClosure().mobility(c, 0.7, 0.058, 0.171, H)
    assert np.array_equal(a, b)


def test_n_less_one_amplifies_gap_contrast():
    """n<1 时 I₁ 的 H 反差被放大：宽/窄流动度比增大。"""
    c = np.full(5, 0.5); H = np.array([0.006, 0.008, 0.010, 0.012, 0.014])
    i_newton = NewtonianClosure().mobility(c, 0.7, 0.058, 0.171, H)
    i_pl = PowerLawGapClosure(0.7).mobility(c, 0.7, 0.058, 0.171, H)
    assert (i_pl[-1] / i_pl[0]) > (i_newton[-1] / i_newton[0])
