"""B-3：幂律间隙一阶修正（n=1 逐位退化；n<1 增强偏心分流）。"""
import warnings

import numpy as np
import pytest

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
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


@pytest.mark.parametrize("bad", [0.0, -0.7, float("nan")])
def test_non_positive_n_rejected(bad):
    """M-1：n ≤ 0（含 NaN）会使 1/n 无定义——构造即拒。"""
    with pytest.raises(ValueError):
        PowerLawGapClosure(bad)


def _tiny_geom(ny=5, nz=4):
    return {"phi": np.linspace(0.0, 1.0, ny), "s": np.linspace(0.0, 30.0, nz),
            "H": np.full((ny, nz), 0.01), "hole_mm": np.full(nz, 215.9),
            "od_mm": np.full(nz, 168.3), "y": np.linspace(0.0, np.pi * 0.1, ny),
            "b": np.full((ny, nz), 0.02)}


def test_bingham_cement_warns_correction_is_inert():
    """I-1：水泥相非幂律（power_law_n is None）时开关静默空转——须一次性告警。"""
    ny, nz = 5, 4
    geom = _tiny_geom(ny, nz)
    mud = FluidSpec("泥浆", FluidRole.MUD, 1200.0, RheologyModel.POWER_LAW,
                    power_law_n=0.8, consistency_k=0.5)
    bingham = FluidSpec("尾浆", FluidRole.TAIL, 1900.0, RheologyModel.BINGHAM,
                        plastic_viscosity_pa_s=0.17, yield_stress_pa=15.0)
    solver = AnnulusD2DGASolver(ny=ny, nz=nz, enable_power_law_gap_correction=True)
    lead = np.full((ny, nz), 0.5); tail = np.full((ny, nz), 0.5)
    w_prev = np.full((ny, nz), 0.5)

    with pytest.warns(UserWarning, match="静默空转"):
        solver._velocity_stream_function(lead, tail, geom, 0.02, w_prev,
                                         mud, bingham, bingham)
    # 第二次调用不再告警（一次性护栏）
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        solver._velocity_stream_function(lead, tail, geom, 0.02, w_prev,
                                         mud, bingham, bingham)
    assert not [r for r in rec if "静默空转" in str(r.message)]


def test_power_law_cement_does_not_warn():
    """幂律水泥相：修正真正生效，不告警。"""
    ny, nz = 5, 4
    geom = _tiny_geom(ny, nz)
    mud = FluidSpec("泥浆", FluidRole.MUD, 1200.0, RheologyModel.POWER_LAW,
                    power_law_n=0.8, consistency_k=0.5)
    cement = FluidSpec("尾浆", FluidRole.TAIL, 1900.0, RheologyModel.POWER_LAW,
                       power_law_n=0.7, consistency_k=0.9)
    solver = AnnulusD2DGASolver(ny=ny, nz=nz, enable_power_law_gap_correction=True)
    lead = np.full((ny, nz), 0.5); tail = np.full((ny, nz), 0.5)
    w_prev = np.full((ny, nz), 0.5)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        solver._velocity_stream_function(lead, tail, geom, 0.02, w_prev,
                                         mud, cement, cement)
    assert not [r for r in rec if "静默空转" in str(r.message)]
