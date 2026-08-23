import numpy as np
from cemdisp.models2d import regime_closure as rc


def test_metzner_reed_re_shape_and_positive():
    w = np.full(4, 0.5); rho = np.full(4, 1500.0); n = np.full(4, 0.7)
    kappa = np.full(4, 0.4); b = np.full(4, 0.02)
    re = rc.metzner_reed_re(w, rho, n, kappa, b)
    assert re.shape == (4,) and np.all(re > 0)


def test_laminar_friction_is_24_over_re_when_no_yield():
    re = np.array([10.0, 100.0, 1000.0]); he = np.zeros(3); n = np.ones(3)
    f = rc.friction_laminar(re, he, n)
    np.testing.assert_allclose(f, 24.0 / re, rtol=1e-6)


def test_dodge_metzner_converges_and_monotone():
    re = np.logspace(3, 6, 50); n = np.full(50, 0.7)
    f = rc.friction_dodge_metzner(re, n)
    assert f.shape == (50,) and np.all(f > 0) and np.all(f < 1.0)
    assert np.all(np.diff(f) <= 1e-6)


def test_drag_weight_layer_mask_and_bounds():
    re = np.array([10.0, 500.0, 5000.0, 50000.0])
    he = np.zeros(4); n = np.full(4, 0.8)
    R, mask = rc.drag_weight(re, he, n, re_crit=np.full(4, 2100.0))
    assert R[0] == 1.0
    assert np.all(R >= 0.3) and np.all(R <= 1.5)
    assert mask.shape == (4,)


def test_hedstrom_zero_for_no_yield():
    tau_y = np.zeros(3); rho = np.full(3, 1500.0); n = np.full(3, 0.7)
    kappa = np.full(3, 0.4); b = np.full(3, 0.02)
    assert np.all(rc.hedstrom_number(tau_y, rho, n, kappa, b) == 0.0)
