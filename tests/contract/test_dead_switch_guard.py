"""A3 惯例：置真但当前路径无消费者的开关必须一次性告警；有消费者时不得告警。"""
import warnings

import pytest

from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver, _dead_switches

OLD_PATH_ONLY = ["enable_regime_split", "enable_true_buoyancy", "enable_power_law_gap_law"]


@pytest.mark.parametrize("name", OLD_PATH_ONLY)
def test_old_path_switch_is_dead_on_new_path(name):
    dead = _dead_switches(enable_stream_function=True, **{name: True})
    assert any(name in d for d in dead)


@pytest.mark.parametrize("name", OLD_PATH_ONLY)
def test_old_path_switch_is_alive_on_old_path(name):
    dead = _dead_switches(enable_stream_function=False, **{name: True})
    assert not any(name in d for d in dead)


def test_yield_gate_dead_until_stream_yield_gate_on():
    dead = _dead_switches(enable_stream_function=True, enable_yield_gate=True,
                          enable_stream_yield_gate=False)
    assert any("enable_yield_gate" in d for d in dead)
    alive = _dead_switches(enable_stream_function=True, enable_yield_gate=True,
                           enable_stream_yield_gate=True)
    assert not any("enable_yield_gate" in d for d in alive)


def test_constructor_warns_once_and_names_switch():
    with pytest.warns(UserWarning, match="enable_regime_split"):
        AnnulusD2DGASolver(total_t=1.0, nz=4, ny=9, enable_regime_split=True)


def test_constructor_quiet_when_all_alive():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        AnnulusD2DGASolver(total_t=1.0, nz=4, ny=9, enable_stream_function=False,
                           enable_regime_split=True)
