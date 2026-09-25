"""A3 惯例：**偏离默认值**且当前路径无消费者的开关必须一次性告警；取默认值/有消费者时不得告警。

R42（2026-09-25 修复轮 1）规则：告警的价值是"用户主动要求了一件不会发生的事"，故
`_dead_switches` 只计"**取值 != `_SWITCH_DEFAULTS` 声明默认值 且 当前路径无消费者**"；
取默认值属基线状态（不是预期落空），默认构造 `AnnulusD2DGASolver()` 必须零告警。
"""
import inspect
import warnings

import pytest

from cemdisp.models2d.annulus_d2dga import (
    AnnulusD2DGASolver,
    _SWITCH_DEFAULTS,
    _dead_switches,
)

OLD_PATH_ONLY = ["enable_regime_split", "enable_true_buoyancy", "enable_power_law_gap_law"]
NEW_PATH_ONLY = ["enable_stream_yield_gate", "enable_power_law_gap_correction"]


# --------------------------------------------------------------------------- #
# 0. 单一真源：签名默认值 == _SWITCH_DEFAULTS（防两处漂移）
# --------------------------------------------------------------------------- #

def test_signature_defaults_match_table():
    """签名实际默认值必须等于表中声明值（禁止签名与守卫各写一份字面量）。"""
    params = inspect.signature(AnnulusD2DGASolver.__init__).parameters
    for name, declared in _SWITCH_DEFAULTS.items():
        assert name in params, f"{name} 应仍是 __init__ 形参"
        assert params[name].default == declared, (
            f"{name} 签名默认值 {params[name].default!r} != _SWITCH_DEFAULTS 声明值 {declared!r}"
        )


# --------------------------------------------------------------------------- #
# 1. 旧路径专属开关：新路径下"偏离默认值"判死 / "取默认值"不判死
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("name", OLD_PATH_ONLY)
def test_old_path_switch_dead_on_new_path_when_deviating(name):
    """偏离默认值 ⇒ 判死（如 enable_regime_split=True 偏离其默认 False）。"""
    dead = _dead_switches(enable_stream_function=True,
                          **{name: not _SWITCH_DEFAULTS[name]})
    assert any(name in d for d in dead)


@pytest.mark.parametrize("name", OLD_PATH_ONLY)
def test_old_path_switch_not_dead_at_default_on_new_path(name):
    """取默认值 ⇒ 属基线空转（非预期落空），不判死。"""
    dead = _dead_switches(enable_stream_function=True,
                          **{name: _SWITCH_DEFAULTS[name]})
    assert not any(name in d for d in dead)


@pytest.mark.parametrize("name", OLD_PATH_ONLY)
@pytest.mark.parametrize("val", [True, False])
def test_old_path_switch_is_alive_on_old_path(name, val):
    dead = _dead_switches(enable_stream_function=False, **{name: val})
    assert not any(name in d for d in dead)


# --------------------------------------------------------------------------- #
# 2. 新路径专属开关（R43 保留的旧路径分支）：偏离默认值判死 / 取默认值不判死
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("name", NEW_PATH_ONLY)
def test_new_path_only_switch_dead_on_old_path_when_deviating(name):
    dead = _dead_switches(enable_stream_function=False,
                          **{name: not _SWITCH_DEFAULTS[name]})
    assert any(name in d for d in dead)


@pytest.mark.parametrize("name", NEW_PATH_ONLY)
def test_new_path_only_switch_not_dead_at_default_on_old_path(name):
    dead = _dead_switches(enable_stream_function=False,
                          **{name: _SWITCH_DEFAULTS[name]})
    assert not any(name in d for d in dead)


# --------------------------------------------------------------------------- #
# 3. 屈服门进算子：只有"输入被关掉"才算死
# --------------------------------------------------------------------------- #

def test_stream_yield_gate_dead_when_its_input_off():
    """用户打开了 stream_yield_gate（偏离默认）却关掉其输入 enable_yield_gate ⇒ 告警。"""
    dead = _dead_switches(enable_stream_function=True, enable_yield_gate=False,
                          enable_stream_yield_gate=True)
    assert any("enable_stream_yield_gate" in d for d in dead)


def test_default_yield_gate_combo_is_not_dead():
    """默认组合（enable_yield_gate=True 且 stream_yield_gate=False）两者皆默认 ⇒ 静默。"""
    dead = _dead_switches(enable_stream_function=True, enable_yield_gate=True,
                          enable_stream_yield_gate=False)
    assert not any("enable_yield_gate" in d for d in dead)


def test_stream_yield_gate_alive_when_input_on():
    dead = _dead_switches(enable_stream_function=True, enable_yield_gate=True,
                          enable_stream_yield_gate=True)
    assert not any("enable_stream_yield_gate" in d for d in dead)


# --------------------------------------------------------------------------- #
# 3b. R50：HB 路径显式丢弃 wall ⇒ stream_yield_gate 同样无消费方
# --------------------------------------------------------------------------- #

def test_stream_yield_gate_dead_on_hb_path():
    """HB 非线性入口无 wall 形参（:1788-1796）⇒ 开了它也无动力学作用。"""
    dead = _dead_switches(enable_stream_function=True, enable_yield_gate=True,
                          enable_stream_yield_gate=True, enable_hb_closure=True)
    assert any("enable_stream_yield_gate" in d for d in dead)
    assert any("enable_hb_closure=True" in d for d in dead)


def test_stream_yield_gate_alive_when_hb_off():
    """HB 关（默认）时，同一组合仍是活配置 ⇒ 不得告警（R50 的反例）。"""
    dead = _dead_switches(enable_stream_function=True, enable_yield_gate=True,
                          enable_stream_yield_gate=True, enable_hb_closure=False)
    assert not any("enable_stream_yield_gate" in d for d in dead)


# --------------------------------------------------------------------------- #
# 3c. R49：enable_local_i3 是唯一"默认关 + 置真后仍可被静默忽略"的开关
#     （消费点 :2282 位于 :2269 的 `enable_d2dga_i3_flux and enable_d2dga` 门内）
# --------------------------------------------------------------------------- #

def test_local_i3_dead_when_i3_flux_off():
    dead = _dead_switches(enable_local_i3=True, enable_d2dga_i3_flux=False,
                          enable_d2dga=True)
    assert any("enable_local_i3" in d for d in dead)
    assert any("enable_d2dga_i3_flux=False" in d for d in dead)
    # 上游开关本身是合法回退，不得被判死（R49 末条）
    assert not any(d.startswith("enable_d2dga") for d in dead)


def test_local_i3_dead_when_d2dga_off():
    dead = _dead_switches(enable_local_i3=True, enable_d2dga_i3_flux=True,
                          enable_d2dga=False)
    assert any("enable_local_i3" in d for d in dead)
    assert any("enable_d2dga=False" in d for d in dead)
    assert not any(d.startswith("enable_d2dga") for d in dead)


def test_local_i3_alive_when_all_three_on():
    """三者全开的反例：不告警。"""
    dead = _dead_switches(enable_local_i3=True, enable_d2dga_i3_flux=True,
                          enable_d2dga=True)
    assert dead == []


def test_local_i3_default_off_is_not_dead():
    """默认关（未偏离默认值）⇒ 静默，即使上游被回退关掉。"""
    dead = _dead_switches(enable_local_i3=False, enable_d2dga_i3_flux=False,
                          enable_d2dga=False)
    assert not any("enable_local_i3" in d for d in dead)


# --------------------------------------------------------------------------- #
# 4. 构造期行为
# --------------------------------------------------------------------------- #

def test_constructor_warns_once_and_names_switch():
    """偏离默认值的死开关仍必告警并点名（招牌用例：regime_split 默认 False）。"""
    with pytest.warns(UserWarning, match="enable_regime_split"):
        AnnulusD2DGASolver(total_t=1.0, nz=4, ny=9, enable_regime_split=True)


def test_constructor_quiet_when_all_alive():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        AnnulusD2DGASolver(total_t=1.0, nz=4, ny=9, enable_stream_function=False,
                           enable_regime_split=True)


def test_default_construction_is_silent():
    """默认构造零告警（R42）：把"默认静默"从偶然变成有意契约。"""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        AnnulusD2DGASolver()
        AnnulusD2DGASolver(nz=4, ny=9, total_t=1.0)


# --------------------------------------------------------------------------- #
# 5. R49/R50/R51 的构造期行为（生产路径实际发生的组合）
# --------------------------------------------------------------------------- #

def test_constructor_warns_when_local_i3_upstream_off():
    """R49：生产脚本（rerun_all_wells_corrected / ht1_004_ablation）传 enable_local_i3=True，
    上游一关就会被静默忽略 ⇒ 必须点名告警。"""
    with pytest.warns(UserWarning, match="enable_local_i3"):
        AnnulusD2DGASolver(total_t=1.0, nz=4, ny=9, enable_local_i3=True,
                           enable_d2dga_i3_flux=False)
    with pytest.warns(UserWarning, match="enable_local_i3"):
        AnnulusD2DGASolver(total_t=1.0, nz=4, ny=9, enable_local_i3=True,
                           enable_d2dga=False)


def test_constructor_quiet_when_local_i3_all_upstream_on():
    """R49 反例：三者全开（=生产 corrected 口径）必须静默。"""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        AnnulusD2DGASolver(total_t=1.0, nz=4, ny=9, enable_local_i3=True)


def test_constructor_warns_when_stream_yield_gate_on_hb_path():
    """R50：构造期语义须与运行期（HB 丢弃 wall 的幂等 UserWarning）一致。"""
    with pytest.warns(UserWarning, match="enable_stream_yield_gate"):
        AnnulusD2DGASolver(total_t=1.0, nz=4, ny=9, enable_stream_yield_gate=True,
                           enable_hb_closure=True)


def test_hb_on_old_path_does_not_claim_replacement_of_b3():
    """R51：旧路径上 HB 根本没运行，不得再宣称"HB 替代 B-3"（与"B-3 走错路径"互相矛盾）；
    但两个开关仍须各自被点名一次。"""
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        AnnulusD2DGASolver(total_t=1.0, nz=4, ny=9, enable_stream_function=False,
                           enable_hb_closure=True, enable_power_law_gap_correction=True)
    msgs = [str(r.message) for r in rec if issubclass(r.category, UserWarning)]
    assert any("enable_hb_closure" in m for m in msgs)
    assert any("enable_power_law_gap_correction" in m for m in msgs)
    assert not any("替代" in m for m in msgs)
