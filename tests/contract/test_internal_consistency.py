# -*- coding: utf-8 -*-
"""自洽校验量（域内效率 / 饥饿份额 / 前缘位置）的契约测试（Task 6）。

口径声明（协调者裁定 #3，务必保留）：
    ``η_E ≡ 1 − 饥饿体积份额`` **只在 c̄ ∈ {0,1}（二值浓度场）时严格成立**。
    生产井浓度场连续，该式是**带阈值的记录**（历史实测偏差 ≤0.007），
    故此处拆成两条断言：二值场要求严格相等、非二值场要求"确实存在空隙"，
    **不得**对连续场写硬相等断言。
"""
import numpy as np
import pytest

from cemdisp.diagnostics.internal_consistency import (
    domain_eta_e,
    front_position_m,
    starved_volume_fraction,
)


def _geom(ny=9, nz=5):
    return {"phi": np.linspace(0, 1, ny), "s": np.linspace(0, 40, nz),
            "md": np.linspace(7000, 6960, nz), "H": np.full((ny, nz), 0.01),
            "b": np.full((ny, nz), 0.02), "y": np.linspace(0, 0.3, ny)}


def test_identity_exact_for_binary_field():
    g = _geom()
    c = np.zeros((9, 5))
    c[:5] = 1.0                       # 二值场 ⇒ 恒等式严格成立
    assert abs((1.0 - domain_eta_e(c, g)) - starved_volume_fraction(c, g)) < 1e-12


def test_identity_is_approximate_for_mixed_field():
    g = _geom()
    c = np.full((9, 5), 0.4)          # 非二值 ⇒ 有"空隙"，须以阈值记录而非硬断言
    assert abs((1.0 - domain_eta_e(c, g)) - starved_volume_fraction(c, g)) > 0.1


def test_front_position_finite_and_monotone():
    g = _geom()
    c = np.zeros((9, 5))
    c[:, :3] = 1.0                    # s 小的一侧（鞋口侧）已到水泥
    assert front_position_m(c, g) > 0.0


# --------------------------------------------------------------------------- #
# 追加：数值安全与退化面（与上面三条同一批交付，防"除零恒零"式假绿）
# --------------------------------------------------------------------------- #

def test_all_cement_binary_gives_full_efficiency_and_zero_starved():
    """全水泥二值场：η_E = 1、饥饿份额 = 0、恒等式偏差 = 0。"""
    g = _geom()
    c = np.ones((9, 5))
    assert domain_eta_e(c, g) == pytest.approx(1.0, abs=1e-12)
    assert starved_volume_fraction(c, g) == pytest.approx(0.0, abs=1e-12)


def test_zero_b_field_returns_finite_values():
    """b ≡ 0（退化几何）不得产生 nan/inf——分母由 1e-12 护栏保护。"""
    g = _geom()
    g["b"] = np.zeros((9, 5))
    c = np.full((9, 5), 0.5)
    assert np.isfinite(domain_eta_e(c, g))
    assert np.isfinite(starved_volume_fraction(c, g))


def test_front_position_falls_back_to_domain_bottom_when_no_cement():
    """全钻井液（未达 level）⇒ 按约定返回域底 md。"""
    g = _geom()
    c = np.zeros((9, 5))
    assert front_position_m(c, g) == pytest.approx(7000.0, abs=1e-12)
