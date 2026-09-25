# -*- coding: utf-8 -*-
"""自洽校验量（域内效率 / 饥饿份额 / 前缘位置）的契约测试（Task 6）。

口径声明（协调者裁定 #3，务必保留）：
    ``η_E ≡ 1 − 饥饿体积份额`` 是**代数恒等**——对任意连续场可严格展开为
    ``Σ_{c≥0.5} b(1−c)/∬b − Σ_{c<0.5} b·c/∬b``（过渡带的代数和，浮点复算一致到 ~1e-16）；
    该代数和**只在 c̄ ∈ {0,1}（二值浓度场）时恒为 0**。生产井浓度场连续，
    故该式是**带阈值的记录**（历史"偏差 ≤0.007"口径出自**基准算例**，非现场井），
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
    # s = 0, 10, 20, 30, 40（自鞋口起）；md = bottom − s（派生标签）
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


def test_front_position_advances_with_cement_advance():
    """s 口径前缘随水泥推进**单调增大**（旧版名为 monotone 却只断言 >0，名不副实）。"""
    g = _geom()                                   # s = 0, 10, 20, 30, 40
    less = np.zeros((9, 5))
    less[:, :3] = 1.0                             # 水泥占 s ≤ 20
    more = np.zeros((9, 5))
    more[:, :4] = 1.0                             # 推进更多：水泥占 s ≤ 30
    f_less, f_more = front_position_m(less, g), front_position_m(more, g)
    assert f_less == pytest.approx(20.0, abs=1e-12)
    assert f_more == pytest.approx(30.0, abs=1e-12)
    assert f_more > f_less > 0.0                  # 推进更多 ⇒ s 前缘更大


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


def test_front_position_returns_zero_when_no_cement():
    """整场无水泥（所有列均值未达 level）⇒ 兜底返回 0.0（与求解器 ``_front`` 同约定）。

    与"水泥只到鞋口 s=0"的正常分支同为 0.0（物理上是同一状态），但与"水泥推进到 s>0"
    的正常分支可区分；且若兜底分支被删（``reached`` 为空仍取 ``.max()``），
    第一条断言会以 ValueError 失败 ⇒ 两分支可辨。
    """
    g = _geom()
    c = np.zeros((9, 5))                          # 整场无水泥 ⇒ 走兜底分支
    assert front_position_m(c, g) == pytest.approx(0.0, abs=1e-12)
    c[:, 0] = 1.0                                 # 仅鞋口列有水泥 ⇒ 同为 0.0
    assert front_position_m(c, g) == pytest.approx(0.0, abs=1e-12)
    c_far = np.zeros((9, 5))
    c_far[:, -1] = 1.0                            # 仅出口侧列有水泥 ⇒ 前缘 = s 域顶（与兜底可辨）
    assert front_position_m(c_far, g) == pytest.approx(40.0, abs=1e-12)
