# -*- coding: utf-8 -*-
"""Task 11 契约测试：屈服门二值 → 连续（2026-09-15）。

文献锚点：停流区判据 τw < f·τy（Pelipenko04 (2.6)-(2.8)）的连续近似
``wall = clip(1 − τw_extrap/(f·τy), 0, 1)``。

旧口径 ``wall = np.where(τw_extrap ≤ f·τy, 1, 0)`` 是二值悬崖：相邻网格元
壁面剪应力跨过阈值时 wall 从 0 跳到 1（差 = 1.0），窄边冻结带呈锯齿状对
排量/浓度微扰过敏。本契约锁定新口径：

1. **连续性**：τw 跨阈值扫掠下相邻元 wall 差 < 0.5（旧二值在阈值处跳变
   1.0，必红）；
2. **极限**：τw→0（窄缝极限）⇒ wall→1（全冻）；τw ≫ f·τy ⇒ wall→0；
3. **R3 不变量①参考元**：每列流动最快元（|w| 最大且 w>0）恒 wall=0，
   即使其自身外推壁剪低于阈值；
4. **R3 不变量②整列冻结**：整列无流动且水泥已到 → 整列 wall=1；
5. **τy=0 零除防护语义**：f·τy = 0 ⇒ 无屈服应力 ⇒ wall=0（且输出无
   NaN/inf——连续式的 0/0 路径必须被防护）。

边界语义（docstring 级约定，此处锁定两侧极限）：τw = f·τy 恰 wall=0
（可流动）——旧二值 ``≤`` 在等号处冻结（wall=1），两口径仅在此零测度集
上不同；Pelipenko04 (2.6)-(2.8) 停流区判据为严格 τw < f·τy，故等号归
可流动。

构造口径：直接单测静态方法 ``_yield_gate_wall``。外推链为
τw_extrap_i = τw_ref·(b_i/b_ref)，τw_ref = mu·6·w_ref/b_ref（参考元取
row0，唯一 w>0），故令 b_i = b_ref·τw_i/τw_ref 即可精确钉住各行外推壁剪。
"""
from __future__ import annotations

import numpy as np
import pytest

from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver

_F_SAFETY = 1.15
_TAU_Y = 8.0
_THRESHOLD = _F_SAFETY * _TAU_Y  # 9.2


def _single_column_b(tau_w_targets_tail, b_ref=0.04, mu=0.1, w_ref=1.0):
    """构造单列 (ny,1) 场：row0 为参考元（唯一 w>0，b=b_ref）。

    参考元外推壁剪 = τw_ref = mu·6·w_ref/b_ref（自身定义）；
    第 i≥1 行外推壁剪 τw_extrap_i = tau_w_targets_tail[i-1]
    （反解 b_i = b_ref·τw_i/τw_ref）。

    Returns:
        (w, b, mu_reg, tau_w_ref)
    """
    tau_w_ref = mu * 6.0 * w_ref / b_ref
    tail = np.asarray(tau_w_targets_tail, dtype=float)
    b_tail = b_ref * tail / tau_w_ref
    b = np.concatenate(([b_ref], b_tail))[:, None]
    ny = b.shape[0]
    w = np.zeros((ny, 1))
    w[0, 0] = w_ref  # 仅参考元流动（外推口径与既有 history 测试同构）
    mu_reg = np.full((ny, 1), mu)
    return w, b, mu_reg, tau_w_ref


def _cement_fields(ny):
    cement_ever = np.ones((ny, 1))
    cement_local = np.full((ny, 1), 0.9)
    return cement_ever, cement_local


def _gate(w, b, mu_reg, tau_y, cement_ever, cement_local, f_safety=_F_SAFETY):
    return AnnulusD2DGASolver._yield_gate_wall(
        w, b, mu_reg, tau_y, cement_ever, cement_local, f_safety)


# ---------------------------------------------------------------------------
# 1. 连续性：跨阈值扫掠无 0/1 跳变（旧二值必红）
# ---------------------------------------------------------------------------
class TestContinuity:
    def test_adjacent_wall_difference_below_half_across_threshold(self):
        """τw 从 1.6·阈值 扫到 0.3·阈值（跨过阈值），相邻元 wall 差全程 < 0.5。

        扫掠必然包含阈值穿越对（旧二值：一侧 wall=0、另一侧 wall=1，
        相邻差恰 = 1.0 ≥ 0.5 → 旧代码本测试必红）。
        """
        tau_w_tail = np.linspace(1.6 * _THRESHOLD, 0.3 * _THRESHOLD, 28)
        w, b, mu_reg, _ = _single_column_b(tau_w_tail)
        cement_ever, cement_local = _cement_fields(b.shape[0])
        tau_y = np.full((b.shape[0], 1), _TAU_Y)
        wall = _gate(w, b, mu_reg, tau_y, cement_ever, cement_local)
        diffs = np.abs(np.diff(wall[:, 0]))
        assert np.all(diffs < 0.5), f"相邻元 wall 差出现 ≥0.5 跳变: max={diffs.max():.4f}"

    def test_threshold_straddle_pair_is_gradual(self):
        """最小穿越对：相邻两元 τw 恰跨阈值（9.5 / 8.9），wall 差必须 < 0.5。

        旧二值门给出 wall(9.5)=0、wall(8.9)=1（差 = 1.0）→ 旧代码必红；
        连续门给出 1−9.5/9.2→0（clip）与 1−8.9/9.2≈0.0326（差 ≈ 0.033）。
        """
        w, b, mu_reg, _ = _single_column_b([9.5, 8.9])
        cement_ever, cement_local = _cement_fields(3)
        tau_y = np.full((3, 1), _TAU_Y)
        wall = _gate(w, b, mu_reg, tau_y, cement_ever, cement_local)
        assert abs(wall[1, 0] - wall[2, 0]) < 0.5, (
            f"阈值两侧相邻元 wall 差 = {abs(wall[1, 0] - wall[2, 0]):.4f}，仍为二值悬崖")


# ---------------------------------------------------------------------------
# 2. 极限行为
# ---------------------------------------------------------------------------
class TestLimits:
    def test_vanishing_tauw_gives_full_freeze(self):
        """τw→0（窄缝极限 b→0）⇒ wall→1（全冻）。

        流动列内 τw_extrap = τw_ref·(b/b_ref) > 0 恒成立，τw→0 只能经
        b→0 极限逼近：b=1e-10 元的外推壁剪 ≈ 4e-8 ≪ 阈值 → wall ≈ 1。
        """
        b = np.array([[0.04], [0.02], [1e-10]])
        ny = b.shape[0]
        w = np.zeros((ny, 1))
        w[0, 0] = 1.0
        mu_reg = np.full((ny, 1), 0.1)
        tau_y = np.full((ny, 1), _TAU_Y)
        cement_ever, cement_local = _cement_fields(ny)
        wall = _gate(w, b, mu_reg, tau_y, cement_ever, cement_local)
        assert wall[2, 0] >= 0.999999, f"τw→0 极限未全冻: wall={wall[2, 0]}"

    def test_tauw_far_above_threshold_gives_no_freeze(self):
        """τw ≫ f·τy ⇒ wall = 0（clip 下界，不出现负值）。"""
        b = np.array([[0.04], [0.030], [0.020]])
        ny = b.shape[0]
        # 参考元高剪切（w=10）→ τw_ref = 150 ≫ 9.2，全列外推壁剪均超阈值
        w = np.zeros((ny, 1))
        w[0, 0] = 10.0
        mu_reg = np.full((ny, 1), 0.1)
        tau_y = np.full((ny, 1), _TAU_Y)
        cement_ever, cement_local = _cement_fields(ny)
        wall = _gate(w, b, mu_reg, tau_y, cement_ever, cement_local)
        assert np.all(wall == 0.0), f"高剪切列不应冻结: {wall[:, 0]}"


# ---------------------------------------------------------------------------
# 3. R3 不变量①：参考元恒 wall=0（即使其外推壁剪低于阈值）
# ---------------------------------------------------------------------------
class TestRefCellInvariant:
    def test_reference_cell_never_frozen_even_below_threshold(self):
        """整列慢速流动（各元外推壁剪均低于阈值）：参考元仍恒 wall=0。

        参考元在定义上正在流动，τw 判据不能冻结参考元自身（R3 裁定），
        否则 w>0 的列失去全部流动通道。
        """
        b = np.full((4, 1), 0.02)
        w = np.full((4, 1), 0.01)  # 均匀小速度 → τw_ref ≈ 0.3 ≪ 9.2
        w[0, 0] = 0.01
        mu_reg = np.full((4, 1), 0.1)
        tau_y = np.full((4, 1), _TAU_Y)
        cement_ever, cement_local = _cement_fields(4)
        wall = _gate(w, b, mu_reg, tau_y, cement_ever, cement_local)
        assert wall[0, 0] == 0.0, f"参考元被冻结: wall={wall[0, 0]}"


# ---------------------------------------------------------------------------
# 4. R3 不变量②：整列无流动且水泥已到 → 整列 wall=1
# ---------------------------------------------------------------------------
class TestColumnFreezeInvariant:
    def test_no_flow_column_with_cement_fully_frozen(self):
        """整列 w=0（has_flow=False）且水泥已到 → 整列 wall=1（恒不变量）。"""
        b = np.full((3, 2), 0.02)
        w = np.zeros((3, 2))
        w[0, 1] = 1.0  # 仅第 1 列有流动；第 0 列整列静止
        mu_reg = np.full((3, 2), 0.1)
        tau_y = np.full((3, 2), _TAU_Y)
        cement_ever = np.ones((3, 2))
        cement_local = np.full((3, 2), 0.9)
        wall = _gate(w, b, mu_reg, tau_y, cement_ever, cement_local)
        assert np.all(wall[:, 0] == 1.0), f"无流动列未整列冻结: {wall[:, 0]}"

    def test_no_flow_column_without_cement_not_frozen(self):
        """整列 w=0 但前锋未到（cement_ever=0）→ 不冻结（前沿可推进）。"""
        b = np.full((3, 2), 0.02)
        w = np.zeros((3, 2))
        w[0, 1] = 1.0
        mu_reg = np.full((3, 2), 0.1)
        tau_y = np.full((3, 2), _TAU_Y)
        cement_ever = np.zeros((3, 2))
        cement_local = np.zeros((3, 2))
        wall = _gate(w, b, mu_reg, tau_y, cement_ever, cement_local)
        assert np.all(wall[:, 0] == 0.0), f"前锋未到列被冻结: {wall[:, 0]}"


# ---------------------------------------------------------------------------
# 5. τy=0 零除防护语义：无屈服应力 ⇒ 无冻结（wall=0），且无 NaN/inf
# ---------------------------------------------------------------------------
class TestZeroYieldGuard:
    def test_zero_yield_stress_column_no_freeze_and_finite(self):
        """τy=0 的流动列：wall 恒 0（无屈服应力即无停流区）且输出有限。

        连续式 wall = 1 − τw/(f·τy) 在 f·τy=0 处为 0 除——语义钉定为
        wall=0（controller 裁定：τy=0 ⇒ 无屈服应力 ⇒ 无冻结），并借此
        消除 0/0 → NaN 对整场的污染路径。
        """
        b = np.full((4, 1), 0.02)
        w = np.zeros((4, 1))
        w[0, 0] = 1.0  # 列内仅参考元流动 → G>0 → 全列 τw_extrap>0
        mu_reg = np.full((4, 1), 0.1)
        tau_y = np.zeros((4, 1))  # 全列零屈服（牛顿/幂律流体混合）
        cement_ever, cement_local = _cement_fields(4)
        wall = _gate(w, b, mu_reg, tau_y, cement_ever, cement_local)
        assert np.all(np.isfinite(wall)), "τy=0 列出现非有限 wall（零除未防护）"
        assert np.all(wall == 0.0), f"τy=0 列不应冻结: {wall[:, 0]}"

    def test_mixed_yield_field_zero_cells_unfrozen(self):
        """混合 τy 场：τy=0 元 wall=0，τy>0 元按连续式取 (0,1) 中间值。

        row2：τw_extrap = τw_ref·(0.02/0.04) = 7.5 < 9.2 → 连续值
        1−7.5/9.2 ≈ 0.1848 ∈ (0,1)（旧二值此处 wall=1.0 → 本断言旧代码必红）。
        """
        b = np.array([[0.04], [0.02], [0.02]])
        w = np.zeros((3, 1))
        w[0, 0] = 1.0  # τw_ref = 0.1·6·1/0.04 = 15 → row1/2 外推壁剪 = 7.5
        mu_reg = np.full((3, 1), 0.1)
        tau_y = np.array([[8.0], [0.0], [8.0]])  # row1 零屈服
        cement_ever, cement_local = _cement_fields(3)
        wall = _gate(w, b, mu_reg, tau_y, cement_ever, cement_local)
        assert np.all(np.isfinite(wall))
        assert wall[1, 0] == 0.0, f"τy=0 元不应冻结: wall={wall[1, 0]}"
        assert 0.0 < wall[2, 0] < 1.0, f"τy>0 元应为连续中间值: wall={wall[2, 0]}"


if __name__ == "__main__":
    import sys

    sys.exit(pytest.main([__file__, "-v"]))
