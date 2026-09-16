"""Task 4.5：``gap_solver`` 批量向量化（批次轴）。

阶段 A 的 ``HBClosure`` 逐格调用 ``solve_fixed_G`` ⇒ 呼101 真实场（10000 点）单次
``mobility()`` 82–116 s（Task 4 实测）。本条把 Uzawa 的内层循环整体 numpy 化到
``(n_points, n_cells)``，各格点仍**互不耦合**（各自 ``c̄_i`` ⇒ 各自的界面掩码、
各自的 ``κ̃/n/τ̃_Y`` 逐格数组）。

锚点（**不重跑同一实现函数**）：

① **一致性**：批量 vs 逐点标量（同参同 ``ny``，含共线/非共线、``τ_Y>0``、
   ``H≠1``、端点 ``c̄=0/1``）⇒ ``rel ≤ 1e-12``（实测见 ``task-4-report.md`` §10）。
   标量路径是**冻结**的既有实现 ⇒ 这是真正的交叉核对。
② **形状与语义契约**：``u`` 在全体轴向时 ``(n_points, ny)``、存在非共线点时为
   ``(n_points, 2, ny)``；``y``/``G``/``iters``/``converged``/``undefined`` 口径。
③ **退化格点掩码**：全场未屈服 / ``G=Gb=0`` 的格点**不再抛异常中断整批**，
   而是进 ``undefined`` 掩码（⇒ ``HBClosure`` 的 R-T1-6 地板可逐格施加，
   见 ``test_hb_closure_tabulated.py`` 与 ``task-4-report.md`` §10）。
④ **失败路径**：``max_iter`` 用尽 ⇒ ``RuntimeError``（不静默，与标量同纪律）；
   形状非法 ⇒ ``ValueError``。

⚠️ 逐位/rel 对比只在同一 numpy 版本内成立（本仓基线 `shenjingwangluo`
py 3.13.7 / numpy 2.3.3）。批量路径的 ``‖·‖``/``dot`` 用逐元素式替代标量路径的
``np.linalg.norm``/``np.dot``（BLAS ⇒ 末位差 ~1 ulp，**轴向时逐位相同**，实测
20000/20000），故契约取 rel ≤ 1e-12 而非逐位。

References
----------
Bararpour & Frigaard (2025), *JFM* **1022**, A15：附录 A.2.1/A.2.3。
Zhang & Frigaard (2022), *JFM* **947**, A32：(4.21a)/(4.21b)（牛顿极限口径）。
"""

import numpy as np
import pytest

from cemdisp.models2d.gap_solver import (
    solve_fixed_G,
    solve_fixed_G_batch,
    closure_integrals_batch,
)

# 代表性样本（全部可解；退化格点单列见 _DEGENERATE）
# 列：c̄, (n₁,n₂), (κ₁,κ₂), (τ_Y,1, τ_Y,2), G, Gb, H
# _VALID_SAMPLE：**全轴向**（G[1]=Gb[1]=0）⇒ 批量 u 形状 (n_points, ny)
_VALID_SAMPLE = [
    (0.0, (1.0, 1.0), (2.0, 2.0), (0.0, 0.0), (1.0, 0.0), (0.0, 0.0), 1.0),      # 单流体牛顿端点
    (1.0, (1.0, 1.0), (2.0, 2.0), (0.0, 0.0), (1.0, 0.0), (0.0, 0.0), 1.0),
    (0.5, (1.0, 1.0), (3.0, 1.0), (0.0, 0.0), (2.0, 0.0), (0.0, 0.0), 1.0),      # 两层牛顿
    (0.5, (1.0, 1.0), (2.0, 2.0), (0.0, 0.0), (1.0, 0.0), (0.0, 0.0), 3.0),      # H≠1
    (0.0, (1.0, 1.0), (2.0, 2.0), (0.4, 0.4), (1.0, 0.0), (0.0, 0.0), 1.0),      # Bingham 单流体
    (0.4, (0.7, 0.8), (1.4, 0.9), (0.0, 0.0), (500.0, 0.0), (0.0, 0.0), 0.008),  # 纯幂律（量纲）
    (0.45, (0.7, 0.8), (1.4, 0.9), (2.0, 0.5), (500.0, 0.0), (0.0, 0.0), 0.008),  # HB 量纲口径
    (0.45, (0.7, 0.8), (1.4, 0.9), (2.0, 0.5), (3.755, 0.0), (0.0, 0.0), 1.0),
    (0.05, (0.6, 0.6), (1.0, 1.0), (0.2, 0.2), (0.8, 0.0), (0.0, 0.0), 1.0),     # 薄中线带
    (0.95, (0.6, 0.6), (1.0, 1.0), (0.2, 0.2), (0.8, 0.0), (0.0, 0.0), 1.0),     # 薄壁面带
    (0.45, (0.7, 0.8), (1.4, 0.9), (2.0, 0.5), (3.755, 0.0), (0.3, 0.0), 1.0),  # Gb≠0（共线）
    (0.5, (1.0, 1.0), (1.0, 1.0), (0.1, 0.1), (1.0, 0.0), (0.3, 0.0), 1.0),      # Gb≠0（共线）
    (0.3, (0.6, 0.6), (1.0, 1.0), (0.1, 0.1), (0.5, 0.0), (0.0, 0.0), 2.0),      # 中等 H
]
# 非共线样本（存在非共线点 ⇒ 批量 u 形状 (n_points, 2, ny)）
# ⚠️ 判据是**归约后**是否轴对齐：``Gb=0`` 时任意 ``G`` 都共线（旋转到 G 方向 ⇒ 轴对齐
# ⇒ u 仍是 1 维），故"G 是向量"本身不构成非共线；必须 ``Gb≠0 且不平行于 G``。
_NONCOLLINEAR_SAMPLE = [
    (0.5, (1.0, 1.0), (1.0, 1.0), (0.1, 0.1), (1.0, 0.4), (0.2, 0.9), 1.0),
    (0.45, (0.7, 0.8), (1.4, 0.9), (2.0, 0.5), (500.0, 150.0), (100.0, 400.0), 0.008),
    (0.35, (0.6, 0.8), (1.2, 1.1), (0.3, 0.15), (2.0, 1.0), (0.4, 0.6), 2.0),    # 非共线 + H≠1
]
# 退化格点：场内**合法**出现的两类（不得中断整批）
_DEGENERATE_SAMPLE = [
    (0.5, (1.0, 1.0), (2.0, 2.0), (5.0, 5.0), (1.0e-3, 0.0), (0.0, 0.0), 0.01),  # 全场未屈服
    (0.5, (0.7, 0.8), (1.4, 0.9), (2.0, 0.5), (0.0, 0.0), (0.0, 0.0), 0.008),    # G=Gb=0
]


def _stack(sample, field):
    """把样本列堆成批量输入（返回 c, n, kappa, tau_y, G, Gb, H）。"""
    idx = dict(c_bar=0, n=1, kappa=2, tau_y=3, G=4, Gb=5, H=6)
    j = idx[field]
    return [row[j] for row in sample]


def _scalar_row(row, ny=201):
    return solve_fixed_G(c_bar=row[0], n=row[1], kappa=row[2], tau_y=row[3],
                         G=row[4], Gb=row[5], H=row[6], ny=ny)


# --------------------------------------------------------------------------- #
# ① 一致性：批量 vs 逐点标量（rel ≤ 1e-12）
# --------------------------------------------------------------------------- #
def test_batch_matches_scalar_on_representative_sample():
    """代表性样本（共线/非共线、τ_Y>0、H≠1、端点 c̄=0/1）批量 vs 标量 rel ≤ 1e-12。"""
    ny = 201
    b = solve_fixed_G_batch(
        _stack(_VALID_SAMPLE, "c_bar"), _stack(_VALID_SAMPLE, "n"),
        _stack(_VALID_SAMPLE, "kappa"), _stack(_VALID_SAMPLE, "tau_y"),
        _stack(_VALID_SAMPLE, "G"), Gb=_stack(_VALID_SAMPLE, "Gb"),
        H=_stack(_VALID_SAMPLE, "H"), ny=ny,
    )
    rows = [_scalar_row(r, ny=ny) for r in _VALID_SAMPLE]
    for field, got in (("I1", b.I1), ("I2", b.I2), ("q0", b.q0)):
        want = np.array([getattr(r, field) for r in rows])
        assert np.allclose(got, want, rtol=1e-12, atol=1e-300), field
    # u：轴向点逐位（见模块 docstring），全体按 rel ≤ 1e-12 钉
    axial = [r[4][1] == 0.0 and r[5][1] == 0.0 for r in _VALID_SAMPLE]
    assert np.all(axial)  # 本条样本全轴向 ⇒ u 为 (n, ny)；非共线形状另测
    want_u = np.stack([np.asarray(r.u) for r in rows])
    assert np.allclose(b.u, want_u, rtol=1e-12, atol=1e-300)
    # iters 与标量逐步一致（同轨迹）
    assert np.array_equal(b.iters, np.array([r.iters for r in rows]))


def test_batch_matches_scalar_including_non_collinear_u_shape():
    """非共线样本：``u`` 形状转 ``(n_points, 2, ny)``，且与标量逐点一致。"""
    sample = _NONCOLLINEAR_SAMPLE
    assert len(sample) >= 3
    ny = 101
    b = solve_fixed_G_batch(
        _stack(sample, "c_bar"), _stack(sample, "n"), _stack(sample, "kappa"),
        _stack(sample, "tau_y"), _stack(sample, "G"), Gb=_stack(sample, "Gb"),
        H=_stack(sample, "H"), ny=ny,
    )
    assert b.u.shape == (len(sample), 2, ny)
    rows = [_scalar_row(r, ny=ny) for r in sample]
    want_u = np.stack([np.asarray(r.u) for r in rows])          # 标量：非共线 ⇒ (2, ny)
    assert want_u.shape == (len(sample), 2, ny)
    assert np.allclose(b.u, want_u, rtol=1e-12, atol=1e-300)
    assert np.allclose(b.I1, [r.I1 for r in rows], rtol=1e-12, atol=1e-300)
    # 壁面无滑移：``u[..., -1]`` 两个分量都为 0（非共线时形状 (n_points,2)）
    assert b.u[..., -1] == pytest.approx(np.zeros((len(sample), 2)), abs=1e-30)


def test_batch_axial_u_wall_index_and_u_bar():
    """``u`` 的末点语义（``u[..., -1]``）与 ``u_bar`` 口径（trapezoid）逐点复核。"""
    ny = 41
    b = solve_fixed_G_batch(
        _stack(_VALID_SAMPLE, "c_bar"), _stack(_VALID_SAMPLE, "n"),
        _stack(_VALID_SAMPLE, "kappa"), _stack(_VALID_SAMPLE, "tau_y"),
        _stack(_VALID_SAMPLE, "G"), Gb=_stack(_VALID_SAMPLE, "Gb"),
        H=_stack(_VALID_SAMPLE, "H"), ny=ny,
    )
    assert b.y.shape == (ny,) and b.y[0] == 0.0 and b.y[-1] == 1.0
    assert np.array_equal(b.u[..., -1], np.zeros(len(_VALID_SAMPLE)))
    assert np.allclose(b.u_bar, np.trapezoid(b.u, b.y, axis=-1), rtol=0, atol=0)


def test_batch_closure_only_entry_matches_scalar_closure():
    """``closure_integrals_batch``（仅闭包量，无 Uzawa）与标量 ``solve_fixed_G`` 的 I₁/I₂/q₀ 一致。

    闭包量由**解析应力场**求积（与 Uzawa 解无关，模块 docstring）⇒ 两条路径必须同值。
    """
    ny = 201
    sample = _VALID_SAMPLE + _NONCOLLINEAR_SAMPLE
    c = closure_integrals_batch(
        _stack(sample, "c_bar"), _stack(sample, "n"),
        _stack(sample, "kappa"), _stack(sample, "tau_y"),
        _stack(sample, "G"), Gb=_stack(sample, "Gb"),
        H=_stack(sample, "H"),
    )
    rows = [_scalar_row(r, ny=ny) for r in sample]
    assert np.allclose(c.I1, [r.I1 for r in rows], rtol=1e-12, atol=1e-300)
    assert np.allclose(c.I2, [r.I2 for r in rows], rtol=1e-12, atol=1e-300)
    assert np.allclose(c.q0, [r.q0 for r in rows], rtol=1e-12, atol=1e-300)
    assert not np.any(c.undefined)


def test_batch_broadcasts_scalar_pairs_and_common_vector():
    """``n/kappa/tau_y`` 允许 ``(2,)`` 广播；``G/Gb`` 允许公共 ``(2,)``。"""
    npts = 4
    c = np.array([0.2, 0.4, 0.6, 0.8])
    b = solve_fixed_G_batch(c, (0.7, 0.8), (1.0, 1.0), (0.1, 0.1),
                            (1.0, 0.0), Gb=(0.0, 0.0), H=1.0, ny=41)
    assert b.I1.shape == (npts,)
    rows = [solve_fixed_G(c_bar=float(ci), n=(0.7, 0.8), kappa=(1.0, 1.0),
                          tau_y=(0.1, 0.1), G=(1.0, 0.0), Gb=(0.0, 0.0), H=1.0, ny=41)
            for ci in c]
    assert np.allclose(b.I1, [r.I1 for r in rows], rtol=1e-12, atol=1e-300)
    # (n_points,) 形状的 G 按第 1 分量口径（标量场）
    b2 = solve_fixed_G_batch(c, (0.7, 0.8), (1.0, 1.0), (0.1, 0.1),
                             np.full(npts, 1.0), H=1.0, ny=41)
    assert np.array_equal(b2.I1, b.I1)


# --------------------------------------------------------------------------- #
# ③ 退化格点掩码（不中断整批）
# --------------------------------------------------------------------------- #
def test_batch_masks_degenerate_points_without_aborting():
    """未屈服 / ``G=Gb=0`` 的格点进 ``undefined`` 掩码；其余格点结果与标量一致。"""
    sample = list(_VALID_SAMPLE[:4]) + list(_DEGENERATE_SAMPLE) + list(_VALID_SAMPLE[4:8])
    ny = 201
    b = solve_fixed_G_batch(
        _stack(sample, "c_bar"), _stack(sample, "n"), _stack(sample, "kappa"),
        _stack(sample, "tau_y"), _stack(sample, "G"), Gb=_stack(sample, "Gb"),
        H=_stack(sample, "H"), ny=ny,
    )
    idx_deg = [4, 5]                                  # _DEGENERATE_SAMPLE 落入的位置
    assert np.array_equal(b.undefined, np.isin(np.arange(len(sample)), idx_deg))
    assert np.all(np.isnan(b.I1[idx_deg])) and np.all(np.isnan(b.u[idx_deg]))
    assert not np.any(b.converged[idx_deg])            # 退化 ⇒ 未收敛标记
    for i, row in enumerate(sample):
        if i in idx_deg:
            with pytest.raises(ValueError):
                _scalar_row(row, ny=ny)                # 标量路径对同一点抛错（交叉核对）
            continue
        ref = _scalar_row(row, ny=ny)
        assert np.isclose(b.I1[i], ref.I1, rtol=1e-12)
        assert b.iters[i] == ref.iters
    assert np.all(b.converged[np.setdiff1d(np.arange(len(sample)), idx_deg)])


def test_batch_all_degenerate_returns_mask_not_exception():
    """整批退化 ⇒ 全掩码（不抛错）——与标量"逐点抛错"的差异在 docstring 里写明。"""
    b = solve_fixed_G_batch(
        _stack(_DEGENERATE_SAMPLE, "c_bar"), _stack(_DEGENERATE_SAMPLE, "n"),
        _stack(_DEGENERATE_SAMPLE, "kappa"), _stack(_DEGENERATE_SAMPLE, "tau_y"),
        _stack(_DEGENERATE_SAMPLE, "G"), Gb=_stack(_DEGENERATE_SAMPLE, "Gb"),
        H=_stack(_DEGENERATE_SAMPLE, "H"), ny=101,
    )
    assert np.all(b.undefined) and not np.any(b.converged)
    assert np.all(np.isnan(b.I1)) and np.all(np.isnan(b.I2)) and np.all(np.isnan(b.q0))


# --------------------------------------------------------------------------- #
# ④ 失败路径与形状校验
# --------------------------------------------------------------------------- #
def test_batch_raises_on_max_iter_exhausted():
    """``max_iter`` 用尽仍不收敛 ⇒ ``RuntimeError``（与标量同纪律，不静默）。"""
    with pytest.raises(RuntimeError, match="未收敛"):
        solve_fixed_G_batch(np.full(3, 0.45), np.tile((0.7, 0.8), (3, 1)),
                            np.tile((1.4, 0.9), (3, 1)), np.tile((2.0, 0.5), (3, 1)),
                            np.full((3, 2), 3.755), H=1.0, ny=201, max_iter=3)


def test_batch_validates_inputs():
    """形状/取值非法 ⇒ 显式 ``ValueError``。"""
    c = np.full(3, 0.5)
    base = dict(n=(0.7, 0.8), kappa=(1.0, 1.0), tau_y=(0.1, 0.1), G=(1.0, 0.0), ny=41)
    solve_fixed_G_batch(c, **base)                                   # 合法基准
    for bad in ({"n": np.tile((0.7, 0.8), (4, 1))},        # 批次长度与 c_bar 不一致
                {"c_bar": np.full(3, 1.5)},
                {"c_bar": np.full(3, -0.1)}, {"H": np.zeros(3)},
                {"G": np.full((3, 3), 1.0)}, {"G": np.full(4, 1.0)},
                {"n": np.tile((0.0, 0.8), (3, 1))},
                {"kappa": np.tile((1.0, -1.0), (3, 1))},
                {"tau_y": np.tile((-1.0, 0.5), (3, 1))},
                {"ny": 2}, {"r": 0.0}):
        with pytest.raises(ValueError):
            solve_fixed_G_batch(**{**dict(c_bar=c, **base), **bad})