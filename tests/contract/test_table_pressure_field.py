# -*- coding: utf-8 -*-
"""Phase 3.1（2026-10-07）契约测试：`TablePressureField`（pressure_field.py 二期追加）。

覆盖（spec §1.3）
----------------
1. 协议兼容：runtime_checkable ``PressureField`` 协议 isinstance。
2. 双线性插值正确性：网格节点精确 / 行内、列内线性 / 面中点手算。
3. oob 三件套**同型 HydrostaticPressureField**：计数条件 = md<0 或 md>bottom_md_m；
   md 钳位端点延拓返回；**时间轴钳位不计数**（对位 Hydrostatic 的"t 被忽略"——本类钳位）。
4. reset_audit / oob_events 字段。
5. 构造校验：形状不符 / 网格非升序 / 长度 <2 ⇒ ValueError。
6. 冒烟数据源 = 靶 `out_annuli_bottom_pressure.csv`（(t, P_bottom) → 2 节点 md 表）：
   井底行精确复现靶列（≤1e-12）、井口行为 0、md 中点=线性。
   ⚠ 仅本阶段落地+单测，**不注入生产路径**（C-14 同型声明）。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

BRANCH_ROOT = Path(__file__).resolve().parents[2]
if str(BRANCH_ROOT) not in sys.path:
    sys.path.insert(0, str(BRANCH_ROOT))

from cemdisp.data.pressure_field import (  # noqa: E402
    PressureField,
    PressureOutOfRangeEvent,
    TablePressureField,
)

SBX = BRANCH_ROOT / "results" / "_probe_matlab靶_2026-10-07" / "sandbox"

# 手算小表：md=[0,1,2]，t=[0,10]s，列0=[0,10,20]，列1=[1,12,25]
MD = [0.0, 1.0, 2.0]
TS = [0.0, 10.0]
TAB = [[0.0, 1.0], [10.0, 12.0], [20.0, 25.0]]


def _mk(**kw) -> TablePressureField:
    return TablePressureField("t_well", MD, TS, TAB, **kw)


# ---------------------------------------------------------------- 1 协议
def test_protocol_compat():
    f = _mk()
    assert isinstance(f, PressureField)  # runtime_checkable 协议
    assert callable(f.P)


# ---------------------------------------------------------------- 2 插值正确性
def test_grid_nodes_exact():
    f = _mk()
    for i, m in enumerate(MD):
        for j, t in enumerate(TS):
            assert f.P(m, t) == pytest.approx(TAB[i][j], abs=1e-12)


def test_row_and_column_linearity():
    f = _mk()
    assert f.P(0.5, 0.0) == pytest.approx(5.0, abs=1e-12)     # md 方向线性（列 0）
    assert f.P(1.5, 10.0) == pytest.approx(18.5, abs=1e-12)   # md 方向线性（列 1）
    assert f.P(0.0, 5.0) == pytest.approx(0.5, abs=1e-12)     # t 方向线性（md=0）
    assert f.P(2.0, 5.0) == pytest.approx(22.5, abs=1e-12)    # t 方向线性（md=2）


def test_bilinear_midpoint_handcalc():
    f = _mk()
    # (md=1, t=5)：列插值 10 与 12 → 11（双线性/矩形网格一致）
    assert f.P(1.0, 5.0) == pytest.approx(11.0, abs=1e-12)


# ---------------------------------------------------------------- 3 oob / 钳位
def test_oob_md_axis_counted_and_clamped():
    f = _mk()
    assert f.P(-3.0, 5.0) == pytest.approx(0.5, abs=1e-12)   # 钳到 md=0 ⇒ t=5 行值
    assert f.P(9.0, 5.0) == pytest.approx(22.5, abs=1e-12)   # 钳到 md=2（bottom_md_m=2）
    assert f.oob_count == 2
    evs = f.oob_events
    assert isinstance(evs, tuple) and all(isinstance(e, PressureOutOfRangeEvent) for e in evs)
    assert evs[0].md_m == -3.0 and evs[0].bottom_md_m == 2.0
    assert evs[1].md_m == 9.0


def test_oob_time_axis_clamped_not_counted():
    f = _mk()
    before = f.oob_count
    assert f.P(1.0, -5.0) == pytest.approx(10.0, abs=1e-12)  # t 钳到 0 ⇒ 列 0 行值
    assert f.P(1.0, 3600.0) == pytest.approx(12.0, abs=1e-12)  # t 钳到 10 ⇒ 列 1
    assert f.oob_count == before  # 时间轴不计数（口径见类 docstring）


def test_reset_audit():
    f = _mk()
    f.P(-1.0, 0.0)
    assert f.oob_count == 1
    f.reset_audit()
    assert f.oob_count == 0 and f.oob_events == ()


def test_bottom_md_override():
    f = _mk(bottom_md_m=1.5)
    f.P(1.6, 0.0)  # 越过自定义底深
    assert f.oob_count == 1


# ---------------------------------------------------------------- 5 构造校验
def test_validation_raises():
    with pytest.raises(ValueError):
        TablePressureField("w", [0.0], [0.0, 1.0], [[1.0, 2.0]])          # md 长度 <2
    with pytest.raises(ValueError):
        TablePressureField("w", [0.0, 1.0], [0.0], [[1.0], [2.0]])        # t 长度 <2
    with pytest.raises(ValueError):
        TablePressureField("w", [0.0, 1.0], [0.0, 10.0], [[1.0]])         # 形状不符
    with pytest.raises(ValueError):
        TablePressureField("w", [1.0, 0.0], [0.0, 10.0], [[1.0, 2.0], [3.0, 4.0]])  # md 非升序
    with pytest.raises(ValueError):
        TablePressureField("w", [0.0, 1.0], [10.0, 0.0], [[1.0, 2.0], [3.0, 4.0]])  # t 非升序


# ---------------------------------------------------------------- 6 靶冒烟
def test_smoke_annuli_bottom_pressure_target():
    """用靶 `out_annuli_bottom_pressure.csv`（[time_min, P_bottom_MPa]）构 2×189 表：
    md=[0, TVD_bottom]，t=time×60；井底行应精确复现靶列，井口行恒 0，md 中点=半值。"""
    path = SBX / "out_annuli_bottom_pressure.csv"
    if not path.is_file():
        pytest.skip("靶件缺失（_probe_matlab靶_2026-10-07）")
    tgt = np.loadtxt(str(path), delimiter=",")
    time_min = tgt[:, 0]
    p_bot = tgt[:, 1]
    bottom = 7656.92  # TVD_bottom（out_summary；摘要级标量允许字面）
    f = TablePressureField("HT1-004", [0.0, bottom], list(time_min * 60.0),
                           np.vstack([np.zeros_like(p_bot), p_bot]))  # 形状 (n_md=2, n_t=189)
    for j in (0, 94, 188):
        t_s = float(time_min[j] * 60.0)
        assert f.P(bottom, t_s) == pytest.approx(float(p_bot[j]), rel=1e-12, abs=1e-12)
        assert f.P(0.0, t_s) == pytest.approx(0.0, abs=1e-12)
        assert f.P(bottom / 2.0, t_s) == pytest.approx(float(p_bot[j]) / 2.0, rel=1e-12, abs=1e-12)
    # 时间列间：井底行 = 相邻时刻线性（双线性退化一致）
    assert f.P(bottom, float((time_min[5] + time_min[6]) / 2.0 * 60.0)) == \
        pytest.approx((p_bot[5] + p_bot[6]) / 2.0, rel=1e-12, abs=1e-12)
    assert f.oob_count == 0
    f.P(bottom + 1.0, 100.0)  # md 越界 ⇒ 计数
    assert f.oob_count == 1
