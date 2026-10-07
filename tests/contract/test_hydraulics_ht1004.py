# -*- coding: utf-8 -*-
"""Phase 3.1（P-2，2026-10-07）契约测试：`cemdisp/transport1d/hydraulics.py` 对 MATLAB 靶逐点复现。

靶权威 = `results/_probe_matlab靶_2026-10-07/sandbox/out_*.csv`（2026-10-07 沙箱 -batch 运行，
wrapper_p30.m 捕获，原脚本零改动；用户裁定 2026-10-07 为权威靶）。判据 = spec §1.4 表：
泵压相对 ≤1e-9、井底 ECD/全场 ECD/走廊/回压绝对 ≤1e-9、conflict 逐位、时间轴逐位。

浮点次序归因（spec §5-1；禁调参、禁放宽容差）
----------------------------------------------
- 实测偏差 1e-13~1e-15 量级，来源两类、均在判据内：
  ① 靶 CSV 由 MATLAB `writematrix` 以 **15 位有效数字**落盘，大数列入读即含
     ~1e-13（对 ~147 MPa 列）量化误差；
  ② `cos/pow` 等超越函数在 MATLAB libm 与 CPython libm 间允许 ±1 ulp 差异，
     沿不动点迭代/递推链传播后仍为 1e-13 量级（地图 §6 纪律 2 允许 ULP 级分歧）。
- 求和顺序镜像（`_sum_seq`/逐字表达式次序）已消除结构性偏差；若某列未来超出
  1e-9，先按上述两因归因再报主会话裁定，**测试内禁放宽判据**。

红线哨兵（spec §1.2 —— 修 bug = 失配 = 本文件红）
--------------------------------------------------
`:1153` tau 混用 / `:1464` safety_margin 不进走廊 / `:733` +1200 L / `n_time=floor=189`
——各由对应列的逐点复现钉住（改动即失配）；另附 `test_s25_events_and_frozen_bp` 钉
冻结事件时刻与 5.117 MPa。
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pytest

BRANCH_ROOT = Path(__file__).resolve().parents[2]
if str(BRANCH_ROOT) not in sys.path:
    sys.path.insert(0, str(BRANCH_ROOT))

from cemdisp.transport1d.hydraulics import run_ht1004_target  # noqa: E402

SBX = BRANCH_ROOT / "results" / "_probe_matlab靶_2026-10-07" / "sandbox"
TOL_ABS = 1e-9  # spec §1.4
TOL_REL = 1e-9
N_TIME = 189  # 摘要级标量：MATLAB :132 floor(189.059)=189（可字面，spec/任务允许）
N_SEG = 333


@pytest.fixture(scope="module")
def res():
    """全链只跑一次（~2 s），各测试共享结果。"""
    assert SBX.is_dir(), f"靶目录缺失：{SBX}"
    return run_ht1004_target(SBX)


def _load(name: str) -> np.ndarray:
    return np.loadtxt(str(SBX / name), delimiter=",")


def _max_abs(mine, tgt) -> float:
    d = float(np.max(np.abs(np.asarray(mine, dtype=float) - np.asarray(tgt, dtype=float))))
    assert np.all(np.isfinite(mine)), "出现 NaN/Inf（关5）"
    return d


def _max_rel(mine, tgt) -> float:
    mine = np.asarray(mine, dtype=float)
    tgt = np.asarray(tgt, dtype=float)
    assert np.all(np.isfinite(mine)) and np.all(np.isfinite(tgt))
    return float(np.max(np.abs(mine - tgt) / np.abs(tgt)))


# --------------------------------------------------------------------------
# 0. 靶件在场与形状
# --------------------------------------------------------------------------

def test_target_files_and_shapes(res):
    rm = res["out_result_matrix_volume_pressure"]
    assert rm.shape == (N_TIME, 5)
    assert res["out_ecd_casing_full"].shape == (N_TIME, N_SEG)
    for name in ("out_result_matrix_volume_pressure.csv", "out_pump_pressure_comparison.csv",
                 "out_bottom_ecd_comparison.csv", "out_four_point_backpressure.csv",
                 "out_four_point_ecd_new.csv", "out_pump_pressure_surface.csv",
                 "out_annuli_bottom_pressure.csv", "out_ecd_casing_full.csv",
                 "out_window_backpressure_ctrl.csv", "out_summary.txt"):
        assert (SBX / name).is_file(), name


def test_time_axis_bitwise(res):
    """时间轴逐位（1..189；MATLAB :1480 `(1:n)'*dt`，dt=1）。"""
    t = np.arange(1, N_TIME + 1, dtype=float)
    for name, col in (("out_pump_pressure_surface.csv", 0),
                      ("out_annuli_bottom_pressure.csv", 0),
                      ("out_result_matrix_volume_pressure.csv", 0),
                      ("out_pump_pressure_comparison.csv", 0),
                      ("out_bottom_ecd_comparison.csv", 0),
                      ("out_four_point_backpressure.csv", 0),
                      ("out_four_point_ecd_new.csv", 0),
                      ("out_window_backpressure_ctrl.csv", 0)):
        assert np.array_equal(_load(name)[:, col], t), name


# --------------------------------------------------------------------------
# 1. 泵压时程（相对 ≤1e-9）
# --------------------------------------------------------------------------

def test_pump_pressure_series_relative(res):
    tgt = _load("out_pump_pressure_surface.csv")
    mine = res["out_pump_pressure_surface"]
    r = _max_rel(mine[:, 1], tgt[:, 1])
    # 归因：15 位落盘量化 + libm pow/cos ±1ulp 经递推链传播 ≈1.6e-14（见模块 docstring）
    assert r <= TOL_REL, f"泵压 max|rel|={r:.3e}"


# --------------------------------------------------------------------------
# 2. 井底 ECD 时程（绝对 ≤1e-9）
# --------------------------------------------------------------------------

def test_bottom_ecd_series_absolute(res):
    tgt = _load("out_bottom_ecd_comparison.csv")
    mine = res["out_bottom_ecd_comparison"]
    assert _max_abs(mine[:, 2], tgt[:, 2]) <= TOL_ABS  # calculated_bottom_ECD（环空井底列）
    # 管内井底 ECD 列（Result_Matrix col3 ≡ out_ecd_casing_full 末列）同源复钉
    full = _load("out_ecd_casing_full.csv")
    assert _max_abs(res["out_result_matrix_volume_pressure"][:, 2], full[:, -1]) <= TOL_ABS


# --------------------------------------------------------------------------
# 3. ECD 全场 189×333（绝对 ≤1e-9）
# --------------------------------------------------------------------------

def test_ecd_casing_full_field(res):
    tgt = _load("out_ecd_casing_full.csv")
    assert _max_abs(res["out_ecd_casing_full"], tgt) <= TOL_ABS


# --------------------------------------------------------------------------
# 4. 四点走廊（conflict 逐位；BP ≤1e-9）
# --------------------------------------------------------------------------

def test_four_point_corridor(res):
    tgt = _load("out_four_point_backpressure.csv")
    mine = res["out_four_point_backpressure"]
    assert np.array_equal(mine[:, 4], tgt[:, 4])  # conflict 逐位（全 1 如实：189/189）
    for j in (1, 2, 3):
        assert _max_abs(mine[:, j], tgt[:, j]) <= TOL_ABS
    # 关3：下界≤上界处不得报冲突；靶全程冲突如实（189/189）
    ok = mine[:, 2] <= mine[:, 3] + 1e-6
    assert not np.any(mine[ok, 4] == 1)
    assert int(mine[:, 4].sum()) == 189


def test_four_point_ecd_new(res):
    tgt = _load("out_four_point_ecd_new.csv")
    assert _max_abs(res["out_four_point_ecd_new"], tgt) <= TOL_ABS


# --------------------------------------------------------------------------
# 5. 井底窗持稳回压 + 状态机事件
# --------------------------------------------------------------------------

def test_bottom_window_columns(res):
    tgt = _load("out_window_backpressure_ctrl.csv")
    mine = res["out_window_backpressure_ctrl"]
    assert _max_abs(mine[:, 1], tgt[:, 1]) <= TOL_ABS  # ECD_base_bottom_win
    assert _max_abs(mine[:, 2], tgt[:, 2]) <= TOL_ABS  # ECD_bottom_ctrl_win（持稳 2.025）
    assert _max_abs(mine[:, 3], tgt[:, 3]) <= TOL_ABS  # ECD_shoe_ctrl_win
    assert _max_abs(mine[:, 4], tgt[:, 4]) <= TOL_ABS  # P_wh_ctrl_win


def test_bottom_window_bp_apply(res):
    """回压施加值 ≤1e-9：BP_apply 非 CSV 直接列，经靶列重构
    bp_target ≈ (ECD_ctrl − ECD_base)·(0.00981·TVD_cum(end))；ECD 15 位量化的重构误差
    ~1e-13（实测 3.7e-13），在判据内。"""
    tgt = _load("out_window_backpressure_ctrl.csv")
    tvd_end = float(np.max(res["tvd_cum_m1b"]))
    den = 0.00981 * tvd_end
    bp_tgt = (tgt[:, 2] - tgt[:, 1]) * den
    assert _max_abs(res["bp_apply_win_MPa"], bp_tgt) <= TOL_ABS


def test_s25_events_and_frozen_bp(res):
    """状态机事件时刻逐位（交点 189 min / 跌破段 1–189）+ 冻结回压 5.117 MPa。

    5.117 = MATLAB fprintf %.3f 显示值（3.0_靶复现判定.md/diary）——摘要级标量，
    按显示精度 ±5e-4 断言；事件时刻/段号逐位。
    """
    s = res["summary"]
    assert s["freeze_segments_1b"] == [(1, N_TIME)]  # 跌破段 1 = 1–189 min（逐位）
    assert s["bp_apply_at_last_min_MPa"] == pytest.approx(5.117, abs=5e-4)
    # 交点后冻结 ⇒ 全列跌破段内 ECD_bottom_ctrl ≡ 2.025
    w = res["out_window_backpressure_ctrl"]
    assert np.max(np.abs(w[:, 2] - 2.025)) <= TOL_ABS


# --------------------------------------------------------------------------
# 6. 对比模块 reference 插值口径（unique-stable+钳位+linear）
# --------------------------------------------------------------------------

def test_reference_interp_design_columns_bitwise(res):
    """S19/S20：设计插值列 ≤1e-9（np 口径与 MATLAB interp1 的 ulp 级次序差远小于判据），
    计算列同源、误差列=计算−插值 一致性钉住 unique-stable+钳位链。"""
    tp = _load("out_pump_pressure_comparison.csv")
    mp = res["out_pump_pressure_comparison"]
    assert _max_abs(mp[:, 3], tp[:, 3]) <= TOL_ABS  # design_pressure_interp
    assert _max_abs(mp[:, 4], tp[:, 4]) <= TOL_ABS  # error
    te = _load("out_bottom_ecd_comparison.csv")
    me = res["out_bottom_ecd_comparison"]
    assert _max_abs(me[:, 3], te[:, 3]) <= TOL_ABS
    assert _max_abs(me[:, 4], te[:, 4]) <= TOL_ABS
    # 误差列内部一致性（MATLAB error = calc − interp 逐字）
    assert np.max(np.abs((me[:, 2] - me[:, 3]) - me[:, 4])) <= TOL_ABS
    assert np.max(np.abs((mp[:, 2] - mp[:, 3]) - mp[:, 4])) <= TOL_ABS


# --------------------------------------------------------------------------
# 7. 结果矩阵 / 体积链闭合（关5）/ 无 NaN
# --------------------------------------------------------------------------

def test_result_matrix_and_volume_chain(res):
    tgt = _load("out_result_matrix_volume_pressure.csv")
    mine = res["out_result_matrix_volume_pressure"]
    assert np.array_equal(mine[:, 1], tgt[:, 1]) or _max_abs(mine[:, 1], tgt[:, 1]) <= TOL_ABS
    assert _max_abs(mine[:, 4], tgt[:, 4]) <= TOL_ABS  # 环空井底压力 MPa
    assert _max_rel(mine[:, 3], tgt[:, 3]) <= TOL_REL  # 泵压列（与 col3 同源复钉）
    # 关5：体积链闭合（累计注入 vs 靶 Cumulative_Volume_m3 逐位级）
    assert _max_abs(res["cumulative_volume_m3"], tgt[:, 1]) <= TOL_ABS
    # 关5：全程无 NaN/Inf
    for k in ("out_result_matrix_volume_pressure", "out_ecd_casing_full",
              "out_four_point_backpressure", "out_four_point_ecd_new",
              "out_window_backpressure_ctrl", "out_pump_pressure_comparison",
              "out_bottom_ecd_comparison", "bp_apply_win_MPa", "cumulative_volume_m3"):
        assert np.all(np.isfinite(res[k])), k


# --------------------------------------------------------------------------
# 8. 方向性锚（关3）+ out_summary 标量
# --------------------------------------------------------------------------

def test_directionality_gate3(res):
    pump = res["out_pump_pressure_surface"][:, 1]
    q = res["pump_Lmin_1b"][1:]
    corr = np.corrcoef(pump, q)[0, 1]
    assert corr > 0.5  # 泵压随排量同向（阶跃接口期允许局部回勾）
    st = res["annulus_static_MPa"][1:, 1:]  # 静液柱随深度非降（ρ、Δz>0）
    assert np.all(np.diff(st, axis=1) >= 0.0)


def test_summary_scalars_against_out_summary(res):
    txt = (SBX / "out_summary.txt").read_text(encoding="utf-8")
    kv = dict(re.findall(r"^(\w+) = ([0-9.eE+-]+)$", txt, flags=re.M))
    s = res["summary"]
    checks = {
        "n_time": s["n_time"],  # 摘要级标量字面允许
        "n_segment": s["n_segment"],
        "TVD_bottom": s["TVD_bottom_m"],
        "max_pump_pressure_MPa": s["max_pump_MPa"],
        "min_pump_pressure_MPa": s["min_pump_MPa"],
        "max_annuli_bottom_MPa": s["max_annuli_bottom_MPa"],
        "max_bottom_ECD_g_cm3": s["max_bottom_ECD_casing_g_cm3"],
        "min_bottom_ECD_g_cm3": s["min_bottom_ECD_casing_g_cm3"],
    }
    for name, mine in checks.items():
        assert name in kv, name
        tv = float(kv[name])
        rel = abs(float(mine) - tv) / abs(tv)
        # out_summary 为 %.10g 打印（10 位有效）⇒ 断言 1e-8（打印量化）
        assert rel <= 1e-8, f"{name}: rel={rel:.3e}"
    assert re.search(r"conflict_steps = 189 / 189", txt)
    assert s["conflict_steps"] == N_TIME
