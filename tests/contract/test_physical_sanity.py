# -*- coding: utf-8 -*-
"""物理合理性闸门契约测试（Task 9）。

纪律（协调者 2026-09-26 附加裁定 #4/#5）
----------------------------------------
1. 7 项判据**每一项**都必须有一条"构造反例 ⇒ 该项 FAIL"的测试；只在恒真情形上断言
   视为缺陷。
2. 条件性判据的"不适用"情形必须判为 ``通过=None``（``说明`` 以 ``未测：`` 开头），
   **不得**默认 True、**不得**改断言别的更容易通过的东西来凑"通过"。
3. 三项判据（P-3/P-4/P-6）在结果对象上不可测（``AnnulusSimulationResult`` 不导出
   速度/通量/流函数场，R155）⇒ 改判测点为**求解器级合成算例**：本文件直接调
   ``solve_stream_function`` / ``velocity_from_stream_function``（P-6 带 ``wall=``）
   做真断言，反例也由真实求解器调用产生（e=0.60 偏心算例 / 破坏 Dirichlet 的 Ψ /
   wall 未进算子），不是手搓的假场。
"""
from __future__ import annotations

import numpy as np
import pytest

from cemdisp.diagnostics.physical_sanity import (
    FLUX_REL_TOL,
    MONOTONE_TOL,
    RATE_RESPONSE_TOL_PP,
    VELOCITY_AMPLIFICATION_MAX,
    sanity_checks,
    synthetic_stream_case,
)
from cemdisp.models2d.stream_function import velocity_from_stream_function

P1 = "P-1 浓度有界"
P2 = "P-2 浓度单调（环空段，无源项平流）"
P3 = "P-3 速度上界"
P4 = "P-4 总通量守恒"
P5 = "P-5 窄边劣势"
P6 = "P-6 冻结区静止"
P7 = "P-7 排量响应方向"


class _Res:
    """最小结果对象替身（只带 sanity_checks 消费的字段）。"""

    def __init__(self, snaps, *, eta_e=None, b_num=None, final=None):
        self.cement_snapshots = tuple(np.asarray(s, dtype=float) for s in snaps)
        self.snapshot_times_s = tuple(float(i) for i in range(len(snaps)))
        self.summary = {}
        if eta_e is not None or b_num is not None:
            final_summary = {}
            if eta_e is not None:
                final_summary["全井段最终有效顶替效率"] = float(eta_e)
            if b_num is not None:
                final_summary["浮力数_b"] = float(b_num)
            self.summary["最终结果"] = final_summary
        if final is not None:
            self.cement_field = np.asarray(final, dtype=float)
        elif snaps:
            self.cement_field = np.asarray(snaps[-1], dtype=float)
        else:
            self.cement_field = np.zeros((0, 0))


def _geom(ny=4, nz=3, b=None):
    return {"b": np.ones((ny, nz)) if b is None else np.asarray(b, dtype=float),
            "H": np.full((ny, nz), 0.01),
            "phi": np.linspace(0, 1, ny), "s": np.linspace(0, 2, nz),
            "y": np.linspace(0, 0.3, ny)}


def _rows(rows):
    return {r["检查项"]: r for r in rows}


def _ratio(case):
    """合成算例的偏心放大比 max|w| / (Q/A_min)。"""
    return float(np.max(np.abs(case["w_m_s"]))) / (case["Q_m3s"] / case["A_min_m2"])


# --------------------------------------------------------------------------- #
# brief 原样的两条（P-2 键名按 R138 更名）
# --------------------------------------------------------------------------- #
def test_bounded_and_monotone_pass_on_valid_field():
    a = np.zeros((4, 3)); a[:, :1] = 1.0
    b = np.zeros((4, 3)); b[:, :2] = 1.0
    geom = {"b": np.ones((4, 3)), "H": np.full((4, 3), 0.01),
            "phi": np.linspace(0, 1, 4), "s": np.linspace(0, 2, 3),
            "y": np.linspace(0, 0.3, 4)}
    rows = {r["检查项"]: r for r in sanity_checks(_Res([a, b]), geom)}
    assert rows[P1]["通过"] is True
    assert rows[P2]["通过"] is True


def test_monotone_fails_on_nonphysical_field():
    a = np.zeros((4, 3)); a[:, :2] = 1.0
    b = np.zeros((4, 3)); b[:, :1] = 1.0        # 浓度"倒退"⇒ 非物理
    geom = {"b": np.ones((4, 3)), "H": np.full((4, 3), 0.01),
            "phi": np.linspace(0, 1, 4), "s": np.linspace(0, 2, 3),
            "y": np.linspace(0, 0.3, 4)}
    rows = {r["检查项"]: r for r in sanity_checks(_Res([a, b]), geom)}
    assert rows[P2]["通过"] is False


# --------------------------------------------------------------------------- #
# 7 条判据各一条"反例 ⇒ FAIL"
# --------------------------------------------------------------------------- #
def test_P1_bounded_fails_on_out_of_range_field():
    """反例：体积分数越界（>1 与 <0 两侧各一例）。"""
    for bad_value in (1.2, -0.05):
        good = np.zeros((4, 3)); good[:, :1] = 1.0
        bad = np.zeros((4, 3)); bad[:, :2] = 1.0
        bad[0, 1] = bad_value
        rows = _rows(sanity_checks(_Res([good, bad]), _geom()))
        assert rows[P1]["通过"] is False, f"越界值 {bad_value} 未被 P-1 抓住"


def test_P2_monotone_fails_on_column_backstep():
    """反例：单列 b 加权浓度回退 0.2（远超 1e-9 容许）。"""
    a = np.zeros((4, 3)); a[:, :1] = 1.0
    b = np.zeros((4, 3)); b[:, :1] = 0.8         # 回退 0.2
    rows = _rows(sanity_checks(_Res([a, b]), _geom()))
    assert rows[P2]["通过"] is False
    assert rows[P2]["判据"].find("1e-09") >= 0 or rows[P2]["判据"].find("1e-9") >= 0


def test_P3_velocity_bound_fails_on_strongly_eccentric_case():
    """反例：e=0.60 强偏心算例 —— 真实求解器调用得到的偏心放大超过 1.5 倍判据。"""
    case = synthetic_stream_case(standoff=1.0 - 0.60)
    rows = _rows(sanity_checks(_Res([np.zeros((4, 3))]), _geom(), stream_case=case))
    assert rows[P3]["通过"] is False
    assert _ratio(case) > VELOCITY_AMPLIFICATION_MAX


def test_P3_passes_on_design_eccentricity_case():
    """正例：设计居中度算例（须真通过，否则 P-3 是恒红噪声）。"""
    case = synthetic_stream_case()
    rows = _rows(sanity_checks(_Res([np.zeros((4, 3))]), _geom(), stream_case=case))
    assert rows[P3]["通过"] is True


def test_P3_matches_analytic_eccentric_amplification():
    """P-3 的数值放大必须等于解析式 (1+e)²/(1+1.5e²)（I₁ ∝ H³ 闭包 + 单位通量 BC）。

    这条把判据与偏心度的关系钉住：1.5 倍的适用边界是 e ≈ 0.31，不是可以随便放宽的旋钮。
    """
    for e in (0.0, 0.22, 0.30):
        case = synthetic_stream_case(standoff=1.0 - e)
        assert _ratio(case) == pytest.approx((1.0 + e) ** 2 / (1.0 + 1.5 * e ** 2), rel=2e-2)


def test_P4_flux_conservation_fails_on_broken_dirichlet_bc():
    """反例：把 φ=1 端 Dirichlet 改成逐列不同 ⇒ 列通量不再守恒（真实 Ψ + 破坏的 BC）。"""
    case = synthetic_stream_case()
    psi = np.array(case["psi"], dtype=float, copy=True)
    psi[-1, :] = 1.0 + 0.1 * np.sin(np.arange(psi.shape[1], dtype=float))
    w_unit, _ = velocity_from_stream_function(psi, case["geom"])
    rows = _rows(sanity_checks(_Res([np.zeros((4, 3))]), _geom(),
                               stream_case={**case, "w_unit": w_unit}))
    assert rows[P4]["通过"] is False


def test_P4_flux_conservation_passes_on_solver_output():
    case = synthetic_stream_case()
    rows = _rows(sanity_checks(_Res([np.zeros((4, 3))]), _geom(), stream_case=case))
    assert rows[P4]["通过"] is True


def test_P5_narrow_disadvantage_fails_when_narrow_quarter_leads():
    """反例：水泥只在窄 1/4 行 ⇒ η_窄 > η_宽，违反偏心环空顶替基本物理。"""
    ny, nz = 8, 5
    final = np.zeros((ny, nz)); final[-(ny // 4):, :] = 1.0
    res = _Res([final], eta_e=0.3, b_num=12.0, final=final)
    rows = _rows(sanity_checks(res, _geom(ny=ny, nz=nz)))
    assert rows[P5]["通过"] is False


def test_P5_passes_when_wide_quarter_leads():
    ny, nz = 8, 5
    final = np.zeros((ny, nz)); final[: ny // 4, :] = 1.0     # 宽边先到（物理方向）
    res = _Res([final], eta_e=0.3, b_num=12.0, final=final)
    rows = _rows(sanity_checks(res, _geom(ny=ny, nz=nz)))
    assert rows[P5]["通过"] is True


def test_P6_frozen_region_fails_when_wall_not_consumed():
    """反例：wall 掩码已算出但未进算子（R42 死开关类缺陷）⇒ 冻结区并不静止。"""
    case = synthetic_stream_case(freeze_rows="auto", apply_wall=False)
    rows = _rows(sanity_checks(_Res([np.zeros((4, 3))]), _geom(), frozen_case=case))
    assert rows[P6]["通过"] is False


def test_P7_rate_response_fails_when_rate_drops_efficiency():
    """反例：排量 ×1.4 反而掉 5 pp（> 1 pp 容许）。"""
    field = np.zeros((4, 3)); field[:, :2] = 1.0
    base = _Res([field], eta_e=0.90, final=field)
    rate = _Res([field], eta_e=0.85, final=field)
    rows = _rows(sanity_checks(base, _geom(), rate_result=rate))
    assert rows[P7]["通过"] is False
    assert float(RATE_RESPONSE_TOL_PP) == 1.0


# --------------------------------------------------------------------------- #
# 未测纪律：条件不满足时必须 None + "未测：<原因>"，绝不默认 True
# --------------------------------------------------------------------------- #
def test_P2_not_measured_without_snapshots():
    rows = _rows(sanity_checks(_Res([]), _geom()))
    assert rows[P2]["通过"] is None
    assert rows[P2]["说明"].startswith("未测：")


def test_P5_not_measured_for_density_unstable_well():
    """密度倒置井（浮力数 b ≤ 0）判据不适用 ⇒ 未测，不得默认 True。"""
    ny, nz = 8, 5
    final = np.zeros((ny, nz)); final[-2:, :] = 1.0
    res = _Res([final], eta_e=0.3, b_num=-3.0, final=final)
    rows = _rows(sanity_checks(res, _geom(ny=ny, nz=nz)))
    assert rows[P5]["通过"] is None
    assert "未测：" in rows[P5]["说明"]


def test_P5_not_measured_without_buoyancy_number():
    ny, nz = 8, 5
    final = np.zeros((ny, nz)); final[-2:, :] = 1.0
    rows = _rows(sanity_checks(_Res([final], final=final), _geom(ny=ny, nz=nz)))
    assert rows[P5]["通过"] is None
    assert rows[P5]["说明"].startswith("未测：")


def test_P7_not_measured_without_rate_run():
    field = np.zeros((4, 3)); field[:, :1] = 1.0
    rows = _rows(sanity_checks(_Res([field], eta_e=0.9, final=field), _geom()))
    assert rows[P7]["通过"] is None
    assert rows[P7]["说明"].startswith("未测：")


def test_P6_not_measured_when_yield_gate_is_off():
    """屈服门未启用（wall=None）⇒ 无冻结单元 ⇒ 未测，不得默认 True。"""
    case = synthetic_stream_case(freeze_rows=None)
    rows = _rows(sanity_checks(_Res([np.zeros((4, 3))]), _geom(), frozen_case=case))
    assert rows[P6]["通过"] is None
    assert rows[P6]["说明"].startswith("未测：")


@pytest.mark.parametrize("key", [P3, P4, P6])
def test_stream_function_checks_not_measured_without_synthetic_case(key):
    """无合成算例（= 只给结果对象）时必须未测，理由是结果对象不导出速度/通量场。"""
    rows = _rows(sanity_checks(_Res([np.zeros((4, 3))]), _geom()))
    assert rows[key]["通过"] is None
    assert rows[key]["说明"].startswith("未测：")
    assert "速度" in rows[key]["说明"] or "流函数" in rows[key]["说明"]


def test_sanity_checks_returns_exactly_seven_rows():
    rows = sanity_checks(_Res([np.zeros((4, 3))]), _geom())
    assert [r["检查项"] for r in rows] == [P1, P2, P3, P4, P5, P6, P7]
    assert all(set(r) == {"检查项", "通过", "实测值", "判据", "说明"} for r in rows)


def test_thresholds_are_the_brief_values():
    """判据常量必须逐字等于 brief 表（防有人为了让台账变绿而放宽）。"""
    assert FLUX_REL_TOL == 1.0e-9
    assert MONOTONE_TOL == 1.0e-9
    assert VELOCITY_AMPLIFICATION_MAX == 1.5
    assert RATE_RESPONSE_TOL_PP == 1.0
