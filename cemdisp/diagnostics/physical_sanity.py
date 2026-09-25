# -*- coding: utf-8 -*-
"""物理合理性闸门：7 项可证伪判据（Task 9）。

目的
----
用户要求"模型要符合物理实际，不只是代码跑通"。本模块把这句话落成 7 条**可证伪**的
物理判据，逐项返回 ``{检查项, 通过, 实测值, 判据, 说明}``，由
``scripts/entrypoints/verify_physical_sanity.py`` 汇总成台账。
判据**一律不做调参**：任何一项不通过都如实记录（台账是诊断输出，不是可调参数）。

判据表（逐字取自计划 Task 9；依据列为其物理来源）
------------------------------------------------
======  ==========================  ==================================================  ==========================
  #      检查项                       判据                                                 依据
======  ==========================  ==================================================  ==========================
 P-1     浓度有界                     ``0 ≤ c ≤ 1``（含全部快照）                          体积分数定义
 P-2     浓度单调（环空段）           每列 b 加权浓度随时间单调不减（容 1e-9 数值回退）      纯平流 + 无源项
 P-3     速度上界                     ``max|w| ≤ 1.5·Q/A_min``                             半环空单位通量 BC + 面积守恒
 P-4     总通量守恒                   每列 ``∫2H·w̄ dφ = const``（相对差 ≤1e-9）            (2.2) 单位通量归一
 P-5     窄边劣势                    密度稳定井（b>0）窄 1/4 行效率 ≤ 宽 1/4 行效率        偏心环空顶替基本物理
 P-6     冻结区静止                   ``wall>0.5`` 处 ``|w| ≤ 0.01·max|w|``                Pelipenko04 (2.6)-(2.8)
 P-7     排量响应方向                 排量 ×1.4 时域内 η_E 不下降超过 1 pp（方向检查）      现场经验方向
======  ==========================  ==================================================  ==========================

测点层级（**两类，说明列逐行标注**，coordinator 裁定 R155/#8）
-------------------------------------------------------------
* **P-1 / P-2 / P-5 / P-7 在结果对象上测**（``AnnulusSimulationResult`` + ``geom``）。
* **P-3 / P-4 / P-6 只能在求解器级合成算例上测**：``AnnulusSimulationResult``
  （``annulus_d2dga.py:305-339``）只导出 ``geom``/五个浓度与壁面场/``metrics``/
  ``depth_profiles``/``summary``/快照，``metric_columns``（``:2464-2482``）**无任何
  速度或通量列、无流函数场**；重建速度需要闭包状态，而闭包状态不落盘。⇒ 这三项经
  :func:`synthetic_stream_case` 直接调 ``solve_stream_function`` /
  ``velocity_from_stream_function`` 测量，说明列注明"单元级合成算例（非八井结果对象）"。
  **不写"未测"**（那是能测却不测），但**也不给生产 metrics 增列**（增列会改变每一次
  生产运行的产物与耗时，且 P-4 本质是离散算子的不变式，最该钉在算子旁边）。

合成算例的偏心放大解析式（P-3 阈值 1.5 的适用边界，**必须公开**）
------------------------------------------------------------------
合成算例用均匀间隙 + 常偏心度 e 的半环空、两层牛顿闭包（Z&F22 (4.21a)，
``I₁ ∝ H³``）。此时 φ-通量守恒给出 ``∂φΨ ∝ H³``，配合 (2.2) 得
``w ∝ H²``，于是

    max|w| / (Q/A_min) = (1+e)² / (1 + 1.5·e²)，  A_min = min_s 2∮b dy（全环空截面积）

（离散网格上实测与该式差 ≤0.2%，见 ``test_P3_matches_analytic_eccentric_amplification``。）
⇒ 判据 ``≤1.5`` 的适用边界是 **e ≲ 0.31**：
设计居中度 0.78–0.83（e=0.17~0.22）⇒ 比值 1.31~1.39 通过；
保守几何代理 0.65~0.70（e=0.30~0.35）⇒ 1.49~1.54，压在阈值上；
呼101 LEGACY 名义剖面 0.38~0.48（e≈0.52~0.62，代码自认"反推情景、不得作验证数字"）
⇒ 1.63~1.71 **超阈**。该边界写进 P-3 说明列，**不**通过调阈值掩盖。

P-6 的实现现状声明（诚实披露，不改判据）
----------------------------------------
``enable_stream_yield_gate`` 的 ``wall`` 自 2026-09-15 Task 11 起是**连续冻结度**
``clip(1 − τw_extrap/(f·τy), 0, 1)``（``annulus_d2dga.py:1096``），而 Pelipenko04
的判据是二值的（τw < f·τy ⇒ 静止）。合成算例取**二值掩码**（wall∈{0,1}），这是
production 的**最有利**配置（连续冻结度只会让 wall>0.5 的过渡带更宽、更不静止）。
即使在最有利配置下，判据在**冻结前锋那一格**也不成立（实测 |w|/max|w| ≈ 9.4e-2，
而冻结区内部 ≈3e-7）——原因是一阶/中心差分的 Ψ 梯度跨在 wall 的 0→1 台阶上。
该结论如实入账（``通过=False``），并在说明列给出前锋格/内部格的分解。

写盘边界：本模块只做纯后处理（零副作用），落盘由调用方负责，且只写
``results/内部自洽加固_2026-09-25/``。本模块**不**加入
``cemdisp/diagnostics/__init__.py`` 导出（R146 有意的不对称）。
"""
from __future__ import annotations

from typing import Any, Mapping, Optional

import numpy as np

# 判据常量（计划 Task 9 判据表逐字取值；禁止为了让台账变绿而放宽）
MONOTONE_TOL = 1.0e-9              # P-2：容许的数值回退
FLUX_REL_TOL = 1.0e-9             # P-4：逐列总通量的相对差上限
VELOCITY_AMPLIFICATION_MAX = 1.5  # P-3：偏心放大上限
FROZEN_SPEED_FRACTION_MAX = 0.01  # P-6：冻结区速度上限（相对 max|w|）
FROZEN_WALL_THRESHOLD = 0.5       # P-6：判"冻结"的 wall 门槛
RATE_RESPONSE_TOL_PP = 1.0        # P-7：η_E 容许下降（百分点）

# 检查项名（P-2 按 R138 更名：其判据落在**环空段**，"前缘"是管内段成对量）
CHECK_BOUNDED = "P-1 浓度有界"
CHECK_MONOTONE = "P-2 浓度单调（环空段，无源项平流）"
CHECK_VELOCITY_BOUND = "P-3 速度上界"
CHECK_FLUX_CONSERVATION = "P-4 总通量守恒"
CHECK_NARROW_DISADVANTAGE = "P-5 窄边劣势"
CHECK_FROZEN_REGION = "P-6 冻结区静止"
CHECK_RATE_RESPONSE = "P-7 排量响应方向"

CHECK_NAMES = (CHECK_BOUNDED, CHECK_MONOTONE, CHECK_VELOCITY_BOUND,
               CHECK_FLUX_CONSERVATION, CHECK_NARROW_DISADVANTAGE,
               CHECK_FROZEN_REGION, CHECK_RATE_RESPONSE)

# 测点层级（进说明列，避免三类判据混淆）
LEVEL_RESULT = "八井结果对象"
LEVEL_SYNTHETIC = "单元级合成算例（非八井结果对象）"

_NA_PREFIX = "未测："

# 合成算例：均匀间隙 + 常偏心度半环空（几何量级取现场尾管井，见 synthetic_stream_case）
_SYNTH_HOLE_MM = 215.9
_SYNTH_OD_MM = 168.3
_SYNTH_STANDOFF = 0.78      # 现场设计中居中度下限（呼2 设计 78%、呼1-003 78%）
_SYNTH_ETA1 = 0.05          # 被顶替液（钻井液）表观黏度 Pa·s
_SYNTH_ETA2 = 0.10          # 顶替液（水泥浆）表观黏度 Pa·s
_SYNTH_M = _SYNTH_ETA1 / _SYNTH_ETA2
_SYNTH_C_BAR = 0.5          # 间隙平均水泥体积分数（两层闭包非退化）
_SYNTH_Q_M3S = 0.55 / 60.0  # 现场排量量级 0.55 m³/min
_SYNTH_NY = 41
_SYNTH_NZ = 21
_SYNTH_LENGTH_M = 900.0


# --------------------------------------------------------------------------- #
# 行构造
# --------------------------------------------------------------------------- #
def _row(name: str, passed: Optional[bool], measured: str, criterion: str,
         note: str) -> dict:
    """按统一 schema 组装一行；``passed=None`` 必须走 :func:`_not_measured`。"""
    if passed is None:
        raise ValueError("通过=None 必须走 _not_measured（说明须以'未测：'开头）")
    return {"检查项": name, "通过": bool(passed), "实测值": measured,
            "判据": criterion, "说明": note}


def _not_measured(name: str, criterion: str, reason: str, measured: str = "未测") -> dict:
    """未测行：``通过=None`` + 说明以 ``未测：`` 开头（**不得**默认为通过）。"""
    return {"检查项": name, "通过": None, "实测值": measured,
            "判据": criterion, "说明": f"{_NA_PREFIX}{reason}"}


def _note(level: str, *parts: str) -> str:
    return "; ".join((f"测点层级={level}",) + tuple(p for p in parts if p))


def _snapshots(result: Any) -> tuple:
    snaps = getattr(result, "cement_snapshots", None)
    return tuple(snaps) if snaps else ()


def _final_result_summary(result: Any) -> Mapping[str, Any]:
    summary = getattr(result, "summary", None) or {}
    block = summary.get("最终结果")
    return block if isinstance(block, Mapping) else {}


def _eta_e_of(result: Any) -> Optional[float]:
    value = _final_result_summary(result).get("全井段最终有效顶替效率")
    try:
        return None if value is None else float(value)
    except (TypeError, ValueError):
        return None


def _column_mean(field: np.ndarray, geom: Mapping[str, Any]) -> np.ndarray:
    """每列（固定 s、沿 φ 求和）的 b 加权平均浓度。"""
    b = np.asarray(geom["b"], dtype=float)
    c = np.asarray(field, dtype=float)
    return (b * c).sum(axis=0) / np.maximum(b.sum(axis=0), 1.0e-30)


def _wide_quarter_efficiency(cement: np.ndarray, geom: Mapping[str, Any]) -> float:
    """宽边 1/4 行效率：与 ``displacement_metrics._narrow_quarter_efficiency`` 同口径，
    仅把行切片从"最后 n_q 行（窄边）"换成"前 n_q 行（宽边，φ=0 侧，h = H(1+e·cos πφ) 最大）"。
    """
    from cemdisp.diagnostics.displacement_metrics import _EPS, _trapz2d  # 延迟导入，防循环导入

    n_q = max(1, int(cement.shape[0]) // 4)
    b_q = np.asarray(geom["b"], dtype=float)[:n_q, :]
    c_q = np.asarray(cement, dtype=float)[:n_q, :]
    y_q = np.asarray(geom["y"], dtype=float)[:n_q]
    num = _trapz2d(b_q * c_q, y_q, geom["s"])
    den = _trapz2d(b_q, y_q, geom["s"])
    return float(np.clip(num / max(den, _EPS), 0.0, 1.0))


# --------------------------------------------------------------------------- #
# 求解器级合成算例（P-3/P-4/P-6 的测点）
# --------------------------------------------------------------------------- #
def synthetic_stream_case(*, standoff: float = _SYNTH_STANDOFF, ny: int = _SYNTH_NY,
                          nz: int = _SYNTH_NZ, hole_mm: float = _SYNTH_HOLE_MM,
                          od_mm: float = _SYNTH_OD_MM, eta1: float = _SYNTH_ETA1,
                          eta2: float = _SYNTH_ETA2, m: float = _SYNTH_M,
                          c_bar: float = _SYNTH_C_BAR, Q_m3s: float = _SYNTH_Q_M3S,
                          freeze_rows: Optional[int] = None,
                          apply_wall: bool = True, banded: bool = True) -> dict:
    """构造单元级合成算例并解一次 (4.22)，供 P-3/P-4/P-6 直接测量。

    几何：均匀间隙半环空，``H(φ) = (hole−od)/4 · (1 + e·cos πφ)``（与生产
    ``annulus_d2dga._build_geom`` 的 ``e = 1 − standoff``、``b = 2H`` 同口径），
    ``e = 1 − standoff``。默认偏心度取**现场设计中居中度下限 0.78**（呼2 设计值 78%、
    呼1-003 78%），即 e=0.22 —— 这是 8 井设计中居中度区间 [0.78, 0.83] 里**最不利**
    的一端（不是为了让判据通过而挑的中间值）。

    Args:
        standoff: 居中度（1.0 = 完全居中）；``e = 1 − standoff``。
        freeze_rows: 冻结的行数（自 φ=1 窄边侧起数）；``None``（默认）⇒ 不构造掩码
            （wall=None，等价 ``enable_stream_yield_gate=False``——P-3/P-4 的基线算例）；
            ``"auto"`` ⇒ ``max(1, ny//4)``（窄 1/4，P-6 的算例）；``0`` ⇒ 掩码全零。
        apply_wall: 是否把 ``wall`` 送进 ``solve_stream_function`` 的算子。
            ``False`` = **反例算例**：掩码已算出但未被消费（R42 死开关类缺陷）。
        banded: 线性求解器选择（与生产默认一致，Task 8）。

    Returns:
        dict，含 ``geom`` / ``psi``（单位通量 Ψ）/ ``w_unit``（模块 w̄）/ ``w_m_s``
        （物理量纲，按 ``w = w_unit·(Q/2)/π``，与 ``_velocity_stream_function`` 同式）/
        ``wall`` / ``wall_in_operator`` / ``Q_m3s`` / ``A_min_m2`` / ``standoff`` /
        ``e`` / ``测点层级``。
    """
    from cemdisp.models2d.stream_function import (  # 延迟导入：与 annulus_d2dga 解耦
        solve_stream_function,
        velocity_from_stream_function,
    )

    e = float(1.0 - float(standoff))
    clearance_m = (float(hole_mm) - float(od_mm)) / 1000.0
    half_gap_mean = clearance_m / 2.0
    mean_radius = ((float(hole_mm) + float(od_mm)) / 4.0) / 1000.0
    y = np.linspace(0.0, np.pi * mean_radius, int(ny))
    phi = y / y[-1]
    H = half_gap_mean * (1.0 + e * np.cos(np.pi * phi))[:, None] * np.ones((1, int(nz)))
    geom = {"phi": phi, "s": np.linspace(0.0, _SYNTH_LENGTH_M, int(nz)), "H": H,
            "b": 2.0 * H, "y": y,
            "hole_mm": np.full(int(nz), float(hole_mm)),
            "od_mm": np.full(int(nz), float(od_mm))}

    if freeze_rows == "auto":
        freeze_rows = max(1, int(ny) // 4)
    if freeze_rows is None:
        wall = None
    else:
        wall = np.zeros((int(ny), int(nz)))
        if int(freeze_rows) > 0:
            wall[-int(freeze_rows):, :] = 1.0

    psi = solve_stream_function(
        geom, np.full((int(ny), int(nz)), float(c_bar)), float(eta1), float(eta2),
        float(m), np.zeros((2, int(ny), int(nz))),
        wall=(wall if (wall is not None and apply_wall) else None),
        banded=banded)
    w_unit, _v_unit = velocity_from_stream_function(psi, geom)
    q_half = float(Q_m3s) / 2.0
    # 与 annulus_d2dga._velocity_stream_function:1732-1734 逐字同式
    w_m_s = w_unit * (q_half / np.pi)
    # 全环空截面积 A(s) = 2·∮b dy（半环空积分 ×2，与 _phase_volume 的因子 2 同口径）
    area_per_s = 2.0 * np.trapezoid(geom["b"], x=y, axis=0)
    return {"geom": geom, "psi": psi, "w_unit": w_unit, "w_m_s": w_m_s, "wall": wall,
            "wall_in_operator": bool(apply_wall and wall is not None),
            "Q_m3s": float(Q_m3s), "A_min_m2": float(np.min(area_per_s)),
            "standoff": float(standoff), "e": e, "测点层级": LEVEL_SYNTHETIC,
            "case": (f"均匀间隙半环空 hole={hole_mm}mm/od={od_mm}mm "
                     f"居中度={standoff:.2f}(e={e:.2f}) {ny}×{nz} "
                     f"c̄={c_bar} η₁={eta1}/η₂={eta2} Q={Q_m3s * 60000:.2f}L/min")}


# --------------------------------------------------------------------------- #
# 逐项检查
# --------------------------------------------------------------------------- #
def check_bounded(result: Any) -> dict:
    """P-1：全部浓度快照都必须落在 [0, 1]（体积分数定义）。"""
    criterion = "0 ≤ c ≤ 1（含全部快照）"
    snaps = _snapshots(result)
    if not snaps:
        return _not_measured(CHECK_BOUNDED, criterion,
                             "结果对象无浓度快照（cement_snapshots 为空）")
    cmin = min(float(np.min(s)) for s in snaps)
    cmax = max(float(np.max(s)) for s in snaps)
    passed = bool(cmin >= 0.0 and cmax <= 1.0)   # NaN 比较为 False ⇒ 一并抓住
    measured = f"min={cmin:.6g}, max={cmax:.6g}（{len(snaps)} 个快照）"
    return _row(CHECK_BOUNDED, passed, measured, criterion,
                _note(LEVEL_RESULT,
                      "NaN 会使比较为假并被判不通过" if not np.isfinite(cmin * cmax) else ""))


def check_monotone(result: Any, geom: Mapping[str, Any]) -> dict:
    """P-2：每列 b 加权浓度随时间单调不减（纯平流 + 无源项；容 1e-9 数值回退）。"""
    criterion = f"每列 b 加权浓度随时间单调不减（容许 {MONOTONE_TOL:.0e} 数值回退）"
    snaps = _snapshots(result)
    if len(snaps) < 2:
        return _not_measured(CHECK_MONOTONE, criterion,
                             f"快照不足（{len(snaps)} 个，单调性需 ≥2 个）")
    cols = [_column_mean(s, geom) for s in snaps]
    worst = 0.0
    for prev, cur in zip(cols[:-1], cols[1:]):
        worst = min(worst, float(np.min(cur - prev)))
    backstep = max(0.0, -worst)
    passed = backstep <= MONOTONE_TOL
    measured = f"最大回退={backstep:.3e}（{len(snaps)} 快照 × {cols[0].size} 列）"
    return _row(CHECK_MONOTONE, passed, measured, criterion,
                _note(LEVEL_RESULT, "环空段量：列均值沿 φ 用 b 加权；窗口止于尾浆入环空时刻，"
                                    "其后继流体尚未进入，故不应出现稀释型回退"))


def check_velocity_bound(case: Optional[Mapping[str, Any]]) -> dict:
    """P-3：``max|w| ≤ 1.5·Q/A_min``（半环空单位通量 BC + 面积守恒）。

    A_min = 全域最小全环空截面积 ``min_s 2·∮b dy``。合成算例的解析放大比见模块
    docstring（``(1+e)²/(1+1.5e²)``，1.5 倍的适用边界 e ≈ 0.31）。
    """
    criterion = f"max|w| ≤ {VELOCITY_AMPLIFICATION_MAX}·Q/A_min"
    if case is None:
        return _not_measured(CHECK_VELOCITY_BOUND, criterion,
                             "未提供求解器级合成算例；结果对象不导出速度/流函数场"
                             "（metric_columns 无速度列，闭包状态不落盘）")
    w = np.asarray(case["w_m_s"], dtype=float)
    q_over_area = float(case["Q_m3s"]) / float(case["A_min_m2"])
    measure = float(np.max(np.abs(w)))
    ratio = measure / q_over_area
    passed = ratio <= VELOCITY_AMPLIFICATION_MAX
    e = float(case.get("e", float("nan")))
    measured = (f"max|w|={measure:.6g} m/s, Q/A_min={q_over_area:.6g} m/s, "
                f"放大比={ratio:.4f}")
    return _row(CHECK_VELOCITY_BOUND, passed, measured, criterion,
                _note(LEVEL_SYNTHETIC, str(case.get("case", "")),
                      f"解析放大比(1+e)²/(1+1.5e²)={((1 + e) ** 2 / (1 + 1.5 * e * e)):.4f}"
                      f"（e={e:.2f}）",
                      f"阈值 {VELOCITY_AMPLIFICATION_MAX} 的适用边界 e≈0.31：设计居中度 "
                      f"0.78~0.83(e=0.17~0.22) 通过、保守代理 0.65~0.70(e=0.30~0.35) 压线、"
                      f"呼101 LEGACY 名义剖面 0.38~0.48(e≈0.52~0.62) 超阈"))


def check_flux_conservation(case: Optional[Mapping[str, Any]]) -> dict:
    """P-4：每列 ``∫2H·w̄ dφ`` 必须为同一常数（(2.2) 单位通量归一 + 差分-梯形恒等式）。

    这里用**单位通量**口径的 w̄ 与 Ψ 同源的几何积分；浮力只在列内重分配轴向通量，
    不改变列总量（``stream_function`` 模块 docstring 的"Q 缩放锚（不变量）"）。
    """
    criterion = f"每列 ∫2H·w̄ dφ 相对差 ≤ {FLUX_REL_TOL:.0e}"
    if case is None:
        return _not_measured(CHECK_FLUX_CONSERVATION, criterion,
                             "未提供求解器级合成算例；结果对象不导出速度/通量/流函数场")
    geom = case["geom"]
    col_flux = np.trapezoid(2.0 * np.asarray(geom["H"], dtype=float)
                            * np.asarray(case["w_unit"], dtype=float),
                            x=np.asarray(geom["phi"], dtype=float), axis=0)
    scale = max(abs(float(np.mean(col_flux))), 1.0e-30)
    rel = float((np.max(col_flux) - np.min(col_flux)) / scale)
    passed = rel <= FLUX_REL_TOL
    measured = (f"列通量 min={np.min(col_flux):.12e}, max={np.max(col_flux):.12e}, "
                f"相对差={rel:.3e}")
    return _row(CHECK_FLUX_CONSERVATION, passed, measured, criterion,
                _note(LEVEL_SYNTHETIC, str(case.get("case", "")),
                      "浮力场为零 ⇒ 本算例只检验算子自身的列守恒；浮力不改变列总量，"
                      "故该项对浮力场无关"))


def check_narrow_disadvantage(result: Any, geom: Mapping[str, Any]) -> dict:
    """P-5：密度稳定井（浮力数 b>0）窄 1/4 行效率 ≤ 宽 1/4 行效率。

    ``b`` 取 ``summary["最终结果"]["浮力数_b"]``（Z&F22 p.8 口径，重驱轻 ⇒ b>0 稳定）。
    b ≤ 0 或缺失 ⇒ 判据不适用 ⇒ **未测**（不得默认通过）。
    """
    criterion = "密度稳定井（b>0）窄 1/4 行效率 ≤ 宽 1/4 行效率"
    cement = getattr(result, "cement_field", None)
    if cement is None or np.asarray(cement).size == 0:
        return _not_measured(CHECK_NARROW_DISADVANTAGE, criterion,
                             "结果对象无最终浓度场（cement_field）")
    b_num = _final_result_summary(result).get("浮力数_b")
    if b_num is None:
        return _not_measured(CHECK_NARROW_DISADVANTAGE, criterion,
                             "summary 无浮力数 b（无法判定密度稳定性）")
    if not float(b_num) > 0.0:
        return _not_measured(CHECK_NARROW_DISADVANTAGE, criterion,
                             f"密度不稳定/倒置井（浮力数 b={float(b_num):+.3f} ≤ 0），"
                             "本判据只适用于密度稳定井（b>0）",
                             measured=f"b={float(b_num):+.3f}")
    from cemdisp.diagnostics.displacement_metrics import _narrow_quarter_efficiency

    cement = np.asarray(cement, dtype=float)
    eta_narrow = float(_narrow_quarter_efficiency(cement, geom))
    eta_wide = _wide_quarter_efficiency(cement, geom)
    passed = eta_narrow <= eta_wide
    measured = (f"η_窄1/4={eta_narrow:.6f}, η_宽1/4={eta_wide:.6f}, "
                f"差={eta_narrow - eta_wide:+.6f}, b={float(b_num):+.3f}")
    return _row(CHECK_NARROW_DISADVANTAGE, passed, measured, criterion,
                _note(LEVEL_RESULT, "行序约定：y[0]=宽边（φ=0，间隙最大）、"
                                    "y[-1]=窄边（φ=1，间隙最小），与求解器 front_narrow 一致"))


def check_frozen_region(case: Optional[Mapping[str, Any]]) -> dict:
    """P-6：``wall>0.5`` 处 ``|w| ≤ 0.01·max|w|``（Pelipenko04 (2.6)-(2.8) 停流区）。

    ``wall=None``（屈服门未启用）或无 ``wall>0.5`` 单元 ⇒ **未测**（不得默认通过）。
    """
    criterion = f"wall>{FROZEN_WALL_THRESHOLD} 处 |w| ≤ {FROZEN_SPEED_FRACTION_MAX}·max|w|"
    if case is None:
        return _not_measured(CHECK_FROZEN_REGION, criterion,
                             "未提供求解器级合成算例；结果对象不导出速度场")
    wall = case.get("wall")
    if wall is None:
        return _not_measured(CHECK_FROZEN_REGION, criterion,
                             "合成算例未启用屈服门（wall=None，等价 "
                             "enable_stream_yield_gate=False）⇒ 无冻结单元")
    wall = np.asarray(wall, dtype=float)
    mask = wall > FROZEN_WALL_THRESHOLD
    if not bool(np.any(mask)):
        return _not_measured(CHECK_FROZEN_REGION, criterion,
                             f"合成算例无 wall>{FROZEN_WALL_THRESHOLD} 的冻结单元")
    w = np.abs(np.asarray(case["w_m_s"], dtype=float))
    peak = float(np.max(w))
    if peak <= 0.0:
        return _not_measured(CHECK_FROZEN_REGION, criterion,
                             "速度场恒为零（无流动），判据退化无判别力")
    frozen_peak = float(np.max(w[mask]))
    ratio = frozen_peak / peak
    passed = ratio <= FROZEN_SPEED_FRACTION_MAX
    front = _front_row_ratio(w, mask)
    measured = (f"冻结区 max|w|={frozen_peak:.3e} m/s, max|w|={peak:.6g} m/s, "
                f"比值={ratio:.3e}")
    return _row(CHECK_FROZEN_REGION, passed, measured, criterion,
                _note(LEVEL_SYNTHETIC, str(case.get("case", "")),
                      "二值掩码（wall∈{0,1}）是 production 连续冻结度（Task 11 起）的"
                      "最有利配置",
                      f"分解：冻结前锋格比值={front[0]:.3e}（{front[1]} 格），"
                      f"其余冻结格最大比值={front[2]:.3e}" if front else ""))


def _front_row_ratio(w: np.ndarray, mask: np.ndarray) -> tuple:
    """拆开"冻结前锋格"与"冻结区内部格"，返回 (前锋最大比值, 前锋格数, 内部最大比值)。

    前锋格 = 与未冻结行相邻的冻结行（差分的 Ψ 梯度跨在 wall 的 0→1 台阶上）；
    内部格 = 冻结行中不与未冻结行相邻者（真正的"静止壁层"应在这里）。
    """
    rows = np.where(mask.any(axis=1))[0]
    if rows.size == 0:
        return ()
    masked_rows = set(int(r) for r in rows)
    n_rows = int(mask.shape[0])
    front_rows = [r for r in masked_rows
                  if (r - 1 >= 0 and (r - 1) not in masked_rows)
                  or (r + 1 < n_rows and (r + 1) not in masked_rows)]
    inner_rows = sorted(masked_rows - set(front_rows))
    peak = float(np.max(w))
    front = mask.copy()
    front[:] = False
    for r in front_rows:
        front[r] = mask[r]
    inner = mask & ~front
    front_vals, inner_vals = w[front], w[inner]
    return (float(np.max(front_vals)) / peak if front_vals.size else 0.0,
            int(front_vals.size),
            float(np.max(inner_vals)) / peak if inner_vals.size else 0.0)


def check_rate_response(result: Any, rate_result: Any) -> dict:
    """P-7：排量 ×1.4 时域内 η_E 不下降超过 1 pp（**方向检查，不设幅值门**）。

    需两次运行（基线 + 排量缩放）；缺对照运行 ⇒ **未测**。
    """
    criterion = f"η_E(排量×1.4) ≥ η_E(基线) − {RATE_RESPONSE_TOL_PP} pp"
    if rate_result is None:
        return _not_measured(CHECK_RATE_RESPONSE, criterion,
                             "未提供排量×1.4 对照运行结果（本判据需两次运行）")
    eta_base = _eta_e_of(result)
    eta_rate = _eta_e_of(rate_result)
    if eta_base is None or eta_rate is None:
        return _not_measured(CHECK_RATE_RESPONSE, criterion,
                             "summary 缺'全井段最终有效顶替效率'，无法比较")
    delta_pp = (eta_rate - eta_base) * 100.0
    passed = delta_pp >= -RATE_RESPONSE_TOL_PP
    measured = (f"η_E(基线)={eta_base:.6f} → η_E(排量×1.4)={eta_rate:.6f}, "
                f"Δ={delta_pp:+.3f} pp")
    return _row(CHECK_RATE_RESPONSE, passed, measured, criterion,
                _note(LEVEL_RESULT, "方向检查：只判趋势方向与 1 pp 容许，不设幅值门"))


# --------------------------------------------------------------------------- #
# 汇总入口
# --------------------------------------------------------------------------- #
def sanity_checks(result: Any, geom: Mapping[str, Any], *,
                  rate_result: Any = None,
                  stream_case: Optional[Mapping[str, Any]] = None,
                  frozen_case: Optional[Mapping[str, Any]] = None) -> list:
    """对一条结果（及可选对照/合成算例）跑 7 项物理判据，返回逐项 dict 列表。

    Args:
        result: ``AnnulusSimulationResult``（或其等价替身）。
        geom: 求解器几何字典（需 ``b``；P-2/P-5 另需 ``y``/``s``/``phi``）。
        rate_result: 排量 ×1.4 的对照结果（P-7）；``None`` ⇒ P-7 未测。
        stream_case: :func:`synthetic_stream_case` 的产物（P-3/P-4 用，**不带冻结带**，
            以便把偏心放大与屈服门分流两个效应分开）；``None`` ⇒ 未测。
        frozen_case: 带 ``wall`` 掩码的合成算例（P-6 用）；``None`` ⇒ P-6 未测。

    Returns:
        长度 7 的列表，元素为 ``{检查项, 通过, 实测值, 判据, 说明}``；
        ``通过`` 取 ``True``/``False``/``None``，``None`` 表示未测且 ``说明`` 以
        ``未测：`` 开头。
    """
    return [check_bounded(result),
            check_monotone(result, geom),
            check_velocity_bound(stream_case),
            check_flux_conservation(stream_case),
            check_narrow_disadvantage(result, geom),
            check_frozen_region(frozen_case),
            check_rate_response(result, rate_result)]
