# -*- coding: utf-8 -*-
"""物理合理性闸门：7 项判据（其中 P-4 为算子恒等式，非独立测量）（Task 9）。

⚠️ **计数口径（终审 I-2，2026-09-26 用户裁定 R195）**：本模块共 **7 项判据**，
但其中 **P-4「总通量守恒」是离散算子的代数恒等式、不是独立测量** ——
``w̄ = ∂φΨ/(2r_aH)`` ⇒ ``2H·w̄ = ∂φΨ/r_a``，而「中心差分 + 梯形」逐列严格等于
``(Ψ(1)−Ψ(0))/r_a``、每列 Dirichlet 值相同 ⇒ **对任意内部解恒成立，无法在求解器
输出上失败**（终审实证：注入强浮力场后列通量相对离散度仍 3.4e-16）。故**可独立
证伪的判据是 6 项**，写作"7 项可证伪判据"属**多算一项**。判据实现与阈值**未**改动
（改 P-4 的实现＝改判据口径，属新工作），只在计数表述上如实标注。

目的
----
用户要求"模型要符合物理实际，不只是代码跑通"。本模块把这句话落成 7 条物理判据
（其中 6 条**可独立证伪**，P-4 为算子恒等式，见上），逐项返回
``{检查项, 通过, 实测值, 判据, 说明}``，由
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
         （算子恒等式，非独立测量——见下方 ⚠️ 计数口径；阈值与实现未改）
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
八井 standoff 剖面**实测**区间 **0.76–0.83**（e=0.17~0.24；呼2 鞋底 0.76、
呼103 设计代理 0.778、呼1-003/004 0.83）⇒ 比值 1.31~1.415 通过
（**但余量薄**：本行算例取 0.78 ⇒ 比值 1.3865 / 阈值 1.5 = **1.08× 余量**，
是全档最薄的一处；0.76 档 1.415 ⇒ 1.06×；降到 0.65 即 1.49 ⇒ 压在阈值上）；
保守几何代理 0.65~0.70（e=0.30~0.35）⇒ 1.49~1.54，压在阈值上；
呼101 LEGACY 名义剖面 0.38~0.48（e≈0.52~0.62，代码自认"反推情景、不得作验证数字"）
⇒ 1.63~1.71 **超阈**。该边界写进 P-3 说明列，**不**通过调阈值掩盖。

P-6 判定落在「完全冻结格」上（判据前提的字面成立处；协调者裁定 2026-09-26）
--------------------------------------------------------------------------
计划 P-6 写作 ``wall>0.5 处 |w| ≤ 0.01·max|w|``，其前提是 **wall 二值**。但
``enable_stream_yield_gate`` 的 ``wall`` 自提交 **``9448572``**「fix(yield): 屈服门
二值→连续，消除阈值悬崖（Pelipenko04 (2.6)-(2.8)）」起是**连续冻结度**
``clip(1 − τw_extrap/(f·τy), 0, 1)``（``annulus_d2dga.py:1087-1096``；该提交的
**目的就是消除**旧 ``np.where(immobile, 1, 0)`` 的 0/1 悬崖）。⇒ ``wall>0.5`` 是
**软**阈值，把过渡带格判成「必须静止」等于要求模型在被刻意连续化的地方按二值行为。
故本模块按裁定改判：

* **判定部分**：``wall ≥ 1 − 1e-6``（τw = 0 的完全冻结格）处 ``|w| ≤ 0.01·max|w|``
  （计划数值阈值逐字沿用，只用在它自己的前提字面成立处）；
* **非空过守卫**：报告完全冻结格数；为 0 ⇒ ``通过=None`` / 未测（给出原因），绝不空过；
* **过渡带（``0.5 < wall < 1−1e-6``）转为披露量**：格数、最高 ``|w|/max|w|``、该格
  自身的 ``wall``，与「相对同算例未加冻结门时同格速度的抑制倍数」一并折进**同一行**。

合成算例的冻结度剖面取 ``wall(φ) = sin²(π·t/2)``（``t = clip((φ−0.625)/0.25, 0, 1)``，
自宽边 wall=0 到窄边 wall=1，两端导数为零的 C¹ 过渡）：production 的
``τw_extrap ∝ b = 2H`` 在 φ 上只有 ~1.56× 展布（且 ``wall=1`` 只在 ``H→0`` 极限或
``col_freeze`` 整列无流动分支出现），跨越阈值时冻结度在整个 φ 上连续变化 ⇒ 单调 C¹ 剖面
是这一连续过渡的最小光滑代表，且使 ``wall≡1`` 的格确实静止。分段线性的过渡剖面会在
**第一个完全冻结格**留下 O(1/Δφ) 的中心差分跨台阶泄漏（实测 1.10e-2 > 0.01）——那反映的是
**离散格式在冻结度突变处的局限**，不是冻结区内部在流动。

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
FROZEN_WALL_THRESHOLD = 0.5       # P-6：判"过渡带"的软阈值（wall>0.5）
# P-6：判"完全冻结"的阈值。`wall` 是**连续**冻结度 `clip(1−τw_extrap/(f·τy),0,1)`
# （自 9448572「屈服门二值→连续，消除阈值悬崖」起），故 `wall=1` ⟺ τw=0 ⟺ 完全冻结；
# `wall>0.5` 是**软**阈值，落在过渡带里的格不能按"必须静止"判。
FROZEN_WALL_FULLY = 1.0 - 1.0e-6
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
_SYNTH_STANDOFF = 0.78      # 设计档下限（呼2 78%、呼1-003 78%）；实测区间下限是 0.76
                            # （呼2 鞋底）⇒ 本档放大比 1.3865/阈值 1.5 = 1.08× 余量（最薄）
_SYNTH_ETA1 = 0.05          # 被顶替液（钻井液）表观黏度 Pa·s
_SYNTH_ETA2 = 0.10          # 顶替液（水泥浆）表观黏度 Pa·s
_SYNTH_M = _SYNTH_ETA1 / _SYNTH_ETA2
_SYNTH_C_BAR = 0.5          # 间隙平均水泥体积分数（两层闭包非退化）
_SYNTH_Q_M3S = 0.55 / 60.0  # 现场排量量级 0.55 m³/min
_SYNTH_NY = 41
_SYNTH_NZ = 21
_SYNTH_LENGTH_M = 900.0
# 合成冻结度剖面的斜坡区间（归一化 φ）：φ≤0.625 wall=0、φ≥0.875 wall=1，
# 中间为 C¹（sin²）过渡。仅作"连续冻结度"的形状代表，见 synthetic_stream_case。
_SYNTH_FREEZE_PHI0 = 0.625
_SYNTH_FREEZE_SPAN = 0.25


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
                          freeze: Optional[str] = None,
                          apply_wall: bool = True, banded: bool = True) -> dict:
    """构造单元级合成算例并解一次 (4.22)，供 P-3/P-4/P-6 直接测量。

    几何：均匀间隙半环空，``H(φ) = (hole−od)/4 · (1 + e·cos πφ)``（与生产
    ``annulus_d2dga._build_geom`` 的 ``e = 1 − standoff``、``b = 2H`` 同口径），
    ``e = 1 − standoff``。默认偏心度取 **0.78**（呼2 设计值 78%、呼1-003 78%），
    即 e=0.22 —— 这是 8 井 standoff 剖面**实测**区间 **[0.76, 0.83]** 里**最不利**的
    一端（不是为了让判据通过而挑的中间值）。⚠️ 该档的放大比 1.3865 相对阈值 1.5
    **只有 1.08× 余量**（全档最薄；0.76 档 1.415 ⇒ 1.06×），见模块 docstring 的
    适用边界段与 P-3 说明列。

    Args:
        standoff: 居中度（1.0 = 完全居中）；``e = 1 − standoff``。
        freeze: 冻结度剖面（P-6 用；P-3/P-4 传 ``None``）：

            * ``None``（默认）⇒ ``wall=None``，等价 ``enable_stream_yield_gate=False``；
            * ``"smooth"`` ⇒ ``wall(φ) = sin²(π·t/2)``、``t = clip((φ−0.625)/0.25, 0, 1)``：
              自宽边（φ=0.625、wall=0）到窄边（φ=0.875、wall=1）的**连续**冻结度，
              两端导数为零（C¹）。既有 ``wall≡1`` 的完全冻结格，也有 ``0.5<wall<1``
              的过渡带格；
            * ``"partial"`` ⇒ ``0.9·sin²(·)``：有过渡带格但**无**完全冻结格
              （用于"非空过"守卫：判据须判为未测）；
            * ``"zero"`` ⇒ 掩码全零（无任何冻结格）。
        apply_wall: 是否把 ``wall`` 送进 ``solve_stream_function`` 的算子。
            ``False`` = **反例算例**：掩码已声明但未进算子（R42 死开关类缺陷）。
        banded: 线性求解器选择（与生产默认一致，Task 8）。

    Returns:
        dict，含 ``geom`` / ``psi``（单位通量 Ψ）/ ``w_unit``（模块 w̄）/ ``w_m_s``
        （物理量纲，按 ``w = w_unit·(Q/2)/π``，与 ``_velocity_stream_function`` 同式）/
        ``wall`` / ``wall_in_operator`` / ``Q_m3s`` / ``A_min_m2`` / ``standoff`` /
        ``e`` / ``测点层级``；另含重解所需的 ``c_bar``/``eta1``/``eta2``/``m``/``banded``
        （供 :func:`case_velocity` 构造反例）。
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

    wall = _freeze_profile(freeze, phi, int(ny), int(nz))
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
            "c_bar": float(c_bar), "eta1": float(eta1), "eta2": float(eta2),
            "m": float(m), "banded": bool(banded),
            "case": (f"均匀间隙半环空 hole={hole_mm}mm/od={od_mm}mm "
                     f"居中度={standoff:.2f}(e={e:.2f}) {ny}×{nz} "
                     f"c̄={c_bar} η₁={eta1}/η₂={eta2} Q={Q_m3s * 60000:.2f}L/min"
                     f" 冻结剖面={freeze}")}


def _freeze_profile(freeze: Optional[str], phi: np.ndarray, ny: int, nz: int):
    """按 ``freeze`` 名构造冻结度剖面 (ny,nz)；``None`` ⇒ 无掩码（返回 None）。"""
    if freeze is None:
        return None
    if freeze == "zero":
        return np.zeros((ny, nz))
    t = np.clip((np.asarray(phi, dtype=float) - _SYNTH_FREEZE_PHI0) / _SYNTH_FREEZE_SPAN,
                0.0, 1.0)
    scale = 0.9 if freeze == "partial" else 1.0
    if freeze not in ("smooth", "partial"):
        raise ValueError(f"freeze 只能是 None/'smooth'/'partial'/'zero'，得到 {freeze!r}")
    return np.repeat((scale * np.sin(np.pi * t / 2.0) ** 2)[:, None], nz, axis=1)


def case_velocity(case: Mapping[str, Any], wall) -> np.ndarray:
    """用**给定的** wall 掩码把同一合成算例重解一次，返回物理轴向速度场 (m/s)。

    用途（反例构造）：把"判据所用的 wall 掩码"与"真正送进算子的 wall"分开，即可复现
    "掩码声明冻结、算子却未冻结"（R42 死开关类）这一类缺陷——被声明为完全冻结的格在
    该情形下仍在流动。缩放式与 :func:`synthetic_stream_case` 逐字同源。
    """
    from cemdisp.models2d.stream_function import (
        solve_stream_function,
        velocity_from_stream_function,
    )

    geom = case["geom"]
    ny, nz = np.asarray(geom["H"]).shape
    psi = solve_stream_function(
        geom, np.full((ny, nz), float(case["c_bar"])), float(case["eta1"]),
        float(case["eta2"]), float(case["m"]), np.zeros((2, ny, nz)),
        wall=wall, banded=bool(case.get("banded", True)))
    w_unit, _v = velocity_from_stream_function(psi, geom)
    return w_unit * (float(case["Q_m3s"]) / 2.0 / np.pi)


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
                      "**本行是上游 clip 的回归/NaN 守卫，不是有界性的独立物理证据**："
                      "快照由求解器以 ``cement = np.clip(lead + tail, 0.0, 1.0)`` 落盘"
                      "（``annulus_d2dga.py:2396``/``:2404`` 处构造、``:2397`` 处 append；"
                      "``lead``/``tail`` 更早亦已 clip），故模型**不可能**越界——"
                      "实测 min=0/max=1 正是该 clip 的签名。只有 NaN（能穿过 ``np.clip``）"
                      "或未来重构移除该 clip 才会让本行转红",
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
                      f"判据适用域（**域内近似**，非全域断言）：放大比 (1+e)²/(1+1.5e²)，阈值 "
                      f"{VELOCITY_AMPLIFICATION_MAX} 只在 e≲0.31（居中度≳0.69）成立 ⇒ "
                      f"八井 standoff 剖面实测区间 [0.76, 0.83] 在域内"
                      f"（0.76=呼2 鞋底、0.778=呼103 设计代理、0.83=呼1-003/004）；"
                      f"本行算例取 0.78（e=0.22），是逐井最小值 0.76 之外**最不利**的一档"
                      f"（0.76 ⇒ e=0.24 ⇒ 放大比 1.415，仍在域内）；"
                      f"⚠️ **余量薄（M-4，2026-09-26 终审）**：本行取 0.78 ⇒ 放大比 "
                      f"1.3865 / 阈值 1.5 = **1.08× 余量**（全档最薄；0.76 档 1.415 ⇒ 1.06×）"
                      f"⇒ 判据对居中度**高度敏感**：降到 0.72 即 1.447、0.70 即 1.466、"
                      f"0.65 即 1.491（压在阈值上）。**不**因此调阈值或改所选取的 standoff——"
                      f"余量如实记录；"
                      f"把居中度降到呼102 保守几何代理 0.65(e=0.35) 会转红——"
                      f"那是**判据适用域的边界，不是模型回归**；呼101 LEGACY 名义剖面 "
                      f"0.38~0.48(e≈0.52~0.62，代码自认「反推情景、不得作验证数字」) 已在域外"))


def check_flux_conservation(case: Optional[Mapping[str, Any]]) -> dict:
    """P-4：每列 ``∫2H·w̄ dφ`` 必须为同一常数（(2.2) 单位通量归一 + 差分-梯形恒等式）。

    ⚠️ **这是算子恒等式，不是守恒律的独立测量**（评审 Fix 2，2026-09-26）：
    ``w̄ = ∂φΨ/(2 r_a H)``（``stream_function.py:548``）⇒ ``2H·w̄ = ∂φΨ/r_a``，而
    "中心差分 + 梯形"这一对离散算子对**逐列**严格等于 ``(Ψ(1)−Ψ(0))/r_a``
    （``stream_function.py:521-522`` 的差分-梯形恒等式），同时求解器在**每一列**
    都同样施加 Ψ(φ=0)=0、Ψ(φ=1)=1（``stream_function.py:316-319`` 与
    ``_assemble_banded_interior`` 的 Dirichlet 行）。⇒ 对**任何**内部解，本行算出的
    逐列常数**恒成立**，**不能**在求解器输出上失败；它检验的只是"φ 两端 Dirichlet 行
    是否被逐列一致地施加"（外加差分/梯形实现是否被改坏），**检测不到**内部通量、
    浮力重分配或冻结门引起的通量误差。反例测试只证明它对**被破坏的边界行**敏感
    （``test_P4_flux_conservation_fails_on_broken_dirichlet_bc``）。
    """
    criterion = (f"每列 ∫2H·w̄ dφ 相对差 ≤ {FLUX_REL_TOL:.0e}"
                 f"（原文：(2.2) 单位通量归一；**限定语**：该式在「中心差分+梯形+每列统一"
                 f" φ-Dirichlet」下是代数恒等式，本行因此只验证 Dirichlet 行被逐列一致施加，"
                 f"不测量内部/浮力/冻结门的通量误差）")
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
                      "**本行是算子恒等式、不是守恒测量**：w̄=∂φΨ/(2r_aH) ⇒ 2H·w̄=∂φΨ/r_a，"
                      "中心差分+梯形逐列严格给出 (Ψ(1)−Ψ(0))/r_a，而每列 Dirichlet 值相同 ⇒ "
                      "任何内部解都恒满足本判据，**无法在求解器输出上失败**",
                      "因此本行**检测不到**内部通量分布、浮力重分配（``b_field``）或 "
                      "``wall`` 冻结门引起的通量误差——那些量在本测点层级不可见",
                      "反例测试 ``test_P4_flux_conservation_fails_on_broken_dirichlet_bc`` "
                      "只证明它对「被破坏的 φ 端 Dirichlet 行」敏感",
                      "浮力场为零 ⇒ 本算例亦不覆盖浮力路径"))


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
    """P-6：**完全冻结格**（``wall ≥ 1 − 1e-6``，即 τw = 0）处 ``|w| ≤ 0.01·max|w|``。

    判据的判定部分只落在"完全冻结"格上——这是计划 P-6 的 "冻结区" 前提**字面成立**之处。
    ``wall`` 是**连续**冻结度 ``clip(1−τw_extrap/(f·τy),0,1)``（自 2026-09-15 提交
    ``9448572``「屈服门二值→连续，消除阈值悬崖」起），故 ``wall>0.5`` 只是**软**阈值：
    把 ``0.5<wall<1`` 的过渡带格判成"必须静止"，等于要求模型在被刻意连续化的地方按二值
    行为。过渡带因此作为**披露量**（格数、最高 ``|w|/max|w|``、该格 ``wall``）折进同一行，
    不进判定。

    ``wall=None``（屈服门未启用）、无完全冻结格（**非空过守卫**）⇒ ``通过=None`` 未测。
    附带报告完全冻结格相对**同算例未加冻结门**时同格速度的抑制倍数（内部静止的直接证据）。
    """
    criterion = (f"完全冻结格（wall ≥ {FROZEN_WALL_FULLY:.6g}，即 τw=0）处 "
                 f"|w| ≤ {FROZEN_SPEED_FRACTION_MAX}·max|w|")
    if case is None:
        return _not_measured(CHECK_FROZEN_REGION, criterion,
                             "未提供求解器级合成算例；结果对象不导出速度场")
    wall = case.get("wall")
    if wall is None:
        return _not_measured(CHECK_FROZEN_REGION, criterion,
                             "合成算例未启用屈服门（wall=None，等价 "
                             "enable_stream_yield_gate=False）⇒ 无冻结单元")
    wall = np.asarray(wall, dtype=float)
    fully = wall >= FROZEN_WALL_FULLY
    band = (wall > FROZEN_WALL_THRESHOLD) & ~fully
    n_fully, n_band = int(np.count_nonzero(fully)), int(np.count_nonzero(band))
    w = np.abs(np.asarray(case["w_m_s"], dtype=float))
    peak = float(np.max(w))
    if peak <= 0.0:
        return _not_measured(CHECK_FROZEN_REGION, criterion,
                             "速度场恒为零（无流动），判据退化无判别力")
    band_ratio, band_wall = _band_stats(w, wall, band, peak)
    if n_fully == 0:
        # 非空过守卫：没有任何格满足"完全冻结"前提 ⇒ 未测，绝不空过
        return _not_measured(
            CHECK_FROZEN_REGION, criterion,
            f"合成算例无完全冻结格（wall ≥ {FROZEN_WALL_FULLY:.6g}）："
            f"wall>{FROZEN_WALL_THRESHOLD} 的过渡带 {n_band} 格最高 "
            f"|w|/max|w|={band_ratio:.3e}，但过渡带不是「必须静止」的判定对象 ⇒ 无从判定",
            measured=(f"完全冻结格=0；过渡带 {n_band} 格最高比值={band_ratio:.3e}"
                      f"（该格 wall={band_wall:.4g}）"))
    frozen_peak = float(np.max(w[fully]))
    ratio = frozen_peak / peak
    passed = ratio <= FROZEN_SPEED_FRACTION_MAX
    suppressed = _suppression(case, w, fully)
    measured = (f"完全冻结格 {n_fully} 个：冻结格内 max|w|={frozen_peak:.3e} m/s, "
                f"全域 max|w|={peak:.6g} m/s, 比值={ratio:.3e}；"
                f"过渡带（{FROZEN_WALL_THRESHOLD}<wall<1−1e-6）{n_band} 格："
                f"最高 |w|/max|w|={band_ratio:.3e}（该格 wall={band_wall:.4g}）")
    return _row(CHECK_FROZEN_REGION, passed, measured, criterion,
                _note(LEVEL_SYNTHETIC, str(case.get("case", "")),
                      f"wall 是**连续**冻结度 clip(1−τw_extrap/(f·τy),0,1)"
                      f"（自提交 9448572「屈服门二值→连续，消除阈值悬崖」起）；"
                      f"wall>{FROZEN_WALL_THRESHOLD} 是**软**阈值，"
                      f"故判定只取 wall≥{FROZEN_WALL_FULLY:.6g}（τw=0）的完全冻结格，"
                      f"过渡带仅披露——把过渡带判成「必须静止」会要求模型"
                      f"在被刻意连续化的地方按二值行为",
                      f"完全冻结格相对同算例未加冻结门时同格速度抑制倍数={suppressed:.3e}"
                      f"（内部静止的直接证据）" if suppressed else "",
                      "**本行是合成场上的算子性质，不是逐井 production 判决**：production 里 "
                      "τw_extrap = τw_ref·b/b_ref 且 b>0"
                      "（``annulus_d2dga.py`` 的 ``_yield_gate_wall`` 体内：``gamma``/``tau_w_field``"
                      "→``tau_w_ref``→``tau_w_extrap``，即 ``:1070-1083``；"
                      "最终 ``wall_new = np.clip(wall_raw, 0.0, 1.0)`` 在 ``:1096``），"
                      "故 wall=1（τw=0）要求该列**流得最快**的参考格静止；而参考格按定义取"
                      "「正在流动」的格 ⇒ production 中 wall≡1 只在 ``col_freeze`` 分支出现"
                      "（``~has_flow & cement_ever>0``、随后 ``wall_new[:, col_freeze] = 1.0``，"
                      "即 ``annulus_d2dga.py:1102-1104``；论证真正依赖的 ``ref_mask`` 置零行在"
                      "``:1101``）。注意 ``:1101`` 是 ``col_freeze`` 之前的 ref_mask 行，"
                      "``col_freeze`` 本身在 ``:1103-1104``。即整列本就静止。"
                      "本行证明的是**算子把 wall≡1 处速度压到零**（该性质由本算例独立证实），"
                      "而非任何一口井在 production 下的冻结行为",
                      f"**判决对冻结度剖面形状敏感**（披露）：完全冻结格相对阈值的余量是 "
                      f"{FROZEN_SPEED_FRACTION_MAX / ratio:.1f}×，但该余量依赖冻结区的"
                      f"φ-展布——本行剖面为 C¹（sin²）过渡。若改用**分段线性**过渡（冻结度在"
                      f"一格内跨完），第一个完全冻结格会因中心差分跨台阶而泄漏到 1.10e-2 "
                      f"（实测，ramp 25→35、ny=41）> {FROZEN_SPEED_FRACTION_MAX}，本行会转红；"
                      f"那是**离散格式在冻结度突变处的局限**，不是冻结区内部在流动",
                      "同源参照：二值掩码下（Task 9 首版算例）冻结前锋格比值 9.38e-2（旧报告节，"
                      "已被 §11 取代）"))


def _band_stats(w: np.ndarray, wall: np.ndarray, band: np.ndarray, peak: float) -> tuple:
    """过渡带（0.5<wall<1−1e-6）的最高 |w|/max|w| 与该格自身的 wall 值。"""
    if not bool(np.any(band)):
        return 0.0, float("nan")
    flat = np.where(band.ravel(), w.ravel(), -1.0)
    k = int(np.argmax(flat))
    return float(flat[k] / peak), float(wall.ravel()[k])


def _suppression(case: Mapping[str, Any], w: np.ndarray, fully: np.ndarray) -> float:
    """完全冻结格在"加冻结门"与"未加冻结门"两种解下的速度之比（越小越静止）。

    需算例提供重解所需的参数（``c_bar``/``eta1``/…）；缺失或异常时返回 0.0（不披露）。
    """
    if not all(k in case for k in ("c_bar", "eta1", "eta2", "m")):
        return 0.0
    try:
        ref = np.abs(case_velocity(case, None))
    except Exception:  # noqa: BLE001 —— 纯披露量，重解失败不影响判定
        return 0.0
    ref_peak = float(np.max(ref[fully])) if bool(np.any(fully)) else 0.0
    if ref_peak <= 0.0:
        return 0.0
    return float(np.max(w[fully]) / ref_peak)


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

    7 项中 P-4 为**算子恒等式**（非独立测量，见模块 docstring 的 ⚠️ 计数口径）。

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
