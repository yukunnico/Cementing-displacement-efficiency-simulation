"""
顶替指标诊断层 (displacement_metrics)
=====================================

本模块从环空二维顶替模拟结果（AnnulusSimulationResult）中提取文献标准的
顶替质量指标，覆盖 Tier 0 路线图的 T0-4 / T0-5 / T0-7 三项
（见 docs/superpowers/specs/2026-07-17-d2dga-full-roadmap-tier0-design.md §2.4/2.5/2.7）。

指标与论文来源
--------------
T0-4 泥浆滞留 + φm/φc 污染指示器（Yang et al. 2021, §4.4 "Analysis of the leftover mud"）
    - mud_retention_fraction = ∫b·mud dA / ∫b dA（最终时刻泥浆滞留体积分数）
    - φm/φc = mud / max(cement, ε) 逐深度比值，即泥浆污染指示器
      （Yang 2021 明确该指示器源自 Zulqarnain and Tyagi, 2016），分级阈值：
      green (φm/φc ≤ 0.05) / yellow (0.05 < φm/φc < 0.1) / red (φm/φc ≥ 0.1)

T0-5 界面长度比 + 窄边效率 η_N（Yang 2021 / Zhang & Frigaard 2023 图 19）
    - interface_length_ratio = |front_wide − front_narrow| / max(mean(fronts), ε)
      偏心环空宽窄边前缘分离程度的归一化度量，0 表示活塞式齐头推进
    - η_N = 窄四分位（方位角最后 1/4 行）的 ∫b·cement/∫b；
      窄边是偏心环空顶替最差区域，物理上应有 η_N ≤ η_E

T0-7 突破时间 t_br（Zhang & Frigaard 2022, Table 3）
    - t_br_s = 水泥前缘（求解器按 c ≥ 0.5 阈值定义）首次到达出口 s_max 的时间
    - 无量纲化约定（Zhang 2022 原文：t_br = t̂_br × ŵ₀/L̂_annu，Table 3 范围 0.33–0.99）：
      t_br_hat = t_br_s · ŵ₀ / L，其中 ŵ₀ 为平均泵速 [m/s]，
      L = geom["s"][-1] − geom["s"][0] 为井段长度 [m]

设计原则（与 quality_proxy.py 一致）
------------------------------------
1. 纯后处理诊断：只消费 AnnulusSimulationResult，不 import models2d 求解逻辑。
2. 数值安全：ε 除零保护、未突破 t_br = inf、窄四分位单行边界退化保护。
3. 约定显式化：无量纲化与 ŵ₀ 估计约定写入 docstring 与结果 notes。

使用示例
--------
>>> from cemdisp.diagnostics.displacement_metrics import compute_displacement_metrics
>>> metrics_result = compute_displacement_metrics(result)  # result 为 AnnulusSimulationResult
>>> print(metrics_result.mud_retention_fraction)
>>> print(metrics_result.quality_zone_counts)
>>> print(metrics_result.t_br_hat)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Dict, Optional

import numpy as np

if TYPE_CHECKING:
    # 仅用于类型注解；运行时不 import models2d，保持诊断层与求解器解耦。
    from cemdisp.models2d.annulus_d2dga import AnnulusSimulationResult


# Yang 2021 §4.4 / Zulqarnain and Tyagi (2016) 的 φm/φc 分级阈值
PHI_GREEN_MAX = 0.05  # φm/φc ≤ 0.05 → green（胶结质量良好）
PHI_RED_MIN = 0.10    # φm/φc ≥ 0.10 → red（泥浆污染严重）；严格介于两者之间 → yellow

_EPS = 1.0e-12  # 全局数值安全阈值


@dataclass(frozen=True)
class DisplacementMetricsResult:
    """顶替指标诊断结果（T0-4 / T0-5 / T0-7 汇总）。

    Attributes:
        mud_retention_fraction: 最终时刻泥浆滞留体积分数 ∫b·mud/∫b ∈ [0, 1]（Yang 2021）
        phi_m_phi_c_profile: 逐深度 φm/φc 比值元组，行序与 depth_profiles 一致（井底→井顶）
        quality_zone_counts: φm/φc 分级段数统计 {"green": …, "yellow": …, "red": …}
        interface_length_ratio: 宽窄边界面长度比 |fw − fn| / mean(fw, fn)（Yang 2021）
        eta_narrow: 窄四分位效率 η_N = ∫b·cement/∫b（Zhang 2023 图 19 口径）
        eta_global: 全局有效顶替效率 η_E（metrics 末行 effective_efficiency）
        t_br_s: 突破时间 [s]；模拟结束时未突破为 inf
        t_br_hat: 无量纲突破时间 t_br·ŵ₀/L；未突破为 inf，ŵ₀ 不可估计时为 nan
        w0_m_s: 实际采用的平均泵速 ŵ₀ [m/s]（外部给定或自动估计值；不可估计为 nan）
        notes: 约定与警告说明列表
    """

    mud_retention_fraction: float
    phi_m_phi_c_profile: tuple[float, ...]
    quality_zone_counts: Dict[str, int]
    interface_length_ratio: float
    eta_narrow: float
    eta_global: float
    t_br_s: float
    t_br_hat: float
    w0_m_s: float
    notes: tuple[str, ...] = ()

    def to_dict(self) -> Dict[str, object]:
        """转为可 JSON 序列化的字典（tuple → list，numpy 标量 → float）。"""
        return {
            "mud_retention_fraction": float(self.mud_retention_fraction),
            "phi_m_phi_c_profile": [float(v) for v in self.phi_m_phi_c_profile],
            "quality_zone_counts": dict(self.quality_zone_counts),
            "interface_length_ratio": float(self.interface_length_ratio),
            "eta_narrow": float(self.eta_narrow),
            "eta_global": float(self.eta_global),
            "t_br_s": float(self.t_br_s),
            "t_br_hat": float(self.t_br_hat),
            "w0_m_s": float(self.w0_m_s),
            "notes": list(self.notes),
        }


def _trapz2d(arr: np.ndarray, y: np.ndarray, s: np.ndarray) -> float:
    """梯形法则二维积分 ∫∫arr dy ds，与 annulus_d2dga._trapez2d 同口径（先 s 后 y）。

    对单行/单列子区域做退化保护（np.trapezoid 单点返回 0，会导致比值 0/0 = nan）：
    单列时退化为方位向积分，单行时退化为轴向积分。
    """
    arr = np.asarray(arr, dtype=float)
    if arr.shape[1] == 1:
        axial = arr[:, 0]
    else:
        axial = np.trapezoid(arr, x=s, axis=1)
    if axial.shape[0] == 1:
        return float(axial[0])
    return float(np.trapezoid(axial, x=y, axis=0))


def _b_weighted_azimuth_mean(field: np.ndarray, b: np.ndarray, y: np.ndarray, eps: float) -> np.ndarray:
    """方位向 b 加权平均 f̄(z) = ∫b·f dy / ∫b dy，返回 (nz,) 数组。

    等价于 depth_profiles 中 "钻井液平均浓度" / "水泥平均浓度" 列的计算口径
    （np.average(field, axis=0, weights=b) 的梯形积分版本）。ny == 1 时退化为该行取值。
    """
    if b.shape[0] == 1:
        den = np.maximum(b[0], eps)
        return (b[0] * field[0]) / den
    den = np.maximum(np.trapezoid(b, x=y, axis=0), eps)
    return np.trapezoid(b * field, x=y, axis=0) / den


def _final_mud_field(result: "AnnulusSimulationResult") -> np.ndarray:
    """最终时刻泥浆体积分数场 mud = clip(1 − cement − spacer, 0, 1)。

    与 annulus_d2dga._depth_profiles 的钻井液定义一致（体积分数闭合反算：
    显式跟踪领浆/尾浆/前置液，钻井液 = 1 − 三相之和）。
    """
    cement = np.clip(np.asarray(result.cement_field, dtype=float), 0.0, 1.0)
    spacer = np.clip(np.asarray(result.spacer_field, dtype=float), 0.0, 1.0)
    return np.clip(1.0 - cement - spacer, 0.0, 1.0)


def _mud_retention_fraction(mud: np.ndarray, geom: Dict[str, np.ndarray]) -> float:
    """泥浆滞留体积分数 = ∫b·mud dA / ∫b dA（Yang 2021 §4.4 mud retention volume）。

    半环空展开积分的因子 2 在比值中约去，结果 ∈ [0, 1]。
    """
    b = np.asarray(geom["b"], dtype=float)
    num = _trapz2d(b * mud, geom["y"], geom["s"])
    den = _trapz2d(b, geom["y"], geom["s"])
    return float(np.clip(num / max(den, _EPS), 0.0, 1.0))


def _phi_m_phi_c_profile(mud: np.ndarray, cement: np.ndarray, geom: Dict[str, np.ndarray], eps: float) -> tuple[float, ...]:
    """逐深度 φm/φc = mud̄(z) / max(cement̄(z), ε)（Yang 2021 §4.4 泥浆污染指示器）。

    mud̄/cement̄ 为方位向 b 加权平均；cement̄ ≈ 0（纯泥浆段）时比值由 ε 保护，
    物理含义为该深度几乎无水泥、污染最严重（判 red）。
    """
    b = np.asarray(geom["b"], dtype=float)
    mud_bar = _b_weighted_azimuth_mean(mud, b, geom["y"], eps)
    cement_bar = _b_weighted_azimuth_mean(cement, b, geom["y"], eps)
    ratio = mud_bar / np.maximum(cement_bar, eps)
    return tuple(float(v) for v in ratio)


def _classify_quality_zones(
    profile: tuple[float, ...],
    green_max: float = PHI_GREEN_MAX,
    red_min: float = PHI_RED_MIN,
) -> Dict[str, int]:
    """按 Yang 2021 §4.4 阈值对逐深度 φm/φc 分级并计数。

    green: φm/φc ≤ 0.05；yellow: 0.05 < φm/φc < 0.1；red: φm/φc ≥ 0.1。
    边界值 0.05 归 green、0.1 归 red，与论文原文不等式方向一致。
    """
    counts = {"green": 0, "yellow": 0, "red": 0}
    for value in profile:
        if value >= red_min:
            counts["red"] += 1
        elif value > green_max:
            counts["yellow"] += 1
        else:
            counts["green"] += 1
    return counts


def _interface_length_ratio(front_wide_m: float, front_narrow_m: float, eps: float) -> float:
    """界面长度比 = |front_wide − front_narrow| / max(mean(fronts), ε)（Yang 2021 核心观察量）。

    0 表示宽窄边前缘齐头推进（活塞式）；越大表示前缘分离越严重。
    两者均为 0（水泥尚未入环空）时由 ε 保护返回 0。
    """
    mean_front = 0.5 * (float(front_wide_m) + float(front_narrow_m))
    return float(abs(float(front_wide_m) - float(front_narrow_m)) / max(mean_front, eps))


def _narrow_quarter_efficiency(cement: np.ndarray, geom: Dict[str, np.ndarray]) -> float:
    """窄边效率 η_N：方位角最后 1/4 行（窄四分位）的 ∫b·cement/∫b（Zhang 2023 图 19 口径）。

    网格约定：geom["y"] 第 0 行为宽边、最后一行为窄边（h = H(1 + e·cos(πφ))，
    φ = 1 处间隙最小），与求解器 front_narrow = cement[-1] 的约定一致。
    行数 n_q = max(1, ny // 4)；ny < 4 时退化为最后 1 行，由 _trapz2d 单行保护。
    """
    ny = int(cement.shape[0])
    n_q = max(1, ny // 4)
    b_q = np.asarray(geom["b"], dtype=float)[-n_q:, :]
    c_q = cement[-n_q:, :]
    y_q = np.asarray(geom["y"], dtype=float)[-n_q:]
    num = _trapz2d(b_q * c_q, y_q, geom["s"])
    den = _trapz2d(b_q, y_q, geom["s"])
    return float(np.clip(num / max(den, _EPS), 0.0, 1.0))


def _breakthrough_time_s(metrics, s_max: float, eps: float) -> float:
    """突破时间 t_br：水泥前缘首次到达出口 s_max 的时间（Zhang 2022 Table 3）。

    前缘位置取 metrics 的 front_wide_m / front_narrow_m 两列
    （求解器按 c ≥ 0.5 阈值沿宽/窄边方位线定义），任一首次 ≥ s_max 即视为突破，
    取两者较早者；均未突破返回 inf。
    """
    t = metrics["time_s"].to_numpy(dtype=float)
    candidates: list[float] = []
    for col in ("front_wide_m", "front_narrow_m"):
        front = metrics[col].to_numpy(dtype=float)
        hit = np.where(front >= s_max - eps)[0]
        if hit.size:
            candidates.append(float(t[hit[0]]))
    return min(candidates) if candidates else float("inf")


def _estimate_w0_from_bulk_fill(metrics, L: float, t_br_s: float, eps: float) -> float:
    """由 bulk_cement_fill 上升段斜率估计平均泵速 ŵ₀ = L · d(fill)/dt。

    物理依据：突破前且纯水泥入井时段，环空内水泥体积以泵排量 Q 线性增长，
    d(fill)/dt = Q/V_annulus = ŵ₀/L（不可压缩、截面平均意义下精确，
    与几何无关）。取逐步正斜率的中位数，对注隔离液平台段（斜率 ≈ 0）
    与数值噪声稳健。拟合窗口限制在 t ≤ t_br（突破后斜率因水泥流出而偏低）。
    估计失败（列缺失 / 无上升段 / L ≤ 0）返回 nan。
    """
    if "bulk_cement_fill" not in metrics.columns or "time_s" not in metrics.columns:
        return float("nan")
    t = metrics["time_s"].to_numpy(dtype=float)
    fill = metrics["bulk_cement_fill"].to_numpy(dtype=float)
    if np.isfinite(t_br_s):
        mask = t <= t_br_s + eps
        if int(mask.sum()) >= 2:
            t, fill = t[mask], fill[mask]
    if t.size < 2 or L <= 0.0:
        return float("nan")
    dt = np.diff(t)
    dfill = np.diff(fill)
    rising = (dt > 0.0) & (dfill > eps)
    if not np.any(rising):
        return float("nan")
    slope = float(np.median(dfill[rising] / dt[rising]))
    return float(L * slope)


def compute_displacement_metrics(
    result: "AnnulusSimulationResult",
    *,
    w0_m_s: Optional[float] = None,
    eps: float = _EPS,
) -> DisplacementMetricsResult:
    """从环空二维顶替结果计算文献标准顶替指标（T0-4 / T0-5 / T0-7）。

    Args:
        result: AnnulusSimulationResult（或具有 geom/cement_field/spacer_field/metrics
            同名属性的对象，便于单元测试 mock）
        w0_m_s: 平均泵速 ŵ₀ [m/s]，用于 t_br 无量纲化；None 时按
            ŵ₀ = L · d(bulk_cement_fill)/dt 从 metrics 上升段自动估计
            （约定见 _estimate_w0_from_bulk_fill docstring）
        eps: 数值安全阈值（除零保护、突破判据容差）

    Returns:
        DisplacementMetricsResult，含泥浆滞留、φm/φc 分级、界面长度比、
        η_N/η_E、突破时间及其无量纲化

    无量纲化约定（Zhang & Frigaard 2022 Table 3 原文定义）：
        t_br_hat = t_br_s · ŵ₀ / L，L = geom["s"][-1] − geom["s"][0]。
        未突破时 t_br_s = t_br_hat = inf；ŵ₀ 不可估计且未外部给定时 t_br_hat = nan。
    """
    geom = result.geom
    s = np.asarray(geom["s"], dtype=float)
    y = np.asarray(geom["y"], dtype=float)
    b = np.asarray(geom["b"], dtype=float)
    cement = np.clip(np.asarray(result.cement_field, dtype=float), 0.0, 1.0)
    mud = _final_mud_field(result)
    metrics = result.metrics
    notes: list[str] = []

    # ---- T0-4：泥浆滞留 + φm/φc 污染指示器（Yang 2021 §4.4）----
    mud_retention = _mud_retention_fraction(mud, geom)
    phi_profile = _phi_m_phi_c_profile(mud, cement, geom, eps)
    zone_counts = _classify_quality_zones(phi_profile)

    # ---- T0-5：界面长度比 + 窄边效率 η_N（Yang 2021 / Zhang 2023 图 19）----
    has_fronts = len(metrics) > 0 and {"front_wide_m", "front_narrow_m"} <= set(metrics.columns)
    if has_fronts:
        final_row = metrics.iloc[-1]
        front_wide = float(final_row["front_wide_m"])
        front_narrow = float(final_row["front_narrow_m"])
    else:
        front_wide = front_narrow = 0.0
        notes.append("metrics 缺少 front_wide_m/front_narrow_m 列或为空，界面长度比按 0 处理。")
    interface_ratio = _interface_length_ratio(front_wide, front_narrow, eps)

    eta_narrow = _narrow_quarter_efficiency(cement, geom)
    if len(metrics) > 0 and "effective_efficiency" in metrics.columns:
        eta_global = float(metrics["effective_efficiency"].iloc[-1])
    else:
        den = _trapz2d(b, y, s)
        eta_global = float(np.clip(_trapz2d(b * cement, y, s) / max(den, eps), 0.0, 1.0))
        notes.append("metrics 缺少 effective_efficiency 列，η_E 由最终水泥场 ∫b·cement/∫b 反算。")
    if eta_narrow - eta_global > 0.02:
        notes.append(
            f"η_N ({eta_narrow:.4f}) > η_E ({eta_global:.4f})，物理上少见"
            "（偏心环空窄边通常顶替最差），请检查方位角方向约定或场数据。"
        )

    # ---- T0-7：突破时间 t_br（Zhang 2022 Table 3）----
    s_max = float(s[-1])
    L = float(s[-1] - s[0])
    if len(metrics) > 0 and {"time_s", "front_wide_m", "front_narrow_m"} <= set(metrics.columns):
        t_br = _breakthrough_time_s(metrics, s_max, eps)
    else:
        t_br = float("inf")

    # ŵ₀：外部给定优先，否则由 bulk_cement_fill 上升段斜率估计
    if w0_m_s is not None:
        w0 = float(w0_m_s)
        w0_source = "外部给定"
    elif len(metrics) > 0:
        w0 = _estimate_w0_from_bulk_fill(metrics, L, t_br, eps)
        w0_source = "bulk_cement_fill 上升段斜率估计（ŵ₀ = L·d(fill)/dt）"
    else:
        w0 = float("nan")
        w0_source = "不可用（metrics 为空）"

    if not np.isfinite(t_br):
        t_br_hat = float("inf")
        notes.append("模拟结束时水泥前缘未到达出口 s_max，t_br = inf（未突破）。")
    elif np.isfinite(w0) and w0 > 0.0 and L > 0.0:
        t_br_hat = float(t_br * w0 / L)
        notes.append(
            f"ŵ₀ = {w0:.6g} m/s（{w0_source}）；"
            f"t̂_br = t_br·ŵ₀/L（Zhang 2022 Table 3 约定），L = {L:.4g} m。"
        )
    else:
        t_br_hat = float("nan")
        notes.append("ŵ₀ 不可估计且未外部给定，t_br_hat = nan；可通过 w0_m_s 参数传入平均泵速。")

    return DisplacementMetricsResult(
        mud_retention_fraction=round(mud_retention, 6),
        phi_m_phi_c_profile=phi_profile,
        quality_zone_counts=zone_counts,
        interface_length_ratio=round(interface_ratio, 6),
        eta_narrow=round(eta_narrow, 6),
        eta_global=round(eta_global, 6),
        t_br_s=t_br,
        t_br_hat=t_br_hat,
        w0_m_s=w0,
        notes=tuple(notes),
    )
