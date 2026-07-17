"""
Regime 分类诊断 (regime_classifiers)
=====================================

本模块提供两个纯后处理诊断（Tier 0，Layer 5），只消费模拟结果与输入数据，
不修改任何求解器状态：

1. **T0-3 浮力数阈值分类**（Zhang & Frigaard 2023, §3.1.3）：
   浮力数 b = Δρ·g·d²/(μ·w₀) 是垂直井顶替的主导无量纲参数，按阈值将工况分为
   forbidden / highly_dispersive / steady_capable / non_dispersive 四类，
   并给出中文设计建议。

2. **T0-6 停泵有限时间衰减诊断**（Moyers-González et al. 2007, 式 3.35/3.40）：
   停泵（SHUTDOWN）后无外加流量（Q=0），若屈服应力足以克服浮力散度驱动的
   有效应力，速度场将在**有限时间**内衰减至零、界面冻结；
   否则停泵期间界面持续密度驱动运移，存在窜流风险。

设计原则
--------
1. 纯诊断：只消费 ``AnnulusSimulationResult`` / ``FluidSpec`` / ``PumpingSchedule``，
   不 import models2d 求解逻辑（类型注解走 TYPE_CHECKING，运行时可接受
   SimpleNamespace 等鸭子类型 mock）。
2. 公式标注论文来源与式号；所有近似在 docstring 中显式声明。
3. 数值安全：除零保护、NaN 守卫、τ_Y,min=0（牛顿流体）时判不满足
   （无屈服应力支撑则速度仅指数衰减，不会有限时间归零）。

近似假设（T0-6，论文为无量纲形式，本模块为有量纲重建）
------------------------------------------------------
* **浮力散度等效应力**（式 3.35 右端 ||∇·f||_∞/2 的有量纲化）：
  论文浮力向量 f = (r_a·cosβ/F², r_a·sin(πφ)·sinβ/F²)（式 2.5b，无量纲），
  其散度 ||∇·f||_∞ = π·r_a·sinβ/F²。以浮力数 b = Δρ·g·d²/(μ·w₀) 表示 1/F²
  并乘以应力尺度 μ·w₀/d，得有量纲等效应力

      S_buoy = π·r_a·Δρ·g·sinβ / 2   [Pa]

  物理意义：密度差 × 重力 × 环空平均半径 = 倾斜界面静水压差的特征尺度。
  垂直井（β=0）时 S_buoy=0，判据在 τ_Y,min>0 时自动满足，与物理一致
  （垂直井停泵后无方位浮力驱动，屈服应力必使界面冻结）。

* **冻结时间**（式 3.40 的有量纲重建）：

      α₁ = μ_ref·(1+e)/(ρ_min·d²)      [1/s]   粘性衰减率
      α₂ = 2(1+e)·Δτ/(ρ_min·π²·d)      [m/s²]  净屈服应力剩余驱动的减速度
      z₀ = w₀·(1+e)/π                  [m/s]   停泵前速度范数（||Ψ||_U(0)≈1 归一化）
      t_s = (1/α₁)·ln(1 + z₀·α₁/α₂)    [s]

  其中 Δτ = τ_Y,min/(1+e) − S_buoy > 0，d 为半间隙，Poincaré 常数 C_Ω=π²
  与 B_Ω=2（论文式 3.37 下）已按论文推导吸入相应定义。
  该重建在 Bingham 流体同心环空极限（e=0, β=0）下与论文式 3.40 后的
  显式特例 t_s = (ρ/κ)·ln(1 + π·κ·||Ψ||_U(0)/(2τ_Y)) 精确一致。

* **停泵前速度 w₀**：由第一个 SHUTDOWN 事件之前最后一个正排量步骤的排量
  除以环空截面积 A = π/4·(D_hole² − D_od²) 估计；缺井径数据时退化为
  典型顶替速度 0.5 m/s；停泵前无任何正排量步骤时 w₀=0（已静止，t_s=0）。

使用示例
--------
>>> from cemdisp.diagnostics.regime_classifiers import (
...     classify_buoyancy_regime, compute_shutdown_decay,
... )
>>> res = classify_buoyancy_regime(50.0)
>>> res.regime
'steady_capable'
>>> decay = compute_shutdown_decay(result, fluids, schedule)
>>> decay.condition_satisfied, decay.freeze_time_s
(True, 4.13)
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Dict, Optional, Sequence, Tuple

import numpy as np

if TYPE_CHECKING:
    from cemdisp.data.fluid_spec import FluidSpec
    from cemdisp.data.pumping_schedule import PumpingSchedule
    from cemdisp.models2d.annulus_d2dga import AnnulusSimulationResult


# ---------------------------------------------------------------------------
# T0-3 浮力数阈值分类（Zhang & Frigaard 2023, §3.1.3）
# ---------------------------------------------------------------------------

#: steady_capable 的 b 下限（0 ≤ b < 20 为 highly_dispersive）
B_STEADY_CAPABLE_MIN: float = 20.0
#: non_dispersive 的 b 下限（20 ≤ b < 80 为 steady_capable）
B_NON_DISPERSIVE_MIN: float = 80.0

_REGIME_FORBIDDEN = "forbidden"
_REGIME_HIGHLY_DISPERSIVE = "highly_dispersive"
_REGIME_STEADY_CAPABLE = "steady_capable"
_REGIME_NON_DISPERSIVE = "non_dispersive"


@dataclass(frozen=True)
class BuoyancyRegimeResult:
    """浮力数阈值分类结果（T0-3）。

    Attributes:
        b_number: 输入浮力数 b = Δρ·g·d²/(μ·w₀)
        regime: 分类标签，四选一：
            "forbidden"（b<0，密度不稳定，禁止施工）、
            "highly_dispersive"（0≤b<20，高度弥散）、
            "steady_capable"（20≤b<80，steady front 可达）、
            "non_dispersive"（b≥80，非弥散）
        design_advice: 中文设计建议
    """

    b_number: float
    regime: str
    design_advice: str

    def to_dict(self) -> Dict[str, Any]:
        """转为字典（可 JSON 序列化）。"""
        return {
            "b_number": self.b_number,
            "regime": self.regime,
            "design_advice": self.design_advice,
        }


def classify_buoyancy_regime(b_number: float) -> BuoyancyRegimeResult:
    """按 Zhang & Frigaard (2023) §3.1.3 的阈值对浮力数 b 分类。

    判据：
        b < 0        → forbidden（密度不稳定：轻浆顶替重浆，禁止施工）
        0 ≤ b < 20   → highly_dispersive（浮力稳定作用不足，界面高度弥散）
        20 ≤ b < 80  → steady_capable（浮力足以形成稳定前缘）
        b ≥ 80       → non_dispersive（浮力主导，界面非弥散，近活塞式顶替）

    Args:
        b_number: 浮力数（通常从 ``result.summary['buoyancy_number']`` 读取，
            由 ``annulus_d2dga.py`` 的 ``_compute_buoyancy_number`` 计算）

    Returns:
        BuoyancyRegimeResult

    Raises:
        ValueError: b_number 为 NaN 时
    """
    b = float(b_number)
    if math.isnan(b):
        raise ValueError(f"b_number 不能为 NaN，当前值：{b_number}")

    if b < 0.0:
        regime = _REGIME_FORBIDDEN
        advice = (
            "浮力数 b<0：密度差方向错误（轻浆顶替重浆），界面密度不稳定，禁止施工；"
            "必须保证重浆顶替轻浆（b≥0），请调整浆体密度序列。"
        )
    elif b < B_STEADY_CAPABLE_MIN:
        regime = _REGIME_HIGHLY_DISPERSIVE
        advice = (
            "0≤b<20：浮力稳定作用不足，界面高度弥散，顶替效率低；"
            "建议提高顶替液与被顶替液的密度差，或降低排量/泥浆粘度以增大 b。"
        )
    elif b < B_NON_DISPERSIVE_MIN:
        regime = _REGIME_STEADY_CAPABLE
        advice = (
            "20≤b<80：浮力足以形成稳定前缘（steady front 可达），但仍存在弥散；"
            "工况可接受，建议结合流动分类诊断（T0-1）确认前缘稳定性。"
        )
    else:
        regime = _REGIME_NON_DISPERSIVE
        advice = (
            "b≥80：浮力主导，界面非弥散，接近活塞式顶替；"
            "注意复核窄边泥浆滞留与失稳指数，避免局部残留。"
        )

    return BuoyancyRegimeResult(b_number=b, regime=regime, design_advice=advice)


def classify_buoyancy_regime_from_result(result: "AnnulusSimulationResult") -> BuoyancyRegimeResult:
    """从模拟结果 summary 读取浮力数并分类（便捷包装）。

    Args:
        result: AnnulusSimulationResult，其 summary 含 'buoyancy_number' 键

    Returns:
        BuoyancyRegimeResult

    Raises:
        KeyError: summary 中无 'buoyancy_number' 键时
    """
    return classify_buoyancy_regime(float(result.summary["buoyancy_number"]))


# ---------------------------------------------------------------------------
# T0-6 停泵有限时间衰减（Moyers-González et al. 2007, 式 3.35/3.40）
# ---------------------------------------------------------------------------

_GRAVITY_M_S2: float = 9.81
#: 参考塑性粘度 fallback（与 annulus_d2dga.py 中 mud 的 None guard 取值一致）
_FALLBACK_VISCOSITY_PA_S: float = 0.05
#: 缺井径数据无法估计停泵前速度时的典型顶替速度 fallback [m/s]
_FALLBACK_VELOCITY_M_S: float = 0.5
#: 缺井径数据时的环空平均半径 fallback [m]（典型 8½in 井眼 × 5½in 套管量级）
_FALLBACK_MEAN_RADIUS_M: float = 0.10


@dataclass(frozen=True)
class ShutdownDecayResult:
    """停泵有限时间衰减诊断结果（T0-6）。

    Attributes:
        condition_satisfied: 式 3.35 判据是否满足（满足则界面有限时间冻结）
        freeze_time_s: 冻结时间 t_s（式 3.40）[s]；判据不满足或诊断不适用时为 inf
        tau_y_min: 各流体最小屈服应力 τ_Y,min [Pa]（无屈服应力字段的流体按 0 计）
        physical_interpretation: 中文物理解释
    """

    condition_satisfied: bool
    freeze_time_s: float
    tau_y_min: float
    physical_interpretation: str

    def to_dict(self) -> Dict[str, Any]:
        """转为字典。注意 freeze_time_s 可能为 float('inf')，
        严格 JSON 序列化时需调用方自行处理。"""
        return {
            "condition_satisfied": self.condition_satisfied,
            "freeze_time_s": self.freeze_time_s,
            "tau_y_min": self.tau_y_min,
            "physical_interpretation": self.physical_interpretation,
        }


def _event_tag_value(tag: object) -> Optional[str]:
    """提取 event_tag 的字符串值（兼容 PumpingStageEvent 枚举与裸字符串 mock）。"""
    if tag is None:
        return None
    return str(getattr(tag, "value", tag))


def _has_shutdown_event(schedule: "PumpingSchedule") -> bool:
    """检查泵注程序是否含 SHUTDOWN 事件（PumpingStageEvent.SHUTDOWN）。"""
    return any(
        _event_tag_value(getattr(step, "event_tag", None)) == "SHUTDOWN"
        for step in schedule.steps
    )


def _min_yield_stress(fluids: Sequence["FluidSpec"]) -> float:
    """各流体最小屈服应力 τ_Y,min [Pa]；yield_stress_pa 缺失（牛顿流体）按 0 计。"""
    values = []
    for fluid in fluids:
        value = getattr(fluid, "yield_stress_pa", None)
        values.append(0.0 if value is None else float(value))
    return min(values) if values else 0.0


def _reference_viscosity(fluids: Sequence["FluidSpec"]) -> float:
    """参考塑性粘度 μ_ref [Pa·s]：取各流体正塑性粘度的最小值；
    全部缺失时退化为 0.05 Pa·s（与求解器 mud None guard 一致）。"""
    values = []
    for fluid in fluids:
        value = getattr(fluid, "plastic_viscosity_pa_s", None)
        if value is not None and float(value) > 0.0:
            values.append(float(value))
    return min(values) if values else _FALLBACK_VISCOSITY_PA_S


def _density_extremes(fluids: Sequence["FluidSpec"]) -> Tuple[float, float]:
    """返回 (ρ_min, Δρ)：最小密度与最大密度差 [kg/m³]。

    停泵后密度分层驱动（重下轻上）由最大密度差决定，故 Δρ = max − min ≥ 0。
    """
    densities = [float(fluid.density_kg_m3) for fluid in fluids]
    rho_min = min(densities)
    return rho_min, max(densities) - rho_min


def _mean_eccentricity(geom: Dict[str, Any]) -> float:
    """平均偏心度 e（geom['e'] 沿深度平均）。"""
    return float(np.mean(np.asarray(geom["e"], dtype=float)))


def _mean_inclination_rad(geom: Dict[str, Any]) -> float:
    """平均井斜角 β [rad]；缺 inc_deg 数据时按垂直井处理（β=0，S_buoy=0）。"""
    inc = geom.get("inc_deg")
    if inc is None:
        return 0.0
    return math.radians(float(np.mean(np.asarray(inc, dtype=float))))


def _mean_radius_m(geom: Dict[str, Any]) -> float:
    """环空平均半径 r_a = mean((D_hole + D_od)/4) [m]；缺井径数据时用 0.10 m fallback。"""
    hole = geom.get("hole_mm")
    od = geom.get("od_mm")
    if hole is None or od is None:
        return _FALLBACK_MEAN_RADIUS_M
    hole_arr = np.asarray(hole, dtype=float)
    od_arr = np.asarray(od, dtype=float)
    return float(np.mean((hole_arr + od_arr) / 4.0) / 1000.0)


def _half_gap_m(geom: Dict[str, Any]) -> float:
    """半间隙 d = mean(geom['b'])/2 [m]（geom['b'] 为局部全间隙），带下界保护。"""
    b_mean = float(np.mean(np.asarray(geom["b"], dtype=float)))
    return max(b_mean / 2.0, 1.0e-6)


def _pre_shutdown_velocity_m_s(schedule: "PumpingSchedule", geom: Dict[str, Any]) -> float:
    """估计停泵前截面平均速度 w₀ [m/s]。

    取第一个 SHUTDOWN 事件之前最后一个正排量步骤的排量 Q，
    除以环空截面积 A = π/4·(D_hole² − D_od²)。

    - 停泵前无任何正排量步骤 → 0.0（速度场已静止，t_s=0）；
    - 缺井径数据或截面积无效 → 典型顶替速度 0.5 m/s fallback。
    """
    rate_m3_min = 0.0
    for step in schedule.steps:
        if _event_tag_value(getattr(step, "event_tag", None)) == "SHUTDOWN":
            break
        rate = getattr(step, "rate_m3_min", 0.0) or 0.0
        if float(rate) > 0.0:
            rate_m3_min = float(rate)
    if rate_m3_min <= 0.0:
        return 0.0

    hole = geom.get("hole_mm")
    od = geom.get("od_mm")
    if hole is None or od is None:
        return _FALLBACK_VELOCITY_M_S
    d_hole = float(np.mean(np.asarray(hole, dtype=float))) / 1000.0
    d_od = float(np.mean(np.asarray(od, dtype=float))) / 1000.0
    area = math.pi / 4.0 * max(d_hole ** 2 - d_od ** 2, 0.0)
    if area <= 1.0e-12:
        return _FALLBACK_VELOCITY_M_S
    return (rate_m3_min / 60.0) / area


def compute_shutdown_decay(
    result: "AnnulusSimulationResult",
    fluids: Sequence["FluidSpec"],
    schedule: "PumpingSchedule",
) -> ShutdownDecayResult:
    """停泵有限时间衰减诊断（Moyers-González et al. 2007, 式 3.35/3.40）。

    物理：停泵后无外加流量（Q=0），若屈服应力足以平衡浮力散度驱动的
    有效应力，则稳态解 Ψs=0 且瞬态解在有限时间 t_s 内衰减到零（界面冻结）；
    否则停泵期间界面持续密度驱动运移（窜流风险）。

    判据（式 3.35 有量纲重建）：
        τ_Y,min/(1+e) ≥ S_buoy = π·r_a·Δρ·g·sinβ/2

    冻结时间（式 3.40 有量纲重建，判据满足时）：
        t_s = (1/α₁)·ln(1 + z₀·α₁/α₂)
        α₁ = μ_ref·(1+e)/(ρ_min·d²)，α₂ = 2(1+e)·Δτ/(ρ_min·π²·d)，
        z₀ = w₀·(1+e)/π，Δτ = τ_Y,min/(1+e) − S_buoy > 0

    触发条件：schedule 含 SHUTDOWN event_tag 时诊断才有意义；
    否则返回 condition_satisfied=False, freeze_time_s=inf 并注明不适用。

    数值安全：τ_Y,min=0（牛顿流体）时判不满足（无屈服则不停滞）；
    半间隙、密度、α₁ 均有下界保护；w₀=0 时 t_s=0（已静止）。

    Args:
        result: AnnulusSimulationResult（使用 geom 的 e/b/inc_deg/hole_mm/od_mm）
        fluids: 流体序列（FluidSpec 或鸭子类型兼容对象），读取
            yield_stress_pa / plastic_viscosity_pa_s / density_kg_m3
        schedule: PumpingSchedule（或鸭子类型兼容对象），检测 SHUTDOWN
            事件并估计停泵前排量

    Returns:
        ShutdownDecayResult
    """
    geom = result.geom
    tau_y_min = _min_yield_stress(fluids)

    # 触发检查：无停泵事件 → 诊断不适用
    if not _has_shutdown_event(schedule):
        return ShutdownDecayResult(
            condition_satisfied=False,
            freeze_time_s=math.inf,
            tau_y_min=tau_y_min,
            physical_interpretation=(
                "泵注程序中未检测到停泵（SHUTDOWN）事件，停泵有限时间衰减诊断不适用；"
                "若施工中存在停泵环节，请在 PumpingScheduleStep 上标注 "
                "event_tag=PumpingStageEvent.SHUTDOWN 后重新评估。"
            ),
        )

    # 牛顿流体守卫：无屈服应力支撑则速度仅指数衰减，不会有限时间归零
    if tau_y_min <= 0.0:
        return ShutdownDecayResult(
            condition_satisfied=False,
            freeze_time_s=math.inf,
            tau_y_min=tau_y_min,
            physical_interpretation=(
                "所有流体屈服应力为零（牛顿流体），停泵后无屈服应力支撑，"
                "速度仅指数衰减而不会有限时间归零，界面无法冻结（式 3.35 恒不满足）。"
            ),
        )

    e = _mean_eccentricity(geom)
    beta = _mean_inclination_rad(geom)
    r_a = _mean_radius_m(geom)
    rho_min, delta_rho = _density_extremes(fluids)

    # 式 3.35：屈服应力支撑 vs 浮力散度等效应力
    s_buoy = math.pi * r_a * delta_rho * _GRAVITY_M_S2 * math.sin(beta) / 2.0
    tau_support = tau_y_min / (1.0 + e)

    if tau_support < s_buoy:
        return ShutdownDecayResult(
            condition_satisfied=False,
            freeze_time_s=math.inf,
            tau_y_min=tau_y_min,
            physical_interpretation=(
                f"式 3.35 不满足：屈服应力支撑 τ_Y,min/(1+e)={tau_support:.2f} Pa "
                f"< 浮力散度等效应力 π·r_a·Δρ·g·sinβ/2={s_buoy:.2f} Pa；"
                "停泵期间界面将持续密度驱动运移，存在窜流风险；"
                "建议提高浆体最小屈服应力，或降低密度差/井斜影响。"
            ),
        )

    # 式 3.40：冻结时间（判据满足，Δτ > 0）
    d = _half_gap_m(geom)
    mu_ref = _reference_viscosity(fluids)
    rho_min_safe = max(rho_min, 1.0e-6)
    w0 = _pre_shutdown_velocity_m_s(schedule, geom)

    delta_tau = tau_support - s_buoy
    alpha1 = max(mu_ref * (1.0 + e) / (rho_min_safe * d * d), 1.0e-12)  # [1/s]
    alpha2 = 2.0 * (1.0 + e) * delta_tau / (rho_min_safe * math.pi ** 2 * d)  # [m/s²]
    z0 = w0 * (1.0 + e) / math.pi  # [m/s]

    if z0 <= 0.0:
        t_s = 0.0  # 停泵前已静止，立即冻结
    else:
        t_s = math.log1p(z0 * alpha1 / alpha2) / alpha1

    interpretation = (
        f"式 3.35 满足：τ_Y,min/(1+e)={tau_support:.2f} Pa ≥ "
        f"浮力散度等效应力 {s_buoy:.2f} Pa；"
        f"停泵后速度将在有限时间 t_s≈{t_s:.2f} s（式 3.40）内衰减至零，"
        "界面冻结，停泵期间不发生密度驱动窜流。"
    )
    if w0 <= 0.0:
        interpretation += "（停泵前排量为零，速度场已静止，t_s=0。）"

    return ShutdownDecayResult(
        condition_satisfied=True,
        freeze_time_s=float(t_s),
        tau_y_min=tau_y_min,
        physical_interpretation=interpretation,
    )
