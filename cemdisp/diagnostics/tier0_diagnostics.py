"""
Tier 0 诊断聚合层 (tier0_diagnostics)
=====================================

本模块聚合 Tier 0 路线图的 4 个诊断模块（见
docs/superpowers/specs/2026-07-17-d2dga-full-roadmap-tier0-design.md §2.8）：

1. ``flow_classification``（T0-1）：流动分类（Zhang & Frigaard 2023 式 3.1-3.6）；
2. ``muskat_regime``（T0-2）：Muskat 三 regime 稳定性（Bararpour & Frigaard 2025 式 3.8-3.13）；
3. ``displacement_metrics``（T0-4/5/7）：泥浆滞留 / 界面长度比 / η_N / 突破时间；
4. ``buoyancy_regime`` + ``shutdown_decay``（T0-3/6）：浮力数阈值分类与停泵有限时间衰减。

设计原则
--------
1. 纯后处理：只消费 ``AnnulusSimulationResult`` 与输入数据（fluids/well_spec/schedule），
   不 import models2d 求解逻辑（类型注解走 TYPE_CHECKING）。
2. 容错聚合：每个子诊断用 try/except 独立包裹，单个失败不拖垮整体；
   失败子项记 None，失败原因追加到 ``notes``。
3. 可序列化：``to_dict()`` 递归转 dict，NaN/Inf → None，tuple → list，
   可直接 ``json.dumps``（严格 JSON，无 NaN/Infinity 字面量）。

使用示例
--------
>>> from cemdisp.diagnostics.tier0_diagnostics import compute_all_tier0_diagnostics
>>> diag = compute_all_tier0_diagnostics(result, fluids=fluids, well_spec=well_spec, schedule=schedule)
>>> print(diag.flow_classification.flow_class)
>>> import json; json.dumps(diag.to_dict(), ensure_ascii=False)
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Dict, Optional, Sequence, Tuple

from cemdisp.diagnostics.displacement_metrics import (
    DisplacementMetricsResult,
    compute_displacement_metrics,
)
from cemdisp.diagnostics.flow_classification import (
    FlowClassificationResult,
    compute_flow_classification,
)
from cemdisp.diagnostics.muskat_regime import (
    MuskatRegimeResult,
    compute_muskat_regime,
)
from cemdisp.diagnostics.regime_classifiers import (
    BuoyancyRegimeResult,
    ShutdownDecayResult,
    classify_buoyancy_regime_from_result,
    compute_shutdown_decay,
)

if TYPE_CHECKING:  # 仅类型注解，避免运行时依赖 models2d 求解逻辑 / data 结构
    from cemdisp.data.fluid_spec import FluidSpec
    from cemdisp.data.pumping_schedule import PumpingSchedule
    from cemdisp.data.well_spec import WellSpec
    from cemdisp.models2d.annulus_d2dga import AnnulusSimulationResult


def _json_safe(value: Any) -> Any:
    """递归转为严格 JSON 可序列化对象：NaN/Inf → None，tuple → list。

    各子结果的 ``to_dict()`` 已处理大部分转换，但部分字段仍可能含
    NaN/Inf（如 t_br_s=inf、freeze_time_s=inf、w0_m_s=nan），此处统一兜底。
    """
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, bool):  # bool 是 int 子类，须先于数值判断
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, (int, str)) or value is None:
        return value
    # numpy 标量等兜底：尽量转 float，失败则转 str
    try:
        as_float = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return str(value)
    return as_float if math.isfinite(as_float) else None


@dataclass(frozen=True)
class Tier0DiagnosticsResult:
    """Tier 0 诊断聚合结果（T0-1 ~ T0-7 汇总）。

    Attributes:
        flow_classification: T0-1 流动分类结果；失败为 None
        muskat_regime: T0-2 Muskat 三 regime 结果；失败为 None
        displacement_metrics: T0-4/5/7 顶替指标结果；失败为 None
        buoyancy_regime: T0-3 浮力数阈值分类结果；失败为 None
        shutdown_decay: T0-6 停泵有限时间衰减结果；失败为 None
        notes: 聚合级备注（各子诊断失败原因、缺失输入说明等）
    """

    flow_classification: Optional[FlowClassificationResult]
    muskat_regime: Optional[MuskatRegimeResult]
    displacement_metrics: Optional[DisplacementMetricsResult]
    buoyancy_regime: Optional[BuoyancyRegimeResult]
    shutdown_decay: Optional[ShutdownDecayResult]
    notes: Tuple[str, ...] = ()

    def to_dict(self) -> Dict[str, Any]:
        """导出为严格 JSON 可序列化字典（NaN/Inf → None，子结果 None → None）。"""
        return _json_safe(
            {
                "flow_classification": (
                    self.flow_classification.to_dict() if self.flow_classification is not None else None
                ),
                "muskat_regime": (
                    self.muskat_regime.to_dict() if self.muskat_regime is not None else None
                ),
                "displacement_metrics": (
                    self.displacement_metrics.to_dict() if self.displacement_metrics is not None else None
                ),
                "buoyancy_regime": (
                    self.buoyancy_regime.to_dict() if self.buoyancy_regime is not None else None
                ),
                "shutdown_decay": (
                    self.shutdown_decay.to_dict() if self.shutdown_decay is not None else None
                ),
                "notes": list(self.notes),
            }
        )


def compute_all_tier0_diagnostics(
    result: "AnnulusSimulationResult",
    *,
    fluids: Optional[Sequence["FluidSpec"]] = None,
    well_spec: Optional["WellSpec"] = None,
    schedule: Optional["PumpingSchedule"] = None,
    w0_m_s: Optional[float] = None,
) -> Tier0DiagnosticsResult:
    """聚合计算全部 Tier 0 诊断（T0-1 ~ T0-7，主入口）。

    每个子诊断独立 try/except 包裹：单个失败不拖垮整体，失败子项记 None
    并将异常信息追加到结果 notes。

    Args:
        result: 环空二维求解结果（AnnulusSimulationResult）
        fluids: 流体规格序列（muskat_regime 与 shutdown_decay 所需；None 时跳过这两项）
        well_spec: 井规格（可选；muskat_regime 的偏心度回退用）
        schedule: 泵注程序（shutdown_decay 所需；None 时跳过该项）
        w0_m_s: 可选平均泵速 [m/s]，透传给 flow_classification 与
            displacement_metrics 的 ŵ₀ 归一化；None 时各模块自动估计

    Returns:
        Tier0DiagnosticsResult
    """
    notes: list[str] = []

    # ---- T0-1 流动分类（仅需 result） ----
    flow_classification: Optional[FlowClassificationResult]
    try:
        flow_classification = compute_flow_classification(result, w0_m_s=w0_m_s)
    except Exception as exc:  # noqa: BLE001 — 聚合层容错，失败记 notes
        flow_classification = None
        notes.append(f"flow_classification 计算失败：{type(exc).__name__}: {exc}")

    # ---- T0-4/5/7 顶替指标（仅需 result） ----
    displacement_metrics: Optional[DisplacementMetricsResult]
    try:
        displacement_metrics = compute_displacement_metrics(result, w0_m_s=w0_m_s)
    except Exception as exc:  # noqa: BLE001
        displacement_metrics = None
        notes.append(f"displacement_metrics 计算失败：{type(exc).__name__}: {exc}")

    # ---- T0-3 浮力数阈值分类（需 result.summary['buoyancy_number']） ----
    buoyancy_regime: Optional[BuoyancyRegimeResult]
    try:
        buoyancy_regime = classify_buoyancy_regime_from_result(result)
    except Exception as exc:  # noqa: BLE001
        buoyancy_regime = None
        notes.append(f"buoyancy_regime 计算失败：{type(exc).__name__}: {exc}")

    # ---- T0-2 Muskat 三 regime（需 fluids；well_spec 可选） ----
    muskat_regime: Optional[MuskatRegimeResult]
    if fluids is None:
        muskat_regime = None
        notes.append("muskat_regime 跳过：未提供 fluids（需含 role='mud' 与 role='lead'/'tail'）。")
    else:
        try:
            muskat_regime = compute_muskat_regime(
                result, fluids, well_spec, mean_velocity_m_s=w0_m_s
            )
        except Exception as exc:  # noqa: BLE001
            muskat_regime = None
            notes.append(f"muskat_regime 计算失败：{type(exc).__name__}: {exc}")

    # ---- T0-6 停泵有限时间衰减（需 fluids + schedule） ----
    shutdown_decay: Optional[ShutdownDecayResult]
    if fluids is None or schedule is None:
        shutdown_decay = None
        missing = [name for name, val in (("fluids", fluids), ("schedule", schedule)) if val is None]
        notes.append(f"shutdown_decay 跳过：未提供 {' 与 '.join(missing)}。")
    else:
        try:
            shutdown_decay = compute_shutdown_decay(result, fluids, schedule)
        except Exception as exc:  # noqa: BLE001
            shutdown_decay = None
            notes.append(f"shutdown_decay 计算失败：{type(exc).__name__}: {exc}")

    return Tier0DiagnosticsResult(
        flow_classification=flow_classification,
        muskat_regime=muskat_regime,
        displacement_metrics=displacement_metrics,
        buoyancy_regime=buoyancy_regime,
        shutdown_decay=shutdown_decay,
        notes=tuple(notes),
    )
