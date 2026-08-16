"""
顶替效率与风险指标诊断层

本模块提供顶替效率计算、风险评估和诊断功能。

主要功能：
- 有效顶替效率计算（考虑泥饼清除后的实际顶替效果）
- 窜槽指数评估（宽窄边水泥前沿差异）
- 混浆指数评估（水泥-钻井液混合程度）
- 失稳指数评估（浮力导致的窄边滑塌风险）

设计原则：
- 区分水力学直接输出（位移效率）和质量解释输出
- 质量响应效率 ≠ 有效顶替效率，不能混为一谈
- 风险指标用于评估施工质量风险，不是效率的直接度量
"""

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
    classify_muskat_regime,
    compute_muskat_regime,
)
from cemdisp.diagnostics.regime_classifiers import (
    BuoyancyRegimeResult,
    ShutdownDecayResult,
    classify_buoyancy_regime,
    classify_buoyancy_regime_from_result,
    compute_shutdown_decay,
)
from cemdisp.diagnostics.tier0_diagnostics import (
    Tier0DiagnosticsResult,
    compute_all_tier0_diagnostics,
)

__all__ = [
    "DisplacementMetricsResult",
    "FlowClassificationResult",
    "compute_displacement_metrics",
    "compute_flow_classification",
    "BuoyancyRegimeResult",
    "ShutdownDecayResult",
    "classify_buoyancy_regime",
    "classify_buoyancy_regime_from_result",
    "compute_shutdown_decay",
    "MuskatRegimeResult",
    "classify_muskat_regime",
    "compute_muskat_regime",
    "Tier0DiagnosticsResult",
    "compute_all_tier0_diagnostics",
]
