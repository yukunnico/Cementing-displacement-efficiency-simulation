"""
顶替效率与风险指标诊断层

本模块提供顶替效率计算、风险评估和诊断功能。

主要功能：
- 有效顶替效率计算（考虑泥饼清除后的实际顶替效果）
- 窜槽指数评估（宽窄边水泥前沿差异）
- 混浆指数评估（水泥-钻井液混合程度）
- 失稳指数评估（浮力导致的窄边滑塌风险）
- CBL 质量风险代理值（基于风险指标与工程因子的独立质量预测）

设计原则：
- 区分水力学直接输出（位移效率）和质量解释输出（CBL代理值）
- 质量响应效率 ≠ 有效顶替效率，不能混为一谈
- 风险指标用于评估施工质量风险，不是效率的直接度量
- CBL 代理值仅用于验证与对比，禁止反向校准求解器
"""

from cemdisp.diagnostics.quality_proxy import (
    CBLQualityProxyResult,
    compute_cbl_quality_proxy,
    compute_from_simulation_metrics,
    summarize_proxy_for_report,
)

__all__ = [
    "CBLQualityProxyResult",
    "compute_cbl_quality_proxy",
    "compute_from_simulation_metrics",
    "summarize_proxy_for_report",
]
