"""
CBL 质量风险代理解释层 (quality_proxy)
=========================================

本模块提供独立于水力核心求解器的质量解释功能，
将流体力学顶替结果转换为 CBL（水泥胶结测井）质量风险代理值。

设计原则
--------
1. **绝不反向校准求解器**：CBL 代理值只能用于后验验证与风险筛查，
   不能用于修改水力模型参数。
2. **明确区分概念**：
   - 水力学顶替效率 = 水泥浆体积占据率（物理过程）
   - CBL 质量代理值 = 基于工程经验的胶结质量预测区间（统计代理）
3. **输出带不确定区间**：给出点估计 + 置信区间，避免单值误导。

影响因素
--------
CBL 测井质量受多因素共同作用，包括但不限于：
- 顶替效率（水力模型直接输出）
- 泥饼残留与清除程度
- 水泥浆收缩率与微环隙
- 井筒温度分布与水化热
- 地层气体窜槽
- 套管居中度与环空偏心

本模块当前基于风险指标加权给出代理区间，后续可扩展为
多物理场耦合的贝叶斯代理模型。

使用示例
--------
>>> from cemdisp.diagnostics.quality_proxy import compute_cbl_quality_proxy
>>> result = compute_cbl_quality_proxy(
...     displacement_efficiency=0.85,
...     channeling_index=0.12,
...     mixing_index=0.05,
...     instability_index=0.08,
...     mud_cake_clearance=0.7,  # 泥饼清除率 70%
...     cement_shrinkage=0.03,   # 水泥收缩率 3%
... )
>>> print(result.point_estimate)   # 0.62
>>> print(result.lower_bound)      # 0.45
>>> print(result.upper_bound)      # 0.78
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np


@dataclass(frozen=True)
class CBLQualityProxyResult:
    """CBL 质量风险代理值结果。

    Attributes:
        point_estimate: 点估计值（最可能代理值）
        lower_bound: 下界（保守估计，如 P10）
        upper_bound: 上界（乐观估计，如 P90）
        confidence_level: 置信水平说明（如 "80% 置信区间"）
        hydraulic_efficiency: 输入的水力学顶替效率
        major_factors: 主要影响因素字典，key 为因素名，value 为影响系数
        notes: 解释性备注列表
    """

    point_estimate: float
    lower_bound: float
    upper_bound: float
    confidence_level: str
    hydraulic_efficiency: float
    major_factors: Dict[str, float]
    notes: tuple[str, ...]


def _compute_quality_factor(
    channeling_index: float,
    mixing_index: float,
    instability_index: float,
    *,
    channeling_weight: float = 0.55,
    mixing_weight: float = 0.35,
    instability_weight: float = 0.25,
    penalty_scale: float = 0.099,
) -> float:
    """基于风险指标计算综合质量折减因子。

    该因子反映水力过程不完善性对 CBL 质量的潜在负面影响，
    不等同于 CBL 实测合格率，仅为工程经验代理。

    Args:
        channeling_index: 窜槽指数 [0, 1]，越大表示宽窄边差异越严重
        mixing_index: 混浆指数 [0, 1]，越大表示界面混合越严重
        instability_index: 失稳指数 [0, 1]，越大表示窄边滑塌风险越高
        channeling_weight: 窜槽权重，默认 0.55
        mixing_weight: 混浆权重，默认 0.35
        instability_weight: 失稳权重，默认 0.25
        penalty_scale: 惩罚缩放系数，默认 0.099（基于 HU102 现场单点标定的历史值，
            仅作默认参考，用户应根据多井数据重新标定或采用不确定性分析）

    Returns:
        质量折减因子 [0, 1]，1.0 表示无折减，0.0 表示完全折减
    """
    raw_penalty = (
        channeling_weight * max(float(channeling_index), 0.0)
        + mixing_weight * max(float(mixing_index), 0.0)
        + instability_weight * max(float(instability_index), 0.0)
    )
    quality_factor = float(np.clip(1.0 - penalty_scale * raw_penalty, 0.0, 1.0))
    return quality_factor


def _apply_engineering_factors(
    base_proxy: float,
    *,
    mud_cake_clearance: Optional[float] = None,
    cement_shrinkage: Optional[float] = None,
    temperature_effect: Optional[float] = None,
    gas_migration_risk: Optional[float] = None,
    formation_quality: Optional[float] = None,
) -> tuple[float, Dict[str, float]]:
    """应用工程因子修正，将基础代理值调整为更全面的 CBL 质量预测。

    各因子均为 [0, 1] 范围内的比例系数：
    - mud_cake_clearance: 泥饼清除率（1.0 = 完全清除，0.0 = 无清除）
    - cement_shrinkage: 水泥收缩补偿率（1.0 = 无收缩问题，0.0 = 严重收缩）
    - temperature_effect: 温度效应系数（1.0 = 温度有利，0.0 = 温度严重不利）
    - gas_migration_risk: 气体窜槽风险控制（1.0 = 无风险，0.0 = 高风险）
    - formation_quality: 地层胶结条件（1.0 = 地层良好，0.0 = 地层极差）

    若因子未提供（None），则默认取中性值 0.85（存在中等不确定性）。

    Returns:
        (修正后代理值, 各因子影响系数字典)
    """
    factors: Dict[str, float] = {}
    adjustments: list[float] = []

    # 泥饼清除影响：泥饼残留会显著降低一界面胶结质量
    if mud_cake_clearance is not None:
        mud_factor = float(np.clip(mud_cake_clearance, 0.0, 1.0))
        factors["泥饼清除率"] = mud_factor
        adjustments.append(0.10 + 0.90 * mud_factor)  # 权重 10%-100%
    else:
        factors["泥饼清除率(默认)"] = 0.85
        adjustments.append(0.865)  # 中性值

    # 水泥收缩影响：收缩产生微环隙，降低 CBL 响应
    if cement_shrinkage is not None:
        shrink_factor = float(np.clip(1.0 - cement_shrinkage, 0.0, 1.0))
        factors["水泥收缩补偿"] = shrink_factor
        adjustments.append(0.05 + 0.95 * shrink_factor)
    else:
        factors["水泥收缩补偿(默认)"] = 0.85
        adjustments.append(0.8575)

    # 温度效应：温度分布不均影响水化程度和胶结质量
    if temperature_effect is not None:
        temp_factor = float(np.clip(temperature_effect, 0.0, 1.0))
        factors["温度效应"] = temp_factor
        adjustments.append(0.15 + 0.85 * temp_factor)
    else:
        factors["温度效应(默认)"] = 0.85
        adjustments.append(0.8725)

    # 气体窜槽风险：气体侵入水泥浆形成窜槽通道
    if gas_migration_risk is not None:
        gas_factor = float(np.clip(gas_migration_risk, 0.0, 1.0))
        factors["气体窜槽控制"] = gas_factor
        adjustments.append(0.20 + 0.80 * gas_factor)
    else:
        factors["气体窜槽控制(默认)"] = 0.85
        adjustments.append(0.88)

    # 地层胶结条件：地层本身疏松或裂缝会降低二界面胶结
    if formation_quality is not None:
        form_factor = float(np.clip(formation_quality, 0.0, 1.0))
        factors["地层胶结条件"] = form_factor
        adjustments.append(0.10 + 0.90 * form_factor)
    else:
        factors["地层胶结条件(默认)"] = 0.85
        adjustments.append(0.865)

    # 综合修正：采用乘积模型（各因子独立作用，联合效应为乘积）
    combined_adjustment = float(np.prod(adjustments) ** (1.0 / len(adjustments)))
    adjusted_proxy = float(np.clip(base_proxy * combined_adjustment, 0.0, 1.0))

    return adjusted_proxy, factors


def compute_cbl_quality_proxy(
    displacement_efficiency: float,
    channeling_index: float,
    mixing_index: float,
    instability_index: float,
    *,
    penalty_scale: float = 0.099,
    channeling_weight: float = 0.55,
    mixing_weight: float = 0.35,
    instability_weight: float = 0.25,
    mud_cake_clearance: Optional[float] = None,
    cement_shrinkage: Optional[float] = None,
    temperature_effect: Optional[float] = None,
    gas_migration_risk: Optional[float] = None,
    formation_quality: Optional[float] = None,
    uncertainty_level: float = 0.15,
) -> CBLQualityProxyResult:
    """计算 CBL 质量风险代理值及其置信区间。

    本函数将水力学顶替结果与工程经验因子结合，输出一个带不确定区间的
    CBL 质量代理预测。明确声明：该代理值不等同于 CBL 实测真值，
    仅供风险筛查和施工方案对比参考。

    Args:
        displacement_efficiency: 水力学有效顶替效率 [0, 1]
        channeling_index: 窜槽指数 [0, 1]
        mixing_index: 混浆指数 [0, 1]
        instability_index: 失稳指数 [0, 1]
        penalty_scale: 风险惩罚缩放系数，默认 0.099
        channeling_weight: 窜槽权重，默认 0.55
        mixing_weight: 混浆权重，默认 0.35
        instability_weight: 失稳权重，默认 0.25
        mud_cake_clearance: 泥饼清除率 [0, 1]，可选
        cement_shrinkage: 水泥收缩率 [0, 1]（注意：输入的是收缩率，
            函数内部会自动转换为补偿系数 = 1 - shrinkage），可选
        temperature_effect: 温度效应系数 [0, 1]，可选
        gas_migration_risk: 气体窜槽风险控制系数 [0, 1]，可选
        formation_quality: 地层胶结条件 [0, 1]，可选
        uncertainty_level: 不确定性水平，默认 0.15（15% 相对不确定性），
            用于计算上下界。若现场数据丰富，可降低至 0.10；
            若数据稀缺，应提高至 0.20-0.25。

    Returns:
        CBLQualityProxyResult，包含点估计、置信区间、主要影响因素和备注

    Raises:
        ValueError: 若输入参数超出合理范围
    """
    # 输入校验
    if not (0.0 <= float(displacement_efficiency) <= 1.0):
        raise ValueError(f"displacement_efficiency 必须在 [0, 1] 范围内，当前值：{displacement_efficiency}")
    if not (0.0 <= float(channeling_index) <= 1.0):
        raise ValueError(f"channeling_index 必须在 [0, 1] 范围内，当前值：{channeling_index}")
    if not (0.0 <= float(mixing_index) <= 1.0):
        raise ValueError(f"mixing_index 必须在 [0, 1] 范围内，当前值：{mixing_index}")
    if not (0.0 <= float(instability_index) <= 1.0):
        raise ValueError(f"instability_index 必须在 [0, 1] 范围内，当前值：{instability_index}")

    # Step 1: 基于水力风险指标计算基础质量因子
    quality_factor = _compute_quality_factor(
        channeling_index=channeling_index,
        mixing_index=mixing_index,
        instability_index=instability_index,
        channeling_weight=channeling_weight,
        mixing_weight=mixing_weight,
        instability_weight=instability_weight,
        penalty_scale=penalty_scale,
    )

    # Step 2: 水力代理值 = 顶替效率 × 质量因子
    hydraulic_proxy = float(np.clip(float(displacement_efficiency) * quality_factor, 0.0, 1.0))

    # Step 3: 应用工程因子修正
    adjusted_proxy, factors = _apply_engineering_factors(
        hydraulic_proxy,
        mud_cake_clearance=mud_cake_clearance,
        cement_shrinkage=cement_shrinkage,
        temperature_effect=temperature_effect,
        gas_migration_risk=gas_migration_risk,
        formation_quality=formation_quality,
    )

    # Step 4: 计算不确定性区间
    # 使用对数正态近似：在代理值两侧按 uncertainty_level 扩展
    # 下界不能小于 0，上界不能大于 1
    delta = float(uncertainty_level) * adjusted_proxy
    lower = float(np.clip(adjusted_proxy - delta, 0.0, 1.0))
    upper = float(np.clip(adjusted_proxy + delta, 0.0, 1.0))

    # 构建备注信息
    notes_list = [
        "CBL 质量代理值基于水力学顶替结果与工程经验因子联合预测，不等同于 CBL 实测真值。",
        f"水力学有效顶替效率 = {float(displacement_efficiency):.4f}，质量折减因子 = {quality_factor:.4f}。",
        f"工程因子综合修正系数 = {float(np.prod(list(factors.values())) ** (1.0 / len(factors))):.4f}。",
    ]
    if any(v is None for v in [mud_cake_clearance, cement_shrinkage, temperature_effect, gas_migration_risk, formation_quality]):
        notes_list.append("部分工程因子未提供，已使用默认值 0.85（存在中等不确定性），建议补充现场实测数据以缩小区间。")
    notes_list.append(f"当前不确定性水平设置为 {float(uncertainty_level)*100:.0f}%，置信区间宽度可通过调整 uncertainty_level 参数控制。")
    notes_list.append("本代理值仅用于风险筛查和方案对比，禁止用于反向校准水力求解器参数。")

    return CBLQualityProxyResult(
        point_estimate=round(adjusted_proxy, 4),
        lower_bound=round(lower, 4),
        upper_bound=round(upper, 4),
        confidence_level=f"{int((1.0 - uncertainty_level) * 100)}% 近似置信区间",
        hydraulic_efficiency=round(float(displacement_efficiency), 4),
        major_factors=factors,
        notes=tuple(notes_list),
    )


def summarize_proxy_for_report(result: CBLQualityProxyResult) -> str:
    """将 CBL 代理结果格式化为中文报告摘要字符串。

    Args:
        result: CBLQualityProxyResult 对象

    Returns:
        Markdown 格式的摘要字符串
    """
    lines = [
        "## CBL 质量风险代理预测",
        "",
        f"- 水力学有效顶替效率：{result.hydraulic_efficiency:.4f}",
        f"- CBL 质量代理点估计：{result.point_estimate:.4f}",
        f"- {result.confidence_level}：[{result.lower_bound:.4f}, {result.upper_bound:.4f}]",
        "",
        "### 主要影响因素",
    ]
    for factor_name, factor_value in result.major_factors.items():
        lines.append(f"- {factor_name}：{factor_value:.4f}")
    lines.append("")
    lines.append("### 备注")
    for note in result.notes:
        lines.append(f"- {note}")
    return "\n".join(lines)


# 便捷函数：从模拟结果 DataFrame 的最后一行快速计算代理值
def compute_from_simulation_metrics(
    final_metrics_row: dict,
    *,
    displacement_efficiency_key: str = "cbl_eval_interval_efficiency",
    channeling_key: str = "channeling_index",
    mixing_key: str = "mixing_index",
    instability_key: str = "instability_index",
    **kwargs,
) -> CBLQualityProxyResult:
    """从模拟结果指标字典（通常为 DataFrame 最后一行）快速计算 CBL 代理值。

    Args:
        final_metrics_row: 包含各指标的字典，如 {"cbl_eval_interval_efficiency": 0.85, ...}
        displacement_efficiency_key: 顶替效率字段名
        channeling_key: 窜槽指数字段名
        mixing_key: 混浆指数字段名
        instability_key: 失稳指数字段名
        **kwargs: 传递给 compute_cbl_quality_proxy 的其他参数

    Returns:
        CBLQualityProxyResult
    """
    return compute_cbl_quality_proxy(
        displacement_efficiency=float(final_metrics_row[displacement_efficiency_key]),
        channeling_index=float(final_metrics_row[channeling_key]),
        mixing_index=float(final_metrics_row[mixing_key]),
        instability_index=float(final_metrics_row[instability_key]),
        **kwargs,
    )
