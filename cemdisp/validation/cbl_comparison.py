"""
CBL 实测数据对比验证接口

本模块提供模型预测结果与现场 CBL（水泥胶结测井）实测数据的对比验证功能。

设计原则：
- 只读对比：读取现场 CBL 报告，与模型水力效率进行定性/半定量对比
- 不反向校准：对比结果仅用于评估模型可信度，禁止用于修改求解器参数
- 输出偏差分析：计算预测值与实测值的绝对偏差、相对偏差和趋势一致性
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple


@dataclass(frozen=True)
class CBLValidationResult:
    """CBL 验证对比结果。

    Attributes:
        well_name: 井名
        cbl_interval_top_m: CBL 评价井段顶深（米）
        cbl_interval_bottom_m: CBL 评价井段底深（米）
        simulated_efficiency: 模型预测的 CBL 评价井段有效顶替效率
        measured_pass_rate: 现场 CBL 实测合格率（如 0.6665 表示 66.65%）
        absolute_delta: 绝对偏差（预测 - 实测）
        relative_delta_pct: 相对偏差（%）
        trend_consistent: 是否与趋势一致（定性判断）
        notes: 解释性备注
    """

    well_name: str
    cbl_interval_top_m: float
    cbl_interval_bottom_m: float
    simulated_efficiency: float
    measured_pass_rate: float | None
    absolute_delta: float | None
    relative_delta_pct: float | None
    trend_consistent: str | None
    notes: tuple[str, ...]


def validate_against_cbl(
    well_name: str,
    cbl_interval_top_m: float,
    cbl_interval_bottom_m: float,
    simulated_efficiency: float,
    measured_pass_rate: float | None,
    *,
    tolerance_absolute: float = 0.15,
    tolerance_relative_pct: float = 20.0,
) -> CBLValidationResult:
    """将模型预测的顶替效率与现场 CBL 实测合格率进行对比验证。

    注意：本函数仅做后验对比，不修改任何模型参数。
    CBL 合格率受多因素影响（泥饼、收缩、温度、气体等），
    模型预测的水力效率通常会高于 CBL 实测合格率，这是正常的物理差异。

    Args:
        well_name: 井名
        cbl_interval_top_m: CBL 评价井段顶深（米）
        cbl_interval_bottom_m: CBL 评价井段底深（米）
        simulated_efficiency: 模型预测的 CBL 评价井段有效顶替效率 [0, 1]
        measured_pass_rate: 现场 CBL 实测合格率 [0, 1]，若 None 表示无实测数据
        tolerance_absolute: 绝对偏差容忍阈值，默认 0.15
        tolerance_relative_pct: 相对偏差容忍阈值（%），默认 20.0%

    Returns:
        CBLValidationResult
    """
    if measured_pass_rate is None:
        notes = (
            "无现场 CBL 实测数据，无法进行定量对比。",
            "建议补充 CBL 测井报告或固井质量评价资料。",
        )
        trend = "无实测数据"
        abs_delta: Optional[float] = None
        rel_delta: Optional[float] = None
    else:
        abs_delta = float(simulated_efficiency - measured_pass_rate)
        rel_delta = (abs_delta / max(measured_pass_rate, 1e-9)) * 100.0 if measured_pass_rate > 0 else None

        # 趋势一致性判断
        if abs_delta is not None:
            if abs(abs_delta) <= tolerance_absolute:
                trend = "基本一致（绝对偏差在容忍范围内）"
            elif abs_delta > 0:
                trend = "水力效率高于 CBL 实测（符合预期，CBL 受多因素降低）"
            else:
                trend = "水力效率低于 CBL 实测（需检查模型参数或实测数据）"
        else:
            trend = "无法判断"

        notes_list = [
            f"水力效率（{simulated_efficiency:.4f}）与 CBL 合格率（{measured_pass_rate:.4f}）的对比：",
            f"绝对偏差 = {abs_delta:.4f}，相对偏差 = {rel_delta:.2f}%（若可计算）。",
        ]
        if abs_delta is not None and abs(abs_delta) > tolerance_absolute:
            notes_list.append(
                f"注意：偏差绝对值（{abs(abs_delta):.4f}）超过容忍阈值（{tolerance_absolute:.4f}），"
                "可能原因包括：泥饼残留、水泥收缩、温度效应、气体窜槽等非水力因素。"
            )
        else:
            notes_list.append("偏差在容忍范围内，模型预测与实测数据基本吻合。")

        notes = tuple(notes_list)

    return CBLValidationResult(
        well_name=well_name,
        cbl_interval_top_m=cbl_interval_top_m,
        cbl_interval_bottom_m=cbl_interval_bottom_m,
        simulated_efficiency=round(float(simulated_efficiency), 4),
        measured_pass_rate=measured_pass_rate,
        absolute_delta=abs_delta,
        relative_delta_pct=rel_delta,
        trend_consistent=trend,
        notes=notes,
    )


def summarize_validation(result: CBLValidationResult) -> str:
    """将验证结果格式化为中文报告摘要字符串。

    Args:
        result: CBLValidationResult 对象

    Returns:
        Markdown 格式的摘要字符串
    """
    lines = [
        f"## {result.well_name} CBL 实测对比验证",
        "",
        f"- CBL 评价井段：{result.cbl_interval_top_m:.2f} - {result.cbl_interval_bottom_m:.2f} m",
        f"- 模型预测水力效率：{result.simulated_efficiency:.4f}",
    ]
    if result.measured_pass_rate is not None:
        lines.append(f"- 现场 CBL 实测合格率：{result.measured_pass_rate:.4f}")
        lines.append(f"- 绝对偏差：{result.absolute_delta:.4f}")
        if result.relative_delta_pct is not None:
            lines.append(f"- 相对偏差：{result.relative_delta_pct:.2f}%")
        lines.append(f"- 趋势一致性：{result.trend_consistent}")
    else:
        lines.append("- 现场 CBL 实测数据：无")
    lines.append("")
    lines.append("### 备注")
    for note in result.notes:
        lines.append(f"- {note}")
    return "\n".join(lines)
