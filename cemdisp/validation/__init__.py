"""
数值与现场对比验证层

本模块提供模型验证和对比分析功能。

主要功能：
1. 数值验证
   - 质量守恒验证（输入体积 vs 输出体积）
   - 时间推进稳定性验证
   - 能量/动量守恒检查

2. 与现场数据对比
   - CBL合格率对比（仅用于后验验证，不反向校准）
   - 水泥面位置对比
   - 泵压曲线对比

3. 敏感性分析
   - 参数敏感性评估
   - 不确定性量化

4. 验证报告生成
   - 可视化对比图
   - 误差分析摘要
   - 改进建议

使用流程：
1. 运行模型获得预测结果
2. 与现场CBL/施工数据对比（仅验证，不校准）
3. 进行误差分析和敏感性评估
4. 生成验证报告

注意：CBL 合格率 ≠ 流体力学顶替效率，两者物理定义不同。
验证结果只反映模型预测与现场实测的差异，不用于修改求解器参数。
"""

from cemdisp.validation.cbl_comparison import (
    CBLValidationResult,
    summarize_validation,
    validate_against_cbl,
)
from cemdisp.validation.mass_balance import CementMassBalanceDiagnostics, validate_cement_mass_balance

__all__ = [
    "CementMassBalanceDiagnostics",
    "validate_cement_mass_balance",
    "CBLValidationResult",
    "validate_against_cbl",
    "summarize_validation",
]
