"""
现场校验资料路径集合

本模块定义了与单口井相关的现场校验资料的文件路径结构，用于后续模型验证和报告生成。

主要字段：
- cbl_profile_path: CBL（水泥胶结测井）剖面原始数据文件路径
- cbl_summary_path: CBL评价汇总报告文件路径
- job_report_path: 固井施工总结报告文件路径
- pump_pressure_series_path: 泵压时序数据文件路径
- returns_report_path: 返排记录报告文件路径
- cement_top_path: 水泥面位置确认文件路径
- notes: 备注信息元组

注意：这些路径指向现场实测资料，用于与模型预测结果进行对比验证。
字段可为None，表示该资料尚未提供或不适用。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Tuple


@dataclass(frozen=True)
class ValidationData:
    """与单井相关的现场校验资料路径。"""

    cbl_profile_path: Optional[Path] = None
    cbl_summary_path: Optional[Path] = None
    # 实测 CBL 合格率（如 0.6665 表示 66.65%），若已知可直接填入，用于后验验证
    cbl_pass_rate: Optional[float] = None
    job_report_path: Optional[Path] = None
    pump_pressure_series_path: Optional[Path] = None
    returns_report_path: Optional[Path] = None
    cement_top_path: Optional[Path] = None
    notes: Tuple[str, ...] = field(default_factory=tuple)
