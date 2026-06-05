"""
cemdisp - 多井通用固井顶替效率模型包

本包旨在构建一个多井通用的固井顶替效率模型。
核心功能包括：
- 标准输入数据结构（井筒规格、流体物性、施工程序）
- 套管内一维输运层（流体前缘追踪、鞋口出流状态）
- 环空二维顶替核心层（偏心环空顶替模拟）
- 诊断与质量解释层（效率计算、风险评估）
- 报告与可视化输出（中文图表）

主要子模块：
- cemdisp.data: 标准输入数据结构
- cemdisp.transport1d: 套管内一维输运层
- cemdisp.models2d: 环空二维顶替核心
- cemdisp.diagnostics: 诊断与指标计算
- cemdisp.reporting: 图表与报告输出
- cemdisp.runners: 各井段模型运行器
- cemdisp.validation: 数值验证
"""

from cemdisp.data import (
    DepthValuePoint,
    EvaluationWindow,
    FluidRole,
    FluidSpec,
    PumpingSchedule,
    PumpingScheduleStep,
    RheologyModel,
    ValidationData,
    WellSpec,
)
from cemdisp.models2d import (
    AnnulusD2DGASolver,
    AnnulusInletState,
    AnnulusSimulationResult,
    pipe_exit_to_annulus_inlet,
)
from cemdisp.transport1d import CasingFlowResult, InterfaceFront, PipeExitState

__all__ = [
    "AnnulusD2DGASolver",
    "AnnulusInletState",
    "AnnulusSimulationResult",
    "CasingFlowResult",
    "DepthValuePoint",
    "EvaluationWindow",
    "FluidRole",
    "FluidSpec",
    "InterfaceFront",
    "PipeExitState",
    "PumpingSchedule",
    "PumpingScheduleStep",
    "RheologyModel",
    "ValidationData",
    "WellSpec",
    "pipe_exit_to_annulus_inlet",
]
