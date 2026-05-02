"""
环空二维顶替核心与边界桥接层

本模块实现偏心环空二维顶替效率计算的核心求解器，以及套管内1D到环空2D的边界桥接。

主要组件：
- AnnulusD2DGASolver: 环空二维D2DGA（偏心环空双液双区）求解器
    * run(): 执行环空二维顶替模拟
    * 输入：井筒规格、流体列表、环空入口状态提供器
    * 输出：AnnulusSimulationResult（水泥场、效率指标、风险指数等）

- AnnulusSimulationResult: 环空二维求解结果
    * geom: 几何参数（网格坐标、偏心度、井径等）
    * cement_field: 水泥浓度场（ny×nz数组）
    * wall_field: 壁面泥饼清除场
    * metrics: 时间序列指标（效率、前沿位置、风险指数）
    * depth_profiles: 深度方向平均剖面

- AnnulusInletState: 环空入口边界状态
    * 包含时间和相分数，供求解器使用

- pipe_exit_to_annulus_inlet(): 鞋口出流→环空入口的直接映射函数
- build_coupled_annulus_inlet_provider(): 1D-2D耦合的边界提供器工厂函数

设计原则：
- 环空求解器可独立运行
- 1D→2D边界桥接接口清晰
- 支持后续扩展套管内对流-弥散，而不影响环空核心
"""

from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver, AnnulusSimulationResult
from cemdisp.models2d.boundary_bridge import (
    AnnulusInletState,
    build_coupled_annulus_inlet_provider,
    pipe_exit_to_annulus_inlet,
)

__all__ = [
    "AnnulusD2DGASolver",
    "AnnulusInletState",
    "AnnulusSimulationResult",
    "build_coupled_annulus_inlet_provider",
    "pipe_exit_to_annulus_inlet",
]
