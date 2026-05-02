"""
逐井数据加载入口

本模块提供逐口井标准数据加载的入口函数和工具函数。

主要功能：
- 统一管理各井的数据加载器注册
- 提供已识别井名列表查询
- 定义井数据加载结果的类型别名

类型别名：
- Hu102LoaderResult: 呼102井数据加载器返回的四元组
    (WellSpec, Tuple[FluidSpec, ...], PumpingSchedule, ValidationData)

- Hu102InletProviderFactory: 环空入口边界提供器工厂函数类型
    接收(施工日程, 流体元组, 边界模式)，返回给定时间返回环空入口状态的函数

当前支持的井：
- 呼102: 通过 hu102_loader 模块加载

参考文档路径约定：
- REFERENCE_DOCS_ROOT: 项目根目录下"参考文档"文件夹
"""

from pathlib import Path
from typing import Callable, Tuple

from cemdisp.data.fluid_spec import FluidSpec
from cemdisp.data.pumping_schedule import PumpingSchedule
from cemdisp.data.validation_data import ValidationData
from cemdisp.data.well_spec import WellSpec
from cemdisp.models2d.boundary_bridge import AnnulusInletState

from cemdisp.data.loaders.hu102_loader import (
    build_hu102_annulus_inlet_provider,
    load_hu102_tailpipe,
)


REFERENCE_DOCS_ROOT = Path(__file__).resolve().parents[3] / "参考文档"


def available_well_names() -> Tuple[str, ...]:
    """返回当前参考文档目录中已识别的井名列表。"""

    if not REFERENCE_DOCS_ROOT.exists():
        return tuple()
    return tuple(sorted(entry.name for entry in REFERENCE_DOCS_ROOT.iterdir() if entry.is_dir()))


Hu102LoaderResult = Tuple[WellSpec, Tuple[FluidSpec, ...], PumpingSchedule, ValidationData]
Hu102InletProviderFactory = Callable[[PumpingSchedule, Tuple[FluidSpec, ...], str], Callable[[float], AnnulusInletState]]

__all__ = [
    "Hu102InletProviderFactory",
    "Hu102LoaderResult",
    "REFERENCE_DOCS_ROOT",
    "available_well_names",
    "build_hu102_annulus_inlet_provider",
    "load_hu102_tailpipe",
]
