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
- 呼101: 通过 hu101_loader 模块加载
- 呼102: 通过 hu102_loader 模块加载
- 呼103: 通过 hu103_loader 模块加载
- 呼探1: 通过 hu1_loader 模块加载
- 呼探1-001（HT1-001）: 通过 ht1_001_loader 模块加载
- 呼探1-002（HT1-002）: 通过 hu2_loader 模块加载

参考文档路径约定：
- REFERENCE_DOCS_ROOT: 项目根目录下"参考文档"文件夹
"""

from pathlib import Path
from typing import Callable

from cemdisp.data.fluid_spec import FluidSpec
from cemdisp.data.pumping_schedule import PumpingSchedule
from cemdisp.data.validation_data import ValidationData
from cemdisp.data.well_spec import WellSpec
from cemdisp.models2d.boundary_bridge import AnnulusInletState

from cemdisp.data.loaders.hu101_loader import (
    build_hu101_annulus_inlet_provider,
    export_hu101_sync_card_markdown,
    load_hu101_tailpipe,
)
from cemdisp.data.loaders.hu102_loader import (
    build_hu102_annulus_inlet_provider,
    export_hu102_sync_card_markdown,
    load_hu102_tailpipe,
)
from cemdisp.data.loaders.hu103_loader import (
    build_hu103_annulus_inlet_provider,
    export_hu103_sync_card_markdown,
    load_hu103_tailpipe,
)
from cemdisp.data.loaders.hu1_loader import (
    build_hu1_annulus_inlet_provider,
    export_hu1_sync_card_markdown,
    load_hu1_tailpipe,
)
from cemdisp.data.loaders.hu2_loader import (
    build_hu2_annulus_inlet_provider,
    export_hu2_sync_card_markdown,
    load_hu2_tailpipe,
)
from cemdisp.data.loaders.ht1_001_loader import (
    build_ht1_001_annulus_inlet_provider,
    export_ht1_001_sync_card_markdown,
    load_ht1_001_tailpipe,
)
from cemdisp.data.loaders.ht1_003_loader import (
    load_ht1_003_tailpipe,
)


REFERENCE_DOCS_ROOT = Path(__file__).resolve().parents[3] / "参考文档"


def available_well_names() -> tuple[str, ...]:
    """返回当前参考文档目录中已识别的井名列表。"""

    if not REFERENCE_DOCS_ROOT.exists():
        return tuple()
    return tuple(sorted(entry.name for entry in REFERENCE_DOCS_ROOT.iterdir() if entry.is_dir()))


Hu101LoaderResult = tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]
Hu102LoaderResult = tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]
Hu103LoaderResult = tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]
Hu1LoaderResult = tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]
Hu2LoaderResult = tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]
Ht1_001LoaderResult = tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]
Ht1_003LoaderResult = tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]
Hu101InletProviderFactory = Callable[[PumpingSchedule, tuple[FluidSpec, ...], str], Callable[[float], AnnulusInletState]]
Hu102InletProviderFactory = Callable[[PumpingSchedule, tuple[FluidSpec, ...], str], Callable[[float], AnnulusInletState]]
Hu103InletProviderFactory = Callable[[PumpingSchedule, tuple[FluidSpec, ...], str], Callable[[float], AnnulusInletState]]
Hu1InletProviderFactory = Callable[[PumpingSchedule, tuple[FluidSpec, ...], str], Callable[[float], AnnulusInletState]]
Hu2InletProviderFactory = Callable[[PumpingSchedule, tuple[FluidSpec, ...], str], Callable[[float], AnnulusInletState]]
Ht1_001InletProviderFactory = Callable[[PumpingSchedule, tuple[FluidSpec, ...], str], Callable[[float], AnnulusInletState]]

__all__ = [
    "Hu101InletProviderFactory",
    "Hu101LoaderResult",
    "Hu102InletProviderFactory",
    "Hu102LoaderResult",
    "Hu103InletProviderFactory",
    "Hu103LoaderResult",
    "Hu1InletProviderFactory",
    "Hu1LoaderResult",
    "Hu2InletProviderFactory",
    "Hu2LoaderResult",
    "Ht1_001InletProviderFactory",
    "Ht1_001LoaderResult",
    "Ht1_003LoaderResult",
    "REFERENCE_DOCS_ROOT",
    "available_well_names",
    "build_hu101_annulus_inlet_provider",
    "build_hu102_annulus_inlet_provider",
    "build_hu103_annulus_inlet_provider",
    "build_hu1_annulus_inlet_provider",
    "build_hu2_annulus_inlet_provider",
    "build_ht1_001_annulus_inlet_provider",
    "export_hu101_sync_card_markdown",
    "export_hu102_sync_card_markdown",
    "export_hu103_sync_card_markdown",
    "export_hu1_sync_card_markdown",
    "export_hu2_sync_card_markdown",
    "export_ht1_001_sync_card_markdown",
    "load_hu101_tailpipe",
    "load_hu102_tailpipe",
    "load_hu103_tailpipe",
    "load_hu1_tailpipe",
    "load_hu2_tailpipe",
    "load_ht1_001_tailpipe",
    "load_ht1_003_tailpipe",
]