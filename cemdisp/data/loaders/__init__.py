"""
逐井数据加载入口

本模块提供逐口井标准数据加载的入口函数和工具函数。

主要功能：
- 统一管理各井的数据加载器注册
- 定义井数据加载结果的类型别名

类型别名：
- Hu102LoaderResult: 呼102井数据加载器返回的四元组
    (WellSpec, Tuple[FluidSpec, ...], PumpingSchedule, ValidationData)

当前支持的井：
- 呼101: 通过 hu101_loader 模块加载
- 呼102: 通过 hu102_loader 模块加载
- 呼103: 通过 hu103_loader 模块加载
- 呼探1: 通过 hu1_loader 模块加载
- 呼探1-001（HT1-001）: 通过 ht1_001_loader 模块加载
- 呼探1-002（HT1-002）: 通过 hu2_loader 模块加载
"""

from cemdisp.data.fluid_spec import FluidSpec
from cemdisp.data.pumping_schedule import PumpingSchedule
from cemdisp.data.validation_data import ValidationData
from cemdisp.data.well_spec import WellSpec

from cemdisp.data.loaders.hu101_loader import (
    export_hu101_sync_card_markdown,
    load_hu101_tailpipe,
)
from cemdisp.data.loaders.hu102_loader import (
    export_hu102_sync_card_markdown,
    load_hu102_tailpipe,
)
from cemdisp.data.loaders.hu103_loader import (
    export_hu103_sync_card_markdown,
    load_hu103_tailpipe,
)
from cemdisp.data.loaders.hu1_loader import (
    export_hu1_sync_card_markdown,
    load_hu1_tailpipe,
)
from cemdisp.data.loaders.hu2_loader import (
    export_hu2_sync_card_markdown,
    load_hu2_tailpipe,
)
from cemdisp.data.loaders.ht1_001_loader import (
    export_ht1_001_sync_card_markdown,
    load_ht1_001_tailpipe,
)
from cemdisp.data.loaders.ht1_003_loader import (
    load_ht1_003_tailpipe,
)


Hu101LoaderResult = tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]
Hu102LoaderResult = tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]
Hu103LoaderResult = tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]
Hu1LoaderResult = tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]
Hu2LoaderResult = tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]
Ht1_001LoaderResult = tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]
Ht1_003LoaderResult = tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]

__all__ = [
    "Hu101LoaderResult",
    "Hu102LoaderResult",
    "Hu103LoaderResult",
    "Hu1LoaderResult",
    "Hu2LoaderResult",
    "Ht1_001LoaderResult",
    "Ht1_003LoaderResult",
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
