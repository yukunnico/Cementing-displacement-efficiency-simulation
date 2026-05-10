"""注入流体来源口径兼容导出层。

真实登记已迁移到 :mod:`cemdisp.data.provenance`。保留本模块是为了兼容
既有 runner / 报告代码的流体来源 API，不在这里维护第二份登记表。
"""

from __future__ import annotations

from cemdisp.data.provenance import (
    FluidProvenance,
    WELL_FLUID_PROVENANCE,
    build_injected_fluid_provenance_summary,
    format_injected_fluid_provenance_markdown,
)


__all__ = [
    "FluidProvenance",
    "WELL_FLUID_PROVENANCE",
    "build_injected_fluid_provenance_summary",
    "format_injected_fluid_provenance_markdown",
]
