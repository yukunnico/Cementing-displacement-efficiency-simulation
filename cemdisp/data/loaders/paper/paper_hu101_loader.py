"""Paper-only wrapper for 呼101."""

from __future__ import annotations

from pathlib import Path

from cemdisp.data.loaders.hu101_loader import load_hu101_tailpipe
from cemdisp.data.loaders.paper.common import (
    PaperLoaderResult,
    PaperWellMetadata,
    validate_paper_loader_result,
    validation_with_cbl_pass_rate,
)


_METADATA = PaperWellMetadata(
    well_id="hu101",
    well_name_cn="呼101",
    paper_version="pead_v1",
    sample_class="A_main_validation",
    loader_version="paper_loader_v1",
    include_in_main_results=True,
    include_in_cbl_metrics=True,
    cbl_pass_rate=0.6277,
    geometry_source="field_measured + equivalent composite liner approximation",
    pumping_schedule_source="field_measured construction sequence with shoe-lag correction inherited from legacy loader",
    fluid_source="field_measured densities and rheology from extracted field report",
    caliper_source="legacy nominal/degenerate profile; adopted as paper v1 equivalent geometry",
    inclination_source="legacy nominal profile",
    standoff_source="model_assumption profile documented in loader notes",
    cbl_source="interpreted official CBL report 100312.PDF",
    notes=(
        "主验证候选井；CBL合格率采用0.6277。",
        "上部168.3mm+下部139.7mm复合尾管按面积守恒等效几何进入当前单外径求解器。",
        "193.7mm回接段资料不进入本 paper loader 主验证口径。",
    ),
)


def get_paper_metadata() -> PaperWellMetadata:
    """Return paper metadata for 呼101."""

    return _METADATA


def load_paper_hu101_tailpipe(*, reference_root: Path | None = None) -> PaperLoaderResult:
    """Load 呼101 using the paper-only provenance wrapper."""

    well_spec, fluids, schedule, validation_data = load_hu101_tailpipe(reference_root=reference_root)
    validation_data = validation_with_cbl_pass_rate(validation_data, _METADATA.cbl_pass_rate)
    return validate_paper_loader_result((well_spec, fluids, schedule, validation_data), _METADATA)
