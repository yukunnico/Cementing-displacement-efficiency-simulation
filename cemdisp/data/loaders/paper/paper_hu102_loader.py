"""Paper-only wrapper for 呼102."""

from __future__ import annotations

from pathlib import Path

from cemdisp.data.loaders.hu102_loader import load_hu102_tailpipe
from cemdisp.data.loaders.paper.common import (
    PaperLoaderResult,
    PaperWellMetadata,
    validate_paper_loader_result,
    validation_with_cbl_pass_rate,
)


_METADATA = PaperWellMetadata(
    well_id="hu102",
    well_name_cn="呼102",
    paper_version="pead_v1",
    sample_class="A_main_validation",
    loader_version="paper_loader_v1",
    include_in_main_results=True,
    include_in_cbl_metrics=True,
    cbl_pass_rate=0.6665,
    geometry_source="field-rich source package; legacy reconstructed loader profile adopted for paper v1",
    pumping_schedule_source="legacy reconstructed schedule from field-rich records",
    fluid_source="field/extracted densities and rheology as represented by legacy loader",
    caliper_source="legacy loader profile pending full extracted CSV rebuild",
    inclination_source="legacy loader profile pending full extracted CSV rebuild",
    standoff_source="model_assumption because continuous field standoff is unavailable",
    cbl_source="interpreted official CBL report 100413.PDF",
    notes=(
        "主验证候选井；CBL合格率采用0.6665。",
        "设计文档指出旧loader需重构；paper v1明确标注当前仍采用legacy reconstructed口径。",
        "CBL评价窗与目标层段必须分别读取，不得混用。",
    ),
)


def get_paper_metadata() -> PaperWellMetadata:
    """Return paper metadata for 呼102."""

    return _METADATA


def load_paper_hu102_tailpipe(*, reference_root: Path | None = None) -> PaperLoaderResult:
    """Load 呼102 using the paper-only provenance wrapper."""

    well_spec, fluids, schedule, validation_data = load_hu102_tailpipe(reference_root=reference_root)
    validation_data = validation_with_cbl_pass_rate(validation_data, _METADATA.cbl_pass_rate)
    return validate_paper_loader_result((well_spec, fluids, schedule, validation_data), _METADATA)
