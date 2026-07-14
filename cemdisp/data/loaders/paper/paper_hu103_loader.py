"""Paper-only wrapper for 呼103."""

from __future__ import annotations

from pathlib import Path

from cemdisp.data.loaders.hu103_loader import load_hu103_tailpipe
from cemdisp.data.loaders.paper.common import (
    PaperLoaderResult,
    PaperWellMetadata,
    validate_paper_loader_result,
    validation_with_cbl_pass_rate,
)


_METADATA = PaperWellMetadata(
    well_id="hu103",
    well_name_cn="呼103",
    paper_version="pead_v1",
    sample_class="A_main_validation",
    loader_version="paper_loader_v1",
    include_in_main_results=True,
    include_in_cbl_metrics=True,
    cbl_pass_rate=0.1206,
    geometry_source="field measured 139.7mm lower liner interval with explicit composite-liner context",
    pumping_schedule_source="confirmed design/actual schedule represented by legacy loader for paper v1",
    fluid_source="field/model-ready fluid table represented by legacy loader",
    caliper_source="field caliper CSV used by legacy loader",
    inclination_source="field inclination CSV used by legacy loader",
    standoff_source="design centralization proxy adjusted by clearance; model_assumption",
    cbl_source="interpreted CBL for 139.7mm interval 7338-7712m",
    notes=(
        "复合尾管建模改进重点井；主CBL验证仅采用139.7mm段7338-7712m合格率0.1206。",
        "168.3mm段CBL 0.0004和整段综合0.0605单独保留在说明中，不混入主验证窗。",
        "当前求解器若按单外径运行，paper metadata必须保留等效/代理说明。",
    ),
)


def get_paper_metadata() -> PaperWellMetadata:
    """Return paper metadata for 呼103."""

    return _METADATA


def load_paper_hu103_tailpipe(
    *,
    caliper_csv_path: Path | None = None,
    reference_root: Path | None = None,
) -> PaperLoaderResult:
    """Load 呼103 using the paper-only provenance wrapper."""

    well_spec, fluids, schedule, validation_data = load_hu103_tailpipe(
        caliper_csv_path=caliper_csv_path,
        reference_root=reference_root,
    )
    validation_data = validation_with_cbl_pass_rate(validation_data, _METADATA.cbl_pass_rate)
    return validate_paper_loader_result((well_spec, fluids, schedule, validation_data), _METADATA)
