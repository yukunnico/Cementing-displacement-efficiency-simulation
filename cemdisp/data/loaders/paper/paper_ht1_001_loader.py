"""Paper-only wrapper for 呼探1-001 / HT1-001."""

from __future__ import annotations

from pathlib import Path

from cemdisp.data.loaders.ht1_001_loader import load_ht1_001_tailpipe
from cemdisp.data.loaders.paper.common import (
    PaperLoaderResult,
    PaperWellMetadata,
    validate_paper_loader_result,
    validation_with_cbl_pass_rate,
)


_METADATA = PaperWellMetadata(
    well_id="ht1_001",
    well_name_cn="呼探1-001井（HT1-001）",
    paper_version="pead_v1",
    sample_class="B_application_pending_validation",
    loader_version="paper_loader_v1",
    include_in_main_results=True,
    include_in_cbl_metrics=False,
    cbl_pass_rate=None,
    geometry_source="legacy/proxy geometry pending extracted-data cleanup",
    pumping_schedule_source="legacy/proxy schedule pending cleaned schedule generation",
    fluid_source="legacy/proxy fluid table pending source freeze",
    caliper_source="proxy/legacy profile",
    inclination_source="proxy/legacy profile",
    standoff_source="model_assumption",
    cbl_source="pending frozen CBL evaluation window",
    notes=(
        "设计文档要求先清洗schedule并冻结CBL窗口后再进入定量验证。",
        "paper v1允许输出应用性主结果，但不纳入CBL误差统计。",
        "几何和施工口径均需在后续资料清洗后升级。",
    ),
)


def get_paper_metadata() -> PaperWellMetadata:
    """Return paper metadata for HT1-001."""

    return _METADATA


def load_paper_ht1_001_tailpipe(*, reference_root: Path | None = None) -> PaperLoaderResult:
    """Load HT1-001 using the paper-only provenance wrapper."""

    well_spec, fluids, schedule, validation_data = load_ht1_001_tailpipe(reference_root=reference_root)
    validation_data = validation_with_cbl_pass_rate(validation_data, _METADATA.cbl_pass_rate)
    return validate_paper_loader_result((well_spec, fluids, schedule, validation_data), _METADATA)
