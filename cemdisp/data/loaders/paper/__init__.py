"""Paper-only loader public API."""

from cemdisp.data.loaders.paper.common import (
    PaperLoaderResult,
    PaperWellMetadata,
    get_window,
    metadata_to_dict,
    validate_paper_loader_result,
    validation_with_cbl_pass_rate,
)
from cemdisp.data.loaders.paper.paper_ht1_001_loader import load_paper_ht1_001_tailpipe
from cemdisp.data.loaders.paper.paper_hu101_loader import load_paper_hu101_tailpipe
from cemdisp.data.loaders.paper.paper_hu102_loader import load_paper_hu102_tailpipe
from cemdisp.data.loaders.paper.paper_hu103_loader import load_paper_hu103_tailpipe
from cemdisp.data.loaders.paper.registry import (
    PAPER_WELL_REGISTRY,
    PaperWellRecord,
    get_paper_metadata,
    get_paper_well_record,
    iter_paper_wells,
    load_paper_well,
)

__all__ = [
    "PAPER_WELL_REGISTRY",
    "PaperLoaderResult",
    "PaperWellMetadata",
    "PaperWellRecord",
    "get_paper_metadata",
    "get_paper_well_record",
    "get_window",
    "iter_paper_wells",
    "load_paper_ht1_001_tailpipe",
    "load_paper_hu101_tailpipe",
    "load_paper_hu102_tailpipe",
    "load_paper_hu103_tailpipe",
    "load_paper_well",
    "metadata_to_dict",
    "validate_paper_loader_result",
    "validation_with_cbl_pass_rate",
]
