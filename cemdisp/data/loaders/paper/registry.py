"""Registry for paper-only well loaders."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable

from cemdisp.data.loaders.paper.common import PaperLoaderResult, PaperWellMetadata
from cemdisp.data.loaders.paper.paper_ht1_001_loader import (
    get_paper_metadata as get_ht1_001_metadata,
    load_paper_ht1_001_tailpipe,
)
from cemdisp.data.loaders.paper.paper_hu101_loader import (
    get_paper_metadata as get_hu101_metadata,
    load_paper_hu101_tailpipe,
)
from cemdisp.data.loaders.paper.paper_hu102_loader import (
    get_paper_metadata as get_hu102_metadata,
    load_paper_hu102_tailpipe,
)
from cemdisp.data.loaders.paper.paper_hu103_loader import (
    get_paper_metadata as get_hu103_metadata,
    load_paper_hu103_tailpipe,
)


@dataclass(frozen=True)
class PaperWellRecord:
    """Registered paper well loader entry."""

    well_id: str
    loader: Callable[[], PaperLoaderResult]
    metadata_factory: Callable[[], PaperWellMetadata]

    @property
    def metadata(self) -> PaperWellMetadata:
        return self.metadata_factory()


PAPER_WELL_REGISTRY: dict[str, PaperWellRecord] = {
    "hu101": PaperWellRecord("hu101", load_paper_hu101_tailpipe, get_hu101_metadata),
    "hu102": PaperWellRecord("hu102", load_paper_hu102_tailpipe, get_hu102_metadata),
    "hu103": PaperWellRecord("hu103", load_paper_hu103_tailpipe, get_hu103_metadata),
    "ht1_001": PaperWellRecord("ht1_001", load_paper_ht1_001_tailpipe, get_ht1_001_metadata),
}


def get_paper_well_record(well_id: str) -> PaperWellRecord:
    """Return the paper registry record for a well id."""

    try:
        return PAPER_WELL_REGISTRY[well_id]
    except KeyError as exc:
        available = ", ".join(PAPER_WELL_REGISTRY)
        raise KeyError(f"Unknown paper well {well_id!r}; available: {available}") from exc


def iter_paper_wells(*, include_pending: bool = False) -> Iterable[PaperWellRecord]:
    """Iterate paper wells in stable paper-table order."""

    for record in PAPER_WELL_REGISTRY.values():
        metadata = record.metadata
        if include_pending or metadata.include_in_cbl_metrics:
            yield record


def load_paper_well(well_id: str) -> PaperLoaderResult:
    """Load one registered paper well."""

    return get_paper_well_record(well_id).loader()


def get_paper_metadata(well_id: str) -> PaperWellMetadata:
    """Return metadata for one registered paper well."""

    return get_paper_well_record(well_id).metadata
