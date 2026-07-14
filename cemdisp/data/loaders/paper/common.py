"""Common helpers for paper-only loader wrappers and outputs."""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any, Iterable, Sequence

from cemdisp.data.fluid_spec import FluidSpec, FluidRole
from cemdisp.data.pumping_schedule import PumpingSchedule
from cemdisp.data.validation_data import ValidationData
from cemdisp.data.well_spec import EvaluationWindow, WellSpec

PaperLoaderResult = tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]


@dataclass(frozen=True)
class PaperWellMetadata:
    """Provenance and paper-sample metadata attached to one paper loader."""

    well_id: str
    well_name_cn: str
    paper_version: str
    sample_class: str
    loader_version: str
    include_in_main_results: bool
    include_in_cbl_metrics: bool
    cbl_pass_rate: float | None
    geometry_source: str
    pumping_schedule_source: str
    fluid_source: str
    caliper_source: str
    inclination_source: str
    standoff_source: str
    cbl_source: str
    notes: tuple[str, ...] = ()


def metadata_to_dict(metadata: PaperWellMetadata) -> dict[str, object]:
    """Convert metadata to a JSON-serializable dict."""

    data = asdict(metadata)
    data["notes"] = list(metadata.notes)
    return data


def validation_with_cbl_pass_rate(
    validation_data: ValidationData,
    cbl_pass_rate: float | None,
) -> ValidationData:
    """Return validation data with the paper CBL pass rate attached."""

    if cbl_pass_rate is None:
        return validation_data
    return replace(validation_data, cbl_pass_rate=cbl_pass_rate)


def get_window(well_spec: WellSpec, window_type: str) -> EvaluationWindow | None:
    """Return the first evaluation window with the requested type."""

    for window in well_spec.evaluation_windows:
        if window.window_type == window_type:
            return window
    return None


def validate_paper_loader_result(
    result: PaperLoaderResult,
    metadata: PaperWellMetadata,
) -> PaperLoaderResult:
    """Validate a paper loader result and return it unchanged.

    The checks intentionally mirror the design document's loader quality gate while
    staying lightweight enough for unit tests.
    """

    well_spec, fluids, schedule, validation_data = result
    if not isinstance(well_spec, WellSpec):
        raise TypeError(f"{metadata.well_id}: loader must return WellSpec as item 0")
    if not all(isinstance(fluid, FluidSpec) for fluid in fluids):
        raise TypeError(f"{metadata.well_id}: loader must return tuple[FluidSpec, ...] as item 1")
    if not isinstance(schedule, PumpingSchedule):
        raise TypeError(f"{metadata.well_id}: loader must return PumpingSchedule as item 2")
    if not isinstance(validation_data, ValidationData):
        raise TypeError(f"{metadata.well_id}: loader must return ValidationData as item 3")
    if well_spec.top_md_m >= well_spec.bottom_md_m:
        raise ValueError(f"{metadata.well_id}: top depth must be shallower than bottom depth")
    if schedule.total_volume_m3 <= 0.0:
        raise ValueError(f"{metadata.well_id}: pumping schedule total volume must be positive")
    fluid_names = {fluid.name for fluid in fluids}
    for step in schedule.steps:
        if step.fluid_name not in fluid_names:
            raise ValueError(f"{metadata.well_id}: schedule fluid {step.fluid_name!r} has no FluidSpec")
        if step.volume_m3 < 0.0:
            raise ValueError(f"{metadata.well_id}: schedule volume must be non-negative")
    roles = {fluid.role for fluid in fluids}
    required_roles = {FluidRole.MUD, FluidRole.SPACER, FluidRole.TAIL, FluidRole.DISPLACEMENT}
    if not required_roles.issubset(roles):
        missing = ", ".join(sorted(role.value for role in required_roles.difference(roles)))
        raise ValueError(f"{metadata.well_id}: missing required fluid roles: {missing}")
    cbl_window = get_window(well_spec, "cbl")
    if cbl_window is not None:
        if cbl_window.top_md_m < well_spec.top_md_m or cbl_window.bottom_md_m > well_spec.bottom_md_m:
            raise ValueError(f"{metadata.well_id}: CBL window exceeds model interval")
    if not metadata.notes:
        raise ValueError(f"{metadata.well_id}: metadata notes must document assumptions/provenance")
    return result


def write_markdown_table(path: Path, rows: Sequence[dict[str, Any]], columns: Sequence[str]) -> Path:
    """Write a simple GitHub-flavored Markdown table."""

    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    lines.append("| " + " | ".join(columns) + " |")
    lines.append("| " + " | ".join("---" for _ in columns) + " |")
    for row in rows:
        values = [str(row.get(column, "")) for column in columns]
        lines.append("| " + " | ".join(value.replace("\n", " ") for value in values) + " |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_csv(path: Path, rows: Iterable[dict[str, Any]], columns: Sequence[str]) -> Path:
    """Write rows to a UTF-8 BOM CSV file with a stable column order."""

    path.parent.mkdir(parents=True, exist_ok=True)
    materialized = list(rows)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(materialized)
    return path
