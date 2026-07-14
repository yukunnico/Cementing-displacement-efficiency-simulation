"""Build CBL validation tables for paper outputs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Sequence

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from cemdisp.data.loaders.paper import get_paper_metadata, iter_paper_wells  # noqa: E402
from cemdisp.data.loaders.paper.common import write_csv, write_markdown_table  # noqa: E402
from scripts.paper_data.run_paper_main_results import RESULTS_ROOT  # noqa: E402

CBL_COLUMNS = [
    "well_id",
    "well_name_cn",
    "predicted_cbl_efficiency",
    "measured_cbl_pass_rate",
    "absolute_error",
    "relative_error_percent",
    "included_in_metrics",
    "notes",
]

TRUTH_COLUMNS = [
    "well_id",
    "well_name_cn",
    "measured_cbl_pass_rate",
    "included_in_metrics",
    "cbl_source",
    "notes",
]


def _round_or_blank(value: float | None, digits: int = 6) -> float | str:
    if value is None:
        return ""
    return round(float(value), digits)


def compute_cbl_error_metrics(
    well_id: str,
    *,
    well_name_cn: str,
    predicted: float | None,
    measured: float | None,
    included: bool,
    notes: str,
) -> dict[str, object]:
    """Compute one CBL error metrics row."""

    absolute_error: float | None = None
    relative_error_percent: float | None = None
    if predicted is not None and measured is not None:
        absolute_error = abs(float(predicted) - float(measured))
        if measured != 0.0:
            relative_error_percent = absolute_error / abs(float(measured)) * 100.0
    return {
        "well_id": well_id,
        "well_name_cn": well_name_cn,
        "predicted_cbl_efficiency": _round_or_blank(predicted),
        "measured_cbl_pass_rate": _round_or_blank(measured),
        "absolute_error": _round_or_blank(absolute_error),
        "relative_error_percent": _round_or_blank(relative_error_percent),
        "included_in_metrics": bool(included),
        "notes": notes,
    }


def _load_summary(root: Path, well_id: str) -> dict[str, Any]:
    summary_path = root / "01_main_results" / well_id / "summary.json"
    if not summary_path.exists():
        return {}
    return json.loads(summary_path.read_text(encoding="utf-8"))


def build_cbl_validation_tables(root: Path = RESULTS_ROOT) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Build CBL truth and error tables from summaries and registry metadata."""

    out_dir = root / "02_cbl_validation"
    out_dir.mkdir(parents=True, exist_ok=True)
    truth_rows: list[dict[str, object]] = []
    metric_rows: list[dict[str, object]] = []
    for record in iter_paper_wells(include_pending=True):
        metadata = record.metadata
        summary = _load_summary(root, record.well_id)
        predicted = summary.get("cbl_eval_interval_efficiency")
        predicted_float = float(predicted) if predicted not in (None, "") else None
        note = "; ".join(metadata.notes)
        truth_rows.append(
            {
                "well_id": metadata.well_id,
                "well_name_cn": metadata.well_name_cn,
                "measured_cbl_pass_rate": _round_or_blank(metadata.cbl_pass_rate),
                "included_in_metrics": bool(metadata.include_in_cbl_metrics),
                "cbl_source": metadata.cbl_source,
                "notes": note,
            }
        )
        metric_rows.append(
            compute_cbl_error_metrics(
                metadata.well_id,
                well_name_cn=metadata.well_name_cn,
                predicted=predicted_float,
                measured=metadata.cbl_pass_rate,
                included=metadata.include_in_cbl_metrics,
                notes=note,
            )
        )
    write_csv(out_dir / "cbl_truth_table.csv", truth_rows, TRUTH_COLUMNS)
    write_csv(out_dir / "cbl_error_metrics.csv", metric_rows, CBL_COLUMNS)
    write_markdown_table(out_dir / "cbl_comparison.md", metric_rows, CBL_COLUMNS)
    return truth_rows, metric_rows


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Build paper CBL validation tables")
    parser.add_argument("--output-root", type=Path, default=RESULTS_ROOT)
    args = parser.parse_args(argv)
    _, rows = build_cbl_validation_tables(args.output_root)
    print(f"Wrote CBL validation rows: {len(rows)}")


if __name__ == "__main__":
    main()
