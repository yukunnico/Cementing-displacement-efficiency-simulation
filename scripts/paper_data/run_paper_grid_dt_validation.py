"""Run paper grid and time-step validation cases."""

from __future__ import annotations

import argparse
import sys
from dataclasses import replace
from pathlib import Path
from typing import Sequence

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from cemdisp.data.loaders.paper.common import write_csv, write_markdown_table  # noqa: E402
from scripts.paper_data.run_paper_main_results import (  # noqa: E402
    RESULTS_ROOT,
    PaperModelConfig,
    default_paper_model_config,
    run_single_paper_well,
)

VALIDATION_COLUMNS = [
    "well_id",
    "nz",
    "ny",
    "dt",
    "final_effective_efficiency",
    "cbl_eval_interval_efficiency",
    "channeling_index",
    "mixing_index",
    "instability_index",
    "runtime_s",
    "reference_case",
    "relative_error_to_reference_percent",
]


def _relative_error(value: object, reference: object) -> float | str:
    try:
        value_f = float(value)
        reference_f = float(reference)
    except (TypeError, ValueError):
        return ""
    if reference_f == 0.0:
        return ""
    return abs(value_f - reference_f) / abs(reference_f) * 100.0


def _row_from_result(well_id: str, config: PaperModelConfig, result: dict[str, object], reference: dict[str, object] | None) -> dict[str, object]:
    reference_value = reference.get("final_effective_efficiency") if reference else result.get("final_effective_efficiency")
    return {
        "well_id": well_id,
        "nz": config.nz,
        "ny": config.ny,
        "dt": config.dt,
        "final_effective_efficiency": result.get("final_effective_efficiency"),
        "cbl_eval_interval_efficiency": result.get("cbl_eval_interval_efficiency"),
        "channeling_index": result.get("channeling_index"),
        "mixing_index": result.get("mixing_index"),
        "instability_index": result.get("instability_index"),
        "runtime_s": result.get("runtime_s"),
        "reference_case": reference is None,
        "relative_error_to_reference_percent": _relative_error(result.get("final_effective_efficiency"), reference_value),
    }


def run_grid_dt_validation(
    wells: Sequence[str],
    *,
    quick: bool = False,
    root: Path = RESULTS_ROOT,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Run grid convergence and timestep sensitivity tables."""

    out_dir = root / "04_grid_dt_validation"
    out_dir.mkdir(parents=True, exist_ok=True)
    base = default_paper_model_config()
    grid_schemes = ((60, 12), (80, 16)) if quick else ((250, 20), (500, 40), (750, 60))
    dt_schemes = (8.0, 10.0) if quick else (2.0, 4.0, 8.0)
    grid_rows: list[dict[str, object]] = []
    timestep_rows: list[dict[str, object]] = []

    for well_id in wells:
        reference_result: dict[str, object] | None = None
        for idx, (nz, ny) in enumerate(grid_schemes):
            config = replace(base, nz=nz, ny=ny, dt=8.0 if quick else 4.0, mode="grid_quick" if quick else "grid")
            result = run_single_paper_well(
                well_id,
                config,
                root=out_dir / "case_outputs" / well_id / f"grid_nz{nz}_ny{ny}_dt{config.dt:g}",
            )
            if idx == 0:
                reference_result = result
                grid_rows.append(_row_from_result(well_id, config, result, None))
            else:
                grid_rows.append(_row_from_result(well_id, config, result, reference_result))

        timestep_reference: dict[str, object] | None = None
        for idx, dt in enumerate(dt_schemes):
            config = replace(base, nz=80 if quick else 500, ny=16 if quick else 40, dt=dt, mode="dt_quick" if quick else "dt")
            result = run_single_paper_well(
                well_id,
                config,
                root=out_dir / "case_outputs" / well_id / f"timestep_nz{config.nz}_ny{config.ny}_dt{dt:g}",
            )
            if idx == 0:
                timestep_reference = result
                timestep_rows.append(_row_from_result(well_id, config, result, None))
            else:
                timestep_rows.append(_row_from_result(well_id, config, result, timestep_reference))

    write_csv(out_dir / "grid_convergence.csv", grid_rows, VALIDATION_COLUMNS)
    write_csv(out_dir / "timestep_sensitivity.csv", timestep_rows, VALIDATION_COLUMNS)
    lines = [
        "# Numerical reliability",
        "",
        f"- quick mode: {quick}",
        f"- grid cases: {len(grid_rows)}",
        f"- timestep cases: {len(timestep_rows)}",
        "",
        "Grid and timestep rows use the first case per well as the reference for relative-error reporting.",
    ]
    (out_dir / "numerical_reliability.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    write_markdown_table(out_dir / "grid_convergence.md", grid_rows, VALIDATION_COLUMNS)
    write_markdown_table(out_dir / "timestep_sensitivity.md", timestep_rows, VALIDATION_COLUMNS)
    return grid_rows, timestep_rows


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Run paper grid/time-step validation")
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--wells", nargs="+", default=["hu101", "hu103"])
    parser.add_argument("--output-root", type=Path, default=RESULTS_ROOT)
    args = parser.parse_args(argv)
    grid_rows, timestep_rows = run_grid_dt_validation(args.wells, quick=args.quick, root=args.output_root)
    print(f"Wrote grid rows: {len(grid_rows)}, timestep rows: {len(timestep_rows)}")


if __name__ == "__main__":
    main()
