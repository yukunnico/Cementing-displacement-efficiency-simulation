"""Run paper D2DGA on/off ablation cases."""

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
    smoke_paper_model_config,
)

ABLATION_COLUMNS = [
    "well_id",
    "enable_d2dga",
    "final_effective_efficiency",
    "cbl_eval_interval_efficiency",
    "channeling_index",
    "mixing_index",
    "instability_index",
    "runtime_s",
]


def run_d2dga_ablation(
    wells: Sequence[str],
    config: PaperModelConfig,
    *,
    root: Path = RESULTS_ROOT,
) -> list[dict[str, object]]:
    """Run on/off D2DGA cases and write comparison tables."""

    rows: list[dict[str, object]] = []
    out_dir = root / "03_d2dga_ablation"
    out_dir.mkdir(parents=True, exist_ok=True)
    for well_id in wells:
        for enabled in (True, False):
            case_config = replace(config, enable_d2dga=enabled)
            case_root = out_dir / "case_outputs" / well_id / ("d2dga_on" if enabled else "d2dga_off")
            result = run_single_paper_well(well_id, case_config, root=case_root)
            rows.append(
                {
                    "well_id": well_id,
                    "enable_d2dga": enabled,
                    "final_effective_efficiency": result.get("final_effective_efficiency"),
                    "cbl_eval_interval_efficiency": result.get("cbl_eval_interval_efficiency"),
                    "channeling_index": result.get("channeling_index"),
                    "mixing_index": result.get("mixing_index"),
                    "instability_index": result.get("instability_index"),
                    "runtime_s": result.get("runtime_s"),
                }
            )
    write_csv(out_dir / "d2dga_on_off_comparison.csv", rows, ABLATION_COLUMNS)
    write_markdown_table(out_dir / "d2dga_on_off_comparison.md", rows, ABLATION_COLUMNS)
    return rows


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Run paper D2DGA ablation")
    parser.add_argument("--mode", choices=["smoke", "formal"], default="smoke")
    parser.add_argument("--wells", nargs="+", default=["hu101", "hu103"])
    parser.add_argument("--output-root", type=Path, default=RESULTS_ROOT)
    args = parser.parse_args(argv)
    config = smoke_paper_model_config() if args.mode == "smoke" else default_paper_model_config()
    rows = run_d2dga_ablation(args.wells, config, root=args.output_root)
    print(f"Wrote D2DGA ablation rows: {len(rows)}")


if __name__ == "__main__":
    main()
