"""理论算例：放大偏心度到 35%/45%，展示 I3/真体力的物理改进能力。
HT1-004 现场偏心固定 17%（standoff 83%），改进效果受限，故设理论算例。（Task 6a Step 2）

运行 R0 和 R3（levels 0, 3）在 e=0.17, 0.35, 0.45 三个偏心度下。
预期：高偏心度下 R3 相对 R0 效率差异应明显大于现场 17% 偏心。
"""
from __future__ import annotations

import csv
import dataclasses
import json
import os
from pathlib import Path

from cemdisp.data.loaders.ht1_004_loader import load_ht1_004_tailpipe
from cemdisp.data.well_spec import DepthValuePoint, WellSpec
from cemdisp.runners.ht1_004_ablation import (
    ABLATION_LEVELS,
    extract_ablation_metrics,
    run_one_level,
)

OUTPUT_DIR = "results/ht1_004_ablation"
CSV_PATH = os.path.join(OUTPUT_DIR, "ablation_summary.csv")
CSV_COLUMNS = [
    "run_id", "ablation_level", "eccentricity", "nz", "dt",
    "effective_efficiency", "channeling_index", "mixing_index",
    "cement_occupation", "instability_index", "buoyancy_number",
]


def _append_csv_row(row: dict, csv_path: str, columns: list[str]) -> None:
    """Append a single row to the CSV; write header if file is new."""
    Path(csv_path).parent.mkdir(parents=True, exist_ok=True)
    write_header = not os.path.exists(csv_path)
    with open(csv_path, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def _dump_json(result, run_id: str, output_dir: str) -> None:
    """Dump per-run metrics JSON."""
    p = Path(output_dir) / f"{run_id}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    payload = extract_ablation_metrics(result)
    payload["run_id"] = run_id
    p.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8",
    )


def _build_well_with_eccentricity(well: WellSpec, eccentricity: float) -> WellSpec:
    """Rebuild WellSpec with standoff = 1 - eccentricity at all depth points."""
    standoff = 1.0 - eccentricity
    depths = [p.depth_md_m for p in well.standoff_profile]
    new_profile = tuple(DepthValuePoint(depth_md_m=d, value=standoff) for d in depths)
    return dataclasses.replace(well, standoff_profile=new_profile)


if __name__ == "__main__":
    nz, dt = 280, 4.0
    eccentricities = [0.17, 0.35, 0.45]
    levels_to_run = [ABLATION_LEVELS[0], ABLATION_LEVELS[3]]  # R0 and R3 only

    # Load well once to get the original depth profile
    well_orig, _, _, _ = load_ht1_004_tailpipe()

    print(f"=== Theoretical Eccentricity Cases (nz={nz}, dt={dt}) ===\n")
    print(f"{'e':<8} {'Level':<6} {'eff':>8} {'channeling':>12} {'mixing':>10} {'b':>10}")
    print("-" * 60)

    for e in eccentricities:
        well_new = _build_well_with_eccentricity(well_orig, e)
        for lv in levels_to_run:
            run_id = f"theoretical_e{int(e*100)}_{lv.name}"
            print(f"  Running {run_id} ...")
            result = run_one_level(
                lv, nz=nz, dt=dt, total_t=None,
                output_dir=OUTPUT_DIR,
                well_spec_override=well_new,
            )
            m = extract_ablation_metrics(result)
            b = m.get("buoyancy_number")
            b_str = f"{b:.2f}" if isinstance(b, (int, float)) else "N/A"
            print(
                f"{e:<8.2f} {lv.name:<6} {m['effective_efficiency']:>8.4f} "
                f"{m['channeling_index']:>12.6f} {m['mixing_index']:>10.6f} {b_str:>10}"
            )

            # Dump JSON
            _dump_json(result, f"theoretical_e{int(e*100)}_{lv.name}", OUTPUT_DIR)

            # Append CSV
            row = {
                "run_id": run_id,
                "ablation_level": lv.name,
                "eccentricity": e,
                "nz": nz,
                "dt": dt,
                **m,
            }
            _append_csv_row(row, CSV_PATH, CSV_COLUMNS)

    print(f"\nCSV appended to: {CSV_PATH}")
    print("Done.")