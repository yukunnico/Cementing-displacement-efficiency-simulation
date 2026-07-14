"""网格收敛 + 时间步敏感性验证（论文图9/表6）。（Task 6a Step 4）

运行 R3 (level 3) at:
- Grid: nz=140, 280, 500 (dt=4.0)
- Time step: dt=2.0, 4.0, 8.0 (nz=280)

预期: nz=280 vs 500 效率变化 <2%; dt=2 vs 4 变化 <1%。
"""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

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

R3 = ABLATION_LEVELS[3]


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


if __name__ == "__main__":
    print("=== Grid Convergence (R3, dt=4.0) ===\n")
    grid_nzs = [140, 280, 500]
    print(f"{'nz':<8} {'eff':>8} {'channeling':>12} {'mixing':>10} {'b':>10}")
    print("-" * 50)
    for nz in grid_nzs:
        run_id = f"convergence_nz{nz}_dt4.0"
        print(f"  Running {run_id} ...")
        result = run_one_level(R3, nz=nz, dt=4.0, total_t=None, output_dir=OUTPUT_DIR)
        m = extract_ablation_metrics(result)
        b = m.get("buoyancy_number")
        b_str = f"{b:.2f}" if isinstance(b, (int, float)) else "N/A"
        print(
            f"{nz:<8} {m['effective_efficiency']:>8.4f} "
            f"{m['channeling_index']:>12.6f} {m['mixing_index']:>10.6f} {b_str:>10}"
        )
        _dump_json(result, run_id, OUTPUT_DIR)
        row = {
            "run_id": run_id,
            "ablation_level": "R3",
            "eccentricity": 0.17,
            "nz": nz,
            "dt": 4.0,
            **m,
        }
        _append_csv_row(row, CSV_PATH, CSV_COLUMNS)

    dts = [2.0, 4.0, 8.0]
    print(f"\n=== Time-Step Convergence (R3, nz=280) ===\n")
    print(f"{'dt':<8} {'eff':>8} {'channeling':>12} {'mixing':>10} {'b':>10}")
    print("-" * 50)
    for dt in dts:
        run_id = f"convergence_nz280_dt{dt}"
        print(f"  Running {run_id} ...")
        result = run_one_level(R3, nz=280, dt=dt, total_t=None, output_dir=OUTPUT_DIR)
        m = extract_ablation_metrics(result)
        b = m.get("buoyancy_number")
        b_str = f"{b:.2f}" if isinstance(b, (int, float)) else "N/A"
        print(
            f"{dt:<8.1f} {m['effective_efficiency']:>8.4f} "
            f"{m['channeling_index']:>12.6f} {m['mixing_index']:>10.6f} {b_str:>10}"
        )
        _dump_json(result, run_id, OUTPUT_DIR)
        row = {
            "run_id": run_id,
            "ablation_level": "R3",
            "eccentricity": 0.17,
            "nz": 280,
            "dt": dt,
            **m,
        }
        _append_csv_row(row, CSV_PATH, CSV_COLUMNS)

    print(f"\nCSV appended to: {CSV_PATH}")
    print("Done.")