"""HT1-004 R0->R3 消融运行入口。产出 results/ht1_004_ablation/ 下的各级摘要JSON + 汇总CSV。（Task 6a）"""
from __future__ import annotations

import os
from pathlib import Path

from cemdisp.runners.ht1_004_ablation import (
    run_full_ablation,
    extract_ablation_metrics,
    ABLATION_LEVELS,
    ABLATION_CSV_COLUMNS,
    append_ablation_csv_row,
)

OUTPUT_DIR = "results/ht1_004_ablation"
CSV_PATH = os.path.join(OUTPUT_DIR, "ablation_summary.csv")


if __name__ == "__main__":
    nz, dt = 500, 4.0
    print(f"=== HT1-004 Full Ablation R0->R3 (nz={nz}, dt={dt}) ===\n")
    results = run_full_ablation(
        nz=nz, dt=dt, total_t=None, output_dir=OUTPUT_DIR,
        run_id_prefix="ht1_004_ablation",
    )

    # Print summary table
    print(f"{'Level':<6} {'eff':>8} {'channeling':>12} {'mixing':>10} {'b':>10}")
    print("-" * 50)
    for lv in ABLATION_LEVELS:
        res = results[lv.name]
        m = extract_ablation_metrics(res)
        b = m.get("buoyancy_number")
        b_str = f"{b:.2f}" if isinstance(b, (int, float)) else "N/A"
        print(
            f"{lv.name:<6} {m['effective_efficiency']:>8.4f} "
            f"{m['channeling_index']:>12.6f} {m['mixing_index']:>10.6f} {b_str:>10}"
        )

        # Append to CSV
        row = {
            "run_id": f"full_R{lv.name}_{nz}_{dt}",
            "ablation_level": lv.name,
            "eccentricity": 0.17,
            "nz": nz,
            "dt": dt,
            **m,
        }
        append_ablation_csv_row(row, CSV_PATH)

    print(f"\nCSV appended to: {CSV_PATH}")
    print("Done.")