"""网格收敛 + 时间步敏感性验证（论文图9/表6）。（Task 6a Step 4）

运行 R3 (level 3) at:
- Grid: nz=140, 280, 500 (dt=4.0)
- Time step: dt=2.0, 8.0 (nz=280)  — dt=4.0 skipped: already covered by grid sweep

预期: nz=280 vs 500 效率变化 <2%; dt=2 vs 4 变化 <1%。
"""
from __future__ import annotations

import os

from cemdisp.runners.ht1_004_ablation import (
    ABLATION_LEVELS,
    ABLATION_CSV_COLUMNS,
    append_ablation_csv_row,
    extract_ablation_metrics,
    run_one_level,
)

OUTPUT_DIR = "results/ht1_004_ablation"
CSV_PATH = os.path.join(OUTPUT_DIR, "ablation_summary.csv")

R3 = ABLATION_LEVELS[3]


if __name__ == "__main__":
    print("=== Grid Convergence (R3, dt=4.0) ===\n")
    grid_nzs = [140, 280, 500]
    print(f"{'nz':<8} {'eff':>8} {'channeling':>12} {'mixing':>10} {'b':>10}")
    print("-" * 50)
    for nz in grid_nzs:
        run_id = f"convergence_nz{nz}_dt4.0"
        print(f"  Running {run_id} ...")
        result = run_one_level(
            R3, nz=nz, dt=4.0, total_t=None, output_dir=OUTPUT_DIR,
            run_id=run_id, eccentricity=0.17,
        )
        m = extract_ablation_metrics(result)
        b = m.get("buoyancy_number")
        b_str = f"{b:.2f}" if isinstance(b, (int, float)) else "N/A"
        print(
            f"{nz:<8} {m['effective_efficiency']:>8.4f} "
            f"{m['channeling_index']:>12.6f} {m['mixing_index']:>10.6f} {b_str:>10}"
        )
        row = {
            "run_id": run_id,
            "ablation_level": "R3",
            "eccentricity": 0.17,
            "nz": nz,
            "dt": 4.0,
            **m,
        }
        append_ablation_csv_row(row, CSV_PATH)

    # dt=4.0 is intentionally skipped — it is already covered by the grid
    # sweep above (nz=280/dt=4.0), avoiding a duplicate CSV row.
    dts = [2.0, 8.0]
    print(f"\n=== Time-Step Convergence (R3, nz=280) ===\n")
    print(f"{'dt':<8} {'eff':>8} {'channeling':>12} {'mixing':>10} {'b':>10}")
    print("-" * 50)
    for dt in dts:
        run_id = f"convergence_nz280_dt{dt}"
        print(f"  Running {run_id} ...")
        result = run_one_level(
            R3, nz=280, dt=dt, total_t=None, output_dir=OUTPUT_DIR,
            run_id=run_id, eccentricity=0.17,
        )
        m = extract_ablation_metrics(result)
        b = m.get("buoyancy_number")
        b_str = f"{b:.2f}" if isinstance(b, (int, float)) else "N/A"
        print(
            f"{dt:<8.1f} {m['effective_efficiency']:>8.4f} "
            f"{m['channeling_index']:>12.6f} {m['mixing_index']:>10.6f} {b_str:>10}"
        )
        row = {
            "run_id": run_id,
            "ablation_level": "R3",
            "eccentricity": 0.17,
            "nz": 280,
            "dt": dt,
            **m,
        }
        append_ablation_csv_row(row, CSV_PATH)

    print(f"\nCSV appended to: {CSV_PATH}")
    print("Done.")