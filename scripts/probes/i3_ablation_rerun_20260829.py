"""I3 修复后 R2 消融重测（2026-08-29）。

背景：T1-3b 已实现 I3 浮力弥散通量物理系数直驱；本脚本在修复后口径重跑
HT1-004 基线消融（R0→R3，CFL off、nz=500、dt=4），确认 R2 vs R1 是否仍近零
（修复前 +2e-10）；并补跑 enable_local_i3=True 的 R2/R3 对比，检验 I3 局部化
后是否生效。

输出：
    results/i3_ablation_rerun_2026-08-29/i3_ablation_rerun.csv
    results/i3_ablation_rerun_2026-08-29/README.md

用法：
    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/i3_ablation_rerun_20260829.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from cemdisp.runners.ht1_004_ablation import (
    ABLATION_LEVELS,
    extract_ablation_metrics,
    run_one_level,
)

OUTPUT_DIR = _PROJECT_ROOT / "results" / "i3_ablation_rerun_2026-08-29"
NZ = 500
DT = 4.0

CSV_COLUMNS = [
    "run_id", "ablation_level", "enable_local_i3", "nz", "dt",
    "effective_efficiency", "channeling_index", "mixing_index",
    "cement_occupation", "instability_index", "buoyancy_number",
    "elapsed_s",
]


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    totals: list[float] = []

    runs = [(lv, False) for lv in ABLATION_LEVELS] + [
        (ABLATION_LEVELS[2], True),   # R2 + local_i3
        (ABLATION_LEVELS[3], True),   # R3 + local_i3
    ]
    for lv, local_i3 in runs:
        run_id = f"ht1_004_{'li3' if local_i3 else 'base'}_{lv.name}"
        print(f"=== 运行 {run_id}（nz={NZ}, dt={DT}, CFL off）===")
        t0 = time.perf_counter()
        result = run_one_level(
            lv,
            nz=NZ,
            dt=DT,
            output_dir=str(OUTPUT_DIR),
            run_id=run_id,
            enable_local_i3=local_i3,
        )
        elapsed = time.perf_counter() - t0
        totals.append(elapsed)
        metrics = extract_ablation_metrics(result)
        row = {
            "run_id": run_id,
            "ablation_level": lv.name,
            "enable_local_i3": local_i3,
            "nz": NZ,
            "dt": DT,
            **metrics,
            "elapsed_s": round(elapsed, 1),
        }
        rows.append(row)
        print(f"--- {run_id}: eta_E={metrics['effective_efficiency']:.6f} "
              f"耗时 {elapsed:.1f} s")

    # 写 CSV
    csv_path = OUTPUT_DIR / "i3_ablation_rerun.csv"
    _write_csv(rows, csv_path)
    print(f"CSV 已写出：{csv_path}")
    print(f"6 次运行总耗时 {sum(totals):.1f} s")
    print("完成。")


def _write_csv(rows: list[dict], csv_path: Path) -> None:
    """写消融结果 CSV（utf-8-sig 便于 Excel 直接打开）。"""
    import csv

    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
