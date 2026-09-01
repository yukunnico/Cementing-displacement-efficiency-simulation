"""hu101 领/尾浆主检 vs 复检双口径敏感性（2026-08-29 校准后续分析项）。

复用 scripts/rerun_all_wells_corrected.py 的官方 corrected 口径 run_one（M1+M3+I3+M4+M2、
nz=250、CFL 自适应），仅切换 load_hu101_tailpipe 的 rheology_source：
  primary = 主检 2011122.pdf（W301-22094，89C）：领 0.844/0.381、尾 0.830/0.352
  recheck = 复检 2011121.doc（93C）：领 0.719/0.815、尾 0.722/0.684
输出 results/校准后重跑_2026-08-29/hu101_流变敏感性_主检vs复检.csv
"""
from __future__ import annotations
import csv
import sys
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

import cemdisp.data.loaders as L
from scripts.rerun_all_wells_corrected import PROJECT_ROOT, run_one

OUT_DIR = PROJECT_ROOT / "results" / "校准后重跑_2026-08-29"


def main() -> None:
    rows = []
    for tag, loader in (
        ("primary", L.load_hu101_tailpipe),
        ("recheck", lambda: L.load_hu101_tailpipe(rheology_source="recheck")),
    ):
        print(f"=== hu101 rheology_source={tag} ===", flush=True)
        row = run_one(f"hu101_{tag}", loader, corrected=True, cfl_on=True)
        row["rheology_source"] = row.get("rheology_source", tag)
        row["well"] = "hu101"
        rows.append(row)
        print(f"  eta_E={row['eta_E']:.4f} eta_N={row['eta_N']:.4f} ({row['elapsed_s']}s)", flush=True)
    out = OUT_DIR / "hu101_流变敏感性_主检vs复检.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    keys = sorted({k for r in rows for k in r})
    with open(out, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"CSV: {out}", flush=True)


if __name__ == "__main__":
    main()
