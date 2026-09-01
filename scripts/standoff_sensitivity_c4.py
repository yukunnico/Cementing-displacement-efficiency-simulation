"""C4: standoff ±0.1 敏感性（修正后配置，nz=250, CFL on）——论文素材。

对 hu101（有实测剖面）、ht1_001、ht1_003（覆盖中低/高 η_E 代表井）各跑三档：
standoff 基线剖面逐点 -0.1 / 0.0 / +0.1（clip 到 [0,1]），输出 η_E/η_N 敏感性区间。

口径：
- 与最终基线同款配置（CORRECTED_KW + nz=250 + CFL on，见 _sensitivity_common）；
- 0.0 档默认复用 A3 基线 JSON（results/最终基线_2026-08-29/cfl_on/单井结果/
  <well>_corrected_on.json），--rerun-zero 可强制重跑验证确定性；
- hu101 用名义 standoff 剖面（0.38–0.48，loader 默认 measured_standoff=None，
  与基线一致），非实测剖面口径。

输出：
  results/最终基线_2026-08-29/standoff_sens_0.1/<well>_delta{-0.1|+0.1}.json  # 断点续跑
  results/最终基线_2026-08-29/standoff_sensitivity_0.1.csv                    # 汇总表
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SCRIPTS_DIR.parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from _sensitivity_common import load_or_run, shift_standoff  # noqa: E402
import cemdisp.data.loaders as L  # noqa: E402

BASELINE_DIR = _PROJECT_ROOT / "results" / "最终基线_2026-08-29"
CASE_DIR = BASELINE_DIR / "standoff_sens_0.1"
CSV_PATH = BASELINE_DIR / "standoff_sensitivity_0.1.csv"

WELLS = ["hu101", "ht1_001", "ht1_003"]
DELTAS = [-0.1, 0.0, 0.1]
LOADERS = {
    "hu101": L.load_hu101_tailpipe,
    "ht1_001": L.load_ht1_001_tailpipe,
    "ht1_003": L.load_ht1_003_tailpipe,
}

CSV_FIELDS = ["well", "standoff_delta", "eta_E", "eta_N", "mixing", "channeling",
              "cement_occ", "elapsed_s", "source"]


def baseline_row(well_id: str) -> dict:
    """读 A3 基线单井 JSON 作为 0.0 档（source=baseline_reuse）。"""
    p = BASELINE_DIR / "cfl_on" / "单井结果" / f"{well_id}_corrected_on.json"
    base = json.loads(p.read_text(encoding="utf-8"))
    keep = ("eta_E", "eta_N", "mixing", "channeling", "cement_occ")
    return {"well": well_id, **{k: float(base[k]) for k in keep},
            "elapsed_s": float(base.get("elapsed_s", 0.0)), "source": "baseline_reuse"}


def main() -> None:
    ap = argparse.ArgumentParser(description="C4: standoff ±0.1 敏感性（3 井×3 档）")
    ap.add_argument("--rerun-zero", action="store_true",
                    help="0.0 档不使用基线复用，强制重跑（确定性验证用）")
    args = ap.parse_args()

    rows: list[dict] = []
    for well_id in WELLS:
        print(f"\n=== {well_id} ===", flush=True)
        well, fluids, schedule, _ = LOADERS[well_id]()  # 基线口径加载（名义 standoff）
        for delta in DELTAS:
            if delta == 0.0 and not args.rerun_zero:
                row = baseline_row(well_id)
                print(f"  [基线复用] delta=0.0 eta_E={row['eta_E']:.4f}", flush=True)
            else:
                case_path = CASE_DIR / f"{well_id}_delta{delta:+.1f}.json"
                row = load_or_run(case_path, well_id,
                                  compute={"well_override": shift_standoff(well, delta)},
                                  describe=f"{well_id} standoff{delta:+.1f}")
            rows.append({**row, "standoff_delta": delta,
                         "source": row.get("source", "computed")})

    with open(CSV_PATH, "w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nCSV: {CSV_PATH.relative_to(_PROJECT_ROOT)}", flush=True)

    # 敏感性区间摘要
    print("\n=== η_E 敏感性区间（standoff ±0.1） ===")
    for well_id in WELLS:
        rs = {r["standoff_delta"]: r for r in rows if r["well"] == well_id}
        e_lo, e0, e_hi = rs[-0.1]["eta_E"], rs[0.0]["eta_E"], rs[0.1]["eta_E"]
        print(f"  {well_id:<8} η_E: {e_lo:.4f} (−0.1) → {e0:.4f} (0) → {e_hi:.4f} (+0.1)"
              f"   区间 [{min(e_lo, e_hi):.4f}, {max(e_lo, e_hi):.4f}]"
              f"   全宽 {abs(e_hi - e_lo):+.4f}")


if __name__ == "__main__":
    main()
