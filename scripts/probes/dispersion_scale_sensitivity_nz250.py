"""NZ=250 弥散敏感性附表（C1 方案 A）：8 井 × dispersion_dt_scale ∈ {0.0, 0.5}。

口径：修正配置（CORRECTED_KW）+ nz=250 生产网格 + CFL on，仅改弥散 dt 归一
标度 dispersion_dt_scale（基线 = 1.0，直接复用 A3 基线 JSON，不重跑）。
scale=0.0 表示关闭弥散项的 dt 归一贡献（κ→数值正则化参数声明见 obsidian
弥散敏感性附表.md）；⚠️ scale=0 时 HT1-004 会塌到很低（冒烟网格曾 1.7%），
如实记录——这正是"弥散主导"证据，论文口径必须注明是 nz=250 生产网格数字。

输出：
  results/最终基线_2026-08-29/dispersion_sens_nz250/<well>_scale{0.0|0.5}.json  # 断点续跑
  results/最终基线_2026-08-29/dispersion_scale_sensitivity_nz250.csv            # 汇总表
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SCRIPTS_DIR.parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from _sensitivity_common import load_or_run  # noqa: E402
from rerun_all_wells_corrected import CORRECTED_KW, WELLS  # noqa: E402

BASELINE_DIR = _PROJECT_ROOT / "results" / "最终基线_2026-08-29"
CASE_DIR = BASELINE_DIR / "dispersion_sens_nz250"
CSV_PATH = BASELINE_DIR / "dispersion_scale_sensitivity_nz250.csv"

SCALES = [0.0, 0.5]  # 1.0 = A3 基线（复用），此处只跑 0.0 与 0.5 两档
CSV_FIELDS = ["well", "dispersion_dt_scale", "eta_E", "eta_N", "mixing", "channeling",
              "cement_occ", "elapsed_s", "source"]


def baseline_row(well_id: str) -> dict:
    """读 A3 基线单井 JSON 作为 scale=1.0 档（source=baseline_reuse）。"""
    p = BASELINE_DIR / "cfl_on" / "单井结果" / f"{well_id}_corrected_on.json"
    base = json.loads(p.read_text(encoding="utf-8"))
    keep = ("eta_E", "eta_N", "mixing", "channeling", "cement_occ")
    return {"well": well_id, **{k: float(base[k]) for k in keep},
            "elapsed_s": float(base.get("elapsed_s", 0.0)), "source": "baseline_reuse"}


def main() -> None:
    rows: list[dict] = []
    for well_id, _loader in WELLS:
        print(f"\n=== {well_id} ===", flush=True)
        # scale=1.0：基线复用
        b = baseline_row(well_id)
        rows.append({**b, "dispersion_dt_scale": 1.0})
        print(f"  [基线复用] scale=1.0 eta_E={b['eta_E']:.4f}", flush=True)
        # scale=0.0 / 0.5：实跑
        for scale in SCALES:
            case_path = CASE_DIR / f"{well_id}_scale{scale:.1f}.json"
            row = load_or_run(
                case_path, well_id,
                compute={"extra_kw": {"dispersion_dt_scale": scale}},
                describe=f"{well_id} dispersion_dt_scale={scale}")
            rows.append({**row, "dispersion_dt_scale": scale,
                         "source": row.get("source", "computed")})

    rows.sort(key=lambda r: (list(w for w, _ in WELLS).index(r["well"]),
                             r["dispersion_dt_scale"]))
    with open(CSV_PATH, "w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nCSV: {CSV_PATH.relative_to(_PROJECT_ROOT)}", flush=True)

    # 汇总打印
    print("\n=== 弥散标度敏感性（nz=250 生产网格） ===")
    print(f"  {'井':<10}{'scale=0.0':>12}{'scale=0.5':>12}{'scale=1.0(基线)':>16}")
    for well_id, _ in WELLS:
        rs = {r["dispersion_dt_scale"]: r for r in rows if r["well"] == well_id}
        print(f"  {well_id:<10}{rs[0.0]['eta_E']:>12.4f}{rs[0.5]['eta_E']:>12.4f}"
              f"{rs[1.0]['eta_E']:>16.4f}")


if __name__ == "__main__":
    main()
