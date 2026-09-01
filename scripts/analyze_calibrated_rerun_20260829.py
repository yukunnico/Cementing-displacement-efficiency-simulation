"""校准后重跑 vs 校准前基线 + CBL 对比分析（2026-08-29）。"""
from __future__ import annotations
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NEW = ROOT / "results" / "校准后重跑_2026-08-29" / "汇总.csv"
OLD = ROOT / "results" / "全井修正前后" / "汇总.csv"
CBL = {"hu101": 0.6277, "hu102": 0.6665, "hu103": 0.1206, "ht1_003": 0.787, "ht1_004": 0.003}
WELLS = ["hu101", "hu102", "hu103", "hu1", "hu2", "ht1_001", "ht1_003", "ht1_004"]


def load(path: Path, corrected: bool = True) -> dict:
    with open(path, encoding="utf-8-sig") as f:
        rows = [r for r in csv.DictReader(f) if (r.get("corrected") == "True") == corrected]
    return {r["well"]: r for r in rows}


def main() -> None:
    old = load(OLD)
    new = load(NEW)
    print("=== 校准后 corrected 结果（vs 校准前 corrected 基线 + CBL） ===")
    hdr = f"{'well':<9}{'eta_E':>8}{'old':>8}{'delta':>9}{'eta_N':>8}{'mixing':>8}{'chan':>7}{'instab':>7}{'wall':>7}{'CBL':>7}"
    print(hdr)
    print("-" * len(hdr))
    for w in WELLS:
        r = new.get(w)
        if not r:
            print(f"{w:<9}missing")
            continue
        o = old.get(w)
        d = float(r["eta_E"]) - float(o["eta_E"]) if o else float("nan")
        c = CBL.get(w)
        cbl_str = f"{c:.4f}" if c is not None else "     -"
        print(f"{w:<9}{float(r['eta_E']):>8.4f}{float(o['eta_E']):>8.4f}{d:>+8.4f}"
              f"{float(r['eta_N']):>8.4f}{float(r['mixing']):>8.4f}{float(r['channeling']):>7.3f}"
              f"{float(r['instability']):>7.3f}{float(r['wall_frac']):>7.3f}{cbl_str:>7}")
    print()
    print("CBL 口径注：hu101=62.77%(5390-7810 测量段)；hu102=66.65%(6840-7665)；"
          "hu103=12.06%(139.7段 7338-7712)；ht1_003=78.7%(数字化总占比口径，官方判不合格)；"
          "ht1_004=0.3%(尾管评价段 5245-7581，官方全井 29.99% 不合格)。hu1/hu2/ht1_001 无定量 CBL。")


if __name__ == "__main__":
    main()
