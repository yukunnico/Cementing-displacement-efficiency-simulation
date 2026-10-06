# -*- coding: utf-8 -*-
"""Phase 0 专项分析（2026-10-06）——本批专属产物，不写入任何既有批次目录。

产出（全部写在本文件所在目录 = 本批产物目录）：
  1. standoff_密度_Ton_vs_Toff.csv/.md
     T-on static（本批）vs T-off（既有批：呼101=09-16、呼1-004=09-27；呼1-003=无）
  2. 跨批复现校验.csv
     本批 vs T2 批同名变体：T-off 组应**逐位**（口径无关）；T-on 组差异由
     Phase 0.0「密度就近取」（Q16）解释，并列 Δ 供核对
  3. 解读数据.json
     供人工撰写「一页解读.md」的机读事实（最优点、单调性、退化档提示）

口径与红线（本脚本只读既有产物，绝不写入）：
  - T-off 批不受 Phase 0.0 影响（T-off 不调用 fluid_at）⇒ 可直接对照；
  - 表格档 r0.6/0.8 在本分析中不出现（本批不含 table 档）。
"""
from __future__ import annotations

import json
from pathlib import Path

PHASE0 = Path(__file__).resolve().parent
ROOT = PHASE0.parents[1]
T2 = ROOT / "results" / "敏感性变体_温压T2_2026-10-01"
TOFF_DIR = {
    "呼101": ROOT / "results" / "敏感性变体_当前口径_2026-09-16",
    "呼1-004": ROOT / "results" / "敏感性变体_呼1-004_2026-09-27",
}  # 呼1-003 无任何既有 T-off 资产（本批首建）

WELLS = ("呼101", "呼1-003", "呼1-004")

# 档位名映射：本批 Ton_static_standoff_<tag>_rate_x1.0 ↔ 09-16 矩阵 standoff_<tag2>
STANDOFF_MAP = (
    ("m0.30", "standoff_m0.30"), ("m0.20", "standoff_m0.20"),
    ("m0.15", "standoff_m0.15"), ("m0.10", "standoff_m0.1"),
    ("m0.05", "standoff_m0.05"), ("p0.05", "standoff_p0.05"),
    ("p0.10", "standoff_p0.1"), ("p0.20", "standoff_p0.20"),
    ("p0.30", "standoff_p0.30"),
)
DENS_MAP = (("p100", "spacer_dens_p100"), ("m100", "spacer_dens_m100"))

# SO=0.83 井上 +0.20 / +0.30 经 np.clip 后同为 1.0（变体退化，见 spec §4）
DEGENERATE = {"呼1-003": [("p0.20", "p0.30")], "呼1-004": [("p0.20", "p0.30")]}


def _final(p: Path):
    if not p.exists():
        return None
    fr = json.loads(p.read_text(encoding="utf-8"))["最终结果"]
    return float(fr["全井段最终有效顶替效率"]), float(fr["窄四分位效率"])


def _extra(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}


def _p0(well: str, variant: str):
    return _final(PHASE0 / f"{well}_{variant}_结果摘要.json")


def _p0x(well: str, variant: str) -> dict:
    return _extra(PHASE0 / f"{well}_{variant}_判别量.json")


rows_cmp: list[dict] = []
facts: dict = {"井": {}}

for well in WELLS:
    tf = TOFF_DIR.get(well)
    well_facts: dict = {"有T_off对照": tf is not None, "档": {}}

    base_off = _p0(well, "Toff_zero")
    base_on = _p0(well, "Ton_static_rate_x1.0")
    well_facts["Toff_zero"] = base_off
    well_facts["Ton_static_rate_x1.0"] = base_on

    for kind, mp in (("standoff", STANDOFF_MAP), ("密度", DENS_MAP)):
        for tag, toff_name in mp:
            v_on = (f"Ton_static_standoff_{tag}_rate_x1.0" if kind == "standoff"
                    else f"Ton_static_spacer_dens_{tag}_rate_x1.0")
            on = _p0(well, v_on)
            off = _final(tf / f"{well}_{toff_name}_结果摘要.json") if tf else None
            row = {
                "井": well, "类型": kind, "档": tag,
                "T_off_η_E": "" if off is None else off[0],
                "T_off_η_N": "" if off is None else off[1],
                "T_on_η_E": "" if on is None else on[0],
                "T_on_η_N": "" if on is None else on[1],
                "ΔηE_pp": "" if (off is None or on is None) else (on[0] - off[0]) * 100.0,
                "ΔηN_pp": "" if (off is None or on is None) else (on[1] - off[1]) * 100.0,
                "T_off_来源": "" if tf is None else tf.name,
                "备注": "",
            }
            rows_cmp.append(row)
            if on is not None:
                well_facts["档"].setdefault(kind, {})[tag] = {
                    "T_off": off, "T_on": on,
                    "ΔηN_pp": None if off is None else (on[1] - off[1]) * 100.0,
                    "ΔηE_pp": None if off is None else (on[0] - off[0]) * 100.0,
                }
            well_facts["档"].setdefault("诊断", {})[f"{kind}:{tag}"] = {
                k: _p0x(well, v_on).get(k)
                for k in ("饥饿份额", "屈服门活化率_b加权", "屈服门_wall占比",
                          "front_narrow_m", "front_wide_m")
            }
    for a, b in DEGENERATE.get(well, []):
        for r in rows_cmp:
            if r["井"] == well and r["档"] == a:
                r["备注"] = f"⚠️ 与本井 {b} 档经 clip 后退化为同一剖面（同为 SO=1.0）"
    facts["井"][well] = well_facts

XCHECK = {
    "呼101": ["Toff_zero", "Ton_const60_rate_x1.0", "Ton_static_rate_x1.0"],
    "呼1-003": [],
    "呼1-004": ["Toff_zero", "Ton_const60_rate_x1.0", "Ton_static_rate_x1.0",
                "Toff_rate_x0.8", "Toff_rate_x1.2",
                "Ton_const60_rate_x0.8", "Ton_const60_rate_x1.2",
                "Ton_static_rate_x0.8", "Ton_static_rate_x1.2"],
}
rows_x: list[dict] = []
for well, variants in XCHECK.items():
    for v in variants:
        a = _p0(well, v)
        b = _final(T2 / f"{well}_{v}_结果摘要.json")
        if a is None:
            continue
        if b is None:
            rows_x.append({"井": well, "变体": v, "T2批": "（T2 批无此档）",
                           "本批η_E": a[0], "本批η_N": a[1],
                           "ΔηE_pp": "", "ΔηN_pp": "", "位级": ""})
            continue
        de, dn = (a[0] - b[0]) * 100.0, (a[1] - b[1]) * 100.0
        bitwise = (a[0] == b[0] and a[1] == b[1])
        rows_x.append({
            "井": well, "变体": v, "T2批": "有",
            "本批η_E": a[0], "本批η_N": a[1],
            "T2η_E": b[0], "T2η_N": b[1],
            "ΔηE_pp": de, "ΔηN_pp": dn,
            "位级": "逐位相同" if bitwise else
                    ("T-off（口径无关，应逐位）" if v.startswith("Toff") else
                     "T-on（Δ 归因 Phase 0.0 密度就近取 Q16）"),
        })


def _write_csv(name: str, rows: list[dict]) -> None:
    import csv
    with (PHASE0 / name).open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


_write_csv("standoff_密度_Ton_vs_Toff.csv", rows_cmp)
_write_csv("跨批复现校验.csv", rows_x)
(PHASE0 / "解读数据.json").write_text(
    json.dumps(facts, ensure_ascii=False, indent=2), encoding="utf-8")

md = [
    "# standoff / 隔离液密度：T-on static vs T-off（2026-10-06 补跑批）",
    "",
    "口径：**T-on static** = 本批 `Ton_static_<档>_rate_x1.0`（09-16 链 + 静温剖面，"
    "**Phase 0.0 密度「就近取」口径**）；**T-off** = 既有批同名档（呼101 ← "
    "`敏感性变体_当前口径_2026-09-16`；呼1-004 ← `敏感性变体_呼1-004_2026-09-27`）。",
    "T-off 不调用 `fluid_at` ⇒ 不受密度口径变更影响，**可直接对照**。",
    "呼1-003 无任何既有 T-off 资产（本批首建）⇒ 该井只有 T-on 列。",
    "",
    "⚠️ 呼1-003 / 呼1-004 的 standoff 剖面恒 0.83，`shift_standoff` 用 `np.clip`，"
    "**`p0.20` 与 `p0.30` 档退化为同一剖面**（同为 1.0）——27 名义格 = 25 个互异剖面。",
    "",
    "| 井 | 类型 | 档 | T-off η_N | T-on η_N | Δη_N/pp | T-off η_E | T-on η_E | Δη_E/pp | 备注 |",
    "|---|---|---|---|---|---|---|---|---|---|",
]
for r in rows_cmp:
    f = lambda v: "—" if v == "" else f"{float(v):.4f}"
    g = lambda v: "—" if v == "" else f"{float(v):+.3f}"
    md.append(
        f"| {r['井']} | {r['类型']} | {r['档']} | {f(r['T_off_η_N'])} | {f(r['T_on_η_N'])} | "
        f"{g(r['ΔηN_pp'])} | {f(r['T_off_η_E'])} | {f(r['T_on_η_E'])} | {g(r['ΔηE_pp'])} | "
        f"{r['备注'] or '—'} |"
    )
(PHASE0 / "standoff_密度_Ton_vs_Toff.md").write_text("\n".join(md), encoding="utf-8")

print("\n".join(md))
print(f"\n跨批复现校验 {len(rows_x)} 行 → 跨批复现校验.csv")
print("DONE")
