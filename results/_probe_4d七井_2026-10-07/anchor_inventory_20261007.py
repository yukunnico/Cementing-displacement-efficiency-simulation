#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Phase 4d 主批前置：七井温度锚点可用性盘点（只读，确定性覆盖）

目的（计划 §4「4d 主批」第 2 条）：三重点井（hu101/ht1_003/ht1_004）的 k 与
锚点族已由用户裁定（k=0.90/0.85/0.85；族=场景 A 电测/作业史）；**其余各井**
须先盘点「锚点可用性」，无锚井一律 Geothermal 回退（口径声明），**不硬造锚**。

本盘点**不做系数取值裁定**：逐行列出各井带温度行的 `md_m / temperature_c /
温度族 / 出处 / notes 原文`，并按 notes 机读提取「明示温度系数」候选区间；
凡不能归族或候选系数越出 [0.80, 0.95] 观测带者，一律如实标出（对应计划
§8-12 锚点无法归族即停、§8-13 k 越带即停）。

输入（只读）：`参考文档/现场资料提取/<目录>/temperature_pressure_profile.csv`
（逐井列名三套 schema + hu1 异构，见 load_rows 内注释）。

输出（本目录，确定性覆盖）：
- 七井锚点盘点表.md   人读汇总（逐井带温行明细 + 归族 + 明示系数 + 可构造性）
- 七井锚点盘点表.csv  长表（井/行号/md/T/族/出处/notes原文/明示系数）

运行：
  PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -I anchor_inventory_20261007.py
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]                      # repo 根（cement model_温压耦合分支）
sys.path.insert(0, str(ROOT))               # 前插分支根：确保 import 分支内 cemdisp

from cemdisp.data.temperature_field import (  # noqa: E402
    GEO_GRAD_C_PER_M,
    GEO_T0_C,
)

# 八井全景（canonical id → 现场资料目录）；4d 批井位裁定见报告 §0
WELLS = {
    "hu101": "hu101_呼101",
    "hu102": "hu102_呼102",
    "hu103": "hu103_呼103",
    "hu1": "hu1_呼探1",
    "hu2": "ht1_002_呼探1-002",              # hu2 ≡ 呼探1-002（HT1-002）
    "ht1_001": "ht1_001_呼探1-001",
    "ht1_003": "ht1_003_呼1-003",
    "ht1_004": "ht1_004_呼1-004",
}
WELL_LABEL = {
    "hu101": "呼101", "hu102": "呼102", "hu103": "呼103", "hu1": "呼探1",
    "hu2": "呼探1-002", "ht1_001": "呼探1-001", "ht1_003": "呼1-003",
    "ht1_004": "呼1-004",
}
# 4d 主批井位（计划 §4-4d-2 括号内 5 名 vs「其余四井」不自洽；按「hu103 退出重点井、
# 新批一律不跑」裁定取 8 井 − hu103 = 七井，见 七井锚点盘点表.md §0）
BATCH_WELLS = ("hu101", "ht1_003", "ht1_004", "hu102", "hu1", "hu2", "ht1_001")
FIELD_DIR = ROOT / "参考文档" / "现场资料提取"

# 明示温度系数正则（notes 原文，三类写法都覆盖）：
#   ① 系数0.85 / 温度系数0.85        ② 0.9x温度系数
#   ③ 循环温度=156x0.9=140C / =150*0.85=127.5C（乘号后紧跟等号/单位者 = 系数）
R_COEF_WORD = re.compile(r"系数\s*(\d\.\d+)")
R_COEF_PRE = re.compile(r"(\d\.\d+)\s*(?:x|×|\*)\s*温度系数")
R_COEF_MUL = re.compile(r"[x×*]\s*(\d\.\d+)\s*(?:=|C|℃)")
# 同行「静止/循环」比值（经验观测，最可辩护的 k 主材）
R_PAIR = re.compile(r"静止/循环(?:温度)?\s*(\d+(?:\.\d+)?)\s*/\s*(\d+(?:\.\d+)?)")
R_STAT = re.compile(r"静止(?:温度)?\s*(\d+(?:\.\d+)?)")
R_CIRC = re.compile(r"循环(?:温度)?\s*(\d+(?:\.\d+)?)")

# 温度族关键词（与既有 probe 同源，按序判定）
R_NEIGHBOR = re.compile(r"邻井")
R_OUTFLOW = re.compile(r"出口")
R_EXPT = re.compile(r"实验")
R_ELEC = re.compile(r"电测")
R_STATIC_WORD = re.compile(r"静止|BHST")
R_MD_LEAD = re.compile(r"^\s*(\d+(?:\.\d+)?)")


def geo_t(md_m: float) -> float:
    return GEO_T0_C + GEO_GRAD_C_PER_M * md_m


def _clean_md(raw: str) -> float | None:
    """hu1 的 md 列含中文括注（如「5694(套管鞋)」「7601(井底)」）⇒ 取前导数值。"""
    m = R_MD_LEAD.match(raw or "")
    return float(m.group(1)) if m else None


def load_rows(well_key: str) -> list[dict]:
    """读单井 CSV。

    逐井列名不统一：标准三套 + hu1 异构（无 well_id/source_location/confidence，
    md 列名 measured_depth_m 且含中文括注，温度列可含非数值串）。
    """
    path = FIELD_DIR / WELLS[well_key] / "temperature_pressure_profile.csv"
    rows: list[dict] = []
    with path.open(encoding="utf-8-sig", newline="") as fh:
        for idx, r in enumerate(csv.DictReader(fh), start=1):
            md_raw = r.get("md_m") or r.get("measured_depth_m") or ""
            src = (r.get("source_description") or "").strip()
            if not src:
                src = (r.get("source_file") or "").strip()
                loc = (r.get("source_location") or "").strip()
                if loc:
                    src = f"{src} §{loc}"
            t_raw = (r.get("temperature_c") or "").strip()
            try:
                t_val = float(t_raw) if t_raw else None
            except ValueError:
                t_val = None          # hu1 的「表面温度」等非数值串
            rows.append(
                {
                    "row": idx,
                    "well": well_key,
                    "md_raw": md_raw,
                    "md_m": _clean_md(md_raw),
                    "T_c": t_val,
                    "T_raw": t_raw,
                    "data_type": (r.get("data_type") or r.get("type") or "").strip(),
                    "source": src,
                    "confidence": (r.get("confidence") or "").strip(),
                    "notes": (r.get("notes") or "").strip(),
                }
            )
    return rows


def static_and_circ_refs(notes: str) -> tuple[list[float], list[float]]:
    pairs = R_PAIR.findall(notes)
    srefs = [float(s) for s, _ in pairs]
    cvals = [float(c) for _, c in pairs]
    stripped = R_PAIR.sub(" ", notes)
    srefs += [float(x) for x in R_STAT.findall(stripped)]
    cvals += [float(x) for x in R_CIRC.findall(stripped)]
    dedup = lambda seq: list(dict.fromkeys(seq))  # noqa: E731
    return dedup(srefs), dedup(cvals)


def classify(row: dict) -> str:
    """按 notes/source 关键词归温度族（按序、确定性）；匹配不上 ⇒ 无法归族（不硬塞）。

    判定序：邻井 > 出口 > 实验 > 静止/BHST > 电测 > 无法归族。
    说明：末两级是「行的读数本身即静温读数」的措辞判据（如 `电测温度@7120m`
    这类列值不在 notes 里复现的行——不靠数值复现，靠措辞）；`无法归族` 一律
    如实保留，由裁定环节决定取舍。
    """
    t, notes = row["T_c"], row["notes"]
    if t is None:
        return "非温度行"
    if R_NEIGHBOR.search(notes) or R_NEIGHBOR.search(row["source"]):
        return "邻井"
    if R_OUTFLOW.search(notes):
        return "出口温度"
    if R_EXPT.search(notes):
        return "静温(实验)"
    if R_STATIC_WORD.search(notes):
        return "静温(静止/设计)"
    if R_ELEC.search(notes):
        return "静温(电测)"
    return "无法归族"


def explicit_coeffs(notes: str) -> list[float]:
    """自 notes 机读「明示温度系数」候选；结果一律过滤到 [0.50, 1.00] 合理带内。"""
    out: list[float] = []
    out += [float(v) for v in R_COEF_WORD.findall(notes)]
    out += [float(v) for v in R_COEF_PRE.findall(notes)]
    out += [float(v) for v in R_COEF_MUL.findall(notes)]
    return [c for c in dict.fromkeys(out) if 0.50 <= c <= 1.00]


def pair_ratios(notes: str) -> list[float]:
    """同行「静止/循环 a/b」的经验比值 b/a（k 的最可辩护主材）。"""
    out = []
    for a, b in R_PAIR.findall(notes):
        af, bf = float(a), float(b)
        if af > 0 and bf > 0 and bf <= af:
            out.append(round(bf / af, 4))
    return list(dict.fromkeys(out))


def main() -> None:
    all_rows: list[dict] = []
    per_well: dict[str, dict] = {}

    for wk in WELLS:
        rows = load_rows(wk)
        for r in rows:
            r["family"] = classify(r)
            r["coeffs"] = explicit_coeffs(r["notes"])
            all_rows.append(r)

        temp_rows = [r for r in rows if r["T_c"] is not None]
        static_rows = [
            r for r in temp_rows
            if r["family"].startswith("静温") and r["family"] != "静温(实验)"
        ]
        expt_rows = [r for r in temp_rows if r["family"] == "静温(实验)"]
        # 场景A 静温锚（去重 md；同 md 异值 = 族混装冲突）
        anchors: dict[float, list[float]] = {}
        for r in static_rows:
            if r["md_m"] is None:
                continue
            anchors.setdefault(float(r["md_m"]), [])
            if r["T_c"] not in anchors[float(r["md_m"])]:
                anchors[float(r["md_m"])].append(r["T_c"])
        conflicts = {md: v for md, v in anchors.items() if len(v) > 1}
        coeffs = sorted({c for r in rows for c in r["coeffs"]})
        ratios = sorted({c for r in rows for c in pair_ratios(r["notes"])})

        per_well[wk] = {
            "n_rows": len(rows),
            "n_temp": len(temp_rows),
            "static_rows": static_rows,
            "expt_rows": expt_rows,
            "anchors": anchors,
            "conflicts": conflicts,
            "coeffs": coeffs,
            "ratios": ratios,
            "constructible": (len(anchors) - len(conflicts)) >= 2,
        }

    # ---------------- md 报告 ----------------
    out = HERE.with_name("七井锚点盘点表.md")
    L: list[str] = []
    L.append("# 七井温度锚点可用性盘点（Phase 4d 主批前置，2026-10-07）\n")
    L.append("> **本表不做系数取值裁定。** 三重点井 k 已由用户裁定（呼101=0.90、"
             "呼1-003=0.85、呼1-004=0.85，族=场景 A）；其余各井的锚点可用性与"
             "候选系数一律**待用户裁定/追认**（计划 §6-6 先报后动、§8-12/§8-13）。\n")
    L.append("## 0. 井位裁定与「七井」口径\n")
    L.append("- **八井全集**（loader + 现场资料目录一一对应）："
             "hu101/hu102/hu103/hu1/**hu2**/ht1_001/ht1_003/ht1_004；"
             "其中 **hu2 ≡ 呼探1-002（HT1-002）**，现场资料目录名为 `ht1_002_呼探1-002`"
             "（`hu2_loader.py` docstring 实证），仓库内无独立「呼2」井。\n")
    L.append("- **计划 §4-4d-2 原文**：「其余四井（hu102/hu103/hu1/hu2/ht1_001）」"
             "——「四井」与括号内 5 名不自洽。按已裁定「呼103 退出重点井——**新批一律不跑**」"
             "（总纲 §1 裁定表），取 **七井 = 八井 − hu103** = "
             f"{'、'.join(WELL_LABEL[w] for w in BATCH_WELLS)}（括号内 hu103 为多列项）。\n")
    L.append("- 依据锚：AnchoredProfileField 静温锚分段线性（k 必给无默认）、"
             "同 md 异 T 构造期拒绝、域外回退统一地温式 "
             f"T(z)={GEO_T0_C}+{GEO_GRAD_C_PER_M}·z。\n")

    L.append("## 1. 归族总览\n")
    L.append("| 井 | 数据行 | 带温度行 | 静温族(场景A)行 | 唯一静温 md | 同md冲突 | 可构造(≥2锚) | 明示系数(notes) |")
    L.append("|---|---|---|---|---|---|---|---|")
    for wk in WELLS:
        d = per_well[wk]
        conf = "、".join(f"{md:.1f}" for md in d["conflicts"]) or "—"
        coef = "、".join(f"{c}" for c in d["coeffs"]) or "—"
        mark = "**批外(hu103)**" if wk not in BATCH_WELLS else ""
        L.append(
            f"| {WELL_LABEL[wk]} ({wk}) {mark} | {d['n_rows']} | {d['n_temp']} | "
            f"{len(d['static_rows'])} | {len(d['anchors'])} | {conf} | "
            f"{'是' if d['constructible'] else '**否**'} | {coef} |"
        )
    L.append("")

    L.append("## 2. 逐井带温度行明细（原文可核）\n")
    for wk in WELLS:
        L.append(f"### {WELL_LABEL[wk]}（{wk}）\n")
        L.append("| 行 | md_m(原文) | T°C | data_type | 族 | 地温式°C | 残差°C | 出处 | notes 原文 |")
        L.append("|---|---|---|---|---|---|---|---|---|")
        for r in all_rows:
            if r["well"] != wk or r["T_c"] is None:
                continue
            md = r["md_m"]
            geo = f"{geo_t(md):.2f}" if md is not None else "—"
            res = f"{r['T_c'] - geo_t(md):+.2f}" if md is not None else "—"
            L.append(
                f"| {r['row']} | {r['md_raw']} | {r['T_c']} | {r['data_type']} | "
                f"{r['family']} | {geo} | {res} | {r['source']} | {r['notes']} |"
            )
        L.append("")
        d = per_well[wk]
        anchors = d["anchors"]
        if anchors:
            L.append("- **场景A 静温锚候选集**："
                     + "、".join(f"{md:.3f} m@{v[0]}°C" for md, v in sorted(anchors.items())
                                 if len(v) == 1))
            L.append("")
        if d["conflicts"]:
            for md, v in sorted(d["conflicts"].items()):
                L.append(f"- ⚠ **同 md 冲突** {md:.3f} m：{' / '.join(f'{x}°C' for x in v)}"
                         "（AnchoredProfileField 构造期将拒绝——族混装须先裁定）")
            L.append("")
        if d["expt_rows"]:
            L.append("- 实验族行（场景 A 不取）："
                     + "；".join(f"{r['md_raw']}@{r['T_c']}°C" for r in d["expt_rows"]))
            L.append("")

    L.append("## 3. 明示温度系数候选（notes 机读，无量纲乘子）与同行经验比值\n")
    L.append("| 井 | notes 明示系数 | 同行 静止/循环 经验比值 | 全落 [0.80,0.95] 观测带 | 依据原文 |")
    L.append("|---|---|---|---|---|")
    for wk in WELLS:
        d = per_well[wk]
        if not d["coeffs"] and not d["ratios"]:
            L.append(f"| {WELL_LABEL[wk]} | — | — | — | （无系数/成对行）|")
            continue
        inband = all(0.80 <= c <= 0.95 for c in d["coeffs"] + d["ratios"])
        refs = "；".join(
            r["notes"] for r in all_rows if r["well"] == wk and (r["coeffs"] or r["notes"])
        )
        L.append(
            f"| {WELL_LABEL[wk]} | {'、'.join(str(c) for c in d['coeffs']) or '—'} | "
            f"{'、'.join(str(c) for c in d['ratios']) or '—'} | "
            f"{'是' if inband else '**否**'} | {refs} |"
        )
    L.append("")

    L.append("## 4. 4d 主批拟用锚点集与候选系数（**待用户裁定/追认**）\n")
    L.append("> 三重点井 = 已裁定行，非裁定项；其余各井的 k 与「是否上锚定场」"
             "= 计划 §6-6 先报后动硬停点，本表只列材料、不代裁。\n")
    L.append("| 井 | 拟用静温锚点集（场景A） | 明示系数 | 经验比值 | 建议 k（待裁） | 无锚回退 |")
    L.append("|---|---|---|---|---|---|")
    for wk in WELLS:
        d = per_well[wk]
        pts = "、".join(f"{md:.3f}@{v[0]}°C" for md, v in sorted(d["anchors"].items())
                        if len(v) == 1) or "—"
        if wk == "hu103":
            L.append(f"| {WELL_LABEL[wk]}（批外） | {pts} | "
                     f"{'、'.join(str(c) for c in d['coeffs']) or '—'} | "
                     f"{'、'.join(str(c) for c in d['ratios']) or '—'} | **不跑（退出重点井）** | — |")
            continue
        if not d["constructible"]:
            L.append(f"| {WELL_LABEL[wk]} | {pts} | — | — | **不可构造（<2 锚）** | "
                     f"**Geothermal 回退** |")
            continue
        kex = "、".join(str(c) for c in d["coeffs"]) or "—"
        kr = "、".join(str(c) for c in d["ratios"]) or "—"
        if wk in ("ht1_003", "ht1_004"):
            sug = "0.85（已裁定）"
        elif wk == "hu101":
            sug = "0.90（已裁定）"
        elif len(d["coeffs"]) == 1:
            sug = f"{d['coeffs'][0]}（notes 明示）"
        elif len(d["coeffs"]) > 1:
            sug = "**多值冲突**（" + "、".join(str(c) for c in d["coeffs"]) + "）"
        else:
            sug = "待裁（无明示系数）"
        L.append(f"| {WELL_LABEL[wk]} | {pts} | {kex} | {kr} | {sug} | — |")
    L.append("")

    L.append("## 5. 复跑\n")
    L.append("```\nPYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -I "
             "results/_probe_4d七井_2026-10-07/anchor_inventory_20261007.py\n```\n")
    out.write_text("\n".join(L), encoding="utf-8")

    # ---------------- csv 长表 ----------------
    csv_path = HERE.with_name("七井锚点盘点表.csv")
    with csv_path.open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["well", "well_label", "row", "md_raw", "md_m", "T_c", "data_type",
                    "family", "geo_t_c", "residual_c", "explicit_coeffs", "source", "notes"])
        for r in all_rows:
            md = r["md_m"]
            w.writerow([
                r["well"], WELL_LABEL[r["well"]], r["row"], r["md_raw"],
                "" if md is None else f"{md:.6f}",
                "" if r["T_c"] is None else r["T_c"],
                r["data_type"], r["family"],
                "" if md is None else f"{geo_t(md):.6f}",
                "" if (md is None or r["T_c"] is None) else f"{r['T_c'] - geo_t(md):.6f}",
                ",".join(str(c) for c in r["coeffs"]), r["source"], r["notes"],
            ])

    print(f"[inventory] 写 {out.name} / {csv_path.name}")
    for wk in WELLS:
        d = per_well[wk]
        print(f"  {WELL_LABEL[wk]:<10} 带温 {d['n_temp']:>2} 行 | 静温锚 "
              f"{len(d['anchors']):>2} | 冲突 {len(d['conflicts'])} | "
              f"可构造={'Y' if d['constructible'] else 'N'} | 明示系数 {d['coeffs']}")


if __name__ == "__main__":
    main()
