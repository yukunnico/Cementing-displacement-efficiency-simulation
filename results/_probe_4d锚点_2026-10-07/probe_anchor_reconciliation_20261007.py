#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Phase 4d 前置对账 probe：三井静温锚点对账（可重跑驱动）

目的（Phase 4 设计规格 §2-4d 第一步）：为「温度系数 0.85–0.90 逐井取值 =
先报后动硬停点」准备裁定材料。**本 probe 不做系数取值裁定**，只产出
按温度族分列的对账表、地温式残差、AnchoredProfileField k 扫描逐锚残差，
以及每井候选 k 区间（建议性质，标注"待用户裁定，禁自行取值"）。

输入（只读）：参考文档/现场资料提取/<hu101_呼101|ht1_003_呼1-003|ht1_004_呼1-004>/
temperature_pressure_profile.csv（四族口径混装，测绘项16；逐行按 notes/source
关键词归族，无法归族者如实列出）。

输出（本目录，确定性覆盖）：
- 三井锚点对账表.md      人读汇总表（族分类/地温残差/k扫描/候选区间/质量注记）
- 三井锚点对账表.csv     锚点长表（井/md/实测T/族/出处/notes原文/地温式/候选静温线A|B|C）
- 三井锚点_ks扫描.csv    井×场景×k×锚点 残差长表

运行：PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -I probe_anchor_reconciliation_20261007.py
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]                     # repo 根（cement model_温压耦合分支）
sys.path.insert(0, str(ROOT))              # 确保 import 分支内 cemdisp（非冻结源目录）

from cemdisp.data.temperature_field import (  # noqa: E402
    GEO_GRAD_C_PER_M,
    GEO_T0_C,
    AnchoredProfileField,
)

WELLS = {
    "hu101": "hu101_呼101",
    "ht1_003": "ht1_003_呼1-003",
    "ht1_004": "ht1_004_呼1-004",
}
K_GRID = [0.80, 0.85, 0.875, 0.90, 0.95, 1.0]
FIELD_DIR = ROOT / "参考文档" / "现场资料提取"


def geo_t(md_m: float) -> float:
    """统一地温式 T(z)=16.006+1.7598e-2·z（09-30 §1 裁定；复用模块单一真源常量）。"""
    return GEO_T0_C + GEO_GRAD_C_PER_M * md_m


# ---------------------------------------------------------------------------
# 读表（列名逐井不统一：hu101 用 source_file+source_location，其余 source_description）
# ---------------------------------------------------------------------------
def load_rows(well_key: str) -> list[dict]:
    path = FIELD_DIR / WELLS[well_key] / "temperature_pressure_profile.csv"
    rows = []
    with path.open(encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh):
            src = r.get("source_description") or ""
            if not src:
                src = (r.get("source_file") or "") + (
                    (" §" + r["source_location"]) if r.get("source_location") else ""
                )
            t = r.get("temperature_c", "").strip()
            md = float(r["md_m"])
            rows.append(
                {
                    "well": well_key,
                    "md_m": md,
                    "T_c": float(t) if t else None,
                    "data_type": r.get("data_type", ""),
                    "source": src,
                    "confidence": r.get("confidence", ""),
                    "notes": r.get("notes", "") or "",
                }
            )
    return rows


# ---------------------------------------------------------------------------
# 温度族归族规则（机读、确定性；关键词 = notes/source 原文）
# ---------------------------------------------------------------------------
R_NEIGHBOR = re.compile(r"邻井")
R_OUTFLOW = re.compile(r"出口")
R_EXPT = re.compile(r"实验温度")
R_ELEC = re.compile(r"电测")
R_WORKH = re.compile(r"作业史")
# 静止↔循环 数值：成对（静止/循环 150/129）或分写（静止温度 152C循环129.2C）
R_PAIR = re.compile(r"静止/循环(?:温度)?\s*(\d+(?:\.\d+)?)\s*/\s*(\d+(?:\.\d+)?)")
R_STAT = re.compile(r"静止(?:温度)?\s*(\d+(?:\.\d+)?)")
R_CIRC = re.compile(r"循环(?:温度)?\s*(\d+(?:\.\d+)?)")
R_TEMPX = re.compile(r"温度\s*(\d+(?:\.\d+)?)\s*(?:℃|C)")
R_ELEC_AFTER = re.compile(r"电测[^\d]{0,15}(\d+(?:\.\d+)?)")
R_ELEC_BEFORE = re.compile(r"(\d+(?:\.\d+)?)[^\d。；]{0,6}电测")


def static_and_circ_refs(notes: str) -> tuple[list[float], list[float]]:
    """从 notes 提取静止参考值与循环值（先摘走 静止/循环 a/b 成对，避免误捕）。"""
    pairs = R_PAIR.findall(notes)
    srefs = [float(s) for s, _ in pairs]
    cvals = [float(c) for _, c in pairs]
    stripped = R_PAIR.sub(" ", notes)
    srefs += [float(x) for x in R_STAT.findall(stripped)]
    cvals += [float(x) for x in R_CIRC.findall(stripped)]
    dedup = lambda seq: list(dict.fromkeys(seq))  # noqa: E731
    return dedup(srefs), dedup(cvals)


def classify(row: dict) -> str:
    """按列值与 notes 关键词归温度族；匹配不上 ⇒ 无法归族（不硬塞）。"""
    t, notes = row["T_c"], row["notes"]
    if t is None:
        return "非温度行"
    if R_NEIGHBOR.search(notes) or R_NEIGHBOR.search(row["source"]):
        return "邻井"
    if R_OUTFLOW.search(notes):
        return "出口温度"
    srefs, _c = static_and_circ_refs(notes)
    markers = []
    if R_EXPT.search(notes) and any(abs(t - s) < 1e-9 for s in srefs):
        markers.append("实验")
    if R_WORKH.search(notes) and any(abs(t - s) < 1e-9 for s in srefs):
        markers.append("作业史")
    if R_ELEC.search(notes):
        elec_vals = [float(x) for x in R_ELEC_AFTER.findall(notes)]
        elec_vals += [float(x) for x in R_ELEC_BEFORE.findall(notes)]
        if any(abs(t - s) < 1e-9 for s in elec_vals):
            markers.append("电测")
    if markers:
        return "静温(" + "+".join(markers) + ")"
    return "无法归族"


# 无法归族行的推断线索（仅注记，不并入列值归族结果；场景C 中显式标“推断”）
INFER_HINT = {
    ("hu101", 5400.0): "notes 空；与 5700 行同源同节(20124.doc 2.9)，疑为回接段电测静止（推断）",
    ("ht1_003", 6500.0): "列值 152 与 5305 行实验区间(5305-6500m 静止 152C)一致，疑为实验静温区间终点（推断）",
    ("ht1_004", 6600.0): "列值 155 与 5241 行实验区间(5241-6600m 静止 155C)一致，疑为实验静温区间终点（推断）",
    ("ht1_004", 7660.0): "notes 为地层压力系数内容（无温度族措辞）；与作业史(八)行同 md 同值，疑为静止族重复行（推断）",
}


def build_anchor_table(well_key: str, rows: list[dict]) -> list[dict]:
    out = []
    for r in rows:
        if r["T_c"] is None:
            continue  # 纯压力/当量行不进锚点对账（计数在别处报告）
        fam = classify(r)
        srefs, cvals = static_and_circ_refs(r["notes"])
        out.append(
            {
                **r,
                "family": fam,
                "geo_T": geo_t(r["md_m"]),
                "resid_vs_geo": (r["T_c"] - geo_t(r["md_m"])),
                "circ_refs": cvals,
                "stat_refs": srefs,
                "infer_hint": (INFER_HINT.get((well_key, r["md_m"]), "")
                                   if fam == "无法归族" else ""),
            }
        )
    return sorted(out, key=lambda a: a["md_m"])


def circ_anchors(anchors: list[dict]) -> list[dict]:
    """循环族次级锚点（列值全为空，值从 notes 提取；md=所在行 md）。"""
    out = []
    for a in anchors:
        for c in a["circ_refs"]:
            out.append(
                {
                    "well": a["well"],
                    "md_m": a["md_m"],
                    "T_c": c,
                    "family": "循环(notes次级)",
                    "row_family": a["family"],
                    "row_T": a["T_c"],
                    "source": a["source"],
                    "notes": a["notes"],
                }
            )
    return out


# ---------------------------------------------------------------------------
# 候选静温锚集（场景 A/B/C —— 族取舍本身也是待裁定项，probe 并列呈现不代裁）
# ---------------------------------------------------------------------------
def scenarios(anchors: list[dict]) -> dict[str, list[tuple[float, float]]]:
    A = [
        (a["md_m"], a["T_c"])
        for a in anchors
        if a["family"].startswith("静温") and "实验" not in a["family"]
    ]
    B = []
    for a in anchors:
        if "实验" in a["family"]:
            B.append((a["md_m"], a["T_c"]))
        elif a["family"] == "无法归族" and "实验静温区间终点" in a["infer_hint"]:
            B.append((a["md_m"], a["T_c"]))  # 区间终点推断行，如实参与场景B
    C = [
        (a["md_m"], a["T_c"])
        for a in anchors
        if a["family"].startswith("静温")
        or (a["family"] == "无法归族" and "静止" in a["infer_hint"])
        or (a["family"] == "无法归族" and "静止" in a["infer_hint"] + a["notes"])
    ]
    # C = 全体静温族（含推断行）；同 md 同值重复行交由 AnchoredProfileField 合并
    ded = {}
    for md, t in C:
        ded.setdefault(md, []).append(t)
    C2 = []
    for md, tv in sorted(ded.items()):
        if len(set(tv)) == 1:
            C2.append((md, tv[0]))
        else:
            C2.extend((md, x) for x in tv)  # 保留冲突原样 → 类构造将实证拒绝
    return {"A电测/作业史静温": A, "B实验静温+区间终点": B, "C全静温并集(含推断)": C2}


R_COEF = re.compile(r"温度系数\s*(0?\.\d+|\d+(?:\.\d+)?)")


def row_ratios(anchors: list[dict]) -> list[dict]:
    """同行 静止↔循环 比值（k 经验观测量，与场景构造无关）。"""
    out = []
    for a in anchors:
        if not a["family"].startswith("静温"):
            continue
        for c in a["circ_refs"]:
            out.append(
                {
                    "well": a["well"],
                    "md_m": a["md_m"],
                    "static_T": a["T_c"],
                    "circ_T": c,
                    "ratio": c / a["T_c"],
                    "design_coef": (R_COEF.findall(a["notes"]) or [""])[0],
                    "family": a["family"],
                    "notes": a["notes"],
                }
            )
    return out


def k_scan(well_key: str, anchors: list[dict]) -> tuple[list[dict], list[dict]]:
    """每场景×k：AnchoredProfileField(circulating) 在全部带值锚点上的逐锚残差。"""
    targets = [(a["md_m"], a["T_c"], a["family"], a["infer_hint"]) for a in anchors]
    for c in circ_anchors(anchors):
        targets.append((c["md_m"], c["T_c"], c["family"], ""))
    rows, implied = [], []
    for name, pts in scenarios(anchors).items():
        if len(pts) < 2:
            rows.append({"well": well_key, "scenario": name, "k": None,
                         "note": f"锚点不足2点（n={len(pts)}），无法构造场",
                         "target_md": None})
            continue
        try:
            line_f = AnchoredProfileField(pts, 1.0)  # k=1 ⇒ 静温线
        except ValueError as exc:
            rows.append({"well": well_key, "scenario": name, "k": None,
                         "note": f"构造拒绝（族混装实证）：{exc}",
                         "target_md": None})
            continue
        for c in targets:
            if c[2] != "循环(notes次级)":
                continue
            lv = line_f.T(c[0], 0.0)
            implied.append({"well": well_key, "scenario": name, "md_m": c[0],
                            "circ_T": c[1], "static_line_T": lv,
                            "implied_k": (c[1] / lv) if lv else float("nan")})
        for k in K_GRID:
            f = AnchoredProfileField(pts, k)
            for md, t, fam, hint in targets:
                rows.append({
                    "well": well_key, "scenario": name, "k": k,
                    "target_md": md, "target_T": t, "target_family": fam,
                    "T_field": f.T(md, 0.0),
                    "residual": f.T(md, 0.0) - t,
                    "note": hint,
                })
    return rows, implied


# ---------------------------------------------------------------------------
# 输出：CSV（机读） + MD（人读）
# ---------------------------------------------------------------------------
def build_line_fields(per_well: dict) -> dict:
    """每井每场景的静温线场（k=1.0）；返回 {well: {scen: field|str}}（不可构造为说明）。"""
    out = {}
    for wk, anchors in per_well.items():
        m = {}
        for name, pts in scenarios(anchors).items():
            try:
                m[name] = AnchoredProfileField(pts, 1.0)
            except ValueError as exc:
                m[name] = f"不可构造：{exc}"
        out[wk] = m
    return out


def write_recon_csv(path: Path, per_well: dict, line_fields: dict) -> list[dict]:
    ratios = []
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "well", "md_m", "T_c", "family", "geo_T", "resid_vs_geo",
            "static_line_A", "static_line_B", "static_line_C",
            "circ_refs_notes", "infer_hint", "source", "notes",
        ])
        for wk, anchors in per_well.items():
            for a in anchors:
                cols = []
                for name in line_fields[wk]:
                    lf = line_fields[wk][name]
                    cols.append(lf.T(a["md_m"], 0.0) if hasattr(lf, "T") else "")
                w.writerow([
                    wk, a["md_m"], a["T_c"], a["family"],
                    round(a["geo_T"], 4), round(a["resid_vs_geo"], 4),
                    *(round(c, 4) if isinstance(c, float) else "" for c in cols),
                    ";".join(str(x) for x in a["circ_refs"]),
                    a["infer_hint"], a["source"], a["notes"],
                ])
        for wk, anchors in per_well.items():
            for c in circ_anchors(anchors):
                w.writerow([
                    c["well"], c["md_m"], c["T_c"], c["family"], "", "",
                    "", "", "", "", f"所在行族={c['row_family']}；所在行列值={c['row_T']}",
                    c["source"], c["notes"],
                ])
            ratios.extend(row_ratios(anchors))
    return ratios


def write_scan_csv(path: Path, scan_rows: list[dict], implied: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["record_type", "well", "scenario", "k", "target_md",
                    "target_T", "target_family", "T_field", "residual",
                    "static_line_T", "implied_k", "note"])
        for r in scan_rows:
            w.writerow(["k_scan", r.get("well"), r.get("scenario"), r.get("k"),
                        r.get("target_md", ""), r.get("target_T", ""),
                        r.get("target_family", ""), r.get("T_field", ""),
                        r.get("residual", ""), "", "", r.get("note", "")])
        for r in implied:
            w.writerow(["implied_k", r["well"], r["scenario"], "", r["md_m"],
                        r["circ_T"], "循环(notes次级)", "", "",
                        round(r["static_line_T"], 4), round(r["implied_k"], 4), ""])


def fmt(x, nd=2):
    return "" if x is None or x == "" else (f"{x:.{nd}f}" if isinstance(x, (int, float)) else str(x))


def scan_summary(scan_rows: list[dict]) -> dict:
    """{(well, scenario): {k: {max静温, mean静温, max循环, mean循环, n静温, n循环}}}"""
    agg: dict = {}
    for r in scan_rows:
        if r.get("k") is None:
            agg.setdefault((r["well"], r["scenario"]), {"note": r["note"]})
            continue
        key = (r["well"], r["scenario"])
        d = agg.setdefault(key, {})
        kd = d.setdefault(r["k"], {"res_s": [], "res_c": []})
        if r["target_family"].startswith("静温"):
            kd["res_s"].append(r["residual"])
        elif r["target_family"] == "循环(notes次级)":
            kd["res_c"].append(r["residual"])
    return agg


def write_md(path: Path, per_well: dict, ratios: list[dict], scan_rows: list[dict],
             implied: list[dict]) -> None:
    L: list[str] = []
    A = L.append
    A("# 三井静温锚点对账表（Phase 4d 前置 probe，2026-10-07）")
    A("")
    A("> **本表不做任何系数取值裁定。**温度系数 k 逐井取值（0.85–0.90 候选）= 先报后动硬停点"
      "（Phase 4 设计规格 §2-4d / §8-6）；以下候选区间均为**建议性质，待用户裁定，禁自行取值**。")
    A("")
    A("## 0. 口径声明")
    A("")
    A("- 输入：`参考文档/现场资料提取/<井>/temperature_pressure_profile.csv`（工作区未跟踪副本，"
      "内容与 415ed26 版本一致；hu101/ht1_003 仅差 BOM）。")
    A("- AnchoredProfileField 语义：`T(md,t) = k · 静温锚分段线`（regime=circulating，"
      "t 缺省不消费）；k=1 ⇒ 纯静温锚。静温↔循环为**乘性系数**口径（notes 实证："
      "152×0.85=129.2；155×0.85≈131.75≈132）。")
    A("- 残差 = 模型/回退线 − 实测（°C）。地温式 = T(z)=16.006+1.7598e-2·z（09-30 §1 统一式）。")
    A("- 场景 A/B/C 的**族取舍本身待裁定**，本表并列呈现、不代裁；C 保留同 md 冲突原样，"
      "由 AnchoredProfileField 构造拒绝以实证族混装校验。")
    A("- k 扫描网格：0.80 / 0.85 / 0.875 / 0.90 / 0.95 / 1.0（1.0=纯静温参照）。")
    A("")
    A("## 1. 归族总览")
    A("")
    A("| 井 | 数据行 | 带温度行 | 静温族 | 出口 | 邻井 | 无法归族 | 循环(notes次级) | 纯压力/当量行 |")
    A("|---|---|---|---|---|---|---|---|---|")
    for wk, dirn in WELLS.items():
        rows = load_rows(wk)
        an = per_well[wk]
        nstat = sum(1 for a in an if a["family"].startswith("静温"))
        A(f"| {wk} | {len(rows)} | {len(an)} | {nstat} | "
          f"{sum(1 for a in an if a['family']=='出口温度')} | "
          f"{sum(1 for a in an if a['family']=='邻井')} | "
          f"{sum(1 for a in an if a['family']=='无法归族')} | "
          f"{len(circ_anchors(an))} | {sum(1 for r in rows if r['T_c'] is None)} |")
    A("")
    A("## 2. 锚点明细（按族分列；含出处列原文）")
    for wk in WELLS:
        A("")
        A(f"### {wk}")
        A("")
        A("| md_m | 实测T°C | 族 | 地温式°C | 残差(实测−地温)°C | 推断注记 | 出处 | notes 原文 |")
        A("|---|---|---|---|---|---|---|---|")
        for a in per_well[wk]:
            A(f"| {a['md_m']:.0f} | {a['T_c']:.0f} | {a['family']} | "
              f"{a['geo_T']:.2f} | {a['resid_vs_geo']:+.2f} | {a['infer_hint']} | "
              f"{a['source']} | {a['notes']} |")
        for c in circ_anchors(per_well[wk]):
            A(f"| {c['md_m']:.0f} | {c['T_c']:.1f} | 循环(notes次级，非列值) | — | — | "
              f"所在行族={c['row_family']}，列值={c['row_T']} | {c['source']} | {c['notes']} |")
    A("")
    A("## 3. 静温族锚 vs 地温式残差（重点列）")
    A("")
    A("| 井 | md_m | 实测T°C | 族 | 地温式°C | 残差°C |")
    A("|---|---|---|---|---|---|")
    for wk in WELLS:
        for a in per_well[wk]:
            if a["family"].startswith("静温"):
                A(f"| {wk} | {a['md_m']:.0f} | {a['T_c']:.0f} | {a['family']} | "
                  f"{a['geo_T']:.2f} | {a['resid_vs_geo']:+.2f} |")
    A("")
    A("## 4. 同行 静止↔循环 比值（k 经验观测，与场景构造无关）")
    A("")
    A("| 井 | md_m | 静止T°C | 循环T°C | 比值=circ/static | notes 明示系数 |")
    A("|---|---|---|---|---|---|")
    for r in ratios:
        A(f"| {r['well']} | {r['md_m']:.0f} | {r['static_T']:.0f} | {r['circ_T']:.1f} | "
          f"{r['ratio']:.3f} | {r['design_coef']} |")
    A("")
    A("## 5. AnchoredProfileField k 扫描逐锚残差（摘要；全量见 三井锚点_ks扫描.csv）")
    A("")
    A("> 目标集=全部带列值行+循环次级锚；**族外锚**（如以电测线衡量实验族锚）的残差同样计入"
      "「静温锚残差」列——族间残差即族选择敏感性的量化，非缺陷。")
    agg = scan_summary(scan_rows)
    for (wk, scen), d in agg.items():
        A("")
        A(f"### {wk} — 场景 {scen}")
        if "note" in d:
            A(f"（{d['note']}）")
            continue
        A("")
        A("| k | 静温锚残差 max/mean °C | 循环锚残差 max/mean °C |")
        A("|---|---|---|")
        for k in K_GRID:
            kd = d.get(k)
            if not kd:
                continue
            rs, rc = kd["res_s"], kd["res_c"]
            s_txt = (f"{max(abs(x) for x in rs):.2f} / {sum(abs(x) for x in rs)/len(rs):.2f}"
                     if rs else "—")
            c_txt = (f"{max(abs(x) for x in rc):.2f} / {sum(abs(x) for x in rc)/len(rc):.2f}"
                     if rc else "—")
            A(f"| {k} | {s_txt} | {c_txt} |")
    A("")
    A("## 6. 每井候选 k 区间（**建议性质——待用户裁定，禁自行取值**）")
    A("")
    A("| 井 | 同行观测比值范围 | notes 明示设计系数 | 场景 implied_k 范围（构造可行者） | 建议候选区间 |")
    A("|---|---|---|---|---|")
    for wk in WELLS:
        rr = [r["ratio"] for r in ratios if r["well"] == wk]
        coefs = sorted({r["design_coef"] for r in ratios if r["well"] == wk and r["design_coef"]})
        im = [x["implied_k"] for x in implied if x["well"] == wk]
        im_txt = (f"{min(im):.3f}–{max(im):.3f}" if im else "—")
        s_txt = (f"[{min(rr):.2f}, {max(rr):.2f}]" if rr else "—")
        cand = "待用户裁定（见下方逐井建议）"
        A(f"| {wk} | {s_txt}（n={len(rr)}） | {'、'.join(coefs) if coefs else '无'} | "
          f"{im_txt} | {cand} |")
    best: dict = {}
    for (wk2, scen), d in agg.items():
        if "note" in d:
            continue
        c2 = [
            (sum(abs(x) for x in kd["res_c"]) / len(kd["res_c"]), k)
            for k in K_GRID if (kd := d.get(k)) and kd["res_c"]
        ]
        if c2:
            best.setdefault(wk2, []).append((scen, min(c2)[1]))
    A("")
    for wk in WELLS:
        rr = [r["ratio"] for r in ratios if r["well"] == wk]
        rng = (f"[{min(rr):.3f}, {max(rr):.3f}]" if rr else "—")
        arg = "、".join(f"{s_} 于 k={k:.3f}" for s_, k in best.get(wk, [])) or "—"
        A(f"- **{wk}**：同行观测区间 {rng}；场景 argmin(循环锚平均残差)：{arg}。"
          f"候选建议={rng}（**建议性质，待用户裁定，禁自行取值**；是否取中值/设计值/分段待定）。")
    A("- 注：场景 implied_k 范围受族配对支配（如实验静止的循环值对电测线量得 k>1/k<0.7 伪影），"
      "裁定主材应以 §4 同行比值与本区间为准。")
    A("")
    A("## 7. 数据质量注记（如实列出，不硬塞）")
    A("")
    A("1. **无法归族行**（有列值但 notes 无族措辞/措辞与列值不一致；§2 有原文）：")
    for wk in WELLS:
        for a in per_well[wk]:
            if a["family"] == "无法归族":
                A(f"   - {wk} md={a['md_m']:.0f} m T={a['T_c']:.0f} °C —— {a['infer_hint'] or '无可靠推断依据'}")
    A("2. **同 md 双族冲突**：ht1_004 md=5241 m 电测静止 124 °C vs 实验静止 155 °C（相差 31 °C）。"
      "场景 C 构造被 AnchoredProfileField 拒绝（族混装构造期校验实证）；"
      "ht1_003 5290(电测123)↔5305(实验152) 15 m 内 29 °C 跳变同源——**静温族选谁直接改变基准线**。")
    A("3. **循环族零列值行**：三井 CSV 的 temperature_c 列无一为循环温度值；§2/§4 循环锚均为 "
      "notes 正则次级提取（静止/循环 a/b 对、循环(温度)x、实验区间分写式），md=所在行 md。")
    A("4. **出口温度 0 m @60 °C**：井口出口（循环洗井两周后）条件，非地层静温，不入静温锚、"
      "不参与 k 扫描目标（仅在 §2 列出备查）。")
    A("5. **邻井行**（ht1_004 内 7746 m@150=HT1-001、7618 m@152=HT1-003，confidence=medium）："
      "属他井域且 7746 m 超过本井 TD 7660 m——一律排除出本井场锚点。")
    A("6. **列名逐井不统一**：hu101 用 source_file+source_location，ht1_003/004 用 source_description；"
      "已统一读为 source 列。hu101/ht1_003 文件带 BOM。")
    A("7. **文件版本**：三 CSV 现均未被 git 跟踪（3664c02 起 参考文档 整树移出版本控制）；"
      "工作区内容与 415ed26 版本一致（hu101/ht1_003 仅差 BOM，ht1_004 逐字节一致），"
      "与测绘项16 行统计（003=8/5、004=12/8 含邻井2、101=7/3）相符。")
    A("8. **notes 明示系数**：ht1_003/ht1_004 实验行均写“领浆/尾浆温度系数 0.85”；hu101 无明示系数行，"
      "仅有 静止150↔循环135（0.90）与 静止123↔循环99（0.805）两组比值对。")
    A("")
    A("## 8. 复跑")
    A("")
    A("```")
    A("PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -I "
      "results/_probe_4d锚点_2026-10-07/probe_anchor_reconciliation_20261007.py")
    A("```")
    A("")
    A("*输出文件：本目录 `三井锚点对账表.md/.csv` + `三井锚点_ks扫描.csv`。"
      "本 probe 由 Phase 4d 前期任务生成；系数取值决定权在用户（先报后动）。*")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> int:
    per_well = {wk: build_anchor_table(wk, load_rows(wk)) for wk in WELLS}
    line_fields = build_line_fields(per_well)
    out = HERE.parent
    ratios = write_recon_csv(out / "三井锚点对账表.csv", per_well, line_fields)
    scan_rows: list[dict] = []
    implied: list[dict] = []
    for wk in WELLS:
        rs, im = k_scan(wk, per_well[wk])
        scan_rows += rs
        implied += im
    write_scan_csv(out / "三井锚点_ks扫描.csv", scan_rows, implied)
    write_md(out / "三井锚点对账表.md", per_well, ratios, scan_rows, implied)
    print("== 三井静温锚点对账 probe 完成 ==")
    for wk in WELLS:
        an = per_well[wk]
        print(f"{wk}: 带温度行 {len(an)} | 静温族 {sum(1 for a in an if a['family'].startswith('静温'))} | "
              f"无法归族 {sum(1 for a in an if a['family']=='无法归族')} | "
              f"循环次级 {len(circ_anchors(an))}")
    print(f"k-scan 行数 {len(scan_rows)}；implied_k 行数 {len(implied)}；"
          f"同行比值 {len(ratios)} 组")
    print(f"输出目录：{out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
