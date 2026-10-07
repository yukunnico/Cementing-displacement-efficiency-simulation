#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Phase 4d 主批收集器：七井温度场对账表 + 汇总表 + 一页解读（确定性覆盖）

输入：本目录 `jobs/*.json`（七井批驱动产物）。
输出：
- `七井温度场对账表.md/.csv`  逐锚点：井 / k / md / 实测T / 锚定场值(k·T) / 地温式 / 残差
- `汇总表_七井温度场.csv/.md`  井 × {Toff, Tstatic_col, Tanchored_col} 的 η_E/η_N + 关5 健康度
- `一页解读.md`                井间分化叙事（接 2e）+ 口径声明 + 对桥结论

运行：PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -I collect_seven_wells_20261007.py
"""

from __future__ import annotations

import csv
import importlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
for _p in (str(ROOT), str(ROOT / "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from cemdisp.data.temperature_field import GEO_GRAD_C_PER_M, GEO_T0_C  # noqa: E402
from entrypoints.run_sensitivity_current_20260916 import (  # noqa: E402
    ANCHORED_NO_ANCHOR_WELLS,
    ANCHORED_WELLS,
)

OUT = HERE.parent                                # results/_probe_4d七井_2026-10-07
JOBS_DIR = OUT / "jobs"
WELLS = ("呼101", "呼1-003", "呼1-004", "呼102", "呼探1", "呼探1-002", "呼探1-001")
WELL_CN = {"呼101": "呼101", "呼1-003": "呼1-003", "呼1-004": "呼1-004",
           "呼102": "呼102", "呼探1": "呼探1", "呼探1-002": "呼探1-002(hu2)",
           "呼探1-001": "呼探1-001"}
LOADERS = {
    "呼101": ("cemdisp.data.loaders.hu101_loader", "load_hu101_tailpipe"),
    "呼1-003": ("cemdisp.data.loaders.ht1_003_loader", "load_ht1_003_tailpipe"),
    "呼1-004": ("cemdisp.data.loaders.ht1_004_loader", "load_ht1_004_tailpipe"),
    "呼102": ("cemdisp.data.loaders.hu102_loader", "load_hu102_tailpipe"),
    "呼探1": ("cemdisp.data.loaders.hu1_loader", "load_hu1_tailpipe"),
    "呼探1-002": ("cemdisp.data.loaders.hu2_loader", "load_hu2_tailpipe"),
    "呼探1-001": ("cemdisp.data.loaders.ht1_001_loader", "load_ht1_001_tailpipe"),
}
GROUP_SUFFIX = (("Toff", "off"), ("Tstatic_col", "static"), ("Tanchored_col", "anchored"))


def geo_t(md: float) -> float:
    return GEO_T0_C + GEO_GRAD_C_PER_M * md


def fmt(v, nd=6):
    return "—" if v is None else f"{v:.{nd}f}"


def fmt_signed(v, nd=3):
    return "—" if v is None else f"{v:+.{nd}f}"


def load_rows() -> dict[str, dict]:
    out: dict[str, dict] = {}
    for p in sorted(JOBS_DIR.glob("*.json")):
        if p.name.startswith("_spec_"):
            continue
        out[p.stem] = json.loads(p.read_text(encoding="utf-8"))
    return out


def well_domain(wk: str) -> dict:
    mod, fn = LOADERS[wk]
    well = getattr(importlib.import_module(mod), fn)()[0]
    return {"top_md_m": float(well.top_md_m), "bottom_md_m": float(well.bottom_md_m),
            "shoe_md_m": float(well.shoe_md_m)}


def field_curve(wk: str, md: float) -> float:
    """`anchored` 档在该井该深度的场值（circulating 语义 = k·静温锚线）。"""
    if wk in ANCHORED_NO_ANCHOR_WELLS:
        return geo_t(md)
    k, anchors, _ = ANCHORED_WELLS[wk]
    return k * float(np.interp(md, [a[0] for a in anchors], [a[1] for a in anchors]))


def main() -> None:
    res = load_rows()
    missing = [f"{wk}_{s}" for wk in WELLS for s, _ in GROUP_SUFFIX
               if f"{wk}_{s}" not in res]
    if missing:
        print(f"[collector] 缺 run：{missing}", flush=True)

    # ---------------- 对账表 ----------------
    csv_rows = []
    for wk in WELLS:
        if wk not in ANCHORED_WELLS:
            continue
        k, anchors, _ = ANCHORED_WELLS[wk]
        for md, t_meas in anchors:
            csv_rows.append({
                "井": wk, "k": k, "md_m": md, "实测T_C": t_meas,
                "锚定场值_C(循环=k·实测)": round(k * t_meas, 4),
                "地温式_C": round(geo_t(md), 4),
                "实测-地温_残差C": round(t_meas - geo_t(md), 4),
                "锚定场-地温_C": round(k * t_meas - geo_t(md), 4),
            })
    with (OUT / "七井温度场对账表.csv").open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(csv_rows[0]))
        w.writeheader()
        w.writerows(csv_rows)

    L: list[str] = []
    L.append("# 七井温度场对账表（Phase 4d 主批，2026-10-07）\n")
    L.append("> 口径：`anchored` 档 = `AnchoredProfileField(regime=\"circulating\")`，"
             "循环准稳态语义 `T = k · 静温锚分段线`（k、锚点集见 spec §1，"
             "出处原文见 `七井锚点盘点表.md` §2）。**无瞬态项**（`transient_tau_s=None`）。"
             "地温式 = 统一线 `T(z)=16.006+1.7598e-2·z`（09-30 §1 裁定）。\n")
    L.append("## 1. 逐井场参数与域内取值\n")
    L.append("| 井 | 场类型 | k | 锚点数 | 域顶 md | 域底 md | T(域顶) | T(域中) | T(域底) |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for wk in WELLS:
        d = well_domain(wk)
        mid = 0.5 * (d["top_md_m"] + d["bottom_md_m"])
        ftype = "Geothermal(无锚回退)" if wk in ANCHORED_NO_ANCHOR_WELLS else "Anchored"
        k = ANCHORED_WELLS[wk][0] if wk in ANCHORED_WELLS else "—"
        n = len(ANCHORED_WELLS[wk][1]) if wk in ANCHORED_WELLS else 0
        L.append(f"| {WELL_CN[wk]} | {ftype} | {k} | {n} | {d['top_md_m']:.1f} | "
                 f"{d['bottom_md_m']:.1f} | {field_curve(wk, d['top_md_m']):.2f} | "
                 f"{field_curve(wk, mid):.2f} | {field_curve(wk, d['bottom_md_m']):.2f} |")
    L.append("")
    L.append("## 2. 逐锚点：实测 vs 锚定场 vs 地温式（**锚点热图对照**）\n")
    L.append("| 井 | k | md_m | 实测T°C | 锚定场值°C(循环=k·实测) | 地温式°C | "
             "实测−地温 残差°C | 锚定场−地温°C |")
    L.append("|---|---|---|---|---|---|---|---|")
    for r in csv_rows:
        L.append(f"| {WELL_CN[r['井']]} | {r['k']} | {r['md_m']} | {r['实测T_C']} | "
                 f"{r['锚定场值_C(循环=k·实测)']} | {r['地温式_C']} | "
                 f"{r['实测-地温_残差C']:+.2f} | {r['锚定场-地温_C']:+.2f} |")
    L.append("")
    L.append("## 3. 跨井锚点矩阵（模型「anchored」场值 °C；`—`=该井无此锚深度）\n")
    uniq_md = sorted({md for wk in WELLS if wk in ANCHORED_WELLS
                      for md, _ in ANCHORED_WELLS[wk][1]})
    L.append("| md_m | " + " | ".join(WELL_CN[w] for w in WELLS if w in ANCHORED_WELLS) + " |")
    L.append("|---" * (1 + sum(1 for w in WELLS if w in ANCHORED_WELLS)) + "|")
    for md in uniq_md:
        cells = []
        for wk in WELLS:
            if wk not in ANCHORED_WELLS:
                continue
            hit = [t for m, t in ANCHORED_WELLS[wk][1] if abs(m - md) < 1e-9]
            cells.append(f"{field_curve(wk, md):.2f}" if hit else "—")
        L.append(f"| {md} | " + " | ".join(cells) + " |")
    L.append("")
    (OUT / "七井温度场对账表.md").write_text("\n".join(L), encoding="utf-8")

    # ---------------- 汇总表 ----------------
    sum_rows = []
    for wk in WELLS:
        row: dict = {"井": wk}
        for suf, mode in GROUP_SUFFIX:
            r = res.get(f"{wk}_{suf}")
            row[f"{mode}_η_E"] = None if r is None else r["η_E"]
            row[f"{mode}_η_N"] = None if r is None else r["η_N"]
        rC, rB, rA = (res.get(f"{wk}_Tanchored_col"), res.get(f"{wk}_Tstatic_col"),
                      res.get(f"{wk}_Toff"))
        row["dEtaN_anchored_minus_static_pp"] = (
            None if not (rC and rB) else (rC["η_N"] - rB["η_N"]) * 100.0)
        row["dEtaN_static_minus_off_pp"] = (
            None if not (rB and rA) else (rB["η_N"] - rA["η_N"]) * 100.0)
        row["stop_t_s"] = None if rC is None else rC["配置"]["stop_t_s"]
        row["步数"] = None if rC is None else rC["health"]["步数"]
        row["col_memo"] = None if rC is None else rC["health"]["col_memo_size"]
        row["col_batches"] = None if rC is None else rC["health"]["col_batches"]
        row["uniform_hits"] = None if rC is None else rC["health"]["uniform_field_hits"]
        row["oob_1d_2d"] = (None if rC is None
                            else f"{rC['health']['oob_1d']}/{rC['health']['oob_2d']}")
        row["oob_col_2d"] = None if rC is None else rC["health"]["oob_col_2d"]
        row["cfl_clip"] = None if rC is None else rC["health"]["cfl_clip_events"]
        row["dt_median_s"] = None if rC is None else rC["health"]["dt_median_s"]
        sum_rows.append(row)

    with (OUT / "汇总表_七井温度场.csv").open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sum_rows[0]))
        w.writeheader()
        w.writerows(sum_rows)

    S: list[str] = []
    S.append("# 汇总表 · 七井温度场（Phase 4d 主批，2026-10-07）\n")
    S.append("> 底座：(F,F) 角 + `pressure_mode=\"off\"` + `nz=250` + CFL 自适应；r1.0。"
             "组 B/C 开 2D 逐列（探针 ctor 直装）；组 D（对桥）关。"
             "**全部后4a 口径**（4a 把 1D 温度消费从域顶单点改逐点现查）。\n")
    S.append("| 井 | η_E(off) | η_N(off) | η_E(static) | η_N(static) | η_E(anchored) | "
             "η_N(anchored) | Δη_N(锚−静)pp | Δη_N(静−关)pp | stop_t s | 步数 |")
    S.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in sum_rows:
        S.append("| " + " | ".join([
            WELL_CN[r["井"]], fmt(r["off_η_E"]), fmt(r["off_η_N"]),
            fmt(r["static_η_E"]), fmt(r["static_η_N"]),
            fmt(r["anchored_η_E"]), fmt(r["anchored_η_N"]),
            fmt_signed(r["dEtaN_anchored_minus_static_pp"]),
            fmt_signed(r["dEtaN_static_minus_off_pp"]),
            fmt(r["stop_t_s"], 1), "—" if r["步数"] is None else str(r["步数"]),
        ]) + " |")
    S.append("")
    S.append("## 关5 数值健康度（anchored 组）\n")
    S.append("| 井 | col_memo | col_batches | uniform_hits | oob(1D/2D) | oob_col(2D) | "
             "cfl_clip | dt中位 s |")
    S.append("|---|---|---|---|---|---|---|---|")
    for r in sum_rows:
        S.append("| " + " | ".join([
            WELL_CN[r["井"]],
            "—" if r["col_memo"] is None else str(r["col_memo"]),
            "—" if r["col_batches"] is None else str(r["col_batches"]),
            "—" if r["uniform_hits"] is None else str(r["uniform_hits"]),
            "—" if r["oob_1d_2d"] is None else str(r["oob_1d_2d"]),
            "—" if r["oob_col_2d"] is None else str(r["oob_col_2d"]),
            "—" if r["cfl_clip"] is None else str(r["cfl_clip"]),
            fmt(r["dt_median_s"], 3),
        ]) + " |")
    S.append("")
    (OUT / "汇总表_七井温度场.md").write_text("\n".join(S), encoding="utf-8")

    # ---------------- 一页解读 ----------------
    lines: list[str] = []
    lines.append("# Phase 4d 主批 · 一页解读（七井温度场，2026-10-07）\n")
    lines.append("## 口径与标签\n")
    lines.append("- **后4a 口径**（本批全部）：4a 已把 1D 温度消费从「域顶单点」改「逐点现查」；"
                 "引用 2e 四角 / 10-02 排量探针 / T2 批必须标「**前4a 口径**」。")
    lines.append("- **底座**：(F,F) 角 + pressure off + nz250 + CFL 自适应 + r1.0；"
                 "组 B/C 开 2D 逐列，组 D（对桥）关 ⇒ **跨组相减不可解释，组内比较有效**。")
    lines.append("- **无锚井**：呼探1 的温度行只有 3 条且同落 md=7601 m（1 个唯一深度，"
                 "`AnchoredProfileField` 需 ≥2 点）⇒ 走 **Geothermal 回退**（口径声明，不硬造锚）；"
                 "该井组 B/C 必然逐位重合，**不得读成「锚定无效」**。\n")
    lines.append("## 三条结论\n")
    dC = [r for r in sum_rows if r["dEtaN_anchored_minus_static_pp"] is not None]
    if dC:
        vals = [r["dEtaN_anchored_minus_static_pp"] for r in dC]
        worst = max(dC, key=lambda r: abs(r["dEtaN_anchored_minus_static_pp"]))
        lines.append(
            f"1. **锚定场 vs 统一地温式（组 C − 组 B）**：{len(dC)} 井可得，"
            f"Δη_N ∈ [{min(vals):+.3f}, {max(vals):+.3f}] pp；"
            f"绝对值最大者 = {WELL_CN[worst['井']]} "
            f"{worst['dEtaN_anchored_minus_static_pp']:+.3f} pp。"
            "Δ 由「锚定场在该井域内相对统一地温式的整体冷/热偏移」驱动。")
    dB = [r for r in sum_rows if r["dEtaN_static_minus_off_pp"] is not None]
    if dB:
        vals = [r["dEtaN_static_minus_off_pp"] for r in dB]
        lines.append(
            f"2. **温度档效应（组 B − 组 A）**：Δη_N ∈ [{min(vals):+.3f}, {max(vals):+.3f}] pp"
            "——温度流变口径替换与温度场的净效应之和（本设计内两者不可分离，"
            "术语见 CONTEXT.md「温度场效应与口径差」）。")
    rows_anch = [r for r in sum_rows if r["anchored_η_N"] is not None]
    if rows_anch:
        hi = max(rows_anch, key=lambda r: r["anchored_η_N"])
        lo = min(rows_anch, key=lambda r: r["anchored_η_N"])
        lines.append(
            f"3. **L3（跨井流变差异传导）**：同一公式分派 + 逐井不同锚定场 ⇒ η_N 跨井分化 "
            f"[{lo['anchored_η_N']:.4f}（{WELL_CN[lo['井']]}）→ "
            f"{hi['anchored_η_N']:.4f}（{WELL_CN[hi['井']]}）]，跨度 "
            f"{(hi['anchored_η_N'] - lo['anchored_η_N']) * 100:.2f} pp；"
            "锚点集与 k 逐井不同（spec §1）是输入侧的直接证据。")
    lines.append("")
    lines.append("## 对桥与表档\n")
    torch = res.get("呼101_Tstatic_nocol_ANCHOR")
    if torch:
        a = torch.get("anchor", {})
        lines.append(f"- **组 D hard 锚**：呼101 × static × 逐列关 ⇒ Δη_E={a.get('Δη_E')}、"
                     f"Δη_N={a.get('Δη_N')}（{'PASS' if a.get('PASS') else 'FAIL'}）"
                     "——与「后4a 口径排量锚」r1.0 点逐位一致。")
    tt = res.get("呼1-004_Ttablext_r1.0_col")
    tb = res.get("呼1-004_Ttable_r1.0_col")
    if tt and tb:
        same = (tt["η_E"] == tb["η_E"] and tt["η_N"] == tb["η_N"])
        lines.append(f"- **表档对桥**：呼1-004 扩展表 r1.0 η_E={tt['η_E']!r} / "
                     f"η_N={tt['η_N']!r}；交付表 r1.0 η_E={tb['η_E']!r} / "
                     f"η_N={tb['η_N']!r} ⇒ "
                     f"{'**逐位一致**（扩展表接线正确性实证）' if same else '**不一致，须排查**'}。")
    for tag in ("呼1-004_Ttablext_r0.6_col", "呼1-004_Ttablext_r0.8_col"):
        r = res.get(tag)
        if r:
            lines.append(f"- **{tag}**：η_E={r['η_E']:.10f}、η_N={r['η_N']:.10f}；"
                         f"coverage={r.get('coverage')}；oob="
                         f"{r['health']['oob_1d']}/{r['health']['oob_2d']}（真覆盖）。")
    lines.append("")
    lines.append("> **必附声明**：表档 r0.6/0.8 解禁 **≠ 洗脱名义排量热史错配**"
                 "（扩展段按末档 0.7 m³/min 外推，r0.6 实际为 0.6× 排量持续循环）；"
                 "**r0.5 维持禁引**（stop_t≈22819 s > 扩展表末 21600 s）。")
    lines.append("")
    (OUT / "一页解读.md").write_text("\n".join(lines), encoding="utf-8")

    print(f"[collector] 对账表 {len(csv_rows)} 锚点行 / 汇总表 {len(sum_rows)} 井 / "
          f"一页解读 已写")
    if missing:
        print(f"[collector] 仍缺 {len(missing)} 个 run")


if __name__ == "__main__":
    main()
