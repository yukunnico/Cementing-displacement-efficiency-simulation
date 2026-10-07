#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Phase 5e：T2 三重点井关键档「后4a 口径」总对账（新旧并列 + 口径升级声明）

输入：
- 旧（前4a 口径，2026-10-06 批）：`results/敏感性补跑_三重点井_20261006/` 的
  `{井}_{变体}_结果摘要.json`
- 新（后4a 口径，本窗口重跑）：`results/敏感性补跑_三重点井_后4a_2026-10-07/` 同名文件

输出（本目录）：
- `5e_T2总对账.md` / `5e_T2总对账.csv`（逐变体新旧 η_E/η_N + Δ + 口径标签）

口径声明（必须随表引用）：
- 旧 = **前4a 口径**（4a 之前：1D 温度消费取「域顶单点」）；新 = **后4a 口径**
  （4a 之后：1D 逐点现查）。**两套数字禁混引**，本表 Δ 即口径升级的净位移。

运行：PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -I collect_5e_20261007.py
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
# 稳健定位仓库根：向上找含 `cemdisp` 的目录（不依赖脚本层级深度）
def _find_root(start: Path) -> Path:
    for d in (start, *start.parents):
        if (d / "cemdisp").is_dir():
            return d
    raise RuntimeError(f"未找到含 cemdisp 的仓库根，起点 {start}")


ROOT = _find_root(HERE)
OLD = ROOT / "results" / "敏感性补跑_三重点井_20261006"
NEW = ROOT / "results" / "敏感性补跑_三重点井_后4a_2026-10-07"
WELLS = ("呼101", "呼1-003", "呼1-004")
_SUF = "_结果摘要.json"


def _final(path: Path) -> dict | None:
    if not path.exists():
        return None
    d = json.loads(path.read_text(encoding="utf-8"))
    f = d.get("最终结果", {})
    return {
        "η_E": float(f["全井段最终有效顶替效率"]),
        "η_N": float(f["窄四分位效率"]),
    }


def main() -> None:
    rows: list[dict] = []
    for w in WELLS:
        old_files = {p.name[len(w) + 1:-len(_SUF)]: p
                     for p in sorted(OLD.glob(f"{w}_*{_SUF}"))}
        new_files = {p.name[len(w) + 1:-len(_SUF)]: p
                     for p in sorted(NEW.glob(f"{w}_*{_SUF}"))}
        for var in sorted(set(old_files) | set(new_files)):
            o = _final(old_files[var]) if var in old_files else None
            n = _final(new_files[var]) if var in new_files else None
            rows.append({
                "井名": w, "变体": var,
                "η_E_前4a": None if o is None else o["η_E"],
                "η_N_前4a": None if o is None else o["η_N"],
                "η_E_后4a": None if n is None else n["η_E"],
                "η_N_后4a": None if n is None else n["η_N"],
                "Δη_E_pp": None if not (o and n) else (n["η_E"] - o["η_E"]) * 100.0,
                "Δη_N_pp": None if not (o and n) else (n["η_N"] - o["η_N"]) * 100.0,
                "状态": "双有" if (o and n) else ("仅前4a" if o else "仅后4a"),
            })

    if rows:
        with (HERE / "5e_T2总对账.csv").open("w", encoding="utf-8-sig", newline="") as fh:
            wtr = csv.DictWriter(fh, fieldnames=list(rows[0]))
            wtr.writeheader()
            wtr.writerows(rows)

    L: list[str] = []
    L.append("# Phase 5e · T2 三重点井关键档总对账（新旧并列）\n")
    L.append("## 口径声明（引用必附）\n")
    L.append("- **旧 = 前4a 口径**（4a 之前：1D 温度消费取「域顶单点」，全程恒定）；"
             "**新 = 后4a 口径**（4a 之后：1D 逐点现查 `T(md,t)`）。")
    L.append("- **两套数字禁混引**；本表 `Δ` 即**口径升级引入的净位移**，"
             "不是模型改进的效果量。旧批目录 `results/敏感性补跑_三重点井_20261006/` "
             "**未被覆盖**（红线：补跑新目录带日期后缀）。")
    L.append("- 底座（新批）：`CORRECTED_KW` + nz250 + CFL 自适应 + `pressure_mode=\"off\"`；"
             "**未开 2D 逐列**（`enable_depthwise_temperature` 默认 False）——"
             "即新批相对旧批的**唯一系统性差异 = 4a 的 1D 逐点现查**。\n")
    dual = [r for r in rows if r["状态"] == "双有"]
    L.append(f"## 汇总：{len(rows)} 行，其中双有（旧+新同变体）**{len(dual)}** 个\n")
    if dual:
        dE = sorted(r["Δη_E_pp"] for r in dual)
        dN = sorted(r["Δη_N_pp"] for r in dual)
        mid = lambda v: v[len(v) // 2]
        L.append(f"- **Δη_E（后4a − 前4a）**：min {dE[0]:+.4f} / median {mid(dE):+.4f} / "
                 f"max {dE[-1]:+.4f} pp")
        L.append(f"- **Δη_N**：min {dN[0]:+.4f} / median {mid(dN):+.4f} / "
                 f"max {dN[-1]:+.4f} pp")
        L.append("- 与既有结论「T-on 数字整体位移 −0.3 ~ −0.8 pp」对照：本表是该结论在"
                 "三重点井全变体上的完整展开。\n")
    L.append("## 逐变体明细\n")
    L.append("| 井 | 变体 | η_E(前4a) | η_N(前4a) | η_E(后4a) | η_N(后4a) | "
             "Δη_E pp | Δη_N pp | 状态 |")
    L.append("|---|---|---|---|---|---|---|---|---|")

    def f6(v):
        return "—" if v is None else f"{v:.6f}"

    def f3(v):
        return "—" if v is None else f"{v:+.3f}"

    for r in rows:
        L.append(f"| {r['井名']} | {r['变体']} | {f6(r['η_E_前4a'])} | {f6(r['η_N_前4a'])} | "
                 f"{f6(r['η_E_后4a'])} | {f6(r['η_N_后4a'])} | {f3(r['Δη_E_pp'])} | "
                 f"{f3(r['Δη_N_pp'])} | {r['状态']} |")
    L.append("")
    (HERE / "5e_T2总对账.md").write_text("\n".join(L), encoding="utf-8")
    print(f"[5e] {len(rows)} 行（双有 {len(dual)}）-> 5e_T2总对账.md/.csv")


if __name__ == "__main__":
    main()
