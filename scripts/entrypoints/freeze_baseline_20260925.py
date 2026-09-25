# -*- coding: utf-8 -*-
"""冻结 2026-09-25 基线：HEAD、环境、8 井当前口径数字。

只读脚本——不跑求解、不写任何 `results/` 目录。产物是
`docs/superpowers/plans/baseline-2026-09-25.md`，供"内部自洽加固"计划后续任务
比对位移（Task 7 的位移台账以它为基线行）。

数据来源：`results/源模型口径重跑_2026-09-14/汇总.csv`（源模型口径重跑的八井
汇总；表头 2026-09-25 核实为
`cement_occ,cfl_mode,channeling,corrected,elapsed_s,eta_E,eta_N,instability,mixing,wall_frac,well`，
UTF-8 BOM ⇒ 用 `utf-8-sig` 读取）。

运行：
    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/entrypoints/freeze_baseline_20260925.py
"""
from __future__ import annotations

import csv
import platform
import subprocess
import sys
from pathlib import Path

import numpy
import scipy

ROOT = Path(__file__).resolve().parents[2]
SUMMARY = ROOT / "results" / "源模型口径重跑_2026-09-14" / "汇总.csv"
OUT = ROOT / "docs" / "superpowers" / "plans" / "baseline-2026-09-25.md"


def _git(*args: str) -> str:
    """在当前仓库执行只读 git 命令并返回去尾空白 stdout。"""
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True,
                          check=True).stdout.strip()


def collect_baseline() -> dict:
    """采集基线：HEAD / 工作树是否脏 / 解释器与依赖版本 / 八井 η_E、η_N。

    井序沿用 `汇总.csv` 自带顺序，不做重排（保证与源表可比）。
    """
    wells = []
    with SUMMARY.open(encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            wells.append({"well": row["well"], "eta_E": float(row["eta_E"]),
                          "eta_N": float(row["eta_N"])})
    return {"head_commit": _git("rev-parse", "HEAD"),
            "dirty": bool(_git("status", "--porcelain")),
            "python": platform.python_version(),
            "numpy": numpy.__version__,
            "scipy": scipy.__version__,
            "wells": wells}


def main() -> int:
    b = collect_baseline()
    lines = ["# 基线冻结（2026-09-25，内部自洽加固前）", "",
             f"- HEAD: `{b['head_commit']}`（工作树{'有' if b['dirty'] else '无'}未提交改动）",
             f"- 环境: Python {b['python']} / numpy {b['numpy']} / scipy {b['scipy']}",
             f"- 口径来源: `results/源模型口径重跑_2026-09-14/汇总.csv`（源模型口径重跑）",
             "", "> 供后续任务比对位移；**加固期间不得修改**（Task 7 位移台账的基线行）。",
             "", "| 井 | η_E | η_N |", "|---|---|---|"]
    lines += [f"| {w['well']} | {w['eta_E']:.6f} | {w['eta_N']:.6f} |" for w in b["wells"]]
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[baseline] {OUT}  ({len(b['wells'])} 井)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
