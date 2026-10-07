#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Phase 4d 主批作业卡生成器（确定性；输出 `jobs_seven_wells_20261007.json`）

矩阵见 spec `docs/superpowers/specs/2026-10-07-phase4d-seven-well-batch-design.md` §3/§3.1：

| 组 | 井 | 温度档 | 2D 逐列 | 说明 |
|---|---|---|---|---|
| A | 七井 | off | — | T-off 对照 |
| B | 七井 | static | 开 | 无锚场对照 |
| C | 七井 | anchored | 开 | 主产物（锚定场逐井不同） |
| D | 呼101 | static / anchored | 关 | static 支带 hard 锚（后4a 基线逐位） |
| E | 呼1-004 | table_ext / table | 开 | 扩展表 r0.6/0.8/1.0 + 交付表 r1.0 对桥 |

统一底座：`(F,F)` 角 + `pressure_mode="off"` + `pressure_caliber="shoe"` + `nz=250`。

运行：PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -I make_jobs_20261007.py
"""

from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

WELLS = ("呼101", "呼1-003", "呼1-004", "呼102", "呼探1", "呼探1-002", "呼探1-001")

BASE = {
    "split": False,          # (F,·) 角：include_yield_term=False
    "gate": False,           # (·,F) 角：enable_stream_yield_gate=False
    "nz": 250,
    "rate": 1.0,
    "pressure_mode": "off",
    "pressure_caliber": "shoe",
}


def job(tag: str, well: str, mode: str, *, depthwise: bool, rate: float = 1.0,
        anchor: dict | None = None) -> dict:
    return {"tag": tag, "well": well, "mode": mode, "depthwise": depthwise,
            "anchor": anchor, **{**BASE, "rate": rate}}


def build() -> list[dict]:
    jobs: list[dict] = []
    # 组 A：T-off 对照（七井）
    for wk in WELLS:
        jobs.append(job(f"{wk}_Toff", wk, "off", depthwise=False))
    # 组 B：static（统一地温式）+ 2D 逐列
    for wk in WELLS:
        jobs.append(job(f"{wk}_Tstatic_col", wk, "static", depthwise=True))
    # 组 C：anchored（七井锚定场；呼探1 ⇒ Geothermal 回退）+ 2D 逐列
    for wk in WELLS:
        jobs.append(job(f"{wk}_Tanchored_col", wk, "anchored", depthwise=True))
    # 组 D：呼101 与后4a 基线对桥（逐列关 ⇒ 与基线同口径）
    jobs.append(job("呼101_Tstatic_nocol_ANCHOR", "呼101", "static", depthwise=False,
                    anchor={"kind": "file", "mode": "hard"}))
    jobs.append(job("呼101_Tanchored_nocol", "呼101", "anchored", depthwise=False))
    # 组 E：呼1-004 表档（扩展表 r0.6/0.8/1.0 + 交付表 r1.0 对桥）
    for r in (0.6, 0.8, 1.0):
        jobs.append(job(f"呼1-004_Ttablext_r{r}_col", "呼1-004", "table_ext",
                        depthwise=True, rate=r))
    jobs.append(job("呼1-004_Ttable_r1.0_col", "呼1-004", "table",
                    depthwise=True, rate=1.0))
    return jobs


def main() -> None:
    jobs = build()
    tags = [j["tag"] for j in jobs]
    assert len(tags) == len(set(tags)), "tag 重复"
    out = HERE / "jobs_seven_wells_20261007.json"
    out.write_text(json.dumps(jobs, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"[make_jobs] {len(jobs)} 作业 -> {out.name}")
    for j in jobs:
        print(f"  {j['tag']:<32} well={j['well']:<10} mode={j['mode']:<10} "
              f"rate={j['rate']} depthwise={j['depthwise']} anchor={bool(j['anchor'])}")


if __name__ == "__main__":
    main()
