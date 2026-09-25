# -*- coding: utf-8 -*-
"""2026-09-25 基线冻结（内部自洽加固 Task 1）。

本文件守护"改动前基线可复算"这一件事：`collect_baseline()` 必须能读出现存
口径表并给出 HEAD、解释器版本与八井 η 三元组。基线数字本身冻结在
`docs/superpowers/plans/baseline-2026-09-25.md`，供后续任务比对位移。

实际表头（2026-09-25 核实，UTF-8 BOM）：
`cement_occ,cfl_mode,channeling,corrected,elapsed_s,eta_E,eta_N,instability,mixing,wall_frac,well`
——`well`/`eta_E`/`eta_N` 三列名与计划假设一致，无需改名（见 Ruling R4）。
"""
import importlib


def test_collect_baseline_shape():
    mod = importlib.import_module("scripts.entrypoints.freeze_baseline_20260925")
    b = mod.collect_baseline()
    assert b["head_commit"] and len(b["head_commit"]) >= 7
    assert b["python"].startswith("3.13")
    assert b["wells"] and all(set(w) >= {"well", "eta_E", "eta_N"} for w in b["wells"])
