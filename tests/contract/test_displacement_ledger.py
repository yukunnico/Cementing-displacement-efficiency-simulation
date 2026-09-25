# -*- coding: utf-8 -*-
"""位移台账汇总（内部自洽加固 Task 7）契约测试。

守护四件事：
1. ``COLUMNS`` 逐字固定（Step 1 的原测试）；
2. 归并行数 = 基线 8 行 + 原始台账各行 + 未测 2 行（R136/R137），井名经规范化；
3. 基线 md 解析失败**必须报错**，不得静默产出空表/短表（R14）；
4. ``Δ*_pp`` 由全精度 修正前/修正后 重算、4 位小数（R143）——2 位小数会抹掉
   Task 8 那行的全部信息（真值 −0.000602/+0.001525 pp 渲染成 −0.00/0.00）；
   超出显示舍入的不一致必须报错，容差内的不一致必须在 notes 里暴露；
5. 原始 ``位移台账.csv`` 跑一次后**逐字节不变**（R137，脚本只读）。
"""
import importlib

import pytest

# 与 results/内部自洽加固_2026-09-25/位移台账.csv 的 Task 8 行逐字一致
# （η_E 真值 −0.000602 pp、η_N 真值 +0.001525 pp；原始 Δ 单元格是 2 位舍入的 -0.00/0.00）。
TASK8_RAW_ROW = {
    "井名": "hu101(锚)",
    "修正项": "线性求解换带状Cholesky",
    "修正前eta_E": "0.9203933244336414",
    "修正后eta_E": "0.9203873022217786",
    "Δeta_E_pp": "-0.00",
    "修正前eta_N": "0.7188836497828229",
    "修正后eta_N": "0.7188988965864732",
    "Δeta_N_pp": "0.00",
    "说明": "锚算例=hu101；本行仅用于测试。",
}

SYNTHETIC_BASELINE = """# 基线冻结（测试夹具）

| 井 | η_E | η_N |
|---|---|---|
| 呼101 | 0.971075 | 0.829436 |
| 呼102 | 0.985278 | 0.886027 |
| 呼103 | 0.997243 | 0.987984 |
| 呼探1 | 0.870732 | 0.790153 |
| 呼探1-002 | 0.986364 | 0.991253 |
| 呼探1-001 | 0.947220 | 0.933968 |
| 呼1-003 | 0.998120 | 0.992673 |
| 呼1-004 | 0.998741 | 0.994083 |
"""

EXPECTED_BASELINE_ORDER = ["hu101", "hu102", "hu103", "hu1", "hu2",
                           "ht1_001", "ht1_003", "ht1_004"]


def _mod():
    return importlib.import_module("scripts.entrypoints.collect_displacement_ledger")


def test_ledger_columns():
    mod = _mod()
    assert mod.COLUMNS == ["井名", "修正项", "修正前eta_E", "修正后eta_E", "Δeta_E_pp",
                           "修正前eta_N", "修正后eta_N", "Δeta_N_pp", "说明"]


def test_aggregation_row_counts_and_well_normalization():
    """基线 8 行 + 原始 1 行 + 未测 2 行 = 11 行；中文井名规范化回内部代号。"""
    mod = _mod()
    baseline_rows = mod.parse_baseline_table(SYNTHETIC_BASELINE)
    assert [r["well"] for r in baseline_rows] == EXPECTED_BASELINE_ORDER
    rows, _ = mod.build_summary_rows(baseline_rows, [dict(TASK8_RAW_ROW)])
    assert len(rows) == 8 + 1 + 2

    names = [r["井名"] for r in rows]
    # ① 基线 8 行用内部代号，顺序沿用基线 md
    assert names[:8] == EXPECTED_BASELINE_ORDER
    # ② 原始行原样并入（井名保持 "hu101(锚)"，不重写）
    assert names[8] == "hu101(锚)"
    # ③ 未测 2 行
    assert names[9:] == ["ht1_001", "hu2"]
    # 基线行自比 Δ=0（R143 裁定 Δ 列一律 4 位小数 ⇒ 渲染为 "0.0000"，非计划原文的 "0.00"）
    for r in rows[:8]:
        assert r["修正项"] == mod.BASELINE_MODIFICATION
        assert r["Δeta_E_pp"] == "0.0000" and r["Δeta_N_pp"] == "0.0000"
        assert r["修正前eta_E"] == r["修正后eta_E"]
        assert "不参与位移归因" in r["说明"]
    # 未测行如实标注、无编造数值（R136）
    for r in rows[9:]:
        assert r["修正前eta_E"] == r["修正后eta_E"] == r["Δeta_E_pp"] == mod.UNMEASURED
        assert "R60" in r["说明"] and "R67/R68" in r["说明"] and "Task 12" in r["说明"]


def test_normalize_well_name_idempotent_and_longest_prefix():
    """规范化幂等；'呼探1-002' 必须命中精确键，不被 '呼探1' 前缀吃掉（两义陷阱）。"""
    mod = _mod()
    for internal in EXPECTED_BASELINE_ORDER:
        assert mod.normalize_well_name(internal) == internal
    assert mod.normalize_well_name("呼101") == "hu101"
    assert mod.normalize_well_name("呼探1") == "hu1"
    assert mod.normalize_well_name("呼探1-002") == "hu2"
    assert mod.normalize_well_name("呼探1-001") == "ht1_001"
    assert mod.normalize_well_name("hu101(锚)") == "hu101"  # 前缀兜底（仅用于识别）
    # 前缀兜底也必须**最长优先**：带后缀的复合名不得被更短的 '呼探1' 吃掉
    # （否则 "呼探1-002(管容)" 会误判成 hu1、"呼探1-001(锚)" 同）
    assert mod.normalize_well_name("呼探1-002(管容)") == "hu2"
    assert mod.normalize_well_name("呼探1-001(锚)") == "ht1_001"
    assert mod.normalize_well_name("呼1-004_尾管") == "ht1_004"
    with pytest.raises(ValueError):
        mod.normalize_well_name("呼999")
    with pytest.raises(ValueError):
        mod.normalize_well_name("   ")


@pytest.mark.parametrize("bad_text, why", [
    ("# 没有表格\n\n纯文字，没有任何 `| 井 | …` 表。\n", "无表头"),
    ("| 井 | η_E | η_N |\n|---|---|---|\n", "有表头无数据行"),
    ("| 井 | η_E | η_N |\n|---|---|---|\n| hu101 | abc | 0.1 |\n", "η 非数值"),
    ("| 井 | η_E | η_N |\n|---|---|---|\n| hu101 | 0.97 | 0.82 |\n", "缺 7 口井"),
    ("| 井 | η_E | η_N |\n|---|---|---|\n"
     "| hu101 | 0.97 | 0.82 |\n| hu101 | 0.97 | 0.82 |\n", "井名重复"),
])
def test_malformed_baseline_fails_loudly(bad_text, why):
    """解析失败必须抛错（R14）——不得静默产出空表/短表。"""
    mod = _mod()
    with pytest.raises(ValueError):
        mod.parse_baseline_table(bad_text)


def test_delta_recomputed_at_four_decimals_and_signal_preserved():
    """原始 Δ 单元格 2 位舍入为 0.00，但真值非零 ⇒ 汇总必须按重算值写、且暴露不一致。"""
    mod = _mod()
    baseline_rows = mod.parse_baseline_table(SYNTHETIC_BASELINE)
    rows, notes = mod.build_summary_rows(baseline_rows, [dict(TASK8_RAW_ROW)])
    fix = rows[8]
    # 原始单元格（若照抄会得到 -0.00 / 0.00，抹掉该列全部信号）
    assert TASK8_RAW_ROW["Δeta_E_pp"] == "-0.00" and TASK8_RAW_ROW["Δeta_N_pp"] == "0.00"
    # 汇总：4 位小数、由全精度两列重算 ⇒ 非零
    assert fix["Δeta_E_pp"] == "-0.0006"
    assert fix["Δeta_N_pp"] == "0.0015"
    # 一致性校验结果必须在 notes 里被暴露（R143）
    assert notes, "原始 Δ 单元格与重算值不一致必须在 notes 里暴露"
    assert any("hu101" in n and "eta_E" in n for n in notes)
    assert any("hu101" in n and "eta_N" in n for n in notes)


def test_consistent_stored_delta_yields_no_note():
    """原始 Δ 单元格与重算值逐位一致时不得产生噪音暴露（真一致，非仅显示舍入）。"""
    mod = _mod()
    raw = dict(TASK8_RAW_ROW)
    # 前后完全相同 ⇒ 真值恰为 0；原始单元格 "0.00" 与重算值逐位一致 ⇒ 无暴露
    raw["修正前eta_E"] = raw["修正后eta_E"] = "0.5000000000000000"
    raw["修正前eta_N"] = raw["修正后eta_N"] = "0.4000000000000000"
    raw["Δeta_E_pp"] = "0.00"
    raw["Δeta_N_pp"] = "0.00"
    rows, notes = mod.build_summary_rows(mod.parse_baseline_table(SYNTHETIC_BASELINE), [raw])
    assert rows[8]["Δeta_E_pp"] == "0.0000" and rows[8]["Δeta_N_pp"] == "0.0000"
    assert notes == []


def test_stored_delta_beyond_display_rounding_raises():
    """原始 Δ 单元格与重算值之差超出显示舍入 ⇒ 报错退出，不静默取其一（R143）。"""
    mod = _mod()
    raw = dict(TASK8_RAW_ROW)
    # 修正后 − 修正前 = +1 pp（远超 2 位显示舍入的 0.005），而原始单元格写 0.00
    raw["修正前eta_E"] = "0.5000000000000000"
    raw["修正后eta_E"] = "0.5100000000000000"
    raw["Δeta_E_pp"] = "0.00"
    with pytest.raises(ValueError):
        mod.build_summary_rows(mod.parse_baseline_table(SYNTHETIC_BASELINE), [raw])


def test_raw_ledger_bytes_unchanged_by_run(tmp_path):
    """跑一次完整归并（含真实基线 md + 真实原始台账）后，原始台账逐字节不变（R137）。"""
    mod = _mod()
    before = mod.RAW_LEDGER.read_bytes()
    rows, _ = mod.collect(baseline_md=mod.BASELINE_MD, raw_ledger=mod.RAW_LEDGER)
    out = mod.write_summary(rows, path=tmp_path / "位移台账_汇总.csv")
    assert out.exists()
    assert mod.RAW_LEDGER.read_bytes() == before, "原始位移台账被改写（违反 R137 只读约束）"
    # 真实产物：基线 8 + 修正（当前 1） + 未测 2
    assert len(rows) == 8 + 1 + 2
    # 表头逐字等于 COLUMNS（无 BOM，UTF-8）
    first_line = out.read_text(encoding="utf-8").splitlines()[0]
    assert first_line == ",".join(mod.COLUMNS)
