# -*- coding: utf-8 -*-
"""位移台账汇总（内部自洽加固 Task 7）。

把三个来源归并成一张汇总台账，写到**新文件**
``results/内部自洽加固_2026-09-25/位移台账_汇总.csv``：

1. **基线行** —— ``docs/superpowers/plans/baseline-2026-09-25.md``（Task 1 冻结）的 8 口井。
   基线数字只有 6 位小数，且产自**另一提交/解释器**（R8/R13）⇒ 只能作历史参照，
   严禁当作任何修正的"修正前"；因此基线行**自比**（修正前=修正后=基线值，Δ=0.00），
   说明列注明"不参与位移归因"。
2. **修正行** —— 原样并入 ``位移台账.csv`` 已有的行（每行是同窗 A/B，R3）。
   例外：``Δ*_pp`` 不复制原行已舍入的单元格，而由全精度 ``修正前/修正后`` **重算**（R143）。
3. **未测行** —— ht1_001 / hu2 的管容位移：用户 2026-09-25 裁定暂缓（R60），
   故按 R136 如实标"未测"，**不编造/外推数值**。

设计约束（附加裁定）：

* 原始 ``位移台账.csv`` 在本脚本里**只读**，逐字节不改（R137）——脚本幂等。
* 基线 md 表格**解析失败必须抛错**（``ValueError``），不得静默产出空表或漏行（R14）。
* 井名一律经 :func:`normalize_well_name` 规范化回内部代号（R9/R12）；
  映射表建在本脚本内，不改 Task 1 产物。
* ``Δ*_pp`` 一律 **4 位小数**（R143，覆盖计划原文的"2 位小数"）；
  并对原始 Δ 单元格做**显示舍入范围内**的一致性校验，超出即报错退出（不静默取其一）。

术语（仓库根 ``CONTEXT.md`` 强制）：本脚本只登记台账，不描述物理量；报告里描述环空
水泥尖端时**不得**称"前缘"（前缘/尾缘成对，且定义在**管内段**）。

运行：
    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/entrypoints/collect_displacement_ledger.py
"""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASELINE_MD = ROOT / "docs" / "superpowers" / "plans" / "baseline-2026-09-25.md"
RAW_LEDGER = ROOT / "results" / "内部自洽加固_2026-09-25" / "位移台账.csv"
SUMMARY_CSV = ROOT / "results" / "内部自洽加固_2026-09-25" / "位移台账_汇总.csv"

#: 汇总台账列名（逐字固定；test_displacement_ledger.py::test_ledger_columns 断言本表）。
COLUMNS = ["井名", "修正项", "修正前eta_E", "修正后eta_E", "Δeta_E_pp",
           "修正前eta_N", "修正后eta_N", "Δeta_N_pp", "说明"]

#: Δ 列显示小数位（R143：4 位，覆盖计划原文的 2 位）。
DELTA_DECIMALS = 4
#: 基线 md 的小数位（Task 1 冻结口径，R8/R13）。
BASELINE_DECIMALS = 6
#: 基线行"修正项"文案（R137）。
BASELINE_MODIFICATION = "基线（09-14 源模型口径重跑）"
#: "未测"占位（R136：如实标注，不编造数值）。
UNMEASURED = "未测"

#: 基线 md 应含的 8 口井（内部代号，R14：缺井/多井都必须报错）。
BASELINE_WELLS = ("hu101", "hu102", "hu103", "hu1", "hu2",
                  "ht1_001", "ht1_003", "ht1_004")

#: 井名规范化映射（R9/R12）。键含内部代号本身 ⇒ 幂等；
#: ⚠️ "呼探1" 在本库有两义（呼探1-001/-002 前缀，以及 hu1 的显示名），
#: 故 "呼探1-002"/"呼探1-001" 单列精确键，并由最长键优先兜底（见 normalize_well_name）。
WELL_NAME_MAP: dict[str, str] = {
    # 内部代号（幂等直通）
    "hu101": "hu101", "hu102": "hu102", "hu103": "hu103",
    "hu1": "hu1", "hu2": "hu2",
    "ht1_001": "ht1_001", "ht1_003": "ht1_003", "ht1_004": "ht1_004",
    # 现场/图件中文名（runner 输出命名口径）
    "呼101": "hu101", "呼102": "hu102", "呼103": "hu103",
    "呼探1": "hu1",            # hu1_tailpipe 输出 "呼探1尾管_…"
    "呼探1-002": "hu2",        # hu2_tailpipe 输出 "呼探1-002尾管_…"
    "呼探1-001": "ht1_001",    # ht1_001_tailpipe 输出 "呼探1-001尾管_…"
    "呼1-003": "ht1_003",      # ht1_003_tailpipe 输出 "呼1-003_…"
    "呼1-004": "ht1_004",      # ht1_004_tailpipe 输出 "呼1-004_…"
}
# 最长键优先：保证 "呼探1-002" 命中精确键，而不是被 "呼探1" 前缀吃掉。
_MAP_KEYS_LONGEST_FIRST = tuple(sorted(WELL_NAME_MAP, key=len, reverse=True))

#: 未测行（R136）。说明文案按 R136 逐字要求写"未测：用户 2026-09-25 裁定暂缓（R60），
#: 前提=常量口径；常量已由 R67/R68 裁定、Task 12 已落地收口"。
DEFERRED_ROWS = (
    ("ht1_001", "管容修正（地面→鞋口口径，R67）"),
    ("hu2", "管容修正（三段偏离纠正，R68）"),
)
_DEFERRED_NOTE = ("未测：用户 2026-09-25 裁定暂缓（R60），前提=常量口径；"
                  "常量已由 R67/R68 裁定、Task 12 已落地收口。"
                  "补跑成本估计见 .superpowers/sdd/2026-09-25-internal-consistency-hardening/"
                  "task-7-report.md")


def normalize_well_name(name: str) -> str:
    """把井名规范化回内部代号（R9/R12）。

    先精确匹配映射表（含内部代号本身 ⇒ 幂等）；否则按**最长键优先**做前缀匹配。
    无法识别时抛 ``ValueError``——绝不静默放行未知井名（R14）。
    """
    key = str(name).strip()
    if not key:
        raise ValueError("井名为空，无法规范化")
    if key in WELL_NAME_MAP:
        return WELL_NAME_MAP[key]
    for cand in _MAP_KEYS_LONGEST_FIRST:
        if key.startswith(cand):
            return WELL_NAME_MAP[cand]
    raise ValueError(f"未知井名 {name!r}：映射表无匹配（不静默放行）")


def parse_baseline_table(text: str, *,
                         expected_wells: tuple[str, ...] | None = BASELINE_WELLS) -> list[dict]:
    """从基线 md 解析 η 表，返回 ``[{"well": 内部代号, "eta_E": float, "eta_N": float}, …]``。

    **解析失败一律抛 ``ValueError``**（R14）：无表头、无数据行、列数不对、η 非数值、
    井名重复、井集与 ``expected_wells`` 不符——任一情形都不静默降级为短表/空表。
    """
    rows: list[dict] = []
    header_seen = False
    for lineno, raw_line in enumerate(text.splitlines(), 1):
        line = raw_line.strip()
        if not line.startswith("|"):
            if header_seen:
                break  # 表格结束（表头之后遇到首个非表格行即止）
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if not header_seen:
            if cells == ["井", "η_E", "η_N"]:
                header_seen = True
            continue
        if all(c and set(c) <= {"-", ":"} for c in cells):
            continue  # 分隔行 |---|---|…|
        if len(cells) != 3:
            raise ValueError(f"基线 md 第 {lineno} 行列数不为 3：{line!r}")
        well_raw, eta_e_s, eta_n_s = cells
        try:
            eta_e = float(eta_e_s)
            eta_n = float(eta_n_s)
        except ValueError as exc:
            raise ValueError(f"基线 md 第 {lineno} 行 η 不是数值：{line!r}") from exc
        rows.append({"well": normalize_well_name(well_raw), "eta_E": eta_e, "eta_N": eta_n})
    if not header_seen:
        raise ValueError("基线 md 未找到 `| 井 | η_E | η_N |` 表头 ⇒ 解析失败，拒绝产出空表（R14）")
    if not rows:
        raise ValueError("基线 md 有表头但无数据行 ⇒ 解析失败，拒绝产出空表（R14）")
    got = [r["well"] for r in rows]
    if len(set(got)) != len(got):
        raise ValueError(f"基线 md 井名重复：{got}")
    if expected_wells is not None and set(got) != set(expected_wells):
        raise ValueError(
            f"基线 md 井集与预期不符（缺井/多井）⇒ 解析失败："
            f"实测 {sorted(got)} != 预期 {sorted(expected_wells)}")
    return rows


def read_raw_ledger(path: Path = RAW_LEDGER) -> list[dict[str, str]]:
    """读原始位移台账（**只读**，绝不写回）。表头必须逐字等于 ``COLUMNS``（R137）。"""
    with path.open(encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        header = reader.fieldnames
        if header != COLUMNS:
            raise ValueError(f"原始位移台账表头与 COLUMNS 不符：{header!r}")
        return [dict(row) for row in reader]


def recompute_delta_pp(before: float, after: float) -> float:
    """``Δ*_pp = (修正后 − 修正前) × 100``（全精度；R143）。"""
    return (after - before) * 100.0


def format_delta_pp(value: float) -> str:
    """Δ 列显示格式：4 位小数（R143）。"""
    return f"{value:.{DELTA_DECIMALS}f}"


def _display_decimals(cell: str) -> int:
    """原始 Δ 单元格的显示小数位（用于推"显示舍入"容差）。"""
    return len(cell.split(".", 1)[1]) if "." in cell else 0


def check_delta_consistency(*, well: str, item: str, axis: str,
                            stored_cell: str, recomputed_pp: float) -> str | None:
    """校验原始 Δ 单元格与重算值一致（R143）。

    容差 = 该单元格**显示分辨率的一半**（末位舍入上界）。超出 ⇒ 抛 ``ValueError``
    （解析/口径错，绝不静默取其一）；在容差内但与重算值不等 ⇒ 返回一条说明供上层暴露
    （"原始单元格是舍入显示值，汇总按全精度重算值写入"）；未测/空单元格返回 ``None``。
    """
    cell = (stored_cell or "").strip()
    if cell in ("", UNMEASURED, "nan", "NaN"):
        return None
    try:
        stored = float(cell)
    except ValueError as exc:
        raise ValueError(f"[{well}/{item}] {axis} 原始 Δ 单元格不可解析：{cell!r}") from exc
    tol = 0.5 * 10.0 ** (-_display_decimals(cell)) + 1e-9
    dev = stored - recomputed_pp
    if abs(dev) > tol:
        raise ValueError(
            f"[{well}/{item}] {axis} 原始 Δ={cell!r} 与由 修正前/修正后 重算的 "
            f"{recomputed_pp:.6f} pp 之差 {dev:+.6g} 超出显示舍入容差 {tol:g} "
            f"⇒ 解析/口径错，报错退出（不静默取其一，R143）")
    if stored == recomputed_pp:
        return None
    return (f"[{well}/{item}] {axis} 原始 Δ 单元格 {cell!r}"
            f"（{_display_decimals(cell)} 位显示舍入）≠ 重算值 {format_delta_pp(recomputed_pp)}"
            f"（全精度 {recomputed_pp:+.6g} pp）；汇总按重算值写入（R143）")


def build_summary_rows(baseline_rows: list[dict],
                       raw_rows: list[dict]) -> tuple[list[dict], list[str]]:
    """归并三条来源 ⇒ ``(rows, notes)``。

    顺序：基线行 → 原样并入的修正行 → 未测行（R137）。``notes`` 记录原始 Δ 单元格与
    重算值在显示舍入内的差异（暴露"复制原始 Δ 会丢信号"，R143）。
    """
    rows: list[dict] = []
    notes: list[str] = []

    # ① 基线行：自比 Δ=0.00；6 位小数只作历史参照（R8/R13）。
    for b in baseline_rows:
        eta_e = f"{b['eta_E']:.{BASELINE_DECIMALS}f}"
        eta_n = f"{b['eta_N']:.{BASELINE_DECIMALS}f}"
        rows.append({
            "井名": b["well"],
            "修正项": BASELINE_MODIFICATION,
            "修正前eta_E": eta_e, "修正后eta_E": eta_e,
            "Δeta_E_pp": f"{0.0:.{DELTA_DECIMALS}f}",
            "修正前eta_N": eta_n, "修正后eta_N": eta_n,
            "Δeta_N_pp": f"{0.0:.{DELTA_DECIMALS}f}",
            "说明": "基线行不参与位移归因：修正前=修正后=基线值（自比），Δ 恒为 0；"
                    "基线 6 位小数产自另一提交/解释器（R8/R13），只作历史参照，"
                    "严禁当作任何修正的“修正前”",
        })

    # ② 修正行：井名/修正项/修正前/修正后/说明 原样并入；Δ 由全精度两列重算（R143）。
    for raw_row in raw_rows:
        out = {col: (raw_row.get(col) or "").strip() for col in COLUMNS}
        well_cn = out["井名"]
        try:
            well = normalize_well_name(well_cn)
        except ValueError as exc:
            raise ValueError(f"原始位移台账行井名无法规范化：{well_cn!r}") from exc
        item = out["修正项"]
        for axis in ("eta_E", "eta_N"):
            before_s, after_s = out[f"修正前{axis}"], out[f"修正后{axis}"]
            try:
                before, after = float(before_s), float(after_s)
            except ValueError as exc:
                raise ValueError(
                    f"[{well}/{item}] {axis} 修正前/修正后 不可解析"
                    f"（{before_s!r}/{after_s!r}）⇒ 无法重算 Δ（R143）") from exc
            recomputed = recompute_delta_pp(before, after)
            note = check_delta_consistency(well=well, item=item, axis=axis,
                                           stored_cell=out[f"Δ{axis}_pp"],
                                           recomputed_pp=recomputed)
            if note:
                notes.append(note)
            out[f"Δ{axis}_pp"] = format_delta_pp(recomputed)
        rows.append(out)

    # ③ 未测行：R136 —— 如实标注，不编造/外推（R60 用户裁定暂缓）。
    for well, item in DEFERRED_ROWS:
        rows.append({
            "井名": well, "修正项": item,
            "修正前eta_E": UNMEASURED, "修正后eta_E": UNMEASURED, "Δeta_E_pp": UNMEASURED,
            "修正前eta_N": UNMEASURED, "修正后eta_N": UNMEASURED, "Δeta_N_pp": UNMEASURED,
            "说明": _DEFERRED_NOTE,
        })

    return rows, notes


def collect(*, baseline_md: Path = BASELINE_MD,
            raw_ledger: Path = RAW_LEDGER) -> tuple[list[dict], list[str]]:
    """读入基线 md + 原始台账 ⇒ ``(汇总行, 暴露说明)``。只读，不写任何文件。"""
    baseline_rows = parse_baseline_table(baseline_md.read_text(encoding="utf-8-sig"))
    raw_rows = read_raw_ledger(raw_ledger)
    return build_summary_rows(baseline_rows, raw_rows)


def write_summary(rows: list[dict], path: Path = SUMMARY_CSV) -> Path:
    """写汇总台账（新文件，R137）。原始 ``位移台账.csv`` 不在任何写路径上。

    编码与行尾沿用同族原始台账（UTF-8 无 BOM + CRLF），使表头**逐字节**等于
    ``COLUMNS`` 的 UTF-8 编码、且与原始台账并排可比。
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def main() -> int:
    rows, notes = collect()
    n_baseline = sum(1 for r in rows if r["修正项"] == BASELINE_MODIFICATION)
    n_deferred = sum(1 for r in rows if r["Δeta_E_pp"] == UNMEASURED)
    n_fix = len(rows) - n_baseline - n_deferred
    path = write_summary(rows)
    print(f"[ledger] 汇总 {len(rows)} 行 → {path}")
    print(f"         （基线 {n_baseline} 行 + 修正 {n_fix} 行 + 未测 {n_deferred} 行）")
    for note in notes:
        print(f"  ⚠ {note}")
    print(f"[ledger] 原始 {RAW_LEDGER.name} 只读，未改动（R137）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
