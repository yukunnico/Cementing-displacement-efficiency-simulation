# -*- coding: utf-8 -*-
"""物理合理性台账：对可跑井跑 7 项物理判据（Task 9）。

产出
----
``results/内部自洽加固_2026-09-25/物理合理性台账.csv``，列为：
``井名, 检查项, 通过, 实测值, 判据, 说明``（``通过`` ∈ {通过, 不通过, 未测}）。

**行结构 = 8 井 × 4 项结果对象级判据（P-1/P-2/P-5/P-7）+ 3 行合成算例判据
（P-3/P-4/P-6，``井名`` 列记 ``合成算例``）= 35 行**。P-3/P-4/P-6 **不复制逐井行**：
它们在结果对象上不可测（R155），把同一份合成算例实测值挂到 8 个井名下会假造
"逐井测过"的印象。

口径声明（复核者必读）
----------------------
1. **短窗非生产数字**：本脚本用 ``nz=60 / ny=12`` 的**短窗**网格（与 Task 6 的
   ``verify_internal_consistency.py`` 同窗同网格，便于横向对照）。其 η 值**不得**当作
   论文/验证数字、**不得**与现场 CBL 比对。生产口径是 ``cemdisp/runners/*_tailpipe.py``
   （nz=250）与 ``scripts/entrypoints/rerun_all_wells_corrected.py``。
2. **1D/2D 口径**：1D 用生产 runner 口径（``enable_gravity`` + T1 三开关
   ``mixing_contact_time`` / ``plug_face_zero_mixing`` / ``has_plug``），2D 用环空
   **默认**开关（不加 ``CORRECTED_KW``）。停算时长 ``total_t = min(泵注总时长 + 1200s,
   尾浆到鞋时刻)``——常量与函数**直接 import 自** ``verify_internal_consistency``，
   避免口径漂移。
3. **P-7 需要两次运行**：基线 + 排量 ×1.4。排量缩放**只改施工程序的泵注步**
   （``PumpingSchedule.steps`` 的 ``rate_m3_min``，体积不变），复用的正是
   本脚本内联的 ``scale_schedule``（语义与 ``scripts/entrypoints/run_sensitivity_current_20260916.py:85``
   逐字一致——只改 ``rate_m3_min``，体积不变；**不**从该脚本导入，以免把
   ``rerun_all_wells_corrected``/``CORRECTED_KW`` 拖进本条刻意使用默认开关的脚本）。
   缩放后两次运行的窗口长度不同（泵注更快 ⇒ ``total_t`` 更短）——这是"排量 ×1.4 时域内"的
   题中之义，逐井把两个 ``total_t`` 写进说明列。
4. **测点层级两类**（协调者裁定 R155）：P-1/P-2/P-5/P-7 在**八井结果对象**上测；
   P-3/P-4/P-6 在**单元级合成算例**上测（``AnnulusSimulationResult`` 不导出速度/通量/
   流函数场），说明列逐行标注层级，台账中这三行 ``井名`` 记为 ``合成算例``。
   P-6 按协调者 2026-09-26 裁定改判：判定只落在**完全冻结格**（``wall ≥ 1−1e-6``，
   即 τw=0），过渡带（``0.5<wall<1−1e-6``）转为同一行的**披露量**；无完全冻结格 ⇒
   未测（非空过守卫）。理由：``wall`` 自提交 ``9448572`` 起是**连续**冻结度。
5. **台账是诊断输出，不是可调参数**：任何一项不通过**如实记录**，并在脚本末尾的
   结论段与 task-9-report.md 写明；**不得**为了让台账"全绿"而放宽判据、改阈值或跳过井。

写盘边界：只写 ``results/内部自洽加固_2026-09-25/``；不碰
``results/<井名>_1D2D耦合模型/``（权威目录），不覆盖任何既有文件。
"""
from __future__ import annotations

import csv
import sys
import time
import traceback
from dataclasses import replace
from pathlib import Path
from typing import Any, Optional

PROJECT_ROOT = Path(__file__).resolve().parents[2]
_SCRIPTS_DIR = PROJECT_ROOT / "scripts"
for _path in (str(PROJECT_ROOT), str(_SCRIPTS_DIR)):
    if _path not in sys.path:
        sys.path.append(_path)

import numpy as np  # noqa: E402

import cemdisp.data.loaders as L  # noqa: E402
from cemdisp.diagnostics.physical_sanity import (  # noqa: E402
    CHECK_BOUNDED,
    CHECK_MONOTONE,
    CHECK_NARROW_DISADVANTAGE,
    CHECK_RATE_RESPONSE,
    sanity_checks,
    synthetic_stream_case,
)
from cemdisp.models2d import AnnulusD2DGASolver  # noqa: E402
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider  # noqa: E402
from cemdisp.transport1d import CasingFlowSolver  # noqa: E402

# 口径唯一来源：与 Task 6 台账共用同一组常量与停算时刻函数（防漂移）。
# ⚠️ 排量缩放**不**从 run_sensitivity_current_20260916 导入：那条 import 会把
# rerun_all_wells_corrected/CORRECTED_KW 拖进本脚本，而本脚本的 2D 求解**刻意**用默认
# 开关，两处口径容易被读成同源（评审 Fix，2026-09-26）。此处内联同一语义的 3 行变换。
from entrypoints.verify_internal_consistency import (  # noqa: E402
    NY,
    NZ,
    TAIL_MARGIN_S,
    stop_time_s,
    total_pump_time_s,
)

OUT_DIR = PROJECT_ROOT / "results" / "内部自洽加固_2026-09-25"
OUT_CSV = OUT_DIR / "物理合理性台账.csv"

COLUMNS = ["井名", "检查项", "通过", "实测值", "判据", "说明"]

# 与 runner / Task 6 台账一致的 8 井（内部代号）
WELLS = [
    ("hu101", L.load_hu101_tailpipe), ("hu102", L.load_hu102_tailpipe),
    ("hu103", L.load_hu103_tailpipe), ("hu1", L.load_hu1_tailpipe),
    ("hu2", L.load_hu2_tailpipe), ("ht1_001", L.load_ht1_001_tailpipe),
    ("ht1_003", L.load_ht1_003_tailpipe), ("ht1_004", L.load_ht1_004_tailpipe),
]

RATE_FACTOR = 1.4          # P-7 的对照排量倍数（现场经验方向检查）
# 结果对象级判据（有逐井行）；P-3/P-4/P-6 只在"合成算例"测点上有行（R155）。
# 名字一律取自 physical_sanity 的常量（单一真源，禁止在本文件里另写一份字面量）。
RESULT_LEVEL_CHECKS = (CHECK_BOUNDED, CHECK_MONOTONE, CHECK_NARROW_DISADVANTAGE,
                       CHECK_RATE_RESPONSE)
SYNTHETIC_LABEL = "合成算例"  # P-3/P-4/P-6 的测点层级标签（非井）


def render_passed(passed: Any) -> str:
    """把 ``通过`` 渲染成台账文本（None ⇒ 未测，绝不写成通过）。"""
    if passed is None:
        return "未测"
    return "通过" if passed else "不通过"


def scale_schedule(schedule: Any, factor: float) -> Any:
    """泵序逐段排量缩放（体积不变，时长随 1/factor 变化）。

    与 ``scripts/entrypoints/run_sensitivity_current_20260916.scale_schedule`` **逐字同语义**
    （只改 ``rate_m3_min``，排量 ≤ 0 的段保持不动），此处内联以免把
    ``rerun_all_wells_corrected``/``CORRECTED_KW`` 一并拖入本脚本的 import 面
    （本脚本的 2D 求解刻意用默认开关）。
    """
    return replace(schedule, steps=tuple(
        replace(st, rate_m3_min=st.rate_m3_min * factor) if st.rate_m3_min > 0 else st
        for st in schedule.steps))


def run_case(loader: Any, factor: float = 1.0) -> tuple:
    """跑一口井的一次完整流水线（1D 生产口径 + 2D 短窗），返回 (result, 诊断信息)。"""
    well, fluids, schedule, _ = loader()
    if factor != 1.0:
        schedule = scale_schedule(schedule, factor)   # 只改泵注步排量，体积不变
    casing = CasingFlowSolver(
        enable_gravity=True,
        mixing_contact_time=True,
        plug_face_zero_mixing=True,
        has_plug=True,
    )
    casing_result = casing.run(well, fluids, schedule)
    t_stop_s = stop_time_s(casing_result, fluids)
    if t_stop_s is None:
        raise RuntimeError(
            "无法确定 2D 停算时刻：cement_end_time_s 缺失且前缘扫描无回退")
    pump_time_s = total_pump_time_s(schedule)
    total_t_s = min(pump_time_s + TAIL_MARGIN_S, t_stop_s)
    inlet = build_coupled_annulus_inlet_provider(
        casing_result, CasingFlowSolver(enable_gravity=True), fluids,
        split_cement_phases=True)
    result = AnnulusD2DGASolver(total_t=total_t_s, nz=NZ, ny=NY).run(
        well, fluids, inlet)
    return result, {"total_t_s": total_t_s, "pump_s": pump_time_s, "stop_s": t_stop_s}


def well_rows(label: str, loader: Any, case: dict, frozen_case: dict) -> list:
    """单井 4 行（P-1/P-2/P-5/P-7）：基线 + 排量 ×1.4 两次运行。

    P-3/P-4/P-6 **没有逐井行**：它们在结果对象上不可测（R155），测点是本台账末尾的
    ``合成算例`` 三行（单元级合成算例）。把同一份合成算例实测值复制到 8 个井名下
    会假造"逐井测过"的印象，故不复制。
    """
    base, base_info = run_case(loader, 1.0)
    rate, rate_info = run_case(loader, RATE_FACTOR)
    rows = [r for r in sanity_checks(base, base.geom, rate_result=rate,
                                     stream_case=case, frozen_case=frozen_case)
            if r["检查项"] in RESULT_LEVEL_CHECKS]
    for row in rows:
        row["井名"] = label
        if row["检查项"] == CHECK_RATE_RESPONSE:
            row["说明"] += (
                f"; 两次运行窗口：基线 total_t={base_info['total_t_s']:.0f}s"
                f"（泵注{base_info['pump_s']:.0f}s/到鞋{base_info['stop_s']:.0f}s）"
                f" vs 排量×{RATE_FACTOR} total_t={rate_info['total_t_s']:.0f}s"
                f"（泵注{rate_info['pump_s']:.0f}s/到鞋{rate_info['stop_s']:.0f}s）"
                f"；排量缩放只改泵注步（本脚本内联的 scale_schedule，语义同 run_sensitivity_current_20260916:85）")
    return rows


def synthetic_rows(case: dict, frozen_case: dict) -> list:
    """合成算例的 3 行（P-3/P-4/P-6）：结果对象侧传入零快照 ⇒ P-1/P-2/P-5/P-7 不参与。"""
    stub = _EmptyResult()
    rows = sanity_checks(stub, case["geom"], stream_case=case, frozen_case=frozen_case)
    wanted = {"P-3 速度上界", "P-4 总通量守恒", "P-6 冻结区静止"}
    out = []
    for row in rows:
        if row["检查项"] in wanted:
            row["井名"] = SYNTHETIC_LABEL
            out.append(row)
    return out


class _EmptyResult:
    """空结果替身：合成算例行只取 P-3/P-4/P-6，其余项在本替身上必然"未测"。"""

    cement_snapshots: tuple = ()
    summary: dict = {}


def failed_rows(label: str, exc: BaseException,
                criteria: Optional[dict] = None) -> list:
    """失败井：结果对象级判据全部如实标记未测并写明失败原因（不允许静默消失）。

    ``criteria``：``{检查项: 判据文本}``，由成功井（或空替身）上跑一次判据收集而得；
    这样失败行也**带着判据文本**，读者不必回头查别的井才知道本行在判什么。
    """
    reason = f"未测：本井运行失败 {type(exc).__name__}: {exc}"
    criteria = criteria or {}
    return [{"井名": label, "检查项": name, "通过": "未测", "实测值": "未测",
             "判据": criteria.get(name, "—"), "说明": reason}
            for name in RESULT_LEVEL_CHECKS]


def write_ledger(rows: list, path: Path = OUT_CSV) -> Path:
    """按 ``COLUMNS`` 写出台账 CSV（utf-8-sig，与仓内其它产物一致）。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({"井名": row["井名"], "检查项": row["检查项"],
                             "通过": render_passed(row["通过"]),
                             "实测值": row["实测值"], "判据": row["判据"],
                             "说明": row["说明"]})
    return path


def main() -> int:
    print(f"物理合理性台账（短窗 nz={NZ}/ny={NY}，仅作自洽诊断）")
    print(f"输出：{OUT_CSV}")
    case = synthetic_stream_case(freeze=None)                  # P-3/P-4 基线（无屈服门）
    frozen_case = synthetic_stream_case(freeze="smooth")       # P-6 连续冻结度算例
    # 判据文本表：在空替身上收集一次，供失败井行原样沿用（判据不进"未测"）
    criteria = {r["检查项"]: r["判据"] for r in sanity_checks(
        _EmptyResult(), case["geom"], stream_case=case, frozen_case=frozen_case)}
    print(f"[合成算例] {case['case']}")
    print(f"[合成算例/冻结] {frozen_case['case']}")

    rows: list = []
    n_failed = 0
    for label, loader in WELLS:
        t0 = time.perf_counter()
        try:
            got = well_rows(label, loader, case, frozen_case)
            status = "OK"
        except Exception as exc:  # noqa: BLE001 —— 逐井隔离，失败必须留痕
            traceback.print_exc()
            got = failed_rows(label, exc, criteria)
            n_failed += 1
            status = "FAIL"
        rows.extend(got)
        write_ledger(rows)      # 逐井落盘，长跑中断也留下已完成部分
        verdict = " ".join(f"{r['检查项'][:3]}={render_passed(r['通过'])}" for r in got)
        print(f"[{status}] {label}: {verdict} ({time.perf_counter() - t0:.1f}s)")

    rows.extend(synthetic_rows(case, frozen_case))
    write_ledger(rows)
    print(f"\n台账已写入 {OUT_CSV}（{len(rows)} 行，失败 {n_failed} 口）")
    return _conclusion(rows)


def _conclusion(rows: list) -> int:
    """结论段：逐条列出不通过与未测（诊断输出，不是验收门；不因不通过而退出非零）。"""
    bad = [r for r in rows if r["通过"] is False]
    na = [r for r in rows if r["通过"] is None]
    print(f"\n=== 结论 ===")
    print(f"总行数 {len(rows)}；通过 {sum(1 for r in rows if r['通过'] is True)}；"
          f"不通过 {len(bad)}；未测 {len(na)}")
    if not bad:
        print("无判据被判不通过。")
    for r in bad:
        print(f"❌ 不通过：{r['井名']} / {r['检查项']}：{r['实测值']}"
              f"（判据 {r['判据']}）")
    for r in na:
        print(f"⚠️ 未测：{r['井名']} / {r['检查项']}：{r['说明']}")
    print("\n声明：本脚本用短窗 nz=60/ny=12 网格，产物仅供自洽诊断，"
          "不得当作论文/验证数字。判据与阈值一律按计划 Task 9 原值，未做任何放宽。")
    return n_failed_count(rows)


def n_failed_count(rows: list) -> int:
    """返回**失败井数**（按 井名 去重；一口失败井有多行，不能把行数当井数）。"""
    return len({r["井名"] for r in rows
                if r["说明"].startswith("未测：本井运行失败")})


if __name__ == "__main__":
    raise SystemExit(main())
