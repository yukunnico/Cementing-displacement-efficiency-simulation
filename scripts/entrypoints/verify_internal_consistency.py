# -*- coding: utf-8 -*-
"""内部自洽校验台账：域内效率 / 饥饿份额恒等式 / 前缘位置 / 质量守恒（Task 6）。

产出
----
``results/内部自洽加固_2026-09-25/一致性台账.csv``，列为：
``井名, 域内eta_E, 饥饿份额, 恒等式偏差, 1D尾浆到鞋时刻_s, 2D前缘位置_m, 质量守恒误差, 说明``

口径声明（复核者必读）
----------------------
1. **短窗非生产数字**：本脚本用 ``nz=60 / ny=12`` 的**短窗**网格，只为让三项自洽量跑得动、
   逐井可比；其 η 值**不得**当作论文/验证数字、**不得**与现场 CBL 比对。生产口径是
   ``cemdisp/runners/*_tailpipe.py``（nz=250）与 ``scripts/entrypoints/rerun_all_wells_corrected.py``。
   ⚠️ **该列不是"误差"，是过渡带的代数和**：对连续场可严格展开为
   ``Σ_{c≥0.5} b(1−c)/∬b − Σ_{c<0.5} b·c/∬b``（与 ``(1−η_E) − 饥饿份额`` 逐位相等，
   实测差 0.0，见 task-6-report.md），二值场时恒为 0。网格加密只会收窄过渡带、
   不会把连续场变二值：实测 hu101 nz=60 → 4.29e-2、nz=120 → 3.27e-2、nz=250 → 2.62e-2；
   而**基准算例**（纯水泥恒定入口 ⇒ 场近乎二值）实测 6.7e-10 ~ 1.0e-3 ⇒ 历史口径
   "偏差 ≤0.007" 的出处是基准算例，**不是**现场井。故本列只用于**同窗横向比较**。
2. **1D/2D 口径**：1D 用生产 runner 口径（`enable_gravity` + T1 三开关
   `mixing_contact_time` / `plug_face_zero_mixing` / `has_plug`），2D 用环空**默认**开关
   （不加 `CORRECTED_KW`），停算时刻 = ``CasingFlowResult.cement_end_time_s``
   （尾浆全部进入环空的时刻；缺失时按 `rerun_all_wells_corrected._stop_t` 的前缘扫描回退）。
3. **质量守恒列写"未测"的理由**（协调者裁定 #4：假设不成立即如实写未测，不得自造公式）：
   ``zhang2022_benchmark.mass_conservation_error``（:430）的口径是
   ``max |V_ann·bulk_cement_fill − Q·t| / (Q·t)``，其成立前提有两条：
   ① **单一恒定排量** ``q_m3s``；② **入口自 t=0 恒为纯水泥**（基准算例 `build_inlet_provider`
   即 ``AnnulusInletState(..., (("tail", 1.0),))`` + 恒定 ``Q̂0``，见 :252）。
   本仓 8 口现场井**两条都不成立**：施工程序为分段变排量（实测每井 2–9 种排量，
   如 hu101 = 0.55/0.6/1.0/1.2/1.5 m³/min），且入口序为
   钻井液/隔离液/领浆/尾浆/替浆 ⇒ ``bulk_cement_fill`` 只计水泥而 ``Q·t`` 计全部注入流体。
   套用该函数会得到无物理含义的数（既非守恒误差也非约定残差），故该列一律写 ``未测``。

写盘边界：只写 ``results/内部自洽加固_2026-09-25/``（本计划自己的输出目录）；
不碰 ``results/<井名>_1D2D耦合模型/``（权威目录），不覆盖任何既有文件。
"""
from __future__ import annotations

import csv
import sys
import time
import traceback
from pathlib import Path
from typing import Any

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import cemdisp.data.loaders as L  # noqa: E402
from cemdisp.data.fluid_spec import FluidRole  # noqa: E402
from cemdisp.diagnostics.internal_consistency import (  # noqa: E402
    domain_eta_e,
    front_position_m,
    starved_volume_fraction,
)
from cemdisp.models2d import AnnulusD2DGASolver  # noqa: E402
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider  # noqa: E402
from cemdisp.transport1d import CasingFlowSolver  # noqa: E402

OUT_DIR = PROJECT_ROOT / "results" / "内部自洽加固_2026-09-25"
OUT_CSV = OUT_DIR / "一致性台账.csv"

# 短窗诊断网格（**非生产口径**，见模块 docstring 第 1 条）
NZ = 60
NY = 12

COLUMNS = ["井名", "域内eta_E", "饥饿份额", "恒等式偏差", "1D尾浆到鞋时刻_s",
           "2D前缘位置_m", "质量守恒误差", "说明"]

# 与 runner 一致的 8 井（内部代号，与 Task 1/3/7 台账口径统一）
WELLS = [
    ("hu101", L.load_hu101_tailpipe), ("hu102", L.load_hu102_tailpipe),
    ("hu103", L.load_hu103_tailpipe), ("hu1", L.load_hu1_tailpipe),
    ("hu2", L.load_hu2_tailpipe), ("ht1_001", L.load_ht1_001_tailpipe),
    ("ht1_003", L.load_ht1_003_tailpipe), ("ht1_004", L.load_ht1_004_tailpipe),
]

_NT = "未测"

# 质量守恒列"未测"的简短理由（完整论证见模块 docstring 第 3 条）
_MC_REASON = ("质量守恒未测：mass_conservation_error 前提(单一恒定排量+入口自t=0恒为纯水泥)"
              "对现场分段变排量/多流体入口不成立，按裁定不自造公式")


def stop_time_s(casing_result: Any, fluids: tuple) -> float:
    """2D 停算时刻 = 尾浆全部进入环空（``cement_end_time_s``），缺失时按前缘扫描回退。

    与 ``scripts/entrypoints/rerun_all_wells_corrected._stop_t`` 同口径（F2 修复 2026-09-01）。
    """
    if casing_result.cement_end_time_s is not None:
        return float(casing_result.cement_end_time_s)
    cement_roles = {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}
    by_name = {f.name: f for f in fluids}
    fronts = sorted(casing_result.fronts, key=lambda f: f.time_s)
    last = next((f.time_s for f in fronts
                 if by_name.get(f.fluid_name) is not None
                 and by_name[f.fluid_name].role in cement_roles), None)
    if last is not None:
        for f in fronts:
            fl = by_name.get(f.fluid_name)
            if fl is None or fl.role in cement_roles:
                continue
            if f.time_s >= last - 1.0e-9:
                return float(f.time_s)
    return float("nan")


def row_from_result(label: str, result: Any, t_stop_s: float) -> dict[str, str]:
    """由 2D 结果组装台账行（三项自洽量 + 模型 s 口径前缘，供说明列对照）。"""
    geom = result.geom
    cement = np.asarray(result.cement_field, dtype=float)
    eta = domain_eta_e(cement, geom)
    starved = starved_volume_fraction(cement, geom)
    deviation = (1.0 - eta) - starved
    front_md = front_position_m(cement, geom)
    last = result.metrics.iloc[-1]
    # 同口径自审：域内积分必须与求解器自算的 effective_efficiency（≡ bulk_cement_fill）逐位一致；
    # 不一致说明归一化/口径漂移，如实写进说明列（是诊断提示，不是验收门）。
    solver_eta = float(last["effective_efficiency"])
    notes = [
        f"短窗nz={NZ}/ny={NY}（非生产数字）",
        ("eta_E与solver末行effective_efficiency逐位一致"
         if eta == solver_eta else
         f"⚠️eta_E与solver末行不一致：{eta:.6f} vs {solver_eta:.6f}"),
        f"模型s口径前缘自鞋口(m)：宽{float(last['front_wide_m']):.0f}"
        f"/中{float(last['front_mid_m']):.0f}/窄{float(last['front_narrow_m']):.0f}；"
        f"域长{float(geom['s'][-1]):.0f}",
        f"md口径前缘≡域底{float(geom['md'][0]):.0f}（定义见 internal_consistency.front_position_m）",
        _MC_REASON,
    ]
    if abs(deviation) > 0.05:
        notes.append(f"⚠️恒等式偏差{deviation:.3f}>0.05（非二值场+粗网格双因素，见报告）")
    return {
        "井名": label,
        "域内eta_E": f"{eta:.6f}",
        "饥饿份额": f"{starved:.6f}",
        "恒等式偏差": f"{deviation:.6f}",
        "1D尾浆到鞋时刻_s": f"{t_stop_s:.1f}",
        "2D前缘位置_m": f"{front_md:.2f}",
        "质量守恒误差": _NT,
        "说明": "; ".join(notes),
    }


def failed_row(label: str, exc: BaseException) -> dict[str, str]:
    """失败行：如实记录异常类型与消息（跑不动的井不允许静默消失）。"""
    return {
        "井名": label, "域内eta_E": _NT, "饥饿份额": _NT, "恒等式偏差": _NT,
        "1D尾浆到鞋时刻_s": _NT, "2D前缘位置_m": _NT, "质量守恒误差": _NT,
        "说明": f"运行失败：{type(exc).__name__}: {exc}",
    }


def compute_row(label: str, loader: Any) -> dict[str, str]:
    """单井完整流程：loader → 1D（生产口径）→ 2D 短窗 → 台账行。"""
    well, fluids, schedule, _ = loader()
    casing_solver = CasingFlowSolver(
        enable_gravity=True,
        # T1 生产口径三开关（与 cemdisp/runners/*_tailpipe.py 一致）
        mixing_contact_time=True,
        plug_face_zero_mixing=True,
        has_plug=True,
    )
    casing_result = casing_solver.run(well, fluids, schedule)
    t_stop_s = stop_time_s(casing_result, fluids)
    inlet = build_coupled_annulus_inlet_provider(
        casing_result, CasingFlowSolver(enable_gravity=True), fluids,
        split_cement_phases=True)
    result = AnnulusD2DGASolver(total_t=t_stop_s, nz=NZ, ny=NY).run(well, fluids, inlet)
    return row_from_result(label, result, t_stop_s)


def write_ledger(rows: list[dict[str, str]], path: Path = OUT_CSV) -> Path:
    """按 ``COLUMNS`` 写出台账 CSV（utf-8-sig，与仓内其它产物一致）。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def main() -> int:
    """逐井跑短窗自洽诊断并写出台账；返回失败井数。"""
    print(f"内部自洽校验台账（短窗 nz={NZ}/ny={NY}，仅作自洽诊断）")
    print(f"输出：{OUT_CSV}")
    rows: list[dict[str, str]] = []
    n_failed = 0
    for label, loader in WELLS:
        t0 = time.perf_counter()
        try:
            row = compute_row(label, loader)
            print(f"[OK]   {label}: eta_E={row['域内eta_E']} 偏差={row['恒等式偏差']} "
                  f"({time.perf_counter() - t0:.1f}s)")
        except Exception as exc:  # noqa: BLE001 —— 逐井隔离，失败必须留痕
            traceback.print_exc()
            row = failed_row(label, exc)
            n_failed += 1
            print(f"[FAIL] {label}: {type(exc).__name__}: {exc}")
        rows.append(row)
        write_ledger(rows)  # 逐井落盘，长跑中断也留下已完成部分
    print(f"\n台账已写入 {OUT_CSV}（{len(rows)} 口井，失败 {n_failed} 口）")
    for row in rows:
        print(f"  {row['井名']:<8} eta_E={row['域内eta_E']:<9} 饥饿={row['饥饿份额']:<9} "
              f"偏差={row['恒等式偏差']:<10} 守恒={row['质量守恒误差']}")
    return n_failed


if __name__ == "__main__":
    raise SystemExit(main())
