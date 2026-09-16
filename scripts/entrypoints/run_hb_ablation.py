"""阶段 B（Task 4）八井 A/B：量化流变口径对 η_E/η_N 的净影响（填补既有未测量缺口）。

被量化的两个开关均为阶段 B 新增、**默认关**、opt-in：

- ``enable_stream_yield_gate``      （B-2）屈服门进 (4.22) 流函数算子；
- ``enable_power_law_gap_correction``（B-3）幂律间隙一阶修正 I₁·(H/H̄)^{1/n−1}。

变体
----
==============  ==========================================================
V0_baseline     默认（两开关关，``enable_stream_function=True``）
V1_yield_gate   仅 ``enable_stream_yield_gate=True``
V2_power_law    仅 ``enable_power_law_gap_correction=True``
V3_both         两开关同开
V4_legacy       旧代数流动度路径（``enable_stream_function=False``）
==============  ==========================================================

读法：**V1/V2/V3 相对 V0 的 Δη 即各开关（及其组合）的净效应**；
**V0 vs V4 是"流函数路径 vs 旧代数路径"的联合效应**（含流变口径与速度场构造）。

口径与边界
----------
- **相对对照，非权威绝对数字**：本脚本用粗网格 ``nz=60``（权威 8 井数字为
  ``nz=250``，入口 ``scripts/entrypoints/rerun_all_wells_corrected.py``）。
  粗网格用于快速横向比较变体间差异；**不得**把本表数字当作论文/验证数字引用。
- 停泵时刻沿权威运行器口径（F2 修复，2026-09-01）：``total_t = min(施工总时长 +
  1200s, 尾浆全部入环空时刻 cement_end_time_s)``——不用固定 3600s，否则各井
  （施工时长 8.0k–20.9k s）会被截断成无意义的近零顶替。
- 入口边界沿权威口径：``CasingFlowSolver(enable_gravity=True)`` → 鞋口时间线 →
  ``build_coupled_annulus_inlet_provider(..., split_cement_phases=True)``。
  逐变体重建（套管 1D 求解开销 ~ms 级），避免边界提供器跨变体复用引入状态污染。
- ⚠️ **ht1_004 的 V2 是开关空转**：该井 LEAD/TAIL 均为 BINGHAM（``power_law_n is
  None``），B-3 运行时一次性告警后 n_rep 回落 1.0、因子恒 1 ⇒ V2 与 V0 逐位相同。
  **该行的 Δη≈0 是"开关空转"而非"物理中性"**——CSV 的 ``pl_idle`` 列逐行标注该
  变体是否触发了空转告警。见 docs/源模型口径与适用域声明.md 声明 3b。

输出
----
``results/流变口径A-B_2026-09-16/ablation_8wells.csv``，列：
``well, variant, eta_E, eta_N, b_number, n_cement, pl_idle``。
（``n_cement`` = 水泥相代表幂律指数，None 记空；``pl_idle`` = 该次运行是否触发
"静默空转"告警。）

用法::

    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/entrypoints/run_hb_ablation.py
"""
from __future__ import annotations

import csv
import time
import warnings
from pathlib import Path
from typing import cast

from cemdisp.data import loaders
from cemdisp.data.fluid_spec import FluidRole
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = PROJECT_ROOT / "results" / "流变口径A-B_2026-09-16"
NZ = 60  # 粗网格快速对照（权威口径为 250）

# 显式列出 8 井（与 rerun_all_wells_corrected.py 同序），不用 dir() 反射——
# 避免把 *_design / *_actual 变体混入（brief 的 dir() 反射会带进这些）。
WELLS = [
    ("hu101", loaders.load_hu101_tailpipe),
    ("hu102", loaders.load_hu102_tailpipe),
    ("hu103", loaders.load_hu103_tailpipe),
    ("hu1", loaders.load_hu1_tailpipe),
    ("hu2", loaders.load_hu2_tailpipe),
    ("ht1_001", loaders.load_ht1_001_tailpipe),
    ("ht1_003", loaders.load_ht1_003_tailpipe),
    ("ht1_004", loaders.load_ht1_004_tailpipe),
]

VARIANTS: dict[str, dict[str, bool]] = {
    "V0_baseline": {},
    "V1_yield_gate": {"enable_stream_yield_gate": True},
    "V2_power_law": {"enable_power_law_gap_correction": True},
    "V3_both": {"enable_stream_yield_gate": True, "enable_power_law_gap_correction": True},
    "V4_legacy": {"enable_stream_function": False},
}

_CEMENT_ROLES = (FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL)
_IDLE_MARK = "静默空转"


def _schedule_total_time_s(schedule) -> float:
    return sum(0.0 if s.rate_m3_min <= 0 else s.volume_m3 / s.rate_m3_min * 60.0
               for s in schedule.steps)


def _cement_representative_n(fluids) -> float | None:
    """水泥相代表幂律指数（lead 优先，缺则取任一水泥相）——与 solver 内部同口径。"""
    cement = [f for f in fluids if f.role in _CEMENT_ROLES]
    for f in cement:
        if f.role == FluidRole.LEAD and f.power_law_n is not None:
            return float(f.power_law_n)
    for f in cement:
        if f.power_law_n is not None:
            return float(f.power_law_n)
    return None


def _annulus_total_t(schedule, casing_result) -> float:
    """停泵时刻（权威口径，F2 修复 2026-09-01）：min(施工总时长+1200s, 尾浆全部入环空)。"""
    total_t = _schedule_total_time_s(schedule)
    stop = (float(casing_result.cement_end_time_s)
            if casing_result.cement_end_time_s is not None else total_t + 1200.0)
    return min(total_t + 1200.0, stop)


def _build_inlet(well, fluids, schedule):
    """按权威口径构建环空入口提供器（每次调用重新求解套管 1D，避免跨变体状态污染）。"""
    casing_result = CasingFlowSolver(enable_gravity=True).run(well, fluids, schedule)
    inlet = build_coupled_annulus_inlet_provider(
        casing_result, CasingFlowSolver(enable_gravity=True), fluids,
        split_cement_phases=True,
    )
    return inlet, _annulus_total_t(schedule, casing_result)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    t_start = time.perf_counter()
    for wname, loader in WELLS:
        well, fluids, schedule, _ = loader()          # 第 4 元是 ValidationData，非 provider
        n_cement = _cement_representative_n(fluids)
        # total_t 各变体共用（同井同边界），保证变体间可比；入口 provider 逐变体重建
        _, total_t = _build_inlet(well, fluids, schedule)
        print(f"\n=== {wname} (nz={NZ}, total_t={total_t:.0f}s, n_cement={n_cement}) ===", flush=True)
        for vname, kw in VARIANTS.items():
            t0 = time.perf_counter()
            inlet_v, _ = _build_inlet(well, fluids, schedule)   # 逐变体重建边界
            solver = AnnulusD2DGASolver(total_t=total_t, nz=NZ, **kw)
            with warnings.catch_warnings(record=True) as rec:
                warnings.simplefilter("always")
                res = solver.run(well, fluids, inlet_v, schedule=schedule)
            idle = any(_IDLE_MARK in str(r.message) for r in rec)
            fr = cast(dict, res.summary["最终结果"])
            row = {
                "well": wname, "variant": vname,
                "eta_E": float(fr["全井段最终有效顶替效率"]),
                "eta_N": float(fr["窄四分位效率"]),
                "b_number": float(fr["浮力数_b"]),
                "n_cement": n_cement,
                "pl_idle": idle,
            }
            rows.append(row)
            print(f"  {wname:8s} {vname:12s} eta_E={row['eta_E']:.4f} "
                  f"eta_N={row['eta_N']:.4f} b={row['b_number']:+.3f}"
                  f"{'  [V2 空转]' if idle else ''}  ({time.perf_counter() - t0:.1f}s)",
                  flush=True)

    csv_path = OUT_DIR / "ablation_8wells.csv"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["well", "variant", "eta_E", "eta_N",
                                          "b_number", "n_cement", "pl_idle"])
        w.writeheader(); w.writerows(rows)

    # ---- 控制台宽表：每井 V0 为基准的 Δ ----
    by = {(r["well"], r["variant"]): r for r in rows}
    others = [v for v in VARIANTS if v != "V0_baseline"]
    print("\n" + "=" * 96)
    hdr = "".join(f"{('Δ' + v.split('_')[0]):>8}" for v in others)
    print(f"{'井':<9}{'V0 η_E':>9}{hdr} |{'V0 η_N':>9}{hdr}")
    for wname, _ in WELLS:
        b0 = by.get((wname, "V0_baseline"))
        if b0 is None:
            print(f"{wname:<9}  [无 V0 基准]"); continue
        de = [by[(wname, v)]["eta_E"] - b0["eta_E"] for v in others]
        dn = [by[(wname, v)]["eta_N"] - b0["eta_N"] for v in others]
        print(f"{wname:<9}{b0['eta_E']:>9.4f}" + "".join(f"{x:>+8.4f}" for x in de)
              + f" |{b0['eta_N']:>9.4f}" + "".join(f"{x:>+8.4f}" for x in dn))
    print("=" * 96)
    print(f"CSV: {csv_path}")
    print(f"总耗时 {time.perf_counter() - t_start:.0f}s（{len(rows)} 次求解）")

    # ---- 开关生效自检（须在 40 次求解之后仍成立）----
    _self_checks(by, rows)


def _self_checks(by: dict, rows: list[dict]) -> None:
    """确认两开关在真实运行中确实生效，且空转井被如实标注。"""
    print("\n[自检]")
    # 1) 屈服门：至少一井 V1 ≠ V0（否则开关未接线/被吞）
    if ("hu101", "V1_yield_gate") in by:
        moved = [w for w, _ in WELLS
                 if (w, "V1_yield_gate") in by and (
                     by[(w, "V1_yield_gate")]["eta_E"] != by[(w, "V0_baseline")]["eta_E"]
                     or by[(w, "V1_yield_gate")]["eta_N"] != by[(w, "V0_baseline")]["eta_N"])]
        print(f"  屈服门生效井数（V1≠V0）: {len(moved)}/{len(WELLS)} {moved}")
        assert moved, "enable_stream_yield_gate 在 8 井上全部无效应——开关疑似未生效"
    # 2) 幂律修正：幂律井 V2 应变化；非幂律井 V2 须逐位等于 V0 且已告警空转
    for wname, _ in WELLS:
        b0, v2 = by.get((wname, "V0_baseline")), by.get((wname, "V2_power_law"))
        if b0 is None or v2 is None:
            continue
        if b0["n_cement"] is None:
            same = (v2["eta_E"] == b0["eta_E"] and v2["eta_N"] == b0["eta_N"])
            print(f"  {wname}: n_cement=None ⇒ V2 空转，逐位等于 V0 = {same}，"
                  f"告警已触发 = {v2['pl_idle']}")
            assert same and v2["pl_idle"], f"{wname} 非幂律但 V2 未空转/未告警——判据失效"
        else:
            print(f"  {wname}: n_cement={b0['n_cement']:.3f} ⇒ V2 "
                  f"Δη_E={v2['eta_E'] - b0['eta_E']:+.4f} Δη_N={v2['eta_N'] - b0['eta_N']:+.4f}")
    print("[自检] 通过")


if __name__ == "__main__":
    main()
