"""呼101井居中度实测值 vs 模型假设值 敏感性对比。

基于居中度检测图（Pipe Standoff）提取的实测剖面，与模型 legacy 假设
剖面（0.38–0.48）对比，定量评估居中度输入不确定对顶替效率的影响。

实测剖面从检测图读取，取两套：
  A. 扶正器之间 (Between Cent.) — 代表大部分管长的实际居中度（偏下限）
  B. 扶正器处 (At Cent.)     — 代表扶正器位置的居中度（偏上限）

模型域：5400–7868 m（尾管段）
运行模式：1D-2D 耦合（与生产口径一致）
网格：nz=250（呼101 legacy 口径）
"""
from __future__ import annotations

import csv
import json
import time
from dataclasses import replace
from pathlib import Path
from typing import cast

import numpy as np

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.loaders import load_hu101_tailpipe
from cemdisp.data.pumping_schedule import PumpingSchedule
from cemdisp.data.well_spec import DepthValuePoint, WellSpec
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.reporting.plots import plot_depth_profiles
from cemdisp.transport1d import CasingFlowResult, CasingFlowSolver


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PROJECT_ROOT / "results" / "呼101_居中度实测对比"

# 模型域深度范围（米）
TOP_MD = 5400.0
BOTTOM_MD = 7868.0

# ============================================================
# 实测居中度剖面（从居中度检测图读取，深度-居中度 对）
# 横坐标 Pipe Standoff 0–1；深度对应左侧 Well 图纵轴
# ============================================================

# 扶正器之间 Between Cent.（红色横线读取值）
MEASURED_BETWEEN_CENT = [
    (5400.0, 0.78),   # 最顶部
    (5700.0, 0.72),
    (6000.0, 0.65),
    (6300.0, 0.58),
    (6600.0, 0.52),
    (6900.0, 0.48),
    (7200.0, 0.42),
    (7500.0, 0.32),
    (7700.0, 0.25),
    (7868.0, 0.22),   # 最底部附近
]

# 扶正器处 At Cent.（蓝色十字读取值）
MEASURED_AT_CENT = [
    (5400.0, 0.88),
    (5700.0, 0.85),
    (6000.0, 0.80),
    (6300.0, 0.76),
    (6600.0, 0.72),
    (6900.0, 0.70),
    (7200.0, 0.68),
    (7500.0, 0.65),
    (7700.0, 0.62),
    (7868.0, 0.60),
]

# 模型假设（legacy）剖面 — 与 hu101_loader.py 完全一致
ASSUMED_PROFILE = [
    (5400.0, 0.45),
    (6100.0, 0.38),
    (6796.0, 0.44),
    (7200.0, 0.48),
    (7600.0, 0.42),
    (7868.0, 0.46),
]


def _schedule_total_time_s(schedule: PumpingSchedule) -> float:
    return sum(
        0.0 if step.rate_m3_min <= 0.0 else step.volume_m3 / step.rate_m3_min * 60.0
        for step in schedule.steps
    )


def _annulus_stop_time_s(casing_result: CasingFlowResult, fluids: tuple[FluidSpec, ...]) -> float:
    cement_roles = {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}
    fluid_by_name = {f.name: f for f in fluids}
    sorted_fronts = sorted(casing_result.fronts, key=lambda f: f.time_s)
    last_cement_time_s: float | None = None
    for front in sorted_fronts:
        fluid = fluid_by_name.get(front.fluid_name)
        if fluid is not None and fluid.role in cement_roles:
            last_cement_time_s = front.time_s
    if last_cement_time_s is not None:
        for front in sorted_fronts:
            fluid = fluid_by_name.get(front.fluid_name)
            if fluid is None or fluid.role in cement_roles:
                continue
            if front.time_s >= last_cement_time_s - 1.0e-9:
                return float(front.time_s)
    return float(casing_result.cement_end_time_s)


def make_well_with_standoff(base_well: WellSpec, profile_points: list[tuple[float, float]]) -> WellSpec:
    profile = tuple(DepthValuePoint(md, val) for md, val in profile_points)
    return replace(base_well, standoff_profile=profile)


def run_case(
    case_name: str,
    well_spec: WellSpec,
    fluids: tuple[FluidSpec, ...],
    schedule: PumpingSchedule,
    nz: int = 250,
) -> dict:
    """运行单个 case，返回指标字典。"""
    print(f"\n{'='*60}")
    print(f"  Case: {case_name}")
    print(f"{'='*60}")

    t0 = time.perf_counter()

    # 1D 套管内
    casing_solver = CasingFlowSolver(enable_gravity=True)
    casing_result = casing_solver.run(well_spec, fluids, schedule)
    annulus_stop = _annulus_stop_time_s(casing_result, fluids)

    # 耦合入口
    inlet_provider = build_coupled_annulus_inlet_provider(
        casing_result, casing_solver, fluids, split_cement_phases=True,
    )

    # 2D 环空
    total_t = _schedule_total_time_s(schedule) + 20.0 * 60.0
    total_t = min(total_t, annulus_stop + 10.0 * 60.0)

    solver = AnnulusD2DGASolver(total_t=total_t, nz=nz)
    result = solver.run(well_spec, fluids, inlet_provider)

    elapsed = time.perf_counter() - t0

    final = cast(dict, result.summary["最终结果"])
    metrics = {
        "case": case_name,
        "effective_efficiency": float(final["全井段最终有效顶替效率"]),
        "cement_occupation": float(final["最终水泥浆占据率"]),
        "channeling_index": float(final["最终窜槽指数"]),
        "mixing_index": float(final["最终混浆指数"]),
        "instability_index": float(final["最终失稳指数"]),
        "elapsed_s": round(elapsed, 1),
    }

    # 剖面均值
    prof = result.depth_profiles
    if "居中度" in prof.columns:
        metrics["mean_standoff"] = round(float(np.mean(prof["居中度"].values)), 4)

    print(f"  有效效率:   {metrics['effective_efficiency']:.4f}")
    print(f"  水泥占据率: {metrics['cement_occupation']:.4f}")
    print(f"  窜槽指数:   {metrics['channeling_index']:.4f}")
    print(f"  混浆指数:   {metrics['mixing_index']:.4f}")
    print(f"  失稳指数:   {metrics['instability_index']:.4f}")
    print(f"  均居中度:   {metrics.get('mean_standoff', 'N/A')}")
    print(f"  耗时:       {elapsed:.1f}s")

    # 保存子目录数据
    case_dir = OUTPUT_DIR / case_name
    case_dir.mkdir(parents=True, exist_ok=True)
    _ = result.metrics.to_csv(case_dir / "时间序列.csv", index=False, encoding="utf-8-sig")
    _ = result.depth_profiles.to_csv(case_dir / "深度剖面.csv", index=False, encoding="utf-8-sig")
    _ = plot_depth_profiles(result, well_spec=well_spec, output_dir=case_dir)

    np.savez(
        case_dir / "2D场数据.npz",
        cement_final=result.cement_field,
        spacer_final=result.spacer_field,
        wall_final=result.wall_field,
        md=result.geom["md"],
        y=result.geom["y"],
    )

    return metrics


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("  呼101井 居中度实测值 vs 模型假设值 对比")
    print("=" * 70)
    print(f"  模型域: {TOP_MD:.0f} – {BOTTOM_MD:.0f} m")
    print(f"  网格:   nz=250")
    print(f"  模式:   1D-2D 耦合（enable_gravity=True）")
    print("=" * 70)

    # 加载基线
    base_well, fluids, schedule, _ = load_hu101_tailpipe()
    print(f"\n基线井加载完成: {base_well.well_name}")
    print(f"  深度范围: {base_well.top_md_m} – {base_well.bottom_md_m} m")

    # 三套剖面统计
    print("\n三套居中度剖面对比（模型域内 200 点插值）:")
    print(f"  {'剖面':<25s}  {'均值':>6s}  {'最小':>6s}  {'最大':>6s}")
    print(f"  {'-'*25}  {'-'*6}  {'-'*6}  {'-'*6}")

    md_ref = np.linspace(TOP_MD, BOTTOM_MD, 200)
    for name, pts in [
        ("模型假设 (legacy)", ASSUMED_PROFILE),
        ("实测-扶正器之间", MEASURED_BETWEEN_CENT),
        ("实测-扶正器处", MEASURED_AT_CENT),
    ]:
        mds = [p[0] for p in pts]
        vals = [p[1] for p in pts]
        interp = np.interp(md_ref, mds, vals)
        print(f"  {name:<25s}  {np.mean(interp):6.3f}  {np.min(interp):6.3f}  {np.max(interp):6.3f}")

    # 构建三个 case
    cases = [
        ("baseline_assumed", "模型假设(legacy)", ASSUMED_PROFILE),
        ("measured_between_cent", "实测-扶正器之间", MEASURED_BETWEEN_CENT),
        ("measured_at_cent", "实测-扶正器处", MEASURED_AT_CENT),
    ]

    all_metrics = []
    for case_id, case_label, profile in cases:
        well = make_well_with_standoff(base_well, profile)
        m = run_case(case_label, well, fluids, schedule, nz=250)
        m["case_id"] = case_id
        all_metrics.append(m)

    # ============================================================
    # 汇总对比
    # ============================================================
    print("\n" + "=" * 70)
    print("  汇总对比")
    print("=" * 70)

    baseline = all_metrics[0]
    keys = ["effective_efficiency", "cement_occupation", "channeling_index", "mixing_index", "instability_index"]
    labels = {
        "effective_efficiency": "有效顶替效率",
        "cement_occupation": "水泥占据率",
        "channeling_index": "窜槽指数",
        "mixing_index": "混浆指数",
        "instability_index": "失稳指数",
    }

    header = f"  {'指标':<12s}  {'模型假设':>10s}  {'实测-扶正间':>12s}  {'实测-扶正处':>12s}  {'间-假设差':>10s}  {'处-假设差':>10s}"
    print(f"\n{header}")
    print(f"  {'-'*12}  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*10}  {'-'*10}")

    for key in keys:
        base_val = baseline[key]
        btw = all_metrics[1][key]
        atc = all_metrics[2][key]
        print(f"  {labels[key]:<12s}  {base_val:10.4f}  {btw:12.4f}  {atc:12.4f}  {btw-base_val:+10.4f}  {atc-base_val:+10.4f}")

    # 相对变化
    print(f"\n  相对变化率 (相对模型假设):")
    print(f"  {'指标':<12s}  {'实测-扶正间':>12s}  {'实测-扶正处':>12s}")
    print(f"  {'-'*12}  {'-'*12}  {'-'*12}")
    for key in keys:
        base_val = baseline[key]
        btw = all_metrics[1][key]
        atc = all_metrics[2][key]
        p_btw = (btw - base_val) / base_val * 100 if base_val != 0 else float("nan")
        p_atc = (atc - base_val) / base_val * 100 if base_val != 0 else float("nan")
        print(f"  {labels[key]:<12s}  {p_btw:+11.2f}%  {p_atc:+11.2f}%")

    # 保存 CSV
    csv_path = OUTPUT_DIR / "居中度对比汇总.csv"
    fieldnames = ["case_id", "case", "effective_efficiency", "cement_occupation",
                  "channeling_index", "mixing_index", "instability_index",
                  "mean_standoff", "elapsed_s"]
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for m in all_metrics:
            w.writerow({k: m.get(k) for k in fieldnames})

    # 保存 JSON
    summary = {
        "description": "呼101井居中度实测值 vs 模型假设值 敏感性对比",
        "model_domain_m": [TOP_MD, BOTTOM_MD],
        "nz": 250,
        "mode": "1D-2D coupled with gravity",
        "baseline_profile_note": "与 hu101_loader.py legacy 假设 0.38-0.48 完全一致",
        "measured_source_note": "从呼101井尾管居中度检测图读取（Pipe Standoff %）",
        "cases": all_metrics,
        "delta_between_vs_assumed": {
            k: round(all_metrics[1][k] - baseline[k], 6) for k in keys
        },
        "delta_at_vs_assumed": {
            k: round(all_metrics[2][k] - baseline[k], 6) for k in keys
        },
    }
    json_path = OUTPUT_DIR / "对比摘要.json"
    json_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"\n结果已保存至: {OUTPUT_DIR}")
    print(f"  - 居中度对比汇总.csv")
    print(f"  - 对比摘要.json")
    print(f"  - 各 case 子目录含深度剖面和2D场")


if __name__ == "__main__":
    main()
