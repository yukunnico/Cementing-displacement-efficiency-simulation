"""偏心度(居中度 standoff)敏感性参数扫描脚本。

对 HT1-004 井做 standoff → 有效顶替效率的定量扫描：
- 固定 R3 三闭包全开（enable_d2dga_auto_m=True, enable_d2dga_i3_flux=True, enable_true_buoyancy=True）
- 扫 standoff ∈ {0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30}
- 每次运行完整 1D-2D 流水线，提取全部指标
- 输出 CSV 到 results/偏心度敏感性扫描/敏感性扫描结果.csv
"""
from __future__ import annotations

import csv
import os
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from cemdisp.data.well_spec import DepthValuePoint, WellSpec
from cemdisp.data.loaders.ht1_004_loader import load_ht1_004_tailpipe
from cemdisp.runners.ht1_004_ablation import (
    AblationLevel,
    run_one_level,
    extract_ablation_metrics,
)

# R3: 三闭包全开
R3 = AblationLevel("R3", True, True, True)

# 扫描参数
STANDOFF_VALUES = [0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30]
NZ = 500          # 生产网格分辨率
DT = 4.0          # 时间步长 (s)
OUTPUT_CSV = Path(
    "D:/users/desktop/research/控压固井项目/cement model/"
    "results/偏心度敏感性扫描/敏感性扫描结果.csv"
)

CSV_COLUMNS = [
    "standoff", "eccentricity",
    "effective_efficiency", "channeling_index", "mixing_index",
    "instability_index", "cement_occupation", "elapsed_s",
]


def make_standoff_override(
    base_well: WellSpec, standoff: float
) -> WellSpec:
    """用恒定居中度覆盖井眼剖面，返回新的 WellSpec。"""
    profile = (
        DepthValuePoint(base_well.top_md_m, standoff),
        DepthValuePoint(base_well.bottom_md_m, standoff),
    )
    return replace(base_well, standoff_profile=profile)


def main() -> None:
    print("=" * 70)
    print("HT1-004 偏心度(居中度)敏感性参数扫描")
    print(f"三闭包: R3 (auto_m + I3 flux + true buoyancy)")
    print(f"网格: nz={NZ}, dt={DT}s")
    print(f"扫描点数: {len(STANDOFF_VALUES)}")
    print("=" * 70)

    # 加载基准数据
    base_well, fluids, schedule, _ = load_ht1_004_tailpipe()

    # 确保输出目录存在
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    # 写 CSV 表头
    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()

    rows = []

    for idx, standoff in enumerate(STANDOFF_VALUES):
        e_raw = 1.0 - standoff
        e_clipped = float(np.clip(e_raw, 0.05, 0.55))
        run_id = f"standoff_scan_{standoff:.2f}"

        print(f"\n--- [{idx + 1}/{len(STANDOFF_VALUES)}] "
              f"standoff={standoff:.2f}  e_raw={e_raw:.3f}  e_clipped={e_clipped:.3f} ---")

        well_override = make_standoff_override(base_well, standoff)
        t0 = time.perf_counter()

        try:
            result = run_one_level(
                R3,
                nz=NZ,
                dt=DT,
                well_spec_override=well_override,
                run_id=run_id,
                eccentricity=e_clipped,
            )
            elapsed = time.perf_counter() - t0
            metrics = extract_ablation_metrics(result)

            row = {
                "standoff": standoff,
                "eccentricity": e_clipped,
                "effective_efficiency": metrics["effective_efficiency"],
                "channeling_index": metrics["channeling_index"],
                "mixing_index": metrics["mixing_index"],
                "instability_index": metrics["instability_index"],
                "cement_occupation": metrics["cement_occupation"],
                "elapsed_s": round(elapsed, 1),
            }
            rows.append(row)
            print(f"  -> 效率={metrics['effective_efficiency']:.4f}  "
                  f"窜槽={metrics['channeling_index']:.4f}  "
                  f"混浆={metrics['mixing_index']:.4f}  "
                  f"失稳={metrics['instability_index']:.4f}  "
                  f"水泥占据={metrics['cement_occupation']:.4f}  "
                  f"耗时={elapsed:.1f}s")

        except Exception as exc:
            elapsed = time.perf_counter() - t0
            print(f"  !! 失败: {exc}")
            row = {
                "standoff": standoff,
                "eccentricity": e_clipped,
                "effective_efficiency": None,
                "channeling_index": None,
                "mixing_index": None,
                "instability_index": None,
                "cement_occupation": None,
                "elapsed_s": round(elapsed, 1),
            }
            rows.append(row)

        # 每完成一个 case 就追加写入 CSV（防中途崩溃丢数据）
        with open(OUTPUT_CSV, "a", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
            writer.writerow(row)

    # ---------- 汇总打印 ----------
    print("\n" + "=" * 70)
    print("扫描结果汇总")
    print("=" * 70)

    header = (
        f"{'standoff':>10s}  {'e':>6s}  {'效率':>8s}  "
        f"{'窜槽':>8s}  {'混浆':>8s}  {'失稳':>8s}  "
        f"{'水泥占据':>10s}  {'耗时(s)':>8s}"
    )
    print(header)
    print("-" * len(header))

    for r in rows:
        def _fmt(v):
            if v is None:
                return "    None"
            return f"{v:8.4f}"

        print(
            f"{r['standoff']:10.2f}  {r['eccentricity']:6.3f}  "
            f"{_fmt(r['effective_efficiency'])}  {_fmt(r['channeling_index'])}  "
            f"{_fmt(r['mixing_index'])}  {_fmt(r['instability_index'])}  "
            f"{_fmt(r['cement_occupation'])}  {r['elapsed_s']:8.1f}"
        )

    print(f"\n结果已保存至: {OUTPUT_CSV}")

    # ---------- 定量敏感性分析 ----------
    valid = [r for r in rows if r["effective_efficiency"] is not None]
    if len(valid) >= 2:
        effs = np.array([r["effective_efficiency"] for r in valid])
        es = np.array([r["eccentricity"] for r in valid])
        chans = np.array([r["channeling_index"] for r in valid])
        mixes = np.array([r["mixing_index"] for r in valid])
        instabs = np.array([r["instability_index"] for r in valid])

        # 线性区间斜率（e <= 0.50 部分，避免裁剪饱和）
        linear_mask = es <= 0.50
        if linear_mask.sum() >= 2:
            slope_eff, _ = np.polyfit(es[linear_mask], effs[linear_mask], 1)
            slope_chan, _ = np.polyfit(es[linear_mask], chans[linear_mask], 1)
            slope_mix, _ = np.polyfit(es[linear_mask], mixes[linear_mask], 1)
            slope_inst, _ = np.polyfit(es[linear_mask], instabs[linear_mask], 1)

            print(f"\n线性区间 (e <= 0.50) 斜率 (每单位偏心度变化):")
            print(f"  有效顶替效率: {slope_eff:+.4f} / e")
            print(f"  窜槽指数:     {slope_chan:+.4f} / e")
            print(f"  混浆指数:     {slope_mix:+.4f} / e")
            print(f"  失稳指数:     {slope_inst:+.4f} / e")

        # 近同心 vs 严重偏心
        near = rows[0]   # standoff=0.90
        far = rows[-1]   # standoff=0.30
        if near["effective_efficiency"] is not None and far["effective_efficiency"] is not None:
            d_eff = near["effective_efficiency"] - far["effective_efficiency"]
            print(f"\n近同心(standoff=0.90) vs 严重偏心(standoff=0.30):")
            print(f"  效率差: {d_eff:+.4f}")
            print(f"  窜槽差: {far['channeling_index'] - near['channeling_index']:+.4f}")
            print(f"  混浆差: {far['mixing_index'] - near['mixing_index']:+.4f}")
            print(f"  失稳差: {far['instability_index'] - near['instability_index']:+.4f}")


if __name__ == "__main__":
    main()