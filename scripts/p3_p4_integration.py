"""
P3 集成验证 + P4 参数标定脚本

功能：
1. 对每口井运行基准配置（M=0, dispersion=False）和改进配置（M=100, dispersion=True）
2. 提取关键指标进行对比
3. 以 hu102 为基准扫描参数空间，找最优 (M, alpha) 组合
4. 用 hu103 做交叉验证

使用方式：
    python scripts/p3_p4_integration.py --mode compare --wells hu101 hu102
    python scripts/p3_p4_integration.py --mode scan --wells hu102
    python scripts/p3_p4_integration.py --mode validate --wells hu103
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

# 确保项目根目录在 PYTHONPATH 中
_PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from cemdisp.data.loaders import (
    load_hu101_tailpipe,
    load_hu102_tailpipe,
    load_hu103_tailpipe,
    load_hu1_tailpipe,
    load_hu2_tailpipe,
    load_ht1_001_tailpipe,
)
from cemdisp.data.validation_data import ValidationData
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import AnnulusInletState, build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver
from cemdisp.transport1d.casing_flow import CasingFlowResult


# 井配置：加载函数、已知 CBL 实测值（如有）
WELL_CONFIGS: dict[str, dict[str, object]] = {
    "hu101": {
        "loader": load_hu101_tailpipe,
        "cbl_measured": 0.6277,  # 62.77%
        "nz": 500,
    },
    "hu102": {
        "loader": load_hu102_tailpipe,
        "cbl_measured": 0.6665,  # 66.65%
        "nz": 500,
    },
    "hu103": {
        "loader": load_hu103_tailpipe,
        "cbl_measured": 0.1206,  # 12.06%
        "nz": 500,
    },
    "hu1": {
        "loader": load_hu1_tailpipe,
        "cbl_measured": None,
        "nz": 500,
    },
    "hu2": {
        "loader": load_hu2_tailpipe,
        "cbl_measured": None,  # HT1-002 定性参考
        "nz": 500,
    },
    "ht1_001": {
        "loader": load_ht1_001_tailpipe,
        "cbl_measured": None,
        "nz": 500,
    },
}


def run_single_config(
    well_key: str,
    *,
    yield_regularization_M: float = 0.0,
    enable_axial_dispersion: bool = False,
    dispersion_alpha: float = 0.2,
) -> dict[str, object]:
    """对单井运行一次配置，返回关键指标和元数据。

    参数：
        yield_regularization_M: 屈服正则化参数，0 表示关闭（基准）
        enable_axial_dispersion: 是否启用管内轴向弥散
        dispersion_alpha: 无量纲弥散系数

    返回：
        dict 包含关键指标、运行时间和配置参数
    """
    config = WELL_CONFIGS[well_key]
    loader: Callable = config["loader"]  # type: ignore[assignment]
    nz: int = config["nz"]  # type: ignore[assignment]

    well_spec, fluids, schedule, validation_data = loader()

    t0 = time.perf_counter()

    # 1D 套管内流动
    casing_solver = CasingFlowSolver(
        enable_gravity=True,
        enable_axial_dispersion=enable_axial_dispersion,
        dispersion_alpha=dispersion_alpha,
    )
    casing_result = casing_solver.run(well_spec, fluids, schedule)

    # 鞋口出流 → 环空入口桥接
    coupled_provider = build_coupled_annulus_inlet_provider(
        casing_result,
        casing_solver,
        fluids,
        split_cement_phases=True,
    )

    # 环空二维顶替停止时间
    total_t_s = float(casing_result.cement_end_time_s)

    # 2D 环空求解
    annulus_solver = AnnulusD2DGASolver(
        total_t=total_t_s,
        nz=nz,
        yield_regularization_M=yield_regularization_M,
    )
    result = annulus_solver.run(well_spec, fluids, coupled_provider)

    elapsed_s = time.perf_counter() - t0

    # 提取最终关键指标
    final = result.metrics.iloc[-1]
    metrics = {
        "effective_efficiency": float(final["effective_efficiency"]),
        "cbl_eval_interval_efficiency": float(final["cbl_eval_interval_efficiency"]),
        "channeling_index": float(final["channeling_index"]),
        "mixing_index": float(final["mixing_index"]),
        "instability_index": float(final["instability_index"]),
        "mean_cement": float(final["mean_cement"]),
        "mean_mud": float(final["mean_mud"]),
    }

    # 提取窄边水泥浓度剖面（最终时刻）
    depth_profiles = result.depth_profiles
    narrow_cement = None
    wide_cement = None
    if "窄边水泥浓度" in depth_profiles.columns:
        narrow_cement = depth_profiles["窄边水泥浓度"].to_list()
    if "宽边水泥浓度" in depth_profiles.columns:
        wide_cement = depth_profiles["宽边水泥浓度"].to_list()

    return {
        "well_key": well_key,
        "config": {
            "yield_regularization_M": yield_regularization_M,
            "enable_axial_dispersion": enable_axial_dispersion,
            "dispersion_alpha": dispersion_alpha,
        },
        "metrics": metrics,
        "elapsed_s": elapsed_s,
        "narrow_cement_profile": narrow_cement,
        "wide_cement_profile": wide_cement,
        "depth_m": depth_profiles["井深_m"].to_list() if "井深_m" in depth_profiles.columns else None,
    }


def compare_baseline_vs_improved(well_key: str) -> dict[str, object]:
    """对比单井的基准配置和改进配置。

    基准：M=0, dispersion=False（等同于旧模型行为）
    改进：M=100, dispersion=True, alpha=0.2
    """
    print(f"\n=== [{well_key}] 基准配置 ===")
    baseline = run_single_config(
        well_key,
        yield_regularization_M=0.0,
        enable_axial_dispersion=False,
        dispersion_alpha=0.2,
    )
    print(f"  耗时: {baseline['elapsed_s']:.1f}s")
    print(f"  CBL评价井段效率: {baseline['metrics']['cbl_eval_interval_efficiency']:.4f}")
    print(f"  窜槽指数: {baseline['metrics']['channeling_index']:.4f}")

    print(f"\n=== [{well_key}] 改进配置 (M=100, dispersion=True) ===")
    improved = run_single_config(
        well_key,
        yield_regularization_M=100.0,
        enable_axial_dispersion=True,
        dispersion_alpha=0.2,
    )
    print(f"  耗时: {improved['elapsed_s']:.1f}s")
    print(f"  CBL评价井段效率: {improved['metrics']['cbl_eval_interval_efficiency']:.4f}")
    print(f"  窜槽指数: {improved['metrics']['channeling_index']:.4f}")

    # 计算变化量
    delta = {
        key: improved["metrics"][key] - baseline["metrics"][key]
        for key in baseline["metrics"]
    }

    print(f"\n=== [{well_key}] 变化量 (改进 - 基准) ===")
    for key, value in delta.items():
        print(f"  {key}: {value:+.4f}")

    cbl_measured = WELL_CONFIGS[well_key].get("cbl_measured")
    if cbl_measured is not None:
        baseline_deviation = baseline["metrics"]["cbl_eval_interval_efficiency"] - cbl_measured
        improved_deviation = improved["metrics"]["cbl_eval_interval_efficiency"] - cbl_measured
        print(f"\n  实测 CBL: {cbl_measured:.4f}")
        print(f"  基准偏差: {baseline_deviation:+.4f}")
        print(f"  改进偏差: {improved_deviation:+.4f}")
        print(f"  偏差改善: {abs(improved_deviation) - abs(baseline_deviation):+.4f}")

    return {
        "well_key": well_key,
        "baseline": baseline,
        "improved": improved,
        "delta": delta,
        "cbl_measured": cbl_measured,
    }


def scan_parameters(
    well_key: str,
    *,
    M_values: tuple[float, ...] = (0.0, 50.0, 100.0, 200.0, 300.0, 500.0),
    alpha_values: tuple[float, ...] = (0.05, 0.1, 0.2, 0.3, 0.5),
) -> dict[str, object]:
    """扫描参数空间，找使预测 CBL 与实测最接近的组合。

    仅适用于已知 CBL 实测值的井。
    """
    cbl_measured = WELL_CONFIGS[well_key].get("cbl_measured")
    if cbl_measured is None:
        raise ValueError(f"{well_key} 没有已知 CBL 实测值，无法做参数扫描")

    results: list[dict[str, object]] = []
    best_deviation = float("inf")
    best_config = None

    total_runs = len(M_values) * len(alpha_values)
    run_idx = 0

    print(f"\n=== [{well_key}] 参数扫描: {total_runs} 组配置 ===")
    for M in M_values:
        for alpha in alpha_values:
            run_idx += 1
            print(f"  [{run_idx}/{total_runs}] M={M}, alpha={alpha} ... ", end="", flush=True)
            t0 = time.perf_counter()
            result = run_single_config(
                well_key,
                yield_regularization_M=M,
                enable_axial_dispersion=(alpha > 0.0),
                dispersion_alpha=alpha,
            )
            elapsed = time.perf_counter() - t0
            cbl_pred = result["metrics"]["cbl_eval_interval_efficiency"]
            deviation = abs(cbl_pred - cbl_measured)
            print(f"CBL_pred={cbl_pred:.4f}, dev={deviation:.4f}, {elapsed:.1f}s")

            row = {
                "M": M,
                "alpha": alpha,
                **result["metrics"],
                "deviation": deviation,
                "elapsed_s": elapsed,
            }
            results.append(row)

            if deviation < best_deviation:
                best_deviation = deviation
                best_config = row.copy()

    print(f"\n=== 最优参数组合 ===")
    if best_config is not None:
        print(f"  M={best_config['M']}, alpha={best_config['alpha']}")
        print(f"  预测 CBL={best_config['cbl_eval_interval_efficiency']:.4f}")
        print(f"  实测 CBL={cbl_measured:.4f}")
        print(f"  绝对偏差={best_config['deviation']:.4f}")

    return {
        "well_key": well_key,
        "cbl_measured": cbl_measured,
        "best_config": best_config,
        "all_results": results,
    }


def export_comparison(results: list[dict[str, object]], output_dir: Path) -> None:
    """导出对比结果到 JSON 和 CSV。"""
    output_dir.mkdir(parents=True, exist_ok=True)

    # JSON 汇总
    json_path = output_dir / "p3_p4_integration_results.json"
    json_path.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")

    # CSV 指标对比
    rows: list[dict[str, object]] = []
    for well_result in results:
        well_key = well_result["well_key"]
        baseline = well_result["baseline"]["metrics"]
        improved = well_result["improved"]["metrics"]
        delta = well_result["delta"]
        cbl_measured = well_result.get("cbl_measured")

        for metric_key in baseline:
            rows.append({
                "井名": well_key,
                "指标": metric_key,
                "基准值": baseline[metric_key],
                "改进值": improved[metric_key],
                "变化量": delta[metric_key],
                "实测CBL": cbl_measured,
            })

    df = pd.DataFrame(rows)
    csv_path = output_dir / "p3_p4_comparison.csv"
    df.to_csv(csv_path, index=False, encoding="utf-8-sig")
    print(f"\n结果已导出到: {output_dir}")


def export_scan_results(scan_result: dict[str, object], output_dir: Path) -> None:
    """导出参数扫描结果。"""
    output_dir.mkdir(parents=True, exist_ok=True)

    well_key = scan_result["well_key"]

    # JSON
    json_path = output_dir / f"{well_key}_parameter_scan.json"
    json_path.write_text(json.dumps(scan_result, ensure_ascii=False, indent=2), encoding="utf-8")

    # CSV
    df = pd.DataFrame(scan_result["all_results"])
    csv_path = output_dir / f"{well_key}_parameter_scan.csv"
    df.to_csv(csv_path, index=False, encoding="utf-8-sig")
    print(f"扫描结果已导出到: {csv_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="P3 集成验证 + P4 参数标定")
    parser.add_argument("--mode", choices=["compare", "scan", "validate"], default="compare",
                        help="运行模式: compare=基准vs改进对比, scan=参数扫描, validate=交叉验证")
    parser.add_argument("--wells", nargs="+", default=["hu101", "hu102"],
                        help="要运行的井名列表")
    parser.add_argument("--output-dir", type=Path, default=_PROJECT_ROOT / "results" / "p3_p4_integration",
                        help="结果输出目录")
    args = parser.parse_args()

    if args.mode == "compare":
        all_results: list[dict[str, object]] = []
        for well_key in args.wells:
            if well_key not in WELL_CONFIGS:
                print(f"警告: 未知井名 {well_key}，跳过")
                continue
            result = compare_baseline_vs_improved(well_key)
            all_results.append(result)
        export_comparison(all_results, args.output_dir)

    elif args.mode == "scan":
        for well_key in args.wells:
            if well_key not in WELL_CONFIGS:
                print(f"警告: 未知井名 {well_key}，跳过")
                continue
            scan_result = scan_parameters(well_key)
            export_scan_results(scan_result, args.output_dir)

    elif args.mode == "validate":
        # 交叉验证：使用已标定的最优参数在未知井上运行
        # TODO: 从 scan 结果中读取最优参数
        print("交叉验证模式: 需要先从扫描结果中获取最优参数")


if __name__ == "__main__":
    main()
