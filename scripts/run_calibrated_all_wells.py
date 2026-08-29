# DEPRECATED: 依赖 scripts/p3_p4_integration 的已删 metrics 列 cbl_eval_interval_efficiency
# （validation/ 于 2026-07-29 删除，运行即 KeyError），且标定参数 M=500/alpha=0.05 来自
# 已废弃的旧参数扫描。参见 scripts/cbl_window_comparison.py（B1 新建）
# 或 scripts/rerun_all_wells_corrected.py（论文/正式 8 井官方口径）替代方案。
"""
用标定参数运行所有井（Step 4.4）
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from scripts.p3_p4_integration import run_single_config, WELL_CONFIGS


def run_all_wells_with_calibrated_params():
    """使用标定后的参数运行所有井。"""
    # 基于 hu101 标定结果：M=500, alpha=0.05
    calibrated_M = 500.0
    calibrated_alpha = 0.05

    results = []
    for well_key in WELL_CONFIGS:
        print(f"\n=== [{well_key}] 标定参数运行 (M={calibrated_M}, alpha={calibrated_alpha}) ===")
        try:
            result = run_single_config(
                well_key,
                yield_regularization_M=calibrated_M,
                enable_axial_dispersion=True,
                dispersion_alpha=calibrated_alpha,
            )
            metrics = result["metrics"]
            cbl_measured = WELL_CONFIGS[well_key].get("cbl_measured")
            deviation = abs(metrics["cbl_eval_interval_efficiency"] - cbl_measured) if cbl_measured else None

            print(f"  CBL预测: {metrics['cbl_eval_interval_efficiency']:.4f}")
            if cbl_measured:
                print(f"  CBL实测: {cbl_measured:.4f}")
                print(f"  偏差: {deviation:.4f}")

            results.append({
                "well_key": well_key,
                "cbl_predicted": metrics["cbl_eval_interval_efficiency"],
                "cbl_measured": cbl_measured,
                "deviation": deviation,
                "effective_efficiency": metrics["effective_efficiency"],
                "channeling_index": metrics["channeling_index"],
                "mean_cement": metrics["mean_cement"],
                "mean_mud": metrics["mean_mud"],
            })
        except Exception as exc:
            print(f"  错误: {exc}")
            results.append({
                "well_key": well_key,
                "error": str(exc),
            })

    # 保存结果
    out_path = _PROJECT_ROOT / "results" / "p3_p4_integration" / "calibrated_all_wells.json"
    out_path.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n结果已保存到: {out_path}")

    # 打印汇总表
    print("\n=== 标定参数全井验证汇总 ===")
    print(f"{'井名':<10} {'预测CBL':<10} {'实测CBL':<10} {'偏差':<10} {'有效效率':<10} {'窜槽指数':<10}")
    print("-" * 70)
    for r in results:
        if "error" in r:
            print(f"{r['well_key']:<10} 错误: {r['error']}")
        else:
            pred = f"{r['cbl_predicted']:.4f}" if r['cbl_predicted'] is not None else "N/A"
            meas = f"{r['cbl_measured']:.4f}" if r['cbl_measured'] else "N/A"
            dev = f"{r['deviation']:.4f}" if r['deviation'] else "N/A"
            eff = f"{r['effective_efficiency']:.4f}" if r['effective_efficiency'] else "N/A"
            ch = f"{r['channeling_index']:.4f}" if r['channeling_index'] else "N/A"
            print(f"{r['well_key']:<10} {pred:<10} {meas:<10} {dev:<10} {eff:<10} {ch:<10}")


if __name__ == "__main__":
    run_all_wells_with_calibrated_params()
