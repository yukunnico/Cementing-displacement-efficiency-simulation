"""三闭包 R0->R3 逐级贡献分解扫描脚本。
额外包含 NONE 级别（enable_d2dga=False 纯平流基线），共 5 级。
产出 results/三闭包贡献分解/三闭包贡献分解结果.csv

注：本脚本非 D2DGA Tier 1 spec 范围，属 R3 模型辅助敏感性分析工具。
"""
from __future__ import annotations

import csv
import os
import time
from pathlib import Path
from typing import Dict

from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver, AnnulusSimulationResult
from cemdisp.transport1d import CasingFlowSolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.runners.ht1_004_tailpipe import annulus_stop_time_s
from cemdisp.data.loaders.ht1_004_loader import load_ht1_004_tailpipe

# ---------------------------------------------------------------------------
# 级别定义：NONE / R0 / R1 / R2 / R3
# ---------------------------------------------------------------------------
LEVELS = [
    # (name, enable_d2dga, enable_d2dga_auto_m, enable_d2dga_i3_flux, enable_true_buoyancy)
    ("NONE", False, False, False, False),
    ("R0",   True,  False, False, False),
    ("R1",   True,  True,  False, False),
    ("R2",   True,  True,  True,  False),
    ("R3",   True,  True,  True,  True),
]

OUTPUT_DIR = Path("results/三闭包贡献分解")
CSV_PATH = OUTPUT_DIR / "三闭包贡献分解结果.csv"
CSV_COLUMNS = [
    "level", "enable_d2dga", "auto_m", "i3_flux", "true_buoyancy",
    "effective_efficiency", "channeling_index", "mixing_index",
    "instability_index", "cement_occupation", "elapsed_s",
]


def extract_metrics(result: AnnulusSimulationResult) -> dict:
    """从 AnnulusSimulationResult 提取扁平化指标字典。"""
    final = result.summary.get("最终结果", {}) if isinstance(result.summary, dict) else {}
    return {
        "effective_efficiency": final.get("全井段最终有效顶替效率"),
        "cement_occupation": final.get("最终水泥浆占据率"),
        "channeling_index": final.get("最终窜槽指数"),
        "mixing_index": final.get("最终混浆指数"),
        "instability_index": final.get("最终失稳指数"),
    }


def run_one_level(
    name: str,
    enable_d2dga: bool,
    enable_d2dga_auto_m: bool,
    enable_d2dga_i3_flux: bool,
    enable_true_buoyancy: bool,
    *,
    nz: int = 500,
    dt: float = 4.0,
) -> tuple[dict, float]:
    """运行一个闭包级别，返回 (指标字典, 耗时秒)。"""
    loaded_well, fluids, schedule, _ = load_ht1_004_tailpipe()
    well_spec = loaded_well

    casing_solver = CasingFlowSolver(enable_gravity=True)
    casing_result = casing_solver.run(well_spec, fluids, schedule)
    provider = build_coupled_annulus_inlet_provider(casing_result, casing_solver, fluids)
    total_t = annulus_stop_time_s(casing_result=casing_result, fluids=fluids)

    solver = AnnulusD2DGASolver(
        dt=dt, nz=nz, ny=40, total_t=total_t,
        enable_d2dga=enable_d2dga,
        enable_d2dga_auto_m=enable_d2dga_auto_m,
        enable_d2dga_i3_flux=enable_d2dga_i3_flux,
        enable_true_buoyancy=enable_true_buoyancy,
        open_outlet=True,
    )

    t0 = time.perf_counter()
    result = solver.run(well_spec, fluids, provider)
    elapsed = time.perf_counter() - t0

    metrics = extract_metrics(result)
    return metrics, elapsed


def save_csv(rows: list[dict], csv_path: Path) -> None:
    """保存结果到 CSV。"""
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def print_results_table(rows: list[dict]) -> None:
    """打印结果表。"""
    header = f"{'Level':<6} {'enable_d2dga':>12} {'auto_m':>7} {'i3_flux':>7} {'true_buoy':>9} {'效率':>8} {'窜槽':>8} {'混浆':>8} {'失稳':>8} {'水泥占据':>8} {'耗时s':>8}"
    sep = "-" * len(header)
    print("\n=== 三闭包贡献分解：全级别结果 ===")
    print(sep)
    print(header)
    print(sep)
    for r in rows:
        print(
            f"{r['level']:<6} {str(r['enable_d2dga']):>12} {str(r['auto_m']):>7} {str(r['i3_flux']):>7} {str(r['true_buoyancy']):>9} "
            f"{r['effective_efficiency']:>8.4f} {r['channeling_index']:>8.4f} {r['mixing_index']:>8.4f} "
            f"{r['instability_index']:>8.4f} {r['cement_occupation']:>8.4f} {r['elapsed_s']:>8.1f}"
        )
    print(sep)


def print_marginal_table(rows: list[dict]) -> None:
    """打印边际贡献表。"""
    # 按 level 名索引
    by_level = {r["level"]: r for r in rows}

    def _delta(hi: str, lo: str, key: str) -> float | None:
        h = by_level.get(hi)
        l = by_level.get(lo)
        if h and l and key in h and key in l:
            return h[key] - l[key]
        return None

    metrics = ["effective_efficiency", "channeling_index", "mixing_index", "instability_index", "cement_occupation"]
    labels = {
        "effective_efficiency": "最终效率",
        "channeling_index": "窜槽指数",
        "mixing_index": "混浆指数",
        "instability_index": "失稳指数",
        "cement_occupation": "水泥占据率",
    }

    print("\n=== 边际贡献表（增量） ===")
    header = f"{'贡献项':<20}" + "".join(f"{labels[m]:>10}" for m in metrics)
    sep = "-" * len(header)
    print(sep)
    print(header)
    print(sep)

    contributions = [
        ("D2DGA整体 (R0-NONE)", "R0", "NONE"),
        ("R1贡献 (R1-R0)", "R1", "R0"),
        ("R2贡献 (R2-R1)", "R2", "R1"),
        ("R3贡献 (R3-R2)", "R3", "R2"),
        ("R3相对R0 (R3-R0)", "R3", "R0"),
    ]

    for label, hi, lo in contributions:
        vals = [_delta(hi, lo, m) for m in metrics]
        val_str = "".join(f"{v:>10.4f}" if v is not None else f"{'N/A':>10}" for v in vals)
        print(f"{label:<20}{val_str}")

    print(sep)


def main() -> None:
    nz, dt = 500, 4.0
    print(f"=== 三闭包贡献分解扫描 (nz={nz}, dt={dt}) ===")
    print(f"共 {len(LEVELS)} 个级别: {', '.join(l[0] for l in LEVELS)}")

    rows: list[dict] = []
    for name, en_d2dga, auto_m, i3_flux, true_buoy in LEVELS:
        print(f"\n--- 运行 {name} (enable_d2dga={en_d2dga}, auto_m={auto_m}, i3_flux={i3_flux}, true_buoyancy={true_buoy}) ---")
        metrics, elapsed = run_one_level(
            name, en_d2dga, auto_m, i3_flux, true_buoy,
            nz=nz, dt=dt,
        )
        row = {
            "level": name,
            "enable_d2dga": en_d2dga,
            "auto_m": auto_m,
            "i3_flux": i3_flux,
            "true_buoyancy": true_buoy,
            **metrics,
            "elapsed_s": round(elapsed, 1),
        }
        rows.append(row)
        print(f"  -> 效率={metrics['effective_efficiency']:.4f}, 耗时={elapsed:.1f}s")

    # 保存 CSV
    save_csv(rows, CSV_PATH)
    print(f"\nCSV 已保存: {CSV_PATH}")

    # 打印结果表
    print_results_table(rows)

    # 打印边际贡献表
    print_marginal_table(rows)


if __name__ == "__main__":
    main()