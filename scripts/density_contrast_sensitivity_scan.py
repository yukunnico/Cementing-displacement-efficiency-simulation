"""密度差（水泥浆密度 vs 钻井液密度）敏感性参数扫描脚本。

对 HT1-004 井做 cement density -> 有效顶替效率的定量扫描：
- 固定 R3 三闭包全开（enable_d2dga_auto_m=True, enable_d2dga_i3_flux=True, enable_true_buoyancy=True）
- 固定泥浆密度 1900 kg/m³，扫描水泥浆密度（lead+tail 同时平移，保持 offset=30）
- 扫 cement_density ∈ {1600, 1750, 1900, 2000, 2100, 2200, 2350} kg/m³
- 每次运行完整 1D-2D 流水线，提取全部指标
- 输出 CSV 到 results/密度差敏感性扫描/密度差敏感性扫描结果.csv

注：本脚本非 D2DGA Tier 1 spec 范围，属 R3 模型辅助敏感性分析工具。
"""
from __future__ import annotations

import csv
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.loaders.ht1_004_loader import load_ht1_004_tailpipe
from cemdisp.runners.ht1_004_ablation import (
    AblationLevel,
    run_one_level,
    extract_ablation_metrics,
)

# R3: 三闭包全开
R3 = AblationLevel("R3", True, True, True)

# 扫描参数
# 泥浆密度固定为 HT1-004 默认值 1900 kg/m³ (1.90 g/cm³)
MUD_DENSITY_KG_M3 = 1900.0
# 水泥浆密度扫描点 (kg/m³) — 即 g/cm³ × 1000
CEMENT_DENSITY_VALUES = [1600, 1750, 1900, 2000, 2100, 2200, 2350]
# 原始 lead-tail offset = 1930 - 1900 = 30 kg/m³，保持此相对关系
LEAD_TAIL_OFFSET = 30

NZ = 500          # 生产网格分辨率
DT = 4.0          # 时间步长 (s)
OUTPUT_CSV = Path(
    "D:/users/desktop/research/控压固井项目/cement model/"
    "results/密度差敏感性扫描/密度差敏感性扫描结果.csv"
)

CSV_COLUMNS = [
    "cement_density_gcc", "lead_density_gcc", "mud_density_gcc",
    "density_contrast_gcc", "buoyancy_number",
    "effective_efficiency", "channeling_index", "mixing_index",
    "instability_index", "cement_occupation", "elapsed_s",
]


def make_density_fluids(
    base_fluids: tuple[FluidSpec, ...], cement_density: float
) -> tuple[FluidSpec, ...]:
    """修改 fluids tuple 中 LEAD 和 TAIL 的密度，保持 lead-tail offset=30 kg/m³。

    cement_density 为尾浆密度 (kg/m³)，领浆 = cement_density + 30。
    """
    new_fluids: list[FluidSpec] = []
    for fluid in base_fluids:
        if fluid.role == FluidRole.LEAD:
            new_fluids.append(replace(fluid, density_kg_m3=cement_density + LEAD_TAIL_OFFSET))
        elif fluid.role == FluidRole.TAIL:
            new_fluids.append(replace(fluid, density_kg_m3=cement_density))
        else:
            new_fluids.append(fluid)
    return tuple(new_fluids)


def compute_buoyancy_number_estimate(
    rho_disp: float, rho_mud: float, gap_m: float = 0.025,
    mu_mud: float = 0.053, velocity: float = 0.5,
) -> float:
    """估算浮力数 b = (Δρ)·g·(gap/2)² / (μ_mud · w₀)。

    使用典型值：gap=0.025m (25mm环空间隙), mu_mud=0.053 Pa·s, w₀=0.5 m/s。
    """
    g = 9.81
    d_half = max(gap_m / 2.0, 1.0e-6)
    delta_rho = rho_disp - rho_mud
    denom = max(mu_mud * max(velocity, 1.0e-6), 1.0e-9)
    return delta_rho * g * d_half ** 2 / denom


def main() -> None:
    print("=" * 70)
    print("HT1-004 密度差（水泥浆 vs 钻井液）敏感性参数扫描")
    print(f"三闭包: R3 (auto_m + I3 flux + true buoyancy)")
    print(f"网格: nz={NZ}, dt={DT}s")
    print(f"泥浆密度固定: {MUD_DENSITY_KG_M3} kg/m3 ({MUD_DENSITY_KG_M3/1000:.2f} g/cm3)")
    print(f"扫描点数: {len(CEMENT_DENSITY_VALUES)}")
    print(f"水泥密度范围: {CEMENT_DENSITY_VALUES[0]/1000:.2f} ~ {CEMENT_DENSITY_VALUES[-1]/1000:.2f} g/cm3")
    print("=" * 70)

    # 加载基准数据
    base_well, base_fluids, _, _ = load_ht1_004_tailpipe()

    # 确保输出目录存在
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    # 写 CSV 表头
    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()

    rows = []

    for idx, cement_density in enumerate(CEMENT_DENSITY_VALUES):
        cement_density_gcc = cement_density / 1000.0
        lead_density = cement_density + LEAD_TAIL_OFFSET
        lead_density_gcc = lead_density / 1000.0
        density_contrast = cement_density_gcc - MUD_DENSITY_KG_M3 / 1000.0

        # 估算浮力数
        b_est = compute_buoyancy_number_estimate(cement_density, MUD_DENSITY_KG_M3)

        run_id = f"density_scan_{cement_density_gcc:.2f}"

        print(f"\n--- [{idx + 1}/{len(CEMENT_DENSITY_VALUES)}] "
              f"cement={cement_density_gcc:.2f} g/cm3  "
              f"lead={lead_density_gcc:.2f}  "
              f"drho={density_contrast:+.2f} g/cm3  "
              f"b_est≈{b_est:+.3f} ---")

        # 修改 fluids 中的水泥密度
        fluids_override = make_density_fluids(base_fluids, cement_density)

        t0 = time.perf_counter()

        try:
            # run_one_level 内部会重新 load_ht1_004_tailpipe()，
            # 因此我们需要直接构建流水线，手动传入 fluids_override
            from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
            from cemdisp.transport1d import CasingFlowSolver
            from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
            from cemdisp.runners.ht1_004_tailpipe import annulus_stop_time_s

            loaded_well, _, schedule, _ = load_ht1_004_tailpipe()
            well_spec = base_well  # 使用默认 well spec

            # 1D casing flow
            casing_solver = CasingFlowSolver(enable_gravity=True)
            casing_result = casing_solver.run(well_spec, fluids_override, schedule)

            # Coupled inlet provider
            provider = build_coupled_annulus_inlet_provider(
                casing_result, casing_solver, fluids_override
            )

            # Annulus stop time
            total_t = annulus_stop_time_s(casing_result=casing_result, fluids=fluids_override)

            # 2D D2DGA solver
            solver = AnnulusD2DGASolver(
                dt=DT,
                nz=NZ,
                ny=40,
                total_t=total_t,
                enable_d2dga=True,
                enable_d2dga_auto_m=R3.enable_d2dga_auto_m,
                enable_d2dga_i3_flux=R3.enable_d2dga_i3_flux,
                enable_true_buoyancy=R3.enable_true_buoyancy,
                open_outlet=True,
            )

            result = solver.run(well_spec, fluids_override, provider)
            elapsed = time.perf_counter() - t0
            metrics = extract_ablation_metrics(result)

            # 从结果中提取实际浮力数
            b_actual = result.summary.get("buoyancy_number") if isinstance(result.summary, dict) else None

            row = {
                "cement_density_gcc": cement_density_gcc,
                "lead_density_gcc": lead_density_gcc,
                "mud_density_gcc": MUD_DENSITY_KG_M3 / 1000.0,
                "density_contrast_gcc": density_contrast,
                "buoyancy_number": round(b_actual, 4) if b_actual is not None else None,
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
                  f"浮力数b={b_actual:.4f}" if b_actual is not None else "浮力数b=N/A",
                  f" 耗时={elapsed:.1f}s")

        except Exception as exc:
            elapsed = time.perf_counter() - t0
            print(f"  !! 失败: {exc}")
            import traceback
            traceback.print_exc()
            row = {
                "cement_density_gcc": cement_density_gcc,
                "lead_density_gcc": lead_density_gcc,
                "mud_density_gcc": MUD_DENSITY_KG_M3 / 1000.0,
                "density_contrast_gcc": density_contrast,
                "buoyancy_number": None,
                "effective_efficiency": None,
                "channeling_index": None,
                "mixing_index": None,
                "instability_index": None,
                "cement_occupation": None,
                "elapsed_s": round(elapsed, 1),
            }
            rows.append(row)

        # 每完成一个 case 就追加写入 CSV
        with open(OUTPUT_CSV, "a", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
            writer.writerow(row)

    # ---------- 汇总打印 ----------
    print("\n" + "=" * 70)
    print("扫描结果汇总")
    print("=" * 70)

    header = (
        f"{'水泥密度':>10s}  {'泥浆密度':>10s}  {'Δρ':>8s}  {'浮力数b':>10s}  "
        f"{'效率':>8s}  {'窜槽':>8s}  {'混浆':>8s}  {'失稳':>8s}  "
        f"{'水泥占据':>10s}  {'耗时(s)':>8s}"
    )
    print(header)
    print("-" * len(header))

    for r in rows:
        def _fmt(v):
            if v is None:
                return "     None"
            return f"{v:8.4f}"

        print(
            f"{r['cement_density_gcc']:10.2f}  {r['mud_density_gcc']:10.2f}  "
            f"{r['density_contrast_gcc']:8.2f}  {_fmt(r['buoyancy_number'])}  "
            f"{_fmt(r['effective_efficiency'])}  {_fmt(r['channeling_index'])}  "
            f"{_fmt(r['mixing_index'])}  {_fmt(r['instability_index'])}  "
            f"{_fmt(r['cement_occupation'])}  {r['elapsed_s']:8.1f}"
        )

    print(f"\n结果已保存至: {OUTPUT_CSV}")

    # ---------- 定量敏感性分析 ----------
    valid = [r for r in rows if r["effective_efficiency"] is not None]
    if len(valid) >= 2:
        dc_arr = np.array([r["density_contrast_gcc"] for r in valid])
        effs = np.array([r["effective_efficiency"] for r in valid])
        chans = np.array([r["channeling_index"] for r in valid])
        mixes = np.array([r["mixing_index"] for r in valid])
        instabs = np.array([r["instability_index"] for r in valid])
        occs = np.array([r["cement_occupation"] for r in valid])

        # 全区间线性斜率
        slope_eff, intercept_eff = np.polyfit(dc_arr, effs, 1)
        slope_chan, _ = np.polyfit(dc_arr, chans, 1)
        slope_mix, _ = np.polyfit(dc_arr, mixes, 1)
        slope_inst, _ = np.polyfit(dc_arr, instabs, 1)
        slope_occ, _ = np.polyfit(dc_arr, occs, 1)

        print(f"\n全区间线性斜率 (每 0.1 g/cm3 密度差变化):")
        print(f"  有效顶替效率: {slope_eff * 0.1:+.4f} / 0.1 g/cm3")
        print(f"  窜槽指数:     {slope_chan * 0.1:+.4f} / 0.1 g/cm3")
        print(f"  混浆指数:     {slope_mix * 0.1:+.4f} / 0.1 g/cm3")
        print(f"  失稳指数:     {slope_inst * 0.1:+.4f} / 0.1 g/cm3")
        print(f"  水泥占据率:   {slope_occ * 0.1:+.4f} / 0.1 g/cm3")

        # 最低 vs 最高密度差
        lo = valid[0]   # 最低密度 (1.60 g/cm³, Δρ = -0.30)
        hi = valid[-1]  # 最高密度 (2.35 g/cm³, Δρ = +0.45)
        if lo["effective_efficiency"] is not None and hi["effective_efficiency"] is not None:
            d_eff = hi["effective_efficiency"] - lo["effective_efficiency"]
            d_chan = hi["channeling_index"] - lo["channeling_index"]
            d_mix = hi["mixing_index"] - lo["mixing_index"]
            d_inst = hi["instability_index"] - lo["instability_index"]
            print(f"\n最低密度差(Δρ={lo['density_contrast_gcc']:+.2f}) vs 最高密度差(Δρ={hi['density_contrast_gcc']:+.2f}):")
            print(f"  效率差: {d_eff:+.4f}")
            print(f"  窜槽差: {d_chan:+.4f}")
            print(f"  混浆差: {d_mix:+.4f}")
            print(f"  失稳差: {d_inst:+.4f}")

        # 敏感性排名
        sensitivities = {
            "有效顶替效率": abs(slope_eff),
            "窜槽指数": abs(slope_chan),
            "混浆指数": abs(slope_mix),
            "失稳指数": abs(slope_inst),
            "水泥占据率": abs(slope_occ),
        }
        ranked = sorted(sensitivities.items(), key=lambda x: x[1], reverse=True)
        print(f"\n指标对密度差的敏感性排名 (线性斜率绝对值):")
        for name, sens in ranked:
            print(f"  {name}: {sens:.4f}")


if __name__ == "__main__":
    main()