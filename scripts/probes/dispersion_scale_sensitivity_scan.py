"""M1 弥散系数 κ 缩放敏感性扫描脚本（CFL on/off 双组存档）。

任务背景（.superpowers/sdd/2026-08-23-annulus-distortion-fix/task-5-brief.md）：
Task 4 (M1) 已将弥散系数改为按 dt 归一：_dt_norm = dispersion_dt_scale * (dt_step / dispersion_dt_ref)。
本脚本扫描 dispersion_dt_scale ∈ {0.0, 0.25, 0.5, 1.0}，CFL on/off 各一组：

- CFL off：固定 dt_step = 4.0 = dispersion_dt_ref，_dt_norm = scale（κ 归一在固定步长下直接生效）
- CFL on ：dt_step 自适应（远小于 4.0），_dt_norm = scale * dt_step/4.0（κ 按实际步长归一）

验收目标（设计稿 §4）：
- mixing 从 0.59 降到 0.2-0.35；
- scale→0 时 mixing 主要由平流数值扩散主导；
- 前沿长度回到米级（2-5m）。

对 HT1-004 井，R3 三闭包全开（enable_d2dga_auto_m / i3_flux / true_buoyancy）。
因 run_one_level 不透传弥散参数，采用 density_contrast_sensitivity_scan.py 的手动流水线：
1D casing → build_coupled_annulus_inlet_provider → annulus_stop_time_s → AnnulusD2DGASolver → run → 指标 → CSV。

每 case 即时追加写 CSV（防崩）。支持冒烟模式（环境变量 M1_SMOKE=1）：
仅 NZ=30, NY=10, scale∈{0.0,1.0} × cfl_on={True,False}，用于端到端自检。
"""
from __future__ import annotations

import csv
import os
import time
from pathlib import Path

import numpy as np

from cemdisp.data.loaders.ht1_004_loader import load_ht1_004_tailpipe
from cemdisp.runners.ht1_004_ablation import extract_ablation_metrics

# ---------------------------------------------------------------------------
# 常量（全量扫描口径，brief Step 1；冒烟模式在 __main__ 里覆盖为小网格）
# ---------------------------------------------------------------------------
SCALE_VALUES = [0.0, 0.25, 0.5, 1.0]
NZ = 500          # 生产网格分辨率（轴向）
NY = 40           # 周向分辨率
DT = 4.0          # 名义时间步长 (s)，等于 dispersion_dt_ref
ENABLE_CFL_ADAPTIVE = True  # 归一化只在 dt_step≠4.0 时真正生效（CFL on 组）

# CSV 输出目录与列
OUTPUT_DIR = Path(
    "D:/users/desktop/research/控压固井项目/cement model/results/m1_dispersion_scale"
)
OUTPUT_CSV = OUTPUT_DIR / "m1_dispersion_scale_scan.csv"

CSV_COLUMNS = [
    "scale", "cfl_on",
    "effective_efficiency", "channeling_index", "mixing_index",
    "instability_index", "cement_occupation", "front_length_m", "elapsed_s",
]

# 前缘长度阈值（等值线位置）
FRONT_C_LOW = 0.02    # 前缘包络浓度阈值
FRONT_C_HIGH = 0.98   # 尾缘包络浓度阈值


def _contour_md(profile: np.ndarray, md: np.ndarray, level: float) -> float:
    """在剖面 profile(md) 上求浓度等值线 profile = level 的测深 [m]（线性插值）。

    口径与 cemdisp/diagnostics/flow_classification._contour_position 一致：
    水泥浆自鞋口（md 大端）向上顶替，profile 总体随 md 递减（1→0）。
    对存在多个交叉的非单调剖面（如弥散尖峰），取**最下游（md 最大）**交叉点，
    对应前缘/尾缘包络，与论文"所有波速 w_f 的并集"定义一致。

    Args:
        profile: 截面平均浓度剖面 c̄(md)，形状 (nz,)
        md: 测深 [m]，形状 (nz,)，单调递增
        level: 浓度阈值，取值 (0, 1)

    Returns:
        等值线测深 [m]；若剖面整体位于 level 同侧（无交叉）则返回 np.nan
    """
    above = profile >= level
    transitions = np.flatnonzero(above[:-1] != above[1:])
    if transitions.size == 0:
        return float("nan")
    j = int(transitions[-1])  # 最下游交叉点
    p0, p1 = float(profile[j]), float(profile[j + 1])
    s0, s1 = float(md[j]), float(md[j + 1])
    dp = p1 - p0
    if abs(dp) < 1.0e-12:  # 零除保护：平台期取左端点
        return s0
    return s0 + (level - p0) * (s1 - s0) / dp


def compute_front_length_m(profiles) -> float:
    """从深度剖面计算水泥浆前沿长度 [m] = MD(c=0.02) 与 MD(c=0.98) 的测深差。

    profiles: result.depth_profiles (pd.DataFrame)，含 "井深_m" 与 "水泥平均浓度" 列。
    任一等值线缺失（np.nan）则返回 nan（不虚构数值）。
    """
    df = profiles.sort_values("井深_m")
    md = np.asarray(df["井深_m"].to_numpy(), dtype=float)
    c = np.asarray(df["水泥平均浓度"].to_numpy(), dtype=float)
    if md.size < 2:
        return float("nan")
    pos_low = _contour_md(c, md, FRONT_C_LOW)
    pos_high = _contour_md(c, md, FRONT_C_HIGH)
    if np.isnan(pos_low) or np.isnan(pos_high):
        return float("nan")
    return abs(pos_low - pos_high)


def run_one_case(scale: float, cfl_on: bool, *, nz: int, ny: int) -> dict:
    """运行单个 (scale, cfl_on) 组合的完整 1D-2D 流水线，返回一行指标。

    手动流水线（density_contrast_sensitivity_scan.py 模式）：
    1D casing → coupled inlet provider → annulus stop time → 2D solver → 指标。
    """
    from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
    from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
    from cemdisp.runners.ht1_004_tailpipe import annulus_stop_time_s
    from cemdisp.transport1d import CasingFlowSolver

    t0 = time.perf_counter()

    loaded_well, fluids, schedule, _ = load_ht1_004_tailpipe()
    well_spec = loaded_well

    # 1D casing flow
    casing_solver = CasingFlowSolver(enable_gravity=True)
    casing_result = casing_solver.run(well_spec, fluids, schedule)

    # Coupled inlet provider
    provider = build_coupled_annulus_inlet_provider(
        casing_result, casing_solver, fluids
    )

    # Annulus stop time
    total_t = annulus_stop_time_s(casing_result=casing_result, fluids=fluids)

    # 2D D2DGA solver：R3 三闭包全开，M1 弥散参数透传
    solver = AnnulusD2DGASolver(
        dt=DT,
        nz=nz,
        ny=ny,
        total_t=total_t,
        enable_d2dga=True,
        enable_d2dga_auto_m=True,
        enable_d2dga_i3_flux=True,
        enable_true_buoyancy=True,
        open_outlet=True,
        enable_cfl_adaptive=cfl_on,
        dispersion_dt_scale=scale,
    )

    result = solver.run(well_spec, fluids, provider)
    elapsed = time.perf_counter() - t0
    metrics = extract_ablation_metrics(result)

    front_length = compute_front_length_m(result.depth_profiles)

    return {
        "scale": scale,
        "cfl_on": cfl_on,
        "effective_efficiency": metrics["effective_efficiency"],
        "channeling_index": metrics["channeling_index"],
        "mixing_index": metrics["mixing_index"],
        "instability_index": metrics["instability_index"],
        "cement_occupation": metrics["cement_occupation"],
        "front_length_m": round(front_length, 3) if not np.isnan(front_length) else None,
        "elapsed_s": round(elapsed, 1),
    }


def main() -> None:
    smoke = os.environ.get("M1_SMOKE", "0") == "1"
    if smoke:
        nz, ny = 30, 10
        scales = [0.0, 1.0]
    else:
        nz, ny = NZ, NY
        scales = SCALE_VALUES
    cfl_values = [True, False]  # CFL on/off 双组（design 方案乙）

    print("=" * 70)
    print("M1 弥散系数 κ 缩放敏感性扫描（dispersion_dt_scale × CFL on/off）")
    print(f"井: HT1-004   R3 三闭包全开")
    print(f"网格: nz={nz}, ny={ny}, dt={DT}s (dispersion_dt_ref={DT}s)")
    print(f"scale 扫描: {scales}")
    print(f"CFL 组: {'SMOKE' if smoke else 'on/off 双组'}")
    print(f"CSV: {OUTPUT_CSV}")
    print("=" * 70)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # 写 CSV 表头
    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()

    rows = []
    n_cases = len(scales) * len(cfl_values)
    idx = 0
    for scale in scales:
        for cfl_on in cfl_values:
            idx += 1
            tag = "CFL_on " if cfl_on else "CFL_off"
            print(f"\n--- [{idx}/{n_cases}] scale={scale}  {tag}  (nz={nz}, ny={ny}) ---")

            try:
                row = run_one_case(scale, cfl_on, nz=nz, ny=ny)
                print(
                    f"  -> 效率={row['effective_efficiency']}  窜槽={row['channeling_index']}  "
                    f"混浆={row['mixing_index']}  失稳={row['instability_index']}  "
                    f"水泥占据={row['cement_occupation']}  前沿长度={row['front_length_m']}m  "
                    f"耗时={row['elapsed_s']}s"
                )
            except Exception as exc:  # noqa: BLE001 —— 单 case 失败不中断全扫描
                import traceback

                traceback.print_exc()
                row = {
                    "scale": scale,
                    "cfl_on": cfl_on,
                    "effective_efficiency": None,
                    "channeling_index": None,
                    "mixing_index": None,
                    "instability_index": None,
                    "cement_occupation": None,
                    "front_length_m": None,
                    "elapsed_s": None,
                }
                print(f"  !! 失败: {exc}")

            rows.append(row)

            # 每完成一个 case 就追加写入 CSV（防崩）
            with open(OUTPUT_CSV, "a", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
                writer.writerow(row)

    # ---------- 汇总打印 ----------
    print("\n" + "=" * 70)
    print("扫描结果汇总")
    print("=" * 70)

    header = (
        f"{'scale':>6s}  {'CFL':>6s}  {'效率':>8s}  {'窜槽':>8s}  {'混浆':>8s}  "
        f"{'失稳':>8s}  {'水泥占据':>10s}  {'前沿长度m':>10s}  {'耗时(s)':>8s}"
    )
    print(header)
    print("-" * len(header))

    for r in rows:
        def _fmt(v):
            if v is None:
                return "     None"
            return f"{v:8.4f}"

        print(
            f"{r['scale']:6.2f}  {'on ' if r['cfl_on'] else 'off':>6s}  "
            f"{_fmt(r['effective_efficiency'])}  {_fmt(r['channeling_index'])}  "
            f"{_fmt(r['mixing_index'])}  {_fmt(r['instability_index'])}  "
            f"{_fmt(r['cement_occupation'])}  {_fmt(r['front_length_m'])}  "
            f"{r['elapsed_s']:8.1f}"
        )

    print(f"\n结果已保存至: {OUTPUT_CSV}")

    # ---------- 各 CFL 组内，指标对 scale 的线性斜率 ----------
    for cfl_on in cfl_values:
        valid = [r for r in rows if r["cfl_on"] == cfl_on and r["mixing_index"] is not None]
        if len(valid) < 2:
            print(f"\nCFL {'on' if cfl_on else 'off'} 组有效 case 不足 2 个，跳过线性斜率")
            continue
        s_arr = np.array([r["scale"] for r in valid])
        mixes = np.array([r["mixing_index"] for r in valid])
        effs = np.array([r["effective_efficiency"] for r in valid])
        chans = np.array([r["channeling_index"] for r in valid])
        fronts = np.array([r["front_length_m"] if r["front_length_m"] is not None else np.nan for r in valid])

        slope_mix, _ = np.polyfit(s_arr, mixes, 1)
        slope_eff, _ = np.polyfit(s_arr, effs, 1)
        slope_chan, _ = np.polyfit(s_arr, chans, 1)
        # 前沿长度可能有 nan，仅对有限值回归
        finite = np.isfinite(fronts)
        slope_front = float("nan")
        if finite.sum() >= 2:
            slope_front, _ = np.polyfit(s_arr[finite], fronts[finite], 1)

        print(f"\nCFL {'on' if cfl_on else 'off'} 组 · 指标对 scale 的线性斜率 (每 0.1 κ):")
        print(f"  混浆指数: {slope_mix * 0.1:+.4f} / 0.1")
        print(f"  有效顶替效率: {slope_eff * 0.1:+.4f} / 0.1")
        print(f"  窜槽指数: {slope_chan * 0.1:+.4f} / 0.1")
        print(f"  前沿长度: {slope_front * 0.1:+.4f} / 0.1")


if __name__ == "__main__":
    main()
