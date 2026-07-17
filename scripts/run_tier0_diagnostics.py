"""Tier 0 诊断独立验证脚本（SA5 集成验证）。

流程（仿 runner，纯调用，不改任何求解器代码）：
    load 井 → CasingFlowSolver → build_coupled_annulus_inlet_provider
    → AnnulusD2DGASolver.run → compute_all_tier0_diagnostics

用法：
    python scripts/run_tier0_diagnostics.py --well hu102
    python scripts/run_tier0_diagnostics.py --well ht1_004
    python scripts/run_tier0_diagnostics.py --well all      # 默认，两井都跑

输出：
    results/tier0_diagnostics/<井标识>_tier0.json  +  控制台中文摘要
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

# 允许直接以脚本方式运行（python scripts/run_tier0_diagnostics.py）
_PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from cemdisp.data.loaders import load_hu102_tailpipe
from cemdisp.data.loaders.ht1_004_loader import load_ht1_004_tailpipe
from cemdisp.diagnostics import compute_all_tier0_diagnostics
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.runners.ht1_004_tailpipe import annulus_stop_time_s as ht1_004_stop_time_s
from cemdisp.runners.hu102_tailpipe import annulus_stop_time_s as hu102_stop_time_s
from cemdisp.transport1d import CasingFlowSolver

OUTPUT_DIR = _PROJECT_ROOT / "results" / "tier0_diagnostics"

# 井标识 → (loader, 环空停止时间函数)
WELL_REGISTRY = {
    "hu102": (load_hu102_tailpipe, hu102_stop_time_s),
    "ht1_004": (load_ht1_004_tailpipe, ht1_004_stop_time_s),
}


def run_well_diagnostics(well_id: str, *, nz: int = 500, dt: float = 4.0) -> dict:
    """对单井跑 1D-2D 耦合求解并聚合 Tier 0 诊断，返回可 JSON 序列化字典。"""
    loader, stop_time_fn = WELL_REGISTRY[well_id]
    well_spec, fluids, schedule, _ = loader()
    print(f"\n{'=' * 70}")
    print(f"井：{well_spec.well_name}（标识 {well_id}）")
    print(f"{'=' * 70}")

    # 1D 套管内流动（启用重力项，与现场模式 runner 一致）
    t0 = time.perf_counter()
    casing_solver = CasingFlowSolver(enable_gravity=True)
    casing_result = casing_solver.run(well_spec, fluids, schedule)
    provider = build_coupled_annulus_inlet_provider(casing_result, casing_solver, fluids)
    total_t_s = float(stop_time_fn(casing_result=casing_result, fluids=fluids))
    print(f"[1/3] 1D 套管内流动完成，环空顶替总时长 total_t = {total_t_s:.1f} s")

    # 2D 环空 D2DGA 求解（参数与现场模式 runner 一致：nz=500，其余默认）
    solver = AnnulusD2DGASolver(dt=dt, nz=nz, total_t=total_t_s)
    result = solver.run(well_spec, fluids, provider)
    elapsed = time.perf_counter() - t0
    print(f"[2/3] 2D 环空顶替求解完成（nz={nz}, dt={dt}），耗时 {elapsed:.1f} s")

    # Tier 0 诊断聚合
    diag = compute_all_tier0_diagnostics(
        result, fluids=fluids, well_spec=well_spec, schedule=schedule
    )
    print("[3/3] Tier 0 诊断聚合完成")

    payload = {
        "well_id": well_id,
        "well_name": well_spec.well_name,
        "solver_config": {"nz": nz, "dt": dt, "total_t_s": total_t_s},
        "tier0_diagnostics": diag.to_dict(),
    }
    return payload


def print_summary(payload: dict) -> None:
    """打印控制台可读中文摘要。"""
    diag = payload["tier0_diagnostics"]
    print(f"\n--- {payload['well_name']} Tier 0 诊断摘要 ---")

    fc = diag.get("flow_classification")
    if fc is not None:
        print(f"[T0-1 流动分类] flow_class = {fc['flow_class']}"
              f"（delta_w_f = {_fmt(fc['delta_w_f'])}，sigma_wr+ = {_fmt(fc['sigma_wr_plus'])}，"
              f"w0 = {_fmt(fc['w0_m_s'])} m/s）")
    else:
        print("[T0-1 流动分类] 失败（见 notes）")

    mr = diag.get("muskat_regime")
    if mr is not None:
        print(f"[T0-2 Muskat] regime = {mr['regime']}"
              f"（m = {_fmt(mr['viscosity_ratio'])}，b = {_fmt(mr['buoyancy_number'])}，"
              f"e = {_fmt(mr['eccentricity'])}，c_critical = {_fmt(mr['c_critical'])}，"
              f"宽边失稳 = {mr['wide_side_unstable']}，窄边失稳 = {mr['narrow_side_unstable']}）")
    else:
        print("[T0-2 Muskat] 失败（见 notes）")

    br = diag.get("buoyancy_regime")
    if br is not None:
        print(f"[T0-3 浮力数分类] b = {_fmt(br['b_number'])} → {br['regime']}")
    else:
        print("[T0-3 浮力数分类] 失败（见 notes）")

    dm = diag.get("displacement_metrics")
    if dm is not None:
        print(f"[T0-4/5/7 顶替指标] 泥浆滞留 = {_fmt(dm['mud_retention_fraction'])}，"
              f"eta_N = {_fmt(dm['eta_narrow'])}，eta_E = {_fmt(dm['eta_global'])}，"
              f"界面长度比 = {_fmt(dm['interface_length_ratio'])}，"
              f"t_br = {_fmt(dm['t_br_s'])} s（t_br_hat = {_fmt(dm['t_br_hat'])}），"
              f"质量分区 green/yellow/red = {dm['quality_zone_counts']}")
    else:
        print("[T0-4/5/7 顶替指标] 失败（见 notes）")

    sd = diag.get("shutdown_decay")
    if sd is not None:
        print(f"[T0-6 停泵衰减] 判据满足 = {sd['condition_satisfied']}，"
              f"冻结时间 = {_fmt(sd['freeze_time_s'])} s，tau_Y,min = {_fmt(sd['tau_y_min'])} Pa")
        print(f"    {sd['physical_interpretation']}")
    else:
        print("[T0-6 停泵衰减] 失败（见 notes）")

    for note in diag.get("notes", []):
        print(f"  [note] {note}")


def _fmt(value: object) -> str:
    """格式化数值（None/非有限值显示为 N/A）。"""
    if value is None:
        return "N/A"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def main() -> None:
    # Windows GBK 控制台兜底：无法编码的字符替换为 "?"，避免 UnicodeEncodeError
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(errors="replace")
    parser = argparse.ArgumentParser(description="Tier 0 诊断独立验证（hu102 / ht1_004）")
    parser.add_argument(
        "--well",
        choices=[*WELL_REGISTRY.keys(), "all"],
        default="all",
        help="选择井标识（默认 all，两井都跑）",
    )
    parser.add_argument("--nz", type=int, default=500, help="轴向网格数（默认 500）")
    parser.add_argument("--dt", type=float, default=4.0, help="时间步长 s（默认 4.0）")
    args = parser.parse_args()

    well_ids = list(WELL_REGISTRY.keys()) if args.well == "all" else [args.well]
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for well_id in well_ids:
        payload = run_well_diagnostics(well_id, nz=args.nz, dt=args.dt)
        json_path = OUTPUT_DIR / f"{well_id}_tier0.json"
        json_path.write_text(
            json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        print_summary(payload)
        print(f"\n诊断 JSON 已写出：{json_path}")

    print("\n全部完成。")


if __name__ == "__main__":
    main()
