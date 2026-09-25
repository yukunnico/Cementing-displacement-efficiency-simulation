"""当前口径敏感性重跑（2026-09-16）——为论文第 4 章取当前模型数字。

背景
----
论文第 3/4 章现有数字全部来自 2026-09-09/09-10 口径（求解器旧路径，HEAD 早于
2026-09-14 源模型口径重构 Task 0–14）。重构把环空二维求解切到流函数新路径
（`enable_stream_function=True` 默认），权威 8 井数字随之改变
（`results/源模型口径重跑_2026-09-14/汇总.csv`）。本脚本在**当前口径**下重跑
同一批单变量变体，供论文数字替换。

口径（与 `results/源模型口径重跑_2026-09-14` 严格同源）
------------------------------------------------------
- 求解器：AnnulusD2DGASolver(total_t, nz=250, enable_cfl_adaptive=True, **CORRECTED_KW)
  —— CORRECTED_KW / NZ / _stop_t / _total_t 直接 import 自
  `scripts/entrypoints/rerun_all_wells_corrected.py`，避免口径漂移。
- 停止时刻：tt = min(总泵注时长 + 1200 s, 尾浆入库时刻)（F2 口径，无 +600 s 尾窗）。
- 套管内 1D：CasingFlowSolver(enable_gravity=True, mixing_contact_time=True,
  plug_face_zero_mixing=True, has_plug=True)（T1 生产口径）。
- 基线：`results/源模型口径重跑_2026-09-14/单井结果/<well>_corrected_on.json`
  （当前口径权威基线，**不**用 results/<井名>_1D2D耦合模型/ 的旧路径数字）。

变体
----
12 个现场杠杆（与 2026-09-10 敏感性组逐个同名，便于新旧逐格对照）：
    rate_x0.6 / rate_x0.8 / rate_x1.2 / rate_x1.4      排量缩放（泵序逐段）
    cement_n_p0.1 / cement_n_m0.1                       水泥浆 n ±0.1
    spacer_dens_p100 / spacer_dens_m100                 隔离液密度 ±100 kg/m³
    mud_pv_x0.7 / mud_pv_x1.3                           钻井液 PV ±30%
    standoff_p0.1 / standoff_m0.1                       standoff 剖面整体 ±0.1
另加 standoff 连续扫描 7 档（±0.05 / ±0.15 / ±0.20 / ±0.30），与 ±0.1 两档合成
9 档连续曲线，用于判定「居中度阈值悬崖」在当前口径下是否仍然存在。

输出
----
results/敏感性变体_当前口径_2026-09-16/（逐变体摘要 JSON + 汇总 CSV/MD）。

用法
----
    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/entrypoints/run_sensitivity_current_20260916.py
"""
from __future__ import annotations

import csv
import importlib
import json
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.pumping_schedule import PumpingSchedule
from cemdisp.data.well_spec import WellSpec
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver

# 口径唯一来源：与权威 8 井重跑共用同一组常量与函数，防止漂移。
_SCRIPTS_DIR = Path(__file__).resolve().parents[1]
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))
from entrypoints.rerun_all_wells_corrected import (  # noqa: E402
    CORRECTED_KW,
    NZ,
    _stop_t,
    _total_t,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = PROJECT_ROOT / "results" / "敏感性变体_当前口径_2026-09-16"
BASELINE_DIR = PROJECT_ROOT / "results" / "源模型口径重跑_2026-09-14" / "单井结果"

# 井名 → (loader 模块, loader 函数, 当前口径基线 JSON 名)
WELLS: dict[str, tuple[str, str, str]] = {
    "呼103": ("cemdisp.data.loaders.hu103_loader", "load_hu103_tailpipe", "hu103_corrected_on.json"),
    "呼101": ("cemdisp.data.loaders.hu101_loader", "load_hu101_tailpipe", "hu101_corrected_on.json"),
}

CEMENT_ROLES = {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}


# ---------------------------------------------------------------- 输入变换
def scale_schedule(schedule: PumpingSchedule, factor: float) -> PumpingSchedule:
    """泵序逐段排量缩放（体积不变，时长随 1/factor 拉长）。"""
    return replace(schedule, steps=tuple(
        replace(s, rate_m3_min=s.rate_m3_min * factor) if s.rate_m3_min > 0 else s
        for s in schedule.steps
    ))


def shift_cement_n(fluids: tuple[FluidSpec, ...], delta: float) -> tuple[FluidSpec, ...]:
    """水泥浆幂律指数 n ± delta（仅水泥三角色且 power_law_n 非空时生效）。"""
    return tuple(
        replace(f, power_law_n=f.power_law_n + delta)
        if f.role in CEMENT_ROLES and f.power_law_n is not None else f
        for f in fluids
    )


def shift_spacer_density(fluids: tuple[FluidSpec, ...], d_kg_m3: float) -> tuple[FluidSpec, ...]:
    """隔离液密度 ± d_kg_m3（仅 SPACER 角色）。"""
    return tuple(
        replace(f, density_kg_m3=f.density_kg_m3 + d_kg_m3) if f.role == FluidRole.SPACER else f
        for f in fluids
    )


def scale_mud_pv(fluids: tuple[FluidSpec, ...], factor: float) -> tuple[FluidSpec, ...]:
    """钻井液塑性黏度 × factor（仅 MUD 角色）。"""
    return tuple(
        replace(f, plastic_viscosity_pa_s=f.plastic_viscosity_pa_s * factor)
        if f.role == FluidRole.MUD and f.plastic_viscosity_pa_s is not None else f
        for f in fluids
    )


def shift_standoff(well: WellSpec, delta: float) -> WellSpec:
    """standoff 剖面整体平移 delta（clip [0,1]），形状与深度点不变。"""
    from cemdisp.data.well_spec import DepthValuePoint

    if not well.standoff_profile:
        raise ValueError(f"{well.well_name}: 无 standoff_profile")
    return replace(well, standoff_profile=tuple(
        DepthValuePoint(p.depth_md_m, float(np.clip(p.value + delta, 0.0, 1.0)))
        for p in well.standoff_profile
    ))


def _identity(x):
    return x


# ---------------------------------------------------------------- 变体表
def build_variants() -> list[tuple[str, object, object, object]]:
    """(变体名, well 变换, fluids 变换, schedule 变换)。"""
    return [
        # 12 杠杆（与 2026-09-10 组同名，逐格可比）
        ("rate_x1.2", _identity, _identity, lambda s: scale_schedule(s, 1.2)),
        ("rate_x0.8", _identity, _identity, lambda s: scale_schedule(s, 0.8)),
        ("rate_x1.4", _identity, _identity, lambda s: scale_schedule(s, 1.4)),
        ("rate_x0.6", _identity, _identity, lambda s: scale_schedule(s, 0.6)),
        ("cement_n_p0.1", _identity, lambda fs: shift_cement_n(fs, +0.1), _identity),
        ("cement_n_m0.1", _identity, lambda fs: shift_cement_n(fs, -0.1), _identity),
        ("spacer_dens_p100", _identity, lambda fs: shift_spacer_density(fs, +100.0), _identity),
        ("spacer_dens_m100", _identity, lambda fs: shift_spacer_density(fs, -100.0), _identity),
        ("mud_pv_x1.3", _identity, lambda fs: scale_mud_pv(fs, 1.3), _identity),
        ("mud_pv_x0.7", _identity, lambda fs: scale_mud_pv(fs, 0.7), _identity),
        ("standoff_p0.1", lambda w: shift_standoff(w, +0.1), _identity, _identity),
        ("standoff_m0.1", lambda w: shift_standoff(w, -0.1), _identity, _identity),
        # standoff 连续扫描补充档
        ("standoff_m0.05", lambda w: shift_standoff(w, -0.05), _identity, _identity),
        ("standoff_p0.05", lambda w: shift_standoff(w, +0.05), _identity, _identity),
        ("standoff_m0.15", lambda w: shift_standoff(w, -0.15), _identity, _identity),
        ("standoff_m0.20", lambda w: shift_standoff(w, -0.20), _identity, _identity),
        ("standoff_p0.20", lambda w: shift_standoff(w, +0.20), _identity, _identity),
        ("standoff_m0.30", lambda w: shift_standoff(w, -0.30), _identity, _identity),
        ("standoff_p0.30", lambda w: shift_standoff(w, +0.30), _identity, _identity),
    ]


def run_variant(loader, well_fn, fluid_fn, sched_fn) -> tuple[dict, dict]:
    """单变体一次完整流水线（1D 重跑 + 环空二维），返回 (最终结果, 完整 summary)。"""
    well, fluids, schedule, _ = loader()
    well2, fluids2, schedule2 = well_fn(well), fluid_fn(fluids), sched_fn(schedule)

    casing = CasingFlowSolver(
        enable_gravity=True,
        mixing_contact_time=True,
        plug_face_zero_mixing=True,
        has_plug=True,
    )
    cr = casing.run(well2, fluids2, schedule2)
    inlet = build_coupled_annulus_inlet_provider(
        cr, casing, fluids2, split_cement_phases=True,
    )
    tt = min(_total_t(schedule2) + 1200.0, _stop_t(cr, fluids2))
    solver = AnnulusD2DGASolver(
        total_t=tt, nz=NZ, enable_cfl_adaptive=True, **CORRECTED_KW,
    )
    res = solver.run(well2, fluids2, inlet, schedule=schedule2)
    return res.summary["最终结果"], res.summary


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    variants = build_variants()
    rows: list[dict] = []

    for well_key, (mod_name, fn_name, baseline_file) in WELLS.items():
        loader = getattr(importlib.import_module(mod_name), fn_name)
        base = json.loads((BASELINE_DIR / baseline_file).read_text(encoding="utf-8"))
        base_eta_e, base_eta_n = float(base["eta_E"]), float(base["eta_N"])
        print(f"\n=== {well_key}  基线（当前口径） η_E={base_eta_e:.4f} η_N={base_eta_n:.4f} ===",
              flush=True)

        for var_name, well_fn, fluid_fn, sched_fn in variants:
            case_json = OUT_DIR / f"{well_key}_{var_name}_结果摘要.json"
            if case_json.exists():                       # 断点续跑
                summary = json.loads(case_json.read_text(encoding="utf-8"))
                final, elapsed = summary["最终结果"], None
                print(f"  [复用] {well_key} × {var_name}", flush=True)
            else:
                t0 = time.perf_counter()
                final, summary = run_variant(loader, well_fn, fluid_fn, sched_fn)
                elapsed = round(time.perf_counter() - t0, 1)
                case_json.write_text(
                    json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8",
                )
                print(f"  [计算] {well_key} × {var_name}: "
                      f"η_E={float(final['全井段最终有效顶替效率']):.4f} "
                      f"η_N={float(final['窄四分位效率']):.4f} ({elapsed}s)", flush=True)

            eta_e = float(final["全井段最终有效顶替效率"])
            eta_n = float(final["窄四分位效率"])
            rows.append({
                "井名": well_key,
                "变体": var_name,
                "η_E": eta_e,
                "η_N": eta_n,
                "窜槽": float(final["最终窜槽指数"]),
                "混浆": float(final["最终混浆指数"]),
                "失稳": float(final["最终失稳指数"]),
                "Δη_E_pp": (eta_e - base_eta_e) * 100.0,
                "Δη_N_pp": (eta_n - base_eta_n) * 100.0,
                "基线η_E": base_eta_e,
                "基线η_N": base_eta_n,
                "耗时_s": elapsed if elapsed is not None else "",
            })

    csv_path = OUT_DIR / "汇总表_敏感性变体_当前口径.csv"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    md = [
        "# 敏感性变体汇总（当前口径 2026-09-16，流函数新路径）",
        "",
        f"口径：nz={NZ}、CFL 自适应 on、CORRECTED_KW、T1 套管生产开关；"
        "与 `results/源模型口径重跑_2026-09-14` 严格同源。",
        "基线：`results/源模型口径重跑_2026-09-14/单井结果/<well>_corrected_on.json`。",
        "",
        "| 井名 | 变体 | η_E | η_N | Δη_E/pp | Δη_N/pp | 耗时/s |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        md.append(f"| {r['井名']} | {r['变体']} | {r['η_E']:.4f} | {r['η_N']:.4f} "
                  f"| {r['Δη_E_pp']:+.2f} | {r['Δη_N_pp']:+.2f} | {r['耗时_s']} |")
    (OUT_DIR / "汇总表_敏感性变体_当前口径.md").write_text("\n".join(md), encoding="utf-8")

    print(f"\n完成：{len(rows)} 个变体 → {OUT_DIR}", flush=True)


if __name__ == "__main__":
    main()
