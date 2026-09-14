"""呼101 低分归因探针（spike，2026-09-11）。

问题：呼101 生产口径 η_E(全井)=0.4919 / η_E(CBL窗)=0.5035，是八井唯一低分
      （其余 0.95–0.998）。这是模型层缺陷，还是忠实反映输入？

已知线索（读制品所得，非本探针产出）：
  - 居中度不是唯一解释：呼102 均 SO=0.491 与呼101(0.438) 同级，η_E 却 0.9718。
  - 呼101 浮力数 b=52.7（八井最高，呼102 仅 2.5），且尾浆 ρ1.90 < 钻井液 ρ1.96。
  - `enable_true_buoyancy=False`（R3→R0/R2 代理）在 正文补齐变体 表里给 +15.6pp，
    而同一开关在呼103(e=0.278)/呼102(e=0.509) 上分别只给 +0.58/+1.12pp（随 e 完美线性），
    e=0.562 外推只该 +1.24pp —— 12 倍离群，疑似阈值分叉或缺陷。

探针设计（14 次 1D+2D 重跑，口径与 runners/hu101_tailpipe.py 逐位一致）：
  A 基线复现
  B 真浮力项关（复核离群值）
  C e_clip_max 扫描 {0.45,0.50,0.60,0.70,0.90}（呼101 原始 e=0.562 正卡在 0.55 截断线上）
  D 弥散系数扫描 {×0.5, ×2, ×3}（量化"标定旋钮"额度）
  E nz=500 网格收敛
  F 交叉诊断：hu102 井挂 hu101 居中度 / hu101 井挂 hu102 居中度
  G 密度倒置消除：hu101 尾浆 1.90→2.10

诊断量：全井/CBL窗 η_E、η_N、宽边/窄边浓度、深度剖面上下段、窜槽/混浆指数。
本脚本为一次性探针，不改任何生产默认值。
"""

from __future__ import annotations

import csv
import json
import time
from dataclasses import replace
from pathlib import Path
from typing import Any, cast

import numpy as np

from cemdisp.data.loaders import load_hu101_tailpipe, load_hu102_tailpipe
from cemdisp.data.well_spec import DepthValuePoint
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.runners.hu101_tailpipe import annulus_stop_time_s
from cemdisp.transport1d import CasingFlowSolver

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PROJECT_ROOT / "results" / "呼101_低分归因探针_2026-09-11"

NZ = 250
CBL_WINDOW = "CBL评价井段(单层套管可评价段)"

CASING_KW = dict(
    enable_gravity=True,
    mixing_contact_time=True,
    plug_face_zero_mixing=True,
    has_plug=True,
)


def _casing_stop(well_spec, fluids, schedule):
    """跑 1D 套管流，返回 (stop_time_s, casing_result, casing_solver)。"""
    casing_solver = CasingFlowSolver(**CASING_KW)  # type: ignore[arg-type]
    casing_result = casing_solver.run(well_spec, fluids, schedule)
    stop = annulus_stop_time_s(casing_result=casing_result, fluids=fluids)
    return stop, casing_result, casing_solver


def _field_stats(cement: np.ndarray) -> dict[str, float]:
    """从最终浓度场提取方位/深度诊断量。cement shape = (ny, nz)。"""
    ny = cement.shape[0]
    q = max(ny // 4, 1)
    nyz = cement.shape[1]
    bands = np.array_split(np.arange(nyz), 4)
    return {
        "wide_q1": float(cement[:q, :].mean()),
        "narrow_q4": float(cement[-q:, :].mean()),
        "depth_bottom_band": float(cement[:, bands[0]].mean()),   # 最深 1/4
        "depth_top_band": float(cement[:, bands[-1]].mean()),     # 最浅 1/4
        "narrow_side_gap_frac": float(np.mean(cement[-q:, :] < 0.05)),
    }


def run_case(
    *,
    label: str,
    group: str,
    well_spec: Any,
    fluids: Any,
    schedule: Any,
    solver_kw: dict[str, Any] | None = None,
    nz: int = NZ,
    casing_cache: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """按生产口径跑一个 case；1D 结果按 (井名, 流体指纹) 缓存复用。"""
    solver_kw = dict(solver_kw or {})
    key = f"{well_spec.well_name}|{len(fluids)}|{hash(tuple((f.name, f.density_kg_m3) for f in fluids))}"
    if casing_cache is not None and key in casing_cache:
        stop, casing_result, casing_solver = casing_cache[key]
    else:
        stop, casing_result, casing_solver = _casing_stop(well_spec, fluids, schedule)
        if casing_cache is not None:
            casing_cache[key] = (stop, casing_result, casing_solver)

    t0 = time.perf_counter()
    inlet = build_coupled_annulus_inlet_provider(
        casing_result, casing_solver, fluids, split_cement_phases=True,
    )
    solver = AnnulusD2DGASolver(total_t=stop, nz=nz, **solver_kw)
    result = solver.run(well_spec, fluids, inlet, schedule=schedule)
    elapsed = time.perf_counter() - t0

    final = cast(dict, result.summary["最终结果"])
    windows = cast(dict, result.summary["评价窗效率"])
    cbl = windows.get(CBL_WINDOW, {})
    row: dict[str, Any] = {
        "case": label,
        "group": group,
        "nz": nz,
        "stop_s": round(stop, 1),
        "eta_E_well": float(final["全井段最终有效顶替效率"]),
        "eta_N_well": float(final["窄四分位效率"]),
        "eta_E_cbl": float(cbl.get("eta_E", float("nan"))),
        "eta_N_cbl": float(cbl.get("eta_N", float("nan"))),
        "channeling_index": float(final["最终窜槽指数"]),
        "mixing_index": float(final["最终混浆指数"]),
        "b_number": float(final["浮力数_b"]),
        "elapsed_s": round(elapsed, 1),
    }
    row.update(_field_stats(np.asarray(result.cement_field, dtype=float)))
    print(
        f"  [{label:<26s}] ηE全井={row['eta_E_well']:.4f} ηE_CBL={row['eta_E_cbl']:.4f} "
        f"ηN全井={row['eta_N_well']:.4f} 宽={row['wide_q1']:.3f}/窄={row['narrow_q4']:.3f} "
        f"深={row['depth_bottom_band']:.3f}/浅={row['depth_top_band']:.3f} ({elapsed:.0f}s)",
        flush=True,
    )
    return row


def _profile_mean(profile: tuple[DepthValuePoint, ...]) -> float:
    """按 200 点等距插值求剖面均值。"""
    vals = np.array([p.value for p in profile], dtype=float)
    mds = np.array([p.depth_md_m for p in profile], dtype=float)
    return float(np.interp(np.linspace(mds.min(), mds.max(), 200), mds, vals).mean())


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    cache: dict[str, Any] = {}

    print("=" * 90)
    print("  呼101 低分归因探针（诊断：模型缺陷 vs 忠实反映输入）")
    print("=" * 90)

    well101, fluids101, sched101, _ = load_hu101_tailpipe()
    well102, fluids102, sched102, _ = load_hu102_tailpipe()

    def add(label, group, **kw):
        rows.append(run_case(label=label, group=group, casing_cache=cache, **kw))

    # ---- A 基线 ----
    print("\n[A] 基线复现（期望 ηE全井=0.49194 / ηE_CBL=0.50347）")
    add("A_基线", "A", well_spec=well101, fluids=fluids101, schedule=sched101)

    # ---- B 真浮力项关 ----
    print("\n[B] 真浮力项关（复核 正文补齐变体 的 +15.6pp 离群值）")
    add("B_真浮力关", "B", well_spec=well101, fluids=fluids101, schedule=sched101,
        solver_kw={"enable_true_buoyancy": False})

    # ---- C e_clip 扫描 ----
    print("\n[C] e_clip_max 扫描（呼101 原始 e=0.562，默认截断 0.55）")
    for ec in (0.45, 0.50, 0.60, 0.70, 0.90):
        add(f"C_eclip={ec:.2f}", "C", well_spec=well101, fluids=fluids101, schedule=sched101,
            solver_kw={"e_clip_max": ec})

    # ---- D 弥散扫描 ----
    print("\n[D] 弥散系数扫描（基线 axial=0.018 / azimuthal=0.015）")
    for scale in (0.5, 2.0, 3.0):
        add(f"D_disp_x{scale:g}", "D", well_spec=well101, fluids=fluids101, schedule=sched101,
            solver_kw={"dispersion_axial": 0.018 * scale, "dispersion_azimuthal": 0.015 * scale})

    # ---- E 网格收敛 ----
    print("\n[E] 网格收敛 nz=500")
    add("E_nz500", "E", well_spec=well101, fluids=fluids101, schedule=sched101, nz=500)

    # ---- F 交叉诊断：把居中度剖面互换 ----
    print("\n[F] 交叉诊断（居中度互换：分离'居中度'与'流体/几何'）")
    mu101 = _profile_mean(well101.standoff_profile)
    mu102 = _profile_mean(well102.standoff_profile)
    add("F_hu101挂hu102居中度", "F",
        well_spec=replace(well101, standoff_profile=well102.standoff_profile, standoff_measured=False),
        fluids=fluids101, schedule=sched101)
    const101 = (
        DepthValuePoint(float(well101.top_md_m), mu101),
        DepthValuePoint(float(well101.bottom_md_m), mu101),
    )
    add(f"F_hu102挂hu101居中度({mu101:.3f})", "F",
        well_spec=replace(well102, standoff_profile=const101, standoff_measured=False),
        fluids=fluids102, schedule=sched102)

    # ---- G 密度倒置消除 ----
    print("\n[G] 密度倒置消除（hu101 尾浆 1.90→2.10，消除 尾浆<钻井液）")
    fluids101_fixed = tuple(
        replace(f, density_kg_m3=2100.0) if f.name == "尾浆" else f for f in fluids101
    )
    add("G_尾浆密度1.90到2.10", "G", well_spec=well101, fluids=fluids101_fixed, schedule=sched101)

    # ---- 汇总 ----
    base = rows[0]
    for r in rows:
        r["d_eta_E_well"] = r["eta_E_well"] - base["eta_E_well"]
        r["d_eta_E_cbl"] = r["eta_E_cbl"] - base["eta_E_cbl"]

    print("\n" + "=" * 90)
    print("  汇总（Δ 相对 A_基线）")
    print("=" * 90)
    print(f"  {'case':<30s}{'ηE全井':>9s}{'Δ全井':>9s}{'ηE_CBL':>9s}{'ΔCBL':>9s}"
          f"{'ηN全井':>9s}{'宽':>7s}{'窄':>7s}{'深':>7s}{'浅':>7s}")
    for r in rows:
        print(f"  {r['case']:<30s}{r['eta_E_well']:>9.4f}{r['d_eta_E_well']:>+9.4f}"
              f"{r['eta_E_cbl']:>9.4f}{r['d_eta_E_cbl']:>+9.4f}{r['eta_N_well']:>9.4f}"
              f"{r['wide_q1']:>7.3f}{r['narrow_q4']:>7.3f}"
              f"{r['depth_bottom_band']:>7.3f}{r['depth_top_band']:>7.3f}")

    fields = ["case", "group", "nz", "stop_s", "eta_E_well", "d_eta_E_well", "eta_E_cbl",
              "d_eta_E_cbl", "eta_N_well", "eta_N_cbl", "channeling_index", "mixing_index",
              "b_number", "wide_q1", "narrow_q4", "depth_bottom_band", "depth_top_band",
              "narrow_side_gap_frac", "elapsed_s"]
    with (OUTPUT_DIR / "归因探针对比表.csv").open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})

    (OUTPUT_DIR / "归因探针摘要.json").write_text(
        json.dumps({
            "description": "呼101 低分归因探针（模型缺陷 vs 忠实反映输入）",
            "caliber": "runners/hu101_tailpipe.py 生产口径（T1三开关+annulus_stop_time+nz=250+schedule）",
            "baseline": base,
            "cases": rows,
        }, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(f"\n结果已落盘：{OUTPUT_DIR}")


if __name__ == "__main__":
    main()
