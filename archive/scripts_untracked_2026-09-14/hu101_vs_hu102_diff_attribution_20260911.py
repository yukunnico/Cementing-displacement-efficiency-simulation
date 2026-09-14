"""呼101 vs 呼102 差异归因探针（spike，2026-09-11 第三轮）。

背景：前两轮已排除
  - 居中度/e_clip：K101_SO 细扫 0.40–0.52 全平（0.48–0.49），SO=0.50（e=0.50≈呼102 的
    0.506）仍在平台；C_eclip 0.45→0.90 仅 ±2pp
  - 密度倒置：G 变体（尾浆 1.90→2.10）仅 +0.02pp
  - 网格：nz500 仅 −1.9pp
  - 出口/屈服门/B1：全部为负向（−1.7 / −5.1 / −3.7 pp）
尚未解释的核心矛盾：
  呼101(e=0.547, 间隙45.7mm) 窄边浓度 0.028；呼102(e=0.506, 间隙26.1mm) 窄边浓度 0.664；
  呼103(e=0.278, 间隙50.1mm) 窄边浓度 0.909。e 相近而结果天壤 ⇒ 另有主控量。

本轮：把呼101 的四个候选差异逐个替换为呼102 的值（单变量），并做反向验证。
  1 井径剖面 → 呼102 尺度（缩小环空）
  2 钻井液流变 → 呼102（PV 0.066 / YP 10.0）
  3 水泥浆流变 → 呼102（n 0.737 / K 0.947，领/尾同值）
  4 井斜 → 呼102 均值常数
  5 井径 + 钻井液流变 组合
  6 反向：呼102 井径 → 呼101 尺度

口径与 runners/hu101_tailpipe.py 逐位一致（1D 随流体/井径变化重跑）。一次性探针。
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
OUTPUT_DIR = PROJECT_ROOT / "results" / "呼101_vs_呼102_差异归因_2026-09-11"
NZ = 250
CBL_WINDOW = "CBL评价井段(单层套管可评价段)"
CASING_KW = dict(enable_gravity=True, mixing_contact_time=True, plug_face_zero_mixing=True, has_plug=True)

# 呼102 目标值（现场 loader 口径）
T102_MUD_PV = 0.066
T102_MUD_YP = 10.0
T102_MUD_RHO = 2020.0
T102_CEM_N = 0.737
T102_CEM_K = 0.947
T102_CEM_RHO = 2100.0


def scale_hole(profile, factor: float):
    """按比例缩放井径剖面。"""
    return tuple(DepthValuePoint(p.depth_md_m, p.value * factor) for p in profile)


def mean_hole(profile) -> float:
    return float(np.mean([p.value for p in profile]))


def mean_gap(well_spec) -> float:
    """平均径向环空间隙 (mm)。"""
    holes = np.array([p.value for p in well_spec.hole_diameter_profile])
    return float((holes.mean() - well_spec.liner_od_mm) / 2.0)


def stats(cement: np.ndarray) -> dict[str, float]:
    ny = cement.shape[0]
    q = max(ny // 4, 1)
    return {
        "wide_q1": float(cement[:q, :].mean()),
        "mid": float(cement[ny // 2, :].mean()),
        "narrow_q4": float(cement[-q:, :].mean()),
    }


def run_case(*, label, well_spec, fluids, schedule, cache=None) -> dict[str, Any]:
    key = (f"{well_spec.well_name}|"
           f"{hash(tuple((f.name, f.density_kg_m3, f.power_law_n, f.plastic_viscosity_pa_s) for f in fluids))}|"
           f"{hash(tuple(round(p.value, 3) for p in well_spec.hole_diameter_profile))}")
    if cache is not None and key in cache:
        stop, cres, csol = cache[key]
    else:
        csol = CasingFlowSolver(**CASING_KW)  # type: ignore[arg-type]
        cres = csol.run(well_spec, fluids, schedule)
        stop = annulus_stop_time_s(casing_result=cres, fluids=fluids)
        if cache is not None:
            cache[key] = (stop, cres, csol)

    t0 = time.perf_counter()
    inlet = build_coupled_annulus_inlet_provider(cres, csol, fluids, split_cement_phases=True)
    result = AnnulusD2DGASolver(total_t=stop, nz=NZ).run(well_spec, fluids, inlet, schedule=schedule)
    elapsed = time.perf_counter() - t0

    final = cast(dict, result.summary["最终结果"])
    cbl = cast(dict, result.summary["评价窗效率"]).get(CBL_WINDOW, {})
    row: dict[str, Any] = {
        "case": label,
        "eta_E_well": float(final["全井段最终有效顶替效率"]),
        "eta_N_well": float(final["窄四分位效率"]),
        "eta_E_cbl": float(cbl.get("eta_E", float("nan"))),
        "gap_mean_mm": round(mean_gap(well_spec), 2),
        "elapsed_s": round(elapsed, 1),
    }
    row.update(stats(np.asarray(result.cement_field, dtype=float)))
    print(f"  [{label:<30s}] ηE={row['eta_E_well']:.4f} ηN={row['eta_N_well']:.4f} "
          f"宽={row['wide_q1']:.3f}/中={row['mid']:.3f}/窄={row['narrow_q4']:.3f} "
          f"(间隙{row['gap_mean_mm']}mm, {elapsed:.0f}s)", flush=True)
    return row


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    cache: dict[str, Any] = {}

    w101, f101, s101, _ = load_hu101_tailpipe()
    w102, f102, s102, _ = load_hu102_tailpipe()

    print("=" * 92)
    print("  呼101 vs 呼102 差异归因（单变量替换）")
    print("=" * 92)
    print(f"  呼101 参考：间隙 {mean_gap(w101):.1f}mm  ηE(已知)=0.4919")
    print(f"  呼102 参考：间隙 {mean_gap(w102):.1f}mm  ηE(已知)=0.9718\n")

    # 0 hu102 基线（本脚本口径复现）
    rows.append(run_case(label="0_呼102基线", well_spec=w102, fluids=f102, schedule=s102, cache=cache))

    h101_mean = mean_hole(w101.hole_diameter_profile)
    h102_mean = mean_hole(w102.hole_diameter_profile)
    factor_to_102 = h102_mean / h101_mean

    # 1 井径 → 呼102 尺度
    rows.append(run_case(
        label=f"1_井径到呼102尺度(x{factor_to_102:.3f})",
        well_spec=replace(w101, hole_diameter_profile=scale_hole(w101.hole_diameter_profile, factor_to_102)),
        fluids=f101, schedule=s101, cache=cache))

    # 2 钻井液流变 → 呼102
    f101_mud102 = tuple(
        replace(f, plastic_viscosity_pa_s=T102_MUD_PV, yield_stress_pa=T102_MUD_YP, density_kg_m3=T102_MUD_RHO)
        if f.role.name == "MUD" else f for f in f101)
    rows.append(run_case(label="2_钻井液流变到呼102", well_spec=w101, fluids=f101_mud102,
                         schedule=s101, cache=cache))

    # 3 水泥浆流变 → 呼102
    f101_cem102 = tuple(
        replace(f, power_law_n=T102_CEM_N, consistency_k=T102_CEM_K, density_kg_m3=T102_CEM_RHO)
        if f.role.name in ("LEAD", "TAIL") else f for f in f101)
    rows.append(run_case(label="3_水泥浆流变到呼102", well_spec=w101, fluids=f101_cem102,
                         schedule=s101, cache=cache))

    # 4 井斜 → 呼102 均值常数
    inc_mean_102 = float(np.mean([p.value for p in w102.inclination_profile]))
    const_inc = (DepthValuePoint(float(w101.top_md_m), inc_mean_102),
                 DepthValuePoint(float(w101.bottom_md_m), inc_mean_102))
    rows.append(run_case(label=f"4_井斜到{inc_mean_102:.2f}deg常数",
                         well_spec=replace(w101, inclination_profile=const_inc),
                         fluids=f101, schedule=s101, cache=cache))

    # 5 井径 + 钻井液流变 组合
    rows.append(run_case(
        label="5_井径+钻井液到呼102",
        well_spec=replace(w101, hole_diameter_profile=scale_hole(w101.hole_diameter_profile, factor_to_102)),
        fluids=f101_mud102, schedule=s101, cache=cache))

    # 6 反向：呼102 井径 → 呼101 尺度
    factor_to_101 = h101_mean / h102_mean
    rows.append(run_case(
        label=f"6_呼102井径到呼101尺度(x{factor_to_101:.3f})",
        well_spec=replace(w102, hole_diameter_profile=scale_hole(w102.hole_diameter_profile, factor_to_101)),
        fluids=f102, schedule=s102, cache=cache))

    print("\n" + "=" * 92)
    print(f"  {'case':<34s}{'ηE全井':>10s}{'ηN全井':>10s}{'宽':>8s}{'中':>8s}{'窄':>8s}{'间隙mm':>9s}")
    for r in rows:
        print(f"  {r['case']:<34s}{r['eta_E_well']:>10.4f}{r['eta_N_well']:>10.4f}"
              f"{r['wide_q1']:>8.3f}{r['mid']:>8.3f}{r['narrow_q4']:>8.3f}{r['gap_mean_mm']:>9.1f}")

    fields = ["case", "eta_E_well", "eta_N_well", "eta_E_cbl", "wide_q1", "mid", "narrow_q4",
              "gap_mean_mm", "elapsed_s"]
    with (OUTPUT_DIR / "差异归因对比表.csv").open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})
    (OUTPUT_DIR / "差异归因摘要.json").write_text(
        json.dumps({"description": "呼101 vs 呼102 单变量替换归因", "cases": rows},
                   ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n结果已落盘：{OUTPUT_DIR}")


if __name__ == "__main__":
    main()
