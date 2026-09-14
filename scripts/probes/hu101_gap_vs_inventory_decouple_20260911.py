"""间隙 vs 库存比 解耦探针（spike，2026-09-11 第四轮）。

第三轮发现：把呼101 井径缩到呼102 尺度 → η_E 0.4919→0.9248；反向 → 0.9718→0.5788。
但缩井径的同时环空体积按 ~0.69 缩放，水泥量不变 ⇒ 库存比 1.068→1.55，
放井径时 1.374→0.955（<1）。因此"间隙"与"水泥相对充裕度"两个因素被混淆。

本轮用四组正交对照把它们分开：
  Q1 呼101 井径×0.831 且水泥量×0.690（库存比锁定 1.068）→ 纯间隙效应
  Q2 呼101 井径不变、水泥量×1.287（库存比 1.068→1.374，对齐呼102）→ 纯库存效应
  Q3 呼102 井径×1.203 且水泥量×1.447（库存比锁定 1.374）→ 反向纯间隙效应
  Q4 呼101 井径×0.831 且水泥量×1.287（间隙+库存都按呼102）→ 对照 Q1/Q2 可加性

口径与 runners/hu101_tailpipe.py 逐位一致（1D 随井径/泵量重跑）。一次性探针。
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
from cemdisp.data.pumping_schedule import PumpingSchedule
from cemdisp.data.well_spec import DepthValuePoint
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.runners.hu101_tailpipe import annulus_stop_time_s
from cemdisp.transport1d import CasingFlowSolver

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = PROJECT_ROOT / "results" / "呼101_间隙vs库存比解耦_2026-09-11"
NZ = 250
CBL_WINDOW = "CBL评价井段(单层套管可评价段)"
CASING_KW = dict(enable_gravity=True, mixing_contact_time=True, plug_face_zero_mixing=True, has_plug=True)

CEMENT_NAMES = ("领浆", "尾浆", "尾管水泥浆")


def scale_hole(profile, factor: float):
    return tuple(DepthValuePoint(p.depth_md_m, p.value * factor) for p in profile)


def scale_cement_schedule(schedule: PumpingSchedule, factor: float) -> PumpingSchedule:
    """按比例缩放所有水泥浆泵注段的体积（隔离液/替浆段不动）。"""
    steps = tuple(
        replace(step, volume_m3=step.volume_m3 * factor) if step.fluid_name in CEMENT_NAMES else step
        for step in schedule.steps
    )
    return replace(schedule, steps=steps)


def annulus_volume_m3(well_spec) -> float:
    """按井径/外径剖面估算物理环空体积。"""
    md = np.array([p.depth_md_m for p in well_spec.hole_diameter_profile], dtype=float)
    hole = np.array([p.value for p in well_spec.hole_diameter_profile], dtype=float)
    od = well_spec.liner_od_mm
    area = np.pi / 4.0 * (hole**2 - od**2) * 1e-6  # m2
    return float(np.trapezoid(area, md))


def cement_volume_m3(schedule: PumpingSchedule) -> float:
    return float(sum(s.volume_m3 for s in schedule.steps if s.fluid_name in CEMENT_NAMES))


def stats(cement: np.ndarray) -> dict[str, float]:
    ny = cement.shape[0]
    q = max(ny // 4, 1)
    bands = np.array_split(np.arange(cement.shape[1]), 4)
    return {
        "wide_q1": float(cement[:q, :].mean()),
        "mid": float(cement[ny // 2, :].mean()),
        "narrow_q4": float(cement[-q:, :].mean()),
        "depth_bottom_band": float(cement[:, bands[0]].mean()),
        "depth_top_band": float(cement[:, bands[-1]].mean()),
    }


def run_case(*, label, well_spec, fluids, schedule, cache=None) -> dict[str, Any]:
    key = (f"{well_spec.well_name}|"
           f"{hash(tuple(round(p.value, 2) for p in well_spec.hole_diameter_profile))}|"
           f"{cement_volume_m3(schedule):.3f}")
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
    holes = np.array([p.value for p in well_spec.hole_diameter_profile])
    V = annulus_volume_m3(well_spec)
    C = cement_volume_m3(schedule)
    row: dict[str, Any] = {
        "case": label,
        "gap_mm": round(float((holes.mean() - well_spec.liner_od_mm) / 2.0), 2),
        "annulus_m3": round(V, 2),
        "cement_m3": round(C, 1),
        "inventory_ratio": round(C / V, 3),
        "eta_E_well": float(final["全井段最终有效顶替效率"]),
        "eta_N_well": float(final["窄四分位效率"]),
        "eta_E_cbl": float(cbl.get("eta_E", float("nan"))),
        "elapsed_s": round(elapsed, 1),
    }
    row.update(stats(np.asarray(result.cement_field, dtype=float)))
    print(f"  [{label:<28s}] ηE={row['eta_E_well']:.4f} ηN={row['eta_N_well']:.4f} "
          f"窄={row['narrow_q4']:.3f} 间隙{row['gap_mm']}mm V={row['annulus_m3']}m3 "
          f"C={row['cement_m3']}m3 库存比={row['inventory_ratio']} ({elapsed:.0f}s)", flush=True)
    return row


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    cache: dict[str, Any] = {}

    w101, f101, s101, _ = load_hu101_tailpipe()
    w102, f102, s102, _ = load_hu102_tailpipe()
    h101 = float(np.mean([p.value for p in w101.hole_diameter_profile]))
    h102 = float(np.mean([p.value for p in w102.hole_diameter_profile]))
    k_down = h102 / h101     # 0.831：呼101 井径缩到呼102
    k_up = h101 / h102       # 1.203：呼102 井径放到呼101

    print("=" * 100)
    print("  间隙 vs 库存比 解耦（呼101 / 呼102）")
    print("=" * 100)

    rows.append(run_case(label="A_呼101基线", well_spec=w101, fluids=f101, schedule=s101, cache=cache))
    rows.append(run_case(label="B_呼102基线", well_spec=w102, fluids=f102, schedule=s102, cache=cache))

    # Q1 纯间隙：缩井径 + 缩水泥量，库存比锁定 1.068
    w101_small = replace(w101, hole_diameter_profile=scale_hole(w101.hole_diameter_profile, k_down))
    s101_small = scale_cement_schedule(s101, k_down**2)
    rows.append(run_case(label="Q1_小间隙+库存比锁定", well_spec=w101_small,
                         fluids=f101, schedule=s101_small, cache=cache))

    # Q2 纯库存：井径不变，水泥量×1.287 → 库存比 1.374
    s101_big = scale_cement_schedule(s101, 1.374 / 1.068)
    rows.append(run_case(label="Q2_原间隙+库存比1.374", well_spec=w101,
                         fluids=f101, schedule=s101_big, cache=cache))

    # Q3 反向纯间隙：呼102 放大井径 + 放大水泥量，库存比锁定 1.374
    w102_big = replace(w102, hole_diameter_profile=scale_hole(w102.hole_diameter_profile, k_up))
    s102_big = scale_cement_schedule(s102, k_up**2)
    rows.append(run_case(label="Q3_呼102大间隙+库存比锁定", well_spec=w102_big,
                         fluids=f102, schedule=s102_big, cache=cache))

    # Q4 间隙+库存都按呼102
    rows.append(run_case(label="Q4_小间隙+库存比1.374", well_spec=w101_small,
                         fluids=f101, schedule=s101_big, cache=cache))

    print("\n" + "=" * 100)
    print(f"  {'case':<30s}{'ηE':>9s}{'ηN':>9s}{'窄':>8s}{'间隙mm':>9s}{'环空m3':>9s}{'水泥m3':>9s}{'库存比':>8s}")
    for r in rows:
        print(f"  {r['case']:<30s}{r['eta_E_well']:>9.4f}{r['eta_N_well']:>9.4f}{r['narrow_q4']:>8.3f}"
              f"{r['gap_mm']:>9.2f}{r['annulus_m3']:>9.2f}{r['cement_m3']:>9.1f}{r['inventory_ratio']:>8.3f}")

    fields = ["case", "gap_mm", "annulus_m3", "cement_m3", "inventory_ratio",
              "eta_E_well", "eta_N_well", "eta_E_cbl", "wide_q1", "mid", "narrow_q4",
              "depth_bottom_band", "depth_top_band", "elapsed_s"]
    with (OUTPUT_DIR / "解耦对比表.csv").open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})
    (OUTPUT_DIR / "解耦摘要.json").write_text(
        json.dumps({"description": "间隙 vs 库存比解耦", "cases": rows}, ensure_ascii=False, indent=2),
        encoding="utf-8")
    print(f"\n结果已落盘：{OUTPUT_DIR}")


if __name__ == "__main__":
    main()
