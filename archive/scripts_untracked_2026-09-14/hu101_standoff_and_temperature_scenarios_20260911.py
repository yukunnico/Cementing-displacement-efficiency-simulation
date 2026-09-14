"""呼101 居中度口径 × 温度口径 全场景探针（spike，2026-09-11 第六轮）。

用户裁定要做三件事：
  1. 把居中度平均"优化"到 0.8；
  2. 做温度口径验证（钻井液流变换成井下温度口径）；
  3. 把口径换成固井设计软件输出的居中度检测图剖面。

⚠️ 性质分级（写进产物，防止误用）：
  ✅ 可辩护：A2/A3（检测图剖面，有出处）、B（温度口径，现场多温度数据外推）
  ❌ 反推场景：A4/A5（把均居中度抬到 0.80）—— 不得作为验证数字写进论文

温度外推：现场 rheometer_readings.csv 只有钻井液 40/50/60/70/80℃ 六速读数（1.4.2 系列）。
  用 log-linear 拟合 ln(PV)、ln(YP) 对 T 外推到 100/130℃；拟合参数打印供审计。
  ⚠️ 100/130℃ 远超实测区间；水泥浆仅 93℃ 单点，无多温度数据 ⇒ 温度档只动钻井液（被顶替相）。

口径与 runners/hu101_tailpipe.py 逐位一致。一次性探针，不改生产默认值。
"""

from __future__ import annotations

import csv
import json
import time
from dataclasses import replace
from pathlib import Path
from typing import Any, cast

import numpy as np

from cemdisp.data.loaders import load_hu101_tailpipe
from cemdisp.data.well_spec import DepthValuePoint
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.runners.hu101_tailpipe import annulus_stop_time_s
from cemdisp.transport1d import CasingFlowSolver

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PROJECT_ROOT / "results" / "呼101_居中度与温度口径场景_2026-09-11"
NZ = 250
CBL_WINDOW = "CBL评价井段(单层套管可评价段)"
CASING_KW = dict(enable_gravity=True, mixing_contact_time=True, plug_face_zero_mixing=True, has_plug=True)

# 现场钻井液多温度六速读数（1.4.2 系列，2011111.doc）：θ600/θ300
# PV = θ600 − θ300 (mPa·s)；YP = θ300 − PV (lb/100ft²)，1 lb/100ft² = 0.4788 Pa
_MUD_SERIES = {
    40: (270, 149), 50: (218, 126), 60: (175, 102), 70: (140, 82), 80: (118, 70),
}
TARGET_MEAN_SO = 0.80


def _fit_loglinear(temps: np.ndarray, values: np.ndarray) -> tuple[float, float]:
    """ln(value) = a + b·T 的最小二乘拟合。"""
    a, b = np.polyfit(temps, np.log(values), 1)
    return float(a), float(b)


def mud_rheology_at(temp_c: float) -> tuple[float, float, dict[str, Any]]:
    """把钻井液流变外推到给定温度，返回 (PV_Pa_s, YP_Pa, 审计信息)。"""
    t = np.array(sorted(_MUD_SERIES), dtype=float)
    pv = np.array([_MUD_SERIES[int(x)][0] - _MUD_SERIES[int(x)][1] for x in t], dtype=float)
    yp = np.array([_MUD_SERIES[int(x)][1] - pv[i] for i, x in enumerate(t)], dtype=float)
    a_pv, b_pv = _fit_loglinear(t, pv)
    a_yp, b_yp = _fit_loglinear(t, yp)
    pv_t = float(np.exp(a_pv + b_pv * temp_c))
    yp_lbf = float(np.exp(a_yp + b_yp * temp_c))
    audit: dict[str, Any] = {
        "measured_temps_c": t.tolist(),
        "measured_pv_mPas": pv.tolist(),
        "measured_yp_lb100ft2": yp.tolist(),
        "fit_ln_pv": {"a": a_pv, "b": b_pv},
        "fit_ln_yp": {"a": a_yp, "b": b_yp},
        "PV_mPas_at_T": pv_t,
        "YP_lb100ft2_at_T": yp_lbf,
        "YP_Pa_at_T": yp_lbf * 0.4788,
        "extrapolated_beyond_measured_range": bool(temp_c > 80.0),
    }
    return pv_t / 1000.0, yp_lbf * 0.4788, audit


def scale_profile_to_mean(profile, target: float, cap: float = 0.98):
    """保形缩放剖面，使其 200 点插值均值为 target（上限 cap）。"""
    mds = np.array([p.depth_md_m for p in profile], dtype=float)
    vals = np.array([p.value for p in profile], dtype=float)
    grid = np.linspace(mds.min(), mds.max(), 200)
    cur = float(np.interp(grid, mds, vals).mean())
    scaled = np.clip(vals * (target / cur), 0.05, cap)
    return tuple(DepthValuePoint(float(m), float(v)) for m, v in zip(mds, scaled)), cur


def uniform(well, so: float):
    return (DepthValuePoint(float(well.top_md_m), so),
            DepthValuePoint(float(well.bottom_md_m), so))


def mean_so(profile) -> float:
    mds = np.array([p.depth_md_m for p in profile], dtype=float)
    vals = np.array([p.value for p in profile], dtype=float)
    return float(np.interp(np.linspace(mds.min(), mds.max(), 200), mds, vals).mean())


def narrow_stats(geom: dict, cement: np.ndarray) -> dict[str, float]:
    b = np.asarray(geom["b"], float)
    s = np.asarray(geom["s"], float)
    y = np.asarray(geom["y"], float)
    c = np.asarray(cement, float)
    half = c.shape[0] // 2
    tot = float(np.trapezoid(np.trapezoid(b, x=s, axis=1), x=y, axis=0))
    idx = np.where(c[-1] >= 0.5)[0]
    tn = float(np.trapezoid(np.trapezoid(b[half:], x=s, axis=1), x=y[half:], axis=0))
    return {
        "starved_vol_frac": float(np.trapezoid(np.trapezoid(b * (c < 0.5), x=s, axis=1), x=y, axis=0)) / tot,
        "eta_narrow_half": float(np.trapezoid(np.trapezoid(b[half:] * c[half:], x=s, axis=1),
                                              x=y[half:], axis=0)) / tn,
        "front_narrow_m": float(s[idx.max()]) if idx.size else 0.0,
        "front_span_m": float(s[-1]),
    }


def run_case(*, label, grade, well_spec, fluids, schedule, cache=None) -> dict[str, Any]:
    key = (f"{well_spec.well_name}|{hash(tuple(round(p.value, 3) for p in well_spec.standoff_profile))}|"
           f"{hash(tuple((f.name, f.density_kg_m3, f.plastic_viscosity_pa_s) for f in fluids))}")
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
    mud = [f for f in fluids if f.role.name == "MUD"][0]
    row: dict[str, Any] = {
        "case": label,
        "grade": grade,
        "mean_standoff": round(mean_so(well_spec.standoff_profile), 4),
        "mud_PV_mPas": round((mud.plastic_viscosity_pa_s or 0) * 1000, 1),
        "mud_YP_Pa": round(mud.yield_stress_pa or 0, 2),
        "eta_E_well": float(final["全井段最终有效顶替效率"]),
        "eta_N_well": float(final["窄四分位效率"]),
        "eta_E_cbl": float(cbl.get("eta_E", float("nan"))),
        "elapsed_s": round(elapsed, 1),
    }
    row.update(narrow_stats(result.geom, np.asarray(result.cement_field, float)))
    print(f"  [{label:<30s}] 均SO={row['mean_standoff']:.3f} PV={row['mud_PV_mPas']:.0f} "
          f"YP={row['mud_YP_Pa']:.1f} ηE={row['eta_E_well']:.4f} ηE_CBL={row['eta_E_cbl']:.4f} "
          f"窄前缘={row['front_narrow_m']:.0f}m ({elapsed:.0f}s)", flush=True)
    return row


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    well, fluids, schedule, _ = load_hu101_tailpipe()
    cache: dict[str, Any] = {}
    rows: list[dict[str, Any]] = []

    print("=" * 108)
    print("  呼101 居中度口径 × 温度口径 全场景")
    print("=" * 108)

    rh_audit: dict[str, Any] = {}
    for T in (100.0, 130.0):
        pv, yp, audit = mud_rheology_at(T)
        rh_audit[f"{T:.0f}C"] = audit
        print(f"  钻井液外推 {T:.0f}℃: PV={pv*1000:.1f} mPa·s, YP={yp:.2f} Pa"
              f"（实测区间 40–80℃，外推={audit['extrapolated_beyond_measured_range']}）")

    def with_mud_temp(fl, T):
        pv, yp, _ = mud_rheology_at(T)
        return tuple(replace(f, plastic_viscosity_pa_s=pv, yield_stress_pa=yp)
                     if f.role.name == "MUD" else f for f in fl)

    # ---- A 居中度口径 ----
    print("\n[A] 居中度口径")
    for label, grade, kw in (
        ("A1_legacy假设(0.38-0.48)", "基线(model_assumption)", {}),
        ("A2_检测图-扶正器间", "可辩护(设计软件)", {"measured_standoff": "between_centralizers"}),
        ("A3_检测图-扶正器处", "可辩护(设计软件)", {"measured_standoff": "at_centralizers"}),
    ):
        w, fl, sch, _ = load_hu101_tailpipe(**kw)
        rows.append(run_case(label=label, grade=grade, well_spec=w, fluids=fl, schedule=sch, cache=cache))

    rows.append(run_case(label="A4_均匀SO=0.80", grade="[反推场景]", well_spec=replace(
        well, standoff_profile=uniform(well, TARGET_MEAN_SO), standoff_measured=True),
        fluids=fluids, schedule=schedule, cache=cache))

    for src, tag in (("at_centralizers", "扶正器处"), ("between_centralizers", "扶正器间")):
        w_src, _, _, _ = load_hu101_tailpipe(measured_standoff=src)
        prof, cur = scale_profile_to_mean(w_src.standoff_profile, TARGET_MEAN_SO)
        rows.append(run_case(
            label=f"A5_{tag}保形到均0.80", grade="[反推场景]",
            well_spec=replace(well, standoff_profile=prof, standoff_measured=True),
            fluids=fluids, schedule=schedule, cache=cache))
        print(f"     （{tag}剖面原均 {cur:.3f} → 缩放后均 {mean_so(prof):.3f}）")

    # ---- B 温度口径（legacy 居中度）----
    print("\n[B] 温度口径（钻井液，legacy 居中度）")
    for T in (65.0, 100.0, 130.0):
        fl = fluids if T == 65.0 else with_mud_temp(fluids, T)
        grade = "基线(现场化验65℃)" if T == 65.0 else f"可辩护(外推{T:.0f}℃)"
        rows.append(run_case(label=f"B_钻井液@{T:.0f}℃", grade=grade, well_spec=well,
                             fluids=fl, schedule=schedule, cache=cache))

    # ---- C 交叉：均0.80 × 温度 ----
    print("\n[C] 交叉（均SO=0.80 × 温度口径）")
    w80 = replace(well, standoff_profile=uniform(well, TARGET_MEAN_SO), standoff_measured=True)
    for T in (100.0, 130.0):
        rows.append(run_case(label=f"C_均0.80+钻井液@{T:.0f}℃", grade="[反推场景]", well_spec=w80,
                             fluids=with_mud_temp(fluids, T), schedule=schedule, cache=cache))

    # ---- 汇总 ----
    base = rows[0]["eta_E_well"]
    print("\n" + "=" * 108)
    print(f"  {'case':<30s}{'性质':<22s}{'均SO':>7s}{'ηE全井':>9s}{'Δ':>9s}{'ηE_CBL':>9s}{'饥饿份额':>10s}{'窄前缘m':>9s}")
    for r in rows:
        print(f"  {r['case']:<30s}{r['grade']:<22s}{r['mean_standoff']:>7.3f}{r['eta_E_well']:>9.4f}"
              f"{r['eta_E_well']-base:>+9.4f}{r['eta_E_cbl']:>9.4f}{r['starved_vol_frac']:>10.3f}"
              f"{r['front_narrow_m']:>9.0f}")

    fields = ["case", "grade", "mean_standoff", "mud_PV_mPas", "mud_YP_Pa", "eta_E_well",
              "eta_N_well", "eta_E_cbl", "starved_vol_frac", "eta_narrow_half",
              "front_narrow_m", "front_span_m", "elapsed_s"]
    with (OUTPUT_DIR / "场景对比表.csv").open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})
    (OUTPUT_DIR / "场景摘要.json").write_text(json.dumps({
        "description": "呼101 居中度口径 × 温度口径 全场景",
        "warning": "grade 标 [反推场景] 的行不得作为验证数字引用",
        "rheology_extrapolation_audit": rh_audit,
        "cases": rows,
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n结果已落盘：{OUTPUT_DIR}")


if __name__ == "__main__":
    main()
