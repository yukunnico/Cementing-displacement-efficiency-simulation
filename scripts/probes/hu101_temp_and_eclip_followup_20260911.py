"""呼101 温度口径 + e_clip 分叉 补测（spike，2026-09-11 第六轮续）。

前一轮（hu101_standoff_and_temperature_scenarios_20260911.py）A 组有效、B/C 组因
`np.polyfit` 返回值顺序写反（polyfit 返回 [slope, intercept]）导致温度外推指数溢出而作废。
本脚本修正拟合方向并补测，同时加做用户裁定的 e_clip 分叉对照。

修正后钻井液外推（log-linear 拟合 40–80℃ 实测 ln(PV)/ln(YP)）：
  斜率 ≈ −0.0231/℃（PV）、−0.0061/℃（YP），100℃/130℃ 落在合理量级。
  ⚠️ 仍为超出实测区间的外推。

case：
  A'/A3b 检测图剖面 + e_clip 保持 0.55（对照 A2/A3 的放开 0.90）
  B100/B130 钻井液外推 100/130℃（legacy 居中度）
  C100/C130 均匀 SO=0.80 + 钻井液 100/130℃
A 组其余结果取自前一轮控制台输出（已核对有效）。一次性探针，不改生产默认值。
"""

from __future__ import annotations

import csv
import json
import math
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

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = PROJECT_ROOT / "results" / "呼101_居中度与温度口径场景_2026-09-11"
NZ = 250
CBL_WINDOW = "CBL评价井段(单层套管可评价段)"
CASING_KW = dict(enable_gravity=True, mixing_contact_time=True, plug_face_zero_mixing=True, has_plug=True)
_MUD_SERIES = {40: (270, 149), 50: (218, 126), 60: (175, 102), 70: (140, 82), 80: (118, 70)}
TARGET_MEAN_SO = 0.80

# A 组（前一轮有效结果，直接并入汇总，不重跑）
A_ROWS: list[dict[str, Any]] = [
    {"case": "A1_legacy假设(0.38-0.48)", "grade": "基线(model_assumption)", "mean_standoff": 0.429,
     "eta_E_well": 0.4919, "eta_E_cbl": 0.5035, "front_narrow_m": 50.0, "source": "前一轮"},
    {"case": "A2_检测图-扶正器间(eclip0.90)", "grade": "可辩护(设计软件)", "mean_standoff": 0.515,
     "eta_E_well": 0.4814, "eta_E_cbl": 0.4941, "front_narrow_m": 40.0, "source": "前一轮"},
    {"case": "A3_检测图-扶正器处(eclip0.90)", "grade": "可辩护(设计软件)", "mean_standoff": 0.734,
     "eta_E_well": 0.5371, "eta_E_cbl": 0.5539, "front_narrow_m": 119.0, "source": "前一轮"},
    {"case": "A4_均匀SO=0.80", "grade": "[反推场景]", "mean_standoff": 0.800,
     "eta_E_well": 0.6805, "eta_E_cbl": 0.7183, "front_narrow_m": 1130.0, "source": "前一轮"},
    {"case": "A5_扶正器处保形到均0.80", "grade": "[反推场景]", "mean_standoff": 0.800,
     "eta_E_well": 0.6333, "eta_E_cbl": 0.6534, "front_narrow_m": 337.0, "source": "前一轮"},
    {"case": "A5b_扶正器间保形到均0.765(上限卡住)", "grade": "[反推场景]", "mean_standoff": 0.765,
     "eta_E_well": 0.5349, "eta_E_cbl": 0.5341, "front_narrow_m": 40.0, "source": "前一轮"},
]


def _fit_loglinear(temps: np.ndarray, values: np.ndarray) -> tuple[float, float]:
    """返回 (intercept, slope)：ln(value) = intercept + slope·T。

    ⚠️ np.polyfit 返回 [slope, intercept]（降幂），此处显式换位。
    """
    slope, intercept = np.polyfit(temps, np.log(values), 1)
    return float(intercept), float(slope)


def mud_rheology_at(temp_c: float) -> tuple[float, float, dict[str, Any]]:
    t = np.array(sorted(_MUD_SERIES), dtype=float)
    pv = np.array([_MUD_SERIES[int(x)][0] - _MUD_SERIES[int(x)][1] for x in t], dtype=float)
    yp = np.array([_MUD_SERIES[int(x)][1] - pv[i] for i, x in enumerate(t)], dtype=float)
    i_pv, s_pv = _fit_loglinear(t, pv)
    i_yp, s_yp = _fit_loglinear(t, yp)
    pv_t = float(np.exp(i_pv + s_pv * temp_c))
    yp_lbf = float(np.exp(i_yp + s_yp * temp_c))
    audit = {
        "measured_temps_c": t.tolist(), "measured_pv_mPas": pv.tolist(),
        "measured_yp_lb100ft2": yp.tolist(),
        "fit": {"pv_intercept": i_pv, "pv_slope_per_C": s_pv,
                "yp_intercept": i_yp, "yp_slope_per_C": s_yp},
        "PV_mPas_at_T": pv_t, "YP_Pa_at_T": yp_lbf * 0.4788,
        "extrapolated_beyond_measured_range": bool(temp_c > 80.0),
    }
    if not (math.isfinite(pv_t) and pv_t > 0 and math.isfinite(yp_lbf) and yp_lbf > 0):
        raise ValueError(f"温度外推得到非物理值 T={temp_c}: PV={pv_t}, YP={yp_lbf}")
    return pv_t / 1000.0, yp_lbf * 0.4788, audit


def mean_so(profile) -> float:
    mds = np.array([p.depth_md_m for p in profile], dtype=float)
    vals = np.array([p.value for p in profile], dtype=float)
    return float(np.interp(np.linspace(mds.min(), mds.max(), 200), mds, vals).mean())


def uniform(well, so: float):
    return (DepthValuePoint(float(well.top_md_m), so), DepthValuePoint(float(well.bottom_md_m), so))


def narrow_stats(geom: dict, cement: np.ndarray) -> dict[str, float]:
    b = np.asarray(geom["b"], float); s = np.asarray(geom["s"], float)
    y = np.asarray(geom["y"], float); c = np.asarray(cement, float)
    tot = float(np.trapezoid(np.trapezoid(b, x=s, axis=1), x=y, axis=0))
    idx = np.where(c[-1] >= 0.5)[0]
    return {
        "starved_vol_frac": float(np.trapezoid(np.trapezoid(b * (c < 0.5), x=s, axis=1), x=y, axis=0)) / tot,
        "front_narrow_m": float(s[idx.max()]) if idx.size else 0.0,
    }


def run_case(*, label, grade, well_spec, fluids, schedule, cache=None) -> dict[str, Any]:
    key = (f"{well_spec.well_name}|{hash(tuple(round(p.value, 3) for p in well_spec.standoff_profile))}|"
           f"{hash(tuple((f.name, f.plastic_viscosity_pa_s) for f in fluids))}")
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
        "case": label, "grade": grade, "mean_standoff": round(mean_so(well_spec.standoff_profile), 4),
        "mud_PV_mPas": round((mud.plastic_viscosity_pa_s or 0) * 1000, 2),
        "mud_YP_Pa": round(mud.yield_stress_pa or 0, 2),
        "eta_E_well": float(final["全井段最终有效顶替效率"]),
        "eta_E_cbl": float(cbl.get("eta_E", float("nan"))),
        "elapsed_s": round(elapsed, 1), "source": "本轮",
    }
    row.update(narrow_stats(result.geom, np.asarray(result.cement_field, float)))
    print(f"  [{label:<34s}] 均SO={row['mean_standoff']:.3f} PV={row['mud_PV_mPas']:.1f}mPa·s "
          f"YP={row['mud_YP_Pa']:.2f}Pa ηE={row['eta_E_well']:.4f} ηE_CBL={row['eta_E_cbl']:.4f} "
          f"窄前缘={row['front_narrow_m']:.0f}m ({elapsed:.0f}s)", flush=True)
    return row


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    well, fluids, schedule, _ = load_hu101_tailpipe()
    cache: dict[str, Any] = {}
    new_rows: list[dict[str, Any]] = []

    print("=" * 112)
    print("  呼101 温度口径 + e_clip 分叉 补测")
    print("=" * 112)
    rh_audit: dict[str, Any] = {}
    for T in (100.0, 130.0):
        pv, yp, audit = mud_rheology_at(T)
        rh_audit[f"{T:.0f}C"] = audit
        print(f"  钻井液外推 {T:.0f}℃: PV={pv*1000:.1f} mPa·s, YP={yp:.2f} Pa "
              f"(拟合斜率 {audit['fit']['pv_slope_per_C']:.5f}/℃, {audit['fit']['yp_slope_per_C']:.5f}/℃)")

    def with_mud_temp(fl, T):
        pv, yp, _ = mud_rheology_at(T)
        return tuple(replace(f, plastic_viscosity_pa_s=pv, yield_stress_pa=yp)
                     if f.role.name == "MUD" else f for f in fl)

    # ---- A': e_clip 保持 0.55 的检测图口径 ----
    print("\n[A'] e_clip 分叉（检测图剖面，保持 e_clip=0.55）")
    for src, tag in (("between_centralizers", "扶正器间"), ("at_centralizers", "扶正器处")):
        w_src, _, _, _ = load_hu101_tailpipe(measured_standoff=src)
        new_rows.append(run_case(
            label=f"A'_检测图-{tag}_eclip0.55", grade="口径裁定对照",
            well_spec=replace(well, standoff_profile=w_src.standoff_profile, standoff_measured=False),
            fluids=fluids, schedule=schedule, cache=cache))

    # ---- B：温度口径 ----
    print("\n[B] 温度口径（钻井液，legacy 居中度）")
    for T in (100.0, 130.0):
        new_rows.append(run_case(label=f"B_钻井液@{T:.0f}C", grade=f"可辩护(外推{T:.0f}℃)",
                                 well_spec=well, fluids=with_mud_temp(fluids, T),
                                 schedule=schedule, cache=cache))

    # ---- C：均0.80 × 温度 ----
    print("\n[C] 交叉（均匀 SO=0.80 × 温度口径）")
    w80 = replace(well, standoff_profile=uniform(well, TARGET_MEAN_SO), standoff_measured=True)
    for T in (100.0, 130.0):
        new_rows.append(run_case(label=f"C_均0.80+钻井液@{T:.0f}C", grade="[反推场景]",
                                 well_spec=w80, fluids=with_mud_temp(fluids, T),
                                 schedule=schedule, cache=cache))

    rows = A_ROWS + new_rows
    base = rows[0]["eta_E_well"]
    print("\n" + "=" * 112)
    print(f"  {'case':<36s}{'性质':<20s}{'均SO':>7s}{'ηE全井':>9s}{'Δ':>9s}{'ηE_CBL':>9s}{'窄前缘m':>9s}")
    for r in rows:
        print(f"  {r['case']:<36s}{r['grade']:<20s}{r['mean_standoff']:>7.3f}{r['eta_E_well']:>9.4f}"
              f"{r['eta_E_well']-base:>+9.4f}{r['eta_E_cbl']:>9.4f}"
              f"{r.get('front_narrow_m', float('nan')):>9.0f}")

    fields = ["case", "grade", "mean_standoff", "mud_PV_mPas", "mud_YP_Pa", "eta_E_well",
              "eta_E_cbl", "starved_vol_frac", "front_narrow_m", "elapsed_s", "source"]
    with (OUTPUT_DIR / "场景对比表.csv").open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    (OUTPUT_DIR / "场景摘要.json").write_text(json.dumps({
        "description": "呼101 居中度口径 × 温度口径 全场景",
        "warning": "grade 标 [反推场景] 的行不得作为验证数字引用",
        "rheology_extrapolation_audit": rh_audit,
        "cases": rows,
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n结果已落盘：{OUTPUT_DIR}")


if __name__ == "__main__":
    main()
