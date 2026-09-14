"""居中度→顶替效率 响应强度调研（spike，2026-09-11 第五轮）。

用户诉求：模型环空段怎么调，能把"居中度对顶替效率的影响"调上去。
物理直觉（用户给出）：居中度低 → 流体走宽边多、窄边少。

先诊断后调参。先验假设：低居中度区间响应不是"弱"而是"饱和"——
窄边已 100% 空，且体积加权 η_E 在该区间反而奖励更高偏心度
（饥饿的窄间隙体积份额随 e 增大而减小）。

A 段（响应曲线形状 + 分解）：
  均匀 SO ∈ {0.35,0.45,0.55,0.65,0.75}，每 case 拆出：
    eta_E（体积加权）、饥饿体积份额 ∫b·1[c<0.5]/∫b、宽半区/窄半区效率、窄边前缘位置。
  判据：若 eta_E ≈ 1 − 饥饿体积份额，则平台是纯体积加权假象。

B 段（旋钮对响应斜率的控制）：5 个候选旋钮 × SO ∈ {0.40,0.55,0.70}
  响应斜率定义 Δ = η_E(0.70) − η_E(0.40)，斜率越大 = 居中度影响越强。
    b1 e_clip_max=0.90（放开截断，暴露真实 e）
    b2 enable_power_law_gap_law=False（通量指数 2.20 → 2.00）
    b3 yield_gate_f_safety=1.6（更强的壁面冻结门槛）
    b4 yield_regularization_M=1000（屈服应力在低剪切区更陡）
    b5 dispersion_azimuthal=0.002（方位弥散≈关，锋面更锐）

口径与 runners/hu101_tailpipe.py 逐位一致。一次性探针，不改生产默认值。
"""

from __future__ import annotations

import csv
import json
import time
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np

from cemdisp.data.loaders import load_hu101_tailpipe
from cemdisp.data.well_spec import DepthValuePoint
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.runners.hu101_tailpipe import annulus_stop_time_s
from cemdisp.transport1d import CasingFlowSolver

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PROJECT_ROOT / "results" / "居中度响应强度调研_2026-09-11"
NZ = 250
CASING_KW = dict(enable_gravity=True, mixing_contact_time=True, plug_face_zero_mixing=True, has_plug=True)
SO_SWEEP = (0.35, 0.45, 0.55, 0.65, 0.75)
SO_KNOB = (0.40, 0.55, 0.70)
KNOBS: dict[str, dict[str, Any]] = {
    "b1_eclip0.90": {"e_clip_max": 0.90, "enable_e_clip_ruling": False},
    "b2_gaplaw_off": {"enable_power_law_gap_law": False},
    "b3_fsafety1.6": {"yield_gate_f_safety": 1.6},
    "b4_M1000": {"yield_regularization_M": 1000.0},
    "b5_dispaz0.002": {"dispersion_azimuthal": 0.002},
}


def uniform(well, so: float):
    return (DepthValuePoint(float(well.top_md_m), so),
            DepthValuePoint(float(well.bottom_md_m), so))


def trapz2d(arr: np.ndarray, s: np.ndarray, y: np.ndarray) -> float:
    return float(np.trapezoid(np.trapezoid(arr, x=s, axis=1), x=y, axis=0))


def decompose(geom: dict, cement: np.ndarray) -> dict[str, float]:
    """按体积权重拆解浓度场：整体/宽半区/窄半区效率 + 饥饿体积份额 + 窄边前缘。"""
    b = np.asarray(geom["b"], dtype=float)
    s = np.asarray(geom["s"], dtype=float)
    y = np.asarray(geom["y"], dtype=float)
    c = np.asarray(cement, dtype=float)
    tot = trapz2d(b, s, y)
    half = c.shape[0] // 2
    tw = trapz2d(b[:half], s, y[:half])
    tn = trapz2d(b[half:], s, y[half:])
    line = c[-1]                                  # 最窄一行（窄边）
    idx = np.where(line >= 0.5)[0]
    return {
        "eta_E": trapz2d(b * c, s, y) / tot,
        "eta_wide_half": trapz2d(b[:half] * c[:half], s, y[:half]) / tw,
        "eta_narrow_half": trapz2d(b[half:] * c[half:], s, y[half:]) / tn,
        "starved_vol_frac": trapz2d(b * (c < 0.5), s, y) / tot,
        "wide_vol_frac": tw / tot,
        "front_narrow_m": float(s[idx.max()]) if idx.size else 0.0,
        "front_span_m": float(s[-1]),
    }


def run_case(*, label, well_spec, fluids, schedule, solver_kw=None, cache=None) -> dict[str, Any]:
    solver_kw = dict(solver_kw or {})
    key = f"{well_spec.well_name}|{hash(tuple(round(p.value, 3) for p in well_spec.standoff_profile))}"
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
    result = AnnulusD2DGASolver(total_t=stop, nz=NZ, **solver_kw).run(
        well_spec, fluids, inlet, schedule=schedule)
    elapsed = time.perf_counter() - t0

    row: dict[str, Any] = {"case": label}
    row.update(decompose(result.geom, np.asarray(result.cement_field, dtype=float)))
    row["elapsed_s"] = round(elapsed, 1)
    print(f"  [{label:<26s}] ηE={row['eta_E']:.4f} 饥饿份额={row['starved_vol_frac']:.3f} "
          f"1-饥饿={1-row['starved_vol_frac']:.3f} 宽半={row['eta_wide_half']:.3f} "
          f"窄半={row['eta_narrow_half']:.3f} 窄前缘={row['front_narrow_m']:.0f}m ({elapsed:.0f}s)", flush=True)
    return row


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    well, fluids, schedule, _ = load_hu101_tailpipe()
    cache: dict[str, Any] = {}
    span = float(well.bottom_md_m - well.top_md_m)

    print("=" * 104)
    print("  居中度 → 顶替效率 响应强度调研（呼101 井）")
    print("=" * 104)
    print(f"  模型域 5400–7868m，轴向前缘跨度 {span:.0f}m\n")

    # ---- A 段：响应曲线 ----
    print("[A] 均匀居中度扫描（含分解）")
    a_rows: list[dict[str, Any]] = []
    for so in SO_SWEEP:
        a_rows.append(run_case(
            label=f"A_SO={so:.2f}",
            well_spec=replace(well, standoff_profile=uniform(well, so), standoff_measured=False),
            fluids=fluids, schedule=schedule, cache=cache))
    a_rows.append(run_case(label="A_实测口径(0.38-0.48)", well_spec=well, fluids=fluids,
                           schedule=schedule, cache=cache))

    # ---- B 段：旋钮 × 居中度 ----
    print("\n[B] 旋钮对响应斜率的控制（斜率 = ηE(0.70) − ηE(0.40)）")
    b_rows: list[dict[str, Any]] = []
    for kname, kw in KNOBS.items():
        for so in SO_KNOB:
            b_rows.append(run_case(
                label=f"{kname}_SO={so:.2f}",
                well_spec=replace(well, standoff_profile=uniform(well, so), standoff_measured=False),
                fluids=fluids, schedule=schedule, solver_kw=kw, cache=cache))

    # ---- 汇总 ----
    print("\n" + "=" * 104)
    print("  [A] 响应曲线分解")
    print(f"  {'case':<24s}{'ηE':>8s}{'饥饿份额':>10s}{'1-饥饿':>9s}{'宽半':>8s}{'窄半':>8s}{'窄前缘m':>9s}")
    for r in a_rows:
        print(f"  {r['case']:<24s}{r['eta_E']:>8.4f}{r['starved_vol_frac']:>10.3f}"
              f"{1-r['starved_vol_frac']:>9.3f}{r['eta_wide_half']:>8.3f}"
              f"{r['eta_narrow_half']:>8.3f}{r['front_narrow_m']:>9.0f}")

    print("\n  [B] 旋钮效应")
    print(f"  {'knob @ SO':<30s}{'ηE':>8s}{'饥饿份额':>10s}{'宽半':>8s}{'窄半':>8s}{'窄前缘m':>9s}")
    slopes: dict[str, float] = {}
    for kname in KNOBS:
        sub = {r["case"].rsplit("_SO=", 1)[1]: r for r in b_rows if r["case"].startswith(kname + "_SO=")}
        for so_s, r in sub.items():
            print(f"  {kname + ' @ SO=' + so_s:<30s}{r['eta_E']:>8.4f}{r['starved_vol_frac']:>10.3f}"
                  f"{r['eta_wide_half']:>8.3f}{r['eta_narrow_half']:>8.3f}{r['front_narrow_m']:>9.0f}")
        if "0.40" in sub and "0.70" in sub:
            slopes[kname] = sub["0.70"]["eta_E"] - sub["0.40"]["eta_E"]
    a04 = next((r for r in a_rows if r["case"] == "A_SO=0.45"), None)
    print("\n  响应斜率对照（ηE@0.70 − ηE@0.40；基线见 A 段内插）：")
    for kname, sl in sorted(slopes.items(), key=lambda kv: -kv[1]):
        print(f"    {kname:<30s}{sl:+.4f}")

    fields = ["case", "eta_E", "eta_wide_half", "eta_narrow_half", "starved_vol_frac",
              "wide_vol_frac", "front_narrow_m", "front_span_m", "elapsed_s"]
    with (OUTPUT_DIR / "响应强度调研对比表.csv").open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in a_rows + b_rows:
            w.writerow({k: r[k] for k in fields})
    (OUTPUT_DIR / "响应强度调研摘要.json").write_text(
        json.dumps({"description": "居中度→ηE 响应强度调研", "curve": a_rows,
                    "knobs": b_rows, "slopes": slopes}, ensure_ascii=False, indent=2),
        encoding="utf-8")
    print(f"\n结果已落盘：{OUTPUT_DIR}")


if __name__ == "__main__":
    main()
