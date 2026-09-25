"""P1 批次：打开流函数路径上被丢弃的屈服门（enable_stream_yield_gate=True）。

背景（2026-09-19 取证）
----------------------
默认生产口径 ``CORRECTED_KW = {enable_yield_gate=True, enable_regime_split=True,
enable_local_i3=True}`` 里 ``enable_yield_gate=True`` 会在 ``annulus_d2dga.py``
每步算出连续冻结度 ``wall``，但 ``_compute_velocity`` 在 ``:1829`` 用
``wall=(wall if self.enable_stream_yield_gate else None)`` 把它**丢掉**——新路径
（``enable_stream_function=True``，默认）的速度场因此**零屈服物理**。
``enable_stream_yield_gate`` 默认 False，且构造函数的"静默无效开关"守卫
（``:561-576``）只检查反方向组合，**不告警**。本批次 = P1。

验收判据（用户 2026-09-19 指令）
--------------------------------
① 八井全井 η_E/η_N 相对 G3 基线（``results/源模型口径重跑_2026-09-14``）的位移；
② **居中度灵敏度不得失灵**——呼101（强偏心）+ 呼103（居中好）各 4 档 standoff
   位移，与 ``results/敏感性变体_当前口径_2026-09-16``（同变体、P1 关）对照。

口径纪律
--------
- 求解器配置唯一来源 = ``rerun_all_wells_corrected.CORRECTED_KW``，本脚本只**叠加**
  一个开关，不复制常量（防漂移）。
- 输出写到**新建目录** ``results/P1屈服门_2026-09-19/``，不触碰权威目录与既有结果目录。
- 断点续跑：已存在的 ``*_结果摘要.json`` 直接复用。

用法
----
    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -u \
        scripts/entrypoints/run_p1_yield_gate_20260919.py [--group wells|standoff|all]
"""
from __future__ import annotations

import argparse
import csv
import importlib
import json
import sys
import time
from pathlib import Path

from cemdisp.data.fluid_spec import FluidRole
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
from entrypoints.run_sensitivity_current_20260916 import shift_standoff  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = PROJECT_ROOT / "results" / "P1屈服门_2026-09-19"
# P1 关的对照基线：八井用 G3 重跑；standoff 档用 09-16 敏感性目录（同变体）。
G3_DIR = PROJECT_ROOT / "results" / "源模型口径重跑_2026-09-14" / "单井结果"
SENS_DIR = PROJECT_ROOT / "results" / "敏感性变体_当前口径_2026-09-16"

# P1：唯一改动——把每步算出的 wall 真正送进流函数算子。
P1_KW = dict(CORRECTED_KW, enable_stream_yield_gate=True)

WELLS: dict[str, tuple[str, str, str]] = {
    "hu101": ("cemdisp.data.loaders.hu101_loader", "load_hu101_tailpipe", "hu101_corrected_on.json"),
    "hu102": ("cemdisp.data.loaders.hu102_loader", "load_hu102_tailpipe", "hu102_corrected_on.json"),
    "hu103": ("cemdisp.data.loaders.hu103_loader", "load_hu103_tailpipe", "hu103_corrected_on.json"),
    "hu1": ("cemdisp.data.loaders.hu1_loader", "load_hu1_tailpipe", "hu1_corrected_on.json"),
    "hu2": ("cemdisp.data.loaders.hu2_loader", "load_hu2_tailpipe", "hu2_corrected_on.json"),
    "ht1_001": ("cemdisp.data.loaders.ht1_001_loader", "load_ht1_001_tailpipe", "ht1_001_corrected_on.json"),
    "ht1_003": ("cemdisp.data.loaders.ht1_003_loader", "load_ht1_003_tailpipe", "ht1_003_corrected_on.json"),
    "ht1_004": ("cemdisp.data.loaders.ht1_004_loader", "load_ht1_004_tailpipe", "ht1_004_corrected_on.json"),
}

# standoff 档：与 09-16 目录中的变体名逐字一致，便于横向对照。
STANDOFF_WELLS = {
    "呼101": ("cemdisp.data.loaders.hu101_loader", "load_hu101_tailpipe"),
    "呼103": ("cemdisp.data.loaders.hu103_loader", "load_hu103_tailpipe"),
}
STANDOFF_DELTAS = [(-0.30, "standoff_m0.30"), (-0.10, "standoff_m0.1"),
                   (+0.10, "standoff_p0.1"), (+0.30, "standoff_p0.30")]

CEMENT_ROLES = {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}


def _run_one(loader, well_fn=None) -> tuple[dict, dict, float]:
    """一次完整流水线（1D 套管 + 环空二维），P1 口径。返回 (最终结果, summary, 耗时s)。"""
    well, fluids, schedule, _ = loader()
    if well_fn is not None:
        well = well_fn(well)

    casing = CasingFlowSolver(
        enable_gravity=True,
        mixing_contact_time=True,
        plug_face_zero_mixing=True,
        has_plug=True,
    )
    cr = casing.run(well, fluids, schedule)
    inlet = build_coupled_annulus_inlet_provider(
        cr, casing, fluids, split_cement_phases=True,
    )
    tt = min(_total_t(schedule) + 1200.0, _stop_t(cr, fluids))
    solver = AnnulusD2DGASolver(total_t=tt, nz=NZ, enable_cfl_adaptive=True, **P1_KW)
    t0 = time.perf_counter()
    res = solver.run(well, fluids, inlet, schedule=schedule)
    return res.summary["最终结果"], res.summary, round(time.perf_counter() - t0, 1)


def _cached(case_json: Path, fn) -> tuple[dict, float | None]:
    if case_json.exists():
        summary = json.loads(case_json.read_text(encoding="utf-8"))
        return summary["最终结果"], None
    final, summary, elapsed = fn()
    case_json.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return final, elapsed


def _row(well_key, var, final, base_e, base_n, elapsed):
    eta_e = float(final["全井段最终有效顶替效率"])
    eta_n = float(final["窄四分位效率"])
    return {
        "井名": well_key, "变体": var, "η_E": eta_e, "η_N": eta_n,
        "窜槽": float(final["最终窜槽指数"]), "混浆": float(final["最终混浆指数"]),
        "失稳": float(final["最终失稳指数"]),
        "Δη_E_pp": (eta_e - base_e) * 100.0, "Δη_N_pp": (eta_n - base_n) * 100.0,
        "基线η_E": base_e, "基线η_N": base_n,
        "耗时_s": elapsed if elapsed is not None else "",
    }


def run_wells() -> list[dict]:
    rows = []
    for well_key, (mod, fn, base_file) in WELLS.items():
        base = json.loads((G3_DIR / base_file).read_text(encoding="utf-8"))
        be, bn = float(base["eta_E"]), float(base["eta_N"])
        loader = getattr(importlib.import_module(mod), fn)
        print(f"\n=== [P1·井] {well_key}  对照基线 G3: η_E={be:.4f} η_N={bn:.4f} ===", flush=True)
        cj = OUT_DIR / f"{well_key}_P1_结果摘要.json"
        final, elapsed = _cached(cj, lambda: _run_one(loader))
        r = _row(well_key, "P1", final, be, bn, elapsed)
        print(f"  [{'复用' if elapsed is None else '计算'}] {well_key}: "
              f"η_E={r['η_E']:.4f} ({r['Δη_E_pp']:+.2f}pp) "
              f"η_N={r['η_N']:.4f} ({r['Δη_N_pp']:+.2f}pp)"
              + (f" ({elapsed}s)" if elapsed else ""), flush=True)
        rows.append(r)
    return rows


def run_standoff() -> list[dict]:
    rows = []
    for well_key, (mod, fn) in STANDOFF_WELLS.items():
        loader = getattr(importlib.import_module(mod), fn)
        for delta, var in STANDOFF_DELTAS:
            ref = SENS_DIR / f"{well_key}_{var}_结果摘要.json"
            if not ref.exists():
                raise FileNotFoundError(f"对照基线缺失：{ref}")
            rs = json.loads(ref.read_text(encoding="utf-8"))
            rf = rs["最终结果"] if "最终结果" in rs else rs
            be = float(rf["全井段最终有效顶替效率"])
            bn = float(rf["窄四分位效率"])
            print(f"\n=== [P1·standoff] {well_key} × {var}  "
                  f"对照(P1关) η_E={be:.4f} η_N={bn:.4f} ===", flush=True)
            cj = OUT_DIR / f"{well_key}_{var}_P1_结果摘要.json"
            final, elapsed = _cached(
                cj, lambda: _run_one(loader, lambda w: shift_standoff(w, delta)))
            r = _row(well_key, var, final, be, bn, elapsed)
            print(f"  [{'复用' if elapsed is None else '计算'}] {well_key}×{var}: "
                  f"η_E={r['η_E']:.4f} ({r['Δη_E_pp']:+.2f}pp) "
                  f"η_N={r['η_N']:.4f} ({r['Δη_N_pp']:+.2f}pp)"
                  + (f" ({elapsed}s)" if elapsed else ""), flush=True)
            rows.append(r)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", choices=["wells", "standoff", "all"], default="all")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("P1 口径 = CORRECTED_KW + enable_stream_yield_gate=True\n"
          f"CORRECTED_KW = {CORRECTED_KW}\n输出 = {OUT_DIR}", flush=True)

    rows: list[dict] = []
    if args.group in ("wells", "all"):
        rows += run_wells()
    if args.group in ("standoff", "all"):
        rows += run_standoff()

    csv_path = OUT_DIR / f"汇总表_P1_{args.group}.csv"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\n完成：{len(rows)} 行 → {csv_path}", flush=True)

    # 验收判据②：居中度灵敏度（相对各井自身 P1 关基线）
    print("\n=== 居中度灵敏度（P1 口径，Δ相对该井自身基线）===", flush=True)
    for well_key in STANDOFF_WELLS:
        sub = [r for r in rows if r["井名"] == well_key]
        if not sub:
            continue
        span_n = max(r["Δη_N_pp"] for r in sub) - min(r["Δη_N_pp"] for r in sub)
        span_e = max(r["Δη_E_pp"] for r in sub) - min(r["Δη_E_pp"] for r in sub)
        print(f"  {well_key}: standoff -0.30→+0.30 跨度  "
              f"Δη_N={span_n:.2f}pp  Δη_E={span_e:.2f}pp", flush=True)


if __name__ == "__main__":
    main()
