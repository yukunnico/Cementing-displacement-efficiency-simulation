# -*- coding: utf-8 -*-
"""P3.1 调试：run_ht1004_target vs sandbox/out_*.csv 逐点偏差（临时件，最后保留供审查）。"""
import sys
import time
from pathlib import Path

import numpy as np

BRANCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BRANCH))
SBX = BRANCH / "results" / "_probe_matlab靶_2026-10-07" / "sandbox"

from cemdisp.transport1d.hydraulics import run_ht1004_target  # noqa: E402

t0 = time.time()
res = run_ht1004_target(SBX)
print(f"[run] done in {time.time()-t0:.1f}s")
s = res["summary"]
print("[summary]", s)


def load(name):
    return np.loadtxt(str(SBX / name), delimiter=",")


def rep(label, mine, tgt, mode="abs"):
    d = np.abs(mine - tgt)
    if mode == "abs":
        mx = float(np.max(d))
        idx = int(np.argmax(d))
        print(f"  {label:46s} max|diff|={mx:.3e} @idx={idx}")
    else:  # rel
        rel = d / np.maximum(np.abs(tgt), 1e-300)
        mx = float(np.max(rel))
        idx = int(np.argmax(rel))
        print(f"  {label:46s} max|rel|={mx:.3e} @idx={idx}")
    if not np.all(np.isfinite(mine)):
        print("    !! mine has non-finite")


print("== out_pump_pressure_surface ==")
tgt = load("out_pump_pressure_surface.csv")
mine = res["out_pump_pressure_surface"]
rep("time", mine[:, 0], tgt[:, 0]); rep("pump_MPa", mine[:, 1], tgt[:, 1], "rel")

print("== out_annuli_bottom_pressure ==")
tgt = load("out_annuli_bottom_pressure.csv")
mine = res["out_annuli_bottom_pressure"]
rep("annuli_bottom_MPa", mine[:, 1], tgt[:, 1])

print("== out_result_matrix_volume_pressure ==")
tgt = load("out_result_matrix_volume_pressure.csv")
mine = res["out_result_matrix_volume_pressure"]
names = ["time", "cumvol_m3", "ECD_casing_bot", "pump_MPa", "annuli_bot_MPa"]
for j, nm in enumerate(names):
    rep(f"col{j} {nm}", mine[:, j], tgt[:, j], "rel" if nm == "pump_MPa" else "abs")

print("== out_ecd_casing_full (189x333) ==")
tgt = load("out_ecd_casing_full.csv")
mine = res["out_ecd_casing_full"]
rep("ECD_casing_full", mine, tgt)

print("== out_four_point_backpressure ==")
tgt = load("out_four_point_backpressure.csv")
mine = res["out_four_point_backpressure"]
for j, nm in enumerate(["time", "required", "BP_lower", "BP_upper", "conflict"]):
    rep(f"col{j} {nm}", mine[:, j], tgt[:, j])
print("  conflict bitwise equal:", bool(np.array_equal(mine[:, 4], tgt[:, 4])))

print("== out_four_point_ecd_new ==")
tgt = load("out_four_point_ecd_new.csv")
mine = res["out_four_point_ecd_new"]
for j in range(5):
    rep(f"col{j}", mine[:, j], tgt[:, j])

print("== out_window_backpressure_ctrl ==")
tgt = load("out_window_backpressure_ctrl.csv")
mine = res["out_window_backpressure_ctrl"]
for j, nm in enumerate(["time", "ECD_base", "ECD_bot_ctrl", "ECD_shoe_ctrl", "P_wh"]):
    rep(f"col{j} {nm}", mine[:, j], tgt[:, j])
tvdb = float(res["summary"]["TVD_bottom_m"])
den = 0.00981 * float(np.max(res["tvd_cum_m1b"]))
bp_tgt = (tgt[:, 2] - tgt[:, 1]) * den
rep("BP_apply_win(recast)", res["bp_apply_win_MPa"], bp_tgt)

print("== out_pump_pressure_comparison ==")
tgt = load("out_pump_pressure_comparison.csv")
mine = res["out_pump_pressure_comparison"]
for j, nm in enumerate(["time", "cumvol", "calc", "design_interp", "error"]):
    rep(f"col{j} {nm}", mine[:, j], tgt[:, j])

print("== out_bottom_ecd_comparison ==")
tgt = load("out_bottom_ecd_comparison.csv")
mine = res["out_bottom_ecd_comparison"]
for j, nm in enumerate(["time", "cumvol", "calc_bot_ECD", "design_interp", "error"]):
    rep(f"col{j} {nm}", mine[:, j], tgt[:, j])
