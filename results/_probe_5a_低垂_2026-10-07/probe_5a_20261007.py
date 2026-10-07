#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Phase 5a 探针：R5 低垂三项出数（①停泵冻结判据 / ④mud_extrapolate 对照 / ⑤缝隙律 A/B）

底座 = Phase 4d 主批 C 组：`anchored` 档 + `(F,F)` 角 + pressure off + nz250 + CFL 自适应
+ 2D 逐列（`enable_depthwise_temperature=True`）+ r1.0 —— 使 ①④⑤ 的数字与 4d 汇总表同底。

作业（5 run）：
- `呼101_Tanchored_col_5a`        基准（①判据 + ④⑤ 的对照底座）
- `呼101_Tanchored_col_mudext`    ④ `mud_extrapolate=True`（钻井液 T>80 °C 按公式外推）
- `呼101_Tanchored_col_gapoff`    ⑤ `enable_power_law_gap_law=False`（缝隙律关，旧路径对照）
- `呼1-003_Tanchored_col_5a`      ① 判据（第二口井）
- `呼1-004_Tanchored_col_5a`      ① 判据（第三口井）

说明：`smooth_break` **不跑**（`fluid_at` 对其显式 `raise NotImplementedError`，占位开关；
忠实两段式是既定口径「100 °C 断点不平滑」）——如实报告，不伪造对照。

运行：PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -I probe_5a_20261007.py --batch
"""

from __future__ import annotations

import argparse
import importlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
BRANCH_ROOT = HERE.parents[2]
for _p in (str(BRANCH_ROOT), str(BRANCH_ROOT / "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from entrypoints.run_sensitivity_current_20260916 import (  # noqa: E402
    annulus_kwargs_from_opts, build_pressure_field, build_temperature_fields,
    casing_kwargs_from_opts, extra_metrics, normalize_run_opts, scale_schedule,
)
from entrypoints.rerun_all_wells_corrected import (  # noqa: E402
    CORRECTED_KW, _stop_t, _total_t,
)
from cemdisp.diagnostics.stop_pump_freeze import (  # noqa: E402
    compute_stop_pump_freeze,
)
from cemdisp.models2d import AnnulusD2DGASolver  # noqa: E402
from cemdisp.models2d.boundary_bridge import (  # noqa: E402
    build_coupled_annulus_inlet_provider,
)
from cemdisp.transport1d import CasingFlowSolver  # noqa: E402

OUT = HERE.parent
WELL_LOADERS = {
    "呼101": ("cemdisp.data.loaders.hu101_loader", "load_hu101_tailpipe"),
    "呼1-003": ("cemdisp.data.loaders.ht1_003_loader", "load_ht1_003_tailpipe"),
    "呼1-004": ("cemdisp.data.loaders.ht1_004_loader", "load_ht1_004_tailpipe"),
}

JOBS = [
    {"tag": "呼101_Tanchored_col_5a", "well": "呼101"},
    {"tag": "呼101_Tanchored_col_mudext", "well": "呼101", "mud_extrapolate": True},
    {"tag": "呼101_Tanchored_col_gapoff", "well": "呼101", "gap_law": False},
    {"tag": "呼1-003_Tanchored_col_5a", "well": "呼1-003"},
    {"tag": "呼1-004_Tanchored_col_5a", "well": "呼1-004"},
]


def _beta_at(well, md: np.ndarray) -> np.ndarray:
    """井斜角剖面 beta(md)（度）——**复用求解器同款口径**（`_profile_to_arrays` +
    `np.interp`，见 `annulus_d2dga._build_geom` 的 inc_deg 行），保证与求解器内部一致。"""
    from cemdisp.models2d.annulus_d2dga import _profile_to_arrays
    prof = getattr(well, "inclination_profile", ())
    if not prof:
        return np.zeros_like(md)
    zs, vs = _profile_to_arrays(prof)
    return np.interp(md, zs, vs)


def run_job(job: dict) -> dict:
    tag, well_key = job["tag"], job["well"]
    mod, fn = WELL_LOADERS[well_key]
    loader = getattr(importlib.import_module(mod), fn)

    opts = normalize_run_opts({
        "enable_temperature_rheology": True,
        "temperature_mode": "anchored",
        "enable_yield_gate": None,
        "mud_extrapolate": bool(job.get("mud_extrapolate", False)),
        "pressure_mode": "off", "pressure_caliber": "shoe",
        "include_yield_term": False, "enable_stream_yield_gate": False,
    })

    t0 = time.perf_counter()
    well, fluids, schedule, _ = loader()
    schedule2 = scale_schedule(schedule, 1.0)
    field_1d, field_2d, note = build_temperature_fields(well_key, "anchored")
    p_field = build_pressure_field(well, fluids, schedule2, "off")

    casing = CasingFlowSolver(
        enable_gravity=True, mixing_contact_time=True, plug_face_zero_mixing=True,
        has_plug=True, pressure_field=p_field, **casing_kwargs_from_opts(opts))
    cr = casing.run(well, fluids, schedule2, temperature_field=field_1d)
    inlet = build_coupled_annulus_inlet_provider(
        cr, casing, fluids, split_cement_phases=True)
    stop_t = float(_stop_t(cr, fluids))
    tt = min(_total_t(schedule2) + 1200.0, stop_t)

    solver_kw = {**CORRECTED_KW, **annulus_kwargs_from_opts(opts),
                 "pressure_field": p_field}
    if job.get("gap_law", True) is False:
        solver_kw["enable_power_law_gap_law"] = False
    solver = AnnulusD2DGASolver(total_t=tt, nz=250, enable_cfl_adaptive=True,
                                enable_depthwise_temperature=True, **solver_kw)

    snap: dict = {}
    orig_refresh = solver._refresh_geom_temperature

    def refresh(geom, t_s):
        snap["geom"] = {k: np.array(geom[k], copy=True)
                        for k in ("b", "standoff", "md") if k in geom}
        return orig_refresh(geom, t_s)

    solver._refresh_geom_temperature = refresh
    orig_cv = solver._compute_velocity

    def cv(*a, **kw):
        out = orig_cv(*a, **kw)
        try:
            snap["rho"] = np.array(out[3], dtype=float, copy=True)
            snap["tau_y"] = np.array(out[8], dtype=float, copy=True)
        except Exception:
            pass
        return out

    solver._compute_velocity = cv
    res = solver.run(well, fluids, inlet, schedule=schedule2,
                     temperature_field=field_2d)
    elapsed = round(time.perf_counter() - t0, 1)

    extra = extra_metrics(res, res.summary, field_1d, field_2d)
    final = res.summary["最终结果"]
    eta_e = float(final["全井段最终有效顶替效率"])
    eta_n = float(final["窄四分位效率"])

    row = {
        "tag": tag, "井名": well_key, "档": "anchored", "note": note,
        "η_E": eta_e, "η_N": eta_n, "elapsed_s": elapsed,
        "配置": {"mud_extrapolate": opts["mud_extrapolate"],
                 "gap_law": job.get("gap_law", True),
                 "depthwise": True, "nz": 250, "stop_t_s": stop_t},
        "extra": extra,
        "health": {"步数": int(np.size(res.dt_history)),
                   "cfl_clip_events": int(res.cfl_clip_events),
                   "col_memo_size": len(getattr(solver, "_col_memo", {}))},
    }

    if "rho" in snap and "tau_y" in snap and "geom" in snap:
        g = snap["geom"]
        b2d = np.asarray(g["b"], dtype=float)
        b1d = b2d.mean(axis=0) if b2d.ndim > 1 else b2d
        md = np.asarray(g["md"], dtype=float)
        so = np.asarray(g["standoff"], dtype=float).ravel()
        fr = compute_stop_pump_freeze(
            tau_y_field=snap["tau_y"], rho_field=snap["rho"],
            gap_m=b1d, standoff=so, beta_deg=_beta_at(well, md))
        row["冻结判据"] = {
            "verdict": fr.verdict,
            "冻结占比": fr.frozen_frac,
            "tau_y_min_pa": {"min": float(fr.tau_y_min_pa.min()),
                             "median": float(np.median(fr.tau_y_min_pa)),
                             "max": float(fr.tau_y_min_pa.max())},
            "eccentricity": {"min": float(fr.eccentricity.min()),
                             "max": float(fr.eccentricity.max())},
            "delta_rho_az_kg_m3_median": float(np.median(fr.delta_rho_az_kg_m3)),
            "driving_pa_median": float(np.median(fr.driving_pa)),
            "driving_pa_max": float(fr.driving_pa.max()),
            "resisting_pa_median": float(np.median(fr.resisting_pa)),
            "ratio_p05": float(np.percentile(np.clip(fr.ratio, 0, 1e12), 5)),
            "ratio_p50": float(np.percentile(np.clip(fr.ratio, 0, 1e12), 50)),
        }
    else:
        row["冻结判据"] = None

    (OUT / "jobs").mkdir(parents=True, exist_ok=True)
    (OUT / "jobs" / f"{tag}.json").write_text(
        json.dumps(row, ensure_ascii=False, indent=1, default=str), encoding="utf-8")
    verdict = row["冻结判据"]["verdict"] if row["冻结判据"] else "N/A"
    print(f"[{tag}] eta_E={eta_e:.10f} eta_N={eta_n:.10f} freeze={verdict} "
          f"elapsed={elapsed}s", flush=True)
    return row


def run_batch(jobs: list[dict], workers: int) -> int:
    import os
    (OUT / "jobs").mkdir(parents=True, exist_ok=True)
    pending = list(jobs); running = []; failed = []; done = []
    env = dict(os.environ); env["PYTHONIOENCODING"] = "utf-8"; env["PYTHONUTF8"] = "1"
    py = sys.executable
    while pending or running:
        while pending and len(running) < workers:
            job = pending.pop(0)
            jf = OUT / "jobs" / f"_spec_{job['tag']}.json"
            jf.write_text(json.dumps(job, ensure_ascii=False), encoding="utf-8")
            logf = open(OUT / "jobs" / f"{job['tag']}.log", "w", encoding="utf-8")
            p = subprocess.Popen([py, str(HERE), "--job", str(jf)],
                                 stdout=logf, stderr=subprocess.STDOUT, env=env,
                                 cwd=str(BRANCH_ROOT))
            running.append((p, job["tag"], logf))
            print(f"[launcher] 启动 {job['tag']} (pid={p.pid})", flush=True)
        time.sleep(2.0)
        still = []
        for p, tag, logf in running:
            rc = p.poll()
            if rc is None:
                still.append((p, tag, logf)); continue
            logf.close()
            if rc == 0:
                done.append(tag)
            else:
                failed.append((tag, rc))
            print(f"[launcher] {'完成' if rc == 0 else 'FAIL'} {tag} exit={rc}", flush=True)
        running = still
    print(f"[launcher] 批完成：{len(done)} 成功 / {len(failed)} 失败 {failed}", flush=True)
    return 0 if not failed else 1


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--job", type=str, default=None)
    ap.add_argument("--batch", action="store_true")
    ap.add_argument("--workers", type=int, default=3)
    args = ap.parse_args()
    if args.job:
        run_job(json.loads(Path(args.job).read_text(encoding="utf-8")))
    elif args.batch:
        sys.exit(run_batch(JOBS, args.workers))
    else:
        ap.error("需 --job 或 --batch")


if __name__ == "__main__":
    main()
