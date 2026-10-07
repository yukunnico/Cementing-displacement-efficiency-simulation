#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Phase 4d 主批驱动：七井温度场（T-off / static / anchored）+ 呼1-004 扩展表档

范式 = `results/_probe_屈服门四角_2026-10-07/probe_four_corner_20261007.py`
（计划 §3-6：装配 = `run_variant_res` 逐行复制 + 纯透传观测包装 + 锚读文件禁手抄 +
hard 锚失败 exit 2 停批）。与本驱动的差异仅三处：

1. 七井井位（WELL_LOADERS）+ 逐作业 `mode`（off / static / anchored / table_ext）；
2. 2D `enable_depthwise_temperature` 由作业字段直装 ctor（计划 §3-2：run_opts 接线留 5d）；
3. `table_ext` 档在 `stop_t` 求得后、环空求解前调 `assert_time_table_coverage`
   （批次级前置门，计划 §8-14 禁静默 clamp），并附加 hard 判据 `oob_count == 0`。

口径（与「后4a 口径排量锚」同底座，`results/_probe_排量锚后4a_2026-10-07/`）：
`CORRECTED_KW` + (F,F) 角（`include_yield_term=False`、`enable_stream_yield_gate=False`）
+ `pressure_mode="off"` + `pressure_caliber="shoe"` + `nz=250` + CFL 自适应。

**hard 锚**（组 D）：呼101 × static × 逐列关 × r1.0 × (F,F) × pressure off
⇒ 逐位复现后4a 基线（读
`results/_probe_排量锚后4a_2026-10-07/jobs/hu101_r1.0_off_FF_post4a.json`，禁手抄）。

用法：
  单作业：python probe_seven_wells_20261007.py --job <job.json>
  批    ：python probe_seven_wells_20261007.py --batch <jobs.json> [--workers 4]
"""

from __future__ import annotations

import argparse
import importlib
import importlib.util
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
BRANCH_ROOT = HERE.parents[2]                       # repo 根（cement model_温压耦合分支）
for _p in (str(BRANCH_ROOT), str(BRANCH_ROOT / "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import cemdisp  # noqa: E402
from entrypoints.run_sensitivity_current_20260916 import (  # noqa: E402
    annulus_kwargs_from_opts,
    build_pressure_field,
    build_temperature_fields,
    casing_kwargs_from_opts,
    extra_metrics,
    normalize_run_opts,
    scale_schedule,
)
from entrypoints.rerun_all_wells_corrected import (  # noqa: E402
    CORRECTED_KW,
    _stop_t,
    _total_t,
)
from cemdisp.data.temperature_field import assert_time_table_coverage  # noqa: E402
from cemdisp.models2d import AnnulusD2DGASolver  # noqa: E402
from cemdisp.models2d import buoyancy  # noqa: E402
from cemdisp.models2d.boundary_bridge import (  # noqa: E402
    build_coupled_annulus_inlet_provider,
)
from cemdisp.transport1d import CasingFlowSolver  # noqa: E402

# ---- 复用四角探针的纯观测辅助件（同族产物，避免重复实现） -----------------
_FOUR_CORNER = (BRANCH_ROOT / "results" / "_probe_屈服门四角_2026-10-07"
                / "probe_four_corner_20261007.py")
_spec = importlib.util.spec_from_file_location("_probe_four_corner", _FOUR_CORNER)
_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_base)
_MuRecorder = _base._MuRecorder
_install_solver_recorders = _base._install_solver_recorders
_conservation_residual = _base._conservation_residual
_provenance = _base._provenance

OUT = HERE.parent                                  # results/_probe_4d七井_2026-10-07
ANCHOR_JSON = (BRANCH_ROOT / "results" / "_probe_排量锚后4a_2026-10-07" / "jobs"
               / "hu101_r1.0_off_FF_post4a.json")
ANCHOR_TOL = 1e-12      # 「逐位」判据：JSON 存全精度 repr，round-trip 精确

# 七井（spec §1：八井 − 呼103）
WELL_LOADERS = {
    "呼101": ("cemdisp.data.loaders.hu101_loader", "load_hu101_tailpipe"),
    "呼1-003": ("cemdisp.data.loaders.ht1_003_loader", "load_ht1_003_tailpipe"),
    "呼1-004": ("cemdisp.data.loaders.ht1_004_loader", "load_ht1_004_tailpipe"),
    "呼102": ("cemdisp.data.loaders.hu102_loader", "load_hu102_tailpipe"),
    "呼探1": ("cemdisp.data.loaders.hu1_loader", "load_hu1_tailpipe"),
    "呼探1-002": ("cemdisp.data.loaders.hu2_loader", "load_hu2_tailpipe"),
    "呼探1-001": ("cemdisp.data.loaders.ht1_001_loader", "load_ht1_001_tailpipe"),
}
_identity = lambda x: x                          # noqa: E731


def load_hard_anchor(job: dict) -> dict | None:
    """组 D hard 锚：读后4a 基线 JSON（禁手抄）。job["anchor"]=None ⇒ None。"""
    if not job.get("anchor"):
        return None
    if not ANCHOR_JSON.exists():
        raise FileNotFoundError(f"hard 锚文件缺失：{ANCHOR_JSON}")
    d = json.loads(ANCHOR_JSON.read_text(encoding="utf-8"))
    return {
        "η_E": float(d["η_E"]),
        "η_N": float(d["η_N"]),
        "source": ANCHOR_JSON.name,
        "mode": job["anchor"].get("mode", "hard"),
    }


def run_job(job: dict) -> dict:
    tag = job["tag"]
    well_key = job["well"]
    split = bool(job.get("split", False))
    gate = job.get("gate", None)          # None ⇒ 不出键（构造默认 False）
    nz = int(job.get("nz", 250))
    rate = float(job.get("rate", 1.0))
    mode = str(job.get("mode", "static"))
    pmode = str(job.get("pressure_mode", "off"))
    depthwise = bool(job.get("depthwise", False))

    mod_name, fn_name = WELL_LOADERS[well_key]
    loader = getattr(importlib.import_module(mod_name), fn_name)

    opts = normalize_run_opts({
        "enable_temperature_rheology": (mode != "off"),
        "temperature_mode": mode,
        "enable_yield_gate": None,        # 沿用 CORRECTED_KW（True）
        "pressure_mode": pmode,
        "pressure_caliber": "shoe",
        "include_yield_term": split,
        "enable_stream_yield_gate": gate,
    })

    t0 = time.perf_counter()
    well, fluids, schedule, _ = loader()
    well2, fluids2 = _identity(well), _identity(fluids)
    schedule2 = scale_schedule(schedule, rate)

    field_1d, field_2d, note = build_temperature_fields(well_key, mode)
    p_field = build_pressure_field(well2, fluids2, schedule2, opts["pressure_mode"])

    casing = CasingFlowSolver(
        enable_gravity=True, mixing_contact_time=True,
        plug_face_zero_mixing=True, has_plug=True,
        pressure_field=p_field,
        **casing_kwargs_from_opts(opts),
    )
    cr = casing.run(well2, fluids2, schedule2, temperature_field=field_1d)
    inlet = build_coupled_annulus_inlet_provider(
        cr, casing, fluids2, split_cement_phases=True,
    )
    stop_t = float(_stop_t(cr, fluids2))
    tt = min(_total_t(schedule2) + 1200.0, stop_t)

    # ---- table_ext 批次级前置门（计划 §8-14：禁静默 clamp）----
    coverage = None
    if mode == "table_ext":
        assert_time_table_coverage([field_1d, field_2d], stop_t, label=tag)
        coverage = {"stop_t_s": stop_t, "table_end_s": float(field_2d.time_s[-1]),
                    "assert": "PASS"}

    solver_kw = {**CORRECTED_KW, **annulus_kwargs_from_opts(opts),
                 "pressure_field": p_field}
    solver = AnnulusD2DGASolver(
        total_t=tt, nz=nz, enable_cfl_adaptive=True,
        enable_depthwise_temperature=depthwise, **solver_kw,
    )

    # 观测包装（纯透传）
    mu_rec = _MuRecorder()
    orig_fav = buoyancy.fluid_apparent_viscosity
    buoyancy.fluid_apparent_viscosity = mu_rec(orig_fav)
    m_rec, wall_series = _install_solver_recorders(solver)
    try:
        res = solver.run(well2, fluids2, inlet, schedule=schedule2,
                         temperature_field=field_2d)
    finally:
        buoyancy.fluid_apparent_viscosity = orig_fav
    elapsed = round(time.perf_counter() - t0, 1)

    extra = extra_metrics(res, res.summary, field_1d, field_2d)
    final = res.summary["最终结果"]
    eta_e = float(final["全井段最终有效顶替效率"])
    eta_n = float(final["窄四分位效率"])

    dt = np.asarray(res.dt_history, dtype=float)
    wall_arr = np.asarray(wall_series, dtype=float) if wall_series else np.empty((0, 2))
    health = {
        "步数": int(dt.size),
        "dt_median_s": float(np.median(dt)) if dt.size else None,
        "cfl_clip_events": int(res.cfl_clip_events),
        "cfl_clip_steps": int(res.cfl_clip_steps),
        "col_memo_size": len(getattr(solver, "_col_memo", {})),
        "col_batches": int(getattr(solver, "_col_batches", 0)),
        "uniform_field_hits": int(getattr(solver, "_uniform_field_hits", 0)),
        "col_audit_counts": dict(getattr(solver, "_col_audit_counts", {})),
        "oob_1d": int(getattr(field_1d, "oob_count", 0)),
        "oob_2d": int(getattr(field_2d, "oob_count", 0)),
        "oob_col_1d": int(getattr(field_1d, "oob_column_clamped_total", 0)),
        "oob_col_2d": int(getattr(field_2d, "oob_column_clamped_total", 0)),
        "守恒残差": _conservation_residual(res),
    }

    row = {
        "tag": tag, "井名": well_key, "变体": f"T{mode}_rate_x{float(rate)}",
        "配置": {
            "mode": mode, "note": note, "depthwise": depthwise,
            "split(include_yield_term)": split,
            "gate(enable_stream_yield_gate)": gate,
            "nz": nz, "ny": int(solver.ny), "rate": rate,
            "pressure_mode": pmode, "pressure_caliber": "shoe",
            "CORRECTED_KW": CORRECTED_KW, "tt_s": tt, "stop_t_s": stop_t,
            "opts": dict(opts),
        },
        "η_E": eta_e, "η_N": eta_n,
        "elapsed_s": elapsed,
        "extra": extra,
        "health": health,
        "coverage": coverage,
        "summary_浮力数_b": final.get("浮力数_b"),
        "provenance": _provenance(job),
    }

    npz_path = OUT / "jobs" / f"{tag}.npz"
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(npz_path, wall_series=wall_arr, dt_history=dt)
    row["npz"] = str(npz_path.relative_to(BRANCH_ROOT))

    # ---- 断言（表档真覆盖 / hard 锚逐位）----
    problems = []
    if mode in ("table", "table_ext") and (health["oob_1d"] or health["oob_2d"]):
        problems.append(
            f"表档 oob≠0（1D={health['oob_1d']} 2D={health['oob_2d']}）⇒ 覆盖不足"
        )
    anchor = load_hard_anchor(job)
    if anchor is not None:
        d_e = abs(eta_e - anchor["η_E"])
        d_n = abs(eta_n - anchor["η_N"])
        row["anchor"] = {**anchor, "Δη_E": d_e, "Δη_N": d_n,
                         "PASS": bool(d_e <= ANCHOR_TOL and d_n <= ANCHOR_TOL)}
        print(f"[{tag}] 锚({anchor['source']}) Δη_E={d_e:.3e} Δη_N={d_n:.3e} "
              f"{'PASS' if row['anchor']['PASS'] else 'FAIL'}", flush=True)
        if not row["anchor"]["PASS"] and anchor["mode"] == "hard":
            problems.append(
                f"hard 锚未复现：实跑 {eta_e!r}/{eta_n!r} vs 锚 "
                f"{anchor['η_E']!r}/{anchor['η_N']!r}"
            )
    row["problems"] = problems

    (OUT / "jobs" / f"{tag}.json").write_text(
        json.dumps(row, ensure_ascii=False, indent=1, default=str), encoding="utf-8")
    print(f"[{tag}] η_E={eta_e:.10f} η_N={eta_n:.10f} "
          f"饥饿={extra.get('饥饿份额')} wall占比={extra.get('屈服门_wall占比')} "
          f"步数={health['步数']} memo={health['col_memo_size']} "
          f"oob={health['oob_1d']}/{health['oob_2d']} clip={health['cfl_clip_events']} "
          f"耗时={elapsed}s", flush=True)
    if problems:
        print(f"致命：{tag} {problems}，停止。", flush=True)
        sys.exit(2)
    return row


# ---------------------------------------------------------------- 批模式
def run_batch(jobs: list[dict], workers: int) -> int:
    """子进程池：每作业独立进程（Windows spawn 安全）；hard 锚失败 ⇒ 停批。"""
    import os
    jobs_dir = OUT / "jobs"
    jobs_dir.mkdir(parents=True, exist_ok=True)
    pending = list(jobs)
    running: list[tuple] = []
    failed, done = [], []
    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "utf-8"
    env["PYTHONUTF8"] = "1"
    py = sys.executable
    while pending or running:
        while pending and len(running) < workers:
            job = pending.pop(0)
            jf = jobs_dir / f"_spec_{job['tag']}.json"
            jf.write_text(json.dumps(job, ensure_ascii=False), encoding="utf-8")
            logf = open(jobs_dir / f"{job['tag']}.log", "w", encoding="utf-8")
            p = subprocess.Popen([py, str(HERE), "--job", str(jf)],
                                 stdout=logf, stderr=subprocess.STDOUT, env=env,
                                 cwd=str(BRANCH_ROOT))
            running.append((p, job["tag"], logf))
            print(f"[launcher] 启动 {job['tag']} (pid={p.pid})，在跑 "
                  f"{len(running)}/{workers}", flush=True)
        time.sleep(2.0)
        still = []
        for p, tag, logf in running:
            rc = p.poll()
            if rc is None:
                still.append((p, tag, logf))
                continue
            logf.close()
            if rc == 0:
                done.append(tag)
                print(f"[launcher] 完成 {tag}", flush=True)
            else:
                failed.append((tag, rc))
                print(f"[launcher] ✗ {tag} exit={rc}"
                      f"{'（hard 锚/覆盖断言失败 ⇒ 停批）' if rc == 2 else ''}",
                      flush=True)
        running = still
        if any(rc == 2 for _, rc in failed):
            for p, tag, logf in running:
                p.terminate()
                logf.close()
            print("[launcher] 停批：hard 锚失败。已完成:", done, "失败:", failed, flush=True)
            return 2
    print(f"[launcher] 批完成。{len(done)} 成功，{len(failed)} 失败: {failed}", flush=True)
    return 0 if not failed else 1


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--job", type=str, default=None)
    ap.add_argument("--batch", type=str, default=None)
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    if args.job:
        run_job(json.loads(Path(args.job).read_text(encoding="utf-8")))
    elif args.batch:
        jobs = json.loads(Path(args.batch).read_text(encoding="utf-8"))
        sys.exit(run_batch(jobs, args.workers))
    else:
        ap.error("需给 --job 或 --batch")


if __name__ == "__main__":
    main()
