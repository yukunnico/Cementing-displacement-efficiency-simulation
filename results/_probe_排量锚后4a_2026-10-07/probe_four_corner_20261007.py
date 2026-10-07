"""Phase 2e 四角矩阵 / nz 收敛 / 排量锚 / f_safety 跑批驱动（2026-10-07）。

登记：Δ-P2-6（四角驱动脚本 = 新文件，改动面外，已登记且用户知悉）。

口径（续作计划 §4.0 统一口径）：
  T-on · temperature_mode="static"（地温式 T(z)=16.006+1.7598e-2·z）· rate 缩放 ·
  pressure_mode ∈ {hydrostatic, off}（默认 hydrostatic）· pressure_caliber="shoe" ·
  CORRECTED_KW · nz 可变（统一口径 250；nz 收敛检查用 100/140）· CFL 自适应 ·
  tt = min(泵总+1200, stop_t)。

装配 = `run_variant_res`（scripts/entrypoints/run_sensitivity_current_20260916.py）
的逐行复制，唯一差异 = ① nz 参数化（run_variant_res 钉死 NZ=250）② extra_kw
（yield_gate_f_safety 档）。含 C-08 压力场另路注入（casing 与 annulus 都吃
pressure_field）——与 run_variant_res 相同，勿退回探针旧形态（旧探针无压力场）。

观测包装（**纯透传记录，不改任何数值路径**；包装器只读实参/返回值并聚合统计）：
  ① `buoyancy.fluid_apparent_viscosity` 模块级包装 → A2.1 逐站点 μ 统计
     （站点标签 = 直接调用方函数名：_froude_squared_at / _buoyancy_number_at /
     _velocity_stream_function；η₁/η₂ 同站点按流体名区分）；
  ② `solver._compute_velocity` 实例包装 → m_field 值域（A2.1「m 范围」）；
  ③ `solver._yield_gate_wall` 实例包装 → wall 占比时程（A2.2；每泵注步一点，
     b 用实参 effective_b 加权）。
A2.2 其余量：dt_history 中位数/顶格占比（顶格 = dt ≥ solver.dt−1e-12，CFL 上限
即 ctor 基准 dt=4.0，见 annulus_d2dga `dt_step = min(dt_cfl, self.dt, ...)`）、
cfl_clip_events/cfl_clip_steps（res 字段）、守恒残差。
**守恒残差驱动级口径（登记于台账）**：末步各相场（lead/tail/spacer/flusher）
0≤φ≤1 越界幅值 + 总相和超 1 幅值 + 水泥/隔离液/flusher 快照全程最大越界——
仓内无现成"守恒残差"实现（tier0 无质量守恒项），此口径为纯诊断量，不进 summary。
A2.3：mu_reg_field / λ_op / γ̇_rep / wall_field 统计 + NPZ 落盘
（(F,T)→(T,T) wall 分布漂移在收集器算）。

锚（**读文件，禁手抄**；hard 锚失败 ⇒ exit 2 停批）：
  - t2csv     : results/敏感性变体_温压T2_2026-10-01/汇总表_温压T2.csv（井×变体行）
  - phase0    : results/敏感性补跑_三重点井_20261006/<井>_Ton_static_rate_x<r>_结果摘要.json
  - probe_gate: results/_probe_屈服门_2026-10-06/探针结果.json（B_gate_on 行）
  - probe_rate: results/_probe_排量敏感性_2026-10-02/探针结果.json（井×变体行）

用法：
  单作业：python probe_four_corner_20261007.py --job <job.json>
  批    ：python probe_four_corner_20261007.py --batch <jobs.json> [--workers 3]
          （子进程池；任一 hard 锚失败 exit 2 ⇒ 停批，不发后续作业）
"""
from __future__ import annotations

import argparse
import csv
import importlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

BRANCH_ROOT = Path(__file__).resolve().parents[2]
for _p in (str(BRANCH_ROOT), str(BRANCH_ROOT / "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# 自检：cemdisp 必须解析到分支树（铁律 §0.1-2）
import cemdisp  # noqa: E402
assert Path(cemdisp.__file__).resolve().is_relative_to(BRANCH_ROOT), \
    f"cemdisp 指向错误树: {cemdisp.__file__}"

from entrypoints.run_sensitivity_current_20260916 import (  # noqa: E402
    _identity, annulus_kwargs_from_opts, build_pressure_field,
    build_temperature_fields, casing_kwargs_from_opts, extra_metrics,
    normalize_run_opts, scale_schedule,
)
from entrypoints.rerun_all_wells_corrected import CORRECTED_KW, _stop_t, _total_t  # noqa: E402
from cemdisp.models2d import AnnulusD2DGASolver  # noqa: E402
from cemdisp.models2d import buoyancy  # noqa: E402
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider  # noqa: E402
from cemdisp.transport1d import CasingFlowSolver  # noqa: E402

OUT = Path(__file__).resolve().parent
RESULTS = BRANCH_ROOT / "results"
T2_CSV = RESULTS / "敏感性变体_温压T2_2026-10-01" / "汇总表_温压T2.csv"
PHASE0_DIR = RESULTS / "敏感性补跑_三重点井_20261006"
PROBE_GATE_JSON = RESULTS / "_probe_屈服门_2026-10-06" / "探针结果.json"
PROBE_RATE_JSON = RESULTS / "_probe_排量敏感性_2026-10-02" / "探针结果.json"

MODE = "static"
WELL_LOADERS = {
    "呼101": ("cemdisp.data.loaders.hu101_loader", "load_hu101_tailpipe"),
    "呼1-003": ("cemdisp.data.loaders.ht1_003_loader", "load_ht1_003_tailpipe"),
    "呼1-004": ("cemdisp.data.loaders.ht1_004_loader", "load_ht1_004_tailpipe"),
}
ANCHOR_TOL = 1e-12   # "逐位"判据：CSV/JSON 均存 17 位有效数字，round-trip 精确


# ---------------------------------------------------------------- 锚读取
def _rate_tag(rate: float) -> str:
    # 锚文件命名 = str(float)（"x1.0"/"x0.6"/"x1.4"），:g 会丢 ".0" 导致失配
    return f"x{float(rate)}"


def load_anchor(spec: dict | None, well: str, rate: float) -> dict | None:
    """按 anchor.kind 读权威文件 → {η_E, η_N, source, mode}。spec=None → None。"""
    if spec is None:
        return None
    kind, mode = spec["kind"], spec.get("mode", "hard")
    variant = f"Ton_{MODE}_rate_{_rate_tag(rate)}"
    if kind == "t2csv":
        with T2_CSV.open(encoding="utf-8-sig", newline="") as fh:
            for r in csv.DictReader(fh):
                if r["井名"] == well and r["变体"] == variant:
                    return {"η_E": float(r["η_E"]), "η_N": float(r["η_N"]),
                            "source": f"t2csv:{variant}", "mode": mode}
        raise RuntimeError(f"T2 CSV 未找到 {well} × {variant}")
    if kind == "phase0":
        p = PHASE0_DIR / f"{well}_Ton_static_rate_{_rate_tag(rate)}_结果摘要.json"
        d = json.loads(p.read_text(encoding="utf-8"))["最终结果"]
        return {"η_E": float(d["全井段最终有效顶替效率"]),
                "η_N": float(d["窄四分位效率"]),
                "source": f"phase0:{p.name}", "mode": mode}
    if kind == "probe_gate":
        rows = json.loads(PROBE_GATE_JSON.read_text(encoding="utf-8"))
        for r in rows:
            if r["井名"] == well and r["run"] == "B_gate_on":
                return {"η_E": r["η_E"], "η_N": r["η_N"],
                        "source": "probe_gate:B_gate_on", "mode": mode}
        raise RuntimeError(f"探针 JSON 未找到 {well} B_gate_on")
    if kind == "probe_rate":
        rows = json.loads(PROBE_RATE_JSON.read_text(encoding="utf-8"))
        for r in rows:
            if r["井名"] == well and r["变体"] == variant:
                return {"η_E": r["η_E"], "η_N": r["η_N"],
                        "source": f"probe_rate:{variant}", "mode": mode}
        raise RuntimeError(f"排量探针 JSON 未找到 {well} × {variant}")
    raise ValueError(f"未知锚类型 {kind!r}")


# ---------------------------------------------------------------- 观测包装
class _MuRecorder:
    """fluid_apparent_viscosity 的聚合记录器（纯透传）。"""

    def __init__(self):
        self.by_site: dict[tuple, dict] = {}
        self.split_min = float("inf")    # include_yield_term=True 的全局值域
        self.split_max = float("-inf")
        self.split_calls = 0

    def __call__(self, fn):
        def wrapper(fluid, shear_rate, *, include_yield_term=False):
            mu = fn(fluid, shear_rate, include_yield_term=include_yield_term)
            try:
                site = sys._getframe(1).f_code.co_name  # wrapper 的直接调用方（站点函数）
                key = (site, str(fluid.name), bool(include_yield_term))
                rec = self.by_site.setdefault(
                    key, {"n": 0, "min": float("inf"), "max": float("-inf"),
                          "sum": 0.0, "gmin": float("inf"), "gmax": float("-inf")})
                rec["n"] += 1
                rec["min"] = min(rec["min"], mu)
                rec["max"] = max(rec["max"], mu)
                rec["sum"] += mu
                g = float(shear_rate)
                rec["gmin"] = min(rec["gmin"], g)
                rec["gmax"] = max(rec["gmax"], g)
                if include_yield_term:
                    self.split_calls += 1
                    self.split_min = min(self.split_min, mu)
                    self.split_max = max(self.split_max, mu)
            except Exception:  # 观测失败不得影响数值路径
                pass
            return mu
        return wrapper

    def report(self) -> dict:
        sites = {}
        for (site, fluid, flag), r in sorted(self.by_site.items()):
            sites[f"{site}|{fluid}|split={'T' if flag else 'F'}"] = {
                "n": r["n"], "mu_min": r["min"], "mu_max": r["max"],
                "mu_mean": r["sum"] / max(r["n"], 1),
                "gamma_min": r["gmin"], "gamma_max": r["gmax"]}
        return {"逐站点": sites,
                "split_μ_min": None if self.split_calls == 0 else self.split_min,
                "split_μ_max": None if self.split_calls == 0 else self.split_max,
                "split_调用数": self.split_calls}


def _install_solver_recorders(solver):
    """②m_field 值域 ③wall 占比时程——实例级纯透传包装。"""
    mrec = {"n": 0, "min": float("inf"), "max": float("-inf"), "sum": 0.0}
    wrec_t = []   # (frac_pos, b_weighted_mean) 每泵注步一点

    orig_cv = solver._compute_velocity

    def cv_wrapper(*a, **kw):
        out = orig_cv(*a, **kw)
        try:
            m_field = out[7]
            if m_field is not None and np.size(m_field):
                m = np.asarray(m_field, dtype=float)
                mrec["n"] += 1
                mrec["min"] = min(mrec["min"], float(np.min(m)))
                mrec["max"] = max(mrec["max"], float(np.max(m)))
                mrec["sum"] += float(np.mean(m))
        except Exception:
            pass
        return out

    solver._compute_velocity = cv_wrapper

    orig_gate = solver._yield_gate_wall

    def gate_wrapper(w, geom_b, mu, tau_y, cement_ever, cement_local, f_safety):
        wall = orig_gate(w, geom_b, mu, tau_y, cement_ever, cement_local, f_safety)
        try:
            wa = np.asarray(wall, dtype=float)
            b = np.asarray(geom_b, dtype=float)
            frac = float(np.mean(wa > 0.0))
            bw = float(np.sum(wa * b) / max(np.sum(b), 1e-12))
            wrec_t.append((frac, bw))
        except Exception:
            pass
        return wall

    solver._yield_gate_wall = gate_wrapper
    return mrec, wrec_t


def _conservation_residual(res) -> dict:
    """驱动级守恒/有界残差（口径见模块 docstring；纯诊断不进 summary）。"""
    def _viol(a):
        if a is None or not np.size(a):
            return 0.0, 0.0
        arr = np.asarray(a, dtype=float)
        return float(max(0.0, -arr.min())), float(max(0.0, arr.max() - 1.0))

    out = {}
    for name in ("lead_field", "tail_field", "spacer_field", "flusher_field"):
        lo, hi = _viol(getattr(res, name, None))
        out[f"末步_{name}"] = {"低于0": lo, "高于1": hi}
    if all(getattr(res, n, None) is not None for n in ("lead_field", "tail_field", "spacer_field")):
        flush = res.flusher_field if res.flusher_field is not None else 0.0
        total = res.lead_field + res.tail_field + res.spacer_field + flush
        out["末步_总相和"] = {"min": float(np.min(total)), "max": float(np.max(total)),
                              "超1幅值": float(max(0.0, np.max(total) - 1.0))}
    snap_max = 0.0
    for snaps in (res.cement_snapshots, res.spacer_snapshots, res.flusher_snapshots):
        for s in snaps or ():
            arr = np.asarray(s, dtype=float)
            if arr.size:
                snap_max = max(snap_max, float(max(0.0, -arr.min())),
                               float(max(0.0, arr.max() - 1.0)))
    out["快照_最大越界"] = snap_max
    return out


# ---------------------------------------------------------------- 单作业
def run_job(job: dict) -> dict:
    tag = job["tag"]
    well_key = job["well"]
    split = bool(job.get("split", False))
    gate = job.get("gate", None)          # None ⇒ 不出键（构造默认 False）
    nz = int(job.get("nz", 250))
    rate = float(job.get("rate", 1.0))
    f_safety = job.get("f_safety", None)
    pmode = str(job.get("pressure_mode", "hydrostatic"))

    mod_name, fn_name = WELL_LOADERS[well_key]
    loader = getattr(importlib.import_module(mod_name), fn_name)

    opts = normalize_run_opts({
        "enable_temperature_rheology": True,
        "temperature_mode": MODE,
        "enable_yield_gate": None,        # 沿用 CORRECTED_KW（True）
        "pressure_mode": pmode,
        "pressure_caliber": "shoe",
        "include_yield_term": split,
        "enable_stream_yield_gate": gate,
    })

    t0 = time.perf_counter()
    well, fluids, schedule, _ = loader()
    well2, fluids2, schedule2 = _identity(well), _identity(fluids), scale_schedule(schedule, rate)

    field_1d, field_2d, note = build_temperature_fields(well_key, MODE)
    # C-08：压力场是对象（需 well/fluids/schedule）⇒ 另路构造，casing 与 annulus 同吃
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
    tt = min(_total_t(schedule2) + 1200.0, _stop_t(cr, fluids2))

    extra_kw = {} if f_safety is None else {"yield_gate_f_safety": float(f_safety)}
    solver_kw = {**CORRECTED_KW, **annulus_kwargs_from_opts(opts),
                 "pressure_field": p_field, **extra_kw}
    solver = AnnulusD2DGASolver(
        total_t=tt, nz=nz, enable_cfl_adaptive=True, **solver_kw,
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

    # ---- A2.1：拆分值域断言（clip 结构常数 [1e-5, 3.0]）----
    a21 = mu_rec.report()
    a21["m_field"] = ({"n": m_rec["n"], "min": m_rec["min"], "max": m_rec["max"],
                       "mean": m_rec["sum"] / max(m_rec["n"], 1)}
                      if m_rec["n"] else None)
    if split and mu_rec.split_calls > 0:
        assert 1.0e-5 - 1e-15 <= mu_rec.split_min <= mu_rec.split_max <= 3.0 + 1e-15, \
            f"A2.1 断言失败：split μ 值域 [{mu_rec.split_min}, {mu_rec.split_max}] ⊄ [1e-5, 3.0]"
        a21["A2.1_断言"] = "PASS（split μ ∈ [1e-5, 3.0]）"
    elif split:
        a21["A2.1_断言"] = "WARN：split=True 但记录器零调用（站点未触达？须排查）"
    else:
        a21["A2.1_断言"] = "N/A（split=False）"

    # ---- A2.2：健康度 ----
    dt = np.asarray(res.dt_history, dtype=float)
    dt_cap = float(solver.dt)
    wall_arr = np.asarray(wall_series, dtype=float) if wall_series else np.empty((0, 2))
    a22 = {
        "步数": int(dt.size),
        "dt_median_s": float(np.median(dt)) if dt.size else None,
        "dt_min_s": float(dt.min()) if dt.size else None,
        "dt_max_s": float(dt.max()) if dt.size else None,
        "dt_p05_s": float(np.percentile(dt, 5)) if dt.size else None,
        "dt_p95_s": float(np.percentile(dt, 95)) if dt.size else None,
        "dt_顶格步占比": float(np.mean(dt >= dt_cap - 1e-12)) if dt.size else None,
        "dt_顶格上限_s": dt_cap,
        "cfl_clip_events": int(res.cfl_clip_events),
        "cfl_clip_steps": int(res.cfl_clip_steps),
        "wall占比时程_点数": int(wall_arr.shape[0]),
        "wall占比时程_frac_min": float(wall_arr[:, 0].min()) if wall_arr.size else None,
        "wall占比时程_frac_max": float(wall_arr[:, 0].max()) if wall_arr.size else None,
        "wall占比时程_bw首": float(wall_arr[0, 1]) if wall_arr.size else None,
        "wall占比时程_bw末": float(wall_arr[-1, 1]) if wall_arr.size else None,
        "守恒残差": _conservation_residual(res),
    }

    # ---- A2.3：共线量 ----
    mu_reg = res.mu_reg_field
    wall_f = res.wall_field
    a23 = {
        "lambda_op_last": res.lambda_op_last,
        "shear_rate_rep_last": res.shear_rate_rep_last,
        "mu_reg_field": ({"min": float(np.min(mu_reg)), "max": float(np.max(mu_reg)),
                          "mean": float(np.mean(mu_reg))}
                         if mu_reg is not None and np.size(mu_reg) else None),
        "wall_field": ({"min": float(np.min(wall_f)), "max": float(np.max(wall_f)),
                        "mean": float(np.mean(wall_f)),
                        "frac_pos": float(np.mean(np.asarray(wall_f) > 0.0))}
                       if wall_f is not None and np.size(wall_f) else None),
    }

    row = {
        "tag": tag, "井名": well_key, "变体": f"Ton_{MODE}_rate_{_rate_tag(rate)}",
        "配置": {"split(include_yield_term)": split, "gate(enable_stream_yield_gate)": gate,
                 "nz": nz, "ny": int(solver.ny), "rate": rate, "f_safety": f_safety,
                 "pressure_mode": pmode, "pressure_caliber": "shoe",
                 "CORRECTED_KW": CORRECTED_KW, "tt_s": tt,
                 "opts": {k: (str(v) if v is not None and not isinstance(v, (bool, int, float, str)) else v)
                          for k, v in opts.items()}},
        "η_E": eta_e, "η_N": eta_n,
        "elapsed_s": elapsed,
        "extra": extra,
        "A2.1": a21, "A2.2": a22, "A2.3": a23,
        "summary_浮力数_b": res.summary["最终结果"].get("浮力数_b"),
        "provenance": _provenance(job),
    }

    # ---- NPZ 落盘（场量 + 时程）----
    npz_path = OUT / "jobs" / f"{tag}.npz"
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        npz_path,
        wall_field=(np.asarray(wall_f) if wall_f is not None else np.empty(0)),
        mu_reg_field=(np.asarray(mu_reg) if mu_reg is not None else np.empty(0)),
        wall_series=wall_arr,
        dt_history=dt,
    )
    row["npz"] = str(npz_path.relative_to(BRANCH_ROOT))

    # ---- 锚比对（JSON 先落盘，hard 失败也保留结果取证）----
    anchor = load_anchor(job.get("anchor"), well_key, rate)
    if anchor is not None:
        d_e = abs(eta_e - anchor["η_E"])
        d_n = abs(eta_n - anchor["η_N"])
        row["anchor"] = {**anchor, "Δη_E": d_e, "Δη_N": d_n,
                         "PASS": bool(d_e <= ANCHOR_TOL and d_n <= ANCHOR_TOL)}
    json_path = OUT / "jobs" / f"{tag}.json"
    json_path.write_text(json.dumps(row, ensure_ascii=False, indent=1, default=str),
                         encoding="utf-8")
    if anchor is not None:
        print(f"[{tag}] 锚({anchor['source']},{anchor['mode']}) "
              f"Δη_E={row['anchor']['Δη_E']:.3e} Δη_N={row['anchor']['Δη_N']:.3e} "
              f"{'PASS' if row['anchor']['PASS'] else 'FAIL'}", flush=True)
        if not row["anchor"]["PASS"] and anchor["mode"] == "hard":
            print(f"致命：hard 锚未复现 A={eta_e!r}/{eta_n!r} "
                  f"锚={anchor['η_E']!r}/{anchor['η_N']!r}，停止。", flush=True)
            sys.exit(2)
    print(f"[{tag}] η_E={eta_e:.10f} η_N={eta_n:.10f} 饥饿={extra.get('饥饿份额'):.6f} "
          f"wall占比={extra.get('屈服门_wall占比')} 步数={a22['步数']} "
          f"dt中位={a22['dt_median_s']} clip={a22['cfl_clip_events']} 耗时={elapsed}s", flush=True)
    return row


def _provenance(job: dict) -> dict:
    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=BRANCH_ROOT,
                              capture_output=True, text=True, timeout=30)
        head = head.stdout.strip()
    except Exception:
        head = "unknown"
    return {"git_HEAD": head, "job": job,
            "time": time.strftime("%Y-%m-%d %H:%M:%S"),
            "python": sys.version.split()[0],
            "cemdisp": cemdisp.__file__}


# ---------------------------------------------------------------- 批模式
def run_batch(jobs: list[dict], workers: int) -> int:
    """子进程池：每作业独立进程（Windows spawn 安全）；hard 锚失败 ⇒ 停批。"""
    import os
    jobs_dir = OUT / "jobs"
    jobs_dir.mkdir(parents=True, exist_ok=True)
    pending = list(jobs)
    running: list[tuple[subprocess.Popen, str, object]] = []
    failed = []
    done = []
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
            p = subprocess.Popen([py, str(Path(__file__).resolve()), "--job", str(jf)],
                                 stdout=logf, stderr=subprocess.STDOUT, env=env,
                                 cwd=str(BRANCH_ROOT))
            running.append((p, job["tag"], logf))
            print(f"[launcher] 启动 {job['tag']} (pid={p.pid})，在跑 {len(running)}/{workers}",
                  flush=True)
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
                      f"{'（hard 锚失败 ⇒ 停批）' if rc == 2 else ''}", flush=True)
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
    ap.add_argument("--workers", type=int, default=3)
    args = ap.parse_args()
    if args.job:
        job = json.loads(Path(args.job).read_text(encoding="utf-8"))
        run_job(job)
    elif args.batch:
        jobs = json.loads(Path(args.batch).read_text(encoding="utf-8"))
        sys.exit(run_batch(jobs, args.workers))
    else:
        ap.error("需 --job 或 --batch")


if __name__ == "__main__":
    main()
