"""Phase 4b 逐列温度 A/B 探针（呼101，2026-10-07）。

口径（= 四角批 §6.1 统一口径）：T-on · static（Geothermal T(z)=16.006+1.7598e-2·z）·
rate×1.0 · pressure_mode="hydrostatic" · pressure_caliber="shoe" · CORRECTED_KW ·
nz=250 · CFL 自适应 · tt=min(泵总+1200, stop_t)。

作业：
  A = enable_depthwise_temperature=False（列化关）——**hard 锚**：逐位复现
      results/_probe_屈服门四角_2026-10-07/jobs/hu101_FF_nz250_hydro_r1.0.json 的
      η_E/η_N（读 JSON 比对，禁手抄；容差 1e-12；失败 ⇒ exit 2）。
  B = enable_depthwise_temperature=True（列化开）——新口径数字 + 观测落盘。

装配 = `run_variant_res`（scripts/entrypoints/run_sensitivity_current_20260916.py）
的逐行复制（同四角驱动范式：extra_kw 直装 ctor kwarg，本波不做 run_opts 接线）。

落盘（关2/关4/关5 + F-3 声明证据）：
  均匀场快捷路径命中计数、列派生批次数、_col_memo 尺寸、fluid_at 审计列批聚合、
  逐站点 μ 记录器（_froude_squared_at/_buoyancy_number_at 应保持标量代表口径
  ——μ 值域逐站点记录为 F²/b_num 标量取值不变证据）、λ_op/γ̇_rep、summary 浮力数_b、
  dt 中位数、η_E/η_N A/B 并列 + Δη（口径升级声明；不设通过线）。

用法：
  PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python probe_4b_columnwise_20261007.py
"""
from __future__ import annotations

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
ANCHOR_JSON = (BRANCH_ROOT / "results/_probe_屈服门四角_2026-10-07/jobs"
               / "hu101_FF_nz250_hydro_r1.0.json")
# HEAD 逐位基线（归因取证产物）：三文件改动 stash 后在**原始 HEAD 2025cbc** 上跑同
# 一 A 作业捕获（attrib_head_20261007.py → attrib_head_result.json，2026-10-07）。
# 背景：四角锚捕获于 e3b5dee；上游提交 847efa6（Phase 4a/4c 1D 逐深温度，耦合入口
# 改变）已使 hu101 (F,F) 在 T-on static 下位移 Δη_E≈−2.6e-3/Δη_N≈−6.8e-3——与 4b
# 无关。关2 红线在 4b 的实际判据 = A 逐位复现【当前 HEAD】；四角旧锚失配如实记录并
# 归因（见"锚_旧四角失配归因"字段），**未做任何调参拟合**。
HEAD_BASELINE_JSON = Path(__file__).resolve().parent / "attrib_head_result.json"
WELL_LOADER = ("cemdisp.data.loaders.hu101_loader", "load_hu101_tailpipe")
ANCHOR_TOL = 1e-12


class _MuRecorder:
    """fluid_apparent_viscosity 聚合记录器（纯透传；同四角驱动范式）。"""

    def __init__(self):
        self.by_site: dict = {}

    def __call__(self, fn):
        def wrapper(fluid, shear_rate, *, include_yield_term=False):
            mu = fn(fluid, shear_rate, include_yield_term=include_yield_term)
            try:
                site = sys._getframe(1).f_code.co_name
                key = (site, str(fluid.name))
                rec = self.by_site.setdefault(key, {"n": 0, "min": float("inf"),
                                                    "max": float("-inf")})
                rec["n"] += 1
                rec["min"] = min(rec["min"], mu)
                rec["max"] = max(rec["max"], mu)
            except Exception:  # 观测失败不得影响数值路径
                pass
            return mu
        return wrapper

    def report(self) -> dict:
        return {f"{site}|{fl}": r for (site, fl), r in sorted(self.by_site.items())}


def run_job(tag: str, depthwise: bool) -> dict:
    t0 = time.perf_counter()
    mod_name, fn_name = WELL_LOADER
    loader = getattr(importlib.import_module(mod_name), fn_name)

    opts = normalize_run_opts({
        "enable_temperature_rheology": True,
        "temperature_mode": "static",
        "enable_yield_gate": None,        # 沿用 CORRECTED_KW（True）
        "pressure_mode": "hydrostatic",
        "pressure_caliber": "shoe",
        "include_yield_term": False,      # 四角 (F,F) 同配置
        "enable_stream_yield_gate": None,
    })

    well, fluids, schedule, _ = loader()
    well2, fluids2, schedule2 = _identity(well), _identity(fluids), scale_schedule(schedule, 1.0)
    field_1d, field_2d, _note = build_temperature_fields("呼101", "static")
    p_field = build_pressure_field(well2, fluids2, schedule2, opts["pressure_mode"])

    casing = CasingFlowSolver(
        enable_gravity=True, mixing_contact_time=True,
        plug_face_zero_mixing=True, has_plug=True,
        pressure_field=p_field,
        **casing_kwargs_from_opts(opts),
    )
    cr = casing.run(well2, fluids2, schedule2, temperature_field=field_1d)
    inlet = build_coupled_annulus_inlet_provider(cr, casing, fluids2, split_cement_phases=True)
    tt = min(_total_t(schedule2) + 1200.0, _stop_t(cr, fluids2))

    # 探针 ctor kwarg 直装（本波不做 run_opts 接线——范式=四角驱动 extra_kw）。
    # A（关）⇒ **不出键**（沿用"默认不产生额外 kwarg"不变量 ⇒ 可在改动前 HEAD 上
    # 复跑同一作业做归因取证）。
    extra_kw = {"enable_depthwise_temperature": True} if depthwise else {}
    solver_kw = {**CORRECTED_KW, **annulus_kwargs_from_opts(opts),
                 "pressure_field": p_field, **extra_kw}
    solver = AnnulusD2DGASolver(total_t=tt, nz=250, enable_cfl_adaptive=True, **solver_kw)

    mu_rec = _MuRecorder()
    orig_fav = buoyancy.fluid_apparent_viscosity
    buoyancy.fluid_apparent_viscosity = mu_rec(orig_fav)
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
    row = {
        "tag": tag, "井名": "呼101", "变体": "Ton_static_rate_x1.0",
        "配置": {"depthwise_temperature": depthwise, "nz": 250, "ny": int(solver.ny),
                 "pressure_mode": "hydrostatic", "pressure_caliber": "shoe",
                 "include_yield_term": False, "enable_stream_yield_gate": False,
                 "CORRECTED_KW": CORRECTED_KW, "tt_s": tt},
        "η_E": eta_e, "η_N": eta_n, "elapsed_s": elapsed,
        "summary_浮力数_b": final.get("浮力数_b"),
        "extra_饥饿份额": extra.get("饥饿份额"),
        "extra_front_narrow_m": extra.get("front_narrow_m"),
        "extra_温度流变审计": extra.get("温度流变审计"),
        "观测": {
            "均匀场快捷路径命中": int(solver._uniform_field_hits),
            "列派生批次数": int(solver._col_batches),
            "列memo尺寸": int(len(solver._col_memo)),
            "列审计聚合": dict(solver._col_audit_counts),
            "lambda_op_last": res.lambda_op_last,
            "shear_rate_rep_last": res.shear_rate_rep_last,
            "dt中位_s": float(np.median(dt)) if dt.size else None,
            "步数": int(dt.size),
            "cfl_clip_events": int(res.cfl_clip_events),
            "逐站点μ(F²/b_num 标量口径证据)": mu_rec.report(),
        },
        "provenance": _provenance({"depthwise": depthwise}),
    }
    (OUT / f"{tag}.json").write_text(
        json.dumps(row, ensure_ascii=False, indent=1, default=str), encoding="utf-8")

    if not depthwise:  # A：hard 锚（读 JSON，禁手抄；失败 exit 2）
        # 主判据（关2 红线 = 列化关逐位复现当前 HEAD 基线）
        hb = json.loads(HEAD_BASELINE_JSON.read_text(encoding="utf-8"))
        h_e, h_n = float(hb["η_E"]), float(hb["η_N"])
        hd_e, hd_n = abs(eta_e - h_e), abs(eta_n - h_n)
        row["anchor_head_baseline"] = {
            "source": str(HEAD_BASELINE_JSON), "η_E": h_e, "η_N": h_n,
            "Δη_E": hd_e, "Δη_N": hd_n,
            "PASS": bool(hd_e <= ANCHOR_TOL and hd_n <= ANCHOR_TOL)}
        # 次判据（spec §2 关4 原文四角锚；捕获于 e3b5dee ⇒ 预期失配，失配归因见
        # attrib_head_20261007.py 取证：未改 4b 三文件的原始 HEAD 跑出同一失配值）
        anchor = json.loads(ANCHOR_JSON.read_text(encoding="utf-8"))
        a_e, a_n = float(anchor["η_E"]), float(anchor["η_N"])
        d_e, d_n = abs(eta_e - a_e), abs(eta_n - a_n)
        row["anchor_旧四角失配归因"] = {
            "source": str(ANCHOR_JSON), "η_E": a_e, "η_N": a_n,
            "Δη_E": d_e, "Δη_N": d_n,
            "PASS": bool(d_e <= ANCHOR_TOL and d_n <= ANCHOR_TOL),
            "说明": ("失配来自上游 847efa6（Phase 4a/4c 1D 逐深温度⇒耦合入口变化），"
                    "非 4b 列化：三文件 stash 后原始 HEAD 复跑逐位得到同一 A 值"
                    "（attrib_head_result.json）。重锚与否=用户裁定项。")}
        print(f"[{tag}] HEAD基线锚 Δη_E={hd_e:.3e} Δη_N={hd_n:.3e} "
              f"{'PASS' if row['anchor_head_baseline']['PASS'] else 'FAIL'}；"
              f"旧四角锚 Δη_E={d_e:.3e} Δη_N={d_n:.3e} "
              f"{'PASS' if row['anchor_旧四角失配归因']['PASS'] else 'FAIL(预期)'}",
              flush=True)
        (OUT / f"{tag}.json").write_text(
            json.dumps(row, ensure_ascii=False, indent=1, default=str), encoding="utf-8")
        if not row["anchor_head_baseline"]["PASS"]:
            print(f"致命：A 未逐位复现当前 HEAD 基线 A={eta_e!r}/{eta_n!r} "
                  f"基线={h_e!r}/{h_n!r} ⇒ 4b 改动的列化关不 inert，停止。", flush=True)
            sys.exit(2)
    print(f"[{tag}] depthwise={depthwise} η_E={eta_e:.10f} η_N={eta_n:.10f} "
          f"耗时={elapsed}s 列批={row['观测']['列派生批次数']} "
          f"memo={row['观测']['列memo尺寸']} 命中={row['观测']['均匀场快捷路径命中']}",
          flush=True)
    return row


def _provenance(job: dict) -> dict:
    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=BRANCH_ROOT,
                              capture_output=True, text=True, timeout=30).stdout.strip()
    except Exception:
        head = "unknown"
    return {"git_HEAD": head, "job": job, "time": time.strftime("%Y-%m-%d %H:%M:%S"),
            "python": sys.version.split()[0], "cemdisp": cemdisp.__file__}


def main() -> None:
    a = run_job("A_col_off", depthwise=False)
    b = run_job("B_col_on", depthwise=True)
    summary = {
        "口径声明": ("Phase 4b 列化 = 口径升级（新旧并列）。A=列化关（=Phase 2e 四角 (F,F) "
                    "同配置），B=列化开（静温梯度场逐列物性）。Δη 如实，不设通过线。"),
        "A": {k: a[k] for k in ("η_E", "η_N", "summary_浮力数_b", "elapsed_s")},
        "B": {k: b[k] for k in ("η_E", "η_N", "summary_浮力数_b", "elapsed_s")},
        "Δ(B−A)": {"η_E_pp": (b["η_E"] - a["η_E"]) * 100.0,
                   "η_N_pp": (b["η_N"] - a["η_N"]) * 100.0,
                   "浮力数_b": (b["summary_浮力数_b"] - a["summary_浮力数_b"])
                   if (a["summary_浮力数_b"] is not None and b["summary_浮力数_b"] is not None)
                   else None},
        "F-3 声明": ("λ_op/F²/b_num/γ̇_rep/诊断维持域均标量口径：B 的 λ_op_last/shear_rate_rep_last/"
                    "浮力数_b 均为标量 float；逐站点 μ 记录器显示 _froude_squared_at/"
                    "_buoyancy_number_at 消费标量代表流体（η 列只进两层闭包与场公式）。"
                    "A/B 的 μ 值域差异来自流场 w 变化对 γ̇_rep 的反馈（合法动力学效应），"
                    "非口径漂移；形态差异=η₁/η₂ 列场 × F² 标量（F-3 代价条款）。"),
        "锚_HEAD基线": a.get("anchor_head_baseline"),
        "锚_旧四角失配归因": a.get("anchor_旧四角失配归因"),
        "观测": {"A": a["观测"], "B": b["观测"]},
    }
    (OUT / "AB结果.json").write_text(json.dumps(summary, ensure_ascii=False, indent=1,
                                               default=str), encoding="utf-8")
    print(json.dumps(summary["Δ(B−A)"], ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
