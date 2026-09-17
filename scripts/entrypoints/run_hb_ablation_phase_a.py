# -*- coding: utf-8 -*-
"""Phase A（B&F25 HB 闭包）Task 7：三分量 A/B 消融矩阵 + 报告。

用户已裁定的矩阵（2026-09-17，controller 转达；不得静默更改/缩小）
----------------------------------------------------------------
**7 行 × 4 井**（井 = hu101 / hu102 / ht1_003 / ht1_004）：

======  ===========================================================
H0      全默认（B-2 off）——L1 锚（与生产 runner 同配置逐位）
H0p     ``enable_stream_yield_gate=True``（B-2 on 基线）
H1      H0p + ``enable_hb_closure=True``（只换闭包，τ_Y2=0、屈服门不动）
H2a_LS  H0p + ``hb_fix_cement_tau_y=True`` + 口径 = Bingham-LS 截距
H2b_HB  H0p + ``hb_fix_cement_tau_y=True`` + 口径 = HB 拟合 τy
H3a_LS  H0p + 两开关全开 + LS 口径
H3b_HB  H0p + 两开关全开 + HB 口径
======  ===========================================================

归因分解（spec §5.1）：``H1−H0p``=闭包效应、``H2−H0p``=屈服门通道效应、
``H3−H0p``=两者、``H3−H1−H2+H0p``=交互项、``H0p−H0``=B-2 算子耦合效应本身。
已知且如实呈现：HB 路径上 B-2 无效（非线性入口无 wall 形参，Task 6 一次性
告警）——H1/H3 行的 B-2 是惰性的。

ht1_004 特例（用户裁定 R-T6-6 = spec 优先）：该井水泥 spec 自带 YP（13/14 Pa）
⇒ 其 a/b 口径变体完全重合——H2/H3 各跑一次（spec 口径），CSV 里两个标签都指向
该结果，``annotation`` 列写明"spec 优先，a/b 重合"。

τy 取数纪律（R4/R7 红线）
------------------------
- 只读 :mod:`cemdisp.data.cement_yield_stress`（Task 0 模块，**绝不修改**）；
- Bingham-LS 口径 = 逐相读 ``CEMENT_YIELD_STRESS[well][phase].fitted_bingham_ls_intercept_pa``；
- HB 口径 = ``yield_stress_by_role(well_key, "fitted_new")``；
- 解析映射（井 → role → 口径 → 数值 + 出处字段）落 manifest 供审计；
- **禁止**在本脚本编造/改写任何 τy 值；缺失（``MISSING`` 哨兵）⇒ 按求解器既有
  语义告警跳过（贡献 0），manifest 记录之。

执行口径
--------
- 每井数据装载**复用生产 runner 的路径**：``load_<well>_tailpipe()`` →
  ``CasingFlowSolver(enable_gravity=True, mixing_contact_time=True,
  plug_face_zero_mixing=True, has_plug=True)`` → 鞋口时序 →
  ``annulus_stop_time_s``（F2 口径）→ ``build_coupled_annulus_inlet_provider(
  ..., split_cement_phases=True)`` → ``AnnulusD2DGASolver(total_t=stop_t, nz=250,
  **变体开关)``——与 ``cemdisp/runners/<well>_tailpipe.py`` 唯一差别 = 注入的
  变体开关，保证 H0 与生产 runner 同配置。
- fallback 纪律：HB 非线性外迭代回退（RuntimeWarning"回退牛顿线性闭包"）
  计数 > 0 的 run 在汇总 CSV 标 ``CONTAMINATED``（其 HB 结果被牛顿回退污染），
  不静默混入。
- 诊断：反求侧 spy（``solve_g_from_mean_velocity_batch`` 包装，只统计不改动）
  记录每时间步最后反求轮的 undefined 数（n_static 口径）与全场 undefined 步数；
  闭包侧 R-T1-6 地板告警（"格点闭包无定义"，每闭包实例至多一次）另行计数并
  解析格点数。
- 可断点续跑：启动读 manifest，``status`` 为 DONE/CONTAMINATED 的 run 跳过。
- 输出：``results/HB闭包A-B_2026-09-17/``（逐 run 子目录 + manifest.json +
  汇总.csv + 归因分解.csv + l1_determinism.json）。

用法::

    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/entrypoints/run_hb_ablation_phase_a.py
    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/entrypoints/run_hb_ablation_phase_a.py --list
    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/entrypoints/run_hb_ablation_phase_a.py --l3
    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/entrypoints/run_hb_ablation_phase_a.py --wells hu101 --variants H0,H0p
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
import warnings
from pathlib import Path
from typing import Any, Callable, Mapping

import numpy as np

from cemdisp.data import loaders
from cemdisp.data.cement_yield_stress import (
    CEMENT_YIELD_STRESS,
    yield_stress_by_role,
)
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.runners.hu101_tailpipe import annulus_stop_time_s
from cemdisp.runners.zhang2022_benchmark import (
    ZHANG2022_CASES,
    annulus_volume_of,
    run_case,
)
from cemdisp.transport1d import CasingFlowSolver

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = PROJECT_ROOT / "results" / "HB闭包A-B_2026-09-17"
MANIFEST_PATH = OUT_DIR / "manifest.json"
SUMMARY_CSV = OUT_DIR / "汇总.csv"
ATTRIBUTION_CSV = OUT_DIR / "归因分解.csv"
L3_DIR = OUT_DIR / "l3_zhang2022"

WELLS: dict[str, Callable[[], tuple]] = {
    "hu101": loaders.load_hu101_tailpipe,
    "hu102": loaders.load_hu102_tailpipe,
    "ht1_003": loaders.load_ht1_003_tailpipe,
    "ht1_004": loaders.load_ht1_004_tailpipe,
}

_ROLE_BY_PHASE = {"lead": "LEAD", "intermediate": "INTERMEDIATE", "tail": "TAIL"}
_CEMENT_ROLE_NAMES = ("LEAD", "INTERMEDIATE", "TAIL")


# ---------------------------------------------------------------------------
# τy 口径取数（只读 Task 0 模块；解析映射落 manifest 供审计，禁止编造）
# ---------------------------------------------------------------------------
def build_tau_y_calibers(well_key: str) -> dict[str, Any]:
    """井 → 两套口径的 role 映射 + 逐 role 出处明细（审计用）。

    返回 ``{"bingham_ls": {...}, "hb_fitted_new": {...}, "detail": {...}}``；
    值为 Task 0 冻结常数原样（``MISSING`` 哨兵原样透传，由求解器告警跳过）。
    """
    hb_map = dict(yield_stress_by_role(well_key, "fitted_new"))
    ls_map: dict[str, Any] = {}
    detail: dict[str, Any] = {
        "bingham_ls": {"values": {}, "provenance": {}},
        "hb_fitted_new": {"values": dict(hb_map), "provenance": {}},
    }
    for phase, role in _ROLE_BY_PHASE.items():
        record = CEMENT_YIELD_STRESS[well_key].get(phase)
        if record is None:
            continue
        ls_map[role] = record.fitted_bingham_ls_intercept_pa
        detail["bingham_ls"]["values"][role] = record.fitted_bingham_ls_intercept_pa
        detail["bingham_ls"]["provenance"][role] = {
            "phase": phase,
            "source_file": record.source_file,
            "source_location": record.source_location,
            "field": "fitted_bingham_ls_intercept_pa",
            "in_rheometer_csv": record.in_rheometer_csv,
        }
        detail["hb_fitted_new"]["provenance"][role] = {
            "phase": phase,
            "source_file": record.source_file,
            "source_location": record.source_location,
            "field": "fitted_hb_tau_y_pa（route=fitted_new）",
            "in_rheometer_csv": record.in_rheometer_csv,
        }
    return {"bingham_ls": ls_map, "hb_fitted_new": hb_map, "detail": detail}


# ---------------------------------------------------------------------------
# 变体矩阵（用户裁定 7 行；ht1_004 按 R-T6-6 收敛为 6 run/井 + H0_repeat）
# ---------------------------------------------------------------------------
def build_run_plan(well_key: str) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """返回该井的 run 计划（有序）与 τy 口径审计块。"""
    calibers = build_tau_y_calibers(well_key)
    spec_priority = well_key == "ht1_004"

    def _gate_only() -> dict[str, Any]:
        return {"enable_stream_yield_gate": True}

    def _closure() -> dict[str, Any]:
        return {"enable_stream_yield_gate": True, "enable_hb_closure": True}

    def _gate_ty(mapping: Mapping[str, Any]) -> dict[str, Any]:
        return {"enable_stream_yield_gate": True, "hb_fix_cement_tau_y": True,
                "cement_tau_y_by_role": dict(mapping)}

    def _both(mapping: Mapping[str, Any]) -> dict[str, Any]:
        return {"enable_stream_yield_gate": True, "enable_hb_closure": True,
                "hb_fix_cement_tau_y": True, "cement_tau_y_by_role": dict(mapping)}

    spec_note = "spec 优先，a/b 重合（R-T6-6）：水泥 spec 自带 YP，常数映射未消费"
    plan: list[dict[str, Any]] = [
        {"run_id": f"{well_key}__H0", "labels": ["H0"], "caliber": "",
         "kwargs": {}, "annotation": "全默认（L1 锚，B-2 off，与生产 runner 同配置）"},
        {"run_id": f"{well_key}__H0p", "labels": ["H0p"], "caliber": "",
         "kwargs": _gate_only(), "annotation": "B-2 on 基线"},
        {"run_id": f"{well_key}__H1", "labels": ["H1"], "caliber": "",
         "kwargs": _closure(),
         "annotation": "H0p+HB 闭包（τ_Y2=0、屈服门不动；B-2 在 HB 路径惰性，一次性告警留痕）"},
    ]
    if spec_priority:
        # 传入映射 = HB 口径（被 spec 覆盖、不消费；两套口径值都落 manifest 审计）。
        plan += [
            {"run_id": f"{well_key}__H2_spec", "labels": ["H2a_LS", "H2b_HB"],
             "caliber": "spec 优先", "kwargs": _gate_ty(calibers["hb_fitted_new"]),
             "annotation": spec_note},
            {"run_id": f"{well_key}__H3_spec", "labels": ["H3a_LS", "H3b_HB"],
             "caliber": "spec 优先", "kwargs": _both(calibers["hb_fitted_new"]),
             "annotation": spec_note + "；HB 路径上 B-2 惰性"},
        ]
    else:
        plan += [
            {"run_id": f"{well_key}__H2a_LS", "labels": ["H2a_LS"],
             "caliber": "bingham_ls", "kwargs": _gate_ty(calibers["bingham_ls"]),
             "annotation": "H0p+水泥 τy 注入（Bingham-LS 截距口径，屈服门通道）"},
            {"run_id": f"{well_key}__H2b_HB", "labels": ["H2b_HB"],
             "caliber": "hb_fitted_new", "kwargs": _gate_ty(calibers["hb_fitted_new"]),
             "annotation": "H0p+水泥 τy 注入（HB 拟合口径，屈服门通道）"},
            {"run_id": f"{well_key}__H3a_LS", "labels": ["H3a_LS"],
             "caliber": "bingham_ls", "kwargs": _both(calibers["bingham_ls"]),
             "annotation": "H0p+两开关全开（LS 口径）；HB 路径上 B-2 惰性"},
            {"run_id": f"{well_key}__H3b_HB", "labels": ["H3b_HB"],
             "caliber": "hb_fitted_new", "kwargs": _both(calibers["hb_fitted_new"]),
             "annotation": "H0p+两开关全开（HB 口径，生产口径候选）；HB 路径上 B-2 惰性"},
        ]
    plan.append({"run_id": f"{well_key}__H0_repeat", "labels": ["H0_repeat"],
                 "caliber": "", "kwargs": {},
                 "annotation": "L1 确定性检查：H0 同配置第二次运行（逐位对比）"})
    return plan, calibers


# ---------------------------------------------------------------------------
# 反求/外迭代 spy（只统计、不改行为）
# ---------------------------------------------------------------------------
class NonlinearSpy:
    """包装 ``annulus_d2dga.solve_stream_function_nonlinear`` 与
    ``stream_function.solve_g_from_mean_velocity_batch``：

    - 每个非线性外迭代调用 = 一个"步"bucket（含冷/热耗时）；
    - 每次批量反求记录 undefined 数（最后反求轮 = 该步 n_static 口径）；
    - 异常（含回退前的 RuntimeError）原样抛出，只留痕。
    """

    def __init__(self) -> None:
        self.steps: list[dict[str, Any]] = []
        self.total_inverse_calls = 0
        self._current: dict[str, Any] | None = None
        self._orig_nl: Any = None
        self._orig_inv: Any = None

    def __enter__(self) -> "NonlinearSpy":
        import cemdisp.models2d.annulus_d2dga as ad
        import cemdisp.models2d.stream_function as sf

        self._ad, self._sf = ad, sf
        self._orig_nl = ad.solve_stream_function_nonlinear
        self._orig_inv = sf.solve_g_from_mean_velocity_batch
        ad.solve_stream_function_nonlinear = self._wrap_nonlinear
        sf.solve_g_from_mean_velocity_batch = self._wrap_inverse
        return self

    def __exit__(self, *exc: object) -> bool:
        self._ad.solve_stream_function_nonlinear = self._orig_nl
        self._sf.solve_g_from_mean_velocity_batch = self._orig_inv
        return False

    def _wrap_nonlinear(self, *args: Any, **kwargs: Any) -> Any:
        bucket: dict[str, Any] = {"inverse_calls": 0, "undefined_final": None,
                                  "undefined_max": 0, "cells": None, "ok": False}
        self._current = bucket
        self.steps.append(bucket)
        t0 = time.perf_counter()
        try:
            out = self._orig_nl(*args, **kwargs)
        finally:
            bucket["dt_s"] = time.perf_counter() - t0
        bucket["ok"] = True
        return out

    def _wrap_inverse(self, *args: Any, **kwargs: Any) -> Any:
        result = self._orig_inv(*args, **kwargs)
        self.total_inverse_calls += 1
        bucket = self._current
        if bucket is not None:
            bucket["inverse_calls"] += 1
            undefined = int(np.asarray(result.undefined).sum())
            bucket["cells"] = int(np.asarray(result.undefined).size)
            bucket["undefined_final"] = undefined
            bucket["undefined_max"] = max(bucket["undefined_max"], undefined)
        return result

    # ---- 聚合 ------------------------------------------------------------
    def aggregate(self) -> dict[str, Any]:
        steps = self.steps
        ok_steps = [s for s in steps if s["ok"]]
        dts = [s["dt_s"] for s in ok_steps]
        finals = [s["undefined_final"] for s in steps if s["undefined_final"] is not None]
        return {
            "nonlinear_steps": len(steps),
            "fallback_steps": sum(1 for s in steps if not s["ok"]),
            "inverse_calls_total": self.total_inverse_calls,
            "cold_step_s": round(dts[0], 3) if dts else None,
            "warm_step_median_s": (round(float(np.median(dts[1:])), 3) if len(dts) > 1 else None),
            "warm_step_max_s": (round(float(np.max(dts[1:])), 3) if len(dts) > 1 else None),
            "max_n_static": max(finals) if finals else None,
            "all_undefined_steps": (
                sum(1 for s in steps
                    if s["cells"] and s["undefined_final"] is not None
                    and s["undefined_final"] == s["cells"])),
        }


def classify_warnings(caught: list[warnings.WarningMessage]) -> dict[str, Any]:
    """告警分类计数：回退 / HB-忽略wall / τy 跳过汇总 / 闭包地板 / 其它。"""
    out: dict[str, Any] = {"fallback_count": 0, "hb_wall_ignored_count": 0,
                           "tau_y_skip_flush_count": 0, "floor_warning_count": 0,
                           "floor_warning_cells": [], "others": []}
    for w in caught:
        message = str(w.message)
        if issubclass(w.category, RuntimeWarning) and "回退牛顿线性闭包" in message:
            out["fallback_count"] += 1
        elif "不消费屈服门冻结度" in message:
            out["hb_wall_ignored_count"] += 1
        elif "hb_fix_cement_tau_y=True：井" in message and "显式跳过项" in message:
            out["tau_y_skip_flush_count"] += 1
            out.setdefault("tau_y_skip_messages", []).append(message)
        elif "格点闭包无定义" in message:
            out["floor_warning_count"] += 1
            digits = "".join(ch if ch.isdigit() else " " for ch in message.split("个格点")[0])
            parts = digits.split()
            if parts:
                out["floor_warning_cells"].append(int(parts[-1]))
        else:
            out["others"].append(f"{w.category.__name__}: {message[:160]}")
    out["others"] = sorted(set(out["others"]))[:20]
    return out


# ---------------------------------------------------------------------------
# 求解执行（复用生产 runner 路径，只注入变体开关）
# ---------------------------------------------------------------------------
def _digest(arr: Any) -> str:
    return hashlib.sha256(
        np.ascontiguousarray(np.asarray(arr, dtype=float)).tobytes()).hexdigest()


def execute_run(run: Mapping[str, Any], well_key: str,
                loader: Callable[[], tuple]) -> dict[str, Any]:
    """跑一个 run（1D 装载 → 2D 求解 → 指标/诊断收集），返回 manifest 记录。"""
    run_dir = OUT_DIR / run["run_id"]
    run_dir.mkdir(parents=True, exist_ok=True)

    well_spec, fluids, schedule, _ = loader()
    casing = CasingFlowSolver(
        enable_gravity=True, mixing_contact_time=True,
        plug_face_zero_mixing=True, has_plug=True)
    t0 = time.perf_counter()
    casing_result = casing.run(well_spec, fluids, schedule)
    casing_1d_wall_s = time.perf_counter() - t0
    stop_t = annulus_stop_time_s(casing_result=casing_result, fluids=fluids)
    provider = build_coupled_annulus_inlet_provider(
        casing_result, casing, fluids, split_cement_phases=True)

    cement_names = {f.name for f in fluids if f.role.name in _CEMENT_ROLE_NAMES}
    v_cement_m3 = float(sum(s.volume_m3 for s in schedule.steps
                            if s.fluid_name in cement_names))
    v_annulus_m3 = annulus_volume_of(well_spec)

    spy = NonlinearSpy()
    with spy:
        solver = AnnulusD2DGASolver(total_t=stop_t, nz=250, **dict(run["kwargs"]))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            t0 = time.perf_counter()
            result = solver.run(well_spec, fluids, provider, schedule=schedule)
            wall_s = time.perf_counter() - t0

    diag = spy.aggregate()
    warn = classify_warnings(list(caught))
    metrics = result.metrics
    dts = metrics["time_s"].diff().dropna().to_numpy(dtype=float)
    summary = result.summary
    final = summary["最终结果"]
    eta_e = float(final["全井段最终有效顶替效率"])
    eta_n = float(summary["eta_narrow"])
    cement_occ = float(final["最终水泥浆占据率"])
    mean_wall = float(np.mean(result.wall_field))
    inventory_ratio = (cement_occ * v_annulus_m3 / v_cement_m3) if v_cement_m3 > 0 else None

    status = "CONTAMINATED" if warn["fallback_count"] > 0 else "DONE"
    record: dict[str, Any] = {
        "run_id": run["run_id"],
        "well": well_key,
        "labels": list(run["labels"]),
        "caliber": run["caliber"],
        "annotation": run["annotation"],
        "solver_kwargs": {k: (v if not isinstance(v, dict) else dict(v))
                          for k, v in run["kwargs"].items()},
        "status": status,
        "total_t_s": round(stop_t, 3),
        "casing_1d_wall_s": round(casing_1d_wall_s, 2),
        "wall_time_s": round(wall_s, 2),
        "n_steps": int(len(metrics)),
        "dt_min_s": float(dts.min()) if dts.size else None,
        "dt_median_s": float(np.median(dts)) if dts.size else None,
        "dt_max_s": float(dts.max()) if dts.size else None,
        "eta_E": eta_e,
        "eta_N": eta_n,
        "cement_occ": cement_occ,
        "mean_wall": mean_wall,
        "channeling_index": float(final["最终窜槽指数"]),
        "mixing_index": float(final["最终混浆指数"]),
        "buoyancy_number": float(final["浮力数_b"]),
        "inventory_ratio": inventory_ratio,
        "annulus_volume_m3": v_annulus_m3,
        "cement_volume_m3": v_cement_m3,
        "diagnostics": diag,
        "warnings": warn,
        "digests": {
            "cement": _digest(result.cement_field),
            "spacer": _digest(result.spacer_field),
            "wall": _digest(result.wall_field),
            "lead": _digest(result.lead_field),
            "tail": _digest(result.tail_field),
        },
        "output_dir": str(run_dir),
        "finished_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    # 逐 run 产物：摘要 JSON（全量）+ 时间序列 CSV（不做云图/NPZ/GIF，控制体积）
    payload = {
        "run": {k: v for k, v in record.items() if k != "output_dir"},
        "solver_summary": summary,
    }
    (run_dir / "摘要.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=1, default=float), encoding="utf-8")
    metrics.to_csv(run_dir / "时间序列结果.csv", index=False, encoding="utf-8-sig")
    return record


# ---------------------------------------------------------------------------
# manifest / 汇总 / 归因 / L1
# ---------------------------------------------------------------------------
def _git_head() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT,
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception as exc:  # pragma: no cover - 环境异常时留痕不中断
        return f"UNAVAILABLE: {exc}"


def new_manifest(wells: list[str]) -> dict[str, Any]:
    wells_meta: dict[str, Any] = {}
    for key in wells:
        plan, calibers = build_run_plan(key)
        loader = WELLS[key]
        well_spec, fluids, _schedule, _ = loader()
        fluids_audit = [{
            "name": f.name, "role": f.role.name,
            "rheology_model": getattr(f.rheology_model, "name", str(f.rheology_model)),
            "power_law_n": f.power_law_n, "consistency_k": f.consistency_k,
            "yield_stress_pa": f.yield_stress_pa,
        } for f in fluids]
        wells_meta[key] = {
            "tau_y_calibers": calibers,
            "run_plan": [{k: r[k] for k in ("run_id", "labels", "caliber", "annotation")}
                         for r in plan],
            "fluids_audit": fluids_audit,
            "well_name": well_spec.well_name,
        }
    return {
        "meta": {
            "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "script": "scripts/entrypoints/run_hb_ablation_phase_a.py",
            "git_head": _git_head(),
            "python": sys.version,
            "numpy": np.__version__,
            "platform": platform.platform(),
            "conda_env": os.environ.get("CONDA_DEFAULT_ENV", ""),
            "matrix": "7 行 × 4 井（用户裁定 2026-09-17）；ht1_004 按 R-T6-6 spec 优先 a/b 重合",
            "runner_reuse": "load_<well>_tailpipe + CasingFlowSolver(T1 生产口径) + "
                            "annulus_stop_time_s(F2) + split_cement_phases=True + nz=250；"
                            "唯一差别 = 变体开关",
        },
        "wells": wells_meta,
        "runs": {},
    }


def load_manifest() -> dict[str, Any] | None:
    if MANIFEST_PATH.exists():
        return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    return None


def save_manifest(manifest: dict[str, Any]) -> None:
    tmp = MANIFEST_PATH.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(manifest, ensure_ascii=False, indent=1, default=float),
                   encoding="utf-8")
    os.replace(tmp, MANIFEST_PATH)


def write_summary_csv(manifest: dict[str, Any]) -> None:
    """逐 (井, 标签) 一行；ht1_004 的双标签行共享 run_id 与指标。"""
    columns = ["well", "label", "run_id", "caliber", "status", "annotation",
               "enable_stream_yield_gate", "enable_hb_closure", "hb_fix_cement_tau_y",
               "tau_y_map", "fallback_count", "fallback_steps", "nonlinear_steps",
               "max_n_static", "all_undefined_steps", "floor_warning_count",
               "n_steps", "total_t_s", "wall_time_s",
               "eta_E", "eta_N", "cement_occ", "mean_wall", "inventory_ratio",
               "buoyancy_number", "digest_cement", "digest_wall", "output_dir"]
    rows: list[dict[str, Any]] = []
    for run_id, rec in manifest["runs"].items():
        kwargs = rec["solver_kwargs"]
        for label in rec["labels"]:
            rows.append({
                "well": rec["well"], "label": label, "run_id": run_id,
                "caliber": rec["caliber"], "status": rec["status"],
                "annotation": rec["annotation"],
                "enable_stream_yield_gate": kwargs.get("enable_stream_yield_gate", False),
                "enable_hb_closure": kwargs.get("enable_hb_closure", False),
                "hb_fix_cement_tau_y": kwargs.get("hb_fix_cement_tau_y", False),
                "tau_y_map": json.dumps(kwargs.get("cement_tau_y_by_role", {}),
                                        ensure_ascii=False),
                "fallback_count": rec["warnings"]["fallback_count"],
                "fallback_steps": rec["diagnostics"]["fallback_steps"],
                "nonlinear_steps": rec["diagnostics"]["nonlinear_steps"],
                "max_n_static": rec["diagnostics"]["max_n_static"],
                "all_undefined_steps": rec["diagnostics"]["all_undefined_steps"],
                "floor_warning_count": rec["warnings"]["floor_warning_count"],
                "n_steps": rec["n_steps"], "total_t_s": rec["total_t_s"],
                "wall_time_s": rec["wall_time_s"],
                "eta_E": rec["eta_E"], "eta_N": rec["eta_N"],
                "cement_occ": rec["cement_occ"], "mean_wall": rec["mean_wall"],
                "inventory_ratio": rec["inventory_ratio"],
                "buoyancy_number": rec["buoyancy_number"],
                "digest_cement": rec["digests"]["cement"],
                "digest_wall": rec["digests"]["wall"],
                "output_dir": rec["output_dir"],
            })
    with SUMMARY_CSV.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


_LABEL_DELTA = {
    "B2_op": ("H0p", "H0"),
    "closure_H1": ("H1", "H0p"),
    "gate_LS_H2a": ("H2a_LS", "H0p"),
    "gate_HB_H2b": ("H2b_HB", "H0p"),
    "both_LS_H3a": ("H3a_LS", "H0p"),
    "both_HB_H3b": ("H3b_HB", "H0p"),
    "interaction_LS": ("H3a_LS", "H1", "H2a_LS", "H0p"),   # H3−H1−H2+H0p
    "interaction_HB": ("H3b_HB", "H1", "H2b_HB", "H0p"),
}
_METRIC_KEYS = ("eta_E", "eta_N", "cement_occ")


def write_attribution_csv(manifest: dict[str, Any]) -> None:
    """spec §5.1 三分量分解（按标签取数；缺 run 记 NaN 不编造）。"""
    columns = ["well", "component", "eta_E", "eta_N", "cement_occ", "missing_labels"]
    rows: list[dict[str, Any]] = []
    wells = sorted({rec["well"] for rec in manifest["runs"].values()})
    for well in wells:
        by_label: dict[str, dict[str, Any]] = {}
        for rec in manifest["runs"].values():
            if rec["well"] != well:
                continue
            for label in rec["labels"]:
                by_label[label] = rec
        for component, labels in _LABEL_DELTA.items():
            values: list[float | None] = []
            missing = [lb for lb in labels if lb not in by_label]
            for metric in _METRIC_KEYS:
                if missing:
                    values.append(None)
                    continue
                terms = [by_label[lb][metric] for lb in labels]
                total = terms[0] - terms[1]
                if len(labels) == 4:  # 交互项 H3−H1−H2+H0p
                    total = terms[0] - terms[1] - terms[2] + terms[3]
                values.append(total)
            rows.append({"well": well, "component": component,
                         "eta_E": values[0], "eta_N": values[1],
                         "cement_occ": values[2], "missing_labels": ";".join(missing)})
    with ATTRIBUTION_CSV.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


# 既有权威结果（信息性对照用；R18 权威目录只读不写）
_AUTHORITATIVE_SUMMARY = {
    "hu101": PROJECT_ROOT / "results" / "呼101尾管_1D2D耦合模型"
             / "呼101尾管_1D2D耦合模型_结果摘要.json",
    "hu102": PROJECT_ROOT / "results" / "呼102尾管_1D2D耦合模型"
             / "呼102尾管_1D2D耦合模型_结果摘要.json",
    "ht1_003": PROJECT_ROOT / "results" / "呼1-003_1D2D耦合模型"
               / "呼1-003_1D2D耦合模型_结果摘要.json",
    "ht1_004": PROJECT_ROOT / "results" / "呼1-004_1D2D耦合模型"
               / "呼1-004_1D2D耦合模型_结果摘要.json",
}


def write_l1_determinism(manifest: dict[str, Any]) -> dict[str, Any]:
    """L1：H0 与 H0_repeat 同配置两次运行逐位对比（shenjingwangluo 域内）。

    并与该井既有权威 runner 摘要做**信息性**对照（只报漂移、不追查——代码默认
    路径逐位已有 Task 6 探针证据）。
    """
    out: dict[str, Any] = {}
    for well in WELLS:
        first = manifest["runs"].get(f"{well}__H0")
        second = manifest["runs"].get(f"{well}__H0_repeat")
        if first is None or second is None:
            out[well] = {"status": "PENDING"}
            continue
        entry: dict[str, Any] = {
            "status": "DONE",
            "bitwise_equal": first["digests"] == second["digests"],
            "digests_first": first["digests"],
            "digests_repeat": second["digests"],
            "eta_E_first": first["eta_E"], "eta_E_repeat": second["eta_E"],
            "eta_N_first": first["eta_N"], "eta_N_repeat": second["eta_N"],
        }
        reference = _AUTHORITATIVE_SUMMARY.get(well)
        if reference is not None and reference.exists():
            ref = json.loads(reference.read_text(encoding="utf-8"))
            ref_final = ref.get("最终结果", {})
            entry["authoritative_reference"] = {
                "path": str(reference),
                "eta_E": ref_final.get("全井段最终有效顶替效率"),
                "eta_N": ref.get("eta_narrow"),
                "drift_eta_E_pp": (round((first["eta_E"] - ref_final["全井段最终有效顶替效率"])
                                         * 100.0, 3)
                                   if ref_final.get("全井段最终有效顶替效率") is not None else None),
                "drift_eta_N_pp": (round((first["eta_N"] - ref["eta_narrow"]) * 100.0, 3)
                                   if ref.get("eta_narrow") is not None else None),
                "note": "信息性对照：权威摘要的产生时点/口径可能与本次 H0 不同，只报漂移不追查",
            }
        out[well] = entry
    path = OUT_DIR / "l1_determinism.json"
    path.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    return out


# ---------------------------------------------------------------------------
# L3：Z&F22 基准算例 H0 与 H3b（最低集）
# ---------------------------------------------------------------------------
def run_l3() -> None:
    """基准流体为纯牛顿（FluidSpec 默认 NEWTONIAN、无 yield_stress）⇒ HB 闭包
    结构性短路为牛顿恒等（R2），L3 实为回归检查——如实写明。

    H3b 口径：``cement_tau_y_by_role={"TAIL": 0.0}`` 是基准水泥相的**真实** τy
    （无屈服），非编造值；R-T6-2 的 τ_Y1 取泥浆 spec YP = None ⇒ 0。
    """
    L3_DIR.mkdir(parents=True, exist_ok=True)
    variants = {
        "H0": {},
        "H3b_HB": {"enable_stream_yield_gate": True, "enable_hb_closure": True,
                   "hb_fix_cement_tau_y": True, "cement_tau_y_by_role": {"TAIL": 0.0}},
    }
    all_rows: dict[str, list[dict[str, Any]]] = {}
    for name, overrides in variants.items():
        rows = []
        for case in ZHANG2022_CASES:
            cid = int(case["case_id"])
            print(f"[L3] {name} case {cid} ...", flush=True)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                row = run_case(cid, solver_overrides=dict(overrides), keep_result=True)
            row["fallback_count"] = sum(
                1 for w in caught
                if issubclass(w.category, RuntimeWarning) and "回退牛顿线性闭包" in str(w.message))
            result = row.pop("_result")
            row["digest_cement"] = _digest(result.cement_field)
            row["digest_wall"] = _digest(result.wall_field)
            row["digest_spacer"] = _digest(result.spacer_field)
            rows.append(row)
        all_rows[name] = rows
        path = L3_DIR / f"对照表_{name}.csv"
        with path.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()),
                                    extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)

    def _thresholds(rows: list[dict[str, Any]]) -> dict[str, Any]:
        # 三门槛判据 = 台账既有口径（progress.md Task 12 行）：
        # ① |Δη_E|≤0.05 命中数（现 5/10）；② mass<0.05 命中数（现原口径 1/10）；
        # ③ e=0.8→0.1 跨度（现 75.31pp = mean(η_E@e=0.1) − min(η_E@e=0.8)，与台账逐位对上）。
        hit_eta = sum(1 for r in rows if abs(r["eta_E_偏差_vs_D2DGA"]) <= 0.05)
        hit_mass = sum(1 for r in rows if r["mass_conservation_error"] < 0.05)
        e_low = [r for r in rows if abs(r["e"] - 0.1) < 1e-9]
        e_high = [r for r in rows if abs(r["e"] - 0.8) < 1e-9]
        span = (float(np.mean([r["eta_E_模型"] for r in e_low])
                      - min(r["eta_E_模型"] for r in e_high))
                if e_low and e_high else float("nan"))
        return {"eta_E_within_0.05_count": hit_eta,
                "mass_lt_0.05_count": hit_mass,
                "span_e0.8_to_e0.1_pp": round(span * 100.0, 2),
                "thresholds_met": int(hit_eta >= 5) + int(hit_mass >= 10) + int(span * 100.0 >= 75.31)}

    payload = {
        "说明": "基准流体纯牛顿 ⇒ HB 闭包结构性短路为牛顿恒等（R2），L3=回归检查；"
                "H3b 的 cement_tau_y_by_role={'TAIL': 0.0} 为基准水泥相真实 τy（无屈服）",
        "H0": {"thresholds": _thresholds(all_rows["H0"]), "rows": all_rows["H0"]},
        "H3b_HB": {"thresholds": _thresholds(all_rows["H3b_HB"]), "rows": all_rows["H3b_HB"]},
        "bitwise_equal_H3b_vs_H0": all(
            a["digest_cement"] == b["digest_cement"]
            and a["digest_spacer"] == b["digest_spacer"]
            and a["digest_wall"] == b["digest_wall"]
            for a, b in zip(all_rows["H0"], all_rows["H3b_HB"])),
    }
    baseline = PROJECT_ROOT / "results" / "基准算例对照_2026-09-14" / "对照表.csv"
    if baseline.exists():
        with baseline.open(encoding="utf-8-sig") as handle:
            base_rows = list(csv.DictReader(handle))
        base_rows = [{**r, "e": float(r["e"]),
                      "eta_E_偏差_vs_D2DGA": (float(r["eta_E_模型"])
                                              - float(r["eta_E_论文D2DGA"])),
                      "mass_conservation_error": float(r["mass_conservation_error"]),
                      "eta_E_模型": float(r["eta_E_模型"])} for r in base_rows]
        payload["baseline_2026-09-14"] = {
            "path": str(baseline),
            "thresholds": _thresholds(base_rows),
            "note": "既有权威基准（λ_op 修正后口径，Task 9/12 线；同一判据函数复算），"
                    "供'不劣于现 1/3'对照",
        }
    (L3_DIR / "l3_summary.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=1, default=float), encoding="utf-8")
    print("[L3] thresholds H0:", json.dumps(payload["H0"]["thresholds"], ensure_ascii=False))
    print("[L3] thresholds H3b:", json.dumps(payload["H3b_HB"]["thresholds"], ensure_ascii=False))
    if "baseline_2026-09-14" in payload:
        print("[L3] baseline 2026-09-14:",
              json.dumps(payload["baseline_2026-09-14"]["thresholds"], ensure_ascii=False))
    print("[L3] H3b ≡ H0 (bitwise on key metrics):", payload["bitwise_equal_H3b_vs_H0"])


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--wells", type=str,
                        default="hu101,hu102,ht1_003,ht1_004",
                        help="逗号分隔井键（注册表键名，非中文井名）")
    parser.add_argument("--variants", type=str, default="",
                        help="逗号分隔标签过滤（如 H0,H0p）；空 = 全部")
    parser.add_argument("--rerun", type=str, default="",
                        help="逗号分隔 run_id：强制重跑（忽略已完成状态）")
    parser.add_argument("--list", action="store_true", help="只打印 run 计划后退出")
    parser.add_argument("--l3", action="store_true",
                        help="只跑 L3（zhang2022 基准 H0+H3b）后退出")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if args.list:
        for well in [w.strip() for w in args.wells.split(",") if w.strip()]:
            plan, calibers = build_run_plan(well)
            print(f"== {well} ==")
            for run in plan:
                print(f"  {run['run_id']:<24} labels={','.join(run['labels']):<16} "
                      f"caliber={run['caliber']:<14} {run['annotation']}")
            print(f"  τy 口径: bingham_ls={calibers['bingham_ls']} "
                  f"hb_fitted_new={calibers['hb_fitted_new']}")
        return
    if args.l3:
        run_l3()
        return

    wells = [w.strip() for w in args.wells.split(",") if w.strip()]
    for well in wells:
        if well not in WELLS:
            raise SystemExit(f"未知井键 {well!r}；可选 {sorted(WELLS)}")
    force = {r.strip() for r in args.rerun.split(",") if r.strip()}
    label_filter = ({v.strip() for v in args.variants.split(",") if v.strip()}
                    or None)

    manifest = load_manifest()
    if manifest is None:
        manifest = new_manifest(wells)
        save_manifest(manifest)
        print(f"[manifest] 新建 {MANIFEST_PATH}", flush=True)
    else:
        print(f"[manifest] 续跑：已有 {len(manifest['runs'])} 个 run 记录", flush=True)

    total = 0
    for well in wells:
        plan, _calibers = build_run_plan(well)
        if well not in manifest["wells"]:
            fresh = new_manifest([well])
            manifest["wells"][well] = fresh["wells"][well]
        for run in plan:
            if label_filter is not None and not (set(run["labels"]) & label_filter):
                continue
            total += 1
            existing = manifest["runs"].get(run["run_id"])
            if existing is not None and existing["status"] in ("DONE", "CONTAMINATED") \
                    and run["run_id"] not in force:
                print(f"[skip] {run['run_id']}（已完成 {existing['status']}）", flush=True)
                continue
            print(f"[run ] {run['run_id']} 开始：{run['annotation']}", flush=True)
            t0 = time.perf_counter()
            try:
                record = execute_run(run, well, WELLS[well])
            except Exception as exc:  # 失败留痕不中断矩阵（续跑可重试）
                manifest["runs"][run["run_id"]] = {
                    "run_id": run["run_id"], "well": well, "labels": run["labels"],
                    "caliber": run["caliber"], "annotation": run["annotation"],
                    "solver_kwargs": {k: (v if not isinstance(v, dict) else dict(v))
                                      for k, v in run["kwargs"].items()},
                    "status": "FAILED", "error": f"{type(exc).__name__}: {exc}",
                    "finished_at": time.strftime("%Y-%m-%d %H:%M:%S"),
                }
                save_manifest(manifest)
                print(f"[FAIL] {run['run_id']}: {type(exc).__name__}: {exc}", flush=True)
                continue
            manifest["runs"][run["run_id"]] = record
            save_manifest(manifest)
            write_summary_csv(manifest)
            write_attribution_csv(manifest)
            write_l1_determinism(manifest)
            print(
                f"[done] {run['run_id']} {record['status']} "
                f"eta_E={record['eta_E']:.4f} eta_N={record['eta_N']:.4f} "
                f"steps={record['n_steps']} fallback={record['warnings']['fallback_count']} "
                f"wall={record['wall_time_s']:.0f}s (总 {time.perf_counter() - t0:.0f}s)",
                flush=True)
    print(f"[plan] 本轮计划 run 数 = {total}", flush=True)


if __name__ == "__main__":
    main()
