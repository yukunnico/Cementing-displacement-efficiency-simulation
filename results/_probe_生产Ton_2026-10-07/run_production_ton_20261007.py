#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""生产口径「最终模拟结果」驱动（D8-A 镜像脚本，**不改 runner 本体**）

定位与红线
----------
用户 2026-09-09 裁定：**权威结果口径 = `cemdisp/runners/` 产出的 `results/<井>_1D2D耦合模型/`**。
D8 裁定（2026-10-06）：runner 链 T-on 走 **A = 镜像脚本试点（零改码）**；「产品化（改 8 井
runner 构造）」= D8-B，**等另裁**。本脚本即 A 路线的产品化延伸：

- **逐行复刻 runner 口径**：`CasingFlowSolver(enable_gravity=True, mixing_contact_time=True,
  plug_face_zero_mixing=True, has_plug=True)` + `build_coupled_annulus_inlet_provider(...,
  split_cement_phases=True)` + `AnnulusD2DGASolver(total_t=annulus_stop_time_s, nz=250)`，
  **不加 CORRECTED_KW**（runner 本身不用）。
- **在其上叠加唯一系统性差异 = 温压耦合开启**（`enable_temperature_rheology=True`、
  `temperature_mode="anchored"`、`pressure_mode="hydrostatic"`、`enable_depthwise_temperature=True`）。
- **产物写新日期目录** `results/<权威目录名>_Ton_2026-10-07/`，**绝不触碰权威目录**。

井位：**七井**（八井 − 呼103）——呼103 按裁定「退出重点井、新批一律不跑」，其 k 亦未裁定，
故不在本批；呼探1 无锚 ⇒ Geothermal 回退（`ANCHORED_NO_ANCHOR_WELLS`）。

用法：PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -I run_production_ton_20261007.py --batch
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
    annulus_kwargs_from_opts,
    build_pressure_field,
    build_temperature_fields,
    casing_kwargs_from_opts,
    extra_metrics,
    normalize_run_opts,
)
from cemdisp.models2d import AnnulusD2DGASolver  # noqa: E402
from cemdisp.models2d.boundary_bridge import (  # noqa: E402
    build_coupled_annulus_inlet_provider,
)
from cemdisp.transport1d import CasingFlowSolver  # noqa: E402

OUT = HERE.parent
OUT_ROOT = BRANCH_ROOT / "results"

# 七井：canonical id → (中文名, loader 模块, loader 函数, runner 模块[取 annulus_stop_time_s],
#                       权威 T-off 目录名)
WELLS = {
    "hu101": ("呼101", "cemdisp.data.loaders.hu101_loader", "load_hu101_tailpipe",
              "cemdisp.runners.hu101_tailpipe", "呼101尾管_1D2D耦合模型"),
    "ht1_003": ("呼1-003", "cemdisp.data.loaders.ht1_003_loader", "load_ht1_003_tailpipe",
                "cemdisp.runners.ht1_003_tailpipe", "呼1-003_1D2D耦合模型"),
    "ht1_004": ("呼1-004", "cemdisp.data.loaders.ht1_004_loader", "load_ht1_004_tailpipe",
                "cemdisp.runners.ht1_004_tailpipe", "呼1-004_1D2D耦合模型"),
    "hu102": ("呼102", "cemdisp.data.loaders.hu102_loader", "load_hu102_tailpipe",
              "cemdisp.runners.hu102_tailpipe", "呼102尾管_1D2D耦合模型"),
    "hu1": ("呼探1", "cemdisp.data.loaders.hu1_loader", "load_hu1_tailpipe",
            "cemdisp.runners.hu1_tailpipe", "呼探1尾管_1D2D耦合模型"),
    "hu2": ("呼探1-002", "cemdisp.data.loaders.hu2_loader", "load_hu2_tailpipe",
            "cemdisp.runners.hu2_tailpipe", "呼探1-002尾管_1D2D耦合模型"),
    "ht1_001": ("呼探1-001", "cemdisp.data.loaders.ht1_001_loader", "load_ht1_001_tailpipe",
                "cemdisp.runners.ht1_001_tailpipe", "呼探1-001尾管_1D2D耦合模型"),
}


def authoritative_teff(wk: str) -> dict | None:
    """读权威目录（runner 产出）里的 T-off 结果摘要（**读文件，禁手抄**）。"""
    d = OUT_ROOT / WELLS[wk][4]
    if not d.is_dir():
        return None
    for p in sorted(d.glob("*_结果摘要.json")):
        try:
            r = json.loads(p.read_text(encoding="utf-8"))["最终结果"]
            return {"η_E": float(r["全井段最终有效顶替效率"]),
                    "η_N": float(r["窄四分位效率"]), "source": p.name}
        except Exception:
            continue
    return None


def run_well(wk: str) -> dict:
    cn, mod, fn, runner_mod, auth_dir = WELLS[wk]
    loader = getattr(importlib.import_module(mod), fn)
    stop_fn = getattr(importlib.import_module(runner_mod), "annulus_stop_time_s")

    opts = normalize_run_opts({
        "enable_temperature_rheology": True,
        "temperature_mode": "anchored",
        "enable_yield_gate": None,          # 沿用 annulus 构造默认（runner 口径）
        "pressure_mode": "hydrostatic",
        "pressure_caliber": "shoe",
        "include_yield_term": False,
        "enable_stream_yield_gate": False,
        "enable_depthwise_temperature": True,
    })

    t0 = time.perf_counter()
    well, fluids, schedule, _ = loader()
    field_1d, field_2d, note = build_temperature_fields(cn, "anchored")
    # 压力场：默认 hydrostatic；若该井井身 md 网格非严格升序（实测呼探1）则**降级 off
    # 并把原因记进结果**（不掩盖、不静默）——温度部分不受影响。
    pressure_note = "hydrostatic"
    p_field = None
    try:
        p_field = build_pressure_field(well, fluids, schedule, "hydrostatic")
    except ValueError as exc:
        pressure_note = f"降级 off（hydrostatic 建场失败：{exc}）"
        opts["pressure_mode"] = "off"
        p_field = None

    # ---- 逐行复刻 runner 的 1D 构造（仅叠加温度/压力场）----
    casing = CasingFlowSolver(
        enable_gravity=True, mixing_contact_time=True,
        plug_face_zero_mixing=True, has_plug=True,
        pressure_field=p_field, **casing_kwargs_from_opts(opts))
    cr = casing.run(well, fluids, schedule, temperature_field=field_1d)
    inlet = build_coupled_annulus_inlet_provider(
        cr, casing, fluids, split_cement_phases=True)
    stop_t = float(stop_fn(casing_result=cr, fluids=fluids))

    # ---- 逐行复刻 runner 的 2D 构造（nz=250 纯默认 + 温度/压力/逐列）----
    # 逐列开关由 `annulus_kwargs_from_opts(opts)` 出键（5d-③ 已接线）——勿再显式传，
    # 否则 `**` 展开报 "multiple values for keyword argument"。
    solver = AnnulusD2DGASolver(
        total_t=stop_t, nz=250,
        pressure_field=p_field, **annulus_kwargs_from_opts(opts))
    res = solver.run(well, fluids, inlet, schedule=schedule,
                     temperature_field=field_2d)
    elapsed = round(time.perf_counter() - t0, 1)

    extra = extra_metrics(res, res.summary, field_1d, field_2d)
    final = res.summary["最终结果"]
    row = {
        "井": cn, "canonical": wk,
        "温度档": opts["temperature_mode"], "场备注": note,
        "压力档": pressure_note, "2D逐列": True,
        "nz": 250, "stop_t_s": stop_t,
        "η_E": float(final["全井段最终有效顶替效率"]),
        "η_N": float(final["窄四分位效率"]),
        "饥饿份额": extra.get("饥饿份额"),
        "屈服门_wall占比": extra.get("屈服门_wall占比"),
        "浮力数_b": final.get("浮力数_b"),
        "步数": int(np.size(res.dt_history)),
        "cfl_clip_events": int(res.cfl_clip_events),
        "col_memo_size": len(getattr(solver, "_col_memo", {})),
        "温度审计_oob": extra.get("温度审计_oob"),
        "耗时_s": elapsed,
        "权威T_off": authoritative_teff(wk),
        "result_json": None,
    }
    d = OUT_ROOT / f"{auth_dir}_Ton_2026-10-07"
    d.mkdir(parents=True, exist_ok=True)
    p = d / f"{cn}_Ton全链_结果摘要.json"
    p.write_text(json.dumps(row, ensure_ascii=False, indent=1, default=str),
                 encoding="utf-8")
    row["result_json"] = str(p.relative_to(BRANCH_ROOT))
    (OUT / "jobs").mkdir(parents=True, exist_ok=True)
    (OUT / "jobs" / f"{wk}.json").write_text(
        json.dumps(row, ensure_ascii=False, indent=1, default=str), encoding="utf-8")
    print(f"[{cn}] etaE={row['η_E']:.6f} etaN={row['η_N']:.6f} "
          f"stop_t={stop_t:.0f}s elapsed={elapsed}s", flush=True)
    return row


def run_batch(which: list[str], workers: int) -> int:
    import os
    (OUT / "jobs").mkdir(parents=True, exist_ok=True)
    pending = list(which); running = []; failed = []; done = []
    env = dict(os.environ); env["PYTHONIOENCODING"] = "utf-8"; env["PYTHONUTF8"] = "1"
    py = sys.executable
    while pending or running:
        while pending and len(running) < workers:
            wk = pending.pop(0)
            jf = OUT / "jobs" / f"_spec_{wk}.json"
            jf.write_text(json.dumps({"well": wk}), encoding="utf-8")
            logf = open(OUT / "jobs" / f"{wk}.log", "w", encoding="utf-8")
            p = subprocess.Popen([py, str(HERE), "--well", wk],
                                 stdout=logf, stderr=subprocess.STDOUT, env=env,
                                 cwd=str(BRANCH_ROOT))
            running.append((p, wk, logf))
            print(f"[launcher] 启动 {wk} (pid={p.pid})", flush=True)
        time.sleep(2.0)
        still = []
        for p, wk, logf in running:
            rc = p.poll()
            if rc is None:
                still.append((p, wk, logf)); continue
            logf.close()
            if rc == 0:
                done.append(wk)
            else:
                failed.append((wk, rc))
            print(f"[launcher] {'完成' if rc == 0 else 'FAIL'} {wk} exit={rc}", flush=True)
        running = still
    print(f"[launcher] 批完成：{len(done)} 成功 / {len(failed)} 失败 {failed}", flush=True)
    return 0 if not failed else 1


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--well", type=str, default=None)
    ap.add_argument("--batch", action="store_true")
    ap.add_argument("--workers", type=int, default=3)
    args = ap.parse_args()
    if args.well:
        run_well(args.well)
    elif args.batch:
        sys.exit(run_batch(list(WELLS), args.workers))
    else:
        ap.error("需 --well 或 --batch")


if __name__ == "__main__":
    main()
