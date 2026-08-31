"""F1+F2 弥散量级与屈服门槛安全系数标定扫描（8 井）。

spec: docs/superpowers/specs/2026-08-31-dispersion-magnitude-calibration.md

背景（spec 侦察结论）：
- 当前 _smooth_dispersion 常数 0.018/0.015 隐含 D_eff ≈ 1e-3–5e-3 m²/s，恰为湍流
  Taylor 量级（Taylor 1954: D_T = 5.05·b·u*），但 8 井湍流网格占比全为 0（全深层流，
  Bingham 泥浆 He=O(1e5–1e7) 把 re_crit 推到 1e5–1e7）→ 弥散被系统性高估 1–3 个量级；
- 层流物理弥散按 Taylor-Aris 平行板（Vedel & Bruus 2012）仅 ≈ 1e-5–1e-3 m²/s，
  湍流支公式 D_T = 5.05·b·u* 仅作适用域声明（当前 8 井不生效）；
- CBL 对照呈"好井偏低/差井偏高"的向均值回归失真，与弥散过强（抹平界面差异）的
  机理指纹一致。

扫描矩阵（CORRECTED_KW 最终基线口径上覆盖，CFL on + nz=250）：
- F1: dispersion_dt_scale ∈ {0.25, 0.5}（基线 1.0 已有；scale=0 已知病态，不扫）；
- F2: yield_gate_f_safety ∈ {1.3}（基线 1.15 已有）。
- 任务拆分：F1 任务固定 yield_gate_f_safety=基线 1.15；F2 任务固定
  dispersion_dt_scale=基线 1.0 —— 保证 F2 与 F1 独立评估（spec 执行清单：
  scale 0.25/0.5 × 8 井 + f_safety 1.3 × 8 井 = 默认 24 任务）。

产出：results/弥散产额标定_2026-08-31/ 下逐井逐配置 JSON（断点续跑）+ 扫描汇总.csv。
评价（CBL 方向一致性收敛度 S、η_N/wall_frac 方向、稳定性）由分析步骤另行完成，见 spec F1。

用法：
  python scripts/calibrate_dispersion_yield.py                     # 默认矩阵全量（24 任务）
  python scripts/calibrate_dispersion_yield.py --workers 3         # 并行进程数（默认 3）
  python scripts/calibrate_dispersion_yield.py --wells ht1_003 --scales 0.5 --f-safeties ""
                                                                   # 单井冒烟（跳过 F2）
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import cast

import numpy as np

_SCRIPTS_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SCRIPTS_DIR.parent
for _p in (str(_SCRIPTS_DIR), str(_PROJECT_ROOT)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from cemdisp.models2d import AnnulusD2DGASolver  # noqa: E402
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider  # noqa: E402
from cemdisp.transport1d import CasingFlowSolver  # noqa: E402
from rerun_all_wells_corrected import CORRECTED_KW, WELLS, _stop_t, _total_t  # noqa: E402

OUT_DIR = _PROJECT_ROOT / "results" / "弥散产额标定_2026-08-31"
JOB_DIR = OUT_DIR / "单井结果"

# 各扫描轴的基线值（另一轴固定在基线上，保证 F1/F2 独立评估）
BASELINE_SCALE = 1.0      # CORRECTED_KW 里的 dispersion_dt_scale
BASELINE_F_SAFETY = 1.15  # AnnulusD2DGASolver 默认 yield_gate_f_safety

CSV_FIELDS = ["well", "dispersion_dt_scale", "yield_gate_f_safety", "eta_E", "eta_N",
              "mixing", "channeling", "cement_occ", "wall_frac", "elapsed_s"]


def run_config(well_id: str, *, dispersion_dt_scale: float,
               yield_gate_f_safety: float) -> dict:
    """跑单井单配置（CORRECTED_KW 基础上覆盖扫描参数），返回汇总行 dict。

    注意：CORRECTED_KW 已含 dispersion_dt_scale=1.0，不能与显式关键字参数同传
    （否则 TypeError: got multiple values for keyword argument），故先拷贝再覆盖。
    """
    loader = next(l for w, l in WELLS if w == well_id)
    t0 = time.perf_counter()
    well, fluids, schedule, _ = loader()
    cr = CasingFlowSolver(enable_gravity=True).run(well, fluids, schedule)
    inlet = build_coupled_annulus_inlet_provider(
        cr, CasingFlowSolver(enable_gravity=True), fluids, split_cement_phases=True)
    tt = min(_total_t(schedule) + 1200.0, _stop_t(cr, fluids) + 600.0)
    kw = dict(CORRECTED_KW)
    kw["dispersion_dt_scale"] = dispersion_dt_scale
    kw["yield_gate_f_safety"] = yield_gate_f_safety
    res = AnnulusD2DGASolver(
        total_t=tt, nz=250, enable_cfl_adaptive=True, **kw
    ).run(well, fluids, inlet)
    fr = cast(dict, res.summary["最终结果"])
    return {
        "well": well_id,
        "dispersion_dt_scale": dispersion_dt_scale,
        "yield_gate_f_safety": yield_gate_f_safety,
        "eta_E": float(fr["全井段最终有效顶替效率"]),
        "eta_N": float(fr.get("窄四分位效率", float("nan"))),
        "mixing": float(fr["最终混浆指数"]),
        "channeling": float(fr["最终窜槽指数"]),
        "cement_occ": float(fr["最终水泥浆占据率"]),
        "wall_frac": float(np.mean(res.wall_field)) if res.wall_field is not None else float("nan"),
        "elapsed_s": round(time.perf_counter() - t0, 1),
    }


def _job_path(well_id: str, scale: float, f_safety: float) -> Path:
    return JOB_DIR / f"{well_id}_s{scale}_f{f_safety}.json"


def _load_job(well_id: str, scale: float, f_safety: float) -> dict | None:
    """已落盘 JSON 直接返回（断点续跑）；不存在或损坏返回 None。"""
    path = _job_path(well_id, scale, f_safety)
    if path.exists():
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            return None
    return None


def main() -> None:
    ap = argparse.ArgumentParser(description="F1+F2 弥散量级/屈服门槛安全系数标定扫描（8 井）")
    ap.add_argument("--scales", default="0.25,0.5", help="F1: dispersion_dt_scale 扫描值，逗号分隔")
    ap.add_argument("--f-safeties", default="1.3", help="F2: yield_gate_f_safety 扫描值，逗号分隔；空串跳过 F2")
    ap.add_argument("--wells", default="", help="井名过滤，逗号分隔；空=全部 8 井")
    ap.add_argument("--workers", type=int, default=3, help="并行进程数（默认 3）")
    args = ap.parse_args()

    scales = [float(s) for s in args.scales.split(",") if s.strip()]
    f_safeties = [float(s) for s in args.f_safeties.split(",") if s.strip()]
    wells = [w.strip() for w in args.wells.split(",") if w.strip()] or [w for w, _ in WELLS]

    # 任务拆分（另一轴固定基线，保证 F1/F2 独立评估）：
    #   F1 任务 = (scale, f_safety=基线 1.15)；F2 任务 = (scale=基线 1.0, f_safety)。
    # 默认 8 井 × {0.25,0.5} + 8 井 × {1.3} = 24 任务。同一 (w,s,f) 去重。
    seen: set[tuple[str, float, float]] = set()
    jobs: list[tuple[str, float, float]] = []
    for w, s, f in ([(w, s, BASELINE_F_SAFETY) for w in wells for s in scales]
                    + [(w, BASELINE_SCALE, f) for w in wells for f in f_safeties]):
        if (w, s, f) not in seen:
            seen.add((w, s, f))
            jobs.append((w, s, f))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    JOB_DIR.mkdir(parents=True, exist_ok=True)

    # 断点续跑：已落盘的任务直接复用
    pending, cached = [], []
    for w, s, f in jobs:
        saved = _load_job(w, s, f)
        (cached if saved is not None else pending).append(saved or (w, s, f))
    print(f"任务 {len(jobs)} 个：缓存 {len(cached)} | 待跑 {len(pending)}（workers={args.workers}）",
          flush=True)

    rows: list = list(cached)
    t_start = time.perf_counter()
    if pending:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {
                pool.submit(run_config, w, dispersion_dt_scale=s, yield_gate_f_safety=f): (w, s, f)
                for w, s, f in pending
            }
            for fut in as_completed(futures):
                w, s, f = futures[fut]
                try:
                    row = fut.result()
                except Exception as exc:  # 单井失败不拖垮整批
                    traceback.print_exc()
                    print(f"[失败] {w} s={s} f={f}: {exc}", flush=True)
                    continue
                _job_path(w, s, f).write_text(
                    json.dumps(row, ensure_ascii=False, indent=2), encoding="utf-8")
                rows.append(row)
                done = len(rows)
                print(f"[{done}/{len(jobs)}] {w} s={s} f={f}: eta_E={row['eta_E']:.4f} "
                      f"eta_N={row['eta_N']:.4f} wall={row['wall_frac']:.3f} "
                      f"({row['elapsed_s']}s) 累计 {time.perf_counter() - t_start:.0f}s",
                      flush=True)

    csv_path = OUT_DIR / "扫描汇总.csv"
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda r: (r["well"], r["dispersion_dt_scale"],
                                                     r["yield_gate_f_safety"])))
    print(f"\n汇总: {csv_path}（{len(rows)} 行）", flush=True)
    for row in sorted(rows, key=lambda r: (r["well"], r["dispersion_dt_scale"],
                                           r["yield_gate_f_safety"])):
        print(f"  {row['well']:<8} s={row['dispersion_dt_scale']:<5} "
              f"f={row['yield_gate_f_safety']:<5} eta_E={row['eta_E']:.4f} "
              f"eta_N={row['eta_N']:.4f} mix={row['mixing']:.4f}", flush=True)


if __name__ == "__main__":
    main()
