"""为 2026-09-19 组会汇报重跑三口关键井并出图。

输出目录：results/汇报出图_2026-09-19/<井名>/
**不触碰权威目录** results/<井名>_1D2D耦合模型/。

口径 = 各 runner 的生产默认口径，与 run_<well>_tailpipe_initial 完全一致：
CasingFlowSolver(enable_gravity=True, mixing_contact_time=True,
                 plug_face_zero_mixing=True, has_plug=True)
+ build_coupled_annulus_inlet_provider(split_cement_phases=True)
+ total_t_s = 该井 annulus_stop_time_s 口径。

唯一差异：关闭 GIF 动画生成（animate_cement_field 置空）。
纯省时，不写任何数值结果、不改任何数值口径。

选择的三口井：
  hu101    强偏心井，CBL 靶值 62.77%，全套里故事最强
  ht1_003  重构前 η_E 最高（0.9973），CBL 靶值 78.7%
  ht1_004  设计井，CBL 靶值 0.3%（界面胶结问题典型）
"""

from __future__ import annotations

import json
import sys
import time
import traceback
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.runners import ht1_003_tailpipe, ht1_004_tailpipe, hu101_tailpipe
from cemdisp.transport1d import CasingFlowSolver

OUT_ROOT = PROJECT_ROOT / "results" / "汇报出图_2026-09-19"

WELLS = [
    ("hu101", hu101_tailpipe, hu101_tailpipe.load_hu101_tailpipe),
    ("ht1_003", ht1_003_tailpipe, ht1_003_tailpipe.load_ht1_003_tailpipe),
    ("ht1_004", ht1_004_tailpipe, ht1_004_tailpipe.load_ht1_004_tailpipe),
]


def _disable_animation(module) -> None:
    """关闭 GIF 动画生成：只影响一个附加产物，不触及任何数值路径。"""
    module.animate_cement_field = lambda *a, **k: None


def run_one(name: str, module, loader) -> dict:
    out_dir = OUT_ROOT / name
    out_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    well_spec, fluids, schedule, _ = loader()

    casing_solver = CasingFlowSolver(
        enable_gravity=True,
        mixing_contact_time=True,
        plug_face_zero_mixing=True,
        has_plug=True,
    )
    casing_result = casing_solver.run(well_spec, fluids, schedule)
    provider = build_coupled_annulus_inlet_provider(
        casing_result, casing_solver, fluids, split_cement_phases=True
    )
    stop_t = module.annulus_stop_time_s(casing_result=casing_result, fluids=fluids)

    if hasattr(module, "_export_casing_flow_timing"):
        module._export_casing_flow_timing(
            output_dir=out_dir,
            schedule=schedule,
            casing_result=casing_result,
            casing_solver=casing_solver,
        )

    module.run_and_export(
        mode_title="1D2D耦合模型",
        output_dir=out_dir,
        inlet_provider=provider,
        total_t_s=stop_t,
    )
    return {"well": name, "status": "OK", "wall_s": round(time.time() - t0, 1),
            "total_t_s": float(stop_t), "out_dir": str(out_dir)}


def main() -> int:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    records = []
    for name, module, loader in WELLS:
        _disable_animation(module)
        print(f"[start] {name}", flush=True)
        try:
            rec = run_one(name, module, loader)
        except Exception:
            rec = {"well": name, "status": "FAILED", "traceback": traceback.format_exc()}
            print(rec["traceback"], flush=True)
        print(f"[done] {json.dumps(rec, ensure_ascii=False)}", flush=True)
        records.append(rec)
    (OUT_ROOT / "rerun_manifest.json").write_text(
        json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("[ALL_DONE]", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
