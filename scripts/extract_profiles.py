"""
快速提取窄边水泥浓度剖面和鞋口出流时序（Step 3.3 & 3.4）

baseline / improved 两组的边界口径（2026-08-22 弥散接线生产后）：
- baseline：屈服正则 M=0、轴向弥散关。弥散关闭时鞋口边界与体积追踪
  单相口径等价（shoe_timeline 与 pipe_exit_state_at 一致）；
- improved：屈服正则 M=100、轴向弥散开。鞋口边界为多相过渡带
  （Taylor-Aris 弥散 + 界面混浆增强），环空入口浓度边界被平滑。
两组差异同时含 M 与弥散两因素，归因单一机制需分组消融。
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from cemdisp.data.loaders import load_hu101_tailpipe, load_hu102_tailpipe
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver


def extract_profiles(well_key: str):
    loader = load_hu101_tailpipe if well_key == "hu101" else load_hu102_tailpipe
    well_spec, fluids, schedule, _ = loader()

    results = {}
    for label, M, enable_disp, alpha in [
        ("baseline", 0.0, False, 0.2),
        ("improved", 100.0, True, 0.2),
    ]:
        casing_solver = CasingFlowSolver(
            enable_gravity=True,
            enable_axial_dispersion=enable_disp,
            dispersion_alpha=alpha,
        )
        casing_result = casing_solver.run(well_spec, fluids, schedule)
        provider = build_coupled_annulus_inlet_provider(
            casing_result, casing_solver, fluids, split_cement_phases=True,
        )

        # 鞋口出流时序
        shoe_events = []
        for ev in casing_result.shoe_timeline.events:
            shoe_events.append({
                "time_s": ev.time_s,
                "kind": ev.kind.name if hasattr(ev.kind, "name") else str(ev.kind),
                "stage_name": ev.stage_name,
                "phase_fractions": ev.phase_fractions,
            })

        total_t = float(casing_result.cement_end_time_s)
        solver = AnnulusD2DGASolver(total_t=total_t, nz=500, yield_regularization_M=M)
        result = solver.run(well_spec, fluids, provider)

        dp = result.depth_profiles
        profile = {
            "depth_m": dp["井深_m"].to_list(),
            "narrow_cement": dp["窄边水泥浓度"].to_list(),
            "wide_cement": dp["宽边水泥浓度"].to_list(),
            "narrow_eff": dp["窄边有效效率"].to_list(),
            "wide_eff": dp["宽边有效效率"].to_list(),
        }

        results[label] = {
            "profile": profile,
            "shoe_events": shoe_events,
            "cbl_eval_eff": float(result.metrics.iloc[-1]["cbl_eval_interval_efficiency"]),
        }

    # 保存
    out_path = _PROJECT_ROOT / "results" / "p3_p4_integration" / f"{well_key}_profiles.json"
    out_path.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Saved {out_path}")


if __name__ == "__main__":
    for well in ["hu101", "hu102"]:
        print(f"\n=== Extracting {well} ===")
        extract_profiles(well)
