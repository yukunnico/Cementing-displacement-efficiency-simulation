"""C4/弥散敏感性脚本共享运行骨架：单 case = 基线同款流水线 + 可覆盖开关/井规格。

与 rerun_all_wells_corrected.run_one 同源（CORRECTED_KW + nz=250 + CFL on +
split_cement_phases + 套管 1D 重力项），仅允许：
- extra_kw：追加/覆盖 AnnulusD2DGASolver 开关（如 dispersion_dt_scale=0.0）；
- well_override：替换 WellSpec（如 standoff 剖面 ±0.1 的 dataclasses.replace 结果）。

注意：套管 1D（CasingFlowSolver）与入口桥接不吃 standoff/弥散开关，case 间
仅环空 2D 求解器输入变化，骨架保持基线同款以防口径漂移。
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import cast

import numpy as np

import sys
_SCRIPTS_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SCRIPTS_DIR.parent
for _p in (str(_SCRIPTS_DIR), str(_PROJECT_ROOT)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from cemdisp.data.well_spec import WellSpec  # noqa: E402
from cemdisp.models2d import AnnulusD2DGASolver  # noqa: E402
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider  # noqa: E402
from cemdisp.transport1d import CasingFlowSolver  # noqa: E402
from rerun_all_wells_corrected import CORRECTED_KW, NZ, WELLS, _stop_t, _total_t  # noqa: E402


def run_case(well_id: str, *, well_override: WellSpec | None = None,
             extra_kw: dict | None = None) -> dict:
    """跑单个敏感性 case，返回与基线 JSON 同构的指标行。"""
    loader = next(l for w, l in WELLS if w == well_id)
    t0 = time.perf_counter()
    well, fluids, schedule, _ = loader()
    if well_override is not None:
        well = well_override
    cr = CasingFlowSolver(enable_gravity=True).run(well, fluids, schedule)
    inlet = build_coupled_annulus_inlet_provider(
        cr, CasingFlowSolver(enable_gravity=True), fluids, split_cement_phases=True)
    tt = min(_total_t(schedule) + 1200.0, _stop_t(cr, fluids) + 600.0)
    kw = dict(total_t=tt, nz=NZ, enable_cfl_adaptive=True, **CORRECTED_KW)
    if extra_kw:
        kw.update(extra_kw)
    res = AnnulusD2DGASolver(**kw).run(well, fluids, inlet)
    fr = cast(dict, res.summary["最终结果"])
    return {
        "well": well_id,
        "eta_E": float(fr["全井段最终有效顶替效率"]),
        "eta_N": float(fr.get("窄四分位效率", float("nan"))),
        "mixing": float(fr["最终混浆指数"]),
        "channeling": float(fr["最终窜槽指数"]),
        "cement_occ": float(fr["最终水泥浆占据率"]),
        "elapsed_s": round(time.perf_counter() - t0, 1),
    }


def load_or_run(case_path: Path, well_id: str, *, compute: dict,
                describe: str) -> dict:
    """断点续跑：case JSON 已存在则复用，否则计算并落盘。compute 传 run_case 参数。"""
    if case_path.exists():
        row = json.loads(case_path.read_text(encoding="utf-8"))
        print(f"  [复用] {case_path.name} eta_E={row['eta_E']:.4f}", flush=True)
        return row
    print(f"  [计算] {describe} ...", flush=True)
    row = run_case(well_id, **compute)
    case_path.parent.mkdir(parents=True, exist_ok=True)
    case_path.write_text(json.dumps(row, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"  [落盘] {case_path.name} eta_E={row['eta_E']:.4f} ({row['elapsed_s']}s)", flush=True)
    return row


def shift_standoff(well: WellSpec, delta: float) -> WellSpec:
    """对 standoff 剖面每点整体平移 delta（clip 到 [0,1]），返回新 WellSpec。

    口径：原剖面（loader 名义/实测）逐点 +delta，保持深度点与形状不变；
    delta=0 即基线原值。
    """
    from dataclasses import replace
    from cemdisp.data.well_spec import DepthValuePoint

    if not well.standoff_profile:
        raise ValueError(f"{well.well_name}: WellSpec 无 standoff_profile，无法平移")
    shifted = tuple(
        DepthValuePoint(p.depth_md_m, float(np.clip(p.value + delta, 0.0, 1.0)))
        for p in well.standoff_profile
    )
    return replace(well, standoff_profile=shifted)
