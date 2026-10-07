"""归因取证（2026-10-07）：hu101 四角 (F,F) 作业在当前工作树【不含 4b 列化 kwarg】
上的 η 复现值——用于判定探针 A hard 锚失配的归因（本波改动 vs 上游 4a/4c 耦合链）。

不引用任何 Phase 4b 属性/kwarg ⇒ 可在 git stash 后的原始 HEAD 树上原样运行。
装配 = probe_4b_columnwise_20261007.run_job 的 A 配置逐行复制。
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

BRANCH_ROOT = Path(__file__).resolve().parents[2]
for _p in (str(BRANCH_ROOT), str(BRANCH_ROOT / "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import cemdisp  # noqa: E402
assert Path(cemdisp.__file__).resolve().is_relative_to(BRANCH_ROOT), \
    f"cemdisp 指向错误树: {cemdisp.__file__}"

from entrypoints.run_sensitivity_current_20260916 import (  # noqa: E402
    _identity, annulus_kwargs_from_opts, build_pressure_field,
    build_temperature_fields, casing_kwargs_from_opts, normalize_run_opts,
    scale_schedule,
)
from entrypoints.rerun_all_wells_corrected import CORRECTED_KW, _stop_t, _total_t  # noqa: E402
from cemdisp.models2d import AnnulusD2DGASolver  # noqa: E402
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider  # noqa: E402
from cemdisp.transport1d import CasingFlowSolver  # noqa: E402

ANCHOR_JSON = (BRANCH_ROOT / "results/_probe_屈服门四角_2026-10-07/jobs"
               / "hu101_FF_nz250_hydro_r1.0.json")


def main() -> None:
    opts = normalize_run_opts({
        "enable_temperature_rheology": True,
        "temperature_mode": "static",
        "enable_yield_gate": None,
        "pressure_mode": "hydrostatic",
        "pressure_caliber": "shoe",
        "include_yield_term": False,
        "enable_stream_yield_gate": None,
    })
    from cemdisp.data.loaders.hu101_loader import load_hu101_tailpipe
    well, fluids, schedule, _ = load_hu101_tailpipe()
    well2, fluids2, schedule2 = _identity(well), _identity(fluids), scale_schedule(schedule, 1.0)
    field_1d, field_2d, _note = build_temperature_fields("呼101", "static")
    p_field = build_pressure_field(well2, fluids2, schedule2, opts["pressure_mode"])
    casing = CasingFlowSolver(
        enable_gravity=True, mixing_contact_time=True,
        plug_face_zero_mixing=True, has_plug=True,
        pressure_field=p_field, **casing_kwargs_from_opts(opts))
    cr = casing.run(well2, fluids2, schedule2, temperature_field=field_1d)
    inlet = build_coupled_annulus_inlet_provider(cr, casing, fluids2, split_cement_phases=True)
    tt = min(_total_t(schedule2) + 1200.0, _stop_t(cr, fluids2))
    solver_kw = {**CORRECTED_KW, **annulus_kwargs_from_opts(opts), "pressure_field": p_field}
    solver = AnnulusD2DGASolver(total_t=tt, nz=250, enable_cfl_adaptive=True, **solver_kw)
    res = solver.run(well2, fluids2, inlet, schedule=schedule2, temperature_field=field_2d)
    final = res.summary["最终结果"]
    eta_e = float(final["全井段最终有效顶替效率"])
    eta_n = float(final["窄四分位效率"])
    anchor = json.loads(ANCHOR_JSON.read_text(encoding="utf-8"))
    out = {"git_head_attribution": True, "cemdisp": cemdisp.__file__,
           "η_E": eta_e, "η_N": eta_n,
           "锚_η_E": float(anchor["η_E"]), "锚_η_N": float(anchor["η_N"]),
           "Δη_E": abs(eta_e - float(anchor["η_E"])),
           "Δη_N": abs(eta_n - float(anchor["η_N"]))}
    print(json.dumps(out, ensure_ascii=False, indent=1), flush=True)
    p = Path(__file__).resolve().parent / "attrib_head_result.json"
    p.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
