"""2026-09-06 三项修复（WASH/SPACER 选相 + e_clip 裁定 + 幂律缝隙律）8 井量化重跑。

对照基线：results/胶塞语义修复后终跑_2026-09-03/汇总.csv（三项修复前）。
口径：
- 全部井用生产口径（纯默认 nz250 ny40 + 三项修复默认开）；
- hu101 额外跑 measured_standoff='between_centralizers'/'at_centralizers' 实测剖面
  口径（e_clip 裁定自动放开到 0.90 + 实测剖面本身更低 standoff）。
输出：results/三项修复重跑_2026-09-06/
"""
from __future__ import annotations

import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib  # noqa: E402

matplotlib.use("Agg")

from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver, _trapez2d  # noqa: E402
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider  # noqa: E402
from cemdisp.transport1d import CasingFlowSolver  # noqa: E402
from scripts.lib.mass_balance_diag import WELLS, build_case, integrate_injection  # noqa: E402

OUT = PROJECT_ROOT / "results" / "三项修复重跑_2026-09-06"
OUT.mkdir(parents=True, exist_ok=True)
NZ = 250


def run_one(tag: str, name: str, loader, **loader_kw) -> dict:
    t0 = time.time()
    well_spec, fluids, schedule, _ = loader(**loader_kw)
    casing_solver = CasingFlowSolver(enable_gravity=True)
    casing_result = casing_solver.run(well_spec, fluids, schedule)
    provider = build_coupled_annulus_inlet_provider(
        casing_result, casing_solver, fluids, split_cement_phases=True)
    stop = float(casing_result.cement_end_time_s)
    inj = integrate_injection(provider, stop)

    solver = AnnulusD2DGASolver(total_t=stop, nz=NZ, ny=40)
    res = solver.run(well_spec, fluids, provider, schedule=schedule)
    fin = res.metrics.iloc[-1]
    g = res.geom
    v_dom = 2.0 * _trapez2d(g["b"] * res.cement_field, g)
    v_full = 2.0 * _trapez2d(g["b"], g)
    s_max = float(g["s"][-1])
    # 记录等效隔离液密度（选相修复生效证据）
    mud_f, lead_f, tail_f, sp_f, _ = solver._pick_fluids(fluids)
    row = {
        "井": tag,
        "η_E": round(float(fin["effective_efficiency"]), 4),
        "η_N": round(float(res.summary.get("eta_narrow", np.nan)), 4),
        "stop_s": round(stop, 1),
        "入环空水泥_m3": round(inj["cement"], 2),
        "域内水泥_m3": round(v_dom, 2),
        "域满体积_m3": round(v_full, 2),
        "守恒率": round(v_dom / max(inj["cement"], 1e-9), 4),
        "宽边前缘m": round(float(fin["front_wide_m"]), 1),
        "窄边前缘m": round(float(fin["front_narrow_m"]), 1),
        "窄边到位率": round(float(fin["front_narrow_m"]) / s_max, 4),
        "混浆指数": round(float(fin["mixing_index"]), 4),
        "壁面冻结占比": round(float(fin["mean_wall_mud"]), 4),
        "等效隔离液密度": (sp_f.density_kg_m3 if sp_f is not None else None),
        "等效隔离液名": (sp_f.name if sp_f is not None else None),
        "standoff实测标记": bool(getattr(well_spec, "standoff_measured", False)),
        "耗时s": round(time.time() - t0, 1),
    }
    print(json.dumps(row, ensure_ascii=False), flush=True)
    (OUT / f"{tag}_摘要.json").write_text(
        json.dumps(res.summary, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    res.metrics.to_csv(OUT / f"{tag}_时间序列.csv", index=False, encoding="utf-8-sig")
    res.depth_profiles.to_csv(OUT / f"{tag}_深度剖面.csv", index=False, encoding="utf-8-sig")
    return row


def main(names: list[str]) -> None:
    rows = []
    for name in names:
        try:
            if name == "hu101":
                rows.append(run_one("hu101_名义", name, WELLS[name]))
                rows.append(run_one("hu101_实测between", name, WELLS[name],
                                    measured_standoff="between_centralizers"))
                rows.append(run_one("hu101_实测at", name, WELLS[name],
                                    measured_standoff="at_centralizers"))
            else:
                rows.append(run_one(name, name, WELLS[name]))
        except Exception as exc:
            rows.append({"井": name, "状态": "error", "错误": f"{type(exc).__name__}: {exc}"})
            traceback.print_exc()
        pd.DataFrame(rows).to_csv(OUT / "汇总.csv", index=False, encoding="utf-8-sig")
    print("\n[done] 全部完成，汇总：")
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    ns = sys.argv[1].split(",") if len(sys.argv) > 1 else list(WELLS)
    main(ns)
