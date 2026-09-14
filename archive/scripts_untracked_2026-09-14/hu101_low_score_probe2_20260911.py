"""呼101 第二次归因探针（spike，2026-09-11）：阈值锐度 + 结构开关复检。

承接 hu101_low_score_attribution_probe_20260911.py 的第一轮。第二轮回答两个问题：

Q1【阈值锐度】呼101 低分是"忠实反映输入卡在动员阈值下"，还是数值人造的刀锋？
   均匀 SO 细扫 {0.40,0.42,0.44,0.46,0.48,0.50,0.52}，
   并以 hu103(低 e) 同口径对照——若同一段 SO 在低 e 井上平滑、在高 e 井上断崖，
   则锐度来自偏心环空屈服门本身的物理分叉，而非 bug。

Q2【结构开关】除了第一轮的浮力/e_clip/弥散，呼101 还有哪些结构性开关有量级？
   - open_outlet=False：B1 报告显示呼101 稳态守恒比仅 0.603–0.645（其余井 1.1–1.24），
     即水泥沿窄边窜顶**流出评价域**。关闭出口是否显著抬升 η_E？
   - f_amp→1（B1 缺陷单元隔离）：复核 B1 报告中呼101 +3.66pp 的虚高。
   - enable_yield_gate=False：壁面冻结层是否就是刀锋的来源？

口径与 runners/hu101_tailpipe.py 逐位一致。一次性探针，不改生产默认值。
"""

from __future__ import annotations

import csv
import json
import time
from dataclasses import replace
from pathlib import Path
from typing import Any, cast

import numpy as np

import cemdisp.models2d.annulus_d2dga as _annulus
from cemdisp.data.loaders import load_hu101_tailpipe, load_hu103_tailpipe
from cemdisp.data.well_spec import DepthValuePoint
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.runners.hu101_tailpipe import annulus_stop_time_s
from cemdisp.transport1d import CasingFlowSolver

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PROJECT_ROOT / "results" / "呼101_低分归因探针2_2026-09-11"

NZ = 250
CBL_WINDOW = "CBL评价井段(单层套管可评价段)"
CASING_KW = dict(
    enable_gravity=True,
    mixing_contact_time=True,
    plug_face_zero_mixing=True,
    has_plug=True,
)


def _field_stats(cement: np.ndarray) -> dict[str, float]:
    ny = cement.shape[0]
    q = max(ny // 4, 1)
    bands = np.array_split(np.arange(cement.shape[1]), 4)
    return {
        "wide_q1": float(cement[:q, :].mean()),
        "narrow_q4": float(cement[-q:, :].mean()),
        "depth_bottom_band": float(cement[:, bands[0]].mean()),
        "depth_top_band": float(cement[:, bands[-1]].mean()),
    }


def run_case(*, label, well_spec, fluids, schedule, solver_kw=None, nz=NZ, cache=None) -> dict[str, Any]:
    solver_kw = dict(solver_kw or {})
    key = f"{well_spec.well_name}|{hash(tuple((f.name, f.density_kg_m3) for f in fluids))}"
    if cache is not None and key in cache:
        stop, cres, csol = cache[key]
    else:
        csol = CasingFlowSolver(**CASING_KW)  # type: ignore[arg-type]
        cres = csol.run(well_spec, fluids, schedule)
        stop = annulus_stop_time_s(casing_result=cres, fluids=fluids)
        if cache is not None:
            cache[key] = (stop, cres, csol)

    t0 = time.perf_counter()
    inlet = build_coupled_annulus_inlet_provider(cres, csol, fluids, split_cement_phases=True)
    result = AnnulusD2DGASolver(total_t=stop, nz=nz, **solver_kw).run(
        well_spec, fluids, inlet, schedule=schedule
    )
    elapsed = time.perf_counter() - t0

    final = cast(dict, result.summary["最终结果"])
    cbl = cast(dict, result.summary["评价窗效率"]).get(CBL_WINDOW, {})
    row: dict[str, Any] = {
        "case": label,
        "well": well_spec.well_name,
        "eta_E_well": float(final["全井段最终有效顶替效率"]),
        "eta_N_well": float(final["窄四分位效率"]),
        "eta_E_cbl": float(cbl.get("eta_E", float("nan"))),
        "channeling_index": float(final["最终窜槽指数"]),
        "elapsed_s": round(elapsed, 1),
    }
    row.update(_field_stats(np.asarray(result.cement_field, dtype=float)))
    print(f"  [{label:<30s}] ηE全井={row['eta_E_well']:.4f} ηE_CBL={row['eta_E_cbl']:.4f} "
          f"ηN={row['eta_N_well']:.4f} 窄={row['narrow_q4']:.3f} ({elapsed:.0f}s)", flush=True)
    return row


def _uniform(well, so: float):
    return (
        DepthValuePoint(float(well.top_md_m), so),
        DepthValuePoint(float(well.bottom_md_m), so),
    )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    cache: dict[str, Any] = {}

    print("=" * 90)
    print("  呼101 归因探针 2：动员阈值锐度 + 结构开关复检")
    print("=" * 90)

    w101, f101, s101, _ = load_hu101_tailpipe()
    w103, f103, s103, _ = load_hu103_tailpipe()

    # ---- 基线 ----
    print("\n[A] 基线")
    rows.append(run_case(label="A_基线", well_spec=w101, fluids=f101, schedule=s101, cache=cache))

    # ---- Q2 结构开关 ----
    print("\n[Q2] 结构开关")
    rows.append(run_case(label="H_open_outlet=False", well_spec=w101, fluids=f101, schedule=s101,
                         solver_kw={"open_outlet": False}, cache=cache))
    rows.append(run_case(label="J_yield_gate=False", well_spec=w101, fluids=f101, schedule=s101,
                         solver_kw={"enable_yield_gate": False}, cache=cache))

    # f_amp≡1（B1 缺陷单元隔离）：只把通量放大因子置 1，其余 D2DGA 分支照常
    _orig_amp = _annulus.d2dga_flux_amplification
    _annulus.d2dga_flux_amplification = lambda *a, **k: 1.0  # type: ignore[assignment]
    try:
        rows.append(run_case(label="I_B1off(f_amp=1)", well_spec=w101, fluids=f101, schedule=s101,
                             cache=cache))
    finally:
        _annulus.d2dga_flux_amplification = _orig_amp  # type: ignore[assignment]

    # ---- Q1 阈值锐度：呼101 细扫 ----
    print("\n[Q1] 均匀 SO 细扫（呼101，高偏心 e≈0.56）")
    for so in (0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52):
        rows.append(run_case(
            label=f"K101_SO={so:.2f}",
            well_spec=replace(w101, standoff_profile=_uniform(w101, so), standoff_measured=False),
            fluids=f101, schedule=s101, cache=cache))

    # ---- Q1 对照：呼103 低偏心同样细扫 ----
    print("\n[Q1] 均匀 SO 细扫（呼103，低偏心 e≈0.28，对照）")
    for so in (0.40, 0.44, 0.48):
        rows.append(run_case(
            label=f"K103_SO={so:.2f}",
            well_spec=replace(w103, standoff_profile=_uniform(w103, so), standoff_measured=False),
            fluids=f103, schedule=s103, cache=cache))

    base = rows[0]
    print("\n" + "=" * 90)
    print(f"  {'case':<30s}{'ηE全井':>10s}{'Δ':>10s}{'ηE_CBL':>10s}{'ηN':>9s}{'窄':>8s}")
    for r in rows:
        print(f"  {r['case']:<30s}{r['eta_E_well']:>10.4f}{r['eta_E_well']-base['eta_E_well']:>+10.4f}"
              f"{r['eta_E_cbl']:>10.4f}{r['eta_N_well']:>9.4f}{r['narrow_q4']:>8.3f}")

    fields = ["case", "well", "eta_E_well", "eta_N_well", "eta_E_cbl", "channeling_index",
              "wide_q1", "narrow_q4", "depth_bottom_band", "depth_top_band", "elapsed_s"]
    with (OUTPUT_DIR / "归因探针2对比表.csv").open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})
    (OUTPUT_DIR / "归因探针2摘要.json").write_text(
        json.dumps({"description": "呼101 阈值锐度 + 结构开关复检", "cases": rows},
                   ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n结果已落盘：{OUTPUT_DIR}")


if __name__ == "__main__":
    main()
