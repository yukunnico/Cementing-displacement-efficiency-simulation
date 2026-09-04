# -*- coding: utf-8 -*-
"""(A) 三口井新旧结果最终占位剖面对比；(B) ht1_004 居中度敏感性决定性实验。
只读仓库数据；输出仅打印到 stdout。行序约定：深度剖面行序=井底→井顶。"""
import dataclasses
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(r"D:\users\desktop\research\控压固井项目\cement model")
NEW = ROOT / "results" / "胶塞语义修复后终跑_2026-09-03"
OLD = {
    "hu101": ROOT / "results" / "呼101尾管_1D2D耦合模型" / "呼101尾管_1D2D耦合模型_深度剖面.csv",
    "ht1_003": ROOT / "results" / "呼1-003_1D2D耦合模型" / "呼1-003_1D2D耦合模型_深度剖面.csv",
    "ht1_004": ROOT / "results" / "呼1-004_1D2D耦合模型" / "呼1-004_1D2D耦合模型_深度剖面.csv",
}


def prof_stats(path, tag):
    df = pd.read_csv(path)
    dcol = next(c for c in df.columns if "深" in c)

    def col(*keys):
        return next((c for c in df.columns if all(k in c for k in keys)), None)

    d = df[dcol].to_numpy(float)
    print(f"  [{tag}] 行={len(df)} 深度{d.min():.0f}~{d.max():.0f}m（行序井底→井顶）")
    for name, c in (
        ("水泥", col("水泥", "浓度")),
        ("前置液", col("前置", "浓度")),
        ("钻井液", col("钻井液", "浓度")),
    ):
        if c is None:
            print(f"    {name}: 缺列")
            continue
        v = df[c].to_numpy(float)
        n = len(v)
        band = np.where(v >= 0.1)[0]
        band_s = f"{d[band].min():.0f}~{d[band].max():.0f}m" if band.size else "无"
        print(
            f"    {name}: 均{v.mean():.3f} 底1/3={v[: n // 3].mean():.2f} "
            f"中1/3={v[n // 3 : 2 * n // 3].mean():.2f} 顶1/3={v[2 * n // 3 :].mean():.2f} "
            f"≥0.1带:{band_s}"
        )
    c = col("水泥", "浓度")
    if c is not None:
        hit = np.where(df[c].to_numpy(float) >= 0.5)[0]
        if hit.size:
            print(f"    水泥顶界(最浅≥0.5): {d[hit].min():.0f}m")


print("########## (A) 三口井最终占位：新(09-03终跑) vs 旧(09-02午 runner) ##########")
for w in ("hu101", "ht1_003", "ht1_004"):
    print("=" * 24, w, "=" * 24)
    prof_stats(NEW / f"{w}_深度剖面.csv", "新·终跑")
    if OLD[w].exists():
        prof_stats(OLD[w], "旧·runner")
    sm = json.loads((NEW / f"{w}_摘要.json").read_text(encoding="utf-8"))
    for k, v in sm.get("评价窗效率", {}).items():
        print(f"    [新窗口] {k}: eta_E={v.get('eta_E'):.4f} eta_N={v.get('eta_N'):.4f}")

print()
print("########## (B) ht1_004 居中度敏感性（nz=120, ny=40, 其余口径同终跑） ##########")
from cemdisp.data.loaders import load_ht1_004_tailpipe
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver

well_spec, fluids, schedule, _ = load_ht1_004_tailpipe()
cs = CasingFlowSolver(enable_gravity=True)
cr = cs.run(well_spec, fluids, schedule)
provider = build_coupled_annulus_inlet_provider(cr, cs, fluids, split_cement_phases=True)
stop = float(cr.cement_end_time_s)
print(f"stop={stop:.1f}s  域顶={well_spec.top_md_m:.1f} 域底={well_spec.bottom_md_m:.1f}")

VARIANTS = [
    ("so=0.83 设计假设(基线)", 0.83, 0.55),
    ("so=0.65", 0.65, 0.55),
    ("so=0.50", 0.50, 0.55),
    ("so=0.40 (e被clip0.55)", 0.40, 0.55),
    ("so=0.40 (e_clip=0.90)", 0.40, 0.90),
]
for label, so, e_clip in VARIANTS:
    t0 = time.time()
    new_prof = tuple(
        dataclasses.replace(p, value=so) for p in well_spec.standoff_profile
    )
    ws = dataclasses.replace(well_spec, standoff_profile=new_prof)
    solver = AnnulusD2DGASolver(total_t=stop, nz=120, ny=40, e_clip_max=e_clip)
    try:
        res = solver.run(ws, fluids, provider, schedule=schedule)
    except Exception as exc:
        print(f"  [{label}] ERROR {type(exc).__name__}: {exc}")
        continue
    fin = res.metrics.iloc[-1]
    s_max = float(res.geom["s"][-1])
    cbl = next(
        (
            f"{k}: eta_E={v.get('eta_E'):.4f}"
            for k, v in res.summary.get("评价窗效率", {}).items()
            if "CBL评价" in k
        ),
        "无CBL窗",
    )
    print(
        f"  [{label}] eta_E={float(fin['effective_efficiency']):.4f} "
        f"eta_N={float(res.summary.get('eta_narrow', float('nan'))):.4f} "
        f"宽/窄前缘={float(fin['front_wide_m']):.0f}/{float(fin['front_narrow_m']):.0f}m "
        f"(域长{s_max:.0f}) {cbl} 耗时{time.time() - t0:.0f}s",
        flush=True,
    )
