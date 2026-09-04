# -*- coding: utf-8 -*-
"""复现 hu1/hu2 摘要中 eta_E=0 但 eta_N≈1 的评价窗（只读调试）。"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\users\desktop\research\控压固井项目\cement model")
sys.path.insert(0, str(ROOT))
RES = ROOT / "results" / "胶塞语义修复后终跑_2026-09-03"

from cemdisp.data.loaders.hu1_loader import load_hu1_tailpipe
from cemdisp.data.loaders.hu2_loader import load_hu2_tailpipe
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver

CASES = {
    "hu1": (load_hu1_tailpipe, ["CBL质量段(差)", "CBL质量段(不合格40m)"]),
    "hu2": (load_hu2_tailpipe, ["高压水层(地层目标)"]),
}

for wname, (loader, targets) in CASES.items():
    well_spec, fluids, schedule = loader()[:3]
    solver = AnnulusD2DGASolver(total_t=12000.0, nz=250, ny=40)
    geom = solver._build_geom(well_spec, None)
    npz = np.load(RES / f"{wname}_2D场数据.npz")
    cement = np.clip(npz["cement_final"], 0.0, 1.0)
    md = geom["md"]
    print(f"\n===== {wname} ===== 域 md {md.min():.1f}–{md.max():.1f}")
    for win in well_spec.evaluation_windows:
        if not any(t in win.name for t in targets):
            continue
        mask = (md >= win.top_md_m) & (md <= win.bottom_md_m)
        k = int(mask.sum())
        b_win = geom["b"][:, mask]
        c_win = cement[:, mask]
        n_q = max(1, 40 // 4)
        cq = c_win[-n_q:, :]
        bq = geom["b"][-n_q:, mask]
        print(f"窗[{win.name}] type={win.window_type} md {win.top_md_m:.1f}–{win.bottom_md_m:.1f} "
              f"命中格点 k={k} 网格步长≈{abs(np.diff(md[:2])[0]) if len(md) > 1 else 0:.1f}m")
        if k:
            print(f"  c_win: min={c_win.min():.6g} max={c_win.max():.6g} mean={c_win.mean():.6g}")
            print(f"  窄1/4 c: min={cq.min():.6g} max={cq.max():.6g} mean={cq.mean():.6g}")
            full = float(np.trapezoid(np.trapezoid(b_win * c_win, x=geom['s'][mask], axis=1),
                                      x=geom['y'], axis=0))
            den = float(np.trapezoid(np.trapezoid(b_win, x=geom['s'][mask], axis=1),
                                     x=geom['y'], axis=0))
            num_q = float(np.trapezoid(np.trapezoid(bq * cq, x=geom['s'][mask], axis=1),
                                       x=geom['y'][-n_q:], axis=0))
            den_q = float(np.trapezoid(np.trapezoid(bq, x=geom['s'][mask], axis=1),
                                       x=geom['y'][-n_q:], axis=0))
            print(f"  复算 eta_E={full / max(den, 1e-12):.10g}  eta_N={num_q / max(den_q, 1e-12):.10g}")
        else:
            print("  mask 为空 → 窗在域外，摘要中不应出现该窗")
    # 对照：摘要里的值
    sm = json.loads((RES / f"{wname}_摘要.json").read_text(encoding="utf-8"))
    for k2, v in sm["评价窗效率"].items():
        if any(t in k2 for t in targets):
            print(f"  摘要[{k2}]: eta_E={v['eta_E']} eta_N={v['eta_N']}")
