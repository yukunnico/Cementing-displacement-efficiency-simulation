"""2026-09-02 阶段4：因子2守恒修正后，按 runner 生产口径(nz=250, 默认开关)重算 8 井。
输出每井 η_E/η_N/质量守恒/前缘/壁面冻结，并与 09-02 修正前 runner 基线对照。"""
from __future__ import annotations
import json, time, sys
from pathlib import Path
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from scripts._mass_balance_diag_20260902 import WELLS, build_case, integrate_injection
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver, _trapez2d

OUT = PROJECT_ROOT / "results" / "_质量平衡取证_2026-09-02"
# 09-02 修正前 runner 基线（nz=250，results/全井runner重跑_停止修复_2026-09-02/汇总.md）
OLD = {"hu101":0.3992,"hu102":0.3814,"hu103":0.1940,"hu1":0.5014,"hu2":0.5690,
       "ht1_001":0.5195,"ht1_003":0.6132,"ht1_004":0.6284}

def main(nz=250, names=None):
    names = names or list(WELLS)
    rows = []
    for name in names:
        well_spec, fluids, schedule, provider, stop, v_cem_design = build_case(WELLS[name])
        v_ann = AnnulusD2DGASolver()._physical_annular_volume(well_spec)
        inj = integrate_injection(provider, stop)
        t0 = time.time()
        solver = AnnulusD2DGASolver(total_t=stop, nz=nz, ny=40)  # 与 runner 完全一致的默认口径
        res = solver.run(well_spec, fluids, provider, schedule=schedule)
        g = res.geom
        v_dom = 2.0*_trapez2d(g["b"]*res.cement_field, g)
        fin = res.metrics.iloc[-1]
        sm = res.summary
        row = {
            "井": name, "修正前η_E": OLD.get(name),
            "修正后η_E": round(float(fin["effective_efficiency"]),4),
            "η_N窄边": round(float(sm.get("eta_narrow", np.nan)),4),
            "入环空水泥_m3": round(inj["cement"],2), "域内水泥_m3": round(v_dom,2),
            "守恒率": round(v_dom/max(inj["cement"],1e-9),4),
            "库存比": round(v_cem_design/v_ann,3),
            "宽边前缘m": round(float(fin["front_wide_m"]),1),
            "窄边前缘m": round(float(fin["front_narrow_m"]),1),
            "域长m": round(float(g["s"][-1]),1),
            "壁面冻结占比": round(float(fin["mean_wall_mud"]),4),
            "混浆指数": round(float(fin.get("mixing_index",np.nan)),4),
            "耗时s": round(time.time()-t0,1),
        }
        rows.append(row)
        print(json.dumps(row, ensure_ascii=False), flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT/f"阶段4_八井守恒修正后_nz{nz}.csv", index=False, encoding="utf-8-sig")
    (OUT/f"阶段4_八井守恒修正后_nz{nz}.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    print("\n", df.to_string(index=False))

if __name__ == "__main__":
    import sys
    nz = int(sys.argv[1]) if len(sys.argv)>1 else 250
    names = sys.argv[2].split(",") if len(sys.argv)>2 else None
    main(nz, names)
