
from __future__ import annotations
import json
from pathlib import Path
import pandas as pd
from cemdisp.data.loaders.hu101_loader import load_hu101_tailpipe
from cemdisp.transport1d.casing_flow import CasingFlowSolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.models2d.annulus_paper_d2dga import AnnulusPaperD2DGASolver
from cemdisp.reporting.plots import plot_all

PROJECT_ROOT = Path(__file__).resolve().parents[2]

def run(output_root: Path | None=None):
    output_root = output_root or PROJECT_ROOT / 'results' / '呼101_1D2D_论文D2DGA耦合模型'
    fig_root = PROJECT_ROOT / 'figures'
    output_root.mkdir(parents=True, exist_ok=True)
    fig_root.mkdir(parents=True, exist_ok=True)
    well, fluids, schedule, validation = load_hu101_tailpipe()
    casing = CasingFlowSolver()
    casing_result = casing.run(well, fluids, schedule)
    provider = build_coupled_annulus_inlet_provider(casing, well, fluids, schedule)
    # 停止时刻：尾浆全部越过鞋口，后续替浆刚到环空入口；用于评价最终水泥放置。
    solver = AnnulusPaperD2DGASolver(total_t=casing_result.cement_end_time_s, nz=420, ny=48, dt=6.0)
    result = solver.run(well, fluids, provider)
    result.metrics.to_csv(output_root/'呼101_1D2D耦合_时间序列.csv',index=False,encoding='utf-8-sig')
    result.depth_profiles.to_csv(output_root/'呼101_1D2D耦合_深度剖面.csv',index=False,encoding='utf-8-sig')
    result.segment_efficiency.to_csv(output_root/'呼101_1D2D耦合_分段效率.csv',index=False,encoding='utf-8-sig')
    (output_root/'呼101_1D2D耦合_结果摘要.json').write_text(json.dumps(result.summary,ensure_ascii=False,indent=2),encoding='utf-8')
    casing_rows=[front.__dict__ for front in casing_result.fronts]
    pd.DataFrame(casing_rows).to_csv(output_root/'呼101_套管段1D前沿追踪.csv',index=False,encoding='utf-8-sig')
    plot_all(result, fig_root)
    return result, casing_result

if __name__ == '__main__':
    res, cas = run()
    print(json.dumps(res.summary, ensure_ascii=False, indent=2))
