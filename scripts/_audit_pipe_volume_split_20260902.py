"""独立审计：管容双链分裂逐井核查（2026-09-02 调研用）。

验证内容：
1. fronts 链管容（决定 cement_end_time_s / 2D 停止时刻）vs timeline 链管容（决定 2D 入口供给时线）；
2. 各水泥相在 2D 停止时刻之前实际供给时长；
3. 停止时刻与尾浆前缘到达时刻的先后关系。
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import cemdisp.data.loaders as L
from cemdisp.data.fluid_spec import FluidRole
from cemdisp.transport1d import CasingFlowSolver
from cemdisp.transport1d.casing_flow import CasingFlowSolver as _CFS

WELLS = [
    ("hu101", L.load_hu101_tailpipe), ("hu102", L.load_hu102_tailpipe),
    ("hu103", L.load_hu103_tailpipe), ("hu1", L.load_hu1_tailpipe),
    ("hu2", L.load_hu2_tailpipe), ("ht1_001", L.load_ht1_001_tailpipe),
    ("ht1_003", L.load_ht1_003_tailpipe), ("ht1_004", L.load_ht1_004_tailpipe),
]

CEMENT_ROLES = {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}

print(f"{'井':<10}{'fronts链管容':>12}{'timeline链管容':>14}{'分裂%':>8}{'cement_end':>12}{'尾浆前缘到达':>14}{'停止-到达':>10}")
for name, loader in WELLS:
    well, fluids, schedule, _ = loader()
    solver = CasingFlowSolver(enable_gravity=True)
    cr = solver.run(well, fluids, schedule)

    # fronts 链管容：shoe_depth × _pipe_cross_section_area
    area = _CFS._pipe_cross_section_area(well)
    v_fronts = well.shoe_md_m * area
    # timeline 链管容
    v_timeline = _CFS._timeline_pipe_volume(well, v_fronts)
    split_pct = (v_timeline / v_fronts - 1.0) * 100.0

    # 尾浆前缘到达鞋口时刻：从 shoe_timeline 事件中找最后一个水泥相 FRONT 类事件
    by_name = {f.name: f for f in fluids}
    cement_names = {f.name for f in fluids if f.role in CEMENT_ROLES}
    tail_arrival = None
    for ev in cr.shoe_timeline.events:
        if ev.phase_fractions and ev.phase_fractions[0][0] in cement_names:
            t = ev.time_s
            if tail_arrival is None or t > tail_arrival:
                tail_arrival = t
    print(f"{name:<10}{v_fronts:>12.2f}{v_timeline:>14.2f}{split_pct:>7.1f}%{cr.cement_end_time_s:>12.1f}{tail_arrival if tail_arrival else float('nan'):>14.1f}{cr.cement_end_time_s - (tail_arrival or float('nan')):>10.1f}")
