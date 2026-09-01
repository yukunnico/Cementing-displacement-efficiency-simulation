"""F2 修复验证（2026-09-01）：逐井对比停止标志修复前后的 1D 时间轴。

修复内容：
- casing_flow.py 尾缘重力修正改按"后继流体(pusher) vs 尾浆"配对（同一物理界面同刻），
  并落实 max 约束（水泥结束时刻 >= 任一水泥浆前缘到达时刻）；
- 8 个 runner 与 rerun_all_wells_corrected.py 的停止口径改为优先取 cement_end_time_s，
  rerun 的 total_t 去掉 +600s 尾窗（停于尾浆全部入库）。

本脚本在修复后的代码上运行：
- old_stop：复现修复前 runner 口径（前缘扫描；fronts 计算未变，逐位复现修复前数值）；
- new_stop：修复后 cement_end_time_s（=尾浆全部进入环空时刻，≡ 现场碰压断面）。
"""

from __future__ import annotations

import cemdisp.data.loaders as L
from cemdisp.data.fluid_spec import FluidRole
from cemdisp.transport1d import CasingFlowSolver

WELLS = [
    ("hu101", L.load_hu101_tailpipe),
    ("hu102", L.load_hu102_tailpipe),
    ("hu103", L.load_hu103_tailpipe),
    ("hu1", L.load_hu1_tailpipe),
    ("hu2", L.load_hu2_tailpipe),
    ("ht1_001", L.load_ht1_001_tailpipe),
    ("ht1_003", L.load_ht1_003_tailpipe),
    ("ht1_004", L.load_ht1_004_tailpipe),
]


def old_scan_stop(cr, fluids) -> float:
    """修复前 runner 口径：末段水泥前缘之后首个非水泥前缘（fronts 计算未变，逐位复现）。"""
    roles = {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}
    by = {f.name: f for f in fluids}
    fs = sorted(cr.fronts, key=lambda f: f.time_s)
    last = None
    for f in fs:
        fl = by.get(f.fluid_name)
        if fl is not None and fl.role in roles:
            last = f.time_s
    if last is not None:
        for f in fs:
            fl = by.get(f.fluid_name)
            if fl is None or fl.role in roles:
                continue
            if f.time_s >= last - 1.0e-9:
                return float(f.time_s)
    return float(cr.cement_end_time_s)


def main() -> None:
    rows = []
    for name, loader in WELLS:
        well, fluids, schedule, _ = loader()
        cr = CasingFlowSolver(enable_gravity=True).run(well, fluids, schedule)
        sched_end = cr.pumping_end_time_s
        old_stop = old_scan_stop(cr, fluids)
        new_stop = cr.cement_end_time_s
        if new_stop is None:
            new_stop = sched_end
        rows.append((name, sched_end, old_stop, new_stop))
        print(
            f"{name:8s} 日程末={sched_end:8.0f}s  旧口径stop={old_stop:8.0f}s  "
            f"新口径stop(尾浆全入库)={new_stop:8.0f}s  Δ={new_stop - old_stop:+7.0f}s"
        )
    lines = ["well,sched_end_s,old_stop_s,new_stop_s,delta_s"] + [
        f"{n},{se:.1f},{o:.1f},{w:.1f},{w - o:.1f}" for n, se, o, w in rows
    ]
    with open("results/数据校验_2026-09-01/stop_fix_对照.csv", "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("saved -> results/数据校验_2026-09-01/stop_fix_对照.csv")


if __name__ == "__main__":
    main()
