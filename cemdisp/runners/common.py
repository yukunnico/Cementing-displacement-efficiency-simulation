"""Runner 通用逻辑与新旧入口对比工具

本模块提取各井 runner 的公共逻辑，避免六口井 runner 中重复出现
相同的边界同步、对比导出和文件命名代码。
"""

from __future__ import annotations

import csv
from collections.abc import Callable
from pathlib import Path

from cemdisp.models2d.boundary_bridge import AnnulusInletState


def export_old_vs_new_inlet_comparison(
    *,
    output_dir: Path,
    old_provider: Callable[[float], AnnulusInletState],
    new_provider: Callable[[float], AnnulusInletState],
    sample_times_s: tuple[float, ...],
    well_name: str,
) -> Path:
    """导出新旧入口边界逐时对比 CSV。

    在每一采样时刻同时查询旧入口假设与新同步入口，
    输出相态与排量的并排结果，便于判断边界同步变更对最终效率的影响。
    """

    comparison_path = output_dir / f"{well_name}_入口边界对比.csv"
    with comparison_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("time_s", "old_phase", "new_phase", "old_rate_m3_s", "new_rate_m3_s"),
        )
        writer.writeheader()
        for time_s in sample_times_s:
            old_state = old_provider(time_s)
            new_state = new_provider(time_s)
            writer.writerow(
                {
                    "time_s": f"{time_s:.3f}",
                    "old_phase": old_state.phase_fractions[0][0] if old_state.phase_fractions else "",
                    "new_phase": new_state.phase_fractions[0][0] if new_state.phase_fractions else "",
                    "old_rate_m3_s": f"{old_state.flow_rate_m3_s:.6f}",
                    "new_rate_m3_s": f"{new_state.flow_rate_m3_s:.6f}",
                }
            )
    return comparison_path
