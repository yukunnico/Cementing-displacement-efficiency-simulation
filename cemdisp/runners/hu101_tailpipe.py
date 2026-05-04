"""
呼101尾管段顶替效率模型运行器

运行两种模式并导出中文命名结果：
1. 硬编码环空入口(sustained_tail)：替浆阶段环空入口按尾浆等效处理
2. 1D-2D耦合：套管内前沿追踪 → 鞋口出流 → 环空入口

输出目录：results/呼101尾管_初版模型/ 和 results/呼101尾管_1D2D耦合模型/
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
import json
from pathlib import Path
from typing import cast

import numpy as np

from cemdisp.data.loaders import build_hu101_annulus_inlet_provider, load_hu101_tailpipe
from cemdisp.data.pumping_schedule import PumpingSchedule
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import AnnulusInletState, build_coupled_annulus_inlet_provider
from cemdisp.reporting.animation import animate_cement_field
from cemdisp.reporting.contour_plots import (
    plot_annulus_snapshots,
    plot_depth_time_contour,
    plot_final_fields_contour,
)
from cemdisp.reporting.plots import (
    plot_depth_profiles,
    plot_efficiency_summary_bar,
    plot_risk_indices,
    plot_time_series,
)
from cemdisp.transport1d import CasingFlowSolver


_CEMDISP_ROOT = Path(__file__).resolve().parents[1]  # cemdisp/
PROJECT_ROOT = _CEMDISP_ROOT.parent                  # cement model/


def _schedule_total_time_s(schedule: PumpingSchedule) -> float:
    """按现场分段排量计算施工总时长。"""

    return sum(0.0 if step.rate_m3_min <= 0.0 else step.volume_m3 / step.rate_m3_min * 60.0 for step in schedule.steps)


def run_and_export(
    *,
    mode_title: str,
    output_dir: Path,
    inlet_provider: Callable[[float], AnnulusInletState],
) -> None:
    """运行呼101环空模型并导出一套中文命名结果。"""

    output_dir.mkdir(parents=True, exist_ok=True)
    well_spec, fluids, schedule, _ = load_hu101_tailpipe()
    # 呼101现场施工存在多段排量：1.2、1.5、1.0、0.55 m³/min，并非 Hu102 的单一平均排量。
    # 因此 total_t 不写死为固定小时数，而按现场施工程序逐段累加得到碰压前总时长，
    # 再增加20分钟窗口，覆盖停泵后早期重力分异与场数据导出快照。
    # 另外沿用 Hu101 legacy 模型口径：
    # 1) nz=500：legacy 收敛检查以 500 作为主计算网格；
    # 2) quality_penalty_scale=0.671：按 Hu101 现场 CBL 合格率 62.77% 做单井校准。
    total_t = _schedule_total_time_s(schedule) + 20.0 * 60.0
    solver = AnnulusD2DGASolver(
        total_t=total_t,
        enable_gravity=True,
        nz=500,
        quality_penalty_scale=0.671,
    )
    result = solver.run(well_spec, fluids, inlet_provider)

    metrics_path = output_dir / f"呼101尾管_{mode_title}_时间序列结果.csv"
    profiles_path = output_dir / f"呼101尾管_{mode_title}_深度剖面.csv"
    _ = result.metrics.to_csv(metrics_path, index=False, encoding="utf-8-sig")  # pyright: ignore[reportUnknownMemberType]
    _ = result.depth_profiles.to_csv(profiles_path, index=False, encoding="utf-8-sig")  # pyright: ignore[reportUnknownMemberType]

    summary_json_path = output_dir / f"呼101尾管_{mode_title}_结果摘要.json"
    _ = summary_json_path.write_text(json.dumps(result.summary, ensure_ascii=False, indent=2), encoding="utf-8")

    summary_md_path = output_dir / f"呼101尾管_{mode_title}_结果摘要.md"
    final_result = cast(Mapping[str, float], result.summary["最终结果"])
    _ = summary_md_path.write_text(
        "\n".join(
            [
                f"# 呼101尾管{mode_title}结果摘要",
                "",
                f"- 模拟对象：{result.summary['模拟对象']}",
                f"- 全井段最终有效顶替效率：{final_result['全井段最终有效顶替效率']:.4f}",
                f"- CBL评价井段模拟有效顶替效率：{final_result['CBL评价井段模拟有效顶替效率']:.4f}",
                f"- 目标层段模拟有效顶替效率：{final_result['目标层段模拟有效顶替效率']:.4f}",
                f"- 最终水泥浆占据率：{final_result['最终水泥浆占据率']:.4f}",
                f"- 最终质量响应效率：{final_result['最终质量响应效率']:.4f}",
                f"- 最终窜槽/混浆/失稳指数：{final_result['最终窜槽指数']:.4f} / {final_result['最终混浆指数']:.4f} / {final_result['最终失稳指数']:.4f}",
            ]
        ),
        encoding="utf-8",
    )

    _ = plot_time_series(result, output_dir=output_dir)
    _ = plot_depth_profiles(result, well_spec=well_spec, output_dir=output_dir)
    _ = plot_risk_indices(result, output_dir=output_dir)
    _ = plot_efficiency_summary_bar(result, output_dir=output_dir)
    _ = plot_depth_time_contour(result, output_dir=output_dir)
    _ = plot_annulus_snapshots(result, output_dir=output_dir)
    _ = plot_final_fields_contour(result, output_dir=output_dir)

    npz_path = output_dir / f"呼101尾管_{mode_title}_2D场数据.npz"
    _ = np.savez(
        npz_path,
        cement_snapshots=np.array(result.cement_snapshots),
        spacer_snapshots=np.array(result.spacer_snapshots),
        wall_snapshots=np.array(result.wall_snapshots),
        gel_strength_snapshots=np.array(result.gel_strength_snapshots),
        mud_cake_snapshots=np.array(result.mud_cake_snapshots),
        reynolds_snapshots=np.array(result.reynolds_snapshots),
        turbulent_viscosity_snapshots=np.array(result.turbulent_viscosity_snapshots),
        snapshot_times_s=np.array(result.snapshot_times_s),
        md=result.geom["md"],
        y=result.geom["y"],
        cement_final=result.cement_field,
        spacer_final=result.spacer_field,
        wall_final=result.wall_field,
        mud_cake_final=result.mud_cake_field,
    )

    _ = animate_cement_field(result, output_dir=output_dir, save_format="gif")

    print(f"\n=== {mode_title} ===")
    print(json.dumps(result.summary, ensure_ascii=False, indent=2))


def run_hu101_tailpipe_initial() -> None:
    """呼101尾管段顶替效率模型完整运行入口。"""

    well_spec, fluids, schedule, _ = load_hu101_tailpipe()

    # 硬编码环空入口模式：保留 52m³ 鞋口滞后，但不再把尾浆无限持续送入环空；
    # 尾浆越鞋后，后续轻泥浆/中置液/井浆按 mud invasion 真实进入环空，贴近现场过替过程。
    hardcoded_provider = build_hu101_annulus_inlet_provider(
        schedule,
        fluids,
        annulus_boundary_mode="field_order_realistic",
        split_cement_phases=True,
    )
    run_and_export(
        mode_title="初版模型",
        output_dir=PROJECT_ROOT / "results" / "呼101尾管_初版模型",
        inlet_provider=hardcoded_provider,
    )

    # 1D-2D耦合模式：由现场分段施工程序先经过套管内前沿追踪，再转成环空入口边界。
    casing_solver = CasingFlowSolver(enable_gravity=True)
    casing_result = casing_solver.run(well_spec, fluids, schedule)
    coupled_provider = build_coupled_annulus_inlet_provider(
        casing_result,
        casing_solver,
        fluids,
        split_cement_phases=True,
    )
    run_and_export(
        mode_title="1D2D耦合模型",
        output_dir=PROJECT_ROOT / "results" / "呼101尾管_1D2D耦合模型",
        inlet_provider=coupled_provider,
    )


if __name__ == "__main__":
    run_hu101_tailpipe_initial()
