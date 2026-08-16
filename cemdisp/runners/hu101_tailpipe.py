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

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.fluid_provenance import build_injected_fluid_provenance_summary, format_injected_fluid_provenance_markdown
from cemdisp.data.loaders import load_hu101_tailpipe
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
from cemdisp.reporting.reference_figures import export_reference_figure_set
from cemdisp.transport1d import CasingFlowResult, CasingFlowSolver


_CEMDISP_ROOT = Path(__file__).resolve().parents[1]  # cemdisp/
PROJECT_ROOT = _CEMDISP_ROOT.parent                  # cement model/


def _schedule_total_time_s(schedule: PumpingSchedule) -> float:
    """按现场分段排量计算施工总时长。"""

    return sum(0.0 if step.rate_m3_min <= 0.0 else step.volume_m3 / step.rate_m3_min * 60.0 for step in schedule.steps)


def annulus_stop_time_s(
    *,
    casing_result: CasingFlowResult,
    fluids: tuple[FluidSpec, ...],
) -> float:
    """返回 Hu101 环空二维顶替应停止的地面累计时间。

    遍历鞋口前缘序列：找到水泥浆之后的首个非水泥流体到达鞋口的时刻，
    此时水泥浆已全部进入环空，停止避免替浆液稀释水泥场。
    """

    cement_roles = {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}
    fluid_by_name = {f.name: f for f in fluids}
    # 按到达时间升序排序后再扫描，对个别井的乱序 fronts 自防御：
    # 不可压缩管流下 front 本应按泵序单调，但历史上当某前缘在泵注结束前
    # 未到达鞋口时会回退为该步骤结束时间而插到更早前缘之间。
    # 排序后取"末段水泥之后的首个非水泥"到达时刻，不提前截断环空顶替。
    sorted_fronts = sorted(casing_result.fronts, key=lambda f: f.time_s)
    last_cement_time_s: float | None = None
    for front in sorted_fronts:
        fluid = fluid_by_name.get(front.fluid_name)
        if fluid is not None and fluid.role in cement_roles:
            last_cement_time_s = front.time_s
    if last_cement_time_s is not None:
        for front in sorted_fronts:
            fluid = fluid_by_name.get(front.fluid_name)
            if fluid is None or fluid.role in cement_roles:
                continue
            if front.time_s >= last_cement_time_s - 1.0e-9:
                return float(front.time_s)
    return float(casing_result.cement_end_time_s)


def run_and_export(
    *,
    mode_title: str,
    output_dir: Path,
    inlet_provider: Callable[[float], AnnulusInletState],
    total_t_s: float | None = None,
) -> None:
    """运行呼101环空模型并导出一套中文命名结果。"""

    output_dir.mkdir(parents=True, exist_ok=True)
    well_spec, fluids, schedule, _ = load_hu101_tailpipe()
    # 呼101现场施工存在多段排量：1.2、1.5、1.0、0.55 m³/min，并非 Hu102 的单一平均排量。
    # 因此 total_t 不写死为固定小时数，而按现场施工程序逐段累加得到碰压前总时长，
    # 再增加20分钟窗口，覆盖停泵后早期重力分异与场数据导出快照。
    # 另外沿用 Hu101 legacy 模型的主计算网格口径：nz=250。
    total_t = total_t_s if total_t_s is not None else _schedule_total_time_s(schedule) + 20.0 * 60.0
    solver = AnnulusD2DGASolver(
        total_t=total_t,
        nz=250,
    )
    result = solver.run(well_spec, fluids, inlet_provider)

    metrics_path = output_dir / f"呼101尾管_{mode_title}_时间序列结果.csv"
    profiles_path = output_dir / f"呼101尾管_{mode_title}_深度剖面.csv"
    _ = result.metrics.to_csv(metrics_path, index=False, encoding="utf-8-sig")  # pyright: ignore[reportUnknownMemberType]
    _ = result.depth_profiles.to_csv(profiles_path, index=False, encoding="utf-8-sig")  # pyright: ignore[reportUnknownMemberType]

    summary_payload = dict(result.summary)

    summary_json_path = output_dir / f"呼101尾管_{mode_title}_结果摘要.json"
    _ = summary_json_path.write_text(json.dumps(summary_payload, ensure_ascii=False, indent=2), encoding="utf-8")

    summary_md_path = output_dir / f"呼101尾管_{mode_title}_结果摘要.md"
    final_result = cast(Mapping[str, float], result.summary["最终结果"])
    _ = summary_md_path.write_text(
        "\n".join(
            [
                f"# 呼101尾管{mode_title}结果摘要",
                "",
                f"- 模拟对象：{result.summary['模拟对象']}",
                f"- 全井段最终有效顶替效率：{final_result['全井段最终有效顶替效率']:.4f}",
                f"- 最终水泥浆占据率：{final_result['最终水泥浆占据率']:.4f}",
                f"- 最终窜槽/混浆/失稳指数：{final_result['最终窜槽指数']:.4f} / {final_result['最终混浆指数']:.4f} / {final_result['最终失稳指数']:.4f}",
            ]
        ),
        encoding="utf-8",
    )

    _ = plot_time_series(result, output_dir=output_dir)
    _ = plot_depth_profiles(result, well_spec=well_spec, output_dir=output_dir)
    _ = plot_risk_indices(result, output_dir=output_dir)
    _ = plot_efficiency_summary_bar(result, output_dir=output_dir)
    _ = export_reference_figure_set(result, well_spec, output_dir=output_dir)
    _ = plot_depth_time_contour(result, output_dir=output_dir)
    _ = plot_annulus_snapshots(result, output_dir=output_dir)
    _ = plot_final_fields_contour(result, output_dir=output_dir)

    npz_path = output_dir / f"呼101尾管_{mode_title}_2D场数据.npz"
    _ = np.savez(
        npz_path,
        cement_snapshots=np.array(result.cement_snapshots),
        spacer_snapshots=np.array(result.spacer_snapshots),
        wall_snapshots=np.array(result.wall_snapshots),
        snapshot_times_s=np.array(result.snapshot_times_s),
        md=result.geom["md"],
        y=result.geom["y"],
        cement_final=result.cement_field,
        spacer_final=result.spacer_field,
        wall_final=result.wall_field,
        lead_snapshots=np.array(result.lead_snapshots),
        tail_snapshots=np.array(result.tail_snapshots),
        lead_final=result.lead_field,
        tail_final=result.tail_field,
    )

    _ = animate_cement_field(result, output_dir=output_dir, save_format="gif")

    print(f"\n=== {mode_title} ===")
    print(json.dumps(summary_payload, ensure_ascii=False, indent=2))


def run_hu101_tailpipe_initial() -> None:
    """呼101尾管段顶替效率模型完整运行入口。"""

    well_spec, fluids, schedule, _ = load_hu101_tailpipe()

    # 1D-2D耦合模式：由现场分段施工程序先经过套管内前沿追踪，再转成环空入口边界。
    casing_solver = CasingFlowSolver(enable_gravity=True)
    casing_result = casing_solver.run(well_spec, fluids, schedule)
    annulus_stop_time_value_s = annulus_stop_time_s(casing_result=casing_result, fluids=fluids)
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
        total_t_s=annulus_stop_time_value_s,
    )


if __name__ == "__main__":
    run_hu101_tailpipe_initial()
