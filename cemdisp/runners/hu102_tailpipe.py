"""
呼102尾管段顶替效率模型运行器

严格按呼102现场主作业运行 1D-2D 耦合模式并导出中文命名结果：
地面开泵 → 套管内前沿追踪 → 鞋口出流 → 环空入口 → 环空二维顶替

输出目录：results/呼102尾管_1D2D耦合模型/
输出文件：CSV(时间序列/深度剖面) + JSON(摘要) + Markdown(摘要) + PNG(静态图) + NPZ(2D场数据) + GIF(动画)

⚠️ 口径声明：本 runner 为基线口径（默认参数）；论文/正式 8 井数字见 scripts/rerun_all_wells_corrected.py（修正后配置）。
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
import csv
import json
from pathlib import Path
from typing import cast

import numpy as np

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.fluid_provenance import build_injected_fluid_provenance_summary, format_injected_fluid_provenance_markdown
from cemdisp.data.loaders import load_hu102_tailpipe
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
from cemdisp.transport1d import CasingFlowSolver
from cemdisp.transport1d.casing_flow import CasingFlowResult


# 项目根目录：从cemdisp/runners向上两级到达cement model根目录
_CEMDISP_ROOT = Path(__file__).resolve().parents[1]  # cemdisp/
PROJECT_ROOT = _CEMDISP_ROOT.parent                  # cement model/


def run_and_export(
    *,
    mode_title: str,
    output_dir: Path,
    inlet_provider: Callable[[float], AnnulusInletState],
    total_t_s: float,
) -> None:
    """运行环空模型并导出一套中文命名结果。

    参数：
        mode_title: 模式标题（如"初版模型"、"1D2D耦合模型"），用于文件名和打印
        output_dir: 结果输出目录
        inlet_provider: 环空入口边界状态提供器
        total_t_s: 环空二维顶替总时长（秒）
    """

    output_dir.mkdir(parents=True, exist_ok=True)
    well_spec, fluids, schedule, _ = load_hu102_tailpipe()
    # Hu102 严格现场模式下，这里的 total_t_s 由 1D 鞋口时序决定：
    # 当替浆液第一次到达鞋口时，代表整段水泥浆已全部进入环空，
    # 环空顶替计算到此结束，不再继续让替浆液入环空稀释既有水泥场。
    solver = AnnulusD2DGASolver(total_t=total_t_s, nz=250)
    # 传入泵注程序：使末尾 Tier0 诊断聚合的 T0-6 停泵衰减诊断可用
    result = solver.run(well_spec, fluids, inlet_provider, schedule=schedule)

    # 导出CSV
    metrics_path = output_dir / f"呼102尾管_{mode_title}_时间序列结果.csv"
    profiles_path = output_dir / f"呼102尾管_{mode_title}_深度剖面.csv"
    _ = result.metrics.to_csv(metrics_path, index=False, encoding="utf-8-sig")  # pyright: ignore[reportUnknownMemberType]
    _ = result.depth_profiles.to_csv(profiles_path, index=False, encoding="utf-8-sig")  # pyright: ignore[reportUnknownMemberType]

    fluid_provenance_summary = build_injected_fluid_provenance_summary(well_spec.well_name, schedule, fluids)
    summary_payload = dict(result.summary)
    summary_payload["注入流体现场符合性检查"] = fluid_provenance_summary

    # 导出JSON摘要
    summary_json_path = output_dir / f"呼102尾管_{mode_title}_结果摘要.json"
    _ = summary_json_path.write_text(
        json.dumps(summary_payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    # 导出Markdown摘要
    summary_md_path = output_dir / f"呼102尾管_{mode_title}_结果摘要.md"
    final_result = cast(Mapping[str, float], result.summary["最终结果"])
    _ = summary_md_path.write_text(
        "\n".join(
            [
                f"# 呼102尾管{mode_title}结果摘要",
                "",
                f"- 模拟对象：{result.summary['模拟对象']}",
                f"- 全井段最终有效顶替效率：{final_result['全井段最终有效顶替效率']:.4f}",
                f"- 最终水泥浆占据率：{final_result['最终水泥浆占据率']:.4f}",
                f"- 最终窜槽/混浆/失稳指数：{final_result['最终窜槽指数']:.4f} / {final_result['最终混浆指数']:.4f} / {final_result['最终失稳指数']:.4f}",
                "",
                *format_injected_fluid_provenance_markdown(fluid_provenance_summary),
            ]
        ),
        encoding="utf-8",
    )

    # 导出静态图表（中文标签+中文文件名）
    _ = plot_time_series(result, output_dir=output_dir)
    _ = plot_depth_profiles(result, well_spec=well_spec, output_dir=output_dir)
    _ = plot_risk_indices(result, output_dir=output_dir)
    _ = plot_efficiency_summary_bar(result, output_dir=output_dir)

    # 导出参考项目风格图件（顶替效率时程、水泥体积分数剖面、宽窄边前沿推进等）
    _ = export_reference_figure_set(result, well_spec, output_dir=output_dir)

    # 导出云图（深度-时间等值线 + 多时刻截面快照 + 最终场三联图）
    _ = plot_depth_time_contour(result, output_dir=output_dir)
    _ = plot_annulus_snapshots(result, output_dir=output_dir)
    _ = plot_final_fields_contour(result, output_dir=output_dir)

    # 导出2D场数据NPZ（水泥/隔离液/壁面快照 + 时间点 + 网格坐标）
    npz_path = output_dir / f"呼102尾管_{mode_title}_2D场数据.npz"
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
    )

    # 导出水泥浓度场时间演化动画（GIF格式）
    _ = animate_cement_field(result, output_dir=output_dir, save_format="gif")

    # 打印摘要
    print(f"\n=== {mode_title} ===")
    print(json.dumps(summary_payload, ensure_ascii=False, indent=2))


def annulus_stop_time_s(
    *,
    casing_result: CasingFlowResult,
    fluids: tuple[FluidSpec, ...],
) -> float:
    """返回 Hu102 环空二维顶替应停止的地面累计时间。

    遍历鞋口前缘序列：找到第一种水泥浆（TAIL 等）之后的首个非水泥流体
    （替浆液）到达鞋口的时刻，此时水泥浆已全部进入环空，停止避免替浆液稀释水泥场。
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


def _export_casing_flow_timing(
    *,
    output_dir: Path,
    schedule: PumpingSchedule,
    casing_result: CasingFlowResult,
    casing_solver: CasingFlowSolver,
) -> None:
    """导出地面泵注与鞋口出流时序，明确环空顶替时间口径。

    该表用于解释 1D-2D 耦合边界：surface_time 从地面开泵算起，
    annulus_time 从水泥浆首次到达鞋口、进入环空后算起。
    """

    output_dir.mkdir(parents=True, exist_ok=True)
    timing_path = output_dir / "呼102尾管_1D2D耦合模型_鞋口出流时序.csv"
    cumulative_volume_m3 = 0.0
    elapsed_time_s = 0.0
    scheduled_windows: list[tuple[str, float, float, float, float, float]] = []
    for step in schedule.steps:
        duration_s = 0.0 if step.rate_m3_min <= 0.0 else step.volume_m3 / step.rate_m3_min * 60.0
        start_time_s = elapsed_time_s
        end_time_s = start_time_s + duration_s
        scheduled_windows.append(
            (
                step.fluid_name,
                step.rate_m3_min,
                start_time_s,
                end_time_s,
                cumulative_volume_m3,
                cumulative_volume_m3 + step.volume_m3,
            )
        )
        elapsed_time_s = end_time_s
        cumulative_volume_m3 += step.volume_m3

    pump_end_time_s = elapsed_time_s
    pipe_volume_m3 = casing_result.shoe_md_m * casing_result.pipe_cross_section_m2

    def _arrival_time_for_volume(target_volume_m3: float) -> float | None:
        """按施工累计体积反算某体积坐标到达鞋口的地面累计时间。"""

        for _, rate_m3_min, start_time_s, end_time_s, volume_start_m3, volume_end_m3 in scheduled_windows:
            if target_volume_m3 <= volume_end_m3 + 1.0e-12:
                if rate_m3_min <= 0.0:
                    return end_time_s
                volume_into_step_m3 = max(target_volume_m3 - volume_start_m3, 0.0)
                return start_time_s + volume_into_step_m3 / rate_m3_min * 60.0
        return None

    arrival_times: list[float] = []
    cement_arrival_time_s: float | None = None
    for fluid_name, _, _, _, volume_start_m3, _ in scheduled_windows:
        arrival_time_s = _arrival_time_for_volume(volume_start_m3 + pipe_volume_m3)
        if arrival_time_s is None:
            continue
        arrival_times.append(arrival_time_s)
        if fluid_name == "尾管水泥浆":
            cement_arrival_time_s = arrival_time_s

    annulus_start_time_s = cement_arrival_time_s if cement_arrival_time_s is not None else pump_end_time_s
    sample_times = sorted({0.0, pump_end_time_s, *arrival_times})

    with timing_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "地面累计时间_s",
                "地面累计时间_min",
                "环空顶替时间_s",
                "环空顶替时间_min",
                "地面施工阶段",
                "鞋口出流流体",
                "环空入口相",
                "排量_m3_min",
            ),
        )
        writer.writeheader()
        for time_s in sample_times:
            state = casing_solver.pipe_exit_state_at(casing_result, time_s)
            inlet_fluid = state.phase_fractions[0][0] if state.phase_fractions else "未知"
            annulus_time_s = max(time_s - annulus_start_time_s, 0.0)
            writer.writerow(
                {
                    "地面累计时间_s": f"{time_s:.3f}",
                    "地面累计时间_min": f"{time_s / 60.0:.3f}",
                    "环空顶替时间_s": f"{annulus_time_s:.3f}",
                    "环空顶替时间_min": f"{annulus_time_s / 60.0:.3f}",
                    "地面施工阶段": state.stage_name,
                    "鞋口出流流体": inlet_fluid,
                    "环空入口相": "水泥相" if inlet_fluid == "尾管水泥浆" else "泥浆相",
                    "排量_m3_min": f"{state.flow_rate_m3_s * 60.0:.6f}",
                }
            )


def run_hu102_tailpipe_initial() -> None:
    """呼102尾管段现场实录 1D-2D 耦合模型运行入口。"""

    well_spec, fluids, schedule, _ = load_hu102_tailpipe()
    output_dir = PROJECT_ROOT / "results" / "呼102尾管_1D2D耦合模型"

    # 严格现场模式只使用 1D-2D 耦合：套管内前沿追踪 → 鞋口出流 → 环空入口。
    # 套管内同样启用重力项，使鞋口边界能反映停泵后的密度分异趋势。
    casing_solver = CasingFlowSolver(enable_gravity=True)
    casing_result = casing_solver.run(well_spec, fluids, schedule)
    # split_cement_phases=True 与其余 7 井口径一致（领浆/中间浆并入前导水泥相，尾浆单独成相）。
    coupled_provider = build_coupled_annulus_inlet_provider(
        casing_result, casing_solver, fluids, split_cement_phases=True
    )
    annulus_stop_time_value_s = annulus_stop_time_s(casing_result=casing_result, fluids=fluids)
    _export_casing_flow_timing(
        output_dir=output_dir,
        schedule=schedule,
        casing_result=casing_result,
        casing_solver=casing_solver,
    )
    run_and_export(
        mode_title="1D2D耦合模型",
        output_dir=output_dir,
        inlet_provider=coupled_provider,
        total_t_s=annulus_stop_time_value_s,
    )


if __name__ == "__main__":
    run_hu102_tailpipe_initial()
