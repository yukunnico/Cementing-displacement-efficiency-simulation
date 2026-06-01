"""
呼1-003井（HT1-003）尾管段顶替效率模型运行器

严格按呼1-003井（HT1-003）现场主作业运行 1D-2D 耦合模式并导出中文命名结果：
地面开泵 → 套管内前沿追踪 → 鞋口出流 → 环空入口 → 环空二维顶替

输出目录：results/呼1-003_1D2D耦合模型/
输出文件：CSV(时间序列/深度剖面) + JSON(摘要) + Markdown(摘要) + PNG(静态图) + NPZ(2D场数据) + GIF(动画)
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
from cemdisp.data.loaders.ht1_003_loader import load_ht1_003_tailpipe
from cemdisp.data.pumping_schedule import PumpingSchedule
from cemdisp.data.validation_data import ValidationData
from cemdisp.diagnostics.quality_proxy import compute_cbl_quality_proxy
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
from cemdisp.validation.cbl_comparison import validate_against_cbl


# 项目根目录：从cemdisp/runners向上两级到达cement model根目录
_CEMDISP_ROOT = Path(__file__).resolve().parents[1]  # cemdisp/
PROJECT_ROOT = _CEMDISP_ROOT.parent                  # cement model/


def run_and_export(
    *,
    mode_title: str,
    output_dir: Path,
    inlet_provider: Callable[[float], AnnulusInletState],
    total_t_s: float,
    validation_data: ValidationData | None = None,
) -> None:
    """运行环空模型并导出一套中文命名结果。

    参数：
        mode_title: 模式标题（如"初版模型"、"1D2D耦合模型"），用于文件名和打印
        output_dir: 结果输出目录
        inlet_provider: 环空入口边界状态提供器
        total_t_s: 环空二维顶替总时长（秒）
        validation_data: 可选的现场验证资料（含 CBL 报告路径等），用于后验对比
    """

    output_dir.mkdir(parents=True, exist_ok=True)
    well_spec, fluids, schedule, _ = load_ht1_003_tailpipe()
    # HT1-003 严格现场模式下，这里的 total_t_s 由 1D 鞋口时序决定：
    # 当替浆液第一次到达鞋口时，代表整段水泥浆已全部进入环空，
    # 环空顶替计算到此结束，不再继续让替浆液入环空稀释既有水泥场。
    solver = AnnulusD2DGASolver(total_t=total_t_s, nz=500)
    result = solver.run(well_spec, fluids, inlet_provider)

    # 导出CSV
    metrics_path = output_dir / f"呼1-003_{mode_title}_时间序列结果.csv"
    profiles_path = output_dir / f"呼1-003_{mode_title}_深度剖面.csv"
    _ = result.metrics.to_csv(metrics_path, index=False, encoding="utf-8-sig")  # pyright: ignore[reportUnknownMemberType]
    _ = result.depth_profiles.to_csv(profiles_path, index=False, encoding="utf-8-sig")  # pyright: ignore[reportUnknownMemberType]

    fluid_provenance_summary = build_injected_fluid_provenance_summary(well_spec.well_name, schedule, fluids)
    summary_payload = dict(result.summary)
    summary_payload["注入流体现场符合性检查"] = fluid_provenance_summary

    # 计算 CBL 质量风险代理值（独立模块，不等同于水力效率）
    final_metrics = result.metrics.iloc[-1]
    cbl_proxy = compute_cbl_quality_proxy(
        displacement_efficiency=float(final_metrics["cbl_eval_interval_efficiency"]),
        channeling_index=float(final_metrics["channeling_index"]),
        mixing_index=float(final_metrics["mixing_index"]),
        instability_index=float(final_metrics["instability_index"]),
    )
    summary_payload["CBL质量风险代理预测"] = {
        "点估计": cbl_proxy.point_estimate,
        "置信区间": [cbl_proxy.lower_bound, cbl_proxy.upper_bound],
        "置信水平": cbl_proxy.confidence_level,
        "水力有效顶替效率": cbl_proxy.hydraulic_efficiency,
    }

    # 若提供了现场 CBL 实测数据，进行后验验证对比
    cbl_validation_md_lines: list[str] = []
    if validation_data is not None and validation_data.cbl_pass_rate is not None:
        cbl_val = validate_against_cbl(
            well_name=result.well_name,
            cbl_interval_top_m=well_spec.evaluation_windows[0].top_md_m,
            cbl_interval_bottom_m=well_spec.evaluation_windows[0].bottom_md_m,
            simulated_efficiency=float(final_metrics["cbl_eval_interval_efficiency"]),
            measured_pass_rate=validation_data.cbl_pass_rate,
        )
        summary_payload["CBL实测对比验证"] = {
            "现场CBL合格率": cbl_val.measured_pass_rate,
            "绝对偏差": cbl_val.absolute_delta,
            "相对偏差_%": cbl_val.relative_delta_pct,
            "趋势一致性": cbl_val.trend_consistent,
        }
        cbl_validation_md_lines = [
            "",
            "## CBL 实测对比验证",
            "",
            f"- 模型预测水力效率：{cbl_val.simulated_efficiency:.4f}",
            f"- 现场 CBL 实测合格率：{cbl_val.measured_pass_rate:.4f}",
            f"- 绝对偏差：{cbl_val.absolute_delta:.4f}",
            f"- 相对偏差：{cbl_val.relative_delta_pct:.2f}%",
            f"- 趋势一致性：{cbl_val.trend_consistent}",
            "",
            "> **注意**：CBL 合格率受泥饼残留、水泥收缩、温度、气体窜槽等多因素影响，",
            "> 通常低于水力学有效顶替效率。偏差在可接受范围内表明模型预测合理。",
        ]

    # 导出JSON摘要（已追加 CBL 代理和验证结果）
    summary_json_path = output_dir / f"呼1-003_{mode_title}_结果摘要.json"
    _ = summary_json_path.write_text(
        json.dumps(summary_payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    # 导出Markdown摘要
    summary_md_path = output_dir / f"呼1-003_{mode_title}_结果摘要.md"
    final_result = cast(Mapping[str, float], result.summary["最终结果"])
    _ = summary_md_path.write_text(
        "\n".join(
            [
                f"# 呼1-003{mode_title}结果摘要",
                "",
                f"- 模拟对象：{result.summary['模拟对象']}",
                f"- 全井段最终有效顶替效率：{final_result['全井段最终有效顶替效率']:.4f}",
                f"- CBL评价井段水力有效顶替效率：{final_result['CBL评价井段水力有效顶替效率']:.4f}",
                f"- 目标层段水力有效顶替效率：{final_result['目标层段水力有效顶替效率']:.4f}",
                f"- 最终水泥浆占据率：{final_result['最终水泥浆占据率']:.4f}",
                f"- 最终窜槽/混浆/失稳指数：{final_result['最终窜槽指数']:.4f} / {final_result['最终混浆指数']:.4f} / {final_result['最终失稳指数']:.4f}",
                "",
                *format_injected_fluid_provenance_markdown(fluid_provenance_summary),
                "",
                "## CBL 质量风险代理预测",
                "",
                f"- 水力学有效顶替效率：{cbl_proxy.hydraulic_efficiency:.4f}",
                f"- CBL 质量代理点估计：{cbl_proxy.point_estimate:.4f}",
                f"- {cbl_proxy.confidence_level}：[{cbl_proxy.lower_bound:.4f}, {cbl_proxy.upper_bound:.4f}]",
                "",
                "> **声明**：CBL 质量代理值基于水力学结果与工程经验因子联合预测，",
                "> 不等同于 CBL 实测真值，仅供风险筛查和方案对比参考。",
                *cbl_validation_md_lines,
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

    # 导出2D场数据NPZ（水泥/隔离液/泥饼/触变/湍流快照 + 时间点 + 网格坐标）
    # 这些数组直接来自求解器结果对象，确保导出数据与模型实际计算场一致。
    npz_path = output_dir / f"呼1-003_{mode_title}_2D场数据.npz"
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
    """返回 HT1-003 环空二维顶替应停止的地面累计时间。

    当替浆液（水泥浆之后的首个非水泥流体）到达鞋口时停止，
    此时水泥浆已完全进入环空，继续泵入会稀释既有水泥场。
    """
    cement_roles = {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}
    fluid_by_name = {f.name: f for f in fluids}
    found_cement = False
    for front in casing_result.fronts:
        fluid = fluid_by_name.get(front.fluid_name)
        if fluid is None:
            continue
        if fluid.role in cement_roles:
            found_cement = True
            continue
        if found_cement:
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
    timing_path = output_dir / "呼1-003_1D2D耦合模型_鞋口出流时序.csv"
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
        if fluid_name in ("领浆", "尾浆", "尾管水泥浆"):
            if cement_arrival_time_s is None:
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


def run_ht1_003_tailpipe_initial() -> None:
    """呼1-003井（HT1-003）尾管段现场实录 1D-2D 耦合模型运行入口。"""

    well_spec, fluids, schedule, validation_data = load_ht1_003_tailpipe()
    output_dir = PROJECT_ROOT / "results" / "呼1-003_1D2D耦合模型"

    # 严格现场模式只使用 1D-2D 耦合：套管内前沿追踪 → 鞋口出流 → 环空入口。
    # 套管内同样启用重力项，使鞋口边界能反映停泵后的密度分异趋势。
    casing_solver = CasingFlowSolver(enable_gravity=True)
    casing_result = casing_solver.run(well_spec, fluids, schedule)
    coupled_provider = build_coupled_annulus_inlet_provider(casing_result, casing_solver, fluids)
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
        validation_data=validation_data,
    )


if __name__ == "__main__":
    run_ht1_003_tailpipe_initial()
