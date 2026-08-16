"""
呼探1-002井（HT1-002）尾管段顶替效率模型运行器

严格按呼探1-002现场主作业运行 1D-2D 耦合模式并导出中文命名结果：
地面开泵 → 套管内前沿追踪 → 鞋口出流 → 环空入口 → 环空二维顶替

输出目录：results/呼探1-002尾管_1D2D耦合模型/
输出文件：CSV(时间序列/深度剖面) + JSON(摘要) + Markdown(摘要) + PNG(静态图) + NPZ(2D场数据) + GIF(动画)
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
import csv
import json
from pathlib import Path
from typing import Protocol, cast

import numpy as np

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.fluid_provenance import build_injected_fluid_provenance_summary, format_injected_fluid_provenance_markdown
from cemdisp.data.loaders.hu2_loader import load_hu2_tailpipe
from cemdisp.data.pumping_schedule import PumpingSchedule
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import AnnulusInletState, build_coupled_annulus_inlet_provider
from cemdisp.reporting.animation import animate_cement_field
from cemdisp.reporting.reference_figures import export_reference_figure_set
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
from cemdisp.transport1d.casing_flow import CasingFlowResult


# 项目根目录：runner 位于 cemdisp/runners/ 下，向上两级到达 cement model 根目录。
PROJECT_ROOT = Path(__file__).resolve().parents[2]


class _CsvWritable(Protocol):
    """表示具备 to_csv 导出能力的表格对象。"""

    def to_csv(self, path_or_buf: Path, *, index: bool, encoding: str) -> None:
        """将表格写入 CSV 文件。"""


def _export_table_csv(table: _CsvWritable, path: Path) -> None:
    """以 UTF-8 BOM 编码导出表格，便于 Windows Excel 正确识别中文。"""

    table.to_csv(path, index=False, encoding="utf-8-sig")


def run_and_export(
    *,
    mode_title: str,
    output_dir: Path,
    inlet_provider: Callable[[float], AnnulusInletState],
    total_t_s: float,
) -> None:
    """运行呼探1-002环空二维模型并导出一套中文命名结果。

    参数：
        mode_title: 模式标题（例如“1D2D耦合模型”），用于中文文件名和摘要标题。
        output_dir: 结果输出目录，默认入口使用 results/呼探1-002尾管_1D2D耦合模型/。
        inlet_provider: 由 1D 鞋口出流桥接得到的环空入口边界状态提供器。
        total_t_s: 环空二维顶替计算总时长（秒），由替浆液到达鞋口时刻确定。
    """

    # 创建输出目录，并从呼2加载器读取井筒、流体和评价窗口等标准输入。
    output_dir.mkdir(parents=True, exist_ok=True)
    well_spec, fluids, schedule, _ = load_hu2_tailpipe()

    # 呼探1-002严格现场模式下，二维求解时长由 1D 套管内前沿追踪确定：
    # 当替浆液第一次到达鞋口时，水泥浆柱已全部进入环空，随后不再让替浆液继续稀释水泥场。
    solver = AnnulusD2DGASolver(total_t=total_t_s, nz=250)
    result = solver.run(well_spec, fluids, inlet_provider)

    # 导出时间序列与深度剖面 CSV；列名由求解器保持中文口径。
    metrics_path = output_dir / f"呼探1-002尾管_{mode_title}_时间序列结果.csv"
    profiles_path = output_dir / f"呼探1-002尾管_{mode_title}_深度剖面.csv"
    _export_table_csv(cast(_CsvWritable, result.metrics), metrics_path)
    _export_table_csv(cast(_CsvWritable, result.depth_profiles), profiles_path)

    fluid_provenance_summary = build_injected_fluid_provenance_summary(well_spec.well_name, schedule, fluids)
    summary_payload = dict(result.summary)
    summary_payload["注入流体现场符合性检查"] = fluid_provenance_summary

    # 导出 JSON 摘要，保留机器可读的完整结果字典。
    summary_json_path = output_dir / f"呼探1-002尾管_{mode_title}_结果摘要.json"
    _ = summary_json_path.write_text(
        json.dumps(summary_payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    # 导出 Markdown 摘要，字段名严格使用中文，避免混淆顶替效率与质量响应效率。
    summary_md_path = output_dir / f"呼探1-002尾管_{mode_title}_结果摘要.md"
    final_result = cast(Mapping[str, float], result.summary["最终结果"])
    _ = summary_md_path.write_text(
        "\n".join(
            [
                f"# 呼探1-002尾管{mode_title}结果摘要",
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

    # 导出静态图表；reporting 模块统一处理中文标题、坐标轴与图例。
    _ = plot_time_series(result, output_dir=output_dir)
    _ = plot_depth_profiles(result, well_spec=well_spec, output_dir=output_dir)
    _ = plot_risk_indices(result, output_dir=output_dir)
    _ = plot_efficiency_summary_bar(result, output_dir=output_dir)

    # 导出参考项目风格图件（顶替效率时程、水泥体积分数剖面、宽窄边前沿推进等）
    _ = export_reference_figure_set(result, well_spec, output_dir=output_dir)

    # 导出二维云图：深度-时间等值线、多时刻截面快照与最终场图。
    _ = plot_depth_time_contour(result, output_dir=output_dir)
    _ = plot_annulus_snapshots(result, output_dir=output_dir)
    _ = plot_final_fields_contour(result, output_dir=output_dir)

    # 导出 2D 场数据 NPZ。呼探1-002为多浆柱体系，必须包含 lead/tail 快照和最终场。
    npz_path = output_dir / f"呼探1-002尾管_{mode_title}_2D场数据.npz"
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

    # 导出水泥浓度场时间演化动画（GIF），用于快速检查环空推进形态。
    _ = animate_cement_field(result, output_dir=output_dir, save_format="gif")

    # 控制台打印完整摘要，便于批处理运行后直接查看关键指标。
    print(f"\n=== {mode_title} ===")
    print(json.dumps(summary_payload, ensure_ascii=False, indent=2))


def annulus_stop_time_s(
    *,
    casing_result: CasingFlowResult,
    fluids: tuple[FluidSpec, ...],
) -> float:
    """确定呼探1-002环空二维顶替应停止的地面累计时间。

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


def _is_cement_slurry_name(fluid_name: str) -> bool:
    """根据中文流体名称识别是否属于水泥浆柱。"""

    # loader 中的角色才是权威数据；本函数仅用于鞋口时序表的中文相名展示。
    cement_keywords = ("领浆", "中间浆", "尾浆", "水泥浆")
    excluded_keywords = ("替浆", "压塞", "平衡", "冲洗", "隔离")
    return any(keyword in fluid_name for keyword in cement_keywords) and not any(
        keyword in fluid_name for keyword in excluded_keywords
    )


def _annulus_phase_label(fluid_name: str) -> str:
    """将鞋口出流流体名称转换为时序表中的环空入口相中文标签。"""

    # 水泥浆、前置液/隔离液与泥浆分开标注，便于人工核查 1D-2D 边界切换。
    if _is_cement_slurry_name(fluid_name):
        return "水泥相"
    if any(keyword in fluid_name for keyword in ("隔离", "冲洗", "平衡")):
        return "前置液/隔离液相"
    return "泥浆相"


def _export_casing_flow_timing(
    *,
    output_dir: Path,
    schedule: PumpingSchedule,
    casing_result: CasingFlowResult,
    casing_solver: CasingFlowSolver,
) -> None:
    """导出地面泵注与鞋口出流时序表。

    该表用于说明 1D-2D 耦合时间口径：地面累计时间从开泵起算，
    环空顶替时间从第一种水泥浆到达鞋口、进入环空后起算。
    """

    # 创建输出目录并固定呼探1-002中文文件名，保持与其他结果文件的命名口径一致。
    output_dir.mkdir(parents=True, exist_ok=True)
    timing_path = output_dir / "呼探1-002尾管_1D2D耦合模型_鞋口出流时序.csv"

    # 将现场施工步骤展开为时间窗，后续按累计体积反算各流体前缘到鞋口的时刻。
    cumulative_volume_m3 = 0.0
    elapsed_time_s = 0.0
    scheduled_windows: list[tuple[str, float, float, float, float, float]] = []
    for step in schedule.steps:
        duration_s = 0.0 if step.rate_m3_min <= 0.0 else step.volume_m3 / step.rate_m3_min * 60.0
        start_time_s = step.start_time_s if step.start_time_s is not None else elapsed_time_s
        end_time_s = step.end_time_s if step.end_time_s is not None else start_time_s + duration_s
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

    # 计算施工结束时间和管内容积，用于把地面累计体积转换成鞋口到达时间。
    pump_end_time_s = elapsed_time_s
    pipe_volume_m3 = casing_result.shoe_md_m * casing_result.pipe_cross_section_m2

    def _arrival_time_for_volume(target_volume_m3: float) -> float | None:
        """按施工累计体积反算指定体积坐标到达鞋口的地面累计时间。"""

        # 顺序扫描施工时间窗，找到目标累计体积落入的施工阶段并线性换算到时间。
        for _, rate_m3_min, start_time_s, end_time_s, volume_start_m3, volume_end_m3 in scheduled_windows:
            if target_volume_m3 <= volume_end_m3 + 1.0e-12:
                if rate_m3_min <= 0.0:
                    return end_time_s
                volume_into_step_m3 = max(target_volume_m3 - volume_start_m3, 0.0)
                return start_time_s + volume_into_step_m3 / rate_m3_min * 60.0
        return None

    # 采样点包含开泵、停泵和所有能到达鞋口的流体前缘，便于审阅边界切换。
    arrival_times: list[float] = []
    cement_arrival_time_s: float | None = None
    for fluid_name, _, _, _, volume_start_m3, _ in scheduled_windows:
        arrival_time_s = _arrival_time_for_volume(volume_start_m3 + pipe_volume_m3)
        if arrival_time_s is None:
            continue
        arrival_times.append(arrival_time_s)
        if cement_arrival_time_s is None and _is_cement_slurry_name(fluid_name):
            cement_arrival_time_s = arrival_time_s

    # 环空时间零点取第一种水泥浆到达鞋口时刻；若无法识别，则退回施工结束时刻。
    annulus_start_time_s = cement_arrival_time_s if cement_arrival_time_s is not None else pump_end_time_s
    sample_times = sorted({0.0, pump_end_time_s, *arrival_times})

    # 使用 UTF-8 BOM 写出 CSV，确保中文列名和中文流体名在 Excel 中正确显示。
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
                    "环空入口相": _annulus_phase_label(inlet_fluid),
                    "排量_m3_min": f"{state.flow_rate_m3_s * 60.0:.6f}",
                }
            )


def run_hu2_tailpipe_initial() -> None:
    """呼探1-002尾管段现场实录 1D-2D 耦合模型运行入口。"""

    # 加载呼2标准数据结构：井筒、流体、施工程序和验证资料由 hu2_loader 统一提供。
    well_spec, fluids, schedule, _ = load_hu2_tailpipe()
    output_dir = PROJECT_ROOT / "results" / "呼探1-002尾管_1D2D耦合模型"

    # 严格现场耦合流程：先在套管内做 1D 前沿追踪，再把鞋口出流桥接到环空入口。
    # 启用重力项，使停泵和密度差对鞋口出流时序的影响能被保留到边界条件中。
    casing_solver = CasingFlowSolver(enable_gravity=True)
    casing_result = casing_solver.run(well_spec, fluids, schedule)
    coupled_provider = build_coupled_annulus_inlet_provider(
        casing_result,
        casing_solver,
        fluids,
        split_cement_phases=True,
    )
    annulus_stop_time_value_s = annulus_stop_time_s(casing_result=casing_result, fluids=fluids)

    # 输出地面施工与鞋口出流的对照表，便于追溯 1D-2D 边界切换口径。
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
    run_hu2_tailpipe_initial()
