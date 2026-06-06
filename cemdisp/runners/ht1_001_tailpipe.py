"""
呼探1-001井尾管段顶替效率模型运行器

严格按呼探1-001井现场主作业运行 1D-2D 耦合模式并导出中文命名结果：
地面开泵 → 套管内前沿追踪 → 鞋口出流 → 环空入口 → 环空二维顶替

注意：呼探1-001井与呼探1井是两口不同的井，本模块专用于呼探1-001井。

输出目录：results/呼探1-001尾管_1D2D耦合模型/
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
from cemdisp.data.loaders.ht1_001_loader import load_ht1_001_tailpipe
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
from cemdisp.transport1d.casing_flow import CasingFlowResult


# 项目根目录：runner 位于 cemdisp/runners/ 下，向上两级到达 cement model 根目录。
PROJECT_ROOT = Path(__file__).resolve().parents[2]


# CSV 导出协议：只约束本模块实际使用的 to_csv 方法，避免强绑定 pandas 类型。
class _CsvWritable(Protocol):
    """表示具备 to_csv 导出能力的表格对象。"""

    def to_csv(self, path_or_buf: Path, *, index: bool, encoding: str) -> None:
        """将表格写入 CSV 文件。"""


# 表格导出工具：统一使用 UTF-8 BOM，保证 Windows Excel 打开中文列名不乱码。
def _export_table_csv(table: _CsvWritable, path: Path) -> None:
    """以 UTF-8 BOM 编码导出表格，便于 Windows Excel 正确识别中文。"""

    table.to_csv(path, index=False, encoding="utf-8-sig")


# 二维环空模型主导出流程：按呼探1-001井中文命名口径导出完整结果资产。
def run_and_export(
    *,
    mode_title: str,
    output_dir: Path,
    inlet_provider: Callable[[float], AnnulusInletState],
    total_t_s: float,
) -> None:
    """运行呼探1-001井环空二维模型并导出一套中文命名结果。

    参数：
        mode_title: 模式标题（例如“1D2D耦合模型”），用于中文文件名和摘要标题。
        output_dir: 结果输出目录，默认入口使用 results/呼探1-001尾管_1D2D耦合模型/。
        inlet_provider: 由 1D 鞋口出流桥接得到的环空入口边界状态提供器。
        total_t_s: 环空二维顶替计算总时长（秒），由替浆液到达鞋口时刻确定。
    """

    # 创建结果目录并加载呼探1-001井尾管段标准数据。
    output_dir.mkdir(parents=True, exist_ok=True)
    well_spec, fluids, schedule, _ = load_ht1_001_tailpipe()

    # 呼探1-001严格现场模式下，二维求解时长由 1D 套管内前沿追踪确定：
    # 当替浆液第一次到达鞋口时，水泥浆柱已全部进入环空，随后不再让替浆液继续稀释水泥场。
    solver = AnnulusD2DGASolver(total_t=total_t_s)
    result = solver.run(well_spec, fluids, inlet_provider)

    # 导出时间序列与深度剖面 CSV；列名由求解器保持中文口径。
    metrics_path = output_dir / f"呼探1-001尾管_{mode_title}_时间序列结果.csv"
    profiles_path = output_dir / f"呼探1-001尾管_{mode_title}_深度剖面.csv"
    _export_table_csv(cast(_CsvWritable, result.metrics), metrics_path)
    _export_table_csv(cast(_CsvWritable, result.depth_profiles), profiles_path)

    fluid_provenance_summary = build_injected_fluid_provenance_summary(well_spec.well_name, schedule, fluids)
    summary_payload = dict(result.summary)
    summary_payload["注入流体现场符合性检查"] = fluid_provenance_summary

    # 导出 JSON 摘要，保留机器可读的完整结果字典。
    summary_json_path = output_dir / f"呼探1-001尾管_{mode_title}_结果摘要.json"
    _ = summary_json_path.write_text(
        json.dumps(summary_payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    # 导出 Markdown 摘要，字段名严格使用中文，避免混淆顶替效率与质量响应效率。
    summary_md_path = output_dir / f"呼探1-001尾管_{mode_title}_结果摘要.md"
    final_result = cast(Mapping[str, float], result.summary["最终结果"])
    _ = summary_md_path.write_text(
        "\n".join(
            [
                f"# 呼探1-001尾管{mode_title}结果摘要",
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

    # 导出二维云图：深度-时间等值线、多时刻截面快照与最终场图。
    _ = plot_depth_time_contour(result, output_dir=output_dir)
    _ = plot_annulus_snapshots(result, output_dir=output_dir)
    _ = plot_final_fields_contour(result, output_dir=output_dir)

    # 导出 2D 场数据 NPZ。呼探1-001为多浆柱体系，必须包含 lead/tail 快照和最终场。
    npz_path = output_dir / f"呼探1-001尾管_{mode_title}_2D场数据.npz"
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


# 停止时间判定：以替浆液首次到达鞋口作为二维环空计算截止时刻。
def annulus_stop_time_s(
    *,
    casing_result: CasingFlowResult,
    fluids: tuple[FluidSpec, ...],
) -> float:
    """确定呼探1-001井环空二维顶替应停止的地面累计时间。"""

    del fluids
    return float(casing_result.cement_end_time_s)


# 水泥浆名称识别：仅用于时序表中文展示，不替代 loader 中的角色定义。
def _is_cement_slurry_name(fluid_name: str) -> bool:
    """根据中文流体名称识别是否属于水泥浆柱。"""

    # loader 中的角色才是权威数据；本函数仅用于鞋口时序表的中文相名展示。
    cement_keywords = ("领浆", "中间浆", "尾浆", "水泥浆")
    excluded_keywords = ("替浆", "压塞", "平衡", "冲洗", "隔离")
    return any(keyword in fluid_name for keyword in cement_keywords) and not any(
        keyword in fluid_name for keyword in excluded_keywords
    )


# 环空入口相标签：将鞋口出流流体名称映射为审阅时序表中的中文相名。
def _annulus_phase_label(fluid_name: str) -> str:
    """将鞋口出流流体名称转换为时序表中的环空入口相中文标签。"""

    # 水泥浆统一显示为水泥相，隔离/冲洗/平衡液统一显示为前置液/隔离液相。
    if _is_cement_slurry_name(fluid_name):
        return "水泥相"
    if any(keyword in fluid_name for keyword in ("隔离", "冲洗", "平衡")):
        return "前置液/隔离液相"
    return "泥浆相"


# 鞋口时序导出：记录地面泵注时间、鞋口出流流体与环空入口相切换关系。
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

    # 创建输出目录，并使用呼探1-001井专用中文文件名导出 UTF-8 BOM CSV。
    output_dir.mkdir(parents=True, exist_ok=True)
    timing_path = output_dir / "呼探1-001尾管_1D2D耦合模型_鞋口出流时序.csv"

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

    # 计算套管内容积和施工结束时刻，用于地面时间到鞋口出流时间的体积坐标转换。
    pump_end_time_s = elapsed_time_s
    pipe_volume_m3 = casing_result.shoe_md_m * casing_result.pipe_cross_section_m2

    # 内部换算函数：按累计泵入体积找到指定体积坐标对应的地面累计时间。
    def _arrival_time_for_volume(target_volume_m3: float) -> float | None:
        """按施工累计体积反算指定体积坐标到达鞋口的地面累计时间。"""

        # 遍历泵注窗口，定位目标体积落在哪一个施工步骤内。
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

    # 环空时间从首个水泥浆到达鞋口起算；若未识别到水泥浆，则保守回退到停泵时刻。
    annulus_start_time_s = cement_arrival_time_s if cement_arrival_time_s is not None else pump_end_time_s
    sample_times = sorted({0.0, pump_end_time_s, *arrival_times})

    # 写出带中文列名的鞋口出流时序 CSV，用于追溯 1D-2D 边界切换。
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


# 公开运行入口：执行呼探1-001井 139.7mm 完井尾管段现场实录 1D-2D 耦合模型。
def run_ht1_001_tailpipe_initial() -> None:
    """呼探1-001井尾管段现场实录 1D-2D 耦合模型运行入口。"""

    # 加载呼探1-001井尾管段井身、流体、泵注程序与附加元数据。
    well_spec, fluids, schedule, _ = load_ht1_001_tailpipe()
    output_dir = PROJECT_ROOT / "results" / "呼探1-001尾管_1D2D耦合模型"

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

    # 运行二维环空模型并导出 CSV、JSON、Markdown、PNG、NPZ 与 GIF 全套结果。
    run_and_export(
        mode_title="1D2D耦合模型",
        output_dir=output_dir,
        inlet_provider=coupled_provider,
        total_t_s=annulus_stop_time_value_s,
    )


# 命令行入口：允许 python -m cemdisp.runners.ht1_001_tailpipe 直接运行模型。
if __name__ == "__main__":
    run_ht1_001_tailpipe_initial()
