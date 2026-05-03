"""
呼102尾管段顶替效率模型运行器

运行两种模式并导出中文命名结果：
1. 硬编码环空入口(sustained_tail)：替浆期间环空入口保持尾浆
2. 1D-2D耦合：套管内前沿追踪 → 鞋口出流 → 环空入口

输出目录：results/呼102尾管_初版模型/ 和 results/呼102尾管_1D2D耦合模型/
输出文件：CSV(时间序列/深度剖面) + JSON(摘要) + Markdown(摘要) + PNG(静态图) + NPZ(2D场数据) + GIF(动画)
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
import json
from pathlib import Path
from typing import cast

import numpy as np

from cemdisp.data.loaders import build_hu102_annulus_inlet_provider, load_hu102_tailpipe
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


# 项目根目录：从cemdisp/runners向上两级到达cement model根目录
_CEMDISP_ROOT = Path(__file__).resolve().parents[1]  # cemdisp/
PROJECT_ROOT = _CEMDISP_ROOT.parent                  # cement model/


def run_and_export(
    *,
    mode_title: str,
    output_dir: Path,
    inlet_provider: Callable[[float], AnnulusInletState],
) -> None:
    """运行环空模型并导出一套中文命名结果。

    参数：
        mode_title: 模式标题（如"初版模型"、"1D2D耦合模型"），用于文件名和打印
        output_dir: 结果输出目录
        inlet_provider: 环空入口边界状态提供器
    """

    output_dir.mkdir(parents=True, exist_ok=True)
    well_spec, fluids, _, _ = load_hu102_tailpipe()
    solver = AnnulusD2DGASolver()
    result = solver.run(well_spec, fluids, inlet_provider)

    # 导出CSV
    metrics_path = output_dir / f"呼102尾管_{mode_title}_时间序列结果.csv"
    profiles_path = output_dir / f"呼102尾管_{mode_title}_深度剖面.csv"
    _ = result.metrics.to_csv(metrics_path, index=False, encoding="utf-8-sig")
    _ = result.depth_profiles.to_csv(profiles_path, index=False, encoding="utf-8-sig")

    # 导出JSON摘要
    summary_json_path = output_dir / f"呼102尾管_{mode_title}_结果摘要.json"
    _ = summary_json_path.write_text(
        json.dumps(result.summary, ensure_ascii=False, indent=2), encoding="utf-8"
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
                f"- CBL评价井段模拟有效顶替效率：{final_result['CBL评价井段模拟有效顶替效率']:.4f}",
                f"- 目标层段模拟有效顶替效率：{final_result['目标层段模拟有效顶替效率']:.4f}",
                f"- 最终水泥浆占据率：{final_result['最终水泥浆占据率']:.4f}",
                f"- 最终质量响应效率：{final_result['最终质量响应效率']:.4f}",
                f"- 最终窜槽/混浆/失稳指数：{final_result['最终窜槽指数']:.4f} / {final_result['最终混浆指数']:.4f} / {final_result['最终失稳指数']:.4f}",
            ]
        ),
        encoding="utf-8",
    )

    # 导出静态图表（中文标签+中文文件名）
    _ = plot_time_series(result, output_dir=output_dir)
    _ = plot_depth_profiles(result, well_spec=well_spec, output_dir=output_dir)
    _ = plot_risk_indices(result, output_dir=output_dir)
    _ = plot_efficiency_summary_bar(result, output_dir=output_dir)

    # 导出云图（深度-时间等值线 + 多时刻截面快照 + 最终场三联图）
    _ = plot_depth_time_contour(result, output_dir=output_dir)
    _ = plot_annulus_snapshots(result, output_dir=output_dir)
    _ = plot_final_fields_contour(result, output_dir=output_dir)

    # 导出2D场数据NPZ（水泥/隔离液/泥饼/触变/湍流快照 + 时间点 + 网格坐标）
    # 这些数组直接来自求解器结果对象，确保导出数据与模型实际计算场一致。
    npz_path = output_dir / f"呼102尾管_{mode_title}_2D场数据.npz"
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
    print(json.dumps(result.summary, ensure_ascii=False, indent=2))


def run_hu102_tailpipe_initial() -> None:
    """呼102尾管段顶替效率模型完整运行入口。"""

    well_spec, fluids, schedule, _ = load_hu102_tailpipe()

    # 硬编码环空入口模式（sustained_tail）：向后兼容，便于对比历史结果。
    hardcoded_provider = build_hu102_annulus_inlet_provider(schedule, fluids)
    run_and_export(
        mode_title="初版模型",
        output_dir=PROJECT_ROOT / "results" / "呼102尾管_初版模型",
        inlet_provider=hardcoded_provider,
    )

    # 1D-2D耦合模式：套管内前沿追踪 → 鞋口出流 → 环空入口。
    casing_solver = CasingFlowSolver()
    casing_result = casing_solver.run(well_spec, fluids, schedule)
    coupled_provider = build_coupled_annulus_inlet_provider(casing_result, casing_solver, fluids)
    run_and_export(
        mode_title="1D2D耦合模型",
        output_dir=PROJECT_ROOT / "results" / "呼102尾管_1D2D耦合模型",
        inlet_provider=coupled_provider,
    )
