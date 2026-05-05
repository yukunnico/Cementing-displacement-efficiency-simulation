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

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
import pandas as pd

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.loaders import build_hu101_annulus_inlet_provider, load_hu101_tailpipe
from cemdisp.data.pumping_schedule import PumpingSchedule
from cemdisp.models2d import AnnulusD2DGASolver, AnnulusSimulationResult
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
from cemdisp.transport1d import CasingFlowResult, CasingFlowSolver


_CEMDISP_ROOT = Path(__file__).resolve().parents[1]  # cemdisp/
PROJECT_ROOT = _CEMDISP_ROOT.parent                  # cement model/
HU101_CBL_PROFILE_PATH = PROJECT_ROOT / "参考文档" / "呼101" / "提取数据" / "100312_CBL剖面_Excel版.csv"


def _schedule_total_time_s(schedule: PumpingSchedule) -> float:
    """按现场分段排量计算施工总时长。"""

    return sum(0.0 if step.rate_m3_min <= 0.0 else step.volume_m3 / step.rate_m3_min * 60.0 for step in schedule.steps)


def annulus_stop_time_s(
    *,
    casing_result: CasingFlowResult,
    fluids: tuple[FluidSpec, ...],
) -> float:
    """返回 Hu101 环空二维顶替应停止的地面累计时间。

    对 Hu101 coupled 模式，停止条件定义为：最后一段水泥浆之后的第一种
    DISPLACEMENT 流体第一次到达鞋口。此时整段水泥浆已全部进入环空，
    不再继续让后续替浆流体进入环空稀释既有水泥场。
    """

    role_by_name = {fluid.name: fluid.role for fluid in fluids}
    for front in casing_result.fronts:
        if role_by_name.get(front.fluid_name) == FluidRole.DISPLACEMENT:
            return float(front.time_s)
    raise ValueError("Hu101 现场耦合模型未找到替浆流体到鞋口时刻，无法确定环空顶替停止时间")


def _load_hu101_cbl_profile() -> pd.DataFrame:
    """读取 Hu101 CBL 剖面代理 CSV，并保留用于同深度对比的数值字段。"""

    columns = [
        "row_id",
        "depth_md_m",
        "cbl_amplitude_pct",
        "quality_proxy",
        "quality_proxy_pct",
        "quality_grade_cn",
        "segment_type_cn",
        "segment_note_cn",
        "is_target_interval_cn",
        "is_double_liner_excluded_cn",
    ]
    last_error: UnicodeDecodeError | None = None
    for encoding in ("utf-8-sig", "gbk", "gb18030"):
        try:
            cbl_profile = pd.read_csv(HU101_CBL_PROFILE_PATH, skiprows=2, names=columns, encoding=encoding)
            break
        except UnicodeDecodeError as exc:
            last_error = exc
    else:
        raise RuntimeError(f"无法读取 Hu101 CBL 剖面 CSV: {HU101_CBL_PROFILE_PATH}") from last_error

    for column in ("depth_md_m", "cbl_amplitude_pct", "quality_proxy", "quality_proxy_pct"):
        cbl_profile[column] = pd.to_numeric(cbl_profile[column], errors="coerce")
    return cbl_profile.dropna(subset=["depth_md_m", "quality_proxy"]).reset_index(drop=True)


def _export_cbl_profile_comparison(result: AnnulusSimulationResult, output_dir: Path) -> None:
    """导出 Hu101 模拟剖面与 CBL 代理剖面的逐深度对比结果。"""

    cbl_profile = _load_hu101_cbl_profile()
    profile = result.depth_profiles.sort_values("井深_m")

    simulated = np.interp(
        cbl_profile["depth_md_m"].to_numpy(dtype=float),
        profile["井深_m"].to_numpy(dtype=float),
        profile["平均有效顶替效率"].to_numpy(dtype=float),
    )
    observed = cbl_profile["quality_proxy"].to_numpy(dtype=float)
    diff = simulated - observed

    comparison = cbl_profile[[
        "depth_md_m",
        "cbl_amplitude_pct",
        "quality_proxy",
        "quality_proxy_pct",
        "quality_grade_cn",
        "segment_type_cn",
        "segment_note_cn",
    ]].copy()
    comparison["simulated_effective_efficiency"] = simulated
    comparison["simulated_minus_cbl_proxy"] = diff
    comparison_path = output_dir / "呼101尾管_模拟剖面_vs_CBL代理剖面.csv"
    comparison.to_csv(comparison_path, index=False, encoding="utf-8-sig")

    severe_mask = (cbl_profile["depth_md_m"] >= 6050.0) & (cbl_profile["depth_md_m"] <= 6210.0)
    anomaly_mask = (cbl_profile["depth_md_m"] >= 6100.0) & (cbl_profile["depth_md_m"] <= 6400.0)
    target_mask = (cbl_profile["depth_md_m"] >= 6153.0) & (cbl_profile["depth_md_m"] <= 7741.0)

    def _mean_for(mask: pd.Series, array: NDArray[np.float64]) -> float:
        mask_array = np.asarray(mask, dtype=bool)
        return float(np.mean(array[mask_array])) if np.any(mask_array) else float("nan")

    diagnostics = {
        "CBL剖面对比点数": float(len(cbl_profile)),
        "模拟_minus_CBL代理_平均差": float(np.mean(diff)),
        "模拟_minus_CBL代理_MAE": float(np.mean(np.abs(diff))),
        "模拟_minus_CBL代理_RMSE": float(np.sqrt(np.mean(diff**2))),
        "CBL代理_全剖面均值": float(np.mean(observed)),
        "模拟剖面_全剖面均值": float(np.mean(simulated)),
        "6050_6210m_CBL代理均值": _mean_for(severe_mask, observed),
        "6050_6210m_模拟均值": _mean_for(severe_mask, simulated),
        "6100_6400m_CBL代理均值": _mean_for(anomaly_mask, observed),
        "6100_6400m_模拟均值": _mean_for(anomaly_mask, simulated),
        "6153_7741m_CBL代理均值": _mean_for(target_mask, observed),
        "6153_7741m_模拟均值": _mean_for(target_mask, simulated),
    }

    diagnostics_json_path = output_dir / "呼101尾管_CBL剖面对比诊断.json"
    diagnostics_json_path.write_text(json.dumps(diagnostics, ensure_ascii=False, indent=2), encoding="utf-8")

    diagnostics_md_path = output_dir / "呼101尾管_CBL剖面对比诊断.md"
    diagnostics_md_path.write_text(
        "\n".join(
            [
                "# 呼101尾管模拟剖面与 CBL 代理剖面对比诊断",
                "",
                f"- CBL剖面对比点数：{diagnostics['CBL剖面对比点数']:.0f}",
                f"- 模拟-CBL平均差：{diagnostics['模拟_minus_CBL代理_平均差']:.4f}",
                f"- MAE：{diagnostics['模拟_minus_CBL代理_MAE']:.4f}",
                f"- RMSE：{diagnostics['模拟_minus_CBL代理_RMSE']:.4f}",
                f"- CBL代理全剖面均值：{diagnostics['CBL代理_全剖面均值']:.4f}",
                f"- 模拟全剖面均值：{diagnostics['模拟剖面_全剖面均值']:.4f}",
                f"- 6050–6210m：模拟 {diagnostics['6050_6210m_模拟均值']:.4f} / CBL {diagnostics['6050_6210m_CBL代理均值']:.4f}",
                f"- 6100–6400m：模拟 {diagnostics['6100_6400m_模拟均值']:.4f} / CBL {diagnostics['6100_6400m_CBL代理均值']:.4f}",
                f"- 6153–7741m：模拟 {diagnostics['6153_7741m_模拟均值']:.4f} / CBL {diagnostics['6153_7741m_CBL代理均值']:.4f}",
            ]
        ),
        encoding="utf-8",
    )

    chart_path = output_dir / "呼101尾管_模拟剖面_vs_CBL代理剖面.png"
    plt.figure(figsize=(7, 8))
    plt.plot(simulated, cbl_profile["depth_md_m"], label="模拟有效顶替效率", linewidth=2.2)
    plt.plot(observed, cbl_profile["depth_md_m"], label="CBL质量代理", linewidth=2.2)
    plt.axhspan(6050.0, 6210.0, color="#EF4444", alpha=0.10, label="6050–6210m 强异常段")
    plt.axhspan(6100.0, 6400.0, color="#F59E0B", alpha=0.08, label="6100–6400m 异常段")
    plt.gca().invert_yaxis()
    plt.xlabel("效率 / 质量代理")
    plt.ylabel("井深 / m")
    plt.title("呼101尾管模拟剖面与 CBL 代理剖面逐深度对比")
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(chart_path, dpi=220)
    plt.close()


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
    # 另外沿用 Hu101 legacy 模型口径：
    # 1) nz=500：legacy 收敛检查以 500 作为主计算网格；
    # 2) quality_penalty_scale=0.671：按 Hu101 现场 CBL 合格率 62.77% 做单井校准。
    total_t = total_t_s if total_t_s is not None else _schedule_total_time_s(schedule) + 20.0 * 60.0
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
    _export_cbl_profile_comparison(result, output_dir)

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
