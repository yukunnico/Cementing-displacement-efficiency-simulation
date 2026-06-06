"""参考项目风格图件输出。"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure

from cemdisp.data.well_spec import WellSpec
from cemdisp.models2d import AnnulusSimulationResult
from cemdisp.reporting.plots import _save_figure

# matplotlib/pandas 的类型桩对动态绘图 API 支持有限；本模块通过运行时数据列约定保证绘图输入正确。
# pyright: reportPrivateUsage=false, reportUnknownMemberType=false, reportUnknownArgumentType=false, reportUnusedCallResult=false, reportUnknownVariableType=false, reportAny=false, reportAttributeAccessIssue=false


def plot_reference_efficiency_timeseries(
    result: AnnulusSimulationResult,
    output_dir: Path | str | None = None,
) -> Figure:
    """绘制参考项目口径的顶替效率时程图。"""
    metrics = result.metrics
    fig, ax = plt.subplots(figsize=(9, 5.4))
    ax.plot(metrics["time_min"], metrics["effective_efficiency"], label="全井段有效顶替效率", color="#2E86AB", linewidth=2.2)
    ax.set_title("顶替效率时程", fontsize=14, fontweight="bold")
    ax.set_xlabel("时间 / min", fontsize=11)
    ax.set_ylabel("顶替效率", fontsize=11)
    ax.set_ylim(0.0, 1.05)
    ax.grid(True, alpha=0.2, linestyle="--")
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    _save_figure(fig, output_dir, "顶替效率时程.png")
    return fig


def plot_reference_depth_profile(
    result: AnnulusSimulationResult,
    output_dir: Path | str | None = None,
) -> Figure:
    """绘制参考项目口径的最终水泥体积分数深度剖面图。"""
    profiles = result.depth_profiles.sort_values("井深_m")
    fig, ax = plt.subplots(figsize=(7.5, 8.0))
    ax.plot(profiles["水泥平均浓度"], profiles["井深_m"], label="周向平均", color="black", linewidth=2.3)
    ax.plot(profiles["宽边水泥浓度"], profiles["井深_m"], label="宽边", color="#2E86AB", linewidth=2.0)
    ax.plot(profiles["窄边水泥浓度"], profiles["井深_m"], label="窄边", color="#C73E1D", linewidth=2.0)
    ax.set_title("最终水泥体积分数深度剖面", fontsize=14, fontweight="bold")
    ax.set_xlabel("水泥体积分数", fontsize=11)
    ax.set_ylabel("井深 / m", fontsize=11)
    ax.set_xlim(0.0, 1.05)
    ax.grid(True, alpha=0.2, linestyle="--")
    ax.legend(loc="best", fontsize=9)
    ax.invert_yaxis()
    fig.tight_layout()
    _save_figure(fig, output_dir, "最终水泥体积分数深度剖面.png")
    return fig


def plot_reference_front_progress(
    result: AnnulusSimulationResult,
    well_spec: WellSpec,
    output_dir: Path | str | None = None,
) -> Figure:
    """绘制参考项目口径的宽窄边前沿推进图。"""
    metrics = result.metrics
    wide_md = well_spec.bottom_md_m - metrics["front_wide_m"]
    mid_md = well_spec.bottom_md_m - metrics["front_mid_m"]
    narrow_md = well_spec.bottom_md_m - metrics["front_narrow_m"]

    fig, ax = plt.subplots(figsize=(9, 5.4))
    ax.plot(metrics["time_min"], wide_md, label="宽边前沿", color="#2E86AB", linewidth=2.2)
    ax.plot(metrics["time_min"], mid_md, label="中线前沿", color="#F18F01", linewidth=2.0)
    ax.plot(metrics["time_min"], narrow_md, label="窄边前沿", color="#C73E1D", linewidth=2.2)
    ax.set_title("宽窄边前沿推进", fontsize=14, fontweight="bold")
    ax.set_xlabel("时间 / min", fontsize=11)
    ax.set_ylabel("井深 / m", fontsize=11)
    ax.grid(True, alpha=0.2, linestyle="--")
    ax.legend(loc="best", fontsize=9)
    ax.invert_yaxis()
    fig.tight_layout()
    _save_figure(fig, output_dir, "宽窄边前沿推进.png")
    return fig


def plot_reference_final_cement_field(
    result: AnnulusSimulationResult,
    output_dir: Path | str | None = None,
) -> Figure:
    """绘制参考项目风格的最终环空二维水泥浓度场。"""
    cement = np.asarray(result.cement_field, dtype=float)
    md = np.asarray(result.geom["md"], dtype=float)
    fig, ax = plt.subplots(figsize=(8.8, 6.4))
    image = ax.imshow(
        cement.T,
        aspect="auto",
        origin="upper",
        extent=(0.0, 1.0, float(md.max()), float(md.min())),
        cmap="viridis",
        vmin=0.0,
        vmax=1.0,
    )
    ax.set_title("最终环空二维水泥浓度场", fontsize=14, fontweight="bold")
    ax.set_xlabel("方位归一化坐标（宽边 → 窄边）", fontsize=11)
    ax.set_ylabel("井深 / m", fontsize=11)
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("水泥浓度", fontsize=10)
    fig.tight_layout()
    _save_figure(fig, output_dir, "最终环空二维水泥浓度场.png")
    return fig


def plot_reference_segment_efficiency(
    result: AnnulusSimulationResult,
    well_spec: WellSpec,
    output_dir: Path | str | None = None,
) -> Figure:
    """绘制参考项目同类型的分段最终顶替效率柱状图。"""
    labels: list[str] = ["全井段"]
    values: list[float] = [_window_mean_efficiency(result, well_spec.top_md_m, well_spec.bottom_md_m)]

    sorted_windows = sorted(well_spec.evaluation_windows, key=lambda window: (window.top_md_m, window.bottom_md_m, window.name))
    for window in sorted_windows:
        labels.append(window.name)
        values.append(_window_mean_efficiency(result, window.top_md_m, window.bottom_md_m))

    fig, ax = plt.subplots(figsize=(max(8.0, 1.8 * len(labels)), 5.4))
    bars = ax.bar(labels, values, color="#2E86AB", edgecolor="#333333", linewidth=0.6)
    ax.set_title("分段最终顶替效率", fontsize=14, fontweight="bold")
    ax.set_ylabel("最终有效顶替效率", fontsize=11)
    ax.set_ylim(0.0, 1.05)
    ax.grid(axis="y", alpha=0.2, linestyle="--")
    ax.tick_params(axis="x", rotation=18)
    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            bar.get_height() + 0.015,
            f"{value:.3f}",
            ha="center",
            va="bottom",
            fontsize=10,
        )
    fig.tight_layout()
    _save_figure(fig, output_dir, "分段最终顶替效率.png")
    return fig


def export_reference_figure_set(
    result: AnnulusSimulationResult,
    well_spec: WellSpec,
    output_dir: Path | str | None = None,
) -> dict[str, Figure]:
    """导出一整套参考项目风格图件。"""
    return {
        "顶替效率时程.png": plot_reference_efficiency_timeseries(result, output_dir),
        "最终水泥体积分数深度剖面.png": plot_reference_depth_profile(result, output_dir),
        "宽窄边前沿推进.png": plot_reference_front_progress(result, well_spec, output_dir),
        "最终环空二维水泥浓度场.png": plot_reference_final_cement_field(result, output_dir),
        "分段最终顶替效率.png": plot_reference_segment_efficiency(result, well_spec, output_dir),
    }


def _window_mean_efficiency(result: AnnulusSimulationResult, top_md_m: float, bottom_md_m: float) -> float:
    """按井深窗口计算平均水泥体积分数。"""
    profiles = result.depth_profiles
    mask = (profiles["井深_m"] >= top_md_m) & (profiles["井深_m"] <= bottom_md_m)
    if not bool(mask.any()):
        return float("nan")
    value_column = "平均有效顶替效率" if "平均有效顶替效率" in profiles.columns else "水泥平均浓度"
    return float(profiles.loc[mask, value_column].mean())


__all__ = [
    "export_reference_figure_set",
    "plot_reference_depth_profile",
    "plot_reference_efficiency_timeseries",
    "plot_reference_final_cement_field",
    "plot_reference_front_progress",
    "plot_reference_segment_efficiency",
]
