"""
固井顶替效率模型的二维云图输出。

本模块提供面向报告与论文插图的二维场可视化函数，包括：
1. 深度-时间顶替效率云图；
2. 环空水泥浓度场多时刻快照图；
3. 最终水泥浓度、有效顶替效率与壁面泥饼场分布图。

所有图表均沿用 :mod:`cemdisp.reporting.plots` 中的中文字体、保存逻辑与安全文件名规则，
并避免在通用绘图函数中硬编码任何单井专用数据。
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from numpy.typing import NDArray

# matplotlib 的类型桩对动态绘图 API 支持有限；本模块通过运行时校验保证数组维度正确。
# pyright: reportAny=false, reportPrivateUsage=false, reportUnknownMemberType=false, reportUnknownVariableType=false, reportUnusedCallResult=false

from cemdisp.models2d import AnnulusSimulationResult
from cemdisp.reporting.plots import _safe_filename_component, _save_figure, _setup_chinese_font

Array = NDArray[np.float64]


def _require_snapshot_data(result: AnnulusSimulationResult) -> tuple[tuple[Array, ...], Array]:
    """读取并校验水泥浓度快照及其时间坐标。"""
    snapshots = tuple(np.asarray(snapshot, dtype=float) for snapshot in result.cement_snapshots)
    if not snapshots:
        raise ValueError("result.cement_snapshots 为空，无法绘制二维场快照云图。")

    # 优先使用专门的快照时间；若旧结果对象未提供，则退回到 time_points_s。
    raw_times = result.snapshot_times_s or result.time_points_s
    if len(raw_times) != len(snapshots):
        raise ValueError("快照数量与快照时间数量不一致，无法建立时间坐标。")

    times_min = np.asarray(raw_times, dtype=float) / 60.0
    return snapshots, times_min


def _depth_coordinates(result: AnnulusSimulationResult) -> Array:
    """从几何字典中读取井深坐标。"""
    if "md" not in result.geom:
        raise KeyError("result.geom 缺少井深坐标字段 'md'。")
    return np.asarray(result.geom["md"], dtype=float)


def _azimuth_coordinates(result: AnnulusSimulationResult) -> Array:
    """从几何字典中读取方位坐标，用于校验二维场尺寸。"""
    if "y" not in result.geom:
        raise KeyError("result.geom 缺少方位坐标字段 'y'。")
    return np.asarray(result.geom["y"], dtype=float)


def _validate_field_shape(field: Array, depth: Array, azimuth: Array, field_name: str) -> None:
    """检查二维场是否符合 ny×nz 的环空网格约定。"""
    expected_shape = (azimuth.size, depth.size)
    if field.shape != expected_shape:
        raise ValueError(f"{field_name} 形状应为 {expected_shape}，实际为 {field.shape}。")


def _field_extent(depth: Array) -> tuple[float, float, float, float]:
    """返回 imshow 使用的绘图范围，方位坐标采用宽边到窄边的归一化显示。"""
    return (float(depth.min()), float(depth.max()), 0.0, 1.0)


def plot_depth_time_contour(
    result: AnnulusSimulationResult,
    output_dir: Path | str | None = None,
) -> Figure:
    """绘制深度-时间有效顶替效率云图。

    每个时间快照先沿方位角方向取平均，得到对应时刻的深度剖面；随后按
    ``(n_depths, n_snapshots)`` 组装为二维矩阵，用于展示水泥浓度随时间和井深的推进过程。

    Args:
        result: 环空二维模拟结果，需包含 ``cement_snapshots`` 与 ``snapshot_times_s``。
        output_dir: 可选输出目录；若提供则自动保存 PNG 图表。

    Returns:
        matplotlib Figure 对象。
    """
    required_columns = {"time_min", "effective_efficiency", "front_wide_m", "front_narrow_m"}
    missing_columns = required_columns.difference(result.metrics.columns)
    if missing_columns:
        raise KeyError(f"metrics 缺少绘图所需字段：{', '.join(sorted(missing_columns))}")

    snapshots, times_min = _require_snapshot_data(result)
    depth = _depth_coordinates(result)
    azimuth = _azimuth_coordinates(result)
    for index, snapshot in enumerate(snapshots):
        _validate_field_shape(snapshot, depth, azimuth, f"cement_snapshots[{index}]")

    # 沿方位角方向平均，得到每个快照的井深方向水泥浓度剖面。
    depth_time_matrix = np.column_stack([np.mean(snapshot, axis=0) for snapshot in snapshots])
    time_grid, depth_grid = np.meshgrid(times_min, depth)

    well_name = _safe_filename_component(result.well_name)
    fig, ax = plt.subplots(figsize=(10.5, 6.4))
    filled = ax.contourf(
        time_grid,
        depth_grid,
        depth_time_matrix,
        levels=np.linspace(0.0, 1.0, 21),
        cmap="RdYlGn",
        vmin=0.0,
        vmax=1.0,
    )
    contour_lines = ax.contour(
        time_grid,
        depth_grid,
        depth_time_matrix,
        levels=[0.5, 0.7, 0.9],
        colors="#333333",
        linewidths=0.8,
    )
    ax.clabel(contour_lines, fmt="%.1f", inline=True, fontsize=9)

    colorbar = fig.colorbar(filled, ax=ax, pad=0.02)
    colorbar.set_label("有效顶替效率")
    ax.set_title(f"{result.well_name} 深度-时间顶替效率云图", fontsize=14, fontweight="bold")
    ax.set_xlabel("时间 / min")
    ax.set_ylabel("井深 / m")
    ax.grid(True, alpha=0.18)

    fig.tight_layout()
    _save_figure(fig, output_dir, f"{well_name}_深度-时间顶替效率云图.png")
    return fig


def plot_annulus_snapshots(
    result: AnnulusSimulationResult,
    output_dir: Path | str | None = None,
    n_panels: int = 6,
) -> Figure:
    """绘制环空水泥浓度场多时刻快照图。

    Args:
        result: 环空二维模拟结果，需包含水泥浓度场快照。
        output_dir: 可选输出目录；若提供则自动保存 PNG 图表。
        n_panels: 展示的快照面板数量，按时间均匀抽取。

    Returns:
        matplotlib Figure 对象。
    """
    if n_panels <= 0:
        raise ValueError("n_panels 必须为正整数。")

    snapshots, times_min = _require_snapshot_data(result)
    depth = _depth_coordinates(result)
    azimuth = _azimuth_coordinates(result)
    for index, snapshot in enumerate(snapshots):
        _validate_field_shape(snapshot, depth, azimuth, f"cement_snapshots[{index}]")

    panel_count = min(n_panels, len(snapshots))
    selected_indices = np.linspace(0, len(snapshots) - 1, panel_count, dtype=int)
    selected_indices = np.unique(selected_indices)
    panel_count = int(selected_indices.size)
    n_cols = min(3, panel_count)
    n_rows = int(np.ceil(panel_count / n_cols))

    well_name = _safe_filename_component(result.well_name)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.8 * n_cols, 3.6 * n_rows), squeeze=False)
    image = None
    for panel_index, snapshot_index in enumerate(selected_indices):
        row, col = divmod(panel_index, n_cols)
        ax = axes[row][col]
        # 原始数组按窄边到宽边存储；翻转后按“宽边→窄边”展示。
        image = ax.imshow(
            snapshots[int(snapshot_index)][::-1, :],
            extent=_field_extent(depth),
            origin="lower",
            aspect="auto",
            cmap="viridis",
            vmin=0.0,
            vmax=1.0,
        )
        ax.set_title(f"t = {times_min[int(snapshot_index)]:.1f} min")
        ax.set_xlabel("井深 / m")
        ax.set_ylabel("方位角位置 (宽边→窄边)")

    for empty_index in range(panel_count, n_rows * n_cols):
        row, col = divmod(empty_index, n_cols)
        axes[row][col].axis("off")

    fig.suptitle(f"{result.well_name} 水泥浓度场演化过程", fontsize=15, fontweight="bold")
    if image is not None:
        colorbar = fig.colorbar(image, ax=axes.ravel().tolist(), pad=0.02, shrink=0.92)
        colorbar.set_label("水泥浓度")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    _save_figure(fig, output_dir, f"{well_name}_水泥浓度场演化过程.png")
    return fig


def plot_final_fields_contour(
    result: AnnulusSimulationResult,
    output_dir: Path | str | None = None,
) -> Figure:
    """绘制最终水泥浓度、有效顶替效率和壁面泥饼场分布图。

    Args:
        result: 环空二维模拟结果，需包含最终水泥浓度场和壁面泥饼场。
        output_dir: 可选输出目录；若提供则自动保存 PNG 图表。

    Returns:
        matplotlib Figure 对象。
    """
    depth = _depth_coordinates(result)
    azimuth = _azimuth_coordinates(result)
    cement_field = np.asarray(result.cement_field, dtype=float)
    wall_field = np.asarray(result.wall_field, dtype=float)
    effective_field = cement_field * (1.0 - wall_field)

    _validate_field_shape(cement_field, depth, azimuth, "cement_field")
    _validate_field_shape(wall_field, depth, azimuth, "wall_field")

    panels = [
        ("水泥浓度", cement_field, "viridis"),
        ("有效顶替效率", effective_field, "RdYlGn"),
        ("壁面泥饼", wall_field, "YlOrRd"),
    ]
    well_name = _safe_filename_component(result.well_name)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), sharex=True, sharey=True)

    for ax, (title, field, colormap) in zip(axes, panels):
        # 按宽边到窄边显示，避免用户读图时误解方位方向。
        image = ax.imshow(
            field[::-1, :],
            extent=_field_extent(depth),
            origin="lower",
            aspect="auto",
            cmap=colormap,
            vmin=0.0,
            vmax=1.0,
        )
        ax.set_title(title, fontweight="bold")
        ax.set_xlabel("井深 / m")
        ax.set_ylabel("方位角 (宽边→窄边)")
        colorbar = fig.colorbar(image, ax=ax, pad=0.02)
        colorbar.set_label(title)

    fig.suptitle(f"{result.well_name} 最终场分布", fontsize=15, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    _save_figure(fig, output_dir, f"{well_name}_最终场分布.png")
    return fig


_setup_chinese_font()


__all__ = [
    "plot_annulus_snapshots",
    "plot_depth_time_contour",
    "plot_final_fields_contour",
]
