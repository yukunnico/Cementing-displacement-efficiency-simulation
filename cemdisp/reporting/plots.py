"""
固井顶替效率模型的通用中文图表输出

本模块提供固井顶替模型的各类图表绘制功能，包括：
1. 时间序列图（plot_time_series）：有效顶替效率与前沿推进
2. 深度剖面图（plot_depth_profiles）：深度方向效率与浓度分布
3. 风险指标图（plot_risk_indices）：窜槽、混浆、失稳风险演变
4. 结果对比图（plot_efficiency_summary_bar）：最终各项指标对比

所有图表均使用中文标签和标题，适配常见 Windows 中文字体。
"""

from __future__ import annotations

import logging
from collections.abc import Iterable
from pathlib import Path
from typing import Optional

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import pandas as pd

from cemdisp.models2d import AnnulusSimulationResult

LOGGER = logging.getLogger(__name__)

# 学术图表配色方案
ACADEMIC_COLORS = {
    "primary": "#2E86AB",      # 主色调（蓝色）
    "success": "#A23B72",      # 成功/正向（紫红色）
    "warning": "#F18F01",      # 警告（橙色）
    "danger": "#C73E1D",       # 危险（红色）
    "info": "#3B1F2B",         # 信息（深色）
    "light": "#44BBA4",        # 浅色（绿色）
}

# 标准图表尺寸
STANDARD_FIGSIZE = {
    "single": (8, 6),          # 单图
    "double": (12, 8),         # 双图
    "wide": (14, 6),           # 宽图
    "tall": (8, 10),           # 高图
}


def _setup_chinese_font() -> None:
    """配置 matplotlib 中文字体与负号显示。

    优先查找 SimHei、Microsoft YaHei、SimSun、FangSong 等常见 Windows 字体。
    若未找到中文字体，回退到 sans-serif 并记录警告日志。
    """
    preferred_fonts = ["SimHei", "Microsoft YaHei", "SimSun", "FangSong"]
    available_fonts = {font.name.lower(): font.name for font in fm.fontManager.ttflist}

    selected_font: Optional[str] = None
    for font_name in preferred_fonts:
        if font_name.lower() in available_fonts:
            selected_font = available_fonts[font_name.lower()]
            break

    if selected_font is None:
        plt.rcParams["font.sans-serif"] = ["sans-serif", "DejaVu Sans", "dejavusans"]
        LOGGER.warning("未找到常见 Windows 中文字体，已回退到 sans-serif/DejaVu Sans；中文可能无法完整显示。")
    else:
        plt.rcParams["font.sans-serif"] = [selected_font, "DejaVu Sans", "sans-serif"]

    plt.rcParams["axes.unicode_minus"] = False


def _safe_filename_component(name: str) -> str:
    """返回适合作为文件名前缀的井名。

    移除文件名不允许的字符（如<>:"/\\|?*），
    确保生成的图表文件名合法。
    """
    invalid_chars = '<>:"/\\|?*'
    cleaned = "".join("_" if char in invalid_chars else char for char in name).strip()
    return cleaned or "固井模型"


def _save_figure(fig: Figure, output_dir: Optional[Path | str], filename: str) -> None:
    """按统一参数保存图片。

    Args:
        fig: matplotlib 图表对象
        output_dir: 输出目录路径
        filename: 文件名
    """
    if output_dir is None:
        return

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path / filename, dpi=300, bbox_inches="tight")


def _require_columns(dataframe: pd.DataFrame, columns: Iterable[str], source_name: str) -> None:
    """检查绘图所需列是否存在。

    若缺少必需字段，抛出 KeyError 并列出缺失的列名。
    """
    missing = [column for column in columns if column not in dataframe.columns]
    if missing:
        raise KeyError(f"{source_name} 缺少绘图所需字段：{', '.join(missing)}")


def plot_time_series(result: AnnulusSimulationResult, output_dir: Optional[Path | str] = None) -> Figure:
    """绘制有效顶替效率与前沿推进时间序列。

    输出两张子图：
    - 上图：全井段有效顶替效率及水泥浆占据率随时间变化
    - 下图：宽边、中线、窄边三个方位的水泥浆前沿推进距离随时间变化

    Args:
        result: 环空二维模拟结果
        output_dir: 可选的输出目录路径，若提供则自动保存图表

    Returns:
        matplotlib Figure 对象
    """
    metrics = result.metrics
    _require_columns(
        metrics,
        [
            "time_min",
            "effective_efficiency",
            "bulk_cement_fill",
            "front_wide_m",
            "front_mid_m",
            "front_narrow_m",
        ],
        "metrics",
    )

    well_name = _safe_filename_component(result.well_name)
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle(f"{result.well_name} 顶替效率与前沿推进时间序列", fontsize=15, fontweight="bold")

    axes[0].plot(metrics["time_min"], metrics["effective_efficiency"], label="全井段有效顶替效率", color=ACADEMIC_COLORS["primary"], linewidth=2.2)
    axes[0].plot(metrics["time_min"], metrics["bulk_cement_fill"], label="水泥浆占据率", color=ACADEMIC_COLORS["success"], linestyle="--", linewidth=2)
    axes[0].set_ylabel("有效顶替效率 / 占据率", fontsize=11)
    axes[0].set_ylim(0, 1.05)
    axes[0].grid(True, alpha=0.2, linestyle="--")
    axes[0].legend(loc="best", fontsize=9)

    axes[1].plot(metrics["time_min"], metrics["front_wide_m"], label="宽边", color=ACADEMIC_COLORS["primary"], linewidth=2.1)
    axes[1].plot(metrics["time_min"], metrics["front_mid_m"], label="中线", color=ACADEMIC_COLORS["warning"], linewidth=2.1)
    axes[1].plot(metrics["time_min"], metrics["front_narrow_m"], label="窄边", color=ACADEMIC_COLORS["danger"], linewidth=2.1)
    axes[1].set_xlabel("时间 / min", fontsize=11)
    axes[1].set_ylabel("前沿推进距离 / m", fontsize=11)
    axes[1].grid(True, alpha=0.2, linestyle="--")
    axes[1].legend(loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=9)

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    _save_figure(fig, output_dir, f"{well_name}_顶替效率时间序列.png")
    return fig


def plot_depth_profiles(
    result: AnnulusSimulationResult,
    well_spec: Optional[object] = None,
    output_dir: Optional[Path | str] = None,
) -> Figure:
    """绘制深度方向有效顶替效率与水泥浓度剖面。"""
    del well_spec  # 预留给后续按井段标注，不在通用函数中硬编码单井边界。
    profiles = result.depth_profiles
    _require_columns(
        profiles,
        [
            "井深_m",
            "宽边有效效率",
            "中线有效效率",
            "窄边有效效率",
            "平均有效顶替效率",
            "宽边水泥浓度",
            "中线水泥浓度",
            "窄边水泥浓度",
            "水泥平均浓度",
        ],
        "depth_profiles",
    )

    well_name = _safe_filename_component(result.well_name)
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle(f"{result.well_name} 深度剖面分布", fontsize=15, fontweight="bold")

    depth = profiles["井深_m"]
    axes[0].plot(depth, profiles["宽边有效效率"], label="宽边有效效率", color=ACADEMIC_COLORS["primary"], linewidth=2)
    axes[0].plot(depth, profiles["中线有效效率"], label="中线有效效率", color=ACADEMIC_COLORS["warning"], linewidth=2)
    axes[0].plot(depth, profiles["窄边有效效率"], label="窄边有效效率", color=ACADEMIC_COLORS["danger"], linewidth=2)
    axes[0].plot(depth, profiles["平均有效顶替效率"], label="平均有效顶替效率", color="black", linewidth=2.4)
    axes[0].set_ylabel("有效顶替效率", fontsize=11)
    axes[0].set_ylim(0, 1.05)
    axes[0].grid(True, alpha=0.2, linestyle="--")
    axes[0].legend(loc="upper left", bbox_to_anchor=(1.02, 1), ncol=1, fontsize=9)

    axes[1].plot(depth, profiles["宽边水泥浓度"], label="宽边水泥浓度", color=ACADEMIC_COLORS["primary"], linewidth=2)
    axes[1].plot(depth, profiles["中线水泥浓度"], label="中线水泥浓度", color=ACADEMIC_COLORS["warning"], linewidth=2)
    axes[1].plot(depth, profiles["窄边水泥浓度"], label="窄边水泥浓度", color=ACADEMIC_COLORS["danger"], linewidth=2)
    axes[1].plot(depth, profiles["水泥平均浓度"], label="水泥平均浓度", color="black", linewidth=2.4)
    axes[1].set_xlabel("井深 / m", fontsize=11)
    axes[1].set_ylabel("水泥浓度", fontsize=11)
    axes[1].set_ylim(0, 1.05)
    axes[1].grid(True, alpha=0.2, linestyle="--")
    axes[1].legend(loc="upper left", bbox_to_anchor=(1.02, 1), ncol=1, fontsize=9)

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    _save_figure(fig, output_dir, f"{well_name}_深度剖面分布.png")
    return fig


def plot_risk_indices(result: AnnulusSimulationResult, output_dir: Optional[Path | str] = None) -> Figure:
    """绘制窜槽、混浆和失稳风险指标时间演变。

    绘制三类风险指数随时间的变化曲线：
    - 窜槽指数：宽窄边水泥前沿差异（偏心导致的顶替不均匀）
    - 混浆指数：水泥-钻井液混合程度（顶替界面的混合带）
    - 失稳指数：浮力导致的窄边水泥滑塌风险

    Args:
        result: 环空二维模拟结果
        output_dir: 可选的输出目录路径，若提供则自动保存图表

    Returns:
        matplotlib Figure 对象
    """
    metrics = result.metrics
    _require_columns(
        metrics,
        ["time_min", "channeling_index", "mixing_index", "instability_index"],
        "metrics",
    )

    well_name = _safe_filename_component(result.well_name)
    fig, ax = plt.subplots(figsize=(10, 5.6))
    ax.plot(metrics["time_min"], metrics["channeling_index"], label="窜槽指数", color=ACADEMIC_COLORS["danger"], linewidth=2.2)
    ax.plot(metrics["time_min"], metrics["mixing_index"], label="混浆指数", color=ACADEMIC_COLORS["warning"], linewidth=2.2)
    ax.plot(metrics["time_min"], metrics["instability_index"], label="失稳指数", color=ACADEMIC_COLORS["info"], linewidth=2.2)
    ax.set_title(f"{result.well_name} 风险指标时间演变", fontsize=14, fontweight="bold")
    ax.set_xlabel("时间 / min", fontsize=11)
    ax.set_ylabel("指数", fontsize=11)
    ax.grid(True, alpha=0.2, linestyle="--")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=9)

    fig.tight_layout()
    _save_figure(fig, output_dir, f"{well_name}_风险指标时间演变.png")
    return fig


def plot_efficiency_summary_bar(
    result: AnnulusSimulationResult,
    output_dir: Optional[Path | str] = None,
) -> Figure:
    """绘制最终顶替指标对比柱状图。

    绘制两个最终指标的柱状对比图：
    - 全井段有效顶替效率
    - 水泥浆占据率（水泥体积占环空体积的比例）

    所有指标均基于全井段计算，不区分 CBL 评价段或目标层段。

    Args:
        result: 环空二维模拟结果
        output_dir: 可选的输出目录路径，若提供则自动保存图表

    Returns:
        matplotlib Figure 对象
    """
    metrics = result.metrics
    _require_columns(
        metrics,
        [
            "effective_efficiency",
            "bulk_cement_fill",
        ],
        "metrics",
    )
    final_row = metrics.iloc[-1]
    labels = [
        "全井段有效顶替效率",
        "水泥浆占据率",
    ]
    values = [
        final_row["effective_efficiency"],
        final_row["bulk_cement_fill"],
    ]
    colors = [ACADEMIC_COLORS["primary"], ACADEMIC_COLORS["success"]]

    well_name = _safe_filename_component(result.well_name)
    fig, ax = plt.subplots(figsize=STANDARD_FIGSIZE["single"])
    bars = ax.bar(labels, values, color=colors, edgecolor="#333333", linewidth=0.6)
    ax.set_title(f"{result.well_name} 最终结果对比", fontsize=14, fontweight="bold")
    ax.set_ylabel("效率 / 占据率", fontsize=11)
    ax.set_ylim(0, 1.08)
    ax.grid(axis="y", alpha=0.2, linestyle="--")
    ax.tick_params(axis="x", rotation=18)

    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.015,
            f"{float(value):.3f}",
            ha="center",
            va="bottom",
            fontsize=10,
        )

    fig.tight_layout()
    _save_figure(fig, output_dir, f"{well_name}_最终结果对比.png")
    return fig


_setup_chinese_font()


__all__ = [
    "ACADEMIC_COLORS",
    "STANDARD_FIGSIZE",
    "_save_figure",
    "_setup_chinese_font",
    "plot_depth_profiles",
    "plot_efficiency_summary_bar",
    "plot_risk_indices",
    "plot_time_series",
]
