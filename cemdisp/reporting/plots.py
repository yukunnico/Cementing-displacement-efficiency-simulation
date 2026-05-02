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
    fig.savefig(output_path / filename, dpi=150, bbox_inches="tight")


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
    - 上图：全井段、CBL评价井段、目标层段的有效顶替效率及水泥浆占据率随时间变化
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
            "cbl_eval_interval_efficiency",
            "target_interval_efficiency",
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

    axes[0].plot(metrics["time_min"], metrics["effective_efficiency"], label="全井段", linewidth=2.2)
    axes[0].plot(metrics["time_min"], metrics["cbl_eval_interval_efficiency"], label="CBL评价井段", linewidth=2)
    axes[0].plot(metrics["time_min"], metrics["target_interval_efficiency"], label="目标层段", linewidth=2)
    axes[0].plot(metrics["time_min"], metrics["bulk_cement_fill"], label="水泥浆占据率", linestyle="--", linewidth=2)
    axes[0].set_ylabel("有效顶替效率 / 占据率")
    axes[0].set_ylim(0, 1.05)
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(loc="best")

    axes[1].plot(metrics["time_min"], metrics["front_wide_m"], label="宽边", linewidth=2.1)
    axes[1].plot(metrics["time_min"], metrics["front_mid_m"], label="中线", linewidth=2.1)
    axes[1].plot(metrics["time_min"], metrics["front_narrow_m"], label="窄边", linewidth=2.1)
    axes[1].set_xlabel("时间 / min")
    axes[1].set_ylabel("前沿推进距离 / m")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(loc="best")

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
    axes[0].plot(depth, profiles["宽边有效效率"], label="宽边有效效率", linewidth=2)
    axes[0].plot(depth, profiles["中线有效效率"], label="中线有效效率", linewidth=2)
    axes[0].plot(depth, profiles["窄边有效效率"], label="窄边有效效率", linewidth=2)
    axes[0].plot(depth, profiles["平均有效顶替效率"], label="平均有效顶替效率", color="black", linewidth=2.4)
    axes[0].set_ylabel("有效顶替效率")
    axes[0].set_ylim(0, 1.05)
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(loc="best", ncol=2)

    axes[1].plot(depth, profiles["宽边水泥浓度"], label="宽边水泥浓度", linewidth=2)
    axes[1].plot(depth, profiles["中线水泥浓度"], label="中线水泥浓度", linewidth=2)
    axes[1].plot(depth, profiles["窄边水泥浓度"], label="窄边水泥浓度", linewidth=2)
    axes[1].plot(depth, profiles["水泥平均浓度"], label="水泥平均浓度", color="black", linewidth=2.4)
    axes[1].set_xlabel("井深 / m")
    axes[1].set_ylabel("水泥浓度")
    axes[1].set_ylim(0, 1.05)
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(loc="best", ncol=2)

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
    ax.plot(metrics["time_min"], metrics["channeling_index"], label="窜槽指数", linewidth=2.2)
    ax.plot(metrics["time_min"], metrics["mixing_index"], label="混浆指数", linewidth=2.2)
    ax.plot(metrics["time_min"], metrics["instability_index"], label="失稳指数", linewidth=2.2)
    ax.set_title(f"{result.well_name} 风险指标时间演变", fontsize=14, fontweight="bold")
    ax.set_xlabel("时间 / min")
    ax.set_ylabel("指数")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")

    fig.tight_layout()
    _save_figure(fig, output_dir, f"{well_name}_风险指标时间演变.png")
    return fig


def plot_efficiency_summary_bar(
    result: AnnulusSimulationResult,
    output_dir: Optional[Path | str] = None,
) -> Figure:
    """绘制最终有效顶替效率、水泥浆占据率与质量响应效率对比柱状图。

    绘制五个最终指标的柱状对比图：
    - 全井段有效顶替效率
    - CBL评价井段有效顶替效率
    - 目标层段有效顶替效率
    - 水泥浆占据率（水泥体积占环空体积的比例）
    - 质量响应效率（基于CBL代理值计算，与有效顶替效率不同）

    图表底部注释说明：质量响应效率 ≠ 有效顶替效率

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
            "cbl_eval_interval_efficiency",
            "target_interval_efficiency",
            "bulk_cement_fill",
            "cbl_quality_proxy",
        ],
        "metrics",
    )
    final_row = metrics.iloc[-1]
    labels = [
        "全井段有效顶替效率",
        "CBL评价井段有效顶替效率",
        "目标层段有效顶替效率",
        "水泥浆占据率",
        "质量响应效率",
    ]
    values = [
        final_row["effective_efficiency"],
        final_row["cbl_eval_interval_efficiency"],
        final_row["target_interval_efficiency"],
        final_row["bulk_cement_fill"],
        final_row["cbl_quality_proxy"],
    ]
    colors = ["#2f80ed", "#2f80ed", "#2f80ed", "#27ae60", "#f2994a"]

    well_name = _safe_filename_component(result.well_name)
    fig, ax = plt.subplots(figsize=(11, 6))
    bars = ax.bar(labels, values, color=colors, edgecolor="#333333", linewidth=0.6)
    ax.set_title(f"{result.well_name} 最终结果对比", fontsize=14, fontweight="bold")
    ax.set_ylabel("效率 / 占据率")
    ax.set_ylim(0, 1.08)
    ax.grid(axis="y", alpha=0.25)
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

    ax.text(
        0.01,
        -0.22,
        "* 质量响应效率 ≠ 有效顶替效率，基于 CBL 代理值计算",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        color="#8a4b08",
    )

    fig.tight_layout()
    _save_figure(fig, output_dir, f"{well_name}_最终结果对比.png")
    return fig


_setup_chinese_font()


__all__ = [
    "plot_depth_profiles",
    "plot_efficiency_summary_bar",
    "plot_risk_indices",
    "plot_time_series",
]
