"""
D2DGA消融实验论文图表（图6-9）

本模块为改进D2DGA三闭包消融论文产出：
- 图6: R0->R3 窜槽/混浆指数对比柱状图
- 图7: I3浮力弥散函数分析曲线（解析形态，替代后验q_phi重建）
- 图8: 真体力场形状与浮力数b（解析形态，替代后验场重建）
- 图9: R0->R3 顶替效率演进 + 理论偏心放大对照

所有图表均复用 cemdisp/reporting/plots.py 的中文字体配置和学术配色方案。
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from cemdisp.models2d.d2dga_flux import d2dga_dispersion_function_I3
from cemdisp.reporting.plots import (
    ACADEMIC_COLORS,
    STANDARD_FIGSIZE,
    _setup_chinese_font,
    _save_figure,
)

LOGGER = logging.getLogger(__name__)

# --- 复用 plots.py 的中文字体配置 ---
_setup_chinese_font()

# --- 输出目录 ---
_DEFAULT_OUTPUT_DIR = Path("results/ht1_004_ablation/figures")


def _ensure_output_dir(output_dir: Path | str | None = None) -> Path:
    """确保输出目录存在并返回 Path 对象。"""
    p = Path(output_dir) if output_dir is not None else _DEFAULT_OUTPUT_DIR
    p.mkdir(parents=True, exist_ok=True)
    return p


# ===========================================================================
# 图6: R0->R3 窜槽/混浆指数对比柱状图
# ===========================================================================


def plot_fig6_ablation_channeling_mixing(
    csv_path: str | Path = "results/ht1_004_ablation/ablation_summary.csv",
    output_dir: str | Path | None = None,
) -> Figure:
    """绘制图6: D2DGA三闭包消融——窜槽/混浆指数对比。

    从 full ablation 行（eccentricity=0.17, nz=500）提取 R0/R1/R2/R3
    的 channeling_index 和 mixing_index，以分组柱状图展示各级闭包的风险变化。

    Args:
        csv_path: ablation_summary.csv 路径。
        output_dir: 输出目录，默认 results/ht1_004_ablation/figures/。

    Returns:
        matplotlib Figure 对象。
    """
    df = pd.read_csv(csv_path)
    # 筛选 full ablation 行
    full = df[(df["eccentricity"] == 0.17) & (df["nz"] == 500)].copy()
    full = full.sort_values("ablation_level")

    levels = full["ablation_level"].tolist()
    channeling = full["channeling_index"].tolist()
    mixing = full["mixing_index"].tolist()

    x = np.arange(len(levels))
    width = 0.35

    fig, ax = plt.subplots(figsize=STANDARD_FIGSIZE["single"])
    bars1 = ax.bar(
        x - width / 2, channeling, width,
        label="窜槽指数 (channeling_index)",
        color=ACADEMIC_COLORS["danger"], edgecolor="#333333", linewidth=0.5,
    )
    bars2 = ax.bar(
        x + width / 2, mixing, width,
        label="混浆指数 (mixing_index)",
        color=ACADEMIC_COLORS["warning"], edgecolor="#333333", linewidth=0.5,
    )

    ax.set_xlabel("消融级别", fontsize=12)
    ax.set_ylabel("指数值", fontsize=12)
    ax.set_title("D2DGA三闭包消融：窜槽/混浆指数对比", fontsize=14, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(levels)
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(axis="y", alpha=0.2, linestyle="--")

    # 数值标注
    for bar, val in zip(bars1, channeling):
        ax.text(
            bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.003,
            f"{val:.4f}", ha="center", va="bottom", fontsize=8, rotation=90,
        )
    for bar, val in zip(bars2, mixing):
        ax.text(
            bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.003,
            f"{val:.4f}", ha="center", va="bottom", fontsize=8, rotation=90,
        )

    fig.tight_layout()
    out_dir = _ensure_output_dir(output_dir)
    _save_figure(fig, out_dir, "图6_消融窜槽混浆指数对比.png")
    return fig


# ===========================================================================
# 图7: I3浮力弥散函数分析曲线（解析形态）
# ===========================================================================


def plot_fig7_i3_buoyancy_dispersion(
    output_dir: str | Path | None = None,
) -> Figure:
    """绘制图7: I3(ḉ, m) 浮力弥散函数曲线（解析形态）。

    绘制 I3(ḉ, m) = ḉ²(1-ḉ)³[4mḉ + 3(1-ḉ)] / {2m[mḉ³ + 1 - ḉ³]} 在
    m ∈ {0.5, 1.0, 2.0} 下的解析曲线，展示 D2DGA 浮力弥散通量对水泥
    浓度和黏度比的依赖关系。

    注：本图为解析形态，替代后验 q_phi 重建。后验 q_phi 需要从 R2 运行的
    完整场数据中提取 ḉ 切片、H、η2、f_phi 等，且需多步近似。解析曲线
    更清晰地展示了 I3 函数的物理本质。

    Args:
        output_dir: 输出目录。

    Returns:
        matplotlib Figure 对象。
    """
    fig, ax = plt.subplots(figsize=STANDARD_FIGSIZE["single"])

    c_bar = np.linspace(0.0, 1.0, 200)
    m_values = [0.5, 1.0, 2.0]
    colors = [ACADEMIC_COLORS["primary"], ACADEMIC_COLORS["warning"], ACADEMIC_COLORS["danger"]]
    styles = ["-", "--", "-."]

    for m, color, style in zip(m_values, colors, styles):
        i3 = d2dga_dispersion_function_I3(c_bar, m)
        ax.plot(c_bar, i3, color=color, linestyle=style, linewidth=2.0,
                label=f"m = {m}")

    ax.set_xlabel("间隙平均水泥浓度 c", fontsize=12)
    ax.set_ylabel("I3(c, m)", fontsize=12)
    ax.set_title("D2DGA浮力弥散函数 I3(c, m) 解析形态", fontsize=14, fontweight="bold")
    ax.legend(loc="upper right", fontsize=10)
    ax.grid(True, alpha=0.2, linestyle="--")
    ax.set_xlim(0, 1)
    ax.set_ylim(bottom=0)

    # 注释：c≈0.5附近达峰
    ax.annotate(
        "c≈0.5 附近达峰\n（最大浮力弥散效应）",
        xy=(0.48, float(d2dga_dispersion_function_I3(0.48, 1.0))),
        xytext=(0.25, 0.85 * float(d2dga_dispersion_function_I3(0.48, 1.0))),
        arrowprops=dict(arrowstyle="->", color="gray", alpha=0.7),
        fontsize=9, color="gray",
    )

    fig.tight_layout()
    out_dir = _ensure_output_dir(output_dir)
    _save_figure(fig, out_dir, "图7_I3浮力弥散函数.png")
    return fig


# ===========================================================================
# 图8: 真体力场形状与浮力数 b
# ===========================================================================


def plot_fig8_true_buoyancy_shape(
    output_dir: str | Path | None = None,
    buoyancy_number: float = 4.02,
    eccentricity: float = 0.17,
) -> Figure:
    """绘制图8: 真体力场形状与浮力数 b。

    绘制 R3 真体力 (2φ-1)·ebar 的归一化方位角形状剖面，并标注浮力数 b。
    注：本图为解析形态，替代后验场重建。后验场重建需要从 R3 运行的完整
    geom + 密度场中提取，涉及多层近似。解析形态直接展示了 R3 真体力的
    数学结构（式2.5b）。

    Args:
        output_dir: 输出目录。
        buoyancy_number: 浮力数 b（来自 full ablation R3）。
        eccentricity: 偏心度 e。

    Returns:
        matplotlib Figure 对象。
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=STANDARD_FIGSIZE["double"])

    # --- 左图: 真体力方位角形状 ---
    phi = np.linspace(0, 1, 100)  # 归一化方位角: 0=宽边, 1=窄边
    # (2φ-1) 形状因子: 在 φ=0(宽边) 为负(-1), φ=1(窄边) 为正(+1)
    shape = 2.0 * phi - 1.0

    ax1.plot(phi, shape, color=ACADEMIC_COLORS["primary"], linewidth=2.5)
    ax1.fill_between(phi, 0, shape, color=ACADEMIC_COLORS["primary"], alpha=0.15)
    ax1.axhline(0, color="black", linewidth=0.5, linestyle="--")
    ax1.set_xlabel("归一化方位角 φ", fontsize=11)
    ax1.set_ylabel("(2φ-1) 形状因子", fontsize=11)
    ax1.set_title("真体力方位角形状 (2φ-1)", fontsize=12, fontweight="bold")
    ax1.grid(True, alpha=0.2, linestyle="--")
    ax1.annotate(
        "宽边\n(φ=0, 轻侧)", xy=(0.02, -0.95),
        fontsize=9, color="gray", ha="left",
    )
    ax1.annotate(
        "窄边\n(φ=1, 重侧)", xy=(0.98, 0.95),
        fontsize=9, color="gray", ha="right",
    )

    # --- 右图: 浮力数 b 标注 + ebar 信息 ---
    ax2.axis("off")
    text_lines = [
        f"R3 真体力参数 (HT1-004)",
        "",
        f"浮力数 b = {buoyancy_number:.2f}",
        f"（b > 0: 密度稳定，重顶替轻）",
        "",
        f"偏心度 e = {eccentricity:.2f}",
        f"（e·cos(πφ) 调制间隙场）",
        "",
        "R3 真体力形式：",
        "f_body ∝ e · (2φ-1) · Δρ/η",
        "",
        "注：本图为解析形态，替代后验场重建。",
        "后验重建需从 R3 全场数据逐层提取",
        "密度场、几何场和体力向量，近似较多。",
        "解析形态已充分展示 R3 真体力的",
        "数学结构（式 2.5b）。",
    ]
    for i, line in enumerate(text_lines):
        y_pos = 0.95 - i * 0.055
        if i == 0:
            ax2.text(0.05, y_pos, line, transform=ax2.transAxes,
                     fontsize=13, fontweight="bold", va="top")
        elif line.startswith("注："):
            ax2.text(0.05, y_pos, line, transform=ax2.transAxes,
                     fontsize=8, color="gray", va="top", style="italic")
        else:
            ax2.text(0.05, y_pos, line, transform=ax2.transAxes,
                     fontsize=10, va="top")

    fig.suptitle("R3 真体力场形状与浮力数", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out_dir = _ensure_output_dir(output_dir)
    _save_figure(fig, out_dir, "图8_真体力场形状与浮力数.png")
    return fig


# ===========================================================================
# 图9: R0->R3 顶替效率演进 + 理论偏心放大对照
# ===========================================================================


def plot_fig9_efficiency_evolution(
    csv_path: str | Path = "results/ht1_004_ablation/ablation_summary.csv",
    output_dir: str | Path | None = None,
) -> Figure:
    """绘制图9: 顶替效率R0->R3演进（现场17% vs 理论35%/45%偏心）。

    图9展示三个子图/面板：
    1. 现场偏心 e=17%: R0->R3 效率演进（主图折线）
    2. 理论偏心 e=35%/45%: R0 vs R3 效率差异（散点/标记）
    3. 网格收敛性: nz=140/280/500 的效率变化

    Args:
        csv_path: ablation_summary.csv 路径。
        output_dir: 输出目录。

    Returns:
        matplotlib Figure 对象。
    """
    df = pd.read_csv(csv_path)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))

    # --- 左图: 现场偏心 e=17% 的 R0->R3 折线 ---
    ax = axes[0]
    full = df[(df["eccentricity"] == 0.17) & (df["nz"] == 500)]
    full = full.sort_values("ablation_level")
    ax.plot(
        full["ablation_level"], full["effective_efficiency"],
        marker="o", color=ACADEMIC_COLORS["primary"], linewidth=2.5, markersize=8,
        label="现场 e=17%, nz=500",
    )
    for _, row in full.iterrows():
        ax.annotate(
            f"{row['effective_efficiency']:.4f}",
            (row["ablation_level"], row["effective_efficiency"]),
            textcoords="offset points", xytext=(0, 10),
            fontsize=8, ha="center",
        )
    ax.set_xlabel("消融级别", fontsize=11)
    ax.set_ylabel("有效顶替效率", fontsize=11)
    ax.set_title("现场 e=17% R0→R3 效率演进", fontsize=11, fontweight="bold")
    ax.grid(True, alpha=0.2, linestyle="--")
    ax.set_ylim(0.94, 0.97)

    # --- 中图: 理论偏心 e=35%/45% 的 R0 vs R3 散点 ---
    ax = axes[1]
    theoretical = df[df["eccentricity"].isin([0.35, 0.45])].copy()
    theoretical = theoretical.sort_values(["eccentricity", "ablation_level"])

    e_colors = {0.35: ACADEMIC_COLORS["warning"], 0.45: ACADEMIC_COLORS["danger"]}
    e_markers = {0.35: "s", 0.45: "^"}

    for e_val in [0.35, 0.45]:
        subset = theoretical[theoretical["eccentricity"] == e_val]
        ax.plot(
            subset["ablation_level"], subset["effective_efficiency"],
            marker=e_markers[e_val], color=e_colors[e_val],
            linewidth=2.0, markersize=9,
            label=f"e={int(e_val*100)}%",
        )
        for _, row in subset.iterrows():
            ax.annotate(
                f"{row['effective_efficiency']:.4f}",
                (row["ablation_level"], row["effective_efficiency"]),
                textcoords="offset points", xytext=(0, -14),
                fontsize=7.5, ha="center",
            )

    # 标注 R3-R0 差距
    for e_val in [0.35, 0.45]:
        subset = theoretical[theoretical["eccentricity"] == e_val]
        r0 = subset[subset["ablation_level"] == "R0"]["effective_efficiency"].values[0]
        r3 = subset[subset["ablation_level"] == "R3"]["effective_efficiency"].values[0]
        ax.annotate(
            f"Δ={r0-r3:.4f}",
            xy=(1, r3), xytext=(1.15, (r0 + r3) / 2),
            arrowprops=dict(arrowstyle="->", color=e_colors[e_val], alpha=0.6),
            fontsize=8, color=e_colors[e_val],
        )

    ax.set_xlabel("消融级别", fontsize=11)
    ax.set_ylabel("有效顶替效率", fontsize=11)
    ax.set_title("理论偏心放大：e=35%/45%", fontsize=11, fontweight="bold")
    ax.legend(loc="lower left", fontsize=8)
    ax.grid(True, alpha=0.2, linestyle="--")

    # --- 右图: 网格收敛性 ---
    ax = axes[2]
    conv = df[df["run_id"].str.startswith("convergence")].copy()
    # nz convergence
    nz_conv = conv[conv["dt"] == 4.0].sort_values("nz")
    ax.plot(
        nz_conv["nz"], nz_conv["effective_efficiency"],
        marker="o", color=ACADEMIC_COLORS["primary"], linewidth=2.0, markersize=8,
        label="nz 收敛 (dt=4s)",
    )
    for _, row in nz_conv.iterrows():
        ax.annotate(
            f"{row['effective_efficiency']:.4f}",
            (row["nz"], row["effective_efficiency"]),
            textcoords="offset points", xytext=(0, 8),
            fontsize=8, ha="center",
        )

    # dt convergence
    dt_conv = conv[conv["nz"] == 280].sort_values("dt")
    ax_twin = ax.twinx()
    # Plot dt convergence on same x-axis but different y scale
    # Actually, let's use a different approach - plot dt on secondary axis
    ax2 = ax.twiny()
    ax2.plot(
        dt_conv["dt"], dt_conv["effective_efficiency"],
        marker="s", color=ACADEMIC_COLORS["danger"], linewidth=2.0, markersize=8,
        linestyle="--", label="dt 收敛 (nz=280)",
    )
    ax2.set_xlabel("dt / s", fontsize=10, color=ACADEMIC_COLORS["danger"])
    ax2.tick_params(axis="x", colors=ACADEMIC_COLORS["danger"])

    ax.set_xlabel("nz", fontsize=11)
    ax.set_ylabel("有效顶替效率", fontsize=11)
    ax.set_title("网格/时间步收敛 (R3)", fontsize=11, fontweight="bold")
    ax.grid(True, alpha=0.2, linestyle="--")

    # Combine legends
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc="lower right", fontsize=8)

    fig.suptitle("顶替效率R0→R3演进（现场17% vs 理论35%/45%偏心）", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out_dir = _ensure_output_dir(output_dir)
    _save_figure(fig, out_dir, "图9_效率演进与理论偏心放大.png")
    return fig


# ===========================================================================
# 批量生成所有图表
# ===========================================================================


def generate_all_ablation_figures(
    csv_path: str | Path = "results/ht1_004_ablation/ablation_summary.csv",
    output_dir: str | Path | None = None,
    buoyancy_number: float = 4.02,
    eccentricity: float = 0.17,
) -> dict[str, Figure]:
    """批量生成图6-9全部四张消融论文图表。

    Args:
        csv_path: ablation_summary.csv 路径。
        output_dir: 输出目录，默认 results/ht1_004_ablation/figures/。
        buoyancy_number: R3 浮力数 b（用于图8）。
        eccentricity: 现场偏心度 e（用于图8）。

    Returns:
        {"fig6": Figure, "fig7": Figure, "fig8": Figure, "fig9": Figure}
    """
    figures: dict[str, Figure] = {}
    LOGGER.info("生成图6: 窜槽/混浆指数对比...")
    figures["fig6"] = plot_fig6_ablation_channeling_mixing(csv_path, output_dir)

    LOGGER.info("生成图7: I3浮力弥散函数...")
    figures["fig7"] = plot_fig7_i3_buoyancy_dispersion(output_dir)

    LOGGER.info("生成图8: 真体力场形状...")
    figures["fig8"] = plot_fig8_true_buoyancy_shape(
        output_dir, buoyancy_number=buoyancy_number, eccentricity=eccentricity,
    )

    LOGGER.info("生成图9: 效率演进与理论偏心放大...")
    figures["fig9"] = plot_fig9_efficiency_evolution(csv_path, output_dir)

    LOGGER.info("图6-9全部生成完毕。")
    return figures


__all__ = [
    "plot_fig6_ablation_channeling_mixing",
    "plot_fig7_i3_buoyancy_dispersion",
    "plot_fig8_true_buoyancy_shape",
    "plot_fig9_efficiency_evolution",
    "generate_all_ablation_figures",
]