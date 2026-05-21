"""
生成专利图2：偏心环空几何示意图

参考 SnowShot 样式，包含三个子图：
  (a) 三维偏心环空透视示意图
  (b) 二维横截面几何定义（含标注）
  (c) 环空间隙厚度 H(φ) 三维曲面图
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc, FancyArrowPatch, Circle
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# ── 字体设置 ──────────────────────────────────────────────────
plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def draw_eccentric_annulus_3d(ax: Axes3D) -> None:
    """绘制三维偏心环空透视示意图（左图）。"""
    # 参数
    xi_max = 4.0
    n_xi = 3
    n_phi = 60
    phi = np.linspace(0, 2 * np.pi, n_phi)
    xi = np.linspace(0, xi_max, n_xi)

    r_o = 1.0          # 外筒半径（井壁）
    r_i = 0.55         # 内筒半径（套管）
    delta_r = 0.22     # 偏心距
    casing_center = (delta_r, 0.0)  # 套管中心偏右

    # 外筒表面
    Phi, Xi = np.meshgrid(phi, xi)
    X_o = r_o * np.cos(Phi)
    Y_o = r_o * np.sin(Phi)
    Z_o = Xi

    ax.plot_surface(X_o, Y_o, Z_o, alpha=0.18, color="steelblue", linewidth=0, antialiased=True)

    # 内筒（套管）表面
    X_i = casing_center[0] + r_i * np.cos(Phi)
    Y_i = casing_center[1] + r_i * np.sin(Phi)
    Z_i = Xi
    ax.plot_surface(X_i, Y_i, Z_i, alpha=0.75, color="indianred", linewidth=0, antialiased=True)

    # 顶部和底部的圆形边界
    for z_val in [0, xi_max]:
        # 外圆
        ax.plot(r_o * np.cos(phi), r_o * np.sin(phi), z_val,
                color="steelblue", linewidth=1.2)
        # 内圆
        ax.plot(casing_center[0] + r_i * np.cos(phi),
                casing_center[1] + r_i * np.sin(phi), z_val,
                color="indianred", linewidth=1.2)

    # ξ 轴箭头
    ax.quiver(-1.3, -1.3, 0, 0, 0, xi_max + 0.6,
              color="black", arrow_length_ratio=0.06, linewidth=1.5)
    ax.text(-1.5, -1.5, xi_max + 0.8, r"$\xi$", fontsize=14, ha="center", va="center")

    # 底部 φ=0 和 φ=1 标注
    # φ=0 在宽边（右边）
    ax.text(r_o + 0.25, 0.0, -0.3, r"$\phi=0$", fontsize=12, ha="left", va="center")
    # φ=1 在窄边（左边）
    ax.text(-r_o - 0.25, 0.0, -0.3, r"$\phi=1$", fontsize=12, ha="right", va="center")

    # 视图
    ax.view_init(elev=22, azim=-55)
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.4, 1.4)
    ax.set_zlim(-0.5, xi_max + 0.8)
    ax.axis("off")


def draw_cross_section(ax: plt.Axes) -> None:
    """绘制二维横截面几何定义图（右上），匹配参考图片样式。"""
    # 参数
    r_o = 1.0
    r_i = 0.55
    delta_r = 0.22
    outer_center = (0.0, 0.0)
    inner_center = (delta_r, 0.0)
    e = delta_r / (r_o - r_i)

    # 外圆
    outer_circle = Circle(outer_center, r_o, fill=False, edgecolor="steelblue", linewidth=2.2)
    ax.add_patch(outer_circle)

    # 内圆（灰色填充）
    inner_circle = Circle(inner_center, r_i, fill=True, facecolor="lightgray",
                          edgecolor="dimgray", linewidth=2.0)
    ax.add_patch(inner_circle)

    # 外圆中心标注
    ax.plot(*outer_center, "ko", markersize=4)
    ax.text(outer_center[0] - 0.12, outer_center[1] - 0.12, r"$O$", fontsize=12,
            ha="center", va="center")

    # 内圆中心标注
    ax.plot(*inner_center, "ko", markersize=4)
    ax.text(inner_center[0] + 0.14, inner_center[1] + 0.14, r"$O'$", fontsize=12,
            ha="center", va="center")

    # r_i 箭头（从内圆中心到内圆边缘）
    ax.annotate("", xy=(inner_center[0] + r_i * 0.85, inner_center[1] + r_i * 0.53),
                xytext=inner_center, fontsize=11,
                arrowprops=dict(arrowstyle="->", color="black", lw=1.2))
    ax.text(inner_center[0] + r_i * 0.45, inner_center[1] + r_i * 0.35,
            r"$r_i$", fontsize=13, ha="center", va="center")

    # r_o 箭头（从中心线到外圆）
    ax.annotate("", xy=(r_o * 0.88, 0.12),
                xytext=(inner_center[0], 0.0), fontsize=11,
                arrowprops=dict(arrowstyle="->", color="steelblue", lw=1.2))
    ax.text(r_o * 0.55, 0.25, r"$r_o$", fontsize=13, ha="center", va="center",
            color="steelblue")

    # Δr 标注 — 两个中心之间的水平距离
    # 参考线
    ax.plot([outer_center[0], outer_center[0]], [-0.85, -0.75], "k-", linewidth=0.8)
    ax.plot([inner_center[0], inner_center[0]], [-0.85, -0.75], "k-", linewidth=0.8)
    ax.annotate("", xy=(inner_center[0], -0.80), xytext=(outer_center[0], -0.80),
                arrowprops=dict(arrowstyle="<->", color="black", lw=1.2))
    ax.text(delta_r / 2, -0.92, r"$\Delta r$", fontsize=13, ha="center", va="top")

    # 偏心度公式 — 右上角
    ax.text(0.72, 0.95,
            r"$e = \frac{\Delta r}{r_o - r_i}$",
            transform=ax.transAxes, fontsize=14,
            ha="center", va="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.9))

    # H(ϕ) — 间隙宽度标注（左上区域）
    phi_angle = np.deg2rad(135)  # ≈135° 位置
    # 外圆上的点
    x_outer = r_o * np.cos(phi_angle)
    y_outer = r_o * np.sin(phi_angle)
    # 内圆上的点 (近似同方向)
    x_inner = inner_center[0] + r_i * np.cos(phi_angle)
    y_inner = inner_center[1] + r_i * np.sin(phi_angle)

    ax.annotate("", xy=(x_outer, y_outer), xytext=(x_inner, y_inner),
                arrowprops=dict(arrowstyle="<->", color="darkgreen", lw=1.2))
    mid_x = (x_outer + x_inner) / 2 - 0.15
    mid_y = (y_outer + y_inner) / 2 + 0.1
    ax.text(mid_x, mid_y, r"$H(\phi)$", fontsize=13, ha="center", va="center",
            color="darkgreen")

    # πϕ 角度标注
    # 从水平线逆时针到某个角度
    angle_arc = np.deg2rad(60)
    arc_radius = 0.35
    theta = np.linspace(0, angle_arc, 40)
    ax.plot(arc_radius * np.cos(theta), arc_radius * np.sin(theta), "k-", linewidth=1.0)
    ax.text(arc_radius * np.cos(angle_arc / 2) + 0.12,
            arc_radius * np.sin(angle_arc / 2) + 0.06,
            r"$\pi\phi$", fontsize=13, ha="center", va="center")

    ax.set_aspect("equal")
    ax.set_xlim(-1.25, 1.35)
    ax.set_ylim(-1.15, 1.25)
    ax.axis("off")


def draw_gap_surface(ax: Axes3D) -> None:
    """绘制环空间隙厚度 H(ϕ) 三维曲面图（右下）。"""
    # 几何参数（与上面一致）
    r_o = 1.0
    r_i = 0.55
    delta_r = 0.22

    # H(ϕ) = sqrt(r_o² - Δr² sin²ϕ) - r_i - Δr cosϕ (近似)
    # 使用精确公式
    phi_vals = np.linspace(0, 1, 50)  # 归一化 ϕ
    phi_rad = phi_vals * np.pi
    xi_vals = np.linspace(0, 1, 10)

    # 间隙宽度 H(ϕ)：外圆半径减内圆半径（沿径向）
    # H(ϕ) = sqrt(r_o² - Δr²·sin²(πϕ)) - r_i - Δr·cos(πϕ)
    H = np.sqrt(r_o**2 - delta_r**2 * np.sin(phi_rad)**2) - r_i - delta_r * np.cos(phi_rad)

    Phi_grid, Xi_grid = np.meshgrid(phi_vals, xi_vals)
    H_grid = np.tile(H, (len(xi_vals), 1))

    ax.plot_surface(Phi_grid, H_grid, Xi_grid, alpha=0.7, color="lightcoral",
                    edgecolor="darkred", linewidth=0.3, antialiased=True)

    # 坐标轴
    ax.set_xlabel(r"$\phi$", fontsize=13, labelpad=8)
    ax.set_ylabel(r"$H(\phi)$", fontsize=13, labelpad=8)
    # ξ 标签手动放置在纵轴左侧（3D zlabel 在部分视角下不渲染）
    ax.text2D(-0.02, 0.52, r"$\xi$", fontsize=13,
              transform=ax.transAxes, ha="right", va="center")

    # ϕ 轴刻度
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(["0", "0.5", "1"])

    # H 轴原点标注
    ax.text(0, 0, 0, "0", fontsize=11, ha="right", va="top")

    # 视角
    ax.view_init(elev=22, azim=-45)


def main() -> None:
    fig = plt.figure(figsize=(13, 9))

    # ── 布局 ────────────────────────────────────────────────
    # 左侧：3D 透视图（占左半部分）
    ax_3d = fig.add_axes([0.02, 0.08, 0.46, 0.88], projection="3d")
    draw_eccentric_annulus_3d(ax_3d)

    # 右上：2D 横截面
    ax_cross = fig.add_axes([0.52, 0.52, 0.46, 0.44])
    draw_cross_section(ax_cross)

    # 右下：H(ϕ) 曲面
    ax_surf = fig.add_axes([0.52, 0.06, 0.46, 0.42], projection="3d")
    draw_gap_surface(ax_surf)

    # ── 子图编号 ─────────────────────────────────────────────
    fig.text(0.02, 0.97, "(a) 三维偏心环空示意图", fontsize=13, fontweight="bold")
    fig.text(0.52, 0.97, "(b) 横截面参数定义", fontsize=13, fontweight="bold")
    fig.text(0.52, 0.49, "(c) 环空间隙 $H(\\phi)$ 曲面", fontsize=13, fontweight="bold")

    # ── 保存 ─────────────────────────────────────────────────
    output_path = (
        r"D:\users\desktop\固井顶替效率专利文件\专利撰写\figures"
        r"\fig2_eccentric_annulus.png"
    )
    fig.savefig(output_path, dpi=200, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"图2已保存至: {output_path}")


if __name__ == "__main__":
    main()
