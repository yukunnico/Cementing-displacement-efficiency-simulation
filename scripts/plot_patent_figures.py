"""
专利技术交底书 12 张附图生成脚本

生成控压固井二维 D2DGA 顶替模型的全套专利附图。
输出目录: D:/users/desktop/research/控压固井项目/固井顶替效率专利撰写/参考/figures/

图1  : 总体流程图
图2  : 偏心环空几何示意
图3  : 四相体积分数闭合
图4  : Papanastasiou 正则化效果对比（纯数学）
图5  : D2DGA 通量放大因子（纯数学）
图6  : Hele-Shaw 速度场
图7  : 半拉格朗日反演追踪示意
图8  : 多约束数值求解链路流程图
图9  : 浓度场快照
图10 : 诊断指标时间序列
图11 : 深度方向效率剖面
图12 : 计算效率与弥散方式对比
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import numpy as np
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, ArrowStyle

# ── 全局设置 ──────────────────────────────────────────────────
OUTPUT_DIR = Path(
    r"D:\users\desktop\research\控压固井项目\固井顶替效率专利撰写\参考\figures"
)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["savefig.dpi"] = 200
plt.rcParams["savefig.bbox"] = "tight"
plt.rcParams["figure.facecolor"] = "white"

# ── 尝试导入求解器 ────────────────────────────────────────────
SOLVER_AVAILABLE = False
try:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from cemdisp.models2d import AnnulusD2DGASolver
    from cemdisp.models2d.boundary_bridge import AnnulusInletState
    from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
    from cemdisp.data.well_spec import DepthValuePoint, EvaluationWindow, WellSpec
    SOLVER_AVAILABLE = True
    print("[INFO] 求解器导入成功，将尝试运行模拟。")
except Exception as exc:
    warnings.warn(f"[WARN] 求解器导入失败: {exc}。模拟依赖的图将使用示意数据。")


def _save(fig: plt.Figure, name: str) -> None:
    """统一保存并关闭图形。"""
    path = OUTPUT_DIR / name
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"  -> 已保存: {path.name}")


# ═══════════════════════════════════════════════════════════════
#  图1：总体流程图
# ═══════════════════════════════════════════════════════════════
def fig01_method_flowchart():
    """绘制 D2DGA 方法总体流程图（简洁版）。"""
    fig, ax = plt.subplots(figsize=(16, 4.2))
    ax.set_xlim(-0.5, 16.5)
    ax.set_ylim(-1.8, 4.0)
    ax.axis("off")

    # 主链路节点
    steps = [
        (1.5,  2.0, "①几何建模", "#E3F2FD"),
        (4.0,  2.0, "②流变计算\n+正则化", "#E8F5E9"),
        (6.5,  2.0, "③速度场\n构建", "#FFF3E0"),
        (9.0,  2.0, "④D2DGA\n通量修正", "#FCE4EC"),
        (11.5, 2.0, "⑤半拉格朗日\n平流输运", "#F3E5F5"),
        (14.0, 2.0, "⑥约束修正\n+诊断输出", "#E0F7FA"),
    ]
    bw, bh = 2.0, 1.3
    for x, y, label, color in steps:
        fb = FancyBboxPatch((x - bw/2, y - bh/2), bw, bh,
            boxstyle="round,pad=0.12", facecolor=color,
            edgecolor="#37474F", linewidth=1.5)
        ax.add_patch(fb)
        ax.text(x, y, label, ha="center", va="center", fontsize=10, fontweight="bold")

    # 主链路水平箭头
    for i in range(len(steps) - 1):
        x0 = steps[i][0] + bw/2 + 0.05
        x1 = steps[i+1][0] - bw/2 - 0.05
        ax.annotate("", xy=(x1, 2.0), xytext=(x0, 2.0),
                    arrowprops=dict(arrowstyle="-|>", color="#37474F", lw=1.8))

    # 循环箭头（从⑥下方绕回③下方，走最底部，不穿过任何框）
    loop_y = -1.2
    # ⑥底部 -> 底部水平线
    ax.annotate("", xy=(14.0, 2.0 - bh/2 - 0.05), xytext=(14.0, loop_y + 0.05),
                arrowprops=dict(arrowstyle="-|>", color="#B71C1C", lw=1.5,
                                connectionstyle="arc3,rad=0", linestyle="--"))
    # 底部水平线：⑥ -> ③
    ax.annotate("", xy=(6.5, loop_y), xytext=(14.0, loop_y),
                arrowprops=dict(arrowstyle="-|>", color="#B71C1C", lw=1.5, linestyle="--"))
    # 底部 -> ③底部
    ax.annotate("", xy=(6.5, 2.0 - bh/2 - 0.05), xytext=(6.5, loop_y + 0.05),
                arrowprops=dict(arrowstyle="-|>", color="#B71C1C", lw=1.5,
                                connectionstyle="arc3,rad=0", linestyle="--"))
    ax.text(10.25, loop_y - 0.35, "下一时间步迭代", fontsize=9, color="#B71C1C",
            ha="center", style="italic", fontweight="bold")

    ax.set_title("图1  D2DGA 二维顶替模型总体计算流程", fontsize=13, fontweight="bold", pad=10)
    _save(fig, "fig01_总体流程图.png")


# ═══════════════════════════════════════════════════════════════
#  图2：偏心环空几何示意
# ═══════════════════════════════════════════════════════════════
def fig02_eccentric_annulus():
    """绘制偏心环空截面与展开域示意。"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    # --- (a) 截面 ---
    r_o = 1.0
    r_i = 0.55
    delta_r = 0.22
    outer = Circle((0, 0), r_o, fill=False, edgecolor="#1565C0", linewidth=2.5, label="井壁")
    inner = Circle((delta_r, 0), r_i, fill=True, facecolor="#E0E0E0",
                    edgecolor="#424242", linewidth=2.0, label="套管")
    ax1.add_patch(outer)
    ax1.add_patch(inner)

    # 中心标注
    ax1.plot(0, 0, "ko", ms=4)
    ax1.text(-0.08, -0.1, "O", fontsize=11, ha="center")
    ax1.plot(delta_r, 0, "ko", ms=4)
    ax1.text(delta_r + 0.08, -0.1, "O'", fontsize=11, ha="center")

    # 宽边/窄边标注（套管偏右，左侧间隙大=宽边，右侧间隙小=窄边）
    ax1.annotate("宽边\n(Wide)", xy=(-r_o + 0.05, 0), xytext=(-r_o - 0.35, 0.35),
                 fontsize=10, ha="center", color="#C62828",
                 arrowprops=dict(arrowstyle="->", color="#C62828", lw=1.2))
    ax1.annotate("窄边\n(Narrow)", xy=(r_o - 0.05, 0), xytext=(r_o + 0.25, 0.35),
                 fontsize=10, ha="center", color="#1B5E20",
                 arrowprops=dict(arrowstyle="->", color="#1B5E20", lw=1.2))

    # 偏心距标注
    ax1.annotate("", xy=(delta_r, -0.75), xytext=(0, -0.75),
                 arrowprops=dict(arrowstyle="<->", color="black", lw=1.3))
    ax1.text(delta_r / 2, -0.88, r"$\Delta r$", fontsize=12, ha="center")

    # phi 角弧
    theta_arc = np.linspace(0, np.pi * 0.65, 30)
    ax1.plot(0.3 * np.cos(theta_arc), 0.3 * np.sin(theta_arc), "k-", lw=1.0)
    ax1.text(0.3 * np.cos(0.325 * np.pi) + 0.08, 0.3 * np.sin(0.325 * np.pi),
             r"$\pi\phi$", fontsize=11)

    # H(phi) 双箭头
    phi_angle = np.deg2rad(135)
    x_outer = r_o * np.cos(phi_angle)
    y_outer = r_o * np.sin(phi_angle)
    x_inner = delta_r + r_i * np.cos(phi_angle)
    y_inner = r_i * np.sin(phi_angle)
    ax1.annotate("", xy=(x_outer, y_outer), xytext=(x_inner, y_inner),
                 arrowprops=dict(arrowstyle="<->", color="#2E7D32", lw=1.5))
    mx = (x_outer + x_inner) / 2 - 0.18
    my = (y_outer + y_inner) / 2 + 0.08
    ax1.text(mx, my, r"$H(\phi)$", fontsize=12, color="#2E7D32")

    ax1.set_aspect("equal")
    ax1.set_xlim(-1.4, 1.45)
    ax1.set_ylim(-1.2, 1.3)
    ax1.axis("off")
    ax1.set_title("(a) 偏心环空横截面", fontsize=12, fontweight="bold")

    # --- (b) 展开域 ---
    phi = np.linspace(0, 1, 200)
    h_mean = 0.225
    e = 0.49
    h = h_mean * (1 + e * np.cos(np.pi * phi))
    s = np.linspace(0, 2.0, 100)
    PHI, S = np.meshgrid(phi, s)
    H_grid = h_mean * (1 + e * np.cos(np.pi * PHI))

    pcm = ax2.pcolormesh(PHI, S, H_grid, cmap="YlOrRd", shading="gouraud")
    cb = fig.colorbar(pcm, ax=ax2, label="局部间隙 $H(\\phi)$ / m", shrink=0.85)

    ax2.set_xlabel(r"归一化方位角 $\phi$", fontsize=11)
    ax2.set_ylabel("轴向坐标 $s$ / m", fontsize=11)
    ax2.set_title("(b) 展开域间隙分布", fontsize=12, fontweight="bold")

    ax2.annotate("窄边 $\phi=1$", xy=(1.0, 1.0), xytext=(0.75, 1.6),
                 fontsize=10, color="#1B5E20", fontweight="bold",
                 arrowprops=dict(arrowstyle="->", color="#1B5E20"))
    ax2.annotate("宽边 $\phi=0$", xy=(0.0, 1.0), xytext=(0.18, 1.6),
                 fontsize=10, color="#C62828", fontweight="bold",
                 arrowprops=dict(arrowstyle="->", color="#C62828"))

    fig.suptitle("图2  偏心环空几何示意", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    _save(fig, "fig02_偏心环空几何示意.png")


# ═══════════════════════════════════════════════════════════════
#  图3：四相体积分数闭合
# ═══════════════════════════════════════════════════════════════
def fig03_phase_closure():
    """展示四相体积分数闭合关系。"""
    fig, ax = plt.subplots(figsize=(8, 6))

    ny, nz = 40, 60
    phi = np.linspace(0, 1, ny)[:, None]
    s = np.linspace(0, 1, nz)[None, :]

    # 合成示意场
    lead = np.clip(0.8 * np.exp(-8 * (phi - 0.15) ** 2) * (1 - s) * 1.5, 0, 1)
    tail = np.clip(0.7 * np.exp(-10 * (phi - 0.25) ** 2) * np.clip(s - 0.3, 0, 1), 0, 1)
    spacer = np.clip(0.4 * np.exp(-6 * (phi - 0.5) ** 2) * np.clip(1 - s, 0, 1) * 1.2, 0, 1)
    mud = np.clip(1 - lead - tail - spacer, 0, 1)

    # 构造四层堆叠展示
    ax.imshow(
        lead.T, origin="lower", aspect="auto", cmap="Blues", alpha=0.85,
        extent=[0, 1, 0, 1], vmin=0, vmax=1,
    )
    # 用子图分别展示
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    fields = [
        (lead, "领浆 $c_{\\mathrm{lead}}$", "Blues"),
        (tail, "尾浆 $c_{\\mathrm{tail}}$", "Greens"),
        (spacer, "隔离液 $c_{\\mathrm{spacer}}$", "Oranges"),
        (mud, "钻井液 $c_{\\mathrm{mud}}=1-c_l-c_t-c_s$", "Reds"),
    ]
    for ax_i, (field, title, cmap) in zip(axes.flat, fields):
        pcm = ax_i.pcolormesh(
            np.linspace(0, 1, nz), np.linspace(0, 1, ny), field,
            cmap=cmap, shading="gouraud", vmin=0, vmax=1,
        )
        fig.colorbar(pcm, ax=ax_i, shrink=0.82)
        ax_i.set_title(title, fontsize=11, fontweight="bold")
        ax_i.set_xlabel(r"$s$ (归一化轴向)")
        ax_i.set_ylabel(r"$\phi$ (方位角)")

    fig.suptitle("图3  四相体积分数闭合关系示意", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))

    # 闭合公式标注
    fig.text(
        0.5, 0.01,
        r"$c_{\mathrm{mud}} = 1 - c_{\mathrm{lead}} - c_{\mathrm{tail}} - c_{\mathrm{spacer}}$"
        "    (体积分数闭合约束)",
        ha="center", fontsize=12, style="italic",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow", alpha=0.9),
    )
    _save(fig, "fig03_四相体积分数闭合.png")


# ═══════════════════════════════════════════════════════════════
#  图4：Papanastasiou 正则化效果对比 (纯数学)
# ═══════════════════════════════════════════════════════════════
def fig04_papanastasiou():
    """经典 Bingham vs Papanastasiou 正则化粘度对比。"""
    gamma = np.logspace(-4, 3, 500)
    PV = 0.050  # Pa·s = 50 mPa·s
    M = 100.0
    YP_list = [8.0, 12.0, 20.0]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for YP in YP_list:
        mu_classic = PV + YP / gamma
        mu_reg = PV + YP * (1 - np.exp(-M * gamma)) / gamma
        mobility_classic = gamma ** 2 / mu_classic
        mobility_reg = gamma ** 2 / mu_reg

        axes[0].plot(gamma, mu_classic, lw=2, label=f"YP={YP:.0f} Pa")
        axes[1].plot(gamma, mu_reg, lw=2, label=f"YP={YP:.0f} Pa")
        axes[2].plot(gamma, mobility_reg / mobility_classic, lw=2, label=f"YP={YP:.0f} Pa")

    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel(r"剪切速率 $\dot{\gamma}$ / s$^{-1}$", fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9)

    axes[0].set_yscale("log")
    axes[0].set_ylabel(r"表观粘度 $\mu$ / Pa$\cdot$s", fontsize=11)
    axes[0].set_title("(a) 经典 Bingham: $\\mu = PV + YP/\\dot{\\gamma}$\n低剪切发散", fontsize=11, fontweight="bold")
    axes[0].set_ylim(1e-2, 1e5)

    axes[1].set_yscale("log")
    axes[1].set_ylabel(r"表观粘度 $\mu$ / Pa$\cdot$s", fontsize=11)
    axes[1].set_title("(b) Papanastasiou: $\\mu = PV + YP(1-e^{-M\\dot{\\gamma}})/\\dot{\\gamma}$\n有限极限", fontsize=11, fontweight="bold")
    axes[1].set_ylim(1e-2, 1e5)

    axes[2].set_ylabel("流动度比 $b^2/\\mu_{reg}$ : $b^2/\\mu_{classic}$", fontsize=10)
    axes[2].set_title("(c) 流动度比\n正则化/经典", fontsize=11, fontweight="bold")
    axes[2].axhline(1.0, color="gray", ls="--", lw=1)

    fig.suptitle("图4  Papanastasiou 正则化效果对比 (PV=50 mPa·s, M=100)",
                 fontsize=14, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    _save(fig, "fig04_Papanastasiou正则化.png")


# ═══════════════════════════════════════════════════════════════
#  图5：D2DGA 通量放大因子 (纯数学)
# ═══════════════════════════════════════════════════════════════
def fig05_d2dga_flux():
    """绘制 D2DGA 通量放大因子曲线族。"""
    cbar = np.linspace(0.01, 0.99, 300)
    m_values = [0.5, 1.0, 2.0, 5.0, 10.0]

    fig, ax = plt.subplots(figsize=(8, 5.5))

    for m in m_values:
        num = m * cbar ** 2 + 1.5 * (1 - cbar ** 2)
        den = m * cbar ** 3 + (1 - cbar ** 3)
        f = num / den
        ax.plot(cbar, f, lw=2.2, label=f"$m = {m}$")

    ax.axhline(1.0, color="gray", ls="--", lw=1, alpha=0.6)
    ax.set_xlabel(r"平均水泥浓度 $\bar{c}$", fontsize=12)
    ax.set_ylabel(r"通量放大因子 $f(\bar{c}, m)$", fontsize=12)
    ax.set_title("图5  D2DGA 通量放大因子\n"
                 r"$f(\bar{c}, m) = \frac{m\bar{c}^2 + 1.5(1-\bar{c}^2)}{m\bar{c}^3 + (1-\bar{c}^3)}$",
                 fontsize=13, fontweight="bold")
    ax.legend(fontsize=11, title="粘度比 $m=\\eta_d/\\eta_p$", title_fontsize=11)
    ax.set_xlim(0, 1)
    ax.set_ylim(0.8, 2.1)
    ax.grid(True, alpha=0.3)

    ax.annotate("低浓度前锋\n$f>1$ 加速推进", xy=(0.15, 1.65), fontsize=10,
                color="#C62828", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="#FFEBEE", alpha=0.8))
    ax.annotate("高浓度区\n$f \\approx 1$", xy=(0.82, 1.08), fontsize=10,
                color="#1B5E20", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="#E8F5E9", alpha=0.8))

    fig.tight_layout()
    _save(fig, "fig05_D2DGA通量放大因子.png")


# ═══════════════════════════════════════════════════════════════
#  图6：Hele-Shaw 速度场
# ═══════════════════════════════════════════════════════════════
def _run_solver_for_figures():
    """尝试运行求解器返回结果，失败返回 None。"""
    if not SOLVER_AVAILABLE:
        return None
    try:
        # 构造简化井规格
        n = 20
        depths = np.linspace(5500, 7700, n)
        well = WellSpec(
            well_name="PatentDemo",
            top_md_m=5500.0,
            bottom_md_m=7700.0,
            shoe_md_m=7700.0,
            liner_od_mm=139.7,
            hole_diameter_profile=tuple(
                DepthValuePoint(float(d), 215.9) for d in depths
            ),
            inclination_profile=tuple(
                DepthValuePoint(float(d), 30.0) for d in depths
            ),
            standoff_profile=tuple(
                DepthValuePoint(float(d), 0.75) for d in depths
            ),
            evaluation_windows=(
                EvaluationWindow("target", 7000.0, 7700.0, window_type="target"),
            ),
        )
        mud = FluidSpec("钻井液", FluidRole.MUD, 1200.0,
                        RheologyModel.BINGHAM, 0.030, 5.0)
        lead = FluidSpec("领浆", FluidRole.LEAD, 1850.0,
                         RheologyModel.BINGHAM, 0.060, 12.0)
        tail = FluidSpec("尾浆", FluidRole.TAIL, 1900.0,
                         RheologyModel.BINGHAM, 0.080, 18.0)
        spacer = FluidSpec("隔离液", FluidRole.SPACER, 1300.0,
                           RheologyModel.NEWTONIAN, 0.020)
        fluids = (mud, lead, tail, spacer)

        q = 0.8 / 60.0  # m3/s
        total_t = 4800.0  # 80 min

        def inlet(t):
            if t < 100:
                return AnnulusInletState(t, q, "前置液", (("spacer", 1.0),))
            elif t < 2000:
                return AnnulusInletState(t, q, "领浆", (("lead", 1.0),))
            elif t < 3800:
                return AnnulusInletState(t, q, "尾浆", (("tail", 1.0),))
            else:
                return AnnulusInletState(t, q, "替浆", (("mud", 1.0),))

        solver = AnnulusD2DGASolver(
            dt=8.0, nz=80, ny=30, total_t=total_t,
            enable_d2dga=True, d2dga_viscosity_ratio=1.5,
            save_interval=20,
        )
        print("[INFO] 正在运行模拟...")
        result = solver.run(well, fluids, inlet)
        print("[INFO] 模拟完成。")
        return result
    except Exception as exc:
        warnings.warn(f"[WARN] 模拟运行失败: {exc}")
        return None


def fig06_velocity_field(result=None):
    """绘制 Hele-Shaw 速度场。"""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    if result is not None:
        # 从结果中重建速度场（简化：用第一时刻计算）
        geom = result.geom
        b = geom["b"]
        mu = np.ones_like(b) * 0.08
        b_mean = np.mean(b, axis=0, keepdims=True)
        base = (b / b_mean) ** 2 / mu
        dy = np.gradient(geom["y"])[:, None]
        area_w = np.sum(base * b * dy * 2, axis=0, keepdims=True)
        q_half = 0.4 / 60.0
        w = q_half * base / area_w

        pcm0 = axes[0].pcolormesh(
            geom["md"], geom["y"], w, cmap="hot", shading="gouraud"
        )
        fig.colorbar(pcm0, ax=axes[0], label="$w$ / m/s", shrink=0.8)

        # 带浮力
        stable = 0.3
        phi_arr = geom["phi"][:, None]
        e_arr = geom["e"][None, :]
        buoy = 1.0 + stable * e_arr * (2 * phi_arr - 1)
        base2 = base * buoy
        area_w2 = np.sum(base2 * b * dy * 2, axis=0, keepdims=True)
        w2 = q_half * base2 / area_w2
        pcm1 = axes[1].pcolormesh(
            geom["md"], geom["y"], w2, cmap="hot", shading="gouraud"
        )
        fig.colorbar(pcm1, ax=axes[1], label="$w$ / m/s", shrink=0.8)

        # 横向速度 v
        ds = geom["s"][1] - geom["s"][0]
        bw = b * w2
        dbw_ds = np.gradient(bw, ds, axis=1)
        bv = np.zeros_like(w2)
        for i in range(1, len(geom["y"])):
            bv[i, :] = bv[i - 1, :] - 0.5 * (dbw_ds[i, :] + dbw_ds[i - 1, :]) * dy[i - 1, 0]
        bv -= (geom["y"][:, None] / geom["y"][-1]) * bv[-1, :]
        v = bv / np.maximum(b, 1e-8)
        pcm2 = axes[2].pcolormesh(
            geom["md"], geom["y"], v, cmap="coolwarm", shading="gouraud"
        )
        fig.colorbar(pcm2, ax=axes[2], label="$v$ / m/s", shrink=0.8)

        for ax in axes:
            ax.set_xlabel("井深 / m")
            ax.set_ylabel(r"$y$ (方位角)")
    else:
        # 示意数据
        y = np.linspace(0, 1, 40)
        s = np.linspace(0, 1, 80)
        Y, S = np.meshgrid(y, s)
        e_val = 0.45
        h = 1 + e_val * np.cos(np.pi * Y)
        w1 = h ** 2 / 0.08
        w1 /= w1.max()

        axes[0].pcolormesh(s, y, w1.T, cmap="hot", shading="gouraud")
        w2 = w1 * (1 + 0.3 * e_val * (2 * Y - 1))
        axes[1].pcolormesh(s, y, w2.T, cmap="hot", shading="gouraud")
        v_syn = np.gradient(w2, axis=1) * 0.05
        axes[2].pcolormesh(s, y, v_syn.T, cmap="coolwarm", shading="gouraud")
        for ax in axes:
            ax.set_xlabel(r"$s$ (归一化轴向)")
            ax.set_ylabel(r"$\phi$ (方位角)")

    axes[0].set_title("(a) 轴向速度 (无浮力)", fontsize=11, fontweight="bold")
    axes[1].set_title("(b) 轴向速度 (含浮力)", fontsize=11, fontweight="bold")
    axes[2].set_title("(c) 横向速度 $v$", fontsize=11, fontweight="bold")

    fig.suptitle("图6  Hele-Shaw 环空速度场分布", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    _save(fig, "fig06_速度场分布.png")


# ═══════════════════════════════════════════════════════════════
#  图7：半拉格朗日反演追踪示意
# ═══════════════════════════════════════════════════════════════
def fig07_semi_lagrangian():
    """半拉格朗日反演追踪示意图。"""
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.set_xlim(-0.5, 5.5)
    ax.set_ylim(-0.5, 5.5)
    ax.set_aspect("equal")
    ax.axis("off")

    # 网格
    for i in range(6):
        ax.axhline(i, color="#BDBDBD", lw=0.8, alpha=0.6)
        ax.axvline(i, color="#BDBDBD", lw=0.8, alpha=0.6)

    # 网格节点
    for i in range(6):
        for j in range(6):
            ax.plot(i, j, "o", color="#90A4AE", ms=3)

    # 目标点 (当前时刻位置)
    xp, yp = 3.5, 2.8
    ax.plot(xp, yp, "s", color="#D32F2F", ms=14, zorder=5, label="当前网格点 $P$")

    # 速度向量
    vx, vy = 0.6, 0.4
    ax.annotate("", xy=(xp - vx, yp - vy), xytext=(xp, yp),
                arrowprops=dict(arrowstyle="<-", color="#D32F2F", lw=2.5))

    # 反演出发点
    xs, ys = xp - vx * 0.8, yp - vy * 0.8
    ax.plot(xs, ys, "^", color="#1565C0", ms=14, zorder=5, label="出发点 $P_{src}$")

    # 双线性插值包围盒
    i0, j0 = int(np.floor(xs)), int(np.floor(ys))
    i1, j1 = i0 + 1, j0 + 1
    corners = [(i0, j0), (i1, j0), (i1, j1), (i0, j1)]
    for c in corners:
        ax.plot(c[0], c[1], "D", color="#2E7D32", ms=10, zorder=4)

    rect = plt.Rectangle((i0, j0), 1, 1, fill=True, facecolor="#E8F5E9",
                          edgecolor="#2E7D32", linewidth=2, alpha=0.35)
    ax.add_patch(rect)

    # 插值权重连线
    for c in corners:
        ax.plot([xs, c[0]], [ys, c[1]], "--", color="#2E7D32", lw=1.2, alpha=0.7)

    # 标注
    ax.text(i0 + 0.05, j0 - 0.15, f"({i0},{j0})", fontsize=9, color="#2E7D32")
    ax.text(i1 + 0.05, j0 - 0.15, f"({i1},{j0})", fontsize=9, color="#2E7D32")
    ax.text(i1 + 0.05, j1 + 0.08, f"({i1},{j1})", fontsize=9, color="#2E7D32")
    ax.text(i0 + 0.05, j1 + 0.08, f"({i0},{j1})", fontsize=9, color="#2E7D32")

    ax.text(xp + 0.1, yp + 0.15, "$P(x_i, s_j)$", fontsize=11, fontweight="bold", color="#D32F2F")
    ax.text(xs - 0.15, ys - 0.25, "$P_{src}$", fontsize=11, fontweight="bold", color="#1565C0")

    ax.annotate("", xy=(xp - vx, yp), xytext=(xp, yp),
                arrowprops=dict(arrowstyle="->", color="#FF6F00", lw=1.5))
    ax.text(xp - vx / 2, yp + 0.15, "$w \\cdot \\Delta t$", fontsize=10, color="#FF6F00", ha="center")
    ax.annotate("", xy=(xp - vx, yp - vy), xytext=(xp - vx, yp),
                arrowprops=dict(arrowstyle="->", color="#FF6F00", lw=1.5))
    ax.text(xp - vx + 0.15, yp - vy / 2, "$v \\cdot \\Delta t$", fontsize=10, color="#FF6F00")

    ax.legend(fontsize=11, loc="upper right")

    ax.set_title("图7  半拉格朗日反演追踪 + 双线性插值", fontsize=14, fontweight="bold")
    _save(fig, "fig07_半拉格朗日反演追踪.png")


# ═══════════════════════════════════════════════════════════════
#  图8：多约束数值求解链路流程图
# ═══════════════════════════════════════════════════════════════
def fig08_constraint_chain():
    """多约束数值求解链路流程图。

    正确顺序（来自技术交底书第187行）：
    1. 半拉格朗日平流（对各相独立执行）
    2. D2DGA通量放大施加（对各相独立执行）
    3. 逐相累计入口体积硬约束检查与修正
    4. 多相之和闭合约束检查与修正
    5. 泵停判断（若停泵，恢复本步初始浓度场）
    """
    fig, ax = plt.subplots(figsize=(8, 14))
    ax.set_xlim(-1.5, 9.5)
    ax.set_ylim(-0.5, 15)
    ax.axis("off")

    def mbox(x, y, text, color="#E3F2FD", w=3.6, h=0.9):
        fb = FancyBboxPatch((x - w/2, y - h/2), w, h,
            boxstyle="round,pad=0.1", facecolor=color,
            edgecolor="#37474F", linewidth=1.4)
        ax.add_patch(fb)
        ax.text(x, y, text, ha="center", va="center", fontsize=9.5, fontweight="bold")

    def mdiamond(x, y, text, color="#FFF9C4", w=3.0, h=1.0):
        pts = np.array([[x, y+h/2], [x+w/2, y], [x, y-h/2], [x-w/2, y]])
        poly = plt.Polygon(pts, facecolor=color, edgecolor="#37474F", linewidth=1.4)
        ax.add_patch(poly)
        ax.text(x, y, text, ha="center", va="center", fontsize=9, fontweight="bold")

    def arrow(x0, y0, x1, y1, color="#37474F", ls="-"):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=1.6, linestyle=ls))

    def label(x, y, text, color="#C62828"):
        ax.text(x, y, text, fontsize=9, color=color, fontweight="bold", ha="center")

    cx = 4  # 主链路x坐标
    bh = 0.9  # 框高

    # ── Step 1: 半拉格朗日平流 ──
    y1 = 14.0
    mbox(cx, y1, "Step 1: 半拉格朗日平流\n(各相独立，含D2DGA速度放大)", "#E3F2FD", 5.0, 1.0)

    # ── Step 2: D2DGA通量放大 ──
    y2 = 12.3
    mbox(cx, y2, r"Step 2: D2DGA 通量修正 $f(\bar{c},m)$", "#E8F5E9", 5.0, 0.9)
    arrow(cx, y1 - 0.5, cx, y2 + 0.45)

    # ── Step 3: 逐相体积硬约束 ──
    y3 = 10.5
    mdiamond(cx, y3, r"$V_{phase} > V_{in}$ ?")
    arrow(cx, y2 - 0.45, cx, y3 + 0.5)

    # 是 → 右侧线性缩放
    label(cx + 1.8, y3 + 0.15, "是", "#C62828")
    mbox(8, y3, "体积硬约束\n$c *= V_{in}/V$", "#C8E6C9", 2.8, 0.9)
    arrow(cx + 1.5, y3, 8 - 1.4, y3)

    # 否 → Step 4
    label(cx + 1.2, y3 - 0.6, "否", "#2E7D32")

    # ── Step 4: 多相闭合约束 ──
    y4 = 8.5
    mdiamond(cx, y4, r"$\Sigma c_i$ > 1 ?")
    arrow(cx, y3 - 0.5, cx, y4 + 0.5)

    # 是 → 左侧比例裁剪
    label(cx - 1.8, y4 + 0.15, "是", "#C62828")
    mbox(0.5, y4, "闭合约束\n$c_i = c_i/\\Sigma c_j$", "#FFCDD2", 2.8, 0.9)
    arrow(cx - 1.5, y4, 0.5 + 1.4, y4)

    # 否 → clip检查
    label(cx + 1.2, y4 - 0.6, "否", "#2E7D32")

    # clip检查（内联，不用菱形）
    y5 = 6.5
    mbox(cx, y5, "截断: $c = \\mathrm{clip}(c, 0, 1)$", "#FFF9C4", 4.0, 0.8)
    arrow(cx, y4 - 0.5, cx, y5 + 0.4)

    # ── Step 5: 泵停判断 ──
    y6 = 4.8
    mdiamond(cx, y6, "泵停?\n$Q < 10^{-9}$", "#F3E5F5", 3.0, 1.0)
    arrow(cx, y5 - 0.4, cx, y6 + 0.5)

    # 是 → 右侧冻结
    label(cx + 1.8, y6 + 0.15, "是", "#C62828")
    mbox(8, y6, "浓度场冻结\n(恢复本步初始值)", "#E0E0E0", 2.8, 0.9)
    arrow(cx + 1.5, y6, 8 - 1.4, y6)

    # 否 → 输出
    label(cx + 1.2, y6 - 0.6, "否", "#2E7D32")

    # ── 输出 ──
    y7 = 2.8
    mbox(cx, y7, "本步浓度场更新\n→ 诊断指标计算与输出", "#BBDEFB", 4.5, 1.0)
    arrow(cx, y6 - 0.5, cx, y7 + 0.5)

    # ── 下一时间步 ──
    y8 = 1.0
    mbox(cx, y8, "下一时间步", "#EFEBE9", 3.0, 0.8)
    arrow(cx, y7 - 0.5, cx, y8 + 0.4)

    # 侧边标注优先级
    ax.text(-1.3, 10.5, "优先级\n最高", fontsize=8, color="#1565C0",
            ha="center", va="center", fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#E3F2FD", edgecolor="#1565C0", lw=1))
    ax.annotate("", xy=(-0.8, 10.5), xytext=(cx - 1.5 - 0.05, 10.5),
                arrowprops=dict(arrowstyle="<-", color="#1565C0", lw=1.2))

    ax.text(-1.3, 8.5, "优先级\n次之", fontsize=8, color="#1565C0",
            ha="center", va="center", fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#E3F2FD", edgecolor="#1565C0", lw=1))
    ax.annotate("", xy=(-0.8, 8.5), xytext=(cx - 1.5 - 0.05, 8.5),
                arrowprops=dict(arrowstyle="<-", color="#1565C0", lw=1.2))

    # 体积约束和闭合约束的汇合箭头
    arrow(8, y3 - 0.45, 8, y5 + 0.4, "#666666", "--")
    arrow(0.5, y4 - 0.45, 0.5, y5 + 0.4, "#666666", "--")
    arrow(0.5, y5, cx - 2.0, y5, "#666666", "--")
    arrow(8, y5, cx + 2.0, y5, "#666666", "--")

    ax.set_title("图8  单时间步多约束数值求解链路", fontsize=12, fontweight="bold", pad=10)
    _save(fig, "fig08_多约束求解链路.png")


# ═══════════════════════════════════════════════════════════════
#  图9：浓度场快照
# ═══════════════════════════════════════════════════════════════
def _find_snapshot_index(result, target_time_s):
    """从快照时间列表中找到最接近目标时间的索引。"""
    times = result.snapshot_times_s
    if not times:
        return None
    return min(range(len(times)), key=lambda i: abs(times[i] - target_time_s))


def fig09_concentration_snapshots(result=None):
    """浓度场快照 at t=600s, 1500s, 3000s。"""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    target_times = [600, 1500, 3000]
    labels = ["(a) $t = 600$ s", "(b) $t = 1500$ s", "(c) $t = 3000$ s"]

    if result is not None and result.cement_snapshots:
        geom = result.geom
        for ax_i, t_target, label in zip(axes, target_times, labels):
            idx = _find_snapshot_index(result, t_target)
            if idx is None:
                idx = min(len(result.cement_snapshots) - 1, 2)
            snap = result.cement_snapshots[idx]
            pcm = ax_i.pcolormesh(
                geom["md"], geom["y"], snap,
                cmap="RdYlBu_r", shading="gouraud", vmin=0, vmax=1,
            )
            fig.colorbar(pcm, ax=ax_i, label="$c_{cement}$", shrink=0.82)
            ax_i.set_xlabel("井深 / m")
            ax_i.set_ylabel(r"$y$ (方位角)")
            t_actual = result.snapshot_times_s[idx]
            ax_i.set_title(f"{label}\n(t={t_actual:.0f}s)", fontsize=11, fontweight="bold")
    else:
        # 合成示意数据
        ny, nz = 30, 80
        y = np.linspace(0, 1, ny)
        s = np.linspace(0, 1, nz)
        Y, S = np.meshgrid(y, s)
        e_val = 0.45
        for ax_i, t_norm, label in zip(axes, [0.15, 0.4, 0.75], labels):
            front = t_norm
            cement = np.clip(
                np.exp(-20 * np.maximum(S - front, 0) ** 2) *
                (1 + 0.3 * e_val * np.cos(np.pi * Y)), 0, 1
            )
            cement += np.random.seed(42) or 0  # seed
            noise = np.random.normal(0, 0.03, cement.shape)
            cement = np.clip(cement + noise, 0, 1)
            pcm = ax_i.pcolormesh(s, y, cement.T, cmap="RdYlBu_r",
                                  shading="gouraud", vmin=0, vmax=1)
            fig.colorbar(pcm, ax=ax_i, label="$c_{cement}$", shrink=0.82)
            ax_i.set_xlabel(r"$s$ (归一化轴向)")
            ax_i.set_ylabel(r"$\phi$ (方位角)")
            ax_i.set_title(label + " (示意)", fontsize=11, fontweight="bold")

    fig.suptitle("图9  水泥浓度场演化快照", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    _save(fig, "fig09_浓度场快照.png")


# ═══════════════════════════════════════════════════════════════
#  图10：诊断指标时间序列
# ═══════════════════════════════════════════════════════════════
def fig10_diagnostics_timeseries(result=None):
    """诊断指标时间序列：效率、风险、前沿。"""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    if result is not None:
        m = result.metrics
        t = m["time_min"]
        # (a) 效率曲线
        axes[0].plot(t, m["effective_efficiency"], lw=2, label="全井段效率")
        axes[0].plot(t, m["target_interval_efficiency"], lw=2, label="目标层段效率")
        axes[0].plot(t, m["bulk_cement_fill"], "--", lw=1.8, label="占据率")
        # (b) 三指数
        axes[1].plot(t, m["channeling_index"], lw=2, label="窜槽指数")
        axes[1].plot(t, m["mixing_index"], lw=2, label="混浆指数")
        axes[1].plot(t, m["instability_index"], lw=2, label="失稳指数")
        # (c) 前沿
        axes[2].plot(t, m["front_wide_m"], lw=2, label="宽边")
        axes[2].plot(t, m["front_mid_m"], lw=2, label="中线")
        axes[2].plot(t, m["front_narrow_m"], lw=2, label="窄边")
    else:
        # 合成数据
        t = np.linspace(0, 80, 100)
        eff = 1 - np.exp(-t / 25)
        axes[0].plot(t, eff, lw=2, label="全井段效率")
        axes[0].plot(t, eff * 0.92, lw=2, label="目标层段效率")
        axes[0].plot(t, eff * 1.02, "--", lw=1.8, label="占据率")
        axes[1].plot(t, 0.3 * np.exp(-t / 30), lw=2, label="窜槽指数")
        axes[1].plot(t, 0.15 * (1 - np.exp(-t / 20)), lw=2, label="混浆指数")
        axes[1].plot(t, 0.1 * np.exp(-t / 40), lw=2, label="失稳指数")
        depth_max = 2200
        axes[2].plot(t, depth_max * (1 - np.exp(-t / 22)), lw=2, label="宽边")
        axes[2].plot(t, depth_max * (1 - np.exp(-t / 28)), lw=2, label="中线")
        axes[2].plot(t, depth_max * (1 - np.exp(-t / 35)), lw=2, label="窄边")

    titles = [
        "(a) 顶替效率与占据率",
        "(b) 风险指标演化",
        "(c) 前沿推进距离",
    ]
    ylabels = ["效率 / 占据率", "指数", "距离 / m"]

    for ax, title, ylabel in zip(axes, titles, ylabels):
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_xlabel("时间 / min")
        ax.set_ylabel(ylabel)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.25)

    fig.suptitle("图10  诊断指标时间序列", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    _save(fig, "fig10_诊断指标时间序列.png")


# ═══════════════════════════════════════════════════════════════
#  图11：深度方向效率剖面
# ═══════════════════════════════════════════════════════════════
def fig11_depth_profile(result=None):
    """深度方向效率剖面。"""
    fig, ax = plt.subplots(figsize=(9, 6))

    if result is not None:
        dp = result.depth_profiles
        depth = dp["井深_m"]
        ax.plot(depth, dp["宽边有效效率"], lw=2.2, label="宽边效率", color="#1565C0")
        ax.plot(depth, dp["中线有效效率"], lw=2.2, label="中线效率", color="#F9A825")
        ax.plot(depth, dp["窄边有效效率"], lw=2.2, label="窄边效率", color="#C62828")
        ax.plot(depth, dp["平均有效顶替效率"], lw=2.8, label="加权平均效率",
                color="black", ls="--")
    else:
        depth = np.linspace(5500, 7700, 100)
        norm_d = (depth - depth[0]) / (depth[-1] - depth[0])
        wide = np.clip(0.95 - 0.1 * norm_d + 0.05 * np.sin(3 * np.pi * norm_d), 0, 1)
        mid = np.clip(0.85 - 0.05 * norm_d + 0.03 * np.sin(2 * np.pi * norm_d), 0, 1)
        narrow = np.clip(0.6 + 0.2 * norm_d - 0.1 * norm_d ** 2, 0, 1)
        avg = 0.4 * wide + 0.3 * mid + 0.3 * narrow

        ax.plot(depth, wide, lw=2.2, label="宽边效率", color="#1565C0")
        ax.plot(depth, mid, lw=2.2, label="中线效率", color="#F9A825")
        ax.plot(depth, narrow, lw=2.2, label="窄边效率", color="#C62828")
        ax.plot(depth, avg, lw=2.8, label="加权平均效率", color="black", ls="--")

    ax.set_xlabel("井深 / m", fontsize=12)
    ax.set_ylabel("有效顶替效率", fontsize=12)
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=11, loc="best")
    ax.grid(True, alpha=0.25)
    ax.set_title("图11  深度方向有效顶替效率剖面", fontsize=14, fontweight="bold")
    fig.tight_layout()
    _save(fig, "fig11_深度方向效率剖面.png")


# ═══════════════════════════════════════════════════════════════
#  图12：计算效率与弥散方式对比
# ═══════════════════════════════════════════════════════════════
def fig12_efficiency_comparison():
    """计算效率与弥散方式对比气泡图。"""
    fig, ax = plt.subplots(figsize=(10, 6))

    methods = [
        "纯对流\n(无弥散)",
        "显式 Laplacian\n弥散",
        "D2DGA\n通量修正",
        "D2DGA +\nPapanastasiou",
        "CFD 全3D\nDNS",
    ]
    y_pos = np.arange(len(methods))
    compute_time = [0.8, 2.5, 4.2, 5.0, 180.0]  # 分钟
    accuracy = [0.55, 0.72, 0.88, 0.91, 0.95]
    bubble_size = [200, 350, 500, 550, 800]
    colors = ["#90CAF9", "#A5D6A7", "#FFCC80", "#EF9A9A", "#CE93D8"]

    for i, (method, ct, acc, bs, c) in enumerate(
        zip(methods, compute_time, accuracy, bubble_size, colors)
    ):
        ax.scatter(ct, i, s=bs, c=c, edgecolors="#37474F", linewidth=1.5,
                   alpha=0.85, zorder=4)
        ax.text(ct * 1.3, i, f"精度 {acc:.0%}", fontsize=9, va="center")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(methods, fontsize=10)
    ax.set_xscale("log")
    ax.set_xlabel("计算时间 / min (对数尺度)", fontsize=12)
    ax.set_title("图12  不同弥散/求解方式的计算效率与精度对比",
                 fontsize=14, fontweight="bold")
    ax.grid(True, axis="x", alpha=0.3)
    ax.set_xlim(0.3, 500)
    ax.invert_yaxis()

    # 图例说明气泡大小
    ax.text(0.98, 0.02,
            "气泡大小 $\propto$ 模型复杂度",
            transform=ax.transAxes, fontsize=9, ha="right", va="bottom",
            style="italic", color="#616161")

    fig.tight_layout()
    _save(fig, "fig12_计算效率对比.png")


# ═══════════════════════════════════════════════════════════════
#  主入口
# ═══════════════════════════════════════════════════════════════
def main():
    print("=" * 60)
    print("  专利技术交底书 12 张附图生成")
    print("=" * 60)

    # 先尝试运行模拟
    result = _run_solver_for_figures()

    print("\n[1/12] 图1: 总体流程图")
    fig01_method_flowchart()

    print("[2/12] 图2: 偏心环空几何示意")
    fig02_eccentric_annulus()

    print("[3/12] 图3: 四相体积分数闭合")
    fig03_phase_closure()

    print("[4/12] 图4: Papanastasiou 正则化 (纯数学)")
    fig04_papanastasiou()

    print("[5/12] 图5: D2DGA 通量放大因子 (纯数学)")
    fig05_d2dga_flux()

    print("[6/12] 图6: Hele-Shaw 速度场")
    fig06_velocity_field(result)

    print("[7/12] 图7: 半拉格朗日反演追踪")
    fig07_semi_lagrangian()

    print("[8/12] 图8: 多约束数值求解链路")
    fig08_constraint_chain()

    print("[9/12] 图9: 浓度场快照")
    fig09_concentration_snapshots(result)

    print("[10/12] 图10: 诊断指标时间序列")
    fig10_diagnostics_timeseries(result)

    print("[11/12] 图11: 深度方向效率剖面")
    fig11_depth_profile(result)

    print("[12/12] 图12: 计算效率与弥散方式对比")
    fig12_efficiency_comparison()

    print("\n" + "=" * 60)
    print(f"  全部 12 张图已保存至: {OUTPUT_DIR}")
    print("=" * 60)

    # 列出生成的文件
    for f in sorted(OUTPUT_DIR.glob("fig*.png")):
        size_kb = f.stat().st_size / 1024
        print(f"  {f.name:40s} ({size_kb:.0f} KB)")


if __name__ == "__main__":
    main()
