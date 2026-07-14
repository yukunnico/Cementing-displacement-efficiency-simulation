"""
各井尾管段顶替效率沿深度分布图（学术风格 · 简洁版）

基于 plot_ht1_004_efficiency_depth_distribution.py 改写为多井参数化版本，
分别绘制：呼探1井、呼1-001、呼1-002、呼1-003 四口井的尾管段顶替效率
沿深度分布图。各井按深度每 300 m 等距分段，计算各段平均顶替效率；
段标签只写深度范围，不起段名。
"""

import csv, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.font_manager as fm
import numpy as np
from matplotlib.patches import Patch

def _setup_font():
    available = {f.name.lower(): f.name for f in fm.fontManager.ttflist}
    # 强制 Microsoft YaHei，找不到则用 SimHei
    if "microsoft yahei" in available:
        chosen = available["microsoft yahei"]
    elif "simhei" in available:
        chosen = available["simhei"]
    else:
        chosen = "sans-serif"
    plt.rcParams["font.sans-serif"] = [chosen, "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    plt.rcParams["axes.titleweight"] = "bold"
    print(f"[字体] {chosen}")
_setup_font()

# ---- 配色 ----
GRADE_C = {"良好": "#2E7D32", "合格": "#EF6C00", "不合格": "#C62828"}
LINE_C = {"avg": "#1565C0", "wide": "#2E7D32", "narrow": "#C62828"}

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# ============================================================
#  各井配置：关键深度线（套管鞋/变径等井身特征）、标题、页脚说明
#  分段不在此配置，由 load() 按深度每 STEP m 等距自动切分
#  key_lines:  左图仅画虚线的深度
#  key_labels: 右图标注 (z, 文本)
# ============================================================
STEP = 300.0  # 每段深度步长 (m)

WELLS = [
    # ---- 呼探1井（HU1）168.3+139.7 双径尾管 ----
    {
        "name": "呼探 1 井",
        "csv": os.path.join(ROOT, "results", "呼探1尾管_1D2D耦合模型",
                            "呼探1尾管_1D2D耦合模型_深度剖面.csv"),
        "out": os.path.join(ROOT, "results", "呼探1尾管_1D2D耦合模型",
                            "呼探1_顶替效率沿深度分布图.png"),
        "key_lines": [7174.938],
        "key_labels": [(7174.938, "变径 7175 m")],
        "footer": "域 5470–7746 m  168.3 + 139.7 mm 双径尾管",
    },
    # ---- 呼1-001（HT1-001）168.3+139.7 双径尾管 ----
    {
        "name": "呼 1-001 井",
        "csv": os.path.join(ROOT, "results", "呼探1-001尾管_1D2D耦合模型",
                            "呼探1-001尾管_1D2D耦合模型_深度剖面.csv"),
        "out": os.path.join(ROOT, "results", "呼探1-001尾管_1D2D耦合模型",
                            "呼探1-001_顶替效率沿深度分布图.png"),
        "key_lines": [7174.938],
        "key_labels": [(7174.938, "变径 7175 m")],
        "footer": "域 5470–7746 m  168.3 + 139.7 mm 双径尾管",
    },
    # ---- 呼1-002（HT1-002）单一 139.7mm 尾管 ----
    {
        "name": "呼 1-002 井",
        "csv": os.path.join(ROOT, "results", "呼探1-002尾管_1D2D耦合模型",
                            "呼探1-002尾管_1D2D耦合模型_深度剖面.csv"),
        "out": os.path.join(ROOT, "results", "呼探1-002尾管_1D2D耦合模型",
                            "呼探1-002_顶替效率沿深度分布图.png"),
        "key_lines": [7355.0],
        "key_labels": [(7355.0, "大肚子 7355 m")],
        "footer": "域 5293–7554 m  139.7 mm 单径尾管",
    },
    # ---- 呼1-003（HT1-003）168.3+139.7 双径尾管 ----
    {
        "name": "呼 1-003 井",
        "csv": os.path.join(ROOT, "results", "呼1-003_1D2D耦合模型",
                            "呼1-003_1D2D耦合模型_深度剖面.csv"),
        "out": os.path.join(ROOT, "results", "呼1-003_1D2D耦合模型",
                            "呼1-003_顶替效率沿深度分布图.png"),
        "key_lines": [5568.0, 7089.576, 7096.0],
        "key_labels": [(5568.0, "套管鞋"),
                       (7089.576, "变径 7090 m"),
                       (7096.0, "井眼变径 7096 m")],
        "footer": "域 5308–7618 m  168.3 + 139.7 mm 双径尾管",
    },
]


# ---- 数据 ----
def load(cfg):
    rows = []
    with open(cfg["csv"], encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            rows.append({
                "d": float(r["井深_m"]), "e": float(r["平均有效顶替效率"]),
                "ca": float(r["水泥平均浓度"]), "cw": float(r["宽边水泥浓度"]),
                "cn": float(r["窄边水泥浓度"]),
            })
    # 从数据域顶起每 STEP m 一段，末段收口到域底
    z_top = min(r["d"] for r in rows)
    z_bot = max(r["d"] for r in rows)
    bounds = []
    z = z_top
    while z < z_bot - 1e-6:
        bounds.append((z, min(z + STEP, z_bot)))
        z += STEP
    n = len(bounds)
    segs = []
    for i, (top, bot) in enumerate(bounds):
        # 末段用闭区间，确保域底端点不被遗漏
        if i == n - 1:
            pts = [r for r in rows if top <= r["d"] <= bot]
        else:
            pts = [r for r in rows if top <= r["d"] < bot]
        if not pts:
            continue
        m = len(pts)
        e = sum(r["e"] for r in pts) / m
        g = "良好" if e >= 0.80 else ("合格" if e >= 0.70 else "不合格")
        segs.append({"top": top, "bot": bot, "mid": (top + bot) / 2,
                     "label": f"{top:.0f} – {bot:.0f} m",
                     "eff": e * 100,
                     "ca": sum(r["ca"] for r in pts) / m * 100,
                     "cw": sum(r["cw"] for r in pts) / m * 100,
                     "cn": sum(r["cn"] for r in pts) / m * 100,
                     "grade": g})
    return segs, rows


# ---- 绘图 ----
def plot(cfg, segs, profile):
    z0 = min(s["top"] for s in segs)
    z1 = max(s["bot"] for s in segs)

    fig = plt.figure(figsize=(22, 13), facecolor="white")
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.2, 1], wspace=0.18,
                           left=0.06, right=0.97, top=0.92, bottom=0.10)

    # ======= 左：条形图 =======
    ax = fig.add_subplot(gs[0, 0])
    ax.set_facecolor("white")

    for s in segs:
        c = GRADE_C[s["grade"]]
        top, bot, mid, val, g = s["top"], s["bot"], s["mid"], s["eff"], s["grade"]
        h = bot - top
        # 柱体
        ax.barh(mid, val, height=h, color=c, edgecolor="#333333",
                linewidth=0.6, alpha=0.85, zorder=3)
        # 效率数值 - 放在柱内右侧
        ax.text(val - 3, mid, f"{val:.1f}%", ha="right", va="center",
                fontsize=14, fontweight="bold", color="white",
                bbox=dict(facecolor=c, edgecolor="none", alpha=0.0, pad=2))
        # 段标签 - 放在柱左侧紧贴（只写深度范围）
        ax.text(2, mid, s["label"], ha="left", va="center",
                fontsize=12, color="#37474F", fontweight="bold")

    # 阈值线
    for th, clr in [(80, GRADE_C["良好"]), (70, GRADE_C["合格"])]:
        ax.axvline(x=th, color=clr, linestyle=(0, (8, 4)), linewidth=1.2, alpha=0.5, zorder=1)

    # 关键深度线（各井不同）
    for z in cfg["key_lines"]:
        ax.axhline(y=z, color="#B0BEC5", linestyle=(0, (3, 3)), linewidth=0.6, alpha=0.4)

    ax.set_ylim(z1 + 40, z0 - 40)
    ax.set_xlim(0, 108)
    ax.set_xlabel("模型计算顶替效率  /  %", fontsize=15, fontweight="bold", labelpad=10)
    ax.set_ylabel("测深  /  m", fontsize=15, fontweight="bold", labelpad=10)
    ax.tick_params(labelsize=12)
    ax.grid(axis="x", color="#CFD8DC", linestyle="--", linewidth=0.4, alpha=0.6)

    # 图例
    ax.legend(handles=[
        Patch(facecolor=GRADE_C["良好"], edgecolor="#333", label="良好  ≥ 80%"),
        Patch(facecolor=GRADE_C["合格"], edgecolor="#333", label="合格  70–80%"),
        Patch(facecolor=GRADE_C["不合格"], edgecolor="#333", label="不合格  < 70%"),
    ], loc="lower right", fontsize=12, framealpha=0.9, edgecolor="#B0BEC5",
       title="顶替效率分级", title_fontsize=13)

    # ======= 右：水泥体积分数剖面 =======
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_facecolor("white")

    Z = np.array([r["d"] for r in profile])
    ax2.plot(np.array([r["cw"] for r in profile]) * 100, Z,
             color=LINE_C["wide"], linewidth=2.0, label="宽边", alpha=0.85)
    ax2.plot(np.array([r["ca"] for r in profile]) * 100, Z,
             color=LINE_C["avg"], linewidth=2.6, label="周向平均", alpha=0.95)
    ax2.plot(np.array([r["cn"] for r in profile]) * 100, Z,
             color=LINE_C["narrow"], linewidth=2.0, label="窄边", alpha=0.85)

    # 分级背景
    ax2.axvspan(0, 70, alpha=0.05, color=GRADE_C["不合格"], zorder=0)
    ax2.axvspan(70, 80, alpha=0.05, color=GRADE_C["合格"], zorder=0)
    ax2.axvspan(80, 108, alpha=0.05, color=GRADE_C["良好"], zorder=0)

    # 关键深度线 + 标注（各井不同）
    for z, lbl in cfg["key_labels"]:
        ax2.axhline(y=z, color="#B0BEC5", linestyle=(0, (3, 3)), linewidth=0.6, alpha=0.4)
        ax2.text(104, z, lbl, fontsize=9, color="#78909C", va="center",
                 bbox=dict(facecolor="white", edgecolor="none", alpha=0.8))

    ax2.set_ylim(z1 + 40, z0 - 40)
    ax2.set_xlim(0, 108)
    ax2.set_xlabel("水泥体积分数  /  %", fontsize=15, fontweight="bold", labelpad=10)
    ax2.set_ylabel("测深  /  m", fontsize=15, fontweight="bold", labelpad=10)
    ax2.tick_params(labelsize=12)
    ax2.grid(axis="x", color="#CFD8DC", linestyle="--", linewidth=0.4, alpha=0.6)
    ax2.legend(loc="lower left", fontsize=12, framealpha=0.9, edgecolor="#B0BEC5")

    # 标题 + 页脚
    overall = sum(s["eff"] * (s["bot"] - s["top"]) for s in segs) / (z1 - z0)
    fig.suptitle(f"{cfg['name']}  尾管段固井顶替效率沿深度分布",
                 fontsize=22, fontweight="bold", y=0.98)
    gN = sum(1 for s in segs if s["grade"] == "良好")
    oN = sum(1 for s in segs if s["grade"] == "合格")
    bN = sum(1 for s in segs if s["grade"] == "不合格")
    fig.text(0.5, 0.03,
             f"加权平均效率  {overall:.1f}%    |    "
             f"良好 {gN} 段  ·  合格 {oN} 段  ·  不合格 {bN} 段    |    "
             f"1D-2D 耦合 D2DGA 环空顶替模型    |    {cfg['footer']}",
             ha="center", fontsize=11, color="#546E7A")

    os.makedirs(os.path.dirname(cfg["out"]), exist_ok=True)
    plt.savefig(cfg["out"], dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"OK  {cfg['out']}")


if __name__ == "__main__":
    for cfg in WELLS:
        segs, profile = load(cfg)
        print(f"\n=== {cfg['name']} ===")
        for s in segs:
            print(f"  {s['top']:.1f}–{s['bot']:.1f}m  {s['eff']:.1f}%  {s['grade']}")
        plot(cfg, segs, profile)
