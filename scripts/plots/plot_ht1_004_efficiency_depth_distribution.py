"""
呼1-004尾管段顶替效率沿深度分布图（学术风格 · 简洁版）
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

# ---- 路径 ----
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV = os.path.join(ROOT, "results", "呼1-004_1D2D耦合模型", "呼1-004_1D2D耦合模型_深度剖面.csv")
OUT = os.path.join(ROOT, "results", "呼1-004_1D2D耦合模型", "呼1-004_顶替效率沿深度分布图.png")

# 按实际井深/工程界限分段，不严格等距
SEGS = [
    (5243.2, 5578.0, "5243 – 5578 m\n套管重叠段"),
    (5578.0, 6000.0, "5578 – 6000 m\n上裸眼段"),
    (6000.0, 6600.0, "6000 – 6600 m\n裸眼段"),
    (6600.0, 7000.0, "6600 – 7000 m\n领浆封固段"),
    (7000.0, 7378.1, "7000 – 7378 m\n领浆封固段"),
    (7378.1, 7521.0, "7378 – 7521 m\n139.7mm 尾管段"),
    (7521.0, 7660.0, "7521 – 7660 m\n139.7mm 尾管段"),
]

# ---- 数据 ----
def load():
    rows = []
    with open(CSV, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            rows.append({
                "d": float(r["井深_m"]), "e": float(r["平均有效顶替效率"]),
                "ca": float(r["水泥平均浓度"]), "cw": float(r["宽边水泥浓度"]),
                "cn": float(r["窄边水泥浓度"]),
            })
    segs = []
    for top, bot, label in SEGS:
        pts = [r for r in rows if top <= r["d"] < bot]
        if not pts: continue
        n = len(pts)
        e = sum(r["e"] for r in pts) / n
        g = "良好" if e >= 0.80 else ("合格" if e >= 0.70 else "不合格")
        segs.append({"top": top, "bot": bot, "mid": (top+bot)/2,
                      "label": label, "eff": e*100,
                      "ca": sum(r["ca"] for r in pts)/n*100,
                      "cw": sum(r["cw"] for r in pts)/n*100,
                      "cn": sum(r["cn"] for r in pts)/n*100,
                      "grade": g})
    return segs, rows

# ---- 绘图 ----
def plot(segs, profile):
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
        # 效率数值 — 放在柱内右侧
        ax.text(val - 3, mid, f"{val:.1f}%", ha="right", va="center",
                fontsize=14, fontweight="bold", color="white",
                bbox=dict(facecolor=c, edgecolor="none", alpha=0.0, pad=2))
        # 段标签 — 放在柱左侧紧贴
        ax.text(2, mid, s["label"], ha="left", va="center",
                fontsize=12, color="#37474F", fontweight="bold",
                linespacing=1.3)

    # 阈值线
    for th, clr in [(80, GRADE_C["良好"]), (70, GRADE_C["合格"])]:
        ax.axvline(x=th, color=clr, linestyle=(0, (8, 4)), linewidth=1.2, alpha=0.5, zorder=1)

    # 关键深度线
    for z in [5578, 7378.1, 7521]:
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
    ax2.plot(np.array([r["cw"] for r in profile])*100, Z,
             color=LINE_C["wide"], linewidth=2.0, label="宽边", alpha=0.85)
    ax2.plot(np.array([r["ca"] for r in profile])*100, Z,
             color=LINE_C["avg"], linewidth=2.6, label="周向平均", alpha=0.95)
    ax2.plot(np.array([r["cn"] for r in profile])*100, Z,
             color=LINE_C["narrow"], linewidth=2.0, label="窄边", alpha=0.85)

    # 分级背景
    ax2.axvspan(0, 70, alpha=0.05, color=GRADE_C["不合格"], zorder=0)
    ax2.axvspan(70, 80, alpha=0.05, color=GRADE_C["合格"], zorder=0)
    ax2.axvspan(80, 108, alpha=0.05, color=GRADE_C["良好"], zorder=0)

    for z, lbl in [(5578, "套管鞋"), (7378.1, "变径 7378 m"), (7521, "变径 7521 m")]:
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

    # 标题
    overall = sum(s["eff"] * (s["bot"] - s["top"]) for s in segs) / (z1 - z0)
    fig.suptitle("呼 1-004 井  尾管段固井顶替效率沿深度分布",
                 fontsize=22, fontweight="bold", y=0.98)
    gN = sum(1 for s in segs if s["grade"] == "良好")
    oN = sum(1 for s in segs if s["grade"] == "合格")
    bN = sum(1 for s in segs if s["grade"] == "不合格")
    fig.text(0.5, 0.03,
             f"加权平均效率  {overall:.1f}%    |    "
             f"良好 {gN} 段  ·  合格 {oN} 段  ·  不合格 {bN} 段    |    "
             f"1D-2D 耦合 D2DGA 环空顶替模型    |    域 5243–7660 m  168.3 + 139.7 mm 双径尾管",
             ha="center", fontsize=11, color="#546E7A")

    plt.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"OK  {OUT}")

if __name__ == "__main__":
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    segs, profile = load()
    for s in segs:
        print(f"  {s['top']:.0f}–{s['bot']:.0f}m  {s['eff']:.1f}%  {s['grade']}")
    plot(segs, profile)
