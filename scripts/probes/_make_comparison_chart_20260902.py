"""2026-09-02 生成八井修正前后 η_E 对比图与汇总 CSV（中文标注，遵循 AGENTS.md 图表规范）。"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

# Windows 中文字体回退：优先微软雅黑/黑体，缺失时回退默认并保证负号正常
for _f in ["Microsoft YaHei", "SimHei", "SimSun", "Arial Unicode MS"]:
    if any(_f in f.name for f in font_manager.fontManager.ttflist):
        plt.rcParams["font.sans-serif"] = [_f]
        break
plt.rcParams["axes.unicode_minus"] = False

OUT = Path(__file__).resolve().parents[2] / "results" / "_质量平衡取证_2026-09-02"
OUT.mkdir(parents=True, exist_ok=True)

# nz=250，三处修复（速度归一化因子2 + 物理τw屈服门 + 双内径停止时刻）后的最终结果
df = pd.DataFrame([
    # 井, 修正前ηE, 修正后ηE, 窄边ηN, 平均居中度
    ("呼101", 0.3992, 0.5348, 0.0305, 0.438),
    ("呼102", 0.3814, 0.9687, 0.8921, 0.491),
    ("呼103", 0.1940, 0.8529, 0.7357, 0.719),
    ("呼1",   0.5014, 0.9494, 0.8295, 0.650),
    ("呼2",   0.5690, 0.9843, 0.9897, 0.778),
    ("HT1-001", 0.5195, 0.9888, 0.9764, 0.804),
    ("HT1-003", 0.6132, 0.9976, 0.9900, 0.830),
    ("HT1-004", 0.6284, 0.9985, 0.9938, 0.830),
], columns=["井", "修正前顶替效率", "修正后顶替效率", "修正后窄边效率ηN", "平均居中度"])
df["提升幅度"] = df["修正后顶替效率"] - df["修正前顶替效率"]
df.to_csv(OUT / "八井顶替效率_修正前后对照_nz250.csv", index=False, encoding="utf-8-sig")

fig, ax = plt.subplots(figsize=(11, 6), dpi=150)
x = np.arange(len(df)); w = 0.38
b1 = ax.bar(x - w/2, df["修正前顶替效率"], w, label="修正前 η_E", color="#C0504D", edgecolor="black", linewidth=0.5)
b2 = ax.bar(x + w/2, df["修正后顶替效率"], w, label="修正后 η_E", color="#4F81BD", edgecolor="black", linewidth=0.5)
ax.plot(x, df["平均居中度"], "D-", color="#2E7D32", lw=1.8, ms=6, label="平均居中度（右轴口径参考）")
for b in list(b1) + list(b2):
    ax.text(b.get_x()+b.get_width()/2, b.get_height()+0.012, f"{b.get_height():.2f}",
            ha="center", va="bottom", fontsize=9)
ax.axhline(0.85, ls="--", c="gray", lw=1.2)
ax.text(len(df)-0.5, 0.862, "文献常见合格区间下限≈0.85", ha="right", fontsize=9, color="gray")
ax.set_xticks(x); ax.set_xticklabels(df["井"])
ax.set_ylim(0, 1.12); ax.set_ylabel("顶替效率 η_E（水泥体积占据率）")
ax.set_title("八口尾管井顶替效率：三处守恒/时序缺陷修正前后对照（nz=250）", fontsize=13, fontweight="bold")
ax.legend(loc="upper left", ncol=3, fontsize=9)
ax.grid(axis="y", ls=":", alpha=0.5)
fig.tight_layout()
fig.savefig(OUT / "八井顶替效率_修正前后对照.png", bbox_inches="tight")
print("saved:", OUT / "八井顶替效率_修正前后对照.png")
print(df.to_string(index=False))
