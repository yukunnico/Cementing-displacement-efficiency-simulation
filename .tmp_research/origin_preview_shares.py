# -*- coding: utf-8 -*-
"""占比-时间图 matplotlib 预览（dataviz 定稿用，先于 Origin 版本验证设计）。
数据：results/<井>_1D2D耦合模型/<井>_1D2D耦合模型_全井深度时间占比_长格式.csv
色板：NPG 校准版（validate_palette.js 通过）mud=#375E99 spacer=#2FA8CC lead=#009463 tail=#E64B35
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

COLORS = {"mud": "#375E99", "spacer": "#2FA8CC", "lead": "#009463", "tail": "#E64B35"}
NAMES = {"mud": "钻井液", "spacer": "隔离液", "lead": "领浆", "tail": "尾浆"}
# 每井 4 个代表深度（悬挂器下 / 上部 / 下部 / 尾管鞋附近）
DEPTHS = {
    "呼1-003": [5310, 6000, 7000, 7580],
    "呼1-004": [5250, 6000, 7000, 7650],
}
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei"]
plt.rcParams["axes.unicode_minus"] = False

fig, axes = plt.subplots(2, 4, figsize=(16, 7.2), dpi=160, sharey=True)
for r, well in enumerate(["呼1-003", "呼1-004"]):
    lg = pd.read_csv(rf"results/{well}_1D2D耦合模型/{well}_1D2D耦合模型_全井深度时间占比_长格式.csv")
    lg = lg[lg.source == "2D_annulus"]
    lg["t_min"] = lg["time_s"] / 60
    all_d = lg.depth_m.unique()
    for c, d in enumerate(DEPTHS[well]):
        ax = axes[r, c]
        dn = min(all_d, key=lambda x: abs(x - d))   # 就近匹配表格深度
        sub = lg[lg.depth_m == dn].pivot_table(index="t_min", columns="fluid", values="share")
        for fl in ["mud", "spacer", "lead", "tail"]:
            if fl in sub:
                ax.plot(sub.index, sub[fl], color=COLORS[fl], lw=2, solid_capstyle="round")
        ax.set_title(f"深度 {d} m", fontsize=11)
        ax.set_xlim(0, sub.index.max())
        ax.set_ylim(-0.02, 1.05)
        ax.grid(True, color="#d9d9d9", lw=0.6, alpha=0.8)
        ax.set_axisbelow(True)
        for s in ["top", "right"]:
            ax.spines[s].set_visible(False)
        if c == 0:
            ax.set_ylabel(f"{well}\n体积占比", fontsize=11)
        if r == 1:
            ax.set_xlabel("时间 (min)", fontsize=10)
        # 最右列线端直接标注（墨色文字，非系列色）
        if c == 3:
            for fl in ["mud", "spacer", "lead", "tail"]:
                if fl in sub and sub[fl].iloc[-1] > 0.04:
                    ax.annotate(NAMES[fl], (sub.index[-1], sub[fl].iloc[-1]),
                                xytext=(4, 0), textcoords="offset points",
                                fontsize=9, color="#333333", va="center")
handles = [plt.Line2D([], [], color=COLORS[f], lw=2) for f in COLORS]
fig.legend(handles, [NAMES[f] for f in COLORS], loc="upper center",
           ncol=4, frameon=False, fontsize=11, bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(rect=(0, 0, 1, 0.94))
out = ".tmp_research/占比时间预览_两井_matplotlib.png"
fig.savefig(out, bbox_inches="tight")
print("saved:", out)
