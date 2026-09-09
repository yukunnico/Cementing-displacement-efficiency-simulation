# -*- coding: utf-8 -*-
"""Origin 版占比-时间图 终版（v5 布局 + add_label 文本对象，绕开 LabTalk 转义）。

产出：.tmp_research/占比时间_Origin工程.opju + .tmp_research/origin_out/占比时间_Origin.png
"""
from pathlib import Path
import pandas as pd
import originpro as op
from originpro.graph import ocolor

ROOT = Path(__file__).resolve().parents[1]
COLORS = {"mud": "#375E99", "spacer": "#2FA8CC", "lead": "#009463", "tail": "#E64B35"}
NAMES = {"mud": "钻井液", "spacer": "隔离液", "lead": "领浆", "tail": "尾浆"}
DEPTHS = {"呼1-003": [5310, 6000, 7000, 7580], "呼1-004": [5250, 6000, 7000, 7650]}
WELLS = ["呼1-003", "呼1-004"]

op.new()
wb = op.new_book("w", lname="占比时间数据")

sheet_meta = []
for well in WELLS:
    lg = pd.read_csv(ROOT / rf"results/{well}_1D2D耦合模型/{well}_1D2D耦合模型_全井深度时间占比_长格式.csv")
    lg = lg[lg.source == "2D_annulus"]
    lg["t_min"] = lg.time_s / 60
    all_d = sorted(lg.depth_m.unique())
    for d in DEPTHS[well]:
        dn = min(all_d, key=lambda x: abs(x - d))
        sub = lg[lg.depth_m == dn].pivot_table(index="t_min", columns="fluid", values="share")
        sh = wb.add_sheet()
        sh.name = f"{well[-3:]}_{int(dn)}m"
        sh.from_df(sub.reset_index()[["t_min", "mud", "spacer", "lead", "tail"]])
        sheet_meta.append((well, int(dn), sh))

gp = op.new_graph(lname="占比时间")
for si in range(7):
    gp.add_layer()

def lt(i, cmd):
    gp.obj.LT_execute(f"page.active={i+1}; " + cmd)

# ── 1. 2×4 布局（页面百分比）──
MARGIN_L, TOP = 11, 10
col_w, col_step = 19.5, 22.0
row_h, row_step = 37.0, 44.0
for i in range(8):
    r, c = divmod(i, 4)
    lt(i, f"layer.left={MARGIN_L + c*col_step}; layer.top={TOP + r*row_step}; "
          f"layer.width={col_w}; layer.height={row_h};")

# ── 2. 每层：曲线 + 轴 + 样式 ──
for si, (well, dn, sh) in enumerate(sheet_meta):
    r, c = divmod(si, 4)
    lay = gp[si]
    for fl, col in [("mud", 1), ("spacer", 2), ("lead", 3), ("tail", 4)]:
        p = lay.add_plot(sh, colx=0, coly=col, type="line")
        p.color = ocolor(COLORS[fl])
        p.set_int("Width", 2.5)
    lay.set_xlim(0, 200)
    lay.set_ylim(0, 1.05)
    lt(si, "layer.x.inc=50; layer.y.inc=0.25;")
    lt(si, 'xb.text$="时间 (min)"' if r == 1 else 'xb.text$=""')
    lt(si, 'yl.text$="体积占比"' if c == 0 else 'yl.text$=""')
    lt(si, "xb.fSize=16; yl.fSize=16;")
    lt(si, "legend.x=-100; legend.y=-100;")   # 藏默认图例
    # 面板深度标题：层内 add_label（数据坐标，层顶中部）
    lab = lay.add_label(f"深度 {dn} m", 100, 0.99)
    if lab:
        lab.set_int("fSize", 18)
        lab.set_int("attach", 2)   # 2 = page 单位？先试

# ── 3. 页顶总图例：第一层 add_label（层内数据坐标 y>1.0 处），四相各一条 ──
lay0 = gp[0]
legend_items = [("mud", "钻井液"), ("spacer", "隔离液"), ("lead", "领浆"), ("tail", "尾浆")]
for k, (fl, nm) in enumerate(legend_items):
    lab = lay0.add_label(nm, 15 + k * 55, 1.13)
    if lab:
        lab.set_int("fSize", 20)

op.save(str(ROOT / ".tmp_research/占比时间_Origin工程.opju"))
out = ROOT / ".tmp_research/origin_out"
out.mkdir(exist_ok=True)
gp.save_fig(str(out / "占比时间_Origin.png"), type="png", width=1920)
print("OK")
op.exit()
