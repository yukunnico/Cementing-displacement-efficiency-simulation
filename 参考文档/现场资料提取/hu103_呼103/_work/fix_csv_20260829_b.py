# -*- coding: utf-8 -*-
"""2026-08-29 hu103 提取包 CSV 校准修正（第二轮：修正列索引，防重复）。"""
import csv
from pathlib import Path

base = Path(r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\hu103_呼103")
log = []

def load(name):
    with open(base / name, encoding="utf-8-sig", newline="") as f:
        return list(csv.reader(f))

def save(name, rows):
    with open(base / name, "w", encoding="utf-8-sig", newline="") as f:
        csv.writer(f, lineterminator="\n").writerows(rows)

# ---- A. cbl_evaluation.csv 综合行（pass_rate=r[4]，区间顶 r[2]）----
rows = load("cbl_evaluation.csv")
hit = False
for r in rows:
    if len(r) > 13 and r[0] == "hu103" and r[4] == "6.05":
        r[4] = "2.12"
        r[2] = "7712.0"
        r[6] = ("168.3mm段(0.04%)+139.7mm段(12.06%)长度加权：按100515/100516测量井段5540-7330.6m(1790.6m)与"
                "7338-7712m(374m)，综合=(0.04*1790.6+12.06*374)/2164.6≈2.12%。"
                "2026-08-29校准：原值6.05%系两段简单平均(0.04+12.06)/2误标为加权，其原注公式实算2.16%，"
                "现按0708 PDF测量井段改为长度加权2.12%。")
        r[9] = "综合计算（168.3mm+139.7mm长度加权，2026-08-29校准）"
        r[13] = "整段尾管综合合格率约2.12%——模型顶替效率应与此对比验证。"
        hit = True
        log.append(("cbl_evaluation.csv", "综合行 pass_rate/区间顶", "6.05(简单平均,原公式实算2.16) / 7770.0", "2.12(长度加权,测量井段1790.6m+374m) / 7712.0"))
        break
if not hit:
    print("WARN: cbl_evaluation 综合行未命中（检查是否已改）")
save("cbl_evaluation.csv", rows)

# ---- C/D. pumping_schedule.csv（volume=r[5]，stage_name=r[2]）----
rows = load("pumping_schedule.csv")
h_c = h_d = False
for r in rows:
    if len(r) > 15 and r[0] == "hu103":
        if r[1] == "10" and r[5] == "46.2":
            r[5] = "45.2"
            r[15] = ("设计排量1.8-0.8m3/min变排量; 2026-08-29校准:46.2->45.2(20313.doc 5.2.4(4)排量表46.2系笔误,"
                     "同文替浆浆柱表45.2/7.1.4顶替量计算理论总量90.2=7+23+15+45.2/7.2工艺流程表20+10+10+5.2三处自洽)。")
            h_c = True
            log.append(("pumping_schedule.csv", "step10 替钻井液(二段) volume", "46.2", "45.2(20313三处自洽,46.2判笔误)"))
        if r[1] == "13" and r[2] == "注水泥(实际)" and r[5] == "":
            r[5] = "83"
            r[15] = ("实际排量0.6m3/min; 泵压16-8MPa; 2026-08-29补充:实际注入量分相领浆22/中间浆35/尾浆26"
                     "(20314施工记录表'实际注入量'm3栏,20315施工流水账同记35),合计83m³;原volume空缺。")
            h_d = True
            log.append(("pumping_schedule.csv", "step13 注水泥(实际) volume", "(空)", "83(分相22/35/26,20314实际注入量)"))
if not h_c: print("WARN: pumping step10 未命中")
if not h_d: print("WARN: pumping step13 未命中")
save("pumping_schedule.csv", rows)

# ---- fluid_properties 防重复检查 ----
rows = load("fluid_properties.csv")
names = [r[1] for r in rows if r]
dup = names.count("驱油隔离液(实测流变)")
print("fluid_properties 隔离液(实测流变)行数:", dup, "(应为1)")
print("fluid_properties 中置液 notes 尾部:", [r[16][-40:] for r in rows if r and r[1] == "中置液"][0])

for t in log:
    print(" -", t[0], "|", t[1], "|", t[2], "->", t[3])
