# -*- coding: utf-8 -*-
"""2026-08-29 hu103 提取包 CSV 校准修正（仅确证错误，保持原格式 utf-8-sig/逗号/行序）。"""
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

# ---- A. cbl_evaluation.csv 综合行：6.05%(简单平均误标加权) -> 2.12%(长度加权, PDF测量井段) ----
rows = load("cbl_evaluation.csv")
for r in rows:
    if len(r) > 4 and r[3] == "6.05" and r[1] == "呼103":
        old = ",".join(r)
        r[3] = "2.12"
        r[2] = "7712.0"  # 综合评价上界=人工井底/阻位7712.41(CBL实测终点), 原7770为完钻井深口径
        r[6] = ("168.3mm段(0.04%)+139.7mm段(12.06%)长度加权：按100515/100516测量井段5540-7330.6m(1790.6m)与"
                "7338-7712m(374m)，综合=(0.04*1790.6+12.06*374)/2164.6≈2.12%。"
                "2026-08-29校准：原值6.05%系两段简单平均(0.04+12.06)/2误标为加权，其原注公式实算2.16%，"
                "现按0708 PDF测量井段改为长度加权2.12%。")
        r[9] = "综合计算（168.3mm+139.7mm长度加权，2026-08-29校准）"
        r[13] = "整段尾管综合合格率约2.12%——模型顶替效率应与此对比验证。"
        log.append(("cbl_evaluation.csv", "综合行 pass_rate", "6.05(简单平均,原公式实算2.16)", "2.12(长度加权,测量井段1790.6m/374m); 区间顶7770->7712"))
        break
save("cbl_evaluation.csv", rows)

# ---- B. rheometer_readings.csv：203111表7"隔离液"行误标为"平衡液" ----
rows = load("rheometer_readings.csv")
for r in rows:
    if len(r) > 4 and r[1] == "平衡液" and r[2] == "140":
        r[1] = "驱油隔离液"
        r[11] = ("140C↘93C驱油隔离液六速; n=0.54 K=2.12Pa.s^n(203111表7)。"
                 "2026-08-29校准:原行名'平衡液'系提取错误,203111.docx表7该行明确为'隔离液';六速数值不变。")
        log.append(("rheometer_readings.csv", "140C行 fluid_name", "平衡液", "驱油隔离液(n=0.54/K=2.12); 数值128/74/54/33/6/5不变"))
        break
save("rheometer_readings.csv", rows)

# ---- C/D. pumping_schedule.csv：46.2->45.2；实际水泥浆步补分相体积 ----
rows = load("pumping_schedule.csv")
for r in rows:
    if len(r) > 10 and r[0] == "hu103":
        if r[1] == "10" and r[6] == "46.2":
            r[6] = "45.2"
            r[15] = ("设计排量1.8-0.8m3/min变排量; 2026-08-29校准:46.2->45.2(20313.doc 5.2.4(4)排量表46.2系笔误,"
                     "同文替浆浆柱表45.2/7.1.4顶替量计算理论总量90.2=7+23+15+45.2/7.2工艺流程表20+10+10+5.2三处自洽)。")
            log.append(("pumping_schedule.csv", "step10 替钻井液(二段) volume", "46.2", "45.2(20313内部45.2/90.2/工艺流程三处自洽,46.2判笔误)"))
        if r[1] == "13" and r[5] == "注水泥(实际)" and (r[6] == "" or r[6] is None):
            r[6] = "83"
            r[15] = ("实际排量0.6m3/min; 泵压16-8MPa; 2026-08-29补充:实际注入量分相领浆22/中间浆35/尾浆26(20314施工记录表"
                     "'实际注入量'm3栏,20315施工流水账同记35),合计83m³;原volume空缺。")
            log.append(("pumping_schedule.csv", "step13 注水泥(实际) volume", "(空)", "83(分相22/35/26,20314实际注入量)"))
save("pumping_schedule.csv", rows)

# ---- E. fluid_properties.csv：中置液notes扩注分歧；新增隔离液实测流变行 ----
rows = load("fluid_properties.csv")
out = []
for r in rows:
    if len(r) > 3 and r[0] == "hu103" and r[1] == "中置液":
        r[16] = (r[16] + "; 2026-08-29校准注:密度存在0708内部口径分歧——7.1.3后置液表/9.2材料清单/203111表3记2.00,"
                 "5.2.4(3)替浆浆柱表记1.95(与配方195%重晶石对应),实际施工记1.90(20314/20315);本表保守维持1.95,未改值。")
    out.append(r)
new_row = ["hu103", "驱油隔离液(实测流变)", "spacer", "2.00", "幂律", "", "", "", "2.12", "0.54", "", "140↘93",
           "design", r"2\203\2031\20311\203111.docx", "Table7流变性能", "high",
           "2026-08-29补充:203111表7隔离液实测流变 n=0.54 K=2.12Pa.s^n(140C↘93C); 六速128/74/54/33/6/5; "
           "设计密度2.00(隔离液1)/1.95(隔离液2); 原fluid_properties漏记此实测行。"]
out.append(new_row)
log.append(("fluid_properties.csv", "隔离液实测流变行", "(缺)", "新增 n=0.54 K=2.12 @140↘93C(203111表7); 中置液行notes扩注口径分歧(值未改)"))
save("fluid_properties.csv", out)

print("修正完成:")
for t in log:
    print(" -", t[0], "|", t[1], "|", t[2], "->", t[3])

# 校验回读
for name in ["cbl_evaluation.csv", "rheometer_readings.csv", "pumping_schedule.csv", "fluid_properties.csv"]:
    rows = load(name)
    print(name, "行数", len(rows), "首列尾行:", rows[-1][:3])
