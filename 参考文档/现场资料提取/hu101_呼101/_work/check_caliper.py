# -*- coding: utf-8 -*-
"""重建 2011113.doc 电测井径/井斜表并与提取包 CSV 逐点比对。"""
import csv, re, sys
from pathlib import Path

WORK = Path(r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\hu101_呼101\_work")
PKG = Path(r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\hu101_呼101")

text = (WORK / "doc2txt" / "2011113.txt").read_text(encoding="utf-8")
lines = text.splitlines()

# 定位电测表区间
start = next(i for i, l in enumerate(lines) if "三、电测井径" in l)
end = next(i for i, l in enumerate(lines) if "四、钻井复杂情况" in l)

# 行模式：|序号|底深|井径|井斜|方位| 重复两组（左右列）
row_re = re.compile(r"^\|(\s*\d+)?\s*\|(\s*\d+)?\s*\|([\d.]+)?\s*\|([\d.]+)?\s*\|([\d.]+)?\s*\|(\s*\d+)?\s*\|(\s*\d+)?\s*\|([\d.]+)?\s*\|([\d.]+)?\s*\|([\d.]+)?\s*\|")
cont_re = re.compile(r"^\|(?:\s*\|){10}$")

pts = {}  # md -> [cal, inc, azi]
pending = None  # 待合并的断行
for l in lines[start:end]:
    m = row_re.match(l)
    if m:
        g = m.groups()
        for base in (0, 5):
            idx, md, cal, inc, azi = g[base:base+5]
            if idx and md and cal:
                md = int(md)
                pts[md] = [float(cal), float(inc), float(azi)]
                pending = (md, 1)  # 记录最后写入的行（先左后右），断行归属最后一个
    elif cont_re.match(l) and pending:
        # 断行 "5" 归属最后一个完整行：检查该行原始断行位置
        # 原始表格单元宽 6；若某数值列在完整行长度<6 且下一行该列为5，则末位为5
        cols = [c.strip() for c in l.strip("|").split("|")]
        md, _ = pending
        if md in pts:
            for ci, cv in enumerate(cols):
                if cv == "5":
                    # 该列对应 pts 索引：列0=序号,1=底深,2=井径,3=井斜,4=方位（左组）
                    if ci in (2, 3, 4):
                        pts[md][ci-2] = float(f"{pts[md][ci-2]:.2f}" + "5")
        pending = None

# 生成重建表
print(f"重建点数: {len(pts)}  深度范围: {min(pts)}-{max(pts)}")
mds = sorted(pts)

def read_csv(p):
    with p.open(encoding="utf-8-sig") as f:
        return {float(r["md_m"]): r for r in csv.DictReader(f)}

cal_csv = read_csv(PKG / "caliper_profile.csv")
inc_csv = read_csv(PKG / "inclination_profile.csv")
print(f"CSV 井径点数: {len(cal_csv)}  深度 {min(cal_csv)}-{max(cal_csv)}")
print(f"CSV 井斜点数: {len(inc_csv)}  深度 {min(inc_csv)}-{max(inc_csv)}")

diffs = []
for md in mds:
    cal, inc, azi = pts[md]
    tag = []
    if md in cal_csv:
        c = float(cal_csv[md]["caliper_mm"])
        if abs(c - round(cal, 2)) > 0.011:
            tag.append(f"井径 CSV={c} vs 0708={cal}")
    else:
        tag.append("CSV缺井径点")
    if md in inc_csv:
        v = float(inc_csv[md]["inclination_deg"])
        if abs(v - round(inc, 2)) > 0.011:
            tag.append(f"井斜 CSV={v} vs 0708={inc}")
    else:
        tag.append("CSV缺井斜点")
    if tag:
        diffs.append((md, cal, inc, "; ".join(tag)))

print(f"\n差异点数: {len(diffs)}")
for md, cal, inc, t in diffs:
    print(f"  md={md} 0708井径={cal} 井斜={inc} | {t}")

# CSV 有而 0708 无的点
extra = [md for md in sorted(set(list(cal_csv) + list(inc_csv))) if md not in pts]
print(f"\nCSV 有而 0708 电测表无的点: {extra}")

# 关键锚点打印
print("\n关键锚点（0708 重建值）:")
for md in (5720, 6060, 6340, 6796, 6800, 6820, 6840, 7048, 7420, 7440, 7500, 7700, 7820, 7840, 7868):
    near = md if md in pts else min(pts, key=lambda x: abs(x - md))
    print(f"  {md}: 实际md={near} 井径={pts[near][0]} 井斜={pts[near][1]}")

# 井斜最大值
mx = max(pts.items(), key=lambda kv: kv[1][1])
print(f"\n井斜最大: md={mx[0]} inc={mx[1][1]}")
