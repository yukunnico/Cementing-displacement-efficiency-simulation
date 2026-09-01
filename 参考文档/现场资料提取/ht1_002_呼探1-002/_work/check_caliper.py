# -*- coding: utf-8 -*-
"""全量核对提取包 64 点井径/井斜 vs 0708 施工设计 1.4.1/1.4.2 电测表（Word COM 导出 txt）。"""
import csv
import re
from pathlib import Path

TXT = Path(r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\ht1_002_呼探1-002\_work\doc_txt\wukai_HT1-002井139.7mm尾管完井固井施工设计（研究中心审核版审批版7.1）(1).txt")
CAL = Path(r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\ht1_002_呼探1-002\caliper_profile.csv")
INC = Path(r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\ht1_002_呼探1-002\inclination_profile.csv")

text = TXT.read_text(encoding="utf-8")
lines = [x.strip() for x in text.splitlines()]

def is_num(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def parse_caliper(lines, start, end):
    """井径表：找 深度行(56xx-75xx) 后紧邻 井径行(188-275)。"""
    pairs = {}
    for i in range(start, end - 1):
        a, b = lines[i], lines[i + 1]
        if is_num(a) and is_num(b):
            d, c = float(a), float(b)
            if 5600 <= d <= 7600 and 185 <= c <= 280 and c not in (float("inf"),):
                # 深度步长 25-45m 才算（30m 网格），排除其他数字误配
                pairs[d] = c
    return pairs

def parse_incl(lines, start, end):
    """井斜表：深度行(56xx-75xx) 后接 井斜(0-5) 再接 方位(0-360)。"""
    pairs = {}
    for i in range(start, end - 2):
        a, b, c = lines[i], lines[i + 1], lines[i + 2]
        if is_num(a) and is_num(b) and is_num(c):
            d, inc, azi = float(a), float(b), float(c)
            if 5600 <= d <= 7600 and 0 < inc <= 5.0 and 0 <= azi <= 360:
                pairs[d] = (inc, azi)
    return pairs

# 1.4.2 井径表范围（txt 行号 686-1050），1.4.1 井斜表范围（行 321-528，按标题定位）
cal_start = next(i for i, l in enumerate(lines) if "电测井径" in l)
inc_start = next(i for i, l in enumerate(lines) if "井斜、方位" in l)
cal_pairs = parse_caliper(lines, cal_start, cal_start + 400)
inc_pairs = parse_incl(lines, inc_start, inc_start + 260)

def read_csv(path, val_col):
    rows = {}
    with path.open("r", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            rows[float(r["md_m"])] = float(r[val_col])
    return rows

csv_cal = read_csv(CAL, "caliper_mm")
csv_inc = read_csv(INC, "inclination_deg")

print(f"设计表解析: 井径 {len(cal_pairs)} 点, 井斜 {len(inc_pairs)} 点")
print(f"CSV: 井径 {len(csv_cal)} 点, 井斜 {len(csv_inc)} 点")

print("\n--- 井径逐点核对（设计 -> CSV）---")
mismatch = 0
for d in sorted(cal_pairs):
    c0 = cal_pairs[d]
    c1 = csv_cal.get(d)
    flag = "" if (c1 is not None and abs(c1 - c0) < 0.05) else "  <== 不一致"
    if flag:
        mismatch += 1
        print(f"{d:.0f}: 设计 {c0} vs CSV {c1}{flag}")
print(f"井径不一致数: {mismatch}")
missing_cal = sorted(set(csv_cal) - set(cal_pairs))
print(f"CSV 多出/设计没有的深度: {missing_cal}")

print("\n--- 井斜逐点核对 ---")
mismatch_i = 0
for d in sorted(inc_pairs):
    i0, a0 = inc_pairs[d]
    i1 = csv_inc.get(d)
    flag = "" if (i1 is not None and abs(i1 - i0) < 0.05) else "  <== 不一致"
    if flag:
        mismatch_i += 1
        print(f"{d:.0f}: 设计 {i0} vs CSV {i1}{flag}")
print(f"井斜不一致数: {mismatch_i}")
