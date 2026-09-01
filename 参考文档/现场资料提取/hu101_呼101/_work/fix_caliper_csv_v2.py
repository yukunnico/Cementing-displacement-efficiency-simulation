# -*- coding: utf-8 -*-
"""基于 2011116.xls Sheet2（电测数据电子版，109 测点）全量重写 caliper/inclination CSV。

v1 修正（依据 antiword 重建 2011113.doc）已处理右列错位，但 antiword 丢失了
序号1（5700）与序号56（6800）两行；Sheet2 电子版 109 点完整且与 antiword 重建
逐点一致（6820/6840/7440/7840/7868 等互证），故以 Sheet2 为权威源全量重写。
"""
import csv, shutil
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path
import xlrd

WORK = Path(r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\hu101_呼101\_work")
PKG = Path(r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\hu101_呼101")

wb = xlrd.open_workbook(str(WORK / "2011116.xls"))
sh = wb.sheet_by_name("Sheet2")
pts = {}  # md -> (cal, inc, azi)
for r in range(sh.nrows):
    v = [sh.cell_value(r, c) for c in range(10)]
    for base in (0, 5):
        idx, md, cal, inc, azi = v[base:base+5]
        if isinstance(idx, float) and idx > 0 and isinstance(md, float) and md > 0:
            pts[int(md)] = (float(cal), float(inc), float(azi))
print(f"Sheet2 点数: {len(pts)}  范围 {min(pts)}-{max(pts)}")
assert len(pts) == 109 and 5700 in pts and 6800 in pts and 7868 in pts

def r2(x):
    return str(Decimal(str(x)).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP))

def g(x):
    s = r2(x)
    return s.rstrip("0").rstrip(".") if "." in s else s

# ---- caliper ----
cal_path = PKG / "caliper_profile.csv"
shutil.copy(cal_path, WORK / "caliper_profile.csv.v1bak")
rows = list(csv.DictReader(cal_path.open(encoding="utf-8-sig")))
fields = list(rows[0].keys())
out = []
n_change = n_add = 0
existing = {float(r["md_m"]): r for r in rows}
for md in sorted(pts):
    cal, inc, azi = pts[md]
    row = existing.get(float(md))
    if row is None:
        row = {k: "" for k in fields}
        row["well_id"] = "hu101"
        row["md_m"] = str(md)
        row["borehole_nominal_diameter_mm"] = "215.9"
        row["data_type"] = "field_measured"
        row["source_file"] = "2/201/2011/20111/201111/2011113.doc + 2011116.xls Sheet2"
        row["source_location"] = "三、电测井径、井斜、方位"
        row["confidence"] = "high"
        row["notes"] = "2026-08-29 校准：按 2011116.xls Sheet2 电子版全量重建（109 点）"
        n_add += 1
    new_cal = g(cal)
    if row["caliper_mm"] != new_cal:
        n_change += 1
    row["caliper_mm"] = new_cal
    row["enlargement_ratio"] = g((float(new_cal) - 215.9) / 215.9 * 100)
    out.append(row)
with cal_path.open("w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(out)
print(f"caliper_profile.csv: {len(out)} 行（新增 {n_add}，改值 {n_change}）")

# ---- inclination ----
inc_path = PKG / "inclination_profile.csv"
shutil.copy(inc_path, WORK / "inclination_profile.csv.v1bak")
rows = list(csv.DictReader(inc_path.open(encoding="utf-8-sig")))
fields = list(rows[0].keys())
out = []
n_change = n_add = 0
existing = {float(r["md_m"]): r for r in rows}
for md in sorted(pts):
    cal, inc, azi = pts[md]
    row = existing.get(float(md))
    if row is None:
        row = {k: "" for k in fields}
        row["well_id"] = "hu101"
        row["md_m"] = str(md)
        row["data_type"] = "field_measured"
        row["source_file"] = "2/201/2011/20111/201111/2011113.doc + 2011116.xls Sheet2"
        row["source_location"] = "三、电测井径、井斜、方位"
        row["confidence"] = "high"
        row["notes"] = "2026-08-29 校准：按 2011116.xls Sheet2 电子版全量重建（109 点）"
        n_add += 1
    new_inc, new_azi = g(inc), g(azi)
    if row["inclination_deg"] != new_inc or row["azimuth_deg"] != new_azi:
        n_change += 1
    row["inclination_deg"] = new_inc
    row["azimuth_deg"] = new_azi
    row["tvd_m"] = g(md * 0.999)  # 提取包原近似口径
    out.append(row)
with inc_path.open("w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(out)
print(f"inclination_profile.csv: {len(out)} 行（新增 {n_add}，改值 {n_change}）")

# 抽样验证
for md in (5700, 5720, 6800, 6820, 7440, 7840, 7868):
    print(md, pts[md])
