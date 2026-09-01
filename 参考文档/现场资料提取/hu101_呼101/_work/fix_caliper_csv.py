# -*- coding: utf-8 -*-
"""修正 caliper_profile.csv / inclination_profile.csv 右列错位。

0708 依据：2011113.doc《固井作业史》三、电测井径井斜方位表。
原 CSV 把右列 53 个测点（序号57-109，md=6820-7868）的值整体前移一个测点存储
（6800←6820 的值，6820←6840 的值，…，7840←7868 的值），再复制出 7868 行。
修正：删除 6800 行（0708 无此测点），6820-7868 按 0708 重建值重写。
"""
import csv, re, shutil
from pathlib import Path

WORK = Path(r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\hu101_呼101\_work")
PKG = Path(r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\hu101_呼101")

text = (WORK / "doc2txt" / "2011113.txt").read_text(encoding="utf-8")
lines = text.splitlines()
start = next(i for i, l in enumerate(lines) if "三、电测井径" in l)
end = next(i for i, l in enumerate(lines) if "四、钻井复杂情况" in l)

row_re = re.compile(r"^\|(\s*\d+)?\s*\|(\s*\d+)?\s*\|([\d.]+)?\s*\|([\d.]+)?\s*\|([\d.]+)?\s*\|(\s*\d+)?\s*\|(\s*\d+)?\s*\|([\d.]+)?\s*\|([\d.]+)?\s*\|([\d.]+)?\s*\|")
cont_re = re.compile(r"^\|(?:\s*\|){10}$")

pts = {}
pending = None  # (左md 或 None, 右md 或 None)
for l in lines[start:end]:
    m = row_re.match(l)
    if m:
        g = m.groups()
        left = right = None
        for base, slot in ((0, "L"), (5, "R")):
            idx, md, cal, inc, azi = g[base:base+5]
            if idx and md and cal:
                md = int(md)
                pts[md] = [float(cal), float(inc), float(azi)]
                if slot == "L":
                    left = md
                else:
                    right = md
        pending = (left, right)
    elif cont_re.match(l) and pending:
        cols = [c.strip() for c in l.strip("|").split("|")]
        for ci, cv in enumerate(cols):
            if cv == "5" and ci in (2, 7):  # 仅井径列有断行
                md = pending[0] if ci == 2 else pending[1]
                if md is not None and md in pts:
                    base = f"{pts[md][ci//5-0 if False else 0]:.2f}"
                    pts[md][0] = float(base + "5")
        pending = None

mds = sorted(pts)
print(f"重建点数 {len(mds)}: {mds[0]}..{mds[-1]}")

# ---- 修正 caliper_profile.csv ----
cal_path = PKG / "caliper_profile.csv"
shutil.copy(cal_path, WORK / "caliper_profile.csv.bak")
rows = list(csv.DictReader(cal_path.open(encoding="utf-8-sig")))
fields = list(rows[0].keys())
out = []
changed = 0
for r in rows:
    md = float(r["md_m"])
    if md == 6800:
        changed += 1
        print(f"删除 6800 行（0708 无此测点；原值 {r['caliper_mm']} 实为 6820 测点值）")
        continue
    if md >= 6820:
        if md not in pts:
            raise SystemExit(f"重建表缺 {md}")
        cal, inc, azi = pts[md]
        new_cal = round(cal, 2)
        if abs(float(r["caliper_mm"]) - new_cal) > 1e-9:
            changed += 1
        r["caliper_mm"] = f"{new_cal:g}"
        r["enlargement_ratio"] = f"{round((new_cal - 215.9) / 215.9 * 100, 2):g}"
    out.append(r)

with cal_path.open("w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(out)
print(f"caliper_profile.csv: {len(out)} 行, 修改/删除 {changed} 处")

# ---- 修正 inclination_profile.csv ----
inc_path = PKG / "inclination_profile.csv"
shutil.copy(inc_path, WORK / "inclination_profile.csv.bak")
rows = list(csv.DictReader(inc_path.open(encoding="utf-8-sig")))
fields = list(rows[0].keys())
out = []
changed = 0
for r in rows:
    md = float(r["md_m"])
    if md == 6800:
        changed += 1
        print(f"删除 6800 行（0708 无此测点；原值 inc={r['inclination_deg']})")
        continue
    if md >= 6820:
        if md not in pts:
            raise SystemExit(f"重建表缺 {md}")
        cal, inc, azi = pts[md]
        new_inc = round(inc, 2)
        new_azi = round(azi, 2)
        if abs(float(r["inclination_deg"]) - new_inc) > 1e-9:
            changed += 1
        r["inclination_deg"] = f"{new_inc:g}"
        r["azimuth_deg"] = f"{new_azi:g}"
        r["tvd_m"] = f"{round(md * 0.999, 2):g}"  # 提取包原近似口径 md*0.999
    out.append(r)

with inc_path.open("w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(out)
print(f"inclination_profile.csv: {len(out)} 行, 修改/删除 {changed} 处")

# 打印修正后右列抽样
for md in (6820, 6840, 6900, 7440, 7700, 7840, 7868):
    print(md, pts[md])
