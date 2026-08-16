# -*- coding: utf-8 -*-
"""精确探测一/二界面色带的 x 边界与条纹宽度"""
import fitz, os
from collections import Counter

SRC_DIR = r"D:\users\desktop\research\控压固井项目\项目相关数据资料（未重新命名）\HT1-003、HT1-004井固声幅图"
SRC_DIR = r"D:\users\desktop\research\控压固井项目\项目相关数据资料（未重新命名）\HT1-003、HT1-004井固井声幅图"
FILES = [
    "HT1-003_固井质量测井评价图_3500-7595_20260622_完井.pdf",
    "HT1-003_固井质量测井评价图_15-3500_20260622_完井.pdf",
]

CYAN = (0.0, 0.6, 0.6); GREEN = (0.0, 0.5, 0.0)

def close(c, ref, tol=0.04):
    return all(abs(a-b) <= tol for a, b in zip(c, ref))

for fn in FILES:
    print("=" * 70)
    print("FILE:", fn)
    doc = fitz.open(os.path.join(SRC_DIR, fn))
    page = doc[0]
    # 青色/绿色填充的 x0/x1 分布
    cyan_x0 = Counter(); cyan_x1 = Counter(); green_x0 = Counter(); green_x1 = Counter()
    cyan_w = Counter(); green_w = Counter()
    for d in page.get_drawings():
        if d.get('type') != 'f' or not d.get('fill'):
            continue
        c = d['fill']; r = d['rect']
        if close(c, CYAN):
            cyan_x0[round(r.x0,1)] += 1; cyan_x1[round(r.x1,1)] += 1
            cyan_w[round(r.width,1)] += 1
        elif close(c, GREEN):
            green_x0[round(r.x0,1)] += 1; green_x1[round(r.x1,1)] += 1
            green_w[round(r.width,1)] += 1
    def show(name, cnt, top=12):
        print(f"  {name}:")
        for k, n in cnt.most_common(top):
            print(f"    {k}: {n}")
    show("CYAN x0", cyan_x0)
    show("CYAN x1", cyan_x1)
    show("CYAN width", cyan_w)
    show("GREEN x0", green_x0)
    show("GREEN x1", green_x1)
    show("GREEN width", green_w)
    doc.close()
