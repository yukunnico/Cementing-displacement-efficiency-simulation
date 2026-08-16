# -*- coding: utf-8 -*-
import fitz, os
SRC_DIR = r"D:\users\desktop\research\控压固井项目\项目相关数据资料（未重新命名）\HT1-003、HT1-004井固声幅图"
print("SRC_DIR repr:", repr(SRC_DIR))
print("EXISTS dir:", os.path.exists(SRC_DIR))
FN = "HT1-003_固井质量测井评价图_3500-7595_20260622_完井.pdf"
full = os.path.join(SRC_DIR, FN)
print("EXISTS file:", os.path.exists(full))
doc = fitz.open(full)
page = doc[0]
draws = page.get_drawings()
print("--- 图例区 y 680-730 所有填充 ---")
for d in draws:
    if d.get('type') != 'f' or not d.get('fill'): continue
    r = d['rect']
    if r.y0 > 680 and r.y1 < 730:
        c = d['fill']
        print(f"  fill={tuple(round(x,2) for x in c)} rect=({r.x0:.1f},{r.y0:.1f},{r.x1:.1f},{r.y1:.1f}) w={r.width:.1f}")
print("--- 图例区 y 740-810 所有填充 ---")
for d in draws:
    if d.get('type') != 'f' or not d.get('fill'): continue
    r = d['rect']
    if r.y0 > 740 and r.y1 < 810:
        c = d['fill']
        print(f"  fill={tuple(round(x,2) for x in c)} rect=({r.x0:.1f},{r.y0:.1f},{r.x1:.1f},{r.y1:.1f}) w={r.width:.1f}")
doc.close()
