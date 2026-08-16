# -*- coding: utf-8 -*-
import fitz, os
SRC_DIR = r"D:\users\desktop\research\控压固井项目\项目相关数据资料（未重新命名）\HT1-003、HT1-004井固井声幅图"
FN = "HT1-003_固井质量测井评价图_3500-7595_20260622_完井.pdf"
doc = fitz.open(os.path.join(SRC_DIR, FN))
page = doc[0]
draws = page.get_drawings()
# 找 '胶结良好' 文字 y
td = page.get_text('dict')
for blk in td['blocks']:
    if blk['type'] != 0: continue
    for line in blk['lines']:
        txt = "".join(s['text'] for s in line['spans']).strip()
        if txt in ('胶结良好','胶结中等','胶结差','空套管'):
            b = line['bbox']
            print(f"text='{txt}' bbox x=({b[0]:.1f},{b[2]:.1f}) y=({b[1]:.1f},{b[3]:.1f})")
print("--- x 140-290, y 665-700 所有绘制 ---")
for d in draws:
    r = d['rect']
    if r.x0 < 290 and r.x1 > 140 and r.y0 < 700 and r.y1 > 665:
        t = d.get('type'); c = d.get('fill'); 
        print(f"  type={t} fill={tuple(round(x,2) for x in c) if c else None} rect=({r.x0:.1f},{r.y0:.1f},{r.x1:.1f},{r.y1:.1f}) w={r.width:.1f}")
doc.close()
