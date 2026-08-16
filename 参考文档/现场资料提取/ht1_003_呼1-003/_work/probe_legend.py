# -*- coding: utf-8 -*-
"""验证图例区色块：找 胶结良好/中等/差/空套管 文字 bbox，检查其附近填充颜色"""
import fitz, os

SRC_DIR = r"D:\users\desktop\research\控压固井项目\项目相关数据资料（未重新命名）\HT1-003、HT1-004井固井声幅图"
FN = "HT1-003_固井质量测井评价图_3500-7595_20260622_完井.pdf"
doc = fitz.open(os.path.join(SRC_DIR, FN))
page = doc[0]

# 找图例文字
td = page.get_text('dict')
for blk in td['blocks']:
    if blk['type'] != 0: continue
    for line in blk['lines']:
        txt = "".join(s['text'] for s in line['spans']).strip()
        if any(k in txt for k in ['胶结良好','胶结中等','胶结差','空套管','图','例']):
            b = line['bbox']
            print(f"text='{txt}'  bbox=({b[0]:.0f},{b[1]:.0f},{b[2]:.0f},{b[3]:.0f})")

# 在文字 y 附近找色块填充（图例色块一般在文字左边 x<300）
draws = page.get_drawings()
def close(c, ref, tol=0.04):
    return all(abs(a-b) <= tol for a,b in zip(c, ref))
print("\n--- 图例附近(y 800-920, x 250-500)的非白填充 ---")
for d in draws:
    if d.get('type') != 'f' or not d.get('fill'): continue
    r = d['rect']
    if r.y0 > 780 and r.y1 < 940 and 250 < r.x0 < 520:
        c = d['fill']
        if close(c, (1,1,1)): continue
        print(f"  fill={tuple(round(x,2) for x in c)} rect=({r.x0:.0f},{r.y0:.0f},{r.x1:.0f},{r.y1:.0f}) w={r.width:.1f}")

# 一界面道内所有非青/绿填充颜色（排查混淆）
print("\n--- 一界面道 688.5-719.7 内所有填充颜色 ---")
from collections import Counter
cc = Counter()
for d in draws:
    if d.get('type') != 'f' or not d.get('fill'): continue
    r = d['rect']
    cx = (r.x0+r.x1)/2
    if 688.5-2 <= cx <= 719.7+2 and r.y1 > 915:
        cc[tuple(round(x,2) for x in d['fill'])] += 1
for c, n in cc.most_common():
    print(f"  {c}: {n}")
doc.close()
