# -*- coding: utf-8 -*-
"""查黄色填充（油气水层段标注）位置与 x~230 列内容"""
import fitz, os, re
SRC_DIR = r"D:\users\desktop\research\控压固井项目\项目相关数据资料（未重新命名）\HT1-003、HT1-004井固井声幅图"
FN = "HT1-003_固井质量测井评价图_3500-7595_20260622_完井.pdf"
YELLOW = (1.0, 1.0, 0.0)
def close(c, ref, tol=0.04):
    return all(abs(a-b) <= tol for a,b in zip(c, ref))
doc = fitz.open(os.path.join(SRC_DIR, FN))
page = doc[0]
td = page.get_text('dict')
pts = []
for blk in td['blocks']:
    if blk['type'] != 0: continue
    for line in blk['lines']:
        for span in line['spans']:
            t = span['text'].strip()
            if re.fullmatch(r'\d{4}', t) and span['bbox'][0] < 70:
                pts.append((int(t), (span['bbox'][1]+span['bbox'][3])/2))
pts.sort()
ys = [p[0] for p in pts]; xs = [p[1] for p in pts]
k = (ys[-1]-ys[0])/(xs[-1]-xs[0]); b = ys[0] - k*xs[0]
def y2d(y): return k*y + b

print("--- 黄色填充（油气水层段标注）---")
for it in page.get_drawings():
    if it.get('type') != 'f' or not it.get('fill'): continue
    c = it['fill']; r = it['rect']
    if close(c, YELLOW):
        print(f"  rect=({r.x0:.1f},{r.y0:.1f},{r.x1:.1f},{r.y1:.1f}) -> depth {y2d(r.y1):.1f}-{y2d(r.y0):.1f}m w={r.width:.1f}")

print("\n--- x 225-260 列的文字内容（解释结论列/深度数字）---")
for blk in td['blocks']:
    if blk['type'] != 0: continue
    for line in blk['lines']:
        for span in line['spans']:
            t = span['text'].strip()
            b = span['bbox']
            if 220 < b[0] < 270 and b[1] > 900:
                if re.fullmatch(r'\d{2,4}', t) or t in ('解释','结论'):
                    print(f"  '{t}' at x={b[0]:.0f} y={b[1]:.0f} -> depth ~{y2d((b[1]+b[3])/2):.0f}m")
doc.close()
