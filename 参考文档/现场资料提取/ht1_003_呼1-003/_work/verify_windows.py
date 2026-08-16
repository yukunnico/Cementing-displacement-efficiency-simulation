# -*- coding: utf-8 -*-
"""复核关键窗口的一界面填充明细（原始矢量，不经过采样合并）"""
import fitz, os, re
SRC_DIR = r"D:\users\desktop\research\控压固井项目\项目相关数据资料（未重新命名）\HT1-003、HT1-004井固井声幅图"
FN = "HT1-003_固井质量测井评价图_3500-7595_20260622_完井.pdf"
CYAN = (0.0, 0.6, 0.6); GREEN = (0.0, 0.5, 0.0)
def close(c, ref, tol=0.04):
    return all(abs(a-b) <= tol for a,b in zip(c, ref))

doc = fitz.open(os.path.join(SRC_DIR, FN))
page = doc[0]
# 标定
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
print(f"calib k={k:.5f} b={b:.2f}")

def y2d(y): return k*y + b

WINDOWS = [(5240, 5360, "红线1窗口：重合段5286-5316m"), (7420, 7600, "红线2窗口：目的井段7442-7618m")]
for dmin, dmax, label in WINDOWS:
    print("\n" + "="*60)
    print(label, f"({dmin}-{dmax}m)")
    ymin = (dmin - b)/k; ymax = (dmax - b)/k
    # 一界面道 688.5-719.7 内所有填充，转深度
    items = []
    for it in page.get_drawings():
        if it.get('type') != 'f' or not it.get('fill'): continue
        c = it['fill']; r = it['rect']
        cx = (r.x0+r.x1)/2
        if 688.5-2 <= cx <= 719.7+2 and r.y1 > ymin-5 and r.y0 < ymax+5:
            if close(c, CYAN): cls = 'cyan'
            elif close(c, GREEN): cls = 'green'
            else: continue
            d0 = y2d(r.y0); d1 = y2d(r.y1)
            items.append((cls, d0, d1, r.width, r.y0, r.y1))
    items.sort(key=lambda t: t[4])
    for cls, d0, d1, w, yy0, yy1 in items:
        if d1 > dmin and d0 < dmax:
            label_c = {'cyan':'青','green':'绿'}[cls]
            solid = '实心' if (cls=='cyan' and w>15) else ('条纹' if cls=='cyan' else '细条')
            print(f"  {d0:7.1f}-{d1:7.1f}m  {label_c}{solid} w={w:.1f}")
doc.close()
