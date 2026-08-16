# -*- coding: utf-8 -*-
"""HT1-004 CBL 评价图矢量数字化：
1. 深度刻度 → 深度-y 标定
2. 一界面/二界面解释道颜色填充 → 分段质量
3. 解释结论列黄色条带+标注数值
4. 声幅道黑描边曲线 → 声幅值（试验）
"""
import fitz, os, re, json
from collections import defaultdict

SRC_DIR = r"D:\users\desktop\research\控压固井项目\项目相关数据资料（未重新命名）\HT1-003、HT1-004井固井声幅图"
OUT = r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\ht1_004_呼1-004\cbl_digitization\_work"

# 颜色定义（与图例实测一致）
CYAN = (0.0, 0.6, 0.6)     # 胶结良好(实心)/胶结中等(条纹)
GREEN = (0.0, 0.5, 0.0)    # 胶结差(条纹)
YELLOW = (1.0, 1.0, 0.0)   # 解释结论列标注条带
RED = (1.0, 0.0, 0.0)
BLUE = (0.0, 0.0, 1.0)

def close(c, ref, tol=0.03):
    return all(abs(a-b) <= tol for a, b in zip(c, ref))

def classify_fill(c):
    """按图例颜色分类填充色"""
    if close(c, CYAN): return 'cyan'
    if close(c, GREEN): return 'green'
    if close(c, YELLOW): return 'yellow'
    if close(c, (1,1,1)): return 'white'
    if close(c, (0,0,0)): return 'black'
    return 'other'

def analyze_pdf(fn, tag):
    doc = fitz.open(os.path.join(SRC_DIR, fn))
    page = doc[0]
    pw, ph = page.rect.width, page.rect.height
    td = page.get_text('dict')

    # 1. 深度刻度
    ticks = []
    for blk in td['blocks']:
        if blk['type']!=0: continue
        for line in blk['lines']:
            for span in line['spans']:
                t = span['text'].strip()
                if re.fullmatch(r'\d{4}', t):
                    v = int(t)
                    b = span['bbox']
                    ticks.append((v, (b[1]+b[3])/2))
    # 取 x<70 的左侧深度刻度
    ticks_l = []
    for blk in td['blocks']:
        if blk['type']!=0: continue
        for line in blk['lines']:
            for span in line['spans']:
                t = span['text'].strip()
                if re.fullmatch(r'\d{4}', t):
                    v = int(t); b = span['bbox']
                    if b[0] < 70:
                        ticks_l.append((v, (b[1]+b[3])/2))
    ticks_l.sort()
    print(f"[{tag}] depth ticks(left): {[v for v,_ in ticks_l][:5]}...{ticks_l[-1] if ticks_l else ''} n={len(ticks_l)}")
    if len(ticks_l) >= 2:
        # 线性拟合
        xs = [y for _, y in ticks_l]; ys = [v for v, _ in ticks_l]
        k = (ys[-1]-ys[0])/(xs[-1]-xs[0])  # m per pt
        b0 = ys[0] - k*xs[0]
        print(f"[{tag}] calib: depth = {k:.5f}*y + {b0:.2f}   (y of tick {ys[0]}={xs[0]:.1f})")

    # 2. 一界面/二界面填充（按图例判断，这里动态找最右两个彩色道）
    #    先收集所有 y>740 的非白填充，按 x 聚类
    fills = []
    for it in page.get_drawings():
        if it.get('type')=='f':
            c = it.get('fill')
            r = it['rect']
            if r.y0 > 740 and r.y0 < 3700 and c and not close(c,(1,1,1)):
                fills.append((r, c))
    # x 聚类
    xgroups = defaultdict(list)
    for r, c in fills:
        xgroups[round(r.x0/10)].append((r, c))
    # 打印 x 位置分布（按填充条带数）
    print(f"[{tag}] 非白填充条带数 by x-bin(10pt):")
    for xbin in sorted(xgroups):
        g = xgroups[xbin]
        cls = defaultdict(int)
        for _, c in g:
            cls[classify_fill(c)] += 1
        x0 = min(r.x0 for r,_ in g); x1 = max(r.x1 for r,_ in g)
        print(f"   x[{x0:.0f}-{x1:.0f}] n={len(g)} {dict(cls)}")
    doc.close()
    return ticks_l

if __name__ == '__main__':
    JOBS = [
        ("HT1-004_固井质量测井评价图_7382-7581_20260721_完井.pdf", "tail"),
        ("HT1-004_固井质量测井评价图_5249-7384_20260721_完井.pdf", "mid"),
        ("HT1-004_固井质量测井评价图_11-5251_20260721_完井.pdf", "shallow"),
    ]
    for fn, tag in JOBS:
        analyze_pdf(fn, tag)
        print()
