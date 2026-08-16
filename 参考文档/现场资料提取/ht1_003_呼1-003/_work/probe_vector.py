# -*- coding: utf-8 -*-
"""HT1-003 矢量探测：检查两个 PDF 的一/二界面解释填充色带"""
import fitz, os

SRC_DIR = r"D:\users\desktop\research\控压固井项目\项目相关数据资料（未重新命名）\HT1-003、HT1-004井固井声幅图"
FILES = [
    "HT1-003_固井质量测井评价图_3500-7595_20260622_完井.pdf",
    "HT1-003_固井质量测井评价图_15-3500_20260622_完井.pdf",
]

for fn in FILES:
    print("=" * 70)
    print("FILE:", fn)
    doc = fitz.open(os.path.join(SRC_DIR, fn))
    page = doc[0]
    pw, ph = page.rect.width, page.rect.height
    print(f"  page: {pw:.1f} x {ph:.1f} pt")
    # 文字层概况
    txt = page.get_text()
    print(f"  text chars: {len(txt)}")
    draws = page.get_drawings()
    print(f"  drawings: {len(draws)}")
    # 统计填充对象
    fills = [d for d in draws if d.get('type') == 'f' and d.get('fill')]
    print(f"  filled paths: {len(fills)}")
    # 颜色直方图（近似到 0.05）
    from collections import Counter
    ccount = Counter()
    for d in fills:
        c = tuple(round(x, 2) for x in d['fill'])
        ccount[c] += 1
    print("  fill colors (count):")
    for c, n in ccount.most_common(15):
        print(f"    {c}: {n}")
    # 非白色填充的 x 分布
    xhist = Counter()
    for d in fills:
        c = tuple(round(x, 2) for x in d['fill'])
        if c in ((1.0, 1.0, 1.0), (1, 1, 1)):
            continue
        r = d['rect']
        cx = round((r.x0 + r.x1) / 2 / 10) * 10
        xhist[cx] += 1
    print("  non-white fill center-x histogram (bin=10pt):")
    for x in sorted(xhist):
        if xhist[x] > 2:
            print(f"    x~{x}: {xhist[x]}")
    doc.close()
