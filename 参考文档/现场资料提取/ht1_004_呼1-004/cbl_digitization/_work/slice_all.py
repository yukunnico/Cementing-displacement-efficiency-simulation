# -*- coding: utf-8 -*-
"""HT1-004 主评价图整条切片（6000px 高，200px 重叠），fitz clip 分块渲染。
以 pt 为单位：6000px@120dpi = 3600pt；重叠 200px = 120pt。
普通版 3 段全切；RCD 版只切顶部 1 片（图头）。"""
import fitz, os

SRC_DIR = r"D:\users\desktop\research\控压固井项目\项目相关数据资料（未重新命名）\HT1-003、HT1-004井固井声幅图"
OUT = r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\ht1_004_呼1-004\cbl_digitization\pages"
os.makedirs(OUT, exist_ok=True)

JOBS = [
    ("HT1-004_固井质量测井评价图_11-5251_20260721_完井.pdf", "seg11_5251"),
    ("HT1-004_固井质量测井评价图_5249-7384_20260721_完井.pdf", "seg5249_7384"),
    ("HT1-004_固井质量测井评价图_7382-7581_20260721_完井.pdf", "seg7382_7581"),
    ("HT1-004_固井质量测井评价图_11-5251_20260721_完井RCD.pdf", "rcd11_5251"),
    ("HT1-004_固井质量测井评价图_5249-7384_20260721_完井RCD.pdf", "rcd5249_7384"),
    ("HT1-004_固井质量测井评价图_7382-7581_20260721_完井RCD.pdf", "rcd7382_7581"),
]

DPI = 120
PX2PT = 72.0 / DPI          # 0.6
SLICE_H_PT = 6000 * PX2PT   # 3600 pt
OVERLAP_PT = 200 * PX2PT    # 120 pt
STEP_PT = SLICE_H_PT - OVERLAP_PT

for fn, tag in JOBS:
    src = os.path.join(SRC_DIR, fn)
    doc = fitz.open(src)
    page = doc[0]
    pw, ph = page.rect.width, page.rect.height
    print(f"=== {fn} page={pw:.0f}x{ph:.0f}pt = {ph*DPI/72:.0f}px@{DPI}dpi")

    # RCD 只取顶部 1 片（图头）
    n_max = 1 if tag.startswith("rcd") else 9999

    slices = []
    y0 = 0.0
    while y0 < ph - 1:
        y1 = min(y0 + SLICE_H_PT, ph)
        slices.append((y0, y1))
        if y1 >= ph - 1 or len(slices) >= n_max:
            break
        y0 += STEP_PT

    print(f"  slices: {len(slices)}")
    for i, (y0_, y1_) in enumerate(slices):
        clip = fitz.Rect(0, y0_, pw, y1_)
        pix = page.get_pixmap(dpi=DPI, clip=clip, alpha=False)
        fn_out = os.path.join(OUT, f"{tag}_slice_{i+1:02d}.png")
        pix.save(fn_out)
        print(f"  {tag}_slice_{i+1:02d}: y_pt[{y0_:.0f},{y1_:.0f}] px[{y0_/PX2PT:.0f},{y1_/PX2PT:.0f}] size={pix.width}x{pix.height}")
    doc.close()

print("done")
