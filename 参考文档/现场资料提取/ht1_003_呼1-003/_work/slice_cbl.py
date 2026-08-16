# -*- coding: utf-8 -*-
"""HT1-003 主评价图整条切片（6000px 高，200px 重叠），fitz clip 分块渲染。
以 pt 为单位：6000px@120dpi = 3600pt；重叠 200px = 120pt。"""
import fitz, os

SRC = r"D:\users\desktop\research\控压固井项目\项目相关数据资料（未重新命名）\HT1-003、HT1-004井固井声幅图\HT1-003_固井质量测井评价图_3500-7595_20260622_完井.pdf"
OUT = r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\ht1_003_呼1-003\cbl_pages"

DPI = 120
PX2PT = 72.0 / DPI          # 0.6
SLICE_H_PT = 6000 * PX2PT   # 3600 pt
OVERLAP_PT = 200 * PX2PT    # 120 pt
STEP_PT = SLICE_H_PT - OVERLAP_PT

doc = fitz.open(SRC)
page = doc[0]
pw, ph = page.rect.width, page.rect.height
print(f"page size pts: {pw:.1f} x {ph:.1f}  (= {ph*DPI/72:.0f} px @{DPI}dpi)")

# 清除旧切片
for f in os.listdir(OUT):
    if f.startswith("slice_") and f.endswith(".png"):
        os.remove(os.path.join(OUT, f))

slices = []
y0 = 0.0
while y0 < ph - 1:
    y1 = min(y0 + SLICE_H_PT, ph)
    slices.append((y0, y1))
    if y1 >= ph - 1:
        break
    y0 += STEP_PT

print(f"total slices: {len(slices)}")
os.makedirs(OUT, exist_ok=True)
for i, (y0, y1) in enumerate(slices):
    clip = fitz.Rect(0, y0, pw, y1)
    pix = page.get_pixmap(dpi=DPI, clip=clip, alpha=False)
    fn = os.path.join(OUT, f"slice_{i:02d}.png")
    pix.save(fn)
    print(f"slice_{i:02d}: y_pt[{y0:.0f},{y1:.0f}] = px[{y0/PX2PT:.0f},{y1/PX2PT:.0f}] size={pix.width}x{pix.height}")
doc.close()
print("done")
