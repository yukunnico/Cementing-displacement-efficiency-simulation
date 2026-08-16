# -*- coding: utf-8 -*-
import fitz, os
p = r"D:\users\desktop\research\控压固井项目\项目相关数据资料（未重新命名）\甲方完成数据\HT1-004井油层尾管上交资料1\HT1-004井油层尾管上交资料（甲方）1\化验报告\HT1-004井油层尾管第三方检测报告.pdf"
outdir = r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\ht1_004_呼1-004\_work\3rdparty_pages"
os.makedirs(outdir, exist_ok=True)
doc = fitz.open(p)
mat = fitz.Matrix(150/72, 150/72)
for i, page in enumerate(doc):
    pix = page.get_pixmap(matrix=mat)
    fp = os.path.join(outdir, f"3rdparty_p{i+1:02d}.png")
    pix.save(fp)
    print("saved", fp, pix.width, pix.height)
