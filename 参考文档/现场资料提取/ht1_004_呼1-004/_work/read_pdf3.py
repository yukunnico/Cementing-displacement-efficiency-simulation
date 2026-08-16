# -*- coding: utf-8 -*-
import fitz
p=r"D:\users\desktop\research\控压固井项目\项目相关数据资料（未重新命名）\甲方完成数据\HT1-004井油层尾管上交资料1\HT1-004井油层尾管上交资料（甲方）1\化验报告\HT1-004井油层尾管第三方检测报告.pdf"
doc=fitz.open(p)
print("pages=",doc.page_count)
total=0
for i in range(doc.page_count):
    t=doc[i].get_text().strip()
    total+=len(t)
    print(f"--p{i} chars={len(t)}")
    if t:
        print(t[:500])
print("TOTAL_CHARS=",total)
