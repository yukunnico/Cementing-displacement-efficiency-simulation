# -*- coding: utf-8 -*-
import fitz, sys
for p in [
 r"D:\users\desktop\research\控压固井项目\项目相关数据资料（未重新命名）\甲方完成数据\HT1-004井油层尾管上交资料1\HT1-004井油层尾管上交资料（甲方）1\XD80502队HT1-004井139.7+168.28㎜上扣扭矩.pdf",
 r"D:\users\desktop\research\控压固井项目\项目相关数据资料（未重新命名）\甲方完成数据\HT1-004井油层尾管上交资料1\HT1-004井油层尾管上交资料（甲方）1\化验报告\HT1-004井油层尾管第三方检测报告.pdf",
]:
    try:
        doc=fitz.open(p)
        print("\n===== PDF:", p.split("\\")[-1], "pages=",doc.page_count)
        total_chars=0
        for i in range(min(doc.page_count,8)):
            t=doc[i].get_text().strip()
            total_chars+=len(t)
            print(f"--p{i} chars={len(t)}")
            if t:
                print(t[:800])
        print("TOTAL_TEXT_CHARS(first8):", total_chars)
    except Exception as e:
        print("ERR", p, e)
