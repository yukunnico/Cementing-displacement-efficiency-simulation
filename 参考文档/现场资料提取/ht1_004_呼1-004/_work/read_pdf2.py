# -*- coding: utf-8 -*-
import fitz
p=r"D:\users\desktop\research\控压固井项目\项目相关数据资料（未重新命名）\甲方完成数据\HT1-004井油层尾管上交资料1\HT1-004井油层尾管上交资料（甲方）1\XD80502队HT1-004井139.7+168.28㎜上扣扭矩.pdf"
doc=fitz.open(p)
full="\n".join(doc[i].get_text() for i in range(doc.page_count))
# find 168.28 section records
import re
idx=[m.start() for m in re.finditer("168.28", full)]
print("168.28 occurrences:", len(idx))
print(full[:0])
# print a 168.28 record snippet
for i in idx[:3]:
    print("---snippet---")
    print(full[max(0,i-300):i+300])
# summary lines
for kw in ["合格总量","采集总量","重新上扣"]:
    for m in re.finditer(kw, full):
        print(full[max(0,m.start()-60):m.start()+80])
        break
# count TP-SFJ (168.3 thread)
print("TP-SFJ count:", full.count("TP-SFJ"))
print("BG-FJU count:", full.count("BG-FJU"))
print("TP-G2 count:", full.count("TP-G2"))
# torque range for 168.3
for m in re.finditer("168.28", full):
    seg=full[m.start():m.start()+400]
    mm=re.search(r"合格\n(\d+)\n", seg)
    if mm: print("168.28 torque sample:", mm.group(1))
    break
