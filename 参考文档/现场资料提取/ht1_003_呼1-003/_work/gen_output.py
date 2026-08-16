# -*- coding: utf-8 -*-
"""生成 HT1-003 井 16 个标准提取文件"""
import csv, os

OUT = r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\ht1_003_呼1-003"
WELL = "ht1_003"
DSGN = "HT1-003井168.3+139.7mm油层尾管控压固井施工设计 (已审批).doc"
REC  = "HT1-003井168.3+139.7mm尾管固井施工记录表.doc"
SUM  = "HT1-003井168.3+139.7mm油层尾管固井总结.doc"
HIST = "HT1-003井油层尾管固井作业史.doc"
LAB  = "化验报告/HT1-003 油层尾管 化验报告.docx"

def w(fname, header, rows):
    p = os.path.join(OUT, fname)
    with open(p, "w", encoding="utf-8-sig", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(header)
        wr.writerows(rows)
    print(f"wrote {fname}: {len(rows)} rows")
