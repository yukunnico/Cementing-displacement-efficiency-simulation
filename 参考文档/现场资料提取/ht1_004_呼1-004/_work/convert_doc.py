# -*- coding: utf-8 -*-
"""用 Word COM 将 .doc 转 .txt"""
import os, sys, glob
import win32com.client

SRC = r"D:\users\desktop\research\控压固井项目\项目相关数据资料（未重新命名）\甲方完成数据\HT1-004井油层尾管上交资料1\HT1-004井油层尾管上交资料（甲方）1"
OUT = r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\ht1_004_呼1-004\_work"

files = [
    "HT1-004井168.3+139.7mm尾管固井施工记录表.doc",
    "HT1-004井168.3+139.7mm油层尾管固井总结.doc",
    "HT1-004井油层尾管固井作业史.doc",
    "HT1-004井168.3+139.7mm油层尾管控压固井施工设计 (已审批) .doc",
]

wdFormatUnicodeText = 7
app = win32com.client.Dispatch("Word.Application")
app.Visible = False
app.DisplayAlerts = 0
try:
    for f in files:
        src = os.path.join(SRC, f)
        name = f.replace(".doc", "").replace(" ", "_")
        dst = os.path.join(OUT, name + "_word.txt")
        doc = app.Documents.Open(src, ReadOnly=True)
        doc.SaveAs(dst, FileFormat=wdFormatUnicodeText)
        doc.Close(False)
        print("OK:", name, os.path.getsize(dst))
finally:
    app.Quit()
