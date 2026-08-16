# -*- coding: utf-8 -*-
"""用 Word COM 把 .doc 转成 .txt，用于 antiword 无法处理的旧 Word 文档"""
import os, sys, glob
import win32com.client

SRC = r"D:\users\desktop\research\控压固井项目\项目相关数据资料（未重新命名）\甲方完成数据\HT1-003井四开尾管上交资料1\HT1-003井四开尾管上交资料（甲方）1"
OUT = r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\ht1_003_呼1-003\_work"

files = [
    os.path.join(SRC, "HT1-003井168.3+139.7mm油层尾管控压固井施工设计 (已审批) .doc"),
    os.path.join(SRC, "HT1-003井油层尾管固井作业史.doc"),
]

word = win32com.client.DispatchEx("Word.Application")
word.Visible = False
word.DisplayAlerts = 0
try:
    for f in files:
        try:
            doc = word.Documents.Open(f, ReadOnly=True, AddToRecentFiles=False)
            txt_path = os.path.join(OUT, os.path.splitext(os.path.basename(f))[0] + "_word.txt")
            doc.SaveAs(txt_path, FileFormat=7)  # 7 = wdFormatEncodedText/Unicode text
            doc.Close(False)
            print("OK:", os.path.basename(f))
        except Exception as e:
            print("FAIL:", os.path.basename(f), repr(e))
finally:
    word.Quit()
print("DONE")
