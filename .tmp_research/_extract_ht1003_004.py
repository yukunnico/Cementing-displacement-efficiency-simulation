# -*- coding: utf-8 -*-
"""临时脚本：提取 HT1-003/004 固井设计 doc 全文（只读）。"""
import win32com.client
import os

OUT_DIR = r"D:\users\desktop\research\控压固井项目\cement model\.tmp_research"

TARGETS = [
    ("HT1-003设计", r"D:\users\desktop\research\控压固井项目\0708\甲方完成数据\HT1-003井四开尾管上交资料1\HT1-003井四开尾管上交资料（甲方）1\HT1-003井168.3+139.7mm油层尾管控压固井施工设计 (已审批) .doc"),
    ("HT1-004设计", r"D:\users\desktop\research\控压固井项目\0708\甲方完成数据\HT1-004井油层尾管上交资料1\HT1-004井油层尾管上交资料（甲方）1\HT1-004井168.3+139.7mm油层尾管控压固井施工设计 (已审批) .doc"),
]

word = win32com.client.Dispatch("Word.Application")
word.Visible = False
try:
    for name, path in TARGETS:
        doc = word.Documents.Open(path, ReadOnly=True)
        text = doc.Content.Text
        doc.Close(False)
        out = os.path.join(OUT_DIR, name + ".txt")
        with open(out, "w", encoding="utf-8") as f:
            f.write(text)
        print(name, "done, len:", len(text))
finally:
    word.Quit()
