# -*- coding: utf-8 -*-
"""临时脚本：批量转 doc 到文本以定位呼探1四开（20441-20444）文档内容（只读）。"""
import win32com.client
import os

OUT_DIR = r"D:\users\desktop\research\控压固井项目\cement model\.tmp_research"

TARGETS = [
    ("20441", r"D:\users\desktop\research\控压固井项目\0708\2\204\2044\20441.doc"),
    ("20442", r"D:\users\desktop\research\控压固井项目\0708\2\204\2044\20442.doc"),
    ("20443", r"D:\users\desktop\research\控压固井项目\0708\2\204\2044\20443.doc"),
    ("20444", r"D:\users\desktop\research\控压固井项目\0708\2\204\2044\20444.doc"),
]

word = win32com.client.Dispatch("Word.Application")
word.Visible = False
try:
    for name, path in TARGETS:
        doc = word.Documents.Open(path, ReadOnly=True)
        text = doc.Content.Text
        ntab = len(doc.Tables)
        doc.Close(False)
        out = os.path.join(OUT_DIR, name + ".txt")
        with open(out, "w", encoding="utf-8") as f:
            f.write(text)
        head = text[:80].replace("\r", " ")
        print(name, "len:", len(text), "tables:", ntab, "|", head)
finally:
    word.Quit()
