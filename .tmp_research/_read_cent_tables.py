# -*- coding: utf-8 -*-
"""临时脚本：读 HT1-002/HT1-001 设计 doc 中扶正器表格结构（只读）。"""
import win32com.client

TARGETS = [
    ("HT1-002设计", r"D:\users\desktop\research\控压固井项目\0708\2\HT1-002井现场资料\HT1-002五开139.7mm尾管上交\HT1-002井139.7mm尾管完井固井施工设计（研究中心审核版审批版7.1）(1).doc"),
    ("HT1-001设计", r"D:\users\desktop\research\控压固井项目\0708\2\HT1-001完井交\HT1-001井139.7+168.3mm尾管完井固井施工设计（审核）.doc"),
]

word = win32com.client.Dispatch("Word.Application")
word.Visible = False
try:
    for name, path in TARGETS:
        doc = word.Documents.Open(path, ReadOnly=True)
        print(f"########## {name}: {len(doc.Tables)} tables")
        for ti in range(1, len(doc.Tables) + 1):
            tb = doc.Tables(ti)
            rows, cols = tb.Rows.Count, tb.Columns.Count
            if rows < 3 or cols > 12:
                continue
            # check if table contains 扶正 keyword
            body_probe = ""
            for r in range(1, min(rows, 6) + 1):
                for c in range(1, min(cols, 10) + 1):
                    try:
                        body_probe += tb.Cell(r, c).Range.Text
                    except Exception:
                        pass
            if "扶正" in body_probe or "居中" in body_probe:
                print(f"--- table {ti}: {rows}x{cols} ---")
                for r in range(1, rows + 1):
                    cells = []
                    for c in range(1, cols + 1):
                        try:
                            t = tb.Cell(r, c).Range.Text.replace("\r\x07", "").replace("\r", "")
                            cells.append(t.strip())
                        except Exception:
                            cells.append("<merged>")
                    print(" | ".join(cells))
        doc.Close(False)
finally:
    word.Quit()
