# -*- coding: utf-8 -*-
"""临时脚本：批量读各井设计 doc 的扶正器/居中度表格结构（只读）。"""
import win32com.client

TARGETS = [
    ("HT1-003设计", r"D:\users\desktop\research\控压固井项目\0708\甲方完成数据\HT1-003井四开尾管上交资料1\HT1-003井四开尾管上交资料（甲方）1\HT1-003井168.3+139.7mm油层尾管控压固井施工设计 (已审批) .doc"),
    ("HT1-004设计", r"D:\users\desktop\research\控压固井项目\0708\甲方完成数据\HT1-004井油层尾管上交资料1\HT1-004井油层尾管上交资料（甲方）1\HT1-004井168.3+139.7mm油层尾管控压固井施工设计 (已审批) .doc"),
    ("呼103设计20313", r"D:\users\desktop\research\控压固井项目\0708\2\203\2031\20313.doc"),
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
            if rows < 2 or cols > 12:
                continue
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
                            cells.append("<m>")
                    print(" | ".join(cells))
        doc.Close(False)
finally:
    word.Quit()
