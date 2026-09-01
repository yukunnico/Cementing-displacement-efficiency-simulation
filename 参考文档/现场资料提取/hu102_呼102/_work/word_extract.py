# -*- coding: utf-8 -*-
"""独立重新解析 0708 原件：Word COM 提取 .doc/.docx 全文；不依赖提取包产物。"""
import sys, os
from pathlib import Path
import win32com.client

WORK = Path(__file__).parent
FILES = {
    "src_10043.doc": "out_10043.txt",
    "src_20211.doc": "out_20211.txt",
    "src_20212.doc": "out_20212.txt",
    "src_20214.doc": "out_20214.txt",
    "src_20216.doc": "out_20216.txt",
    "src_20234.doc": "out_20234.txt",
    "src_202391.docx": "out_202391.txt",
    "src_20251.doc": "out_20251.txt",
    "src_20223.doc": "out_20223.txt",
}

def main():
    word = win32com.client.DispatchEx("Word.Application")
    word.Visible = False
    word.DisplayAlerts = 0
    try:
        for src, dst in FILES.items():
            src_p = WORK / src
            dst_p = WORK / dst
            if not src_p.exists():
                print(f"[MISS] {src}")
                continue
            try:
                doc = word.Documents.Open(str(src_p), ReadOnly=True, AddToRecentFiles=False)
                text = doc.Content.Text
                # 表格文本也包含在 Content 中；再补 tables 展开（防丢失）
                if doc.Tables.Count > 0:
                    tb = []
                    for i in range(1, doc.Tables.Count + 1):
                        t = doc.Tables(i)
                        try:
                            tb.append(f"\n===== TABLE {i} =====\n" + t.Range.Text)
                        except Exception as e:
                            tb.append(f"\n===== TABLE {i} (err {e}) =====\n")
                    text = text + "\n".join(tb)
                doc.Close(False)
                dst_p.write_text(text, encoding="utf-8", errors="replace")
                print(f"[OK] {src} -> {dst} ({len(text)} chars)")
            except Exception as e:
                print(f"[ERR] {src}: {e}")
    finally:
        word.Quit()

if __name__ == "__main__":
    main()
