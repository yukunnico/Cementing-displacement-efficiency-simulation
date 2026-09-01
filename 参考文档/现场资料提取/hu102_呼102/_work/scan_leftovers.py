# -*- coding: utf-8 -*-
"""遗留未识别文件快速排查：Word COM 批量提取 + PDF 文字层 + xls/xlsx sheet 概览。
只输出每个文件的前若干字符与关键词命中，判定是否含模型相关数据。"""
import sys, re
from pathlib import Path

WORK = Path(r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\hu102_呼102\_work")
OUT = WORK / "leftovers"
OUT.mkdir(exist_ok=True)

DOCS = {
    "10044": r"D:\users\desktop\research\控压固井项目\0708\1\1004\10044.doc",
    "10046": r"D:\users\desktop\research\控压固井项目\0708\1\1004\10046.doc",
    "10047": r"D:\users\desktop\research\控压固井项目\0708\1\1004\10047.doc",
    "10048": r"D:\users\desktop\research\控压固井项目\0708\1\1004\10048.doc",
    "10049": r"D:\users\desktop\research\控压固井项目\0708\1\1004\10049.doc",
    "100491": r"D:\users\desktop\research\控压固井项目\0708\1\1004\100491.doc",
    "20213": r"D:\users\desktop\research\控压固井项目\0708\2\202\2021\20213.doc",
    "20231": r"D:\users\desktop\research\控压固井项目\0708\2\202\2023\20231.doc",
    "20233": r"D:\users\desktop\research\控压固井项目\0708\2\202\2023\20233.doc",
    "20236": r"D:\users\desktop\research\控压固井项目\0708\2\202\2023\20236.doc",
    "20237": r"D:\users\desktop\research\控压固井项目\0708\2\202\2023\20237.doc",
    "20239": r"D:\users\desktop\research\控压固井项目\0708\2\202\2023\20239.doc",
    "202392": r"D:\users\desktop\research\控压固井项目\0708\2\202\2023\202392.doc",
    "202399": r"D:\users\desktop\research\控压固井项目\0708\2\202\2023\202399.doc",
    "20247": r"D:\users\desktop\research\控压固井项目\0708\2\202\2024\20247.doc",
    "20248": r"D:\users\desktop\research\控压固井项目\0708\2\202\2024\20248.doc",
    "20249": r"D:\users\desktop\research\控压固井项目\0708\2\202\2024\20249.doc",
    "202492": r"D:\users\desktop\research\控压固井项目\0708\2\202\2024\202492.doc",
    "202493": r"D:\users\desktop\research\控压固井项目\0708\2\202\2024\202493.doc",
    "202496": r"D:\users\desktop\research\控压固井项目\0708\2\202\2024\202496.doc",
    "202497": r"D:\users\desktop\research\控压固井项目\0708\2\202\2024\202497.doc",
    "20232": r"D:\users\desktop\research\控压固井项目\0708\2\202\2023\20232.docx",
    "202394": r"D:\users\desktop\research\控压固井项目\0708\2\202\2023\202394.docx",
    "20261": r"D:\users\desktop\research\控压固井项目\0708\2\202\2026\20261.docx",
    "20266": r"D:\users\desktop\research\控压固井项目\0708\2\202\2026\20266.docx",
    "20267": r"D:\users\desktop\research\控压固井项目\0708\2\202\2026\20267.doc",
    "20268": r"D:\users\desktop\research\控压固井项目\0708\2\202\2026\20268.doc",
    "202691": r"D:\users\desktop\research\控压固井项目\0708\2\202\2026\202691.docx",
    "202692": r"D:\users\desktop\research\控压固井项目\0708\2\202\2026\202692.doc",
    "202693": r"D:\users\desktop\research\控压固井项目\0708\2\202\2026\202693.doc",
    "202696": r"D:\users\desktop\research\控压固井项目\0708\2\202\2026\202696.docx",
    "202697": r"D:\users\desktop\research\控压固井项目\0708\2\202\2026\202697.doc",
}

PDFS = {
    "100419": r"D:\users\desktop\research\控压固井项目\0708\1\1004\10041\100419.PDF",
    "100420": r"D:\users\desktop\research\控压固井项目\0708\1\1004\10041\100420.PDF",
    "100421": r"D:\users\desktop\research\控压固井项目\0708\1\1004\10041\100421.PDF",
    "100422": r"D:\users\desktop\research\控压固井项目\0708\1\1004\10041\100422.PDF",
    "100423": r"D:\users\desktop\research\控压固井项目\0708\1\1004\10041\100423.PDF",
    "100424": r"D:\users\desktop\research\控压固井项目\0708\1\1004\10041\100424.PDF",
    "20225": r"D:\users\desktop\research\控压固井项目\0708\2\202\2022\20225.pdf",
    "20238": r"D:\users\desktop\research\控压固井项目\0708\2\202\2023\20238.pdf",
    "202491": r"D:\users\desktop\research\控压固井项目\0708\2\202\2024\202491.pdf",
    "202495": r"D:\users\desktop\research\控压固井项目\0708\2\202\2024\202495.pdf",
}

XLS = {
    "20241.xls": r"D:\users\desktop\research\控压固井项目\0708\2\202\2024\20241.xls",
    "202494.xlsx": r"D:\users\desktop\research\控压固井项目\0708\2\202\2024\202494.xlsx",
    "202397.xls": r"D:\users\desktop\research\控压固井项目\0708\2\202\2023\202397.xls",
    "202398.xlsx": r"D:\users\desktop\research\控压固井项目\0708\2\202\2023\202398.xlsx",
    "202393.xlsm": r"D:\users\desktop\research\控压固井项目\0708\2\202\2023\202393.xlsm",
    "20269.xls": r"D:\users\desktop\research\控压固井项目\0708\2\202\2026\20269.xls",
}

KEYWORDS = ["井身结构", "套管", "尾管", "水泥浆", "隔离液", "泵压", "排量", "替浆", "顶替",
            "密度", "流变", "粘度", "声幅", "CBL", "VDL", "固井质量", "扶正器", "注水泥",
            "悬挂器", "井径", "井斜", "压塞液", "平衡液", "前置液"]

def summarize(tag, text):
    text = text or ""
    hits = {k: text.count(k) for k in KEYWORDS if text.count(k) > 0}
    head = re.sub(r"\s+", " ", text[:220])
    print(f"[{tag}] chars={len(text)}")
    print(f"    head: {head}")
    top = sorted(hits.items(), key=lambda x: -x[1])[:10]
    print(f"    kw: {top}")
    (OUT / f"{tag}.txt").write_text(text, encoding="utf-8", errors="replace")

def main():
    import win32com.client
    word = win32com.client.DispatchEx("Word.Application")
    word.Visible = False
    word.DisplayAlerts = 0
    try:
        for tag, p in DOCS.items():
            try:
                doc = word.Documents.Open(p, ReadOnly=True, AddToRecentFiles=False)
                text = doc.Content.Text
                if doc.Tables.Count > 0:
                    for i in range(1, min(doc.Tables.Count, 6) + 1):
                        try:
                            text += "\n[TABLE %d]\n" % i + doc.Tables(i).Range.Text
                        except Exception:
                            pass
                doc.Close(False)
                summarize(tag, text)
            except Exception as e:
                print(f"[{tag}] WORD ERR: {e}")
    finally:
        word.Quit()

    import pdfplumber
    for tag, p in PDFS.items():
        try:
            with pdfplumber.open(p) as pdf:
                texts = []
                for pg in pdf.pages[:4]:
                    texts.append(pg.extract_text() or "")
                summarize(tag + "_pdf", "\n".join(texts))
        except Exception as e:
            print(f"[{tag}] PDF ERR: {e}")

    import pandas as pd
    for tag, p in XLS.items():
        try:
            xl = pd.ExcelFile(p)
            print(f"[{tag}] sheets={xl.sheet_names[:8]}")
            for sn in xl.sheet_names[:3]:
                df = xl.parse(sn, header=None, nrows=6)
                vals = [str(v)[:18] for v in df.values.flatten() if str(v) != "nan"][:14]
                print(f"    {sn}: {vals}")
        except Exception as e:
            print(f"[{tag}] XLS ERR: {e}")

if __name__ == "__main__":
    main()
