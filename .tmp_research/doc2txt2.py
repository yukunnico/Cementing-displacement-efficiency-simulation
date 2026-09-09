# -*- coding: utf-8 -*-
import sys, os, io
sys.stdout.reconfigure(encoding='utf-8')
SRC = r"D:\users\desktop\research\控压固井项目\论文构思及草稿\0708"
OUT = r"D:\users\desktop\research\控压固井项目\cement model\.tmp_research\0708_txt"
FILES = [
    r"2\202\2021\20211.doc",
    r"2\202\2021\20212.doc",
    r"2\202\2021\20213.doc",
    r"2\202\2021\20216.doc",
    r"2\205\2051\20519.xlsx",
    r"2\205\2051\20517\205171.docx",
    r"1\1003\10033.doc",
    r"1\1004\10041\100413.PDF",
]
import win32com.client, pythoncom
pythoncom.CoInitialize()
word = win32com.client.DispatchEx("Word.Application")
word.Visible = False; word.DisplayAlerts = 0
for rel in FILES:
    src = os.path.join(SRC, rel)
    if not os.path.exists(src):
        print("MISSING", rel); continue
    name = os.path.splitext(os.path.basename(rel))[0]
    dst = os.path.join(OUT, name + ".txt")
    try:
        doc = word.Documents.Open(src, ReadOnly=True, AddToRecentFiles=False)
        doc.SaveAs2(dst, FileFormat=7)
        doc.Close(False)
        raw = open(dst,'rb').read()
        try: t = raw.decode('gb18030')
        except: t = raw.decode('utf-16', errors='replace')
        io.open(dst,'w',encoding='utf-8').write(t)
        print("OK", name, os.path.getsize(dst))
    except Exception as e:
        print("FAIL", name, repr(e)[:120])
word.Quit()
