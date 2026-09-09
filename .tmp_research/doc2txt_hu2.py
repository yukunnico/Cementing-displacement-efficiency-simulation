# -*- coding: utf-8 -*-
import sys, os, io, glob
sys.stdout.reconfigure(encoding='utf-8')
SRC = r"D:\users\desktop\research\控压固井项目\论文构思及草稿\0708\2\HT1-002井现场资料"
OUT = r"D:\users\desktop\research\控压固井项目\cement model\.tmp_research\0708_txt"
import win32com.client, pythoncom
pythoncom.CoInitialize()
word = win32com.client.DispatchEx("Word.Application")
word.Visible = False; word.DisplayAlerts = 0
targets = [p for p in glob.glob(os.path.join(SRC, "**", "*.doc"), recursive=True) if ("139.7mm尾管" in p or "作业史" in p)]
for src in targets:
    name = os.path.splitext(os.path.basename(src))[0]
    dst = os.path.join(OUT, "HU2_" + name + ".txt")
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
        print("FAIL", name, repr(e)[:150])
word.Quit()
