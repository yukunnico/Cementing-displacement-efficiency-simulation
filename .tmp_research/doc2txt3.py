# -*- coding: utf-8 -*-
import sys, os, io
sys.stdout.reconfigure(encoding='utf-8')
JOBS = [
 (r"D:\users\desktop\research\控压固井项目\0708\甲方完成数据\HT1-003井四开尾管上交资料1\HT1-003井四开尾管上交资料（甲方）1\HT1-003井168.3+139.7mm油层尾管固井总结.doc", "HT1-003固井总结"),
 (r"D:\users\desktop\research\控压固井项目\0708\甲方完成数据\HT1-003井四开尾管上交资料1\HT1-003井四开尾管上交资料（甲方）1\HT1-003井油层尾管固井作业史.doc", "HT1-003作业史"),
 (r"D:\users\desktop\research\控压固井项目\0708\甲方完成数据\HT1-004井油层尾管上交资料1\HT1-004井油层尾管上交资料（甲方）1\HT1-004井168.3+139.7mm油层尾管固井总结.doc", "HT1-004固井总结"),
 (r"D:\users\desktop\research\控压固井项目\0708\甲方完成数据\HT1-004井油层尾管上交资料1\HT1-004井油层尾管上交资料（甲方）1\HT1-004井油层尾管固井作业史.doc", "HT1-004作业史"),
 (r"D:\users\desktop\research\控压固井项目\论文构思及草稿\0708\1\1006\10068.doc", "呼探1_10068"),
 (r"D:\users\desktop\research\控压固井项目\论文构思及草稿\0708\2\204\2041\20411\204111.doc", "呼探1_204111设计"),
]
OUT = r"D:\users\desktop\research\控压固井项目\cement model\.tmp_research\0708_txt"
import win32com.client, pythoncom
pythoncom.CoInitialize()
word = win32com.client.DispatchEx("Word.Application")
word.Visible = False; word.DisplayAlerts = 0
for src, name in JOBS:
    dst = os.path.join(OUT, name + ".txt")
    if os.path.exists(dst):
        print("skip", name); continue
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
