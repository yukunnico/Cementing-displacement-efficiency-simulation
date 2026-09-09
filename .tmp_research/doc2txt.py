# -*- coding: utf-8 -*-
"""0708 老 .doc 批量转 txt：优先 Word COM，失败则 OLE 粗提取。"""
import sys, os, io
sys.stdout.reconfigure(encoding='utf-8')
SRC = r"D:\users\desktop\research\控压固井项目\论文构思及草稿\0708"
OUT = r"D:\users\desktop\research\控压固井项目\cement model\.tmp_research\0708_txt"
FILES = [
    r"2\202\2023\20233.doc",
    r"2\203\2031\20318.doc",
    r"2\204\2041\20415\204151.doc",
    r"2\HT1-001完井交\HT1-001井139.7mm尾管固井施工小结（写实）.doc",
    r"HT1-002井现场资料\HT1-002五开139.7mm尾管上交\HT1-002井139.7mm尾管固井施工小结（写实）.doc",
    r"2\205\2051\20513\205131.doc",
    r"2\205\2051\205191.doc",
    r"2\205\2051\205192.doc",
    r"2\206\2061\20611\206116.doc",
    r"2\206\2061\20611\206117.doc",
    r"2\206\2061\20611\206119.doc",
    r"2\206\2061\20611\2061191.doc",
    r"2\206\2061\20611\2061192.doc",
]
os.makedirs(OUT, exist_ok=True)
use_com = True
word = None
try:
    import win32com.client, pythoncom
    pythoncom.CoInitialize()
    word = win32com.client.DispatchEx("Word.Application")
    word.Visible = False
    word.DisplayAlerts = 0
except Exception as e:
    print("COM init failed:", e)
    use_com = False

def ole_extract(path):
    """OLE WordDocument 流粗提取：UTF-16LE 与 GBK 连续可读段。"""
    import olefile, re
    ole = olefile.OleFileIO(path)
    data = ole.openstream("WordDocument").read()
    ole.close()
    chunks = []
    # UTF-16LE 扫描：连续中日韩/ASCII 可读
    txt = data.decode('utf-16-le', errors='ignore')
    for m in re.finditer(r'[\u4e00-\u9fffA-Za-z0-9，。：:；;（）()\.\-—~、\s/%×Φφ‰→↑↓]{12,}', txt):
        s = m.group(0)
        if sum('\u4e00' <= c <= '\u9fff' for c in s) >= 2:
            chunks.append(s)
    return '\n'.join(chunks)

for rel in FILES:
    src = os.path.join(SRC, rel)
    name = os.path.splitext(os.path.basename(rel))[0]
    dst = os.path.join(OUT, name + ".txt")
    if os.path.exists(dst):
        print("skip", name); continue
    ok = False
    if use_com:
        try:
            doc = word.Documents.Open(src, ReadOnly=True, AddToRecentFiles=False)
            doc.SaveAs2(dst, FileFormat=7, Encoding=65001)  # 7=wdFormatEncodedText, 65001=UTF-8
            doc.Close(False)
            ok = os.path.exists(dst) and os.path.getsize(dst) > 100
        except Exception as e:
            print("COM fail", name, repr(e)[:120])
    if not ok:
        try:
            txt = ole_extract(src)
            with io.open(dst, 'w', encoding='utf-8') as f:
                f.write(txt)
            print("OLE", name, os.path.getsize(dst))
        except Exception as e:
            print("OLE fail", name, repr(e)[:120])
    else:
        print("COM", name, os.path.getsize(dst))
if word is not None:
    word.Quit()
