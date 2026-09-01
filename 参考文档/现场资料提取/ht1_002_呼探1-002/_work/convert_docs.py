# -*- coding: utf-8 -*-
"""用 Word COM 直接读 doc.Content.Text（绕开 SaveAs 编码坑），批量导出 UTF-8 txt。"""
import sys
from pathlib import Path

import win32com.client

JOBS = [
    (r"D:\users\desktop\research\控压固井项目\0708\2\HT1-002井现场资料\HT1-002五开139.7mm尾管上交", [
        "HT1-002井139.7mm尾管固井施工小结（写实）.doc",
        "HT1-002井139.7mm尾管完井固井施工设计（研究中心审核版审批版7.1）(1).doc",
        "HT1-002井139.7mm尾管油气井固井作业史.doc",
        "HT1-002井139.7mm油层尾管技术总结.doc",
        "HT1-002井完井尾管固井施工记录表.doc",
    ]),
    (r"D:\users\desktop\research\控压固井项目\0708\2\HT1-002井现场资料\HT1-002井193.7+168.3mm回接上交", [
        "HT1-002井193.7+168.3mm回接固井作业史（终）.doc",
        "HT1-002井回接固井施工设计7.9审批改.doc",
        "HT1-002井油套回接固井施工记录表.doc",
        "HT1-002井油套回接技术总结.doc",
        "探井固井工程量清单（HT1-002 尾管回接）.doc",
    ]),
    (r"D:\users\desktop\research\控压固井项目\0708\2\206\2061\20611\206118", [
        "呼探1-002井219.1mm技术尾管固井施工设计.doc",
        "呼探1-002井219.1mm技术尾管总结.doc",
        "呼探1-002井固井施工记录表.doc",
        "呼探1-002井219.1mm技术尾管固井作业史.doc",
        "固井施工技术要求(供参考).doc",
        "探井固井工程量清单（呼探1-002）.doc",
    ]),
    (r"D:\users\desktop\research\控压固井项目\0708\2\206\2061\20611", [
        "206116.doc", "206117.doc", "206119.doc", "2061191.doc",
        "2061192.doc", "2061194.doc", "2061197.doc",
    ]),
    (r"D:\users\desktop\research\控压固井项目\0708\2\206\2062", [
        "20622.doc", "20623.doc", "20624.doc", "20625.doc",
    ]),
    (r"D:\users\desktop\research\控压固井项目\0708\2\206\2063", [
        "20636.doc", "20637.doc", "20638.doc", "20639.doc",
    ]),
]

OUT_DIR = Path(r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\ht1_002_呼探1-002\_work\doc_txt")
PREFIX = {
    r"D:\users\desktop\research\控压固井项目\0708\2\HT1-002井现场资料\HT1-002五开139.7mm尾管上交": "wukai_",
    r"D:\users\desktop\research\控压固井项目\0708\2\HT1-002井现场资料\HT1-002井193.7+168.3mm回接上交": "huijie_",
    r"D:\users\desktop\research\控压固井项目\0708\2\206\2061\20611\206118": "t219_",
    r"D:\users\desktop\research\控压固井项目\0708\2\206\2061\20611": "L2061x_",
    r"D:\users\desktop\research\控压固井项目\0708\2\206\2062": "L2062x_",
    r"D:\users\desktop\research\控压固井项目\0708\2\206\2063": "L2063x_",
}


def main() -> None:
    word = win32com.client.DispatchEx("Word.Application")
    word.Visible = False
    word.DisplayAlerts = 0
    try:
        for root, names in JOBS:
            for name in names:
                src = Path(root) / name
                if not src.exists():
                    print(f"[MISS] {src}")
                    continue
                doc = word.Documents.Open(str(src), ReadOnly=True, AddToRecentFiles=False)
                # 统计内嵌图形（判断是否扫描件）
                n_inline = doc.InlineShapes.Count
                n_shapes = doc.Shapes.Count
                text = doc.Content.Text
                doc.Close(False)
                # 单元格结束符 \x07 换成制表符，段落 \r 换行
                text = text.replace("\x07", "\t").replace("\r", "\n")
                prefix = PREFIX.get(root, "")
                out = OUT_DIR / (prefix + name.rsplit(".", 1)[0] + ".txt")
                out.write_text(text, encoding="utf-8")
                print(f"[OK] {out.name}  {len(text)} chars  img={n_inline}+{n_shapes}")
    finally:
        word.Quit()


if __name__ == "__main__":
    sys.exit(main())
