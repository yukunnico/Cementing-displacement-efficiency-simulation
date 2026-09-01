# -*- coding: utf-8 -*-
"""将 _work 下 Word 转出的 txt（UTF-16LE）统一转为 UTF-8，便于后续读取。"""
from pathlib import Path

work = Path(__file__).parent
for p in sorted(work.glob("*.txt")):
    raw = p.read_bytes()
    enc_used = None
    text = None
    for enc in ("utf-16", "gbk", "utf-8-sig", "utf-8"):
        try:
            candidate = raw.decode(enc)
            if "\x00" in candidate:
                continue
            text, enc_used = candidate, enc
            break
        except Exception:
            continue
    if text is None:
        print("FAIL", p.name)
        continue
    p.write_text(text, encoding="utf-8")
    print(f"{p.name}\t{enc_used}\t{len(text)}")
