import sys, fitz
sys.stdout.reconfigure(encoding='utf-8')
d = r"D:\users\desktop\research\控压固井项目\0708\HT1-003、HT1-004井固井声幅图"
import glob, os
for f in sorted(glob.glob(d + r"\HT1-004_*.pdf")):
    doc = fitz.open(f)
    print("="*100); print("FILE:", os.path.basename(f), "pages:", len(doc))
    t = doc[0].get_text()
    print(t[:2400])
    doc.close()
