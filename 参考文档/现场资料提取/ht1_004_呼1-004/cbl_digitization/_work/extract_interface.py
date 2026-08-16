# -*- coding: utf-8 -*-
"""HT1-004 CBL 矢量数字化核心：一界面/二界面解释填充 → 深度区间 + 质量分级 + 合格率。"""
import fitz, os, re, json
from collections import defaultdict

SRC_DIR = r"D:\users\desktop\research\控压固井项目\项目相关数据资料（未重新命名）\HT1-003、HT1-004井固井声幅图"
OUT = r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\ht1_004_呼1-004\cbl_digitization\_work"

CYAN = (0.0, 0.6, 0.6)
GREEN = (0.0, 0.5, 0.0)
WHITE = (1.0, 1.0, 1.0)

def close(c, ref, tol=0.04):
    return all(abs(a-b) <= tol for a, b in zip(c, ref))

def get_calib(page):
    """深度刻度 → (k, b) depth = k*y + b，取左侧 x<70 刻度"""
    td = page.get_text('dict')
    pts = []
    for blk in td['blocks']:
        if blk['type'] != 0: continue
        for line in blk['lines']:
            for span in line['spans']:
                t = span['text'].strip()
                if re.fullmatch(r'\d{4}', t) and span['bbox'][0] < 70:
                    pts.append((int(t), (span['bbox'][1]+span['bbox'][3])/2))
    pts.sort()
    ys = [p[0] for p in pts]; xs = [p[1] for p in pts]
    k = (ys[-1]-ys[0])/(xs[-1]-xs[0]); b = ys[0] - k*xs[0]
    return k, b

def y_to_depth(y, k, b): return k*y + b

def extract_interface(page, track_x0, track_x1, k, b, y_data_start=915.0):
    """对一道，逐 y 采样分类；返回 [(depth_top, depth_bot, class)] 合并区间。"""
    fills = []
    for it in page.get_drawings():
        if it.get('type') == 'f':
            c = it.get('fill'); r = it['rect']
            if not c: continue
            # 道内：中心 x 在范围内
            cx = (r.x0+r.x1)/2
            if track_x0-2 <= cx <= track_x1+2 and r.y1 > y_data_start:
                if close(c, CYAN):
                    cls = 'cyan'
                elif close(c, GREEN):
                    cls = 'green'
                else:
                    continue
                fills.append((cls, r.y0, r.y1, r.width))
    # 逐 y 采样（步长 0.5pt≈0.035m）
    ymin = min((f[1] for f in fills), default=None)
    ymax = max((f[2] for f in fills), default=None)
    if ymin is None:
        return []  # 全空
    # 生成采样分类
    samples = []  # (y, cls)
    y = ymin
    step = 0.5
    while y <= ymax:
        # 该 y 处覆盖的填充
        covering = [f for f in fills if f[1] <= y <= f[2]]
        cls = 'empty'
        for c, y0_, y1_, w in covering:
            if c == 'green':
                cls = 'green'; break
            elif c == 'cyan' and w > 15:
                cls = 'cyan_solid'  # 良好
            elif c == 'cyan' and cls in ('empty','cyan_strip'):
                cls = 'cyan_strip'
        samples.append((y, cls))
        y += step
    # 合并连续同类
    segs = []
    cur_cls = None; cur_y0 = None
    for y, cls in samples:
        if cls == cur_cls:
            continue
        if cur_cls is not None:
            segs.append((cur_y0, y, cur_cls))
        cur_cls = cls; cur_y0 = y
    if cur_cls is not None:
        segs.append((cur_y0, samples[-1][0]+step, cur_cls))
    # 深度
    out = [(y_to_depth(a,k,b), y_to_depth(b2,k,b), m) for a,b2,m in segs]
    return out

def classify_display(cls):
    return {'green':'胶结差','cyan_strip':'胶结中等','cyan_solid':'胶结良好','empty':'空套管'}[cls]

def main():
    JOBS = [
        ("HT1-004_固井质量测井评价图_7382-7581_20260721_完井.pdf", "tail", 7382, 7581),
        ("HT1-004_固井质量测井评价图_5249-7384_20260721_完井.pdf", "mid", 5249, 7384),
        ("HT1-004_固井质量测井评价图_11-5251_20260721_完井.pdf", "shallow", 11, 5251),
    ]
    all_rows = []
    for fn, tag, dmin, dmax in JOBS:
        doc = fitz.open(os.path.join(SRC_DIR, fn))
        page = doc[0]
        k, b = get_calib(page)
        print(f"=== {tag} ({dmin}-{dmax}m) calib k={k:.5f} b={b:.2f}")
        # 一界面 688.5-719.7, 二界面 719.7-750.9
        for name, x0, x1 in [("一界面", 688.5, 719.7), ("二界面", 719.7, 750.9)]:
            segs = extract_interface(page, x0, x1, k, b)
            # 裁剪到测量井段
            segs = [(a,b2,m) for a,b2,m in segs if b2 > dmin and a < dmax]
            # 合并类名
            print(f"  {name}:")
            tot = 0; good = 0
            for a, b2, m in segs:
                md = max(a, dmin); mdx = min(b2, dmax)
                if mdx <= md: continue
                length = mdx - md
                tot += length
                if m in ('cyan_strip','cyan_solid'):
                    good += length
                print(f"    {md:.1f}-{mdx:.1f}m  {classify_display(m)}  len={length:.1f}m")
                all_rows.append((tag, name, round(md,2), round(mdx,2), m, round(length,2)))
            if tot > 0:
                print(f"    合计={tot:.1f}m  中等及以上占比={good/tot*100:.1f}%")
        doc.close()
    # 存 JSON
    with open(os.path.join(OUT, "interface_segments.json"), "w", encoding="utf-8") as f:
        json.dump(all_rows, f, ensure_ascii=False, indent=1)
    print("saved interface_segments.json")

if __name__ == "__main__":
    main()
