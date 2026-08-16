# -*- coding: utf-8 -*-
"""修复版连续段检查"""
import json, os
WORK = r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\ht1_003_呼1-003\_work"
with open(os.path.join(WORK, "interface_segments.json"), encoding="utf-8") as f:
    segs = json.load(f)
one = {}
for tag, name, md, mdx, cls, ln in segs:
    if name != "一界面": continue
    for m in range(int(round(md)), int(round(mdx))):
        one[m] = cls  # 'green'/'cyan_strip'/'cyan_solid'

def find_runs(mp, m0, m1):
    runs = []; cur = None; s = None
    for m in range(m0, m1+1):
        c = mp.get(m, 'empty')
        if c in ('cyan_strip','cyan_solid'):
            if cur is None: s = m; cur = c
        else:
            if cur is not None: runs.append((s, m-1, m-s)); cur = None
    if cur is not None: runs.append((s, m1, m1-s+1))
    return runs

print("全井 15-7595 一界面 连续中等及以上段:")
runs = find_runs(one, 15, 7595)
long_runs = [(a,b,l) for a,b,l in runs if l >= 25]
print(f"  总数={len(runs)}, >=25m 的={len(long_runs)}")
for a,b,l in sorted(long_runs, key=lambda x:-x[2])[:12]:
    print(f"    {a}-{b}m  len={l}m")
total_mid = sum(l for _,_,l in runs)
print(f"  中等及以上总长度={total_mid}m")

print("\n关键窗口（官方红线1 候选）:")
for name, a, b in [("5543-5568 技术套管鞋5568以上25m", 5543, 5568), ("5286-5316 回接鞋/悬挂器附近", 5286, 5316)]:
    rr = find_runs(one, a, b)
    print(f"  {name}: 连续段={[(x,y,l) for x,y,l in rr]} 最长={max([l for _,_,l in rr]+[0])}m")

print("\n油气水层段内连续中等及以上（8段）:")
oil_gas = [(7352.9,7357.7),(7380.2,7383.2),(7398.0,7407.4),(7483.2,7489.8),
           (7499.4,7502.8),(7504.0,7508.6),(7526.6,7529.2),(7548.4,7559.6)]
for a,b in oil_gas:
    rr = find_runs(one, int(round(a)), int(round(b)))
    print(f"  {a:.1f}-{b:.1f}m: {[(x,y,l) for x,y,l in rr]} 最长={max([l for _,_,l in rr]+[0])}m")
