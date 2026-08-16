# -*- coding: utf-8 -*-
"""HT1-003 CBL 补充统计：二界面分窗口 + 连续25m中等及以上检查 + 关键窗口"""
import json, os

WORK = r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\ht1_003_呼1-003\_work"
with open(os.path.join(WORK, "interface_segments.json"), encoding="utf-8") as f:
    segs = json.load(f)

one = {}; two = {}
for tag, name, md, mdx, cls, ln in segs:
    for m in range(int(round(md)), int(round(mdx))):
        if name == "一界面": one[m] = cls
        else: two[m] = cls

def cn(c): return {'green':'差','cyan_strip':'中等','cyan_solid':'良好','empty':'空套管'}[c]

def stat(mp, rng, label):
    cnt = {'良好':0,'中等':0,'差':0,'空套管':0}
    for m in rng:
        cnt[cn(mp.get(m,'empty'))] += 1
    tot = sum(cnt.values()); mid = cnt['中等']+cnt['良好']
    if tot:
        print(f"  [{label}] 良好={cnt['良好']} 中等={cnt['中等']} 差={cnt['差']} 空={cnt['空套管']} | 评价={tot}m 中等及以上={mid}m ({mid/tot*100:.1f}%)")
    return mid/tot*100 if tot else None

print("===== 二界面分窗口 =====")
stat(two, range(5568, 7596), "二界面 5568-7595(单层套管评价段)")
stat(two, range(5308, 7515), "二界面 尾管评价段 5307.54-7514.21")
stat(two, range(7442, 7596), "二界面 目的井段 7442-7595")
stat(two, range(7353, 7561), "二界面 油气水层段 7353-7560")

print("\n===== 连续25m中等及以上（一界面）检查 =====")
# 找出所有 连续良好/中等 段（run），并统计长度>=25m的
runs = []; cur = None; start = None
for m in range(15, 7596):
    c = one.get(m, 'empty')
    if c in ('良好','中等'):
        if cur is None: start = m; cur = c
    else:
        if cur is not None:
            runs.append((start, m-1))
            cur = None
if cur is not None: runs.append((start, 7595))
long_runs = [(a,b,b-a+1) for a,b in runs if b-a+1 >= 25]
print(f"  连续良好/中等段总数={len(runs)}, 其中 >=25m 的={len(long_runs)}")
print("  最长10段:")
for a,b,l in sorted(long_runs, key=lambda x:-x[2])[:10]:
    print(f"    {a}-{b}m  len={l}m")

print("\n===== 官方红线1关键窗口 5543-5568m（上层套管鞋5568以上25m, 假设）一界面 =====")
stat(one, range(5543, 5569), "5543-5568 一界面")
# 检查该窗口内最长连续中等及以上
r = []; cur = None; s = None
for m in range(5543, 5569):
    c = one.get(m, 'empty')
    if c in ('良好','中等'):
        if cur is None: s = m; cur = c
    else:
        if cur is not None: r.append((s, m-1)); cur = None
if cur is not None: r.append((s, 5568))
print(f"  该窗口内连续中等及以上段: {[(a,b,b-a+1) for a,b in r]}")

print("\n===== 红线1另一候选 5286-5316m（悬挂器/回接套管鞋附近 25m）一界面 =====")
stat(one, range(5286, 5317), "5286-5316 一界面")
r2 = []; cur = None; s = None
for m in range(5286, 5317):
    c = one.get(m, 'empty')
    if c in ('良好','中等'):
        if cur is None: s = m; cur = c
    else:
        if cur is not None: r2.append((s, m-1)); cur = None
if cur is not None: r2.append((s, 5316))
print(f"  该窗口内连续中等及以上段: {[(a,b,b-a+1) for a,b in r2]}")

# 油气水层段长度
oil_gas = [(7352.9,7357.7),(7380.2,7383.2),(7398.0,7407.4),(7483.2,7489.8),
           (7499.4,7502.8),(7504.0,7508.6),(7526.6,7529.2),(7548.4,7559.6)]
print("\n===== 油气水层段长度（均<25m?） =====")
for a,b in oil_gas:
    print(f"  {a:.1f}-{b:.1f}m  len={b-a:.1f}m")
