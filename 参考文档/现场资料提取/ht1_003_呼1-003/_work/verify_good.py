# -*- coding: utf-8 -*-
import json, os
WORK = r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\ht1_003_呼1-003\_work"
with open(os.path.join(WORK, "interface_segments.json"), encoding="utf-8") as f:
    segs = json.load(f)
one = {}
for tag, name, md, mdx, cls, ln in segs:
    if name != "一界面": continue
    for m in range(int(round(md)), int(round(mdx))+1):
        one[m] = cls
def cnt(rng, cls='cyan_solid'):
    return sum(1 for m in rng if one.get(m) == cls)
print("尾管段 5308-7514 良好:", cnt(range(5308,7515)))
print("目的段 7442-7595 良好:", cnt(range(7442,7596)))
print("7514-7595 良好:", cnt(range(7515,7596)))
print("7329-7442 良好:", cnt(range(7329,7443)))
print("7442-7514 良好:", cnt(range(7442,7515)))
print("交叉验证 7329-7442 + 7442-7514 =", cnt(range(7329,7443))+cnt(range(7442,7515)))
