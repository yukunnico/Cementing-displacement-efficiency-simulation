# -*- coding: utf-8 -*-
"""HT1-003 CBL 数字化汇总：interface_segments.json -> 逐米CSV + 各口径合格率统计"""
import json, csv, os

WORK = r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\ht1_003_呼1-003\_work"
OUT_DIR = r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\ht1_003_呼1-003"

with open(os.path.join(WORK, "interface_segments.json"), encoding="utf-8") as f:
    segs = json.load(f)  # (tag, name, md, mdx, cls, len)

def cls_cn(c):
    return {'green':'差','cyan_strip':'中等','cyan_solid':'良好','empty':'空套管'}[c]

# 逐米: 一界面 class (m -> cls)，notes 附二界面
# 主图 main: 3500-7595, 浅层 shallow: 15-3500
one = {}   # m -> cls  (一界面)
two = {}   # m -> cls  (二界面)
for tag, name, md, mdx, cls, ln in segs:
    m0 = int(round(md)); m1 = int(round(mdx))
    for m in range(m0, m1+1):
        if name == "一界面":
            one[m] = cls
        else:
            two[m] = cls

# 覆盖深度
dmin = 15; dmax = 7595
rows = []
for m in range(dmin, dmax+1):
    c1 = one.get(m, 'empty')
    c2 = two.get(m, None)  # None = 不评价
    if m < 5568:
        note2 = "二界面不评价(多层套管)"
    else:
        note2 = f"二界面{cls_cn(c2)}" if c2 else "二界面未评价"
    notes = f"一界面{cls_cn(c1)};{note2}"
    rows.append({
        'well_id':'ht1_003','md_m':f"{m}.0",'amp_pct':'',
        'quality_class':cls_cn(c1),'data_type':'interpreted',
        'source_description':'HT1-003_固井质量测井评价图(15-3500/3500-7595) 图头+一界面解释填充(矢量LEAD5.0) 20260622 完井',
        'confidence':'medium','notes':notes
    })

csv_path = os.path.join(OUT_DIR, "cbl_digitization.csv")
with open(csv_path, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=['well_id','md_m','amp_pct','quality_class','data_type','source_description','confidence','notes'])
    w.writeheader()
    w.writerows(rows)
print(f"CSV written: {csv_path}  rows={len(rows)}  cover={dmin}-{dmax}m")

# ===== 统计 =====
def stat(depths_list, one_map, two_map=None, label=""):
    """给定深度列表，统计一界面/二界面各类别长度与中等及以上占比"""
    res = {}
    for key, mp in [("一界面", one_map)]:
        cnt = {'良好':0,'中等':0,'差':0,'空套管':0}
        for m in depths_list:
            cnt[cls_cn(mp.get(m,'empty'))] += 1
        tot = sum(cnt.values()); mid = cnt['中等']+cnt['良好']
        res[key] = (cnt, tot, mid)
        print(f"  [{label}] {key}: 良好={cnt['良好']}m 中等={cnt['中等']}m 差={cnt['差']}m 空={cnt['空套管']}m | 评价={tot}m 中等及以上={mid}m ({mid/tot*100:.1f}%)" if tot else f"  [{label}] {key}: 无数据")
    return res

ranges = {
    "全井测量段 3500-7595": range(3500, 7596),
    "尾管评价段 5307.54-7514.21(水泥面-人工井底)": range(5308, 7515),
    "尾管段 5316.036-7595": range(5316, 7596),
    "目的井段 7442-7595(图中)": range(7442, 7596),
    "油气水层段标注 7353-7560": range(7353, 7561),
}
print("\n===== 合格率统计（中等及以上占比） =====")
for label, rng in ranges.items():
    stat(list(rng), one, two, label)

# 油气水层段逐个（黄色标注 8 段）
oil_gas = [(7352.9,7357.7),(7380.2,7383.2),(7398.0,7407.4),(7483.2,7489.8),
           (7499.4,7502.8),(7504.0,7508.6),(7526.6,7529.2),(7548.4,7559.6)]
print("\n===== 油气水层段逐段一界面质量 =====")
for a, b in oil_gas:
    ms = list(range(int(round(a)), int(round(b))+1))
    cnt = {'良好':0,'中等':0,'差':0,'空套管':0}
    for m in ms:
        cnt[cls_cn(one.get(m,'empty'))] += 1
    mid = cnt['中等']+cnt['良好']
    print(f"  {a:.1f}-{b:.1f}m: 良好={cnt['良好']} 中等={cnt['中等']} 差={cnt['差']} | 中等及以上={mid}/{len(ms)} ({mid/len(ms)*100:.0f}%)")
