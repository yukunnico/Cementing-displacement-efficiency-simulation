# -*- coding: utf-8 -*-
"""汇总 interface_segments.json → 各评价窗口合格率 + 生成逐深度 CSVD。"""
import json, os, csv

OUT = r"D:\users\desktop\research\控压固井项目\cement model\参考文档\现场资料提取\ht1_004_呼1-004\cbl_digitization"
WORK = os.path.join(OUT, "_work")
rows = json.load(open(os.path.join(WORK, "interface_segments.json"), encoding="utf-8"))
# rows: (tag, interface, top_m, bot_m, class_key, len)

def frac(pairs, a, b):
    """pairs: (top,bot,cls_key)；返回区间[a,b]内长度分布与中等及以上占比"""
    length = {'良好':0.0,'中等':0.0,'差':0.0,'空':0.0,'不评':0.0}
    clsmap = {'cyan_solid':'良好','cyan_strip':'中等','green':'差','empty':'空'}
    tot = 0.0; good = 0.0; evald = 0.0
    for top, bot, k in pairs:
        md = max(top, a); mx = min(bot, b)
        if mx <= md: continue
        ln = mx - md
        c = clsmap.get(k, '空')
        length[c] += ln
        tot += ln
        if c in ('良好','中等'):
            good += ln; evald += ln
        elif c == '差':
            evald += ln
    return length, tot, good, evald

def rate(pairs, a, b):
    length, tot, good, evald = frac(pairs, a, b)
    r = (good/evald*100) if evald>0 else None
    return length, tot, good, evald, r

seg = {'tail':(7382,7581), 'mid':(5249,7384), 'shallow':(11,5251)}
inter1 = [ (t,b,k) for (tag,iface,t,b,k,l) in rows if iface=='一界面']
inter2 = [ (t,b,k) for (tag,iface,t,b,k,l) in rows if iface=='二界面']

WINDOWS = [
    ("全井一界面 11-7581", inter1, 11, 7581),
    ("尾管评价段一界面 5245-7581", inter1, 5245, 7581),
    ("尾管段一界面 5316-7581(任务初估)", inter1, 5316, 7581),
    ("目的井段一界面 7495-7550", inter1, 7495, 7550),
    ("油气水层段一界面 7482-7560", inter1, 7482, 7560),
    ("全井二界面 5578-7581", inter2, 5578, 7581),
    ("尾管评价段二界面 5245-7581", inter2, 5245, 7581),
]
print("=== HT1-004 CBL 合格率估计（中等及以上占比，基于矢量解释填充）===")
for name, pairs, a, b in WINDOWS:
    length, tot, good, evald, r = rate(pairs, a, b)
    rstr = f"{r:.1f}%" if r is not None else "N/A(无评价)"
    print(f"{name:28s} 区间[{a}-{b}m] 评价长度={evald:.0f}m 良好={length['良好']:.0f}m 中等={length['中等']:.0f}m 差={length['差']:.0f}m → 中等及以上={good:.0f}m = {rstr}")

# 一界面逐 500m 分段
print("\n=== 一界面每500m分段 中等及以上占比 ===")
for a in range(0, 8000, 500):
    b = a+500
    if b <= 11: continue
    length, tot, good, evald, r = rate(inter1, a, b)
    if evald > 0 and r is not None:
        print(f"  {a:5d}-{b:5d}m  评价={evald:6.0f}m  良好={length['良好']:5.0f} 中等={length['中等']:5.0f} 差={length['差']:6.0f} → {r:5.1f}%")

# 导出逐 1m 深度序列 CSV（解释结论）
print("\n=== 导出 cbl_digitization.csv（逐米一界面解释结论）===")
out_rows = []
md = 11.0
while md < 7581:
    # 一界面
    cls_i = '空'
    for top, bot, k in inter1:
        if top <= md < bot:
            cls_i = {'cyan_solid':'良好','cyan_strip':'中等','green':'差','empty':'空'}[k]
            break
    cls_ii = '不评价'
    for top, bot, k in inter2:
        if top <= md < bot:
            cls_ii = {'cyan_solid':'良好','cyan_strip':'中等','green':'差','empty':'空'}[k]
            break
    out_rows.append((md, cls_i, cls_ii))
    md += 1.0

# 写 CSV（每米一行）
csvpath = os.path.join(OUT, "cbl_digitization.csv")
with open(csvpath, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["well_id","md_m","amp_pct","quality_class","data_type","source_description","confidence","notes"])
    for md, ci, cii in out_rows:
        q = ci if ci != '空' else ('空套管' if md < 7382 else '未评价')
        amp = ""  # 声幅未从矢量曲线提取（见 notes）
        w.writerow(["ht1_004", f"{md:.1f}", amp, q, "interpreted",
                    "HT1-004_固井质量测井评价图 图头+一界面解释填充(矢量) 20260721 完井",
                    "medium", f"一界面{ci};二界面{cii}"])
print("written", csvpath, "rows", len(out_rows))
