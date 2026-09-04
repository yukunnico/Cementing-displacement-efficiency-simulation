# 深度剖面窄边残留分析 + hu101 npz 窄边通道证据 + CBL 窗 1e-13 搜索
import json, sys
import numpy as np
import pandas as pd
sys.stdout.reconfigure(encoding='utf-8')
BASE = r"D:/users/desktop/research/控压固井项目/cement model/results/胶塞语义修复后终跑_2026-09-03"
WELLS = ["hu1","hu2","hu101","hu102","hu103","ht1_001","ht1_003","ht1_004"]

# ---- 1. 逐深度 e（偏心度指标）触顶检查 ----
print('=== 深度剖面 e（偏心度指标）分布 ===')
print('井\te域均\t触顶0.55段占比\te>0.5占比\t最小e\t最大e\tstandoff<0.5占比(独立复算)')
for w in WELLS:
    df = pd.read_csv(rf"{BASE}/{w}_深度剖面.csv")
    e = df['偏心度指标'].to_numpy()
    so = df['居中度'].to_numpy()
    print(f"{w}\t{e.mean():.4f}\t{(e >= 0.55 - 1e-9).mean():.3f}\t{(e > 0.5).mean():.3f}\t{e.min():.3f}\t{e.max():.3f}\t{(so < 0.5).mean():.3f}")

# ---- 2. hu101 / hu1 / hu102 窄边残留沿深度 ----
for w in ['hu101', 'hu1', 'hu102']:
    df = pd.read_csv(rf"{BASE}/{w}_深度剖面.csv")
    n = len(df)
    bands = [('底部0-10%', 0, int(n*0.1)), ('10-50%', int(n*0.1), int(n*0.5)),
             ('50-90%', int(n*0.5), int(n*0.9)), ('顶部90-100%', int(n*0.9), n)]
    print(f'\n=== {w} 深度分带（窄边 vs 宽边）===')
    print('带\t窄边效率均\t宽边效率均\t窄边水泥浓度均\t宽边水泥浓度均\t钻井液浓度均\t环空间隙m(均)')
    for name, a, b in bands:
        seg = df.iloc[a:b]
        print(f"{name}\t{seg['窄边有效效率'].mean():.3f}\t{seg['宽边有效效率'].mean():.3f}\t"
              f"{seg['窄边水泥浓度'].mean():.3f}\t{seg['宽边水泥浓度'].mean():.3f}\t"
              f"{seg['钻井液平均浓度'].mean():.3f}\t{seg['环空间隙_m'].mean():.4f}")
    # 窄边效率 < 0.05 的深度点占比（复算低尾指标口径）
    print(f"窄边效率<0.05 深度点占比: {(df['窄边有效效率'] < 0.05).mean():.4f}；"
          f"窄边效率中位数: {df['窄边有效效率'].median():.4f}")

# ---- 3. hu101 npz 窄边通道证据 ----
print('\n=== hu101 2D场数据.npz 窄边通道证据 ===')
d = np.load(rf"{BASE}/hu101_2D场数据.npz", allow_pickle=True)
y, md = d['y'], d['md']
wf = d['wall_final']; cf = d['cement_final']
print('y 首末:', y[0], y[-1], '单调增:', bool(np.all(np.diff(y) > 0)))
# 判断哪端是窄边：深度剖面 窄边水泥浓度 最低的深度附近，看 wall_final 高值集中在 y 的哪端
dp = pd.read_csv(rf"{BASE}/hu101_深度剖面.csv")
# 深度剖面 251 点 vs md 250 点：取前 250 对齐（方向需核对）
mdc = dp['井深_m'].to_numpy()
print('深度剖面井深首末:', mdc[0], mdc[-1], ' npz md 首末:', md[0], md[-1])
# y 两端的 wall 浓度（最终时刻）
print('wall_final y=0 端均值: %.4f | y=末端均值: %.4f' % (wf[:, 0].mean(), wf[:, -1].mean()))
print('cement_final y=0 端均值: %.4f | y=末端均值: %.4f' % (cf[:, 0].mean(), cf[:, -1].mean()))
# 中部泥浆通道：wall>0.5 的 (y,md) 网格占比，按 y 分带
w5 = (wf > 0.5).mean(axis=1)  # 每条 y 线上 wall>0.5 的深度占比
print('wall>0.5 沿 y 的深度占比（y 分 4 带）:')
for i, name in enumerate(['y带0(端A)', 'y带1', 'y带2', 'y带3(端B)']):
    sl = slice(i*10, (i+1)*10)
    print(f'  {name}: {w5[sl].mean():.4f}')
# 窄边判定后给出：若 y 末端 wall 高 → 末端为窄边
# 泥浆主要残留在哪条 y：argmax
iy = int(np.argmax(w5))
print(f'wall>0.5 占比最高的 y 线索引={iy}/40, y值={y[iy]:.4f}, 占比={w5[iy]:.4f}')

# ---- 4. CBL 窗 1e-13 搜索 ----
print('\n=== 评价窗效率中的机器精度数字（|eta|<1e-6 且 !=0）===')
found = []
for w in WELLS:
    with open(rf"{BASE}/{w}_摘要.json", encoding='utf-8') as f:
        d2 = json.load(f)
    for name, info in d2.get('评价窗效率', {}).items():
        for k in ('eta_E', 'eta_N'):
            v = info.get(k)
            if isinstance(v, float) and 0 < abs(v) < 1e-6:
                found.append((w, name, k, v))
for f_ in found:
    print(f_)
print('共', len(found), '条机器精度数字')
# 全摘要文本搜索 e-13
import re
for w in WELLS:
    txt = open(rf"{BASE}/{w}_摘要.json", encoding='utf-8').read()
    hits = re.findall(r'\d?\.?\d*[eE]-1[0-9]', txt)
    if hits:
        print(w, '含小数字:', set(hits))
