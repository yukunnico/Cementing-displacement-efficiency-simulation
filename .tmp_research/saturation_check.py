# 收尾计算：scipy Spearman p 值 + η_E 恒等式 + 漂移结构 + 入库完成度
import json, sys
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
BASE = r"D:/users/desktop/research/控压固井项目/cement model/results/胶塞语义修复后终跑_2026-09-03"
WELLS = ["hu1","hu2","hu101","hu102","hu103","ht1_001","ht1_003","ht1_004"]

from scipy.stats import spearmanr
print('scipy OK')

merged = json.load(open(r"D:/users/desktop/research/控压固井项目/cement model/.tmp_research/merged_metrics.json", encoding='utf-8'))
# 阶段1 json 的设计水泥量（与 09-03 build_case 同源）
stage1 = json.load(open(r"D:/users/desktop/research/控压固井项目/cement model/results/_质量平衡取证_2026-09-02/阶段1_库存核算.json", encoding='utf-8'))

print('\n=== η_E 恒等式核验：η_E ?= 域内水泥/物理环空体积 ===')
print('井\tη_E\t域内/环空\t偏差\t入库完成度(09-03入环空/设计)\t漂移/环空')
identity_err = {}
for w in WELLS:
    r = merged[w]
    ratio = r['cement_in_domain'] / r['annulus_vol']
    design = stage1[w]['设计水泥量_m3']
    intake_frac = r['cement_in'] / design
    drift_frac = r['field_drift'] / r['annulus_vol']
    identity_err[w] = abs(r['eta_E'] - ratio)
    print(f"{w}\t{r['eta_E']:.4f}\t{ratio:.4f}\t{abs(r['eta_E']-ratio):.2e}\t{intake_frac:.4f}\t{drift_frac:.4f}")
print('恒等式最大偏差: %.2e' % max(identity_err.values()))

print('\n=== Spearman p 值 ===')
def corr(wl, xk, yk):
    xs = [merged[w][xk] for w in wl]
    ys = [merged[w][yk] for w in wl]
    return spearmanr(xs, ys)
med6 = ['hu1','hu2','hu103','ht1_001','ht1_003','ht1_004']
med7 = med6 + ['hu102']
for gname, wl in [('6井中等偏心', med6), ('7井(除hu101)', med7), ('8井全体', WELLS)]:
    for yk in ['eta_E', 'eta_N', 'narrow_arrival', 'narrow_lt005']:
        res = corr(wl, 'e', yk)
        print(f'{gname}: e vs {yk}: rho={res.statistic:.3f} p={res.pvalue:.3f}')
# 漂移/环空 vs η_E（7 井）
res = corr(med7, 'field_drift', 'eta_E')
print(f"7井: 漂移 vs η_E: rho={res.statistic:.3f} p={res.pvalue:.3f}")
res = corr(WELLS, 'inv_ratio', 'field_drift')
print(f"8井: 库存比 vs 漂移: rho={res.statistic:.3f} p={res.pvalue:.3f}")

print('\n=== 潜在窄边残留井识别（e 中等偏高 且 库存比低）===')
cands = sorted([w for w in WELLS if w != 'hu101'], key=lambda w: (merged[w]['e'], -merged[w]['inv_ratio']), reverse=True)
for w in cands:
    r = merged[w]
    print(f"{w}: e={r['e']:.3f} 库存比={r['inv_ratio']:.3f} η_E={r['eta_E']:.4f} η_N={r['eta_N']:.4f} 窄边<0.05域占比={r['narrow_lt005']:.4f}")

print('\n=== hu101 时间序列：窄边前缘推进与停泵时刻 ===')
import pandas as pd
ts = pd.read_csv(rf"{BASE}/hu101_时间序列.csv")
last = ts.iloc[-1]
print('末行: t=%.1fs stage=%s front_wide=%.1f front_narrow=%.1f eta_E=%.4f' % (
    last['time_s'], last['stage'], last['front_wide_m'], last['front_narrow_m'], last['effective_efficiency']))
# 窄边前缘何时开始停滞
fn = ts['front_narrow_m'].to_numpy()
tt = ts['time_s'].to_numpy()
print('窄边前缘最终值 %.1f m（域长 2468m，到位率 %.4f）' % (fn[-1], fn[-1]/2468))
# 前缘到达 >49m 的最后时刻
idx = np.argmax(fn > 49.6)
print('窄边前缘超过最终值的时间点: %.1fs' % tt[idx] if fn.max() > 49.5 else 'n/a')
