# 修正低尾指标提取 + 居中度-效率关系量化（Spearman）+ 库存比-η_E 结构
import json, sys
sys.stdout.reconfigure(encoding='utf-8')
BASE = r"D:/users/desktop/research/控压固井项目/cement model/results/胶塞语义修复后终跑_2026-09-03"
WELLS = ["hu1","hu2","hu101","hu102","hu103","ht1_001","ht1_003","ht1_004"]

summaries = {}
for w in WELLS:
    with open(rf"{BASE}/{w}_摘要.json", encoding='utf-8') as f:
        summaries[w] = json.load(f)
with open(rf"{BASE}/汇总.json", encoding='utf-8') as f:
    sums = {r['井']: r for r in json.load(f)}

print('=== 低尾指标（直接键访问）===')
for w in WELLS:
    lw = summaries[w].get('低尾指标', {})
    print(w, json.dumps(lw, ensure_ascii=False))

print('\n=== 合并表（摘要+汇总）===')
hdr = ['井','e','库存比','η_E','η_N','窜槽','混浆','窄边到位率','宽边到位率','standoff<0.5占比','窄边<0.05域占比','守恒率','场漂移_m3','域长_m','物理环空_m3','域满_m3','入环空水泥_m3','域内水泥_m3','2D入库积分_m3','stop_s','stop_B3前_s']
print('\t'.join(hdr))
data = {}
for w in WELLS:
    d = summaries[w]
    s = sums[w]
    lw = d.get('低尾指标', {})
    rec = dict(
        e=d['tier0_diagnostics']['muskat_regime']['eccentricity'],
        inv_ratio=s['库存比'],
        eta_E=d['effective_efficiency'],
        eta_N=d['eta_narrow'],
        channeling=d['channeling_index'],
        mixing=d['mixing_index'],
        narrow_arrival=s['窄边到位率'],
        wide_arrival=s['宽边到位率'],
        standoff_lt05=lw.get('standoff低于0.5段占比'),
        narrow_lt005=lw.get('窄边效率低于0.05域占比'),
        conservation=s['守恒率'],
        field_drift=s['场体积漂移_m3'],
        domain_len=s['域长m'],
        annulus_vol=s['物理环空体积_m3'],
        domain_full=s['域满体积_m3'],
        cement_in=s['入环空水泥_m3'],
        cement_in_domain=s['域内水泥_m3'],
        intake_2d=s['2D入库积分_m3'],
        stop_s=s['stop_s'],
        stop_preB3=s['stop_B3前_s'],
    )
    data[w] = rec
    print('\t'.join([w] + [f'{rec[k]:.4g}' if isinstance(rec[k], float) else str(rec[k]) for k in list(rec)[1:]]))

with open(r"D:/users/desktop/research/控压固井项目/cement model/.tmp_research/merged_metrics.json", 'w', encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=1)

# ---- Spearman ----
try:
    from scipy.stats import spearmanr
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False

def rank(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0]*len(v)
    i = 0
    while i < len(v):
        j = i
        while j+1 < len(v) and v[order[j+1]] == v[order[i]]:
            j += 1
        avg = (i + j)/2 + 1
        for k in range(i, j+1):
            r[order[k]] = avg
        i = j+1
    return r

def pearson(a, b):
    n = len(a)
    ma, mb = sum(a)/n, sum(b)/n
    num = sum((x-ma)*(y-mb) for x, y in zip(a, b))
    da = sum((x-ma)**2 for x in a)**0.5
    db = sum((y-mb)**2 for y in b)**0.5
    return num/(da*db) if da*db > 0 else float('nan')

def spearman(x, y):
    return pearson(rank(x), rank(y))

med = [w for w in WELLS if w != 'hu101']
med6 = ['hu1','hu2','hu103','ht1_001','ht1_003','ht1_004']  # e 0.17-0.35，不含 hu101/hu102
print('\n=== Spearman 相关（n 很小，仅供参考）===')
groups = {'7井(除hu101)': med, '6井中等偏心(e0.17-0.35)': med6, '8井全体': WELLS}
for gname, wl in groups.items():
    es = [data[w]['e'] for w in wl]
    for mk in ['eta_E','eta_N','narrow_arrival','narrow_lt005','channeling']:
        vs = [data[w][mk] for w in wl]
        rho = spearman(es, vs)
        print(f'{gname}: e vs {mk}: rho={rho:.3f} (n={len(wl)})')

print('\n=== 7 井（除hu101）η_E/η_N 极差 ===')
for mk in ['eta_E','eta_N','narrow_arrival','narrow_lt005']:
    vals = sorted((data[w][mk], w) for w in med)
    print(f'{mk}: min={vals[0][0]:.4f}({vals[0][1]}) max={vals[-1][0]:.4f}({vals[-1][1]}) range={vals[-1][0]-vals[0][0]:.4f}')

print('\n=== 库存比 vs η_E / 场漂移（饱和结构检验）===')
print('井\t库存比\tη_E\tη_N\t场漂移_m3\t窄边<0.05域占比')
for w in sorted(WELLS, key=lambda w: data[w]['inv_ratio']):
    r = data[w]
    print(f"{w}\t{r['inv_ratio']:.3f}\t{r['eta_E']:.4f}\t{r['eta_N']:.4f}\t{r['field_drift']:.2f}\t{r['narrow_lt005']:.4f}")

print('\n=== e vs 窄边指标 排序核对 ===')
for w in sorted(WELLS, key=lambda w: data[w]['e']):
    r = data[w]
    print(f"e={r['e']:.3f} {w}\tη_E={r['eta_E']:.4f}\tη_N={r['eta_N']:.4f}\t窄边到位={r['narrow_arrival']:.4f}\t窄边<0.05占比={r['narrow_lt005']:.4f}")
