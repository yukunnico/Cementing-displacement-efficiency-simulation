# 提取 8 井摘要.json + 汇总.json 的全部指标，输出总表
import json, csv, os, sys
sys.stdout.reconfigure(encoding='utf-8')

BASE = r"D:/users/desktop/research/控压固井项目/cement model/results/胶塞语义修复后终跑_2026-09-03"
WELLS = ["hu1","hu2","hu101","hu102","hu103","ht1_001","ht1_003","ht1_004"]

def full_keys(o, prefix=''):
    keys = []
    if isinstance(o, dict):
        for k, v in o.items():
            p = prefix + '.' + k if prefix else k
            if isinstance(v, dict):
                keys.append((p, type(v).__name__, len(v)))
                keys.extend(full_keys(v, p))
            elif isinstance(v, list):
                keys.append((p, f'list[{len(v)}]', None))
            else:
                keys.append((p, type(v).__name__, v))
    return keys

sums = {}
with open(os.path.join(BASE, '汇总.json'), encoding='utf-8') as f:
    for row in json.load(f):
        sums[row['井']] = row

summaries = {}
for w in WELLS:
    with open(os.path.join(BASE, f'{w}_摘要.json'), encoding='utf-8') as f:
        summaries[w] = json.load(f)

# --- 键名清单（以 hu101 为代表，但检查 8 井键集合差异）---
print('=== 摘要.json 顶层键（hu101）===')
for k, v in summaries['hu101'].items():
    print(f'  {k}: {type(v).__name__}' + (f' len={len(v)}' if isinstance(v,(dict,list)) else f' = {v}'))

def leaf_paths(d, prefix=''):
    out = set()
    if isinstance(d, dict):
        for k, v in d.items():
            p = f'{prefix}.{k}' if prefix else k
            if isinstance(v, dict):
                out |= leaf_paths(v, p)
            else:
                out.add(p)
    return out

all_paths = {}
for w in WELLS:
    all_paths[w] = leaf_paths(summaries[w])
union = set().union(*all_paths.values())
only_hu101 = all_paths['hu101'] - set().union(*[all_paths[w] for w in WELLS if w != 'hu101'])
missing_in_some = {w: union - all_paths[w] for w in WELLS}
print('\n=== 叶子路径总数（每井）===')
for w in WELLS:
    print(f'  {w}: {len(all_paths[w])}')
print('仅某井独有的路径:', {w: sorted(v) for w, v in missing_in_some.items() if v} or '无（8 井叶子路径一致）')

def get_path(d, path):
    cur = d
    for part in path.split('.'):
        if isinstance(cur, dict) and part in cur:
            cur = cur[part]
        else:
            return None
    return cur

INTEREST = [
    'effective_efficiency', 'eta_narrow', 'channeling_index', 'mixing_index',
    '最终结果.最终失稳指数', '最终结果.最终失稳指数_线性', '最终结果.浮力数_b',
    '最终结果.最终水泥浆占据率',
    'tier0_diagnostics.muskat_regime.eccentricity',
    'tier0_diagnostics.muskat_regime.regime',
    'tier0_diagnostics.muskat_regime.viscosity_ratio',
    'tier0_diagnostics.muskat_regime.shock_detected',
    'tier0_diagnostics.flow_classification.flow_class',
    'tier0_diagnostics.flow_classification.w0_m_s',
    'tier0_diagnostics.displacement_metrics.mud_retention_fraction',
    'tier0_diagnostics.displacement_metrics.interface_length_ratio',
    'tier0_diagnostics.displacement_metrics.t_br_s',
    'tier0_diagnostics.displacement_metrics.t_br_hat',
    'tier0_diagnostics.displacement_metrics.eta_narrow',
    'tier0_diagnostics.displacement_metrics.eta_global',
    'tier0_diagnostics.displacement_metrics.mean_flusher',
    'tier0_diagnostics.buoyancy_regime.b_number',
    'tier0_diagnostics.buoyancy_regime.regime',
    '低尾指标.standoff低于0.5段占比',
    '低尾指标.窄边效率低于0.05域占比',
    '物理环空体积_m3', '井段_m', '模拟对象',
]

print('\n=== 井×指标总表（摘要.json）===')
hdr = ['井'] + [i.split('.')[-1] if i.count('.')>=1 else i for i in INTEREST]
print('\t'.join(hdr))
rows = {}
for w in WELLS:
    d = summaries[w]
    vals = [get_path(d, i) for i in INTEREST]
    rows[w] = dict(zip(INTEREST, vals))
    print('\t'.join([w] + [(f'{v:.4g}' if isinstance(v, float) else str(v)) for v in vals]))

# 评价窗效率
print('\n=== 评价窗效率 ===')
win_names = set()
for w in WELLS:
    win_names |= set(summaries[w].get('评价窗效率', {}).keys())
win_names = sorted(win_names)
print('窗口名:', win_names)
print('井\t' + '\t'.join(f'{n}' for n in win_names))
for w in WELLS:
    vals = []
    for n in win_names:
        info = summaries[w].get('评价窗效率', {}).get(n)
        if info is None:
            vals.append('N/A')
        else:
            vals.append(f"E={info['eta_E']:.4g},N={info['eta_N']:.4g}")
    print(w + '\t' + '\t'.join(vals))

# 汇总.json 全列
print('\n=== 汇总.json 全列 ===')
cols = list(sums[WELLS[0]].keys())
print('\t'.join(cols))
for w in WELLS:
    print('\t'.join(str(sums[w].get(c, '')) for c in cols))

# 存成 json 供后续脚本用
out = {'rows': rows, 'windows': {w: summaries[w].get('评价窗效率', {}) for w in WELLS},
       'sums': sums}
with open(r"D:/users/desktop/research/控压固井项目/cement model/.tmp_research/extracted_metrics.json", 'w', encoding='utf-8') as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print('\nsaved extracted_metrics.json')
