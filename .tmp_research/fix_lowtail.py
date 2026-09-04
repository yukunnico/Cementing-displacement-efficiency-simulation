# 修正低尾指标提取 + 居中度-效率关系量化（Spearman）
import json, sys
sys.stdout.reconfigure(encoding='utf-8')
BASE = r"D:/users/desktop/research/控压固井项目/cement model/results/胶塞语义修复后Tab终跑_2026-09-03"
BASE = r"D:/users/desktop/research/控压固井项目/cement model/results/胶塞语义修复后终跑_2026-09-03"
WELLS = ["hu1","hu2","hu101","hu102","hu103","ht1_001","ht1_003","ht1_004"]

summaries = {}
for w in WELLS:
    with open(rf"{BASE}/{w}_摘要.json", encoding='utf-8') as f:
        summaries[w] = json.load(f)

print('=== 低尾指标（直接键访问）===')
for w in WELLS:
    lw = summaries[w].get('低尾指标', {})
    print(w, lw)

print('\n=== 逐井核对的会话结论数字 ===')
print('井\te\t库存比\tη_E\tη_N\t窄边到位率\t窜槽\tstandoff<0.5占比\t窄边<0.05域占比\t守恒率\t场漂移_m3\t域长m')
for w in WELLS:
    d = summaries[w]
    e = d['tier0_diagnostics']['muskat_regime']['eccentricity']
    lw = d.get('低尾指标', {})
    print(f"{w}\t{e:.4f}\t{sums_ := d.get('库存比', None)}\t{d['effective_efficiency']:.4f}\t{d['eta_narrow']:.4f}\t"
          f"{lw.get('窄边效率低于0.05域占比', float('nan')):.4f}")
