"""定位 σ_off=0.0 band：84s 处的 band 是 REAR_EXIT 同刻残片"""
import sys; sys.path.insert(0, 'tests')
from test_casing_mixing_contact_time import _dispersion_bands
# flag=False 84s band: center=84, sigma=0 —— 5子事件中 t_sub 被 max(...,0.0) 压缩？
# 检查 80.81→87.19 band: sigma=0.797 ✓ 真过渡带。84s 的 sigma=0 band：
# 事件序: 80.8132, 82.4066, [REAR_EXIT 84], 84 FRONT(尾浆0.5), [RATE_SWITCH 84], 85.5934, 87.1868, 87.1868
# _dispersion_bands 扫描：遇 80.8132 frac<1 启动 band，while 扫到 85.5934（中间夹 REAR_EXIT/RATE_SWITCH
# 不匹配 kind==FRONT_ARRIVAL 条件 → j 停在 82.4066 后的 84.0 FRONT? 不，while 条件要求连续 FRONT_ARRIVAL
# 且同相名—— REAR_EXIT@84 打断连续性 → band1=[80.8132,82.4066] sigma=(82.4066-80.8132)/2=0.797?
# 但输出显示 band center=81.6099 sigma=0.797 ✓（80.8132..82.4066 跨度1.5934/2=0.797——这其实只扫到2个子事件）
# band2: 84.0 frac=0.5<1 启动 → while 扫 [84.0] 只有1个（85.5934 也匹配！）→ 84.0..87.1868? 
# 实际 band2 sigma=0 center=84：84..84? 
print("band 提取器缺陷：REAR_EXIT/RATE_SWITCH 事件夹在过渡带中间打断连续段，")
print("导致同一物理过渡带被切成多个 band，sigma 提取错误。")
print("且 84.0 处 frac=0.5 子事件与 80.8132 处 frac=0.0786 子事件属同一 5 子事件序列,")
print("但被 REAR_EXIT@84 切断。")
