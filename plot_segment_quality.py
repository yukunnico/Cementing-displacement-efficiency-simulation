import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib

# 设置中文字体
matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
matplotlib.rcParams['axes.unicode_minus'] = False

# 读取深度剖面数据
csv_path = Path("results/呼1-003_1D2D耦合模型/呼1-003_1D2D耦合模型_深度剖面.csv")
df = pd.read_csv(csv_path)

# 找到水泥浆前缘位置（从井底往上，水泥平均浓度>0.01的最近位置）
# 由于数据是从井底到井口排列的，找到第一个水泥浓度接近0的位置
frontier_idx = None
for i in range(len(df)):
    if df['水泥平均浓度'].iloc[i] < 0.01:
        frontier_idx = i
        break

if frontier_idx is None:
    frontier_idx = len(df) - 1

# 水泥浆占据的井段范围
bottom_md = df['井深_m'].iloc[0]  # 井底
top_md = df['井深_m'].iloc[frontier_idx]  # 前缘位置

print(f"井底深度: {bottom_md:.1f}m")
print(f"水泥浆前缘深度: {top_md:.1f}m")
print(f"水泥浆占据段长: {bottom_md - top_md:.1f}m")

# 按200m分段
segment_length = 200.0
n_segments = int(np.ceil((bottom_md - top_md) / segment_length))

segments = []
for i in range(n_segments):
    seg_bottom = bottom_md - i * segment_length
    seg_top = max(seg_bottom - segment_length, top_md)
    
    # 提取该段数据
    mask = (df['井深_m'] <= seg_bottom) & (df['井深_m'] > seg_top)
    seg_data = df[mask]
    
    if len(seg_data) > 0:
        avg_cement = seg_data['水泥平均浓度'].mean()
        avg_efficiency = seg_data['平均有效顶替效率'].mean()
        
        # 判断质量等级
        if avg_cement > 0.80:
            quality = "良好"
            color = "#2ecc71"  # 绿色
        elif avg_cement > 0.70:
            quality = "合格"
            color = "#f39c12"  # 橙色
        else:
            quality = "不合格"
            color = "#e74c3c"  # 红色
        
        segments.append({
            '段号': i + 1,
            '井段': f"{seg_bottom:.0f}-{seg_top:.0f}m",
            '段长_m': seg_bottom - seg_top,
            '平均水泥占据率': avg_cement,
            '平均顶替效率': avg_efficiency,
            '质量等级': quality,
            '颜色': color
        })

# 创建DataFrame
seg_df = pd.DataFrame(segments)
print("\n各段顶替质量评价:")
print(seg_df[['段号', '井段', '平均水泥占据率', '质量等级']].to_string(index=False))

# 绘制图表
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), gridspec_kw={'height_ratios': [3, 1]})

# 上图：分段柱状图
y_pos = np.arange(len(seg_df))
bars = ax1.barh(y_pos, seg_df['平均水泥占据率'] * 100, 
                color=seg_df['颜色'], edgecolor='black', linewidth=0.5)

# 添加段号标签
for i, (bar, row) in enumerate(zip(bars, seg_df.itertuples())):
    width = bar.get_width()
    ax1.text(width + 1, bar.get_y() + bar.get_height()/2, 
             f"{row.井段}\n{width:.1f}% ({row.质量等级})", 
             va='center', fontsize=10, fontweight='bold')

# 添加合格线
ax1.axvline(x=70, color='red', linestyle='--', linewidth=2, label='合格线 (70%)')
ax1.axvline(x=80, color='green', linestyle='--', linewidth=2, label='良好线 (80%)')

ax1.set_yticks(y_pos)
ax1.set_yticklabels([f"第{s['段号']}段" for s in segments])
ax1.set_xlabel('水泥浆占据率 (%)', fontsize=12)
ax1.set_title('呼1-003井（HT1-003）尾管段水泥浆顶替质量分段评价\n（每段200m，从井底向上）', fontsize=14, fontweight='bold')
ax1.set_xlim(0, 105)
ax1.legend(loc='lower right')
ax1.grid(axis='x', alpha=0.3)

# 反转y轴，使井底在上方
ax1.invert_yaxis()

# 下图：深度-水泥占据率连续剖面
ax2.plot(df['井深_m'], df['水泥平均浓度'] * 100, 'b-', linewidth=2, label='水泥占据率')
ax2.axhline(y=70, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
ax2.axhline(y=80, color='green', linestyle='--', linewidth=1.5, alpha=0.7)
ax2.fill_between(df['井深_m'], 0, df['水泥平均浓度'] * 100, 
                  where=(df['水泥平均浓度'] > 0.8), alpha=0.3, color='green', label='良好区')
ax2.fill_between(df['井深_m'], 0, df['水泥平均浓度'] * 100, 
                  where=((df['水泥平均浓度'] > 0.7) & (df['水泥平均浓度'] <= 0.8)), 
                  alpha=0.3, color='orange', label='合格区')
ax2.fill_between(df['井深_m'], 0, df['水泥平均浓度'] * 100, 
                  where=(df['水泥平均浓度'] <= 0.7), alpha=0.3, color='red', label='不合格区')

ax2.set_xlabel('井深 (m)', fontsize=12)
ax2.set_ylabel('水泥占据率 (%)', fontsize=12)
ax2.set_title('水泥浆占据率深度剖面', fontsize=12)
ax2.legend(loc='upper left')
ax2.grid(alpha=0.3)
ax2.set_xlim(top_md - 100, bottom_md + 100)

plt.tight_layout()
output_path = Path("results/呼1-003_1D2D耦合模型/呼1-003井_分段顶替质量评价.png")
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"\n图表已保存: {output_path}")

# 保存CSV
seg_df[['段号', '井段', '段长_m', '平均水泥占据率', '平均顶替效率', '质量等级']].to_csv(
    "results/呼1-003_1D2D耦合模型/呼1-003井_分段顶替质量评价.csv", 
    index=False, encoding='utf-8-sig'
)

# 统计
excellent_count = len(seg_df[seg_df['质量等级'] == '良好'])
qualified_count = len(seg_df[seg_df['质量等级'] == '合格'])
unqualified_count = len(seg_df[seg_df['质量等级'] == '不合格'])
total_segments = len(seg_df)

print(f"\n统计结果:")
print(f"  总段数: {total_segments}")
print(f"  良好 (>80%): {excellent_count}段 ({excellent_count/total_segments*100:.1f}%)")
print(f"  合格 (70-80%): {qualified_count}段 ({qualified_count/total_segments*100:.1f}%)")
print(f"  不合格 (<70%): {unqualified_count}段 ({unqualified_count/total_segments*100:.1f}%)")
