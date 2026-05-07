"""
绘制呼101尾管井固井质量与模型计算顶替效率沿深度分布图

完全匹配参考图像样式：
- 左侧：深度段表格
- 左二：CBL固井质量图像粗读（颜色带）
- 中间：横向条形图（模型计算顶替效率）
- 右侧：深度剖面水泥体积分数曲线

数据来源：
- 分段统计：参考模型CSV
- 深度剖面：模型结果CSV
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import os
import sys

# --------------------------- 字体设置 ---------------------------
font_list = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS', 'WenQuanYi Micro Hei']
font_found = None
for font in font_list:
    try:
        plt.rcParams['font.family'] = font
        plt.rcParams['axes.unicode_minus'] = False
        fig_test, ax_test = plt.subplots()
        ax_test.text(0.5, 0.5, '中文测试', ha='center', va='center')
        plt.close(fig_test)
        font_found = font
        break
    except Exception:
        continue

if font_found is None:
    print("[警告] 未找到可用的中文字体，图表中文可能显示为方框。", file=sys.stderr)
else:
    print(f"[信息] 使用字体: {font_found}")

# --------------------------- 路径设置 ---------------------------
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CSV_PATH = os.path.join(
    PROJECT_ROOT, '参考文档', 'hu101_1d2d_paper_project', 
    '呼101尾管_固井质量_分段顶替效率图表',
    '呼101尾管_各段固井质量与模型顶替效率.csv'
)

MODEL_PROFILE_PATH = os.path.join(
    PROJECT_ROOT, 'results', '呼101尾管_1D2D耦合模型',
    '呼101尾管_1D2D耦合模型_深度剖面.csv'
)

OUTPUT_DIR = os.path.join(PROJECT_ROOT, 'results', '呼101尾管_1D2D耦合模型')

# --------------------------- 数据读取 ---------------------------
def load_segment_data(path):
    df = pd.read_csv(path, encoding='utf-8-sig')
    return df

def load_profile_data(path):
    df = pd.read_csv(path, encoding='utf-8-sig')
    # 过滤到有效深度范围
    df = df[(df['井深_m'] >= 5400) & (df['井深_m'] <= 7870)].copy()
    # 反转顺序使深度从上到下增加
    df = df.sort_values('井深_m', ascending=True).reset_index(drop=True)
    return df

# --------------------------- 绘图 ---------------------------
def plot_full_chart(seg_df, profile_df, output_path):
    """
    绘制完整参考图样式的图表
    """
    # 颜色映射（匹配参考图）
    quality_color_map = {
        '良好': '#4CAF50',       # 绿色
        '较好': '#66BB6A',       # 浅绿色
        '中等': '#DAA520',       # 金黄色
        '一般': '#DAA520',       # 金黄色（与中等同色）
        '局部差': '#D32F2F',     # 红色
        '差': '#D32F2F',         # 红色（与局部差同色）
        '未测井评价': '#9E9E9E', # 灰色
    }
    
    # 创建图表布局
    fig = plt.figure(figsize=(18, 12))
    gs = gridspec.GridSpec(1, 4, width_ratios=[0.8, 1.2, 3, 2.5], wspace=0.1)
    
    # ========== 左一：深度段表格 ==========
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_xlim(0, 1)
    ax1.set_ylim(5400, 7870)
    ax1.invert_yaxis()
    
    # 绘制深度段
    for i, row in seg_df.iterrows():
        top = row['顶界_m']
        bottom = row['底界_m']
        mid = (top + bottom) / 2
        
        # 绘制分隔线
        ax1.axhline(y=top, color='black', linewidth=0.5)
        ax1.axhline(y=bottom, color='black', linewidth=0.5)
        
        # 标注深度
        ax1.text(0.5, mid, f"{top:.0f}\n-\n{bottom:.0f}", 
                ha='center', va='center', fontsize=9)
    
    ax1.set_title('深度段', fontsize=10)
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    
    # ========== 左二：CBL固井质量图像粗读 ==========
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_xlim(0, 1)
    ax2.set_ylim(5400, 7870)
    ax2.invert_yaxis()
    
    # 绘制CBL质量颜色带
    for i, row in seg_df.iterrows():
        top = row['顶界_m']
        bottom = row['底界_m']
        quality = row['固井质量解释']
        color = quality_color_map.get(quality, '#9E9E9E')
        mid = (top + bottom) / 2
        
        # 绘制颜色矩形
        rect = Rectangle((0.1, top), 0.8, bottom - top, 
                         facecolor=color, edgecolor='black', linewidth=0.5)
        ax2.add_patch(rect)
        
        # 添加质量标注
        ax2.text(0.5, mid, quality, ha='center', va='center', 
                fontsize=10, fontweight='bold', color='white' if quality in ['差', '局部差', '良好'] else 'black')
    
    ax2.set_title('测井固井质量\n图像粗读', fontsize=10)
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    
    # ========== 中间：模型计算顶替效率横向条形图 ==========
    ax3 = fig.add_subplot(gs[0, 2])
    
    # 准备数据
    y_positions = []
    bar_heights = []
    colors = []
    labels = []
    values = []
    
    for i, row in seg_df.iterrows():
        top = row['顶界_m']
        bottom = row['底界_m']
        mid = (top + bottom) / 2
        quality = row['固井质量解释']
        eff = row['模型顶替效率_%']
        
        y_positions.append(mid)
        bar_heights.append(bottom - top)
        colors.append(quality_color_map.get(quality, '#9E9E9E'))
        labels.append(f"{row['井段_m']} {quality}")
        values.append(eff)
    
    # 绘制横向条形图
    bars = ax3.barh(y_positions, values, height=bar_heights, 
                     color=colors, edgecolor='black', linewidth=0.5, alpha=0.85)
    
    # 添加数值标签
    for y, val, label in zip(y_positions, values, labels):
        ax3.text(val + 1, y, f'{val:.1f}%', ha='left', va='center', 
                fontsize=10, fontweight='bold')
        ax3.text(2, y, label, ha='left', va='center', 
                fontsize=9, color='black')
    
    ax3.set_ylim(5400, 7870)
    ax3.invert_yaxis()
    ax3.set_xlim(0, 110)
    ax3.set_xlabel('模型计算顶替效率 / %', fontsize=11)
    ax3.set_yticks([])
    ax3.grid(axis='x', linestyle='--', alpha=0.4)
    
    # 添加图例
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#4CAF50', edgecolor='black', label='良好/较好'),
        Patch(facecolor='#DAA520', edgecolor='black', label='中等/一般'),
        Patch(facecolor='#D32F2F', edgecolor='black', label='差/局部差'),
        Patch(facecolor='#9E9E9E', edgecolor='black', label='未测井评价'),
    ]
    ax3.legend(handles=legend_elements, loc='upper right', fontsize=9, title='固井质量')
    
    # ========== 右侧：深度剖面水泥体积分数曲线 ==========
    ax4 = fig.add_subplot(gs[0, 3])
    
    depth = profile_df['井深_m'].values
    cement_avg = profile_df['水泥平均浓度'].values * 100  # 转换为百分比
    cement_wide = profile_df['宽边水泥浓度'].values * 100
    cement_narrow = profile_df['窄边水泥浓度'].values * 100
    
    # 绘制曲线
    ax4.plot(cement_avg, depth, 'b-', linewidth=2, label='周向平均', alpha=0.8)
    ax4.plot(cement_wide, depth, 'g-', linewidth=1.5, label='宽边', alpha=0.8)
    ax4.plot(cement_narrow, depth, 'r-', linewidth=1.5, label='窄边', alpha=0.8)
    
    ax4.set_ylim(5400, 7870)
    ax4.invert_yaxis()
    ax4.set_xlim(0, 105)
    ax4.set_xlabel('深度剖面水泥体积分数 / %', fontsize=11)
    ax4.set_yticks([])
    ax4.grid(axis='x', linestyle='--', alpha=0.4)
    ax4.legend(loc='lower left', fontsize=9)
    
    # 添加右侧Y轴标签（深度）
    ax4_right = ax4.twinx()
    ax4_right.set_ylim(5400, 7870)
    ax4_right.invert_yaxis()
    ax4_right.set_ylabel('井深 / m', fontsize=11)
    
    # 整体标题
    fig.suptitle('呼101井尾管段固井质量与模型计算顶替效率沿深度分布图', 
                 fontsize=16, fontweight='bold', y=0.98)
    
    # 添加说明文字
    fig.text(0.5, 0.93, '固井质量分段依据：100312.PDF右侧固井质量评价色带的视觉抽象；效率值来自上一版1D-2D耦合+论文D2DGA环空模型结果。',
             ha='center', fontsize=10, style='italic', color='gray')
    
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"[信息] 图表已保存: {output_path}")

# --------------------------- 主程序 ---------------------------
if __name__ == '__main__':
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("[信息] 读取分段统计数据 ...")
    seg_df = load_segment_data(CSV_PATH)
    print(f"       分段数据行数: {len(seg_df)}")
    
    print("[信息] 读取模型深度剖面数据 ...")
    profile_df = load_profile_data(MODEL_PROFILE_PATH)
    print(f"       剖面数据行数: {len(profile_df)}")
    
    output_path = os.path.join(OUTPUT_DIR, '呼101_固井质量与模型效率沿深度分布图.png')
    plot_full_chart(seg_df, profile_df, output_path)
