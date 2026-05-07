"""
绘制呼101尾管井CBL固井质量与模型顶替效率分段对比图

数据来源：参考模型分段统计表
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
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
# 脚本位于 scripts/ 下，项目根目录为其父目录
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CSV_PATH = os.path.join(
    PROJECT_ROOT, '参考文档', 'hu101_1d2d_paper_project', 
    '呼101尾管_固井质量_分段顶替效率图表',
    '呼101尾管_各段固井质量与模型顶替效率.csv'
)

OUTPUT_DIR = os.path.join(PROJECT_ROOT, 'results', '呼101尾管_1D2D耦合模型')

# --------------------------- 数据读取 ---------------------------
def load_data(path):
    df = pd.read_csv(path, encoding='utf-8-sig')
    return df

# --------------------------- 绘图 ---------------------------
def plot_comparison(df, output_path):
    """
    分组柱状图：
    - 每组两个柱子：模型效率（蓝色） vs CBL固井质量等级（按等级着色）
    """
    # 颜色映射：CBL 固井质量解释 -> 颜色（匹配PDF右侧实际分布）
    quality_color_map = {
        '良好': '#008000',       # 绿色
        '较好': '#66BB6A',       # 浅绿色（介于良好和中等之间）
        '中等': '#009696',       # 青绿色
        '一般': '#FFB300',       # 琥珀/黄色（一般质量）
        '局部差': '#FF5722',     # 橙红色（局部差）
        '差': '#FF0000',         # 红色（差）
        '未测井评价': '#9E9E9E', # 灰色
    }

    x_labels = df['井段_m'].tolist()
    model_vals = df['模型顶替效率_%'].values
    
    # CBL侧：我们使用100%作为基准，但用颜色表示质量等级
    # 这样可以在图上同时显示模型效率和质量等级颜色
    cbl_vals = np.ones(len(df)) * 100.0  # 用100%作为高度，颜色表示质量
    cbl_colors = [quality_color_map.get(q, '#9E9E9E') for q in df['固井质量解释']]

    n = len(x_labels)
    x = np.arange(n)
    width = 0.35

    fig, ax = plt.subplots(figsize=(14, 8))

    # 模型效率柱
    bars_model = ax.bar(x - width/2, model_vals, width, label='模型顶替效率', color='#1976D2', edgecolor='black', linewidth=0.5)
    
    # CBL 固井质量柱（高度100%，颜色表示质量等级）
    bars_cbl = ax.bar(x + width/2, cbl_vals, width, label='CBL固井质量等级', color=cbl_colors, edgecolor='black', linewidth=0.5, alpha=0.8)

    # 数值标签 - 模型效率
    for bar in bars_model:
        height = bar.get_height()
        if not np.isnan(height):
            ax.annotate(f'{height:.1f}%',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=9, fontweight='bold')

    # CBL柱上标注质量等级文字
    for bar, quality in zip(bars_cbl, df['固井质量解释']):
        height = bar.get_height()
        ax.annotate(f'{quality}',
                    xy=(bar.get_x() + bar.get_width() / 2, height/2),
                    ha='center', va='center', fontsize=9, fontweight='bold',
                    color='white' if quality in ['差', '局部差', '良好', '较好'] else 'black')

    # 添加宽窄边效率差异的次坐标轴或注释
    for i, row in df.iterrows():
        diff = row['宽窄边差值_%']
        if diff > 20:  # 高差异标注
            ax.text(x[i], model_vals[i] + 8, f'宽窄差:{diff:.1f}%', 
                   ha='center', va='bottom', fontsize=7, color='#D32F2F', style='italic')

    # 轴标签与标题
    ax.set_ylabel('百分比 (%)', fontsize=12)
    ax.set_xlabel('井深分段', fontsize=12)
    ax.set_title('呼101尾管井 CBL固井质量与模型顶替效率分段对比', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, rotation=30, ha='right', fontsize=10)
    ax.set_ylim(0, 120)
    ax.axhline(y=100, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)

    # 图例
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#1976D2', edgecolor='black', label='模型顶替效率'),
        Patch(facecolor='#008000', edgecolor='black', label='CBL良好（绿色）'),
        Patch(facecolor='#66BB6A', edgecolor='black', label='CBL较好（浅绿）'),
        Patch(facecolor='#009696', edgecolor='black', label='CBL中等（青绿）'),
        Patch(facecolor='#FFB300', edgecolor='black', label='CBL一般（黄色）'),
        Patch(facecolor='#FF5722', edgecolor='black', label='CBL局部差（橙红）'),
        Patch(facecolor='#FF0000', edgecolor='black', label='CBL差（红色）'),
        Patch(facecolor='#9E9E9E', edgecolor='black', label='未测井评价（灰色）'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=9, ncol=2)

    ax.grid(axis='y', linestyle='--', alpha=0.4)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"[信息] 图表已保存: {output_path}")

# --------------------------- 主程序 ---------------------------
if __name__ == '__main__':
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("[信息] 读取参考模型分段数据 ...")
    df = load_data(CSV_PATH)
    print(f"       数据行数: {len(df)}")
    print(df[['井段_m', '固井质量解释', '模型顶替效率_%']].to_string(index=False))

    output_path = os.path.join(OUTPUT_DIR, '呼101_CBL固井质量与模型效率分段对比_参考模型.png')
    plot_comparison(df, output_path)
