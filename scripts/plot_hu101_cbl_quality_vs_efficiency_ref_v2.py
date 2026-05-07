"""
绘制呼101尾管井CBL固井质量与模型顶替效率分段对比图

数据来源：参考模型分段统计表
包含：宽边效率、模型平均效率、窄边效率 + CBL固井质量等级
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
    - 每组三个柱子：宽边效率、模型平均效率、窄边效率
    - 用背景色/柱体边框色表示CBL固井质量等级
    """
    # 颜色映射：CBL 固井质量解释 -> 颜色（匹配PDF右侧实际分布）
    quality_color_map = {
        '良好': '#008000',       # 绿色
        '较好': '#66BB6A',       # 浅绿色
        '中等': '#009696',       # 青绿色
        '一般': '#FFB300',       # 琥珀/黄色
        '局部差': '#FF5722',     # 橙红色
        '差': '#FF0000',         # 红色
        '未测井评价': '#9E9E9E', # 灰色
    }

    x_labels = df['井段_m'].tolist()
    wide_vals = df['宽边效率_%'].values
    model_vals = df['模型顶替效率_%'].values
    narrow_vals = df['窄边效率_%'].values
    cbl_quality = df['固井质量解释'].tolist()
    
    n = len(x_labels)
    x = np.arange(n)
    width = 0.25

    fig, ax = plt.subplots(figsize=(16, 9))

    # 绘制三组柱子
    bars_wide = ax.bar(x - width, wide_vals, width, label='宽边效率', color='#4FC3F7', edgecolor='black', linewidth=0.5)
    bars_model = ax.bar(x, model_vals, width, label='模型平均效率', color='#1976D2', edgecolor='black', linewidth=0.8)
    bars_narrow = ax.bar(x + width, narrow_vals, width, label='窄边效率', color='#90A4AE', edgecolor='black', linewidth=0.5)

    # 为每个组添加CBL质量等级的背景色带
    for i, quality in enumerate(cbl_quality):
        color = quality_color_map.get(quality, '#9E9E9E')
        # 在x位置添加半透明背景矩形
        ax.axvspan(i - 0.5, i + 0.5, alpha=0.15, color=color, zorder=0)

    # 数值标签
    for bar in bars_wide:
        height = bar.get_height()
        if not np.isnan(height):
            ax.annotate(f'{height:.1f}',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=8, color='#0277BD')

    for bar in bars_model:
        height = bar.get_height()
        if not np.isnan(height):
            ax.annotate(f'{height:.1f}',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=9, fontweight='bold', color='#0D47A1')

    for bar in bars_narrow:
        height = bar.get_height()
        if not np.isnan(height):
            ax.annotate(f'{height:.1f}',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=8, color='#455A64')

    # 在x轴下方添加CBL质量等级标注
    for i, quality in enumerate(cbl_quality):
        color = quality_color_map.get(quality, '#9E9E9E')
        ax.text(i, -8, quality, ha='center', va='top', fontsize=10, fontweight='bold', color=color)

    # 轴标签与标题
    ax.set_ylabel('顶替效率 (%)', fontsize=12)
    ax.set_xlabel('井深分段', fontsize=12)
    ax.set_title('呼101尾管井 各段固井质量与模型顶替效率对比\n（宽边/平均/窄边效率 + CBL固井质量等级）', 
                 fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, rotation=30, ha='right', fontsize=10)
    ax.set_ylim(-15, 110)
    ax.axhline(y=100, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)

    # 图例
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#4FC3F7', edgecolor='black', label='宽边效率'),
        Patch(facecolor='#1976D2', edgecolor='black', label='模型平均效率'),
        Patch(facecolor='#90A4AE', edgecolor='black', label='窄边效率'),
        Patch(facecolor='#008000', edgecolor='black', label='CBL良好'),
        Patch(facecolor='#66BB6A', edgecolor='black', label='CBL较好'),
        Patch(facecolor='#009696', edgecolor='black', label='CBL中等'),
        Patch(facecolor='#FFB300', edgecolor='black', label='CBL一般'),
        Patch(facecolor='#FF5722', edgecolor='black', label='CBL局部差'),
        Patch(facecolor='#FF0000', edgecolor='black', label='CBL差'),
        Patch(facecolor='#9E9E9E', edgecolor='black', label='未测井评价'),
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
    print(df[['井段_m', '固井质量解释', '模型顶替效率_%', '宽边效率_%', '窄边效率_%']].to_string(index=False))

    output_path = os.path.join(OUTPUT_DIR, '呼101_CBL固井质量与模型效率分段对比_参考模型_v2.png')
    plot_comparison(df, output_path)
