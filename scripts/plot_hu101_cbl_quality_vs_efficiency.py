"""
绘制呼101尾管井CBL胶结质量与模型顶替效率分段对比图

本脚本读取模型深度剖面与CBL测井数据，按统一井段分段统计：
- 模型平均有效顶替效率
- CBL实测胶结质量（以 quality_proxy_pct 作为胶结质量定量代表）
- CBL胶结质量等级（胶结良好 / 胶结中等 / 差/空套风险）用于颜色标注

输出要求：
- 中文图例、标题、坐标轴
- 兼容Windows中文字体
- 文件名使用中文描述
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import sys

# --------------------------- 字体设置 ---------------------------
# 优先尝试常见Windows中文字体，并启用负号正常显示
font_list = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS', 'WenQuanYi Micro Hei']
font_found = None
for font in font_list:
    try:
        plt.rcParams['font.family'] = font
        plt.rcParams['axes.unicode_minus'] = False
        # 快速测试
        fig_test, ax_test = plt.subplots()
        ax_test.text(0.5, 0.5, '中文测试', ha='center', va='center')
        plt.close(fig_test)
        font_found = font
        break
    except Exception:
        continue

if font_found is None:
    # 若均不可用，使用默认字体但显式警告
    print("[警告] 未找到可用的中文字体，图表中文可能显示为方框。", file=sys.stderr)
else:
    print(f"[信息] 使用字体: {font_found}")

# --------------------------- 路径设置 ---------------------------
# 脚本位于 scripts/ 下，项目根目录为其父目录
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

MODEL_CSV = os.path.join(
    PROJECT_ROOT, 'results', '呼101尾管_1D2D耦合模型',
    '呼101尾管_1D2D耦合模型_深度剖面.csv'
)
CBL_CSV = os.path.join(
    PROJECT_ROOT, '参考文档', '呼101', '提取数据',
    '100312_CBL剖面_Excel版.csv'
)
OUTPUT_DIR = os.path.join(PROJECT_ROOT, 'results', '呼101尾管_1D2D耦合模型')

# --------------------------- 分段定义 ---------------------------
# 与模型结果摘要保持一致的分段口径
SEGMENTS = [
    (5402.85, 5700, '5403-5700m'),
    (5700, 6153, '5700-6153m'),
    (6153, 6800, '6153-6800m'),
    (6800, 7200, '6800-7200m'),
    (7200, 7492, '7200-7492m'),
    (7492, 7600, '7492-7600m'),
    (7600, 7735, '7600-7735m'),
    (7735, 7810, '7735-7810m'),
    (7810, 7868, '7810-7868m'),
]

# --------------------------- 数据读取 ---------------------------
def load_model_data(path):
    df = pd.read_csv(path, encoding='utf-8-sig')
    # 列名：井深_m, 平均有效顶替效率, ...
    return df

def load_cbl_data(path):
    # 第一行为 Excel 导出标记 sep=...，需跳过
    df = pd.read_csv(path, encoding='gbk', skiprows=1)
    # 过滤到有效测井深度范围
    df = df[(df['depth_md_m'] >= 5400) & (df['depth_md_m'] <= 7870)].copy()
    return df

# --------------------------- 分段统计 ---------------------------
def compute_segment_stats(model_df, cbl_df, segments):
    records = []
    for start, end, label in segments:
        # 模型效率：该深度区间内平均有效顶替效率的平均值
        m_mask = (model_df['井深_m'] >= start) & (model_df['井深_m'] < end)
        m_sub = model_df.loc[m_mask, '平均有效顶替效率']
        model_eff = float(m_sub.mean()) if len(m_sub) > 0 else np.nan

        # CBL 数据
        c_mask = (cbl_df['depth_md_m'] >= start) & (cbl_df['depth_md_m'] < end)
        c_sub = cbl_df.loc[c_mask]
        if len(c_sub) == 0:
            cbl_proxy = np.nan
            dominant_grade = '无数据'
            grade_counts = {}
        else:
            cbl_proxy = float(c_sub['quality_proxy_pct'].mean())
            grade_counts = c_sub['质量等级_代理'].value_counts().to_dict()
            dominant_grade = max(grade_counts, key=grade_counts.get)

        records.append({
            'segment_label': label,
            'depth_start': start,
            'depth_end': end,
            'model_efficiency_pct': model_eff * 100.0 if not np.isnan(model_eff) else np.nan,
            'cbl_quality_pct': cbl_proxy,
            'dominant_grade': dominant_grade,
            'grade_counts': grade_counts,
        })
    return pd.DataFrame(records)

# --------------------------- 绘图 ---------------------------
def plot_comparison(seg_df, output_path):
    """
    分组柱状图：
    - 每组两个柱子：模型效率（固定蓝色） vs CBL胶结质量（按主导等级着色）
    - 添加数值标签
    """
    # 颜色映射：CBL 主导等级 -> 颜色（匹配PDF右侧固井胶结质量道实际分布）
    # 根据100312.PDF右侧实际颜色分布：
    # - 胶结良好 = 绿色 (RGB ~0,128,0)
    # - 胶结中等 = 青绿色/蓝绿色 (RGB ~0,150,150)
    # - 胶结差 = 红色 (RGB ~255,0,0)
    # - 空套管 = 白色/空白
    grade_color_map = {
        '胶结良好': '#008000',   # 绿色 (实际分布中的胶结良好色)
        '胶结中等': '#009696',   # 青绿色/蓝绿色 (实际分布中的胶结中等色)
        '差/空套风险': '#FF0000', # 红色 (实际分布中的胶结差色)
        '无数据': '#9E9E9E',     # 灰色
    }

    x_labels = seg_df['segment_label'].tolist()
    model_vals = seg_df['model_efficiency_pct'].values
    cbl_vals = seg_df['cbl_quality_pct'].values
    cbl_colors = [grade_color_map.get(g, '#9E9E9E') for g in seg_df['dominant_grade']]

    n = len(x_labels)
    x = np.arange(n)
    width = 0.35

    fig, ax = plt.subplots(figsize=(12, 7))

    # 模型效率柱
    bars_model = ax.bar(x - width/2, model_vals, width, label='模型有效顶替效率', color='#1976D2', edgecolor='black', linewidth=0.5)
    # CBL 胶结质量柱
    bars_cbl = ax.bar(x + width/2, cbl_vals, width, label='CBL实测胶结质量', color=cbl_colors, edgecolor='black', linewidth=0.5)

    # 数值标签
    for bar in bars_model:
        height = bar.get_height()
        if not np.isnan(height):
            ax.annotate(f'{height:.1f}%',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=9)

    for bar, grade in zip(bars_cbl, seg_df['dominant_grade']):
        height = bar.get_height()
        if not np.isnan(height):
            ax.annotate(f'{height:.1f}%',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=9)

    # 对含多种等级的井段，在x轴下方添加小字注释
    for i, row in seg_df.iterrows():
        counts = row['grade_counts']
        if len(counts) > 1:
            # 构造如 "良好35/差11" 的简短注释
            parts = []
            for g, c in counts.items():
                short = g.replace('胶结', '').replace('/空套风险', '')
                parts.append(f'{short}{int(c)}')
            note = '/'.join(parts)
            ax.text(x[i], -8, note, ha='center', va='top', fontsize=7, color='#424242')

    # 轴标签与标题
    ax.set_ylabel('百分比 (%)', fontsize=12)
    ax.set_xlabel('井深分段', fontsize=12)
    ax.set_title('呼101尾管井 CBL胶结质量与模型顶替效率分段对比', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, rotation=30, ha='right', fontsize=10)
    ax.set_ylim(0, 115)
    ax.axhline(y=100, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)

    # 图例：除默认两个标签外，再补充 CBL 等级颜色说明（匹配PDF右侧实际分布色）
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#1976D2', edgecolor='black', label='模型有效顶替效率'),
        Patch(facecolor='#008000', edgecolor='black', label='CBL胶结良好（绿色）'),
        Patch(facecolor='#009696', edgecolor='black', label='CBL胶结中等（青绿色）'),
        Patch(facecolor='#FF0000', edgecolor='black', label='CBL差/空套风险（红色）'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=10)

    ax.grid(axis='y', linestyle='--', alpha=0.4)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"[信息] 图表已保存: {output_path}")

# --------------------------- 主程序 ---------------------------
if __name__ == '__main__':
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("[信息] 读取模型深度剖面 ...")
    model_df = load_model_data(MODEL_CSV)
    print(f"       模型数据行数: {len(model_df)}")

    print("[信息] 读取CBL测井数据 ...")
    cbl_df = load_cbl_data(CBL_CSV)
    print(f"       CBL数据行数: {len(cbl_df)}")

    print("[信息] 计算分段统计 ...")
    seg_df = compute_segment_stats(model_df, cbl_df, SEGMENTS)
    print(seg_df[['segment_label', 'model_efficiency_pct', 'cbl_quality_pct', 'dominant_grade']].to_string(index=False))

    output_path = os.path.join(OUTPUT_DIR, '呼101_CBL胶结质量与模型效率分段对比.png')
    plot_comparison(seg_df, output_path)
