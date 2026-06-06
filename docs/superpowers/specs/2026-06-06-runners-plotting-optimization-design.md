# Runners 绘图代码优化设计

## 概述

优化 `cemdisp/runners/` 目录中各井的绘图部分代码，重点是：
1. 删除所有 CBL 相关代码
2. 统一学术论文风格
3. 全面审查绘图逻辑，确保无数据错误、坐标轴错误、图例重叠等问题

## 设计原则

- **各井代码保持独立**：不进行代码去重，每口井的 runner 保持独立文件
- **学术论文风格**：300dpi、规范字体、清晰图例、不重叠
- **先审查后优化**：先全面检查绘图逻辑，确认正确后再调整样式

## 审查范围

### 1. reporting 模块（绘图实现）

| 文件 | 功能 | 检查重点 |
|------|------|----------|
| `reporting/plots.py` | 通用静态图 | 字段名匹配、坐标轴标签、图例位置 |
| `reporting/contour_plots.py` | 二维云图 | 色标范围、坐标轴方向、标签 |
| `reporting/animation.py` | GIF 动画 | 帧率、尺寸、中文标题 |
| `reporting/reference_figures.py` | 参考风格图 | 数据来源、图例、布局 |

### 2. runners 模块（绘图调用）

| 文件 | 对应井 | 检查重点 |
|------|--------|----------|
| `hu1_tailpipe.py` | 呼探1井 | 绘图调用参数、输出路径 |
| `hu2_tailpipe.py` | 呼探1-002井 | 同上 |
| `hu101_tailpipe.py` | 呼101 | CBL 相关代码删除、绘图调用 |
| `hu102_tailpipe.py` | 呼102 | CBL 相关代码删除 |
| `hu103_tailpipe.py` | 呼103 | 绘图调用参数 |
| `ht1_001_tailpipe.py` | 呼探1-001井 | 绘图调用参数 |
| `ht1_003_tailpipe.py` | 呼1-003井 | CBL 相关代码删除 |

### 3. data loaders（数据源）

检查各井 loader 输出的数据字段是否与绘图函数期望的列名一致。

## 检查清单

### 数据正确性
- [ ] 数据字段名称与绘图函数 `expected columns` 匹配
- [ ] 效率值范围在 0-1 之间
- [ ] 深度值单位正确（米）
- [ ] 时间值单位正确（分钟）

### 坐标轴与标签
- [ ] X/Y 轴标签正确（中文、单位）
- [ ] 坐标轴范围合理
- [ ] 刻度值格式正确

### 图例与布局
- [ ] 图例不与数据线重叠
- [ ] 图例位置合理（best/upper right/lower left 等）
- [ ] 子图之间不重叠
- [ ] `tight_layout()` 正确使用
- [ ] `bbox_inches="tight"` 保存时使用

### 样式统一
- [ ] 配色使用 `ACADEMIC_COLORS`
- [ ] 线宽：主线 2.0-2.5，辅助线 1.0-1.5
- [ ] 字号：标题 14-15pt，轴标签 11-12pt，图例 9-10pt
- [ ] 网格：alpha=0.2，虚线

## 删除 CBL 相关代码

### hu101_tailpipe.py

删除以下内容：
- `HU101_CBL_PROFILE_PATH` 常量
- `_load_hu101_cbl_profile()` 函数
- `_export_cbl_profile_comparison()` 函数
- `run_and_export()` 中调用 `_export_cbl_profile_comparison(result, output_dir)` 的行
- `import pandas as pd`（如果仅用于 CBL）
- `import matplotlib.pyplot as plt`（如果仅用于 CBL）

### hu102_tailpipe.py

删除以下内容：
- `compute_cbl_quality_proxy` 相关代码
- `validate_against_cbl` 相关代码
- `summary_payload["CBL质量风险代理预测"]` 相关代码
- `summary_payload["CBL实测对比验证"]` 相关代码
- Markdown 摘要中的 CBL 相关段落
- `from cemdisp.diagnostics.quality_proxy import compute_cbl_quality_proxy`
- `from cemdisp.validation.cbl_comparison import validate_against_cbl`
- `from cemdisp.data.validation_data import ValidationData`
- `run_and_export()` 的 `validation_data` 参数

### ht1_003_tailpipe.py

删除内容同 hu102_tailpipe.py。

### 其他 runner

检查是否有隐式 CBL 引用，如有则删除。

## 学术风格参数

```python
# 配色方案（已在 plots.py 中定义）
ACADEMIC_COLORS = {
    "primary": "#2E86AB",      # 主色调（蓝色）
    "success": "#A23B72",      # 成功/正向（紫红色）
    "warning": "#F18F01",      # 警告（橙色）
    "danger": "#C73E1D",       # 危险（红色）
    "info": "#3B1F2B",         # 信息（深色）
    "light": "#44BBA4",        # 浅色（绿色）
}

# 图表参数
DPI = 300
FONT_SIZE_TITLE = 14
FONT_SIZE_LABEL = 11
FONT_SIZE_LEGEND = 9
LINE_WIDTH_MAIN = 2.2
LINE_WIDTH_AUX = 1.5
GRID_ALPHA = 0.2
GRID_STYLE = "--"
```

## 实施计划

### 阶段 1：全面审查（并行）

使用子 agent 并行检查：
1. reporting 模块的绘图实现
2. 各井 runner 的绘图调用
3. data loaders 的数据输出

### 阶段 2：删除 CBL 代码（并行）

使用子 agent 并行修改：
1. hu101_tailpipe.py
2. hu102_tailpipe.py
3. ht1_003_tailpipe.py

### 阶段 3：样式统一

修改 `reporting/plots.py`、`reporting/contour_plots.py` 中的样式参数。

### 阶段 4：验证

使用子 agent 逐井检查修改后的代码，确保：
- 无语法错误
- 绘图逻辑正确
- 图例不重叠
- 输出文件路径正确

## 输出文件

- 设计文档：`docs/superpowers/specs/2026-06-06-runners-plotting-optimization-design.md`
- 修改的文件：各 runner 文件、reporting 模块文件
