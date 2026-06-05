# 固井模型图表输出优化设计文档

> **日期**：2026-06-05
> **版本**：v1.1
> **状态**：已确认（更新：删除CBL/目标层段相关，聚焦全井分布）

---

## 1. 项目背景

### 1.1 当前图表模块状态

**reporting/plots.py** - 已实现功能：
- 时间序列图（plot_time_series）：有效顶替效率与前沿推进
- 深度剖面图（plot_depth_profiles）：深度方向效率与浓度分布
- 风险指标图（plot_risk_indices）：窜槽、混浆、失稳风险演变
- 结果对比图（plot_efficiency_summary_bar）：最终各项指标对比

**reporting/contour_plots.py** - 已实现功能：
- 深度-时间顶替效率云图（plot_depth_time_contour）
- 环空水泥与隔离液浓度场多时刻快照图（plot_annulus_snapshots）
- 最终场分布图（plot_final_fields_contour）

**reporting/animation.py** - 已实现功能：
- 浓度场演化动画（animate_cement_field）

**reporting/reference_figures.py** - 已实现功能：
- 参考项目风格图件输出

### 1.2 已发现的问题

#### 1.2.1 数据计算问题

| 问题 | 位置 | 说明 |
|------|------|------|
| 旧版本兼容逻辑 | contour_plots.py | `_require_snapshot_data`、`_optional_spacer_snapshots`、`_optional_spacer_field` 中有旧版本兼容代码 |
| 无效计算 | contour_plots.py | `effective_field = cement_field * (1.0 - wall_field)` 中 wall_field 当前恒为零 |
| 旧版本兼容 | plots.py | `plot_efficiency_summary_bar` 中有 `cbl_quality_proxy` 兼容逻辑 |
| 冗余字段 | plots.py, contour_plots.py | CBL评价段、目标层段相关字段和图表需要删除 |

#### 1.2.2 样式问题

| 问题 | 位置 | 说明 |
|------|------|------|
| 等值线标注重叠 | contour_plots.py | 时间 125-150 min 之间，等值线标注拥挤 |
| 图例遮挡数据 | plots.py | 图例遮挡了部分起始数据点 |
| 曲线震荡 | plots.py | 井深 7500m 之后曲线剧烈震荡，需要平滑处理 |
| X轴范围 | contour_plots.py | X轴刻度止于 175，但数据延伸到约 187.5 |

#### 1.2.3 布局问题

| 问题 | 位置 | 说明 |
|------|------|------|
| 多面板排列 | contour_plots.py | 多面板图表的排列可以优化 |
| colorbar | contour_plots.py | colorbar 的位置和标签可以调整 |
| 图表尺寸 | plots.py, contour_plots.py | 图表尺寸可以统一 |

---

## 2. 优化目标

1. **中文显示优化**：统一设置全局图表样式
2. **数据计算修复**：删除旧版本兼容逻辑，修复无效计算，删除CBL/目标层段相关
3. **样式调整**：修复等值线标注、图例位置、曲线平滑等问题
4. **布局优化**：优化多面板图表的排列，调整 colorbar 位置
5. **新增图表类型**：湍流区域分布图、雷诺数云图、流态判断图
6. **聚焦全井分布**：所有图表只展示全井段数据，不区分CBL评价段和目标层段

---

## 3. 设计方案

### 3.1 全局样式设置

#### 3.1.1 新增 `setup_academic_style()` 函数

在 `plots.py` 中新增全局学术风格设置函数：

```python
def setup_academic_style() -> None:
    """设置学术论文风格的全局图表样式。

    包括：
    - 中文字体设置（SimHei、Microsoft YaHei 等）
    - LaTeX 数学公式支持
    - 全局网格样式
    - 坐标轴样式
    - 线宽、标记大小等默认值
    """
    # 中文字体设置
    _setup_chinese_font()

    # LaTeX 支持
    plt.rcParams["text.usetex"] = False  # 不使用系统 LaTeX，使用 matplotlib 内置
    plt.rcParams["mathtext.fontset"] = "cm"  # 使用 Computer Modern 字体

    # 全局网格样式
    plt.rcParams["axes.grid"] = True
    plt.rcParams["grid.alpha"] = 0.25
    plt.rcParams["grid.linestyle"] = "--"
    plt.rcParams["grid.linewidth"] = 0.5

    # 坐标轴样式
    plt.rcParams["axes.linewidth"] = 0.8
    plt.rcParams["axes.edgecolor"] = "#333333"
    plt.rcParams["axes.labelsize"] = 12
    plt.rcParams["axes.titlesize"] = 14
    plt.rcParams["axes.titleweight"] = "bold"

    # 刻度样式
    plt.rcParams["xtick.labelsize"] = 10
    plt.rcParams["ytick.labelsize"] = 10
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"

    # 图例样式
    plt.rcParams["legend.fontsize"] = 10
    plt.rcParams["legend.framealpha"] = 0.8
    plt.rcParams["legend.edgecolor"] = "#333333"

    # 线宽和标记大小
    plt.rcParams["lines.linewidth"] = 1.8
    plt.rcParams["lines.markersize"] = 6

    # 图表尺寸
    plt.rcParams["figure.figsize"] = (10, 6)
    plt.rcParams["figure.dpi"] = 150
```

#### 3.1.2 学术颜色方案

定义学术论文风格的颜色方案：

```python
# 学术颜色方案
ACADEMIC_COLORS = {
    "primary": "#2171b5",      # 主色调（蓝色）
    "secondary": "#6baed6",    # 次色调（浅蓝色）
    "accent": "#cb181d",       # 强调色（红色）
    "success": "#238b45",      # 成功色（绿色）
    "warning": "#fe9929",      # 警告色（橙色）
    "neutral": "#636363",      # 中性色（灰色）
    "light": "#f0f0f0",        # 浅色背景
    "dark": "#252525",         # 深色文字
}

# 多线颜色方案（用于多条曲线）
MULTI_LINE_COLORS = ["#2171b5", "#cb181d", "#238b45", "#fe9929", "#6a3d9a", "#ff7f00"]

# 云图颜色方案
CONTOUR_CMAP = {
    "efficiency": "RdYlGn",      # 顶替效率（红-黄-绿）
    "cement": "viridis",          # 水泥浓度
    "spacer": "coolwarm",         # 隔离液浓度
    "reynolds": "plasma",         # 雷诺数
    "turbulence": "YlOrRd",       # 湍流强度
}
```

### 3.2 数据计算修复

#### 3.2.1 删除旧版本兼容逻辑

**contour_plots.py 中要删除的代码**：

| 函数 | 删除内容 |
|------|----------|
| `_require_snapshot_data` | 旧的 `snapshot_times_s` 兼容逻辑 |
| `_optional_spacer_snapshots` | 旧版本单流体兼容逻辑 |
| `_optional_spacer_field` | 旧版本兼容逻辑 |

**plots.py 中要删除的代码**：

| 函数 | 删除内容 |
|------|----------|
| `plot_efficiency_summary_bar` | `cbl_quality_proxy` 兼容逻辑 |

#### 3.2.2 删除 CBL 和目标层段相关

**删除以下内容**：

1. **plots.py 中的 `plot_time_series`**：
   - 删除 CBL 评价段效率曲线
   - 删除目标层段效率曲线
   - 只保留全井段有效顶替效率

2. **plots.py 中的 `plot_efficiency_summary_bar`**：
   - 删除 CBL 评价井段有效顶替效率柱状
   - 删除目标层段有效顶替效率柱状
   - 只保留全井段有效顶替效率和水泥浆占据率

3. **AnnulusD2DGASolver.run() 中的指标计算**：
   - 删除 `target_efficiency` 计算
   - 删除 `cbl_efficiency` 计算
   - 只保留 `effective_efficiency`

4. **AnnulusSimulationResult 中的字段**：
   - 删除 `target_interval_efficiency` 相关
   - 删除 `cbl_eval_interval_efficiency` 相关

#### 3.2.3 修复无效计算

**contour_plots.py 中的 `plot_final_fields_contour`**：

```python
# 修改前（无效计算）
effective_field = cement_field * (1.0 - wall_field)

# 修改后（直接使用水泥浓度）
effective_field = cement_field  # wall_field 当前恒为零，直接使用水泥浓度作为有效顶替效率
```

### 3.3 样式调整

#### 3.3.1 修复等值线标注重叠

**contour_plots.py 中的 `plot_depth_time_contour`**：

```python
# 修改前
contour_lines = ax.contour(
    time_grid,
    depth_grid,
    depth_time_matrix,
    levels=[0.5, 0.7, 0.9],
    colors="#333333",
    linewidths=0.8,
)
ax.clabel(contour_lines, fmt="%.1f", inline=True, fontsize=9)

# 修改后（减少标注密度，使用自动避让）
contour_lines = ax.contour(
    time_grid,
    depth_grid,
    depth_time_matrix,
    levels=[0.5, 0.7, 0.9],
    colors="#333333",
    linewidths=0.8,
)
ax.clabel(contour_lines, fmt="%.1f", inline=True, fontsize=9, manual=False)
```

#### 3.3.2 优化图例位置

**plots.py 中的 `plot_depth_profiles`**：

```python
# 修改前
axes[0].legend(loc="best", ncol=2)

# 修改后（将图例放在绘图区域之外）
axes[0].legend(
    loc="upper left",
    bbox_to_anchor=(1.02, 1.0),
    fontsize=10,
    framealpha=0.8,
    edgecolor="#333333",
    ncol=1,
)
```

#### 3.3.3 曲线平滑处理

**plots.py 中的 `plot_depth_profiles`**：

```python
# 添加平滑函数
from scipy.signal import savgol_filter

def _smooth_curve(data: np.ndarray, window_length: int = 11, polyorder: int = 3) -> np.ndarray:
    """使用 Savitzky-Golay 滤波器平滑曲线。

    Args:
        data: 原始数据
        window_length: 窗口长度（必须为奇数）
        polyorder: 多项式阶数

    Returns:
        平滑后的数据
    """
    if len(data) < window_length:
        return data
    return savgol_filter(data, window_length, polyorder)

# 使用示例
axes[0].plot(
    depth,
    _smooth_curve(profiles["宽边有效效率"]),
    label="宽边有效效率",
    linewidth=2,
)
```

#### 3.3.4 调整 X 轴范围

**contour_plots.py 中的 `plot_depth_time_contour`**：

```python
# 修改前
ax.set_xlabel("时间 / min")

# 修改后（设置合适的 X 轴范围）
ax.set_xlabel("时间 / min")
ax.set_xlim(left=0)  # 从 0 开始
```

#### 3.3.5 统一坐标轴标签格式

```python
# 坐标轴标签格式
axes[0].set_ylabel("有效顶替效率 / 占据率", fontsize=12)
axes[0].set_xlabel("时间 / min", fontsize=12)
```

#### 3.3.6 优化图例样式

```python
# 图例样式
axes[0].legend(
    loc="best",
    fontsize=10,
    framealpha=0.8,
    edgecolor="#333333",
    ncol=1,
)
```

#### 3.3.7 调整线宽和标记

```python
# 线宽和标记
axes[0].plot(
    metrics["time_min"],
    metrics["effective_efficiency"],
    label="全井段有效顶替效率",
    linewidth=2.0,
    color=ACADEMIC_COLORS["primary"],
    marker="o",
    markersize=4,
    markevery=10,  # 每10个点标记一次
)
```

### 3.4 布局优化

#### 3.4.1 统一图表尺寸

```python
# 标准图表尺寸
STANDARD_FIGSIZE = {
    "single": (10, 6),       # 单面板
    "double": (10, 8),       # 双面板
    "triple": (12, 10),      # 三面板
    "quad": (12, 8),         # 四面板
    "wide": (14, 6),         # 宽幅图
}
```

#### 3.4.2 优化 colorbar 位置

```python
# colorbar 位置优化
colorbar = fig.colorbar(
    filled,
    ax=ax,
    pad=0.02,
    shrink=0.8,
    aspect=30,
)
colorbar.set_label("有效顶替效率", fontsize=11)
colorbar.ax.tick_params(labelsize=9)
```

### 3.5 新增图表类型

#### 3.5.1 湍流区域分布图

基于雷诺数绘制湍流区域分布图：

```python
def plot_turbulence_region(
    result: AnnulusSimulationResult,
    output_dir: Path | str | None = None,
) -> Figure:
    """绘制湍流区域分布图。

    基于雷诺数判断流态，标记层流和湍流区域：
    - 层流区域（Re < Re_critical）：蓝色
    - 湍流区域（Re > Re_critical）：红色

    Args:
        result: 环空二维模拟结果，需包含 reynolds_snapshots
        output_dir: 可选输出目录

    Returns:
        matplotlib Figure 对象
    """
    # 获取最终时刻的雷诺数场
    reynolds_field = result.reynolds_snapshots[-1] if result.reynolds_snapshots else None
    if reynolds_field is None:
        raise ValueError("结果中不包含雷诺数快照，无法绘制湍流区域分布图。")

    # 判断流态
    re_critical = 2100.0  # 临界雷诺数
    turbulence_mask = reynolds_field > re_critical

    # 绘图
    fig, ax = plt.subplots(figsize=(10, 6))
    depth = result.geom["md"]

    # 绘制湍流区域
    ax.imshow(
        turbulence_mask.astype(float),
        extent=_field_extent(depth),
        origin="lower",
        aspect="auto",
        cmap="RdBu_r",
        vmin=0.0,
        vmax=1.0,
        alpha=0.6,
    )

    ax.set_title(f"{result.well_name} 湍流区域分布", fontsize=14, fontweight="bold")
    ax.set_xlabel("井深 / m")
    ax.set_ylabel("方位角 (宽边→窄边)")

    # 添加 colorbar
    colorbar = fig.colorbar(ax.images[0], ax=ax, pad=0.02)
    colorbar.set_label("流态 (0=层流, 1=湍流)")
    colorbar.set_ticks([0.25, 0.75])
    colorbar.set_ticklabels(["层流", "湍流"])

    fig.tight_layout()
    _save_figure(fig, output_dir, f"{well_name}_湍流区域分布.png")
    return fig
```

#### 3.5.2 雷诺数云图

```python
def plot_reynolds_contour(
    result: AnnulusSimulationResult,
    output_dir: Path | str | None = None,
) -> Figure:
    """绘制雷诺数云图。

    展示最终时刻的雷诺数分布，用于分析流态：
    - Re < 2100：层流区域
    - 2100 < Re < 4000：过渡区域
    - Re > 4000：湍流区域

    Args:
        result: 环空二维模拟结果，需包含 reynolds_snapshots
        output_dir: 可选输出目录

    Returns:
        matplotlib Figure 对象
    """
    # 获取最终时刻的雷诺数场
    reynolds_field = result.reynolds_snapshots[-1] if result.reynolds_snapshots else None
    if reynolds_field is None:
        raise ValueError("结果中不包含雷诺数快照，无法绘制雷诺数云图。")

    # 绘图
    fig, ax = plt.subplots(figsize=(10, 6))
    depth = result.geom["md"]

    # 绘制雷诺数云图
    image = ax.imshow(
        reynolds_field,
        extent=_field_extent(depth),
        origin="lower",
        aspect="auto",
        cmap="plasma",
        norm=plt.cm.colors.LogNorm(vmin=100, vmax=10000),
    )

    ax.set_title(f"{result.well_name} 雷诺数分布", fontsize=14, fontweight="bold")
    ax.set_xlabel("井深 / m")
    ax.set_ylabel("方位角 (宽边→窄边)")

    # 添加 colorbar
    colorbar = fig.colorbar(image, ax=ax, pad=0.02)
    colorbar.set_label("雷诺数 Re")

    # 添加临界雷诺数参考线
    ax.axhline(y=2100, color="white", linestyle="--", linewidth=1.5, label="Re_critical = 2100")
    ax.legend(loc="upper right")

    fig.tight_layout()
    _save_figure(fig, output_dir, f"{well_name}_雷诺数分布.png")
    return fig
```

#### 3.5.3 流态判断图

```python
def plot_flow_regime(
    result: AnnulusSimulationResult,
    output_dir: Path | str | None = None,
) -> Figure:
    """绘制流态判断图。

    展示不同井深处的流态分布：
    - 宽边、中线、窄边的雷诺数随深度变化
    - 临界雷诺数参考线
    - 层流/湍流区域标记

    Args:
        result: 环空二维模拟结果，需包含 reynolds_snapshots
        output_dir: 可选输出目录

    Returns:
        matplotlib Figure 对象
    """
    # 获取最终时刻的雷诺数场
    reynolds_field = result.reynolds_snapshots[-1] if result.reynolds_snapshots else None
    if reynolds_field is None:
        raise ValueError("结果中不包含雷诺数快照，无法绘制流态判断图。")

    depth = result.geom["md"]
    ny, nz = reynolds_field.shape

    # 提取宽边、中线、窄边的雷诺数
    re_wide = reynolds_field[0, :]           # 宽边
    re_mid = reynolds_field[ny // 2, :]      # 中线
    re_narrow = reynolds_field[-1, :]        # 窄边

    # 绘图
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.semilogy(depth, re_wide, label="宽边", color=ACADEMIC_COLORS["primary"], linewidth=2)
    ax.semilogy(depth, re_mid, label="中线", color=ACADEMIC_COLORS["neutral"], linewidth=2)
    ax.semilogy(depth, re_narrow, label="窄边", color=ACADEMIC_COLORS["accent"], linewidth=2)

    # 添加临界雷诺数参考线
    ax.axhline(y=2100, color="black", linestyle="--", linewidth=1.5, label=r"$Re_{critical}$ = 2100")

    # 填充层流和湍流区域
    ax.fill_between(depth, 100, 2100, alpha=0.1, color="blue", label="层流区域")
    ax.fill_between(depth, 2100, 10000, alpha=0.1, color="red", label="湍流区域")

    ax.set_title(f"{result.well_name} 流态判断图", fontsize=14, fontweight="bold")
    ax.set_xlabel("井深 / m", fontsize=12)
    ax.set_ylabel("雷诺数 Re", fontsize=12)
    ax.set_ylim(100, 10000)
    ax.legend(loc="best", fontsize=10)
    ax.grid(True, alpha=0.25)

    fig.tight_layout()
    _save_figure(fig, output_dir, f"{well_name}_流态判断图.png")
    return fig
```

---

## 4. 实施顺序

### 第一阶段：数据清理

1. 删除 contour_plots.py 中的旧版本兼容逻辑
2. 删除 plots.py 中的旧版本兼容逻辑
3. 删除 CBL 评价段和目标层段相关代码
4. 修复 `effective_field` 计算
5. 更新 AnnulusD2DGASolver.run() 中的指标计算
6. 更新 AnnulusSimulationResult 字段

### 第二阶段：样式修复

1. 修复等值线标注重叠问题
2. 优化图例位置，避免遮挡数据
3. 添加曲线平滑处理
4. 调整 X 轴范围

### 第三阶段：全局样式设置

1. 新增 `setup_academic_style()` 函数
2. 定义学术颜色方案
3. 统一坐标轴标签格式
4. 优化图例样式
5. 调整线宽和标记

### 第四阶段：布局优化

1. 统一图表尺寸
2. 优化 colorbar 位置
3. 调整多面板排列

### 第五阶段：新增图表类型

1. 实现湍流区域分布图
2. 实现雷诺数云图
3. 实现流态判断图

---

## 5. 验证标准

### 5.1 功能验证

- [ ] 所有图表能正常生成
- [ ] 中文标签、图例、标题显示正确
- [ ] 数据计算逻辑正确
- [ ] 新增图表类型功能正常
- [ ] CBL 和目标层段相关代码已完全删除
- [ ] 只展示全井段数据

### 5.2 样式验证

- [ ] 等值线标注不重叠
- [ ] 图例不遮挡数据
- [ ] 曲线平滑处理正确
- [ ] X 轴范围合适
- [ ] 颜色方案符合学术论文风格
- [ ] 坐标轴标签格式统一

### 5.3 布局验证

- [ ] 图表尺寸统一
- [ ] colorbar 位置合理
- [ ] 多面板排列整齐

---

## 6. 参考文献

1. Matplotlib 官方文档：https://matplotlib.org/stable/
2. 学术论文图表规范：Nature、Science 等期刊的图表要求
3. 中文字体设置：SimHei、Microsoft YaHei 等字体的使用规范
