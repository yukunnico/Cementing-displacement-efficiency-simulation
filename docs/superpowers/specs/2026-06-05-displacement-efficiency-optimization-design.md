# 固井顶替效率模型优化设计文档

> **日期**：2026-06-05
> **版本**：v1.0
> **状态**：已确认

---

## 1. 项目背景

### 1.1 当前模型状态

**套管段模型 (CasingFlowSolver)** - 已实现功能：
- 体积推进法追踪多流体前缘
- 重力沉降修正（密度差驱动，支持井斜角投影）
- 屈服应力修正（抑制沉降）
- 停泵沉降增强模型（凝胶强度发展 + 屈服应力效应）
- 轴向弥散（Taylor-Aris 型，将尖锐前缘转为平滑 S 形过渡带）
- 鞋口出流事件时间轴
- 双径向井支持

**环空段模型 (AnnulusD2DGASolver)** - 已实现功能：
- 偏心环空几何展开
- 基于局部流动度的轴向/方位角平均速度场
- D2DGA 通量放大修正（Zhang & Frigaard 2022）
- 半拉格朗日平流输运（双线性插值）
- 显式扩散（轴向 + 方位角 Laplacian）
- Papanastasiou 正则化屈服应力模型
- 四相跟踪（领浆、尾浆、前置液/隔离液、钻井液）
- 泵停冻结
- 开放边界（允许水泥浆流出到重叠段）

### 1.2 已禁用功能（代码中标注为"论文口径核心不再使用"）

- 泥饼清除模型
- 温度修正
- 凝胶强度影响
- 湍流修正
- CBL 质量惩罚

### 1.3 现场数据情况

- **井身结构数据**：有，但提取不全面
- **流体物性数据**：有，但提取不全面
- **施工过程数据**：有，但提取不全面
- **CBL 验证数据**：暂时不在模型中添加，后续再提取

---

## 2. 优化目标

1. **添加湍流修正**：基于雷诺数自动判断流态，湍流条件下修正粘度
2. **清理 legacy 兼容代码**：删除所有不再使用的旧字段、旧函数、兼容参数
3. **优化数据文件**：提取公共逻辑，各井数据分离，统一接口格式
4. **优化代码注释**：整个模型代码的注释都要优化，删掉老版本注释

---

## 3. 设计方案

### 3.1 湍流修正模型

#### 3.1.1 物理背景

参考文献：
- Maleki & Frigaard (2019) - "Comparing laminar and turbulent primary cementing flows"
- Zhang & Frigaard (2022) - JFM Vol.947, A32

关键发现：
- 湍流可以显著提高偏心环空窄边的顶替效率
- 层流条件下偏心度对顶替效率影响更大
- 流态判断基于广义雷诺数

#### 3.1.2 实现方案

**新增参数**（`AnnulusD2DGASolver.__init__`）：

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `enable_turbulence` | `bool` | `True` | 是否启用湍流修正 |
| `re_critical` | `float` | `2100.0` | 临界雷诺数 |
| `turbulence_enhancement_factor` | `float` | `0.2` | 湍流增强系数 α |

**新增方法**：

```python
def _compute_reynolds_number(
    self,
    rho: Array,
    w: Array,
    b: Array,
    mu: Array,
) -> Array:
    """计算局部广义雷诺数。

    公式：Re_gen = ρ * V * D_h / μ_apparent
    其中：
    - ρ：混合流体密度 (kg/m³)
    - V：轴向速度 (m/s)
    - D_h：水力直径 = 2×局部间隙宽度b (m)
    - μ_apparent：表观粘度 (Pa·s)

    Args:
        rho: 密度场 (ny×nz)，单位 kg/m³
        w: 轴向速度场 (ny×nz)，单位 m/s
        b: 局部间隙宽度 (ny×nz)，单位 m
        mu: 表观粘度场 (ny×nz)，单位 Pa·s

    Returns:
        雷诺数场 (ny×nz)，无量纲
    """
    D_h = 2.0 * b
    return rho * np.abs(w) * D_h / np.maximum(mu, 1.0e-6)


def _turbulent_viscosity_correction(
    self,
    mu: Array,
    Re: Array,
) -> Array:
    """湍流粘度修正。

    当 Re > Re_critical 时，使用湍流增强因子修正粘度：
    μ_effective = μ_laminar × (1 + α × max(0, Re/Re_critical - 1))

    物理意义：
    - 层流区（Re < Re_critical）：粘度不变
    - 湍流区（Re > Re_critical）：粘度增大，反映湍流脉动对动量传递的增强

    Args:
        mu: 层流表观粘度场 (ny×nz)，单位 Pa·s
        Re: 雷诺数场 (ny×nz)，无量纲

    Returns:
        修正后的有效粘度场 (ny×nz)，单位 Pa·s
    """
    excess_ratio = np.maximum(Re / self.re_critical - 1.0, 0.0)
    enhancement = 1.0 + self.turbulence_enhancement_factor * excess_ratio
    return mu * enhancement
```

**修改 `_compute_velocity` 方法**：

在计算完层流粘度后，添加湍流修正：

```python
# 计算雷诺数
Re = self._compute_reynolds_number(rho_kg_m3, w_prev, b, mu_reg)

# 湍流粘度修正
if self.enable_turbulence:
    mu_reg = self._turbulent_viscosity_correction(mu_reg, Re)
```

**诊断输出**：

保留 `reynolds_snapshots` 字段，用于输出雷诺数场快照，便于分析流态分布。

#### 3.1.3 临界雷诺数说明

- 牛顿流体：Re_critical = 2100（经典值）
- 非牛顿流体：使用 Metzner-Reed 广义雷诺数，临界值约 2100
- 实际应用中可根据现场数据调整

---

### 3.2 Legacy 兼容代码清理

#### 3.2.1 删除范围

**AnnulusSimulationResult 中要删除的字段**：

| 字段名 | 说明 | 删除原因 |
|--------|------|----------|
| `gel_strength_snapshots` | 凝胶强度快照 | 不再填充数据 |
| `mud_cake_field` | 泥饼场 | 不再填充数据 |
| `mud_cake_snapshots` | 泥饼快照 | 不再填充数据 |
| `turbulent_viscosity_snapshots` | 湍流粘度快照 | 不再填充数据 |
| `snapshot_times_s` | 快照时间 | 与 `time_points_s` 重复 |

**AnnulusD2DGASolver 中要删除/简化的代码**：

| 内容 | 说明 |
|------|------|
| `_update_effective_gap` 方法 | 当前实现只是赋值，可删除并在调用处使用 `geom["b"]` |
| `_build_geom` 中的 `mud_cake_thickness` 参数 | 已不再使用 |
| `_apparent_viscosity` 中的 `gel_strength` 和 `temperature_correction` 参数 | 已不再使用 |
| `_compute_props` 中的 `gel_strength` 和 `temperature_correction` 参数 | 已不再使用 |
| `_compute_velocity` 中的 `gel_strength` 和 `temperature_correction` 参数 | 已不再使用 |
| `_limit_phase_volume` 函数 | 当前实现只是裁剪，可内联 |

**boundary_bridge.py 中要删除的代码**：

| 内容 | 说明 |
|------|------|
| `build_coupled_annulus_inlet_provider` 的旧调用方式 | arg1=casing_result, arg2=casing_solver |
| `_legacy_phase_fractions_for_fluid` 函数 | 旧的相分数映射 |
| `_legacy_provider` 函数 | 旧的边界提供器 |

**各 loader 中要删除的代码**：

| 内容 | 说明 |
|------|------|
| `build_xxx_annulus_inlet_provider` 函数 | 已标记为 deprecated |
| `annulus_boundary_mode` 参数 | 旧的边界模式 |
| `export_xxx_sync_card_markdown` 函数 | 如果不再使用 |

#### 3.2.2 保持不变的部分

- `reynolds_snapshots`：湍流修正后重新使用
- `lead_field`、`tail_field`：当前仍在使用
- 套管段模型 (`CasingFlowSolver`)：不修改
- `cement_snapshots`、`lead_snapshots`、`tail_snapshots`、`spacer_snapshots`：当前仍在使用

#### 3.2.3 下游代码同步更新

删除 legacy 字段后，需要同步更新以下文件：
- `reporting/contour_plots.py`
- `reporting/plots.py`
- `reporting/reference_figures.py`
- `reporting/animation.py`
- 各 runner 文件

---

### 3.3 数据文件优化

#### 3.3.1 新增 base.py

在 `cemdisp/data/loaders/base.py` 中提取公共逻辑：

```python
"""数据加载器公共基类。

本模块提供各井数据加载器的公共方法，包括：
- CSV 文件读取
- 居中度计算
- WellSpec 构建
- FluidSpec 列表构建
- PumpingSchedule 构建

使用方式：
    各井 loader 继承或调用本模块的公共方法，只定义井特有的数据。
"""
```

**公共方法**：

| 方法 | 说明 |
|------|------|
| `read_profile_csv(csv_path)` | 读取井径/井斜 CSV |
| `calculate_standoff(depth, hole_diameter, liner_od, design_standoff)` | 计算居中度 |
| `build_well_spec(well_name, params, profile_rows)` | 构建 WellSpec |
| `build_fluids(fluid_definitions)` | 构建 FluidSpec 列表 |
| `build_pumping_schedule(steps_definitions)` | 构建 PumpingSchedule |

#### 3.3.2 各井 loader 重构

**重构后的 loader 结构**：

```python
"""呼103尾管段标准数据加载器。

本模块定义呼103井的特有数据，并调用 base.py 的公共方法组装标准数据对象。

数据来源标注：
- [实测] 现场实测数据
- [代理] 首版代理值，后续需要更新
- [估算] 根据经验公式估算
"""
```

**各井 loader 中只保留**：
- 井特有的常量定义（井段参数、流体物性、施工程序）
- 调用基类方法组装数据
- 井特有的特殊处理逻辑（如有）
- 数据来源标注

#### 3.3.3 接口统一

所有 loader 统一返回格式：

```python
def load_xxx_tailpipe() -> tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载 xxx 井尾管段完整输入。

    Returns:
        (井筒参数, 流体参数列表, 施工程序, 验证数据)
    """
```

---

### 3.4 代码注释优化

#### 3.4.1 注释规范

1. **模块级注释**：说明模块功能、主要类、使用方式
2. **类级注释**：说明类的作用、主要属性、使用示例
3. **方法级注释**：说明方法功能、参数、返回值、物理意义
4. **行内注释**：说明关键算法步骤、物理公式、取值依据

#### 3.4.2 删除旧注释

- 删除标注为"兼容旧接口"的注释
- 删除标注为"首版"但已不再适用的注释
- 删除重复或冗余的注释

#### 3.4.3 注释语言

- 使用中文注释
- 专业术语保留英文原文（如 Reynolds number、D2DGA）
- 公式使用 LaTeX 或清晰的文本表示

---

## 4. 实施顺序

### 第一阶段：湍流修正模型

1. 修改 `AnnulusD2DGASolver.__init__`，添加湍流参数
2. 新增 `_compute_reynolds_number` 方法
3. 新增 `_turbulent_viscosity_correction` 方法
4. 修改 `_compute_velocity`，集成湍流修正
5. 更新 `run` 方法，保存雷诺数快照

### 第二阶段：Legacy 代码清理

1. 清理 `AnnulusSimulationResult` 字段
2. 清理 `AnnulusD2DGASolver` 方法和参数
3. 清理 `boundary_bridge.py`
4. 清理各 loader 文件
5. 更新下游代码（reporting、runners）

### 第三阶段：数据文件优化

1. 新增 `base.py`，提取公共逻辑
2. 重构各 loader，调用基类方法
3. 添加数据来源标注
4. 统一接口格式

### 第四阶段：代码注释优化

1. 更新模块级注释
2. 更新类级注释
3. 更新方法级注释
4. 删除旧注释

---

## 5. 验证标准

### 5.1 功能验证

- [ ] 湍流修正启用时，雷诺数大于临界值的区域粘度增大
- [ ] 湍流修正禁用时，模型行为与修改前一致
- [ ] 所有井的 loader 能正常加载数据
- [ ] 所有 runner 能正常运行
- [ ] 所有 reporting 函数能正常输出图表

### 5.2 代码质量验证

- [ ] 无 legacy 兼容代码残留
- [ ] 所有注释清晰、准确、无冗余
- [ ] 各 loader 接口统一
- [ ] 公共逻辑无重复

---

## 6. 参考文献

1. Maleki, A., & Frigaard, I. A. (2019). Comparing laminar and turbulent primary cementing flows. Journal of Petroleum Science and Engineering, 175, 392-405.
2. Zhang, Z., & Frigaard, I. A. (2022). Cementing in eccentric annuli: A two-dimensional modelling approach. Journal of Fluid Mechanics, 947, A32.
3. Romero, J., & Carter, R. (1999). Gravity settling in inclined wells. SPE 55927.
4. Shah, S. N., & Sutton, D. L. (1990). Yield stress effects on settling. SPE 18036.
5. Taylor, G. I. (1953). Dispersion of soluble matter in solvent flowing slowly through a tube. Proc. R. Soc. A, 219, 186-203.
