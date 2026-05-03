# cemdisp 包文件结构参考文档

> **维护规则**：每当 `cemdisp/` 下新增或修改 `.py` 文件，必须同步更新此文档对应章节，保持文档与代码一致。

---

## 目录总览

```
cemdisp/
├── __init__.py                  # 顶层包导出
├── PACKAGE_REFERENCE.md         # 本文档
├── data/                        # 数据结构与井资料加载层
│   ├── __init__.py
│   ├── well_spec.py             # 井身结构数据类
│   ├── fluid_spec.py            # 流体物性数据类
│   ├── pumping_schedule.py      # 施工程序数据类
│   ├── validation_data.py       # 验证数据（CBL、压力等）
│   └── loaders/                 # 单井加载器
│       ├── __init__.py
│       └── hu102_loader.py      # 呼102尾管段加载器
├── models2d/                    # 环空2D顶替核心层
│   ├── __init__.py
│   ├── annulus_d2dga.py         # D2DGA偏心环空求解器
│   └── boundary_bridge.py       # 1D→2D鞋口边界桥接
├── transport1d/                 # 套管内1D输运层
│   ├── __init__.py
│   ├── casing_flow.py           # 套管内前沿追踪求解器
│   ├── interface_tracking.py    # 多流体界面追踪
│   └── pipe_exit_state.py       # 鞋口出流状态数据类
├── reporting/                   # 图表与摘要输出层
│   ├── __init__.py
│   ├── plots.py                 # 静态图表（时间序列/深度剖面/风险指标/柱状对比）
│   ├── contour_plots.py         # 云图（深度-时间等值线/截面快照/最终场三联图）
│   └── animation.py             # 动画（水泥浓度场时间演化GIF/MP4）
├── runners/                      # 各井段模型运行器
│   ├── __init__.py
│   └── hu102_tailpipe.py         # 呼102尾管段运行器（初版+1D2D耦合）
├── correlations/                # 经验相关式与快速判别层（占位）
│   └── __init__.py
├── diagnostics/                 # 顶替效率与风险指标诊断层（占位）
│   └── __init__.py
├── utils/                       # 通用工具层（占位）
│   └── __init__.py
└── validation/                  # 数值与现场对比验证层（占位）
    └── __init__.py
```

---

## 1. 顶层包 `cemdisp/__init__.py`

**功能**：统一导出包内所有公开符号，使外部可直接 `from cemdisp import WellSpec, AnnulusD2DGASolver` 等。

**导出符号（17个）**：

| 来源模块 | 导出符号 |
|----------|----------|
| `data.well_spec` | `DepthValuePoint`, `EvaluationWindow`, `WellSpec` |
| `data.fluid_spec` | `FluidRole`, `FluidSpec`, `RheologyModel` |
| `data.pumping_schedule` | `PumpingSchedule`, `PumpingScheduleStep` |
| `data.validation_data` | `ValidationData` |
| `data.loaders` | `load_hu102_tailpipe`, `build_hu102_annulus_inlet_provider`, `available_well_names` |
| `models2d.annulus_d2dga` | `AnnulusD2DGASolver`, `AnnulusSimulationResult` |
| `models2d.boundary_bridge` | `AnnulusInletState`, `build_coupled_annulus_inlet_provider` |
| `transport1d.casing_flow` | `CasingFlowSolver`, `CasingFlowResult` |
| `transport1d.interface_tracking` | `InterfaceTracker`, `InterfaceFront` |
| `transport1d.pipe_exit_state` | `PipeExitState` |

---

## 2. 数据层 `cemdisp/data/`

### 2.1 `cemdisp/data/__init__.py`

**功能**：导出数据层的8个公开数据类。

**导出符号**：`DepthValuePoint`, `EvaluationWindow`, `WellSpec`, `FluidRole`, `FluidSpec`, `RheologyModel`, `PumpingSchedule`, `PumpingScheduleStep`, `ValidationData`

---

### 2.2 `cemdisp/data/well_spec.py`

**功能**：定义井身结构相关的三个冻结（frozen）数据类，用于表达井的几何参数和评价窗口。

#### `DepthValuePoint`

- **类型**：`@dataclass(frozen=True)`
- **字段**：
  - `depth_md_m: float` — 测深（米）
  - `value: float` — 对应测深处的数值（如井径、井斜角、居中度等）
- **验证**：`__post_init__` 要求 depth_md_m 和 value 均为有限数值

#### `EvaluationWindow`

- **类型**：`@dataclass(frozen=True)`
- **字段**：
  - `name: str` — 窗口名称（如"CBL评价段"、"目的层段"）
  - `top_md_m: float` — 窗口顶界测深（米）
  - `bottom_md_m: float` — 窗口底界测深（米）
  - `window_type: str` — 窗口类型标识（如"target"、"cbl"、"full"）
- **验证**：`__post_init__` 要求 top < bottom，且两者均为正有限值

#### `WellSpec`

- **类型**：`@dataclass(frozen=True)`
- **字段**：
  - `well_name: str` — 井号名称（如"呼102"）
  - `top_md_m: float` — 模拟段顶界测深（米）
  - `bottom_md_m: float` — 模拟段底界测深（米）
  - `shoe_md_m: float` — 鞋口测深（米）
  - `hanger_md_m: float` — 悬挂器测深（米）
  - `casing_id_mm: float` — 外层套管内径（毫米）
  - `liner_od_mm: float` — 尾管外径（毫米）
  - `liner_id_mm: float` — 尾管内径（毫米）
  - `hole_diameter_profile: list[DepthValuePoint]` — 井径随深度变化剖面
  - `inclination_profile: list[DepthValuePoint]` — 井斜角随深度变化剖面（度）
  - `standoff_profile: list[DepthValuePoint]` — 居中度随深度变化剖面（百分比）
  - `evaluation_windows: list[EvaluationWindow]` — 评价窗口列表
  - `reference_root: str = ""` — 井资料来源标识
  - `notes: str = ""` — 备注
- **验证**：`__post_init__` 逐项校验深度范围合理性、管径正值、profile 非空

**内部辅助函数**：
- `_require_finite(name, value)` — 校验数值有限
- `_require_positive(name, value)` — 校验数值为正

---

### 2.3 `cemdisp/data/fluid_spec.py`

**功能**：定义流体角色枚举、流变模型枚举和流体物性数据类。

#### `RheologyModel`

- **类型**：`StrEnum`
- **成员**：
  - `NEWTONIAN` — 牛顿流体
  - `BINGHAM` — 宾汉流体（需 PV + YP）
  - `POWER_LAW` — 幂律流体（需 n + K）
  - `HERSCHEL_BULKLEY` — HB流体（需 n + K + YP）

#### `FluidRole`

- **类型**：`StrEnum`
- **成员**：
  - `MUD` — 钻井液/泥浆
  - `WASH` — 冲洗液
  - `SPACER` — 隔离液
  - `LEAD` — 领浆
  - `TAIL` — 尾浆
  - `DISPLACEMENT` — 替浆液
  - `OTHER` — 其他

#### `FluidSpec`

- **类型**：`@dataclass(frozen=True)`
- **字段**：
  - `name: str` — 流体名称
  - `role: FluidRole` — 流体角色
  - `density_kg_m3: float` — 密度（kg/m³）
  - `rheology_model: RheologyModel` — 流变模型类型
  - `plastic_viscosity_pa_s: float = 0.0` — 塑性粘度（Pa·s，宾汉/HB用）
  - `yield_stress_pa: float = 0.0` — 屈服应力（Pa，宾汉/HB用）
  - `power_law_n: float = 1.0` — 幂律指数（幂律/HB用）
  - `consistency_k: float = 0.001` — 幂律稠度系数（Pa·s^n，幂律/HB用）
- **验证**：`__post_init__` 校验密度正值、流变参数匹配模型类型

---

### 2.4 `cemdisp/data/pumping_schedule.py`

**功能**：定义施工程序步骤和整体施工时序数据类。

#### `PumpingScheduleStep`

- **类型**：`@dataclass(frozen=True)`
- **字段**：
  - `step_name: str` — 步骤名称（如"注水泥尾浆"、"替浆"）
  - `fluid_name: str` — 该步骤注入的流体名称
  - `volume_m3: float` — 注入体积（m³）
  - `rate_m3_min: float` — 注入排量（m³/min）
  - `start_time_s: float = 0.0` — 步骤开始时间（秒），可自动计算或手动指定
  - `end_time_s: float = 0.0` — 步骤结束时间（秒），可自动计算或手动指定
  - `remarks: str = ""` — 备注说明
- **验证**：`__post_init__` 校验体积正值、排量正值

#### `PumpingSchedule`

- **类型**：`@dataclass(frozen=True)`
- **字段**：
  - `steps: list[PumpingScheduleStep]` — 施工步骤列表
  - `notes: str = ""` — 备注
- **验证**：`__post_init__` 自动计算各步骤的 start_time_s 和 end_time_s（若未手动指定），确保时序连续
- **方法**：
  - `total_volume_m3() -> float` — 返回所有步骤注入体积之和

**内部辅助函数**：
- `_require_finite_non_negative(name, value)` — 校验非负有限值
- `_require_finite_positive(name, value)` — 校验正有限值

---

### 2.5 `cemdisp/data/validation_data.py`

**功能**：定义验证数据数据类，用于存储现场实测数据（CBL曲线、压力曲线等）以对比模型输出。

#### `ValidationData`

- **类型**：`@dataclass(frozen=True)`
- **字段**：
  - `cbl_csv_path: str = ""` — CBL测井数据CSV路径
  - `pressure_csv_path: str = ""` — 压力数据CSV路径
  - `temperature_csv_path: str = ""` — 温度数据CSV路径
  - `cbl_qualified_rate: float = 0.0` — CBL合格率（百分比，如66.65）
  - `notes: str = ""` — 备注

---

### 2.6 `cemdisp/data/loaders/__init__.py`

**功能**：导出加载器的公开函数。

**导出符号**：
- `load_hu102_tailpipe` — 加载呼102尾管段完整输入
- `build_hu102_annulus_inlet_provider` — 构建呼102环空入口状态提供器
- `available_well_names` — 返回当前可用井号列表

---

### 2.7 `cemdisp/data/loaders/hu102_loader.py`

**功能**：呼102井尾管段（139.70mm尾管，6823.10–7735.00m）的专用加载器，将硬编码井数据组装成标准数据类对象。此文件是**唯一允许包含单井硬编码数据**的地方。

**主要函数**：

| 函数 | 返回类型 | 说明 |
|------|----------|------|
| `load_hu102_tailpipe()` | `tuple[WellSpec, list[FluidSpec], PumpingSchedule, ValidationData]` | 返回呼102尾管段的完整四元组输入数据 |
| `build_hu102_annulus_inlet_provider(schedule, fluids, mode)` | `Callable[[float], AnnulusInletState]` | 构建环空入口状态提供器；`mode`参数决定模式：`"sustained_tail"`（硬编码持续注入水泥浆）或`"coupled"`（接收1D鞋口出流状态）|
| `available_well_names()` | `list[str]` | 返回 `["呼102"]` |

**硬编码数据清单**（仅存于此文件）：
- 井号："呼102"
- 尾管段：6823.10–7735.00m，鞋口7735m，悬挂器6823.1m
- 管径：套管内径219.10mm，尾管外径139.70mm，尾管内径108.10mm
- 井径剖面：实测CSV数据，7120–7735m段取自20215.xlsx
- 井斜剖面：实测CSV数据，7120–7735m段取自20215.xlsx
- 居中度剖面：由井径+管径+扶正器代理构造
- 流体（3种，按现场记录）：
  - 钻井液（Bingham, ρ=2.02g/cm³, PV=80mPa·s, YP=15Pa, MUD角色）
  - 替浆液（Bingham, ρ=2.02g/cm³, PV=80mPa·s, YP=15Pa, DISPLACEMENT角色）
  - 尾管水泥浆（Power-Law, ρ=2.10g/cm³, n=0.722, K=0.684, TAIL角色）
- 可选补充流体（0708邻井代理，暂不强制注入）：
  - 冲洗液（Bingham, ρ=1.88g/cm³, PV=25mPa·s, YP=1.5Pa, WASH角色）— 呼103邻井
  - 驱油隔离液（Bingham, ρ=1.85g/cm³, PV=35mPa·s, YP=8Pa, SPACER角色）— 呼103邻井
- 施工程序：注入尾管水泥浆16.67m³/1.30m³/min → 替浆液推进74.0m³/1.30m³/min
- 验证数据：CBL合格率66.65%
- 现场记录来源：10042.xlsx Row 26 (2022-11-22)，注水泥35t/2.10g/cm³，替浆液74m³/2.02g/cm³

---

## 3. 环空2D顶替核心层 `cemdisp/models2d/`

### 3.1 `cemdisp/models2d/__init__.py`

**功能**：导出求解器、结果类、入口状态和耦合函数。

**导出符号**：`AnnulusD2DGASolver`, `AnnulusSimulationResult`, `AnnulusInletState`, `build_coupled_annulus_inlet_provider`

---

### 3.2 `cemdisp/models2d/annulus_d2dga.py`

**功能**：偏心环空2D顶替求解器（D2DGA算法），基于 Zhang & Frigaard (2022, JFM) 的二维深度-周向网格模型。这是整个包的核心求解器。

#### `AnnulusSimulationResult`

- **类型**：`@dataclass(frozen=False)`（可变，因包含 numpy 数组）
- **字段**：
  - `well_name: str` — 井号
  - `geom: dict` — 环空几何信息字典（含 `y`, `phi`, `h`, `b`, `e`, `standoff`, `inc_deg`, `hole`, `od_mm`, `clearance`, `half_gap_mean`, `mean_radius`）
  - `cement_field: ndarray` — 水泥浆浓度场（shape: ny×nz，值0~1）
  - `wall_field: ndarray` — 壁面泥饼残余率场（shape: ny×nz，值0~1）
  - `metrics: DataFrame` — 逐时刻指标表（含时间、排量、效率等列）
  - `depth_profiles: dict[str, ndarray]` — 深度剖面字典（含"水泥浆占据率"、"壁面泥饼残余率"、"顶替效率"、"有效顶替效率"）
  - `summary: dict[str, float]` — 最终汇总指标（含各效率值、风险指标等）
  - `time_points_s: list[float]` — 时间节点列表（秒）
  - `cement_snapshots: Tuple[ndarray, ...]` — 水泥浓度场快照序列（每元素shape: ny×nz），按 `save_interval` 间隔保存
  - `wall_snapshots: Tuple[ndarray, ...]` — 壁面泥饼场快照序列（每元素shape: ny×nz）
  - `snapshot_times_s: Tuple[float, ...]` — 快照时间点列表（秒），与快照序列一一对应
  - `notes: list[str]` — 运行备注

#### `AnnulusD2DGASolver`

- **类型**：普通类（非dataclass）
- **构造参数**（`__init__`）：
  - `dt: float = 4.0` — 时间步长（秒）
  - `nz: int = 140` — 周向网格数
  - `ny: int = 40` — 深度网格数
  - `alpha_clean: float = 0.085` — 壁面清洁率基准系数
  - `total_t: float = 6600.0` — 总模拟时间（秒）
  - `quality_penalty_scale: float = 0.099` — 质量响应惩罚缩放因子（由CBL合格率66.65%校准）
  - `channeling_penalty_weight: float = 0.35` — 窜槽惩罚权重
  - `mixing_penalty_weight: float = 0.30` — 混浆惩罚权重
  - `instability_penalty_weight: float = 0.35` — 失稳惩罚权重
  - `instability_decay_scale: float = 0.001` — 失稳衰减尺度
  - `save_interval: int = 60` — 快照保存间隔（每隔N个时间步保存一次2D场快照）

- **核心方法**：

| 方法 | 参数 | 返回 | 说明 |
|------|------|------|------|
| `run(well_spec, fluids, inlet_state_provider)` | 井身、流体列表、入口状态函数 | `AnnulusSimulationResult` | 执行完整时间推进模拟 |
| `_build_geom(well_spec)` | 井身数据 | `dict` | 构建环空几何网格（深度y、周向phi、间隙h、宽度b） |
| `_physical_annular_volume(well_spec)` | 井身数据 | `float` | 计算物理环空体积（m³） |
| `_pick_fluids(fluids)` | 流体列表 | `tuple[FluidSpec, FluidSpec]` | 从流体列表中选出泥浆和水泥浆 |
| `_compute_props(cement, w_prev, geom, mud_fluid, cement_fluid)` | 浓度场、几何、流体 | `tuple[ndarray, ndarray, ndarray]` | 计算各网格点的表观粘度μ、密度ρ、泥浆分数 |
| `_compute_velocity(cement, geom, q_m3s, w_prev, mud_fluid, cement_fluid)` | 浓度场、几何、排量、流体 | `ndarray` | 计算周向流速场v（含浮力修正和窄间隙加速） |
| `_apparent_viscosity(fluid, gamma)` | 流体、剪切速率 | `float` | 计算单流体表观粘度 |
| `_depth_profiles(geom, cement, wall, eff)` | 几何、浓度场、壁面场、效率 | `dict` | 计算深度方向平均剖面 |

- **求解流程**（`run`方法内部）：
  1. 构建几何网格和初始场
  2. 时间循环：每个时间步调用 `inlet_state_provider(current_time)` 获取入口状态
  3. 停泵检测：若 `flow_rate < 1e-9 m³/s`，冻结水泥浆和壁面场，跳过平流/扩散/壁面清洁
  4. 计算流速场 `_compute_velocity`
  5. 半拉格朗日平流推进（bilinear插值）
  6. 显式扩散（轴向+周向Laplacian）
  7. 壁面泥饼清除（剪切驱动 + 0.45基准清除率）
  8. 计算效率指标和风险指标
  9. 最终汇总：计算质量响应效率（quality_proxy = bulk_fill × quality_factor）

- **关键指标列**（metrics DataFrame）：

| 列名 | 含义 |
|------|------|
| `时间_s` | 模拟时间（秒） |
| `排量_m3_min` | 当前排量（m³/min） |
| `泵状态` | 泵是否运行 |
| `水泥浆前沿_宽侧_m` | 宽侧水泥浆前沿位置 |
| `水泥浆前沿_窄侧_m` | 窄侧水泥浆前沿位置 |
| `水泥浆前沿_中间_m` | 中间位置前沿 |
| `窜槽指数` | 宽窄侧前沿差值 |
| `混浆指数` | 0.1~0.9过渡区面积占比 |
| `宽侧流度比` | 宽侧μ_mud/μ_cement |
| `窄侧流度比` | 窄侧μ_mud/μ_cement |
| `失稳代理值` | 浮力-粘性力比值 |
| `顶替效率` | 水泥浆占据率（bulk_fill） |
| `有效顶替效率` | 窜槽/混浆/失稳修正后的效率 |
| `质量响应效率` | 经CBL校准缩放的最终代理值 |

**内部辅助函数**：
- `_profile_to_arrays(points)` — 将 `DepthValuePoint` 列表转为深度和数值数组
- `_window_mask(well_spec, md, window_type)` — 生成评价窗口布尔掩码
- `_phase_fraction(inlet_state, phase_name)` — 从入口状态提取相分数
- `_trapez2d(arr, geom)` — 二维梯形积分
- `_bilinear_interp(field, ysrc, ssrc, geom, inlet_value)` — 双线性插值

---

### 3.3 `cemdisp/models2d/boundary_bridge.py`

**功能**：实现1D套管输运到2D环空求解的鞋口边界桥接，将套管鞋口出流状态转换为环空入口边界条件。

#### `AnnulusInletState`

- **类型**：`@dataclass(frozen=True)`
- **字段**：
  - `cement_fraction: float` — 水泥浆体积分数（0~1）
  - `flow_rate_m3_s: float` — 出口排量（m³/s）
  - `fluid_name: str` — 当前出口流体名称
  - `fluid_role: FluidRole` — 当前出口流体角色

#### `pipe_exit_to_annulus_inlet(exit_state)`

- **参数**：`PipeExitState`（来自transport1d层）
- **返回**：`AnnulusInletState`
- **映射规则**：
  - MUD / DISPLACEMENT → `cement_fraction=0.0`
  - LEAD / TAIL → `cement_fraction=1.0`
  - WASH / SPACER → `cement_fraction=0.0`（中间流体视为非水泥浆）
  - OTHER → 保持原值
  - flow_rate 和 fluid_name 直接传递

#### `build_coupled_annulus_inlet_provider(casing_result, fluids)`

- **参数**：`CasingFlowResult`, `list[FluidSpec]`
- **返回**：`Callable[[float], AnnulusInletState]` — 可传入求解器的入口状态函数
- **功能**：将1D套管求解结果包装为2D求解器可调用的时间→入口状态函数，实现1D-2D耦合

---

## 4. 套管内1D输运层 `cemdisp/transport1d/`

### 4.1 `cemdisp/transport1d/__init__.py`

**功能**：导出输运层的所有公开符号。

**导出符号**：`CasingFlowSolver`, `CasingFlowResult`, `InterfaceTracker`, `InterfaceFront`, `PipeExitState`

---

### 4.2 `cemdisp/transport1d/casing_flow.py`

**功能**：套管内1D前沿追踪求解器，基于体积守恒的解析方法追踪多流体界面在套管内的下行位置。

#### `CasingFlowResult`

- **类型**：`@dataclass(frozen=False)`
- **字段**：
  - `fronts: list[InterfaceFront]` — 各流体前沿信息列表
  - `schedule_steps: list[_ScheduledStep]` — 施工步骤时间序列
  - `pipe_cross_section_m2: float` — 套管截面积（m²）
  - `shoe_md_m: float` — 鞋口深度（米）
  - `notes: list[str]` — 备注

#### `_ScheduledStep`（内部类）

- **类型**：`@dataclass(frozen=False)`
- **字段**：
  - `step: PumpingScheduleStep` — 原始施工步骤
  - `start_time_s: float` — 实际开始时间
  - `end_time_s: float` — 实际结束时间
  - `cumulative_volume_start_m3: float` — 步骤开始时的累计注入体积
  - `cumulative_volume_end_m3: float` — 步骤结束时的累计注入体积

#### `CasingFlowSolver`

- **类型**：普通类
- **构造参数**：
  - `dt: float = 1.0` — 内部时间步长（秒），用于前沿到达时间计算

- **核心方法**：

| 方法 | 参数 | 返回 | 说明 |
|------|------|------|------|
| `run(well_spec, fluids, schedule)` | 井身、流体、施工程序 | `CasingFlowResult` | 执行1D前沿追踪，计算各流体到达鞋口的时间和位置 |
| `pipe_exit_state_at(result, time_s)` | 求解结果、查询时间 | `PipeExitState` | 查询某时刻鞋口出流状态（哪个流体、排量、累计体积）|

- **求解逻辑**：
  1. 计算套管截面积和管内容积
  2. 将施工程序转换为 `_ScheduledStep` 序列（含时间轴和累计体积轴）
  3. 对每个流体前沿，解析计算到达鞋口的时间（需累计注入体积超过管内容积）
  4. 记录各前沿到达时刻和最终位置

**内部辅助方法**：
- `_pipe_cross_section_area(well_spec)` — 从井身数据计算套管截面积
- `_build_scheduled_steps(schedule)` — 将PumpingSchedule转为时间轴序列
- `_initial_fluid_name(fluids, schedule)` — 确定管内初始流体
- `_front_arrival_time(front_step, scheduled_steps, pipe_volume_m3)` — 解析计算前沿到达鞋口时间
- `_active_step_at(scheduled_steps, time_s)` — 查询某时刻的活跃施工步骤
- `_cumulative_volume_at(scheduled_steps, time_s)` — 查询某时刻的累计注入体积
- `_fluid_by_injected_volume(scheduled_steps, volume_m3)` — 查询某注入体积对应的流体
- `_scheduled_steps_for_result(result)` — 从结果中提取施工步骤序列

---

### 4.3 `cemdisp/transport1d/interface_tracking.py`

**功能**：多流体界面追踪器，记录各流体前沿在套管内的位置和到达时间。

#### `InterfaceFront`

- **类型**：`@dataclass(frozen=False)`
- **字段**：
  - `fluid_name: str` — 流体名称
  - `arrival_time_s: float` — 到达鞋口的时间（秒，-1表示未到达）
  - `final_distance_m: float` — 最终下行距离（米）

#### `InterfaceTracker`

- **类型**：普通类
- **功能**：管理多个流体前沿的追踪，提供添加前沿和查询接口

---

### 4.4 `cemdisp/transport1d/pipe_exit_state.py`

**功能**：定义鞋口出流状态数据类，作为1D→2D边界桥接的数据载体。

#### `PipeExitState`

- **类型**：`@dataclass(frozen=True)`
- **字段**：
  - `time_s: float` — 查询时间点（秒）
  - `fluid_name: str` — 当前鞋口出流流体名称
  - `fluid_role: FluidRole` — 当前出流流体角色
  - `flow_rate_m3_s: float` — 当前出流排量（m³/s）
  - `cumulative_volume_m3: float` — 当前累计注入体积（m³）

---

## 5. 图表与摘要输出层 `cemdisp/reporting/`

### 5.1 `cemdisp/reporting/__init__.py`

**功能**：仅包含模块级docstring "图表、表格和摘要输出层"，暂无导出。

---

### 5.2 `cemdisp/reporting/plots.py`

**功能**：四类绘图函数，全部使用中文标签、中文图例、中文文件名，并处理中文字体兼容性问题。

**导出**：`__all__` 包含 `plot_time_series`, `plot_depth_profiles`, `plot_risk_indices`, `plot_efficiency_summary_bar`

#### 内部辅助函数

| 函数 | 说明 |
|------|------|
| `_setup_chinese_font()` | 设置中文字体：优先SimHei→Microsoft YaHei→系统默认；修复负号显示 |
| `_safe_filename_component(name)` | 清理文件名中的非法字符 |
| `_save_figure(fig, output_dir, filename)` | 保存图片到指定目录，dpi=150，bbox_inches=tight |
| `_require_columns(dataframe, columns, source_name)` | 校验DataFrame包含所需列 |

#### `plot_time_series(result, output_dir)`

- **参数**：`AnnulusSimulationResult`, 输出目录路径
- **输出文件**：`{井号}_时间序列结果.png`
- **内容**：双面板图表
  - 上图：排量时间曲线（含泵状态标注）
  - 下图：有效顶替效率、目标层段效率、CBL段效率随时间变化
- **脚注**：区分"有效顶替效率"与"质量响应效率"

#### `plot_depth_profiles(result, well_spec, output_dir)`

- **参数**：`AnnulusSimulationResult`, `WellSpec`, 输出目录路径
- **输出文件**：`{井号}_深度效率剖面对比.png`
- **内容**：双面板图表
  - 上图：水泥浆占据率、有效顶替效率随深度变化
  - 下图：壁面泥饼残余率随深度变化
- **评价窗口**：用阴影标注CBL段和目标层段

#### `plot_risk_indices(result, output_dir)`

- **参数**：`AnnulusSimulationResult`, 输出目录路径
- **输出文件**：`{井号}_风险指标时间序列.png`
- **内容**：单面板图表
  - 窜槽指数、混浆指数、失稳代理值随时间变化
- **脚注**：各指标含义说明

#### `plot_efficiency_summary_bar(result, output_dir)`

- **参数**：`AnnulusSimulationResult`, 输出目录路径
- **输出文件**：`{井号}_效率汇总柱状图.png`
- **内容**：单面板柱状图
  - 四根柱：顶替效率、有效顶替效率、目标层段效率、质量响应效率
  - 颜色区分：蓝/绿/橙/红
- **脚注**：明确区分各效率概念，标注"质量响应效率≠有效顶替效率真值"

---

### 5.3 `cemdisp/reporting/contour_plots.py`

**功能**：三类云图绘制函数，展示环空顶替过程的二维场分布。

**导出**：`__all__` 包含 `plot_depth_time_contour`, `plot_annulus_snapshots`, `plot_final_fields_contour`

**依赖**：复用 `plots.py` 中的 `_setup_chinese_font`, `_save_figure`, `_safe_filename_component`

#### `plot_depth_time_contour(result, output_dir)`

- **参数**：`AnnulusSimulationResult`, 输出目录路径
- **输出文件**：`{井号}_深度-时间顶替效率云图.png`
- **内容**：
  - X轴：时间（min），Y轴：深度（m），颜色：水泥浓度
  - 使用 `contourf` + `RdYlGn` 配色（红=低效，绿=高效）
  - 添加0.5/0.7/0.9等值线及标签
  - 数据来源：`cement_snapshots` 沿方位角平均后构建深度-时间矩阵

#### `plot_annulus_snapshots(result, output_dir, n_panels=6)`

- **参数**：`AnnulusSimulationResult`, 输出目录路径, 面板数
- **输出文件**：`{井号}_水泥浓度场演化过程.png`
- **内容**：
  - 多面板并排展示不同时刻的环空截面浓度场
  - 每面板：X轴=井深(m)，Y轴=方位角(宽边→窄边)，颜色=水泥浓度
  - 使用 `viridis` 配色，每面板标题 "t = XXX.X min"
  - 从 `cement_snapshots` 中均匀选取 `n_panels` 个时刻

#### `plot_final_fields_contour(result, output_dir)`

- **参数**：`AnnulusSimulationResult`, 输出目录路径
- **输出文件**：`{井号}_最终场分布.png`
- **内容**：
  - 三联图：水泥浓度(viridis) + 有效效率(RdYlGn) + 泥饼残余(YlOrRd)
  - 展示最终时刻的完整二维场分布

---

### 5.4 `cemdisp/reporting/animation.py`

**功能**：水泥浓度场时间演化动画生成。

**导出**：`__all__` 包含 `animate_cement_field`

**依赖**：复用 `plots.py` 中的 `_setup_chinese_font`, `_safe_filename_component`

#### `animate_cement_field(result, output_dir, interval_ms=200, fps=10, save_format='gif')`

- **参数**：`AnnulusSimulationResult`, 输出目录路径, 帧间隔(ms), 帧率, 保存格式(gif/mp4)
- **输出文件**：`{井号}_顶替过程动画.gif` 或 `.mp4`
- **内容**：
  - 单面板动画，每帧展示一时刻的水泥浓度场
  - X轴=井深(m)，Y轴=方位角(宽边→窄边)
  - 标题每帧更新："水泥浓度场 — t = XXX.X min"
  - 使用 `viridis` 配色，`FuncAnimation` 生成
  - GIF使用 `pillow` 写入器，MP4使用 `ffmpeg` 写入器
  - 空快照时提前返回并记录警告

---

## 6. 运行器子包 `cemdisp/runners/`

运行器负责将"加载→求解→导出→打印"全流程封装为可调函数，每个模块对应一口井的一个固井段。

### 6.1 `cemdisp/runners/__init__.py`

**功能**：子包声明与docstring，说明运行器职责范围。

### 6.2 `cemdisp/runners/hu102_tailpipe.py`

**功能**：呼102尾管段顶替效率模型运行器，执行硬编码入口和1D-2D耦合两种模式的模拟并导出中文命名结果。

**核心函数**：

| 函数 | 作用 |
|------|------|
| `run_and_export(mode_title, output_dir, inlet_provider)` | 封装加载→求解→导出→打印全流程 |
| `run_hu102_tailpipe_initial()` | 完整运行入口，依次执行初版和耦合两种模式 |

**关键常量**：

| 常量 | 值 | 说明 |
|------|----|------|
| `PROJECT_ROOT` | `Path(__file__).resolve().parents[1].parent` | cement model根目录，用于定位results/输出 |

**输出**（每模式13个文件）：
- `呼102尾管_{模式}_时间序列结果.csv` — 逐时刻指标CSV
- `呼102尾管_{模式}_深度剖面.csv` — 深度剖面CSV
- `呼102尾管_{模式}_结果摘要.json` — 完整结果JSON
- `呼102尾管_{模式}_结果摘要.md` — 结果摘要Markdown
- `呼102尾管_{模式}_2D场数据.npz` — 二维场数据（水泥浓度快照+壁面泥饼快照+网格坐标）
- 4张静态PNG图表（时间序列、深度剖面、风险指标、效率汇总柱状图）
- 3张云图PNG（深度-时间等值线、水泥浓度场演化过程、最终场分布三联图）
- 1张动画GIF（水泥浓度场时间演化）

---

## 7. 占位模块（暂无实现）

以下模块仅含docstring占位，后续按需求逐步填充：

| 模块 | docstring | 规划用途 |
|------|-----------|----------|
| `cemdisp/correlations/__init__.py` | "经验相关式与快速判别层" | 流变转换、临界排量、流态判别等经验公式 |
| `cemdisp/diagnostics/__init__.py` | "顶替效率与风险指标诊断层" | 实时诊断、风险预警、质量解释 |
| `cemdisp/utils/__init__.py` | "通用工具层" | 物理常数、单位转换、数值安全 |
| `cemdisp/validation/__init__.py` | "数值与现场对比验证层" | 质量守恒验证、CBL对比、多井切换验证 |

---

## 8. 外部脚本与测试

### 8.1 `scripts/run_hu102_tailpipe_initial.py`

**功能**：薄wrapper脚本，仅导入并调用 `cemdisp.runners.hu102_tailpipe.run_hu102_tailpipe_initial()`。

**核心逻辑已迁移至**：`cemdisp/runners/hu102_tailpipe.py`（见第6章）。

### 8.2 `tests/` 目录

| 测试文件 | 测试数 | 覆盖内容 |
|----------|--------|----------|
| `test_hu102_loader.py` | 4 | WellSpec构造、流体列表、施工程序、入口状态提供器 |
| `test_hu102_initial_model.py` | 1 | 初版模型端到端冒烟测试 |
| `test_casing_flow.py` | 5 | 套管求解器构造、步骤序列、前沿到达、出流查询、流体识别 |
| `test_boundary_bridge_coupling.py` | 3 | 入口状态映射、耦合提供器构建、流体角色→相分数 |

---

## 9. 关键设计决策备忘

| 决策 | 内容 | 原因 |
|------|------|------|
| 冻结dataclass | 所有数据类使用 `frozen=True` | 防止意外修改输入数据，确保求解器收到不可变输入 |
| 停泵冻结 | `flow_rate < 1e-9` 时冻结水泥浆和壁面场 | 水泥浆静态胶凝强度在短停泵期抵抗浮力滑落；防止v_eff和壁面清洁基准项导致水泥"蒸发" |
| 单相水泥 | 初版仅追踪水泥浆（TAIL）浓度场（0/1） | 缺乏冲洗液/隔离液/领浆的直接数据证据，待0708资料确认后扩展 |
| 1D纯前沿追踪 | 套管内不考虑管内混合/扩散 | 首版简化，后续可加入对流-弥散而无需重写环空核心 |
| 入口提供器 | 函数式 `Callable[[float], AnnulusInletState]` | 解耦1D与2D，支持硬编码和耦合两种模式切换 |
| 质量响应效率 | `quality_proxy = bulk_fill × quality_factor`，α=0.099 | 由Hu102 CBL合格率66.65%校准，明确标注为代理值而非顶替效率真值 |

---

> **文档版本**：2026-05-03 | 基于 cemdisp/ 当前25个 .py 文件编写（含新增 runners/ 子包、contour_plots.py、animation.py）