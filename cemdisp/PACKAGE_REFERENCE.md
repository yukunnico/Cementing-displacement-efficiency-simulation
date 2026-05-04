# cemdisp 包文件结构参考文档

> **维护规则**：每当 `cemdisp/` 下新增或修改 `.py` 文件，必须同步更新此文档对应章节，保持文档与代码一致。

---

## 目录总览

```
cemdisp/
├── __init__.py                  # 顶层包导出
├── data/                        # 数据结构与井资料加载层
│   ├── __init__.py
│   ├── well_spec.py             # 井身结构数据类
│   ├── fluid_spec.py            # 流体物性数据类
│   ├── pumping_schedule.py      # 施工程序数据类
│   ├── validation_data.py       # 验证数据（CBL、压力等）
│   └── loaders/                 # 单井加载器
│       ├── __init__.py
│       ├── hu101_loader.py      # 呼101尾管段加载器
│       └── hu102_loader.py      # 呼102尾管段加载器
├── transport1d/                 # 套管内1D输运层
│   ├── __init__.py
│   ├── casing_flow.py           # 套管内前沿追踪求解器
│   ├── interface_tracking.py    # 多流体界面追踪
│   └── pipe_exit_state.py       # 鞋口出流状态数据类
├── models2d/                    # 环空2D顶替核心层
│   ├── __init__.py
│   ├── annulus_d2dga.py         # D2DGA偏心环空求解器
│   └── boundary_bridge.py       # 1D→2D鞋口边界桥接
├── reporting/                   # 图表与摘要输出层
│   ├── __init__.py
│   ├── plots.py                 # 静态图表（时间序列/深度剖面/风险指标/柱状对比）
│   ├── contour_plots.py         # 二维云图（深度-时间/快照/最终场）
│   └── animation.py             # 浓度场演化动画
├── runners/                     # 井段模型运行器
│   ├── __init__.py
│   ├── hu101_tailpipe.py        # 呼101尾管段运行器
│   └── hu102_tailpipe.py        # 呼102尾管段运行器
├── correlations/                # 经验相关式与快速判别层
│   └── __init__.py
├── diagnostics/                 # 顶替效率与风险指标诊断层
│   └── __init__.py
├── utils/                       # 通用工具层
│   └── __init__.py
└── validation/                  # 数值与现场对比验证层
    └── __init__.py
```

---

## 1. 顶层包 `cemdisp/__init__.py`

**功能**：统一导出包内所有公开符号，使外部可直接 `from cemdisp import WellSpec, AnnulusD2DGASolver` 等。

**主要子模块说明**：
- `cemdisp.data`: 标准输入数据结构（井筒规格、流体物性、施工程序）
- `cemdisp.transport1d`: 套管内一维输运层（流体前缘追踪、鞋口出流状态）
- `cemdisp.models2d`: 环空二维顶替核心层（偏心环空顶替模拟）
- `cemdisp.diagnostics`: 顶替效率与风险指标诊断层
- `cemdisp.reporting`: 图表与报告输出（中文图表、云图、动画）
- `cemdisp.correlations`: 经验相关式
- `cemdisp.utils`: 通用工具
- `cemdisp.runners`: 井段模型运行器
- `cemdisp.validation`: 数值验证

**导出符号**：

| 来源模块 | 导出符号 |
|----------|----------|
| `data.well_spec` | `DepthValuePoint`, `EvaluationWindow`, `WellSpec` |
| `data.fluid_spec` | `FluidRole`, `FluidSpec`, `RheologyModel` |
| `data.pumping_schedule` | `PumpingSchedule`, `PumpingScheduleStep` |
| `data.validation_data` | `ValidationData` |
| `models2d.annulus_d2dga` | `AnnulusD2DGASolver`, `AnnulusSimulationResult` |
| `models2d.boundary_bridge` | `AnnulusInletState`, `pipe_exit_to_annulus_inlet` |
| `transport1d.casing_flow` | `CasingFlowResult` |
| `transport1d.interface_tracking` | `InterfaceFront` |
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

**功能**：定义现场校验资料路径集合，用于存储与单井相关的现场实测资料路径。

#### `ValidationData`

- **类型**：`@dataclass(frozen=True)`
- **字段**：
  - `cbl_profile_path: Path | None` — CBL（水泥胶结测井）剖面原始数据文件路径
  - `cbl_summary_path: Path | None` — CBL评价汇总报告文件路径
  - `job_report_path: Path | None` — 固井施工总结报告文件路径
  - `pump_pressure_series_path: Path | None` — 泵压时序数据文件路径
  - `returns_report_path: Path | None` — 返排记录报告文件路径
  - `cement_top_path: Path | None` — 水泥面位置确认文件路径
  - `notes: Tuple[str, ...]` — 备注信息元组

- **说明**：这些路径指向现场实测资料，用于与模型预测结果进行对比验证。字段可为 `None`，表示该资料尚未提供或不适用。

---

### 2.6 `cemdisp/data/loaders/__init__.py`

**功能**：导出加载器的公开函数。

**导出符号**：
- `load_hu101_tailpipe` — 加载呼101尾管段完整输入
- `load_hu102_tailpipe` — 加载呼102尾管段完整输入
- `build_hu101_annulus_inlet_provider` — 构建呼101环空入口状态提供器
- `build_hu102_annulus_inlet_provider` — 构建呼102环空入口状态提供器
- `available_well_names` — 返回当前可用井号列表
- `REFERENCE_DOCS_ROOT` — 参考文档根目录路径
- `Hu101LoaderResult`, `Hu102LoaderResult` — 加载结果类型别名
- `Hu101InletProviderFactory`, `Hu102InletProviderFactory` — 入口提供器工厂类型别名

---

### 2.7 `cemdisp/data/loaders/hu101_loader.py`

**功能**：呼101井尾管段（上部168.3mm尾管+下部139.70mm尾管，5400–7868.00m）的专用加载器，将硬编码井数据组装成标准数据类对象。此文件是**唯一允许包含单井硬编码数据**的地方。

**主要函数**：

| 函数 | 返回类型 | 说明 |
|------|----------|------|
| `load_hu101_tailpipe(reference_root=None)` | `Tuple[WellSpec, Tuple[FluidSpec, ...], PumpingSchedule, ValidationData]` | 返回呼101尾管段的完整四元组输入数据 |
| `build_hu101_annulus_inlet_provider(schedule, fluids, annulus_boundary_mode="field_order_realistic", *, split_cement_phases=False)` | `Callable[[float], AnnulusInletState]` | 构建环空入口状态提供器；支持多种边界模式：`"sustained_tail"`、`"volume_limited"`、`"tail_then_mud"`、`"field_order_realistic"`；支持`split_cement_phases`参数分离领浆/尾浆相 |

**硬编码数据清单**（仅存于此文件）：
- 井号："呼101"
- 尾管段：5400–7868m，鞋口7868m，悬挂器5407.46m
- 管径：技术套管当量内径273.10mm，上部尾管外径168.30mm（上部井段：5400–6796m），下部尾管外径139.70mm（下部井段：6796–7868m），尾管壁厚15.80mm，尾管等效内径由52m³鞋口滞后反推
- 井径剖面：由实测数据和段近似构造，上部段用等效井径保持环空面积
- 井斜剖面：由实测数据点构造，5400–7868m
- 居中度剖面：由实测数据点构造
- 流体（8种，按现场记录）：
  - 钻井液（Bingham, ρ=1960kg/m³, PV=58mPa·s, YP=5Pa, MUD角色）
  - 平衡液（Bingham, ρ=1850kg/m³, PV=30mPa·s, YP=3Pa, WASH角色）
  - 驱油隔离液（Bingham, ρ=2000kg/m³, PV=30mPa·s, YP=5Pa, SPACER角色）
  - 领浆（Power-Law, ρ=2100kg/m³, n=0.719, K=0.815, LEAD角色）
  - 尾浆（Power-Law, ρ=1900kg/m³, n=0.722, K=0.684, TAIL角色）
  - 轻泥浆（Bingham, MUD/DISPLACEMENT角色）
  - 中置液（Bingham, DISPLACEMENT角色）
  - 井浆（Bingham, DISPLACEMENT角色）
- 施工程序（9步，按现场记录）：
  1. 注平衡液：25.0m³ / 1.20m³/min
  2. 注驱油隔离液：25.0m³ / 1.20m³/min
  3. 注领浆：47.0m³ / 1.20m³/min
  4. 注尾浆：23.0m³ / 1.20m³/min
  5. 注后置液：2.0m³ / 0.60m³/min
  6. 注轻泥浆：26.0m³ / 1.50m³/min
  7. 注中置液：10.0m³ / 1.20m³/min
  8. 井浆快替：40.0m³ / 1.00m³/min
  9. 井浆慢替：23.4m³ / 0.55m³/min
- 验证数据：
  - CBL合格率62.77%，评价井段5700–7810m

---

### 2.8 `cemdisp/data/loaders/hu102_loader.py`

**功能**：呼102井尾管段（139.70mm尾管，6823.10–7735.00m）的专用加载器，将硬编码井数据组装成标准数据类对象。此文件是**唯一允许包含单井硬编码数据**的地方。

**主要函数**：

| 函数 | 返回类型 | 说明 |
|------|----------|------|
| `load_hu102_tailpipe()` | `Tuple[WellSpec, Tuple[FluidSpec, ...], PumpingSchedule, ValidationData]` | 返回呼102尾管段的完整四元组输入数据 |
| `build_hu102_annulus_inlet_provider(schedule, fluids, mode)` | `Callable[[float], AnnulusInletState]` | 构建环空入口状态提供器；`mode`参数决定模式：`"sustained_tail"`（替浆期环空入口保持尾浆）、`"volume_limited"`（管内推进期间环空入口保持尾浆，排量设为0）、`"tail_then_mud"`（替浆期结束后环空入口切换为钻井液）|

**硬编码数据清单**（仅存于此文件）：
- 井号："呼102"
- 尾管段：6823.10–7735.00m，鞋口7735m，悬挂器6823.1m
- 管径：套管内径219.10mm，尾管外径139.70mm，尾管壁厚15.80mm，尾管内径108.10mm
- 井径剖面：实测CSV数据，7120–7735m段取自20215.xlsx
- 井斜剖面：实测CSV数据，7120–7735m段取自20215.xlsx
- 居中度剖面：由井径+管径按分段规则构造
- 流体（2种，按现场记录）：
  - 钻井液（Bingham, ρ=2.02g/cm³, PV=80mPa·s, YP=15Pa, MUD角色）
  - 尾管水泥浆（Power-Law, ρ=2.10g/cm³, n=0.722, K=0.684, TAIL角色）
- 施工程序（首版简化两步）：
  1. 注入尾管水泥浆：16.67m³ / 1.30m³/min
  2. 替浆推进：74.0m³ / 1.30m³/min
- 验证数据：
  - CBL合格率66.65%，评价井段6840–7665m
  - 水泥浆35t/平均密度2.10g/cm³
  - 替浆量74m³

---

## 3. 环空2D顶替核心层 `cemdisp/models2d/`

### 3.1 `cemdisp/models2d/__init__.py`

**功能**：导出求解器、结果类、入口状态和边界桥接函数。

**导出符号**：`AnnulusD2DGASolver`, `AnnulusSimulationResult`, `AnnulusInletState`, `pipe_exit_to_annulus_inlet`, `build_coupled_annulus_inlet_provider`

**子模块说明**：
- 环空二维求解器可独立运行，无需依赖套管内1D层
- 1D→2D边界桥接接口清晰，支持后续扩展

---

### 3.2 `cemdisp/models2d/annulus_d2dga.py`

**功能**：多井通用环空二维D2DGA（双液双区面积）求解器首版实现。

核心物理过程：
1. **平流输运**：水泥浆在偏心环空中随平均流场运移
2. **弥散效应**：横向和垂向扩散（与偏心度、流体性质相关）
3. **壁面清除**：泥饼在水泥浆剪切清除作用下逐渐移除
4. **浮力效应**：密度差导致的窄边优先顶替和失稳风险

网格系统：
- `s`方向：井深方向（从鞋口到悬挂器）
- `y`方向：方位角方向（0=窄边，π=宽边）
- `φ`方向：归一化方位角（0~1对应窄边到宽边）

#### `AnnulusSimulationResult`

- **类型**：`@dataclass(frozen=True)`
- **字段**：
  - `well_name: str` — 井号
  - `geom: Dict[str, Array]` — 几何参数字典（含s、md、y、phi、H、b、e、standoff、inc_deg、hole_mm、od_mm等）
  - `cement_field: Array` — 水泥浓度场（shape: ny×nz，值0~1）
  - `wall_field: Array` — 壁面泥饼清除场（shape: ny×nz，值0~1）
  - `metrics: pd.DataFrame` — 时间序列指标DataFrame
  - `depth_profiles: pd.DataFrame` — 深度方向平均剖面DataFrame
  - `summary: Dict[str, object]` — 最终结果摘要字典
  - `time_points_s: Tuple[float, ...]` — 时间点序列
  - `notes: Tuple[str, ...]` — 运行备注

#### `AnnulusD2DGASolver`

- **类型**：普通类
- **构造参数**（`__init__`）：
  - `dt: float = 4.0` — 时间步长（秒）
  - `nz: int = 140` — 井深方向网格数
  - `ny: int = 40` — 方位角方向网格数
  - `alpha_clean: float = 0.085` — 泥饼清除系数
  - `total_t: float = 6600.0` — 总模拟时间（秒），默认6600秒（110分钟）
  - `quality_penalty_scale: float = 0.099` — 质量惩罚因子
  - `channeling_penalty_weight: float = 0.55` — 窜槽风险权重
  - `mixing_penalty_weight: float = 0.35` — 混浆风险权重
  - `instability_penalty_weight: float = 0.25` — 失稳风险权重
  - `instability_decay_scale: float = 5.0` — 失稳指数衰减标度

- **核心方法**：

| 方法 | 参数 | 返回 | 说明 |
|------|------|------|------|
| `run(well_spec, fluids, inlet_state_provider)` | 井身、流体列表、入口状态函数 | `AnnulusSimulationResult` | 执行完整时间推进模拟 |
| `_build_geom(well_spec)` | 井身数据 | `dict` | 构建环空几何网格 |
| `_physical_annular_volume(well_spec)` | 井身数据 | `float` | 计算物理环空体积 |
| `_pick_fluids(fluids)` | 流体列表 | `tuple[FluidSpec, FluidSpec]` | 选取钻井液和水泥浆 |
| `_apparent_viscosity(fluid, gamma)` | 流体、剪切速率 | `Array` | 计算单流体表观粘度 |
| `_compute_props(cement, w_prev, geom, mud_fluid, cement_fluid)` | 浓度场、几何、流体 | `tuple[Array, Array, Array]` | 计算混合物系表观粘度、密度、钻井液分数 |
| `_compute_velocity(cement, geom, q_m3s, w_prev, mud_fluid, cement_fluid)` | 浓度场、几何、排量、流体 | `tuple[Array, Array, Array, Array, Array]` | 计算环空速度场 |
| `_depth_profiles(geom, cement, wall)` | 几何、浓度场 | `pd.DataFrame` | 计算深度方向平均剖面 |

- **求解流程**（`run`方法内部）：
  1. 构建几何网格和初始场（初始水泥浓度=0，壁面泥饼=1）
  2. 时间循环：每个时间步调用 `inlet_state_provider(current_time)` 获取入口状态
  3. 停泵检测：若 `flow_rate < 1e-9 m³/s`，冻结水泥浆和壁面场
  4. 计算流速场（井深方向速度w、方位角方向速度v）
  5. 半拉格朗日平流推进（双线性插值）
  6. 显式扩散（轴向+方位角Laplacian）
  7. 壁面泥饼清除（剪切驱动）
  8. 计算效率指标和风险指标
  9. 最终汇总

- **关键指标列**（metrics DataFrame）：

| 列名 | 含义 |
|------|------|
| `time_s` | 模拟时间（秒） |
| `time_min` | 模拟时间（分钟） |
| `stage` | 当前施工阶段名称 |
| `bulk_cement_fill` | 水泥浆占据率 |
| `effective_efficiency` | 全井段有效顶替效率 |
| `target_interval_efficiency` | 目标层段有效顶替效率 |
| `cbl_eval_interval_efficiency` | CBL评价井段有效顶替效率 |
| `cbl_quality_proxy` | 质量响应效率代理值 |
| `front_wide_m` | 宽边水泥浆前沿位置 |
| `front_narrow_m` | 窄边水泥浆前沿位置 |
| `front_mid_m` | 中线水泥浆前沿位置 |
| `channeling_index` | 窜槽指数 |
| `mixing_index` | 混浆指数 |
| `instability_proxy` | 失稳代理值 |
| `instability_index` | 失稳指数 |
| `mean_wall_mud` | 壁面平均泥饼残余率 |
| `mean_cement` | 平均水泥浓度 |
| `mean_mud` | 平均钻井液浓度 |

**内部辅助函数**：
- `_profile_to_arrays(points)` — 将 `DepthValuePoint` 列表转为深度和数值数组
- `_window_mask(well_spec, md, window_type)` — 生成评价窗口布尔掩码
- `_phase_fraction(inlet_state, phase_name)` — 从入口状态提取相分数
- `_trapez2d(arr, geom)` — 二维梯形积分
- `_bilinear_interp(field, ysrc, ssrc, geom, inlet_value)` — 双线性插值

---

### 3.3 `cemdisp/models2d/boundary_bridge.py`

**功能**：实现1D套管输运到2D环空求解的鞋口边界桥接。

#### `AnnulusInletState`

- **类型**：`@dataclass(frozen=True)`
- **字段**：
  - `time_s: float` — 当前时间（秒）
  - `flow_rate_m3_s: float` — 排量（立方米/秒）
  - `stage_name: str` — 施工阶段名称
  - `phase_fractions: tuple[tuple[str, float], ...]` — 各相体积分数元组

#### `pipe_exit_to_annulus_inlet(pipe_exit_state)`

- **功能**：直接将鞋口出流状态映射为环空入口边界状态
- **参数**：`PipeExitState`（来自transport1d层）
- **返回**：`AnnulusInletState`

#### `build_coupled_annulus_inlet_provider(casing_result, casing_solver, fluids)`

- **功能**：构建1D-2D耦合的环空入口边界提供器
- **参数**：`CasingFlowResult`, `CasingFlowSolver`, `tuple[FluidSpec, ...]`
- **返回**：`Callable[[float], AnnulusInletState]`
- **说明**：将套管鞋口出流状态转换为环空入口，自动将流体名称映射为环空两相（cement/mud）

---

## 4. 套管内1D输运层 `cemdisp/transport1d/`

### 4.1 `cemdisp/transport1d/__init__.py`

**功能**：导出输运层的所有公开符号。

**导出符号**：`CasingFlowSolver`, `CasingFlowResult`, `InterfaceTracker`, `InterfaceFront`, `PipeExitState`

**子模块说明**：
- 套管内1D输运层提供更真实的鞋口出流边界，不污染环空2D核心
- 环空求解器可独立运行，无需依赖套管内1D层
- 1D→2D边界桥接接口清晰，支持后续扩展

---

### 4.2 `cemdisp/transport1d/casing_flow.py`

**功能**：套管内1D前沿追踪求解器，基于体积守恒的解析方法追踪多流体界面在套管内的下行位置。

#### `CasingFlowResult`

- **类型**：`@dataclass(frozen=True)`
- **字段**：
  - `fronts: tuple[InterfaceFront, ...]` — 各流体前缘元组
  - `schedule_steps: tuple[PumpingScheduleStep, ...]` — 施工步骤记录
  - `pipe_cross_section_m2: float` — 管内截面积（平方米）
  - `shoe_md_m: float` — 鞋口深度（米）
  - `notes: tuple[str, ...]` — 运行备注

#### `_ScheduledStep`（内部类）

- **类型**：`@dataclass(frozen=True)`
- **字段**：
  - `step: PumpingScheduleStep` — 原始施工步骤
  - `start_time_s: float` — 实际开始时间（秒）
  - `end_time_s: float` — 实际结束时间（秒）
  - `cumulative_volume_start_m3: float` — 步骤开始时的累计注入体积
  - `cumulative_volume_end_m3: float` — 步骤结束时的累计注入体积

#### `CasingFlowSolver`

- **类型**：普通类
- **构造参数**：
  - `dt: float = 2.0` — 时间步长（秒），仅用于内部时间查询容差判断

- **核心方法**：

| 方法 | 参数 | 返回 | 说明 |
|------|------|------|------|
| `run(well_spec, fluids, schedule)` | 井身、流体、施工程序 | `CasingFlowResult` | 执行1D前沿追踪 |
| `pipe_exit_state_at(result, time_s)` | 求解结果、查询时间 | `PipeExitState` | 查询某时刻鞋口出流状态 |

**模型假设**：
- 套管内流体以施工程序给定排量注入
- 流体前缘按体积推进法追踪，首版不考虑管内扩散
- 鞋口深度为从地面到鞋口的总测深
- 前缘到达鞋口后，对应流体从鞋口进入环空

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

**功能**：套管内流体前缘追踪数据结构和轻量追踪器。

#### `InterfaceFront`

- **类型**：`@dataclass(frozen=True)`
- **字段**：
  - `fluid_name: str` — 流体名称
  - `distance_m: float` — 从地面算起的前缘深度（米）
  - `time_s: float` — 到达该位置的时间（秒）

#### `InterfaceTracker`

- **类型**：普通类
- **构造参数**：
  - `shoe_depth_m: float` — 鞋口深度（米）
  - `pipe_area_m2: float` — 管内截面积（平方米）

- **核心方法**：
  - `advance_front(fluid_name, rate_m3_s, dt)` — 推进指定流体的前缘
  - `fronts_snapshot(time_s)` — 获取所有前缘的当前状态快照

- **设计说明**：
  - 主求解器（CasingFlowSolver）当前采用解析计算
  - InterfaceTracker 保留用于可视化或非恒定排量细化场景
  - 位移增量 = 排量(m³/s) × 时间步长(s) / 管内截面积(m²)

---

### 4.4 `cemdisp/transport1d/pipe_exit_state.py`

**功能**：定义鞋口出流状态数据结构。

#### `PipeExitState`

- **类型**：`@dataclass(frozen=True)`
- **字段**：
  - `time_s: float` — 当前时刻（秒）
  - `flow_rate_m3_s: float` — 当前排量（立方米/秒）
  - `stage_name: str` — 当前施工阶段名称
  - `phase_fractions: Tuple[Tuple[str, float], ...]` — 各相流体体积分数元组

- **用途**：作为套管内1D输运层向环空2D层传递的边界信息

---

## 5. 井段模型运行器 `cemdisp/runners/`

### 5.1 `cemdisp/runners/__init__.py`

**功能**：提供各井段顶替效率模型的执行与导出功能。运行器负责加载井段输入数据、配置求解器与边界条件、执行模拟并导出结果。

**导出符号**：
- `run_hu101_tailpipe_initial` — 运行呼101尾管段初始模拟
- `run_hu102_tailpipe_initial` — 运行呼102尾管段初始模拟

### 5.2 `cemdisp/runners/hu101_tailpipe.py`

**功能**：呼101井尾管段模型运行器，封装从数据加载到模拟执行再到结果导出的完整流程。

**主要函数**：

| 函数 | 返回类型 | 说明 |
|------|----------|------|
| `run_hu101_tailpipe_initial(output_dir=None)` | `AnnulusSimulationResult` | 运行呼101尾管段初始模拟并返回结果，支持指定输出目录 |

### 5.3 `cemdisp/runners/hu102_tailpipe.py`

**功能**：呼102井尾管段模型运行器，封装从数据加载到模拟执行再到结果导出的完整流程。

**主要函数**：

| 函数 | 返回类型 | 说明 |
|------|----------|------|
| `run_hu102_tailpipe_initial(output_dir=None)` | `AnnulusSimulationResult` | 运行呼102尾管段初始模拟并返回结果，支持指定输出目录 |

---

## 6. 图表与摘要输出层 `cemdisp/reporting/`

### 6.1 `cemdisp/reporting/__init__.py`

**功能**：图表、表格和摘要输出层，提供固井顶替模型的通用中文图表输出功能。

**字体配置**：
- 优先使用 SimHei（黑体）、Microsoft YaHei（微软雅黑）
- 回退到 SimSun（宋体）或系统默认 sans-serif
- 负号显示正常（axes.unicode_minus=False）

---

### 6.2 `cemdisp/reporting/plots.py`

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

- **参数**：`AnnulusSimulationResult`, 可选的输出目录路径
- **输出文件**：`{井号}_顶替效率时间序列.png`
- **内容**：双面板图表
  - 上图：全井段、CBL评价井段、目标层段的有效顶替效率及水泥浆占据率
  - 下图：宽边、中线、窄边三个方位的水泥浆前沿推进距离
- **脚注**：区分"有效顶替效率"与"质量响应效率"

#### `plot_depth_profiles(result, well_spec, output_dir)`

- **参数**：`AnnulusSimulationResult`, `WellSpec`（预留）, 可选的输出目录路径
- **输出文件**：`{井号}_深度剖面分布.png`
- **内容**：双面板图表
  - 上图：宽边、中线、窄边及平均有效顶替效率沿井深分布
  - 下图：宽边、中线、窄边及平均水泥浓度沿井深分布

#### `plot_risk_indices(result, output_dir)`

- **参数**：`AnnulusSimulationResult`, 可选的输出目录路径
- **输出文件**：`{井号}_风险指标时间演变.png`
- **内容**：单面板图表
  - 窜槽指数：宽窄边水泥前沿差异
  - 混浆指数：水泥-钻井液混合程度
  - 失稳指数：浮力导致的窄边水泥滑塌风险

#### `plot_efficiency_summary_bar(result, output_dir)`

- **参数**：`AnnulusSimulationResult`, 可选的输出目录路径
- **输出文件**：`{井号}_最终结果对比.png`
- **内容**：单面板柱状图
  - 五根柱：全井段有效顶替效率、CBL评价井段有效顶替效率、目标层段有效顶替效率、水泥浆占据率、质量响应效率
- **脚注**：明确标注"质量响应效率≠有效顶替效率"

---

### 6.3 `cemdisp/reporting/contour_plots.py`

**功能**：提供面向报告与论文插图的二维场可视化函数，包括深度-时间顶替效率云图、环空水泥与隔离液浓度场多时刻快照图、最终水泥浓度/隔离液浓度/有效顶替效率/壁面泥饼场分布图。全部使用中文标签、中文图例、中文文件名，并复用 plots.py 中的中文字体与文件名清理逻辑。

**导出**：`__all__` 包含 `plot_depth_time_contour`, `plot_annulus_snapshots`, `plot_final_fields_contour`

#### `plot_depth_time_contour(result, output_dir)`

- **参数**：`AnnulusSimulationResult`, 可选的输出目录路径
- **输出文件**：`{井号}_深度-时间顶替效率云图.png`
- **内容**：深度-时间二维云图
  - 颜色表示有效顶替效率，按时间快照沿方位角平均后的深度剖面组装
  - 等高线标注关键效率阈值（0.5, 0.7, 0.9）
  - x轴：时间（分钟），y轴：井深（米）

#### `plot_annulus_snapshots(result, output_dir, n_panels=6)`

- **参数**：`AnnulusSimulationResult`, 可选的输出目录路径, 可选的快照面板数
- **输出文件**：`{井号}_水泥-隔离液浓度场演化过程.png`
- **内容**：多时刻快照对照图
  - 上排：各时刻水泥浓度场（viridis色图）
  - 下排：各时刻隔离液浓度场（coolwarm色图）
  - 每个快照显示归一化方位角（0=宽边→1=窄边）与井深
  - 支持旧版单流体结果兼容（无隔离液快照时显示零浓度）

#### `plot_final_fields_contour(result, output_dir)`

- **参数**：`AnnulusSimulationResult`, 可选的输出目录路径
- **输出文件**：`{井号}_最终场分布.png`
- **内容**：四面板最终场分布图
  - 面板1：水泥浓度场
  - 面板2：隔离液浓度场
  - 面板3：有效顶替效率场
  - 面板4：壁面泥饼场
  - 各场独立色标，避免混淆

---

### 6.4 `cemdisp/reporting/animation.py`

**功能**：提供面向用户报告的动画生成功能，用于展示环空二维模型中水泥与隔离液浓度场随施工时间推进的演化过程。全部使用中文标签、中文图例、中文文件名，并复用 plots.py 中的中文字体与文件名清理逻辑。

**导出**：`__all__` 包含 `animate_cement_field`

#### `animate_cement_field(result, output_dir=None, interval_ms=200, fps=10, save_format="gif")`

- **参数**：`AnnulusSimulationResult`, 可选的输出目录路径, 帧间隔（毫秒）, 帧率, 保存格式（"gif"或"mp4"）
- **输出文件**：`{井号}_顶替过程动画.gif` 或 `{井号}_顶替过程动画.mp4`
- **内容**：左右双面板动态热力图
  - 左面板：水泥浓度场演化（viridis色图）
  - 右面板：隔离液浓度场演化（coolwarm色图）
  - 标题显示当前时间
  - 支持旧版单流体结果兼容（无隔离液快照时显示零浓度）
