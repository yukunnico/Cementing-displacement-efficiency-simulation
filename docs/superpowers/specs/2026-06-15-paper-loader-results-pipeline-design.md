# 论文专用 Paper Loader 与结果流水线设计

## 0. 文档目的

本文档设计一套论文专用数据加载与模拟结果流水线，用于支持 `cement model` 项目中《石油勘探与开发》版本论文的数据生成、模型验证和图表输出。

该设计的核心目标是：

> 将论文采用的现场数据、模型输入、模拟结果和论文图表从历史 loader 与历史 results 中隔离出来，形成统一、可追溯、可复现的论文数据口径。

本文档只描述设计，不实施代码修改。

---

## 1. 背景与问题

当前 `cemdisp/data/loaders` 中已有多个井的 loader，包括：

- `hu101_loader.py`
- `hu102_loader.py`
- `hu103_loader.py`
- `hu1_loader.py`
- `ht1_001_loader.py`
- `hu2_loader.py`
- `ht1_003_loader.py`
- `ht1_004_loader.py`

这些 loader 支撑了历史模拟、调试、敏感性分析和部分现场案例。但根据现场提取数据与 loader 核对结果，现有 loader 存在以下共性问题：

1. **现场值、模型代理值、优化值混杂**。部分参数来自现场资料，部分来自等效几何或经验假设，部分来自优化分段。
2. **旧 results 目录口径混杂**。历史结果可能来自不同 loader、不同时间步、不同网格和不同边界条件。
3. **复合尾管处理不统一**。例如 hu101、hu103 现场为复合尾管，但部分 loader 使用单外径或等效几何。
4. **CBL评价窗、目标层段和模型关注窗语义混用**。这会影响论文验证口径。
5. **部分井缺独立现场核验**。例如 ht1_003、ht1_004 当前更适合作为应用或敏感性案例，不宜写成现场定量验证井。

如果直接修改旧 loader 并继续输出到旧 `results/`，论文结果将难以追溯，也容易混用不同数据口径。

---

## 2. 设计目标

本设计要实现以下目标：

1. **隔离论文输入**：新增 paper loader 层，保留旧 loader 不动。
2. **隔离论文结果**：新增 `results_paper/pead_v1/`，不混用历史 `results/`。
3. **统一论文参数**：所有论文主结果使用同一模型配置。
4. **保留数据来源**：每口井记录现场值、解释值、模型假设和优化值。
5. **支持多类结果**：主结果、CBL验证、D2DGA消融、网格/时间步验证分目录输出。
6. **服务论文写作**：输出论文可直接引用的 CSV、Markdown 表格和图件索引。

---

## 3. 非目标

本设计暂不做以下事情：

1. 不删除或覆盖历史 `results/`。
2. 不直接修改旧 loader 的历史口径。
3. 不把所有现场提取CSV自动转为模型输入。
4. 不解决完整三维CFD或热-化学-力学耦合问题。
5. 不将CBL合格率等同于水力顶替效率。

---

## 4. 推荐总体架构

建议新增一条论文专用流水线：

```text
现场资料提取 CSV
        ↓
论文采用值对照表 adopted values
        ↓
paper loader
        ↓
paper runner
        ↓
results_paper/pead_v1
        ↓
论文数据表、图件和正文
```

旧流水线保留：

```text
legacy loader → legacy runner → results/
```

新旧结果的关系：

| 层级 | 用途 | 是否进论文主表 |
|---|---|---|
| legacy loader/results | 历史结果、开发参考、对比排查 | 默认不进入 |
| paper loader/results | 论文主数据口径 | 进入 |
| paper ablation results | 方法对比 | 进入方法验证 |
| paper grid/dt results | 数值可靠性 | 进入可靠性验证 |

---

## 5. 新增目录设计

### 5.1 Paper loader 层

推荐新增：

```text
cemdisp/data/loaders/paper/
├── __init__.py
├── common.py
├── registry.py
├── paper_hu101_loader.py
├── paper_hu102_loader.py
├── paper_hu103_loader.py
└── paper_ht1_001_loader.py
```

后续视资料补齐情况再加入：

```text
paper_hu2_loader.py
paper_ht1_003_loader.py
paper_ht1_004_loader.py
```

### 5.2 Paper runner 脚本层

推荐新增：

```text
scripts/paper_data/
├── run_paper_main_results.py
├── run_paper_d2dga_ablation.py
├── run_paper_grid_dt_validation.py
├── build_cbl_validation_tables.py
├── build_paper_tables.py
└── build_paper_figures.py
```

### 5.3 Paper results 层

推荐新增：

```text
results_paper/pead_v1/
├── 00_run_config/
│   ├── paper_model_config.json
│   ├── adopted_loader_manifest.csv
│   ├── well_sample_classification.csv
│   └── run_log.md
├── 01_main_results/
│   ├── all_wells_main_results.csv
│   ├── all_wells_main_results.md
│   ├── hu101/
│   ├── hu102/
│   ├── hu103/
│   └── ht1_001/
├── 02_cbl_validation/
│   ├── cbl_truth_table.csv
│   ├── cbl_error_metrics.csv
│   ├── cbl_comparison.md
│   └── figures/
├── 03_d2dga_ablation/
│   ├── d2dga_on_off_comparison.csv
│   ├── d2dga_on_off_comparison.md
│   └── figures/
├── 04_grid_dt_validation/
│   ├── grid_convergence.csv
│   ├── timestep_sensitivity.csv
│   ├── numerical_reliability.md
│   └── figures/
├── 05_figures_for_paper/
│   ├── fig1_model_schematic.png
│   ├── fig2_coupling_flowchart.png
│   ├── fig3_typical_field_snapshot.png
│   ├── fig4_cbl_comparison.png
│   ├── fig5_d2dga_ablation.png
│   └── fig6_grid_dt_check.png
└── 99_archive/
```

---

## 6. 数据采用值设计

### 6.1 为什么需要 adopted values

现场提取数据、旧 loader 值和论文采用值并不总是一致。例如：

- hu101 的上部复合尾管需要等效几何；
- hu102 的旧 loader 剖面明显过度简化；
- hu103 的真实悬挂器和 loader 中的 top/hanger 语义不同；
- ht1_004 的输入更像优化参数，不是现场实测。

因此，论文流水线需要一张采用值对照表。

### 6.2 推荐文件

```text
参考文档/现场资料提取/01_总表/loader优化采用值对照表.csv
```

字段：

```csv
well_id,field_group,field_name,field_value,loader_value,adopted_value,unit,adopted_source_type,confidence,action,notes
```

### 6.3 数据类型

`adopted_source_type` 建议使用：

```text
field_measured
interpreted
model_assumption
optimized_input
legacy_value
```

定义：

| 类型 | 含义 | 示例 |
|---|---|---|
| field_measured | 现场实测或施工记录 | 泵注体积、密度、井斜测点 |
| interpreted | 测井解释或图像识别 | CBL合格率、CBL分段评价 |
| model_assumption | 模型假设或等效值 | standoff、等效井径 |
| optimized_input | 优化参数化输入 | ht1_004多段优化施工程序 |
| legacy_value | 旧loader遗留值 | 未核验旧参数 |

---

## 7. Paper loader 接口设计

### 7.1 返回值保持兼容

每个 paper loader 应继续返回现有模型需要的四元组：

```python
WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData
```

这样可以复用现有 runner、2D solver、报告和图表工具。

### 7.2 增加论文元数据

paper loader 应额外提供元数据函数，例如：

```python
def get_paper_metadata() -> dict:
    return {
        "well_id": "hu101",
        "paper_version": "pead_v1",
        "sample_class": "A_main_validation",
        "geometry_source": "field_measured + equivalent composite liner approximation",
        "pumping_schedule_source": "field_measured",
        "fluid_source": "field_measured",
        "standoff_source": "model_assumption",
        "cbl_source": "official interpreted CBL report",
    }
```

如果不想改变函数接口，也可以把这些信息输出到 `adopted_input_summary.json`。

---

## 8. 各井 paper loader 设计

### 8.1 hu101

定位：主验证候选井。

设计策略：

- 保留已匹配的现场施工和流体参数；
- 明确复合尾管采用等效几何；
- 结构化写入 CBL合格率 `0.6277`；
- standoff 继续作为模型假设；
- 排除193.7mm回接段资料。

关键说明：

```text
hu101 paper loader 使用 field-measured pumping/fluid/CBL，使用 equivalent composite liner geometry approximation。
```

### 8.2 hu102

定位：现场数据丰富，但旧 loader 需重构后才能进入论文主结果。

设计策略：

- 使用现场井径/井斜剖面替代旧占位剖面；
- 重建完整泵注程序；
- 增加独立 lead slurry、tail slurry、spacer、wash、displacement；
- 拆分 CBL评价窗和地层目标窗；
- 标注缺失的连续 standoff 为模型假设。

### 8.3 hu103

定位：复合尾管建模改进重点井。

设计策略：

- 显式记录真实悬挂器、168.3mm段和139.7mm段；
- 当前求解器若仍单外径，应在 paper loader 中明确等效规则；
- 使用现场 caliper/inclination；
- 泵注程序区分 design schedule 与 actual schedule，论文主结果采用经过确认的一套；
- 补结构化 CBL：`0.1206` 对应 `7338–7712m` 的139.7mm段；
- 168.3mm段 CBL `0.0004` 和整段综合 `0.0605` 单独保存，不混入主验证窗。

### 8.4 ht1_001

定位：待清洗后判断是否纳入验证。

设计策略：

- 不直接使用原始采集CSV；
- 先生成 cleaned schedule；
- 标注几何代理值；
- CBL窗口冻结后再进入 paper loader。

### 8.5 hu2 / HT1-002

定位：暂不作为现场验证井。

设计策略：

- 暂不建立主验证 paper loader；
- 可保留 application loader；
- 元数据标注 `loader_only_without_extracted_field_verification`。

### 8.6 ht1_003

定位：应用或敏感性井。

设计策略：

- 标注固定 `standoff=0.83` 为模型假设；
- 标注井径来源为等效/计算值；
- 不纳入主CBL误差统计。

### 8.7 ht1_004

定位：优化或敏感性案例。

设计策略：

- 标注为 `optimized_input`；
- 不写成现场实测验证；
- 若后续补齐现场日志，再升级为 field-verified loader。

---

## 9. Paper runner 设计

### 9.1 主结果 runner

脚本：

```text
scripts/paper_data/run_paper_main_results.py
```

职责：

1. 从 `paper.registry` 获取论文井清单；
2. 加载 paper loader；
3. 使用统一主模型参数；
4. 运行1D-2D耦合模拟；
5. 输出单井结果和全井汇总表。

主模型参数建议：

```json
{
  "dt": 4.0,
  "nz": 500,
  "ny": 40,
  "enable_d2dga": true,
  "yield_regularization_M": 100,
  "open_outlet": true
}
```

### 9.2 D2DGA消融 runner

脚本：

```text
scripts/paper_data/run_paper_d2dga_ablation.py
```

职责：

- 同一 paper loader；
- 同一网格和时间步；
- 分别运行 `enable_d2dga=True` 和 `False`；
- 输出效率、窜槽、混浆、失稳等对比表。

### 9.3 网格和时间步验证 runner

脚本：

```text
scripts/paper_data/run_paper_grid_dt_validation.py
```

建议代表井：

- hu101
- hu103

网格方案：

```text
nz/ny = 250/20
nz/ny = 500/40
nz/ny = 750/60
```

时间步方案：

```text
dt = 2 s
dt = 4 s
dt = 8 s
```

---

## 10. 单井结果输出设计

每口井主结果目录：

```text
results_paper/pead_v1/01_main_results/<well_id>/
├── run_config.json
├── adopted_input_summary.json
├── time_series.csv
├── depth_profile.csv
├── field_2d.npz
├── summary.json
├── summary.md
├── final_field.png
├── depth_efficiency.png
└── provenance.md
```

### 10.1 run_config.json

记录运行参数：

```json
{
  "paper_version": "pead_v1",
  "loader_version": "paper_loader_v1",
  "well_id": "hu101",
  "dt": 4.0,
  "nz": 500,
  "ny": 40,
  "enable_d2dga": true,
  "yield_regularization_M": 100,
  "open_outlet": true
}
```

### 10.2 adopted_input_summary.json

记录数据来源：

```json
{
  "geometry_source": "field_measured + equivalent composite liner approximation",
  "pumping_schedule_source": "field_measured",
  "fluid_source": "field_measured",
  "caliper_source": "field_measured",
  "inclination_source": "field_measured",
  "standoff_source": "model_assumption",
  "cbl_source": "interpreted official CBL report"
}
```

### 10.3 provenance.md

说明：

- 使用哪个 paper loader；
- 与旧 loader 的主要差异；
- 哪些字段来自现场；
- 哪些字段是模型假设；
- 是否可纳入论文定量验证。

---

## 11. 汇总表设计

### 11.1 全井主结果

文件：

```text
results_paper/pead_v1/01_main_results/all_wells_main_results.csv
```

字段：

```csv
well_id,well_name_cn,sample_class,top_md_m,bottom_md_m,annulus_volume_m3,final_effective_efficiency,final_bulk_cement_fill,cbl_eval_interval_efficiency,target_interval_efficiency,channeling_index,mixing_index,instability_index,summary_json_path,time_series_csv_path,depth_profile_csv_path
```

### 11.2 CBL验证表

文件：

```text
results_paper/pead_v1/02_cbl_validation/cbl_error_metrics.csv
```

字段：

```csv
well_id,well_name_cn,predicted_cbl_efficiency,measured_cbl_pass_rate,absolute_error,relative_error_percent,included_in_metrics,notes
```

### 11.3 D2DGA对比表

文件：

```text
results_paper/pead_v1/03_d2dga_ablation/d2dga_on_off_comparison.csv
```

字段：

```csv
well_id,enable_d2dga,final_effective_efficiency,cbl_eval_interval_efficiency,channeling_index,mixing_index,instability_index,runtime_s
```

### 11.4 网格和时间步表

文件：

```text
results_paper/pead_v1/04_grid_dt_validation/grid_convergence.csv
results_paper/pead_v1/04_grid_dt_validation/timestep_sensitivity.csv
```

字段：

```csv
well_id,nz,ny,dt,final_effective_efficiency,cbl_eval_interval_efficiency,channeling_index,mixing_index,instability_index,runtime_s,reference_case,relative_error_to_reference_percent
```

---

## 12. 错误处理与质量控制

### 12.1 Loader 质量检查

每个 paper loader 必须通过：

1. 能返回 `WellSpec, FluidSpec, PumpingSchedule, ValidationData`；
2. 深度窗口合法；
3. 泵注体积非负；
4. 流体角色完整；
5. CBL窗口不超出模拟段；
6. target window 与 CBL window 不混淆；
7. notes 或 metadata 写明假设。

### 12.2 结果质量检查

每次 paper simulation 必须检查：

1. 是否生成 `run_config.json`；
2. 是否生成 `adopted_input_summary.json`；
3. 是否生成时间序列和深度剖面；
4. 是否有质量守恒或体积异常；
5. 是否有负效率或超过1的异常效率；
6. 是否记录运行耗时和参数。

---

## 13. 测试策略

### 13.1 单元测试

建议新增：

```text
tests/paper_data/test_paper_loaders.py
tests/paper_data/test_paper_result_schema.py
tests/paper_data/test_cbl_metrics.py
```

测试内容：

- paper loader 返回对象类型正确；
- sample class 正确；
- CBL窗口在模拟段内；
- pumping schedule 总量大于0；
- CBL误差计算公式正确。

### 13.2 Smoke simulation

正式高分辨率运行前，先用低分辨率跑通：

```text
dt = 8 s
nz = 120
ny = 24
```

目标：

- 不崩溃；
- 输出字段完整；
- 效率范围合理；
- 质量守恒无明显异常。

---

## 14. 实施阶段建议

### 阶段1： adopted values 对照表

输出：

```text
参考文档/现场资料提取/01_总表/loader优化采用值对照表.csv
```

优先井：

- hu101
- hu102
- hu103

### 阶段2： paper loader 最小实现

优先实现：

- `paper_hu101_loader.py`
- `paper_hu102_loader.py`
- `paper_hu103_loader.py`

### 阶段3： paper main runner

实现：

- `scripts/paper_data/run_paper_main_results.py`

先用 smoke 参数跑通。

### 阶段4： CBL验证与主结果表

实现：

- `build_cbl_validation_tables.py`
- `build_paper_tables.py`

### 阶段5： D2DGA和网格验证

实现：

- `run_paper_d2dga_ablation.py`
- `run_paper_grid_dt_validation.py`

### 阶段6：正式论文主参数运行

使用：

```text
dt = 4 s
nz = 500
ny = 40
```

---

## 15. 设计决策

推荐采用：

> 新建 paper loader 层 + 新建 paper runner + 新建 results_paper 结果目录。

不推荐直接修改旧 loader，原因是：

1. 旧 loader 已被历史 runner 和 results 依赖；
2. 直接修改会破坏历史可追溯性；
3. 论文需要更严格的数据来源标注；
4. 新旧结果必须隔离，才能避免论文表格和图件混用口径。

---

## 16. 待确认事项

进入实施前，需要确认：

1. 第一批 paper loader 是否只做 hu101、hu102、hu103；
2. ht1_001 是否先清洗数据再进入 paper loader；
3. ht1_003、ht1_004 是否只作为敏感性/优化案例；
4. paper results 根目录是否采用 `results_paper/pead_v1/`；
5. 论文主模型参数是否固定为 `dt=4, nz=500, ny=40`。
