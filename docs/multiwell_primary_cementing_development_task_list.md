# 多井通用固井顶替效率模型开发任务清单

## 1. 文档目的

本文档将《多井通用固井顶替效率模型技术方案》进一步拆解为可执行的开发任务清单，重点回答：

1. 先改哪些文件；
2. 先建哪些数据结构；
3. 先实现哪个模块；
4. 每个阶段完成后的验收标准是什么。

本清单面向当前项目代码结构，默认以 `cemdisp/` 为主开发目录，旧的单井脚本视为可迁移的参考资产，而不是长期核心。

---

## 2. 总体开发策略

采用“**先统一输入，再泛化环空核心，再接入套管内输运，最后做解释层和批量化**”的路线。

### 总体原则

1. **先重构，再增强**：先把单井硬编码拆掉，再新增上游建模能力；
2. **保留现有有效资产**：Hu102 的二维环空模型逻辑尽量迁移而非废弃；
3. **每一阶段都可运行**：不能做成只有最终阶段才能验证的大重构；
4. **多井复用优先于单井调参**：所有新增能力都应优先做成“可配置”的。

---

## 3. 第一批必须修改/新增的文件

## 3.1 第一优先级：数据结构与输入层

### 需要新增的文件

| 文件 | 作用 | 优先级 |
|---|---|---|
| `cemdisp/data/well_spec.py` | 定义井筒几何、井段窗口、鞋口/悬挂器等统一数据结构 | P0 |
| `cemdisp/data/fluid_spec.py` | 定义流体物性统一数据结构 | P0 |
| `cemdisp/data/pumping_schedule.py` | 定义施工程序、分段注入、停泵/续泵结构 | P0 |
| `cemdisp/data/validation_data.py` | 定义 CBL、施工日报等校验数据结构 | P1 |
| `cemdisp/data/loaders/hu102_loader.py` | 把 Hu102 现有资料装配成标准数据包 | P0 |
| `cemdisp/data/loaders/__init__.py` | 统一数据加载接口 | P0 |

### 需要调整的文件

| 文件 | 调整内容 |
|---|---|
| `cemdisp/geometry.py` | 从“当前几何工具”提升为与 `WellSpec` 配套使用的通用几何层 |
| `cemdisp/fluid.py` | 与新 `FluidSpec` 对齐，区分输入规范与求解内部表示 |

---

## 3.2 第二优先级：环空二维核心泛化

### 需要新增的文件

| 文件 | 作用 | 优先级 |
|---|---|---|
| `cemdisp/models2d/annulus_d2dga.py` | 多井通用的环空二维顶替核心求解器 | P0 |
| `cemdisp/models2d/annulus_state.py` | 定义二维状态对象（浓度场、泥饼场、速度场） | P0 |
| `cemdisp/models2d/evaluation_windows.py` | 统一处理 CBL段、目标层段等窗口计算 | P1 |

### 需要参考迁移的旧文件

| 来源文件 | 迁移内容 |
|---|---|
| `hu102model/hu102_tail_d2dga_model.py` | `build_geom()`、`compute_props()`、`compute_velocity()`、`simulate()`、效率计算逻辑 |
| `hu101model/hu101_d2dga_model.py` | 图表与多井输出风格参考 |

### 需要剥离的硬编码

重点去掉：

- 井号硬编码；
- 井深常量；
- 固定井段掩码；
- `FIELD_RATE_M3_MIN` 之类施工量常量；
- `TRACKED = ["balance", "spacer", "lead", "tail"]` 这种固定流体序列限制。

---

## 3.3 第三优先级：套管内 1D 输运层

### 需要新增的文件

| 文件 | 作用 | 优先级 |
|---|---|---|
| `cemdisp/transport1d/casing_flow.py` | 套管内 1D 输运主模块 | P0 |
| `cemdisp/transport1d/interface_tracking.py` | 前缘/段塞追踪 | P0 |
| `cemdisp/transport1d/pipe_exit_state.py` | 生成鞋口处实时出流状态 | P0 |
| `cemdisp/transport1d/dispersion.py` | 套管内对流-弥散扩展（后期） | P2 |

### 模块顺序建议

先实现：

1. `interface_tracking.py`
2. `casing_flow.py`
3. `pipe_exit_state.py`

最后再接：

4. `dispersion.py`

原因是先把“从地面开泵到鞋口边界”的**基础容积时序**打通，比先做复杂混合更重要。

---

## 3.4 第四优先级：1D → 2D 边界桥接

### 需要新增的文件

| 文件 | 作用 | 优先级 |
|---|---|---|
| `cemdisp/models2d/boundary_bridge.py` | 把鞋口出流转成环空二维入口边界 | P0 |

### 替代目标

该模块最终要替代当前 Hu102 脚本里 `boundary_state()` 的硬编码逻辑。

桥接接口建议统一为：

```python
get_annulus_inlet_state(t) -> {
    "vec": ...,
    "q_m3s": ...,
    "stage": ...,
}
```

这样现有二维求解器只需依赖统一入口接口，不直接感知套管内模块细节。

---

## 3.5 第五优先级：诊断与解释层

### 需要新增的文件

| 文件 | 作用 | 优先级 |
|---|---|---|
| `cemdisp/diagnostics/efficiency.py` | 实时与最终效率计算 | P0 |
| `cemdisp/diagnostics/channeling.py` | 窜槽指标 | P1 |
| `cemdisp/diagnostics/mixing.py` | 混浆指标 | P1 |
| `cemdisp/diagnostics/instability.py` | 失稳指标 | P1 |
| `cemdisp/diagnostics/quality_proxy.py` | 质量响应效率映射 | P1 |

### 需要拆分的旧逻辑

当前这些内容散落在 `hu102_tail_d2dga_model.py` 的 `simulate()` 和 `save_outputs()` 中，后续应拆成单独模块，避免主求解器同时承担“算流场 + 做评价 + 画图”三类职责。

---

## 3.6 第六优先级：验证与输出层

### 需要新增的文件

| 文件 | 作用 | 优先级 |
|---|---|---|
| `cemdisp/validation/mass_balance.py` | 质量守恒验证 | P0 |
| `cemdisp/validation/field_compare.py` | 与现场摘要数据对比 | P1 |
| `cemdisp/validation/cbl_compare.py` | 与 CBL 剖面对比 | P1 |
| `cemdisp/reporting/summary.py` | 汇总 JSON / Markdown 摘要 | P1 |
| `cemdisp/reporting/plots.py` | 统一图表输出 | P1 |
| `cemdisp/reporting/tables.py` | 统一 CSV / 表格输出 | P1 |

---

## 4. 第一批必须先建的数据结构

开发顺序上，数据结构必须先于求解器泛化。

## 4.1 `WellSpec`

必须先建，原因：

- 当前 Hu102 和 Hu101 的几何常量都写死在脚本里；
- 不先抽象井筒结构，后面的“多井复用”无法成立。

### 最小字段建议

- `well_name`
- `top_md_m`
- `bottom_md_m`
- `shoe_md_m`
- `hanger_md_m`
- `casing_id_mm`
- `liner_od_mm`
- `liner_id_mm`
- `hole_profile`
- `inclination_profile`
- `standoff_profile`
- `evaluation_windows`

## 4.2 `FluidSpec`

必须先建，原因：

- 当前不同井脚本把流体参数定义方式混在求解器里；
- 后续若支持更多井型和更多流体段，必须配置化。

### 最小字段建议

- `name`
- `role`
- `density`
- `rheology_model`
- `plastic_viscosity`
- `yield_stress`
- `power_law_n`
- `consistency_k`

## 4.3 `PumpingScheduleStep`

必须先建，原因：

- 模型要从地面开泵开始；
- 没有标准施工步骤对象，就无法建立套管内输运层。

### 最小字段建议

- `step_name`
- `fluid_name`
- `volume_m3`
- `rate_m3_min`
- `start_time_s`
- `end_time_s`
- `remarks`

## 4.4 `EvaluationWindow`

必须先建，原因：

- 当前 Hu102 / Hu101 对井段评价掩码都写死；
- 后续多井需要灵活切换 CBL段、目标层段、异常段。

### 最小字段建议

- `name`
- `top_md_m`
- `bottom_md_m`
- `window_type`

---

## 5. 先实现哪个模块

### 结论

**第一个真正要实现的模块，不是套管内 1D 输运，而是“多井通用环空二维核心”的参数化版本。**

原因：

1. 当前最成熟的已有资产在环空层；
2. 如果环空层还是单井脚本，就算加了套管内模块，也没法多井复用；
3. 套管内 1D 层只是入口边界的上游增强，而环空 2D 层才是最终效率计算的核心。

### 推荐顺序

1. **先实现通用数据结构**
2. **再实现通用环空 2D 核心**
3. **再实现套管内 1D 输运层**
4. **再实现边界桥接层**
5. **最后补质量解释层和报告层**

---

## 6. 分阶段开发任务清单

## Phase 0：文档与输入口径统一

### 任务 0.1

建立标准输入规范文档。

**完成标准**：

- 明确井数据、流体数据、施工程序、校验数据四类输入结构；
- 给出 Hu102 样例字段说明。

### 任务 0.2

创建 `cemdisp/data/` 目录及基础数据类文件。

**涉及文件**：

- `cemdisp/data/well_spec.py`
- `cemdisp/data/fluid_spec.py`
- `cemdisp/data/pumping_schedule.py`
- `cemdisp/data/__init__.py`

---

## Phase 1：Hu102 环空核心剥离与泛化

### 任务 1.1

把 `hu102_tail_d2dga_model.py` 中的几何、流体、效率评价逻辑按职责拆分。

**优先迁移函数**：

- `build_geom()`
- `compute_props()`
- `compute_velocity()`
- `interp2()`
- `simulate()` 中与求解有关的部分

### 任务 1.2

在 `cemdisp/models2d/annulus_d2dga.py` 中建立通用求解器接口。

**完成标准**：

- 不再依赖 Hu102 文件名或固定井段；
- 由 `WellSpec`、`FluidSpec`、`EvaluationWindow` 输入驱动。

### 任务 1.3

将 Hu102 作为第一口标准样例井接入新接口。

**涉及文件**：

- `cemdisp/data/loaders/hu102_loader.py`
- `cemdisp/models2d/annulus_d2dga.py`

**完成标准**：

- 新接口对 Hu102 跑出的核心效率曲线与旧脚本趋势一致；
- 质量守恒误差可接受。

---

## Phase 2：套管内 1D 输运 MVP

### 任务 2.1

实现 `PumpingSchedule` 到“套管内前缘传播”的基础容积追踪。

**涉及文件**：

- `cemdisp/transport1d/interface_tracking.py`

### 任务 2.2

实现套管内 1D 输运主模块。

**涉及文件**：

- `cemdisp/transport1d/casing_flow.py`

**完成标准**：

- 输入地面程序；
- 输出任意时刻各流体前缘在套管内的位置。

### 任务 2.3

实现鞋口出流状态计算。

**涉及文件**：

- `cemdisp/transport1d/pipe_exit_state.py`

**完成标准**：

- 给出 `vec(t)`、`q(t)`、`stage(t)`；
- 可用于替代 Hu102 现有 `boundary_state()`。

---

## Phase 3：1D → 2D 耦合

### 任务 3.1

建立桥接层，将套管内输出接入环空模型。

**涉及文件**：

- `cemdisp/models2d/boundary_bridge.py`

### 任务 3.2

修改二维求解器，使其不再依赖脚本内硬编码边界函数，而是依赖统一边界接口。

**涉及文件**：

- `cemdisp/models2d/annulus_d2dga.py`

### 任务 3.3

用 Hu102 跑通“地面开泵 → 套管内输运 → 环空顶替”的第一版全流程模拟。

**完成标准**：

- 能输出鞋口出流时序；
- 能输出实时顶替效率与最终顶替效率；
- 结果曲线可与当前环空入口起算版对比。

---

## Phase 4：诊断与输出重构

### 任务 4.1

拆分效率与质量解释逻辑。

**涉及文件**：

- `cemdisp/diagnostics/efficiency.py`
- `cemdisp/diagnostics/quality_proxy.py`

### 任务 4.2

统一输出模块。

**涉及文件**：

- `cemdisp/reporting/summary.py`
- `cemdisp/reporting/plots.py`
- `cemdisp/reporting/tables.py`

### 任务 4.3

重新生成 Hu102 全部图表和报告，验证新架构未丢失现有可交付物能力。

---

## Phase 5：第二口井接入与多井复用验证

### 任务 5.1

选一口第二口井，整理输入数据包。

### 任务 5.2

编写对应 loader。

**涉及文件**：

- `cemdisp/data/loaders/<well_name>_loader.py`

### 任务 5.3

验证无需修改求解器，只通过替换输入数据即可运行。

**完成标准**：

- 新井可运行；
- 输出结构一致；
- 评价井段可配置。

---

## Phase 6：套管内混合扩展（条件触发）

### 触发条件

当以下任一情况成立时启动本阶段：

- 鞋口纯流体边界与现场数据偏差明显；
- 前置液/隔离液体积较小且界面混合显著；
- 施工存在频繁停泵/续泵；
- 长套管内输运导致明显界面退化。

### 任务 6.1

实现 1D 对流-弥散模块。

**涉及文件**：

- `cemdisp/transport1d/dispersion.py`

### 任务 6.2

比较“无混合”和“有混合”两套鞋口边界对环空最终效率的影响。

---

## 7. 第一版测试清单

## 7.1 必须先建的测试文件

| 文件 | 测试内容 |
|---|---|
| `tests/test_well_spec.py` | 井数据结构初始化与字段合法性 |
| `tests/test_pumping_schedule.py` | 施工程序排序、体积、时长一致性 |
| `tests/test_casing_flow.py` | 套管内前缘传播与体积守恒 |
| `tests/test_pipe_exit_state.py` | 鞋口边界生成正确性 |
| `tests/test_annulus_d2dga.py` | 环空二维质量守恒与基本趋势 |
| `tests/test_boundary_bridge.py` | 1D→2D 耦合接口正确性 |

## 7.2 第一批验收项

1. Hu102 在新架构下可跑通；
2. 不改井数据结构即可换第二口井；
3. 可输出从地面开泵开始的实时效率曲线；
4. 可输出鞋口出流历史；
5. 可输出最终顶替效率和评价井段效率；
6. 质量守恒检查通过；
7. 图表与摘要文件能自动生成。

---

## 8. 开发顺序总结

如果只看“现在开始第一周应该干什么”，建议顺序如下：

### 第一步

先建输入数据结构：

- `well_spec.py`
- `fluid_spec.py`
- `pumping_schedule.py`

### 第二步

把 Hu102 现有二维环空模型迁到：

- `cemdisp/models2d/annulus_d2dga.py`

### 第三步

实现 Hu102 标准 loader：

- `hu102_loader.py`

### 第四步

实现套管内 1D 前缘追踪：

- `interface_tracking.py`
- `casing_flow.py`
- `pipe_exit_state.py`

### 第五步

实现边界桥接层：

- `boundary_bridge.py`

### 第六步

拆出诊断与报告模块。

---

## 9. 最终建议

后续开发不要再沿着“再做一份 Hu103、再做一份 Hu104 脚本”的路径继续走。正确方向是：

> **先把 Hu102 现有二维环空模型泛化成通用核心，再加地面开泵到鞋口的 1D 输运模块，最后把输入数据包和评价窗口全部配置化。**

这样做的好处是：

1. 现有 Hu102 资产保住了；
2. 地面真实施工流程接进来了；
3. 后续换井不需要复制模型；
4. 模型复杂度是渐进增加的，而不是一开始失控。

这就是当前项目最稳妥、最可落地的开发路线。
