# cement model 所需现场井尾管固井数据清单与提取规划

## 0. 文档目的

本文档只说明 `cement model` 固井顶替模型及《石油勘探与开发》论文计算所需要的**现场井数据类型、字段、优先级和提取口径**。

本文档不要求执行者逐个分析某个目录中的具体文件名，也不要求在文档中绑定具体文件。执行AI或人工助手只需根据本文档列出的数据需求，到现场资料中查找、提取、整理对应数据。

核心目标：

> 把现场尾管固井资料整理成 `cement model` 可用的标准输入数据，并为论文中的井资料完整性审查、模型计算、CBL验证和结果分析提供依据。

---

## 1. 模型需要的现场数据总览

`cement model` 当前是一个面向尾管固井顶替效率评价的 1D—2D耦合模型。模型需要的现场数据可分为8类：

| 类别 | 用途 | 是否核心必需 |
|---|---|---|
| A. 井基本信息 | 确认井号、井名、区块、资料对应关系 | 必需 |
| B. 井身结构与管柱数据 | 建立套管/尾管/环空几何 | 必需 |
| C. 井径数据 | 建立真实环空间隙、井径扩大率 | 高优先级 |
| D. 井斜/轨迹数据 | 考虑井斜、重力分异、模拟段空间位置 | 高优先级 |
| E. 居中度/扶正器/偏心数据 | 建立偏心环空宽边/窄边差异 | 高优先级 |
| F. 固井施工泵注程序 | 驱动1D套管输运和2D环空入口边界 | 必需 |
| G. 流体性能数据 | 描述钻井液、前置液、隔离液、水泥浆、顶替液 | 必需 |
| H. CBL/VDL固井质量评价数据 | 论文后验验证，不参与反标定 | 论文验证必需 |

此外，若现场资料中存在目标层段、油气水层、温度、压力、施工异常等信息，可作为论文解释和讨论的辅助数据。

---

## 2. 当前目录表反映出的资料类型规划

根据现场资料目录表截图，现场资料大致包含以下类型：

1. 固井设计书、固井施工设计；
2. 水泥浆设计、实验报告、浆体性能资料；
3. 套管柱结构、尾管结构、管串数据；
4. 井斜、井径、测井曲线资料；
5. CBL/VDL固井质量评价资料；
6. 油气水层、目的层、解释层段资料；
7. 施工总结、施工记录、施工参数表；
8. 图片、扫描件、图片型PDF、Word内嵌图；
9. 各井历史数据表、截图、图像曲线。

因此提取规划应按“模型字段需求”组织，而不是按目录文件名组织。执行者看到现场资料后，应优先寻找能填充下列标准数据表的内容。

---

## 3. 待覆盖现场井

优先覆盖以下现场井。若现场资料中发现其他相关井，可补充。

| well_id | 现场井名 | 模型/论文用途 |
|---|---|---|
| hu101 | 呼101 | 主验证候选井 / 多井计算 |
| hu102 | 呼102 | 主验证候选井 / 多井计算 |
| hu103 | 呼103 | 主验证候选井 / 多井计算 |
| hu1 | 呼探1 | 多井计算 / 验证候选 |
| ht1_001 | 呼探1-001 / HT1-001 | 多井计算 / 验证候选 |
| hu2 | 呼探1-002 / HT1-002 | 多井计算 / 验证候选 |
| ht1_003 | 呼1-003 / HT1-003 | 敏感性分析 / 多井计算 |
| ht1_004 | 呼1-004 / HT1-004 | 敏感性分析 / 多井计算 |

---

## 4. A类数据：井基本信息

### 4.1 用途

用于确认资料归属，避免井名、井号、别名混淆。

### 4.2 必需字段

建议输出到：

```text
well_basic_info.csv
```

字段：

```csv
well_id,well_name_cn,alias,block,field,well_type,cementing_type,data_source,confidence,notes
```

### 4.3 需要提取的数据

| 字段 | 含义 | 说明 |
|---|---|---|
| well_id | 模型内部井ID | 如 hu101、ht1_004 |
| well_name_cn | 中文井名 | 如 呼101、呼1-004 |
| alias | 别名 | 如 HT1-004、呼探1-002 |
| block | 区块 | 若资料有则提取 |
| field | 气田/油田 | 若资料有则提取 |
| well_type | 井型 | 直井、定向井、水平井等 |
| cementing_type | 固井类型 | 尾管固井、技术套管固井等 |

---

## 5. B类数据：井身结构与管柱数据

### 5.1 用途

这是模型建立几何结构的核心输入，决定套管内1D输运体积、环空体积、模拟段顶底深、鞋口位置等。

### 5.2 建议输出文件

```text
well_geometry.csv
casing_liner_string.csv
```

### 5.3 `well_geometry.csv`字段

```csv
well_id,item,value,unit,top_md_m,bottom_md_m,data_type,source_description,confidence,notes
```

### 5.4 必须提取的井身结构字段

| item | 含义 | 模型用途 | 优先级 |
|---|---|---|---|
| simulation_top_md | 模拟段顶深 | 2D环空模拟上边界 | 必需 |
| simulation_bottom_md | 模拟段底深 | 2D环空模拟下边界 | 必需 |
| openhole_top_md | 裸眼段顶深 | 确定裸眼环空段 | 高 |
| openhole_bottom_md | 裸眼段底深 | 确定裸眼环空段 | 高 |
| liner_hanger_md | 尾管悬挂器深度 | 环空顶端/尾管重叠段判断 | 必需 |
| liner_shoe_md | 尾管鞋深度 | 1D到2D耦合入口 | 必需 |
| casing_inner_diameter | 上部套管内径 | 1D输运截面积 | 必需 |
| casing_outer_diameter | 上部套管外径 | 环空几何参考 | 中 |
| liner_inner_diameter | 尾管内径 | 1D尾管输运截面积 | 必需 |
| liner_outer_diameter | 尾管外径 | 环空间隙计算 | 必需 |
| borehole_nominal_diameter | 名义井眼直径 | 环空基准尺寸 | 必需 |
| drill_bit_size | 钻头尺寸 | 校核名义井径 | 高 |
| float_collar_md | 浮箍深度 | 施工解释/附件位置 | 中 |
| float_shoe_md | 浮鞋深度 | 鞋口位置校核 | 中 |

### 5.5 `casing_liner_string.csv`字段

```csv
well_id,component,top_md_m,bottom_md_m,outer_diameter_mm,inner_diameter_mm,weight_or_grade,quantity,spacing_m,data_type,source_description,confidence,notes
```

### 5.6 需要提取的管柱项目

| component | 含义 |
|---|---|
| surface_casing | 表层套管，如有 |
| intermediate_casing | 技术套管，如有 |
| production_casing | 生产套管，如有 |
| liner | 尾管 |
| liner_hanger | 尾管悬挂器 |
| float_collar | 浮箍 |
| float_shoe | 浮鞋 |
| centralizer | 扶正器 |
| plug | 胶塞/碰压塞 |
| other_accessory | 其他附件 |

---

## 6. C类数据：井径数据

### 6.1 用途

井径数据用于计算真实环空间隙和井径扩大率，是偏心环空2D顶替模型的重要输入。

### 6.2 建议输出文件

```text
caliper_profile.csv
```

字段：

```csv
well_id,md_m,caliper_mm,borehole_nominal_diameter_mm,enlargement_ratio,data_type,source_description,confidence,notes
```

### 6.3 提取要求

| 情况 | 处理方式 |
|---|---|
| 有连续井径曲线数据 | 提取深度—井径序列 |
| 有测井图像但无表格 | 进行曲线数字化，标记为 interpreted |
| 只有分段平均井径 | 提取分段顶底深和平均井径 |
| 只有钻头尺寸 | 只能作为名义井径，不能伪造成实际井径曲线 |
| 无井径资料 | 标记 missing |

### 6.4 模型最小需求

最低限度需要：

- 模拟段名义井径；
- 或模拟段平均井径；
- 或可用于估计环空体积的井径数据。

若无实际井径曲线，后续模型只能使用名义井径或假设井径，论文中必须标记该井资料质量较低。

---

## 7. D类数据：井斜/井轨迹数据

### 7.1 用途

井斜影响重力分异、偏心趋势和宽窄边顶替特征。

### 7.2 建议输出文件

```text
inclination_profile.csv
```

字段：

```csv
well_id,md_m,inclination_deg,azimuth_deg,tvd_m,data_type,source_description,confidence,notes
```

### 7.3 提取要求

| 数据 | 是否必需 | 说明 |
|---|---|---|
| MD测深 | 必需 | 井斜测点深度 |
| 井斜角 | 必需 | 模型主要使用 |
| 方位角 | 可选 | 若模型后续拓展可用 |
| TVD | 可选 | 用于深度校核 |

若只有“最大井斜”“井斜范围”，应记录为摘要信息，不得伪造连续测斜表。

---

## 8. E类数据：居中度、扶正器和偏心数据

### 8.1 用途

用于建立偏心环空几何，区分宽边、窄边和中间区域。

### 8.2 建议输出文件

```text
centralization_profile.csv
casing_liner_string.csv
```

### 8.3 `centralization_profile.csv`字段

```csv
well_id,md_m,standoff_ratio,centralization_percent,eccentricity_ratio,data_type,source_description,confidence,notes
```

### 8.4 需要提取的数据

| 数据 | 含义 | 说明 |
|---|---|---|
| 扶正器位置 | 每个扶正器深度 | 可用于估算居中度 |
| 扶正器数量 | 总数或分段数量 | 用于判断偏心风险 |
| 扶正器间距 | 相邻扶正器距离 | 用于估算居中度 |
| 扶正器型号 | 刚性/弹性/螺旋等 | 用于论文解释 |
| 居中度 | 现场计算或设计值 | 若有则优先提取 |
| 偏心度 | 模型直接可用 | 若有则优先提取 |
| standoff | 离壁比 | 若有则优先提取 |

### 8.5 注意事项

- 现场资料中很多井可能没有连续居中度曲线；
- 若只有扶正器设计，不能直接写成连续居中度；
- 可以记录“具备扶正器信息，但需模型估算居中度”；
- 估算值必须标记为 `model_assumption`，不能标为现场实测。

---

## 9. F类数据：固井施工泵注程序

### 9.1 用途

这是1D套管输运模型的直接驱动输入。模型需要知道各流体段的泵入顺序、体积、排量和密度。

### 9.2 建议输出文件

```text
pumping_schedule.csv
```

字段：

```csv
well_id,step_index,stage_name,fluid_role,fluid_name,volume_m3,rate_m3_min,start_time_min,end_time_min,density_g_cm3,event_tag,data_type,source_description,confidence,notes
```

### 9.3 必须提取的施工阶段

| fluid_role | 说明 | 是否必需 |
|---|---|---|
| mud | 井内原钻井液 | 高 |
| wash | 冲洗液/前置液 | 有则提取 |
| spacer | 隔离液 | 有则提取 |
| lead_cement | 领浆 | 有则提取 |
| intermediate_cement | 中间浆 | 有则提取 |
| tail_cement | 尾浆 | 必需，若有水泥浆 |
| displacement | 顶替液 | 必需 |
| plug | 胶塞/碰压 | 有则提取 |
| shutdown | 停泵/候凝事件 | 有则提取 |
| rate_change | 变排量事件 | 有则提取 |

### 9.4 每个施工阶段至少需要

| 字段 | 是否必需 | 说明 |
|---|---|---|
| 阶段顺序 | 必需 | step_index |
| 流体名称 | 必需 | 如隔离液、尾浆、清水等 |
| 流体角色 | 必需 | 统一到fluid_role |
| 体积 | 必需 | m³ |
| 排量 | 高 | m³/min |
| 密度 | 高 | g/cm³ |
| 开始/结束时间 | 有则提取 | 可用于校核排量 |
| 碰压/停泵事件 | 有则提取 | 用于解释施工异常 |

### 9.5 特别注意

- 现场施工设计和实际施工记录可能不同，应优先记录实际施工值；
- 若设计值和实际值均存在，应分别记录，并在notes中说明；
- 若只有体积和总时间，可计算平均排量，但必须说明是计算值；
- 顶替液末段补量、返出、漏失、停泵等异常必须记录。

---

## 10. G类数据：流体性能数据

### 10.1 用途

流体性能用于计算非牛顿流体表观黏度、流度比、顶替稳定性和D2DGA通量修正。

### 10.2 建议输出文件

```text
fluid_properties.csv
```

字段：

```csv
well_id,fluid_name,fluid_role,density_g_cm3,rheology_model,pv_mpa_s,yp_pa,yield_stress_pa,consistency_k_pa_sn,flow_index_n,viscosity_mpa_s,test_temperature_c,data_type,source_description,confidence,notes
```

### 10.3 需要提取的流体

| fluid_role | 流体 | 是否必需 |
|---|---|---|
| mud | 钻井液 | 必需 |
| wash | 冲洗液/前置液 | 有则提取 |
| spacer | 隔离液 | 有则提取 |
| lead_cement | 领浆 | 有则提取 |
| intermediate_cement | 中间浆 | 有则提取 |
| tail_cement | 尾浆 | 必需 |
| displacement | 顶替液 | 必需 |

### 10.4 需要提取的参数

| 参数 | 含义 | 优先级 |
|---|---|---|
| density_g_cm3 | 密度 | 必需 |
| rheology_model | 流变模型 | 高 |
| pv_mpa_s | 塑性黏度 | 高 |
| yp_pa | 动切力 | 高 |
| yield_stress_pa | 屈服应力 | 高 |
| consistency_k_pa_sn | 稠度系数K | 中 |
| flow_index_n | 流性指数n | 中 |
| viscosity_mpa_s | 牛顿黏度或表观黏度 | 中 |
| test_temperature_c | 测试温度 | 中 |

### 10.5 六速旋转黏度计读数

若现场资料中提供六速旋转黏度计读数，应单独保留原始读数，建议附加表：

```text
rheometer_readings.csv
```

字段：

```csv
well_id,fluid_name,temperature_c,theta_600,theta_300,theta_200,theta_100,theta_6,theta_3,source_description,confidence,notes
```

然后再根据需要换算PV、YP、n、K等参数。

### 10.6 特别注意

- 不得把 `yield_regularization_M` 当作现场流体参数；
- `yield_regularization_M` 是数值正则化参数，只属于模型设置；
- 水泥浆密度、钻井液密度、隔离液密度必须区分；
- 领浆和尾浆不能混写；
- 设计流变参数和现场复测流变参数如有差异，应分别记录。

---

## 11. H类数据：CBL/VDL固井质量评价数据

### 11.1 用途

CBL/VDL数据用于论文后验验证，即判断模型预测的水力顶替效率与现场固井质量评价是否一致。

CBL不参与模型反标定，不等同于顶替效率。

### 11.2 建议输出文件

```text
cbl_evaluation.csv
```

字段：

```csv
well_id,well_name_cn,cbl_top_md_m,cbl_bottom_md_m,cbl_pass_rate,cbl_quality_class,interpretation_summary,include_in_validation,data_type,source_description,confidence,notes
```

### 11.3 需要提取的数据

| 字段 | 含义 | 是否必需 |
|---|---|---|
| cbl_top_md_m | CBL评价段顶深 | 验证必需 |
| cbl_bottom_md_m | CBL评价段底深 | 验证必需 |
| cbl_pass_rate | CBL合格率 | 定量验证必需 |
| cbl_quality_class | 固井质量等级 | 有则提取 |
| interpretation_summary | 解释结论 | 有则提取 |
| 差胶结井段 | 差胶结区间 | 有则提取 |
| 优质胶结井段 | 优质区间 | 有则提取 |

### 11.4 处理原则

| 情况 | 处理方式 |
|---|---|
| 有明确合格率 | 可用于定量验证候选 |
| 只有分段质量等级 | 用于定性验证 |
| 只有CBL/VDL图像 | 提取评价段，必要时低置信度数字化 |
| 无CBL资料 | 不纳入CBL定量误差统计 |

---

## 12. 辅助数据：目标层段、油气水层、温压和施工异常

这些数据不是所有模型计算都必需，但对论文解释很有价值。

### 12.1 目标层段 `target_intervals.csv`

字段：

```csv
well_id,interval_name,top_md_m,bottom_md_m,purpose,source_description,data_type,confidence,notes
```

用途：

- 统计目标层段顶替效率；
- 解释CBL评价重点井段；
- 论文中展示目标区间风险。

### 12.2 温度压力数据 `temperature_pressure_profile.csv`

字段：

```csv
well_id,md_m,temperature_c,pressure_mpa,data_type,source_description,confidence,notes
```

用途：

- 解释流变参数适用性；
- 后续可扩展温度耦合模型；
- 当前主模型可作为背景信息。

### 12.3 施工异常 `construction_events.csv`

字段：

```csv
well_id,event_time_min,event_md_m,event_type,description,related_stage,source_description,confidence,notes
```

常见异常：

- 漏失；
- 返出异常；
- 停泵；
- 变排量；
- 碰压异常；
- 替量不足；
- 施工压力异常。

---

## 13. 图片、扫描PDF和Word内嵌图片中的数据要求

现场资料目录表显示现场资料可能包含大量图片、扫描件或截图。执行AI必须检查这些图像资料。

### 13.1 必须从图片中识别的数据类型

| 图片内容 | 应提取的数据 |
|---|---|
| 井身结构图 | 套管/尾管尺寸、鞋深、悬挂器深度、井段 |
| 管柱结构图 | 尾管、浮箍、浮鞋、扶正器、附件位置 |
| 施工程序截图 | 泵注顺序、体积、排量、密度、时间 |
| 水泥浆实验表截图 | 密度、PV、YP、稠化时间、流变读数 |
| 井径曲线 | 深度—井径数据或分段井径 |
| 井斜曲线 | 深度—井斜角数据 |
| CBL/VDL图 | 评价段、差胶结段、合格率或质量等级 |
| 目标层图 | 目标层段顶底深 |

### 13.2 图片数据记录要求

每口井应建立：

```text
image_extraction_index.csv
```

字段：

```csv
well_id,image_id,source_description,image_type,detected_content,contains_useful_data,extracted_to_file,ocr_or_vision_method,confidence,needs_human_review,notes
```

### 13.3 曲线数字化要求

若只有曲线图，应记录：

```csv
well_id,curve_type,x_axis_name,x_axis_unit,y_axis_name,y_axis_unit,digitization_method,estimated_error,confidence,notes
```

原则：

- 曲线读数一般标记为 `interpreted`；
- 置信度通常为 `medium` 或 `low`；
- 不得把曲线粗读数写成高精度实测数据；
- 坐标轴不清的曲线只能列入人工复核清单。

---

## 14. 数据完整性分级

提取完成后，应按井给出完整性初评。

| 类别 | 条件 | 论文用途 |
|---|---|---|
| A_main_validation | 几何、施工、流体、CBL定量资料均完整 | 主验证井 |
| B_application | 几何、施工、流体基本完整，CBL不足 | 多井应用计算 |
| C_sensitivity | 输入较完整，适合参数扰动 | 敏感性分析 |
| D_pending | 关键资料缺失或图片/OCR未复核 | 暂不纳入主论文 |

建议输出：

```text
well_data_completeness.csv
```

字段：

```csv
well_id,well_name_cn,has_geometry,has_liner_string,has_caliper,has_inclination,has_centralization_or_centralizer,has_pumping_schedule,has_fluid_density,has_fluid_rheology,has_cbl_interval,has_cbl_pass_rate,has_target_interval,has_image_checked,suggested_class,missing_key_data,notes
```

---

## 15. 单井标准交付文件

每口井建议最终形成：

```text
well_basic_info.csv
well_geometry.csv
casing_liner_string.csv
caliper_profile.csv
inclination_profile.csv
centralization_profile.csv
pumping_schedule.csv
fluid_properties.csv
rheometer_readings.csv
cbl_evaluation.csv
target_intervals.csv
temperature_pressure_profile.csv
construction_events.csv
image_extraction_index.csv
source_trace.csv
extraction_notes.md
```

如果某类数据不存在，仍应保留文件或在说明中明确写明 `missing`。

---

## 16. 最小可用数据集定义

### 16.1 能进入模型主计算的最低要求

一口井至少需要：

1. 模拟段顶底深；
2. 尾管鞋深；
3. 尾管外径/内径；
4. 井眼直径或井径估计；
5. 泵注程序：体积 + 排量 + 顺序；
6. 主要流体密度；
7. 水泥浆和钻井液基本流变参数，或明确的模型假设；
8. 井斜数据或井斜范围。

若缺少以上关键项，不能直接作为主模型计算井。

### 16.2 能进入CBL定量验证的最低要求

除主计算要求外，还需要：

1. CBL评价段顶底深；
2. CBL合格率或可量化固井质量指标；
3. CBL资料来源明确；
4. CBL评价段与模型模拟段有重叠。

---

## 17. 给执行AI模型的简短指令模板

可直接把下面指令交给其他AI：

```text
请按《cement model 所需现场井尾管固井数据清单与提取规划》提取指定井的现场尾管固井数据。

只需要根据模型和论文需求提取字段，不要运行模型，不要修改cement model代码，不要整理历史模拟结果。

必须检查现场资料中的文字、表格、图片、扫描PDF、Word内嵌图片和曲线图。所有提取值必须记录来源、单位、置信度和是否需要人工复核。

最终按单井标准交付文件输出，并给出该井能否进入模型主计算、CBL定量验证或敏感性分析的建议分类。
```

---

## 18. 后续衔接

当所有井现场数据提取完成后，再进入论文执行流程：

1. 正式井资料完整性审查；
2. A/B/C/D样本分类确认；
3. 主模型参数固化；
4. 多井主结果生成；
5. CBL验证；
6. D2DGA消融；
7. 网格和时间步验证；
8. 论文Markdown初稿；
9. 《石油勘探与开发》格式Word终稿。
