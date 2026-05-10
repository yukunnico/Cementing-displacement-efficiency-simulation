# 0708资料整理与多井数据接入实施计划（改进版）

> **版本**: 1.1
> **修订日期**: 2026-05-10
> **修订说明**: 针对初版计划进行四处关键补强——补充呼101井、补全双径向字段模板、增加提取数据回写闭环、明确证据等级与provenance状态映射。

**Goal:** 把 0708 资料按"证据化提取 → 标准字段映射 → 参考文档沉淀 → Loader回写 → Provenance升级"的完整链路落地，先以呼103为样板井跑通四类提取表模板，再逐步扩展到呼101、呼102和呼探1系列。

**Architecture:** 本计划覆盖资料治理层与代码更新层。核心产物是：四类 CSV 提取模板、六井样板提取数据、`参考文档/` 沉淀结构、以及提取数据回写到 `cemdisp/data/loaders/` 和 `cemdisp/data/provenance.py` 的闭环。

**Tech Stack:** Markdown、CSV、pandoc（docx 提取）、主代理 PDF 直读、Python（回写验证）

---

## 与初版计划相比的改进点

| 序号 | 改进项 | 初版状态 | 改进后 | 影响范围 |
|------|--------|---------|--------|----------|
| 1 | 呼101井缺失 | 未覆盖 | 新增 Task 9（呼101资料整理） | 六井全覆盖 |
| 2 | WellSpec模板缺双径向字段 | 无上段内径/鞋口滞后体积字段 | 补充 `upper_section_bottom_md_m` 等5个字段 | 边界同步精度 |
| 3 | 提取→代码闭环断裂 | 到缺口汇总即结束 | 新增 Task 11：Loader回写与Provenance升级 | 资料真正生效 |
| 4 | 证据等级与provenance无映射 | 未定义A/B/C→field/partial/proxy的映射 | 明确映射规则，并在Task 11批量升级 | 可追溯性 |

---

## 文件结构总览

本计划涉及的文件：

| 文件 | 职责 |
|---|---|
| `参考文档/提取模板/well_spec_template.csv` | WellSpec 提取表模板（含双径向字段） |
| `参考文档/提取模板/fluid_spec_template.csv` | FluidSpec 提取表模板 |
| `参考文档/提取模板/pumping_schedule_template.csv` | PumpingSchedule 提取表模板 |
| `参考文档/提取模板/validation_data_template.csv` | ValidationData 提取表模板 |
| `参考文档/提取模板/原始资料索引_template.md` | 原始资料索引模板 |
| `参考文档/提取模板/README.md` | 模板使用说明 |
| `参考文档/呼101/原始资料索引.md` | 呼101 资料来源清单（新增） |
| `参考文档/呼101/提取数据/*.csv` | 呼101 四类提取表（新增） |
| `参考文档/呼103/原始资料索引.md` | 呼103 资料来源清单 |
| `参考文档/呼103/提取数据/*.csv` | 呼103 四类提取表 |
| `参考文档/呼103/提取数据/extraction_notes.md` | 呼103 提取过程备注 |
| `参考文档/呼102/原始资料索引.md` | 呼102 资料来源清单 |
| `参考文档/呼102/提取数据/*.csv` | 呼102 四类提取表 |
| `参考文档/呼探1/原始资料索引.md` | 呼探1 资料来源清单 |
| `参考文档/呼探1/提取数据/*.csv` | 呼探1 四类提取表 |
| `参考文档/呼探1-001/原始资料索引.md` | 呼探1-001 资料来源清单 |
| `参考文档/呼探1-001/提取数据/*.csv` | 呼探1-001 四类提取表 |
| `参考文档/呼探1-002/原始资料索引.md` | 呼探1-002 资料来源清单 |
| `参考文档/呼探1-002/提取数据/*.csv` | 呼探1-002 四类提取表 |
| `参考文档/多井资料缺口汇总.md` | 缺口汇总与 fallback 策略 |
| `docs/superpowers/plans/2026-05-09-0708-data-mapping-implementation-v1.1.md` | 本计划文档 |

---

## Task 1: 创建四类提取表模板（改进）

**Files:**
- Create/Modify: `参考文档/提取模板/well_spec_template.csv`
- Create: `参考文档/提取模板/fluid_spec_template.csv`
- Create: `参考文档/提取模板/pumping_schedule_template.csv`
- Create: `参考文档/提取模板/validation_data_template.csv`
- Create: `参考文档/提取模板/原始资料索引_template.md`
- Create: `参考文档/提取模板/README.md`

### Step 1: 创建 `well_spec_template.csv`（改进版）

```csv
字段名,字段值,单位,适用井段/套管段,来源文件,页码/图号,证据等级,备注
well_name,,,-,-,-,A,井号
top_md_m,,m,-,-,-,A,井段顶部深度
bottom_md_m,,m,-,-,-,A,井段底部深度
shoe_md_m,,m,-,-,-,A,鞋深（测深）
hanger_md_m,,m,-,-,-,A/B,悬挂器位置
casing_od_mm,,mm,-,-,-,A,套管外径
casing_id_mm,,mm,-,-,-,A,套管内径
liner_od_mm,,mm,-,-,-,A,尾管外径
liner_id_mm,,mm,-,-,-,A,尾管内径
liner_wall_thickness_mm,,mm,-,-,-,A/B,尾管壁厚
upper_section_bottom_md_m,,m,-,-,-,A/B,上段尾管底界（双径井必填）
upper_liner_od_mm,,mm,-,-,-,A/B,上段尾管外径（双径井必填）
upper_liner_id_mm,,mm,-,-,-,A/B,上段尾管内径（双径井必填）
shoe_lag_volume_m3,,m3,-,-,-,A/B,鞋口滞后体积（替代管内截面积×鞋深粗算）
hole_diameter_mm,,mm,-,-,-,A/B,井眼直径
packer_depth_m,,m,-,-,-,A/B,分级箍位置
evaluation_window_top_m,,m,-,-,-,A/B,评价窗口顶部
evaluation_window_bottom_m,,m,-,-,-,A/B,评价窗口底部
cementing_method,,,-,-,-,A,固井方式（单级/分级）
```

> **改进说明**：补充了 `upper_section_bottom_md_m`、`upper_liner_od_mm`、`upper_liner_id_mm`、`shoe_lag_volume_m3`、`liner_wall_thickness_mm` 五个字段。这些字段直接影响 Task 3–5 中 `CasingFlowSolver._timeline_pipe_volume()` 的计算精度。

### Step 2: 创建 `fluid_spec_template.csv`

```csv
流体名称,流体角色,密度(g/cm3),流变模型,n,K(Pa·sn),PV(mPa·s),YP(Pa),配方摘要,稳定性(g/cm3),冲洗效率(%),API失水(mL),24h强度(MPa),稠化时间(min),相容性摘要,来源文件,页码/表号,证据等级,备注
```

### Step 3: 创建 `pumping_schedule_template.csv`

```csv
步骤序号,步骤名,流体名称,设计体积(m3),实际体积(m3),设计排量(m3/min),实际排量(m3/min),设计密度(g/cm3),实际密度(g/cm3),开始时间,结束时间,泵压(MPa),备注,来源文件,页码/图号,证据等级
```

### Step 4: 创建 `validation_data_template.csv`

```csv
资料类型,文件路径,测量井段顶(m),测量井段底(m),完钻井深(m),固井方式,固井质量结论,水泥面预计(m),水泥面实测(m),固井日期,测井日期,解释日期,井液密度(g/cm3),水泥密度(g/cm3),来源文件,页码,证据等级,备注
```

### Step 5: 创建 `原始资料索引_template.md`

```markdown
# <井号> 原始资料索引

## 资料清单

| 序号 | 来源文件路径 | 文件类型 | 主要用途 | 证据等级上限 | 是否已完成抽取 |
|---|---|---|---|---|---|
| 1 | | | | | |

## 说明

- 主要用途分类：几何 / 流体 / 施工 / 验证
- 证据等级：A（原始正文） / B（二次汇总） / C（类比兜底）
- 证据等级与代码 provenance 状态映射：
  - A → `field`（现场符合）
  - B → `partial`（部分符合）
  - C → `proxy`（代理/暂定）
```

### Step 6: 创建 `README.md` 使用说明

写明模板用途、各列含义、证据等级定义、使用注意事项，以及**证据等级与 `cemdisp.data.provenance` 状态字段的映射规则**。

### Step 7: 验证模板可读性

用文本编辑器打开所有 CSV，确认列对齐、无乱码、中文表头正确。

---

## Task 2: 建立呼103原始资料索引

**Files:**
- Create: `参考文档/呼103/原始资料索引.md`

步骤与初版一致，但在说明部分增加证据等级映射规则。

---

## Task 3: 呼103 WellSpec 提取

**Files:**
- Create: `参考文档/呼103/提取数据/well_spec.csv`

步骤与初版一致，但需按**改进版模板**填写，包含双径向字段（呼103上段内径为代理值，证据等级标注为 B）。

---

## Task 4: 呼103 FluidSpec 提取

**Files:**
- Create: `参考文档/呼103/提取数据/fluid_spec.csv`

步骤与初版一致。

---

## Task 5: 呼103 PumpingSchedule 提取

**Files:**
- Create: `参考文档/呼103/提取数据/pumping_schedule.csv`

步骤与初版一致。

---

## Task 6: 呼103 ValidationData 提取

**Files:**
- Create: `参考文档/呼103/提取数据/validation_data.csv`

步骤与初版一致。

---

## Task 7: 呼103 提取过程备注

**Files:**
- Create: `参考文档/呼103/提取数据/extraction_notes.md`

步骤与初版一致，但增加一节：
- **证据等级映射记录**：列出哪些字段已满足 A 级可升级为 `field`，哪些仍为 B/C 级保持 `partial`/`proxy`。

---

## Task 8: 呼102 资料整理

**Files:**
- Create: `参考文档/呼102/原始资料索引.md`
- Create: `参考文档/呼102/提取数据/well_spec.csv`
- Create: `参考文档/呼102/提取数据/fluid_spec.csv`
- Create: `参考文档/呼102/提取数据/pumping_schedule.csv`
- Create: `参考文档/呼102/提取数据/validation_data.csv`

步骤与初版一致，但 `well_spec.csv` 需补充 `shoe_lag_volume_m3` 字段（呼102为单径井，该字段可留空或按单径计算值填入）。

---

## Task 9: 呼探1系列资料整理

**Files:**
- Create: `参考文档/呼探1/原始资料索引.md`
- Create: `参考文档/呼探1/提取数据/*.csv`
- Create: `参考文档/呼探1-001/原始资料索引.md`
- Create: `参考文档/呼探1-001/提取数据/*.csv`
- Create: `参考文档/呼探1-002/原始资料索引.md`
- Create: `参考文档/呼探1-002/提取数据/*.csv`

步骤与初版一致。呼探1和呼探1-001为双径向井，`well_spec.csv` 必须填写 `upper_section_bottom_md_m`、`upper_liner_od_mm`、`upper_liner_id_mm`。

---

## Task 10: 呼101 资料整理（新增）

**Files:**
- Create: `参考文档/呼101/原始资料索引.md`
- Create: `参考文档/呼101/提取数据/well_spec.csv`
- Create: `参考文档/呼101/提取数据/fluid_spec.csv`
- Create: `参考文档/呼101/提取数据/pumping_schedule.csv`
- Create: `参考文档/呼101/提取数据/validation_data.csv`

- [ ] **Step 1: 建立呼101原始资料索引**

从 `0708/` 中定位呼101相关资料。呼101为设计文档中明确列出的双径/变径等效型代表井，具有 `shoe_lag_volume_m3=52.0` 的特殊参数。

- [ ] **Step 2: 提取呼101四类提取表**

重点提取：
- 双径向几何：上段 `168.3mm`、下段 `139.7mm`、变径位置
- `shoe_lag_volume_m3` 的实测或计算依据
- 多前置液 + 多段替浆的详细施工步骤

- [ ] **Step 3: 标注缺口**

呼101施工程序最复杂（平衡液 → 驱油隔离液 → 领浆 → 尾浆 → 后置液 → 轻泥浆 → 中置液 → 井浆快替/慢替），需特别注意各段体积和排量的提取完整性。

---

## Task 11: 提取数据回写 Loader 与 Provenance 升级（新增）

**Files:**
- Modify: `cemdisp/data/loaders/hu101_loader.py`
- Modify: `cemdisp/data/loaders/hu102_loader.py`
- Modify: `cemdisp/data/loaders/hu103_loader.py`
- Modify: `cemdisp/data/loaders/hu1_loader.py`
- Modify: `cemdisp/data/loaders/hu2_loader.py`
- Modify: `cemdisp/data/loaders/ht1_001_loader.py`
- Modify: `cemdisp/data/provenance.py`

### Step 1: 将 CSV 提取表数据映射到 loader 常量

- 用提取表中的实测值替换 loader 中的代理值
- 保留原始代理值作为 fallback 注释（加 `# 旧代理值: ...`）
- 对双径向井，确保 `upper_section_bottom_md_m`、`upper_liner_id_mm` 等字段从提取表同步到 `WellSpec` 构造

### Step 2: 更新 provenance 状态标签

证据等级 → provenance 状态映射规则：

| 证据等级 | 含义 | provenance.status | 典型场景 |
|---------|------|-------------------|---------|
| A | 原始正文/实测 | `field` | CBL原始PDF、实验报告流变数据 |
| B | 二次汇总/设计值 | `partial` | 施工设计、汇总表交叉校验 |
| C | 类比兜底/文献值 | `proxy` | 邻井代理、历史脚本占位值 |
| 未找到 | — | `unknown` | 当前缺失，需后续补录 |

- 将已补齐的流体从 `"proxy"` / `"partial"` 升级为 `"field"`
- 在 `FluidProvenance.note` 中注明数据来源文件（如 `"密度与幂律参数来自 203211.docx 实验报告，A级提取"`）
- 对仍为代理值的字段，更新 `note` 说明当前 fallback 依据

### Step 3: 验证回写一致性

- 运行 loader 生成 `WellSpec`、`FluidSpec`、`PumpingSchedule`
- 与 CSV 提取表逐项比对，确保无遗漏
- 运行 `CasingFlowSolver`，检查 `shoe_timeline` 事件与提取表中的施工步骤是否一致

### Step 4: 运行端到端测试

```bash
python -m unittest tests.test_six_well_integration -v
```

确认所有六口井的同步画像卡正常生成，且 `proxy` 警告数量减少。

---

## Task 12: 缺口汇总与 fallback 策略文档

**Files:**
- Create: `参考文档/多井资料缺口汇总.md`

步骤与初版一致，但增加：
- [ ] **Step 4: 汇总 provenance 升级状态**

列出每口井各流体的当前状态分布（field / partial / proxy / unknown），与初版计划对比，量化资料补强进度。

---

## 验收标准

### 阶段 1 验收（Task 1–7：呼103样板化）

- [ ] 四类 CSV 模板创建完成，列定义与改进版 spec 一致（含双径向字段）
- [ ] 呼103 原始资料索引覆盖已确认的高价值资料
- [ ] 呼103 四类提取表填写完成，每条记录有来源文件和证据等级
- [ ] 呼103 extraction_notes.md 记录了提取过程和缺口
- [ ] 所有 CSV 文件可被文本编辑器正常打开，无乱码

### 阶段 2 验收（Task 8–10：呼102、呼探1系列、呼101）

- [ ] 呼102、呼101、呼探1、呼探1-001、呼探1-002 原始资料索引创建完成
- [ ] 五口井四类提取表填写完成
- [ ] 五口井缺口标注完成
- [ ] 双径向井的 `upper_section_bottom_md_m`、`upper_liner_id_mm` 等字段已提取

### 阶段 3 验收（Task 11：回写闭环）

- [ ] Loader 常量已从提取表更新，代理值保留注释
- [ ] `provenance.py` 中各井流体状态已按证据等级升级
- [ ] 六口井的 `shoe_timeline` 与提取表施工步骤一致
- [ ] `python -m unittest discover -s tests` 全部通过

### 阶段 4 验收（Task 12：缺口汇总）

- [ ] 多井资料缺口汇总文档创建完成
- [ ] fallback 策略评估完成
- [ ] provenance 升级状态汇总完成
- [ ] 后续补录优先级建议完成

---

## 与主计划的衔接

本计划（第二阶段：逐井资料补强）完成后，应回到主计划 `2026-05-09-multiwell-boundary-synchronization-design.md` 的第三阶段：

1. 用更新后的 loader 重新运行六口井的边界同步
2. 对比新旧 `ShoeTimeline` 的关键事件时间差异
3. 判断边界同步精度提升后，最终效率是否更贴近现场
4. 若仍有显著偏差，再评估是否需要增强 1D（对流-弥散）或 2D 核心

---

## 附录：呼101井关键参数备忘

来自设计文档 `2026-05-09-multiwell-boundary-synchronization-design.md` 第 5.1 节：

- **结构特征**：双径/变径等效型
- **程序特征**：平衡液 → 驱油隔离液 → 领浆 → 尾浆 → 后置液 → 轻泥浆 → 中置液 → 井浆快替/慢替
- **设计重点**：
  - 52m³ 鞋口滞后口径需要保留
  - lead / tail 分相要能继续支持
  - 上段 `168.3mm` 与下段 `139.7mm` 的等效几何不能在边界同步时丢掉

这些参数必须在 `well_spec.csv` 和 `pumping_schedule.csv` 中完整体现。
