# 模型数字地基收口（P0）+ Track B 物理改进路线图 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 石油学报论文优先 —— 收口模型"数字地基"（P0，约 1–1.5 周，纯正确性/口径，无物理风险），让论文数字可信可追溯；同时给出 Track B 物理改进的可行性评估与阶段路线（英文线备用）。

**依据:** `obsidian/固井顶替效率改进/尾管顶替效率模型现状与改进空间_综合调研_2026-08-29.md` §7（以下称"综合报告"）。四地基问题：①8 井 runner 与 scripts 口径分裂；②弥散"钉死"结果；③居中度敏感被六层机制栈掩盖；④CFL 口径未裁定 + CBL 对照链已断。

> **⚠️ 可行性审查修正（2026-08-29，执行前并入，优先于下文对应 Task 原文）**：
> - **R1(A1)**：`ht1_004 实算 0.8677` 为 T1-7 时代旧配置软证据（当前可追溯实算：修正配置 0.1418 / runner 存量 0.6884），**不得作现口径引用**；A1 以 Step 2 实测（CFL on vs 固定 dt）为准，0.8677 仅可作历史注脚。
> - **R2(B1/Phase B)**：8 井存量 summary 均无"评价窗效率"键（生成早于 M0），npz 无 b 场无法离线重建 → **模型侧窗口数字只能来自 A3 重跑后的新 summary**；Phase B"可先用 08-16 基线试跑"删除，B1→B2 排期按"A3 完成后"计。
> - **R3(A2/全局)**：新增 **Task 0（前置）**：先提交工作区未提交改动（`cemdisp/models2d/boundary_bridge.py`、`cemdisp/transport1d/casing_flow.py`、`cemdisp/PACKAGE_REFERENCE.md`），保证 `adopted_config.json` 的 commit 追溯与运行代码一致。A2 Step 3 中"12/9 元组签名修正"无对应物，改为修正 PACKAGE_REFERENCE.md 两处真实错误：①`:275` loader 签名缺 `measured_standoff` 参数；②`:444-445` 仍列出已删除的 metrics 列 `target_interval_efficiency`/`cbl_eval_interval_efficiency`。
> - **R4(C2/D1)**：C2 需先扩展 `ht1_004_ablation.py::run_one_level` 转发修正配置与 `enable_local_i3`；D1 为"`run()` 加可选 `schedule=None` + 8 井 runner 接线"（非单行），并顺手修正该文件头过时 docstring。
> - **R5(全局/验收)**：pytest 验收须 `--ignore=tests/paper_data`（孤儿测试既有结论）；C4 standoff 扫描重估：8 井×多点约 1.5–4 小时，收缩为 **hu101（有实测剖面）+ 1–2 口代表井**。

**Architecture:** 本计划分两类工作——
- **P0 收口**：改 runner/脚本/文档口径，不新增物理机制。所有"数字变更"均为**预期行为变更**（区别于失真修正计划的"开关默认关=逐位复现"约束）；涉及求解器默认值的修改仍须走开关化+默认保持。
- **Track B 物理**：仅概述与可行性评估，不在此展开 Task 级（待石油学报投稿后立项，另出规格）。

**Tech Stack:** Python 3、NumPy、pytest。conda 环境 `shenjingwangluo`。测试命令 `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/ -v`。

**基线:** commit `8e0e5ae`（环空失真修正合并）为主基线；`results/全井修正前后/汇总.csv`（nz=250）为 8 井修正后数字基线。

## Global Constraints

- **口径唯一性**：论文引用的每个数字必须能回溯到"哪个脚本+哪个配置"；不允许同一数字存在两个来源口径。
- **决策需落盘**：每个裁定（CFL 口径、runner 双轨、I3 口径、弥散系数）必须在文档中记录"选什么+为什么+何时重估"，禁止只在对话中口头裁定。
- **不新增物理**：P0 阶段不引入新常数/新闭包；弥散溯源（C1）若需改系数，必须走开关化+默认保持逐位复现（M1 同款）。
- **不删生产代码**：stale 脚本修复优先，无法修复的归档标注（加 `DEPRECATED` 头注释+README 说明），不直接删。
- **中文注释、UTF-8**；每个 Task 一个 commit，信息格式 `fix(data): <描述>` 或 `docs(paper): <描述>`，附验收结果。
- **论文联动**：每个 Task 产出的论文素材（表/段落/声明）统一收集到 `论文撰写/石油学报数据准备_2026-08-29/`（obsidian），供写作直接取用。

### 文件结构（本计划产出）

| 文件 | 责任 | 动作 |
|------|------|------|
| `cemdisp/runners/hu102_tailpipe.py` | 补 `split_cement_phases=True`（7/8 井统一） | 修改 |
| `cemdisp/runners/*_tailpipe.py`（其余 7 井） | 视 A2 裁定：显式化配置 or 标注双轨 | 修改/注释 |
| `scripts/rerun_all_wells_corrected.py` | 官方 8 井结果入口（修正后配置） | 完善（输出口径声明 JSON） |
| `scripts/cbl_window_comparison.py` | CBL 窗口口径对照模块（新建） | 新建 |
| `cemdisp/models2d/annulus_d2dga.py` | Tier0 内嵌调用补传 schedule（T0-6 接线） | 修改 |
| `scripts/p3_p4_integration.py` 等 4 个 stale 脚本 | 修复 or 归档标注 | 修改/标注 |
| `docs/superpowers/specs/2026-08-29-cfl-caliber-ruling.md` | CFL 口径裁定落盘（方案乙） | 新建 |
| `docs/superpowers/plans/2026-08-29-model-foundation-closure-p0-and-trackb-roadmap.md` | 本计划 | 新建 |
| obsidian `论文撰写/石油学报数据准备_2026-08-29/` | 论文素材收集（表/声明/说明段） | 新建 |

### 阶段总览与依赖

```
Task 0  前置：提交工作区未提交改动（boundary_bridge/casing_flow/PACKAGE_REFERENCE）
Phase A 口径裁定+runner 收口 (A1→A2→A3)        ← 数字地基，最先
Phase B CBL 窗口对照闭环 (B1→B2→B3)            ← 模型侧数字依赖 A3 重跑（R2）
Phase C 弥散溯源+论文核验 (C1/C2/C3 并行→C4)   ← 与 A/B 独立，可并行
Phase D 诊断与遗留修复 (D1/D2/D3 独立小修)      ← 随时可做
```

**硬约束：A3（8 井最终重跑）必须先于 B2（窗口对照表定稿）与任何论文数字引用；A1（CFL 裁定）必须先于 A3。**

**工期估计（合计约 1–1.5 周）**：Phase A 2–3 天；Phase B 2–3 天；Phase C 2–3 天；Phase D 1 天。

---

## Phase A：口径裁定与 runner 收口（数字地基）

## Task A1: CFL 口径裁定落盘（方案乙）

**Files:**
- New: `docs/superpowers/specs/2026-08-29-cfl-caliber-ruling.md`
- Data: 固定 dt=4s 扫描附表（复用 `scripts/rerun_all_wells_corrected.py` 或独立小脚本）

**Interfaces:**
- Produces: 口径裁定文档（选"方案乙"：主结果 CFL on + 固定 dt=4s 扫描附表），含三类数字的官方引用规则：
  1. 主结果表（CFL on，nz=250，修正后配置）
  2. 附表（固定 dt=4s，同配置）——证明结论对时间步口径不敏感
  3. 消融表（CFL off，nz=500）——仅用于机制归因，**不得与主结果同表混引**

- [ ] **Step 1:** 写口径裁定文档初稿：引用综合报告 §2.3 的 30pp 差异、08-23 设计"方案乙"裁定、ht1_004 实算 0.8677。
- [ ] **Step 2:** 跑固定 dt=4s 附表（8 井，修正后配置），确认与 CFL on 主结果的差异是否<2%（同 nz）；若>2%，记录并在文档中给出解释（弥散 dt 归一后的残差）。
- [ ] **Step 3:** 落盘最终裁定：明确论文正文数字=CFL on 主结果，附表=固定 dt；消融口径单独声明。
- [ ] **Step 4:** 在 obsidian `论文撰写/石油学报数据准备_2026-08-29/` 写"CFL 口径声明"段落素材。

**Acceptance:** 裁定文档含"三类数字引用规则"；8 井 CFL on vs 固定 dt 差异表已产出；论文可据此写口径声明。

## Task A2: 8 井 runner 口径统一（hu102 split_cement_phases + 双轨声明）

**Files:**
- Modify: `cemdisp/runners/hu102_tailpipe.py`（补 `split_cement_phases=True`）
- Modify: `cemdisp/runners/*_tailpipe.py`（其余 7 井）：视裁定加配置显式化注释
- Modify: `scripts/rerun_all_wells_corrected.py`：输出 `adopted_config.json`（记录启用开关清单+版本）

**Interfaces:**
- Produces: 8 井 runner 与 corrected 脚本口径一致（hu102 不再例外）
- Produces: `adopted_config.json` 含 solver 开关快照（M1/M2/M3/I3/e_clip_max/split_cement_phases/CFL），随结果存档

- [ ] **Step 1:** 修 hu102 runner 补 `split_cement_phases=True`；跑 hu102 对比修前修后差异，记录在 commit 说明（预期 η_E 变化，因 lead/tail 入口分配改变）。
- [ ] **Step 2:** 裁定 runner 定位（二选一，写进计划附录）：
  - (a) runner 改为"修正后配置"默认（论文=runner 口径，脚本仅留隔离分析）
  - (b) runner 保持默认基线，`rerun_all_wells_corrected.py` 为官方论文入口（论文=scripts 口径），runner 文件头加醒目双轨声明
  - **建议 (b)**：风险最小（不改 8 个 runner 默认行为，回归测试不受影响），论文口径集中在官方脚本。
- [ ] **Step 3:** 按裁定落地：若 (b)，给 8 个 runner 文件头加"当前为基线口径，论文数字见 rerun_all_wells_corrected.py"声明 + `PACKAGE_REFERENCE.md` 同步更新（含 12/9 元组签名修正）。
- [ ] **Step 4:** `rerun_all_wells_corrected.py` 增加 `adopted_config.json` 导出（含 git commit 号，`git rev-parse HEAD`）。

**Acceptance:** hu102 与其他井同口径；任何 8 井数字都能通过 `adopted_config.json` 追溯配置与代码版本。

## Task A3: 8 井最终基线重跑 + 新旧数字对照表（弃用旧数字）

**Files:**
- Data: `results/最终基线_2026-08-29/`（新目录，避免覆盖历史）
- Data: `results/最终基线_2026-08-29/新旧数字对照表.md`（8 井：旧论文数字/08-16 基线/08-24 修正/最终基线 四列）

**Interfaces:**
- Produces: 论文用 8 井新基线表（η_E/η_N/mixing/wall%/CFL 口径/配置版本）
- Consumes: Task A1（CFL 裁定）、A2（口径统一）结果

- [ ] **Step 1:** 用官方入口（A2 裁定后的）重跑 8 井，nz=250，修正后配置，输出到 `results/最终基线_2026-08-29/`。
- [ ] **Step 2:** 生成新旧数字对照表：四列对照 + "弃用理由"列（旧数字基于病态冻结值/旧口径）。
- [ ] **Step 3:** 更新 obsidian `论文撰写/石油学报数据准备_2026-08-29/8井最终基线表.md`（含 CFL 口径声明引用）。
- [ ] **Step 4:** 全量 pytest 回归（确保 runner 改动不破坏测试）。

**Acceptance:** 8 井最终数字齐备且可追溯；对照表明确标注哪些旧数字作废（hu103/hu1≈0、67% 阈值、"20 秒"卖点等）；pytest 全绿。

---

## Phase B：CBL 窗口对照闭环

## Task B1: CBL 窗口口径对照模块

**Files:**
- New: `scripts/cbl_window_comparison.py`
- Modify（若需）: `cemdisp/data/well_spec.py` 的 EvaluationWindow（确认四口径 cbl/cbl_quality/formation_target/model_focus 可被脚本直接消费）

**Interfaces:**
- Consumes: 各井 `WellSpec.EvaluationWindow`（M0 已算 per-window η_E/η_N）、`参考文档/现场资料提取/<井>/cbl_evaluation.csv`
- Produces: 每井"窗口对照表"（窗口 η_E/η_N vs 窗口 CBL 合格率 vs 窗口长度 vs 可对比性标记）

- [ ] **Step 1:** 写脚本骨架：读 loader → 取 EvaluationWindow 列表 → 对每个窗口输出模型 η_E/η_N（来自 summary 或结果 npz）→ 读 cbl_evaluation.csv 对应窗口 → 合并成表。
- [ ] **Step 2:** 处理"窗口错配"：hu101 模型域 5400–7868 vs CBL 窗 5699.8–7810 → 交集窗口口径；ht1_004 模型域含大量 CBL 不可评价段 → 双口径输出（全域/可评价段）。
- [ ] **Step 3:** 修复 4 个 stale 脚本（`p3_p4_integration.py:142` 等 KeyError）——优先修，修不动则加 `DEPRECATED` 头注释 + 指向 cbl_window_comparison.py。

**Acceptance:** 8 井窗口对照表可复现产出；stale 脚本不再报 KeyError（修复或已标注归档）。

## Task B2: 8 井窗口对照表 + 结构性错配说明段（论文素材）

**Files:**
- Data: `results/最终基线_2026-08-29/CBL窗口对照总表.csv`
- Doc: obsidian `论文撰写/石油学报数据准备_2026-08-29/CBL对照说明段.md`

- [ ] **Step 1:** 汇总 B1 各井窗口表为总表，标注"可比/部分可比/不可比"三档（依据窗口长度覆盖与判据口径）。
- [ ] **Step 2:** 写论文 CBL 对照说明段素材：含三重结构性错配声明（全井均值 vs 评价段、泵停 vs 候凝 24h+、体积分数 vs 胶结等级）+ 各井偏差归因（hu103 漏失、ht1_004 域错配、hu101 待窗口复算）。
- [ ] **Step 3:** 明确"模型不对标 CBL、CBL 仅作对比"的表述（三条铁律之一），写入论文讨论段素材。

**Acceptance:** 论文 §CBL 对照段落素材齐备（表+说明段+归因）。

## Task B3: 呼101 −10.6pp 窗口口径复算

- [ ] **Step 1:** 用 B1 模块对 hu101 做 CBL 窗口（5699.8–7810）口径复算：模型窗口 η_E/η_N vs CBL 62.77%。
- [ ] **Step 2:** 记录复算结果与旧基线"65.04 vs 62.77 (+2.3pp)"叙事的关系：确认窗口口径下偏差是否收窄或反转，写进 B2 说明段。
- [ ] **Step 3:** 若窗口口径下 hu101 偏差仍大，列出可解释项（standoff 假设、流变单点、漏失注灰施工异常）供讨论段引用。

**Acceptance:** hu101 窗口口径数字产出并有归因说明（不再停留在全井均值 −10.6pp 单一口径）。

---

## Phase C：弥散溯源与论文数字核验

## Task C1: 弥散系数溯源（0.018/0.015 vs Zhang22 间隙尺度弥散）

**Files:**
- Doc: `docs/superpowers/research/dispersion_coefficient_origin_2026-08-29.md`
- 若需改系数：`cemdisp/models2d/annulus_d2dga.py` 走开关化（M1 同款），默认保持逐位复现

- [ ] **Step 1:** 精读 Zhang22 的间隙尺度弥散推导（§4 弥散项），核对当前 0.018/0.015 是否对应原文某参数化；查 Z23 实验的弥散量级。
- [ ] **Step 2:** 给出结论：系数有据（改/不改）/ 无据（必须处理）。若无据，量化"弥散主导"风险（引用 HT1-004 scale 扫描 1.7%↔66.5%），给出最小处理方案（如：系数公开声明为标定参数+敏感性附表，或改走物理推导式）。
- [ ] **Step 3:** 产出论文可用的"弥散处理声明"段落素材（含敏感性附表或溯源结论）。

**Acceptance:** 弥散系数来源结论落盘；论文有对应声明素材；若改动，走开关化且默认逐位复现。

## Task C2: I3 修复后 R2 消融重测

- [ ] **Step 1:** 跑 HT1-004 R0→R3 消融（`ht1_004_ablation.py`），确认 R2 是否仍 ≈R1（修复前 +2e-10），记录新值。
- [ ] **Step 2:** 补跑 `enable_local_i3=True` 对比（eta2/Δρ 局部化后 I3 是否生效）。
- [ ] **Step 3:** 结果写入论文消融表素材（若 I3 仍近零，如实报告并归因：垂直井 f_phi→0 使 I3 退化为纯轴向分量 + 弥散主导掩盖）。

**Acceptance:** R2 消融新数据产出；I3 效果结论更新（有效/仍无效+原因）。

## Task C3: I3 双口径核对签字

- [ ] **Step 1:** 核对代码 I3 公式 = Zhang22 式(4.26)；对照 Bararpour25 式(2.27) 差异（分子 3(1−c̄) vs 3(1−c̄²)、分母 2m[…] vs 12[…]）。
- [ ] **Step 2:** 裁定论文口径（选一个为主+注明归一化差异）；产出"严禁写相差 6/m 因子"的核对签字记录（obsidian 数据准备目录）。
- [ ] **Step 3:** 论文方法章 I3 公式段落素材定稿。

**Acceptance:** I3 口径签字记录落盘；方法章公式段落素材定稿。

## Task C4: 论文数据准备清单收口（standoff ±0.1、T1-10/T1-11 决策记录、η_N 三件套）

- [ ] **Step 1:** standoff ±0.1 区间：对 8 井输出 η_E/η_N 对 standoff±0.1 的敏感性区间（复用 `standoff_sensitivity_scan.py` 思路，修正后配置）。
- [ ] **Step 2:** T1-10/T1-11 接线决策记录：写"混浆/Atwood 已接线生产+接线日期+消融口径（弥散 −1pp、重力 −6.5pp）"入论文素材。
- [ ] **Step 3:** η_N 三件套论证素材：η_N + 评价窗效率 + 窜槽指数 作为居中度敏感性主论证（引用 Foolad 2021"体积效率欺骗性"）；η_E 降级为辅助指标——写入论文方法/讨论素材。
- [ ] **Step 4:** 汇总所有素材到 obsidian `论文撰写/石油学报数据准备_2026-08-29/`，生成 `清单核对.md`（逐项打勾）。

**Acceptance:** 论文数据准备清单全绿；每项素材可定位。

---

## Phase D：诊断与遗留修复

## Task D1: Tier0 T0-6 传 schedule 接线

**Files:**
- Modify: `cemdisp/models2d/annulus_d2dga.py:1351`（内嵌诊断调用补传 schedule）

- [ ] **Step 1:** 在求解器 run() 内嵌诊断调用处传入 schedule（从 self 或 runner 传入），使 T0-6 停泵衰减不再恒跳过。
- [ ] **Step 2:** 单测：跑任一井确认 `summary["tier0_diagnostics"]` 含 shutdown_decay 非 None。
- [ ] **Step 3:** 确认开关关闭路径（无 schedule 时）仍优雅降级（记 notes 而非报错）。

**Acceptance:** T0-6 在 runner 结果中产出；回归测试全绿。

## Task D2: stale 脚本清理

- [ ] **Step 1:** 列出引用已删除路径的脚本（`p3_p4_integration.py`、`run_calibrated_all_wells.py` 等，已知 4 个）。
- [ ] **Step 2:** 逐脚本判定：修复（指向新模块）或加 `DEPRECATED` 头注释+原因（引用旧 cbl 指标路径已删）。
- [ ] **Step 3:** `scripts/README.md` 或 `cemdisp/PACKAGE_REFERENCE.md` 更新脚本清单状态。

**Acceptance:** 无脚本运行即 KeyError；每个 stale 脚本有明确状态标注。

## Task D3: M2 re_crit 标定决策

- [ ] **Step 1:** 记录 M2 现状：`re_crit=2100(1+0.1He)` 临时标定钮；对现 8 井中性（Bingham 泥浆 He 高→re_crit 高→全层流 R=1）。
- [ ] **Step 2:** 裁定：石油学报阶段不标定（对现 8 井无影响），文档声明"启用前须三重回归（Walton&Bittleston / Z22 Table3 / Foolad 定性）"；标注在 Track B 前置条件。
- [ ] **Step 3:** 在 M2 代码注释与 spec 中加"未标定"醒目标注。

**Acceptance:** M2 标定状态与启用条件有明确文档记录。

---

## Track B：物理改进路线图（概述 + 可行性评估，英文线备用）

> 石油学报投稿后立项。每个项目另出 spec + Task 级计划（TDD）。以下为可行性评估与排序。

| 项 | 内容 | 可行性 | 前置 | 工作量 | 论文杠杆 | 风险 |
|---|---|---|---|---|---|---|
| **B1** | 弥散 TVD / 间隙尺度弥散（去人工弥散主导） | 高 | C1 溯源结论 | ~1 周 | **高**（解"钉死"，结果不再被 0.018/0.015 主导） | 需防数值扩散回潮；TVD 限幅器选型 |
| **B2** | 真实横流（椭圆流函数求解 v，替换连续性回推） | 中 | M5 设计；需重排速度求解 | 3–5 天 | 高（窄边"饿死"正反馈的物理解，I3 才有载体） | 破坏基线逐位复现；需回归锚定 |
| **B3** | static wall 物理化（Bararpour25 c̄min 判据替换 c_min=0.05） | 高 | M3 基础（yield gate 已就位） | 待估 | 中高（文献级门槛落地） | 判据公式需逐字符核对 |
| **B4** | 窄缝穿透门槛 | 中 | B2 | 待估 | 中 | 物理边界待文献核实 |
| **B5** | Uzawa 精确屈服（augmented Lagrangian） | 中 | 求解器重构（Papanastasiou→Uzawa） | 2–3 周 | 中高（真 plug/static layer，η_N 可信化） | 周期长；数值收敛 |
| **B6** | 瞬态（停泵滑移/凝胶发展） | 低 | 架构不支持（浓度场冻结） | 2–3 周 | 低（周期长，石油学报不用） | 需重设计停泵分支 |
| **B7** | 序列优化器（浆柱结构优化） | 高 | — | 1–2 周 | 中（工程应用价值） | 与论文创新点重叠风险 |

**排序建议（若做英文线）**：**B1 → B2 → B3**（与综合报告 §7.3 杠杆排序一致）；B5 视时间；B6/B7 缓。
**硬约束**：B2/B3 依赖 A4 已完成的 e 去截断（M4）与 M3 就位——**已满足**（8e0e5ae）；B1 必须基于 C1 溯源结论，禁止凭空改系数。

---

## 新方向评估（不进入近期计划，仅记录）

| 方向 | 文献 | 评估 |
|---|---|---|
| 控压固井温压-流变耦合 | 刘金璐 2025 石油科学通报/2024 钻采工艺 | **与石油学报控压固井卖点直接相关**：可先作为论文对照引用；耦合实现放 M6（流变温压修正） |
| 套管旋转 | PoF 2024/2025 | 远期；控压固井中旋转套管提高顶替的量化 |
| 泡沫水泥 | Moreira 2025 | 低密度体系扩展，与控压固井低密度相关但工作量大 |
| PINN 代理 | EAAI 2025 | 8 井快速重算/敏感性代理，纯工具价值 |

---

## 里程碑与验收

| 里程碑 | 内容 | 完成判据 |
|---|---|---|
| **M1** | Phase A 完成 | 数字地基裁定落盘（CFL/runner 双轨/8 井最终基线）；论文可开始引用数字 |
| **M2** | Phase B 完成 | CBL 窗口对照闭环；论文 §CBL 对照素材齐备 |
| **M3** | Phase C+D 完成 | 论文数字核验清单全绿（弥散/I3/standoff/η_N 三件套素材齐备）；T0-6 接线；无 stale 脚本 |

**总验收**：obsidian `论文撰写/石油学报数据准备_2026-08-29/清单核对.md` 全勾；`results/最终基线_2026-08-29/` 数字可追溯；pytest 全绿；每个 Task 一个 commit。

---

## 附录：P0 与石油学报论文审计 P0 的对照

综合报告 §6.2 + 记忆"论文框架审计 2026-08-16"的 P0 项 → 本计划映射：

| 论文审计 P0 项 | 本计划 Task |
|---|---|
| CFL 口径裁定 | A1 |
| 表 4 窜槽列错配 / 65.04% 指标删除 / R3 公式不符 | A3（新旧数字对照）+ C3（I3 口径） |
| 旧论文数字（hu103/hu1≈0、"20 秒"卖点）弃用 | A3 |
| CBL 段效率未实现 | B1/B2/B3 |
| standoff ±0.1 区间 | C4 |
| T1-10/T1-11 接线决策记录 | C4 |
| I3 双口径核对签字 | C3 |
