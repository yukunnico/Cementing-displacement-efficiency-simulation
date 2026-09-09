# 套管段混浆带路线 B（接触时间积分 + 胶塞面零掺混）Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 把套管内混浆带宽度从"单一标定量（alpha=0.25·U·R 出口后处理 + 全程行程时间近似）"升级为"界面真实接触时间历程积分 × 既有 D_eff 闭包 + 胶塞面零掺混"。混浆带从"可调参数"变为"半预测量"，同时修复 09-03 胶塞语义裁定的一个实现遗漏。这是用户裁定的套管段实施方向（路线 B 主线先行）；路线 A（1D c(z,t) PDE）为二期预留，本计划的接触时间与 D_eff 闭包概念在 A 中整体复用。

**背景调研:** obsidian《套管段管内浓度剖面真解_调研与实施方案_2026-09-08》（§0 结论速览 / §二 三路线对比 / §六 实施卡）；探针实测 `.tmp_research/probe_mixing_width_20260908.py` + `_results.json`。执行者开工前必读该报告 §一/§二。

**Tech Stack:** Python 3、NumPy、pytest。conda 环境 `shenjingwangluo`；测试命令 `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/ -v`。

## Global Constraints

- **基线逐位复现**：双开关 `mixing_contact_time=False` / `plug_face_zero_mixing=False`（默认）时，所有井 shoe_timeline 逐位等于当前 main（c471ab7 之后工作区）基线。每个任务的测试都要验证。
- **验收红线**：开关全开时，8 井 η_E 全井 |Δ|<0.5pp、stop（cement_end_time_s）逐位不变；混浆指数允许有可感知变化（09-06 消融已证当前库存比≥1 工况下混浆带对 η 是弱通道）。若某井超线 → 停，回退分析（见 Risk）。
- **不加新物理常数**：D_eff 闭包完全复用 `_compute_dispersion_coefficient` + `_interface_instability_factor`；本计划只换 σ_t 里的时间项与跳过判据。
- **代码风格**：中文注释、关键物理逻辑注释保留、不做无关重构；每 Task 一个 commit（`feat(transport1d): <描述>（开关=xx 默认关）`）。
- **权威基线数字** = `results/三项修复重跑_2026-09-06/汇总.csv`（勿引用 runner 目录数字，口径分裂见 obsidian 09-08 三问裁定报告）。

### 文件结构（本计划产出）

| 文件 | 责任 | 动作 |
|------|------|------|
| `cemdisp/transport1d/casing_flow.py` | 1D 求解器：+2 构造参数、+`_contact_time_integrated_sigma`、`_apply_dispersion_to_timeline` σ_t 来源替换与胶塞面跳过 | 修改 |
| `tests/test_casing_mixing_contact_time.py` | 逐位回归 + 接触时间单调性 + 胶塞面跳过 + 分段注入边界 | 新建 |
| `.tmp_research/probe_mixing_width_20260908.py` | （已存在，只读探针）Task 0 对照复核用 | 只读 |
| `results/浓度剖面真解路线B_2026-09-XX/` | 8 井双开关对照表 | 新建 |
| obsidian《套管段管内浓度剖面真解…2026-09-08.md》 | 实施结果附记 | 修改 |

### 代码锚点（已核实，main 工作区）

| 锚点 | 内容 |
|------|------|
| `casing_flow.py:562-671` | `_apply_dispersion_to_timeline`：σ_t 生成主体（prev_fluid 查找 :598-609、混浆增强 :612-619、σ_t 计算 :626-632、5 子事件+F4 收尾 :637-671） |
| `casing_flow.py:626-631` | 现状 σ_t 公式：`t_travel = shoe_md_m / U`（全程近似）；`sigma_t = sqrt(2·D_eff·t_travel)/U`；防御上限 `min(σ, 0.5·t_travel)`、下限 `max(σ, dt)` |
| `casing_flow.py:903-918` | `_front_arrival_time`：体积坐标→时刻反解（target_volume = cumulative_start + pipe_volume；按步线性内插）——t_inject 反解复用此逻辑 |
| `casing_flow.py:74-84` | `_ScheduledStep`（step/start_time_s/end_time_s/cumulative_volume_start_m3/cumulative_volume_end_m3） |
| `casing_flow.py:225-245` | `run()`：`_timeline_pipe_volume`（双内径感知管容）→ `_displacement_sequence_cutoff`（胶塞 RESTART 截断） |
| `casing_flow.py:144/220` | `has_plug` 构造参数（09-03 胶塞语义开关，生产 has_plug=False、唯一 hu102 True） |
| `casing_flow.py:730` | `_apply_dispersion_to_timeline` 唯一调用处（`_build_shoe_timeline` 内） |
| `cemdisp/data/fluid_spec.py:50-60` | `FluidRole` 枚举：**压塞液在各 loader 全部标 `FluidRole.OTHER` + 名"压塞液"**（ht1_001:404/ht1_003:437/ht1_004:464,504/hu2:334）；无独立 PLUG 角色 |
| `cemdisp/models2d/boundary_bridge.py:72-97` | `_phase_fractions_for_fluid`：OTHER→mud 相（压塞液与替浆同相——胶塞面语义在 1D 层处理，不动桥接） |

---

## Task 0：基线冻结锚（0.5h）

- [ ] 跑 `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python .tmp_research/probe_mixing_width_20260908.py` 确认探针在当前工作区可复现（结果应与 `probe_mixing_width_20260908_results.json` 一致：hu101 过渡带 1.477 m³/占评价域 2.254%、hu102 13.5%、ht1_003 3.8%）。不一致 → 停，先对齐工作区状态。
- [ ] 跑一次全量 `pytest tests/`（363+ 过）记录基线。
- [ ] 保存基线锚：用 hu103（无 RESTART、界面齐全）构造 `CasingFlowSolver(enable_gravity=True).run(...)` 输出的 `shoe_timeline.events` 元组，把逐事件字段（time_s/kind/flow_rate_m3_s/stage_name/phase_fractions）写入 `tests/_baseline_shoe_timeline_hu103.json` 供 Task 1 逐位回归。

## Task 1：接触时间积分 σ_t（S1，~1 天）

- [ ] **测试先行**（`tests/test_casing_mixing_contact_time.py`）：
  1. `test_default_off_bitwise`：双开关默认，hu103 shoe_timeline 逐事件等于基线锚；
  2. `test_contact_time_monotone`：开关开 → 对 hu103 每个界面计算 t_contact = t_arrival − t_inject，断言 (a) 所有 t_contact ≤ t_travel，(b) 晚注入界面（尾浆前缘）σ_t < 基线 σ_t，(c) 首界面（钻井液→平衡液，注入≈开泵时刻）σ_t ≈ 基线（相对差 <5%）；
  3. `test_segmented_injection_boundary`：同流体分两段注入（如"尾浆 RATE_SWITCH + 尾浆 FRONT_ARRIVAL"）时界面取**该前缘所在注入步**的出发时刻（与 prev_fluid 跳过同名事件逻辑对齐，不串到第一段）。
- [ ] **实现**：`CasingFlowSolver.__init__` 增 `mixing_contact_time: bool = False`；新增静态方法：
  ```
  _inject_start_time(front_step, scheduled_steps_full) -> float
      # 界面（front_step 流体）从地面注入的时刻：front_step.cumulative_volume_start_m3
      # 用 _front_arrival_time 同款体积→时刻反解在未截断全序列上查（注意不是管容偏移）。
      # 特例：首步流体 t_inject = scheduled_steps_full[0].start_time_s（≈0）。
  _contact_time_integrated_sigma(t_arrival, t_inject, t_travel, D_eff, U, dt) -> float
      # t_contact = max(t_arrival - t_inject, 0.0)   # 防负值（尾浆滞留井 t_arrival 可能被截断）
      # sigma_t = sqrt(2·D_eff·t_contact)/U；防御上限 min(σ, 0.5·t_travel)、下限 max(σ, dt) 与现状一致
  ```
- [ ] `_apply_dispersion_to_timeline` 内：`mixing_contact_time=True` 时 σ_t 改由 `_contact_time_integrated_sigma` 给出（D_eff/混浆增强/防御上下限全部不动）；`t_travel` 仍计算（防御上限用）。开关关 → 走原路径逐位不变。
- [ ] 跑 Task 1 测试 + 全量回归。commit：`feat(transport1d): 接触时间积分σ_t（mixing_contact_time 默认关，逐位复现）`。

## Task 2：胶塞面零掺混（S2，~0.5 天）

- [ ] **测试先行**：
  4. `test_plug_face_no_dispersion`：ht1_003（碰压成功、loader 生产口径 has_plug=False——单测直接构造 `solver.has_plug=True` + 开关开），断言"尾浆→压塞液"界面（prev=尾浆名, next=压塞液名）**无过渡子事件**（只有原阶跃事件）、其他界面过渡带不受影响；
  5. `test_plug_face_off_bitwise`：开关关 → ht1_003 时间线逐位等于基线；
  6. `test_hu102_unaffected`：hu102（唯一生产 has_plug=True 井、未碰压、替浆步被截断）开关开 → 输出不受本开关影响（保护性断言：其界面集无"尾浆→压塞液"异物对）。
- [ ] **实现**：`__init__` 增 `plug_face_zero_mixing: bool = False`；`_apply_dispersion_to_timeline` 生成过渡带前判据：
  ```
  if (self.plug_face_zero_mixing and self.has_plug
          and self._is_plug_release_fluid(fluid_name, fluids)):
      dispersed_events.append(event); continue   # 机械隔离界面，零掺混
  ```
  `_is_plug_release_fluid`：流体 role==`FluidRole.OTHER` **且名含"压塞液"**（各 loader 现状=OTHER+名字，不新增角色枚举、不动 FluidSpec——改动最小且向后兼容；docstring 注明这是权宜判据，若未来引入 PLUG 角色则切换）。
- [ ] 跑 Task 2 测试 + 全量回归。commit：`feat(transport1d): 胶塞面零掺混（plug_face_zero_mixing 默认关）`。

## Task 3：8 井双开关验证（S4，~1 天；可与前置消融合并）

- [ ] 新建 `scripts/rerun8_mixing_contact_time_20260908.py`（复用 `scripts/rerun8_three_fixes_20260906.py` 的**纯默认 solver 口径**——勿用 `_sensitivity_common.py` 的 CORRECTED_KW 骨架，口径不同）：
  - 变体矩阵：8 井 × {基线(双关), contact_time 开, plug_face 开, 双开}；
  - 输出 `results/浓度剖面真解路线B_2026-09-08/`：每井 η_E/η_N/混浆指数/stop/窄边前缘 对照 CSV + 汇总表。
- [ ] **验收核对**（红线）：双开 vs 基线——η_E 全井 |Δ|<0.5pp、stop 逐位不变。预期效应集中在混浆指数（下降）与 hu102（过渡带占评价域 13.5% 最敏感井，η_E 预期 <0.3pp 变化）。
- [ ] （可选前置，建议先做）**alpha 敏感性消融**：hu103/hu102 双井 × alpha∈{0.04,0.12,0.25,0.50}（基线开关口径），量化"混浆带宽→η"灵敏度上界——若 alpha 减半 η 都不动 <0.2pp，则路线 B 的 η 效应预期同样微小，论文叙事按"混浆带物理化修正"而非"结果改变"写。
- [ ] obsidian 调研报告追加"§七 实施结果"（对照表 + 验收结论 + 遗留）。commit：`feat(results): 路线B 8井双开关验证落盘`。

## Risk 与回退

| 风险 | 概率 | 缓解 |
|------|------|------|
| t_inject 反解在分段注入/RESTART 截断井取错步 | 中 | 测试 3 覆盖；实现沿 prev_fluid 查找逻辑取"前缘流体最后一段注入步"；hu102/hu101（尾浆滞留井）t_arrival 可能被 pumping_end 截断——t_contact 用 `max(t_arrival - t_inject, 0.0)` 防负值 |
| 8 井 η 扰动超 0.5pp | 低 | 双开关关=逐位回退；退回做 alpha 消融定位灵敏度，向用户报告后裁定 |
| σ_t 收窄后过渡带 <2 子事件 | 低 | 保留 n_sub=5 与 F4 收尾结构不动，只变 σ_t 数值；σ_t 下限 max(σ,dt) 防退化 |
| 胶塞面判据误伤（未来 loader 改名） | 低 | `_is_plug_release_fluid` 单测 + docstring 权宜声明；role=OTHER 是稳定现状（5 loader 一致） |

## 路线 A 过渡预留（二期，不在本计划）

- 接触时间积分是 A（c(z,t) PDE）中"界面追踪 + D_eff(z,t) 输入"的降维形式；B 的 D_eff 闭包函数、混浆增强因子、胶塞面语义在 A 中全部复用；
- A 开工信号：出现"多界面掺混带重叠"井（8 井现况探针实测无重叠——最薄隔离液段 10 m³ >> 相邻带 0.5 m³）或论文审稿要求管内剖面真解；
- 本计划双开关届时退役为 A 的边界条件校验器（B 的 σ_t 解析解 vs A 的网格解一致性检查）。
