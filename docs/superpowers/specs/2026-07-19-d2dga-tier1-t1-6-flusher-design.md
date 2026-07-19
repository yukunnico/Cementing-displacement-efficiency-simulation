# D2DGA Tier 1 T1-6 FLUSHER 流体角色 + 独立浓度场设计

> **日期**：2026-07-19
> **范围**：Tier 1 剩余首项 T1-6（FLUSHER 独立浓度场，被动平流相）
> **关联**：`2026-07-18-d2dga-tier1-design.md`（T1-1~T1-5 已完成）、`docs/顶替效率模型参考文献调研_改进方向.md`、Yang 2021
> **执行模式**：主 agent 统筹，子 agent 执行（spec §6 含 Task 分解，作 plan 用）
> **验证约束**：⚠️ 不对标 CBL，机制就绪验证

---

## 0. 背景与设计约束

### 0.1 前置状态

- **T1-1~T1-5 已完成**（commits `5cb03b0..44c0868`，236 测试通过）：去通量限幅、I3 物理化、体力向量注入、两层黏度闭包、static wall c_min。D2DGA 物理闭包已物理化。
- **T1-6 是 Tier 1 剩余首项**：加 FLUSHER（冲洗液，Yang 2021 mud-spacer-flusher-cement 序列 cement 前独立相）独立浓度场。

### 0.2 WASH vs FLUSHER 语义（关键澄清）

- **WASH = 清洗液**（`fluid_spec.py:16` 注释），代码映射→spacer 相。**保留不动**。
- **FLUSHER = 冲洗液**（Yang 2021），mud-spacer-flusher-cement 序列中 cement 前的独立相（Case 7 = 89.67%，比三相 +3.79%；Case 8 flusher 在 spacer 前 = -4.41%）。**新增独立相**。
- 现有 loader 把中文"冲洗液/先导浆"塞 WASH 是**语义错配**（把冲洗液当清洗液）。T1-6 逐井重映射。

### 0.3 设计原则

1. **最小独立浓度场**：flusher 作被动平流相，独立 `flusher_field`，**不参与** D2DGA 通量放大/I3/体力注入/两层闭包（T1-1~T1-5 成果不动）。
2. **五相体积闭合**：`mud = 1 - lead - tail - spacer - flusher`。
3. **后向兼容**：`AnnulusSimulationResult` 新字段默认 None/空，下游脚本/测试不破坏。
4. **机制就绪验证**（不对标 CBL）：flusher 独立追踪 + 五相闭合 + 不崩溃；Yang Case 7/8 仅"序列可表达"检查，不要求复现 +3.79%；为 T2-5 序列优化器铺路。
5. **移植而非原创**：FLUSHER 角色与序列来自 Yang 2021 已发表。

---

## 1. 改动总览

| 层 | 文件 | 改动 | 风险 |
|----|------|------|------|
| L1 | `cemdisp/data/fluid_spec.py` | `FluidRole` 加 `FLUSHER="flusher"` | 低 |
| L3 | `cemdisp/models2d/boundary_bridge.py` | `_phase_fractions_for_fluid` 加 `FLUSHER→"flusher"` 映射 | 低 |
| L4 | `cemdisp/models2d/annulus_d2dga.py` | `AnnulusSimulationResult` 加 `flusher_field/snapshots` + `run` 推进 + 五相闭合 + 入口 + 输出 | 中 |
| L2 | `cemdisp/transport1d/casing_flow.py` | 无改动（flusher 作 FluidSpec 自动追踪） | — |
| loaders | `cemdisp/data/loaders/hu102_loader.py` / `ht1_004_loader.py` | 逐井重映射"冲洗液/先导浆"→FLUSHER | 中 |

依赖顺序：L1 fluid_spec → L3 boundary_bridge → L4 annulus_d2dga → loaders。L2 无改动。

---

## 2. 各层详细设计

### 2.1 L1 `fluid_spec.py`：`FluidRole.FLUSHER`

- `FluidRole` 枚举（行 47-58）加 `FLUSHER = "flusher"`，注释"冲洗液（Yang 2021 mud-spacer-flusher-cement 序列 cement 前独立相）"。
- `FluidSpec` 不加专属字段——`role=FluidRole.FLUSHER` 即标识冲洗液，物性复用现有 `density/rheology_model/...` 字段。
- WASH 注释保持"清洗液"不动。

### 2.2 L3 `boundary_bridge.py`：FLUSHER→"flusher" 相名映射

- `_phase_fractions_for_fluid`（行 70-93）加分支：`role==FLUSHER → "flusher"`（独立相名，**不并入 spacer**）。
- `AnnulusInletState.phase_fractions`（行 45-57，`tuple[tuple[str, float], ...]`）支持 `"flusher"` 名。
- `_phase_fractions_from_state`（行 96-112）合并 flusher 分数（与其他相并列，不再合并到 spacer）。
- 入口相分数仍**均匀施加**（无方位加权，T1-9 另做）。

### 2.3 L4 `annulus_d2dga.py`：flusher 独立浓度场 + 五相闭合

**`AnnulusSimulationResult`（行 116-152）**：加 `flusher_field: Array` + `flusher_snapshots: tuple`，默认 None/空（后向兼容）。

**`run` 循环**：
- 初始化 `flusher = np.zeros((ny, nz))`（与 lead/tail/spacer 同构，行 689-691 区）。
- 平流推进：与 lead/tail/spacer 同构（行 757-762 区），用 `inlet_flusher_fraction`。
- 入口加 `inlet_flusher_fraction = _phase_fraction(inlet_state, "flusher")`（行 716-719 区）。
- 快照保存：`flusher_snapshots` 与 cement/spacer 同构（行 814 区）。

**`_compute_props`（行 425-481）**：体积闭合改 `mud = 1 - lead - tail - spacer - flusher`（行 438-440），`np.clip(mud, 0, 1)` 守卫。**flusher 不参与** `c_bar=lead+tail`、`m_field`、`eta1/eta2`、体力、I3、两层闭包（被动相）。

**`_depth_profiles`（行 650-677）**：加 `"冲洗液平均浓度": np.average(flusher, axis=0, weights=geom["b"])`。

**`metrics`（行 878-915）**：加 `mean_flusher`。

**`_compute_velocity`**：flusher 不影响流动度（被动，不进 pref/base）。

### 2.4 L2 `casing_flow.py`：无改动

flusher 作为 `FluidSpec`（role=FLUSHER）加入 `PumpingSchedule`，`CasingFlowSolver` 按流体名自动追踪（行 478 `phase_fractions = ((fluid_name, 1.0),)`）。无需改动。

### 2.5 loaders：逐井重映射

- `hu102_loader.py` / `ht1_004_loader.py`：现有"冲洗液/先导浆"（WASH）按物理判断：
  - mud-spacer-flusher-cement 序列中 **cement 前的冲洗液** → `role=FluidRole.FLUSHER`，加 flusher `FluidSpec` + schedule step。
  - 套管清洗液 → 留 `WASH`。
- 若某井无真实冲洗液，加合成 FLUSHER 流体验证序列可表达。

---

## 3. 集成与数据流

- **`AnnulusSimulationResult` 后向兼容**：`flusher_field`/`flusher_snapshots` 默认 None/空，下游脚本（`quality_proxy`/`tier0_diagnostics`/ablation/plot）不引用这两个字段则不破坏；引用处加 `getattr` 守卫。
- **五相闭合同步**：所有 `mud = 1 - ...` 处（`_compute_props` 行 438-440、`_depth_profiles` 行 650-677 的 mud 列）同步改五相。
- **flusher 被动**：不进 `c_bar`/`m_field`/`eta1/eta2`/体力/I3/两层闭包/static wall。D2DGA 物理闭包（T1-1~T1-5）零改动。
- **Tier0 诊断复跑**：flusher 不影响 muskat_regime/浮力分类/η_N（因不参与通量/体力），诊断结果对 flusher 引入应不变（仅多 mean_flusher 指标）。
- **入口相名**：`run` 内 `"cement"/"lead"/"tail"/"spacer"` 硬编码处（行 716-719）加 `"flusher"`。

---

## 4. 测试策略

| 层 | 测试 | 文件 |
|----|------|------|
| L1 | `FluidRole.FLUSHER` 枚举值；`FluidSpec(role=FLUSHER)` 构造不报错 | `tests/test_fluid_spec.py`（或现有 well_spec 测试） |
| L3 | `FLUSHER→"flusher"` 映射；`AnnulusInletState` 含 flusher 分数；不并入 spacer | `tests/test_boundary_bridge.py` |
| L4 | flusher 独立追踪（注入 flusher 后 flusher_field 增、其他相不混入）；五相闭合 sum(lead+tail+spacer+flusher+mud)=1；平流推进；`AnnulusSimulationResult` 默认 None 后向兼容 | `tests/test_improved_d2dga_annulus.py` |
| 序列 | 合成 mud-spacer-flusher-cement 序列可表达（Yang Case 7 结构）；不要求复现 +3.79% | `tests/test_improved_d2dga_annulus.py` |
| 回归 | R0-R3 消融不崩溃 + 六井集成（不对标 CBL） | `tests/test_improved_d2dga_annulus.py` + `test_six_well_integration.py` |
| 诊断 | hu102/ht1_004 实跑 + Tier0 诊断复跑（mean_flusher 新增，其他不变） | `scripts/run_tier0_diagnostics.py` |

**R0-R3 回归基线**：T1-6 不改 D2DGA 物理闭包，R0-R3 效率值应与 T1-5 后一致（flusher 被动，不影响 mud-cement 顶替）；若 loader 重映射引入 flusher 流体，则实跑结果含 flusher 相。

---

## 5. 验证门（无 CBL，机制就绪）

| 验证类型 | 内容 |
|---------|------|
| 机制就绪 | flusher 独立追踪（不与 mud/cement/spacer 混）、五相闭合 sum=1、平流推进、不崩溃 |
| 序列可表达 | 合成 mud-spacer-flusher-cement（Yang Case 7 结构）可构造+运行；Case 8（flusher 在 spacer 前）亦可；不要求复现 +3.79% |
| 回归 | R0-R3 消融 + 六井集成，不崩溃+物理合理 |
| 诊断复跑 | hu102/ht1_004 Tier0 诊断：mean_flusher 新增，muskat/浮力/η_N 不变 |

---

## 6. 执行计划（子 agent 分工 + Task 分解，作 plan 用）

### 6.1 子 agent 分工

| 步骤 | 子 agent | 任务 | 产出 |
|------|---------|------|------|
| SA1 | L1+L3 | `fluid_spec.py` FLUSHER + `boundary_bridge.py` 映射 + 测试 | 2 模块 + 测试 |
| SA2 | L4 | `annulus_d2dga.py` flusher_field + run 推进 + 五相闭合 + 入口 + 输出 + 测试 | 模块 + 测试 |
| SA3 | loaders | hu102/ht1_004 逐井重映射 + 测试 | 2 loader + 测试 |
| 主 agent | 统筹 | 审查 + 全量测试 + R0-R3/六井/Tier0 诊断/Yang 序列验证 | — |

### 6.2 Task 分解（bite-sized TDD，按 Task 1 模板：写失败测试→跑确认 FAIL→实现→跑确认 PASS→commit）

- **Task 1**（L1）：`fluid_spec.py` `FluidRole.FLUSHER="flusher"` + 注释；测试 `FluidRole.FLUSHER.value=="flusher"` + `FluidSpec(role=FLUSHER)` 构造。
- **Task 2**（L3，依赖 Task1）：`boundary_bridge.py` `_phase_fractions_for_fluid` 加 `FLUSHER→"flusher"` + `_phase_fractions_from_state` 合并；测试 FLUSHER 映射独立、不并入 spacer、`AnnulusInletState` 含 flusher。
- **Task 3**（L4-a，依赖 Task2）：`AnnulusSimulationResult` 加 `flusher_field/flusher_snapshots`（默认 None/空）；`run` 初始化 `flusher=np.zeros` + 平流推进（用 `inlet_flusher_fraction`）+ 入口加 `inlet_flusher_fraction` + 快照；测试 flusher 独立追踪 + 平流推进 + 后向兼容（None 默认）。
- **Task 4**（L4-b，依赖 Task3）：`_compute_props` 体积闭合改五相 `mud=1-lead-tail-spacer-flusher` + clip 守卫；`_depth_profiles` 加 flusher 平均浓度；`metrics` 加 `mean_flusher`；`run` 内 `1-...` 处同步五相；测试五相闭合 sum=1 + flusher 不混入 mud/cement。
- **Task 5**（loaders，依赖 Task1-4）：hu102/ht1_004 逐井重映射"冲洗液/先导浆"→FLUSHER（按物理判断）+ flusher FluidSpec + schedule step；测试 loader 构造含 FLUSHER。
- **Task 6**（集成验证）：全量 `pytest tests/` + R0-R3 消融（不崩溃）+ 六井集成 + Tier0 诊断复跑（hu102/ht1_004，mean_flusher 新增）+ 合成 mud-spacer-flusher-cement 序列可表达（Yang Case 7 结构）+ commit。

### 6.3 子 agent 统一约束

- 中文 docstring + 类型注解 + 公式/来源标注。
- 数值安全：`np.clip` 守卫五相闭合（sum 不超 1）、NaN 守卫。
- 后向兼容：`flusher_field` 默认 None，下游引用处 `getattr` 守卫。
- 环境：`conda activate shenjingwangluo`；`/d/apps/Anaconda/python` 加 `PYTHONIOENCODING=utf-8 PYTHONUTF8=1`；`pip install -e .`。
- TDD：每 Task 先写失败测试，跑确认 FAIL，再实现，跑确认 PASS，commit。

---

## 7. 风险与缓解

| 风险 | 缓解 |
|------|------|
| 体积闭合五相未同步（`1-...` 遗漏）→ 负泥浆或越界 | Task 4 集中改所有 `1-...` 处 + `np.clip(mud,0,1)` 守卫 + sum=1 测试 |
| 入口相名硬编码缺 `"flusher"` → flusher 被丢弃 | Task 3 入口加 `inlet_flusher_fraction` + 测试 |
| `AnnulusSimulationResult` 新字段破坏下游 | 默认 None/空 + 下游 `getattr` 守卫 + 后向兼容测试 |
| WASH vs FLUSHER loader 逐井误判 | Task 5 按物理判断（mud-spacer-flusher-cement 序列 cement 前→FLUSHER；套管清洗→WASH）+ 测试 |
| flusher 被动无法复现 Yang +3.79% | 已接受（机制就绪验证，为 T2-5 铺路） |

---

## 8. 后续衔接

- T1-6 完成 → FLUSHER 独立相就绪，为 **T2-5 泵注序列自动优化器**（Yang 2021）提供相基础。
- Tier 1 剩余：T1-7 CFL 自适应 / T1-8 1D Taylor 弥散 / T1-9 入口方位加权（各自子 spec）。
- 若后续需 flusher 物理参与（复现 Yang 效率），扩展范围 B/C（flusher 进黏度比/浮力），单独 spec。


