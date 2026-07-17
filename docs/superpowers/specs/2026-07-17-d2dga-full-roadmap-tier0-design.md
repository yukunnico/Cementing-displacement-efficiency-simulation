# D2DGA 模型全路线图改进 — 主路线图 + Tier 0 详细设计

> **日期**：2026-07-17
> **类型**：模型驱动（独立改进路线），非论文驱动
> **范围**：全路线图（Tier 0/1/2）分解 + Tier 0 详细设计
> **关联文档**：`docs/顶替效率模型参考文献调研_改进方向.md`（差距分析源头）、`docs/模型求解流程与代码结构详解报告.md`（代码地图）、Obsidian `改进措施文献调研表.md`（前期文献映射）
> **执行模式**：主 agent 统筹，子 agent 执行
> **验证约束**：⚠️ **暂不对标 CBL**。验证以理论对标为主（Zhang 2023 Table 3 / Bararpour 2025 regime / Yang 2021 序列效率）。`cbl_comparison.py` 与 `quality_proxy.py` 在本路线图验证门中**不使用**。

---

## 0. 背景与设计约束

用户明确三点：
1. **模型驱动**：目标是提升模型物理保真度与诊断可信度本身，不为特定论文服务。
2. **全层改进**：5 层（数据/1D输运/鞋口桥接/2D核心/诊断报告）的代码都要改进，不只是 2D 核心。
3. **暂不对标 CBL**：验证改用理论对标。

设计原则：
- **移植而非原创**：所有改进都有 Frigaard 课题组已发表、已同行评议的公式支撑（见调研报告），技术风险低。
- **增量验证**：每个 Tier 完成后过验证门，再进下一 Tier。
- **YAGNI**：Tier 2 细节不在本 spec 写死，到达时再各自出详细 spec。

---

## 1. 全路线图分解（3 阶段 × 5 层）

### 1.1 阶段总览

| 阶段 | 名称 | 覆盖层 | 风险 | 预估 | 验证门 |
|------|------|--------|------|------|--------|
| **Tier 0** | 诊断可信度提升 | Layer 5（诊断）| 低（纯后处理，不改求解器）| ~2 周 | 理论对标 + 单元测试 |
| **Tier 1** | 核心物理修正 | Layer 4（2D核心）+ Layer 2（1D）+ Layer 1/3（FLUSHER）| 中（改速度场/通量）| ~3-4 周 | R0-R3 消融回归 + 理论对标 |
| **Tier 2** | 根本性重构 | Layer 4 + 新优化模块 | 高（重构求解器）| 长期 | 理论对标 + 3D 基准（未来）|

### 1.2 Tier 0 改进项（本 spec 详细设计，见 §2）

| # | 改进项 | 来源 | 层 |
|---|--------|------|-----|
| T0-1 | 流动分类判据（steady/unsteady × dispersive/non-dispersive）| Zhang 2023 式 3.1-3.6 | L5 |
| T0-2 | Muskat 三 regime（stable/partial/unstable）| Bararpour 2025 式 3.8-3.13 / Pelipenko 2004 式 3.33-3.34 | L5 |
| T0-3 | 浮力数 b 设计阈值分类 | Zhang 2023 §3.1.3 | L5 |
| T0-4 | 泥浆滞留体积 + φm/φc 污染指示器 | Yang 2021 | L5 |
| T0-5 | 界面长度比 + 窄边效率 η_N | Yang 2021 / Zhang 2023 图 19 | L5 |
| T0-6 | 停泵有限时间衰减诊断 | Moyers-González 2007 式 3.35/3.40 | L5 |
| T0-7 | 突破时间 t_br | Zhang 2022 Table 3 | L5 |

### 1.3 Tier 1 改进项（路线图级，到达时出详细 spec）

| # | 改进项 | 来源 | 层 |
|---|--------|------|-----|
| T1-1 | 移交通量放大人工限幅 `[0.5,2.0]` | Zhang 2022 式 4.28 | L4 `d2dga_flux.py:51-52,78` |
| T1-2 | I3 通量系数物理化（`flux_strength=0.05`→`ΔρH³/(6η₂)`）| Zhang 2022 式 4.25 | L4 `annulus_d2dga.py:775` |
| T1-3 | 激活浮力体力向量（含 I2/I1 分层修正，式 4.24）| Zhang 2022 式 4.13/4.21/4.24 | L4 `annulus_d2dga.py:466-487,587-604` |
| T1-4 | 壁面剪切两层闭包 `1/η_mix=c̄³/η₂+(1-c̄³)/η₁` | Zhang 2022 式 4.23 | L4 `annulus_d2dga.py:566-570` |
| T1-5 | static wall layer 判据激活 wall 场（c_min 判据）| Bararpour 2025 式 2.35-2.41 | L4 `annulus_d2dga.py:668,923` |
| T1-6 | FLUSHER 流体角色 + 独立浓度场 | Yang 2021 | L1 `fluid_spec.py:47` + L3 `boundary_bridge.py:91` + L4 |
| T1-7 | CFL 自适应时间步 | 数值稳定性 | L4 `annulus_d2dga.py:654-686` |
| T1-8 | 1D HB/Bingham 弥散替换 `α·U·R` 为屈服应力 Taylor 弥散 | E1 JNNFM 2016（见文献调研表）| L2 `casing_flow.py:311-313` |
| T1-9 | 入口相分数按间隙宽度方位加权（当前均匀施加）| 调研报告 §2.3 | L3 `boundary_bridge.py` + L4 |

### 1.4 Tier 2 改进项（路线图级，长期）

| # | 改进项 | 来源 | 层 |
|---|--------|------|-----|
| T2-1 | HB D2DGA 闭包查表（I1/I2/I3/q0 数值积分）| Bararpour 2025 §2 | L4 `d2dga_flux.py` |
| T2-2 | 1D 间隙方向 augmented Lagrangian 求解器 | Bararpour 2025 Appendix A | L4 新建 `gap_solver.py` |
| T2-3 | TVD 激波捕捉替换拉普拉斯平滑 | Bararpour 2025 / Zhang 2022 FCT | L4 `annulus_d2dga.py:380-410` |
| T2-4 | 瞬态流函数演化方程（捕捉脉动）| Moyers-González 2007 式 2.21 | L4 速度场求解 |
| T2-5 | 泵注序列自动优化器 | Yang 2021 | 新建 `cemdisp/optimization/` |

### 1.5 依赖关系

```
T0（全部独立，纯 L5 诊断，可并行）
  └─ 解锁：诊断可信度基线，为 T1/T2 提供验证工具
T1-6（FLUSHER）─→ Yang 2021 序列验证 + T2-5（序列优化器）
T2-2（augmented Lagrangian）─→ T2-1（HB 闭包查表）
T1-3（浮力体力向量）─→ 水平井场景重新验证
T1-1/T1-2（去限幅）─→ 应在 T1-3/T1-4 之前（先恢复物理通量形式）
```

### 1.6 验证策略（无 CBL）

| 验证类型 | 内容 | 用于 |
|---------|------|------|
| **理论对标** | Zhang 2023 Table 3（10 算例 t_br/η_E）、Bararpour 2025 regime 边界（牛顿 m>1.5 失稳、幂律 b~70 稳定）、Yang 2021 序列效率（89.67% vs 85.88%）| Tier 0/1/2 |
| **单元测试** | 每个诊断/闭包对解析值或已知极限 | 全部 |
| **回归测试** | R0-R3 消融（`test_improved_d2dga_annulus.py`）+ 六井集成（`test_six_well_integration.py`），仅检查不崩溃 + 物理合理，**不对标 CBL** | Tier 1/2 |
| **通量纯函数测试** | `test_d2dga_flux.py` | Tier 1 改通量时 |

---

## 2. Tier 0 详细设计

### 2.0 核心架构决策：纯 Layer 5 诊断，零求解器改动

**关键发现**：`AnnulusSimulationResult` 已保存 `cement_snapshots`（tuple of ny×nz）、`snapshot_times_s`、`geom['b']`、`lead/tail/spacer_snapshots`、`metrics`、`depth_profiles`、`summary`。Tier 0 全部 7 项可作为**纯诊断函数**消费现有结果，**无需改 `annulus_d2dga.py`**。

流动分类所需的截面平均浓度场 `c̄(s,t)` 在诊断层从快照重建：
```python
c_bar_st = np.stack([np.average(snap, axis=0, weights=geom['b']) for snap in result.cement_snapshots], axis=1)  # (nz, n_snapshots)
```

**模块组织**：新建 4 个诊断模块于 `cemdisp/diagnostics/`，遵循现有 `quality_proxy.py` 模式（frozen dataclass 结果 + 纯计算函数 + 不 import models2d 核心，只消费 `AnnulusSimulationResult`）。

| 模块 | 覆盖项 | 结果 dataclass |
|------|--------|---------------|
| `flow_classification.py` | T0-1 | `FlowClassificationResult` |
| `muskat_regime.py` | T0-2 | `MuskatRegimeResult` |
| `displacement_metrics.py` | T0-4, T0-5, T0-7 | `DisplacementMetricsResult` |
| `regime_classifiers.py` | T0-3, T0-6 | `BuoyancyRegimeResult`, `ShutdownDecayResult` |

### 2.1 T0-1 流动分类判据（`flow_classification.py`）

**物理**：Zhang 2023 用相似变换 `c̄(ẑ/t̂)` 判断对流主导，定义归一化前缘速度 `w_f(c̄)=(ẑ/t̂)/ŵ₀`。

**计算**：
- `Δw_f = w_f(c̄=0.3) - w_f(c̄=0.7)`（式 3.1）；`Δw_f ≤ 0.1` → steady front
- 相对速度 `w_r(c̄)=w_f(c̄)-1`，分正负部 `w_r±`
- `σ_{w_r+}=sqrt(mean(w_r+²))`（式 3.2）；`|w̄_{r+}|=∫w_r+dc̄`（式 3.4）
- dispersive 需同时 `σ_{w_r+}>0.08` 且 `|w̄_{r+}|>0.05`（式 3.6a,b）
- 残余流体 `|w̄_{r-}|=-∫w_r-dc̄`（式 3.5），随偏心度增大

**输出** `FlowClassificationResult`：`delta_w_f`, `sigma_wr_plus`, `abs_wr_plus`, `abs_wr_minus`, `is_steady`, `is_dispersive`, `flow_class`（三态：unsteady_dispersive / dispersive_steady / non_dispersive_steady）, `c_bar_st`, `w0`。

**测试**：合成 `c̄(s,t)` 场——(a) 活塞流（steady non-dispersive，w_f≡1）；(b) 线性扩散前沿（dispersive）。验证判据正确分类。

### 2.2 T0-2 Muskat 三 regime（`muskat_regime.py`）

**物理**：Bararpour 2025 式 3.8-3.13。指进速度 vs 前缘速度比较。

**计算**：
- `w_finger(c̄) = I1(1)/I1(c̄) + b·I1(1)·[I2(c̄)/I1(c̄)+c̄-1]`（式 3.12）
- `w_f(c̄) = q0'(c̄) + b·I3'(c̄)`（式 3.8，含激波 Rankine-Hugoniot 式 3.9）
- `Δw(c̄)=w_finger-w_f`（式 3.13）
- 判据：全正→unstable；存在 c_critical 使 [0,c_critical] 为负且 c*<c_critical→stable；否则 partial_penetration
- 牛顿极限：I1/I2/I3/q0 用解析式（`d2dga_flux.py` 已有 q0 放大因子与 I3）；A_w/A_n 由 Pelipenko 2004 式 3.33-3.34 补充宽/窄边判据

**输出** `MuskatRegimeResult`：`regime`（stable/partial_penetration/unstable）, `c_critical`, `c_star`, `delta_w_profile`, `wide_side_unstable`, `narrow_side_unstable`。

**测试**：Bararpour 2025 已知边界——牛顿等密度 `m>1.5`→unstable；`b=10`→partial；`b=100`→stable；幂律 `b~70`→stable。

### 2.3 T0-3 浮力数阈值分类（`regime_classifiers.py`）

**物理**：Zhang 2023 §3.1.3。浮力数 b 是垂直井主导参数。

**判据**：
- `b < 0`：密度不稳定，**禁止**（forbidden）
- `0 ≤ b < 20`：高度弥散（highly_dispersive）
- `20 ≤ b < 80`：steady front 可达（steady_capable）
- `b ≥ 80~100`：non-dispersive（non_dispersive）

**输出** `BuoyancyRegimeResult`：`b_number`, `regime`, `design_advice`（中文建议字符串）。

**测试**：边界值 -5/10/50/100 分类正确。b 从 `result.summary['buoyancy_number']` 读（已算，`annulus_d2dga.py:892`）。

### 2.4 T0-4 泥浆滞留 + φm/φc（`displacement_metrics.py`）

**物理**：Yang 2021。偏心导致窄边泥浆滞留，φm/φc 指示胶结质量。

**计算**：
- `mud_retention_fraction = ∫b·mud dA / ∫b dA`（最终时刻）
- `phi_m_phi_c = mud / max(cement, ε)` 逐深度，按 0.05/0.1 分红/黄/绿
- 提升到 `summary` + `depth_profiles` 加列

**输出**：并入 `DisplacementMetricsResult`。

### 2.5 T0-5 界面长度比 + η_N（`displacement_metrics.py`）

**计算**：
- `interface_length_ratio = |front_wide - front_narrow| / max(mean(front_wide, front_narrow), ε)`（Yang 2021 核心观察量）
- `η_N`：窄四分位（azimuth 最后 1/4）的 `∫b·cement/∫b`（Zhang 2023 图 19）

### 2.6 T0-6 停泵有限时间衰减（`regime_classifiers.py`）

**物理**：Moyers-González 2007 式 3.35/3.40。停泵后屈服应力+小密度差→速度有限时间归零→界面冻结。

**判据**：`τ_{Y,min}/(1+e) ≥ ||div(f)||_∞/2`（式 3.35）满足则输出冻结时间 `t_s`（式 3.40）。

**触发**：`PumpingSchedule` 有 SHUTDOWN event_tag 时。

**输出** `ShutdownDecayResult`：`condition_satisfied`, `freeze_time_s`, `physical_interpretation`。

**测试**：解析值——构造满足/不满足条件的参数，验证判据与 t_s 量级。

### 2.7 T0-7 突破时间 t_br（`displacement_metrics.py`）

**计算**：`t_br` = 水泥前缘（c̄≥0.5）首次到达出口 `s_max` 的时间。从 `metrics` 的 `front_wide/front_narrow` 或 `c_bar_st` 提取。与 Zhang 2022 Table 3 的 t_br 范围（0.33-0.99 无量纲）对照需无量纲化 `t̂_br=t_br·ŵ₀/L`。

### 2.8 集成与报告

- 各诊断结果汇入 `result.summary`（新增 `flow_classification`/`muskat_regime`/`buoyancy_regime`/`displacement_metrics` 子键），**不改 `AnnulusSimulationResult` 结构**（用 summary dict 扩展，避免破坏现有 schema 测试）。
- 新增诊断图（可选，Tier 0 后期）：流动分类 `c̄(ẑ/t̂)` 坍缩图、Muskat `Δw(c̄)` 剖面图。
- runner 接线：在 `run_and_export` 中调用诊断聚合函数 `compute_all_tier0_diagnostics(result, well_spec, fluids, schedule)`，输出 JSON + 并入 markdown 摘要。

---

## 3. Tier 0 执行计划（子 agent 分工）

### 3.1 执行结构

由于 Tier 0 全部纯 Layer 5、无求解器改动、模块间独立，采用**并行子 agent**：

| 步骤 | 子 agent | 任务 | 产出 |
|------|---------|------|------|
| **SA1** | 诊断-流动分类 | `flow_classification.py`（T0-1）+ 单元测试 | 模块 + `tests/test_flow_classification.py` |
| **SA2** | 诊断-Muskat | `muskat_regime.py`（T0-2）+ 单元测试 | 模块 + `tests/test_muskat_regime.py` |
| **SA3** | 诊断-顶替指标 | `displacement_metrics.py`（T0-4/5/7）+ 单元测试 | 模块 + `tests/test_displacement_metrics.py` |
| **SA4** | 诊断-regime分类 | `regime_classifiers.py`（T0-3/6）+ 单元测试 | 模块 + `tests/test_regime_classifiers.py` |
| **SA5** | 集成 | 聚合函数 + summary 接线 + runner 接线 + `__init__.py` 导出 | 集成代码 |
| **主agent** | 统筹 | 审查各 SA 产出、跑全部测试、解决冲突、验证 | — |

### 3.2 执行顺序

1. **并行启动 SA1-SA4**（4 个诊断模块，相互独立，无依赖）。
2. 主 agent 审查 SA1-SA4 产出，跑各自单元测试。
3. **SA5 集成**（依赖 SA1-SA4 的模块接口，需等它们完成）。
4. 主 agent 跑全量测试（`pytest tests/`）+ 在 hu102/ht1_004 上实跑验证（仅检查诊断产出物理合理，**不对标 CBL**）。

### 3.3 子 agent 统一约束（写进每个 SA 的 prompt）

- 遵循现有 `cemdisp/diagnostics/quality_proxy.py` 模式：frozen dataclass 结果 + 纯计算函数 + 类型注解 + docstring。
- **不 import `models2d` 核心**；只消费 `AnnulusSimulationResult`（从 `cemdisp.models2d.annulus_d2dga` import 类型用于注解即可，不调用求解）。
- 公式必须标注论文来源与式号（注释）。
- 每个模块配 `tests/test_*.py`，用合成数据测试已知极限。
- 数值安全：除零保护、`np.clip`、NaN 守卫。
- 中文 docstring，与现有代码风格一致。
- 环境：`conda activate shenjingwangluo`，`pip install -e .`，测试 `python -m pytest tests/test_<module>.py -v`。

### 3.4 验收标准（Tier 0 完成定义）

- [ ] 4 个诊断模块 + 4 个测试文件，全部 `pytest` 通过。
- [ ] 集成后 `pytest tests/` 全量通过（不破坏现有测试）。
- [ ] hu102 + ht1_004 实跑，诊断 JSON 产出，各指标物理合理（如 b>0、flow_class 非空、η_N≤η_E、t_br>0）。
- [ ] Muskat regime 对牛顿极限 `m>1.5`→unstable 复现。
- [ ] 浮力数分类边界正确。
- [ ] summary 新增诊断子键，markdown 摘要含新指标。

---

## 4. 后续 Tier 衔接

Tier 0 完成后：
1. 诊断工具就绪，可作为 Tier 1/2 的验证手段（如用 Muskat regime 验证 T1-3 浮力体力向量激活后的稳定性预测改善）。
2. 出 **Tier 1 详细 spec**（重点：T1-1/T1-2 去限幅 → T1-3 激活浮力体力向量 → T1-5 static wall layer → T1-6 FLUSHER）。
3. Tier 1 需跑 R0-R3 消融回归 + 六井集成（不对标 CBL）。

---

## 5. 风险与缓解

| 风险 | 缓解 |
|------|------|
| 快照时间分辨率不足（save_interval=60×dt=4s=240s）影响流动分类精度 | 诊断层对 `c_bar_st` 做时间插值；必要时在 runner 调小 `save_interval`（配置改动，非求解器改动）|
| Muskat regime 的 I1/I2 对 HB 流体无数值积分（Tier 0 用牛顿解析式）| Tier 0 明确标注"牛顿近似"，HB 精确版留给 Tier 2（T2-1/T2-2）|
| 子 agent 接口不一致 | 主 agent 在 SA5 集成前统一定义 4 个结果 dataclass 的字段名（本 spec §2 已定）|
| 无量纲化约定（t̂、ŵ₀、m 定义）| 每个模块 docstring 写明所用约定，与 `d2dga_flux.py` 的 m=μ_mud/μ_cement 保持一致 |
| I3 分母 `2m(...)` vs Bararpour 式 2.27 `12(...)` 差异 | Muskat 模块内注明采用 `d2dga_flux.py` 现有 I3（与 Zhang 2022 式 4.26 一致），差异留待 Tier 1 核实 |
