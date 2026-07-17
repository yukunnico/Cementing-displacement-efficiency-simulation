# D2DGA Tier 1 详细设计 — 核心物理修正（T1-1~T1-5）

> **日期**：2026-07-18
> **范围**：Tier 1 的 L4 核心子集 T1-1~T1-5（不跨层）
> **关联**：`2026-07-17-d2dga-full-roadmap-tier0-design.md`、`docs/顶替效率模型参考文献调研_改进方向.md`、`docs/模型求解流程与代码结构详解报告.md`
> **执行模式**：主 agent 统筹，子 agent 执行
> **验证约束**：⚠️ 暂不对标 CBL，沿用 Tier 0 理论对标（Zhang 2023 Table 3 / Bararpour 2025 regime / Yang 2021 序列效率）

---

## 0. 背景与设计约束

### 0.1 前置状态

- **Tier 0 已完成**（commit `83c4a43`）：4 诊断模块 + 81 单元测试通过 + hu102/ht1_004 实跑诊断 JSON 产出。诊断工具就绪，可作 Tier 1 验证手段。
- **R0-R3 消融已实现**（commits `31066ee`/`ab70a3e`/`7f9506e`）：`enable_d2dga_auto_m`（R1，m 场）、`enable_d2dga_i3_flux`（R2，I3 弥散通量）、`enable_true_buoyancy`（R3，真浮力）三开关。

### 0.2 Tier 1 定位（关键澄清）

R0-R3 消融是 Tier 1 的**开关版前序**——做了 m 场、I3 公式、浮力方位角因子，但带三处人工限幅/简化：

| 位置 | 人工限幅/简化 | 物理本应 |
|------|--------------|---------|
| `d2dga_flux.py:78` | `np.clip(amplification, 0.5, 2.0)` | `f(c,m)` 自然取值 |
| `annulus_d2dga.py:775` | `flux_strength = 0.05` 硬缩放 | `ΔρH³/(6η₂)·I3` 物理系数直接驱动 |
| `annulus_d2dga.py:587-599` | `(2φ−1)` 密度对比简化代理 | 式 2.5b 体力向量注入流动度 + I1/I2 分层修正 |

Tier 1 = 在 R0-R3 开关版上**去掉这三处限幅/简化 + 补两层闭包 + 补 static wall**。二者不重复，Tier 1 是深化。

### 0.3 设计原则

1. **纯 L4 改动**：只动 `annulus_d2dga.py` + `d2dga_flux.py`，不跨层（T1-6/T1-8/T1-9 拆后续子 spec），不改 `AnnulusSimulationResult` 结构（用 summary 扩展）。
2. **直接物理化删限幅，不加新开关**：R0-R3 三开关保留原语义做消融对照；限幅版不可回退（已接受）。
3. **数值安全靠公式内保护，不靠开关**：`c_safe=clip(c,0.01,0.99)`（防除零，非物理限幅）+ 局部 CFL 裁剪（防单步越界，非全局限幅）。
4. **增量验证**：每项配单元测试，全部完成后过 R0-R3 消融回归 + 六井集成 + Tier 0 诊断复跑。
5. **移植而非原创**：所有公式有 Frigaard 课题组已发表支撑（见调研报告）。

---

## 1. 改进项总览

| # | 标题 | 来源式号 | 文件:行 | 现状 | 风险 |
|---|------|---------|---------|------|------|
| T1-1 | 去通量放大限幅 [0.5,2.0] | 式 4.28 | `d2dga_flux.py:78` | 人工限幅 | 低 |
| T1-2 | I3 物理系数接线（去 0.05）| 式 4.25 | `annulus_d2dga.py:775` | 开关版半成 | 中 |
| T1-3 | 体力向量注入流动度 + I1/I2 分层修正 | 式 2.5b/4.13/4.21/4.24 | `annulus_d2dga.py:466-487,587-604` | 开关版半成 | 高 |
| T1-4 | 两层黏度闭包 1/η_mix=c̄³/η₂+(1−c̄³)/η₁ | 式 4.23 | `annulus_d2dga.py:566-604,424-464` | 未实现 | 高 |
| T1-5 | static wall layer c_min 判据 | 式 2.35-2.41 | `annulus_d2dga.py:668,627,923` | 未实现（零场） | 中 |

依赖：T1-3 的 I1/I2 须先于体力注入；T1-4 的相黏度返回须先于两层闭包；T1-2/T1-3 共用 `_buoyancy_force_vector`。

---

## 2. 各项详细设计

### 2.1 T1-1 去通量放大限幅（`d2dga_flux.py:78`）

**现状**：`d2dga_flux_amplification` 末行 `amplification = np.clip(numerator/denominator, 0.5, 2.0)` 硬裁剪。R0-R3 的 `auto_m` 只改 m 输入，不触及此 clip。

**设计**：
- 删 `np.clip(..., 0.5, 2.0)`，让 `f(c,m)=[m·c²+1.5(1−c²)]/[m·c³+(1−c³)]` 自然取值。
- 保留 `c_safe=np.clip(c, 0.01, 0.99)`（数值防除零，非物理限幅）。
- 物理边界：c→0 趋 1.5、c→1 趋 1/m，无需硬裁剪。
- 保留 `min/max_amplification` 参数签名但默认设为 `±inf`（向后兼容 API，不再实际裁剪）。

**测试**（`tests/test_d2dga_flux.py`）：
- `f(0)=1.5`、`f(1)=1/m`（m=2 → 0.5）端点。
- 无 clip 边界突变：c 从 0.99→1.0 平滑过渡。
- m 场（数组）输入形状兼容。

### 2.2 T1-2 I3 物理系数接线（`annulus_d2dga.py:775`）

**现状**：`d2dga_buoyancy_flux`（`d2dga_flux.py:119`）已实现物理系数 `ΔρH³/(6η₂)·I3`，返回 `(q_phi, q_xi)`；但 `annulus_d2dga.py:775` 更新 `lead/tail` 时乘 `flux_strength=0.05` 限幅，物理系数算出却被盖掉。

**设计**：
- 删 `flux_strength = 0.05`（775），让 `lead/tail -= div_q · frac · dt` 直接用物理系数散度驱动（`div_q` 由 `q_phi/q_xi` 散度算出，763-772 现有逻辑保留）。
- `Δρ`/`η₂`/`H`/`m` 已局部化（749-762），无需改。
- **数值保护**（非物理限幅）：对 `div_q·dt` 加局部 CFL 裁剪 `|div_q|·dt ≤ α·Δs`（α=0.5，仅防单步越界，不做全局限幅）。
- `enable_d2dga_i3_flux` 开关语义不变（True=物理化 I3，False=关）。

**测试**：
- 构造已知 `Δρ/H/η₂` 解析工况，验证更新量量级 = 物理系数·梯度。
- 加"去 0.05 后通量幅度提升"回归断言（更新量比限幅版大）。

### 2.3 T1-3 体力向量注入流动度 + I1/I2 分层修正（最高风险）

**现状**：
- `_buoyancy_force_vector`（466-487）返回 `f_phi/f_xi`，只喂 R2 I3 通量（752），**未注入流动度**。
- `enable_true_buoyancy` 分支（587-599）实际用简化代理 `buoyancy_shape = 1 + clip(2·(ρ−ρ_displaced), −0.35, 0.45)·ebar·(2φ−1)`，未用完整体力向量。
- `_compute_velocity` 流动度链：`base = (b/mean(b))² / μ_reg`（570）→ 经 `buoyancy_shape` 修正得 `pref`（604）→ `w = q_half·pref / area_weight`（609）。
- `d2dga_flux.py` **无 I1/I2 函数**，仅有 `d2dga_flux_amplification`/`d2dga_dispersion_function_I3`/`d2dga_buoyancy_flux`。

**设计**：
1. **新建 I1/I2 函数**于 `d2dga_flux.py`（Zhang 2022 式 4.13/4.21）：
   - `d2dga_dispersion_I1(c_bar, m)`：对流放大分层修正项。
   - `d2dga_dispersion_I2(c_bar, m)`：方位再分配分层修正项。
   - 端点性质：ḉ=0/1 处 I1/I2 取已知极限值（实现时对照论文式号核实）。
2. **体力向量注入流动度**（式 2.5b）：在 `_compute_velocity` 的 `pref` 乘子处，用真体力向量 `f=(f_phi, f_xi)` 替换 `(2φ−1)` 简化代理：
   - `pref = base · (1 + Δρ · f · I2(c̄, m))` 形式注入方位角流动度重分配（重顶替轻→窄边流动度提升）。
   - `f` 由 `_buoyancy_force_vector(geom, beta_deg)` 给出（已存在，复用）。
3. **I1 分层修正**（式 4.24）：对流项（`base` 中的放大因子）乘 `I1(c̄, m)` 体现水泥/泥浆分层。
4. **保留 `enable_true_buoyancy` 开关**：True=真体力注入+I1/I2（式 2.5b）；False=回退 `(2φ−1)` 简化代理（R0-R2 消融可复现）。
5. **数值安全**：体力项 `Δρ·f·I2` 幅度局部 CFL 裁剪（同 T1-2，非全局限幅）。

**测试**（`tests/test_d2dga_flux.py` + `tests/test_improved_d2dga_annulus.py`）：
- I1/I2 端点性质（ḉ=0/1）。
- 合成均匀密度差+已知井斜工况，验证体力向量方位角重分配方向（重顶替轻→窄边 `pref` 提升）。
- `enable_true_buoyancy=False` 时回退 `(2φ−1)` 简化（R0-R2 可复现）。

### 2.4 T1-4 两层黏度闭包（`annulus_d2dga.py:424-464,566-604`，式 4.23）

**现状**：
- `_compute_props`（424-464）按相体积分数加权混合黏度，内部已算 `μ_mud_field`（453）、`μ_cement_field`（455-457）但仅用于构造 `m_field`，**未作为相黏度场返回**。
- `_compute_velocity` 正则化为 `μ_reg`（549-559），流动度 `base=(b/mean(b))²/μ_reg`（570），单相混合黏度模型。

**设计**：
1. `_compute_props` 额外返回 `η1=μ_mud_field`、`η2=μ_cement_field`（相黏度场）。返回签名 5→7 元组（`mu, rho, mud, tau_y, m_field, eta1, eta2`）。
2. **注意测试同步**：`tests/test_improved_d2dga_annulus.py::TestMFieldFromProps::test_m_field_returned_and_shape` 期望 `len(out)==5`，需同步更新为 7。
3. `_compute_velocity` 用两层闭包替换 `μ_reg`：`1/η_mix = c̄³/η₂ + (1−c̄³)/η₁`，`base=(b/mean(b))²/η_mix`。
4. 数值安全：`η_mix` 下限 epsilon（如 1e-9）防除零。

**测试**：
- `c̄=0 → η_mix=η1`（纯泥浆）、`c̄=1 → η_mix=η2`（纯水泥）端点。
- 混合值量级介于 η1/η2 之间。
- `_compute_props` 返回 7 元组的形状断言。

### 2.5 T1-5 static wall layer c_min 判据（`annulus_d2dga.py:668,627,923`，式 2.35-2.41）

**现状**：
- `wall` 场在 `run` 中初始化为零（668），注释说明论文口径下不做壁面泥饼清洗，恒为零场。
- 无 `c_min` 判据，`wall_field` 仅作结果字段输出（923）+ 快照（814/932），`_depth_profiles` 中 `del wall`（627）丢弃，**不参与计算**。

**设计**：
1. 实现 `c_min` 判据（Bararpour 2025 式 2.35-2.41）：局部水泥浓度 `c < c_min` 处壁面层泥浆滞留不流动，`wall=1`；否则 `wall=0`。
2. `wall` 反馈流动度：`wall=1` 处 `pref *= (1−wall)`（壁面静止层，速度归零）。
3. 去掉 `del wall`（627），`wall_field` 参与 `_compute_velocity`。
4. `c_min` 作求解器参数（默认 0.05，可配）。

**测试**（`tests/test_improved_d2dga_annulus.py`）：
- `c > c_min` 全场 → `wall` 全零。
- `c < c_min` 区域 → `wall=1` 且该处速度=0。
- `c_min` 参数可配。

---

## 3. 集成与数据流

- **`AnnulusSimulationResult` 结构不改**：诊断/指标仍走 `summary` dict 扩展（沿用 Tier 0 模式）。
- **`_compute_props` 返回签名变更**（5→7 元组）是本 spec 唯一破坏性接口改动，影响所有调用点（`run` 主循环、`_compute_velocity`）。
- **三开关语义保留**：`enable_d2dga_auto_m`/`enable_d2dga_i3_flux`/`enable_true_buoyancy` 做 R0-R3 消融对照；去限幅是无条件的（开关只控 R0-R3 的物理项开关，不限幅）。
- **入口/出口/快照/时间步不动**：`save_interval`、`dt`、`nz` 沿用现场模式。
- **Tier 0 诊断复跑**：`compute_all_tier0_diagnostics` 不改，但实跑结果会变（去限幅后物理量变）。

---

## 4. 测试策略

| 层 | 内容 | 文件 |
|----|------|------|
| 通量纯函数 | I1/I2 端点、f(c,m) 去 clip 边界、`d2dga_buoyancy_flux` 量级 | `tests/test_d2dga_flux.py` |
| 求解器单元 | `_compute_props` 7 元组、两层闭包端点、`wall` c_min 判据、体力向量方位重分配 | `tests/test_improved_d2dga_annulus.py` |
| 消融回归 | R0-R3 四级（不崩溃+物理合理，**不对标 CBL**） | `tests/test_improved_d2dga_annulus.py` + `runners/ht1_004_ablation.py` |
| 六井集成 | 不崩溃+物理合理 | `tests/test_six_well_integration.py` |
| 实跑+诊断 | hu102/ht1_004 实跑 + Tier 0 诊断复跑（Muskat regime/浮力分类应改善） | `scripts/run_tier0_diagnostics.py` |

**R0-R3 回归基线处理**：测试无硬编码效率期望值（现状），去限幅后只需保证"不崩溃 + 物理合理"（m>1、η_N≤η_E、t_br>0、flow_class 非空）。R3 效率值会变，论文图6-9 需重跑（已接受）。

---

## 5. 验证门（无 CBL）

| 验证类型 | 内容 | 用于 |
|---------|------|------|
| 理论对标 | Zhang 2023 Table 3（t_br/η_E）、Bararpour 2025 regime 边界、Yang 2021 序列效率（89.67% vs 85.88%）| T1-1~T1-5 |
| 单元测试 | 每项对解析值/已知极限 | 全部 |
| 消融回归 | R0-R3 四级 + 六井集成，仅检查不崩溃+物理合理 | T1 全部 |
| Tier 0 诊断复跑 | Muskat regime/浮力分类/η_N 应随物理化改善 | T1-2/T1-3 |

---

## 6. 执行计划（子 agent 分工）

### 6.1 执行结构

| 步骤 | 子 agent | 任务 | 产出 |
|------|---------|------|------|
| **SA1** | 通量物理化 | `d2dga_flux.py`：T1-1 去 clip + T1-3 新建 I1/I2 函数 + 单元测试 | 模块 + `tests/test_d2dga_flux.py` 增量 |
| **SA2** | 求解器物理化 | `annulus_d2dga.py`：T1-2 去 0.05 + T1-3 体力注入 + T1-4 两层闭包 + T1-5 static wall + 单元测试 | 模块 + `tests/test_improved_d2dga_annulus.py` 增量 |
| **主 agent** | 统筹 | 审查 SA1/SA2 产出、跑全量测试、R0-R3 消融回归 + 六井集成 + 实跑诊断、解决冲突 | — |

### 6.2 执行顺序

1. **SA1 先行**：I1/I2 函数 + 去 clip（T1-1 + T1-3 的 d2dga_flux 部分）。SA2 的体力注入依赖 I1/I2。
2. **SA2 接力**：T1-2 去 0.05 → T1-4 相黏度返回+两层闭包 → T1-3 体力注入（用 SA1 的 I1/I2）→ T1-5 static wall。
3. 主 agent 跑 `pytest tests/` 全量 + R0-R3 消融（`run_one_level` 四级）+ 六井集成 + `run_tier0_diagnostics.py` 实跑（hu102/ht1_004，不对标 CBL）。

### 6.3 子 agent 统一约束（写进每个 SA prompt）

- 遵循现有 `d2dga_flux.py`/`annulus_d2dga.py` 代码风格：类型注解 + 中文 docstring + 公式标注论文式号。
- 数值安全：除零保护、`np.clip` 仅用于防除零（非物理限幅）、NaN 守卫。
- 每项配单元测试，用合成数据测已知极限。
- 不跨层：只动 `d2dga_flux.py`/`annulus_d2dga.py` + 对应测试。
- 环境：`conda activate shenjingwangluo`，`pip install -e .`，测试 `python -m pytest tests/test_<module>.py -v`。

---

## 7. 风险与缓解

| 风险 | 缓解 |
|------|------|
| 去限幅后数值发散（大排量/大密度差）| 公式内 `c_safe` clip + 局部 CFL 裁剪（`|div_q|·dt ≤ α·Δs`）；R0-R3 消融逐级验证 |
| `_compute_props` 返回签名 5→7 破坏测试 | 同步更新 `test_m_field_returned_and_shape` 及所有调用点 |
| 体力向量注入改变速度场结构（T1-3 最高风险）| `enable_true_buoyancy` 开关保留 False 回退；R0-R3 消融逐步验证稳定性 |
| I1/I2 公式需查论文核实 | 实现时对照 Zhang 2022 式 4.13/4.21/4.24，端点性质作为测试锚点 |
| R3 效率值变化影响论文图6-9 | 已接受，重跑图6-9（`run_ht1_004_ablation.py`） |
| static wall c_min 阈值经验性 | 默认 0.05 可配，敏感性留待 Tier 0 诊断复跑后调 |

---

## 8. 后续 Tier 衔接（本 spec 排除项）

本 spec 排除的 T1-6~T1-9 各自后续出子 spec：
- **T1-6 FLUSHER 流体角色 + 独立浓度场**：跨 L1/L3/L4，数据结构改动最大，**优先单独 spec**。
- **T1-7 CFL 自适应时间步**：L4 数值稳定性，去限幅后可能需提前（作为 T1-1~T1-5 的数值保护后备）。
- **T1-8 1D HB/Bingham 屈服 Taylor 弥散**：L2 `casing_flow.py` 独立。
- **T1-9 入口相分数方位加权**：L3 `boundary_bridge.py` + L4。

Tier 1 完成后进 Tier 2（HB D2DGA 闭包查表 / augmented Lagrangian / TVD 激波捕捉 / 瞬态流函数 / 泵注优化器），各出详细 spec。

---

## 9. 验收标准（Tier 1 完成定义）

- [ ] T1-1~T1-5 五项全部实现，`pytest tests/` 全量通过（含更新后的 `_compute_props` 7 元组测试）。
- [ ] R0-R3 四级消融回归通过（不崩溃 + 物理合理，不对标 CBL）。
- [ ] 六井集成（`test_six_well_integration.py`）通过。
- [ ] hu102 + ht1_004 实跑，Tier 0 诊断复跑：Muskat regime/浮力分类/η_N 随物理化改善（b>0、flow_class 非空、η_N≤η_E、t_br>0）。
- [ ] `d2dga_flux.py` 新增 I1/I2 函数 + 端点测试。
- [ ] `_compute_props` 返回相黏度场 η1/η2，两层闭包替换 μ_reg。
- [ ] `wall` 场按 c_min 判据激活，参与 `_compute_velocity`。
- [ ] 论文图6-9 重跑（R3 数据更新）。






