# 内部自洽加固 · 设计规格

> **日期**：2026-09-25
> **代码锚点**：`feat/hb-closure-phase-a` @ `bc155b0` + 工作树未提交改动
> **触发**：用户裁定验收标准 = **内部自洽优先**（"至少逻辑是通的，能采信这个模型"）
> **上位调研**：obsidian《模型模块现状与环空段模拟调研_2026-09-25》《三项需求方案调研_敏感性_套管段二维_深度浓度追踪_2026-09-25》
> **性质**：本规格只定义"让代码实际行为 = 其声明行为，并建立可复核的自洽校验闸门"。**不改变任何物理主张、不引入新物理**。

---

## 1. 目标与非目标

### 1.1 目标
1. **消除静默失效**：任何"置真却在当前路径上无消费者"的开关都必须告警（当前有 5 类空转开关**零告警**）。
2. **消除输入口径的静默偏差**：已存在的管容锚点常量必须真正进入 WellSpec；停泵沉降半径必须用本井内径而非硬编码 0.05 m。
3. **消除调用面腐烂**：形参已删却仍被传的调用点（当前 4 处，必抛 `TypeError`）必须修复，并加**自动闸门**防复发。
4. **建立可复核的自洽证据**：质量守恒、`η_E ≡ 1 − 饥饿体积份额`、前缘到达时刻 vs 解析解、死开关不变性——产出可归档的台账。
5. **记录位移**：上述修正前后，逐井记录 η 变化（**接受数字变化，包括变差**）。

### 1.2 非目标（明确排除，另案裁定）
- 接 (4.25) 通量形式（把 f_amp 从速度层移到通量层）——物理变更，属"物理保真优先"路线。
- HB 路径提速、套管段二维化、指标换主验证量（η_N / 1/η_N）。
- 调整任何**有量纲标定钮**（`dispersion_alpha=0.25`、`yield_gate_f_safety=1.15`、`K_AXIAL=1/3`、`re_crit`、`μ` clip 上限 3.0 Pa·s、`w_prev` 初值 0.45 m/s）。
- 重跑论文数字、更新论文正文。

---

## 2. 设计决定

| # | 决定 | 理由 |
|---|---|---|
| **D1** | **死开关守卫白名单化 + 双向完备**：覆盖 `enable_regime_split` / `enable_true_buoyancy` / `enable_power_law_gap_law` 在 `enable_stream_function=True` 时的空转，以及 `enable_yield_gate=True` 但 `enable_stream_yield_gate=False` 时"wall 算而不用"。**一次性 `UserWarning`，不改任何数值** | 现状守卫只检查反方向组合（`annulus_d2dga.py:561-576`），默认路径的空转零告警 |
| **D2** | **q₀ 不改接线，改为声明一致化**：`two_layer.isotropic_flux_q0` 与 `docs/源模型口径与适用域声明.md` §1 声明 5 明确标注"默认路径无消费方，仅作 `gap_solver` 回归锚" | 接线属物理变更（非目标）；但"声明说弥散由 q₀+I₃ 承载"与代码不符，必须改文字 |
| **D3** | **管容锚点接入**：`ht1_001_loader` / `hu2_loader` 把既有常量传入 `WellSpec.shoe_lag_volume_m3`，与 `ht1_003_loader.py:395` 同构 | 常量已存在且注释完整，只是从未接线；改的是 1D 输入口径 |
| **D4** | **停泵半径传参**：`_effective_pipe_radius_m(well_spec)` 由调用方传入本井 `liner_id_mm`；并**证明**该路径在生产时间窗内不可达（预期零生产影响） | `casing_flow.py:1371` 无参调用 ⇒ 恒 0.05 m，真实 0.0539–0.0556 |
| **D5** | **调用签名闸门**：新增 AST 扫描器 + 契约测试，冻结 `AnnulusD2DGASolver` / `CasingFlowSolver` 的全部生产调用点 | 仓内无此类工具；形参删除后腐烂 4 处且无人察觉 |
| **D6** | **官方口径脚本与生产 runner 的 4 处口径差：本次只出对照报告，不改** | 属口径决定，须单独裁定（见 §4 R-A） |

---

## 3. 验收判据

1. 新增契约测试全绿；`tests/contract/` 除既存 1 红（`test_hb_newtonian_limit.py::test_frozen_snapshot_stream_function_defaults`）外全绿。
2. **不变性可执行证明**：对每个"死开关"，`on/off` 两次运行结果**逐位相同**（由既有逐位锚测试 `test_six_well_integration.py` 承担）。
3. 一致性台账产出全部可跑井的：守恒误差、`(1−η_E) − 饥饿份额`、前缘到达时刻（2D）vs 解析（1D）。
4. `results/内部自洽加固_2026-09-25/位移台账.csv` 逐井含 `井名, 修正项, 修正前eta_E, 修正后eta_E, Δeta_E_pp, 修正前eta_N, 修正后eta_N, Δeta_N_pp, 说明`。
5. 每条结论附 `文件:行号` 或产物路径，可被第三人复核。

---

## 4. 需用户裁定（不影响前 5 个任务开工）

| # | 分歧点 | 建议默认 |
|---|---|---|
| **R-A** | 官方口径脚本（`rerun_all_wells_corrected.py`）是否收口到生产 runner 口径（1D 三开关 + `schedule=` + 统一 `total_t`） | **收口**（生产口径为真相），但会再次移动 G3 八井数字 |
| **R-B** | P1（`enable_stream_yield_gate=True`）是否进生产口径 | 先**不**进，作为 A/B 记录（它只对 1/8 井有量级，但强化居中度杠杆） |
| **R-C** | 位移台账的接受阈值（Δη_E 超过多少必须停下复查） | 建议 **>2 pp 逐案复查**，≤0.5 pp 记录即可 |
| **R-D** | `save_interval` 是否由"步数"改"物理时间"（影响浓度-时间曲线口径） | 另案，不在本计划内 |

---

## 5. 风险

| 风险 | 缓解 |
|---|---|
| 管容修正使 ht1_001/hu2 前缘时刻后移 ⇒ 其 η 位移可能不小 | 单井 A/B + 位移台账；先跑这两口井 |
| 守卫白名单写错 ⇒ 误报噪音，掩盖真告警 | 每条守卫规则都配"不告警"的反例测试 |
| 一致性闸门暴露更多问题（例如守恒率超阈） | 这是**目的**：如实记录，不修 |
| 用户后续选择物理路线 ⇒ 本计划的基线再次作废 | 本计划只锚定"当前口径 + 修正项"，与物理路线解耦 |

---

## 6. 交付物

| 产物 | 路径 |
|---|---|
| 基线冻结 | `docs/superpowers/plans/baseline-2026-09-25.md` |
| 死开关守卫 + 测试 | `cemdisp/models2d/annulus_d2dga.py`、`tests/contract/test_dead_switch_guard.py` |
| 管容接入 + 测试 | `cemdisp/data/loaders/{ht1_001,hu2}_loader.py`、`tests/contract/test_shoe_lag_wiring.py` |
| 半径传参 + 测试 | `cemdisp/transport1d/casing_flow.py`、`tests/contract/test_settled_pipe_radius.py` |
| 签名闸门 + 测试 | `scripts/entrypoints/check_call_signatures.py`、`tests/contract/test_entrypoint_signatures.py` |
| 一致性闸门 | `cemdisp/diagnostics/internal_consistency.py`、`scripts/entrypoints/verify_internal_consistency.py`、`tests/contract/test_internal_consistency.py` |
| 位移台账 | `results/内部自洽加固_2026-09-25/位移台账.csv` 等 |
