# Phase 1 设计规格：「fluid_at 签名波」= R1 参数化 + P-1 静压（2026-10-06）

> **上游**：《详细执行计划_真温压响应_2026-10-06》§3 Phase 1 + §4 五关 + §6 先报后动 + §7 红线；
> 对抗评审《评审A_R1-R4设计对抗评审_2026-10-06》§1 缺陷总表 + §4 补充断言清单 A1.x / A3.x；
> 实施面映射《实施面映射C_四大空缺与敏感性资产_2026-10-05》§2（空缺二 Phase P）。
> **前置**：Phase 0.0（Q16 密度就近取）与 Phase 0（三重点井补跑）已验收。
> **基线**：HEAD `c039915`；全量 `pytest tests/ -q` = **17 failed / 1006 passed / 241 subtests**
> （2026-10-06 实测 664.90s，失败清单与基线逐条一致）。
>
> **评审重排依据**：R1 与 P-1 同动 `fluid_at` 签名 / 调用点 / memo 键，合并一波避免二次返工
> （计划 §3 Phase 1 导语）；评审 D-05 裁定 **P-1 必须先于 R2**（隔离液 τy 经混合 τy 场
> `annulus:1321-1327` 直接喂 `_yield_gate_wall`，P-1 使 τy +17~29% ⇒ wall 通道增强）。

---

## 0. 改动面总览（两 commit）

| commit | 内容 | 文件 |
|---|---|---|
| **C1** `feat(temp): R1 温变流变公式系数参数化` | `RheologyFormulaParams` + 全链下传 + solver kwarg + run_opts 键 | `cemdisp/data/rheology_vs_temperature.py`、`cemdisp/data/__init__.py`、`cemdisp/models2d/annulus_d2dga.py`、`cemdisp/transport1d/casing_flow.py`、`scripts/entrypoints/run_sensitivity_current_20260916.py`、`scripts/entrypoints/run_sensitivity_temperature_t2_20261001.py`、`tests/contract/test_rheology_formula_params.py`（新） |
| **C2** `feat(temp): P-1 静压场接入` | `PressureField` + `HydrostaticPressureField` + 三井对账 | `cemdisp/data/pressure_field.py`（新）、`cemdisp/data/__init__.py`、两个 solver、两个敏感性脚本、`scripts/entrypoints/verify_temperature_coupling.py`、`tests/contract/test_pressure_field.py`（新）、`tests/contract/test_pressure_wiring.py`（新）、ADR-0001、CONTEXT.md |

**红线自检（逐条对 §7）**：`_SWITCH_DEFAULTS:192-220` **一字不动**；一切走运行侧 opt-in；
关1 锚 `results/_baseline_T_off/基线摘要_T_off.json` 实值（η_E 0.9996663176608406 /
η_N 0.9983246031489231，**读 JSON 比对，禁硬编码**）；预存 17F 勿动；源目录冻结。

---

## 1. R1：温变流变公式系数参数化

### 1.1 `RheologyFormulaParams`（frozen dataclass，置于 `rheology_vs_temperature.py` 内）

字段分组、默认值 = 现行字面量，逐字搬运（**不改任何数值、不改任何表达式次序**）：

| 组 | 字段（含默认） |
|---|---|
| 断点 | `break_T=100.0` |
| 组A 水泥 | `ca_tau0_q=(0.002793, -0.391702, 21.6378)`；`ca_tau0_l=(0.624883, -51.7949)`；`ca_mup=(0.53129, -0.012495, -0.00204)`；`ca_T_lo=20.0`；`ca_T_hi=170.0` |
| 组B 水泥 | `cb_tau0_q=(4.433e-4, -0.07843, 5.962)`；`cb_tau0_l=(0.03206, 1.979)`；`cb_mup=(1.355e-5, -6.570e-4, 0.08215)`；`cb_T_lo=20.0`；`cb_T_hi=200.0` |
| 隔离液 | `sp_1p95_ty`、`sp_1p95_mp`、`sp_2p05_ty`、`sp_2p05_mp`（各 6 元组，逐字）；`sp_T_lo=20.0`、`sp_T_hi=200.0`、`sp_P_lo=0.1`、`sp_P_hi=200.0`、`sp_p_default=0.1`；`sp_lo=1.90`、`sp_hi=2.10`、`sp_mid=2.000`、`sp_anchor_lo=1.95`、`sp_anchor_hi=2.05` |
| 水泥密度域 | `cm_lo=1.88`、`cm_hi=2.12`、`cm_mid=2.000`、`cm_anchor_a=2.10`、`cm_anchor_b=1.90` |
| 钻井液 | `mud_mup=(0.297154, -0.02310525)`；`mud_tauy=(21.346710, -0.00830631)`；`mud_T_lo=40.0`、`mud_T_hi=80.0` |

**逐位约束（关2 硬红线）**：表达式**保持现有运算次序**。
`a*T*T - b*T + c` 只能落成 `q0*T*T + q1*T + q2`（`q1 = -b`，IEEE754 下 `x-y ≡ x+(-y)`）；
**禁止 Horner 化**——实测（210 温点网格）：加法形式 **0 处变位**，Horner 形式 **109/210 处变位**。

### 1.2 全链下传

```
fluid_at(fluid, T_c, P_mpa=None, *, mud_extrapolate=False, smooth_break=False,
         params: RheologyFormulaParams | None = None)
   ├─ _mud_values(T, extrapolate, name, params)
   ├─ _spacer_values(fluid, T, P_mpa, params)
   ├─ _cement_values(fluid, T, params)
   └─ _route(fluid)                      ← **不动**（路由是结构非系数）
```

`params is None` ⇒ 模块级 `_DEFAULT_PARAMS` 单例。`__all__` += `RheologyFormulaParams`。
`_cement_a_tau0 / _cement_a_mup / _cement_b_tau0 / _cement_b_mup / _mud_mup / _mud_tauy`
私有签名加 `params`。
**仓内无任何测试/脚本直接 import 这些私有函数**（已核实：外部仅引用 `_route/fluid_at/get_audit/reset_audit`）。

### 1.3 solver 侧

| 文件 | 改动 |
|---|---|
| `annulus_d2dga.py` ctor 尾（`:453-454` 之后） | `rheology_formula_params`、`mud_extrapolate` 两个新形参（默认 `None` / `False`） |
| `annulus_d2dga.py` `__init__` 体（`:699-711` 邻近） | `self.rheology_formula_params = ...`、`self.mud_extrapolate = ...`（**不建场、零痕迹**） |
| `annulus_d2dga.py` `_phase_props` `:1193` | `fluid_at(fluid, T, params=..., mud_extrapolate=...)` |
| `annulus_d2dga.py` `run()` 构造层 `:2395/:2398/:2401/:2404` | 同参（**D-06：构造层与逐步同口径，不得只改一处**） |
| `casing_flow.py` ctor（`:169` 之后） | 同两个 kwarg |
| `casing_flow.py` `_phase_props` `:533` | 同参 |

### 1.4 敏感性链 run_opts（**白名单两处同改 + 映射两处**）

真闸 = `normalize_run_opts` 默认 dict（`:203-207`），校验比较对象是它（`:209-211`），**不是**
`RUN_OPTS_KEYS`（`:195`）。但评审 **A1.2** 要求新增单测把两处钉在一起 ⇒ **两处都改**：

| 位置 | 改动 |
|---|---|
| `RUN_OPTS_KEYS` `:195` | += `"rheology_formula"`, `"mud_extrapolate"` |
| `normalize_run_opts` 默认 dict `:203-207` | += `rheology_formula: None`, `mud_extrapolate: False` + 类型校验（dict｜None / bool） |
| `casing_kwargs_from_opts` `:237-239` | 映射到 `rheology_formula_params` / `mud_extrapolate` |
| `annulus_kwargs_from_opts` `:242-255` | 同上 |
| T2 脚本 `_opts()` `:181-187` | 加两个同名形参（默认 `None`/`False` ⇒ 现有调用逐位不变） |

`rheology_formula` 载荷（dict）：`{"scale": {...}}` 形式——按字段名给出乘子，经
`_params_from_spec()` 转 `RheologyFormulaParams`（`dataclasses.replace`）。**默认 `None` ⇒ 默认参。**

### 1.5 验收（= **L2 流变响应达成**）

| 关 | 判据 |
|---|---|
| 关1 | T-off 路径**零改动**（`fluid_at` 根本不被调用）⇒ 关1 锚逐位 |
| 关2 | `fluid_at(f,T,P)` ≡ `fluid_at(f,T,P,params=None)` 逐位；关2 审计事件**按 kind 分桶计数**逐键对齐（评审 D-12） |
| 关4-A1.1 | 全域网格逐位：全族 × T∈[20,170]×5℃ × P∈{None,0.1,70,140} 输出 `(τy,μp)` 逐位等；**并断言全域 τy≥0 ∧ μp>0**（评审 D-09 负物性敞口；扰动档同样闸） |
| 关4-A1.2 | `RUN_OPTS_KEYS == tuple(normalize_run_opts({}))`；每键经 kwargs 映射**到达 solver 属性**（构造后 `getattr` 断言） |
| 关4-A1.3 | 同算例重构前后 `temperature_rheology_audit` 按 kind 计数全等 |
| 关4-锚点表 | 隔离液 2.05 式 @20°C,0.1MPa：τy=`10.27960477731`、μp=`0.10536316320405001`（实测值，**禁手抄**）；泥浆 τy(80°C)=`10.98348553218639`；扰动组A τ₀ 系数 ×1.1 → `fluid_at` 输出 = 手算 |
| 排量响应锚 | 呼101 static ×0.6/1.0/1.4 η_N = 0.8264/0.7777/0.7302 单调（本阶段应**不动**——R1 默认参不改数） |

---

## 2. P-1：静压场接入

### 2.1 新建 `cemdisp/data/pressure_field.py`（同型于 `temperature_field.py`）

- `PressureField` Protocol：`P(md_m: float, t_s: float) -> float`（**MPa**）
- `ConstantPressureField(p_mpa)`：常数场（关2 自洽 / 对照档）
- `HydrostaticPressureField`：
  - **TVD 口径**：`TVD(md) = ∫₀^md cosθ(md') dmd'`，θ 取自 **`well.inclination_profile`**
    （`DepthValuePoint(depth_md_m, value=inclination_deg)`，与 loader 读的 `inclination_profile.csv` **同源**）；
    剖面首点以上按 0°（直井）外推
  - 密度剖面：均匀单值 **或** 分层 `((top_md, rho_kg_m3), ...)`
  - `P(md, t) = g · ∫₀^{TVD(md)} ρ dTVD' / 1e6`，`g = 9.80665`（与停用绘图段
    `plot_jieti_timing_figs_20260928.py:353-466` 同源；MATLAB 参照 `p_jaifang1.m:701-715` 用 9.81，
    差 0.03%，口径注记）
  - `from_well(well, rho_kg_m3)` / `from_layers(well, layers)`；`oob_count`/`reset_audit` 同型
- `insitu_column_density(fluids, schedule)`：**在场相密度按设计泵注体积加权均值**（Q9/D5 口径 ③；
  素材C「相密度用当前在场相常数」）；排除 `volume==0` 或不在 fluids 表中的步骤

### 2.2 调用点补传 P（**三处，且必须同口径**）

| 位置 | 口径 |
|---|---|
| `annulus._phase_props` `:1171-1195` | `P = self._representative_pressure_mpa(geom, t)`；`pressure_caliber="shoe"` ⇒ `P(max(geom["md"]), t)`；`"mean"` ⇒ `mean(P(md_i, t))` |
| `annulus.run()` 构造层 `:2395-2404` | **同传**（D-06 裁定：诊断物性 ≡ 逐步物性；`_temp_rheo_fluids` 与求解同源） |
| `casing._phase_props` `:514-535` | 鞋深单点同口径 |

**memo 键 (fluid, T) → (fluid, T, P)**（annulus `:1189`、casing `:529-532`）——否则缓存污染（红线）。
**HB memo key 温度维 `:1949-1953` 不改**（HB 默认关；R4 期另裁 —— 评审 D-04）。

### 2.3 solver kwargs

两个 solver 各加 `pressure_field`（默认 `None`）、`pressure_caliber`（默认 `"shoe"`）。
`pressure_field is None` ⇒ `P=None` ⇒ `fluid_at` 走 `p_default` 审计 ⇒ **逐位 = 旧行为**（关2）。

### 2.4 敏感性链

run_opts 再加 `pressure_mode ∈ {"off","hydrostatic"}`（默认 `"off"`）、
`pressure_caliber ∈ {"shoe","mean"}`（默认 `"shoe"`）；09-16 脚本加
`build_pressure_field(well, fluids, schedule, mode)`（仿 `build_temperature_fields`）。

### 2.5 验收

| 关 | 判据 |
|---|---|
| 关2 | 无 `pressure_field` ⇒ 审计 `p_default` 语义逐位保留；关1 锚逐位 |
| 关4-A3.1 | T-on 首步 `_temp_rheo_fluids[phase] is _phase_props(phase, geom, t0)`（同 memo 路径 ⇒ 同一对象） |
| 关4-A3.2 | 三井 vs ρgh 参照**逐节点 diff=0**；**呼1-003 P(shoe) 对设计锚 145.96 MPa 在 ±2% 内**；T-on+传 P 的 run 内 `p_default` 审计计数 **必须为 0** |
| 关4-A3.3 | 同 run 内 (fluid,T,P1)/(fluid,T,P2) 两次查询返回**不同对象**且各自 `_quad6` 手算对（键未含 P 则必红） |
| 关4-口径差表 v2 | `verify_temperature_coupling.py` 加 P 口径 → **新目录**（v1 作废，红线） |
| 关3 | 方向：隔离液 τy 随 P 增（1.95 式 +17%／2.05 式 +29% @70MPa,80°C）；μp 符号随温翻转 |

**独立复算预置（本 spec 编制时实测，供验收比对）**：

| 井 | TVD(shoe) | ρ̄（体积加权） | P(shoe) g=9.80665 | 现场／设计锚 | 偏差 |
|---|---|---|---|---|---|
| 呼101 | 7863.16 m | 1975.58 | **152.34 MPa** | 152.16 MPa（**实测**，notes 文本列 7868 m） | **+0.12%** |
| 呼1-003 | 7617.24 m | 1949.53 | **145.63 MPa** | 145.96 MPa（**设计值**，notes 文本列 7618 m） | **−0.23%** |
| 呼1-004 | 7656.80 m | 1881.05 | 141.24 MPa | 无鞋深压力锚；地层压力当量 1.933 → 145.2 MPa，**低于**该值 ⇒ 与 P-4 的 MPD 回压／ECD 属同一缺口，如实记录 | — |

⚠ 计划 §3 原文「ρ≈1.94 时 ρgh≈145.6」**算术有误**（源出素材C 误用呼1-004 井深 7660 m）：
呼1-003 井深 7618 m 下 ρ=1.940 → **144.98 MPa**；得 145.96 需 ρ≈**1.9533**。本 spec 以上表为准。

---

## 3. 偏离登记（§6-1，执行窗口已记录）

| # | 计划原文 | 实际 | 性质 |
|---|---|---|---|
| Δ1 | 「run_opts 白名单**两处**同改（元组 :195 + normalize 默认 dict）」 | `RUN_OPTS_KEYS` 全仓**零引用**，改它无效；真闸只有 normalize 默认 dict。但评审 A1.2 要求两者钉一起 ⇒ **两处都改 + 新增 A1.2 单测** | 计划半实，已按评审补强 |
| Δ2 | 改动面未列 annulus `run()` 构造层 4 个 `fluid_at` 调用点 | 已补列（D-06 裁定必须同传） | 计划漏列 |
| Δ3 | 「呼1-003 环空静压 145.96 MPa **实测**锚」 | 该数仅在 `temperature_pressure_profile.csv` 的 **notes 文本列**（`pressure_mpa` 列全空），自标「(设计)」⇒ **设计锚** | 计划措辞失真 |
| Δ4 | 「MD→TVD 用 `inclination_profile.csv`」 | 改用 `well.inclination_profile`（**同一 CSV** 经 loader 解析；呼1-003 该 CSV 的 `tvd_m` 列全空、井斜 ≤4.12° 须自行积分） | 实现路径简化，同源 |
| Δ5 | 计划未给 P 的 run_opts 键 | 新增 `pressure_mode`／`pressure_caliber`（无此键 P-1 在敏感性链中为**死代码**，Phase 2「P-1 后口径」批无法执行） | **计划缺口**，须用户追认 |
| Δ6 | Phase 1.5 的 `enable_stream_yield_gate` 称「第 5 键」 | 加 P 两键后成为第 8 键 | 键序漂移，无语义影响 |
| Δ7 | 计划 §8/§10 台账路径 `docs/superpowers/progress.md` | 不存在；真实台账 = `.superpowers/sdd/2026-10-06-thermobaric-true-response/progress.md` | 计划路径错（已在台账登记） |
| Δ8 | 计划未提 TVD 数据源二选一 | hu101/ht1_004 的 `inclination_profile.csv` **有**权威 `tvd_m` 列（比自积分准 3.15 m ≈ 0.04%）；呼1-003 该列**全空**，必须自积分。本阶段**统一走 `well.inclination_profile` 自积分**（单一路径、三井一致、误差 ≤0.05%），CSV 直读列为 P-2 细化项 | 实现选择，误差已量化 |
| Δ9 | 计划未裁「在场相密度」的角色口径 | 默认 = **全泵注体积加权**（含管内替浆/压塞液），实测两口径对呼1-003 差 0.02%（1.9495 vs 1.9492）；提供 `roles=` 供细化。**须用户追认** | 口径裁量，登记 |

---

## 4. 回退

两 commit 均为**加法式**（新参默认 = 旧行为）。`git revert` 单 commit 即回到 HEAD 口径；
新产物目录可整目录删除，不影响任何既有产物。

---

## 6. Phase 1.5 追加：激活管道骨架（2026-10-06 同窗口执行）

计划 §3 Phase 1.5「独立小 PR」。用户 2026-10-06 裁定 Δ5/Δ9/Δ10/C-14 全部按推荐口径追认，
本波继续。

### 改动面

| 位置 | 改动 |
|---|---|
| `RUN_OPTS_KEYS` | += `"enable_stream_yield_gate"`（第 8 键） |
| `normalize_run_opts` 默认 dict | += `"enable_stream_yield_gate": None`（**None ⇒ 不给键**，沿用 annulus 构造默认 False）；类型闸 bool｜None |
| `annulus_kwargs_from_opts` | `is not None` 时出键 |
| `casing_kwargs_from_opts` | **禁出此键**（casing 无该形参；出键=运行期 TypeError，见 C-01 不对称版） |
| T2 `_opts()` | += `stream_yield_gate: bool \| None = None` |

### 语义区分（计划 §3 Phase 1.5 点名要求写入 docstring）

| 开关 | 端 | 含义 | 默认 | T2 既有对照用的是 |
|---|---|---|---|---|
| `enable_yield_gate` | **生产端** | wall 场**算不算** | True | ✅ 本键（gateoff 对照） |
| `enable_stream_yield_gate` | **消费端** | wall **进不进**流函数算子 | False | ✗ |

二者组合的死活判定见 `annulus._dead_switches`（置真但 wall 到不了算子 ⇒ 判"死开关"告警）。

### 执行窗口发现的语义要点（**重要**）

`enable_stream_yield_gate` **不是 T-on 专属开关**——它把 wall 接进算子，**在 T-off 下同样改变
流场**（实测：合成小算例 T-off 下置真 ⇒ `effective_efficiency` 变化）。故契约测试的判据是
两段式：

- **不给键 / 给键=默认值** ⇒ 与基线**逐位相同**（计划 §3 Phase 1.5 的"无此键时行为逐位不变"）；
- **给键=True** ⇒ 结果**必须变** —— 这正是"确实到达算子"的证据（反制评审 A1.2 的
  "添加了 ≠ 在运行"）。

### 验收

| 项 | 结果 |
|---|---|
| `test_signature_defaults_match_table`（计划点名） | **绿**（本波不动 `_SWITCH_DEFAULTS`） |
| 全 8 键穿透（A1.2） | `test_every_key_reaches_solver_attribute` —— 构造后逐属性断言 |
| casing 不对称 | `test_casing_kwargs_never_emit_stream_yield_gate`（`inspect` 证明 casing 无该形参） |
| 关2 无键逐位 | `test_no_key_is_bitwise_zero_trace_and_key_true_is_live` |
| 两门不互代偿 | `test_two_gates_are_distinct_and_not_substitutable` |
| 全量 | `tests/contract` **940 passed / 1 failed**（唯一失败=预存冻结锚，在 17F 基线内） |

新增契约测试文件 `tests/contract/test_run_opts_wiring.py`（8 条）。

⚠️ 该文件**未**写真调 `CasingFlowSolver(enable_stream_yield_gate=...)` 的反例：
`test_entrypoint_signatures.py` 的「豁免集冻结为**唯一**一处蓄意反例」是刻意收窄的守卫，
不应为一条断言扩容——改用 `inspect` 证明，效力等价。


---

## 5. 执行窗口对抗核查吸收（2026-10-06，7 路透镜 + 完整性批评者）

核查方法：6 路独立透镜（改动面完整性 / 关2 逐位风险 / P-1 口径算术 / run_opts 装配面 /
测试与回归面 / 既有资产复用）+ 1 名完整性批评者；**证伪优先**，每条须给 file:line 与实跑证据。
共 49 条发现（BLOCKER 0 / MAJOR 16 / MINOR 16 / INFO 17）。逐条处置：

| # | 发现 | 处置 |
|---|---|---|
| C-01 | **run_opts 键 `rheology_formula` 与 solver 形参 `rheology_formula_params` 不同名**，而 `*_kwargs_from_opts` 的返回 dict 被 `**` 直接展开 ⇒ 运行期 `TypeError`；仓库签名闸门 `check_call_signatures.find_bad_kwargs` 对"函数返回的 dict"是**盲区**（实测返回空） | **已修**：两 kwargs 函数内**显式改名**；新增 `test_kwargs_from_opts_keys_are_valid_solver_params` 断言返回键集 ⊆ `__init__` 形参集（补闸） |
| C-02 | **符号陷阱**：`a*T*T - b*T + c` 的 `-b` 若被搬成正数再代入 `+(-b)*T`，线性项符号静默翻转（实测 τy 偏 5.7×） | **已修**：字段存**带符号**系数，函数体一律 `q0*T*T + q1*T + q2`；新增 `test_field_values_are_bit_identical_to_original_literals`（struct.pack 位级）钉死 |
| C-03 | 计划/素材C 的「ρ≈1.94 → ρgh≈145.6」是**串井**算术（用呼1-004 井深 7660 m）；呼1-003 井深 7618 m 下 ρ=1.940 → 144.98 MPa；得 145.96 需 ρ≈1.9533 | **已更正**：§2.5 表以实测为准；并实测 **g=9.80665 → 145.63（−0.23%）／g=9.81 → 145.68（−0.19%）**，§2.1 明确 g 为可注入参数 |
| C-04 | 「呼1-003 井斜 ≤4°」不准：实测**最大 4.12°**（7450 m）；`inclination_profile.csv` 的 `tvd_m` 列**全空** | **已更正**（Δ4）；首点以上按 0° 的误差 ≤0.05%（评审实测） |
| C-05 | `verify_temperature_coupling.py` 的 `OUT_DIR` **硬编码**（`:81`）+ 正文硬编码路径（`:18/:44/:554/:667`），v1 已成 **tracked** 文件 ⇒ 就地重跑=覆盖既有产物（违红线） | **已列入 C2 改动面**：OUT_DIR 参数化（`--out-dir`）+ 正文路径自 OUT_DIR 派生 |
| C-06 | 该脚本 **6 处 `fluid_at`** 未分类：`:171/:176/:186`=关2「参照构造法」（须与 solver 内同 params，**且不得接 P**——关2 靠 P=None）；`:308/:391-396`=公式侧口径差表（须接 params + P） | **已列入 C2 改动面** |
| C-07 | 仓内**无 ADR 目录/模板**（`find -iname '*adr*'` 零命中） | **已列入 C2**：新建 `docs/superpowers/adr/ADR-0001-*.md` |
| C-08 | `pressure_mode/pressure_caliber` 加进 run_opts 后**无转发点**（kwargs 函数只收 `opts`，而 `pressure_field` 是需 well/fluids/schedule 的对象）⇒ 第二种死代码形态 | **已列入 C2**：落点 = `run_variant_res` 内 `solver_kw['pressure_field'] = build_pressure_field(...)` |
| C-09 | `_pressure_at` 无场时若返回 **0.0**（而非 `None`）会走 clamp 而非 `p_default` ⇒ 关2 语义漂移 | **已定**：无场 ⇒ **返回 `None`**（第一分支，早于任何计算） |
| C-10 | `memo` 键 `(fluid,T,P)` 未含 **params** 维 ⇒ 同实例换参档会吃陈旧派生对象 | **已修**：键扩为 `(fluid, T, P, params)` |
| C-11 | `insitu_column_density` 口径含混：全 schedule 会把**管内**流体（替钻井液/压塞液）计入环空柱 | **已裁（口径声明）**：默认 = **全泵注体积加权**（Q9/D5 口径③「在场相密度」不加角色限定），显式声明其为近似；实测两口径对呼1-003 差 0.02%（1.9495 vs 1.9492）。提供 `roles=` 供 P-2 细化。**该口径须用户追认** |
| C-12 | 2D 合成隔离液（呼1-004 ρ≈1.81，域外 borrow）叠真实 P = **双重外推**，且呼1-004 无鞋深压力锚 | **已列入 C2 验收**：保留 borrow 审计断言；口径差表**单列**呼1-004 合成相 τy/μp 增量 |
| C-13 | P-1 抬高隔离液 τy 经混合 τy 场灌入 `_yield_gate_wall`（`:1339`）的**耦合面无验收** | **已列入 C2 验收**：`pressure_mode` off vs hydrostatic 的 wall 活化率/冻结度差 + spacer τy 增量的口径差分解 |
| C-14 | `cemdisp/runners/*` **全 11 文件零温度接线**（`enable_temperature_rheology` 取默认 False）⇒ R1/P-1 在**权威 runner 路径上是死代码** | **声明**：R1/P-1 = **敏感性/探针口径，不进生产 runner**（§6-8：runner 产品化需另裁）。Phase 2 的"P-1 后口径"批走敏感性链，**不指望** runners 产出的权威目录 |
| C-15 | annulus **构造层无 geom**（geom 在循环内才建），无法用 `_representative_pressure_mpa(geom,t)` | **已定**：构造层用**井级标量** `P(鞋深 md, 0.0)`；静压场与 t 无关 ⇒ 与循环内逐步值同 memo 同对象（A3.1 成立） |
| C-16 | casing `_phase_props(self, fluid)` **无 geom/t**（与 annulus 三参签名不同） | **已定**：casing 侧 `run()` 内一次性算 `self._step_P = P(鞋深, 0.0)`，`_phase_props` 内取该标量；**不照抄** annulus 形态 |
| C-17 | `pressure_mode != 'off'` 而温度开关关 ⇒ 静默空转（同 `temperature_mode` 陷阱） | **已定**：`normalize_run_opts` 增交叉守卫，响亮报错 |
| C-18 | 现有关键测试 `test_temperature_phase_props.py:282-308` 硬钉 SPACER 派生 == `fluid_at(base, 60.0)`（P=None 口径） | **落实为红线**：无场 ⇒ `P=None` 是第一分支 |
| C-19 | `pressure_field.py` 应补 `__all__` 与 oob 三件套（与 `temperature_field.py` 数据层同型） | **已列入 C2 改动面** |
| C-20 | `tvd_m` 权威列在 hu101/ht1_004 的 CSV 中存在（比自积分准 3.15 m ≈ 0.04%） | **登记 Δ8**：本阶段统一走 `well.inclination_profile` 自积分（呼1-003 该列全空，必须积分），误差 ≤0.05%；CSV 直读列为 P-2 细化项 |
| C-21 | 关2 比对须对 `git show HEAD:...` 而非工作树 | **已执行**：`probe_equiv.py` 即为 HEAD 版 vs 新版的 5250 组合 + 2953 审计事件比对 |
| C-22 | `_step_T_c` / HB memo 前向注记 | **已记**：HB memo 温度维本阶段不动（HB 默认关）；若 R4 让水泥式吃 P，HB 键须同步加 P 维 |
| C-23 | memo 仅**步内**有效（每步 `clear()`）⇒ P 为标量无性能风险 | **已记**：P-1 不变量 = **每步单标量 P**；若 Phase P 需逐格 `P(z,t)`，memo 键须改粒度 |

**核查的独立确认（无需改动）**：`fluid_at` 生产调用点恰 6 处（annulus 5 + casing 1），无漏网；
`_phase_props` 唯一定义各 1 处；`_phase_memo` 无键结构断言；无测试 import 私有系数函数；
无第二处公式系数副本（无分叉）；两 solver 的 `__init__` 形参**全为 KEYWORD_ONLY 且无 VAR_KEYWORD**
⇒ 新 kwarg 不可能造成位置错位；`test_dead_switch_guard::test_signature_defaults_match_table` 只遍历
`_SWITCH_DEFAULTS` ⇒ 新 kwarg 不影响；ht1_004 的 TVD 双源（inclination CSV 的 `tvd_m` 与细井身结构 CSV
的 `vertical_depth_for_logging_m_`）逐位相同；`RUN_OPTS_KEYS` 确为全仓零引用的装饰常量。

**全程实测的关键数字**：R1 参数化 vs HEAD：**5250 组合 0 失配 + 审计事件 2953 条逐条全等**；
加法形式 210 温点 0 变位、Horner 109/210 变位（故禁 Horner）；
`test_default_path_bitwise_anchor`（关1）**1 passed**；`tests/contract/` 相关 7 文件 **185 passed**；
新增 R1 契约测试 **16 passed**。
