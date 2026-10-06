# Phase 2 设计规格：R2 屈服门接回 + μp/τy 真拆分（2026-10-06）

> **上游**：《详细执行计划_真温压响应_2026-10-06》§3 Phase 2 + §4 五关 + §6 先报后动 + §7 红线；
> 《评审A_R1-R4设计对抗评审_2026-10-06》D-02/D-03/D-05 + §4 A2.1/A2.2/A2.3；
> 《文献复核B2_屈服门正则化可辩护性_2026-10-06》（谱系定位与引用锚）；
> 《温压耦合真响应改进方案_2026-10-05》§R2。
> **前置**：Phase 0/0.0/1/1.5 已验收（HEAD `4b63ba4`）。评审 D-05 裁定 **P-1 必须先于 R2** —— 已满足。
> **基线**：全量 `pytest tests/ -q` = **17 failed / 1050 passed / 241 subtests**（2026-10-06 复跑 366.38s 实测）。
> **预检**：§8（行号 8 条 + 实质 5 条）。
> **对抗核查**：§10（8 agent / 133 万 token / 44 条发现：3 BLOCKER / 8 MAJOR / 18 MINOR / 15 INFO），
> 其中 **BLOCKER-1/2 推翻了本 spec 初稿的核心论断**（见 §2.3 与 §8 Δ-P2-3）。

---

## 0. 改动面总览

| 组 | 文件 | 内容 |
|---|---|---|
| **2a 拆分本体** | `cemdisp/models2d/buoyancy.py` | `fluid_apparent_viscosity(..., *, include_yield_term=False)`；`_YIELD_TERM_MU_CLIP_PA_S = 3.0` |
| | `cemdisp/models2d/annulus_d2dga.py` | ctor kwarg `include_yield_term`；4 个站点传开关；观测字段（§4） |
| **2b 管道** | `scripts/entrypoints/run_sensitivity_current_20260916.py` | `RUN_OPTS_KEYS` + `normalize_run_opts` + `annulus_kwargs_from_opts`（`casing_kwargs_from_opts` **禁出**） |
| | `scripts/entrypoints/run_sensitivity_temperature_t2_20261001.py` | `_opts()` **形参 + 返回 dict 双改**（见 §3 陷阱） |
| **2c 熔断** | `cemdisp/models2d/annulus_d2dga.py` | dt 时程 / α_cfl 裁剪计数 / `mu_reg` 场 / `λ_op` / `γ̇_rep`（**不进 summary**，§4） |
| **2d 测试** | `tests/contract/test_yield_term_split.py`（**新**）、`test_run_opts_wiring.py`（扩）、`test_buoyancy.py`（扩，已核无冲突）、`test_rheology_formula_params.py`（两端同改即自动通过） | §5 |
| **2e 驱动** | `results/_probe_屈服门四角_<日期>/probe_four_corner_<日期>.py`（**新文件**，按 §6-1 登记为计划改动面外） | §6 |
| **文档** | `docs/superpowers/adr/ADR-0002-wall-channel-primary.md` | wall 主通道 + 两族互为对照 |

**红线自检**：`_SWITCH_DEFAULTS:193-219` 一字不动；一切走运行侧 opt-in；关 1 锚读 JSON 比对；
预存 17F 勿动；源目录 `cement model\` 冻结；新产物目录带日期后缀。

---

## 1. 依据章（每条计算方法改动必须挂锚）

> 用户 2026-10-06：模型代码改动（尤其温压耦合相关）**必须有依据**（文献 + 其他资料）。
> 本章是本阶段全部改动的依据来源；**依据强度如实分级**，不足处显式声明。

### 1.1 本构层：μ_app = μp + τy/γ̇

| 锚 | 强度 |
|---|---|
| **仓内既有实现**（主依据）：`annulus_d2dga.py:1121` Bingham `mu = PV + τy/gamma`；`:1130` HB `mu = τy/gamma + K·gamma**(n-1)`；`:1133` `np.clip(mu,1e-5,3.0)`。已生产运行（T-off 全批） | **强（内部一致性）** |
| **现场资料**：六速表 `notes` 列逐字记 AV/PV/YP（如 "AV136mPa.s PV114Pa YP22Pa"；ht1_003/004 施工设计"一.1.3.4 低温性能"节同式） | **强（现场口径）** |
| **外部文献**：Bingham (1916) *An Investigation of the Laws of Plastic Flow*, Bull. Bur. Stand. 13:309；API RP 13D *Rheology and Hydraulics of Oil-Well Drilling Fluids* | ⚠ **本阶段素材库内无著录**（§10 E-3）⇒ **须补正式著录锚**，补前不得以之充当唯一外部依据 |

**结论**：R2「真拆分」= 把**标量口径**补成与**场口径**同一本构，消除"同一模型两套黏度口径"。
核心式的外部文献锚**偏弱**，如实声明：本构式依据 = 场口径既有实现 + 现场 AV/PV/YP 口径；外部文献仅支撑**正则化/上限**（§1.2 来源四）。

### 1.2 上限 3.0 Pa·s

**来源一（口径唯一性，主依据，强）**：`:1133` 的 `np.clip(mu, 1e-5, 3.0)` 是仓内首个提交 `0f78af9`
起的既有结构常数，喂 `mu_reg` → `_yield_gate_wall`。标量口径若另取正则化 ⇒ 出现**第三套黏度口径**。

**来源二（实测边界，支撑性）**：八井六速反算（Fann 35 R1B1：τ[Pa]=0.511·θ，γ̇[s⁻¹]=1.703·N）。
⚠ **2026-10-06 对抗核查（§10 E-1/E-2）发现原始 CSV 存在列/字段错位**（hu1 表头 13 列、12 行 14 字段；
ht1_004 第三方复检行的 `notes` 自带真实转速-读数对，与列名冲突）。**统计口径改为：notes 有
`<转速>:<读数>` 标注时以 note 为准重建，否则用列名；错位不可复位者剔除并注明。** 复算数见 §10 E-1 处置后回填。

**来源三（工作域不生效，强）**：拆分后的 μ_app 在**代表剪切率** γ̇_rep = 6⟨|w|⟩/⟨b⟩ 上取值。
按 q/A（三井井身结构 + 排量 0.55–1.5 m³/min）估 **γ̇_rep ≈ 63–185 s⁻¹**（**q/A 是 ⟨|w|⟩ 的上界**：
真实场含近静止格 ⇒ 实际更低、μ_app 更高，量级仍远低于上限）。
上限触发条件 γ̇ < τy/(3.0 − μp)：

| 流体 | 阈值 γ̇ |
|---|---|
| 组 B 水泥 @150 °C（μp=0.28848, τy=6.788） | **2.50 s⁻¹** |
| 钻井液式 @T clamp 80 °C（μp=0.046797, τy=10.9835） | **3.72 s⁻¹** |
| 组 A 水泥 @150 °C（μp=0.07952, τy=41.938） | **14.36 s⁻¹** |

⇒ 上限只在 **γ̇ ≲ 2.5–14.4 s⁻¹** 生效，即**六速最低转速 5.11 s⁻¹ 之下或刚跨过**的近静止端。
各地 μ_app 量级（γ̇∈[70,185]）：**钻井液式 0.11–0.20**、**组 A 水泥 0.31–0.68**、**组 B 水泥 ~0.31–0.31**
—— **均远低于 3.0**。

**来源四（文献，正则化必要性，中）**：γ̇→0 的 τy/γ̇ 奇异性必须有界正则化——Saramito & Wachs (2017)
*Rheol. Acta* 56(3):211–230；Papanastasiou (1987) *J. Rheol.* 31:385；Balmforth, Frigaard & Ovarlez (2014)
*Annu. Rev. Fluid Mech.* 46:121–146（"replacing the plugs with **very viscous** fluid"）；
Putz, Frigaard & Martinez (2009) *JNNFM* 163:62–77。谱系定位与 B 级评级见《文献复核B2》。
本仓另有 Papanastasiou 路径（`yield_regularization_M`，`annulus:2195-2200`）同族、不冲突。

**裁定**：上限 = **3.0 Pa·s**，**结构常数、禁作标定钮**（f_safety = 1.15 同型纪律）。

### 1.3 γ̇ floor 保持 1e-8 不动

`buoyancy.py:66` 的 floor 位于默认分支幂律路径之前；改 1e-6 会破关 2。
且有了 3.0 上限后 floor 不再起作用 ⇒ 保留原值，只在新分支加 clip。

### 1.4 数值健康度门槛（A2.2）的依据 —— **本阶段升级为载荷项**

评审 §4 A2.2「最危险的单一验收盲区」；**对抗核查实测坐实**：拆分轴 (T,F) 角 Δη_N 在
nz=100→140 摆动 **10 pp**（−34.60 → −44.15）。⇒ 熔断指标**必须落盘**，且四角报告须先过 nz 收敛检查。

### 1.5 通道共线量落盘（A2.3）的依据

评审 D-03：四组合数学上不可干净分解 ⇒ 报告写「**四角总量 + 交互残差**」；落盘每角
`(mu_reg 场, γ̇_rep, λ_op, wall 场)`。

---

## 2. 2a 拆分本体设计

### 2.1 `buoyancy.fluid_apparent_viscosity`（唯一改动点）

```python
_YIELD_TERM_MU_CLIP_PA_S = 3.0   # 结构常数；与 annulus._apparent_viscosity:1133 同值同义；禁作标定钮

def _clip_yield(mu):             # 与 annulus:1133 逐字对齐
    return min(max(float(mu), 1.0e-5), _YIELD_TERM_MU_CLIP_PA_S)

def fluid_apparent_viscosity(fluid, shear_rate, *, include_yield_term: bool = False) -> float:
    g = max(float(shear_rate), 1e-8)                      # §1.3：保留原值（关 2 红线）
    if fluid.rheology_model in (POWER_LAW, HERSCHEL_BULKLEY):
        if fluid.consistency_k is not None and fluid.power_law_n is not None:
            mu = float(fluid.consistency_k) * g ** (float(fluid.power_law_n) - 1.0)
            if include_yield_term and fluid.yield_stress_pa:
                mu = float(fluid.yield_stress_pa) / g + mu     # HB 补 τy/γ̇
            if fluid.plastic_viscosity_pa_s:
                mu = max(mu, float(fluid.plastic_viscosity_pa_s))   # 不变
            return _clip_yield(mu) if include_yield_term else mu
    if fluid.plastic_viscosity_pa_s is not None:
        mu = float(fluid.plastic_viscosity_pa_s)
        if include_yield_term:
            if fluid.yield_stress_pa:
                mu += float(fluid.yield_stress_pa) / g          # Bingham 补 τy/γ̇
            return _clip_yield(mu)          # ← 无条件 clip（§10 A-1：与 :1133 同构）
        return mu
    raise ValueError(...)
```

- `include_yield_term=False` ⇒ 两条 `return` 与 HEAD **位级相同** ⇒ 关 2 天然成立（§10 B-1 实跑 104 例 0 mismatch）。
- **clip 在 `include_yield_term=True` 时无条件生效**（含 τy=0 的 Bingham、含 HB 带 PV 的情形）——对抗核查 A-1/F-4 指出初稿把 clip 挂在 τy 判据下，与 `:1133` 不同构。
- **与场口径的已知差异（显式声明，不作掩盖）**：`_apparent_viscosity` 的 γ floor = **1e-6**（本处 1e-8），
  且 HB 分支**无** `max(mu, PV)` 地板。⇒ T2 的"逐位同值"**只在 γ̇ ≥ 1e-6 且 HB 流体 PV=None 时成立**。

### 2.2 solver 侧（4 站点**同传**）

| 站点 | 位置 | 所在方法 |
|---|---|---|
| `_froude_squared_at` | `:1515` | 方法本身（被 `:1894` 流函数段、`:2245` 旧路径、`:2655` run 循环 I3 块调用） |
| `_buoyancy_number_at` | `:1563` | 方法本身（被 `:2249` 旧路径、`~:2869` summary 调用） |
| η₁ | `:1885` | **`_velocity_stream_function`**（§10 C-3：初稿误标 `_compute_velocity`） |
| η₂ | `:1888` | 同上 |

ctor：`AnnulusD2DGASolver.__init__`（`:442` 邻近）加 `include_yield_term: bool = False`（keyword-only）。
**`__init__` 同时初始化** `self._lambda_op_last = None` / `self._shear_rate_rep_last = None`（§10 C-1：消除早退分支 AttributeError 面）。

⚠ **「4 站点同传」只在生产新路径（`enable_stream_function=True`）下可满足**：η₁/η₂ 两站点位于
`_velocity_stream_function`，旧代数路径不经过；旧路径下若开 True，只有 F²/b 变、η 不变 ⇒
**非对称**。⇒ 登记为"仅新路径支持；旧路径开此键为**非对称对照档**，须加注"。

### 2.3 拆分的作用面（**BLOCKER 更正版**）

> ⚠ 本节初稿断言「椭圆算子尺度在单位通量 BC 下完全抵消 ⇒ 流函数路径免疫」——**已证伪**（§10 L3-1/2、L6-1）。
> 错因：把 BC 当成事后归一化。实际 ∇·(αI₁∇Ψ) = −∇·(βb) 在 **Ψ(1)−Ψ(0)=1 的硬约束**下，
> Ψ = (β/α)Ψ₀ **违反 BC**，必须叠加修正项 ⇒ **形状改变**。

| 链路 | 对拆分的依赖 | 实测 |
|---|---|---|
| 流函数源 `b_field` | `λ_op/F²` **精确相消**（该不变量成立，实测偏离 ~1e-16）；但 `χ` 经 `m` 变 | b_field **非逐位不变**：maxrelΔ **1.43%**、meanrelΔ 0.20%（§10 L6-1）；χ 符号为 **正**，+0.05%~+0.11%（hu101，§10 L3-4） |
| **椭圆算子** `a_cell = 1/(2I₁)`，`I₁ ∝ 1/√(η₁η₂)` | **主导通道** | 隔离实验：冻结 b_field 只改 η₁ ⇒ **relmax\|Δw\| = 16–28%**（§10 L3-1） |
| **I3 浮力弥散通量** | `q ∝ f_φ ∝ 1/μ̂₁`（经 `:2655` 的 F² 站；**I3 的 η₂ 系数与 m 场是拆分不变量**） | 实测仅 **+0.6 ~ +1.5 pp**（η_N）——小量 |
| 诊断量（summary `浮力数_b`、旧代数路径 mobility/F²） | `∝ 1/μ̂₁` | 数值变化，属口径升级 |

**整机实跑（hu101，monkeypatch 精确复刻 §2.1）**：拆分单开 ⇒ **Δη_N = −44.15 pp（nz=140）/ −34.60 pp（nz=100）**，
Δη_E = −8.48 / −1.16 pp。⇒ **拆分轴是主导轴，量级 ≈ gate 轴的 3–4 倍、方向相反，
且带 10 pp 级网格摆动（A2.2 载荷）**。

⇒ §6.2 四角期望值据此重填；**跑前必须先做 nz 收敛与 dt 中位漂移检查**，确认非数值污染后
才可写方向性结论（§6-11）。

---

## 3. 2b 管道接线（Option B，用户 2026-10-06 裁定）

| 位置 | 改动 |
|---|---|
| `RUN_OPTS_KEYS` | += `"include_yield_term"`（第 9 键；**仅测试镜像用**，生产零引用） |
| `normalize_run_opts` 默认 dict | += `"include_yield_term": False` + `bool` 类型闸 |
| `annulus_kwargs_from_opts` | **仅当为 True 时出键**（沿用 `enable_stream_yield_gate` 的"默认不产生额外 kwarg"模式，§10 B-3，保既有契约测试 docstring 语义） |
| `casing_kwargs_from_opts` | **禁出此键**（`CasingFlowSolver` 无该形参；出键=运行期 TypeError） |
| T2 `_opts()` | ⚠ **形参 + 返回 dict 双改**：`include_yield_term: bool = False` **且** 返回 dict 增 `"include_yield_term": include_yield_term`。**只加形参不加键 ⇒ 静默回退 False ⇒ 四角 (T,·) 角实际按 off 跑**（§10 W-1，MAJOR） |

**新增契约测试（§5 W1）**：`normalize_run_opts(_opts("static", include_yield_term=True))["include_yield_term"] is True`
且 `annulus_kwargs_from_opts(该 opts)["include_yield_term"] is True` —— 把「标签=拆分 on ⇔ 真到达 solver」钉死。
`_opts` 为 keyword-only，全部既有 40+ 调用逐位安全（§10 W-4）。

---

## 4. 2c 数值熔断（关 5 实现）

新增观测（**不进 summary dict**；关 1 锚比对对象为 summary JSON ⇒ 位级成立，§10 B-2）：

| 量 | 落点 | 说明 |
|---|---|---|
| `dt_history: Tuple[float, ...]` | `AnnulusSimulationResult`（带默认值字段） | **三个收集点**：`:2548`（CFL-off 循环顶）、`:2604`（CFL-on 泵注）、`:2740`（CFL-on 泵停）。⚠ `:2606`/`:2742` 是 CFL-off 的**同值重复赋值**，**禁止**作收集点（否则双计，§10 A-2/B-5）。更稳：在循环体汇合处（`:2745` 前）**统一 append 一次** |
| `cfl_clip_events: int` / `cfl_clip_steps: int` | 同上 | I3 块 `:2694-2695`：`abs(div_q) > step_limit` 的**单元-步**计数与发生裁剪的**步**数 |
| `mu_reg_field: Array or None` | 同上 | 末步 `mu_reg` 场（定义在 `:2201`） |
| `lambda_op_last: float or None` | 同上 | `:1933`（**早退分支不更新** ⇒ 语义 = "最后一次流函数解的活跃步"，None 表示全程无解） |
| `shear_rate_rep_last: float or None` | 同上 | `:1883`，语义同 |

solver summary **零改动**。**不静默**：`run()` 末尾若 dt 中位数相对批次基线漂移 > 2×，发**一次性
`UserWarning`**（不进 summary ⇒ 不破关 1/关 2；§10 A-3 采纳方案 (a)，保留计划"不静默"语义）。
`extra_metrics`（敏感性层）仅**透传** res 新字段并派生 `dt_中位_s / dt_顶格步占比 / alpha_cfl_裁剪事件数 /
mu_reg_均值·max / λ_op / γ̇_rep / wall 分位`；落表走 Phase 2 **新目录的 sidecar JSON**（不改 T2
`CSV_COLUMNS`/`_row_of`，不破关 2 逐位；`EXTRA_SCHEMA_VERSION` 新目录下无需 bump，§10 X-4）。

**判据（A2.2）**：任一角相对 (0,0) 的 `median(dt_step)` 漂移 > 2× ⇒ 标「数值污染候选」，
**不得进方向性结论**。

---

## 5. 2d 专项测试（**逐条按 §10 更正**）

| # | 测试 | 更正后的判据 | 依据 |
|---|---|---|---|
| W1 | **run_opts 端到端**（新） | `_opts(..., include_yield_term=True)` ⇒ normalize ⇒ kwargs ⇒ **构造后 `s.include_yield_term is True`**；`False`/缺键 ⇒ 属性 False 且**不产生额外 kwarg** | §10 W-1/W-2 |
| T1 | 默认路径逐位 | 全族 × γ 网格，`include_yield_term=False` vs HEAD 实现**位级**一致 | 关 2 |
| T2 | clip 边界 + 场口径对齐（**域限定**） | ①flag=True ⇒ 输出恒 ∈[1e-5,3.0]，γ̇→0 收敛到 3.0；②与 `_apparent_viscosity` **逐位同值仅在「γ̇ ≥ 1e-6 且（Bingham 且 τy ≥ (3−PV)·1e-6，或 HB 且 PV=None）」**；③**负例单列**：Bingham τy=0 且 PV=5 ⇒ 两者均 3.0；HB 带 PV ⇒ **已知失配**，断言写在"已声明的差异"测试里 | §10 F-4/A-1/B-2 |
| T3 | Bingham/HB 手算对照 | `PV + τy/γ̇`、`τy/γ̇ + Kγ̇^(n−1)`（**τy 项在前，与 `:1130` 字面对齐**）与 numpy 手算逐位 | §1.1 |
| T4 | **拆成两条** | (1) **解析解对照**：复用 `tests/contract/test_gap_solver.py:55` 同构（`solve_fixed_G` + Bingham 单流体 + `u_exact` 逐位）；(2) **地板残余**：严格 `wall≡1` 单元（**排除过渡界面元**，界面元实测保留 26–36% 基线速度）上 `I₁_eff/I₁ ≤ 1e-5`；**不得**用"<0.1%/步"作单测判据 | §10 F-3 |
| T5 | **连续屈服门自解析极限**（替换 Balmforth） | 单列构造（仿 `tests/contract/test_yield_gate_continuous.py:46`）：τw < f·τy ⇒ wall→1；τw ≫ f·τy ⇒ wall→0；参考元恒 wall=0；**线性内插** `wall_i = clip(1 − τw_extrap/(f·τy), 0, 1)`。Balmforth 坡面式本模型无自由面，**不适用**，不得作 `_yield_gate_wall` 的解析真值 | §10 F-2 |
| T6 | **λ_op·F² 相消 + 算子不免疫（回归锚）** | (a) **乘积不变量**：4 站点同传放大 η ⇒ `λ_op·(1/F²)` **逐位不变**；(b) **负对照**：只传 F² 站（漏 λ_op 侧）⇒ b 相对变化 **≈ −55%**；(c) **回归锚**：同算例记录 `a_cell=1/(2I₁)` 与最终 `w` 的相对变化量级（实测 1.8–3.0× / 16–28%），把"不免疫"钉成锚。⚠ `b_field` 未入 `AnnulusSimulationResult` ⇒ 测试须 monkeypatch `solve_stream_function` 捕获其 `b_field` 实参。**不得**断言 `b_field` 逐位不变 | §10 L3-3/L6-1 |
| T7 | **网格收敛（判据重定义）** | "冻结过渡带" = `res.wall_field` 中 `0 < wall < 1−δ` 的单元计数（字段已存在）；**pin 网格阶梯**（ny,nz ∈ {(20,80),(30,120),(40,160)}）且 gate-on/off **分别**做；判据 = **两级间 \|Δη_N\| 单调递减**（趋势断言），**不写固定 <1%**（实测同实现下可 +24% → +0.22%，固定阈值不能区分） | §10 F-5 |
| T8 | 两门不互代偿 | 生产端 vs 消费端语义断言 | Phase 1.5 |
| T9 | f_safety ±15% | 1.0/1.15/1.3 三档呼101（**批量**） | D2-A |

**测试红线（§10 A-6）**：新测试**不得**写求解器构造的 `pytest.raises(TypeError)` 坏键反例
（会扩容 `test_entrypoint_signatures.py:283` 的冻结豁免集 `{('test_no_invented_dispersion.py',57)}`）；
证明 casing 无该形参照 `test_run_opts_wiring.py:112-113` 用 `inspect.signature`。

---

## 6. 2e 跑批（四角 + 三井 gate）

### 6.1 统一口径
`T-on` · `static`（T(z)=16.006+1.7598e-2·z）· `rate×1.0` · `pressure_mode="hydrostatic"` ·
`pressure_caliber="shoe"` · `CORRECTED_KW` · `nz=250` · CFL 自适应 · `tt=min(泵总+1200, stop_t)`。
> **已核对（§10 X-3）**：呼101 在 off / hydrostatic 下 η_E/η_N **逐位相同（Δ = 0.0）** ⇒ §6.2 的
> P = 0.1 MPa 参考值在 hydrostatic 口径下同样成立。呼1-003/呼1-004 **无此保证**，各自加跑一次 off vs hydrostatic。

### 6.2 四角矩阵（3 井 × 4 角 = 12 run）
轴 1 = 拆分 `include_yield_term` ∈ {F, T}；轴 2 = 消费端 `enable_stream_yield_gate` ∈ {F, T}。

| 角 | 拆分 | stream_gate | 呼101 期望（**已按 §10 L3-1/2 重填**） |
|---|---|---|---|
| (F,F) | off | off | 基准；η_N 应复现 **0.7777300041326185**（探针 A 已逐位验证） |
| (T,F) | on | off | ⚠ **算子通道主导**，实测量级 **Δη_N −35 ~ −44 pp**（nz 100/140）；**先过 nz 收敛与 A2.2 才写方向性结论** |
| (F,T) | off | on | wall 通道，探针实测 **Δη_N −10.4975 pp** |
| (T,T) | on | on | 两通道同开（**交互残差**单列，禁写"分解"） |

### 6.3 三井 gate 验收批（6 run，其中 3 run 与 (T,F) 复用）

| 井 | 角色 | 判据 |
|---|---|---|
| **呼101** | L1 主判 | 显著非零 + 窄边滞留↑ + 窄边前缘 −500 m 量级（探针 −545 m） |
| **呼1-003** | 首建 | wall 占比未知——如实记录；≈0 则效应 ≈0 属正常 |
| **呼1-004** | 交叉 | gate 轴 ≈ 0（探针 −0.0008 pp）；**拆分轴 ≠ 0**（§10 X-3 提示）⇒ 分列 |

### 6.4 追加批
- **关 1**：T-off 锚，读 `results/_baseline_T_off/基线摘要_T_off.json` 比对（禁手抄）。
- **关 2 排量锚回归（§10 X-2）**：呼101 **gate-off** ×0.6/×1.0/×1.4 复跑，与
  `results/_probe_排量敏感性_2026-10-02/探针结果.json`（0.8264188630 / 0.7777300041 / 0.7301758287）
  **近位对账**（CONTEXT.md「排量响应锚」是每阶段必测项，初稿只排了 gate-on）。
- **排量锚（gate-on 新观测）**：呼101 ×0.6/1.0/1.4 —— §6-5 检查点：**若破坏单调性 ⇒ 停、如实报告**。
- **f_safety ±15%**：1.0/1.15/1.3（呼101）。

**机时**：≈ 24 run；探针实测呼101 78–82 s / 呼1-004 47 s ⇒ 估 **35–85 min**。

### 6.5 产物
- `results/_probe_屈服门四角_<日期>/`：驱动脚本 + 四角判别量 + 共线量（A2.3）+ 健康度（A2.2）+ nz 收敛表。
- `results/敏感性补跑_Phase2_<日期>/`：三井 gate 批 + 排量锚（gate-off/gate-on）+ f ±15%。
- 报告：三井 gate 报告；**ADR-0002**。
- **措辞红线**：只写「**四角总量 + 交互残差**」，**禁写「贡献分解」**；wall 通道**禁写"精确满足屈服"**。

---

## 7. 五关验收

| 关 | 判据 |
|---|---|
| **关 1** | T-off 与 `基线摘要_T_off.json` **读 JSON 比对**（η_E 0.9996663176608406 / η_N 0.9983246031489231）；预存 17F 不增 |
| **关 2** | 缺键/False ⇒ 不产生额外 kwarg 且与 T2 对应行**逐位**；T1 位级；`_SWITCH_DEFAULTS` 未动 |
| **关 3** | 呼101：gate-on ⇒ 窄边滞留↑ / η_N↓；拆分通道 ⇒ 记录主导链路（算子）与量级 |
| **关 4** | 四角总量 + 交互残差；T2–T5 锚；T6 相消+不免疫锚；T7 趋势收敛 |
| **关 5** | A2.2：任一角 dt 中位漂移 > 2× 或 **nz 摆动 > 2 pp** ⇒ 标「数值污染候选」，不得进方向性结论；红灯而关 3 绿灯 ⇒ 停 + 报（§6-11） |
| **两锚** | 排量响应锚（**gate-off 回归 + gate-on 新观测**）；**L1 流变响应达成判定** |

---

## 8. 预检：计划偏差

### 8.1 行号（按符号检索复核，基准 HEAD `4b63ba4`）

| # | 计划/评审写法 | 实测 | 性质 |
|---|---|---|---|
| A1 | `buoyancy.py:54-75` | 吻合 | 准确 |
| A2 | 调用点 `:1820/:1823`、`:1450`、`:1498` | `:1885/:1888`、`:1515`、`:1563` | 漂移 +65；站点数 4 正确 |
| A3 | 场口径 clip `:1108` | Bingham 公式 `:1121`、**HB 公式 `:1130`**、clip `:1133`；γ floor **1e-6**（`:1110`） | 漂移 |
| A4 | 唯一扯断点 `:2155` | `:2220` | 漂移 +65 |
| A5 | `_yield_gate_wall :1211-1269` | `:1276-1333` | 漂移 |
| A6 | `stream_function.py:421-435` | `if wall is None:` `:430`，块 `:419-436`；地板 `:182` | 漂移 |
| A7 | `_SWITCH_DEFAULTS:192-220` | `:193-219` | 基本吻合 |
| A8 | `_dead_switches:234-262` | `:235-289` | 上界 +27 |

### 8.2 实质偏差

- **Δ-P2-1**（计划缺口，已裁定）：「拆分」无 run_opts 键名 ⇒ 新增 `include_yield_term`（Option B）。
- **Δ-P2-2**（改动面漏列，已裁定）：A2.2/A2.3 观测字段无现成落盘位 ⇒ 新增带默认值字段（§4）。
- **Δ-P2-3（初稿论断被证伪，已改写）**：初稿称「λ_op 与 F² 精确相消 ⇒ 流函数路径免疫、唯一作用面是 I3」。
  **(a) 对的部分**：`λ_op·(1/F²)` 乘积确实逐位不变（实测 ~1e-16）。
  **(b) 错的部分**：由此推出"算子免疫"是错的——算子 `a_cell = 1/(2I₁)` 随拆分变 1.8–3.0×，
  实测速度场变 16–28%、整机 Δη_N −35~−44 pp。**错因 = 把单位通量 BC 当成事后归一化**（§10 L3-1/2）。
- **Δ-P2-4（作用面排序改写）**：初稿称 I3 是"唯一实质作用面"。实测 I3 仅 **+0.6~+1.5 pp**；
  主导是**椭圆算子通道**。I3 的 η₂ 系数与 m 场是拆分**不变量**，作用面只在经 `:2655` F² 站的 `f_φ`。
- **Δ-P2-5**：计划 §3 Phase 2e「呼1-004 交叉（预期≈0）」只对 gate 轴成立；拆分轴 ≠0。
- **Δ-P2-6（新增，§6-1 登记）**：2e 四角**驱动脚本为新文件**（`results/_probe_屈服门四角_<日期>/*.py`），
  越出计划改动面清单 ⇒ 已登记（§0）。
- **Δ-P2-7（新增，§6-1 登记）**：观测字段 + 运行期一次性告警，与计划 §3 2c「告警**进 summary**」
  不同（改为 warning + 关 5 gate），理由 = 保关 1/关 2 位级红线。

---

## 9. 回退

2a/2b/2c 全为**加法式**（新参默认=旧行为、新字段带默认值、summary 零改动）。
`git revert` 单 commit 即回 HEAD 口径；新产物目录可整目录删除。ADR-0002 与 spec 为纯文档。

---

## 10. 对抗核查吸收表（2026-10-06，8 agent / 133 万 token / 44 条）

**方法**：7 路独立透镜（改动面完整性 / 关 2 逐位 / 相消算术 / run_opts 装配 / 熔断可实现性 /
测试可证伪性 / 依据面）+ 1 名完整性批评者；**证伪优先**，每条须给 file:line + 实跑证据。
汇总：**BLOCKER 3 / MAJOR 8 / MINOR 18 / INFO 15**。

### L3 相消算术（BLOCKER ×2）
| # | 发现 | 处置 |
|---|---|---|
| L3-1 | **BLOCKER**：`I₁ ∝ 1/√(η₁η₂)` 并未在单位通量 BC 下抵消；冻结 b_field 只改 η₁ ⇒ 速度场变 16–28%（隔离实验）+ 一致缩放实验 relmax\|Δw\| = 0.482/0.743 | **已改写 §2.3/§8 Δ-P2-3**；算子列标为主导通道 |
| L3-2 | **BLOCKER**：(T,F) 角期望写反：实测 Δη_N −44.15 pp(nz140)/−34.60 pp(nz100)，非"小量" | **已重填 §6.2**；加 nz 收敛前置 |
| L3-3 | MAJOR：T6 只测 b_field（真），测不到算子不抵消，会绿着放行错误前提 | **T6 已改**：乘积不变量 + 负对照 + 不免疫回归锚 |
| L3-4 | MINOR：χ 符号写错（实为 **+0.05%~+0.11%**） | 已更正 |
| L3-5 | INFO：I3 倍率应逐井现算（呼101 实际 ÷1.9–3.3，非通用 ÷2.3–4.4） | 已改为"按井现算" |

### L6 测试可证伪性（BLOCKER ×1 / MAJOR ×4）
| # | 发现 | 处置 |
|---|---|---|
| L6-1 | **BLOCKER**：T6 的 `b_field` 逐位不变**被证伪**（maxrelΔ 1.43%，源自 χ(m)） | T6 改判据（见上） |
| L6-2 | MAJOR：T5 牛头不对马嘴（`_yield_gate_wall` 自 2026-09-15 起是连续口径，无二值分支；Balmforth 坡面式与本构造不同构） | **T5 已替换**为连续门自解析极限 |
| L6-3 | MAJOR：T4 不可执行（仓内塞流解析解在 gap_solver 路径；"地板残余"实测界元 26–36%、全带 8.6–12%；总通量由 BC 守恒） | **T4 已拆两条** |
| L6-4 | MAJOR：T2 在 HB 带 PV、Bingham τy=0 两类输入上失配 | **T2 已加域限定 + 负例**；实现侧 Bingham 改无条件 clip |
| L6-5 | MAJOR：T7 判据欠定义，实测随网格对摆动 +24% ↔ +0.22%，固定 <1% 无区分力 | **T7 已改趋势断言 + pin 网格** |

### L1 改动面完整性（MINOR ×2 / INFO ×4）
| # | 发现 | 处置 |
|---|---|---|
| A-1 | MINOR：Bingham 分支 clip 被挂在 τy 判据下，与 `:1133` 不同构 | **已修**（§2.1 无条件 clip） |
| A-2 | MINOR：dt 赋值实为 **5** 处；`:2606/:2742` 不得作收集点 | **已注明**（§4） |
| A-3 | INFO：初稿"不进 summary"与计划"告警进 summary"相左 | **已登记 Δ-P2-7**；采纳一次性 UserWarning |
| A-4 | INFO：wiring 测试只自动覆盖"合法性"，不覆盖"到达" | **已加 W1**（§5） |
| A-5 | INFO：2e 驱动脚本未列 | **已登记 Δ-P2-6**（§0/§6.5） |
| A-6 | INFO：新测试不得写 `pytest.raises(TypeError)` 坏键反例（污染冻结豁免集） | **已列测试红线**（§5） |

### L2 关 2 逐位（MINOR ×3 / INFO ×2）
| # | 发现 | 处置 |
|---|---|---|
| B-1 | MINOR：HB 分支"次序与 :1131 一致"措辞与代码矛盾（=`:`1130`，且书写次序相反；IEEE-754 下等价） | **已改写 §2.1**（τy 项在前，字面对齐 `:1130`） |
| B-2 | MINOR：T2"逐位同值"在 γ̇<1e-6 与 HB 带 PV 时不成立 | **已加域限定**（§2.1/§5 T2） |
| B-3 | MINOR：`annulus_kwargs_from_opts` 恒出键与既有"新键默认值不产生额外 kwarg"不变量冲突 | **已改为仅 True 时出键**（§3） |
| B-4 | INFO：旧代数路径"4 站点同传"物理上不可满足 | **已声明**（§2.2） |
| B-5 | INFO：dt 双计风险 | **已注明**（§4） |

### L5 熔断可实现性（MINOR ×2 / INFO ×1）
| # | 发现 | 处置 |
|---|---|---|
| C-1 | MINOR：`_lambda_op_last` 在 `q≤0` 早退分支不更新；类无初始化 ⇒ AttributeError 面 | **已定**：`__init__` 置 None + 语义声明（§2.2/§4） |
| C-2 | MINOR：改 `extra_metrics` 的 schema 版本与落表方式未定 | **已定**：新目录 sidecar JSON，不动 T2 CSV（§4） |
| C-3 | INFO：`:1885/:1888` 在 `_velocity_stream_function`（非 `_compute_velocity`）；λ_op `:1933`、mu_reg `:2201` | **已校正**（§2.2/§4） |

### L7 依据面（MAJOR ×1 / MINOR ×6 / INFO ×2）
| # | 发现 | 处置 |
|---|---|---|
| E-1 | **MAJOR**：§1.2 三个头条实测值由 **CSV 列错位**产生（呼1-004 第三方复检行 notes 自带真实转速对；hu1 字段数不齐），"3.0 低于个别流体 6 rpm 实测值"不成立 | **已声明清洗口径**（§1.2 来源二）；重算数**待回填** |
| E-2 | MINOR：有效读数 389 不可复现（独立重算 399） | **已声明按显式口径统计**；数待回填 |
| E-3 | MINOR：Bingham (1916) / API RP 13D 在素材库**无著录** | **已标"须补锚"**（§1.1），并如实声明核心式外部锚偏弱 |
| E-4 | MINOR：Balmforth 引文漏 "viscous"（语义反转） | **已更正**（§1.2 来源四） |
| E-5 | MINOR：HB 行号 `:1131` → `:1130` | **已校正**（§8.1 A3） |
| E-6 | MINOR：阈值区间漏组 B 的 **2.50**（下界应为 2.5 非 3.7） | **已补齐三档表**（§1.2 来源三） |
| E-7 | MINOR：三重点井 3 rpm 最大为 **1.0002** 非 1.00；μ_app 0.14–0.16 只是钻井液（水泥 0.31–0.68） | **已更正**（§1.2 来源三） |
| E-8 | INFO：q/A 是 ⟨\|w\|⟩ 上界，下界实测 **63**（非 70）；方向应向更低 γ̇/更高 μ_app | **已披露**（§1.2 来源三） |
| E-9 | INFO：八井六速表列数/列义不齐 | **已声明**（§1.2 来源二） |

### L8 完整性批评者（MINOR ×3 / INFO ×3）
| # | 发现 | 处置 |
|---|---|---|
| X-1 | MINOR：§2.3 把 I3 作用面引到 `d2dga_flux.py:123-125`（该处 η₂ 与 m 是拆分不变量），与 §8 自相矛盾 | **已更正**为 `annulus:2655`（经 F² 站）（§2.3） |
| X-2 | MINOR：CONTEXT.md「排量响应锚」要求每阶段 **gate-off** 三点回归，初稿只排 gate-on | **已补**（§6.4/§7） |
| X-3 | INFO：§6.1 强制 P-hydrostatic，而 §6.2 参考值来自 P=0.1 MPa 探针；仅因呼101 对压力完全不敏感才一致（实跑 Δ=0.0） | **已声明**（§6.1）；另两井各加一次 off-vs-hydro 检查 |
| X-4 | MINOR：反驳 L5 的"必须 bump EXTRA_SCHEMA_VERSION"（Phase 2 用新目录，无旧产物可复用） | **已采纳**（§4） |
| X-5 | INFO：用户约束②风险——核心本构式只有内部一致性依据，外部文献锚偏弱 | **已如实声明**（§1.1） |
| X-6 | INFO：`tests/contract/test_buoyancy.py` 已核**无冲突**（默认 NEWTONIAN，断言 `mu_displaced==PV` 在拆分默认关下仍绿） | **已登记**（§0/§5） |

### 实测关键数字（全部来自本轮对抗核查实跑）
- 隔离实验：冻结 b_field 只改 η₁ ⇒ relmax\|Δw\| = **0.277/0.206/0.158** @ γ̇=70/123/185。
- 整机（hu101，T-on static，monkeypatch 复刻 §2.1）：拆分单开 Δη_N = **−44.15 pp（nz=140）/ −34.60 pp（nz=100）**；
  分解：算子 −46.68/−35.60 pp、I3 **+0.6~+1.5 pp**。
- b_field 相对变化：maxrelΔ **1.433e-02**、meanrelΔ 1.978e-03（全部来自 χ(m)）；`λ_op·(1/F²)` 不变量偏离 ~1e-16。
- χ 变化：**+0.05% ~ +0.11%**（hu101，符号为正）。
- 默认路径位级：104 例 **0 mismatch**。
- 呼101 off vs hydrostatic：**Δη_E = Δη_N = 0.0**。
- 现状既有测试：`test_run_opts_wiring.py + test_rheology_formula_params.py` = 24 passed；
  `test_dead_switch_guard + test_buoyancy + test_run_opts_wiring + test_entrypoint_signatures` = **75 passed**。
