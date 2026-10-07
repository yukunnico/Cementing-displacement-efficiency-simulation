# Phase 4（逐深温度）改动面测绘调研（2026-10-07）

> 基线 = `feat/temperature-coupling @ 2ab0514`（Phase 3.1 之后）。总纲 Phase 4 原行号基于 ca8bba7，全部已重定位。
> 本调研只读；行号均为当前 HEAD 实测（`file:line` + 符号 + 一行摘录）。
> **全局现状**：8 个生产 runner **均未开温变流变**（`cemdisp/runners/*.py` 无 `enable_temperature_rheology` 传参），
> T-on 链路目前只活在 scripts/entrypoints（`run_sensitivity_current_20260916.py` / `run_sensitivity_temperature_t2_20261001.py` / `verify_temperature_coupling.py`）。
> Phase 4 的"生产波及面"= 这几个 entrypoints + 契约测试；runner 接线本身是独立决定。

---

## 4a（1D `cemdisp/transport1d/casing_flow.py`）

### 1. `_step_T_c` 单点定标（原 :339）

- **当前位**：`casing_flow.py:315` `self._step_T_c: float | None = None`（init）；`:319-320` 构造层恒温场回退；
  **定标点 `:362`** `self._step_T_c = self._temperature_field.T(well_spec.top_md_m, 0.0)`（run() 内，域顶 t=0，一次定标非逐步）；
  消费点 `:557` `T = self._step_T_c`（`_phase_props` 内）；`:363-368` P-1 同型标量 `self._step_P_mpa = pressure_field.P(shoe_md, 0)`。
- **核实**：成立（仍是单点标量语义）。变化点 = Phase 1 加了 P 维（`:313/:365-368`），4a 存 field 引用逐深/逐时查询时 **T 与 P 须一并处理**（`_step_P_mpa` 现为鞋深单点）。场对象引用其实已保存在 `self._temperature_field`（`:314`），改动位是"消费方式"而非"取不到场"。
- **波及面**：`run(temperature_field=...)` 形参（`:334`）已有注入通道；无外部读 `_step_T_c` 的 casing 侧代码（仅 `:557`）。
- **风险/建议**：4a 若改 `T(md,t)` 逐深，`_step_T_c` 标量对 HB memo/审计的语义（见项 5）在 1D 侧无 HB 消费，风险小；`_step_P_mpa` 是否同步逐深取决于 P-1 口径（现静压场与 t 无关，逐深=P(md)，改动成本低）。

### 2. `_phase_props` 现签名与调用点（原 :514-535）

- **当前位**：**`:543`** `def _phase_props(self, fluid: FluidSpec) -> FluidSpec:`；memo 键 **`:560`** `key = (fluid, T, P, self.rheology_formula_params)`；派生 `:564-566` `fluid_at(fluid, T, P, params=…, mud_extrapolate=…)`。
- **全部调用方（1D 内仅 3 处，模块外无调用——私有方法）**：
  1. `:601`（`_compute_dispersion_coefficient` 弥散入口）
  2. `:674`（`_effective_viscosity` 有效黏度入口）
  3. `:1351`（`_gravity_corrected_arrival_time` 屈服应力侧）
  （另 `:175` 注释对比 annulus 三参签名、`:1288+` 消费。）
- **核实**：签名**未变**（仍单参 `fluid`，annulus 侧才带 geom/t）；原计划改 `(fluid, md_m, t_s)` 的波及面 = 上述 3 调用点 + 键 `:560` + `:557-558` 取温/取压 + 项 4 新增 2 处包点。
- **风险/建议**：1D 的调用点在时间线**构建期**（非逐步循环），md/t 由体积链反推（项 3），签名扩为 `(fluid, md_m, t_s)` 后 memo 键从 `(fluid,T,P,params)` 变 `(fluid,T(md,t),P(md),params)`——1D 每 run 派生次数≈事件数（数十），无爆炸问题；键里放连续浮点 T 前建议量化/就近取档（沿用 Q16 口径）。

### 3. 三入口现状（原弥散 :568/:783、有效黏度 :641/:715-716、重力 :1316）

- **弥散**：`_compute_dispersion_coefficient` def **`:570`**（_phase_props@`:601`，Fan&Wang τy 消费 `:624/:629/:633`）；被 `_apply_dispersion_to_timeline`（def **`:763`**）在 **`:816`** `D_eff = self._compute_dispersion_coefficient(pipe_radius_m, fluid, U)` 调用。
  **体积链反推 (time, depth) 可行性 = 成立**：事件时刻 `event.time_s`（`:858` `t_arrival = event.time_s`）；`scheduled_steps` 携 `cumulative_volume_start/end_m3 + start/end_time_s`（`:102-103`），`_cumulative_volume_at`（`:522/:1043`）可逆算 `md_front(t) = V(t)/pipe_area`（活塞流）；`_inject_start_time`（`:1205`）与 `_contact_time_integrated_sigma`（`:1243-1261`，σ_t=√(2·D_eff·t_contact)/U）已按体积链配对 t_inject/t_arrival——逐深化只需把"单点 D_eff×总接触时间"扩成"D_eff(T(md(t),t)) 路径积分"，数据全成。
- **有效黏度**：`_effective_viscosity` def **`:652`**（_phase_props@`:674`）；其**唯一消费方**是 `_interface_instability_factor`（def `:709`）的 `:748-749`，而该函数在 **`:737-738` `if self.has_plug: return 1.0` 短路在前**。instability 调用点仅 `:848`（弥散带增强）。
  **核实：死路仍成立**——8 井 runner 全 `has_plug=True`（`runners/ht1_001:364、ht1_003:309、ht1_004:327、hu101:196、hu102:305、hu103:308、hu1:347、hu2:353`）⇒ `_effective_viscosity` 的 T-hook 在生产口径永不执行。4a spec 应显式声明该入口"仅探针/非胶塞井可达"，勿把它写进生产逐深承诺。
- **重力**：`_gravity_corrected_arrival_time` def **`:1288`**；鞋深 = `well_spec.shoe_md_m`（经 `_effective_pipe_radius_m` 等）、到达时刻 = 形参 `arrival_time_s`（`:1290`）——"现成"成立；`_phase_props`@`:1351`、τy 消费 `:1352-1358`。调用点：`:449`（`_ordered_front_arrival_times` 内）、`:990`（`_build_shoe_timeline` 内）；另有 `:521` 停泵判据旁路。
- **风险**：弥散/重力两个活入口的 (md,t) 各自含义不同（弥散=沿程路径、重力=鞋口终点单点），4a spec 要分开口径，不要共用一个"前缘深度"定义。

### 4. B8 两路 raw τy 包 `_phase_props`（原 :1431-1432/:1444-1451、:1571-1573/:1580-1588）

- **路径①（legacy 重力）**：`_legacy_gravity_correction` def **`:1441`**；raw 取流体 `:1464` `fluid = next((f for f in fluids if f.name == fluid_name), None)`；raw τy 消费 **`:1477` `if fluid is not None and fluid.yield_stress_pa is not None …`** 与 **`:1482` `yield_ratio = min(fluid.yield_stress_pa / tau_critical, 1.0)`**——**无 _phase_props 包**，核实成立。
- **路径②（停泵增强沉降）**：`_settled_exit_fluid_name_enhanced` def **`:1559`**；raw 取流体 `:1595` `current_fluid = next(…)`；raw τy 消费 **`:1614`/`:1618`**——**无包**，核实成立。
- **对照**：现代重力路径已包（`:1351`，Phase 1 成果），现代弥散/有效黏度已包（`:601/:674`）。
- **波及/建议**：包法仍是"取到 fluid 后补一行 `fluid = self._phase_props(fluid)`"（`:1464→1465` 后、`:1595→1596` 后），各 ~1-5 行；注意 4a 改签名后这两处也要给出各自 (md,t)（legacy 重力=鞋深/arrival_time；停泵沉降=停泵时段）。**触发面小**：路径①仅 `enable_buoyancy_physics=False` 时进入（`:1327-1329`），路径②仅停泵段（`:1603-1605`）——包它是口径一致性修补，不是数值主链；spec 应把它们与逐深 T 的优先级分开。

### 5. HB memo 温度维 + 逐列化键爆炸评估（原 annulus :1949-1953）

- **当前位（2D）**：`_hb_closure_for` def `annulus_d2dga.py:2034`；基础键 **`:2060`** `key = (n1, n2, k1, k2, tau_y1, tau_y2)`；**T 维 `:2061-2063`** `if self.enable_temperature_rheology: key = key + (self._step_T_c,)`（Phase 1 已落地，标量代表温度）；缓存读写 `:2064/:2072`；`self._hb_closure_memo = {}` 建于 **`:777`**，**整个 run 不清**（grep 无 clear）。
- `_phase_memo`（派生流体缓存）键 `:1287` `key = (fluid, T, P, self.rheology_formula_params)`；按 run 清（`:2510`）+ **按步清（`:2616`）**。
- **memo 命中路径数（2D `_phase_props` 调用点）**：`_compute_props` 4 相（`:1392-1395`）、`_velocity_stream_function` 3 相（`:1913-1915`）、run 循环每步重派生 4 相（`:2618-2621`）、诊断 4 相（`:3026`）⇒ **~15 次/步/运行路径**，步内去重后每步 4 个唯一键。
- **唯一 T 量级估算**：生产 `nz=250`（`runners/ht1_004:86`、`ht1_003:72`、`hu101:118`）；T-on 表格场每步代表温度（域均值）≈ 一个新值 ⇒ 现键上界 ≈ 4 流体×步数（total_t/dt，CFL 自适应 dt≈亚秒级，量级 10³–10⁴ 步）≈ **10⁴–4×10⁴ 派生对象/run**。
  **逐列化后**：键 (fluid, **T 列向量**) —— `np.ndarray` **不可哈希**，现键结构直接失效；若拆成逐列标量键 ⇒ 上界 = nz×步数×4 流体 ≈ **2.5×10⁶ 级**，HB 闭包更重（每个 `HBClosure(ny_gap=201)` 带查表，`:2070-2071`，逐列重建=250×步数，不可行）。
- **对策建议（供 spec 裁定）**：
  1. **静温类场**（Constant/Geothermal/AnchoredProfile，T 与 t 无关）：逐列派生**每 run 一次**，缓存 `(fluid, col) → FluidSpec`（4×250=10³ 对象，步间零重建）。
  2. **表格场**：按表时间节点预计算参数数组 `(fluid, col, t_node)`（4×250×200=2×10⁵，一次/run），逐步按 t 索引/插值——键问题消解为数组索引。
  3. `_hb_closure_memo` 需改为 per-column 闭包数组（或 LRU）；若 Phase 4 不打算开 HB 路径逐列（生产 hb_closure 默认关，`_SWITCH_DEFAULTS:214`），可声明"HV memo 仍按代表标量 T 键、HB 路径不逐列"作过渡。
  4. **审计溢出**：`rheology_vs_temperature.py:26` `fluid_at` 审计定长 10000（超限丢最旧，接线后按步清空）——逐列派生每步 1000 次 `fluid_at` 必然溢出 ⇒ 审计改抽样或仅记录 clamp/borrow 事件，spec 必须处理。

---

## 4b（2D `cemdisp/models2d/`）

### 6. `_representative_temperature` + 均匀场快捷路径（原 :1147-1169）

- **当前位**：def `annulus_d2dga.py:1207`；**返回仍是标量 float**。
- **均匀场快捷路径（原 :1166-1168）现为 `:1225-1228`**：
  `:1225` `values = np.asarray(arr, dtype=float).ravel()` → `:1226` `first = float(values[0])` → `:1227` `if np.all(values == first):` → `:1228` `return first`；`:1229` `return float(np.mean(values))` 为非均匀均值分支。
- **核实**：**关 2 红线位置确认**——逐列化改造时 `:1226-1228`（均匀⇒取该值、免均值浮点漂移）必须原样保留（Constant 场下 `np.mean` 与 `values[0]` 不保证逐位相等，此短路是"与手写同参参照自洽逐位"测试的实现要件）；另缺省回查分支 `:1217-1224`（geom 无 "T" 键⇒域中部现查）。
- **建议**：逐列化后该函数的归宿 = 保留为 F²/b_num/`_step_T_c` 标量口径的**代表值**来源（与逐列场并存），或退役并由列场均值显式定义——spec 需裁定并守住关 2 测试。

### 7. annulus 侧 `_phase_props` 列版派生波及面（原 :1192 前后）

- **当前位**：def **`:1266-1295`** `def _phase_props(self, fluid, geom, t)`；内部 `:1282` `T = self._representative_temperature(geom, t)`、`:1283` 写 `_step_T_c`、`:1284` P、`:1287` 键、`:1291-1293` `fluid_at`。
- **全部调用点**：`:1392-1395`（`_compute_props`）、`:1913-1915`（`_velocity_stream_function`，不派生 spacer）、`:2618-2621`（run 循环逐步重派生，含 spacer）、`:3026`（诊断 `_temp_rheo_fluids` 收派生后流体，T1-2 裁定 `:3020`）。
- **同族旁路（不经 memo 的构造层直调）**：`run()` 构造层 `:2513` `t_c = field.T(top_md, 0)` + `:2519-2532` 四次 `fluid_at(…, t_c, p_c, …)`；geom 初始化 `:922` `geom["T"] = self._temperature_array(md, 0.0)`；每步刷新 `_refresh_geom_temperature` `:1203-1205`（调用位 `:2498/2614` 注释链）。
- **核实**：列版签名（返回 `(nz,)` 列参数而非单 FluidSpec）波及 = 12 个调用点 + 4 个构造层直调 + 诊断暴露语义（`_temp_rheo_fluids` 现为标量派生 4 相）+ 键结构（项 5）。
- **风险**：`_phase_props` 是"物性唯一入口"的设计支柱（消灭双算，注释 `:47-48/:1270-1272`）——列版若绕过它各自现算，将破坏"构造层=逐步派生"（评审 A3.1）不变量；spec 应保持单入口并把列场做成该入口的返回值形态。

### 8. `_compute_props` 场公式吃列参（原 :1292-1295/:1321-1327）

- **当前位**：def **`:1371-1456`**。吃流体参数的场公式（全部已是 numpy 广播，缺的只是"流体参数从标量变 (nz,) 列"）：
  - `:1400-1406` `mu = mud * self._apparent_viscosity(mud_fluid, gamma)`（四相黏度场）
  - `:1407-1413` `rho` 密度场（当前密度无 T 依赖）
  - `:1421-1427` `tau_y` 场（`_phase_cement_tau_y` `:1297` 水泥相公式 τy + `_fluid_yield_stress`）
  - `:1432-1443` `m_field = mu_mud/mu_cement` 场 + `np.clip(0.1,10)` `:1443`
  - `:1445-1446` `eta1=mu_mud_field, eta2=mu_cement_field`（**已是 (ny,nz) 场返回**）
  - `:1449-1455` `n_mix/κ_mix`（体积加权 n、对数加权 K，标量参数×场权重）
- **核实**：原假设成立——场算子侧几乎零改动，改动集中在上游"每相一套标量参数 ⇒ 每列一套参数"；`_apparent_viscosity/_phase_power_law_params/_fluid_yield_stress` 需接受列化参数（或其调用改传列 FluidSpec 数组）。
- **风险**：`gamma = 6|w|/b`（`:1399`）与列参数组合时注意 `max(k,1e-12)` 类标量防御（`:1453` `max(k_mud,1e-12)` 是 **python max 对 float**——列化后必须改 `np.maximum`，否则数组比较报错）；逐列 clip 的 `0.1/10.0` 口径保留。

### 9. 标量端 m/η₁/η₂ 逐列化位置 + Phase 2 四站点交互（原 :1818-1826）

- **当前位（`_velocity_stream_function`，def `:1777`）**：shear_rate `:1921-1923`（γ̇=6⟨|w|⟩/⟨b⟩，全场标量）；**eta1 `:1924-1926`**、**eta2 `:1927-1930`**（`fluid_apparent_viscosity(…, include_yield_term=self.include_yield_term)`）、**m_ratio `:1932`** `eta1/eta2`；`λ_op ∝ eta1` `:1973-1979`；`mobility_i1/i2` 消费 `:1952-1955`。
- **Phase 2 `include_yield_term` 四站点当前行号**：**`:1552`**（`_froude_squared_at`）、**`:1601`**（`_buoyancy_number_at`）、**`:1926`**、**`:1930`**（均在 `_velocity_stream_function`）；init `:462`/self `:674`。
- **核实/交互**：逐列化只作用于 `:1924-1932`（→列场，走项 10 数组化）；`:1552/:1601` 属**全域代表标量**口径（F²/b 进 (4.14) 标量动力学）——逐列后这两处的 μ̂₁ 定义（取鞋口列？域代表？）是**必须裁定的口径**，且四站点同传约束（Phase 2 拆 μp/τy）要求列化时同一 `include_yield_term` 语义贯穿，不能一处列一处标量各吃一套参数。
- **风险**：η₁ 同喂 `λ_op`（`:1975`）——逐列化后 λ_op 是 (nz,) 乘 (2,ny,nz) b_field（广播可行），但 λ_op×F² 的 μ̂₁ 精确相消不变量（Phase 2 注释 `:1976-1978`）在逐列后只对"同列同参数"成立，spec 需重申该不变量防归因错读。

### 10. `two_layer.mobility_i1/i2` 数组化可行性（原 :57-129）

- **当前位**：`mobility_i1` def **`two_layer.py:57-91`**、`mobility_i2` def **`:94-128`**（行号几乎零漂移）。`c_bar` 已场化（`:86/:123` `np.asarray`）；**`m/eta1/eta2/H` 中仅 H 是场**；**`math.sqrt(m)` `:87`/`:124`**、**`math.sqrt(eta1*eta2)` `:90`/`:127`** 为数组化唯二卡点。
- **核实/评估**：`math.sqrt` 对 float64 与 `np.sqrt` 同为 IEEE-754 正确舍入平方根 ⇒ **对标量逐位同值**，替换数值风险≈0（建议契约测试钉死标量档逐位）。真实风险在返回类型分支 `:91/:128` `float(out) if np.isscalar(c_bar) else out.astype(float)`：若 m 数组化而 c_bar 标量，`float(out)` 在数组上直接 TypeError——重排分支（按"任一输入为数组⇒场"判定）；形参注解 `m: float`（`:59/:96`）放宽。
- **建议**：改动极小（4 行 sqrt + 2 处返回分支 + 注解），先行做"数组化但全场等值"回归（逐位=HEAD）再放真列参。

### 11. `solve_stream_function` 接受数组 eta1/eta2 的现状（原 :324-326/:421-435）

- **当前位**：def **`stream_function.py:324-326`**（零漂移）`solve_stream_function(geom, c_bar, eta1, eta2, m: float, b_field, *, closure, wall, ny, nz, banded)`；docstring `:334-335` 已写"标量或全场一致数组"。
- **`_scalar_viscosity` `:227-236`：是——强制标量**：`np.allclose(arr, first, rtol=1e-12)` 否则 raise（`:234-235`），消费点 **`:396-397`**；`m = float(m)` **`:398`**。
- **b_field 场入算子先例 = 确认**：(2,ny,nz) 全向量，形状检查 `:409-416`；annulus 装配 `annulus_d2dga.py:1959-1961`（φ/ξ 槽）+ `:1979`（λ_op 乘子）——算子本体（a_cell/c_cell/调和面系数 `:438-447`）**已全程 (ny,nz) 逐点广播**，非均匀流动度（wall 通道 `:428-433` I₁_eff=I₁·conductance）已在生产路径证明场化可行。
- **结论**：数组化 eta1/eta2 的硬卡点只有三处——`:396-397`（放宽为真场：把逐点 μ 传给 mobility 闭包/直接算 I₁ 场）、`:398`、`two_layer` sqrt（项 10）；`ClosureProvider.mobility(self, c_bar, m, eta1, eta2, H)`（`:577-586`，含 `NewtonianClosure`/`PowerLawGapClosure`/`HBClosure`）全部签名 `m/eta*: float` 需同步放宽；**HB 非线性路径** `_rheology_from_closure` `:758` + `hb_closure.mobility(c, m, e1, e2, H)` `:785` 是最大遗留（查表闭包按标量 m/B 构造，逐列 ⇒ 列版闭包族，见项 5）。
- **风险**：`_scalar_viscosity` 的"禁止静默平均"注释是 Z&F22 (4.21) 每流体单黏度前提的守卫——放宽即离开原闭式适用域，spec 须写明"列版两层闭包 = 逐列独立 (4.21) + 列间耦合仍由椭圆算子"的合法性论证。

### 12. 4c 同源化验收位：呼101 鞋口 1D/2D 代表温差 21.7°C（原 4c 断裂）

- **1D 取值点**：`casing_flow.py:362` —— `T(well_spec.top_md_m, 0.0)`（**域顶、t=0、全程不变**）。
- **2D 取值点**：`annulus_d2dga.py:1225-1229` —— 每步 `geom["T"]` 剖面：**均匀→`values[0]`、非均匀→深度均值**；构造层 `:2513` 与 1D 同（top_md,0），入循环后改口径。
- **核实（算术复现）**：成立。呼101 域 `top=5400 / bottom=7868`（`hu101_loader.py:79-80`），静温档 Geothermal 线 `16.006+0.017598·md`（`temperature_field.py:112-113`）：1D=111.04°C（5400），2D=132.75°C（均值 md=6634），差 = 0.017598×1234 ≈ **21.72°C** ✓。断裂就在"1D 域顶单点标量"与"2D 循环内深度均值"两语义之间。
- **建议**：4c 同源化的验收探针应三方对账：`casing_flow:362` vs `annulus:2513`（构造层，本已同源）vs `annulus:1225-1229`（逐步）；逐深化后前两者被列场取代，断裂自然消解，但**关 2 的 Constant 场逐位红线测试锚在 `:1226-1228`**，4c 改口径不得顺手"归一化"掉均值/首值之分。

### 13. `AnnulusInletState.temperature_c` 占位（T1-4 接口位）

- **当前位**：`boundary_bridge.py:62` `temperature_c: float = float("nan")`（doc `:12`）。
- **核实**：仍为**纯占位**——全仓 grep 无任何消费方（`cement_yield_stress.py` 的同名字段是出处标注字符串，无关）。4c/4d 走场注入路径，不经 inlet state ⇒ **不冗余**，spec 声明"保留占位"或明确接线即可。

---

## 4d（`cemdisp/data/temperature_field.py`）

### 14. 场类清单 / 协议 / `AnchoredProfileField` 新增位 / oob 三件套

- **当前类**：`ClampEvent:70`、`ConstantTemperatureField:82`、`GeothermalTemperatureField:116`（常量单一真源 `:112-113`，注释含 09-30 §1 裁定与 333 点拟合残差 ≤0.05°C）、`TableTemperatureField:160`、`load_delivered_pair:444`。
- **协议**：**无 formal Protocol 类**——duck 接口 `T(md_m, t_s)->float` + 审计三件套（`oob_count`/`oob_events`/`reset_audit`，Constant `:92/:97/:100`、Geothermal `:141/:146/:149`、Table `:232/:237/:241`）；对照 `pressure_field.py:73` 已有 `class PressureField(Protocol)`（4d spec 可顺手补对称 Protocol，成本低）。
- **`AnchoredProfileField` 新增位**：`GeothermalTemperatureField` 类尾 **`:155/157` 之间**插入（静温同型零变查询；若锚点含静止/循环双族或分段非单调，需决定是否启用真 oob 审计——三件套形状照抄 Geothermal+Table 混合）；导出加 `__all__ :35-48`。
- **风险**：oob 审计为**逐查询记账**（Table `:256-263` append 每条）——逐深×逐步×多井消费会爆内存/列表，spec 给去重（按 (md_bin,t_bin) 聚合计次）或定长环策略；`load_delivered_pair:482-488` 的井底行同值校验是 Table 独有约定，Anchored 类无对应校验 ⇒ 锚点冲突（见项 16 静态/循环两族）需在类内显式裁定或抛错。

### 15. 呼1-004 温度表时程 / `stop_t>表末` 现状断言位

- **实测（`参考文档/温压耦合数据、/T_out.npz` 缓存）**：shape **(333, 200)** ✓；`time_s` 0 → **11940.0 s**（=199×60，col0=初始时刻，**200 列≠12000 s**——轴末是 11940 s）；深度轴 30→7660 m；全域温幅 16.53–150.76°C。常量：`TIME_STEP_S=60.0 :60`、`EXPECTED_SHAPE (333,200) :59`。
- **stop_t 关系**：`runners/ht1_004_tailpipe.py:333/:346` `total_t_s = annulus_stop_time_s(...)`；现口径 stop_t ≤ 11940 s——该断言不在 cemdisp（库内仅 clamp+审计 `temperature_field.py:253-263`，越界**静默延伸**），而在探针脚本：`scripts/entrypoints/verify_temperature_coupling.py:684`（"时间轴 0–11940 s 覆盖 stop_t"）+ `:686/:903-904`（oob_in/oob_out 计数，0=查询全在表域）。
- **>11940 s 扩展需求**：Phase 4 若把查询推到停泵后/候凝延时（或 1D 管容链时刻>表末），现行为=clamp 到末列 + 逐条 oob 审计 ⇒ spec 需（a）在 4d 加显式 `t>表末` 门（raise 或声明口径），或（b）把表尾延拓（末列常值）写进 Anchored/Table 语义，且与 `EXPECTED_SHAPE` 硬校验（`:184-188`——扩表列数即失配报错，需同步放开）。

### 16. 三重点井 `temperature_pressure_profile.csv` 锚点可用性（静温锚点对账 probe 数据源）

实测行统计（`参考文档/现场资料提取/<井>/temperature_pressure_profile.csv`）：

| 井 | 数据行 | 带温度行 | 域内温度锚点（去重 md） | 备注 |
|---|---|---|---|---|
| ht1_003 呼1-003 | 8 | **5** | 0(60)、5290(123)、5305(152)、6500(152)、7618(150) | 另 3 行纯压力/当量（7463/6880/7618）；5 个唯一域深 |
| ht1_004 呼1-004 | 12 | 8 | 0(60)、5241(124∥155)、6600(155)、7660(155×2) | 5241 **静止/循环同md两值**；含 2 行邻井（7746@150=HT1-001、7618@152=HT1-003，medium）；域内唯一 md=4 |
| hu101 呼101 | 7 | **3** | 5400(117)、5700(123)、7868(150) | 5700/5400 为回接段电测；另 4 行纯压力/漏层 |

- **核实**：总纲"呼1-003 锚点最富"**部分成立需修正**——按带温度行数 004(8)>003(5)>101(3)，但 004 的 8 行含 2 邻井 + 同 md 双族重复；按域内唯一 md 数 003=5≈004=4-5>101=3。**口径陷阱**：三井 CSV 混装 出口循环/悬挂器电测/实验静止/邻井 四类温度族，`AnchoredProfileField` probe 必须先逐井裁定吃哪一族（建议与 09-30 §1 统一地温线对账时只用静止族，notes 列含出处可机读）。

---

## Phase 4 spec 必须处理的冲突清单

1. **逐列化 × memo 键结构（1D 小、2D 致命）**：`(fluid,T,P,params)` 键的 T 一旦变列向量即不可哈希；`_hb_closure_memo`（`:777`，run 内不清）逐列重建不可行 ⇒ 采用"静温类逐列一次派生 + 表格类按表节点预计算参数数组"的缓存形态（项 5），并保留步内去重语义（"物性唯一入口、消灭双算"不变量）。
2. **逐列化 × 标量代表温度语义链**：`_representative_temperature`（`:1207-1229`）的均匀场短路 `:1226-1228` = 关 2 红线，须原样保留；但 F²/b_num/`include_yield_term` 四站点中的 `:1552/:1601` 仍要标量代表值 ⇒ 逐列后"代表标量怎么取"是新增口径决定（不能默认沿用均值，21.7°C 教训=语义差会被读成物理效应）。
3. **`_phase_props` 单入口 vs 构造层旁路**：2D 构造层 `:2513-2532` 直调 `fluid_at`（P-1 裁定"构造层=逐步派生同 memo 键"）——列版入口若改返回值形态（列场），构造层四直调、诊断 `_temp_rheo_fluids`（`:3026`）、契约测试 `test_temperature_phase_props.py`（钉 `fluid_at(base,60.0)`）同批重排。
4. **1D 三入口中"有效黏度"是死路**：`has_plug` 短路（`:737-738`）× 8 井全 has_plug ⇒ 该入口的 T-hook 生产不可达；spec 不得把三入口并列承诺"全部逐深化生效"，应标注可达性。
5. **`two_layer`/`stream_function` 数组化守卫链**：`_scalar_viscosity`（`:227-236`）的"全场一致否则 raise"与 Z&F22 每流体单黏度前提冲突——放宽前要有逐列独立 (4.21) 的合法性论证；`math.sqrt→np.sqrt` 逐位安全但返回类型分支（`two_layer.py:91/:128`）会炸。
6. **HB 非线性闭包未列化**：`HBClosure/hb_groups` 按标量 (m,B,τy1,τy2) 构造（`annulus:2053-2071`、`stream_function:758-785`）——若 Phase 4 不逐列 HB，需显式声明"HB 路径保持标量代表口径"，否则 T2/Phase2e 的 HB 变体与逐深 T 互斥。
7. **审计溢出**：`fluid_at` 审计定长 10000（rheology_vs_temperature.py:26）+ Table oob 逐查询 append（temperature_field.py:256-263）——逐列×逐步必然溢出 ⇒ 抽样/聚合/按步清，三处一起改。
8. **表末 11940 s ≠ 12000 s**：时间轴末=199×60=11940 s；`EXPECTED_SHAPE=(333,200)` 硬校验（`:184-188`）挡扩展；>11940 s 现仅静默 clamp+审计，库内无断言（断言活在 verify 脚本 `:684`）——4d 需把"stop_t/总时长 vs 表末"变成正式门或声明口径。
9. **锚点 CSV 混合温度族**：静/循环/出口/邻井四族同表混装（项 16），`AnchoredProfileField` probe 的选族规则不裁定 = 制造第二个"21.7°C 型"语义断裂。
