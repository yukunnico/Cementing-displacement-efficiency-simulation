# 环空顶替模型失真修正：设计规格（M0–M4 + 附加项）

> **日期**：2026-08-23
> **状态**：设计稿，待用户审阅后决定实施单元
> **范围**：只覆盖环空段失真修正（A 方案），不含 Track B 深化（椭圆流函数横流、TVD、瞬态、Uzawa 精确屈服）
> **决策记录**：①范围=A（M0/M1/M3/M2/M4 + 2 附加项）；②CFL=方案乙（主结果 CFL 自适应开，动手前先做 on/off 对比存档）；③重跑节奏=合并制（M0 一轮 / M1 一轮 / M3+M2 一轮 / M4 收尾一轮，开发期用 hu101 三剖面快检）；④实现路线=路线2（增量开关，新机制默认关）+路线3 轻量基准；⑤实施顺序=**M0→M1→M3→M2→M4→附加项**（2026-08-23 核查修订8：M4 单独开会重演 529:1，须在 M3 屈服门槛就位后再放开 e_clip_max；开发期 M3 用默认 e=0.55 调参）
> **代码基线**：commit 2832a08；`annulus_d2dga.py`/`d2dga_flux.py`/`displacement_metrics.py` 与 HEAD 一致；`boundary_bridge.py`/`casing_flow.py` 有未提交改动（不涉及本计划对象）
> **依据**：紧邻四篇 08-23/08-20 诊断报告 + 代码全景调研（C1）+ obsidian 约束调研（C2）+ 7 篇文献全文精读（Maleki 2017/2019、Foroushan 2020、Escudier 2002、Roustaei 2015 P1/P2、Jung&Frigaard 2022）+ 主会话流态计算。所有 file:line 以 C1 代码调研为准。

> **⚠️ 代码核查修订（2026-08-23，4 路独立核查后）**：以下条目修正原稿与代码不符之处，**以本节为准**。动手实施前必须先消化：
> 1. **【I3 会崩，必改一行】** `d2dga_flux.py:136` 是内建 `max(eta2, 1.0e-9)`，对 numpy 数组抛 `ValueError: truth value of an array is ambiguous`（已实测复现）。原稿§8"C1 已核 numpy 广播 OK"是**错的**。I3 局部化必须先把该行改为 `np.maximum(eta2, 1.0e-9)`；I3 公式本体（:110-111）不动。
> 2. **【扩元是 3 处不是 1 处】** `_compute_velocity` 8 元组扩元要同步改：`annulus_d2dga.py:805`（泵注解包）、`:933`（停泵解包）、外加原稿漏列的 **`tests/test_yield_deadzone.py:209`**（精确 8 元组解包，会炸）。
> 3. **【M3 停泵"自然冻结"不成立】** γ̇ 在 :590 用 `w_prev`，停泵期 `w_prev` 不更新且保留上一步非零值 → 外推 G≠0，会在停泵期把冻结元误判为可流动。原稿风险表"q=0→G=0→自然冻结（C1 已核无冲突）"**与代码不符**。必须二选一：① wall 更新显式门控在 `pump_active`（:801）内（推荐，停泵期保持上一步 wall）；② 改用当前步 w（但 q=0 时 wall 会全 1，改变 wall_field/`mean_wall_mud` 输出语义，不推荐）。
> 4. **【M3 宽边参考元规则必须显式定义】** "用流动中的宽边"不够：若 row 0（宽边）自身冻结 → τ_w,ref=0 → G=0 → 整片永久冻结无法解冻。参考元须取"该深度 |w| 最大且 >0 的元"（或 w>0 的最宽元），并对全切片无流动元时回退为全冻结。
> 5. **【M3 无屈服流体的质变须兜底】** 泥浆/水泥 YP≈0 时 `τ_w ≤ f·0` 恒假 → wall 层整体消失、b³ 放大全开。原稿"微量残泥下限"必须在实现前定形（建议保留一个小的几何/体积残泥底，或对 τ_y_local 设数值下限），否则是行为反转而非参数微调。
> 6. **【M3 振荡风险】** 每步重算 wall + pref→w→τ_w→wall 反馈可能逐帧翻转；现有 wall 作参数传入天然滞后一步，保留这一步滞后或加欠松弛，勿在同一步内迭代到收敛。
> 7. **【M2 漏列隐藏工作量】** 代码只有分相表观黏度 + 体积分数加权混合 μ，**无混合物局部 n/κ_p**，而 Metzner-Reed Re/He（式58/59）必须要它们。需新增"从混合 μ(γ̇) 对数斜率拟合局部 n/κ"逻辑；固定点迭代 2-3 轮的收敛/欠松弛准则也未定义。原稿"2-3 天"偏乐观，按 4-6 天估。
> 8. **【M4 排序裁定】** 决策记录⑤"M4 提前"与 §12 风险表"M4 排在 M3/M2 之后"矛盾。**裁定：M4 在 M3（屈服门槛）就位后再开**——M4 单独开、M3 未就位时 e=0.90 会重演 A4 的 529:1 拉爆。实施顺序：M0→M1→M3→M2→M4→附加项。开发期用 e_clip_max=0.55（默认）调 M3，不得把 e=0.90 结果当结论。
> 9. **【测试类名/行号】** 要重写的测试实际类名是 **`TestStaticWallLayer`**（`tests/test_improved_d2dga_annulus.py:557`），不是原稿写的 `TestCminCriteria`。各断言行号微漂：`test_wall_consistency_after_run` 实际 :571-596、`test_constructor_has_cmin_default` :560-564、`test_cmin_0_3_more_wall...` :624-647、`test_wall_zeros_velocity_in_wall_cells` :598-622。
> 10. **【"只加键不清键"措辞修正】** runner 层（8 个 tailpipe）是**直接索引**既有中文键（`result.summary["最终结果"]["全井段最终有效顶替效率"]` 等），不是 `.get`——改键会 KeyError，加键安全。`extract_ablation_metrics` 才是全 `.get`（改键静默返回 None，不炸）。reporting 层不读 summary，只对 depth_profiles/metrics 做"必需列存在性"检查，加列安全、删列炸。核心安全性结论不变。
> 11. **【M0 Tier0 except 勿收窄】** 原稿"收窄 try/except"会把静默跳过变成 run() 崩溃。已改为：保持宽 `except Exception`，但把错误写入 `result.summary["tier0_diagnostics_error"]`（2026-08-23 已落地）。M0 新键计算放在 try 块之外。

---

## 1. 目标与成功判据

**目标**：让环空段模拟"不失真"——恢复套管居中度对顶替效率的真实敏感性（从"死区 +0.32%"变为"阈值型响应"），并让主指标能看见窄边窜槽（η_N）。**不追求**把 8 井效率拟合到 CBL（口径结构性错配，见 §12）。

**成功判据**（三层）：
1. **机制可复现**：每个改造单元开启后，hu101 三剖面回归的 η_N 随均居中度单调；between_cent 不再出现 A4 的 −1.41% 反物理反转；mixing 从 0.59 降至 0.2–0.35（M1）。
2. **文献量级锚定**（M-F19 精确靶）：e=0.3 → η_N≥95%、e=0.6 → η_N≈30–35%；e=0.6/0.8 最宽/最窄 1/4 体积比 3.25/6.15。
3. **零回归**：所有未开启单元的既有 summary 键数值与基线逐位一致（证明"关=与 2832a08 完全一样"，归因可信）。

**诚实边界（C2 约束，写入验收原则）**：不得主张三闭包显著提效（既有消融 R3−R0≈−1.6pp）；不向 66.7%/67% 规范值校准（它是工程惯例，Jung&Frigaard 模拟显示 standoff 25% 仍可能有效顶替）；CFL 消融"数据集②（−1.56pp）"仅存 review 未落盘，不得作为依据。

---

## 2. 方案结构：开关体系（路线 2 骨架）

所有改动集中于 `AnnulusD2DGASolver`（`cemdisp/models2d/annulus_d2dga.py`）。每个机制独立 `__init__` 开关，**新机制默认关**，生产 runner 显式开启并做消融：

| 机制 | 开关参数（拟名） | 默认值 | 含义 |
|------|----------------|--------|------|
| M1 弥散重标定 | `dispersion_dt_scale` | `1.0` | =1.0 时与现行行为一致（系数乘 `dt_step/4.0` 的开关） |
| M4 e 护栏 | `e_clip_max` | `0.55` | =0.55 时与现行行为一致 |
| M3 屈服门槛 | `enable_yield_gate` | `False` | 开启后取代 c_min 冻结判据 |
| M2 流态修正 | `enable_regime_split` | `False` | 开启后启用 Maleki 2017 闭包 |
| 附加 I3 局部化 | `enable_local_i3` | `False` | eta2/Δρ 局部化 |

配套改造：
- **只加键不清键**：summary/depth_profiles 新增键，不动既有中文键（`最终结果` 等）；runner、`ht1_004_ablation.extract_ablation_metrics`、reporting 层均已核实只依赖既有键或 `.get` 语义，改键会炸、加键安全（C1）。
- **`_compute_velocity` 返回值扩元**：现返回 8 元组（:673），扩为含 `tau_y`、`eta2`（可选 `G` 场）的元组，一处透传供 M3/M2/I3 复用，避免重复计算。

---

## 3. M0 指标与口径层（输出层，约 0.5–1 天；零物理扰动）

| # | 改动 | 落点（C1 已核） | 内容 |
|---|------|----------------|------|
| 1 | η_N 进 summary | `:1124-1130` 收窄 try/except + `:1079-1096` 构建区 | 内联调用 `displacement_metrics._narrow_quarter_efficiency`（:211-225；无循环依赖已核实），写 `summary["最终结果"]["窄四分位效率"]`+英文别名，与 η_E 并列 |
| 2 | 失稳指数去饱和 | `:988-991` | 原 `instability_index` 保留；新增 `失稳指数_线性`（proxy 原值）+`失稳指数_对数`（log10(1+proxy)） |
| 3 | 评价窗效率 | `well_spec.py:70-93`；`solver.run()` 已持有 well_spec | 按 **2832a08 四口径 `window_type`**（cbl/cbl_quality/cbl_digitization/formation_target/model_focus/oil_gas_show/custom，已废除 target）枚举遍历 EvaluationWindow，每窗输出 b 加权 2D 积分的 η_E 与 η_N（不用 depth_profiles 单线平均）；呼101 cbl 窗 5699.8–7810、ht1_003 数字化窗与官方口径分开、ht1_001 只用 model_focus 不得冒充 cbl |
| 4 | 低尾指标 | `geom["standoff"]`（:325）可用 | `standoff低于0.5段占比` + `窄边效率低于0.05域占比` |

**验收**：①既有全部 summary 键数值与基线逐位一致；②hu101 三 case η_N=0.057/0.090/0.180 复现；③呼101 cbl 窗效率可算且"窄边拉开 25 倍（0.034→0.850）"论文论点不破坏；④`test_displacement_metrics`/`test_regime_classifiers` mock 驱动安全；`TestExtractAblationMetrics` 用 `.get` 安全。

**测试**：新增 `test_m0_metrics`（η_N 进 summary + 四口径窗遍历 + 低尾键）。

---

## 4. M1 弥散重标定（约 1–2 天标定，不换格式；TVD 属 Track B 后置）

**现状**：`_smooth_dispersion`（:410-440）每步固定幅值；调用点 `:860-863` 硬编码 lead/tail=0.018/0.015、spacer/flusher=0.012/0.012；CFL 自适应 dt 4→0.118s，单位物理时间弥散放大 ~34 倍 → mixing≈0.59 钉死基线（C2 确认硬编码名与名义 dt 标定时点）。

**改法（分步，最小改动优先）**：
```python
# 步骤A（先恢复量纲正确性）：
axial_eff   = self.dispersion_axial   * (dt_step / dt_ref)   # dt_ref=4.0（标定时名义 dt）
azimuth_eff = self.dispersion_azimuth * (dt_step / dt_ref)
lead   = self._smooth_dispersion(lead,   axial=axial_eff,   azimuthal=azimuth_eff)
tail   = self._smooth_dispersion(tail,   axial=axial_eff,   azimuthal=azimuth_eff)
spacer = self._smooth_dispersion(spacer, axial=axial_eff*0.667, azimuthal=azimuth_eff*0.8)  # 保持 spacer/flusher 与 lead 的 0.012/0.018 相对比例
flusher= self._smooth_dispersion(flusher, axial=axial_eff*0.667, azimuthal=azimuth_eff*0.8)
# dispersion_axial=0.018 / dispersion_azimuth=0.015 / dt_ref=4.0 均提为 __init__ 参数；dispersion_dt_scale 控制开关
```
- 保持 spacer/flusher 与 lead/tail 的相对系数比（0.012/0.018≈0.667、0.012/0.015=0.8），不引入新的无依据常数。
- 只在校正后（:857-871 之间）泵注分支调用；停泵分支不调用（:931-952，C1 已核）。

**验收（轻量基准）**：
1. 固定物理时长的 κ 扫描（scale∈{0, 0.25, 0.5, 1.0}）看 mixing/η_E/channeling 收敛性；预期 mixing 0.59→0.2–0.35；
2. 网格双收敛（nz 140/250、ny 20/40）：mixing 平台值应由 κ 主导而非 dt_step（证明已随 dt 归一）；
3. 前沿长度：从"半个环空"回到米级（现场参照：石油钻采工艺 2021 界面长度 2–5m）；
4. 既有消融表口径（ht1_004_ablation.py:104 现为 CFL 关固定 dt）——重跑时与主结果口径统一声明。

**论文数字影响（C2 明确，M1 提交说明必附"受影响数字对照表"）**：表4（8 井效率/窜槽/混浆）、表5（消融）、表6（网格收敛）、67% 居中度阈值、呼101"65.04% vs 62.77% +2.27pp"（孤儿数字）、"20 秒"卖点（4 处）。**M0 与 M1 分开提交**。

**风险**：M1 同向移动 8 井 η_E 与既有消融结论；CFL 口径（P0-1）未裁定前不能落地——本设计按方案乙裁定，动手前先跑一次 CFL on/off 对比存档。

---

## 5. M4 e 护栏（约 0.5 天，几何层单行）

**现状**：`e = np.clip(1.0 - standoff, 0.05, 0.55)`（:303），e=0.55 硬天花板使 standoff≤0.45 全部饱和；下游 e 用于 :314（b 分布）、:636（fallback 浮力分支）、:730（depth_profiles 偏心度列）。

**改法**：
```python
e = np.clip(1.0 - standoff, 0.05, self.e_clip_max)   # e_clip_max 默认 0.55，生产跑道显式设为 0.90
```
- 提为 `__init__` 参数 `e_clip_max`（默认 0.55 = 现行）。
- **不做 b_narrow ≥5mm 下限**（hu101 实测窄边 8.4mm 本就 >5mm，下限会掩盖真实风险）。
- 依据：Pelipenko 几何全域 e∈[0,1)、Z23 实验到 0.8、Foolad 到 1.0；A4 证伪"0.55 截断是解药"但完全放开 0.999 会在 M2/M3 未就位时重演 529:1 拉爆。

**验收**：开放 0.90 后 hu101 三剖面几何真正变化（between_cent 下部 7055–7868m 段 e 从 0.55→0.78），M3/M2 才能在该几何上调参；`_diag` 校核体积守恒（scale :330-336 在 e 变化后仍正确）。

---

## 6. M3 屈服门槛替代 c_min 冻结（约 2–3 天；新主因，先于 M2）

**现状**：c_min=0.05（默认 :208）→ wall 判据 `:923-929`（`cement_ever>0` 门控 + `cement<c_min`）→ `:653-655` `pref *= (1-wall)` 永久冻结，"饿死→冻结→永远饿死"。

**改法（经代码调研修正的实现路径）**：
1. **可重启的屈服判据替代 c_min**：用局部壁面剪应力判断"泥浆能否被动员"：
   - `immobile = τ_w_eff ≤ f_safety · τ_y_local`
   - **τ_w_eff 计算**（关键：不能直接反推 dp/ds，速度归一化含任意因子，C1 核实不可靠）：用同切片**流动中的宽边壁剪 × 间隙比外推**（并联槽流等压力梯度 G 假设：τ_w(φ) ≈ G·g(φ)/2，G=2τ_w,ref/b_ref）。**核查修订4/5**：参考元取"该深度 |w| 最大且 >0 的元"（row 0 自身冻结会致 G=0 全片永冻）；平行槽流 τ_w=G·g/2 与 μ 无关，μ 只进速度剖面，无需"μ 差异修正"；全切片无流动元时回退全冻结。
   - `γ̇ = 6|w|/b`（:590 已有，标准槽流壁面剪切率；b 为体积校正后的物理间隙，C1 已确认）；`μ_reg`（:599-600）、`τ_y` 场（`_compute_props`:476-482）均在 `_compute_velocity` 作用域内。
   - `f_safety` 安全裕度（建议 1.1–1.2，防 Roustaei P2 惯性警告：局部 Re~O(10²) 时纯壁剪判据会过度解冻；保留微量残泥下限）。
   - **每步重算、可解冻**：去掉 `cement_ever` 一次性触发；排量增大（快替 1.0/1.5 段）或 Δρ 增大时窄边可重启。
2. **水泥进入门槛**（McLean 第二条）：`cement_blocked = (τ_w_eff < τ_y_cement_local)`，或改用 B25 判据 `c̄min=[H/(a·τY,1)]^(1/(2n₂+1))`；两方向都保留，用 `enable_yield_gate` 下细分子判据可配。
3. **浮力项**（Lockyear Eq.2 含浮力部分）：**标注为模型扩展项**（Roustaei 两篇均无浮力项、P2 明言留给 3D，论文不引这两篇作浮力依据）；浮力通过"宽边壁剪外推"中间接体现（Δρ 影响 G），不作为独立判据项。

**数值带入（hu101，主会话前期计算）**：无浮力门槛 standoff≈0.72；含 Δρ=140 浮力≈0.24；YP 温度减半后 0.12–0.36——**必须用井下温压修正后流变**（联动 M6 附加项/后续数据层），否则门槛位置不可信；实现前先做流变温度敏感性分析。

**测试（会踩的现有断言，需重写，C1 已核）**：`tests/test_improved_d2dga_annulus.py::TestCminCriteria`（:558-667）——
- `test_wall_consistency_after_run`（:592-596，断言 wall=1⟺cement<c_min）→ 按新判据重写；
- `test_cmin_0_3_more_wall_than_cmin_0_05`（:624-646）→ 新判据下改/删；
- `test_constructor_has_cmin_default`（:561-564）→ c_min 参数保留但语义改为 fallback，断言更新；
- `test_wall_zeros_velocity_in_wall_cells`（:598-622）→ 若 wall 更新移入 `_compute_velocity` 则调用面变。

**验收**：hu101 三剖面 between_cent 底部段窄边呈"阶段相关"启动/静止（快替段尝试重启），而非全程冻结；窜槽不再被 A4 拉到 0.992 反物理。

---

## 7. M2 b³ 流态修正（约 2–3 天；次要但必做）

**现状**：pref ∝ (b/b̄)²/η_mix·I1（:619, :652），分流 ∝ b³（:659-660）；全井默认两侧层流；Re 只在 :602-604 诊断不参与。

**改法（Maleki 2017 闭包，2026-08-23 精读）**：
1. 在 :652（pref 构造后）与 :653-655（wall 零化）之间插入局部流态修正：
```python
# 固定点迭代 2-3 轮（M-F17 附录A）：
# Rep = 6ρ·w²/[κp·(γ̇N)^n]        （式58，Metzner-Reed；γ̇N=6|w|/(2Ĥ)... 用局部 w、b 与局部 n/κp/τy）
# He  = τY·(ρ^n (2Ĥ)^(2n)/κp²)^(1/(2−n))   （式59，Hedström）
# 摩擦因子 Hw↔Rep：层流式60-61（含塞流核 yY=He/Hw） / 湍流 Dodge-Metzner 式65（对 Hw 迭代） / 过渡区摩擦因子对数插值（连续化，禁硬切换）
# 屈服静止（式63）：|S|≤ra·τY/H ⟹ |∇Ψ|=0  —— 与 M3 判据天然衔接
# 用新阻力律重算 pref 归一化分流 → 迭代至收敛（防"宽边自证湍流"虚假正反馈）
```
2. **关键约束（C1 精读汇总）**：局部 Re 用 Metzner-Reed+He（不能挂牛顿 2100）；过渡区连续插值；必须固定点自洽迭代；保持 Hele-Shaw 椭圆耦合与共同压力梯度（**不退回独立 1D 切片**，Foroushan 批评正是此）；层流侧**不降指数**（幂律指数 2+1/n≈3.39>b³=3，"b³ 过陡"仅对宽边湍流成立）；把"湍流救窄边"作为不变量排除（M-F19：e=0.6 全湍流 η_N 最多 35%、等密度湍流 η_N=0）。
3. 开关注册 `enable_regime_split`；湍流阻力律与 M3 的屈服静止共用同一 `_compute_velocity` 扩展返回值。

**验收**：e=0.78 混合律分配比落在 36–530（自研回归：Walton&Bittleston 解析解 + Z22 Table3 + Foolad 定性基准）；at_cent（e≤0.4，宽边偏过渡）结果基本保留；不做下降指数。M-F19 量级靶：e=0.3→η_N≥95%、e=0.6→30–35%。

---

## 8. 附加项 A：I3 局部化（约 0.5 天）

**现状**：`eta2 = float(np.mean(mu))`（:880，混合场均值，语义错）；`delta_rho = (rho.mean()-mud)·1000`（:882，全场均值，非局部 Δρ(c̄)）。

**改法**：
- 从 `_compute_velocity` 扩展返回值透传 `eta2` 场（`_compute_props` :499 已产出，:673 未返回）。
- **⚠️ 必改（核查修订1）**：`d2dga_buoyancy_flux`（`d2dga_flux.py:136`）现用内建 `max(eta2, 1e-9)`，对数组会崩（已实测 `ValueError`）；必须先改为 `np.maximum(eta2, 1e-9)`，eta2 场才能广播透传。
- **扩元解包是 3 处**（核查修订2）：`_compute_velocity` 8 元组扩元同步改 `annulus_d2dga.py:805`、`:933`，以及漏列的 `tests/test_yield_deadzone.py:209`（精确解包，不改会炸）。
- Δρ 局部化：`Δρ(c̄) = c̄·ρ₂ + (1−c̄)·ρ₁ − ρ₁`（两层密度）。
- **不动公式分母**：I3 公式保持 Zhang22 式4.26 逐字符一致（`I3=c²(1−c)³[4mc+3(1−c)]/{2m[mc³+1−c³]}`，d2dga_flux.py:110-111）；只改 eta2/Δρ 取值，不引发论文侧"6/m 换算"争议扩大。
- 开关注册 `enable_local_i3`。

**验收**：R2 消融（I3 单独）效率变化不再 +2e-10（现近乎无效的证据）；修复后 R2 行重测，论文表 5/6 联动。

---

## 9. 附加项 B：hu101 实测 standoff 接回 loader（约 0.5 天）

**现状**：`hu101_loader.py:222-234` 为 model_assumption 剖面（0.38–0.48，注释明示无实测）；实测剖面硬编码在 `scripts/hu101_standoff_measured_vs_assumed.py:48-83`（扶正器间 0.22–0.78、扶正器处 0.60–0.88）。

**改法**：用 `dataclasses.replace(base_well, standoff_profile=实测剖面)` 把实测值接回 `hu101_loader.py` 的 `WellSpec.standoff_profile`（脚本已验证此法可行）；其余 7 井维持假设值 + 标注 `model_assumption`。**不改默认**（保持"假设值=默认运行"），实测剖面作为可选 profile 提供，避免未经确认就改默认输入导致旧基线不可复现。数据层伴随：standoff ±0.1 区间标注（A5；依据：长江大学理论 vs 成像测井误差 ~11%、埕海实测裸眼段仅 33.4%）。

---

## 10. 验收基准集（轻量路线 3）

每个单元内建验收步骤，不另建完整试验场：

| 单元 | 验收锚 |
|------|--------|
| M0 | 既有键逐位一致 + η_N 三 case 复现 + 四口径窗遍历 + "25 倍"论点不破坏 |
| M1 | κ 扫描收敛 + 网格/Δt 双收敛（mixing 由 κ 主导）+ 前沿长度米级 |
| M4 | hu101 三剖面几何真变 + 体积守恒校核 |
| M3 | between_cent 底部阶段相关启动/静止 + 窜槽不反物理 + 流变温度敏感性先行 |
| M2 | e=0.78 分配比 36–530 + at_cent 保留 + M-F19 量级靶 |
| I3 | R2 消融不再 +2e-10 |
| standoff | loader 可选 profile 生效 + ±0.1 区间 |

开发期回归脚本：`scripts/hu101_standoff_measured_vs_assumed.py`（三剖面）——判据：η_N 单调、between 不反转、mixing 收敛、M-F19 量级靶。

---

## 11. 重跑与论文数字影响

**重跑节奏（Q3-c）**：M0 一轮（零物理，验证即过）→ M1 一轮（独立归因）→ M4+M3+M2 合并一轮；全量 8 井。开发期全程 hu101 三剖面快检。

**论文数字对照表（M1 提交说明必附，C2 已核）**：表4（8 井效率/窜槽/混浆）、表5（消融）、表6（网格收敛）、67% 居中度阈值、呼101"65.04% vs 62.77% +2.27pp"、"20 秒"卖点（4 处）。**规避**：M0 与 M1 分开提交；改表时同步更新论文框架/撰写指南相应条目；"窄边拉开 25 倍"核心论点不得被 M1 破坏。

---

## 12. 风险清单与回退

| 风险 | 影响 | 缓解/回退 |
|------|------|----------|
| M2 自研参数化无文献闭合式 | 审稿可击穿 | 配 Walton&Bittleston 解析解 + Z22 Table3 + Foolad 三重回归；不改层流指数 |
| M3 浮力项无文献背书 | 论文表述被质疑 | 标注"模型扩展项"，不引 Roustaei 作浮力依据 |
| M3 门槛数值不可定标（0.12–0.75 摆动） | 门槛位置不可信 | 流变温压修正先行 + 敏感性分析；f_safety 裕度留 P2 惯性 |
| M1 移动全部基线 | 论文数字失效 | 与 M0 分开、提交附数字对照表、CFL on/off 存档 |
| e 开放到 0.90 重演 529:1 | 窜槽拉爆 | M4 排在 M3/M2 之后（先门槛后护栏） |
| I3 局部化改变 R2 消融 | 论文表5/6 连锁 | I3 单独开关、R2 重测、公式分母不动 |
| wall 判据移进 `_compute_velocity` 改变停泵行为 | 停泵期错冻结 | **核查修订3：原稿"q=0→G=0→自然冻结"不成立**——γ̇ 用 w_prev，停泵期 w_prev 不更新且非零，会误解冻。wall 更新须显式门控在 `pump_active`（:801）内，停泵期保持上一步 wall；并保留一步滞后防振荡（核查修订6） |
| 体积校正×M2 pref 变化 | 体积守恒疑云 | M2 只改 pref 权重，总流量由 area_weight 归一（:660），体积闭合由 overfilled(:850-855)+`_limit_phase_volume`(:918-921) 独立保证 |

**回退**：所有单元默认关=与 2832a08 完全一致；任一机制出险，关对应开关即可，无级联破坏。

---

## 13. 反方论点（Devil's Advocate）

1. 湍流不能"救"窄边（M-F19/Foolad 2021 一致）→ M2 预期是"窜槽量级温和修正"，η_E 大变化来自 M3/M5，勿把 M2 当提效手段。
2. A4 的 between 惩罚可能部分正确（底部 8.4mm 窄边泥浆静止有物理依据）→ 修复后若 between 效率大涨、窜槽大降，要检查是否过度矫正（残层 25–55% 与高体积效率并存是真实实验现象）。
3. 66.7% 是工程惯例非物理阈值（Jung&Frigaard 模拟 standoff 25% 仍有效）→ 敏感性扫描放宽到 e∈[0,1]、关注 e→1 质变，不向 66.7% 校准。
4. M2 无闭合公式 → 三重回归必须落地，否则不可辩护。
5. Roustaei P2"动了≠清走" → M0 指标应区分 flowing area 与 displaced area（本设计 M0 用卷加权 η_N 已部分覆盖；深化留后置）。

---

## 14. 提交与测试约定

- 每个单元独立 commit，信息格式 `feat(annulus): <Mx 描述>（开关=xx 默认关）`；每 commit 附单元验收结果。
- 测试：M0 新增 test_m0_metrics；M3 重写 TestCminCriteria 三断言；M2 新增局部 Re/闭包单测；M4 无既有断言冲突。`tests/paper_data/test_paper_result_schema.py` 为孤儿测试（import 不存在的模块），ignore，不新增风险。

---

## 附录 A：file:line 速查

| 位置 | 内容 |
|------|------|
| `annulus_d2dga.py:303` | e clip [0.05,0.55]（M4） |
| `:314-336` | h/b 构造 + 体积校正 scale（b=物理间隙，C1 已核） |
| `:589-604` | γ̇=6w/b、μ_reg、Re（诊断用）（M2/M3） |
| `:619,652-660` | pref ∝b³、wall 零化、面积归一（M2/M3） |
| `:673` | `_compute_velocity` 返回 8 元组（扩元供 M2/M3/I3） |
| `:860-863` | `_smooth_dispersion` 调用（硬编码系数）（M1） |
| `:880,882` | eta2=mean(mu)、delta_rho=全场均值（I3 局部化） |
| `:923-929` | c_min wall 判据（M3） |
| `:988-991` | 失稳指数饱和（M0） |
| `:1079-1096,1124-1130` | summary 构建、Tier0 try/except（M0） |
| `d2dga_flux.py:110-111,120-139` | I3 公式（不动）、buoyancy_flux（接数组） |
| `displacement_metrics.py:211-225` | `_narrow_quarter_efficiency` 定义（M0） |
| `well_spec.py:70-93` | EvaluationWindow（M0 窗遍历） |
| `hu101_loader.py:222-234` | 假设 standoff 剖面（附加B） |
| `scripts/hu101_standoff_measured_vs_assumed.py:48-83` | 实测剖面硬编码（附加B） |
| `tests/test_improved_d2dga_annulus.py:558-667` | TestCminCriteria（M3 重写） |

## 附录 B：关键文献锚点（精读确认 2026-08-23）

- **Maleki & Frigaard 2017**（J. Eng. Math. 107:201-230）：式58 局部幂律 Rep、式59 Hedström 数、式60-66 层流/屈服/湍流 Dodge-Metzner 闭式 + 广义指数、附录A.3 过渡区摩擦因子对数插值；"宽边湍流窄边层流甚至静止可同截面并存"。
- **Maleki & Frigaard 2019**：Couturier 四规则（①顶替液密度≥被顶替+10% ②摩阻梯度≥+20% ③窄边剪应力>被顶替屈服 ④宽边顶替速度≤窄边被顶替速度）；e=0.3→η_N≥95%、e=0.6→30–35%、体积比 3.25/6.15。
- **Foroushan et al. 2020**（SPE JDC 35(2):297-316）：常数表观黏度使窄边速度高估 43%（e=0.5）/68.6%（e=0.75）；批评逐切片+同心动量方程。
- **Roustaei et al. 2015 P1/P2**（JNNFM 220:87-98 / 226:1-15）：临界压差动员判据；无屈服→无静态残泥；P2 惯性警告（中 Re 段静区反而变大；"动了≠清走"）；**无浮力项**。
- **Escudier et al. 2002**：无 b³ 直接评价（间接背书双层流模型）；e=0.8 宽边几乎完美、中缝 -9%。
- **Jung & Frigaard 2022**（doi:10.1016/j.petrol.2021.109622）：67% 为工程惯例（源 Lee 1986）；模拟显示 standoff 25% 仍可有效顶替；模型=Maleki 2017 架构 + stiff-string 居中 + η_N/CBL 数字化验证。

## 附录 C：本文档未覆盖（后置/不在范围）

椭圆流函数横流（Z22 式2.2-2.3）、瞬态流函数+停泵有限时间冻结（MG07）、TVD minmod 输运、Uzawa 精确屈服（B25 附录A）、8 井 CBL 窗口对照闭环模块、流变井下温压修正——均属 Track B / 数据层后置，本设计不含。
