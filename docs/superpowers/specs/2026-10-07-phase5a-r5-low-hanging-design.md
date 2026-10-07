# Phase 5a 设计规格：R5 低垂五项（2026-10-07）

> 上游 = 《详细执行计划_4d主批起_续作_2026-10-07.md》§4 Phase 5「5a R5 低垂」+
> 《温压耦合真响应改进方案_2026-10-05》§「R5 配套修复（各 ≤1 天，可穿插并行）」。
> 起点 HEAD `ccdb893`（Phase 4d 主批已收口）；测试基线 17F/1218P/241sub。
> 依据纪律（§0.2）：每项须挂锚（文献 / 仓内既有实现 / 现场资料）；**不发明物理**。

---

## 5a-① P1-1 停泵衰减判据后处理（新诊断模块）

**判据（引用）**：Moyers-González & Frigaard —— 停泵后 Q=0，界面冻结当且仅当
屈服抗力压得住浮力体力散度：`τ_{Y,min}/(1+e) ≥ ‖div f‖∞ / 2`。与 R2（屈服门进流场）正交：
R2 讲**流动中**的通道导纳，本判据讲**静止后**的冻结判据，二者共用同一 τy 场但作用时段不同。

**实现面（新增 `cemdisp/diagnostics/stop_pump_freeze.py`，纯后处理，零求解器改动）**：

| 量 | 取法（本模型口径） |
|---|---|
| `τ_{Y,min}(z)` | 混合屈服应力场逐深最小值：`τy(z) = Σ_phase c_phase(z)·τy_phase`（与 `_compute_props` 的 τy 场同一混合规则）；相 τy 取温变派生态（T-on）或静态 spec（T-off） |
| `e` | 偏心度 = `1 − standoff`（逐深剖面取界面带内最大值） |
| `f` | 浮力体力向量 = `_buoyancy_force_vector` 的式 (2.5b) 口径 `f = (r_a·cosβ, r_a·sin(πφ)·sinβ)/F²`；**本模块按同一式重建** |
| `‖div f‖∞` | 沿**间隙坐标 x** 对方位分量求散度：`∂f_x/∂x`，逐深取 ∞-范数 |

**尺寸口径（本次映射假设，**待追认**——计划原文写的是无量纲判据的简写）**：
`f` 在式 (2.5b) 中为无量纲体力（被 `F²` 归一），而 `τ_Y` 为 Pa。为使两侧同量纲，
判别式取**带间隙因子**的形式
`τ_{Y,min}/(1+e) ≥ (b(z)/2) · ‖∂f_x/∂x‖∞ / 2`
——即浮力体力散度乘以特征长度 `b/2` 才还原为应力尺度的 Pa。
依据：槽流中体力梯度作用在半个间隙上的等效壁面剪应力 = `(∇·f)·(b/2)`（与
`_yield_gate_wall` 的 `τw = G·b/2` 外推同型）。**该因子是本次实现的显式映射假设**，
与计划原文简写不一致处**以本行为准并标注待追认**。

**输出**：逐深判据值 / 冻结布尔场 / 界面带冻结占比 / 判定（冻结/未冻结/部分），
落 `results/_probe_5a_低垂_2026-10-07/`。
**与既有 wall 通道的关系**：本模块不改 `_yield_gate_wall`，只作停泵后的**后处理判据**，
结果不进 summary（关1 红线）。

## 5a-② I₃ 窄口回落 hb_fix（T-on 未覆盖水泥相丢回落）

**病灶**（方案 §4「存疑需修复」）：`_phase_cement_tau_y` 在 T-on 下**无条件**返回
`_fluid_yield_stress(fluid)`，绕过了 `_cement_phase_yield_stress` 的 `hb_fix_cement_tau_y`
常数机制 ⇒ **公式未覆盖的水泥相**（密度越出 `[cm_lo, cm_hi]` 捕获区）在 T-on 下丢失常数回落。

**修复**：T-on 分支加覆盖判定——
`if self.enable_temperature_rheology and is_replaced(base, self.rheology_formula_params): 公式 τy;
 else: self._cement_phase_yield_stress(base)`（常数机制/静态 spec）。

**组合语义裁定**：**覆盖优先、未覆盖回落**——"公式绝对替换"只作用于被覆盖相；
未覆盖相按 T-off 同一路径取常数机制（`hb_fix_cement_tau_y=True`）或静态 spec（`False`）。

**关2 依据（逐位）**：`hb_fix_cement_tau_y=False`（生产默认）时
`_cement_phase_yield_stress == _fluid_yield_stress` ⇒ **两支同值，默认路径逐位不变**；
仅 `hb_fix_cement_tau_y=True` 时未覆盖相由「丢静态值/0」改为「常数回落」（正是修复目标）。

**新增单一真源**：`rheology_vs_temperature.is_replaced(fluid, params=None) -> bool`
（与 `fluid_at` 的替换判据同源：`_route` 族 + 水泥密度档 `cm_lo ≤ d ≤ cm_hi`）。

## 5a-③ 1D 宾汉丢 PV 地雷

**病灶**：`casing_flow._effective_viscosity` 对 Bingham 用
`k_cons = fluid.consistency_k or 0.01; n = 1.0` ⇒ 返回 `τy/γ̇ + 0.01`，
**PV（`plastic_viscosity_pa_s`）被 0.01 顶掉**（Dai 2024 A.11 的 `k` 对宾汉应为 PV、n=1）。

**修复（实施中改为 opt-in，见下）**：Bingham 分支显式用 PV：`μ_eff = τy/γ̇ + PV`
（PV 为 None 时才回退 0.01，沿用原缺值回退语义）；HB / PowerLaw 分支不动。

**⚠ 实施更正（2026-10-07，实测驱动）**：原判断「生产不可达 ⇒ 零影响（逐位）」
**不完整**——该入口在 `has_plug=False` 的配置下**可达**（经
`_interface_instability_factor` → 弥散带增强）。直接改值会位移
`tests/history/test_casing_mixing_contact_time.py::TestMixingContactTime::test_default_off_bitwise`
的 **hu103 默认路径冻结锚**（实测 15 个事件失配 + 1 个 FAILED ⇒ 新增失败 16 处），
即触碰「缺省/不传 ⇒ 逐位 = HEAD」红线（计划 §8-16 同型）。
⇒ 按项目铁律**一切 opt-in**：新增 ctor 开关 `CasingFlowSolver(bingham_pv_fix=...)`，
**默认 False = 沿用 τy/γ̇ + 0.01（逐位=HEAD）**；置真才用 PV。
生产 8 井 `has_plug=True` ⇒ 该入口不可达，两种取值对生产数字均无影响。
**新增改动面**：`casing_flow.CasingFlowSolver.__init__` 加 keyword 形参
`bingham_pv_fix: bool = False`（超出原 §2 清单，故在此补记）。

## 5a-④ smooth_break / mud_extrapolate 对照档

- `mud_extrapolate=True`（钻井液 T>80 °C 按公式外推而非 clamp）对照档一跑。
- **`smooth_break` 不可跑**：`rheology_vs_temperature.fluid_at` 对其显式
  `raise NotImplementedError`（占位开关，忠实两段式是既定口径——红线「100 °C 断点不平滑」）。
  ⇒ 本项**如实声明为不可执行**并记录原因，不伪造对照（**范围缩减，按实报告**）。

## 5a-⑤ 缝隙律单因子 A/B

`enable_power_law_gap_law` on/off 单因子对照（旧路径只作对照、**不进生产**），
补「机制 A 净贡献」不可证伪缺口（方案 §「风险」第 8 条：论文不得把速度场
「只吃标量表观黏度」表述为「非牛顿顶替模拟」，须给 on/off 单因子数字）。

---

## 验收（计划 §5 五关）

| 关 | 判据 |
|---|---|
| 关1 | T-off 逐位：③ 不可达路径 + ② 默认 `hb_fix=False` 两支同值 ⇒ 冻结锚 JSON 不动；全量 17F/1218P/241sub 不增，FAILED 逐条一致 |
| 关2 | ②/③ 默认档逐位不变（新测试钉住：覆盖相走公式、未覆盖相回落、默认两支同值；Bingham μ_eff = τy/γ̇+PV） |
| 关3 | 方向性：② 在 `hb_fix=True` + T-on + 未覆盖相下 τy 由静态值 → 常数（方向=抬升）；⑤ 缝隙律 on/off 单因子差 |
| 关4 | ① 判据对停泵末态给出逐井裁决 + 间隙因子映射假设显式声明；⑤ 单因子数字落表 |
| 关5 | 对照档 cfl_clip=0、健康度与 4d 同型 |

## 回退

②③ 为分支内小改（`git revert` 可回）；① 为**新增独立模块**（零耦合）；
④⑤ 只产 `results/` 日期目录，不覆盖既有产物。
