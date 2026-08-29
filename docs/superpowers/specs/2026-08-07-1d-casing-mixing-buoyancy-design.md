# 套管段1D模型：混浆段增强分散 + 浮力滑移物理化

> **日期**：2026-08-07
> **状态**：已实施（T1-10/T1-11 已落地并有单测；2026-08-22 弥散/混浆已接线生产 boundary_bridge，U型管 T1-12 已放弃并删除钩子）
> **任务编号**：T1-10（混浆段）+ T1-11（浮力滑移物理化）
> **叙事目标**：论文"入口边界精度提升" -- 从阶跃切换到物理化混浆+浮力修正
> **工作量估计**：T1-10 ~1.5天 + T1-11 ~1天 = ~2.5天

---

## 实施后记（2026-08-22）

- **实现位置**：`cemdisp/transport1d/casing_flow.py`（T1-10 混浆增强分散、T1-11 Atwood 浮力修正均已落地）。
- **测试**：`tests/test_casing_mixing_buoyancy.py`，覆盖浮力物理化、混浆增强、默认参数与接线集成。
- **生产接线**：`cemdisp/models2d/boundary_bridge.py` 已同步接通弥散/混浆路径（由另一 agent 改造）。
- **T1-12 U 型管**：决定放弃，`enable_utube` 参数、`_utube_corrected_arrival_time` 空钩子及 `_build_shoe_timeline` 中的调用已全部删除。
- **弥散量级修复（2026-08-22 接线当日）**：接线后首次端到端实测暴露 Taylor-Aris 渐近公式在固井 Pe~10^8 下失效（D_eff~10^5 m²/s、σ_t~10^4 s，hu101 效率崩塌至 3%）。已在 `_compute_dispersion_coefficient` 中以对流尺度混合上限 `dispersion_alpha·U·R` 截断（湍流弥散文献量级），并修复时间线自弥散 bug（prev_fluid 跳过同名流体）与负时间事件（子事件 clamp≥0、σ_t≤0.5·t_travel 防御上限）。修复后 hu101 消融：弥散贡献约 -1.0pp（0.5260→0.5161），量级回到物理合理区间；重力（Atwood）贡献约 -6.5pp。

---

## 1. 背景与动机

### 1.1 当前模型的问题

套管段1D体积追踪模型（`cemdisp/transport1d/casing_flow.py`）存在两个物理缺陷：

1. **重力修正基准错误**（`casing_flow.py:670`）：
   ```python
   gravity_factor = self.settling_velocity_factor * (fluid_density - 1000.0) / 1000.0
   ```
   以水密度 1000 kg/m³ 为基准，而非被顶替流体密度。对于密度 2100 kg/m³ 的水泥浆，修正幅度仅 0.165%，远小于实际浮力效应。

2. **混浆仅靠 Taylor-Aris 分散**：当前 `_apply_dispersion_to_timeline` 的过渡带宽度完全由分子扩散驱动的 Taylor-Aris 分散决定，未考虑密度不稳定（重驱轻）导致的界面混合增强。文献（Dai 2024 Petroleum Research）证明忽略管内混浆会导致顶替效率被显著高估。

### 1.2 文献支撑

| 文献 | 关键贡献 |
|------|---------|
| Dai et al. (2024) Petroleum Research | Atwood 数 `At=(ρh-ρl)/(ρh+ρl)` 表征密度差；3D VOF 模型证明同心模型高估效率 |
| Skadsem et al. (2022) J. Energy Res. Tech. | 1D 垂直井顶替模型，Froude 数判据 Fr≈0.94 |
| Zhang & Frigaard (2022) JFM 947 A32 | 浮力数是垂直井固井最关键无量纲参数 |
| Maleki & Frigaard (2018) Phys. Fluids | 屈服应力抑制浮力效应 |

### 1.3 已确认的决策

- 论文叙事：**入口边界精度提升**（论文发表优先）
- 实施范围：**T1-10 混浆段 + T1-11 浮力滑移物理化**
- U型管：**本次不实现**，只留 `enable_utube` 接口钩子
- 混浆方案：**方案A -- 基于失稳判据的增强分散**（复用 Dai 2024 物理判据）

---

## 2. 代码现状（已核实）

### 2.1 修改文件

仅修改 `cemdisp/transport1d/casing_flow.py`。不涉及 `boundary_bridge.py`、`shoe_timeline.py`、`pipe_exit_state.py`。

### 2.2 关键方法位置

| 方法 | 行号 | 现状 |
|------|------|------|
| `__init__` | L113-164 | 8个参数 |
| `_gravity_corrected_arrival_time` | L634-688 | **L670 bug**：以水为基准 |
| `_compute_dispersion_coefficient` | L277-344 | T1-8 已实现三分支 |
| `_apply_dispersion_to_timeline` | L346-416 | **L378** 混浆增强注入点 |
| `_build_shoe_timeline` | L418-469 | L440/445 重力修正调用点 |
| `_get_fluid_density` | 已存在 | 按名称查密度 |
| `_effective_pipe_radius_m` | 已存在 | 计算等效管半径 |

### 2.3 数据结构

`FluidSpec`（`cemdisp/data/fluid_spec.py`）已有字段：
- `name`, `role`, `density_kg_m3`, `rheology_model`
- `plastic_viscosity_pa_s`, `yield_stress_pa`, `power_law_n`, `consistency_k`

`FluidRole` 枚举：MUD, WASH, SPACER, LEAD, INTERMEDIATE, TAIL, DISPLACEMENT, FLUSHER, OTHER

---

## 3. T1-11 浮力滑移物理化

### 3.1 目标

将 `_gravity_corrected_arrival_time` 中的经验乘子替换为基于 Atwood 数的物理化浮力修正，修正密度差基准（从水改为被顶替流体）。

### 3.2 新增构造函数参数

```python
enable_buoyancy_physics: bool = True      # 启用物理化浮力修正（替代旧经验乘子）
buoyancy_correction_factor: float = 1.0   # 浮力修正标定系数（审查修正：0.1 时修正量与旧公式几乎相同，初值 1.0，需六井数据反标定）
```

当 `enable_buoyancy_physics=True` 时，使用新的 Atwood 数公式；当 `False` 时，回退到旧经验乘子（保持向后兼容）。

### 3.3 方法签名变更

`_gravity_corrected_arrival_time` 新增参数 `displaced_fluid_name: str`：

```python
def _gravity_corrected_arrival_time(
    self,
    arrival_time_s: float,
    fluid_name: str,
    displaced_fluid_name: str,   # 新增
    fluids: tuple[FluidSpec, ...],
    well_spec: WellSpec | None = None,
) -> float:
```

### 3.4 修正逻辑

```python
if not self.enable_buoyancy_physics:
    # 回退到旧经验乘子（向后兼容）
    # 将当前 L668-688 的旧逻辑提取为 _legacy_gravity_correction 私有方法
    return self._legacy_gravity_correction(arrival_time_s, fluid_name, fluids, well_spec)

rho_fluid = self._get_fluid_density(fluid_name, fluids)
rho_displaced = self._get_fluid_density(displaced_fluid_name, fluids)

# Atwood 数（无量纲密度差）
at = abs(rho_fluid - rho_displaced) / (rho_fluid + rho_displaced)

# 浮力修正因子
gravity_scale = self.g_constant / 9.81
gravity_factor = self.buoyancy_correction_factor * at * gravity_scale

# 井斜角投影（保留原有逻辑）
if well_spec is not None and well_spec.inclination_profile:
    avg_inclination_rad = self._average_inclination_rad(well_spec)
    gravity_factor *= max(math.cos(avg_inclination_rad), 0.0)

# 屈服应力抑制（保留原有逻辑）
fluid = next((f for f in fluids if f.name == fluid_name), None)
if fluid is not None and fluid.yield_stress_pa and fluid.yield_stress_pa > 0.0:
    pipe_radius_m = self._effective_pipe_radius_m(well_spec)
    delta_rho = abs(rho_fluid - rho_displaced)
    tau_critical = delta_rho * self.g_constant * pipe_radius_m
    if tau_critical > 1e-6:
        yield_ratio = min(fluid.yield_stress_pa / tau_critical, 1.0)
        gravity_factor *= (1.0 - 0.8 * yield_ratio)

# 方向判断：重驱轻 → 加速；轻驱重 → 减速
if rho_fluid > rho_displaced:
    return max(arrival_time_s * (1.0 - gravity_factor), 0.0)
else:
    return arrival_time_s * (1.0 + gravity_factor)
```

### 3.5 调用点更新

`_gravity_corrected_arrival_time` 有 3 处调用，均需传入被顶替流体名：

1. **`run()` L190**：迭代 `scheduled_steps` 时，被顶替流体为上一步的流体（或初始泥浆）
2. **`_build_shoe_timeline()` L440**：FRONT_ARRIVAL 修正
3. **`_build_shoe_timeline()` L445**：REAR_EXIT 修正

**循环改造说明**：`_build_shoe_timeline` L435 和 `run()` L187 当前循环为 `for scheduled in scheduled_steps:`（无索引），需改为 `for i, scheduled in enumerate(scheduled_steps):` 以获取 `current_index` 传入 `_displaced_fluid_name`。

新增辅助方法确定被顶替流体：

```python
@staticmethod
def _displaced_fluid_name(
    scheduled_steps: tuple[_ScheduledStep, ...],
    current_index: int,
    initial_fluid: str,
) -> str:
    """返回当前步骤的被顶替流体名（上一步流体，或初始泥浆）。"""
    if current_index <= 0:
        return initial_fluid
    return scheduled_steps[current_index - 1].step.fluid_name
```

---

## 4. T1-10 混浆段增强分散

### 4.1 目标

在流体界面处，根据密度不稳定判据（Atwood 数 + Reynolds 数）对 Taylor-Aris 分散系数施加增强因子，使过渡带宽度反映物理混浆效应。

### 4.2 新增构造函数参数

```python
enable_mixing_enhancement: bool = True       # 启用混浆增强分散
mixing_enhancement_factor: float = 5.0       # 增强系数 k_mix
max_mixing_enhancement: float = 10.0         # 增强因子上限
has_plug: bool = False                        # 是否有胶塞（有则抑制混浆）
```

### 4.3 新增方法：`_interface_instability_factor`

```python
def _interface_instability_factor(
    self,
    fluid_next: FluidSpec,
    fluid_prev: FluidSpec,
    pipe_radius_m: float,
    mean_velocity_m_s: float,
) -> float:
    """界面失稳增强因子（垂直井适配）。

    返回 >1 表示失稳增强分散，=1 表示稳定（仅 Taylor-Aris）。
    受 Dai 2024 启发的垂直井适配判据。注意：Dai 2024 原判据（eq. A.10/A.14/A.15）
    依赖 cos(β) 和 Fr，专为斜井设计；垂直井（β≈0）时 vt→0、Fr→∞，原判据不可直接用。
    此处简化为 At>0 AND Re>100 的密度不稳定判据。

    物理依据：
    - 重驱轻（At > 0）→ Rayleigh-Taylor 型界面失稳 → 混合增强
    - 黏度差越大 → 界面越不稳定
    - 有胶塞时 → 胶塞刮拭阻止混合 → 增强因子=1
    """
    if self.has_plug:
        return 1.0

    rho_h = max(fluid_next.density_kg_m3, fluid_prev.density_kg_m3)
    rho_l = min(fluid_next.density_kg_m3, fluid_prev.density_kg_m3)
    at = (rho_h - rho_l) / (rho_h + rho_l)  # Atwood number

    if at < 1e-9:
        return 1.0  # 等密度，无失稳

    # 有效黏度（几何平均，Dai 2024 eq. A.12）
    mu_next = self._effective_viscosity(fluid_next, mean_velocity_m_s, pipe_radius_m)
    mu_prev = self._effective_viscosity(fluid_prev, mean_velocity_m_s, pipe_radius_m)
    mu_mean = math.sqrt(max(mu_next * mu_prev, 1e-12))

    # Reynolds number
    rho_avg = (fluid_next.density_kg_m3 + fluid_prev.density_kg_m3) / 2.0
    D = 2.0 * pipe_radius_m
    Re = rho_avg * mean_velocity_m_s * D / max(mu_mean, 1e-12)

    # 垂直井密度不稳定判据：重驱轻 + Re 足够大 → 失稳
    if at > 0 and Re > 100.0:
        enhancement = 1.0 + self.mixing_enhancement_factor * at * math.sqrt(Re / 100.0)
        return min(enhancement, self.max_mixing_enhancement)
    return 1.0
```

> **注**：增强公式 `1 + k_mix·At·√(Re/100)` 是**经验公式**，不来自特定论文，系数 k_mix=5.0 需六井数据标定。HU102 水泥-泥浆界面验证：At=0.0194, Re≈1950 → enhancement≈1.43，过渡带宽度增加 ~20%，物理上合理但偏保守。

### 4.4 新增辅助方法：`_effective_viscosity`

```python
def _effective_viscosity(
    self,
    fluid: FluidSpec,
    mean_velocity_m_s: float,
    pipe_radius_m: float,
) -> float:
    """计算流体在给定剪切率下的表观黏度。

    使用 Dai 2024 eq. A.11: μ = τ₀/(u/D) + k·(u/D)^(n-1)
    """
    if mean_velocity_m_s < 1e-9 or pipe_radius_m < 1e-9:
        return fluid.plastic_viscosity_pa_s or 0.01

    shear_rate = 8.0 * mean_velocity_m_s / (2.0 * pipe_radius_m)  # 壁面剪切率近似
    shear_rate = max(shear_rate, 1e-6)

    if fluid.rheology_model == RheologyModel.NEWTONIAN:
        return fluid.plastic_viscosity_pa_s or 0.01

    tau_y = fluid.yield_stress_pa or 0.0
    k_cons = fluid.consistency_k or 0.01
    n = fluid.power_law_n or 1.0

    return tau_y / shear_rate + k_cons * shear_rate ** (n - 1.0)
```

### 4.5 在 `_apply_dispersion_to_timeline` 中注入

在 L378（`D_eff = self._compute_dispersion_coefficient(...)` 之后）：

```python
D_eff = self._compute_dispersion_coefficient(pipe_radius_m, fluid, U)

# T1-10: 混浆增强分散
if self.enable_mixing_enhancement:
    prev_fluid_spec = fluid_by_name.get(prev_fluid)
    if prev_fluid_spec is not None:
        instability = self._interface_instability_factor(
            fluid, prev_fluid_spec, pipe_radius_m, U
        )
        D_eff *= instability
```

---

## 5. U型管预留接口

### 5.1 新增构造函数参数

```python
enable_utube: bool = False  # U型管修正开关（默认关闭）
```

### 5.2 新增钩子方法

```python
def _utube_corrected_arrival_time(
    self,
    arrival_time_s: float,
    fluid_name: str,
    fluids: tuple[FluidSpec, ...],
    well_spec: WellSpec | None = None,
) -> float:
    """U型管/自由下落修正钩子。当前为空实现，预留未来扩展。

    未来实现方向（T1-12）：
    1. 计算套管-环空静压平衡
    2. 检测自由下落条件（ΔP_hydro > ΔP_friction）
    3. 修正到达时间（加速效应）
    """
    if not self.enable_utube:
        return arrival_time_s
    # 未来：填入 U 型管物理计算
    return arrival_time_s
```

### 5.3 调用点

在 `_build_shoe_timeline` 中，重力修正之后调用：

```python
if self.enable_gravity:
    front_time_s = self._gravity_corrected_arrival_time(...)
front_time_s = self._utube_corrected_arrival_time(front_time_s, ...)  # 新增
```

---

## 6. 参数汇总

| 参数 | 默认值 | 所属任务 | 说明 |
|------|--------|---------|------|
| `enable_buoyancy_physics` | `True` | T1-11 | 启用 Atwood 数浮力修正 |
| `buoyancy_correction_factor` | `1.0` | T1-11 | 浮力修正标定系数（需六井标定） |
| `enable_mixing_enhancement` | `True` | T1-10 | 启用混浆增强分散 |
| `mixing_enhancement_factor` | `5.0` | T1-10 | 增强系数 k_mix |
| `max_mixing_enhancement` | `10.0` | T1-10 | 增强因子上限 |
| `has_plug` | `False` | T1-10 | 胶塞开关 |
| `enable_utube` | `False` | 预留 | U型管开关 |

---

## 7. 测试计划

### 7.1 单元测试

| 测试项 | 预期 |
|--------|------|
| T1-11：Δρ=0（等密度流体） | 修正为零，到达时间不变 |
| T1-11：重驱轻（水泥驱泥浆） | 到达时间缩短 |
| T1-11：轻驱重（泥浆驱水泥） | 到达时间延长 |
| T1-11：τ_y→0 极限 | 回到无屈服抑制的浮力修正 |
| T1-11：`enable_buoyancy_physics=False` | 回退到旧经验乘子 |
| T1-10：At=0（等密度界面） | 增强因子=1 |
| T1-10：重驱轻+高Re | 增强因子>1 |
| T1-10：has_plug=True | 增强因子=1 |
| T1-10：enable_mixing_enhancement=False | 不增强 |
| U型管：enable_utube=False | 不影响结果 |

### 7.2 回归测试

- 六井全量回归：HU102, HU101, HU103, HU2, HT1-001, HT1-003, HT1-004
- 验证顶替效率变化在合理范围（预期 ±5% 以内）
- 验证 `PipeExitState.phase_fractions` 多相分数正确传递到环空入口

---

## 8. 实施顺序

```
Step 1: T1-11 浮力滑移物理化（~1天）
  ├── 新增 _displaced_fluid_name 辅助方法
  ├── 修改 _gravity_corrected_arrival_time 签名和逻辑
  ├── 更新 3 处调用点
  └── 单元测试

Step 2: T1-10 混浆段增强分散（~1.5天）
  ├── 新增 _effective_viscosity 辅助方法
  ├── 新增 _interface_instability_factor 方法
  ├── 修改 _apply_dispersion_to_timeline 注入增强因子
  └── 单元测试

Step 3: U型管预留接口（~0.5天）
  ├── 新增 _utube_corrected_arrival_time 钩子
  └── 在 _build_shoe_timeline 中调用

Step 4: 六井回归测试（~0.5天）
```

---

## 9. 风险与缓解

| 风险 | 概率 | 影响 | 缓解 |
|------|------|------|------|
| buoyancy_correction_factor 标定不准 | 中 | 中 | 用六井数据反标定，初值 1.0（审查修正） |
| mixing_enhancement_factor 过大导致过渡带不物理 | 中 | 低 | 上限保护 max_mixing_enhancement=10 |
| 被顶替流体判断错误（首步无前置流体） | 低 | 高 | 默认回退到初始泥浆 |
| 胶塞信息不在 WellSpec 中 | 低 | 低 | has_plug 默认 False，手动设置 |
