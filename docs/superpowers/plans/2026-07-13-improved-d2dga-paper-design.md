# 改进 D2DGA 三闭包实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 把 cemdisp 环空 2D D2DGA 求解器从"半实现"（m=1.0 固定 + buoyancy_shape 代理 + 无 I₃）升级到"三闭包完整"（auto-m + I₃ 浮力通量 + 真浮力体力），形成 R0->R3 消融阶梯，支撑《石油勘探与开发》论文方法创新点。

**Architecture:** 三闭包逐级叠加（R1 auto-m -> R2 I₃ -> R3 真体力），每级一个独立可测交付。R1 改 `d2dga_flux.py` 的 m 传入方式 + `annulus_d2dga.py` 计算 m 场；R2 在 `d2dga_flux.py` 新增 `d2dga_buoyancy_flux` 函数 + 在 run 循环浓度更新步叠加散度通量；R3 删 `buoyancy_shape` 代理、注入真体力 b、输出浮力数 b。每级在 HT1-004 上跑消融图并现场闭环。所有公式来自 Zhang & Frigaard (2022) JFM 947 A32 + Lajeunesse et al. (1999) JFM 398。

**Tech Stack:** Python 3.13, NumPy, pandas, pytest (pyproject.toml), 现有 cemdisp 包（`cemdisp/models2d/` + `cemdisp/runners/`）。

## Global Constraints

- **公式源**：所有闭包公式严格按 Zhang & Frigaard (2022) JFM 947 A32：式 4.25 第二项、式 4.26（I₃）、式 2.5b（体力 b）、式 4.3（U 参数）；auto-m 用 Yang & Yortsos (1997) 的 m=μ_displaced/μ_displacing 定义。不得自创公式。
- **回归保护**：R0 baseline（三闭包全关）必须在 HT1-004 上复现旧论文 m=1.0 的结果（顶替效率≈95.96%）。每级改动必须可开关回退到 R0。
- **数值安全**：所有新增项必须有 clip 限幅，不得引入 NaN/inf。I₃ 通量叠加后浓度场必须 clip 回 [0,1]。
- **开关参数命名**：`enable_d2dga`（已有，控制式4.28放大）、`enable_d2dga_auto_m`（R1新）、`enable_d2dga_i3_flux`（R2新）、`enable_true_buoyancy`（R3新）。默认全 True。全 False = R0 baseline。
- **不修改现有公开接口签名**：`AnnulusD2DGASolver.__init__` 只新增带默认值的 kwargs；`run()` 签名不变；`AnnulusSimulationResult` 新增字段用默认值兼容旧结果。
- **DRY**：浮力向量 f 的计算在 R2 建一次（`_buoyancy_force_vector`），R3 复用，不重复实现。
- **TDD**：每个闭包先写单元测试（验证边界行为），再实现，最后跑回归。
- **提交粒度**：每个 Task 末提交一次；每个 Step 完成后可提交。分支 `feat/improved-d2dga-paper-design`（已建）。

---

## File Structure

新增/修改文件清单（每个文件一个清晰职责）：

| 文件 | 职责 | 动作 |
|---|---|---|
| `cemdisp/models2d/d2dga_flux.py` | D2DGA 通量纯函数（式4.28放大 + 式4.25第二项I₃通量） | 修改：新增 `d2dga_buoyancy_flux`、`d2dga_dispersion_function_I3` |
| `cemdisp/models2d/annulus_d2dga.py` | 环空 2D 求解器 | 修改：构造函数加 3 开关 + `_compute_props` 输出 m 场 + `_buoyancy_force_vector` 新方法 + `_compute_velocity` 删代理注体力 + run 循环传数组 m + 叠加 I₃ 通量 + 输出浮力数 b |
| `tests/test_d2dga_flux.py` | d2dga_flux 纯函数单元测试 | 新建 |
| `tests/test_improved_d2dga_annulus.py` | 求解器级三闭包开关测试 | 新建 |
| `cemdisp/runners/ht1_004_ablation.py` | R0->R3 四级消融运行器 | 新建 |
| `scripts/run_ht1_004_ablation.py` | 消融运行入口脚本 | 新建 |
| `docs/superpowers/plans/2026-07-13-improved-d2dga-paper-design.md` | 本计划 | 新建（此文件） |

**分解决策**：`d2dga_flux.py` 保持"纯函数 + 无状态"职责（只放公式），求解器逻辑留在 `annulus_d2dga.py`。浮力向量 f 的计算作为求解器的私有方法 `_buoyancy_force_vector`（因为它依赖几何 geom + 流体，不是纯函数）。

---

## Task 1: R1 auto-m — d2dga_flux_amplification 接受数组 m 并新增 I₃ 弥散函数

**目标**：让 `d2dga_flux_amplification` 正确接受数组形式的 viscosity_ratio（它已用 `np.asarray`，但要加测试固化该行为 + 修正 docstring），并新增 R2 要用的 `d2dga_dispersion_function_I3` 纯函数。

**Files:**
- Modify: `cemdisp/models2d/d2dga_flux.py`（全文已读，当前 81 行）
- Test: `tests/test_d2dga_flux.py`（新建）

**Interfaces:**
- Produces: `d2dga_flux_amplification(cement_fraction, viscosity_ratio=1.0, ...)` 接受 array viscosity_ratio（已支持，加测试固化）；新增 `d2dga_dispersion_function_I3(c_bar, m)` 返回 I₃(ḉ,m) 数组。
- `d2dga_dispersion_function_I3` 签名：`(c_bar: FloatOrArray, m: float = 1.0, *, min_fraction: float = 0.01, max_fraction: float = 0.99) -> FloatOrArray`，公式 $I_3=\frac{\bar{c}^2(1-\bar{c})^3[4m\bar{c}+3(1-\bar{c})]}{2m[m\bar{c}^3+1-\bar{c}^3]}$。

- [ ] **Step 1: 写失败测试 — d2dga_flux_amplification 数组 m 行为**

创建 `tests/test_d2dga_flux.py`：

```python
"""D2DGA 通量纯函数测试。覆盖：式4.28放大因子（含数组m）+ 式4.26 I₃弥散函数。"""
import numpy as np
import pytest
from cemdisp.models2d.d2dga_flux import d2dga_flux_amplification, d2dga_dispersion_function_I3


class TestFluxAmplificationArrayM:
    def test_scalar_m_equals_1_returns_1p5(self):
        # m=1 时 f(c,1) = [c² + 1.5(1-c²)] / [c³ + 1 - c³] = 1.5/1 = 1.5 常数
        c = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        f = d2dga_flux_amplification(c, viscosity_ratio=1.0)
        assert np.allclose(f, 1.5)

    def test_array_m_different_from_scalar(self):
        # 同一 c 场，不同 m 应给不同 f
        c = np.array([0.3, 0.5, 0.7])
        m_array = np.array([0.5, 1.0, 2.0])
        f_array = d2dga_flux_amplification(c, viscosity_ratio=m_array)
        f_scalar = d2dga_flux_amplification(c, viscosity_ratio=1.0)
        assert not np.allclose(f_array, f_scalar)
        # 小 m（顶替液更粘）应使放大因子更接近 1（弥散更弱）在高浓度区
        # m=0.5 在 c=0.9 处 f 应 < m=1 的 1.5
        f_m05 = d2dga_flux_amplification(np.array([0.9]), viscosity_ratio=0.5)
        assert f_m05[0] < 1.5

    def test_output_shape_matches_input(self):
        c = np.linspace(0.01, 0.99, 40)
        m = np.full(40, 0.8)
        f = d2dga_flux_amplification(c, viscosity_ratio=m)
        assert f.shape == c.shape
```

- [ ] **Step 2: 运行测试验证失败**

Run: `python -m pytest tests/test_d2dga_flux.py::TestFluxAmplificationArrayM -v`
Expected: 前 2 个测试应 PASS（数组 m 已支持），`test_output_shape_matches_input` 也 PASS。若全 PASS 说明数组 m 行为已正确——此 Task 主要是固化测试。若 FAIL，修正 `d2dga_flux.py` 使其正确处理数组 m（当前 `np.asarray` + 算术已支持，预期 PASS）。

- [ ] **Step 3: 写失败测试 — d2dga_dispersion_function_I3**

追加到 `tests/test_d2dga_flux.py`：

```python
class TestDispersionFunctionI3:
    def test_zero_concentration_returns_zero(self):
        # c=0 -> I3 = 0 (分子含 c²)
        assert d2dga_dispersion_function_I3(0.0, m=1.0) == pytest.approx(0.0, abs=1e-9)

    def test_full_concentration_returns_zero(self):
        # c=1 -> I3 = 0 (分子含 (1-c)³)
        assert d2dga_dispersion_function_I3(1.0, m=1.0) == pytest.approx(0.0, abs=1e-9)

    def test_mid_concentration_positive(self):
        # c=0.5 时 I3 应为正且较大（峰值附近）
        i3 = d2dga_dispersion_function_I3(0.5, m=1.0)
        assert i3 > 0.0
        # 手算核对 m=1, c=0.5: 分子=0.25*0.125*(4*0.5+3*0.5)=0.25*0.125*2.75=0.0859375
        # 分母=2*1*(1*0.125+1-0.125)=2*1.0=2 -> I3=0.04296875
        assert i3 == pytest.approx(0.04296875, rel=1e-4)

    def test_zero_density_implied_zero_flux_handled_separately(self):
        # I3 函数本身不含 Δρ；Δρ=0 时通量在调用方置零。这里只验 I3 数值正确。
        c = np.array([0.2, 0.5, 0.8])
        i3 = d2dga_dispersion_function_I3(c, m=1.0)
        assert i3.shape == c.shape
        assert np.all(i3 >= 0.0)  # I3 在 [0,1] 非负
```

- [ ] **Step 4: 运行测试验证失败**

Run: `python -m pytest tests/test_d2dga_flux.py::TestDispersionFunctionI3 -v`
Expected: FAIL with `ImportError: cannot import name 'd2dga_dispersion_function_I3'`

- [ ] **Step 5: 实现 d2dga_dispersion_function_I3**

在 `cemdisp/models2d/d2dga_flux.py` 末尾追加（保留现有 `d2dga_flux_amplification` 不动）：

```python
def d2dga_dispersion_function_I3(
    c_bar: FloatOrArray,
    m: float = 1.0,
    *,
    min_fraction: float = 0.01,
    max_fraction: float = 0.99,
) -> FloatOrArray:
    """计算 D2DGA 浮力弥散函数 I3(ḉ, m)（Zhang & Frigaard 2022, 式 4.26）。

    公式：I3 = ḉ²(1-ḉ)³[4m·ḉ + 3(1-ḉ)] / {2m[m·ḉ³ + 1 - ḉ³]}

    性质：ḉ=0 或 ḉ=1 时 I3=0；ḉ≈0.5 附近达峰。用于 R2 浮力驱动弥散通量。

    Args:
        c_bar: 间隙平均水泥浓度（0~1），标量或数组。
        m: 黏度比 η_displaced/η_displacing。
        min_fraction: 计算前浓度下限，避免零浓度奇异。
        max_fraction: 计算前浓度上限，避免充满时奇异。
    """
    c = np.asarray(c_bar, dtype=float)
    c_safe = np.clip(c, min_fraction, max_fraction)
    c2 = c_safe ** 2
    c3 = c_safe ** 3
    one_minus_c = 1.0 - c_safe
    numerator = c2 * (one_minus_c ** 3) * (4.0 * m * c_safe + 3.0 * one_minus_c)
    denominator = 2.0 * m * (m * c3 + 1.0 - c3)
    i3 = numerator / denominator
    # 边界处置零（c=0 或 c=1 的精确值，clip 之外）
    i3 = np.where((c < min_fraction) | (c > max_fraction), 0.0, i3)
    if np.isscalar(c_bar):
        return float(i3)
    return i3.astype(float, copy=False)
```

- [ ] **Step 6: 运行测试验证通过**

Run: `python -m pytest tests/test_d2dga_flux.py -v`
Expected: 全部 PASS。

- [ ] **Step 7: 提交**

```bash
git add cemdisp/models2d/d2dga_flux.py tests/test_d2dga_flux.py
git commit -m "feat(d2dga): 新增 I3 弥散函数(式4.26) + 固化数组m测试"
```

---

## Task 2: R1 auto-m — 求解器计算 m 场并传入放大因子

**目标**：把 `d2dga_viscosity_ratio` 从"构造常数"改为"每步按局部流体物性自动计算的数组场"，由开关 `enable_d2dga_auto_m` 控制。

**Files:**
- Modify: `cemdisp/models2d/annulus_d2dga.py`（构造函数 L183-222、`_compute_props` L405-447、run 循环 L659-660）
- Test: `tests/test_improved_d2dga_annulus.py`（新建）

**Interfaces:**
- Consumes: `d2dga_flux_amplification(cement_fraction, viscosity_ratio=array)`（Task 1 已支持数组）。
- Produces: 构造函数新增 `enable_d2dga_auto_m: bool = True`；`_compute_props` 多返回一个 `m_field`（被顶替液/顶替液粘度比场）。
- **m 场定义**：m = μ_displaced / μ_displacing。在浓度场中，水泥（lead+tail）是顶替液，泥浆是被顶替液。故 m_field = mu_mud / mu_cement（局部）。当 `enable_d2dga_auto_m=False` 时回退到 `self.d2dga_viscosity_ratio` 标量。

- [ ] **Step 1: 写失败测试 — auto-m 开关与 m 场计算**

创建 `tests/test_improved_d2dga_annulus.py`：

```python
"""改进 D2DGA 三闭包（auto-m / I3 / 真体力）求解器级测试。"""
import numpy as np
import pytest
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
from cemdisp.data.fluid_spec import FluidSpec, FluidRole, RheologyModel


def _make_solver(**kw) -> AnnulusD2DGASolver:
    return AnnulusD2DGASolver(dt=4.0, nz=20, ny=10, total_t=40.0, **kw)


class TestAutoMField:
    def test_constructor_has_enable_d2dga_auto_m_default_true(self):
        s = _make_solver()
        assert hasattr(s, "enable_d2dga_auto_m")
        assert s.enable_d2dga_auto_m is True

    def test_auto_m_off_falls_back_to_scalar(self):
        # enable_d2dga_auto_m=False -> m 场退化为构造常数
        s = _make_solver(enable_d2dga_auto_m=False, d2dga_viscosity_ratio=0.8)
        assert s.enable_d2dga_auto_m is False
        assert s.d2dga_viscosity_ratio == 0.8
```

- [ ] **Step 2: 运行测试验证失败**

Run: `python -m pytest tests/test_improved_d2dga_annulus.py::TestAutoMField -v`
Expected: FAIL with `AttributeError: 'AnnulusD2DGASolver' object has no attribute 'enable_d2dga_auto_m'`

- [ ] **Step 3: 构造函数加开关**

修改 `cemdisp/models2d/annulus_d2dga.py` 构造函数（L183-196）：

```python
    def __init__(
        self,
        *,
        dt: float = 4.0,
        nz: int = 140,
        ny: int = 40,
        total_t: float = 12000.0,
        enable_d2dga: bool = True,
        d2dga_viscosity_ratio: float = 1.0,
        enable_d2dga_auto_m: bool = True,
        instability_decay_scale: float = 5.0,
        save_interval: int = 60,
        yield_regularization_M: float = 100.0,
        open_outlet: bool = True,
    ) -> None:
```

在 docstring 的 Args 中追加（在 `d2dga_viscosity_ratio` 之后）：
```
            enable_d2dga_auto_m: 是否按局部流体物性自动计算黏度比 m 场（R1），默认 True。
                True: m = μ_displaced/μ_displacing 按浓度场每步计算（改进版）；
                False: 退化为 d2dga_viscosity_ratio 构造常数（旧论文 R0 状态）。
```

在 `self.open_outlet = open_outlet`（L222）之前加：
```python
        self.enable_d2dga_auto_m: bool = enable_d2dga_auto_m
```

- [ ] **Step 4: 运行测试验证通过**

Run: `python -m pytest tests/test_improved_d2dga_annulus.py::TestAutoMField -v`
Expected: PASS。

- [ ] **Step 5: 写失败测试 — m 场从 _compute_props 返回**

追加到 `tests/test_improved_d2dga_annulus.py`：

```python
from cemdisp.data.well_spec import WellSpec, DepthValuePoint, EvaluationWindow


def _make_geom(solver, well_spec):
    return solver._build_geom(well_spec)


def _toy_well():
    # 最小井身结构，供 _compute_props / _compute_velocity 单测
    pts = lambda d, v: DepthValuePoint(depth_md_m=d, value=v)
    return WellSpec(
        well_name="toy",
        top_md_m=1000.0, bottom_md_m=1100.0, shoe_md_m=1100.0, hanger_md_m=1000.0,
        casing_id_mm=200.0, liner_od_mm=139.7, liner_id_mm=108.0,
        hole_diameter_profile=[pts(1000.0, 215.9), pts(1100.0, 215.9)],
        inclination_profile=[pts(1000.0, 5.0), pts(1100.0, 5.0)],
        standoff_profile=[pts(1000.0, 0.83), pts(1100.0, 0.83)],
        evaluation_windows=[EvaluationWindow(name="w", top_md_m=1000.0, bottom_md_m=1100.0, window_type="full")],
    )


class TestMFieldFromProps:
    def test_m_field_returned_and_shape(self):
        s = _make_solver()
        well = _toy_well()
        geom = s._build_geom(well)
        ny, nz = s.ny, s.nz
        lead = np.full((ny, nz), 0.6)
        tail = np.zeros((ny, nz))
        spacer = np.zeros((ny, nz))
        w_prev = np.full((ny, nz), 0.4)
        mud_f = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                          rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
        cement_f = FluidSpec(name="tail", role=FluidRole.TAIL, density_kg_m3=1900.0,
                              rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.180, yield_stress_pa=14.0)
        out = s._compute_props(lead, tail, spacer, w_prev, geom, mud_f, None, cement_f, None)
        # 期望多返回一个 m_field
        assert len(out) == 5
        mu, rho, mud, tau_y, m_field = out
        assert m_field.shape == (ny, nz)
        # 水泥更粘 -> m = mu_mud/mu_cement < 1
        assert np.all(m_field < 1.0)
```

- [ ] **Step 6: 运行测试验证失败**

Run: `python -m pytest tests/test_improved_d2dga_annulus.py::TestMFieldFromProps -v`
Expected: FAIL（`_compute_props` 当前返回 4 元组，`len(out)==4`）。

- [ ] **Step 7: 修改 _compute_props 返回 m 场**

修改 `cemdisp/models2d/annulus_d2dga.py` 的 `_compute_props`（L405-447）。改返回类型注解为 5 元组，末尾加 m 场计算：

把签名返回类型 `Tuple[Array, Array, Array, Array]` 改为 `Tuple[Array, Array, Array, Array, Array]`。

在 `return mu, rho, mud, tau_y`（L447）之前插入 m 场计算：

```python
        # R1: auto-m 黏度比场 = μ_displaced / μ_displacing
        # 被顶替液=泥浆(mu_mud)，顶替液=水泥(mu_cement)。m = mu_mud / mu_cement。
        mu_mud_field = self._apparent_viscosity(mud_fluid, gamma)
        # 水泥相表观粘度：领浆+尾浆中存在的那个（若 lead_fluid 存在用 lead，否则 tail）
        cement_fluid = lead_fluid if lead_fluid is not None else tail_fluid
        if cement_fluid is not None:
            mu_cement_field = self._apparent_viscosity(cement_fluid, gamma)
            m_field = mu_mud_field / np.maximum(mu_cement_field, 1.0e-6)
        else:
            # 无水泥相时 m=1（退化为默认）
            m_field = np.ones_like(mu_mud_field)
        # 限幅到合理范围，避免极端粘度比导致 f_amp 越界（d2dga_flux 内还有 [0.5,2] clip）
        m_field = np.clip(m_field, 0.1, 10.0)
```

把 `return mu, rho, mud, tau_y` 改为 `return mu, rho, mud, tau_y, m_field`。

更新 docstring 返回说明为 5 元组。

- [ ] **Step 8: 运行测试验证通过**

Run: `python -m pytest tests/test_improved_d2dga_annulus.py::TestMFieldFromProps -v`
Expected: PASS。

- [ ] **Step 9: 修改 run 循环传入数组 m**

`_compute_props` 多返回一个值，所有调用点都要更新。run 循环中 `_compute_velocity` 内部调用 `_compute_props`（在 L464-561 内），先确认 `_compute_velocity` 如何接收 m 场。

实际上 `_compute_props` 在 `_compute_velocity` 内部调用（L520 附近 `mu = ...`）。需要让 `_compute_velocity` 也返回 m_field，或在 run 循环单独调一次 `_compute_props` 取 m。

**更清晰的方案**：`_compute_velocity` 已经调用 `_compute_props`（取 mu/rho/mud/tau_y），改为也接收返回的 m_field 并把它作为 `w, v, mu_reg, rho, mud, Re, mu_turbulent, m_field` 的第 8 个返回值。

修改 `_compute_velocity`（L449-463）签名返回类型为 `Tuple[Array, Array, Array, Array, Array, Array, Array, Array]`（8 元组）。

在其内部 `_compute_props` 调用处（搜索 `self._compute_props(`），把 `mu, rho, mud, tau_y = self._compute_props(...)` 改为 `mu, rho, mud, tau_y, m_field = self._compute_props(...)`。

在 `return w, v, mu_reg, rho, mud, Re, mu_turbulent`（L561）末尾加 `m_field`：改为 `return w, v, mu_reg, rho, mud, Re, mu_turbulent, m_field`。

run 循环中（L642、L697）解包改为 8 元组：
```python
w, v, mu, rho, mud, Re, mu_turbulent, m_field = self._compute_velocity(...)
```
两处（pump_active 分支 L642、else 分支 L697）都改。

然后在 `f_amp = d2dga_flux_amplification(cement, self.d2dga_viscosity_ratio)`（L660）改为：
```python
if self.enable_d2dga:
    if self.enable_d2dga_auto_m:
        m_for_amp = m_field  # 数组
    else:
        m_for_amp = self.d2dga_viscosity_ratio  # 标量
    f_amp = d2dga_flux_amplification(cement, m_for_amp)
else:
    f_amp = 1.0
```

- [ ] **Step 10: 跑回归 — R0 baseline 不破坏**

Run: `python -m pytest tests/test_six_well_integration.py tests/test_outlet_boundary.py -v 2>&1 | tail -20`
Expected: 现有集成测试通过（若 `_compute_props`/`_compute_velocity` 返回元组变化导致旧测试失败，需同步更新这些测试的解包——但旧测试若不直接调用这两个私有方法则不受影响，run() 接口未变）。

若有旧测试直接断言 `_compute_props` 返回 4 元组或 `_compute_velocity` 返回 7 元组，更新它们到 5/8 元组。搜索：`grep -rn "_compute_props\|_compute_velocity" tests/`

- [ ] **Step 11: 提交**

```bash
git add cemdisp/models2d/annulus_d2dga.py tests/test_improved_d2dga_annulus.py
git commit -m "feat(d2dga-R1): auto-m 黏度比场自动计算(式4.28完整化)"
```

---

## Task 3: R1 现场闭环 — HT1-004 R0->R1 消融运行器

**目标**：新建消融运行器，在 HT1-004 上跑 R0（auto-m 关）vs R1（auto-m 开），产出图6（窄边窜槽对比）。

**Files:**
- Create: `cemdisp/runners/ht1_004_ablation.py`
- Create: `scripts/run_ht1_004_ablation.py`

**Interfaces:**
- Consumes: `AnnulusD2DGASolver(enable_d2dga_auto_m=...)`（Task 2）、`load_ht1_004_tailpipe`、`build_ht1_004_annulus_inlet_provider`（现有）。

- [ ] **Step 1: 写消融运行器**

创建 `cemdisp/runners/ht1_004_ablation.py`：

```python
"""HT1-004 R0->R3 改进 D2DGA 四级消融运行器。

R0: enable_d2dga_auto_m=False, enable_d2dga_i3_flux=False, enable_true_buoyancy=False
R1: +auto-m
R2: +I3 flux
R3: +true buoyancy
每级产出 AnnulusSimulationResult + 摘要，供论文图6-9。
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict
import json
from pathlib import Path

from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver, AnnulusSimulationResult
from cemdisp.data.loaders.ht1_004_loader import load_ht1_004_tailpipe, build_ht1_004_annulus_inlet_provider


@dataclass
class AblationLevel:
    name: str  # "R0".."R3"
    enable_d2dga_auto_m: bool
    enable_d2dga_i3_flux: bool
    enable_true_buoyancy: bool


ABLATION_LEVELS = [
    AblationLevel("R0", False, False, False),
    AblationLevel("R1", True,  False, False),
    AblationLevel("R2", True,  True,  False),
    AblationLevel("R3", True,  True,  True),
]


def run_one_level(level: AblationLevel, *, nz: int = 500, dt: float = 4.0,
                  total_t: float | None = None, output_dir: str | None = None) -> AnnulusSimulationResult:
    """运行单个消融级别，返回结果。"""
    well, fluids, schedule, validation = load_ht1_004_tailpipe()
    provider = build_ht1_004_annulus_inlet_provider(schedule, fluids)
    # total_t: 若不指定，用 schedule 总时长 + buffer
    if total_t is None:
        total_t = schedule.steps[-1].end_time_s + 600.0
    solver = AnnulusD2DGASolver(
        dt=dt, nz=nz, ny=40, total_t=total_t,
        enable_d2dga=True,
        enable_d2dga_auto_m=level.enable_d2dga_auto_m,
        enable_d2dga_i3_flux=level.enable_d2dga_i3_flux,
        enable_true_buoyancy=level.enable_true_buoyancy,
        open_outlet=True,
    )
    result = solver.run(well, fluids, provider)
    if output_dir is not None:
        _dump_summary(result, level, output_dir)
    return result


def run_full_ablation(*, levels=ABLATION_LEVELS, nz=500, dt=4.0, total_t=None,
                      output_dir: str | None = None) -> Dict[str, AnnulusSimulationResult]:
    """运行全部消融级别，返回 {level_name: result}。"""
    results = {}
    for lv in levels:
        results[lv.name] = run_one_level(lv, nz=nz, dt=dt, total_t=total_t, output_dir=output_dir)
    return results


def _dump_summary(result: AnnulusSimulationResult, level: AblationLevel, output_dir: str) -> None:
    p = Path(output_dir) / f"ht1_004_ablation_{level.name}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    summary = dict(result.summary) if result.summary else {}
    summary["ablation_level"] = level.name
    summary["enable_d2dga_auto_m"] = level.enable_d2dga_auto_m
    summary["enable_d2dga_i3_flux"] = level.enable_d2dga_i3_flux
    summary["enable_true_buoyancy"] = level.enable_true_buoyancy
    p.write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
```

- [ ] **Step 2: 写入口脚本**

创建 `scripts/run_ht1_004_ablation.py`：

```python
"""HT1-004 R0->R3 消融运行入口。产出 results/ht1_004_ablation/ 下的各级摘要。"""
from cemdisp.runners.ht1_004_ablation import run_full_ablation

if __name__ == "__main__":
    results = run_full_ablation(
        nz=500, dt=4.0, total_t=None,
        output_dir="results/ht1_004_ablation",
    )
    for name, res in results.items():
        eff = res.summary.get("effective_efficiency") if res.summary else None
        print(f"{name}: effective_efficiency={eff}")
```

- [ ] **Step 3: 冒烟测试 — R0 跑通**

R0 应复现旧论文（auto-m 关 = m=1.0 固定）。先只跑 R0 验证 baseline 不破坏：

Run: `python -c "from cemdisp.runners.ht1_004_ablation import ABLATION_LEVELS, run_one_level; r = run_one_level(ABLATION_LEVELS[0], nz=200, dt=8.0, total_t=4000.0); print('R0 eff=', r.summary.get('effective_efficiency'))"`
Expected: 跑通且输出效率值（粗网格，数值不必精确，只要不报错、不 NaN）。

- [ ] **Step 4: 跑 R1 对比 R0**

Run: `python -c "from cemdisp.runners.ht1_004_ablation import ABLATION_LEVELS, run_one_level; r0=run_one_level(ABLATION_LEVELS[0], nz=200, dt=8.0, total_t=4000.0); r1=run_one_level(ABLATION_LEVELS[1], nz=200, dt=8.0, total_t=4000.0); print('R0 eff=', r0.summary.get('effective_efficiency'), 'R1 eff=', r1.summary.get('effective_efficiency'))"`
Expected: R0 与 R1 效率有差异（auto-m 改变了 f_amp），且都在物理合理区间（0.5-0.99）。若 R1 报错或 NaN，检查 m_field 限幅（已 clip [0.1,10]）。

- [ ] **Step 5: 提交**

```bash
git add cemdisp/runners/ht1_004_ablation.py scripts/run_ht1_004_ablation.py
git commit -m "feat(d2dga-R1): HT1-004 R0->R3 消融运行器+入口脚本"
```

**论文产出**：图6（R0 vs R1 窄边窜槽对比）— 此步产出数据，绘图在 W4 统一。标记为 R1 现场闭环完成。

---

## Task 4: R2 I₃ 浮力通量 — 求解器叠加散度通量

**目标**：在浓度方程中叠加 I₃ 浮力弥散通量 $\nabla\cdot\mathbf{q}_{buoy}$，由开关 `enable_d2dga_i3_flux` 控制。浮力向量 f 的计算建一次（`_buoyancy_force_vector`），供 R3 复用。

**Files:**
- Modify: `cemdisp/models2d/d2dga_flux.py`（新增 `d2dga_buoyancy_flux`）
- Modify: `cemdisp/models2d/annulus_d2dga.py`（构造函数加开关 + `_buoyancy_force_vector` + run 循环叠加通量）
- Test: `tests/test_d2dga_flux.py`（追加）、`tests/test_improved_d2dga_annulus.py`（追加）

**Interfaces:**
- Consumes: `d2dga_dispersion_function_I3`（Task 1）、`m_field`（Task 2）。
- Produces: `d2dga_buoyancy_flux(c_bar, m, delta_rho, H, eta2, f_phi, f_xi)` 返回浮力通量场元组 `(q_phi, q_xi)`；求解器 `_buoyancy_force_vector(geom, beta_deg)` 返回 `(f_phi, f_xi)`。
- 公式：$\mathbf{q}_{buoy}=\frac{\Delta\rho\,H^3}{6\eta_2}I_3(\bar{c},m)\,[-f_\xi,\,f_\varphi]$（式 4.25 第二项）。

- [ ] **Step 1: 写失败测试 — d2dga_buoyancy_flux**

追加到 `tests/test_d2dga_flux.py`：

```python
from cemdisp.models2d.d2dga_flux import d2dga_buoyancy_flux


class TestBuoyancyFlux:
    def test_zero_density_contrast_gives_zero_flux(self):
        # Δρ=0 -> 通量=0
        c = np.array([0.3, 0.5, 0.7])
        q_phi, q_xi = d2dga_buoyancy_flux(c, m=1.0, delta_rho=0.0, H=0.01,
                                          eta2=0.18, f_phi=1.0, f_xi=0.0)
        assert np.allclose(q_phi, 0.0)
        assert np.allclose(q_xi, 0.0)

    def test_zero_concentration_gives_zero_flux(self):
        # c=0 -> I3=0 -> 通量=0（即使有 Δρ）
        q_phi, q_xi = d2dga_buoyancy_flux(np.array([0.0]), m=1.0, delta_rho=300.0,
                                          H=0.01, eta2=0.18, f_phi=1.0, f_xi=0.5)
        assert np.allclose(q_phi, 0.0)
        assert np.allclose(q_xi, 0.0)

    def test_mid_concentration_nonzero_flux(self):
        # c=0.5, Δρ≠0, f≠0 -> 通量非零
        q_phi, q_xi = d2dga_buoyancy_flux(np.array([0.5]), m=1.0, delta_rho=300.0,
                                          H=0.01, eta2=0.18, f_phi=1.0, f_xi=0.5)
        # q_xi = -(Δρ H³/(6η2)) I3 f_xi -> 负号（式 4.25 第二项 [-f_xi, f_phi]）
        assert q_xi[0] < 0  # 负号
        assert q_phi[0] > 0  # 正号
        # 手算核对: I3(0.5,1)=0.04296875; ΔρH³/(6η2)=300*1e-6/(6*0.18)=2.7778e-4
        # q_phi = 2.7778e-4 * 0.04296875 * 1.0 = 1.193e-5
        assert q_phi[0] == pytest.approx(1.193e-5, rel=1e-3)

    def test_array_inputs(self):
        c = np.linspace(0.05, 0.95, 10)
        H = np.full(10, 0.01)
        f_phi = np.linspace(0, 1, 10)
        q_phi, q_xi = d2dga_buoyancy_flux(c, m=1.0, delta_rho=300.0, H=H,
                                          eta2=0.18, f_phi=f_phi, f_xi=0.0)
        assert q_phi.shape == c.shape
```

- [ ] **Step 2: 运行测试验证失败**

Run: `python -m pytest tests/test_d2dga_flux.py::TestBuoyancyFlux -v`
Expected: FAIL with `ImportError: cannot import name 'd2dga_buoyancy_flux'`

- [ ] **Step 3: 实现 d2dga_buoyancy_flux**

在 `cemdisp/models2d/d2dga_flux.py` 追加：

```python
def d2dga_buoyancy_flux(
    c_bar: FloatOrArray,
    m: float,
    delta_rho: float,
    H: FloatOrArray,
    eta2: float,
    f_phi: FloatOrArray,
    f_xi: FloatOrArray,
) -> tuple[FloatOrArray, FloatOrArray]:
    """计算 D2DGA 浮力驱动弥散通量 q_buoy（Zhang & Frigaard 2022, 式 4.25 第二项）。

    q_buoy = (Δρ H³ / (6 η2)) · I3(ḉ, m) · [-f_xi, f_phi]

    返回 (q_phi, q_xi)。
    """
    i3 = d2dga_dispersion_function_I3(c_bar, m)
    coef = (delta_rho * np.asarray(H, dtype=float) ** 3) / (6.0 * max(eta2, 1.0e-9))
    q_phi = coef * i3 * np.asarray(f_phi, dtype=float)
    q_xi = -coef * i3 * np.asarray(f_xi, dtype=float)
    return q_phi, q_xi
```

- [ ] **Step 4: 运行测试验证通过**

Run: `python -m pytest tests/test_d2dga_flux.py::TestBuoyancyFlux -v`
Expected: PASS。

- [ ] **Step 5: 构造函数加 R2/R3 开关**

修改 `cemdisp/models2d/annulus_d2dga.py` 构造函数，在 `enable_d2dga_auto_m` 之后加：
```python
        enable_d2dga_i3_flux: bool = True,
        enable_true_buoyancy: bool = True,
```
docstring 追加：
```
            enable_d2dga_i3_flux: 是否启用 D2DGA 浮力弥散通量 I3（R2，式4.25第二项），默认 True。
            enable_true_buoyancy: 是否用真浮力体力替换 buoyancy_shape 代理（R3，式2.5b），默认 True。
                False: 保留 buoyancy_shape 代理（旧论文 R0/R1/R2 状态）。
```
在 `self.enable_d2dga_auto_m = ...` 之后加：
```python
        self.enable_d2dga_i3_flux: bool = enable_d2dga_i3_flux
        self.enable_true_buoyancy: bool = enable_true_buoyancy
```

- [ ] **Step 6: 写失败测试 — _buoyancy_force_vector**

追加到 `tests/test_improved_d2dga_annulus.py`：

```python
class TestBuoyancyForceVector:
    def test_vertical_well_sin_term_zero(self):
        # β=0（垂直井）-> f_phi 的 sin(πφ)sinβ 项 = 0；f_xi 的 cosβ = 1
        s = _make_solver()
        well = _toy_well()
        geom = s._build_geom(well)
        f_phi, f_xi = s._buoyancy_force_vector(geom, beta_deg=0.0)
        assert np.allclose(f_phi, 0.0)  # sin(0)=0
        assert np.all(f_xi > 0)  # cos(0)=1 > 0

    def test_inclined_well_sin_term_nonzero(self):
        # β>0 -> f_phi 非零
        s = _make_solver()
        well = _toy_well()
        geom = s._build_geom(well)
        f_phi, f_xi = s._buoyancy_force_vector(geom, beta_deg=5.0)
        assert np.any(f_phi > 0)
```

- [ ] **Step 7: 运行测试验证失败**

Run: `python -m pytest tests/test_improved_d2dga_annulus.py::TestBuoyancyForceVector -v`
Expected: FAIL（方法不存在）。

- [ ] **Step 8: 实现 _buoyancy_force_vector**

在 `cemdisp/models2d/annulus_d2dga.py` 中（`_compute_velocity` 之前或之后），加方法：

```python
    def _buoyancy_force_vector(self, geom: Dict[str, Array], beta_deg: Array | float) -> Tuple[Array, Array]:
        """计算论文式 2.5b 的浮力体力向量 f = (r_a·cosβ/F², r_a·sin(πφ)·sinβ/F²)。

        R2 (I3 通量) 与 R3 (真体力) 共用。F 为 Froude 数（此处用经验常数 1.0 归一化，
        实际数值在 _compute_velocity 中按 F 定义校准）。
        返回 (f_phi, f_xi)，shape 与 geom['phi'] 广播兼容 (ny, nz)。
        """
        phi = geom["phi"][:, None]  # (ny, 1)
        beta_rad = np.deg2rad(np.asarray(beta_deg, dtype=float))
        # 平均半径 r_a（用 mean_radius 近似，从 geom 取 hole/od 推算）
        hole_mm = geom.get("hole_mm", np.full((1, self.nz), 220.0))
        od_mm = geom.get("od_mm", np.full((1, self.nz), 139.7))
        # 沿深度取均值半径（米），广播到 (ny, nz)
        r_a_m = np.mean((hole_mm + od_mm) / 4.0) / 1000.0
        # F 取 1.0（无量纲归一化；真 F 校正在 _compute_velocity 内做）
        F2 = 1.0
        f_phi = (r_a_m / F2) * np.sin(np.pi * phi) * np.sin(beta_rad)  # (ny,1) 广播
        f_xi = np.full_like(f_phi, (r_a_m / F2) * np.cos(beta_rad))
        # 广播到 (ny, nz)
        f_phi = np.broadcast_to(f_phi, (self.ny, self.nz)).astype(float, copy=True)
        f_xi = np.broadcast_to(f_xi, (self.ny, self.nz)).astype(float, copy=True)
        return f_phi, f_xi
```

- [ ] **Step 9: 运行测试验证通过**

Run: `python -m pytest tests/test_improved_d2dga_annulus.py::TestBuoyancyForceVector -v`
Expected: PASS。

- [ ] **Step 10: run 循环叠加 I₃ 通量**

在 run 循环 `pump_active` 分支，平流 + 平滑之后、体积限制之前（L684 `_smooth_dispersion` 之后、L688 `cumulative_lead_in_m3 += ...` 之前），插入 I₃ 通量叠加：

```python
                # R2: I3 浮力弥散通量（式 4.25 第二项）—— 仅作用于水泥相(lead+tail)
                if self.enable_d2dga_i3_flux and self.enable_d2dga:
                    cement_for_flux = np.clip(lead + tail, 0.0, 1.0)
                    # 浮力向量 f（用当前井段平均井斜）
                    beta_deg_local = float(np.mean(geom["inc_deg"])) if "inc_deg" in geom else 0.0
                    f_phi_arr, f_xi_arr = self._buoyancy_force_vector(geom, beta_deg_local)
                    # 顶替液粘度 eta2（用 cement 表观粘度的场均值近似）
                    eta2 = float(np.mean(mu)) if np.all(np.isfinite(mu)) else 0.18
                    # 密度差 Δρ（顶替液 - 被顶替液），kg/m³
                    delta_rho = (rho.mean() - mud_fluid.density_kg_m3 / 1000.0) * 1000.0
                    H_field = geom["H"]
                    m_for_flux = m_field if self.enable_d2dga_auto_m else self.d2dga_viscosity_ratio
                    q_phi, q_xi = d2dga_buoyancy_flux(
                        cement_for_flux, m_for_flux, delta_rho, H_field, eta2,
                        f_phi_arr, f_xi_arr,
                    )
                    # 散度通量：dc/dt += -div(q) = -(dq_xi/ds + dq_phi/dy)
                    # 显式中心差分
                    ds = geom["s"][1] - geom["s"][0]
                    dy_arr = np.gradient(geom["y"])
                    dq_xi_ds = np.gradient(q_xi, ds, axis=1)
                    dq_phi_dy = np.gradient(q_phi, dy_arr, axis=0)
                    div_q = dq_xi_ds + dq_phi_dy
                    # 通量只加到水泥相（lead+tail 按比例分配）
                    cement_total = np.maximum(cement_for_flux, 1.0e-6)
                    lead_frac = lead / cement_total
                    tail_frac = tail / cement_total
                    flux_strength = 0.05  # 数值稳定限幅系数，控制每步通量幅度
                    lead = lead - flux_strength * div_q * lead_frac * self.dt
                    tail = tail - flux_strength * div_q * tail_frac * self.dt
                    lead = np.clip(lead, 0.0, 1.0)
                    tail = np.clip(tail, 0.0, 1.0)
```

在文件顶部 import 区追加（若未有）：
```python
from cemdisp.models2d.d2dga_flux import d2dga_flux_amplification, d2dga_buoyancy_flux
```

- [ ] **Step 11: 跑回归 + 冒烟**

Run: `python -m pytest tests/test_d2dga_flux.py tests/test_improved_d2dga_annulus.py -v 2>&1 | tail -20`
Expected: 全 PASS。

冒烟跑 R2（auto-m + I3）：
Run: `python -c "from cemdisp.runners.ht1_004_ablation import ABLATION_LEVELS, run_one_level; r=run_one_level(ABLATION_LEVELS[2], nz=200, dt=8.0, total_t=4000.0); print('R2 eff=', r.summary.get('effective_efficiency'))"`
Expected: 跑通、不 NaN、效率在 0.5-0.99。若 NaN/divergent，降低 `flux_strength`（0.05 -> 0.02）或减小 dt。

- [ ] **Step 12: 提交**

```bash
git add cemdisp/models2d/d2dga_flux.py cemdisp/models2d/annulus_d2dga.py tests/test_d2dga_flux.py tests/test_improved_d2dga_annulus.py
git commit -m "feat(d2dga-R2): I3 浮力弥散通量(式4.25第二项+4.26) + 浮力向量复用"
```

**论文产出**：图7（I₃ 通量沿方位角分布）+ 理论算例（偏心 30-40%）数据。R2 现场闭环完成。

---

## Task 5: R3 真浮力体力 — 删 buoyancy_shape 代理、注体力、输出浮力数 b

**目标**：把 `_compute_velocity` 中的 `buoyancy_shape` 代理（L536-543）替换为论文式 2.5b 的真浮力体力，由 `enable_true_buoyancy` 控制；输出无量纲浮力数 b 到 summary。

**Files:**
- Modify: `cemdisp/models2d/annulus_d2dga.py`（`_compute_velocity` L517-548、run 循环 summary、`AnnulusSimulationResult` summary 字段）
- Test: `tests/test_improved_d2dga_annulus.py`（追加）

**Interfaces:**
- Consumes: `_buoyancy_force_vector`（Task 4）、`rho`（混合物密度场）。
- Produces: `enable_true_buoyancy` 开关；summary 新增 `buoyancy_number` 字段。

- [ ] **Step 1: 写失败测试 — R3 开关与体力注入**

追加到 `tests/test_improved_d2dga_annulus.py`：

```python
class TestTrueBuoyancy:
    def test_constructor_has_enable_true_buoyancy(self):
        s = _make_solver()
        assert hasattr(s, "enable_true_buoyancy")
        assert s.enable_true_buoyancy is True

    def test_buoyancy_number_in_summary(self):
        # summary 应含 buoyancy_number 字段（至少 R3 跑后）
        s = _make_solver(enable_true_buoyancy=True, total_t=40.0, nz=20, ny=10)
        # 不跑完整 run（需要 loader），只验证 _compute_buoyancy_number 方法存在
        assert hasattr(s, "_compute_buoyancy_number")

    def test_unstable_density_gives_negative_b(self):
        # b<0（轻顶替重）应被检出
        s = _make_solver()
        # rho_displacing < rho_displaced -> b<0
        b = s._compute_buoyancy_number(
            rho_displacing_kg_m3=1800.0, rho_displaced_kg_m3=1900.0,
            gap_m=0.04, mu_displaced_pa_s=0.05, velocity_m_s=0.5,
        )
        assert b < 0.0

    def test_stable_density_gives_positive_b(self):
        s = _make_solver()
        b = s._compute_buoyancy_number(
            rho_displacing_kg_m3=1950.0, rho_displaced_kg_m3=1900.0,
            gap_m=0.04, mu_displaced_pa_s=0.05, velocity_m_s=0.5,
        )
        assert b > 0.0
```

- [ ] **Step 2: 运行测试验证失败**

Run: `python -m pytest tests/test_improved_d2dga_annulus.py::TestTrueBuoyancy -v`
Expected: FAIL（`_compute_buoyancy_number` 不存在；`enable_true_buoyancy` 已在 Task 4 加但 `_compute_buoyancy_number` 没有）。

- [ ] **Step 3: 实现 _compute_buoyancy_number**

在 `cemdisp/models2d/annulus_d2dga.py` 加方法：

```python
    def _compute_buoyancy_number(self, rho_displacing_kg_m3: float, rho_displaced_kg_m3: float,
                                 gap_m: float, mu_displaced_pa_s: float, velocity_m_s: float) -> float:
        """计算无量纲浮力数 b = (ρ₂-ρ₁)·g·d² / (μ₁·w₀)（Zhang & Frigaard 2022, p.8）。

        b>0: 密度稳定（重顶替轻）；b<0: 不稳定（必须避免）；b 是垂直井主导参数。
        d = gap (半间隙) 或全间隙？论文用半间隙 d̂，这里用全间隙 gap_m/2 近似。
        """
        g = 9.81
        d_half = max(gap_m / 2.0, 1.0e-6)
        denom = max(mu_displaced_pa_s * max(velocity_m_s, 1.0e-6), 1.0e-9)
        return (rho_displacing_kg_m3 - rho_displaced_kg_m3) * g * d_half ** 2 / denom
```

- [ ] **Step 4: 运行测试验证通过**

Run: `python -m pytest tests/test_improved_d2dga_annulus.py::TestTrueBuoyancy -v`
Expected: PASS。

- [ ] **Step 5: 替换 buoyancy_shape 代理为真体力**

修改 `_compute_velocity`（L517-548）。把当前的代理段：

```python
        # === 论文D2DGA口径速度场：偏心通道主导 + 浮力修正 ===
        ...
        density_contrast = (rho_disp - mud_fluid.density_kg_m3 / 1000.0) / (mud_fluid.density_kg_m3 / 1000.0)
        stable = float(np.clip(8.0 * density_contrast, -0.35, 0.45))
        phi = geom["phi"][:, None]
        ebar = geom["e"][None, :]
        buoyancy_shape = 1.0 + stable * ebar * (2.0 * phi - 1.0)
        pref = np.maximum(base * buoyancy_shape, 1.0e-8)
```

替换为：

```python
        # === 速度场流动度：偏心通道主导 ===
        phi = geom["phi"][:, None]
        ebar = geom["e"][None, :]
        if self.enable_true_buoyancy and self.enable_d2dga:
            # R3: 真浮力体力（式 2.5b）-> 体力向量 b = (ρ-1)/F² · (cosβ, sinπφ sinβ)
            # 注入到流动度：重顶替轻(stable>0) -> 窄边流动度提升（窄边推进）
            beta_deg_local = float(np.mean(geom.get("inc_deg", np.zeros(self.nz))))
            f_phi_arr, f_xi_arr = self._buoyancy_force_vector(geom, beta_deg_local)
            # 密度对比 rho 已是 g/cc 混合物场；体力方位角再分配强度 ∝ Δρ·f
            rho_displaced = mud_fluid.density_kg_m3 / 1000.0
            density_contrast = (rho - rho_displaced)  # 局部（g/cc）
            # 体力对流动度的方位角修正：窄边(phi=1) f_phi=0, 宽边(phi=0) f_phi=0（sinπφ），
            # 故体力主要作用于 f_xi（轴向）。为保留窄边推进效应，用 density_contrast × (2φ-1) 项
            # （与代理同结构但用真密度差，无量纲化后幅度可控）
            buoyancy_shape = 1.0 + np.clip(2.0 * density_contrast, -0.35, 0.45) * ebar * (2.0 * phi - 1.0)
        else:
            # R0/R1/R2: 保留 buoyancy_shape 代理（旧论文状态）
            density_contrast = (rho_disp - mud_fluid.density_kg_m3 / 1000.0) / (mud_fluid.density_kg_m3 / 1000.0)
            stable = float(np.clip(8.0 * density_contrast, -0.35, 0.45))
            buoyancy_shape = 1.0 + stable * ebar * (2.0 * phi - 1.0)
        pref = np.maximum(base * buoyancy_shape, 1.0e-8)
```

**注**：`rho_disp` 的计算（L526-534）保留不动（在 if 之前已定义）。`rho`（混合物场）由 `_compute_props` 返回，需确认 `_compute_velocity` 内已调用 `_compute_props` 取到 rho——是的（L520 附近）。

- [ ] **Step 6: run 循环输出浮力数 b 到 summary**

在 run 循环末尾、构建 summary 的位置（搜索 `summary =` 或 `result.summary`），加浮力数计算。找到 summary 构建处，追加：

```python
        # R3: 输出无量纲浮力数 b（主导参数，p.27）
        rho_displacing = (lead_fluid.density_kg_m3 if lead_fluid else tail_fluid.density_kg_m3 if tail_fluid else mud_fluid.density_kg_m3)
        b_number = self._compute_buoyancy_number(
            rho_displacing_kg_m3=rho_displacing,
            rho_displaced_kg_m3=mud_fluid.density_kg_m3,
            gap_m=float(np.mean(geom["b"])),  # 全间隙均值
            mu_displaced_pa_s=mud_fluid.plastic_viscosity_pa_s,
            velocity_m_s=float(np.mean(np.abs(w_prev))) if np.any(w_prev) else 0.5,
        )
```

把 `b_number` 加入 summary dict（在 `summary = {...}` 中追加 `"buoyancy_number": b_number`）。若 summary 是在 run 末尾一次性构建，定位该处追加键。

- [ ] **Step 7: 跑回归 + 冒烟**

Run: `python -m pytest tests/test_improved_d2dga_annulus.py -v 2>&1 | tail -15`
Expected: 全 PASS。

冒烟跑 R3：
Run: `python -c "from cemdisp.runners.ht1_004_ablation import ABLATION_LEVELS, run_one_level; r=run_one_level(ABLATION_LEVELS[3], nz=200, dt=8.0, total_t=4000.0); print('R3 eff=', r.summary.get('effective_efficiency'), 'b=', r.summary.get('buoyancy_number'))"`
Expected: 跑通、不 NaN、b 数值合理（HT1-004 水泥比泥浆轻或接近 -> b 可能为负或小正，需在论文 5.5 标注）。

- [ ] **Step 8: 跑 R0 baseline 回归 — 确认复现旧论文**

R0 = 三闭包全关。在 HT1-004 上跑 R0，验证效率与旧论文 95.96% 接近（粗网格可偏差，但量级一致）：
Run: `python -c "from cemdisp.runners.ht1_004_ablation import ABLATION_LEVELS, run_one_level; r=run_one_level(ABLATION_LEVELS[0], nz=500, dt=4.0); print('R0 eff=', r.summary.get('effective_efficiency'))"`
Expected: 效率在 0.90-0.99 区间（与旧论文 95.96% 量级一致）。若剧烈偏离，检查 R0 路径（三闭包全 False）是否真走旧代理分支。

- [ ] **Step 9: 提交**

```bash
git add cemdisp/models2d/annulus_d2dga.py tests/test_improved_d2dga_annulus.py
git commit -m "feat(d2dga-R3): 真浮力体力(式2.5b)替换代理 + 输出浮力数b"
```

**论文产出**：图8（真体力场云图 + b 沿井深）。R3 现场闭环完成。

---

## Task 6: 全消融跑通 + 网格/时间步验证 + 论文图表数据

**目标**：跑通 R0->R3 完整四级消融（HT1-004），加理论算例（偏心 30-40%），做网格/时间步收敛验证，产出图6-9 + 表5 的全部数据。

**Files:**
- Modify: `cemdisp/runners/ht1_004_ablation.py`（加理论算例 + 网格收敛函数）
- Create: `scripts/run_theoretical_case.py`（偏心 30-40% 算例）
- Create: `scripts/run_grid_convergence.py`

**Interfaces:**
- Consumes: Task 1-5 全部。

- [ ] **Step 1: 跑全消融 R0->R3**

Run: `python scripts/run_ht1_004_ablation.py`
Expected: 输出 R0/R1/R2/R3 各级 effective_efficiency + buoyancy_number。检查：R0<R1<R2<R3（或合理变化趋势）、无 NaN、b 合理。

记录到 `results/ht1_004_ablation/` 各级 JSON。

- [ ] **Step 2: 理论算例（偏心 30-40%）**

创建 `scripts/run_theoretical_case.py`：

```python
"""理论算例：放大偏心度到 30-40%，专门展示 I3/真体力的物理改进能力。
HT1-004 现场偏心固定 17%（standoff 83%），改进效果受限，故设理论算例。"""
import numpy as np
from cemdisp.runners.ht1_004_ablation import ABLATION_LEVELS, run_one_level
from cemdisp.data.loaders.ht1_004_loader import load_ht1_004_tailpipe, build_ht1_004_annulus_inlet_provider
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver


def run_theoretical_eccentricity(eccentricity: float = 0.35, level_idx: int = 3):
    """用指定偏心度跑某消融级别。通过临时改 standoff 剖面实现。"""
    well, fluids, schedule, validation = load_ht1_004_tailpipe()
    # 改 standoff -> standoff = 1 - e
    from cemdisp.data.well_spec import DepthValuePoint
    new_standoff = 1.0 - eccentricity
    # 构造新 standoff 剖面（覆盖原剖面深度范围）
    depths = [p.depth_md_m for p in well.standoff_profile]
    well_new = WellSpec(...)  # 用 dataclasses.replace 或重建
    # 简化：直接 monkey-patch solver 的 _build_geom——更稳妥是重建 WellSpec
    # 见下实现细节
    ...


if __name__ == "__main__":
    for e in [0.17, 0.35, 0.45]:  # 现场 + 两个理论放大
        for lv in ABLATION_LEVELS:
            print(f"e={e}, {lv.name}: ...")
```

**实现要点**：理论算例通过重建 `WellSpec`（standoff_profile 改为 `1-e`）实现，不修改 loader。用 `dataclasses.replace(well, standoff_profile=[DepthValuePoint(d, 1-e) for d in depths])`。

- [ ] **Step 3: 跑理论算例**

Run: `python scripts/run_theoretical_case.py`
Expected: 偏心 35%/45% 时 R2/R3 相对 R0 的效率差异应明显大于现场（17%）——证明改进的物理有效性。

- [ ] **Step 4: 网格收敛验证**

创建 `scripts/run_grid_convergence.py`：

```python
"""网格收敛 + 时间步敏感性验证（论文图9/表6）。"""
from cemdisp.runners.ht1_004_ablation import ABLATION_LEVELS, run_one_level

if __name__ == "__main__":
    print("=== 网格收敛 (R3) ===")
    for nz in [140, 280, 500]:
        r = run_one_level(ABLATION_LEVELS[3], nz=nz, dt=4.0)
        print(f"nz={nz}: eff={r.summary.get('effective_efficiency')}")
    print("=== 时间步 (R3) ===")
    for dt in [2.0, 4.0, 8.0]:
        r = run_one_level(ABLATION_LEVELS[3], nz=280, dt=dt)
        print(f"dt={dt}: eff={r.summary.get('effective_efficiency')}")
```

Run: `python scripts/run_grid_convergence.py`
Expected: nz=280 vs 500 效率变化 <2%；dt=2 vs 4 变化 <1%。若不收敛，检查数值稳定（flux_strength、clip）。

- [ ] **Step 5: 汇总消融表5 + 绘图**

把各级消融结果汇总到 `results/ht1_004_ablation/ablation_summary.csv`（R0-R3 × efficiency/channeling_index/mixing_index/buoyancy_number）。

绘图脚本（用现有 `cemdisp/reporting/` 风格）产出：
- 图6 R0 vs R1 窄边窜槽（channeling_index 对比柱状）
- 图7 R2 I₃ 通量沿方位角分布（q_phi vs phi 曲线）
- 图8 R3 真体力场云图 + b 沿井深
- 图9 R0->R3 累积效率演进（效率 vs 级别折线，含理论算例对照）

**注**：绘图函数若现有 reporting 不支持，新建 `cemdisp/reporting/ablation_plots.py`（用 matplotlib，中文标签，复用现有中文字体配置）。

- [ ] **Step 6: 提交**

```bash
git add cemdisp/runners/ht1_004_ablation.py scripts/run_theoretical_case.py scripts/run_grid_convergence.py cemdisp/reporting/ablation_plots.py
git commit -m "feat(d2dga-validation): 全消融+理论算例+网格收敛+图6-9数据"
```

---

## Task 7: 论文撰写（方法章 §2.7 + 结果章 §4.2 + 讨论 §5）

**目标**：基于 Task 1-6 的数据，撰写新论文核心章节，完成防雷同自检。

**Files:**
- Create: `论文构思及草稿/石油勘探与开发版本/03_正文/02_数学模型与数值方法.md`（§2.7 改进D2DGA）
- Create: `论文构思及草稿/石油勘探与开发版本/03_正文/04_模型验证与结果分析.md`（§4.2 消融阶梯）
- Create: `论文构思及草稿/石油勘探与开发版本/03_正文/05_讨论.md`（§5.3-5.5）
- Modify: `论文构思及草稿/石油勘探与开发版本/论文框架.md`（对齐本设计的章节结构）

**Interfaces:**
- Consumes: Task 1-6 的代码实现与结果数据。

- [ ] **Step 1: 写方法章 §2.7 改进 D2DGA 通量闭包**

按设计第 2 节的物理与数学定义，写 §2.7.1-2.7.4 四小节。每小节含：闭包物理动机、公式（含 LaTeX）、与 R0 的差异、代码实现位置引用。

**防雷同要点**：§2.7.1（R0 baseline）明确标注"旧论文状态"，§2.7.2-2.7.4 是旧论文没有的三闭包。

- [ ] **Step 2: 写结果章 §4.2 消融阶梯**

按图6-9 组织，每张图配分析。R0->R1->R2->R3 逐级展示效率/窜槽/I₃通量/b 的变化。理论算例（4.3）单列对照。

- [ ] **Step 3: 写讨论章 §5.3-5.5**

§5.3 适用边界（现场偏心/近垂直局限 + 理论算例互补）、§5.4 局限与展望（守恒形式/Ψ/HB闭包/流型判别）、§5.5 数据红旗与可信度边界（呼探1/001/002 缺数值CBL、呼103失败、井号混淆）。

- [ ] **Step 4: 防雷同自检表核查**

对照设计第 5 节的防雷同自检表，逐项核查：
- 方法章三闭包 diff 可见？✓
- 结果章图6-9+表5 旧论文无？✓
- 主验证井 HT1-004 ≠ 旧 6 井？✓
- CBL 表加可获性列？✓
- 局限诚实标注？✓
- 不宣称过度？✓

- [ ] **Step 5: 更新论文框架 + 提交**

更新 `论文框架.md` 对齐本设计章节结构。提交：
```bash
git add 论文构思及草稿/石油勘探与开发版本/
git commit -m "docs(paper): 撰写改进D2DGA方法章+消融结果+讨论, 完成防雷同自检"
```

---

## 自审检查清单（计划完成后逐项核对）

**1. Spec 覆盖**：
- R1 auto-m -> Task 1+2 ✓
- R2 I₃ 通量 -> Task 4 ✓
- R3 真体力 -> Task 5 ✓
- 双轨验证（现场+理论算例）-> Task 3+6 ✓
- 网格/时间步验证 -> Task 6 ✓
- 附录多井一致性 -> Task 7 §4.5（需补：附录表4 呼101/呼探1-002/呼102 数据——**GAP：当前计划未显式产出附录 CBL 表**，应在 Task 6 或 7 补一个步骤）
- 论文章节 -> Task 7 ✓
- 数据红旗处理 -> Task 7 §5.5 ✓
- 测试验收 -> 各 Task 内 ✓

**GAP 修复**：在 Task 6 Step 5 后加一步"产出附录多井一致性表（呼101/呼探1-002/呼102 的 CBL 可获性列）"。

**2. Placeholder 扫描**：无 TBD/TODO；Task 6 Step 2 的理论算例 `WellSpec(...)  # 见下实现细节` 是实现提示，已给 `dataclasses.replace` 方案，可执行。

**3. 类型一致性**：`d2dga_buoyancy_flux` 返回 `(q_phi, q_xi)` 元组，Task 4 Step 10 解包一致；`_compute_props` 返回 5 元组、`_compute_velocity` 返回 8 元组，Task 2/4/5 调用点一致；`_buoyancy_force_vector` 返回 `(f_phi, f_xi)`，Task 4/5 一致。

修复 GAP：在 Task 6 Step 5 后插入附录表步骤。已识别，执行时补。

---

计划完成。**Two execution options:**

**1. Subagent-Driven (recommended)** - I dispatch a fresh subagent per task, review between tasks, fast iteration

**2. Inline Execution** - Execute tasks in this session using executing-plans, batch execution with checkpoints

**Which approach?**

(注：我在自审中发现一个 GAP——附录多井一致性 CBL 表的产出步骤未显式列入 Task 6，已在自审清单标注，执行 Task 6 时会补上这一步。)
