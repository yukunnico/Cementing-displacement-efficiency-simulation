# D2DGA Tier 1 (T1-1~T1-5) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在 R0-R3 开关版上去掉三处人工限幅/简化、补 I1/I2 分层修正、两层黏度闭包与 static wall，完成 D2DGA 核心物理修正（T1-1~T1-5）。

**Architecture:** 纯 L4 改动，只动 `cemdisp/models2d/d2dga_flux.py` + `cemdisp/models2d/annulus_d2dga.py` 及对应测试；不改 `AnnulusSimulationResult` 结构；R0-R3 三开关保留原语义做消融对照；数值安全靠 `c_safe` 防除零 + 局部 CFL 裁剪，不靠全局限幅。

**Tech Stack:** Python (numpy), pytest, cemdisp 包。

## Global Constraints

- Python 解释器：`/d/apps/Anaconda/python`，运行含中文输出脚本前设 `PYTHONIOENCODING=utf-8 PYTHONUTF8=1`。
- conda 环境：`shenjingwangluo`；安装 `pip install -e .`。
- 代码风格：中文 docstring + 类型注解 + 公式标注论文式号注释，与现有 `d2dga_flux.py`/`annulus_d2dga.py` 一致。
- 验证：⚠️ **不对标 CBL**，用 R0-R3 消融回归 + 六井集成 + Tier 0 诊断复跑。
- 默认 `nz=500, dt=4.0`（生产网格），单元测试用小网格 `nz=20, ny=10, total_t=40`。

---

## Task 1: T1-1 去通量放大人工限幅 [0.5,2.0]

**Files:**
- Modify: `cemdisp/models2d/d2dga_flux.py`（`d2dga_flux_amplification` 默认参数，行 27-28、40-41、50-51、78）
- Test: `tests/test_d2dga_flux.py`（新增 `TestNoClipDefault` 类）

**Interfaces:**
- Produces: `d2dga_flux_amplification(cement_fraction, viscosity_ratio=1.0, *, min_fraction=0.01, max_fraction=0.99, min_amplification=-np.inf, max_amplification=np.inf)` — 默认不物理裁剪，仅 `c_safe=clip(c,0.01,0.99)` 防除零；显式传 `min/max_amplification` 仍生效（向后兼容，muskat_regime 已用此路径）。

- [ ] **Step 1: Write the failing test**

在 `tests/test_d2dga_flux.py` 末尾追加：

```python
class TestNoClipDefault:
    """T1-1: 去通量放大限幅 [0.5,2.0]，默认不物理裁剪（仅 c_safe 防除零）。"""

    def test_large_m_high_c_not_clipped_low(self):
        # m=10, c=0.99: f -> 1/m = 0.1 (远小于 0.5)，旧默认 clip 到 0.5，新默认应 ≈0.1
        f = d2dga_flux_amplification(np.array([0.99]), viscosity_ratio=10.0)
        assert f[0] < 0.5  # 不被 clip 到 0.5
        assert f[0] == pytest.approx(0.1, abs=0.02)

    def test_c_zero_returns_near_1p5(self):
        # c=0 -> c_safe clip 到 0.01, f(0.01,1) ≈ 1.497（趋势趋向 1.5）
        f = d2dga_flux_amplification(0.0, viscosity_ratio=1.0)
        assert f == pytest.approx(1.5, abs=0.02)

    def test_explicit_clip_still_available(self):
        # 向后兼容：显式传 min/max_amplification 仍生效
        f = d2dga_flux_amplification(np.array([0.99]), viscosity_ratio=10.0,
                                     min_amplification=0.5, max_amplification=2.0)
        assert f[0] == 0.5  # 显式 clip 仍生效
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/test_d2dga_flux.py::TestNoClipDefault -v`
Expected: FAIL — `test_large_m_high_c_not_clipped_low` 断言 `f[0] < 0.5` 失败（旧默认 clip 到 0.5）。

- [ ] **Step 3: Write minimal implementation**

修改 `cemdisp/models2d/d2dga_flux.py`，把 `d2dga_flux_amplification` 的默认 `min_amplification`/`max_amplification` 从 `0.5`/`2.0` 改为 `-np.inf`/`np.inf`（两个 `@overload` 签名行 27-28、40-41 与主函数签名行 50-51 共三处）。行 78 的 `amplification = np.clip(numerator / denominator, min_amplification, max_amplification)` 保留（默认 ±inf 实际不裁剪，显式传值仍生效）。

行 50-51 改为：
```python
    min_amplification: float = -np.inf,
    max_amplification: float = np.inf,
```
（两个 overload 的对应行同样改。）

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/test_d2dga_flux.py -v`
Expected: PASS — 全部（含现有 `TestFluxAmplificationArrayM` 等，因现有用例 f∈[1.0,1.5] 不触 clip）。

- [ ] **Step 5: Commit**

```bash
git add cemdisp/models2d/d2dga_flux.py tests/test_d2dga_flux.py
git commit -m "feat(d2dga-T1-1): 去通量放大人工限幅[0.5,2.0],默认物理化(式4.28)"
```

---

## Task 2: T1-3a 提升 I1/I2 到 d2dga_flux.py 公共函数（DRY）

**Files:** Modify `cemdisp/models2d/d2dga_flux.py`（加 `d2dga_dispersion_I1`/`d2dga_dispersion_I2` + `import math`）；Modify `cemdisp/diagnostics/muskat_regime.py`（删私有 `_mean_mobility_I1`/`_buoyant_mobility_I2`，改从 `d2dga_flux` import，`_finger_velocity` 调用替换）；Test `tests/test_d2dga_flux.py`（加 `TestDispersionI1I2`）。

**Interfaces:** Produces `d2dga_dispersion_I1(c_bar, m)`=`[√m·c³+(1−c³)/√m]/3`（Bararpour 2025 式 2.24）；`d2dga_dispersion_I2(c_bar, m)`=`[2√m·c³(1−c)+c(1−c)²(1+2c)/√m]/6`（式 2.25）。标量→`float`，数组→`Array`。

**TDD 五步**（按 Task 1 模板：写失败测试→跑确认 FAIL→实现→跑确认 PASS→commit）：

- 测试 `TestDispersionI1I2`：`I1(0,m=4)=1/6`、`I1(1,m=4)=2/3`；`I2(0)=I2(1)=0`；`I2(0.5,1)>0`；数组形状兼容。
- 实现（`d2dga_flux.py` 顶部加 `import math`，末尾加两函数）：

```python
def d2dga_dispersion_I1(c_bar: FloatOrArray, m: float) -> FloatOrArray:
    """牛顿平均流动度 I1(ḉ,m)（Bararpour 2025 式 2.24，H³ 归一化）。"""
    c = np.asarray(c_bar, dtype=float)
    sq_m = math.sqrt(m)
    out = (sq_m * c**3 + (1.0 - c**3) / sq_m) / 3.0
    return float(out) if np.isscalar(c_bar) else out.astype(float, copy=False)


def d2dga_dispersion_I2(c_bar: FloatOrArray, m: float) -> FloatOrArray:
    """牛顿浮力流动度 I2(ḉ,m)（Bararpour 2025 式 2.25，H⁴ 归一化）。"""
    c = np.asarray(c_bar, dtype=float)
    sq_m = math.sqrt(m)
    out = (2.0 * sq_m * c**3 * (1.0 - c) + c * (1.0 - c)**2 * (1.0 + 2.0 * c) / sq_m) / 6.0
    return float(out) if np.isscalar(c_bar) else out.astype(float, copy=False)
```

- `muskat_regime.py`：删 `_mean_mobility_I1`/`_buoyant_mobility_I2`（194-211），import 处加 `d2dga_dispersion_I1, d2dga_dispersion_I2`，`_finger_velocity`（243-253）内 `_mean_mobility_I1`→`d2dga_dispersion_I1`、`_buoyant_mobility_I2`→`d2dga_dispersion_I2`。
- 跑：`PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/test_d2dga_flux.py tests/test_muskat_regime.py -v`，预期全过（muskat 现有测试不破坏）。
- commit：`git commit -m "feat(d2dga-T1-3a): 提升I1/I2到d2dga_flux公共函数(DRY,式2.24/2.25)"`

---

## Task 3: T1-2 I3 物理系数接线（去 flux_strength=0.05）

**Files:** Modify `cemdisp/models2d/annulus_d2dga.py`（775-779 I3 flux 更新段）；Test `tests/test_improved_d2dga_annulus.py`（加 `TestI3FluxPhysicalUpdate`）。

**Interfaces:** `lead/tail -= div_q · frac · dt`，`div_q` 由 `d2dga_buoyancy_flux` 物理系数 `ΔρH³/(6η₂)·I3`（式 4.25）直驱；局部 CFL 裁剪 `|div_q|·dt ≤ α·Δs`（α=0.5，防单步越界，非全局限幅）。`enable_d2dga_i3_flux` 开关语义不变。

**TDD 五步**（按 Task 1 模板）：

- 测试 `TestI3FluxPhysicalUpdate`：构造已知 `Δρ=300/H=0.01/η₂=0.18` 工况（参考 `test_d2dga_flux.py::TestBuoyancyFlux`），断言更新量量级 = 物理系数·梯度（`ΔρH³/(6η₂)·I3 ≈ 2.78e-4`）；断言"去 0.05 后单步更新量"大于"0.05 限幅版"。
- 实现（改 775-779）：

```python
# T1-2: 去人工限幅 flux_strength=0.05；物理系数 ΔρH³/(6η₂)·I3 直驱（式 4.25）
# 局部 CFL 裁剪防单步越界（α=0.5，非全局限幅）
alpha_cfl = 0.5
ds = geom["s"][1] - geom["s"][0]
step_limit = alpha_cfl * ds / max(self.dt, 1.0e-9)
div_q_clipped = np.clip(div_q, -step_limit, step_limit)
lead = lead - div_q_clipped * lead_frac * self.dt
tail = tail - div_q_clipped * tail_frac * self.dt
lead = np.clip(lead, 0.0, 1.0)
tail = np.clip(tail, 0.0, 1.0)
```

（原 `flux_strength = 0.05` 行删除；`div_q` 散度计算 763-772 保留。）
- 跑：`PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/test_improved_d2dga_annulus.py -v`，预期全过。
- commit：`git commit -m "feat(d2dga-T1-2): I3物理系数接线,去flux_strength=0.05限幅(式4.25)"`

---

## Task 4: T1-4 两层黏度闭包（式 4.23）+ 相黏度返回

**Files:** Modify `cemdisp/models2d/annulus_d2dga.py` `_compute_props`（424-464，返回 5→7 元组）+ `_compute_velocity`（566-570 `base`）；Test `tests/test_improved_d2dga_annulus.py`（更新 `test_m_field_returned_and_shape` + 加 `TestTwoLayerViscosity`）。

**Interfaces:** `_compute_props` 返回 `(mu, rho, mud, tau_y, m_field, eta1, eta2)`，`eta1=μ_mud_field`（泥浆相，约 453 行已算）、`eta2=μ_cement_field`（水泥相，约 455-457 行已算）。`_compute_velocity` 用 `1/η_mix = c̄³/η₂ + (1−c̄³)/η₁` 替换 `μ_reg` 作 `base` 分母。

**TDD 五步**（按 Task 1 模板）：

- 测试：`_compute_props` 返回 `len==7`（更新 `test_m_field_returned_and_shape`）；`TestTwoLayerViscosity`：`c̄=0 → η_mix=η1`、`c̄=1 → η_mix=η2`、混合值介于两者间。
- 实现 `_compute_props` 末尾 return（eta1/eta2 用已算的相黏度场）：

```python
return mu, rho, mud, tau_y, m_field, eta1, eta2
# eta1=μ_mud_field, eta2=μ_cement_field（约 453/455-457 已算，原仅用于 m_field，现额外返回）
```

- 实现 `_compute_velocity`（解构改 7 元组；`base` 行替换约 566-570，`mu_reg` 仍用于 Re/正则化保留）：

```python
mu, rho, mud, tau_y, m_field, eta1, eta2 = self._compute_props(...)  # 7 元组
...
b_mean = np.mean(b, axis=0, keepdims=True)
# T1-4: 两层黏度闭包 1/η_mix = c̄³/η₂ + (1−c̄³)/η₁（式 4.23）替换单相 μ_reg
c_bar = np.clip(lead + tail, 0.0, 1.0)  # 局部水泥浓度
eta_mix = 1.0 / (c_bar**3 / np.maximum(eta2, 1.0e-9)
                 + (1.0 - c_bar**3) / np.maximum(eta1, 1.0e-9))
base = (b / np.maximum(b_mean, 1.0e-12)) ** 2 / np.maximum(eta_mix, 1.0e-9)
```

- 同步：`run` 主循环内解构 `_compute_props` 的位置一并改 7 元组。
- 跑：`PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/test_improved_d2dga_annulus.py -v`，预期全过。
- commit：`git commit -m "feat(d2dga-T1-4): 两层黏度闭包1/η_mix(式4.23)+相黏度返回"`

---

## Task 5: T1-3b 体力向量注入流动度（式 2.5b/4.24）

**Files:** Modify `cemdisp/models2d/annulus_d2dga.py` `enable_true_buoyancy` 分支（587-599）+ 顶部 import；Test `tests/test_improved_d2dga_annulus.py`（加 `TestBuoyancyForceInjection`）。

**Interfaces:** `enable_true_buoyancy=True` 时 `pref = base·(1 + Δρ·f_phi·I2(c̄,m)/I1(c̄,m))`（式 4.24 形式），`f_phi` 由 `_buoyancy_force_vector(geom, beta_deg)` 给出；`False` 回退 `(2φ−1)` 简化代理。

**TDD 五步**（按 Task 1 模板）：

- 测试 `TestBuoyancyForceInjection`：合成均匀密度差（重顶替轻，`delta_rho>0`）+ 已知井斜，断言窄边（`phi` 大）`pref` 高于宽边（方位重分配方向正确）；`enable_true_buoyancy=False` 时回退 `(2φ−1)` 简化（与旧 R0-R2 一致）。
- 顶部 import：`from cemdisp.models2d.d2dga_flux import d2dga_dispersion_I1, d2dga_dispersion_I2`（追加到现有 `d2dga_flux` import 行）。
- 实现（改 587-599 `enable_true_buoyancy` 分支）：

```python
if self.enable_true_buoyancy and self.enable_d2dga:
    # T1-3b: 真浮力体力向量注入流动度（式 2.5b/4.24），替换 (2φ−1) 简化代理
    beta_deg_local = float(np.mean(geom.get("inc_deg", np.zeros(self.nz))))
    f_phi_arr, _ = self._buoyancy_force_vector(geom, beta_deg_local)
    rho_displaced = mud_fluid.density_kg_m3 / 1000.0
    delta_rho = (rho - rho_displaced)  # g/cc 局部密度差
    m_local = float(np.mean(m_field)) if np.all(np.isfinite(m_field)) else self.d2dga_viscosity_ratio
    c_bar = np.clip(lead + tail, 0.0, 1.0)
    i1 = d2dga_dispersion_I1(c_bar, m_local)
    i2 = d2dga_dispersion_I2(c_bar, m_local)
    # 式 4.24: 方位修正 = Δρ·f_phi·(I2/I1)；重顶替轻→窄边(f_phi 大) pref 提升
    correction = np.clip(delta_rho * (i2 / np.maximum(i1, 1.0e-12)), -0.5, 0.5)
    buoyancy_shape = 1.0 + correction * f_phi_arr
else:
    # R0/R1/R2 回退 (2φ−1) 简化代理
    density_contrast = (rho_disp - mud_fluid.density_kg_m3) / mud_fluid.density_kg_m3
    stable = float(np.clip(8.0 * density_contrast, -0.35, 0.45))
    buoyancy_shape = 1.0 + stable * ebar * (2.0 * phi - 1.0)
pref = np.maximum(base * buoyancy_shape, 1.0e-8)
```

- 跑：`PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/test_improved_d2dga_annulus.py -v`，预期全过。
- commit：`git commit -m "feat(d2dga-T1-3b): 体力向量注入流动度(式2.5b/4.24),替(2φ-1)简化"`

---

## Task 6: T1-5 static wall layer c_min 判据（式 2.35-2.41）

**Files:** Modify `cemdisp/models2d/annulus_d2dga.py`（`__init__` 加 `c_min` 约 190、`run` 内 wall 更新约 668、`_compute_velocity` 加 wall 反馈、去 `del wall` 约 627）；Test `tests/test_improved_d2dga_annulus.py`（加 `TestStaticWallLayer`）。

**Interfaces:** 求解器加 `c_min: float = 0.05`；`wall` 场按 `c < c_min` 处置 1（壁面静止层）；`_compute_velocity` 内 `pref *= (1−wall)`；`_depth_profiles` 去 `del wall`。

**TDD 五步**（按 Task 1 模板）：

- 测试 `TestStaticWallLayer`：`c > c_min` 全场 → `wall` 全零；`c < c_min` 区域 → `wall=1` 且该处 `w≈0`；`c_min=0.3` 时更多区域 `wall=1`。
- 实现 `__init__`（约 190）：加 `c_min: float = 0.05,` + `self.c_min: float = c_min`。
- 实现 `run` 内 wall 更新（668 区 `wall = np.zeros(...)` 替换）：

```python
# T1-5: static wall layer c_min 判据（Bararpour 2025 式 2.35-2.41）
# 局部水泥浓度 c < c_min 处壁面层泥浆滞留不流动 → wall=1
cement_local = np.clip(lead + tail, 0.0, 1.0)
wall = (cement_local < self.c_min).astype(float)
```

- 实现 `_compute_velocity` 加 `wall` 参数 + 反馈（`pref` 行后）：

```python
pref = np.maximum(base * buoyancy_shape, 1.0e-8)
pref = pref * (1.0 - wall)  # T1-5: wall=1 处壁面静止层，流动度归零
```

（`_compute_velocity` 签名加 `wall` 参数，`run` 传入更新后的 wall；`_depth_profiles` 删 `del wall` 行。）
- 跑：`PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/test_improved_d2dga_annulus.py -v`，预期全过。
- commit：`git commit -m "feat(d2dga-T1-5): static wall c_min判据激活wall场(式2.35-2.41)"`

---

## Task 7: 集成验证（R0-R3 消融 + 六井 + Tier 0 诊断复跑）

**Files:** 无新代码；运行验证。

- [ ] **Step 1: 全量单元测试**

Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/ -q`
Expected: PASS（含更新后的 `test_m_field_returned_and_shape` len==7、I1/I2、体力注入、static wall 测试，不破坏现有 Tier 0 诊断/六井测试）。

- [ ] **Step 2: R0-R3 消融回归**（不对标 CBL）

Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -c "from cemdisp.runners.ht1_004_ablation import run_one_level, ABLATION_LEVELS; [run_one_level(l, nz=500, dt=4.0) for l in ABLATION_LEVELS]"`
Expected: R0/R1/R2/R3 四级不崩溃；物理合理（`m>1`、η_N≤η_E、`t_br>0`、`flow_class` 非空）。R3 效率值会变（去限幅后，接受）。

- [ ] **Step 3: 六井集成**

Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/test_six_well_integration.py -v`
Expected: PASS（不崩溃 + 物理合理）。

- [ ] **Step 4: Tier 0 诊断复跑**

Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python scripts/run_tier0_diagnostics.py --well all`
Expected: `results/tier0_diagnostics/hu102_tier0.json` / `ht1_004_tier0.json` 重生成；Muskat regime/浮力分类/η_N 随物理化改善（b>0、flow_class 非空、η_N≤η_E、t_br>0）。

- [ ] **Step 5: 论文图6-9 重跑（R3 数据更新）**

Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python scripts/run_ht1_004_ablation.py`
Expected: 消融 CSV 重生成（R3 效率值更新）；重绘图6-9。

- [ ] **Step 6: 最终 commit**

```bash
git add -A
git commit -m "test(d2dga-T1): 集成验证通过(R0-R3消融+六井+Tier0诊断复跑,不对标CBL)"
```

---

## Self-Review（writing-plans 自检）

- **Spec 覆盖**：T1-1（Task1）/T1-2（Task3）/T1-3（Task2 I1/I2 + Task5 体力注入）/T1-4（Task4）/T1-5（Task6）全覆盖；集成验证（Task7）对应 spec §9 验收标准。
- **占位符**：无 TBD/TODO；每 Task 含实际测试代码 + 实现代码 + 命令。
- **类型一致**：`d2dga_dispersion_I1`/`d2dga_dispersion_I2` 在 Task2 定义、Task5 调用一致；`_compute_props` 7 元组在 Task4 定义、Task5 用 `m_field` 一致；`c_min` 在 Task6 定义一致。
- **已知不确定点**：`_compute_props` 内相黏度场变量名（`μ_mud_field`/`μ_cement_field` 约 453/455-457）需实现时按实际行号核实；`_buoyancy_force_vector` 返回 `f_phi` 形状需与 `buoyancy_shape` 广播兼容（实现时验证）。
