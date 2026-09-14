# 控压固井模型 · 源模型口径重构 + 仓库清理 实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 清理仓库脚本/测试/历史产物混乱，并按 Zhang & Frigaard (2022) D2DGA **原模型方法**重建环空 2D 求解器的速度—压力—浮力链路；**方法学取原文口径，数值输入紧贴现场 HB/斜井实际工况**。

**Architecture:**
- 速度场从"局部代数流动度 + Q 归一"改为**解流函数椭圆方程** `∇a·[(r_a/2I₁)∇aΨ + b] = 0`（Z&F22 (4.22)）。**关键简化**：I₁/I₂ 只依赖 c̄ 与 m（本步已知）⇒ 每步是**变系数线性 Poisson**，用 `scipy.sparse.linalg.spsolve` 直接解；论文的增广拉格朗日仅在 HB 非线性 `S(|∇Ψ|)` 时才需要，本计划用外层 Picard 迭代处理。
- 浮力拆成两阶：**领先阶轴向浮力数 b**（Z&F22 (4.14)/(4.19)，进动力学）+ **O(ε) 方位浮力**（(2.5b)，按 (2.6) 的 F² 定标）。
- 弥散**删除自创拉普拉斯项**，改用文献的 `q₀(c̄,m)` + `I₃(c̄,m)` 分层通量闭合。
- 新增三个聚焦模块：`stream_function.py`、`two_layer.py`、`buoyancy.py`；`annulus_d2dga.py` 保留为编排层。

**Tech Stack:** Python 3.13 / numpy / scipy.sparse / pandas / pytest；conda 环境 `shenjingwangluo`；Windows 需 `PYTHONIOENCODING=utf-8 PYTHONUTF8=1`。

**Spec（必读，执行者需同时读这两份）：**
- `D:\obsidian\obsidian-storage\1.科研\新疆油田控压固井项目\固井顶替效率改进\半环空设定与源文献适用域审查_2026-09-14.md`
- `D:\obsidian\obsidian-storage\1.科研\新疆油田控压固井项目\固井顶替效率改进\模型全面调研_早期版本对照与敏感性根因_2026-09-14.md`

## Global Constraints

- 代码注释/用户可见输出用**中文**；模块/函数/类名用**英文**（`AGENTS.md`）。
- **每个改动必须能指到文献式号**；指不到的一律不得进入论文"方法学依据"。
- 已确认正确的实现**不得改动**：`I₁/I₂` 闭式、`I₃` 用 **Z&F22 (4.26)**、`q₀` 放大因子 (4.28)、半环空 + `sinπφ` 方位浮力、`η_E/η_N` 定义、`b^{1+1/n}` 幂律缝隙律（出处 Walton & Bittleston 1991，**不是** Z&F22）。
- **权威结果目录 `results/<井名>_1D2D耦合模型/` 在 Task 13 之前不得改动**；其余 results 子目录只做归档。
- **删除一律改为归档**（用户裁定）：移到 `archive/` 或 obsidian，不做 `rm`。
- **禁用 `scripts/entrypoints/smoke_all_wells.py` 作为冒烟手段**（controller 裁定 R8）：该脚本会写入 `results/<井名>_1D2D耦合模型/`（权威目录），触碰「Task 13 前不得改动权威目录」。冒烟统一改用只读的 `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest tests/test_six_well_integration.py -q`（Task 2 完成后路径变为 `tests/contract/test_six_well_integration.py`）。
- 测试基线：**396 passed**（2026-09-14 实测 69.44s）。任何任务结束必须保持全绿（`tests/history/` 的"预期变红"除外，须逐条标注）。
- 每完成一个 Task 提交一次 commit（conventional commits，中文正文）。
- **不引入新第三方依赖**（scipy 已在环境中）。

---

# Phase 0 — 安全网与仓库清理

### Task 0: 建分支 + 冻结测试基线

**Files:**
- Create: `docs/superpowers/plans/baseline-2026-09-14.md`

**Interfaces:**
- Produces: 分支 `refactor/d2dga-source-fidelity`；基线数字（供后续任务回归比对）

- [ ] **Step 1: 记录既有改动**

Run:
```bash
cd "D:/users/desktop/research/控压固井项目/cement model" && git status --porcelain | head -20
```
预期：约 41 未跟踪 / 31 删除 / 37 修改（本会话前既有）。**记录，不回滚用户改动。**

- [ ] **Step 2: 建分支**

```bash
git checkout -b refactor/d2dga-source-fidelity
```
预期：`Switched to a new branch 'refactor/d2dga-source-fidelity'`

- [ ] **Step 3: 跑测试基线**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest tests/ -q --tb=short 2>&1 | tail -5
```
预期：`396 passed in ~70s`

- [ ] **Step 4: 记录基线**

写 `docs/superpowers/plans/baseline-2026-09-14.md`：commit hash、测试数、8 井当前 η_E（从 `results/<井名>_1D2D耦合模型/*_结果摘要.json` 的 `最终结果.全井段最终有效顶替效率` 提取）、基准算例 `mass_conservation_error`（`results/基准算例对照_2026-09-10/对照表.csv`）。

- [ ] **Step 5: Commit**

```bash
git add docs/superpowers/plans/baseline-2026-09-14.md
git commit -m "chore(baseline): 冻结 2026-09-14 测试与结果基线"
```

---

### Task 1: scripts/ 归档与索引

**Files:**
- Create: `scripts/README.md`、`scripts/{entrypoints,probes,reruns,plots,lib}/`
- Move: 见映射表
- Modify: 所有 `from scripts.xxx import ...` 与 `sys.path.insert` 调用点

**Interfaces:**
- Produces: `scripts/README.md` 单页索引（脚本名 | 用途 | 状态 | 输出目录 | 依赖）

**归档映射：**

| 目标 | 内容 |
|---|---|
| `entrypoints/` | `rerun_all_wells_corrected.py`、`run_tier0_diagnostics.py`、`run_theoretical_case.py`、`run_grid_convergence.py`、`export_depth_time_concentration.py`、`cbl_window_comparison.py`、`verify_loader_against_field_v1.py`、`calibrate_dispersion_yield.py`、`closure_contribution_scan.py`、`run_ablation_variants.py`、`run_gap_fill_variants.py`、`run_sensitivity_variants.py`、`test_all_wells.py`→改名 `smoke_all_wells.py` |
| `lib/` | `_mass_balance_diag_20260902.py`→`mass_balance_diag.py`、`_sensitivity_common.py`→`sensitivity_common.py` |
| `probes/` | 全部 `hu101_*`(10)、`test_dual_factor_ht1_004.py`、`test_geometry_rate_factor.py`、`test_standoff_reverse_ht1_004.py`、`test_standoff_sensitivity.py`、全部 `_diag_*`/`_audit_*`/`_pure_plug_*`/`_wallfreeze_*`/`_yieldgate_*`/`bisect_*`/`c_verify_*`/`corrected_ref_*`/`debug_*`/`verify_fix_*`/`_make_comparison_*`/`_rerun8_after_conservation_fix_*`、`dump_*`、`isolate_fix_mechanisms.py`、`analyze_profiles.py`、`extract_profiles.py`、`compare_old_new_model.py`、`hu101_rheology_sensitivity_20260829.py`、`analyze_calibrated_rerun_20260829.py`、`i3_ablation_rerun_20260829.py`、`dispersion_scale_sensitivity_scan.py`、`density_contrast_sensitivity_scan.py`、`standoff_sensitivity_scan.py`、`ht1_003_sensitivity.py`、`ht1_004_sensitivity.py`、`analyze_distortion_fix.py`、`dispersion_scale_sensitivity_nz250.py`、`standoff_sensitivity_c4.py` |
| `reruns/` | 全部 `rerun8_*`(6)、`rerun_stop_fix_20260901.py` |
| `plots/` | `plot_patent_figures.py`、`plot_ht1_001_*`、`plot_ht1_004_*`、`plot_wells_*`、`plot_hu101_cbl_*`、`plot_hu101_full_depth_*` |
| **改为 `*.DEPRECATED.py`（不移动，仅改名以断 import）** | `p3_p4_integration.py`、`run_calibrated_all_wells.py` |

- [ ] **Step 1: 记录移动前的 import 关系**

```bash
grep -rn "from scripts\.\|import scripts\.\|sys.path.insert" --include=*.py . | grep -v __pycache__
```
预期：`_mass_balance_diag_20260902` 被 6 处 import、`_sensitivity_common` 被 3 处。**全部记录。**

- [ ] **Step 2: 执行归档**

```bash
mkdir -p scripts/{entrypoints,probes,reruns,plots,lib} archive/scripts_untracked_2026-09-14
# 未跟踪文件先备份（不可逆）
cp scripts/hu101_*.py scripts/run_ablation_variants.py scripts/run_gap_fill_variants.py \
   scripts/run_sensitivity_variants.py archive/scripts_untracked_2026-09-14/
# 已跟踪用 git mv（逐个按映射表执行）
git mv scripts/rerun_all_wells_corrected.py scripts/entrypoints/
git mv scripts/_mass_balance_diag_20260902.py scripts/lib/mass_balance_diag.py
git mv scripts/_sensitivity_common.py scripts/lib/sensitivity_common.py
# ...其余按映射表逐条 git mv
```
预期：`ls scripts/` 只剩 5 个子目录

- [ ] **Step 3: 修正 import 路径**

把 `from scripts._mass_balance_diag_20260902 import X` → `from scripts.lib.mass_balance_diag import X`；`_sensitivity_common` 同理。检查各脚本的 `sys.path.insert(0, ...)` 仍指向 repo 根。

- [ ] **Step 4: 验证无断裂导入**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -c "
import pathlib, ast
bad=[]
for p in pathlib.Path('scripts').rglob('*.py'):
    try: ast.parse(p.read_text(encoding='utf-8'))
    except Exception as e: bad.append((str(p), str(e)))
print('语法错误:', bad or '无')
"
```
预期：`语法错误: 无`

- [ ] **Step 5: 写 `scripts/README.md`**（活跃面只列 `entrypoints/` 13 个 + `lib/` 2 个）

- [ ] **Step 6: 冒烟 + Commit**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest tests/test_six_well_integration.py -q 2>&1 | tail -5
# 注意：禁用 smoke_all_wells.py（R8：它会覆写权威结果目录）
git add -A scripts/ archive/ .gitignore
git commit -m "chore(scripts): 归档为 entrypoints/probes/reruns/plots/lib + README 索引"
```

---

### Task 2: tests/ 分层 + results/ 归档

**Files:**
- Create: `tests/contract/`、`tests/history/`、`tests/README.md`
- Move: tests 文件（见映射）；`results/*` 历史目录 → `archive/results_2026-09-14/`

**Interfaces:**
- Produces: `tests/README.md` 说明"契约测试"（守护当前正确性）vs"历史锁定测试"（改模型预期变红）

**分层映射：**

| 目标 | 测试文件 |
|---|---|
| `tests/contract/` | `test_well_spec.py`、`test_fluid_spec.py`、`test_pumping_schedule.py`、`test_provenance.py`、`test_evaluation_windows.py`、`test_d2dga_flux.py`、`test_displacement_metrics.py`、`test_m0_metrics.py`、`test_sync_cards.py`、`test_export_depth_time_shares.py`、`test_tier0_schedule_wiring.py`、`test_shoe_timeline.py`、`test_casing_flow.py`、`test_boundary_bridge.py`、`test_outlet_boundary.py`、`test_regime_closure.py`、`test_regime_classifiers.py`、`test_muskat_regime.py`、`test_flow_classification.py`、`test_axial_dispersion.py`、`test_yield_deadzone.py`、`test_six_well_integration.py` |
| `tests/history/` | `test_improved_d2dga_annulus.py`、`test_three_fixes_20260906.py`、`test_plug_semantics_restart.py`、`test_pipe_capacity_chain_fix.py`、`test_casing_mixing_contact_time.py`、`test_casing_mixing_buoyancy.py`、`test_m1_dispersion.py`、`test_hu101_loader_standoff.py`、`test_zhang2022_benchmark.py` |

- [ ] **Step 1: 建目录 + `git mv`**（按上表逐条）
- [ ] **Step 2: 修 pytest 发现路径** —— `pyproject.toml`/`pytest.ini` 检查 `testpaths`；必要时加 `--import-mode=importlib` 或给子目录加 `__init__.py`

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest tests/ -q --tb=no 2>&1 | tail -3
```
预期：`396 passed`（数量不变）

- [ ] **Step 3: 在 `tests/history/` 各文件头插入标注**

```python
# STATUS: history —— 本文件锁定"重构前"的历史行为契约。
# 源模型口径重构（2026-09-14）后，相关断言预期变红；红灯不等同回归失败，
# 需按 docs/superpowers/plans/2026-09-14-d2dga-source-fidelity-and-repo-cleanup.md 逐条复核。
```

- [ ] **Step 4: results/ 归档（只动非权威目录）**

**保留原地**：`results/<井名>_1D2D耦合模型/`（8 个）、`results/基准算例对照_2026-09-10/`。
**其余全部**移到 `archive/results_2026-09-14/`（写显式清单，避免误伤）。完成后 `ls results/ | wc -l` 预期为 10。

- [ ] **Step 5: 写 `tests/README.md` + Commit**

```bash
git add -A tests/ results/ archive/
git commit -m "chore(tests,results): tests 分层(contract/history) + results 历史目录归档"
```

---

# Phase 1 — 浮力道：口径统一 + 定标 + 进动力学

> 文献：(4.14) `(Hv̄, Hw̄) = −I₁·G + I₂·(Δρ/(H r_a))·(−f_ξ, f_φ)`；(2.5b) `b = r_a(ρ−1)/F²·(cosβ, sinπφ sinβ)`；(2.6) `F = sqrt(τ̂₀/(ρ̂₁ĝδ₀r̂ₐ*))` ⇒ **`F² = τ̂₀/(ρ̂₁ĝδ₀r̂ₐ*)`**（不是它的倒数），`τ̂₀ = μ̂₁ŵ₀/d̂`。

### Task 3: 统一浮力口径（新建 `buoyancy.py`）

**Files:**
- Create: `cemdisp/models2d/buoyancy.py`
- Test: `tests/contract/test_buoyancy.py`
- Modify: `cemdisp/models2d/annulus_d2dga.py:1412-1424`

**Interfaces:**
- Produces:
  - `displacing_density_kg_m3(lead_fluid, tail_fluid, mud_fluid) -> float` —— **全仓唯一**顶替液密度定义
  - `fluid_apparent_viscosity(fluid, shear_rate) -> float` —— 幂律/HB 用 `K·γ̇^(n−1)`，**不静默回退**
  - `buoyancy_number(rho_displacing, rho_displaced, half_gap_m, mu_displaced, w0_mps) -> float`
  - `froude_squared(mu_displaced, w0_mps, half_gap_m, rho_displaced, gap_scale_m, mean_radius_m) -> float`

- [ ] **Step 1: 写失败测试**

```python
# tests/contract/test_buoyancy.py
import pytest
from cemdisp.models2d.buoyancy import (buoyancy_number, froude_squared,
                                       fluid_apparent_viscosity)
from cemdisp.data.fluid_spec import FluidSpec, FluidRole, RheologyModel

def test_buoyancy_number_formula_and_sign():
    b = buoyancy_number(2100.0, 2020.0, 0.0257, 0.066, 0.563)
    assert b == pytest.approx(80 * 9.81 * 0.0257 ** 2 / (0.066 * 0.563), rel=1e-12)
    assert b > 0

def test_density_inversion_is_negative():
    assert buoyancy_number(1900.0, 1960.0, 0.045, 0.058, 0.764) < 0

def test_power_law_viscosity_not_silent_fallback():
    pl = FluidSpec(name="泥浆", role=FluidRole.MUD, density_kg_m3=1960.0,
                   rheology_model=RheologyModel.POWER_LAW,
                   power_law_n=0.72, consistency_k=0.381)
    v = fluid_apparent_viscosity(pl, shear_rate=10.0)
    assert v == pytest.approx(0.381 * 10.0 ** (0.72 - 1.0), rel=1e-12)
    assert v != 0.05

def test_froude_squared_matches_manual_scale():
    f2 = froude_squared(mu_displaced=0.058, w0_mps=0.764, half_gap_m=0.0458,
                        rho_displaced=1200.0, gap_scale_m=0.0917,
                        mean_radius_m=0.1071)
    assert 1e-4 < f2 < 1.0     # 实测量级 ~4.7e-3，绝不是 1.0
```

- [ ] **Step 2: 跑测试确认失败** → `ModuleNotFoundError: cemdisp.models2d.buoyancy`

- [ ] **Step 3: 实现**

```python
# cemdisp/models2d/buoyancy.py
"""浮力口径统一模块（Zhang & Frigaard 2022）。

全仓唯一的浮力定义入口，消除 annulus_d2dga.py 内并存的两套顶替液口径。
文献锚点：
- 浮力数 b = (ρ̂₂−ρ̂₁)·ĝ·d̂²/(μ̂₁·ŵ₀)        —— Z&F22 p.8
- Froude 数 F² = ρ̂₁·ĝ·δ₀·r̂ₐ*/(μ̂₁·ŵ₀/d̂)     —— Z&F22 (2.6)，τ̂₀ = μ̂₁ŵ₀/d̂
"""
from __future__ import annotations

from cemdisp.data.fluid_spec import FluidSpec, RheologyModel

G = 9.81
LEAD_WEIGHT = 0.67   # 体积加权口径，与 annulus_d2dga._compute_velocity 一致


def displacing_density_kg_m3(lead_fluid, tail_fluid, mud_fluid) -> float:
    """顶替液代表密度（0.67×领浆 + 0.33×尾浆）。全仓唯一口径。"""
    if lead_fluid is not None and tail_fluid is not None:
        return LEAD_WEIGHT * lead_fluid.density_kg_m3 + (1 - LEAD_WEIGHT) * tail_fluid.density_kg_m3
    if lead_fluid is not None:
        return float(lead_fluid.density_kg_m3)
    if tail_fluid is not None:
        return float(tail_fluid.density_kg_m3)
    return float(mud_fluid.density_kg_m3)


def fluid_apparent_viscosity(fluid: FluidSpec, shear_rate: float) -> float:
    """表观黏度。幂律/HB 用 K·γ̇^(n−1)；回退到 PV。缺参数时抛错，不静默回退。"""
    g = max(float(shear_rate), 1e-8)
    if fluid.rheology_model in (RheologyModel.POWER_LAW, RheologyModel.HERSCHEL_BULKLEY):
        if fluid.consistency_k is not None and fluid.power_law_n is not None:
            mu = float(fluid.consistency_k) * g ** (float(fluid.power_law_n) - 1.0)
            if fluid.plastic_viscosity_pa_s:
                mu = max(mu, float(fluid.plastic_viscosity_pa_s))
            return mu
    if fluid.plastic_viscosity_pa_s is not None:
        return float(fluid.plastic_viscosity_pa_s)
    raise ValueError(f"流体 {fluid.name} 缺少可用流变参数，禁止静默回退")


def buoyancy_number(rho_displacing: float, rho_displaced: float, half_gap_m: float,
                    mu_displaced: float, w0_mps: float) -> float:
    """无量纲浮力数 b（Z&F22 p.8）。b>0 密度稳定；b<0 密度倒置（文献警告严格避免）。"""
    d = max(float(half_gap_m), 1e-9)
    denom = max(float(mu_displaced) * max(float(w0_mps), 1e-9), 1e-12)
    return (float(rho_displacing) - float(rho_displaced)) * G * d * d / denom


def froude_squared(mu_displaced: float, w0_mps: float, half_gap_m: float,
                   rho_displaced: float, gap_scale_m: float, mean_radius_m: float) -> float:
    """Froude 数平方（Z&F22 (2.6)）。替代此前硬编码的 F2 = 1.0。"""
    d = max(float(half_gap_m), 1e-9)
    tau0 = float(mu_displaced) * max(float(w0_mps), 1e-9) / d
    denom = max(float(rho_displaced) * G * max(float(gap_scale_m), 1e-9)
                * max(float(mean_radius_m), 1e-9), 1e-12)
    return tau0 / denom
```

- [ ] **Step 4: 跑测试确认通过** → `4 passed`

- [ ] **Step 5: `annulus_d2dga.py:1412-1424` 改为调用新模块**（顶替液密度用 `displacing_density_kg_m3`；泥浆黏度用 `fluid_apparent_viscosity`；`gap_m` 改传几何半间隙而非经体积 scale 的 `mean(geom["b"])`）

- [ ] **Step 6: 回归 + Commit**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest tests/ -q --tb=short 2>&1 | tail -5
git add cemdisp/models2d/buoyancy.py tests/contract/test_buoyancy.py cemdisp/models2d/annulus_d2dga.py
git commit -m "feat(buoyancy): 统一浮力口径模块——消除双口径/静默回退/末步速度"
```

---

### Task 4: `F²` 按 (2.6) 定标

**Files:**
- Modify: `cemdisp/models2d/annulus_d2dga.py:759-780`
- Test: `tests/contract/test_buoyancy.py`（追加）

**Interfaces:**
- Consumes: `buoyancy.froude_squared`
- Produces: `_buoyancy_force_vector(geom, beta_deg, f2: float)` —— 新增 `f2` 形参，删除内部 `F2 = 1.0`

- [ ] **Step 1: 写失败测试**

```python
def test_buoyancy_force_vector_uses_froude_scale():
    """f_φ 必须带 1/F² 标定。f2 的取值按 (2.6) 用呼101 实参**现算**，
    不要硬抄——`F² = τ̂₀/(ρ̂₁ĝδ₀r̂ₐ*)`，`τ̂₀ = μ̂₁ŵ₀/d̂`，`d̂ = mean(geom["H"])`。
    参考量级 O(1e-2)。"""
    import numpy as np
    from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
    s = AnnulusD2DGASolver(ny=40, nz=2)
    geom = {"phi": np.linspace(0, 1, 40), "hole_mm": np.full((1, 2), 260.0),
            "od_mm": np.full((1, 2), 168.3)}
    f2 = 1.0e-2   # ← 用 (2.6) 现算后替换；仅作为量级示例
    f_unit, _ = s._buoyancy_force_vector(geom, 1.9, f2=1.0)
    f_phys, _ = s._buoyancy_force_vector(geom, 1.9, f2=f2)
    assert f_phys.max() > 50 * f_unit.max()
```

- [ ] **Step 2: 跑测试确认失败** → `TypeError: unexpected keyword argument 'f2'`

- [ ] **Step 3: 实现** —— `_buoyancy_force_vector` 加形参 `f2`，删 `F2 = 1.0`，两处 `/F2` 改 `/max(f2, 1e-12)`；两个调用点（`_compute_velocity` 内、run 循环 I3 通量处）按 `froude_squared(...)` 计算实参传入。

- [ ] **Step 4: 跑测试确认通过**

- [ ] **Step 5: 记录量级变化** —— 跑 hu101 求解，记录 `buoyancy_shape` min/max（预期从 `1±0.00004` → `1±0.01` 量级），写入 `baseline-2026-09-14.md`

- [ ] **Step 6: Commit** `fix(buoyancy): F² 按 Z&F22 (2.6) 定标，废止硬编码 F2=1.0`

---

### Task 5: 领先阶轴向浮力数接进动力学

**Files:**
- Modify: `cemdisp/models2d/annulus_d2dga.py:900-930`（`_compute_velocity` 的 buoyancy 段）
- Test: `tests/contract/test_buoyancy_axial.py`

**Interfaces:**
- Consumes: `buoyancy.buoyancy_number`、`d2dga_flux.d2dga_dispersion_I2`
- Produces: `_mobility_profile(c_bar, b_num, geom, i1_base, m_local, beta_deg) -> Array`（抽出的纯函数，便于单测）；`_compute_velocity` 的 `pref` 含 (4.14) 的 `I₂·(Δρ/(H·r_a))` 浮力项

- [ ] **Step 1: 写失败测试**

```python
# tests/contract/test_buoyancy_axial.py
import numpy as np
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver

def _geom(ny=40, nz=6):
    phi = np.linspace(0, 1, ny)
    H = np.broadcast_to(0.0458 * (1 + 0.3 * np.cos(np.pi * phi))[:, None], (ny, nz)).copy()
    return {"phi": phi, "H": H, "b": 2 * H, "y": np.linspace(0, np.pi * 0.1071, ny),
            "hole_mm": np.full((1, nz), 260.0), "od_mm": np.full((1, nz), 168.3),
            "s": np.linspace(0, 10.0, nz)}

def test_axial_buoyancy_raises_narrow_side_mobility():
    """b>0（重顶替轻）应抬高窄边流动度份额。"""
    s = AnnulusD2DGASolver(ny=40, nz=6)
    g = _geom()
    c = np.full_like(g["H"], 0.5)
    p0 = s._mobility_profile(c, b_num=0.0, geom=g, i1_base=0.44, m_local=0.5, beta_deg=0.0)
    p1 = s._mobility_profile(c, b_num=50.0, geom=g, i1_base=0.44, m_local=0.5, beta_deg=0.0)
    ny = g["H"].shape[0]
    narrow = slice(3 * ny // 4, ny)
    assert p1[narrow].mean() > p0[narrow].mean()

def test_zero_buoyancy_is_identity():
    s = AnnulusD2DGASolver(ny=40, nz=6)
    g = _geom()
    c = np.full_like(g["H"], 0.5)
    base = s._mobility_profile(c, b_num=0.0, geom=g, i1_base=0.44, m_local=0.5, beta_deg=0.0)
    assert np.allclose(base, s._mobility_profile(c, b_num=0.0, geom=g,
                                                 i1_base=0.44, m_local=0.5, beta_deg=0.0))
```

- [ ] **Step 2: 跑测试确认失败** → `AttributeError: _mobility_profile`

- [ ] **Step 3: 实现** —— 从 `_compute_velocity` 抽出流动度构造为 `_mobility_profile`，并按 (4.14)/(2.5b) 加入轴向浮力项：

```python
def _mobility_profile(self, c_bar, b_num, geom, i1_base, m_local, beta_deg):
    """局部流动度 pref（Z&F22 (4.14) + (2.5b)）。

    基础项：幂律缝隙律 b^(1+1/n)/η_mix（Walton & Bittleston 1991）
    浮力项：(4.14) 的 I₂·(Δρ/(H r_a))·(−f_ξ, f_φ)
      —— 以无量纲浮力数 b_num 为幅值、I₂/I₁ 为浓度权重、cosβ 为轴向投影。
      K_AXIAL 由 pref 与 b_num 同无量纲这一事实取 1（推导写入 docstring）。
    """
    ...
```

- [ ] **Step 4: 跑测试确认通过**

- [ ] **Step 5: 8 井快速扫描**（记录每井 b 与 η_E 相对重构前的位移）写入 `baseline-2026-09-14.md`

- [ ] **Step 6: Commit** `feat(buoyancy): 领先阶轴向浮力数 b 接入动力学（Z&F22 (4.14)）`

---

# Phase 2 — 两层闭包（文献口径）

### Task 6: 新建 `two_layer.py`

**Files:**
- Create: `cemdisp/models2d/two_layer.py`
- Modify: `cemdisp/models2d/d2dga_flux.py`（保留为 re-export 薄层，避免破坏既有 import）
- Test: `tests/contract/test_two_layer.py`

**Interfaces:**
- Produces:
  - `mobility_i1(c_bar, m, eta1, eta2, H) -> Array`（Z&F22 (4.21a)）
  - `mobility_i2(c_bar, m, eta1, eta2, H) -> Array`（(4.21b)）
  - `buoyancy_flux_distribution_i3(c_bar, m) -> Array`（**(4.26)**；**禁止**换成 B&F25 (2.27)）
  - `isotropic_flux_q0(c_bar, m) -> Array`（(4.28)）
  - `layer_thickness_fraction(c_bar) -> float`（`c̄ = y_i/H`）

- [ ] **Step 1: 写失败测试** —— 用 `m ∈ {0.2,0.5,1,2,5}` × `c̄ ∈ {0.25,0.5,0.75}` 全网格验证 (4.26) 与独立解析复算一致到 `1e-12`；验证 `I₁(0)=I₁(1)` 的牛顿退化；验证 `q0(0)=0, q0(1)=1`
- [ ] **Step 2: 跑测试确认失败**
- [ ] **Step 3: 实现**（迁移现有实现 + 新增 `q0` 与 `layer_thickness_fraction`）
- [ ] **Step 4: 跑测试确认通过**
- [ ] **Step 5: 全量回归** → `396 passed`
- [ ] **Step 6: Commit** `refactor(closure): two_layer.py 收敛 I1/I2/I3/q0 文献闭包`

---

### Task 7: 删自创弥散，改用分层通量闭合

**Files:**
- Modify: `cemdisp/models2d/annulus_d2dga.py:658-687`（删 `_smooth_dispersion`）、`:1192-1202`（删两处调用）、`__init__` 的 `dispersion_*` 形参
- Test: `tests/contract/test_no_invented_dispersion.py`

- [ ] **Step 1: 写失败测试**

```python
def test_no_laplacian_dispersion_in_concentration_update():
    """浓度更新中不得出现自创拉普拉斯平滑（Z&F22 p.11：no diffusive terms）。"""
    import inspect
    import cemdisp.models2d.annulus_d2dga as m
    assert "_smooth_dispersion" not in inspect.getsource(m.AnnulusD2DGASolver.run)
```

- [ ] **Step 2: 跑测试确认失败**
- [ ] **Step 3: 删除 `_smooth_dispersion` 与两处调用**；`dispersion_axial/azimuthal/dt_ref/dt_scale` 保留形参但置 `None` 并 `warnings.warn`（避免破坏 8 个 runner）；**通量闭合改由 Task 6 的 `q0`/`I3` 承担**
- [ ] **Step 4: 跑测试确认通过**
- [ ] **Step 5: 更正模块 docstring `:198`**"自创弥散=论文口径"的错误表述
- [ ] **Step 6: Commit** `fix(core): 删除自创拉普拉斯弥散，改用 q0/I3 分层通量（Z&F22）`

---

# Phase 3 — 流函数椭圆方程（原模型方法核心）

### Task 8: 新建 `stream_function.py`

**Files:**
- Create: `cemdisp/models2d/stream_function.py`
- Test: `tests/contract/test_stream_function.py`

**Interfaces:**
- Produces:
  - `solve_stream_function(geom, c_bar, eta1, eta2, m, b_field, *, ny=None, nz=None) -> Array`
    解 `∇a·[(r_a/2I₁)∇aΨ + b] = 0`；BC `Ψ(φ=0)=0`、`Ψ(φ=1)=1`（单位通量），调用方按 Q 缩放
  - `velocity_from_stream_function(psi, geom) -> tuple[Array, Array]`
    `w̄ = ∂φΨ/(2 r_a H)`、`v̄ = −∂ξΨ/(2 H)`（Z&F22 (2.2)）

- [ ] **Step 1: 写失败测试**

```python
# tests/contract/test_stream_function.py
import numpy as np
import pytest
from cemdisp.models2d.stream_function import (solve_stream_function,
                                              velocity_from_stream_function)

def _geom(ny=40, nz=5):
    y = np.linspace(0, np.pi * 0.1071, ny)
    phi = y / y[-1]
    H = np.broadcast_to(0.0458 * (1 + 0.3 * np.cos(np.pi * phi))[:, None], (ny, nz)).copy()
    return {"y": y, "phi": phi, "H": H, "b": 2 * H, "s": np.linspace(0, 10.0, nz),
            "hole_mm": np.full((1, nz), 260.0), "od_mm": np.full((1, nz), 168.3)}

def test_boundary_conditions_and_monotonicity():
    g = _geom()
    psi = solve_stream_function(g, np.zeros_like(g["H"]), np.ones_like(g["H"]),
                                np.ones_like(g["H"]), 1.0, np.zeros_like(g["H"]))
    assert psi[0, :] == pytest.approx(0.0, abs=1e-12)
    assert psi[-1, :] == pytest.approx(1.0, rel=1e-9)
    assert np.all(np.diff(psi[:, 0]) >= -1e-12)

def test_flux_is_depth_independent():
    """同心/无浮力时每列通量必须相同。"""
    g = _geom()
    c = np.full_like(g["H"], 0.5)
    psi = solve_stream_function(g, c, np.ones_like(c), np.ones_like(c), 1.0, np.zeros_like(c))
    w, _ = velocity_from_stream_function(psi, g)
    flux = np.trapezoid(2.0 * g["H"] * w, x=g["phi"], axis=0)
    assert np.allclose(flux, flux[0], rtol=1e-6)

def test_buoyancy_shifts_flow_to_narrow_side():
    g = _geom()
    c = np.full_like(g["H"], 0.5)
    z = np.zeros_like(c)
    w0, _ = velocity_from_stream_function(
        solve_stream_function(g, c, np.ones_like(c), np.ones_like(c), 1.0, z), g)
    w1, _ = velocity_from_stream_function(
        solve_stream_function(g, c, np.ones_like(c), np.ones_like(c), 1.0, z + 5.0), g)
    ny = g["H"].shape[0]
    narrow = slice(3 * ny // 4, ny)
    assert w1[narrow].mean() > w0[narrow].mean()
```

- [ ] **Step 2: 跑测试确认失败** → `ModuleNotFoundError`

- [ ] **Step 3: 实现** —— 变系数线性 Poisson，5 点格式，行主序编号 `idx = i*nz + j`，φ 方向两端为 Dirichlet（`Ψ=0` / `Ψ=1`），ξ 方向两端零梯度（一阶单侧）；`scipy.sparse.linalg.spsolve` 求解

```python
# cemdisp/models2d/stream_function.py
"""D2DGA 流函数椭圆方程求解（Zhang & Frigaard 2022, (2.3)/(4.22)）。

∇a·[ S + b ] = 0，  S = (r_a/2I₁)·∇aΨ，  ∇a = (1/r_a ∂φ, ∂ξ)

关键性质：I₁、I₂ 只依赖 c̄ 与 m（本时间步已知），故方程对 Ψ 线性——每步只需解
变系数线性 Poisson，不需论文的增广拉格朗日（那用于 HB 的 S(|∇aΨ|) 非线性）。
"""
```

- [ ] **Step 4: 跑测试确认通过**（3 个）
- [ ] **Step 5: 性能检查（ny=40,nz=250 单次 <0.15s）**
- [ ] **Step 6: Commit** `feat(stream): 流函数椭圆方程求解器（Z&F22 (2.3)/(4.22)，变系数线性 Poisson）`

---

### Task 9: 求解器接入 Ψ，替换代数流动度

**Files:**
- Modify: `cemdisp/models2d/annulus_d2dga.py:794-985`（`_compute_velocity`）、`:1162-1169`（f_amp 速度乘子）
- Test: `tests/history/test_zhang2022_benchmark.py`（预期部分变红，逐条复核）

**Interfaces:**
- Consumes: `stream_function.solve_stream_function`、`velocity_from_stream_function`
- Produces: `_compute_velocity` 的 `(w, v)` 由 Ψ 导出；新增开关 `enable_stream_function: bool = True`（False 走旧代数路径，保证可回退）

- [ ] **Step 1: 加开关 `enable_stream_function`（默认 True）**
- [ ] **Step 2: 新路径：`solve_stream_function` → `velocity_from_stream_function` → 按 `Q/2` 缩放 Ψ（半环空）**
- [ ] **Step 3: 废止 `w_d2dga = w * f_amp`（`:1168-1169`）** —— `f_amp` 只保留在通量层（由 Task 6 的 `q0` 承担）
- [ ] **Step 4: 基准算例跑一次，记录偏心度响应跨度**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest tests/history/test_zhang2022_benchmark.py -q --tb=short 2>&1 | tail -20
```

- [ ] **Step 5: 质量守恒检查** —— `mass_conservation_error` 目标从 0.33–0.48 降到 `<0.05`
- [ ] **Step 6: Commit** `feat(core): 速度场改用流函数椭圆方程解，废止 f_amp 速度乘子`

---

# Phase 4 — 几何与屈服

### Task 10: 消 `e_clip` 截断

**Files:**
- Modify: `cemdisp/models2d/annulus_d2dga.py:395-401`
- Test: `tests/contract/test_eccentricity_no_clip.py`

- [ ] **Step 1: 写失败测试** —— 断言 `standoff=0.35` 与 `0.45` 给出**不同**的 `e`（文献 `e∈[0,1)`，Pelipenko04 (2.1)）
- [ ] **Step 2: 跑测试确认失败**（现状 `e=clip(1−SO,0.05,0.55)` 使两者同为 0.55）
- [ ] **Step 3: 实现** —— `e = np.clip(1.0 - standoff, 1e-6, 1.0 - 1e-6)`；`e_clip_max`/`e_clip_measured_max`/`enable_e_clip_ruling` 保留形参并 `warnings.warn`
- [ ] **Step 4: 跑测试确认通过**
- [ ] **Step 5: 8 井居中度响应扫描**（standoff 0.35→0.80，确认 η_E 单调递增）
- [ ] **Step 6: Commit** `fix(geometry): 消 e_clip 截断，采用文献 e∈[0,1)（Pelipenko04 (2.1)）`

---

### Task 11: 屈服门 → 连续化

**Files:**
- Modify: `cemdisp/models2d/annulus_d2dga.py:624-655`（`_yield_gate_wall`）
- Test: `tests/contract/test_yield_gate_continuous.py`

- [ ] **Step 1: 写失败测试** —— 断言相邻 τw 的 `wall` 差 `< 0.5`（连续性），而非 0/1 跳变
- [ ] **Step 2: 跑测试确认失败**（现状 `np.where(immobile, 1.0, 0.0)` 二值）
- [ ] **Step 3: 实现** —— `wall = np.clip(1.0 - tau_w_extrap / (f_safety * tau_y), 0.0, 1.0)`（Pelipenko04 (2.6)-(2.8) 停流区判据的连续近似）
- [ ] **Step 4: 跑测试确认通过**
- [ ] **Step 5: 全量回归** —— `tests/history/` 预期多个红，逐条标注是否属"预期变红"
- [ ] **Step 6: Commit** `fix(yield): 屈服门二值→连续，消除阈值悬崖（Pelipenko04 (2.6)-(2.8)）`

---

# Phase 5 — 验证

### Task 12: 基准算例复现（数值正确性）

**Files:**
- Modify: `cemdisp/runners/zhang2022_benchmark.py` —— ⚠️ **该文件从未被 git 跟踪**（2026-09-14 核实），Task 2 已随自洽性提交补入版本控制；若仍为未跟踪，Task 12 须先 `git add`
- Create: `results/基准算例对照_2026-09-14/对照表.csv`

- [ ] **Step 1: 重跑 10 算例**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m cemdisp.runners.zhang2022_benchmark --out-dir results/基准算例对照_2026-09-14
```

- [ ] **Step 2: 逐项对照 Z&F22 Table 3 并填验收表**

| case | e | m | b | η_E 模型 | η_E 论文 | Δ | t_br 模型 | t_br 论文 | Δ | mass_err |
|---|---|---|---|---|---|---|---|---|---|---|

- [ ] **Step 3: 判定门槛**
  - `|Δη_E| ≤ 0.05` 的算例 ≥ **8/10**（当前 2/10）
  - `mass_conservation_error < 0.05`（当前 0.33–0.48）
  - **e 从 0.8→0.1 的 η_E 响应跨度 ≥ 20pp**（当前 1.3pp，论文 31pp）

- [ ] **Step 4: 未达标则写归因报告并回 Phase 3 迭代（最多 2 轮）**
- [ ] **Step 5: Commit** `test(benchmark): 源模型口径重构后基准算例复现度验收`

---

### Task 13: 8 井重跑 + 现场对照

**Files:**
- Create: `results/源模型口径重跑_2026-09-14/`（**先独立目录，经用户确认后再替换权威目录**）
- Create: `results/源模型口径重跑_2026-09-14/CBL窗口对照总表.csv`

- [ ] **Step 1: 8 井重跑（nz=250）**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/entrypoints/rerun_all_wells_corrected.py --out-dir results/源模型口径重跑_2026-09-14
```

- [ ] **Step 2: CBL 窗口对照**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/entrypoints/cbl_window_comparison.py --results results/源模型口径重跑_2026-09-14
```

- [ ] **Step 3: 与现场 CBL 的秩相关**（**用 CBL 评价窗口口径，不是全井口径**）
  靶值：呼101 62.77% / 呼102 66.65% / 呼103 12.06% / 呼1-003 78.7% / 呼1-004 0.3%。目标 `ρ ≥ +0.7`。

- [ ] **Step 4: 与重构前对照**

| 井 | 重构前 η_E | 重构后 η_E | 重构前 η_N | 重构后 η_N | CBL 靶值 | 窗口 η_E（新） | 偏差 pp |
|---|---|---|---|---|---|---|---|

- [ ] **Step 5: 更新权威目录**（**默认不替换，等用户确认**）
- [ ] **Step 6: Commit** `feat(results): 源模型口径 8 井重跑 + CBL 窗口对照`

---

# Phase 6 — 论文口径

### Task 14: 边界声明与配置表同步

**Files:**
- Create: `docs/源模型口径与适用域声明.md`
- Modify: `cemdisp/models2d/annulus_d2dga.py` 模块 docstring

- [ ] **Step 1: 写六条边界声明**（每条形影对应文献）
  1. 半环空 = 人为 Ψ 对称条件（Z&F22 p.25/p.32），抹掉 3-D 方位非对称
  2. 窄边=低边假设（Pelipenko04 p.4）
  3. 源模型严格适用域 = 竖直井 + 牛顿（Z&F23 p.34）；本项目 HB/斜井为**扩展应用，未获外部验证**
  4. 8 井 `b ∈ [−27, +74]`，**无一口达 `b ≳ 80`**；n 口密度倒置（`b<0`，文献警告区）
  5. 弥散为 `q₀/I₃` 解析闭包，**无自由弥散参数**
  6. `e ∈ [0,1)`，无截断

- [ ] **Step 2: 同步方法章公式出处表**（式号 → 代码锚点）
- [ ] **Step 3: Commit** `docs: 源模型口径与适用域声明 + 公式-代码对照表`

---

# 自检（Self-Review）

**Spec 覆盖**：§2 半环空 → T14；§3 适用域 → T12/13/14；§4 浮力定位 → T3/4/5；§5.1 做对的 → Global Constraints 禁止改动；§5.2 偏离 1 椭圆方程 → T8/9；偏离 2 浮力缺席 → T3/4/5；偏离 3 自创弥散 → T7；偏离 4 屈服门 → T11；偏离 5 e_clip → T10；偏离 6 两层 → T6；§7 路线 C → 全计划；§6 早期逻辑 → T3/7/10/11。

**显式范围外（下一轮）**：偏离 7（润滑小参数 ε 校验）、P3-1（HB 数值制表）、P2-2（M2 标定）、P2-3（两层 `y_i` 作为独立状态量）。

**类型一致性**：`buoyancy.py` → `buoyancy_number / froude_squared / fluid_apparent_viscosity / displacing_density_kg_m3`；`two_layer.py` → `mobility_i1 / mobility_i2 / buoyancy_flux_distribution_i3 / isotropic_flux_q0 / layer_thickness_fraction`；`stream_function.py` → `solve_stream_function / velocity_from_stream_function`。T4/5/8/9 的调用签名与之一致。

---

# 执行方式

**Subagent-Driven**：每 Task 派一个 fresh subagent；Task 完成后我做两阶段复核（代码审查 + 测试验证）再进下一个。

**并行化**：T1/T2 可并行；T3/T4 可并行；**T5/T7/T9/T10/T11 都改 `annulus_d2dga.py` —— 必须串行**，每个开始前确认上一 Task 已提交。

**风险闸门**：T9 与 T12 是最大风险点（速度场替换 + 基准验收）。T9 若使 `tests/history/` 大面积变红属预期；T12 若未达门槛，回 T8/T9 迭代最多 2 轮，仍不达标则如实报告并停在可回退点（`enable_stream_function=False` 即回旧口径）。
