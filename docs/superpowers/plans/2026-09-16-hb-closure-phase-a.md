# Phase A 实施计划：完整 B&F25 HB 闭包

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 把阶段 B 的 `ClosureProvider` 协议真正实现为 HB 闭包——一维间隙弱解（含塞流/static wall layer）→ `Ī₁/Ī₂/q₀` 查表 → (4.22) 椭圆方程非线性外迭代，使 2D 环空速度场承载完整 B&F25 非牛顿流变。

**Architecture:** 新建 `gap_solver.py` 实现 B&F25 附录 A 的增广拉格朗日/Uzawa 求解器（纯函数、无状态）；`hb_closure.HBClosure` 实现阶段 B 冻结的 `ClosureProvider` 协议（G 由外迭代注入，不改签名）；`stream_function` 增固定点外迭代；`annulus_d2dga` 经 `enable_hb_closure`（**默认 False**）接线。

**Tech Stack:** Python 3.13 / numpy / scipy / pytest；conda 环境 `shenjingwangluo`。

**Spec:** `docs/superpowers/specs/2026-09-16-hb-closure-phase-a-design.md`（含 B&F25 全式号锚点表）
**上游**：`docs/superpowers/plans/2026-09-16-hb-closure-route2-phase-b.md`（阶段 B，已合并 `main`）
**调查**：`docs/superpowers/research/2026-09-16-yield-gate-inertia.md`（屈服门根因，§5.1 依赖它）

## Global Constraints

- 运行测试/脚本必须设 `PYTHONIOENCODING=utf-8 PYTHONUTF8=1`。
- **L1 硬约束**：`n=1` 且 `τ_Y=0`（牛顿极限）时，全链路与 HEAD **逐位一致**；阶段 B 的牛顿极限护栏（`tests/contract/test_hb_newtonian_limit.py`）全程必须绿。
- **阶段 B 冻结的接口不得改**：`ClosureProvider.mobility/buoyant_mobility(c_bar, m, eta1, eta2, H)` 签名；`solve_stream_function(..., closure=, wall=)`。
- 新开关 `enable_hb_closure` **默认 False**；关闭时逐位 = HEAD。
- **`Ī₃` 仍用 Z&F22 (4.26)**（`two_layer.buoyancy_flux_distribution_i3` 不动）——B&F25 (2.27) 印刷式不自洽。
- **不得引 Z&F22/23 作 HB/斜井依据**；A 属"扩展应用，未获外部验证"。
- 公式**先写解析极限测试再实现**（阶段 B 的 C1 教训）。Uzawa 固定 `ρ=r=1`（满足 A31）。
- 每 Task 结束提交一次，message 用 `feat(hb):` / `test(hb):` 前缀（中文正文）。

---

### Task 1: `gap_solver.py` 骨架 + 牛顿/Bingham 解析锚（A-1a）

**Files:**
- Create: `cemdisp/models2d/gap_solver.py`
- Test: `tests/contract/test_gap_solver.py`

**Interfaces:**
- Produces:
  - `solve_fixed_G(c_bar: float, n: tuple[float,float], kappa: tuple[float,float], tau_y: tuple[float,float], G, Gb=(0.0,0.0), *, ny: int=201, r: float=1.0, tol: float=1e-10, max_iter: int=2000) -> GapSolution`
  - `GapSolution`（dataclass）：`I1, I2, q0: float`、`u: Array(ny,)`、`y: Array(ny,)`、`G`、`iters: int`、`converged: bool`
  - 网格：`y ∈ [0,1]` 归一化（论文 ỹ∈[0,1]，c̄ = 界面位置）

- [ ] **Step 1: 写解析极限测试（先失败）**

```python
# tests/contract/test_gap_solver.py
"""A-1：一维间隙弱解。先钉三个解析极限，再实现。"""
import numpy as np
import pytest

from cemdisp.models2d.gap_solver import solve_fixed_G


def test_newtonian_single_fluid_matches_slot_solution():
    """单流体牛顿（c̄=0）⇒ Ī₁=1/(3η)，槽流 u(y)=G(1−y²)/(2η)。"""
    sol = solve_fixed_G(c_bar=0.0, n=(1.0, 1.0), kappa=(2.0, 2.0),
                        tau_y=(0.0, 0.0), G=(1.0, 0.0))
    assert sol.converged
    assert sol.I1 == pytest.approx(1.0 / (3.0 * 2.0), rel=1e-8)
    y = sol.y
    assert np.allclose(sol.u, (1.0 - y**2) / (2.0 * 2.0), rtol=1e-6, atol=1e-8)


def test_newtonian_two_layer_matches_two_layer_closure():
    """牛顿两层 ⇒ Ī₁ 与 two_layer.mobility_i1 的归一化闭式一致。"""
    from cemdisp.models2d.two_layer import mobility_i1
    c, e1, e2 = 0.4, 3.0, 1.0
    m = e1 / e2
    sol = solve_fixed_G(c_bar=c, n=(1.0, 1.0), kappa=(e1, e2),
                        tau_y=(0.0, 0.0), G=(1.0, 0.0))
    expect = float(np.asarray(mobility_i1(c, m, eta1=e1, eta2=e2, H=1.0)))
    assert sol.I1 == pytest.approx(expect, rel=1e-7)


def test_bingham_plug_exists_and_matches_slot_solution():
    """Bingham 单流体：|τ|≤τY 区 du/dy=0（塞流），塞流起于 y=1−τY/|G|。"""
    G, tauY, kappa = 10.0, 3.0, 1.0
    sol = solve_fixed_G(c_bar=0.0, n=(1.0, 1.0), kappa=(kappa, kappa),
                        tau_y=(tauY, tauY), G=(G, 0.0))
    assert sol.converged
    y_plug = 1.0 - tauY / G
    du = np.gradient(sol.u, sol.y)
    mask = sol.y > y_plug + 1e-6
    assert np.max(np.abs(du[mask])) < 1e-6
```

- [ ] **Step 2: 跑测试确认失败**

Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest tests/contract/test_gap_solver.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'cemdisp.models2d.gap_solver'`

- [ ] **Step 3: 实现骨架（A.2.1 定 G 模式的 Uzawa，A16—A21）**

按 B&F25 附录 A.2.1：
- 交错网格：`u` 在格边（ny+1 点）、`q`/`λ` 在格心（ny 点）；
- `λ0`（A17）平衡 G/Gb；`λ̃^k = λ^k − λ0`；
- 迭代：① 由 (A18) 从壁面 `u(1)=0` 往回积分得 `u^{k+1}`；② 由 (A19)/(A20) 二分求 `θ` 得 `q^{k+1}`（`|m|≤τ_Y ⇒ q=0`）；③ (A21) 更新 `λ̃`；
- 收敛判据 `‖u^{k+1}−u^k‖_p`、`‖q^{k+1}−q^k‖_p`、`‖λ̃^k‖_p`，`p = 1+min(n₁,n₂)`。

```python
# cemdisp/models2d/gap_solver.py （结构示意——离散细节按 A.2.3 补全）
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class GapSolution:
    I1: float; I2: float; q0: float
    u: np.ndarray; y: np.ndarray; G: tuple; iters: int; converged: bool

def solve_fixed_G(c_bar, n, kappa, tau_y, G, Gb=(0.0, 0.0), *,
                  ny=201, r=1.0, tol=1e-10, max_iter=2000) -> GapSolution:
    """B&F25 A.2.1：固定 G/Gb 求 ũ，再算 Ī₁/Ī₂/q₀。"""
    ...
```

- [ ] **Step 4: 跑测试确认通过** — Expected: PASS（3 passed）
- [ ] **Step 5: 提交**

```bash
git add cemdisp/models2d/gap_solver.py tests/contract/test_gap_solver.py
git commit -m "feat(hb): A-1a 一维间隙 HB 弱解骨架（Uzawa 定G）+ 牛顿/Bingham 解析锚"
```

---

### Task 2: 幂律标度验收 + Uzawa 收敛性（A-1b）

**Files:** Test: `tests/contract/test_gap_solver.py`（追加）

- [ ] **Step 1: 写失败测试（幂律标度 + 收敛）**

```python
@pytest.mark.parametrize("n_exp,expect", [(1.0, 3.0), (0.7, 2 + 1/0.7), (0.4, 2 + 1/0.4)])
def test_single_fluid_power_law_I1_scales_as_H_pow_2_plus_1_over_n(n_exp, expect):
    """阶段 B 已数值核验：I₁ ∝ H^{2+1/n}。固定 G 扫描 H 验斜率。"""
    Hs = np.array([0.006, 0.008, 0.010, 0.012, 0.014])
    I1s = [solve_fixed_G(0.0, (n_exp, n_exp), (1.0, 1.0), (0.0, 0.0),
                         G=(1.0 / H, 0.0)).I1 * H**3 for H in Hs]
    slope = np.polyfit(np.log(Hs), np.log(np.array(I1s)), 1)[0]
    assert slope == pytest.approx(expect, abs=1e-3)


def test_uzawa_converges_and_reports_iters():
    sol = solve_fixed_G(0.5, (0.7, 0.7), (1.0, 1.0), (0.0, 0.0), G=(5.0, 0.0))
    assert sol.converged and 0 < sol.iters <= 2000
```

- [ ] **Step 2-4:** 跑失败 → 修实现（网格/H 标度/二分）→ 跑通过
- [ ] **Step 5: 复现 B&F25 图 6（L2 文献锚）**：探针脚本输出 `Ī₁(c̄)`/`q₀(c̄)` 随 `n=1/0.8/0.6/0.4/0.2`（m=1、w̄=1、b=0），与图 6 目视比对，结论写进报告。**不阻塞 Task 3，但必须做。**
- [ ] **Step 6: 提交**

```bash
git commit -am "feat(hb): A-1b 幂律 H 标度验收 + B&F25 图6 复现"
```

---

### Task 3: 定均速模式（A.2.2）——D2DGA 实际需要

**Files:** Modify: `cemdisp/models2d/gap_solver.py`；Test: 同上（追加）

**Interfaces:** Produces `solve_fixed_mean_velocity(c_bar, n, kappa, tau_y, u_bar, Gb=(0.0,0.0), **kw) -> GapSolution`（返回含解出的 `G`）

- [ ] **Step 1: 写失败测试（往返自洽）**

```python
def test_fixed_mean_velocity_is_self_consistent():
    """给定 ū* 反求 G，再用 G 正算 ū 必须回到 ū*。"""
    from cemdisp.models2d.gap_solver import solve_fixed_G, solve_fixed_mean_velocity
    args = dict(c_bar=0.45, n=(0.7, 0.8), kappa=(1.4, 0.9), tau_y=(2.0, 0.5))
    star = (0.30, 0.0)
    s = solve_fixed_mean_velocity(u_bar=star, **args)
    assert s.converged
    chk = solve_fixed_G(G=s.G, **args)
    assert chk.u.mean() == pytest.approx(star[0], rel=1e-6)
```

- [ ] **Step 2-4:** 跑失败 → 按 (A22)—(A30) 实现（`u=u_I+u_P`，α 由 (A28) 定 G）→ 跑通过
- [ ] **Step 5: 提交** — `git commit -am "feat(hb): A-1c 定均速模式（B&F25 A.2.2）——D2DGA 接口"`

---

### Task 4: HB 群换算 + `HBClosure`（A-2）

**Files:** Modify `two_layer.py`（增 `hb_groups(...)`，牛顿闭式不动）、`hb_closure.py`（增 `HBClosure`）；Test: `tests/contract/test_hb_closure_tabulated.py`

**Interfaces:**
- `two_layer.hb_groups(kappa1,kappa2,n1,n2,tauY1,tauY2,gamma0) -> dict`（(2.29)—(2.34) 的 `mu_e/tau_0/m/B/kappa_k/tau_Yk`）
- `hb_closure.HBClosure(n, kappa, tau_y, m, B, *, ny_gap=201)` 实现 `ClosureProvider`；`set_pressure_gradient(G)` 注入当前 G 场

- [ ] **Step 1: 写失败测试（群换算退化 + 牛顿极限）**

```python
def test_hb_groups_reduce_to_newtonian():
    g = hb_groups(kappa1=0.1, kappa2=0.05, n1=1.0, n2=1.0,
                  tauY1=0.0, tauY2=0.0, gamma0=100.0)
    assert g["B"] == pytest.approx(0.0)
    assert g["m"] == pytest.approx(0.1 / 0.05)


def test_hb_closure_newtonian_limit():
    from cemdisp.models2d.hb_closure import HBClosure, NewtonianClosure
    c = np.linspace(0, 1, 7); H = np.full(7, 0.012)
    hb = HBClosure(n=(1.0, 1.0), kappa=(1.0, 1.0), tau_y=(0.0, 0.0), m=1.0, B=0.0)
    assert np.allclose(hb.mobility(c, 1.0, 1.0, 1.0, H),
                       NewtonianClosure().mobility(c, 1.0, 1.0, 1.0, H), rtol=1e-6)
```

- [ ] **Step 2-5:** 跑失败 → 实现群换算 + `HBClosure`（内部调 `gap_solver` 定均速模式；按 8 井参数范围做 `(c̄, B)` 二维表 + 插值）→ 跑通过 → 提交

```bash
git commit -am "feat(hb): A-2 HB 无量纲群 (2.29)-(2.34) + HBClosure 查表（牛顿极限退回）"
```

---

### Task 5: 椭圆方程非线性外迭代（A-3a）

**Files:** Modify `stream_function.py`；Test `tests/contract/test_stream_function_nonlinear.py`

**Interfaces:** Produces `solve_stream_function_nonlinear(geom, c_bar, hb_closure, b_field, *, omega=0.5, tol=1e-6, max_outer=50) -> Array`

- [ ] **Step 1: 写失败测试（牛顿极限逐位）**

```python
def test_nonlinear_outer_loop_newtonian_limit_bitwise():
    """HBClosure 取牛顿极限 ⇒ 与线性 solve_stream_function 逐位一致。"""
    ...
    psi_lin = solve_stream_function(g, c, e1, e2, m, b)
    psi_nl = solve_stream_function_nonlinear(g, c, newtonian_limit_hb_closure, b)
    assert np.array_equal(psi_lin, psi_nl)
```

- [ ] **Step 2-5:** 跑失败 → 实现（冻结 Ī₁ → 解线性 Ψ → 由 ∇aΨ 经定均速反求局部 G → 更新 Ī₁ → 欠松弛 ω）→ 跑通过 → 提交

> ⚠️ **不收敛 ⇒ 抛错/告警并回退牛顿闭包**，不得静默。

---

### Task 6: `annulus_d2dga` 接线 + `enable_hb_closure`（A-3b）

**Files:** Modify `annulus_d2dga.py`；Test `tests/contract/test_hb_closure_wiring.py`

**Interfaces:** 求解器形参 `enable_hb_closure: bool = False`；可选第二开关 `hb_fix_cement_tau_y: bool = False`（供 Task 7 归因分离）

- [ ] **Step 1: 写失败测试（默认关逐位 = HEAD）**

```python
def test_default_off_is_bitwise_head():
    """enable_hb_closure 默认 False ⇒ 与不传该参数的求解器结果逐位一致。"""
```

- [ ] **Step 2-5:** 跑失败 → 接线（`FluidSpec` → HB 群换算；**水泥相 τ_y 不再被丢弃**，并由 `hb_fix_cement_tau_y` 控制可分离）→ 跑通过 → 提交

---

### Task 7: 三分量 A/B + 报告（**归因必须分开**）

**Files:** Create `scripts/entrypoints/run_hb_ablation_phase_a.py`、`results/HB闭包A-B_<date>/`

**依据（spec §5.1）**：A 的"换闭包"与"补水泥 τy"会同时出现，**必须分开报告**。

- [ ] **Step 1: 定义三分量变体**

```python
VARIANTS = {
  "H0_baseline":       {},                                                        # 阶段 B 默认
  "H1_closure_only":   {"enable_hb_closure": True,  "hb_fix_cement_tau_y": False}, # 只换 Ī₁
  "H2_cement_ty_only": {"enable_hb_closure": False, "hb_fix_cement_tau_y": True},  # 只补 τy
  "H3_both":           {"enable_hb_closure": True,  "hb_fix_cement_tau_y": True},  # 生产口径
}
```
（若两变化在实现上不可分离，须在报告中**明确说明并改用算子层分离**。）

- [ ] **Step 2:** 跑八井 × 4 变体，产出 CSV
- [ ] **Step 3: 验收**：**L1** `H0` 与阶段 B 默认逐位一致；**L3** Z&F22 Table 3 三门槛**不劣于现 1/3**；记录 `H3−H0` 的 Δη 并**按 §5.1 三分解**
- [ ] **Step 4:** 增补 `docs/源模型口径与适用域声明.md`：HB 闭包为"扩展应用、未获外部验证"
- [ ] **Step 5: 提交**

---

## Self-Review 记录

- **Spec 覆盖**：A-1→Task1-3、A-2→Task4、A-3→Task5-6；spec §4 四层锚落于 Task1(解析)/Task2(图6)/Task7(L3)；spec §5.1 的"三分量分开报告"落在 Task 7。
- **占位符扫描**：Task 3/4/5 的中间步骤以"跑失败→实现→跑通过→提交"表述并给出关键代码；离散细节以式号锚点给出，**执行者须对照 B&F25 原文**（`参考文档/_extracted_pdfs/...txt` 已入库）。
- **类型一致性**：`GapSolution` 字段、`solve_fixed_G/solve_fixed_mean_velocity` 签名、`HBClosure(n, kappa, tau_y, m, B)` 在 Task 1/3/4 间一致；`ClosureProvider` 沿用阶段 B 冻结签名。
- **未决**：`hb_fix_cement_tau_y` 的可分离性（Task 6/7）——若不可分离，Task 7 Step 1 备注给出降级方案。
