# HB 闭包改造 · 阶段 B 实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 把 2D 环空默认速度场路径的闭包解耦成可插拔协议，接入屈服门（static wall layer 一阶物理）与幂律间隙一阶修正，并建立牛顿极限逐位回归护栏 + 八井 A/B 量化。

**Architecture:** 新增 `models2d/hb_closure.py` 定义 `ClosureProvider` 协议与牛顿实现；`solve_stream_function` 增加 `closure=`/`wall=` 两个可选注入点（缺省值 ⇒ 逐位等于 HEAD）；`annulus_d2dga` 通过两个 **opt-in** 开关（默认 `False`）暴露新物理，保证既有权威结果与全部测试不受影响。Phase A 的 HB 闭包将实现同一协议接入。

**Tech Stack:** Python 3.13 / numpy / scipy（`scipy.sparse.linalg.spsolve`）/ pytest；conda 环境 `shenjingwangluo`。

**Spec:** `docs/superpowers/specs/2026-09-16-hb-closure-route2-design.md`

## Global Constraints

- 运行测试/脚本必须设 `PYTHONIOENCODING=utf-8 PYTHONUTF8=1`。
- 基线 HEAD = `594a132`；**任何新增行为必须 opt-in 且默认关闭**，关闭时与 HEAD **逐位一致**。
- **L1 硬约束**：`n=1` 且 `τ_y=0`（牛顿极限）时，所有新路径与 HEAD 逐位一致。
- 新增开关命名：`enable_stream_yield_gate`、`enable_power_law_gap_correction`，默认均 `False`。
- 不改 `two_layer.py` 的既有函数签名/数值（仅新增可调用者）。
- 公式出处：`I₁` 幂律 H 指数 `2+1/n`（B&F25 (2.14)；幂律槽流 Walton & Bittleston 1991）；屈服判据 `wall`（Pelipenko04 (2.6)-(2.8)）。**不得引 Z&F22 作 HB 依据**。
- 每个 Task 结束提交一次；提交信息用 `feat(closure):` / `test(closure):` 前缀（中文正文）。

---

### Task 1: 闭包协议抽象 + 牛顿实现（零行为变更）

**Files:**
- Create: `cemdisp/models2d/hb_closure.py`
- Modify: `cemdisp/models2d/stream_function.py:169-298`（增 `closure=` 形参）
- Test: `tests/contract/test_hb_closure_protocol.py`

**Interfaces:**
- Consumes: `cemdisp.models2d.two_layer.mobility_i1/mobility_i2`（现有）。
- Produces:
  - `ClosureProvider`（`Protocol`）：`mobility(c_bar, m, eta1, eta2, H) -> Array`、`buoyant_mobility(c_bar, m, eta1, eta2, H) -> Array`
  - `NewtonianClosure()`：实现上述协议，逐位等于 `mobility_i1/i2`
  - `solve_stream_function(geom, c_bar, eta1, eta2, m, b_field, *, closure=None, wall=None, ny=None, nz=None) -> Array`

- [ ] **Step 1: 写失败测试**

```python
# tests/contract/test_hb_closure_protocol.py
"""B-1：闭包协议与牛顿实现（逐位等价 + 求解器注入无扰）。"""
import numpy as np

from cemdisp.models2d.hb_closure import ClosureProvider, NewtonianClosure
from cemdisp.models2d.two_layer import mobility_i1, mobility_i2


def test_newtonian_closure_implements_protocol():
    assert isinstance(NewtonianClosure(), ClosureProvider)


def test_newtonian_closure_bitwise_equals_two_layer():
    c = np.linspace(0.0, 1.0, 7)
    H = np.full(7, 0.0123)
    n = NewtonianClosure()
    assert np.array_equal(n.mobility(c, 0.7, 0.058, 0.171, H),
                          np.asarray(mobility_i1(c, 0.7, eta1=0.058, eta2=0.171, H=H)))
    assert np.array_equal(n.buoyant_mobility(c, 0.7, 0.058, 0.171, H),
                          np.asarray(mobility_i2(c, 0.7, eta1=0.058, eta2=0.171, H=H)))


def _tiny_geom(ny=5, nz=4):
    phi = np.linspace(0.0, 1.0, ny)
    s = np.linspace(0.0, 30.0, nz)
    H = np.full((ny, nz), 0.01)
    return {"phi": phi, "s": s, "H": H, "hole_mm": np.full(nz, 215.9),
            "od_mm": np.full(nz, 168.3), "y": np.linspace(0.0, np.pi * 0.1, ny),
            "b": np.full((ny, nz), 0.02)}


def test_solve_stream_default_equals_explicit_newtonian():
    """closure=None 与 closure=NewtonianClosure() 逐位一致。"""
    from cemdisp.models2d.stream_function import solve_stream_function
    g = _tiny_geom()
    c = np.full((5, 4), 0.5)
    b = np.zeros((2, 5, 4))
    psi_a = solve_stream_function(g, c, 0.058, 0.171, 0.34, b)
    psi_b = solve_stream_function(g, c, 0.058, 0.171, 0.34, b, closure=NewtonianClosure())
    assert np.array_equal(psi_a, psi_b)
```

- [ ] **Step 2: 跑测试确认失败**

Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest tests/contract/test_hb_closure_protocol.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'cemdisp.models2d.hb_closure'`

- [ ] **Step 3: 写实现**

```python
# cemdisp/models2d/hb_closure.py
"""流变闭包提供者协议（Phase B-1）。

把 stream_function 求解器与具体闭包解耦：求解器只依赖 ClosureProvider 协议。
NewtonianClosure 逐位等于 two_layer.mobility_i1/i2（Z&F22 (4.21a)/(4.21b)）。
Phase A 的 HB 查表闭包（HBClosure）实现同一协议，额外依赖 G（由外层迭代注入）。

注意（适用域）：源模型 Z&F22/23 严格适用于竖直井 + 牛顿流体；本协议的幂律/HB
实现属扩展应用，未获外部验证（见 docs/源模型口径与适用域声明.md 声明 3）。
"""
from __future__ import annotations

from typing import Protocol, runtime_checkable

import numpy as np
from numpy.typing import NDArray

from cemdisp.models2d.two_layer import mobility_i1, mobility_i2

Array = NDArray[np.float64]


@runtime_checkable
class ClosureProvider(Protocol):
    """两层闭包协议：返回间隙平均流动度 I₁ 与浮力流动度 I₂。

    约定（Z&F22 (4.21a)/(4.21b) 口径）：c_bar 为 (ny,nz) 场或标量；m、eta1、
    eta2 为标量；H 为 (ny,nz) 场或标量；返回与 c_bar/H 广播后的 float 数组。
    """

    def mobility(self, c_bar, m: float, eta1: float, eta2: float, H) -> Array: ...

    def buoyant_mobility(self, c_bar, m: float, eta1: float, eta2: float, H) -> Array: ...


class NewtonianClosure:
    """牛顿两层闭包（Z&F22 (4.21a)/(4.21b)）——逐位等于 two_layer 原函数。"""

    def mobility(self, c_bar, m: float, eta1: float, eta2: float, H) -> Array:
        return np.asarray(mobility_i1(c_bar, m, eta1=eta1, eta2=eta2, H=H), dtype=float)

    def buoyant_mobility(self, c_bar, m: float, eta1: float, eta2: float, H) -> Array:
        return np.asarray(mobility_i2(c_bar, m, eta1=eta1, eta2=eta2, H=H), dtype=float)
```

在 `stream_function.py`：
1. 模块 import 区加 `from cemdisp.models2d.hb_closure import ClosureProvider, NewtonianClosure`
2. 改签名（`:169`）为 `def solve_stream_function(geom, c_bar, eta1, eta2, m, b_field, *, closure=None, wall=None, ny=None, nz=None)`（`wall` 本步未用，Task 2 接入）
3. 改 `:224-225` 的 `I1 = np.asarray(mobility_i1(...))` 为：
```python
    if closure is None:
        closure = NewtonianClosure()
    I1 = np.asarray(closure.mobility(c, m, e1, e2, H), dtype=float)
```
（`e1/e2` 沿用已算好的标量 `_scalar_viscosity` 结果，保证逐位。）

- [ ] **Step 4: 跑测试确认通过 + 全量回归**

Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest tests/contract/test_hb_closure_protocol.py -v`
Expected: PASS（3 passed）

Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest tests/ -v`
Expected: 与改动前**同样的通过/预期红集合**（无新增失败）。

- [ ] **Step 5: 提交**

```bash
git add cemdisp/models2d/hb_closure.py cemdisp/models2d/stream_function.py tests/contract/test_hb_closure_protocol.py
git commit -m "feat(closure): B-1 闭包协议抽象 + NewtonianClosure（逐位等价，零行为变更）"
```

---

### Task 2: 屈服门接入椭圆算子（opt-in）

**Files:**
- Modify: `cemdisp/models2d/stream_function.py`（`conductance` 进 `a_cell`/`c_cell`）
- Modify: `cemdisp/models2d/annulus_d2dga.py`（`:263-299` 加形参；`:1159` `_velocity_stream_function` 增 `wall`；`:1437-1441` 传参）
- Test: `tests/contract/test_stream_yield_gate.py`

**Interfaces:**
- Consumes: Task 1 的 `solve_stream_function(..., wall=...)`；现有 `_yield_gate_wall`（`annulus_d2dga.py:714-772`）。
- Produces: 求解器形参 `enable_stream_yield_gate: bool = False`；冻结区电导地板常量 `_WALL_CONDUCTANCE_FLOOR = 1e-6`。

- [ ] **Step 1: 写失败测试**

```python
# tests/contract/test_stream_yield_gate.py
"""B-2：屈服门接入流函数算子（wall=None 逐位无扰；wall>0 冻结窄边）。"""
import numpy as np

from cemdisp.models2d.stream_function import solve_stream_function, _WALL_CONDUCTANCE_FLOOR


def _geom(ny=9, nz=4):
    return {"phi": np.linspace(0.0, 1.0, ny), "s": np.linspace(0.0, 30.0, nz),
            "H": np.full((ny, nz), 0.01), "hole_mm": np.full(nz, 215.9),
            "od_mm": np.full(nz, 168.3), "y": np.linspace(0.0, np.pi * 0.1, ny),
            "b": np.full((ny, nz), 0.02)}


def test_wall_none_is_bitwise_unchanged():
    g = _geom(); c = np.full((9, 4), 0.5); b = np.zeros((2, 9, 4))
    a = solve_stream_function(g, c, 0.058, 0.171, 0.34, b)
    z = solve_stream_function(g, c, 0.058, 0.171, 0.34, b, wall=np.zeros((9, 4)))
    assert np.array_equal(a, z)


def test_wall_freezes_narrow_side():
    """窄边（φ→1）wall=1 ⇒ 该侧轴向速度显著下降。"""
    from cemdisp.models2d.stream_function import velocity_from_stream_function
    g = _geom(); c = np.full((9, 4), 0.5); b = np.zeros((2, 9, 4))
    wall = np.zeros((9, 4)); wall[-3:, :] = 1.0
    psi0 = solve_stream_function(g, c, 0.058, 0.171, 0.34, b)
    psi1 = solve_stream_function(g, c, 0.058, 0.171, 0.34, b, wall=wall)
    w0, _ = velocity_from_stream_function(psi0, g)
    w1, _ = velocity_from_stream_function(psi1, g)
    assert np.mean(np.abs(w1[-3:])) < np.mean(np.abs(w0[-3:]))
    assert np.all(np.isfinite(psi1))


def test_conductance_floor_prevents_singular():
    g = _geom(); c = np.full((9, 4), 0.5); b = np.zeros((2, 9, 4))
    psi = solve_stream_function(g, c, 0.058, 0.171, 0.34, b, wall=np.ones((9, 4)))
    assert np.all(np.isfinite(psi))
```

- [ ] **Step 2: 跑测试确认失败**

Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest tests/contract/test_stream_yield_gate.py -v`
Expected: FAIL — `ImportError: cannot import name '_WALL_CONDUCTANCE_FLOOR'`

- [ ] **Step 3: 写实现**

`stream_function.py` 模块级常量（import 区后）：
```python
# B-2：屈服门冻结区的电导地板——防止 wall≡1 时系数全零导致矩阵奇异。
_WALL_CONDUCTANCE_FLOOR = 1.0e-6
```
求解器内，把 `:240-241` 的两行改为：
```python
    # B-2：屈服门电导因子（wall=None ⇒ 1.0，逐位无扰）
    if wall is None:
        conductance = 1.0
    else:
        w_arr = np.asarray(wall, dtype=float)
        if w_arr.shape != H.shape:
            raise ValueError("wall 形状须与 geom['H'] 相同 (ny,nz)")
        conductance = np.maximum(1.0 - w_arr, _WALL_CONDUCTANCE_FLOOR)
    # ⚠️ 方向（2026-09-16 Task 2 实测订正）：a_cell=1/(2I₁) 是导度，冻结须放大它
    # （I₁_eff→0 ⇒ a_cell→∞ ⇒ Ψ 平直 ⇒ w→0）。严禁写成 a_cell*(1−wall)。
    I1_eff = I1 * conductance
    a_cell = 1.0 / (2.0 * I1_eff)        # φ-槽系数 1/(2I₁_eff)
    c_cell = r_a / (2.0 * I1_eff)        # ξ-槽系数 r_a/(2I₁_eff)
```

`annulus_d2dga.py`：
1. `__init__` 形参表（紧随 `enable_power_law_gap_law: bool = True,`）加：
```python
        enable_stream_yield_gate: bool = False,  # B-2 opt-in：屈服门进流函数算子（默认关=HEAD 逐位）
```
并在存 `self.enable_power_law_gap_law`（`:440`）附近存属性：
```python
        self.enable_stream_yield_gate = enable_stream_yield_gate
```
2. `_velocity_stream_function` 签名（`:1159`）末尾加 `wall: Array | None = None`；`solve_stream_function(...)` 调用处（`:1338`）加 `wall=wall`。
3. `_compute_velocity` 新路径分支（`:1437-1441`）改为：
```python
        if self.enable_stream_function:
            w, v = self._velocity_stream_function(
                lead, tail, geom, q_m3s, w_prev, mud_fluid, lead_fluid, tail_fluid,
                wall=(wall if self.enable_stream_yield_gate else None),
            )
            return w, v, mu_reg, rho, mud, Re, mu_turbulent, m_field, tau_y, eta2, n_mix, kappa_mix
```

- [ ] **Step 4: 跑测试确认通过 + 全量回归**

Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest tests/contract/test_stream_yield_gate.py -v`
Expected: PASS（3 passed）

Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest tests/ -v`
Expected: 无新增失败（默认关 ⇒ 逐位=HEAD）。

- [ ] **Step 5: 提交**

```bash
git add cemdisp/models2d/stream_function.py cemdisp/models2d/annulus_d2dga.py tests/contract/test_stream_yield_gate.py
git commit -m "feat(closure): B-2 屈服门接入流函数算子（opt-in enable_stream_yield_gate，默认关逐位=HEAD）"
```

---

### Task 3: 幂律间隙一阶修正（opt-in）

**Files:**
- Modify: `cemdisp/models2d/hb_closure.py`（增 `PowerLawGapClosure`）
- Modify: `cemdisp/models2d/annulus_d2dga.py`（形参 + `_velocity_stream_function` 选闭包）
- Test: `tests/contract/test_power_law_gap_closure.py`

**Interfaces:**
- Consumes: Task 1 的 `ClosureProvider`、`NewtonianClosure`。
- Produces: `PowerLawGapClosure(n, base=None)`；求解器形参 `enable_power_law_gap_correction: bool = False`。

**指数依据**：`I₁ ∝ H^{2+1/n}G^{1/n−1}`（B&F25 (2.14) + 幂律槽流）；相对牛顿 `H³` 的场修正 `(H/H̄)^{1/n−1}`；`G` 纯标量因子在 flux 归一化中消去。`n=1` 时因子恒 1。

- [ ] **Step 1: 写失败测试**

```python
# tests/contract/test_power_law_gap_closure.py
"""B-3：幂律间隙一阶修正（n=1 逐位退化；n<1 增强偏心分流）。"""
import numpy as np

from cemdisp.models2d.hb_closure import PowerLawGapClosure, NewtonianClosure


def test_n_equals_one_is_bitwise_newtonian():
    c = np.linspace(0.0, 1.0, 7); H = np.linspace(0.008, 0.014, 7)
    a = PowerLawGapClosure(1.0).mobility(c, 0.7, 0.058, 0.171, H)
    b = NewtonianClosure().mobility(c, 0.7, 0.058, 0.171, H)
    assert np.array_equal(a, b)


def test_n_less_one_amplifies_gap_contrast():
    """n<1 时 I₁ 的 H 反差被放大：宽/窄流动度比增大。"""
    c = np.full(5, 0.5); H = np.array([0.006, 0.008, 0.010, 0.012, 0.014])
    i_newton = NewtonianClosure().mobility(c, 0.7, 0.058, 0.171, H)
    i_pl = PowerLawGapClosure(0.7).mobility(c, 0.7, 0.058, 0.171, H)
    assert (i_pl[-1] / i_pl[0]) > (i_newton[-1] / i_newton[0])
```

- [ ] **Step 2: 跑测试确认失败**

Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest tests/contract/test_power_law_gap_closure.py -v`
Expected: FAIL — `ImportError: cannot import name 'PowerLawGapClosure'`

- [ ] **Step 3: 写实现**

`hb_closure.py` 追加：
```python
class PowerLawGapClosure:
    """牛顿闭包 + 幂律间隙一阶修正（B-3，**近似**，非 B&F25 闭包口径）。

    I₁ → I₁ · (H/H̄)^{1/n − 1}（H̄ = 全场均值）。依据：单流体幂律槽流
    ū ∝ H^{1+1/n}G^{1/n} ⇒ I₁ = Hū/G ∝ H^{2+1/n}G^{1/n−1}；相对牛顿 H³ 的
    场修正即 (H/H̄)^{1/n−1}。G 的纯标量因子在椭圆解 flux 归一化中消去，
    故只有 H 依赖进入算子形状。n=1 时因子恒为 1（逐位退化为牛顿）。

    ⚠️ 适用域：一阶近似，不得引 B&F25 作为方法学依据（见设计规格 §1.5）。
    """

    def __init__(self, n: float, base: "ClosureProvider | None" = None) -> None:
        self._n = float(n)
        self._base = base if base is not None else NewtonianClosure()

    def mobility(self, c_bar, m: float, eta1: float, eta2: float, H) -> Array:
        I1 = np.asarray(self._base.mobility(c_bar, m, eta1, eta2, H), dtype=float)
        if abs(self._n - 1.0) < 1e-12:
            return I1
        Harr = np.asarray(H, dtype=float)
        Hbar = float(np.mean(Harr))
        return I1 * (Harr / Hbar) ** (1.0 / self._n - 1.0)

    def buoyant_mobility(self, c_bar, m: float, eta1: float, eta2: float, H) -> Array:
        return self._base.buoyant_mobility(c_bar, m, eta1, eta2, H)
```

`annulus_d2dga.py`：
1. `__init__` 形参（Task 2 新参之后）加：
```python
        enable_power_law_gap_correction: bool = False,  # B-3 opt-in：幂律间隙一阶修正（默认关）
```
存属性 `self.enable_power_law_gap_correction = enable_power_law_gap_correction`。
2. `_velocity_stream_function` 内，在算完 `m_ratio`（`:1292`）之后、`solve_stream_function` 调用（`:1338`）之前构造闭包，并把该调用改为：
```python
        # B-3：代表幂律指数（顶替液/水泥相）；非幂律/HB 时 n=1（无修正）
        if self.enable_power_law_gap_correction:
            n_rep = 1.0
            if cement_fluid is not None and getattr(cement_fluid, "power_law_n", None):
                n_rep = float(cement_fluid.power_law_n)
            from cemdisp.models2d.hb_closure import PowerLawGapClosure
            closure = PowerLawGapClosure(n_rep)
        else:
            closure = None
        psi = solve_stream_function(geom, c_bar, eta1, eta2, m_ratio, b_field,
                                    ny=ny, nz=nz, closure=closure, wall=wall)
```

> 注：Task 2 已把 `wall=wall` 加进该调用；Task 3 只再加 `closure=closure`。

- [ ] **Step 4: 跑测试确认通过 + 全量回归**

Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest tests/contract/test_power_law_gap_closure.py -v`
Expected: PASS（2 passed）

Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest tests/ -v`
Expected: 无新增失败。

- [ ] **Step 5: 提交**

```bash
git add cemdisp/models2d/hb_closure.py cemdisp/models2d/annulus_d2dga.py tests/contract/test_power_law_gap_closure.py
git commit -m "feat(closure): B-3 幂律间隙一阶修正（opt-in，默认关；n=1 逐位退化）"
```

---

### Task 4: 牛顿极限回归护栏 + 八井 A/B 量化

**Files:**
- Test: `tests/contract/test_hb_newtonian_limit.py`
- Create: `scripts/entrypoints/run_hb_ablation.py`
- Create: `results/流变口径A-B_2026-09-16/`（脚本产出 CSV）
- Modify: `docs/源模型口径与适用域声明.md`（增补声明）

**Interfaces:**
- Consumes: Task 2/3 的两个开关；`cemdisp.data.loaders`（8 井 loader）；`cemdisp.models2d.AnnulusD2DGASolver`。
- Produces: 对照 CSV `ablation_8wells.csv`（列：`well, variant, eta_E, eta_N, b_number`）。

- [ ] **Step 1: 写测试（牛顿极限护栏）**

```python
# tests/contract/test_hb_newtonian_limit.py
"""B-4：L1 硬约束——牛顿极限（τ_y=0、n=1）下新开关不改变结果。"""
import numpy as np

from cemdisp.models2d.stream_function import solve_stream_function


def _geom_arrays(ny=9, nz=4):
    return {"phi": np.linspace(0.0, 1.0, ny), "s": np.linspace(0.0, 30.0, nz),
            "H": np.full((ny, nz), 0.01), "hole_mm": np.full(nz, 215.9),
            "od_mm": np.full(nz, 168.3), "y": np.linspace(0.0, np.pi * 0.1, ny),
            "b": np.full((ny, nz), 0.02), "inc_deg": np.zeros(nz)}


def test_zero_yield_wall_is_bitwise_unchanged():
    """τ_y=0 ⇒ wall 恒 0 ⇒ 与不传 wall 逐位一致。"""
    g = _geom_arrays(); c = np.full((9, 4), 0.5); b = np.zeros((2, 9, 4))
    base = solve_stream_function(g, c, 0.058, 0.171, 0.34, b)
    gated = solve_stream_function(g, c, 0.058, 0.171, 0.34, b, wall=np.zeros((9, 4)))
    assert np.array_equal(base, gated)
```

- [ ] **Step 2: 跑测试确认通过**

Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest tests/contract/test_hb_newtonian_limit.py -v`
Expected: PASS（Task 2 完成后）。若 FAIL，回到 Task 2 修 `conductance` 的 wall=None 分支。

- [ ] **Step 3: 写 A/B 脚本**

```python
# scripts/entrypoints/run_hb_ablation.py
"""阶段 B 八井 A/B：量化流变口径对 η_E/η_N 的净影响（填补 R34 缺口）。

变体：V0 默认 / V1 屈服门 / V2 幂律修正 / V3 两者 / V4 旧代数路径。
用法：PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/entrypoints/run_hb_ablation.py
输出：results/流变口径A-B_2026-09-16/ablation_8wells.csv
"""
from __future__ import annotations

import csv
from pathlib import Path

from cemdisp.data import loaders
from cemdisp.models2d import AnnulusD2DGASolver

WELL_LOADERS = {n[5:-8]: getattr(loaders, n) for n in dir(loaders)
                if n.startswith("load_") and n.endswith("_tailpipe")}
VARIANTS = {
    "V0_baseline": {},
    "V1_yield_gate": {"enable_stream_yield_gate": True},
    "V2_power_law": {"enable_power_law_gap_correction": True},
    "V3_both": {"enable_stream_yield_gate": True, "enable_power_law_gap_correction": True},
    "V4_legacy_algebraic": {"enable_stream_function": False},
}


def main() -> None:
    out = Path("results/流变口径A-B_2026-09-16")
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for wname, loader in sorted(WELL_LOADERS.items()):
        well, fluids, schedule, provider = loader()
        for vname, kw in VARIANTS.items():
            solver = AnnulusD2DGASolver(total_t=3600.0, nz=60, **kw)  # 粗网格快速对照
            res = solver.run(well, fluids, provider, schedule=schedule)
            fin = res.summary["最终结果"]
            row = {"well": wname, "variant": vname,
                   "eta_E": fin["全井段最终有效顶替效率"],
                   "eta_N": fin["窄四分位效率"], "b_number": fin["浮力数_b"]}
            rows.append(row)
            print(f"{wname:8s} {vname:22s} eta_E={row['eta_E']:.4f} eta_N={row['eta_N']:.4f}")
    with (out / "ablation_8wells.csv").open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["well", "variant", "eta_E", "eta_N", "b_number"])
        w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
```

> ⚠️ **实施前核对**：先跑
> `python -c "import cemdisp.data.loaders as m; print([n for n in dir(m) if n.startswith('load_')])"`
> 确认 loader 函数名与返回签名；不符即按实际导出名修正 `WELL_LOADERS` 构造。

- [ ] **Step 4: 跑脚本产出对照表**

Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/entrypoints/run_hb_ablation.py`
Expected: 打印 8×5=40 行；生成 `ablation_8wells.csv`。**V0 与 V4 的差异**即"流变口径 + 路径"的联合效应；**V1/V2 相对 V0 的 Δη** 即两开关的各自效应。

- [ ] **Step 5: 增补适用域声明**

在 `docs/源模型口径与适用域声明.md` §1 声明 3 之后追加：
```markdown
### 声明 3b：阶段 B 的幂律间隙修正是一阶近似（非 B&F25 口径）

`enable_power_law_gap_correction` 用 `I₁·(H/H̄)^{1/n−1}` 近似幂律的间隙标度
（依据：单流体幂律槽流 ū∝H^{1+1/n}G^{1/n}；B&F25 (2.14)）。它**不是** B&F25 的
数值闭包，**不得**在论文中引 B&F25 作为其方法学依据。完整 HB 闭包见阶段 A。
```

- [ ] **Step 6: 提交**

```bash
git add tests/contract/test_hb_newtonian_limit.py scripts/entrypoints/run_hb_ablation.py "results/流变口径A-B_2026-09-16" docs/源模型口径与适用域声明.md
git commit -m "feat(closure): B-4 牛顿极限护栏 + 八井 A/B 量化 + 适用域声明增补"
```

---

## Self-Review 记录

- **Spec 覆盖**：B-1→Task1、B-2→Task2、B-3→Task3、B-4→Task4；规格 §8 文件清单齐；§6 的 L1 锚落在 Task4 Step1 + Task2 Step1。
- **占位符扫描**：无 TBD/TODO；所有代码步均给出可运行代码。
- **类型一致性**：`ClosureProvider.mobility(c_bar, m, eta1, eta2, H)` 在 Task1 定义、Task3 沿用；`_WALL_CONDUCTANCE_FLOOR` Task2 定义并被测试；`enable_stream_yield_gate`/`enable_power_law_gap_correction` 命名全程一致。
- **偏离规格项（需用户知悉）**：两个新开关默认 `False`（规格暗示 True）——为保 8 井权威结果与既有测试逐位不变；翻转默认属阶段 B 验收后的独立决策。
