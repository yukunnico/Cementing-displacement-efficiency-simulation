# 环空顶替模型失真修正 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 恢复环空段套管居中度对顶替效率的真实（阈值型）敏感性，让主指标能看见窄边窜槽（η_N），同时修复 CFL 自适应步长下的弥散放大、b³ 分流高 Re 过陡、I3 浮力通量近乎无效等失真。

**Architecture:** 在 `AnnulusD2DGASolver`（`cemdisp/models2d/annulus_d2dga.py`）内以独立 `__init__` 开关叠加 6 个机制，新机制默认关，关闭时逐位复现 commit `4893a86` 基线。`_compute_velocity` 返回值由 8 元组固定扩为 12 元组（新增 tau_y/eta2/n_mix/kappa_mix 四场，始终计算），统一服务 M2/M3/I3。M2 的 Maleki 摩擦闭包抽到独立纯函数模块 `cemdisp/models2d/regime_closure.py`。

**Tech Stack:** Python 3、NumPy（全场向量化，无 Python 网格循环）、pytest。conda 环境 `shenjingwangluo`。

**Spec:** `docs/superpowers/specs/2026-08-23-annulus-distortion-fix-design.md`（含 2026-08-23 代码核查修订 11 条；本计划基于该修订后的规格）。执行者必须同时读规格与本计划。

## Global Constraints

- **基线逐位复现**：所有新开关默认关闭时，summary/depth_profiles/wall/metrics 数值必须与 commit `4893a86` 逐位一致（浮点误差内）。这是归因可信的前提，每个任务的测试都要验证。
- **固定 12 元组**：`_compute_velocity` 返回长度固定为 12，不用条件长度；`_compute_props` 由 7 扩为 9。新增场始终计算（加权运算轻），开关只控制是否参与 pref/通量。
- **代码简洁**：不新增无依据的常数；保留关键物理注释；不做无关重构。
- **中文注释、UTF-8**：运行脚本前设 `PYTHONIOENCODING=utf-8 PYTHONUTF8=1`。
- **测试命令**：`PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/ -v`（或单文件）。
- **提交粒度**：每个 Task 一个 commit，信息格式 `feat(annulus): <Mx 描述>（开关=xx 默认关）`，附验收结果。
- **不删生产代码**：legacy 脚本归档另行处理，不在本计划范围。

### 文件结构（本计划产出）

| 文件 | 责任 | 动作 |
|------|------|------|
| `cemdisp/models2d/annulus_d2dga.py` | 核心求解器：M0 指标、M1 弥散、M3 屈服门槛、M2 编排、M4 e护栏、I3 局部化 | 修改 |
| `cemdisp/models2d/regime_closure.py` | M2 纯函数：Metzner-Reed Re、Hedström、层流/湍流/过渡摩擦因子、阻力权重 | 新建 |
| `cemdisp/models2d/d2dga_flux.py` | I3 浮力通量：`max`→`np.maximum` 修数组崩溃 + 类型注解 | 修改（1行+注解） |
| `tests/test_m0_metrics.py` | M0 指标层 solver 级测试 | 新建 |
| `tests/test_regime_closure.py` | M2 摩擦闭包纯函数单测 | 新建 |
| `tests/test_improved_d2dga_annulus.py` | M3 门开新用例、扩元解包同步、`_compute_props` 7→9 断言更新 | 修改 |
| `tests/test_yield_deadzone.py` | `_compute_velocity` 8→12 元组解包同步 | 修改 |
| `tests/test_d2dga_flux.py` | I3 数组 eta2 回归测试 | 修改 |
| `scripts/dispersion_scale_sensitivity_scan.py` | M1 κ 扫描（scale∈{0,0.25,0.5,1.0}） | 新建 |
| `docs/superpowers/plans/2026-08-23-annulus-distortion-fix.md` | 本计划 | 新建 |

### 实施顺序（阶段）

阶段0 M0 → 阶段1 M1 → 阶段2 M3 → 阶段3 M2 → 阶段4 M4 → 阶段5 I3 → 阶段6 附加B（可选）。
**硬约束：M4（Task 10）不得在 M3（Task 8-9）就位前开启 e_clip_max>0.55**，否则重演 A4 的 529:1 窜槽拉爆。

---

## Task 1: 扩元基础设施（_compute_props 7→9，_compute_velocity 8→12，3 处解包同步）

**Files:**
- Modify: `cemdisp/models2d/annulus_d2dga.py`（`_compute_props` 返回 :500；`_compute_velocity` 返回 :673；解包 :805、:933）
- Test: `tests/test_improved_d2dga_annulus.py:56`（`len(out)==7`）
- Test: `tests/test_yield_deadzone.py:209`（精确 8 元组解包）；同文件 `test_compute_props_returns_seven_values`（约 :223-250，断言 7 元组）

**Interfaces:**
- Produces: `_compute_props(...) -> Tuple[mu, rho, mud, tau_y, m_field, eta1, eta2, n_mix, kappa_mix]`（9 元组，新增 `n_mix`/`kappa_mix`，形状均 `(ny,nz)`）。
- Produces: `_compute_velocity(...) -> Tuple[w, v, mu_reg, rho, mud, Re, mu_turbulent, m_field, tau_y, eta2, n_mix, kappa_mix]`（12 元组）。
- Consumes: 后续 Task 2/5/8/9/11 直接从解包取 `tau_y`/`eta2`/`n_mix`/`kappa_mix`，不再各自重算。

- [ ] **Step 1: 先加会失败的断言（扩元长度契约）**

在 `tests/test_improved_d2dga_annulus.py` 的 `TestMFieldFromProps` 类里，把现有 `test_m_field_returned_and_shape`（约 :41-62）末尾的 `assert len(out) == 7` 改为：

```python
        # Task1: _compute_props 扩为 9 元组（+ n_mix/kappa_mix）
        assert len(out) == 9
        n_mix, kappa_mix = out[7], out[8]
        assert n_mix.shape == (ny, nz)
        assert kappa_mix.shape == (ny, nz)
        assert np.all(n_mix > 0.0) and np.all(kappa_mix > 0.0)
```

再在该文件末尾新增一个测试类验证 `_compute_velocity` 12 元组：

```python
class TestComputeVelocityTuple12:
    def test_returns_twelve_tuple(self):
        from cemdisp.data.fluid_spec import FluidSpec, FluidRole, RheologyModel
        s = _make_solver()
        well = _toy_well()
        geom = s._build_geom(well)
        ny, nz = s.ny, s.nz
        lead = np.full((ny, nz), 0.3); tail = np.zeros((ny, nz))
        spacer = np.zeros((ny, nz)); flusher = np.zeros((ny, nz))
        w_prev = np.full((ny, nz), 0.3); wall = np.zeros((ny, nz))
        mud_f = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                          rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
        lead_f = FluidSpec(name="lead", role=FluidRole.LEAD, density_kg_m3=1900.0,
                           rheology_model=RheologyModel.POWER_LAW, power_law_n=0.7, consistency_k=0.4)
        out = s._compute_velocity(well, geom, mud_f, None, lead_f, None, None,
                                  None, 0.01, w_prev, 0.0, wall=wall)
        assert len(out) == 12
        assert out[8].shape == (ny, nz)   # tau_y
        assert out[9].shape == (ny, nz)   # eta2
        assert out[10].shape == (ny, nz)  # n_mix
        assert out[11].shape == (ny, nz)  # kappa_mix
```

> 注意：`_compute_velocity` 的精确形参以当前签名为准（读 :537-551）。若实参顺序与上例不符，按当前签名调整，保持"传最小合法输入"即可。

- [ ] **Step 2: 运行测试确认失败**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/test_improved_d2dga_annulus.py::TestMFieldFromProps tests/test_improved_d2dga_annulus.py::TestComputeVelocityTuple12 -v
```
Expected: FAIL（当前 7 元组/8 元组，长度不匹配）。

- [ ] **Step 3: 实现 _phase_power_law_params 助手 + n_mix/kappa_mix 加权**

在 `annulus_d2dga.py` 中给 `AnnulusD2DGASolver` 加一个静态方法（放在 `_apparent_viscosity` 附近，约 :408 之后）。先确认文件已 `from cemdisp.data.fluid_spec import RheologyModel`（若无则加）：

```python
    @staticmethod
    def _phase_power_law_params(fluid) -> tuple[float, float]:
        """把任意流变模型映射为 (幂律指数 n, 稠度 K[Pa·s^n])，供 M2 混合 n/k 加权。

        NEWTONIAN/BINGHAM -> n=1, K=plastic_viscosity_pa_s；
        POWER_LAW/HERSCHEL_BULKLEY -> (power_law_n, consistency_k)。
        HB 的屈服应力由 tau_y 场单独携带，不在这里折进 K。
        """
        if fluid is None:
            return 1.0, 1.0e-6
        rm = fluid.rheology_model
        if rm == RheologyModel.POWER_LAW or rm == RheologyModel.HERSCHEL_BULKLEY:
            return float(fluid.power_law_n), float(fluid.consistency_k)
        return 1.0, float(fluid.plastic_viscosity_pa_s)
```

在 `_compute_props`（:442-500）里，与 mu/tau_y 加权同模式（:461-482）新增 n/k 加权。在四相 mu 累加块之后、`return`（:500）之前插入：

```python
        # Task1/M2: 混合物幂律参数（体积分数加权 n，对数加权 K），供 Metzner-Reed Re/He
        n_mud, k_mud = self._phase_power_law_params(mud_fluid)
        n_lead, k_lead = self._phase_power_law_params(lead_fluid)
        n_tail, k_tail = self._phase_power_law_params(tail_fluid)
        n_sp, k_sp = self._phase_power_law_params(spacer_fluid)
        n_mix = mud * n_mud + lead * n_lead + tail * n_tail + spacer * n_sp
        log_k_mix = (mud * np.log(max(k_mud, 1e-12)) + lead * np.log(max(k_lead, 1e-12))
                     + tail * np.log(max(k_tail, 1e-12)) + spacer * np.log(max(k_sp, 1e-12)))
        kappa_mix = np.exp(log_k_mix)
```

把 `return`（:500）改为 9 元组：

```python
        return mu, rho, mud, tau_y, m_field, eta1, eta2, n_mix, kappa_mix
```

- [ ] **Step 4: _compute_velocity 解包 9 元组并返回 12 元组**

在 `_compute_velocity` 内，把 `_compute_props` 的解包（约 :577，现为 7 个变量）改为接收 n_mix/kappa_mix：

```python
        mu, rho, mud, tau_y, m_field, eta1, eta2, n_mix, kappa_mix = self._compute_props(...)
```

把返回语句（:673）改为：

```python
        return w, v, mu_reg, rho, mud, Re, mu_turbulent, m_field, tau_y, eta2, n_mix, kappa_mix
```

同步该函数返回类型注解（约 :551）为 12 个 `Array`。

- [ ] **Step 5: 同步 3 处解包点**

`annulus_d2dga.py` 内两处（:805 泵注、:933 停泵），把现 8 变量解包末尾加 `, tau_y, eta2, n_mix, kappa_mix`。变量名若与局部已存变量冲突，解包时用 `_tau_y, _eta2, _n_mix, _kappa_mix`（这些场在 Task 8/11 才消费；Task 1 只保证解包不炸）。

`tests/test_yield_deadzone.py:209` 把：
```python
        w, v, mu_reg, rho, mud_frac, Re, mu_turbulent, m_field = solver._compute_velocity(
```
改为：
```python
        w, v, mu_reg, rho, mud_frac, Re, mu_turbulent, m_field, _tau_y, _eta2, _n_mix, _kappa_mix = solver._compute_velocity(
```

同文件 `test_compute_props_returns_seven_values`（约 :223-250）：把 `len(result) == 7` 改 `== 9`，解包语句补 `, n_mix, kappa_mix` 两个占位变量，并加 `assert n_mix.shape == (...)`（按该测试现有形状断言风格）。

- [ ] **Step 6: 运行全量测试确认零回归**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/ -q
```
Expected: 全部 PASS（基线数值不变，因为新场未被任何 pref/通量消费）。

- [ ] **Step 7: Commit**

```bash
git add cemdisp/models2d/annulus_d2dga.py tests/test_improved_d2dga_annulus.py tests/test_yield_deadzone.py
git commit -m "feat(annulus): _compute_props/_compute_velocity 扩元(9/12元组)供M2/M3/I3复用"
```

---

## Task 2: M0 — η_N 窄四分位效率进 summary

**Files:**
- Modify: `cemdisp/models2d/annulus_d2dga.py`（顶部 import；summary 构建 :1079-1096）
- Test: `tests/test_m0_metrics.py`（新建）

**Interfaces:**
- Consumes: `cemdisp.diagnostics.displacement_metrics._narrow_quarter_efficiency(cement: ndarray, geom: dict) -> float`、`_trapez2d(field, geom) -> float`（两者均为模块级，无循环依赖；`displacement_metrics` 运行时仅 import numpy）。
- Produces: `summary["最终结果"]["窄四分位效率"]: float`、`summary["eta_narrow"]: float`。

- [ ] **Step 1: 写失败测试**

新建 `tests/test_m0_metrics.py`：

```python
"""M0 指标层 solver 级测试：η_N 进 summary、评价窗、失稳指数、低尾指标。"""
import numpy as np
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
from cemdisp.data.well_spec import WellSpec, DepthValuePoint, EvaluationWindow


def _solver(**kw):
    return AnnulusD2DGASolver(dt=4.0, nz=30, ny=12, total_t=40.0, **kw)


def _multi_window_well():
    pts = lambda d, v: DepthValuePoint(depth_md_m=d, value=v)
    return WellSpec(
        well_name="m0", top_md_m=1000.0, bottom_md_m=1100.0,
        shoe_md_m=1100.0, hanger_md_m=1000.0,
        casing_id_mm=200.0, liner_od_mm=139.7, liner_id_mm=108.0,
        hole_diameter_profile=[pts(1000.0, 215.9), pts(1100.0, 215.9)],
        inclination_profile=[pts(1000.0, 5.0), pts(1100.0, 5.0)],
        standoff_profile=[pts(1000.0, 0.83), pts(1100.0, 0.83)],
        evaluation_windows=[
            EvaluationWindow(name="cbl窗", top_md_m=1020.0, bottom_md_m=1080.0, window_type="cbl"),
            EvaluationWindow(name="model窗", top_md_m=1000.0, bottom_md_m=1100.0, window_type="model_focus"),
        ],
    )


def _run_minimal():
    """跑一个最小 solver；run() 实参以当前签名为准，参照 tests/test_six_well_integration.py 装配。"""
    s = _solver()
    return s.run(_multi_window_well(), fluids=None, inlet_state_provider=None)


def test_eta_narrow_in_summary_and_in_range():
    from cemdisp.diagnostics.displacement_metrics import compute_displacement_metrics
    res = _run_minimal()
    fr = res.summary["最终结果"]
    assert "窄四分位效率" in fr
    assert 0.0 <= fr["窄四分位效率"] <= 1.0
    assert 0.0 <= res.summary["eta_narrow"] <= 1.0
    dm = compute_displacement_metrics(res)
    assert abs(fr["窄四分位效率"] - dm.eta_narrow) < 1e-9
```

> 注意：`run()` 的实参以当前签名为准。若 `fluids=None`/`inlet_state_provider=None` 不被接受，参照 `tests/test_six_well_integration.py` 或现有 solver 测试的最小流体/入口装配；测试目标是断言 summary 键。

- [ ] **Step 2: 运行确认失败**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/test_m0_metrics.py::test_eta_narrow_in_summary_and_in_range -v
```
Expected: FAIL（KeyError: '窄四分位效率'）。

- [ ] **Step 3: 实现**

顶部 import 区（:34-39 附近）加（先 grep `_trapez2d` 确认 solver 内是否已有同名模块级函数，有则只 import `_narrow_quarter_efficiency`，避免重复）：
```python
from cemdisp.diagnostics.displacement_metrics import _narrow_quarter_efficiency
```

在 summary 字典的 `"最终结果"` 块（:1084-1091）内，与 `"全井段最终有效顶替效率"` 并列加：
```python
            "窄四分位效率": _narrow_quarter_efficiency(cement, geom),
```
在顶层英文别名区（:1092-1095）加：
```python
        "eta_narrow": _narrow_quarter_efficiency(cement, geom),
```
`cement`（:1053）、`geom`（:744）在该作用域均为局部变量。

- [ ] **Step 4: 运行测试通过**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/test_m0_metrics.py -v
```
Expected: PASS。

- [ ] **Step 5: Commit**

```bash
git add cemdisp/models2d/annulus_d2dga.py tests/test_m0_metrics.py
git commit -m "feat(annulus): M0 窄四分位效率η_N进summary（零物理扰动）"
```

---

## Task 3: M0 — 失稳指数去饱和 + 评价窗效率 + 低尾指标

**Files:**
- Modify: `cemdisp/models2d/annulus_d2dga.py`（summary 构建 :1079-1096）
- Test: `tests/test_m0_metrics.py`（追加）

**Interfaces:**
- Produces:
  - `summary["最终结果"]["最终失稳指数_线性"]: float`（= metrics 末行 `instability_proxy`）
  - `summary["最终结果"]["最终失稳指数_对数"]: float`（= log10(1+proxy)）
  - `summary["评价窗效率"]: dict[str, dict]`，每窗 `{window_type, eta_E, eta_N}`
  - `summary["低尾指标"]: {"standoff低于0.5段占比": float, "窄边效率低于0.05域占比": float}`

- [ ] **Step 1: 写失败测试**

在 `tests/test_m0_metrics.py` 追加：

```python
def test_instability_index_linear_and_log():
    res = _run_minimal()
    fr = res.summary["最终结果"]
    assert "最终失稳指数_线性" in fr and "最终失稳指数_对数" in fr
    proxy = float(res.metrics["instability_proxy"].iloc[-1])
    assert abs(fr["最终失稳指数_线性"] - proxy) < 1e-9
    assert abs(fr["最终失稳指数_对数"] - np.log10(1.0 + proxy)) < 1e-9


def test_evaluation_window_efficiencies():
    res = _run_minimal()
    we = res.summary["评价窗效率"]
    assert {"cbl窗", "model窗"} <= set(we.keys())
    for name, d in we.items():
        assert d["window_type"] in ("cbl", "model_focus")
        assert 0.0 <= d["eta_E"] <= 1.0
        assert 0.0 <= d["eta_N"] <= 1.0


def test_low_tail_indicators():
    res = _run_minimal()
    lt = res.summary["低尾指标"]
    assert 0.0 <= lt["standoff低于0.5段占比"] <= 1.0
    assert 0.0 <= lt["窄边效率低于0.05域占比"] <= 1.0


def test_m0_does_not_change_existing_keys():
    res = _run_minimal()
    fr = res.summary["最终结果"]
    for k in ["全井段最终有效顶替效率", "最终水泥浆占据率", "最终窜槽指数", "最终混浆指数", "最终失稳指数"]:
        assert k in fr
    assert "effective_efficiency" in res.summary
```

- [ ] **Step 2: 运行确认失败**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/test_m0_metrics.py -k "instability or window or low_tail or existing_keys" -v
```
Expected: FAIL（缺键）。

- [ ] **Step 3: 实现失稳指数两键**

在 summary 构建之前（`final = metrics.iloc[-1]` 之后，约 :1055），加：
```python
    _inst_lin = float(final["instability_proxy"])
    _inst_log = float(np.log10(1.0 + _inst_lin))
```
在 `"最终结果"` 字典内（与 `"最终失稳指数"` 并列）加：
```python
            "最终失稳指数_线性": _inst_lin,
            "最终失稳指数_对数": _inst_log,
```

- [ ] **Step 4: 实现评价窗效率（模块级纯函数）**

在 `annulus_d2dga.py` 模块级（类外，`_trapez2d` 附近）定义：
```python
def _evaluation_window_efficiencies(well_spec, geom, cement):
    """对每个 EvaluationWindow 做 b 加权 2D 积分，返回 {窗名: {window_type, eta_E, eta_N}}。"""
    out = {}
    md = geom["md"]            # (nz,)，md = bottom - s
    b_full = geom["b"]         # (ny,nz)
    s_full = geom["s"]         # (nz,)
    for w in well_spec.evaluation_windows:
        mask = (md >= w.top_md_m) & (md <= w.bottom_md_m)
        if not bool(mask.any()):
            continue
        b_win = b_full[:, mask]
        c_win = cement[:, mask]
        geom_win = {**geom, "b": b_win, "s": s_full[mask]}
        denom = _trapez2d(b_win, geom_win)
        eta_e = float(_trapez2d(b_win * c_win, geom_win) / max(denom, 1e-12))
        eta_n = float(_narrow_quarter_efficiency(c_win, geom_win))
        out[w.name] = {"window_type": w.window_type, "eta_E": eta_e, "eta_N": eta_n}
    return out
```
> 关键：`_narrow_quarter_efficiency` 按行取最后 ny//4（窄边），窗口掩码只切列，传子网格正确；但必须传切片后的 `geom_win`，否则 b 形状与 c 不匹配会崩。若 solver 内 `_trapez2d` 名字不同（grep 确认，现有有效率积分在 :976-977），用现有的。

在 summary 字典内（`"最终结果"` 块之后、顶层别名之前）加：
```python
        "评价窗效率": _evaluation_window_efficiencies(well_spec, geom, cement),
```

- [ ] **Step 5: 实现低尾指标（模块级纯函数）**

```python
def _low_tail_indicators(geom, cement, ny):
    b_full = geom["b"]
    so_frac = float(np.mean(geom["standoff"] < 0.5))
    n_q = max(1, ny // 4)
    b_q = b_full[-n_q:, :]
    low = (cement[-n_q:, :] < 0.05).astype(float)
    geom_q = {**geom, "b": b_q}
    tail_frac = float(_trapez2d(b_q * low, geom_q) / max(_trapez2d(b_q, geom_q), 1e-12))
    return {"standoff低于0.5段占比": so_frac, "窄边效率低于0.05域占比": tail_frac}
```
summary 字典内加：
```python
        "低尾指标": _low_tail_indicators(geom, cement, self.ny),
```

- [ ] **Step 6: 全量测试 + 零回归确认**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/ -q
```
Expected: 全部 PASS。既有键数值不变（只加键）。

- [ ] **Step 7: 顺手删死 wall 形参（可选但推荐）**

`_depth_profiles`（约 :705）的 `wall` 形参函数体未用。grep 全仓 `_depth_profiles` 确认无外部位置调用后，删除该形参与调用处（约 :1054）的 `wall` 实参。跑全量测试。若有外部位置传参，保留形参仅加 `# 死参数，保留兼容` 注释。

- [ ] **Step 8: Commit**

```bash
git add cemdisp/models2d/annulus_d2dga.py tests/test_m0_metrics.py
git commit -m "feat(annulus): M0 失稳指数去饱和+评价窗效率+低尾指标（只加键零物理扰动）"
```

---

## Task 4: M1 — 弥散系数按 dt 归一（恢复量纲正确性）

**Files:**
- Modify: `cemdisp/models2d/annulus_d2dga.py`（`__init__` :188-259 加参数；弥散调用 :860-863）
- Test: `tests/test_m1_dispersion.py`（新建）

**Interfaces:**
- Produces（`__init__` 新参数，默认值=基线）：
  - `dispersion_axial: float = 0.018`
  - `dispersion_azimuthal: float = 0.015`
  - `dispersion_dt_ref: float = 4.0`（名义/唯一运行 dt，系数 fa44ace 引入时即 dt=4.0）
  - `dispersion_dt_scale: float = 1.0`（=1.0 时固定 dt 模式逐位复现基线；CFL 下用 dt_step/dt_ref 归一）

- [ ] **Step 1: 写失败测试**

新建 `tests/test_m1_dispersion.py`：

```python
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver


def test_dispersion_defaults_match_baseline_constants():
    s = AnnulusD2DGASolver(dt=4.0, nz=20, ny=10, total_t=40.0)
    assert s.dispersion_axial == 0.018
    assert s.dispersion_azimuthal == 0.015
    assert s.dispersion_dt_ref == 4.0
    assert s.dispersion_dt_scale == 1.0


def test_dispersion_eff_scales_with_dt():
    s = AnnulusD2DGASolver(dt=4.0, nz=20, ny=10, total_t=40.0)
    # dt_step=dt_ref 且 scale=1 时有效系数=原系数；dt 缩小则同比缩小
    assert abs(s.dispersion_axial * (4.0 / s.dispersion_dt_ref) - 0.018) < 1e-12
    assert abs(s.dispersion_axial * (0.118 / 4.0) - 0.018 * 0.0295) < 1e-6
```

- [ ] **Step 2: 运行确认失败**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/test_m1_dispersion.py -v
```
Expected: FAIL（属性不存在）。

- [ ] **Step 3: __init__ 加参数并存 self**

在 `__init__` 签名（dt 相关参数附近，:188-208）加上述 4 个参数及默认值。在 `self.` 赋值区（:259 附近）加：
```python
        self.dispersion_axial = dispersion_axial
        self.dispersion_azimuthal = dispersion_azimuthal
        self.dispersion_dt_ref = dispersion_dt_ref
        self.dispersion_dt_scale = dispersion_dt_scale
```

- [ ] **Step 4: 改弥散调用点（:860-863）**

把硬编码：
```python
            lead = self._smooth_dispersion(lead, axial=0.018, azimuthal=0.015)
            tail = self._smooth_dispersion(tail, axial=0.018, azimuthal=0.015)
            spacer = self._smooth_dispersion(spacer, axial=0.012, azimuthal=0.012)
            flusher = self._smooth_dispersion(flusher, axial=0.012, azimuthal=0.012)
```
改为（dt_step 在 :822-825 已算好；核实该作用域变量名）：
```python
            _dt_norm = self.dispersion_dt_scale * (dt_step / self.dispersion_dt_ref)
            _ax = self.dispersion_axial * _dt_norm
            _az = self.dispersion_azimuthal * _dt_norm
            # spacer/flusher 保持与 lead 的相对比 0.012/0.018=0.667、0.012/0.015=0.8
            _ax_sf = self.dispersion_axial * 0.667 * _dt_norm
            _az_sf = self.dispersion_azimuthal * 0.8 * _dt_norm
            lead = self._smooth_dispersion(lead, axial=_ax, azimuthal=_az)
            tail = self._smooth_dispersion(tail, axial=_ax, azimuthal=_az)
            spacer = self._smooth_dispersion(spacer, axial=_ax_sf, azimuthal=_az_sf)
            flusher = self._smooth_dispersion(flusher, axial=_ax_sf, azimuthal=_az_sf)
```
> 关键：固定 dt 模式 dt_step==self.dt==4.0 且 scale==1.0 → 系数=0.018/0.015/0.012，与基线逐位一致。

- [ ] **Step 5: 全量测试确认零回归**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/ -q
```
Expected: 全部 PASS。

- [ ] **Step 6: Commit**

```bash
git add cemdisp/models2d/annulus_d2dga.py tests/test_m1_dispersion.py
git commit -m "feat(annulus): M1 弥散系数按dt归一（开关默认1=逐位复现基线）"
```

---

## Task 5: M1 — κ 扫描验收脚本 + CFL on/off 存档

**Files:**
- Create: `scripts/dispersion_scale_sensitivity_scan.py`

**Interfaces:**
- 复用 `scripts/density_contrast_sensitivity_scan.py:136-169` 的手动流水线模式（1D casing → build_coupled_annulus_inlet_provider → AnnulusD2DGASolver → run → 指标 → CSV），因为 `run_one_level` 不透传弥散参数。

- [ ] **Step 1: 写脚本**

照抄 density 脚本骨架，常量改为：
```python
SCALE_VALUES = [0.0, 0.25, 0.5, 1.0]
NZ, NY = 500, 40
DT = 4.0
ENABLE_CFL_ADAPTIVE = True   # 归一化只在 dt_step≠4.0 时真正生效
```
循环里构造 solver 时传：
```python
solver = AnnulusD2DGASolver(
    dt=DT, nz=NZ, ny=NY, total_t=annulus_stop_time_s,
    enable_cfl_adaptive=ENABLE_CFL_ADAPTIVE,
    dispersion_dt_scale=scale,
)
```
CSV 列：`scale, cfl_on, effective_efficiency, channeling_index, mixing_index, instability_index, cement_occupation, front_length_m, elapsed_s`。

`front_length_m`：从 `result.depth_profiles` 的"水泥平均浓度"列，插值 c=0.02 与 c=0.98 的深度差；可复用 `cemdisp/diagnostics/flow_classification.py` 等值线位置思路（grep `_contour_position`），拿不到就用 numpy.interp 手写。

每 case 即时追加写 CSV（防崩），末尾打印各 scale 下 mixing/η_E/channeling/front_length 表 + 线性斜率。

- [ ] **Step 2: CFL on/off 存档前置（设计稿方案乙要求）**

先跑 CFL off（`enable_cfl_adaptive=False`，固定 dt=4.0）scale=1.0 基线，再跑 CFL on。两组存同一 CSV（`cfl_on` 列区分）到 `results/m1_dispersion_scale/`。这是 M1 归因基线。

- [ ] **Step 3: 运行脚本（手动，耗时）**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python scripts/dispersion_scale_sensitivity_scan.py
```
验收（设计稿 §4）：①mixing 从 0.59 降到 0.2–0.35；②scale→0 时 mixing 主要由平流数值扩散主导，平台值由 κ 主导而非 dt_step；③前沿长度回到米级（2–5m）。若 mixing 未降到区间，不要调常数凑数——记录实际值回报（可能平流数值扩散占主导，属 Track B TVD 范畴）。

- [ ] **Step 4: Commit**

```bash
git add scripts/dispersion_scale_sensitivity_scan.py
git commit -m "test(scripts): M1 弥散κ扫描脚本(scale∈0/0.25/0.5/1.0)+CFL on/off存档"
```

---

## Task 6: M2 — regime_closure 纯函数模块 + 单测

**Files:**
- Create: `cemdisp/models2d/regime_closure.py`
- Test: `tests/test_regime_closure.py`（新建）

**Interfaces（本 Task 只建纯函数，不接 solver；Task 7 才接）：**
```python
def metzner_reed_re(w, rho, n, kappa, b) -> np.ndarray
def hedstrom_number(tau_y, rho, n, kappa, b) -> np.ndarray
def friction_laminar(re_p, he, n) -> np.ndarray
def friction_dodge_metzner(re_p, n, f0=0.01, tol=1e-4, it_max=15) -> np.ndarray
def friction_transition(re_p, re_crit, re_turb, f_lam_cr, f_turb) -> np.ndarray
def drag_weight(re_p, he, n, re_crit, re_turb_ratio=1.8) -> tuple[np.ndarray, np.ndarray]
```
所有函数数组友好（标量或 ndarray 逐元运算），返回 ndarray。

- [ ] **Step 1: 写失败测试**

新建 `tests/test_regime_closure.py`：

```python
import numpy as np
from cemdisp.models2d import regime_closure as rc


def test_metzner_reed_re_shape_and_positive():
    w = np.full(4, 0.5); rho = np.full(4, 1500.0); n = np.full(4, 0.7)
    kappa = np.full(4, 0.4); b = np.full(4, 0.02)
    re = rc.metzner_reed_re(w, rho, n, kappa, b)
    assert re.shape == (4,) and np.all(re > 0)


def test_laminar_friction_is_24_over_re_when_no_yield():
    re = np.array([10.0, 100.0, 1000.0]); he = np.zeros(3); n = np.ones(3)
    f = rc.friction_laminar(re, he, n)
    np.testing.assert_allclose(f, 24.0 / re, rtol=1e-6)


def test_dodge_metzner_converges_and_monotone():
    re = np.logspace(3, 6, 50); n = np.full(50, 0.7)
    f = rc.friction_dodge_metzner(re, n)
    assert f.shape == (50,) and np.all(f > 0) and np.all(f < 1.0)
    assert np.all(np.diff(f) <= 1e-6)


def test_drag_weight_layer_mask_and_bounds():
    re = np.array([10.0, 500.0, 5000.0, 50000.0])
    he = np.zeros(4); n = np.full(4, 0.8)
    R, mask = rc.drag_weight(re, he, n, re_crit=np.full(4, 2100.0))
    assert R[0] == 1.0
    assert np.all(R >= 0.3) and np.all(R <= 1.5)
    assert mask.shape == (4,)


def test_hedstrom_zero_for_no_yield():
    tau_y = np.zeros(3); rho = np.full(3, 1500.0); n = np.full(3, 0.7)
    kappa = np.full(3, 0.4); b = np.full(3, 0.02)
    assert np.all(rc.hedstrom_number(tau_y, rho, n, kappa, b) == 0.0)
```

- [ ] **Step 2: 运行确认失败**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/test_regime_closure.py -v
```
Expected: FAIL（模块不存在）。

- [ ] **Step 3: 实现 regime_closure.py**

新建 `cemdisp/models2d/regime_closure.py`：

```python
"""M2 局部流态修正的纯函数闭包（Maleki & Frigaard 2017 式58-66）。

所有函数数组友好：标量或 np.ndarray 逐元运算，返回 ndarray。
几何约定：b=effective_b（全间隙，与求解器 γ̇=6|w|/b 自洽）。
"""
import numpy as np


def metzner_reed_re(w, rho, n, kappa, b):
    """式58：Re_p = 6ρw²/(κ γ̇_N^n)，γ̇_N=6|w|/b。"""
    w = np.asarray(w, dtype=float)
    b = np.maximum(np.asarray(b, dtype=float), 1e-9)
    gamma = np.maximum(6.0 * np.abs(w) / b, 1e-6)
    tau_app = np.maximum(np.asarray(kappa, dtype=float) * gamma ** np.asarray(n, dtype=float), 1e-12)
    return 6.0 * np.asarray(rho, dtype=float) * w * w / tau_app


def hedstrom_number(tau_y, rho, n, kappa, b):
    """式59：He = τy [ρ^n b^(2n)/κ²]^(1/(2-n))。"""
    tau_y = np.asarray(tau_y, dtype=float)
    rho = np.asarray(rho, dtype=float)
    n = np.asarray(n, dtype=float)
    kappa = np.maximum(np.asarray(kappa, dtype=float), 1e-12)
    b = np.maximum(np.asarray(b, dtype=float), 1e-9)
    inner = rho ** n * b ** (2.0 * n) / kappa ** 2
    return tau_y * np.power(np.maximum(inner, 1e-300), 1.0 / (2.0 - n))


def friction_laminar(re_p, he, n):
    """式60-61：层流 f=24/Re_p，含塞流核修正 yY=He/Hw。

    无屈服(he=0)时退化为 24/Re_p。塞流核幂次以 Maleki 式60-61 为准，
    实现后用文献算例复核；此处给标准槽流形式 (1-yY)^2。
    """
    re_p = np.maximum(np.asarray(re_p, dtype=float), 1e-9)
    he = np.asarray(he, dtype=float)
    hw = 24.0 / re_p
    yY = np.clip(he / np.maximum(hw, 1e-9), 0.0, 0.95)
    return 24.0 * (1.0 - yY) ** 2 / re_p


def friction_dodge_metzner(re_p, n, f0=0.01, tol=1e-4, it_max=15):
    """式65：1/√f = (4/n^0.75) log10(Re f^(1-n/2)) - 0.4/n^1.2。固定点迭代。"""
    re_p = np.asarray(re_p, dtype=float)
    n = np.asarray(n, dtype=float)
    f = np.full_like(re_p, f0, dtype=float)
    for _ in range(it_max):
        rhs = ((4.0 / n ** 0.75) * np.log10(np.maximum(re_p * f ** (1.0 - n / 2.0), 1e-30))
               - 0.4 / n ** 1.2)
        f_new = 1.0 / np.maximum(rhs, 1e-9) ** 2
        if np.all(np.abs(f_new - f) / np.maximum(f, 1e-12) < tol):
            f = f_new
            break
        f = f_new
    return f


def friction_transition(re_p, re_crit, re_turb, f_lam_cr, f_turb):
    """过渡区 log-Re 空间线性插值（禁硬切换）。"""
    log_re = np.log10(np.maximum(np.asarray(re_p, dtype=float), 1e-9))
    lo = np.log10(np.maximum(np.asarray(re_crit, dtype=float), 1e-9))
    hi = np.log10(np.maximum(np.asarray(re_turb, dtype=float), 1e-9))
    t = np.clip((log_re - lo) / np.maximum(hi - lo, 1e-9), 0.0, 1.0)
    log_f = (np.log10(np.maximum(f_lam_cr, 1e-12))
             + t * (np.log10(np.maximum(f_turb, 1e-12)) - np.log10(np.maximum(f_lam_cr, 1e-12))))
    return 10.0 ** log_f


def drag_weight(re_p, he, n, re_crit, re_turb_ratio=1.8):
    """返回 (R, regime_mask)。R 乘到 pref 上；mask: 0 层流/1 过渡/2 湍流。

    层流 R=1（保持现有 b² 幂律，不降指数）；过渡/湍流 R=clip(f_lam_cr/f_eff,0.3,1.5)。
    总面积归一保证总流量守恒。Re_crit 由调用方按 He 给出（屈服推迟转捩）。
    """
    re_p = np.asarray(re_p, dtype=float)
    he = np.asarray(he, dtype=float)
    n = np.asarray(n, dtype=float)
    re_crit = np.asarray(re_crit, dtype=float)
    re_turb = re_crit * re_turb_ratio
    f_lam_cr = friction_laminar(re_crit, he, n)
    f_turb = friction_dodge_metzner(np.maximum(re_turb, 1e-9), n)
    f_eff = friction_transition(re_p, re_crit, re_turb, f_lam_cr, f_turb)
    R = np.clip(f_lam_cr / np.maximum(f_eff, 1e-9), 0.3, 1.5)
    R = np.where(re_p <= re_crit, 1.0, R)
    mask = np.where(re_p <= re_crit, 0, np.where(re_p >= re_turb, 2, 1)).astype(int)
    return R, mask
```
> `friction_laminar` 的塞流核 yY 是简化联立；Task 7 三重回归（Walton&Bittleston/Z22/Foolad）发现偏差时按 Maleki 式60-61 精确修正，接口不变。

- [ ] **Step 4: 运行单测通过**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/test_regime_closure.py -v
```
Expected: 5 个 PASS。

- [ ] **Step 5: Commit**

```bash
git add cemdisp/models2d/regime_closure.py tests/test_regime_closure.py
git commit -m "feat(annulus): M2 regime_closure纯函数(Metzner-Reed/He/层流湍流摩擦/阻力权重)"
```

---

## Task 7: M2 — 在 _compute_velocity 内接固定点迭代（开关默认关）

**Files:**
- Modify: `cemdisp/models2d/annulus_d2dga.py`（`__init__` 加开关；`:652-660` 包裹固定点循环）
- Test: `tests/test_regime_closure.py`（加 solver 级门开门关对比）

**Interfaces:**
- Produces（`__init__` 新参数）：
  - `enable_regime_split: bool = False`
  - `regime_relax_alpha: float = 0.5`
  - `regime_max_iter: int = 6`
  - `regime_tol_rel: float = 1e-3`
  - `regime_re_turb_ratio: float = 1.8`
- 关键：关闭时 `:652-660` 原代码逐字节不变（零回归）。

- [ ] **Step 1: 写失败测试（门关门开层流极限一致 + 守恒）**

在 `tests/test_regime_closure.py` 追加 solver 级测试（参照 Task 1 `TestComputeVelocityTuple12` 的最小装配，低速使 Re<re_crit 全层流）：

```python
def test_regime_split_off_equals_on_in_laminar_limit():
    import numpy as np
    from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
    from cemdisp.data.fluid_spec import FluidSpec, FluidRole, RheologyModel
    from cemdisp.data.well_spec import WellSpec, DepthValuePoint
    # 用极小排量构造层流极限；门开门关 w 应一致
    ...
    # w_off = solver(enable_regime_split=False); w_on = solver(enable_regime_split=True)
    # np.testing.assert_allclose(w_on, w_off, atol=1e-9)


def test_regime_split_conserves_total_flow():
    # area_weight 归一 => Σw·b·dy·2 == q_half（任意 R 分布）
    ...
```
> 装配复用 `tests/test_improved_d2dga_annulus.py` 里现成的 `_make_solver()`/`_toy_well()` 与流体构造；目标是断言层流极限一致与流量守恒。

- [ ] **Step 2: 运行确认失败**（开关不存在 → AttributeError）

- [ ] **Step 3: __init__ 加 5 参数存 self**（默认值见上）

- [ ] **Step 4: 包裹 :652-660 固定点迭代**

先读当前 :619-660，记下 `dy`、sum 轴/keepdims、`q_half`、`b_mean`、`b_shape`、`i1_base` 的确切变量名与维度。当前结构：
```python
pref = np.maximum(base * b_shape, 1e-8)
pref = pref * (1.0 - wall)
area_weight = 2.0 * np.sum(pref * b * dy, ...)
w = q_half * pref / area_weight
```
改为：
```python
        pref = np.maximum(base * b_shape, 1.0e-8)
        pref *= (1.0 - wall)
        if self.enable_regime_split:
            from cemdisp.models2d import regime_closure as rc
            he = rc.hedstrom_number(tau_y, rho, n_mix, kappa_mix, b)
            re_crit = 2100.0 * (1.0 + 0.1 * he)   # 屈服推迟转捩；标定钮
            w_k = w_prev.copy()
            R = np.ones_like(w_k)
            for _ in range(self.regime_max_iter):
                re_p = rc.metzner_reed_re(w_k, rho, n_mix, kappa_mix, b)
                R_new, _ = rc.drag_weight(re_p, he, n_mix, re_crit, self.regime_re_turb_ratio)
                pref_k = np.maximum(base * b_shape * R_new, 1.0e-8) * (1.0 - wall)
                area_w = 2.0 * np.sum(pref_k * b * dy)  # 轴/keepdims 按现有对齐
                w_raw = q_half * pref_k / np.maximum(area_w, 1e-12)
                w_new = self.regime_relax_alpha * w_raw + (1.0 - self.regime_relax_alpha) * w_k
                if (np.max(np.abs(w_new - w_k))
                        < self.regime_tol_rel * max(np.max(np.abs(w_k)), 1e-12)):
                    w_k = w_new; R = R_new; break
                w_k = w_new; R = R_new
            w = w_k
            area_weight = 2.0 * np.sum(
                np.maximum(base * b_shape * R, 1.0e-8) * (1.0 - wall) * b * dy)
        else:
            area_weight = 2.0 * np.sum(pref * b * dy)   # 原代码逐字复制
            w = q_half * pref / np.maximum(area_weight, 1e-12)
```
> 门关分支必须是原代码逐字复制，保证零回归。浓度相关量（base/b_shape/wall/n_mix/kappa_mix/tau_y/rho/b）在迭代外算好（缓存，黏度保持 w_prev 一步滞后），迭代只重算 Re_p/R/pref/area_weight/w。

- [ ] **Step 5: 全量测试零回归 + 门开层流极限一致**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/ -q
```
Expected: 全 PASS；门关数值与基线逐位一致；门开层流极限与门关一致 <1e-9；总流量守恒。

- [ ] **Step 6: Commit**

```bash
git add cemdisp/models2d/annulus_d2dga.py tests/test_regime_closure.py
git commit -m "feat(annulus): M2 局部流态修正固定点迭代（enable_regime_split默认关，门关逐位复现）"
```

---

## Task 8: M3 — 屈服门槛纯逻辑 + 单测

**Files:**
- Modify: `cemdisp/models2d/annulus_d2dga.py`（新增静态方法 `_yield_gate_wall`；`__init__` 加开关）
- Test: `tests/test_improved_d2dga_annulus.py`（追加门开单测）

**Interfaces:**
- Produces（`__init__` 新参数）：
  - `enable_yield_gate: bool = False`
  - `yield_gate_f_safety: float = 1.15`
  - `yield_gate_c_min_residual: float = 0.01`
- Produces（静态方法，纯函数）：
```python
@staticmethod
def _yield_gate_wall(w, b, mu_reg, tau_y, cement_ever, cement_local,
                     f_safety, c_min_residual) -> np.ndarray  # wall(ny,nz) 0/1
```
- Consumes: Task 1 已让 `_compute_velocity` 返回 `tau_y`；wall 组装留 run() 泵注分支（Task 9），本 Task 只写纯逻辑。

- [ ] **Step 1: 写失败测试**

在 `tests/test_improved_d2dga_annulus.py` 追加：

```python
class TestYieldGateWall:
    def test_reference_row_pick_avoids_frozen_wide_side(self):
        ny, nz = 4, 3
        w = np.array([[0.0, 0.8, 0.0],   # 宽边 row0：第1列冻结、第2列流、第3列冻结
                      [0.0, 0.5, 0.0],
                      [0.0, 0.2, 0.0],
                      [0.0, 0.0, 0.0]])  # 窄边全冻
        b = np.full((ny, nz), 0.02); mu_reg = np.full((ny, nz), 0.1)
        tau_y = np.full((ny, nz), 5.0)
        cement_ever = np.ones((ny, nz)); cement_local = np.full((ny, nz), 0.9)
        wall = AnnulusD2DGASolver._yield_gate_wall(
            w, b, mu_reg, tau_y, cement_ever, cement_local, 1.15, 0.01)
        assert wall.shape == (ny, nz)
        # 第1/3列无流动参考元 -> 整列冻结
        assert wall[0, 0] == 1.0 and wall[0, 2] == 1.0

    def test_residual_wall_keeps_front_gate(self):
        ny, nz = 2, 2
        wall = AnnulusD2DGASolver._yield_gate_wall(
            np.zeros((ny, nz)), np.full((ny, nz), 0.02), np.full((ny, nz), 0.1),
            np.zeros((ny, nz)), np.zeros((ny, nz)), np.zeros((ny, nz)), 1.15, 0.01)
        # cement_ever=0 前锋未到，不得全域冻结
        assert np.all(wall == 0.0)
```

- [ ] **Step 2: 运行确认失败**（静态方法不存在）

- [ ] **Step 3: 实现静态方法**

在 `AnnulusD2DGASolver` 内（`_apparent_viscosity`/`_phase_power_law_params` 附近）加：

```python
    @staticmethod
    def _yield_gate_wall(w, b, mu_reg, tau_y, cement_ever, cement_local,
                         f_safety, c_min_residual):
        """M3 可重启屈服门槛：每深度列以流动最宽元为参考外推壁剪，
        immobile = τw_extrap ≤ f·τy；OR 残泥下限(cement_ever>0 且 cement<c_min_residual)。
        全列无流动 -> 整列冻结。停泵期不调用（run() 泵注分支门控）。"""
        b = np.maximum(b, 1e-12)
        gamma = np.maximum(6.0 * np.abs(w) / b, 1e-6)
        tau_w_field = mu_reg * gamma
        # 每列参考元：|w| 最大且 w>0；非流动元罚为 -1
        w_rank = np.where(w > 0.0, np.abs(w), -1.0)
        ref_row = np.argmax(w_rank, axis=0)            # (nz,)
        has_flow = np.any(w > 0.0, axis=0)            # (nz,)
        col = np.arange(w.shape[1])
        ref_row_safe = np.where(has_flow, ref_row, 0)
        tau_w_ref = tau_w_field[ref_row_safe, col]
        b_ref = b[ref_row_safe, col]
        G = 2.0 * tau_w_ref / np.maximum(b_ref, 1e-12)
        tau_w_extrap = G[None, :] * b / 2.0
        immobile = tau_w_extrap <= f_safety * tau_y
        residual_wall = (cement_ever > 0.0) & (cement_local < c_min_residual)
        wall_new = np.where(immobile | residual_wall, 1.0, 0.0)
        wall_new[:, ~has_flow] = 1.0
        return wall_new.astype(float)
```
> 平行槽流 τw=G·b/2 与 μ 无关，外推不需 μ 差异修正。残泥门保留 `cement_ever>0` 前锋门控（否则水泥未到全域 wall=1 堵死）。

- [ ] **Step 4: __init__ 加 3 参数存 self**（默认值见上）

- [ ] **Step 5: 运行单测通过**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/test_improved_d2dga_annulus.py::TestYieldGateWall -v
```
Expected: 2 PASS。

- [ ] **Step 6: Commit**

```bash
git add cemdisp/models2d/annulus_d2dga.py tests/test_improved_d2dga_annulus.py
git commit -m "feat(annulus): M3 屈服门槛纯逻辑_yield_gate_wall+单测（参考元外推/全冻结回退/残泥兜底）"
```

---

## Task 9: M3 — 接入 run() 泵注分支（停泵分支零改动）

**Files:**
- Modify: `cemdisp/models2d/annulus_d2dga.py`（run() wall 更新 :923-929；解包 :805 取 w/mu/tau_y）
- Test: `tests/test_improved_d2dga_annulus.py`（门开集成 + 停泵不更新）

**Interfaces:**
- Consumes: `_compute_velocity` 在 :805 返回的 `w`(第0)、`mu_reg`(第2)、`tau_y`(第8)；`cement_local`(:927)、`cement_ever`(:928)、`b`（geom["effective_b"]，与 :574-575 一致）。

- [ ] **Step 1: 写失败测试（门关=旧行为；停泵不更新 wall）**

在 `tests/test_improved_d2dga_annulus.py` 追加 `TestYieldGateIntegration`：

```python
class TestYieldGateIntegration:
    def test_gate_off_preserves_cmin_behavior(self):
        # enable_yield_gate=False 下 wall 仍由 c_min 决定（与基线一致）
        ...  # 参照 TestStaticWallLayer.test_wall_consistency_after_run 装配
    def test_gate_on_wall_can_reopen(self):
        # 低排量窄边冻结 -> 高排量后同格 wall 能变 0（可重启）
        ...
    def test_gate_does_not_update_during_shutdown(self):
        # 停泵阶段 wall 保持泵注末值，不因 w_prev 非零误解冻
        ...
```
> 装配参照现有 `TestStaticWallLayer`（:557-667）；两段排量用 mock/inlet_state_provider 返回不同 flow_rate；若太重，至少用 solver 内 `pump_active` 分支单测覆盖停泵不重算。

- [ ] **Step 2: 运行确认失败**

- [ ] **Step 3: 改 wall 更新（:923-929）**

现有：
```python
cement_local = np.clip(lead + tail, 0, 1)
cement_ever = np.maximum(cement_ever, cement_local)
wall = np.where(cement_ever > 0, (cement_local < self.c_min).astype(float), 0.0)
```
改为：
```python
cement_local = np.clip(lead + tail, 0, 1)
cement_ever = np.maximum(cement_ever, cement_local)
if self.enable_yield_gate:
    wall = self._yield_gate_wall(
        w, geom["effective_b"], mu, tau_y, cement_ever, cement_local,
        self.yield_gate_f_safety, self.yield_gate_c_min_residual)
else:
    wall = np.where(cement_ever > 0, (cement_local < self.c_min).astype(float), 0.0)
```
> `w`/`mu`/`tau_y` 是本步泵注分支 :805 解包的当前步值（`mu` 即第2个 mu_reg，按解包别名）。停泵分支 :931-952 **不动**，wall 自然保持上一泵注步值（核查修订3）。

- [ ] **Step 4: 更新 TestStaticWallLayer 断言**

- `test_wall_consistency_after_run`（:571-596）：门关路径保留原断言（wall=1 ⟺ cement<c_min）。
- `test_cmin_0_3_more_wall_than_cmin_0_05`（:624-647）：门关保留；门开下 c_min 仅残泥下限，单调性可能变——新增门开用例断言"τ_y 越大/f_safety 越大 → wall 越多"。
- `test_constructor_has_cmin_default`（:560-564）：追加 `assert s.enable_yield_gate is False`、`assert s.yield_gate_f_safety == 1.15`、`assert s.yield_gate_c_min_residual == 0.01`；保留 `s.c_min == 0.05`（参数保留为门关 fallback）。
- `test_wall_zeros_velocity_in_wall_cells`（:598-622）：保留（wall 机制不变，只改生产方式）。
- `test_wall_zero_before_cement_arrival`（:649-667）：保留；门开回归由 Task 8 单测覆盖。

- [ ] **Step 5: 全量测试零回归 + 门开集成**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/ -q
```
Expected: 全 PASS；门关（默认）逐位复现基线。

- [ ] **Step 6: hu101 三剖面快检（手动）**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python scripts/hu101_standoff_measured_vs_assumed.py
```
验收（设计稿 §6）：between_cent 底部窄边阶段相关启动/静止（快替段尝试重启），非全程冻结；窜槽不拉到 0.992。流变做 ±50% YP 敏感性（无温压数据，仅井口常数）。不达标记录回报而非调常数。

- [ ] **Step 7: Commit**

```bash
git add cemdisp/models2d/annulus_d2dga.py tests/test_improved_d2dga_annulus.py
git commit -m "feat(annulus): M3 屈服门槛接入泵注分支（停泵零改动，enable_yield_gate默认关）"
```

---

## Task 10: M4 — e 护栏（M3 就位后才开 e_clip_max>0.55）

**Files:**
- Modify: `cemdisp/models2d/annulus_d2dga.py`（`__init__` 加参数；`:303` 一行）
- Test: `tests/test_improved_d2dga_annulus.py`（追加）

**Interfaces:**
- Produces: `e_clip_max: float = 0.55`（默认=现行；生产跑道显式设 0.90）

- [ ] **Step 1: 写失败测试**

```python
def test_e_clip_max_default_preserves_055_ceiling():
    s = _make_solver(); assert s.e_clip_max == 0.55
def test_e_clip_max_configurable():
    s = _make_solver(e_clip_max=0.90); assert s.e_clip_max == 0.90
```

- [ ] **Step 2: 运行确认失败**

- [ ] **Step 3: __init__ 加参数；改 :303**

```python
e = np.clip(1.0 - standoff, 0.05, self.e_clip_max)
```
存 self。**不做 b_narrow≥5mm 下限**（设计稿 §5）。

- [ ] **Step 4: 全量测试 + 体积守恒校核**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/ -q
```
手动核 hu101 三剖面（`e_clip_max=0.90`，M3 已开）：between_cent 下部 7055–7868m e 从 0.55→~0.78，体积校正 scale（:330-336）自动重算、half_volume 守恒。**硬约束：本 Task 在 Task 9 之后执行。**

- [ ] **Step 5: Commit**

```bash
git add cemdisp/models2d/annulus_d2dga.py tests/test_improved_d2dga_annulus.py
git commit -m "feat(annulus): M4 e护栏提为e_clip_max参数（默认0.55=现行，生产跑道显式0.90）"
```

---

## Task 11: I3 — eta2 数组透传 + Δρ 局部化

**Files:**
- Modify: `cemdisp/models2d/d2dga_flux.py`（:136 一行 + 2 注解）
- Modify: `cemdisp/models2d/annulus_d2dga.py`（`__init__` 加开关；:879-882 分支；eta2 从 :805 解包取）
- Test: `tests/test_d2dga_flux.py`（加数组 eta2 回归）

**Interfaces:**
- Produces: `enable_local_i3: bool = False`（默认关=基线均值语义）

- [ ] **Step 1: 写失败测试（数组 eta2 不再崩）**

`tests/test_d2dga_flux.py` 追加：
```python
def test_buoyancy_flux_accepts_array_eta2():
    import numpy as np
    from cemdisp.models2d.d2dga_flux import d2dga_buoyancy_flux
    c = np.full(5, 0.5); H = np.full(5, 0.02)
    eta2 = np.linspace(0.05, 0.3, 5)   # 数组，曾触发 max() ValueError
    qp, qx = d2dga_buoyancy_flux(c, 2.0, 100.0, H, eta2, np.ones(5), np.ones(5))
    assert qp.shape == (5,) and np.all(np.isfinite(qp))
```

- [ ] **Step 2: 运行确认失败**（`ValueError: truth value of an array...`）

- [ ] **Step 3: 修 d2dga_flux.py:136**

```python
    coef = (delta_rho * np.asarray(H, dtype=float) ** 3) / (6.0 * np.maximum(eta2, 1.0e-9))
```
把函数签名（:120-128）里 `delta_rho: float`、`eta2: float` 注解改为 `float | np.ndarray`（若文件有 `FloatOrArray` 别名就用）。公式本体 :110-111 不动。

- [ ] **Step 4: __init__ 加 enable_local_i3；改 :879-882 分支**

```python
if self.enable_local_i3:
    eta2_flux = eta2 if np.all(np.isfinite(eta2)) else float(np.mean(eta2))
    delta_rho_flux = (rho - mud_density_gcc) * 1000.0
else:
    eta2_flux = float(np.mean(mu)) if np.all(np.isfinite(mu)) else 0.18
    delta_rho_flux = (rho.mean() - mud_density_gcc) * 1000.0
```
`eta2` 从 :805 解包第 9 个取；`mud_density_gcc` 在 :758 已有；`rho` 第 3 个。把浮力通量调用（:885-888）的实参换成 `eta2_flux`/`delta_rho_flux`。停泵分支不调浮力通量（核实 :931-952），无需改。

- [ ] **Step 5: 全量测试 + R2 消融**

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 /d/apps/Anaconda/python -m pytest tests/ -q
```
门关逐位复现（既有标量 eta2 测试逐位不变）；门开数组测试通过；R2 消融（I3 单独）效率变化不再 +2e-10。

- [ ] **Step 6: Commit**

```bash
git add cemdisp/models2d/d2dga_flux.py cemdisp/models2d/annulus_d2dga.py tests/test_d2dga_flux.py
git commit -m "feat(annulus): I3局部化（max→np.maximum修数组崩溃+eta2场透传+Δρ局部化，enable_local_i3默认关）"
```

---

## Task 12: 附加 B（可选）— hu101 实测 standoff 可选 profile

**Files:**
- Modify: `cemdisp/data/loaders/hu101_loader.py`（:222-234 附近）
- Reference: `scripts/hu101_standoff_measured_vs_assumed.py:48-83`（实测硬编码剖面）

> 独立可选项，不阻塞主线。不改默认（假设值=默认运行），实测剖面用 `dataclasses.replace(base_well, standoff_profile=...)` 作为可选 profile。伴随 ±0.1 区间标注。

- [ ] **Step 1:** 在 hu101_loader.py 暴露可选实测 standoff profile（`dataclasses.replace`，扶正器间 0.22–0.78、扶正器处 0.60–0.88）。
- [ ] **Step 2:** 加 loader 测试（默认仍是 model_assumption；可选 profile 生效且形状匹配）。
- [ ] **Step 3:** 全量测试 + Commit。

---

## Task 13: 全量 8 井重跑 + 论文数字对照表

**Files:**
- Output: `results/<井名>_1D2D耦合模型/`（8 个 runner 输出目录）
- Doc: 提交说明附"受影响论文数字对照表"

- [ ] **Step 1: 重跑节奏**（设计稿 §11）
  1. M0 完成后一轮（零物理，验证即过）
  2. M1 完成后一轮（独立归因）
  3. M3+M2 合并一轮
  4. M4 收尾一轮
  每轮全程 hu101 三剖面快检 + 末轮全量 8 井（`scripts/test_all_wells.py` 列 7 井，ht1_004 单独跑其 runner）。
- [ ] **Step 2: M1 提交必附对照表**：表4（8井效率/窜槽/混浆）、表5（消融）、表6（网格收敛）、67%居中度阈值、呼101"65.04% vs 62.77% +2.27pp"、"20秒"卖点（4处）。同步更新 obsidian 论文框架/撰写指南。
- [ ] **Step 3: 量级锚验收**（设计稿 §1）：M-F19 靶 e=0.3→η_N≥95%、e=0.6→η_N≈30–35%；e=0.6/0.8 体积比 3.25/6.15；hu101 三剖面 η_N 单调、between_cent 不反转、mixing 收敛。
- [ ] **Step 4: 诚实边界**：不得主张三闭包显著提效；不向 66.7% 校准；CFL 消融数据集②未落盘不得引用。流变井下温压修正属数据层后置，结论标注"井口常数流变"。

---

## Self-Review 记录

- **Spec coverage**：M0(Task2-3)、M1(Task4-5)、M3(Task8-9)、M2(Task6-7)、M4(Task10)、I3(Task11)、附加B(Task12)、重跑(Task13) 全覆盖；扩元基础设施 Task1。核查修订 1-11 逐条落入：修1(I3 max)→Task11、修2(3处解包)→Task1、修3(停泵门控)→Task9、修4(参考元)→Task8、修5(无屈服兜底)→Task8、修6(滞后)→Task9 天然、修7(混合n/k)→Task1/6/7、修8(M4排序)→Task10 硬约束、修9(类名/行号)→Task9、修10(runner直接索引)→Global Constraints、修11(Tier0宽except)→已在基线 commit 47f6806 落地。
- **Placeholder scan**：无 TBD/TODO；Task 5/9/12 标"手动/参照现有装配"处给出了具体文件行与判据，非空泛。
- **Type consistency**：12 元组顺序 `..., tau_y(8), eta2(9), n_mix(10), kappa_mix(11)` 在 Task1/7/8/9/11 统一；`_compute_props` 9 元组 `..., n_mix(7), kappa_mix(8)` 在 Task1 统一。开关名 `enable_regime_split`/`enable_yield_gate`/`enable_local_i3`/`e_clip_max`/`dispersion_*` 全文一致。
