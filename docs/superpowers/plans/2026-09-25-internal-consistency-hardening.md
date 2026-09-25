# 内部自洽加固 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 让模型"代码实际行为 = 其声明行为"，并建立可复核的自洽校验台账——不接受任何为让数字好看而做的调整。

**Architecture:** 只做四类改动：① 静默失效的开关改为告警（零数值影响）；② 已存在但未接线的输入常量接进 WellSpec；③ 无参调用导致的硬编码半径改为传参；④ 形参已删的调用点修复 + 自动闸门。随后建立三项自洽校验（守恒 / `η_E` 恒等式 / 前缘解析对照）并产出位移台账。

**Tech Stack:** Python 3.13（conda env `shenjingwangluo`）、numpy 2.3.3、scipy 1.16.2、pytest。

**Spec:** `docs/superpowers/specs/2026-09-25-internal-consistency-hardening-design.md`

## Global Constraints

- 所有 pytest / 脚本运行前必须 `PYTHONIOENCODING=utf-8 PYTHONUTF8=1`（Windows 控制台 GBK）。
- 解释器一律用显式路径 `D:/apps/Anaconda/envs/shenjingwangluo/python.exe`（避免 R-T3-7 陷阱）。
- 代码注释与提交信息用中文；`cemdisp/` 模块名、函数名保持英文。
- **不得改动** `hu101model/`、`hu102model/`（legacy）；**不得改动** 任何有量纲标定钮（`dispersion_alpha`、`yield_gate_f_safety`、`K_AXIAL`、`re_crit`、`μ` clip 上限、`w_prev` 初值）。
- 涉及结果口径的修正（T3 管容）必须**先报后动**，并在位移台账留痕。
- 每个任务独立提交；提交信息前缀 `fix(consistency):` 或 `test(consistency):`。
- 权威结果目录 `results/<井名>_1D2D耦合模型/` 的写入必须走既有防覆写守卫（参见 `tests/contract/test_zhang2022_benchmark_outdir_guard.py`）。

## Review Focus

1. **守卫误报**：白名单写错会让每天都出现的告警变成噪音，最终掩盖真告警——期望"只在真正的空转组合下告警，且一次运行最多一次"。
2. **管容修正的位移方向**：`ht1_001`/`hu2` 管容由 70.88/69.12 升到 94.5/81.13 m³，前缘到鞋时刻**推后**——期望 η 位移被逐井记录而不是被平均掉。
3. **零影响证明**：停泵半径修正**不应该**改变任何生产井结果（该路径在 2D 时间窗内不可达）——期望有一条测试或台账证据证明"零位移"，而不是"没检查"。
4. **恒等式的近似性**：`η_E ≡ 1 − 饥饿份额` 只在 `c̄∈{0,1}` 时严格成立——期望闸门是带阈值的**记录**而非硬断言"相等"。
5. **同名不同义**：`_effective_pipe_radius_m`（管内半径）与 `shoe_lag/shoe_md` 推出的**等效管容半径**是两个量（差 42%）——期望测试只改前者，不动后者。

---

### Task 1: 基线冻结（改动前）

**Files:**
- Create: `docs/superpowers/plans/baseline-2026-09-25.md`
- Create: `scripts/entrypoints/freeze_baseline_20260925.py`
- Test: `tests/contract/test_baseline_freeze_20260925.py`

**Interfaces:**
- Produces: markdown 基线文件（供后续任务比对位移）；`collect_baseline() -> dict` 供测试调用。

- [ ] **Step 1: Write the failing test**

```python
# tests/contract/test_baseline_freeze_20260925.py
import importlib


def test_collect_baseline_shape():
    mod = importlib.import_module("scripts.entrypoints.freeze_baseline_20260925")
    b = mod.collect_baseline()
    assert b["head_commit"] and len(b["head_commit"]) >= 7
    assert b["python"].startswith("3.13")
    assert b["wells"] and all(set(w) >= {"well", "eta_E", "eta_N"} for w in b["wells"])
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 D:/apps/Anaconda/envs/shenjingwangluo/python.exe -m pytest tests/contract/test_baseline_freeze_20260925.py -v`
Expected: FAIL — `ModuleNotFoundError: scripts.entrypoints.freeze_baseline_20260925`

- [ ] **Step 3: Write minimal implementation**

```python
# scripts/entrypoints/freeze_baseline_20260925.py
"""冻结 2026-09-25 基线：HEAD、环境、8 井当前口径数字。只读，不改任何结果目录。"""
from __future__ import annotations

import csv
import platform
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SUMMARY = ROOT / "results" / "源模型口径重跑_2026-09-14" / "汇总.csv"
OUT = ROOT / "docs" / "superpowers" / "plans" / "baseline-2026-09-25.md"


def _git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True,
                          check=True).stdout.strip()


def collect_baseline() -> dict:
    wells = []
    with SUMMARY.open(encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            wells.append({"well": row["well"], "eta_E": float(row["eta_E"]),
                          "eta_N": float(row["eta_N"])})
    return {"head_commit": _git("rev-parse", "HEAD"),
            "dirty": bool(_git("status", "--porcelain")),
            "python": platform.python_version(),
            "numpy": __import__("numpy").__version__,
            "scipy": __import__("scipy").__version__,
            "wells": wells}


def main() -> int:
    b = collect_baseline()
    lines = ["# 基线冻结（2026-09-25，内部自洽加固前）", "",
             f"- HEAD: `{b['head_commit']}`（工作树{'有' if b['dirty'] else '无'}未提交改动）",
             f"- 环境: Python {b['python']} / numpy {b['numpy']} / scipy {b['scipy']}",
             "", "| 井 | η_E | η_N |", "|---|---|---|"]
    lines += [f"| {w['well']} | {w['eta_E']:.6f} | {w['eta_N']:.6f} |" for w in b["wells"]]
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[baseline] {OUT}  ({len(b['wells'])} 井)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run test + generate baseline**

Run: `... -m pytest tests/contract/test_baseline_freeze_20260925.py -v` → PASS
Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 D:/apps/Anaconda/envs/shenjingwangluo/python.exe scripts/entrypoints/freeze_baseline_20260925.py`
Expected: 打印 `[baseline] …baseline-2026-09-25.md (8 井)`；文件内 8 行且与 `results/源模型口径重跑_2026-09-14/汇总.csv` 逐位一致

- [ ] **Step 5: Commit**

```bash
git add docs/superpowers/plans/baseline-2026-09-25.md scripts/entrypoints/freeze_baseline_20260925.py tests/contract/test_baseline_freeze_20260925.py
git commit -m "test(consistency): 冻结 2026-09-25 基线（HEAD/环境/8 井数字）"
```

---

### Task 2: 死开关守卫完备化（零数值影响）

**Files:**
- Modify: `cemdisp/models2d/annulus_d2dga.py`（替换 `:561-576` 的"静默无效开关"守卫块）
- Modify: `cemdisp/models2d/two_layer.py:172` 附近（`isotropic_flux_q0` docstring，D2）
- Modify: `docs/源模型口径与适用域声明.md` §1 声明 5（D2）
- Test: `tests/contract/test_dead_switch_guard.py`（新建）

**Interfaces:**
- Produces: 模块级纯函数 `_dead_switches(**switches) -> list[str]`；构造期告警文本含每个空转开关名

- [ ] **Step 1: Write the failing test**

```python
# tests/contract/test_dead_switch_guard.py
"""A3 惯例：置真但当前路径无消费者的开关必须一次性告警；有消费者时不得告警。"""
import warnings

import pytest

from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver, _dead_switches

OLD_PATH_ONLY = ["enable_regime_split", "enable_true_buoyancy", "enable_power_law_gap_law"]


@pytest.mark.parametrize("name", OLD_PATH_ONLY)
def test_old_path_switch_is_dead_on_new_path(name):
    dead = _dead_switches(enable_stream_function=True, **{name: True})
    assert any(name in d for d in dead)


@pytest.mark.parametrize("name", OLD_PATH_ONLY)
def test_old_path_switch_is_alive_on_old_path(name):
    dead = _dead_switches(enable_stream_function=False, **{name: True})
    assert not any(name in d for d in dead)


def test_yield_gate_dead_until_stream_yield_gate_on():
    dead = _dead_switches(enable_stream_function=True, enable_yield_gate=True,
                          enable_stream_yield_gate=False)
    assert any("enable_yield_gate" in d for d in dead)
    alive = _dead_switches(enable_stream_function=True, enable_yield_gate=True,
                           enable_stream_yield_gate=True)
    assert not any("enable_yield_gate" in d for d in alive)


def test_constructor_warns_once_and_names_switch():
    with pytest.warns(UserWarning, match="enable_regime_split"):
        AnnulusD2DGASolver(total_t=1.0, nz=4, ny=9, enable_regime_split=True)


def test_constructor_quiet_when_all_alive():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        AnnulusD2DGASolver(total_t=1.0, nz=4, ny=9, enable_stream_function=False,
                           enable_regime_split=True)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `... -m pytest tests/contract/test_dead_switch_guard.py -v`
Expected: FAIL — `ImportError: cannot import name '_dead_switches'`

- [ ] **Step 3: Write minimal implementation**

在 `cemdisp/models2d/annulus_d2dga.py` 模块层（`_limit_phase_volume` 之前）新增纯函数：

```python
_OLD_PATH_ONLY_SWITCHES = (
    ("enable_regime_split", "M2 局部流态修正"),
    ("enable_true_buoyancy", "真浮力体力"),
    ("enable_power_law_gap_law", "幂律缝隙律"),
)


def _dead_switches(**switches) -> list[str]:
    """返回"置真但在当前路径上无消费者"的开关名（A3 惯例，纯函数便于测试）。

    规则（与代码实际消费点一一对应；改消费点时必须同步改这里）：
      * 旧代数路径专属开关：仅在 ``enable_stream_function=False`` 时被消费
        （``_compute_velocity`` 在新路径早退，见 :1826-1831）；
      * ``enable_yield_gate``：新路径下 ``wall`` 只写诊断量，需
        ``enable_stream_yield_gate=True`` 才进流函数算子（:1829）；
      * ``K_AXIAL`` 为模块常量、非构造形参，仅旧路径消费，附在告警文本里说明。
    """
    dead: list[str] = []
    new_path = bool(switches.get("enable_stream_function", True))
    if new_path:
        for name, desc in _OLD_PATH_ONLY_SWITCHES:
            if switches.get(name, False):
                dead.append(f"{name}（{desc}：仅旧代数路径 enable_stream_function=False 消费）")
        if switches.get("enable_yield_gate", True) and not switches.get(
                "enable_stream_yield_gate", False):
            dead.append("enable_yield_gate（wall 已算出但不进算子：需 "
                        "enable_stream_yield_gate=True 才有动力学作用，当前仅写诊断量）")
    return dead
```

把 `:561-576` 的原守卫块整体替换为：

```python
        _dead = _dead_switches(
            enable_stream_function=enable_stream_function,
            enable_yield_gate=enable_yield_gate,
            enable_stream_yield_gate=enable_stream_yield_gate,
            enable_regime_split=enable_regime_split,
            enable_true_buoyancy=enable_true_buoyancy,
            enable_power_law_gap_law=enable_power_law_gap_law,
        )
        if _dead:
            warnings.warn(
                "AnnulusD2DGASolver 的开关 " + "；".join(_dead)
                + " 在当前配置下无效（新路径不消费；K_AXIAL=1/3 同属旧路径专属）。"
                "请修正配置或移除这些开关。", UserWarning, stacklevel=2)
```

保留 `:534-556` 的 `_hb_dead` 块与 HB×B-3 互斥告警**不动**。

同一步内完成 D2 的文字一致化：在 `two_layer.isotropic_flux_q0` docstring 与 `docs/源模型口径与适用域声明.md` §1 声明 5 各加一句："默认路径（`enable_stream_function=True`）下 q₀ 无消费方，仅作 `gap_solver` 回归锚；弥散由 I₃ 单独承载。"

- [ ] **Step 4: Run test + 全量契约 + 不变性证明**

Run: `... -m pytest tests/contract/test_dead_switch_guard.py tests/contract/test_stream_yield_gate.py tests/contract/test_yield_deadzone.py -v` → PASS
Run: `... -m pytest tests/contract/ -q` → 除既存 1 红（`test_hb_newtonian_limit.py::test_frozen_snapshot_stream_function_defaults`）外全绿
Run（逐位不变性）：`... -m pytest tests/contract/test_six_well_integration.py -v` → PASS（该测试逐位锚定默认路径 ⇒ 证明本次改动零数值影响）

- [ ] **Step 5: Commit**

```bash
git add cemdisp/models2d/annulus_d2dga.py cemdisp/models2d/two_layer.py docs/源模型口径与适用域声明.md tests/contract/test_dead_switch_guard.py
git commit -m "fix(consistency): 死开关守卫完备化 + q₀ 声明一致化（零数值影响）"
```

---

### Task 3: 管容锚点接入（ht1_001 / hu2）★需先报后动

**Files:**
- Modify: `cemdisp/data/loaders/ht1_001_loader.py:383-400`（WellSpec 构造）
- Modify: `cemdisp/data/loaders/hu2_loader.py:292-318`（WellSpec 构造）
- Test: `tests/contract/test_shoe_lag_wiring.py`（新建）

**Interfaces:**
- Consumes: `HT1_001_SHOE_LAG_VOLUME_M3`（`:137`）、`HU2_SHOE_LAG_VOLUME_M3`（`:101`）
- Produces: `load_ht1_001_tailpipe()[0].shoe_lag_volume_m3 == 94.5`；`load_hu2_tailpipe()[0].shoe_lag_volume_m3 ≈ 81.13`

- [ ] **Step 1: Write the failing test**

```python
# tests/contract/test_shoe_lag_wiring.py
"""T1（管容死常量）：常量必须真正进入 WellSpec，否则 1D 管容回退到单一内径。"""
import math

import pytest

from cemdisp.data.loaders import ht1_001_loader, hu2_loader


def test_ht1_001_shoe_lag_volume_wired():
    well, _, _, _ = ht1_001_loader.load_ht1_001_tailpipe()
    assert well.shoe_lag_volume_m3 is not None, "常量未接入 WellSpec（T1 缺陷）"
    assert well.shoe_lag_volume_m3 == pytest.approx(
        ht1_001_loader.HT1_001_SHOE_LAG_VOLUME_M3, rel=1e-9)


def test_hu2_shoe_lag_volume_wired():
    well, _, _, _ = hu2_loader.load_hu2_tailpipe()
    assert well.shoe_lag_volume_m3 is not None, "常量未接入 WellSpec（T1 缺陷）"
    assert well.shoe_lag_volume_m3 == pytest.approx(
        hu2_loader.HU2_SHOE_LAG_VOLUME_M3, rel=1e-9)


def test_tube_volume_larger_than_single_id_fallback():
    """修正后管容必须显著大于旧的单一内径口径（旧值 ht1_001≈70.88 / hu2≈69.12 m³）。"""
    well, _, _, _ = ht1_001_loader.load_ht1_001_tailpipe()
    old = math.pi * (well.liner_id_mm / 2000.0) ** 2 * well.shoe_md_m
    assert well.shoe_lag_volume_m3 > old * 1.20      # 旧口径低约 25%
```

- [ ] **Step 2: Run test to verify it fails**

Run: `... -m pytest tests/contract/test_shoe_lag_wiring.py -v`
Expected: FAIL — `assert None is not None` 或 `AttributeError: 'WellSpec' object has no attribute 'shoe_lag_volume_m3'`

- [ ] **Step 3: Write minimal implementation**

`ht1_001_loader.py` 的 `WellSpec(...)`（`:383`）中，在 `liner_id_mm=HT1_001_LINER_ID_MM,` 之后加一行（与 `ht1_003_loader.py:395` 同构）：

```python
        shoe_lag_volume_m3=HT1_001_SHOE_LAG_VOLUME_M3,
```

`hu2_loader.py` 的 `WellSpec(...)`（`:292`）同位置加：

```python
        shoe_lag_volume_m3=HU2_SHOE_LAG_VOLUME_M3,
```

字段名以 `ht1_003_loader.py:395` 为准（该处已在用，字段必然存在）。

- [ ] **Step 4: Run test + 对照旧口径**

Run: `... -m pytest tests/contract/test_shoe_lag_wiring.py tests/contract/test_casing_flow.py -v` → PASS
Run（打印新旧管容对照）：
`PYTHONIOENCODING=utf-8 PYTHONUTF8=1 D:/apps/Anaconda/envs/shenjingwangluo/python.exe -c "import math;from cemdisp.data.loaders import ht1_001_loader as a, hu2_loader as b;[print(m.__name__, w.shoe_lag_volume_m3, math.pi*(w.liner_id_mm/2000)**2*w.shoe_md_m) for m,f in ((a,a.load_ht1_001_tailpipe),(b,b.load_hu2_tailpipe)) for w,_,_,_ in [f()]]"`
Expected: 94.5 vs 70.88（−25%）、81.13 vs 69.12（−15%）

- [ ] **Step 5: 先报后动 → 重跑两口井 → 记录位移**

**此步开工前必须向用户报告并取得同意**（用户全局规则：涉及数据口径的决定先报后动）。
Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 D:/apps/Anaconda/envs/shenjingwangluo/python.exe scripts/entrypoints/rerun_all_wells_corrected.py`（若该脚本不支持子集，复制其 `run_one` 到临时探针脚本，仅跑 ht1_001/hu2；**不得**覆盖 `results/<井名>_1D2D耦合模型/`）
把两口井修正前后 η_E/η_N 写入 `results/内部自洽加固_2026-09-25/位移台账.csv`。

- [ ] **Step 6: Commit**

```bash
git add cemdisp/data/loaders/ht1_001_loader.py cemdisp/data/loaders/hu2_loader.py tests/contract/test_shoe_lag_wiring.py results/内部自洽加固_2026-09-25/位移台账.csv
git commit -m "fix(consistency): ht1_001/hu2 管容锚点接入 WellSpec（T1 死常量）"
```

---

### Task 4: 停泵沉降半径传参（并证明生产零影响）

**Files:**
- Modify: `cemdisp/transport1d/casing_flow.py:1371`（`_settled_exit_fluid_name_enhanced` 内的调用）
- Test: `tests/contract/test_settled_pipe_radius.py`（新建）

**Interfaces:**
- Consumes: `CasingFlowSolver._effective_pipe_radius_m(well_spec=None) -> float`（`:1303-1316`）
- Produces: 无新接口；仅调用点补传 `well_spec`

- [ ] **Step 1: Write the failing test**

```python
# tests/contract/test_settled_pipe_radius.py
"""T5：停泵沉降的 τ_c 必须用本井内径（R = liner_id/2），不得用硬编码 0.05 m。"""
import inspect

from cemdisp.transport1d.casing_flow import CasingFlowSolver


class _WellStub:
    liner_id_mm = 111.16


def test_effective_radius_uses_well_spec():
    s = CasingFlowSolver(enable_gravity=True)
    assert s._effective_pipe_radius_m(_WellStub()) == 111.16 / 2000.0
    assert s._effective_pipe_radius_m(None) == 0.05


def test_settled_path_passes_well_spec():
    src = inspect.getsource(CasingFlowSolver._settled_exit_fluid_name_enhanced)
    assert "_effective_pipe_radius_m()" not in src, (
        "停泵沉降仍以无参调用取半径（恒 0.05 m，T5 缺陷）")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `... -m pytest tests/contract/test_settled_pipe_radius.py -v`
Expected: `test_effective_radius_uses_well_spec` PASS；`test_settled_path_passes_well_spec` FAIL

- [ ] **Step 3: Write minimal implementation**

在 `_settled_exit_fluid_name_enhanced` 内把 `:1371` 的

```python
            pipe_radius_m = self._effective_pipe_radius_m()
```

改为

```python
            pipe_radius_m = self._effective_pipe_radius_m(well_spec)
```

若该方法签名中无 `well_spec`，加形参 `well_spec: WellSpec | None = None`，并同步其唯一调用方 `pipe_exit_state_at`（保持向后兼容）。

- [ ] **Step 4: Run test + 零影响证明**

Run: `... -m pytest tests/contract/test_settled_pipe_radius.py tests/contract/test_casing_flow.py tests/contract/test_shoe_timeline.py -v` → PASS
Run（零影响）：`... -m pytest tests/contract/test_six_well_integration.py -v` → PASS
若该测试失败 ⇒ 说明该路径其实可达，**停下并上报**（不得调阈值掩盖）。

- [ ] **Step 5: Commit**

```bash
git add cemdisp/transport1d/casing_flow.py tests/contract/test_settled_pipe_radius.py
git commit -m "fix(consistency): 停泵沉降半径改传 well_spec（T5 硬编码 0.05 m）"
```

---

### Task 5: 调用签名闸门 + 修复 4 处坏调用点

**Files:**
- Create: `scripts/entrypoints/check_call_signatures.py`
- Modify: `cemdisp/runners/ht1_004_ablation.py:175-181`
- Modify: `scripts/entrypoints/closure_contribution_scan.py:73`
- Modify: `scripts/probes/density_contrast_sensitivity_scan.py`
- Modify: `scripts/probes/dispersion_scale_sensitivity_scan.py`
- Test: `tests/contract/test_entrypoint_signatures.py`（新建）

**Interfaces:**
- Produces: `find_bad_kwargs(root: Path) -> list[tuple[str, int, str, str]]`（文件、行号、类名、未知形参名）

- [ ] **Step 1: Write the failing test**

```python
# tests/contract/test_entrypoint_signatures.py
"""形参被删却仍被传 ⇒ 运行即 TypeError。本闸门冻结全部生产调用点。"""
from pathlib import Path

from scripts.entrypoints.check_call_signatures import find_bad_kwargs

ROOT = Path(__file__).resolve().parents[2]


def test_no_bad_keyword_arguments_anywhere():
    bad = find_bad_kwargs(ROOT)
    assert bad == [], "存在未知形参调用点：" + "; ".join(
        f"{f}:{n}:{c}.{k}" for f, n, c, k in bad)


def test_scanner_detects_a_known_bad_call(tmp_path):
    p = tmp_path / "sample.py"
    p.write_text("from cemdisp.models2d import AnnulusD2DGASolver\n"
                 "s = AnnulusD2DGASolver(enable_d2dga_auto_m=True)\n", encoding="utf-8")
    assert "enable_d2dga_auto_m" in {k for *_, k in find_bad_kwargs(tmp_path)}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `... -m pytest tests/contract/test_entrypoint_signatures.py -v`
Expected: FAIL — `ModuleNotFoundError: scripts.entrypoints.check_call_signatures`

- [ ] **Step 3: Write minimal implementation**

```python
# scripts/entrypoints/check_call_signatures.py
"""扫描全仓对 AnnulusD2DGASolver / CasingFlowSolver 的调用关键字，比对真实签名。

用途：形参删除后（如 enable_d2dga_auto_m、dispersion_dt_scale）静默腐烂的调用点
只在运行期抛 TypeError；本扫描器在测试期就抓住它。
"""
from __future__ import annotations

import ast
import inspect
import sys
from pathlib import Path

TARGETS = {"AnnulusD2DGASolver", "CasingFlowSolver"}
SKIP_DIRS = {".git", "__pycache__", "archive", "hu101model", "hu102model", ".tmp_research"}


def _signatures() -> dict[str, set[str]]:
    from cemdisp.models2d import AnnulusD2DGASolver
    from cemdisp.transport1d import CasingFlowSolver

    out: dict[str, set[str]] = {}
    for cls in (AnnulusD2DGASolver, CasingFlowSolver):
        params = inspect.signature(cls.__init__).parameters
        out[cls.__name__] = {n for n in params if n != "self"}
    return out


def find_bad_kwargs(root: Path) -> list[tuple[str, int, str, str]]:
    sigs = _signatures()
    bad: list[tuple[str, int, str, str]] = []
    for path in Path(root).rglob("*.py"):
        if any(part in SKIP_DIRS for part in path.parts):
            continue
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"))
        except (SyntaxError, UnicodeDecodeError):
            continue
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            fn = node.func
            name = fn.id if isinstance(fn, ast.Name) else getattr(fn, "attr", None)
            if name not in TARGETS:
                continue
            for kw in node.keywords:
                if kw.arg and kw.arg not in sigs[name]:
                    bad.append((str(path), node.lineno, name, kw.arg))
    return bad


def main() -> int:
    bad = find_bad_kwargs(Path(__file__).resolve().parents[2])
    for f, n, c, k in bad:
        print(f"{f}:{n}: {c}({k}=...) 形参不存在")
    print(f"[signatures] {len(bad)} 处问题")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: 按扫描结果修 4 处**

删除各调用点已被移除的形参（`enable_d2dga_auto_m` / `dispersion_dt_scale`）；若某构造参数确实仍需要，改从现有 `level` 对象的既有字段读取，**不得**重新引入已删形参。

Run: `... -m pytest tests/contract/test_entrypoint_signatures.py -v` → PASS
Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 D:/apps/Anaconda/envs/shenjingwangluo/python.exe scripts/entrypoints/check_call_signatures.py` → `[signatures] 0 处问题`

- [ ] **Step 5: Commit**

```bash
git add scripts/entrypoints/check_call_signatures.py tests/contract/test_entrypoint_signatures.py cemdisp/runners/ht1_004_ablation.py scripts/entrypoints/closure_contribution_scan.py scripts/probes/density_contrast_sensitivity_scan.py scripts/probes/dispersion_scale_sensitivity_scan.py
git commit -m "fix(consistency): 修复 4 处已删形参调用点 + 新增调用签名闸门"
```

---

### Task 6: 自洽校验闸门（恒等式 / 前缘对照 / 守恒）

**Files:**
- Create: `cemdisp/diagnostics/internal_consistency.py`
- Create: `scripts/entrypoints/verify_internal_consistency.py`
- Test: `tests/contract/test_internal_consistency.py`（新建）

**Interfaces:**
- Consumes: `cemdisp.models2d.annulus_d2dga._trapez2d`；`cemdisp.runners.zhang2022_benchmark.mass_conservation_error`（`:430`）
- Produces:
  - `domain_eta_e(cement, geom) -> float` = `∬b·c / ∬b`
  - `starved_volume_fraction(cement, geom, threshold=0.5) -> float` = `∬b·1[c<thr] / ∬b`
  - `front_position_m(cement, geom, level=0.5) -> float`

- [ ] **Step 1: Write the failing test**

```python
# tests/contract/test_internal_consistency.py
import numpy as np

from cemdisp.diagnostics.internal_consistency import (
    domain_eta_e, front_position_m, starved_volume_fraction)


def _geom(ny=9, nz=5):
    return {"phi": np.linspace(0, 1, ny), "s": np.linspace(0, 40, nz),
            "md": np.linspace(7000, 6960, nz), "H": np.full((ny, nz), 0.01),
            "b": np.full((ny, nz), 0.02), "y": np.linspace(0, 0.3, ny)}


def test_identity_exact_for_binary_field():
    g = _geom()
    c = np.zeros((9, 5))
    c[:5] = 1.0                       # 二值场 ⇒ 恒等式严格成立
    assert abs((1.0 - domain_eta_e(c, g)) - starved_volume_fraction(c, g)) < 1e-12


def test_identity_is_approximate_for_mixed_field():
    g = _geom()
    c = np.full((9, 5), 0.4)          # 非二值 ⇒ 有"空隙"，须以阈值记录而非硬断言
    assert abs((1.0 - domain_eta_e(c, g)) - starved_volume_fraction(c, g)) > 0.1


def test_front_position_finite_and_monotone():
    g = _geom()
    c = np.zeros((9, 5))
    c[:, :3] = 1.0                    # s 小的一侧（鞋口侧）已到水泥
    assert front_position_m(c, g) > 0.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `... -m pytest tests/contract/test_internal_consistency.py -v`
Expected: FAIL — `ModuleNotFoundError: cemdisp.diagnostics.internal_consistency`

- [ ] **Step 3: Write minimal implementation**

```python
# cemdisp/diagnostics/internal_consistency.py
"""内部自洽校验量：域内效率、饥饿体积份额、前缘位置。

⚠️ 恒等式 ``η_E = 1 − 饥饿份额`` **只在 c̄∈{0,1} 时严格成立**；现场浓度场连续，
故生产口径下该式是**带阈值的记录**，不是硬断言（历史实测偏差 ≤0.007）。
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from cemdisp.models2d.annulus_d2dga import _trapez2d


def domain_eta_e(cement: NDArray, geom: dict) -> float:
    """域内体积加权效率 η_E = ∬b·c / ∬b（与 `_evaluation_window_efficiencies` 同口径）。"""
    b = geom["b"]
    return float(_trapez2d(b * np.asarray(cement, float), geom)
                 / max(_trapez2d(b, geom), 1e-12))


def starved_volume_fraction(cement: NDArray, geom: dict, threshold: float = 0.5) -> float:
    """饥饿体积份额 ∬b·1[c<threshold] / ∬b（口径同 2026-09-11 探针脚本）。"""
    b = geom["b"]
    ind = (np.asarray(cement, float) < threshold).astype(float)
    return float(_trapez2d(b * ind, geom) / max(_trapez2d(b, geom), 1e-12))


def front_position_m(cement: NDArray, geom: dict, level: float = 0.5) -> float:
    """前缘位置（m）：b 加权列均值首达 ``level`` 的最深 md；未达则返回域底 md。"""
    b = np.asarray(geom["b"], float)
    col = (b * np.asarray(cement, float)).sum(axis=0) / np.maximum(b.sum(axis=0), 1e-30)
    reached = np.where(col >= level)[0]
    md = np.asarray(geom["md"], float)
    return float(md[reached.min()]) if reached.size else float(md.max())
```

`scripts/entrypoints/verify_internal_consistency.py`：对每口井跑一次**短窗**（`nz=60`，仅作自洽诊断，**产物不得当论文数字**），输出 `results/内部自洽加固_2026-09-25/一致性台账.csv`，列：
`井名, 域内eta_E, 饥饿份额, 恒等式偏差, 1D尾浆到鞋时刻_s, 2D前缘位置_m, 质量守恒误差, 说明`

- `1D尾浆到鞋时刻_s` 取 `CasingFlowResult.cement_end_time_s`；
- `质量守恒误差` 调 `zhang2022_benchmark.mass_conservation_error`——**先读其 docstring 确认参数语义**，按其文档口径调用（参照 `zhang2022_benchmark.py:488`）；若其假设对生产井不成立，该列写 `未测` 并在 `说明` 列写原因（**不得**自造公式）。

- [ ] **Step 4: Run test + 生成台账**

Run: `... -m pytest tests/contract/test_internal_consistency.py -v` → PASS
Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 D:/apps/Anaconda/envs/shenjingwangluo/python.exe scripts/entrypoints/verify_internal_consistency.py`
Expected: 生成 CSV；逐井检查 `恒等式偏差` 是否与历史量级（≤0.01）相符。**若某井偏差 >0.05，如实记录并在 `说明` 列标注**（诊断输出，不是验收门）。

- [ ] **Step 5: Commit**

```bash
git add cemdisp/diagnostics/internal_consistency.py scripts/entrypoints/verify_internal_consistency.py tests/contract/test_internal_consistency.py "results/内部自洽加固_2026-09-25/一致性台账.csv"
git commit -m "feat(consistency): 自洽校验闸门（恒等式/前缘/守恒）+ 井台账"
```

---

### Task 7: 位移台账汇总 + 官方口径对照报告（不改口径）

**Files:**
- Create: `scripts/entrypoints/collect_displacement_ledger.py`
- Create: `results/内部自洽加固_2026-09-25/口径对照_官方脚本vs生产runner.md`
- Test: `tests/contract/test_displacement_ledger.py`

**Interfaces:**
- Consumes: `docs/superpowers/plans/baseline-2026-09-25.md`（Task 1）、`results/内部自洽加固_2026-09-25/位移台账.csv`（Task 3）
- Produces: 汇总后的 `位移台账.csv`（含基线行 + 修正行）

- [ ] **Step 1: Write the failing test**

```python
# tests/contract/test_displacement_ledger.py
import importlib


def test_ledger_columns():
    mod = importlib.import_module("scripts.entrypoints.collect_displacement_ledger")
    assert mod.COLUMNS == ["井名", "修正项", "修正前eta_E", "修正后eta_E", "Δeta_E_pp",
                           "修正前eta_N", "修正后eta_N", "Δeta_N_pp", "说明"]
```

- [ ] **Step 2: Run test to verify it fails** → FAIL（模块不存在）

- [ ] **Step 3: Write minimal implementation**

读基线 md 表格 + 现有位移 csv，按 `COLUMNS` 归并写出；`Δ*_pp = (后 − 前) × 100`，保留 2 位小数。

- [ ] **Step 4: 写官方口径对照报告**（只写文档，不改脚本）

固定三部分：① 4 处口径差逐条列出（`rerun_all_wells_corrected.py:73` 不传 1D 三开关、`:79` 不传 `schedule=`、`total_t` 公式不同、额外 `CORRECTED_KW`）；② 已知量级（1D 三开关 ≤0.07 pp、hu2 −0.66 pp；日程与 `total_t` 差**未量化**）；③ 结论：**G3 八井数字非生产 runner 口径，两者不可横比**，是否收口待裁定（Spec §4 R-A）。

- [ ] **Step 5: Run + Commit**

```bash
git add scripts/entrypoints/collect_displacement_ledger.py tests/contract/test_displacement_ledger.py "results/内部自洽加固_2026-09-25/"
git commit -m "test(consistency): 位移台账汇总 + 官方口径对照报告（口径不改，待裁定）"
```

---

## Self-Review

**Spec coverage**：D1→Task 2；D2→Task 2 Step 3 末段；D3→Task 3；D4→Task 4；D5→Task 5；D6→Task 7 Step 4；验收判据 §3.1/3.2→Task 2 Step 4；§3.3→Task 6；§3.4→Task 3 Step 5 + Task 7；§3.5 覆盖全任务。

**Placeholder scan**：无 TBD/TODO。Task 6 Step 3 关于守恒列要求"先读 docstring 再按其口径调用，否则写未测"是**明确的失败处理路径**，不是占位。

**Type consistency**：`_dead_switches`（Task 2）；`find_bad_kwargs` 返回四元组（Task 5，测试解包一致）；`domain_eta_e` / `starved_volume_fraction` / `front_position_m`（Task 6，测试与实现命名一致）。

**Review Focus 覆盖**：①→Task 2 `test_constructor_quiet_when_all_alive`；②→Task 3 Step 5 台账；③→Task 4 Step 4；④→Task 6 `test_identity_is_approximate_for_mixed_field`；⑤→Task 4 的 stub 只测 `_effective_pipe_radius_m`，不触碰管容半径。

---

## 追加任务（用户 2026-09-25 追加要求："时间太长是不行的" + "模型要符合物理实际，不只是代码跑通"）

### Task 8: 线性求解换带状 Cholesky（**纯数值等价**，~2×）

**背景**：默认路径单步 30.4 ms，其中 `scipy.sparse.linalg.spsolve` 占 **64%**（19.45 ms，`stream_function.py:395`）。该算子为**均匀网格 5 点变系数 Poisson**、内部子块精确对称正定（实测 κ=3.16e4），各向异性 10⁶–10⁷ ⇒ Krylov/AMG/ILU 均实测负收益，唯一正收益路线是**节点重排 + 带状 Cholesky**（实测 17.46 → 4.38 ms = 4.0×，与 spsolve 解相对差 5e-14）。**离散系统与物理口径完全不变**。

**Files:**
- Modify: `cemdisp/models2d/stream_function.py:341-398`（矩阵装配与求解段）
- Modify: `cemdisp/models2d/annulus_d2dga.py`（新增构造形参 `enable_banded_solve: bool = True` 并透传）
- Test: `tests/contract/test_banded_solve_equivalence.py`（新建）

**Interfaces:**
- Produces: `solve_stream_function(..., banded: bool = True)`；同函数在 `banded=False` 时**逐位等于改动前**（供历史锚测试使用）
- Consumes: `a_face`/`c_face`（`:341-343`）、Dirichlet 行（`i=0` 为 0、`i=ny-1` 为 1）

- [ ] **Step 1: Write the failing test**

```python
# tests/contract/test_banded_solve_equivalence.py
"""带状 Cholesky 与 spsolve 解同一离散系统 ⇒ Ψ 相对差应为舍入级。"""
import numpy as np
import pytest

from cemdisp.models2d.stream_function import solve_stream_function


def _geom(ny=40, nz=250):
    rng = np.random.default_rng(0)
    phi = np.linspace(0.0, 1.0, ny)
    H = 0.01 * (1.0 + 0.3 * np.cos(np.pi * phi))[:, None] * np.ones((1, nz))
    H = 0.01 + 0.002 * rng.random((ny, nz)) + H * 0.0
    return {"phi": phi, "s": np.linspace(0.0, 2468.0, nz), "H": H,
            "hole_mm": np.full(nz, 215.9), "od_mm": np.full(nz, 168.3),
            "y": phi * np.pi * 0.1, "b": 2.0 * H}


@pytest.mark.parametrize("ny,nz", [(9, 8), (40, 60)])
def test_banded_matches_sparse(ny, nz):
    g = _geom(ny, nz)
    c = np.clip(0.5 + 0.4 * np.sin(np.linspace(0, 3, nz))[None, :] * np.ones((ny, 1)), 0, 1)
    b = np.zeros((2, ny, nz))
    b[0] = 0.1 * np.cos(np.pi * g["phi"])[:, None]
    psi_sp = solve_stream_function(g, c, 0.058, 0.171, 0.34, b, banded=False)
    psi_bd = solve_stream_function(g, c, 0.058, 0.171, 0.34, b, banded=True)
    assert np.all(np.isfinite(psi_bd))
    rel = np.max(np.abs(psi_bd - psi_sp)) / max(np.max(np.abs(psi_sp)), 1e-30)
    assert rel < 1e-12
    # 端点 Dirichlet 口径不变
    assert psi_bd[0, :] == pytest.approx(0.0, abs=1e-14)
    assert psi_bd[-1, :] == pytest.approx(1.0, abs=1e-14)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 D:/apps/Anaconda/envs/shenjingwangluo/python.exe -m pytest tests/contract/test_banded_solve_equivalence.py -v`
Expected: FAIL — `TypeError: solve_stream_function() got an unexpected keyword argument 'banded'`

- [ ] **Step 3: Write minimal implementation**

在 `solve_stream_function` 内，把现有"5 点矩阵装配 + `spsolve`"保留为 `banded=False` 分支；新增 `banded=True`（默认）分支：

```python
def _solve_banded_interior(a_face, c_face, src, *, ny: int, nz: int,
                           dphi: float, dxi: float, r_a: float):
    """消去 φ 两端 Dirichlet 后，用带状 Cholesky 解内部 SPD 系统。

    节点重排 m = j·(ny-2) + i（ξ 外层、φ 内层）⇒ 带宽 = ny-2（生产 ny=40 ⇒ 38）。
    列 j 的 φ 两端 Dirichlet 值：i=0 → 0、i=ny-1 → 1（单位通量 BC）。
    """
    from scipy.linalg import cho_solve_banded, cholesky_banded

    n_int = ny - 2
    band = n_int            # 上/下带宽
    ab = np.zeros((band + 1, n_int * nz))   # LAPACK 带状下三角存法
    rhs = np.zeros(n_int * nz)
    for j in range(nz):
        for i in range(1, ny - 1):
            m = j * n_int + (i - 1)
            a_w = a_face[i - 1, j] / dphi ** 2          # φ⁻ 面
            a_e = a_face[i, j] / dphi ** 2              # φ⁺ 面
            c_s = (c_face[i, j - 1] / dxi ** 2) if j > 0 else 0.0
            c_n = (c_face[i, j] / dxi ** 2) if j < nz - 1 else 0.0
            ab[0, m] = a_w + a_e + c_s + c_n            # 对角
            if i - 1 > 0:
                ab[1, m - 1] = -a_w                     # 下一条（φ⁻ 邻格）
            if j > 0:
                ab[band, m - band] = -c_s               # 最远一条（ξ⁻ 邻格）
            rhs[m] = src[i, j]
        # 边界贡献已含在 Dirichlet 反代：i=1 处无 φ⁻ 项（i=0 已知 0），
        # i=ny-2 处 φ⁺ 项作用于 i=ny-1（值 1）⇒ rhs 加 a_e·1
        rhs[j * n_int + (n_int - 1)] += a_face[ny - 2, j] / dphi ** 2 * 1.0
    factor = cholesky_banded(ab, lower=True)
    sol = cho_solve_banded((factor, True), rhs)
    psi = np.zeros((ny, nz))
    for j in range(nz):
        psi[1:-1, j] = sol[j * n_int:j * n_int + n_int]
        psi[-1, j] = 1.0
    return psi
```

实现要点：① `ab` 的带状索引严格按 LAPACK `?pbtrf`（`ab[0]` 对角、`ab[k]` 第 k 条下对角线）；② 每步只需重建 `ab` 与 `rhs`（系数随 I₁ 变），不做符号分解复用；③ `banded=False` 分支保持原 `spsolve` 代码**逐字不动**以保住历史锚。在 `annulus_d2dga.py` 加形参 `enable_banded_solve: bool = True` 并在 `_velocity_stream_function` / `_solve_stream_function_hb` 两处透传。

- [ ] **Step 4: Run tests + 计时**

Run: `... -m pytest tests/contract/test_banded_solve_equivalence.py tests/contract/test_stream_function.py tests/contract/test_stream_function_solver_integration.py -v` → PASS
Run（计时对照）：`... -m pytest tests/contract/test_banded_solve_equivalence.py -v -s` 中打印两个 ny/nz 组合的单次 `spsolve` vs 带状耗时，人工确认 ≥3×。
Run: `... -m pytest tests/contract/ -q` → 除既存 1 红全绿；**`test_six_well_integration.py` 若失败属预期**（它逐位锚定 spsolve 口径）⇒ 按 Step 5 重锚。

- [ ] **Step 5: 重锚历史锚测试并留位移**

把 `test_six_well_integration.py` 的逐位锚改为 `enable_banded_solve=False`（保留历史可比性），另**新增**一条 `enable_banded_solve=True` 的容差锚（η_E/η_N 相对差 ≤1e-3 pp）；并在 `results/内部自洽加固_2026-09-25/位移台账.csv` 追加"线性求解替换"行（呼1-003 或任一口可跑井）。

- [ ] **Step 6: Commit**

```bash
git add cemdisp/models2d/stream_function.py cemdisp/models2d/annulus_d2dga.py tests/contract/test_banded_solve_equivalence.py tests/contract/test_six_well_integration.py "results/内部自洽加固_2026-09-25/位移台账.csv"
git commit -m "perf(consistency): 流函数线性求解换带状 Cholesky（纯数值等价，~2×）"
```

---

### Task 9: 物理合理性闸门（"符合物理实际，不只是代码跑通"）

**Files:**
- Create: `cemdisp/diagnostics/physical_sanity.py`
- Create: `scripts/entrypoints/verify_physical_sanity.py`
- Test: `tests/contract/test_physical_sanity.py`（新建）

**Interfaces:**
- Consumes: `AnnulusSimulationResult`（`cement_snapshots`/`snapshot_times_s`/`summary`）、`geom`
- Produces: `sanity_checks(result, geom) -> list[dict]`，每项 `{检查项, 通过, 实测值, 判据, 说明}`

**物理判据（每条都必须可证伪、并写明来源）**：

| # | 检查项 | 判据 | 依据 |
|---|---|---|---|
| P-1 | 浓度有界 | `0 ≤ c ≤ 1`（含全部快照） | 体积分数定义 |
| P-2 | 前缘单调 | 每列 b 加权浓度随时间**单调不减**（容许 1e-9 数值回退） | 纯平流 + 无源项 |
| P-3 | 速度上界 | `max\|w\| ≤ 1.5·Q/A_min`（间隙平均速度的偏心放大有界） | 半环空单位通量 BC + 面积守恒 |
| P-4 | 总通量守恒 | 每列 `∫2H·w̄ dφ = const`（相对差 ≤1e-9） | (2.2) 单位通量归一 |
| P-5 | 窄边劣势 | 密度稳定井（b>0）窄 1/4 行效率 ≤ 宽 1/4 行效率 | 偏心环空顶替基本物理 |
| P-6 | 冻结区静止 | `wall>0.5` 处 `\|w\| ≤ 0.01·max\|w\|`（仅 `enable_stream_yield_gate=True` 时检查） | Pelipenko04 (2.6)-(2.8) |
| P-7 | 排量响应方向 | 排量 ×1.4 时域内 η_E 不下降超过 1 pp（趋势方向） | 现场经验方向（**方向检查，不设幅值门**） |

- [ ] **Step 1: Write the failing test**

```python
# tests/contract/test_physical_sanity.py
import numpy as np

from cemdisp.diagnostics.physical_sanity import sanity_checks


class _Res:
    def __init__(self, snaps):
        self.cement_snapshots = tuple(snaps)
        self.snapshot_times_s = tuple(float(i) for i in range(len(snaps)))
        self.summary = {}


def test_bounded_and_monotone_pass_on_valid_field():
    a = np.zeros((4, 3)); a[:, :1] = 1.0
    b = np.zeros((4, 3)); b[:, :2] = 1.0
    geom = {"b": np.ones((4, 3)), "H": np.full((4, 3), 0.01),
            "phi": np.linspace(0, 1, 4), "s": np.linspace(0, 2, 3),
            "y": np.linspace(0, 0.3, 4)}
    rows = {r["检查项"]: r for r in sanity_checks(_Res([a, b]), geom)}
    assert rows["P-1 浓度有界"]["通过"] is True
    assert rows["P-2 前缘单调"]["通过"] is True


def test_monotone_fails_on_nonphysical_field():
    a = np.zeros((4, 3)); a[:, :2] = 1.0
    b = np.zeros((4, 3)); b[:, :1] = 1.0        # 浓度"倒退"⇒ 非物理
    geom = {"b": np.ones((4, 3)), "H": np.full((4, 3), 0.01),
            "phi": np.linspace(0, 1, 4), "s": np.linspace(0, 2, 3),
            "y": np.linspace(0, 0.3, 4)}
    rows = {r["检查项"]: r for r in sanity_checks(_Res([a, b]), geom)}
    assert rows["P-2 前缘单调"]["通过"] is False
```

- [ ] **Step 2: Run test to verify it fails** → FAIL（模块不存在）

- [ ] **Step 3: Write minimal implementation**

`physical_sanity.py`：按上表实现 7 项检查，返回逐项结果（未提供判据所需输入时该项 `通过=None` 且 `说明="未测：<原因>"`，**不得**默认为通过）。`scripts/entrypoints/verify_physical_sanity.py`：对可跑井跑短窗（`nz=60`）并写 `results/内部自洽加固_2026-09-25/物理合理性台账.csv`。

- [ ] **Step 4: Run test + 生成台账**

Run: `... -m pytest tests/contract/test_physical_sanity.py -v` → PASS
Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 D:/apps/Anaconda/envs/shenjingwangluo/python.exe scripts/entrypoints/verify_physical_sanity.py`
Expected: 生成 CSV；**任何一项不通过都必须如实记录并在结论段写明**（这是诊断输出，不是可调参数）。

- [ ] **Step 5: Commit**

```bash
git add cemdisp/diagnostics/physical_sanity.py scripts/entrypoints/verify_physical_sanity.py tests/contract/test_physical_sanity.py "results/内部自洽加固_2026-09-25/物理合理性台账.csv"
git commit -m "feat(consistency): 物理合理性闸门（7 项可证伪判据）+ 台账"
```

---

### Task 10: 端到端默认路径逐位锚（由 Task 1 评审 R10 追认）

**背景（计划缺陷修正）**：原计划四处引用 `tests/contract/test_six_well_integration.py` 作为"逐位锚定默认路径 ⇒ 证明零数值影响"。**评审已核实该测试全文仅为同步画像卡冒烟测试、零数值断言**，且全仓不存在端到端默认路径的逐位锚。⇒ Task 2/4 的"零数值影响"证明在此之前**无证据可依**，故先建锚。

**Files:**
- Create: `tests/contract/test_default_path_bitwise_anchor.py`
- Create: `tests/contract/_default_path_anchor_hu101.json`（锚文件，首次运行生成后由人工确认并提交）

**Interfaces:**
- Produces: 可被 Task 2 / Task 4 / Task 8 复用的"默认路径逐位不变"证明；`_run_default_case() -> dict` 供测试与后续任务调用

- [ ] **Step 1: Write the failing test**

```python
# tests/contract/test_default_path_bitwise_anchor.py
"""默认路径端到端逐位锚（R10）。

口径：hu101 生产 loader + T1 生产 1D 三开关 + 默认环空开关（不含 CORRECTED_KW），
微网格 nz=30/ny=12/total_t=200s（只求逐位稳定，不求数值精度）。
锁：cement/lead/tail/spacer/wall 五场 sha256 + η_E/η_N/窜槽/混浆/失稳五个标量。
用途：任何声称"零数值影响"的改动都必须让本测试保持全绿。
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from cemdisp.data.loaders import hu101_loader
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver

ANCHOR = Path(__file__).with_name("_default_path_anchor_hu101.json")


def _run_default_case() -> dict:
    well, fluids, schedule, _ = hu101_loader.load_hu101_tailpipe()
    cr = CasingFlowSolver(enable_gravity=True, mixing_contact_time=True,
                          plug_face_zero_mixing=True, has_plug=True).run(
        well, fluids, schedule)
    inlet = build_coupled_annulus_inlet_provider(
        cr, CasingFlowSolver(enable_gravity=True), fluids, split_cement_phases=True)
    solver = AnnulusD2DGASolver(total_t=200.0, nz=30, ny=12)
    res = solver.run(well, fluids, inlet, schedule=schedule)
    out = {}
    for name in ("cement", "lead", "tail", "spacer", "wall"):
        arr = np.asarray(getattr(res, f"{name}_final"), dtype=float)
        out[f"sha_{name}"] = hashlib.sha256(arr.tobytes()).hexdigest()
    fr = res.summary["最终结果"]
    for key, field in (("eta_E", "全井段最终有效顶替效率"),
                       ("eta_N", "窄四分位效率"),
                       ("channeling", "最终窜槽指数"),
                       ("mixing", "最终混浆指数"),
                       ("instability", "最终失稳指数")):
        out[key] = float(fr[field])
    return out


def test_default_path_matches_bitwise_anchor():
    got = _run_default_case()
    if not ANCHOR.exists():
        pytest.skip(f"锚文件不存在，已输出实测供确认：{json.dumps(got)[:400]}")
    want = json.loads(ANCHOR.read_text(encoding="utf-8"))
    assert got == want, "默认路径数值发生位移（若为有意改动，须走重锚并记录位移台账）"
```

- [ ] **Step 2: Run it, confirm the skip + capture actual values**

Run: `... -m pytest tests/contract/test_default_path_bitwise_anchor.py -v -s`
Expected: `SKIPPED`，输出含实测五场 sha256 与五个标量

- [ ] **Step 3: 生成锚文件并人工确认**

把 Step 2 输出的 JSON 写入 `tests/contract/_default_path_anchor_hu101.json`（格式化、UTF-8），**再跑两次** `_run_default_case()` 确认两次结果**完全一致**（同机同版本下逐位可复现）。若两次不一致 ⇒ **停下上报**（说明默认路径本身不确定，"逐位锚"路线不成立）。

- [ ] **Step 4: Run test to verify it passes**

Run: `... -m pytest tests/contract/test_default_path_bitwise_anchor.py -v` → PASS

- [ ] **Step 5: Commit**

```bash
git add tests/contract/test_default_path_bitwise_anchor.py tests/contract/_default_path_anchor_hu101.json
git commit -m "test(consistency): 新建默认路径端到端逐位锚（R10 修正计划缺陷）"
```

---

### Task 11: 锚加固（防假绿 + 防误改 + 来源指纹）——由 Task 10 评审 R17/R18/R19 立案

**背景**：Task 10 建成的默认路径逐位锚是后续所有"零数值影响"结论的唯一证据载体，但评审查出两处缺口：① **假绿**——锚文件缺失时 `pytest.skip`，CI 变绿而零验证（`pyproject.toml` 无 `addopts`、全仓无 `conftest.py` 兜底）；② **无键集守卫**——锚 JSON 被手工改成少键时 `got == want` 仍 PASS。**且最危险的方向不是改测试而是改 JSON 数值**（可让任何位移静默通过，不触发任何检查），这在 Task 8 换线性求解器时必然发生位移的情形下尤其致命。

**Files:**
- Modify: `tests/contract/test_default_path_bitwise_anchor.py`
- Modify: `tests/contract/_default_path_anchor_hu101.json`（只加来源注释键，**不动 5 个数值**）
- Create: `tests/contract/test_anchor_integrity.py`

**Interfaces:**
- Produces: `expected_keys() -> set[str]`（锚必须覆盖的键集，供本测试与 Task 8 复用）

- [ ] **Step 1: Write the failing test**

```python
# tests/contract/test_anchor_integrity.py
"""锚完整性守卫（R17/R18）：缺锚必须失败而非跳过；键集必须齐备；必须有来源指纹。"""
import json
from pathlib import Path

import pytest

from tests.contract.test_default_path_bitwise_anchor import (
    ANCHOR, expected_keys, _FIELD_ATTRS, _SCALAR_FIELDS)

FINGERPRINT_KEYS = {"_note", "_generated_from", "_env"}


def test_anchor_file_must_exist():
    assert ANCHOR.is_file(), (
        f"逐位锚文件缺失：{ANCHOR} —— 缺锚会让默认路径的零数值影响证明静默失效，"
        "严禁用 skip 掩盖")


def test_anchor_covers_exactly_expected_keys():
    want = json.loads(ANCHOR.read_text(encoding="utf-8"))
    keys = set(want) - FINGERPRINT_KEYS
    assert keys == expected_keys(), (
        f"锚键集不匹配：缺 {expected_keys() - keys}；多 {keys - expected_keys()}")


def test_expected_keys_derived_not_hardcoded():
    """键集必须由字段表推导，避免两处各写一份而漂移。"""
    assert expected_keys() == {f"sha_{n}" for n in _FIELD_ATTRS} | set(_SCALAR_FIELDS)


def test_anchor_carries_provenance_note():
    want = json.loads(ANCHOR.read_text(encoding="utf-8"))
    assert "_note" in want and "_generated_from" in want and "_env" in want, (
        "锚必须带来源指纹：微网格非生产数字的警告 + 生成条件 + 环境版本")


def test_anchor_run_does_not_skip():
    """锚存在时主测试不得走 skip 分支（用 skip 计数断言）。"""
    src = Path(ANCHOR.with_name("test_default_path_bitwise_anchor.py")).read_text(
        encoding="utf-8")
    assert "pytest.skip" not in src, "主测试仍可在缺锚时静默跳过（R17 假绿路径）"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONIOENCODING=utf-8 PYTHONUTF8=1 D:/apps/Anaconda/envs/shenjingwangluo/python.exe -m pytest tests/contract/test_anchor_integrity.py -v`
Expected: 至少 2 条 FAIL（`expected_keys` 不存在 / 锚无 `_note`），`test_anchor_run_does_not_skip` FAIL

- [ ] **Step 3: Write minimal implementation**

① 在 `test_default_path_bitwise_anchor.py` 内把字段表提升为模块常量并新增推导函数：

```python
_FIELD_ATTRS = {"cement": "cement_field", "lead": "lead_field", "tail": "tail_field",
                "spacer": "spacer_field", "wall": "wall_field"}
_SCALAR_FIELDS = ("eta_E", "eta_N", "channeling", "mixing", "instability")


def expected_keys() -> set[str]:
    """锚必须覆盖的键集（由字段表推导，禁止与断言行各写一份）。"""
    return {f"sha_{n}" for n in _FIELD_ATTRS} | set(_SCALAR_FIELDS)
```

② 把缺锚分支从 `pytest.skip` 改为**显式失败**：

```python
    want = json.loads(ANCHOR.read_text(encoding="utf-8"))
    assert set(want) >= expected_keys(), (
        f"锚键集不齐：缺 {expected_keys() - set(want)}（R18）")
    got = _run_default_case()
    mismatched = {k for k in expected_keys() if got.get(k) != want.get(k)}
    assert not mismatched, (
        f"默认路径数值位移：{sorted(mismatched)}（若为有意改动，须走重锚："
        "与改动同批次、并在 results/内部自洽加固_2026-09-25/位移台账.csv 留迹）")
```

③ 锚 JSON 增加来源指纹（**不动任何数值**）：

```json
{
  "_note": "微网格锚（nz=30/ny=12/total_t=14000s），非生产数字，不得用于论文或与现场 CBL 比对",
  "_generated_from": "HEAD bc155b0 + 4 个未提交 HB 改动（工作树脏态）；生成命令见 task-10-report.md",
  "_env": "Python 3.13.7 / numpy 2.3.3 / scipy 1.16.2（对账源：docs/superpowers/plans/baseline-2026-09-25.md）",
  "...": "原有 5 个 sha_* 与 5 个标量逐字保留"
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `... -m pytest tests/contract/test_anchor_integrity.py tests/contract/test_default_path_bitwise_anchor.py -v` → 全 PASS
Run（哨兵验证）：临时把锚 JSON 里 `eta_N` 改 1e-9 ⇒ `test_default_path_bitwise_anchor` 必须 FAIL；恢复 byte-identical

- [ ] **Step 5: Commit**

```bash
git add tests/contract/test_default_path_bitwise_anchor.py tests/contract/_default_path_anchor_hu101.json tests/contract/test_anchor_integrity.py
git commit -m "test(consistency): 锚加固——缺锚即失败/键集守卫/来源指纹（R17-R19）"
```
