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

---

### Task 12: 管容常量口径落地（用户裁定 R67/R68）

**Files:**
- Modify: `cemdisp/data/loaders/ht1_001_loader.py`（常量值 + 注释链）
- Modify: `cemdisp/data/loaders/hu2_loader.py`（分段链纠正 + 注释链）
- Modify: `tests/contract/test_shoe_lag_wiring.py`（等值断言随新值更新）
- Test: `tests/contract/test_shoe_lag_wiring.py`

**Interfaces:** Consumes 证据报告 `docs/superpowers/research/2026-09-25-shoe-lag-volume-caliber-evidence.md`（含一手资料出处与逐字摘录）

- [ ] **Step 1: ht1_001 常量改为 94.439**，并在注释链写清推导：`93.4934（7.1.4 表四段和，止于阻位 7642.674m）+ 0.95（7.1.1 尾浆(下塞)行，阻位→鞋 103.326m）= 94.439`
- [ ] **Step 2: hu2 三段偏离纠正**（按设计 7.1.4 表口径）：149.2 钻杆改用设计表的 `12.91 L/m`（≈ID128.2）；补 `5292.5→5308.564` 的 16.064 m；四段和算至**阻位 7438.9**。**Step 2 的产出 = 两个数值**：`V_to_stop_collar`（止于阻位）与 `V_to_shoe = V_to_stop_collar + (阻位→鞋段容积)`
- [ ] **Step 3: 取 `V_to_shoe` 作为 `HU2_SHOE_LAG_VOLUME_M3`**（与 R67 的"严格地面→鞋"一致），注释里**并列写明 `V_to_stop_collar` 的值**并注明"若改为阻位口径，只需把常量换成该值"
- [ ] **Step 4: 更新 `tests/contract/test_shoe_lag_wiring.py` 的等值断言**为新值（`pytest.approx(rel=1e-9)`），并保留 R62 的 1.05 量级哨兵与交叉判别测试
- [ ] **Step 5: 更正注释链里的段落号错标**（`6.4.1` → `7.1.4`，并注明 `6.4.1` 实为"套管基础数据"）；**不要**改 `pumping_schedule.csv`/`well_geometry.csv`（另案）
- [ ] **Step 6: 原始输出**：两口井新旧值对照打印 + 本任务测试绿 + contract 仅 1 条既存红
- [ ] **Step 7: Commit**（**不**执行 Task 3 的 Step 5；不写 `results/`）

---

### Task 13: 停止时刻逻辑（用户裁定 R69；**实施前须先出方案报用户**）

**背景**：`casing_flow.py:1033`/`:1109` 的到达判据 `target = 步累计体积 + 管容` 超出末步累计即 `return None`，使 `annulus_stop_time_s()` 回退为"泵注结束"。证据显示该边界**零余量**（hu2 边界 79.0、ht1_001 边界 93.7）⇒ 结构性脆弱。用户裁定"修停止时刻逻辑"。

**Files:** 待方案确定（预计 `cemdisp/transport1d/casing_flow.py`、各 runner 的 `annulus_stop_time_s` 调用点）

- [ ] **Step 1: 出方案报用户**（至少三个候选：① 显式工艺停泵时刻 + 尾缘迟到量；② 停泵时刻 = 泵注结束 + 尾缘到鞋的预期滞后（解析外推）；③ 保留判据但在回退时强制告警并写入摘要），每条附：物理依据、对 8 井的波及面、实现量、是否改数值口径
- [ ] **Step 2: 用户选定后再拆实现步骤**

---

### Task 14: dict→splat 调用面治理（由 Task 5 评审 R98 立案）

**背景**：Task 5 建成的调用签名闸门只抓**直接关键字调用**；而**至少 9 个脚本**用 `dict(...)` + `**kw`/`extra_kw` 把已删形参 `dispersion_dt_scale` 传进求解器 ⇒ **运行期仍 `TypeError`，闸门却报绿**。闸门"绿而无用"与本计划的目的直接冲突。

**Files:**
- Modify: `scripts/entrypoints/check_call_signatures.py`（扩展 dict→splat 数据流校验）
- Modify/Delete: 涉及的 9 个脚本（按 `scripts/entrypoints/rerun_all_wells_corrected.py` 为样板清理）
- Test: `tests/contract/test_entrypoint_signatures.py`（补相应断言）

**Interfaces:** Consumes `find_bad_kwargs(root) -> list[tuple[str,int,str,str]]`（Task 5 产物）

- [ ] **Step 1: 先出清单**：用 ast 或 grep 列出全部 `dict(...)`+`**` 进入 `AnnulusD2DGASolver`/`CasingFlowSolver` 的调用点，并逐个判定其 dict 中是否有**已删形参**；把清单写进报告（数量以实测为准，不采信"9"这个数）
- [ ] **Step 2: 清理**：从各 dict 中删掉已删形参键（**不得**重新引入形参）；若某脚本已整体作废（如 `dispersion_scale_sensitivity_scan.py`，R102 已记"8 次仿真逐位相同"），**优先归档/删除**并说明
- [ ] **Step 3: 扩展闸门**：在 `check_call_signatures.py` 内增加**同函数作用域**的 `name = {...}` → `Call(**name)` 键集校验（够用即可，不必做跨函数的全量数据流）
- [ ] **Step 4: 测试钉死**：加一条"坏 dict + splat 必须被抓到"的用例（构造临时文件）
- [ ] **Step 5: 判据**：本任务测试绿 + 锚两文件绿 + 既存红不新增；**按 R104 用 pathspec 提交**

---

## 附录 A：影响论文的结构性发现（受版本控制留档，2026-09-25）

### A.1 `R0` 与 `R1` 消融级结构性恒等 ⇒ 不得再生产 R0/R1 差异类图表

**事实（Task 5 评审三条独立证据）**：
1. 级别定义仅 `enable_d2dga_auto_m` 一项不同（`cemdisp/runners/ht1_004_ablation.py:47-52` 的 `ABLATION_LEVELS`；`scripts/entrypoints/closure_contribution_scan.py:26-32` 的 `LEVELS` 同理），而该形参**已自 `d8917b7`（2026-09-08，"R0 兼容分支删除——auto-m 恒开"）起从 `__init__` 移除** ⇒ 两者构造出的 kwargs 逐字相同。
2. 求解器确定性（无 RNG，CFL 自适应为确定式）⇒ 同 kwargs ⇒ 输出**逐位相同**。
3. `inspect.signature` 实测 `has VAR_KEYWORD=False`（无 `**kwargs`）⇒ 删参前那 4 处调用**确实**抛 `TypeError`；即两个消融脚本自 09-08 起**已 17 天不可调用**。

**受害物（规划层）**：`docs/superpowers/plans/2026-07-13-improved-d2dga-paper-design.md` §Task 3 规划「**图6（R0 vs R1 窄边窜槽对比）**」并写明 "Expected: R0 与 R1 效率有差异"。该图现**结构性恒平（差值 ≡ 0）**。

**约束（写作红线）**：
- **不得**再生产任何"R0 vs R1 有差异"的图表或表格。
- 若论文/汇报材料已引用该类对比，须核实其来源：要么**早于 2026-09-08**（当时 auto_m 仍是有效开关），要么是**退化产物**（差值恒 0），后者必须撤除。
- `closure_contribution_scan.py` 输出的 CSV 仍保留 `auto_m` 列（R0=0/R1=1），而指标列将逐位相同 ⇒ 该列是**纯标签、零信息**，却看起来像做过 auto_m 消融。引用该表时须显式说明。

**状态**：**未处理**（需先裁定：删除 R0 级、或用真实区分量替代 R0/R1 的对照维度）。已并入 Task 14 的范围讨论。

---

## 附录 B：退化对照总表（Task 14 清理轮留档，2026-09-26，受版本控制）

### B.0 这份表解决什么问题

附录 A 留档了 **R0≡R1 结构性恒等**（`enable_d2dga_auto_m` 形参已删 ⇒ 两级 kwargs 逐字相同）。
Task 14 的清理又暴露**同一类问题的另一批实例**：一批"看起来做过消融/敏感性对照、实际差值恒 0"
的行。这类行**不能当对照读**——引用它会产出"某机制无影响"的假结论。

本附录给出**逐行可审计**的清单：变体名 / 脚本 / 它等于哪一行 / **判定的方法** / 失去意义的图或表。
审计论文表格的人不需要重做本文的工作：只要拿本表的"变体名 + 脚本"去对账即可。

### B.0.1 两种方法，不得混用

| 方法 | 含义 | 需要的依据 |
|---|---|---|
| **结构性** | 两行的 kwargs **逐字相同**（或差异仅在"等于形参默认值"的位上）⇒ 同一求解器 + 确定性 ⇒ 输出逐位相同 | 读代码给出 kwargs 逐字对比；"等于默认值"须给出 `inspect.signature` 的默认值 |
| **实测** | kwargs **不同**，但**跑出来逐位相同** | 给出实际跑的算例、比对口径与结果（哈希/逐位断言） |

**两条前提**（下面所有"结构性"行的公共前提）：

1. **求解器确定性（实测）**：`tests/contract/test_default_path_bitwise_anchor.py` 的冻结锚
   是在**另一个提交/进程**里生成的，当下运行仍**逐位复现**（该测试长期全绿，且不做浮点容差）。
   本轮另做同进程双跑对照：同 kwargs 两次构造运行，五场 sha256 + 五标量**完全一致**
   （见 B.0.2 的 A 行与 F 行）。
2. **`has VAR_KEYWORD = False`**：`AnnulusD2DGASolver.__init__` 无 `**kwargs`
   （闸门 `_var_keyword_targets()` 每次运行都实测为空）。因此"已删形参"只会触发 `TypeError`，
   不会被静默吸收。

### B.0.2 实测：默认路径上哪些开关真的改变输出

口径（与冻结锚同款，**微网格、非生产数字**）：hu101 生产 loader + T1 生产 1D 三开关 +
**默认环空配置**（不含 `CORRECTED_KW`），`nz=30 / ny=12 / total_t=14000s`（覆盖整个顶替序列）。
指纹 = 五场（cement/lead/tail/spacer/wall）sha256 + 五个标量（η_E/η_N/窜槽/混浆/失稳）拼接后取前 16 位。

| kwargs | 指纹 | 与基线 |
|---|---|---|
| `{}`（基线） | `a24abcfffac3b0a7` | — |
| `{}` 重跑 | `a24abcfffac3b0a7` | **逐位相同**（确定性对照） |
| `enable_regime_split=True` | `a24abcfffac3b0a7` | **逐位相同 ⇒ 该开关在默认路径上无消费者** |
| `enable_true_buoyancy=False` | `a24abcfffac3b0a7` | **逐位相同 ⇒ 同上** |
| `enable_power_law_gap_law=False` | `a24abcfffac3b0a7` | **逐位相同 ⇒ 同上** |
| `e_clip_max=0.90` | `a24abcfffac3b0a7` | **逐位相同 ⇒ 同上（已弃用）** |
| `enable_e_clip_ruling=False` | `a24abcfffac3b0a7` | **逐位相同 ⇒ 同上（已弃用）** |
| `enable_local_i3=True` | `4972b49f78c2e756` | 不同（活） |
| `yield_gate_f_safety=1.6` | `b97836ac3f47ea3c` | 不同（活） |
| `yield_regularization_M=1000.0` | `50b93341c7754516` | 不同（活） |
| `enable_yield_gate=False` | `07e91da61a97401e` | 不同（活） |
| `enable_stream_yield_gate=True` | `c78a6a763d2a9e50` | 不同（活） |
| `enable_d2dga_i3_flux=False`（对照） | `71a0c623c6b9e5ae` | 不同（**证明比对有牙**） |
| `enable_d2dga=False`（对照） | `71a0c623c6b9e5ae` | 不同；**与上一行同指纹**——两者闸的是同一个 I3 通量块（与 `_dead_switches` 的注释一致） |
| `enable_power_law_gap_correction=True`（对照） | `cd10aaf804decc41` | 不同（对照） |

**形参默认值（`inspect.signature` 实测）**：`enable_yield_gate=True`、`enable_regime_split=False`、
`enable_local_i3=False`、`enable_true_buoyancy=True`、`e_clip_max=0.55`、`enable_e_clip_ruling=True`。
⇒ **显式传"等于默认值"的开关与不传等价**（Python 语义，确定性前提同上）。

**代码侧独立印证（结构性）**：`cemdisp/models2d/annulus_d2dga.py` 的 `_dead_switches`
（A3 惯例，R42）与 `_OLD_PATH_ONLY_SWITCHES` 已**在代码里声明**：`enable_regime_split` /
`enable_true_buoyancy` / `enable_power_law_gap_law` 属"仅旧代数路径
（`enable_stream_function=False`）消费"，默认路径偏离即告警；`e_clip_*` 三形参由构造器
DeprecationWarning 声明"显式传值不再生效"。**实测与代码声明一致**，故下文的"实测"行
既可由微网格复现，也有代码侧的路径级依据（路径属性与网格无关，故可外推到脚本的 nz=80/250）。

### B.1 本次清理**直接造成**的退化行（一行一条）

删键前这些行会 `TypeError`（脚本根本跑不起来）；删键后它们能跑，但**与某一行逐位相同**。

| # | 变体名 | 脚本 | 等于哪一行 | 方法（依据） | 失去意义的图/表 |
|---|---|---|---|---|---|
| B1-1 | `dispersion_zero` | `scripts/entrypoints/run_ablation_variants.py` | **同脚本的 `i3_localized`**（两者只差一个 `enable_local_i3`；见下方 kwarg 并排） — ⚠️ **不等于**表内任何"基线"行，也不等于"求解器全默认"（见"方法"栏的边界说明） | **结构性（限于"不再覆盖任何参数"这一条）**：删键后 `VARIANTS["dispersion_zero"]` 的 dict 为空，构造实参只剩字面量 `total_t`/`nz`（`:104-105`）⇒ 与同脚本 `i3_localized` 相比**仅少 `enable_local_i3`**。<br>**未实测**：本行**不**声明"与表内某基线行逐位相同" —— 表内基线的数值来自 `_read_baseline()` 读的**权威摘要 JSON**（`:99`），那是**另一个 runner**（`cemdisp/runners/*`，带它自己的 `CORRECTED_KW`）写的，kwargs 与 `total_t`/`nz` 都无法在本表内并排核对 ⇒ 按 B.0.1 的要求，该断言在此**不成立、不作声明** | 该脚本自产的 `results/消融变体_runner口径_2026-09-09/`（CSV/MD/逐变体 JSON）里 `dispersion_zero` 行：**与同表"求解器全默认"行互为重复**（同脚本内两者 kwargs 相同）；按上栏的边界，**不得**把它读成"与权威基线无差别" |
| B1-2 | `M1 only` | `scripts/probes/isolate_fix_mechanisms.py` | `BASELINE(全关)` | **结构性**：两行 kwargs 均为 `{}`，逐字相同 | 09-02 呼101 机制归因表；**M1 独立归因失效**（对应 2026-08-23 计划 Task 5 / 设计稿 §4 的"M1 一轮独立归因"） |
| B1-3 | `M1+M3` | 同上 | `M3 only` | **结构性**：两行 kwargs 均为 `{enable_yield_gate: True}`，逐字相同 | 同上 |
| B1-4 | `ALL corrected` | 同上 | **不再与同表任何行重合** | — | 我 Task 14 首轮报告把 `ALL corrected` 也列为退化行，**此处更正：该判断过强**。它只是少了已删的 M1 维，仍与 `M3+M2+M4` 不同（多 `enable_local_i3`），标签与内容仍相符 |
| B1-5 | `b5_dispaz0.002` | `scripts/probes/hu101_standoff_response_tuning_survey_20260911.py` | **同 SO 的"无旋钮"配置**（A 段在 `SO=0.55` 有同名行；`SO=0.40/0.70` 无 A 段对应行，但其 B 段斜率即无旋钮基线斜率） | **结构性**：`solver_kw={}` 与 `solver_kw=None` 都归约为 `{}`；同 `well_spec`/`fluids`/`schedule`/缓存键 | 该脚本 B 段"旋钮对响应斜率的控制"表的 `b5` 行；**无已知论文图/表受害者**（09-11 调研报告的"方位弥散空杠杆"结论**方向仍对**，但今天不再由该脚本的对照给出） |
| B1-6 | `D_disp_x0.5` / `D_disp_x2` / `D_disp_x3` | `scripts/probes/hu101_low_score_attribution_probe_20260911.py` | `A_基线` | **结构性**：同 `well_spec`/`fluids`/`schedule`/缓存键，`solver_kw={}` 与缺省等价 | 该脚本 D 段表（3 行 = `A_基线` 的 3 份副本）；无已知论文图/表受害者 |
| B1-7 | `scale ∈ {0.0, 0.5}` 全部档 | `scripts/probes/dispersion_scale_sensitivity_nz250.py`（**已归档**） | 彼此逐位相同；与外部基线 JSON 档（`scale=1.0`，`source=baseline_reuse`）的关系**未实测** | **结构性**（彼此：kwargs 逐字相同）。**未实测**：脚本断点续跑会复用旧 JSON；"与外部 JSON 也相同"只有脚本自述的同源口径，我没有跑 | 2026-08-23 计划 Task 5 的 κ 扫描验收（①mixing 0.59→0.2–0.35 ②scale→0 平台由 κ 主导 ③前沿米级）；已归档 |
| B1-8 | 全部 8 档（`scale∈{0,0.25,0.5,1.0}` × CFL on/off） | `scripts/probes/dispersion_scale_sensitivity_scan.py`（**已归档**） | 彼此逐位相同 | **实测（R102 记录，非我本次实测）**：`.superpowers/sdd/2026-09-25-internal-consistency-hardening/progress.md:268` 记"8 次仿真输出逐位相同"；本次未复跑（脚本已归档） | 同上（M1 κ 扫描验收）；已归档 |
| B1-9 | `关壁面冻结` | `scripts/probes/_wallfreeze_grid_mechanism_20260902.py` | `BASE` | **结构性**：两行 kwargs 均为 `dict()`；`run_variant(name, ...)` 的 `name` 只进输出标签，不进构造 | `results/_质量平衡取证_2026-09-02/阶段5_壁面冻结网格机制.json` 的两行；B2 屈服门取证链 |
| B1-10 | `关壁面冻结+关弥散` | 同上 | `BASE`（**经 B1-9 传递** —— 直接的同一行是 `关壁面冻结`） | **结构性**：与 `关壁面冻结` 同为 `dict()`（经 B1-9 传递 ⇒ 亦等于 `BASE`）；`run_variant(name, ...)` 的 `name` 只进输出标签，不进构造。⚠️ 本行的**直接**同一行是 B1-9，不是 `BASE` ——"等于 `BASE`"只是传递结论，审计时应先对 B1-9 | 同上 |

**B1-1 的 kwarg 并排**（唯一可当场核对的比对；均引自 `run_ablation_variants.py`）：

```python
# :104-105 —— 两个变体走同一条构造语句，只有 **overrides 不同
solver = AnnulusD2DGASolver(total_t=total_t_s, nz=250, **overrides)

# :40  i3_localized     → {"enable_local_i3": True}
# :42  m2_regime_split  → {"enable_regime_split": True}
# :43  dispersion_zero  → {}            ← 删键后为空：total_t/nz 之外不再传任何东西
```

⇒ `dispersion_zero` 与 `i3_localized` 的差别**只有一个 `enable_local_i3`**；与"求解器全默认"的差别只有字面量 `total_t`/`nz`（这两项是脚本的时间/网格口径，不是开关）。
**表内基线的数值不可用于本行比对**：它由 `_read_baseline()`（`:99`）从权威摘要 JSON 读出，写它的是 `cemdisp/runners/*`（带自己的 `CORRECTED_KW`），kwargs 与 `total_t`/`nz` 都不公开在本表内。

| B1-11 | `对照_完全无壁面层` | `scripts/probes/_yieldgate_verify_20260902.py` | `对照_浓度冻结基线` | **结构性**：两行 kwargs 均为 `dict()`，逐字相同 | `results/_质量平衡取证_2026-09-02/阶段6_物理屈服门验证.json` 的两行 |

### B.2 同一批脚本里**本次清理之外**发现的同类退化（另行发现，非本次清理引入）

这些行的 kwargs **并未**被本次清理改动；它们的退化来自"死开关 / 显式传默认值"，
依据是 B.0.2 的**实测**与代码侧 `_dead_switches` 声明。**按行类列**（成员逐条列在"行成员"里，
以便逐行对账）。

| 类 | 退化的开关（依据） | 行成员（脚本 · 变体名） | 失去意义的图/表 |
|---|---|---|---|
| B2-1 | `enable_regime_split=True` 在默认路径无消费者（**实测**，`_OLD_PATH_ONLY_SWITCHES` 声明） | `run_ablation_variants` · `m2_regime_split`；`isolate_fix_mechanisms` · `M2 only`、`M3+M2+M4`（该位）、`ALL corrected`（该位）；`CORRECTED_KW` 的该位（`bisect_hu103_20260902` / `c_verify_convergence_20260902` / `corrected_ref_hu1_hu103_20260902` / `debug_hu1_hu103_eta_zero_20260902` / `dump_tailwindow_2d_v1` / `rerun_stop_fix_20260901` / `analyze_distortion_fix` / `sensitivity_common` / `mass_balance_diag`） | `run_ablation_variants` 消融表 `m2` 行（与 `dispersion_zero` 互为重复）；09-02 各脚本的 `+enable_regime_split` 逐开关行 |
| B2-2 | `enable_true_buoyancy=False` 在默认路径无消费者（**实测**，同上声明） | `hu101_low_score_attribution_probe_20260911` · `B_真浮力关`；`run_gap_fill_variants` · `r2_i3`（唯一键）、`r0_base`（该位）；`_wallfreeze_grid_mechanism_20260902` · `关壁面冻结+关弥散+关D2DGA`（该位） | 09-11 调研报告的"真浮力关 +15.6pp"类结论（**注**：该数字来自**旧代数路径**口径；此处只声明"在**当前 HEAD** 上该脚本无法复现该对照"，不否认历史口径的成立） |
| B2-3 | `enable_power_law_gap_law=False` 在默认路径无消费者（**实测**） | `hu101_standoff_response_tuning_survey_20260911` · `b2_gaplaw_off` | 该脚本 B 段旋钮表 `b2` 行（5 个旋钮里 2 个是假对照：`b1`/`b2`） |
| B2-4 | `e_clip_max` / `e_clip_measured_max` / `enable_e_clip_ruling` 已弃用、显式传值不再生效（构造器 DeprecationWarning 声明 + **实测**） | `hu101_standoff_response_tuning_survey_20260911` · `b1_eclip0.90`；`hu101_low_score_attribution_probe_20260911` · 整个 `C` 段（5 档 `C_eclip=…`）；`isolate_fix_mechanisms` · `M4 only(e=.90)`、`M3+M4`、`M3+M2+M4`、`ALL corrected`（该位）；`_wallfreeze_grid_mechanism_20260902` · `关壁面冻结+关弥散+近同心`；`CORRECTED_KW` 的该位（同 B2-1 的脚本清单） | 该两脚本的旋钮/扫描表中 `b1` 与整个 `C` 段；`isolate_fix_mechanisms` 的 `M4` 归因 |
| B2-5 | 显式传"等于形参默认值"的开关 ≡ 不传（`inspect.signature` **实测**默认值：`enable_yield_gate=True`） | `isolate_fix_mechanisms` · `M3 only`、`M1+M3`、`M3+M4`、`M3+M2+M4`；`_yieldgate_verify_20260902` · `物理屈服门` | `isolate_fix_mechanisms` 的 `M3` 归因；`_yieldgate_verify_20260902` 的**整张表**（3 行全部等价 ⇒ 零信息） |
| B2-6 | `dispersion_axial` / `dispersion_azimuthal` 已删（Task 7 删除自创拉普拉斯弥散；删前 `dispersion=False` 分支传的也是**默认值 0.0**） | `scripts/probes/_pure_plug_conservation_20260902.py` · `纯水泥_0.7_同心_无弥散`（`阶段3_纯水泥守恒实验.json` **第 1 行**）vs `纯水泥_0.7_同心_有弥散`（**第 5 行**） | 该脚本 docstring 明写"无弥散 vs 默认弥散"是要扫描的维度之一，`:94` 也确实是唯一的 `dispersion=True` 调用 —— 但它**没有对照物**：**实测**（本任务修复轮 2，`run_one(..., dispersion=False, nz=10)` vs `dispersion=True`）两次输出字典**除两个标签键外逐字节相同**（九个数值字段全等：守恒率 0.9962、η_E 0.6973、域内水泥 8.63 …）⇒ **第 5 行不是弥散对照**，读它会得出"弥散对纯水泥守恒无影响"的假结论 |

**B2 的合并后果（逐脚本）**：

- `run_ablation_variants.py`：3 个变体里 **2 个**（`m2_regime_split`、`dispersion_zero`）**不再覆盖任何参数**，
  只有 `i3_localized` 触及活开关 ⇒ 该消融表在 M2 与"弥散"两维上零信息。
  ⚠️ 注意口径（见 B1-1 的边界说明）：本脚本的表内"基线"数值来自**另一个 runner** 写的权威摘要 JSON，
  其 kwargs 无法在本表内并排 ⇒ **不得**把 `dispersion_zero` 读成"与权威基线逐位相同"；可核对的只是"同脚本内不再覆盖参数"。
- `isolate_fix_mechanisms.py`：10 行塌成 **2 个等价类** ——
  `{I3 only, ALL corrected}`（活的部分是 `enable_local_i3`）与 `{其余 8 行}`。
- `hu101_low_score_attribution_probe_20260911.py`：**B/C/D 三段共 9 行全部退化为 `A_基线` 副本**
  （B：1 行真浮力关；C：5 档 e_clip；D：3 档弥散），只有 `E_nz500`、`F_*`、`G_*` 仍有区分力。
- `_yieldgate_verify_20260902.py`：整表 3 行等价。
- `hu101_standoff_response_tuning_survey_20260911.py`：B 段 5 旋钮里 `b1`/`b2`/`b5` 三个是假对照，
  只有 `b3_fsafety1.6`、`b4_M1000` 有区分力（两者**实测为活**，与 09-11 旧路径口径下"M1000 ±0.0000 死"的结论**不同**——
  那是旧路径的结论，不得跨路径引用）。

### B.3 状态与需裁定

- **B.1（本次清理造成）**：已按 Task 14 的硬约束处理 —— **只删键、不改口径**，标签保留，
  在代码内加注释披露，并在本表逐行留档。**不需要**再改脚本，除非裁定"删行"。
- **B.2（本次清理之外）**：**未处理**。建议裁定方向（三选一，可组合）：
  ① 删掉无区分力的行；② 用有真实区分力的量替换该维度；③ 保留行但在输出表/图注里显式标注"该行与基线同构、零信息"。
  无论选哪条，都应同时裁定**是否做全仓同类扫描**（本附录只覆盖 Task 14 触碰过的 18 个脚本，
  `enable_true_buoyancy=False` / `e_clip_max` 在**其它**脚本（如 `scripts/entrypoints/run_gap_fill_variants.py`
  以外的跑批、`cemdisp/runners/*`）的用法**未**逐一枚举）。
- **不在本附录判断范围**：2026-08-23 设计稿 §4 列的"表4（8 井效率/窜槽/混浆）、表5（消融）、
  表6（网格收敛）、67% 居中度阈值、呼101 '65.04% vs 62.77% +2.27pp'、'20 秒'卖点"——那是
  **M1 是否落地**对论文数字的影响清单，与"灵敏度扫描有无判别力"是两件事。附录 A 已覆盖 R0/R1
  的图6/图9/表5 受害者，本附录不重复声明。
- **时效口径（重要）**：B.0.2 与 B.2 的"默认值恒等"判定按**当前 HEAD** 的形参默认值给出；
  我**未**追溯历史默认值，因此这些行在**历史运行**时是否同构**不在此表结论内**——
  本表只保证"**今天重跑**时同构"。

### B.4 复核方式（可复现）

```bash
# 1) 闸门：全仓调用面（应 0 处问题）
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/entrypoints/check_call_signatures.py

# 2) 形参默认值（B2-5 的结构性依据）
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -c "import inspect; \
from cemdisp.models2d import AnnulusD2DGASolver as S; \
p=inspect.signature(S.__init__).parameters; \
print({k: p[k].default for k in ('enable_yield_gate','enable_regime_split','enable_local_i3',\
'enable_true_buoyancy','e_clip_max','enable_e_clip_ruling')})"

# 3) B.0.2 的开关实测：按 tests/contract/test_default_path_bitwise_anchor.py 的
#    _run_default_case() 同款流水线（hu101 loader + T1 三开关 + 默认环空配置，
#    nz=30/ny=12/total_t=14000s），对每个 kwarg 组合跑一次，比五场 sha256 + 五标量。
#    对照项（enable_d2dga_i3_flux=False / enable_d2dga=False /
#    enable_power_law_gap_correction=True）必须"不同"，否则说明该算例无判别力。

# 4) 确定性前提
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest \
  tests/contract/test_default_path_bitwise_anchor.py tests/contract/test_anchor_integrity.py -q
```
