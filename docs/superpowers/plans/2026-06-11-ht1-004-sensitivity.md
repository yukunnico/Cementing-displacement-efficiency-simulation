# HT1-004 Sensitivity Script Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `scripts/ht1_004_sensitivity.py` to run 呼1-004 tail-slurry sensitivity analysis with the same parameter ranges and outputs as the existing 呼1-003 script, while using HT1-004 input data and HT1-004 runner stop-time semantics with `nz=140`.

**Architecture:** Implement a standalone script by copying the proven structure of `scripts/ht1_003_sensitivity.py`, changing only well-specific imports, names, output paths, and the annulus stop-time helper. The script stays isolated from existing loaders/runners and does not refactor shared utilities.

**Tech Stack:** Python, pandas, numpy, matplotlib Agg backend, tqdm, existing `cemdisp` 1D/2D solver APIs.

---

## File Structure

- Create: `scripts/ht1_004_sensitivity.py`
  - Standalone HT1-004 sensitivity-analysis CLI.
  - Loads `load_ht1_004_tailpipe()`.
  - Mutates only `尾浆` fluid properties and `注入尾浆` pumping rate per parameter combo.
  - Uses local `_annulus_stop_time_s()` copied semantically from `cemdisp/runners/ht1_004_tailpipe.py`, but keeps `AnnulusD2DGASolver(..., nz=140)`.
- Existing reference only: `scripts/ht1_003_sensitivity.py`
  - Source template; do not modify.
- Existing reference only: `cemdisp/runners/ht1_004_tailpipe.py`
  - Stop-time behavior source; do not modify.
- Existing reference only: `cemdisp/data/loaders/ht1_004_loader.py`
  - Loader source; do not modify.

## Task 1: Create HT1-004 sensitivity script

**Files:**
- Create: `scripts/ht1_004_sensitivity.py`
- Reference: `scripts/ht1_003_sensitivity.py`
- Reference: `cemdisp/runners/ht1_004_tailpipe.py:152-174`

- [ ] **Step 1: Copy the HT1-003 script as the starting point**

Copy `scripts/ht1_003_sensitivity.py` to `scripts/ht1_004_sensitivity.py` using the editor or filesystem copy. Do not modify `scripts/ht1_003_sensitivity.py`.

- [ ] **Step 2: Update module docstring and imports**

In `scripts/ht1_004_sensitivity.py`, replace the header docstring with:

```python
"""
呼1-004井尾浆参数敏感性分析脚本

对尾浆的注入排量、屈服应力、塑性黏度进行全因子敏感性分析，
评估各参数对全井有效顶替效率的影响。

参数空间：
- 排量: 0.30~1.60 m^3/min, 步长0.10, 共14个水平
- 屈服应力: 5~25 Pa, 步长1, 共21个水平
- 塑性黏度: 50~350 mPa·s, 步长20, 共16个水平
- 全因子组合: 14 × 21 × 16 = 4704次模拟

输出：
- CSV汇总表: results/呼1-004_敏感性分析/呼1-004_敏感性分析结果.csv
- 单因素曲线图 ×3
- 两两交互热力图 ×3

使用方式：
    python scripts/ht1_004_sensitivity.py
    python scripts/ht1_004_sensitivity.py --resume  # 断点续跑
    python scripts/ht1_004_sensitivity.py --plot-only  # 仅生成图表
"""
```

Replace the HT1-003 loader import:

```python
from cemdisp.data.loaders.ht1_003_loader import load_ht1_003_tailpipe
```

with:

```python
from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.loaders.ht1_004_loader import load_ht1_004_tailpipe
```

If the copied file already imports `FluidSpec`, update that existing line to import both `FluidRole` and `FluidSpec` instead of adding a duplicate import.

- [ ] **Step 3: Update output paths and strings**

Set these constants exactly:

```python
OUTPUT_DIR = _PROJECT_ROOT / "results" / "呼1-004_敏感性分析"
CSV_PATH = OUTPUT_DIR / "呼1-004_敏感性分析结果.csv"
CHECKPOINT_PATH = OUTPUT_DIR / "呼1-004_敏感性分析_中间结果.csv"
```

Replace all visible user-facing `呼1-003` labels in the copied script with `呼1-004`, including plot filenames, argparse description, banner text, and final output messages.

Keep parameter constants unchanged:

```python
RATE_MIN = 0.30
RATE_MAX = 1.60
RATE_STEP = 0.10
RATE_BASELINE = 1.20
YP_MIN = 5
YP_MAX = 25
YP_STEP = 1
YP_BASELINE = 14.0
PV_MIN = 50
PV_MAX = 350
PV_STEP = 20
PV_BASELINE = 180.0
```

- [ ] **Step 4: Add HT1-004 annulus stop-time helper**

Add this helper near `_run_single()` before it is used:

```python
def _annulus_stop_time_s(casing_result, fluids: tuple[FluidSpec, ...]) -> float:
    """返回呼1-004环空二维顶替应停止的地面累计时间。

    对齐 cemdisp.runners.ht1_004_tailpipe.annulus_stop_time_s：
    当水泥浆之后的首个非水泥流体到达鞋口时停止；若未找到，
    回退到 cement_end_time_s。
    """
    cement_roles = {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}
    fluid_by_name = {fluid.name: fluid for fluid in fluids}
    found_cement = False
    for front in casing_result.fronts:
        fluid = fluid_by_name.get(front.fluid_name)
        if fluid is None:
            continue
        if fluid.role in cement_roles:
            found_cement = True
            continue
        if found_cement:
            return float(front.time_s)
    return float(casing_result.cement_end_time_s)
```

- [ ] **Step 5: Update data loading and single-run stop time**

In `run_sensitivity_analysis()`, replace:

```python
well_spec, base_fluids, base_schedule, _ = load_ht1_003_tailpipe()
```

with:

```python
well_spec, base_fluids, base_schedule, _ = load_ht1_004_tailpipe()
```

In `_run_single()`, replace the copied HT1-003 stop-time block:

```python
# 环空二维顶替（停止时间 = 水泥浆全部进入环空后替浆液到鞋口）
# 敏感性分析关注相对趋势，用 cement_end_time_s 即可，速度快
# 网格分辨率 nz=140：效率与 nz=500 差异<0.4%，速度快3倍
total_t_s = float(casing_result.cement_end_time_s)
annulus_solver = AnnulusD2DGASolver(total_t=total_t_s, nz=140)
```

with:

```python
# 环空二维顶替：停止时间对齐呼1-004 runner 的严格现场口径，
# 即水泥之后的首个非水泥流体到达鞋口时停止；批量敏感性仍用 nz=140 加速。
total_t_s = _annulus_stop_time_s(casing_result, fluids)
annulus_solver = AnnulusD2DGASolver(total_t=total_t_s, nz=140)
```

- [ ] **Step 6: Check all references to HT1-003 are gone from the new script**

Run a text search in `scripts/ht1_004_sensitivity.py` for these strings:

```text
ht1_003
load_ht1_003_tailpipe
呼1-003
HT1-003
```

Expected: no matches, except none at all is preferred. The new script should use only HT1-004 labels and `load_ht1_004_tailpipe`.

## Task 2: Verify syntax and one-combo execution

**Files:**
- Run only: `scripts/ht1_004_sensitivity.py`
- Generated output: `results/呼1-004_敏感性分析/`

- [ ] **Step 1: Compile the new script**

Run:

```powershell
python -m py_compile "scripts/ht1_004_sensitivity.py"
```

Expected: command exits with status 0 and no output.

- [ ] **Step 2: Run a one-combo smoke test**

Run:

```powershell
python "scripts/ht1_004_sensitivity.py" --limit 1
```

Expected:

- Prints `呼1-004井尾浆参数敏感性分析`.
- Prints `限制模式: 均匀采样 1 组覆盖整个参数空间`.
- Completes one simulation without traceback.
- Writes `results/呼1-004_敏感性分析/呼1-004_敏感性分析结果.csv`.
- Generates sensitivity/heatmap PNG files, though many heatmaps may have sparse one-row data due to `--limit 1`.

- [ ] **Step 3: Inspect generated CSV header**

Open or inspect the first row of:

```text
results/呼1-004_敏感性分析/呼1-004_敏感性分析结果.csv
```

Expected header exactly:

```text
排量_m3min,屈服应力_Pa,黏度_mPas,全井有效顶替效率,耗时_s
```

Expected data row count for `--limit 1`: one result row after the header.

- [ ] **Step 4: Report verification result**

Summarize whether both commands passed. If the smoke test is slow but still running, report that honestly and do not claim completion until it finishes.

## Self-Review

- Spec coverage: The plan creates the requested HT1-004 script, preserves the 003 parameter ranges/steps, uses the HT1-004 loader, keeps output/CLI behavior, and implements HT1-004 stop-time semantics with `nz=140`.
- Placeholder scan: No placeholders remain; every code step provides exact replacement content.
- Type consistency: `_annulus_stop_time_s(casing_result, fluids: tuple[FluidSpec, ...])` uses `FluidRole` imported from `cemdisp.data.fluid_spec`; `_run_single()` already receives `fluids`, so no new call-site data is needed.
