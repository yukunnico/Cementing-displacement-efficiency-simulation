# Paper Loader Results Pipeline Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a paper-only loader/results pipeline that isolates paper data provenance, model runs, result tables, and validation outputs under `results_paper/pead_v1/` without changing legacy loader behavior.

**Architecture:** Add a narrow `cemdisp.data.loaders.paper` package that wraps existing well loaders with paper metadata, source classification, and quality checks while preserving the existing `(WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData)` contract. Add reusable paper runner utilities in `scripts/paper_data/` so main results, CBL validation, D2DGA ablation, and grid/time-step validation all share the same registry and CSV/Markdown schemas.

**Tech Stack:** Python 3.9+, dataclasses, pathlib, csv/json, numpy, pandas, matplotlib-backed existing reporting helpers, pytest/unittest.

---

## File Structure

- Create `参考文档/现场资料提取/01_总表/loader优化采用值对照表.csv`: paper adopted-values control table with one row per adopted field/source decision.
- Create `cemdisp/data/loaders/paper/__init__.py`: public exports for paper loaders and registry helpers.
- Create `cemdisp/data/loaders/paper/common.py`: metadata dataclasses, JSON/CSV helpers, loader validation, cbl/window helper functions.
- Create `cemdisp/data/loaders/paper/registry.py`: authoritative paper well registry for `hu101`, `hu102`, `hu103`, and `ht1_001`.
- Create `cemdisp/data/loaders/paper/paper_hu101_loader.py`: paper wrapper for `load_hu101_tailpipe()` with CBL pass rate 0.6277 and paper metadata.
- Create `cemdisp/data/loaders/paper/paper_hu102_loader.py`: paper wrapper for `load_hu102_tailpipe()` and paper metadata.
- Create `cemdisp/data/loaders/paper/paper_hu103_loader.py`: paper wrapper for `load_hu103_tailpipe()` with CBL pass rate 0.1206 and explicit window notes.
- Create `cemdisp/data/loaders/paper/paper_ht1_001_loader.py`: paper wrapper for `load_ht1_001_tailpipe()` marked pending validation / not included in CBL metrics.
- Create `scripts/paper_data/__init__.py`: make scripts importable in tests.
- Create `scripts/paper_data/run_paper_main_results.py`: shared simulation runner and main results CLI.
- Create `scripts/paper_data/build_cbl_validation_tables.py`: CBL error metrics and Markdown comparison builder.
- Create `scripts/paper_data/build_paper_tables.py`: aggregate paper tables from per-well summaries.
- Create `scripts/paper_data/run_paper_d2dga_ablation.py`: run `enable_d2dga=True/False` comparisons.
- Create `scripts/paper_data/run_paper_grid_dt_validation.py`: run grid/time-step sensitivity cases.
- Create `tests/paper_data/test_paper_loaders.py`: loader contract and metadata tests.
- Create `tests/paper_data/test_paper_result_schema.py`: summary/table schema tests without long simulations.
- Create `tests/paper_data/test_cbl_metrics.py`: CBL formula tests.

Known baseline note: before implementation, `python -m pytest tests/test_loader_deprecation.py tests/test_well_spec.py tests/test_pumping_schedule.py -q` has 3 pre-existing failures in `tests/test_well_spec.py::TestWellSpecPositionalBackwardCompatibility`; do not treat those as new paper-pipeline regressions unless this task edits `cemdisp/data/well_spec.py`.

---

### Task 1: Paper adopted-values manifest

**Files:**
- Create: `参考文档/现场资料提取/01_总表/loader优化采用值对照表.csv`

- [ ] **Step 1: Write the adopted-values CSV**

Use this exact header:

```csv
well_id,field_group,field_name,field_value,loader_value,adopted_value,unit,adopted_source_type,confidence,action,notes
```

Include at least these rows:

```csv
hu101,cbl,cbl_pass_rate,0.6277,,0.6277,fraction,interpreted,high,adopt,"100312.PDF CBL合格率62.77%; paper CBL评价窗5700-7810m"
hu101,geometry,upper_composite_liner,"168.3mm liner in 260.35mm hole","equivalent hole diameter with 139.7mm OD reference","area-preserving equivalent geometry",rule,model_assumption,medium,adopt,"paper loader keeps explicit equivalent-geometry note"
hu101,program,lead_volume_m3,47.0,47.0,47.0,m3,field_measured,high,adopt,"现场抽取施工量"
hu101,program,tail_volume_m3,23.0,23.0,23.0,m3,field_measured,high,adopt,"现场抽取施工量"
hu102,cbl,cbl_pass_rate,0.6665,0.6665,0.6665,fraction,interpreted,high,adopt,"100413.PDF CBL合格率66.65%; CBL评价窗6840-7665m"
hu102,geometry,caliper_profile,"field profile available","legacy loader profile","legacy loader current profile pending full CSV rebuild",profile,legacy_value,medium,adopt_for_v1,"paper metadata marks hu102 as field-rich but loader-rebuild sensitive"
hu102,program,pumping_schedule,"field-rich construction records","legacy reconstructed schedule","legacy reconstructed schedule",schedule,legacy_value,medium,adopt_for_v1,"future improvement should rebuild from extracted pumping CSV"
hu103,cbl,cbl_pass_rate,0.1206,,0.1206,fraction,interpreted,high,adopt,"139.7mm段7338-7712m; do not mix 168.3mm段0.0004或综合0.0605"
hu103,geometry,lower_liner_interval,"139.7mm tail liner 7330.6-7770m","139.7mm single OD",139.7,mm,field_measured,high,adopt,"paper main validation window is 7338-7712m"
hu103,program,tail_volume_m3,19.0,19.0,19.0,m3,field_measured,medium,adopt,"含下塞1.18m3; intermediate slurry remains proxy"
ht1_001,sample_class,main_validation_status,"pending cleaned schedule and frozen CBL window","legacy/proxy loader","application_pending_validation",class,model_assumption,low,defer,"not included in CBL error metrics until cleaned schedule and CBL window are frozen"
```

- [ ] **Step 2: Verify CSV parses**

Run:

```bash
python - <<'PY'
from pathlib import Path
import csv
p = Path('参考文档/现场资料提取/01_总表/loader优化采用值对照表.csv')
rows = list(csv.DictReader(p.open(encoding='utf-8-sig')))
assert rows, 'adopted values table is empty'
assert {'well_id','field_group','field_name','adopted_source_type','action'}.issubset(rows[0])
print(len(rows))
PY
```

Expected: prints a row count greater than 0.

---

### Task 2: Paper loader package and registry

**Files:**
- Create: `cemdisp/data/loaders/paper/common.py`
- Create: `cemdisp/data/loaders/paper/registry.py`
- Create: `cemdisp/data/loaders/paper/__init__.py`
- Create: `cemdisp/data/loaders/paper/paper_hu101_loader.py`
- Create: `cemdisp/data/loaders/paper/paper_hu102_loader.py`
- Create: `cemdisp/data/loaders/paper/paper_hu103_loader.py`
- Create: `cemdisp/data/loaders/paper/paper_ht1_001_loader.py`
- Test: `tests/paper_data/test_paper_loaders.py`

- [ ] **Step 1: Write failing tests**

Add tests that import `cemdisp.data.loaders.paper`, call all four paper loaders, assert return types, assert `metadata.sample_class`, assert CBL windows are within `WellSpec.top_md_m..bottom_md_m`, assert total pumping volume is positive, assert hu101/hu102/hu103 have measured CBL pass rates, and assert `iter_paper_wells(include_pending=True)` returns four records.

- [ ] **Step 2: Run tests to verify failure**

Run:

```bash
python -m pytest tests/paper_data/test_paper_loaders.py -q
```

Expected before implementation: import failure for `cemdisp.data.loaders.paper`.

- [ ] **Step 3: Implement `common.py`**

Define:

```python
@dataclass(frozen=True)
class PaperWellMetadata:
    well_id: str
    well_name_cn: str
    paper_version: str
    sample_class: str
    loader_version: str
    include_in_main_results: bool
    include_in_cbl_metrics: bool
    cbl_pass_rate: float | None
    geometry_source: str
    pumping_schedule_source: str
    fluid_source: str
    caliper_source: str
    inclination_source: str
    standoff_source: str
    cbl_source: str
    notes: tuple[str, ...] = ()
```

Also implement `metadata_to_dict()`, `validation_with_cbl_pass_rate()`, `get_window()`, `validate_paper_loader_result()`, and `write_markdown_table()`.

- [ ] **Step 4: Implement paper wrapper loaders**

Each wrapper exports `load_paper_<well>_tailpipe()` and `get_paper_metadata()`. The loader calls the legacy loader, replaces `ValidationData.cbl_pass_rate` where the paper metadata has one, validates the result, and returns the same four-tuple.

- [ ] **Step 5: Implement registry**

Define `PaperWellRecord` with `well_id`, `loader`, and `metadata_factory`. Export `PAPER_WELL_REGISTRY`, `get_paper_well_record(well_id)`, `iter_paper_wells(include_pending=False)`, and `load_paper_well(well_id)`.

- [ ] **Step 6: Run tests**

Run:

```bash
python -m pytest tests/paper_data/test_paper_loaders.py -q
```

Expected: all tests in the file pass.

---

### Task 3: Paper main result runner and result schema

**Files:**
- Create: `scripts/paper_data/__init__.py`
- Create: `scripts/paper_data/run_paper_main_results.py`
- Test: `tests/paper_data/test_paper_result_schema.py`

- [ ] **Step 1: Write schema tests**

Test a pure helper `summary_row_from_files(summary_path, metadata, output_root)` with a minimal temporary `summary.json`. Assert it emits the exact spec fields:

```python
[
    'well_id','well_name_cn','sample_class','top_md_m','bottom_md_m','annulus_volume_m3',
    'final_effective_efficiency','final_bulk_cement_fill','cbl_eval_interval_efficiency',
    'target_interval_efficiency','channeling_index','mixing_index','instability_index',
    'summary_json_path','time_series_csv_path','depth_profile_csv_path'
]
```

- [ ] **Step 2: Implement runner helpers**

Implement `PaperModelConfig`, `default_paper_model_config()`, `smoke_paper_model_config()`, `ensure_results_tree()`, `run_single_paper_well()`, `write_single_well_outputs()`, `summary_row_from_files()`, `write_all_wells_main_results()`, and CLI `main()`.

- [ ] **Step 3: Use shared simulation flow**

`run_single_paper_well()` must follow the existing flow in `scripts/p3_p4_integration.py`: paper loader → `CasingFlowSolver` → `build_coupled_annulus_inlet_provider(split_cement_phases=True)` → `AnnulusD2DGASolver(dt=config.dt, nz=config.nz, ny=config.ny, total_t=casing_result.cement_end_time_s, enable_d2dga=config.enable_d2dga, yield_regularization_M=config.yield_regularization_M, open_outlet=config.open_outlet)`.

- [ ] **Step 4: Write required per-well outputs**

For each well directory under `results_paper/pead_v1/01_main_results/<well_id>/`, write `run_config.json`, `adopted_input_summary.json`, `time_series.csv`, `depth_profile.csv`, `field_2d.npz`, `summary.json`, `summary.md`, and `provenance.md`. Figure creation is optional but, if plotting succeeds, write `final_field.png` and `depth_efficiency.png`; plotting failures must be recorded in `summary.json['warnings']` instead of failing the run.

- [ ] **Step 5: Run tests**

Run:

```bash
python -m pytest tests/paper_data/test_paper_result_schema.py -q
```

Expected: pass.

---

### Task 4: CBL validation and paper table builders

**Files:**
- Create: `scripts/paper_data/build_cbl_validation_tables.py`
- Create: `scripts/paper_data/build_paper_tables.py`
- Test: `tests/paper_data/test_cbl_metrics.py`

- [ ] **Step 1: Write CBL formula tests**

Test `compute_cbl_error_metrics('hu101', predicted=0.60, measured=0.50, included=True, notes='x')` returns `absolute_error == 0.10` and `relative_error_percent == 20.0`.

- [ ] **Step 2: Implement CBL builder**

Read per-well `summary.json` files, combine with registry metadata, and write:

```text
results_paper/pead_v1/02_cbl_validation/cbl_truth_table.csv
results_paper/pead_v1/02_cbl_validation/cbl_error_metrics.csv
results_paper/pead_v1/02_cbl_validation/cbl_comparison.md
```

Rows for ht1_001 must set `included_in_metrics=False` and blank measured/error values.

- [ ] **Step 3: Implement main table builder**

`build_paper_tables.py` should rebuild `01_main_results/all_wells_main_results.csv` and `.md` from per-well summaries using the same field order as Task 3.

- [ ] **Step 4: Run tests**

Run:

```bash
python -m pytest tests/paper_data/test_cbl_metrics.py -q
```

Expected: pass.

---

### Task 5: D2DGA ablation and grid/time-step validation runners

**Files:**
- Create: `scripts/paper_data/run_paper_d2dga_ablation.py`
- Create: `scripts/paper_data/run_paper_grid_dt_validation.py`

- [ ] **Step 1: Implement D2DGA ablation CLI**

Use `run_single_paper_well()` twice per requested well with `enable_d2dga=True` and `False`. Write `results_paper/pead_v1/03_d2dga_ablation/d2dga_on_off_comparison.csv` with fields:

```csv
well_id,enable_d2dga,final_effective_efficiency,cbl_eval_interval_efficiency,channeling_index,mixing_index,instability_index,runtime_s
```

Also write `d2dga_on_off_comparison.md`.

- [ ] **Step 2: Implement grid/time-step CLI**

Use representative wells defaulting to `hu101` and `hu103`. Grid schemes are `250/20`, `500/40`, `750/60`; timestep schemes are `2`, `4`, `8`. Provide `--quick` to run a tiny validation set (`nz/ny=60/12` and `80/16`, `dt=8` and `dt=10`) for output checks. Write:

```text
results_paper/pead_v1/04_grid_dt_validation/grid_convergence.csv
results_paper/pead_v1/04_grid_dt_validation/timestep_sensitivity.csv
results_paper/pead_v1/04_grid_dt_validation/numerical_reliability.md
```

- [ ] **Step 3: Run import smoke tests**

Run:

```bash
python - <<'PY'
import scripts.paper_data.run_paper_d2dga_ablation as a
import scripts.paper_data.run_paper_grid_dt_validation as g
print(a.__name__)
print(g.__name__)
PY
```

Expected: both module names print.

---

### Task 6: Execute, validate outputs, and review

**Files:**
- Output: `results_paper/pead_v1/**`

- [ ] **Step 1: Run focused paper tests**

Run:

```bash
python -m pytest tests/paper_data -q
```

Expected: pass.

- [ ] **Step 2: Run smoke main results**

Run:

```bash
python scripts/paper_data/run_paper_main_results.py --mode smoke --wells hu101 hu102 hu103 ht1_001
```

Expected: creates per-well outputs and `all_wells_main_results.csv` under `results_paper/pead_v1/01_main_results/`.

- [ ] **Step 3: Build validation tables**

Run:

```bash
python scripts/paper_data/build_cbl_validation_tables.py
python scripts/paper_data/build_paper_tables.py
```

Expected: writes CBL and main Markdown/CSV tables.

- [ ] **Step 4: Run quick method-validation outputs**

Run:

```bash
python scripts/paper_data/run_paper_d2dga_ablation.py --mode smoke --wells hu101 hu103
python scripts/paper_data/run_paper_grid_dt_validation.py --quick --wells hu101 hu103
```

Expected: writes D2DGA and grid/time-step CSV/Markdown outputs.

- [ ] **Step 5: Run output schema/range check**

Run:

```bash
python - <<'PY'
from pathlib import Path
import csv, json
root = Path('results_paper/pead_v1')
required = [
    root/'00_run_config'/'paper_model_config.json',
    root/'00_run_config'/'adopted_loader_manifest.csv',
    root/'01_main_results'/'all_wells_main_results.csv',
    root/'02_cbl_validation'/'cbl_error_metrics.csv',
    root/'03_d2dga_ablation'/'d2dga_on_off_comparison.csv',
    root/'04_grid_dt_validation'/'grid_convergence.csv',
    root/'04_grid_dt_validation'/'timestep_sensitivity.csv',
]
missing = [str(p) for p in required if not p.exists()]
assert not missing, missing
rows = list(csv.DictReader((root/'01_main_results'/'all_wells_main_results.csv').open(encoding='utf-8-sig')))
assert rows, 'main results empty'
for row in rows:
    for key in ('final_effective_efficiency','final_bulk_cement_fill','cbl_eval_interval_efficiency'):
        value = row.get(key, '')
        if value:
            f = float(value)
            assert 0.0 <= f <= 1.0, (row['well_id'], key, f)
print('checked', len(rows), 'main-result rows')
PY
```

Expected: prints checked row count.

- [ ] **Step 6: If time allows, run formal main configuration**

Run:

```bash
python scripts/paper_data/run_paper_main_results.py --mode formal --wells hu101 hu102 hu103 ht1_001
```

Expected: updates main results with `dt=4`, `nz=500`, `ny=40`. If runtime is too long for the session, keep smoke outputs and report formal run as skipped/not completed.

---

## Self-Review

- Spec coverage: The plan creates the adopted-values table, paper loader layer, registry, paper results root, main runner, CBL validation builder, D2DGA ablation runner, grid/time validation runner, and paper table builders requested by the design doc.
- Placeholder scan: No task uses TBD/TODO/fill-later language; pending status in metadata is an explicit sample classification for ht1_001, not a plan placeholder.
- Type consistency: Loader return type remains the existing four-tuple; metadata is a separate dataclass converted to dict for JSON outputs; result table fields match the design document.
- Scope control: hu101/hu102/hu103 are main validation wells; ht1_001 is supported but marked pending validation; ht1_003/ht1_004 remain out of this first paper pipeline because the approved design says they are sensitivity/optimization cases unless field verification is added.
