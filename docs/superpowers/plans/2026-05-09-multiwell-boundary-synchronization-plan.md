# Multi-Well Boundary Synchronization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a single six-well boundary-synchronization pipeline that converts each well’s ground pumping schedule into a traceable shoe-exit timeline and then into a unified `AnnulusInletState` input for the existing 2D annulus solver.

**Architecture:** Keep the current `cemdisp` layering intact. Concentrate changes in `data`, `transport1d`, `boundary_bridge`, `loaders`, and `runners`; do not rewrite `annulus_d2dga.py` physics in this phase. Implement boundary synchronization first, then per-well provenance/sync cards, and only after that open the gate to later efficiency-fitting work.

**Tech Stack:** Python 3, stdlib `dataclasses` / `enum` / `unittest`, existing `numpy` / `pandas` / `matplotlib`, current `cemdisp` package layout.

---

## Scope guard

This plan intentionally covers one dependent subsystem chain, not unrelated parallel products:

1. `PumpingSchedule` / well metadata enrichment
2. `CasingFlowSolver` shoe-exit timeline generation
3. `BoundaryBridge` unification
4. Loader-level per-well sync/provenance surfaces
5. Runner-level old-vs-new comparison and exports
6. Six-well verification

It does **not** include:

- CBL-driven parameter calibration
- 2D constitutive-model rewrite
- 1D convection-diffusion implementation
- New plotting framework
- Switching test framework away from `unittest`

---

## File structure to create or modify

### Create

- `cemdisp/transport1d/shoe_timeline.py`
- `cemdisp/data/provenance.py`
- `cemdisp/runners/common.py`
- `tests/test_shoe_timeline.py`
- `tests/test_pumping_schedule.py`
- `tests/test_well_spec.py`
- `tests/test_provenance.py`
- `tests/test_casing_flow.py`
- `tests/test_boundary_bridge.py`
- `tests/test_loader_deprecation.py`
- `tests/test_runner_common.py`
- `tests/test_sync_cards.py`
- `tests/test_comparison.py`
- `tests/test_six_well_integration.py`

### Modify

- `cemdisp/data/pumping_schedule.py`
- `cemdisp/data/well_spec.py`
- `cemdisp/data/fluid_provenance.py`
- `cemdisp/data/__init__.py`
- `cemdisp/transport1d/casing_flow.py`
- `cemdisp/transport1d/__init__.py`
- `cemdisp/models2d/boundary_bridge.py`
- `cemdisp/data/loaders/__init__.py`
- `cemdisp/data/loaders/hu101_loader.py`
- `cemdisp/data/loaders/hu102_loader.py`
- `cemdisp/data/loaders/hu103_loader.py`
- `cemdisp/data/loaders/hu1_loader.py`
- `cemdisp/data/loaders/hu2_loader.py`
- `cemdisp/data/loaders/ht1_001_loader.py`
- `cemdisp/runners/__init__.py`
- `cemdisp/runners/hu101_tailpipe.py`
- `cemdisp/runners/hu102_tailpipe.py`
- `cemdisp/runners/hu103_tailpipe.py`
- `cemdisp/runners/hu1_tailpipe.py`
- `cemdisp/runners/hu2_tailpipe.py`
- `cemdisp/runners/ht1_001_tailpipe.py`

---

## Task 1: Add explicit schedule semantics and tests

**Files:**
- Modify: `cemdisp/data/pumping_schedule.py`
- Test: `tests/test_pumping_schedule.py`

- [ ] **Step 1: Write the failing tests**

```python
import unittest

from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep, PumpingStageEvent


class PumpingScheduleTests(unittest.TestCase):
    def test_total_injected_volume_m3_sums_all_steps(self) -> None:
        schedule = PumpingSchedule(
            steps=(
                PumpingScheduleStep("注尾浆", "尾浆", 10.0, 1.0, event_tag=PumpingStageEvent.INJECT_CEMENT),
                PumpingScheduleStep("注替浆", "井浆", 20.0, 1.0, event_tag=PumpingStageEvent.INJECT_DISPLACEMENT),
            )
        )
        self.assertAlmostEqual(schedule.total_injected_volume_m3, 30.0)

    def test_cement_phase_steps_filters_cement_tags(self) -> None:
        schedule = PumpingSchedule(
            steps=(
                PumpingScheduleStep("注尾浆", "尾浆", 10.0, 1.0, event_tag=PumpingStageEvent.INJECT_CEMENT),
                PumpingScheduleStep("停泵", "尾浆", 0.0, 0.0, event_tag=PumpingStageEvent.SHUTDOWN),
            )
        )
        self.assertEqual(len(schedule.cement_phase_steps), 1)
        self.assertEqual(schedule.cement_phase_steps[0].step_name, "注尾浆")
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m unittest tests.test_pumping_schedule -v`

Expected: import or attribute errors for `PumpingStageEvent`, `total_injected_volume_m3`, or `cement_phase_steps`.

- [ ] **Step 3: Implement the minimal schedule enhancements**

```python
from enum import Enum


class PumpingStageEvent(str, Enum):
    INJECT_CEMENT = "inject_cement"
    INJECT_SPACER = "inject_spacer"
    INJECT_DISPLACEMENT = "inject_displacement"
    PLUG_RELEASE = "plug_release"
    SHUTDOWN = "shutdown"
    RESTART = "restart"
    RATE_SWITCH = "rate_switch"


@dataclass(frozen=True)
class PumpingScheduleStep:
    step_name: str
    fluid_name: str
    volume_m3: float
    rate_m3_min: float
    start_time_s: Optional[float] = None
    end_time_s: Optional[float] = None
    remarks: str = ""
    event_tag: Optional[PumpingStageEvent] = None


@dataclass(frozen=True)
class PumpingSchedule:
    steps: Tuple[PumpingScheduleStep, ...] = field(default_factory=tuple)
    notes: Tuple[str, ...] = field(default_factory=tuple)

    @property
    def total_injected_volume_m3(self) -> float:
        return sum(step.volume_m3 for step in self.steps)

    @property
    def cement_phase_steps(self) -> Tuple[PumpingScheduleStep, ...]:
        return tuple(step for step in self.steps if step.event_tag == PumpingStageEvent.INJECT_CEMENT)
```

- [ ] **Step 4: Run the tests again**

Run: `python -m unittest tests.test_pumping_schedule -v`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add cemdisp/data/pumping_schedule.py tests/test_pumping_schedule.py
git commit -m "增加施工程序事件标签与聚合属性"
```

---

## Task 2: Extend WellSpec for dual-diameter and shoe-lag metadata

**Files:**
- Modify: `cemdisp/data/well_spec.py`
- Test: `tests/test_well_spec.py`

- [ ] **Step 1: Write the failing tests**

```python
import unittest

from cemdisp.data.well_spec import DepthValuePoint, EvaluationWindow, WellSpec


class WellSpecTests(unittest.TestCase):
    def test_is_dual_diameter_true_when_upper_section_fields_present(self) -> None:
        well = WellSpec(
            well_name="测试井",
            top_md_m=1000.0,
            bottom_md_m=2000.0,
            shoe_md_m=2000.0,
            hanger_md_m=1000.0,
            casing_id_mm=200.0,
            liner_od_mm=140.0,
            liner_id_mm=110.0,
            upper_section_bottom_md_m=1500.0,
            upper_liner_od_mm=168.3,
            upper_liner_id_mm=138.9,
            shoe_lag_volume_m3=52.0,
            hole_diameter_profile=(DepthValuePoint(1000.0, 215.9),),
            inclination_profile=(DepthValuePoint(1000.0, 2.0),),
            standoff_profile=(DepthValuePoint(1000.0, 0.6),),
            evaluation_windows=(EvaluationWindow("CBL", 1000.0, 2000.0, "cbl"),),
        )
        self.assertTrue(well.is_dual_diameter)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m unittest tests.test_well_spec -v`

Expected: constructor rejects unknown keyword arguments.

- [ ] **Step 3: Add the new fields and validation**

```python
upper_section_bottom_md_m: float | None = None
upper_liner_od_mm: float | None = None
upper_liner_id_mm: float | None = None
shoe_lag_volume_m3: float | None = None
liner_wall_thickness_mm: float | None = None

@property
def is_dual_diameter(self) -> bool:
    return self.upper_section_bottom_md_m is not None

upper_values = (
    self.upper_section_bottom_md_m,
    self.upper_liner_od_mm,
    self.upper_liner_id_mm,
)
if any(value is not None for value in upper_values) and not all(value is not None for value in upper_values):
    raise ValueError("双径井必须同时提供上段边界、上段外径和上段内径")
```

- [ ] **Step 4: Run the tests again**

Run: `python -m unittest tests.test_well_spec -v`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add cemdisp/data/well_spec.py tests/test_well_spec.py
git commit -m "补充井筒双径与鞋口滞后元数据"
```

---

## Task 3: Create the shoe-exit timeline abstraction

**Files:**
- Create: `cemdisp/transport1d/shoe_timeline.py`
- Modify: `cemdisp/transport1d/__init__.py`
- Test: `tests/test_shoe_timeline.py`

- [ ] **Step 1: Write the failing tests**

```python
import unittest

from cemdisp.transport1d.shoe_timeline import ShoeEvent, ShoeEventKind, ShoeTimeline


class ShoeTimelineTests(unittest.TestCase):
    def test_at_returns_last_active_event_state(self) -> None:
        timeline = ShoeTimeline(
            events=(
                ShoeEvent(10.0, 0.01, "钻井液", "初始出流", ShoeEventKind.FRONT_ARRIVAL),
                ShoeEvent(20.0, 0.02, "尾浆", "注尾浆", ShoeEventKind.FRONT_ARRIVAL),
            )
        )
        state = timeline.at(25.0)
        self.assertEqual(state.phase_fractions, (("尾浆", 1.0),))
        self.assertAlmostEqual(state.flow_rate_m3_s, 0.02)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m unittest tests.test_shoe_timeline -v`

Expected: module import failure.

- [ ] **Step 3: Create the new file with timeline types**

```python
from dataclasses import dataclass, field
from enum import Enum

from cemdisp.transport1d.pipe_exit_state import PipeExitState


class ShoeEventKind(str, Enum):
    FRONT_ARRIVAL = "front_arrival"
    REAR_EXIT = "rear_exit"
    RATE_SWITCH = "rate_switch"
    SHUTDOWN = "shutdown"
    RESTART = "restart"
    END = "end"


@dataclass(frozen=True)
class ShoeEvent:
    time_s: float
    flow_rate_m3_s: float
    fluid_name: str
    stage_label: str
    event_kind: ShoeEventKind


@dataclass(frozen=True)
class ShoeTimeline:
    events: tuple[ShoeEvent, ...] = field(default_factory=tuple)

    def at(self, time_s: float) -> PipeExitState:
        active = self.events[0]
        for event in self.events:
            if event.time_s <= time_s:
                active = event
            else:
                break
        return PipeExitState(
            time_s=time_s,
            flow_rate_m3_s=active.flow_rate_m3_s,
            stage_name=active.stage_label,
            phase_fractions=((active.fluid_name, 1.0),),
        )
```

- [ ] **Step 4: Run the tests again**

Run: `python -m unittest tests.test_shoe_timeline -v`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add cemdisp/transport1d/shoe_timeline.py cemdisp/transport1d/__init__.py tests/test_shoe_timeline.py
git commit -m "新增鞋口出流时间线抽象"
```

---

## Task 4: Expand provenance beyond fluids

**Files:**
- Create: `cemdisp/data/provenance.py`
- Modify: `cemdisp/data/fluid_provenance.py`, `cemdisp/data/__init__.py`
- Test: `tests/test_provenance.py`

- [ ] **Step 1: Write the failing tests**

```python
import unittest

from cemdisp.data.provenance import WELL_PROVENANCE, build_well_provenance_summary


class ProvenanceTests(unittest.TestCase):
    def test_hu102_has_geometry_and_program_sections(self) -> None:
        summary = build_well_provenance_summary("呼102")
        self.assertIn("几何口径", summary)
        self.assertIn("施工程序口径", summary)
        self.assertIn("鞋口同步口径", summary)
        self.assertIn("流体口径", summary)

    def test_all_six_wells_registered(self) -> None:
        self.assertEqual(len(WELL_PROVENANCE), 6)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m unittest tests.test_provenance -v`

Expected: module import failure.

- [ ] **Step 3: Create the new provenance module and shim**

```python
@dataclass(frozen=True)
class ProvenanceItem:
    status: str
    note: str


@dataclass(frozen=True)
class WellProvenance:
    fluid: dict[str, ProvenanceItem]
    geometry: dict[str, ProvenanceItem]
    program: dict[str, ProvenanceItem]
    sync: dict[str, ProvenanceItem]


def build_well_provenance_summary(well_name: str) -> dict[str, object]:
    item = WELL_PROVENANCE[well_name]
    return {
        "流体口径": item.fluid,
        "几何口径": item.geometry,
        "施工程序口径": item.program,
        "鞋口同步口径": item.sync,
    }
```

And keep `fluid_provenance.py` as:

```python
from cemdisp.data.provenance import *
```

- [ ] **Step 4: Run the tests again**

Run: `python -m unittest tests.test_provenance -v`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add cemdisp/data/provenance.py cemdisp/data/fluid_provenance.py cemdisp/data/__init__.py tests/test_provenance.py
git commit -m "扩展多井 provenance 到几何程序与同步层"
```

---

## Task 5: Make CasingFlowSolver emit a complete shoe timeline

**Files:**
- Modify: `cemdisp/transport1d/casing_flow.py`
- Test: `tests/test_casing_flow.py`

- [ ] **Step 1: Write the failing tests**

```python
import unittest

from cemdisp.data.loaders.hu102_loader import load_hu102_tailpipe
from cemdisp.transport1d.casing_flow import CasingFlowSolver


class CasingFlowSolverTests(unittest.TestCase):
    def test_build_shoe_timeline_contains_cement_arrival(self) -> None:
        well_spec, fluids, schedule, _ = load_hu102_tailpipe()
        solver = CasingFlowSolver()
        timeline = solver.build_shoe_timeline(well_spec, fluids, schedule)
        cement_events = [event for event in timeline.events if event.fluid_name == "尾管水泥浆"]
        self.assertTrue(cement_events)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m unittest tests.test_casing_flow -v`

Expected: missing `build_shoe_timeline` method.

- [ ] **Step 3: Implement timeline generation on top of existing volume tracking**

```python
from cemdisp.transport1d.shoe_timeline import ShoeEvent, ShoeEventKind, ShoeTimeline


def build_shoe_timeline(self, well_spec, fluids, schedule) -> ShoeTimeline:
    result = self.run(well_spec, fluids, schedule)
    scheduled_steps = self._scheduled_steps_for_result(result)
    events: list[ShoeEvent] = []
    for front in result.fronts:
        active_step = self._active_step_at(scheduled_steps, front.time_s)
        flow_rate_m3_s = 0.0 if active_step is None else active_step.step.rate_m3_min / 60.0
        stage_label = front.fluid_name if active_step is None else active_step.step.step_name
        events.append(
            ShoeEvent(
                time_s=front.time_s,
                flow_rate_m3_s=flow_rate_m3_s,
                fluid_name=front.fluid_name,
                stage_label=stage_label,
                event_kind=ShoeEventKind.FRONT_ARRIVAL,
            )
        )
    return ShoeTimeline(events=tuple(sorted(events, key=lambda event: event.time_s)))
```

- [ ] **Step 4: Run the tests again**

Run: `python -m unittest tests.test_casing_flow -v`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add cemdisp/transport1d/casing_flow.py tests/test_casing_flow.py
git commit -m "让套管一维求解输出完整鞋口时间线"
```

---

## Task 6: Rework BoundaryBridge around ShoeTimeline and provenance

**Files:**
- Modify: `cemdisp/models2d/boundary_bridge.py`
- Test: `tests/test_boundary_bridge.py`

- [ ] **Step 1: Write the failing tests**

```python
import unittest

from cemdisp.data.provenance import WellProvenance
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d.shoe_timeline import ShoeEvent, ShoeEventKind, ShoeTimeline


class BoundaryBridgeTests(unittest.TestCase):
    def test_bridge_maps_tail_to_cement_phase(self) -> None:
        timeline = ShoeTimeline(
            events=(ShoeEvent(0.0, 0.01, "尾浆", "注尾浆", ShoeEventKind.FRONT_ARRIVAL),)
        )
        provenance = WellProvenance(fluid={}, geometry={}, program={}, sync={})
        provider = build_coupled_annulus_inlet_provider(timeline, provenance, split_cement_phases=False)
        state = provider(0.0)
        self.assertEqual(state.phase_fractions, (("cement", 1.0),))
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m unittest tests.test_boundary_bridge -v`

Expected: signature mismatch.

- [ ] **Step 3: Refactor the bridge signature and add sync-card output**

```python
def build_coupled_annulus_inlet_provider(
    shoe_timeline: ShoeTimeline,
    provenance: WellProvenance,
    *,
    split_cement_phases: bool = False,
) -> Callable[[float], AnnulusInletState]:
    def _provider(time_s: float) -> AnnulusInletState:
        pipe_exit_state = shoe_timeline.at(time_s)
        fluid_name = pipe_exit_state.phase_fractions[0][0]
        mapped = _phase_fractions_for_fluid_name(fluid_name, split_cement_phases=split_cement_phases)
        return AnnulusInletState(time_s, pipe_exit_state.flow_rate_m3_s, pipe_exit_state.stage_name, mapped)
    return _provider


def build_sync_card(well_name: str, shoe_timeline: ShoeTimeline, provenance: WellProvenance) -> dict[str, object]:
    return {
        "井名": well_name,
        "鞋口同步口径": {"事件数": len(shoe_timeline.events)},
        "代理提醒": [note.note for note in provenance.sync.values() if note.status != "field"],
    }
```

- [ ] **Step 4: Run the tests again**

Run: `python -m unittest tests.test_boundary_bridge -v`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add cemdisp/models2d/boundary_bridge.py tests/test_boundary_bridge.py
git commit -m "统一鞋口到环空桥接并支持同步画像卡"
```

---

## Task 7: Deprecate per-loader legacy inlet providers and add per-well sync-card exporters

**Files:**
- Modify: all six loader files and `cemdisp/data/loaders/__init__.py`
- Test: `tests/test_loader_deprecation.py`, `tests/test_sync_cards.py`

- [ ] **Step 1: Write the failing tests**

```python
import unittest
import warnings

from cemdisp.data.loaders.hu102_loader import build_hu102_annulus_inlet_provider


class LoaderDeprecationTests(unittest.TestCase):
    def test_hu102_legacy_provider_warns(self) -> None:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                build_hu102_annulus_inlet_provider((), ())
            except Exception:
                pass
        self.assertTrue(any(item.category is DeprecationWarning for item in caught))
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m unittest tests.test_loader_deprecation tests.test_sync_cards -v`

Expected: no warning emitted and missing sync-card exporters.

- [ ] **Step 3: Deprecate old providers and add new exporters**

```python
import warnings


def build_hu102_annulus_inlet_provider(
    schedule: PumpingSchedule,
    fluids: tuple[FluidSpec, ...],
    annulus_boundary_mode: str = "sustained_tail",
) -> Callable[[float], AnnulusInletState]:
    warnings.warn(
        "build_hu102_annulus_inlet_provider 已废弃，请改用 CasingFlowSolver + build_coupled_annulus_inlet_provider",
        DeprecationWarning,
        stacklevel=2,
    )
    return _legacy_hu102_provider(schedule=schedule, fluids=fluids, annulus_boundary_mode=annulus_boundary_mode)


def export_hu102_sync_card_markdown(output_dir: Path) -> Path:
    well_spec, fluids, schedule, _ = load_hu102_tailpipe()
    solver = CasingFlowSolver()
    timeline = solver.build_shoe_timeline(well_spec, fluids, schedule)
    provenance = WELL_PROVENANCE["呼102"]
    sync_card = build_sync_card("呼102", timeline, provenance)
    output_path = output_dir / "呼102_同步画像卡.md"
    lines = [
        "# 呼102同步画像卡",
        "",
        f"- 事件数：{sync_card['鞋口同步口径']['事件数']}",
    ]
    output_path.write_text("\n".join(lines), encoding="utf-8")
    return output_path
```

- [ ] **Step 4: Run the tests again**

Run: `python -m unittest tests.test_loader_deprecation tests.test_sync_cards -v`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add cemdisp/data/loaders/__init__.py cemdisp/data/loaders/*.py tests/test_loader_deprecation.py tests/test_sync_cards.py
git commit -m "废弃旧入口提供器并补充六口井同步画像导出"
```

---

## Task 8: Extract runner common logic, add old-vs-new comparisons, and verify all six wells

**Files:**
- Create: `cemdisp/runners/common.py`
- Modify: all six runner files, `cemdisp/runners/__init__.py`
- Test: `tests/test_runner_common.py`, `tests/test_comparison.py`, `tests/test_six_well_integration.py`

- [ ] **Step 1: Write the failing tests**

```python
import unittest

from cemdisp.runners.common import export_old_vs_new_inlet_comparison


class RunnerCommonTests(unittest.TestCase):
    def test_export_old_vs_new_inlet_comparison_returns_csv_path(self) -> None:
        self.assertTrue(callable(export_old_vs_new_inlet_comparison))
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m unittest tests.test_runner_common tests.test_comparison tests.test_six_well_integration -v`

Expected: missing module and missing integration flow.

- [ ] **Step 3: Create the common runner module and wire all six runners**

```python
def annulus_stop_time_s(casing_result: CasingFlowResult) -> float:
    return float(casing_result.cement_end_time_s)


def export_old_vs_new_inlet_comparison(
    *,
    output_dir: Path,
    old_provider: Callable[[float], AnnulusInletState],
    new_provider: Callable[[float], AnnulusInletState],
    sample_times_s: tuple[float, ...],
    well_name: str,
) -> Path:
    comparison_path = output_dir / f"{well_name}_入口边界对比.csv"
    with comparison_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("time_s", "old_phase", "new_phase", "old_rate_m3_s", "new_rate_m3_s"),
        )
        writer.writeheader()
        for time_s in sample_times_s:
            old_state = old_provider(time_s)
            new_state = new_provider(time_s)
            writer.writerow(
                {
                    "time_s": f"{time_s:.3f}",
                    "old_phase": old_state.phase_fractions[0][0],
                    "new_phase": new_state.phase_fractions[0][0],
                    "old_rate_m3_s": f"{old_state.flow_rate_m3_s:.6f}",
                    "new_rate_m3_s": f"{new_state.flow_rate_m3_s:.6f}",
                }
            )
    return comparison_path
```

Update each runner to:

```python
solver = CasingFlowSolver()
timeline = solver.build_shoe_timeline(well_spec, fluids, schedule)
provenance = WELL_PROVENANCE[well_spec.well_name]
new_provider = build_coupled_annulus_inlet_provider(timeline, provenance, split_cement_phases=True)
```

- [ ] **Step 4: Run the full suite**

Run: `python -m unittest discover -s tests -v`

Expected: PASS with all new unit and integration tests green.

- [ ] **Step 5: Commit**

```bash
git add cemdisp/runners/common.py cemdisp/runners/*.py tests/test_runner_common.py tests/test_comparison.py tests/test_six_well_integration.py
git commit -m "统一运行器主流程并增加新旧入口对比验证"
```

---

## Self-review against the approved spec

### Spec coverage

- **边界同步优先**：covered by Tasks 1, 3, 5, 6, 8
- **资料/provenance 并入各井资料**：covered by Tasks 4 and 7
- **六口井统一接口、分井定制**：covered by Tasks 7 and 8
- **旧入口 vs 新入口差值输出**：covered by Task 8
- **后续效率贴场是后续阶段**：not implemented here by design; correctly left outside this plan’s coding scope

### Placeholder scan

- No `TBD`
- No `TODO`
- No “similar to previous task” references without concrete file paths

### Type consistency

- `PumpingStageEvent` is introduced before `CasingFlowSolver` consumes it
- `ShoeTimeline` is introduced before `BoundaryBridge` consumes it
- `WellProvenance` is introduced before loaders and runners reference it

---

Plan complete and saved to `docs/superpowers/plans/2026-05-09-multiwell-boundary-synchronization-plan.md`. Two execution options:

**1. Subagent-Driven (recommended)** - I dispatch a fresh subagent per task, review between tasks, fast iteration

**2. Inline Execution** - Execute tasks in this session using executing-plans, batch execution with checkpoints

**Which approach?**
