# -*- coding: utf-8 -*-
"""Task 9（T1-4）契约测试：`AnnulusInletState.temperature_c` 温度接口位。

覆盖（task-9-brief 需求逐条）：

1. **默认 NaN=未知**：不传 `temperature_c` 构造，读出为 NaN（`math.isnan`）。
2. **旧构造签名向后兼容**：位置参数四参旧签名（`AnnulusInletState(t, q, stage, fracs)`）
   与既有 keyword 构造直接可跑；既有桥接函数（`pipe_exit_to_annulus_inlet` /
   `build_coupled_annulus_inlet_provider` 两分支）输出默认 NaN——它们不传温度，
   零改动即兼容。
3. **新字段可传可读**：显式传值原样读回；`dataclasses.replace` 可派生新实例。
4. **frozen 不可变**：对实例赋 `temperature_c` 抛 `FrozenInstanceError`。
"""
from __future__ import annotations

import math
from dataclasses import FrozenInstanceError, replace

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.provenance import FluidProvenance, SectionProvenance, WellProvenance
from cemdisp.data.well_spec import DepthValuePoint, EvaluationWindow, WellSpec
from cemdisp.models2d.boundary_bridge import (
    AnnulusInletState,
    build_coupled_annulus_inlet_provider,
    pipe_exit_to_annulus_inlet,
)
from cemdisp.transport1d.pipe_exit_state import PipeExitState
from cemdisp.transport1d.shoe_timeline import ShoeEvent, ShoeEventKind, ShoeTimeline

_Q_M3S = 1.0 / 60.0


# --------------------------------------------------------------------------- #
# 1. 默认 NaN
# --------------------------------------------------------------------------- #

def test_default_temperature_is_nan() -> None:
    """不传 temperature_c ⇒ 默认 NaN=未知。"""

    state = AnnulusInletState(
        time_s=0.0, flow_rate_m3_s=_Q_M3S, stage_name="注尾浆",
        phase_fractions=(("cement", 1.0),),
    )
    assert math.isnan(state.temperature_c)


def test_minimal_positional_defaults_to_nan() -> None:
    """最简 keyword 构造（连 phase_fractions 都走默认）也读出 NaN。"""

    state = AnnulusInletState(time_s=0.0, flow_rate_m3_s=_Q_M3S, stage_name="s")
    assert math.isnan(state.temperature_c)
    assert state.phase_fractions == ()


# --------------------------------------------------------------------------- #
# 2. 旧构造签名向后兼容（既有构造点零改动）
# --------------------------------------------------------------------------- #

def test_legacy_positional_signature_still_works() -> None:
    """四位置参数旧签名直接跑，新字段取默认 NaN。"""

    state = AnnulusInletState(10.0, _Q_M3S, "前置液", (("spacer", 1.0),))
    assert state.time_s == 10.0
    assert state.flow_rate_m3_s == _Q_M3S
    assert state.stage_name == "前置液"
    assert state.phase_fractions == (("spacer", 1.0),)
    assert math.isnan(state.temperature_c)


def test_pipe_exit_mapping_defaults_to_nan() -> None:
    """pipe_exit_to_annulus_inlet 不传温度 ⇒ 输出默认 NaN（零改动兼容）。"""

    exit_state = PipeExitState(
        time_s=5.0,
        flow_rate_m3_s=_Q_M3S,
        stage_name="注水泥",
        phase_fractions=(("水泥", 1.0),),
    )
    state = pipe_exit_to_annulus_inlet(exit_state)
    assert state.time_s == 5.0
    assert math.isnan(state.temperature_c)


def test_timeline_provider_defaults_to_nan() -> None:
    """新分支 build_coupled_annulus_inlet_provider 输出默认 NaN。"""

    timeline = ShoeTimeline(events=(
        ShoeEvent(
            time_s=0.0,
            kind=ShoeEventKind.FRONT_ARRIVAL,
            flow_rate_m3_s=_Q_M3S,
            stage_name="注尾浆",
            phase_fractions=(("尾浆", 1.0),),
        ),
    ))
    fluids = (
        FluidSpec(name="尾浆", role=FluidRole.TAIL, density_kg_m3=1900.0,
                  plastic_viscosity_pa_s=0.08),
    )
    provenance = WellProvenance(
        well_name="测试井",
        fluid={},
        geometry=SectionProvenance("field", ""),
        program=SectionProvenance("field", ""),
        sync=SectionProvenance("field", ""),
    )
    provider = build_coupled_annulus_inlet_provider(timeline, provenance, fluids)
    state = provider(0.0)
    assert state.phase_fractions == (("cement", 1.0),)
    assert math.isnan(state.temperature_c)


def test_annulus_solver_tolerates_default_nan_inlet() -> None:
    """最小 2D 求解消费默认 NaN 入口可跑完（既有消费方不破坏）。

    本轮温度字段无 2D 消费方（接口位先行）——NaN 不进任何数值路径，
    求解照常完成即为向后兼容证据。
    """

    from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver

    top, bottom = 100.0, 300.0
    well = WellSpec(
        well_name="T14_inlet",
        top_md_m=top,
        bottom_md_m=bottom,
        shoe_md_m=bottom,
        hole_diameter_profile=(DepthValuePoint(top, 260.0), DepthValuePoint(bottom, 250.0)),
        inclination_profile=(DepthValuePoint(top, 0.0), DepthValuePoint(bottom, 0.0)),
        standoff_profile=(DepthValuePoint(top, 0.8), DepthValuePoint(bottom, 0.8)),
        liner_od_mm=139.7,
        liner_id_mm=124.3,
        evaluation_windows=(EvaluationWindow("target", top, bottom),),
    )
    fluids = (
        FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1200.0,
                  plastic_viscosity_pa_s=0.02),
        FluidSpec(name="tail", role=FluidRole.TAIL, density_kg_m3=1900.0,
                  plastic_viscosity_pa_s=0.08),
    )

    def provider(t: float) -> AnnulusInletState:
        # 不传 temperature_c ⇒ 默认 NaN（正是要喂给 2D 的旧口径状态）
        return AnnulusInletState(t, _Q_M3S, "tail", (("tail", 1.0),))

    solver = AnnulusD2DGASolver(nz=8, ny=6, dt=4.0, total_t=8.0,
                                enable_cfl_adaptive=False, open_outlet=True)
    result = solver.run(well, fluids, provider)
    assert result is not None


# --------------------------------------------------------------------------- #
# 3. 新字段可传可读
# --------------------------------------------------------------------------- #

def test_explicit_temperature_roundtrip() -> None:
    """显式传值原样读回。"""

    state = AnnulusInletState(
        time_s=0.0, flow_rate_m3_s=_Q_M3S, stage_name="注尾浆",
        phase_fractions=(("cement", 1.0),),
        temperature_c=45.5,
    )
    assert state.temperature_c == 45.5
    assert not math.isnan(state.temperature_c)


def test_explicit_temperature_zero_is_not_unknown() -> None:
    """显式 0.0 是合法温度，不与默认 NaN 混淆。"""

    state = AnnulusInletState(0.0, _Q_M3S, "s", (), temperature_c=0.0)
    assert state.temperature_c == 0.0


def test_replace_can_set_temperature() -> None:
    """dataclasses.replace 派生新实例可写入温度（frozen 下的合法改法）。"""

    base = AnnulusInletState(0.0, _Q_M3S, "s", (("mud", 1.0),))
    warm = replace(base, temperature_c=60.0)
    assert math.isnan(base.temperature_c)          # 原实例不动
    assert warm.temperature_c == 60.0
    assert warm.time_s == base.time_s


# --------------------------------------------------------------------------- #
# 4. frozen 不可变
# --------------------------------------------------------------------------- #

def test_temperature_field_is_immutable() -> None:
    """frozen dataclass：实例上直接赋 temperature_c 抛 FrozenInstanceError。"""

    state = AnnulusInletState(0.0, _Q_M3S, "s", ())
    try:
        state.temperature_c = 30.0  # type: ignore[misc]
    except FrozenInstanceError:
        pass
    else:
        raise AssertionError("frozen dataclass 应拒绝 temperature_c 赋值")
    assert math.isnan(state.temperature_c)
