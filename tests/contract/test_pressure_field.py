# -*- coding: utf-8 -*-
"""P-1（2026-10-06）契约测试：静液柱压力场 + 求解器接线。

上游
----
- ``docs/superpowers/specs/2026-10-06-phase1-fluid-at-signature-wave-design.md`` §2
- 评审A §4 断言 **A3.1 / A3.2 / A3.3** + 缺陷 D-06（构造层口径分裂）
- 执行窗口对抗核查 MAJOR：**C-09**（无场须返回 None 而非 0.0）、**C-16**（casing 无 geom）、
  **C-13**（P-1 经混合 τy 场灌入屈服门的耦合面）、**C-12**（合成隔离液双重外推）

覆盖
----
1. 均匀柱：``P(md) == g·ρ̄·TVD(md)/1e6`` 解析对照（浮点容差 <1e-12）
2. 分层柱：介于两端均匀柱之间 + 与手算分段累加一致
3. ``ConstantPressureField`` + oob 三件套
4. 三井真实数据对账：呼1-003 设计锚 145.96 MPa **±2%**、呼101 现场锚 152.16 MPa ±2%
5. 关2：无场 ⇒ ``fluid_at`` 收 P=None ⇒ ``p_default`` 审计保留（且不得走 clamp）
6. A3.1：构造层物性 ≡ 逐步派生（同 memo 键 ⇒ 同一 P 档）
7. A3.2：传场 run 内 ``p_default`` 计数 **== 0**
8. A3.3：memo 键含 P ⇒ 同 (fluid,T) 不同 P 返回**不同对象**且各自手算对（键漏 P 必红）
9. C-13：spacer τy 随 P 增 ⇒ 混合屈服应力场/屈服门输入随之改变（方向锚）
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[2]
for _p in (str(_ROOT), str(_ROOT / "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel  # noqa: E402
from cemdisp.data.pressure_field import (  # noqa: E402
    ANNULUS_ROLES,
    ConstantPressureField,
    HydrostaticPressureField,
    insitu_column_density,
)
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep  # noqa: E402
from cemdisp.data.rheology_vs_temperature import (  # noqa: E402
    fluid_at,
    get_audit,
    reset_audit,
)
from cemdisp.data.well_spec import DepthValuePoint, WellSpec  # noqa: E402
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver  # noqa: E402
from cemdisp.models2d.boundary_bridge import AnnulusInletState  # noqa: E402
from cemdisp.transport1d.casing_flow import CasingFlowSolver  # noqa: E402

WELL_LOADERS = {
    "呼101": ("cemdisp.data.loaders.hu101_loader", "load_hu101_tailpipe"),
    "呼1-003": ("cemdisp.data.loaders.ht1_003_loader", "load_ht1_003_tailpipe"),
    "呼1-004": ("cemdisp.data.loaders.ht1_004_loader", "load_ht1_004_tailpipe"),
}
# 静压锚：呼101 = 现场实测（notes 文本列 7868 m）；呼1-003 = 设计值（notes 文本列 7618 m）
ANCHORS_MPA = {"呼101": (152.16, "field_measured"), "呼1-003": (145.96, "design_value")}


@pytest.fixture(autouse=True)
def _clean_audit():
    reset_audit()
    yield
    reset_audit()


# --------------------------------------------------------------------------- #
# 0. 合成装配（小算例，供接线测试用）
# --------------------------------------------------------------------------- #

def _well_spec() -> WellSpec:
    top, bottom = 100.0, 400.0
    return WellSpec(
        well_name="P1_wiring",
        top_md_m=top, bottom_md_m=bottom, shoe_md_m=bottom,
        liner_id_mm=139.7, casing_id_mm=168.3,
        hole_diameter_profile=(DepthValuePoint(top, 260.0), DepthValuePoint(bottom, 250.0)),
        liner_od_profile=(DepthValuePoint(top, 168.3), DepthValuePoint(bottom, 168.3)),
        inclination_profile=(DepthValuePoint(top, 2.0), DepthValuePoint(bottom, 5.0)),
        standoff_profile=(DepthValuePoint(top, 0.75), DepthValuePoint(bottom, 0.60)),
    )


def _fluids() -> tuple:
    mud = FluidSpec("mud", FluidRole.MUD, 1050.0, RheologyModel.BINGHAM,
                    plastic_viscosity_pa_s=0.022, yield_stress_pa=6.0)
    spacer = FluidSpec("spacer", FluidRole.SPACER, 2050.0, RheologyModel.POWER_LAW,
                       power_law_n=0.7, consistency_k=0.3)
    lead = FluidSpec("lead", FluidRole.LEAD, 1890.0, RheologyModel.POWER_LAW,
                     power_law_n=0.75, consistency_k=0.55)
    tail = FluidSpec("tail", FluidRole.TAIL, 1920.0, RheologyModel.POWER_LAW,
                     power_law_n=0.70, consistency_k=0.90)
    return mud, spacer, lead, tail


def _schedule(fluids: tuple) -> PumpingSchedule:
    steps = tuple(
        PumpingScheduleStep(step_name=f"注{f.name}", fluid_name=f.name,
                            volume_m3=10.0, rate_m3_min=1.0)
        for f in fluids
    )
    return PumpingSchedule(steps=steps)


def _solver(**overrides) -> AnnulusD2DGASolver:
    params = dict(nz=24, ny=12, dt=2.0, total_t=240.0, enable_cfl_adaptive=False,
                  open_outlet=True)
    params.update(overrides)
    return AnnulusD2DGASolver(**params)


def _provider(t: float) -> AnnulusInletState:
    if t < 40.0:
        return AnnulusInletState(t, 0.01, "spacer", (("spacer", 1.0),))
    if t < 120.0:
        return AnnulusInletState(t, 0.01, "lead", (("lead", 1.0),))
    return AnnulusInletState(t, 0.01, "tail", (("tail", 1.0),))


def _load(well: str):
    mod, fn = WELL_LOADERS[well]
    t = getattr(importlib.import_module(mod), fn)()
    return t[0], t[1], t[2]


# --------------------------------------------------------------------------- #
# 1–3. 场本体
# --------------------------------------------------------------------------- #

def test_uniform_column_matches_analytic_rho_g_h():
    """均匀柱：``P(md) == g·ρ·TVD(md)/1e6``（解析对照，容差 <1e-12）。"""
    well = _well_spec()
    rho = 1900.0
    f = HydrostaticPressureField.from_well(well, rho)
    grid, cum = HydrostaticPressureField._tvd_axis(well)
    for md, tvd in zip(grid, cum):
        want = 9.80665 * rho * tvd / 1e6
        assert abs(f.P(float(md), 0.0) - want) < 1e-12, md
    assert f.P(400.0, 0.0) == f.P(400.0, 12345.0), "静压口径与 t 无关"
    assert f.P(0.0, 0.0) == 0.0


def test_layered_column_is_sandwiched_and_matches_hand_calc():
    """分层柱：与手算分段累加一致，且介于两端均匀柱之间。"""
    well = _well_spec()
    f = HydrostaticPressureField.from_layers(well, ((0.0, 1000.0), (250.0, 2000.0)))
    grid, cum = HydrostaticPressureField._tvd_axis(well)
    tvd_mid = float(np.interp(250.0, grid, cum))
    tvd_bot = float(np.interp(400.0, grid, cum))
    want = 9.80665 * (1000.0 * tvd_mid + 2000.0 * (tvd_bot - tvd_mid)) / 1e6
    assert abs(f.P(400.0, 0.0) - want) < 1e-12
    lo = HydrostaticPressureField.from_well(well, 1000.0).P(400.0, 0.0)
    hi = HydrostaticPressureField.from_well(well, 2000.0).P(400.0, 0.0)
    assert lo < f.P(400.0, 0.0) < hi


def test_layers_must_start_at_zero():
    with pytest.raises(ValueError):
        HydrostaticPressureField.from_layers(_well_spec(), ((10.0, 1000.0),))


def test_constant_field_and_oob_audit():
    f = ConstantPressureField(70.0, bottom_md_m=1000.0)
    assert f.P(500.0, 0.0) == 70.0 and f.P(500.0, 999.0) == 70.0
    assert f.oob_count == 0
    f.P(-1.0, 0.0)
    f.P(1001.0, 0.0)
    assert f.oob_count == 2 and len(f.oob_events) == 2
    f.reset_audit()
    assert f.oob_count == 0


def test_hydrostatic_oob_audit_and_endpoint_extrapolation():
    well = _well_spec()
    f = HydrostaticPressureField.from_well(well, 1900.0)
    p_bottom = f.P(well.bottom_md_m, 0.0)
    assert f.oob_count == 0
    assert f.P(well.bottom_md_m + 50.0, 0.0) == p_bottom, "越界按端点延拓"
    assert f.oob_count == 1
    assert f.P(-5.0, 0.0) == 0.0
    assert f.oob_count == 2


# --------------------------------------------------------------------------- #
# 4. 三井真实数据对账（A3.2）
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("well", ["呼101", "呼1-003", "呼1-004"])
def test_three_wells_insitu_density_and_rho_g_h_consistency(well):
    """口径③：``insitu_column_density`` 与手算体积加权逐位同；P(shoe) 与 ρ̄gh 解析一致。"""
    w, fluids, sched = _load(well)
    rho = insitu_column_density(fluids, sched)
    by_name = {f.name.strip(): f for f in fluids}
    vol = acc = 0.0
    for s in sched.steps:
        f = by_name.get(str(s.fluid_name).strip())
        if f is None or s.volume_m3 is None or float(s.volume_m3) <= 0.0:
            continue
        vol += float(s.volume_m3)
        acc += float(s.volume_m3) * float(f.density_kg_m3)
    assert rho == acc / vol, "必须与手算体积加权逐位相同"
    fld = HydrostaticPressureField.from_well(w, rho)
    _, cum = HydrostaticPressureField._tvd_axis(w)
    want = 9.80665 * rho * float(cum[-1]) / 1e6
    assert abs(fld.P(w.shoe_md_m, 0.0) - want) < 1e-9
    assert fld.oob_count == 0


@pytest.mark.parametrize("well", ["呼101", "呼1-003"])
def test_field_pressure_anchors_within_tolerance(well):
    """现场/设计静压锚 ±2%（A3.2）。两锚均只存在于 CSV 的 notes 文本列。"""
    w, fluids, sched = _load(well)
    fld = HydrostaticPressureField.from_well(w, insitu_column_density(fluids, sched))
    got = fld.P(w.shoe_md_m, 0.0)
    want, kind = ANCHORS_MPA[well]
    rel = abs(got - want) / want
    assert rel < 0.02, f"{well}({kind}): P(shoe)={got:.3f} vs 锚 {want} ⇒ {rel*100:.2f}%"


def test_ht1_003_design_anchor_exact_values_documented():
    """把呼1-003 的对账数字**显式钉住**（防口径悄悄漂移）。"""
    w, fluids, sched = _load("呼1-003")
    rho = insitu_column_density(fluids, sched)
    fld = HydrostaticPressureField.from_well(w, rho)
    got = fld.P(w.shoe_md_m, 0.0)
    assert abs(rho - 1949.5347) < 0.01, rho
    assert abs(got - 145.6294) < 0.002, got
    assert -0.005 < (got - 145.96) / 145.96 < 0.0, "设计锚**负偏**，量级 ~0.23%"


def test_annulus_roles_caliber_differs_slightly():
    """C-11：环空角色口径与全泵注口径的差异**已量化**（呼1-003 <0.05%）。"""
    _w, fluids, sched = _load("呼1-003")
    full = insitu_column_density(fluids, sched)
    ann = insitu_column_density(fluids, sched, roles=ANNULUS_ROLES)
    assert full != ann
    assert abs(full - ann) / full < 5e-4, (full, ann)


# --------------------------------------------------------------------------- #
# 5–9. 求解器接线
# --------------------------------------------------------------------------- #

def test_guan2_no_field_keeps_p_default_audit():
    """关2：无压力场 ⇒ `_representative_pressure_mpa` 返回 **None**（不是 0.0）。"""
    solver = _solver(enable_temperature_rheology=True)
    geom = solver._build_geom(_well_spec())
    assert solver._representative_pressure_mpa(geom, 0.0) is None
    _, spacer, _, _ = _fluids()
    reset_audit()
    derived = solver._phase_props(spacer, geom, 0.0)
    kinds = [e["kind"] for e in get_audit()]
    assert "p_default" in kinds, kinds
    assert kinds.count("clamp") == 0, "不得把 None 当 0.0 走 clamp"
    assert derived == fluid_at(spacer, 60.0)


def test_field_injection_changes_spacer_only():
    """传场后：隔离液 τy 随 P 变（含压二次曲面），其余族不受 P 影响。"""
    solver0 = _solver(enable_temperature_rheology=True)
    solver1 = _solver(enable_temperature_rheology=True,
                      pressure_field=ConstantPressureField(70.0))
    geom = solver0._build_geom(_well_spec())
    mud, spacer, lead, tail = _fluids()
    s0 = solver0._phase_props(spacer, geom, 0.0)
    s1 = solver1._phase_props(spacer, geom, 0.0)
    assert s1 is not s0 and s1.yield_stress_pa > s0.yield_stress_pa, \
        "隔离液 τy 随 P 增（C-13 方向锚）"
    for f in (mud, lead, tail):
        assert solver0._phase_props(f, geom, 0.0) == solver1._phase_props(f, geom, 0.0), \
            "水泥/钻井液公式纯温度 ⇒ P 不改变其派生"


def test_A3_3_memo_key_contains_P():
    """A3.3 memo 反污染：同 (fluid,T)、不同 P ⇒ **不同对象**，且各自等于手算。

    若 memo 键漏了 P，第二次查询会命中第一次的缓存 ⇒ 两值相等 ⇒ 必红。
    """
    p1, p2 = 10.0, 140.0
    solver = _solver(enable_temperature_rheology=True,
                     pressure_field=ConstantPressureField(p1))
    geom = solver._build_geom(_well_spec())
    _, spacer, _, _ = _fluids()
    a = solver._phase_props(spacer, geom, 0.0)
    assert a == fluid_at(spacer, 60.0, p1)
    solver.pressure_field = ConstantPressureField(p2)
    b = solver._phase_props(spacer, geom, 0.0)
    assert b == fluid_at(spacer, 60.0, p2)
    assert a is not b and a.yield_stress_pa != b.yield_stress_pa
    assert len([k for k in solver._phase_memo if k[0] is spacer]) == 2


def test_A3_1_construction_layer_equals_step_derivation():
    """A3.1：T-on + 传场的一次短 run 内，构造层暴露物性 ≡ 逐步派生（同一 P 档）。"""
    mud, spacer, lead, tail = _fluids()
    well = _well_spec()
    sched = _schedule((mud, spacer, lead, tail))
    fld = HydrostaticPressureField.from_well(well, 1900.0)
    solver = _solver(enable_temperature_rheology=True, pressure_field=fld)
    res = solver.run(well, (mud, spacer, lead, tail), _provider, schedule=sched)
    exposed = solver._temp_rheo_fluids
    p_expect = fld.P(float(well.shoe_md_m), 0.0)
    for name, f in (("mud", mud), ("spacer", spacer), ("lead", lead), ("tail", tail)):
        assert exposed[name] == fluid_at(f, 60.0, p_expect), name
    assert res.summary["effective_efficiency"] > 0.0
    p_vals = {k[2] for k in solver._phase_memo if k[0] is spacer}
    assert p_vals == {p_expect}, p_vals


def test_A3_2_run_with_field_has_zero_p_default():
    """A3.2：传场的 run 内 ``p_default`` 审计计数必须为 **0**。"""
    mud, spacer, lead, tail = _fluids()
    well = _well_spec()
    sched = _schedule((mud, spacer, lead, tail))
    fld = HydrostaticPressureField.from_well(well, 1900.0)
    solver = _solver(enable_temperature_rheology=True, pressure_field=fld)
    solver.run(well, (mud, spacer, lead, tail), _provider, schedule=sched)
    counts = solver._temp_rheo_audit_counts
    assert counts.get("p_default", 0) == 0, counts
    assert sum(counts.values()) > 0, "审计应非空（否则本守卫形同虚设）"


def test_casing_wiring_uses_shoe_scalar():
    """C-16：casing 侧 P 在 run 内定标为鞋深标量（`_phase_props(self, fluid)` 无 geom/t）。"""
    mud, spacer, lead, tail = _fluids()
    well = _well_spec()
    sched = _schedule((mud, spacer, lead, tail))
    fld = HydrostaticPressureField.from_well(well, 1900.0)
    solver = CasingFlowSolver(enable_temperature_rheology=True, pressure_field=fld)
    solver.run(well, (mud, spacer, lead, tail), sched)
    assert solver._step_P_mpa == fld.P(float(well.shoe_md_m), 0.0)
    assert solver._phase_props(spacer) == fluid_at(spacer, 60.0, solver._step_P_mpa)
    plain = CasingFlowSolver(enable_temperature_rheology=True)
    plain.run(well, (mud, spacer, lead, tail), sched)
    assert plain._step_P_mpa is None
    assert plain._phase_props(spacer) == fluid_at(spacer, 60.0)


def test_pressure_caliber_validation_and_mean_mode():
    """口径开关：非法值响亮报错；"mean" 档给出与 "shoe" 档**不同**的 P。"""
    with pytest.raises(ValueError):
        AnnulusD2DGASolver(dt=2.0, nz=8, ny=4, total_t=10.0, pressure_caliber="bad")
    with pytest.raises(ValueError):
        CasingFlowSolver(pressure_caliber="bad")
    hf = HydrostaticPressureField.from_well(_well_spec(), 1900.0)
    s_shoe = _solver(enable_temperature_rheology=True, pressure_field=hf)
    s_mean = _solver(enable_temperature_rheology=True, pressure_field=hf,
                     pressure_caliber="mean")
    geom = s_shoe._build_geom(_well_spec())
    p_shoe = s_shoe._representative_pressure_mpa(geom, 0.0)
    p_mean = s_mean._representative_pressure_mpa(geom, 0.0)
    assert p_shoe > p_mean > 0.0


def test_t_off_path_unaffected_by_pressure_field():
    """关2：T-off 时传压力场也**零痕迹**（`_phase_props` 恒等返回，不写 memo）。"""
    fld = ConstantPressureField(70.0)
    solver = _solver(pressure_field=fld)   # 默认 T-off
    geom = solver._build_geom(_well_spec())
    _, spacer, _, _ = _fluids()
    assert solver._phase_props(spacer, geom, 0.0) is spacer
    assert solver._phase_memo == {}
    assert "T" not in geom
    mud, spacer2, lead, tail = _fluids()
    well = _well_spec()
    solver2 = _solver(pressure_field=fld)
    res = solver2.run(well, (mud, spacer2, lead, tail), _provider,
                      schedule=_schedule((mud, spacer2, lead, tail)))
    assert np.isfinite(res.summary["effective_efficiency"])
