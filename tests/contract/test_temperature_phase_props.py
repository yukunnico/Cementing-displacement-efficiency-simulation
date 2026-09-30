# -*- coding: utf-8 -*-
"""Task 7（T1-2）契约测试：2D 温度挂接——`_phase_props` 物性唯一入口、
geom["T"] 每步刷新、屈服门 τy(T)、HB memo key 温度维、Tier0 派生口径。

覆盖（task-7-brief 需求逐条）：

1. **`_phase_props` 两调用点一致性**：`_compute_props`（场口径，含旧代数旁路
   依赖的同一函数）与 `_velocity_stream_function`（标量代表剪切率口径）同走
   `_phase_props`；同输入 ⇒ 同对象同输出（memo 消灭双算），且派生值 ==
   手写 ``fluid_at(base, T)``。
2. **memo key 隔离（T 维）**：T-on 时 HB memo key 追加本步代表温度，
   T1≠T2 ⇒ 不同 key；T-off 时 key 为 6 元组（不追加温度 ⇒ 序列与 HEAD 一致）。
3. **屈服门 τy 口径**：T-on 吃公式 τy（``fluid_at`` 派生的 ``yield_stress_pa``，
   绝对替换优先于 ``cement_yield_stress`` 常数机制）；T-off 照旧吃常数机制
   （``hb_fix_cement_tau_y=True`` + ``cement_tau_y_by_role`` 映射）。
4. **geom["T"]**：T-on 有键且每步刷新（注入场对象生效）；T-off 无键（零痕迹）。
5. **Tier0 口径**：T-on 收派生后流体（诊断所见 = 求解所用）；T-off 收原元组
   （同一对象）。

T-off 逐位红线由 `.superpowers/sdd/温压耦合改进计划_2026-09-30/task-7-toff-*`
探针（byte 级 diff）承担；本文件只做 T-on/接口级断言。
"""
from __future__ import annotations

import numpy as np

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.rheology_vs_temperature import fluid_at
from cemdisp.data.well_spec import DepthValuePoint, WellSpec
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import AnnulusInletState

# 生产级排量（同 test_temperature_rheology_switch）。
_Q_M3S = 1.0 / 60.0


# --------------------------------------------------------------------------- #
# 公共装配（与 Task 6 契约测试同源，便于跨任务对照）
# --------------------------------------------------------------------------- #

def _well_spec() -> WellSpec:
    top, bottom = 100.0, 400.0
    return WellSpec(
        well_name="T12_phase_props",
        top_md_m=top,
        bottom_md_m=bottom,
        shoe_md_m=bottom,
        hole_diameter_profile=(DepthValuePoint(top, 260.0), DepthValuePoint(bottom, 250.0)),
        liner_od_profile=(DepthValuePoint(top, 168.3), DepthValuePoint(bottom, 168.3)),
        inclination_profile=(DepthValuePoint(top, 2.0), DepthValuePoint(bottom, 5.0)),
        standoff_profile=(DepthValuePoint(top, 0.75), DepthValuePoint(bottom, 0.60)),
    )


def _fluids() -> tuple:
    """Bingham 泥浆/隔离液 + POWER_LAW 水泥（spec 无 τy ⇒ 公式/常数可分辨）。"""
    mud = FluidSpec("mud", FluidRole.MUD, 1050.0, RheologyModel.BINGHAM,
                    plastic_viscosity_pa_s=0.022, yield_stress_pa=6.0)
    spacer = FluidSpec("spacer", FluidRole.SPACER, 1120.0, RheologyModel.BINGHAM,
                       plastic_viscosity_pa_s=0.030, yield_stress_pa=9.0)
    lead = FluidSpec("lead", FluidRole.LEAD, 1890.0, RheologyModel.POWER_LAW,
                     power_law_n=0.75, consistency_k=0.55)
    tail = FluidSpec("tail", FluidRole.TAIL, 1920.0, RheologyModel.POWER_LAW,
                     power_law_n=0.70, consistency_k=0.90)
    return mud, spacer, lead, tail


def _provider(t: float) -> AnnulusInletState:
    if t < 40.0:
        return AnnulusInletState(t, _Q_M3S, "spacer", (("spacer", 1.0),))
    if t < 420.0:
        return AnnulusInletState(t, _Q_M3S, "lead", (("lead", 1.0),))
    return AnnulusInletState(t, _Q_M3S, "tail", (("tail", 1.0),))


def _solver(**overrides) -> AnnulusD2DGASolver:
    params = dict(nz=24, ny=12, dt=2.0, total_t=240.0, enable_cfl_adaptive=False,
                  open_outlet=True)
    params.update(overrides)
    return AnnulusD2DGASolver(**params)


def _zeros(solver: AnnulusD2DGASolver, fill: float = 0.0) -> np.ndarray:
    return np.full((solver.ny, solver.nz), fill, dtype=float)


class _StubField:
    """斜坡温度场（°C = 10 + 0.01·md + 0.001·t）：验证注入接口与逐步刷新。"""

    def T(self, md_m: float, t_s: float) -> float:
        return 10.0 + 0.01 * float(md_m) + 0.001 * float(t_s)


# --------------------------------------------------------------------------- #
# 1. `_phase_props` 两调用点一致性（同输入同输出 + 同走同一函数 + 旁路同步）
# --------------------------------------------------------------------------- #

def test_phase_props_off_is_identity():
    """T-off：恒等返回入参，不派生、不写温度状态（off 红线）。"""
    solver = _solver()  # 默认关
    mud, spacer, lead, tail = _fluids()
    geom = solver._build_geom(_well_spec())
    assert "T" not in geom, "T-off 的 geom 不得出现 T 键（零痕迹）"
    assert solver._phase_props(lead, geom, 0.0) is lead
    assert solver._phase_props(None, geom, 0.0) is None
    assert solver._step_T_c is None, "T-off 不得写入派生温度"
    assert solver._phase_memo == {}


def test_phase_props_on_matches_manual_fluid_at():
    """T-on + 恒温 60°C：派生值逐字段等于手写 ``fluid_at(base, 60)``（自洽）。"""
    solver = _solver(enable_temperature_rheology=True)
    mud, spacer, lead, tail = _fluids()
    geom = solver._build_geom(_well_spec())  # T-on ⇒ geom["T"] ≡ 60（恒温场）
    assert np.all(geom["T"] == 60.0)
    got = solver._phase_props(lead, geom, 0.0)
    assert got == fluid_at(lead, 60.0)
    assert got is not lead
    assert solver._step_T_c == 60.0, "派生温度须记入 _step_T_c（HB memo key 消费）"


def test_two_call_sites_share_phase_props():
    """`_compute_props`（场口径）与 `_velocity_stream_function`（标量口径）
    同走 `_phase_props`，且同输入 ⇒ **同一派生对象**（memo 消灭双算）。"""
    solver = _solver(enable_temperature_rheology=True)
    mud, spacer, lead, tail = _fluids()
    geom = solver._build_geom(_well_spec())

    returned: list = []
    orig = solver._phase_props

    def spy(fluid, g, t):
        out = orig(fluid, g, t)
        returned.append((fluid, out))
        return out

    solver._phase_props = spy
    zero = _zeros(solver)
    w_prev = _zeros(solver, 0.45)

    # 调用点 A：场口径（_compute_props）
    solver._compute_props(zero, zero, zero, w_prev, geom, mud, lead, tail, spacer)
    n_after_a = len(returned)
    assert n_after_a >= 4, f"_compute_props 应逐相走 _phase_props，实际 {n_after_a} 次"

    # 调用点 B：流函数标量口径
    solver._velocity_stream_function(zero, zero, geom, _Q_M3S, w_prev, mud, lead, tail)
    assert len(returned) > n_after_a, "流函数段未走 _phase_props"

    # 同输入同输出：两调用点对同一 base lead 得到同一对象，且 == 手写参照
    lead_outs = [out for (f, out) in returned if f is lead]
    assert len(lead_outs) >= 2
    assert all(out is lead_outs[0] for out in lead_outs), (
        "两调用点对同一 (fluid, T) 必须复用同一派生对象（memo）")
    assert lead_outs[0] == fluid_at(lead, 60.0)


def test_old_path_bypass_also_routes_through_phase_props():
    """旁路（enable_stream_function=False 旧代数路径）同步走 `_phase_props`
    ——经 `_compute_props` 同一入口（含旁路必须同步的裁定）。"""
    solver = _solver(enable_temperature_rheology=True, enable_stream_function=False)
    mud, spacer, lead, tail = _fluids()
    geom = solver._build_geom(_well_spec())

    calls: list = []
    orig = solver._phase_props

    def spy(fluid, g, t):
        calls.append(fluid)
        return orig(fluid, g, t)

    solver._phase_props = spy
    zero = _zeros(solver)
    w_prev = _zeros(solver, 0.45)
    solver._compute_velocity(zero, zero, zero, geom, _Q_M3S, w_prev,
                             mud, lead, tail, spacer)
    assert len(calls) >= 4, "旁路路径的 _compute_props 未走 _phase_props"


# --------------------------------------------------------------------------- #
# 2. HB memo key 温度维
# --------------------------------------------------------------------------- #

def test_hb_memo_key_isolated_by_temperature():
    """T-on：key 追加本步代表温度 ⇒ T1≠T2 不串缓存；T-off：key 不追加（6 元组）。"""
    mud, spacer, lead, tail = _fluids()

    on = _solver(enable_temperature_rheology=True, enable_hb_closure=True)
    on._step_T_c = 60.0
    on._hb_closure_for(mud, lead, 10.0)
    on._step_T_c = 90.0
    on._hb_closure_for(mud, lead, 10.0)
    keys_on = list(on._hb_closure_memo)
    assert len(keys_on) == 2, "不同 T 必须落在不同 memo 条目（不得串缓存）"
    assert all(len(k) == 7 for k in keys_on), f"T-on key 应为 6 参数 + 1 温度，得到 {[len(k) for k in keys_on]}"
    assert {k[-1] for k in keys_on} == {60.0, 90.0}

    # 同 T 复用
    on._step_T_c = 90.0
    on._hb_closure_for(mud, lead, 10.0)
    assert len(on._hb_closure_memo) == 2

    off = _solver(enable_hb_closure=True)  # 默认关
    off._hb_closure_for(mud, lead, 10.0)
    keys_off = list(off._hb_closure_memo)
    assert len(keys_off) == 1
    assert len(keys_off[0]) == 6, (
        f"T-off key 序列必须与 HEAD 一致（6 元组无温度维），得到 {keys_off[0]}")


# --------------------------------------------------------------------------- #
# 3. 屈服门 τy 口径：T-on 公式绝对替换 / T-off 常数机制
# --------------------------------------------------------------------------- #

def test_cement_tauy_formula_on_constant_off():
    """水泥相 τy：T-on 走公式（fluid_at 派生 yield_stress_pa，绝对替换优先）；
    T-off + hb_fix_cement_tau_y=True 照旧吃 cement_tau_y_by_role 常数。"""
    mapping = {"LEAD": 5.0, "TAIL": 7.0}
    mud, spacer, lead, tail = _fluids()
    expected_formula = fluid_at(lead, 60.0).yield_stress_pa
    assert expected_formula is not None and abs(expected_formula - 5.0) > 1e-6, (
        "装配前提：公式 τy(60°C) 须与常数映射可分辨")

    # ---- T-off：常数机制（hb_fix=True ⇒ cement_tau_y_by_role 映射） ----
    off = _solver(hb_fix_cement_tau_y=True, cement_tau_y_by_role=mapping)
    geom_off = off._build_geom(_well_spec())
    zero, one = _zeros(off), _zeros(off, 1.0)
    _, _, _, tauy_off, *_ = off._compute_props(
        one, zero, zero, _zeros(off, 0.45), geom_off, mud, lead, tail, None)
    assert np.allclose(tauy_off, 5.0), (
        f"T-off 应吃常数映射 5.0，得到 {float(np.mean(tauy_off))}")

    # ---- T-on：公式绝对替换（同一常数映射被旁路） ----
    on = _solver(hb_fix_cement_tau_y=True, cement_tau_y_by_role=mapping,
                 enable_temperature_rheology=True)
    geom_on = on._build_geom(_well_spec())
    _, _, _, tauy_on, *_ = on._compute_props(
        one, zero, zero, _zeros(on, 0.45), geom_on, mud, lead, tail, None)
    assert np.allclose(tauy_on, expected_formula), (
        f"T-on 应吃公式 τy={expected_formula}，得到 {float(np.mean(tauy_on))}")
    assert not np.allclose(tauy_on, 5.0), "T-on 不得退回 cement_yield_stress 常数机制"


# --------------------------------------------------------------------------- #
# 4. geom["T"]：初始化 + 每步刷新 + 注入场对象 + T-off 无键
# --------------------------------------------------------------------------- #

def test_geom_T_initialized_and_refreshed():
    solver = _solver(enable_temperature_rheology=True)
    solver._temperature_field = _StubField()  # 注入场对象（run() 同款接口）
    geom = solver._build_geom(_well_spec())
    md = geom["md"]
    assert "T" in geom and geom["T"].shape == md.shape
    assert np.allclose(geom["T"], 10.0 + 0.01 * md), "t=0 初始化应等于场的初态剖面"

    solver._refresh_geom_temperature(geom, 600.0)
    assert np.allclose(geom["T"], 10.0 + 0.01 * md + 0.6), "每步刷新未生效"
    # 代表温度 = 非均匀剖面的深度均值
    rep = solver._representative_temperature(geom, 600.0)
    assert abs(rep - float(np.mean(geom["T"]))) < 1e-12


def test_run_writes_geom_T_only_when_on():
    """run 级：T-on 结果 geom 含 T 剖面；T-off 结果 geom 无 T 键（零痕迹）。"""
    off = _solver()
    res_off = off.run(_well_spec(), _fluids(), _provider)
    assert "T" not in res_off.geom

    on = _solver(enable_temperature_rheology=True)
    res_on = on.run(_well_spec(), _fluids(), _provider, temperature_field=_StubField())
    assert "T" in res_on.geom and res_on.geom["T"].shape == (on.nz,)
    # 末步刷新时刻 = total_t（固定 dt：末次迭代 current_time_s = 240）
    expected = 10.0 + 0.01 * res_on.geom["md"] + 0.001 * 240.0
    assert np.allclose(res_on.geom["T"], expected), (
        "geom['T'] 应为注入温度场在末步时刻的剖面（逐步刷新接线）")


# --------------------------------------------------------------------------- #
# 5. Tier0 诊断口径：T-on 收派生后流体 / T-off 收原元组
# --------------------------------------------------------------------------- #

def test_tier0_fluids_follow_temperature_switch(monkeypatch):
    import cemdisp.diagnostics.tier0_diagnostics as t0mod

    captured: list = []
    orig = t0mod.compute_all_tier0_diagnostics

    def spy(result, **kw):
        captured.append(kw.get("fluids"))
        return orig(result, **kw)

    monkeypatch.setattr(t0mod, "compute_all_tier0_diagnostics", spy)

    fluids_off = _fluids()
    _solver().run(_well_spec(), fluids_off, _provider)
    assert captured[-1] is fluids_off, "T-off 必须收原 fluids 元组（同一对象，零变化）"

    fluids_on = _fluids()
    _solver(enable_temperature_rheology=True).run(_well_spec(), fluids_on, _provider)
    diag_on = captured[-1]
    assert diag_on is not fluids_on, "T-on 不得收原始 fluids 元组"
    assert len(diag_on) == len(fluids_on)
    base_by_name = {f.name: f for f in fluids_on}
    for f in diag_on:
        if f.role in (FluidRole.LEAD, FluidRole.TAIL, FluidRole.MUD, FluidRole.SPACER):
            assert f is not base_by_name[f.name], f"{f.name} 应为派生后流体"
            assert f == fluid_at(base_by_name[f.name], 60.0), (
                f"{f.name} 应与求解同口径的 fluid_at(base, 60°C) 一致")


# --------------------------------------------------------------------------- #
# 6. 开关开 = 进入动力学（回归护栏：T1-2 改动不得把开关改空转）
# --------------------------------------------------------------------------- #

def test_on_still_changes_dynamics_and_off_stays_clean():
    res_off = _solver().run(_well_spec(), _fluids(), _provider)
    res_on = _solver(enable_temperature_rheology=True).run(_well_spec(), _fluids(), _provider)
    assert not hasattr(_solver(), "_temp_rheo_fluids")  # 关=零痕迹（新实例）
    assert "temperature_rheology_audit" not in res_off.summary
    assert res_off.cement_field.tobytes() != res_on.cement_field.tobytes(), (
        "开/关浓度场必须不同（开关真被消费）")
