# -*- coding: utf-8 -*-
"""Task 6（A-3b）契约测试：`annulus_d2dga` 接线 `enable_hb_closure` + `hb_fix_cement_tau_y` 双开关。

覆盖（裁定逐条）：

1. **默认关逐位 = HEAD**（L1 硬约束）：三个新形参默认 False/None；显式传默认值与
   不传的求解器全结果**逐位一致**（float64 位级 sha256 摘要 + summary 逐项）。
2. **H1/H2/H3 语义可分离**（R-T6-1）：
   - H1（只开闭包）⇒ 混合 τy 场与 H0 **逐位一致**（屈服门不动）+ 非线性入口被消费；
   - H2（只补 τy）⇒ 混合 τy 场注入水泥常数 + 非线性入口**不**被消费（闭包保持牛顿）；
   - H3（两处同时）⇒ τy 场与 H2 一致 + 非线性入口被消费 + 闭包 τ_Y2 与屈服门通道
     **同源同值**（`_cement_phase_yield_stress` ≡ `_hb_closure_cement_tau_y`）。
3. **回退路径**（R-T5-3）：非线性外迭代 ``RuntimeError`` ⇒ 显式告警 + 回退牛顿线性
   路径，结果与纯牛顿路径**逐位一致**，且有告警证据。
4. **velocity_scale 接线**（R-T5-1 REVISED）：非线性入口收到 ``ŵ = q_half/π``、
   ``Gb=(0,0)``、**物理口径** (n, κ, τ_Y) 闭包（推荐路径，不走 Option B 等价式）。
5. **R4（τy 只受显式传入值）**：``hb_fix_cement_tau_y=True`` 无映射 ⇒ 构造期报错；
   映射缺相/MISSING 哨兵 ⇒ **告警跳过**（不静默用 0）；水泥 spec 自带 yield_stress_pa
   （Bingham，如 ht1_004）⇒ 以 spec 为准（常数不覆盖 spec）；非法值构造期拒绝。
6. **R7（红线）**：接线不触碰水泥相 ``FluidSpec`` 的流变模型/屈服字段（Task 0
   常数只经映射显式注入，水泥保持 ``POWER_LAW`` + ``yield_stress_pa is None``）。

测试用水泥相一律 ``POWER_LAW``（生产 8 井口径，spec 无 τy）——这正是 R-T6-1(a)
"水泥相贡献由 0 改为常数"的前提形态。常数用**测试内合成值**（接线测试不消费
Task 0 数据；该模块有自己的 provenance 契约测试）。
"""
from __future__ import annotations

import hashlib
import warnings

import numpy as np
import pytest

import cemdisp.models2d.annulus_d2dga as ann_mod
from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.well_spec import DepthValuePoint, WellSpec
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import AnnulusInletState


# --------------------------------------------------------------------------- #
# 公共装配（field-like 小域：Bingham 泥浆/隔离液 + POWER_LAW 水泥，屈服门默认开）
# --------------------------------------------------------------------------- #

# 生产级排量（呼探1 现场口径 1.0 m³/min ⇒ 环空均速 ~0.54 m/s，memory 09-11 实测锚）。
# ⚠️ 排量不能取小：Bingham 泥浆（PV=0.022/YP=6 Pa）在低流速下整层低于屈服门槛
# （τw < 1.5·τY），纯泥浆格点落在屈服悬崖上 ⇒ 反求不收敛（gap_solver 割线在
# 求积不连续处 2-循环）⇒ 触发 R-T5-3 回退。生产流速下泥浆格点 τw≈10 Pa > 9 Pa，
# 稳定越过悬崖——本文件判别测试因此用生产级排量。
_Q_M3S = 1.0 / 60.0


def _well_spec() -> WellSpec:
    top, bottom = 100.0, 400.0
    return WellSpec(
        well_name="T6_hb_wiring",
        top_md_m=top,
        bottom_md_m=bottom,
        shoe_md_m=bottom,
        hole_diameter_profile=(DepthValuePoint(top, 260.0), DepthValuePoint(bottom, 250.0)),
        liner_od_profile=(DepthValuePoint(top, 168.3), DepthValuePoint(bottom, 168.3)),
        inclination_profile=(DepthValuePoint(top, 2.0), DepthValuePoint(bottom, 5.0)),
        standoff_profile=(DepthValuePoint(top, 0.75), DepthValuePoint(bottom, 0.60)),
    )


def _fluids() -> tuple:
    """Bingham 泥浆/隔离液（带 YP，喂既有屈服门）+ POWER_LAW 水泥（spec 无 τy）。"""
    mud = FluidSpec("mud", FluidRole.MUD, 1050.0, RheologyModel.BINGHAM,
                    plastic_viscosity_pa_s=0.022, yield_stress_pa=6.0)
    spacer = FluidSpec("spacer", FluidRole.SPACER, 1120.0, RheologyModel.BINGHAM,
                       plastic_viscosity_pa_s=0.030, yield_stress_pa=9.0)
    lead = FluidSpec("lead", FluidRole.LEAD, 1890.0, RheologyModel.POWER_LAW,
                     power_law_n=0.75, consistency_k=0.55)
    tail = FluidSpec("tail", FluidRole.TAIL, 1920.0, RheologyModel.POWER_LAW,
                     power_law_n=0.70, consistency_k=0.90)
    return mud, spacer, lead, tail


def _provider():
    def provider(t: float) -> AnnulusInletState:
        if t < 40.0:
            return AnnulusInletState(t, _Q_M3S, "spacer", (("spacer", 1.0),))
        if t < 420.0:
            return AnnulusInletState(t, _Q_M3S, "lead", (("lead", 1.0),))
        return AnnulusInletState(t, _Q_M3S, "tail", (("tail", 1.0),))

    return provider


def _run_solver(**overrides):
    params = dict(nz=24, ny=12, dt=2.0, total_t=240.0, enable_cfl_adaptive=False,
                  open_outlet=True)
    params.update(overrides)
    solver = AnnulusD2DGASolver(**params)
    return solver.run(_well_spec(), _fluids(), _provider())


def _digest(arr: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(arr, dtype=np.float64).tobytes()).hexdigest()[:16]


def _result_digests(res) -> dict:
    s = res.summary["最终结果"]
    return {
        "cement": _digest(res.cement_field),
        "spacer": _digest(res.spacer_field),
        "wall": _digest(res.wall_field),
        "eta_E": float(s["全井段最终有效顶替效率"]),
        "eta_N": float(s["窄四分位效率"]),
        "mean_wall": float(res.metrics["mean_wall_mud"].iloc[-1]),
    }


_ORIGINAL_YIELD_GATE = AnnulusD2DGASolver._yield_gate_wall  # 原始 staticmethod 函数


def _patch_yield_gate_capture(monkeypatch, store: list) -> None:
    """捕获每步进屈服门的混合 τy 场（R-T6-1(a) 注入位），行为不改。

    用 staticmethod 包装避免 self 绑定；原函数在模块加载时捕获一次，
    防同测试内链式 patch 把上一个 spy 当成原函数。
    """
    def spy(w, b, mu_reg, tau_y, cement_ever, cement_local, f_safety):
        store.append(np.array(tau_y, dtype=float, copy=True))
        return _ORIGINAL_YIELD_GATE(w, b, mu_reg, tau_y, cement_ever, cement_local, f_safety)

    monkeypatch.setattr(AnnulusD2DGASolver, "_yield_gate_wall", staticmethod(spy))


def _patch_nonlinear_spy(monkeypatch, store: list) -> None:
    """包装真实非线性入口：记录调用与闭包对象并透传（判别「非线性入口是否被消费」）。"""
    real = ann_mod.solve_stream_function_nonlinear

    def spy(geom, c_bar, closure, b_field, *args, **kwargs):
        store.append({"closure": closure, "kwargs": dict(kwargs)})
        return real(geom, c_bar, closure, b_field, *args, **kwargs)

    monkeypatch.setattr(ann_mod, "solve_stream_function_nonlinear", spy)


# --------------------------------------------------------------------------- #
# 1. 默认关逐位 = HEAD（L1 硬约束）
# --------------------------------------------------------------------------- #

def test_switch_defaults_are_off():
    solver = AnnulusD2DGASolver()
    assert solver.enable_hb_closure is False
    assert solver.hb_fix_cement_tau_y is False
    assert solver.cement_tau_y_by_role is None


def test_default_off_is_bitwise_head():
    """三个新形参显式传默认值 ⇒ 与完全不传的求解器全结果逐位一致。"""
    base = _run_solver(total_t=600.0)
    explicit = _run_solver(total_t=600.0, enable_hb_closure=False,
                           hb_fix_cement_tau_y=False, cement_tau_y_by_role=None)
    assert _result_digests(base) == _result_digests(explicit)


# --------------------------------------------------------------------------- #
# 2. H1/H2/H3 语义可分离（R-T6-1）
# --------------------------------------------------------------------------- #

def _tau_y_field(solver, geom, ny, nz):
    """给定 solver 配置，在固定浓度场状态下调 `_compute_props`，返回混合 τy 场。

    τy 场由相体积加权而来（annulus_d2dga.py `_compute_props`），同一输入下
    「注入是否改变了 τy」可以确定性判别（跨 run 对比会受动力学差异污染）。
    返回序：mu, rho, mud, tau_y, ...（tau_y = 下标 3）。
    """
    lead = np.full((ny, nz), 0.3)
    tail = np.full((ny, nz), 0.2)
    spacer = np.full((ny, nz), 0.1)
    w_prev = np.full((ny, nz), 0.54)
    mud, spacer_f, lead_f, tail_f = _fluids()
    return solver._compute_props(lead, tail, spacer, w_prev, geom,
                                 mud, lead_f, tail_f, spacer_f)[3]


def test_h1_only_closure_leaves_yield_gate_field_bitwise_unchanged(monkeypatch):
    """H1（只开闭包）：同一浓度场状态下混合 τy 场与 H0 **逐位一致**（屈服门不动，
    R-T6-1）+ 非线性入口被消费。"""
    ny, nz = 12, 24
    solver_h0 = AnnulusD2DGASolver(nz=nz, ny=ny)
    solver_h1 = AnnulusD2DGASolver(nz=nz, ny=ny, enable_hb_closure=True)
    geom = solver_h0._build_geom(_well_spec())
    tau0 = _tau_y_field(solver_h0, geom, ny, nz)
    tau1 = _tau_y_field(solver_h1, geom, ny, nz)
    assert np.array_equal(tau0, tau1), "H1 不得改动混合 τy 场（屈服门不动，R-T6-1）"
    calls_h1: list = []
    _patch_nonlinear_spy(monkeypatch, calls_h1)
    _run_solver(total_t=120.0, enable_hb_closure=True)
    assert calls_h1, "H1 必须消费非线性外迭代入口（solve_stream_function_nonlinear）"
    assert "velocity_scale" in calls_h1[0]["kwargs"]


def test_h2_only_tau_y_injects_constant_and_keeps_newtonian_closure(monkeypatch):
    """H2（只补 τy）：同一浓度场状态下 τy 场注入水泥常数（≠H0，增量 = 份额×常数）
    + 非线性入口不被消费（闭包保持牛顿，R-T6-1）。"""
    ny, nz = 12, 24
    solver_h0 = AnnulusD2DGASolver(nz=nz, ny=ny)
    solver_h2 = AnnulusD2DGASolver(nz=nz, ny=ny, hb_fix_cement_tau_y=True,
                                   cement_tau_y_by_role={"LEAD": 8.0, "TAIL": 11.0})
    geom = solver_h0._build_geom(_well_spec())
    tau0 = _tau_y_field(solver_h0, geom, ny, nz)
    tau2 = _tau_y_field(solver_h2, geom, ny, nz)
    assert np.any(tau2 != tau0), "H2 必须改变混合 τy 场"
    lead = np.full((ny, nz), 0.3)
    tail = np.full((ny, nz), 0.2)
    np.testing.assert_allclose(tau2 - tau0, lead * 8.0 + tail * 11.0, rtol=0, atol=1e-12)
    calls_h2: list = []
    _patch_nonlinear_spy(monkeypatch, calls_h2)
    _run_solver(total_t=120.0, hb_fix_cement_tau_y=True,
                cement_tau_y_by_role={"LEAD": 8.0, "TAIL": 11.0})
    assert not calls_h2, "H2 闭包必须保持牛顿（非线性入口不得被消费，R-T6-1）"


def test_h3_both_channels_same_source_same_value():
    """H3（两处同时）：闭包 τ_Y2 与屈服门通道的水泥 τy **同源同值**（逐 role）。"""
    solver = AnnulusD2DGASolver(enable_hb_closure=True, hb_fix_cement_tau_y=True,
                                cement_tau_y_by_role={"LEAD": 8.0, "TAIL": 11.0})
    mud, spacer, lead, tail = _fluids()
    for fluid, expect in ((lead, 8.0), (tail, 11.0)):
        gate_value = solver._cement_phase_yield_stress(fluid)
        closure_value = solver._hb_closure_cement_tau_y(fluid)
        assert gate_value == expect
        assert closure_value == expect, "两处注入必须同源同值（R-T6-1）"


def test_h1_changes_dynamics_via_closure():
    """H1 的 HB 闭包（n₂≠1、τ_Y1=泥浆YP）改变算子 ⇒ 动力学结果与 H0 可分。

    同时验证 R-T6-2 验收：YP 口径 + 生产排量的整跑**不触发回退**（外迭代
    在 warm-start + 两段式冷启动 + 悬崖滞回下全程收敛）。
    """
    res_h0 = _run_solver(total_t=120.0)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        res_h1 = _run_solver(total_t=120.0, enable_hb_closure=True)
    fallbacks = [r for r in rec
                 if issubclass(r.category, RuntimeWarning) and "回退" in str(r.message)]
    assert not fallbacks, f"YP 口径生产工况不应回退：{[str(r.message)[:80] for r in fallbacks]}"
    d0, d1 = _result_digests(res_h0), _result_digests(res_h1)
    assert d0["cement"] != d1["cement"], "H1（HB 闭包）应改变浓度场"


def test_h2_reaches_dynamics_through_yield_gate_channel():
    """H2 经既有屈服门通道接通动力学：B-2（enable_stream_yield_gate）开启时 wall 进
    线性算子 ⇒ H2 与 H0 的浓度场可分（τy 注入真的到了动力学，非仅诊断量）。"""
    res_h0 = _run_solver(total_t=120.0, enable_stream_yield_gate=True)
    res_h2 = _run_solver(total_t=120.0, enable_stream_yield_gate=True,
                         hb_fix_cement_tau_y=True,
                         cement_tau_y_by_role={"LEAD": 8.0, "TAIL": 11.0})
    d0, d2 = _result_digests(res_h0), _result_digests(res_h2)
    assert d0["mean_wall"] != d2["mean_wall"] or d0["cement"] != d2["cement"], (
        "H2 的 τy 注入应经屈服门通道改变 wall/浓度场（B-2 开启时）")


# --------------------------------------------------------------------------- #
# 3. 回退路径（R-T5-3）：不收敛 ⇒ 告警 + 回退牛顿，结果与牛顿路径逐位一致
# --------------------------------------------------------------------------- #

def test_nonlinear_failure_falls_back_to_newtonian_with_warning(monkeypatch):
    """构造不收敛工况（非线性入口抛 RuntimeError）⇒ 显式告警 + 回退线性闭包；
    回退结果与纯牛顿路径逐位一致。"""
    calls = {"n": 0}

    def _boom(*args, **kwargs):
        calls["n"] += 1
        raise RuntimeError("非线性外迭代在 max_outer=50 轮内未收敛（模拟工况）")

    monkeypatch.setattr(ann_mod, "solve_stream_function_nonlinear", _boom)
    with pytest.warns(RuntimeWarning, match="回退"):
        res_hb = _run_solver(total_t=40.0, enable_hb_closure=True)
    assert calls["n"] >= 1, "回退前必须真的尝试过非线性入口"
    res_lin = _run_solver(total_t=40.0)
    assert _result_digests(res_hb) == _result_digests(res_lin), (
        "每步回退牛顿后，整跑结果必须与纯牛顿路径逐位一致")


# --------------------------------------------------------------------------- #
# 4. velocity_scale 接线（R-T5-1 REVISED）：ŵ = q_half/π + 物理 (n, κ, τ_Y)
# --------------------------------------------------------------------------- #

def test_velocity_scale_is_q_half_over_pi_and_closure_params_physical(monkeypatch):
    calls: list = []
    _patch_nonlinear_spy(monkeypatch, calls)
    _run_solver(total_t=40.0, enable_hb_closure=True)
    assert calls, "HB 路径未调用非线性入口"
    kw = calls[0]["kwargs"]
    # ŵ = q_half/π = (Q/2)/π（生产换算锚 annulus_d2dga w = w_unit*(q_half/np.pi)）
    assert kw["velocity_scale"] == pytest.approx((_Q_M3S / 2.0) / np.pi, rel=1e-12)
    # Gb 保持 (0,0)：Phase A 浮力全部走 b_field（不得两处算同一份浮力）
    closure = calls[0]["closure"]
    assert np.all(np.asarray(closure._Gb) == 0.0)
    # 物理 (n, κ, τ_Y)（推荐路径；H1 ⇒ τ_Y2=0；τ_Y1=泥浆 spec YP；κ/n 来自 FluidSpec）
    assert closure.n == (1.0, 0.75)          # 泥浆 Bingham→n=1；水泥 power_law_n
    assert closure.kappa == (0.022, 0.55)    # 泥浆 PV；水泥 consistency_k
    assert closure.tau_y == (6.0, 0.0)       # τ_Y1=泥浆 spec YP（R-T6-2）；τ_Y2=0（H1）


# --------------------------------------------------------------------------- #
# 5. R4：τy 只受显式传入值（缺失显式告警跳过 / 非法拒绝 / spec 优先）
# --------------------------------------------------------------------------- #

def test_hb_fix_requires_explicit_constants():
    """hb_fix_cement_tau_y=True 而无常数映射 ⇒ 构造期 ValueError（禁止代码编造）。"""
    with pytest.raises(ValueError, match="cement_tau_y_by_role"):
        AnnulusD2DGASolver(hb_fix_cement_tau_y=True)


def test_missing_constant_warns_and_skips_not_silent():
    """映射缺相 ⇒ 贡献 0 但**显式告警**（一次性汇总），不得静默。"""
    solver = AnnulusD2DGASolver(hb_fix_cement_tau_y=True, cement_tau_y_by_role={})
    _, _, lead, tail = _fluids()
    assert solver._cement_phase_yield_stress(lead) == 0.0
    with pytest.warns(UserWarning, match="LEAD"):
        solver._hb_flush_tau_y_skips()
    # 一次性：再次 flush 不再告警
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        solver._hb_flush_tau_y_skips()
    assert not [r for r in rec if issubclass(r.category, UserWarning)]


def test_missing_sentinel_string_is_skipped_with_warning():
    """MISSING 哨兵（字符串）⇒ 同样告警跳过（不静默、不抛错、不跨路线借值）。"""
    solver = AnnulusD2DGASolver(hb_fix_cement_tau_y=True,
                                cement_tau_y_by_role={"LEAD": "MISSING", "TAIL": 11.0})
    _, _, lead, tail = _fluids()
    assert solver._cement_phase_yield_stress(tail) == 11.0
    assert solver._cement_phase_yield_stress(lead) == 0.0
    with pytest.warns(UserWarning, match="LEAD"):
        solver._hb_flush_tau_y_skips()


def test_spec_yield_stress_takes_priority_over_constants():
    """水泥 spec 自带 yield_stress_pa（Bingham，如 ht1_004）⇒ 以 spec 为准，
    常数不覆盖 spec（且留痕告警）；H1 语义下闭包 τ_Y2 恒 0（裁定字面）。"""
    solver = AnnulusD2DGASolver(hb_fix_cement_tau_y=True,
                                cement_tau_y_by_role={"TAIL": 5.0})
    bingham_tail = FluidSpec("tail", FluidRole.TAIL, 1920.0, RheologyModel.BINGHAM,
                             plastic_viscosity_pa_s=0.17, yield_stress_pa=13.0)
    assert solver._cement_phase_yield_stress(bingham_tail) == 13.0
    with pytest.warns(UserWarning, match="自带"):
        solver._hb_flush_tau_y_skips()
    solver_h1 = AnnulusD2DGASolver(enable_hb_closure=True)
    assert solver_h1._hb_closure_cement_tau_y(bingham_tail) == 0.0


def test_invalid_constant_value_rejected():
    """映射值非负有限校验：负值/NaN 构造期拒绝。"""
    with pytest.raises(ValueError):
        AnnulusD2DGASolver(hb_fix_cement_tau_y=True, cement_tau_y_by_role={"TAIL": -1.0})
    with pytest.raises(ValueError):
        AnnulusD2DGASolver(hb_fix_cement_tau_y=True,
                           cement_tau_y_by_role={"TAIL": float("nan")})


def test_r7_cement_fluidspec_untouched_by_wiring():
    """R7 红线：接线不触碰水泥相 FluidSpec（流变模型/屈服字段与构造时逐项一致）。"""
    solver = AnnulusD2DGASolver(enable_hb_closure=True, hb_fix_cement_tau_y=True,
                                cement_tau_y_by_role={"LEAD": 8.0, "TAIL": 11.0})
    _, _, lead, tail = _fluids()
    assert lead.rheology_model is RheologyModel.POWER_LAW and lead.yield_stress_pa is None
    assert tail.rheology_model is RheologyModel.POWER_LAW and tail.yield_stress_pa is None
    solver._cement_phase_yield_stress(lead)
    solver._hb_closure_cement_tau_y(tail)
    assert lead.rheology_model is RheologyModel.POWER_LAW and lead.yield_stress_pa is None
    assert tail.rheology_model is RheologyModel.POWER_LAW and tail.yield_stress_pa is None


# --------------------------------------------------------------------------- #
# 6. 死开关告警（项目 A3 惯例）
# --------------------------------------------------------------------------- #

def test_dead_switch_warnings():
    """enable_hb_closure 需流函数路径；常数映射需 hb_fix；hb_fix 需屈服门或闭包。"""
    with pytest.warns(UserWarning, match="enable_hb_closure"):
        AnnulusD2DGASolver(enable_hb_closure=True, enable_stream_function=False)
    with pytest.warns(UserWarning, match="cement_tau_y_by_role"):
        AnnulusD2DGASolver(cement_tau_y_by_role={"TAIL": 5.0})
    with pytest.warns(UserWarning, match="hb_fix_cement_tau_y"):
        AnnulusD2DGASolver(hb_fix_cement_tau_y=True, cement_tau_y_by_role={"TAIL": 5.0},
                           enable_yield_gate=False)
    # 合法组合不告警
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        AnnulusD2DGASolver(enable_hb_closure=True, hb_fix_cement_tau_y=True,
                           cement_tau_y_by_role={"TAIL": 5.0})
    assert not [r for r in rec if issubclass(r.category, UserWarning)]


# --------------------------------------------------------------------------- #
# 7. R-T6-2 修复轮 1：warm-start 默认逐位 + 汇总尾缀分流（评审 Important-1）
# --------------------------------------------------------------------------- #

def test_warm_start_default_is_bitwise_cold_start():
    """solve_stream_function_nonlinear 的 initial_G 默认 None ⇒ 与不传参逐位一致。"""
    from cemdisp.models2d.hb_closure import HBClosure
    from cemdisp.models2d.stream_function import (
        solve_stream_function, solve_stream_function_nonlinear)
    y = np.linspace(0.0, np.pi * 0.1071, 12)
    phi = y / y[-1]
    nz = 6
    H = np.broadcast_to(0.0458 * (1 + 0.3 * np.cos(np.pi * phi))[:, None], (12, nz)).copy()
    geom = {"y": y, "phi": phi, "H": H, "b": 2 * H,
            "s": np.linspace(0.0, 10.0, nz),
            "hole_mm": np.full((1, nz), 260.0), "od_mm": np.full((1, nz), 168.3)}
    c = 0.5 * np.ones((12, nz))
    b = np.zeros((2, 12, nz))
    b[0] = 8.0e4 * np.sin(np.pi * phi)[:, None] * c
    hb = HBClosure(n=(1.0, 0.8), kappa=(0.022, 0.55), tau_y=(6.0, 0.0),
                   m=0.4, B=0.3, ny_gap=41)
    hb.set_pressure_gradient(300.0)  # 冷启动语义：调用方注入初值
    psi_a = solve_stream_function_nonlinear(geom, c, hb, b, omega=0.5,
                                            tol=1.0, max_outer=3)
    hb2 = HBClosure(n=(1.0, 0.8), kappa=(0.022, 0.55), tau_y=(6.0, 0.0),
                    m=0.4, B=0.3, ny_gap=41)
    hb2.set_pressure_gradient(300.0)
    psi_b = solve_stream_function_nonlinear(geom, c, hb2, b, initial_G=None,
                                            omega=0.5, tol=1.0, max_outer=3)
    assert np.array_equal(psi_a, psi_b)


def test_flush_summary_suffix_split_by_kind():
    """汇总告警尾缀按 skip 类型分流：missing ⇒ 贡献 0 提示；spec ⇒ 贡献非 0 声明。"""
    _, _, lead, tail = _fluids()
    solver_missing = AnnulusD2DGASolver(hb_fix_cement_tau_y=True,
                                        cement_tau_y_by_role={})
    solver_missing._cement_phase_yield_stress(lead)
    with pytest.warns(UserWarning) as rec_missing:
        solver_missing._hb_flush_tau_y_skips()
    msg_missing = str(rec_missing[0].message)
    assert "LEAD" in msg_missing
    assert "贡献为 0" in msg_missing
    assert "贡献非 0" not in msg_missing

    solver_spec = AnnulusD2DGASolver(hb_fix_cement_tau_y=True,
                                     cement_tau_y_by_role={"TAIL": 5.0})
    bingham_tail = FluidSpec("tail", FluidRole.TAIL, 1920.0, RheologyModel.BINGHAM,
                             plastic_viscosity_pa_s=0.17, yield_stress_pa=13.0)
    solver_spec._cement_phase_yield_stress(bingham_tail)
    with pytest.warns(UserWarning) as rec_spec:
        solver_spec._hb_flush_tau_y_skips()
    msg_spec = str(rec_spec[0].message)
    assert "TAIL" in msg_spec
    assert "贡献非 0" in msg_spec
    assert "贡献为 0" not in msg_spec
