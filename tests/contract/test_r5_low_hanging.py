# -*- coding: utf-8 -*-
"""Phase 5a 契约测试：R5 低垂三项代码修复（2026-10-07）

上游：`docs/superpowers/specs/2026-10-07-phase5a-r5-low-hanging-design.md`

覆盖
----
1. **5a-②** `rheology_vs_temperature.is_replaced`：覆盖判定与 `fluid_at` 的
   「原样返回」判据**同源**（no_replace 族 / 水泥密度越界 ⇒ False；其余 ⇒ True），
   并以对象恒等 `fluid_at(f, T) is f` 交叉钉住。
2. **5a-②** `_phase_cement_tau_y` 组合语义：**覆盖优先、未覆盖回落**；
   `hb_fix_cement_tau_y=False`（生产默认）时两支**同值**（关2 逐位依据）。
3. **5a-③** `CasingFlowSolver._effective_viscosity`：宾汉返回 `τy/γ̇ + PV`
   （PV 被 0.01 顶掉的地雷已修）；PV 缺省时才回退 0.01。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[2]
for _p in (str(_ROOT), str(_ROOT / "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel  # noqa: E402
from cemdisp.data.rheology_vs_temperature import fluid_at, is_replaced  # noqa: E402
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver  # noqa: E402
from cemdisp.transport1d.casing_flow import CasingFlowSolver  # noqa: E402


def _fluid(name: str, role: FluidRole, *, rho: float = 1900.0,
           tauy: float | None = 10.0, pv: float | None = 0.05) -> FluidSpec:
    return FluidSpec(
        name=name, role=role, density_kg_m3=rho,
        rheology_model=RheologyModel.BINGHAM,
        plastic_viscosity_pa_s=pv, yield_stress_pa=tauy,
    )


# --------------------------------------------------------------------------- #
# 1. is_replaced 与 fluid_at 同源
# --------------------------------------------------------------------------- #
def test_is_replaced_matches_fluid_at_identity():
    """`is_replaced(f) is False` ⇔ `fluid_at(f, T) is f`（同源判据，逐例交叉）。"""
    cases = [
        _fluid("钻井液", FluidRole.MUD, rho=1200.0),            # mud ⇒ True
        _fluid("隔离液", FluidRole.SPACER, rho=2050.0),          # spacer ⇒ True
        _fluid("领浆", FluidRole.LEAD, rho=1900.0),              # cement 档内 ⇒ True
        _fluid("尾浆", FluidRole.TAIL, rho=2400.0),              # cement 档外 ⇒ False
        _fluid("冲洗液", FluidRole.FLUSHER, rho=1000.0),         # no_replace ⇒ False
        _fluid("压塞液", FluidRole.DISPLACEMENT, rho=1000.0),    # no_replace ⇒ False
    ]
    for f in cases:
        replaced = is_replaced(f)
        assert (fluid_at(f, 60.0) is f) is (not replaced), f"{f.name} 判据不同源"


def test_is_replaced_cement_density_gate():
    """水泥族只看密度捕获区：域内 True、域外 False（与 _cement_values 的 None 条件同）。"""
    assert is_replaced(_fluid("领浆", FluidRole.LEAD, rho=1900.0)) is True
    assert is_replaced(_fluid("领浆", FluidRole.LEAD, rho=2100.0)) is True
    assert is_replaced(_fluid("领浆", FluidRole.LEAD, rho=999.0)) is False
    assert is_replaced(_fluid("领浆", FluidRole.LEAD, rho=3000.0)) is False


def test_is_replaced_returns_bool():
    assert isinstance(is_replaced(_fluid("钻井液", FluidRole.MUD, rho=1200.0)), bool)


# --------------------------------------------------------------------------- #
# 2. _phase_cement_tau_y 组合语义
# --------------------------------------------------------------------------- #
def _solver(**kw) -> AnnulusD2DGASolver:
    return AnnulusD2DGASolver(total_t=10.0, nz=10, ny=6,
                              enable_temperature_rheology=True, **kw)


def test_uncovered_cement_falls_back_and_default_branches_agree():
    """未覆盖水泥相：T-on 回落分支；`hb_fix=False`（默认）时两支**同值**（关2）。"""
    s = _solver()                                   # hb_fix_cement_tau_y 默认 False
    out_of_tier = _fluid("尾浆", FluidRole.TAIL, rho=2400.0, tauy=17.5)
    assert is_replaced(out_of_tier) is False
    got = s._phase_cement_tau_y(out_of_tier)
    # 默认 hb_fix=False ⇒ `_cement_phase_yield_stress == _fluid_yield_stress`
    assert got == s._fluid_yield_stress(out_of_tier) == 17.5


def test_covered_cement_still_uses_formula_value():
    """被覆盖相：仍走公式派生态的 τy（绝对替换裁定不变）。"""
    s = _solver()
    derived = _fluid("领浆", FluidRole.LEAD, rho=1900.0, tauy=42.0)
    assert is_replaced(derived) is True
    assert s._phase_cement_tau_y(derived) == 42.0


def test_uncovered_cement_uses_constant_when_hb_fix_on():
    """`hb_fix=True` + T-on + **未覆盖**相 ⇒ 常数回落（修复目标；原实现会丢）。"""
    s = _solver(hb_fix_cement_tau_y=True, cement_tau_y_by_role={"TAIL": 33.0})
    # spec 无 yield_stress_pa（FluidSpec 对 Bingham 强制要求 τy ⇒ 用幂律相表达"无该字段"）
    out_of_tier = FluidSpec(name="尾浆", role=FluidRole.TAIL, density_kg_m3=2400.0,
                            rheology_model=RheologyModel.POWER_LAW,
                            power_law_n=0.7, consistency_k=0.5)
    assert is_replaced(out_of_tier) is False
    assert s._phase_cement_tau_y(out_of_tier) == pytest.approx(33.0)


# --------------------------------------------------------------------------- #
# 3. 1D 宾汉丢 PV 地雷
# --------------------------------------------------------------------------- #
def _casing(**kw) -> CasingFlowSolver:
    return CasingFlowSolver(enable_gravity=True, has_plug=True, **kw)


def test_bingham_pv_fix_default_off_is_bitwise_head():
    """**默认关**：`bingham_pv_fix=False` ⇒ 沿用 τy/γ̇ + 0.01（逐位=HEAD 红线）。

    该默认值不是笔误：本入口在 `has_plug=False` 配置下可达（经
    `_interface_instability_factor` → 弥散带增强），直接改值会位移
    `tests/history/test_casing_mixing_contact_time.py::test_default_off_bitwise`
    的 hu103 默认路径冻结锚。故按铁律「一切 opt-in」，修复默认关。
    """
    s = _casing()
    f = _fluid("水泥", FluidRole.TAIL, tauy=20.0, pv=0.25)
    U, R = 1.0, 0.1
    gamma = 8.0 * U / (2.0 * R)
    assert s.bingham_pv_fix is False
    assert s._effective_viscosity(f, U, R) == pytest.approx(20.0 / gamma + 0.01)


def test_bingham_effective_viscosity_keeps_pv_when_enabled():
    """开关置真：宾汉 μ_eff = τy/γ̇ + PV（γ̇ = 8U/(2R)）；PV 不再被 0.01 顶掉。"""
    s = _casing(bingham_pv_fix=True)
    f = _fluid("水泥", FluidRole.TAIL, tauy=20.0, pv=0.25)
    U, R = 1.0, 0.1
    gamma = 8.0 * U / (2.0 * R)
    assert s._effective_viscosity(f, U, R) == pytest.approx(20.0 / gamma + 0.25)


def test_hb_effective_viscosity_unchanged():
    """HB 支不动：μ_eff = τy/γ̇ + K·γ̇^(n-1)（仍走 consistency_k，不误用 PV）。"""
    s = _casing()
    f = FluidSpec(name="尾浆", role=FluidRole.TAIL, density_kg_m3=1900.0,
                  rheology_model=RheologyModel.HERSCHEL_BULKLEY,
                  yield_stress_pa=20.0, plastic_viscosity_pa_s=0.9,
                  consistency_k=0.4, power_law_n=0.8)
    U, R = 1.0, 0.1
    gamma = 8.0 * U / (2.0 * R)
    assert s._effective_viscosity(f, U, R) == pytest.approx(20.0 / gamma + 0.4 * gamma ** (-0.2))


def test_power_law_effective_viscosity_unchanged():
    """幂律支不动：μ_eff = K·γ̇^(n-1)（τy=0）。"""
    s = _casing()
    f = FluidSpec(name="钻井液", role=FluidRole.MUD, density_kg_m3=1200.0,
                  rheology_model=RheologyModel.POWER_LAW, power_law_n=0.6,
                  consistency_k=0.8)
    U, R = 1.0, 0.1
    gamma = 8.0 * U / (2.0 * R)
    assert s._effective_viscosity(f, U, R) == pytest.approx(0.8 * gamma ** (-0.4))


# --------------------------------------------------------------------------- #
# 4. 5a-① 停泵衰减判据（Moyers-González & Frigaard）
# --------------------------------------------------------------------------- #
from cemdisp.diagnostics.stop_pump_freeze import (  # noqa: E402
    StopPumpFreezeResult,
    compute_stop_pump_freeze,
    fields_from_state,
)


def test_freeze_verdict_frozen_when_yield_dominates():
    """屈服抗力 >> 浮力驱动 ⇒ 冻结（判据方向）。"""
    ny, nz = 4, 3
    tau = np.full((ny, nz), 50.0)          # 50 Pa
    rho = np.zeros((ny, nz))
    rho[0, :] = 1500.0                     # 方位向差 100 kg/m³
    rho[-1, :] = 1600.0
    r = compute_stop_pump_freeze(
        tau_y_field=tau, rho_field=rho, gap_m=np.full(nz, 0.02),
        standoff=np.full(nz, 0.8), beta_deg=30.0)
    assert isinstance(r, StopPumpFreezeResult)
    assert r.frozen.all() and r.verdict == "完全冻结"


def test_freeze_verdict_unfrozen_when_driving_dominates():
    """屈服抗力 << 浮力驱动 ⇒ 未冻结（窄缝 + 大密度差 + 大井斜）。"""
    ny, nz = 4, 3
    tau = np.full((ny, nz), 1e-3)
    rho = np.zeros((ny, nz))
    rho[0, :] = 1000.0
    rho[-1, :] = 2000.0
    r = compute_stop_pump_freeze(
        tau_y_field=tau, rho_field=rho, gap_m=np.full(nz, 0.05),
        standoff=np.full(nz, 0.5), beta_deg=90.0)
    assert (~r.frozen).all() and r.verdict == "未冻结"


def test_freeze_zero_driving_means_frozen():
    """水平段（或方位向无密度差）⇒ 驱动为 0 ⇒ 判为静止稳定（冻结、ratio=inf）。"""
    ny, nz = 3, 2
    r = compute_stop_pump_freeze(
        tau_y_field=np.zeros((ny, nz)), rho_field=np.full((ny, nz), 1500.0),
        gap_m=np.full(nz, 0.02), standoff=0.7, beta_deg=0.0)
    assert (r.driving_pa == 0.0).all()
    assert r.frozen.all() and np.isinf(r.ratio).all()


def test_freeze_eccentricity_raises_bar():
    """同 τy/驱动下，偏心度越大（standoff 越小）⇒ 抗力侧 /(1+e) 越低（阈值更难过）。"""
    ny, nz = 3, 1
    tau = np.full((ny, nz), 30.0)
    rho = np.zeros((ny, nz)); rho[0, 0] = 1400.0; rho[-1, 0] = 1600.0
    lo_e = compute_stop_pump_freeze(tau_y_field=tau, rho_field=rho,
                                    gap_m=np.full(nz, 0.03), standoff=0.99,
                                    beta_deg=45.0)
    hi_e = compute_stop_pump_freeze(tau_y_field=tau, rho_field=rho,
                                    gap_m=np.full(nz, 0.03), standoff=0.10,
                                    beta_deg=45.0)
    assert hi_e.resisting_pa[0] < lo_e.resisting_pa[0]


def test_freeze_fields_from_state_mixes_by_volume_fraction():
    """逐相份额加权：纯水泥相 ⇒ 密度/τy 等于该相值。"""
    ny, nz = 2, 2
    comp = {"tail": np.ones((ny, nz))}
    r = fields_from_state(
        composition=comp,
        phase_rho_kg_m3={"tail": 1900.0, "mud": 1200.0},
        phase_tau_y_pa={"tail": 40.0, "mud": 5.0},
        gap_m=np.full(nz, 0.02), standoff=0.9, beta_deg=20.0)
    assert np.allclose(r.tau_y_min_pa, 40.0)


def test_freeze_shape_mismatch_raises():
    with pytest.raises(ValueError):
        compute_stop_pump_freeze(
            tau_y_field=np.zeros((3, 4)), rho_field=np.zeros((3, 4)),
            gap_m=np.zeros(5), standoff=0.5, beta_deg=10.0)
