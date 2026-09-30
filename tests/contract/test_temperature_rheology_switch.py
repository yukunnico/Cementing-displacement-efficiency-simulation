# -*- coding: utf-8 -*-
"""Task 6（T1-1）契约测试：`annulus_d2dga` 接线 `enable_temperature_rheology`
温变流变总开关 + `fluid_at` 构造层接线（默认关 ⇒ 逐位=HEAD）。

覆盖（task-6-brief 需求逐条）：

1. **默认关**：`_SWITCH_DEFAULTS["enable_temperature_rheology"] is False`，且与
   `__init__` 签名默认值一致（单一真源，由 test_dead_switch_guard 泛化守护；此处
   再钉一个直接断言）。
2. **T-off 逐位（L1 硬约束）**：显式传 `enable_temperature_rheology=False` 与完全
   不传的求解器全结果**逐位一致**（浓度场/wall float64 sha256 + summary 逐项）。
3. **关=零痕迹**：关时不产出任何温变派生痕迹（summary 无审计键、无派生流体暴露）。
4. **开 + ConstantTemperatureField(60)**（接口级）：四相经 `fluid_at` 绝对替换派生，
   派生流体是 Bingham 且 `yield_stress_pa > 0`；原 FluidSpec 不被就地修改
   （`fluid_at` 走 `dataclasses.replace`，frozen）。
5. **审计接线（最小）**：开时 `fluid_at` 的 clamp/borrow/no_replace/p_default 计数
   进 summary（`temperature_rheology_audit`）+ `[D2DGA]` 日志行可见。
6. **开关真被消费**：开 vs 关的浓度场摘要不同（构造层接线确实进入动力学）。

逐步刷新 geom["T"] / props 合并 / 屈服门 / memo key 是 Task 7 的事，本文件不覆盖。
"""
from __future__ import annotations

import hashlib

import numpy as np

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.well_spec import DepthValuePoint, WellSpec
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver, _SWITCH_DEFAULTS
from cemdisp.models2d.boundary_bridge import AnnulusInletState

# 生产级排量（同 test_hb_closure_wiring：低流速下 Bingham 泥浆整层低于屈服门槛）。
_Q_M3S = 1.0 / 60.0


# --------------------------------------------------------------------------- #
# 公共装配（field-like 小域；与 test_hb_closure_wiring 同源，便于跨任务对照）
# --------------------------------------------------------------------------- #

def _well_spec() -> WellSpec:
    top, bottom = 100.0, 400.0
    return WellSpec(
        well_name="T11_temp_rheo_switch",
        top_md_m=top,
        bottom_md_m=bottom,
        shoe_md_m=bottom,
        hole_diameter_profile=(DepthValuePoint(top, 260.0), DepthValuePoint(bottom, 250.0)),
        liner_od_profile=(DepthValuePoint(top, 168.3), DepthValuePoint(bottom, 168.3)),
        inclination_profile=(DepthValuePoint(top, 2.0), DepthValuePoint(bottom, 5.0)),
        standoff_profile=(DepthValuePoint(top, 0.75), DepthValuePoint(bottom, 0.60)),
    )


def _fluids() -> tuple:
    """Bingham 泥浆/隔离液 + POWER_LAW 水泥（生产 8 井口径，spec 无 τy）。

    密度刻意落在 fluid_at 分派表的可替换域：泥浆走公式、隔离液 1.12 g/cm³ 触发
    borrow（域外就近端）、领浆 1.89 / 尾浆 1.92 落水泥档。
    """
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


def _run_solver(**overrides) -> tuple:
    """返回 (solver, result)——solver 用于读派生流体暴露面，result 用于摘要对照。"""
    params = dict(nz=24, ny=12, dt=2.0, total_t=240.0, enable_cfl_adaptive=False,
                  open_outlet=True)
    params.update(overrides)
    solver = AnnulusD2DGASolver(**params)
    return solver, solver.run(_well_spec(), _fluids(), _provider)


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
        "b": float(s["浮力数_b"]),
    }


# --------------------------------------------------------------------------- #
# 1. 默认关（单一真源）
# --------------------------------------------------------------------------- #

def test_switch_default_is_off():
    assert _SWITCH_DEFAULTS["enable_temperature_rheology"] is False
    assert AnnulusD2DGASolver().enable_temperature_rheology is False
    assert AnnulusD2DGASolver().temperature_rheology_t_c == 60.0


# --------------------------------------------------------------------------- #
# 2. T-off 逐位（L1 硬约束）
# --------------------------------------------------------------------------- #

def test_default_off_is_bitwise_head():
    """显式传默认值与不传的求解器全结果逐位一致。"""
    _, base = _run_solver()
    _, explicit = _run_solver(enable_temperature_rheology=False,
                              temperature_rheology_t_c=60.0)
    assert _result_digests(base) == _result_digests(explicit)


def test_off_leaves_no_trace():
    """关=零痕迹：无派生流体暴露面、summary 无审计键（输出风格与 HEAD 同型）。"""
    solver, res = _run_solver()
    assert not hasattr(solver, "_temp_rheo_fluids")
    assert "temperature_rheology_audit" not in res.summary


# --------------------------------------------------------------------------- #
# 3. 开 + ConstantTemperatureField(60)：派生流体是 Bingham 且 τy > 0（接口级）
# --------------------------------------------------------------------------- #

def test_on_derives_bingham_at_60c():
    solver, _ = _run_solver(enable_temperature_rheology=True,
                            temperature_rheology_t_c=60.0)
    derived = solver._temp_rheo_fluids
    assert set(derived) == {"mud", "lead", "tail", "spacer"}
    for name, fluid in derived.items():
        assert fluid.rheology_model is RheologyModel.BINGHAM, (
            f"{name} 派生后应为 Bingham（绝对替换），得到 {fluid.rheology_model}")
        assert fluid.yield_stress_pa is not None and fluid.yield_stress_pa > 0.0, (
            f"{name} 派生 τy 应 >0，得到 {fluid.yield_stress_pa}")
        assert fluid.plastic_viscosity_pa_s is not None and fluid.plastic_viscosity_pa_s > 0.0
        assert fluid.power_law_n is None and fluid.consistency_k is None


def test_on_does_not_mutate_original_specs():
    """fluid_at 走 dataclasses.replace（frozen）：原 FluidSpec 逐项不变。"""
    mud0, spacer0, lead0, tail0 = _fluids()
    _run_solver(enable_temperature_rheology=True, temperature_rheology_t_c=60.0)
    assert lead0.rheology_model is RheologyModel.POWER_LAW
    assert lead0.power_law_n == 0.75 and lead0.yield_stress_pa is None
    assert tail0.rheology_model is RheologyModel.POWER_LAW
    assert tail0.power_law_n == 0.70 and tail0.yield_stress_pa is None
    assert mud0.yield_stress_pa == 6.0 and mud0.plastic_viscosity_pa_s == 0.022
    assert spacer0.yield_stress_pa == 9.0


def test_on_is_consumed_into_dynamics():
    """开关真被消费：开 vs 关的浓度场摘要不同（构造层接线进入动力学）。"""
    _, off = _run_solver()
    _, on = _run_solver(enable_temperature_rheology=True,
                        temperature_rheology_t_c=60.0)
    d_off, d_on = _result_digests(off), _result_digests(on)
    assert d_off["cement"] != d_on["cement"], (
        "开/关浓度场摘要必须不同，否则构造层接线未进入求解（开关空转）")


# --------------------------------------------------------------------------- #
# 4. 审计接线（最小）：clamp/borrow/no_replace 计数进 summary + [D2DGA] 日志行
# --------------------------------------------------------------------------- #

def test_audit_counts_in_summary_and_log(capsys):
    solver, res = _run_solver(enable_temperature_rheology=True,
                              temperature_rheology_t_c=60.0)
    out = capsys.readouterr().out
    assert "温变流变" in out, f"[D2DGA] 日志行应含温变流变审计，得到：\n{out[-500:]}"
    audit = res.summary["temperature_rheology_audit"]
    assert isinstance(audit, dict) and audit
    # 本装配的可预期事件：隔离液 1.12 g/cm³ 域外借用端点公式 + P 缺省常压。
    assert audit.get("borrow", 0) >= 1, audit
    assert audit.get("p_default", 0) >= 1, audit
    # 暴露面与 summary 同源（同一 run 的计数）
    assert solver._temp_rheo_audit_counts == audit
