# -*- coding: utf-8 -*-
"""Phase 1.5（2026-10-06）契约测试：run_opts 装配面全键穿透 + 两门语义区分。

上游
----
- 《详细执行计划_真温压响应_2026-10-06》§3 Phase 1.5「激活管道骨架」
- 评审A §4 断言 **A1.2**（每键必须「到达 solver 属性」——防"添加了≠在运行"）
- 执行窗口对抗核查 **C-01**（run_opts 键名 ≠ solver 形参名 ⇒ `**` 展开 TypeError；
  仓库签名闸门对"函数返回的 dict"是盲区 ⇒ 本文件补闸）

覆盖
----
1. **全 8 键穿透**：每个 run_opts 键经 kwargs 映射**到达**对应 solver 的构造属性
2. **casing 侧禁出 `enable_stream_yield_gate`**：该消费端开关只在 `AnnulusD2DGASolver` 上
   存在；casing kwargs 同样被 `**` 展开 ⇒ 出键即运行期 TypeError（C-01 同型陷阱的**不对称**版）
3. **默认不变量**：`normalize_run_opts(None)` 与 `_opts("off", on=False)` 映射出的 dict
   **逐键等于**扩展前的四键口径（关2：无新键 ⇒ 行为逐位不变）
4. **两门语义区分**：`enable_yield_gate`（生产端=wall 算不算）与
   `enable_stream_yield_gate`（消费端=wall 进不进算子）互不代偿
5. `build_pressure_field` 的 off/hydrostatic 两档 + 非法档响亮报错
"""
from __future__ import annotations

import inspect
import sys
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[2]
for _p in (str(_ROOT), str(_ROOT / "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from cemdisp.data.pressure_field import HydrostaticPressureField  # noqa: E402
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver  # noqa: E402
from cemdisp.transport1d.casing_flow import CasingFlowSolver  # noqa: E402
from entrypoints.run_sensitivity_current_20260916 import (  # noqa: E402
    RUN_OPTS_KEYS,
    annulus_kwargs_from_opts,
    build_pressure_field,
    casing_kwargs_from_opts,
    normalize_run_opts,
)
from entrypoints.run_sensitivity_temperature_t2_20261001 import _opts  # noqa: E402

_ANNULUS_PARAMS = set(inspect.signature(AnnulusD2DGASolver.__init__).parameters)
_CASING_PARAMS = set(inspect.signature(CasingFlowSolver.__init__).parameters)

# 扩展前的四键口径（关2 逐键对照基准）
_BASE_KEYS_BEFORE = ("enable_temperature_rheology", "rheology_formula_params",
                     "mud_extrapolate", "pressure_caliber")

_ALL_ON = {
    "enable_temperature_rheology": True,
    "temperature_mode": "static",
    "enable_yield_gate": True,
    "rheology_formula": {"scale": {"ca_tau0_q": 1.1}},
    "mud_extrapolate": True,
    "pressure_mode": "hydrostatic",
    "pressure_caliber": "mean",
    "enable_stream_yield_gate": True,
}


# --------------------------------------------------------------------------- #
# 1. 全键穿透（A1.2）
# --------------------------------------------------------------------------- #

def test_run_opts_keys_equal_normalize_default_keys():
    """`RUN_OPTS_KEYS` 必须与 `normalize_run_opts` 默认键集一致（防两处漂移）。"""
    assert tuple(normalize_run_opts(None)) == tuple(RUN_OPTS_KEYS)


def test_every_key_reaches_solver_attribute():
    """A1.2 全键穿透：八键全开 ⇒ 每个键都**到达** solver 构造属性（非"添加了没运行"）。"""
    opts = normalize_run_opts(_ALL_ON)
    ak = annulus_kwargs_from_opts(opts)
    ck = casing_kwargs_from_opts(opts)
    assert set(ak) <= _ANNULUS_PARAMS, f"坏键：{set(ak) - _ANNULUS_PARAMS}"
    assert set(ck) <= _CASING_PARAMS, f"坏键：{set(ck) - _CASING_PARAMS}"

    s = AnnulusD2DGASolver(dt=2.0, nz=4, ny=4, total_t=1.0, **ak)
    assert s.enable_temperature_rheology is True
    assert s.enable_yield_gate is True                    # 生产端
    assert s.enable_stream_yield_gate is True             # 消费端（Phase 1.5 新穿透位）
    assert s.mud_extrapolate is True
    assert s.pressure_caliber == "mean"
    assert s.rheology_formula_params is not None
    assert s.rheology_formula_params.ca_tau0_q[0] == pytest.approx(0.002793 * 1.1)

    c = CasingFlowSolver(**ck)
    assert c.enable_temperature_rheology is True
    assert c.mud_extrapolate is True
    assert c.pressure_caliber == "mean"
    assert c.rheology_formula_params is not None


def test_casing_kwargs_never_emit_stream_yield_gate():
    """C-01 不对称版：`enable_stream_yield_gate` **只在 annulus 上存在**。

    casing kwargs 同样被 `**` 展开进 `CasingFlowSolver(...)` ⇒ 出键即运行期 TypeError。
    本测试钉住"为什么必须分函数映射"。

    故意**不**写真调 `CasingFlowSolver(enable_stream_yield_gate=...)` 的反例：
    `tests/contract/test_entrypoint_signatures.py` 的「豁免集冻结为**唯一**一处蓄意反例」
    （`test_no_invented_dispersion.py:57`）是刻意收窄的守卫，不应为一条断言扩容；
    这里用 `inspect` 直接证明 `CasingFlowSolver.__init__` 没有该形参，效力等价。
    """
    opts = normalize_run_opts(_ALL_ON)
    assert "enable_stream_yield_gate" in annulus_kwargs_from_opts(opts)
    assert "enable_stream_yield_gate" not in casing_kwargs_from_opts(opts)
    assert "enable_stream_yield_gate" not in _CASING_PARAMS
    assert "enable_stream_yield_gate" in _ANNULUS_PARAMS


# --------------------------------------------------------------------------- #
# 2. 默认不变量（关2）：无新键 ⇒ 映射出的 dict 逐键等于扩展前
# --------------------------------------------------------------------------- #

def test_default_and_t_off_opts_map_to_legacy_four_key_caliber():
    """关2：`normalize_run_opts(None)` 与 `_opts("off", on=False)` 映射出的 dict
    **键集**等于 Phase 1 前的四键口径（新键默认值不产生额外 kwarg）。"""
    for opts in (normalize_run_opts(None),
                 normalize_run_opts(_opts("off", on=False)),
                 normalize_run_opts(_opts("static"))):
        ck = casing_kwargs_from_opts(opts)
        assert tuple(sorted(ck)) == tuple(sorted(_BASE_KEYS_BEFORE)), sorted(ck)
        assert ck["rheology_formula_params"] is None
        assert ck["mud_extrapolate"] is False
        assert ck["pressure_caliber"] == "shoe"
        ak = annulus_kwargs_from_opts(opts)
        assert "enable_stream_yield_gate" not in ak
        assert "enable_yield_gate" not in ak, "None ⇒ 不给键（沿用调用方口径）"


def test_no_key_is_bitwise_zero_trace_and_key_true_is_live():
    """关2「无此键时行为逐位不变」+ Phase 1.5 穿透证明（两者由同一算例给出）。

    ⚠️ **语义要点**：`enable_stream_yield_gate` 是**路径无关**的物理开关——
    它把 wall 场接进流函数算子，**在 T-off 下同样改变流场**（不是 T-on 专属）。
    因此正确判据是：
      - **不给键** / 给键=False（= 构造默认）⇒ 与基线**逐位相同**（关2 红线）；
      - **给键=True** ⇒ 结果**必须变**——这正是"确实到达算子"的证据
        （反制评审 A1.2 的"添加了 ≠ 在运行"）。
    同例的 `mud_extrapolate=True` 在 T-off 下为零痕迹（T-off 不调用 `fluid_at`）。
    """
    from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
    from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep
    from cemdisp.data.well_spec import DepthValuePoint, WellSpec
    from cemdisp.models2d.boundary_bridge import AnnulusInletState

    well = WellSpec(
        well_name="p15_default", top_md_m=100.0, bottom_md_m=400.0, shoe_md_m=400.0,
        hole_diameter_profile=(DepthValuePoint(100.0, 260.0), DepthValuePoint(400.0, 250.0)),
        liner_od_profile=(DepthValuePoint(100.0, 168.3), DepthValuePoint(400.0, 168.3)),
        inclination_profile=(DepthValuePoint(100.0, 2.0), DepthValuePoint(400.0, 5.0)),
        standoff_profile=(DepthValuePoint(100.0, 0.75), DepthValuePoint(400.0, 0.60)),
    )
    fluids = (
        FluidSpec("mud", FluidRole.MUD, 1050.0, RheologyModel.BINGHAM,
                  plastic_viscosity_pa_s=0.022, yield_stress_pa=6.0),
        FluidSpec("lead", FluidRole.LEAD, 1890.0, RheologyModel.POWER_LAW,
                  power_law_n=0.75, consistency_k=0.55),
        FluidSpec("tail", FluidRole.TAIL, 1920.0, RheologyModel.POWER_LAW,
                  power_law_n=0.70, consistency_k=0.90),
    )
    sched = PumpingSchedule(steps=tuple(
        PumpingScheduleStep(step_name=f"注{f.name}", fluid_name=f.name,
                            volume_m3=10.0, rate_m3_min=1.0) for f in fluids))

    def provider(t: float) -> AnnulusInletState:
        if t < 60.0:
            return AnnulusInletState(t, 0.01, "lead", (("lead", 1.0),))
        return AnnulusInletState(t, 0.01, "tail", (("tail", 1.0),))

    base = dict(dt=2.0, nz=24, ny=12, total_t=200.0,
                enable_cfl_adaptive=False, open_outlet=True)

    def run(**extra):
        return AnnulusD2DGASolver(**dict(base, **extra)).run(
            well, fluids, provider, schedule=sched)

    ref = run()                                              # 不给任何新键
    same = run(mud_extrapolate=True, enable_stream_yield_gate=False)
    for k in ("effective_efficiency", "eta_narrow"):
        assert ref.summary[k] == same.summary[k], f"{k}：无键/默认键必须逐位不变（关2）"

    gated = run(enable_stream_yield_gate=True)
    assert gated.summary["effective_efficiency"] != ref.summary["effective_efficiency"], \
        "给键 True 必须改变流场 —— 否则该键是死代码（未到达算子）"


# --------------------------------------------------------------------------- #
# 3. 两门语义区分
# --------------------------------------------------------------------------- #

def test_two_gates_are_distinct_and_not_substitutable():
    """`enable_yield_gate`（生产端）与 `enable_stream_yield_gate`（消费端）不是同一个键。"""
    o_prod = normalize_run_opts(_opts("static", yield_gate=False))
    o_cons = normalize_run_opts(_opts("static", stream_yield_gate=False))
    assert annulus_kwargs_from_opts(o_prod).get("enable_yield_gate") is False
    assert "enable_stream_yield_gate" not in annulus_kwargs_from_opts(o_prod)
    assert annulus_kwargs_from_opts(o_cons).get("enable_stream_yield_gate") is False
    assert "enable_yield_gate" not in annulus_kwargs_from_opts(o_cons)


def test_stream_yield_gate_type_validation():
    """类型闸：非 bool|None 响亮报错（防 'yes' 这类静默真值）。"""
    with pytest.raises(TypeError):
        normalize_run_opts({"enable_stream_yield_gate": "yes"})


# --------------------------------------------------------------------------- #
# 4. build_pressure_field 两档
# --------------------------------------------------------------------------- #

def test_build_pressure_field_modes():
    """off ⇒ None；hydrostatic ⇒ `HydrostaticPressureField`；非法档 ⇒ ValueError。"""
    from cemdisp.data.loaders.ht1_003_loader import load_ht1_003_tailpipe
    well, fluids, sched, _ = load_ht1_003_tailpipe()
    assert build_pressure_field(well, fluids, sched, "off") is None
    fld = build_pressure_field(well, fluids, sched, "hydrostatic")
    assert isinstance(fld, HydrostaticPressureField)
    assert np.isfinite(fld.P(well.shoe_md_m, 0.0))
    with pytest.raises(ValueError):
        build_pressure_field(well, fluids, sched, "nope")
