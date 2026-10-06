# -*- coding: utf-8 -*-
"""Phase 2（2026-10-06）契约测试：R2 μp/τy 真拆分（`include_yield_term`）。

上游
----
- 《详细执行计划_真温压响应_2026-10-06》§3 Phase 2「R2 屈服门接回 + μp/τy 真拆分」
- 评审A §4 断言 **D-02**（拆分式必带 η 上限）/ **A2.1**（η ∈ [1e-5, 3.0]）
- 执行窗口对抗核查 44 条（Phase 2 spec §10）：**L3-1/2**（算子不免疫——推翻"流函数免疫"）、
  **L6-1**（`b_field` 非逐位不变）、**F-4**（场口径对齐的域限定）、**W-1**（T2 `_opts` 陷阱）

覆盖（spec §5 的 T1/T2/T3/T6）
----
* **T1 默认路径逐位**：`include_yield_term=False` 与 HEAD 公式**位级**一致（关 2 红线）
* **T2 上限与场口径对齐**：flag=True ⇒ 输出恒 ∈ [1e-5, 3.0]；γ̇→0 收敛到 3.0；
  与 `AnnulusD2DGA._apparent_viscosity` 在**已声明域**内逐位同值，域外显式记录已知差异
* **T3 手算对照**：Bingham `PV + τy/γ̇`、HB `τy/γ̇ + K·γ̇^(n−1)`
* **T6 相消与"不免疫"回归锚**：4 站点同传时 ①剪切率输入不变 ②`λ_op·(1/F²)` 乘积不变量
  ③速度场**确实改变**（把 L3-1 的"不免疫"钉成锚，防回归到错误前提）

⚠ 测试红线（spec §5，对抗核查 A-6）：本文件**不得**写求解器构造的
``pytest.raises(TypeError)`` 坏键反例——那会扩容
``test_entrypoint_signatures.py::test_exemption_set_is_frozen_to_the_one_deliberate_negative``
的冻结豁免集。需要证明"casing 无此形参"时用 ``inspect``（见 test_run_opts_wiring.py）。
"""
from __future__ import annotations

import inspect
import struct

import numpy as np
import pytest

from cemdisp.data.fluid_spec import FluidSpec, FluidRole, RheologyModel
from cemdisp.data.well_spec import DepthValuePoint, EvaluationWindow, WellSpec
from cemdisp.models2d import annulus_d2dga as annulus_mod
from cemdisp.models2d import buoyancy
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import AnnulusInletState

CLIP = buoyancy._YIELD_TERM_MU_CLIP_PA_S


def _bits(x: float) -> bytes:
    return struct.pack("<d", float(x))


# --------------------------------------------------------------------------- #
# 测试流体族（覆盖四类本构 + 两个退化态）
# --------------------------------------------------------------------------- #
def _newton() -> FluidSpec:
    return FluidSpec("牛顿", FluidRole.MUD, 1000.0, RheologyModel.NEWTONIAN, 0.05, None)


def _bingham() -> FluidSpec:
    return FluidSpec("宾汉泥浆", FluidRole.MUD, 1960.0, RheologyModel.BINGHAM, 0.058, 9.2)


def _bingham_zero_yield() -> FluidSpec:
    """τy=0 的 Bingham（FluidSpec 明确允许；仓内 tests/contract/test_axial_dispersion.py:96 有此例）。"""
    return FluidSpec("宾汉零屈服", FluidRole.MUD, 1960.0, RheologyModel.BINGHAM, 5.0, 0.0)


def _power_law() -> FluidSpec:
    return FluidSpec("幂律", FluidRole.LEAD, 1900.0, RheologyModel.POWER_LAW,
                     None, None, power_law_n=0.719, consistency_k=0.815)


def _hb() -> FluidSpec:
    return FluidSpec("HB", FluidRole.SPACER, 2000.0, RheologyModel.HERSCHEL_BULKLEY,
                     None, 13.0, power_law_n=0.597, consistency_k=1.622)


def _hb_with_pv() -> FluidSpec:
    return FluidSpec("HB带PV", FluidRole.SPACER, 2000.0, RheologyModel.HERSCHEL_BULKLEY,
                     0.03, 13.0, power_law_n=0.597, consistency_k=1.622)


ALL_FLUIDS = [_newton(), _bingham(), _bingham_zero_yield(), _power_law(), _hb(), _hb_with_pv()]

# 覆盖：亚地板、地板附近、六速实测最低转速(5.11)、代表剪切率区(63–185)、高剪切
GAMMAS = [0.0, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-3, 0.1, 1.0, 3.0,
          5.11, 10.22, 63.0, 79.0, 123.0, 185.0, 1021.8, 1e4, 1e6]


# --------------------------------------------------------------------------- #
# T1 默认路径逐位（关 2 红线）
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("fluid", ALL_FLUIDS, ids=lambda f: f.name)
def test_t1_default_path_is_bitwise_legacy(fluid):
    """`include_yield_term=False`（含不传该 kwarg）必须与 HEAD 公式**位级**一致。"""
    for gv in GAMMAS:
        g = max(float(gv), 1e-8)
        if fluid.rheology_model in (RheologyModel.POWER_LAW, RheologyModel.HERSCHEL_BULKLEY) \
                and fluid.consistency_k is not None and fluid.power_law_n is not None:
            ref = float(fluid.consistency_k) * g ** (float(fluid.power_law_n) - 1.0)
            if fluid.plastic_viscosity_pa_s:
                ref = max(ref, float(fluid.plastic_viscosity_pa_s))
        else:
            assert fluid.plastic_viscosity_pa_s is not None
            ref = float(fluid.plastic_viscosity_pa_s)
        a = buoyancy.fluid_apparent_viscosity(fluid, gv)
        b = buoyancy.fluid_apparent_viscosity(fluid, gv, include_yield_term=False)
        assert _bits(a) == _bits(b) == _bits(ref), (
            f"{fluid.name} @γ̇={gv}: got {a!r}/{b!r}, ref {ref!r}")


# --------------------------------------------------------------------------- #
# T2 上限与场口径对齐（含域限定与负例）
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("fluid", ALL_FLUIDS, ids=lambda f: f.name)
def test_t2_output_is_clipped_into_unit_interval(fluid):
    """A2.1：flag=True ⇒ 输出恒 ∈ [1e-5, 3.0]（含 τy=0 与牛顿——**无条件** clip）。"""
    for gv in GAMMAS:
        v = buoyancy.fluid_apparent_viscosity(fluid, gv, include_yield_term=True)
        assert 1.0e-5 <= v <= CLIP, f"{fluid.name} @γ̇={gv} 越界：{v!r}"


def test_t2_zero_shear_saturates_at_cap_not_infinity():
    """D-02 核心：γ̇→0 时 τy/γ̇ 发散必须被上限吃掉（HEAD 无 clip 会给 ~1e9 Pa·s）。"""
    f = _bingham()
    for gv in (0.0, 1e-12, 1e-9, 1e-8):
        assert buoyancy.fluid_apparent_viscosity(f, gv, include_yield_term=True) == CLIP
    assert buoyancy.fluid_apparent_viscosity(f, 1e-9) == float(f.plastic_viscosity_pa_s)


def test_t2_zero_yield_bingham_still_clipped():
    """对抗核查 F-4 负例：τy=0 的 Bingham，若把 clip 挂在 τy 判据下会漏 clip（PV=5 > 3）。"""
    f = _bingham_zero_yield()
    assert f.yield_stress_pa == 0.0
    for gv in (1e-6, 1.0, 79.0, 1e6):
        assert buoyancy.fluid_apparent_viscosity(f, gv, include_yield_term=True) == CLIP


@pytest.mark.parametrize("fluid,gamma", [
    (_bingham(), 79.0), (_bingham(), 5.11), (_hb(), 79.0),
    (_hb_with_pv(), 79.0), (_power_law(), 79.0), (_newton(), 79.0),
])
def test_t2_matches_field_caliber_within_declared_domain(fluid, gamma):
    """与场口径 `AnnulusD2DGA._apparent_viscosity` **逐位同值**——**域限定**（spec §2.1）：
    ① γ̇ ≥ 1e-6（本处 floor 1e-8 vs 场口径 1e-6）；② HB 需 `plastic_viscosity_pa_s is None`。"""
    s = AnnulusD2DGASolver(ny=4, nz=2)
    field_v = float(s._apparent_viscosity(fluid, np.array([gamma]))[0])
    scalar_v = buoyancy.fluid_apparent_viscosity(fluid, gamma, include_yield_term=True)
    assert _bits(field_v) == _bits(scalar_v), (
        f"{fluid.name} @γ̇={gamma}: 场口径 {field_v!r} vs 标量口径 {scalar_v!r}")


# --------------------------------------------------------------------------- #
# T3 手算对照（本构式）
# --------------------------------------------------------------------------- #
def test_t3_bingham_manual():
    f = _bingham()
    for gv in (5.11, 63.0, 79.0, 123.0, 185.0):
        expect = float(f.plastic_viscosity_pa_s) + float(f.yield_stress_pa) / gv
        got = buoyancy.fluid_apparent_viscosity(f, gv, include_yield_term=True)
        assert got == pytest.approx(min(max(expect, 1e-5), CLIP), rel=0, abs=0)


def test_t3_hb_manual():
    f = _hb()
    for gv in (5.11, 63.0, 79.0, 123.0, 185.0):
        expect = (float(f.yield_stress_pa) / gv
                  + float(f.consistency_k) * gv ** (float(f.power_law_n) - 1.0))
        got = buoyancy.fluid_apparent_viscosity(f, gv, include_yield_term=True)
        assert got == pytest.approx(min(max(expect, 1e-5), CLIP), rel=0, abs=0)


def test_t3_split_multiplier_matches_reported_magnitude():
    """spec §1.2 量级自证：钻井液式 @T clamp 80°C 的倍率 = 1 + τy/(μp·γ̇)；
    γ̇∈[63,185] ⇒ 倍率 1.9–4.8（对抗核查 L3-5：按流体现算，不用通用区间）。"""
    f = FluidSpec("钻井液式", FluidRole.MUD, 1960.0, RheologyModel.BINGHAM, 0.046797, 10.9835)
    for gv, lo, hi in ((63.0, 4.0, 4.8), (79.0, 3.4, 4.0), (123.0, 2.6, 3.2), (185.0, 2.1, 2.6)):
        r = (buoyancy.fluid_apparent_viscosity(f, gv, include_yield_term=True)
             / buoyancy.fluid_apparent_viscosity(f, gv))
        assert lo <= r <= hi, f"γ̇={gv}: 倍率 {r:.3f} 超出 [{lo},{hi}]"


# --------------------------------------------------------------------------- #
# T6 相消不变量 + 「算子不免疫」回归锚
# --------------------------------------------------------------------------- #
_TOY_HOLE_MM, _TOY_OD_MM = 215.9, 139.7


def _toy_well() -> WellSpec:
    pts = lambda d, v: DepthValuePoint(depth_md_m=d, value=v)  # noqa: E731
    return WellSpec(
        well_name="toy", top_md_m=1000.0, bottom_md_m=1100.0,
        shoe_md_m=1100.0, hanger_md_m=1000.0,
        casing_id_mm=200.0, liner_od_mm=_TOY_OD_MM, liner_id_mm=108.0,
        hole_diameter_profile=[pts(1000.0, _TOY_HOLE_MM), pts(1100.0, _TOY_HOLE_MM)],
        inclination_profile=[pts(1000.0, 5.0), pts(1100.0, 5.0)],
        standoff_profile=[pts(1000.0, 0.60), pts(1100.0, 0.60)],
        evaluation_windows=[EvaluationWindow(name="w", top_md_m=1000.0,
                                             bottom_md_m=1100.0, window_type="full")],
    )


def _toy_fluids():
    mud = FluidSpec("toy泥浆", FluidRole.MUD, 1960.0, RheologyModel.BINGHAM, 0.058, 9.2)
    tail = FluidSpec("toy尾浆", FluidRole.TAIL, 1900.0, RheologyModel.HERSCHEL_BULKLEY,
                     None, 14.0, power_law_n=0.869, consistency_k=0.669)
    return mud, tail


def _inlet_for(q: float = 0.02):
    def _inlet(t: float) -> AnnulusInletState:
        return AnnulusInletState(time_s=t, flow_rate_m3_s=q, stage_name="pump",
                                 phase_fractions=(("cement", 1.0), ("tail", 1.0)))
    return _inlet


def _run_toy(include_yield_term: bool):
    """跑玩具算例；返回 (调用记录, 捕获的 b_field 列表, 结果)。"""
    calls: list[tuple[str, str, float, bool]] = []
    b_fields: list[np.ndarray] = []
    real_visc = buoyancy.fluid_apparent_viscosity
    real_solve = annulus_mod.solve_stream_function

    def _spy_visc(fluid, shear_rate, *, include_yield_term=False):
        calls.append((inspect.currentframe().f_back.f_code.co_name,
                      fluid.name, float(shear_rate), bool(include_yield_term)))
        return real_visc(fluid, shear_rate, include_yield_term=include_yield_term)

    def _spy_solve(geom, c_bar, eta1, eta2, m, b_field, **kw):
        b_fields.append(np.array(b_field, dtype=float, copy=True))
        return real_solve(geom, c_bar, eta1, eta2, m, b_field, **kw)

    mud, tail = _toy_fluids()
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(buoyancy, "fluid_apparent_viscosity", _spy_visc)
        mp.setattr(annulus_mod, "solve_stream_function", _spy_solve)
        solver = AnnulusD2DGASolver(dt=1.0, nz=20, ny=8, total_t=20.0,
                                    enable_cfl_adaptive=False,
                                    include_yield_term=include_yield_term)
        res = solver.run(_toy_well(), (mud, tail), _inlet_for())
    return calls, b_fields, res


def test_t6_wiring_and_shear_rate_inputs_are_identical():
    """4 站点同传 + **剪切率输入逐位不变**（这正是 λ_op 与 F² 共用同一个 μ̂₁ 的前提）。

    固定 dt（`enable_cfl_adaptive=False`，dt=1.0 / total_t=20）⇒ 两次 run 步数相同、
    逐步状态可比（若用 CFL 自适应，拆分改变流场 ⇒ dt 分岔 ⇒ 步数不同，无法对位比较）。
    """
    from collections import Counter
    calls_off, b_off, _ = _run_toy(False)
    calls_on, b_on, _ = _run_toy(True)

    names_off = {c[0] for c in calls_off}
    names_on = {c[0] for c in calls_on}
    assert "_froude_squared_at" in names_off and "_buoyancy_number_at" in names_off
    assert "_velocity_stream_function" in names_off
    assert names_on == names_off, f"站点集变化：off={sorted(names_off)} on={sorted(names_on)}"

    # 站点实例数：每步 `_velocity_stream_function` 内 η₁+η₂ 各一次；`_froude_squared_at`
    # 每步两次（流函数段一次 + run 循环 I3 块一次）。两者计数相等是结构不变量。
    cnt_off = Counter(c[0] for c in calls_off)
    assert cnt_off["_velocity_stream_function"] == 2 * 0 + cnt_off["_froude_squared_at"], cnt_off
    assert cnt_off["_velocity_stream_function"] % 2 == 0, cnt_off

    # ① 逐步：把**泥浆**的调用按剪切率（位级）切成连续段；每段内必须同时出现
    #    `_velocity_stream_function`（η₁ 站）与 `_froude_squared_at`（F² 站）
    #    —— 这是 λ_op 与 F² 共用同一个 μ̂₁ 的直接证据。
    steps: list[list[tuple]] = []
    core = calls_off
    for _i, _c in enumerate(calls_off):        # 截掉末尾 summary 段的 _buoyancy_number_at
        if _c[0] == "_buoyancy_number_at":
            core = calls_off[:_i]
            break
    # ① η₁ 站（VSF 内）与其**紧邻的下一次** F² 站必须读到同一位级剪切率
    #    —— λ_op（∝η₁）与 F²（∝μ̂₁）相消的前提，正是"两处读同一个 μ̂₁"。
    n_pair = 0
    for i, c in enumerate(core):
        if c[0] != "_velocity_stream_function" or c[1] != "toy泥浆":
            continue
        nxt = next((d for d in core[i + 1:] if d[0] == "_froude_squared_at"), None)
        assert nxt is not None, "VSF 后未见 F² 站"
        assert _bits(c[2]) == _bits(nxt[2]), (
            f"η₁ 站与 F² 站剪切率不一致：{c[2]!r} vs {nxt[2]!r} ⇒ 相消前提被破坏")
        n_pair += 1
    assert n_pair >= 2, f"配对样本不足：{n_pair}"

    # ② 跨 run：步数相同（固定 dt）；且**第一步**的 η₁ 站剪切率逐位相同
    #    （同初始态 ⇒ 拆分不得改变第一条剪切率读数）。更深步的剪切率是上一步流场的
    #    **输出**，拆分本就会改变它——故**不**断言全序列跨 run 相等（那是错误的判据）。
    assert len(b_off) == len(b_on) > 0
    mud_off = [c for c in calls_off if c[1] == "toy泥浆"]
    mud_on = [c for c in calls_on if c[1] == "toy泥浆"]
    assert _bits(mud_off[0][2]) == _bits(mud_on[0][2]), "第一步剪切率即被拆分改变 ⇒ 异常"

    # ③ 开关确实到达全部站点
    assert all(c[3] is False for c in calls_off)
    assert all(c[3] is True for c in calls_on)


def test_t6_bfield_amplitude_is_insensitive_but_not_bitwise():
    """`b_field` 幅度对拆分**近乎不敏感**（λ_op·(1/F²) 相消），但**非逐位不变**。

    对抗核查 L6-1 实测 maxrelΔ = 1.43%（全部来自 χ(m)）。阈值 5% 作**回归锚**
    （"只传部分站点"的负对照见下一条，那种情形变化 ~50%+）。"""
    _, b_off, _ = _run_toy(False)
    _, b_on, _ = _run_toy(True)
    assert len(b_off) == len(b_on) > 0
    # 逐步比对（固定 dt ⇒ 步数相同）；取全时程最大相对变化
    rel_max = 0.0
    for a, b in zip(b_off, b_on):
        denom = max(float(np.max(np.abs(a))), 1e-30)
        rel_max = max(rel_max, float(np.max(np.abs(b - a))) / denom)
    assert rel_max < 0.05, f"b_field 相对变化 {rel_max:.4g} 过大：λ_op·(1/F²) 相消可能被破坏"
    assert not np.array_equal(b_off[-1], b_on[-1]), (
        "b_field 逐位相同——若真实如此，L6-1 的 χ(m) 结论需复核（勿据此改断言）")


def test_t6_negative_control_partial_wiring_breaks_cancellation():
    """负对照（对抗核查 L3-3）：只把开关传给 F² 站而**不传** λ_op 侧 η₁ ⇒
    `b ∝ 1/μ̂₁` 直接掉 ~半 ⇒ 相对变化远大于 5%。证明上面的阈值确有区分力。"""
    real_visc = buoyancy.fluid_apparent_viscosity
    mud, tail = _toy_fluids()

    def _b_last(split_at_froude_only: bool) -> np.ndarray:
        b_fields: list[np.ndarray] = []
        real_solve = annulus_mod.solve_stream_function

        def _spy_visc(fluid, shear_rate, *, include_yield_term=False):
            on = (include_yield_term
                  and inspect.currentframe().f_back.f_code.co_name == "_froude_squared_at")
            return real_visc(fluid, shear_rate,
                             include_yield_term=bool(on and split_at_froude_only))

        def _spy_solve(geom, c_bar, eta1, eta2, m, b_field, **kw):
            b_fields.append(np.array(b_field, dtype=float, copy=True))
            return real_solve(geom, c_bar, eta1, eta2, m, b_field, **kw)

        with pytest.MonkeyPatch.context() as mp:
            mp.setattr(buoyancy, "fluid_apparent_viscosity", _spy_visc)
            mp.setattr(annulus_mod, "solve_stream_function", _spy_solve)
            solver = AnnulusD2DGASolver(dt=1.0, nz=20, ny=8, total_t=20.0,
                                        enable_cfl_adaptive=False,
                                        include_yield_term=True)
            solver.run(_toy_well(), (mud, tail), _inlet_for())
        return b_fields[-1]

    base = _b_last(False)          # 两站点对称（都不开）
    asym = _b_last(True)           # 仅 F² 站开 ⇒ 相消被破坏
    denom = max(float(np.max(np.abs(base))), 1e-30)
    rel = float(np.max(np.abs(asym - base))) / denom
    assert rel > 0.05, (
        f"非对称接线未造成 >5% 的 b_field 变化（实测 {rel:.4g}）——阈值无区分力")
