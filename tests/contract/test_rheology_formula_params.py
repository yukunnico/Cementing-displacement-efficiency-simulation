# -*- coding: utf-8 -*-
"""R1（2026-10-06）契约测试：温变流变公式系数参数化。

上游
----
- ``docs/superpowers/specs/2026-10-06-phase1-fluid-at-signature-wave-design.md`` §1
- 评审A《R1-R4 设计对抗评审》§4 断言 **A1.1 / A1.2 / A1.3** + 缺陷 **D-09**（负物性敞口）
- 执行窗口对抗核查（2026-10-06）MAJOR：**M4 符号陷阱**、**M5/M7 run_opts 键名映射**
  （仓库签名闸门 ``check_call_signatures`` 对"函数返回的 dict"是盲区，本文件补闸）

覆盖
----
1. 字段值 == 2026-10-01 原字面量（**struct.pack 位级**）——防搬运笔误/符号翻转
2. 默认参输出逐位 == 原表达式（断点两侧 + 全族温压网格）——防 Horner/重排（A1.1）
3. ``params=None`` ≡ ``RheologyFormulaParams()``（同一默认单例路径）
4. 扰动真消费：组 A τ₀ 系数 ×1.1 → 输出 == 手算
5. 物理可达域 τy≥0 ∧ μp>0；**已知负值角落**（2.05 式 T≳172°C 且低压）显式登记为
   HEAD 同型（本仓改造前即存在）——非本次引入
6. 审计：``params=None`` 时事件序列与"不传参"逐条相同（A1.3）
7. run_opts 装配：``RUN_OPTS_KEYS`` == normalize 默认键；kwargs 映射**改名正确**且
   返回键集 ⊆ 两 solver ``__init__`` 形参集（A1.2 + M5/M7）
8. solver 接线：ctor 新参**到达** ``_phase_props``（A1.2 的"添加了≠在运行"）
"""
from __future__ import annotations

import inspect
import math
import struct
import sys
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[2]
for _p in (str(_ROOT), str(_ROOT / "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel  # noqa: E402
from cemdisp.data.rheology_vs_temperature import (  # noqa: E402
    RheologyFormulaParams,
    fluid_at,
    get_audit,
    reset_audit,
)
from cemdisp.data.well_spec import DepthValuePoint, WellSpec  # noqa: E402
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver  # noqa: E402
from cemdisp.transport1d.casing_flow import CasingFlowSolver  # noqa: E402
from entrypoints.run_sensitivity_current_20260916 import (  # noqa: E402
    RUN_OPTS_KEYS,
    annulus_kwargs_from_opts,
    casing_kwargs_from_opts,
    normalize_run_opts,
    params_from_spec,
)

# 原字面量复核用（与本模块字段对照，独立重述）
_SP_1P95_TY_REF = (7.379374, -0.001857, 0.0225351,
                   -0.000157779, -0.00000615055, -0.0000993563)
_SP_1P95_MP_REF = (0.0912039, -0.000719380, 0.000167466,
                   0.00000277807, -0.00000322380, 0.00000119469)
_SP_2P05_TY_REF = (10.787391, -0.0206255, 0.0320167,
                   -0.000247286, 0.000219522, -0.000253669)
_SP_2P05_MP_REF = (0.125694, -0.00108616, 0.000323974,
                   0.00000341573, -0.00000316667, 0.000000714405)


def _bits(x: float) -> str:
    """IEEE754 位级指纹（`==` 对 -0.0/0.0 与 NaN 不足，位级才可信）。"""
    return struct.pack("<d", float(x)).hex()


@pytest.fixture(autouse=True)
def _clean_audit():
    reset_audit()
    yield
    reset_audit()


# --------------------------------------------------------------------------- #
# 0. 流体与几何工厂
# --------------------------------------------------------------------------- #

def _spacer(name: str = "隔离液1", rho_g_cm3: float = 1.95) -> FluidSpec:
    return FluidSpec(name=name, role=FluidRole.SPACER, density_kg_m3=rho_g_cm3 * 1000.0,
                     rheology_model=RheologyModel.POWER_LAW,
                     power_law_n=0.7, consistency_k=0.3)


def _cement(name: str = "领浆", rho_g_cm3: float = 2.10,
            role: FluidRole = FluidRole.LEAD) -> FluidSpec:
    return FluidSpec(name=name, role=role, density_kg_m3=rho_g_cm3 * 1000.0,
                     rheology_model=RheologyModel.POWER_LAW,
                     power_law_n=0.8, consistency_k=0.5)


def _mud(name: str = "钻井液", role: FluidRole = FluidRole.MUD) -> FluidSpec:
    return FluidSpec(name=name, role=role, density_kg_m3=1900.0,
                     rheology_model=RheologyModel.BINGHAM,
                     plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)


_ALL_FLUIDS = (
    _mud(), _mud("平衡液", FluidRole.WASH), _mud("先导浆", FluidRole.WASH),
    _spacer("隔离液1", 1.95), _spacer("隔离液1", 2.05), _spacer("隔离液1", 2.00),
    _spacer("隔离液2", 1.75),
    _cement("领浆", 2.10), _cement("领浆", 1.93), _cement("尾浆", 1.90, FluidRole.TAIL),
    _mud("替钻井液", FluidRole.DISPLACEMENT), _mud("冲洗液（FLUSHER）", FluidRole.FLUSHER),
)

_T_GRID = (20.0, 40.0, 60.0, 80.0, 99.999, 100.0, 100.001, 120.0, 150.0, 170.0)
_P_GRID = (None, 0.1, 70.0, 140.0, 200.0)


def _well_spec() -> WellSpec:
    top, bottom = 100.0, 400.0
    return WellSpec(
        well_name="R1_params_wiring",
        top_md_m=top, bottom_md_m=bottom, shoe_md_m=bottom,
        hole_diameter_profile=(DepthValuePoint(top, 260.0), DepthValuePoint(bottom, 250.0)),
        liner_od_profile=(DepthValuePoint(top, 168.3), DepthValuePoint(bottom, 168.3)),
        inclination_profile=(DepthValuePoint(top, 2.0), DepthValuePoint(bottom, 5.0)),
        standoff_profile=(DepthValuePoint(top, 0.75), DepthValuePoint(bottom, 0.60)),
    )


# --------------------------------------------------------------------------- #
# 1. 字段值 == 原字面量（位级）——防符号翻转 / 搬运笔误（评审 M4）
# --------------------------------------------------------------------------- #

def test_field_values_are_bit_identical_to_original_literals():
    """每个字段必须与 2026-10-01 原字面量**位级相同**。

    评审 M4「符号陷阱」：``a*T*T - b*T + c`` 的 ``-b`` 若被搬成正数再把字段值
    当 ``b`` 代入 ``+(-b)*T``，线性项符号静默翻转（实测 τy 偏差 5.7×）。
    本测试把符号约定钉死在**字段值本身**上。
    """
    p = RheologyFormulaParams()
    expected = {
        "break_T": 100.0,
        "ca_tau0_q": (0.002793, -0.391702, 21.6378),
        "ca_tau0_l": (0.624883, -51.7949),
        "ca_mup": (0.53129, -0.012495, -0.00204),
        "ca_T_lo": 20.0, "ca_T_hi": 170.0,
        "cb_tau0_q": (4.433e-4, -0.07843, 5.962),
        "cb_tau0_l": (0.03206, 1.979),
        "cb_mup": (1.355e-5, -6.570e-4, 0.08215),
        "cb_T_lo": 20.0, "cb_T_hi": 200.0,
        "sp_1p95_ty": _SP_1P95_TY_REF,
        "sp_1p95_mp": _SP_1P95_MP_REF,
        "sp_2p05_ty": _SP_2P05_TY_REF,
        "sp_2p05_mp": _SP_2P05_MP_REF,
        "sp_T_lo": 20.0, "sp_T_hi": 200.0, "sp_P_lo": 0.1, "sp_P_hi": 200.0,
        "sp_p_default": 0.1,
        "sp_lo": 1.90, "sp_hi": 2.10, "sp_mid": 2.000,
        "sp_anchor_lo": 1.95, "sp_anchor_hi": 2.05,
        "cm_lo": 1.88, "cm_hi": 2.12, "cm_mid": 2.000,
        "cm_anchor_a": 2.10, "cm_anchor_b": 1.90,
        "mud_mup": (0.297154, -0.02310525),
        "mud_tauy": (21.346710, -0.00830631),
        "mud_T_lo": 40.0, "mud_T_hi": 80.0,
    }
    assert set(expected) == set(RheologyFormulaParams.__dataclass_fields__), \
        "字段增删必须同步本断言（防参数化面悄悄扩大）"
    for name, want in expected.items():
        got = getattr(p, name)
        if isinstance(want, tuple):
            assert len(got) == len(want), name
            for g, w in zip(got, want):
                assert _bits(g) == _bits(w), f"{name}: {g!r} != {w!r}（位级）"
        else:
            assert _bits(got) == _bits(want), f"{name}: {got!r} != {want!r}（位级）"


# --------------------------------------------------------------------------- #
# 2. 默认参输出逐位 == 原表达式（A1.1；防 Horner / 重排）
# --------------------------------------------------------------------------- #

def _orig_tauy_a(T: float) -> float:
    """2026-10-01 原字面量表达式（对照组，独立重述）。"""
    if T <= 100.0:
        return 0.002793 * T * T - 0.391702 * T + 21.6378
    return 0.624883 * T - 51.7949


def _orig_mup_a(T: float) -> float:
    return 0.53129 * math.exp(-0.012495 * T) - 0.00204


def _orig_tauy_b(T: float) -> float:
    if T <= 100.0:
        return 4.433e-4 * T * T - 0.07843 * T + 5.962
    return 0.03206 * T + 1.979


def _orig_mup_b(T: float) -> float:
    return 1.355e-5 * T * T - 6.570e-4 * T + 0.08215


def _orig_mud_tauy(T: float) -> float:
    return 21.346710 * math.exp(-0.00830631 * T)


def _orig_mud_mup(T: float) -> float:
    return 0.297154 * math.exp(-0.02310525 * T)


def _orig_quad6(c, T, P):
    return c[0] + c[1] * T + c[2] * P + c[3] * T * T + c[4] * T * P + c[5] * P * P


def test_default_params_bitwise_equal_original_expressions():
    """默认参 ⇒ 各公式函数输出与原字面量表达式**位级**相同（含断点两侧）。"""
    from cemdisp.data import rheology_vs_temperature as R
    p = RheologyFormulaParams()
    for T in _T_GRID:
        assert _bits(R._cement_a_tau0(T, p)) == _bits(_orig_tauy_a(T)), T
        assert _bits(R._cement_a_mup(T, p)) == _bits(_orig_mup_a(T)), T
        assert _bits(R._cement_b_tau0(T, p)) == _bits(_orig_tauy_b(T)), T
        assert _bits(R._cement_b_mup(T, p)) == _bits(_orig_mup_b(T)), T
        assert _bits(R._mud_tauy(T, p)) == _bits(_orig_mud_tauy(T)), T
        assert _bits(R._mud_mup(T, p)) == _bits(_orig_mud_mup(T)), T
    for T in _T_GRID:
        for P in (0.1, 70.0, 140.0):
            assert _bits(R._quad6(p.sp_1p95_ty, T, P)) == \
                _bits(_orig_quad6(_SP_1P95_TY_REF, T, P)), (T, P)
            assert _bits(R._quad6(p.sp_2p05_ty, T, P)) == \
                _bits(_orig_quad6(_SP_2P05_TY_REF, T, P)), (T, P)


def test_default_params_fluid_at_bitwise_matches_reference_grid():
    """`fluid_at` 全族 × 温压网格：默认参输出 == 手工重述的旧口径输出（位级）。"""
    from cemdisp.data import rheology_vs_temperature as R
    for f in _ALL_FLUIDS:
        d = f.density_kg_m3 / 1000.0
        family, _note = R._route(f)
        for T in _T_GRID:
            for P in _P_GRID:
                if family == "no_replace":
                    assert fluid_at(f, T, P) is f
                    continue
                got = fluid_at(f, T, P)
                if family == "mud":
                    Tm = min(max(T, 40.0), 80.0)
                    exp_ty, exp_mp = _orig_mud_tauy(Tm), _orig_mud_mup(Tm)
                elif family == "spacer":
                    Tc = min(max(T, 20.0), 200.0)
                    Pc = min(max(0.1 if P is None else P, 0.1), 200.0)
                    c = (_SP_1P95_TY_REF, _SP_1P95_MP_REF) if d < 2.000 \
                        else (_SP_2P05_TY_REF, _SP_2P05_MP_REF)
                    exp_ty, exp_mp = _orig_quad6(c[0], Tc, Pc), _orig_quad6(c[1], Tc, Pc)
                else:  # cement
                    if not (1.88 <= d <= 2.12):
                        assert got is f
                        continue
                    if d < 2.000:
                        Tc = min(max(T, 20.0), 200.0)
                        exp_ty, exp_mp = _orig_tauy_b(Tc), _orig_mup_b(Tc)
                    else:
                        Tc = min(max(T, 20.0), 170.0)
                        exp_ty, exp_mp = _orig_tauy_a(Tc), _orig_mup_a(Tc)
                if exp_ty < 0.0:
                    # 已知负值角落（D-09 登记）：两侧同抛，见 test_known_negative_corner
                    continue
                assert _bits(got.yield_stress_pa) == _bits(exp_ty), (f.name, T, P)
                assert _bits(got.plastic_viscosity_pa_s) == _bits(exp_mp), (f.name, T, P)


def test_params_none_is_default_singleton_path():
    """``params=None`` 与显式默认实例逐位同。"""
    for f in _ALL_FLUIDS:
        for T in _T_GRID:
            for P in (None, 0.1, 70.0):
                a = fluid_at(f, T, P)
                b = fluid_at(f, T, P, params=RheologyFormulaParams())
                if a is f:
                    assert b is f
                else:
                    assert _bits(a.yield_stress_pa) == _bits(b.yield_stress_pa)
                    assert _bits(a.plastic_viscosity_pa_s) == _bits(b.plastic_viscosity_pa_s)


# --------------------------------------------------------------------------- #
# 3. 扰动真消费（L2 流变响应：公式系数是敏感性轴）
# --------------------------------------------------------------------------- #

def test_perturbation_is_consumed_and_matches_hand_calc():
    """组 A τ₀ 系数 ×1.1 ⇒ `fluid_at` 输出 == 手算（证明参**真被消费**）。"""
    p = params_from_spec({"scale": {"ca_tau0_q": 1.1}})
    assert p is not None
    assert p.ca_tau0_q == (0.002793 * 1.1, -0.391702 * 1.1, 21.6378 * 1.1)
    assert p.ca_tau0_l == RheologyFormulaParams().ca_tau0_l, "未列字段不得被改动"
    T = 60.0
    out = fluid_at(_cement("领浆", 2.10), T, params=p)
    q0, q1, q2 = p.ca_tau0_q
    assert out.yield_stress_pa == q0 * T * T + q1 * T + q2
    base = fluid_at(_cement("领浆", 2.10), T)
    assert out.yield_stress_pa != base.yield_stress_pa


def test_params_from_spec_rejects_bad_payload():
    """载荷错字段 / 空载荷 / 非 dict 一律响亮报错（防"配错了静默按默认跑"）。"""
    with pytest.raises(KeyError):
        params_from_spec({"scale": {"not_a_field": 1.1}})
    with pytest.raises(KeyError):
        params_from_spec({"unknown_sub": {}})
    with pytest.raises(ValueError):
        params_from_spec({})
    with pytest.raises(TypeError):
        params_from_spec([1, 2])


# --------------------------------------------------------------------------- #
# 4. 正物性域（A1.1）+ 已知负值角落登记（D-09）
# --------------------------------------------------------------------------- #

def test_positive_properties_over_physical_domain():
    """物理可达域（T ≤ 170°C，三重点井 T≤155°C）内 τy≥0 ∧ μp>0，全族全压。"""
    for f in _ALL_FLUIDS:
        for T in (20.0, 60.0, 100.0, 100.001, 120.0, 150.0, 170.0):
            for P in (None, 0.1, 70.0, 140.0, 200.0):
                out = fluid_at(f, T, P, mud_extrapolate=True)
                if out is f:
                    continue
                assert out.yield_stress_pa >= 0.0, (f.name, T, P, out.yield_stress_pa)
                assert out.plastic_viscosity_pa_s > 0.0, (f.name, T, P)


def test_known_negative_corner_is_pre_existing_and_unchanged():
    """**已知负值角落**（D-09）：2.05 式隔离液 T≳172°C 且低压下 τy<0。

    该角落**在本仓改造前（HEAD）即已存在**（同一 `_quad6`、同一系数），R1 不引入、
    不修复——本测试把它**显式登记**下来，防后续有人误判为 R1 回归；
    也防"全域 τy≥0"被当成无条件的验收判据（A1.1 须限定在物理可达域）。
    """
    from cemdisp.data import rheology_vs_temperature as R
    p = RheologyFormulaParams()
    assert R._quad6(p.sp_2p05_ty, 200.0, 0.1) < 0.0, "角落应仍为负（系数未动）"
    assert R._quad6(p.sp_2p05_ty, 170.0, 0.1) > 0.0, "170°C 仍为正"
    with pytest.raises(ValueError):
        # FluidSpec 构造校验拦截（负 yield_stress_pa）⇒ 不会被静默放进求解器
        fluid_at(_spacer("隔离液1", 2.05), 200.0)
    assert R._quad6(p.sp_2p05_ty, 155.0, 150.0) > 0.0, "三重点井物理域安全"


# --------------------------------------------------------------------------- #
# 5. 审计恒等（A1.3）
# --------------------------------------------------------------------------- #

def test_audit_sequence_identical_for_params_none_vs_omitted():
    """``params=None`` 与"不传 params"的审计事件**逐条相同**（kind/fluid/detail/extra）。"""
    for f in _ALL_FLUIDS:
        for T in (20.0, 60.0, 100.0, 150.0):
            for P in (None, 0.1, 70.0):
                reset_audit()
                fluid_at(f, T, P)
                a = get_audit()
                reset_audit()
                fluid_at(f, T, P, params=None)
                b = get_audit()
                assert a == b, (f.name, T, P)
                reset_audit()
                fluid_at(f, T, P, params=RheologyFormulaParams())
                c = get_audit()
                assert a == c, (f.name, T, P)


# --------------------------------------------------------------------------- #
# 6. run_opts 装配面（A1.2 + 评审 M5/M7）
# --------------------------------------------------------------------------- #

def test_run_opts_keys_equal_normalize_default_keys():
    """A1.2：``RUN_OPTS_KEYS`` 必须与 `normalize_run_opts` 的默认键集一致（防两处漂移）。"""
    assert tuple(normalize_run_opts(None)) == tuple(RUN_OPTS_KEYS)


def test_kwargs_from_opts_keys_are_valid_solver_params():
    """A1.2 + M5/M7：kwargs 映射返回的**键名**必须逐一存在于对应 solver 的 `__init__`。

    仓库签名闸门 `check_call_signatures.find_bad_kwargs` 对 ``**f(opts)``（函数返回的
    dict）是盲区 ⇒ 本测试补闸。这里同时钉死 run_opts 键
    ``rheology_formula`` → 求解器形参 ``rheology_formula_params`` 的**显式改名**。
    """
    opts = normalize_run_opts({"enable_temperature_rheology": True,
                               "temperature_mode": "static"})
    ck = casing_kwargs_from_opts(opts)
    ak = annulus_kwargs_from_opts(opts)
    assert "rheology_formula_params" in ck and "rheology_formula" not in ck
    assert "rheology_formula_params" in ak and "rheology_formula" not in ak
    casing_params = set(inspect.signature(CasingFlowSolver.__init__).parameters)
    annulus_params = set(inspect.signature(AnnulusD2DGASolver.__init__).parameters)
    assert set(ck) <= casing_params, f"坏键：{set(ck) - casing_params}"
    assert set(ak) <= annulus_params, f"坏键：{set(ak) - annulus_params}"
    assert ck["rheology_formula_params"] is None
    assert ck["mud_extrapolate"] is False
    assert ak["rheology_formula_params"] is None


def test_normalize_rejects_formula_without_temperature_switch():
    """交叉守卫：温度关时给 `rheology_formula` ⇒ 响亮报错（防静默空转/标签失真）。"""
    with pytest.raises(ValueError):
        normalize_run_opts({"rheology_formula": {"scale": {"ca_tau0_q": 1.1}}})
    with pytest.raises(TypeError):
        normalize_run_opts({"enable_temperature_rheology": True,
                            "temperature_mode": "static", "mud_extrapolate": "yes"})


# --------------------------------------------------------------------------- #
# 7. solver 接线："添加了 ≠ 在运行"（A1.2）
# --------------------------------------------------------------------------- #

def test_ctor_params_reach_phase_props_annulus():
    """ctor 新参**到达** `_phase_props`：扰动档经求解器物性唯一入口真生效。"""
    p = params_from_spec({"scale": {"ca_tau0_q": 1.1}})
    solver = AnnulusD2DGASolver(dt=2.0, nz=24, ny=12, total_t=240.0,
                                enable_temperature_rheology=True,
                                rheology_formula_params=p)
    assert solver.rheology_formula_params is p
    geom = solver._build_geom(_well_spec())
    lead = _cement("领浆", 2.10)
    got = solver._phase_props(lead, geom, 0.0)
    assert got == fluid_at(lead, 60.0, params=p)
    plain = AnnulusD2DGASolver(dt=2.0, nz=24, ny=12, total_t=240.0,
                               enable_temperature_rheology=True)
    assert plain.rheology_formula_params is None
    assert plain._phase_props(lead, geom, 0.0) == fluid_at(lead, 60.0)


def test_ctor_params_reach_phase_props_casing():
    """1D 同款接线断言（两 solver 必须同口径）。"""
    p = params_from_spec({"scale": {"ca_tau0_q": 1.1}})
    solver = CasingFlowSolver(enable_temperature_rheology=True,
                              rheology_formula_params=p)
    assert solver.rheology_formula_params is p
    lead = _cement("领浆", 2.10)
    assert solver._phase_props(lead) == fluid_at(lead, 60.0, params=p)
    plain = CasingFlowSolver(enable_temperature_rheology=True)
    assert plain._phase_props(lead) == fluid_at(lead, 60.0)


def test_mud_extrapolate_reaches_phase_props():
    """`mud_extrapolate` 经 ctor 到达派生：T>80°C 时外推档与 clamp 档**不同**。"""
    mud = _mud()
    solver_on = AnnulusD2DGASolver(dt=2.0, nz=24, ny=12, total_t=240.0,
                                   enable_temperature_rheology=True,
                                   temperature_rheology_t_c=150.0,
                                   mud_extrapolate=True)
    geom = solver_on._build_geom(_well_spec())
    got = solver_on._phase_props(mud, geom, 0.0)
    assert got == fluid_at(mud, 150.0, mud_extrapolate=True)
    assert got.yield_stress_pa != fluid_at(mud, 150.0).yield_stress_pa


def test_t_off_path_bitwise_unchanged_by_new_kwargs():
    """关2：T-off 时新参**零痕迹**（`_phase_props` 恒等返回，memo 空）。"""
    for kw in ({}, {"rheology_formula_params": params_from_spec(
        {"scale": {"ca_tau0_q": 1.1}})}, {"mud_extrapolate": True}):
        solver = AnnulusD2DGASolver(dt=2.0, nz=24, ny=12, total_t=240.0, **kw)
        geom = solver._build_geom(_well_spec())
        assert "T" not in geom
        lead = _cement("领浆", 2.10)
        assert solver._phase_props(lead, geom, 0.0) is lead
        assert solver._phase_memo == {}
