# STATUS: history —— 本文件锁定“重构前”的历史行为契约。
# 源模型口径重构（2026-09-14）后，相关断言预期变红；红灯不等同回归失败，
# 需按 docs/superpowers/plans/2026-09-14-d2dga-source-fidelity-and-repo-cleanup.md 逐条复核。
"""Zhang & Frigaard (2022, JFM 947:A32) 基准算例对照测试。

对照目标：论文 Table 1（物性）、Table 2（无量纲参数）、Table 3（突破时间与顶替效率）。

诚实定位
--------
论文 Table 3 的 "D2DGA" 列本身是 FCT 有限体积数值解，本模型是半拉格朗日平流 +
经验拉普拉斯弥散（``dispersion_axial=0.018`` / ``dispersion_azimuthal=0.015``）。
因此本测试只能支撑 **同源算例半定量交叉验证**，不能声称"数值正确性验证通过"。
容差先验：η_E ±0.05、t_br ±0.10（论文量纲一单位）。
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import cemdisp.models2d.annulus_d2dga as _ann_mod
from cemdisp.runners.zhang2022_benchmark import (
    ZHANG2022_CASES,
    ZHANG2022_CASE_BY_ID,
    ZHANG2022_D_HALF_M,
    ZHANG2022_G,
    ZHANG2022_INNER_RADIUS_M,
    ZHANG2022_L_M,
    ZHANG2022_OUTER_RADIUS_M,
    build_case_fluids,
    build_case_solver,
    build_case_well_spec,
    build_inlet_provider,
    case_annulus_volume_m3,
    case_buoyancy_number,
    case_mean_velocity_m_s,
    case_viscosity_ratio,
    clip_trigger_bound,
    mass_conservation_error,
    paper_eta_e,
    paper_t_br,
    peak_volume_ratio,
    run_case,
)


# --------------------------------------------------------------------------
# Table 2 / Table 1 逐位核对
# --------------------------------------------------------------------------
def test_case_table_matches_paper():
    """Table 2 的 10 个算例 (e, Re, m, b) 逐位核对论文。"""
    assert len(ZHANG2022_CASES) == 10
    assert ZHANG2022_CASES[0]["e"] == 0.8
    assert ZHANG2022_CASES[0]["b"] == -50.0
    assert ZHANG2022_CASES[3]["b"] == 1000.0


def test_table2_all_rows_match_paper():
    """Table 2 全 10 行 (e, Re, m, b) 逐位。"""
    paper_rows = [
        (1, 0.8, 20, 0.2, -50.0),
        (2, 0.8, 20, 0.2, 100.0),
        (3, 0.6, 20, 0.2, 10.0),
        (4, 0.6, 20, 0.2, 1000.0),
        (5, 0.4, 20, 0.2, 100.0),
        (6, 0.4, 20, 5.0, 100.0),
        (7, 0.2, 20, 0.5, 10.0),
        (8, 0.2, 20, 2.0, 10.0),
        (9, 0.1, 100, 0.2, 100.0),
        (10, 0.1, 1000, 0.2, 100.0),
    ]
    got = [(c["case_id"], c["e"], c["Re"], c["m"], c["b"]) for c in ZHANG2022_CASES]
    assert got == paper_rows


def test_table1_matches_table2_definition():
    """Table 1 物性 + 论文定义式 → 逐位重现 Table 2 的 b 与 m。"""
    for case in ZHANG2022_CASES:
        assert math.isclose(case_viscosity_ratio(case), case["m"], rel_tol=1e-9, abs_tol=1e-12)
        b = case_buoyancy_number(case)
        # 论文 b 保留 2~3 位有效数字，物性也是四舍五入后的打印值 → abs 容差 0.05
        assert math.isclose(b, case["b"], rel_tol=1e-3, abs_tol=0.05), (case["case_id"], b)


def test_table1_q0_is_w0_times_annulus_area():
    """Table 1 的 Q̂0 必须等于 ŵ0 × 环空面积。

    case 9 例外：论文打印的 Q̂0 = 2.38e-4 与 ŵ0 = 0.04 差 10 倍，且与 Table 2 的
    b = 100、Re = 100 矛盾（按该 Q̂0 反推 b = 10）→ 判定为论文 Table 1 笔误，
    runner 统一取 Q0 = ŵ0·A_annulus。
    """
    area = case_annulus_volume_m3() / ZHANG2022_L_M
    for case in ZHANG2022_CASES:
        implied = case_mean_velocity_m_s(case) * area
        assert math.isclose(case["Q0"], implied, rel_tol=1e-12)
        if case["case_id"] == 9:
            assert math.isclose(case["Q0_printed"], 10.0 * implied, rel_tol=1e-2)
        else:
            assert math.isclose(case["Q0_printed"], implied, rel_tol=2e-3), case["case_id"]


def test_table3_reference_values_present():
    """Table 3 的 t_br / η_E 两列（D2DGA 与 3-D）逐位记录。"""
    expected = {
        1: (0.44, 0.33, 0.66, 0.61),
        2: (0.95, 0.95, 0.95, 0.95),
        3: (0.79, 0.67, 0.92, 0.91),
        4: (0.99, 0.98, 1.00, 0.98),
        5: (0.95, 0.93, 0.97, 0.95),
        6: (0.93, 0.94, 0.93, 0.96),
        7: (0.78, 0.70, 0.90, 0.90),
        8: (0.78, 0.70, 0.84, 0.86),
        9: (0.97, 0.96, 0.97, 0.95),
        10: (0.97, 0.96, 0.97, 0.95),
    }
    for cid, (t_d2dga, t_3d, e_d2dga, e_3d) in expected.items():
        case = ZHANG2022_CASE_BY_ID[cid]
        assert (case["t_br_d2dga"], case["t_br_3d"], case["eta_e_d2dga"], case["eta_e_3d"]) == (
            t_d2dga, t_3d, e_d2dga, e_3d,
        )


# --------------------------------------------------------------------------
# 几何
# --------------------------------------------------------------------------
def test_geometry_is_lab_scale():
    ws = build_case_well_spec(ZHANG2022_CASES[0])
    assert np.isclose(ws.bottom_md_m - ws.top_md_m, 4.8, atol=1e-9)


def test_geometry_gap_matches_paper_radii():
    """几何半间隙 d 必须等于 (r_o − r_i)/2 = 2.385 mm。"""
    assert math.isclose(ZHANG2022_D_HALF_M, (ZHANG2022_OUTER_RADIUS_M - ZHANG2022_INNER_RADIUS_M) / 2.0)
    assert math.isclose(ZHANG2022_D_HALF_M, 2.385e-3, rel_tol=1e-12)
    solver = build_case_solver(ZHANG2022_CASES[0])
    geom = solver._build_geom(build_case_well_spec(ZHANG2022_CASES[0]))
    # 求解器内部 geom["b"] 为体校正后的全环隙（= 2d），均值应等于 4.77 mm
    assert np.isclose(float(np.mean(geom["b"])), 2.0 * ZHANG2022_D_HALF_M, rtol=1e-9)


def test_geometry_annulus_volume_matches_paper():
    """环空体积 = π(r_o²−r_i²)·L。"""
    v_expected = math.pi * (ZHANG2022_OUTER_RADIUS_M**2 - ZHANG2022_INNER_RADIUS_M**2) * ZHANG2022_L_M
    assert math.isclose(case_annulus_volume_m3(), v_expected, rel_tol=1e-12)
    solver = build_case_solver(ZHANG2022_CASES[0])
    v_model = solver._physical_annular_volume(build_case_well_spec(ZHANG2022_CASES[0]))
    assert math.isclose(v_model, v_expected, rel_tol=1e-9)


def test_geometry_vertical_and_eccentricity():
    for case in ZHANG2022_CASES:
        ws = build_case_well_spec(case)
        assert all(p.value == 0.0 for p in ws.inclination_profile)  # 直井
        assert np.isclose(ws.standoff_profile[0].value, 1.0 - case["e"], atol=1e-12)


# --------------------------------------------------------------------------
# e_clip 解锁
# --------------------------------------------------------------------------
def test_e_clip_released():
    """基准算例必须解除 e_clip 0.55 截断，否则 e>=0.6 的算例几何被篡改。"""
    solver = build_case_solver(ZHANG2022_CASES[0])
    assert solver.e_clip_max == 1.0


def test_e_clip_released_geometry_keeps_eccentricity():
    """e=0.8 算例在解锁后 geom["e"] 必须仍为 0.8（默认 0.55 会被截断）。"""
    case = ZHANG2022_CASE_BY_ID[1]
    ws = build_case_well_spec(case)
    released = build_case_solver(case, nz=20)._build_geom(ws)["e"][0]
    assert np.isclose(released, 0.8, atol=1e-12)
    from cemdisp.models2d import AnnulusD2DGASolver

    default_cap = AnnulusD2DGASolver(nz=20)._build_geom(ws)["e"][0]
    assert np.isclose(default_cap, 0.55, atol=1e-12)  # 说明截断确实存在


# --------------------------------------------------------------------------
# inlet / fluids
# --------------------------------------------------------------------------
def test_inlet_provider_is_pure_cement_at_case_flow_rate():
    case = ZHANG2022_CASES[0]
    provider = build_inlet_provider(case)
    state = provider(0.0)
    assert state.flow_rate_m3_s == case["Q0"]
    assert dict(state.phase_fractions)["tail"] == 1.0


def test_fluids_are_newtonian_with_table1_properties():
    case = ZHANG2022_CASE_BY_ID[6]
    mud, cem = build_case_fluids(case)
    assert mud.density_kg_m3 == 1000.0
    assert mud.plastic_viscosity_pa_s == 0.005
    assert cem.density_kg_m3 == 1358.41
    assert cem.plastic_viscosity_pa_s == 0.001


# --------------------------------------------------------------------------
# 裁剪统计界
# --------------------------------------------------------------------------
def test_clip_bound_is_inactive_for_all_cases():
    """论文式 4.24 的 ±0.5 硬裁剪在这些算例上理论上不应触发（max|Δρ·I2/I1| << 0.5）。"""
    for case in ZHANG2022_CASES:
        ub = clip_trigger_bound(case)
        assert ub < 0.5, (case["case_id"], ub)


# --------------------------------------------------------------------------
# 端到端（用最快的算例 10：w0=0.4 → 2L/w0 = 24 s）
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def case10_run():
    return run_case(ZHANG2022_CASE_BY_ID[10], nz=140, ny=40, keep_result=True)


@pytest.fixture(scope="module")
def case5_run():
    return run_case(ZHANG2022_CASE_BY_ID[5], nz=140, ny=40)


def test_run_produces_paper_time_window(case10_run):
    """模拟终点 = 2 个环空体积注入时间 2L/ŵ0（论文口径）。"""
    res = case10_run["_result"]
    case = ZHANG2022_CASE_BY_ID[10]
    t_end = 2.0 * ZHANG2022_L_M / case_mean_velocity_m_s(case)
    assert res.metrics["time_s"].iloc[-1] == pytest.approx(t_end, rel=1e-9)


def test_paper_eta_e_in_unit_interval(case10_run):
    case = ZHANG2022_CASE_BY_ID[10]
    eta = paper_eta_e(case10_run["_result"], case10_run["_well_spec"], case_mean_velocity_m_s(case))
    assert 0.0 <= eta <= 1.0


def test_mass_conservation_is_violated_by_d2dga_flux_amplification(case10_run):
    """质量守恒检验（论文 §2.4 的硬锚点）——**本模型不满足**，且缺口可定位。

    目标 <1e-6，实测：突破前环空内水泥体积比累计注入体积**多出 ~28–33%**。
    归因（对照实验）：把 ``enable_d2dga=False`` 关掉通量放大后，同一算例的缺口降到
    ~3%，说明缺口主要来自 ``w_d2dga = w · f_amp(c)``（``annulus_d2dga.py:1168``）——
    把论文的**通量**放大函数当成**速度**乘子用在半拉格朗日平流里（f_amp 在 c→0 处
    趋于 1.5），于是前锋以 1.5·ŵ0 推进、凭空造出体积。论文口径下守恒律
    ∂c/∂t + ∂F(c)/∂ξ = 0 的 Rankine-Hugoniot 前锋速度应为 F(1)−F(0) = ŵ0。

    本测试锁定当前实测量级（characterization test），不是"通过"的断言；
    D2DGA 通量实现一旦修好，此测试会失败并提示更新。
    """
    case = ZHANG2022_CASE_BY_ID[10]
    res, ws = case10_run["_result"], case10_run["_well_spec"]
    w0 = case_mean_velocity_m_s(case)
    t_br, _ = paper_t_br(res, ws, w0)
    err_on = mass_conservation_error(res, ws, case["Q0"], t_end_s=0.9 * t_br)
    assert err_on > 0.1, f"预期 D2DGA 通量放大破坏守恒（>10%），实测 {err_on:.4f}"
    # 硬判据：环空内水泥体积峰值超过累计注入体积 ⇒ 凭空造出水泥
    peak_on = peak_volume_ratio(res, ws, case["Q0"])
    assert peak_on > 1.2, f"预期体积创造峰值 >1.2，实测 {peak_on:.3f}"

    # 对照实验：关闭 D2DGA 通量放大 → 缺口应大幅收窄、且不再造体积
    row_off = run_case(case, nz=140, ny=40, solver_overrides={"enable_d2dga": False},
                       keep_result=True)
    err_off = row_off["mass_conservation_error"]
    assert err_off < 0.05, f"关闭 D2DGA 后仍不守恒 {err_off:.4f}"
    assert peak_volume_ratio(row_off["_result"], row_off["_well_spec"], case["Q0"]) <= 1.0
    assert err_on > 5.0 * err_off, (err_on, err_off)


def test_clip_instrument_works_and_never_triggers(case10_run):
    """式 4.24 的 ±0.5 硬裁剪：插桩生效（每步一次调用）且 10 个算例均未触发。"""
    assert case10_run["clip_调用次数"] > 100
    assert case10_run["clip_triggered"] == 0
    assert case10_run["clip_最大绝对值"] < 0.5


def test_case5_eta_e_within_prior_tolerance(case5_run):
    """容差先验抽查（characterization）：case 5（e=0.4, b=100, m=0.2）。

    η_E 命中先验容差 ±0.05；t_br 超出先验容差 ±0.10（模型 0.58 vs 论文 0.95）。
    整表命中率见 results/基准算例对照_2026-09-10/对照明细.csv。
    """
    assert abs(case5_run["eta_E_模型"] - 0.97) < 0.05
    assert abs(case5_run["t_br_模型"] - 0.95) > 0.10


def test_grid_convergence_nz_140_to_500(case10_run, monkeypatch):
    """网格收敛（任务目标）：nz 140 → 500 时 η_E 变化必须 < 0.02。

    CFL 自适应使 dt 随 ds 同步细化（case 10：dt 中位 0.033 s → 0.020 s），
    故该检验同时覆盖网格与时间步收敛。

    STATUS（Task 5 机械适配，2026-09-14）：轴向浮力数接线（K_AXIAL，Z&F22 (4.14)）
    改变了 t_br 的网格敏感度（nz 140→500 的 |Δt_br| 由 <0.10 变 ~0.12）。本测试属
    "重构前"行为契约，故在测试内把 K_AXIAL 钉回 0（= Task 4 状态），并在本测试内
    重跑 nz=140 基准（共享 fixture case10_run 现为新物理口径，不能跨物理比较），
    **断言未改动**；新行为的网格收敛由 Task 12 基准算例重跑复核。
    """
    monkeypatch.setattr(_ann_mod, "K_AXIAL", 0.0)
    base = run_case(ZHANG2022_CASE_BY_ID[10], nz=140, ny=40)
    coarse = run_case(ZHANG2022_CASE_BY_ID[10], nz=500, ny=40)
    assert abs(coarse["eta_E_模型"] - base["eta_E_模型"]) < 0.02
    assert abs(coarse["t_br_模型"] - base["t_br_模型"]) < 0.10


def test_breakthrough_time_is_reported(case10_run):
    case = ZHANG2022_CASE_BY_ID[10]
    res, ws = case10_run["_result"], case10_run["_well_spec"]
    w0 = case_mean_velocity_m_s(case)
    assert case10_run["t_br_模型"] > 0.0
    t_br_raw, t_br_dimensionless = paper_t_br(res, ws, w0)
    assert t_br_raw > 0.0
    assert t_br_dimensionless == pytest.approx(t_br_raw * w0 / ZHANG2022_L_M, rel=1e-12)


def test_vertical_well_buoyancy_number_does_not_change_breakthrough(monkeypatch):
    """结构性缺口（characterization test）：竖直井 + 默认 ``enable_true_buoyancy=True`` 时
    ``_buoyancy_force_vector`` 的方位分量 ``f_phi ∝ sinβ ≡ 0``（β = 井斜 = 0），
    于是 ``buoyancy_shape ≡ 1``，b 完全不进入方位速度分布——case 1（b = −50）与
    case 2（b = +100）的 t_br / η_E 逐位相同。

    论文 Table 3：case 1 t_br = 0.44、case 2 t_br = 0.95（相差 0.51）。

    对照实验：切到旧代理路径 ``enable_true_buoyancy=False``
    （``buoyancy_shape = 1 + stable·e·(2φ−1)``）后模型恢复对 b 符号的区分能力
    （0.51 vs 0.64），但仍远弱于论文 ⇒ 缺口定位到 R3 真浮力路径在竖直井上的退化。

    STATUS（Task 5 机械适配，2026-09-14）：本测试钉住的"结构性缺口"（b 不进动力学）
    正是 Task 5 按 Z&F22 (4.14)/(4.22) 接线轴向浮力数所要修复的——接线后 case 1/2 的
    t_br 有意分离（实测 |Δ| ≈ 0.33）。本测试属"重构前"行为契约，故在测试内把
    K_AXIAL 钉回 0（= Task 4 状态）继续锁定旧路径，**断言未改动**；
    新行为（b 分离突破时间）的契约由 tests/contract/test_buoyancy_axial.py 锁定。
    """
    monkeypatch.setattr(_ann_mod, "K_AXIAL", 0.0)
    r1 = run_case(ZHANG2022_CASE_BY_ID[1], nz=140, ny=40)
    r2 = run_case(ZHANG2022_CASE_BY_ID[2], nz=140, ny=40)
    assert abs(r2["t_br_模型"] - r1["t_br_模型"]) < 0.05, (r1["t_br_模型"], r2["t_br_模型"])

    off = {"enable_true_buoyancy": False}
    o1 = run_case(ZHANG2022_CASE_BY_ID[1], nz=140, ny=40, solver_overrides=off)
    o2 = run_case(ZHANG2022_CASE_BY_ID[2], nz=140, ny=40, solver_overrides=off)
    assert o2["t_br_模型"] - o1["t_br_模型"] > 0.10, (o1["t_br_模型"], o2["t_br_模型"])
