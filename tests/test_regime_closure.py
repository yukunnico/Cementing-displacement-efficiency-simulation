import numpy as np
from cemdisp.models2d import regime_closure as rc
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
from cemdisp.data.fluid_spec import FluidSpec, FluidRole, RheologyModel
from cemdisp.data.well_spec import WellSpec, DepthValuePoint, EvaluationWindow


def _regime_solver(**kw) -> AnnulusD2DGASolver:
    return AnnulusD2DGASolver(dt=4.0, nz=20, ny=10, total_t=40.0, **kw)


def _regime_well() -> WellSpec:
    pts = lambda d, v: DepthValuePoint(depth_md_m=d, value=v)
    return WellSpec(
        well_name="toy",
        top_md_m=1000.0, bottom_md_m=1100.0, shoe_md_m=1100.0, hanger_md_m=1000.0,
        casing_id_mm=200.0, liner_od_mm=139.7, liner_id_mm=108.0,
        hole_diameter_profile=[pts(1000.0, 215.9), pts(1100.0, 215.9)],
        inclination_profile=[pts(1000.0, 5.0), pts(1100.0, 5.0)],
        standoff_profile=[pts(1000.0, 0.83), pts(1100.0, 0.83)],
        evaluation_windows=[EvaluationWindow(name="w", top_md_m=1000.0, bottom_md_m=1100.0, window_type="full")],
    )


def _call_velocity(solver: AnnulusD2DGASolver, q_m3s: float, w_prev_val: float = 0.3,
                   wall: np.ndarray | None = None, bingham: bool = False):
    """最小装配直调 _compute_velocity，返回 w 场（Task 7 solver 级测试）。"""
    well = _regime_well()
    geom = solver._build_geom(well)
    ny, nz = solver.ny, solver.nz
    lead = np.full((ny, nz), 0.6)
    tail = np.zeros((ny, nz))
    spacer = np.zeros((ny, nz))
    flusher = np.zeros((ny, nz))
    w_prev = np.full((ny, nz), w_prev_val)
    if bingham:
        mud_f = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                          rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
        lead_f = FluidSpec(name="lead", role=FluidRole.LEAD, density_kg_m3=1900.0,
                           rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.180, yield_stress_pa=14.0)
    else:
        mud_f = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                          rheology_model=RheologyModel.POWER_LAW, power_law_n=0.7, consistency_k=0.4)
        lead_f = FluidSpec(name="lead", role=FluidRole.LEAD, density_kg_m3=1900.0,
                           rheology_model=RheologyModel.POWER_LAW, power_law_n=0.7, consistency_k=0.4)
    out = solver._compute_velocity(
        lead, tail, spacer, flusher, geom, q_m3s, w_prev, mud_f, lead_f, None, None, wall=wall,
    )
    return out[0]  # w


class TestRegimeSplitSolver:
    """Task7 M2: _compute_velocity 内接流态固定点迭代（开关默认关）。"""

    def test_constructor_defaults_exist_and_off(self):
        s = _regime_solver()
        assert s.enable_regime_split is False
        assert s.regime_relax_alpha == 0.5
        assert s.regime_max_iter == 24
        assert s.regime_tol_rel == 1e-3
        assert s.regime_re_turb_ratio == 1.8

    def test_off_equals_on_in_laminar_limit(self):
        """极小排量全层流（Re<re_crit，R=1）：门开门关 w 逐位一致（层流元/壁面层不变）。"""
        s_off = _regime_solver()
        # 收敛到固定点需要足够迭代+紧容差（欠松弛 α=0.5 线性收敛）
        s_on = _regime_solver(enable_regime_split=True, regime_max_iter=200, regime_tol_rel=1e-14)
        wall = np.zeros((s_off.ny, s_off.nz))
        wall[:, : s_off.nz // 2] = 1.0  # 一半深度列带壁面静止层
        w_off = _call_velocity(s_off, q_m3s=1e-4, wall=wall)
        w_on = _call_velocity(s_on, q_m3s=1e-4, wall=wall)
        np.testing.assert_allclose(w_on, w_off, atol=1e-9)
        # 壁面层元两边都是精确零（层流元不变）
        assert np.all(np.abs(w_on[:, : s_off.nz // 2]) < 1e-12)
        assert np.all(w_on[:, s_off.nz // 2:] > 0)

    def test_off_equals_on_in_laminar_limit_defaults_no_crash(self):
        """门开默认参数（max_iter=24, tol=1e-3）在层流极限不崩，结果有限且物理合理。
        注：欠松弛(α=0.5)下 24 次迭代从 w_prev 起步未必收敛到 1e-9（见收敛版测试），
        此处仅验证默认参数路径可运行、守恒仍由面积归一保证。"""
        s_on = _regime_solver(enable_regime_split=True)
        w_on = _call_velocity(s_on, q_m3s=1e-4, wall=None)
        assert np.all(np.isfinite(w_on))
        assert np.all(w_on > 0)          # 无壁面层，全为正

    def test_off_equals_on_in_laminar_limit_wall_none(self):
        """wall=None 分支：门开门关层流极限一致（ON 分支必须处理 wall=None）。"""
        s_off = _regime_solver()
        s_on = _regime_solver(enable_regime_split=True, regime_max_iter=200, regime_tol_rel=1e-14)
        w_off = _call_velocity(s_off, q_m3s=1e-4, wall=None)
        w_on = _call_velocity(s_on, q_m3s=1e-4, wall=None)
        np.testing.assert_allclose(w_on, w_off, atol=1e-9)

    def test_split_conserves_total_flow(self):
        """门开高排量（部分元 Re>re_crit 非层流）：逐列 2·Σw·b·dy == q_half（area_weight 归一）。"""
        s_on = _regime_solver(enable_regime_split=True, regime_max_iter=200, regime_tol_rel=1e-14)
        q = 0.05
        w = _call_velocity(s_on, q_m3s=q, w_prev_val=2.0, wall=None, bingham=False)
        geom = s_on._build_geom(_regime_well())
        b = geom["effective_b"]
        dy = np.gradient(geom["y"])[:, None]
        flow_per_col = np.sum(w * b * dy * 2.0, axis=0)
        np.testing.assert_allclose(flow_per_col, q / 2.0, atol=1e-9)

    def test_split_changes_w_when_turbulent(self):
        """高排量使部分元 Re>re_crit：R≠1 令 w 相对层流归一发生可感变化（证明开关在工作）。"""
        s_off = _regime_solver()
        s_on = _regime_solver(enable_regime_split=True, regime_max_iter=200, regime_tol_rel=1e-14)
        q = 0.15
        w_off = _call_velocity(s_off, q_m3s=q, w_prev_val=2.0, wall=None)
        w_on = _call_velocity(s_on, q_m3s=q, w_prev_val=2.0, wall=None)
        assert not np.allclose(w_on, w_off, rtol=1e-3)

    def test_default_params_conserves_total_flow(self):
        """门开默认参数（max_iter=24, tol=1e-3, α=0.5）命中过渡/湍流 Re：逐列精确守恒。
        性质由 fix #1 保证——迭代只收敛 R，最终 w 以该 R 直接归一（2·Σw·b·dy == q_half），
        对任意迭代步数成立；旧 w_k（欠松弛阻尼迭代量）返回方式在默认 24 步下仍有
        相对 ~1e-3 量级（实测 8.7e-4）的瞬态守恒误差，此测试会将其钉死。"""
        s_on = _regime_solver(enable_regime_split=True)  # 全默认：max_iter=24, tol=1e-3, alpha=0.5
        assert s_on.regime_max_iter == 24
        q = 0.05
        w = _call_velocity(s_on, q_m3s=q, w_prev_val=2.0, wall=None, bingham=False)
        geom = s_on._build_geom(_regime_well())
        b = geom["effective_b"]
        dy = np.gradient(geom["y"])[:, None]
        flow_per_col = 2.0 * np.sum(w * b * dy, axis=0)
        np.testing.assert_allclose(flow_per_col, q / 2.0, rtol=1e-9)  # q_half


def test_metzner_reed_re_shape_and_positive():
    w = np.full(4, 0.5); rho = np.full(4, 1500.0); n = np.full(4, 0.7)
    kappa = np.full(4, 0.4); b = np.full(4, 0.02)
    re = rc.metzner_reed_re(w, rho, n, kappa, b)
    assert re.shape == (4,) and np.all(re > 0)


def test_laminar_friction_is_24_over_re_when_no_yield():
    re = np.array([10.0, 100.0, 1000.0]); he = np.zeros(3); n = np.ones(3)
    f = rc.friction_laminar(re, he, n)
    np.testing.assert_allclose(f, 24.0 / re, rtol=1e-6)


def test_dodge_metzner_converges_and_monotone():
    re = np.logspace(3, 6, 50); n = np.full(50, 0.7)
    f = rc.friction_dodge_metzner(re, n)
    assert f.shape == (50,) and np.all(f > 0) and np.all(f < 1.0)
    assert np.all(np.diff(f) <= 1e-6)


def test_drag_weight_layer_mask_and_bounds():
    re = np.array([10.0, 500.0, 5000.0, 50000.0])
    he = np.zeros(4); n = np.full(4, 0.8)
    R, mask = rc.drag_weight(re, he, n, re_crit=np.full(4, 2100.0))
    assert R[0] == 1.0
    assert np.all(R >= 0.3) and np.all(R <= 1.5)
    assert mask.shape == (4,)


def test_hedstrom_zero_for_no_yield():
    tau_y = np.zeros(3); rho = np.full(3, 1500.0); n = np.full(3, 0.7)
    kappa = np.full(3, 0.4); b = np.full(3, 0.02)
    assert np.all(rc.hedstrom_number(tau_y, rho, n, kappa, b) == 0.0)


def test_laminar_friction_independent_of_he_and_24_over_re():
    # 屈服在层流分支不反向削减摩擦；f=24/Re 对任意 he 成立
    re = np.array([10.0, 100.0, 2100.0])
    for he_val in (0.0, 2.0, 10.0, 76.0):
        he = np.full_like(re, he_val)
        f = rc.friction_laminar(re, he, n=np.ones_like(re))
        np.testing.assert_allclose(f, 24.0 / re, rtol=1e-12)


def test_drag_weight_not_pinned_to_floor_for_yield_fluid():
    # he>0 的生产工况下 R 不应恒触 0.3 下限（回归旧 bug）
    re = np.array([2100.0, 5000.0, 50000.0])
    he = np.array([10.0, 10.0, 10.0]); n = np.full(3, 0.7)
    re_crit = 2100.0 * (1.0 + 0.1 * he)
    R, mask = rc.drag_weight(re, he, n, re_crit=re_crit)
    assert R[0] == 1.0                       # 层流元精确为 1
    assert np.all(R > 0.3)                   # 没有任何元被钉在 0.3 下限
    assert not np.allclose(R[1:], 0.3)       # 湍流/过渡元对 Re 敏感非常数
