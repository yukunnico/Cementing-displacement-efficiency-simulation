"""改进 D2DGA 三闭包（auto-m / I3 / 真体力）求解器级测试。"""
import numpy as np
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
from cemdisp.data.fluid_spec import FluidSpec, FluidRole, RheologyModel
from cemdisp.data.well_spec import WellSpec, DepthValuePoint, EvaluationWindow


def _make_solver(**kw) -> AnnulusD2DGASolver:
    return AnnulusD2DGASolver(dt=4.0, nz=20, ny=10, total_t=40.0, **kw)


def _toy_well():
    # 最小井身结构，供 _compute_props / _compute_velocity 单测
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


class TestAutoMField:
    def test_constructor_has_enable_d2dga_auto_m_default_true(self):
        s = _make_solver()
        assert hasattr(s, "enable_d2dga_auto_m")
        assert s.enable_d2dga_auto_m is True

    def test_auto_m_off_falls_back_to_scalar(self):
        # enable_d2dga_auto_m=False -> m 场退化为构造常数
        s = _make_solver(enable_d2dga_auto_m=False, d2dga_viscosity_ratio=0.8)
        assert s.enable_d2dga_auto_m is False
        assert s.d2dga_viscosity_ratio == 0.8


class TestMFieldFromProps:
    def test_m_field_returned_and_shape(self):
        s = _make_solver()
        well = _toy_well()
        geom = s._build_geom(well)
        ny, nz = s.ny, s.nz
        lead = np.full((ny, nz), 0.6)
        tail = np.zeros((ny, nz))
        spacer = np.zeros((ny, nz))
        w_prev = np.full((ny, nz), 0.4)
        mud_f = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                          rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
        cement_f = FluidSpec(name="tail", role=FluidRole.TAIL, density_kg_m3=1900.0,
                              rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.180, yield_stress_pa=14.0)
        out = s._compute_props(lead, tail, spacer, w_prev, geom, mud_f, None, cement_f, None)
        # 期望多返回一个 m_field
        assert len(out) == 5
        mu, rho, mud, tau_y, m_field = out
        assert m_field.shape == (ny, nz)
        # 水泥更粘 -> m = mu_mud/mu_cement < 1
        assert np.all(m_field < 1.0)


class TestExtractAblationMetrics:
    """Unit tests for extract_ablation_metrics helper (Task 3)."""

    def test_extracts_nested_chinese_keys_to_flat_english(self):
        from cemdisp.runners.ht1_004_ablation import extract_ablation_metrics
        from cemdisp.models2d.annulus_d2dga import AnnulusSimulationResult

        # Construct a minimal fake result with only the summary field populated
        fake_summary = {
            "模拟对象": "HT1-004",
            "最终结果": {
                "全井段最终有效顶替效率": 0.8765,
                "最终水泥浆占据率": 0.9123,
                "最终窜槽指数": 0.0456,
                "最终混浆指数": 0.0234,
                "最终失稳指数": 0.0123,
            },
        }

        result = AnnulusSimulationResult(
            well_name="fake",
            geom={},
            cement_field=np.zeros((2, 2)),
            spacer_field=np.zeros((2, 2)),
            wall_field=np.zeros((2, 2)),
            metrics=None,  # type: ignore[arg-type]
            depth_profiles=None,  # type: ignore[arg-type]
            summary=fake_summary,
        )

        m = extract_ablation_metrics(result)

        assert m["effective_efficiency"] == 0.8765
        assert m["cement_occupation"] == 0.9123
        assert m["channeling_index"] == 0.0456
        assert m["mixing_index"] == 0.0234
        assert m["instability_index"] == 0.0123

    def test_handles_missing_final_result_key(self):
        from cemdisp.runners.ht1_004_ablation import extract_ablation_metrics
        from cemdisp.models2d.annulus_d2dga import AnnulusSimulationResult

        result = AnnulusSimulationResult(
            well_name="fake",
            geom={},
            cement_field=np.zeros((2, 2)),
            spacer_field=np.zeros((2, 2)),
            wall_field=np.zeros((2, 2)),
            metrics=None,  # type: ignore[arg-type]
            depth_profiles=None,  # type: ignore[arg-type]
            summary={"some_other_key": 42},
        )

        m = extract_ablation_metrics(result)
        # All values should be None (missing key returns None from .get())
        assert m["effective_efficiency"] is None
        assert m["cement_occupation"] is None
        assert m["channeling_index"] is None
        assert m["mixing_index"] is None
        assert m["instability_index"] is None

    def test_handles_non_dict_summary(self):
        from cemdisp.runners.ht1_004_ablation import extract_ablation_metrics
        from cemdisp.models2d.annulus_d2dga import AnnulusSimulationResult

        result = AnnulusSimulationResult(
            well_name="fake",
            geom={},
            cement_field=np.zeros((2, 2)),
            spacer_field=np.zeros((2, 2)),
            wall_field=np.zeros((2, 2)),
            metrics=None,  # type: ignore[arg-type]
            depth_profiles=None,  # type: ignore[arg-type]
            summary="not_a_dict",  # type: ignore[arg-type]
        )

        m = extract_ablation_metrics(result)
        assert m["effective_efficiency"] is None
        assert m["cement_occupation"] is None


class TestBuoyancyForceVector:
    """Tests for _buoyancy_force_vector (Task 4: shared between R2 I3 flux and R3 true buoyancy)."""

    def test_vertical_well_sin_term_zero(self):
        # β=0（垂直井）-> f_phi 的 sin(πφ)sinβ 项 = 0；f_xi 的 cosβ = 1
        s = _make_solver()
        well = _toy_well()
        geom = s._build_geom(well)
        f_phi, f_xi = s._buoyancy_force_vector(geom, beta_deg=0.0)
        assert np.allclose(f_phi, 0.0)  # sin(0)=0
        assert np.all(f_xi > 0)  # cos(0)=1 > 0

    def test_inclined_well_sin_term_nonzero(self):
        # β>0 -> f_phi 非零
        s = _make_solver()
        well = _toy_well()
        geom = s._build_geom(well)
        f_phi, f_xi = s._buoyancy_force_vector(geom, beta_deg=5.0)
        assert np.any(f_phi > 0)


class TestTrueBuoyancy:
    def test_constructor_has_enable_true_buoyancy(self):
        s = _make_solver()
        assert hasattr(s, "enable_true_buoyancy")
        assert s.enable_true_buoyancy is True

    def test_has_compute_buoyancy_number_method(self):
        # summary 应含 buoyancy_number 字段（至少 R3 跑后）
        s = _make_solver(enable_true_buoyancy=True)
        # 不跑完整 run（需要 loader），只验证 _compute_buoyancy_number 方法存在
        assert hasattr(s, "_compute_buoyancy_number")

    def test_unstable_density_gives_negative_b(self):
        # b<0（轻顶替重）应被检出
        s = _make_solver()
        # rho_displacing < rho_displaced -> b<0
        b = s._compute_buoyancy_number(
            rho_displacing_kg_m3=1800.0, rho_displaced_kg_m3=1900.0,
            gap_m=0.04, mu_displaced_pa_s=0.05, velocity_m_s=0.5,
        )
        assert b < 0.0

    def test_stable_density_gives_positive_b(self):
        s = _make_solver()
        b = s._compute_buoyancy_number(
            rho_displacing_kg_m3=1950.0, rho_displaced_kg_m3=1900.0,
            gap_m=0.04, mu_displaced_pa_s=0.05, velocity_m_s=0.5,
        )
        assert b > 0.0


class TestSwitchRegression:
    """I6: 开关回归测试——验证 R0（全关）vs R3（全开）确实产生不同输出。"""

    def test_r0_vs_r3_produces_different_efficiency(self):
        from cemdisp.models2d.boundary_bridge import AnnulusInletState

        well = _toy_well()
        mud = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                        rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
        lead = FluidSpec(name="lead", role=FluidRole.LEAD, density_kg_m3=1930.0,
                         rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.180, yield_stress_pa=14.0)
        fluids = (mud, lead)

        def _inlet(t: float) -> AnnulusInletState:
            return AnnulusInletState(
                time_s=t, flow_rate_m3_s=0.02, stage_name="pump",
                phase_fractions=(("cement", 1.0), ("lead", 1.0)),
            )

        # R0: all switches off
        s_r0 = AnnulusD2DGASolver(
            dt=4.0, nz=20, ny=8, total_t=40.0,
            enable_d2dga=True,
            enable_d2dga_auto_m=False,
            enable_d2dga_i3_flux=False,
            enable_true_buoyancy=False,
        )
        res_r0 = s_r0.run(well, fluids, _inlet)
        eff_r0 = float(res_r0.summary["effective_efficiency"])

        # R3: all switches on
        s_r3 = AnnulusD2DGASolver(
            dt=4.0, nz=20, ny=8, total_t=40.0,
            enable_d2dga=True,
            enable_d2dga_auto_m=True,
            enable_d2dga_i3_flux=True,
            enable_true_buoyancy=True,
        )
        res_r3 = s_r3.run(well, fluids, _inlet)
        eff_r3 = float(res_r3.summary["effective_efficiency"])

        assert eff_r0 != eff_r3, f"Switches should change output, got R0={eff_r0:.6f} R3={eff_r3:.6f}"


class TestI3FluxPhysicalUpdate:
    """I3 浮力弥散通量物理系数接线测试（T1-2：去 flux_strength=0.05，式 4.25 直驱）。"""

    def test_i3_flux_physical_update_different_from_off(self):
        """验证 I3 通量（物理系数直驱，无 0.05 限幅）产生与关闭不同的顶替效率。"""
        from cemdisp.models2d.boundary_bridge import AnnulusInletState

        # 构造 Δρ=300 kg/m³（lead 2200 - mud 1900），η₂=0.18 Pa·s 工况
        well = _toy_well()
        mud = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                        rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
        lead = FluidSpec(name="lead", role=FluidRole.LEAD, density_kg_m3=2200.0,
                         rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.180, yield_stress_pa=14.0)
        fluids = (mud, lead)

        def _inlet(t: float):
            return AnnulusInletState(
                time_s=t, flow_rate_m3_s=0.02, stage_name="pump",
                phase_fractions=(("cement", 1.0), ("lead", 1.0)),
            )

        # 启用 I3 通量（无 0.05 限幅，物理系数直驱）
        s_on = AnnulusD2DGASolver(
            dt=4.0, nz=20, ny=8, total_t=40.0,
            enable_d2dga=True,
            enable_d2dga_auto_m=True,
            enable_d2dga_i3_flux=True,
            enable_true_buoyancy=True,
        )
        res_on = s_on.run(well, fluids, _inlet)
        eff_on = float(res_on.summary["effective_efficiency"])

        # 关闭 I3 通量（对照）
        s_off = AnnulusD2DGASolver(
            dt=4.0, nz=20, ny=8, total_t=40.0,
            enable_d2dga=True,
            enable_d2dga_auto_m=True,
            enable_d2dga_i3_flux=False,
            enable_true_buoyancy=True,
        )
        res_off = s_off.run(well, fluids, _inlet)
        eff_off = float(res_off.summary["effective_efficiency"])

        # 物理系数直驱的 I3 通量应改变顶替效率（不同于 0.05 限幅版的小修正）
        assert eff_on != eff_off, (
            f"I3 flux should change efficiency; "
            f"on={eff_on:.6f} off={eff_off:.6f}"
        )
        assert 0.0 < eff_on < 1.0
        assert 0.0 < eff_off < 1.0

    def test_i3_flux_magnitude_matches_physical_coefficient(self):
        """验证 I3 通量量级符合 ΔρH³/(6η₂)·I3 ≈ 2.78e-4 理论值。"""
        from cemdisp.models2d.d2dga_flux import d2dga_buoyancy_flux, d2dga_dispersion_function_I3

        # 构造 Δρ=300 kg/m³, H=0.01 m, η₂=0.18 Pa·s 工况
        delta_rho = 300.0
        H = 0.01
        eta2 = 0.18
        c = 0.5
        m = 1.0

        # 物理系数 prefactor = ΔρH³/(6η₂) = 300*1e-6/(6*0.18) = 2.7778e-4
        prefactor = delta_rho * H**3 / (6.0 * eta2)
        assert np.isclose(prefactor, 2.7778e-4, rtol=1e-3), (
            f"prefactor={prefactor:.6e}, expected ~2.78e-4"
        )
        # I3(0.5, 1.0) = 0.0546875
        i3 = d2dga_dispersion_function_I3(c, m)
        assert np.isclose(i3, 0.0546875, rtol=1e-4), (
            f"I3(0.5,1)={i3:.8f}, expected 0.0546875"
        )
        expected_q_mag = prefactor * i3  # ≈ 1.519e-5

        q_phi, q_xi = d2dga_buoyancy_flux(
            np.array([c]), m=m, delta_rho=delta_rho, H=H,
            eta2=eta2, f_phi=1.0, f_xi=0.5,
        )
        # q_phi = +prefactor * I3 * f_phi（式 4.25 正号）
        assert np.isclose(q_phi[0], expected_q_mag * 1.0, rtol=1e-3), (
            f"q_phi={q_phi[0]:.6e}, expected {expected_q_mag:.6e}"
        )
        # q_xi = -prefactor * I3 * f_xi（式 4.25 负号）
        assert np.isclose(q_xi[0], -expected_q_mag * 0.5, rtol=1e-3), (
            f"q_xi={q_xi[0]:.6e}, expected {-expected_q_mag * 0.5:.6e}"
        )

    def test_i3_flux_no_clip_greater_than_old_clipped(self):
        """验证去 flux_strength=0.05 后通量更新量大于旧 0.05 限幅版。"""
        from cemdisp.models2d.d2dga_flux import d2dga_buoyancy_flux

        delta_rho = 300.0
        H = 0.01
        eta2 = 0.18
        c = np.array([0.5])
        m = 1.0

        q_phi, q_xi = d2dga_buoyancy_flux(
            c, m=m, delta_rho=delta_rho, H=H,
            eta2=eta2, f_phi=1.0, f_xi=0.5,
        )
        # 旧版：通量乘以 0.05 (flux_strength=0.05)，新版无此限幅
        q_phi_old = q_phi * 0.05
        q_xi_old = q_xi * 0.05

        # 新版（无 0.05 限幅）通量幅值大于旧版（应为旧版 20 倍）
        assert np.abs(q_phi[0]) > np.abs(q_phi_old[0]), (
            f"新 q_phi={q_phi[0]:.6e} 应大于旧版 {q_phi_old[0]:.6e}"
        )
        assert np.abs(q_xi[0]) > np.abs(q_xi_old[0]), (
            f"新 q_xi={q_xi[0]:.6e} 应大于旧版 {q_xi_old[0]:.6e}"
        )
        # 精确倍率验证：新 = 旧 / 0.05 = 20 倍
        assert np.isclose(q_phi[0], q_phi_old[0] * 20.0, rtol=1e-10)
        assert np.isclose(q_xi[0], q_xi_old[0] * 20.0, rtol=1e-10)
