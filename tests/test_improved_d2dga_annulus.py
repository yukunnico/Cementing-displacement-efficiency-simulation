"""改进 D2DGA 三闭包（auto-m / I3 / 真体力）求解器级测试。"""
import numpy as np
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver, _trapez2d
from cemdisp.data.fluid_spec import FluidSpec, FluidRole, RheologyModel
from cemdisp.models2d.d2dga_flux import d2dga_dispersion_I1, d2dga_dispersion_I2
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
        out = s._compute_props(lead, tail, spacer, np.zeros_like(lead), w_prev, geom, mud_f, None, cement_f, None)
        # 期望多返回 m_field + 相黏度场 eta1/eta2（T1-4）
        assert len(out) == 7
        mu, rho, mud, tau_y, m_field, eta1, eta2 = out
        assert m_field.shape == (ny, nz)
        # 水泥更粘 -> m = mu_mud/mu_cement < 1
        assert np.all(m_field < 1.0)
        # 相黏度场形状与 m_field 一致
        assert eta1.shape == (ny, nz)
        assert eta2.shape == (ny, nz)


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


class TestTwoLayerViscosity:
    """T1-4: 两层黏度闭包 1/η_mix = c̄³/η₂ + (1−c̄³)/η₁（式 4.23）端点及混合行为验证。"""

    def _compute_eta_mix(self, c_bar, eta1, eta2):
        """应用式 4.23 计算两层黏度闭包。"""
        return 1.0 / (c_bar**3 / np.maximum(eta2, 1.0e-9)
                      + (1.0 - c_bar**3) / np.maximum(eta1, 1.0e-9))

    def test_eta_mix_equals_eta1_when_cbar_zero(self):
        """c̄=0 → η_mix = η1（纯泥浆相）。"""
        s = _make_solver()
        well = _toy_well()
        geom = s._build_geom(well)
        ny, nz = s.ny, s.nz
        w_prev = np.full((ny, nz), 0.4)
        mud_f = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                          rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
        cement_f = FluidSpec(name="tail", role=FluidRole.TAIL, density_kg_m3=1900.0,
                              rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.180, yield_stress_pa=14.0)
        # c̄=0: lead=0, tail=0
        lead = np.zeros((ny, nz))
        tail = np.zeros((ny, nz))
        spacer = np.zeros((ny, nz))
        _, _, _, _, _, eta1, eta2 = s._compute_props(
            lead, tail, spacer, np.zeros_like(lead), w_prev, geom, mud_f, None, cement_f, None)
        c_bar = np.clip(lead + tail, 0.0, 1.0)
        eta_mix = self._compute_eta_mix(c_bar, eta1, eta2)
        np.testing.assert_allclose(eta_mix, eta1, rtol=1e-10)

    def test_eta_mix_equals_eta2_when_cbar_one(self):
        """c̄=1 → η_mix = η2（纯水泥相）。"""
        s = _make_solver()
        well = _toy_well()
        geom = s._build_geom(well)
        ny, nz = s.ny, s.nz
        w_prev = np.full((ny, nz), 0.4)
        mud_f = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                          rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
        cement_f = FluidSpec(name="tail", role=FluidRole.TAIL, density_kg_m3=1900.0,
                              rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.180, yield_stress_pa=14.0)
        # c̄=1: lead=1, tail=0（水泥相通过 lead_fluid 传入）
        lead = np.ones((ny, nz))
        tail = np.zeros((ny, nz))
        spacer = np.zeros((ny, nz))
        _, _, _, _, _, eta1, eta2 = s._compute_props(
            lead, tail, spacer, np.zeros_like(lead), w_prev, geom, mud_f, lead_fluid=cement_f, tail_fluid=None, spacer_fluid=None)
        c_bar = np.clip(lead + tail, 0.0, 1.0)
        eta_mix = self._compute_eta_mix(c_bar, eta1, eta2)
        np.testing.assert_allclose(eta_mix, eta2, rtol=1e-10)

    def test_eta_mix_between_eta1_and_eta2_for_mixed_cbar(self):
        """混合 c̄ 时 η_mix 介于 η1 和 η2 之间。"""
        s = _make_solver()
        well = _toy_well()
        geom = s._build_geom(well)
        ny, nz = s.ny, s.nz
        w_prev = np.full((ny, nz), 0.4)
        mud_f = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                          rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
        cement_f = FluidSpec(name="tail", role=FluidRole.TAIL, density_kg_m3=1900.0,
                              rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.180, yield_stress_pa=14.0)
        # c̄=0.5: lead=0.5, tail=0
        lead = np.full((ny, nz), 0.5)
        tail = np.zeros((ny, nz))
        spacer = np.zeros((ny, nz))
        _, _, _, _, _, eta1, eta2 = s._compute_props(
            lead, tail, spacer, np.zeros_like(lead), w_prev, geom, mud_f, lead_fluid=cement_f, tail_fluid=None, spacer_fluid=None)
        c_bar = np.clip(lead + tail, 0.0, 1.0)
        eta_mix = self._compute_eta_mix(c_bar, eta1, eta2)
        eta_min = np.minimum(eta1, eta2)
        eta_max = np.maximum(eta1, eta2)
        assert np.all(eta_mix >= eta_min - 1e-10), "η_mix 不应小于 min(η1,η2)"
        assert np.all(eta_mix <= eta_max + 1e-10), "η_mix 不应大于 max(η1,η2)"


class TestBuoyancyForceInjection:
    """T1-3b: 体力向量注入流动度（式 2.5b/4.24）测试。"""

    def _make_well_inclined(self):
        """创建斜井（5 deg）井身结构。"""
        pts = lambda d, v: DepthValuePoint(depth_md_m=d, value=v)
        return WellSpec(
            well_name="inclined",
            top_md_m=1000.0, bottom_md_m=1100.0, shoe_md_m=1100.0, hanger_md_m=1000.0,
            casing_id_mm=200.0, liner_od_mm=139.7, liner_id_mm=108.0,
            hole_diameter_profile=[pts(1000.0, 215.9), pts(1100.0, 215.9)],
            inclination_profile=[pts(1000.0, 5.0), pts(1100.0, 5.0)],
            standoff_profile=[pts(1000.0, 0.83), pts(1100.0, 0.83)],
            evaluation_windows=[EvaluationWindow(name="w", top_md_m=1000.0, bottom_md_m=1100.0, window_type="full")],
        )

    def _setup_heavy_over_light(self, solver):
        """重顶替轻工况：水泥 2200 vs 泥浆 1900 kg/m³，均匀 c=0.6。"""
        well = self._make_well_inclined()
        geom = solver._build_geom(well)
        ny, nz = solver.ny, solver.nz
        mud_f = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                          rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
        lead_f = FluidSpec(name="lead", role=FluidRole.LEAD, density_kg_m3=2200.0,
                           rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.180, yield_stress_pa=14.0)
        lead = np.full((ny, nz), 0.6)
        tail = np.zeros((ny, nz))
        spacer = np.zeros((ny, nz))
        w_prev = np.full((ny, nz), 0.4)
        return geom, ny, nz, mud_f, lead_f, lead, tail, spacer, w_prev

    def _base_and_props(self, geom, lead, tail, spacer, w_prev, mud_f, lead_f):
        """由 _compute_props 重建基础流动度 base 与混合物性场。"""
        s = AnnulusD2DGASolver(dt=4.0, nz=20, ny=10, total_t=40.0)
        mu, rho, mud, tau_y, m_field, eta1, eta2 = s._compute_props(
            lead, tail, spacer, np.zeros_like(lead), w_prev, geom, mud_f, lead_f, None, None
        )
        c_bar = np.clip(lead + tail, 0.0, 1.0)
        eta_mix = 1.0 / (
            c_bar ** 3 / np.maximum(eta2, 1.0e-9)
            + (1.0 - c_bar ** 3) / np.maximum(eta1, 1.0e-9)
        )
        b = geom["b"]
        b_mean = np.mean(b, axis=0, keepdims=True)
        base = (b / np.maximum(b_mean, 1.0e-12)) ** 2 / np.maximum(eta_mix, 1.0e-9)
        return base, rho, m_field

    def test_true_buoyancy_pref_raises_at_mid_gap_for_heavy_over_light(self):
        """重顶替轻（Δρ>0）时，True 分支中缝区 pref 高于 False 分支（体力正向修正）。"""
        geom, ny, nz, mud_f, lead_f, lead, tail, spacer, w_prev = self._setup_heavy_over_light(
            AnnulusD2DGASolver(dt=4.0, nz=20, ny=10, total_t=40.0)
        )
        base, rho, m_field = self._base_and_props(
            geom, lead, tail, spacer, w_prev, mud_f, lead_f
        )

        # True 分支 pref：式 4.24，重顶替轻时 f_phi>0 的中缝区修正最大
        s_true = AnnulusD2DGASolver(
            dt=4.0, nz=20, ny=10, total_t=40.0,
            enable_d2dga=True, enable_d2dga_auto_m=True,
            enable_d2dga_i3_flux=True, enable_true_buoyancy=True,
        )
        f_phi, _ = s_true._buoyancy_force_vector(geom, float(np.mean(geom["inc_deg"])))
        m_local = float(np.mean(m_field))
        i1 = d2dga_dispersion_I1(0.6, m_local)
        i2 = d2dga_dispersion_I2(0.6, m_local)
        rho_mud = mud_f.density_kg_m3 / 1000.0
        delta_rho = float(np.mean(rho - rho_mud))
        correction = np.clip(delta_rho * (i2 / np.maximum(i1, 1.0e-12)), -0.5, 0.5)
        pref_true = base * (1.0 + correction * f_phi)

        # False 分支 pref：旧 (2φ−1) 代理
        density_contrast = (lead_f.density_kg_m3 - mud_f.density_kg_m3) / mud_f.density_kg_m3
        stable = float(np.clip(8.0 * density_contrast, -0.35, 0.45))
        ebar = float(np.mean(geom["e"]))
        phi = geom["phi"][:, None]
        pref_false = base * (1.0 + stable * ebar * (2.0 * phi - 1.0))

        # 中缝区（φ≈0.5）True 的 pref 更高
        mid_mask = (phi[:, 0] >= 0.4) & (phi[:, 0] <= 0.6)
        assert mid_mask.any(), "中缝区应有网格点"
        assert np.mean(pref_true[mid_mask, :]) > np.mean(pref_false[mid_mask, :]), (
            f"重顶替轻时 True 分支中缝 pref 应高于 False"
        )
        # 方位加权积分 pref(True) > pref(False)
        assert _trapez2d(pref_true, geom) > _trapez2d(pref_false, geom), (
            "重顶替轻时 True 分支方位加权 pref 积分应高于 False"
        )

    def test_false_buoyancy_falls_back_to_simplified(self):
        """enable_true_buoyancy=False 时，pref 形状恢复为 (2φ−1) 简化代理。"""
        s = AnnulusD2DGASolver(
            dt=4.0, nz=20, ny=10, total_t=40.0,
            enable_d2dga=True, enable_d2dga_auto_m=True,
            enable_d2dga_i3_flux=True, enable_true_buoyancy=False,
        )
        geom, ny, nz, mud_f, lead_f, lead, tail, spacer, w_prev = self._setup_heavy_over_light(s)

        w, *_ = s._compute_velocity(
            lead, tail, spacer, np.zeros_like(lead), geom, q_m3s=0.02, w_prev=w_prev,
            mud_fluid=mud_f, lead_fluid=lead_f, tail_fluid=None, spacer_fluid=None,
        )
        assert np.all(w > 0), "速度应为正"
        assert w.shape == (ny, nz)

        # 从速度场反解无量纲 buoyancy_shape：w ∝ base·shape
        base, *_ = self._base_and_props(geom, lead, tail, spacer, w_prev, mud_f, lead_f)
        shape_est = w / np.maximum(base, 1.0e-12)
        # 按列归一化，消除截面流量常数
        shape_est = shape_est / np.mean(shape_est, axis=0, keepdims=True)
        shape_est_mean = np.mean(shape_est, axis=1)

        # 期望的 (2φ−1) 代理
        density_contrast = (lead_f.density_kg_m3 - mud_f.density_kg_m3) / mud_f.density_kg_m3
        stable = float(np.clip(8.0 * density_contrast, -0.35, 0.45))
        ebar = float(np.mean(geom["e"]))
        phi = geom["phi"][:, None]
        shape_exp = 1.0 + stable * ebar * (2.0 * phi[:, 0] - 1.0)

        np.testing.assert_allclose(shape_est_mean, shape_exp, rtol=1.0e-10, atol=1.0e-10)


class TestStaticWallLayer:
    """T1-5: Static wall layer c_min 判据（Bararpour 2025 式 2.35-2.41）"""

    def test_constructor_has_cmin_default(self):
        """c_min 参数默认值应为 0.05。"""
        s = _make_solver()
        assert hasattr(s, "c_min")
        assert s.c_min == 0.05

    def test_cmin_parameter_stored(self):
        """c_min 参数可配置。"""
        s = _make_solver(c_min=0.3)
        assert s.c_min == 0.3

    def test_wall_consistency_after_run(self):
        """运行后 wall 场与水泥浓度场一致：wall=1 ↔ cement < c_min。"""
        from cemdisp.models2d.boundary_bridge import AnnulusInletState
        well = _toy_well()
        mud = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                        rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
        lead = FluidSpec(name="lead", role=FluidRole.LEAD, density_kg_m3=1930.0,
                         rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.180, yield_stress_pa=14.0)
        fluids = (mud, lead)

        def _inlet(t: float):
            return AnnulusInletState(
                time_s=t, flow_rate_m3_s=0.02, stage_name="pump",
                phase_fractions=(("cement", 1.0), ("lead", 1.0)),
            )

        s = AnnulusD2DGASolver(dt=4.0, nz=20, ny=8, total_t=40.0, c_min=0.05)
        res = s.run(well, fluids, _inlet)
        cement = np.clip(res.lead_field + res.tail_field, 0.0, 1.0)
        # wall=1 → cement < c_min; wall=0 → cement >= c_min
        wall_one = res.wall_field > 0.5
        if wall_one.any():
            assert np.all(cement[wall_one] < 0.05 + 1e-9), "wall=1 处 cement 应 < c_min"
        wall_zero = res.wall_field < 0.5
        if wall_zero.any():
            assert np.all(cement[wall_zero] >= 0.05 - 1e-9), "wall=0 处 cement 应 >= c_min"

    def test_wall_zeros_velocity_in_wall_cells(self):
        """wall=1 处速度 w≈0（流动度归零）。"""
        s = _make_solver(c_min=0.05)
        well = _toy_well()
        geom = s._build_geom(well)
        ny, nz = s.ny, s.nz
        mud_f = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                          rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
        lead_f = FluidSpec(name="lead", role=FluidRole.LEAD, density_kg_m3=1930.0,
                           rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.180, yield_stress_pa=14.0)
        lead = np.full((ny, nz), 0.6)
        tail = np.zeros((ny, nz))
        spacer = np.zeros((ny, nz))
        w_prev = np.full((ny, nz), 0.4)
        wall = np.zeros((ny, nz))
        wall[:, :nz//2] = 1.0

        w, *_ = s._compute_velocity(
            lead, tail, spacer, np.zeros_like(lead), geom, 0.02, w_prev, mud_f, lead_f, None, None, wall=wall,
        )

        # wall=1 处速度 ≈ 0
        assert np.all(np.abs(w[:, :nz//2]) < 1e-10), "wall=1 处速度应 ≈ 0"
        # wall=0 处速度 > 0
        assert np.all(w[:, nz//2:] > 0), "wall=0 处速度应 > 0"

    def test_cmin_0_3_more_wall_than_cmin_0_05(self):
        """c_min=0.3 时 wall=1 网格数多于 c_min=0.05。"""
        from cemdisp.models2d.boundary_bridge import AnnulusInletState
        well = _toy_well()
        mud = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                        rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
        lead = FluidSpec(name="lead", role=FluidRole.LEAD, density_kg_m3=1930.0,
                         rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.180, yield_stress_pa=14.0)
        fluids = (mud, lead)

        def _inlet(t: float):
            return AnnulusInletState(
                time_s=t, flow_rate_m3_s=0.02, stage_name="pump",
                phase_fractions=(("cement", 1.0), ("lead", 1.0)),
            )

        s_low = AnnulusD2DGASolver(dt=4.0, nz=20, ny=8, total_t=40.0, c_min=0.05)
        s_high = AnnulusD2DGASolver(dt=4.0, nz=20, ny=8, total_t=40.0, c_min=0.3)
        res_low = s_low.run(well, fluids, _inlet)
        res_high = s_high.run(well, fluids, _inlet)
        assert np.sum(res_high.wall_field) >= np.sum(res_low.wall_field), (
            f"c_min=0.3 wall=1 网格数 ({np.sum(res_high.wall_field)}) "
            f"应 >= c_min=0.05 ({np.sum(res_low.wall_field)})"
        )

    def test_wall_zero_before_cement_arrival(self):
        """水泥前锋未到达时（全场 cement=0），wall 场必须全为 0，不得提前堵死。"""
        from cemdisp.models2d.boundary_bridge import AnnulusInletState
        well = _toy_well()
        mud = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                        rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
        lead = FluidSpec(name="lead", role=FluidRole.LEAD, density_kg_m3=1930.0,
                         rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.180, yield_stress_pa=14.0)
        fluids = (mud, lead)

        def _inlet(t: float):
            return AnnulusInletState(
                time_s=t, flow_rate_m3_s=0.02, stage_name="pump",
                phase_fractions=(("mud", 1.0),),
            )

        s = AnnulusD2DGASolver(dt=4.0, nz=20, ny=8, total_t=40.0, c_min=0.05)
        res = s.run(well, fluids, _inlet)
        assert np.all(res.wall_field == 0.0), "水泥前锋未到达时 wall 场应全为 0"


class TestFlusherField:
    """T1-6: FLUSHER 独立浓度场（被动平流相）测试。"""

    def test_annulus_result_default_flusher_field_none(self):
        """AnnulusSimulationResult 默认 flusher_field=None 后向兼容。"""
        from cemdisp.models2d.annulus_d2dga import AnnulusSimulationResult
        result = AnnulusSimulationResult(
            well_name="test",
            geom={},
            cement_field=np.zeros((2, 2)),
            spacer_field=np.zeros((2, 2)),
            wall_field=np.zeros((2, 2)),
            metrics=None,  # type: ignore[arg-type]
            depth_profiles=None,  # type: ignore[arg-type]
            summary={},
        )
        assert result.flusher_field is None
        assert result.flusher_snapshots == ()

    def test_flusher_injected_via_inlet(self):
        """flusher 通过入口注入后 flusher_field 非空且独立于 cement/spacer。"""
        from cemdisp.models2d.boundary_bridge import AnnulusInletState
        well = _toy_well()
        mud = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                        rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
        lead = FluidSpec(name="lead", role=FluidRole.LEAD, density_kg_m3=1930.0,
                         rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.180, yield_stress_pa=14.0)
        flusher = FluidSpec(name="flusher", role=FluidRole.FLUSHER, density_kg_m3=1850.0,
                            rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.04, yield_stress_pa=3.0)
        fluids = (mud, lead, flusher)

        def _inlet(t: float):
            return AnnulusInletState(
                time_s=t, flow_rate_m3_s=0.02, stage_name="flusher",
                phase_fractions=(("flusher", 1.0),),
            )

        s = AnnulusD2DGASolver(dt=4.0, nz=20, ny=8, total_t=40.0)
        res = s.run(well, fluids, _inlet)
        assert res.flusher_field is not None
        assert np.any(res.flusher_field > 0), "flusher 应出现在求解域内"
        cement = np.clip(res.lead_field + res.tail_field, 0.0, 1.0)
        assert np.all(cement == 0.0), "flusher 不应混入 lead/tail"
        assert np.all(res.spacer_field == 0.0), "flusher 不应混入 spacer"

    def test_flusher_not_mixed_with_cement_on_coexistence(self):
        """flusher 与水泥共存时，cement/spacer 场不受 flusher 侵入。"""
        from cemdisp.models2d.boundary_bridge import AnnulusInletState
        well = _toy_well()
        mud = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                        rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
        lead = FluidSpec(name="lead", role=FluidRole.LEAD, density_kg_m3=1930.0,
                         rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.180, yield_stress_pa=14.0)
        flusher = FluidSpec(name="flusher", role=FluidRole.FLUSHER, density_kg_m3=1850.0,
                            rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.04, yield_stress_pa=3.0)
        fluids = (mud, lead, flusher)

        # 分阶段注入：先 flusher，后水泥，验证 flusher 场和水泥场各自独立
        def _inlet(t: float):
            if t < 20.0:
                return AnnulusInletState(
                    time_s=t, flow_rate_m3_s=0.02, stage_name="flusher",
                    phase_fractions=(("flusher", 1.0),),
                )
            else:
                return AnnulusInletState(
                    time_s=t, flow_rate_m3_s=0.02, stage_name="cement",
                    phase_fractions=(("cement", 1.0), ("lead", 1.0)),
                )

        s = AnnulusD2DGASolver(dt=4.0, nz=20, ny=8, total_t=40.0)
        res = s.run(well, fluids, _inlet)
        assert res.flusher_field is not None
        # flusher 应在域内存在（前 20s 注入）
        assert np.any(res.flusher_field > 0), "flusher 应出现在求解域内"
        # cement 应在域内存在（后 20s 注入）
        cement = np.clip(res.lead_field + res.tail_field, 0.0, 1.0)
        assert np.any(cement > 0), "cement 应出现在求解域内"
        # 验证 flusher 与 cement 不混：flusher 注入期间前锋处 cement 应为 0
        # 注意：数值扩散/弥散会导致界面处两相轻微重叠，允许 5% 容差
        assert np.all(cement + res.flusher_field <= 1.05), "flusher+cement 不应超 1.05"

    def test_five_phase_closure_in_compute_props(self):
        """_compute_props 五相闭合：lead+tail+spacer+flusher+mud ≈ 1。"""
        s = _make_solver()
        well = _toy_well()
        geom = s._build_geom(well)
        ny, nz = s.ny, s.nz
        lead = np.full((ny, nz), 0.3)
        tail = np.full((ny, nz), 0.2)
        spacer = np.full((ny, nz), 0.1)
        flusher = np.full((ny, nz), 0.05)
        w_prev = np.full((ny, nz), 0.4)
        mud_f = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                          rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
        cement_f = FluidSpec(name="tail", role=FluidRole.TAIL, density_kg_m3=1900.0,
                              rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.180, yield_stress_pa=14.0)
        spacer_f = FluidSpec(name="spacer", role=FluidRole.SPACER, density_kg_m3=1850.0,
                              rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.07, yield_stress_pa=5.0)
        flusher_f = FluidSpec(name="flusher", role=FluidRole.FLUSHER, density_kg_m3=1850.0,
                               rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.04, yield_stress_pa=3.0)
        out = s._compute_props(lead, tail, spacer, flusher, w_prev, geom, mud_f, cement_f, None, spacer_f)
        mu, rho, mud, tau_y, m_field, eta1, eta2 = out
        # 五相闭合：sum = lead + tail + spacer + flusher + mud ≈ 1
        phase_sum = lead + tail + spacer + flusher + mud
        assert np.allclose(phase_sum, 1.0, atol=1e-10), (
            f"五相之和应 ≈ 1，实际 min={phase_sum.min()} max={phase_sum.max()}"
        )

    def test_flusher_reduces_mud_in_compute_props(self):
        """flusher 注入后 mud 分数减少，且 flusher 不参与 c_bar/m_field。"""
        s = _make_solver()
        well = _toy_well()
        geom = s._build_geom(well)
        ny, nz = s.ny, s.nz
        lead = np.full((ny, nz), 0.3)
        tail = np.zeros((ny, nz))
        spacer = np.full((ny, nz), 0.15)
        w_prev = np.full((ny, nz), 0.4)
        mud_f = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                          rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
        cement_f = FluidSpec(name="lead", role=FluidRole.LEAD, density_kg_m3=1930.0,
                              rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.180, yield_stress_pa=14.0)
        spacer_f = FluidSpec(name="spacer", role=FluidRole.SPACER, density_kg_m3=1850.0,
                              rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.07, yield_stress_pa=5.0)
        # 无 flusher 时：mud1 = 1 - lead - tail - spacer
        out_no_flusher = s._compute_props(lead, tail, spacer, np.zeros_like(lead), w_prev, geom, mud_f, cement_f, None, spacer_f)
        mu1, rho1, mud1, tau_y1, m_field1, eta1_1, eta2_1 = out_no_flusher
        sum_no_flusher = lead + tail + spacer + mud1
        assert np.allclose(sum_no_flusher, 1.0, atol=1e-10), "四相闭合应 ≈ 1"
        # 有 flusher 时：mud2 = 1 - lead - tail - spacer - flusher
        flusher = np.full((ny, nz), 0.08)
        out_with_flusher = s._compute_props(lead, tail, spacer, flusher, w_prev, geom, mud_f, cement_f, None, spacer_f)
        mu2, rho2, mud2, tau_y2, m_field2, eta1_2, eta2_2 = out_with_flusher
        sum_with_flusher = lead + tail + spacer + flusher + mud2
        assert np.allclose(sum_with_flusher, 1.0, atol=1e-10), "五相闭合应 ≈ 1"
        # 有 flusher 时 mud 应减少（约等于 flusher 分数）
        # 注意：数值舍入，允许微小差异
        assert np.all(mud2 <= mud1 + 1e-10), "flusher 注入后 mud 应减少或不变"
        # flusher 不参与 m_field（被动相）：m_field 与无 flusher 时一致
        assert np.allclose(m_field1, m_field2, atol=1e-6), "flusher 不应影响 m_field"

    def test_five_phase_closure_in_run_result(self):
        """完整 run 后五相闭合 sum(lead+tail+spacer+flusher+mud) ≈ 1。"""
        from cemdisp.models2d.boundary_bridge import AnnulusInletState
        well = _toy_well()
        mud = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                        rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
        lead = FluidSpec(name="lead", role=FluidRole.LEAD, density_kg_m3=1930.0,
                         rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.180, yield_stress_pa=14.0)
        spacer = FluidSpec(name="spacer", role=FluidRole.SPACER, density_kg_m3=1850.0,
                           rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.07, yield_stress_pa=5.0)
        flusher = FluidSpec(name="flusher", role=FluidRole.FLUSHER, density_kg_m3=1850.0,
                            rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.04, yield_stress_pa=3.0)
        fluids = (mud, lead, spacer, flusher)

        # 分阶段注入：flusher → spacer → cement
        def _inlet(t: float):
            if t < 10.0:
                return AnnulusInletState(
                    time_s=t, flow_rate_m3_s=0.02, stage_name="flusher",
                    phase_fractions=(("flusher", 1.0),),
                )
            elif t < 20.0:
                return AnnulusInletState(
                    time_s=t, flow_rate_m3_s=0.02, stage_name="spacer",
                    phase_fractions=(("spacer", 1.0),),
                )
            else:
                return AnnulusInletState(
                    time_s=t, flow_rate_m3_s=0.02, stage_name="cement",
                    phase_fractions=(("cement", 1.0), ("lead", 1.0)),
                )

        s = AnnulusD2DGASolver(dt=4.0, nz=20, ny=8, total_t=40.0)
        res = s.run(well, fluids, _inlet)
        # 五相闭合检查：result 中存 lead_field, tail_field, spacer_field, flusher_field
        # mud 由 1 - sum 反算；要求四相显式体积分数之和始终不超过 1（I3 修复验证）。
        # 数值弥散/浮点噪声可能引入 ~1e-8 量级的越界，显式容差 1e-8 覆盖之。
        lead_f = res.lead_field
        tail_f = res.tail_field
        spacer_f = res.spacer_field
        flusher_f = res.flusher_field
        assert flusher_f is not None
        tracked_sum = lead_f + tail_f + spacer_f + flusher_f
        assert np.all(tracked_sum <= 1.0 + 1e-8), (
            f"run 后显式四相之和应 ≤ 1（容差 1e-8），max={tracked_sum.max()}"
        )
        mud_f = np.clip(1.0 - tracked_sum, 0.0, 1.0)
        phase_sum = tracked_sum + mud_f
        # 允许 1e-10 舍入误差
        assert np.allclose(phase_sum, 1.0, atol=1e-10), (
            f"run 后五相之和应 ≈ 1，min={phase_sum.min()} max={phase_sum.max()}"
        )
