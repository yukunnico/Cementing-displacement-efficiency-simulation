"""改进 D2DGA 三闭包（auto-m / I3 / 真体力）求解器级测试。"""
import numpy as np
import pytest
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
from cemdisp.data.fluid_spec import FluidSpec, FluidRole, RheologyModel
from cemdisp.data.well_spec import WellSpec, DepthValuePoint, EvaluationWindow


def _make_solver(**kw) -> AnnulusD2DGASolver:
    return AnnulusD2DGASolver(dt=4.0, nz=20, ny=10, total_t=40.0, **kw)


def _make_geom(solver, well_spec):
    return solver._build_geom(well_spec)


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