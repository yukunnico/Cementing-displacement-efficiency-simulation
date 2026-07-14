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