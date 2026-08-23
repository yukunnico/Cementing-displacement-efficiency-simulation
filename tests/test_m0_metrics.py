"""M0 指标层 solver 级测试：η_N 进 summary、评价窗、失稳指数、低尾指标。"""
import numpy as np
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
from cemdisp.data.well_spec import WellSpec, DepthValuePoint, EvaluationWindow
from cemdisp.data.fluid_spec import FluidSpec, FluidRole, RheologyModel
from cemdisp.models2d.boundary_bridge import AnnulusInletState


def _solver(**kw):
    return AnnulusD2DGASolver(dt=4.0, nz=30, ny=12, total_t=40.0, **kw)


def _multi_window_well():
    pts = lambda d, v: DepthValuePoint(depth_md_m=d, value=v)
    return WellSpec(
        well_name="m0", top_md_m=1000.0, bottom_md_m=1100.0,
        shoe_md_m=1100.0, hanger_md_m=1000.0,
        casing_id_mm=200.0, liner_od_mm=139.7, liner_id_mm=108.0,
        hole_diameter_profile=[pts(1000.0, 215.9), pts(1100.0, 215.9)],
        inclination_profile=[pts(1000.0, 5.0), pts(1100.0, 5.0)],
        standoff_profile=[pts(1000.0, 0.83), pts(1100.0, 0.83)],
        evaluation_windows=[
            EvaluationWindow(name="cbl窗", top_md_m=1020.0, bottom_md_m=1080.0, window_type="cbl"),
            EvaluationWindow(name="model窗", top_md_m=1000.0, bottom_md_m=1100.0, window_type="model_focus"),
        ],
    )


def _fluids():
    mud = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                    rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
    lead = FluidSpec(name="lead", role=FluidRole.LEAD, density_kg_m3=1930.0,
                     rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.180, yield_stress_pa=14.0)
    return (mud, lead)


def _inlet(t: float):
    return AnnulusInletState(
        time_s=t, flow_rate_m3_s=0.02, stage_name="pump",
        phase_fractions=(("cement", 1.0), ("lead", 1.0)),
    )


def _run_minimal():
    """跑一个最小 solver；run() 实参以当前签名为准，参照 TestStaticWallLayer 装配。"""
    s = _solver()
    return s.run(_multi_window_well(), _fluids(), _inlet)


def test_eta_narrow_in_summary_and_in_range():
    from cemdisp.diagnostics.displacement_metrics import compute_displacement_metrics
    res = _run_minimal()
    fr = res.summary["最终结果"]
    assert "窄四分位效率" in fr
    assert 0.0 <= fr["窄四分位效率"] <= 1.0
    assert 0.0 <= res.summary["eta_narrow"] <= 1.0
    dm = compute_displacement_metrics(res)
    # dm 层按设计 round(eta_narrow, 6)（displacement_metrics.py），
    # 故 summary 全精度值与 dm 值允许 6 位小数舍入误差（< 5e-7）。
    assert abs(fr["窄四分位效率"] - dm.eta_narrow) < 1e-5
