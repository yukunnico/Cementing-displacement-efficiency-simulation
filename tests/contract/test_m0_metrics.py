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


def test_instability_index_linear_and_log():
    res = _run_minimal()
    fr = res.summary["最终结果"]
    assert "最终失稳指数_线性" in fr and "最终失稳指数_对数" in fr
    proxy = float(res.metrics["instability_proxy"].iloc[-1])
    assert abs(fr["最终失稳指数_线性"] - proxy) < 1e-9
    assert abs(fr["最终失稳指数_对数"] - np.log10(1.0 + proxy)) < 1e-9


def test_evaluation_window_efficiencies():
    res = _run_minimal()
    we = res.summary["评价窗效率"]
    assert {"cbl窗", "model窗"} <= set(we.keys())
    for name, d in we.items():
        assert d["window_type"] in ("cbl", "model_focus")
        assert 0.0 <= d["eta_E"] <= 1.0
        assert 0.0 <= d["eta_N"] <= 1.0


def test_low_tail_indicators():
    res = _run_minimal()
    lt = res.summary["低尾指标"]
    assert 0.0 <= lt["standoff低于0.5段占比"] <= 1.0
    assert 0.0 <= lt["窄边效率低于0.05域占比"] <= 1.0


def test_m0_does_not_change_existing_keys():
    res = _run_minimal()
    fr = res.summary["最终结果"]
    for k in ["全井段最终有效顶替效率", "最终水泥浆占据率", "最终窜槽指数", "最终混浆指数", "最终失稳指数"]:
        assert k in fr
    assert "effective_efficiency" in res.summary
