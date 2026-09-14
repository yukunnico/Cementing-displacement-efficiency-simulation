"""浮力口径统一模块契约测试。

守护 `cemdisp.models2d.buoyancy` 的规格（Task 3）：
1. 浮力数 b 的定义与符号（Z&F22 p.8）；
2. 密度倒置 → b<0；
3. 幂律/HB 表观黏度用 K·γ̇^(n−1)，**不得静默回退**到硬编码 0.05；
4. Froude 数平方的量级（Z&F22 (2.6)），不是硬编码 1.0；
5. **全仓唯一**顶替液密度口径的四态退化（M1）；
6. **半间隙 d̂ = (r_o−r_i)/2 = (井径−外径)/4 的语义钉死**（M1/C1），
   以及求解器 summary 段的调用点确实按该口径传参（C1 回归锚）。
"""

import numpy as np
import pytest

from cemdisp.data.fluid_spec import FluidSpec, FluidRole, RheologyModel
from cemdisp.data.well_spec import DepthValuePoint, EvaluationWindow, WellSpec
from cemdisp.models2d import buoyancy
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import AnnulusInletState
from cemdisp.models2d.buoyancy import (buoyancy_number, froude_squared,
                                       fluid_apparent_viscosity)


def test_buoyancy_number_formula_and_sign():
    b = buoyancy_number(2100.0, 2020.0, 0.0257, 0.066, 0.563)
    assert b == pytest.approx(80 * 9.81 * 0.0257 ** 2 / (0.066 * 0.563), rel=1e-12)
    assert b > 0


def test_density_inversion_is_negative():
    assert buoyancy_number(1900.0, 1960.0, 0.045, 0.058, 0.764) < 0


def test_power_law_viscosity_not_silent_fallback():
    pl = FluidSpec(name="泥浆", role=FluidRole.MUD, density_kg_m3=1960.0,
                   rheology_model=RheologyModel.POWER_LAW,
                   power_law_n=0.72, consistency_k=0.381)
    v = fluid_apparent_viscosity(pl, shear_rate=10.0)
    assert v == pytest.approx(0.381 * 10.0 ** (0.72 - 1.0), rel=1e-12)
    assert v != 0.05


def test_froude_squared_matches_manual_scale():
    f2 = froude_squared(mu_displaced=0.058, w0_mps=0.764, half_gap_m=0.0458,
                        rho_displaced=1200.0, gap_scale_m=0.0917,
                        mean_radius_m=0.1071)
    assert 1e-4 < f2 < 1.0     # 实测量级 ~4.7e-3，绝不是 1.0


# ---------------------------------------------------------------------------
# M1: displacing_density_kg_m3 —— 全仓唯一口径的四态退化
# ---------------------------------------------------------------------------
def _fluid(name: str, role: FluidRole, rho: float) -> FluidSpec:
    return FluidSpec(name=name, role=role, density_kg_m3=rho,
                     plastic_viscosity_pa_s=0.05)


def test_displacing_density_four_states():
    """(领浆,尾浆) 都在 → 0.67/0.33 加权；单在 → 取该相；都无 → 取泥浆。"""
    lead = _fluid("领浆", FluidRole.LEAD, 2100.0)
    tail = _fluid("尾浆", FluidRole.TAIL, 1900.0)
    mud = _fluid("泥浆", FluidRole.MUD, 1960.0)

    assert buoyancy.displacing_density_kg_m3(lead, tail, mud) == pytest.approx(
        0.67 * 2100.0 + 0.33 * 1900.0, rel=1e-12)
    assert buoyancy.displacing_density_kg_m3(lead, None, mud) == pytest.approx(2100.0, rel=1e-12)
    assert buoyancy.displacing_density_kg_m3(None, tail, mud) == pytest.approx(1900.0, rel=1e-12)
    assert buoyancy.displacing_density_kg_m3(None, None, mud) == pytest.approx(1960.0, rel=1e-12)


# ---------------------------------------------------------------------------
# M1/C1: 半间隙语义钉死 —— d̂ = (r_o−r_i)/2 = (井径−外径)/4
# ---------------------------------------------------------------------------
_TOY_HOLE_MM = 215.9
_TOY_OD_MM = 139.7


def _toy_well() -> WellSpec:
    pts = lambda d, v: DepthValuePoint(depth_md_m=d, value=v)  # noqa: E731
    return WellSpec(
        well_name="toy",
        top_md_m=1000.0, bottom_md_m=1100.0, shoe_md_m=1100.0, hanger_md_m=1000.0,
        casing_id_mm=200.0, liner_od_mm=_TOY_OD_MM, liner_id_mm=108.0,
        hole_diameter_profile=[pts(1000.0, _TOY_HOLE_MM), pts(1100.0, _TOY_HOLE_MM)],
        inclination_profile=[pts(1000.0, 5.0), pts(1100.0, 5.0)],
        standoff_profile=[pts(1000.0, 0.83), pts(1100.0, 0.83)],
        evaluation_windows=[EvaluationWindow(name="w", top_md_m=1000.0,
                                             bottom_md_m=1100.0, window_type="full")],
    )


def test_geom_H_semantics_is_quarter_of_diametral_clearance():
    """钉死半间隙语义：geom["H"] 是 Z&F22 的 d̂ = (r_o−r_i)/2 = (井径−外径)/4，
    输运用的 geom["b"] 是它的两倍（b=2H 逐格成立，即全环隙 2d̂）。

    这条是 C1 的语义锚：把 d̂ 误当 (井径−外径)/2（即 2d̂）会让 b 放大 4×。
    """
    solver = AnnulusD2DGASolver(dt=4.0, nz=20, ny=8, total_t=40.0)
    geom = solver._build_geom(_toy_well())
    d_half_paper = (_TOY_HOLE_MM - _TOY_OD_MM) / 4.0 / 1000.0   # (r_o−r_i)/2，单位 m
    assert d_half_paper == pytest.approx(0.01905, rel=1e-12)

    mean_h = float(np.mean(geom["H"]))
    mean_b = float(np.mean(geom["b"]))
    assert mean_h == pytest.approx(d_half_paper, rel=0.02)      # 差值为体积 scale 残差
    assert mean_b == pytest.approx(2.0 * mean_h, rel=1e-12)     # b 恒为 2H（逐格）
    # 反证：2d̂ 口径（C1 的错误值）会把 b 放大 4×，必须被这条测试挡住
    assert mean_h != pytest.approx(2.0 * d_half_paper, rel=0.02)


def test_summary_call_site_uses_unique_density_and_paper_half_gap(monkeypatch):
    """求解器 summary 段的调用点口径锚（C1 回归锚）。

    用 spy 捕获 `buoyancy.buoyancy_number` 的实参，断言：
    顶替液密度走唯一口径、半间隙是 geom["H"]（= d̂，而非 2d̂）、w₀ = q/A。
    """
    captured: dict[str, float] = {}
    real = buoyancy.buoyancy_number

    def _spy(*args, **kwargs):
        captured.update(kwargs)
        return real(*args, **kwargs)

    monkeypatch.setattr(buoyancy, "buoyancy_number", _spy)

    well = _toy_well()
    mud = _fluid("泥浆", FluidRole.MUD, 1960.0)
    tail = _fluid("尾浆", FluidRole.TAIL, 1900.0)
    fluids = (mud, tail)
    q_m3s = 0.02

    def _inlet(t: float) -> AnnulusInletState:
        return AnnulusInletState(time_s=t, flow_rate_m3_s=q_m3s, stage_name="pump",
                                 phase_fractions=(("cement", 1.0), ("tail", 1.0)))

    solver = AnnulusD2DGASolver(dt=4.0, nz=20, ny=8, total_t=40.0)
    result = solver.run(well, fluids, _inlet)

    geom = solver._build_geom(well)
    assert captured, "summary 段未调用 buoyancy.buoyancy_number"
    # ① 顶替液密度：唯一口径
    assert captured["rho_displacing"] == pytest.approx(
        buoyancy.displacing_density_kg_m3(None, tail, mud), rel=1e-12)
    # ② 半间隙：geom["H"]（= d̂），且 ≈ (井径−外径)/4（而非其 2 倍）
    assert captured["half_gap_m"] == pytest.approx(float(np.mean(geom["H"])), rel=1e-12)
    assert captured["half_gap_m"] == pytest.approx(
        (_TOY_HOLE_MM - _TOY_OD_MM) / 4.0 / 1000.0, rel=0.02)
    # ③ 被顶替液黏度：Bingham 走 PV 而非静默回退
    assert captured["mu_displaced"] == pytest.approx(mud.plastic_viscosity_pa_s, rel=1e-12)
    # ④ w₀：截面平均速度 q/A
    area = float(np.pi / 4.0 * ((_TOY_HOLE_MM / 1000.0) ** 2 - (_TOY_OD_MM / 1000.0) ** 2))
    assert captured["w0_mps"] == pytest.approx(q_m3s / area, rel=1e-12)
    # summary 的值即该调用返回值
    assert float(result.summary["buoyancy_number"]) == pytest.approx(
        real(**captured), rel=1e-12)


# ---------------------------------------------------------------------------
# Task 4: _buoyancy_force_vector 的 F² 标定（Z&F22 (2.6)/(2.5b)）
# ---------------------------------------------------------------------------
def test_buoyancy_force_vector_uses_froude_scale():
    """f_φ 必须带 1/F² 标定。f2 的取值按 (2.6) 用呼101 实参**现算**，
    不要硬抄——`F² = τ̂₀/(ρ̂₁ĝδ₀r̂ₐ*)`，`τ̂₀ = μ̂₁ŵ₀/d̂`，`d̂ = mean(geom["H"])`。
    参考量级 O(1e-2)。"""
    import numpy as np
    from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
    s = AnnulusD2DGASolver(ny=40, nz=2)
    geom = {"phi": np.linspace(0, 1, 40), "hole_mm": np.full((1, 2), 260.0),
            "od_mm": np.full((1, 2), 168.3)}
    f_unit, _ = s._buoyancy_force_vector(geom, 1.9, f2=1.0)
    f_phys, _ = s._buoyancy_force_vector(geom, 1.9, f2=1.0e-2)
    # 断言**正比关系**（f ∝ 1/F²），而非任意阈值——物理 F² 取 O(10⁻²)，
    # 固定"50×"阈值会随 F² 取值失效。
    assert f_phys.max() == pytest.approx(f_unit.max() / 1.0e-2, rel=1e-12)


def test_solver_call_sites_pass_physical_f2(monkeypatch):
    """两个调用点必须传入按 (2.6) 现算的 F²（量级 O(10⁻²)），不得硬编码 1.0。

    覆盖：`_compute_velocity` 的 R3 真体力段（式 2.5b）与 run 循环的
    R2 I3 浮力弥散通量段（式 4.25 第二项）——两处共用同一浮力向量。
    """
    captured: list[float] = []
    real = AnnulusD2DGASolver._buoyancy_force_vector

    def _spy(self, geom, beta_deg, f2):
        captured.append(float(f2))
        return real(self, geom, beta_deg, f2)

    monkeypatch.setattr(AnnulusD2DGASolver, "_buoyancy_force_vector", _spy)
    well = _toy_well()
    fluids = (_fluid("泥浆", FluidRole.MUD, 1960.0), _fluid("尾浆", FluidRole.TAIL, 1900.0))
    q_m3s = 0.02

    def _inlet(t: float) -> AnnulusInletState:
        return AnnulusInletState(time_s=t, flow_rate_m3_s=q_m3s, stage_name="pump",
                                 phase_fractions=(("cement", 1.0), ("tail", 1.0)))

    solver = AnnulusD2DGASolver(dt=4.0, nz=20, ny=8, total_t=40.0)
    solver.run(well, fluids, _inlet)

    assert captured, "两个调用点都未调用 _buoyancy_force_vector"
    for f2 in captured:
        assert f2 != 1.0
        assert 1e-4 < f2 < 1.0
