"""领先阶轴向浮力数 b 接进动力学的契约测试（Task 5，Z&F22 (4.14)/(4.22)）。

守护：
1. b_num>0（重顶替轻）应抬高窄边流动度——(4.22) 分层浮力项 ``−Δρ·I₂/(H·I₁)·f_φ``
   在本项目 Δρ≡ρ̂₂−ρ̂₁ 口径下为正，且 (4.14) 的逐通道半间隙分母 H(φ) 给出窄边放大；
2. b_num=0 时 pref 与**不含浮力项的基础流动度**逐位一致（重构锚：
   轴向浮力项关闭 ⇔ Task 4 行为，测试独立复算基础项而非同参自比）；
3. R26：`_mobility_profile` 的输出必须随 `f2` 变化（方位 f_φ = r_a·sinπφ·sinβ/F²
   通路）——防 Task 4 的 F² 定标在抽函数时二次静默消失；
4. R27：现场浮力数域（八井 b∈[−2.9, +20.7]）内 buoyancy_shape 恒正
   （K_AXIAL = 1/𝒢 推导的自洽下界，禁止 clip）。
"""

import numpy as np

from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver


def _geom(ny=40, nz=6):
    """合成几何（与 brief 一致）：H 随方位 ±30% 变化，b=2H 逐格成立。"""
    phi = np.linspace(0, 1, ny)
    H = np.broadcast_to(0.0458 * (1 + 0.3 * np.cos(np.pi * phi))[:, None], (ny, nz)).copy()
    return {"phi": phi, "H": H, "b": 2 * H, "y": np.linspace(0, np.pi * 0.1071, ny),
            "hole_mm": np.full((1, nz), 260.0), "od_mm": np.full((1, nz), 168.3),
            "s": np.linspace(0, 10.0, nz)}


def test_axial_buoyancy_raises_narrow_side_mobility():
    """b>0（重顶替轻）应抬高窄边流动度份额（(4.22) 分层项 + 1/H̃ 窄边放大）。"""
    s = AnnulusD2DGASolver(ny=40, nz=6)
    g = _geom()
    c = np.full_like(g["H"], 0.5)
    p0 = s._mobility_profile(c, b_num=0.0, geom=g, i1_base=0.44, m_local=0.5,
                             beta_deg=0.0, f2=1.0)
    p1 = s._mobility_profile(c, b_num=50.0, geom=g, i1_base=0.44, m_local=0.5,
                             beta_deg=0.0, f2=1.0)
    ny = g["H"].shape[0]
    narrow = slice(3 * ny // 4, ny)
    assert p1[narrow].mean() > p0[narrow].mean()


def test_zero_buoyancy_is_identity():
    """b_num=0 时 pref 与不含浮力项的基础流动度逐位一致（有判别力版）。

    基础流动度独立复算（默认 eta1=eta2=1、n_mix=1 ⇒ 幂律缝隙律退化牛顿：
    base = (b/b̄)²·I₁），不做同参数自比——任何 b_num=0 下残留的浮力项都会红。
    """
    s = AnnulusD2DGASolver(ny=40, nz=6)
    g = _geom()
    c = np.full_like(g["H"], 0.5)
    pref = s._mobility_profile(c, b_num=0.0, geom=g, i1_base=0.44, m_local=0.5,
                               beta_deg=0.0, f2=1.0)
    b_mean = np.mean(g["b"], axis=0, keepdims=True)
    base_expected = (g["b"] / np.maximum(b_mean, 1.0e-12)) ** 2 * 0.44
    assert np.array_equal(pref, base_expected)


def test_axial_buoyancy_is_additive_in_b_num():
    """轴向浮力项随 b_num 线性：pref(b)−pref(0) ∝ b（(4.14) 项的线性结构）。"""
    s = AnnulusD2DGASolver(ny=40, nz=6)
    g = _geom()
    c = np.full_like(g["H"], 0.5)
    kw = dict(geom=g, i1_base=0.44, m_local=0.5, beta_deg=0.0, f2=1.0)
    p0 = s._mobility_profile(c, b_num=0.0, **kw)
    p1 = s._mobility_profile(c, b_num=10.0, **kw)
    p2 = s._mobility_profile(c, b_num=30.0, **kw)
    delta1 = p1 - p0
    delta2 = p2 - p0
    assert np.allclose(delta2, 3.0 * delta1, rtol=1e-12)


def test_mobility_profile_output_varies_with_f2():
    """R26：输出必须随 f2 变化——方位 f_φ = r_a·sinπφ·sinβ/F² 通路（β>0、Δρ≠0）。"""
    s = AnnulusD2DGASolver(ny=40, nz=6)
    g = _geom()
    c = np.full_like(g["H"], 0.5)
    kw = dict(b_num=10.0, geom=g, i1_base=0.44, m_local=0.5, beta_deg=5.0,
              delta_rho=0.15)
    p_unit = s._mobility_profile(c, f2=1.0, **kw)
    p_phys = s._mobility_profile(c, f2=1.0e-2, **kw)
    assert not np.allclose(p_unit, p_phys)


def test_buoyancy_shape_stays_positive_over_field_b_range():
    """R27：现场浮力数域内（含密度倒置 hu1 b=−2.85）buoyancy_shape 恒正，禁止 clip。"""
    s = AnnulusD2DGASolver(ny=40, nz=6)
    g = _geom()
    # 八井实测 b 域端点 + 密度倒置井（Task 3 模型实测：hu1 −2.85 / hu2 −0.94）
    for b_num in (-2.9, -0.94, 0.0, 1.96, 9.06, 18.89, 20.7):
        pref = s._mobility_profile(np.full_like(g["H"], 0.5), b_num=b_num, geom=g,
                                   i1_base=0.44, m_local=0.5, beta_deg=0.0, f2=1.0)
        assert pref.min() > 0.0, f"b_num={b_num} 时 pref 出现非正值（min={pref.min()}）"


def test_compute_velocity_passthrough_uses_buoyancy_number(monkeypatch):
    """_compute_velocity 必须以 summary 同口径把 b_num 传入 _mobility_profile（接线锚）。

    T9（2026-09-15）钉旧代数路径（enable_stream_function=False）：b_num→
    `_mobility_profile` 是旧路径的 (4.14) 分层浮力接线；新路径（默认）浮力经
    (4.22) b 向量完整进入、不消费 b_num（双重计入排查见
    `_velocity_stream_function` docstring），其接线由
    tests/contract/test_stream_function_solver_integration.py 锁定。
    断言未改动。

    spy `buoyancy.buoyancy_number` 并检查**调用栈链**：Task 5 前它只在 run() 末尾的
    summary 块被调用；接线后动力学时间步（`_compute_velocity` → `_buoyancy_number_at`）
    也必须调用（浮力数进动力学而非只进 summary）。按调用栈断言，避免"只数次数"对步数敏感。
    """
    import inspect

    from cemdisp.models2d import buoyancy as buoyancy_mod

    from_dynamics: list[bool] = []
    real = buoyancy_mod.buoyancy_number

    def _spy(*args, **kwargs):
        frame = inspect.currentframe().f_back
        chain: list[str] = []
        f = frame
        while f is not None:
            chain.append(f.f_code.co_name)
            f = f.f_back
        from_dynamics.append("_compute_velocity" in chain)
        return real(*args, **kwargs)

    monkeypatch.setattr(buoyancy_mod, "buoyancy_number", _spy)

    # 合成几何 + 幂律泥浆/尾浆的最小 run（不需要现场数据）
    from cemdisp.data.fluid_spec import FluidSpec, FluidRole
    from cemdisp.models2d.boundary_bridge import AnnulusInletState

    mud = FluidSpec(name="泥浆", role=FluidRole.MUD, density_kg_m3=1960.0,
                    plastic_viscosity_pa_s=0.05)
    tail = FluidSpec(name="尾浆", role=FluidRole.TAIL, density_kg_m3=1900.0,
                     plastic_viscosity_pa_s=0.05)
    well = _toy_well()
    q_m3s = 0.02

    def _inlet(t: float) -> AnnulusInletState:
        return AnnulusInletState(time_s=t, flow_rate_m3_s=q_m3s, stage_name="pump",
                                 phase_fractions=(("cement", 1.0), ("tail", 1.0)))

    solver = AnnulusD2DGASolver(dt=4.0, nz=20, ny=8, total_t=40.0,
                                enable_stream_function=False)
    solver.run(well, (mud, tail), _inlet)

    assert any(from_dynamics), (
        "动力学时间步未调用 buoyancy.buoyancy_number（b_num 未接进 _compute_velocity）")
    assert not all(from_dynamics), "summary 段的 buoyancy_number 调用丢失（应两条路径都有）"


def _toy_well():
    """最小井规格（与 tests/contract/test_buoyancy.py 的 toy 井同口径）。"""
    from cemdisp.data.well_spec import DepthValuePoint, EvaluationWindow, WellSpec

    pts = lambda d, v: DepthValuePoint(depth_md_m=d, value=v)  # noqa: E731
    return WellSpec(
        well_name="toy",
        top_md_m=1000.0, bottom_md_m=1100.0, shoe_md_m=1100.0, hanger_md_m=1000.0,
        casing_id_mm=200.0, liner_od_mm=139.7, liner_id_mm=108.0,
        hole_diameter_profile=[pts(1000.0, 215.9), pts(1100.0, 215.9)],
        inclination_profile=[pts(1000.0, 5.0), pts(1100.0, 5.0)],
        standoff_profile=[pts(1000.0, 0.83), pts(1100.0, 0.83)],
        evaluation_windows=[EvaluationWindow(name="w", top_md_m=1000.0,
                                             bottom_md_m=1100.0, window_type="full")],
    )
