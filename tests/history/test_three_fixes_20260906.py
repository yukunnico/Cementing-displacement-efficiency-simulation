# STATUS: history —— 本文件锁定“重构前”的历史行为契约。
# 源模型口径重构（2026-09-14）后，相关断言预期变红；红灯不等同回归失败，
# 需按 docs/superpowers/plans/2026-09-14-d2dga-source-fidelity-and-repo-cleanup.md 逐条复核。
"""2026-09-06 三项修复的专项测试。

1. WASH/SPACER 选相修复：多种 WASH/SPACER 并存时用体积加权等效流体，
   不再取第一个（旧口径使真实隔离液物性在 2D 闭包中失效）。
2. e_clip 上限裁定：实测居中度井（standoff_measured=True）放开到 e_clip_measured_max，
   设计/代理值井维持保守默认；enable_e_clip_ruling=False 退回旧口径。
3. 幂律缝隙律：速度场指数 1+1/n（牛顿 n=1 退化为 b² 逐位兼容旧口径；
   剪切变稀 n<1 窄边分流比 b³ 更极端）；关闭开关回退 b²。
"""
import numpy as np
import pytest

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.well_spec import DepthValuePoint, WellSpec
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver


def _make_solver(**kw) -> AnnulusD2DGASolver:
    return AnnulusD2DGASolver(dt=4.0, nz=20, ny=10, total_t=40.0, **kw)


def _toy_well(standoff_measured=False, standoff=0.45):
    pts = lambda d, v: DepthValuePoint(depth_md_m=d, value=v)
    return WellSpec(
        well_name="toy",
        top_md_m=1000.0, bottom_md_m=1100.0, shoe_md_m=1100.0, hanger_md_m=1000.0,
        casing_id_mm=200.0, liner_od_mm=139.7, liner_id_mm=108.0,
        hole_diameter_profile=[pts(1000.0, 215.9), pts(1100.0, 215.9)],
        inclination_profile=[pts(1000.0, 5.0), pts(1100.0, 5.0)],
        standoff_profile=[pts(1000.0, standoff), pts(1100.0, standoff)],
        standoff_measured=standoff_measured,
        evaluation_windows=[],
    )


# ---------------------------------------------------------------------------
# 1. WASH/SPACER 选相修复
# ---------------------------------------------------------------------------
class TestCompositeSpacer:
    def test_single_spacer_passthrough(self):
        """单一 WASH/SPACER 时直接透传（不合成）。"""
        f = FluidSpec("隔离液", FluidRole.SPACER, 2050.0, RheologyModel.BINGHAM,
                      plastic_viscosity_pa_s=0.03, yield_stress_pa=8.0)
        mud = FluidSpec("泥浆", FluidRole.MUD, 1900.0, RheologyModel.BINGHAM,
                        plastic_viscosity_pa_s=0.05, yield_stress_pa=9.0)
        lead = FluidSpec("领浆", FluidRole.LEAD, 2100.0, RheologyModel.POWER_LAW,
                         power_law_n=0.8, consistency_k=0.4)
        s = _make_solver()
        picked = s._pick_fluids((mud, f, lead))
        assert picked[3] is f  # 同一对象透传

    def test_composite_bingham_weighted(self):
        """两种 Bingham 合成：PV/YP/密度均线性加权。"""
        a = FluidSpec("平衡液", FluidRole.WASH, 1750.0, RheologyModel.BINGHAM,
                      plastic_viscosity_pa_s=0.03, yield_stress_pa=3.0)
        b = FluidSpec("隔离液", FluidRole.SPACER, 2050.0, RheologyModel.BINGHAM,
                      plastic_viscosity_pa_s=0.05, yield_stress_pa=9.0)
        comp = AnnulusD2DGASolver._composite_spacer_fluid([a, b], [40.0, 20.0])
        assert comp.density_kg_m3 == pytest.approx((1750.0 * 2 + 2050.0) / 3.0)
        assert comp.plastic_viscosity_pa_s == pytest.approx((0.03 * 2 + 0.05) / 3.0)
        assert comp.yield_stress_pa == pytest.approx((3.0 * 2 + 9.0) / 3.0)
        assert comp.rheology_model == RheologyModel.BINGHAM

    def test_composite_power_law_log_k(self):
        """两种幂律合成：n 线性、K 对数加权（与 n_mix/kappa_mix 口径一致）。"""
        a = FluidSpec("隔离液1", FluidRole.SPACER, 2000.0, RheologyModel.POWER_LAW,
                      power_law_n=0.6, consistency_k=1.0)
        b = FluidSpec("隔离液2", FluidRole.SPACER, 1950.0, RheologyModel.POWER_LAW,
                      power_law_n=0.8, consistency_k=2.0)
        comp = AnnulusD2DGASolver._composite_spacer_fluid([a, b], [1.0, 1.0])
        assert comp.rheology_model == RheologyModel.POWER_LAW
        assert comp.power_law_n == pytest.approx(0.7)
        assert comp.consistency_k == pytest.approx(np.sqrt(1.0 * 2.0))

    def test_composite_mixed_models_bingham_fallback(self):
        """Bingham+幂律混合：统一映射 Bingham（幂律折算 γ_ref=20 等效黏度）。"""
        a = FluidSpec("平衡液", FluidRole.WASH, 1850.0, RheologyModel.BINGHAM,
                      plastic_viscosity_pa_s=0.041, yield_stress_pa=9.2)
        b = FluidSpec("驱油隔离液", FluidRole.SPACER, 2000.0, RheologyModel.POWER_LAW,
                      power_law_n=0.54, consistency_k=2.12)
        comp = AnnulusD2DGASolver._composite_spacer_fluid([a, b], [25.0, 25.0])
        assert comp.rheology_model == RheologyModel.BINGHAM
        eff_k = 2.12 * 20.0 ** (0.54 - 1.0)
        assert comp.plastic_viscosity_pa_s == pytest.approx((0.041 + eff_k) / 2.0)
        # 密度仍是严格体积加权
        assert comp.density_kg_m3 == pytest.approx(1925.0)

    def test_volume_weights_from_schedule(self):
        """_wash_spacer_volume_weights：按泵注体积取权，泵注名不匹配时权重 0。"""
        wash = FluidSpec("平衡液", FluidRole.WASH, 1850.0, RheologyModel.BINGHAM,
                         plastic_viscosity_pa_s=0.04, yield_stress_pa=9.0)
        spacer = FluidSpec("隔离液", FluidRole.SPACER, 2000.0, RheologyModel.BINGHAM,
                           plastic_viscosity_pa_s=0.03, yield_stress_pa=5.0)
        steps = [
            type("Step", (), {"fluid_name": "平衡液", "volume_m3": 25.0, "rate_m3_min": 1.5})(),
            type("Step", (), {"fluid_name": "驱油隔离液", "volume_m3": 25.0, "rate_m3_min": 1.5})(),
        ]
        sched = type("S", (), {"steps": steps})()
        weights = AnnulusD2DGASolver._wash_spacer_volume_weights([wash, spacer], sched)
        assert weights == [25.0, 0.0]

    def test_pick_fluids_composite_name_and_density(self):
        """_pick_fluids：两种 WASH/SPACER 并存时合成，名称与密度可观测。"""
        mud = FluidSpec("泥浆", FluidRole.MUD, 1900.0, RheologyModel.BINGHAM,
                        plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
        wash = FluidSpec("平衡液", FluidRole.WASH, 1750.0, RheologyModel.BINGHAM,
                         plastic_viscosity_pa_s=0.03, yield_stress_pa=3.0)
        spacer_f = FluidSpec("隔离液", FluidRole.SPACER, 2050.0, RheologyModel.BINGHAM,
                             plastic_viscosity_pa_s=0.05, yield_stress_pa=9.0)
        tail = FluidSpec("尾浆", FluidRole.TAIL, 1900.0, RheologyModel.BINGHAM,
                         plastic_viscosity_pa_s=0.18, yield_stress_pa=14.0)
        s = _make_solver()
        mud_f, lead_f, tail_f, sp_f, _ = s._pick_fluids((mud, wash, spacer_f, tail))
        assert sp_f.name == "平衡液+隔离液"
        assert sp_f.density_kg_m3 == pytest.approx((1750.0 + 2050.0) / 2.0)


# ---------------------------------------------------------------------------
# 2. e_clip 上限裁定
# ---------------------------------------------------------------------------
class TestEClipRuling:
    def test_assumed_well_keeps_default_cap(self):
        """设计值井（standoff_measured=False）：e 上限保持 0.55。

        STATUS（Task 10 预期红，2026-09-15）：e_clip 硬截断已移除（Pelipenko04
        (2.1) e∈[0,1) 文献口径），e=1−standoff 直取 → 本井 e=0.95 ≠ 0.55。
        截断语义由 tests/contract/test_eccentricity_no_clip.py 锁定。断言未改动。
        """
        s = _make_solver()
        well = _toy_well(standoff_measured=False, standoff=0.05)  # e 期望 0.55
        geom = s._build_geom(well)
        assert geom["e"].max() == pytest.approx(0.55)

    def test_measured_well_raises_to_090(self):
        """实测井（standoff_measured=True）：e 上限放开到 0.90。

        STATUS（Task 10 预期红，2026-09-15）：截断移除后 e=0.95 ≠ 0.90；
        2026-09-06 "实测井放开 0.90" 裁定随截断一同退役（文献口径无按数据
        来源选上限概念）。断言未改动。
        """
        s = _make_solver()
        well = _toy_well(standoff_measured=True, standoff=0.05)
        geom = s._build_geom(well)
        assert geom["e"].max() == pytest.approx(0.90)

    def test_ruling_off_falls_back(self):
        """enable_e_clip_ruling=False：实测井也维持 0.55（旧口径）。

        STATUS（Task 10 预期红，2026-09-15）：截断移除后 e=0.95 ≠ 0.55；
        enable_e_clip_ruling 形参弃用（偏离 legacy 默认传值仅触发
        DeprecationWarning）。断言未改动。
        """
        s = _make_solver(enable_e_clip_ruling=False)
        well = _toy_well(standoff_measured=True, standoff=0.05)
        geom = s._build_geom(well)
        assert geom["e"].max() == pytest.approx(0.55)

    def test_explicit_override_respected_for_assumed(self):
        """显式 e_clip_max=0.90 + 设计值井：上限用显式值。

        STATUS（Task 10 预期红，2026-09-15）：截断移除后显式 e_clip_max 不再
        生效，e=0.95 ≠ 0.90。断言未改动。
        """
        s = _make_solver(e_clip_max=0.90)
        well = _toy_well(standoff_measured=False, standoff=0.05)
        geom = s._build_geom(well)
        assert geom["e"].max() == pytest.approx(0.90)

    def test_hu101_assumed_still_flagged_false(self):
        """hu101 名义口径（measured_standoff=None）standoff_measured=False。"""
        from cemdisp.data.loaders import load_hu101_tailpipe
        well, _f, _s, _v = load_hu101_tailpipe()
        assert well.standoff_measured is False

    def test_hu101_measured_flagged_true(self):
        """hu101 实测口径（measured_standoff='between_centralizers'）标记 True。"""
        from cemdisp.data.loaders import load_hu101_tailpipe
        well, _f, _s, _v = load_hu101_tailpipe(measured_standoff="between_centralizers")
        assert well.standoff_measured is True
        assert well.standoff_profile[0].value == pytest.approx(0.78)


# ---------------------------------------------------------------------------
# 3. 幂律缝隙律
# ---------------------------------------------------------------------------
class TestPowerLawGapLaw:
    @staticmethod
    def _fields_with_lead(solver, well, n_lead):
        """构造均匀水泥场（lead=0.6），返回所需输入与 n_mix。"""
        geom = solver._build_geom(well)
        ny, nz = solver.ny, solver.nz
        lead = np.full((ny, nz), 0.6)
        tail = np.zeros((ny, nz))
        w_prev = np.full((ny, nz), 0.4)
        mud_f = FluidSpec("mud", FluidRole.MUD, 1900.0, RheologyModel.BINGHAM,
                          plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
        lead_f = FluidSpec("lead", FluidRole.LEAD, 2100.0, RheologyModel.POWER_LAW,
                           power_law_n=n_lead, consistency_k=0.4)
        _, _, _, _, _, _, _, n_mix, _ = solver._compute_props(
            lead, tail, np.zeros_like(lead), w_prev, geom,
            mud_f, lead_f, None, None)
        return geom, lead, tail, w_prev, mud_f, lead_f, n_mix

    def test_shear_thinning_steepens_profile(self):
        """剪切变稀（n_mix<1）：窄边/宽边速度比应低于旧口径（b²）。

        STATUS（Task 9 预期红，2026-09-15）：幂律缝隙律 (b/b̄)^(1+1/n) 是旧代数
        路径的 base 构造（`_mobility_base`，enable_power_law_gap_law 消费）；
        T9 新路径（默认）为 Z&F22 (4.21) 牛顿两层闭包——流变经标量表观黏度
        η₁/η₂/m 进入，幂律指数不进椭圆算子 ⇒ 两侧都走新路径时 ratio_new 与
        ratio_old 相同，断言 ratio_new < ratio_old 不再成立。旧路径行为由
        enable_stream_function=False 保留（幂律缝隙律语义不变）。断言未改动。
        """
        well = _toy_well(standoff=0.45)
        s_new = _make_solver()
        s_old = _make_solver(enable_power_law_gap_law=False)
        geom, lead, tail, w_prev, mud_f, lead_f, n_mix = self._fields_with_lead(s_new, well, 0.54)
        assert np.mean(n_mix) < 1.0
        w_new, *_ = s_new._compute_velocity(lead, tail, np.zeros_like(lead),
                                            geom, q_m3s=0.02, w_prev=w_prev,
                                            mud_fluid=mud_f, lead_fluid=lead_f,
                                            tail_fluid=None, spacer_fluid=None)
        w_old, *_ = s_old._compute_velocity(lead, tail, np.zeros_like(lead),
                                            geom, q_m3s=0.02, w_prev=w_prev,
                                            mud_fluid=mud_f, lead_fluid=lead_f,
                                            tail_fluid=None, spacer_fluid=None)
        ratio_new = w_new[-1, :].mean() / w_new[0, :].mean()
        ratio_old = w_old[-1, :].mean() / w_old[0, :].mean()
        assert ratio_new < ratio_old

    def test_flux_conservation_preserved(self):
        """幂律缝隙律下每深度截面通量仍守恒（B1 口径：每列 Σw·b·dy = q_half）。

        STATUS（Task 9 预期红，2026-09-15）：旧路径按构造使矩形和
        Σw·b·dy = q_half 逐位成立；T9 新路径的列通量在**梯形求积**下精确
        = q_half（差分-梯形恒等式，与模型体积核算 _trapez2d 同口径），矩形和
        与梯形积差 O(边界半权) —— 本测试用 np.sum（矩形）+ rtol=1e-9 断言在
        新路径下不成立。新路径守恒契约由
        tests/contract/test_stream_function_solver_integration.py::
        TestStreamFunctionPathWiring::test_new_path_column_flux_conservation
        锁定。断言未改动。
        """
        well = _toy_well()
        s_new = _make_solver()
        geom, lead, tail, w_prev, mud_f, lead_f, _ = self._fields_with_lead(s_new, well, 0.54)
        w, *_ = s_new._compute_velocity(lead, tail, np.zeros_like(lead),
                                        geom, q_m3s=0.02, w_prev=w_prev,
                                        mud_fluid=mud_f, lead_fluid=lead_f,
                                        tail_fluid=None, spacer_fluid=None)
        dy = geom["y"][1] - geom["y"][0]
        col_flux = np.sum(w * geom["b"], axis=0) * dy  # 每深度截面半环空通量
        assert np.allclose(col_flux, 0.01, rtol=1e-9)  # = q_half = Q/2

    def test_pure_mud_uniform_n1_bitwise_matches_old(self):
        """全场纯泥浆（无水泥，n_mix≡1）：两口径速度场逐位一致。"""
        well = _toy_well()
        s_new = _make_solver()
        s_old = _make_solver(enable_power_law_gap_law=False)
        geom = s_new._build_geom(well)
        ny, nz = s_new.ny, s_new.nz
        lead = np.zeros((ny, nz))
        w_prev = np.full((ny, nz), 0.4)
        mud_f = FluidSpec("mud", FluidRole.MUD, 1900.0, RheologyModel.BINGHAM,
                          plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
        w_new, *_ = s_new._compute_velocity(lead, lead, lead,
                                            geom, q_m3s=0.02, w_prev=w_prev,
                                            mud_fluid=mud_f, lead_fluid=None,
                                            tail_fluid=None, spacer_fluid=None)
        w_old, *_ = s_old._compute_velocity(lead, lead, lead,
                                            geom, q_m3s=0.02, w_prev=w_prev,
                                            mud_fluid=mud_f, lead_fluid=None,
                                            tail_fluid=None, spacer_fluid=None)
        np.testing.assert_allclose(w_new, w_old, rtol=1e-12)
