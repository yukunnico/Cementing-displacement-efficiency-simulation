"""
T0-2 温变流变公式全集契约测试（cemdisp.data.rheology_vs_temperature）

覆盖 brief 需求（seam=fluid_at 公共接口；期望值一律取自 §3.7 回归锚点，
禁止用实现公式自算）：

- 水泥组 A / 组 B 锚点：7 行 × 2（含 100−=T=100、100+=T=100.1 两个独立断点用例）
- 钻井液锚点（T=60）+ clamp 关系断言（T=20→40、T=100→80）+ 审计计数 +1
- 隔离液锚点（T=20, P=0.1）：1.95 → τy≈7.28/μp≈0.0779；
  2.05 → τy≈10.28/μp≈0.1054（2026-09-30 Q1① 裁定修正锚，原 10.75/0.1233 废止）
- 插值单调性：ρ=1.90/1.93/1.95/2.05/2.10 序列 τ₀(60°C) 严格递增
- 分派表逐行一例：钻井液/先导浆/平衡液/替浆链不替换/隔离液插值/隔离液域外借用/
  领尾中间浆三密度档/冲洗液不替换/不匹配不替换
- Bingham 绝对替换切换、smooth_break 占位、clamp/borrow/审计 API
"""

from __future__ import annotations

import pytest

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.rheology_vs_temperature import _route, fluid_at, get_audit, reset_audit

# 锚点断言容差（brief：τ ±0.01 Pa、μp ±0.5%）
TAU_TOL = 0.01
MUP_REL = 0.005


@pytest.fixture(autouse=True)
def _clean_audit():
    """每个用例前后清空模块级审计，保证计数断言独立。"""
    reset_audit()
    yield
    reset_audit()


# ---------------------------------------------------------------------------
# 流体工厂（输入侧任意合法模型；绝对替换后应全部切 Bingham）
# ---------------------------------------------------------------------------
def _cement(name: str, rho_g_cm3: float, role: FluidRole = FluidRole.LEAD) -> FluidSpec:
    return FluidSpec(
        name=name,
        role=role,
        density_kg_m3=rho_g_cm3 * 1000.0,
        rheology_model=RheologyModel.POWER_LAW,
        power_law_n=0.8,
        consistency_k=0.5,
    )


def _spacer(name: str = "隔离液1", rho_g_cm3: float = 1.95) -> FluidSpec:
    return FluidSpec(
        name=name,
        role=FluidRole.SPACER,
        density_kg_m3=rho_g_cm3 * 1000.0,
        rheology_model=RheologyModel.POWER_LAW,
        power_law_n=0.7,
        consistency_k=0.3,
    )


def _mud(name: str = "钻井液", role: FluidRole = FluidRole.MUD) -> FluidSpec:
    return FluidSpec(
        name=name,
        role=role,
        density_kg_m3=1200.0,
        rheology_model=RheologyModel.BINGHAM,
        plastic_viscosity_pa_s=0.05,
        yield_stress_pa=5.0,
    )


def _assert_tau(actual_pa: float, expected_pa: float) -> None:
    assert actual_pa == pytest.approx(expected_pa, abs=TAU_TOL)


def _assert_mup(actual_pa_s: float, expected_mpas: float) -> None:
    """μp：水泥/mud 锚点以 mPa·s 给出，按 ÷1000 换算为 Pa·s 后比相对误差。"""
    assert actual_pa_s == pytest.approx(expected_mpas / 1000.0, rel=MUP_REL)


def _assert_mup_pas(actual_pa_s: float, expected_pa_s: float) -> None:
    """μp：隔离液锚点以 Pa·s 给出（brief 注明单位 Pa、Pa·s），直接比相对误差。"""
    assert actual_pa_s == pytest.approx(expected_pa_s, rel=MUP_REL)


def _clamp_count() -> int:
    return sum(1 for e in get_audit() if e["kind"] == "clamp")


# §3.7 回归锚点（T, τ₀_A, μp_A(mPa·s), τ₀_B, μp_B(mPa·s)）；100+=T=100.1（§四 裁定）
CEMENT_ANCHORS = [
    pytest.param(20.0, 14.921, 411.8, 4.571, 74.4, id="20"),
    pytest.param(60.0, 8.190, 249.0, 2.852, 91.5, id="60"),
    pytest.param(100.0, 10.398, 150.3, 2.552, 151.9, id="100-"),
    pytest.param(100.1, 10.756, 150.1, 5.188, 152.2, id="100+"),
    pytest.param(120.0, 23.191, 116.6, 5.826, 198.4, id="120"),
    pytest.param(155.0, 45.062, 74.6, 6.948, 305.9, id="155"),
    pytest.param(170.0, 54.435, 61.5, 7.429, 362.1, id="170"),
]


# ---------------------------------------------------------------------------
# 水泥组 A（ρ=2.10 → cement_A_2p1 档）
# ---------------------------------------------------------------------------
class TestCementAnchorA:
    @pytest.mark.parametrize("T,tau_a,mpa_a,tau_b,mpb_b", CEMENT_ANCHORS)
    def test_anchor_A(self, T, tau_a, mpa_a, tau_b, mpb_b):
        out = fluid_at(_cement("领浆", 2.10), T)
        _assert_tau(out.yield_stress_pa, tau_a)
        _assert_mup(out.plastic_viscosity_pa_s, mpa_a)


# ---------------------------------------------------------------------------
# 水泥组 B（ρ=1.90 → cement_B_1p9 档）
# ---------------------------------------------------------------------------
class TestCementAnchorB:
    @pytest.mark.parametrize("T,tau_a,mpa_a,tau_b,mpb_b", CEMENT_ANCHORS)
    def test_anchor_B(self, T, tau_a, mpa_a, tau_b, mpb_b):
        out = fluid_at(_cement("尾浆", 1.90), T)
        _assert_tau(out.yield_stress_pa, tau_b)
        _assert_mup(out.plastic_viscosity_pa_s, mpb_b)


# ---------------------------------------------------------------------------
# 钻井液锚点 + clamp 关系断言
# ---------------------------------------------------------------------------
class TestMudAnchorAndClamp:
    def test_anchor_60(self):
        out = fluid_at(_mud(), 60.0)
        _assert_tau(out.yield_stress_pa, 12.97)
        _assert_mup(out.plastic_viscosity_pa_s, 74.3)

    def test_clamp_T20_equals_T40_and_audit(self):
        out20 = fluid_at(_mud(), 20.0)
        out40 = fluid_at(_mud(), 40.0)
        assert out20.yield_stress_pa == out40.yield_stress_pa
        assert out20.plastic_viscosity_pa_s == out40.plastic_viscosity_pa_s
        assert _clamp_count() == 1
        ev = [e for e in get_audit() if e["kind"] == "clamp"]
        assert ev[0]["fluid"] == "钻井液"
        assert ev[0]["param"] == "T"
        assert ev[0]["requested"] == pytest.approx(20.0)
        assert ev[0]["clamped"] == pytest.approx(40.0)

    def test_clamp_T100_equals_T80_and_audit(self):
        out100 = fluid_at(_mud(), 100.0)
        out80 = fluid_at(_mud(), 80.0)
        assert out100.yield_stress_pa == out80.yield_stress_pa
        assert out100.plastic_viscosity_pa_s == out80.plastic_viscosity_pa_s
        assert _clamp_count() == 1
        ev = [e for e in get_audit() if e["kind"] == "clamp"]
        assert ev[0]["clamped"] == pytest.approx(80.0)

    def test_in_domain_no_clamp(self):
        fluid_at(_mud(), 60.0)
        assert _clamp_count() == 0

    def test_mud_extrapolate_switch(self):
        clamped = fluid_at(_mud(), 20.0)
        reset_audit()
        ext = fluid_at(_mud(), 20.0, mud_extrapolate=True)
        # 外推值 ≠ clamp 到 40°C 的端值，且不记 clamp、记 extrapolate
        assert ext.yield_stress_pa != clamped.yield_stress_pa
        assert ext.plastic_viscosity_pa_s != clamped.plastic_viscosity_pa_s
        assert _clamp_count() == 0
        assert any(e["kind"] == "extrapolate" for e in get_audit())


# ---------------------------------------------------------------------------
# 隔离液锚点（T=20, P=0.1；2.05 档为 Q1① 修正锚）
# ---------------------------------------------------------------------------
class TestSpacerAnchor:
    def test_1p95(self):
        out = fluid_at(_spacer("隔离液1", 1.95), 20.0, 0.1)
        _assert_tau(out.yield_stress_pa, 7.28)
        _assert_mup_pas(out.plastic_viscosity_pa_s, 0.0779)

    def test_2p05(self):
        out = fluid_at(_spacer("隔离液1", 2.05), 20.0, 0.1)
        _assert_tau(out.yield_stress_pa, 10.28)
        _assert_mup_pas(out.plastic_viscosity_pa_s, 0.1054)


# ---------------------------------------------------------------------------
# 中间密度插值单调性
# ---------------------------------------------------------------------------
class TestCementInterpMonotonic:
    def test_tau0_increases_with_density(self):
        rhos = [1.90, 1.93, 1.95, 2.05, 2.10]
        taus = [fluid_at(_cement("领浆", r), 60.0).yield_stress_pa for r in rhos]
        assert all(b > a for a, b in zip(taus, taus[1:])), taus


# ---------------------------------------------------------------------------
# 分派表逐行一例
# ---------------------------------------------------------------------------
class TestDispatchTable:
    def test_row_mud_drilling_fluid(self):
        """钻井液 → 钻井液式（以 mud 锚点证明路由）。"""
        out = fluid_at(_mud("钻井液", FluidRole.MUD), 60.0)
        _assert_tau(out.yield_stress_pa, 12.97)
        _assert_mup(out.plastic_viscosity_pa_s, 74.3)
        assert not any(e["kind"] == "model_assumption" for e in get_audit())

    @pytest.mark.parametrize(
        "name,role",
        [
            pytest.param("先导浆", FluidRole.WASH, id="先导浆"),
            pytest.param("平衡液", FluidRole.WASH, id="平衡液"),
        ],
    )
    def test_row_assumption_mud(self, name, role):
        """先导浆、平衡液 → 钻井液式 + model_assumption 标注。"""
        f = FluidSpec(
            name=name, role=role, density_kg_m3=1100.0,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=0.02, yield_stress_pa=3.0,
        )
        out = fluid_at(f, 60.0)
        mud = fluid_at(_mud(), 60.0)
        assert out.yield_stress_pa == mud.yield_stress_pa
        assert out.plastic_viscosity_pa_s == mud.plastic_viscosity_pa_s
        assert any(
            e["kind"] == "model_assumption" and e["fluid"] == name
            for e in get_audit()
        )

    def test_row_hu1_balance_alias(self):
        """hu1「平衡液(先导泥浆)」(role=WASH) → 钻井液式 + model_assumption。

        review Important：全等匹配漏掉该别名会落 unmatched 不替换；
        改子串匹配后必须命中行2（回归用例）。
        """
        f = FluidSpec(
            name="平衡液(先导泥浆)", role=FluidRole.WASH, density_kg_m3=1100.0,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=0.02, yield_stress_pa=3.0,
        )
        out = fluid_at(f, 60.0)
        mud = fluid_at(_mud(), 60.0)
        assert out.yield_stress_pa == mud.yield_stress_pa
        assert out.plastic_viscosity_pa_s == mud.plastic_viscosity_pa_s
        assert out.rheology_model is RheologyModel.BINGHAM
        assert any(
            e["kind"] == "model_assumption" and e["fluid"] == "平衡液(先导泥浆)"
            for e in get_audit()
        )
        assert not any(e["kind"] == "no_replace" for e in get_audit())

    def test_row_hu103_balance_mud_role_still_assumption(self):
        """hu103「平衡液」role=MUD → 名字档优先于 role（review minor② 裁定）。

        公式同为钻井液式，差异仅在审计标注：必须带 model_assumption。
        """
        f = FluidSpec(
            name="平衡液", role=FluidRole.MUD, density_kg_m3=1050.0,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=0.02, yield_stress_pa=2.0,
        )
        out = fluid_at(f, 60.0)
        mud = fluid_at(_mud(), 60.0)
        assert out.yield_stress_pa == mud.yield_stress_pa
        assert out.plastic_viscosity_pa_s == mud.plastic_viscosity_pa_s
        assert any(
            e["kind"] == "model_assumption" and e["fluid"] == "平衡液"
            for e in get_audit()
        )

    @pytest.mark.parametrize(
        "name,role",
        [
            pytest.param("压塞液", FluidRole.OTHER, id="压塞液"),
            pytest.param("替钻井液", FluidRole.DISPLACEMENT, id="替钻井液"),
            pytest.param("井浆", FluidRole.DISPLACEMENT, id="井浆"),
            pytest.param("基液", FluidRole.DISPLACEMENT, id="基液"),
            pytest.param("保护液", FluidRole.DISPLACEMENT, id="保护液"),
        ],
    )
    def test_row_replacement_chain_no_replace(self, name, role):
        """替浆链（压塞液/替钻井液/井浆/基液/保护液）→ 不替换，原对象返回。"""
        f = FluidSpec(
            name=name, role=role, density_kg_m3=1200.0,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=0.05, yield_stress_pa=5.0,
        )
        out = fluid_at(f, 60.0)
        assert out is f
        assert any(
            e["kind"] == "no_replace" and e["fluid"] == name
            for e in get_audit()
        )

    def test_row_spacer_density_interp(self):
        """隔离液 ρ∈[1.95,2.05] → 六系数密度插值（ρ=2.00 = 两端锚点均值，值线性）。"""
        out = fluid_at(_spacer("隔离液1", 2.00), 20.0, 0.1)
        _assert_tau(out.yield_stress_pa, (7.28 + 10.28) / 2.0)
        _assert_mup_pas(out.plastic_viscosity_pa_s, (0.0779 + 0.1054) / 2.0)
        assert not any(e["kind"] == "borrow" for e in get_audit())

    @pytest.mark.parametrize(
        "rho,ref_rho",
        [
            pytest.param(1.75, 1.95, id="1.75就近借1.95"),
            pytest.param(1.92, 1.95, id="1.92就近借1.95"),
            pytest.param(2.10, 2.05, id="2.10就近借2.05"),
        ],
    )
    def test_row_spacer_out_of_range_borrow(self, rho, ref_rho):
        """隔离液 ρ 域外 → clamp 到最近端公式 + model_assumption 借用审计。"""
        name = "隔离液2" if rho == 1.75 else "驱油隔离液"
        f = _spacer(name, rho)
        out = fluid_at(f, 20.0, 0.1)
        ref = fluid_at(_spacer("参照", ref_rho), 20.0, 0.1)
        # 系数不外推：借用=与端点完全同一式
        assert out.yield_stress_pa == ref.yield_stress_pa
        assert out.plastic_viscosity_pa_s == ref.plastic_viscosity_pa_s
        borrow_events = [
            e for e in get_audit()
            if e["kind"] == "borrow" and e["fluid"] == name
        ]
        assert borrow_events, "缺 borrow 审计事件"
        assert borrow_events[0].get("note") == "model_assumption"

    # ------------------------------------------------------------------ #
    # C1 回归（2026-10-01 终审修复）：role==SPACER 判定优先于名字子串
    # ------------------------------------------------------------------ #
    @staticmethod
    def _synthetic_spacer() -> FluidSpec:
        """2D 合成等效隔离液（annulus_d2dga._composite_spacer_fluid 同型）。

        名字含子串「先导浆」+ role=SPACER + ρ≈1.82（域外）——修复前名字子串
        先于 role 判定 ⇒ 误落泥浆式（τy≈10.98、T 被误 clamp [40,80]）。
        """
        return FluidSpec(
            name="先导浆+隔离液1+隔离液2", role=FluidRole.SPACER,
            density_kg_m3=1820.0, rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=0.05, yield_stress_pa=8.0,
        )

    def test_c1_synthetic_spacer_routes_to_spacer_family(self):
        """合成相（名含「先导浆」、role=SPACER）→ 隔离液族 + borrow 审计。

        路由必须 role==SPACER 优先于名字子串：ρ=1.82 域外就近借 1.95 式，
        与隔离液族参照逐位同值；不得出现 model_assumption 泥浆式标注、
        也不得有泥浆域 clamp（[40,80]）事件。
        """
        f = self._synthetic_spacer()
        assert _route(f) == ("spacer", None)

        out = fluid_at(f, 60.0)                       # P 缺省=0.1（p_default 审计）
        ref = fluid_at(_spacer("参照", 1.95), 60.0, 0.1)
        assert out.yield_stress_pa == ref.yield_stress_pa
        assert out.plastic_viscosity_pa_s == ref.plastic_viscosity_pa_s
        # 与泥浆式（同 T）必须不同值——修复前两者相等即为缺陷
        assert out.yield_stress_pa != fluid_at(_mud(), 60.0).yield_stress_pa

        events = get_audit()
        assert any(e["kind"] == "borrow" and e["fluid"] == f.name for e in events)
        assert not any(
            e["kind"] == "model_assumption" and e["fluid"] == f.name for e in events
        )
        assert not any(
            e["kind"] == "clamp" and e["fluid"] == f.name for e in events
        )

    def test_c1_standalone_lead_mud_still_mud_family(self):
        """单独「先导浆」(role=WASH) 仍走泥浆式——子串档只对非 SPACER 生效。"""
        f = FluidSpec(
            name="先导浆", role=FluidRole.WASH, density_kg_m3=1750.0,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=0.058, yield_stress_pa=9.8,
        )
        assert _route(f) == ("mud", "model_assumption")
        out = fluid_at(f, 60.0)
        mud = fluid_at(_mud(), 60.0)
        assert out.yield_stress_pa == mud.yield_stress_pa
        assert out.plastic_viscosity_pa_s == mud.plastic_viscosity_pa_s
        assert any(
            e["kind"] == "model_assumption" and e["fluid"] == "先导浆"
            for e in get_audit()
        )

    def test_row_cement_B_band(self):
        """领/尾/中间浆 ρ∈[1.88,1.92] → cement_B_1p9（以 B 锚点证明路由）。"""
        out = fluid_at(_cement("尾浆", 1.90), 60.0)
        _assert_tau(out.yield_stress_pa, 2.852)
        _assert_mup(out.plastic_viscosity_pa_s, 91.5)

    def test_row_cement_A_band(self):
        """领/尾/中间浆 ρ∈[2.08,2.12] → cement_A_2p1（以 A 锚点证明路由）。"""
        out = fluid_at(_cement("领浆", 2.10), 60.0)
        _assert_tau(out.yield_stress_pa, 8.190)
        _assert_mup(out.plastic_viscosity_pa_s, 249.0)

    def test_row_cement_interp_band(self):
        """领/尾/中间浆 ρ∈(1.92,2.08) → cement_interp(ρ)：期望由 A/B 锚点线性组合。"""
        out = fluid_at(_cement("领浆", 1.93), 60.0)
        _assert_tau(out.yield_stress_pa, 2.852 + 0.15 * (8.190 - 2.852))
        _assert_mup(out.plastic_viscosity_pa_s, 91.5 + 0.15 * (249.0 - 91.5))

    @pytest.mark.parametrize(
        "rho,expected_tau",
        [
            pytest.param(1.90, 2.852, id="中间浆1.90→B"),
            pytest.param(2.05, 2.852 + 0.75 * (8.190 - 2.852), id="中间浆2.05→interp"),
        ],
    )
    def test_row_intermediate_merged_into_cement(self, rho, expected_tau):
        """中间浆 INTERMEDIATE 并入水泥档，按 ρ 分派（补裁③）。"""
        out = fluid_at(
            _cement("中间浆", rho, role=FluidRole.INTERMEDIATE), 60.0
        )
        _assert_tau(out.yield_stress_pa, expected_tau)

    def test_row_flusher_no_replace(self):
        """冲洗液（FLUSHER）→ 不替换 + 审计标记，原对象返回。"""
        f = FluidSpec(
            name="冲洗液（FLUSHER）", role=FluidRole.FLUSHER,
            density_kg_m3=1050.0, rheology_model=RheologyModel.NEWTONIAN,
            plastic_viscosity_pa_s=0.001,
        )
        out = fluid_at(f, 60.0)
        assert out is f
        assert any(
            e["kind"] == "no_replace" and e["fluid"] == f.name
            for e in get_audit()
        )

    def test_row_unmatched_no_replace(self):
        """LEAD/TAIL/INTERMEDIATE 之外不匹配任何档 → 不替换 + 审计标记。"""
        f = FluidSpec(
            name="特殊液", role=FluidRole.OTHER, density_kg_m3=1200.0,
            rheology_model=RheologyModel.NEWTONIAN,
            plastic_viscosity_pa_s=0.001,
        )
        out = fluid_at(f, 60.0)
        assert out is f
        assert any(
            e["kind"] == "no_replace" and "不匹配" in e["detail"]
            for e in get_audit()
        )


# ---------------------------------------------------------------------------
# Bingham 绝对替换切换
# ---------------------------------------------------------------------------
class TestBinghamAbsoluteReplacement:
    def test_switch_and_clear_power_law(self):
        src = _cement("领浆", 1.93)  # 输入为幂律
        out = fluid_at(src, 60.0)
        assert out.rheology_model is RheologyModel.BINGHAM
        assert out.yield_stress_pa is not None and out.yield_stress_pa > 0
        assert out.plastic_viscosity_pa_s is not None and out.plastic_viscosity_pa_s > 0
        assert out.power_law_n is None
        assert out.consistency_k is None
        # 密度/名字/角色不变；原对象（frozen）未被改动
        assert out.density_kg_m3 == src.density_kg_m3
        assert out.name == src.name
        assert out.role == src.role
        assert src.rheology_model is RheologyModel.POWER_LAW
        assert src.power_law_n == 0.8
        assert src.consistency_k == 0.5


# ---------------------------------------------------------------------------
# smooth_break 占位（默认关=忠实两段式）
# ---------------------------------------------------------------------------
class TestSmoothBreakPlaceholder:
    def test_default_faithful_two_segment(self):
        out = fluid_at(_cement("尾浆", 1.90), 100.0, smooth_break=False)
        _assert_tau(out.yield_stress_pa, 2.552)  # T=100 走低温式

    def test_true_raises_notimplemented(self):
        with pytest.raises(NotImplementedError):
            fluid_at(_cement("尾浆", 1.90), 100.0, smooth_break=True)


# ---------------------------------------------------------------------------
# 审计 API + 压力语义
# ---------------------------------------------------------------------------
class TestAuditAndPressure:
    def test_get_audit_and_reset(self):
        f = FluidSpec(
            name="冲洗液（FLUSHER）", role=FluidRole.FLUSHER,
            density_kg_m3=1050.0, rheology_model=RheologyModel.NEWTONIAN,
            plastic_viscosity_pa_s=0.001,
        )
        fluid_at(f, 60.0)
        events = get_audit()
        assert isinstance(events, list) and len(events) == 1
        assert isinstance(events[0], dict)
        assert events[0]["kind"] == "no_replace"
        reset_audit()
        assert get_audit() == []

    def test_cement_T_clamp(self):
        """水泥域越界 clamp：A 域 [20,170]（180→170 与 10→20 各 +1 审计）。"""
        out180 = fluid_at(_cement("领浆", 2.10), 180.0)
        out170 = fluid_at(_cement("领浆", 2.10), 170.0)
        assert out180.yield_stress_pa == out170.yield_stress_pa
        assert _clamp_count() == 1
        reset_audit()
        out10 = fluid_at(_cement("领浆", 2.10), 10.0)
        out20 = fluid_at(_cement("领浆", 2.10), 20.0)
        assert out10.yield_stress_pa == out20.yield_stress_pa
        assert _clamp_count() == 1

    def test_spacer_P_clamp(self):
        out = fluid_at(_spacer(), 20.0, 250.0)
        ref = fluid_at(_spacer(), 20.0, 200.0)
        assert out.yield_stress_pa == ref.yield_stress_pa
        assert _clamp_count() == 1
        ev = [e for e in get_audit() if e["kind"] == "clamp"]
        assert ev[0]["param"] == "P"
        assert ev[0]["clamped"] == pytest.approx(200.0)

    def test_spacer_P_default_constant_0p1(self):
        """P_mpa=None → 常压 0.1 MPa 求值 + p_default 审计。"""
        out = fluid_at(_spacer(), 20.0)
        ref = fluid_at(_spacer(), 20.0, 0.1)
        assert out.yield_stress_pa == ref.yield_stress_pa
        assert out.plastic_viscosity_pa_s == ref.plastic_viscosity_pa_s
        assert any(e["kind"] == "p_default" for e in get_audit())

    def test_non_spacer_ignores_P(self):
        """非隔离液相忽略 P（传入也不 clamp/不审计）。"""
        fluid_at(_cement("尾浆", 1.90), 60.0, 999.0)
        assert get_audit() == []
