"""BoundaryBridge  ShoeTimeline 与 provenance 重构测试

测试目标：
- 新 build_coupled_annulus_inlet_provider 直接接受 ShoeTimeline 与 WellProvenance；
- 流体角色到环空相分数的映射保持旧口径；
- build_sync_card 输出包含时间轴统计与 provenance 代理提醒。
"""

import unittest

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.provenance import (
    FluidProvenance,
    SectionProvenance,
    WellProvenance,
)
from cemdisp.models2d.boundary_bridge import (
    AnnulusInletState,
    build_coupled_annulus_inlet_provider,
    build_sync_card,
)
from cemdisp.transport1d.shoe_timeline import ShoeEvent, ShoeEventKind, ShoeTimeline


def _timeline_with_events(*events: ShoeEvent) -> ShoeTimeline:
    """辅助工厂：按时间排序构造 ShoeTimeline。"""

    return ShoeTimeline(events=tuple(sorted(events, key=lambda e: e.time_s)))


def _fluids() -> tuple[FluidSpec, ...]:
    """测试用流体，覆盖 mud/spacer/lead/intermediate/tail 角色。"""

    return (
        FluidSpec(name="泥浆", role=FluidRole.MUD, density_kg_m3=1200.0, plastic_viscosity_pa_s=0.02),
        FluidSpec(name="隔离液", role=FluidRole.SPACER, density_kg_m3=1100.0, plastic_viscosity_pa_s=0.01),
        FluidSpec(name="领浆", role=FluidRole.LEAD, density_kg_m3=1800.0, plastic_viscosity_pa_s=0.05),
        FluidSpec(name="中间浆", role=FluidRole.INTERMEDIATE, density_kg_m3=1850.0, plastic_viscosity_pa_s=0.06),
        FluidSpec(name="尾浆", role=FluidRole.TAIL, density_kg_m3=1900.0, plastic_viscosity_pa_s=0.08),
    )


def _provenance(well_name: str = "测试井") -> WellProvenance:
    """测试用 provenance，sync 状态设为非 field 以验证代理提醒。"""

    return WellProvenance(
        well_name=well_name,
        fluid={"泥浆": FluidProvenance("field", "现场实测")},
        geometry=SectionProvenance("field", "几何口径完整"),
        program=SectionProvenance("field", "施工程序完整"),
        sync=SectionProvenance("proxy", "鞋口滞后体积按双径等效计算，含代理成分。"),
    )


class TestBoundaryBridgeTimeline(unittest.TestCase):
    """测试基于 ShoeTimeline 的新桥接入口。"""

    def test_bridge_maps_tail_to_cement_phase(self) -> None:
        """尾浆角色在不拆分水泥相时映射为 cement。"""

        timeline = _timeline_with_events(
            ShoeEvent(
                time_s=0.0,
                kind=ShoeEventKind.FRONT_ARRIVAL,
                flow_rate_m3_s=0.01,
                stage_name="注尾浆",
                phase_fractions=(("尾浆", 1.0),),
            ),
        )
        provider = build_coupled_annulus_inlet_provider(
            timeline, _provenance(), _fluids(), split_cement_phases=False
        )
        state = provider(0.0)
        self.assertEqual(state.phase_fractions, (("cement", 1.0),))

    def test_bridge_maps_spacer_to_spacer_phase(self) -> None:
        """隔离液角色映射为 spacer。"""

        timeline = _timeline_with_events(
            ShoeEvent(
                time_s=10.0,
                kind=ShoeEventKind.FRONT_ARRIVAL,
                flow_rate_m3_s=0.02,
                stage_name="注隔离液",
                phase_fractions=(("隔离液", 1.0),),
            ),
        )
        provider = build_coupled_annulus_inlet_provider(
            timeline, _provenance(), _fluids(), split_cement_phases=False
        )
        state = provider(10.0)
        self.assertEqual(state.phase_fractions, (("spacer", 1.0),))

    def test_bridge_maps_mud_to_mud_phase(self) -> None:
        """泥浆角色映射为 mud。"""

        timeline = _timeline_with_events(
            ShoeEvent(
                time_s=5.0,
                kind=ShoeEventKind.FRONT_ARRIVAL,
                flow_rate_m3_s=0.015,
                stage_name="初始出流",
                phase_fractions=(("泥浆", 1.0),),
            ),
        )
        provider = build_coupled_annulus_inlet_provider(
            timeline, _provenance(), _fluids(), split_cement_phases=False
        )
        state = provider(5.0)
        self.assertEqual(state.phase_fractions, (("mud", 1.0),))

    def test_split_cement_phases_true_maps_lead_to_lead(self) -> None:
        """split_cement_phases=True 时领浆映射为 lead。"""

        timeline = _timeline_with_events(
            ShoeEvent(
                time_s=20.0,
                kind=ShoeEventKind.FRONT_ARRIVAL,
                flow_rate_m3_s=0.03,
                stage_name="注领浆",
                phase_fractions=(("领浆", 1.0),),
            ),
        )
        provider = build_coupled_annulus_inlet_provider(
            timeline, _provenance(), _fluids(), split_cement_phases=True
        )
        state = provider(20.0)
        self.assertEqual(state.phase_fractions, (("lead", 1.0),))

    def test_split_cement_phases_true_maps_tail_to_tail(self) -> None:
        """split_cement_phases=True 时尾浆映射为 tail。"""

        timeline = _timeline_with_events(
            ShoeEvent(
                time_s=30.0,
                kind=ShoeEventKind.FRONT_ARRIVAL,
                flow_rate_m3_s=0.04,
                stage_name="注尾浆",
                phase_fractions=(("尾浆", 1.0),),
            ),
        )
        provider = build_coupled_annulus_inlet_provider(
            timeline, _provenance(), _fluids(), split_cement_phases=True
        )
        state = provider(30.0)
        self.assertEqual(state.phase_fractions, (("tail", 1.0),))

    def test_provider_interpolates_between_events(self) -> None:
        """provider 在事件之间保持上一个事件的状态。"""

        timeline = _timeline_with_events(
            ShoeEvent(
                time_s=0.0,
                kind=ShoeEventKind.FRONT_ARRIVAL,
                flow_rate_m3_s=0.01,
                stage_name="初始出流",
                phase_fractions=(("泥浆", 1.0),),
            ),
            ShoeEvent(
                time_s=60.0,
                kind=ShoeEventKind.FRONT_ARRIVAL,
                flow_rate_m3_s=0.02,
                stage_name="注尾浆",
                phase_fractions=(("尾浆", 1.0),),
            ),
        )
        provider = build_coupled_annulus_inlet_provider(
            timeline, _provenance(), _fluids(), split_cement_phases=False
        )
        before = provider(30.0)
        after = provider(60.0)
        self.assertEqual(before.phase_fractions, (("mud", 1.0),))
        self.assertEqual(after.phase_fractions, (("cement", 1.0),))

    def test_unknown_fluid_defaults_to_mud(self) -> None:
        """时间轴中出现未知流体名称时，默认映射为 mud。"""

        timeline = _timeline_with_events(
            ShoeEvent(
                time_s=0.0,
                kind=ShoeEventKind.FRONT_ARRIVAL,
                flow_rate_m3_s=0.01,
                stage_name="未知阶段",
                phase_fractions=(("未知流体", 1.0),),
            ),
        )
        provider = build_coupled_annulus_inlet_provider(
            timeline, _provenance(), _fluids(), split_cement_phases=False
        )
        state = provider(0.0)
        self.assertEqual(state.phase_fractions, (("mud", 1.0),))


class TestSyncCard(unittest.TestCase):
    """测试 build_sync_card 输出。"""

    def test_sync_card_contains_event_count(self) -> None:
        """画像卡包含时间轴事件数与首尾时间。"""

        timeline = _timeline_with_events(
            ShoeEvent(
                time_s=10.0,
                kind=ShoeEventKind.FRONT_ARRIVAL,
                flow_rate_m3_s=0.01,
                stage_name="阶段1",
                phase_fractions=(("泥浆", 1.0),),
            ),
            ShoeEvent(
                time_s=120.0,
                kind=ShoeEventKind.END,
                flow_rate_m3_s=0.0,
                stage_name="结束",
                phase_fractions=(),
            ),
        )
        provenance = _provenance("呼102")
        card = build_sync_card("呼102", timeline, provenance)

        self.assertEqual(card["井名"], "呼102")
        sync_info = card["鞋口同步口径"]
        self.assertEqual(sync_info["事件数"], 2)
        self.assertEqual(sync_info["首事件时间_s"], 10.0)
        self.assertEqual(sync_info["末事件时间_s"], 120.0)

    def test_sync_card_shows_proxy_note_when_sync_not_field(self) -> None:
        """sync 状态非 field 时，代理提醒应包含 provenance 说明。"""

        timeline = ShoeTimeline(events=())
        provenance = _provenance()
        card = build_sync_card("测试井", timeline, provenance)

        self.assertIn("代理", card["代理提醒"])

    def test_sync_card_empty_proxy_note_when_sync_is_field(self) -> None:
        """sync 状态为 field 时，代理提醒应为空字符串。"""

        timeline = ShoeTimeline(events=())
        provenance = WellProvenance(
            well_name="测试井",
            fluid={},
            geometry=SectionProvenance("field", ""),
            program=SectionProvenance("field", ""),
            sync=SectionProvenance("field", "全部现场实测"),
        )
        card = build_sync_card("测试井", timeline, provenance)

        self.assertEqual(card["代理提醒"], "")


class TestFlusherPhaseMapping(unittest.TestCase):
    """测试 FLUSHER 流体→"flusher" 相名映射（T1-6 Task 2）。"""

    def test_flusher_maps_to_flusher_phase(self) -> None:
        """FLUSHER 角色映射为独立 flusher 相，不并入 spacer。"""

        timeline = _timeline_with_events(
            ShoeEvent(
                time_s=0.0,
                kind=ShoeEventKind.FRONT_ARRIVAL,
                flow_rate_m3_s=0.01,
                stage_name="注冲洗液",
                phase_fractions=(("冲洗液", 1.0),),
            ),
        )
        fluids = _fluids() + (
            FluidSpec(name="冲洗液", role=FluidRole.FLUSHER, density_kg_m3=1500.0, plastic_viscosity_pa_s=0.03),
        )
        provider = build_coupled_annulus_inlet_provider(
            timeline, _provenance(), fluids, split_cement_phases=False
        )
        state = provider(0.0)
        self.assertEqual(state.phase_fractions, (("flusher", 1.0),))

    def test_flusher_does_not_merge_into_spacer(self) -> None:
        """FLUSHER 与 SPACER 同时存在时，flusher 保持独立不合并到 spacer。"""

        timeline = _timeline_with_events(
            ShoeEvent(
                time_s=0.0,
                kind=ShoeEventKind.FRONT_ARRIVAL,
                flow_rate_m3_s=0.01,
                stage_name="混合出流",
                phase_fractions=(("隔离液", 0.5), ("冲洗液", 0.5)),
            ),
        )
        fluids = _fluids() + (
            FluidSpec(name="冲洗液", role=FluidRole.FLUSHER, density_kg_m3=1500.0, plastic_viscosity_pa_s=0.03),
        )
        provider = build_coupled_annulus_inlet_provider(
            timeline, _provenance(), fluids, split_cement_phases=False
        )
        state = provider(0.0)
        mapped = dict(state.phase_fractions)
        self.assertIn("flusher", mapped, "FLUSHER 应独立存在于 phase_fractions")
        self.assertIn("spacer", mapped, "SPACER 应独立存在于 phase_fractions")
        self.assertAlmostEqual(mapped["spacer"] + mapped["flusher"], 1.0,
                               msg="spacer 与 flusher 分数之和应为 1.0")

    def test_annulus_inlet_state_supports_flusher(self) -> None:
        """AnnulusInletState 直接构造含 flusher 分数不报错。"""

        state = AnnulusInletState(
            time_s=0.0,
            flow_rate_m3_s=0.01,
            stage_name="测试",
            phase_fractions=(("flusher", 1.0),),
        )
        mapped = dict(state.phase_fractions)
        self.assertIn("flusher", mapped)
        self.assertEqual(mapped["flusher"], 1.0)

    def test_flusher_with_split_cement_phases(self) -> None:
        """split_cement_phases=True 时 FLUSHER 仍映射为独立 flusher。"""

        timeline = _timeline_with_events(
            ShoeEvent(
                time_s=0.0,
                kind=ShoeEventKind.FRONT_ARRIVAL,
                flow_rate_m3_s=0.01,
                stage_name="注冲洗液",
                phase_fractions=(("冲洗液", 1.0),),
            ),
        )
        fluids = _fluids() + (
            FluidSpec(name="冲洗液", role=FluidRole.FLUSHER, density_kg_m3=1500.0, plastic_viscosity_pa_s=0.03),
        )
        provider = build_coupled_annulus_inlet_provider(
            timeline, _provenance(), fluids, split_cement_phases=True
        )
        state = provider(0.0)
        self.assertEqual(state.phase_fractions, (("flusher", 1.0),))


if __name__ == "__main__":
    _ = unittest.main()
