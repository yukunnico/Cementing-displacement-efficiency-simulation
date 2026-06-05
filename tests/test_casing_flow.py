"""
CasingFlowSolver 鞋口时间轴输出测试

测试目标：
- run() 保持原有 CasingFlowResult 字段行为，同时附带 ShoeTimeline；
- ShoeTimeline 按现有体积追踪结果给出鞋口出流状态；
- 双径向井上段内径参与新时间轴的鞋口迟到体积计算，不回改旧字段口径。
"""

import math
import unittest

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.pumping_schedule import (
    PumpingSchedule,
    PumpingScheduleStep,
    PumpingStageEvent,
)
from cemdisp.data.well_spec import WellSpec
from cemdisp.transport1d.casing_flow import CasingFlowSolver
from cemdisp.transport1d.shoe_timeline import ShoeEventKind, ShoeTimeline


def _liner_id_for_area(area_m2: float) -> float:
    """按目标截面积反算内径，便于测试中直接控制管内容积。"""

    return math.sqrt(4.0 * area_m2 / math.pi) * 1000.0


def _fluids() -> tuple[FluidSpec, ...]:
    """构建测试用流体，名称与施工步骤保持一致。"""

    return (
        FluidSpec(name="泥浆", role=FluidRole.MUD, density_kg_m3=1200.0, plastic_viscosity_pa_s=0.02),
        FluidSpec(name="隔离液", role=FluidRole.SPACER, density_kg_m3=1100.0, plastic_viscosity_pa_s=0.01),
        FluidSpec(name="尾浆", role=FluidRole.TAIL, density_kg_m3=1900.0, plastic_viscosity_pa_s=0.08),
        FluidSpec(name="顶替液", role=FluidRole.DISPLACEMENT, density_kg_m3=1050.0, plastic_viscosity_pa_s=0.01),
    )


def _well(*, upper_area_m2: float | None = None) -> WellSpec:
    """构建鞋深 100 m 的测试井；默认单一内径对应 1 m³ 管内容积。"""

    if upper_area_m2 is None:
        return WellSpec(
            well_name="测试井",
            top_md_m=1.0,
            bottom_md_m=120.0,
            shoe_md_m=100.0,
            liner_id_mm=_liner_id_for_area(0.01),
        )
    return WellSpec(
        well_name="测试井",
        top_md_m=1.0,
        bottom_md_m=120.0,
        shoe_md_m=100.0,
        liner_id_mm=_liner_id_for_area(0.01),
        upper_section_bottom_md_m=50.0,
        upper_liner_od_mm=200.0,
        upper_liner_id_mm=_liner_id_for_area(upper_area_m2),
    )


def _schedule() -> PumpingSchedule:
    """构建含隔离液、水泥浆、停泵、重启顶替的显式施工程序。"""

    return PumpingSchedule(
        steps=(
            PumpingScheduleStep(
                step_name="注入隔离液",
                fluid_name="隔离液",
                volume_m3=1.0,
                rate_m3_min=1.0,
                start_time_s=0.0,
                end_time_s=60.0,
                event_tag=PumpingStageEvent.INJECT_SPACER,
            ),
            PumpingScheduleStep(
                step_name="注入尾浆",
                fluid_name="尾浆",
                volume_m3=1.0,
                rate_m3_min=1.0,
                start_time_s=60.0,
                end_time_s=120.0,
                event_tag=PumpingStageEvent.INJECT_CEMENT,
            ),
            PumpingScheduleStep(
                step_name="停泵候凝",
                fluid_name="顶替液",
                volume_m3=0.0,
                rate_m3_min=0.0,
                start_time_s=120.0,
                end_time_s=180.0,
                event_tag=PumpingStageEvent.SHUTDOWN,
            ),
            PumpingScheduleStep(
                step_name="重启顶替",
                fluid_name="顶替液",
                volume_m3=2.0,
                rate_m3_min=2.0,
                start_time_s=180.0,
                end_time_s=240.0,
                event_tag=PumpingStageEvent.RESTART,
            ),
        )
    )


class TestCasingFlowShoeTimeline(unittest.TestCase):
    """测试 CasingFlowSolver 输出 ShoeTimeline 的新增路径。"""

    def test_run_emits_shoe_timeline_without_changing_legacy_result_fields(self) -> None:
        """run() 返回 ShoeTimeline，且单径井查询结果与旧 pipe_exit_state_at 口径一致。"""

        solver = CasingFlowSolver(enable_gravity=False, enable_axial_dispersion=False)
        result = solver.run(_well(), _fluids(), _schedule())

        self.assertIsInstance(result.shoe_timeline, ShoeTimeline)
        self.assertAlmostEqual(result.pipe_cross_section_m2, 0.01)
        self.assertEqual(result.shoe_md_m, 100.0)
        self.assertEqual(result.schedule_steps, _schedule().steps)

        kinds = [event.kind for event in result.shoe_timeline.events]
        self.assertIn(ShoeEventKind.FRONT_ARRIVAL, kinds)
        self.assertIn(ShoeEventKind.REAR_EXIT, kinds)
        self.assertIn(ShoeEventKind.SHUTDOWN, kinds)
        self.assertIn(ShoeEventKind.RESTART, kinds)
        self.assertEqual(result.shoe_timeline.events[-1].kind, ShoeEventKind.END)

        for time_s in (0.0, 60.0, 150.0, 210.0, 240.0):
            legacy_state = solver.pipe_exit_state_at(result, time_s)
            timeline_state = result.shoe_timeline.at(time_s)
            self.assertEqual(timeline_state.flow_rate_m3_s, legacy_state.flow_rate_m3_s)
            self.assertEqual(timeline_state.stage_name, legacy_state.stage_name)
            self.assertEqual(timeline_state.phase_fractions, legacy_state.phase_fractions)

    def test_dual_diameter_metadata_delays_timeline_without_rewriting_legacy_fronts(self) -> None:
        """双径向井时间轴使用上段内径累计体积；旧 fronts 仍保持单径算法结果。"""

        solver = CasingFlowSolver(enable_gravity=False, enable_axial_dispersion=False)
        result = solver.run(_well(upper_area_m2=0.02), _fluids(), _schedule())

        # 旧结果仍按 result.pipe_cross_section_m2 * shoe_md_m = 1 m³ 计算首个前缘。
        self.assertAlmostEqual(result.fronts[0].time_s, 60.0)

        before_arrival = result.shoe_timeline.at(89.0)
        at_arrival = result.shoe_timeline.at(90.0)

        self.assertEqual(before_arrival.phase_fractions, (("泥浆", 1.0),))
        self.assertEqual(at_arrival.phase_fractions, (("隔离液", 1.0),))

    def test_gravity_corrected_front_time_is_used_by_timeline(self) -> None:
        """启用重力修正时，时间轴前缘事件必须与旧 fronts 的校正时间一致。"""

        solver = CasingFlowSolver(enable_gravity=True, settling_velocity_factor=1.0, enable_axial_dispersion=False)
        result = solver.run(_well(), _fluids(), _schedule())

        corrected_front_time_s = result.fronts[0].time_s
        self.assertAlmostEqual(corrected_front_time_s, 54.0)

        front_events = [
            event for event in result.shoe_timeline.events
            if event.kind == ShoeEventKind.FRONT_ARRIVAL and event.phase_fractions == (("隔离液", 1.0),)
        ]
        self.assertTrue(front_events)
        self.assertAlmostEqual(front_events[0].time_s, corrected_front_time_s)

        timeline_state = result.shoe_timeline.at(corrected_front_time_s)
        self.assertEqual(timeline_state.phase_fractions, (("隔离液", 1.0),))


if __name__ == "__main__":
    _ = unittest.main()
