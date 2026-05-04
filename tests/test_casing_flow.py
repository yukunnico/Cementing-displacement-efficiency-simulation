"""套管内 1D 前沿追踪求解器测试。"""

from __future__ import annotations

import math
import unittest
from typing import cast

from cemdisp.data.fluid_spec import FluidSpec
from cemdisp.data.loaders import load_hu102_tailpipe
from cemdisp.data.pumping_schedule import PumpingSchedule
from cemdisp.data.well_spec import WellSpec
from cemdisp.transport1d.casing_flow import CasingFlowResult
from cemdisp.transport1d import CasingFlowSolver


class CasingFlowSolverTestCase(unittest.TestCase):
    """验证套管内前沿追踪的核心体积推进逻辑。"""

    well_spec: WellSpec = cast(WellSpec, None)
    fluids: tuple[FluidSpec, ...] = cast(tuple[FluidSpec, ...], None)
    schedule: PumpingSchedule = cast(PumpingSchedule, None)
    solver: CasingFlowSolver = cast(CasingFlowSolver, None)
    result: CasingFlowResult = cast(CasingFlowResult, None)

    def setUp(self) -> None:
        super().setUp()
        self.well_spec, self.fluids, self.schedule, _ = load_hu102_tailpipe(include_wash_spacer=False)
        self.solver = CasingFlowSolver()
        self.result = self.solver.run(self.well_spec, self.fluids, self.schedule)

    def test_casing_solver_pipe_area(self) -> None:
        expected_area = math.pi * ((self.well_spec.liner_id_mm or 0.0) / 1000.0) ** 2 / 4.0

        self.assertAlmostEqual(self.result.pipe_cross_section_m2, expected_area, places=8)
        self.assertAlmostEqual(self.result.pipe_cross_section_m2, 0.00918, places=4)

    def test_casing_solver_front_tracking(self) -> None:
        pipe_volume_m3 = self.result.shoe_md_m * self.result.pipe_cross_section_m2
        cement_volume_m3 = self.schedule.steps[0].volume_m3

        # 水泥浆注入量（~16.67 m³）远不足以填满套管到鞋口（管容积 ~71 m³），
        # 因此水泥前缘最终距离 = 水泥体积 / 管截面积 ≈ 1816 m，而非鞋口 7735 m。
        expected_cement_distance = cement_volume_m3 / self.result.pipe_cross_section_m2
        cement_front = self.result.fronts[0]
        displacement_front = self.result.fronts[1]

        self.assertEqual(cement_front.fluid_name, "尾管水泥浆")
        # 水泥浆自身注入量不足以到达鞋口，前缘停在约 1816 m 处。
        self.assertAlmostEqual(cement_front.distance_m, expected_cement_distance, places=1)
        self.assertGreater(pipe_volume_m3, cement_volume_m3)
        self.assertAlmostEqual(expected_cement_distance, 1816.0, delta=20.0)
        # 替浆流体体积大（74 m³ > 管容积 ~71 m³），前缘已穿过鞋口。
        self.assertAlmostEqual(displacement_front.distance_m, self.result.shoe_md_m, places=6)

    def test_casing_solver_shoe_arrival_time(self) -> None:
        pipe_volume_m3 = self.result.shoe_md_m * self.result.pipe_cross_section_m2
        expected_cement_arrival_s = pipe_volume_m3 / self.schedule.steps[0].rate_m3_min * 60.0
        expected_displacement_arrival_s = (
            self.schedule.steps[0].volume_m3 + pipe_volume_m3
        ) / self.schedule.steps[0].rate_m3_min * 60.0

        cement_front = self.result.fronts[0]
        displacement_front = self.result.fronts[1]

        self.assertAlmostEqual(cement_front.time_s, expected_cement_arrival_s, delta=2.0)
        self.assertAlmostEqual(cement_front.time_s, 11272.0, delta=20.0)
        self.assertAlmostEqual(displacement_front.time_s, expected_displacement_arrival_s, delta=2.0)

    def test_pipe_exit_state_before_arrival(self) -> None:
        state = self.solver.pipe_exit_state_at(self.result, 3000.0)

        self.assertAlmostEqual(state.flow_rate_m3_s, 0.378 / 60.0)
        self.assertEqual(state.phase_fractions, (("钻井液", 1.0),))

    def test_pipe_exit_state_after_arrival(self) -> None:
        state = self.solver.pipe_exit_state_at(self.result, 12000.0)

        self.assertAlmostEqual(state.flow_rate_m3_s, 0.378 / 60.0)
        self.assertEqual(state.phase_fractions, (("尾管水泥浆", 1.0),))


if __name__ == "__main__":
    _ = unittest.main()
