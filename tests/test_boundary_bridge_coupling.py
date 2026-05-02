"""1D-2D 耦合边界桥接测试。"""

from __future__ import annotations

import unittest
from collections.abc import Callable
from typing import cast

from cemdisp.data.fluid_spec import FluidSpec
from cemdisp.data.loaders import build_hu102_annulus_inlet_provider, load_hu102_tailpipe
from cemdisp.data.pumping_schedule import PumpingSchedule
from cemdisp.data.well_spec import WellSpec
from cemdisp.models2d.boundary_bridge import AnnulusInletState
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowResult, CasingFlowSolver


class BoundaryBridgeCouplingTestCase(unittest.TestCase):
    """验证套管 1D 出流可正确驱动环空二维入口。"""

    well_spec: WellSpec = cast(WellSpec, None)
    fluids: tuple[FluidSpec, ...] = cast(tuple[FluidSpec, ...], None)
    schedule: PumpingSchedule = cast(PumpingSchedule, None)
    casing_solver: CasingFlowSolver = cast(CasingFlowSolver, None)
    casing_result: CasingFlowResult = cast(CasingFlowResult, None)
    provider: Callable[[float], AnnulusInletState] = cast(Callable[[float], AnnulusInletState], None)

    def setUp(self) -> None:
        super().setUp()
        self.well_spec, self.fluids, self.schedule, _ = load_hu102_tailpipe()
        self.casing_solver = CasingFlowSolver()
        self.casing_result = self.casing_solver.run(self.well_spec, self.fluids, self.schedule)
        self.provider = build_coupled_annulus_inlet_provider(
            self.casing_result,
            self.casing_solver,
            self.fluids,
        )

    def test_coupled_provider_before_shoe_arrival(self) -> None:
        state = self.provider(3000.0)

        self.assertAlmostEqual(state.flow_rate_m3_s, 1.30 / 60.0)
        self.assertEqual(state.phase_fractions, (("mud", 1.0),))

    def test_coupled_provider_after_shoe_arrival(self) -> None:
        state = self.provider(3400.0)

        self.assertAlmostEqual(state.flow_rate_m3_s, 1.30 / 60.0)
        self.assertEqual(state.phase_fractions, (("cement", 1.0),))

    def test_coupled_provider_matches_hardcoded_sustained_tail(self) -> None:
        hardcoded_provider = build_hu102_annulus_inlet_provider(self.schedule, self.fluids, "sustained_tail")

        coupled_state = self.provider(3400.0)
        hardcoded_state = hardcoded_provider(3400.0)

        self.assertEqual(coupled_state.phase_fractions, hardcoded_state.phase_fractions)
        self.assertAlmostEqual(coupled_state.flow_rate_m3_s, hardcoded_state.flow_rate_m3_s)


if __name__ == "__main__":
    _ = unittest.main()
