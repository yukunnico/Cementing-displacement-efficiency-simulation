"""Hu101 尾管段 1D-2D 停止口径测试。"""

from __future__ import annotations

from collections.abc import Mapping
import unittest
from typing import cast

from cemdisp.data.fluid_spec import FluidRole
from cemdisp.data.loaders import load_hu101_tailpipe
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.runners.hu101_tailpipe import annulus_stop_time_s
from cemdisp.transport1d import CasingFlowSolver


class Hu101InitialModelStopTestCase(unittest.TestCase):
    def test_coupled_model_runs_until_surface_pumping_ends(self) -> None:
        well_spec, fluids, schedule, _ = load_hu101_tailpipe()
        casing_solver = CasingFlowSolver(enable_gravity=True)
        casing_result = casing_solver.run(well_spec, fluids, schedule)
        stop_time_s = annulus_stop_time_s(casing_result=casing_result, fluids=fluids)
        displacement_front_time_s = next(
            front.time_s
            for front in casing_result.fronts
            if next(fluid.role for fluid in fluids if fluid.name == front.fluid_name) == FluidRole.DISPLACEMENT
        )
        self.assertAlmostEqual(stop_time_s, casing_result.pumping_end_time_s, places=6)
        self.assertGreaterEqual(stop_time_s, displacement_front_time_s)

        inlet_provider = build_coupled_annulus_inlet_provider(
            casing_result,
            casing_solver,
            fluids,
            split_cement_phases=True,
        )
        solver = AnnulusD2DGASolver(dt=4.0, nz=24, ny=12, total_t=stop_time_s)
        result = solver.run(well_spec, fluids, inlet_provider)

        final_time_s = float(result.metrics["time_s"].iloc[-1])
        self.assertLess(final_time_s, stop_time_s)
        self.assertNotIn("cbl_quality_proxy", result.metrics.columns)
        final_summary = cast(Mapping[str, object], result.summary["最终结果"])
        self.assertNotIn("最终质量响应效率", final_summary)
