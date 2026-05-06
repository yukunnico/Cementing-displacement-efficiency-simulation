"""Hu102 尾管段首版环空模型冒烟测试。"""

from __future__ import annotations

from collections.abc import Mapping
import unittest
from typing import cast

from cemdisp.data.fluid_spec import FluidRole
from cemdisp.data.loaders import build_hu102_annulus_inlet_provider, load_hu102_tailpipe
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.runners.hu102_tailpipe import annulus_stop_time_s
from cemdisp.transport1d import CasingFlowSolver


class Hu102InitialModelSmokeTestCase(unittest.TestCase):
    def test_solver_runs_initial_model(self) -> None:
        well_spec, fluids, schedule, _ = load_hu102_tailpipe()
        inlet_provider = build_hu102_annulus_inlet_provider(schedule, fluids)
        solver = AnnulusD2DGASolver(dt=20.0, nz=24, ny=12, total_t=200.0)

        result = solver.run(well_spec, fluids, inlet_provider)

        self.assertEqual(result.well_name, "呼102")
        self.assertFalse(result.metrics.empty)
        self.assertIn("effective_efficiency", result.metrics.columns)
        self.assertNotIn("cbl_quality_proxy", result.metrics.columns)
        self.assertIn("最终结果", result.summary)
        self.assertGreaterEqual(float(result.metrics["effective_efficiency"].iloc[-1]), 0.0)
        self.assertLessEqual(float(result.metrics["effective_efficiency"].iloc[-1]), 1.0)
        self.assertEqual(result.cement_field.shape, (12, 24))
        self.assertEqual(result.wall_field.shape, (12, 24))
        self.assertTrue((result.wall_field == 0.0).all())
        final_summary = cast(Mapping[str, object], result.summary["最终结果"])
        self.assertNotIn("最终质量响应效率", final_summary)

    def test_coupled_model_runs_until_surface_pumping_ends(self) -> None:
        well_spec, fluids, schedule, _ = load_hu102_tailpipe()
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

        inlet_provider = build_coupled_annulus_inlet_provider(casing_result, casing_solver, fluids)
        solver = AnnulusD2DGASolver(dt=4.0, nz=24, ny=12, total_t=stop_time_s)
        result = solver.run(well_spec, fluids, inlet_provider)

        final_time_s = float(result.metrics["time_s"].iloc[-1])
        final_efficiency = float(result.metrics["effective_efficiency"].iloc[-1])
        peak_efficiency = float(result.metrics["effective_efficiency"].max())
        self.assertLess(final_time_s, stop_time_s)
        self.assertLess(final_efficiency, peak_efficiency)
