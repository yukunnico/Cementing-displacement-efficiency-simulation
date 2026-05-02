"""Hu102 尾管段首版环空模型冒烟测试。"""

from __future__ import annotations

import unittest

from cemdisp.data.loaders import build_hu102_annulus_inlet_provider, load_hu102_tailpipe
from cemdisp.models2d import AnnulusD2DGASolver


class Hu102InitialModelSmokeTestCase(unittest.TestCase):
    def test_solver_runs_initial_model(self) -> None:
        well_spec, fluids, schedule, _ = load_hu102_tailpipe()
        inlet_provider = build_hu102_annulus_inlet_provider(schedule, fluids)
        solver = AnnulusD2DGASolver(dt=20.0, nz=24, ny=12, total_t=200.0)

        result = solver.run(well_spec, fluids, inlet_provider)

        self.assertEqual(result.well_name, "呼102")
        self.assertFalse(result.metrics.empty)
        self.assertIn("effective_efficiency", result.metrics.columns)
        self.assertIn("最终结果", result.summary)
        self.assertGreaterEqual(float(result.metrics["effective_efficiency"].iloc[-1]), 0.0)
        self.assertLessEqual(float(result.metrics["effective_efficiency"].iloc[-1]), 1.0)
        self.assertEqual(result.cement_field.shape, (12, 24))
        self.assertEqual(result.wall_field.shape, (12, 24))
