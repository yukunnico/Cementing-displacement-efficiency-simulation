"""Hu101 尾管段数据加载器测试。"""

from __future__ import annotations

import unittest

from cemdisp.data.fluid_spec import FluidRole, RheologyModel
from cemdisp.data.loaders import build_hu101_annulus_inlet_provider, load_hu101_tailpipe
from cemdisp.data.validation_data import ValidationData
from cemdisp.data.well_spec import WellSpec


class Hu101LoaderTestCase(unittest.TestCase):
    def test_loader_returns_expected_structures(self) -> None:
        well_spec, fluids, schedule, validation_data = load_hu101_tailpipe()

        self.assertIsInstance(well_spec, WellSpec)
        self.assertEqual(well_spec.well_name, "呼101")
        self.assertAlmostEqual(well_spec.top_md_m, 5400.0)
        self.assertAlmostEqual(well_spec.bottom_md_m, 7868.0)
        self.assertAlmostEqual(well_spec.shoe_md_m, 7868.0)
        self.assertAlmostEqual(well_spec.liner_od_mm or 0.0, 139.70)
        self.assertGreaterEqual(len(well_spec.hole_diameter_profile), 5)
        self.assertEqual(len(well_spec.hole_diameter_profile), len(well_spec.inclination_profile))
        self.assertEqual(len(well_spec.hole_diameter_profile), len(well_spec.standoff_profile))
        self.assertEqual(len(well_spec.evaluation_windows), 2)

        fluids_by_role = {fluid.role: fluid for fluid in fluids}
        self.assertEqual(len(fluids), 8)
        self.assertAlmostEqual(fluids_by_role[FluidRole.MUD].density_kg_m3, 1960.0)
        self.assertAlmostEqual(fluids_by_role[FluidRole.WASH].density_kg_m3, 1850.0)
        self.assertAlmostEqual(fluids_by_role[FluidRole.SPACER].density_kg_m3, 2000.0)
        self.assertAlmostEqual(fluids_by_role[FluidRole.LEAD].density_kg_m3, 2100.0)
        self.assertAlmostEqual(fluids_by_role[FluidRole.TAIL].density_kg_m3, 1900.0)
        self.assertEqual(fluids_by_role[FluidRole.MUD].rheology_model, RheologyModel.BINGHAM)
        self.assertEqual(fluids_by_role[FluidRole.LEAD].rheology_model, RheologyModel.POWER_LAW)
        self.assertEqual(fluids_by_role[FluidRole.TAIL].rheology_model, RheologyModel.POWER_LAW)

        self.assertEqual(len(schedule.steps), 9)
        self.assertEqual(schedule.steps[0].fluid_name, "平衡液")
        self.assertEqual(schedule.steps[1].fluid_name, "驱油隔离液")
        self.assertEqual(schedule.steps[2].fluid_name, "领浆")
        self.assertEqual(schedule.steps[3].fluid_name, "尾浆")
        self.assertAlmostEqual(schedule.steps[0].volume_m3, 25.0)
        self.assertAlmostEqual(schedule.steps[2].volume_m3, 47.0)
        self.assertAlmostEqual(schedule.steps[3].volume_m3, 23.0)
        self.assertAlmostEqual(schedule.steps[-1].rate_m3_min, 0.55)
        self.assertAlmostEqual(sum(step.volume_m3 for step in schedule.steps[4:]), 101.4)

        self.assertIsInstance(validation_data, ValidationData)
        self.assertIsNotNone(validation_data.cbl_summary_path)
        self.assertIsNotNone(validation_data.job_report_path)

    def test_annulus_inlet_provider_applies_52m3_shoe_lag(self) -> None:
        _, fluids, schedule, _ = load_hu101_tailpipe()
        provider = build_hu101_annulus_inlet_provider(
            schedule,
            fluids,
            annulus_boundary_mode="field_order_realistic",
            split_cement_phases=True,
        )

        early_state = provider(1.0)
        lead_arrival_time_s = (52.0 + 25.0 + 25.0 + 1.0) / 1.2 * 60.0
        lead_state = provider(lead_arrival_time_s)
        tail_arrival_time_s = (52.0 + 25.0 + 25.0 + 47.0 + 1.0) / 1.2 * 60.0
        tail_state = provider(tail_arrival_time_s)

        self.assertEqual(early_state.phase_fractions, (("mud", 1.0),))
        self.assertIn("鞋口前仍为钻井液", early_state.stage_name)
        self.assertAlmostEqual(early_state.flow_rate_m3_s, 1.2 / 60.0)
        self.assertEqual(lead_state.phase_fractions, (("lead", 1.0),))
        self.assertIn("注领浆", lead_state.stage_name)
        self.assertEqual(tail_state.phase_fractions, (("tail", 1.0),))
        self.assertIn("注尾浆", tail_state.stage_name)

    def test_annulus_inlet_provider_tracks_spacer_phase_before_cement(self) -> None:
        _, fluids, schedule, _ = load_hu101_tailpipe()
        provider = build_hu101_annulus_inlet_provider(
            schedule,
            fluids,
            annulus_boundary_mode="field_order_realistic",
            split_cement_phases=True,
        )

        balance_arrival_time_s = (52.0 + 1.0) / 1.2 * 60.0
        spacer_arrival_time_s = (52.0 + 25.0 + 1.0) / 1.2 * 60.0

        self.assertEqual(provider(balance_arrival_time_s).phase_fractions, (("spacer", 1.0),))
        self.assertIn("注平衡液", provider(balance_arrival_time_s).stage_name)
        self.assertEqual(provider(spacer_arrival_time_s).phase_fractions, (("spacer", 1.0),))
        self.assertIn("注驱油隔离液", provider(spacer_arrival_time_s).stage_name)

    def test_annulus_inlet_provider_switches_to_mud_after_tail_invasion(self) -> None:
        _, fluids, schedule, _ = load_hu101_tailpipe()
        provider = build_hu101_annulus_inlet_provider(
            schedule,
            fluids,
            annulus_boundary_mode="field_order_realistic",
            split_cement_phases=True,
        )

        mud_invasion_time_s = (
            (25.0 + 25.0 + 47.0 + 23.0) / 1.2
            + 2.0 / 0.6
            + 26.0 / 1.5
            + 10.0 / 1.2
            + 15.0 / 1.0
        ) * 60.0
        state = provider(mud_invasion_time_s)

        self.assertEqual(state.phase_fractions, (("mud", 1.0),))
        self.assertIn("注后置液", state.stage_name)

    def test_standoff_profile_range(self) -> None:
        well_spec, _, _, _ = load_hu101_tailpipe()
        values = [point.value for point in well_spec.standoff_profile]

        self.assertTrue(all(0.30 <= value <= 0.82 for value in values))
        self.assertLess(min(values), max(values))
