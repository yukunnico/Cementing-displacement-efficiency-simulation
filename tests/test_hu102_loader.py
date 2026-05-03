"""Hu102 尾管段数据加载器测试。"""

from __future__ import annotations

import unittest

from cemdisp.data.fluid_spec import FluidRole, RheologyModel
from cemdisp.data.loaders import build_hu102_annulus_inlet_provider, load_hu102_tailpipe
from cemdisp.data.validation_data import ValidationData
from cemdisp.data.well_spec import WellSpec


class Hu102LoaderTestCase(unittest.TestCase):
    def test_loader_returns_expected_structures(self) -> None:
        well_spec, fluids, schedule, validation_data = load_hu102_tailpipe()

        self.assertIsInstance(well_spec, WellSpec)
        self.assertEqual(well_spec.well_name, "呼102")
        self.assertAlmostEqual(well_spec.top_md_m, 6823.10)
        self.assertAlmostEqual(well_spec.bottom_md_m, 7735.00)
        self.assertAlmostEqual(well_spec.shoe_md_m, 7735.00)
        self.assertAlmostEqual(well_spec.hanger_md_m or 0.0, 6823.10)
        self.assertAlmostEqual(well_spec.casing_id_mm or 0.0, 219.10)
        self.assertAlmostEqual(well_spec.liner_od_mm or 0.0, 139.70)
        self.assertGreaterEqual(len(well_spec.hole_diameter_profile), 60)
        self.assertEqual(len(well_spec.hole_diameter_profile), len(well_spec.inclination_profile))
        self.assertEqual(len(well_spec.hole_diameter_profile), len(well_spec.standoff_profile))
        self.assertEqual(len(well_spec.evaluation_windows), 3)

        mud = next(fluid for fluid in fluids if fluid.role == FluidRole.MUD)
        displacement = next(fluid for fluid in fluids if fluid.role == FluidRole.DISPLACEMENT)
        cement = next(fluid for fluid in fluids if fluid.role == FluidRole.TAIL)
        wash = next(fluid for fluid in fluids if fluid.role == FluidRole.WASH)
        spacer = next(fluid for fluid in fluids if fluid.role == FluidRole.SPACER)
        self.assertEqual(len(fluids), 5)
        self.assertEqual(mud.rheology_model, RheologyModel.BINGHAM)
        self.assertEqual(displacement.rheology_model, RheologyModel.BINGHAM)
        self.assertEqual(cement.rheology_model, RheologyModel.POWER_LAW)
        self.assertEqual(wash.rheology_model, RheologyModel.BINGHAM)
        self.assertEqual(spacer.rheology_model, RheologyModel.BINGHAM)
        self.assertAlmostEqual(mud.density_kg_m3, 2020.0)
        self.assertAlmostEqual(displacement.density_kg_m3, 2020.0)
        self.assertAlmostEqual(cement.density_kg_m3, 2100.0)
        self.assertAlmostEqual(wash.density_kg_m3, 1880.0)
        self.assertAlmostEqual(spacer.density_kg_m3, 1850.0)

        self.assertEqual(len(schedule.steps), 2)
        self.assertAlmostEqual(schedule.steps[0].volume_m3, 16.6666666667, places=2)
        self.assertAlmostEqual(schedule.steps[1].volume_m3, 74.0)
        self.assertEqual(schedule.steps[0].fluid_name, "尾管水泥浆")
        self.assertEqual(schedule.steps[1].fluid_name, "替浆液")
        self.assertAlmostEqual(schedule.steps[0].rate_m3_min, 1.30)
        self.assertAlmostEqual(schedule.steps[1].rate_m3_min, 1.30)

        self.assertIsInstance(validation_data, ValidationData)
        self.assertIsNotNone(validation_data.cbl_summary_path)
        self.assertIsNotNone(validation_data.job_report_path)

    def test_optional_wash_spacer_schedule(self) -> None:
        _, fluids, schedule, _ = load_hu102_tailpipe(include_wash_spacer=True)

        self.assertEqual(len(schedule.steps), 4)
        self.assertEqual(schedule.steps[0].fluid_name, "冲洗液")
        self.assertEqual(schedule.steps[1].fluid_name, "隔离液")
        self.assertEqual(schedule.steps[2].fluid_name, "尾管水泥浆")
        self.assertEqual(schedule.steps[3].fluid_name, "替浆液")
        self.assertAlmostEqual(schedule.steps[0].volume_m3, 3.0)
        self.assertAlmostEqual(schedule.steps[1].volume_m3, 5.0)

        provider = build_hu102_annulus_inlet_provider(schedule, fluids, "tail_then_mud")
        wash_state = provider(1.0)
        spacer_time_s = schedule.steps[0].volume_m3 / schedule.steps[0].rate_m3_min * 60.0 + 1.0
        spacer_state = provider(spacer_time_s)
        displacement_time_s = sum(
            step.volume_m3 / step.rate_m3_min * 60.0 for step in schedule.steps[:3]
        ) + 1.0
        displacement_state = provider(displacement_time_s)

        self.assertEqual(wash_state.phase_fractions, (("spacer", 1.0),))
        self.assertEqual(spacer_state.phase_fractions, (("spacer", 1.0),))
        self.assertEqual(displacement_state.phase_fractions, (("mud", 1.0),))

    def test_standoff_profile_range(self) -> None:
        well_spec, _, _, _ = load_hu102_tailpipe()
        values = [point.value for point in well_spec.standoff_profile]
        self.assertTrue(all(0.30 <= value <= 0.82 for value in values))

    def test_annulus_inlet_provider_modes(self) -> None:
        _, fluids, schedule, _ = load_hu102_tailpipe()

        sustained_provider = build_hu102_annulus_inlet_provider(schedule, fluids, "sustained_tail")
        tail_then_mud_provider = build_hu102_annulus_inlet_provider(schedule, fluids, "tail_then_mud")

        cement_state = sustained_provider(10.0)
        self.assertAlmostEqual(cement_state.flow_rate_m3_s, 1.30 / 60.0)
        self.assertEqual(cement_state.phase_fractions, (("cement", 1.0),))

        displacement_time_s = schedule.steps[0].volume_m3 / schedule.steps[0].rate_m3_min * 60.0 + 1.0
        sustained_state = sustained_provider(displacement_time_s)
        mud_state = tail_then_mud_provider(displacement_time_s)
        self.assertEqual(sustained_state.phase_fractions, (("cement", 1.0),))
        self.assertEqual(mud_state.phase_fractions, (("mud", 1.0),))
