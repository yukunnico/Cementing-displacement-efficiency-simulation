"""Tests for paper-only loader wrappers."""

from __future__ import annotations

import unittest

from cemdisp.data.fluid_spec import FluidSpec
from cemdisp.data.loaders.paper import (
    get_paper_metadata,
    iter_paper_wells,
    load_paper_well,
)
from cemdisp.data.pumping_schedule import PumpingSchedule
from cemdisp.data.validation_data import ValidationData
from cemdisp.data.well_spec import WellSpec


class TestPaperLoaders(unittest.TestCase):
    """Validate paper loader contract and metadata."""

    def test_registry_returns_all_wells_when_pending_included(self) -> None:
        records = tuple(iter_paper_wells(include_pending=True))
        self.assertEqual([record.well_id for record in records], ["hu101", "hu102", "hu103", "ht1_001"])

    def test_default_registry_excludes_pending_validation_well(self) -> None:
        records = tuple(iter_paper_wells())
        self.assertEqual([record.well_id for record in records], ["hu101", "hu102", "hu103"])

    def test_loaders_keep_legacy_four_tuple_contract(self) -> None:
        for well_id in ("hu101", "hu102", "hu103", "ht1_001"):
            with self.subTest(well_id=well_id):
                well_spec, fluids, schedule, validation_data = load_paper_well(well_id)
                metadata = get_paper_metadata(well_id)

                self.assertIsInstance(well_spec, WellSpec)
                self.assertTrue(all(isinstance(fluid, FluidSpec) for fluid in fluids))
                self.assertIsInstance(schedule, PumpingSchedule)
                self.assertIsInstance(validation_data, ValidationData)
                self.assertEqual(metadata.well_id, well_id)
                self.assertEqual(metadata.paper_version, "pead_v1")
                self.assertGreater(schedule.total_volume_m3, 0.0)

                for window in well_spec.evaluation_windows:
                    self.assertGreaterEqual(window.top_md_m, well_spec.top_md_m)
                    self.assertLessEqual(window.bottom_md_m, well_spec.bottom_md_m)

    def test_main_validation_wells_have_cbl_pass_rate(self) -> None:
        expected = {"hu101": 0.6277, "hu102": 0.6665, "hu103": 0.1206}
        for well_id, measured in expected.items():
            with self.subTest(well_id=well_id):
                _, _, _, validation_data = load_paper_well(well_id)
                metadata = get_paper_metadata(well_id)
                self.assertEqual(metadata.sample_class, "A_main_validation")
                self.assertTrue(metadata.include_in_cbl_metrics)
                self.assertAlmostEqual(metadata.cbl_pass_rate or -1.0, measured, places=4)
                self.assertAlmostEqual(validation_data.cbl_pass_rate or -1.0, measured, places=4)

    def test_ht1_001_is_pending_application_sample(self) -> None:
        metadata = get_paper_metadata("ht1_001")
        self.assertEqual(metadata.sample_class, "B_application_pending_validation")
        self.assertFalse(metadata.include_in_cbl_metrics)
        self.assertIsNone(metadata.cbl_pass_rate)


if __name__ == "__main__":
    unittest.main()
