"""Tests for paper CBL validation metrics."""

from __future__ import annotations

import unittest

from scripts.paper_data.build_cbl_validation_tables import compute_cbl_error_metrics


class TestCBLMetrics(unittest.TestCase):
    """Validate CBL error formulae used for paper tables."""

    def test_compute_absolute_and_relative_error(self) -> None:
        row = compute_cbl_error_metrics(
            "hu101",
            well_name_cn="呼101",
            predicted=0.60,
            measured=0.50,
            included=True,
            notes="x",
        )
        self.assertAlmostEqual(row["absolute_error"], 0.10)
        self.assertAlmostEqual(row["relative_error_percent"], 20.0)
        self.assertTrue(row["included_in_metrics"])

    def test_missing_measured_value_leaves_error_blank(self) -> None:
        row = compute_cbl_error_metrics(
            "ht1_001",
            well_name_cn="呼探1-001井（HT1-001）",
            predicted=0.70,
            measured=None,
            included=False,
            notes="pending",
        )
        self.assertEqual(row["measured_cbl_pass_rate"], "")
        self.assertEqual(row["absolute_error"], "")
        self.assertFalse(row["included_in_metrics"])


if __name__ == "__main__":
    unittest.main()
