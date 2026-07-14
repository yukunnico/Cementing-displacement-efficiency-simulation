"""Schema tests for paper result table helpers."""

from __future__ import annotations

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from cemdisp.data.loaders.paper import get_paper_metadata
from scripts.paper_data.run_paper_main_results import MAIN_RESULT_COLUMNS, summary_row_from_files


class TestPaperResultSchema(unittest.TestCase):
    """Validate stable paper result schemas without running simulations."""

    def test_summary_row_uses_exact_main_result_columns(self) -> None:
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            well_dir = root / "01_main_results" / "hu101"
            well_dir.mkdir(parents=True)
            summary_path = well_dir / "summary.json"
            summary_path.write_text(
                json.dumps(
                    {
                        "well_id": "hu101",
                        "top_md_m": 5400.0,
                        "bottom_md_m": 7868.0,
                        "annulus_volume_m3": 321.0,
                        "final_effective_efficiency": 0.71,
                        "final_bulk_cement_fill": 0.68,
                        "cbl_eval_interval_efficiency": 0.60,
                        "target_interval_efficiency": 0.75,
                        "channeling_index": 0.10,
                        "mixing_index": 0.20,
                        "instability_index": 0.30,
                    },
                    ensure_ascii=False,
                ),
                encoding="utf-8",
            )
            (well_dir / "time_series.csv").write_text("time_s\n0\n", encoding="utf-8")
            (well_dir / "depth_profile.csv").write_text("井深_m\n5400\n", encoding="utf-8")

            row = summary_row_from_files(summary_path, get_paper_metadata("hu101"), root)

        expected = [
            "well_id",
            "well_name_cn",
            "sample_class",
            "top_md_m",
            "bottom_md_m",
            "annulus_volume_m3",
            "final_effective_efficiency",
            "final_bulk_cement_fill",
            "cbl_eval_interval_efficiency",
            "target_interval_efficiency",
            "channeling_index",
            "mixing_index",
            "instability_index",
            "summary_json_path",
            "time_series_csv_path",
            "depth_profile_csv_path",
        ]
        self.assertEqual(MAIN_RESULT_COLUMNS, expected)
        self.assertEqual(list(row.keys()), expected)
        self.assertEqual(row["well_id"], "hu101")
        self.assertEqual(row["sample_class"], "A_main_validation")
        self.assertEqual(row["summary_json_path"], "01_main_results/hu101/summary.json")


if __name__ == "__main__":
    unittest.main()
