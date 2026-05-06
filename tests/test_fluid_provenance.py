from __future__ import annotations

from typing import cast
import unittest

from cemdisp.data.fluid_provenance import (
    build_injected_fluid_provenance_summary,
    format_injected_fluid_provenance_markdown,
)
from cemdisp.data.loaders import load_hu101_tailpipe, load_hu102_tailpipe
from cemdisp.data.loaders.hu2_loader import load_hu2_tailpipe


class FluidProvenanceSummaryTestCase(unittest.TestCase):
    def test_hu102_summary_marks_cement_as_partial(self) -> None:
        well_spec, fluids, schedule, _ = load_hu102_tailpipe()
        summary = build_injected_fluid_provenance_summary(well_spec.well_name, schedule, fluids)

        items = cast(list[dict[str, str]], summary["明细"])
        details = {item["流体名称"]: item for item in items}
        self.assertFalse(bool(summary["是否全部现场符合"]))
        self.assertEqual(details["尾管水泥浆"]["来源口径"], "部分符合")

    def test_hu2_summary_mentions_field_density_for_plug_and_buffer(self) -> None:
        well_spec, fluids, schedule, _ = load_hu2_tailpipe()
        summary = build_injected_fluid_provenance_summary(well_spec.well_name, schedule, fluids)

        items = cast(list[dict[str, str]], summary["明细"])
        details = {item["流体名称"]: item for item in items}
        self.assertIn("1.90g/cm³", details["压塞液"]["说明"])
        self.assertIn("1.90g/cm³", details["中置液"]["说明"])

    def test_markdown_formatter_outputs_expected_heading(self) -> None:
        well_spec, fluids, schedule, _ = load_hu101_tailpipe()
        summary = build_injected_fluid_provenance_summary(well_spec.well_name, schedule, fluids)
        markdown_lines = format_injected_fluid_provenance_markdown(summary)

        self.assertEqual(markdown_lines[0], "## 注入流体现场符合性检查")
        self.assertTrue(any("轻泥浆" in line for line in markdown_lines))


if __name__ == "__main__":
    unittest.main()
