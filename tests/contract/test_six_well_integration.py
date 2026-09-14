"""六井集成测试

验证所有六口井的 loader 都能正常生成 ShoeTimeline 与同步画像卡，
确保 boundary-synchronization 链条端到端可运行。
"""

import tempfile
import unittest
from pathlib import Path

from cemdisp.data.loaders import (
    export_hu101_sync_card_markdown,
    export_hu102_sync_card_markdown,
    export_hu103_sync_card_markdown,
    export_hu1_sync_card_markdown,
    export_hu2_sync_card_markdown,
    export_ht1_001_sync_card_markdown,
)


class TestSixWellIntegration(unittest.TestCase):
    """六口井同步画像卡端到端集成测试。"""

    def test_all_six_wells_export_sync_card(self) -> None:
        """每口井的 export_*_sync_card_markdown 都能生成非空 Markdown 文件。"""

        exporters = [
            ("呼101", export_hu101_sync_card_markdown),
            ("呼102", export_hu102_sync_card_markdown),
            ("呼103", export_hu103_sync_card_markdown),
            ("呼探1井", export_hu1_sync_card_markdown),
            ("呼探1-002井（HT1-002）", export_hu2_sync_card_markdown),
            ("呼探1-001井（HT1-001）", export_ht1_001_sync_card_markdown),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            for well_name, exporter in exporters:
                with self.subTest(well=well_name):
                    path = exporter(output_dir)
                    self.assertTrue(path.exists(), f"{well_name} sync card not created")
                    content = path.read_text(encoding="utf-8")
                    self.assertIn(well_name, content)
                    self.assertIn("鞋口同步事件数", content)


if __name__ == "__main__":
    _ = unittest.main()
