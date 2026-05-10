"""Sync-card 导出测试

验证各井 loader 的 export_*_sync_card_markdown 函数能正确生成 Markdown 文件。
"""

import tempfile
import unittest
from pathlib import Path

from cemdisp.data.loaders.hu102_loader import export_hu102_sync_card_markdown


class TestSyncCardExport(unittest.TestCase):
    """测试同步画像卡 Markdown 导出。"""

    def test_hu102_sync_card_md_created(self) -> None:
        """export_hu102_sync_card_markdown 应生成非空 Markdown 文件。"""

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            output_path = export_hu102_sync_card_markdown(output_dir)

            self.assertTrue(output_path.exists())
            content = output_path.read_text(encoding="utf-8")
            self.assertIn("# 呼102 同步画像卡", content)
            self.assertIn("鞋口同步事件数", content)
            self.assertIn("井名", content)


if __name__ == "__main__":
    _ = unittest.main()
