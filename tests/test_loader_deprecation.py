"""Loader 废弃警告测试

验证旧的 build_*_annulus_inlet_provider 函数在调用时发出 DeprecationWarning。
"""

import unittest
import warnings

from cemdisp.data.loaders.hu102_loader import build_hu102_annulus_inlet_provider
from cemdisp.data.loaders.hu101_loader import build_hu101_annulus_inlet_provider


class TestLoaderDeprecation(unittest.TestCase):
    """测试 legacy loader provider 发出废弃警告。"""

    def test_hu102_legacy_provider_warns(self) -> None:
        """build_hu102_annulus_inlet_provider 调用时应发出 DeprecationWarning。"""

        from cemdisp.data.loaders.hu102_loader import load_hu102_tailpipe

        _, fluids, schedule, _ = load_hu102_tailpipe()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                _ = build_hu102_annulus_inlet_provider(schedule, fluids)
            except Exception:
                pass  # 即使运行时出错，警告也应在函数入口发出
        self.assertTrue(
            any(item.category is DeprecationWarning for item in caught),
            f"Expected DeprecationWarning, got: {[w.message for w in caught]}",
        )

    def test_hu101_legacy_provider_warns(self) -> None:
        """build_hu101_annulus_inlet_provider 调用时应发出 DeprecationWarning。"""

        from cemdisp.data.loaders.hu101_loader import load_hu101_tailpipe

        _, fluids, schedule, _ = load_hu101_tailpipe()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                _ = build_hu101_annulus_inlet_provider(schedule, fluids)
            except Exception:
                pass
        self.assertTrue(
            any(item.category is DeprecationWarning for item in caught),
            f"Expected DeprecationWarning, got: {[w.message for w in caught]}",
        )


if __name__ == "__main__":
    _ = unittest.main()
