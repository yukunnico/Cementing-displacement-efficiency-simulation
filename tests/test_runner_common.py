"""Runner 通用逻辑测试

验证 export_old_vs_new_inlet_comparison 能正确生成 CSV。
"""

import csv
import tempfile
import unittest
from pathlib import Path

from cemdisp.models2d.boundary_bridge import AnnulusInletState
from cemdisp.runners.common import export_old_vs_new_inlet_comparison


def _constant_provider(phase: str, rate: float):
    """返回恒定状态的测试 provider。"""

    def _provider(time_s: float) -> AnnulusInletState:
        return AnnulusInletState(
            time_s=time_s,
            flow_rate_m3_s=rate,
            stage_name="测试阶段",
            phase_fractions=((phase, 1.0),),
        )

    return _provider


class TestExportOldVsNewComparison(unittest.TestCase):
    """测试 export_old_vs_new_inlet_comparison CSV 导出。"""

    def test_csv_created_with_expected_columns(self) -> None:
        """导出的 CSV 应包含 time_s, old_phase, new_phase, old_rate_m3_s, new_rate_m3_s 列。"""

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            old_provider = _constant_provider("mud", 0.01)
            new_provider = _constant_provider("cement", 0.02)
            sample_times = (0.0, 60.0, 120.0)

            path = export_old_vs_new_inlet_comparison(
                output_dir=output_dir,
                old_provider=old_provider,
                new_provider=new_provider,
                sample_times_s=sample_times,
                well_name="测试井",
            )

            self.assertTrue(path.exists())
            rows = list(csv.DictReader(path.read_text(encoding="utf-8-sig").splitlines()))
            self.assertEqual(len(rows), 3)
            self.assertEqual(rows[0]["old_phase"], "mud")
            self.assertEqual(rows[0]["new_phase"], "cement")
            self.assertEqual(rows[0]["old_rate_m3_s"], "0.010000")
            self.assertEqual(rows[0]["new_rate_m3_s"], "0.020000")


if __name__ == "__main__":
    _ = unittest.main()
