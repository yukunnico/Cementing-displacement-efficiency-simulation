"""
出口边界条件测试。

验证 AnnulusD2DGASolver 的 open_outlet 参数行为：
- open_outlet=True（默认）：_limit_phase_volume 不限制体积，仅裁剪到 [0, 1]
- open_outlet=False：_limit_phase_volume 按累计入环空体积限制场量
"""

import unittest
import numpy as np

from cemdisp.models2d.annulus_d2dga import (
    AnnulusD2DGASolver,
    _limit_phase_volume,
    _trapez2d,
)


class TestLimitPhaseVolume(unittest.TestCase):
    """测试 _limit_phase_volume 函数的开放/封闭出口行为。"""

    def _make_geom(self, ny: int = 5, nz: int = 10) -> dict:
        """构造最小几何参数字典。"""
        y = np.linspace(0.0, 1.0, ny)
        s = np.linspace(0.0, 100.0, nz)
        b = np.ones((ny, nz)) * 0.02
        return {"y": y, "s": s, "b": b}

    def test_open_outlet_clips_to_unit_interval(self) -> None:
        """开放出口时，仅裁剪到 [0, 1]，不缩放场量。"""
        geom = self._make_geom()
        field = np.full((5, 10), 0.8)
        result = _limit_phase_volume(field, geom, target_volume_m3=0.0, open_outlet=True)
        np.testing.assert_array_equal(result, field)

    def test_open_outlet_ignores_target_volume(self) -> None:
        """开放出口时，即使 target_volume=0 也不清零场量。"""
        geom = self._make_geom()
        field = np.full((5, 10), 0.5)
        result = _limit_phase_volume(field, geom, target_volume_m3=0.0, open_outlet=True)
        np.testing.assert_array_equal(result, field)

    def test_open_outlet_clips_negative_values(self) -> None:
        """开放出口时，负值被裁剪到 0。"""
        geom = self._make_geom()
        field = np.full((5, 10), -0.1)
        result = _limit_phase_volume(field, geom, target_volume_m3=1.0, open_outlet=True)
        np.testing.assert_array_equal(result, np.zeros((5, 10)))

    def test_open_outlet_clips_values_above_one(self) -> None:
        """开放出口时，大于 1 的值被裁剪到 1。"""
        geom = self._make_geom()
        field = np.full((5, 10), 1.5)
        result = _limit_phase_volume(field, geom, target_volume_m3=1.0, open_outlet=True)
        np.testing.assert_array_equal(result, np.ones((5, 10)))

    def test_closed_outlet_zero_target_clears_field(self) -> None:
        """封闭出口且 target_volume=0 时，清零场量。"""
        geom = self._make_geom()
        field = np.full((5, 10), 0.5)
        result = _limit_phase_volume(field, geom, target_volume_m3=0.0, open_outlet=False)
        np.testing.assert_array_equal(result, np.zeros((5, 10)))

    def test_closed_outlet_scales_down_when_overfilled(self) -> None:
        """封闭出口且当前体积超过目标时，按比例缩小。"""
        geom = self._make_geom()
        field = np.full((5, 10), 0.8)
        target = 0.001  # 远小于当前体积
        result = _limit_phase_volume(field, geom, target_volume_m3=target, open_outlet=False)
        # 结果应小于原始值
        self.assertTrue(np.all(result <= field))
        # 结果非负
        self.assertTrue(np.all(result >= 0.0))

    def test_closed_outlet_no_change_when_within_limit(self) -> None:
        """封闭出口且当前体积未超目标时，不做缩放。"""
        geom = self._make_geom()
        field = np.full((5, 10), 0.001)  # 非常小的值
        target = 100.0  # 远大于当前体积
        result = _limit_phase_volume(field, geom, target_volume_m3=target, open_outlet=False)
        np.testing.assert_allclose(result, field, rtol=1.0e-6)

    def test_default_open_outlet_is_false(self) -> None:
        """_limit_phase_volume 默认 open_outlet=False（与函数签名一致）。"""
        geom = self._make_geom()
        field = np.full((5, 10), 0.5)
        result_default = _limit_phase_volume(field, geom, target_volume_m3=0.0)
        result_explicit = _limit_phase_volume(field, geom, target_volume_m3=0.0, open_outlet=False)
        np.testing.assert_array_equal(result_default, result_explicit)


class TestAnnulusD2DGASolverOpenOutlet(unittest.TestCase):
    """测试 AnnulusD2DGASolver 的 open_outlet 参数。"""

    def test_default_open_outlet_is_true(self) -> None:
        """默认 open_outlet=True（与 __init__ 签名一致）。"""
        solver = AnnulusD2DGASolver()
        self.assertTrue(solver.open_outlet)

    def test_open_outlet_false_explicit(self) -> None:
        """显式设置 open_outlet=False 时正确保存。"""
        solver = AnnulusD2DGASolver(open_outlet=False)
        self.assertFalse(solver.open_outlet)

    def test_open_outlet_true_explicit(self) -> None:
        """显式设置 open_outlet=True 时正确保存。"""
        solver = AnnulusD2DGASolver(open_outlet=True)
        self.assertTrue(solver.open_outlet)

    def test_backward_compatibility_no_open_outlet_arg(self) -> None:
        """不传 open_outlet 参数时，行为与 open_outlet=True 一致（向后兼容）。"""
        solver_default = AnnulusD2DGASolver()
        solver_explicit = AnnulusD2DGASolver(open_outlet=True)
        self.assertEqual(solver_default.open_outlet, solver_explicit.open_outlet)
        self.assertEqual(solver_default.dt, solver_explicit.dt)
        self.assertEqual(solver_default.nz, solver_explicit.nz)
        self.assertEqual(solver_default.ny, solver_explicit.ny)


if __name__ == "__main__":
    unittest.main()
