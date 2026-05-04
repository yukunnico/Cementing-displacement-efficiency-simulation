from __future__ import annotations

import unittest

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cemdisp.models2d import AnnulusSimulationResult
from cemdisp.reporting.contour_plots import plot_annulus_snapshots, plot_final_fields_contour


class ReportingOrientationTestCase(unittest.TestCase):
    def _make_result(self) -> AnnulusSimulationResult:
        cement = np.array([[1.0, 0.8], [0.2, 0.0]], dtype=float)
        spacer = np.array([[0.0, 0.1], [0.3, 0.4]], dtype=float)
        wall = np.zeros_like(cement)
        return AnnulusSimulationResult(
            well_name="测试井",
            geom={"md": np.array([100.0, 90.0], dtype=float), "y": np.array([0.0, 1.0], dtype=float)},
            cement_field=cement,
            spacer_field=spacer,
            wall_field=wall,
            metrics=pd.DataFrame(),
            depth_profiles=pd.DataFrame(),
            summary={},
            snapshot_times_s=(0.0,),
            cement_snapshots=(cement.copy(),),
            spacer_snapshots=(spacer.copy(),),
            wall_snapshots=(wall.copy(),),
        )

    def test_snapshot_plot_keeps_wide_to_narrow_row_order(self) -> None:
        result = self._make_result()
        figure = plot_annulus_snapshots(result, output_dir=None, n_panels=1)
        cement_image_array = np.asarray(figure.axes[0].images[0].get_array())
        spacer_image_array = np.asarray(figure.axes[1].images[0].get_array())

        np.testing.assert_allclose(cement_image_array, result.cement_snapshots[0])
        np.testing.assert_allclose(spacer_image_array, result.spacer_snapshots[0])
        plt.close(figure)

    def test_final_field_plot_keeps_wide_to_narrow_row_order(self) -> None:
        result = self._make_result()
        figure = plot_final_fields_contour(result, output_dir=None)
        cement_image_array = np.asarray(figure.axes[0].images[0].get_array())
        spacer_image_array = np.asarray(figure.axes[1].images[0].get_array())
        effective_image_array = np.asarray(figure.axes[2].images[0].get_array())

        np.testing.assert_allclose(cement_image_array, result.cement_field)
        np.testing.assert_allclose(spacer_image_array, result.spacer_field)
        np.testing.assert_allclose(effective_image_array, result.cement_field * (1.0 - result.wall_field))
        plt.close(figure)


if __name__ == "__main__":
    unittest.main()
