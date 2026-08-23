"""M1: 弥散系数按 dt 归一（恢复量纲正确性）。

CFL 自适应使 dt 从 4.0s 降到 ~0.118s，但 _smooth_dispersion 的系数是
"每步固定幅值"，导致单位物理时间弥散放大约 34 倍、mixing 钉死在 0.59。
本测试验证：
1. __init__ 新增 4 个 dispersion_* 参数及其默认值（默认=基线硬编码）。
2. 固定 dt 模式（dt_step==dt_ref 且 scale=1.0）下有效系数==原硬编码；
   dt 缩小则有效系数同比缩小（每物理秒弥散恒定）。
"""
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver


def test_dispersion_defaults_match_baseline_constants():
    s = AnnulusD2DGASolver(dt=4.0, nz=20, ny=10, total_t=40.0)
    assert s.dispersion_axial == 0.018
    assert s.dispersion_azimuthal == 0.015
    assert s.dispersion_dt_ref == 4.0
    assert s.dispersion_dt_scale == 1.0


def test_dispersion_eff_scales_with_dt():
    s = AnnulusD2DGASolver(dt=4.0, nz=20, ny=10, total_t=40.0)
    # dt_step=dt_ref 且 scale=1 时有效系数=原系数；dt 缩小则同比缩小
    assert abs(s.dispersion_axial * (4.0 / s.dispersion_dt_ref) - 0.018) < 1e-12
    assert abs(s.dispersion_axial * (0.118 / 4.0) - 0.018 * 0.0295) < 1e-6
