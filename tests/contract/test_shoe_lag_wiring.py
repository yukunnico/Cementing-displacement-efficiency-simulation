# -*- coding: utf-8 -*-
"""T1（管容死常量）：常量必须真正进入 WellSpec，否则 1D 管容回退到单一内径。

背景（内部自洽加固 Task 3）：`ht1_001` / `hu2` 两口井的"鞋口滞后体积"常量早已
定义且注释完整，但从未传入 `WellSpec`，导致 `CasingFlowSolver._timeline_pipe_volume`
与 `_pipe_cross_section_area` 双双回退到"单一内径 × 鞋深"兜底口径，管容实测偏低
约 −25%（ht1_001）/ −14.8%（hu2），前缘到鞋时刻系统性提前。

本测试是纯 loader 断言（R16）：只查常量是否进入 WellSpec、管容是否大于旧兜底口径，
不涉及数值求解路径。
"""
import math

import pytest

from cemdisp.data.loaders import ht1_001_loader, hu2_loader


def test_ht1_001_shoe_lag_volume_wired():
    well, _, _, _ = ht1_001_loader.load_ht1_001_tailpipe()
    assert well.shoe_lag_volume_m3 is not None, "常量未接入 WellSpec（T1 缺陷）"
    assert well.shoe_lag_volume_m3 == pytest.approx(
        ht1_001_loader.HT1_001_SHOE_LAG_VOLUME_M3, rel=1e-9)


def test_hu2_shoe_lag_volume_wired():
    well, _, _, _ = hu2_loader.load_hu2_tailpipe()
    assert well.shoe_lag_volume_m3 is not None, "常量未接入 WellSpec（T1 缺陷）"
    assert well.shoe_lag_volume_m3 == pytest.approx(
        hu2_loader.HU2_SHOE_LAG_VOLUME_M3, rel=1e-9)


def test_tube_volume_larger_than_single_id_fallback():
    """修正后管容必须显著大于旧的单一内径口径（旧值 ht1_001≈70.88 / hu2≈69.12 m³）。"""
    well, _, _, _ = ht1_001_loader.load_ht1_001_tailpipe()
    old = math.pi * (well.liner_id_mm / 2000.0) ** 2 * well.shoe_md_m
    assert well.shoe_lag_volume_m3 > old * 1.20      # 旧口径低约 25%


def test_hu2_tube_volume_larger_than_single_id_fallback():
    """hu2 旧口径低约 14.8%，阈值取 1.15（81.13 > 69.12×1.15 = 79.49）。"""
    well, _, _, _ = hu2_loader.load_hu2_tailpipe()
    old = math.pi * (well.liner_id_mm / 2000.0) ** 2 * well.shoe_md_m
    assert well.shoe_lag_volume_m3 > old * 1.15
