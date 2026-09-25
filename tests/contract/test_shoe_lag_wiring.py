# -*- coding: utf-8 -*-
"""T1（管容死常量）：常量必须真正进入 WellSpec，否则 1D 管容回退到单一内径。

背景（内部自洽加固 Task 3）：`ht1_001` / `hu2` 两口井的"鞋口滞后体积"常量早已
定义且注释完整，但从未传入 `WellSpec`，导致 `CasingFlowSolver._timeline_pipe_volume`
与 `_pipe_cross_section_area` 双双回退到"单一内径 × 鞋深"兜底口径，管容实测偏低
约 −25%（ht1_001）/ −14.8%（hu2），前缘到鞋时刻系统性提前。

本测试是纯 loader 断言（R16）：只查常量是否进入 WellSpec、管容是否大于旧兜底口径，
不涉及数值求解路径。

修复轮 1（评审 Important-1 / R62）：
- 量级阈值的分母是"旧兜底口径"，比值本身随**常量口径**移动，故阈值必须按候选口径的
  下限来定，而不是按 loader 现值。hu2 的 loader 注释并列三套口径（设计 ECD 模拟 77 /
  实际替浆 79 / 组件累计 81.1），原阈值 1.15 在口径取 77 或 79 时会**假红**
  （1.15×69.1245 = 79.49 > 77），因而会误导即将进行的口径裁定；下调到 1.05 后
  三套口径的比值（1.114 / 1.143 / 1.174）全部通过，而"退回单一内径"比值恰为 1.000，
  真实回归依旧必红。ht1_001 侧候选口径（91.7 纯替浆 / 93.7 含压塞液 / 94.5 loader 现值
  / 95 全管内腔）对应比值 1.294–1.340，原阈值 1.20 对最不利候选仍留 7.7% 余量，
  不脆弱，保持不变。
- 量级阈值**不防张冠李戴**：把 ht1_001 的 94.5 误接到 hu2，94.5 > 旧兜底×1.05 照样通过。
  真正的交叉判别力在"等值断言"里（错接时 `well.shoe_lag_volume_m3` 等于**对方**常量、
  不等于本模块常量），故新增 `test_two_wells_do_not_borrow_each_others_constant`
  把这一性质显式固化。
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


def test_two_wells_do_not_borrow_each_others_constant():
    """防张冠李戴：两口井各自引用**自己的**常量（ht1_001 的 94.5 不得落到 hu2）。

    量级阈值抓不到这类错接（94.5 远大于 hu2 旧兜底口径的 1.05 倍），等值断言才是真守卫；
    本测试把"等值断言确实具备交叉判别力"这一性质显式固化，顺带锁定其前提
    （两井常量本身不同——若将来口径裁定使二者相等，本条会红，提示交叉判别已失效）。
    """
    ht1_001_well, _, _, _ = ht1_001_loader.load_ht1_001_tailpipe()
    hu2_well, _, _, _ = hu2_loader.load_hu2_tailpipe()

    # 前提：两井常量不同；否则"等于自己的常量"与"等于对方的常量"不可区分。
    assert hu2_loader.HU2_SHOE_LAG_VOLUME_M3 != ht1_001_loader.HT1_001_SHOE_LAG_VOLUME_M3

    # 错接即红：hu2 拿到 ht1_001 的常量 / ht1_001 拿到 hu2 的常量。
    assert hu2_well.shoe_lag_volume_m3 != ht1_001_loader.HT1_001_SHOE_LAG_VOLUME_M3
    assert ht1_001_well.shoe_lag_volume_m3 != hu2_loader.HU2_SHOE_LAG_VOLUME_M3


def test_tube_volume_larger_than_single_id_fallback():
    """修正后管容必须显著大于旧的单一内径口径（旧值 ht1_001≈70.88 m³）。

    阈值 1.20 的余量按**最不利候选口径**校核：91.7/70.8814 = 1.294（最紧），
    93.7 → 1.322，94.5 → 1.333，95 → 1.340，故 1.20 在任一裁定下都不假红。
    """
    well, _, _, _ = ht1_001_loader.load_ht1_001_tailpipe()
    old = math.pi * (well.liner_id_mm / 2000.0) ** 2 * well.shoe_md_m
    assert well.shoe_lag_volume_m3 > old * 1.20      # 旧口径低约 25%


def test_hu2_tube_volume_larger_than_single_id_fallback():
    """hu2 旧口径低约 14.8%（81.13 vs 69.12）；阈值取 1.05（修复轮 1 由 1.15 下调）。

    1.15 的余量仅 2.06%（81.1296 vs 79.4932），而 loader 注释并列的三套口径里
    77 与 79 都会假红 ⇒ 会误导即将进行的口径裁定。1.05 对最不利候选 77 仍留 6.1%
    余量（77/69.1245 = 1.114），79 → 1.143、81.13 → 1.174 亦通过；而"退回单一内径"
    的比值恰为 1.000，依旧必红。
    """
    well, _, _, _ = hu2_loader.load_hu2_tailpipe()
    old = math.pi * (well.liner_id_mm / 2000.0) ** 2 * well.shoe_md_m
    assert well.shoe_lag_volume_m3 > old * 1.05
