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

Task 12（用户裁定 R67 / R68，2026-09-25）：管容常量口径落地，等值断言同时钉死**字面值**，
不再只与模块常量自比（自比无法发现常量本身被改错）。裁定口径：
- ht1_001 = 94.439m³（严格"地面→鞋口"）= 7.1.4 表四段和 93.4934（止于阻位 7642.674m）
  + 7.1.1"尾浆（下塞）"行 103.326m@9.15L/m（0.9454）⇒ 未取整 94.43884664，按裁定取 94.439。
- hu2 = 80.0858m³ = 纠正三段偏离后的四段和 79.0324（止于阻位 7438.9m）
  + 阻位→鞋 115.125m@9.15L/m（1.0534）⇒ 未取整 80.08580375。
  止于阻位的变体 79.0324 另立常量 `HU2_SHOE_LAG_V_TO_STOP_COLLAR_M3` 并单测。
"""
import math

import pytest

from cemdisp.data.loaders import ht1_001_loader, hu2_loader

# ---- 裁定口径的字面值（评审/用户复核用；与 loader 注释链逐一对应）----
HT1_001_SHOE_LAG_LITERAL_M3 = 94.439          # 93.4934 + 0.9454（未取整 94.43884664）
HU2_V_TO_STOP_COLLAR_LITERAL_M3 = 79.0324     # 设计 7.1.4 四段和，止于阻位 7438.9
HU2_V_TO_SHOE_LITERAL_M3 = 80.0858            # + 阻位→鞋 115.125m@9.15L/m（未取整 80.08580375）

# hu2 设计 7.1.4 表逐行（长度 m, 单位容积 L/m）；独立复算用，不引用 loader 常量。
_HU2_DESIGN_714_SEGMENTS = (
    (3489.0, 12.91),   # 149.2 钻杆 0→3489（等效 ID≈128.2，非 loader 旧用的 129.9）
    (1803.5, 7.42),    # 114.3 钻杆 3489→5292.5（等效 ID≈97.2）
    (1759.8, 9.70),    # 139.7/14.27 尾管 5292.5→7052.3（含 5292.5→5308.564 的 16.064m）
    (386.6, 9.15),     # 139.7/15.88 尾管 7052.3→7438.9（阻位）
)
_HU2_STOP_COLLAR_TO_SHOE_M = 7554.0 - 7438.875  # = 115.125，阻位（作业史 7438.875）→鞋 7554
_HU2_SHOE_UNIT_VOLUME_L_PER_M = 9.15            # 同段 139.7/15.88 尾管单位容积


def test_ht1_001_shoe_lag_volume_wired():
    well, _, _, _ = ht1_001_loader.load_ht1_001_tailpipe()
    assert well.shoe_lag_volume_m3 is not None, "常量未接入 WellSpec（T1 缺陷）"
    assert well.shoe_lag_volume_m3 == pytest.approx(
        ht1_001_loader.HT1_001_SHOE_LAG_VOLUME_M3, rel=1e-9)
    # R67 裁定字面值：不再只与模块常量自比，否则常量本身被改错也看不出来。
    assert well.shoe_lag_volume_m3 == pytest.approx(HT1_001_SHOE_LAG_LITERAL_M3, rel=1e-9)


def test_ht1_001_caliber_chain_reproduces_constant():
    """R67 口径链可复现：7.1.4 四段（止于阻位 7642.674）+ 7.1.1 尾浆（下塞）103.326m。

    四段 = 钻杆149.2 3548.575m@12.91 + 钻杆127 1911.584m@9.11 + 套管168.3 1716.128m@15.15
    + 套管139.7 466.387m@9.15 = 93.49341374；加下塞 103.326m@9.15 = 0.9454329
    ⇒ 94.43884664，按裁定取 94.439（3 位小数）。故此处用 abs 容差校核"取整关系"。
    """
    four_segments_m3 = (
        3548.575 * 12.91 + 1911.584 * 9.11 + 1716.128 * 15.15 + 466.387 * 9.15
    ) / 1000.0
    assert four_segments_m3 == pytest.approx(93.4934, abs=1e-4)
    assert ht1_001_loader.HT1_001_SHOE_LAG_VOLUME_M3 == pytest.approx(
        four_segments_m3 + 103.326 * 9.15 / 1000.0, abs=5e-4)
    assert ht1_001_loader.HT1_001_SHOE_LAG_VOLUME_M3 == pytest.approx(
        HT1_001_SHOE_LAG_LITERAL_M3, rel=1e-9)


def test_hu2_shoe_lag_volume_wired():
    well, _, _, _ = hu2_loader.load_hu2_tailpipe()
    assert well.shoe_lag_volume_m3 is not None, "常量未接入 WellSpec（T1 缺陷）"
    assert well.shoe_lag_volume_m3 == pytest.approx(
        hu2_loader.HU2_SHOE_LAG_VOLUME_M3, rel=1e-9)
    # R68 裁定字面值（严格"地面→鞋口"）。常量是未取整的算术值，故用 rel=1e-6 钉 80.0858。
    assert well.shoe_lag_volume_m3 == pytest.approx(HU2_V_TO_SHOE_LITERAL_M3, rel=1e-6)


def test_hu2_caliber_chain_reproduces_constants():
    """R68 口径链可复现：设计 7.1.4 四段（止于阻位 7438.9）与"止于阻位"变体常量。

    独立复算（不引用 loader 分段），同时锁定"两个候选值"：
    V_to_stop_collar = 79.03241；V_to_shoe = 79.03241 + 115.125m@9.15L/m = 80.08580375。
    """
    stop_collar_m3 = sum(
        length_m * unit_l_per_m for length_m, unit_l_per_m in _HU2_DESIGN_714_SEGMENTS
    ) / 1000.0
    assert stop_collar_m3 == pytest.approx(HU2_V_TO_STOP_COLLAR_LITERAL_M3, rel=1e-6)

    assert hu2_loader.HU2_SHOE_LAG_V_TO_STOP_COLLAR_M3 == pytest.approx(
        stop_collar_m3, rel=1e-9)
    assert hu2_loader.HU2_SHOE_LAG_VOLUME_M3 == pytest.approx(
        stop_collar_m3 + _HU2_STOP_COLLAR_TO_SHOE_M * _HU2_SHOE_UNIT_VOLUME_L_PER_M / 1000.0,
        rel=1e-9)
    assert hu2_loader.HU2_SHOE_LAG_VOLUME_M3 == pytest.approx(
        HU2_V_TO_SHOE_LITERAL_M3, rel=1e-6)


def test_two_wells_do_not_borrow_each_others_constant():
    """防张冠李戴：两口井各自引用**自己的**常量（ht1_001 的 94.439 不得落到 hu2）。

    量级阈值抓不到这类错接（94.439 远大于 hu2 旧兜底口径的 1.05 倍），等值断言才是真守卫；
    本测试把"等值断言确实具备交叉判别力"这一性质显式固化，顺带锁定其前提
    （两井常量本身不同——若将来口径裁定使二者相等，本条会红，提示交叉判别已失效）。
    """
    ht1_001_well, _, _, _ = ht1_001_loader.load_ht1_001_tailpipe()
    hu2_well, _, _, _ = hu2_loader.load_hu2_tailpipe()

    # 前提：两井常量不同（Task 12 裁定后为 94.439 vs 80.0858）；否则"等于自己的常量"
    # 与"等于对方的常量"不可区分。
    assert hu2_loader.HU2_SHOE_LAG_VOLUME_M3 != ht1_001_loader.HT1_001_SHOE_LAG_VOLUME_M3

    # 错接即红：hu2 拿到 ht1_001 的常量 / ht1_001 拿到 hu2 的常量。
    assert hu2_well.shoe_lag_volume_m3 != ht1_001_loader.HT1_001_SHOE_LAG_VOLUME_M3
    assert ht1_001_well.shoe_lag_volume_m3 != hu2_loader.HU2_SHOE_LAG_VOLUME_M3


def test_tube_volume_larger_than_single_id_fallback():
    """修正后管容必须显著大于旧的单一内径口径（旧值 ht1_001≈70.88 m³）。

    阈值 1.20 的余量按**最不利候选口径**校核（Task 12 裁定后常量已定，此处保留候选扫描
    记录以便追溯）：91.7/70.8814 = 1.294（最紧），93.7 → 1.322，94.439 → 1.333，
    95 → 1.340，故 1.20 在任一裁定下都不假红。
    """
    well, _, _, _ = ht1_001_loader.load_ht1_001_tailpipe()
    old = math.pi * (well.liner_id_mm / 2000.0) ** 2 * well.shoe_md_m
    assert well.shoe_lag_volume_m3 > old * 1.20      # 旧口径低约 25%


def test_hu2_tube_volume_larger_than_single_id_fallback():
    """hu2 旧口径低约 13.7%（Task 12 裁定后 80.0858 vs 69.1245）；阈值取 1.05（修复轮 1 下调）。

    1.15 的余量仅 2.06%（口径裁定前 81.1296 vs 79.4932），而 loader 注释并列的三套口径里
    77 与 79 都会假红 ⇒ 会误导口径裁定。1.05 对最不利候选 77 仍留 6.1% 余量
    （77/69.1245 = 1.114），79 → 1.143、80.0858 → 1.159 亦通过；而"退回单一内径"
    的比值恰为 1.000，依旧必红。
    """
    well, _, _, _ = hu2_loader.load_hu2_tailpipe()
    old = math.pi * (well.liner_id_mm / 2000.0) ** 2 * well.shoe_md_m
    assert well.shoe_lag_volume_m3 > old * 1.05
