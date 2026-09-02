"""2026-09-02 管容链现场核实修复测试（0708 原件取证，TDD 先行）。

背景：casing_flow._timeline_pipe_volume 与 _pipe_cross_section_area 双链均优先消费
well_spec.shoe_lag_volume_m3（cemdisp/transport1d/casing_flow.py:706-753），直接设置
该字段可让两条管容链同时统一到现场真值。

三井修复口径（0708 原件核实）：
- hu103: 88.55 m³（20315 灌水+流量计双证实测；20313.txt 行490 设计表 90.2 m³ 并存注记）
- hu101: 102.6 m³（0708 分段真实内径链 149.2 钻杆 ID129.9 0-5397.21m + 尾管两段；
  2011114.txt 行157"理论需102.2方碰压"/行158"累计泵冲到量碰压…=101.8方"呼应）；
  liner_id=91.73（由 52m³ 无出处滞后反推，2026-08-29 查证 missing）退役，
  改用 139.7mm 段厚壁真实内径 108.10mm（=139.7-2×15.8，与壁厚常量自洽）。
- hu102: 88.9 m³（149.2 钻杆 0-5101.25 + 114.3 钻杆 ID97.18 5101.25-6819.37 + 尾管链）；
  实际替浆仅 74 m³（泵冲 72，单流阀失效、到量未碰压——现场真实），分歧注记保留。
其余 5 井（hu1/hu2/ht1_001/ht1_003/ht1_004）不动。
"""

import pytest

from cemdisp.data.loaders import (
    load_hu101_tailpipe,
    load_hu102_tailpipe,
    load_hu103_tailpipe,
    load_hu103_tailpipe_actual,
)
from cemdisp.transport1d.casing_flow import CasingFlowSolver


# ---------------------------------------------------------------------------
# hu103（P0）：shoe_lag_volume_m3 = 88.55
# ---------------------------------------------------------------------------


def test_hu103_shoe_lag_volume_matches_field_chain():
    """hu103 全井管容实测链 88.55 m³（20315 灌水+流量计双证），注记双口径出处。"""
    well, _, _, _ = load_hu103_tailpipe()
    assert well.shoe_lag_volume_m3 == pytest.approx(88.55)
    notes_text = "\n".join(well.notes)
    assert "88.55" in notes_text, "notes 应含实测值 88.55"
    assert "20315" in notes_text, "notes 应注记 20315（灌水+流量计双证）出处"
    assert "20313" in notes_text, "notes 应注记 20313（设计表 90.2/90.7）并存口径"
    # 实际版共用同一 WellSpec
    well_actual, _, _, _ = load_hu103_tailpipe_actual()
    assert well_actual.shoe_lag_volume_m3 == pytest.approx(88.55)


def test_hu103_pipe_capacity_chain_uses_field_volume():
    """hu103 双链（截面积+时间轴迟到体积）统一由 88.55 驱动。"""
    well, fluids, schedule, _ = load_hu103_tailpipe()
    solver = CasingFlowSolver(enable_gravity=False)
    result = solver.run(well, fluids, schedule)
    # _pipe_cross_section_area 优先级：shoe_lag_volume_m3 / shoe_md_m
    assert result.pipe_cross_section_m2 == pytest.approx(88.55 / well.shoe_md_m)


def test_hu103_cement_end_before_pumping_end():
    """hu103 尾浆尾缘 146+88.55=234.55 < 总注入 236.2 → 在泵注结束（碰压）前过鞋口。

    设计版：前段 1.6m³/min（63m³）+ 后段 1.8m³/min。
    尾浆尾缘过鞋口时刻 = 63/1.6 + 83/1.8 + 88.55/1.8 = 134.68min ≈ 8081s；
    泵注结束 = 63/1.6 + 173.2/1.8 = 135.60min ≈ 8136s；差约 55s（≈碰压时刻）。
    """
    well, fluids, schedule, _ = load_hu103_tailpipe()
    solver = CasingFlowSolver(enable_gravity=False)
    result = solver.run(well, fluids, schedule)
    assert result.cement_end_time_s < result.pumping_end_time_s
    assert result.pumping_end_time_s - result.cement_end_time_s == pytest.approx(55.0, abs=90.0)
    assert result.cement_end_time_s == pytest.approx(8080.8, abs=120.0)


# ---------------------------------------------------------------------------
# hu101（P1）：shoe_lag_volume_m3 = 102.6，liner_id 91.73 退役
# ---------------------------------------------------------------------------


def test_hu101_shoe_lag_volume_matches_field_chain():
    """hu101 真实管容链 ≈102.6 m³（2011114 行157 理论碰压 102.2/行158 实泵 101.8 呼应）。"""
    well, _, _, _ = load_hu101_tailpipe()
    assert well.shoe_lag_volume_m3 == pytest.approx(102.6)
    notes_text = "\n".join(well.notes)
    assert "102.6" in notes_text, "notes 应含真值 102.6 与出处注记"
    assert "2011114" in notes_text, "notes 应注记 2011114 出处"


def test_hu101_liner_id_legacy_inference_retired():
    """liner_id=91.73（52m³ 无出处滞后反推）退役 → 139.7mm 段厚壁真实内径 108.10。"""
    from cemdisp.data.loaders.hu101_loader import (
        HU101_LINER_ID_MM,
        HU101_LINER_WALL_THICKNESS_MM,
        HU101_LOWER_LINER_OD_MM,
        HU101_SHOE_LAG_VOLUME_M3,
    )

    assert HU101_SHOE_LAG_VOLUME_M3 == pytest.approx(102.6)
    assert HU101_LINER_ID_MM == pytest.approx(108.10)
    # 与壁厚常量自洽：139.7 - 2×15.8
    assert HU101_LINER_ID_MM == pytest.approx(
        HU101_LOWER_LINER_OD_MM - 2.0 * HU101_LINER_WALL_THICKNESS_MM
    )
    assert HU101_LINER_ID_MM != pytest.approx(91.73)
    well, _, _, _ = load_hu101_tailpipe()
    assert well.liner_id_mm == pytest.approx(108.10)


def test_hu101_cement_end_at_pumping_end():
    """hu101 尾浆尾缘 120+102.6=222.6 > 总注入 221.4 → 尾缘到不了鞋口，停泵末兜底。"""
    well, fluids, schedule, _ = load_hu101_tailpipe()
    solver = CasingFlowSolver(enable_gravity=False)
    result = solver.run(well, fluids, schedule)
    assert result.cement_end_time_s == pytest.approx(result.pumping_end_time_s)


# ---------------------------------------------------------------------------
# hu102（P2）：shoe_lag_volume_m3 = 88.9（无送入段链 → 有送入段链）
# ---------------------------------------------------------------------------


def test_hu102_shoe_lag_volume_matches_field_chain():
    """hu102 真实管容链 ≈88.9 m³；74/72（实际替浆/泵冲到量）与 88.9 分歧注记保留。"""
    well, _, _, _ = load_hu102_tailpipe()
    assert well.shoe_lag_volume_m3 == pytest.approx(88.9)
    notes_text = "\n".join(well.notes)
    assert "88.9" in notes_text, "notes 应含真值 88.9 与链组成注记"
    assert "74" in notes_text and "72" in notes_text, (
        "notes 应注记实际替浆 74/泵冲 72 与 88.9 的分歧（单流阀失效、到量未碰压）"
    )


def test_hu102_displacement_front_cannot_reach_shoe():
    """hu102 尾浆尾缘 52+88.9=140.9 > 替浆步累计 139 → 替浆步内尾缘到不了鞋口。

    schedule 末尾的"循环排混浆"（次日后处理，RESTART，41m³）继续注入，
    尾浆尾缘 140.9 将落在该步时间窗内——模型按步骤串行如实反映；
    本测试只锁定 0708 核实事实：替浆到量（139m³）时尾浆尾缘仍未过鞋口。
    """
    well, fluids, schedule, _ = load_hu102_tailpipe()
    solver = CasingFlowSolver(enable_gravity=False)
    result = solver.run(well, fluids, schedule)
    displacement_end_cum = sum(s.volume_m3 for s in schedule.steps[:8])
    # 替浆步结束时累计 139 < 尾浆尾缘过鞋口所需 52+88.9=140.9
    assert displacement_end_cum == pytest.approx(139.0)
    assert 52.0 + well.shoe_lag_volume_m3 > displacement_end_cum
    # 尾浆尾缘过鞋口时刻落在循环排混浆步内：晚于替浆步结束时刻
    circulation_volume = schedule.steps[8].volume_m3
    circulation_rate_m3_s = schedule.steps[8].rate_m3_min / 60.0
    circulation_start_s = result.pumping_end_time_s - circulation_volume / circulation_rate_m3_s
    assert result.cement_end_time_s > circulation_start_s
