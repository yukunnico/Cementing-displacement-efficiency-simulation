# STATUS: history —— 本文件锁定“重构前”的历史行为契约。
# 源模型口径重构（2026-09-14）后，相关断言预期变红；红灯不等同回归失败，
# 需按 docs/superpowers/plans/2026-09-14-d2dga-source-fidelity-and-repo-cleanup.md 逐条复核。
"""Task 12: hu101 实测居中度剖面可选项测试。"""
import pytest

from cemdisp.data.loaders import load_hu101_tailpipe


def test_default_standoff_is_assumed_profile():
    """默认（measured_standoff=None）剖面。

    ⚠️ 2026-09-11 用户指令临时改动：默认由 legacy 名义剖面 0.38–0.48（均 0.429）
    改为全井常数 0.80（反推情景，与现场"居中度下降"记录方向相反）。
    LEGACY 断言：min==0.38 / max==0.48 / mean<0.50。
    """
    well, _, _, _ = load_hu101_tailpipe()
    vals = [p.value for p in well.standoff_profile]
    assert min(vals) == pytest.approx(0.80)
    assert max(vals) == pytest.approx(0.80)


def test_measured_between_centralizers_profile():
    """'between_centralizers' 用扶正器间实测剖面（偏下限，含 0.22 低值）。"""
    well, _, _, _ = load_hu101_tailpipe(measured_standoff="between_centralizers")
    vals = [p.value for p in well.standoff_profile]
    assert min(vals) == pytest.approx(0.22)
    assert max(vals) == pytest.approx(0.78)
    assert len(vals) == 10


def test_measured_at_centralizers_profile():
    """'at_centralizers' 用扶正器处实测剖面（偏上限，全 >=0.60）。"""
    well, _, _, _ = load_hu101_tailpipe(measured_standoff="at_centralizers")
    vals = [p.value for p in well.standoff_profile]
    assert min(vals) >= 0.60
    assert max(vals) == pytest.approx(0.88)
    assert len(vals) == 10


def test_invalid_measured_standoff_raises():
    with pytest.raises(ValueError, match="measured_standoff"):
        load_hu101_tailpipe(measured_standoff="bogus")
