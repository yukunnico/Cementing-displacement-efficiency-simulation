"""Task 12: hu101 实测居中度剖面可选项测试。"""
import numpy as np
import pytest

from cemdisp.data.loaders import load_hu101_tailpipe


def test_default_standoff_is_assumed_profile():
    """默认（measured_standoff=None）保持 0.38-0.48 名义剖面。"""
    well, _, _, _ = load_hu101_tailpipe()
    vals = [p.value for p in well.standoff_profile]
    assert min(vals) == pytest.approx(0.38)
    assert max(vals) == pytest.approx(0.48)
    assert np.mean(vals) < 0.50


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
