# -*- coding: utf-8 -*-
"""8 井验证窗口（EvaluationWindow）一致性与现场值校验。

窗口四口径约定见 cemdisp/data/well_spec.py EvaluationWindow docstring
（核对汇总 §10.2，2026-08-17 全仓统一）：
cbl / cbl_quality / cbl_digitization / formation_target / model_focus / oil_gas_show / custom。
"""

from __future__ import annotations

import pytest

from cemdisp.data.loaders.hu101_loader import load_hu101_tailpipe
from cemdisp.data.loaders.hu102_loader import load_hu102_tailpipe
from cemdisp.data.loaders.hu103_loader import load_hu103_tailpipe
from cemdisp.data.loaders.hu1_loader import load_hu1_tailpipe
from cemdisp.data.loaders.hu2_loader import load_hu2_tailpipe
from cemdisp.data.loaders.ht1_001_loader import load_ht1_001_tailpipe
from cemdisp.data.loaders.ht1_003_loader import load_ht1_003_tailpipe
from cemdisp.data.loaders.ht1_004_loader import load_ht1_004_tailpipe

ALLOWED_WINDOW_TYPES = {
    "cbl",
    "cbl_quality",
    "cbl_digitization",
    "formation_target",
    "model_focus",
    "oil_gas_show",
    "custom",
}

# 井名 -> 默认加载入口
LOADERS = {
    "hu101": load_hu101_tailpipe,
    "hu102": load_hu102_tailpipe,
    "hu103": load_hu103_tailpipe,
    "hu1": load_hu1_tailpipe,
    "hu2": load_hu2_tailpipe,
    "ht1_001": load_ht1_001_tailpipe,
    "ht1_003": load_ht1_003_tailpipe,
    "ht1_004": load_ht1_004_tailpipe,
}


def _windows(well: str):
    well_spec, _fluids, _schedule, _validation = LOADERS[well]()
    return well_spec, list(well_spec.evaluation_windows)


def test_window_types_unified() -> None:
    """全仓不再使用旧 'target' 标签，且取值均在约定集合内。"""
    for well in LOADERS:
        _, windows = _windows(well)
        for w in windows:
            assert w.window_type in ALLOWED_WINDOW_TYPES, (
                f"{well} 窗口 {w.name} 的 window_type={w.window_type!r} 不在约定集合"
            )
            assert w.window_type != "target"


def test_windows_within_model_domain() -> None:
    """所有窗口必须落在模型域 [top_md_m, bottom_md_m] 内。"""
    for well in LOADERS:
        well_spec, windows = _windows(well)
        for w in windows:
            assert w.top_md_m >= well_spec.top_md_m - 1e-6, f"{well} 窗口 {w.name} 顶深越界"
            assert w.bottom_md_m <= well_spec.bottom_md_m + 1e-6, f"{well} 窗口 {w.name} 底深越界"


def test_hu101_official_cbl_window_and_formation_targets() -> None:
    """hu101：可评价段 5699.8-7810 为 cbl 窗（正式测量段 5390-7810 中双层套管段不评价）；>=3 地层目标。"""
    well_spec, windows = _windows("hu101")
    cbl = [w for w in windows if w.window_type == "cbl"]
    assert any(
        abs(w.top_md_m - 5699.8) < 1e-6 and abs(w.bottom_md_m - 7810.0) < 1e-6 for w in cbl
    )
    formation = [w for w in windows if w.window_type == "formation_target"]
    assert len(formation) >= 3


def test_hu103_formation_targets_added() -> None:
    """hu103：补齐气层 7507-7557 等 4 个地层目标窗。"""
    _, windows = _windows("hu103")
    formation = [w for w in windows if w.window_type == "formation_target"]
    assert len(formation) >= 4
    assert any(
        abs(w.top_md_m - 7507.0) < 1e-6 and abs(w.bottom_md_m - 7557.0) < 1e-6 for w in formation
    )


def test_hu1_cbl_quality_segments() -> None:
    """hu1：CBL 定性 16 段中最差段 7474-7514（不合格 40m）已结构化。"""
    _, windows = _windows("hu1")
    quality = [w for w in windows if w.window_type == "cbl_quality"]
    assert any(
        abs(w.top_md_m - 7474.0) < 1e-6 and abs(w.bottom_md_m - 7514.0) < 1e-6 for w in quality
    )


def test_hu2_formation_targets() -> None:
    """hu2：油气显示层 7402-2500 + 两个目的层 + 高压水层，均为 formation_target。"""
    _, windows = _windows("hu2")
    formation = [w for w in windows if w.window_type == "formation_target"]
    assert len(formation) >= 4
    assert any(abs(w.top_md_m - 7402.0) < 1e-6 for w in formation)


def test_ht1_001_no_cbl_window_but_formation() -> None:
    """ht1_001：CBL 定量缺失（无 cbl 窗），地层目标/录井显示窗 >=5。"""
    _, windows = _windows("ht1_001")
    assert not any(w.window_type == "cbl" for w in windows)
    formation = [w for w in windows if w.window_type == "formation_target"]
    assert len(formation) >= 5
    assert any(abs(w.top_md_m - 7460.0) < 1e-6 for w in formation)


def test_ht1_003_digitization_window() -> None:
    """ht1_003：数字化窗 5307.54-7514.21 与 cbl 窗分口径并存。"""
    _, windows = _windows("ht1_003")
    digit = [w for w in windows if w.window_type == "cbl_digitization"]
    assert any(
        abs(w.top_md_m - 5307.54) < 1e-6 and abs(w.bottom_md_m - 7514.21) < 1e-6 for w in digit
    )
    assert any(w.window_type == "cbl" for w in windows)


def test_ht1_004_cbl_window() -> None:
    """ht1_004：CBL 窗 5245-7581（数字化 0.3% 对应）。"""
    _, windows = _windows("ht1_004")
    cbl = [w for w in windows if w.window_type == "cbl"]
    assert any(
        abs(w.top_md_m - 5245.0) < 1e-6 and abs(w.bottom_md_m - 7581.0) < 1e-6 for w in cbl
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
