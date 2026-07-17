"""displacement_metrics 单元测试（T0-4/5/7）。

用 SimpleNamespace mock AnnulusSimulationResult（不跑真实求解器），覆盖：
- 已知极限场：全水泥 / 全泥浆 / 全隔离液
- φm/φc 分级边界（Yang 2021 §4.4：0.05 归 green、0.1 归 red）
- 界面长度比（front_wide == front_narrow → 0）
- 窄四分位 η_N（含 ny 不整除 4 的单行边界）
- 突破时间 t_br（首次到达 / 未突破 → inf）与无量纲化 t_br_hat
- ŵ₀ 由 bulk_cement_fill 上升段斜率自动估计（含隔离液平台段）
"""
from dataclasses import FrozenInstanceError
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from cemdisp.diagnostics.displacement_metrics import (
    DisplacementMetricsResult,
    _classify_quality_zones,
    compute_displacement_metrics,
)

L_M = 100.0  # 井段长度 [m]


def _make_geom(ny: int = 8, nz: int = 5, L: float = L_M, b_value: float = 0.04) -> dict:
    """均匀间隙的合成几何（b 与 y、s 无关，便于手算核对积分比值）。"""
    s = np.linspace(0.0, L, nz)
    y = np.linspace(0.0, np.pi * 0.1, ny)
    b = np.full((ny, nz), b_value)
    return {"s": s, "y": y, "b": b, "md": L - s}


def _make_metrics(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows)


def _make_result(cement: np.ndarray, spacer: np.ndarray, geom: dict, metrics: pd.DataFrame) -> SimpleNamespace:
    return SimpleNamespace(
        well_name="mock_well",
        geom=geom,
        cement_field=cement,
        spacer_field=spacer,
        metrics=metrics,
    )


def _metrics_row(t, fw, fn, eff, fill):
    return {
        "time_s": float(t),
        "front_wide_m": float(fw),
        "front_narrow_m": float(fn),
        "effective_efficiency": float(eff),
        "bulk_cement_fill": float(fill),
    }


class TestKnownLimitFields:
    """已知场极限：全水泥 / 全泥浆 / 全隔离液。"""

    def test_all_cement_field(self):
        geom = _make_geom()
        cement = np.ones((8, 5))
        spacer = np.zeros((8, 5))
        metrics = _make_metrics([
            _metrics_row(0.0, 0.0, 0.0, 0.0, 0.0),
            _metrics_row(50.0, 100.0, 100.0, 1.0, 1.0),
        ])
        res = compute_displacement_metrics(_make_result(cement, spacer, geom, metrics))

        assert res.mud_retention_fraction == pytest.approx(0.0, abs=1e-9)
        assert res.eta_narrow == pytest.approx(1.0, abs=1e-9)
        assert res.eta_global == pytest.approx(1.0, abs=1e-9)
        # φm/φc = 0/max(1,ε) = 0 → 全部 green
        assert all(v == pytest.approx(0.0, abs=1e-9) for v in res.phi_m_phi_c_profile)
        assert res.quality_zone_counts == {"green": 5, "yellow": 0, "red": 0}
        # 宽窄边前缘齐头 → 界面长度比 0
        assert res.interface_length_ratio == pytest.approx(0.0, abs=1e-9)
        # t_br = 50 s；ŵ₀ 自动估计 = L·d(fill)/dt = 100·(1/50) = 2 m/s → t̂ = 50·2/100 = 1.0
        assert res.t_br_s == pytest.approx(50.0)
        assert res.w0_m_s == pytest.approx(2.0)
        assert res.t_br_hat == pytest.approx(1.0)

    def test_all_mud_field(self):
        geom = _make_geom()
        cement = np.zeros((8, 5))
        spacer = np.zeros((8, 5))
        metrics = _make_metrics([
            _metrics_row(0.0, 0.0, 0.0, 0.0, 0.0),
            _metrics_row(50.0, 0.0, 0.0, 0.0, 0.0),
        ])
        res = compute_displacement_metrics(_make_result(cement, spacer, geom, metrics))

        assert res.mud_retention_fraction == pytest.approx(1.0, abs=1e-9)
        assert res.eta_narrow == pytest.approx(0.0, abs=1e-9)
        # cement = 0 → φm/φc = 1/ε 巨大 → 全部 red
        assert res.quality_zone_counts == {"green": 0, "yellow": 0, "red": 5}
        # 前缘从未到达出口 → 未突破
        assert np.isinf(res.t_br_s)
        assert np.isinf(res.t_br_hat)
        assert any("未突破" in note for note in res.notes)

    def test_all_spacer_field_numerical_safety(self):
        """全隔离液：mud = 0 且 cement = 0，φm/φc = 0/ε = 0 → green，不产生 nan/inf。"""
        geom = _make_geom()
        cement = np.zeros((8, 5))
        spacer = np.ones((8, 5))
        metrics = _make_metrics([_metrics_row(0.0, 0.0, 0.0, 0.0, 0.0)])
        res = compute_displacement_metrics(_make_result(cement, spacer, geom, metrics))

        assert res.mud_retention_fraction == pytest.approx(0.0, abs=1e-9)
        assert np.all(np.isfinite(res.phi_m_phi_c_profile))
        assert res.quality_zone_counts == {"green": 5, "yellow": 0, "red": 0}
        assert res.eta_narrow == pytest.approx(0.0, abs=1e-9)


class TestPhiMPhiCClassification:
    """φm/φc 分级边界（Yang 2021 §4.4：≤0.05 green，(0.05,0.1) yellow，≥0.1 red）。"""

    def test_zone_boundaries(self):
        counts = _classify_quality_zones((0.0, 0.04, 0.05, 0.07, 0.099, 0.1, 0.2))
        # 0.0/0.04/0.05 → green；0.07/0.099 → yellow；0.1/0.2 → red
        assert counts == {"green": 3, "yellow": 2, "red": 2}

    def test_profile_from_fields(self):
        """逐列构造已知比值：均匀 b 下 φm/φc = mud_col/cement_col。

        注意：精确边界值（0.05/0.1）的分级由 test_zone_boundaries 用精确浮点覆盖；
        本测试经 1−cement−spacer 反算 mud，存在 ~1e-17 浮点噪声，故用区间内部值。
        """
        ny, nz = 4, 6
        geom = _make_geom(ny=ny, nz=nz)
        ratios = [0.04, 0.03, 0.07, 0.12, 0.2, 0.0]
        cement = np.zeros((ny, nz))
        spacer = np.zeros((ny, nz))
        for j, r in enumerate(ratios):
            cement[:, j] = 0.5
            mud_col = 0.5 * r
            spacer[:, j] = 1.0 - 0.5 - mud_col  # mud = 1 − cement − spacer = 0.5·r
        metrics = _make_metrics([_metrics_row(0.0, 0.0, 0.0, 0.5, 0.5)])
        res = compute_displacement_metrics(_make_result(cement, spacer, geom, metrics))

        assert res.phi_m_phi_c_profile == pytest.approx(tuple(ratios), rel=1e-9)
        # 0.04/0.03/0.0 → green；0.07 → yellow；0.12/0.2 → red
        assert res.quality_zone_counts == {"green": 3, "yellow": 1, "red": 2}


class TestInterfaceLengthRatio:
    def test_equal_fronts_gives_zero(self):
        geom = _make_geom()
        metrics = _make_metrics([_metrics_row(10.0, 80.0, 80.0, 0.8, 0.8)])
        res = compute_displacement_metrics(
            _make_result(np.ones((8, 5)), np.zeros((8, 5)), geom, metrics)
        )
        assert res.interface_length_ratio == pytest.approx(0.0, abs=1e-9)

    def test_known_ratio(self):
        """|90 − 60| / mean(90, 60) = 30/75 = 0.4。"""
        geom = _make_geom()
        metrics = _make_metrics([_metrics_row(10.0, 90.0, 60.0, 0.8, 0.8)])
        res = compute_displacement_metrics(
            _make_result(np.ones((8, 5)), np.zeros((8, 5)), geom, metrics)
        )
        assert res.interface_length_ratio == pytest.approx(0.4, rel=1e-9)

    def test_both_fronts_zero_no_nan(self):
        """水泥尚未入环空（两前缘均 0）时由 ε 保护返回 0，不产生 nan。"""
        geom = _make_geom()
        metrics = _make_metrics([_metrics_row(0.0, 0.0, 0.0, 0.0, 0.0)])
        res = compute_displacement_metrics(
            _make_result(np.zeros((8, 5)), np.zeros((8, 5)), geom, metrics)
        )
        assert res.interface_length_ratio == pytest.approx(0.0, abs=1e-9)
        assert np.isfinite(res.interface_length_ratio)


class TestNarrowQuarterEfficiency:
    def test_narrow_quarter_empty(self):
        """ny=8 → 窄四分位为最后 2 行；窄边水泥为 0 → η_N = 0，全局 η_E = 6/8 = 0.75。"""
        geom = _make_geom(ny=8)
        cement = np.ones((8, 5))
        cement[6:, :] = 0.0
        metrics = _make_metrics([_metrics_row(10.0, 100.0, 40.0, 0.75, 0.75)])
        res = compute_displacement_metrics(_make_result(cement, np.zeros((8, 5)), geom, metrics))

        assert res.eta_narrow == pytest.approx(0.0, abs=1e-9)
        assert res.eta_global == pytest.approx(0.75, abs=1e-9)
        assert res.eta_narrow <= res.eta_global  # 物理合理：窄边不优于全局

    def test_narrow_quarter_full(self):
        geom = _make_geom(ny=8)
        cement = np.ones((8, 5))
        metrics = _make_metrics([_metrics_row(10.0, 100.0, 100.0, 1.0, 1.0)])
        res = compute_displacement_metrics(_make_result(cement, np.zeros((8, 5)), geom, metrics))
        assert res.eta_narrow == pytest.approx(1.0, abs=1e-9)

    def test_single_row_quarter_boundary(self):
        """ny=5 → n_q = max(1, 5//4) = 1 行；单行梯形积分退化保护，不应得 nan/0 误判。"""
        geom = _make_geom(ny=5)
        cement = np.ones((5, 5))
        cement[4, :] = 0.5  # 最后一行（窄边）水泥 0.5
        metrics = _make_metrics([_metrics_row(10.0, 100.0, 50.0, 0.9, 0.9)])
        res = compute_displacement_metrics(_make_result(cement, np.zeros((5, 5)), geom, metrics))

        assert np.isfinite(res.eta_narrow)
        assert res.eta_narrow == pytest.approx(0.5, abs=1e-9)


class TestBreakthroughTime:
    def test_first_breakthrough_from_front_series(self):
        """front_wide 在 t=20 首次到达 s_max=100；front_narrow 未到达 → t_br = 20 s。"""
        geom = _make_geom()
        metrics = _make_metrics([
            _metrics_row(0.0, 0.0, 0.0, 0.0, 0.0),
            _metrics_row(10.0, 50.0, 30.0, 0.3, 0.3),
            _metrics_row(20.0, 100.0, 60.0, 0.5, 0.5),
            _metrics_row(30.0, 100.0, 95.0, 0.6, 0.6),
        ])
        res = compute_displacement_metrics(
            _make_result(np.ones((8, 5)), np.zeros((8, 5)), geom, metrics),
            w0_m_s=0.5,
        )
        assert res.t_br_s == pytest.approx(20.0)
        # t̂_br = t_br·ŵ₀/L = 20·0.5/100 = 0.1
        assert res.t_br_hat == pytest.approx(0.1)
        assert res.w0_m_s == pytest.approx(0.5)

    def test_narrow_side_earlier_breakthrough(self):
        """取宽窄边较早者：front_narrow 在 t=15 先到达 → t_br = 15 s。"""
        geom = _make_geom()
        metrics = _make_metrics([
            _metrics_row(0.0, 0.0, 0.0, 0.0, 0.0),
            _metrics_row(15.0, 80.0, 100.0, 0.4, 0.4),
            _metrics_row(30.0, 100.0, 100.0, 0.6, 0.6),
        ])
        res = compute_displacement_metrics(
            _make_result(np.ones((8, 5)), np.zeros((8, 5)), geom, metrics),
            w0_m_s=1.0,
        )
        assert res.t_br_s == pytest.approx(15.0)
        assert res.t_br_hat == pytest.approx(15.0 * 1.0 / 100.0)

    def test_no_breakthrough_gives_inf(self):
        geom = _make_geom()
        metrics = _make_metrics([
            _metrics_row(0.0, 0.0, 0.0, 0.0, 0.0),
            _metrics_row(30.0, 90.0, 60.0, 0.5, 0.5),
        ])
        res = compute_displacement_metrics(
            _make_result(np.ones((8, 5)), np.zeros((8, 5)), geom, metrics),
            w0_m_s=1.0,
        )
        assert np.isinf(res.t_br_s)
        assert np.isinf(res.t_br_hat)
        assert any("未突破" in note for note in res.notes)


class TestW0Estimation:
    def test_w0_from_bulk_fill_slope_with_spacer_plateau(self):
        """fill = [0, 0, 0.2, 0.4]（前 10 s 注隔离液平台段）→ 正斜率中位数 0.02 /s
        → ŵ₀ = L·slope = 100·0.02 = 2.0 m/s；t_br = 30 s → t̂ = 30·2/100 = 0.6。"""
        geom = _make_geom()
        metrics = _make_metrics([
            _metrics_row(0.0, 0.0, 0.0, 0.0, 0.0),
            _metrics_row(10.0, 0.0, 0.0, 0.0, 0.0),
            _metrics_row(20.0, 50.0, 40.0, 0.2, 0.2),
            _metrics_row(30.0, 100.0, 90.0, 0.4, 0.4),
        ])
        res = compute_displacement_metrics(
            _make_result(np.ones((8, 5)), np.zeros((8, 5)), geom, metrics)
        )
        assert res.t_br_s == pytest.approx(30.0)
        assert res.w0_m_s == pytest.approx(2.0)
        assert res.t_br_hat == pytest.approx(0.6)

    def test_w0_not_estimable_gives_nan_t_br_hat(self):
        """fill 恒为 0（无上升段）但前缘到达出口（假设性数据）→ t_br_hat = nan + 警告备注。"""
        geom = _make_geom()
        metrics = _make_metrics([
            _metrics_row(0.0, 0.0, 0.0, 0.0, 0.0),
            _metrics_row(10.0, 100.0, 100.0, 0.0, 0.0),
        ])
        res = compute_displacement_metrics(
            _make_result(np.ones((8, 5)), np.zeros((8, 5)), geom, metrics)
        )
        assert res.t_br_s == pytest.approx(10.0)
        assert np.isnan(res.t_br_hat)
        assert any("w0_m_s" in note for note in res.notes)


class TestResultStructure:
    def test_to_dict_keys_and_types(self):
        geom = _make_geom()
        metrics = _make_metrics([_metrics_row(10.0, 100.0, 100.0, 1.0, 1.0)])
        res = compute_displacement_metrics(
            _make_result(np.ones((8, 5)), np.zeros((8, 5)), geom, metrics)
        )
        d = res.to_dict()
        expected_keys = {
            "mud_retention_fraction", "phi_m_phi_c_profile", "quality_zone_counts",
            "interface_length_ratio", "eta_narrow", "eta_global",
            "t_br_s", "t_br_hat", "w0_m_s", "notes",
        }
        assert set(d.keys()) == expected_keys
        assert isinstance(d["phi_m_phi_c_profile"], list)
        assert isinstance(d["quality_zone_counts"], dict)
        assert isinstance(d["notes"], list)

    def test_result_is_frozen(self):
        geom = _make_geom()
        metrics = _make_metrics([_metrics_row(10.0, 100.0, 100.0, 1.0, 1.0)])
        res = compute_displacement_metrics(
            _make_result(np.ones((8, 5)), np.zeros((8, 5)), geom, metrics)
        )
        assert isinstance(res, DisplacementMetricsResult)
        with pytest.raises(FrozenInstanceError):
            res.mud_retention_fraction = 0.5  # type: ignore[misc]
