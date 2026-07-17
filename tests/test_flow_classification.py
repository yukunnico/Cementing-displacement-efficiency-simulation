"""流动分类判据单元测试（Zhang & Frigaard 2023 式 3.1-3.6）。

用合成 c̄(s,t) 场构造 mock AnnulusSimulationResult（SimpleNamespace），
不跑真实求解器，验证三类流动 regime 的判据分类与 NaN 守卫：

(a) 活塞流：w_f ≡ 1 → steady + non-dispersive；
(b) 随时间线性展宽的前沿（自相似 L∝t）：w_f(c̄) 线性展开 → unsteady + dispersive；
(c) 陡峭主前沿 + 低浓度前缘尖峰 + 高浓度残留尾 → steady + dispersive；
(d) 快照数 < 3 → NaN 守卫（flow_class="insufficient_data"）。
"""
from types import SimpleNamespace

import numpy as np
import pytest

from cemdisp.diagnostics.flow_classification import (
    FlowClassificationResult,
    compute_flow_classification,
)


def _make_mock_result(c_bar_st: np.ndarray, s: np.ndarray, times, ny: int = 8):
    """由截面平均浓度场 c_bar_st (nz, nt) 反构造 mock result。

    快照在 y 方向取常数（snap[y, j] = c_bar_st[j, i]），使 b 加权截面平均
    精确还原 c_bar_st；b 用非常数场以走真实加权路径。
    """
    s = np.asarray(s, dtype=float)
    times = tuple(float(t) for t in times)
    nz = s.size
    nt = len(times)
    assert c_bar_st.shape == (nz, nt)
    y = np.linspace(0.0, 1.0, ny)
    b = 0.02 * (1.0 + 0.3 * np.sin(np.pi * y))[:, None] * np.ones((1, nz))  # (ny, nz)
    snapshots = tuple(np.repeat(c_bar_st[:, i][None, :], ny, axis=0) for i in range(nt))
    geom = {"s": s, "b": b}
    return SimpleNamespace(
        cement_snapshots=snapshots,
        snapshot_times_s=times,
        geom=geom,
        metrics=None,
        summary={},
    )


class TestPistonFlow:
    """合成场 (a)：活塞流 c̄(s,t) = H(vt − s)，w_f ≡ 1。"""

    def _piston_result(self):
        s = np.linspace(0.0, 100.0, 400)
        times = [10.0, 20.0, 30.0, 40.0, 50.0]
        v = 1.0
        c_bar_st = (s[:, None] <= v * np.asarray(times)[None, :]).astype(float)
        return _make_mock_result(c_bar_st, s, times)

    def test_steady_non_dispersive(self):
        res = compute_flow_classification(self._piston_result())
        assert res.is_steady is True
        assert res.is_dispersive is False
        assert res.flow_class == "non_dispersive_steady"

    def test_criteria_magnitudes(self):
        res = compute_flow_classification(self._piston_result())
        # 阶梯剖面所有浓度等值线重合 → Δw_f ≈ 0，w_r ≈ 0
        assert res.delta_w_f <= 0.1
        assert res.sigma_wr_plus <= 0.08
        assert res.abs_wr_plus <= 0.05
        # ŵ₀ 由 c̄=0.5 等值线估计，应接近真实前缘速度 v = 1 m/s
        assert res.w0_m_s == pytest.approx(1.0, abs=0.01)

    def test_w0_override(self):
        # 显式传入 ŵ₀=2.0：w_f 减半 → w_r ≈ −0.5，仍 steady non-dispersive
        res = compute_flow_classification(self._piston_result(), w0_m_s=2.0)
        assert res.w0_m_s == pytest.approx(2.0)
        assert res.is_steady is True
        assert res.is_dispersive is False
        assert res.flow_class == "non_dispersive_steady"
        assert res.abs_wr_minus > 0.4  # 残留侧面积 ≈ 0.5


class TestDispersiveFlow:
    """合成场 (b)：线性展宽前沿 c̄ = clip((vt + L − s)/(2L), 0, 1)，L = αt。

    自相似：z(c̄,t)/t = v + α(1−2c̄) → w_f(c̄) = 1 + (α/v)(1−2c̄)。
    取 α/v = 0.3：Δw_f = 0.24 > 0.1（unsteady），
    σ_{w_r+} ≈ 0.12 > 0.08，|w̄_{r+}| ≈ 0.075 > 0.05（dispersive）。
    """

    def _dispersive_result(self):
        s = np.linspace(0.0, 100.0, 400)
        times = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
        v, alpha = 1.0, 0.3
        half_width = alpha * times
        c_bar_st = np.clip(
            (v * times[None, :] + half_width[None, :] - s[:, None])
            / (2.0 * half_width[None, :]),
            0.0,
            1.0,
        )
        return _make_mock_result(c_bar_st, s, times)

    def test_unsteady_dispersive(self):
        res = compute_flow_classification(self._dispersive_result())
        assert res.is_dispersive is True
        assert res.is_steady is False
        assert res.flow_class == "unsteady_dispersive"

    def test_criteria_magnitudes(self):
        res = compute_flow_classification(self._dispersive_result())
        # 线性剖面的等值线位置精确（线性插值无弦误差）
        assert res.delta_w_f == pytest.approx(0.24, abs=0.02)
        assert res.sigma_wr_plus == pytest.approx(0.12, abs=0.02)
        assert res.abs_wr_plus == pytest.approx(0.075, abs=0.015)
        assert res.abs_wr_minus == pytest.approx(0.075, abs=0.015)
        assert res.w0_m_s == pytest.approx(1.0, abs=1.0e-6)


class TestDispersiveSteadyFlow:
    """合成场 (c)：自相似剖面 c̄ = F((s−vt)/t)，主前沿陡峭 + 前缘尖峰 + 残留尾。

    F 折点 (u, F)：(-0.40, 1.0) → (-0.04, 0.7) → (0.04, 0.3) → (0.50, 0.0)。
    w_f(c̄) = 1 + u(c̄)：Δw_f = 0.08 ≤ 0.1（steady），
    σ_{w_r+} ≈ 0.16 > 0.08，|w̄_{r+}| ≈ 0.09 > 0.05（dispersive）。
    """

    def _result(self):
        s = np.linspace(0.0, 100.0, 1600)
        times = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
        v = 1.0
        u_pts = np.array([-0.40, -0.04, 0.04, 0.50])
        f_pts = np.array([1.00, 0.70, 0.30, 0.00])
        u = (s[:, None] - v * times[None, :]) / times[None, :]
        c_bar_st = np.interp(u.ravel(), u_pts, f_pts).reshape(u.shape)
        return _make_mock_result(c_bar_st, s, times)

    def test_dispersive_steady(self):
        res = compute_flow_classification(self._result())
        assert res.is_steady is True
        assert res.is_dispersive is True
        assert res.flow_class == "dispersive_steady"

    def test_criteria_magnitudes(self):
        res = compute_flow_classification(self._result())
        assert res.delta_w_f == pytest.approx(0.08, abs=0.02)
        assert res.sigma_wr_plus > 0.08
        assert res.abs_wr_plus > 0.05
        # 残留尾（w_r− 部）面积应显著：|w̄_{r−}| ≈ ∫₀.₅¹ |u| dc̄ ≈ 0.09
        assert res.abs_wr_minus > 0.05


class TestInsufficientSnapshots:
    """NaN 守卫：快照数 < 3 时全部指标为 NaN，flow_class="insufficient_data"。"""

    def test_two_snapshots_return_nan_guard(self):
        s = np.linspace(0.0, 100.0, 100)
        times = [10.0, 20.0]
        c_bar_st = (s[:, None] <= np.asarray(times)[None, :]).astype(float)
        res = compute_flow_classification(_make_mock_result(c_bar_st, s, times))
        assert res.flow_class == "insufficient_data"
        assert res.is_steady is False
        assert res.is_dispersive is False
        for value in (
            res.delta_w_f,
            res.sigma_wr_plus,
            res.abs_wr_plus,
            res.abs_wr_minus,
            res.w0_m_s,
        ):
            assert np.isnan(value)

    def test_zero_snapshots_return_nan_guard(self):
        s = np.linspace(0.0, 100.0, 100)
        c_bar_st = np.zeros((s.size, 0))
        res = compute_flow_classification(_make_mock_result(c_bar_st, s, []))
        assert res.flow_class == "insufficient_data"
        assert np.isnan(res.delta_w_f)

    def test_t0_snapshot_filtered_leaving_too_few(self):
        # 3 个快照但含 t=0 初始快照时仍可工作（过滤后剩 2 个有效时刻）；
        # 仅 t=0 + 1 个有效时刻 → 不足 2 个 → NaN 守卫
        s = np.linspace(0.0, 100.0, 100)
        times = [0.0, 10.0, 20.0]
        c_bar_st = (s[:, None] <= 0.5 * np.asarray(times)[None, :]).astype(float)
        res = compute_flow_classification(_make_mock_result(c_bar_st, s, times))
        assert res.flow_class != "insufficient_data"  # 2 个有效时刻足够
        times_short = [0.0, 10.0]
        c_bar_short = (s[:, None] <= 0.5 * np.asarray(times_short)[None, :]).astype(float)
        res_short = compute_flow_classification(_make_mock_result(c_bar_short, s, times_short))
        assert res_short.flow_class == "insufficient_data"

    def test_mismatched_lengths_raise(self):
        s = np.linspace(0.0, 100.0, 100)
        times = [10.0, 20.0, 30.0]
        c_bar_st = (s[:, None] <= np.asarray(times)[None, :]).astype(float)
        mock = _make_mock_result(c_bar_st, s, times)
        mock.snapshot_times_s = (10.0, 20.0)  # 人为制造不一致
        with pytest.raises(ValueError):
            compute_flow_classification(mock)


class TestResultInterface:
    """结果接口：frozen dataclass 与 to_dict()。"""

    def test_to_dict_keys_and_frozen(self):
        s = np.linspace(0.0, 100.0, 400)
        times = [10.0, 20.0, 30.0, 40.0, 50.0]
        c_bar_st = (s[:, None] <= np.asarray(times)[None, :]).astype(float)
        res = compute_flow_classification(_make_mock_result(c_bar_st, s, times))
        assert isinstance(res, FlowClassificationResult)
        d = res.to_dict()
        assert set(d.keys()) == {
            "delta_w_f",
            "sigma_wr_plus",
            "abs_wr_plus",
            "abs_wr_minus",
            "is_steady",
            "is_dispersive",
            "flow_class",
            "w0_m_s",
        }
        assert d["flow_class"] == res.flow_class
        with pytest.raises(AttributeError):  # frozen dataclass 不可变
            res.flow_class = "tampered"
