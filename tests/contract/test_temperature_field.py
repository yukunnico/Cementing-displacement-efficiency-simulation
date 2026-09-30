"""
T0-1 温度场表数据层契约测试（cemdisp.data.temperature_field）

覆盖 brief 需求：
- 抽查节点：T(z[行], t[列]) 精确命中表内值（含独立硬编码锚点，防解析器同错）
- 双线性：两节点中点 = 两节点均值
- 边界：z/t 越界 clamp 到端值 + 审计计数/位置可查询
- 构造校验：形状错 / 深度乱序 / NaN 抛 ValueError
- ConstantTemperatureField：任意 (z,t) 返回常数
- npz 缓存：二次加载免解析 xlsx
- 井底（最深行）两表共享同值
"""

from __future__ import annotations

import numpy as np
import pytest

from cemdisp.data.temperature_field import (
    ConstantTemperatureField,
    TableTemperatureField,
    load_delivered_pair,
)


# ---------------------------------------------------------------------------
# 合成小表（校验与插值单测用；生产默认形状为 333×200）
# ---------------------------------------------------------------------------
def _synthetic_field() -> TableTemperatureField:
    table = np.arange(12, dtype=float).reshape(3, 4)  # 3 深度 × 4 时间
    depth_m = np.array([0.0, 10.0, 20.0])
    time_s = np.arange(4, dtype=float) * 60.0
    return TableTemperatureField(table, depth_m, time_s, expected_shape=(3, 4))


# ---------------------------------------------------------------------------
# ConstantTemperatureField
# ---------------------------------------------------------------------------
class TestConstantTemperatureField:
    def test_any_query_returns_constant(self):
        f = ConstantTemperatureField(T_c=42.5)
        assert f.T(0.0, 0.0) == pytest.approx(42.5)
        assert f.T(7660.0, 1e9) == pytest.approx(42.5)
        assert f.T(-1.0, -1.0) == pytest.approx(42.5)
        # 恒温场永不越界
        assert f.oob_count == 0

    def test_integer_input_still_float(self):
        f = ConstantTemperatureField(T_c=20)
        assert isinstance(f.T(100, 60), float)


# ---------------------------------------------------------------------------
# 构造校验
# ---------------------------------------------------------------------------
class TestConstructionValidation:
    def test_shape_mismatch_raises(self):
        # 默认期望 333×200，给 100×200 → 报形状错并说清期望/实际
        with pytest.raises(ValueError, match="333"):
            TableTemperatureField(np.zeros((100, 200)), np.arange(100.0))

    def test_depth_length_mismatch_raises(self):
        with pytest.raises(ValueError, match="333"):
            TableTemperatureField(np.zeros((333, 200)), np.arange(100.0))

    def test_depth_not_monotonic_raises(self):
        depth = np.array([0.0, 20.0, 10.0])  # 乱序
        with pytest.raises(ValueError, match="单调"):
            TableTemperatureField(np.zeros((3, 4)), depth, expected_shape=(3, 4))

    def test_nan_in_table_raises(self):
        table = np.zeros((3, 4))
        table[1, 2] = np.nan
        with pytest.raises(ValueError, match="NaN"):
            TableTemperatureField(
                table, np.array([0.0, 10.0, 20.0]), expected_shape=(3, 4)
            )

    def test_nan_in_depth_raises(self):
        depth = np.array([0.0, np.nan, 20.0])
        with pytest.raises(ValueError, match="NaN"):
            TableTemperatureField(np.zeros((3, 4)), depth, expected_shape=(3, 4))

    def test_time_length_mismatch_raises(self):
        with pytest.raises(ValueError, match="时间轴"):
            TableTemperatureField(
                np.zeros((3, 4)),
                np.array([0.0, 10.0, 20.0]),
                time_s=np.arange(3, dtype=float),
                expected_shape=(3, 4),
            )


# ---------------------------------------------------------------------------
# 插值：节点精确命中（交付真表）
# ---------------------------------------------------------------------------
@pytest.fixture(scope="session")
def pair(tmp_path_factory):
    """交付两表（管内, 环空），npz 缓存落在临时目录，不污染参考文档。"""
    cache_dir = tmp_path_factory.mktemp("tcache")
    return load_delivered_pair(cache_dir=cache_dir)


class TestNodeExactHit:
    # 独立硬编码锚点（直接从 xlsx 底层单元格取值，防"读表错+查询错"同错）
    ANCHORS = [
        ("T_in", 0, 0, 16.527999998701688),
        ("T_in", 100, 50, 48.33138788001854),
        ("T_in", 332, 199, 102.22497569040068),
        ("T_out", 100, 50, 75.57940602446222),
        ("T_out", 332, 199, 102.22497569040068),
    ]
    # ≥5 个 (行,列) 抽查点
    SAMPLE_POINTS = [(0, 0), (100, 50), (332, 199), (200, 0), (50, 199), (150, 75)]

    def test_axes_from_delivered_files(self, pair):
        t_in, _t_out = pair
        assert t_in.table.shape == (333, 200)
        assert len(t_in.depth_m) == 333
        assert t_in.depth_m[0] == pytest.approx(30.0)
        assert t_in.depth_m[-1] == pytest.approx(7660.0)
        # 时间轴：col0=初始时刻，步长 1 min
        assert t_in.time_s[0] == 0.0
        assert t_in.time_s[1] == 60.0
        assert t_in.time_s[-1] == pytest.approx(199 * 60.0)

    def test_sample_points_hit_table_values(self, pair):
        t_in, t_out = pair
        for field, table in ((t_in, t_in.table), (t_out, t_out.table)):
            for row, col in self.SAMPLE_POINTS:
                got = field.T(field.depth_m[row], field.time_s[col])
                assert got == pytest.approx(table[row, col], abs=1e-9), (
                    f"节点 ({row},{col}) 未精确命中: {got} vs {table[row, col]}"
                )

    @pytest.mark.parametrize("which,row,col,expected", ANCHORS)
    def test_hardcoded_anchors(self, pair, which, row, col, expected):
        t_in, t_out = pair
        field = t_in if which == "T_in" else t_out
        got = field.T(field.depth_m[row], field.time_s[col])
        assert got == pytest.approx(expected, abs=1e-9)


# ---------------------------------------------------------------------------
# 双线性：线性轴上中点 = 两节点均值
# ---------------------------------------------------------------------------
class TestBilinear:
    def test_depth_midpoint_is_mean(self):
        f = _synthetic_field()
        z_mid = 0.5 * (f.depth_m[0] + f.depth_m[1])
        expected = 0.5 * (f.table[0, 2] + f.table[1, 2])
        assert f.T(z_mid, f.time_s[2]) == pytest.approx(expected, abs=1e-12)

    def test_time_midpoint_is_mean(self):
        f = _synthetic_field()
        t_mid = 0.5 * (f.time_s[1] + f.time_s[2])
        expected = 0.5 * (f.table[1, 1] + f.table[1, 2])
        assert f.T(f.depth_m[1], t_mid) == pytest.approx(expected, abs=1e-12)

    def test_corner_blend_on_synthetic(self):
        f = _synthetic_field()
        # z=5（0~10 中点）、t=30s（0~60 中点）→ 四角均值
        got = f.T(5.0, 30.0)
        corners = [f.table[0, 0], f.table[0, 1], f.table[1, 0], f.table[1, 1]]
        assert got == pytest.approx(float(np.mean(corners)), abs=1e-12)


# ---------------------------------------------------------------------------
# 越界 clamp + 审计
# ---------------------------------------------------------------------------
class TestClampAndAudit:
    def test_depth_below_clamps_and_audits(self):
        f = _synthetic_field()
        got = f.T(-1000.0, f.time_s[2])
        assert got == pytest.approx(f.table[0, 2], abs=1e-12)
        assert f.oob_count == 1
        ev = f.oob_events[0]
        assert ev.md_m == -1000.0
        assert ev.md_clamped_m == pytest.approx(f.depth_m[0])

    def test_depth_above_clamps_and_audits(self):
        f = _synthetic_field()
        got = f.T(99999.0, f.time_s[1])
        assert got == pytest.approx(f.table[-1, 1], abs=1e-12)
        assert f.oob_count == 1
        assert f.oob_events[0].md_clamped_m == pytest.approx(f.depth_m[-1])

    def test_time_clamps_and_audits(self):
        f = _synthetic_field()
        assert f.T(f.depth_m[1], -5.0) == pytest.approx(f.table[1, 0], abs=1e-12)
        assert f.T(f.depth_m[1], 1e12) == pytest.approx(f.table[1, -1], abs=1e-12)
        assert f.oob_count == 2
        assert f.oob_events[1].t_clamped_s == pytest.approx(f.time_s[-1])

    def test_in_domain_query_no_audit(self):
        f = _synthetic_field()
        f.T(f.depth_m[1], f.time_s[1])
        f.T(5.0, 30.0)
        assert f.oob_count == 0
        assert f.oob_events == ()

    def test_reset_audit(self):
        f = _synthetic_field()
        f.T(-1.0, -1.0)
        assert f.oob_count == 1
        f.reset_audit()
        assert f.oob_count == 0
        assert f.oob_events == ()


# ---------------------------------------------------------------------------
# npz 缓存
# ---------------------------------------------------------------------------
class TestNpzCache:
    def test_cache_written_and_second_load_skips_xlsx(
        self, tmp_path, monkeypatch
    ):
        cache = tmp_path / "T_in.npz"
        f1 = TableTemperatureField.from_files(cache_path=cache)
        assert cache.exists()

        def _boom(*_a, **_k):
            raise AssertionError("二次加载不应再解析 xlsx")

        monkeypatch.setattr(
            "cemdisp.data.temperature_field._read_xlsx_numeric", _boom
        )
        f2 = TableTemperatureField.from_files(cache_path=cache)
        np.testing.assert_array_equal(f1.table, f2.table)
        np.testing.assert_array_equal(f1.depth_m, f2.depth_m)
        np.testing.assert_array_equal(f1.time_s, f2.time_s)

    def test_use_cache_false_forces_reparse(self, tmp_path):
        cache = tmp_path / "T_in.npz"
        TableTemperatureField.from_files(cache_path=cache, use_cache=False)
        assert not cache.exists()  # 关缓存则不落盘
        f = TableTemperatureField.from_files(cache_path=cache, use_cache=False)
        assert f.table.shape == (333, 200)


# ---------------------------------------------------------------------------
# 交付表一致性：井底两表共享同值
# ---------------------------------------------------------------------------
class TestDeliveredPair:
    def test_bottom_row_shared(self, tmp_path):
        t_in, t_out = load_delivered_pair(cache_dir=tmp_path)
        assert t_in.table.shape == t_out.table.shape == (333, 200)
        np.testing.assert_allclose(
            t_in.table[-1], t_out.table[-1], atol=1e-9, rtol=0
        )
        # 两表非同一份数据（管内≠环空，中段应有差异）
        assert not np.allclose(t_in.table[100], t_out.table[100], atol=1e-9)
