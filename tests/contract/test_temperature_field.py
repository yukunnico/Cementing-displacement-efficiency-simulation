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
    GEO_GRAD_C_PER_M,
    GEO_T0_C,
    ConstantTemperatureField,
    GeothermalTemperatureField,
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
# GeothermalTemperatureField（T2-1a：地温静温剖面，统一式单一真源）
# ---------------------------------------------------------------------------
class TestGeothermalTemperatureField:
    def test_module_constants_are_canonical(self):
        """模块常量即统一式 T(z)=16.006+1.7598e-2·z（单一真源）。"""
        assert GEO_T0_C == 16.006
        assert GEO_GRAD_C_PER_M == 1.7598e-2

    def test_linear_profile_formula(self):
        f = GeothermalTemperatureField()
        for md in (0.0, 30.0, 1000.0, 5241.0, 7660.0):
            assert f.T(md, 0.0) == pytest.approx(16.006 + 1.7598e-2 * md)

    def test_custom_params(self):
        f = GeothermalTemperatureField(T0_c=20.0, grad_c_per_m=0.03)
        assert f.T(100.0, 0.0) == pytest.approx(23.0)

    def test_time_invariance(self):
        """静温与时刻无关：任意 t 查询同值（时间维不消费）。"""
        f = GeothermalTemperatureField()
        for md in (0.0, 3000.0, 7660.0):
            assert f.T(md, 0.0) == f.T(md, 60.0) == f.T(md, 1e9)

    def test_same_interface_as_constant(self):
        """与 ConstantTemperatureField 完全同型：oob 恒空、reset_audit 空操作。"""
        f = GeothermalTemperatureField()
        assert f.oob_count == 0
        assert f.oob_events == ()
        assert isinstance(f.T(100, 60), float)  # 整型入参仍出 float
        f.T(-1.0, -1.0)
        f.T(1e12, 1e12)
        assert f.oob_count == 0 and f.oob_events == ()
        f.reset_audit()  # 空操作不报错
        assert f.oob_count == 0

    @staticmethod
    def _load_table_or_skip(tmp_path):
        """交付表可加载则返回 t_in，不可加载则 skip 并注明（呼应 brief 要求）。"""
        try:
            t_in, _t_out = load_delivered_pair(cache_dir=tmp_path)
        except Exception as exc:  # noqa: BLE001 —— 交付件缺失时跳过而非失败
            pytest.skip(f"交付温度表不可加载，跳过 col0 对照：{type(exc).__name__}: {exc}")
        return t_in

    def test_matches_table_col0_within_0_05c(self, tmp_path):
        """静温线 vs 交付表首列（col0=初始时刻）在抽查深度差 ≤0.05 °C。

        呼应拟合残差 ~0.04 °C（温压耦合改进计划 §1）；表不可加载则 skip。
        """
        t_in = self._load_table_or_skip(tmp_path)
        f = GeothermalTemperatureField()
        z = t_in.depth_m
        for row in (0, 50, 100, 166, 250, 332):
            got = f.T(float(z[row]), 0.0)
            assert abs(got - float(t_in.table[row, 0])) <= 0.05, (
                f"行 {row}（z={z[row]:.1f} m）：静温 {got:.4f} vs 表 col0 "
                f"{t_in.table[row, 0]:.4f} 超出 0.05 °C"
            )


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
# T_column（F-4 批量列聚合入口；Phase 4d 前置）
# ---------------------------------------------------------------------------
class TestColumnQuery:
    def test_column_matches_scalar_bitwise(self):
        """逐列批量值与逐点标量调用**逐位**同值（向量化不改运算次序）。"""
        f = _synthetic_field()
        mds = np.array([0.0, 2.5, 5.0, 9.9, 10.0, 15.0, 20.0])
        for t in (0.0, 30.0, 61.5, 180.0):
            got = f.T_column(mds, t)
            want = np.array([f.T(float(m), t) for m in mds])
            assert np.array_equal(got, want), f"t={t} 逐位失配"
        assert f.oob_count == 0  # 域内查询不记事件

    def test_column_shape_preserved(self):
        f = _synthetic_field()
        assert f.T_column([1.0, 2.0, 3.0], 60.0).shape == (3,)
        assert f.T_column(np.zeros((2, 3)), 60.0).shape == (2, 3)
        assert f.T_column(np.array([5.0]), 0.0).shape == (1,)

    def test_column_returns_float_even_with_int_inputs(self):
        f = _synthetic_field()
        got = f.T_column(np.array([0, 10, 20]), 0)
        assert got.dtype == np.float64

    def test_column_one_event_per_batch(self):
        """一批内 100 个域外点 ⇒ 事件列表只 +1，域外点数进聚合计数。"""
        f = _synthetic_field()
        mds = np.concatenate([np.full(50, -100.0), np.full(50, 1e6)])
        f.T_column(mds, 60.0)
        assert f.oob_count == 1
        assert f.oob_column_clamped_total == 100
        ev = f.oob_events[0]
        # 代表事件 = 偏移边界最远的域外点（此处 1e6 侧更远）
        assert ev.md_m == pytest.approx(1e6)
        assert ev.md_clamped_m == pytest.approx(f.depth_m[-1])

    def test_column_repeated_batches_grow_linearly(self):
        """逐批聚合：n 批 ⇒ 事件 n 条（而非 n×nz）——防逐列×逐步 append 爆表。"""
        f = _synthetic_field()
        mds = np.full(30, -5.0)
        for _ in range(7):
            f.T_column(mds, 60.0)
        assert f.oob_count == 7
        assert f.oob_column_clamped_total == 210

    def test_column_time_oob_single_event_no_total(self):
        """时间越界：每批 1 条代表事件，但**不计入**域外点数（同 Anchored 口径——
        该计数只统计深度域外）。"""
        f = _synthetic_field()
        f.T_column(np.array([0.0, 10.0, 20.0]), 1e9)
        assert f.oob_count == 1
        assert f.oob_column_clamped_total == 0
        assert f.oob_events[0].t_clamped_s == pytest.approx(f.time_s[-1])

    def test_column_in_domain_batch_no_audit(self):
        f = _synthetic_field()
        f.T_column(np.array([0.0, 5.0, 20.0]), 120.0)
        assert f.oob_count == 0
        assert f.oob_events == ()
        assert f.oob_column_clamped_total == 0

    def test_reset_audit_clears_column_total(self):
        f = _synthetic_field()
        f.T_column(np.array([-1.0, 30.0]), 0.0)
        assert f.oob_column_clamped_total == 2
        f.reset_audit()
        assert f.oob_count == 0
        assert f.oob_column_clamped_total == 0

    def test_scalar_path_unchanged(self):
        """批量入口不改标量路径的既有语义（关 1/关 2 红线：逐位保持）。"""
        f = _synthetic_field()
        f.T_column(np.array([-1000.0]), 0.0)   # 批量越界
        assert f.oob_count == 1
        assert f.T(-1000.0, 0.0) == pytest.approx(f.table[0, 0], abs=1e-12)
        assert f.oob_count == 2                 # 标量仍逐查询计数
        ev = f.oob_events[1]
        assert ev.md_m == -1000.0 and ev.md_clamped_m == pytest.approx(f.depth_m[0])


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
