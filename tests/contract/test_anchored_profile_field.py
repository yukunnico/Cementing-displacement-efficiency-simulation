"""4d AnchoredProfileField 契约测试（静温锚点剖面场；Phase 4 设计规格 §1 项 5）。

覆盖：
- 分段线性插值正确性（手算锚：锚点命中 / 段内线性 / 域外地温回退）
- 温度系数端点行为：k=1 等于纯静温锚；k 必给、无默认（用户硬停点）
- oob 三件套同型 + F-4 列批量聚合计数（逐批至多 1 条事件、域外点数另计）
- 与既有三场协议同型（同一调用形态 T/oob_count/oob_events/reset_audit）
- 族混装构造期拒绝（同 md 异 T 抛 ValueError；同 md 同值合并）
- 瞬态接口位（t=0 静温、t 大收敛准稳态、关闭时时间维不消费）
- 既有场类零改动旁证（模块常量/形状/既有行为锚点复查；全量见
  test_temperature_field.py 同批重跑）
"""

from __future__ import annotations

import inspect

import numpy as np
import pytest

from cemdisp.data.temperature_field import (
    EXPECTED_SHAPE,
    GEO_GRAD_C_PER_M,
    GEO_T0_C,
    AnchoredProfileField,
    ClampEvent,
    ConstantTemperatureField,
    GeothermalTemperatureField,
    TableTemperatureField,
)

# 手算锚（°C）：段 1000→3000 斜率 0.02，段 3000→5000 斜率 0.015
ANCHORS = [(1000.0, 50.0), (3000.0, 90.0), (5000.0, 120.0)]


def _static(k: float = 1.0, **kw) -> AnchoredProfileField:
    return AnchoredProfileField(ANCHORS, k, regime="static", **kw)


def _circ(k: float, **kw) -> AnchoredProfileField:
    return AnchoredProfileField(ANCHORS, k, regime="circulating", **kw)


# ---------------------------------------------------------------------------
# 分段线性插值（手算锚）
# ---------------------------------------------------------------------------
class TestPiecewiseLinear:
    def test_anchor_nodes_exact_hit(self):
        f = _static()
        assert f.T(1000.0, 0.0) == 50.0
        assert f.T(3000.0, 0.0) == 90.0
        assert f.T(5000.0, 0.0) == 120.0

    def test_segment_midpoints_hand_computed(self):
        """段内线性手算：T(2000)=50+0.5*40=70；T(4500)=90+0.75*30=112.5。"""
        f = _static()
        assert f.T(2000.0, 0.0) == pytest.approx(70.0, abs=1e-12)
        assert f.T(4500.0, 0.0) == pytest.approx(112.5, abs=1e-12)

    def test_unsorted_input_is_sorted(self):
        f = AnchoredProfileField(list(reversed(ANCHORS)), 1.0, regime="static")
        assert f.T(2000.0, 0.0) == pytest.approx(70.0, abs=1e-12)

    def test_return_is_float_even_with_int_inputs(self):
        f = _static()
        got = f.T(2000, 0)
        assert isinstance(got, float)


# ---------------------------------------------------------------------------
# 锚域外回退地温式（参数化开关，默认声明）
# ---------------------------------------------------------------------------
class TestGeothermalFallback:
    def test_default_fallback_is_geothermal_line(self):
        """默认 fallback_geothermal=True：域外与 GeothermalTemperatureField 同值。"""
        f = _static()
        geo = GeothermalTemperatureField()
        assert f.fallback_geothermal is True
        for md in (0.0, 500.0, 6000.0, 7868.0):
            assert f.T(md, 0.0) == geo.T(md, 0.0)  # 同一算式，逐位
        assert f.T(0.0, 0.0) == pytest.approx(GEO_T0_C, abs=1e-12)
        assert f.T(6000.0, 0.0) == pytest.approx(
            GEO_T0_C + GEO_GRAD_C_PER_M * 6000.0, abs=1e-12
        )

    def test_endpoint_clamp_when_fallback_off(self):
        f = _static(fallback_geothermal=False)
        assert f.T(0.0, 0.0) == 50.0      # 钳首锚值
        assert f.T(7868.0, 0.0) == 120.0  # 钳末锚值

    def test_fallback_still_audits_oob(self):
        """域外回退/钳位都记越界审计（两开关下都记）。"""
        for fb in (True, False):
            f = _static(fallback_geothermal=fb)
            f.T(0.0, 0.0)
            f.T(7868.0, 0.0)
            f.T(2000.0, 0.0)  # 域内不记
            assert f.oob_count == 2, f"fallback_geothermal={fb}"
            ev = f.oob_events[0]
            assert isinstance(ev, ClampEvent)
            assert ev.md_m == 0.0 and ev.md_clamped_m == 1000.0
            assert f.oob_events[1].md_clamped_m == 5000.0


# ---------------------------------------------------------------------------
# 温度系数 k：端点行为 + 必给无默认（取值 = 用户硬停点）
# ---------------------------------------------------------------------------
class TestTemperatureFactor:
    def test_k_is_required_no_default(self):
        """构造器不得自带默认系数：缺参即 TypeError，签名无默认值。"""
        with pytest.raises(TypeError):
            AnchoredProfileField(ANCHORS)
        sig = inspect.signature(AnchoredProfileField.__init__)
        p = sig.parameters["temperature_factor_k"]
        assert p.default is inspect.Parameter.empty
        assert p.kind is p.POSITIONAL_OR_KEYWORD

    def test_k_one_equals_pure_static_anchor(self):
        """k=1 ⇒ 循环准稳态退化为纯静温锚（同值集合逐位相等）。"""
        circ = _circ(1.0)
        stat = _static(1.0)
        for md in (0.0, 1000.0, 2000.0, 3000.0, 4500.0, 5000.0, 7868.0):
            assert circ.T(md, 0.0) == stat.T(md, 0.0)

    def test_k_multiplies_static_line(self):
        """乘性口径：T = k · 静温线（k=0.85 ⇒ T(2000)=0.85*70=59.5）。"""
        f = _circ(0.85)
        assert f.T(2000.0, 0.0) == pytest.approx(59.5, abs=1e-12)
        assert f.T(1000.0, 0.0) == pytest.approx(42.5, abs=1e-12)
        # 域外同样先回退地温线再乘 k
        geo = GeothermalTemperatureField()
        assert f.T(7868.0, 0.0) == pytest.approx(0.85 * geo.T(7868.0, 0.0), abs=1e-12)

    @pytest.mark.parametrize("bad_k", [0.0, -0.5, 1.0001, 2.0, float("nan"), float("inf")])
    def test_k_outside_unit_interval_raises(self, bad_k):
        with pytest.raises(ValueError, match="温度系数"):
            AnchoredProfileField(ANCHORS, bad_k)


# ---------------------------------------------------------------------------
# 瞬态接口位（一阶集总 Ramey/Hasan-Kabir 型混合）
# ---------------------------------------------------------------------------
class TestTransientInterface:
    def test_tau_none_time_dimension_not_consumed(self):
        """缺省关瞬态：任意 t（含负值）同值、不审计——与 Geothermal 同型。"""
        f = _circ(0.85)
        assert f.T(2000.0, 0.0) == f.T(2000.0, 3600.0) == f.T(2000.0, 1e9)
        assert f.T(2000.0, -5.0) == f.T(2000.0, 0.0)
        assert f.oob_count == 0

    def test_tau_starts_at_static_and_relaxes_to_quasi_steady(self):
        stat = _static()
        f = _circ(0.85, transient_tau_s=100.0)
        assert f.T(2000.0, 0.0) == pytest.approx(stat.T(2000.0, 0.0), abs=1e-9)
        t_qs = f.T(2000.0, 0.0) * 0.85
        assert f.T(2000.0, 5000.0) == pytest.approx(t_qs, rel=1e-6)
        mid = f.T(2000.0, 100.0)
        assert t_qs < mid < stat.T(2000.0, 0.0)  # 单调衰减
        assert mid == pytest.approx(
            t_qs + (59.5 / 0.85 - t_qs) * np.exp(-1.0), abs=1e-9
        )

    def test_negative_time_clamps_and_audits_only_when_transient_active(self):
        f = _circ(0.85, transient_tau_s=100.0)
        assert f.T(2000.0, -30.0) == f.T(2000.0, 0.0)
        assert f.oob_count == 1
        assert f.oob_events[0].t_clamped_s == 0.0

    def test_bad_tau_raises(self):
        for bad in (0.0, -1.0, float("nan")):
            with pytest.raises(ValueError, match="transient_tau_s"):
                AnchoredProfileField(ANCHORS, 0.85, transient_tau_s=bad)


# ---------------------------------------------------------------------------
# oob 审计：标量逐查询（同 Table 型）+ F-4 列批量聚合计数
# ---------------------------------------------------------------------------
class TestOobAuditAndColumnAggregate:
    def test_scalar_path_appends_per_query(self):
        f = _circ(0.85)
        f.T(0.0, 0.0)
        f.T(-500.0, 0.0)
        f.T(7868.0, 0.0)
        assert f.oob_count == 3
        assert f.oob_column_clamped_total == 0

    def test_column_path_one_event_per_batch(self):
        """逐查询 append 爆表对策：一批 4 个域外点只追加 1 条代表事件。"""
        f = _circ(0.85)
        col = f.T_column([-2000.0, -500.0, 2000.0, 7868.0, 9000.0, 100000.0])
        assert col.shape == (6,)
        assert f.oob_count == 1                      # 每批至多 1 条
        assert f.oob_column_clamped_total == 5       # 域外点数另计
        ev = f.oob_events[0]
        assert isinstance(ev, ClampEvent)
        assert ev.md_m == 100000.0                   # 代表 = 偏移边界最远点
        assert ev.md_clamped_m == 5000.0

    def test_column_values_match_scalar_path(self):
        """列批量与逐点标量同值（含域外回退段与代表标量形状）。"""
        f = _circ(0.875)
        mds = np.array([-100.0, 0.0, 1000.0, 2500.0, 3000.0, 5000.0, 6500.0, 7868.0])
        col = f.T_column(mds, t_s=0.0)
        per_pt = np.array([f.T(m, 0.0) for m in mds])
        np.testing.assert_allclose(col, per_pt, rtol=0, atol=0)  # 要求逐位

    def test_column_in_domain_batch_no_audit(self):
        f = _circ(0.85)
        f.T_column([1000.0, 3000.0, 5000.0, 2000.0])
        assert f.oob_count == 0 and f.oob_column_clamped_total == 0

    def test_reset_audit_clears_both_channels(self):
        f = _circ(0.85)
        f.T(0.0, 0.0)
        f.T_column([0.0, 6000.0])
        f.reset_audit()
        assert f.oob_count == 0
        assert f.oob_events == ()
        assert f.oob_column_clamped_total == 0

    def test_column_preserves_input_shape(self):
        f = _circ(0.85)
        col = f.T_column([[1500.0, 2500.0], [3500.0, 4500.0]])
        assert col.shape == (2, 2)


# ---------------------------------------------------------------------------
# 族混装构造期拒绝（测绘项16 / 冲突清单 9）
# ---------------------------------------------------------------------------
class TestFamilyCollision:
    def test_same_md_different_T_raises(self):
        """同 md 两族（如呼1-004 5241 m 电测124 vs 实验155）必须拒绝。"""
        with pytest.raises(ValueError, match="锚点冲突"):
            AnchoredProfileField([(5241.0, 124.0), (5241.0, 155.0), (7660.0, 155.0)], 0.85)

    def test_same_md_same_T_merged(self):
        """跨文档同值重复行（如呼1-004 7660 m 155 ×2）自动合并不报错。"""
        f = AnchoredProfileField(
            [(5241.0, 155.0), (7660.0, 155.0), (7660.0, 155.0)], 0.85
        )
        assert len(f.md_anchor_m) == 2
        assert f.T(7660.0, 0.0) == pytest.approx(0.85 * 155.0)

    def test_fewer_than_two_anchors_raises(self):
        with pytest.raises(ValueError, match="至少需 2 点"):
            AnchoredProfileField([(7660.0, 155.0)], 0.85)
        with pytest.raises(ValueError):
            AnchoredProfileField([], 0.85)

    def test_nonfinite_anchor_raises(self):
        with pytest.raises(ValueError, match="NaN/Inf"):
            AnchoredProfileField([(1000.0, float("nan")), (3000.0, 90.0)], 0.85)

    def test_bad_regime_raises(self):
        with pytest.raises(ValueError, match="regime"):
            AnchoredProfileField(ANCHORS, 0.85, regime="transient")


# ---------------------------------------------------------------------------
# 协议同型：与既有三场同一调用形态（无 formal Protocol，duck 同型——测绘项14）
# ---------------------------------------------------------------------------
def _synthetic_table() -> TableTemperatureField:
    table = np.arange(12, dtype=float).reshape(3, 4)
    depth_m = np.array([0.0, 10.0, 20.0])
    time_s = np.arange(4, dtype=float) * 60.0
    return TableTemperatureField(table, depth_m, time_s, expected_shape=(3, 4))


def _all_fields():
    return [
        ConstantTemperatureField(T_c=42.5),
        GeothermalTemperatureField(),
        _synthetic_table(),
        AnchoredProfileField(ANCHORS, 0.85),
    ]


class TestProtocolIsomorphism:
    @pytest.mark.parametrize("field", _all_fields())
    def test_duck_interface_identical_shape(self, field):
        assert callable(field.T)
        assert isinstance(field.T(1500.0, 60.0), float)
        assert isinstance(field.oob_count, int)
        assert isinstance(field.oob_events, tuple)
        assert field.reset_audit() is None
        for ev in field.oob_events:
            assert isinstance(ev, ClampEvent)

    @pytest.mark.parametrize("field", _all_fields())
    def test_acceptable_by_solver_injection_shape(self, field):
        """模拟求解器消费形态：T(md, t) 每步标量查询 + 审计读取，不炸。"""
        for md in (0.0, 1500.0, 6000.0):
            for t in (0.0, 120.0):
                val = field.T(md, t)
                assert np.isfinite(val)
        assert field.oob_count >= 0


# ---------------------------------------------------------------------------
# 既有场类零改动旁证（全量回归 = test_temperature_field.py 同批重跑）
# ---------------------------------------------------------------------------
class TestExistingUnchanged:
    def test_module_surface_untouched(self):
        assert EXPECTED_SHAPE == (333, 200)
        assert GEO_T0_C == 16.006
        assert GEO_GRAD_C_PER_M == 1.7598e-2
        import cemdisp.data.temperature_field as tf

        for name in ("ConstantTemperatureField", "GeothermalTemperatureField",
                     "TableTemperatureField", "AnchoredProfileField",
                     "ClampEvent", "load_delivered_pair"):
            assert name in tf.__all__

    def test_constant_field_still_constant(self):
        f = ConstantTemperatureField(T_c=42.5)
        assert f.T(0.0, 0.0) == 42.5 and f.T(7660.0, 1e9) == 42.5
        assert f.oob_count == 0

    def test_geothermal_field_formula_still_canonical(self):
        f = GeothermalTemperatureField()
        assert f.T(5400.0, 0.0) == pytest.approx(16.006 + 1.7598e-2 * 5400.0)
        assert f.oob_count == 0 and f.oob_events == ()

    def test_table_field_scalar_clamp_behavior(self):
        f = _synthetic_table()
        assert f.T(-1000.0, f.time_s[2]) == pytest.approx(f.table[0, 2], abs=1e-12)
        assert f.oob_count == 1  # Table 标量路径仍逐查询记账（未变）
        # Phase 4d 主批 §4-1：Table 已补 F-4 聚合入口（原断言"无 T_column"是 4d 前期
        # 的状态快照，随该前置落地翻转为"存在且与标量路径逐位一致"）
        assert hasattr(f, "T_column")
        assert f.T_column([f.depth_m[1]], f.time_s[2])[0] == f.table[1, 2]
