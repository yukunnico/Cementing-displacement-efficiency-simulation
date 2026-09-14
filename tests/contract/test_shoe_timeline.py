"""
ShoeTimeline 抽象与测试

测试目标：
- ShoeEventKind 枚举：鞋口事件类型
- ShoeEvent 数据类：鞋口处单个事件快照
- ShoeTimeline 类：按时间轴管理鞋口事件，提供 .at(time_s) 查询

设计说明：
- ShoeTimeline.at() 在任意时刻返回最近的 PipeExitState
- 时间轴仅允许单调递增的事件（不支持乱序插入）
- 空时间轴在任意时刻返回零排量零相分数的默认状态
"""

import unittest
from cemdisp.transport1d.shoe_timeline import (
    ShoeEventKind,
    ShoeEvent,
    ShoeTimeline,
)
from cemdisp.transport1d.pipe_exit_state import PipeExitState


class TestShoeEventKind(unittest.TestCase):
    """测试 ShoeEventKind 枚举定义。"""

    def test_enum_values(self):
        """枚举包含 FRONT_ARRIVAL / REAR_EXIT / RATE_SWITCH / SHUTDOWN / RESTART / END。"""
        self.assertEqual(ShoeEventKind.FRONT_ARRIVAL.value, "FRONT_ARRIVAL")
        self.assertEqual(ShoeEventKind.REAR_EXIT.value, "REAR_EXIT")
        self.assertEqual(ShoeEventKind.RATE_SWITCH.value, "RATE_SWITCH")
        self.assertEqual(ShoeEventKind.SHUTDOWN.value, "SHUTDOWN")
        self.assertEqual(ShoeEventKind.RESTART.value, "RESTART")
        self.assertEqual(ShoeEventKind.END.value, "END")

    def test_enum_count(self):
        """枚举共有 6 个值。"""
        self.assertEqual(len(ShoeEventKind), 6)


class TestShoeEvent(unittest.TestCase):
    """测试 ShoeEvent 数据类。"""

    def test_fields(self):
        """ShoeEvent 包含 time_s, kind, flow_rate_m3_s, stage_name, phase_fractions。"""
        event = ShoeEvent(
            time_s=100.0,
            kind=ShoeEventKind.FRONT_ARRIVAL,
            flow_rate_m3_s=0.02,
            stage_name="注入领浆",
            phase_fractions=(("领浆", 1.0),),
        )
        self.assertEqual(event.time_s, 100.0)
        self.assertEqual(event.kind, ShoeEventKind.FRONT_ARRIVAL)
        self.assertEqual(event.flow_rate_m3_s, 0.02)
        self.assertEqual(event.stage_name, "注入领浆")
        self.assertEqual(event.phase_fractions, (("领浆", 1.0),))

    def test_phase_fractions_default_empty(self):
        """phase_fractions 默认为空元组。"""
        event = ShoeEvent(
            time_s=50.0,
            kind=ShoeEventKind.REAR_EXIT,
            flow_rate_m3_s=0.0,
            stage_name="停泵",
        )
        self.assertEqual(event.phase_fractions, ())


class TestShoeTimelineConstruction(unittest.TestCase):
    """测试 ShoeTimeline 构造与基本属性。"""

    def test_empty_timeline(self):
        """空时间轴 events 列表为空。"""
        tl = ShoeTimeline(events=[])
        self.assertEqual(tl.events, ())

    def test_single_event(self):
        """单事件时间轴。"""
        event = ShoeEvent(
            time_s=0.0,
            kind=ShoeEventKind.FRONT_ARRIVAL,
            flow_rate_m3_s=0.01,
            stage_name="启动",
            phase_fractions=(),
        )
        tl = ShoeTimeline(events=[event])
        self.assertEqual(len(tl.events), 1)
        self.assertEqual(tl.events[0].time_s, 0.0)

    def test_ascending_events_only(self):
        """乱序插入（时间递减）应抛出 ValueError。"""
        events = [
            ShoeEvent(time_s=10.0, kind=ShoeEventKind.FRONT_ARRIVAL, flow_rate_m3_s=0.0, stage_name=""),
            ShoeEvent(time_s=5.0, kind=ShoeEventKind.FRONT_ARRIVAL, flow_rate_m3_s=0.0, stage_name=""),
        ]
        with self.assertRaises(ValueError):
            ShoeTimeline(events=events)


class TestShoeTimelineAt(unittest.TestCase):
    """测试 ShoeTimeline.at() 查询行为。"""

    def _make_simple_timeline(self):
        """构建一个含三个事件的时间轴，测试边界行为。"""
        events = [
            ShoeEvent(
                time_s=0.0,
                kind=ShoeEventKind.FRONT_ARRIVAL,
                flow_rate_m3_s=0.0,
                stage_name="初始",
                phase_fractions=(),
            ),
            ShoeEvent(
                time_s=60.0,
                kind=ShoeEventKind.FRONT_ARRIVAL,
                flow_rate_m3_s=0.02,
                stage_name="注入领浆",
                phase_fractions=(("领浆", 1.0),),
            ),
            ShoeEvent(
                time_s=180.0,
                kind=ShoeEventKind.REAR_EXIT,
                flow_rate_m3_s=0.0,
                stage_name="停泵",
                phase_fractions=(),
            ),
        ]
        return ShoeTimeline(events=events)

    def test_at_before_first_event(self):
        """查询时刻早于首个事件：返回首个事件的状态（最近过去）。"""
        tl = self._make_simple_timeline()
        state = tl.at(time_s=-10.0)
        self.assertEqual(state.time_s, 0.0)
        self.assertEqual(state.flow_rate_m3_s, 0.0)
        self.assertEqual(state.stage_name, "初始")

    def test_at_first_event_exact(self):
        """查询时刻等于首个事件时间。"""
        tl = self._make_simple_timeline()
        state = tl.at(time_s=0.0)
        self.assertEqual(state.time_s, 0.0)
        self.assertEqual(state.flow_rate_m3_s, 0.0)

    def test_at_middle(self):
        """查询时刻在两个事件之间：返回最近过去事件的状态。"""
        tl = self._make_simple_timeline()
        state = tl.at(time_s=100.0)
        # 最近过去事件是 time_s=60.0 的事件
        self.assertEqual(state.time_s, 60.0)
        self.assertEqual(state.flow_rate_m3_s, 0.02)
        self.assertEqual(state.stage_name, "注入领浆")

    def test_at_last_event_exact(self):
        """查询时刻等于最后事件时间。"""
        tl = self._make_simple_timeline()
        state = tl.at(time_s=180.0)
        self.assertEqual(state.time_s, 180.0)
        self.assertEqual(state.flow_rate_m3_s, 0.0)

    def test_at_after_last_event(self):
        """查询时刻晚于最后事件：返回最后事件的状态。"""
        tl = self._make_simple_timeline()
        state = tl.at(time_s=300.0)
        self.assertEqual(state.time_s, 180.0)
        self.assertEqual(state.flow_rate_m3_s, 0.0)
        self.assertEqual(state.stage_name, "停泵")

    def test_at_exact_middle_event(self):
        """查询时刻恰好落在中间事件的时间点。"""
        tl = self._make_simple_timeline()
        state = tl.at(time_s=60.0)
        self.assertEqual(state.time_s, 60.0)
        self.assertEqual(state.flow_rate_m3_s, 0.02)


class TestShoeTimelineAtEmptyTimeline(unittest.TestCase):
    """测试空时间轴的 .at() 行为。"""

    def test_at_empty_timeline(self):
        """空时间轴任意时刻返回零排量默认 PipeExitState。"""
        tl = ShoeTimeline(events=[])
        state = tl.at(time_s=100.0)
        self.assertIsInstance(state, PipeExitState)
        self.assertEqual(state.time_s, 0.0)
        self.assertEqual(state.flow_rate_m3_s, 0.0)
        self.assertEqual(state.stage_name, "")


class TestShoeTimelinePhaseFractionsPassedThrough(unittest.TestCase):
    """测试 phase_fractions 原样从 ShoeEvent 传递到 PipeExitState。"""

    def test_phase_fractions_preserved(self):
        """ShoeEvent.phase_fractions 完整传递到 .at() 返回值。"""
        events = [
            ShoeEvent(
                time_s=0.0,
                kind=ShoeEventKind.FRONT_ARRIVAL,
                flow_rate_m3_s=0.01,
                stage_name="注入水泥",
                phase_fractions=(("水泥", 0.6), ("水", 0.4)),
            ),
        ]
        tl = ShoeTimeline(events=events)
        state = tl.at(time_s=50.0)
        self.assertEqual(state.phase_fractions, (("水泥", 0.6), ("水", 0.4)))


if __name__ == "__main__":
    unittest.main()