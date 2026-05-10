"""
PumpingSchedule 扩展功能测试

测试新增功能：
- PumpingStageEvent 枚举：施工作业阶段标签
- PumpingScheduleStep.event_tag: Optional[PumpingStageEvent] 类型字段
- PumpingSchedule.total_injected_volume_m3: 所有步骤注入体积之和
- PumpingSchedule.cement_phase_steps: 按 INJECT_CEMENT 过滤的步骤元组
"""

import unittest
from cemdisp.data.pumping_schedule import (
    PumpingScheduleStep,
    PumpingSchedule,
    PumpingStageEvent,
)


class TestPumpingStageEventEnum(unittest.TestCase):
    """测试 PumpingStageEvent 枚举定义。"""

    def test_enum_values(self):
        """枚举包含所有要求的事件类型。"""
        self.assertEqual(PumpingStageEvent.INJECT_CEMENT.value, "INJECT_CEMENT")
        self.assertEqual(PumpingStageEvent.INJECT_SPACER.value, "INJECT_SPACER")
        self.assertEqual(PumpingStageEvent.INJECT_DISPLACEMENT.value, "INJECT_DISPLACEMENT")
        self.assertEqual(PumpingStageEvent.PLUG_RELEASE.value, "PLUG_RELEASE")
        self.assertEqual(PumpingStageEvent.SHUTDOWN.value, "SHUTDOWN")
        self.assertEqual(PumpingStageEvent.RESTART.value, "RESTART")
        self.assertEqual(PumpingStageEvent.RATE_SWITCH.value, "RATE_SWITCH")

    def test_enum_count(self):
        """枚举共有 7 个值。"""
        self.assertEqual(len(PumpingStageEvent), 7)


class TestPumpingScheduleStepEventTag(unittest.TestCase):
    """测试 PumpingScheduleStep.event_tag 类型为 Optional[PumpingStageEvent]。"""

    def test_event_tag_optional(self):
        """event_tag 不提供时为 None。"""
        step = PumpingScheduleStep(
            step_name="注入领浆",
            fluid_name="领浆",
            volume_m3=5.0,
            rate_m3_min=0.5,
        )
        self.assertIsNone(step.event_tag)

    def test_event_tag_inject_cement(self):
        """event_tag 可接受 INJECT_CEMENT 枚举值。"""
        step = PumpingScheduleStep(
            step_name="注入领浆",
            fluid_name="领浆",
            volume_m3=5.0,
            rate_m3_min=0.5,
            event_tag=PumpingStageEvent.INJECT_CEMENT,
        )
        self.assertEqual(step.event_tag, PumpingStageEvent.INJECT_CEMENT)

    def test_event_tag_inject_spacer(self):
        """event_tag 可接受 INJECT_SPACER 枚举值。"""
        step = PumpingScheduleStep(
            step_name="注入前置液",
            fluid_name="前置液",
            volume_m3=2.0,
            rate_m3_min=0.3,
            event_tag=PumpingStageEvent.INJECT_SPACER,
        )
        self.assertEqual(step.event_tag, PumpingStageEvent.INJECT_SPACER)

    def test_event_tag_none_explicit(self):
        """event_tag 可显式传 None（表示无阶段标签）。"""
        step = PumpingScheduleStep(
            step_name="替泥浆",
            fluid_name="泥浆",
            volume_m3=10.0,
            rate_m3_min=1.0,
            event_tag=None,
        )
        self.assertIsNone(step.event_tag)

    def test_event_tag_invalid_string_raises(self):
        """event_tag 传入普通字符串（非 PumpingStageEvent 枚举值）应抛出 ValueError。"""
        # 使用 kwargs 字典构造以绕过静态类型检查，确保触发 __post_init__ 中的运行时校验
        kwargs = {
            "step_name": "注入领浆",
            "fluid_name": "领浆",
            "volume_m3": 5.0,
            "rate_m3_min": 0.5,
            "event_tag": "lead_cement",  # plain string, not allowed — tests runtime path
        }
        with self.assertRaises(ValueError):
            PumpingScheduleStep(**kwargs)


class TestPumpingScheduleTotalInjectedVolume(unittest.TestCase):
    """测试 PumpingSchedule.total_injected_volume_m3 —— 所有步骤体积之和。"""

    def test_total_injected_volume_m3_all_steps(self):
        """total_injected_volume_m3 是所有步骤 volume_m3 之和，不做任何流体过滤。"""
        schedule = PumpingSchedule(
            steps=(
                PumpingScheduleStep(
                    step_name="注入领浆",
                    fluid_name="领浆",
                    volume_m3=10.0,
                    rate_m3_min=1.0,
                    event_tag=PumpingStageEvent.INJECT_CEMENT,
                ),
                PumpingScheduleStep(
                    step_name="替泥浆",
                    fluid_name="泥浆",
                    volume_m3=20.0,
                    rate_m3_min=2.0,
                    event_tag=None,
                ),
                PumpingScheduleStep(
                    step_name="注入尾管浆",
                    fluid_name="尾管浆",
                    volume_m3=5.0,
                    rate_m3_min=0.5,
                    event_tag=PumpingStageEvent.INJECT_CEMENT,
                ),
            )
        )
        # 10 + 20 + 5 = 35（全部体积相加，无过滤）
        self.assertEqual(schedule.total_injected_volume_m3, 35.0)

    def test_total_injected_volume_m3_single_step(self):
        """单步骤时等于该步骤体积。"""
        schedule = PumpingSchedule(
            steps=(
                PumpingScheduleStep(
                    step_name="注入水泥浆",
                    fluid_name="水泥浆",
                    volume_m3=8.0,
                    rate_m3_min=0.8,
                    event_tag=PumpingStageEvent.INJECT_CEMENT,
                ),
            )
        )
        self.assertEqual(schedule.total_injected_volume_m3, 8.0)


class TestPumpingScheduleCementPhaseSteps(unittest.TestCase):
    """测试 PumpingSchedule.cement_phase_steps —— 按 INJECT_CEMENT 过滤。"""

    def test_cement_phase_steps_by_event_tag(self):
        """cement_phase_steps 只返回 event_tag == INJECT_CEMENT 的步骤。"""
        schedule = PumpingSchedule(
            steps=(
                PumpingScheduleStep(
                    step_name="注入领浆",
                    fluid_name="领浆",
                    volume_m3=10.0,
                    rate_m3_min=1.0,
                    event_tag=PumpingStageEvent.INJECT_CEMENT,
                ),
                PumpingScheduleStep(
                    step_name="替泥浆",
                    fluid_name="泥浆",
                    volume_m3=20.0,
                    rate_m3_min=2.0,
                    event_tag=None,
                ),
                PumpingScheduleStep(
                    step_name="注入尾管浆",
                    fluid_name="尾管浆",
                    volume_m3=5.0,
                    rate_m3_min=0.5,
                    event_tag=PumpingStageEvent.INJECT_CEMENT,
                ),
            )
        )
        cement_steps = schedule.cement_phase_steps
        self.assertEqual(len(cement_steps), 2)
        self.assertEqual(cement_steps[0].step_name, "注入领浆")
        self.assertEqual(cement_steps[1].step_name, "注入尾管浆")

    def test_cement_phase_steps_no_cement(self):
        """没有 INJECT_CEMENT 步骤时返回空元组。"""
        schedule = PumpingSchedule(
            steps=(
                PumpingScheduleStep(
                    step_name="替泥浆",
                    fluid_name="泥浆",
                    volume_m3=20.0,
                    rate_m3_min=2.0,
                    event_tag=None,
                ),
                PumpingScheduleStep(
                    step_name="盐水顶替",
                    fluid_name="盐水",
                    volume_m3=30.0,
                    rate_m3_min=3.0,
                    event_tag=PumpingStageEvent.INJECT_DISPLACEMENT,
                ),
            )
        )
        cement_steps = schedule.cement_phase_steps
        self.assertEqual(len(cement_steps), 0)

    def test_cement_phase_steps_all_cement(self):
        """全部为 INJECT_CEMENT 步骤时返回全部。"""
        schedule = PumpingSchedule(
            steps=(
                PumpingScheduleStep(
                    step_name="注入领浆",
                    fluid_name="领浆",
                    volume_m3=10.0,
                    rate_m3_min=1.0,
                    event_tag=PumpingStageEvent.INJECT_CEMENT,
                ),
                PumpingScheduleStep(
                    step_name="注入尾管浆",
                    fluid_name="尾管浆",
                    volume_m3=5.0,
                    rate_m3_min=0.5,
                    event_tag=PumpingStageEvent.INJECT_CEMENT,
                ),
            )
        )
        cement_steps = schedule.cement_phase_steps
        self.assertEqual(len(cement_steps), 2)

    def test_cement_phase_steps_ignores_fluid_name_with_none_tag(self):
        """
        cement_phase_steps 仅按 event_tag 过滤。

        即使 fluid_name 含"浆"（如领浆），若 event_tag=None 也不应被包含。
        这证明过滤基于 event_tag 而非 fluid_name 字符串heuristic。
        """
        schedule = PumpingSchedule(
            steps=(
                PumpingScheduleStep(
                    step_name="注入领浆",
                    fluid_name="领浆",          # 水泥浆名称，但无阶段标签
                    volume_m3=10.0,
                    rate_m3_min=1.0,
                    event_tag=None,
                ),
                PumpingScheduleStep(
                    step_name="盐水顶替",
                    fluid_name="盐水",          # 非水泥浆名称，但有 INJECT_CEMENT 标签
                    volume_m3=20.0,
                    rate_m3_min=2.0,
                    event_tag=PumpingStageEvent.INJECT_CEMENT,
                ),
            )
        )
        cement_steps = schedule.cement_phase_steps
        # 只有 INJECT_CEMENT 的步骤应被包含，与 fluid_name 无关
        self.assertEqual(len(cement_steps), 1)
        self.assertEqual(cement_steps[0].fluid_name, "盐水")
        self.assertEqual(cement_steps[0].step_name, "盐水顶替")

    def test_cement_phase_steps_non_cement_fluid_name_with_inject_cement_tag(self):
        """
        cement_phase_steps 基于 event_tag 而非 fluid_name。

        fluid_name 为"位移液"（非"浆"类），但 event_tag=INJECT_CEMENT，
        因此应被包含。这证明过滤机制是 event_tag 而非 fluid_name 字符串判断。
        """
        schedule = PumpingSchedule(
            steps=(
                PumpingScheduleStep(
                    step_name="位移液驱替",
                    fluid_name="位移液",
                    volume_m3=15.0,
                    rate_m3_min=1.5,
                    event_tag=PumpingStageEvent.INJECT_CEMENT,
                ),
                PumpingScheduleStep(
                    step_name="替泥浆",
                    fluid_name="泥浆",
                    volume_m3=20.0,
                    rate_m3_min=2.0,
                    event_tag=None,
                ),
            )
        )
        cement_steps = schedule.cement_phase_steps
        # 位移液有 INJECT_CEMENT 标签，应被包含；泥浆无标签，不应被包含
        self.assertEqual(len(cement_steps), 1)
        self.assertEqual(cement_steps[0].fluid_name, "位移液")


if __name__ == "__main__":
    unittest.main()
