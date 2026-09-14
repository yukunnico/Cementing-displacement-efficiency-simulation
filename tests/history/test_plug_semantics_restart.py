# STATUS: history —— 本文件锁定“重构前”的历史行为契约。
# 源模型口径重构（2026-09-14）后，相关断言预期变红；红灯不等同回归失败，
# 需按 docs/superpowers/plans/2026-09-14-d2dga-source-fidelity-and-repo-cleanup.md 逐条复核。
"""胶塞语义修复测试（2026-09-03 用户现场工艺裁定，TDD 先行）。

现场工艺事实（用户裁定）：
- 尾管固井尾浆之后有胶塞：替浆顶胶塞、胶塞驱动尾浆进环空，替浆本身不进环空；
- 次日后处理步（首个 RESTART 步起，如 hu102 循环排混浆）与水泥顶替过程解耦，
  不推动管内界面。

修复语义（cemdisp/transport1d/casing_flow.py）：
- "顶替序列" = 首个 RESTART 步之前的全部注入步骤；
- "顶替序列终点" = 首个 RESTART 步的 start_time（无 RESTART 井 = 泵注结束）；
- 界面推进（累计体积）截断到顶替序列之内：RESTART 步及之后的体积不计入累计，
  累计封顶后替浆前缘的到达体积坐标恒超出封顶值 → 替浆（及其后流体）永不到达
  鞋口，"替浆不进环空"由时间线构造自动保证；
- pumping_end_time_s（外推上界）= 顶替序列终点；
- cement_end_time_s = min(尾浆尾缘过鞋口时刻, 顶替序列终点)——替浆不足（未碰压）
  时停在替浆步末，尾浆尾段滞留管内；
- RESTART 步不生成鞋口事件（时间轴止于顶替序列终点）。

该截断是无条件工艺事实（不依赖 has_plug 开关）；has_plug 仅承载混浆增强因子=1
的语义。
"""

import math
import unittest

import pytest

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.pumping_schedule import (
    PumpingSchedule,
    PumpingScheduleStep,
    PumpingStageEvent,
)
from cemdisp.data.well_spec import WellSpec
from cemdisp.transport1d.casing_flow import CasingFlowSolver
from cemdisp.transport1d.shoe_timeline import ShoeEventKind


def _liner_id_for_area(area_m2: float) -> float:
    """按目标截面积反算内径，便于测试中直接控制管内容积。"""

    return math.sqrt(4.0 * area_m2 / math.pi) * 1000.0


def _fluids() -> tuple[FluidSpec, ...]:
    """构建测试用流体，名称与施工步骤保持一致。"""

    return (
        FluidSpec(name="泥浆", role=FluidRole.MUD, density_kg_m3=1200.0, plastic_viscosity_pa_s=0.02),
        FluidSpec(name="前置液", role=FluidRole.SPACER, density_kg_m3=1100.0, plastic_viscosity_pa_s=0.01),
        FluidSpec(name="尾浆", role=FluidRole.TAIL, density_kg_m3=1900.0, plastic_viscosity_pa_s=0.08),
        FluidSpec(name="压塞液", role=FluidRole.DISPLACEMENT, density_kg_m3=1000.0, plastic_viscosity_pa_s=0.01),
        FluidSpec(name="替浆液", role=FluidRole.DISPLACEMENT, density_kg_m3=1050.0, plastic_viscosity_pa_s=0.01),
        FluidSpec(name="钻井液", role=FluidRole.MUD, density_kg_m3=1150.0, plastic_viscosity_pa_s=0.02),
    )


def _well() -> WellSpec:
    """构建鞋深 100 m 的测试井；单一内径对应 1 m³ 管内容积。"""

    return WellSpec(
        well_name="胶塞语义测试井",
        top_md_m=1.0,
        bottom_md_m=120.0,
        shoe_md_m=100.0,
        liner_id_mm=_liner_id_for_area(0.01),
    )


def _schedule_insufficient_with_restart() -> PumpingSchedule:
    """替浆不足 + RESTART 步：[前置0.8 → 尾浆0.7 → 压塞0.5 → 替浆0.4 → RESTART循环1.0]。

    管内容积 1.0 m³。尾浆尾缘过鞋口需 1.5+1.0=2.5 m³ > 顶替序列累计 2.4 m³
    （差 0.1 m³）→ 未碰压，尾浆尾段 0.1 m³ 滞留管内。
    """

    return PumpingSchedule(
        steps=(
            PumpingScheduleStep(
                step_name="注入前置液", fluid_name="前置液", volume_m3=0.8,
                rate_m3_min=1.0, start_time_s=0.0, end_time_s=48.0,
                event_tag=PumpingStageEvent.INJECT_SPACER,
            ),
            PumpingScheduleStep(
                step_name="注入尾浆", fluid_name="尾浆", volume_m3=0.7,
                rate_m3_min=1.0, start_time_s=48.0, end_time_s=90.0,
                event_tag=PumpingStageEvent.INJECT_CEMENT,
            ),
            PumpingScheduleStep(
                step_name="注入压塞液", fluid_name="压塞液", volume_m3=0.5,
                rate_m3_min=1.0, start_time_s=90.0, end_time_s=120.0,
            ),
            PumpingScheduleStep(
                step_name="替浆", fluid_name="替浆液", volume_m3=0.4,
                rate_m3_min=1.0, start_time_s=120.0, end_time_s=144.0,
                event_tag=PumpingStageEvent.INJECT_DISPLACEMENT,
            ),
            PumpingScheduleStep(
                step_name="循环排混浆", fluid_name="钻井液", volume_m3=1.0,
                rate_m3_min=1.0, start_time_s=144.0, end_time_s=204.0,
                event_tag=PumpingStageEvent.RESTART,
            ),
        )
    )


def _schedule_sufficient_no_restart() -> PumpingSchedule:
    """替浆充足、无 RESTART：[前置0.8 → 尾浆0.7 → 压塞0.5 → 替浆0.9]。

    尾浆尾缘 2.5 m³ ≤ 累计 2.9 m³ → 尾浆尾缘在替浆步内（150 s）过鞋口。
    """

    return PumpingSchedule(
        steps=(
            PumpingScheduleStep(
                step_name="注入前置液", fluid_name="前置液", volume_m3=0.8,
                rate_m3_min=1.0, start_time_s=0.0, end_time_s=48.0,
                event_tag=PumpingStageEvent.INJECT_SPACER,
            ),
            PumpingScheduleStep(
                step_name="注入尾浆", fluid_name="尾浆", volume_m3=0.7,
                rate_m3_min=1.0, start_time_s=48.0, end_time_s=90.0,
                event_tag=PumpingStageEvent.INJECT_CEMENT,
            ),
            PumpingScheduleStep(
                step_name="注入压塞液", fluid_name="压塞液", volume_m3=0.5,
                rate_m3_min=1.0, start_time_s=90.0, end_time_s=120.0,
            ),
            PumpingScheduleStep(
                step_name="替浆", fluid_name="替浆液", volume_m3=0.9,
                rate_m3_min=1.0, start_time_s=120.0, end_time_s=174.0,
                event_tag=PumpingStageEvent.INJECT_DISPLACEMENT,
            ),
        )
    )


def _schedule_sufficient_with_restart() -> PumpingSchedule:
    """替浆充足 + RESTART 步：替浆 1.0（120-180）→ RESTART 循环 0.5（180-210）。

    尾浆尾缘 2.5 m³ 在替浆步内（150 s）过鞋口；替浆液前缘恰好于替浆步末（180 s）
    到达鞋口。 cement_end 必须与旧行为（150 s）逐位一致。
    """

    return PumpingSchedule(
        steps=(
            PumpingScheduleStep(
                step_name="注入前置液", fluid_name="前置液", volume_m3=0.8,
                rate_m3_min=1.0, start_time_s=0.0, end_time_s=48.0,
                event_tag=PumpingStageEvent.INJECT_SPACER,
            ),
            PumpingScheduleStep(
                step_name="注入尾浆", fluid_name="尾浆", volume_m3=0.7,
                rate_m3_min=1.0, start_time_s=48.0, end_time_s=90.0,
                event_tag=PumpingStageEvent.INJECT_CEMENT,
            ),
            PumpingScheduleStep(
                step_name="注入压塞液", fluid_name="压塞液", volume_m3=0.5,
                rate_m3_min=1.0, start_time_s=90.0, end_time_s=120.0,
            ),
            PumpingScheduleStep(
                step_name="替浆", fluid_name="替浆液", volume_m3=1.0,
                rate_m3_min=1.0, start_time_s=120.0, end_time_s=180.0,
                event_tag=PumpingStageEvent.INJECT_DISPLACEMENT,
            ),
            PumpingScheduleStep(
                step_name="循环排混浆", fluid_name="钻井液", volume_m3=0.5,
                rate_m3_min=1.0, start_time_s=180.0, end_time_s=210.0,
                event_tag=PumpingStageEvent.RESTART,
            ),
        )
    )


class TestPlugSemanticsRestart(unittest.TestCase):
    """胶塞语义：RESTART 后处理步解耦 + 替浆不进环空。"""

    def setUp(self) -> None:
        self.solver = CasingFlowSolver(
            enable_gravity=False, enable_axial_dispersion=False
        )

    # ------------------------------------------------------------------
    # a) 替浆不足 + RESTART 步
    # ------------------------------------------------------------------

    def test_insufficient_displacement_restarts_do_not_push_interface(self) -> None:
        """替浆不足时 cement_end=替浆步末；RESTART 步不产生鞋口事件、不推进界面。"""

        result = self.solver.run(_well(), _fluids(), _schedule_insufficient_with_restart())

        # 顶替序列终点 = RESTART 步 start = 144 s = 替浆步末
        self.assertEqual(result.pumping_end_time_s, 144.0)
        # 尾浆尾缘需 2.5 m³ > 顶替序列累计 2.4 m³ → 未过鞋口，cement_end 停在替浆步末
        self.assertEqual(result.cement_end_time_s, 144.0)
        # 界面推进截断：fronts 只含顶替序列 4 步（RESTART 步不建前缘）
        self.assertEqual(len(result.fronts), 4)
        self.assertEqual([f.time_s for f in result.fronts], [60.0, 108.0, 144.0, 144.0])

        # RESTART 步不生成鞋口事件：时间轴止于顶替序列终点
        events = result.shoe_timeline.events
        self.assertEqual(events[-1].time_s, 144.0)
        self.assertEqual(events[-1].kind, ShoeEventKind.END)
        self.assertTrue(all(e.time_s <= 144.0 + 1e-9 for e in events))
        self.assertNotIn(ShoeEventKind.RESTART, [e.kind for e in events])
        self.assertFalse(any("循环排混浆" in e.stage_name for e in events))

        # RESTART 窗口内（144-204 s）查询：界面停滞——零排量、鞋口仍为尾浆（尾段滞留）
        state = result.shoe_timeline.at(150.0)
        self.assertEqual(state.flow_rate_m3_s, 0.0)
        self.assertEqual(state.phase_fractions, (("尾浆", 1.0),))
        self.assertEqual(state.stage_name, "施工结束后保持")
        # 与 pipe_exit_state_at（legacy 查询）同口径
        legacy = self.solver.pipe_exit_state_at(result, 150.0)
        self.assertEqual(legacy.phase_fractions, state.phase_fractions)
        self.assertEqual(legacy.flow_rate_m3_s, state.flow_rate_m3_s)

    def test_restart_decoupling_is_unconditional_regardless_of_has_plug(self) -> None:
        """RESTART 解耦是无条件工艺事实：has_plug 开/关结果一致（不引入新开关）。"""

        well, fluids, schedule = _well(), _fluids(), _schedule_insufficient_with_restart()
        r_off = CasingFlowSolver(enable_gravity=False, enable_axial_dispersion=False, has_plug=False).run(well, fluids, schedule)
        r_on = CasingFlowSolver(enable_gravity=False, enable_axial_dispersion=False, has_plug=True).run(well, fluids, schedule)
        self.assertEqual(r_off.cement_end_time_s, r_on.cement_end_time_s)
        self.assertEqual(r_off.pumping_end_time_s, r_on.pumping_end_time_s)

    # ------------------------------------------------------------------
    # b) 替浆充足（回归保护：与旧行为逐位一致）
    # ------------------------------------------------------------------

    def test_sufficient_displacement_no_restart_bit_identical_to_legacy(self) -> None:
        """无 RESTART 井（碰压成功）：cement_end/pumping_end/fronts 与旧行为逐位一致。"""

        result = self.solver.run(_well(), _fluids(), _schedule_sufficient_no_restart())

        # 旧口径（全日程最大结束时刻）：无 RESTART → 顶替序列终点 = 泵注结束
        self.assertEqual(result.pumping_end_time_s, 174.0)
        # 尾浆尾缘 2.5 m³ 在替浆步内（offset 0.5 m³ @1.0 m³/min）→ 150 s 过鞋口
        self.assertEqual(result.cement_end_time_s, 150.0)
        self.assertEqual([f.time_s for f in result.fronts], [60.0, 108.0, 150.0, 174.0])

        events = result.shoe_timeline.events
        self.assertEqual(events[-1].time_s, 174.0)
        self.assertEqual(events[-1].kind, ShoeEventKind.END)
        # 尾浆尾缘 REAR_EXIT 存在且在 cement_end 时刻
        rear_events = [e for e in events if e.kind == ShoeEventKind.REAR_EXIT]
        self.assertTrue(any(e.time_s == pytest.approx(150.0) for e in rear_events))

    def test_sufficient_displacement_with_restart_keeps_cement_end(self) -> None:
        """替浆充足 + RESTART：cement_end 与旧行为逐位一致（150 s），仅时间轴尾部截断。"""

        result = self.solver.run(_well(), _fluids(), _schedule_sufficient_with_restart())

        self.assertEqual(result.cement_end_time_s, 150.0)
        self.assertEqual(result.pumping_end_time_s, 180.0)  # 顶替序列终点（旧口径为 210）
        # 前缘与旧口径逐位一致（顶替序列内累计体积未变）
        self.assertEqual([f.time_s for f in result.fronts], [60.0, 108.0, 150.0, 180.0])

        events = result.shoe_timeline.events
        self.assertNotIn(ShoeEventKind.RESTART, [e.kind for e in events])
        self.assertEqual(events[-1].time_s, 180.0)
        # 顶替序列内的鞋口状态不受影响
        self.assertEqual(result.shoe_timeline.at(150.0).phase_fractions, (("压塞液", 1.0),))

    # ------------------------------------------------------------------
    # c) hu102 真实 loader 回归
    # ------------------------------------------------------------------

    def test_hu102_restart_circulation_does_not_push_interface(self) -> None:
        """hu102：循环排混浆（RESTART，41 m³）不再把尾浆尾缘推出鞋口。

        0708 核实事实：替浆序列累计 139 m³ < 尾浆尾缘过鞋口所需 52+88.9=140.9 m³
        （差 1.9 m³，单流阀失效、到量未碰压）→ cement_end = 替浆步末（= 循环排混浆
        start 时刻），尾浆尾段 1.9 m³ 滞留管内。
        """

        from cemdisp.data.loaders import load_hu102_tailpipe

        well, fluids, schedule, _ = load_hu102_tailpipe()
        solver = CasingFlowSolver(enable_gravity=False)
        result = solver.run(well, fluids, schedule)

        disp_steps = schedule.steps[:8]  # 前 8 步 = 顶替序列（第 9 步循环排混浆为 RESTART）
        disp_total_m3 = sum(s.volume_m3 for s in disp_steps)
        disp_end_s = sum(s.volume_m3 / s.rate_m3_min * 60.0 for s in disp_steps)
        assert disp_total_m3 == pytest.approx(139.0)
        assert 52.0 + well.shoe_lag_volume_m3 - disp_total_m3 == pytest.approx(1.9)

        assert result.pumping_end_time_s == pytest.approx(disp_end_s, abs=1e-6)
        assert result.cement_end_time_s == pytest.approx(disp_end_s, abs=1e-6)

        # RESTART 步不生成鞋口事件；时间轴止于替浆步末
        events = result.shoe_timeline.events
        assert all(e.time_s <= disp_end_s + 1e-9 for e in events)
        assert not any(e.kind == ShoeEventKind.RESTART for e in events)
        assert not any("循环排混浆" in e.stage_name for e in events)

        # RESTART 窗口内查询：界面停滞（尾浆滞留鞋口，零排量）
        state = result.shoe_timeline.at(disp_end_s + 100.0)
        assert state.flow_rate_m3_s == pytest.approx(0.0)
        assert state.phase_fractions == (("尾管水泥浆", 1.0),)

        # 默认重力口径（生产 runner 口径）同样停在替浆步末
        result_g = CasingFlowSolver().run(well, fluids, schedule)
        assert result_g.cement_end_time_s == pytest.approx(disp_end_s, abs=1e-6)
        assert result_g.pumping_end_time_s == pytest.approx(disp_end_s, abs=1e-6)
        # 修复注记写入 result.notes
        assert any("胶塞" in note for note in result_g.notes)


if __name__ == "__main__":
    _ = unittest.main()
