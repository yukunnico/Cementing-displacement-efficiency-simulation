# -*- coding: utf-8 -*-
"""管内段前缘到达次序保序（甲-2）回归契约。

不变量：连续泵入下，**界面到达次序 ≡ 注入次序**。

缺陷背景（2026-09-26 实测，HEAD d6fa344）
----------------------------------------
``_gravity_corrected_arrival_time`` 是**逐界面独立**的乘法修正，量级 ±1%~±8%，
足以让相邻界面对调：

1. ht1_004 隔离液2 修正后 7003.08 s 排到 领浆 修正后 6762.12 s 之后（原始活塞流
   时刻 6670.00 / 7070.00 本是递增）——鞋口时间线出现"后注入的流体先到"的非物理
   序列，环空入口浓度回退，物理合理性闸门 P-2 实测最大回退 1.000e+00（不通过）。
2. 同类第二处（哨兵口径，非保序约束）：ht1_004 压塞液的重力外推**延迟**
   11996.98 s 越过泵注结束时刻 11409.61 s，而"未到鞋口"的哨兵值原样取泵注
   结束时刻 ⇒ 哨兵反而排到压塞液之前。修法是**抬高哨兵**（取 max(泵注结束,
   已到达前缘最大值)），不压掉任何重力修正结果。

本测试把该不变量钉成**可证伪**的契约
------------------------------------
撤掉 ``_ordered_front_arrival_times`` 的保序约束（恢复"逐界面独立修正"，
即 :meth:`_Case.unconstrained` 走的那条老路径）后：

* ``test_ht1_004_lead_slurry_no_longer_overtakes_spacer2`` 转红；
* ``test_ht1_004_shoe_timeline_front_events_follow_injection_order`` 转红。

而 ``test_unaffected_wells_bitwise_equal_unconstrained`` 仍绿——约束在那些井上
本来就是恒等（该测试同时是"未受影响井逐位不变"的常驻守卫）。
``test_fixture_still_discriminates`` 反向钉住算例的判别力：若将来 ht1_004 的
loader/流体密度漂移到"对调消失"，该测试转红，避免保序契约变成零判别力的空过。
"""
import unittest

from cemdisp.data.loaders import (
    load_ht1_001_tailpipe,
    load_ht1_003_tailpipe,
    load_ht1_004_tailpipe,
    load_hu1_tailpipe,
    load_hu2_tailpipe,
    load_hu101_tailpipe,
    load_hu102_tailpipe,
    load_hu103_tailpipe,
)
from cemdisp.data.pumping_schedule import PumpingSchedule
from cemdisp.transport1d import CasingFlowSolver
from cemdisp.transport1d.shoe_timeline import ShoeEventKind

WELLS = (
    ("hu101", load_hu101_tailpipe),
    ("hu102", load_hu102_tailpipe),
    ("hu103", load_hu103_tailpipe),
    ("hu1", load_hu1_tailpipe),
    ("hu2", load_hu2_tailpipe),
    ("ht1_001", load_ht1_001_tailpipe),
    ("ht1_003", load_ht1_003_tailpipe),
    ("ht1_004", load_ht1_004_tailpipe),
)
# 与 ht1_004 无对调、保序约束应为恒等的七口井（不变量面）
UNAFFECTED = tuple((n, l) for n, l in WELLS if n != "ht1_004")


def _production_solver(*, enable_axial_dispersion: bool = True) -> CasingFlowSolver:
    """生产 1D 口径（T1 三开关 + 重力修正），与 verify_physical_sanity / runner 一致。"""

    return CasingFlowSolver(
        enable_gravity=True,
        mixing_contact_time=True,
        plug_face_zero_mixing=True,
        has_plug=True,
        enable_axial_dispersion=enable_axial_dispersion,
    )


class _Case:
    """一口井在**生产截断口径**下求解前缘所需的全套中间量。"""

    def __init__(self, solver: CasingFlowSolver, well, fluids, schedule) -> None:
        self.solver = solver
        self.well = well
        self.fluids = fluids
        self.schedule = schedule
        self.pipe_volume_m3 = solver._timeline_pipe_volume(
            well, well.shoe_md_m * solver._pipe_cross_section_area(well)
        )
        self.steps, self.pumping_end_time_s = solver._displacement_sequence_cutoff(
            solver._build_scheduled_steps(schedule)
        )
        self.initial_fluid = solver._initial_fluid_name(
            fluids, PumpingSchedule(steps=tuple(s.step for s in self.steps))
        ) if self.steps else ""
        self.raw_times = tuple(
            solver._front_arrival_time(s, self.steps, self.pipe_volume_m3)
            for s in self.steps
        )

    def ordered(self) -> tuple:
        """生产路径：重力修正 + 严格保序 + 泵注结束上限（唯一来源）。"""

        return self.solver._ordered_front_arrival_times(
            self.steps, self.pipe_volume_m3, self.fluids, self.initial_fluid, self.well,
        )

    def unconstrained(self) -> tuple:
        """撤掉保序约束的对照路径 = 修复前代码的逐界面独立重力修正。"""

        out = []
        for i, raw in enumerate(self.raw_times):
            if raw is None:
                out.append(None)
                continue
            out.append(self.solver._gravity_corrected_arrival_time(
                raw,
                self.steps[i].step.fluid_name,
                self.solver._displaced_fluid_name(self.steps, i, self.initial_fluid),
                self.fluids,
                self.well,
            ))
        return tuple(out)

    def fluid_names(self) -> tuple:
        return tuple(s.step.fluid_name for s in self.steps)

    def run_result(self):
        return self.solver.run(self.well, self.fluids, self.schedule)


def _case(loader) -> _Case:
    well, fluids, schedule, _ = loader()
    return _Case(_production_solver(), well, fluids, schedule)


def _crossings(times: tuple) -> list:
    """相邻注入步之间的到达次序违例（None = 未到达鞋口，不参与比较）。"""

    return [
        (i, a, b)
        for i, (a, b) in enumerate(zip(times, times[1:]))
        if a is not None and b is not None and a > b + 1e-12
    ]


class TestFrontOrderPreserving(unittest.TestCase):
    """八井生产口径下的到达次序不变量。"""

    def test_all_eight_wells_have_no_arrival_order_violation(self):
        """八井前缘按注入次序非减，且相邻界面间隔 ≥ dt（不得同刻退化）。"""

        total = 0
        for name, loader in WELLS:
            case = _case(loader)
            ordered = case.ordered()
            self.assertEqual([], _crossings(ordered), f"{name}: 注入次序出现对调")
            for i, (a, b) in enumerate(zip(ordered, ordered[1:])):
                if a is None or b is None:
                    continue
                self.assertGreaterEqual(
                    b - a, case.solver.dt - 1e-9,
                    f"{name}: 第 {i}/{i + 1} 界面间隔 {b - a} < dt={case.solver.dt}"
                    "（两界面同刻到达的退化）")
                total += 1
            # run() 暴露的 fronts 序列（含"未到鞋口"上界标记）不得递减
            fronts = tuple(f.time_s for f in case.run_result().fronts)
            self.assertEqual([], _crossings(fronts),
                             f"{name}: result.fronts 相邻到达次序违例")
        self.assertGreater(total, 0, "无任何相邻界面被检查（测试零判别力）")

    def test_ht1_004_lead_slurry_no_longer_overtakes_spacer2(self):
        """ht1_004：领浆被压到隔离液2 之后恰好 Δ_min = dt（保序约束生效）。"""

        case = _case(load_ht1_004_tailpipe)
        names = case.fluid_names()
        i_spacer2 = names.index("隔离液2")
        i_lead = i_spacer2 + 1
        self.assertEqual("领浆", names[i_lead])
        ordered = case.ordered()
        self.assertLess(ordered[i_spacer2], ordered[i_lead])
        self.assertAlmostEqual(
            ordered[i_lead] - ordered[i_spacer2], case.solver.dt, delta=1e-9,
            msg="领浆应被压到 隔离液2 + dt（失去全部修正提前量后仍保留最小间隔）")

    def test_unarrived_front_marks_never_precede_real_arrivals(self):
        """八井：未到鞋口的哨兵值 = max(泵注结束, 已到达前缘最大值)，fronts 非减。"""

        for name, loader in WELLS:
            case = _case(loader)
            real = [t for t in case.ordered() if t is not None]
            bound = max([case.pumping_end_time_s] + real)
            fronts = tuple(f.time_s for f in case.run_result().fronts)
            for i, t in enumerate(fronts):
                self.assertLessEqual(
                    t, bound + 1e-9, f"{name}: 第 {i} 个前缘 {t} 超出上界 {bound}")
            self.assertEqual(bound, max(fronts))
            self.assertEqual([], _crossings(fronts), f"{name}: fronts 出现次序违例")
        # 判别力守卫：ht1_004 是唯一"已到达前缘外推过泵注结束"的井；
        # 没有它，本测试对"哨兵固定取泵注结束时刻"这一旧口径无区分度。
        case = _case(load_ht1_004_tailpipe)
        real_max = max(t for t in case.ordered() if t is not None)
        self.assertGreater(real_max, case.pumping_end_time_s,
                           "ht1_004 已无前缘外推过泵注结束 ⇒ 本测试失去判别力")

    def test_unaffected_wells_bitwise_equal_unconstrained(self):
        """七口无对调的井：保序约束是恒等（输出与未约束重建**逐位相同**）。"""

        for name, loader in UNAFFECTED:
            case = _case(loader)
            self.assertEqual(
                case.unconstrained(), case.ordered(),
                f"{name}: 保序约束不应改动无对调井的任何一个比特")

    def test_fixture_still_discriminates(self):
        """算例判别力守卫：撤掉约束后 ht1_004 的对调必须**真的**复现。

        若将来 loader/流体密度漂移到对调消失，本测试转红——避免保序契约退化成
        零判别力的空过（vacuous green）。
        """

        case = _case(load_ht1_004_tailpipe)
        names = case.fluid_names()
        i_spacer2 = names.index("隔离液2")
        t_spacer2, t_lead = case.unconstrained()[i_spacer2:i_spacer2 + 2]
        self.assertGreater(t_spacer2, t_lead,
                           "ht1_004 未约束重建已无对调 ⇒ 本文件其余断言失去判别力")

    def test_ht1_004_shoe_timeline_front_events_follow_injection_order(self):
        """ht1_004 鞋口时间线的前缘事件次序 ≡ 注入次序（关弥散以单事件对单步）。"""

        well, fluids, schedule, _ = load_ht1_004_tailpipe()
        solver = _production_solver(enable_axial_dispersion=False)
        result = solver.run(well, fluids, schedule)
        case = _Case(solver, well, fluids, schedule)
        finite = [
            case.steps[i].step.fluid_name
            for i, t in enumerate(case.ordered()) if t is not None
        ]
        arrivals = [e for e in result.shoe_timeline.events
                    if e.kind == ShoeEventKind.FRONT_ARRIVAL]
        self.assertEqual(
            finite,
            [e.phase_fractions[0][0] for e in arrivals],
            "鞋口时间线的前缘事件未按注入次序排列（关弥散后应一个前缘一个事件）")
        times = [e.time_s for e in arrivals]
        self.assertEqual(sorted(times), times, "鞋口时间线前缘事件时刻未递增")

    def test_gravity_off_path_is_bitwise_raw(self):
        """关重力修正时，保序包装必须逐位返回原始活塞流时刻（零影响）。"""

        for name, loader in WELLS:
            well, fluids, schedule, _ = loader()
            solver = _production_solver()
            solver.enable_gravity = False
            case = _Case(solver, well, fluids, schedule)
            self.assertEqual(case.raw_times, case.ordered(),
                             f"{name}: 关重力路径被保序包装改动")


if __name__ == "__main__":
    unittest.main()
