"""
接触时间积分 σ_t（路线 B Task 1）测试

mixing_contact_time 开关（默认 False）：
- 关：σ_t 走原路径（全程行程时间近似），hu103 shoe_timeline 逐事件等于基线锚；
- 开：σ_t 改由界面真实接触时间历程积分给出
  （t_contact = t_arrival − t_inject，σ_t = sqrt(2·D_eff·t_contact)/U），
  D_eff/混浆增强/防御上下限/n_sub=5/F4 收尾结构全部不变。

测试项：
1. test_default_off_bitwise        —— 双开关默认，逐位复现基线锚；
2. test_constant_rate_anchor       —— 恒排量+单内径+等密度合成井，σ_on == σ_off
                                       （规格渐近一致性锚点：x_inj=0 首界面接触
                                       时间=全程，恒速井上全界面逐位一致）；
3. test_variable_rate_directions   —— (i) 减速合成井 σ_on < σ_off（受控收窄）；
                                       (ii) hu103 真实井全界面 σ_on > σ_off
                                       （双内径体积链差被接触时间修正暴露），
                                       尾浆界面 pin σ_on/σ_off == sqrt(88.55/V_liner)；
4. test_segmented_injection_boundary —— 同流体分两段注入，界面归着到
                                       该前缘所在注入步的出发时刻。

Task 2（胶塞面零掺混 plug_face_zero_mixing，默认 False）：
5. test_plug_face_no_dispersion    —— ht1_003（loader 生产口径 has_plug=False，
                                       单测属性覆盖 has_plug=True）开关开：
                                       "尾浆→压塞液"界面无过渡子事件（仅原
                                       阶跃单相事件），其余 5 界面过渡带逐位
                                       不受影响；
6. test_plug_face_off_bitwise      —— ht1_003 生产口径（has_plug 默认 False）
                                       仅开 plug_face_zero_mixing：与全关参照
                                       逐位相同（has_plug 门空转，零效应）；
7. test_hu102_unaffected           —— hu102（has_plug=True + 开关开）：其
                                       界面集无压塞液界面（压塞液前缘被
                                       RESTART 截断挡在鞋口上游），开关开 vs 关
                                       shoe_timeline 逐位相同（保护性断言）。
"""

import json
import math
import os
import unittest

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.pumping_schedule import (
    PumpingSchedule,
    PumpingScheduleStep,
    PumpingStageEvent,
)
from cemdisp.data.well_spec import WellSpec
from cemdisp.transport1d.casing_flow import CasingFlowSolver
from cemdisp.transport1d.shoe_timeline import ShoeEventKind

_BASELINE_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "_baseline_shoe_timeline_hu103.json",
)


def _liner_id_for_area(area_m2: float) -> float:
    """按目标截面积反算内径，便于测试中直接控制管内容积。"""

    return math.sqrt(4.0 * area_m2 / math.pi) * 1000.0


def _hu103_case():
    """加载 hu103 真实井（设计版）。"""

    from cemdisp.data.loaders import load_hu103_tailpipe

    return load_hu103_tailpipe()


def _ht1_003_case():
    """加载 ht1_003 真实井（实际施工版，loader 生产口径 has_plug=False）。"""

    from cemdisp.data.loaders import load_ht1_003_tailpipe

    return load_ht1_003_tailpipe()


def _hu102_case():
    """加载 hu102 真实井（唯一含 RESTART 步的重建现场程序）。"""

    from cemdisp.data.loaders import load_hu102_tailpipe

    return load_hu102_tailpipe()


def _synthetic_fluids(mud_rho: float = 1200.0, spacer_rho: float | None = None,
                      cement_rho: float = 1900.0, disp_rho: float = 1050.0):
    """合成井流体：泥浆/隔离液/尾浆/顶替液；spacer_rho=None 时与泥浆等密度。"""

    return (
        FluidSpec(name="泥浆", role=FluidRole.MUD, density_kg_m3=mud_rho,
                  plastic_viscosity_pa_s=0.02),
        FluidSpec(name="隔离液", role=FluidRole.SPACER,
                  density_kg_m3=spacer_rho if spacer_rho is not None else mud_rho,
                  plastic_viscosity_pa_s=0.01),
        FluidSpec(name="尾浆", role=FluidRole.TAIL, density_kg_m3=cement_rho,
                  plastic_viscosity_pa_s=0.08),
        FluidSpec(name="顶替液", role=FluidRole.DISPLACEMENT, density_kg_m3=disp_rho,
                  plastic_viscosity_pa_s=0.01),
    )


def _synthetic_well(shoe_md_m: float = 100.0, area_m2: float = 0.01) -> WellSpec:
    """合成单内径井：无 shoe_lag_volume_m3、无双内径字段，
    _timeline_pipe_volume 回退到截面口径（shoe_md × area）。"""

    return WellSpec(
        well_name="合成井",
        top_md_m=1.0,
        bottom_md_m=shoe_md_m + 20.0,
        shoe_md_m=shoe_md_m,
        liner_id_mm=_liner_id_for_area(area_m2),
    )


def _dispersion_bands(events):
    """提取弥散过渡带：[(next_fluid, prev_fluid, sigma_t, t_center), ...]。

    过渡带由 5 个子事件（t_arr−σ → t_arr+σ）+ F4 收尾事件（与最后子事件同刻、
    frac=1.0）构成；σ_t = (该相组末子事件时刻 − 首子事件时刻)/2。

    按 (next, prev) 相组跨全事件列表分组（同一过渡带的 5 子事件可能被
    REAR_EXIT/RATE_SWITCH 等同时刻事件夹断，不要求连续）；只统计 frac<1.0
    的子事件（F4 收尾 frac=1.0 与末子事件同刻，排除后首末跨度=2σ）。
    """

    by_pair: dict[tuple[str, str], list[float]] = {}
    for e in events:
        if e.kind != ShoeEventKind.FRONT_ARRIVAL or len(e.phase_fractions) != 2:
            continue
        if e.phase_fractions[0][1] >= 1.0:
            continue  # F4 收尾/正体事件
        key = (e.phase_fractions[0][0], e.phase_fractions[1][0])
        by_pair.setdefault(key, []).append(e.time_s)
    bands = []
    for (nxt, prv), times in by_pair.items():
        times.sort()
        sigma = (times[-1] - times[0]) / 2.0
        center = (times[-1] + times[0]) / 2.0
        bands.append((nxt, prv, sigma, center))
    return bands


class TestMixingContactTime(unittest.TestCase):
    """mixing_contact_time 开关行为测试。"""

    # ------------------------------------------------------------------
    # 1. 默认关：逐位复现基线锚
    # ------------------------------------------------------------------
    def test_default_off_bitwise(self) -> None:
        """双开关默认（mixing_contact_time 不传），hu103 shoe_timeline
        逐事件等于基线锚 JSON（time_s/kind/flow_rate/stage/phase_fractions
        全字段 == 级对比），cement_end_time_s 亦逐位一致。"""

        with open(_BASELINE_PATH, encoding="utf-8") as f:
            anchor = json.load(f)

        well, fluids, schedule, _ = _hu103_case()
        solver = CasingFlowSolver(enable_gravity=True)
        self.assertFalse(solver.mixing_contact_time)  # 默认关
        result = solver.run(well, fluids, schedule)

        events = result.shoe_timeline.events
        base = anchor["events"]
        self.assertEqual(anchor["_meta"]["event_count"], len(base))
        self.assertEqual(len(base), len(events))
        for i, (ev, ref) in enumerate(zip(events, base)):
            with self.subTest(event_index=i):
                self.assertEqual(ev.time_s, ref["time_s"])
                self.assertEqual(ev.kind.value, ref["kind"])
                self.assertEqual(ev.flow_rate_m3_s, ref["flow_rate_m3_s"])
                self.assertEqual(ev.stage_name, ref["stage_name"])
                self.assertEqual(
                    [list(pair) for pair in ev.phase_fractions], ref["phase_fractions"]
                )
        self.assertEqual(result.cement_end_time_s, anchor["_meta"]["cement_end_time_s"])

    # ------------------------------------------------------------------
    # 2. 恒排量锚点：σ_on == σ_off（规格渐近一致性）
    # ------------------------------------------------------------------
    def test_constant_rate_anchor(self) -> None:
        """恒排量 + 单内径 + 等密度隔离液合成井：开关开 vs 关，全部界面
        σ_on == σ_off（相对差 ≤1e-12）。

        规格 §146 锚点：x_inj=0 首界面接触时间=全程行程时间，与旧公式在恒速
        井上逐位一致；恒速 + 单内径（_timeline_pipe_volume 回退截面口径，
        到达链与 t_travel 链同体积）下该一致性推广到全部界面。
        重力修正关闭（等密度+单内径合成井，排除密度差修正的时间缩放干扰，
        使 t_arrival 与 t_inject 链在两条路径下严格同源）。
        """

        shoe_md_m = 100.0
        area_m2 = 0.01  # 管容 = 1.0 m³
        well = _synthetic_well(shoe_md_m=shoe_md_m, area_m2=area_m2)
        fluids = _synthetic_fluids(spacer_rho=None)  # 等密度隔离液
        schedule = PumpingSchedule(
            steps=(
                PumpingScheduleStep(
                    step_name="注入隔离液", fluid_name="隔离液", volume_m3=0.4,
                    rate_m3_min=1.0, start_time_s=0.0, end_time_s=24.0,
                    event_tag=PumpingStageEvent.INJECT_SPACER,
                ),
                PumpingScheduleStep(
                    step_name="注入尾浆", fluid_name="尾浆", volume_m3=1.0,
                    rate_m3_min=1.0, start_time_s=24.0, end_time_s=84.0,
                    event_tag=PumpingStageEvent.INJECT_CEMENT,
                ),
                PumpingScheduleStep(
                    step_name="替浆", fluid_name="顶替液", volume_m3=2.0,
                    rate_m3_min=1.0, start_time_s=84.0, end_time_s=204.0,
                    event_tag=PumpingStageEvent.INJECT_DISPLACEMENT,
                ),
            )
        )

        # dt=0.05（非默认 2.0）：默认 dt 下"隔离液←泥浆"界面 raw σ≈1.008s 会被
        # max(σ, dt) 下限 clamp 成 2.0 平凡相等（无区分度）；dt=0.05 下 3 界面
        # raw σ 全部 > dt，断言对所有界面都有区分度（Task 1 审查 F1）。
        solver_off = CasingFlowSolver(enable_gravity=False, mixing_contact_time=False, dt=0.05)
        solver_on = CasingFlowSolver(enable_gravity=False, mixing_contact_time=True, dt=0.05)
        events_off = solver_off.run(well, fluids, schedule).shoe_timeline.events
        events_on = solver_on.run(well, fluids, schedule).shoe_timeline.events

        self.assertEqual(len(events_off), len(events_on))
        bands_off = _dispersion_bands(events_off)
        bands_on = _dispersion_bands(events_on)
        self.assertTrue(bands_off, "恒速合成井应产生弥散过渡带")
        self.assertEqual(len(bands_off), len(bands_on))
        for (nxt_o, prv_o, s_off, c_off), (nxt_n, prv_n, s_on, c_on) in zip(bands_off, bands_on):
            with self.subTest(interface=f"{nxt_n}<-{prv_n}"):
                self.assertEqual((nxt_o, prv_o), (nxt_n, prv_n))
                # 事件骨架逐位一致（D_eff/U/t_arrival 链完全同源）
                self.assertEqual(c_off, c_on)
                self.assertGreater(s_off, 0.0)
                rel = abs(s_on - s_off) / s_off
                self.assertLessEqual(rel, 1e-12)

    # ------------------------------------------------------------------
    # 3. 变排量方向：减速井收窄 / hu103 变宽（含尾浆 pin）
    # ------------------------------------------------------------------
    def test_variable_rate_directions(self) -> None:
        """(i) 减速合成井：水泥段 2.0 m³/min → 后继替浆段 0.5 m³/min，
        尾浆界面 σ_on < σ_off（收窄方向受控验证）；
        (ii) hu103 真实井：全部界面 σ_on > σ_off（t_travel 分母 71.101 m³
        liner 圆柱体积 vs 到达链 88.55 m³ 双内径体积的链差被接触时间修正
        暴露——规格承认的修正暴露，非 bug）；尾浆界面 pin：
        σ_on/σ_off == sqrt(shoe_lag_volume / (shoe_md·π·r_liner²))，
        t_inject/t_arrival 链路在测试内独立重算做守卫。"""

        # ---------- (i) 减速合成井 ----------
        shoe_md_m = 100.0
        area_m2 = 0.01
        well = _synthetic_well(shoe_md_m=shoe_md_m, area_m2=area_m2)
        fluids = _synthetic_fluids(spacer_rho=None)
        schedule = PumpingSchedule(
            steps=(
                PumpingScheduleStep(
                    step_name="注入尾浆", fluid_name="尾浆", volume_m3=0.9,
                    rate_m3_min=2.0, start_time_s=0.0, end_time_s=27.0,
                    event_tag=PumpingStageEvent.INJECT_CEMENT,
                ),
                PumpingScheduleStep(
                    step_name="替浆", fluid_name="顶替液", volume_m3=2.0,
                    rate_m3_min=0.5, start_time_s=27.0, end_time_s=267.0,
                    event_tag=PumpingStageEvent.INJECT_DISPLACEMENT,
                ),
            )
        )
        solver_off = CasingFlowSolver(enable_gravity=False)
        solver_on = CasingFlowSolver(enable_gravity=False, mixing_contact_time=True)
        bands_off = _dispersion_bands(solver_off.run(well, fluids, schedule).shoe_timeline.events)
        bands_on = _dispersion_bands(solver_on.run(well, fluids, schedule).shoe_timeline.events)
        self.assertEqual(len(bands_off), len(bands_on))
        tail_off = next(b for b in bands_off if b[0] == "尾浆")
        tail_on = next(b for b in bands_on if b[0] == "尾浆")
        self.assertLess(tail_on[2], tail_off[2], "减速井尾浆界面 σ_on 应小于 σ_off")

        # ---------- (ii) hu103 真实井 ----------
        well_hu, fluids_hu, schedule_hu, _ = _hu103_case()
        solver_hu_off = CasingFlowSolver(enable_gravity=True)
        solver_hu_on = CasingFlowSolver(enable_gravity=True, mixing_contact_time=True)
        r_off = solver_hu_off.run(well_hu, fluids_hu, schedule_hu)
        r_on = solver_hu_on.run(well_hu, fluids_hu, schedule_hu)

        bands_hu_off = _dispersion_bands(r_off.shoe_timeline.events)
        bands_hu_on = _dispersion_bands(r_on.shoe_timeline.events)
        self.assertEqual(len(bands_hu_off), len(bands_hu_on))
        self.assertEqual(len(bands_hu_off), 7, "hu103 应有 7 个界面过渡带")
        for (nxt_o, prv_o, s_off, _), (nxt_n, prv_n, s_on, _) in zip(bands_hu_off, bands_hu_on):
            self.assertEqual((nxt_o, prv_o), (nxt_n, prv_n))
            with self.subTest(interface=f"{nxt_n}<-{prv_n}"):
                self.assertGreater(s_on, s_off, "hu103 全界面 σ_on 应大于 σ_off（双内径链差暴露）")

        # 尾浆界面 pin：σ_on/σ_off == sqrt(shoe_lag_volume / (shoe_md·π·r_liner²))
        tail_band_off = next(b for b in bands_hu_off if b[0] == "尾浆" and b[1] == "中间浆")
        tail_band_on = next(b for b in bands_hu_on if b[0] == "尾浆" and b[1] == "中间浆")
        shoe_lag_volume_m3 = well_hu.shoe_lag_volume_m3
        self.assertIsNotNone(shoe_lag_volume_m3)
        r_liner_m = well_hu.liner_id_mm / 2000.0
        v_liner_m3 = well_hu.shoe_md_m * math.pi * r_liner_m ** 2
        expected_ratio = math.sqrt(shoe_lag_volume_m3 / v_liner_m3)  # sqrt(88.55/71.101)≈1.115980
        actual_ratio = tail_band_on[2] / tail_band_off[2]
        self.assertLessEqual(
            abs(actual_ratio - expected_ratio) / expected_ratio, 1.0e-2,
            f"尾浆 σ 比值 {actual_ratio} 应等于 sqrt(体积链比) {expected_ratio}",
        )

        # 守卫：t_inject/t_arrival 链路在测试内独立重算，防实现悄悄改口径。
        # 到达链：截断序列上 _front_arrival_time + 重力修正（对齐 _build_shoe_timeline）；
        # 注入链：尾浆步出发体积坐标经 _front_arrival_time 同款反解。
        steps_full = solver_hu_on._build_scheduled_steps(schedule_hu)
        steps, _ = solver_hu_on._displacement_sequence_cutoff(steps_full)
        initial_fluid = solver_hu_on._initial_fluid_name(fluids_hu, schedule_hu)
        area_hu = solver_hu_on._pipe_cross_section_area(well_hu)
        pipe_vol = solver_hu_on._timeline_pipe_volume(well_hu, well_hu.shoe_md_m * area_hu)
        tail_idx = next(i for i, s in enumerate(steps) if s.step.fluid_name == "尾浆")
        t_arrival = solver_hu_on._front_arrival_time(steps[tail_idx], steps, pipe_vol)
        self.assertIsNotNone(t_arrival)
        t_arrival = solver_hu_on._gravity_corrected_arrival_time(
            t_arrival, "尾浆",
            solver_hu_on._displaced_fluid_name(steps, tail_idx, initial_fluid),
            fluids_hu, well_hu,
        )
        t_inject = solver_hu_on._inject_start_time(steps[tail_idx], steps)
        band_q_m3_s = _band_flow_rate(r_off.shoe_timeline.events, tail_band_off)
        t_travel = well_hu.shoe_md_m / (
            band_q_m3_s / (math.pi * r_liner_m ** 2)
        )
        # σ 比值独立重算：sqrt(t_contact / t_travel)
        d_off = solver_hu_on._compute_dispersion_coefficient(
            r_liner_m,
            next(f for f in fluids_hu if f.name == "尾浆"),
            band_q_m3_s / (math.pi * r_liner_m ** 2),
        )
        u_rate = _band_flow_rate(r_off.shoe_timeline.events, tail_band_off) / (math.pi * r_liner_m ** 2)
        sigma_off_recomputed = math.sqrt(2.0 * d_off * t_travel) / u_rate
        sigma_on_recomputed = solver_hu_on._contact_time_integrated_sigma(
            t_arrival, t_inject, t_travel, d_off, u_rate, solver_hu_on.dt
        )
        # 两条重算链与实测 band 比值一致（防实现漂移）
        recomputed_ratio = sigma_on_recomputed / sigma_off_recomputed
        self.assertLessEqual(
            abs(recomputed_ratio - actual_ratio) / actual_ratio, 1.0e-9,
            f"独立重算比值 {recomputed_ratio} 应与实测 band 比值 {actual_ratio} 一致",
        )

    # ------------------------------------------------------------------
    # 4. 分段注入归着
    # ------------------------------------------------------------------
    def test_segmented_injection_boundary(self) -> None:
        """同流体分两段注入（中间隔另一流体）时，第二段前缘界面的 t_inject
        取**第二段**出发体积对应时刻（_front_arrival_time 同款反解），
        不是第一段的；与 prev_fluid 跳过同名事件逻辑对齐，不串段。"""

        shoe_md_m = 100.0
        area_m2 = 0.01
        well = _synthetic_well(shoe_md_m=shoe_md_m, area_m2=area_m2)
        fluids = _synthetic_fluids(spacer_rho=None)
        # 尾浆第一段 0→30s（0.5m³ @1.0），隔离液 30→60s（0.5m³ @1.0），
        # 尾浆第二段 60→120s（1.0m³ @1.0），顶替液 120→240s。
        # 管容 1.0m³：尾浆第一段前缘在 60s 到达（0.5+管容1.0=1.5m³ 坐标在隔离液段内）。
        schedule = PumpingSchedule(
            steps=(
                PumpingScheduleStep(
                    step_name="尾浆一段", fluid_name="尾浆", volume_m3=0.5,
                    rate_m3_min=1.0, start_time_s=0.0, end_time_s=30.0,
                    event_tag=PumpingStageEvent.INJECT_CEMENT,
                ),
                PumpingScheduleStep(
                    step_name="注入隔离液", fluid_name="隔离液", volume_m3=0.5,
                    rate_m3_min=1.0, start_time_s=30.0, end_time_s=60.0,
                    event_tag=PumpingStageEvent.INJECT_SPACER,
                ),
                PumpingScheduleStep(
                    step_name="尾浆二段", fluid_name="尾浆", volume_m3=1.0,
                    rate_m3_min=1.0, start_time_s=60.0, end_time_s=120.0,
                    event_tag=PumpingStageEvent.INJECT_CEMENT,
                ),
                PumpingScheduleStep(
                    step_name="替浆", fluid_name="顶替液", volume_m3=2.0,
                    rate_m3_min=1.0, start_time_s=120.0, end_time_s=240.0,
                    event_tag=PumpingStageEvent.INJECT_DISPLACEMENT,
                ),
            )
        )

        solver = CasingFlowSolver(enable_gravity=False, mixing_contact_time=True)
        steps_full = solver._build_scheduled_steps(schedule)
        steps, _ = solver._displacement_sequence_cutoff(steps_full)
        # 第二段（index 2）前缘的归着注入步必须是 steps[2] 本身
        tail_steps = [i for i, s in enumerate(steps) if s.step.fluid_name == "尾浆"]
        self.assertEqual(tail_steps, [0, 2])
        t_inject_seg2 = solver._inject_start_time(steps[2], steps)
        # 第二段出发体积 = 1.0m³ 坐标 → 时刻 60.0s（第二段 start_time_s），
        # 而非第一段出发时刻 0.0s
        self.assertEqual(t_inject_seg2, 60.0)
        self.assertNotEqual(t_inject_seg2, 0.0)
        t_inject_seg1 = solver._inject_start_time(steps[0], steps)
        self.assertEqual(t_inject_seg1, 0.0)

        # 端到端：验证两段前缘的到达时刻与接触时间取值链。
        # 第一段前缘到达：0.0m³ 出发坐标 + 管容 1.0 → 1.0m³ 坐标 → 60s（隔离液段末）；
        # 第二段前缘到达：1.0m³ 出发坐标 + 管容 1.0 → 2.0m³ 坐标 → 120s（尾浆二段末）。
        t_arr1 = solver._front_arrival_time(steps[0], steps, 1.0)
        self.assertIsNotNone(t_arr1)
        t_arr2 = solver._front_arrival_time(steps[2], steps, 1.0)
        self.assertIsNotNone(t_arr2)
        self.assertEqual(t_arr1, 60.0)
        self.assertEqual(t_arr2, 120.0)
        # 接触时间：段1 t_inject=0 → t_contact=60s；段2 t_inject=60（第二段出发）→
        # t_contact=60s。若归着串段（段2 误取第一段 t_inject=0），t_contact 会是 120s
        # ——σ 比值 sqrt(2) 的差异可区分；正确归着下两段 σ 相等。
        d_eff = solver._compute_dispersion_coefficient(
            _liner_id_for_area(area_m2) / 2000.0,
            next(f for f in fluids if f.name == "尾浆"),
            (1.0 / 60.0) / (math.pi * (_liner_id_for_area(area_m2) / 2000.0) ** 2),
        )
        u_rate = (1.0 / 60.0) / (math.pi * (_liner_id_for_area(area_m2) / 2000.0) ** 2)
        t_travel = shoe_md_m / u_rate
        sigma_seg1 = solver._contact_time_integrated_sigma(t_arr1, t_inject_seg1, t_travel, d_eff, u_rate, solver.dt)
        sigma_seg2 = solver._contact_time_integrated_sigma(t_arr2, t_inject_seg2, t_travel, d_eff, u_rate, solver.dt)
        self.assertEqual(sigma_seg2 / sigma_seg1, 1.0)
        # 反证：若串到第一段出发时刻（t_inject=0），段2 σ 会是 sqrt(2) 倍
        sigma_seg2_wrong = solver._contact_time_integrated_sigma(t_arr2, t_inject_seg1, t_travel, d_eff, u_rate, solver.dt)
        self.assertAlmostEqual(sigma_seg2_wrong / sigma_seg2, math.sqrt(2.0), places=9)


class TestPlugFaceZeroMixing(unittest.TestCase):
    """胶塞面零掺混开关（plug_face_zero_mixing，路线 B Task 2）行为测试。"""

    # ------------------------------------------------------------------
    # 5. 开关开：尾浆→压塞液界面零过渡，其余界面不受影响
    # ------------------------------------------------------------------
    def test_plug_face_no_dispersion(self) -> None:
        """ht1_003（碰压成功井；loader 生产口径 has_plug=False——单测直接
        属性覆盖 solver.has_plug=True + plug_face_zero_mixing=True）：
        "尾浆→压塞液"界面（prev=尾浆, next=压塞液）无过渡子事件——
        弥散函数保留原单相阶跃 FRONT_ARRIVAL（("压塞液",1.0)），不再生成
        双相 (next,prev) 事件对；其余 5 个界面对（钻井液→平衡液/平衡液→
        隔离液1/隔离液1→隔离液2/隔离液2→领浆/领浆→尾浆）的过渡带子事件数
        与开关关时逐位相同（双相事件数不变、band 中心与 σ 不变）。
        """

        well, fluids, schedule, _ = _ht1_003_case()

        solver_off = CasingFlowSolver(enable_gravity=True)
        solver_off.has_plug = True  # 属性覆盖（不重跑 loader），混浆增强=1 口径
        result_off = solver_off.run(well, fluids, schedule)
        events_off = result_off.shoe_timeline.events

        solver_on = CasingFlowSolver(enable_gravity=True, plug_face_zero_mixing=True)
        solver_on.has_plug = True
        result_on = solver_on.run(well, fluids, schedule)
        events_on = result_on.shoe_timeline.events

        # 压塞液界面在开关关时确有过渡带（探针实证 5 子事件）——守卫前置
        bands_off = _dispersion_bands(events_off)
        self.assertIn(("压塞液", "尾浆"), {(n, p) for n, p, _, _ in bands_off})

        # 开关开：无任何 (压塞液, *) 双相子事件
        for e in events_on:
            if e.kind == ShoeEventKind.FRONT_ARRIVAL and len(e.phase_fractions) == 2:
                self.assertNotEqual(
                    e.phase_fractions[0][0], "压塞液",
                    "开关开后不应存在 next=压塞液 的双相过渡子事件",
                )

        # 原单相阶跃事件保留：("压塞液", 1.0) 的 FRONT_ARRIVAL 恰好 1 个
        single_plug_on = [
            e for e in events_on
            if e.kind == ShoeEventKind.FRONT_ARRIVAL
            and e.phase_fractions == (("压塞液", 1.0),)
        ]
        self.assertEqual(len(single_plug_on), 1)
        single_plug_off = [
            e for e in events_off
            if e.kind == ShoeEventKind.FRONT_ARRIVAL
            and e.phase_fractions == (("压塞液", 1.0),)
        ]
        # 开关关时该时刻只存在过渡子事件（无单相阶跃）——探针实证
        self.assertEqual(len(single_plug_off), 0)
        # 保留的阶跃事件与开关关时压塞液界面的到达时刻（band 中心）同刻：
        # 阶跃事件本身即原 FRONT_ARRIVAL，其 time_s 与过渡带中心一致
        plug_band_off = next(b for b in bands_off if b[0] == "压塞液" and b[1] == "尾浆")
        self.assertEqual(single_plug_on[0].time_s, plug_band_off[3])
        self.assertEqual(single_plug_on[0].flow_rate_m3_s,
                         _band_flow_rate(events_off, plug_band_off))

        # 其余 5 个界面对：过渡带逐位不受影响（双相事件数、band 中心与 σ）。
        # _dispersion_bands 键约定为 (next, prev)——与探针实证的
        # ('压塞液','尾浆') 同向，注意别写成 (prev, next)。
        other_pairs = {
            ("平衡液", "钻井液"), ("隔离液1", "平衡液"), ("隔离液2", "隔离液1"),
            ("领浆", "隔离液2"), ("尾浆", "领浆"),
        }
        bands_on = _dispersion_bands(events_on)
        pairs_off = {(n, p): (s, c) for n, p, s, c in bands_off}
        pairs_on = {(n, p): (s, c) for n, p, s, c in bands_on}
        self.assertEqual(set(pairs_on), other_pairs, "开关开后应恰余 5 个界面过渡带")
        for pair in other_pairs:
            with self.subTest(interface=f"{pair[0]}<-{pair[1]}"):
                self.assertIn(pair, pairs_off)
                self.assertEqual(pairs_on[pair], pairs_off[pair])

        # 体积链不受本开关影响（仅时间线后处理）
        self.assertEqual(result_off.cement_end_time_s, result_on.cement_end_time_s)

    # ------------------------------------------------------------------
    # 6. 开关关：ht1_003 生产口径逐位回归（fresh 参照，不依赖冻结文件）
    # ------------------------------------------------------------------
    def test_plug_face_off_bitwise(self) -> None:
        """ht1_003（has_plug 默认 False，生产口径）：仅 plug_face_zero_mixing=True
        vs 开关全关 fresh 参照，shoe_timeline 逐事件逐位相同（time_s/kind/
        flow_rate/stage/phase_fractions == 级对比）。

        生产口径 has_plug=False 使判据第二项（has_plug 门）恒假 → 开关空转、
        零效应；这是"开关关逐位回归"的等价强断言（等价于 Task 1 的
        test_default_off_bitwise 锚思路，但不重复建 ht1_003 冻结文件）。
        """

        well, fluids, schedule, _ = _ht1_003_case()

        solver_ref = CasingFlowSolver(enable_gravity=True)
        self.assertFalse(solver_ref.has_plug)  # loader 生产口径，默认 False
        self.assertFalse(solver_ref.plug_face_zero_mixing)  # 默认关
        events_ref = solver_ref.run(well, fluids, schedule).shoe_timeline.events

        solver_on = CasingFlowSolver(enable_gravity=True, plug_face_zero_mixing=True)
        self.assertFalse(solver_on.has_plug)
        events_on = solver_on.run(well, fluids, schedule).shoe_timeline.events

        self.assertEqual(len(events_ref), len(events_on))
        for i, (ev, ref) in enumerate(zip(events_on, events_ref)):
            with self.subTest(event_index=i):
                self.assertEqual(ev.time_s, ref.time_s)
                self.assertEqual(ev.kind, ref.kind)
                self.assertEqual(ev.flow_rate_m3_s, ref.flow_rate_m3_s)
                self.assertEqual(ev.stage_name, ref.stage_name)
                self.assertEqual(ev.phase_fractions, ref.phase_fractions)

    # ------------------------------------------------------------------
    # 7. hu102 保护性断言：界面集无压塞液界面，开关零效应
    # ------------------------------------------------------------------
    def test_hu102_unaffected(self) -> None:
        """hu102（生产 has_plug=True 语义井、替浆步被 RESTART 截断）+
        plug_face_zero_mixing=True：shoe_timeline 与开关关时逐位相同。

        保护性断言依据：hu102 顶替序列被 RESTART（循环排混浆）截断，压塞液/
        后置液/替浆液前缘均未到达鞋口，其界面集无"尾浆→压塞液"异物对
        （探针实证：过渡带对仅 钻井液→平衡液/隔离液→平衡液/领浆→隔离液/
        尾管水泥浆→领浆）——开关不应有任何效应。has_plug=True 直传构造参数
        （混浆增强=1 语义口径）。"""

        well, fluids, schedule, _ = _hu102_case()

        solver_off = CasingFlowSolver(enable_gravity=True, has_plug=True)
        r_off = solver_off.run(well, fluids, schedule)
        solver_on = CasingFlowSolver(enable_gravity=True, has_plug=True,
                                     plug_face_zero_mixing=True)
        r_on = solver_on.run(well, fluids, schedule)

        # 守卫：确认界面集确无压塞液界面（若未来 hu102 程序变更使压塞液
        # 前缘到达鞋口，本测试的"零效应"前提失效，须重裁定断言）
        events_off = r_off.shoe_timeline.events
        plug_next_events = [
            e for e in events_off
            if e.kind == ShoeEventKind.FRONT_ARRIVAL
            and e.phase_fractions and e.phase_fractions[0][0] == "压塞液"
        ]
        self.assertEqual(plug_next_events, [])

        # 开关开 vs 关：逐事件逐位相同
        events_on = r_on.shoe_timeline.events
        self.assertEqual(len(events_off), len(events_on))
        for i, (ev, ref) in enumerate(zip(events_on, events_off)):
            with self.subTest(event_index=i):
                self.assertEqual(ev.time_s, ref.time_s)
                self.assertEqual(ev.kind, ref.kind)
                self.assertEqual(ev.flow_rate_m3_s, ref.flow_rate_m3_s)
                self.assertEqual(ev.stage_name, ref.stage_name)
                self.assertEqual(ev.phase_fractions, ref.phase_fractions)
        self.assertEqual(r_off.cement_end_time_s, r_on.cement_end_time_s)


def _band_flow_rate(events, band) -> float:
    """取过渡带首事件的排量（m³/s）。band = (next, prev, sigma, center)。"""

    center = band[3]
    for e in events:
        if (e.kind == ShoeEventKind.FRONT_ARRIVAL and len(e.phase_fractions) == 2
                and e.phase_fractions[0][0] == band[0] and e.phase_fractions[1][0] == band[1]
                and abs(e.time_s - center) < 1.0e-9):
            return e.flow_rate_m3_s
    raise AssertionError(f"band {band[:2]} 未找到对应事件")


if __name__ == "__main__":
    _ = unittest.main()
