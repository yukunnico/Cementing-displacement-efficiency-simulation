# -*- coding: utf-8 -*-
"""审查探针：验证 Task 1 测试的区分度与实现口径。

探针1：恒速锚测试的 σ 是否被 max(σ, dt) 下限吞掉（平凡通过风险）。
探针2：hu103 尾浆 pin 的 sqrt 比值数学——为何无重力修正的期望值能匹配实测？
探针3：归着游标在"事件数 > 注入步数"（异常日程）下的防御回退路径。
"""
import sys, math
sys.path.insert(0, ".")

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.pumping_schedule import (
    PumpingSchedule, PumpingScheduleStep, PumpingStageEvent,
)
from cemdisp.data.well_spec import WellSpec
from cemdisp.transport1d.casing_flow import CasingFlowSolver
from cemdisp.transport1d.shoe_timeline import ShoeEventKind


def _liner_id_for_area(area_m2: float) -> float:
    return math.sqrt(4.0 * area_m2 / math.pi) * 1000.0


def _synthetic_fluids():
    return (
        FluidSpec(name="泥浆", role=FluidRole.MUD, density_kg_m3=1200.0,
                  plastic_viscosity_pa_s=0.02),
        FluidSpec(name="隔离液", role=FluidRole.SPACER, density_kg_m3=1200.0,
                  plastic_viscosity_pa_s=0.01),
        FluidSpec(name="尾浆", role=FluidRole.TAIL, density_kg_m3=1900.0,
                  plastic_viscosity_pa_s=0.08),
        FluidSpec(name="顶替液", role=FluidRole.DISPLACEMENT, density_kg_m3=1050.0,
                  plastic_viscosity_pa_s=0.01),
    )


def probe1_constant_rate_anchor_discrimination():
    """恒速锚测试：复现测试日程，检查 σ_off 与 dt 下限的相对大小。"""
    shoe_md_m, area_m2 = 100.0, 0.01
    well = WellSpec(
        well_name="合成井", top_md_m=1.0, bottom_md_m=shoe_md_m + 20.0,
        shoe_md_m=shoe_md_m, liner_id_mm=_liner_id_for_area(area_m2),
    )
    fluids = _synthetic_fluids()
    schedule = PumpingSchedule(steps=(
        PumpingScheduleStep(step_name="注入隔离液", fluid_name="隔离液",
                            volume_m3=0.4, rate_m3_min=1.0, start_time_s=0.0,
                            end_time_s=24.0, event_tag=PumpingStageEvent.INJECT_SPACER),
        PumpingScheduleStep(step_name="注入尾浆", fluid_name="尾浆",
                            volume_m3=1.0, rate_m3_min=1.0, start_time_s=24.0,
                            end_time_s=84.0, event_tag=PumpingStageEvent.INJECT_CEMENT),
        PumpingScheduleStep(step_name="替浆", fluid_name="顶替液",
                            volume_m3=2.0, rate_m3_min=1.0, start_time_s=84.0,
                            end_time_s=204.0, event_tag=PumpingStageEvent.INJECT_DISPLACEMENT),
    ))
    solver = CasingFlowSolver(enable_gravity=False)
    # U = Q/A；Q = 1/60 m³/s, A = 0.01 m² → U = 1.6667 m/s
    U = (1.0 / 60.0) / area_m2
    R = _liner_id_for_area(area_m2) / 2000.0
    fluid_tail = next(f for f in fluids if f.name == "尾浆")
    D_eff = solver._compute_dispersion_coefficient(R, fluid_tail, U)
    t_travel = shoe_md_m / U  # 60 s
    sigma_raw = math.sqrt(2.0 * D_eff * t_travel) / U
    print(f"[probe1] U={U:.4f} m/s  R={R:.5f} m  D_eff={D_eff:.4e} m2/s")
    print(f"[probe1] t_travel={t_travel:.1f} s  sigma_raw={sigma_raw:.6f} s  dt=2.0")
    print(f"[probe1] sigma_raw/dt = {sigma_raw/2.0:.3f}  → {'CLAMPED by max(sigma,dt)!' if sigma_raw < 2.0 else 'not clamped, anchor has discrimination'}")
    # 也算隔离液界面（前一界面，prev=泥浆? 首界面被跳过？）
    fluid_spacer = next(f for f in fluids if f.name == "隔离液")
    D_sp = solver._compute_dispersion_coefficient(R, fluid_spacer, U)
    sigma_sp = math.sqrt(2.0 * D_sp * t_travel) / U
    print(f"[probe1] spacer D_eff={D_sp:.4e}  sigma_sp={sigma_sp:.6f}  clamped={sigma_sp < 2.0}")


def probe2_hu103_tail_pin_math():
    """hu103 尾浆 pin：检查 σ_on/σ_off == sqrt(V_lag/V_liner) 的数学成立条件。

    关键疑问：t_arrival 与 t_inject 都含重力修正吗？
    t_inject 来自 _inject_start_time（无重力修正，纯体积反解）。
    t_arrival 来自 _front_arrival_time + 重力修正。
    t_contact = t_arrival_gravity - t_inject。
    旧 σ_off ∝ sqrt(t_travel) = sqrt(shoe_md/U)。
    新 σ_on ∝ sqrt(t_contact)。
    pin 声称 σ_on/σ_off == sqrt(V_lag/V_liner) = sqrt(t_contact/t_travel)
    → 即 t_contact/t_travel == V_lag/V_liner == 1.2454？
    但 t_contact = t_arrival_grav - t_inject ≠ V_lag/U 除非 U 恒定且无重力。
    hu103 变排量+重力 → pin 精确匹配到 2e-13 很可疑，需查明。
    """
    from cemdisp.data.loaders import load_hu103_tailpipe
    well, fluids, schedule, _ = load_hu103_tailpipe()
    solver_on = CasingFlowSolver(enable_gravity=True, mixing_contact_time=True)
    solver_off = CasingFlowSolver(enable_gravity=True)
    r_off = solver_off.run(well, fluids, schedule)
    r_on = solver_on.run(well, fluids, schedule)

    def bands(events):
        by_pair = {}
        for e in events:
            if e.kind != ShoeEventKind.FRONT_ARRIVAL or len(e.phase_fractions) != 2:
                continue
            if e.phase_fractions[0][1] >= 1.0:
                continue
            key = (e.phase_fractions[0][0], e.phase_fractions[1][0])
            by_pair.setdefault(key, []).append(e.time_s)
        out = []
        for (nxt, prv), times in by_pair.items():
            times.sort()
            out.append((nxt, prv, (times[-1] - times[0]) / 2.0,
                        (times[-1] + times[0]) / 2.0))
        return out

    bo, bn = bands(r_off.shoe_timeline.events), bands(r_on.shoe_timeline.events)
    tail_off = next(b for b in bo if b[0] == "尾浆" and b[1] == "中间浆")
    tail_on = next(b for b in bn if b[0] == "尾浆" and b[1] == "中间浆")
    actual = tail_on[2] / tail_off[2]

    steps_full = solver_on._build_scheduled_steps(schedule)
    steps, _ = solver_on._displacement_sequence_cutoff(steps_full)
    initial_fluid = solver_on._initial_fluid_name(fluids, schedule)
    area = solver_on._pipe_cross_section_area(well)
    pipe_vol = solver_on._timeline_pipe_volume(well, well.shoe_md_m * area)
    tail_idx = next(i for i, s in enumerate(steps) if s.step.fluid_name == "尾浆")
    t_arr_pure = solver_on._front_arrival_time(steps[tail_idx], steps, pipe_vol)
    t_arr_grav = solver_on._gravity_corrected_arrival_time(
        t_arr_pure, "尾浆",
        solver_on._displaced_fluid_name(steps, tail_idx, initial_fluid),
        fluids, well)
    t_inject = solver_on._inject_start_time(steps[tail_idx], steps)
    t_contact = t_arr_grav - t_inject

    Q = tail_off[3] and None  # band center 时刻的排量
    # 找 band center 时刻的事件排量（测试用 _band_flow_rate）
    for e in r_off.shoe_timeline.events:
        if (e.kind == ShoeEventKind.FRONT_ARRIVAL and len(e.phase_fractions) == 2
                and e.phase_fractions[0][0] == "尾浆" and e.phase_fractions[1][0] == "中间浆"
                and abs(e.time_s - tail_off[3]) < 1e-9):
            Q = e.flow_rate_m3_s
            break
    r_liner = well.liner_id_mm / 2000.0
    U = Q / (math.pi * r_liner ** 2)
    t_travel = well.shoe_md_m / U
    print(f"[probe2] hu103: t_arr_pure={t_arr_pure:.4f}  t_arr_grav={t_arr_grav:.4f}")
    print(f"[probe2] t_inject={t_inject:.4f}  t_contact={t_contact:.4f}")
    print(f"[probe2] Q(center)={Q:.8f} m3/s  U={U:.6f} m/s  t_travel={t_travel:.4f}")
    print(f"[probe2] V_lag={well.shoe_lag_volume_m3:.4f}  V_liner={well.shoe_md_m*math.pi*r_liner**2:.4f}")
    print(f"[probe2] t_contact/t_travel = {t_contact/t_travel:.10f}")
    print(f"[probe2] V_lag/V_liner     = {well.shoe_lag_volume_m3/(well.shoe_md_m*math.pi*r_liner**2):.10f}")
    print(f"[probe2] sqrt both: {math.sqrt(t_contact/t_travel):.10f} vs {math.sqrt(well.shoe_lag_volume_m3/(well.shoe_md_m*math.pi*r_liner**2)):.10f}")
    print(f"[probe2] actual band ratio = {actual:.10f}")
    print(f"[probe2] gravity contributes? grav_time={t_arr_grav - t_arr_pure:.4f} s (t_contact 中的占比 {(t_arr_grav-t_arr_pure)/t_contact*100:.4f}%)")


def probe3_cursor_defense_path():
    """异常日程：事件数 > 该流体注入步数时游标耗尽 → t_contact=0 → σ=dt。

    检查该路径是否真的可达：FRONT_ARRIVAL 事件来自 _build_shoe_timeline 里
    每步一个 front 事件（_front_arrival_time 非 None），且 prev_fluid 异名。
    事件数 ≤ 步数恒成立？何种日程会耗尽？
    """
    # 构造一个游标耗尽场景：同一流体两步，但第一步 front 事件被 prev_fluid
    # 同名跳过逻辑吞掉 —— 实际上 dispense 只在 FRONT_ARRIVAL 分支消耗游标。
    # 每个非 None front 都会消耗。事件数=步数。但若某步 front is None（体积不够），
    # 弥散缓存里仍用 end_time_s 标记 → 缓存长度=步数。
    # 耗尽场景：EVENT 数 > 缓存长度不可能（每事件消耗一次，缓存每步一项）。
    # 除非：fluid_name 相同但事件数 > 步数——不可能，因为每步至多一个 front 事件。
    # 但注意：缓存键是 fluid_name，事件也按 fluid_name 匹配——一致。
    # 结论：正常构造下耗尽不可达；只有"手动构造重复事件"才触发。
    # 验证一下守卫逻辑本身：
    solver = CasingFlowSolver(enable_gravity=False)
    # 静态检查：cursor 越界分支 → t_inject = t_arrival → t_contact = 0 → σ = sqrt(0)=0 → max(σ,dt)=dt
    sigma = solver._contact_time_integrated_sigma(100.0, 100.0, 60.0, 1e-4, 1.5, 2.0)
    print(f"[probe3] t_contact=0 → sigma={sigma} (期望 =dt=2.0)")
    assert sigma == 2.0
    print("[probe3] 游标耗尽防御路径：σ 落到 dt 下限（2s），不会崩，但会静默产出极窄过渡带")


if __name__ == "__main__":
    probe1_constant_rate_anchor_discrimination()
    print()
    probe2_hu103_tail_pin_math()
    print()
    probe3_cursor_defense_path()
