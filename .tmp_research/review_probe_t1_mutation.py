# -*- coding: utf-8 -*-
"""审查探针 2：变异测试验证恒速锚/减速井测试区分度 + dt 调小后的数学等价。

M1: 恒速锚在 dt=2.0（现状）下，把实现的时间项"故意错化"（t_contact→t_arrival），
    看 rel 断言是否仍通过（区分度为零的证据）。
M2: 恒速锚在 dt=0.05 下，正确实现 vs 错化实现的 rel（区分度恢复的证据）。
M3: 减速井 σ_off/σ_on raw 值与 clamp 状态复核。
"""
import math, sys
sys.path.insert(0, ".")

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep, PumpingStageEvent
from cemdisp.data.well_spec import WellSpec
from cemdisp.transport1d.casing_flow import CasingFlowSolver
from cemdisp.transport1d.shoe_timeline import ShoeEventKind


def _liner_id_for_area(area_m2):
    return math.sqrt(4.0 * area_m2 / math.pi) * 1000.0


def _fluids():
    return (
        FluidSpec(name="泥浆", role=FluidRole.MUD, density_kg_m3=1200.0, plastic_viscosity_pa_s=0.02),
        FluidSpec(name="隔离液", role=FluidRole.SPACER, density_kg_m3=1200.0, plastic_viscosity_pa_s=0.01),
        FluidSpec(name="尾浆", role=FluidRole.TAIL, density_kg_m3=1900.0, plastic_viscosity_pa_s=0.08),
        FluidSpec(name="顶替液", role=FluidRole.DISPLACEMENT, density_kg_m3=1050.0, plastic_viscosity_pa_s=0.01),
    )


def bands(events):
    by_pair = {}
    for e in events:
        if e.kind != ShoeEventKind.FRONT_ARRIVAL or len(e.phase_fractions) != 2:
            continue
        if e.phase_fractions[0][1] >= 1.0:
            continue
        by_pair.setdefault((e.phase_fractions[0][0], e.phase_fractions[1][0]), []).append(e.time_s)
    out = []
    for (nxt, prv), ts in by_pair.items():
        ts.sort()
        out.append((nxt, prv, (ts[-1] - ts[0]) / 2.0, (ts[-1] + ts[0]) / 2.0))
    return out


CONST_SCHEDULE = PumpingSchedule(steps=(
    PumpingScheduleStep(step_name="注入隔离液", fluid_name="隔离液", volume_m3=0.4,
                        rate_m3_min=1.0, start_time_s=0.0, end_time_s=24.0,
                        event_tag=PumpingStageEvent.INJECT_SPACER),
    PumpingScheduleStep(step_name="注入尾浆", fluid_name="尾浆", volume_m3=1.0,
                        rate_m3_min=1.0, start_time_s=24.0, end_time_s=84.0,
                        event_tag=PumpingStageEvent.INJECT_CEMENT),
    PumpingScheduleStep(step_name="替浆", fluid_name="顶替液", volume_m3=2.0,
                        rate_m3_min=1.0, start_time_s=84.0, end_time_s=204.0,
                        event_tag=PumpingStageEvent.INJECT_DISPLACEMENT),
))


def run_const(dt):
    well = WellSpec(well_name="合成井", top_md_m=1.0, bottom_md_m=120.0,
                    shoe_md_m=100.0, liner_id_mm=_liner_id_for_area(0.01))
    off = CasingFlowSolver(enable_gravity=False, dt=dt).run(well, _fluids(), CONST_SCHEDULE)
    on = CasingFlowSolver(enable_gravity=False, dt=dt, mixing_contact_time=True).run(well, _fluids(), CONST_SCHEDULE)
    return bands(off.shoe_timeline.events), bands(on.shoe_timeline.events)


def M_const_anchor():
    print("=== M1/M2: 恒速锚区分度 ===")
    for dt in (2.0, 0.05):
        bo, bn = run_const(dt)
        print(f"[dt={dt}]")
        for (no, po, so, co), (nn, pn, sn, cn) in zip(bo, bn):
            raw_est = None
            rel = abs(sn - so) / so
            print(f"  {nn}<-{po}: sigma_off={so:.6f} sigma_on={sn:.6f} rel={rel:.3e} "
                  f"{'CLAMPED(both<=dt)' if so <= dt and sn <= dt else ''}")

    # 错化变异：σ_wrong ∝ sqrt(t_arrival)（尾浆 t_arr=84, 隔离液 t_arr=60）
    R = _liner_id_for_area(0.01) / 2000.0
    U = (1.0 / 60.0) / 0.01
    D = 0.25 * U * R  # 对流 cap 主导
    dt = 2.0
    s_off = max(math.sqrt(2 * D * 60.0) / U, dt)   # t_travel=60
    s_correct = max(math.sqrt(2 * D * 60.0) / U, dt)  # t_contact=60（恒速=等价）
    s_wrong_arr = max(math.sqrt(2 * D * 84.0) / U, dt)  # 错化: t_contact→t_arrival=84
    s_wrong_inj = max(math.sqrt(2 * D * 24.0) / U, dt)  # 错化: t_contact→t_inject=24
    print(f"[dt=2.0 变异分析] σ_off={s_off:.4f} σ_correct={s_correct:.4f} "
          f"σ_wrong(t_arrival)={s_wrong_arr:.4f} σ_wrong(t_inject)={s_wrong_inj:.4f}")
    print(f"  → 错化实现 rel vs off: {abs(s_wrong_arr-s_off)/s_off:.3e} / {abs(s_wrong_inj-s_off)/s_off:.3e}"
          f"  （均 ≤1e-12? {'是→锚无区分度' if max(abs(s_wrong_arr-s_off),abs(s_wrong_inj-s_off))/s_off<=1e-12 else '否→有区分度'}）")
    dt = 0.05
    s_off5 = max(math.sqrt(2 * D * 60.0) / U, dt)
    s_wrong5 = math.sqrt(2 * D * 84.0) / U
    print(f"[dt=0.05 变异分析] σ_off={s_off5:.4f} σ_wrong(t_arrival)={s_wrong5:.4f} "
          f"rel={abs(s_wrong5-s_off5)/s_off5:.3e}（应≈0.183→区分度恢复）")


def M_decel():
    print("=== M3: 减速井 clamp 状态 ===")
    well = WellSpec(well_name="合成井", top_md_m=1.0, bottom_md_m=120.0,
                    shoe_md_m=100.0, liner_id_mm=_liner_id_for_area(0.01))
    fluids = _fluids()
    sched = PumpingSchedule(steps=(
        PumpingScheduleStep(step_name="注入尾浆", fluid_name="尾浆", volume_m3=0.9,
                            rate_m3_min=2.0, start_time_s=0.0, end_time_s=27.0,
                            event_tag=PumpingStageEvent.INJECT_CEMENT),
        PumpingScheduleStep(step_name="替浆", fluid_name="顶替液", volume_m3=2.0,
                            rate_m3_min=0.5, start_time_s=27.0, end_time_s=267.0,
                            event_tag=PumpingStageEvent.INJECT_DISPLACEMENT),
    ))
    off = CasingFlowSolver(enable_gravity=False).run(well, fluids, sched)
    on = CasingFlowSolver(enable_gravity=False, mixing_contact_time=True).run(well, fluids, sched)
    bo, bn = bands(off.shoe_timeline.events), bands(on.shoe_timeline.events)
    t_off = next(b for b in bo if b[0] == "尾浆")
    t_on = next(b for b in bn if b[0] == "尾浆")
    print(f"尾浆 σ_off={t_off[2]:.6f} (dt=2.0, {'CLAMPED' if t_off[2] <= 2.0+1e-12 else 'raw'}) "
          f"σ_on={t_on[2]:.6f} ({'CLAMPED=dt' if abs(t_on[2]-2.0)<1e-9 else 'raw'}) "
          f"通过? {t_on[2] < t_off[2]}")
    # raw 值重算
    Q_tail = t_off[3] and next(e.flow_rate_m3_s for e in off.shoe_timeline.events
                               if e.kind == ShoeEventKind.FRONT_ARRIVAL and len(e.phase_fractions) == 2
                               and e.phase_fractions[0][0] == "尾浆" and abs(e.time_s - t_off[3]) < 1e-9)
    U2 = Q_tail / (math.pi * R ** 2)
    D2 = 0.25 * U2 * R
    solver = CasingFlowSolver(enable_gravity=False)
    steps_full = solver._build_scheduled_steps(sched)
    steps, _ = solver._displacement_sequence_cutoff(steps_full)
    vol = solver._timeline_pipe_volume(well, 100.0 * math.pi * R ** 2)
    t_arr = solver._front_arrival_time(steps[0], steps, vol)
    t_inj = solver._inject_start_time(steps[0], steps)
    print(f"t_arr={t_arr:.4f} t_inj={t_inj:.4f} t_contact={t_arr-t_inj:.4f} t_travel={100.0/U2:.4f}")
    print(f"raw σ_off={math.sqrt(2*D2*(100.0/U2))/U2:.4f} raw σ_on={math.sqrt(2*D2*(t_arr-t_inj))/U2:.4f} (dt=2.0)")


if __name__ == "__main__":
    M_const_anchor()
    print()
    M_decel()
