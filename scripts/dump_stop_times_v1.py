"""dump_stop_times_v1: 顶替结束标志取证脚本（1D 时间轴实算，只读不改任何仓库数据）。

对 8 井输出：
- 模型 1D 管容 / 日程总泵量 / 日程总时长
- 各注入步前缘到达鞋口时刻（活塞流 vs 重力修正）
- 末段水泥尾缘到达时刻（修正前后）
- runner 口径 annulus_stop_time_s（复刻 cemdisp/runners/*_tailpipe.py）
- RR 口径 tt = min(日程总时长+1200s, stop+600s)（复刻 scripts/rerun_all_wells_corrected.py:67）
- F2 判定：末段水泥尾缘（修正后）是否晚于 stop（尾缘错配+重排截流判定）+ 截流体积
- 三种评价时刻（stop / stop+600 / 日程总时长+1200）的已出鞋口水泥体积与 eta_E 采样
- 超替检查：替浆步体积 vs 模型管容
- hu101 停泵后浮力交换流特征时间标度锚点（量级估计，非模拟）

用法：PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/dump_stop_times_v1.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import cemdisp.data.loaders as L
from cemdisp.data.fluid_spec import FluidRole
from cemdisp.transport1d import CasingFlowSolver

WELLS = [
    ("hu101", L.load_hu101_tailpipe),
    ("hu102", L.load_hu102_tailpipe),
    ("hu103", L.load_hu103_tailpipe),
    ("hu1", L.load_hu1_tailpipe),
    ("hu2", L.load_hu2_tailpipe),
    ("ht1_001", L.load_ht1_001_tailpipe),
    ("ht1_003", L.load_ht1_003_tailpipe),
    ("ht1_004", L.load_ht1_004_tailpipe),
]

CEMENT_ROLES = {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}


def runner_stop(cr, fluids) -> float:
    """复刻 runner/RR 脚本的 annulus_stop_time_s（含乱序防御排序）。"""
    by = {f.name: f for f in fluids}
    fs = sorted(cr.fronts, key=lambda f: f.time_s)
    last_cem = None
    for f in fs:
        fl = by.get(f.fluid_name)
        if fl is not None and fl.role in CEMENT_ROLES:
            last_cem = f.time_s
    if last_cem is not None:
        for f in fs:
            fl = by.get(f.fluid_name)
            if fl is None or fl.role in CEMENT_ROLES:
                continue
            if f.time_s >= last_cem - 1e-9:
                return float(f.time_s)
    return float(cr.cement_end_time_s)


def _cement_out_volume(cement_steps, cum_pumped, pipe_vol):
    """活塞流口径：累计泵入 cum_pumped 时已出鞋口的水泥体积（m3）。

    水泥步 i 的前缘在泵入 pipe_vol+cum0_i 时出鞋口，尾缘在 pipe_vol+cum1_i 时出完。
    """
    out = 0.0
    for c in cement_steps:
        front_v = pipe_vol + c["cum0"]
        rear_v = pipe_vol + c["cum1"]
        a = max(cum_pumped - front_v, 0.0)
        b = min(cum_pumped, rear_v) - front_v
        out += max(min(a, b), 0.0)
    return out


def main() -> None:
    print("=" * 100)
    for label, loader in WELLS:
        well, fluids, schedule, _ = loader()
        by = {f.name: f for f in fluids}

        solver_g = CasingFlowSolver(enable_gravity=True)
        cr_g = solver_g.run(well, fluids, schedule)
        cr_p = CasingFlowSolver(enable_gravity=False).run(well, fluids, schedule)

        pipe_area = cr_g.pipe_cross_section_m2
        pipe_vol = cr_g.shoe_md_m * pipe_area
        sched = list(schedule.steps)
        step_rows = []
        cum = 0.0
        t0 = 0.0
        pumping_end = 0.0
        total_vol = 0.0
        for s in sched:
            dur = 0.0 if s.rate_m3_min <= 0.0 else s.volume_m3 / s.rate_m3_min * 60.0
            step_rows.append({
                "name": s.step_name,
                "fluid": s.fluid_name,
                "role": by[s.fluid_name].role.name,
                "vol": s.volume_m3,
                "rate": s.rate_m3_min,
                "cum0": cum,
                "cum1": cum + s.volume_m3,
                "t0": t0,
                "t1": t0 + dur,
            })
            cum += s.volume_m3
            total_vol += s.volume_m3
            t0 += dur
            pumping_end = max(pumping_end, t0)

        # 每步前缘到达时刻：活塞流用 cr_p.fronts（同序），重力用 cr_g.fronts
        fronts_p = list(cr_p.fronts)
        fronts_g = list(cr_g.fronts)
        for i, row in enumerate(step_rows):
            row["front_piston"] = fronts_p[i].time_s if i < len(fronts_p) else float("nan")
            row["front_grav"] = fronts_g[i].time_s if i < len(fronts_g) else float("nan")

        # 末段水泥步（最后一个水泥角色步）
        cement_idx = [i for i, r in enumerate(step_rows) if r["role"] in ("LEAD", "INTERMEDIATE", "TAIL")]
        last_cem_i = cement_idx[-1] if cement_idx else None

        # 尾缘（重力修正）：复刻 run() 内逻辑 target = cum1 + pipe_vol
        def rear_arrival_piston(idx: int) -> float | None:
            target = step_rows[idx]["cum1"] + pipe_vol
            for r in step_rows:
                if target <= r["cum1"] + 1e-12:
                    if r["rate"] <= 0.0:
                        return r["t1"]
                    return r["t0"] + max(target - r["cum0"], 0.0) / r["rate"] * 60.0
            return None

        rear_p = rear_arrival_piston(last_cem_i) if last_cem_i is not None else None
        rear_g = cr_g.cement_end_time_s if last_cem_i is not None else None

        stop = runner_stop(cr_g, fluids)
        total_t_sched = pumping_end
        rr_tt = min(total_t_sched + 1200.0, stop + 600.0)

        # 三评价时刻的已出鞋口水泥体积（活塞流口径）——先算，供 F2 截流口径 (b) 使用
        cement_steps = [r for r in step_rows if r["role"] in ("LEAD", "INTERMEDIATE", "TAIL")]
        cement_total = sum(r["vol"] for r in cement_steps)
        ann_vol = 0.0
        # 环空体积：hole_diameter_profile 与 liner_od 沿深度积分，截面积 = pi/4 (D_hole^2 - D_od^2)
        od = well.liner_od_mm / 1000.0
        pts = well.hole_diameter_profile
        for k in range(len(pts) - 1):
            d0, d1 = pts[k].depth_md_m, pts[k + 1].depth_md_m
            a0 = math.pi / 4.0 * ((pts[k].value / 1000.0) ** 2 - od**2)
            a1 = math.pi / 4.0 * ((pts[k + 1].value / 1000.0) ** 2 - od**2)
            if a0 <= 0 or a1 <= 0:
                continue
            ann_vol += 0.5 * (a0 + a1) * (d1 - d0)
        eta_samples = {}
        for tag, tt in (("stop", stop), ("stop+600", stop + 600.0), ("sched+1200", total_t_sched + 1200.0)):
            # 重力修正时间轴 -> 累计泵量：按泵注步骤线性插值
            v = total_vol
            for r in step_rows:
                if tt <= r["t1"] + 1e-9:
                    if r["t1"] > r["t0"]:
                        v = r["cum0"] + (tt - r["t0"]) / (r["t1"] - r["t0"]) * r["vol"]
                    else:
                        v = r["cum1"]
                    break
            cement_out = _cement_out_volume(cement_steps, v, pipe_vol)
            eta_samples[tag] = (tt, v, cement_out, min(cement_out, ann_vol) / ann_vol if ann_vol > 0 else float("nan"))

        # F2 判定：末段水泥尾缘(修正后) vs stop
        f2_trig = (rear_g is not None) and (rear_g > stop + 1e-6)
        q_last_cem = step_rows[last_cem_i]["rate"] if last_cem_i is not None else 0.0
        # 截流体积两个口径：
        #  (a) 时间口径: (rear_g - stop) × Q_尾浆（重力修正时间轴上的错位量）
        #  (b) 体积口径(更硬): 水泥总量 - stop 时刻活塞流已出鞋水泥 = 未入库尾浆
        trapped_t = (rear_g - stop) * q_last_cem / 60.0 if f2_trig else 0.0
        stop_vol = eta_samples["stop"][1]
        trapped_v = cement_total - _cement_out_volume(cement_steps, stop_vol, pipe_vol)
        f2_real = trapped_v > 1e-3

        # 替浆步（DISPLACEMENT）总体积与模型管容
        disp_vol = sum(r["vol"] for r in step_rows if r["role"] == "DISPLACEMENT")
        front_noncem_after_cem = stop  # stop 本身即首个非水泥前缘（重力修正口径）
        # 首个非水泥前缘对应的累计泵量（活塞流下 = pipe_vol + cum0 首个非水泥步）
        first_disp_i = next((i for i, r in enumerate(step_rows) if r["role"] == "DISPLACEMENT"), None)

        print(f"\n### {label}  {well.well_name}")
        print(f"  模型管容={pipe_vol:.2f} m3 | 日程总泵量={total_vol:.2f} m3 | 日程总时长={total_t_sched:.0f}s | 环空体积={ann_vol:.1f} m3")
        print(f"  {'步':<16}{'角色':<14}{'vol':>6}{'rate':>6}{'cum0':>8}{'cum1':>8}{'前缘(活塞)s':>12}{'前缘(重力)s':>12}{'尾缘(活塞)s':>12}{'尾缘(重力)s':>12}")
        for i, r in enumerate(step_rows):
            is_last_cem = (i == last_cem_i)
            rear_p_s = f"{rear_p:.0f}" if (is_last_cem and rear_p is not None) else ""
            rear_g_s = f"{rear_g:.0f}" if is_last_cem else ""
            print(f"  {r['name']:<16}{r['role']:<14}{r['vol']:>6.1f}{r['rate']:>5.2f}{r['cum0']:>8.1f}{r['cum1']:>8.1f}"
                  f"{r['front_piston']:>12.0f}{r['front_grav']:>12.0f}{rear_p_s:>12}{rear_g_s:>12}")
        print(f"  -> stop(runner口径)={stop:.0f}s | cement_end_time_s(尾缘重力)={rear_g:.0f}s"
              f" | RR tt=min({total_t_sched:.0f}+1200, {stop:.0f}+600)={rr_tt:.0f}s")
        print(f"  -> F2 时间错位={rear_g - stop:+.0f}s (尾缘{rear_g:.0f}s vs stop{stop:.0f}s, Q={q_last_cem:.2f}m3/min)"
              f" -> 时间口径截流={trapped_t:.2f} m3")
        print(f"  -> F2 实际截流(体积口径)=水泥总{cement_total:.1f} - stop已出鞋{eta_samples['stop'][2]:.1f} = {trapped_v:.2f} m3"
              f" (占水泥{trapped_v / cement_total * 100:.1f}%)" if cement_total > 0 else "")
        print(f"  -> stop 时累计泵量={eta_samples['stop'][1]:.1f} m3 | 距日程结束还差={total_vol - eta_samples['stop'][1]:.1f} m3"
              f" | 替浆步总体积={disp_vol:.1f} m3 vs 模型管容={pipe_vol:.1f} m3")
        print("  -> 三评价时刻已出鞋口水泥体积 / eta_E 采样:")
        for tag in ("stop", "stop+600", "sched+1200"):
            tt, v, cem_out, eta = eta_samples[tag]
            print(f"     {tag:<11} t={tt:8.0f}s  泵量={v:7.1f}m3  出鞋水泥={cem_out:7.1f}m3  eta_E近似={eta:.4f}")

        # hu101 停泵后浮力交换流锚点
        if label == "hu101":
            mud = by.get("钻井液")
            lead = by.get("领浆")
            tail = by.get("尾浆")
            if mud is not None and lead is not None:
                # 重驱轻不稳定界面：环空内水泥(领浆 2.10) vs 泥浆(1.96)
                delta_rho = (lead.density_kg_m3 - mud.density_kg_m3) / 1000.0  # g/cc
                # 窄边间隙：139.7 段 hole 215.9，standoff 取 loader 名义剖面 0.42-0.46 的代表 0.44
                hole_mm = 215.9
                H = (hole_mm - well.liner_od_mm) / 2000.0  # m
                standoff = 0.44
                b_narrow = 2.0 * H * standoff
                # 有效黏度：尾浆幂律 K,n 在窄边剪切率量级 1/s 下 -> mu_eff ~ K * gamma^n, gamma~1-10
                n = tail.power_law_n
                K = tail.consistency_k
                mu_low = K * 1.0**n    # gamma=1
                mu_high = K * 10.0**n  # gamma=10
                L_cem = 2468.0  # 尾管封固段长度 5400-7868m
                g = 9.81
                print("  -> [5bis-4] 停泵后浮力交换流特征时间锚点（标度估计，非模拟）:")
                print(f"     Δρ={delta_rho:.2f} g/cc | H(单边间隙)={H*1000:.1f}mm | standoff={standoff} | b_narrow={b_narrow*1000:.1f}mm")
                print(f"     尾浆幂律 n={n:.3f} K={K:.3f} Pa.s^n -> mu_eff(γ=1)={mu_low:.2f} Pa.s | mu_eff(γ=10)={mu_high:.2f} Pa.s")
                for mu in (mu_low, mu_high):
                    # 窄槽 Poiseuille 标度: u ~ Δρ g b^2 / (12 μ) (单位密度差驱动梯度 Δρ g)
                    u = delta_rho * 1000.0 * g * b_narrow**2 / (12.0 * mu)
                    t_ex = L_cem / u if u > 0 else float("inf")
                    print(f"     mu={mu:.2f} Pa.s -> u_buoy~{u:.2e} m/s -> t_exchange=L/u~{t_ex:.3e} s ({t_ex/3600:.1f} h)")
        _ = first_disp_i  # 保留变量以备查
    print("=" * 100)


if __name__ == "__main__":
    main()
