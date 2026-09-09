"""开关开验证：hu103 七界面 σ_on > σ_off + 尾浆 pin"""
import math
from cemdisp.data.loaders import load_hu103_tailpipe
from cemdisp.transport1d.casing_flow import CasingFlowSolver
from cemdisp.transport1d.shoe_timeline import ShoeEventKind

well, fluids, schedule, _ = load_hu103_tailpipe()
off = CasingFlowSolver(enable_gravity=True)
on = CasingFlowSolver(enable_gravity=True, mixing_contact_time=True)
r_off = off.run(well, fluids, schedule)
r_on = on.run(well, fluids, schedule)

def first_fracs(events):
    # 收集每个过渡带的 (next_fluid, prev_fluid, sigma)：以 frac 从 <1 首次升序段定位
    out = []
    i = 0
    while i < len(events):
        e = events[i]
        if e.kind == ShoeEventKind.FRONT_ARRIVAL and len(e.phase_fractions) == 2 and e.phase_fractions[0][1] < 1.0:
            nxt, prv = e.phase_fractions[0][0], e.phase_fractions[1][0]
            # 子事件间隔 = 2σ/(n_sub-1) = σ/2 → σ = t5-t1 /2 ... 直接用同流体段跨度/2
            j = i
            while j < len(events) and events[j].kind == ShoeEventKind.FRONT_ARRIVAL and events[j].phase_fractions[:1] and events[j].phase_fractions[0][0] == nxt and events[j].phase_fractions[1][0] == prv:
                j += 1
            t0, t1 = events[i].time_s, events[j-1].time_s
            out.append((nxt, prv, (t1 - t0) / 2.0, t0))
            i = j
        else:
            i += 1
    return out

fo, fn = first_fracs(r_off.shoe_timeline.events), first_fracs(r_on.shoe_timeline.events)
print(f"bands off={len(fo)} on={len(fn)}")
for (n_o, p_o, s_o, _), (n_n, p_n, s_n, _) in zip(fo, fn):
    assert (n_o, p_o) == (n_n, p_n)
    print(f"{n_o:>6}<-{p_o:<6} sigma_off={s_o:.4f} sigma_on={s_n:.4f} ratio={s_n/s_o:.6f} {'WIDER' if s_n>s_o else 'narrower'}")
# 尾浆 pin
r_liner = (well.liner_id_mm or 100.0)/2000.0
expect = math.sqrt(88.55/(well.shoe_md_m*math.pi*r_liner**2))
tail = [(n,p,s) for n,p,s,_ in fn if n=="尾浆" and p=="中间浆"][0]
tail_o = [(n,p,s) for n,p,s,_ in fo if n=="尾浆" and p=="中间浆"][0]
got = tail[2]/tail_o[2]
print(f"pin: sigma_on/sigma_off={got:.6f} expected sqrt(88.55/71.101)={expect:.6f} rel={abs(got-expect)/expect:.2e}")
