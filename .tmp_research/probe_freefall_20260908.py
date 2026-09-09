# -*- coding: utf-8 -*-
"""自由下落（U 型管）效应对 8 井的量级评估探针（只读，不改任何模型代码）。

== 目的 ==
回答："忽略自由下落到底会错多少"——对 8 口尾管井逐一估算
  (a) 静压盈余上界 Δp_static：管内水泥柱 vs 环空同深度液柱的静压差上界
      （"若无摩阻、无背压"的全自由下落驱动力上界）；
  (b) 流动摩阻 Δp_fric(Q)：设计排量下管内+环空 Darcy-Weisbach 粗估；
  (c) 自由下落附加流量 Q_ff：把静压盈余全部消耗在"排量从 Q_pump 增到
      Q_pump+Q_ff 的附加摩阻"上，数值反解 Q_ff（自由下落窗口内的流量放大上界，
      Beirute 1984 口径："free-fall rate is larger than the surface pump rate"）；
  (d) 现场控压值 P_choke（0708 原件施工记录）对自由下落窗口的抑制判定：
      Δp_surplus_eff = Δp_static − P_choke，≤0 → 控压完全抑制（停泵即静止，
      环空闭环下返出必须顶开背压）；
  (e) 现有 Atwood 重力修正（casing_flow._gravity_corrected_arrival_time，
      enable_gravity=True 生产默认）的实际到达时刻偏移量：
      跑 1D solver 开/关重力修正对比 cement_end_time_s（=碰压断面 stop）；
  (f) 理论到达时刻偏移上界：若 (d) 判定自由下落窗口存在，
      界面到达提前量 ~ (1 − Q_pump/(Q_pump+Q_ff)) × 活塞流到达时间。

== 全部公式（可复核）==

1. 管内截面积 A_pipe = shoe_lag_volume_m3 / shoe_md_m（与生产
   CasingFlowSolver._pipe_cross_section_area 同口径；shoe_lag 缺失时用
   liner_id_mm 圆面积）。等效管内径 D_pipe = sqrt(4·A_pipe/π)。

2. 管内水泥柱最大长度 h_cem_max = V_cement_total / A_pipe
   （水泥浆总体积全部进入管内、前缘尚未出鞋口的时刻，静压盈余最大）。

3. 静压盈余上界：
      Δp_static = max_i( |ρ_cem_i − ρ_mud0| ) · g · h_cem_max
   其中 ρ_mud0 = 初始钻井液密度（被顶替液），ρ_cem_i 为各水泥浆（含隔离液
   保守不计——隔离液密度差更小）。取 |Δρ| 是因为轻驱重（尾浆轻于泥浆）时
   U 型管反向（管内减速、到达滞后），量级同阶。

4. Darcy-Weisbach 摩阻（牛顿等效黏度粗估）：
      管内：v = Q/A_pipe；μ_eff = μ_p + τ_y·D_pipe/(6·v)（Bingham 有效黏度，
            幂律取 μ_eff = K·γ̇^(n−1)，γ̇ = 8v/D·(3n+1)/(4n)）；
            Re = ρ·v·D_pipe/μ_eff；f = 64/Re（层流假设，量级评估足够）；
            Δp_pipe = f·(L_pipe/D_pipe)·(ρ·v²/2)，L_pipe = shoe_md_m（全管长）。
      环空：评价段 hanger→shoe，D_out = 裸眼平均井径（hole_diameter_profile
            在该段的算术平均），D_in = liner_od_mm；
            A_ann = π(D_out²−D_in²)/4；v = Q/A_ann；
            水力直径 D_h = D_out − D_in（环形空间，层流 f=96/Re 用缝宽近似）；
            Δp_ann = 96·μ_eff·L_ann·v/( (D_out−D_in)² · ρ )·ρ … 直接用层流平板近似：
            Δp_ann = 12·μ_eff·L_ann·v / ((D_out−D_in)/2)²（无限平板缝 Poiseuille，
            量级评估足够；含 ρ 的形式见下）。
      μ_eff 取管内被顶泥浆与环空泥浆同参数（保守取泥浆 Bingham 参数）。
   局限声明：未做流态判别（Re>2000 应转紊流 f=0.316/Re^0.25，本探针统一层流
   上包络——层流摩阻随 Q 线性、紊流随 Q²，反解 Q_ff 用层流口径会低估紊流工况的
   摩阻约束、高估 Q_ff，因此结论按"上界"解读）。

5. 自由下落附加流量反解：
      Δp_fric(Q) = Δp_pipe(Q) + Δp_ann(Q)（层流下线性于 Q）
      Q_ff: Δp_fric(Q_pump + Q_ff) − Δp_fric(Q_pump) = Δp_static − P_choke
      层流线性 → 闭式解 Q_ff = Q_pump · (Δp_static − P_choke)/Δp_fric(Q_pump)
      （Δp_static−P_choke ≤ 0 → Q_ff = 0，自由下落被抑制）。
   注意 P_choke（井口环空背压）同时抬高管内泵压与环空静压，但 U 型管驱动是
   "管内液柱 vs 环空液柱+井口背压"的失衡，背压直接吃掉静压盈余（闭环 RCD 下
   任何额外返出必须顶开节流阀压差）。

6. 到达时刻偏移：
      现有 Atwood 修正实测偏移 = |t_gravity_on − t_gravity_off|（1D 开关对比）；
      理论上界 = t_piston · Q_ff/(Q_pump+Q_ff)（若自由下落窗口全程存在——
      这是最激进上界，实际窗口只覆盖水泥在管内行进的时段）。

== 输出 ==
results/_自由下落量级_2026-09-08/freefall_magnitude.csv   每井一行全指标
results/_自由下落量级_2026-09-08/README.md                小结
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_ROOT / "results" / "_自由下落量级_2026-09-08"

G = 9.81


# ---------------------------------------------------------------- 现场控压值
# 0708 原件取证（2026-09-08）：施工期（顶替窗口内）有直接"控压值"记录的井仅
# hu102（20213 施工记录表控压值列 1.0→15→30 MPa）与 hu103（20313"控压至
# 6120m ECD 2.05 控压值 3.4 MPa"）。hu1 204121"环空加压 3.2-11MPa/30h"与
# 20326/202392 系的 14.6/12 MPa 均为**候凝期**加压（碰压后），不在顶替窗口内。
# 判据取最保守口径：所有井施工控压按 0（不靠背压抑制）；有据值单列展示。
FIELD_CHOKES_MPA = {}          # 判据用：全部按无控压（最保守）
FIELD_CHOKES_EVIDENCE = {      # 展示用：0708 原件有据的控压/加压值（注明窗口）
    "hu1":      "候凝 3.2–11 MPa/30h（204121）；施工原则'停泵控压、开泵不控压'",
    "hu102":    "施工期控压值列 1.0→15→30 MPa（20213）；候凝环空加压",
    "hu103":    "施工期控压 3.4 MPa（20313）；碰压后环空憋压 6.9；候凝 5.3 MPa/48h（20314）",
    "ht1_004":  "候凝环空加压 14.6 MPa/72h（20326/20327，井归属待核）",
    "hu101":    "无施工控压值记录（loader/0708 未见）",
    "hu2":      "无施工控压值记录",
    "ht1_001":  "无施工控压值记录",
    "ht1_003":  "无施工控压值记录（202392 候凝加压 12 MPa 井归属待核）",
}

CEMENT_ROLES = {"lead", "tail", "intermediate", "cement"}


def _fluid_mu_eff_pa_s(fluid, v_m_s: float, d_m: float) -> float:
    """牛顿等效黏度（量级评估口径）。Bingham: μ_p + τ_y·D/(6v)；
    幂律: K·γ̇^(n−1)；缺流变回退 30 mPa·s。"""
    if fluid is None:
        return 0.03
    if fluid.rheology_model.value == "bingham" and fluid.plastic_viscosity_pa_s is not None:
        mu_p = fluid.plastic_viscosity_pa_s
        tau_y = fluid.yield_stress_pa or 0.0
        v = max(v_m_s, 1e-3)
        return mu_p + tau_y * d_m / (6.0 * v)
    if fluid.rheology_model.value == "power_law" and fluid.power_law_n is not None:
        n = fluid.power_law_n
        k = fluid.consistency_k or 0.1
        v = max(v_m_s, 1e-3)
        gamma = 8.0 * v / d_m * (3.0 * n + 1.0) / (4.0 * n)
        return k * gamma ** (n - 1.0)
    return 0.03


def _mean_hole_diameter_m(well_spec, from_md: float, to_md: float) -> float | None:
    pts = [p.value for p in well_spec.hole_diameter_profile if from_md <= p.depth_md_m <= to_md]
    if not pts:
        return None
    return sum(pts) / len(pts) / 1000.0


def analyze_well(tag: str, loader_kwargs: dict | None = None) -> dict:
    from cemdisp.data import loaders as L
    from cemdisp.transport1d.casing_flow import CasingFlowSolver
    from cemdisp.data.fluid_spec import FluidRole

    fn = getattr(L, f"load_{tag}_tailpipe")
    try:
        well, fluids, schedule, _v = fn(**(loader_kwargs or {}))
    except TypeError:
        well, fluids, schedule, _v = fn()

    fluid_by_name = {f.name: f for f in fluids}
    role_of = {f.name: f.role.value for f in fluids}

    # --- 管内几何（与生产 _pipe_cross_section_area 同口径） ---
    if well.shoe_lag_volume_m3:
        a_pipe = well.shoe_lag_volume_m3 / well.shoe_md_m
    elif well.liner_id_mm:
        a_pipe = math.pi * (well.liner_id_mm / 1000.0) ** 2 / 4.0
    else:
        a_pipe = math.pi * 0.1079 ** 2 / 4.0
    d_pipe = math.sqrt(4.0 * a_pipe / math.pi)

    # --- 泥浆基准密度（初始被顶替液 = MUD 角色第一个；无则取密度最小泥浆系） ---
    mud = next((f for f in fluids if f.role.value == "mud"), None)
    rho_mud = mud.density_kg_m3 if mud else min(f.density_kg_m3 for f in fluids)

    # --- 水泥浆总量与最大密度差 ---
    cement_vols = [(s.volume_m3, fluid_by_name[s.fluid_name]) for s in schedule.steps
                   if role_of.get(s.fluid_name) in CEMENT_ROLES]
    v_cem_total = sum(v for v, _ in cement_vols)
    if cement_vols:
        max_drho = max(abs(f.density_kg_m3 - rho_mud) for _v_, f in cement_vols)
    else:
        max_drho = 0.0

    h_cem_max = v_cem_total / a_pipe if a_pipe > 0 and v_cem_total > 0 else 0.0
    dp_static = max_drho * G * h_cem_max  # Pa（激进上界：全长 × 最大密度差）

    # --- 时序感知静压盈余峰值（更物理的口径） ---
    # 管内液柱按泵注顺序分层：注入体积坐标 V 时，管内占据 (V−V_pipe, V] 的流体层。
    # 静压盈余(V) = Σ_层 (ρ_层 − ρ_mud)·g·层长 = (p_casing(shoe) − ρ_mud·g·L)。
    # 真实自由下落窗口 = 该盈余 > 背压+摩阻 的时段（U 型管驱动力逐时刻平衡式，
    # Beirute 1984 口径）。峰值取 max_V。
    def dp_surplus_at(V: float) -> float:
        lo = max(V - a_pipe * well.shoe_md_m, 0.0)  # 管内窗口下沿（体积坐标）
        surplus = 0.0
        cum = 0.0
        for s in schedule.steps:
            if s.rate_m3_min <= 0:
                continue  # RESTART/停泵步不推界面（与生产胶塞语义一致）
            seg_hi = min(cum + s.volume_m3, V)
            seg_lo = max(cum, lo)
            if seg_hi > seg_lo:
                fl = fluid_by_name.get(s.fluid_name)
                rho = fl.density_kg_m3 if fl else rho_mud
                surplus += (rho - rho_mud) * G * (seg_hi - seg_lo) / a_pipe
            cum += s.volume_m3
            if cum >= V:
                break
        return surplus

    v_total = sum(s.volume_m3 for s in schedule.steps if s.rate_m3_min > 0)
    dp_static_seq = 0.0
    n_scan = 200
    for i in range(1, n_scan + 1):
        v_scan = v_total * i / n_scan
        dp_static_seq = max(dp_static_seq, dp_surplus_at(v_scan))
    dp_static = max(dp_static, 0.0)
    dp_static_peak = dp_static_seq  # 时序感知峰值

    # --- 设计（最大）排量 ---
    q_pump = max(s.rate_m3_min for s in schedule.steps if s.rate_m3_min > 0) / 60.0  # m3/s

    # --- 管内摩阻 @Q_pump（泥浆在被顶段的等效黏度） ---
    v_pipe = q_pump / a_pipe
    mu_pipe = _fluid_mu_eff_pa_s(mud, v_pipe, d_pipe)
    re_pipe = (rho_mud or 2000.0) * v_pipe * d_pipe / mu_pipe
    f_darcy = 64.0 / max(re_pipe, 1.0)
    dp_pipe = f_darcy * (well.shoe_md_m / d_pipe) * (rho_mud or 2000.0) * v_pipe ** 2 / 2.0

    # --- 环空摩阻 @Q_pump（平板缝层流近似） ---
    l_ann = well.shoe_md_m - (well.hanger_md_m or well.top_md_m)
    d_hole = _mean_hole_diameter_m(well, well.hanger_md_m or well.top_md_m, well.shoe_md_m)
    d_od = (well.liner_od_mm or 139.7) / 1000.0
    if d_hole and d_hole > d_od:
        gap = (d_hole - d_od) / 2.0
        a_ann = math.pi * (d_hole ** 2 - d_od ** 2) / 4.0
        v_ann = q_pump / a_ann
        mu_ann = _fluid_mu_eff_pa_s(mud, v_ann, d_hole - d_od)
        dp_ann = 12.0 * mu_ann * l_ann * v_ann / gap ** 2
    else:
        a_ann = v_ann = dp_ann = 0.0

    dp_fric = dp_pipe + dp_ann

    # --- 自由下落窗口判据（施工控压按 0 最保守；有据控压值见 EVIDENCE） ---
    p_choke = FIELD_CHOKES_MPA.get(tag, 0.0) * 1e6
    dp_surplus_noc = dp_static_peak - dp_fric            # 无控压窗口判据（>0 → 会自由下落）
    dp_surplus_mpd = dp_static_peak - p_choke - dp_fric  # 同上（本探针两口径合一）
    # 层流线性闭式反解 Q_ff（控压口径）
    if dp_surplus_mpd > 0 and dp_fric > 0:
        q_ff = q_pump * dp_surplus_mpd / dp_fric
    else:
        q_ff = 0.0
    q_ratio = q_ff / q_pump if q_pump > 0 else 0.0

    # --- 现有 Atwood 修正实测偏移（1D 开关对比） ---
    res_on = CasingFlowSolver(enable_gravity=True).run(well, fluids, schedule)
    res_off = CasingFlowSolver(enable_gravity=False).run(well, fluids, schedule)
    t_stop_on = res_on.cement_end_time_s
    t_stop_off = res_off.cement_end_time_s
    stop_shift_s = t_stop_on - t_stop_off
    # 活塞流到达时刻基准（off 口径）与理论上界
    t_piston = t_stop_off
    t_shift_upper = t_piston * q_ratio / (1.0 + q_ratio) if q_ratio > 0 else 0.0

    return {
        "井号": tag,
        "井深_m": well.shoe_md_m,
        "管容_m3": round(a_pipe * well.shoe_md_m, 1),
        "水泥浆总量_m3": round(v_cem_total, 1),
        "最大密度差_kgm3": round(max_drho, 0),
        "管内水泥柱长上界_m": round(h_cem_max, 0),
        "静压盈余上界_MPa": round(dp_static / 1e6, 3),
        "时序感知静压峰值_MPa": round(dp_static_peak / 1e6, 3),
        "现场控压证据": FIELD_CHOKES_EVIDENCE.get(tag, "无"),
        "控压值_MPa": round(p_choke / 1e6, 1),
        "摩阻管内_MPa": round(dp_pipe / 1e6, 3),
        "摩阻环空_MPa": round(dp_ann / 1e6, 3),
        "摩阻合计_MPa": round(dp_fric / 1e6, 3),
        "无控压盈余_MPa": round(dp_surplus_noc / 1e6, 3),
        "控压后盈余_MPa": round(dp_surplus_mpd / 1e6, 3),
        "自由下落附加流量_m3min": round(q_ff * 60.0, 3),
        "流量放大倍数_QpumpQff": round(1.0 + q_ratio, 3),
        "stop_无重力_s": round(t_stop_off, 0),
        "stop_有重力_s": round(t_stop_on, 0),
        "Atwood修正偏移_s": round(stop_shift_s, 0),
        "Atwood修正偏移_pct": round(stop_shift_s / t_stop_off * 100.0, 2) if t_stop_off else 0.0,
        "理论上界偏移_s": round(t_shift_upper, 0),
        "环空流速峰值低估因子": round(q_ratio, 3),
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for tag in ["hu1", "hu101", "hu102", "hu103", "hu2", "ht1_001", "ht1_003", "ht1_004"]:
        try:
            rows.append(analyze_well(tag))
        except Exception as exc:  # 单井失败不阻断
            rows.append({"井号": tag, "错误": repr(exc)})
            raise

    cols = list(rows[0].keys())
    with open(OUT_DIR / "freefall_magnitude.csv", "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    lines = [
        "# 自由下落（U 型管）量级评估小结 — 8 井（2026-09-08 探针）",
        "",
        "口径：静压盈余两个口径——(A) 激进上界 = max|Δρ|·g·(水泥总量/管内截面积)；",
        "(B) 时序感知峰值 = 按泵注顺序逐时刻积分管内实际分层密度 vs 初始泥浆的静压盈余（物理自洽）。",
        "摩阻 = Darcy-Weisbach 层流下包络（管内 f=64/Re + 环空平板缝 12μLv/gap²）。",
        "判据：Δp_surplus = Δp_static − P_choke − Δp_fric ≤ 0 → 自由下落窗口关闭。",
        "控压判据取最保守（全部按 P_choke=0）；0708 原件有据控压值单列在 csv。",
        "双向保守声明：盈余按上界（环空全泥浆+忽略凝胶抗力）、摩阻按下界（层流，",
        "紊流真实摩阻更高）——即便如此 8 井窗口仍全关，结论对口径误差稳健。",
        "",
        "| 井号 | 时序感知静压峰值 MPa | 摩阻合计 MPa | 无控压盈余 MPa | 判定 |",
        "|---|---|---|---|---|",
    ]
    for r in rows:
        if "错误" in r:
            lines.append(f"| {r['井号']} | 错误 {r['错误']} | | | |")
            continue
        verdict = "窗口关闭" if r["无控压盈余_MPa"] <= 0 else "存在自由下落窗口（需复查）"
        lines.append(
            f"| {r['井号']} | {r['时序感知静压峰值_MPa']} | "
            f"{r['摩阻合计_MPa']} | {r['无控压盈余_MPa']} | {verdict} |"
        )
    lines += ["", "## 到达时刻偏移（现有 Atwood 修正 vs 理论上界）", "",
              "| 井号 | stop无重力 s | stop有重力 s | Atwood偏移 s (%) | 理论上界偏移 s | 流量放大倍数 |",
              "|---|---|---|---|---|---|"]
    for r in rows:
        if "错误" in r:
            continue
        lines.append(
            f"| {r['井号']} | {r['stop_无重力_s']} | {r['stop_有重力_s']} | "
            f"{r['Atwood修正偏移_s']} ({r['Atwood修正偏移_pct']}%) | {r['理论上界偏移_s']} | "
            f"{r['流量放大倍数_QpumpQff']} |"
        )
    (OUT_DIR / "README.md").write_text("\n".join(lines), encoding="utf-8")
    print("written:", OUT_DIR)
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
