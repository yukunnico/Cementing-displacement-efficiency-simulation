# -*- coding: utf-8 -*-
"""Phase 3.2（P-4 产品层）驱动：三井泵压/ECD/控压窗报告 + 稠化窗校核（D7-A）数据层。

依据 = docs/superpowers/specs/2026-10-07-phase3-hydraulics-p2-p4-design.md §2 +
总纲《详细执行计划_真温压响应_2026-10-06》§3 Phase 3.2。

A) 呼1-004：run_ht1004_target（靶权威 = 2026-10-07 沙箱 -batch 运行，wrapper_p30.m 捕获，
   原脚本零改动）→ 全链时程 CSV + summary；
B) 呼1-003/呼101：fluid_table_from/build_pump_schedule_custom/build_segments_custom +
   run_well_forward 正演（两井无逐时控压窗设计数据 ⇒ 只出泵压/ECD + 对可得当量口径余量）；
C) 稠化窗校核（D7-A）数值层；D) 作业史对账数值层（notes 泵压文字抽取 + 内核段映射）。

输出目录 = results/P4产品_压力全链_20261007/。只读输入、只写本目录，不触碰生产链。
"""
from __future__ import annotations

import csv
import io
import json
import math
import re
from pathlib import Path

import numpy as np

from cemdisp.transport1d import hydraulics as H

ROOT = Path(__file__).resolve().parents[2]
SANDBOX = ROOT / "results/_probe_matlab靶_2026-10-07/sandbox"
FIELD = ROOT / "参考文档/现场资料提取"
OUT = ROOT / "results/P4产品_压力全链_20261007"

STRUCT_COLS = list(H._STRUCTURE_COLUMNS)


def t_static_c(z_m):
    """静温剖面口径（任务书/总纲给定）：T(z)=16.006+1.7598e-2*z [degC]。"""
    return 16.006 + 1.7598e-2 * z_m


def fnum(x):
    if x is None:
        return None
    s = str(x).strip()
    if s == "":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def read_csv(path):
    return list(csv.DictReader(io.open(str(path), encoding="utf-8-sig")))


def write_series_csv(path, cols_units, rows):
    with io.open(str(path), "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols_units)
        for r in rows:
            w.writerow(r)


# --------------------------------------------------------------------------
# A) 呼1-004 全链（靶配置）
# --------------------------------------------------------------------------

def part_ht1004():
    res = H.run_ht1004_target(SANDBOX)
    s = res["summary"]
    nt = int(s["n_time"])
    t = np.arange(1, nt + 1, dtype=float)
    pps = res["out_pump_pressure_surface"]
    abp = res["out_annuli_bottom_pressure"]
    ecd_bot = res["out_result_matrix_volume_pressure"][:, 2]
    cumvol = res["cumulative_volume_m3"]
    write_series_csv(OUT / "ht1004_泵压_井底压力_ECD_时程.csv",
                     ["time_min", "pump_pressure_MPa", "annulus_bottom_pressure_MPa",
                      "bottom_ECD_g_cm3", "cumulative_volume_m3"],
                     zip(t, pps[:, 1], abp[:, 1], ecd_bot, cumvol))
    fp = res["out_four_point_backpressure"]
    fe = res["out_four_point_ecd_new"]
    write_series_csv(OUT / "ht1004_四点走廊与控压后ECD.csv",
                     ["time_min", "required_backpressure_MPa", "BP_lower_MPa", "BP_upper_MPa",
                      "conflict_flag", "post_ctrl_ECD_5578shoe_g_cm3",
                      "post_ctrl_ECD_6600_g_cm3", "post_ctrl_ECD_7498_g_cm3",
                      "post_ctrl_ECD_bottom_g_cm3"],
                     zip(fp[:, 0], fp[:, 1], fp[:, 2], fp[:, 3], fp[:, 4],
                         fe[:, 1], fe[:, 2], fe[:, 3], fe[:, 4]))
    wbc = res["out_window_backpressure_ctrl"]
    bpa = res["bp_apply_win_MPa"]
    write_series_csv(OUT / "ht1004_井底窗持稳回压.csv",
                     ["time_min", "ECD_bottom_base_g_cm3", "ECD_bottom_ctrl_g_cm3",
                      "ECD_shoe_ctrl_g_cm3", "P_wh_MPa", "bp_apply_MPa"],
                     zip(wbc[:, 0], wbc[:, 1], wbc[:, 2], wbc[:, 3], wbc[:, 4], bpa))
    return res, s


# --------------------------------------------------------------------------
# B) 通用井输入构造（合成井表 + 参数化装配）
# --------------------------------------------------------------------------

def interp(rows, xcol, ycol, xq):
    pts = sorted((float(r[xcol]), float(r[ycol])) for r in rows
                 if fnum(r.get(xcol)) is not None and fnum(r.get(ycol)) is not None)
    if not pts:
        return None
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    return float(np.interp(xq, xs, ys)), pts[0][0], pts[-1][0]


def build_grid_md(bottom_m, thresholds, step=30.0):
    bset = {0.0, float(bottom_m)}
    cur = step
    while cur < bottom_m:
        bset.add(round(cur, 6))
        cur += step
    bset.update(round(float(x), 6) for x in thresholds)
    inner = sorted(b for b in bset if 0.0 < b < bottom_m)
    edges = [0.0] + inner + [float(bottom_m)]
    return [(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]


def seg_index_of(md, thresholds):
    for j, t in enumerate(thresholds):
        if md <= t:
            return j
    return len(thresholds)


def make_structure(well, bottom_m, thresholds, bore_const_mm, pipe_od_mm, pipe_wall_mm,
                   dev_rows, cal_rows):
    """合成呼1-004 同款 8 列井表并加载为 StructureTable。

    口径声明（写入报告）：
    - bore_const_mm[k]=None 段 = 裸眼段，井径 = caliper_profile 线性插值（端点钳位）；
      常数段 = 套管内径/名义井眼，**直径口径**（legacy "radius 列实为直径"陷阱已在参数名消歧）；
    - deg = inclination_profile 插值（段末 MD），测斜首点以上 = 0.0（设计书直井段口径）；
    - sq/vol 列 = 几何式 (pi/4)(bore^2 - out^2)*len，与靶表生成口径一致（靶行 1 反推核验：
      245.37/149.2 -> 2.980256 dm2 -> 894.0768 L）；tvd_csv = 几何 cumsum（非电测 TVD 表，声明）。
    """
    dev0 = interp(dev_rows, "md_m", "inclination_deg", 0.0)
    dev_first = dev0[1] if dev0 else float("inf")
    grid = build_grid_md(bottom_m, thresholds)
    cd, tvd_geo, lens, dia_cm, vol, deg_arr, sq = [], [], [], [], [], [], []
    acc_tvd = 0.0
    for (a, b) in grid:
        ln = b - a
        k = seg_index_of(b, thresholds)
        if bore_const_mm[k] is not None:
            bore_mm = float(bore_const_mm[k])
        else:
            got = interp(cal_rows, "md_m", "caliper_mm", b)
            bore_mm = got[0] if got is not None else float(pipe_od_mm[k]) + 20.0
        gotd = interp(dev_rows, "md_m", "inclination_deg", b)
        dv = 0.0 if (gotd is None or b < dev_first) else gotd[0]
        bore_m = bore_mm / 1000.0
        out_m = float(pipe_od_mm[k]) / 1000.0
        sq_dm2 = math.pi * (bore_m**2 - out_m**2) / 4.0 * 100.0
        acc_tvd += ln * math.cos(math.radians(dv))
        cd.append(b)
        tvd_geo.append(acc_tvd)
        lens.append(ln)
        dia_cm.append(bore_mm / 10.0)
        vol.append(sq_dm2 * ln * 10.0)
        deg_arr.append(dv)
        sq.append(sq_dm2)
    p = OUT / ("synth_井身结构_" + well + ".csv")
    with io.open(str(p), "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(STRUCT_COLS)
        for i in range(len(cd)):
            w.writerow([i + 1, "%.6g" % cd[i], "%.6g" % tvd_geo[i], "%.6g" % lens[i],
                        "%.9g" % dia_cm[i], "%.9g" % vol[i], "%.9g" % deg_arr[i],
                        "%.9g" % sq[i]])
    st = H.load_structure(p)
    segs = H.build_segments_custom(st, thresholds, bore_const_mm, pipe_od_mm, pipe_wall_mm)
    return st, segs, p


def mu_from_kn(k_pa_sn, n, gamma_ref=100.0):
    """声明换算：幂律 K,n -> 等效塑性粘度 mu_p [mPa.s] @ gamma_ref=100 1/s；tau_y=0。
    仅用于无 PV/YP 实测的浆（换算=建模口径，非实测值，报告中标注）。"""
    return k_pa_sn * (gamma_ref ** (n - 1.0)) * 1000.0


HT1003 = dict(
    well="ht1003_呼1-003", bottom_m=7618.0,
    # 变段分界：套管鞋 1404.746/3905/5568（geometry/casing_liner_string）、
    # 钻杆变径 4025.73（借靶 HT1-004 钻杆组合口径——缺失，声明）、尾管顶 5316.036、
    # 尾管变径 7091.016、钻头变径 7096（geometry borehole_nominal 两档）
    thresholds=[1404.746, 3905.0, 4025.73, 5316.036, 5568.0, 7091.016, 7096.0],
    bore=[475.74, 337.36, 245.37, 245.37, 245.37, None, None, None],
    pipe_od=[149.2, 149.2, 127.0, 127.0, 168.28, 168.28, 139.7, 139.7],
    pipe_wall=[9.65, 9.65, 9.65, 9.65, 14.7, 14.7, 15.88, 15.88],
    # 设计口径泵序（pumping_schedule design_ 行；前 5 = 环空流体系）
    fluid_names=["先导浆", "隔离液1", "隔离液2", "领浆", "尾浆", "压塞液",
                 "钻井液", "保护液", "钻井液", "钻井液", "钻井液", "钻井液"],
    rates=[1.2, 1.2, 1.2, 1.2, 1.4, 1.0, 1.6, 1.4, 1.2, 1.0, 0.8, 0.7],
    vols=[28, 16, 10, 39, 28, 2, 25, 14, 8, 14, 16, 12.9],
    rou=[1.75, 2.05, 1.95, 2.05, 1.95, 1.95, 1.95, 1.95, 1.95, 1.95, 1.95, 1.95],
    miu=[55.0, mu_from_kn(1.245, 0.668), mu_from_kn(1.245, 0.668),
         mu_from_kn(1.622, 0.597), mu_from_kn(1.673, 0.585), 30.0, 51.0, 30.0,
         51.0, 51.0, 51.0, 51.0],
    tau=[9.2, 0.0, 0.0, 0.0, 0.0, 8.0, 10.0, 8.2, 10.0, 10.0, 10.0, 10.0],
    mud=(1.95, 51.0, 10.0),
)

HU101 = dict(
    well="hu101_呼101", bottom_m=7868.0,
    # 变段分界：钻杆变径 4025.73（借靶口径——缺失，声明）、尾管顶 5402.885、
    # 裸眼顶 5700（套管鞋 5699.8）、尾管变径 6796.329、变扣 7048.83
    thresholds=[4025.73, 5402.885, 5700.0, 6796.329, 7048.83],
    # 0-5700 取名义井眼 241.3（centralizer 行记载井眼 241.3；回接套管/中间套管 ID
    # 未提取=数据缺口，且该段仅泥浆柱+宽环空摩擦，影响≈0，声明）
    bore=[241.3, 241.3, 241.3, None, None, None],
    pipe_od=[149.2, 127.0, 168.3, 168.3, 139.7, 139.7],
    pipe_wall=[9.65, 9.65, 13.0, 13.0, 14.27, 15.8],
    # 现场口径泵序（pumping_schedule 行 1-8；替浆末段 63.4m3@1.0-0.55 取均值 0.9——声明）
    fluid_names=["平衡液", "驱油隔离液", "领浆", "尾浆", "压塞液",
                 "轻泥浆(替浆)", "中置液", "钻井液"],
    rates=[1.2, 1.0, 1.0, 1.0, 0.6, 1.5, 1.0, 0.9],
    vols=[25, 25, 47, 23, 2, 26, 10, 63.4],
    rou=[1.85, 2.00, 2.10, 1.90, 2.00, 1.85, 2.00, 1.96],
    miu=[58.0, 30.0, mu_from_kn(0.815, 0.719), mu_from_kn(0.684, 0.722), 30.0,
         58.0, 30.0, 58.0],
    tau=[9.2, 0.0, 0.0, 0.0, 0.0, 9.2, 0.0, 9.2],
    mud=(1.96, 58.0, 9.2),
)


def part_generic(spec, wdir):
    dev = read_csv(FIELD / wdir / "inclination_profile.csv")
    cal = read_csv(FIELD / wdir / "caliper_profile.csv")
    st, segs, _ = make_structure(spec["well"], spec["bottom_m"], spec["thresholds"],
                                 spec["bore"], spec["pipe_od"], spec["pipe_wall"],
                                 dev, cal)
    fluids = H.fluid_table_from(
        [r * 1000.0 for r in spec["rates"]], [v * 1000.0 for v in spec["vols"]],
        spec["rou"], spec["miu"], spec["tau"], *spec["mud"])
    sched = H.build_pump_schedule_custom(fluids)
    res = H.run_well_forward(st, segs, fluids, sched)
    s = res["summary"]
    nt = int(s["n_time"])
    t = np.arange(1, nt + 1, dtype=float)
    write_series_csv(OUT / (spec["well"] + "_泵压_井底压力_ECD_时程.csv"),
                     ["time_min", "pump_pressure_MPa", "annulus_bottom_pressure_MPa",
                      "bottom_ECD_annulus_g_cm3", "bottom_ECD_casing_g_cm3",
                      "cumulative_volume_m3"],
                     zip(t, res["pump_MPa"], res["ann_bottom_MPa"],
                         res["ecd_bottom_annulus_g_cm3"], res["ecd_bottom_casing_g_cm3"],
                         res["cumvol_m3"]))
    # ECD 剖面选时：起泵、尾浆泵毕节点、末步
    n = int(s["n_segment"])
    nodes = res["sched_nodes_min"]
    idx_t = sorted({1, int(min(nt, math.floor(nodes[4]))), nt})
    md = st.c_depth_m
    with io.open(str(OUT / (spec["well"] + "_ECD剖面_选时.csv")), "w",
                 encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["md_m"] + ["ECD_annulus_g_cm3@t%02dmin" % tt for tt in idx_t])
        for i in range(n):
            w.writerow(["%.6g" % md[i]] +
                        ["%.6g" % res["ecd_ann_full"][tt - 1, i] for tt in idx_t])
    res["fluid_names"] = spec["fluid_names"]
    res["nodes_min"] = nodes
    return res


# --------------------------------------------------------------------------
# C) 稠化窗校核（D7-A）：泵注程序（呼1-004 设计口径 190 min）vs 井下 T,P
# --------------------------------------------------------------------------

def part_d7a(res_a):
    st = H.load_structure(SANDBOX / "呼1-004井身结构.csv")
    segs = H.build_segments(st)
    cap_L = segs.capacity_pipe_L
    fl = H.build_fluid_table()
    sc = H.build_pump_schedule(fl)
    nt = sc.n_time
    nodes = sc.nodes_min
    vol_all = sc.vol_all_L_1b

    def first_t_ge(thr_L):
        for t in range(1, nt + 1):
            if vol_all[t] >= thr_L:
                return float(t)
        return float("inf")

    t_lead_shoe = first_t_ge(cap_L + (fl.vols_L[0] + fl.vols_L[1] + fl.vols_L[2]))
    t_tail_shoe = first_t_ge(cap_L + sum(fl.vols_L[:4]))
    t_tail_clear = first_t_ge(cap_L + sum(fl.vols_L[:5]))
    p_bot = res_a["out_result_matrix_volume_pressure"][:, 4]
    tb = int(st.n_segment)
    p_shoe_series = None
    ecd_arr = res_a["ecd_annulus_g_cm3_full"]  # (nt,333) 0基
    t_bot_c = t_static_c(float(st.tvd_csv_m[-1]))
    t_shoe_c = t_static_c(5578.0)
    z110 = (110.0 - 16.006) / 1.7598e-2
    # 现场稠化（fluid_properties.csv 注记，70Bc 判据；132C 实验口径）
    fp_rows = read_csv(FIELD / "ht1_004_呼1-004" / "fluid_properties.csv")
    consist = {}
    for r in fp_rows:
        m = re.search(r"稠化(\d+)min/70Bc", r.get("notes", ""))
        if m and r["fluid_role"] in ("lead_cement", "tail_cement"):
            consist.setdefault(r["fluid_role"], int(m.group(1)))
    lead_c = consist.get("lead_cement", 489)
    tail_c = consist.get("tail_cement", 256)
    mix_lead = nodes[2]      # 领浆泵始（=混始近似，声明）
    mix_tail = nodes[3]      # 尾浆泵始
    bump_ref = 260.0         # 现场碰压（pumping_schedule 行 7 end_time_min）
    rows = [
        ["领浆", lead_c, lead_c / 3.0, "%.2f" % mix_lead, "%.1f" % t_lead_shoe,
         "%.0f" % bump_ref, "现场碰压",
         lead_c - (bump_ref - mix_lead), lead_c / 3.0 - (bump_ref - mix_lead),
         "%.1f" % t_shoe_c, "%.2f" % float(np.max(p_bot[int(mix_lead):min(nt, int(nodes[4]))]))],
        ["尾浆", tail_c, tail_c / 3.0, "%.2f" % mix_tail, "%.1f" % t_tail_shoe,
         "%.0f" % bump_ref, "现场碰压",
         tail_c - (bump_ref - mix_tail), tail_c / 3.0 - (bump_ref - mix_tail),
         "%.1f" % t_bot_c, "%.2f" % float(np.max(p_bot[int(mix_tail):nt]))],
    ]
    write_series_csv(OUT / "稠化窗校核_D7A_表.csv",
                     ["浆", "现场稠化_min@70Bc_高压132C143MPa口径", "D7A_1/3极端折算_min",
                      "混始_内核min", "浆前缘到鞋_内核min", "对照终点_min", "对照点说明",
                      "余量_field_min", "余量_extreme_min", "静温剖面T_at_鞋|底°C",
                      "内核环空井底P_MPa(max in窗)"], rows)
    return dict(cap_L=cap_L, nodes=nodes, t_lead_shoe=t_lead_shoe, t_tail_shoe=t_tail_shoe,
                 t_tail_clear=t_tail_clear, lead_c=lead_c, tail_c=tail_c,
                 T_bottom_static=t_bot_c, T_5578_static=t_shoe_c, z_T110_m=z110,
                 P_bottom_max=float(np.max(p_bot)))


# --------------------------------------------------------------------------
# D) 作业史对账：pumping_schedule notes 泵压文字抽取
# --------------------------------------------------------------------------

def part_notes(wdir):
    rows = read_csv(FIELD / wdir / "pumping_schedule.csv")
    out = []
    allv = []
    pat = re.compile(r"泵压([0-9.]+(?:[~\-–][0-9.]+)*)MPa")
    for r in rows:
        note = r.get("notes", "")
        hits = []
        for m in pat.finditer(note):
            for piece in re.split(r"[~\-–]", m.group(1)):
                if piece:
                    hits.append(float(piece))
        allv.extend(hits)
        out.append([r.get("step_index"), r.get("stage_name"), r.get("fluid_name"),
                    r.get("volume_m3"), r.get("rate_m3_min"), r.get("start_time_min"),
                    r.get("end_time_min"), r.get("density_g_cm3"),
                    "|".join("%.1f" % v for v in hits), note])
    return out, (min(allv) if allv else None), (max(allv) if allv else None)


# --------------------------------------------------------------------------
# 主流程
# --------------------------------------------------------------------------

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    summary = {"generated_utc": "2026-10-07", "design_ref":
               "docs/superpowers/specs/2026-10-07-phase3-hydraulics-p2-p4-design.md §2"}
    # A 呼1-004
    res_a, s_a = part_ht1004()
    pump_a = res_a["out_pump_pressure_surface"][:, 1]
    ecd_a = res_a["out_result_matrix_volume_pressure"][:, 2]
    cmp_p = res_a["out_pump_pressure_comparison"]
    cmp_e = res_a["out_bottom_ecd_comparison"]
    summary["ht1004"] = dict(s_a)
    summary["ht1004"]["ecd_ref_err_min_max_g_cm3"] = [float(np.min(cmp_e[:, 4])),
                                                       float(np.max(cmp_e[:, 4]))]
    summary["ht1004"]["pump_vs_hist_err_MPa_min_max_mean_std"] = [
        float(np.min(cmp_p[:, 4])), float(np.max(cmp_p[:, 4])),
        float(np.mean(cmp_p[:, 4])), float(np.std(cmp_p[:, 4]))]
    wbc = res_a["out_window_backpressure_ctrl"]
    bpa = res_a["bp_apply_win_MPa"]
    summary["ht1004"]["bp_apply_min_max_MPa"] = [float(np.min(bpa)), float(np.max(bpa))]
    summary["ht1004"]["shoe_ecd_ctrl_range"] = [float(np.min(wbc[1:, 3])),
                                                 float(np.max(wbc[1:, 3]))]
    # B 两井
    for spec, wdir in ((HT1003, "ht1_003_呼1-003"), (HU101, "hu101_呼101")):
        res_g = part_generic(spec, wdir)
        summary[spec["well"]] = dict(res_g["summary"])
        summary[spec["well"]]["nodes_min"] = [round(x, 3) for x in res_g["nodes_min"]]
        summary[spec["well"]]["fluid_names"] = res_g["fluid_names"]
    # C 稠化窗
    summary["d7a_ht1004"] = part_d7a(res_a)
    # D 作业史 notes 泵压
    for spec, wdir in ((HT1003, "ht1_003_呼1-003"), (HU101, "hu101_呼101"),
                       (None, "ht1_004_呼1-004")):
        well = wdir if spec is None else spec["well"]
        rows, lo, hi = part_notes(wdir)
        write_series_csv(OUT / ("作业史对账_notes泵压_%s.csv" % well),
                         ["step_index", "stage_name", "fluid_name", "volume_m3",
                          "rate_m3_min", "start_time_min", "end_time_min",
                          "density_g_cm3", "notes_pump_MPa_list", "notes"], rows)
        summary.setdefault("notes_pump_MPa", {})[well] = {"min": lo, "max": hi}
    with io.open(str(OUT / "run_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2, default=float)
    print(json.dumps(summary["ht1004"], ensure_ascii=False, indent=1, default=float)[:800])
    for k in ("ht1003_呼1-003", "hu101_呼101"):
        print(k, json.dumps(summary[k], ensure_ascii=False, default=float))
    print("D7A:", json.dumps({kk: (round(v, 2) if isinstance(v, float) else v)
                              for kk, v in summary["d7a_ht1004"].items()},
                             ensure_ascii=False, default=str))
    print("NOTES:", json.dumps(summary["notes_pump_MPa"], ensure_ascii=False))


if __name__ == "__main__":
    main()
