"""
呼1-004 几何/泵速影响因素分析 + 呼1-002 反向对照。

第一部分（呼1-002 反向对照）：把呼1-002 居中度提到 0.83 + 领/尾浆换 Bingham，
看全井效率能否从 70% 跳到 90%。
  - 居中度 0.58–0.75 -> 均匀 0.83（dataclasses.replace）
  - 领/尾浆 Power-Law -> Bingham（取呼1-004 的领/尾 PV/YP，保留 002 的密度）
  边界：仅改居中度与领/尾浆流变，002 原有的几何(26.7mm窄间隙)/泵速/体积/流体序列不变。

第二部分（呼1-004 几何/泵速）：在呼1-004 基准(0.83+Bingham)上分别/叠加施加：
  - 几何：缩小等效井径 -> 模拟呼1-002 的窄环空（间隙 38mm -> ~26mm）
  - 泵速：等比降低所有 PumpingScheduleStep.rate_m3_min -> 模拟呼1-002 的低泵速
  目的：隔离几何与泵速对呼1-004 高效率的贡献，看能否从 96% 掉到 70%。

不修改模型源代码。结果写到 results/_test_standoff/，不覆盖原始结果。
"""

import csv, json, math, os, sys
from dataclasses import replace
from pathlib import Path

import matplotlib; matplotlib.use("Agg")

ROOT = Path(__file__).resolve().parents[1]
os.chdir(ROOT)
sys.path.insert(0, str(ROOT))

from cemdisp.data.fluid_spec import FluidSpec, FluidRole, RheologyModel
from cemdisp.data.well_spec import DepthValuePoint
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep
from cemdisp.transport1d import CasingFlowSolver
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.runners.ht1_004_tailpipe import annulus_stop_time_s as s004
from cemdisp.runners.hu2_tailpipe import annulus_stop_time_s as s002
from cemdisp.data.loaders.ht1_004_loader import load_ht1_004_tailpipe
from cemdisp.data.loaders.hu2_loader import load_hu2_tailpipe

# 呼1-004 领/尾浆 Bingham 参数（用于给 002 换 Bingham，取自 ht1_004_loader）
HT004_LEAD_PV, HT004_LEAD_YP = 0.170, 13.0
HT004_TAIL_PV, HT004_TAIL_YP = 0.180, 14.0
HT1_004_LOWER_LINER_OD_MM = 139.7


# ====================== 通用工具 ======================
def _const_standoff(well_spec, value):
    top, bot = well_spec.top_md_m, well_spec.bottom_md_m
    prof = (DepthValuePoint(depth_md_m=top, value=value),
            DepthValuePoint(depth_md_m=bot, value=value))
    return replace(well_spec, standoff_profile=prof)


def _swap_lead_tail_to_bingham(fluids, lead_pv, lead_yp, tail_pv, tail_yp):
    """领/尾浆 Power-Law -> Bingham(PV/YP)，保留各自 density/role/name。"""
    out = []
    for f in fluids:
        if f.role == FluidRole.LEAD:
            out.append(replace(f, rheology_model=RheologyModel.BINGHAM,
                               plastic_viscosity_pa_s=lead_pv, yield_stress_pa=lead_yp,
                               power_law_n=None, consistency_k=None))
        elif f.role == FluidRole.TAIL:
            out.append(replace(f, rheology_model=RheologyModel.BINGHAM,
                               plastic_viscosity_pa_s=tail_pv, yield_stress_pa=tail_yp,
                               power_law_n=None, consistency_k=None))
        else:
            out.append(f)
    return tuple(out)


def _narrow_annulus(well_spec, target_gap_mm, liner_od_mm=HT1_004_LOWER_LINER_OD_MM):
    """把等效井径整体缩小，使下段平均环空间隙 ~ target_gap_mm。

    实现：每个节点 hole_d -> liner_od + 2*( (hole_d - liner_od)/2 ) * (target/nominal_ratio)。
    为简洁起见，直接把 (hole_d - liner_od) 缩放到目标间隙比例。这里用一个简单做法：
    把每个节点 hole_d 减小，使径向间隙(hole_d-liner_od)/2 按比例缩小到 target_gap。
    但各点间隙不同，统一用一个缩放因子 r = target_gap/nominal_gap 更稳健。
    """
    # 先估当前下段平均间隙
    pts = [(p.depth_md_m, p.value) for p in well_spec.hole_diameter_profile]
    # 用下半段估 nominal gap
    zmid = (well_spec.top_md_m + well_spec.bottom_md_m) / 2.0
    lower = [d for z, d in pts if z >= zmid]
    if not lower:
        lower = [d for _, d in pts]
    nominal_gap = (sum(lower) / len(lower) - liner_od_mm) / 2.0  # 径向间隙 mm
    if nominal_gap <= 0:
        return well_spec
    r = max(target_gap_mm / nominal_gap, 0.1)  # 间隙缩放比
    new_pts = tuple(
        DepthValuePoint(depth_md_m=z, value=liner_od_mm + (d - liner_od_mm) * r)
        for z, d in pts
    )
    return replace(well_spec, hole_diameter_profile=new_pts)


def _scale_schedule_rates(schedule, factor):
    """所有 PumpingScheduleStep.rate_m3_min 乘 factor，体积不变。"""
    new_steps = tuple(
        replace(s, rate_m3_min=s.rate_m3_min * factor) if s.rate_m3_min else s
        for s in schedule.steps
    )
    return replace(schedule, steps=new_steps)


def _profile_rows(result):
    return result.depth_profiles.to_dict(orient="records")


def _stats(rows):
    if not rows:
        return {}
    keys = list(rows[0].keys())
    def col(*need):
        k = next((k for k in keys if all(t in k for t in need)), None)
        if not k:
            return None
        v = [float(r[k]) for r in rows]
        return {"mean": sum(v) / len(v), "min": min(v), "max": max(v)}
    return {"平均效率": col("平均", "顶替"),
            "窄边效率": col("窄边", "有效"),
            "宽边效率": col("宽边", "有效"),
            "居中度": col("居中度")}


def _run(label, loader_fn, stop_fn, standoff, use_bingham_for_002, geom_gap, rate_factor, out_dir):
    """通用运行。geom_gap/rate_factor 仅对 004 有效；002 仅用 standoff+bingham。"""
    print(f"\n===== {label} =====")
    well_spec, fluids, schedule, _ = loader_fn()

    # 002: 换居中度 + (可选)领尾浆换Bingham
    if standoff is not None:
        well_spec = _const_standoff(well_spec, standoff)
    if use_bingham_for_002:
        fluids = _swap_lead_tail_to_bingham(fluids, HT004_LEAD_PV, HT004_LEAD_YP,
                                             HT004_TAIL_PV, HT004_TAIL_YP)
    # 004: 几何 + 泵速
    if geom_gap is not None:
        well_spec = _narrow_annulus(well_spec, geom_gap)
    if rate_factor is not None:
        schedule = _scale_schedule_rates(schedule, rate_factor)

    casing_solver = CasingFlowSolver(enable_gravity=True)
    casing_result = casing_solver.run(well_spec, fluids, schedule)
    coupled_provider = build_coupled_annulus_inlet_provider(
        casing_result, casing_solver, fluids, split_cement_phases=True)
    stop_s = stop_fn(casing_result=casing_result, fluids=fluids)
    solver = AnnulusD2DGASolver(total_t=stop_s)
    result = solver.run(well_spec, fluids, coupled_provider)

    final = result.summary.get("最终结果", {})
    overall = final.get("全井段最终有效顶替效率", float("nan"))
    st = _stats(_profile_rows(result))
    print(f"  全井加权 = {overall*100:.2f}%")
    for lab, s in st.items():
        if s:
            print(f"  {lab}: mean={s['mean']:.4f} min={s['min']:.4f} max={s['max']:.4f}")
    # 存深度剖面
    rows = _profile_rows(result); keys = list(rows[0].keys())
    out_csv = out_dir / f"{label}_深度剖面.csv"
    with open(out_csv, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
        for r in rows: w.writerow(r)
    return {"overall": overall, "stats": st}


def _read_csv_stats(rel):
    p = ROOT / rel
    if not p.exists():
        return None
    with open(p, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    return _stats(rows), len(rows)


# ====================== 第一部分：呼1-002 反向对照 ======================
if __name__ == "__main__":
    out_dir = ROOT / "results" / "_test_standoff"
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = {}

    print("\n" + "#" * 70)
    print("# 第一部分：呼1-002 反向对照（0.83 + Bingham，能否 70% -> 90%）")
    print("#" * 70)
    # 002 原始
    o002_stats, o002_n = _read_csv_stats(
        "results/呼探1-002尾管_1D2D耦合模型/呼探1-002尾管_1D2D耦合模型_深度剖面.csv")
    o002_json = ROOT / "results/呼探1-002尾管_1D2D耦合模型/呼探1-002尾管_1D2D耦合模型_结果摘要.json"
    o002_overall = json.loads(o002_json.read_text(encoding="utf-8")).get("最终结果", {}).get(
        "全井段最终有效顶替效率") if o002_json.exists() else float("nan")
    print(f"\n002 原始(0.58-0.75+PowerLaw): 全井加权={o002_overall*100:.2f}%  (n={o002_n})")
    for lab, s in o002_stats.items():
        if s:
            print(f"  {lab}: mean={s['mean']:.4f}")
    summary["002_原始"] = {"overall": o002_overall, "stats": o002_stats}

    # 002 只提居中度
    r_002a = _run("002_A_只提居中度0.83", load_hu2_tailpipe, s002,
                  0.83, False, None, None, out_dir)
    summary["002_A_只提居中度0.83"] = r_002a
    # 002 只换Bingham浆体
    r_002b = _run("002_B_只换Bingham浆体", load_hu2_tailpipe, s002,
                  None, True, None, None, out_dir)
    summary["002_B_只换Bingham浆体"] = r_002b
    # 002 双因素：0.83 + Bingham
    r_002c = _run("002_C_双因素0.83+Bingham", load_hu2_tailpipe, s002,
                  0.83, True, None, None, out_dir)
    summary["002_C_双因素0.83+Bingham"] = r_002c

    # ====================== 第二部分：呼1-004 几何/泵速 ======================
    print("\n" + "#" * 70)
    print("# 第二部分：呼1-004 几何/泵速影响因素（基准 0.83+Bingham=95.96%）")
    print("#" * 70)
    # 004 基准
    o004_json = ROOT / "results/呼1-004_1D2D耦合模型/优化参数/呼1-004_1D2D耦合模型_结果摘要.json"
    o004_overall = json.loads(o004_json.read_text(encoding="utf-8")).get("最终结果", {}).get(
        "全井段最终有效顶替效率") if o004_json.exists() else float("nan")
    o004_stats, o004_n = _read_csv_stats(
        "results/呼1-004_1D2D耦合模型/优化参数/呼1-004_1D2D耦合模型_深度剖面.csv")
    print(f"\n004 基准(0.83+Bingham): 全井加权={o004_overall*100:.2f}%  (n={o004_n})")
    summary["004_基准"] = {"overall": o004_overall, "stats": o004_stats}

    # 004 只缩几何（间隙 -> 26mm，模拟002窄环空）
    r_004g = _run("004_D_只缩几何26mm", load_ht1_004_tailpipe, s004,
                  None, False, 26.0, None, out_dir)
    summary["004_D_只缩几何26mm"] = r_004g
    # 004 只降泵速（因子0.5，模拟002的低泵速）
    r_004p = _run("004_E_只降泵速x0.5", load_ht1_004_tailpipe, s004,
                  None, False, None, 0.5, out_dir)
    summary["004_E_只降泵速x0.5"] = r_004p
    # 004 双因素：缩几何 + 降泵速
    r_004gp = _run("004_F_缩几何+降泵速", load_ht1_004_tailpipe, s004,
                  None, False, 26.0, 0.5, out_dir)
    summary["004_F_缩几何+降泵速"] = r_004gp
    # 004 三因素：缩几何 + 降泵速 + 降居中度0.60（看能否掉到70%）
    r_004all = _run("004_G_缩几何+降泵速+降居中0.60", load_ht1_004_tailpipe, s004,
                    0.60, False, 26.0, 0.5, out_dir)
    summary["004_G_缩几何+降泵速+降居中0.60"] = r_004all

    # ====================== 汇总 ======================
    print("\n" + "=" * 80)
    print("呼1-002 反向对照")
    print("=" * 80)
    print(f"{'档':<30} | {'居中度':>8} | {'领尾浆':>9} | {'全井加权':>8} | {'窄边均值':>8}")
    print("-" * 80)
    def line(tag, standoff, rheo, d):
        st = d["stats"]; o = d["overall"]
        nar = st["窄边效率"]["mean"] if st.get("窄边效率") else float("nan")
        print(f"{tag:<30} | {standoff:>8} | {rheo:>9} | {o*100:7.2f}% | {nar*100:7.2f}%")
    line("002 原始", "0.58-0.75", "PowerLaw", summary["002_原始"])
    line("002 A 只提居中度", "0.83", "PowerLaw", summary["002_A_只提居中度0.83"])
    line("002 B 只换Bingham", "0.58-0.75", "Bingham", summary["002_B_只换Bingham浆体"])
    line("002 C 双因素", "0.83", "Bingham", summary["002_C_双因素0.83+Bingham"])

    print("\n" + "=" * 80)
    print("呼1-004 几何/泵速影响（基准 95.96%）")
    print("=" * 80)
    print(f"{'档':<32} | {'几何间隙':>8} | {'泵速':>6} | {'居中度':>6} | {'全井加权':>8} | {'窄边均值':>8}")
    print("-" * 80)
    def line2(tag, gap, rate, cen, d):
        st = d["stats"]; o = d["overall"]
        nar = st["窄边效率"]["mean"] if st.get("窄边效率") else float("nan")
        print(f"{tag:<32} | {gap:>8} | {rate:>6} | {cen:>6} | {o*100:7.2f}% | {nar*100:7.2f}%")
    line2("004 基准", "38mm", "1x", "0.83", summary["004_基准"])
    line2("004 D 只缩几何", "26mm", "1x", "0.83", summary["004_D_只缩几何26mm"])
    line2("004 E 只降泵速", "38mm", "0.5x", "0.83", summary["004_E_只降泵速x0.5"])
    line2("004 F 缩几何+降泵速", "26mm", "0.5x", "0.83", summary["004_F_缩几何+降泵速"])
    line2("004 G 三因素(+降居中)", "26mm", "0.5x", "0.60", summary["004_G_缩几何+降泵速+降居中0.60"])

    json.dump(summary,
              open(out_dir / "几何泵速_双井对照.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=2, default=float)
    print(f"\n结果已写到: {out_dir}")
