"""
双因素对照实验：在呼1-004 上分离"居中度"与"流体体系（领/尾浆流变模型）"各自贡献。

核心问题：把呼1-004 同时降到 0.60 居中度 + 换成 Power-Law 浆体，
全井效率是否真的从 88% 掉到 70%？

为干净分离变量，做四档阶梯对照（呼1-004 井，其余全部不变）：

  A. 基准       : 0.83 居中 + 原始 Bingham 领/尾浆        （= 原始结果，95.96%）
  B. 只降居中度 : 0.60 居中 + 原始 Bingham 领/尾浆        （已测，88.64%）
  C. 只换浆体   : 0.83 居中 + Power-Law 领/尾浆           （隔离纯流体体系贡献）
  D. 双因素     : 0.60 居中 + Power-Law 领/尾浆           （用户问的点）

"换浆体"的边界（单一变量）：
  - 仅把【领浆】和【尾浆】的流变模型从 Bingham(PV/YP) 换成 Power-Law(n/K)，
    参数取自呼1-002 现场（hu2_loader.py）：
      领浆 POWER_LAW n=0.811, K=0.876  ；尾浆 n=0.886, K=0.453
  - 保留呼1-004 原有的：领/尾浆密度、体积、泵速、双隔离液、替浆序列、其余流体（钻井液/
    平衡液/隔离液/压塞液/替浆液等）全部不变。
  - 这样 C、D 相对 A、B 的唯一差异就是领/尾浆流变模型，可隔离"流体体系"贡献。

不修改任何模型源代码：well_spec 用 dataclasses.replace 改 standoff_profile；
fluids 是 tuple，重建一个替换领/尾浆 FluidSpec 的新 tuple。
结果写到 results/_test_standoff/，不覆盖原始结果。
"""

import csv, json, os, sys
from dataclasses import replace
from pathlib import Path

import matplotlib; matplotlib.use("Agg")

ROOT = Path(__file__).resolve().parents[1]
os.chdir(ROOT)
sys.path.insert(0, str(ROOT))

from cemdisp.data.fluid_spec import FluidSpec, FluidRole, RheologyModel
from cemdisp.data.well_spec import DepthValuePoint
from cemdisp.transport1d import CasingFlowSolver
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.runners.ht1_004_tailpipe import annulus_stop_time_s as ht1_004_stop
from cemdisp.data.loaders.ht1_004_loader import load_ht1_004_tailpipe

# ---- 呼1-002 现场领/尾浆 Power-Law 参数（用于替换呼1-004 的领/尾浆流变模型）----
# 来源: cemdisp/data/loaders/hu2_loader.py:96-101, 233-235
HU2_LEAD_N, HU2_LEAD_K = 0.811, 0.876
HU2_TAIL_N, HU2_TAIL_K = 0.886, 0.453


def _const_standoff(well_spec, value):
    top, bot = well_spec.top_md_m, well_spec.bottom_md_m
    prof = (DepthValuePoint(depth_md_m=top, value=value),
           DepthValuePoint(depth_md_m=bot, value=value))
    return replace(well_spec, standoff_profile=prof)


def _swap_lead_tail_to_power_law(fluids):
    """把领浆/尾浆的流变模型从 Bingham 换成 Power-Law（n/K 取自呼1-002），
    保留各自的 density、role、name。其余 fluid 不变。返回新的 tuple。"""
    new = []
    for f in fluids:
        if f.role == FluidRole.LEAD:
            new.append(replace(
                f, rheology_model=RheologyModel.POWER_LAW,
                plastic_viscosity_pa_s=None, yield_stress_pa=None,
                power_law_n=HU2_LEAD_N, consistency_k=HU2_LEAD_K))
        elif f.role == FluidRole.TAIL:
            new.append(replace(
                f, rheology_model=RheologyModel.POWER_LAW,
                plastic_viscosity_pa_s=None, yield_stress_pa=None,
                power_law_n=HU2_TAIL_N, consistency_k=HU2_TAIL_K))
        else:
            new.append(f)
    return tuple(new)


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


def run(tag, standoff_val, use_power_law, out_dir):
    print(f"\n===== {tag}  居中度={standoff_val*100:.0f}%  领尾浆={'PowerLaw' if use_power_law else 'Bingham'} =====")
    well_spec, fluids, schedule, _ = load_ht1_004_tailpipe()
    spec = _const_standoff(well_spec, standoff_val)
    flds = _swap_lead_tail_to_power_law(fluids) if use_power_law else fluids

    casing_solver = CasingFlowSolver(enable_gravity=True)
    casing_result = casing_solver.run(spec, flds, schedule)
    coupled_provider = build_coupled_annulus_inlet_provider(
        casing_result, casing_solver, flds, split_cement_phases=True)
    stop_s = ht1_004_stop(casing_result=casing_result, fluids=flds)
    solver = AnnulusD2DGASolver(total_t=stop_s)
    result = solver.run(spec, flds, coupled_provider)

    final = result.summary.get("最终结果", {})
    overall = final.get("全井段最终有效顶替效率", float("nan"))
    occ = final.get("最终水泥浆占据率", float("nan"))
    print(f"  全井加权 = {overall:.4f}  ({overall*100:.2f}%)   占据率 = {occ:.4f}")
    st = _stats(_profile_rows(result))
    for label, s in st.items():
        if s:
            print(f"  {label}: mean={s['mean']:.4f} min={s['min']:.4f} max={s['max']:.4f}")

    # 保存深度剖面
    rows = _profile_rows(result)
    keys = list(rows[0].keys())
    out_csv = out_dir / f"呼1-004_{tag}_深度剖面.csv"
    with open(out_csv, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
        for r in rows: w.writerow(r)
    print(f"  -> {out_csv}")
    return {"overall": overall, "occupancy": occ, "stats": st}


def _orig():
    """读原始 0.83 结果（避免重算）。"""
    p = ROOT / "results/呼1-004_1D2D耦合模型/优化参数/呼1-004_1D2D耦合模型_深度剖面.csv"
    with open(p, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    j = ROOT / "results/呼1-004_1D2D耦合模型/优化参数/呼1-004_1D2D耦合模型_结果摘要.json"
    d = json.loads(j.read_text(encoding="utf-8")) if j.exists() else {}
    return {"overall": d.get("最终结果", {}).get("全井段最终有效顶替效率", float("nan")),
            "stats": _stats(rows)}


if __name__ == "__main__":
    out_dir = ROOT / "results" / "_test_standoff"
    out_dir.mkdir(parents=True, exist_ok=True)

    A = _orig()                                      # 基准 0.83 + Bingham
    B = run("B_只降居中度", 0.60, False, out_dir)   # 0.60 + Bingham（重测，便于同行打印）
    C = run("C_只换浆体",   0.83, True,  out_dir)   # 0.83 + PowerLaw
    D = run("D_双因素",     0.60, True,  out_dir)   # 0.60 + PowerLaw  ← 用户问的点

    # 汇总
    print("\n" + "=" * 76)
    print("呼1-004 双因素阶梯对照汇总")
    print("=" * 76)
    print(f"{'档':<22} | {'居中度':>6} | {'领尾浆':>8} | {'全井加权':>8} | {'剖面均值':>8} | {'窄边均值':>8} | {'宽边均值':>8}")
    print("-" * 76)
    def row(tag, standoff, rheo, d):
        st = d["stats"]; o = d["overall"]
        avg = st["平均效率"]["mean"] if st.get("平均效率") else float("nan")
        nar = st["窄边效率"]["mean"] if st.get("窄边效率") else float("nan")
        wid = st["宽边效率"]["mean"] if st.get("宽边效率") else float("nan")
        print(f"{tag:<22} | {standoff*100:5.0f}% | {rheo:>8} | {o*100:7.2f}% | {avg*100:7.2f}% | {nar*100:7.2f}% | {wid*100:7.2f}%")
    row("A 基准",            0.83, "Bingham",  A)
    row("B 只降居中度",       0.60, "Bingham",  B)
    row("C 只换浆体(PL)",     0.83, "PowerLaw", C)
    row("D 双因素(降+PL)",    0.60, "PowerLaw", D)

    print("\n边际贡献（相对基准A的全井加权变化）:")
    print(f"  纯居中度 B-A  = {(B['overall']-A['overall'])*100:+.2f} pp")
    print(f"  纯浆体   C-A  = {(C['overall']-A['overall'])*100:+.2f} pp")
    print(f"  双因素   D-A  = {(D['overall']-A['overall'])*100:+.2f} pp")
    print(f"  叠加预期 B+C-A = {(B['overall']+C['overall']-A['overall'])*100:+.2f} pp"
          f"  (若 ≈ D，则两因素近似独立可加)")

    json.dump({"A_基准": A, "B_只降居中度": B, "C_只换浆体": C, "D_双因素": D},
              open(out_dir / "呼1-004_双因素对照.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=2, default=float)
    print(f"\n结果已写到: {out_dir}")
