"""
反对照实验：把呼1-004 的居中度从 0.83 降到 0.60，看顶替效率是否暴跌。

逻辑：若居中度是首要驱动因素，则 0.83 -> 0.60（跌破 0.67 阈值）应导致效率大幅下降；
若效率几乎不变，则说明居中度不是主因，前面"居中度主导"的结论需要修正。

不修改任何模型源代码，仅以 dataclasses.replace 在内存中替换 well_spec.standoff_profile。
结果写到 results/_test_standoff/，不覆盖原始结果。
"""

import csv, json, os, sys
from dataclasses import replace
from pathlib import Path

import matplotlib; matplotlib.use("Agg")

ROOT = Path(__file__).resolve().parents[1]
os.chdir(ROOT)
sys.path.insert(0, str(ROOT))

from cemdisp.data.well_spec import DepthValuePoint
from cemdisp.transport1d import CasingFlowSolver
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.runners.ht1_004_tailpipe import annulus_stop_time_s as ht1_004_stop
from cemdisp.data.loaders.ht1_004_loader import load_ht1_004_tailpipe


def _const_standoff(well_spec, value):
    """用全井段恒定居中度替换 standoff_profile。"""
    top, bot = well_spec.top_md_m, well_spec.bottom_md_m
    prof = (
        DepthValuePoint(depth_md_m=top, value=value),
        DepthValuePoint(depth_md_m=bot, value=value),
    )
    return replace(well_spec, standoff_profile=prof)


def _profile_rows(result):
    """从 solver result.depth_profiles (DataFrame) 取行 dict。"""
    dp = result.depth_profiles
    return dp.to_dict(orient="records")


def _stats(rows):
    keys = list(rows[0].keys()) if rows else []
    def col(name_pat):
        k = next((k for k in keys if all(t in k for t in name_pat)), None)
        if not k:
            return None
        v = [float(r[k]) for r in rows]
        return {"mean": sum(v) / len(v), "min": min(v), "max": max(v)}
    return {
        "平均效率": col(["平均", "顶替"]),
        "窄边效率": col(["窄边", "有效"]),
        "宽边效率": col(["宽边", "有效"]),
        "居中度": col(["居中度"]),
    }


def run_with_standoff(standoff_val, out_dir, tag):
    """用指定恒定居中度运行呼1-004 全流程，返回摘要+深度剖面统计。"""
    print(f"\n===== 呼1-004  居中度={standoff_val*100:.1f}%  ({tag}) =====")
    well_spec, fluids, schedule, _ = load_ht1_004_tailpipe()
    spec = _const_standoff(well_spec, standoff_val)

    casing_solver = CasingFlowSolver(enable_gravity=True)
    casing_result = casing_solver.run(spec, fluids, schedule)
    coupled_provider = build_coupled_annulus_inlet_provider(
        casing_result, casing_solver, fluids, split_cement_phases=True
    )
    stop_s = ht1_004_stop(casing_result=casing_result, fluids=fluids)
    solver = AnnulusD2DGASolver(total_t=stop_s)
    result = solver.run(spec, fluids, coupled_provider)

    final = result.summary.get("最终结果", {})
    overall = final.get("全井段最终有效顶替效率", float("nan"))
    occ = final.get("最终水泥浆占据率", float("nan"))
    print(f"  全井段最终有效顶替效率 = {overall:.4f}  ({overall*100:.2f}%)")
    print(f"  最终水泥浆占据率       = {occ:.4f}")

    rows = _profile_rows(result)
    st = _stats(rows)
    for label, s in st.items():
        if s:
            print(f"  {label}: mean={s['mean']:.4f} min={s['min']:.4f} max={s['max']:.4f}")

    # 保存深度剖面
    keys = list(rows[0].keys())
    out_csv = out_dir / f"呼1-004_{tag}_standoff{standoff_val*100:.0f}_深度剖面.csv"
    with open(out_csv, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"  -> {out_csv}")

    return {"overall": overall, "occupancy": occ, "stats": st}


def read_original_csv(csv_rel):
    p = ROOT / csv_rel
    if not p.exists():
        return None
    with open(p, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    return _stats(rows), len(rows), p


if __name__ == "__main__":
    out_dir = ROOT / "results" / "_test_standoff"
    out_dir.mkdir(parents=True, exist_ok=True)

    summary = {}

    # 原始：0.83（读已有结果 CSV，避免重算）
    orig_csv = "results/呼1-004_1D2D耦合模型/优化参数/呼1-004_1D2D耦合模型_深度剖面.csv"
    orig_stats, orig_n, orig_path = read_original_csv(orig_csv)
    print(f"\n##### 原始(0.83) 读自: {orig_path}  (n={orig_n}) #####")
    if orig_stats:
        for label, s in orig_stats.items():
            if s:
                print(f"  {label}: mean={s['mean']:.4f} min={s['min']:.4f} max={s['max']:.4f}")
    summary["原始_0.83"] = {"source": str(orig_path), "stats": orig_stats}

    # 反对照：降到 0.60
    r60 = run_with_standoff(0.60, out_dir, "降到060")
    summary["降到0.60"] = r60

    # 顺带：降到 0.70（阈值附近），观察是否单调
    r70 = run_with_standoff(0.70, out_dir, "降到070")
    summary["降到0.70"] = r70

    # 汇总
    print("\n" + "=" * 64)
    print("呼1-004 居中度反对照汇总")
    print("=" * 64)
    print(f"{'居中度':>8} | {'全井加权':>8} | {'剖面均值':>8} | {'窄边均值':>8} | {'宽边均值':>8}")
    print("-" * 64)

    def line(tag, overall, stats):
        avg = stats["平均效率"]["mean"] if stats.get("平均效率") else float("nan")
        nar = stats["窄边效率"]["mean"] if stats.get("窄边效率") else float("nan")
        wid = stats["宽边效率"]["mean"] if stats.get("宽边效率") else float("nan")
        print(f"{tag:>8} | {overall*100:7.2f}% | {avg*100:7.2f}% | {nar*100:7.2f}% | {wid*100:7.2f}%")

    s0 = summary["原始_0.83"]["stats"]
    line("0.83原", 0.0 if not s0 else 0.0, s0)  # 原始 CSV 无 summary overall，用占位
    # 原始用摘要 JSON 的加权值
    import json as _json
    orig_json = ROOT / "results/呼1-004_1D2D耦合模型/优化参数/呼1-004_1D2D耦合模型_结果摘要.json"
    orig_overall = None
    if orig_json.exists():
        d = _json.loads(orig_json.read_text(encoding="utf-8"))
        orig_overall = d.get("最终结果", {}).get("全井段最终有效顶替效率")
    if orig_overall is not None:
        line("0.83原", orig_overall, s0)
    else:
        line("0.83原", float("nan"), s0)
    line("0.70", summary["降到0.70"]["overall"], summary["降到0.70"]["stats"])
    line("0.60", summary["降到0.60"]["overall"], summary["降到0.60"]["stats"])

    (out_dir / "呼1-004_反对照汇总.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, default=float),
        encoding="utf-8",
    )
    print(f"\n结果已写到: {out_dir}")
