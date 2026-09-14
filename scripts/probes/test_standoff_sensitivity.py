"""
测试：将呼探1-001、呼探1-002 的居中度改为现场设计文档第6.3节给出的设计模拟值，
对比改前/改后的顶替效率。

- 呼探1-001：文档居中度 80.4%（centralization_profile.csv 第2行）
- 呼探1-002：文档居中度 78%（斯伦贝谢模拟，centralization_profile.csv 第2行）

仅运行求解器管线，结果写到 results/_test_standoff/，不覆盖原始结果。
"""

import csv, json, os, sys
from dataclasses import replace
from pathlib import Path

import matplotlib; matplotlib.use("Agg")  # 防止部分导入链弹窗
import matplotlib.font_manager as fm

ROOT = Path(__file__).resolve().parents[1]
os.chdir(ROOT)
sys.path.insert(0, str(ROOT))

from cemdisp.data.well_spec import DepthValuePoint
from cemdisp.transport1d import CasingFlowSolver
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.runners.ht1_001_tailpipe import annulus_stop_time_s as ht1_001_stop
from cemdisp.runners.hu2_tailpipe import annulus_stop_time_s as hu2_stop
from cemdisp.data.loaders.ht1_001_loader import load_ht1_001_tailpipe
from cemdisp.data.loaders.hu2_loader import load_hu2_tailpipe


def _const_standoff(well_spec, value):
    """用全井段恒定居中度替换 standoff_profile。"""
    top, bot = well_spec.top_md_m, well_spec.bottom_md_m
    prof = (
        DepthValuePoint(depth_md_m=top, value=value),
        DepthValuePoint(depth_md_m=bot, value=value),
    )
    return replace(well_spec, standoff_profile=prof)


def _depth_profile_rows(result):
    """从 solver result 取深度剖面行（dict）。"""
    # result.depth_profiles 暴露 to_csv；优先用 pandas 行；否则用底层列
    dp = result.depth_profiles
    # 兼容两种形态：DataFrame-like 或列字典
    if hasattr(dp, "iterrows"):
        return [dict(r) for _, r in dp.iterrows()]
    if hasattr(dp, "to_dict"):
        d = dp.to_dict(orient="records") if hasattr(dp, "columns") else None
        if d:
            return d
    # 列字典形态
    keys = list(dp.keys()) if isinstance(dp, dict) else []
    n = len(dp[keys[0]]) if keys else 0
    return [{k: dp[k][i] for k in keys} for i in range(n)]


def run_well(name, loader_fn, stop_fn, standoff_val, out_dir):
    print(f"\n===== {name}  居中度={standoff_val*100:.1f}% =====")
    well_spec, fluids, schedule, _ = loader_fn()
    spec_patched = _const_standoff(well_spec, standoff_val)

    casing_solver = CasingFlowSolver(enable_gravity=True)
    casing_result = casing_solver.run(spec_patched, fluids, schedule)
    coupled_provider = build_coupled_annulus_inlet_provider(
        casing_result, casing_solver, fluids, split_cement_phases=True
    )
    stop_s = stop_fn(casing_result=casing_result, fluids=fluids)
    solver = AnnulusD2DGASolver(total_t=stop_s)
    result = solver.run(spec_patched, fluids, coupled_provider)

    # 摘要
    final = result.summary.get("最终结果", {})
    overall = final.get("全井段最终有效顶替效率", float("nan"))
    occ = final.get("最终水泥浆占据率", float("nan"))
    print(f"  全井段最终有效顶替效率 = {overall:.4f}  ({overall*100:.2f}%)")
    print(f"  最终水泥浆占据率       = {occ:.4f}")

    # 深度剖面统计：平均/窄边/宽边
    rows = _depth_profile_rows(result)
    if rows:
        keys = list(rows[0].keys())
        avg_key = next((k for k in keys if "平均" in k and "顶替" in k), None)
        narrow_key = next((k for k in keys if "窄边" in k and "有效" in k), None)
        wide_key = next((k for k in keys if "宽边" in k and "有效" in k), None)
        cen_key = next((k for k in keys if "居中度" in k), None)
        if avg_key:
            vals = [float(r[avg_key]) for r in rows]
            print(f"  深度剖面 平均效率 mean={sum(vals)/len(vals):.4f} "
                  f"min={min(vals):.4f} max={max(vals):.4f}")
        if narrow_key:
            vals = [float(r[narrow_key]) for r in rows]
            print(f"  窄边效率 mean={sum(vals)/len(vals):.4f} "
                  f"min={min(vals):.4f}")
        if wide_key:
            vals = [float(r[wide_key]) for r in rows]
            print(f"  宽边效率 mean={sum(vals)/len(vals):.4f}")
        if cen_key:
            vals = [float(r[cen_key]) for r in rows]
            print(f"  居中度(校验) mean={sum(vals)/len(vals):.4f}")

        # 保存深度剖面 CSV 便于核对
        out_csv = out_dir / f"{name}_standoff{standoff_val*100:.1f}_深度剖面.csv"
        with open(out_csv, "w", encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in rows:
                w.writerow(r)
        print(f"  -> {out_csv}")

    return {"overall": overall, "occupancy": occ}


def read_original(name, csv_rel):
    """读原始结果 CSV 的统计。"""
    p = ROOT / csv_rel
    if not p.exists():
        return None
    with open(p, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    avg = [float(r["平均有效顶替效率"]) for r in rows]
    narrow = [float(r["窄边有效效率"]) for r in rows]
    wide = [float(r["宽边有效效率"]) for r in rows]
    return {
        "csv_mean": sum(avg) / len(avg),
        "narrow_mean": sum(narrow) / len(narrow),
        "wide_mean": sum(wide) / len(wide),
        "n": len(rows),
    }


if __name__ == "__main__":
    out_dir = ROOT / "results" / "_test_standoff"
    out_dir.mkdir(parents=True, exist_ok=True)

    results = {}

    # 呼探1-001：文档居中度 80.4%
    orig_001 = read_original("呼探1-001",
                             "results/呼探1-001尾管_1D2D耦合模型/呼探1-001尾管_1D2D耦合模型_深度剖面.csv")
    new_001 = run_well("呼探1-001", load_ht1_001_tailpipe, ht1_001_stop, 0.804, out_dir)
    results["呼探1-001"] = {"original": orig_001, "patched_0.804": new_001}

    # 呼探1-002：文档居中度 78%
    orig_002 = read_original("呼探1-002",
                             "results/呼探1-002尾管_1D2D耦合模型/呼探1-002尾管_1D2D耦合模型_深度剖面.csv")
    new_002 = run_well("呼探1-002", load_hu2_tailpipe, hu2_stop, 0.78, out_dir)
    results["呼探1-002"] = {"original": orig_002, "patched_0.78": new_002}

    # 汇总
    print("\n" + "=" * 60)
    print("汇总对比（改前 loader代理值 -> 改后 文档设计模拟值）")
    print("=" * 60)
    for well, d in results.items():
        o = d["original"]
        n = d["patched_0.804"] if "patched_0.804" in d else d["patched_0.78"]
        if o:
            print(f"\n{well}:")
            print(f"  改前(loader代理): CSV均值={o['csv_mean']:.4f}  "
                  f"窄边均值={o['narrow_mean']:.4f}  宽边均值={o['wide_mean']:.4f}")
            print(f"  改后(文档居中度): 全井加权={n['overall']:.4f}  "
                  f"占据率={n['occupancy']:.4f}")
        else:
            print(f"\n{well}: 原始CSV未找到")

    # 写 JSON
    (out_dir / "对比汇总.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2, default=float),
        encoding="utf-8",
    )
    print(f"\n结果已写到: {out_dir}")
