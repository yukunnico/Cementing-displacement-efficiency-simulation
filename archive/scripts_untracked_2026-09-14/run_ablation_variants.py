"""runner 单变量变体实验（消融组）—— 2026-09-09 用户裁定执行。

不改 cemdisp/runners/ 任何代码：loader 无参调用 + CasingFlowSolver(T1 生产参数)
+ AnnulusD2DGASolver 构造参数覆盖实现变体，链路与 8 井 runner 同构。

变体：i3_localized(enable_local_i3=True) / m2_regime_split(enable_regime_split=True)
     / dispersion_zero(弥散两系数=0)
基线 = results/<井名>_1D2D耦合模型/*_结果摘要.json（权威口径，不重跑）。
输出 = results/消融变体_runner口径_2026-09-09/（CSV/MD 汇总 + 每变体摘要 JSON）。
"""

from __future__ import annotations

import csv
import json
import time
from importlib import import_module
from pathlib import Path

from cemdisp.data.fluid_spec import FluidRole
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_ROOT / "results" / "消融变体_runner口径_2026-09-09"

WELLS: dict[str, tuple[str, str, str]] = {
    "呼探1": ("cemdisp.data.loaders.hu1_loader", "load_hu1_tailpipe", "呼探1尾管_1D2D耦合模型"),
    "呼101": ("cemdisp.data.loaders.hu101_loader", "load_hu101_tailpipe", "呼101尾管_1D2D耦合模型"),
    "呼102": ("cemdisp.data.loaders.hu102_loader", "load_hu102_tailpipe", "呼102尾管_1D2D耦合模型"),
    "呼103": ("cemdisp.data.loaders.hu103_loader", "load_hu103_tailpipe", "呼103尾管_1D2D耦合模型"),
    "呼探1-001": ("cemdisp.data.loaders.ht1_001_loader", "load_ht1_001_tailpipe", "呼探1-001尾管_1D2D耦合模型"),
    "呼探1-002": ("cemdisp.data.loaders.hu2_loader", "load_hu2_tailpipe", "呼探1-002尾管_1D2D耦合模型"),
    "呼1-003": ("cemdisp.data.loaders.ht1_003_loader", "load_ht1_003_tailpipe", "呼1-003_1D2D耦合模型"),
    "呼1-004": ("cemdisp.data.loaders.ht1_004_loader", "load_ht1_004_tailpipe", "呼1-004_1D2D耦合模型"),
}

VARIANTS: dict[str, dict[str, float | bool]] = {
    "i3_localized": {"enable_local_i3": True},
    "m2_regime_split": {"enable_regime_split": True},
    "dispersion_zero": {"dispersion_axial": 0.0, "dispersion_azimuthal": 0.0},
}


def _read_baseline(well_dir: str) -> dict[str, float]:
    """读权威基线摘要（排除 hu101 的"初版模型"目录产物）。"""
    import glob

    hits = sorted(glob.glob(str(PROJECT_ROOT / "results" / well_dir / "*结果摘要.json")))
    hits = [h for h in hits if "初版" not in h]
    d = json.load(open(hits[0], encoding="utf-8"))
    return dict(d["最终结果"])


def _stop_time_s(casing_result, fluids) -> float:
    """停止时刻：cement_end_time_s 优先（F2 口径≡碰压断面），前缘扫描回退+排序防御。

    复制自 runners/hu103_tailpipe.py annulus_stop_time_s（8 井通用防御版）。
    """
    if casing_result.cement_end_time_s is not None:
        return float(casing_result.cement_end_time_s)
    cement_roles = {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}
    fluid_by_name = {f.name: f for f in fluids}
    sorted_fronts = sorted(casing_result.fronts, key=lambda f: f.time_s)
    last_cement_time_s: float | None = None
    for front in sorted_fronts:
        fluid = fluid_by_name.get(front.fluid_name)
        if fluid is not None and fluid.role in cement_roles:
            last_cement_time_s = front.time_s
    if last_cement_time_s is not None:
        for front in sorted_fronts:
            fluid = fluid_by_name.get(front.fluid_name)
            if fluid is None or fluid.role in cement_roles:
                continue
            if front.time_s >= last_cement_time_s - 1.0e-9:
                return float(front.time_s)
    return float(casing_result.cement_end_time_s)


def main() -> None:
    out_dir = PROJECT_ROOT / "results" / "消融变体_runner口径_2026-09-09"
    out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for well_key, (module_name, func_name, well_dir) in WELLS.items():
        loader = getattr(import_module(module_name), func_name)
        well_spec, fluids, schedule, _ = loader()
        casing_solver = CasingFlowSolver(
            enable_gravity=True,
            mixing_contact_time=True,
            plug_face_zero_mixing=True,
            has_plug=True,
        )
        casing_result = casing_solver.run(well_spec, fluids, schedule)
        provider = build_coupled_annulus_inlet_provider(
            casing_result, casing_solver, fluids, split_cement_phases=True,
        )
        total_t_s = _stop_time_s(casing_result, fluids)
        baseline = _read_baseline(well_dir)
        for variant_name, overrides in VARIANTS.items():
            t0 = time.perf_counter()
            solver = AnnulusD2DGASolver(total_t=total_t_s, nz=250, **overrides)
            result = solver.run(well_spec, fluids, provider, schedule=schedule)
            elapsed_s = time.perf_counter() - t0
            final = result.summary["最终结果"]
            base_e = float(baseline["全井段最终有效顶替效率"])
            base_n = float(baseline["窄四分位效率"])
            rows.append({
                "井名": well_key,
                "变体": variant_name,
                "η_E": float(final["全井段最终有效顶替效率"]),
                "η_N": float(final["窄四分位效率"]),
                "窜槽": float(final["最终窜槽指数"]),
                "混浆": float(final["最终混浆指数"]),
                "失稳": float(final["最终失稳指数"]),
                "浮力数_b": float(final.get("浮力数_b", 0.0)),
                "Δη_E_vs_基线": float(final["全井段最终有效顶替效率"]) - base_e,
                "Δη_N_vs_基线": float(final["窄四分位效率"]) - base_n,
                "耗时_s": round(elapsed_s, 1),
            })
            summary_json = out_dir / f"{well_key}_{variant_name}_结果摘要.json"
            summary_json.write_text(
                json.dumps(result.summary, ensure_ascii=False, indent=2), encoding="utf-8",
            )
            print(
                f"[OK] {well_key} × {variant_name}: "
                f"η_E={float(final['全井段最终有效顶替效率']):.4f} "
                f"η_N={float(final['窄四分位效率']):.4f} ({elapsed_s:.0f}s)",
                flush=True,
            )
    csv_path = out_dir / "汇总表_消融变体.csv"
    fieldnames = list(rows[0].keys())
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    md_lines = [
        "# 消融变体汇总（runner 口径单变量，2026-09-09）",
        "",
        "变体：i3_localized=I3通量局部化开 / m2_regime_split=M2流态修正开 / dispersion_zero=弥散置零",
        "基线 = results/<井名>_1D2D耦合模型（唯一口径，不重跑）",
        "",
        "| 井名 | 变体 | η_E | η_N | Δη_E vs 基线 | Δη_N vs 基线 | 耗时_s |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in rows:
        md_lines.append(
            f"| {row['井名']} | {row['变体']} | {row['η_E']:.4f} | {row['η_N']:.4f} "
            f"| {row['Δη_E_vs_基线']:+.4f} | {row['Δη_N_vs_基线']:+.4f} | {row['耗时_s']} |"
        )
    (out_dir / "汇总表_消融变体.md").write_text("\n".join(md_lines), encoding="utf-8")
    print(f"完成：{len(rows)} 个变体 → {out_dir}")


if __name__ == "__main__":
    main()
