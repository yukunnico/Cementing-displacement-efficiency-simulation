"""正文补齐变体实验（runner 口径）—— 2026-09-10 用户裁定"补齐剩余文章内容"。

四组变体（全部不改泵序/排量，施工时长不变，无稠化时间窗风险）：
  A. R 阶梯消融（8 井）：r0_base(I3+浮力关) / r2_i3(仅浮力关)；R3=生产基线(用现成 JSON)
     R 阶梯定义源 scripts/entrypoints/closure_contribution_scan.py:22-30；auto_m 已恒开(d8917b7)故 R1 并入
  B. 固定 dt 稳健性（8 井）：dt_fixed(enable_cfl_adaptive=False, dt=4.0) vs CFL 自适应基线
  C. standoff 连续扫描（呼103）：±0.05/±0.15/±0.2/±0.3 七新档（±0.1 已在敏感性目录）
  D. 隔离液密度大范围（呼103）：±0.2/±0.3 g/cc（±0.1 已在敏感性目录）
基线 = results/<井名>_1D2D耦合模型/*_结果摘要.json（唯一口径）。
输出 = results/正文补齐变体_runner口径_2026-09-10/。
"""

from __future__ import annotations

import csv
import json
import time
from dataclasses import replace
from importlib import import_module
from pathlib import Path

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.well_spec import WellSpec
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = PROJECT_ROOT / "results" / "正文补齐变体_runner口径_2026-09-10"

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

CEMENT_ROLES = {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}
T1_CASING_KW = dict(
    enable_gravity=True,
    mixing_contact_time=True,
    plug_face_zero_mixing=True,
    has_plug=True,
)


def _stop_time_s(casing_result, fluids) -> float:
    """停止时刻：cement_end_time_s 优先（F2 口径≡碰压），前缘扫描回退+排序防御（hu103 runner 同构）。"""
    if casing_result.cement_end_time_s is not None:
        return float(casing_result.cement_end_time_s)
    fluid_by_name = {f.name: f for f in fluids}
    fronts = sorted(casing_result.fronts, key=lambda f: f.time_s)
    last_cement = None
    for front in fronts:
        fluid = fluid_by_name.get(front.fluid_name)
        if fluid is not None and fluid.role in CEMENT_ROLES:
            last_cement = front.time_s
    if last_cement is not None:
        for front in fronts:
            fluid = fluid_by_name.get(front.fluid_name)
            if fluid is None or fluid.role in CEMENT_ROLES:
                continue
            if front.time_s >= last_cement - 1.0e-9:
                return float(front.time_s)
    return float(casing_result.cement_end_time_s)


def _read_baseline(well_dir: str) -> dict[str, float]:
    import glob
    hits = sorted(glob.glob(str(PROJECT_ROOT / "results" / well_dir / "*结果摘要.json")))
    hits = [h for h in hits if "初版" not in h]
    d = json.load(open(hits[0], encoding="utf-8"))
    return dict(d["最终结果"])


def shift_standoff(well: WellSpec, delta: float) -> WellSpec:
    """standoff 剖面整体平移 delta（clip [0,1]）。"""
    import numpy as np
    from cemdisp.data.well_spec import DepthValuePoint

    shifted = tuple(
        DepthValuePoint(p.depth_md_m, float(np.clip(p.value + delta, 0.0, 1.0)))
        for p in well.standoff_profile
    )
    return replace(well, standoff_profile=shifted)


def shift_spacer_density(fluids: tuple[FluidSpec, ...], delta_kg_m3: float) -> tuple[FluidSpec, ...]:
    """隔离液密度平移（仅 SPACER 角色）。"""
    return tuple(
        replace(f, density_kg_m3=f.density_kg_m3 + delta_kg_m3)
        if f.role == FluidRole.SPACER else f
        for f in fluids
    )


def run_case(well_key: str, tag: str, *, solver_kw: dict | None = None,
             well_fn=None, fluid_fn=None, out_dir: Path | None = None) -> dict:
    """单 case：load → (输入变换) → 1D(T1) → 耦合入口 → 停止时刻 → 环空(参数覆盖) → summary。"""
    module_name, func_name, well_dir = WELLS[well_key]
    loader = getattr(import_module(module_name), func_name)
    t0 = time.perf_counter()
    well, fluids, schedule, _ = loader()
    if well_fn is not None:
        well = well_fn(well)
    if fluid_fn is not None:
        fluids = fluid_fn(fluids)
    casing_solver = CasingFlowSolver(**T1_CASING_KW)
    casing_result = casing_solver.run(well, fluids, schedule)
    provider = build_coupled_annulus_inlet_provider(
        casing_result, casing_solver, fluids, split_cement_phases=True,
    )
    total_t_s = _stop_time_s(casing_result, fluids)
    kw = dict(total_t=total_t_s, nz=250)
    if solver_kw:
        kw.update(solver_kw)
    result = AnnulusD2DGASolver(**kw).run(well, fluids, provider, schedule=schedule)
    elapsed_s = time.perf_counter() - t0
    final = result.summary["最终结果"]
    row = {
        "井名": well_key,
        "变体": tag,
        "η_E": float(final["全井段最终有效顶替效率"]),
        "η_N": float(final["窄四分位效率"]),
        "窜槽": float(final["最终窜槽指数"]),
        "混浆": float(final["最终混浆指数"]),
        "失稳": float(final["最终失稳指数"]),
        "耗时_s": round(elapsed_s, 1),
    }
    assert out_dir is not None
    json_path = out_dir / f"{well_key}_{tag}_结果摘要.json"
    json_path.write_text(
        json.dumps(result.summary, ensure_ascii=False, indent=2), encoding="utf-8",
    )
    print(
        f"[OK] {well_key} × {tag}: η_E={row['η_E']:.4f} η_N={row['η_N']:.4f} ({elapsed_s:.0f}s)",
        flush=True,
    )
    return row


def main() -> None:
    out_dir = OUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []

    # A+B：8 井 × {r0_base, r2_i3, dt_fixed}
    for well_key in WELLS:
        for tag, kw in (
            ("r0_base", {"enable_d2dga_i3_flux": False, "enable_true_buoyancy": False}),
            ("r2_i3", {"enable_true_buoyancy": False}),
            ("dt_fixed", {"enable_cfl_adaptive": False, "dt": 4.0}),
        ):
            rows.append(run_case(well_key, tag, solver_kw=kw, out_dir=out_dir))

    # C：呼103 standoff 连续档（±0.1 已在敏感性目录）
    for delta in (-0.3, -0.2, -0.15, -0.05, 0.05, 0.2, 0.3):
        tag = f"standoff_{delta:+.2f}"
        rows.append(run_case(
            "呼103", tag, well_fn=lambda w, d=delta: shift_standoff(w, d), out_dir=out_dir,
        ))

    # D：呼103 隔离液密度大范围（±0.1 已在敏感性目录）
    for delta in (-300.0, -200.0, 200.0, 300.0):
        tag = f"spacer_dens_{delta:+.0f}"
        rows.append(run_case(
            "呼103", tag, fluid_fn=lambda fs, d=delta: shift_spacer_density(fs, d), out_dir=out_dir,
        ))

    # 汇总（附基线对比列）
    for row in rows:
        baseline = _read_baseline(WELLS[row["井名"]][2])
        row["Δη_E_vs_基线"] = row["η_E"] - float(baseline["全井段最终有效顶替效率"])
        row["Δη_N_vs_基线"] = row["η_N"] - float(baseline["窄四分位效率"])
    csv_path = out_dir / "汇总表_正文补齐变体.csv"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    md_lines = [
        "# 正文补齐变体汇总（runner 口径，2026-09-10）",
        "",
        "A. R阶梯消融：r0_base=I3+浮力关 / r2_i3=仅浮力关（R3=生产基线；定义源 closure_contribution_scan.py，auto_m 恒开）",
        "B. dt_fixed=固定4s步长（CFL自适应对照）",
        "C. 呼103 standoff 连续档（±0.05~±0.3；±0.1 见敏感性目录）",
        "D. 呼103 隔离液密度 ±0.2/±0.3 g/cc（±0.1 见敏感性目录）",
        "全部变体不改泵序/排量，施工时长不变",
        "",
        "| 井名 | 变体 | η_E | η_N | Δη_E vs 基线 | Δη_N vs 基线 | 耗时_s |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in rows:
        md_lines.append(
            f"| {row['井名']} | {row['变体']} | {row['η_E']:.4f} | {row['η_N']:.4f} "
            f"| {row['Δη_E_vs_基线']:+.4f} | {row['Δη_N_vs_基线']:+.4f} | {row['耗时_s']} |"
        )
    (out_dir / "汇总表_正文补齐变体.md").write_text("\n".join(md_lines), encoding="utf-8")
    print(f"完成：{len(rows)} 个变体 → {out_dir}")


if __name__ == "__main__":
    main()
