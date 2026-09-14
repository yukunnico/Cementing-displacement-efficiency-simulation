"""runner 口径敏感性组实验（现场杠杆单变量）—— 2026-09-10 用户裁定执行。

与 run_ablation_variants.py 同骨架，差异：杠杆全在输入侧（well/fluids/schedule），
每变体 1D 重跑（排量/流体变化传导到时序与停止时刻）。链路与 8 井 runner 同构：
load → 输入变换 → CasingFlowSolver(T1) → 耦合入口 → 停止时刻 → 环空默认参数 run。

杠杆（档位与 09-07 矩阵一致，用于复现方向、产出唯一口径数字）：
    rate_x1.2 / rate_x0.8 / rate_x1.4 / rate_x0.6   排量缩放（泵序 steps 逐段 ×factor）
    cement_n_p0.1 / cement_n_m0.1                    水泥浆 n ±0.1（LEAD/INTERMEDIATE/TAIL）
    spacer_dens_p100 / spacer_dens_m100              隔离液密度 ±0.10 g/cc（SPACER ±100 kg/m³）
    mud_pv_x1.3 / mud_pv_x0.7                        泥浆 PV ±30%（MUD）
    standoff_p0.1 / standoff_m0.1                    standoff 剖面整体 ±0.1（clip [0,1]）
井：hu103（低偏心 e=0.278）+ hu101（强偏心 e=0.547）两极代表井。
基线 = results/<井名>_1D2D耦合模型/*_结果摘要.json（唯一口径，不新跑）。
输出 = results/敏感性变体_runner口径_2026-09-10/（CSV/MD 汇总 + 每变体摘要 JSON）。
"""

from __future__ import annotations

import csv
import json
import time
from dataclasses import replace
from importlib import import_module
from pathlib import Path

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.pumping_schedule import PumpingSchedule
from cemdisp.data.well_spec import WellSpec
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_ROOT / "results" / "敏感性变体_runner口径_2026-09-10"

WELLS: dict[str, tuple[str, str, str]] = {
    "呼103": ("cemdisp.data.loaders.hu103_loader", "load_hu103_tailpipe", "呼103尾管_1D2D耦合模型"),
    "呼101": ("cemdisp.data.loaders.hu101_loader", "load_hu101_tailpipe", "呼101尾管_1D2D耦合模型"),
}

CEMENT_ROLES = {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}


def scale_schedule(schedule: PumpingSchedule, factor: float) -> PumpingSchedule:
    """泵序逐段排量缩放（体积不变，时长随 1/factor 拉长）。"""
    new_steps = tuple(
        replace(s, rate_m3_min=s.rate_m3_min * factor) if s.rate_m3_min > 0 else s
        for s in schedule.steps
    )
    return replace(schedule, steps=new_steps)


def shift_cement_n(fluids: tuple[FluidSpec, ...], delta: float) -> tuple[FluidSpec, ...]:
    """水泥浆幂律指数 n ± delta（仅水泥三角色且 power_law_n 非空时生效）。"""
    return tuple(
        replace(f, power_law_n=f.power_law_n + delta)
        if f.role in CEMENT_ROLES and f.power_law_n is not None else f
        for f in fluids
    )


def shift_spacer_density(fluids: tuple[FluidSpec, ...], delta_kg_m3: float) -> tuple[FluidSpec, ...]:
    """隔离液密度 ± delta_kg_m3（仅 SPACER 角色）。"""
    return tuple(
        replace(f, density_kg_m3=f.density_kg_m3 + delta_kg_m3) if f.role == FluidRole.SPACER else f
        for f in fluids
    )


def scale_mud_pv(fluids: tuple[FluidSpec, ...], factor: float) -> tuple[FluidSpec, ...]:
    """泥浆塑性黏度 × factor（仅 MUD 角色）。"""
    return tuple(
        replace(f, plastic_viscosity_pa_s=f.plastic_viscosity_pa_s * factor)
        if f.role == FluidRole.MUD and f.plastic_viscosity_pa_s is not None else f
        for f in fluids
    )


def shift_standoff(well: WellSpec, delta: float) -> WellSpec:
    """standoff 剖面整体平移 delta（clip [0,1]），与 09-07 shift_standoff 同实现。"""
    from cemdisp.data.well_spec import DepthValuePoint
    import numpy as np

    if not well.standoff_profile:
        raise ValueError(f"{well.well_name}: 无 standoff_profile")
    shifted = tuple(
        DepthValuePoint(p.depth_md_m, float(np.clip(p.value + delta, 0.0, 1.0)))
        for p in well.standoff_profile
    )
    return replace(well, standoff_profile=shifted)


def _stop_time_s(casing_result, fluids) -> float:
    """停止时刻：cement_end_time_s 优先（F2 口径），前缘扫描回退+排序防御。"""
    if casing_result.cement_end_time_s is not None:
        return float(casing_result.cement_end_time_s)
    fluid_by_name = {f.name: f for f in fluids}
    sorted_fronts = sorted(casing_result.fronts, key=lambda f: f.time_s)
    last_cement_time_s = None
    for front in sorted_fronts:
        fluid = fluid_by_name.get(front.fluid_name)
        if fluid is not None and fluid.role in CEMENT_ROLES:
            last_cement_time_s = front.time_s
    if last_cement_time_s is not None:
        for front in sorted_fronts:
            fluid = fluid_by_name.get(front.fluid_name)
            if fluid is None or fluid.role in CEMENT_ROLES:
                continue
            if front.time_s >= last_cement_time_s - 1.0e-9:
                return float(front.time_s)
    return float(casing_result.cement_end_time_s)


def _read_baseline(well_dir: str) -> dict[str, float]:
    import glob
    hits = sorted(glob.glob(str(PROJECT_ROOT / "results" / well_dir / "*结果摘要.json"))
                  )
    hits = [h for h in hits if "初版" not in h]
    d = json.load(open(hits[0], encoding="utf-8"))
    return dict(d["最终结果"])


def _id(w):
    return w


def _idf(fs):
    return fs


def _ids(s):
    return s


def main() -> None:
    out_dir = PROJECT_ROOT / "results" / "敏感性变体_runner口径_2026-09-10"
    out_dir.mkdir(parents=True, exist_ok=True)
    variants: list[tuple[str, object, object, object]] = [
        ("rate_x1.2", _id, _idf, lambda s: scale_schedule(s, 1.2)),
        ("rate_x0.8", _id, _idf, lambda s: scale_schedule(s, 0.8)),
        ("rate_x1.4", _id, _idf, lambda s: scale_schedule(s, 1.4)),
        ("rate_x0.6", _id, _idf, lambda s: scale_schedule(s, 0.6)),
        ("cement_n_p0.1", _id, lambda fs: shift_cement_n(fs, +0.1), _ids),
        ("cement_n_m0.1", _id, lambda fs: shift_cement_n(fs, -0.1), _ids),
        ("spacer_dens_p100", _id, lambda fs: shift_spacer_density(fs, +100.0), _ids),
        ("spacer_dens_m100", _id, lambda fs: shift_spacer_density(fs, -100.0), _ids),
        ("mud_pv_x1.3", _id, lambda fs: scale_mud_pv(fs, 1.3), _ids),
        ("mud_pv_x0.7", _id, lambda fs: scale_mud_pv(fs, 0.7), _ids),
        ("standoff_p0.1", lambda w: shift_standoff(w, +0.1), _idf, _ids),
        ("standoff_m0.1", lambda w: shift_standoff(w, -0.1), _idf, _ids),
    ]
    rows: list[dict[str, object]] = []
    for well_key, (module_name, func_name, well_dir) in WELLS.items():
        loader = getattr(import_module(module_name), func_name)
        for var_name, well_fn, fluid_fn, sched_fn in variants:
            t0 = time.perf_counter()
            well, fluids, schedule, _ = loader()
            well2 = well_fn(well)
            fluids2 = fluid_fn(fluids)
            schedule2 = sched_fn(schedule)
            casing_solver = CasingFlowSolver(
                enable_gravity=True,
                mixing_contact_time=True,
                plug_face_zero_mixing=True,
                has_plug=True,
            )
            casing_result = casing_solver.run(well2, fluids2, schedule2)
            provider = build_coupled_annulus_inlet_provider(
                casing_result, casing_solver, fluids2, split_cement_phases=True,
            )
            total_t_s = _stop_time_s(casing_result, fluids2)
            solver = AnnulusD2DGASolver(total_t=total_t_s, nz=250)
            result = solver.run(well2, fluids2, provider, schedule=schedule2)
            elapsed_s = time.perf_counter() - t0
            final = result.summary["最终结果"]
            baseline = _read_baseline(well_dir)
            rows.append({
                "井名": well_key,
                "变体": var_name,
                "η_E": float(final["全井段最终有效顶替效率"]),
                "η_N": float(final["窄四分位效率"]),
                "窜槽": float(final["最终窜槽指数"]),
                "混浆": float(final["最终混浆指数"]),
                "失稳": float(final["最终失稳指数"]),
                "Δη_E_vs_基线": float(final["全井段最终有效顶替效率"]) - float(baseline["全井段最终有效顶替效率"]),
                "Δη_N_vs_基线": float(final["窄四分位效率"]) - float(baseline["窄四分位效率"]),
                "耗时_s": round(elapsed_s, 1),
            })
            json_path = out_dir / f"{well_key}_{var_name}_结果摘要.json"
            json_path.write_text(
                json.dumps(result.summary, ensure_ascii=False, indent=2), encoding="utf-8",
            )
            print(
                f"[OK] {well_key} × {var_name}: "
                f"η_E={float(final['全井段最终有效顶替效率']):.4f} "
                f"η_N={float(final['窄四分位效率']):.4f} ({elapsed_s:.0f}s)",
                flush=True,
            )
    csv_path = out_dir / "汇总表_敏感性变体.csv"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    md_lines = [
        "# 敏感性变体汇总（runner 口径单变量，2026-09-10）",
        "",
        "杠杆：排量±20/40% / 水泥浆n±0.1 / 隔离液密度±0.10g/cc / 泥浆PV±30% / standoff±0.1",
        "井：呼103（低偏心 e=0.278）+ 呼101（强偏心 e=0.547）两极代表井；每变体 1D 重跑",
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
    (out_dir / "汇总表_敏感性变体.md").write_text("\n".join(md_lines), encoding="utf-8")
    print(f"完成：{len(rows)} 个变体 → {out_dir}")


if __name__ == "__main__":
    main()
