"""HT1-004 R0->R3 改进 D2DGA 四级消融运行器。

R0: enable_d2dga_auto_m=False, enable_d2dga_i3_flux=False, enable_true_buoyancy=False
R1: +auto-m
R2: +I3 flux
R3: +true buoyancy
每级产出 AnnulusSimulationResult + 摘要，供论文图6-9。

注意：当前只有 enable_d2dga_auto_m 开关存在（Task 2）。i3_flux 和 true_buoyancy
开关由 Task 4/5 添加后会接入此运行器；在此之前 R2/R3 行为与 R1 相同。
"""
from __future__ import annotations

import csv
import json
import os
from dataclasses import dataclass
from typing import Dict, Optional
from pathlib import Path

from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver, AnnulusSimulationResult
from cemdisp.transport1d import CasingFlowSolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.runners.ht1_004_tailpipe import annulus_stop_time_s
from cemdisp.data.loaders.ht1_004_loader import load_ht1_004_tailpipe
from cemdisp.data.well_spec import WellSpec


@dataclass
class AblationLevel:
    """消融级别定义——三个开关的布尔组合。

    当前只有 enable_d2dga_auto_m 在求解器上存在（Task 2）。
    i3_flux 和 true_buoyancy 字段保留供 Task 4/5/6 前向使用。
    """
    name: str  # "R0".."R3"
    enable_d2dga_auto_m: bool
    enable_d2dga_i3_flux: bool
    enable_true_buoyancy: bool


ABLATION_LEVELS = [
    AblationLevel("R0", False, False, False),
    AblationLevel("R1", True,  False, False),
    AblationLevel("R2", True,  True,  False),
    AblationLevel("R3", True,  True,  True),
]

# ---------------------------------------------------------------------------
# Shared CSV / JSON helpers — used by ablation, theoretical, and convergence
# scripts to avoid copy-paste duplication.
# ---------------------------------------------------------------------------

ABLATION_CSV_COLUMNS = [
    "run_id", "ablation_level", "eccentricity", "nz", "dt",
    "effective_efficiency", "channeling_index", "mixing_index",
    "cement_occupation", "instability_index", "buoyancy_number",
]


def append_ablation_csv_row(
    row: dict, csv_path: str, columns: list[str] | None = None
) -> None:
    """Append a single row to the ablation CSV; write header if file is new."""
    if columns is None:
        columns = ABLATION_CSV_COLUMNS
    Path(csv_path).parent.mkdir(parents=True, exist_ok=True)
    write_header = not os.path.exists(csv_path)
    with open(csv_path, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def extract_ablation_metrics(result: AnnulusSimulationResult) -> dict:
    """Read nested Chinese summary keys -> flat English-keyed dict (robust to missing keys).

    AnnulusSimulationResult.summary nests metrics under summary["最终结果"].
    This helper extracts them into a flat English-keyed dict so callers don't
    need to know the Chinese key layout.
    """
    final = result.summary.get("最终结果", {}) if isinstance(result.summary, dict) else {}
    return {
        "effective_efficiency": final.get("全井段最终有效顶替效率"),
        "cement_occupation": final.get("最终水泥浆占据率"),
        "channeling_index": final.get("最终窜槽指数"),
        "mixing_index": final.get("最终混浆指数"),
        "instability_index": final.get("最终失稳指数"),
        "buoyancy_number": result.summary.get("buoyancy_number") if isinstance(result.summary, dict) else None,
    }


def run_one_level(
    level: AblationLevel,
    *,
    nz: int = 500,
    dt: float = 4.0,
    total_t: float | None = None,
    output_dir: str | None = None,
    well_spec_override: WellSpec | None = None,
    run_id: str | None = None,
    eccentricity: float = 0.17,
) -> AnnulusSimulationResult:
    """Run a single ablation level, returning the full simulation result.

    Pipeline (matches ht1_004_tailpipe.py):
    1. Load well data
    2. Run 1D casing flow solver
    3. Build coupled annulus inlet provider
    4. Compute annulus stop time (real, from casing result)
    5. Run 2D D2DGA solver with the level's switches
    6. Optionally dump summary JSON to a unique filename

    All three switches are forwarded to the solver. enable_true_buoyancy
    (Task 5) has no logic yet and is a placeholder.

    Parameters
    ----------
    well_spec_override : WellSpec | None
        If provided, use this well spec instead of the loaded one.
        Enables theoretical eccentricity cases without modifying the loader.
    run_id : str | None
        Unique identifier for this run. Used as the JSON filename stem
        (e.g. "ht1_004_ablation_R0", "theoretical_e35_R3").
        If None, falls back to ``level.name``.
    eccentricity : float
        Casing eccentricity (0..1). Default 0.17 matches HT1-004 field data.
    """
    loaded_well, fluids, schedule, _ = load_ht1_004_tailpipe()
    well_spec = well_spec_override if well_spec_override is not None else loaded_well

    # 1D casing flow
    casing_solver = CasingFlowSolver(enable_gravity=True)
    casing_result = casing_solver.run(well_spec, fluids, schedule)

    # Coupled inlet provider
    provider = build_coupled_annulus_inlet_provider(
        casing_result, casing_solver, fluids
    )

    # Annulus stop time (real, from casing result)
    if total_t is None:
        total_t = annulus_stop_time_s(casing_result=casing_result, fluids=fluids)

    # 2D D2DGA solver — only pass the switches that exist on the constructor
    solver = AnnulusD2DGASolver(
        dt=dt,
        nz=nz,
        ny=40,
        total_t=total_t,
        enable_d2dga=True,
        enable_d2dga_auto_m=level.enable_d2dga_auto_m,
        enable_d2dga_i3_flux=level.enable_d2dga_i3_flux,
        enable_true_buoyancy=level.enable_true_buoyancy,
        open_outlet=True,
    )

    result = solver.run(well_spec, fluids, provider)

    if output_dir is not None:
        _dump_summary(
            result, level, output_dir,
            run_id=run_id, nz=nz, dt=dt, eccentricity=eccentricity,
        )

    return result


def run_full_ablation(
    *,
    levels: list[AblationLevel] | None = None,
    nz: int = 500,
    dt: float = 4.0,
    total_t: float | None = None,
    output_dir: str | None = None,
    run_id_prefix: str = "ht1_004_ablation",
    eccentricity: float = 0.17,
) -> Dict[str, AnnulusSimulationResult]:
    """Run all ablation levels, returning {level_name: result}.

    Each level gets a unique run_id = ``{run_id_prefix}_{level.name}`` so the
    dumped JSONs are never overwritten by other scripts sharing the same
    output_dir.
    """
    if levels is None:
        levels = ABLATION_LEVELS
    results: Dict[str, AnnulusSimulationResult] = {}
    for lv in levels:
        run_id = f"{run_id_prefix}_{lv.name}"
        results[lv.name] = run_one_level(
            lv, nz=nz, dt=dt, total_t=total_t, output_dir=output_dir,
            run_id=run_id, eccentricity=eccentricity,
        )
    return results


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _dump_summary(
    result: AnnulusSimulationResult,
    level: AblationLevel,
    output_dir: str,
    *,
    run_id: str | None = None,
    nz: int,
    dt: float,
    eccentricity: float,
) -> None:
    """Write a per-level summary JSON combining ablation metrics + run metadata.

    The JSON is written to ``{output_dir}/{run_id or level.name}.json`` so
    every caller gets a unique file — no more filename collisions between
    ablation / theoretical / convergence scripts.
    """
    file_stem = run_id if run_id else level.name
    p = Path(output_dir) / f"{file_stem}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "run_id": run_id or f"ht1_004_ablation_{level.name}",
        "ablation_level": level.name,
        "eccentricity": eccentricity,
        "nz": nz,
        "dt": dt,
        "enable_d2dga_auto_m": level.enable_d2dga_auto_m,
        "enable_d2dga_i3_flux": level.enable_d2dga_i3_flux,
        "enable_true_buoyancy": level.enable_true_buoyancy,
        **extract_ablation_metrics(result),
    }
    p.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8",
    )