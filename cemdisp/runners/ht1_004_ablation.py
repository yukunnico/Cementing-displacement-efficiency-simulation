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

from dataclasses import dataclass
from typing import Dict
import json
from pathlib import Path

from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver, AnnulusSimulationResult
from cemdisp.transport1d import CasingFlowSolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.runners.ht1_004_tailpipe import annulus_stop_time_s
from cemdisp.data.loaders.ht1_004_loader import load_ht1_004_tailpipe


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
    }


def run_one_level(
    level: AblationLevel,
    *,
    nz: int = 500,
    dt: float = 4.0,
    total_t: float | None = None,
    output_dir: str | None = None,
) -> AnnulusSimulationResult:
    """Run a single ablation level, returning the full simulation result.

    Pipeline (matches ht1_004_tailpipe.py):
    1. Load well data
    2. Run 1D casing flow solver
    3. Build coupled annulus inlet provider
    4. Compute annulus stop time (real, from casing result)
    5. Run 2D D2DGA solver with the level's switches
    6. Optionally dump summary JSON

    All three switches are forwarded to the solver. enable_true_buoyancy
    (Task 5) has no logic yet and is a placeholder.
    """
    well_spec, fluids, schedule, _ = load_ht1_004_tailpipe()

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
        _dump_summary(result, level, output_dir)

    return result


def run_full_ablation(
    *,
    levels: list[AblationLevel] | None = None,
    nz: int = 500,
    dt: float = 4.0,
    total_t: float | None = None,
    output_dir: str | None = None,
) -> Dict[str, AnnulusSimulationResult]:
    """Run all ablation levels, returning {level_name: result}."""
    if levels is None:
        levels = ABLATION_LEVELS
    results: Dict[str, AnnulusSimulationResult] = {}
    for lv in levels:
        results[lv.name] = run_one_level(
            lv, nz=nz, dt=dt, total_t=total_t, output_dir=output_dir
        )
    return results


def _dump_summary(
    result: AnnulusSimulationResult, level: AblationLevel, output_dir: str
) -> None:
    """Write a per-level summary JSON combining ablation metrics + level flags."""
    p = Path(output_dir) / f"ht1_004_ablation_{level.name}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    payload = extract_ablation_metrics(result)
    payload["ablation_level"] = level.name
    payload["enable_d2dga_auto_m"] = level.enable_d2dga_auto_m
    payload["enable_d2dga_i3_flux"] = level.enable_d2dga_i3_flux
    payload["enable_true_buoyancy"] = level.enable_true_buoyancy
    p.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8",
    )