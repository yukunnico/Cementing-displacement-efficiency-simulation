"""Run paper-only main simulation results and write paper schemas."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from cemdisp.data.loaders.paper import (  # noqa: E402
    PaperWellMetadata,
    get_paper_metadata,
    iter_paper_wells,
    load_paper_well,
    metadata_to_dict,
)
from cemdisp.data.loaders.paper.common import get_window, write_csv, write_markdown_table  # noqa: E402
from cemdisp.models2d import AnnulusD2DGASolver  # noqa: E402
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider  # noqa: E402
from cemdisp.reporting.contour_plots import plot_final_fields_contour  # noqa: E402
from cemdisp.reporting.plots import plot_depth_profiles  # noqa: E402
from cemdisp.transport1d import CasingFlowSolver  # noqa: E402

PAPER_VERSION = "pead_v1"
RESULTS_ROOT = _PROJECT_ROOT / "results_paper" / PAPER_VERSION
ADOPTED_VALUES_PATH = _PROJECT_ROOT / "参考文档" / "现场资料提取" / "01_总表" / "loader优化采用值对照表.csv"

MAIN_RESULT_COLUMNS = [
    "well_id",
    "well_name_cn",
    "sample_class",
    "top_md_m",
    "bottom_md_m",
    "annulus_volume_m3",
    "final_effective_efficiency",
    "final_bulk_cement_fill",
    "cbl_eval_interval_efficiency",
    "target_interval_efficiency",
    "channeling_index",
    "mixing_index",
    "instability_index",
    "summary_json_path",
    "time_series_csv_path",
    "depth_profile_csv_path",
]


@dataclass(frozen=True)
class PaperModelConfig:
    """Numerical configuration for paper simulations."""

    dt: float = 4.0
    nz: int = 500
    ny: int = 40
    enable_d2dga: bool = True
    yield_regularization_M: float = 100.0
    open_outlet: bool = True
    enable_axial_dispersion: bool = True
    dispersion_alpha: float = 0.2
    mode: str = "formal"

    def to_dict(self) -> dict[str, object]:
        data = asdict(self)
        data["paper_version"] = PAPER_VERSION
        return data


def default_paper_model_config() -> PaperModelConfig:
    """Return the formal paper configuration from the design document."""

    return PaperModelConfig()


def smoke_paper_model_config() -> PaperModelConfig:
    """Return the low-cost smoke configuration used before formal runs."""

    return PaperModelConfig(dt=8.0, nz=120, ny=24, mode="smoke")


def ensure_results_tree(root: Path = RESULTS_ROOT) -> Path:
    """Create the paper results directory tree and copy the adopted-values manifest."""

    for relative in (
        "00_run_config",
        "01_main_results",
        "02_cbl_validation/figures",
        "03_d2dga_ablation/figures",
        "04_grid_dt_validation/figures",
        "05_figures_for_paper",
        "99_archive",
    ):
        (root / relative).mkdir(parents=True, exist_ok=True)
    manifest_dest = root / "00_run_config" / "adopted_loader_manifest.csv"
    if ADOPTED_VALUES_PATH.exists():
        shutil.copyfile(ADOPTED_VALUES_PATH, manifest_dest)
    elif not manifest_dest.exists():
        manifest_dest.write_text(
            "well_id,field_group,field_name,field_value,loader_value,adopted_value,unit,adopted_source_type,confidence,action,notes\n",
            encoding="utf-8-sig",
        )
    return root


def _json_default(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return str(value)


def _safe_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(result):
        return None
    return result


def _depth_interval_mean(depth_profiles: Any, top_md_m: float, bottom_md_m: float, column: str = "平均有效顶替效率") -> float | None:
    if column not in depth_profiles.columns or "井深_m" not in depth_profiles.columns:
        return None
    depths = depth_profiles["井深_m"]
    mask = (depths >= top_md_m) & (depths <= bottom_md_m)
    if not bool(mask.any()):
        return None
    return _safe_float(depth_profiles.loc[mask, column].mean())


def _result_metrics(result: Any, well_spec: Any) -> dict[str, float | None]:
    final = result.metrics.iloc[-1]
    cbl_window = get_window(well_spec, "cbl")
    target_window = get_window(well_spec, "target")
    cbl_eff = None
    target_eff = None
    if cbl_window is not None:
        cbl_eff = _depth_interval_mean(result.depth_profiles, cbl_window.top_md_m, cbl_window.bottom_md_m)
    if target_window is not None:
        target_eff = _depth_interval_mean(result.depth_profiles, target_window.top_md_m, target_window.bottom_md_m)
    return {
        "final_effective_efficiency": _safe_float(final.get("effective_efficiency")),
        "final_bulk_cement_fill": _safe_float(final.get("bulk_cement_fill")),
        "cbl_eval_interval_efficiency": cbl_eff if cbl_eff is not None else _safe_float(final.get("effective_efficiency")),
        "target_interval_efficiency": target_eff,
        "channeling_index": _safe_float(final.get("channeling_index")),
        "mixing_index": _safe_float(final.get("mixing_index")),
        "instability_index": _safe_float(final.get("instability_index")),
    }


def run_single_paper_well(
    well_id: str,
    config: PaperModelConfig,
    *,
    root: Path = RESULTS_ROOT,
) -> dict[str, object]:
    """Run one paper well and write its outputs."""

    ensure_results_tree(root)
    metadata = get_paper_metadata(well_id)
    well_spec, fluids, schedule, validation_data = load_paper_well(well_id)

    t0 = time.perf_counter()
    casing_solver = CasingFlowSolver(
        enable_gravity=True,
        enable_axial_dispersion=config.enable_axial_dispersion,
        dispersion_alpha=config.dispersion_alpha,
    )
    casing_result = casing_solver.run(well_spec, fluids, schedule)
    inlet_provider = build_coupled_annulus_inlet_provider(
        casing_result,
        casing_solver,
        fluids,
        split_cement_phases=True,
    )
    solver = AnnulusD2DGASolver(
        dt=config.dt,
        nz=config.nz,
        ny=config.ny,
        total_t=float(casing_result.cement_end_time_s),
        enable_d2dga=config.enable_d2dga,
        yield_regularization_M=config.yield_regularization_M,
        open_outlet=config.open_outlet,
    )
    result = solver.run(well_spec, fluids, inlet_provider)
    runtime_s = time.perf_counter() - t0

    return write_single_well_outputs(
        well_id=well_id,
        metadata=metadata,
        well_spec=well_spec,
        validation_data=validation_data,
        result=result,
        config=config,
        runtime_s=runtime_s,
        root=root,
    )


def write_single_well_outputs(
    *,
    well_id: str,
    metadata: PaperWellMetadata,
    well_spec: Any,
    validation_data: Any,
    result: Any,
    config: PaperModelConfig,
    runtime_s: float,
    root: Path = RESULTS_ROOT,
) -> dict[str, object]:
    """Write all required per-well paper result files."""

    well_dir = root / "01_main_results" / well_id
    well_dir.mkdir(parents=True, exist_ok=True)
    warnings: list[str] = []

    run_config = config.to_dict()
    run_config.update({"well_id": well_id, "loader_version": metadata.loader_version})
    (well_dir / "run_config.json").write_text(json.dumps(run_config, ensure_ascii=False, indent=2), encoding="utf-8")
    (well_dir / "adopted_input_summary.json").write_text(
        json.dumps(metadata_to_dict(metadata), ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    result.metrics.to_csv(well_dir / "time_series.csv", index=False, encoding="utf-8-sig")
    result.depth_profiles.to_csv(well_dir / "depth_profile.csv", index=False, encoding="utf-8-sig")
    np.savez_compressed(
        well_dir / "field_2d.npz",
        cement_final=result.cement_field,
        spacer_final=result.spacer_field,
        wall_final=result.wall_field,
        lead_final=result.lead_field,
        tail_final=result.tail_field,
        md=result.geom.get("md"),
        y=result.geom.get("y"),
    )

    metrics = _result_metrics(result, well_spec)
    summary: dict[str, object] = {
        "paper_version": PAPER_VERSION,
        "well_id": well_id,
        "well_name_cn": metadata.well_name_cn,
        "sample_class": metadata.sample_class,
        "include_in_cbl_metrics": metadata.include_in_cbl_metrics,
        "measured_cbl_pass_rate": metadata.cbl_pass_rate,
        "top_md_m": well_spec.top_md_m,
        "bottom_md_m": well_spec.bottom_md_m,
        "annulus_volume_m3": result.summary.get("物理环空体积_m3"),
        "runtime_s": runtime_s,
        "validation_notes": list(getattr(validation_data, "notes", ())),
        "warnings": warnings,
        **metrics,
    }
    (well_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=_json_default), encoding="utf-8")
    summary_lines = [
        f"# {metadata.well_name_cn} paper result summary",
        "",
        f"- well_id: {well_id}",
        f"- sample_class: {metadata.sample_class}",
        f"- final_effective_efficiency: {summary['final_effective_efficiency']}",
        f"- cbl_eval_interval_efficiency: {summary['cbl_eval_interval_efficiency']}",
        f"- runtime_s: {runtime_s:.2f}",
    ]
    (well_dir / "summary.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    provenance_lines = [
        f"# {metadata.well_name_cn} paper provenance",
        "",
        f"- Paper loader: `{well_id}` / `{metadata.loader_version}`",
        f"- Geometry source: {metadata.geometry_source}",
        f"- Pumping schedule source: {metadata.pumping_schedule_source}",
        f"- Fluid source: {metadata.fluid_source}",
        f"- Standoff source: {metadata.standoff_source}",
        f"- CBL source: {metadata.cbl_source}",
        f"- Include in CBL metrics: {metadata.include_in_cbl_metrics}",
        "",
        "## Notes",
        *[f"- {note}" for note in metadata.notes],
    ]
    (well_dir / "provenance.md").write_text("\n".join(provenance_lines) + "\n", encoding="utf-8")

    try:
        fig = plot_final_fields_contour(result, output_dir=well_dir)
        fig.savefig(well_dir / "final_field.png", dpi=300, bbox_inches="tight")
    except Exception as exc:  # plotting is optional for batch robustness
        warnings.append(f"final_field plot skipped: {exc}")
    try:
        fig = plot_depth_profiles(result, well_spec=well_spec, output_dir=well_dir)
        fig.savefig(well_dir / "depth_efficiency.png", dpi=300, bbox_inches="tight")
    except Exception as exc:
        warnings.append(f"depth profile plot skipped: {exc}")
    if warnings:
        (well_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=_json_default), encoding="utf-8")

    return {"well_id": well_id, "summary_path": well_dir / "summary.json", "runtime_s": runtime_s, **metrics}


def _relative(path: Path, root: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.as_posix()


def summary_row_from_files(summary_path: Path, metadata: PaperWellMetadata, output_root: Path = RESULTS_ROOT) -> dict[str, object]:
    """Build one stable main-result row from a per-well summary file."""

    data = json.loads(summary_path.read_text(encoding="utf-8"))
    well_dir = summary_path.parent
    row = {
        "well_id": data.get("well_id", metadata.well_id),
        "well_name_cn": metadata.well_name_cn,
        "sample_class": metadata.sample_class,
        "top_md_m": data.get("top_md_m"),
        "bottom_md_m": data.get("bottom_md_m"),
        "annulus_volume_m3": data.get("annulus_volume_m3"),
        "final_effective_efficiency": data.get("final_effective_efficiency"),
        "final_bulk_cement_fill": data.get("final_bulk_cement_fill"),
        "cbl_eval_interval_efficiency": data.get("cbl_eval_interval_efficiency"),
        "target_interval_efficiency": data.get("target_interval_efficiency"),
        "channeling_index": data.get("channeling_index"),
        "mixing_index": data.get("mixing_index"),
        "instability_index": data.get("instability_index"),
        "summary_json_path": _relative(summary_path, output_root),
        "time_series_csv_path": _relative(well_dir / "time_series.csv", output_root),
        "depth_profile_csv_path": _relative(well_dir / "depth_profile.csv", output_root),
    }
    return {column: row.get(column, "") for column in MAIN_RESULT_COLUMNS}


def write_all_wells_main_results(
    well_ids: Iterable[str],
    *,
    root: Path = RESULTS_ROOT,
) -> list[dict[str, object]]:
    """Write all-well main CSV and Markdown tables from per-well summaries."""

    rows: list[dict[str, object]] = []
    for well_id in well_ids:
        summary_path = root / "01_main_results" / well_id / "summary.json"
        if not summary_path.exists():
            continue
        rows.append(summary_row_from_files(summary_path, get_paper_metadata(well_id), root))
    write_csv(root / "01_main_results" / "all_wells_main_results.csv", rows, MAIN_RESULT_COLUMNS)
    write_markdown_table(root / "01_main_results" / "all_wells_main_results.md", rows, MAIN_RESULT_COLUMNS)
    return rows


def _well_ids_from_args(wells: Sequence[str] | None, *, include_pending: bool = True) -> list[str]:
    if wells:
        return list(wells)
    return [record.well_id for record in iter_paper_wells(include_pending=include_pending)]


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Run paper-only main results")
    parser.add_argument("--mode", choices=["smoke", "formal"], default="smoke")
    parser.add_argument("--wells", nargs="+", default=None)
    parser.add_argument("--output-root", type=Path, default=RESULTS_ROOT)
    args = parser.parse_args(argv)

    config = smoke_paper_model_config() if args.mode == "smoke" else default_paper_model_config()
    ensure_results_tree(args.output_root)
    config_path = args.output_root / "00_run_config" / "paper_model_config.json"
    config_path.write_text(json.dumps(config.to_dict(), ensure_ascii=False, indent=2), encoding="utf-8")
    well_ids = _well_ids_from_args(args.wells)
    run_log = args.output_root / "00_run_config" / "run_log.md"
    log_lines = [f"# Paper run log", "", f"- mode: {args.mode}", f"- wells: {', '.join(well_ids)}"]
    for well_id in well_ids:
        result = run_single_paper_well(well_id, config, root=args.output_root)
        log_lines.append(f"- {well_id}: runtime_s={result['runtime_s']:.2f}")
    write_all_wells_main_results(well_ids, root=args.output_root)
    run_log.write_text("\n".join(log_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
