# -*- coding: utf-8 -*-
"""全井「深度 × 时间 × 四通道份额」逐分钟导出（2026-09-26 交付物 6，第二次修订）。

与既有 export_fullwell_depth_time_shares 的三处不同（均为 2026-09-26 用户裁定）：
1. **相名统一为四通道**（既有脚本 1D 侧吐中文流体名，与 2D 侧四通道无法并表）；
2. **两个 source 并存**：套管内覆盖井口→井底全长，环空仅覆盖尾管段；
3. **按流体拆分、时间逐分钟**：每个通道一个文件，行=深度、列=逐分钟（宽格式）。

口径：
- 深度轴 = **井身结构表的 `depth_well_logging_m_` 列**（严格按表，不另造网格）。
  呼1-004 用用户提供的 333 点版本（放在输出目录内），其余井回退到
  `参考文档/<井>/<井名>井身结构.csv`。
- 时间轴 = **逐分钟**：0, 1, …, floor(t_max/60) min。管内侧可在任意时刻精确重建；
  环空侧 NPZ 只有快照帧，**逐分钟值由快照时间线性插值得到**（须向用户声明）。
- 管内侧：``build_casing_depth_profile`` 直接建在表的深度轴与分钟时间轴上（无插值）。
- 环空侧：NPZ 快照 → 方位角算术平均 → 先沿时间线性插值到分钟轴，再沿深度线性
  插值到表深度轴（仅尾管段内的深度）。
- 输出目录须由调用方指定，且拒绝覆盖已存在文件，绝不写入权威结果目录。

用法：
    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/entrypoints/export_unified_depth_time.py \
        --well ht1_004 --out results/深度时间统一口径_2026-09-26
"""

from __future__ import annotations

import argparse
import csv
import io
from pathlib import Path

import numpy as np

from cemdisp.transport1d.casing_flow import CasingFlowSolver
from cemdisp.transport1d.casing_depth_profile import (
    CHANNELS,
    build_casing_depth_profile,
    profile_channels,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
_FIELD_ROOT = PROJECT_ROOT.parent / "参考文档"
_CHANNEL_LABELS = {"lead": "领浆", "tail": "尾浆", "spacer": "隔离液", "mud": "钻井液"}

WELL_CONFIGS = {
    "ht1_003": {
        "well_label": "呼1-003",
        "loader": "cemdisp.data.loaders.ht1_003_loader:load_ht1_003_tailpipe",
        "loader_kwargs": {
            "reference_root": _FIELD_ROOT / "呼1-003" / "新",
            "caliper_csv_path": _FIELD_ROOT / "现场资料提取" / "ht1_003_呼1-003" / "caliper_profile.csv",
            "inclination_csv_path": _FIELD_ROOT / "现场资料提取" / "ht1_003_呼1-003" / "inclination_profile.csv",
        },
    },
    "ht1_004": {
        "well_label": "呼1-004",
        "loader": "cemdisp.data.loaders.ht1_004_loader:load_ht1_004_tailpipe",
        "loader_kwargs": {
            "reference_root": _FIELD_ROOT / "呼1-004",
            "caliper_csv_path": _FIELD_ROOT / "现场资料提取" / "ht1_004_呼1-004" / "caliper_profile.csv",
            "inclination_csv_path": _FIELD_ROOT / "现场资料提取" / "ht1_004_呼1-004" / "inclination_profile.csv",
        },
    },
}


def read_well_structure_depths(path: Path) -> np.ndarray:
    """井身结构表 → ``depth_well_logging_m_`` 列（升序）。"""

    rows = list(csv.reader(io.StringIO(Path(path).read_bytes().decode("utf-8-sig"))))
    if not rows:
        raise ValueError(f"井身结构表为空：{path}")
    col = next((i for i, name in enumerate(rows[0]) if "depth_well_logging" in name), None)
    if col is None:
        raise ValueError(f"井身结构表无 depth_well_logging 列：{path} 表头={rows[0]}")
    depths = np.array(
        [float(row[col]) for row in rows[1:] if len(row) > col and row[col].strip()],
        dtype=float,
    )
    if depths.size == 0:
        raise ValueError(f"井身结构表无深度数据：{path}")
    return np.sort(depths)


def resolve_well_structure(well_key: str, out_dir: Path, override: Path | None) -> Path:
    """井身结构表取用顺序：显式参数 > 输出目录内的用户版 > 参考文档标准版。"""

    if override is not None:
        return Path(override)
    label = WELL_CONFIGS[well_key]["well_label"]
    user_version = Path(out_dir) / f"{label}井身结构.csv"
    if user_version.exists():
        return user_version
    return _FIELD_ROOT / label / f"{label}井身结构.csv"


def _load_well(well_key: str):
    import importlib

    cfg = WELL_CONFIGS[well_key]
    module_name, func_name = cfg["loader"].split(":")
    loader = getattr(importlib.import_module(module_name), func_name)
    well_spec, fluids, schedule, _meta = loader(**cfg["loader_kwargs"])
    return cfg, well_spec, fluids, schedule


def _annulus_channels_on_grid(
    npz_path: Path, minute_s: np.ndarray, target_depths: np.ndarray
) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray]:
    """环空 NPZ → 分钟时间轴 × 表深度 的四通道份额。

    返回 ``(域内深度, {通道: (n_min, n_depth)}, 快照时刻)``。时间与深度两处均为
    **线性插值**（NPZ 只有快照帧与自身 ``md`` 网格），端点夹取、不外推。
    """

    data = np.load(npz_path, allow_pickle=True)
    snap_t = np.asarray(data["snapshot_times_s"], dtype=float)
    md = np.asarray(data["md"], dtype=float)
    order = np.argsort(md)
    md = md[order]

    def _mean(attr: str) -> np.ndarray | None:
        if attr not in data.files:
            return None
        return np.mean(np.asarray(data[attr], dtype=float), axis=1)[:, order]

    lead, tail = _mean("lead_snapshots"), _mean("tail_snapshots")
    if lead is None or tail is None:
        cement = _mean("cement_snapshots")
        if cement is None:
            raise KeyError(f"NPZ 既无 lead/tail 也无 cement 快照：{npz_path}")
        lead, tail = np.zeros_like(cement), cement
    spacer = _mean("spacer_snapshots")
    if spacer is None:
        spacer = np.zeros_like(lead)
    raw = {
        "lead": lead,
        "tail": tail,
        "spacer": spacer,
        "mud": np.clip(1.0 - lead - tail - spacer, 0.0, 1.0),
    }

    tol = 1.0e-6
    mask = (target_depths >= md[0] - tol) & (target_depths <= md[-1] + tol)
    in_depths = target_depths[mask]

    out: dict[str, np.ndarray] = {}
    for name, arr in raw.items():
        # (1) 时间插值：逐深度列插到分钟轴
        on_time = np.empty((minute_s.size, md.size), dtype=float)
        for k in range(md.size):
            on_time[:, k] = np.interp(
                minute_s, snap_t, arr[:, k], left=arr[0, k], right=arr[-1, k]
            )
        # (2) 深度插值：逐分钟行插到表深度轴（仅域内深度）
        out[name] = np.stack(
            [np.interp(in_depths, md, on_time[i]) for i in range(minute_s.size)]
        )
    return in_depths, out, snap_t


def export_unified_tables(
    well_key: str, out_dir: Path, *, well_structure: Path | None = None
) -> dict:
    """单井导出：每个通道一个「深度 × 逐分钟」宽表。返回摘要。"""

    cfg, well_spec, fluids, schedule = _load_well(well_key)
    label = cfg["well_label"]
    out_dir = Path(out_dir)
    npz_path = (
        PROJECT_ROOT / "results" / f"{label}_1D2D耦合模型"
        / f"{label}_1D2D耦合模型_2D场数据.npz"
    )

    ws_path = resolve_well_structure(well_key, out_dir, well_structure)
    depths = read_well_structure_depths(ws_path)

    snap_t = np.asarray(
        np.load(npz_path, allow_pickle=True)["snapshot_times_s"], dtype=float
    )
    n_min = int(np.floor(float(snap_t[-1]) / 60.0))
    minute_s = np.arange(n_min + 1, dtype=float) * 60.0

    solver = CasingFlowSolver()
    result = solver.run(well_spec, fluids, schedule)
    profile = build_casing_depth_profile(
        solver, well_spec, fluids, schedule, result,
        depths_m=depths, times_s=minute_s, mixing_band=True,
    )
    casing = profile_channels(profile)

    in_depths, annulus, _snap = _annulus_channels_on_grid(npz_path, minute_s, depths)

    out_dir.mkdir(parents=True, exist_ok=True)
    header = ["深度_m", "source"] + [f"{m}min" for m in range(n_min + 1)]
    written: list[str] = []
    for channel in CHANNELS:
        csv_path = out_dir / f"{label}_{_CHANNEL_LABELS[channel]}_逐分钟_深度时间表.csv"
        if csv_path.exists():
            raise FileExistsError(f"目标已存在，拒绝覆盖：{csv_path}")
        with csv_path.open("w", encoding="utf-8-sig", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(header)
            for j, z in enumerate(depths):
                writer.writerow([round(float(z), 3), "1D_casing"]
                                + [round(float(casing[channel][t, j]), 6)
                                   for t in range(minute_s.size)])
            for j, z in enumerate(in_depths):
                writer.writerow([round(float(z), 3), "2D_annulus"]
                                + [round(float(annulus[channel][t, j]), 6)
                                   for t in range(minute_s.size)])
        written.append(str(csv_path))

    return {
        "well": well_key,
        "label": label,
        "well_structure": str(ws_path),
        "casing_depth_range_m": [float(depths[0]), float(depths[-1])],
        "n_casing_rows": int(depths.size),
        "annulus_depth_range_m": [float(in_depths[0]), float(in_depths[-1])],
        "n_annulus_rows": int(in_depths.size),
        "n_minutes": n_min + 1,
        "t_max_s": float(snap_t[-1]),
        "n_snapshots": int(snap_t.size),
        "files": written,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="全井深度×时间四通道逐分钟导出")
    parser.add_argument("--well", required=True, choices=sorted(WELL_CONFIGS))
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--well-structure", type=Path, default=None)
    args = parser.parse_args()
    summary = export_unified_tables(args.well, args.out, well_structure=args.well_structure)
    for key, value in summary.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
