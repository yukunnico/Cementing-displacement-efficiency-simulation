# -*- coding: utf-8 -*-
"""全井「深度 × 时间 × 四通道份额」统一口径导出（2026-09-26 交付物 6）。

与既有 `export_depth_time_concentration.export_fullwell_depth_time_shares` 的区别
只有一个，但是关键的：**两侧相名统一为四通道**。既有脚本的 1D 侧直接吐出中文
流体名（钻井液/先导浆/隔离液1/…），与 2D 侧的四通道（lead/tail/spacer/mud）
无法并表；本脚本用 `unify_phase_channel` 把 1D 侧也归到四通道。

口径：
- 时间栅格 = 该井 NPZ 的 `snapshot_times_s`（物理时间，非步数）；
- 深度 < 2D 评价域上界 → `source=1D_casing`，取管内剖面重建（带 erf 混浆带，
  2026-09-26 Q11）的四通道聚合份额；
- 深度 ≥ 上界 → `source=2D_annulus`，取 NPZ 快照的方位角算术平均（与既有导出
  同口径；NPZ 不含 b 场故不做 b 加权），mud = clip(1 − 其余三通道, 0, 1)；
- 输出目录须由调用方指定，且**拒绝覆盖已存在文件**，绝不写入权威结果目录。

用法：
    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/entrypoints/export_unified_depth_time.py \
        --well ht1_004 --out results/深度时间统一口径_2026-09-26
"""

from __future__ import annotations

import argparse
import csv
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

# 只有提供 pipe_id_profile 的井可做全井导出（另 5 井缺该字段，属数据层另案）
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

_HEADER = ("time_s", "time_min", "depth_m", "fluid", "share", "source")


def _load_well(well_key: str):
    import importlib

    cfg = WELL_CONFIGS[well_key]
    module_name, func_name = cfg["loader"].split(":")
    loader = getattr(importlib.import_module(module_name), func_name)
    well_spec, fluids, schedule, _meta = loader(**cfg["loader_kwargs"])
    return cfg, well_spec, fluids, schedule


def _annulus_channels(npz_path: Path) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray]:
    """NPZ → ({四通道: (n_t, nz)}, 时刻 (n_t,), 深度升序 (nz,))，四通道逐格闭合。"""

    data = np.load(npz_path, allow_pickle=True)
    times = np.asarray(data["snapshot_times_s"], dtype=float)
    md = np.asarray(data["md"], dtype=float)

    def _mean(attr: str) -> np.ndarray:
        return np.mean(np.asarray(data[attr], dtype=float), axis=1) if attr in data.files \
            else None

    lead, tail = _mean("lead_snapshots"), _mean("tail_snapshots")
    if lead is None or tail is None:
        cement = _mean("cement_snapshots")
        if cement is None:
            raise KeyError(f"NPZ 既无 lead/tail 也无 cement 快照：{npz_path}")
        lead, tail = np.zeros_like(cement), cement
    spacer = _mean("spacer_snapshots")
    if spacer is None:
        spacer = np.zeros_like(lead)
    channels = {
        "lead": lead,
        "tail": tail,
        "spacer": spacer,
        "mud": np.clip(1.0 - lead - tail - spacer, 0.0, 1.0),
    }

    order = np.argsort(md)                      # 深度升序
    channels = {name: arr[:, order] for name, arr in channels.items()}
    return channels, times, md[order]


def export_unified_long_table(well_key: str, out_dir: Path, *, n_casing_depths: int = 200) -> dict:
    """单井导出统一口径长表；返回统计摘要。"""

    cfg, well_spec, fluids, schedule = _load_well(well_key)
    label = cfg["well_label"]
    results_dir = PROJECT_ROOT / "results" / f"{label}_1D2D耦合模型"
    npz_path = results_dir / f"{label}_1D2D耦合模型_2D场数据.npz"

    ann_channels, times, md_asc = _annulus_channels(npz_path)
    domain_top_m = float(md_asc[0])

    solver = CasingFlowSolver()
    result = solver.run(well_spec, fluids, schedule)

    # 管内段深度网格：地面 → 2D 评价域上界（不含端点，避免与 2D 侧重复）
    casing_depths = np.linspace(0.0, domain_top_m, n_casing_depths + 1)[:-1]
    profile = build_casing_depth_profile(
        solver, well_spec, fluids, schedule, result,
        depths_m=casing_depths, times_s=times, mixing_band=True,
    )
    casing_channels = profile_channels(profile)

    depths = np.concatenate([casing_depths, md_asc])
    n_casing = casing_depths.size
    rows = []
    for t_idx in range(times.size):
        ts = round(float(times[t_idx]), 3)
        tm = round(ts / 60.0, 4)
        for j, z in enumerate(depths):
            zd = round(float(z), 3)
            if j < n_casing:
                for ch in CHANNELS:
                    rows.append((ts, tm, zd, ch, round(float(casing_channels[ch][t_idx, j]), 6),
                                 "1D_casing"))
            else:
                k = j - n_casing
                for ch in CHANNELS:
                    rows.append((ts, tm, zd, ch, round(float(ann_channels[ch][t_idx, k]), 6),
                                 "2D_annulus"))

    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"{label}_全井深度时间四通道_长格式.csv"
    if csv_path.exists():
        raise FileExistsError(f"目标已存在，拒绝覆盖：{csv_path}")
    with csv_path.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(_HEADER)
        writer.writerows(rows)

    return {
        "well": well_key,
        "label": label,
        "csv": str(csv_path),
        "n_rows": len(rows),
        "n_times": int(times.size),
        "n_depths": int(depths.size),
        "n_casing_depths": int(n_casing),
        "domain_top_m": domain_top_m,
        "shoe_md_m": float(well_spec.shoe_md_m),
        "t_min_s": float(times[0]),
        "t_max_s": float(times[-1]),
        "channels": list(CHANNELS),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="全井深度×时间四通道统一口径导出")
    parser.add_argument("--well", required=True, choices=sorted(WELL_CONFIGS))
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    summary = export_unified_long_table(args.well, args.out)
    for key, value in summary.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
