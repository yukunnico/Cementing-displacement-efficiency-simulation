# -*- coding: utf-8 -*-
"""环空各深度-各流体浓度随时间变化表导出。

从 2D 场数据 NPZ 或 AnnulusSimulationResult 中提取方位角平均浓度，
导出为 CSV（长格式：每行 = 一个时刻 × 一个深度）。

用法:
    # 方式1: 从 NPZ 文件（不必重新跑求解器）
    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/export_depth_time_concentration.py \
        results/呼1-004_1D2D耦合模型/呼1-004_1D2D耦合模型_2D场数据.npz

    # 方式2: 在 runner 中导入，求解完成后自动导出
    from scripts.export_depth_time_concentration import export_depth_time_csv
    export_depth_time_csv(result, output_dir, well_name="呼1-004")
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from cemdisp.models2d.annulus_d2dga import AnnulusSimulationResult


# ── 方位角平均 ─────────────────────────────────────────────────
def _azimuthal_mean(field: np.ndarray, geom: dict) -> np.ndarray:
    """对方位角方向做 b(半间隙)加权平均 → (n_snapshots, nz)。"""
    n_snapshots, ny, nz = field.shape
    b = geom.get("b", np.ones((ny, nz)))
    if b.ndim == 1:
        b = b[:, np.newaxis]

    result = np.zeros((n_snapshots, nz))
    for t_idx in range(n_snapshots):
        numerator = np.sum(field[t_idx] * b, axis=0)
        denominator = np.sum(b, axis=0)
        denominator = np.where(denominator < 1e-12, 1.0, denominator)
        result[t_idx] = numerator / denominator
    return result


# ── 主导出函数 ──────────────────────────────────────────────────
def export_depth_time_csv(
    result: "AnnulusSimulationResult",
    output_dir: Path,
    *,
    well_name: str = "",
    mode_title: str = "1D2D耦合模型",
) -> Path:
    """从求解结果导出「环空各深度 × 各时刻 × 各流体浓度」CSV 表。

    导出列:
        time_s, time_min, depth_m,
        [lead, tail,] cement_total, spacer, flusher, mud
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    md = result.geom["md"]
    nz = len(md)
    times_s = np.array(result.snapshot_times_s)
    n_snapshots = len(times_s)

    # ── 检测是否拆分了领浆/尾浆 ──
    has_split = (
        hasattr(result, "lead_snapshots")
        and result.lead_snapshots is not None
        and len(result.lead_snapshots) > 0
    )

    profiles: dict[str, np.ndarray] = {}
    if has_split:
        profiles["lead"] = _azimuthal_mean(np.array(result.lead_snapshots), result.geom)
        profiles["tail"] = _azimuthal_mean(np.array(result.tail_snapshots), result.geom)
        profiles["cement_total"] = profiles["lead"] + profiles["tail"]
    else:
        profiles["cement_total"] = _azimuthal_mean(
            np.array(result.cement_snapshots), result.geom
        )

    # spacer
    profiles["spacer"] = (
        _azimuthal_mean(np.array(result.spacer_snapshots), result.geom)
        if result.spacer_snapshots and len(result.spacer_snapshots) > 0
        else np.zeros((n_snapshots, nz))
    )

    # flusher (T1-6)
    flusher_snaps = getattr(result, "flusher_snapshots", None)
    profiles["flusher"] = (
        _azimuthal_mean(np.array(flusher_snaps), result.geom)
        if flusher_snaps is not None and len(flusher_snaps) > 0
        else np.zeros((n_snapshots, nz))
    )

    # mud = 1 - 其余相（五相闭合）
    mud_source = profiles["cement_total"] + profiles["spacer"] + profiles["flusher"]
    profiles["mud"] = np.clip(1.0 - mud_source, 0.0, 1.0)

    # ── 构建列序 ──
    fieldnames = ["time_s", "time_min", "depth_m"]
    if has_split:
        fieldnames.extend(["lead", "tail"])
    fieldnames.extend(["cement_total", "spacer", "flusher", "mud"])

    # ── 写 CSV ──
    prefix = f"{well_name}_" if well_name else ""
    csv_path = output_dir / f"{prefix}{mode_title}_深度时间浓度表.csv"

    with csv_path.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for t_idx in range(n_snapshots):
            ts = round(float(times_s[t_idx]), 3)
            tm = round(ts / 60.0, 4)
            for z_idx in range(nz):
                row: dict[str, float] = {
                    "time_s": ts,
                    "time_min": tm,
                    "depth_m": round(float(md[z_idx]), 3),
                }
                if has_split:
                    row["lead"] = round(float(profiles["lead"][t_idx, z_idx]), 6)
                    row["tail"] = round(float(profiles["tail"][t_idx, z_idx]), 6)
                row["cement_total"] = round(float(profiles["cement_total"][t_idx, z_idx]), 6)
                row["spacer"] = round(float(profiles["spacer"][t_idx, z_idx]), 6)
                row["flusher"] = round(float(profiles["flusher"][t_idx, z_idx]), 6)
                row["mud"] = round(float(profiles["mud"][t_idx, z_idx]), 6)
                writer.writerow(row)

    n_rows = n_snapshots * nz
    print(f"[浓度表] {n_snapshots}时刻×{nz}深度={n_rows}行 → {csv_path}")
    return csv_path


# ── NPZ 直读模式（不需要重新跑求解器）──────────────────────────
def export_from_npz(npz_path: str | Path, output_dir: str | Path | None = None) -> Path:
    """从已有 NPZ 文件直接导出浓度表。

    参数:
        npz_path: 2D 场数据 NPZ 文件路径
        output_dir: 输出目录，默认同 NPZ 所在目录
    """
    npz_path = Path(npz_path)
    output_dir = Path(output_dir) if output_dir else npz_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    data = np.load(npz_path, allow_pickle=True)
    md = data["md"]
    times_s = data["snapshot_times_s"]
    nz, n_snapshots = len(md), len(times_s)

    # 算术平均（无几何加权，NPZ 不含 b 场）
    def _mean(field):
        return np.mean(field, axis=1)

    has_lead = "lead_snapshots" in data
    has_tail = "tail_snapshots" in data

    if has_lead and has_tail:
        lead = _mean(data["lead_snapshots"])
        tail = _mean(data["tail_snapshots"])
        cement = lead + tail
        has_split = True
    elif "cement_snapshots" in data:
        cement = _mean(data["cement_snapshots"])
        has_split = False
    else:
        raise KeyError(f"NPZ 缺少 cement_snapshots，可用字段: {list(data.keys())}")

    spacer = _mean(data["spacer_snapshots"]) if "spacer_snapshots" in data else np.zeros((n_snapshots, nz))
    flusher = _mean(data["flusher_snapshots"]) if "flusher_snapshots" in data else np.zeros((n_snapshots, nz))
    mud = np.clip(1.0 - cement - spacer - flusher, 0.0, 1.0)

    # CSV 列
    fieldnames = ["time_s", "time_min", "depth_m"]
    if has_split:
        fieldnames.extend(["lead", "tail"])
    fieldnames.extend(["cement_total", "spacer", "flusher", "mud"])

    csv_path = output_dir / f"{npz_path.stem}_深度时间浓度表.csv"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for t_idx in range(n_snapshots):
            ts = round(float(times_s[t_idx]), 3)
            tm = round(ts / 60.0, 4)
            for z_idx in range(nz):
                row: dict[str, float] = {
                    "time_s": ts,
                    "time_min": tm,
                    "depth_m": round(float(md[z_idx]), 3),
                }
                if has_split:
                    row["lead"] = round(float(lead[t_idx, z_idx]), 6)
                    row["tail"] = round(float(tail[t_idx, z_idx]), 6)
                row["cement_total"] = round(float(cement[t_idx, z_idx]), 6)
                row["spacer"] = round(float(spacer[t_idx, z_idx]), 6)
                row["flusher"] = round(float(flusher[t_idx, z_idx]), 6)
                row["mud"] = round(float(mud[t_idx, z_idx]), 6)
                writer.writerow(row)

    n_rows = n_snapshots * nz
    print(f"[浓度表] {n_snapshots}时刻×{nz}深度={n_rows}行 → {csv_path}")
    return csv_path


# ── CLI ─────────────────────────────────────────────────────────
if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("用法: python export_depth_time_concentration.py <npz_path> [output_dir]")
        print("示例: python export_depth_time_concentration.py results/呼1-004_1D2D耦合模型/呼1-004_1D2D耦合模型_2D场数据.npz")
        sys.exit(1)

    npz_file = Path(sys.argv[1])
    out_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else None
    if not npz_file.exists():
        print(f"错误: 文件不存在 - {npz_file}")
        sys.exit(1)

    export_from_npz(npz_file, out_dir)
