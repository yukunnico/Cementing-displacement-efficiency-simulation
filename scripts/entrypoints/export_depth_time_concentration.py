# -*- coding: utf-8 -*-
"""环空各深度-各流体浓度随时间变化表导出。

从 2D 场数据 NPZ 或 AnnulusSimulationResult 中提取方位角平均浓度，
导出为 CSV（长格式：每行 = 一个时刻 × 一个深度）。

用法:
    # 方式1: 从 NPZ 文件（不必重新跑求解器）
    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/entrypoints/export_depth_time_concentration.py \
        results/呼1-004_1D2D耦合模型/呼1-004_1D2D耦合模型_2D场数据.npz

    # 方式2: 在 runner 中导入，求解完成后自动导出
    from scripts.entrypoints.export_depth_time_concentration import export_depth_time_csv
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


# ── 全井深度×时间×流体占比导出（环空2D插值 + 套管1D sharp重建）──────────────
#
# 口径（2026-09-09 用户批准）：
# - 环空段（2D 评价域 = NPZ md min→max，悬挂器→尾管鞋）：NPZ 快照算术周向平均
#   （np.mean axis=1，与 export_from_npz 同口径）+ 沿深度线性插值到井身结构表深度；
#   四通道 lead/tail/spacer/mud，mud = clip(1 − 其余, 0, 1) 闭合。
# - 套管内段（评价域外深度）：1D 解析重建 sharp 口径（无混浆带）。
#   体积账：V(t) = 分段线性累计泵入体积（停泵段不增；碰压后冻结在 S_尾浆+Vp_shoe）；
#   深度 z 的管容坐标 v(z) = Vp(z)（pipe_id_profile 梯形积分，鞋口处标定到 1D 求解器
#   管容）；流经体积坐标 u = V(t) − v(z)；u < 0 → 初始钻井液，S_{j-1} ≤ u < S_j →
#   泵序第 j 步流体（半开区间，与 CasingFlowSolver._fluid_by_injected_volume 同口径）。
#   碰压后 V 冻结 → 胶塞面（尾浆/压塞液界面）停在鞋口，替浆/压塞液保持在胶塞上方
#   不进环空（2026-09-03 胶塞工艺裁定）。时间语义 = 地面累计自开泵（与 2D
#   snapshot_times_s 同轴）。
# - 自校验不过即报错拒绝导出：①各流体前缘到达鞋口时刻（体积反解，与
#   CasingFlowSolver._front_arrival_time 同口径）vs 权威 鞋口出流时序.csv 逐项对比；
#   ②重建 t_bump vs NPZ 最后快照时刻（报告差值）。

_FULL_WELL_CHANNELS = ("lead", "tail", "spacer", "mud")
_DOMAIN_TOL_M = 1e-6     # 2D 评价域边界容差
_ARRIVAL_TOL_S = 1.0     # 到达时刻对比容差（鞋口 CSV 保留 3 位小数，秒级口径）


def map_fluid_to_channel(fluid_name: str) -> str:
    """泵序流体名 → 四通道（lead/tail/spacer/mud）。未知名报错，不静默归 mud。"""
    name = fluid_name.strip()
    if "领浆" in name or "中间浆" in name:
        return "lead"
    if "尾浆" in name or "尾管水泥浆" in name:
        return "tail"
    if ("隔离液" in name or "平衡液" in name or "先导浆" in name or "先导液" in name
            or "冲洗" in name or "WASH" in name.upper()):
        return "spacer"
    if ("钻井液" in name or "泥浆" in name or "压塞液" in name or "保护液" in name
            or "基液" in name or "井浆" in name):
        return "mud"
    raise ValueError(f"无法映射泵序流体名到四通道: {fluid_name!r}")


def build_pipe_volume_profile(
    profile_depths,
    profile_ids_mm,
    shoe_md_m: float,
    pipe_volume_target: float,
):
    """管容 Vp(z)：pipe_id_profile 梯形积分，鞋口处标定到 1D 求解器管容。

    标定系数 k = pipe_volume_target / 剖面梯形积分（管容双链对齐：
    pipe_id_profile 剖面积分 vs shoe_lag_volume_m3 现场口径，呼1-003/004 k≈1）。
    返回 (升序深度数组, Vp 数组, k)。
    """
    depths = np.asarray(profile_depths, dtype=float)
    ids_mm = np.asarray(profile_ids_mm, dtype=float)
    order = np.argsort(depths)
    depths, ids_mm = depths[order], ids_mm[order]
    # 裁剪到 [0, shoe]（与 CasingFlowSolver._pipe_cross_section_area 的 mask 口径一致）
    keep = (depths >= -_DOMAIN_TOL_M) & (depths <= shoe_md_m + _DOMAIN_TOL_M)
    depths, ids_mm = depths[keep], ids_mm[keep]
    if len(depths) < 2:
        raise ValueError("pipe_id_profile 在 [0, shoe] 内不足 2 个点，无法积分管容")
    if depths[0] > _DOMAIN_TOL_M:
        depths = np.concatenate(([0.0], depths))
        ids_mm = np.concatenate(([ids_mm[0]], ids_mm))
    areas = np.pi * (ids_mm / 1000.0) ** 2 / 4.0
    seg = (areas[:-1] + areas[1:]) * 0.5 * np.diff(depths)
    vp = np.concatenate(([0.0], np.cumsum(seg)))
    raw_shoe = float(vp[-1])
    if raw_shoe <= 0.0:
        raise ValueError("pipe_id_profile 积分管容为零，请检查内径剖面")
    k = float(pipe_volume_target) / raw_shoe
    return depths, vp * k, k


def invert_volume_to_time(scheduled_steps, target_volume_m3: float) -> float | None:
    """体积坐标 → 时刻（与 CasingFlowSolver._front_arrival_time 同款分段线性反解）。"""
    for scheduled in scheduled_steps:
        if target_volume_m3 <= scheduled.cumulative_volume_end_m3 + 1.0e-12:
            if scheduled.step.rate_m3_min <= 0.0:
                return scheduled.end_time_s
            volume_into_step = max(target_volume_m3 - scheduled.cumulative_volume_start_m3, 0.0)
            return scheduled.start_time_s + volume_into_step / scheduled.step.rate_m3_min * 60.0
    return None


def make_volume_of_time_fn(scheduled_steps, v_frozen: float):
    """时刻 → 累计泵入体积 V(t)：分段线性（停泵段不增），冻结在 v_frozen。

    复用 CasingFlowSolver._cumulative_volume_at 保证与 1D 求解器同一积分口径。
    """
    from cemdisp.transport1d.casing_flow import CasingFlowSolver

    def volume_of_time(t_s: float) -> float:
        return min(CasingFlowSolver._cumulative_volume_at(scheduled_steps, t_s), v_frozen)

    return volume_of_time


def fluid_name_at_volume_coordinate(u: float, scheduled_steps, initial_fluid_name: str) -> str:
    """sharp 界面查表：管内容积坐标 u（= V(t) − v(z)）处的流体名。

    u < 0 → 初始钻井液；否则半开区间 [S_{j-1}, S_j) 查表（与
    CasingFlowSolver._fluid_by_injected_volume 同口径，u ≥ S_末 → 末步流体）。
    """
    if u < 0.0:
        return initial_fluid_name
    from cemdisp.transport1d.casing_flow import CasingFlowSolver

    return CasingFlowSolver._fluid_by_injected_volume(scheduled_steps, u)


def interp_annulus_profiles(mean_fields: dict, md_descending, table_depths, tol_m: float = _DOMAIN_TOL_M):
    """2D 周向平均场沿深度线性插值到井身结构表深度（仅域内，不外推）。

    返回 (域内掩码, 域内深度数组, {通道: (n_t, n_in) 数组})。
    """
    md_asc = np.flipud(np.asarray(md_descending, dtype=float))
    md_min, md_max = float(md_asc[0]), float(md_asc[-1])
    table = np.asarray(table_depths, dtype=float)
    in_mask = (table >= md_min - tol_m) & (table <= md_max + tol_m)
    in_depths = table[in_mask]
    in_vals: dict[str, np.ndarray] = {}
    for ch, field in mean_fields.items():
        # 场与 md 同步翻转（md 降序=井底在前 → 升序=浅在前）：若只翻 md 不翻场，
        # np.interp 会对不齐的 x/y 静默产出镜像剖面（2026-09-09 实锤缺陷）。
        # 场形状 (n_t, nz)：深度在 axis 1，flipud 翻的是 axis 0（时间维）——必须用列反转。
        field_asc = np.asarray(field, dtype=float)[:, ::-1]
        in_vals[ch] = np.stack([np.interp(in_depths, md_asc, field_asc[k])
                                for k in range(field_asc.shape[0])])
    return in_mask, in_depths, in_vals


def check_no_overwrite(paths) -> None:
    """防覆盖守卫：任一目标文件已存在即拒绝（权威目录不覆盖现有文件）。"""
    existing = [Path(p) for p in paths if Path(p).exists()]
    if existing:
        raise RuntimeError("拒绝覆盖已存在文件: " + ", ".join(str(p) for p in existing))


def validate_shoe_arrivals(mine_rows, csv_rows, tol_s: float = _ARRIVAL_TOL_S) -> list[dict]:
    """重建到达时刻 vs 鞋口出流时序 CSV 逐项对比。

    mine_rows: [(流体名, 到达时刻|None)] 按泵序；csv_rows: [(地面累计时间_s, 鞋口出流流体)]。
    同名多步按时间序游标消耗配对；t≈0 行（初始管内流体）与泵注结束行（迟到体积封顶
    的出流快照）视为无害残余。返回 mismatch 列表（空 = 通过）。
    """
    csv_max = max((t for t, _ in csv_rows), default=0.0)
    pool: dict[str, list[float]] = {}
    for t, name in csv_rows:
        if t > tol_s:  # t≈0 行是初始流体状态，不是前缘到达
            pool.setdefault(name, []).append(t)
    cursors: dict[str, int] = {}
    mismatches: list[dict] = []
    for fluid_name, t_mine in mine_rows:
        times = pool.get(fluid_name, [])
        i = cursors.get(fluid_name, 0)
        if t_mine is None:
            for t in times[i:]:
                if abs(t - csv_max) > tol_s:
                    mismatches.append({"fluid": fluid_name, "mine_s": None,
                                       "csv_s": t, "delta_s": None})
            cursors[fluid_name] = len(times)
        else:
            if i >= len(times) or abs(times[i] - t_mine) > tol_s:
                mismatches.append({
                    "fluid": fluid_name,
                    "mine_s": t_mine,
                    "csv_s": times[i] if i < len(times) else None,
                    "delta_s": None if i >= len(times) else times[i] - t_mine,
                })
            else:
                cursors[fluid_name] = i + 1
    return mismatches


def write_csv_utf8sig(path, header, rows) -> None:
    """utf-8-sig（带 BOM）写 CSV，rows 为元组迭代器。"""
    with Path(path).open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)


def _read_table_depths(table_csv_path) -> list[float]:
    """读井身结构表深度列（第 2 列 depth_well_logging_m_，utf-8-sig，段底累计深度）。"""
    with Path(table_csv_path).open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        col = next((i for i, name in enumerate(header) if "depth_well_logging" in name), 1)
        depths = [float(row[col]) for row in reader if row and row[col].strip()]
    if not depths:
        raise ValueError(f"井身结构表无深度数据: {table_csv_path}")
    return depths


def _read_shoe_timeline_rows(shoe_csv_path) -> list[tuple[float, str]]:
    """读权威鞋口出流时序 CSV → [(地面累计时间_s, 鞋口出流流体)]。"""
    rows: list[tuple[float, str]] = []
    with Path(shoe_csv_path).open("r", encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            rows.append((float(row["地面累计时间_s"]), row["鞋口出流流体"]))
    if not rows:
        raise ValueError(f"鞋口出流时序 CSV 为空: {shoe_csv_path}")
    return rows


def _last_cement_step(scheduled_steps, fluids):
    """定位末段水泥浆步骤（LEAD/INTERMEDIATE/TAIL 角色，名称兜底）→ S_尾浆。"""
    from cemdisp.data.fluid_spec import FluidRole

    role_by_name = {f.name: f.role for f in fluids}
    cement_roles = {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}
    cement = [s for s in scheduled_steps if role_by_name.get(s.step.fluid_name) in cement_roles]
    if not cement:
        cement = [s for s in scheduled_steps
                  if ("尾浆" in s.step.fluid_name or "水泥" in s.step.fluid_name
                      or "领浆" in s.step.fluid_name)]
    if not cement:
        raise ValueError("泵序中未识别到水泥浆步骤（LEAD/INTERMEDIATE/TAIL），无法定位 S_尾浆")
    return cement[-1]


def export_fullwell_depth_time_shares(
    npz_path,
    table_csv_path,
    output_dir,
    shoe_timeline_csv_path,
    *,
    well_label: str,
    schedule,
    fluids,
    pipe_profile_depths,
    pipe_profile_ids_mm,
    shoe_md_m: float,
    pipe_volume_target: float,
    mode_title: str = "1D2D耦合模型",
) -> dict:
    """全井深度×时间×流体占比导出：环空段 2D 插值 + 套管段 1D sharp 重建。

    零重跑：仅读权威 NPZ + loader 泵序（与 runner 传给 1D 求解器同源）+ 井身结构表。
    自校验不过即报错拒绝导出。每井产出 6 个文件（长格式 + 4 宽格式 + 说明），
    输出到权威结果目录顶层，不覆盖任何现有文件。
    """
    from cemdisp.transport1d.casing_flow import CasingFlowSolver

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # ── 输入加载 ──
    data = np.load(npz_path, allow_pickle=True)
    md = np.asarray(data["md"], dtype=float)
    times_s = np.asarray(data["snapshot_times_s"], dtype=float)
    n_t, nz = len(times_s), len(md)
    mean_fields = {}
    for ch, key in (("lead", "lead_snapshots"), ("tail", "tail_snapshots"),
                    ("spacer", "spacer_snapshots")):
        if key not in data:
            raise KeyError(f"NPZ 缺少 {key}，无法构建四通道: {list(data.keys())}")
        mean_fields[ch] = np.mean(np.asarray(data[key], dtype=float), axis=1)  # 算术周向平均
    table_depths = _read_table_depths(table_csv_path)
    n_depths = len(table_depths)

    # ── 管容 / 泵序 / 碰压 ──
    vp_depths, vp_vals, k_scale = build_pipe_volume_profile(
        pipe_profile_depths, pipe_profile_ids_mm, shoe_md_m, pipe_volume_target)
    steps_full = CasingFlowSolver._build_scheduled_steps(schedule)
    steps, _ = CasingFlowSolver._displacement_sequence_cutoff(steps_full)
    initial_fluid = CasingFlowSolver._initial_fluid_name(fluids, schedule)
    cement_step = _last_cement_step(steps, fluids)
    s_tail = cement_step.cumulative_volume_end_m3
    v_total = steps[-1].cumulative_volume_end_m3
    v_bump = s_tail + pipe_volume_target          # 尾浆全部出鞋口（≡ 现场碰压断面）
    v_frozen = min(v_bump, v_total)               # 未碰压井冻结在序列终点
    t_bump = invert_volume_to_time(steps, s_tail + pipe_volume_target)
    volume_of_time = make_volume_of_time_fn(steps, v_frozen)

    # ── 自校验①：前缘到达鞋口时刻 vs 权威鞋口 CSV ──
    mine_rows = [(s.step.fluid_name, CasingFlowSolver._front_arrival_time(s, steps, pipe_volume_target))
                 for s in steps]
    shoe_rows = _read_shoe_timeline_rows(shoe_timeline_csv_path)
    arrival_mismatches = validate_shoe_arrivals(mine_rows, shoe_rows, tol_s=_ARRIVAL_TOL_S)
    if arrival_mismatches:
        lines = [f"  {m['fluid']}: 重建={m['mine_s']} CSV={m['csv_s']} Δ={m['delta_s']}"
                 for m in arrival_mismatches]
        raise RuntimeError("1D 重建到达时刻与鞋口出流时序.csv 不符（容差 "
                           f"{_ARRIVAL_TOL_S}s），拒绝导出:\n" + "\n".join(lines))

    # ── 自校验②：t_bump vs NPZ 最后快照（报告差值，不拦截）──
    last_snapshot_s = float(times_s[-1])
    t_bump_delta_s = (t_bump - last_snapshot_s) if t_bump is not None else None

    # ── 深度域判定与 2D 插值 ──
    in_mask, in_depths, in_vals = interp_annulus_profiles(mean_fields, md, table_depths)
    mud_in = np.clip(1.0 - in_vals["lead"] - in_vals["tail"] - in_vals["spacer"], 0.0, 1.0)
    n_in, n_out = int(in_mask.sum()), int((~in_mask).sum())

    # ── 套管段 sharp 状态：每个域外深度的 Vp 坐标 ──
    out_indices = [j for j in range(n_depths) if not in_mask[j]]
    out_vps = [float(np.interp(min(table_depths[j], shoe_md_m), vp_depths, vp_vals))
               for j in out_indices]

    def casing_fluid_at(vp_z: float, t_idx: int) -> str:
        u = volume_of_time(times_s[t_idx]) - vp_z
        return fluid_name_at_volume_coordinate(u, steps, initial_fluid)

    # ── 防覆盖 ──
    prefix = f"{well_label}_{mode_title}"
    long_csv_path = output_dir / f"{prefix}_全井深度时间占比_长格式.csv"
    doc_path = output_dir / f"{prefix}_全井深度时间占比_说明.md"
    wide_paths = {ch: output_dir / f"{prefix}_宽格式_{ch}.csv" for ch in _FULL_WELL_CHANNELS}
    check_no_overwrite([long_csv_path, doc_path, *wide_paths.values()])

    # ── 长格式：环空段四通道全输出；套管段仅输出 share>0 相（sharp 通常 1 个）──
    def _long_rows():
        for t_idx in range(n_t):
            ts = round(float(times_s[t_idx]), 3)
            tm = round(ts / 60.0, 4)
            it_idx = 0
            for j, z in enumerate(table_depths):
                zd = round(float(z), 3)
                if in_mask[j]:
                    for ch, arr in (("lead", in_vals["lead"]), ("tail", in_vals["tail"]),
                                    ("spacer", in_vals["spacer"]), ("mud", mud_in)):
                        yield (ts, tm, zd, ch, round(float(arr[t_idx, it_idx]), 6), "2D_annulus")
                    it_idx += 1
                else:
                    fluid = casing_fluid_at(out_vps[out_indices.index(j)], t_idx)
                    yield (ts, tm, zd, fluid, 1.0, "1D_casing")

    long_header = ("time_s", "time_min", "depth_m", "fluid", "share", "source")
    write_csv_utf8sig(long_csv_path, long_header, _long_rows())
    n_long_rows = n_t * (n_in * 4 + n_out)

    # ── 宽格式 ×4：行 = 全部表格深度（域外行用 1D 四通道映射值填充，无空值）──
    wide_values = {ch: np.zeros((n_depths, n_t)) for ch in _FULL_WELL_CHANNELS}
    for j_local, j in enumerate(range(n_depths)):
        if in_mask[j]:
            i_in = int(np.count_nonzero(in_mask[: j + 1])) - 1
            for t_idx in range(n_t):
                for ch in ("lead", "tail", "spacer"):
                    wide_values[ch][j, t_idx] = in_vals[ch][t_idx, i_in]
                wide_values["mud"][j, t_idx] = mud_in[t_idx, i_in]
        else:
            vp_z = out_vps[out_indices.index(j)]
            for t_idx in range(n_t):
                fluid = casing_fluid_at(vp_z, t_idx)
                wide_values[map_fluid_to_channel(fluid)][j, t_idx] = 1.0
    wide_header = ["深度_m"] + [f"{float(t):.3f}_s" for t in times_s]
    for ch, path in wide_paths.items():
        rows = ((round(float(z), 3), *(round(float(wide_values[ch][j, t]), 6)
                                       for t in range(n_t)))
                for j, z in enumerate(table_depths))
        write_csv_utf8sig(path, wide_header, rows)

    # ── 说明文档 ──
    csv_arrival_by_name: dict[str, float] = {}
    for t, name in shoe_rows:
        if t > _ARRIVAL_TOL_S and name not in csv_arrival_by_name:
            csv_arrival_by_name[name] = t
    arrival_lines = [
        f"| {name} | {'—' if t is None else f'{t:.3f}'} | "
        f"{'—' if csv_arrival_by_name.get(name) is None else f'{csv_arrival_by_name[name]:.3f}'} |"
        for name, t in mine_rows
    ]
    # 同名流体去重（保持首现顺序）
    seen: set[str] = set()
    unique_srcs = [s for s in ([x.step.fluid_name for x in steps] + [initial_fluid])
                   if not (s in seen or seen.add(s))]
    channel_map_lines = "\n".join(f"| {src} | {map_fluid_to_channel(src)} |" for src in unique_srcs)
    doc = f"""# {well_label} 全井深度×时间×流体占比 导出说明

## 数据来源（零重跑）
- 环空 2D 场：`{Path(npz_path).name}`（权威结果目录 NPZ，快照 {n_t} 个 × 模型深度 {nz} 层）
- 泵序/管容：`{well_label}` loader 泵序（与 runner 传给 1D 求解器同源），共 {len(steps)} 步，总泵入 {v_total:.3f} m³
- 深度轴：井身结构表 `{Path(table_csv_path).name}`（{n_depths} 点）
- 校验基准：`{Path(shoe_timeline_csv_path).name}`

## 口径声明
1. **环空段（source=2D_annulus）**：2D 评价域 = NPZ md min→max（悬挂器→尾管鞋），
   本井域 [{md.min():.3f}, {md.max():.3f}] m；周向**算术平均**（NPZ 无 b 场，
   `np.mean(axis=1)`，与现有 `export_from_npz` 同口径）+ 沿深度**线性插值**到表格深度；
   四通道 lead/tail/spacer/mud，mud = clip(1 − 其余, 0, 1) 闭合。
   域内表格深度点数：**{n_in}/{n_depths}**（域外 {n_out} 点由 1D 重建覆盖）。
   NPZ 中另有 wall 场（壁面相）未计入四通道闭合（与现有导出口径一致）。
2. **套管内段（source=1D_casing）**：1D 解析重建 **sharp 界面**口径（无混浆带）。
   体积账：V(t) = 泵序分段线性累计泵入体积（停泵段不增；**碰压后冻结**在
   S_尾浆 + Vp_shoe = {s_tail:.3f} + {pipe_volume_target:.3f} = {v_bump:.3f} m³）；
   流经体积坐标 u = V(t) − v(z)，v(z) = Vp(z) 为管容坐标（井口→z 梯形积分，
   鞋口处标定到 1D 求解器管容，标定系数 k = {k_scale:.6f}）；
   u < 0 → 初始钻井液（{initial_fluid}）；S_(j-1) ≤ u < S_j → 泵序第 j 步流体（半开区间，
   与 CasingFlowSolver._fluid_by_injected_volume 同口径）。
   碰压后 V 冻结 → 胶塞面（尾浆/压塞液界面）停在鞋口，替浆/压塞液保持在胶塞上方
   **不进环空**（2026-09-03 用户胶塞工艺裁定）；未碰压井则尾浆尾段滞留管内。
   **套管内 = 管内流体（非环空）**，与环空段在悬挂器处不连续属物理事实。
3. **时间语义**：地面累计自开泵（秒），与 2D 快照时刻轴同轴（{n_t} 个快照时刻）。
4. **四通道映射**（1D 套管段流体名 → 宽格式通道）：

| 泵序流体名 | 通道 |
|---|---|
{channel_map_lines}

   映射规则：隔离液*/平衡液/先导浆(液)/WASH→spacer；领浆/中间浆→lead；尾浆→tail；
   钻井液*/替钻井液/压塞液/保护液/基液/井浆→mud。

## 自校验结果
### ① 前缘到达鞋口时刻：1D 重建 vs 鞋口出流时序.csv（容差 {_ARRIVAL_TOL_S}s）

| 泵序流体 | 1D 重建 (s) | 权威 CSV (s) |
|---|---|---|
{chr(10).join(arrival_lines)}

结果：**{"全部一致" if not arrival_mismatches else "存在不符（不应发生，导出已拒绝）"}**。
### ② 碰压时刻 t_bump vs NPZ 最后快照
- 重建 t_bump = {f'{t_bump:.3f} s' if t_bump is not None else '未碰压（V 冻结在序列终点）'}
- NPZ 最后快照 = {last_snapshot_s:.3f} s
- 差值 = {f'{t_bump_delta_s:+.3f} s' if t_bump_delta_s is not None else '—'}（2D 模拟终点 ≡ 1D 碰压断面，量级秒级即视为一致）

## 产出文件
- `{long_csv_path.name}`：长格式，{n_long_rows} 行
  （= {n_t} 时刻 × ({n_in} 域内深度 × 4 通道 + {n_out} 域外深度 × 1 sharp 相)）
- `{{...}}_宽格式_{{lead,tail,spacer,mud}}.csv` ×4：行 = 全部 {n_depths} 表格深度，
  列 = {n_t} 个快照时刻（列名 = time_s 保留 3 位小数 + "_s"），无空值
"""
    doc_path.write_text(doc, encoding="utf-8")

    return {
        "well_label": well_label,
        "long_csv": str(long_csv_path),
        "wide_csvs": {ch: str(p) for ch, p in wide_paths.items()},
        "doc_md": str(doc_path),
        "n_long_rows": n_long_rows,
        "n_in_domain": n_in,
        "n_out_domain": n_out,
        "n_snapshots": n_t,
        "k_scale": k_scale,
        "pipe_volume_target_m3": float(pipe_volume_target),
        "s_tail_m3": float(s_tail),
        "v_bump_m3": float(v_bump),
        "v_frozen_m3": float(v_frozen),
        "t_bump_s": t_bump,
        "last_snapshot_s": last_snapshot_s,
        "t_bump_minus_last_snapshot_s": t_bump_delta_s,
        "arrival_table": [{"fluid": n, "reconstructed_s": t} for n, t in mine_rows],
        "arrival_mismatches": arrival_mismatches,
    }


# ── 全井导出：两井配置与入口 ─────────────────────────────────────

PROJECT_ROOT = Path(__file__).resolve().parents[2]  # cement model/

_FULL_WELL_CONFIGS = {
    "ht1_003": {
        "well_label": "呼1-003",
        "loader_module": "cemdisp.data.loaders.ht1_003_loader",
        "loader_func": "load_ht1_003_tailpipe",
        "loader_kwargs": {
            # loader 默认路径基于 parents[3] 指向仓库根，但参考数据实际在上一级
            # 控压固井项目/参考文档/，此处显式传参（几何/流变与 runner 同源）
            "reference_root": PROJECT_ROOT.parent / "参考文档" / "呼1-003" / "新",
            "caliper_csv_path": PROJECT_ROOT.parent / "参考文档" / "现场资料提取" / "ht1_003_呼1-003" / "caliper_profile.csv",
            "inclination_csv_path": PROJECT_ROOT.parent / "参考文档" / "现场资料提取" / "ht1_003_呼1-003" / "inclination_profile.csv",
        },
        "table_csv": r"D:\users\desktop\research\控压固井项目\参考文档\呼1-003\呼1-003井身结构.csv",
    },
    "ht1_004": {
        "well_label": "呼1-004",
        "loader_module": "cemdisp.data.loaders.ht1_004_loader",
        "loader_func": "load_ht1_004_tailpipe",
        "loader_kwargs": {
            "reference_root": PROJECT_ROOT.parent / "参考文档" / "呼1-004",
            "caliper_csv_path": PROJECT_ROOT.parent / "参考文档" / "现场资料提取" / "ht1_004_呼1-004" / "caliper_profile.csv",
            "inclination_csv_path": PROJECT_ROOT.parent / "参考文档" / "现场资料提取" / "ht1_004_呼1-004" / "inclination_profile.csv",
        },
        "table_csv": r"D:\users\desktop\research\控压固井项目\参考文档\呼1-004\呼1-004井身结构.csv",
    },
}


def run_fullwell_export(well_key: str) -> dict:
    """单井全井导出入口：从 loader 取泵序/管容（与 runner 同源），零重跑。"""
    import importlib as _importlib

    from cemdisp.transport1d.casing_flow import CasingFlowSolver

    if well_key not in _FULL_WELL_CONFIGS:
        raise KeyError(f"未配置的井: {well_key}，可选: {list(_FULL_WELL_CONFIGS)}")
    cfg = _FULL_WELL_CONFIGS[well_key]
    loader_fn = getattr(_importlib.import_module(cfg["loader_module"]), cfg["loader_func"])
    well_spec, fluids, schedule, _ = loader_fn(**cfg["loader_kwargs"])
    results_dir = PROJECT_ROOT / "results" / f"{cfg['well_label']}_1D2D耦合模型"
    npz_path = results_dir / f"{cfg['well_label']}_1D2D耦合模型_2D场数据.npz"
    shoe_csv = results_dir / f"{cfg['well_label']}_1D2D耦合模型_鞋口出流时序.csv"
    profile = [(float(p.depth_md_m), float(p.value)) for p in well_spec.pipe_id_profile]
    return export_fullwell_depth_time_shares(
        npz_path, cfg["table_csv"], results_dir, shoe_csv,
        well_label=cfg["well_label"], schedule=schedule, fluids=fluids,
        pipe_profile_depths=[d for d, _ in profile],
        pipe_profile_ids_mm=[i for _, i in profile],
        shoe_md_m=float(well_spec.shoe_md_m),
        # 1D 求解器管容（复用同款优先级链：shoe_lag > pipe_id_profile > liner_id）
        pipe_volume_target=float(well_spec.shoe_md_m)
        * CasingFlowSolver._pipe_cross_section_area(well_spec),
    )


# ── CLI ─────────────────────────────────────────────────────────
if __name__ == "__main__":
    import sys

    def _cli_flag(name: str, default: str | None = None) -> str | None:
        return sys.argv[sys.argv.index(name) + 1] if name in sys.argv else default

    if "--fullwell" in sys.argv:
        # 全井模式：python scripts/entrypoints/export_depth_time_concentration.py --fullwell --well ht1_003
        well_key = _cli_flag("--well")
        if well_key is None:
            print("用法: python export_depth_time_concentration.py --fullwell --well <ht1_003|ht1_004>")
            sys.exit(1)
        summary = run_fullwell_export(well_key)
        print(f"[全井导出] {summary['well_label']}: 长格式 {summary['n_long_rows']} 行 "
              f"(域内 {summary['n_in_domain']} / 域外 {summary['n_out_domain']} 深度点, "
              f"{summary['n_snapshots']} 快照)")
        print(f"  管容标定 k={summary['k_scale']:.6f}, Vp_shoe={summary['pipe_volume_target_m3']:.3f} m³, "
              f"V_bump={summary['v_bump_m3']:.3f} m³")
        if summary["t_bump_s"] is not None:
            print(f"  t_bump={summary['t_bump_s']:.3f}s vs 最后快照 {summary['last_snapshot_s']:.3f}s "
                  f"(Δ={summary['t_bump_minus_last_snapshot_s']:+.3f}s)")
        else:
            print("  未碰压：V 冻结在顶替序列终点")
        for row in summary["arrival_table"]:
            print(f"  到达鞋口: {row['fluid']} @ {row['reconstructed_s']}")
        print(f"  长格式: {summary['long_csv']}")
        print(f"  说明: {summary['doc_md']}")
    else:
        if len(sys.argv) < 2:
            print("用法: python export_depth_time_concentration.py <npz_path> [output_dir]")
            print("      python export_depth_time_concentration.py --fullwell --well <ht1_003|ht1_004>")
            print("示例: python export_depth_time_concentration.py results/呼1-004_1D2D耦合模型/呼1-004_1D2D耦合模型_2D场数据.npz")
            sys.exit(1)

        npz_file = Path(sys.argv[1])
        out_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else None
        if not npz_file.exists():
            print(f"错误: 文件不存在 - {npz_file}")
            sys.exit(1)

        export_from_npz(npz_file, out_dir)
