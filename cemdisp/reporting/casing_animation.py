# -*- coding: utf-8 -*-
"""管内段 + 环空段拼接的顶替过程可视化（2026-09-26 设计规格 §3 修订版）。

消费 :class:`~cemdisp.transport1d.casing_depth_profile.CasingDepthProfile`
（整井剖面 + 放大窗剖面各一份）与环空模拟结果，产出：

1. ``{井名}_管内环空拼接顶替动画.gif`` —— 管内全相着色 + 环空四通道 + 局部放大窗
2. ``{井名}_管内环空拼接快照.png``   —— 论文用多帧拼贴

布局口径（2026-09-26 用户裁定）：
- **遮挡修复**：每根色标占用独立 GridSpec 轴，不使用跨轴 ``fig.colorbar(ax=[...])``；
- **管内 = 全部流体**逐相追踪，按相名分配可区分颜色，混合区按份额做颜色线性混合；
- **环空 = 四通道**：领浆 / 尾浆 / 隔离液 / 钻井液（与导出脚本同口径）；
- **局部放大窗**：整井尺度上带宽仅占域长 0.3% 左右（纸面亚毫米），放大窗是
  唯一能让混浆带与带宽被看见、被检验的载体。

图内强制口径声明（Q14 = 审稿人 + 汇报）：径向均匀非求解结果、带宽由 α 定标、
胶塞面为工艺锐界面。
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from cemdisp.reporting.plots import _safe_filename_component, _setup_chinese_font
from cemdisp.transport1d.casing_depth_profile import (
    band_contract_violations,
    band_width_summary,
)

if TYPE_CHECKING:  # pragma: no cover - 仅类型检查
    from cemdisp.models2d.annulus_d2dga import AnnulusSimulationResult
    from cemdisp.transport1d.casing_depth_profile import CasingDepthProfile

LOGGER = logging.getLogger(__name__)

# 环空四通道（与 export_depth_time_concentration.py 的 _FULL_WELL_CHANNELS 同口径）
_ANNULUS_CHANNELS = ("领浆", "尾浆", "隔离液", "钻井液")
_ANNULUS_CMAPS = {"领浆": "YlGn", "尾浆": "inferno", "隔离液": "Blues", "钻井液": "BrBG"}

# 相名关键词 → 基色（RGB，0-1）。同名组多次出现时按次序做明度扰动以示区分。
_PHASE_KEYWORDS: tuple[tuple[tuple[str, ...], tuple[float, float, float]], ...] = (
    (("压塞液",), (0.85, 0.12, 0.15)),
    (("尾浆", "尾管水泥浆", "水泥"), (0.95, 0.55, 0.08)),
    (("领浆", "中间浆"), (0.97, 0.87, 0.30)),
    (("先导浆", "隔离液", "平衡液", "冲洗"), (0.16, 0.42, 0.82)),
    (("替钻井液", "顶替液"), (0.10, 0.70, 0.68)),
    (("保护液",), (0.58, 0.25, 0.75)),
    (("基液",), (0.20, 0.64, 0.28)),
    (("井浆", "钻井液", "泥浆"), (0.44, 0.35, 0.26)),
)
_PHASE_FALLBACK = (0.55, 0.55, 0.55)

_CAPTION_LINES = (
    "口径声明（不得省略）：",
    "① 管内为井筒竖条：横向=管内径，径向按均匀处理，非求解结果；",
    "② 混浆带宽度由标定参数 α=0.25 经 σ_t 决定，非独立物理预测；",
    "③ 红色实线=胶塞面（工艺锐界面：替浆顶胶塞驱动尾浆）；",
    "④ 放大窗自带深度分辨率，带宽可见性以放大窗为准。",
)


# ── 相配色 ────────────────────────────────────────────────────────────


def phase_colors(fluid_names) -> dict[str, tuple[float, float, float]]:
    """按相名分配可区分颜色；同名关键词组重复出现时做明度扰动。"""

    colors: dict[str, tuple[float, float, float]] = {}
    used: dict[str, int] = {}
    for name in fluid_names:
        base = _PHASE_FALLBACK
        key = "?"
        for keywords, rgb in _PHASE_KEYWORDS:
            if any(token in name for token in keywords):
                base, key = rgb, keywords[0]
                break
        seen = used.get(key, 0)
        used[key] = seen + 1
        if seen == 0:
            colors[name] = base
            continue
        # 第 2 次起向亮/暗交替偏移，保证同关键词组内颜色可区分
        shift = 0.22 * ((seen + 1) // 2)
        sign = 1.0 if seen % 2 == 1 else -1.0
        colors[name] = tuple(
            float(min(max(c + sign * shift * (1.0 - c if sign > 0 else c), 0.0), 1.0))
            for c in base
        )
    return colors


def blend_rgb(shares: np.ndarray, fluid_names, colors: dict) -> np.ndarray:
    """份额场 (n_z, n_f) → 混合 RGB (n_z, 3)；Σ份额=1 保证结果落在色域内。"""

    rgb = np.zeros((shares.shape[0], 3), dtype=float)
    for idx, name in enumerate(fluid_names):
        rgb += shares[:, idx : idx + 1] * np.asarray(colors[name], dtype=float)[None, :]
    return np.clip(rgb, 0.0, 1.0)


# ── 环空四通道 ────────────────────────────────────────────────────────


def annulus_channels(
    annulus_result: "AnnulusSimulationResult",
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """环空结果 → ({通道: (n_t,ny,nz)}, 深度 (nz,))，四通道与导出脚本同口径。"""

    def _stack(attr: str) -> np.ndarray | None:
        raw = getattr(annulus_result, attr, None)
        if raw is None or len(raw) == 0:
            return None
        return np.asarray(raw, dtype=float)

    lead, tail = _stack("lead_snapshots"), _stack("tail_snapshots")
    cement = _stack("cement_snapshots")
    if lead is not None and tail is not None:
        lead_arr, tail_arr = lead, tail
    elif cement is not None:
        lead_arr, tail_arr = np.zeros_like(cement), cement
    else:
        raise ValueError("环空结果既无 lead/tail 也无 cement 快照，无法拼接")

    spacer = _stack("spacer_snapshots")
    if spacer is None:
        spacer = np.zeros_like(lead_arr)
    channels = {
        "领浆": lead_arr,
        "尾浆": tail_arr,
        "隔离液": spacer,
        "钻井液": np.clip(1.0 - lead_arr - tail_arr - spacer, 0.0, 1.0),
    }
    md = np.asarray(annulus_result.geom["md"], dtype=float)
    return channels, md


# ── 几何与标记 ────────────────────────────────────────────────────────


def _casing_strip(rgb: np.ndarray) -> np.ndarray:
    """(n_z, 3) 混合色 → (n_z, 2, 3) 竖条，使井筒在图上具有宽度。"""

    return np.repeat(np.asarray(rgb, dtype=float)[:, np.newaxis, :], 2, axis=1)


def _casing_bore_extent(profile: "CasingDepthProfile") -> tuple[float, float]:
    """管内面板横轴范围：以平均管内半径为准，对应真实管内径。"""

    area = np.asarray(profile.pipe_area_m2, dtype=float)
    radius = float(np.sqrt(max(float(np.mean(area)), 0.0) / math.pi))
    return (-radius, radius)


def _front_lines(profile: "CasingDepthProfile", time_idx: int) -> list[tuple[float, str, bool]]:
    """当前帧各界面位置 → [(深度 m, 标签, 是否胶塞面), ...]。

    界面位置 = **已到达的最深**深度。到达时刻 τ(z) 随深度单调递增，故在时刻 t
    已通过的深度集合是浅端的一个前缀，前缘取该前缀最深处（``flatnonzero[-1]``）。
    用 ``argmax(crossed)`` 会取到第 0 号格（井口）—— 2026-09-26 实测缺陷。
    """

    depths = np.asarray(profile.depths_m, dtype=float)
    lines: list[tuple[float, str, bool]] = []
    for k in range(profile.boundary_arrival_s.shape[0]):
        arrival = profile.boundary_arrival_s[k]
        crossed = np.flatnonzero(np.isfinite(arrival) & (arrival <= profile.times_s[time_idx]))
        if crossed.size == 0:
            continue
        label = f"{profile.fluid_names[k]} → {profile.fluid_names[k + 1]}"
        lines.append((float(depths[int(crossed[-1])]), label, bool(profile.plug_boundary[k])))
    return lines


def _zoom_center(zoom_profile: "CasingDepthProfile") -> float:
    depths = np.asarray(zoom_profile.depths_m, dtype=float)
    return float(0.5 * (depths[0] + depths[-1]))


def _nearest_boundary(profile: "CasingDepthProfile", depth_m: float) -> int | None:
    """离给定深度最近的「已出现」界面下标（用于放大窗标注带宽）。

    界面位置同样取已到达的最深处（与 :func:`_front_lines` 同口径）。
    """

    best, best_gap = None, float("inf")
    for k in range(profile.boundary_arrival_s.shape[0]):
        reached = np.flatnonzero(np.isfinite(profile.boundary_arrival_s[k]))
        if reached.size == 0:
            continue
        gap = abs(float(profile.depths_m[int(reached[-1])]) - depth_m)
        if gap < best_gap:
            best, best_gap = k, gap
    return best


def _apply_front_marks(ax, profile: "CasingDepthProfile", time_idx: int) -> None:
    """叠加前缘轨迹、胶塞面与数值锐化标记（Q8-① / Q16-C）。"""

    for line in list(ax.lines):
        line.remove()
    for text in list(ax.texts):
        text.remove()
    depths = np.asarray(profile.depths_m, dtype=float)
    cell = float(abs(depths[1] - depths[0])) if depths.size > 1 else 1.0
    plug_depths: list[float] = []
    for depth, _label, is_plug in _front_lines(profile, time_idx):
        if is_plug:
            ax.axhline(depth, color="crimson", lw=2.0)
            plug_depths.append(depth)
        else:
            ax.axhline(depth, color="white", lw=0.9, linestyle="--")
    if plug_depths:
        ax.text(0.03, plug_depths[0], " 胶塞面（工艺）", color="crimson", fontsize=7,
                va="bottom", ha="left", transform=ax.get_yaxis_transform())

    sharp = profile.sharp_cells()
    if sharp.shape[1] > 0 and time_idx < sharp.shape[0]:
        for cell_idx in np.flatnonzero(sharp[time_idx]):
            depth = float(depths[cell_idx])
            if any(abs(depth - plug_depth) <= 2.0 * cell for plug_depth in plug_depths):
                continue
            ax.axhline(depth, color="orange", lw=0.7, linestyle=":")
        if np.any(sharp[time_idx]):
            ax.text(0.97, 0.02, "橙点线=数值锐化", color="orange", fontsize=7,
                    va="bottom", ha="right", transform=ax.transAxes)


# ── 图幅 ─────────────────────────────────────────────────────────────


def _build_figure(
    profile: "CasingDepthProfile",
    zoom_profile: "CasingDepthProfile",
    annulus_result: "AnnulusSimulationResult",
    *,
    figsize: tuple[float, float] = (17.0, 10.0),
):
    """搭建图幅：整井 5 面板 + 独立色标行 + 放大窗与声明面板 + 相名图例。

    色标各占独立 GridSpec 轴（第 1 行），与任何数据轴都不重叠——这是 2026-09-26
    「右侧色标遮挡图表」的修复口径。
    """

    channels, md = annulus_channels(annulus_result)
    depths = np.asarray(profile.depths_m, dtype=float)
    depth_min, depth_max = float(depths[0]), float(depths[-1])
    x_lo, x_hi = _casing_bore_extent(profile)
    colors = phase_colors(profile.fluid_names)

    fig = plt.figure(figsize=figsize)
    grid = fig.add_gridspec(3, 5, height_ratios=(3.4, 0.16, 1.25), hspace=0.30, wspace=0.20)
    casing_ax = fig.add_subplot(grid[0, 0])
    ann_axes = [fig.add_subplot(grid[0, i], sharey=casing_ax) for i in range(1, 5)]
    cax_list = [fig.add_subplot(grid[1, i]) for i in range(1, 5)]
    zoom_ax = fig.add_subplot(grid[2, 0:3])
    info_ax = fig.add_subplot(grid[2, 3:5])

    casing_ax.imshow(
        _casing_strip(blend_rgb(profile.shares[0], profile.fluid_names, colors)),
        aspect="auto", extent=(x_lo, x_hi, depth_max, depth_min), origin="upper",
    )
    casing_ax.set_xlabel("管内（井筒竖条）")
    casing_ax.set_ylabel("井深 / m")
    casing_ax.set_title("管内 全相组成", fontsize=12, fontweight="bold")
    casing_ax.set_ylim(depth_max, depth_min)

    ann_images: dict[str, object] = {}
    for ax, cax, name in zip(ann_axes, cax_list, _ANNULUS_CHANNELS):
        image = ax.imshow(
            np.flipud(channels[name][0].T), vmin=0.0, vmax=1.0, cmap=_ANNULUS_CMAPS[name],
            aspect="auto", extent=(0.0, 1.0, float(np.max(md)), float(np.min(md))),
            origin="upper",
        )
        ann_images[name] = image
        ax.set_xlabel("方位角（宽边→窄边）")
        ax.set_title(f"环空 {name}", fontsize=12, fontweight="bold")
        ax.set_ylim(depth_max, depth_min)
        # 每根色标独占一行轴，杜绝遮挡
        fig.colorbar(image, cax=cax, orientation="horizontal").set_label(
            f"{name} 浓度", fontsize=8.5)

    zd = np.asarray(zoom_profile.depths_m, dtype=float)
    for ax in (casing_ax, *ann_axes):
        ax.axhspan(float(zd[0]), float(zd[-1]), facecolor="none",
                   edgecolor="darkorange", lw=1.2, zorder=5)

    handles = [Patch(facecolor=colors[n], label=n) for n in profile.fluid_names]
    handles.append(Line2D([], [], color="crimson", lw=2.0, label="胶塞面（工艺）"))
    handles.append(Line2D([], [], color="orange", lw=0.9, linestyle=":", label="数值锐化"))
    fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=8.5,
               frameon=False, bbox_to_anchor=(0.5, 0.002))

    zoom_ax.set_xlabel("相份额")
    zoom_ax.set_ylabel("井深 / m")
    zoom_ax.set_xlim(0.0, 1.0)
    zoom_ax.set_ylim(float(zd[-1]), float(zd[0]))
    zoom_ax.grid(alpha=0.25, linestyle=":")

    declared = band_contract_violations(
        zoom_profile, cell_size_m=float(np.median(np.diff(zd)))
    )
    summary = band_width_summary(profile)
    info_ax.axis("off")
    info_ax.text(0.0, 1.0, "\n".join(_CAPTION_LINES), va="top", ha="left",
                 fontsize=8.5, color="0.2", transform=info_ax.transAxes)
    info_ax.text(
        0.0, 0.28,
        f"整井尺度带宽诊断（报告，不拦截）：\n"
        f"  界面数 {summary['n_interfaces']}，中位带宽 {summary['median_band_m']:.1f} m"
        f" = 域长 {float(summary['median_domain_frac']) * 100:.2f}%\n"
        f"  即在整井纵轴上约 {6.5 * 25.4 * float(summary['median_domain_frac']):.2f} mm，"
        f"肉眼不可见——故须以放大窗为准\n"
        f"放大窗尺度契约（通过条件，3 格）："
        + ("通过" if not declared else f"违规 {len(declared)} 条"),
        va="top", ha="left", fontsize=8.5, color="0.2", transform=info_ax.transAxes,
    )

    return {
        "fig": fig,
        "casing_ax": casing_ax,
        "ann_axes": ann_axes,
        "ann_images": ann_images,
        "channels": channels,
        "zoom_ax": zoom_ax,
        "zoom_center": _zoom_center(zoom_profile),
        "colors": colors,
        "extent": (x_lo, x_hi, depth_max, depth_min),
    }


def _draw_zoom(
    parts: dict,
    profile: "CasingDepthProfile",
    zoom_profile: "CasingDepthProfile",
    time_idx: int,
) -> None:
    """在放大窗里画当前帧的相份额—深度曲线，并标注最近界面的混浆带宽度。"""

    ax = parts["zoom_ax"]
    for line in list(ax.lines):
        line.remove()
    for text in list(ax.texts):
        text.remove()
    for collection in list(ax.collections):
        collection.remove()

    colors = parts["colors"]
    shares = zoom_profile.shares[time_idx]
    for idx, name in enumerate(zoom_profile.fluid_names):
        column = shares[:, idx]
        if float(np.max(column)) <= 1.0e-6:
            continue
        ax.plot(column, zoom_profile.depths_m, color=colors[name], lw=1.6, label=name)

    k = _nearest_boundary(zoom_profile, parts["zoom_center"])
    note = ""
    if k is not None:
        active = zoom_profile.sigma_t_s[k] > 0.0
        if np.any(active):
            d = zoom_profile.depths_m[active]
            ax.axhspan(float(np.min(d)), float(np.max(d)), color="darkorange",
                       alpha=0.10, zorder=0)
        is_plug = bool(zoom_profile.plug_boundary[k])
        band = float(np.max(zoom_profile.band_width_m[k]))
        note = (
            f"最近界面：{zoom_profile.fluid_names[k]} → {zoom_profile.fluid_names[k + 1]}"
            + ("（胶塞面：工艺锐界面，σ≡0）" if is_plug else f"，带宽 ≈ {band:.1f} m")
        )
    ax.set_title(
        f"放大窗 {float(zoom_profile.depths_m[0]):.0f}–{float(zoom_profile.depths_m[-1]):.0f} m"
        f" — t = {profile.times_s[time_idx] / 60.0:.1f} min\n{note}",
        fontsize=10,
    )
    ax.set_ylim(float(zoom_profile.depths_m[-1]), float(zoom_profile.depths_m[0]))


# ── 对外入口 ─────────────────────────────────────────────────────────


def animate_stitched_displacement(
    profile: "CasingDepthProfile",
    zoom_profile: "CasingDepthProfile",
    annulus_result: "AnnulusSimulationResult",
    output_dir: Path | str,
    *,
    well_name: str | None = None,
    interval_ms: int = 200,
    fps: int = 8,
) -> Path:
    """生成「管内（全相）+ 环空（四通道）+ 局部放大窗」拼接的顶替过程 GIF。

    ``profile`` 与 ``zoom_profile`` 必须用同一组时刻构建（= 环空
    ``snapshot_times_s``），使三处同帧（Q22-b-i）；放大窗可另用更细的深度网格。

    Raises:
        ValueError: 三侧帧数不一致。
    """

    _setup_chinese_font()
    channels, _md = annulus_channels(annulus_result)
    n_frames = len(profile.times_s)
    n_annulus = channels["尾浆"].shape[0]
    if n_annulus != n_frames:
        raise ValueError(f"管内剖面帧数 {n_frames} 与环空快照帧数 {n_annulus} 不一致")
    if len(zoom_profile.times_s) != n_frames:
        raise ValueError(
            f"放大窗帧数 {len(zoom_profile.times_s)} 与整井剖面帧数 {n_frames} 不一致"
        )

    parts = _build_figure(profile, zoom_profile, annulus_result)
    fig = parts["fig"]
    colors = parts["colors"]
    x_lo, x_hi, depth_max, depth_min = parts["extent"]
    title = fig.suptitle("", fontsize=14, fontweight="bold")

    def update_frame(frame_idx: int):
        rgb = blend_rgb(profile.shares[frame_idx], profile.fluid_names, colors)
        parts["casing_ax"].images[0].set_data(_casing_strip(rgb))
        parts["casing_ax"].images[0].set_extent((x_lo, x_hi, depth_max, depth_min))
        _apply_front_marks(parts["casing_ax"], profile, frame_idx)
        for name, ax in zip(_ANNULUS_CHANNELS, parts["ann_axes"]):
            parts["ann_images"][name].set_data(np.flipud(channels[name][frame_idx].T))
        _draw_zoom(parts, profile, zoom_profile, frame_idx)
        title.set_text(
            f"管内—环空顶替过程拼接 — t = {profile.times_s[frame_idx] / 60.0:.1f} min"
        )
        return ()

    animation = FuncAnimation(fig, update_frame, frames=n_frames,
                              interval=interval_ms, blit=False)
    target = Path(output_dir)
    target.mkdir(parents=True, exist_ok=True)
    label = well_name or annulus_result.well_name
    path = target / f"{_safe_filename_component(label)}_管内环空拼接顶替动画.gif"
    fig.subplots_adjust(top=0.93, bottom=0.16, left=0.05, right=0.985)
    try:
        animation.save(path, writer="pillow", fps=fps)
    finally:
        plt.close(fig)
    LOGGER.info("拼接动画已写出：%s", path)
    return path


def plot_stitched_snapshots(
    profile: "CasingDepthProfile",
    zoom_profile: "CasingDepthProfile",
    annulus_result: "AnnulusSimulationResult",
    output_dir: Path | str,
    *,
    well_name: str | None = None,
    n_panels: int = 4,
) -> Path:
    """生成论文用多帧拼贴：每帧一列「管内全相竖条 + 该时刻放大窗剖面」。"""

    _setup_chinese_font()
    channels, _md = annulus_channels(annulus_result)
    n_frames = len(profile.times_s)
    if channels["尾浆"].shape[0] != n_frames or len(zoom_profile.times_s) != n_frames:
        raise ValueError("管内 / 放大窗 / 环空 帧数不一致，无法拼贴")
    indices = np.unique(np.linspace(0, n_frames - 1, min(n_panels, n_frames)).astype(int))

    depths = np.asarray(profile.depths_m, dtype=float)
    depth_min, depth_max = float(depths[0]), float(depths[-1])
    x_lo, x_hi = _casing_bore_extent(profile)
    colors = phase_colors(profile.fluid_names)
    zd = np.asarray(zoom_profile.depths_m, dtype=float)

    fig, axes = plt.subplots(2, len(indices), figsize=(3.0 * len(indices), 9.0),
                             squeeze=False,
                             gridspec_kw={"height_ratios": (2.4, 1.0), "hspace": 0.35})
    for col, frame_idx in enumerate(indices):
        top, bottom = axes[0][col], axes[1][col]
        top.imshow(_casing_strip(blend_rgb(profile.shares[frame_idx],
                                           profile.fluid_names, colors)),
                   aspect="auto", extent=(x_lo, x_hi, depth_max, depth_min), origin="upper")
        top.set_ylim(depth_max, depth_min)
        top.axhspan(float(zd[0]), float(zd[-1]), facecolor="none",
                    edgecolor="darkorange", lw=1.2, zorder=5)
        _apply_front_marks(top, profile, int(frame_idx))
        top.set_title(f"t = {profile.times_s[frame_idx] / 60.0:.0f} min", fontsize=10)
        top.set_xlabel("管内（井筒竖条）")
        if col == 0:
            top.set_ylabel("井深 / m")

        for idx, name in enumerate(zoom_profile.fluid_names):
            column = zoom_profile.shares[frame_idx][:, idx]
            if float(np.max(column)) <= 1.0e-6:
                continue
            bottom.plot(column, zoom_profile.depths_m, color=colors[name], lw=1.4, label=name)
        bottom.set_xlim(0.0, 1.0)
        bottom.set_ylim(float(zd[-1]), float(zd[0]))
        bottom.grid(alpha=0.25, linestyle=":")
        bottom.set_xlabel("相份额")
        if col == 0:
            bottom.set_ylabel("放大窗井深 / m")

    handles = [Patch(facecolor=colors[n], label=n) for n in profile.fluid_names]
    handles.append(Line2D([], [], color="crimson", lw=2.0, label="胶塞面（工艺）"))
    # 相名图例最多 3 行，底部须留足空间，避免与数据轴或声明文字相互遮挡
    fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=8.5,
               frameon=False, bbox_to_anchor=(0.5, 0.003))
    fig.suptitle(
        f"管内全相组成与放大窗剖面 — {well_name or annulus_result.well_name}",
        fontsize=13, fontweight="bold",
    )
    fig.text(0.5, 0.115, " ｜ ".join(_CAPTION_LINES[1:4]), ha="center", va="bottom",
             fontsize=7.5, color="0.3")

    target = Path(output_dir)
    target.mkdir(parents=True, exist_ok=True)
    label = well_name or annulus_result.well_name
    path = target / f"{_safe_filename_component(label)}_管内环空拼接快照.png"
    fig.subplots_adjust(top=0.90, bottom=0.20, left=0.06, right=0.99)
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path


__all__ = [
    "animate_stitched_displacement",
    "annulus_channels",
    "blend_rgb",
    "phase_colors",
    "plot_stitched_snapshots",
]
