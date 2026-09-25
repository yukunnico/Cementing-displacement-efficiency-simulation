# -*- coding: utf-8 -*-
"""管内段 + 环空段拼接的顶替过程可视化（2026-09-26 设计规格 §3）。

本模块消费 :class:`~cemdisp.transport1d.casing_depth_profile.CasingDepthProfile`
与 ``AnnulusSimulationResult``，产出：

1. ``{井名}_管内环空拼接顶替动画.gif`` —— 浓度场 + 前缘轨迹 + 胶塞面标记（Q8-①/③）
2. ``{井名}_管内环空拼接快照.png``   —— 论文用多帧拼贴（Q4-A）
3. ``{井名}_管内段速度剪切率动画.gif`` —— 同一批数据二次渲染（Q8-④ / Q20-4）

口径（不得在图中隐去，Q14 = 审稿人 + 汇报）：
- 管内面板为**井筒竖条**，横向为管内径，**径向均匀**，不是求解结果（Q18-C/A）；
- 管内混浆带宽度由标定参数 ``dispersion_alpha``（默认 0.25）经 ``σ_t`` 决定，
  **不是独立物理预测**；
- 胶塞面为**工艺锐界面**（替浆顶胶塞驱动尾浆），是混合带口径的显式例外（Q16）；
- 环空面板沿用现有约定：横轴归一化方位角（宽边→窄边）、纵轴井深（浅→深）。
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

from cemdisp.reporting.plots import _safe_filename_component, _setup_chinese_font

if TYPE_CHECKING:  # pragma: no cover - 仅类型检查
    from cemdisp.models2d.annulus_d2dga import AnnulusSimulationResult
    from cemdisp.transport1d.casing_depth_profile import CasingDepthProfile

LOGGER = logging.getLogger(__name__)

_CMAP_CEMENT = "viridis"
_CMAP_SPACER = "coolwarm"

# 图中必须出现的口径声明（Q18/Q14：审稿人一定会问）
_CAPTION_RADIAL = "管内为井筒竖条：横向=管内径，径向按均匀处理（非求解结果）"
_CAPTION_BAND = "混浆带宽度由标定参数 α=0.25 经 σ_t 决定，非独立物理预测"
_CAPTION_PLUG = "红色实线=胶塞面（工艺锐界面：替浆顶胶塞驱动尾浆）"
_CAPTION_VELOCITY = "活塞流口径：全管截面同速 U(t,z)=Q(t)/A(z)，径向均匀；剪切率取 8U/(2R)"


def _casing_strip(field: np.ndarray) -> np.ndarray:
    """把 (n_z,) 的管内剖面摊成 (n_z, 2) 竖条，供 imshow 画出有宽度的管柱。"""

    return np.repeat(np.asarray(field, dtype=float)[:, np.newaxis], 2, axis=1)


def _casing_bore_extent(profile: "CasingDepthProfile") -> tuple[float, float]:
    """管内面板横轴范围：以平均管内半径为准，使竖条宽度对应真实管内径。"""

    area = np.asarray(profile.pipe_area_m2, dtype=float)
    radius = float(np.sqrt(max(float(np.mean(area)), 0.0) / np.pi))
    return (-radius, radius)


def _front_lines(profile: "CasingDepthProfile", time_idx: int) -> list[tuple[float, str, bool]]:
    """当前帧各界面位置 → [(深度 m, 标签, 是否胶塞面), ...]。

    界面位置取已越过的最深格点（即 frac ≥ 0.5 的深度端），未到达则跳过。
    """

    depths = np.asarray(profile.depths_m, dtype=float)
    lines: list[tuple[float, str, bool]] = []
    for k in range(profile.boundary_arrival_s.shape[0]):
        arrival = profile.boundary_arrival_s[k]
        crossed = np.isfinite(arrival) & (arrival <= profile.times_s[time_idx])
        if not np.any(crossed):
            continue
        label = f"{profile.fluid_names[k]} → {profile.fluid_names[k + 1]}"
        lines.append((float(depths[int(np.argmax(crossed))]), label, bool(profile.plug_boundary[k])))
    return lines


def _annulus_profiles(
    annulus_result: "AnnulusSimulationResult",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """环空快照 → (水泥 (n_t,ny,nz), 隔离液 (n_t,ny,nz), 深度 (nz,))。"""

    cement = np.asarray(annulus_result.cement_snapshots, dtype=float)
    if cement.size == 0:
        raise ValueError("环空结果不含浓度场快照，无法拼接")
    raw_spacer = getattr(annulus_result, "spacer_snapshots", ())
    spacer = (
        np.asarray(raw_spacer, dtype=float)
        if raw_spacer is not None and len(raw_spacer) > 0
        else np.zeros_like(cement)
    )
    md = np.asarray(annulus_result.geom["md"], dtype=float)
    return cement, spacer, md


def _cement_channel(profile: "CasingDepthProfile") -> str:
    """管内“水泥”代表通道：领浆/尾浆中份额峰值更大者。"""

    cement_names = [n for n in profile.fluid_names if "领浆" in n or "尾浆" in n or "水泥" in n]
    if not cement_names:
        return profile.fluid_names[-1]
    return max(cement_names, key=lambda n: float(np.max(profile.channel_shares(n))))


def _spacer_channel(profile: "CasingDepthProfile") -> str:
    """管内“隔离液”通道：名称含隔离液/平衡液/先导浆/冲洗者。"""

    for name in profile.fluid_names:
        if any(key in name for key in ("隔离液", "平衡液", "先导浆", "冲洗")):
            return name
    raise ValueError("管内剖面中没有可识别的隔离液通道，无法拼接出图")


def _apply_front_marks(ax, profile: "CasingDepthProfile", time_idx: int) -> None:
    """在图中叠加前缘轨迹与胶塞面标记，并**同时标出数值锐化**（Q16-C）。"""

    for line in list(ax.lines):
        line.remove()
    for text in list(ax.texts):
        text.remove()
    for depth, _label, is_plug in _front_lines(profile, time_idx):
        if is_plug:
            ax.axhline(depth, color="crimson", lw=2.0)
            ax.text(0.02, depth, " 胶塞面（工艺）", color="crimson", fontsize=7,
                    va="bottom", ha="left", transform=ax.get_yaxis_transform())
        else:
            ax.axhline(depth, color="white", lw=0.9, linestyle="--")

    # 数值锐化：带模式下去检测相邻深度格 ≥ _SHARP_JUMP 的跳变，且不落在胶塞面
    sharp = profile.sharp_cells()
    if sharp.shape[1] > 0 and time_idx < sharp.shape[0]:
        plug_depths = [d for d, _l, p in _front_lines(profile, time_idx) if p]
        depths = np.asarray(profile.depths_m, dtype=float)
        for cell in np.flatnonzero(sharp[time_idx]):
            depth = float(depths[cell])
            if any(abs(depth - pd) <= 2.0 * abs(depths[1] - depths[0]) for pd in plug_depths):
                continue  # 已由胶塞面标记覆盖
            ax.axhline(depth, color="orange", lw=0.7, linestyle=":")
        if np.any(sharp[time_idx]):
            ax.text(0.98, 0.02, "橙点线=数值锐化", color="orange", fontsize=7,
                    va="bottom", ha="right", transform=ax.transAxes)


def _stitched_figure(
    profile: "CasingDepthProfile",
    annulus_result: "AnnulusSimulationResult",
    *,
    figsize: tuple[float, float] = (13.0, 9.0),
):
    """搭建四联图：管内水泥 / 管内隔离液 / 环空水泥 / 环空隔离液，共用深度轴。"""

    cement_c, spacer_c, md = _annulus_profiles(annulus_result)
    depths = np.asarray(profile.depths_m, dtype=float)
    depth_min, depth_max = float(depths[0]), float(depths[-1])
    x_lo, x_hi = _casing_bore_extent(profile)

    fig = plt.figure(figsize=figsize)
    grid = fig.add_gridspec(1, 4, width_ratios=(1.6, 1.6, 9.0, 9.0), wspace=0.12)
    axes = [fig.add_subplot(grid[0, i]) for i in range(4)]
    cement_ax, spacer_ax, ann_cement_ax, ann_spacer_ax = axes

    cement_ax.imshow(
        _casing_strip(profile.channel_shares(_cement_channel(profile))[0]),
        vmin=0.0, vmax=1.0, cmap=_CMAP_CEMENT, aspect="auto",
        extent=(x_lo, x_hi, depth_max, depth_min), origin="upper",
    )
    spacer_ax.imshow(
        _casing_strip(profile.channel_shares(_spacer_channel(profile))[0]),
        vmin=0.0, vmax=1.0, cmap=_CMAP_SPACER, aspect="auto",
        extent=(x_lo, x_hi, depth_max, depth_min), origin="upper",
    )
    ann_cement_image = ann_cement_ax.imshow(
        np.flipud(cement_c[0].T), vmin=0.0, vmax=1.0, cmap=_CMAP_CEMENT, aspect="auto",
        extent=(0.0, 1.0, float(np.max(md)), float(np.min(md))), origin="upper",
    )
    ann_spacer_image = ann_spacer_ax.imshow(
        np.flipud(spacer_c[0].T), vmin=0.0, vmax=1.0, cmap=_CMAP_SPACER, aspect="auto",
        extent=(0.0, 1.0, float(np.max(md)), float(np.min(md))), origin="upper",
    )

    for ax in axes:
        ax.set_ylim(depth_max, depth_min)  # 共用真实深度轴；环空域外自然留白
    for ax in (cement_ax, spacer_ax):
        ax.set_xlabel("管内（井筒竖条）")
    ann_cement_ax.set_xlabel("方位角（宽边→窄边）")
    ann_spacer_ax.set_xlabel("方位角（宽边→窄边）")
    cement_ax.set_ylabel("井深 / m")
    cement_ax.set_title("管内 水泥浓度", fontsize=12, fontweight="bold")
    spacer_ax.set_title("管内 隔离液浓度", fontsize=12, fontweight="bold")
    ann_cement_ax.set_title("环空 水泥浓度", fontsize=12, fontweight="bold")
    ann_spacer_ax.set_title("环空 隔离液浓度", fontsize=12, fontweight="bold")

    fig.colorbar(ann_cement_image, ax=[cement_ax, ann_cement_ax], shrink=0.8).set_label("水泥浓度")
    fig.colorbar(ann_spacer_image, ax=[spacer_ax, ann_spacer_ax], shrink=0.8).set_label("隔离液浓度")

    fig.text(0.5, 0.012, f"{_CAPTION_RADIAL}\n{_CAPTION_BAND}\n{_CAPTION_PLUG}",
             ha="center", va="bottom", fontsize=8.5, color="0.25")

    return fig, axes, ann_cement_image, ann_spacer_image, (cement_c, spacer_c)


def animate_stitched_displacement(
    profile: "CasingDepthProfile",
    annulus_result: "AnnulusSimulationResult",
    output_dir: Path | str,
    *,
    well_name: str | None = None,
    interval_ms: int = 200,
    fps: int = 8,
) -> Path:
    """生成「管内 + 环空」拼接的顶替过程 GIF。

    时间栅格取 ``profile.times_s``，调用方须先用环空 ``snapshot_times_s`` 构建管内
    剖面，使两侧同帧（Q22-b-i）。

    Raises:
        ValueError: 管内剖面帧数与环空快照帧数不一致。
    """

    _setup_chinese_font()
    cement_c, _spacer_c, _md = _annulus_profiles(annulus_result)
    n_casing = len(profile.times_s)
    n_annulus = cement_c.shape[0]
    if n_casing != n_annulus:
        raise ValueError(
            f"管内剖面帧数 {n_casing} 与环空快照帧数 {n_annulus} 不一致，"
            "请改用环空 snapshot_times_s 构建管内剖面（Q22-b-i）"
        )

    fig, axes, ann_cement_image, ann_spacer_image, (cem, spa) = _stitched_figure(
        profile, annulus_result
    )
    cement_name = _cement_channel(profile)
    spacer_name = _spacer_channel(profile)
    depth_max, depth_min = float(profile.depths_m[-1]), float(profile.depths_m[0])
    x_lo, x_hi = _casing_bore_extent(profile)
    title = fig.suptitle("", fontsize=14, fontweight="bold")

    def update_frame(frame_idx: int):
        for ax, name in ((axes[0], cement_name), (axes[1], spacer_name)):
            ax.images[0].set_data(_casing_strip(profile.channel_shares(name)[frame_idx]))
            ax.images[0].set_extent((x_lo, x_hi, depth_max, depth_min))
            _apply_front_marks(ax, profile, frame_idx)
        ann_cement_image.set_data(np.flipud(cem[frame_idx].T))
        ann_spacer_image.set_data(np.flipud(spa[frame_idx].T))
        title.set_text(
            f"管内—环空顶替过程拼接 — t = {profile.times_s[frame_idx] / 60.0:.1f} min"
        )
        return ()

    animation = FuncAnimation(fig, update_frame, frames=n_casing,
                              interval=interval_ms, blit=False)
    target = Path(output_dir)
    target.mkdir(parents=True, exist_ok=True)
    label = well_name or annulus_result.well_name
    path = target / f"{_safe_filename_component(label)}_管内环空拼接顶替动画.gif"
    # 本图含 gridspec + 跨轴 colorbar，tight_layout 不兼容，改用显式边距
    fig.subplots_adjust(top=0.92, bottom=0.13, left=0.06, right=0.985)
    try:
        animation.save(path, writer="pillow", fps=fps)
    finally:
        plt.close(fig)
    LOGGER.info("拼接动画已写出：%s", path)
    return path


def plot_stitched_snapshots(
    profile: "CasingDepthProfile",
    annulus_result: "AnnulusSimulationResult",
    output_dir: Path | str,
    *,
    well_name: str | None = None,
    n_panels: int = 6,
) -> Path:
    """生成论文用多帧拼贴 PNG（静帧；如需矢量可在此基础上另存 SVG/PDF）。"""

    _setup_chinese_font()
    cement_c, _spacer_c, _md = _annulus_profiles(annulus_result)
    n_frames = len(profile.times_s)
    if n_frames != cement_c.shape[0]:
        raise ValueError("管内剖面帧数与环空快照帧数不一致，无法拼贴")
    indices = np.unique(np.linspace(0, n_frames - 1, min(n_panels, n_frames)).astype(int))

    depths = np.asarray(profile.depths_m, dtype=float)
    depth_min, depth_max = float(depths[0]), float(depths[-1])
    x_lo, x_hi = _casing_bore_extent(profile)
    cement_name = _cement_channel(profile)

    fig, axes = plt.subplots(1, len(indices), figsize=(2.2 * len(indices), 7.0),
                             sharey=True, squeeze=False)
    column = axes[0]
    for ax, frame_idx in zip(column, indices):
        ax.imshow(_casing_strip(profile.channel_shares(cement_name)[frame_idx]),
                  vmin=0.0, vmax=1.0, cmap=_CMAP_CEMENT, aspect="auto",
                  extent=(x_lo, x_hi, depth_max, depth_min), origin="upper")
        ax.set_ylim(depth_max, depth_min)
        _apply_front_marks(ax, profile, int(frame_idx))
        ax.set_title(f"t = {profile.times_s[frame_idx] / 60.0:.0f} min", fontsize=10)
        ax.set_xlabel("管内")
    column[0].set_ylabel("井深 / m")
    fig.suptitle(
        f"管内段水泥浓度演化（井筒竖条）— {well_name or annulus_result.well_name}",
        fontsize=13, fontweight="bold",
    )
    fig.text(0.5, 0.008, f"{_CAPTION_RADIAL}\n{_CAPTION_BAND}\n{_CAPTION_PLUG}",
             ha="center", va="bottom", fontsize=8.0, color="0.25")

    target = Path(output_dir)
    target.mkdir(parents=True, exist_ok=True)
    label = well_name or annulus_result.well_name
    path = target / f"{_safe_filename_component(label)}_管内环空拼接快照.png"
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path


def animate_casing_velocity(
    profile: "CasingDepthProfile",
    output_dir: Path | str,
    *,
    well_name: str,
    interval_ms: int = 200,
    fps: int = 8,
) -> Path:
    """管内段速度 / 壁面剪切率二次渲染动画（Q8-④ / Q20-4，同一批数据）。"""

    _setup_chinese_font()
    depths = np.asarray(profile.depths_m, dtype=float)
    depth_min, depth_max = float(depths[0]), float(depths[-1])
    x_lo, x_hi = _casing_bore_extent(profile)
    v_max = max(float(profile.velocity_m_s.max()), 1.0e-9)
    s_max = max(float(profile.shear_rate_s.max()), 1.0e-9)

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 9.0))
    v_img = axes[0].imshow(_casing_strip(profile.velocity_m_s[0]),
                           vmin=0.0, vmax=v_max, cmap="magma", aspect="auto",
                           extent=(x_lo, x_hi, depth_max, depth_min), origin="upper")
    s_img = axes[1].imshow(_casing_strip(profile.shear_rate_s[0]),
                           vmin=0.0, vmax=s_max, cmap="cividis", aspect="auto",
                           extent=(x_lo, x_hi, depth_max, depth_min), origin="upper")
    fig.colorbar(v_img, ax=axes[0], shrink=0.8).set_label("截面平均流速 / (m/s)")
    fig.colorbar(s_img, ax=axes[1], shrink=0.8).set_label("壁面剪切率 / (1/s)")
    axes[0].set_title("管内流速", fontsize=12, fontweight="bold")
    axes[1].set_title("管内壁面剪切率", fontsize=12, fontweight="bold")
    for ax in axes:
        ax.set_xlabel("管内（井筒竖条）")
        ax.set_ylabel("井深 / m")
    fig.text(0.5, 0.012, _CAPTION_VELOCITY, ha="center", va="bottom",
             fontsize=8.5, color="0.25")
    title = fig.suptitle("", fontsize=14, fontweight="bold")

    def update_frame(frame_idx: int):
        v_img.set_data(_casing_strip(profile.velocity_m_s[frame_idx]))
        s_img.set_data(_casing_strip(profile.shear_rate_s[frame_idx]))
        title.set_text(f"管内流速与剪切率 — t = {profile.times_s[frame_idx] / 60.0:.1f} min")
        return ()

    animation = FuncAnimation(fig, update_frame, frames=len(profile.times_s),
                              interval=interval_ms, blit=False)
    target = Path(output_dir)
    target.mkdir(parents=True, exist_ok=True)
    path = target / f"{_safe_filename_component(well_name)}_管内段速度剪切率动画.gif"
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    try:
        animation.save(path, writer="pillow", fps=fps)
    finally:
        plt.close(fig)
    return path


__all__ = [
    "animate_casing_velocity",
    "animate_stitched_displacement",
    "plot_stitched_snapshots",
]
