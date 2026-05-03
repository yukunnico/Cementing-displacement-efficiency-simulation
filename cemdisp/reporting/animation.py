"""固井顶替过程的水泥浓度场动画输出。

本模块提供面向用户报告的动画生成功能，用于展示环空二维模型中
水泥浓度场随施工时间推进的演化过程。所有图表标签、标题和输出
文件名均使用中文，并复用 reporting.plots 中的中文字体与文件名清理逻辑。
"""

from __future__ import annotations

# pyright: reportAny=false, reportPrivateUsage=false, reportUnknownMemberType=false, reportUnusedCallResult=false

import logging
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from cemdisp.models2d import AnnulusSimulationResult
from cemdisp.reporting.plots import _safe_filename_component, _setup_chinese_font

LOGGER = logging.getLogger(__name__)


def animate_cement_field(
    result: AnnulusSimulationResult,
    output_dir: Path | str | None = None,
    interval_ms: int = 200,
    fps: int = 10,
    save_format: str = "gif",
) -> None:
    """生成水泥浓度场随时间演化的 GIF 或 MP4 动画。

    动画采用单面板热力图展示每个快照时刻的水泥浓度场：横轴为井深，
    纵轴为归一化方位角（0 表示宽边，1 表示窄边），颜色表示水泥浓度。

    Args:
        result: 环空二维模拟结果，需包含 cement_snapshots、snapshot_times_s 和 geom["md"]。
        output_dir: 动画输出目录；为 None 时保存到当前工作目录。
        interval_ms: 动画播放时相邻帧之间的间隔，单位 ms。
        fps: 保存文件时使用的帧率。
        save_format: 保存格式，仅支持 "gif" 或 "mp4"。

    Raises:
        KeyError: 当 result.geom 缺少 "md" 井深数组时抛出。
        ValueError: 当 save_format 不是 "gif" 或 "mp4" 时抛出。
    """
    _setup_chinese_font()

    snapshots = result.cement_snapshots
    if not snapshots:
        message = f"{result.well_name} 未包含水泥浓度场快照，跳过顶替过程动画生成。"
        LOGGER.warning(message)
        warnings.warn(message, RuntimeWarning, stacklevel=2)
        return

    normalized_format = save_format.lower().lstrip(".")
    if normalized_format not in {"gif", "mp4"}:
        raise ValueError("save_format 仅支持 'gif' 或 'mp4'")

    if "md" not in result.geom:
        raise KeyError("result.geom 缺少动画绘制所需的井深数组 'md'")

    # 井深数组来自通用几何输入，不在动画模块中硬编码任何单井深度。
    depth = result.geom["md"]
    depth_min = float(depth[0])
    depth_max = float(depth[-1])
    frame_count = len(snapshots)
    times_s = result.snapshot_times_s

    fig, ax = plt.subplots(figsize=(10, 6))
    image = ax.imshow(
        snapshots[0],
        vmin=0,
        vmax=1,
        cmap="viridis",
        aspect="auto",
        extent=(depth_min, depth_max, 0.0, 1.0),
        origin="lower",
    )
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("水泥浓度")

    ax.set_xlabel("井深 / m")
    ax.set_ylabel("方位角 (宽边→窄边)")
    ax.set_title("水泥浓度场 — t = 0.0 min", fontsize=14, fontweight="bold")

    def update_frame(frame_index: int):
        """刷新单帧数据与标题，并输出生成进度。"""
        cement_field = snapshots[frame_index]
        time_s = times_s[frame_index] if frame_index < len(times_s) else 0.0
        time_min = float(time_s) / 60.0

        image.set_data(cement_field)
        ax.set_title(f"水泥浓度场 — t = {time_min:.1f} min", fontsize=14, fontweight="bold")
        print(f"正在生成动画帧 {frame_index + 1}/{frame_count}")
        return (image,)

    animation = FuncAnimation(
        fig,
        update_frame,
        frames=frame_count,
        interval=interval_ms,
        blit=False,
    )

    output_path = Path("." if output_dir is None else output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    well_name = _safe_filename_component(result.well_name)
    file_path = output_path / f"{well_name}_顶替过程动画.{normalized_format}"

    writer = "pillow" if normalized_format == "gif" else "ffmpeg"
    fig.tight_layout()
    try:
        animation.save(file_path, writer=writer, fps=fps)
    finally:
        plt.close(fig)


__all__ = ["animate_cement_field"]
