# -*- coding: utf-8 -*-
"""管内—环空拼接可视化契约测试（2026-09-26 修订版）。

不跑 2D 求解器：用**鸭子类型桩**提供 ``lead/tail/spacer/cement_snapshots`` /
``snapshot_times_s`` / ``geom["md"]`` / ``well_name``（与 ``AnnulusSimulationResult``
同接口），合成井跑真 1D 求解得到整井与放大窗两份管内剖面。

覆盖 2026-09-26 用户裁定：
- **遮挡回归**：色标轴与数据轴在图形坐标下不得相交（这是被报的 bug）；
- 管内=全部流体、按相名分配可区分颜色、混合区做凸组合着色；
- 环空=领浆/尾浆/隔离液/钻井液 四通道且闭合；
- 放大窗存在且带宽契约按放大窗分辨率判定；
- 三侧帧数对齐守卫。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pytest

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.pumping_schedule import (
    PumpingSchedule,
    PumpingScheduleStep,
    PumpingStageEvent,
)
from cemdisp.data.well_spec import WellSpec
from cemdisp.reporting.casing_animation import (
    _build_figure,
    _front_lines,
    animate_stitched_displacement,
    annulus_channels,
    blend_rgb,
    phase_colors,
    plot_stitched_snapshots,
)
from cemdisp.transport1d.casing_flow import CasingFlowSolver
from cemdisp.transport1d.casing_depth_profile import (
    band_contract_violations,
    build_casing_depth_profile,
)

SHOE_M = 100.0
HANGER_M = 40.0          # 环空 2D 评价域上界（悬挂器）
AREA_M2 = 0.01
PIPE_ID_MM = math.sqrt(4.0 * AREA_M2 / math.pi) * 1000.0
N_TIMES = 6
ZOOM_CENTER_M = 70.0
ZOOM_HALF_M = 15.0


@dataclass
class _AnnulusStub:
    """AnnulusSimulationResult 的最小鸭子类型桩（仅本模块读取的字段）。"""

    well_name: str
    cement_snapshots: tuple
    spacer_snapshots: tuple
    snapshot_times_s: tuple
    geom: dict = field(default_factory=dict)
    lead_snapshots: tuple = ()
    tail_snapshots: tuple = ()


def _fluids() -> tuple[FluidSpec, ...]:
    return (
        FluidSpec("井浆", FluidRole.MUD, 1000.0, RheologyModel.NEWTONIAN,
                  plastic_viscosity_pa_s=0.01),
        FluidSpec("平衡液", FluidRole.SPACER, 1200.0, RheologyModel.NEWTONIAN,
                  plastic_viscosity_pa_s=0.02),
        FluidSpec("尾浆", FluidRole.TAIL, 1900.0, RheologyModel.NEWTONIAN,
                  plastic_viscosity_pa_s=0.05),
        FluidSpec("压塞液", FluidRole.OTHER, 1000.0, RheologyModel.NEWTONIAN,
                  plastic_viscosity_pa_s=0.01),
        FluidSpec("替钻井液", FluidRole.DISPLACEMENT, 1100.0, RheologyModel.NEWTONIAN,
                  plastic_viscosity_pa_s=0.01),
    )


def _well() -> WellSpec:
    return WellSpec(well_name="合成井", top_md_m=1.0, bottom_md_m=120.0,
                    shoe_md_m=SHOE_M, liner_id_mm=PIPE_ID_MM)


def _schedule() -> PumpingSchedule:
    def step(name, fluid, volume, rate, t0, t1, tag):
        return PumpingScheduleStep(
            step_name=name, fluid_name=fluid, volume_m3=volume, rate_m3_min=rate,
            start_time_s=t0, end_time_s=t1, event_tag=tag,
        )

    return PumpingSchedule(steps=(
        step("注入平衡液", "平衡液", 0.5, 0.5, 0.0, 60.0, PumpingStageEvent.INJECT_SPACER),
        step("注入尾浆", "尾浆", 0.5, 0.5, 60.0, 120.0, PumpingStageEvent.INJECT_CEMENT),
        step("注入压塞液", "压塞液", 0.2, 0.2, 120.0, 180.0, PumpingStageEvent.INJECT_CEMENT),
        step("顶替", "替钻井液", 1.0, 0.5, 180.0, 300.0, PumpingStageEvent.INJECT_DISPLACEMENT),
    ))


@pytest.fixture(scope="module")
def scene(tmp_path_factory):
    """(整井剖面, 放大窗剖面, 环空桩, 输出目录)：三侧时刻严格同源（Q22-b-i）。"""

    well, fluids, schedule = _well(), _fluids(), _schedule()
    solver = CasingFlowSolver()
    result = solver.run(well, fluids, schedule)
    times = np.linspace(0.0, 300.0, N_TIMES)

    profile = build_casing_depth_profile(
        solver, well, fluids, schedule, result,
        depths_m=np.linspace(0.0, SHOE_M, 301), times_s=times, mixing_band=True,
    )
    zoom_profile = build_casing_depth_profile(
        solver, well, fluids, schedule, result,
        depths_m=np.linspace(ZOOM_CENTER_M - ZOOM_HALF_M, ZOOM_CENTER_M + ZOOM_HALF_M, 401),
        times_s=times, mixing_band=True,
    )

    ny, nz = 6, 12
    md = np.linspace(SHOE_M, HANGER_M, nz)          # 与 2D 一致：降序（深→浅）
    rng = np.random.default_rng(20260926)
    lead = tuple(rng.uniform(0.0, 0.4, size=(ny, nz)) for _ in range(N_TIMES))
    tail = tuple(rng.uniform(0.0, 0.4, size=(ny, nz)) for _ in range(N_TIMES))
    spacer = tuple(rng.uniform(0.0, 0.2, size=(ny, nz)) for _ in range(N_TIMES))
    cement = tuple(np.clip(lead[i] + tail[i], 0.0, 1.0) for i in range(N_TIMES))
    annulus = _AnnulusStub(
        well_name="合成井",
        cement_snapshots=cement,
        spacer_snapshots=spacer,
        snapshot_times_s=tuple(float(t) for t in times),
        geom={"md": md},
        lead_snapshots=lead,
        tail_snapshots=tail,
    )
    return profile, zoom_profile, annulus, tmp_path_factory.mktemp("casing_anim")


def _is_gif(path: Path) -> bool:
    """GIF89a/GIF87a 魔数校验，避免「文件存在但是空的」假通过。"""

    return path.read_bytes()[:6] in (b"GIF89a", b"GIF87a")


def _overlaps(a, b) -> bool:
    """两个 matplotlib 轴的图形坐标包围盒是否相交。"""

    return not (a.x1 <= b.x0 or b.x1 <= a.x0 or a.y1 <= b.y0 or b.y1 <= a.y0)


# ── 1. 遮挡回归（用户报的 bug）────────────────────────────────────────


def test_colorbar_axes_do_not_overlap_data_axes(scene):
    """色标轴必须与所有数据轴分离——「右侧色标挡住图表」的几何回归。"""

    import matplotlib.pyplot as plt

    profile, zoom_profile, annulus, _out = scene
    parts = _build_figure(profile, zoom_profile, annulus)
    fig = parts["fig"]
    try:
        data_axes = [parts["casing_ax"], *parts["ann_axes"], parts["zoom_ax"]]
        cbar_axes = [ax for ax in fig.axes
                     if ax not in data_axes and len(ax.get_images()) == 0]
        assert len(cbar_axes) >= 4, f"独立色标轴不足：{len(cbar_axes)}"
        for cax in cbar_axes:
            for dax in data_axes:
                assert not _overlaps(cax.get_position(), dax.get_position()), (
                    f"色标轴 {cax.get_position().bounds} 与数据轴 "
                    f"{dax.get_title()!r} {dax.get_position().bounds} 重叠"
                )
    finally:
        plt.close(fig)


# ── 2. 相配色与混合 ──────────────────────────────────────────────────


def test_phase_colors_distinguish_all_phases(scene):
    """管内每一相都要有颜色，且同名关键词组内的多相必须可区分。"""

    profile, _zoom, _annulus, _out = scene
    colors = phase_colors(profile.fluid_names)
    assert set(colors) == set(profile.fluid_names), "有相未分配到颜色"
    dup_names = ("隔离液1", "隔离液2")
    colors_dup = phase_colors(dup_names)
    assert colors_dup[dup_names[0]] != colors_dup[dup_names[1]], (
        "同组的两个隔离液分到同一颜色，无法区分"
    )


def test_blend_rgb_is_convex_combination(scene):
    """纯相位置的混合色必须等于该相颜色；混合区落在色域内。"""

    profile, _zoom, _annulus, _out = scene
    names = profile.fluid_names
    colors = phase_colors(names)
    shares = np.zeros((3, len(names)))
    shares[0, names.index(names[0])] = 1.0
    shares[1, names.index(names[-1])] = 1.0
    shares[2, :] = 1.0 / len(names)
    rgb = blend_rgb(shares, names, colors)
    assert np.allclose(rgb[0], colors[names[0]], atol=1.0e-12)
    assert np.allclose(rgb[1], colors[names[-1]], atol=1.0e-12)
    assert np.all(rgb[2] >= 0.0) and np.all(rgb[2] <= 1.0)


# ── 3. 环空四通道 ────────────────────────────────────────────────────


def test_annulus_channels_are_four_and_closed(scene):
    """环空必须是四通道（领浆/尾浆/隔离液/钻井液）且逐格闭合到 1。"""

    _profile, _zoom, annulus, _out = scene
    channels, md = annulus_channels(annulus)
    assert set(channels) == {"领浆", "尾浆", "隔离液", "钻井液"}
    total = sum(channels.values())
    assert np.allclose(total, 1.0, atol=1.0e-12), "四通道未闭合"
    assert md.size == channels["尾浆"].shape[2]


def test_annulus_channels_fall_back_to_cement(scene):
    """无 lead/tail 时退化为「尾浆=cement」，不报错也不静默丢相。"""

    _profile, _zoom, annulus, _out = scene
    legacy = _AnnulusStub(
        well_name="合成井", cement_snapshots=annulus.cement_snapshots,
        spacer_snapshots=annulus.spacer_snapshots,
        snapshot_times_s=annulus.snapshot_times_s, geom=annulus.geom,
    )
    channels, _md = annulus_channels(legacy)
    assert np.all(channels["领浆"] == 0.0)
    assert np.allclose(channels["尾浆"], np.asarray(legacy.cement_snapshots))


# ── 4. 产物与守卫 ────────────────────────────────────────────────────


def test_stitched_animation_is_a_real_gif(scene):
    profile, zoom_profile, annulus, out = scene
    path = animate_stitched_displacement(profile, zoom_profile, annulus, out,
                                         well_name="合成井", fps=5)
    assert path.exists() and path.stat().st_size > 1024
    assert _is_gif(path)
    assert "拼接" in path.name


def test_stitched_animation_rejects_frame_mismatch(scene):
    """环空帧数不齐必须报错（Q22-b-i）。"""

    profile, zoom_profile, annulus, out = scene
    short = _AnnulusStub(
        well_name="合成井",
        cement_snapshots=annulus.cement_snapshots[:-1],
        spacer_snapshots=annulus.spacer_snapshots[:-1],
        snapshot_times_s=annulus.snapshot_times_s[:-1],
        geom=annulus.geom,
        lead_snapshots=annulus.lead_snapshots[:-1],
        tail_snapshots=annulus.tail_snapshots[:-1],
    )
    with pytest.raises(ValueError, match="不一致"):
        animate_stitched_displacement(profile, zoom_profile, short, out, well_name="合成井")


def test_stitched_animation_rejects_zoom_frame_mismatch(scene):
    """放大窗与整井帧数不齐也必须报错，不允许「放大窗停在前一帧」。"""

    profile, _zoom_profile, annulus, out = scene
    solver = CasingFlowSolver()
    well, fluids, schedule = _well(), _fluids(), _schedule()
    result = solver.run(well, fluids, schedule)
    bad_zoom = build_casing_depth_profile(
        solver, well, fluids, schedule, result,
        depths_m=np.linspace(ZOOM_CENTER_M - ZOOM_HALF_M, ZOOM_CENTER_M + ZOOM_HALF_M, 401),
        times_s=np.linspace(0.0, 300.0, N_TIMES - 1), mixing_band=True,
    )
    with pytest.raises(ValueError, match="放大窗"):
        animate_stitched_displacement(profile, bad_zoom, annulus, out, well_name="合成井")


def test_snapshot_png_written(scene):
    profile, zoom_profile, annulus, out = scene
    path = plot_stitched_snapshots(profile, zoom_profile, annulus, out,
                                   well_name="合成井", n_panels=4)
    assert path.exists() and path.stat().st_size > 1024
    assert path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n", "产出不是 PNG"


def test_zoom_window_resolves_the_band(scene):
    """(C)：放大窗深度分辨率下带宽必须 ≥ 3 格（整井粗网格则报违规）。

    注意合成小井的带宽被 ``max(σ, dt)`` 抬到约 3.9 m，因此"整井不够密"必须用
    真正的粗网格体现（11 点，Δz=10 m），不能用 301 点——否则断言的是错的东西。
    """

    _profile, zoom_profile, _annulus, _out = scene
    solver = CasingFlowSolver()
    well, fluids, schedule = _well(), _fluids(), _schedule()
    result = solver.run(well, fluids, schedule)
    coarse = build_casing_depth_profile(
        solver, well, fluids, schedule, result,
        depths_m=np.linspace(0.0, SHOE_M, 11), times_s=np.linspace(0.0, 300.0, N_TIMES),
        mixing_band=True,
    )
    assert band_contract_violations(coarse), "整井 11 点网格（Δz=10 m）应报带宽不足"
    zoom_cell = float(np.median(np.diff(zoom_profile.depths_m)))
    assert band_contract_violations(zoom_profile, cell_size_m=zoom_cell) == [], (
        "放大窗分辨率下带宽仍不足 3 格，放大窗失去意义"
    )


def test_plug_face_is_marked_among_front_lines(scene):
    profile, _zoom, _annulus, _out = scene
    lines = _front_lines(profile, N_TIMES - 1)
    assert lines, "末帧没有任何界面，前缘轨迹为空"
    assert len([line for line in lines if line[2]]) == 1, "胶塞面标记数应为 1"


def test_front_line_depth_is_deepest_reached(scene):
    """前缘深度必须取「已到达的最深处」，且随时间单调加深。

    回归 2026-09-26 实测缺陷：``np.argmax(crossed)`` 取到第 0 号格（井口），
    导致胶塞面标记停在深度 0、放大窗被定位到井口。
    """

    profile, _zoom, _annulus, _out = scene

    def cement_front(time_idx: int) -> float:
        depths = [depth for depth, _label, _plug in _front_lines(profile, time_idx)]
        assert depths, f"第 {time_idx} 帧没有前缘"
        return max(depths)

    fronts = [cement_front(idx) for idx in range(N_TIMES)]
    assert fronts[-1] > 0.0, f"末帧前缘仍停在井口（深度 {fronts[-1]}），取到的是第 0 号格"
    assert all(b >= a - 1.0e-12 for a, b in zip(fronts, fronts[1:])), (
        f"前缘深度未随时间单调加深：{fronts}"
    )
    assert fronts[-1] > fronts[0], f"末帧前缘未推进：首 {fronts[0]} 末 {fronts[-1]}"
