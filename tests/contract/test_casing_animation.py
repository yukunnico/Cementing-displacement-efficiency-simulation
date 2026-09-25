# -*- coding: utf-8 -*-
"""管内—环空拼接可视化 smoke 测试（2026-09-26 设计规格 §3 交付物 1/2/3/4）。

不跑 2D 求解器：用**鸭子类型桩**提供 ``cement_snapshots`` / ``spacer_snapshots`` /
``snapshot_times_s`` / ``geom["md"]`` / ``well_name``（与
``AnnulusSimulationResult`` 同接口），合成井跑真 1D 求解得到管内剖面。

断言的是**交付物的机器可验判据**（规格 §3）与**帧对齐守卫**，不校验图像像素。
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
    _cement_channel,
    _front_lines,
    _spacer_channel,
    animate_casing_velocity,
    animate_stitched_displacement,
    plot_stitched_snapshots,
)
from cemdisp.transport1d.casing_flow import CasingFlowSolver
from cemdisp.transport1d.casing_depth_profile import build_casing_depth_profile

SHOE_M = 100.0
HANGER_M = 40.0          # 环空 2D 评价域上界（悬挂器）
AREA_M2 = 0.01
PIPE_ID_MM = math.sqrt(4.0 * AREA_M2 / math.pi) * 1000.0
N_TIMES = 6


@dataclass
class _AnnulusStub:
    """AnnulusSimulationResult 的最小鸭子类型桩（仅本模块读取的字段）。"""

    well_name: str
    cement_snapshots: tuple
    spacer_snapshots: tuple
    snapshot_times_s: tuple
    geom: dict = field(default_factory=dict)


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
    """(profile, 环空桩, 输出目录)：时间栅格两侧严格同源（Q22-b-i）。"""

    well, fluids, schedule = _well(), _fluids(), _schedule()
    solver = CasingFlowSolver()
    result = solver.run(well, fluids, schedule)

    times = np.linspace(0.0, 300.0, N_TIMES)
    depths = np.linspace(0.0, SHOE_M, 301)
    profile = build_casing_depth_profile(
        solver, well, fluids, schedule, result,
        depths_m=depths, times_s=times, mixing_band=True,
    )

    ny, nz = 6, 12
    md = np.linspace(SHOE_M, HANGER_M, nz)          # 与 2D 一致：降序（深→浅）
    rng = np.random.default_rng(20260926)
    cement = tuple(rng.uniform(0.0, 1.0, size=(ny, nz)) for _ in range(N_TIMES))
    spacer = tuple(rng.uniform(0.0, 0.5, size=(ny, nz)) for _ in range(N_TIMES))
    annulus = _AnnulusStub(
        well_name="合成井",
        cement_snapshots=cement,
        spacer_snapshots=spacer,
        snapshot_times_s=tuple(float(t) for t in times),
        geom={"md": md},
    )
    return profile, annulus, tmp_path_factory.mktemp("casing_anim")


def _is_gif(path: Path) -> bool:
    """GIF89a/GIF87a 魔数校验，避免「文件存在但是空的」假通过。"""

    return path.read_bytes()[:6] in (b"GIF89a", b"GIF87a")


def test_stitched_animation_is_a_real_gif(scene):
    """交付物 1：拼接 GIF 落盘且是真 GIF。"""

    profile, annulus, out = scene
    path = animate_stitched_displacement(profile, annulus, out, well_name="合成井", fps=5)
    assert path.exists() and path.stat().st_size > 1024
    assert _is_gif(path), "产出的不是 GIF（魔数不符）"
    assert "拼接" in path.name


def test_stitched_animation_rejects_frame_mismatch(scene):
    """帧不对齐必须报错，不允许两侧时间栅格各说各话（Q22-b-i）。"""

    profile, annulus, out = scene
    short = _AnnulusStub(
        well_name="合成井",
        cement_snapshots=annulus.cement_snapshots[:-1],
        spacer_snapshots=annulus.spacer_snapshots[:-1],
        snapshot_times_s=annulus.snapshot_times_s[:-1],
        geom=annulus.geom,
    )
    with pytest.raises(ValueError, match="不一致"):
        animate_stitched_displacement(profile, short, out, well_name="合成井")


def test_snapshot_png_written(scene):
    """交付物 3：论文用多帧拼贴 PNG。"""

    profile, annulus, out = scene
    path = plot_stitched_snapshots(profile, annulus, out, well_name="合成井", n_panels=4)
    assert path.exists() and path.stat().st_size > 1024
    assert path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n", "产出不是 PNG"


def test_velocity_animation_is_a_real_gif(scene):
    """交付物 4：管内段速度/剪切率二次渲染 GIF。"""

    profile, _annulus, out = scene
    path = animate_casing_velocity(profile, out, well_name="合成井", fps=5)
    assert path.exists() and _is_gif(path)
    assert "速度" in path.name


def test_velocity_and_shear_are_physical(scene):
    """速度/剪切率场必须是正量，且剪切率 = 8U/(2R) 自洽。"""

    profile, _annulus, _out = scene
    assert np.all(profile.velocity_m_s >= 0.0)
    assert np.all(profile.shear_rate_s >= 0.0)
    area = np.asarray(profile.pipe_area_m2, dtype=float)
    radius = np.sqrt(area / np.pi)
    expected = 8.0 * profile.velocity_m_s / (2.0 * radius[np.newaxis, :])
    assert np.allclose(profile.shear_rate_s, expected, rtol=1.0e-12)


def test_channel_picks_cement_and_spacer(scene):
    """代表通道识别：水泥取领浆/尾浆，隔离液取平衡液。"""

    profile, _annulus, _out = scene
    assert _cement_channel(profile) == "尾浆"
    assert _spacer_channel(profile) == "平衡液"


def test_plug_face_is_marked_among_front_lines(scene):
    """前缘轨迹里必须恰有一条被标为胶塞面（Q16-C 的标记来源）。"""

    profile, _annulus, _out = scene
    lines = _front_lines(profile, N_TIMES - 1)
    assert lines, "末帧没有任何界面，前缘轨迹为空"
    plug = [line for line in lines if line[2]]
    assert len(plug) == 1, f"胶塞面标记数应为 1，实得 {len(plug)}"
