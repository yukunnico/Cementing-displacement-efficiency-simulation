# -*- coding: utf-8 -*-
"""管内段「深度 × 时间 × 相份额」重建契约测试（2026-09-26 设计规格 §3/§4）。

合成井夹具，不依赖现场资料：

- 锐界面模式与 ``scripts/entrypoints/export_depth_time_concentration.py`` 的
  体积账重建**逐位一致**（这是 Q9「冻结求解器」下唯一可用的强回归锚）；
- 相份额望远镜求和 Σ = 1（含 erf 带模式）；
- 胶塞面 σ ≡ 0（工艺锐界面，Q16），其余界面进入 erf 混浆带（Q11/Q17）；
- 带宽随接触时间增长（深度方向单调）；
- 混浆带自洽下限契约（Q17b + Q21A）在粗网格上如实报违规、在细网格上通过；
- 深度域 / 鞋深口径守卫；
- 重建不改变求解器与结果对象（冻结证明的最小形式）。
"""

from __future__ import annotations

import importlib.util
import math
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
from cemdisp.transport1d.casing_flow import CasingFlowSolver
from cemdisp.transport1d import casing_depth_profile
from cemdisp.transport1d.casing_depth_profile import (
    CHANNELS,
    band_contract_violations,
    band_width_summary,
    profile_channels,
    unify_phase_channel,
    build_casing_depth_profile,
)

_REPO_ROOT = Path(__file__).resolve().parents[2]
_SCRIPT_PATH = _REPO_ROOT / "scripts" / "entrypoints" / "export_depth_time_concentration.py"

SHOE_M = 100.0
AREA_M2 = 0.01           # 管内截面积 → 管容 1.0 m³
PIPE_ID_MM = math.sqrt(4.0 * AREA_M2 / math.pi) * 1000.0


@pytest.fixture(scope="module")
def export_mod():
    """按文件路径加载导出脚本模块（scripts/ 非 package）。"""

    spec = importlib.util.spec_from_file_location("_export_dtc_test", _SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _fluids() -> tuple[FluidSpec, ...]:
    """井浆（初始）+ 平衡液 / 尾浆 / 压塞液（胶塞释放液）/ 替钻井液。"""

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
    """鞋深 100 m 的合成井，管容 1.0 m³（无 pipe_id_profile → 均匀管径分支）。"""

    return WellSpec(
        well_name="合成井",
        top_md_m=1.0,
        bottom_md_m=120.0,
        shoe_md_m=SHOE_M,
        liner_id_mm=PIPE_ID_MM,
    )


def _schedule() -> PumpingSchedule:
    """平衡液 0.5 → 尾浆 0.5 → 压塞液 0.2 → 替钻井液 1.0（累计 2.2 m³）。"""

    def step(name: str, fluid: str, volume: float, rate: float,
             t0: float, t1: float, tag: PumpingStageEvent) -> PumpingScheduleStep:
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
def solved():
    """跑一次合成井 1D 求解，全模块复用（只读）。"""

    well, fluids, schedule = _well(), _fluids(), _schedule()
    solver = CasingFlowSolver()
    result = solver.run(well, fluids, schedule)
    return solver, well, fluids, schedule, result


def _depths(n: int) -> np.ndarray:
    return np.linspace(0.0, SHOE_M, n)


def _times(n: int, horizon_s: float = 300.0) -> np.ndarray:
    return np.linspace(0.0, horizon_s, n)


# ── 1. 锐界面模式与导出脚本逐位一致（冻结下的强回归锚）──────────────────────


def test_sharp_mode_matches_export_volume_ledger(solved, export_mod):
    """mixing_band=False 时，主导相必须与导出脚本的体积账重建逐点同名。"""

    solver, well, fluids, schedule, result = solved
    depths, times = _depths(41), _times(31)

    profile = build_casing_depth_profile(
        solver, well, fluids, schedule, result,
        depths_m=depths, times_s=times, mixing_band=False,
    )

    # 复现导出脚本的重建口径：V(t) 冻结在 S_尾浆 + v(鞋口)
    steps_full = CasingFlowSolver._build_scheduled_steps(schedule)
    steps, _ = CasingFlowSolver._displacement_sequence_cutoff(steps_full)
    last_cement = next(s for s in reversed(steps) if s.step.fluid_name == "尾浆")
    vp_depths, vp_vals, _k = export_mod.build_pipe_volume_profile(
        [0.0, SHOE_M], [PIPE_ID_MM, PIPE_ID_MM], SHOE_M, SHOE_M * AREA_M2,
    )
    v_frozen = last_cement.cumulative_volume_end_m3 + float(vp_vals[-1])
    initial = CasingFlowSolver._initial_fluid_name(fluids, schedule)

    # 界面阈值体积坐标 T_k：用于判定不一致点是否恰好落在界面重合点上
    thresholds = casing_depth_profile._boundary_thresholds(
        casing_depth_profile._merged_segments(steps)
    )

    mismatches: list[tuple] = []
    for t_idx, t_s in enumerate(times):
        v_now = min(CasingFlowSolver._cumulative_volume_at(steps, float(t_s)), v_frozen)
        for z_idx, z_m in enumerate(depths):
            vp_z = float(np.interp(z_m, vp_depths, vp_vals))
            u = v_now - vp_z
            expected = export_mod.fluid_name_at_volume_coordinate(u, steps, initial)
            got = profile.fluid_names[int(np.argmax(profile.shares[t_idx, z_idx]))]
            if got != expected:
                tie_distance = float(np.min(np.abs(thresholds - u)))
                mismatches.append((float(t_s), float(z_m), got, expected, tie_distance))

    # 允许的差异只有一类：u 恰好等于某界面阈值 T_k（测度零的重合点，标签取决于
    # 半开区间约定与浮点最后一位）。任何"远离界面却不同名"都是真实口径错误。
    off_interface = [m for m in mismatches if m[4] > 1.0e-9]
    assert not off_interface, (
        f"锐界面模式与导出脚本口径不一致（且不落在界面重合点），共 {len(off_interface)} 处："
        f"{off_interface[:5]}"
    )


def test_sharp_mode_is_binary(solved):
    """锐界面模式份额只能取 0 或 1（带宽恒为 0）。"""

    solver, well, fluids, schedule, result = solved
    profile = build_casing_depth_profile(
        solver, well, fluids, schedule, result,
        depths_m=_depths(41), times_s=_times(31), mixing_band=False,
    )
    interior = (profile.shares > 1.0e-12) & (profile.shares < 1.0 - 1.0e-12)
    assert not interior.any(), "锐界面模式出现了中间份额"


# ── 2. 份额闭合 ───────────────────────────────────────────────────────


@pytest.mark.parametrize("mixing_band", [True, False])
def test_shares_sum_to_one(solved, mixing_band):
    """每个（时刻, 深度）上各相份额之和恒为 1。"""

    solver, well, fluids, schedule, result = solved
    profile = build_casing_depth_profile(
        solver, well, fluids, schedule, result,
        depths_m=_depths(61), times_s=_times(41), mixing_band=mixing_band,
    )
    total = profile.shares.sum(axis=2)
    assert np.allclose(total, 1.0, atol=1.0e-12), (
        f"份额和偏离 1：max |Σ−1| = {np.abs(total - 1).max():.3e}"
    )
    assert np.all(profile.shares >= 0.0), "出现负份额"


# ── 3. 混浆带 ────────────────────────────────────────────────────────


def test_mixing_band_creates_interior_shares(solved):
    """erf 带模式在界面附近产生 0 < share < 1 的格点（Q11 的核心诉求）。"""

    solver, well, fluids, schedule, result = solved
    profile = build_casing_depth_profile(
        solver, well, fluids, schedule, result,
        depths_m=_depths(201), times_s=_times(61), mixing_band=True,
    )
    interior = (profile.shares > 0.01) & (profile.shares < 0.99)
    assert int(interior.sum()) >= 20, f"混浆带格点数只有 {int(interior.sum())}，带未形成"


def test_plug_face_stays_sharp(solved):
    """胶塞面（压塞液前缘）σ ≡ 0，是混合带口径的显式例外（Q16）。"""

    solver, well, fluids, schedule, result = solved
    profile = build_casing_depth_profile(
        solver, well, fluids, schedule, result,
        depths_m=_depths(101), times_s=_times(41), mixing_band=True,
    )

    plug_indices = [k for k, is_plug in enumerate(profile.plug_boundary) if is_plug]
    assert plug_indices == [2], f"胶塞面定位错误：{plug_indices}（期望唯一的压塞液前缘）"
    k = plug_indices[0]
    assert np.all(profile.sigma_t_s[k] == 0.0), "胶塞面 σ 不为零"
    assert np.all(profile.band_width_m[k] == 0.0), "胶塞面产生了带宽"

    for other in range(len(profile.plug_boundary)):
        if other != k:
            assert np.any(profile.sigma_t_s[other] > 0.0), f"界面 {other} 未生成混浆带"


def test_band_width_grows_with_contact_time():
    """物理量级真实的长井上，带宽随深度（接触时间）增长，且**不**由数值下限支配。

    合成小井（100 m / 高排量）上接触时间只有数十秒，σ 被 ``max(σ, dt)`` 这个
    数值下限夹住 → 带宽由 dt 决定而非弥散物理。本用例改用长井 + 低排量，使
    σ ≫ dt，从而检验真正的物理行为。
    """

    shoe_m, area_m2 = 3000.0, 0.01          # 管容 30 m³
    rate = 0.06                              # m³/min → U = 0.1 m/s
    pipe_id_mm = math.sqrt(4.0 * area_m2 / math.pi) * 1000.0
    well = WellSpec(well_name="长井", top_md_m=1.0, bottom_md_m=3050.0,
                    shoe_md_m=shoe_m, liner_id_mm=pipe_id_mm)

    def step(name, fluid, volume, t0, tag):
        return PumpingScheduleStep(
            step_name=name, fluid_name=fluid, volume_m3=volume, rate_m3_min=rate,
            start_time_s=t0, end_time_s=t0 + volume / rate * 60.0, event_tag=tag,
        )

    schedule = PumpingSchedule(steps=(
        step("注入平衡液", "平衡液", 5.0, 0.0, PumpingStageEvent.INJECT_SPACER),
        step("注入尾浆", "尾浆", 15.0, 5000.0, PumpingStageEvent.INJECT_CEMENT),
        step("注入压塞液", "压塞液", 1.0, 20000.0, PumpingStageEvent.INJECT_CEMENT),
        step("顶替", "替钻井液", 20.0, 21000.0, PumpingStageEvent.INJECT_DISPLACEMENT),
    ))
    fluids = _fluids()
    solver = CasingFlowSolver()
    result = solver.run(well, fluids, schedule)

    depths = np.linspace(0.0, shoe_m, 3001)
    times = np.linspace(0.0, 41000.0, 21)
    profile = build_casing_depth_profile(
        solver, well, fluids, schedule, result,
        depths_m=depths, times_s=times, mixing_band=True,
    )

    k = 0  # 平衡液前缘
    sigma = profile.sigma_t_s[k]
    widths = profile.band_width_m[k]

    # 1) 深部 σ 必须显著大于数值下限 dt —— 否则"带宽"测的是离散化而不是弥散
    assert float(np.max(sigma)) > 10.0 * solver.dt, (
        f"σ 最大值 {float(np.max(sigma)):.3f} s 未显著超过 dt={solver.dt} s，"
        "带宽被数值下限支配，无法检验弥散物理"
    )

    # 2) 未被下限夹住的深部，带宽随深度单调增长
    deep = sigma > solver.dt * 1.001
    assert int(deep.sum()) >= 100, "未被下限夹住的深度点太少"
    assert np.all(np.diff(widths[deep]) >= -1.0e-12), (
        "深部带宽未随深度单调增长（接触时间越长带应越宽）"
    )
    # 3) 最深处带宽应显著大于最浅处的有效带宽
    assert float(widths[deep][-1]) > 2.0 * float(widths[deep][0]), (
        "带宽沿深度的增长幅度不足"
    )


# ── 4. 带宽自洽下限契约 ────────────────────────────────────────────────


def test_band_contract_flags_coarse_grid(solved):
    """粗深度网格下带宽低于 3Δz，契约必须如实报违规（不静默放行）。"""

    solver, well, fluids, schedule, result = solved
    profile = build_casing_depth_profile(
        solver, well, fluids, schedule, result,
        depths_m=_depths(11), times_s=_times(21), mixing_band=True,
    )
    violations = band_contract_violations(profile)
    assert violations, "11 点深度网格下应报带宽违规（3Δz=30 m 远超实际带宽）"
    assert all("最小带宽" in v for v in violations)


def test_band_contract_passes_on_fine_grid(solved):
    """深度网格足够密时契约通过。"""

    solver, well, fluids, schedule, result = solved
    profile = build_casing_depth_profile(
        solver, well, fluids, schedule, result,
        depths_m=_depths(2001), times_s=_times(21), mixing_band=True,
    )
    assert band_contract_violations(profile) == []


def test_band_contract_uses_explicit_cell_size(solved):
    """(C)：放大窗尺度以显式分辨率判定——整井网格报违规，放大窗分辨率通过。"""

    solver, well, fluids, schedule, result = solved
    profile = build_casing_depth_profile(
        solver, well, fluids, schedule, result,
        depths_m=_depths(61), times_s=_times(21), mixing_band=True,
    )
    assert band_contract_violations(profile), "61 点整井网格（Δz=1.67 m）应报违规"
    assert band_contract_violations(profile, cell_size_m=0.1) == [], (
        "放大窗分辨率 0.1 m 下带宽 ≥ 3 格，应通过"
    )


def test_band_width_summary_reports_visibility(solved):
    """带宽诊断摘要必须如实给出中位带宽与占域长比例（不美化）。"""

    solver, well, fluids, schedule, result = solved
    profile = build_casing_depth_profile(
        solver, well, fluids, schedule, result,
        depths_m=_depths(1001), times_s=_times(21), mixing_band=True,
    )
    summary = band_width_summary(profile)
    assert int(summary["n_interfaces"]) >= 1
    assert 0.0 < float(summary["median_band_m"]) < float(summary["domain_m"])
    assert 0.0 < float(summary["median_domain_frac"]) < 1.0
    assert float(summary["min_band_m"]) <= float(summary["median_band_m"])


def test_band_contract_exempts_plug_face(solved):
    """胶塞面为工艺锐界面，豁免带宽契约。"""

    solver, well, fluids, schedule, result = solved
    profile = build_casing_depth_profile(
        solver, well, fluids, schedule, result,
        depths_m=_depths(2001), times_s=_times(21), mixing_band=True,
    )
    k = profile.plug_boundary.index(True)
    assert np.all(profile.sigma_t_s[k] == 0.0)
    assert not any(f"界面 {k}（" in v for v in band_contract_violations(profile))


# ── 5. 口径守卫 ──────────────────────────────────────────────────────


def test_rejects_depth_beyond_shoe(solved):
    """目标深度超出鞋深必须报错，不外推。"""

    solver, well, fluids, schedule, result = solved
    with pytest.raises(ValueError, match="超出鞋深"):
        build_casing_depth_profile(
            solver, well, fluids, schedule, result,
            depths_m=np.array([0.0, SHOE_M + 5.0]), times_s=_times(5),
        )


def test_rejects_non_monotonic_depths(solved):
    solver, well, fluids, schedule, result = solved
    with pytest.raises(ValueError, match="严格升序"):
        build_casing_depth_profile(
            solver, well, fluids, schedule, result,
            depths_m=np.array([0.0, 50.0, 40.0]), times_s=_times(5),
        )


def test_rejects_shoe_mismatch(solved):
    """result 与 well_spec 鞋深不一致时拒绝重建（防止张冠李戴）。"""

    solver, well, fluids, schedule, result = solved
    other = WellSpec(
        well_name="另一口井", top_md_m=1.0, bottom_md_m=220.0,
        shoe_md_m=200.0, liner_id_mm=PIPE_ID_MM,
    )
    with pytest.raises(ValueError, match="不一致"):
        build_casing_depth_profile(
            solver, other, fluids, schedule, result,
            depths_m=np.array([0.0, 100.0, 190.0]), times_s=_times(5),
        )


# ── 6. 冻结：重建不得改动求解器与结果 ──────────────────────────────────


def test_rebuild_does_not_mutate_solver_or_result(solved):
    """Q9「冻结」的最小形式：重建前后 result 全字段逐位不变。"""

    solver, well, fluids, schedule, result = solved

    def snapshot(res):
        return {
            "fronts": [(f.fluid_name, float(f.distance_m), float(f.time_s))
                       for f in res.fronts],
            "pipe_cross_section_m2": float(res.pipe_cross_section_m2),
            "shoe_md_m": float(res.shoe_md_m),
            "pumping_end_time_s": float(res.pumping_end_time_s),
            "cement_end_time_s": float(res.cement_end_time_s),
            "timeline_len": len(res.shoe_timeline.events),
        }

    before = snapshot(result)
    build_casing_depth_profile(
        solver, well, fluids, schedule, result,
        depths_m=_depths(301), times_s=_times(51), mixing_band=True,
    )
    assert snapshot(result) == before, "重建改动了 1D 求解结果（违反冻结约束）"
    assert float(solver.dispersion_alpha) == 0.25, "重建改动了 dispersion_alpha"


# ── 7. 相名统一（交付物 6 的必做前置）────────────────────────────────


# 呼1-004 实际出现的 11 个泵序流体名（取自真实 loader）
_REAL_FLUID_NAMES = ("钻井液", "先导浆", "隔离液1", "隔离液2", "领浆", "尾浆",
                     "压塞液", "替钻井液", "保护液", "基液", "井浆")


def test_unify_phase_channel_matches_export_script(export_mod):
    """相名统一必须与导出脚本 map_fluid_to_channel 逐名一致（两侧同口径）。"""

    for name in _REAL_FLUID_NAMES:
        assert unify_phase_channel(name) == export_mod.map_fluid_to_channel(name), (
            f"{name!r} 映射不一致：包内 {unify_phase_channel(name)} vs "
            f"导出脚本 {export_mod.map_fluid_to_channel(name)}"
        )


def test_unify_phase_channel_covers_four_channels(export_mod):
    """真井 11 个相名必须全部落到四通道，且四通道都被用到。"""

    mapped = {unify_phase_channel(name) for name in _REAL_FLUID_NAMES}
    assert mapped == set(CHANNELS), f"四通道覆盖不全：{mapped}"


def test_unify_phase_channel_rejects_unknown_name():
    """未知名必须报错，**不静默归并**（静默归并会伪造浓度）。"""

    with pytest.raises(ValueError, match="不静默归并"):
        unify_phase_channel("某种没见过的流体")


def test_profile_channels_close_to_one(solved):
    """四通道聚合份额必须逐格闭合到 1，与相份额总和一致。"""

    solver, well, fluids, schedule, result = solved
    profile = build_casing_depth_profile(
        solver, well, fluids, schedule, result,
        depths_m=_depths(101), times_s=_times(21), mixing_band=True,
    )
    channels = profile_channels(profile)
    assert set(channels) == set(CHANNELS)
    total = sum(channels.values())
    assert np.allclose(total, 1.0, atol=1.0e-12), (
        f"四通道未闭合，max |Σ−1| = {np.abs(total - 1).max():.3e}"
    )


def test_unify_phase_channel_handles_hu101_zhongzhiye():
    """呼101 专有相名「中置液」必须归 mud（hu101_loader.py:408 标注其现场名为保护液）。"""

    assert unify_phase_channel("中置液") == "mud"


def test_unify_phase_channel_still_rejects_truly_unknown():
    """补名不得把「不静默归并」的红线一起放宽。"""

    with pytest.raises(ValueError, match="不静默归并"):
        unify_phase_channel("某种真的没见过的流体")
