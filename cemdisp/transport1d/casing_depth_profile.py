# -*- coding: utf-8 -*-
"""管内段「深度 × 时间 × 相份额」重建（erf 混浆带口径）。

本模块是**纯后处理**：不改变 ``CasingFlowSolver`` 的任何数值行为，只消费
``CasingFlowResult`` 与同一套体积账原语，把管内一维解重建为随时间演化的深度
剖面，供可视化层使用。求解器冻结是该模块的硬约束（见 2026-09-26 设计规格）。

口径（2026-09-26 裁定）：
- 管容坐标 ``v(z)``：``pipe_id_profile`` 梯形积分，鞋口标定到 1D 求解器管容
  （与 ``scripts/entrypoints/export_depth_time_concentration.py`` 同口径）；
  ``pipe_id_profile`` 缺失时退化为均匀管径。
- 体积账：``u(z,t) = V(t) − v(z)``；``V(t)`` 为累计泵入体积，碰压后冻结在
  ``S_末段水泥 + v(鞋口)``（胶塞工艺：替浆不进环空，2026-09-03 裁定）。
- 混浆带：界面阈值 ``T_k`` 的到达时刻 ``τ_k(z)`` 由 V 反解；
  ``A_k(z,t) = 0.5·(1+erf((t−τ_k)/σ_k(z)))``。
  ``σ_k(z) = √(2·D_eff·t_contact)/U``，夹取沿用 ``_apply_dispersion_to_timeline``：
  ``min(σ, 0.5·t_travel)``、``max(σ, dt)``。
- 相份额由界面到达分数望远镜求和：``share_泥浆 = 1 − A_0``、
  ``share_段k = A_k − A_{k+1}``、``share_末段 = A_{M−1}``；各界面 σ 不同导致的
  负值先 clip 再归一化到 Σ = 1（2026-09-26 Q17 裁定）。
- 胶塞面 σ ≡ 0（工艺锐界面，Q16 裁定），不施加任何夹取。
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Sequence

import numpy as np
from scipy.special import erf

if TYPE_CHECKING:  # pragma: no cover - 仅类型检查
    from cemdisp.data.fluid_spec import FluidSpec
    from cemdisp.data.pumping_schedule import PumpingSchedule
    from cemdisp.data.well_spec import WellSpec
    from cemdisp.transport1d.casing_flow import CasingFlowResult, CasingFlowSolver

# erf⁻¹(0.99)：|A − 0.5| ≤ 0.49 对应的 |t − τ|/σ，用于把「带宽」定义成可测量
_BAND_Z_HALF = 1.1630871536766740
# 混浆带自洽下限契约（2026-09-26 Q17b + Q21A，纯数值口径，不含现场标定）
_BAND_MIN_CELLS = 3.0
_BAND_MIN_DOMAIN_FRAC = 1.0e-3
# 锐界面自动判据：相邻深度格某相份额跳变阈值
_SHARP_JUMP = 0.98
_DOMAIN_TOL_M = 1.0e-6


@dataclass(frozen=True)
class CasingDepthProfile:
    """管内段重建结果。

    所有数组按 ``(时间, 深度)`` 或 ``(界面, 深度)`` 排列，深度轴升序（浅→深）。
    """

    times_s: np.ndarray              # (n_t,)
    depths_m: np.ndarray             # (n_z,) 升序
    fluid_names: tuple[str, ...]     # (n_f,)，index 0 = 开泵前管内初始流体
    shares: np.ndarray               # (n_t, n_z, n_f)，Σ_f = 1
    sigma_t_s: np.ndarray            # (n_b, n_z) 界面过渡时间尺度（秒）
    band_width_m: np.ndarray         # (n_b, n_z) 解析混浆带宽度（米）
    boundary_arrival_s: np.ndarray   # (n_b, n_z)，不可达为 np.inf
    plug_boundary: tuple[bool, ...]  # (n_b,) 该界面是否胶塞面
    mixing_band: bool
    velocity_m_s: np.ndarray         # (n_t, n_z) 截面平均流速 Q(t)/A(z)
    shear_rate_s: np.ndarray         # (n_t, n_z) 壁面剪切率 8U/(2R)
    pipe_area_m2: np.ndarray         # (n_z,) 管内截面积（绘制井筒竖条宽度用）

    def channel_shares(self, fluid_name: str) -> np.ndarray:
        """取某一相的份额场 (n_t, n_z)。"""

        if fluid_name not in self.fluid_names:
            raise KeyError(
                f"未参与本次顶替的流体: {fluid_name!r}，可选: {list(self.fluid_names)}"
            )
        return self.shares[:, :, self.fluid_names.index(fluid_name)]

    def sharp_cells(self) -> np.ndarray:
        """锐界面标记 (n_t, n_z−1)：相邻深度格存在 ≥ _SHARP_JUMP 的相份额跳变。"""

        if self.shares.shape[1] < 2:
            return np.zeros((self.shares.shape[0], 0), dtype=bool)
        jump = np.max(np.abs(np.diff(self.shares, axis=1)), axis=2)
        return jump >= _SHARP_JUMP


# ── 体积账原语（与 CasingFlowSolver / 导出脚本同口径）──────────────────────


def _build_steps(schedule: "PumpingSchedule"):
    """构建胶塞截断后的注入步骤序列（与 CasingFlowSolver.run 同口径）。"""

    from cemdisp.transport1d.casing_flow import CasingFlowSolver

    full = CasingFlowSolver._build_scheduled_steps(schedule)
    steps, _ = CasingFlowSolver._displacement_sequence_cutoff(full)
    return steps


def _merged_segments(steps) -> list[tuple[str, float, float]]:
    """把相邻同名流体的注入步合并为段 → [(流体名, 起始体积, 结束体积), ...]。"""

    segments: list[tuple[str, float, float]] = []
    for scheduled in steps:
        name = scheduled.step.fluid_name
        if segments and segments[-1][0] == name:
            segments[-1] = (name, segments[-1][1], scheduled.cumulative_volume_end_m3)
        else:
            segments.append(
                (name, scheduled.cumulative_volume_start_m3, scheduled.cumulative_volume_end_m3)
            )
    return segments


def _boundary_thresholds(segments) -> np.ndarray:
    """界面阈值体积坐标 T_k：T_0 = 0（首个流体抵达），T_k = 段 k−1 的结束体积。

    末段的结束体积不构成界面——越过它仍是末段流体（与
    ``CasingFlowSolver._fluid_by_injected_volume`` 的「u ≥ S_末 → 末步流体」同口径）。
    """

    return np.array(
        [0.0 if k == 0 else segments[k - 1][2] for k in range(len(segments))], dtype=float
    )


def _pipe_volume_profile(
    well_spec: "WellSpec", depths_m: np.ndarray
) -> tuple[np.ndarray, np.ndarray, float]:
    """管容坐标 v(z) 与局部截面积 A(z)：pipe_id_profile 梯形积分 + 鞋口标定。

    返回 (vp_m3, area_m2, k_scale)；``pipe_id_profile`` 缺失或不可用时退化为均匀管径。
    """

    from cemdisp.transport1d.casing_flow import CasingFlowSolver

    shoe_m = float(well_spec.shoe_md_m)
    mean_area_m2 = float(CasingFlowSolver._pipe_cross_section_area(well_spec))
    target_m3 = shoe_m * mean_area_m2

    profile = tuple(getattr(well_spec, "pipe_id_profile", ()) or ())
    if len(profile) >= 2:
        p_depths = np.array([p.depth_md_m for p in profile], dtype=float)
        p_ids = np.array([p.value for p in profile], dtype=float)
        order = np.argsort(p_depths)
        p_depths, p_ids = p_depths[order], p_ids[order]
        keep = (p_depths >= -_DOMAIN_TOL_M) & (p_depths <= shoe_m + _DOMAIN_TOL_M)
        p_depths, p_ids = p_depths[keep], p_ids[keep]
        if len(p_depths) >= 2:
            if p_depths[0] > _DOMAIN_TOL_M:
                p_depths = np.concatenate(([0.0], p_depths))
                p_ids = np.concatenate(([p_ids[0]], p_ids))
            p_areas = math.pi * (p_ids / 1000.0) ** 2 / 4.0
            seg = (p_areas[:-1] + p_areas[1:]) * 0.5 * np.diff(p_depths)
            cum = np.concatenate(([0.0], np.cumsum(seg)))
            if cum[-1] > 0.0:
                k_scale = target_m3 / float(cum[-1])
                return (
                    np.interp(depths_m, p_depths, cum * k_scale),
                    np.interp(depths_m, p_depths, p_areas),
                    k_scale,
                )

    # 退化：均匀等效截面积（无 pipe_id_profile 的井走这里）
    return depths_m * mean_area_m2, np.full_like(depths_m, mean_area_m2), 1.0


def _invert_volume_to_time(steps, target_volume_m3: float) -> float | None:
    """体积坐标 → 时刻（与 _front_arrival_time 同款分段线性反解）。"""

    for scheduled in steps:
        if target_volume_m3 <= scheduled.cumulative_volume_end_m3 + 1.0e-12:
            if scheduled.step.rate_m3_min <= 0.0:
                return scheduled.end_time_s
            volume_into = max(target_volume_m3 - scheduled.cumulative_volume_start_m3, 0.0)
            return scheduled.start_time_s + volume_into / scheduled.step.rate_m3_min * 60.0
    return None


def _rate_at_time(steps, time_s: float) -> float:
    """某时刻的活动排量 [m³/s]；不在任何注入步内（停泵）时返回 0。"""

    for scheduled in steps:
        if scheduled.start_time_s <= time_s < scheduled.end_time_s - 1.0e-12:
            return scheduled.step.rate_m3_min / 60.0
    return 0.0


def _last_cement_step(steps, fluids: "tuple[FluidSpec, ...]"):
    """末段水泥浆步骤（LEAD/INTERMEDIATE/TAIL，名称兜底）→ 碰压冻结体积口径。"""

    from cemdisp.data.fluid_spec import FluidRole

    role_by_name = {fluid.name: fluid.role for fluid in fluids}
    cement_roles = {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}
    cement = [s for s in steps if role_by_name.get(s.step.fluid_name) in cement_roles]
    if not cement:
        cement = [
            s
            for s in steps
            if ("尾浆" in s.step.fluid_name or "水泥" in s.step.fluid_name
                or "领浆" in s.step.fluid_name)
        ]
    return cement[-1] if cement else None


# ── 主入口 ────────────────────────────────────────────────────────────


def build_casing_depth_profile(
    solver: "CasingFlowSolver",
    well_spec: "WellSpec",
    fluids: "tuple[FluidSpec, ...]",
    schedule: "PumpingSchedule",
    result: "CasingFlowResult",
    *,
    depths_m: Sequence[float],
    times_s: Sequence[float],
    mixing_band: bool = True,
) -> CasingDepthProfile:
    """把管内一维解重建为「深度 × 时间 × 相份额」。

    Args:
        solver: 已运行的 CasingFlowSolver（只读其弥散系数与体积账原语）。
        well_spec: 井筒规格（提供 shoe_md_m / pipe_id_profile / liner_id_mm）。
        fluids: 流体规格，用于弥散系数与胶塞面判据。
        schedule: 施工程序。
        result: 1D 求解结果，仅用于鞋深口径一致性校验。
        depths_m: 目标深度网格（米，升序，应落在 [0, shoe]）。
        times_s: 目标时刻网格（秒）。
        mixing_band: False 时退化为锐界面（σ ≡ 0），用于与 09-09 口径对照。

    Returns:
        CasingDepthProfile。

    Raises:
        ValueError: 施工程序为空、泵序流体名不在 fluids 中、或鞋深/深度域口径不一致。
    """

    depths = np.asarray(depths_m, dtype=float)
    times = np.asarray(times_s, dtype=float)
    if depths.ndim != 1 or times.ndim != 1:
        raise ValueError("depths_m 与 times_s 必须是一维序列")
    if depths.size < 2 or times.size < 1:
        raise ValueError("depths_m 至少 2 个点、times_s 至少 1 个点")
    if np.any(np.diff(depths) <= 0.0):
        raise ValueError("depths_m 必须严格升序")

    # 口径一致性校验：深度域以 well_spec 为准，result 必须与之同源
    if result is not None and hasattr(result, "shoe_md_m"):
        if abs(float(result.shoe_md_m) - float(well_spec.shoe_md_m)) > _DOMAIN_TOL_M:
            raise ValueError(
                f"result.shoe_md_m={result.shoe_md_m} 与 well_spec.shoe_md_m="
                f"{well_spec.shoe_md_m} 不一致，拒绝重建"
            )
    if float(depths[-1]) > float(well_spec.shoe_md_m) + _DOMAIN_TOL_M:
        raise ValueError(
            f"最深目标深度 {float(depths[-1])} m 超出鞋深 {well_spec.shoe_md_m} m"
        )

    steps = _build_steps(schedule)
    if not steps:
        raise ValueError("施工程序截断后为空，无法重建管内剖面")

    fluid_by_name = {fluid.name: fluid for fluid in fluids}
    segments = _merged_segments(steps)
    for name, _, _ in segments:
        if name not in fluid_by_name:
            raise ValueError(f"泵序流体 {name!r} 不在 fluids 中，无法计算弥散系数")

    initial_fluid = solver._initial_fluid_name(fluids, schedule)
    fluid_names = (initial_fluid,) + tuple(name for name, _, _ in segments)

    vp_m3, area_m2, _k_scale = _pipe_volume_profile(well_spec, depths)
    radius_m = np.sqrt(np.maximum(area_m2, 0.0) / math.pi)

    # 碰压冻结体积：末段水泥结束 + 鞋口管容（替浆不进环空，2026-09-03 裁定）
    last_cement = _last_cement_step(steps, fluids)
    v_shoe_m3 = float(vp_m3[-1])
    v_frozen_m3 = (
        last_cement.cumulative_volume_end_m3 + v_shoe_m3
        if last_cement is not None
        else float("inf")
    )

    thresholds = _boundary_thresholds(segments)
    n_boundary = len(segments)
    boundary_arrival = np.full((n_boundary, depths.size), np.inf, dtype=float)
    sigma_t = np.zeros((n_boundary, depths.size), dtype=float)
    band_width = np.zeros((n_boundary, depths.size), dtype=float)
    plug_boundary: list[bool] = []

    for k in range(n_boundary):
        arriving_name = segments[k][0]
        arriving_fluid = fluid_by_name[arriving_name]
        # 胶塞面：穿过该阈值的流体是胶塞释放液（压塞液）→ 工艺锐界面
        is_plug = bool(solver._is_plug_release_fluid(arriving_name, fluids))
        plug_boundary.append(is_plug)

        threshold_m3 = float(thresholds[k])
        t_inject = _invert_volume_to_time(steps, threshold_m3)
        t_inject = 0.0 if t_inject is None else t_inject

        for j in range(depths.size):
            target_volume = float(vp_m3[j]) + threshold_m3
            if target_volume > v_frozen_m3 + 1.0e-12:
                continue  # 冻结体积之外：该界面永不出现 → A ≡ 0
            tau = _invert_volume_to_time(steps, target_volume)
            if tau is None:
                continue
            boundary_arrival[k, j] = tau
            if is_plug or not mixing_band:
                continue  # 胶塞面 / 锐界面模式：σ ≡ 0，不施加夹取
            rate_m3_s = _rate_at_time(steps, tau)
            if rate_m3_s <= 0.0 or area_m2[j] <= 0.0:
                continue
            u_local = rate_m3_s / area_m2[j]
            d_eff = solver._compute_dispersion_coefficient(
                float(radius_m[j]), arriving_fluid, float(u_local)
            )
            t_travel = float(depths[j]) / u_local
            sigma = solver._contact_time_integrated_sigma(
                tau, t_inject, t_travel, d_eff, u_local, solver.dt
            )
            sigma = min(sigma, 0.5 * t_travel)
            sigma = max(sigma, solver.dt)
            sigma_t[k, j] = sigma
            # 带内判据 |A − 0.5| ≤ 0.49 ⇒ |t − τ| ≤ 1.163σ，换算到深度方向 ×U
            band_width[k, j] = 2.0 * _BAND_Z_HALF * sigma * u_local

    shares = _compose_shares(times, boundary_arrival, sigma_t, n_boundary)

    # 管内速度/剪切率：活塞流口径下全管截面同速 U(t,z) = Q(t)/A(z)（径向均匀，
    # 本模型不含速度剖面）；壁面剪切率沿用求解器的 8U/(2R) 约定。
    q_t = np.array([_rate_at_time(steps, float(t)) for t in times], dtype=float)
    safe_area = np.where(area_m2 > 0.0, area_m2, 1.0)
    velocity = q_t[:, np.newaxis] / safe_area[np.newaxis, :]
    safe_radius = np.where(radius_m > 0.0, radius_m, 1.0)
    shear_rate = 8.0 * velocity / (2.0 * safe_radius[np.newaxis, :])

    return CasingDepthProfile(
        times_s=times,
        depths_m=depths,
        fluid_names=fluid_names,
        shares=shares,
        sigma_t_s=sigma_t,
        band_width_m=band_width,
        boundary_arrival_s=boundary_arrival,
        plug_boundary=tuple(plug_boundary),
        mixing_band=bool(mixing_band),
        velocity_m_s=velocity,
        shear_rate_s=shear_rate,
        pipe_area_m2=area_m2,
    )


def _compose_shares(
    times_s: np.ndarray,
    boundary_arrival_s: np.ndarray,
    sigma_t_s: np.ndarray,
    n_boundary: int,
) -> np.ndarray:
    """由界面到达分数合成相份额 (n_t, n_z, n_f)，并归一化到 Σ = 1。

    望远镜求和：``share_初始 = 1 − A_0``、``share_段k = A_k − A_{k+1}``、
    ``share_末段 = A_{M−1}``。各界面 σ 不同时中间项可能为负，先 clip 再归一化。
    """

    n_t = times_s.size
    n_z = boundary_arrival_s.shape[1]

    arrival = boundary_arrival_s[:, np.newaxis, :]   # (n_b, 1, n_z)
    sigma = sigma_t_s[:, np.newaxis, :]              # (n_b, 1, n_z)
    t = times_s[np.newaxis, :, np.newaxis]           # (1, n_t, 1)

    reachable = np.isfinite(arrival)
    safe_sigma = np.where(sigma > 0.0, sigma, 1.0)
    frac_erf = 0.5 * (1.0 + erf((t - arrival) / safe_sigma))
    frac_step = np.where(t >= arrival, 1.0, 0.0)
    frac = np.where(sigma > 0.0, frac_erf, frac_step)
    frac = np.where(reachable, frac, 0.0)            # (n_b, n_t, n_z)

    raw = np.empty((n_boundary + 1, n_t, n_z), dtype=float)
    raw[0] = 1.0 - frac[0]
    for k in range(n_boundary - 1):
        raw[k + 1] = frac[k] - frac[k + 1]
    raw[n_boundary] = frac[n_boundary - 1]

    np.clip(raw, 0.0, None, out=raw)
    total = raw.sum(axis=0, keepdims=True)
    raw /= np.where(total > 0.0, total, 1.0)
    return np.transpose(raw, (1, 2, 0))


def band_contract_violations(profile: CasingDepthProfile) -> list[str]:
    """混浆带自洽下限契约（Q17b + Q21A）：带宽 ≥ max(3Δz, 0.1% 域长)。

    仅对**实际存在混浆带**（σ > 0）的界面逐深度检查；胶塞面为工艺锐界面，豁免。
    返回违规描述列表，空列表代表通过。
    """

    violations: list[str] = []
    depths = profile.depths_m
    if depths.size < 2:
        return ["深度网格不足 2 个点，无法评估带宽契约"]
    dz = float(np.median(np.diff(depths)))
    domain = float(depths[-1] - depths[0])
    floor = max(_BAND_MIN_CELLS * dz, _BAND_MIN_DOMAIN_FRAC * domain)

    for k, is_plug in enumerate(profile.plug_boundary):
        if is_plug:
            continue
        active = profile.sigma_t_s[k] > 0.0
        if not np.any(active):
            continue
        width_min = float(np.min(profile.band_width_m[k][active]))
        if width_min < floor:
            violations.append(
                f"界面 {k}（{profile.fluid_names[k]}→{profile.fluid_names[k + 1]}）"
                f"最小带宽 {width_min:.4g} m < 下限 {floor:.4g} m"
                f"（{_BAND_MIN_CELLS:g}Δz 与 {_BAND_MIN_DOMAIN_FRAC:.1%} 域长的较大者）"
            )
    return violations


__all__ = ["CasingDepthProfile", "build_casing_depth_profile", "band_contract_violations"]
