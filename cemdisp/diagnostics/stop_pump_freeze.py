# -*- coding: utf-8 -*-
"""停泵衰减判据后处理（Phase 5a-①）：Q=0 后界面是否冻结

判据来源
--------
Moyers-González & Frigaard —— 停泵后 Q=0，界面冻结当且仅当屈服抗力压得住浮力驱动：
``τ_{Y,min}/(1+e) ≥ ‖div f‖∞ / 2``。与 R2（屈服门进流场）**正交**：R2 讲流动中的通道
导纳，本模块讲静止后的冻结判据；二者共用同一 τy 场，但作用时段不同。

实现口径（**显式映射假设，待追认**）
-----------------------------------
计划原文给的是无量纲判据的简写。本实现取**全尺寸（dimensional）形式**，避免跨
Z&F22 无量纲化的歧义，且与仓内既有 wall 通道的应力尺度同型
（``_yield_gate_wall`` 用 ``τw = G·b/2`` 外推，本模块用同一个「体力梯度 × 半间隙」尺度）：

- **驱动侧**（浮力体力散度 → 壁面应力尺度）：
  槽流中体力梯度 ``G = ∂f/∂x`` 作用在半间隙 ``b/2`` 上的等效壁面剪应力 = ``G·b/2``；
  即 ``‖div f‖∞/2 · (b/2)``，其中 ``f`` 为**方位向浮力体力** ``Δρ_az(z)·g·sinβ`` [Pa/m]，
  ``div`` 沿间隙坐标 ``x``。取 ``‖∂f/∂x‖∞ ≈ Δρ_az(z)·g·sinβ / b(z)``
  （方位向密度差在一个间隙宽度内完成变化）
  ⇒ ``driving = Δρ_az·g·sinβ·(b/2)/2`` [Pa]。
- **抗力侧**：``τ_{Y,min}(z)/(1+e(z))``，``e = 1 − standoff``（逐深取界面带内最大值）。
- **冻结**：``resisting ≥ driving``。

边界与已知局限
--------------
- ``Δρ_az`` 为**方位向**密度差（同一深度、宽窄边之间），非径向/纵向密度差——驱动 横向窜流的是前者。
- 本判据**只看最终态**（停泵时刻的场），不含时间演化：不预测冻结所需时间。
- ``‖∂f/∂x‖`` 用 ``Δρ_az/b`` 差分近似（一阶），不做间隙内多点数重构。
- 结果**不进 summary**（关1 红线：诊断字段一律不进 summary dict）。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

__all__ = ["StopPumpFreezeResult", "compute_stop_pump_freeze", "fields_from_state"]

_G = 9.81  # m/s²


@dataclass(frozen=True)
class StopPumpFreezeResult:
    """停泵冻结判据结果（逐深数组 + 裁决）。"""

    tau_y_min_pa: np.ndarray        # (nz,) 混合屈服应力场逐深最小值
    eccentricity: np.ndarray        # (nz,) e = 1 − standoff
    gap_m: np.ndarray               # (nz,) 局部间隙
    delta_rho_az_kg_m3: np.ndarray  # (nz,) 方位向密度差
    driving_pa: np.ndarray          # (nz,) 驱动侧应力尺度
    resisting_pa: np.ndarray        # (nz,) 抗力侧 τ_Y,min/(1+e)
    ratio: np.ndarray               # (nz,) resisting / driving
    frozen: np.ndarray              # (nz,) bool
    frozen_frac: float
    verdict: str


def _as_depth_vector(v, nz: int, name: str) -> np.ndarray:
    a = np.asarray(v, dtype=float)
    if a.ndim == 0:
        return np.full(nz, float(a))
    a = a.ravel()
    if a.size != nz:
        raise ValueError(f"{name} 长度须为 {nz} 或标量，实际 {a.size}")
    return a


def _azimuthal_spread(arr: np.ndarray) -> np.ndarray:
    """(ny,nz) → (nz,) 逐深方位向极差 max−min。"""
    a = np.asarray(arr, dtype=float)
    if a.ndim == 1:
        return np.zeros(a.shape, dtype=float)
    return a.max(axis=0) - a.min(axis=0)


def compute_stop_pump_freeze(
    *,
    tau_y_field: np.ndarray,
    rho_field: np.ndarray,
    gap_m: np.ndarray,
    standoff,
    beta_deg,
    g: float = _G,
) -> StopPumpFreezeResult:
    """按 Moyers-González & Frigaard 判据给出逐深冻结裁决。

    Args:
        tau_y_field: (ny,nz) 混合屈服应力场 [Pa]（与 `_compute_props` 的 τy 场同规则）。
        rho_field: (ny,nz) 密度场 [kg/m³]（逐相体积分数 × 相密度之和）。
        gap_m: (nz,) 局部间隙 b [m]。
        standoff: (nz,) 或标量，居中度 [0,1]。
        beta_deg: (nz,) 或标量，井斜角 [°]。
        g: 重力加速度 [m/s²]。
    """
    ty = np.asarray(tau_y_field, dtype=float)
    rho = np.asarray(rho_field, dtype=float)
    b = np.asarray(gap_m, dtype=float).ravel()
    if b.ndim != 1 or b.size == 0:
        raise ValueError("gap_m 须为非空一维数组")
    nz = b.size
    if ty.shape[-1] != nz or rho.shape[-1] != nz:
        raise ValueError(
            f"τy/ρ 场末维须等于 nz={nz}，实际 {ty.shape} / {rho.shape}"
        )
    if np.any(b <= 0.0):
        raise ValueError("间隙 b 须为正（m）")
    so = _as_depth_vector(standoff, nz, "standoff")
    beta = _as_depth_vector(beta_deg, nz, "beta_deg")

    tau_y_min = ty.min(axis=0) if ty.ndim > 1 else ty
    ecc = 1.0 - so
    delta_rho = _azimuthal_spread(rho)

    driving = delta_rho * g * np.abs(np.sin(np.deg2rad(beta))) * (b / 2.0) / 2.0
    resisting = tau_y_min / (1.0 + ecc)

    # 零驱动（水平段 / 方位向无密度差）⇒ 无浮力驱动，判为静止稳定（冻结）
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(driving > 0.0, resisting / driving, np.inf)
    frozen = ratio >= 1.0
    frac = float(np.mean(frozen))
    if frac >= 0.999:
        verdict = "完全冻结"
    elif frac <= 0.001:
        verdict = "未冻结"
    else:
        verdict = f"部分冻结({frac * 100:.1f}%)"

    return StopPumpFreezeResult(
        tau_y_min_pa=tau_y_min, eccentricity=ecc, gap_m=b,
        delta_rho_az_kg_m3=delta_rho, driving_pa=driving, resisting_pa=resisting,
        ratio=ratio, frozen=frozen, frozen_frac=frac, verdict=verdict,
    )


def fields_from_state(
    *,
    composition: dict,
    phase_rho_kg_m3: dict,
    phase_tau_y_pa: dict,
    gap_m: np.ndarray,
    standoff,
    beta_deg,
    g: float = _G,
    flusher: Optional[np.ndarray] = None,
) -> StopPumpFreezeResult:
    """由「末态相份额 + 逐相物性」重建 (ρ 场, τy 场) 再判据（便利包装）。

    密度/τy 取**体积分数加权**（与 `_compute_props` 的 τy 混合规则同型）；
    泥浆相 = ``1 − Σ已知相``（残差兜底，保证权重和为 1）。
    """
    if not composition:
        raise ValueError("composition 不能为空")
    ny, nz = np.asarray(next(iter(composition.values()))).shape
    rho = np.zeros((ny, nz), dtype=float)
    ty = np.zeros((ny, nz), dtype=float)
    total = np.zeros((ny, nz), dtype=float)
    for n, c in composition.items():
        cc = np.asarray(c, dtype=float)
        if cc.shape != (ny, nz):
            raise ValueError(f"相 {n} 形状 {cc.shape} 与 ({ny},{nz}) 不符")
        total += cc
        rho += cc * float(phase_rho_kg_m3[n])
        ty += cc * float(phase_tau_y_pa[n])
    if flusher is not None:
        total = total + np.asarray(flusher, dtype=float)
    residual = np.clip(1.0 - total, 0.0, 1.0)
    if "mud" in phase_rho_kg_m3:
        rho += residual * float(phase_rho_kg_m3["mud"])
        ty += residual * float(phase_tau_y_pa["mud"])
    return compute_stop_pump_freeze(
        tau_y_field=ty, rho_field=rho, gap_m=gap_m,
        standoff=standoff, beta_deg=beta_deg, g=g,
    )
