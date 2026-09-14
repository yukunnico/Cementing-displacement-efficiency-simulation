"""浮力口径统一模块（Zhang & Frigaard 2022）。

全仓唯一的浮力定义入口，消除 `annulus_d2dga.py` 内并存的两套顶替液口径。
文献锚点：

- 浮力数 ``b = (ρ̂₂ − ρ̂₁)·ĝ·d̂²/(μ̂₁·ŵ₀)`` —— Z&F22 p.8
  （``ρ̂₁``/``μ̂₁`` = **被顶替液**（钻井液），``ρ̂₂`` = **顶替液**（水泥浆）；
  ``d̂`` = 半间隙 = ``(r_o − r_i)/2 = (井径 − 外径)/4``；``ŵ₀`` = 截面平均轴向速度。
  下标约定与 ``τ̂₀ = μ̂₁ŵ₀/d̂`` 一致：μ̂₁ 与被顶替液配对。）
- Froude 数 ``F = √(τ̂₀/(ρ̂₁·ĝ·δ₀·r̂ₐ*))``，即
  ``F² = τ̂₀/(ρ̂₁·ĝ·δ₀·r̂ₐ*)`` —— Z&F22 (2.6)，其中 ``τ̂₀ = μ̂₁ŵ₀/d̂``

单位口径（本模块内一律 SI，调用方负责换算）：

- 密度 kg/m³；半间隙 m；黏度 Pa·s；速度 m/s；剪切率 1/s；半径 m。
- ``b`` 与 ``F²`` 均为无量纲数。

修复的三处旧口径缺陷（Task 3）：

1. **双口径**：`_compute_velocity` 用 ``0.67×领浆 + 0.33×尾浆``，而 summary 段的
   浮力数只用领浆 → ``b`` 符号可能翻转。本模块 ``displacing_density_kg_m3`` 为
   **全仓唯一**口径。
2. **静默回退**：幂律/HB 泥浆被 ``plastic_viscosity_pa_s or 0.05`` 回退到硬编码
   0.05 Pa·s。本模块 ``fluid_apparent_viscosity`` 对幂律/HB 一律用 ``K·γ̇^(n−1)``，
   缺参数直接抛错。
3. **末步速度**：``w₀`` 原取“最后一步”速度场均值。调用方改用截面平均速度
   ``q/A``（末态泵注排量 / 环形截面积）。
"""

from __future__ import annotations

from cemdisp.data.fluid_spec import FluidSpec, RheologyModel

G = 9.81
LEAD_WEIGHT = 0.67   # 体积加权口径，与 annulus_d2dga._compute_velocity 一致


def displacing_density_kg_m3(lead_fluid, tail_fluid, mud_fluid) -> float:
    """顶替液代表密度（0.67×领浆 + 0.33×尾浆）。全仓唯一口径。

    与 `annulus_d2dga._compute_velocity` 的 ``rho_disp`` 体积加权口径一致；
    领浆/尾浆缺失时逐级退化，二者皆无则取泥浆密度。
    返回单位 kg/m³（Z&F22 ``ρ̂₂`` = 顶替液；被顶替液为 ``ρ̂₁``）。
    """
    if lead_fluid is not None and tail_fluid is not None:
        return LEAD_WEIGHT * lead_fluid.density_kg_m3 + (1 - LEAD_WEIGHT) * tail_fluid.density_kg_m3
    if lead_fluid is not None:
        return float(lead_fluid.density_kg_m3)
    if tail_fluid is not None:
        return float(tail_fluid.density_kg_m3)
    return float(mud_fluid.density_kg_m3)


def fluid_apparent_viscosity(fluid: FluidSpec, shear_rate: float) -> float:
    """表观黏度。幂律/HB 用 ``K·γ̇^(n−1)``；回退到 PV。缺参数时抛错，不静默回退。

    Args:
        fluid: 流体规格（`FluidSpec`）。幂律/HB 必须带 ``consistency_k`` 与
            ``power_law_n``；牛顿/Bingham 用 ``plastic_viscosity_pa_s``。
        shear_rate: 剪切率 γ̇，单位 1/s（调用方约定与 ``_compute_props`` 一致：
            ``γ̇ = 6|w|/b``）。非正值被夹到 1e-8 避免幂律奇异。

    Returns:
        表观黏度，单位 Pa·s（Z&F22 ``μ̂₁``）。
    """
    g = max(float(shear_rate), 1e-8)
    if fluid.rheology_model in (RheologyModel.POWER_LAW, RheologyModel.HERSCHEL_BULKLEY):
        if fluid.consistency_k is not None and fluid.power_law_n is not None:
            mu = float(fluid.consistency_k) * g ** (float(fluid.power_law_n) - 1.0)
            if fluid.plastic_viscosity_pa_s:
                mu = max(mu, float(fluid.plastic_viscosity_pa_s))
            return mu
    if fluid.plastic_viscosity_pa_s is not None:
        return float(fluid.plastic_viscosity_pa_s)
    raise ValueError(f"流体 {fluid.name} 缺少可用流变参数，禁止静默回退")


def buoyancy_number(rho_displacing: float, rho_displaced: float, half_gap_m: float,
                    mu_displaced: float, w0_mps: float) -> float:
    """无量纲浮力数 b（Z&F22 p.8）。b>0 密度稳定；b<0 密度倒置（文献警告严格避免）。

    ``b = (ρ_displacing − ρ_displaced)·g·d²/(μ_displaced·w₀)``，``d`` 为半间隙。

    Args:
        rho_displacing: 顶替液密度 ρ̂₂，kg/m³。
        rho_displaced: 被顶替液（泥浆）密度 ρ̂₁，kg/m³。
        half_gap_m: 半间隙 d̂，m（= 全间隙/2 = (井径−外径)/4）。
        mu_displaced: 被顶替液表观黏度 μ̂₁，Pa·s。
        w0_mps: 截面平均轴向速度 ŵ₀，m/s。
    """
    d = max(float(half_gap_m), 1e-9)
    denom = max(float(mu_displaced) * max(float(w0_mps), 1e-9), 1e-12)
    return (float(rho_displacing) - float(rho_displaced)) * G * d * d / denom


def froude_squared(mu_displaced: float, w0_mps: float, half_gap_m: float,
                   rho_displaced: float, gap_scale_m: float, mean_radius_m: float) -> float:
    """Froude 数平方（Z&F22 (2.6)）。替代此前硬编码的 F2 = 1.0。

    ``F = √(τ̂₀/(ρ̂₁·ĝ·δ₀·r̂ₐ*))`` ⇒ ``F² = τ̂₀/(ρ̂₁·ĝ·δ₀·r̂ₐ*)``，其中
    ``τ̂₀ = μ̂₁·ŵ₀/d̂`` 为黏性应力尺度（注意 F² 是比值 τ̂₀/(浮力尺度)，不是其倒数）。

    Args:
        mu_displaced: 被顶替液表观黏度 μ̂₁，Pa·s。
        w0_mps: 截面平均轴向速度 ŵ₀，m/s。
        half_gap_m: 半间隙 d̂，m（= 全间隙/2 = (井径−外径)/4）。
        rho_displaced: 被顶替液密度 ρ̂₁，kg/m³。
        gap_scale_m: 环空间隙尺度 δ₀，m（全间隙 = 2d̂）。
        mean_radius_m: 平均环空半径 r̂ₐ*，m（= (井径+外径)/4）。
    """
    d = max(float(half_gap_m), 1e-9)
    tau0 = float(mu_displaced) * max(float(w0_mps), 1e-9) / d
    denom = max(float(rho_displaced) * G * max(float(gap_scale_m), 1e-9)
                * max(float(mean_radius_m), 1e-9), 1e-12)
    return tau0 / denom
