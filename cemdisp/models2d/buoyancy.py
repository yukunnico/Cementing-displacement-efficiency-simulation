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


_YIELD_TERM_MU_CLIP_PA_S = 3.0
"""屈服项表观黏度的**上限**（Pa·s）——**结构常数，禁作标定钮**。

与 ``AnnulusD2DGA._apparent_viscosity`` 的 ``np.clip(mu, 1.0e-5, 3.0)``（场口径）**同值同义**：
两条口径共用同一上限，避免出现第三套黏度口径。自仓内首个提交 ``0f78af9`` 起 3.0 即为场口径的
既有常数；本模块只是把它引到标量口径上。

依据（Phase 2 spec §1.2）：①口径唯一性（主）；②八井六速实测边界——最低转速（3 rpm，
γ̇ = 5.11 s⁻¹）实测最大 μ_app ≈ 2.60 Pa·s，3.0 与最后一个可信低转速实测点同量级；
③在模型工作域内不生效——触发条件 γ̇ < τy/(3.0−μp) ≈ 2.5~14.4 s⁻¹，落在六速最低转速之下或刚跨过；
④文献：屈服应力流体 τy/γ̇ 在 γ̇→0 的奇异性必须有界正则化（Saramito & Wachs 2017；Papanastasiou 1987）。
与 ``f_safety=1.15`` 同型纪律：只做区间披露，**不参与任何反标定**。
"""


def _clip_yield_term(mu: float) -> float:
    """``include_yield_term=True`` 分支的上限夹逼——与 ``annulus._apparent_viscosity`` 逐字对齐。"""
    return min(max(float(mu), 1.0e-5), _YIELD_TERM_MU_CLIP_PA_S)


def fluid_apparent_viscosity(fluid: FluidSpec, shear_rate: float, *,
                             include_yield_term: bool = False) -> float:
    """表观黏度。幂律/HB 用 ``K·γ̇^(n−1)``；回退到 PV。缺参数时抛错，不静默回退。

    Args:
        fluid: 流体规格（`FluidSpec`）。幂律/HB 必须带 ``consistency_k`` 与
            ``power_law_n``；牛顿/Bingham 用 ``plastic_viscosity_pa_s``。
        shear_rate: 剪切率 γ̇，单位 1/s（调用方约定与 ``_compute_props`` 一致：
            ``γ̇ = 6|w|/b``）。非正值被夹到 1e-8 避免幂律奇异。
        include_yield_term: **R2 μp/τy 真拆分开关**（Phase 2，默认 ``False`` = 旧语义逐位不变）。

            * ``False``（默认）——HEAD 行为：Bingham/牛顿只返 ``PV``；幂律/HB 只返
              ``K·γ̇^(n−1)``（带 ``max(·, PV)`` 地板）。**本分支的两条 return 与 HEAD 位级相同**（关 2 红线）。
            * ``True``——补上屈服应力贡献，使**标量口径与场口径同构**
              （``annulus._apparent_viscosity:1121/:1130`` 的同一本构）：

              * Bingham / 牛顿：``μ = PV + τy/γ̇``
              * HB：``μ = τy/γ̇ + K·γ̇^(n−1)``（τy 项在前，与 ``:1130`` 字面对齐）
              * 幂律：无 τy ⇒ 与 ``False`` 分支同式

              并**无条件**施加上限 ``_clip_yield_term``（与场口径 ``:1133`` 同构）。

            **本构依据**：Bingham (1916) 塑性本构 / API RP 13D 钻井液表观黏度定义
            ``μ_app = τ/γ̇ = μp + τy/γ̇``——现场六速表即以 AV/PV/YP 报告。

            **已知与场口径的差异（显式声明，不掩盖）**：场口径 γ 地板 = 1e-6（本处 1e-8），
            且 HB 分支**无** ``max(μ, PV)`` 地板 ⇒ 两者**逐位同值仅在**「γ̇ ≥ 1e-6 且
            （Bingham 且 τy ≥ (3−PV)·1e-6，或 HB 且 ``plastic_viscosity_pa_s is None``）」时成立。

    Returns:
        表观黏度，单位 Pa·s（Z&F22 ``μ̂₁``）。
    """
    g = max(float(shear_rate), 1e-8)
    if fluid.rheology_model in (RheologyModel.POWER_LAW, RheologyModel.HERSCHEL_BULKLEY):
        if fluid.consistency_k is not None and fluid.power_law_n is not None:
            mu = float(fluid.consistency_k) * g ** (float(fluid.power_law_n) - 1.0)
            if include_yield_term and fluid.yield_stress_pa:
                # τy 项在前（与 _apparent_viscosity:1130 书写次序一致；IEEE-754 下加法可交换）
                mu = float(fluid.yield_stress_pa) / g + mu
            if fluid.plastic_viscosity_pa_s:
                mu = max(mu, float(fluid.plastic_viscosity_pa_s))
            return _clip_yield_term(mu) if include_yield_term else mu
    if fluid.plastic_viscosity_pa_s is not None:
        mu = float(fluid.plastic_viscosity_pa_s)
        if include_yield_term:
            # 无条件 clip（含 τy=0 的 Bingham 与牛顿流体）：与场口径 :1133 同构
            if fluid.yield_stress_pa:
                mu += float(fluid.yield_stress_pa) / g
            return _clip_yield_term(mu)
        return mu
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
    """Froude 数平方（Z&F22 (2.6)）。Task 4 起替代 `_buoyancy_force_vector` 内硬编码的 F2 = 1.0。

    ``F = √(τ̂₀/(ρ̂₁·ĝ·δ₀·r̂ₐ*))`` ⇒ ``F² = τ̂₀/(ρ̂₁·ĝ·δ₀·r̂ₐ*)``，其中
    ``τ̂₀ = μ̂₁·ŵ₀/d̂`` 为被顶替液中的黏性应力尺度（注意 F² 是比值 τ̂₀/(浮力尺度)，
    不是其倒数）。

    **δ₀ 与 r̂ₐ* 的取值约定**（论文出处：§2.1 与 (2.6)）：

    - ``r̂ₐ*`` = 沿环空流道平均的代表性半径（论文："The mean radius ``r̂ₐ*`` is defined
      by averaging along the annular flow path"），量纲 m，即本函数的 ``mean_radius_m``。
    - ``δ₀`` = **无量纲**参考间隙比。论文 (2.1) 的窄间隙参数**本身就是**
      ``δ = d̂/r̂ₐ*``——原文写作 ``d̂/(πr̂ₐ*) = δ/π ≪ 1``（π 只出现在 ``δ/π`` 的写法里，
      不在 ``δ`` 自身上）；(2.6) 的 ``δ₀`` 取其参考值 ⇒ ``δ₀ = d̂/r̂ₐ*``（无量纲），
      等价于 ``δ₀·r̂ₐ* = d̂``。
      该取值使 (2.5b)/(2.6) 的 ``|b| ≈ (ρ−1)/F²``（论文 p.8："b 即浮力向量 b 的大小"）
      与论文 p.8 的浮力数 ``b = Δρ·ĝ·d̂²/(μ̂₁ŵ₀)`` 精确对齐，联立给出
      ``F²·b = Δρ/ρ̂₁``（Atwood 数），等价于 ``F² = μ̂₁ŵ₀/(ρ̂₁·ĝ·d̂²)``。
      ⇒ 调用方应传 ``gap_scale_m = half_gap_m / mean_radius_m``。
      ``gap_scale_m`` 与 ``half_gap_m`` **不是同一个量**（前者无量纲、后者长度），
      但二者乘积恒等于 ``d̂``，故在本式分母中只以乘积 ``δ₀·r̂ₐ* = d̂`` 起作用。
      ⚠️ **不要再给 δ₀ 乘或除 π**：论文的 δ 就是 ``d̂/r̂ₐ*``，加 π 会让 F² 整体差 π 倍。

    Args:
        mu_displaced: 被顶替液（钻井液）表观黏度 μ̂₁，Pa·s。
        w0_mps: 截面平均轴向速度 ŵ₀ = q/A，m/s。
        half_gap_m: 半间隙 d̂，m（= (r_o−r_i)/2 = (井径−外径)/4）。
        rho_displaced: 被顶替液密度 ρ̂₁，kg/m³。
        gap_scale_m: 无量纲参考间隙比 δ₀ = d̂/r̂ₐ*（= half_gap_m/mean_radius_m）。
            ⚠️ 形参名保留 ``_m`` 后缀是 Task 3 的历史遗留，**语义已改为无量纲**。
        mean_radius_m: 沿程平均环空半径 r̂ₐ*，m（= mean((井径+外径)/4)）。

    Returns:
        F²（无量纲）。呼101 实测 **O(10⁻³)（2.0×10⁻³ ~ 5.5×10⁻³）**；
        八井整体跨度 **[1.2×10⁻³, 3.1×10⁻²]**（最大为呼102）。
    """
    d = max(float(half_gap_m), 1e-9)
    tau0 = float(mu_displaced) * max(float(w0_mps), 1e-9) / d
    denom = max(float(rho_displaced) * G * max(float(gap_scale_m), 1e-9)
                * max(float(mean_radius_m), 1e-9), 1e-12)
    return tau0 / denom
