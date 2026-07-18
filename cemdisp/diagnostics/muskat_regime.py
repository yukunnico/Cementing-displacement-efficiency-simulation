"""
Muskat 三 regime 稳定性诊断 (muskat_regime)
=========================================

本模块实现 Bararpour & Frigaard (2025, JFM 1009, A15) §3.3 的 Muskat
指进稳定性分析（式 3.8-3.13），将环空顶替的弥散前缘分类为三种 regime：

1. **stable**（稳定）：指进不能穿透前缘到达以平均泵速弥散的浓度位置；
2. **partial_penetration**（部分穿透）：指进穿透至少 50% 的前缘，
   但在某临界浓度 c_critical 处停止，无法穿到纯被顶替液层；
3. **unstable**（不稳定）：指进在全部浓度区间都快于前缘，完全穿透。

物理方法
--------
比较假想窄指进（充满顶替液/水泥浆）速度 w_finger 与弥散前缘的运动学
波速 w_f（论文式号）：

- 式 3.8：w_f(c̄) = q0'(c̄) + b·I3'(c̄)（可微区间）；
- 式 3.9：激波区间用 Rankine-Hugoniot 条件 q(c̄M)−q(c̄m) = w_s(c̄M−c̄m)。
  本模块以通量函数 F(c̄)=q0+b·I3 的**上凹包络**实现：包络与 F 重合段为
  弥散（稀疏波）区（w_f=F'），弦段为激波区（w_f=弦斜率=w_s，即等面积规则）；
- 式 3.12：w_finger(c̄) = I1(1)/I1(c̄) + b·I1(1)·[I2(c̄)/I1(c̄) + c̄ − 1]；
- 式 3.13：Δw(c̄) = w_finger(c̄) − w_f(c̄)。

判据（论文 §3.3）：

- Δw(c̄) > 0 对全部 c̄∈[0,1] → unstable（不存在 c_critical）；
- 存在 c_critical（Δw=0 的最大浓度），Δw<0 于 [0, c_critical] 且
  c* < c_critical → stable（c* 为 w_f=1 的浓度，即以平均泵速弥散的浓度，
  是体积加权意义下前缘的中点）；
- 其余情形 → partial_penetration。

牛顿解析闭包（Tier 0 近似）
--------------------------
I1/I2/q0 采用两牛顿流体解析式（Bararpour 2025 式 2.24-2.26，几何平均
黏度标度 η1=m^0.5、η2=m^−0.5；Zhang 2022 式 4.21a,b 同构）：

- I1(c̄,m) = [√m·c̄³ + (1−c̄³)/√m]/3                （式 2.24，归一化掉 H³）
- I2(c̄,m) = [2√m·c̄³(1−c̄) + c̄(1−c̄)²(1+2c̄)/√m]/6  （式 2.25，归一化掉 H⁴）
- q0(c̄,m) = c̄·[m·c̄²+1.5(1−c̄²)]/[m·c̄³+1−c̄³]      （式 2.26，复用
  ``d2dga_flux_amplification`` 并解除 [0.5,2] 人工限幅）

I3 复用 ``d2dga_flux.d2dga_dispersion_function_I3``（Zhang & Frigaard 2022
式 4.26：分母 2m·[m·c̄³+1−c̄³]）。注意 Bararpour 2025 式 2.27 的分母为
12·[m·c̄³+1−c̄³]、分子含 3(1−c̄²)，两者存在标度差异；按 Tier 0 spec 约定
采用求解器现有 I3，差异留待 Tier 1 核实。

**m 约定**：黏度比 m = μ_displaced/μ_displacing = μ_mud/μ_cement，与
``d2dga_flux.py`` 及求解器一致（m>1 表示泥浆比水泥浆更黏，指进风险高；
牛顿等密度极限 m > m_c ≈ 1.5 时前缘失稳，论文图 15a）。

**b 约定**：浮力数 b 从 ``result.summary['buoyancy_number']`` 读取
（Zhang 2022 p.8 定义 b=(ρ₂−ρ₁)·g·d²/(μ₁·ŵ₀)，以被顶替液黏度 μ₁=μ_mud
标度；b>0 为密度稳定，即重顶替轻）。Bararpour 闭包严格采用几何平均
黏度标度（b_B = b·√m）；Tier 0 不做此重标定（m≈1 时两者一致），
作为已知近似在 notes 中披露。

宽/窄边补充判据（Pelipenko & Frigaard 2004 式 3.33/3.34）
------------------------------------------------------
对牛顿流体 Hele-Shaw 闭包，宽边（H_w=1+e）水泥指进速度与窄边
（H_n=1−e）泥浆（回退）指进速度可解析求出：

- 宽边失稳（水泥指进快于宽边泥浆基流 1）：m − b·√m·(1+e)²/3 > 1
- 窄边失稳（泥浆指进慢于窄边水泥基流 1，滞留/窜槽风险）：
  1/m + b·(1−e)²/(3√m) < 1

e→0 时与 Pelipenko 2004 式 4.12/4.13 的牛顿极限（(1−δ)/λ>1 与 λ+δ<1，
其中 λ=1/m、δ=b/(3√m)）精确一致；有限 e 的 (1±e)² 因子是本模块按
Hele-Shaw 浮力速度 ∝H² 缩放的扩展（Pelipenko 的 O(e) 摄动在 e→0 时
与该扩展一致）。

数值安全
--------
- c̄→0/1 端点：I1(c̄) 恒正（≥min(√m,1/√m)/3），所有除法以 1e-12 下限
  保护；q0/I3 在端点有良好定义（q0(0)=0、q0(1)=1、I3(0)=I3(1)=0）；
- 激波检测（上凹包络）失败时回退到纯数值微分 w_f=q0'+b·I3' 并在 notes 披露；
- 凹包络过 (0,0)、(1,1) 且 ≥F，故 w_f(1)≤1、Δw(1)=1−w_f(1)>0 恒成立，
  保证 Δw(0)<0 时必存在 c_critical。

局限性（Tier 0 牛顿近似）
------------------------
- HB 流体（屈服应力/幂律）的 I1/I2/I3/q0 需数值积分闭包，留 Tier 2
  （T2-1/T2-2）；对 HB 流体本模块以塑性黏度（或一致性系数）作等效牛顿
  近似，结果仅作筛查用途；
- b 的黏度标度差异（μ_mud vs 几何平均）未重标定；
- 激波处理假设左态 c̄=1、右态 c̄=0 的 Riemann 结构（前缘下游无预置
  混合带）。

使用示例
--------
>>> from cemdisp.diagnostics.muskat_regime import compute_muskat_regime
>>> mr = compute_muskat_regime(result, fluids, well_spec)
>>> print(mr.regime)          # "stable"
>>> print(mr.c_critical)      # 0.83
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Dict, List, Optional, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

from cemdisp.models2d.d2dga_flux import (
    d2dga_dispersion_function_I3,
    d2dga_dispersion_I1,
    d2dga_dispersion_I2,
    d2dga_flux_amplification,
)

if TYPE_CHECKING:  # 仅类型注解，避免运行时依赖 models2d 求解逻辑 / data 结构
    from cemdisp.data.fluid_spec import FluidSpec
    from cemdisp.data.well_spec import WellSpec
    from cemdisp.models2d.annulus_d2dga import AnnulusSimulationResult


Array = NDArray[np.float64]

REGIME_STABLE = "stable"
REGIME_PARTIAL = "partial_penetration"
REGIME_UNSTABLE = "unstable"

_GRAVITY_M_S2: float = 9.81
_EPS: float = 1.0e-12
_MIN_GRID: int = 51


@dataclass(frozen=True)
class MuskatRegimeResult:
    """Muskat 三 regime 判据结果（Bararpour & Frigaard 2025 式 3.8-3.13）。

    Attributes:
        regime: "stable" / "partial_penetration" / "unstable"
        c_critical: Δw=0 的最大浓度（指进停止穿透的临界浓度）；
            unstable（无根）时为 NaN（to_dict 中转 None）
        c_star: w_f(c̄)=1 的浓度（以平均泵速弥散的浓度，前缘体积中点）
        wide_side_unstable: 宽边指进失稳判据（Pelipenko 2004 式 3.33 牛顿形式）
        narrow_side_unstable: 窄边指进失稳判据（Pelipenko 2004 式 3.34 牛顿形式）
        delta_w_min: Δw 最小值（<0 表示存在指进无法穿透的浓度区间）
        delta_w_max: Δw 最大值（恒 >0，因 Δw(1)=1−w_f(1)>0）
        viscosity_ratio: 黏度比 m=μ_mud/μ_cement（与 d2dga_flux 约定一致）
        buoyancy_number: 浮力数 b（Zhang 2022 p.8 标度，b>0 密度稳定）
        eccentricity: 偏心度 e（geom['e'] 深度平均或 standoff 换算）
        shock_detected: 通量函数 F=q0+b·I3 是否存在凸区间（已按式 3.9 处理激波）
        c_bar_grid: 浓度网格（供 Δw 剖面绘图；include_profiles=False 时为空）
        delta_w_profile: Δw(c̄) 剖面（式 3.13）
        notes: 近似与披露备注
    """

    regime: str
    c_critical: float
    c_star: float
    wide_side_unstable: bool
    narrow_side_unstable: bool
    delta_w_min: float
    delta_w_max: float
    viscosity_ratio: float
    buoyancy_number: float
    eccentricity: float
    shock_detected: bool
    c_bar_grid: Tuple[float, ...] = ()
    delta_w_profile: Tuple[float, ...] = ()
    notes: Tuple[str, ...] = ()

    def to_dict(self) -> Dict[str, object]:
        """导出为可 JSON 序列化字典（NaN/Inf → None，tuple → list）。"""

        def _finite_or_none(x: float) -> Optional[float]:
            return x if math.isfinite(x) else None

        return {
            "regime": self.regime,
            "c_critical": _finite_or_none(self.c_critical),
            "c_star": _finite_or_none(self.c_star),
            "wide_side_unstable": self.wide_side_unstable,
            "narrow_side_unstable": self.narrow_side_unstable,
            "delta_w_min": _finite_or_none(self.delta_w_min),
            "delta_w_max": _finite_or_none(self.delta_w_max),
            "viscosity_ratio": _finite_or_none(self.viscosity_ratio),
            "buoyancy_number": _finite_or_none(self.buoyancy_number),
            "eccentricity": _finite_or_none(self.eccentricity),
            "shock_detected": self.shock_detected,
            "c_bar_grid": list(self.c_bar_grid),
            "delta_w_profile": [_finite_or_none(v) for v in self.delta_w_profile],
            "notes": list(self.notes),
        }


# ---------------------------------------------------------------------------
# 牛顿解析闭包（Bararpour & Frigaard 2025 式 2.24-2.26；I3 用 Zhang 2022 式 4.26）
# I1/I2 已提升至 d2dga_flux.d2dga_dispersion_I1/I2（T1-3a DRY）。
# 向后兼容别名（测试中仍 import _mean_mobility_I1 / _buoyant_mobility_I2）。
_mean_mobility_I1 = d2dga_dispersion_I1
_buoyant_mobility_I2 = d2dga_dispersion_I2
# ---------------------------------------------------------------------------


def _isotropic_flux_q0(c_bar: Array, m: float) -> Array:
    """各向同性通量 q0(c̄,m)（Bararpour 2025 式 2.26）= c̄·f(c̄,m)。

    复用 ``d2dga_flux_amplification`` 的 f(c,m)=[m·c²+1.5(1−c²)]/[m·c³+(1−c³)]，
    但**解除**求解器的 [0.5,2] 人工限幅与端点裁剪（Muskat 分析需要全区间
    真实通量；q0(0)=0、q0(1)=1 端点良定义）。
    """
    c = np.asarray(c_bar, dtype=float)
    f = d2dga_flux_amplification(
        c,
        m,
        min_fraction=0.0,
        max_fraction=1.0,
        min_amplification=-np.inf,
        max_amplification=np.inf,
    )
    return c * np.asarray(f, dtype=float)


def _front_flux(c_bar: Array, m: float, b: float) -> Array:
    """前缘通量函数 F(c̄) = q0(c̄) + b·I3(c̄)（式 3.7 的通量）。

    I3 用 ``d2dga_dispersion_function_I3``（Zhang 2022 式 4.26），解除端点
    裁剪（I3(0)=I3(1)=0 端点良定义）。
    """
    i3 = d2dga_dispersion_function_I3(c_bar, m, min_fraction=0.0, max_fraction=1.0)
    return _isotropic_flux_q0(c_bar, m) + b * np.asarray(i3, dtype=float)


def _finger_velocity(c_bar: Array, m: float, b: float) -> Array:
    """指进速度 w_finger(c̄)（式 3.12）。

    w_finger = I1(1)/I1(c̄) + b·I1(1)·[I2(c̄)/I1(c̄) + c̄ − 1]；
    利用了 I2(1)=0。I1 除零保护（I1 恒正，保护仅为防御性）。
    """
    c = np.asarray(c_bar, dtype=float)
    i1 = np.maximum(d2dga_dispersion_I1(c, m), _EPS)
    i1_at_1 = max(float(d2dga_dispersion_I1(1.0, m)), _EPS)
    i2 = d2dga_dispersion_I2(c, m)
    return i1_at_1 / i1 + b * i1_at_1 * (i2 / i1 + c - 1.0)


# ---------------------------------------------------------------------------
# 前缘速度（式 3.8 + 式 3.9 激波处理）
# ---------------------------------------------------------------------------


def _concave_majorant(c_bar: Array, flux: Array) -> Array:
    """上凹包络（Rankine-Hugoniot 等面积激波处理，式 3.9）。

    Andrew 单调链取点集 {(c̄_i, F_i)} 的上凸包：包络与 F 重合段对应弥散
    （稀疏波）区，弦段对应激波区间 [c̄m, c̄M]，弦斜率即激波速度
    w_s = [F(c̄M)−F(c̄m)]/(c̄M−c̄m)（式 3.9）。包络过端点 (0,F(0)) 与
    (1,F(1))，故 ∫₀¹ w_f dc̄ = F(1)−F(0) = 1 在激波处理后保持。
    """
    n = int(c_bar.size)
    hull: List[int] = []
    for i in range(n):
        while len(hull) >= 2:
            i0, i1 = hull[-2], hull[-1]
            cross = (c_bar[i1] - c_bar[i0]) * (flux[i] - flux[i0]) - (
                flux[i1] - flux[i0]
            ) * (c_bar[i] - c_bar[i0])
            if cross > 0.0:  # 非顺时针转折：中间点低于弦，弹出以保留上凹顶点
                hull.pop()
            else:
                break
        hull.append(i)
    idx = np.asarray(hull, dtype=int)
    return np.interp(c_bar, c_bar[idx], flux[idx]).astype(float, copy=False)


def _front_speed(c_bar: Array, flux: Array) -> Tuple[Array, bool, bool]:
    """前缘速度 w_f(c̄)（式 3.8，激波区间按式 3.9 处理）。

    Returns:
        (w_f, shock_detected, fell_back)：w_f 为上凹包络的导数（非递增）；
        shock_detected 表示通量存在凸区间（包络与 F 分离）；fell_back 表示
        包络构造失败已回退到纯数值微分。
    """
    try:
        majorant = _concave_majorant(c_bar, flux)
        tol = 1.0e-9 * max(1.0, float(np.max(np.abs(flux))))
        shock = bool(np.any(majorant - flux > tol))
        w_f = np.gradient(majorant, c_bar, edge_order=1)
        if not np.all(np.isfinite(w_f)):
            raise FloatingPointError("上凹包络导数含非有限值")
        return w_f.astype(float, copy=False), shock, False
    except Exception:
        w_f = np.gradient(flux, c_bar, edge_order=1)
        return w_f.astype(float, copy=False), False, True


def _mean_speed_concentration(c_bar: Array, w_f: Array) -> float:
    """c*：w_f(c̄)=1 的浓度（w_f 经激波处理后非递增，取最大此类浓度并线性插值）。"""
    above = np.nonzero(w_f >= 1.0)[0]
    if above.size == 0:
        return 0.0
    i = int(above[-1])
    if i >= c_bar.size - 1:
        return float(c_bar[-1])
    w0, w1 = float(w_f[i]), float(w_f[i + 1])
    denom = w0 - w1
    frac = 0.0 if abs(denom) < _EPS else (w0 - 1.0) / denom
    return float(c_bar[i] + float(np.clip(frac, 0.0, 1.0)) * (c_bar[i + 1] - c_bar[i]))


def _largest_zero_crossing(c_bar: Array, delta_w: Array) -> float:
    """Δw=0 的最大浓度根（符号变化处线性插值）；无根返回 NaN。"""
    roots: List[float] = []
    for i in range(int(c_bar.size) - 1):
        d0, d1 = float(delta_w[i]), float(delta_w[i + 1])
        if d0 == 0.0:
            roots.append(float(c_bar[i]))
        elif d0 * d1 < 0.0:
            roots.append(float(c_bar[i] - d0 * (c_bar[i + 1] - c_bar[i]) / (d1 - d0)))
    if not roots:
        return float("nan")
    return max(roots)


# ---------------------------------------------------------------------------
# Pelipenko & Frigaard 2004 式 3.33/3.34（牛顿 Hele-Shaw 形式）
# ---------------------------------------------------------------------------


def _wide_narrow_side_criteria(m: float, b: float, e: float) -> Tuple[bool, bool]:
    """宽/窄边局部指进判据（Pelipenko 2004 式 3.33/3.34 的牛顿形式）。

    - 宽边（H_w=1+e）：水泥指进速度 w_fw = m − b·√m·(1+e)²/3 > 1 → 失稳；
    - 窄边（H_n=1−e）：泥浆指进速度 w_fn = 1/m + b·(1−e)²/(3√m) < 1 → 失稳
      （泥浆被滞留，静态窜槽风险）。

    e→0 时退化为 Pelipenko 式 4.12/4.13 牛顿极限（(1−δ)/λ>1、λ+δ<1，
    λ=1/m、δ=b/(3√m)）。
    """
    e_c = float(np.clip(e, 0.0, 0.95))
    sq_m = math.sqrt(m)
    w_finger_wide = m - b * sq_m * (1.0 + e_c) ** 2 / 3.0
    w_finger_narrow = 1.0 / m + b * (1.0 - e_c) ** 2 / (3.0 * sq_m)
    return (w_finger_wide > 1.0, w_finger_narrow < 1.0)


# ---------------------------------------------------------------------------
# 主判据
# ---------------------------------------------------------------------------


def classify_muskat_regime(
    viscosity_ratio: float,
    buoyancy_number: float,
    eccentricity: float = 0.0,
    *,
    n_grid: int = 2001,
    include_profiles: bool = True,
) -> MuskatRegimeResult:
    """Muskat 三 regime 分类（纯函数核心，Bararpour 2025 式 3.8-3.13）。

    Args:
        viscosity_ratio: 黏度比 m=μ_mud/μ_cement（>0）
        buoyancy_number: 浮力数 b（Zhang 2022 p.8 标度；b>0 密度稳定）
        eccentricity: 偏心度 e∈[0,1)，仅用于宽/窄边补充判据
        n_grid: 浓度网格点数（≥51）
        include_profiles: 是否在结果中携带 c̄ 网格与 Δw 剖面

    Returns:
        MuskatRegimeResult

    Raises:
        ValueError: 输入参数非法（m≤0、b 非有限、网格过粗）
    """
    m = float(viscosity_ratio)
    b = float(buoyancy_number)
    e = float(eccentricity)
    if not math.isfinite(m) or m <= 0.0:
        raise ValueError(f"viscosity_ratio 必须为正的有限数值，当前值：{viscosity_ratio}")
    if not math.isfinite(b):
        raise ValueError(f"buoyancy_number 必须为有限数值，当前值：{buoyancy_number}")
    if not math.isfinite(e):
        raise ValueError(f"eccentricity 必须为有限数值，当前值：{eccentricity}")
    if int(n_grid) < _MIN_GRID:
        raise ValueError(f"n_grid 至少为 {_MIN_GRID}，当前值：{n_grid}")

    c = np.linspace(0.0, 1.0, int(n_grid))
    flux = _front_flux(c, m, b)
    w_f, shock_detected, fell_back = _front_speed(c, flux)
    w_finger = _finger_velocity(c, m, b)
    delta_w = w_finger - w_f  # 式 3.13

    delta_w_min = float(np.min(delta_w))
    delta_w_max = float(np.max(delta_w))
    c_star = _mean_speed_concentration(c, w_f)
    zero_tol = 1.0e-6 * max(1.0, float(np.max(np.abs(delta_w))))

    notes: List[str] = [
        "牛顿近似：I1/I2/q0 用 Bararpour 2025 式 2.24-2.26 解析闭包；HB 精确闭包留 Tier 2（T2-1/T2-2）。",
        "I3 采用 d2dga_flux 现有版本（Zhang 2022 式 4.26，分母 2m·[...]）；与 Bararpour 式 2.27（分母 12·[...]）的标度差异留待 Tier 1 核实。",
        "b 为 Zhang 2022 标度（μ_mud 归一）；Bararpour 几何平均黏度标度需 ×√m，Tier 0 未重标定。",
    ]
    if fell_back:
        notes.append("激波检测（上凹包络）失败，已回退到纯数值微分 w_f=q0'+b·I3'（式 3.8）。")
    elif shock_detected:
        notes.append("通量函数 F=q0+b·I3 存在凸区间：激波已按 Rankine-Hugoniot（式 3.9，上凹包络弦）处理。")

    # --- 三 regime 判据（论文 §3.3） ---
    if delta_w_min > zero_tol:
        # Δw(c̄)>0 对全部 c̄∈[0,1]：无 c_critical，指进穿透整个前缘
        regime = REGIME_UNSTABLE
        c_critical = float("nan")
    else:
        c_critical = _largest_zero_crossing(c, delta_w)
        if not math.isfinite(c_critical):
            # 退化情形：Δw 与零相切但未穿越（网格未捕获符号变化）
            c_critical = float(c[int(np.argmin(delta_w))])
            notes.append("Δw 与零相切但未穿越（退化），c_critical 取 Δw 最小值点。")
        below = c <= c_critical + _EPS
        positive_below = bool(np.any(delta_w[below] > zero_tol))
        if positive_below:
            # c_critical 以下仍存在 Δw>0 区间（多次穿越），不满足稳定定义
            regime = REGIME_PARTIAL
        elif c_star < c_critical:
            regime = REGIME_STABLE
        else:
            regime = REGIME_PARTIAL

    wide_side_unstable, narrow_side_unstable = _wide_narrow_side_criteria(m, b, e)

    return MuskatRegimeResult(
        regime=regime,
        c_critical=c_critical,
        c_star=c_star,
        wide_side_unstable=wide_side_unstable,
        narrow_side_unstable=narrow_side_unstable,
        delta_w_min=delta_w_min,
        delta_w_max=delta_w_max,
        viscosity_ratio=m,
        buoyancy_number=b,
        eccentricity=e,
        shock_detected=shock_detected,
        c_bar_grid=tuple(float(v) for v in c) if include_profiles else (),
        delta_w_profile=tuple(float(v) for v in delta_w) if include_profiles else (),
        notes=tuple(notes),
    )


# ---------------------------------------------------------------------------
# 从模拟结果提取输入（Layer 5 诊断入口）
# ---------------------------------------------------------------------------


def _role_name(fluid: object) -> str:
    """提取流体角色名字符串（兼容 FluidRole 枚举与纯字符串 mock）。"""
    role = getattr(fluid, "role", "")
    return str(getattr(role, "value", role)).strip().lower()


def _pick_mud_and_cement(fluids: Sequence["FluidSpec"]) -> Tuple["FluidSpec", "FluidSpec"]:
    """从流体列表选取钻井液（role=mud）与水泥浆（优先 lead，其次 tail）。

    与 ``annulus_d2dga._pick_fluids`` 的角色约定一致；顶替液密度/黏度取
    领浆优先（与求解器 buoyancy_number 计算中 rho_displacing 的选取一致）。
    """
    mud = next((f for f in fluids if _role_name(f) == "mud"), None)
    lead = next((f for f in fluids if _role_name(f) == "lead"), None)
    tail = next((f for f in fluids if _role_name(f) == "tail"), None)
    cement = lead if lead is not None else tail
    if mud is None:
        raise ValueError("fluids 中未找到 role='mud' 的钻井液")
    if cement is None:
        raise ValueError("fluids 中未找到 role='lead'/'tail' 的水泥浆")
    return mud, cement


def _effective_viscosity_pa_s(fluid: "FluidSpec", label: str) -> float:
    """等效牛顿黏度：优先塑性黏度，其次一致性系数（幂律/HB 的牛顿近似）。"""
    mu = getattr(fluid, "plastic_viscosity_pa_s", None)
    if mu is not None and math.isfinite(float(mu)) and float(mu) > 0.0:
        return float(mu)
    k = getattr(fluid, "consistency_k", None)
    if k is not None and math.isfinite(float(k)) and float(k) > 0.0:
        return float(k)
    raise ValueError(f"{label}缺少可用黏度（plastic_viscosity_pa_s 与 consistency_k 均不可用）")


def _extract_buoyancy_number(
    result: "AnnulusSimulationResult",
    mud: "FluidSpec",
    cement: "FluidSpec",
    mu_mud: float,
    mean_velocity_m_s: Optional[float],
) -> Tuple[float, str]:
    """浮力数 b：优先读 result.summary['buoyancy_number']，否则按定义重算。

    重算公式与 ``annulus_d2dga._compute_buoyancy_number`` 一致（Zhang 2022
    p.8）：b = (ρ_cement−ρ_mud)·g·(gap/2)²/(μ_mud·ŵ₀)，gap 取 geom['b'] 均值。
    """
    summary = getattr(result, "summary", None)
    if isinstance(summary, dict):
        b_val = summary.get("buoyancy_number")
        if b_val is not None and math.isfinite(float(b_val)):
            return float(b_val), "result.summary['buoyancy_number']"
    geom = getattr(result, "geom", None)
    if mean_velocity_m_s is None or not isinstance(geom, dict) or "b" not in geom:
        raise ValueError(
            "result.summary 中无有限 'buoyancy_number'，且缺少回退重算所需输入"
            "（需 result.geom['b'] 与 mean_velocity_m_s 参数）"
        )
    gap_m = float(np.mean(np.asarray(geom["b"], dtype=float)))
    d_half = max(gap_m / 2.0, 1.0e-6)
    w0 = max(float(mean_velocity_m_s), 1.0e-6)
    rho_mud = float(getattr(mud, "density_kg_m3"))
    rho_cement = float(getattr(cement, "density_kg_m3"))
    b = (rho_cement - rho_mud) * _GRAVITY_M_S2 * d_half**2 / max(mu_mud * w0, 1.0e-9)
    return float(b), "由流体/几何重算（Zhang 2022 p.8 定义）"


def _extract_eccentricity(
    result: "AnnulusSimulationResult",
    well_spec: Optional["WellSpec"],
) -> Tuple[float, str]:
    """偏心度 e：优先 geom['e'] 深度平均，其次 well_spec.standoff_profile 换算。"""
    geom = getattr(result, "geom", None)
    if isinstance(geom, dict) and "e" in geom:
        return float(np.mean(np.asarray(geom["e"], dtype=float))), "result.geom['e'] 深度平均"
    if well_spec is not None:
        profile = getattr(well_spec, "standoff_profile", None) or ()
        values = [float(p.value) for p in profile if math.isfinite(float(getattr(p, "value", float("nan"))))]
        if values:
            e = float(np.clip(1.0 - float(np.mean(values)), 0.0, 0.95))
            return e, "well_spec.standoff_profile 均值换算 e=1−standoff"
    return 0.0, "无偏心信息，取 e=0（同心近似）"


def compute_muskat_regime(
    result: "AnnulusSimulationResult",
    fluids: Sequence["FluidSpec"],
    well_spec: Optional["WellSpec"] = None,
    *,
    n_grid: int = 2001,
    mean_velocity_m_s: Optional[float] = None,
    include_profiles: bool = True,
) -> MuskatRegimeResult:
    """从环空二维顶替结果计算 Muskat 三 regime 诊断（T0-2 主入口）。

    纯后处理：只消费 ``AnnulusSimulationResult`` 的 summary/geom 与流体物性，
    不调用任何求解逻辑。

    Args:
        result: AnnulusSimulationResult，其 summary 含 'buoyancy_number' 键、
            geom 含 'e'（偏心度）与 'b'（间隙宽度，仅回退重算 b 时用）
        fluids: 流体规格序列（需含 role='mud' 与 role='lead'/'tail'）
        well_spec: 井规格（可选；仅用于 geom['e'] 缺失时的 standoff 回退）
        n_grid: 浓度网格点数
        mean_velocity_m_s: 平均泵速 [m/s]（可选；仅 summary 缺 buoyancy_number
            时用于回退重算 b）
        include_profiles: 是否携带 Δw 剖面

    Returns:
        MuskatRegimeResult

    Raises:
        ValueError: 流体角色缺失、黏度不可用或 b 无法获得时
    """
    mud, cement = _pick_mud_and_cement(fluids)
    mu_mud = _effective_viscosity_pa_s(mud, "钻井液")
    mu_cement = _effective_viscosity_pa_s(cement, "水泥浆")
    m = mu_mud / mu_cement
    b, b_source = _extract_buoyancy_number(result, mud, cement, mu_mud, mean_velocity_m_s)
    e, e_source = _extract_eccentricity(result, well_spec)

    out = classify_muskat_regime(m, b, e, n_grid=n_grid, include_profiles=include_profiles)
    extract_note = (
        f"输入提取：m=μ_mud/μ_cement={m:.4f}（μ_mud={mu_mud:.4g} Pa·s，"
        f"μ_cement={mu_cement:.4g} Pa·s）；b={b:.6g}（{b_source}）；e={e:.4f}（{e_source}）。"
    )
    return replace(out, notes=out.notes + (extract_note,))
