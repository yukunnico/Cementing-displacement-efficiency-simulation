"""Zhang & Frigaard (2022, JFM 947:A32) 基准算例对照运行器。

对照对象
--------
Zhang R. & Frigaard I. A., "Primary cementing of vertical wells: displacement and
dispersion in narrow eccentric annuli", J. Fluid Mech. 947 (2022) A32.
DOI: 10.1017/jfm.2022.626

* 几何（§2.3）：环空长 L = 4.8 m，外半径 r_o = 22.23 mm，内半径 r_i = 17.46 mm，
  ⇒ 半间隙 d = (r_o − r_i)/2 = 2.385 mm；竖直井。
* Table 1：10 个算例的 (ŵ0, Q̂0, ρ̂1, μ̂1, ρ̂2, μ̂2) 物性；
* Table 2：10 个算例的 (e, Re, m, b) 无量纲参数；
* Table 3：10 个算例的 t_br 与 η_E（D2DGA 列 / 3-D 列）。
* 口径（§5）：t_br = t̂_br·ŵ0/L，η_E = t = 1.2·L/ŵ0 时刻的环空已顶替体积分数。

诚实定位（务必随结果一起报告）
------------------------------
论文 Table 3 的 "D2DGA" 列**本身也是数值解**（FCT 有限体积），本模型是半拉格朗日
平流（无人工扩散项——Z&F22 p.11 "we have no diffusive terms"，2026-09-14 Task 7
已删除旧 ``dispersion_*`` 拉普拉斯弥散；弥散由 q₀ + I₃ 分层通量闭合承载），
二者共享同一个理论源头（Zhang & Frigaard 2022 的 D2DGA 闭包）但离散格式不同。
因此本对照是 **同源算例的半定量交叉验证**，不能写成"数值正确性验证通过"。
容差先验：η_E ±0.05、t_br ±0.10。

额外记录项（论文未有的实现差异）
--------------------------------
1. 偏心度口径：``e = 1 − standoff`` 按文献口径 e∈[0,1) 直取（Pelipenko04 (2.1)）。
   旧 ``e_clip`` 硬截断已于 2026-09-15 Task 10 移除，本 runner 不再传弃用形参
   ``e_clip_max``（Task 12 起，构造零弃用警告），e ≥ 0.6 的 case 1-4 几何按
   论文原值进入求解器。
2. ``annulus_d2dga.py:921`` ``correction = np.clip(Δρ·I2/I1, −0.5, 0.5)``：论文式 4.24
   无此裁剪。本 runner 用透明 numpy 代理统计其触发情况（见 ``_clip_probe``），
   并把计数写入结果表 ``clip_triggered`` 列。
3. 论文用 ŵ0 无量纲化；本模型按 Q 驱动，环空截面平均速度 = Q/A。二者在 Table 1 中
   逐例一致（case 9 除外，见 ``ZHANG2022_CASES`` 脚注）。
"""

from __future__ import annotations

import csv
import math
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Callable, Iterator, Mapping

import numpy as np

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.well_spec import DepthValuePoint, WellSpec
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import AnnulusInletState
from cemdisp.models2d.d2dga_flux import d2dga_dispersion_I1, d2dga_dispersion_I2


# ---------------------------------------------------------------------------
# 论文常量
# ---------------------------------------------------------------------------
ZHANG2022_L_M: float = 4.8
ZHANG2022_OUTER_RADIUS_M: float = 22.23e-3
ZHANG2022_INNER_RADIUS_M: float = 17.46e-3
ZHANG2022_D_HALF_M: float = (ZHANG2022_OUTER_RADIUS_M - ZHANG2022_INNER_RADIUS_M) / 2.0
ZHANG2022_G: float = 9.81

_CEMDISP_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = _CEMDISP_ROOT.parent
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "results" / "基准算例对照_2026-09-10"


def _header(labels: tuple[str, ...]) -> str:
    """生成 Markdown 表头分隔行。"""
    return "| " + " | ".join(labels) + " |\n|" + "|".join("---" for _ in labels) + "|"


# ---------------------------------------------------------------------------
# Table 1 + Table 2 + Table 3（逐字核对 参考文档/_extracted_pdfs/
#                  Zhang和Frigaard_-_2022_-_Primary_cementing_of_vertical_wells_.txt）
# ---------------------------------------------------------------------------
# Table 1 列：ŵ0 (m/s), Q̂0 (m³/s), ρ̂1 (kg/m³), μ̂1 (Pa·s), ρ̂2 (kg/m³), μ̂2 (Pa·s)
# Table 2 列：e, Re, m, b
# Table 3 列：t_br(D2DGA), t_br(3-D), η_E(D2DGA), η_E(3-D)
#
# ⚠️ 论文 Table 1 case 9 的 Q̂0 = 2.38×10⁻⁴ 与其 ŵ0 = 0.04 不自洽（ŵ0·A = 2.38×10⁻⁵，
#    差 10 倍），且与 Table 2 给定的 b = 100（= Δρ·g·d²/(μ1·ŵ0)）矛盾——按 Q̂0 反推
#    会得到 b = 10。判定为论文表格笔误，故 ``Q0`` 统一取 ŵ0·A_annulus（自洽值），
#    原始打印值另存 ``Q0_printed`` 备查。
ZHANG2022_CASES: tuple[dict[str, Any], ...] = (
    {"case_id": 1, "e": 0.8, "Re": 20, "m": 0.2, "b": -50.0,
     "w0": 0.032, "Q0_printed": 1.90e-5, "rho1": 1000.0, "mu1": 0.004,
     "rho2": 885.31, "mu2": 0.02,
     "t_br_d2dga": 0.44, "t_br_3d": 0.33, "eta_e_d2dga": 0.66, "eta_e_3d": 0.61},
    {"case_id": 2, "e": 0.8, "Re": 20, "m": 0.2, "b": 100.0,
     "w0": 0.032, "Q0_printed": 1.90e-5, "rho1": 1000.0, "mu1": 0.004,
     "rho2": 1229.38, "mu2": 0.02,
     "t_br_d2dga": 0.95, "t_br_3d": 0.95, "eta_e_d2dga": 0.95, "eta_e_3d": 0.95},
    {"case_id": 3, "e": 0.6, "Re": 20, "m": 0.2, "b": 10.0,
     "w0": 0.032, "Q0_printed": 1.90e-5, "rho1": 1000.0, "mu1": 0.004,
     "rho2": 1022.94, "mu2": 0.02,
     "t_br_d2dga": 0.79, "t_br_3d": 0.67, "eta_e_d2dga": 0.92, "eta_e_3d": 0.91},
    {"case_id": 4, "e": 0.6, "Re": 20, "m": 0.2, "b": 1000.0,
     "w0": 0.008, "Q0_printed": 4.76e-6, "rho1": 1000.0, "mu1": 0.001,
     "rho2": 1143.37, "mu2": 0.005,
     "t_br_d2dga": 0.99, "t_br_3d": 0.98, "eta_e_d2dga": 1.00, "eta_e_3d": 0.98},
    {"case_id": 5, "e": 0.4, "Re": 20, "m": 0.2, "b": 100.0,
     "w0": 0.032, "Q0_printed": 1.90e-5, "rho1": 1000.0, "mu1": 0.004,
     "rho2": 1229.38, "mu2": 0.02,
     "t_br_d2dga": 0.95, "t_br_3d": 0.93, "eta_e_d2dga": 0.97, "eta_e_3d": 0.95},
    {"case_id": 6, "e": 0.4, "Re": 20, "m": 5.0, "b": 100.0,
     "w0": 0.04, "Q0_printed": 2.38e-5, "rho1": 1000.0, "mu1": 0.005,
     "rho2": 1358.41, "mu2": 0.001,
     "t_br_d2dga": 0.93, "t_br_3d": 0.94, "eta_e_d2dga": 0.93, "eta_e_3d": 0.96},
    {"case_id": 7, "e": 0.2, "Re": 20, "m": 0.5, "b": 10.0,
     "w0": 0.032, "Q0_printed": 1.90e-5, "rho1": 1000.0, "mu1": 0.004,
     "rho2": 1022.94, "mu2": 0.008,
     "t_br_d2dga": 0.78, "t_br_3d": 0.70, "eta_e_d2dga": 0.90, "eta_e_3d": 0.90},
    {"case_id": 8, "e": 0.2, "Re": 20, "m": 2.0, "b": 10.0,
     "w0": 0.08, "Q0_printed": 4.76e-5, "rho1": 1000.0, "mu1": 0.010,
     "rho2": 1143.37, "mu2": 0.005,
     "t_br_d2dga": 0.78, "t_br_3d": 0.70, "eta_e_d2dga": 0.84, "eta_e_3d": 0.86},
    {"case_id": 9, "e": 0.1, "Re": 100, "m": 0.2, "b": 100.0,
     "w0": 0.04, "Q0_printed": 2.38e-4, "rho1": 1000.0, "mu1": 0.001,
     "rho2": 1071.68, "mu2": 0.005,
     "t_br_d2dga": 0.97, "t_br_3d": 0.96, "eta_e_d2dga": 0.97, "eta_e_3d": 0.95},
    {"case_id": 10, "e": 0.1, "Re": 1000, "m": 0.2, "b": 100.0,
     "w0": 0.4, "Q0_printed": 2.38e-4, "rho1": 1000.0, "mu1": 0.001,
     "rho2": 1716.83, "mu2": 0.005,
     "t_br_d2dga": 0.97, "t_br_3d": 0.96, "eta_e_d2dga": 0.97, "eta_e_3d": 0.95},
)

# Q0 = ŵ0 · A_annulus（见表头 case 9 脚注）。A 由论文半径按 π(r_o²−r_i²) 计算。
ZHANG2022_ANNULUS_VOLUME_M3: float = (
    math.pi * (ZHANG2022_OUTER_RADIUS_M**2 - ZHANG2022_INNER_RADIUS_M**2) * ZHANG2022_L_M
)
ZHANG2022_CASES = tuple(
    {**case, "Q0": case["w0"] * ZHANG2022_ANNULUS_VOLUME_M3 / ZHANG2022_L_M}
    for case in ZHANG2022_CASES
)
ZHANG2022_CASE_BY_ID: Mapping[int, dict[str, Any]] = {c["case_id"]: c for c in ZHANG2022_CASES}

# 论文中 e 被用作 "偏心度"，模型内部 standoff = 1 − e。
# 环空入口深度须 > 0（WellSpec 约束），故井段取 [1.0, 1.0 + L]。
_TOP_MD_M: float = 1.0
# 求解器内部把 hole/od 视作**直径**（clearance = hole − od，体积按 π(Do²−Di²)/4），
# 故这里把论文半径换算成直径输入。
_HOLE_DIAMETER_MM: float = 2.0 * ZHANG2022_OUTER_RADIUS_M * 1000.0   # 44.46 mm
_LINER_OD_MM: float = 2.0 * ZHANG2022_INNER_RADIUS_M * 1000.0        # 34.92 mm

# 基准尺度下（w0 最大 0.4 m/s、nz=140 → ds≈0.0345 m）求解器默认 dt_min=0.1 s 会
# 使 CFL = w_max·dt/ds > 1（case 10 实测 ~1.4）。基准统一放宽到 0.02 s 让 CFL 控制步长。
_BENCHMARK_DT_MIN: float = 0.02


# ---------------------------------------------------------------------------
# 算例包装
# ---------------------------------------------------------------------------
def _coerce(case: int | Mapping[str, Any]) -> dict[str, Any]:
    """接受 case_id 或 case dict，返回 case dict。"""
    if isinstance(case, int):
        return ZHANG2022_CASE_BY_ID[case]
    return dict(case)  # type: ignore[arg-type]


def case_mean_velocity_m_s(case: int | Mapping[str, Any]) -> float:
    """论文 Table 1 的 ŵ0（环空截面平均速度，m/s）。"""
    return float(_coerce(case)["w0"])


def case_annulus_volume_m3() -> float:
    """论文几何的环空体积 π(r_o²−r_i²)·L（m³）。"""
    return ZHANG2022_ANNULUS_VOLUME_M3


def case_viscosity_ratio(case: int | Mapping[str, Any]) -> float:
    """论文定义 m = μ̂1/μ̂2。"""
    c = _coerce(case)
    return float(c["mu1"]) / float(c["mu2"])


def case_buoyancy_number(case: int | Mapping[str, Any]) -> float:
    """论文定义 b = (ρ̂2 − ρ̂1)·ĝ·d̂²/(μ̂1·ŵ0)，d̂ = 半间隙 2.385 mm。"""
    c = _coerce(case)
    w0 = case_mean_velocity_m_s(c)
    return float(
        (c["rho2"] - c["rho1"]) * ZHANG2022_G * ZHANG2022_D_HALF_M**2 / (float(c["mu1"]) * w0)
    )


def build_case_well_spec(case: int | Mapping[str, Any]) -> WellSpec:
    """构造论文几何的 WellSpec（竖直、等径、恒定偏心度）。"""
    c = _coerce(case)
    top, bottom = _TOP_MD_M, _TOP_MD_M + ZHANG2022_L_M
    standoff = 1.0 - float(c["e"])
    return WellSpec(
        well_name=f"Zhang2022_case{int(c['case_id'])}",
        top_md_m=top,
        bottom_md_m=bottom,
        shoe_md_m=bottom,
        hole_diameter_profile=(DepthValuePoint(top, _HOLE_DIAMETER_MM),
                               DepthValuePoint(bottom, _HOLE_DIAMETER_MM)),
        liner_od_profile=(DepthValuePoint(top, _LINER_OD_MM),
                          DepthValuePoint(bottom, _LINER_OD_MM)),
        inclination_profile=(DepthValuePoint(top, 0.0), DepthValuePoint(bottom, 0.0)),
        standoff_profile=(DepthValuePoint(top, standoff), DepthValuePoint(bottom, standoff)),
        notes=(f"Zhang & Frigaard (2022) JFM 947:A32 基准算例 case {int(c['case_id'])}",),
    )


def build_case_fluids(case: int | Mapping[str, Any]) -> tuple[FluidSpec, FluidSpec]:
    """构造 Table 1 的牛顿流体对：(钻井液 ρ̂1/μ̂1, 水泥浆 ρ̂2/μ̂2)。

    论文中流体 2（顶替液，水泥）顶替流体 1（被顶替液，钻井液）；
    m = μ̂1/μ̂2 与 b = (ρ̂2−ρ̂1)·ĝ·d̂²/(μ̂1·ŵ0) 均按此约定。
    """
    c = _coerce(case)
    mud = FluidSpec(f"case{int(c['case_id'])}_mud", FluidRole.MUD,
                    float(c["rho1"]), plastic_viscosity_pa_s=float(c["mu1"]))
    cement = FluidSpec(f"case{int(c['case_id'])}_cement", FluidRole.TAIL,
                       float(c["rho2"]), plastic_viscosity_pa_s=float(c["mu2"]))
    return mud, cement


def build_case_solver(
    case: int | Mapping[str, Any],
    *,
    nz: int = 140,
    ny: int = 40,
    total_t: float | None = None,
    **overrides: Any,
) -> AnnulusD2DGASolver:
    """构造基准算例求解器。

    显式设置：
    * ``dt_min = 0.02``：基准速度尺度下让 CFL 而非 dt_min 控制步长；
    * ``total_t = 2·L/ŵ0``：论文"泵注两个环空体积"的停止时刻。

    偏心度不再传 ``e_clip_max``（Task 10 已移除 e_clip 硬截断，e = 1−standoff 按
    Pelipenko04 (2.1) 文献口径 e∈[0,1) 直取；Task 12 起停止传该弃用形参，
    构造零弃用警告，论文算例 e 最大 0.8 按原值进入几何）。
    """
    c = _coerce(case)
    if total_t is None:
        total_t = 2.0 * ZHANG2022_L_M / case_mean_velocity_m_s(c)
    kwargs: dict[str, Any] = {
        "nz": nz,
        "ny": ny,
        "total_t": total_t,
        "dt_min": _BENCHMARK_DT_MIN,
    }
    kwargs.update(overrides)
    return AnnulusD2DGASolver(**kwargs)


def build_inlet_provider(case: int | Mapping[str, Any]) -> Callable[[float], AnnulusInletState]:
    """环空入口恒为纯水泥浆、恒定排量 Q̂0（绕开 1D 套管内耦合求解器）。"""
    q0 = float(_coerce(case)["Q0"])

    def _provider(time_s: float) -> AnnulusInletState:
        return AnnulusInletState(time_s, q0, "bench", (("tail", 1.0),))

    return _provider


# ---------------------------------------------------------------------------
# 论文式 4.24 硬裁剪的统计（annulus_d2dga.py:921）
# ---------------------------------------------------------------------------
class _ClipCounter:
    """``np.clip(x, −0.5, 0.5)`` 调用统计容器。"""

    def __init__(self) -> None:
        self.calls: int = 0
        self.elements: int = 0
        self.clipped_elements: int = 0
        self.max_abs: float = 0.0


class _NumpyClipProbe:
    """透明代理 numpy：仅拦截 ``np.clip(x, −0.5, 0.5)``（论文式 4.24 之外的硬裁剪）。

    代理对其它属性/调用完全透传，返回值与未插桩时逐位一致——只统计、不改行为。
    """

    def __init__(self, real_np: Any, counter: _ClipCounter) -> None:
        self._np = real_np
        self._counter = counter

    def __getattr__(self, name: str) -> Any:
        return getattr(self._np, name)

    def clip(self, a: Any, a_min: Any = None, a_max: Any = None, *args: Any, **kwargs: Any) -> Any:
        if (
            not args
            and not kwargs
            and isinstance(a_min, float)
            and isinstance(a_max, float)
            and a_min == -0.5
            and a_max == 0.5
        ):
            arr = self._np.asarray(a, dtype=float)
            self._counter.calls += 1
            self._counter.elements += int(arr.size)
            if arr.size:
                self._counter.clipped_elements += int(self._np.count_nonzero(self._np.abs(arr) > 0.5))
                self._counter.max_abs = max(self._counter.max_abs, float(self._np.max(self._np.abs(arr))))
        return self._np.clip(a, a_min, a_max, *args, **kwargs)


@contextmanager
def _clip_probe(module: Any) -> Iterator[_ClipCounter]:
    """临时替换模块级 ``np`` 以统计硬裁剪触发情况。"""
    counter = _ClipCounter()
    real_np = module.np
    module.np = _NumpyClipProbe(real_np, counter)
    try:
        yield counter
    finally:
        module.np = real_np


def clip_trigger_bound(case: int | Mapping[str, Any]) -> float:
    """式 4.24 硬裁剪触发的理论上界 max_c |Δρ(g/cc)·I2(c,m)/I1(c,m)|。

    与具体浓度场无关，故可作为"这 10 个算例是否可能触发 ±0.5 裁剪"的判据。
    """
    c = _coerce(case)
    m = case_viscosity_ratio(c)
    delta_rho_gcc = abs(float(c["rho2"]) - float(c["rho1"])) / 1000.0
    grid = np.linspace(0.0, 1.0, 2001)
    i1 = np.asarray(d2dga_dispersion_I1(grid, m), dtype=float)
    i2 = np.asarray(d2dga_dispersion_I2(grid, m), dtype=float)
    return float(np.max(np.abs(delta_rho_gcc * i2 / np.maximum(i1, 1e-12))))


# ---------------------------------------------------------------------------
# 论文口径的指标换算
# ---------------------------------------------------------------------------
def annulus_volume_of(well_spec: WellSpec) -> float:
    """按 well_spec 井径/衬管剖面算物理环空体积（π(D_o²−D_i²)/4 沿深度积分）。"""
    hole_md = np.array([p.depth_md_m for p in well_spec.hole_diameter_profile], dtype=float)
    hole = np.array([p.value for p in well_spec.hole_diameter_profile], dtype=float)
    if well_spec.liner_od_profile:
        od_md = np.array([p.depth_md_m for p in well_spec.liner_od_profile], dtype=float)
        od = np.array([p.value for p in well_spec.liner_od_profile], dtype=float)
    else:
        od_md = hole_md
        od = np.full_like(hole, float(well_spec.liner_od_mm or 0.0))
    area = np.pi * ((hole / 1000.0) ** 2 - (np.interp(hole_md, od_md, od) / 1000.0) ** 2) / 4.0
    return float(np.trapezoid(area, x=hole_md))


def peak_volume_ratio(
    result: Any,
    well_spec: WellSpec,
    q_m3s: float,
    warmup_frac: float = 0.02,
) -> float:
    """max_{t>warmup} (V_ann·bulk_fill) / (Q·t)：环空内水泥体积对累计注入体积的峰值比。

    >1 即"凭空造出水泥"——不受出口边界/突破口径影响，是质量守恒的硬判据。
    """
    v_ann = annulus_volume_of(well_spec)
    times = result.metrics["time_s"].to_numpy(dtype=float)
    fill = result.metrics["bulk_cement_fill"].to_numpy(dtype=float)
    t_end = float(times[-1])
    mask = (times > warmup_frac * t_end) & (times > 0.0)
    if not bool(mask.any()):
        return float("nan")
    return float(np.max(v_ann * fill[mask] / (q_m3s * times[mask])))


def paper_eta_e(result: Any, well_spec: WellSpec, w0: float) -> float:
    """论文 η_E：t = 1.2·L/ŵ0 时刻环空已顶替体积分数。

    ``result.metrics['bulk_cement_fill']`` ≡ ∫∫b·c / ∫∫b，正是"已顶替体积 / 环空体积"，
    故直接在该时刻做线性插值。
    """
    length = float(well_spec.bottom_md_m - well_spec.top_md_m)
    t_meas = 1.2 * length / w0
    times = result.metrics["time_s"].to_numpy(dtype=float)
    fill = result.metrics["bulk_cement_fill"].to_numpy(dtype=float)
    return float(np.interp(t_meas, times, fill))


def paper_t_br(
    result: Any,
    well_spec: WellSpec,
    w0: float,
    q_m3s: float | None = None,
    method: str = "front",
    tol_frac: float = 1e-3,
) -> tuple[float | None, float | None]:
    """突破时间：顶替液前缘首次到达环空出口的时刻。

    返回 ``(t_br_秒, t_br_量纲一)``；窗口内未突破时返回 ``(None, None)``。

    两种口径：

    * ``method="front"``（默认，对应论文"顶替液首次流出环空"）：取 metrics 的
      ``front_wide_m/front_narrow_m/front_mid_m``（求解器 ``_front``，浓度 0.5 等值线）
      中最早到达出口者。不需要 Q。
    * ``method="deficit"``（论文图 4 的 V₂(t) 脱离理论注入直线口径）：累计注入体积
      Q·t 与环空内水泥体积之差首次超过 ``tol_frac·V_ann``。需要 ``q_m3s``。

    ⚠️ 本模型 ``method="deficit"`` 口径**被 D2DGA 通量放大的体积不守恒污染**
    （见 ``mass_conservation_error``）：deficit 先为负（多出 28% 体积、后转为正），
    拐点被推迟到环空几乎填满之后，故默认用 front 口径。
    """
    length = float(well_spec.bottom_md_m - well_spec.top_md_m)
    times = result.metrics["time_s"].to_numpy(dtype=float)
    if method == "front":
        fronts = result.metrics[["front_wide_m", "front_narrow_m", "front_mid_m"]].to_numpy(dtype=float)
        reached = np.max(fronts, axis=1) >= length - 1e-9
        idx = np.flatnonzero(reached)
        if idx.size == 0:
            return None, None
        t_br = float(times[idx[0]])
        return t_br, t_br * w0 / length
    if method == "deficit":
        if q_m3s is None:
            raise ValueError("deficit 口径需要 q_m3s")
        v_ann = annulus_volume_of(well_spec)
        fill = result.metrics["bulk_cement_fill"].to_numpy(dtype=float)
        deficit = q_m3s * times - v_ann * fill
        idx = np.flatnonzero(deficit > tol_frac * v_ann)
        if idx.size == 0:
            return None, None
        t_br = float(times[idx[0]])
        return t_br, t_br * w0 / length
    raise ValueError(f"未知突破时间口径 {method!r}（可选 'front' / 'deficit'）")


def mass_conservation_error(
    result: Any,
    well_spec: WellSpec,
    q_m3s: float,
    t_end_s: float,
    warmup_frac: float = 0.02,
) -> float:
    """[0, t_end_s] 窗口内 |V_ann·bulk_fill − Q·t| / (Q·t) 的最大值。

    论文 §2.4 的质量守恒锚点：突破前环空内顶替液体积必须等于累计注入体积。

    ``warmup_frac`` 用于跳过起步瞬态：半拉格朗日入口在第一步把整格浓度置为入口值，
    当 Q·dt 不足一格时会在首个时间步给出 O(10%) 的假缺口。默认跳过窗口前 2%。
    """
    v_ann = annulus_volume_of(well_spec)
    times = result.metrics["time_s"].to_numpy(dtype=float)
    fill = result.metrics["bulk_cement_fill"].to_numpy(dtype=float)
    mask = (times > warmup_frac * t_end_s) & (times <= t_end_s)
    if not bool(mask.any()):
        return float("nan")
    injected = q_m3s * times[mask]
    contained = v_ann * fill[mask]
    return float(np.max(np.abs(contained - injected) / np.maximum(injected, 1e-30)))


# ---------------------------------------------------------------------------
# 单算例运行
# ---------------------------------------------------------------------------
def run_case(
    case: int | Mapping[str, Any],
    *,
    nz: int = 140,
    ny: int = 40,
    total_t: float | None = None,
    solver_overrides: Mapping[str, Any] | None = None,
    eta_tol: float = 0.05,
    t_br_tol: float = 0.10,
    keep_result: bool = False,
) -> dict[str, Any]:
    """运行单个基准算例并返回对照行（含额外诊断字段）。"""
    import cemdisp.models2d.annulus_d2dga as _ann

    c = _coerce(case)
    well_spec = build_case_well_spec(c)
    fluids = build_case_fluids(c)
    solver = build_case_solver(c, nz=nz, ny=ny, total_t=total_t, **dict(solver_overrides or {}))
    provider = build_inlet_provider(c)
    w0 = case_mean_velocity_m_s(c)

    with _clip_probe(_ann) as counter:
        result = solver.run(well_spec, fluids, provider)

    eta_model = paper_eta_e(result, well_spec, w0)
    t_br_s, t_br_model = paper_t_br(result, well_spec, w0, method="front")
    t_br_deficit_s, t_br_deficit = paper_t_br(
        result, well_spec, w0, method="deficit", q_m3s=float(c["Q0"])
    )
    if t_br_s is not None:
        mc_err = mass_conservation_error(result, well_spec, float(c["Q0"]), t_end_s=0.9 * t_br_s)
    else:
        # 窗口内未突破 → 全程都应守恒
        mc_err = mass_conservation_error(result, well_spec, float(c["Q0"]),
                                         t_end_s=float(result.metrics["time_s"].iloc[-1]))

    row: dict[str, Any] = {
        # ---- 论文要求的 12 列 ----
        "case_id": int(c["case_id"]),
        "e": float(c["e"]),
        "Re": int(c["Re"]),
        "m": float(c["m"]),
        "b": float(c["b"]),
        "t_br_模型": float("nan") if t_br_model is None else round(t_br_model, 4),
        "t_br_论文D2DGA": float(c["t_br_d2dga"]),
        "t_br_论文3D": float(c["t_br_3d"]),
        "eta_E_模型": round(eta_model, 4),
        "eta_E_论文D2DGA": float(c["eta_e_d2dga"]),
        "mass_conservation_error": mc_err,
        "clip_triggered": int(counter.clipped_elements),
        # ---- 诊断列（明细表用） ----
        "体积创造峰值": round(peak_volume_ratio(result, well_spec, float(c["Q0"])), 4),
        "eta_E_论文3D": float(c["eta_e_3d"]),
        "eta_E_偏差_vs_D2DGA": round(eta_model - float(c["eta_e_d2dga"]), 4),
        "eta_E_偏差_vs_3D": round(eta_model - float(c["eta_e_3d"]), 4),
        "eta_E_超容差": bool(abs(eta_model - float(c["eta_e_d2dga"])) > eta_tol),
        "t_br_偏差_vs_D2DGA": (
            float("nan") if t_br_model is None else round(t_br_model - float(c["t_br_d2dga"]), 4)
        ),
        "t_br_超容差": (
            bool(t_br_model is None or abs(t_br_model - float(c["t_br_d2dga"])) > t_br_tol)
        ),
        "t_br_原始_s": float("nan") if t_br_s is None else round(t_br_s, 3),
        "t_br_v2口径": float("nan") if t_br_deficit is None else round(t_br_deficit, 4),
        "t_br_v2口径_s": float("nan") if t_br_deficit_s is None else round(t_br_deficit_s, 3),
        "b_反算_论文式": round(case_buoyancy_number(c), 4),
        "b_模型summary": round(float(result.summary["buoyancy_number"]), 4),
        "Q0_用": float(c["Q0"]),
        "Q0_论文打印": float(c["Q0_printed"]),
        "网格_nz": int(nz),
        "网格_ny": int(ny),
        "时间步数": int(len(result.metrics)),
        "clip_调用次数": int(counter.calls),
        "clip_最大绝对值": round(counter.max_abs, 6),
        "clip_理论上界": round(clip_trigger_bound(c), 6),
    }
    if keep_result:
        row["_result"] = result
        row["_solver"] = solver
        row["_well_spec"] = well_spec
    return row


# 对照表列（Task 12 验收口径）：case, e, m, b, η_E 模型/论文, Δη_E,
# t_br 模型/论文, Δt_br, mass_err。更全的诊断字段见 DETAIL_COLUMNS（对照明细.csv）。
COMPARISON_COLUMNS: tuple[str, ...] = (
    "case_id", "e", "m", "b",
    "eta_E_模型", "eta_E_论文D2DGA", "eta_E_偏差_vs_D2DGA",
    "t_br_模型", "t_br_论文D2DGA", "t_br_偏差_vs_D2DGA",
    "mass_conservation_error",
)

DETAIL_COLUMNS: tuple[str, ...] = (
    "case_id", "e", "Re", "m", "b",
    "eta_E_模型", "eta_E_论文D2DGA", "eta_E_论文3D", "eta_E_偏差_vs_D2DGA", "eta_E_偏差_vs_3D",
    "eta_E_超容差",
    "t_br_模型", "t_br_论文D2DGA", "t_br_论文3D", "t_br_偏差_vs_D2DGA", "t_br_超容差",
    "t_br_原始_s", "t_br_v2口径", "t_br_v2口径_s",
    "mass_conservation_error", "体积创造峰值", "clip_triggered",
    "clip_调用次数", "clip_最大绝对值",
    "clip_理论上界", "b_反算_论文式", "b_模型summary", "Q0_用", "Q0_论文打印",
    "网格_nz", "网格_ny", "时间步数",
)


def _write_csv(path: Path, rows: list[Mapping[str, Any]], columns: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in columns})


def run_all_cases(
    *,
    nz: int = 140,
    ny: int = 40,
    output_dir: Path | None = None,
    cases: tuple[int, ...] | None = None,
) -> list[dict[str, Any]]:
    """跑完 10 个基准算例，写出 对照表.csv / 对照明细.csv。"""
    selected = cases if cases is not None else tuple(c["case_id"] for c in ZHANG2022_CASES)
    rows: list[dict[str, Any]] = []
    for cid in selected:
        print(f"[基准算例] case {cid} 开始 ...")
        row = run_case(cid, nz=nz, ny=ny)
        rows.append(row)
        print(
            f"[基准算例] case {cid} 完成: η_E={row['eta_E_模型']:.4f}"
            f"(论文 D2DGA {row['eta_E_论文D2DGA']:.2f}) "
            f"t_br={row['t_br_模型']}(论文 {row['t_br_论文D2DGA']:.2f}) "
            f"守恒误差={row['mass_conservation_error']:.2e} "
            f"体积创造峰值={row['体积创造峰值']:.3f} "
            f"clip={row['clip_triggered']}"
        )

    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        _write_csv(output_dir / "对照表.csv", rows, COMPARISON_COLUMNS)
        _write_csv(output_dir / "对照明细.csv", rows, DETAIL_COLUMNS)
        print(f"[基准算例] 已写出 {output_dir / '对照表.csv'}")
    return rows


# ---------------------------------------------------------------------------
# 网格 / 时间步收敛检验
# ---------------------------------------------------------------------------
def convergence_scan(
    case: int | Mapping[str, Any],
    *,
    nz_values: tuple[int, ...] = (140, 250, 500),
    ny: int = 40,
) -> list[dict[str, float]]:
    """对单个算例做 nz 细化扫描，返回每档 (nz, dt中位, η_E, t_br)。"""
    out: list[dict[str, float]] = []
    for nz in nz_values:
        row = run_case(case, nz=nz, ny=ny, keep_result=True)
        result = row["_result"]
        dt = result.metrics["time_s"].diff().dropna()
        out.append({
            "nz": float(nz),
            "eta_E_模型": float(row["eta_E_模型"]),
            "t_br_模型": float(row["t_br_模型"]),
            "dt_中位_s": float(dt.median()),
            "时间步数": float(row["时间步数"]),
        })
        print(f"[收敛扫描] nz={nz}: η_E={out[-1]['eta_E_模型']:.4f} t_br={out[-1]['t_br_模型']:.4f}")
    return out


def format_report(rows: list[Mapping[str, Any]]) -> str:
    """把对照行渲染成便于阅读的 Markdown 文本。"""
    lines = [
        "## Zhang & Frigaard (2022) 基准算例对照",
        "",
        "| case | e | Re | m | b | t_br 模型 | t_br 论文D2DGA | t_br 论文3-D "
        "| η_E 模型 | η_E 论文D2DGA | 守恒误差 | 体积创造峰值 | clip |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        lines.append(
            f"| {r['case_id']} | {r['e']} | {r['Re']} | {r['m']} | {r['b']} "
            f"| {r['t_br_模型']} | {r['t_br_论文D2DGA']} | {r['t_br_论文3D']} "
            f"| {r['eta_E_模型']} | {r['eta_E_论文D2DGA']} "
            f"| {r['mass_conservation_error']:.2e} | {r['体积创造峰值']:.3f} "
            f"| {r['clip_triggered']} |"
        )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> None:
    """命令行入口：跑 10 个算例并写出对照表。

    用法：``python -m cemdisp.runners.zhang2022_benchmark --out-dir DIR``。
    ⚠️ ``--out-dir`` 必填：缺省时一律报错退出，绝不回落权威目录
    ``results/基准算例对照_2026-09-10``（R18 权威目录不许覆写；终审 I-2 修复，
    2026-09-15）。如确需覆写权威目录，须显式传 ``--force-authoritative-dir``。
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="Zhang & Frigaard (2022, JFM 947:A32) 基准算例对照运行器",
    )
    parser.add_argument(
        "--out-dir", type=Path, default=None,
        help="输出目录（必填。权威目录 results/基准算例对照_2026-09-10 默认拒写，"
             "验收重跑请显式传新目录）",
    )
    parser.add_argument(
        "--force-authoritative-dir", action="store_true",
        help="显式豁免：允许在未传 --out-dir 时写权威目录 "
             "results/基准算例对照_2026-09-10（自担覆写责任）",
    )
    args = parser.parse_args(argv)
    if args.out_dir is not None:
        output_dir = args.out_dir
    elif args.force_authoritative_dir:
        output_dir = DEFAULT_OUTPUT_DIR
    else:
        parser.error(
            "必须显式传 --out-dir 指定输出目录：缺省会回落覆写权威目录 "
            f"{DEFAULT_OUTPUT_DIR}（R18 不许覆写）。如确需写权威目录，"
            "请加 --force-authoritative-dir。"
        )
        return  # 不可达：parser.error 抛 SystemExit(2)
    rows = run_all_cases(output_dir=output_dir)
    print()
    print(format_report(rows))


if __name__ == "__main__":
    main()
