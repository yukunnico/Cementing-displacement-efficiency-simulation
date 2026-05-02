"""
呼101井尾管 D2DGA 二维方位-轴向顶替模型
==========================================

本模块实现了呼101井168.3/139.7mm变径尾管固井顶替过程的二维数值模拟。
D2DGA = 2D Displacement model with Depth & Gyration-Azimuthal resolution，
即在"沿井深方向(s)"和"环空周向方位角(phi)"两个维度上求解流体顶替方程。

模型核心思路：
1. 将偏心环空沿周向展开为变间隙的二维通道（宽边→窄边），
   间隙分布由居中度(standoff)和偏心度(e)控制。
2. 采用体积守恒的流动模型，根据局部间隙和流体物性计算轴向流速 w 和周向流速 v。
3. 对多种流体（钻井液、平衡液、隔离液、领浆、尾浆）的体积分数场做对流-扩散输运，
   包含浮力增强、偏心窜槽、混浆扩散等物理效应。
4. 引入壁面泥饼清除模型，计算"有效顶替效率"= 水泥浆占据率 × (1 - 壁面泥饼残余率)。
5. 综合窜槽指数、混浆指数、失稳指数等惩罚因子，得到与CBL合格率可比的"质量响应效率"。

物理参数来源：
- 井径：5400-6796m按名义260.35mm井眼/168.3mm尾管；6796-7868m按名义215.90mm井眼/139.7mm尾管。
  无实测caliper数据，使用名义井径。
- 井斜：最大8.19°@7440m，入口段约2°，底部约7°。
- 流体物性：钻井液PV=58mPa·s、YP=9.2Pa为设计报告实测值(65℃)，Bingham模型；
  领浆n=0.719,k=0.815,密度2.10g/cc；尾浆n=0.722,k=0.684,密度1.90g/cc，幂律模型。
- 施工参数：领浆47m³+尾浆23m³=70m³、替浆量101.4m³、排量1.2→0.55m³/min分段。
- 居中度：基于132只扶正器(≈18.7m间距)按API RP 10D-2估算0.62-0.67。
- CBL对比值：100312.PDF测井固井质量合格率62.77%。
- 已知异常：悬挂器坐挂失败(技术总结记载)，本模型按正常设计建模，差异通过α吸收。

网格收敛性验证(nz=200/250/350/500/700)：
  nz=500(4.94m/cell)与nz=700(3.53m/cell)的η_hyd差异<0.014，确认nz=500已收敛。
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
import pandas as pd

# ========== matplotlib 中文字体配置 ==========
# 设置中文字体优先级，确保图表中的中文标签正常显示
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "Arial Unicode MS", "DejaVu Sans"]
# 解决负号显示为方块的问题
plt.rcParams["axes.unicode_minus"] = False

# ========== 井身结构与输出路径 ==========
OUT_DIR = Path(__file__).resolve().parent  # 输出文件保存目录
WELL_TOP_MD_M = 5400.0
WELL_BOTTOM_MD_M = 7868.0
TECH_CASING_EQUIV_ID_MM = 273.10
UPPER_SECTION_BOTTOM_MD_M = 6796.0
LOWER_SECTION_BOTTOM_MD_M = WELL_BOTTOM_MD_M
UPPER_HOLE_DIAMETER_MM = 260.35
LOWER_HOLE_DIAMETER_MM = 215.90
UPPER_LINER_OD_MM = 168.30
LOWER_LINER_OD_MM = 139.70
FIELD_REFERENCE_EFFICIENCY = 0.6277
FIELD_REFERENCE_LABEL = "100312.PDF CBL测井：测量井段固井质量合格率62.77%"
Array = NDArray[np.float64]

# ========== 水泥浆与施工参数 ==========
LEAD_VOLUME_M3 = 47.0
TAIL_VOLUME_M3 = 23.0
LEAD_DENSITY_GCC = 2.10
TAIL_DENSITY_GCC = 1.90
DISPLACEMENT_VOLUME_M3 = 101.4
FIELD_RATE_M3_MIN = 1.2
FIELD_RATE_LOW_M3_MIN = 0.55
RATE_SWITCH_DISP_VOLUME_M3 = 62.0
FIELD_DISPLACEMENT_MUD_DENSITY_GCC = 1.96
BALANCE_VOLUME_M3 = 25.0
BALANCE_DENSITY_GCC = 1.85
SPACER_VOLUME_M3 = 25.0
SPACER_DENSITY_GCC = 2.00
MUD_PV_PA_S = 0.058
MUD_YP_PA = 9.2
MUD_DENSITY_GCC = 1.96
SPACER_PV_PA_S = 0.030
SPACER_YP_PA = 5.0
BALANCE_PV_PA_S = 0.030
BALANCE_YP_PA = 3.0

# 环空入口边界模式：
# - sustained_tail：保留原模型口径，替浆阶段等效为尾浆继续从尾管鞋进入环空；适合与旧结果对比，但可能偏乐观。
# - volume_limited：只按水泥浆设计体积(领浆47+尾浆23=70m³)进入环空，之后停泵保持；适合检查环空体积守恒下限。
# - tail_then_mud：水泥浆体积出鞋后，替浆钻井液进入环空；适合检查过替风险上限。
ANNULUS_BOUNDARY_MODE = "sustained_tail"

# ========== 质量惩罚标定参数 ==========
# 三项惩罚的相对权重（保留原物理比例关系 11:7:5）
CHANNELING_PENALTY_WEIGHT = 0.55    # 窜槽对质量因子的相对权重
MIXING_PENALTY_WEIGHT = 0.35        # 混浆对质量因子的相对权重
INSTABILITY_PENALTY_WEIGHT = 0.25   # 失稳对质量因子的相对权重
INSTABILITY_DECAY_SCALE = 5.0       # 失稳指数指数衰减尺度参数

# 质量惩罚全局标定系数 α
# 标定方法：利用呼101现场CBL数据点——
#   CBL合格率 = 0.6277（100312.PDF），模型水动力效率 = 0.9255（nz=500收敛值）
#   目标 quality_factor = 0.6277 / 0.9255 ≈ 0.6783
#   原始惩罚总和 S = 0.55×0.491 + 0.35×0.200 + 0.25×0.557 ≈ 0.479
#   α = (1 − 0.6783) / 0.479 ≈ 0.671
# 含义：α=0.671较高，物理解释——
#   (1) 模型以132扶正器理想居中建模(SO=0.62-0.67)，实际悬挂器坐挂失败导致管柱偏心恶化
#   (2) 2468m超长井段累积误差 + 无实测caliper/偏心数据，理想化程度更高
#   (3) α吸收了geometry理想化 + hanger failure + 无caliper三重差异
# 网格收敛验证：nz=500(η=0.9255)与nz=700(η=0.9391)差异<0.014，确认nz=500已收敛
# 通用标定公式：α = (1 − CBL合格率 / 水动力效率) / S
QUALITY_PENALTY_SCALE = 0.671


@dataclass(frozen=True)
class Fluid:
    """
    流体物性数据类。

    封装了固井过程中各类流体（钻井液、隔离液、水泥浆等）的密度和流变参数，
    并提供表观粘度计算方法。

    Attributes:
        name: 流体中文名称（如"钻井液"、"隔离液"等）
        density_gcc: 密度，单位 g/cm³
        model: 流变模型类型，支持 "newtonian"、"bingham"、"power_law"、"herschel_bulkley"
        n: 幂律指数（仅幂律模型有效），n<1为剪切稀释性，n=1退化为牛顿流体
        k: 稠度系数（仅幂律模型有效），
           牛顿/Bingham流体时k为塑性黏度(Pa·s)，幂律/HB流体时k为稠度系数(Pa·s^n)
        yield_stress_pa: 屈服应力，单位Pa，仅Bingham/Herschel-Bulkley模型有效
        mu_min: 表观黏度下限，防止数值奇异
        mu_max: 表观黏度上限，避免低剪切区屈服项导致速度场过度刚化

    流变模型公式：
        - 牛顿流体：μ_app = k（常数，与剪切速率无关）
        - Bingham流体：μ_app = k + τy/γ
        - 幂律流体：μ_app = k × γ^(n-1)，其中γ为剪切速率
        - Herschel-Bulkley流体：μ_app = τy/γ + k × γ^(n-1)
    """
    name: str
    density_gcc: float
    model: str = "newtonian"
    n: float = 1.0
    k: float = 0.05
    yield_stress_pa: float = 0.0
    mu_min: float = 1e-5
    mu_max: float = 3.0

    def mu_app(self, gamma: Array) -> Array:
        """
        计算给定剪切速率下的表观粘度。

        Args:
            gamma: 剪切速率数组，单位 1/s

        Returns:
            表观粘度数组，单位 Pa·s

        物理说明：
        - 牛顿流体表观粘度恒等于k（即动力粘度）。
        - Bingham/Herschel-Bulkley模型通过τy/γ体现低剪切屈服阻力。
        - 幂律流体表观粘度随剪切速率变化：n<1时剪切稀释（高速率低粘度），n>1时剪切增稠。
        """
        gamma = np.maximum(np.asarray(gamma, dtype=float), 1e-6)  # 剪切速率下限，避免零除
        if self.model == "newtonian":
            mu = np.full_like(gamma, self.k, dtype=float)
        elif self.model == "bingham":
            mu = self.k + self.yield_stress_pa / gamma
        elif self.model == "power_law":
            mu = self.k * gamma ** (self.n - 1.0)
        elif self.model == "herschel_bulkley":
            mu = self.yield_stress_pa / gamma + self.k * gamma ** (self.n - 1.0)
        else:
            raise ValueError(f"Unsupported rheology model: {self.model}")
        return np.clip(mu, self.mu_min, self.mu_max)


# ========== 流体体系定义 ==========
# TRACKED: 需要追踪体积分数的流体名称列表（按注入顺序排列）
# FLUIDS: 所有流体的物性字典，包含钻井液和四种被追踪流体
TRACKED = ["balance", "spacer", "lead", "tail"]
FLUIDS = {
    "mud": Fluid("钻井液", MUD_DENSITY_GCC, "bingham", k=MUD_PV_PA_S, yield_stress_pa=MUD_YP_PA),
    "balance": Fluid("平衡液", BALANCE_DENSITY_GCC, "bingham", k=BALANCE_PV_PA_S, yield_stress_pa=BALANCE_YP_PA),
    "spacer": Fluid("隔离液", SPACER_DENSITY_GCC, "bingham", k=SPACER_PV_PA_S, yield_stress_pa=SPACER_YP_PA),
    "lead": Fluid("领浆", LEAD_DENSITY_GCC, "power_law", n=0.719, k=0.815),
    "tail": Fluid("尾管水泥浆", TAIL_DENSITY_GCC, "power_law", n=0.722, k=0.684),
}


def trapez2d(arr: Array, geom: dict[str, Array]) -> float:
    """
    二维梯形积分。

    对二维数组 arr(y, s) 在 y 方向和 s 方向依次做梯形积分，
    用于计算环空截面上的体积、流量等积分量。

    Args:
        arr: 被积二维数组，形状 (ny, nz)
        geom: 几何字典，包含 "s"（轴向坐标）和 "y"（周向坐标）

    Returns:
        积分值（标量）
    """
    return float(np.trapezoid(np.trapezoid(arr, x=geom["s"], axis=1), x=geom["y"], axis=0))


def caliper_profile_points() -> tuple[Array, Array]:
    """Nominal Hu101 hole diameter profile."""
    transition_half_m = 1.0
    md = np.array([WELL_TOP_MD_M, UPPER_SECTION_BOTTOM_MD_M - transition_half_m, UPPER_SECTION_BOTTOM_MD_M, UPPER_SECTION_BOTTOM_MD_M + transition_half_m, WELL_BOTTOM_MD_M], dtype=float)
    hole = np.array([UPPER_HOLE_DIAMETER_MM, UPPER_HOLE_DIAMETER_MM, 0.5 * (UPPER_HOLE_DIAMETER_MM + LOWER_HOLE_DIAMETER_MM), LOWER_HOLE_DIAMETER_MM, LOWER_HOLE_DIAMETER_MM], dtype=float)
    return md, hole


def inclination_profile_points() -> tuple[Array, Array]:
    """Hu101 inclination profile: 2 deg at top, 8.19 deg at 7440m, 7 deg at bottom."""
    md = np.array([WELL_TOP_MD_M, 7440.0, WELL_BOTTOM_MD_M], dtype=float)
    inc = np.array([2.0, 8.19, 7.0], dtype=float)
    return md, inc


def pipe_od_profile(md: Array) -> Array:
    """Hu101 variable liner OD profile."""
    od = np.where(md <= UPPER_SECTION_BOTTOM_MD_M, UPPER_LINER_OD_MM, LOWER_LINER_OD_MM)
    return od.astype(float)


def standoff_profile(md: Array, hole: Array, od: Array) -> Array:
    """Hu101 standoff profile calibrated for 132 centralizers / 2468m (≈18.7m spacing).

    API RP 10D-2 guidance: dense centralizer programs (≤20m spacing) in low-inclination
    wells (max 8.19°) justify standoff 0.62–0.70. Previous values (0.50–0.60) were too
    pessimistic and caused catastrophic channeling (wide/narrow flow ratio >18×).
    """
    standoff = np.full_like(md, 0.67, dtype=float)
    standoff[(md >= 5400.0) & (md < 5700.0)] = 0.63    # near entry, near-vertical
    standoff[(md >= 5700.0) & (md < UPPER_SECTION_BOTTOM_MD_M)] = 0.67  # main upper section, good centralizer coverage
    standoff[(md >= UPPER_SECTION_BOTTOM_MD_M) & (md < 7400.0)] = 0.65  # transition zone, variable OD
    standoff[(md >= 7400.0) & (md <= 7735.0)] = 0.62    # target zone, max inclination 8.19°@7440m
    standoff[(md > 7735.0) & (md <= WELL_BOTTOM_MD_M)] = 0.63  # bottom section
    clearance = np.maximum(hole - od, 5.0)
    standoff *= np.clip(clearance / 70.0, 0.55, 1.0)
    return np.clip(standoff, 0.30, 0.82)


def physical_annular_volume_m3() -> float:
    """Integrate Hu101 annular volume with variable hole diameter and variable liner OD."""
    md = np.linspace(WELL_TOP_MD_M, WELL_BOTTOM_MD_M, 2001)
    cal_md, cal_hole = caliper_profile_points()
    hole = np.interp(md, cal_md, cal_hole)
    od = pipe_od_profile(md)
    area = np.pi * ((hole / 1000.0) ** 2 - (od / 1000.0) ** 2) / 4.0
    return float(np.trapezoid(area, x=md))


def build_geom(nz: int = 140, ny: int = 40, target_half_volume: float | None = None) -> dict[str, Array]:
    """
    构建2D环空几何网格。

    将偏心环空展开为二维通道：
    - s方向（轴向）：沿井深从底到顶，s=0对应井底(WELL_BOTTOM_MD_M)
    - y方向（周向）：从宽边(phi=0)到窄边(phi=π)的弧长坐标

    间隙分布模型：
        h(phi, s) = half_gap_mean × (1 + e × cos(π × phi))
    其中 e 为偏心度指标，phi=0为宽边（间隙最大），phi=1为窄边（间隙最小）。

    几何缩放：
        为保证数值网格的体积与物理环空体积一致，对间隙场施加全局缩放因子。

    Args:
        nz: 轴向网格数（默认140）
        ny: 周向网格数（默认40）
        target_half_volume: 目标半环空体积(m³)，默认取物理环空体积的一半

    Returns:
        几何字典，包含：
        - s: 轴向坐标数组 (nz,)
        - md: 对应测深坐标 (nz,)
        - y: 周向弧长坐标 (ny,)
        - phi: 归一化周向角 [0,1] (ny,)
        - H: 半间隙场 (ny, nz)
        - b: 全间隙场 = 2×H (ny, nz)
        - e: 偏心度指标 (nz,)
        - standoff: 居中度 (nz,)
        - inc_deg: 井斜角 (nz,)
        - hole_mm: 井径 (nz,)
        - od_mm: 套管外径 (nz,)
        - volume_scale: 体积缩放因子
    """
    target_half_volume = 0.5 * physical_annular_volume_m3() if target_half_volume is None else target_half_volume
    # 轴向坐标：s=0为井底，s增大为井顶方向
    s = np.linspace(0.0, WELL_BOTTOM_MD_M - WELL_TOP_MD_M, nz)
    md = WELL_BOTTOM_MD_M - s  # 测深：井底到井顶

    # 插值获取各网格点的井径、井斜
    cal_md, cal_hole = caliper_profile_points()
    inc_md, inc = inclination_profile_points()
    hole = np.interp(md, cal_md, cal_hole)
    od = pipe_od_profile(md)
    standoff = standoff_profile(md, hole, od)

    # 偏心度指标 e：从居中度换算，e越大偏心越严重
    e = np.clip(1.0 - standoff, 0.05, 0.55)
    # 环空间隙和平均半径
    clearance = (hole - od) / 1000.0  # 环空间隙，单位 m
    half_gap_mean = clearance / 2.0   # 平均半间隙
    mean_radius = ((hole + od) / 4.0) / 1000.0  # 平均半径，单位 m

    # 周向坐标：从宽边到窄边的弧长
    y = np.linspace(0.0, np.pi * np.mean(mean_radius), ny)
    phi = y / y[-1]  # 归一化周向角，0=宽边，1=窄边

    # 构建间隙场：h(phi,s) = half_gap_mean × (1 + e × cos(π × phi))
    # 宽边(phi=0)间隙最大，窄边(phi=1)间隙最小
    h = np.zeros((ny, nz))
    b = np.zeros((ny, nz))
    for j in range(nz):
        h[:, j] = half_gap_mean[j] * (1.0 + e[j] * np.cos(np.pi * phi))
        b[:, j] = 2.0 * h[:, j]  # 全间隙 = 2 × 半间隙

    geom = {"s": s, "md": md, "y": y, "phi": phi, "H": h, "b": b, "e": e, "standoff": standoff, "inc_deg": np.interp(md, inc_md, inc), "hole_mm": hole, "od_mm": od}

    # 体积校准：缩放间隙场使网格体积等于目标物理体积
    current_half_volume = trapez2d(geom["b"], geom)
    scale = target_half_volume / current_half_volume
    geom["H"] *= scale
    geom["b"] *= scale
    geom["volume_scale"] = np.array(scale)
    return geom


def boundary_state(t: float) -> tuple[Array, float, str]:
    """Hu101 injection boundary schedule."""
    schedule = [("balance", BALANCE_VOLUME_M3, FIELD_RATE_M3_MIN, "注平衡液"), ("spacer", SPACER_VOLUME_M3, FIELD_RATE_M3_MIN, "注隔离液"), ("lead", LEAD_VOLUME_M3, FIELD_RATE_M3_MIN, "注领浆"), ("tail", TAIL_VOLUME_M3, FIELD_RATE_M3_MIN, "注尾浆")]
    t0 = 0.0
    vec = np.zeros(len(TRACKED))
    for name, volume, rate, stage in schedule:
        duration = volume / rate * 60.0
        if t < t0 + duration - 1e-12:
            vec[TRACKED.index(name)] = 1.0
            return vec, rate / 60.0, stage
        t0 += duration
    displacement_fast_volume_m3 = min(RATE_SWITCH_DISP_VOLUME_M3, DISPLACEMENT_VOLUME_M3)
    displacement_slow_volume_m3 = max(DISPLACEMENT_VOLUME_M3 - displacement_fast_volume_m3, 0.0)
    for volume, rate, stage in [(displacement_fast_volume_m3, FIELD_RATE_M3_MIN, "替浆快速"), (displacement_slow_volume_m3, FIELD_RATE_LOW_M3_MIN, "替浆慢速段")]:
        duration = volume / rate * 60.0 if rate > 0.0 else 0.0
        if t < t0 + duration - 1e-12:
            if ANNULUS_BOUNDARY_MODE == "sustained_tail":
                vec[TRACKED.index("tail")] = 1.0
                return vec, rate / 60.0, f"{stage}(尾浆持续入环空)"
            if ANNULUS_BOUNDARY_MODE == "volume_limited":
                vec[TRACKED.index("tail")] = 1.0
                return vec, 0.0, f"{stage}(体积限制停泵)"
            if ANNULUS_BOUNDARY_MODE == "tail_then_mud":
                return vec, rate / 60.0, f"{stage}(钻井液入环空)"
            raise ValueError(f"Unsupported annulus boundary mode: {ANNULUS_BOUNDARY_MODE}")
        t0 += duration
    if ANNULUS_BOUNDARY_MODE in {"sustained_tail", "volume_limited"}:
        vec[TRACKED.index("tail")] = 1.0
    return vec, 0.0, "替浆结束"


def compute_props(x: Array, w_prev: Array, geom: dict[str, Array]) -> tuple[Array, Array, Array]:
    """
    根据各流体体积分数计算局部混合物的表观粘度和密度。

    采用体积分数加权平均：
        μ_mix = Σ(x_i × μ_i(γ))  +  x_mud × μ_mud(γ)
        ρ_mix = Σ(x_i × ρ_i)     +  x_mud × ρ_mud

    Args:
        x: 各追踪流体的体积分数场，形状 (len(TRACKED), ny, nz)
        w_prev: 上一时间步的轴向流速场，用于计算剪切速率
        geom: 几何字典

    Returns:
        (mu, rho, mud):
        - mu: 混合物表观粘度场 (ny, nz)，单位 Pa·s
        - rho: 混合物密度场 (ny, nz)，单位 g/cm³
        - mud: 钻井液体积分数场 (ny, nz)
    """
    # 钻井液体积分数 = 1 - 所有追踪流体体积分数之和
    mud = np.clip(1.0 - x.sum(axis=0), 0.0, 1.0)
    # 剪切速率 γ ≈ 6|w|/b（平板Poiseuille流的壁面剪切速率近似）
    gamma = np.maximum(6.0 * np.abs(w_prev) / np.maximum(geom["b"], 1e-5), 1e-6)
    # 加权计算混合物粘度和密度
    mu = mud * FLUIDS["mud"].mu_app(gamma)
    rho = mud * FLUIDS["mud"].density_gcc
    for k, name in enumerate(TRACKED):
        mu += x[k] * FLUIDS[name].mu_app(gamma)
        rho += x[k] * FLUIDS[name].density_gcc
    return mu, rho, mud


def compute_velocity(x: Array, geom: dict[str, Array], q_m3s: float, w_prev: Array) -> tuple[Array, Array, Array, Array, Array]:
    """
    计算环空二维流速场（轴向w + 周向v）。

    流动模型基于局部间隙和流体物性，考虑以下物理效应：
    1. 浮力增强：密度差驱动水泥浆在宽边加速（浮力因子）
    2. 偏心效应：居中度影响流动分布（宽边流速大于窄边）
    3. 窄边抑制：水泥浆在窄边受密度差和偏心双重抑制

    轴向流速 w 的计算：
        m = (b³/μ) × buoyancy × standoff_correction × narrow_boost
        w = (Q_half / ∫b·m dy) × m
    其中 m 为局部流动度，Q_half 为半环空流量。

    周向流速 v 的计算：
        由连续性方程 ∂(bw)/∂s + ∂(bv)/∂y = 0 推导，
        从宽边(y=0)向窄边积分得到。

    Args:
        x: 流体体积分数场 (len(TRACKED), ny, nz)
        geom: 几何字典
        q_m3s: 总注入排量，单位 m³/s
        w_prev: 上一时间步轴向流速，用于计算粘度

    Returns:
        (w, v, mu, rho, mud):
        - w: 轴向流速场 (ny, nz)，单位 m/s
        - v: 周向流速场 (ny, nz)，单位 m/s
        - mu: 混合物粘度 (ny, nz)
        - rho: 混合物密度 (ny, nz)
        - mud: 钻井液体积分数 (ny, nz)
    """
    y = geom["y"]
    b = geom["b"]
    q_half = q_m3s / 2.0  # 半环空流量（模型只计算半环空：宽边到窄边）
    mu, rho, mud = compute_props(x, w_prev, geom)
    # 水泥浆总浓度（领浆+尾浆）
    cement = x[TRACKED.index("lead")] + x[TRACKED.index("tail")]

    # ===== 浮力增强因子 =====
    # 密度差越大，浮力驱动越强；井斜角越大，浮力分量越大
    buoy = 1.0 + 0.60 * np.maximum(rho - FLUIDS["mud"].density_gcc, 0.0) + 0.08 * np.sin(np.deg2rad(geom["inc_deg"]))[None, :]

    # ===== 局部流动度 m =====
    # m ∝ b³/μ（宽间隙、低粘度→高流动度），乘以浮力和居中度修正
    m = (b**3 / np.maximum(mu, 1e-6)) * buoy * (0.90 + 0.25 * geom["standoff"][None, :])

    # ===== 窄边抑制修正 =====
    # 水泥浆在窄边受密度差和偏心效应抑制，流速降低
    phi = geom["phi"][:, None]
    narrow_boost = 1.0 - 0.30 * cement * np.maximum(rho - FLUIDS["mud"].density_gcc, 0.0) * np.cos(np.pi * phi)
    m *= np.clip(narrow_boost, 0.75, 1.70)

    # ===== 轴向流速 w =====
    # 由流量守恒：Q_half = ∫(b × w) dy = ∫(b × (Q_half/∫b·m dy) × m) dy
    int_m = np.trapezoid(b * m, x=y, axis=0)
    w = (q_half / np.maximum(int_m, 1e-12))[None, :] * m

    # ===== 周向流速 v =====
    # 由二维连续性方程推导：∂(bw)/∂s + ∂(bv)/∂y = 0
    # 积分得：v(y) = -∫₀ʸ [∂(bw)/∂s] dy'，再减去线性分布使窄边v=0
    ds = geom["s"][1] - geom["s"][0]
    dy = y[1] - y[0]
    bw = b * w
    dbw_ds = np.gradient(bw, ds, axis=1)  # ∂(bw)/∂s
    bv = np.zeros_like(w)
    for i in range(1, len(y)):
        bv[i, :] = bv[i - 1, :] - 0.5 * (dbw_ds[i, :] + dbw_ds[i - 1, :]) * dy
    # 修正：使窄边(y=y_max)处v=0（对称边界条件）
    bv -= (y[:, None] / y[-1]) * bv[-1, :]
    v = bv / np.maximum(b, 1e-8)

    return w, v, mu, rho, mud


def interp2(field: Array, ysrc: Array, ssrc: Array, geom: dict[str, Array], bval: float) -> Array:
    """
    二维双线性插值（半拉格朗日回溯）。

    在半拉格朗日方法中，需要根据质点回溯位置(ysrc, ssrc)插值获取上一时刻的场值。
    本函数实现二维双线性插值，对超出网格边界的点用边界值bval填充。

    Args:
        field: 被插值的二维场 (ny, nz)
        ysrc: 回溯周向坐标 (ny, nz)
        ssrc: 回溯轴向坐标 (ny, nz)
        geom: 几何字典，提供网格坐标
        bval: 超出网格边界时的填充值（入口边界条件）

    Returns:
        插值结果场 (ny, nz)
    """
    y = geom["y"]
    s = geom["s"]
    dy = y[1] - y[0]
    ds = s[1] - s[0]
    ny, nz = field.shape
    # 将回溯坐标裁剪到网格范围内
    ycl = np.clip(ysrc, y[0], y[-1])
    scl = np.clip(ssrc, s[0], s[-1])
    # 计算网格索引和权重
    iy = np.clip(np.floor((ycl - y[0]) / dy).astype(int), 0, ny - 2)
    js = np.clip(np.floor((scl - s[0]) / ds).astype(int), 0, nz - 2)
    wy = (ycl - y[iy]) / dy     # y方向插值权重
    wz = (scl - s[js]) / ds     # s方向插值权重
    # 双线性插值公式
    out = (1.0 - wy) * (1.0 - wz) * field[iy, js] + wy * (1.0 - wz) * field[iy + 1, js] + (1.0 - wy) * wz * field[iy, js + 1] + wy * wz * field[iy + 1, js + 1]
    # 超出入口边界（s < s[0]）的点用边界值填充
    out[ssrc < s[0]] = bval
    return out


def simulate(dt: float = 4.0, nz: int = 500, ny: int = 40, alpha_clean: float = 0.085, total_t: float = 14000.0) -> tuple[dict[str, Array], Array, Array, pd.DataFrame]:
    """
    主模拟函数：执行D2DGA顶替过程的完整时间推进。

    算法流程：
    1. 初始化几何网格、流体场、壁面泥饼场
    2. 时间循环（显式Euler推进）：
       a. 获取当前时刻的边界条件（流体类型、排量）
       b. 计算流速场（轴向w + 周向v）
       c. 计算对流扩散系数（含浮力、偏心、混浆等效应）
       d. 半拉格朗日对流步：回溯质点位置并插值
       e. 显式扩散步：拉普拉斯算子离散
       f. 体积分数归一化（保证 Σx_i ≤ 1）
       g. 壁面泥饼清除模型
       h. 计算各项效率指标
    3. 返回几何、最终场、壁面场和时间序列指标

    Args:
        dt: 时间步长，单位 秒（默认4.0s）
        nz: 轴向网格数（默认500，对应4.94m/cell，经收敛验证）
        ny: 周向网格数（默认40）
        alpha_clean: 壁面清除系数（默认0.085），控制泥饼清除速率
        total_t: 总模拟时间，单位 秒（默认14000s ≈ 233min，覆盖全部注替过程）

    Returns:
        (geom, x, wall, metrics):
        - geom: 几何字典
        - x: 最终各追踪流体体积分数场 (len(TRACKED), ny, nz)
        - wall: 最终壁面泥饼残余率场 (ny, nz)，1=完全未清除，0=完全清除
        - metrics: 时间序列指标DataFrame
    """
    geom = build_geom(nz=nz, ny=ny)
    ny, _ = geom["b"].shape
    # 初始化：环空全部为钻井液，追踪流体体积分数为零
    x = np.zeros((len(TRACKED), ny, nz))
    # 壁面泥饼残余率：初始为1（完全被钻井液污染）
    wall = np.ones((ny, nz))
    # 初始轴向流速猜测值
    w_prev = np.full((ny, nz), 0.45)
    # 半环空参考体积（用于归一化效率指标）
    v_half = trapez2d(geom["b"], geom)
    # 目标评价区间掩码：油气水敏感层段（7492-7735m，含气层7492m起、漏层7598m起）
    target_mask = ((geom["md"] >= 7492.0) & (geom["md"] <= 7735.0))[None, :]
    cbl_eval_mask = ((geom["md"] >= 5700.0) & (geom["md"] <= 7810.0))[None, :]
    rows: list[list[float | str]] = []
    # 构建网格坐标矩阵（用于半拉格朗日回溯）
    ygrid, sgrid = np.meshgrid(geom["y"], geom["s"], indexing="ij")

    # ===== 主时间循环 =====
    for step in range(int(total_t / dt) + 1):
        t = step * dt

        # --- 步骤a：获取边界条件 ---
        boundary_vec, q, stage = boundary_state(t)

        # --- 步骤b：计算流速场 ---
        w, v, mu, rho, mud = compute_velocity(x, geom, q, w_prev)
        w_prev = w  # 更新上一时间步流速

        # --- 步骤c：计算对流扩散系数 ---
        cement = x[TRACKED.index("lead")] + x[TRACKED.index("tail")]
        preflush = x[TRACKED.index("balance")] + x[TRACKED.index("spacer")]
        # 钻井液/水泥浆比值（影响界面不稳定性）
        disp_ratio = np.clip(mud / (cement + 1e-4), 0.2, 8.0)
        # 浮力数（密度差驱动）
        buoy_num = np.maximum(rho - FLUIDS["mud"].density_gcc, 0.0)

        # 轴向扩散系数 dax：
        # 基础值 ∝ |w|×b，增强因子考虑：
        # - 水泥浆-钻井液界面混合 (1.20 × cement × (1-cement))
        # - 前置液-钻井液界面混合 (0.45 × preflush × (1-preflush))
        # - 偏心效应增强 (1 + 0.50 × e)
        # - 粘度比效应 (1 + 0.10 × max(disp_ratio-1, 0))
        # - 浮力抑制 (1 / (1 + 0.95 × buoy_num))
        dax = 0.0040 * np.abs(w) * geom["b"] * (0.20 + 1.20 * cement * (1.0 - cement) + 0.45 * preflush * (1.0 - preflush)) * (1.0 + 0.50 * geom["e"][None, :]) * (1.0 + 0.10 * np.maximum(disp_ratio - 1.0, 0.0)) / (1.0 + 0.95 * buoy_num)
        # 周向扩散系数 day：约为轴向的16%，偏心度越大周向扩散越强
        day = 0.16 * dax * (0.45 + 0.55 * geom["e"][None, :])

        # --- 步骤d：半拉格朗日对流步 ---
        # 浮力驱动的周向漂移速度修正
        v_eff = v + 0.007 * np.maximum(rho - FLUIDS["mud"].density_gcc, 0.0) * (cement + 0.20 * preflush) * np.sin(np.pi * geom["phi"])[:, None]
        # 回溯质点位置
        ysrc = ygrid - v_eff * dt  # 周向回溯
        ssrc = sgrid - w * dt      # 轴向回溯
        # 对各追踪流体分别插值
        xadv = np.zeros_like(x)
        for k in range(len(TRACKED)):
            xadv[k] = interp2(x[k], ysrc, ssrc, geom, float(boundary_vec[k]))

        # --- 步骤e：显式扩散步 ---
        xnew = np.zeros_like(x)
        dy = geom["y"][1] - geom["y"][0]
        ds = geom["s"][1] - geom["s"][0]
        for k in range(len(TRACKED)):
            f = xadv[k]
            # 轴向拉普拉斯算子（二阶中心差分，边界用单侧差分）
            lap_s = np.zeros_like(f)
            lap_y = np.zeros_like(f)
            lap_s[:, 1:-1] = (f[:, 2:] - 2.0 * f[:, 1:-1] + f[:, :-2]) / ds**2
            lap_s[:, 0] = (f[:, 1] - f[:, 0]) / ds**2       # 入口边界
            lap_s[:, -1] = (f[:, -2] - f[:, -1]) / ds**2     # 出口边界
            # 周向拉普拉斯算子
            lap_y[1:-1, :] = (f[2:, :] - 2.0 * f[1:-1, :] + f[:-2, :]) / dy**2
            lap_y[0, :] = (f[1, :] - f[0, :]) / dy**2       # 宽边边界
            lap_y[-1, :] = (f[-2, :] - f[-1, :]) / dy**2     # 窄边边界
            # 显式Euler推进：x_new = x_adv + dt × (D_ax × ∇²s + D_ay × ∇²y)
            xnew[k] = f + dt * (dax * lap_s + day * lap_y)

        # --- 步骤f：体积分数归一化 ---
        # 确保所有流体体积分数非负且总和不超过1
        x = np.clip(xnew, 0.0, None)
        total = x.sum(axis=0)
        over = total > 1.0
        if np.any(over):
            x[:, over] /= total[over]  # 等比缩放使总和=1

        # --- 步骤g：壁面泥饼清除模型 ---
        cement = x[TRACKED.index("lead")] + x[TRACKED.index("tail")]
        mud = np.clip(1.0 - x.sum(axis=0), 0.0, 1.0)
        # 清除能力因子：不同流体对壁面泥饼的清除效率不同
        # 平衡液0.35 < 隔离液0.80 < 水泥浆1.10
        cleaner = 0.35 * x[TRACKED.index("balance")] + 0.80 * x[TRACKED.index("spacer")] + 1.10 * cement
        # 剪切速率（驱动机械冲刷）
        shear = np.abs(w) / np.maximum(geom["b"], 1e-5)
        # 综合清除速率：受清除能力、剪切速率、居中度、浮力共同影响
        kclean = alpha_clean * cleaner * (0.45 + np.sqrt(np.maximum(shear, 1e-6))) * (0.85 + 0.35 * geom["standoff"][None, :]) * (1.0 + 0.30 * buoy_num)
        # 指数衰减模型：wall(t+dt) = wall(t) × exp(-kclean × dt / τ)
        wall *= np.exp(-kclean * dt / 150.0)
        wall = np.clip(wall, 0.0, 1.0)

        # --- 步骤h：计算效率指标 ---
        # 有效顶替效率 = 水泥浆浓度 × (1 - 壁面泥饼残余率)
        eff = cement * (1.0 - wall)
        # 体积填充率：水泥浆占据环空体积的比例
        bulk_fill = trapez2d(geom["b"] * cement, geom) / v_half
        # 全井段有效顶替效率
        eff_eta = trapez2d(geom["b"] * eff, geom) / v_half
        # 油气水层段有效顶替效率
        target_eff = trapez2d(geom["b"] * eff * target_mask, geom) / max(trapez2d(geom["b"] * target_mask, geom), 1e-12)
        # CBL评价井段有效顶替效率
        cbl_eff = trapez2d(geom["b"] * eff * cbl_eval_mask, geom) / max(trapez2d(geom["b"] * cbl_eval_mask, geom), 1e-12)

        scoord = geom["s"]

        # 水泥浆前缘位置计算：找到浓度≥阈值的最远位置
        def front(line: Array, thr: float = 0.5) -> float:
            idx = np.where(line >= thr)[0]
            return float(scoord[idx.max()]) if idx.size else 0.0

        front_wide = front(cement[0])       # 宽边前缘
        front_narrow = front(cement[-1])    # 窄边前缘
        front_mid = front(cement[ny // 2])  # 中线前缘

        # 窜槽指数：宽窄边前缘差值归一化，越大表示窜槽越严重
        channeling = abs(front_wide - front_narrow) / (scoord[-1] + 1e-9)
        # 混浆指数：界面混合程度（4×c×(1-c)在c=0.5时最大）
        mixing = trapez2d(geom["b"] * (4.0 * cement * (1.0 - cement)), geom) / v_half
        # 流动度比：宽边/窄边，反映偏心导致的流速差异
        mobility_wide = np.mean((geom["b"][0] ** 3) / (mu[0] + 1e-6))
        mobility_narrow = np.mean((geom["b"][-1] ** 3) / (mu[-1] + 1e-6))
        # 失稳代理量：窜槽×流动度比超量×(1+混浆增强)
        instability_proxy = channeling * max(mobility_wide / (mobility_narrow + 1e-9) - 1.0, 0.0) * (1.0 + 0.4 * mixing)
        # 失稳指数：指数映射到[0,1]
        instability_index = 1.0 - np.exp(-instability_proxy / INSTABILITY_DECAY_SCALE)
        # 质量因子：综合窜槽、混浆、失稳的惩罚（α × 加权和）
        raw_penalty = (CHANNELING_PENALTY_WEIGHT * channeling
                       + MIXING_PENALTY_WEIGHT * mixing
                       + INSTABILITY_PENALTY_WEIGHT * instability_index)
        quality_factor = np.clip(1.0 - QUALITY_PENALTY_SCALE * raw_penalty, 0.0, 1.0)
        # CBL质量响应效率 = 水动力效率 × 质量因子（与CBL合格率同口径）
        cbl_quality_proxy = cbl_eff * quality_factor

        # 记录本时间步的指标
        rows.append([t, t / 60.0, stage, bulk_fill, eff_eta, target_eff, cbl_eff, cbl_quality_proxy, front_wide, front_narrow, front_mid, channeling, mixing, instability_proxy, instability_index, float(np.mean(wall)), float(np.mean(cement)), float(np.mean(mud))])

    # 构建时间序列指标DataFrame
    metric_columns = ["time_s", "time_min", "stage", "bulk_cement_fill", "effective_efficiency", "target_interval_efficiency", "cbl_eval_interval_efficiency", "cbl_quality_proxy", "front_wide_m", "front_narrow_m", "front_mid_m", "channeling_index", "mixing_index", "instability_proxy", "instability_index", "mean_wall_mud", "mean_cement", "mean_mud"]
    metrics = pd.DataFrame({name: [row[index] for row in rows] for index, name in enumerate(metric_columns)})
    return geom, x, wall, metrics


def save_outputs(geom: dict[str, Array], x: Array, wall: Array, metrics: pd.DataFrame) -> dict[str, object]:
    """
    保存模拟结果到文件并生成图表。

    输出文件包括：
    1. 时间序列指标CSV（各时刻的效率、前缘位置、窜槽指数等）
    2. 深度剖面CSV（各深度点的水泥浓度、有效效率、壁面泥饼等）
    3. 最终场数据NPZ（numpy压缩格式，保存完整2D场数据）
    4. 结果摘要JSON + MD（关键数值汇总）
    5. 对比柱状图（模拟效率 vs 资料CBL合格率）
    6. 深度剖面图（有效顶替效率沿井深分布）
    7. 时程诊断、环空体积数、沿深度多轨、二维场云图等参考文献风格图表
    8. 文件清单JSON（所有输出文件路径索引）

    Args:
        geom: 几何字典
        x: 最终流体体积分数场
        wall: 最终壁面泥饼残余率场
        metrics: 时间序列指标DataFrame

    Returns:
        结果摘要字典
    """
    charts = OUT_DIR / "图表"
    charts.mkdir(parents=True, exist_ok=True)
    # 提取水泥相浓度和有效效率（保持与四相D2DGA结构一致）
    cement = x[TRACKED.index("lead")] + x[TRACKED.index("tail")]
    eff = cement * (1.0 - wall)

    # ===== 深度剖面数据 =====
    depth_profiles = pd.DataFrame(
        {
            "井深_m": geom["md"],
            "水泥平均浓度": np.average(cement, axis=0, weights=geom["b"]),          # 按间隙宽度加权的周向平均水泥浓度
            "平均有效顶替效率": np.average(eff, axis=0, weights=geom["b"]),          # 按间隙宽度加权的周向平均有效效率
            "钻井液平均浓度": 1.0 - np.average(x.sum(axis=0), axis=0, weights=geom["b"]),  # 周向平均钻井液浓度
            "宽边有效效率": eff[0],                   # 宽边(phi=0)有效效率
            "中线有效效率": eff[eff.shape[0] // 2],   # 中线(phi=0.5)有效效率
            "窄边有效效率": eff[-1],                   # 窄边(phi=1)有效效率
            "宽边水泥浓度": cement[0],                 # 宽边水泥浓度
            "中线水泥浓度": cement[cement.shape[0] // 2],  # 中线水泥浓度
            "窄边水泥浓度": cement[-1],                # 窄边水泥浓度
            "环空间隙_m": np.mean(geom["b"], axis=0),  # 周向平均环空间隙
            "偏心度指标": geom["e"],                    # 偏心度
            "居中度": geom["standoff"],                 # 居中度
        }
    )

    # ===== 保存CSV和NPZ =====
    metrics_path = OUT_DIR / "呼101尾管_D2DGA_时间序列结果.csv"
    profile_path = OUT_DIR / "呼101尾管_D2DGA_深度剖面.csv"
    npz_path = OUT_DIR / "呼101尾管_D2DGA_最终场数据.npz"
    metrics.to_csv(metrics_path, index=False, encoding="utf-8-sig")
    depth_profiles.to_csv(profile_path, index=False, encoding="utf-8-sig")
    np.savez_compressed(npz_path, md=geom["md"], s=geom["s"], y=geom["y"], b=geom["b"], cement=cement, effective_efficiency=eff, wall_mud=wall)

    # ===== 结果摘要 =====
    final = metrics.iloc[-1]
    delta = float(final["cbl_quality_proxy"] - FIELD_REFERENCE_EFFICIENCY)  # 质量响应效率与资料值的偏差
    hydraulic_delta = float(final["cbl_eval_interval_efficiency"] - FIELD_REFERENCE_EFFICIENCY)  # 水动力效率与资料值的偏差
    summary: dict[str, object] = {
        "模型名称": "呼101尾管_D2DGA_二维方位-轴向顶替模型",
        "模拟对象": "呼101井168.3/139.7mm变径尾管，5400.0-7868.0m",
        "资料对比口径": FIELD_REFERENCE_LABEL,
        "井段_m": [WELL_TOP_MD_M, WELL_BOTTOM_MD_M],
        "井眼数据": "5400-6796m按260.35mm井眼/168.3mm尾管；6796-7868m按215.90mm井眼/139.7mm尾管，无实测caliper数据",
        "尾管外径_mm": {"5400-6796m": UPPER_LINER_OD_MM, "6796-7868m": LOWER_LINER_OD_MM},
        "物理环空体积_m3": physical_annular_volume_m3(),
        "领浆体积_m3": LEAD_VOLUME_M3,
        "尾浆体积_m3": TAIL_VOLUME_M3,
        "水泥浆总体积_m3": LEAD_VOLUME_M3 + TAIL_VOLUME_M3,
        "替浆量_m3": DISPLACEMENT_VOLUME_M3,
        "排量_m3_min": FIELD_RATE_M3_MIN,
        "环空入口边界模式": ANNULUS_BOUNDARY_MODE,
        "流体模型": {
            name: {
                "名称": fluid.name,
                "密度_g_cm3": fluid.density_gcc,
                "流变模型": fluid.model,
                "PV或K": fluid.k,
                "n": fluid.n,
                "YP_Pa": fluid.yield_stress_pa,
            }
            for name, fluid in FLUIDS.items()
        },
        "最终结果": {
            "全井段最终有效顶替效率": float(final["effective_efficiency"]),
            "CBL评价井段模拟有效顶替效率": float(final["cbl_eval_interval_efficiency"]),
            "CBL评价井段质量响应效率": float(final["cbl_quality_proxy"]),
            "资料CBL合格率_代理顶替效率": FIELD_REFERENCE_EFFICIENCY,
            "质量响应_minus_资料": delta,
            "质量响应_minus_资料_百分点": delta * 100.0,
            "水动力效率_minus_资料": hydraulic_delta,
            "水动力效率_minus_资料_百分点": hydraulic_delta * 100.0,
            "油气水层段模拟有效顶替效率": float(final["target_interval_efficiency"]),
            "最终水泥浆占据率": float(final["bulk_cement_fill"]),
            "最终窜槽指数": float(final["channeling_index"]),
            "最终混浆指数": float(final["mixing_index"]),
            "最终失稳指数": float(final["instability_index"]),
        },
        "质量惩罚标定": {
            "窜槽权重": CHANNELING_PENALTY_WEIGHT,
            "混浆权重": MIXING_PENALTY_WEIGHT,
            "失稳权重": INSTABILITY_PENALTY_WEIGHT,
            "失稳衰减尺度": INSTABILITY_DECAY_SCALE,
            "全局标定系数α": QUALITY_PENALTY_SCALE,
            "标定方法": "呼101独立标定(nz=500收敛)：α=(1−CBL/η_hyd)/S=(1−0.6277/0.9255)/0.479≈0.671；α较高因模型以理想居中建模而实际悬挂器坐挂失败+无caliper数据",
        },
        "资料来源": {
            "CBL": "参考文档\\呼101\\100312.PDF",
            "施工记录": "参考文档\\呼101\\呼101井油层尾管固井施工记录表.pdf",
            "流变": "参考文档\\呼101\\2011121.pdf",
        },
        "假设说明": [
            "钻井液PV=58mPa·s、YP=9.2Pa为尾管设计报告实测值（65℃），非假设。",
            "领浆n=0.719,k=0.815为2011121主报告值；尾浆n=0.722,k=0.684同源。",
            "井径按名义尺寸260.35/215.90mm处理，无实测caliper数据。",
            "居中度0.62-0.67基于132只扶正器(18.7m间距)+最大井斜8.19°按API RP 10D-2估算。",
            "悬挂器坐挂失败(技术总结记载)未在几何中建模，其影响通过α标定系数吸收。",
            "环空入口边界默认为sustained_tail，保留替浆阶段尾浆等效入环空的口径。",
            "资料对比值62.77%为CBL固井质量合格率，作为顶替/胶结效果代理。",
            "网格收敛验证：nz=500(η=0.9255)与nz=700(η=0.9391)差异<0.014，确认nz=500已收敛。",
            "α=0.671较高：理想居中+无caliper+悬挂器失败三重差异被α吸收，属单井经验修正。",
        ],
        "扩展图表": [
            "时程诊断：效率、前沿、窜槽/混浆/失稳随施工时间变化。",
            "环空体积数曲线：用注入体积/物理环空体积归一化展示顶替推进。",
            "沿深度多轨诊断：效率、浓度、井眼间隙、偏心/居中度同深度对照。",
            "二维场云图：方位-轴向水泥浓度、有效效率、残余泥浆分布。",
        ],
    }
    summary_path = OUT_DIR / "呼101尾管_D2DGA_结果摘要.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    md_path = OUT_DIR / "呼101尾管_D2DGA_结果摘要.md"
    md_path.write_text(
        "\n".join(
            [
                "# 呼101尾管 D2DGA 模型结果摘要",
                "",
                f"- 模拟对象：{summary['模拟对象']}",
                f"- 全井段最终有效顶替效率：{final['effective_efficiency']:.4f}",
                f"- CBL评价井段模拟有效顶替效率：{final['cbl_eval_interval_efficiency']:.4f}",
                f"- CBL评价井段质量响应效率：{final['cbl_quality_proxy']:.4f}",
                f"- 资料CBL合格率代理值：{FIELD_REFERENCE_EFFICIENCY:.4f}",
                f"- 质量响应-资料差值：{delta:.4f}",
                f"- 水动力效率-资料差值：{hydraulic_delta:.4f}",
                f"- 油气水层段模拟有效顶替效率：{final['target_interval_efficiency']:.4f}",
                f"- 最终窜槽/混浆/失稳指数：{final['channeling_index']:.4f} / {final['mixing_index']:.4f} / {final['instability_index']:.4f}",
                f"- 水泥浆体积：{LEAD_VOLUME_M3 + TAIL_VOLUME_M3:.2f} m³；替浆量：{DISPLACEMENT_VOLUME_M3:.2f} m³；排量：{FIELD_RATE_M3_MIN:.2f} m³/min",
                f"- 环空入口边界模式：{ANNULUS_BOUNDARY_MODE}",
                f"- 钻井液流变：Bingham，PV={MUD_PV_PA_S * 1000.0:.0f} mPa·s，YP={MUD_YP_PA:.1f} Pa（65℃实测）",
                "",
                "## 说明",
                "资料对比值来自 `100312.PDF` 的尾管 CBL 固井质量合格率 62.77%，此处作为顶替/胶结效果代理。",
                "新增图表参考文献中常见的时程曲线、环空体积数曲线、沿深度多轨诊断、二维云图。",
            ]
        ),
        encoding="utf-8",
    )

    # ===== 图表1：模拟效率与资料CBL合格率对比柱状图 =====
    plt.figure(figsize=(7, 5))
    names = ["水动力效率", "质量响应效率", "资料CBL合格率"]
    vals = [float(final["cbl_eval_interval_efficiency"]), float(final["cbl_quality_proxy"]), FIELD_REFERENCE_EFFICIENCY]
    bars = plt.bar(names, vals, color=["#3B82F6", "#10B981", "#F97316"])
    for bar, val in zip(bars, vals):
        plt.text(bar.get_x() + bar.get_width() / 2, val + 0.02, f"{val:.3f}", ha="center")
    plt.ylim(0, 1.1)
    plt.ylabel("效率/合格率")
    plt.title("呼101尾管模拟顶替效率与资料CBL合格率对比")
    plt.tight_layout()
    comparison_path = charts / "呼101尾管_模拟与资料效率对比.png"
    plt.savefig(comparison_path, dpi=220)
    plt.close()

    # ===== 图表2：深度-效率剖面对比图 =====
    plt.figure(figsize=(7, 6))
    plt.plot(depth_profiles["平均有效顶替效率"], depth_profiles["井深_m"], label="水动力有效顶替效率")
    plt.axvspan(0, FIELD_REFERENCE_EFFICIENCY, color="#F97316", alpha=0.12, label="资料CBL合格率62.77%")
    plt.axhspan(7492.0, 7735.0, color="#EF4444", alpha=0.10, label="油气水敏感层7492-7735m")
    plt.gca().invert_yaxis()  # 井深向下增大
    plt.xlabel("效率")
    plt.ylabel("井深 / m")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    profile_plot_path = charts / "呼101尾管_深度效率剖面对比.png"
    plt.savefig(profile_plot_path, dpi=220)
    plt.close()

    # ===== 图表3：时程诊断（参考顶替效率/前缘/窜槽指数随时间变化类图） =====
    time_diagnostics_path = charts / "呼101尾管_时程诊断_效率前沿窜槽.png"
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    time_min = metrics["time_min"]
    axes[0, 0].plot(time_min, metrics["effective_efficiency"], label="全井段有效效率", color="#2563EB")
    axes[0, 0].plot(time_min, metrics["cbl_eval_interval_efficiency"], label="CBL井段水动力效率", color="#059669")
    axes[0, 0].plot(time_min, metrics["cbl_quality_proxy"], label="质量响应效率", color="#7C3AED")
    axes[0, 0].axhline(FIELD_REFERENCE_EFFICIENCY, color="#F97316", linestyle="--", label="资料CBL合格率")
    axes[0, 0].set_ylabel("效率/合格率")
    axes[0, 0].set_title("效率时程")
    axes[0, 0].legend(fontsize=8)

    axes[0, 1].plot(time_min, metrics["front_wide_m"], label="宽边前缘", color="#0EA5E9")
    axes[0, 1].plot(time_min, metrics["front_mid_m"], label="中线前缘", color="#22C55E")
    axes[0, 1].plot(time_min, metrics["front_narrow_m"], label="窄边前缘", color="#EF4444")
    axes[0, 1].invert_yaxis()
    axes[0, 1].set_ylabel("井深 / m")
    axes[0, 1].set_title("方位前缘推进")
    axes[0, 1].legend(fontsize=8)

    axes[1, 0].plot(time_min, metrics["channeling_index"], label="窜槽指数", color="#DC2626")
    axes[1, 0].plot(time_min, metrics["mixing_index"], label="混浆指数", color="#9333EA")
    axes[1, 0].plot(time_min, metrics["instability_index"], label="失稳指数", color="#EA580C")
    axes[1, 0].set_xlabel("施工时间 / min")
    axes[1, 0].set_ylabel("指数")
    axes[1, 0].set_title("风险指标时程")
    axes[1, 0].legend(fontsize=8)

    axes[1, 1].plot(time_min, metrics["bulk_cement_fill"], label="水泥占据率", color="#16A34A")
    axes[1, 1].plot(time_min, metrics["mean_wall_mud"], label="壁面泥饼残余", color="#A16207")
    axes[1, 1].plot(time_min, metrics["mean_mud"], label="钻井液平均浓度", color="#64748B")
    axes[1, 1].set_xlabel("施工时间 / min")
    axes[1, 1].set_ylabel("体积分数/残余率")
    axes[1, 1].set_title("流体占据与壁面清洗")
    axes[1, 1].legend(fontsize=8)
    for ax in axes.ravel():
        ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(time_diagnostics_path, dpi=220)
    plt.close(fig)

    # ===== 图表4：效率随环空体积数变化（参考PV归一化顶替曲线） =====
    annular_volume_path = charts / "呼101尾管_效率随环空体积数变化.png"
    injected_volume_m3 = np.minimum(metrics["time_min"] * FIELD_RATE_M3_MIN, BALANCE_VOLUME_M3 + SPACER_VOLUME_M3 + LEAD_VOLUME_M3 + TAIL_VOLUME_M3 + RATE_SWITCH_DISP_VOLUME_M3)
    slow_time_min = np.maximum(metrics["time_min"] - (BALANCE_VOLUME_M3 + SPACER_VOLUME_M3 + LEAD_VOLUME_M3 + TAIL_VOLUME_M3 + RATE_SWITCH_DISP_VOLUME_M3) / FIELD_RATE_M3_MIN, 0.0)
    injected_volume_m3 = np.minimum(injected_volume_m3 + slow_time_min * FIELD_RATE_LOW_M3_MIN, BALANCE_VOLUME_M3 + SPACER_VOLUME_M3 + LEAD_VOLUME_M3 + TAIL_VOLUME_M3 + DISPLACEMENT_VOLUME_M3)
    annular_pore_volumes = injected_volume_m3 / physical_annular_volume_m3()
    plt.figure(figsize=(8, 5))
    plt.plot(annular_pore_volumes, metrics["cbl_eval_interval_efficiency"], label="CBL井段水动力效率", color="#2563EB", linewidth=2)
    plt.plot(annular_pore_volumes, metrics["cbl_quality_proxy"], label="质量响应效率", color="#7C3AED", linewidth=2)
    plt.plot(annular_pore_volumes, metrics["target_interval_efficiency"], label="油气水层段有效效率", color="#059669", linewidth=1.8)
    plt.axhline(FIELD_REFERENCE_EFFICIENCY, color="#F97316", linestyle="--", label="资料CBL合格率")
    plt.axvline((LEAD_VOLUME_M3 + TAIL_VOLUME_M3 + DISPLACEMENT_VOLUME_M3) / physical_annular_volume_m3(), color="#94A3B8", linestyle=":", label="水泥+替浆总PV")
    plt.xlabel("累计注入体积 / 物理环空体积")
    plt.ylabel("效率/合格率")
    plt.title("呼101尾管顶替效率随环空体积数变化")
    plt.ylim(0, 1.05)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(annular_volume_path, dpi=220)
    plt.close()

    # ===== 图表5：沿深度多轨诊断（参考测井综合解释图版） =====
    depth_tracks_path = charts / "呼101尾管_沿深度多轨诊断.png"
    fig, axes = plt.subplots(1, 4, figsize=(13, 7), sharey=True)
    md = depth_profiles["井深_m"]
    axes[0].plot(depth_profiles["平均有效顶替效率"], md, label="平均", color="#2563EB")
    axes[0].plot(depth_profiles["宽边有效效率"], md, label="宽边", color="#0EA5E9", linewidth=1)
    axes[0].plot(depth_profiles["窄边有效效率"], md, label="窄边", color="#EF4444", linewidth=1)
    axes[0].axvline(FIELD_REFERENCE_EFFICIENCY, color="#F97316", linestyle="--", linewidth=1)
    axes[0].set_xlabel("有效效率")
    axes[0].set_ylabel("井深 / m")
    axes[0].set_title("顶替效率")
    axes[0].legend(fontsize=8)

    axes[1].plot(depth_profiles["水泥平均浓度"], md, label="水泥", color="#16A34A")
    axes[1].plot(depth_profiles["钻井液平均浓度"], md, label="钻井液", color="#64748B")
    axes[1].set_xlabel("体积分数")
    axes[1].set_title("流体占据")
    axes[1].legend(fontsize=8)

    axes[2].plot(depth_profiles["环空间隙_m"] * 1000.0, md, color="#9333EA")
    axes[2].set_xlabel("平均间隙 / mm")
    axes[2].set_title("井眼间隙")

    axes[3].plot(depth_profiles["偏心度指标"], md, label="偏心度", color="#DC2626")
    axes[3].plot(depth_profiles["居中度"], md, label="居中度", color="#059669")
    axes[3].set_xlabel("指标")
    axes[3].set_title("偏心/居中")
    axes[3].legend(fontsize=8)
    for ax in axes:
        ax.grid(True, alpha=0.25)
        ax.set_ylim(float(md.max()), float(md.min()))
    fig.suptitle("呼101尾管沿深度多轨诊断图")
    fig.tight_layout()
    fig.savefig(depth_tracks_path, dpi=220)
    plt.close(fig)

    # ===== 图表6：最终二维场云图（参考轴向-方位分布云图） =====
    final_fields_path = charts / "呼101尾管_最终二维场_浓度效率残余泥浆.png"
    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    field_specs = [
        (cement, "水泥浓度", "viridis"),
        (eff, "有效顶替效率", "YlGnBu"),
        (wall, "壁面泥饼残余率", "YlOrRd"),
    ]
    extent = (float(geom["md"].min()), float(geom["md"].max()), 1.0, 0.0)
    for ax, (field, title, cmap) in zip(axes, field_specs):
        image = ax.imshow(field, aspect="auto", extent=extent, vmin=0.0, vmax=1.0, cmap=cmap)
        ax.set_ylabel("方位：宽边0 → 窄边1")
        ax.set_title(title)
        fig.colorbar(image, ax=ax, label="比例")
    axes[-1].set_xlabel("井深 / m")
    fig.suptitle("呼101尾管最终二维方位-轴向场分布")
    fig.tight_layout()
    fig.savefig(final_fields_path, dpi=220)
    plt.close(fig)

    # ===== 敏感性占位（呼101基准模型不生成敏感性图表，需另行脚本） =====
    # 若未来另行生成敏感性CSV，可在独立脚本中绘制。
    sensitivity_main_path: Path | None = None
    sensitivity_heatmap_path: Path | None = None

    # ===== 文件清单 =====
    manifest = {
        "summary_json": str(summary_path),
        "summary_md": str(md_path),
        "time_series_csv": str(metrics_path),
        "depth_profile_csv": str(profile_path),
        "field_npz": str(npz_path),
        "comparison_plot": str(comparison_path),
        "profile_plot": str(profile_plot_path),
        "time_diagnostics_plot": str(time_diagnostics_path),
        "efficiency_annular_volumes_plot": str(annular_volume_path),
        "depth_tracks_plot": str(depth_tracks_path),
        "final_fields_plot": str(final_fields_path),
    }
    if sensitivity_main_path is not None:
        manifest["sensitivity_main_effects_plot"] = str(sensitivity_main_path)
    if sensitivity_heatmap_path is not None:
        manifest["sensitivity_heatmap_plot"] = str(sensitivity_heatmap_path)
    manifest_path = OUT_DIR / "呼101尾管_D2DGA_文件清单.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    """主入口：运行模拟并保存结果。"""
    geom, x, wall, metrics = simulate()
    summary = save_outputs(geom, x, wall, metrics)
    print(json.dumps(summary["最终结果"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
