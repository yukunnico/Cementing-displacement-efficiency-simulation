"""
多井通用环空二维D2DGA求解器首版实现

本模块实现偏心环空二维顶替效率计算的D2DGA（多液双区面积）模型。

核心物理过程：
1. 平流输运：水泥浆与前置/隔离液在偏心环空中随平均流场运移
2. 弥散效应：横向和垂向扩散（与偏心度、流体性质相关）
3. 壁面清除：泥饼在水泥浆和前置/隔离液剪切清除作用下逐渐移除
4. 浮力效应：密度差导致的窄边优先顶替和失稳风险

网格系统：
- s方向：井深方向（从鞋口到悬挂器）
- y方向：方位角方向（0=窄边，π=宽边）
- φ方向：归一化方位角（0~1对应窄边到宽边）

主要类：
- AnnulusD2DGASolver: 环空二维顶替求解器
    * __init__(): 可配置时间步长、网格密度、风险权重等参数
    * run(): 执行完整的环空二维模拟

- AnnulusSimulationResult: 求解结果容器
    * geom: 几何参数字典（网格坐标、偏心度等）
    * cement_field: 水泥浓度场数组
    * wall_field: 壁面泥饼清除场
    * metrics: 时间序列DataFrame（效率、风险指数等）
    * depth_profiles: 深度方向平均剖面

求解参数说明：
- dt: 时间步长（秒），默认4秒
- nz: 井深方向网格数，默认140
- ny: 方位角方向网格数，默认40
- alpha_clean: 泥饼清除系数，默认0.085
- total_t: 总模拟时间（秒），默认12000秒（200分钟）
- enable_d2dga: 是否启用D2DGA通量修正（Zhang & Frigaard 2022），默认开启
- d2dga_viscosity_ratio: D2DGA粘度比 m = η_displaced/η_displacing，默认1.0
- quality_penalty_scale: 质量惩罚因子，默认0.099
- channeling_penalty_weight: 窜槽风险权重，默认0.55
- mixing_penalty_weight: 混浆风险权重，默认0.35
- instability_penalty_weight: 失稳风险权重，默认0.25
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, Tuple

import numpy as np
from numpy.typing import NDArray
import pandas as pd

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.well_spec import DepthValuePoint, WellSpec
from cemdisp.models2d.boundary_bridge import AnnulusInletState


Array = NDArray[np.float64]
BoolArray = NDArray[np.bool_]


def _profile_to_arrays(points: Tuple[DepthValuePoint, ...]) -> Tuple[Array, Array]:
    """将剖面数据点列表转换为NumPy数组。"""
    depths = np.array([float(point.depth_md_m) for point in points], dtype=float)
    values = np.array([float(point.value) for point in points], dtype=float)
    return depths, values


def _window_mask(well_spec: WellSpec, md: Array, window_type: str) -> BoolArray:
    """根据窗口类型和井段规格生成井深掩码数组。"""
    mask = np.zeros_like(md, dtype=bool)
    for window in well_spec.evaluation_windows:
        if window.window_type == window_type:
            mask |= (md >= window.top_md_m) & (md <= window.bottom_md_m)
    return mask[None, :]


def _phase_fraction(inlet_state: AnnulusInletState, phase_name: str) -> float:
    """从环空入口状态中提取指定相的体积分数。"""
    return float(sum(fraction for name, fraction in inlet_state.phase_fractions if name == phase_name))


def _trapez2d(arr: Array, geom: Dict[str, Array]) -> float:
    """使用梯形法则计算二维数组在网格上的积分。"""
    return float(np.trapezoid(np.trapezoid(arr, x=geom["s"], axis=1), x=geom["y"], axis=0))


def _bilinear_interp(field: Array, ysrc: Array, ssrc: Array, geom: Dict[str, Array], inlet_value: float) -> Array:
    """双线性插值，用于平流输运的反演追踪。"""
    y = geom["y"]
    s = geom["s"]
    dy = y[1] - y[0]
    ds = s[1] - s[0]
    ny, nz = field.shape
    ycl = np.clip(ysrc, y[0], y[-1])
    scl = np.clip(ssrc, s[0], s[-1])
    iy = np.clip(np.floor((ycl - y[0]) / dy).astype(int), 0, ny - 2)
    js = np.clip(np.floor((scl - s[0]) / ds).astype(int), 0, nz - 2)
    wy = (ycl - y[iy]) / dy
    wz = (scl - s[js]) / ds
    out = (
        (1.0 - wy) * (1.0 - wz) * field[iy, js]
        + wy * (1.0 - wz) * field[iy + 1, js]
        + (1.0 - wy) * wz * field[iy, js + 1]
        + wy * wz * field[iy + 1, js + 1]
    )
    out[ssrc < s[0]] = inlet_value
    return out


@dataclass(frozen=True)
class AnnulusSimulationResult:
    """环空二维求解结果。

    包含完整的环空二维模拟输出，包括：
    - geom: 几何参数（网格坐标、偏心度、井径等）
    - cement_field: 水泥浓度场（ny×nz数组）
    - spacer_field: 前置液/隔离液浓度场（ny×nz数组）
    - wall_field: 壁面泥饼清除场
    - metrics: 时间序列指标DataFrame
    - depth_profiles: 深度方向平均剖面DataFrame
    - summary: 最终结果摘要字典
    - time_points_s: 时间点序列
    """

    well_name: str
    geom: Dict[str, Array]
    cement_field: Array
    spacer_field: Array
    wall_field: Array
    metrics: pd.DataFrame
    depth_profiles: pd.DataFrame
    summary: Dict[str, object]
    time_points_s: Tuple[float, ...] = field(default_factory=tuple)
    cement_snapshots: Tuple[Array, ...] = field(default_factory=tuple)
    lead_snapshots: Tuple[Array, ...] = field(default_factory=tuple)
    tail_snapshots: Tuple[Array, ...] = field(default_factory=tuple)
    spacer_snapshots: Tuple[Array, ...] = field(default_factory=tuple)
    wall_snapshots: Tuple[Array, ...] = field(default_factory=tuple)
    gel_strength_snapshots: Tuple[Array, ...] = field(default_factory=tuple)
    mud_cake_field: Array = field(default_factory=lambda: np.empty((0, 0), dtype=float))
    mud_cake_snapshots: Tuple[Array, ...] = field(default_factory=tuple)
    reynolds_snapshots: Tuple[Array, ...] = field(default_factory=tuple)
    turbulent_viscosity_snapshots: Tuple[Array, ...] = field(default_factory=tuple)
    snapshot_times_s: Tuple[float, ...] = field(default_factory=tuple)
    notes: Tuple[str, ...] = field(default_factory=tuple)
    lead_field: Array = field(default_factory=lambda: np.empty((0, 0), dtype=float))
    tail_field: Array = field(default_factory=lambda: np.empty((0, 0), dtype=float))


class AnnulusD2DGASolver:
    """环空二维D2DGA顶替求解器。

    实现偏心环空中水泥浆、前置/隔离液和钻井液的多相体积分数模拟。

    主要物理过程：
    1. 平流输运：水泥浆与前置/隔离液随速度场向下游运移
    2. D2DGA分散：基于Zhang & Frigaard (2022)的通量修正，捕捉间隙尺度速度梯度分散
    3. 弥散效应：横向(day)和垂向(dax)扩散，与偏心度和流体性质相关
    4. 壁面清除：泥饼在剪切作用下被水泥浆与前置/隔离液逐渐清除
    5. 浮力效应：密度差导致窄边优先顶替和潜在失稳风险

    D2DGA通量修正（核心改进）：
    - 假设替浆液占据间隙中心（流速快），被替液贴近壁面（流速慢）
    - 通量放大因子 f(c̄, m) = [m·c̄² + 1.5·(1-c̄²)] / [m·c̄³ + (1-c̄³)]
    - 效果：低浓度时f>1（替浆前锋跑得快），高浓度时f≈1（通量不变）
    - 参考文献：Zhang & Frigaard (2022), JFM Vol.947, A32

    模型特点：
    - 网格：ny×nz（方位角×井深），使用双线性插值实现反演追踪
    - 流变：支持牛顿、Bingham、幂律、Herschel-Bulkley四种模型
    - 泵停处理：泵停后冻结水泥场，避免浮力流"蒸发"已就位水泥
    - 风险指标：计算窜槽、混浆、失稳三个风险指数

    使用示例：
        solver = AnnulusD2DGASolver(dt=4.0, nz=140, ny=40)
        result = solver.run(well_spec, fluids, inlet_state_provider)
    """

    def __init__(
        self,
        *,
        dt: float = 4.0,
        nz: int = 140,
        ny: int = 40,
        alpha_clean: float = 0.085,
        total_t: float = 12000.0,
        enable_d2dga: bool = True,
        d2dga_viscosity_ratio: float = 1.0,
        quality_penalty_scale: float = 0.099,
        channeling_penalty_weight: float = 0.55,
        mixing_penalty_weight: float = 0.35,
        instability_penalty_weight: float = 0.25,
        instability_decay_scale: float = 5.0,
        save_interval: int = 60,
        gel_growth_rate: float = 0.001,
        gel_max_pa: float = 50.0,
        gel_break_threshold: float = 100.0,
        T_surface: float = 20.0,
        geothermal_gradient: float = 0.03,
        T_ref: float = 80.0,
        alpha_T: float = 0.01,
        enable_temperature_coupling: bool = False,
        enable_mud_cake: bool = False,
        initial_mud_cake_mm: float = 3.0,
        k_erosion: float = 0.001,
        enable_turbulence: bool = False,
        Re_critical: float = 2100.0,
        turbulence_coefficient: float = 0.16,
        enable_gravity: bool = False,
        g_constant: float = 9.81,
        gravity_yield_factor: float = 0.5,
    ) -> None:
        """初始化环空二维求解器参数。

        Args:
            dt: 时间步长（秒），默认4秒
            nz: 井深方向网格数，默认140
            ny: 方位角方向网格数，默认40
            alpha_clean: 泥饼清除系数，默认0.085
            total_t: 总模拟时间（秒），默认12000秒（200分钟）
            enable_d2dga: 是否启用D2DGA通量修正（Zhang & Frigaard 2022），默认开启
            d2dga_viscosity_ratio: D2DGA粘度比 m = η_displaced/η_displacing，默认1.0
            quality_penalty_scale: 质量惩罚因子，默认0.099
            channeling_penalty_weight: 窜槽风险权重，默认0.55
            mixing_penalty_weight: 混浆风险权重，默认0.35
            instability_penalty_weight: 失稳风险权重，默认0.25
            instability_decay_scale: 失稳指数衰减标度，默认5.0
            save_interval: 二维场快照保存步长，默认每60个时间步保存一次
            gel_growth_rate: 泵停静置时钻井液凝胶强度增长系数（1/s），默认0.001
            gel_max_pa: 钻井液凝胶强度上限（Pa），默认50Pa
            gel_break_threshold: 泵启后凝胶剪切破坏的剪切速率阈值（1/s），默认100
            T_surface: 地面温度（°C），用于线性地温梯度温度场，默认20°C
            geothermal_gradient: 地温梯度（°C/m），默认0.03°C/m
            T_ref: 流变黏度修正参考温度（°C），默认80°C
            alpha_T: 温度-黏度线性修正系数（1/°C），默认0.01
            enable_temperature_coupling: 是否启用温度-流变耦合；默认关闭以保持历史算例输出不变
            enable_mud_cake: 是否启用泥饼残余层模型；默认关闭以保持历史算例输出不变
            initial_mud_cake_mm: 初始泥饼厚度（mm），典型值2-5mm，默认3mm
            k_erosion: 泥饼侵蚀系数（m/s per 1/s shear rate），默认0.001
            enable_turbulence: 是否启用湍流/混合流修正；默认关闭以保持层流算例历史输出不变
            Re_critical: 临界雷诺数，超过该值时按湍流单元处理，默认2100
            turbulence_coefficient: 湍流系数，等价于混合长度模型中的0.4²，默认0.16
            enable_gravity: 是否启用显式井轴重力体力项；默认关闭以保持历史输出不变
            g_constant: 重力加速度（m/s²），默认9.81
            gravity_yield_factor: 重力对流需克服的屈服应力比例，默认取屈服应力的50%
        """
        self.dt = dt
        self.nz = nz
        self.ny = ny
        self.alpha_clean = alpha_clean
        self.total_t = total_t
        self.enable_d2dga = enable_d2dga
        self.d2dga_viscosity_ratio = d2dga_viscosity_ratio
        self.quality_penalty_scale = quality_penalty_scale
        self.channeling_penalty_weight = channeling_penalty_weight
        self.mixing_penalty_weight = mixing_penalty_weight
        self.instability_penalty_weight = instability_penalty_weight
        self.instability_decay_scale = instability_decay_scale
        self.save_interval: int = save_interval
        # 触变模型参数：仅作用于钻井液相，默认初始凝胶强度为0，不改变连续泵注算例。
        self.gel_growth_rate = gel_growth_rate
        self.gel_max_pa = gel_max_pa
        self.gel_break_threshold = gel_break_threshold
        # 温度模型参数：默认保留但不启用耦合，避免改变既有等温测试结果。
        self.T_surface = T_surface
        self.geothermal_gradient = geothermal_gradient
        self.T_ref = T_ref
        self.alpha_T = alpha_T
        self.enable_temperature_coupling = enable_temperature_coupling
        # 泥饼残余层模型参数：默认关闭，关闭时有效间隙等于原始环空间隙，确保向后兼容。
        self.enable_mud_cake: bool = enable_mud_cake
        self.initial_mud_cake_mm: float = initial_mud_cake_mm
        self.k_erosion: float = k_erosion
        # 湍流模型参数：参考 Maleki & Frigaard (2017) 的混合/湍流流态处理思路。
        # 默认关闭，确保既有层流求解路径与测试输出保持完全兼容。
        self.enable_turbulence: bool = enable_turbulence
        self.Re_critical: float = Re_critical
        self.turbulence_coefficient: float = turbulence_coefficient
        # 显式重力模型参数：默认关闭，避免改变既有经验浮力修正路径和历史测试结果。
        self.enable_gravity: bool = enable_gravity
        self.g_constant: float = g_constant
        self.gravity_yield_factor: float = gravity_yield_factor

    def _build_geom(self, well_spec: WellSpec, mud_cake_thickness: Array | None = None) -> Dict[str, Array]:
        """根据井筒规格构建环空二维网格几何参数。

        构建步骤：
        1. 创建井深方向网格s和深度坐标md
        2. 从井径剖面插值得到井径、井斜、偏心度数据
        3. 构建方位角方向网格y和归一化方位角phi
        4. 计算半间隙h和局部环空间隙b
        5. 校正体积使其等于物理环空体积
        6. 可选扣除泥饼厚度，得到流动计算使用的有效环空间隙

        返回几何参数字典，包含：
        - s: 井深坐标（从鞋口算起）
        - md: 测深坐标
        - y: 方位角坐标
        - phi: 归一化方位角
        - H: 半间隙数组
        - b: 局部环空间隙数组
        - e: 偏心度数组
        - standoff: 居中度数组
        - inc_deg: 井斜角度数组
        - hole_mm: 井径数组
        - od_mm: 尾管外径数组
        - effective_b: 有效环空间隙数组；泥饼模型关闭时等于b
        """
        s = np.linspace(0.0, well_spec.bottom_md_m - well_spec.top_md_m, self.nz)
        md = well_spec.bottom_md_m - s

        cal_md, cal_hole = _profile_to_arrays(well_spec.hole_diameter_profile)
        inc_md, inc_values = _profile_to_arrays(well_spec.inclination_profile)
        standoff_md, standoff_values = _profile_to_arrays(well_spec.standoff_profile)

        hole = np.interp(md, cal_md, cal_hole)
        inc_deg = np.interp(md, inc_md, inc_values)
        standoff = np.interp(md, standoff_md, standoff_values)
        od_mm = np.full_like(md, float(well_spec.liner_od_mm or 0.0), dtype=float)

        e = np.clip(1.0 - standoff, 0.05, 0.55)
        clearance = (hole - od_mm) / 1000.0
        half_gap_mean = clearance / 2.0
        mean_radius = ((hole + od_mm) / 4.0) / 1000.0

        y = np.linspace(0.0, np.pi * np.mean(mean_radius), self.ny)
        phi = y / y[-1]

        h = np.zeros((self.ny, self.nz), dtype=float)
        b = np.zeros((self.ny, self.nz), dtype=float)
        for j in range(self.nz):
            h[:, j] = half_gap_mean[j] * (1.0 + e[j] * np.cos(np.pi * phi))
            b[:, j] = 2.0 * h[:, j]

        geom = {
            "s": s,
            "md": md,
            "y": y,
            "phi": phi,
            "H": h,
            "b": b,
            "e": e,
            "standoff": standoff,
            "inc_deg": inc_deg,
            "hole_mm": hole,
            "od_mm": od_mm,
        }
        current_half_volume = _trapez2d(geom["b"], geom)
        target_half_volume = 0.5 * self._physical_annular_volume(well_spec)
        scale = target_half_volume / current_half_volume
        geom["H"] *= scale
        geom["b"] *= scale
        geom["volume_scale"] = np.array(scale)
        if mud_cake_thickness is None:
            geom["effective_b"] = geom["b"].copy()
        else:
            # 在构建几何参数时考虑泥饼厚度。
            # 有效间隙 = 原始间隙 - 泥饼厚度；当前二维截面厚度场已表示两壁合计等效缩减量。
            effective_b = geom["b"] - mud_cake_thickness
            geom["effective_b"] = np.maximum(effective_b, 0.001)  # 最小有效间隙1mm
        return geom

    def _update_effective_gap(self, geom: Dict[str, Array], mud_cake_thickness: Array) -> None:
        """根据当前泥饼厚度更新有效环空间隙。

        泥饼附着会占据流动通道，使水泥浆实际可通过的环空间隙变小；该字段只用于
        速度、剪切和侵蚀计算，不改变原始几何体积与既有效率指标积分口径。
        """
        if self.enable_mud_cake:
            # 有效间隙 = 原始间隙 - 泥饼厚度；下限1mm防止速度/剪切计算奇异。
            geom["effective_b"] = np.maximum(geom["b"] - mud_cake_thickness, 0.001)
        else:
            geom["effective_b"] = geom["b"]

    def _physical_annular_volume(self, well_spec: WellSpec) -> float:
        """计算井段物理环空体积（用于体积校正）。"""
        cal_md, cal_hole = _profile_to_arrays(well_spec.hole_diameter_profile)
        od = np.full_like(cal_md, float(well_spec.liner_od_mm or 0.0), dtype=float)
        area = np.pi * ((cal_hole / 1000.0) ** 2 - (od / 1000.0) ** 2) / 4.0
        return float(np.trapezoid(area, x=cal_md))

    def _pick_fluids(
        self,
        fluids: Tuple[FluidSpec, ...],
    ) -> Tuple[FluidSpec, FluidSpec | None, FluidSpec | None, FluidSpec | None]:
        """从流体列表中选取钻井液、领浆、尾浆和可选前置/隔离液。"""
        mud = next((fluid for fluid in fluids if fluid.role == FluidRole.MUD), None)
        lead = next((fluid for fluid in fluids if fluid.role == FluidRole.LEAD), None)
        tail = next((fluid for fluid in fluids if fluid.role == FluidRole.TAIL), None)
        spacer = next((fluid for fluid in fluids if fluid.role in {FluidRole.WASH, FluidRole.SPACER}), None)
        if mud is None or (lead is None and tail is None):
            raise ValueError("需要钻井液和至少一个水泥浆流体")
        return mud, lead, tail, spacer

    def _apparent_viscosity(
        self,
        fluid: FluidSpec,
        gamma: Array,
        gel_strength: Array | None = None,
        temperature_correction: Array | None = None,
    ) -> Array:
        """根据流变模型计算流体的表观粘度。

        Args:
            fluid: 流体规格
            gamma: 剪切速率数组
            gel_strength: 钻井液静置凝胶强度场；非钻井液相应保持为None
            temperature_correction: 温度导致的黏度修正系数场；None表示等温不修正

        Returns:
            表观粘度数组
        """
        gamma = np.maximum(np.asarray(gamma, dtype=float), 1.0e-6)
        if fluid.rheology_model == fluid.rheology_model.NEWTONIAN:
            assert fluid.plastic_viscosity_pa_s is not None
            mu = np.full_like(gamma, fluid.plastic_viscosity_pa_s, dtype=float)
        elif fluid.rheology_model == fluid.rheology_model.BINGHAM:
            assert fluid.plastic_viscosity_pa_s is not None
            assert fluid.yield_stress_pa is not None
            mu = fluid.plastic_viscosity_pa_s + fluid.yield_stress_pa / gamma
        elif fluid.rheology_model == fluid.rheology_model.POWER_LAW:
            assert fluid.power_law_n is not None
            assert fluid.consistency_k is not None
            mu = fluid.consistency_k * gamma ** (fluid.power_law_n - 1.0)
        elif fluid.rheology_model == fluid.rheology_model.HERSCHEL_BULKLEY:
            assert fluid.yield_stress_pa is not None
            assert fluid.power_law_n is not None
            assert fluid.consistency_k is not None
            mu = fluid.yield_stress_pa / gamma + fluid.consistency_k * gamma ** (fluid.power_law_n - 1.0)
        else:
            raise ValueError(f"Unsupported rheology model: {fluid.rheology_model}")
        if gel_strength is not None:
            # 触变贡献：泵停形成的凝胶强度以 τ_gel/γ 的形式提高钻井液表观粘度。
            # 该项只在钻井液相调用时传入，避免改变水泥浆/隔离液本身的流变模型。
            mu = mu + np.asarray(gel_strength, dtype=float) / gamma
        if temperature_correction is not None:
            # 温度-流变耦合：温度升高时黏度降低，温度降低时黏度升高；所有流体使用同一修正场。
            mu = mu * np.asarray(temperature_correction, dtype=float)
        return np.clip(mu, 1.0e-5, 3.0)

    def _compute_props(
        self,
        lead: Array,
        tail: Array,
        spacer: Array,
        w_prev: Array,
        geom: Dict[str, Array],
        mud_fluid: FluidSpec,
        lead_fluid: FluidSpec | None,
        tail_fluid: FluidSpec | None,
        spacer_fluid: FluidSpec | None,
        gel_strength: Array | None = None,
        temperature_correction: Array | None = None,
    ) -> Tuple[Array, Array, Array]:
        """计算四相混合物系的表观粘度、密度和钻井液分数。"""
        # 四相体积分数闭合：显式跟踪领浆、尾浆和前置/隔离液，钻井液由守恒关系反算。
        mud = np.clip(1.0 - lead - tail - spacer, 0.0, 1.0)
        effective_b = geom.get("effective_b", geom["b"])
        gamma = np.maximum(6.0 * np.abs(w_prev) / np.maximum(effective_b, 1.0e-5), 1.0e-6)
        # 凝胶强度仅随钻井液相贡献到混合表观粘度，水泥浆和隔离液不受该触变状态影响。
        mu = mud * self._apparent_viscosity(mud_fluid, gamma, gel_strength, temperature_correction)
        if lead_fluid is not None:
            mu += lead * self._apparent_viscosity(lead_fluid, gamma, temperature_correction=temperature_correction)
        if tail_fluid is not None:
            mu += tail * self._apparent_viscosity(tail_fluid, gamma, temperature_correction=temperature_correction)
        if spacer_fluid is not None:
            mu += spacer * self._apparent_viscosity(spacer_fluid, gamma, temperature_correction=temperature_correction)
        rho = mud * (mud_fluid.density_kg_m3 / 1000.0)
        if lead_fluid is not None:
            rho += lead * (lead_fluid.density_kg_m3 / 1000.0)
        if tail_fluid is not None:
            rho += tail * (tail_fluid.density_kg_m3 / 1000.0)
        if spacer_fluid is not None:
            rho += spacer * (spacer_fluid.density_kg_m3 / 1000.0)
        return mu, rho, mud

    def _compute_velocity(
        self,
        lead: Array,
        tail: Array,
        spacer: Array,
        geom: Dict[str, Array],
        q_m3s: float,
        w_prev: Array,
        mud_fluid: FluidSpec,
        lead_fluid: FluidSpec | None,
        tail_fluid: FluidSpec | None,
        spacer_fluid: FluidSpec | None,
        gel_strength: Array | None = None,
        temperature_correction: Array | None = None,
    ) -> Tuple[Array, Array, Array, Array, Array, Array, Array]:
        """计算环空速度场。

        计算步骤：
        1. 根据密度差计算浮力修正因子
        2. 计算流度场m（考虑间隙、偏心度、居中度）
        3. 对窄边进行流度增强修正
        4. 由总排量积分得到井深方向速度w
        5. 通过连续性方程求解方位角方向速度v

        Returns:
            w: 井深方向速度（轴向速度）
            v: 方位角方向速度（横向速度）
            mu: 有效表观粘度场；启用湍流时为层流粘度+湍流粘度，否则等于层流粘度
            rho: 密度场
            mud: 钻井液分数场
            Re: 雷诺数场
            mu_turbulent: 湍流粘度修正场
        """
        y = geom["y"]
        # 在速度计算中使用有效间隙；泥饼模型关闭时 effective_b 与原始 b 完全一致。
        effective_b = geom.get("effective_b", geom["b"])
        b = effective_b
        q_half = q_m3s / 2.0
        mu, rho, mud = self._compute_props(
            lead,
            tail,
            spacer,
            w_prev,
            geom,
            mud_fluid,
            lead_fluid,
            tail_fluid,
            spacer_fluid,
            gel_strength,
            temperature_correction,
        )
        mud_density_gcc = mud_fluid.density_kg_m3 / 1000.0
        cement = np.clip(lead + tail, 0.0, 1.0)

        shear_rate = np.maximum(6.0 * np.abs(w_prev) / np.maximum(effective_b, 1.0e-5), 1.0e-6)
        # Reynolds number calculation
        # Re = ρ * v * D_h / μ
        # 其中 D_h 是水力直径，对于环空 D_h = 2 * (r_outer - r_inner)
        D_h = 2.0 * geom["b"]  # 水力直径 (m)
        rho_kg_m3 = rho * 1000.0  # 模型内部rho以g/cm³保存，Re与湍流粘度计算需换算为kg/m³。
        Re = rho_kg_m3 * np.abs(w_prev) * D_h / np.maximum(mu, 1.0e-6)

        if self.enable_turbulence:
            # 湍流粘度修正：基于混合长度模型
            # μ_turb = ρ * l_m^2 * |du/dy|
            # 其中 l_m 是混合长度，对于环空 l_m ≈ 0.4 * y (y是到壁面的距离)
            # 简化为：μ_turb = ρ * (0.4 * b)^2 * shear_rate
            # 总粘度 = 层流粘度 + 湍流粘度
            # turbulence_coefficient 默认0.16，对应 (0.4)^2；仅在 Re 超临界时生效。
            turbulent_mask = Re > self.Re_critical
            mu_turbulent_raw = rho_kg_m3 * self.turbulence_coefficient * geom["b"]**2 * shear_rate
            mu_turbulent = np.where(turbulent_mask, mu_turbulent_raw, 0.0)
            mu_total = np.where(turbulent_mask, mu + mu_turbulent, mu)
        else:
            # 兼容模式：仍计算Re供诊断快照使用，但不改变速度场粘度。
            mu_turbulent = np.zeros_like(mu)
            mu_total = mu
        # 在计算速度时使用总粘度（层流+湍流）
        mu_effective = mu_total

        buoy = 1.0 + 0.60 * np.maximum(rho - mud_density_gcc, 0.0) + 0.08 * np.sin(np.deg2rad(geom["inc_deg"]))[None, :]
        m = (b**3 / np.maximum(mu_effective, 1.0e-6)) * buoy * (0.90 + 0.25 * geom["standoff"][None, :])

        phi = geom["phi"][:, None]
        narrow_boost = 1.0 - 0.30 * cement * np.maximum(rho - mud_density_gcc, 0.0) * np.cos(np.pi * phi)
        m *= np.clip(narrow_boost, 0.75, 1.70)

        int_m = np.trapezoid(b * m, x=y, axis=0)
        w = (q_half / np.maximum(int_m, 1.0e-12))[None, :] * m

        if self.enable_gravity:
            # 在斜井/水平井段中，重力沿井轴分量为 g·sin(θ)，θ为从垂直方向测量的井斜角。
            # 对Hele-Shaw窄缝流动，体力驱动速度近似为 (b²/12μ)·ρ·g·sin(θ)。
            # 这里保留原有排量归一化速度，同时叠加密度相关的显式重力体力项。
            theta = np.deg2rad(geom["inc_deg"])[None, :]
            w_gravity = (
                b**2
                / (12.0 * np.maximum(mu_effective, 1.0e-6))
                * rho_kg_m3
                * self.g_constant
                * np.sin(theta)
            )
            # 减去面积加权平均重力速度，使叠加项不改变每个井深截面的总排量约束。
            gravity_flux = np.trapezoid(b * w_gravity, x=y, axis=0)
            area_flux = np.trapezoid(b, x=y, axis=0)
            w_gravity -= (gravity_flux / np.maximum(area_flux, 1.0e-12))[None, :]
            w = w + w_gravity

        ds = geom["s"][1] - geom["s"][0]
        dy = y[1] - y[0]
        bw = b * w
        dbw_ds = np.gradient(bw, ds, axis=1)
        bv = np.zeros_like(w)
        for i in range(1, len(y)):
            bv[i, :] = bv[i - 1, :] - 0.5 * (dbw_ds[i, :] + dbw_ds[i - 1, :]) * dy
        bv -= (y[:, None] / y[-1]) * bv[-1, :]
        v = bv / np.maximum(b, 1.0e-8)
        return w, v, mu_effective, rho, mud, Re, mu_turbulent

    def _depth_profiles(self, geom: Dict[str, Array], lead: Array, tail: Array, spacer: Array, wall: Array) -> pd.DataFrame:
        """计算深度方向的平均剖面数据。"""
        cement = np.clip(lead + tail, 0.0, 1.0)
        eff = cement * (1.0 - wall)
        mud = np.clip(1.0 - lead - tail - spacer, 0.0, 1.0)
        return pd.DataFrame(
            {
                "井深_m": geom["md"],
                "领浆平均浓度": np.average(lead, axis=0, weights=geom["b"]),
                "尾浆平均浓度": np.average(tail, axis=0, weights=geom["b"]),
                "水泥平均浓度": np.average(cement, axis=0, weights=geom["b"]),
                "前置液隔离液平均浓度": np.average(spacer, axis=0, weights=geom["b"]),
                "平均有效顶替效率": np.average(eff, axis=0, weights=geom["b"]),
                "钻井液平均浓度": np.average(mud, axis=0, weights=geom["b"]),
                "宽边有效效率": eff[0],
                "中线有效效率": eff[eff.shape[0] // 2],
                "窄边有效效率": eff[-1],
                "宽边水泥浓度": cement[0],
                "中线水泥浓度": cement[cement.shape[0] // 2],
                "窄边水泥浓度": cement[-1],
                "宽边前置液隔离液浓度": spacer[0],
                "中线前置液隔离液浓度": spacer[spacer.shape[0] // 2],
                "窄边前置液隔离液浓度": spacer[-1],
                "环空间隙_m": np.mean(geom["b"], axis=0),
                "偏心度指标": geom["e"],
                "居中度": geom["standoff"],
            }
        )

    def run(
        self,
        well_spec: WellSpec,
        fluids: Tuple[FluidSpec, ...],
        inlet_state_provider: Callable[[float], AnnulusInletState],
    ) -> AnnulusSimulationResult:
        """运行 Hu102 尾管段首版环空二维顶替求解。"""

        mud_fluid, lead_fluid, tail_fluid, spacer_fluid = self._pick_fluids(fluids)
        geom = self._build_geom(well_spec)
        lead = np.zeros((self.ny, self.nz), dtype=float)
        tail = np.zeros((self.ny, self.nz), dtype=float)
        spacer = np.zeros((self.ny, self.nz), dtype=float)
        wall = np.ones((self.ny, self.nz), dtype=float)
        # 触变状态场：记录钻井液在泵停静置期间形成的凝胶强度，初始为0以保持正常泵注结果不变。
        gel_strength = np.zeros((self.ny, self.nz), dtype=float)
        # 在 run() 开始时初始化温度场。
        # 假设线性地温梯度：T = T_surface + gradient * depth。
        # 这里的 depth 采用测深 md；温度场在当前简化模型中为准稳态场，随井深变化并在整个时间推进中复用。
        temperature = self.T_surface + self.geothermal_gradient * geom["md"][None, :]
        temperature = np.broadcast_to(temperature, (self.ny, self.nz)).astype(float)
        # 温度对粘度的影响：Arrhenius型关系
        # μ(T) = μ_ref * exp(Ea/R * (1/T - 1/T_ref))
        # 简化为线性近似：μ(T) = μ_ref * (1 + alpha_T * (T_ref - T))
        if self.enable_temperature_coupling:
            viscosity_correction = 1.0 + self.alpha_T * (self.T_ref - temperature)
            viscosity_correction = np.clip(viscosity_correction, 0.5, 2.0)
        else:
            # 默认等温兼容模式：修正系数恒为1，确保历史测试和输出数值不变。
            viscosity_correction = np.ones_like(temperature, dtype=float)
        w_prev = np.full((self.ny, self.nz), 0.45, dtype=float)
        half_volume = _trapez2d(geom["b"], geom)

        target_mask = _window_mask(well_spec, geom["md"], "target")
        cbl_mask = _window_mask(well_spec, geom["md"], "cbl")
        ygrid, sgrid = np.meshgrid(geom["y"], geom["s"], indexing="ij")
        rows: list[list[float | str]] = []
        mud_density_gcc = mud_fluid.density_kg_m3 / 1000.0
        cement_snapshots: list[Array] = []
        lead_snapshots: list[Array] = []
        tail_snapshots: list[Array] = []
        spacer_snapshots: list[Array] = []
        wall_snapshots: list[Array] = []
        gel_strength_snapshots: list[Array] = []
        mud_cake_snapshots: list[Array] = []
        reynolds_snapshots: list[Array] = []
        turbulent_viscosity_snapshots: list[Array] = []
        snapshot_times: list[float] = []

        # 在 run() 开始时初始化泥饼厚度场。
        # 初始泥饼厚度：假设均匀分布，典型值 2-5mm；关闭模型时置零以保持既有输出不变。
        initial_mud_cake_m = self.initial_mud_cake_mm / 1000.0 if self.enable_mud_cake else 0.0
        mud_cake_thickness: Array = np.full((self.ny, self.nz), initial_mud_cake_m, dtype=float)  # 转换为m
        self._update_effective_gap(geom, mud_cake_thickness)

        final_step_index = int(self.total_t / self.dt)
        for step_index in range(final_step_index + 1):
            current_time_s = step_index * self.dt
            inlet_state = inlet_state_provider(current_time_s)
            inlet_cement_fraction = _phase_fraction(inlet_state, "cement")
            inlet_lead_fraction = _phase_fraction(inlet_state, "lead")
            inlet_tail_fraction = _phase_fraction(inlet_state, "tail") + inlet_cement_fraction
            inlet_spacer_fraction = _phase_fraction(inlet_state, "spacer")

            # 泵停判断：排量低于阈值时认为泵已停止。
            # 泵停后水泥场冻结——不再平流、扩散或壁面清除，
            # 因为水泥静凝胶强度在短时间（<2h）内足以抵抗浮力滑塌。
            pump_active = inlet_state.flow_rate_m3_s > 1.0e-9
            self._update_effective_gap(geom, mud_cake_thickness)
            effective_b = geom["effective_b"]

            if pump_active:
                # === 正常泵注阶段：执行全部物理过程 ===
                w, v, mu, rho, mud, Re, mu_turbulent = self._compute_velocity(
                    lead,
                    tail,
                    spacer,
                    geom,
                    inlet_state.flow_rate_m3_s,
                    w_prev,
                    mud_fluid,
                    lead_fluid,
                    tail_fluid,
                    spacer_fluid,
                    gel_strength,
                    viscosity_correction,
                )
                # 泵启/泵注时：钻井液凝胶结构被剪切破坏。
                # 剪切速率超过阈值的位置快速破胶，低剪切区域保留部分弱凝胶结构。
                shear_rate = np.maximum(6.0 * np.abs(w) / np.maximum(effective_b, 1.0e-5), 1.0e-6)
                gel_strength *= np.where(shear_rate > self.gel_break_threshold, 0.1, 0.5)
                w_prev = w

                cement = np.clip(lead + tail, 0.0, 1.0)
                buoy_num = np.maximum(rho - mud_density_gcc, 0.0)
                disp_ratio = np.clip(mud / (cement + 1.0e-4), 0.2, 8.0)
                dax = 0.0040 * np.abs(w) * effective_b * (0.20 + 1.20 * cement * (1.0 - cement)) * (1.0 + 0.50 * geom["e"][None, :]) * (1.0 + 0.10 * np.maximum(disp_ratio - 1.0, 0.0)) / (1.0 + 0.95 * buoy_num)
                day = 0.16 * dax * (0.45 + 0.55 * geom["e"][None, :])

                # === D2DGA通量修正（Zhang & Frigaard 2022, JFM Vol.947, A32）===
                # 论文核心：假设替浆液占据间隙中心（流速快），被替液贴近壁面（流速慢），
                # 不同浓度层以不同速度运动，产生分散效应。
                # 通量放大因子 f(c̄, m) = [m·c̄² + 1.5·(1-c̄²)] / [m·c̄³ + (1-c̄³)]
                # 其中 c̄ = 间隙平均浓度（水泥浆），m = 粘度比（被替/替）
                # 效果：c̄小→f>1（替浆前锋跑得快），c̄大→f≈1（通量不变）
                if self.enable_d2dga:
                    m_ratio = self.d2dga_viscosity_ratio
                    c_safe = np.clip(cement, 0.01, 0.99)  # 避免除零
                    f_amp = (m_ratio * c_safe**2 + 1.5 * (1.0 - c_safe**2)) / (m_ratio * c_safe**3 + (1.0 - c_safe**3))
                    f_amp = np.clip(f_amp, 0.5, 2.0)  # 限制放大范围，保持数值稳定
                else:
                    f_amp = 1.0

                v_eff = v + 0.007 * np.maximum(rho - mud_density_gcc, 0.0) * cement * np.sin(np.pi * geom["phi"])[:, None]
                # D2DGA修正：用放大因子调整有效速度场
                w_d2dga = w * f_amp
                v_eff_d2dga = v_eff * f_amp
                ysrc = ygrid - v_eff_d2dga * self.dt
                ssrc = sgrid - w_d2dga * self.dt
                # 领浆、尾浆与前置/隔离液使用同一速度场输运；只改变被输运的浓度变量。
                lead_adv = _bilinear_interp(lead, ysrc, ssrc, geom, inlet_lead_fraction)
                tail_adv = _bilinear_interp(tail, ysrc, ssrc, geom, inlet_tail_fraction)
                spacer_adv = _bilinear_interp(spacer, ysrc, ssrc, geom, inlet_spacer_fraction)

                dy = geom["y"][1] - geom["y"][0]
                ds = geom["s"][1] - geom["s"][0]
                lead_lap_s = np.zeros_like(lead_adv)
                lead_lap_y = np.zeros_like(lead_adv)
                tail_lap_s = np.zeros_like(tail_adv)
                tail_lap_y = np.zeros_like(tail_adv)
                spacer_lap_s = np.zeros_like(spacer_adv)
                spacer_lap_y = np.zeros_like(spacer_adv)
                lead_lap_s[:, 1:-1] = (lead_adv[:, 2:] - 2.0 * lead_adv[:, 1:-1] + lead_adv[:, :-2]) / ds**2
                lead_lap_s[:, 0] = (lead_adv[:, 1] - lead_adv[:, 0]) / ds**2
                lead_lap_s[:, -1] = (lead_adv[:, -2] - lead_adv[:, -1]) / ds**2
                lead_lap_y[1:-1, :] = (lead_adv[2:, :] - 2.0 * lead_adv[1:-1, :] + lead_adv[:-2, :]) / dy**2
                lead_lap_y[0, :] = (lead_adv[1, :] - lead_adv[0, :]) / dy**2
                lead_lap_y[-1, :] = (lead_adv[-2, :] - lead_adv[-1, :]) / dy**2
                tail_lap_s[:, 1:-1] = (tail_adv[:, 2:] - 2.0 * tail_adv[:, 1:-1] + tail_adv[:, :-2]) / ds**2
                tail_lap_s[:, 0] = (tail_adv[:, 1] - tail_adv[:, 0]) / ds**2
                tail_lap_s[:, -1] = (tail_adv[:, -2] - tail_adv[:, -1]) / ds**2
                tail_lap_y[1:-1, :] = (tail_adv[2:, :] - 2.0 * tail_adv[1:-1, :] + tail_adv[:-2, :]) / dy**2
                tail_lap_y[0, :] = (tail_adv[1, :] - tail_adv[0, :]) / dy**2
                tail_lap_y[-1, :] = (tail_adv[-2, :] - tail_adv[-1, :]) / dy**2
                spacer_lap_s[:, 1:-1] = (spacer_adv[:, 2:] - 2.0 * spacer_adv[:, 1:-1] + spacer_adv[:, :-2]) / ds**2
                spacer_lap_s[:, 0] = (spacer_adv[:, 1] - spacer_adv[:, 0]) / ds**2
                spacer_lap_s[:, -1] = (spacer_adv[:, -2] - spacer_adv[:, -1]) / ds**2
                spacer_lap_y[1:-1, :] = (spacer_adv[2:, :] - 2.0 * spacer_adv[1:-1, :] + spacer_adv[:-2, :]) / dy**2
                spacer_lap_y[0, :] = (spacer_adv[1, :] - spacer_adv[0, :]) / dy**2
                spacer_lap_y[-1, :] = (spacer_adv[-2, :] - spacer_adv[-1, :]) / dy**2
                lead = np.clip(
                    lead_adv + self.dt * (dax * lead_lap_s + day * lead_lap_y),
                    0.0,
                    1.0,
                )
                tail = np.clip(
                    tail_adv + self.dt * (dax * tail_lap_s + day * tail_lap_y),
                    0.0,
                    1.0,
                )
                spacer = np.clip(
                    spacer_adv + self.dt * (dax * spacer_lap_s + day * spacer_lap_y),
                    0.0,
                    1.0,
                )
                # 数值扩散可能使显式两相之和略超1；按比例压回可行域，保持泥浆分数非负。
                tracked_total = lead + tail + spacer
                overfilled = tracked_total > 1.0
                lead[overfilled] /= tracked_total[overfilled]
                tail[overfilled] /= tracked_total[overfilled]
                spacer[overfilled] /= tracked_total[overfilled]
                cement = np.clip(lead + tail, 0.0, 1.0)

                # 清洗能力由水泥浆主导，隔离液/冲洗液也参与泥饼清除但权重较低。
                cleaner = 1.10 * (cement + 0.8 * spacer)
                shear = np.abs(w) / np.maximum(effective_b, 1.0e-5)
                kclean = self.alpha_clean * cleaner * (0.45 + np.sqrt(np.maximum(shear, 1.0e-6))) * (0.85 + 0.35 * geom["standoff"][None, :]) * (1.0 + 0.30 * buoy_num)
                wall *= np.exp(-kclean * self.dt / 150.0)
                wall = np.clip(wall, 0.0, 1.0)

                if self.enable_mud_cake:
                    # 泥饼侵蚀模型：水泥浆剪切作用下泥饼逐渐被清除。
                    # 侵蚀速率 = k_erosion * shear_rate * cement_concentration。
                    shear_rate = np.maximum(6.0 * np.abs(w) / np.maximum(effective_b, 1.0e-5), 1.0e-6)
                    erosion_rate = self.k_erosion * shear_rate * cement
                    mud_cake_thickness -= erosion_rate * self.dt
                    mud_cake_thickness = np.maximum(mud_cake_thickness, 0.0)  # 不能为负
                    self._update_effective_gap(geom, mud_cake_thickness)

            else:
                # === 泵停阶段：冻结水泥场与壁面场，仅记录指标 ===
                # 泵停期间：钻井液触变凝胶强度随静置时间增长。
                # 典型钻井液凝胶强度增长模型：τ_gel = τ_gel0 * (1 + α * t_gel)
                # 这里以显式增量近似累计静置效应，α由 gel_growth_rate 控制，并限制最大凝胶强度。
                gel_strength += self.gel_growth_rate * self.dt
                gel_strength = np.clip(gel_strength, 0.0, self.gel_max_pa)
                # 使用上一步的 mu/rho/mud 计算 mobility 等指标（用于风险指数），
                # 但不更新 cement 或 wall，避免浮力流"蒸发"已就位的水泥。
                w, v, mu, rho, mud, Re, mu_turbulent = self._compute_velocity(
                    lead,
                    tail,
                    spacer,
                    geom,
                    0.0,
                    w_prev,
                    mud_fluid,
                    lead_fluid,
                    tail_fluid,
                    spacer_fluid,
                    gel_strength,
                    viscosity_correction,
                )
                # 注意：w_prev 保持泵停前最后的值，不再更新，
                # 因泵停后速度场已冻结。
                if self.enable_gravity:
                    # 停泵期间Q=0，密度差形成的井轴重力分量成为唯一可能的驱动力。
                    # 驱动剪切应力近似为 |Δρ|·g·sin(θ)·b/2；若小于屈服应力阈值，则认为凝胶/屈服应力锁住流体。
                    theta = np.deg2rad(geom["inc_deg"])[None, :]
                    density_delta_kg_m3 = (rho - mud_density_gcc) * 1000.0
                    gravity_stress = (
                        np.abs(density_delta_kg_m3)
                        * self.g_constant
                        * np.sin(theta)
                        * effective_b
                        / 2.0
                    )
                    if mud_fluid.yield_stress_pa is not None:
                        yield_threshold = mud_fluid.yield_stress_pa * self.gravity_yield_factor
                    else:
                        yield_threshold = 0.0
                    flow_mask = gravity_stress > yield_threshold
                    # 用Hele-Shaw体力速度估算停泵沉降/上浮速度；正负号由密度差决定。
                    w_settle = (
                        effective_b**2
                        / (12.0 * np.maximum(mu, 1.0e-6))
                        * density_delta_kg_m3
                        * self.g_constant
                        * np.sin(theta)
                    )
                    w_settle = np.where(flow_mask, w_settle, 0.0)
                    # 停泵对流是缓慢再分布过程，限制最大速度以避免显式反演追踪数值不稳定。
                    w_settle = np.clip(w_settle, -0.01, 0.01)
                    # 仅做轴向平流，不额外加入扩散或壁面清除，避免把停泵重力滑移误当成继续顶替效率增长。
                    lead = _bilinear_interp(lead, ygrid, sgrid - w_settle * self.dt, geom, inlet_lead_fraction)
                    tail = _bilinear_interp(tail, ygrid, sgrid - w_settle * self.dt, geom, inlet_tail_fraction)
                    spacer = _bilinear_interp(spacer, ygrid, sgrid - w_settle * self.dt, geom, inlet_spacer_fraction)
                    lead = np.clip(lead, 0.0, 1.0)
                    tail = np.clip(tail, 0.0, 1.0)
                    spacer = np.clip(spacer, 0.0, 1.0)
                    tracked_total = lead + tail + spacer
                    overfilled = tracked_total > 1.0
                    lead[overfilled] /= tracked_total[overfilled]
                    tail[overfilled] /= tracked_total[overfilled]
                    spacer[overfilled] /= tracked_total[overfilled]
                    cement = np.clip(lead + tail, 0.0, 1.0)

            # 在物理场更新后、指标计算前保存快照，确保快照与本步指标使用同一状态。
            # 使用 copy() 固化二维场，避免后续时间步原地更新影响已保存结果。
            if step_index % self.save_interval == 0 or step_index == int(self.total_t / self.dt):
                cement = np.clip(lead + tail, 0.0, 1.0)
                cement_snapshots.append(cement.copy())
                lead_snapshots.append(lead.copy())
                tail_snapshots.append(tail.copy())
                spacer_snapshots.append(spacer.copy())
                wall_snapshots.append(wall.copy())
                gel_strength_snapshots.append(gel_strength.copy())
                mud_cake_snapshots.append(mud_cake_thickness.copy())
                # 保存雷诺数和湍流粘度二维场，供后续图表诊断流态分布与湍流修正强度。
                reynolds_snapshots.append(Re.copy())
                turbulent_viscosity_snapshots.append(mu_turbulent.copy())
                snapshot_times.append(current_time_s)

            cement = np.clip(lead + tail, 0.0, 1.0)
            eff = cement * (1.0 - wall)
            bulk_fill = _trapez2d(geom["b"] * cement, geom) / half_volume
            effective_efficiency = _trapez2d(geom["b"] * eff, geom) / half_volume
            target_efficiency = _trapez2d(geom["b"] * eff * target_mask, geom) / max(_trapez2d(geom["b"] * target_mask, geom), 1.0e-12)
            cbl_efficiency = _trapez2d(geom["b"] * eff * cbl_mask, geom) / max(_trapez2d(geom["b"] * cbl_mask, geom), 1.0e-12)

            def _front(line: Array, threshold: float = 0.5) -> float:
                idx = np.where(line >= threshold)[0]
                return float(geom["s"][idx.max()]) if idx.size else 0.0

            front_wide = _front(cement[0])
            front_narrow = _front(cement[-1])
            front_mid = _front(cement[self.ny // 2])
            channeling = abs(front_wide - front_narrow) / (geom["s"][-1] + 1.0e-9)
            mixing = _trapez2d(geom["b"] * (4.0 * cement * (1.0 - cement)), geom) / half_volume
            mobility_wide = np.mean((geom["b"][0] ** 3) / (mu[0] + 1.0e-6))
            mobility_narrow = np.mean((geom["b"][-1] ** 3) / (mu[-1] + 1.0e-6))
            instability_proxy = channeling * max(mobility_wide / (mobility_narrow + 1.0e-9) - 1.0, 0.0) * (1.0 + 0.4 * mixing)
            instability_index = 1.0 - np.exp(-instability_proxy / self.instability_decay_scale)
            raw_penalty = (
                self.channeling_penalty_weight * channeling
                + self.mixing_penalty_weight * mixing
                + self.instability_penalty_weight * instability_index
            )
            quality_factor = np.clip(1.0 - self.quality_penalty_scale * raw_penalty, 0.0, 1.0)
            quality_proxy = cbl_efficiency * quality_factor

            rows.append(
                [
                    current_time_s,
                    current_time_s / 60.0,
                    inlet_state.stage_name,
                    bulk_fill,
                    effective_efficiency,
                    target_efficiency,
                    cbl_efficiency,
                    quality_proxy,
                    front_wide,
                    front_narrow,
                    front_mid,
                    channeling,
                    mixing,
                    instability_proxy,
                    instability_index,
                    float(np.mean(wall)),
                    float(np.mean(cement)),
                    float(np.mean(mud)),
                ]
            )

        metric_columns = [
            "time_s",
            "time_min",
            "stage",
            "bulk_cement_fill",
            "effective_efficiency",
            "target_interval_efficiency",
            "cbl_eval_interval_efficiency",
            "cbl_quality_proxy",
            "front_wide_m",
            "front_narrow_m",
            "front_mid_m",
            "channeling_index",
            "mixing_index",
            "instability_proxy",
            "instability_index",
            "mean_wall_mud",
            "mean_cement",
            "mean_mud",
        ]
        metrics = pd.DataFrame(data=rows, columns=pd.Index(metric_columns))
        cement = np.clip(lead + tail, 0.0, 1.0)
        depth_profiles = self._depth_profiles(geom, lead, tail, spacer, wall)
        final = metrics.iloc[-1]
        summary: Dict[str, object] = {
            "模型名称": "通用尾管段环空二维顶替模型",
            "模拟对象": f"{well_spec.well_name} 尾管段 {well_spec.top_md_m:.2f}-{well_spec.bottom_md_m:.2f}m",
            "井段_m": [well_spec.top_md_m, well_spec.bottom_md_m],
            "物理环空体积_m3": self._physical_annular_volume(well_spec),
            "最终结果": {
                "全井段最终有效顶替效率": float(final["effective_efficiency"]),
                "CBL评价井段模拟有效顶替效率": float(final["cbl_eval_interval_efficiency"]),
                "目标层段模拟有效顶替效率": float(final["target_interval_efficiency"]),
                "最终水泥浆占据率": float(final["bulk_cement_fill"]),
                "最终质量响应效率": float(final["cbl_quality_proxy"]),
                "最终窜槽指数": float(final["channeling_index"]),
                "最终混浆指数": float(final["mixing_index"]),
                "最终失稳指数": float(final["instability_index"]),
            },
        }
        return AnnulusSimulationResult(
            well_name=well_spec.well_name,
            geom=geom,
            cement_field=cement,
            spacer_field=spacer,
            wall_field=wall,
            metrics=metrics,
            depth_profiles=depth_profiles,
            summary=summary,
            time_points_s=tuple(float(value) for value in metrics["time_s"].to_list()),
            cement_snapshots=tuple(cement_snapshots),
            lead_snapshots=tuple(lead_snapshots),
            tail_snapshots=tuple(tail_snapshots),
            spacer_snapshots=tuple(spacer_snapshots),
            wall_snapshots=tuple(wall_snapshots),
            gel_strength_snapshots=tuple(gel_strength_snapshots),
            mud_cake_field=mud_cake_thickness,
            mud_cake_snapshots=tuple(mud_cake_snapshots),
            reynolds_snapshots=tuple(reynolds_snapshots),
            turbulent_viscosity_snapshots=tuple(turbulent_viscosity_snapshots),
            snapshot_times_s=tuple(snapshot_times),
            notes=(
                "当前显式跟踪领浆、尾浆与前置液/隔离液三类入环空相，钻井液由体积分数闭合反算。",
                "效率与质量指标按总水泥浆相（领浆+尾浆）计算，前置液/隔离液仅影响流变混合、输运占据和壁面清洗。",
                "触变凝胶强度仅作用于钻井液表观粘度：泵停增长，泵启剪切破坏，不改变水泥浆输运/扩散公式。",
                "温度-流变耦合使用线性地温梯度与黏度线性修正；默认等温兼容模式下修正系数为1。",
                "泥饼残余层模型默认关闭；启用后泥饼厚度会缩小有效环空间隙，并随水泥浆剪切侵蚀逐步降低。",
                "湍流/混合流态修正基于临界雷诺数与混合长度粘度；默认关闭，开启后仅在Re超过阈值处增加有效粘度。",
            ),
            lead_field=lead,
            tail_field=tail,
        )
