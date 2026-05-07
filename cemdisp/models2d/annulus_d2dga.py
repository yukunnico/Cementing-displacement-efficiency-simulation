"""
多井通用环空二维 D2DGA 求解器。

本模块实现面向 ``cemdisp`` 主链路的偏心环空二维 D2DGA 核心，口径尽量
收敛到 Zhang & Frigaard (2022) 所强调的环空层流顶替主过程：

1. 偏心窄环空几何展开；
2. 基于局部流动度的轴向/方位角平均速度场；
3. D2DGA 通量放大修正，近似捕捉间隙尺度分散；
4. 仅输出求解域内的顶替效率与浓度场。

注意：
- 本模块不再把泥饼、温度、凝胶强度、湍流修正、CBL 质量惩罚等工程扩展
  项作为核心求解的一部分；
- 为兼顾下游脚本兼容性，这些旧参数与旧快照字段仍保留接口或占位输出，
  但不再影响求解结果；
- 现场 CBL/合格率只允许用于验证与对比，不允许反向校准求解器。
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
from cemdisp.models2d.d2dga_flux import d2dga_flux_amplification


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


def _phase_volume(field: Array, geom: Dict[str, Array]) -> float:
    """计算某一相在全环空中的实际占据体积。

    当前求解器在 ``y`` 方向按半环空展开，因此这里对半环空积分结果乘以 2，
    与真实全环空体积口径保持一致。
    """
    return 2.0 * _trapez2d(geom["b"] * np.clip(field, 0.0, 1.0), geom)


def _limit_phase_volume(field: Array, geom: Dict[str, Array], target_volume_m3: float) -> Array:
    """按累计入环空体积限制场量，避免数值扩散凭空放大相体积。"""
    if target_volume_m3 <= 0.0:
        return np.zeros_like(field)
    current_volume_m3 = _phase_volume(field, geom)
    if current_volume_m3 <= target_volume_m3 + 1.0e-9 or current_volume_m3 <= 1.0e-12:
        return np.clip(field, 0.0, 1.0)
    return np.clip(field * (target_volume_m3 / current_volume_m3), 0.0, 1.0)


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
    1. 平流输运：水泥浆与前置/隔离液随平均速度场向下游运移；
    2. D2DGA 分散：基于 Zhang & Frigaard (2022) 的通量放大修正，近似捕捉间隙尺度分散；
    3. 浮力相关横向再分布：以 Hele-Shaw 风格的密度差横向速度近似，保持密度差影响仍在核心层；
    4. 仅输出求解域指标，不在 solver 内叠加现场质量惩罚或 CBL 校准。

    D2DGA通量修正（核心改进）：
    - 假设替浆液占据间隙中心（流速快），被替液贴近壁面（流速慢）
    - 通量放大因子 f(c̄, m) = [m·c̄² + 1.5·(1-c̄²)] / [m·c̄³ + (1-c̄³)]
    - 效果：低浓度时f>1（替浆前锋跑得快），高浓度时f≈1（通量不变）
    - 参考文献：Zhang & Frigaard (2022), JFM Vol.947, A32

    模型特点：
    - 网格：ny×nz（方位角×井深），使用双线性插值实现反演追踪；
    - 流变：继续兼容牛顿、Bingham、幂律、Herschel-Bulkley 四种现场输入；
    - 泵停处理：泵停后冻结浓度场，不再加入停泵重力滑移等额外扩展；
    - 风险指标：仅作为后验诊断输出，不反向影响顶替结果。

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
            total_t: 总模拟时间（秒），默认12000秒（200分钟）
            enable_d2dga: 是否启用D2DGA通量修正（Zhang & Frigaard 2022），默认开启
            d2dga_viscosity_ratio: D2DGA粘度比 m = η_displaced/η_displacing，默认1.0
            quality_penalty_scale: 兼容旧脚本保留，但论文口径核心中不再使用。
            channeling_penalty_weight: 兼容旧脚本保留，仅用于后验风险指标。
            mixing_penalty_weight: 兼容旧脚本保留，仅用于后验风险指标。
            instability_penalty_weight: 兼容旧脚本保留，仅用于后验风险指标。
            instability_decay_scale: 后验失稳指数缩放，默认5.0。
            save_interval: 二维场快照保存步长，默认每60个时间步保存一次
            gel_growth_rate: 兼容旧脚本保留，核心中不再使用。
            gel_max_pa: 兼容旧脚本保留，核心中不再使用。
            gel_break_threshold: 兼容旧脚本保留，核心中不再使用。
            T_surface: 兼容旧脚本保留，核心中不再使用。
            geothermal_gradient: 兼容旧脚本保留，核心中不再使用。
            T_ref: 兼容旧脚本保留，核心中不再使用。
            alpha_T: 兼容旧脚本保留，核心中不再使用。
            enable_temperature_coupling: 兼容旧脚本保留，核心中不再使用。
            enable_mud_cake: 兼容旧脚本保留，核心中不再使用。
            initial_mud_cake_mm: 兼容旧脚本保留，核心中不再使用。
            k_erosion: 兼容旧脚本保留，核心中不再使用。
            enable_turbulence: 兼容旧脚本保留，核心中不再使用。
            Re_critical: 兼容旧脚本保留，仅用于诊断雷诺数阈值注释。
            turbulence_coefficient: 兼容旧脚本保留，核心中不再使用。
            enable_gravity: 兼容旧脚本保留，核心中不再使用。
            g_constant: 重力常数，供论文口径浮力横向速度近似使用。
            gravity_yield_factor: 兼容旧脚本保留，核心中不再使用。
        """
        self.dt = dt
        self.nz = nz
        self.ny = ny
        self.total_t = total_t
        self.enable_d2dga = enable_d2dga
        self.d2dga_viscosity_ratio = d2dga_viscosity_ratio
        self.quality_penalty_scale = quality_penalty_scale
        self.channeling_penalty_weight = channeling_penalty_weight
        self.mixing_penalty_weight = mixing_penalty_weight
        self.instability_penalty_weight = instability_penalty_weight
        self.instability_decay_scale = instability_decay_scale
        self.save_interval: int = save_interval
        self.gel_growth_rate = gel_growth_rate
        self.gel_max_pa = gel_max_pa
        self.gel_break_threshold = gel_break_threshold
        self.T_surface = T_surface
        self.geothermal_gradient = geothermal_gradient
        self.T_ref = T_ref
        self.alpha_T = alpha_T
        self.enable_temperature_coupling = enable_temperature_coupling
        self.enable_mud_cake: bool = enable_mud_cake
        self.initial_mud_cake_mm: float = initial_mud_cake_mm
        self.k_erosion: float = k_erosion
        self.enable_turbulence: bool = enable_turbulence
        self.Re_critical: float = Re_critical
        self.turbulence_coefficient: float = turbulence_coefficient
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
        6. 构造核心求解使用的有效环空间隙（当前与原始环空间隙一致）

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
        - effective_b: 有效环空间隙数组；论文口径核心中等于 b
        """
        del mud_cake_thickness  # 兼容旧接口；论文口径核心不再修改有效环空间隙。
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
        geom["effective_b"] = geom["b"].copy()
        return geom

    def _update_effective_gap(self, geom: Dict[str, Array], mud_cake_thickness: Array) -> None:
        """兼容旧接口：论文口径核心中有效环空间隙始终等于原始环空间隙。"""
        del mud_cake_thickness
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
            gel_strength: 兼容旧接口保留，论文口径核心中不再使用。
            temperature_correction: 兼容旧接口保留，论文口径核心中不再使用。

        Returns:
            表观粘度数组
        """
        del gel_strength, temperature_correction
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
        return np.clip(mu, 1.0e-5, 3.0)

    def _smooth_dispersion(
        self,
        field: Array,
        *,
        axial: float = 0.018,
        azimuthal: float = 0.015,
    ) -> Array:
        """显式小系数拉普拉斯平滑，模拟D2DGA间隙尺度弥散。

        论文版采用显式二阶差分在轴向和方位角方向添加小系数弥散：
        - 轴向弥散系数通常 0.012–0.020；
        - 方位角弥散系数通常 0.012–0.018；
        - 边界处用一阶差分保持单侧稳定性。

        Args:
            field: 二维浓度场 (ny, nz)
            axial: 轴向弥散系数，默认 0.018
            azimuthal: 方位角弥散系数，默认 0.015

        Returns:
            平滑后的浓度场，裁剪到 [0, 1]
        """
        f = field.copy()
        # 轴向平滑（井深方向）：内部用二阶中心差分
        f[:, 1:-1] += axial * (field[:, 2:] - 2.0 * field[:, 1:-1] + field[:, :-2])
        # 方位角平滑（宽边→窄边方向）：内部用二阶中心差分
        f[1:-1, :] += azimuthal * (field[2:, :] - 2.0 * field[1:-1, :] + field[:-2, :])
        # 边界处理：用一阶差分避免越界
        f[0, :] += azimuthal * (field[1, :] - field[0, :])
        f[-1, :] += azimuthal * (field[-2, :] - field[-1, :])
        return np.clip(f, 0.0, 1.0)

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
        del gel_strength, temperature_correction
        mu = mud * self._apparent_viscosity(mud_fluid, gamma)
        if lead_fluid is not None:
            mu += lead * self._apparent_viscosity(lead_fluid, gamma)
        if tail_fluid is not None:
            mu += tail * self._apparent_viscosity(tail_fluid, gamma)
        if spacer_fluid is not None:
            mu += spacer * self._apparent_viscosity(spacer_fluid, gamma)
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
        """计算环空速度场（论文D2DGA口径）。

        采用 Zhang & Frigaard (2022) 的Hele-Shaw风格速度场：
        1. 计算局部混合流体的表观粘度与密度；
        2. 以 ``b²/μ`` 构造偏心通道主导局部流动度；
        3. 根据密度差（顶替液 vs 被顶替液）计算浮力稳定系数；
        4. 用浮力修正项调整宽边/窄边速度分配；
        5. 由截面排量约束得到轴向速度 ``w``。

        Returns:
            w: 井深方向速度（轴向速度）
            v: 方位角方向速度（横向速度）
            mu: 有效表观粘度场
            rho: 密度场
            mud: 钻井液分数场
            Re: 雷诺数场（仅作诊断，不参与湍流修正）
            mu_turbulent: 占位零场，保留旧结果对象兼容性
        """
        y = geom["y"]
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

        shear_rate = np.maximum(6.0 * np.abs(w_prev) / np.maximum(effective_b, 1.0e-5), 1.0e-6)
        D_h = 2.0 * geom["b"]
        rho_kg_m3 = rho * 1000.0
        Re = rho_kg_m3 * np.abs(w_prev) * D_h / np.maximum(mu, 1.0e-6)
        mu_turbulent = np.zeros_like(mu)

        # === 论文D2DGA口径速度场：偏心通道主导 + 浮力修正 ===
        # 基础流动度：偏心通道主导 (b/mean(b))^2 / mu
        b_mean = np.mean(b, axis=0, keepdims=True)
        base = (b / np.maximum(b_mean, 1.0e-12)) ** 2 / np.maximum(mu, 1.0e-6)

        # 浮力修正：基于顶替液与被顶替液的整体密度对比
        # density_contrast > 0 表示顶替液更重（水泥重 vs 泥浆轻），有助于窄边推进
        # density_contrast < 0 表示顶替液更轻，加剧宽边窜流
        if lead_fluid is not None and tail_fluid is not None:
            # 用领浆和尾浆的加权平均密度作为顶替液密度
            rho_disp = lead_fluid.density_kg_m3 * 0.67 + tail_fluid.density_kg_m3 * 0.33
        elif lead_fluid is not None:
            rho_disp = lead_fluid.density_kg_m3
        elif tail_fluid is not None:
            rho_disp = tail_fluid.density_kg_m3
        else:
            rho_disp = mud_fluid.density_kg_m3

        density_contrast = (rho_disp - mud_fluid.density_kg_m3) / mud_fluid.density_kg_m3
        # stable 系数：正值表示密度稳定（重顶替液推轻泥浆），负值表示不稳定
        stable = float(np.clip(8.0 * density_contrast, -0.35, 0.45))
        phi = geom["phi"][:, None]
        ebar = geom["e"][None, :]
        # buoyancy_shape: 在窄边(phi=1)处 = 1 + stable*ebar，在宽边(phi=0)处 = 1 - stable*ebar
        buoyancy_shape = 1.0 + stable * ebar * (2.0 * phi - 1.0)
        pref = np.maximum(base * buoyancy_shape, 1.0e-8)

        # 由截面排量约束得到轴向速度 w
        dy = np.gradient(geom["y"])[:, None]
        area_weight = np.sum(pref * b * dy * 2.0, axis=0, keepdims=True)
        w = q_half * pref / np.maximum(area_weight, 1.0e-12)

        # 横向速度 v 由连续性方程求解（简化处理，论文版未显式计算 v）
        ds = geom["s"][1] - geom["s"][0]
        bw = b * w
        dbw_ds = np.gradient(bw, ds, axis=1)
        bv = np.zeros_like(w)
        for i in range(1, len(y)):
            bv[i, :] = bv[i - 1, :] - 0.5 * (dbw_ds[i, :] + dbw_ds[i - 1, :]) * dy[i - 1, 0]
        bv -= (y[:, None] / y[-1]) * bv[-1, :]
        v = bv / np.maximum(b, 1.0e-8)

        return w, v, mu, rho, mud, Re, mu_turbulent

    def _depth_profiles(self, geom: Dict[str, Array], lead: Array, tail: Array, spacer: Array, wall: Array) -> pd.DataFrame:
        """计算深度方向的平均剖面数据。"""
        cement = np.clip(lead + tail, 0.0, 1.0)
        del wall
        eff = cement
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
        """运行论文口径的环空二维顶替求解。"""

        mud_fluid, lead_fluid, tail_fluid, spacer_fluid = self._pick_fluids(fluids)
        geom = self._build_geom(well_spec)
        lead = np.zeros((self.ny, self.nz), dtype=float)
        tail = np.zeros((self.ny, self.nz), dtype=float)
        spacer = np.zeros((self.ny, self.nz), dtype=float)
        # 为兼容旧结果结构仍保留 wall 字段，但论文口径核心不再做壁面泥饼清洗，故恒为零。
        wall = np.zeros((self.ny, self.nz), dtype=float)
        gel_strength = np.zeros((self.ny, self.nz), dtype=float)
        viscosity_correction = np.ones((self.ny, self.nz), dtype=float)
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
        cumulative_lead_in_m3 = 0.0
        cumulative_tail_in_m3 = 0.0
        cumulative_spacer_in_m3 = 0.0

        mud_cake_thickness: Array = np.zeros((self.ny, self.nz), dtype=float)
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
                # === 正常泵注阶段：仅执行论文口径核心平流 + D2DGA 通量修正 ===
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
                w_prev = w

                cement = np.clip(lead + tail, 0.0, 1.0)
                if self.enable_d2dga:
                    f_amp = d2dga_flux_amplification(cement, self.d2dga_viscosity_ratio)
                else:
                    f_amp = 1.0

                w_d2dga = w * f_amp
                v_d2dga = v * f_amp
                ysrc = ygrid - v_d2dga * self.dt
                ssrc = sgrid - w_d2dga * self.dt
                lead_adv = _bilinear_interp(lead, ysrc, ssrc, geom, inlet_lead_fraction)
                tail_adv = _bilinear_interp(tail, ysrc, ssrc, geom, inlet_tail_fraction)
                spacer_adv = _bilinear_interp(spacer, ysrc, ssrc, geom, inlet_spacer_fraction)
                lead = np.clip(lead_adv, 0.0, 1.0)
                tail = np.clip(tail_adv, 0.0, 1.0)
                spacer = np.clip(spacer_adv, 0.0, 1.0)
                # 数值扩散可能使显式两相之和略超1；按比例压回可行域，保持泥浆分数非负。
                tracked_total = lead + tail + spacer
                overfilled = tracked_total > 1.0
                lead[overfilled] /= tracked_total[overfilled]
                tail[overfilled] /= tracked_total[overfilled]
                spacer[overfilled] /= tracked_total[overfilled]

                # D2DGA间隙尺度弥散：在低浓度前锋更强，模拟间隙尺度分散效应
                lead = self._smooth_dispersion(lead, axial=0.018, azimuthal=0.015)
                tail = self._smooth_dispersion(tail, axial=0.018, azimuthal=0.015)
                spacer = self._smooth_dispersion(spacer, axial=0.012, azimuthal=0.012)

                # D2DGA 通量放大会改变前锋形态，但不应让各相总量超过累计入环空体积。
                # 这里按入口累计体积对领浆、尾浆和前置液/隔离液分别做体积上限约束。
                cumulative_lead_in_m3 += inlet_state.flow_rate_m3_s * inlet_lead_fraction * self.dt
                cumulative_tail_in_m3 += inlet_state.flow_rate_m3_s * inlet_tail_fraction * self.dt
                cumulative_spacer_in_m3 += inlet_state.flow_rate_m3_s * inlet_spacer_fraction * self.dt
                lead = _limit_phase_volume(lead, geom, cumulative_lead_in_m3)
                tail = _limit_phase_volume(tail, geom, cumulative_tail_in_m3)
                spacer = _limit_phase_volume(spacer, geom, cumulative_spacer_in_m3)

            else:
                # === 泵停阶段：冻结浓度场，仅记录指标 ===
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
                # 泵停后保持上一时刻浓度场，不再引入停泵滑移或额外壁面过程。

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
                reynolds_snapshots.append(Re.copy())
                turbulent_viscosity_snapshots.append(mu_turbulent.copy())
                snapshot_times.append(current_time_s)

            cement = np.clip(lead + tail, 0.0, 1.0)
            eff = cement
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

            rows.append(
                [
                    current_time_s,
                    current_time_s / 60.0,
                    inlet_state.stage_name,
                    bulk_fill,
                    effective_efficiency,
                    target_efficiency,
                    cbl_efficiency,
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
                "CBL评价井段水力有效顶替效率": float(final["cbl_eval_interval_efficiency"]),
                "目标层段水力有效顶替效率": float(final["target_interval_efficiency"]),
                "最终水泥浆占据率": float(final["bulk_cement_fill"]),
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
                "solver 核心仅保留论文口径的几何、流动度、D2DGA 通量修正与求解域效率输出。",
                "泥饼、温度、凝胶强度、湍流修正与 CBL 质量惩罚不再影响求解结果；相关字段仅作兼容占位。",
            ),
            lead_field=lead,
            tail_field=tail,
        )
