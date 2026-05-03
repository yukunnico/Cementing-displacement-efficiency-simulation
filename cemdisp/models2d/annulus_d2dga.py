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
- total_t: 总模拟时间（秒），默认6600秒（110分钟）
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
    spacer_snapshots: Tuple[Array, ...] = field(default_factory=tuple)
    wall_snapshots: Tuple[Array, ...] = field(default_factory=tuple)
    snapshot_times_s: Tuple[float, ...] = field(default_factory=tuple)
    notes: Tuple[str, ...] = field(default_factory=tuple)


class AnnulusD2DGASolver:
    """环空二维D2DGA顶替求解器。

    实现偏心环空中水泥浆、前置/隔离液和钻井液的多相体积分数模拟。

    主要物理过程：
    1. 平流输运：水泥浆与前置/隔离液随速度场向下游运移
    2. 弥散效应：横向(day)和垂向(dax)扩散，与偏心度和流体性质相关
    3. 壁面清除：泥饼在剪切作用下被水泥浆与前置/隔离液逐渐清除
    4. 浮力效应：密度差导致窄边优先顶替和潜在失稳风险

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
        total_t: float = 6600.0,
        quality_penalty_scale: float = 0.099,
        channeling_penalty_weight: float = 0.55,
        mixing_penalty_weight: float = 0.35,
        instability_penalty_weight: float = 0.25,
        instability_decay_scale: float = 5.0,
        save_interval: int = 60,
    ) -> None:
        """初始化环空二维求解器参数。

        Args:
            dt: 时间步长（秒），默认4秒
            nz: 井深方向网格数，默认140
            ny: 方位角方向网格数，默认40
            alpha_clean: 泥饼清除系数，默认0.085
            total_t: 总模拟时间（秒），默认6600秒（110分钟）
            quality_penalty_scale: 质量惩罚因子，默认0.099
            channeling_penalty_weight: 窜槽风险权重，默认0.55
            mixing_penalty_weight: 混浆风险权重，默认0.35
            instability_penalty_weight: 失稳风险权重，默认0.25
            instability_decay_scale: 失稳指数衰减标度，默认5.0
            save_interval: 二维场快照保存步长，默认每60个时间步保存一次
        """
        self.dt = dt
        self.nz = nz
        self.ny = ny
        self.alpha_clean = alpha_clean
        self.total_t = total_t
        self.quality_penalty_scale = quality_penalty_scale
        self.channeling_penalty_weight = channeling_penalty_weight
        self.mixing_penalty_weight = mixing_penalty_weight
        self.instability_penalty_weight = instability_penalty_weight
        self.instability_decay_scale = instability_decay_scale
        self.save_interval: int = save_interval

    def _build_geom(self, well_spec: WellSpec) -> Dict[str, Array]:
        """根据井筒规格构建环空二维网格几何参数。

        构建步骤：
        1. 创建井深方向网格s和深度坐标md
        2. 从井径剖面插值得到井径、井斜、偏心度数据
        3. 构建方位角方向网格y和归一化方位角phi
        4. 计算半间隙h和局部环空间隙b
        5. 校正体积使其等于物理环空体积

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
        return geom

    def _physical_annular_volume(self, well_spec: WellSpec) -> float:
        """计算井段物理环空体积（用于体积校正）。"""
        cal_md, cal_hole = _profile_to_arrays(well_spec.hole_diameter_profile)
        od = np.full_like(cal_md, float(well_spec.liner_od_mm or 0.0), dtype=float)
        area = np.pi * ((cal_hole / 1000.0) ** 2 - (od / 1000.0) ** 2) / 4.0
        return float(np.trapezoid(area, x=cal_md))

    def _pick_fluids(self, fluids: Tuple[FluidSpec, ...]) -> Tuple[FluidSpec, FluidSpec, FluidSpec | None]:
        """从流体列表中选取钻井液、水泥浆和可选前置/隔离液。"""
        mud = next((fluid for fluid in fluids if fluid.role == FluidRole.MUD), None)
        cement = next((fluid for fluid in fluids if fluid.role in {FluidRole.LEAD, FluidRole.TAIL}), None)
        spacer = next((fluid for fluid in fluids if fluid.role in {FluidRole.WASH, FluidRole.SPACER}), None)
        if mud is None or cement is None:
            raise ValueError("需要钻井液和至少一个水泥浆流体")
        return mud, cement, spacer

    def _apparent_viscosity(self, fluid: FluidSpec, gamma: Array) -> Array:
        """根据流变模型计算流体的表观粘度。

        Args:
            fluid: 流体规格
            gamma: 剪切速率数组

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
        return np.clip(mu, 1.0e-5, 3.0)

    def _compute_props(
        self,
        cement: Array,
        spacer: Array,
        w_prev: Array,
        geom: Dict[str, Array],
        mud_fluid: FluidSpec,
        cement_fluid: FluidSpec,
        spacer_fluid: FluidSpec | None,
    ) -> Tuple[Array, Array, Array]:
        """计算三相混合物系的表观粘度、密度和钻井液分数。"""
        # 三相体积分数闭合：显式跟踪水泥浆和前置/隔离液，钻井液由守恒关系反算。
        mud = np.clip(1.0 - cement - spacer, 0.0, 1.0)
        gamma = np.maximum(6.0 * np.abs(w_prev) / np.maximum(geom["b"], 1.0e-5), 1.0e-6)
        mu = mud * self._apparent_viscosity(mud_fluid, gamma)
        mu += cement * self._apparent_viscosity(cement_fluid, gamma)
        if spacer_fluid is not None:
            mu += spacer * self._apparent_viscosity(spacer_fluid, gamma)
        rho = mud * (mud_fluid.density_kg_m3 / 1000.0)
        rho += cement * (cement_fluid.density_kg_m3 / 1000.0)
        if spacer_fluid is not None:
            rho += spacer * (spacer_fluid.density_kg_m3 / 1000.0)
        return mu, rho, mud

    def _compute_velocity(
        self,
        cement: Array,
        spacer: Array,
        geom: Dict[str, Array],
        q_m3s: float,
        w_prev: Array,
        mud_fluid: FluidSpec,
        cement_fluid: FluidSpec,
        spacer_fluid: FluidSpec | None,
    ) -> Tuple[Array, Array, Array, Array, Array]:
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
            mu: 表观粘度场
            rho: 密度场
            mud: 钻井液分数场
        """
        y = geom["y"]
        b = geom["b"]
        q_half = q_m3s / 2.0
        mu, rho, mud = self._compute_props(cement, spacer, w_prev, geom, mud_fluid, cement_fluid, spacer_fluid)
        mud_density_gcc = mud_fluid.density_kg_m3 / 1000.0

        buoy = 1.0 + 0.60 * np.maximum(rho - mud_density_gcc, 0.0) + 0.08 * np.sin(np.deg2rad(geom["inc_deg"]))[None, :]
        m = (b**3 / np.maximum(mu, 1.0e-6)) * buoy * (0.90 + 0.25 * geom["standoff"][None, :])

        phi = geom["phi"][:, None]
        narrow_boost = 1.0 - 0.30 * cement * np.maximum(rho - mud_density_gcc, 0.0) * np.cos(np.pi * phi)
        m *= np.clip(narrow_boost, 0.75, 1.70)

        int_m = np.trapezoid(b * m, x=y, axis=0)
        w = (q_half / np.maximum(int_m, 1.0e-12))[None, :] * m

        ds = geom["s"][1] - geom["s"][0]
        dy = y[1] - y[0]
        bw = b * w
        dbw_ds = np.gradient(bw, ds, axis=1)
        bv = np.zeros_like(w)
        for i in range(1, len(y)):
            bv[i, :] = bv[i - 1, :] - 0.5 * (dbw_ds[i, :] + dbw_ds[i - 1, :]) * dy
        bv -= (y[:, None] / y[-1]) * bv[-1, :]
        v = bv / np.maximum(b, 1.0e-8)
        return w, v, mu, rho, mud

    def _depth_profiles(self, geom: Dict[str, Array], cement: Array, spacer: Array, wall: Array) -> pd.DataFrame:
        """计算深度方向的平均剖面数据。"""
        eff = cement * (1.0 - wall)
        mud = np.clip(1.0 - cement - spacer, 0.0, 1.0)
        return pd.DataFrame(
            {
                "井深_m": geom["md"],
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

        mud_fluid, cement_fluid, spacer_fluid = self._pick_fluids(fluids)
        geom = self._build_geom(well_spec)
        cement = np.zeros((self.ny, self.nz), dtype=float)
        spacer = np.zeros((self.ny, self.nz), dtype=float)
        wall = np.ones((self.ny, self.nz), dtype=float)
        w_prev = np.full((self.ny, self.nz), 0.45, dtype=float)
        half_volume = _trapez2d(geom["b"], geom)

        target_mask = _window_mask(well_spec, geom["md"], "target")
        cbl_mask = _window_mask(well_spec, geom["md"], "cbl")
        ygrid, sgrid = np.meshgrid(geom["y"], geom["s"], indexing="ij")
        rows: list[list[float | str]] = []
        mud_density_gcc = mud_fluid.density_kg_m3 / 1000.0
        cement_snapshots: list[Array] = []
        spacer_snapshots: list[Array] = []
        wall_snapshots: list[Array] = []
        snapshot_times: list[float] = []

        final_step_index = int(self.total_t / self.dt)
        for step_index in range(final_step_index + 1):
            current_time_s = step_index * self.dt
            inlet_state = inlet_state_provider(current_time_s)
            inlet_cement_fraction = _phase_fraction(inlet_state, "cement")
            inlet_spacer_fraction = _phase_fraction(inlet_state, "spacer")

            # 泵停判断：排量低于阈值时认为泵已停止。
            # 泵停后水泥场冻结——不再平流、扩散或壁面清除，
            # 因为水泥静凝胶强度在短时间（<2h）内足以抵抗浮力滑塌。
            pump_active = inlet_state.flow_rate_m3_s > 1.0e-9

            if pump_active:
                # === 正常泵注阶段：执行全部物理过程 ===
                w, v, mu, rho, mud = self._compute_velocity(
                    cement,
                    spacer,
                    geom,
                    inlet_state.flow_rate_m3_s,
                    w_prev,
                    mud_fluid,
                    cement_fluid,
                    spacer_fluid,
                )
                w_prev = w

                buoy_num = np.maximum(rho - mud_density_gcc, 0.0)
                disp_ratio = np.clip(mud / (cement + 1.0e-4), 0.2, 8.0)
                dax = 0.0040 * np.abs(w) * geom["b"] * (0.20 + 1.20 * cement * (1.0 - cement)) * (1.0 + 0.50 * geom["e"][None, :]) * (1.0 + 0.10 * np.maximum(disp_ratio - 1.0, 0.0)) / (1.0 + 0.95 * buoy_num)
                day = 0.16 * dax * (0.45 + 0.55 * geom["e"][None, :])

                v_eff = v + 0.007 * np.maximum(rho - mud_density_gcc, 0.0) * cement * np.sin(np.pi * geom["phi"])[:, None]
                ysrc = ygrid - v_eff * self.dt
                ssrc = sgrid - w * self.dt
                # 水泥浆与前置/隔离液使用同一速度场输运；只改变被输运的浓度变量。
                cement_adv = _bilinear_interp(cement, ysrc, ssrc, geom, inlet_cement_fraction)
                spacer_adv = _bilinear_interp(spacer, ysrc, ssrc, geom, inlet_spacer_fraction)

                dy = geom["y"][1] - geom["y"][0]
                ds = geom["s"][1] - geom["s"][0]
                lap_s = np.zeros_like(cement_adv)
                lap_y = np.zeros_like(cement_adv)
                spacer_lap_s = np.zeros_like(spacer_adv)
                spacer_lap_y = np.zeros_like(spacer_adv)
                lap_s[:, 1:-1] = (cement_adv[:, 2:] - 2.0 * cement_adv[:, 1:-1] + cement_adv[:, :-2]) / ds**2
                lap_s[:, 0] = (cement_adv[:, 1] - cement_adv[:, 0]) / ds**2
                lap_s[:, -1] = (cement_adv[:, -2] - cement_adv[:, -1]) / ds**2
                lap_y[1:-1, :] = (cement_adv[2:, :] - 2.0 * cement_adv[1:-1, :] + cement_adv[:-2, :]) / dy**2
                lap_y[0, :] = (cement_adv[1, :] - cement_adv[0, :]) / dy**2
                lap_y[-1, :] = (cement_adv[-2, :] - cement_adv[-1, :]) / dy**2
                spacer_lap_s[:, 1:-1] = (spacer_adv[:, 2:] - 2.0 * spacer_adv[:, 1:-1] + spacer_adv[:, :-2]) / ds**2
                spacer_lap_s[:, 0] = (spacer_adv[:, 1] - spacer_adv[:, 0]) / ds**2
                spacer_lap_s[:, -1] = (spacer_adv[:, -2] - spacer_adv[:, -1]) / ds**2
                spacer_lap_y[1:-1, :] = (spacer_adv[2:, :] - 2.0 * spacer_adv[1:-1, :] + spacer_adv[:-2, :]) / dy**2
                spacer_lap_y[0, :] = (spacer_adv[1, :] - spacer_adv[0, :]) / dy**2
                spacer_lap_y[-1, :] = (spacer_adv[-2, :] - spacer_adv[-1, :]) / dy**2
                cement = np.clip(cement_adv + self.dt * (dax * lap_s + day * lap_y), 0.0, 1.0)
                spacer = np.clip(spacer_adv + self.dt * (dax * spacer_lap_s + day * spacer_lap_y), 0.0, 1.0)
                # 数值扩散可能使显式两相之和略超1；按比例压回可行域，保持泥浆分数非负。
                tracked_total = cement + spacer
                overfilled = tracked_total > 1.0
                cement[overfilled] /= tracked_total[overfilled]
                spacer[overfilled] /= tracked_total[overfilled]

                # 清洗能力由水泥浆主导，隔离液/冲洗液也参与泥饼清除但权重较低。
                cleaner = 1.10 * (cement + 0.8 * spacer)
                shear = np.abs(w) / np.maximum(geom["b"], 1.0e-5)
                kclean = self.alpha_clean * cleaner * (0.45 + np.sqrt(np.maximum(shear, 1.0e-6))) * (0.85 + 0.35 * geom["standoff"][None, :]) * (1.0 + 0.30 * buoy_num)
                wall *= np.exp(-kclean * self.dt / 150.0)
                wall = np.clip(wall, 0.0, 1.0)

            else:
                # === 泵停阶段：冻结水泥场与壁面场，仅记录指标 ===
                # 使用上一步的 mu/rho/mud 计算 mobility 等指标（用于风险指数），
                # 但不更新 cement 或 wall，避免浮力流"蒸发"已就位的水泥。
                w, v, mu, rho, mud = self._compute_velocity(
                    cement, spacer, geom, 0.0, w_prev, mud_fluid, cement_fluid, spacer_fluid
                )
                # 注意：w_prev 保持泵停前最后的值，不再更新，
                # 因泵停后速度场已冻结。

            # 在物理场更新后、指标计算前保存快照，确保快照与本步指标使用同一状态。
            # 使用 copy() 固化二维场，避免后续时间步原地更新影响已保存结果。
            if step_index % self.save_interval == 0 or step_index == int(self.total_t / self.dt):
                cement_snapshots.append(cement.copy())
                spacer_snapshots.append(spacer.copy())
                wall_snapshots.append(wall.copy())
                snapshot_times.append(current_time_s)

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
        depth_profiles = self._depth_profiles(geom, cement, spacer, wall)
        final = metrics.iloc[-1]
        summary: Dict[str, object] = {
            "模型名称": "Hu102尾管段首版环空顶替模型",
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
            spacer_snapshots=tuple(spacer_snapshots),
            wall_snapshots=tuple(wall_snapshots),
            snapshot_times_s=tuple(snapshot_times),
            notes=(
                "当前显式跟踪水泥浆相（领浆/尾浆）与前置液/隔离液相，钻井液由体积分数闭合反算。",
                "效率与质量指标仍按水泥浆相计算，前置液/隔离液仅影响流变混合、输运占据和壁面清洗。",
            ),
        )
