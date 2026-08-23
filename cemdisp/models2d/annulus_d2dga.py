"""
多井通用环空二维 D2DGA 求解器。

本模块实现面向 ``cemdisp`` 主链路的偏心环空二维 D2DGA 核心，口径尽量
收敛到 Zhang & Frigaard (2022) 所强调的环空层流顶替主过程：

1. 偏心窄环空几何展开；
2. 基于局部流动度的轴向/方位角平均速度场；
3. D2DGA 通量放大修正，近似捕捉间隙尺度分散；
4. 仅输出求解域内的顶替效率与浓度场。

出口边界条件：
- 开放出口（open_outlet=True，默认）：允许水泥浆流出求解域到重叠段，适用于只模拟裸眼段；
- 封闭出口（open_outlet=False）：按累计入环空体积限制场量，适用于模拟整个环空。

注意：
- 本模块不再把泥饼、温度、凝胶强度、湍流修正等工程扩展项作为核心求解的一部分；
- 为兼顾下游脚本兼容性，旧参数与旧快照字段仍保留接口或占位输出，但不再影响求解结果。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, Tuple

import numpy as np
from numpy.typing import NDArray
import pandas as pd

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.well_spec import DepthValuePoint, WellSpec
from cemdisp.diagnostics.displacement_metrics import _narrow_quarter_efficiency
from cemdisp.models2d.boundary_bridge import AnnulusInletState
from cemdisp.models2d.d2dga_flux import (
    d2dga_buoyancy_flux,
    d2dga_dispersion_I1,
    d2dga_dispersion_I2,
    d2dga_flux_amplification,
)


Array = NDArray[np.float64]


def _profile_to_arrays(points: Tuple[DepthValuePoint, ...]) -> Tuple[Array, Array]:
    """将剖面数据点列表转换为NumPy数组。"""
    depths = np.array([float(point.depth_md_m) for point in points], dtype=float)
    values = np.array([float(point.value) for point in points], dtype=float)
    return depths, values


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


def _evaluation_window_efficiencies(well_spec: "WellSpec", geom: Dict[str, Array], cement: Array) -> Dict[str, dict]:
    """对每个 EvaluationWindow 做 b 加权 2D 积分，返回 {窗名: {window_type, eta_E, eta_N}}。"""
    out: Dict[str, dict] = {}
    md = geom["md"]            # (nz,)，md = bottom - s
    b_full = geom["b"]         # (ny,nz)
    s_full = geom["s"]         # (nz,)
    for w in well_spec.evaluation_windows:
        mask = (md >= w.top_md_m) & (md <= w.bottom_md_m)
        if not bool(mask.any()):
            continue
        b_win = b_full[:, mask]
        c_win = cement[:, mask]
        geom_win = {**geom, "b": b_win, "s": s_full[mask]}
        denom = _trapez2d(b_win, geom_win)
        eta_e = float(_trapez2d(b_win * c_win, geom_win) / max(denom, 1e-12))
        eta_n = float(_narrow_quarter_efficiency(c_win, geom_win))
        out[w.name] = {"window_type": w.window_type, "eta_E": eta_e, "eta_N": eta_n}
    return out


def _low_tail_indicators(geom: Dict[str, Array], cement: Array, ny: int) -> Dict[str, float]:
    """低尾指标：standoff<0.5 段占比 + 窄边（最后 ny//4 行）cement<0.05 域体积占比。"""
    b_full = geom["b"]
    so_frac = float(np.mean(geom["standoff"] < 0.5))
    n_q = max(1, ny // 4)
    b_q = b_full[-n_q:, :]
    low = (cement[-n_q:, :] < 0.05).astype(float)
    # 行切片后必须同步切 y，否则 _trapez2d 对 axis=0 的 x=geom["y"] 长度与行数不匹配
    geom_q = {**geom, "b": b_q, "y": geom["y"][-n_q:]}
    tail_frac = float(_trapez2d(b_q * low, geom_q) / max(_trapez2d(b_q, geom_q), 1e-12))
    return {"standoff低于0.5段占比": so_frac, "窄边效率低于0.05域占比": tail_frac}


def _limit_phase_volume(field: Array, geom: Dict[str, Array], target_volume_m3: float, open_outlet: bool = False) -> Array:
    """按累计入环空体积限制场量，避免数值扩散凭空放大相体积。

    Args:
        field: 浓度场（ny×nz 数组）
        geom: 几何参数字典
        target_volume_m3: 目标体积（累计入环空体积，立方米）
        open_outlet: 是否开放出口边界。True 时不限制体积，允许水泥浆流出到重叠段。

    Returns:
        限制后的浓度场，裁剪到 [0, 1]
    """
    if open_outlet:
        return np.clip(field, 0.0, 1.0)
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
    snapshot_times_s: Tuple[float, ...] = field(default_factory=tuple)
    notes: Tuple[str, ...] = field(default_factory=tuple)
    lead_field: Array = field(default_factory=lambda: np.empty((0, 0), dtype=float))
    tail_field: Array = field(default_factory=lambda: np.empty((0, 0), dtype=float))
    # T1-6: FLUSHER 独立浓度场（被动平流相，不参与 D2DGA 闭包）
    flusher_field: Array | None = field(default=None)
    flusher_snapshots: Tuple[Array, ...] = field(default_factory=tuple)


class AnnulusD2DGASolver:
    """环空二维D2DGA顶替求解器。

    实现偏心环空中水泥浆、前置/隔离液和钻井液的多相体积分数模拟。

    主要物理过程：
    1. 平流输运：水泥浆与前置/隔离液随平均速度场向下游运移；
    2. D2DGA 分散：基于 Zhang & Frigaard (2022) 的通量放大修正，近似捕捉间隙尺度分散；
    3. 浮力相关横向再分布：以 Hele-Shaw 风格的密度差横向速度近似，保持密度差影响仍在核心层；
    4. 出口边界：支持开放边界（允许水泥浆流出到重叠段）或封闭边界（限制体积）；
    5. 仅输出求解域指标，不在 solver 内叠加现场质量惩罚或 CBL 校准。

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
        # 开放出口边界（默认，适用于只模拟裸眼段）
        solver = AnnulusD2DGASolver(open_outlet=True)

        # 封闭出口边界（适用于模拟整个环空）
        solver = AnnulusD2DGASolver(open_outlet=False)

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
        enable_d2dga_auto_m: bool = True,
        enable_d2dga_i3_flux: bool = True,
        enable_local_i3: bool = False,
        enable_true_buoyancy: bool = True,
        instability_decay_scale: float = 5.0,
        save_interval: int = 60,
        yield_regularization_M: float = 100.0,
        enable_regime_split: bool = False,
        regime_relax_alpha: float = 0.5,
        regime_max_iter: int = 24,
        regime_tol_rel: float = 1e-3,
        regime_re_turb_ratio: float = 1.8,
        open_outlet: bool = True,
        alpha_cfl: float = 0.5,
        enable_cfl_adaptive: bool = True,
        cfl_number: float = 0.5,
        dt_min: float = 0.1,
        c_min: float = 0.05,
        e_clip_max: float = 0.55,
        enable_yield_gate: bool = False,
        yield_gate_f_safety: float = 1.15,
        yield_gate_c_min_residual: float = 0.01,
        dispersion_axial: float = 0.018,
        dispersion_azimuthal: float = 0.015,
        dispersion_dt_ref: float = 4.0,
        dispersion_dt_scale: float = 1.0,
    ) -> None:
        """初始化环空二维求解器参数。

        Args:
            dt: 时间步长（秒），默认4秒
            nz: 井深方向网格数，默认140
            ny: 方位角方向网格数，默认40
            total_t: 总模拟时间（秒），默认12000秒（200分钟）
            enable_d2dga: 是否启用D2DGA通量修正（Zhang & Frigaard 2022），默认开启
            d2dga_viscosity_ratio: D2DGA粘度比 m = η_displaced/η_displacing，默认1.0
            enable_d2dga_auto_m: 是否按局部流体物性自动计算黏度比 m 场（R1），默认 True。
                True: m = μ_displaced/μ_displacing 按浓度场每步计算（改进版）；
                False: 退化为 d2dga_viscosity_ratio 构造常数（旧论文 R0 状态）。
            enable_d2dga_i3_flux: 是否启用 D2DGA 浮力弥散通量 I3（R2，式4.25第二项），默认 True。
            enable_local_i3: I3 通量局部化开关，默认 False（不改变既有行为）。
                False: eta2 用 cement 表观粘度场均值、Δρ 用全场均值（基线逐位复现）；
                True: eta2 透传水泥相黏度场 _eta2、Δρ 用 (rho-mud_density_gcc)*1000 局部场。
            enable_true_buoyancy: 是否用真浮力体力替换 buoyancy_shape 代理（R3，式2.5b），默认 True。
                False: 保留 buoyancy_shape 代理（旧论文 R0/R1/R2 状态）。
            instability_decay_scale: 后验失稳指数缩放，默认5.0。
            save_interval: 二维场快照保存步长，默认每60个时间步保存一次
            yield_regularization_M: Papanastasiou正则化参数，控制屈服应力在低剪切区的平滑过渡，默认100.0
            enable_regime_split: M2 局部流态修正固定点迭代开关（Maleki & Frigaard 2017 式58-66），默认 False。
                False: 原 b²/μ 流动度分配逐字节复现基线；True: 迭代 Re_p/阻力权重 R 修正 pref 分配。
                黏度场保持 w_prev 一步滞后（既有约定），仅流态权重 R 参与迭代。
            regime_relax_alpha: 固定点迭代欠松弛系数 α，默认 0.5。
            regime_max_iter: 固定点迭代最大迭代次数，默认 24。
            regime_tol_rel: w 相对变化收敛容差，默认 1e-3。
            regime_re_turb_ratio: 湍流起始 Re 相对 re_crit 的倍数 re_turb = re_crit·ratio，默认 1.8。
            open_outlet: 是否开放出口边界（允许水泥浆流出到重叠段），默认True。
                True: 开放出口，不限制体积，适用于只模拟裸眼段；
                False: 封闭出口，按累计入环空体积限制场量，适用于模拟整个环空。
            alpha_cfl: CFL 裁剪系数，默认 0.5。控制单步 I3 通量散度裁剪上限
                |div_q|·dt ≤ alpha_cfl·min(ds)，其中 ds 取轴向网格最小间距。
            enable_cfl_adaptive: 是否启用全局 CFL 自适应时间步，默认 True。
                True: 每步按 CFL 条件动态调整 dt_step；
                False: 固定 dt 复现基线。
            cfl_number: 全局 CFL 数（半拉格朗日保守估计），默认 0.5。
                独立于 T1-2 alpha_cfl（局部 I3 通量裁剪）。
            dt_min: 自适应时间步下限（秒），默认 0.1。防 CFL 过小步数爆炸。
            c_min: 壁面静止层浓度阈值（Bararpour 2025 式 2.35-2.41），默认 0.05。
                局部水泥浓度 c < c_min 时该处壁面层泥浆滞留不流动（wall=1）。
                屈服门槛关闭（enable_yield_gate=False）时作为 OFF 路径兜底。
            e_clip_max: M4 偏心度 e 硬截断上限，默认 0.55（逐位复现基线）。
                由 e = clip(1-standoff, 0.05, e_clip_max) 构造几何；
                生产跑道（Task 13 重跑阶段）显式设 0.90 放宽截断。
                体积校正（_build_geom 末尾 scale）每个 run 重算，half_volume 守恒。
            enable_yield_gate: M3 屈服门槛开关，默认 False（不改变既有行为）。
                True 时 pump 分支用 _yield_gate_wall 重建壁面冻结层（Task 9 接线）。
            yield_gate_f_safety: 屈服门槛安全系数 f，默认 1.15。
                immobile 判定：外推壁剪 τw_extrap ≤ f·τy。
            yield_gate_c_min_residual: 残泥下限浓度，默认 0.01。
                水泥已到（cement_ever>0）但局部浓度 < 该值 -> 视为残泥冻结。
            dispersion_axial: D2DGA 间隙尺度弥散轴向系数（每 dt_ref 秒），默认 0.018。
            dispersion_azimuthal: D2DGA 间隙尺度弥散方位角系数（每 dt_ref 秒），默认 0.015。
            dispersion_dt_ref: 弥散系数的名义/参考时间步（秒），默认 4.0。
                系数 fa44ace 引入时即 dt=4.0，故 dt_ref=4.0 使固定 dt 模式逐位复现基线。
            dispersion_dt_scale: 弥散系数按 dt 归一开关，默认 1.0。
                =1.0 时固定 dt 模式（dt_step==dt_ref）逐位复现基线；
                CFL 模式下按 dt_step/dt_ref 同比缩放，使每物理秒弥散恒定。
        """
        self.dt = dt
        self.nz = nz
        self.ny = ny
        self.total_t = total_t
        self.enable_d2dga = enable_d2dga
        self.d2dga_viscosity_ratio = d2dga_viscosity_ratio
        self.instability_decay_scale = instability_decay_scale
        self.save_interval: int = save_interval
        self.yield_regularization_M: float = yield_regularization_M
        self.enable_regime_split: bool = enable_regime_split
        self.regime_relax_alpha: float = regime_relax_alpha
        self.regime_max_iter: int = regime_max_iter
        self.regime_tol_rel: float = regime_tol_rel
        self.regime_re_turb_ratio: float = regime_re_turb_ratio
        self.enable_d2dga_auto_m: bool = enable_d2dga_auto_m
        self.enable_d2dga_i3_flux: bool = enable_d2dga_i3_flux
        self.enable_local_i3: bool = enable_local_i3
        self.enable_true_buoyancy: bool = enable_true_buoyancy
        self.open_outlet: bool = open_outlet
        self.alpha_cfl: float = alpha_cfl
        self.enable_cfl_adaptive: bool = enable_cfl_adaptive
        self.cfl_number: float = cfl_number
        self.dt_min: float = dt_min
        self.c_min: float = c_min
        self.e_clip_max: float = e_clip_max
        self.enable_yield_gate: bool = enable_yield_gate
        self.yield_gate_f_safety: float = yield_gate_f_safety
        self.yield_gate_c_min_residual: float = yield_gate_c_min_residual
        self.dispersion_axial = dispersion_axial
        self.dispersion_azimuthal = dispersion_azimuthal
        self.dispersion_dt_ref = dispersion_dt_ref
        self.dispersion_dt_scale = dispersion_dt_scale

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
        if well_spec.liner_od_profile:
            od_md, od_values = _profile_to_arrays(well_spec.liner_od_profile)
            od_mm = np.interp(md, od_md, od_values)
        else:
            od_mm = np.full_like(md, float(well_spec.liner_od_mm or 0.0), dtype=float)

        e = np.clip(1.0 - standoff, 0.05, self.e_clip_max)
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

    def _physical_annular_volume(self, well_spec: WellSpec) -> float:
        """计算井段物理环空体积（用于体积校正）。"""
        cal_md, cal_hole = _profile_to_arrays(well_spec.hole_diameter_profile)
        if well_spec.liner_od_profile:
            od_md, od_values = _profile_to_arrays(well_spec.liner_od_profile)
            od = np.interp(cal_md, od_md, od_values)
        else:
            od = np.full_like(cal_md, float(well_spec.liner_od_mm or 0.0), dtype=float)
        area = np.pi * ((cal_hole / 1000.0) ** 2 - (od / 1000.0) ** 2) / 4.0
        return float(np.trapezoid(area, x=cal_md))

    def _pick_fluids(
        self,
        fluids: Tuple[FluidSpec, ...],
    ) -> Tuple[FluidSpec, FluidSpec | None, FluidSpec | None, FluidSpec | None, FluidSpec | None]:
        """从流体列表中选取钻井液、领浆、尾浆、可选前置/隔离液和可选冲洗液。

        自 T1-6 起返回 5 元组 (mud, lead, tail, spacer, flusher)，
        其中 lead/tail/spacer/flusher 均可为 None。
        """
        mud = next((fluid for fluid in fluids if fluid.role == FluidRole.MUD), None)
        lead = next((fluid for fluid in fluids if fluid.role == FluidRole.LEAD), None)
        tail = next((fluid for fluid in fluids if fluid.role == FluidRole.TAIL), None)
        spacer = next((fluid for fluid in fluids if fluid.role in {FluidRole.WASH, FluidRole.SPACER}), None)
        flusher = next((fluid for fluid in fluids if fluid.role == FluidRole.FLUSHER), None)
        if mud is None or (lead is None and tail is None):
            raise ValueError("需要钻井液和至少一个水泥浆流体")
        return mud, lead, tail, spacer, flusher

    @staticmethod
    def _fluid_yield_stress(fluid: FluidSpec) -> float:
        """返回流体的屈服应力。幂律和牛顿流体返回 0。"""
        if fluid.yield_stress_pa is not None:
            return fluid.yield_stress_pa
        return 0.0

    def _apparent_viscosity(
        self,
        fluid: FluidSpec,
        gamma: Array,
    ) -> Array:
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

    @staticmethod
    def _phase_power_law_params(fluid) -> tuple[float, float]:
        """把任意流变模型映射为 (幂律指数 n, 稠度 K[Pa·s^n])，供 M2 混合 n/k 加权。

        NEWTONIAN/BINGHAM -> n=1, K=plastic_viscosity_pa_s；
        POWER_LAW/HERSCHEL_BULKLEY -> (power_law_n, consistency_k)。
        HB 的屈服应力由 tau_y 场单独携带，不在这里折进 K。
        """
        if fluid is None:
            return 1.0, 1.0e-6
        rm = fluid.rheology_model
        if rm == RheologyModel.POWER_LAW or rm == RheologyModel.HERSCHEL_BULKLEY:
            return float(fluid.power_law_n), float(fluid.consistency_k)
        return 1.0, float(fluid.plastic_viscosity_pa_s)

    @staticmethod
    def _yield_gate_wall(w, b, mu_reg, tau_y, cement_ever, cement_local,
                         f_safety, c_min_residual):
        """M3 可重启屈服门槛：每深度列以流动最宽元为参考外推壁剪，
        immobile = τw_extrap ≤ f·τy；OR 残泥下限(cement_ever>0 且 cement<c_min_residual)。
        全列无流动 -> 整列冻结。停泵期不调用（run() 泵注分支门控）。"""
        b = np.maximum(b, 1e-12)
        gamma = np.maximum(6.0 * np.abs(w) / b, 1e-6)
        tau_w_field = mu_reg * gamma
        # 每列参考元：|w| 最大且 w>0；非流动元罚为 -1
        w_rank = np.where(w > 0.0, np.abs(w), -1.0)
        ref_row = np.argmax(w_rank, axis=0)            # (nz,)
        has_flow = np.any(w > 0.0, axis=0)            # (nz,)
        col = np.arange(w.shape[1])
        ref_row_safe = np.where(has_flow, ref_row, 0)
        tau_w_ref = tau_w_field[ref_row_safe, col]
        b_ref = b[ref_row_safe, col]
        G = 2.0 * tau_w_ref / np.maximum(b_ref, 1e-12)
        tau_w_extrap = G[None, :] * b / 2.0
        immobile = (tau_w_extrap <= f_safety * tau_y) & (cement_ever > 0.0)
        residual_wall = (cement_ever > 0.0) & (cement_local < c_min_residual)
        wall_new = np.where(immobile | residual_wall, 1.0, 0.0)
        # 整列冻结仅限水泥已到达的列（前锋未到无流动是正常状态，不得冻结，
        # 否则 pre-cement 全域 wall=1 堵死，门不可重启）
        col_freeze = ~has_flow & np.any(cement_ever > 0.0, axis=0)
        wall_new[:, col_freeze] = 1.0
        return wall_new.astype(float)

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
        flusher: Array,
        w_prev: Array,
        geom: Dict[str, Array],
        mud_fluid: FluidSpec,
        lead_fluid: FluidSpec | None,
        tail_fluid: FluidSpec | None,
        spacer_fluid: FluidSpec | None,
    ) -> Tuple[Array, Array, Array, Array, Array, Array, Array, Array, Array]:
        """计算混合物系的表观粘度、密度、钻井液分数、混合屈服应力、黏度比 m 场、相黏度场（η1=泥浆相, η2=水泥相）及混合幂律参数（n_mix, kappa_mix）。
        T1-6: flusher 仅参与体积闭合（五相），不参与 D2DGA 闭包。"""
        # 五相体积分数闭合：显式跟踪领浆、尾浆、前置/隔离液和冲洗液，钻井液由守恒关系反算。
        mud = np.clip(1.0 - lead - tail - spacer - flusher, 0.0, 1.0)
        effective_b = geom.get("effective_b", geom["b"])
        gamma = np.maximum(6.0 * np.abs(w_prev) / np.maximum(effective_b, 1.0e-5), 1.0e-6)
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
        # 新增：混合屈服应力（相体积加权）
        tau_y = mud * self._fluid_yield_stress(mud_fluid)
        if lead_fluid is not None:
            tau_y += lead * self._fluid_yield_stress(lead_fluid)
        if tail_fluid is not None:
            tau_y += tail * self._fluid_yield_stress(tail_fluid)
        if spacer_fluid is not None:
            tau_y += spacer * self._fluid_yield_stress(spacer_fluid)
        # R1: auto-m 黏度比场 = μ_displaced / μ_displacing
        # 被顶替液=泥浆(mu_mud)，顶替液=水泥(mu_cement)。m = mu_mud / mu_cement。
        mu_mud_field = self._apparent_viscosity(mud_fluid, gamma)
        # 水泥相表观粘度：领浆+尾浆中存在的那个（若 lead_fluid 存在用 lead，否则 tail）
        cement_fluid = lead_fluid if lead_fluid is not None else tail_fluid
        if cement_fluid is not None:
            mu_cement_field = self._apparent_viscosity(cement_fluid, gamma)
            m_field = mu_mud_field / np.maximum(mu_cement_field, 1.0e-6)
        else:
            # 无水泥相时 m=1（退化为默认），η2 退化为 η1
            mu_cement_field = mu_mud_field
            m_field = np.ones_like(mu_mud_field)
        # 限幅到合理范围，避免极端粘度比导致 f_amp 越界（d2dga_flux 内还有 [0.5,2] clip）
        m_field = np.clip(m_field, 0.1, 10.0)
        # T1-4: 返回相黏度场 η1=泥浆相, η2=水泥相（两层黏度闭包用）
        eta1 = mu_mud_field
        eta2 = mu_cement_field
        # Task1/M2: 混合物幂律参数（体积分数加权 n，对数加权 K），供 Metzner-Reed Re/He
        n_mud, k_mud = self._phase_power_law_params(mud_fluid)
        n_lead, k_lead = self._phase_power_law_params(lead_fluid)
        n_tail, k_tail = self._phase_power_law_params(tail_fluid)
        n_sp, k_sp = self._phase_power_law_params(spacer_fluid)
        n_mix = mud * n_mud + lead * n_lead + tail * n_tail + spacer * n_sp
        log_k_mix = (mud * np.log(max(k_mud, 1e-12)) + lead * np.log(max(k_lead, 1e-12))
                     + tail * np.log(max(k_tail, 1e-12)) + spacer * np.log(max(k_sp, 1e-12)))
        kappa_mix = np.exp(log_k_mix)
        return mu, rho, mud, tau_y, m_field, eta1, eta2, n_mix, kappa_mix

    def _buoyancy_force_vector(self, geom: Dict[str, Array], beta_deg: Array | float) -> Tuple[Array, Array]:
        """计算论文式 2.5b 的浮力体力向量 f = (r_a·cosβ/F², r_a·sin(πφ)·sinβ/F²)。

        R2 (I3 通量) 与 R3 (真体力) 共用。F 为 Froude 数（此处用经验常数 1.0 归一化，
        实际数值在 _compute_velocity 中按 F 定义校准）。
        返回 (f_phi, f_xi)，shape 与 geom['phi'] 广播兼容 (ny, nz)。
        """
        phi = geom["phi"][:, None]  # (ny, 1)
        beta_rad = np.deg2rad(np.asarray(beta_deg, dtype=float))
        # 平均半径 r_a（用 mean_radius 近似，从 geom 取 hole/od 推算）
        hole_mm = geom.get("hole_mm", np.full((1, self.nz), 220.0))
        od_mm = geom.get("od_mm", np.full((1, self.nz), 139.7))
        # 沿深度取均值半径（米），广播到 (ny, nz)
        r_a_m = np.mean((hole_mm + od_mm) / 4.0) / 1000.0
        # F 取 1.0（无量纲归一化；真 F 校正在 _compute_velocity 内做）
        F2 = 1.0
        f_phi = (r_a_m / F2) * np.sin(np.pi * phi) * np.sin(beta_rad)  # (ny,1) 广播
        f_xi = np.full_like(f_phi, (r_a_m / F2) * np.cos(beta_rad))
        # 广播到 (ny, nz)
        f_phi = np.broadcast_to(f_phi, (self.ny, self.nz)).astype(float, copy=True)
        f_xi = np.broadcast_to(f_xi, (self.ny, self.nz)).astype(float, copy=True)
        return f_phi, f_xi

    def _compute_buoyancy_number(self, rho_displacing_kg_m3: float, rho_displaced_kg_m3: float,
                                 gap_m: float, mu_displaced_pa_s: float, velocity_m_s: float) -> float:
        """计算无量纲浮力数 b = (ρ₂-ρ₁)·g·d² / (μ₁·w₀)（Zhang & Frigaard 2022, p.8）。

        b>0: 密度稳定（重顶替轻）；b<0: 不稳定（必须避免）；b 是垂直井主导参数。
        d = gap (半间隙) 或全间隙？论文用半间隙 d̂，这里用全间隙 gap_m/2 近似。
        """
        g = 9.81
        d_half = max(gap_m / 2.0, 1.0e-6)
        denom = max(mu_displaced_pa_s * max(velocity_m_s, 1.0e-6), 1.0e-9)
        return (rho_displacing_kg_m3 - rho_displaced_kg_m3) * g * d_half ** 2 / denom

    def _compute_velocity(
        self,
        lead: Array,
        tail: Array,
        spacer: Array,
        flusher: Array,
        geom: Dict[str, Array],
        q_m3s: float,
        w_prev: Array,
        mud_fluid: FluidSpec,
        lead_fluid: FluidSpec | None,
        tail_fluid: FluidSpec | None,
        spacer_fluid: FluidSpec | None,
        wall: Array | None = None,
    ) -> Tuple[Array, Array, Array, Array, Array, Array, Array, Array, Array, Array, Array, Array]:
        """计算环空速度场（论文D2DGA口径）。

        采用 Zhang & Frigaard (2022) 的Hele-Shaw风格速度场：
        1. 计算局部混合流体的表观粘度与密度；
        2. 以 ``b²/μ`` 构造偏心通道主导局部流动度；
        3. 根据密度差（顶替液 vs 被顶替液）计算浮力稳定系数；
        4. 用浮力修正项调整宽边/窄边速度分配；
        5. 由截面排量约束得到轴向速度 ``w``。

        T1-6: flusher 仅参与体积闭合，不参与 D2DGA 闭包/两层黏度/浮力修正。

        Returns:
            w: 井深方向速度（轴向速度）
            v: 方位角方向速度（横向速度）
            mu: 有效表观粘度场
            rho: 密度场
            mud: 钻井液分数场
            Re: 雷诺数场（仅作诊断，不参与湍流修正）
            mu_turbulent: 占位零场，保留旧结果对象兼容性
            m_field: 黏度比场 m = μ_displaced/μ_displacing（R1 auto-m）
            tau_y: 混合屈服应力场（Task1 起随 _compute_velocity 返回，Task 8/11 消费）
            eta2: 水泥相黏度场（两层黏度闭包用）
            n_mix: 混合物幂律指数场（Task1 起随速度场返回，供 M2 消费）
            kappa_mix: 混合物稠度场（Task1 起随速度场返回，供 M2 消费）
        """
        y = geom["y"]
        effective_b = geom.get("effective_b", geom["b"])
        b = effective_b
        q_half = q_m3s / 2.0
        mu, rho, mud, tau_y, m_field, eta1, eta2, n_mix, kappa_mix = self._compute_props(
            lead,
            tail,
            spacer,
            flusher,
            w_prev,
            geom,
            mud_fluid,
            lead_fluid,
            tail_fluid,
            spacer_fluid,
        )

        shear_rate = np.maximum(6.0 * np.abs(w_prev) / np.maximum(effective_b, 1.0e-5), 1.0e-6)
        # 新增：Papanastasiou 正则化屈服应力模型
        # 在低剪切区（窄边死区）显著增大有效黏度，使局部流度趋近于零
        # 注意：_apparent_viscosity 对 Bingham/HB 已包含 tau_y/gamma 项，
        # 正则化应替换此项而非叠加，故先减去纯屈服应力贡献
        gamma_safe = np.maximum(shear_rate, 1.0e-8)
        regularization_factor = (1.0 - np.exp(-self.yield_regularization_M * shear_rate)) / gamma_safe
        # 减去 _apparent_viscosity 中的纯屈服应力贡献（tau_y / gamma），
        # 再用 Papanastasiou 正则化项替代；对牛顿/幂律流体 tau_y=0，不影响
        mu_shear = np.maximum(mu - tau_y / gamma_safe, 1.0e-6)
        mu_reg = mu_shear + tau_y * regularization_factor

        D_h = 2.0 * geom["b"]
        rho_kg_m3 = rho * 1000.0
        Re = rho_kg_m3 * np.abs(w_prev) * D_h / np.maximum(mu_reg, 1.0e-6)
        mu_turbulent = np.zeros_like(mu_reg)

        # === 论文D2DGA口径速度场：偏心通道主导 + 浮力修正 ===
        # T1-4: 两层黏度闭包 1/η_mix = c̄³/η₂ + (1−c̄³)/η₁（式 4.23）替换单相 μ_reg
        # 基础流动度：偏心通道主导 (b/mean(b))^2 / η_mix
        # mu_reg 保留用于 Re 诊断
        b_mean = np.mean(b, axis=0, keepdims=True)
        c_bar = np.clip(lead + tail, 0.0, 1.0)  # 局部水泥浓度
        eta_mix = 1.0 / (c_bar**3 / np.maximum(eta2, 1.0e-9)
                         + (1.0 - c_bar**3) / np.maximum(eta1, 1.0e-9))
        # T1-3b: I₁(c̄,m) 乘子（Zhang 2022 式 4.22，S ∝ 1/(2I₁)）
        # 牛顿极限 m→1 时 I₁=1/3，不改变 base 形状；m≠1 时修正方位分布
        m_local = float(np.mean(m_field)) if np.all(np.isfinite(m_field)) else self.d2dga_viscosity_ratio
        i1_base = d2dga_dispersion_I1(c_bar, m_local)
        base = (b / np.maximum(b_mean, 1.0e-12)) ** 2 / np.maximum(eta_mix, 1.0e-9)
        base = base * i1_base  # I₁ 乘子

        # === 速度场流动度：偏心通道主导 + 浮力修正 ===
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

        phi = geom["phi"][:, None]
        ebar = geom["e"][None, :]
        if self.enable_true_buoyancy and self.enable_d2dga:
            # T1-3b: 体力向量注入流动度（式 2.5b/4.24），替换 (2φ−1) 简化代理
            beta_deg_local = float(np.mean(geom.get("inc_deg", np.zeros(self.nz))))
            f_phi_arr, _ = self._buoyancy_force_vector(geom, beta_deg_local)
            rho_displaced = mud_fluid.density_kg_m3 / 1000.0
            delta_rho = (rho - rho_displaced)  # g/cc 局部密度差
            i2 = d2dga_dispersion_I2(c_bar, m_local)
            # 式 4.24: 方位修正 = Δρ·f_phi·(I2/I1)；重顶替轻→窄边(f_phi 大) pref 提升
            correction = np.clip(delta_rho * (i2 / np.maximum(i1_base, 1.0e-12)), -0.5, 0.5)
            buoyancy_shape = 1.0 + correction * f_phi_arr
        else:
            # R0/R1/R2: 保留 buoyancy_shape 代理（旧论文状态）
            density_contrast = (rho_disp - mud_fluid.density_kg_m3) / mud_fluid.density_kg_m3
            stable = float(np.clip(8.0 * density_contrast, -0.35, 0.45))
            buoyancy_shape = 1.0 + stable * ebar * (2.0 * phi - 1.0)
        pref = np.maximum(base * buoyancy_shape, 1.0e-8)
        # T1-5: wall=1 处壁面静止层，流动度归零（式 2.35-2.41）
        if wall is not None:
            pref = pref * (1.0 - wall)

        # 由截面排量约束得到轴向速度 w
        dy = np.gradient(geom["y"])[:, None]
        if self.enable_regime_split:
            # M2: 局部流态修正固定点迭代（Maleki & Frigaard 2017 式58-66）
            # 浓度相关量（base/buoyancy_shape/wall/黏度/密度/b）在迭代外缓存；
            # 黏度场保持 w_prev 一步滞后（既有约定），迭代只重算 Re_p/R/pref/area_weight/w。
            # rho_kg_m3 与 kappa_mix(Pa·s^n)/tau_y(Pa) 同单位系，保证 Re_p/He 无量纲正确。
            from cemdisp.models2d import regime_closure as rc
            wall_factor = np.ones_like(base) if wall is None else (1.0 - wall)
            he = rc.hedstrom_number(tau_y, rho_kg_m3, n_mix, kappa_mix, b)
            re_crit = 2100.0 * (1.0 + 0.1 * he)  # 屈服推迟转捩；re_crit(He) 标定钮（provisional）
            w_k = w_prev.copy()
            R = np.ones_like(w_k)
            for _ in range(self.regime_max_iter):
                re_p = rc.metzner_reed_re(w_k, rho_kg_m3, n_mix, kappa_mix, b)
                R_new, _ = rc.drag_weight(re_p, he, n_mix, re_crit, self.regime_re_turb_ratio)
                pref_k = np.maximum(base * buoyancy_shape * R_new, 1.0e-8) * wall_factor
                area_w = np.sum(pref_k * b * dy * 2.0, axis=0, keepdims=True)
                w_raw = q_half * pref_k / np.maximum(area_w, 1.0e-12)
                w_new = self.regime_relax_alpha * w_raw + (1.0 - self.regime_relax_alpha) * w_k
                if (np.max(np.abs(w_new - w_k))
                        < self.regime_tol_rel * max(np.max(np.abs(w_k)), 1e-12)):
                    w_k = w_new
                    R = R_new
                    break
                w_k = w_new
                R = R_new
            # 欠松弛只是求解 R 的迭代手段；报告的 w 必须是以最终 R 直接归一的结果，
            # 使 2·Σw·b·dy = q_half 对任意迭代步数精确成立（真正的不动点），
            # 消除欠松弛迭代返回 w_k 时 ~1.6% 的瞬态守恒误差。
            pref_final = np.maximum(base * buoyancy_shape * R, 1.0e-8) * wall_factor
            area_final = np.sum(pref_final * b * dy * 2.0, axis=0, keepdims=True)
            w = q_half * pref_final / np.maximum(area_final, 1.0e-12)
            # area_weight 与最终 w 构成一致配对（同一 area_final）
            area_weight = area_final
        else:
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

        # 返回正则化后的黏度 mu_reg，确保下游 mobility 指标反映屈服死区效应
        return w, v, mu_reg, rho, mud, Re, mu_turbulent, m_field, tau_y, eta2, n_mix, kappa_mix

    def _compute_cfl_dt_step(self, w: Array, v: Array, geom: Dict[str, Array], current_time_s: float) -> float:
        """根据 CFL 条件计算自适应时间步长。

        CFL 条件：dt_step * (max|w|/Δs + max|v|/Δy) ≤ cfl_number

        Args:
            w: 轴向速度场 (ny×nz)
            v: 横向速度场 (ny×nz)
            geom: 几何参数字典
            current_time_s: 当前时间（秒）

        Returns:
            dt_step: 自适应时间步长，满足 dt_min ≤ dt_step ≤ min(dt, total_t - current_time_s)
        """
        ds = float(np.min(np.diff(geom["s"])))
        dy_arr = np.gradient(geom["y"])[:, None]
        dy_min = float(np.min(dy_arr)) if dy_arr.size else ds
        denom = max(
            float(np.max(np.abs(w))) / max(ds, 1e-12)
            + float(np.max(np.abs(v))) / max(dy_min, 1e-12),
            1e-12,
        )
        dt_cfl = self.cfl_number / denom
        dt_step = min(dt_cfl, self.dt, self.total_t - current_time_s)
        dt_step = max(dt_step, self.dt_min)
        # 最后步裁剪精确到 total_t
        if current_time_s + dt_step >= self.total_t:
            dt_step = self.total_t - current_time_s
        return dt_step

    def _depth_profiles(self, geom: Dict[str, Array], lead: Array, tail: Array, spacer: Array, flusher: Array) -> pd.DataFrame:
        """计算深度方向的平均剖面数据。T1-6: 含冲洗液平均浓度。"""
        cement = np.clip(lead + tail, 0.0, 1.0)
        eff = cement
        mud = np.clip(1.0 - lead - tail - spacer - flusher, 0.0, 1.0)
        return pd.DataFrame(
            {
                "井深_m": geom["md"],
                "领浆平均浓度": np.average(lead, axis=0, weights=geom["b"]),
                "尾浆平均浓度": np.average(tail, axis=0, weights=geom["b"]),
                "水泥平均浓度": np.average(cement, axis=0, weights=geom["b"]),
                "前置液隔离液平均浓度": np.average(spacer, axis=0, weights=geom["b"]),
                "冲洗液平均浓度": np.average(flusher, axis=0, weights=geom["b"]),
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

        mud_fluid, lead_fluid, tail_fluid, spacer_fluid, flusher_fluid = self._pick_fluids(fluids)
        geom = self._build_geom(well_spec)
        lead = np.zeros((self.ny, self.nz), dtype=float)
        tail = np.zeros((self.ny, self.nz), dtype=float)
        spacer = np.zeros((self.ny, self.nz), dtype=float)
        flusher = np.zeros((self.ny, self.nz), dtype=float)  # T1-6: FLUSHER 独立浓度场
        # wall 场初始化为零；T1-5 后按 c < c_min 判据在泵注阶段动态更新（式 2.35-2.41）
        wall = np.zeros((self.ny, self.nz), dtype=float)
        # 水泥前锋到达标记：c_min 壁面判据只在前锋已到达的网格生效
        cement_ever = np.zeros((self.ny, self.nz), dtype=float)
        w_prev = np.full((self.ny, self.nz), 0.45, dtype=float)
        half_volume = _trapez2d(geom["b"], geom)

        ygrid, sgrid = np.meshgrid(geom["y"], geom["s"], indexing="ij")
        rows: list[list[float | str]] = []
        mud_density_gcc = mud_fluid.density_kg_m3 / 1000.0
        cement_snapshots: list[Array] = []
        lead_snapshots: list[Array] = []
        tail_snapshots: list[Array] = []
        spacer_snapshots: list[Array] = []
        flusher_snapshots: list[Array] = []  # T1-6
        wall_snapshots: list[Array] = []
        snapshot_times: list[float] = []
        cumulative_lead_in_m3 = 0.0
        cumulative_tail_in_m3 = 0.0
        cumulative_spacer_in_m3 = 0.0
        cumulative_flusher_in_m3 = 0.0  # T1-6

        # T1-7: CFL 自适应 → while 循环（每步 current_time_s += dt_step），
        # 固定 dt → for 仿真（current_time_s = step_index * dt，复现基线）
        # 统一循环体，仅在循环入口/出口按 enable_cfl_adaptive 分支
        step_index = 0
        current_time_s = 0.0
        _progress_next_pct = 0.0  # 下一个需要打印的进度百分比
        print(f"[D2DGA] 开始环空二维模拟 total_t={self.total_t:.0f}s nz={self.nz} ny={self.ny}")
        if not self.enable_cfl_adaptive:
            final_step_index = int(self.total_t / self.dt)

        # 单 while 共享循环体（CFL 分支在步进后判 break，固定 dt 分支在步前判 break）。
        # 若分拆两循环则 ~250 行泵注/停泵/指标逻辑重复，等价于 spec 条件 current_time_s < self.total_t。
        while True:
            if not self.enable_cfl_adaptive:
                # 用 while+break 而非 for step_index in range(...)：因循环体与 CFL 分支共享，
                # 拆分两循环将重复 ~250 行逻辑（行 798-1008），保持单循环体减少重复。
                if step_index > final_step_index:
                    break
                current_time_s = step_index * self.dt
                dt_step = self.dt
            inlet_state = inlet_state_provider(current_time_s)
            inlet_cement_fraction = _phase_fraction(inlet_state, "cement")
            inlet_lead_fraction = _phase_fraction(inlet_state, "lead")
            inlet_tail_fraction = _phase_fraction(inlet_state, "tail") + inlet_cement_fraction
            inlet_spacer_fraction = _phase_fraction(inlet_state, "spacer")
            inlet_flusher_fraction = _phase_fraction(inlet_state, "flusher")  # T1-6

            # 泵停判断：排量低于阈值时认为泵已停止。
            # 泵停后水泥场冻结——不再平流、扩散或壁面清除，
            # 因为水泥静凝胶强度在短时间（<2h）内足以抵抗浮力滑塌。
            pump_active = inlet_state.flow_rate_m3_s > 1.0e-9

            if pump_active:
                # === 正常泵注阶段：仅执行论文口径核心平流 + D2DGA 通量修正 ===
                w, v, mu, rho, mud, Re, mu_turbulent, m_field, _tau_y, _eta2, _n_mix, _kappa_mix = self._compute_velocity(
                    lead,
                    tail,
                    spacer,
                    flusher,  # T1-6
                    geom,
                    inlet_state.flow_rate_m3_s,
                    w_prev,
                    mud_fluid,
                    lead_fluid,
                    tail_fluid,
                    spacer_fluid,
                    wall=wall,
                )
                w_prev = w

                # T1-7: CFL 自适应时间步（dt 计算，循环结构 Task 2 改）
                if self.enable_cfl_adaptive:
                    dt_step = self._compute_cfl_dt_step(w, v, geom, current_time_s)
                else:
                    dt_step = self.dt

                cement = np.clip(lead + tail, 0.0, 1.0)
                if self.enable_d2dga:
                    if self.enable_d2dga_auto_m:
                        m_for_amp = m_field  # 数组
                    else:
                        m_for_amp = self.d2dga_viscosity_ratio  # 标量
                    f_amp = d2dga_flux_amplification(cement, m_for_amp)
                else:
                    f_amp = 1.0

                w_d2dga = w * f_amp
                v_d2dga = v * f_amp
                ysrc = ygrid - v_d2dga * dt_step
                ssrc = sgrid - w_d2dga * dt_step
                lead_adv = _bilinear_interp(lead, ysrc, ssrc, geom, inlet_lead_fraction)
                tail_adv = _bilinear_interp(tail, ysrc, ssrc, geom, inlet_tail_fraction)
                spacer_adv = _bilinear_interp(spacer, ysrc, ssrc, geom, inlet_spacer_fraction)
                flusher_adv = _bilinear_interp(flusher, ysrc, ssrc, geom, inlet_flusher_fraction)  # T1-6
                lead = np.clip(lead_adv, 0.0, 1.0)
                tail = np.clip(tail_adv, 0.0, 1.0)
                spacer = np.clip(spacer_adv, 0.0, 1.0)
                flusher = np.clip(flusher_adv, 0.0, 1.0)  # T1-6
                # 数值扩散可能使显式相之和略超1；按比例压回可行域，保持泥浆分数非负。
                tracked_total = lead + tail + spacer + flusher  # T1-6: 五相过填修正
                overfilled = tracked_total > 1.0
                lead[overfilled] /= tracked_total[overfilled]
                tail[overfilled] /= tracked_total[overfilled]
                spacer[overfilled] /= tracked_total[overfilled]
                flusher[overfilled] /= tracked_total[overfilled]  # T1-6

                # D2DGA间隙尺度弥散：在低浓度前锋更强，模拟间隙尺度分散效应。
                # 数值弥散可能使显式相之和略超 1；后续两次 overfilled 修正将其压回可行域，
                # 允许不超过 1e-12 的数值扩散容差。
                # M1: 弥散系数按 dt 归一（恢复量纲正确性）。CFL 自适应使 dt 降到 ~0.118s，
                # 旧硬编码是"每步固定幅值"→单位物理时间弥散放大 dt_ref/dt_step≈34 倍。
                # _dt_norm = scale * dt_step/dt_ref：固定 dt 模式 dt_step==dt_ref 且 scale=1
                # 时系数==基线硬编码（0.018/0.015/0.012），逐位复现；CFL 下每物理秒弥散恒定。
                _dt_norm = self.dispersion_dt_scale * (dt_step / self.dispersion_dt_ref)
                _ax = self.dispersion_axial * _dt_norm
                _az = self.dispersion_azimuthal * _dt_norm
                # spacer/flusher 基础弥散系数为 0.012/0.012（轴向/方位角同值，独立于 lead/tail）。
                # 必须用字面量 0.012，而非 0.018*0.667（=0.012006）或 0.015*0.8：
                # 默认(fixed dt=4, scale=1)下要求 0.012*1.0==0.012 与基线硬编码逐位复现。
                _ax_sf = 0.012 * _dt_norm
                _az_sf = 0.012 * _dt_norm
                lead = self._smooth_dispersion(lead, axial=_ax, azimuthal=_az)
                tail = self._smooth_dispersion(tail, axial=_ax, azimuthal=_az)
                spacer = self._smooth_dispersion(spacer, axial=_ax_sf, azimuthal=_az_sf)
                flusher = self._smooth_dispersion(flusher, axial=_ax_sf, azimuthal=_az_sf)
                # T1-6: 弥散后再次执行五相过填修正，防止 _smooth_dispersion 数值扩散
                # 使 lead+tail+spacer+flusher 再次超过 1，破坏体积分数闭合。
                tracked_total = lead + tail + spacer + flusher
                overfilled = tracked_total > 1.0
                lead[overfilled] /= tracked_total[overfilled]
                tail[overfilled] /= tracked_total[overfilled]
                spacer[overfilled] /= tracked_total[overfilled]
                flusher[overfilled] /= tracked_total[overfilled]

                # R2: I3 浮力弥散通量（式 4.25 第二项）—— 仅作用于水泥相(lead+tail)
                if self.enable_d2dga_i3_flux and self.enable_d2dga:
                    cement_for_flux = np.clip(lead + tail, 0.0, 1.0)
                    # 浮力向量 f（用当前井段平均井斜）
                    beta_deg_local = float(np.mean(geom["inc_deg"])) if "inc_deg" in geom else 0.0
                    f_phi_arr, f_xi_arr = self._buoyancy_force_vector(geom, beta_deg_local)
                    # 顶替液粘度 eta2 + 密度差 Δρ（顶替液 - 被顶替液），kg/m³。
                    # 默认关（enable_local_i3=False）：全场均值，逐位复现基线；
                    # 开启后：eta2 透传水泥相黏度场 _eta2、Δρ 用局部混合密度场，实现 I3 局部化。
                    if self.enable_local_i3:
                        eta2 = _eta2 if np.all(np.isfinite(_eta2)) else float(np.mean(_eta2))
                        delta_rho = (rho - mud_density_gcc) * 1000.0
                    else:
                        eta2 = float(np.mean(mu)) if np.all(np.isfinite(mu)) else 0.18
                        delta_rho = (rho.mean() - mud_density_gcc) * 1000.0
                    H_field = geom["H"]
                    m_for_flux = m_field if self.enable_d2dga_auto_m else self.d2dga_viscosity_ratio
                    q_phi, q_xi = d2dga_buoyancy_flux(
                        cement_for_flux, m_for_flux, delta_rho, H_field, eta2,
                        f_phi_arr, f_xi_arr,
                    )
                    # 散度通量：dc/dt += -div(q) = -(dq_xi/ds + dq_phi/dy)
                    # 注意：np.gradient(f, x, axis) 中 x 是坐标数组（不是间距）。
                    # 若用 np.gradient(geom["y"]) 返回均匀间距数组，再作为坐标传入会导致
                    # 除零间距 -> inf/NaN。此处直接用 geom["s"] (nz,) 和 geom["y"] (ny,)
                    # 作为坐标参数，与 q_xi.shape[1] 和 q_phi.shape[0] 匹配。
                    dq_xi_ds = np.gradient(q_xi, geom["s"], axis=1)    # geom["s"]: 轴向坐标 (nz,)
                    dq_phi_dy = np.gradient(q_phi, geom["y"], axis=0)   # geom["y"]: 方位角坐标 (ny,)
                    div_q = dq_xi_ds + dq_phi_dy
                    # 通量只加到水泥相（lead+tail 按比例分配）
                    cement_total = np.maximum(cement_for_flux, 1.0e-6)
                    lead_frac = lead / cement_total
                    tail_frac = tail / cement_total
                    # T1-2: 去人工限幅 flux_strength=0.05；物理系数 ΔρH³/(6η₂)·I3 直驱（式 4.25）
                    # 局部 CFL 裁剪防单步越界（非全局限幅）
                    # ds 取轴向网格最小间距，保证非均匀网格下 CFL 条件保守
                    ds = float(np.min(np.diff(geom["s"])))
                    step_limit = self.alpha_cfl * ds / max(dt_step, 1.0e-9)
                    div_q_clipped = np.clip(div_q, -step_limit, step_limit)
                    lead = lead - div_q_clipped * lead_frac * dt_step
                    tail = tail - div_q_clipped * tail_frac * dt_step
                    lead = np.clip(lead, 0.0, 1.0)
                    tail = np.clip(tail, 0.0, 1.0)

                # D2DGA 通量放大会改变前锋形态，但不应让各相总量超过累计入环空体积。
                # 这里按入口累计体积对领浆、尾浆、前置液/隔离液和冲洗液分别做体积上限约束。
                cumulative_lead_in_m3 += inlet_state.flow_rate_m3_s * inlet_lead_fraction * dt_step
                cumulative_tail_in_m3 += inlet_state.flow_rate_m3_s * inlet_tail_fraction * dt_step
                cumulative_spacer_in_m3 += inlet_state.flow_rate_m3_s * inlet_spacer_fraction * dt_step
                cumulative_flusher_in_m3 += inlet_state.flow_rate_m3_s * inlet_flusher_fraction * dt_step  # T1-6
                lead = _limit_phase_volume(lead, geom, cumulative_lead_in_m3, self.open_outlet)
                tail = _limit_phase_volume(tail, geom, cumulative_tail_in_m3, self.open_outlet)
                spacer = _limit_phase_volume(spacer, geom, cumulative_spacer_in_m3, self.open_outlet)
                flusher = _limit_phase_volume(flusher, geom, cumulative_flusher_in_m3, self.open_outlet)  # T1-6

                # T1-5: static wall layer c_min 判据（Bararpour 2025 式 2.35-2.41）
                # 局部水泥浓度 c < c_min 处壁面层泥浆滞留不流动 → wall=1
                # 注意：c_min 判据只在水前锋已到达的网格生效（cement_ever > 0），
                # 避免前锋到达前全局 wall=1 堵塞速度场。
                cement_local = np.clip(lead + tail, 0.0, 1.0)
                cement_ever = np.maximum(cement_ever, cement_local)
                # M3: 屈服门槛（Bararpour 2025）—— 开启时用 w/mu_reg/tau_y 重建壁面冻结层，
                # 使壁面层在排量回升后可以解冻（c_min 判据做不到）。默认关闭，else 分支为
                # 原始 T1-5 c_min 判据，逐位复现基线。
                if self.enable_yield_gate:
                    wall = self._yield_gate_wall(
                        w, geom["effective_b"], mu, _tau_y, cement_ever, cement_local,
                        self.yield_gate_f_safety, self.yield_gate_c_min_residual)
                else:
                    wall = np.where(cement_ever > 0, (cement_local < self.c_min).astype(float), 0.0)

            else:
                # === 泵停阶段：冻结浓度场，仅记录指标 ===
                w, v, mu, rho, mud, Re, mu_turbulent, m_field, _tau_y, _eta2, _n_mix, _kappa_mix = self._compute_velocity(
                    lead,
                    tail,
                    spacer,
                    flusher,  # T1-6
                    geom,
                    0.0,
                    w_prev,
                    mud_fluid,
                    lead_fluid,
                    tail_fluid,
                    spacer_fluid,
                    wall=wall,
                )
                # 泵停后保持上一时刻浓度场，不再引入停泵滑移或额外壁面过程。
                # T1-7: 泵停分支也补 dt_step 计算（while 循环需 dt_step 推进时间）
                if self.enable_cfl_adaptive:
                    dt_step = self._compute_cfl_dt_step(w, v, geom, current_time_s)
                else:
                    dt_step = self.dt

            # T1-7: 统一步后时间（current_time_s + dt_step），CFL/固定 dt 分支语义一致
            record_time = min(current_time_s + dt_step, self.total_t)
            _last_step = (
                current_time_s + dt_step >= self.total_t - 1e-9
                if self.enable_cfl_adaptive
                else step_index == int(self.total_t / self.dt)
            )

            # 在物理场更新后、指标计算前保存快照，确保快照与本步指标使用同一状态。
            # 使用 copy() 固化二维场，避免后续时间步原地更新影响已保存结果。
            if step_index % self.save_interval == 0 or _last_step:
                cement = np.clip(lead + tail, 0.0, 1.0)
                cement_snapshots.append(cement.copy())
                lead_snapshots.append(lead.copy())
                tail_snapshots.append(tail.copy())
                spacer_snapshots.append(spacer.copy())
                flusher_snapshots.append(flusher.copy())  # T1-6
                wall_snapshots.append(wall.copy())
                snapshot_times.append(record_time)

            cement = np.clip(lead + tail, 0.0, 1.0)
            eff = cement
            bulk_fill = _trapez2d(geom["b"] * cement, geom) / half_volume
            effective_efficiency = _trapez2d(geom["b"] * eff, geom) / half_volume

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
                    record_time,
                    record_time / 60.0,
                    inlet_state.stage_name,
                    bulk_fill,
                    effective_efficiency,
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
                    float(np.mean(flusher)),  # T1-6
                ]
            )

            # T1-7: 时间步进（while 分支：current_time_s += dt_step；for 分支：step_index 自增）
            if self.enable_cfl_adaptive:
                current_time_s += dt_step
                # 进度输出：每 10% 打印一次
                pct = current_time_s / self.total_t * 100.0
                if pct >= _progress_next_pct:
                    print(f"  [D2DGA] 进度 {pct:.0f}% ({current_time_s:.0f}s/{self.total_t:.0f}s, dt={dt_step:.3f}s)")
                    _progress_next_pct += 10.0
                if current_time_s >= self.total_t - 1e-9:
                    break
            else:
                # 固定 dt 模式的进度输出
                pct = current_time_s / self.total_t * 100.0
                if pct >= _progress_next_pct:
                    print(f"  [D2DGA] 进度 {pct:.0f}% ({current_time_s:.0f}s/{self.total_t:.0f}s)")
                    _progress_next_pct += 10.0
            step_index += 1

        metric_columns = [
            "time_s",
            "time_min",
            "stage",
            "bulk_cement_fill",
            "effective_efficiency",
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
            "mean_flusher",  # T1-6
        ]
        # T1-7: CFL 自适应模式下 dt 缩小 ~34 倍（4s→~0.118s），metrics 行数相应膨胀 ~34 倍，
        # 属预期行为（数值扩散锐减的代价）；后续可加采样降频优化，非阻塞。
        metrics = pd.DataFrame(data=rows, columns=pd.Index(metric_columns))
        cement = np.clip(lead + tail, 0.0, 1.0)
        depth_profiles = self._depth_profiles(geom, lead, tail, spacer, flusher)  # T1-6
        final = metrics.iloc[-1]

        # M0: 失稳指数去饱和——线性代理与对数代理（log10(1+proxy)）进 summary
        _inst_lin = float(final["instability_proxy"])
        _inst_log = float(np.log10(1.0 + _inst_lin))

        # Extract values into locals so both nested Chinese keys and top-level English
        # aliases use one source, avoiding drift.
        eff_efficiency = float(final["effective_efficiency"])
        cement_occ = float(final["bulk_cement_fill"])
        chan_idx = float(final["channeling_index"])
        mix_idx = float(final["mixing_index"])
        inst_idx = float(final["instability_index"])

        # R3: 输出无量纲浮力数 b（主导参数，p.27）
        rho_displacing = (
            lead_fluid.density_kg_m3 if lead_fluid
            else tail_fluid.density_kg_m3 if tail_fluid
            else mud_fluid.density_kg_m3
        )
        b_number = self._compute_buoyancy_number(
            rho_displacing_kg_m3=rho_displacing,
            rho_displaced_kg_m3=mud_fluid.density_kg_m3,
            gap_m=float(np.mean(geom["b"])),
            mu_displaced_pa_s=(mud_fluid.plastic_viscosity_pa_s or 0.05),  # None guard for non-Bingham muds (power-law/HB); 0.05 Pa·s fallback
            velocity_m_s=float(np.mean(np.abs(w_prev))),
        )

        summary: Dict[str, object] = {
            "模型名称": "通用尾管段环空二维顶替模型",
            "模拟对象": f"{well_spec.well_name} 尾管段 {well_spec.top_md_m:.2f}-{well_spec.bottom_md_m:.2f}m",
            "井段_m": [well_spec.top_md_m, well_spec.bottom_md_m],
            "物理环空体积_m3": self._physical_annular_volume(well_spec),
            "最终结果": {
                "全井段最终有效顶替效率": eff_efficiency,
                "窄四分位效率": _narrow_quarter_efficiency(cement, geom),
                "最终水泥浆占据率": cement_occ,
                "最终窜槽指数": chan_idx,
                "最终混浆指数": mix_idx,
                "最终失稳指数": inst_idx,
                "最终失稳指数_线性": _inst_lin,
                "最终失稳指数_对数": _inst_log,
                "浮力数_b": b_number,
            },
            "评价窗效率": _evaluation_window_efficiencies(well_spec, geom, cement),
            "低尾指标": _low_tail_indicators(geom, cement, self.ny),
            "effective_efficiency": eff_efficiency,
            "eta_narrow": _narrow_quarter_efficiency(cement, geom),
            "channeling_index": chan_idx,
            "mixing_index": mix_idx,
            "buoyancy_number": b_number,
        }
        result = AnnulusSimulationResult(
            well_name=well_spec.well_name,
            geom=geom,
            cement_field=cement,
            spacer_field=spacer,
            flusher_field=flusher,  # T1-6
            wall_field=wall,
            metrics=metrics,
            depth_profiles=depth_profiles,
            summary=summary,
            time_points_s=tuple(float(value) for value in metrics["time_s"].to_list()),
            cement_snapshots=tuple(cement_snapshots),
            lead_snapshots=tuple(lead_snapshots),
            tail_snapshots=tuple(tail_snapshots),
            spacer_snapshots=tuple(spacer_snapshots),
            flusher_snapshots=tuple(flusher_snapshots),  # T1-6
            wall_snapshots=tuple(wall_snapshots),
            snapshot_times_s=tuple(snapshot_times),
            notes=(
                "当前显式跟踪领浆、尾浆与前置液/隔离液/冲洗液五类入环空相，钻井液由体积分数闭合反算。",
                "solver 核心仅保留论文口径的几何、流动度、D2DGA 通量修正与求解域效率输出。",
                "泥饼、温度、凝胶强度、湍流修正与 CBL 质量惩罚不再影响求解结果；相关字段仅作兼容占位。",
            ),
            lead_field=lead,
            tail_field=tail,
        )

        # Tier 0 诊断聚合：纯后处理，注入 result.summary（dict 可安全追加）
        # 诊断失败不影响主求解，但记录错误原因（不再静默吞掉），便于排查接线问题
        try:
            from cemdisp.diagnostics.tier0_diagnostics import compute_all_tier0_diagnostics
            tier0 = compute_all_tier0_diagnostics(result, fluids=fluids, well_spec=well_spec)
            result.summary["tier0_diagnostics"] = tier0.to_dict()  # type: ignore[index]
        except Exception as exc:
            result.summary["tier0_diagnostics_error"] = f"{type(exc).__name__}: {exc}"

        return result
