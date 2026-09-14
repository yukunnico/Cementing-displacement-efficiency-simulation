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

import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Callable, Dict, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray
import pandas as pd

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.well_spec import DepthValuePoint, WellSpec
from cemdisp.diagnostics.displacement_metrics import _narrow_quarter_efficiency
from cemdisp.models2d import buoyancy
from cemdisp.models2d.boundary_bridge import AnnulusInletState
from cemdisp.models2d.d2dga_flux import (
    d2dga_buoyancy_flux,
    d2dga_dispersion_I1,
    d2dga_dispersion_I2,
    d2dga_flux_amplification,
)

if TYPE_CHECKING:  # 仅类型注解，运行时不引入 data.pumping_schedule 依赖
    from cemdisp.data.pumping_schedule import PumpingSchedule


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
        e_clip_max: float = 0.55,
        enable_yield_gate: bool = True,  # 2026-09-02 默认启用可逆τw物理屈服门（替代非物理永久浓度冻结，结果网格收敛）
        yield_gate_f_safety: float = 1.15,
        dispersion_axial: float = 0.018,
        dispersion_azimuthal: float = 0.015,
        dispersion_dt_ref: float = 4.0,
        dispersion_dt_scale: float = 1.0,
        enable_e_clip_ruling: bool = True,
        e_clip_measured_max: float = 0.90,
        enable_power_law_gap_law: bool = True,
    ) -> None:
        """初始化环空二维求解器参数。

        Args:
            dt: 时间步长（秒），默认4秒
            nz: 井深方向网格数，默认140
            ny: 方位角方向网格数，默认40
            total_t: 总模拟时间（秒），默认12000秒（200分钟）
            enable_d2dga: 是否启用D2DGA通量修正（Zhang & Frigaard 2022），默认开启
                2026-09-07 R0 分支删除：黏度比 m 恒由 auto-m 自动计算（旧
                d2dga_viscosity_ratio 标量路径为旧论文 R0 状态，已移除）。
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
            e_clip_max: M4 偏心度 e 硬截断上限，默认 0.55（逐位复现基线）。
                由 e = clip(1-standoff, 0.05, e_clip_max) 构造几何；
                生产跑道（Task 13 重跑阶段）显式设 0.90 放宽截断。
                体积校正（_build_geom 末尾 scale）每个 run 重算，half_volume 守恒。
            enable_e_clip_ruling: e_clip 按井数据来源自动裁定开关（2026-09-06），默认 True。
                True: well_spec.standoff_measured=True（实测居中度剖面井，如呼101实测口径）
                时截断上限取 e_clip_measured_max=0.90，否则维持 e_clip_max（设计/代理
                值井保持保守截断——在假设输入上放大模型响应会制造伪敏感性）；
                False: 恒用 e_clip_max（旧口径）。
                显式传入 e_clip_max 时该井仍按裁定选择上限（显式值覆盖设计值井路径）。
            enable_yield_gate: M3 屈服门槛开关，默认 True（2026-09-02 起物理屈服门）。
                True 时 pump 分支用 _yield_gate_wall 重建壁面冻结层；
                False（2026-09-07 起无 c_min 兜底轨）= wall 恒零（无壁面静止层）。
            yield_gate_f_safety: 屈服门槛安全系数 f，默认 1.15。
                immobile 判定：外推壁剪 τw_extrap ≤ f·τy。
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
        self.instability_decay_scale = instability_decay_scale
        self.save_interval: int = save_interval
        self.yield_regularization_M: float = yield_regularization_M
        self.enable_regime_split: bool = enable_regime_split
        self.regime_relax_alpha: float = regime_relax_alpha
        self.regime_max_iter: int = regime_max_iter
        self.regime_tol_rel: float = regime_tol_rel
        self.regime_re_turb_ratio: float = regime_re_turb_ratio
        self.enable_d2dga_i3_flux: bool = enable_d2dga_i3_flux
        self.enable_local_i3: bool = enable_local_i3
        self.enable_true_buoyancy: bool = enable_true_buoyancy
        self.open_outlet: bool = open_outlet
        self.alpha_cfl: float = alpha_cfl
        self.enable_cfl_adaptive: bool = enable_cfl_adaptive
        self.cfl_number: float = cfl_number
        self.dt_min: float = dt_min
        self.e_clip_max: float = e_clip_max
        self.enable_yield_gate: bool = enable_yield_gate
        self.yield_gate_f_safety: float = yield_gate_f_safety
        self.dispersion_axial = dispersion_axial
        self.dispersion_azimuthal = dispersion_azimuthal
        self.dispersion_dt_ref = dispersion_dt_ref
        self.dispersion_dt_scale = dispersion_dt_scale
        # 2026-09-06 e_clip 裁定 + 幂律缝隙律（构造参数见 docstring）
        self.enable_e_clip_ruling = enable_e_clip_ruling
        self.e_clip_measured_max = e_clip_measured_max
        self.enable_power_law_gap_law = enable_power_law_gap_law

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

        # 2026-09-06 e_clip 裁定：实测居中度井（standoff_measured=True）放开截断到
        # e_clip_measured_max（默认 0.90）；设计/代理值井维持 e_clip_max（默认 0.55），
        # 避免在假设输入上放大模型响应。enable_e_clip_ruling=False 退回旧口径。
        e_cap = self.e_clip_max
        if self.enable_e_clip_ruling and getattr(well_spec, "standoff_measured", False):
            e_cap = max(self.e_clip_max, self.e_clip_measured_max)
        e = np.clip(1.0 - standoff, 0.05, e_cap)
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

        2026-09-06 选相修复：当 WASH 与 SPACER 角色并存时（如"平衡液"+"驱油隔离液"），
        返回体积加权等效代表流体（_composite_spacer_fluid），不再取第一个 WASH/SPACER——
        旧口径使真实隔离液（高密度幂律）物性在 2D 闭包中失效，密度场/浮力项失真。
        """
        mud = next((fluid for fluid in fluids if fluid.role == FluidRole.MUD), None)
        lead = next((fluid for fluid in fluids if fluid.role == FluidRole.LEAD), None)
        tail = next((fluid for fluid in fluids if fluid.role == FluidRole.TAIL), None)
        wash_or_spacer = [fluid for fluid in fluids
                          if fluid.role in {FluidRole.WASH, FluidRole.SPACER}]
        spacer: FluidSpec | None = None
        if len(wash_or_spacer) == 1:
            spacer = wash_or_spacer[0]
        elif len(wash_or_spacer) > 1:
            # 多种 WASH/SPACER 并存：等权合成占位，run() 内按入库体积权重重建（见下）。
            spacer = self._composite_spacer_fluid(wash_or_spacer)
        flusher = next((fluid for fluid in fluids if fluid.role == FluidRole.FLUSHER), None)
        if mud is None or (lead is None and tail is None):
            raise ValueError("需要钻井液和至少一个水泥浆流体")
        return mud, lead, tail, spacer, flusher

    @staticmethod
    def _wash_spacer_volume_weights(
        wash_spacer_fluids: Sequence[FluidSpec],
        schedule: "PumpingSchedule | None",
    ) -> list[float] | None:
        """从泵注程序提取各 WASH/SPACER 流体的设计体积权重。

        schedule 为 None 或某流体未出现在泵注序列中时，该流体权重记 0；
        全部权重为 0 时返回 None（退化等权）。仅统计正向泵注（排量>0）步骤。
        """
        if schedule is None:
            return None
        name_by_norm = {f.name.strip(): i for i, f in enumerate(wash_spacer_fluids)}
        volumes = [0.0] * len(wash_spacer_fluids)
        for step in schedule.steps:
            idx = name_by_norm.get(step.fluid_name.strip())
            if idx is None:
                continue
            if float(step.rate_m3_min) > 0.0:
                volumes[idx] += float(step.volume_m3)
        if all(v <= 0.0 for v in volumes):
            return None
        return volumes

    @staticmethod
    def _composite_spacer_fluid(
        fluids: Sequence[FluidSpec],
        volume_fractions: Sequence[float] | None = None,
    ) -> FluidSpec:
        """把多种 WASH/SPACER 流体合成为单一等效代表流体。

        合成规则（体积加权）：
        - 密度：线性加权 ρ_mix = Σfᵢ·ρᵢ；
        - Bingham：PV/YP 线性加权（与 _compute_props 对相分数做体积加权混合的口径一致）；
        - 幂律：n 线性加权、K 对数加权（与 _compute_props 的 n_mix/kappa_mix 同口径）；
        - 组分流变模型不一致时统一映射为 Bingham：幂律/HB 折算 γ_ref=20 s⁻¹ 等效黏度
          参与线性加权（现场泵排量对应剪切速率量级），YP 线性加权。
        全部组分流变一致时保留原模型类型。
        """
        if volume_fractions is None:
            weights = np.full(len(fluids), 1.0 / len(fluids), dtype=float)
        else:
            weights = np.asarray(volume_fractions, dtype=float)
            total = float(weights.sum())
            weights = weights / total if total > 0.0 else np.full(len(fluids), 1.0 / len(fluids))
        models = {f.rheology_model for f in fluids}
        density = float(sum(w * f.density_kg_m3 for w, f in zip(weights, fluids)))
        name = "+".join(f.name for f in fluids)
        role = FluidRole.SPACER
        if len(models) == 1:
            model = models.pop()
            if model == RheologyModel.BINGHAM:
                pv = float(sum(w * (f.plastic_viscosity_pa_s or 0.0) for w, f in zip(weights, fluids)))
                yp = float(sum(w * (f.yield_stress_pa or 0.0) for w, f in zip(weights, fluids)))
                return FluidSpec(name, role, density, model, plastic_viscosity_pa_s=pv,
                                 yield_stress_pa=yp)
            if model == RheologyModel.POWER_LAW:
                n_mix = float(sum(w * (f.power_law_n or 1.0) for w, f in zip(weights, fluids)))
                log_k = float(sum(w * math.log(max(f.consistency_k or 1e-12, 1e-12))
                                  for w, f in zip(weights, fluids)))
                return FluidSpec(name, role, density, model, power_law_n=n_mix,
                                 consistency_k=math.exp(log_k))
            if model == RheologyModel.NEWTONIAN:
                pv = float(sum(w * (f.plastic_viscosity_pa_s or 0.0) for w, f in zip(weights, fluids)))
                return FluidSpec(name, role, density, model, plastic_viscosity_pa_s=pv)
            # HB：屈服+幂律参数全保留
            yp = float(sum(w * (f.yield_stress_pa or 0.0) for w, f in zip(weights, fluids)))
            n_mix = float(sum(w * (f.power_law_n or 1.0) for w, f in zip(weights, fluids)))
            log_k = float(sum(w * math.log(max(f.consistency_k or 1e-12, 1e-12))
                              for w, f in zip(weights, fluids)))
            return FluidSpec(name, role, density, model,
                             yield_stress_pa=yp, power_law_n=n_mix, consistency_k=math.exp(log_k))
        # 混合模型：统一映射为 Bingham（幂律/HB 折算 γ_ref 等效黏度线性加权，YP 线性加权）
        gamma_ref = 20.0
        eff_mu = []
        for f in fluids:
            if f.rheology_model == RheologyModel.POWER_LAW:
                eff_mu.append(float(f.consistency_k or 0.0) * gamma_ref ** (float(f.power_law_n or 1.0) - 1.0))
            elif f.rheology_model == RheologyModel.HERSCHEL_BULKLEY:
                eff_mu.append(float(f.yield_stress_pa or 0.0) / gamma_ref
                              + float(f.consistency_k or 0.0) * gamma_ref ** (float(f.power_law_n or 1.0) - 1.0))
            else:
                eff_mu.append(float(f.plastic_viscosity_pa_s or 0.0))
        pv_mix = float(sum(w * mu for w, mu in zip(weights, eff_mu)))
        yp_mix = float(sum(w * (f.yield_stress_pa or 0.0) for w, f in zip(weights, fluids)))
        return FluidSpec(name, role, density, RheologyModel.BINGHAM,
                         plastic_viscosity_pa_s=max(pv_mix, 1e-6), yield_stress_pa=yp_mix)

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
                         f_safety):
        """M3 可重启屈服门槛：每深度列以该列流动最快元（|w| 最大且 w>0）为参考，
        按平行槽流 τw=G·b/2 外推各元壁面剪应力。immobile = τw_extrap ≤ f·τy。

        关键不变量：参考元（正在流动）本身永不冻结——它在定义上可流动；只有壁面
        剪应力低于 f·τy 的更窄/更慢元才冻结。若某列完全无流动（has_flow=False），
        且水泥已到达，则整列冻结（无法外推 G）；前锋未到列不冻结。
        停泵期不调用（run() 泵注分支门控）。
        2026-09-07 精简：c_min_residual 形参删除（B2 判据不消费）。"""
        b = np.maximum(b, 1e-12)
        gamma = 6.0 * np.abs(w) / b               # w=0 → τw=0 即真实静止，不加 floor
        tau_w_field = mu_reg * gamma
        ny, nz = w.shape
        # 每列参考元：|w| 最大且 w>0；非流动元罚为 -1
        w_rank = np.where(w > 0.0, np.abs(w), -1.0)
        ref_row = np.argmax(w_rank, axis=0)            # (nz,)
        has_flow = np.any(w > 0.0, axis=0)            # (nz,)
        col = np.arange(nz)
        ref_row_safe = np.where(has_flow, ref_row, 0)
        tau_w_ref = tau_w_field[ref_row_safe, col]
        b_ref = b[ref_row_safe, col]
        G = 2.0 * tau_w_ref / np.maximum(b_ref, 1e-12)
        tau_w_extrap = G[None, :] * b / 2.0           # (ny,nz)
        # 参考元掩码：正在流动的最快元永不冻结（它确实在流，τw 判据不能冻结参考元自身）
        ref_mask = np.zeros_like(w, dtype=bool)
        ref_mask[ref_row_safe[has_flow], col[has_flow]] = True
        immobile = (tau_w_extrap <= f_safety * tau_y) & (~ref_mask) & (cement_ever > 0.0)
        # 2026-09-02 删除非物理的 residual_wall 永久浓度冻结（见下行注释），静泥层只由 immobile 决定
        wall_new = np.where(immobile, 1.0, 0.0)  # 仅可逆τw判据；残余泥膜由浓度场c<1计入ηE，不再清零速度
        # 整列无流动且水泥已到 -> 整列冻结（无法定义参考 G）
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
        w_prev: Array,
        geom: Dict[str, Array],
        mud_fluid: FluidSpec,
        lead_fluid: FluidSpec | None,
        tail_fluid: FluidSpec | None,
        spacer_fluid: FluidSpec | None,
    ) -> Tuple[Array, Array, Array, Array, Array, Array, Array, Array, Array]:
        """计算混合物系的表观粘度、密度、钻井液分数、混合屈服应力、黏度比 m 场、相黏度场（η1=泥浆相, η2=水泥相）及混合幂律参数（n_mix, kappa_mix）。
        2026-09-07 精简：flusher 相降级——8 井生产口径入库 flusher 恒为 0（hu102/ht1_004
        实测），2D 不再为 FLUSHER 建独立浓度场；泥浆由 1−lead−tail−spacer 闭合。"""
        # 四相体积分数闭合：显式跟踪领浆、尾浆、前置/隔离液，钻井液由守恒关系反算。
        mud = np.clip(1.0 - lead - tail - spacer, 0.0, 1.0)
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

    def _buoyancy_force_vector(self, geom: Dict[str, Array], beta_deg: Array | float,
                               f2: float) -> Tuple[Array, Array]:
        """计算论文式 (2.5b) 的浮力体力向量 f = (r_a·cosβ, r_a·sin(πφ)·sinβ)/F²。

        R2 (I3 通量) 与 R3 (真体力) 共用。**F² 必须由调用方按 Z&F22 (2.6) 现算传入**
        （见 `_froude_squared_at`，物理量级 O(10⁻²)）：本方法原先把 ``F2 = 1.0`` 写死在
        内部并声称"真 F 校正在 _compute_velocity 内做"，而那次校正从未存在，导致
        方位浮力修正幅度只剩 ~4×10⁻⁵——浮力在动力学里实际缺席。2026-09-14 Task 4 废止。

        Args:
            geom: 几何字典（用 ``phi``；``hole_mm``/``od_mm`` 推平均半径 r_a）。
            beta_deg: 井斜角 β，度。
            f2: Froude 数平方 F²（无量纲，Z&F22 (2.6)）；非正/极小值被夹到 1e-12 防除零。

        Returns:
            (f_phi, f_xi)，shape 与 ``geom['phi']`` 广播兼容，即 (ny, nz)。
        """
        phi = geom["phi"][:, None]  # (ny, 1)
        beta_rad = np.deg2rad(np.asarray(beta_deg, dtype=float))
        # 平均半径 r_a（用 mean_radius 近似，从 geom 取 hole/od 推算）
        hole_mm = geom.get("hole_mm", np.full((1, self.nz), 220.0))
        od_mm = geom.get("od_mm", np.full((1, self.nz), 139.7))
        # 沿深度取均值半径（米），广播到 (ny, nz)
        r_a_m = np.mean((hole_mm + od_mm) / 4.0) / 1000.0
        # Task 4: 1/F² 由外部按 (2.6) 传入（原先硬编码 F2 = 1.0）；下限防除零。
        f2_safe = max(float(f2), 1.0e-12)
        f_phi = (r_a_m / f2_safe) * np.sin(np.pi * phi) * np.sin(beta_rad)  # (ny,1) 广播
        f_xi = np.full_like(f_phi, (r_a_m / f2_safe) * np.cos(beta_rad))
        # 广播到 (ny, nz)
        f_phi = np.broadcast_to(f_phi, (self.ny, self.nz)).astype(float, copy=True)
        f_xi = np.broadcast_to(f_xi, (self.ny, self.nz)).astype(float, copy=True)
        return f_phi, f_xi

    def _froude_squared_at(self, geom: Dict[str, Array], q_m3s: float, w_field: Array,
                           mud_fluid: FluidSpec) -> float:
        """按 Z&F22 (2.6) 现算当前时间步的 F²（式 2.5b 浮力体力向量的标定分母）。

        ``F = √(τ̂₀/(ρ̂₁·ĝ·δ₀·r̂ₐ*))`` ⇒ ``F² = τ̂₀/(ρ̂₁·ĝ·δ₀·r̂ₐ*)``，其中
        ``τ̂₀ = μ̂₁·ŵ₀/d̂`` 为**被顶替液（钻井液）**中的黏性应力尺度。

        取值口径（与 Task 3 的 :mod:`cemdisp.models2d.buoyancy`、summary 段一致）：

        - ``μ̂₁ = buoyancy.fluid_apparent_viscosity(mud, γ̇)``，剪切率 ``γ̇ = 6|w|/b``
          与 `_compute_props` 同约定；单位 Pa·s。
        - ``ŵ₀ = q/A``：截面平均轴向速度（**不是** ``mean(|w|)``），m/s。
        - ``ρ̂₁ = mud_fluid.density_kg_m3``（被顶替液密度），kg/m³。
        - ``d̂ = mean(geom["H"])``：Z&F22 半间隙（``= (r_o−r_i)/2 = (井径−外径)/4``），m。
        - ``r̂ₐ* = mean((hole+od)/4)/1000``：沿程平均半径，m（论文中 ``r̂ₐ*`` 沿环空流道平均）。
        - ``δ₀ = d̂/r̂ₐ*``：**无量纲**参考间隙比。论文 (2.1) 推导处把窄间隙参数写作
          ``δ = d̂/(πr̂ₐ*)``；(2.6) 的 ``δ₀`` 是同一量的参考值，本实现由
          ``δ₀·r̂ₐ* = d̂`` 钉定。该取值的依据：把 (2.5b)/(2.6) 的
          ``|b| ≈ (ρ−1)/F²``（论文 p.8"b 即浮力向量的大小"）与论文 p.8 的浮力数
          ``b = Δρ·ĝ·d̂²/(μ̂₁ŵ₀)`` 联立，得 ``F²·b = Δρ/ρ̂₁``（Atwood 数），
          即 ``F² = μ̂₁ŵ₀/(ρ̂₁·ĝ·d̂²)``，等价于 ``δ₀·r̂ₐ* = d̂``。

        传给 `buoyancy.froude_squared` 时 ``gap_scale_m = half_gap_m / mean_radius_m``
        ——二者是**不同的量**（前者无量纲间隙比、后者长度），但乘积恰为 ``d̂``。

        两个调用点（`_compute_velocity` 的 R3 真体力、run 循环的 R2 I3 弥散通量）
        各自按本方法现算，同一时间步内几何/物性口径一致。
        """
        half_gap_m = float(np.mean(geom["H"])) if "H" in geom else 0.5 * float(np.mean(geom["b"]))
        # 退化口径：合成几何的单元测试可能只给 y/phi/b（无 hole/od/H）。
        # 因本式只以乘积 ``δ₀·r̂ₐ* = d̂`` 起作用，此时取 r̂ₐ* = d̂（即 δ₀ = 1）不影响 F²；
        # 环形截面积改取模型自身的 ``∫ b dy``（与体积 scale 依赖的
        # ``∫∫ b dy ds = 物理环空体积`` 恒等式同源），仍保持 ŵ₀ = q/A 的口径。
        if "hole_mm" in geom and "od_mm" in geom:
            mean_radius_m = float(np.mean((geom["hole_mm"] + geom["od_mm"]) / 4.0)) / 1000.0
            annulus_area_m2 = float(
                np.mean(np.pi / 4.0 * ((geom["hole_mm"] / 1000.0) ** 2
                                       - (geom["od_mm"] / 1000.0) ** 2))
            )
        else:
            mean_radius_m = half_gap_m
            annulus_area_m2 = float(np.trapezoid(np.mean(geom["b"], axis=1), x=geom["y"]))
        if float(q_m3s) > 0.0 and annulus_area_m2 > 0.0:
            w0_mps = float(q_m3s) / annulus_area_m2
        else:
            # 退让口径：无泵注排量（异常输入）时回落到末步速度场均值
            w0_mps = float(np.mean(np.abs(w_field)))
        # 泥浆表观黏度的剪切率约定与 _compute_props 一致：γ̇ = 6|w|/b
        shear_rate_mud = (6.0 * float(np.mean(np.abs(w_field)))
                          / max(float(np.mean(geom["b"])), 1e-12))
        return buoyancy.froude_squared(
            mu_displaced=buoyancy.fluid_apparent_viscosity(mud_fluid, shear_rate_mud),
            w0_mps=w0_mps,
            half_gap_m=half_gap_m,
            rho_displaced=float(mud_fluid.density_kg_m3),
            gap_scale_m=half_gap_m / max(mean_radius_m, 1.0e-12),
            mean_radius_m=mean_radius_m,
        )

    def _compute_buoyancy_number(self, rho_displacing_kg_m3: float, rho_displaced_kg_m3: float,
                                 gap_m: float, mu_displaced_pa_s: float, velocity_m_s: float) -> float:
        """无量纲浮力数 b 的薄委托（2026-09-14 Task 3 起）。

        公式与口径已统一搬到 :mod:`cemdisp.models2d.buoyancy`（Z&F22 p.8）：
        ``b = (ρ_displacing − ρ_displaced)·g·d²/(μ_displaced·w₀)``，d 为半间隙。

        本方法仅为兼容既有调用方/测试保留（签名不变，语义不变）：
        ``gap_m`` 是全间隙，内部按半间隙 ``gap_m/2`` 传入。
        新代码请直接调用 ``buoyancy.buoyancy_number``。
        """
        return buoyancy.buoyancy_number(
            rho_displacing=rho_displacing_kg_m3,
            rho_displaced=rho_displaced_kg_m3,
            half_gap_m=gap_m / 2.0,
            mu_displaced=mu_displaced_pa_s,
            w0_mps=velocity_m_s,
        )

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
        wall: Array | None = None,
    ) -> Tuple[Array, Array, Array, Array, Array, Array, Array, Array, Array, Array, Array, Array]:
        """计算环空速度场（论文D2DGA口径）。

        采用 Zhang & Frigaard (2022) 的Hele-Shaw风格速度场：
        1. 计算局部混合流体的表观粘度与密度；
        2. 以 ``b²/μ`` 构造偏心通道主导局部流动度；
        3. 根据密度差（顶替液 vs 被顶替液）计算浮力稳定系数；
        4. 用浮力修正项调整宽边/窄边速度分配；
        5. 由截面排量约束得到轴向速度 ``w``。

        2026-09-07 精简：flusher 形参删除（相降级，见 _compute_props 注）。

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
        # 2026-09-07 R0 分支删除：m 恒由 _compute_props 自动计算（enable_d2dga_auto_m
        # 恒 True，标量 d2dga_viscosity_ratio 路径为旧论文 R0 状态已移除）。
        m_local = float(np.mean(m_field))
        i1_base = d2dga_dispersion_I1(c_bar, m_local)
        if self.enable_power_law_gap_law:
            # 2026-09-06 幂律缝隙律修正：层流偏心环空各通道流量份额
            # q ∝ b^(2+1/n)·G^(1/n)（Walton & Bittleston 1991 JFM 222:39-60 窄隙
            #   Bingham/幂律槽流；Pelipenko & Frigaard 2004c JFM 520:343-377；
            #   Maleki & Frigaard 2017 式60-61 闭式）。
            # 速度场取缝隙平均速度口径 w = q/b ∝ b^(1+1/n)：
            #   牛顿 n=1 → w ∝ b²、通量 w·b ∝ b³（Poiseuille，与旧口径逐位一致）；
            #   剪切变稀 n<1 → 指数 1+1/n > 2，窄边分流比 b³ 更极端（b³ 低估通道化）。
            # 压降梯度项 G^(1/n) 全截面同值，被归一化分母吸收，不影响方位分配。
            # 混合物 n 用 _compute_props 的 n_mix 场（体积分数加权，Bingham→n=1）。
            n_safe = np.clip(n_mix, 0.25, 1.5)
            gap_exponent = 1.0 + 1.0 / n_safe
            base = (b / np.maximum(b_mean, 1.0e-12)) ** gap_exponent / np.maximum(eta_mix, 1.0e-9)
        else:
            # 旧口径：Hele-Shaw b² 流动度（逐位复现基线）
            base = (b / np.maximum(b_mean, 1.0e-12)) ** 2 / np.maximum(eta_mix, 1.0e-9)
        base = base * i1_base  # I₁ 乘子

        # === 速度场流动度：偏心通道主导 + 浮力修正 ===
        # density_contrast > 0 表示顶替液更重（水泥重 vs 泥浆轻），有助于窄边推进
        # density_contrast < 0 表示顶替液更轻，加剧宽边窜流
        # Task 3: 顶替液密度统一走 buoyancy.displacing_density_kg_m3（全仓唯一口径，
        # 0.67×领浆 + 0.33×尾浆），消除与 summary 段浮力数口径不一致导致的 b 符号翻转。
        rho_disp = buoyancy.displacing_density_kg_m3(lead_fluid, tail_fluid, mud_fluid)

        phi = geom["phi"][:, None]
        ebar = geom["e"][None, :]
        if self.enable_true_buoyancy and self.enable_d2dga:
            # T1-3b: 体力向量注入流动度（式 2.5b/4.24），替换 (2φ−1) 简化代理
            beta_deg_local = float(np.mean(geom.get("inc_deg", np.zeros(self.nz))))
            # Task 4: F² 按 Z&F22 (2.6) 现算（原先内部硬编码 1.0 → 浮力缺席动力学），
            # 量级 O(10⁻²)；ŵ₀ = q/A、μ̂₁ 取泥浆、d̂ = mean(geom["H"])。
            f2_local = self._froude_squared_at(geom, q_m3s, w_prev, mud_fluid)
            f_phi_arr, _ = self._buoyancy_force_vector(geom, beta_deg_local, f2_local)
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

        # 由截面排量约束得到轴向速度 w（2026-09-02 守恒修正：截面权重去掉多余因子2，
        # 使同心极限 w=Q/A、全环空通量 2∫w b dy=Q；旧口径仅输运 Q/2 导致 η_E 系统偏低）。
        dy = np.gradient(geom["y"])[:, None]
        if self.enable_regime_split:
            # M2: 局部流态修正固定点迭代（Maleki & Frigaard 2017 式58-66）
            # 浓度相关量（base/buoyancy_shape/wall/黏度/密度/b）在迭代外缓存；
            # 黏度场保持 w_prev 一步滞后（既有约定），迭代只重算 Re_p/R/pref/area_weight/w。
            # rho_kg_m3 与 kappa_mix(Pa·s^n)/tau_y(Pa) 同单位系，保证 Re_p/He 无量纲正确。
            from cemdisp.models2d import regime_closure as rc
            wall_factor = np.ones_like(base) if wall is None else (1.0 - wall)
            he = rc.hedstrom_number(tau_y, rho_kg_m3, n_mix, kappa_mix, b)
            # ⚠️ 未标定（临时公式）：re_crit = 2100(1+0.1·He) 为 provisional 标定钮（屈服推迟转捩）。
            # 启用 enable_regime_split 前须三重回归锚定（Walton&Bittleston / Z22 Table3 / Foolad 定性）。
            # 对现 8 井中性（Bingham 泥浆 He 高 → 全层流 R=1，此式不影响基线结果）。
            re_crit = 2100.0 * (1.0 + 0.1 * he)
            w_k = w_prev.copy()
            R = np.ones_like(w_k)
            for _ in range(self.regime_max_iter):
                re_p = rc.metzner_reed_re(w_k, rho_kg_m3, n_mix, kappa_mix, b)
                R_new, _ = rc.drag_weight(re_p, he, n_mix, re_crit, self.regime_re_turb_ratio)
                pref_k = np.maximum(base * buoyancy_shape * R_new, 1.0e-8) * wall_factor
                area_w = np.sum(pref_k * b * dy, axis=0, keepdims=True)
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
            # 使全环空通量 2·Σw·b·dy = 2·q_half = Q 精确成立（2026-09-02 守恒修正口径）。
            # 消除欠松弛迭代返回 w_k 时 ~1.6% 的瞬态守恒误差。
            pref_final = np.maximum(base * buoyancy_shape * R, 1.0e-8) * wall_factor
            area_final = np.sum(pref_final * b * dy, axis=0, keepdims=True)
            w = q_half * pref_final / np.maximum(area_final, 1.0e-12)
            # area_weight 与最终 w 构成一致配对（同一 area_final）
            area_weight = area_final
        else:
            area_weight = np.sum(pref * b * dy, axis=0, keepdims=True)
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

    def _depth_profiles(self, geom: Dict[str, Array], lead: Array, tail: Array, spacer: Array) -> pd.DataFrame:
        """计算深度方向的平均剖面数据。2026-09-07 精简：flusher 相降级，
        "冲洗液平均浓度"列保留恒 0（下游 CSV 列位稳定）。"""
        cement = np.clip(lead + tail, 0.0, 1.0)
        eff = cement
        mud = np.clip(1.0 - lead - tail - spacer, 0.0, 1.0)
        return pd.DataFrame(
            {
                "井深_m": geom["md"],
                "领浆平均浓度": np.average(lead, axis=0, weights=geom["b"]),
                "尾浆平均浓度": np.average(tail, axis=0, weights=geom["b"]),
                "水泥平均浓度": np.average(cement, axis=0, weights=geom["b"]),
                "前置液隔离液平均浓度": np.average(spacer, axis=0, weights=geom["b"]),
                "冲洗液平均浓度": np.zeros(geom["md"].shape),  # flusher 相已降级，恒 0 保列位
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
        schedule: "PumpingSchedule | None" = None,
    ) -> AnnulusSimulationResult:
        """运行论文口径的环空二维顶替求解。

        Args:
            well_spec: 井身结构。
            fluids: 流体序列（mud/lead/tail/spacer/flusher 等）。
            inlet_state_provider: 环空入口边界状态提供器（1D-2D 耦合时由
                ``build_coupled_annulus_inlet_provider`` 构造）。
            schedule: 泵注程序（可选，默认 None）。仅用于末尾 Tier0 诊断聚合：
                提供后 T0-6 停泵有限时间衰减诊断（shutdown_decay）可用；
                为 None 时诊断层优雅降级（记 notes "未提供 schedule"），
                不影响求解结果与既有调用方（向后兼容）。
        """

        mud_fluid, lead_fluid, tail_fluid, spacer_fluid, flusher_fluid = self._pick_fluids(fluids)
        # 2026-09-06 选相修复：多种 WASH/SPACER 并存（如平衡液+驱油隔离液）时，
        # 按泵注程序中各流体的设计体积加权重建等效代表流体——进入环空的 spacer 相
        # 由这些流体按入库体积混合而成，物性（密度/黏度/屈服）应取入库加权而非
        # "第一个 WASH/SPACER"（旧口径使真实隔离液物性在 2D 闭包中失效）。
        _wash_spacer_fluids = [f for f in fluids
                               if f.role in {FluidRole.WASH, FluidRole.SPACER}]
        if len(_wash_spacer_fluids) > 1:
            _ws_weights = self._wash_spacer_volume_weights(_wash_spacer_fluids, schedule)
            spacer_fluid = self._composite_spacer_fluid(_wash_spacer_fluids, _ws_weights)
        # 诊断暴露：最近一次 run 实际使用的等效隔离液（rerun/报告脚本读取）
        self._active_spacer_fluid = spacer_fluid
        geom = self._build_geom(well_spec)
        lead = np.zeros((self.ny, self.nz), dtype=float)
        tail = np.zeros((self.ny, self.nz), dtype=float)
        spacer = np.zeros((self.ny, self.nz), dtype=float)
        # 2026-09-07 精简：flusher 相降级，不再建独立浓度场（入库恒 0，见 _compute_props 注）
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
        wall_snapshots: list[Array] = []
        snapshot_times: list[float] = []
        cumulative_lead_in_m3 = 0.0
        cumulative_tail_in_m3 = 0.0
        cumulative_spacer_in_m3 = 0.0
        # Task 3: 最近一次有效泵注排量（m³/s），供 summary 段按截面平均速度口径
        # 计算浮力数 b 的 w₀ = q/A（不再用"最后一步"速度场均值）。
        last_pump_rate_m3s = 0.0

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

            # 泵停判断：排量低于阈值时认为泵已停止。
            # 泵停后水泥场冻结——不再平流、扩散或壁面清除，
            # 因为水泥静凝胶强度在短时间（<2h）内足以抵抗浮力滑塌。
            pump_active = inlet_state.flow_rate_m3_s > 1.0e-9

            if pump_active:
                last_pump_rate_m3s = float(inlet_state.flow_rate_m3_s)  # Task 3: w₀ = q/A 用
                # === 正常泵注阶段：仅执行论文口径核心平流 + D2DGA 通量修正 ===
                w, v, mu, rho, mud, Re, mu_turbulent, m_field, _tau_y, _eta2, _n_mix, _kappa_mix = self._compute_velocity(
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
                    # 2026-09-07 R0 分支删除：m 恒用 auto-m 场
                    f_amp = d2dga_flux_amplification(cement, m_field)
                else:
                    f_amp = 1.0

                w_d2dga = w * f_amp
                v_d2dga = v * f_amp
                ysrc = ygrid - v_d2dga * dt_step
                ssrc = sgrid - w_d2dga * dt_step
                lead_adv = _bilinear_interp(lead, ysrc, ssrc, geom, inlet_lead_fraction)
                tail_adv = _bilinear_interp(tail, ysrc, ssrc, geom, inlet_tail_fraction)
                spacer_adv = _bilinear_interp(spacer, ysrc, ssrc, geom, inlet_spacer_fraction)
                lead = np.clip(lead_adv, 0.0, 1.0)
                tail = np.clip(tail_adv, 0.0, 1.0)
                spacer = np.clip(spacer_adv, 0.0, 1.0)
                # 数值扩散可能使显式相之和略超1；按比例压回可行域，保持泥浆分数非负。
                tracked_total = lead + tail + spacer  # 四相过填修正（flusher 相已降级）
                overfilled = tracked_total > 1.0
                lead[overfilled] /= tracked_total[overfilled]
                tail[overfilled] /= tracked_total[overfilled]
                spacer[overfilled] /= tracked_total[overfilled]

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
                # spacer 基础弥散系数为 0.012/0.012（轴向/方位角同值，独立于 lead/tail）。
                # 必须用字面量 0.012，而非 0.018*0.667（=0.012006）或 0.015*0.8：
                # 默认(fixed dt=4, scale=1)下要求 0.012*1.0==0.012 与基线硬编码逐位复现。
                _ax_sf = 0.012 * _dt_norm
                _az_sf = 0.012 * _dt_norm
                lead = self._smooth_dispersion(lead, axial=_ax, azimuthal=_az)
                tail = self._smooth_dispersion(tail, axial=_ax, azimuthal=_az)
                spacer = self._smooth_dispersion(spacer, axial=_ax_sf, azimuthal=_az_sf)
                # T1-6: 弥散后再次执行四相过填修正，防止 _smooth_dispersion 数值扩散
                # 使 lead+tail+spacer 再次超过 1，破坏体积分数闭合。
                tracked_total = lead + tail + spacer
                overfilled = tracked_total > 1.0
                lead[overfilled] /= tracked_total[overfilled]
                tail[overfilled] /= tracked_total[overfilled]
                spacer[overfilled] /= tracked_total[overfilled]

                # R2: I3 浮力弥散通量（式 4.25 第二项）—— 仅作用于水泥相(lead+tail)
                if self.enable_d2dga_i3_flux and self.enable_d2dga:
                    cement_for_flux = np.clip(lead + tail, 0.0, 1.0)
                    # 浮力向量 f（用当前井段平均井斜）
                    beta_deg_local = float(np.mean(geom["inc_deg"])) if "inc_deg" in geom else 0.0
                    # Task 4: F² 按 Z&F22 (2.6) 现算，与 _compute_velocity 内同一口径
                    # （排量取最近一次有效泵注 q，速度取本步 w），量级 O(10⁻²)。
                    f2_local = self._froude_squared_at(
                        geom, last_pump_rate_m3s, w_prev, mud_fluid)
                    f_phi_arr, f_xi_arr = self._buoyancy_force_vector(
                        geom, beta_deg_local, f2_local)
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
                    q_phi, q_xi = d2dga_buoyancy_flux(
                        cement_for_flux, m_field, delta_rho, H_field, eta2,
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
                # 这里按入口累计体积对领浆、尾浆、前置液/隔离液分别做体积上限约束。
                cumulative_lead_in_m3 += inlet_state.flow_rate_m3_s * inlet_lead_fraction * dt_step
                cumulative_tail_in_m3 += inlet_state.flow_rate_m3_s * inlet_tail_fraction * dt_step
                cumulative_spacer_in_m3 += inlet_state.flow_rate_m3_s * inlet_spacer_fraction * dt_step
                lead = _limit_phase_volume(lead, geom, cumulative_lead_in_m3, self.open_outlet)
                tail = _limit_phase_volume(tail, geom, cumulative_tail_in_m3, self.open_outlet)
                spacer = _limit_phase_volume(spacer, geom, cumulative_spacer_in_m3, self.open_outlet)

                # T1-5/M3: 壁面冻结层（Bararpour 2025）。
                # 2026-09-07 精简：enable_yield_gate=False 的 c_min 浓度兜底轨已删除
                # （09-02 取证 B2 物理屈服门取代之；yield_gate_c_min_residual 在
                # B2 分支本就不消费）。关闭屈服门 = wall 恒零（无壁面静止层）。
                cement_local = np.clip(lead + tail, 0.0, 1.0)
                cement_ever = np.maximum(cement_ever, cement_local)
                if self.enable_yield_gate:
                    wall = self._yield_gate_wall(
                        w, geom["effective_b"], mu, _tau_y, cement_ever, cement_local,
                        self.yield_gate_f_safety)

            else:
                # === 泵停阶段：冻结浓度场，仅记录指标 ===
                w, v, mu, rho, mud, Re, mu_turbulent, m_field, _tau_y, _eta2, _n_mix, _kappa_mix = self._compute_velocity(
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
                wall_snapshots.append(wall.copy())
                snapshot_times.append(record_time)

            cement = np.clip(lead + tail, 0.0, 1.0)
            # η_E ≡ 水泥库存占据率恒等式（η_E = bulk_fill = 域内水泥/域满，09-06 取证
            # 8 井偏差 ≤2.1e-4）。两列保留（下游 30+ 处消费两键），注释声明即可。
            eff = cement
            bulk_fill = _trapez2d(geom["b"] * cement, geom) / half_volume
            effective_efficiency = bulk_fill  # 同值列（恒等式），不再重复全场积分

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
                    0.0,  # mean_flusher：flusher 相已降级，恒 0 保列位
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
        depth_profiles = self._depth_profiles(geom, lead, tail, spacer)
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

        # R3 → Task 3: 无量纲浮力数 b 统一走 cemdisp.models2d.buoyancy（全仓唯一口径，
        # 文献锚点 Z&F22 p.8）。修复三处旧口径缺陷：
        #   ① 顶替液密度双口径：此处旧版只用领浆，而 _compute_velocity 用
        #      0.67×领浆+0.33×尾浆 → b 可能符号翻转；现统一 displacing_density_kg_m3。
        #   ② 幂律/HB 泥浆被 `plastic_viscosity_pa_s or 0.05` 静默回退到 0.05 Pa·s
        #      → 改用 fluid_apparent_viscosity（K·γ̇^(n−1)，缺参数抛错不回退）。
        #   ③ w₀ 旧取"最后一步"速度场均值 → 改为截面平均速度 q/A（末态泵注排量/环形截面积）。
        rho_displacing = buoyancy.displacing_density_kg_m3(lead_fluid, tail_fluid, mud_fluid)
        # 泥浆表观黏度的剪切率约定与 _compute_props 一致：γ̇ = 6|w|/b
        shear_rate_mud = (
            6.0 * float(np.mean(np.abs(w_prev))) / max(float(np.mean(geom["b"])), 1e-12)
        )
        mu_displaced = buoyancy.fluid_apparent_viscosity(mud_fluid, shear_rate_mud)
        # 几何半间隙 d̂ = (r_o−r_i)/2 = (hole−od)/4，直接取模型自己的半间隙场
        # geom["H"]（`_build_geom` 中 H=2b 一半，与输运用的 geom["b"]=2d̂ 严格自洽：
        # b=2H 逐格成立，故 mean(b)/2 ≡ mean(H)）。用 hole/od 手推 (hole−od)/4 亦可，
        # 但会引入 ≤2% 的体积 scale 残差差异，故以 geom["H"] 为准。
        half_gap_m = float(np.mean(geom["H"]))
        # 截面平均轴向速度 w₀ = q/A（环形截面积由 hole/od 算）
        annulus_area_m2 = float(
            np.mean(np.pi / 4.0 * ((geom["hole_mm"] / 1000.0) ** 2
                                   - (geom["od_mm"] / 1000.0) ** 2))
        )
        if last_pump_rate_m3s > 0.0 and annulus_area_m2 > 0.0:
            w0_mps = last_pump_rate_m3s / annulus_area_m2
        else:
            # 退让口径：全程无泵注排量记录（异常输入）时回落到末步速度场均值
            w0_mps = float(np.mean(np.abs(w_prev)))
        b_number = buoyancy.buoyancy_number(
            rho_displacing=rho_displacing,
            rho_displaced=mud_fluid.density_kg_m3,
            half_gap_m=half_gap_m,
            mu_displaced=mu_displaced,
            w0_mps=w0_mps,
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
            # 2026-09-07 精简：flusher 相降级，不再赋值（默认 None/()，向后兼容）
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
            tier0 = compute_all_tier0_diagnostics(
                result, fluids=fluids, well_spec=well_spec, schedule=schedule
            )
            result.summary["tier0_diagnostics"] = tier0.to_dict()  # type: ignore[index]
        except Exception as exc:
            result.summary["tier0_diagnostics_error"] = f"{type(exc).__name__}: {exc}"

        return result
