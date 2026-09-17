"""
多井通用环空二维 D2DGA 求解器。

本模块实现面向 ``cemdisp`` 主链路的偏心环空二维 D2DGA 核心，口径按
Zhang & Frigaard (2022)（JFM 947 A32）源模型重构完成（2026-09-14/15，Task 0–14）。
**方法学边界声明（半环空对称条件、适用域、浮力数区间、弥散口径等六条）与
式号→代码锚点对照表见 ``docs/源模型口径与适用域声明.md``**。核心口径：

1. 偏心窄环空几何展开；求解域为**半环空**（φ∈[0,1]，φ=0 低边）——Z&F22
   p.25/32 "the imposition of symmetry"：人为施加的 Ψ 对称条件，非物理必然；
2. 速度场双路径（Task 9，2026-09-15）：``enable_stream_function=True``
   （默认）由 (4.22) 流函数椭圆方程解 + (2.2) 换算给出 (w, v)，浮力完整经
   b = (ρ − Δρ·I₂/(H·I₁))·f 向量进入（(4.22)/(4.13)/(2.5b) 字面分组）；
   ``False`` 为旧代数流动度路径（`_mobility_profile` + 截面归一，可回退）；
3. 浮力两阶接线（Task 4/5）：O(ε) 方位浮力经 F² 按 (2.6) 现算真实进入
   （八井 F² ∈ [1.2e-3, 3.1e-2]，替换写死的 F2=1.0）；领先阶轴向浮力在
   新路径经 (4.22) b 向量通道 A（ρ·f）+ 通道 B（−Δρ·I₂/(H·I₁)·f）完整进入，
   旧路径经 K_AXIAL·b_num·cosβ·(I₂/I₁)（H⁰ 形状，R28 裁定，provisional 1/𝒢）；
4. I₃ 浮力通量：弥散由 q₀ + I₃ 分层通量闭合承载
   （式 4.25/4.26/4.28）；本模型**无人工扩散项**（Z&F22 p.11
   "we have no diffusive terms"，2026-09-14 Task 7 删除自创拉普拉斯弥散）；
5. e = 1−standoff ∈ [0,1) **无截断**（Pelipenko04 (2.1) 文献口径；Task 10
   移除 e_clip 硬截断 0.55/0.90，弃用警告见 `__init__`）；
6. 屈服门（M3）连续化（Task 11，Pelipenko04 (2.6)-(2.8) 停流判据连续近似
   wall = clip(1−τw/(f·τy), 0, 1)，`_yield_gate_wall`）：新旧路径统一计算，
   旧**代数路径**的 mobility 直接消费 wall；**流函数新路径**（Task 9）的
   两层闭包自身无屈服项，但可经 ``enable_stream_yield_gate=True``（B-2，
   默认关）把 wall 送进 ``solve_stream_function`` 的算子（见
   `_velocity_stream_function` 旧路径差异声明）；
7. 仅输出求解域内的顶替效率与浓度场。
8. HB 闭包接线（A-3b，2026-09-17 Task 6，opt-in）：``enable_hb_closure=True``
   时 (4.22) 算子经 ``solve_stream_function_nonlinear`` 非线性外迭代（B&F25 HB
   两层闭包，物理参数 + ``velocity_scale = q_half/π``，R-T5-1 REVISED 推荐路径；
   不收敛 ⇒ 显式告警 + 回退牛顿线性闭包，R-T5-3）；``hb_fix_cement_tau_y=True``
   （R-T6-1）把水泥 τy 注入两处（屈服门混合 τy 场 + 闭包 τ_Y2，同源同值，
   ``cement_tau_y_by_role`` 显式给值）。双开关默认全关 ⇒ 逐位 = HEAD。
   Phase A 边界：只改质量流动度 I₁/I₂，通量函数 q₀ 仍走牛顿闭式（见
   ``docs/源模型口径与适用域声明.md`` §2.5）。

出口边界条件：
- 开放出口（open_outlet=True，默认）：允许水泥浆流出求解域到重叠段，适用于只模拟裸眼段；
- 封闭出口（open_outlet=False）：按累计入环空体积限制场量，适用于模拟整个环空。

注意：
- 本模块不再把泥饼、温度、凝胶强度、湍流修正等工程扩展项作为核心求解的一部分；
- 为兼顾下游脚本兼容性，旧参数与旧快照字段仍保留接口或占位输出，但不再影响求解结果。
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Callable, Dict, Mapping, Sequence, Tuple

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
from cemdisp.models2d.stream_function import (
    _VELOCITY_COMPONENT_RATIO,
    solve_stream_function,
    solve_stream_function_nonlinear,
    velocity_from_stream_function,
)
from cemdisp.models2d.hb_closure import HBClosure
from cemdisp.models2d.two_layer import hb_groups, mobility_i1, mobility_i2

if TYPE_CHECKING:  # 仅类型注解，运行时不引入 data.pumping_schedule 依赖
    from cemdisp.data.pumping_schedule import PumpingSchedule


Array = NDArray[np.float64]

# Task 5 (R2/R27/R28，provisional)：轴向浮力修正的幅值系数 K_AXIAL = 1/𝒢。
# 𝒢 为 (4.6) 无量纲修正压力梯度的参考量级：把 ∂p/∂ξ 用 τ̂₀/d̂ = μ̂₁ŵ₀/d̂² 无量纲化后，
# 牛顿参考缝隙流（半间隙 d̂、平均速度 ŵ₀ 标度）w̄ = Ĝd̂²/(3μ̂₁) 给出 𝒢 = Ĝd̂²/(μ̂₁ŵ₀) = 3。
# 即 (4.14)/(4.22) 分层浮力项相对压力驱动项的相对修正 = b_num·cosβ·(I₂/I₁)/𝒢（H⁰ 读数，R28）。
# 完整推导与 H⁰ 裁定的三重证据链见 `_mobility_profile` docstring；最终取值由 R27 八井扫描
# 与 Task 12 基准算例裁定，**不得用 clip 兜底**。
# ⚠️ T9（2026-09-15）：本系数仅旧代数路径（enable_stream_function=False）消费；
# 新路径浮力经 (4.22) b 向量完整进入（通道 A+B），再叠加 K_AXIAL 即通道 B 双重计入
# （推导见 `_velocity_stream_function` docstring）。
K_AXIAL = 1.0 / 3.0


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
    1. 平流输运：水泥浆与前置/隔离液随平均速度场向下游运移（q₀ 平均通量层）；
    2. 分层通量闭合：q₀ 平均平流 + I₃ 浮力通量（式 4.25 第二项 / (4.26)）。
       Z&F22 p.11 明言 "we have no diffusive terms"——原模型无任何人工扩散项，
       本模块同样**无人工拉普拉斯弥散**（自创弥散函数已于 2026-09-14 Task 7 删除）；
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
        # ⚠️ 2026-09-14 Task 7 弃用形参（默认 None）：自创拉普拉斯弥散已删除。
        dispersion_axial: float | None = None,
        dispersion_azimuthal: float | None = None,
        dispersion_dt_ref: float | None = None,
        dispersion_dt_scale: float | None = None,
        enable_e_clip_ruling: bool = True,
        e_clip_measured_max: float = 0.90,
        enable_power_law_gap_law: bool = True,
        enable_stream_yield_gate: bool = False,  # B-2 opt-in：屈服门进流函数算子（默认关=HEAD 逐位）
        enable_power_law_gap_correction: bool = False,  # B-3 opt-in：幂律间隙一阶修正（默认关）
        enable_stream_function: bool = True,
        # ⚠️ 2026-09-17 A-3b（Task 6）：HB 闭包接线双开关（opt-in，默认全关 ⇒ 逐位=HEAD）。
        enable_hb_closure: bool = False,
        hb_fix_cement_tau_y: bool = False,
        cement_tau_y_by_role: Mapping[str, float] | None = None,
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
            e_clip_max: ⚠️ 已弃用（2026-09-15 Task 10），默认 0.55。e_clip 硬截断
                已移除，e = clip(1−standoff, 1e-6, 1−1e-6) 按 Pelipenko04 (2.1)
                文献口径 e∈[0,1) 直取（standoff 0.35 与 0.45 不再同入 0.55 死区，
                强偏心不再封顶）。保留形参仅为既有 runner 兼容；偏离 legacy 默认
                （0.55）的传值触发一次性 DeprecationWarning 且不再生效。
            e_clip_measured_max: ⚠️ 已弃用（2026-09-15 Task 10），语义同 e_clip_max
                （legacy 默认 0.90），偏离传值触发 DeprecationWarning 且不再生效。
            enable_e_clip_ruling: ⚠️ 已弃用（2026-09-15 Task 10），语义同 e_clip_max
                （legacy 默认 True），偏离传值触发 DeprecationWarning 且不再生效。
                2026-09-06 裁定（实测井放开 0.90）随截断一同退役——文献口径无
                "按数据来源选上限"概念。
            enable_yield_gate: M3 屈服门槛开关，默认 True（2026-09-02 起物理屈服门）。
                True 时 pump 分支用 _yield_gate_wall 重建壁面冻结层
                （2026-09-15 Task 11 起为连续冻结度 wall ∈ [0,1]，非二值）；
                False（2026-09-07 起无 c_min 兜底轨）= wall 恒零（无壁面静止层）。
            yield_gate_f_safety: 屈服门槛安全系数 f，默认 1.15。
                连续冻结度判据（Pelipenko04 (2.6)-(2.8) 停流区判据的连续近似）：
                wall = clip(1 − τw_extrap/(f·τy), 0, 1)；τw ≥ f·τy 可流动
                （wall=0），τw→0 全冻（wall→1）。
            dispersion_axial: ⚠️ 已弃用（2026-09-14 Task 7），默认 None。自创拉普拉斯
                弥散已删除（Z&F22 p.11 "we have no diffusive terms"），弥散由
                q₀ + I₃ 分层通量闭合承载（(4.25)/(4.26)/(4.28)）。保留形参仅为
                既有 runner 兼容；非 None 传值触发 DeprecationWarning 且不再生效。
            dispersion_azimuthal: ⚠️ 已弃用，语义同 dispersion_axial。
            dispersion_dt_ref: ⚠️ 已弃用，语义同 dispersion_axial。
            dispersion_dt_scale: ⚠️ 已弃用，语义同 dispersion_axial。
            enable_stream_yield_gate: B-2 opt-in 开关，默认 False。
                True 且 enable_stream_function=True 时，把 `_yield_gate_wall` 的
                连续冻结度 wall 传进 `solve_stream_function`
                （I₁_eff = I₁·max(1−wall, _WALL_CONDUCTANCE_FLOOR)，地板防
                wall≡1 矩阵奇异；该常量定义于 `stream_function`，勿在此处硬编码）
                ——冻结区流动度→0、Ψ 局部趋于常数 ⇒ 该处轴向速度→0（static wall
                layer，Pelipenko04 (2.6)-(2.8)），总通量守恒、只重新分配到活跃区。
                False（默认）= 不传 wall，逐位等于 HEAD（B-1 的 L1 硬约束）。
                需 enable_yield_gate=True 才会算出非零 wall（否则静默无效，见构造告警）。
            enable_power_law_gap_correction: B-3 opt-in 开关，默认 False。
                True 且 enable_stream_function=True 时，把 `PowerLawGapClosure`
                注入 `solve_stream_function` 的 closure 形参：I₁ → I₁·(H/H̄)^{1/n−1}
                （n = 水泥相 power_law_n，缺省 1.0），使剪切变稀的偏心间隙放大
                进入 2D 速度场。⚠️ 一阶近似，非 B&F25 闭包口径，不得作方法学依据。
                False（默认）= closure 不注入，逐位等于 HEAD（B-1 的 L1 硬约束）。
            enable_stream_function: 速度场路径开关（2026-09-15 Task 9），默认 True。
                True: (w, v) 由 Z&F22 (4.22) 流函数椭圆方程解经 (2.2) 换算得到
                （`_velocity_stream_function`）——浮力（平均密度 ρ·f + 分层
                −Δρ·I₂/(H·I₁)·f）完整经 (4.22) 的 b 向量进入，废止
                ``w·f_amp`` 速度乘子（B1 缺陷，f_amp 是 (4.28) 通量函数不是速度）；
                False: 旧代数流动度路径（`_mobility_profile`/代理 + 截面归一），
                逐位复现 76a91c1 行为（R7 冻结锚护栏，可回退）。
            enable_hb_closure: A-3b（Phase A Task 6）HB 闭包开关，默认 False。
                True 且 ``enable_stream_function=True`` 时，(4.22) 算子经
                ``solve_stream_function_nonlinear`` 非线性外迭代求解（B&F25 HB 两层
                闭包，物理口径参数 + ``velocity_scale = q_half/π``，R-T5-1 REVISED
                推荐路径）；外迭代不收敛（``RuntimeError``）⇒ 显式告警 + **回退牛顿
                线性闭包**继续该时间步（R-T5-3，不静默）。闭包参数：被顶替液（泥浆）
                取其 ``FluidSpec`` 的幂律等价（Bingham/牛顿 → n₁=1/κ₁=PV），
                **τ_Y1 = 泥浆 spec ``yield_stress_pa``**（R-T6-2，2026-09-17 用户
                裁定：泥浆 spec YP 在闭包中真实可用；收敛性由时间步间 warm-start
                等接线级工程保障）；
                顶替液（水泥）n₂/κ₂ 取 ``power_law_n``/``consistency_k``，τ_Y2 由
                ``hb_fix_cement_tau_y`` 门控（H1 ⇒ 0，R-T6-1）。互斥：
                ``enable_power_law_gap_correction`` 被本开关替代（不生效，构造告警）。
                ⚠️ **耗时量级（2026-09-17 实测）**：HB 口径单井 **1.5–3.3 h**
                （hu102 环空容积小、流速高 ⇒ 约 9 h/run），是牛顿口径（2–8.5 min/井）
                的 **45–50×**。瓶颈在闭包求值**不在 CFD**：``solve_g_from_mean_velocity_batch``
                78% / ``HBClosure._evaluate`` 20% / 线性 Poisson 仅 0.8%。出处
                ``docs/superpowers/research/2026-09-17-hb-path-timing-and-reduction-options.md``。
            hb_fix_cement_tau_y: A-3b 水泥 τy 注入开关，默认 False（R-T6-1）。
                True ⇒ 水泥相 τy 注入**两处**（同源同值）：(a) 既有屈服门的混合 τy 场
                （相体积加权处，水泥相贡献由 0 改为常数）；(b) ``enable_hb_closure=True``
                时 ``HBClosure`` 的水泥层 τ_Y2。常数来源优先级：水泥相 ``FluidSpec``
                自带 ``yield_stress_pa``（Bingham/HB，如 ht1_004）以 spec 为准；
                否则取 ``cement_tau_y_by_role`` 映射值。缺相/``MISSING`` 哨兵 ⇒
                **显式告警跳过**（R4，不静默用 0）。
            cement_tau_y_by_role: 水泥相 τy 常数映射 ``{FluidRole 名: Pa}``（键
                "LEAD"/"INTERMEDIATE"/"TAIL"，与 ``cement_yield_stress.yield_stress_by_role``
                的输出同形），默认 None。``hb_fix_cement_tau_y=True`` 时**必须显式提供**
                （R4：τy 只受显式传入值，禁止代码编造；生产口径 Task 7 出数后由用户
                裁定，本仓不内定默认）。两套可用口径的构建方式（Task 0 访问器）：
                HB 拟合 = ``yield_stress_by_role(well_key, "fitted_new")``；
                Bingham-LS 截距 = 逐相读 ``CEMENT_YIELD_STRESS[well][phase]
                .fitted_bingham_ls_intercept_pa``。``hb_fix_cement_tau_y=False`` 时
                提供本映射不消费（构造告警）。
        """
        # ⚠️ 2026-09-14 Task 7 弃用检查：dispersion_* 任一非 None 即弃用警告。
        # （显式传 None 视同默认，不警告——保证未传参的 runner/脚本行为无感。）
        _dispersion_given = (
            dispersion_axial is not None or dispersion_azimuthal is not None
            or dispersion_dt_ref is not None or dispersion_dt_scale is not None
        )
        if _dispersion_given:
            warnings.warn(
                "AnnulusD2DGASolver 的 dispersion_axial/azimuthal/dt_ref/dt_scale 形参已弃用："
                "自创拉普拉斯弥散已删除（Z&F22 p.11 \"we have no diffusive terms\"），"
                "弥散由 q₀ + I₃ 分层通量闭合承载（式 4.25/4.26/4.28）。"
                "显式传值不再生效，请从调用方移除这些参数。",
                DeprecationWarning,
                stacklevel=2,
            )
        # ⚠️ 2026-09-15 Task 10 弃用检查：e_clip 三形参任一偏离 legacy 默认
        # （0.55/0.90/True）即一次性弃用警告。e_clip 硬截断已移除，e = 1−standoff
        # 按 Pelipenko04 (2.1) 文献口径 e∈[0,1) 直取（仅 1e-6 浮点护栏）。
        # 默认/legacy 等值传参完全静默——调用方无感获得文献口径（本 Task 目的）。
        _e_clip_given = (
            e_clip_max != 0.55
            or e_clip_measured_max != 0.90
            or enable_e_clip_ruling is not True
        )
        if _e_clip_given:
            warnings.warn(
                "AnnulusD2DGASolver 的 e_clip_max/e_clip_measured_max/enable_e_clip_ruling "
                "形参已弃用：e_clip 硬截断已移除，e = 1−standoff 按文献口径 e∈[0,1) "
                "（Pelipenko04 (2.1)）直取，不再按数据来源裁剪上限。"
                "显式传值不再生效，请从调用方移除这些参数。",
                DeprecationWarning,
                stacklevel=2,
            )
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
        # ⚠️ 2026-09-15 Task 10：e_clip 截断已移除（Pelipenko04 (2.1) e∈[0,1)），
        # e_clip_max 形参弃用——仅保留属性以兼容旧脚本读值，不再被 _build_geom 消费。
        self.e_clip_max: float = e_clip_max
        self.enable_yield_gate: bool = enable_yield_gate
        self.yield_gate_f_safety: float = yield_gate_f_safety
        # ⚠️ 2026-09-14 Task 7：弥散形参已弃用——仅保留属性以兼容旧脚本读值，
        # 不再被求解过程消费（显式传值已在上方触发 DeprecationWarning）。
        self.dispersion_axial = dispersion_axial
        self.dispersion_azimuthal = dispersion_azimuthal
        self.dispersion_dt_ref = dispersion_dt_ref
        self.dispersion_dt_scale = dispersion_dt_scale
        # 2026-09-06 e_clip 裁定已随 Task 10 截断移除一同退役（形参弃用，仅存属性）
        self.enable_e_clip_ruling = enable_e_clip_ruling
        self.e_clip_measured_max = e_clip_measured_max
        # 幂律缝隙律（构造参数见 docstring）
        self.enable_power_law_gap_law = enable_power_law_gap_law
        # 2026-09-16 B-2：屈服门进流函数算子（opt-in，默认 False ⇒ 逐位=HEAD）。
        # 仅 enable_stream_function=True 路径消费（见 _compute_velocity/_velocity_stream_function）。
        self.enable_stream_yield_gate = enable_stream_yield_gate
        # 2026-09-16 B-3：幂律间隙一阶修正（opt-in，默认 False ⇒ 逐位=HEAD）。
        # 仅 enable_stream_function=True 路径消费（PowerLawGapClosure 注入
        # solve_stream_function 的 closure 形参）。
        self.enable_power_law_gap_correction = enable_power_law_gap_correction
        # 2026-09-15 Task 9：速度场路径开关（True = (4.22) 流函数椭圆方程，
        # False = 旧代数流动度，逐位复现 76a91c1——R7 冻结锚护栏）
        self.enable_stream_function = enable_stream_function
        # B-3 I-1：运行时一次性告警标志（水泥相非幂律时开关静默空转，见
        # _velocity_stream_function）。__init__ 拿不到 fluids，故不能在此判定。
        self._power_law_gap_correction_warned = False

        # ------------------------------------------------------------------ #
        # A-3b（Phase A Task 6，2026-09-17）：HB 闭包接线双开关（opt-in，
        # 默认全关 ⇒ 逐位 = HEAD；判别测试见 tests/contract/test_hb_closure_wiring.py）。
        # ------------------------------------------------------------------ #
        # R4（τy 只受显式传入值，禁止代码编造）：注入常数必须显式给出映射。
        # 键 = FluidRole 名（"LEAD"/"INTERMEDIATE"/"TAIL"，与 Task 0
        # cement_yield_stress.yield_stress_by_role 的输出同形，键统一大写归一）；
        # 值 = Pa（非负有限）或 MISSING 哨兵（字符串 ⇒ 运行期显式告警跳过）。
        self.enable_hb_closure = enable_hb_closure
        self.hb_fix_cement_tau_y = hb_fix_cement_tau_y
        if hb_fix_cement_tau_y and cement_tau_y_by_role is None:
            raise ValueError(
                "hb_fix_cement_tau_y=True 时必须显式提供 cement_tau_y_by_role"
                "（R4：τy 只受显式传入值，禁止代码编造；两套口径 report_given/"
                "fitted_new 由调用方经 cement_yield_stress 访问器显式选择，"
                "本仓不内定默认口径）"
            )
        if cement_tau_y_by_role is not None:
            normalized: Dict[str, float | str] = {}
            for key, val in cement_tau_y_by_role.items():
                key = str(key).upper()
                if isinstance(val, str):
                    # MISSING 哨兵等字符串：留到运行期按 R4 显式告警跳过。
                    normalized[key] = val
                    continue
                v = float(val)
                if not (math.isfinite(v) and v >= 0.0):
                    raise ValueError(
                        f"cement_tau_y_by_role[{key!r}] 须为非负有限值（Pa），得到 {val!r}"
                    )
                normalized[key] = v
            self.cement_tau_y_by_role: Dict[str, float | str] | None = normalized
        else:
            self.cement_tau_y_by_role = None
        # 运行时状态（每 run 重置；见 run() 开头）
        self._active_well_name: str = ""
        self._hb_tau_y_skips: list[str] = []
        self._hb_tau_y_skips_reported = False
        self._hb_wall_on_hb_path_warned = False
        self._hb_closure_memo: dict = {}
        self._hb_prev_G = None  # 上一时间步收敛的 G 场（warm-start；每 run 重置）

        # A-3b 死开关告警（项目 A3 惯例：置真却无可消费路径/输入 ⇒ 一次性告警，
        # 防「死开关被当活杠杆」）。
        _hb_dead: list[str] = []
        if enable_hb_closure and not enable_stream_function:
            _hb_dead.append(
                "enable_hb_closure（仅 enable_stream_function=True 的流函数路径被消费）")
        if cement_tau_y_by_role is not None and not hb_fix_cement_tau_y:
            _hb_dead.append(
                "cement_tau_y_by_role（仅 hb_fix_cement_tau_y=True 时被消费）")
        if hb_fix_cement_tau_y and not enable_yield_gate and not enable_hb_closure:
            _hb_dead.append(
                "hb_fix_cement_tau_y（enable_yield_gate=False 且 enable_hb_closure=False"
                " ⇒ 注入的 τy 场无任何动力学消费点）")
        if _hb_dead:
            warnings.warn(
                f"AnnulusD2DGASolver 的开关 {_hb_dead[0]}"
                + (f"（同批无效：{'; '.join(_hb_dead[1:])}）" if len(_hb_dead) > 1 else "")
                + " 在当前配置下无效。请修正配置或移除这些开关。",
                UserWarning,
                stacklevel=2,
            )
        if enable_hb_closure and enable_power_law_gap_correction:
            warnings.warn(
                "enable_hb_closure=True 与 enable_power_law_gap_correction=True 互斥："
                "HB 闭包（B&F25 (2.13)-(2.15)）替代幂律一阶修正闭包，后者不生效。",
                UserWarning,
                stacklevel=2,
            )

        # ⚠️ 2026-09-16 A3（Task 2 评审）：静默无效开关告警。``enable_stream_yield_gate``
        # 与 ``enable_power_law_gap_correction`` 均只被新路径（``enable_stream_function=True``）
        # 消费；此外屈服门还需 ``enable_yield_gate=True`` 才会算出非零 wall。置真却
        # 无可消费路径时开关静默失效——按项目惯例一次性告警，防「死开关被当活杠杆」。
        _stream_switches = {
            "enable_stream_yield_gate": enable_stream_yield_gate,
            "enable_power_law_gap_correction": enable_power_law_gap_correction,
        }
        _dead = [name for name, on in _stream_switches.items() if on and not enable_stream_function]
        if enable_stream_yield_gate and enable_stream_function and not enable_yield_gate:
            _dead.append("enable_stream_yield_gate(enable_yield_gate=False)")
        if _dead:
            warnings.warn(
                f"AnnulusD2DGASolver 的开关 {', '.join(_dead)} 在当前配置下无效："
                "流函数新路径开关（enable_stream_yield_gate/enable_power_law_gap_correction）"
                "仅在 enable_stream_function=True 时被消费；enable_stream_yield_gate 还需 "
                "enable_yield_gate=True 才会算出非零冻结度 wall。请修正配置或移除这些开关。",
                UserWarning,
                stacklevel=2,
            )

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

        # 2026-09-15 Task 10：e_clip 硬截断移除——e = 1−standoff 按 Pelipenko04
        # (2.1) 文献口径 e∈[0,1) 直取，仅以 1e-6 浮点护栏夹取开区间端点
        # （standoff→0 防 e=1 退化、standoff→1 防零间隙除零）。旧口径
        # clip(1−SO, 0.05, 0.55) 的 0.55 死区与 2026-09-06 实测井 0.90 裁定
        # 一同退役；e_clip 三形参弃用（构造时偏离 legacy 默认即警告）。
        e = np.clip(1.0 - standoff, 1.0e-6, 1.0 - 1.0e-6)
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

    # ------------------------------------------------------------------ #
    # A-3b（Phase A Task 6）：水泥 τy 注入（R-T6-1）+ HB 闭包接线（R-T5-1/3）
    # ------------------------------------------------------------------ #
    def _cement_phase_yield_stress(self, fluid: FluidSpec) -> float:
        """水泥相（lead/tail/intermediate）在**混合 τy 场**中的屈服贡献（R-T6-1(a)）。

        - ``hb_fix_cement_tau_y=False``（默认）⇒ 逐位等于 ``_fluid_yield_stress``
          （HEAD 行为，硬约束）；
        - ``True`` ⇒ 水泥相贡献由 0 改为常数：spec 自带 ``yield_stress_pa``
          （Bingham/HB 水泥）以 spec 为准（HEAD 行为本就非零，"0 改常数"不适用），
          否则取 ``cement_tau_y_by_role`` 映射值；缺相/``MISSING`` 哨兵 ⇒ 显式告警
          跳过（R4，贡献 0 但不静默）。
        """
        if not self.hb_fix_cement_tau_y:
            return self._fluid_yield_stress(fluid)
        return self._hb_resolved_cement_tau_y(fluid)

    def _hb_closure_cement_tau_y(self, fluid: FluidSpec) -> float:
        """HBClosure 水泥层 τ_Y2（R-T6-1(b)，与 (a) 同源同值）。

        ``hb_fix_cement_tau_y=False``（H1，只开闭包）⇒ **恒 0**（裁定字面：即使
        水泥 spec 自带 τy 也不进闭包——H1 的语义就是"只动闭包、不动屈服门"）；
        ``True`` ⇒ 与 `_cement_phase_yield_stress` 同一解析函数（同源同值）。
        """
        if not self.hb_fix_cement_tau_y:
            return 0.0
        return self._hb_resolved_cement_tau_y(fluid)

    def _hb_resolved_cement_tau_y(self, fluid: FluidSpec) -> float:
        """R-T6-1 的水泥 τy 解析（(a)/(b) 两处共用的唯一来源）：spec 优先，
        其次常数映射；缺失 ⇒ 记入显式跳过清单（一次性汇总告警，R4 不静默）。"""
        if fluid.yield_stress_pa is not None:
            self._hb_note_tau_y_skip("spec", (
                f"{fluid.role.name}: 水泥相自带 yield_stress_pa="
                f"{float(fluid.yield_stress_pa):.6g} Pa，常数映射未消费（以 spec 为准）"))
            return float(fluid.yield_stress_pa)
        mapping = self.cement_tau_y_by_role or {}
        value = mapping.get(fluid.role.name.upper())
        if value is None or isinstance(value, str):
            self._hb_note_tau_y_skip("missing", (
                f"{fluid.role.name}: 常数映射缺该水泥相（或为 MISSING 哨兵）——"
                "本相 τy 贡献显式置 0（R4 告警跳过，非静默）"))
            return 0.0
        return float(value)

    def _hb_note_tau_y_skip(self, kind: str, reason: str) -> None:
        """记录一条 τy 跳过留痕（kind: "missing"/"spec"），同一 run 内一次性汇总。"""
        self._hb_tau_y_skips.append((kind, reason))

    def _hb_flush_tau_y_skips(self) -> None:
        """把 τy 跳过留痕一次性汇总告警（每 run 一次；R4：不得静默）。

        ⚠️ 两类 skip 的物理后果不同，汇总尾缀**分流**（评审 Important-1）：
        - ``missing``：该相贡献显式置 0（与 HEAD 同），提示核对映射相覆盖；
        - ``spec``：该相以 spec 自带 yield_stress_pa 为准，**贡献非 0**（与
          HEAD 同值），与映射覆盖无关。
        """
        if self.hb_fix_cement_tau_y and self._hb_tau_y_skips and not self._hb_tau_y_skips_reported:
            unique = "；".join(dict.fromkeys(reason for _k, reason in self._hb_tau_y_skips))
            kinds = {k for k, _r in self._hb_tau_y_skips}
            tail = []
            if "missing" in kinds:
                tail.append("missing 类跳过相的 τy 贡献为 0（与 HEAD 同），"
                            "请核对 cement_tau_y_by_role 的相覆盖")
            if "spec" in kinds:
                tail.append("spec 类跳过相以水泥相自带 yield_stress_pa 为准"
                            "（贡献非 0、与 HEAD 同值），与常数映射覆盖无关")
            warnings.warn(
                f"hb_fix_cement_tau_y=True：井 {self._active_well_name!r} 的水泥 τy 注入"
                f"存在显式跳过项（{unique}）。" + "；".join(tail) + "。",
                UserWarning,
                stacklevel=2,
            )
            self._hb_tau_y_skips_reported = True

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
        """M3 可重启屈服门槛（连续化，2026-09-15 Task 11）：每深度列以该列流动
        最快元（|w| 最大且 w>0）为参考，按平行槽流 τw=G·b/2 外推各元壁面剪应力，
        冻结度 ``wall = clip(1 − τw_extrap/(f·τy), 0, 1)`` 连续取值。

        文献锚点：Pelipenko04 (2.6)-(2.8) 停流区判据 τw < f·τy 的连续近似——
        τw→0（窄缝极限）⇒ wall→1（全冻）；τw ≥ f·τy ⇒ wall=0（可流动）。
        **边界语义**：τw = f·τy 恰 wall=0（可流动）。旧二值判据
        ``τw_extrap ≤ f·τy``（含等号）在阈值处 0/1 跳变（悬崖），两口径仅在
        等号这一零测度集上不同；连续式消除悬崖，窄边冻结带对排量/浓度微扰
        连续响应。

        关键不变量（R3，连续化必须保留）：
        ①参考元（正在流动）本身永不冻结（wall=0）——它在定义上可流动；只有
        壁面剪应力低于 f·τy 的更窄/更慢元才部分冻结；
        ②某列完全无流动（has_flow=False）且水泥已到达，则整列冻结（无法外推
        G）；前锋未到列（cement_ever=0）不冻结。
        τy=0（或 f·τy=0）零除防护语义：无屈服应力 ⇒ 无停流区，wall=0
        （同时消除连续式 0/0 → NaN 的污染路径）。
        停泵期不调用（run() 泵注分支门控）。
        2026-09-07 精简：c_min_residual 形参删除（B2 判据不消费）。
        2026-09-15 Task 11：二值 np.where(immobile,1,0) → 连续冻结度。"""
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
        # 屈服门连续化（2026-09-15 Task 11，Pelipenko04 (2.6)-(2.8) 停流区判据
        # τw < f·τy 的连续近似）：wall = clip(1 − τw_extrap/(f·τy), 0, 1)。
        # τw→0 ⇒ wall→1（全冻）；τw ≥ f·τy ⇒ wall=0（可流动）——阈值处连续，
        # 消除旧二值 np.where(immobile,1,0) 的 0/1 悬崖。
        # τy=0（或 f·τy=0）零除防护：无屈服应力 ⇒ 无冻结，wall=0
        # （掩码外分母置 1 仅作占位，除后按掩码覆写，同时消除 0/0 → NaN）。
        tau_ref = f_safety * tau_y
        tau_ref_safe = np.where(tau_ref > 0.0, tau_ref, 1.0)
        wall_raw = np.where(tau_ref > 0.0, 1.0 - tau_w_extrap / tau_ref_safe, 0.0)
        wall_new = np.clip(wall_raw, 0.0, 1.0)  # 连续冻结度 ∈ [0,1]
        # 前锋未到列不冻结（旧 immobile 掩码的 cement_ever 门，不变量保持）；
        # 残余泥膜由浓度场 c<1 计入 ηE，不再清零速度（2026-09-02 裁定不变）
        wall_new = np.where(cement_ever > 0.0, wall_new, 0.0)
        # 参考元不变量：正在流动的最快元恒 wall=0（R3①）
        wall_new = np.where(ref_mask, 0.0, wall_new)
        # 整列无流动且水泥已到 -> 整列冻结（无法定义参考 G，R3②）
        col_freeze = ~has_flow & np.any(cement_ever > 0.0, axis=0)
        wall_new[:, col_freeze] = 1.0
        return wall_new.astype(float)

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
        # 新增：混合屈服应力（相体积加权）。
        # A-3b（R-T6-1(a)）：水泥相（lead/tail）贡献经 `_cement_phase_yield_stress`——
        # 默认（hb_fix_cement_tau_y=False）逐位等于 `_fluid_yield_stress`（HEAD）；
        # True 时水泥相贡献由 0 改为常数（spec 优先，其次 cement_tau_y_by_role）。
        tau_y = mud * self._fluid_yield_stress(mud_fluid)
        if lead_fluid is not None:
            tau_y += lead * self._cement_phase_yield_stress(lead_fluid)
        if tail_fluid is not None:
            tau_y += tail * self._cement_phase_yield_stress(tail_fluid)
        if spacer_fluid is not None:
            tau_y += spacer * self._fluid_yield_stress(spacer_fluid)
        # A-3b：R4 显式跳过项一次性汇总告警（幂等；默认关时为空操作）。
        self._hb_flush_tau_y_skips()
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
        （见 `_froude_squared_at`，物理量级 O(10⁻³)，八井 [1.2e-3, 3.1e-2]）：
        本方法原先把 ``F2 = 1.0`` 写死在
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
        - ``δ₀ = d̂/r̂ₐ*``：**无量纲**参考间隙比。论文 (2.1) 的窄间隙参数**本身就是**
          ``δ = d̂/r̂ₐ*``——原文写作 ``d̂/(πr̂ₐ*) = δ/π ≪ 1``（π 只在 ``δ/π`` 的写法里，
          不在 ``δ`` 自身上）；(2.6) 的 ``δ₀`` 取其参考值 ⇒ ``δ₀ = d̂/r̂ₐ*``。
          该取值的依据：把 (2.5b)/(2.6) 的
          ``|b| ≈ (ρ−1)/F²``（论文 p.8"b 即浮力向量的大小"）与论文 p.8 的浮力数
          ``b = Δρ·ĝ·d̂²/(μ̂₁ŵ₀)`` 联立，得 ``F²·b = Δρ/ρ̂₁``（Atwood 数），
          即 ``F² = μ̂₁ŵ₀/(ρ̂₁·ĝ·d̂²)``，等价于 ``δ₀·r̂ₐ* = d̂``。
          ⚠️ **不要再给 δ₀ 乘或除 π**（论文的 δ 不含 π）；该取值已用论文 Table 1/2
          反证：复算 b = −50/100/1000/100/1000 逐例吻合，且 ``F²·b`` 与 Atwood 逐位相等。

        传给 `buoyancy.froude_squared` 时 ``gap_scale_m = half_gap_m / mean_radius_m``
        ——二者是**不同的量**（前者无量纲间隙比、后者长度），但乘积恰为 ``d̂``。

        两个调用点（`_compute_velocity` 的 R3 真体力、run 循环的 R2 I3 弥散通量）
        各自按本方法现算，同一时间步内几何/物性口径一致。
        """
        half_gap_m = float(np.mean(geom["H"])) if "H" in geom else 0.5 * float(np.mean(geom["b"]))
        # 退化口径：合成几何的单元测试可能只给 y/phi/b（无 hole/od/H）。
        # 因本式只以乘积 ``δ₀·r̂ₐ* = d̂`` 起作用，此时取 r̂ₐ* = d̂（即 δ₀ = 1）不影响 F²；
        # 环形截面积改取模型自身的 ``∫ b dy``（模型域为**半环空**，与体积 scale 依赖的
        # ``_trapez2d(b) = 0.5·物理环空体积`` 恒等式同源），仍保持 ŵ₀ = q/A 的口径。
        if "hole_mm" in geom and "od_mm" in geom:
            mean_radius_m = float(np.mean((geom["hole_mm"] + geom["od_mm"]) / 4.0)) / 1000.0
            annulus_area_m2 = float(
                np.mean(np.pi / 4.0 * ((geom["hole_mm"] / 1000.0) ** 2
                                       - (geom["od_mm"] / 1000.0) ** 2))
            )
        else:
            mean_radius_m = half_gap_m
            # ⚠️ 半环空面积须 ×2 才是全环空口径（与上面的真实分支同口径）：
            # 模型域是半环空（run 循环 q_half = q/2、体积 scale 用 _trapez2d(b) =
            # 0.5·V_phys），不乘 2 会让 ŵ₀ 大 2 倍 ⇒ F² 大 2 倍 ⇒ f_φ 小 2 倍。
            annulus_area_m2 = 2.0 * float(
                np.trapezoid(np.mean(geom["b"], axis=1), x=geom["y"]))
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

    def _buoyancy_number_at(self, geom: Dict[str, Array], q_m3s: float, w_field: Array,
                            mud_fluid: FluidSpec,
                            lead_fluid: FluidSpec | None,
                            tail_fluid: FluidSpec | None) -> float:
        """按 summary 段同口径现算当前时间步的无量纲浮力数 b（Z&F22 p.8）。

        Task 5 起 ``_compute_velocity`` 逐时间步调用它，把领先阶轴向浮力数接进
        (4.14) 动力学；run() 末尾 summary 块亦复用本方法，保证两处口径**逐位一致**：

        - 顶替液密度 ρ̂₂ = ``buoyancy.displacing_density_kg_m3``（0.67×领浆+0.33×尾浆，
          全仓唯一口径）；
        - 被顶替液黏度 μ̂₁ = ``fluid_apparent_viscosity(泥浆, γ̇ = 6⟨|w|⟩/⟨b⟩)``，
          剪切率约定与 `_compute_props` 一致；缺参数抛错不静默回退；
        - 半间隙 d̂ = ``mean(geom["H"])``（R20 口径，b=2H 逐格自洽）；
        - ŵ₀ = q/A（截面平均速度；环形截面积由 hole/od 计算，退化几何回退 ∫b dy×2）。
        """
        rho_displacing = buoyancy.displacing_density_kg_m3(lead_fluid, tail_fluid, mud_fluid)
        # 泥浆表观黏度的剪切率约定与 _compute_props 一致：γ̇ = 6|w|/b
        shear_rate_mud = (
            6.0 * float(np.mean(np.abs(w_field))) / max(float(np.mean(geom["b"])), 1e-12)
        )
        mu_displaced = buoyancy.fluid_apparent_viscosity(mud_fluid, shear_rate_mud)
        # 几何半间隙 d̂ = (r_o−r_i)/2 = (hole−od)/4，直接取模型自己的半间隙场 geom["H"]
        half_gap_m = float(np.mean(geom["H"])) if "H" in geom else 0.5 * float(np.mean(geom["b"]))
        # 截面平均轴向速度 w₀ = q/A（环形截面积由 hole/od 算）
        if "hole_mm" in geom and "od_mm" in geom:
            annulus_area_m2 = float(
                np.mean(np.pi / 4.0 * ((geom["hole_mm"] / 1000.0) ** 2
                                       - (geom["od_mm"] / 1000.0) ** 2))
            )
        else:
            # 退化口径（合成几何单元测试）：半环空面积 ×2 与真实分支同口径
            annulus_area_m2 = 2.0 * float(
                np.trapezoid(np.mean(geom["b"], axis=1), x=geom["y"]))
        if float(q_m3s) > 0.0 and annulus_area_m2 > 0.0:
            w0_mps = float(q_m3s) / annulus_area_m2
        else:
            # 退让口径：无泵注排量记录（异常输入）时回落到速度场均值
            w0_mps = float(np.mean(np.abs(w_field)))
        return buoyancy.buoyancy_number(
            rho_displacing=rho_displacing,
            rho_displaced=float(mud_fluid.density_kg_m3),
            half_gap_m=half_gap_m,
            mu_displaced=mu_displaced,
            w0_mps=w0_mps,
        )

    def _mobility_base(
        self,
        c_bar: Array,
        geom: Dict[str, Array],
        i1_base: Array | float,
        eta1: Array | float,
        eta2: Array | float,
        n_mix: Array | float,
    ) -> Array:
        """基础流动度 base（无浮力形状）——旧代数路径与 `_mobility_profile` 共享。

        幂律缝隙律（Walton & Bittleston 1991 JFM 222:39-60 槽流）：
        ``base = (b/b̄)^(1+1/n)/η_mix · I₁(c̄,m)``，η_mix 为 (4.23) 两层闭包
        ``1/(c̄³/η₂+(1−c̄³)/η₁)``；牛顿 n=1 退化 w ∝ b²（与 b²/μ 逐位一致）。

        T9 DRY 收敛（controller 第 5 条）：此前该 ~15 行构造在
        `_compute_velocity` 旧代数分支与 `_mobility_profile` 内逐字重复，
        收敛到本方法；两处消费点运算顺序逐字保持（False 路径逐位不变约束）。
        """
        b = geom.get("effective_b", geom["b"])
        b_mean = np.mean(b, axis=0, keepdims=True)
        # (4.23) 两层黏度闭包
        eta_mix = 1.0 / (c_bar ** 3 / np.maximum(eta2, 1.0e-9)
                         + (1.0 - c_bar ** 3) / np.maximum(eta1, 1.0e-9))
        if self.enable_power_law_gap_law:
            # 幂律缝隙律（Walton & Bittleston 1991）：w ∝ b^(1+1/n)，与
            # _compute_velocity 既有口径逐位一致（牛顿 n=1 → b²）
            n_safe = np.clip(n_mix, 0.25, 1.5)
            gap_exponent = 1.0 + 1.0 / n_safe
            base = (b / np.maximum(b_mean, 1.0e-12)) ** gap_exponent / np.maximum(eta_mix, 1.0e-9)
        else:
            # 旧口径：Hele-Shaw b² 流动度
            base = (b / np.maximum(b_mean, 1.0e-12)) ** 2 / np.maximum(eta_mix, 1.0e-9)
        return base * i1_base  # I₁ 乘子（Z&F22 (4.21a)/(4.22) 的 S ∝ 1/(2I₁)）

    def _mobility_profile(
        self,
        c_bar: Array,
        b_num: float,
        geom: Dict[str, Array],
        i1_base: Array | float,
        m_local: float,
        beta_deg: float,
        f2: float,
        eta1: Array | float = 1.0,
        eta2: Array | float = 1.0,
        n_mix: Array | float = 1.0,
        delta_rho: Array | float = 0.0,
    ) -> Array:
        """局部流动度 pref = 基础项 × 浮力形状（Z&F22 (4.14)/(4.22)/(4.23)）。

        ⚠️ **T9（2026-09-15）起仅旧代数路径消费本方法**
        （``enable_stream_function=False``）：新路径速度由 (4.22) 流函数椭圆
        方程直接给出（`_velocity_stream_function`），浮力经 b 向量完整进入，
        本方法的浮力修正是 (4.14) 第二项（通道 B 分层浮力）的代数近似——
        新路径再叠加即双重计入（推导见 `_velocity_stream_function` docstring）。
        保留本方法用于旧路径回退与 K_AXIAL 语义测试。


        **基础项**（幂律缝隙律，Walton & Bittleston 1991 JFM 222:39-60 槽流；
        Pelipenko & Frigaard 2004c）：
        ``base = (b/b̄)^(1+1/n) / η_mix · I₁(c̄,m)``，η_mix 为 (4.23) 两层闭包
        ``1/(c̄³/η₂+(1−c̄³)/η₁)``；牛顿 n=1 退化 w ∝ b²（与 b²/μ 逐位一致）。

        **浮力形状** = 方位项（既有，T1-3b/Task 4）+ 轴向项（Task 5 新增）：

        - 方位项（式 2.5b/4.24）：``clip(Δρ·(I₂/I₁), ±0.5)·f_φ``，
          ``f_φ = r_a·sin(πφ)·sinβ/F²`` 由 `_buoyancy_force_vector` 按 (2.6) 的
          F² 现算（``f2`` 形参透传，R26）。
        - 轴向项（式 (4.14) 第二项 + (4.22) 浮力向量的分层部分）：

          ``(4.14)``: ``(H v̄, H w̄) = −I₁(·, ·) + I₂·Δρ/(H·r_a)·(−f_ξ, f_φ)``，
          轴向分量的浮力通量 = ``I₂·Δρ·f_φ/(H·r_a)``；``(4.13)`` 给
          ``f_φ = r_a·cosβ/F²``，故该项 = ``I₂·Δρ·cosβ/(H·F²)``。
          论文的 Δρ 取 ρ₁−ρ₂（被顶替−顶替），本项目 buoyancy 模块全仓口径为
          Δρ ≡ ρ̂₂−ρ̂₁（顶替−被顶替，b_num 与之同号），两者差一个符号；
          (4.22) 把同一物理写成浮力向量 ``b = (ρ − Δρ_paper·I₂/(H·I₁))·f``，
          其分层部分 ``−Δρ_paper·(I₂/(H·I₁))·f_φ`` 在本项目口径下恰为
          ``+Δρ·(I₂/(H·I₁))·f_φ``。代入 ``F²·b_num = Δρ/ρ̂₁``（Task 4 钉定的
          Atwood 关系，Δρ 无量纲化）与 f_φ，F²、r_a 恰好消去：

            轴向浮力通量 ∝ b_num·cosβ·I₂/(H·I₁)（b>0 重顶替轻 ⇒ 正）。

          以压力驱动项 ``I₁·𝒢`` 为基准的相对修正（𝒢 = 无量纲修正压力梯度）：

            correction(φ) = K_AXIAL·b_num·cosβ·(I₂/I₁)(c̄,m)　—— **H⁰ 形状（R28 裁定）**。

          **K_AXIAL 的量纲推导（R2，provisional）**：pref 与 correction 同为
          无量纲场，量纲只约束到"K_AXIAL 无量纲"；其数值由压力驱动基准给出——
          把 (4.6) 的 ∂p/∂ξ 用 τ̂₀/d̂ = μ̂₁ŵ₀/d̂² 无量纲化后，牛顿参考缝隙流
          （半间隙 d̂、截面平均速度 ŵ₀ 标度）``w̄ = Ĝd̂²/(3μ̂₁)`` 给
          ``𝒢 = Ĝd̂²/(μ̂₁ŵ₀) = 3``，故 ``K_AXIAL = 1/𝒢 = 1/3``（模块常量）。
          局部黏度比使 𝒢 在水泥富集通道自动增大、修正自我衰减；剩余不确定性
          O(1) 由 R27 八井扫描（shape>0）与 Task 12 基准算例裁定，**禁止 clip**。

          **H⁰ 读数（controller 裁定 R28，2026-09-14 fix round 1；documented
          deviation）**：计划 brief 字面 ``I₂·(Δρ/(H·r_a))`` 中的 1/H 判为文本滞后于
          文献分组——首版实现曾按字面乘 ``1/H̃ = d̂/H(φ)``（窄边放大），审查发现三重
          矛盾后裁定删除：

          1. **量纲闭合**：两个闭式 ``d2dga_dispersion_I1/I2(c̄, m)`` 只吃 (c̄, m)、
             不含几何，是 H³/H⁴ 归一化形状；论文 (4.22) 分组 ``I₂/(H·I₁)`` 用全量纲
             闭式（I₁∝H³、I₂∝H⁴）求值时 H⁴/H³ 恰与字面 1/H 相消 ⇒ **H⁰**——修正
             只随 c̄（界面位置）变化，与 H(φ) 无关。
          2. **兄弟项一致**：同一函数、同一 (4.14) 通量结构的方位项（式 2.5b/4.24，
             Task 4 已验收形态）不乘任何 H 因子；轴向项独加 1/H̃ 构成兄弟项间不一致。
          3. **通量归一化**：乘性结构下 base 的 (b/b̄)^(1+1/n) 已把浮力速度恢复到
             论文的 H² 标度；且 ``w = q·pref/Σ(pref·b·dy)`` 使 pref 绝对量级消去、
             只有 φ/z 形状起作用——1/H̃ 是纯形状失真（窄边响应被 d̂/H(φ) 虚假放大，
             并经 H(z) 剖面污染 z 向再分配）。

          方向不受影响（b>0 重顶替轻 ⇒ 正，窄边流动度份额仍被抬高——经 c̄(φ) 与
          通量重分配；只是幅值不再被 1/H̃ 虚假放大）。

        ``b_num=0`` 时轴向项恒为 0，输出与不含浮力项的基础流动度逐位一致（重构锚）。

        Args:
            c_bar: 局部水泥浓度场 c̄ = clip(lead+tail, 0, 1)，(ny, nz)。
            b_num: 无量纲浮力数 b（Z&F22 p.8，`_buoyancy_number_at` 同口径现算）。
            geom: 几何字典（用 ``b``/``effective_b``；H⁰ 读数下轴向项不用 H，R28）。
            i1_base: I₁(c̄,m) 闭式场（Z&F22 (4.21a) 的 H³ 归一化形状）。
            m_local: 黏度比 m = μ_displaced/μ_displacing（标量，R1 auto-m 口径）。
            beta_deg: 井斜角 β（度）；轴向投影 cosβ 来自 (4.13) 的 f_φ = r_a·cosβ/F²。
            f2: Froude 数平方 F²（Z&F22 (2.6)，R26 透传，供方位 f_φ 使用）。
            eta1: 泥浆相黏度场 η₁（(4.23) 闭包输入；默认 1 供纯函数单测）。
            eta2: 水泥相黏度场 η₂（默认 1）。
            n_mix: 混合物幂律指数场（缝隙律 1+1/n；默认 1 = 牛顿）。
            delta_rho: 局部密度差场（g/cc，rho − ρ̂₁/1000；默认 0 关闭方位项）。

        Returns:
            pref 未饱和乘积 base·buoyancy_shape（无 max/wall——调用方按
            M2/壁面分支自行施加 ``np.maximum(·, 1e-8)`` 与 ``(1−wall)``）。
        """
        # T9 DRY 收敛：base 构造收敛到 `_mobility_base`（运算逐字保持，逐位等价）。
        base = self._mobility_base(c_bar, geom, i1_base, eta1, eta2, n_mix)

        # --- 浮力形状 ---
        beta_rad = np.deg2rad(float(beta_deg))
        cos_beta = float(np.cos(beta_rad))
        i_ratio = d2dga_dispersion_I2(c_bar, m_local) / np.maximum(i1_base, 1.0e-12)
        # 方位项（T1-3b/Task 4 既有路径，式 2.5b/4.24）：f_φ 带 1/F² 标定（R26 通路）
        f_phi_arr, _ = self._buoyancy_force_vector(geom, beta_deg, f2)
        az_correction = np.clip(delta_rho * i_ratio, -0.5, 0.5)
        # 轴向项（Task 5 新增，式 (4.14) 第二项/(4.22) 分层部分）：
        # correction = K_AXIAL·b_num·cosβ·(I₂/I₁)——H⁰ 形状（R28 裁定：闭式 H⁴/H³
        # 与字面 1/H 相消；首版的 1/H̃ = d̂/H 因子已删，见 docstring 三重证据链）。
        axial_correction = K_AXIAL * float(b_num) * cos_beta * i_ratio
        buoyancy_shape = 1.0 + az_correction * f_phi_arr + axial_correction
        return base * buoyancy_shape

    def _velocity_stream_function(
        self,
        lead: Array,
        tail: Array,
        geom: Dict[str, Array],
        q_m3s: float,
        w_prev: Array,
        mud_fluid: FluidSpec,
        lead_fluid: FluidSpec | None,
        tail_fluid: FluidSpec | None,
        wall: Array | None = None,
    ) -> Tuple[Array, Array]:
        """新路径速度场：Z&F22 (4.22) 流函数椭圆方程 + (2.2) 换算（2026-09-15 Task 9）。

        **文献结构（推导全文，controller 裁定第 2 条）**
        ------------------------------------------------
        Z&F22 的速度场由且仅由两个关系确定：(4.22) 椭圆方程
        ``0 = ∇a·[S + b]``、``S = (r_a/(2I₁))∇aΨ``（(4.23)）与 (2.2) 速度换算
        ``w̄ = ∂φΨ/(2r_aH)``、``v̄ = −∂ξΨ/(2H)``。密度**只**通过浮力向量进入流动：

            b = (ρ − Δρ·I₂/(H·I₁))·f，  f = (r_a·cosβ/F², r_a·sin(πφ)·sinβ/F²)
                （(4.22) + (4.13)；(2.5b) 概览式同分组：φ-槽配 cosβ、ξ-槽配 sinπφ·sinβ）

        论文 p.A32-20 明言其两部分（(4.24) 段）："the first part of ∇a·b
        represents changes in the mean density in the direction f. These
        buoyancy gradients are exactly those considered in the models of
        Bittleston et al. (2002) and Maleki & Frigaard (2017). The second
        part results specifically from the layered flow."——即

        - **通道 A（平均密度）**：``ρ·f``，ρ = (1−c̄)ρ₁ + c̄ρ₂（无量纲，ρ₁=1）。
          它就是 (4.14) 压力槽里的 "modified with the mean static pressure
          gradient" 项 ``∂ξp + ρf_φ/r_a``——领先阶轴向重力（cosβ 投影）经此进入；
          源项只含其 φ/ξ-梯度（均匀密度不驱动流动，物理正确）。
        - **通道 B（分层流动）**：``−Δρ·I₂/(H·I₁)·f``，Δρ = ρ₁−ρ₂（论文符号）。
          它就是 (4.14) 第二项 ``I₂·(Δρ/(H·r_a))·(−f_ξ, f_φ)``。

        **浮力双通道整合结论（本方法的核心设计裁定）**：旧代数路径的
        ``_mobility_profile`` 轴向项 ``K_AXIAL·b_num·cosβ·(I₂/I₁)`` 是 (4.14)
        第二项/**通道 B** 的代数近似（见其 docstring），通道 A 在旧路径中缺席。
        新路径里 (4.22) 的 b 向量已**完整承载通道 A+B**，故 `_mobility_profile`
        的浮力修正**不得再叠加**（通道 B 双重计入），新路径速度 = 纯
        (4.22)+(2.2) 解，`_mobility_profile`/截面归一完全不消费。
        排查论证：①若新路径再叠加 K_AXIAL 项，通道 B 以两种离散形式进入两次
        （mobility 乘子 + 椭圆源项），而论文只有一次；②若以为 b 向量只承载
        O(ε) 方位浮力、把领先阶轴向浮力留给 mobility，则与 (4.22) 推导矛盾——
        压力在交叉微分中已被消去，b 是密度的唯一入口（(2.3) 是唯一动量关系），
        轴向重力经 ρ·f（通道 A）的 φ-梯度进入竖直井机制（A32-22 页
        "the elliptic equation ... is driven by gradients in b"，须 χ 的
        φ-梯度即偏心前缘，均匀 c̄ 下源项为零是论文自身的机制属性）；
        ③b_num 的幅值信息不经独立通路进入：1/F² = b_num·ρ̂₁/Δρ̂（Task 4 钉定的
        Atwood 关系），f = r_a/F² 已内含 b_num 的量级——不再需要、也不再消费
        K_AXIAL/b_num（K_AXIAL 仅属旧代数路径）。

        **λ_op 修正（量纲→无量纲装配匹配，T9 推导 D3）**
        ------------------------------------------------
        (4.22) 在论文中是无量纲方程（§2/§4.1：径向标度 d̂、方位与轴向标度
        πr̂ₐ*、黏度 μ̂₁、密度 ρ̂₁、应力 τ̂₀=μ̂₁ŵ₀/d̂，(4.5)-(4.8) 的重力项
        ρ/F² 与压力梯度同槽平衡）。本仓 ``solve_stream_function`` 按量纲 I₁
        （m³/(Pa·s)）、米制 r_a/ξ 装配同一形状的方程。将量纲装配算子与论文
        无量纲算子逐槽匹配（φ-槽与 ξ-槽同时成立）：

            op_模块 = λ_op·op_论文，  λ_op = μ̂₁/(d̂³·π·r̂ₐ*)，

        其中 r_a^nd = r̂ₐ*/(πr̂ₐ*) = 1/π 由双槽各向异性匹配钉定（与论文
        "azimuthal and axial lengths have been scaled with πr̂ₐ*" 一致）；
        源项两侧逐位相等（φ-槽的 r_a 在 (1/r_a)∂φ(χ·r_a/F²) 中消去）。
        故源项必须乘 λ_op 才保持论文的源/算子比——缺此因子浮力响应被压低
        λ_op⁻¹·… 倍（case5 约 5×10⁶）。交叉验证：K_AXIAL 的 𝒢=Ĝd̂²/(μ̂₁ŵ₀)
        槽流压力标度、(4.29) 的 U=3mη₂F²/(ρH³r_a) 一致性均不矛盾。

        **Ψ 的物理缩放（度量换算推导，T9 推导 D2；修正 brief 字面 "Q/2"）**
        ------------------------------------------------
        模块单位 BC 不变量 ``∫₀¹ 2r_aH·w̄ dφ = 1`` 是 **φ-度量**通量；模型输运/
        体积层用的是**弧长度量**（geom["y"] = 半周长弧、geom["b"] 经 _build_geom
        体积标定），同一速度场的弧长列通量 ``∫w·b·dy = π·∫2r_aH·w̄ dφ``
        （弧元 dy = π·r_a·dφ）——数值验证逐位列比恰为 π（task9_probe/
        verify_pi_scaling.py，2026-09-15）。故物理半环空通量 = Q/2 要求

            w = w̄_unit·(Q/2)/π = w̄_unit·q_half/π。

        同理模块 (w,v) 对满足 ``π·∂y(bv) + ∂s(bw) = 0`` 而非物理连续性
        ``∂y(bv) + ∂s(bw) = 0``（同一 ξ 量纲混用：论文 ξ 以 πr̂ₐ* 标定、模块
        ξ 用米）——v 需再乘 π 方与 w 构成物理连续对，即

            v = v̄_unit·q_half。

        两缩放联立后数值验证：列通量逐位 = q_half、物理连续性残差 ~8e-15
        （同上探针脚本 [4][5] 段）。stream_function 模块 docstring 的
        "整体乘以 Q_half" 表述已按本推导修正（T9，注释级改动）。

        **旧路径差异声明（新路径不承载的旧机制）**：①屈服门 wall/`f_safety`
        （两层牛顿闭包无屈服项，壁面带慢速由闭包自身体现——documented
        deviation；B-2 起可经 ``enable_stream_yield_gate=True`` 把 wall 送进
        本路径的算子，默认关）；②M2 流态修正（默认关）；③幂律缝隙律 (b/b̄)^(1+1/n)
        （新路径默认牛顿两层闭包 (4.21)，流变经标量表观黏度 η₁/η₂/m 进入，
        m 无 clip——论文无 clip；B-3 起可经 ``enable_power_law_gap_correction=True``
        叠加一阶间隙修正 I₁·(H/H̄)^{1/n−1}，默认关）；④f_amp 速度乘子（B1 缺陷，
        见 run()）。
        均可经 ``enable_stream_function=False`` 回退到旧路径。

        Args:
            lead/tail: 相浓度场 (ny,nz)。
            geom: 几何字典（`_build_geom` 产物，含 phi/H/s/y/hole_mm/od_mm）。
            q_m3s: 全环空排量（m³/s）；≤0（泵停）直接返回零场——与旧路径
                q_half=0 的零场逐位等价，且省一次椭圆解。
            w_prev: 上一步轴向速度场（代表性剪切率与 F² 口径的输入）。
            mud_fluid: 被顶替液（ρ̂₁/μ̂₁ 口径流体）。
            lead_fluid/tail_fluid: 顶替液（水泥）——η₂ 取 lead（缺则 tail）的
                表观黏度；ρ̂₂ 取全仓唯一口径 ``displacing_density_kg_m3``
                （0.67×领浆+0.33×尾浆，与 b_num/summary 一致）。
            wall: (ny,nz) 屈服门冻结度 ∈ [0,1]（B-2 opt-in，Pelipenko04
                (2.6)-(2.8)），``None`` ⇒ 不进算子（逐位=HEAD）。非 None 时经
                ``solve_stream_function(wall=...)`` 把冻结区流动度压向 0 ⇒ 该处
                速度→0、通量重分配到活跃区。调用点由
                ``enable_stream_yield_gate`` 门控（默认 False）。

        Returns:
            (w, v)：物理量纲速度场 (ny,nz)，轴向/方位间隙平均速度（m/s）。
        """
        ny, nz = self.ny, self.nz
        if not (float(q_m3s) > 0.0):
            # 泵停：旧代数路径在 q_half=0 时 w=v=0（精确零场），逐位等价。
            return (np.zeros((ny, nz)), np.zeros((ny, nz)))
        if "H" not in geom:
            raise ValueError(
                "enable_stream_function=True 需要 geom['H']（Z&F22 半隙场，"
                "(4.22)/(2.2) 的核心几何）；生产 _build_geom 恒提供。合成测试几何"
                "请按 H = b/2 补齐，或改用 enable_stream_function=False。"
            )

        c_bar = np.clip(lead + tail, 0.0, 1.0)

        # ---- 标量黏度口径（与 F²/b_num 同一 Representative 剪切率）------------
        # γ̇ = 6⟨|w|⟩/⟨b⟩：_froude_squared_at/_buoyancy_number_at/_compute_props 同约定。
        shear_rate = (6.0 * float(np.mean(np.abs(w_prev)))
                      / max(float(np.mean(geom["b"])), 1e-12))
        eta1 = buoyancy.fluid_apparent_viscosity(mud_fluid, shear_rate)
        cement_fluid = lead_fluid if lead_fluid is not None else tail_fluid
        if cement_fluid is not None:
            eta2 = buoyancy.fluid_apparent_viscosity(cement_fluid, shear_rate)
        else:
            eta2 = eta1  # 无水泥相：两层退化为均一（m=1，χ 的 φ-梯度仍在）
        m_ratio = eta1 / eta2

        # ---- F²（(2.6)，与旧路径同口径现算）-----------------------------------
        f2 = self._froude_squared_at(geom, q_m3s, w_prev, mud_fluid)
        beta_deg_local = float(np.mean(geom.get("inc_deg", np.zeros(nz))))
        # _buoyancy_force_vector 返回序 = 代码口径 (方位形状, 轴向形状)；
        # (4.13) 的 f_φ=轴向(r_a·cosβ/F²) 进 b 的 φ-槽、f_ξ=方位(sinπφ·sinβ/F²)
        # 进 ξ-槽（(2.5b)/(4.22) 字面分组，R29）。命名陷阱：代码 f_phi=方位，
        # 论文 f_φ=轴向——此处按形状显式重命名防错位。
        f_azimuthal, f_axial = self._buoyancy_force_vector(geom, beta_deg_local, f2)

        # ---- χ = ρ − Δρ·I₂/(H·I₁)（(4.22) 浮力向量的标量因子，(ny,nz) 场）------
        # ρ₁ = 1（被顶替液密度标定）、ρ₂ = ρ̂₂/ρ̂₁（无量纲顶替液密度）；
        # Δρ = ρ₁ − ρ₂（论文符号，重顶替轻时为负）；I₁/I₂ 用全量纲闭式 (4.21a,b)
        # （I₂∝H⁴、I₁·H∝H⁴ ⇒ I₂/(H·I₁) 为 H⁰ 无量纲场，R28 同款消去）。
        rho1_kg_m3 = float(mud_fluid.density_kg_m3)
        rho2_kg_m3 = float(buoyancy.displacing_density_kg_m3(lead_fluid, tail_fluid, mud_fluid))
        rho2_nd = rho2_kg_m3 / rho1_kg_m3
        rho_nd = 1.0 + c_bar * (rho2_nd - 1.0)                     # (1−c̄)·ρ₁ + c̄·ρ₂
        i1_field = np.asarray(mobility_i1(c_bar, m_ratio, eta1=eta1, eta2=eta2,
                                          H=geom["H"]), dtype=float)
        i2_field = np.asarray(mobility_i2(c_bar, m_ratio, eta1=eta1, eta2=eta2,
                                          H=geom["H"]), dtype=float)
        chi = rho_nd + (rho2_nd - 1.0) * i2_field / np.maximum(geom["H"] * i1_field, 1.0e-30)

        # ---- b 向量装配（(2,ny,nz)，(4.22) 字面分组唯一口径，R29）--------------
        b_field = np.empty((2, ny, nz), dtype=float)
        b_field[0] = chi * f_axial       # φ-槽：(4.13) f_φ = r_a·cosβ/F²（轴向重力）
        b_field[1] = chi * f_azimuthal   # ξ-槽：(4.13) f_ξ = r_a·sinπφ·sinβ/F²（方位重力）

        # ---- λ_op 修正（T9 推导 D3，2026-09-15）--------------------------------
        # 论文 (4.22) 是**无量纲**方程（径向标度 d̂、方位/轴向标度 πr̂ₐ*、黏度 μ̂₁、
        # 应力 τ̂₀=μ̂₁ŵ₀/d̂，Z&F22 §2/§4.1 (4.4)-(4.8)）；本仓 solve_stream_function
        # 按量纲 I₁（m³/(Pa·s)）与米制 r_a/ξ 装配。把量纲装配式与论文无量纲式逐槽
        # 匹配（φ-槽与 ξ-槽**同时**成立，r_a^nd = r̂ₐ*/(πr̂ₐ*) = 1/π 由双槽各向异性
        # 匹配钉定），得 op_模块 = λ_op·op_论文、src_模块 = src_论文，故源项须乘
        #     λ_op = μ̂₁/(d̂³·π·r̂ₐ*)
        # 才等价于论文无量纲源/算子比。缺此因子浮力响应被压低 ~μ̂₁/(d̂³πr̂ₐ*) 倍
        # （基准算例 ~5×10⁶——2026-09-15 首跑实测 case1/case2 t_br 差仅 0.004，
        # 论文 0.51，即此根因）。推导全文见 task-9-report 与本方法 docstring。
        half_gap_m = float(np.mean(geom["H"]))            # d̂（R20 口径）
        r_a_m = float(np.mean((geom["hole_mm"] + geom["od_mm"]) / 4.0)) / 1000.0
        lambda_op = eta1 / (half_gap_m ** 3 * np.pi * r_a_m)
        b_field = b_field * lambda_op

        # ---- 椭圆解 + (2.2) 换算 + 物理缩放（推导见 docstring）-----------------
        # B-3（opt-in，默认关）：代表幂律指数取自水泥相（lead 优先，缺则 tail）；
        # 非幂律/HB 或指数未给时 n_rep=1.0 ⇒ PowerLawGapClosure 因子恒 1、
        # mobility 逐位退化为 NewtonianClosure（默认关时 closure=None，
        # solve_stream_function 内部同样落到 NewtonianClosure ⇒ 逐位=HEAD）。
        # ⚠️ M-2 已知不一致（一阶近似包络内）：上方 χ 的 i2_field/i1_field 与
        # solve_stream_function 的 b 场装配用的是**未修正**牛顿闭式 I₁/I₂
        # （I₂/(H·I₁) 为 H⁰ 的 H 幂消去）；B-3 开启时椭圆算子用**已修正** I₁，
        # 两者口径不同。属一阶近似包络内的已知不一致，不单独修正。
        if self.enable_power_law_gap_correction:
            n_rep = 1.0
            n_cement = getattr(cement_fluid, "power_law_n", None) if cement_fluid is not None else None
            if n_cement:
                n_rep = float(n_cement)
            elif not self._power_law_gap_correction_warned:
                # I-1：开关置真但水泥相非幂律（power_law_n is None，如 Bingham/牛顿）
                # ⇒ n_rep 回落 1.0、因子恒 1、**静默空转**——一次性告警，防下游
                # （如 Task 4 变体矩阵）把 Δη=0 误读成「物理中性」。__init__ 拿不到
                # fluids，故判据只能落在运行时（本处），这是唯一能精确定位流体之处。
                self._power_law_gap_correction_warned = True
                warnings.warn(
                    "enable_power_law_gap_correction=True 但水泥相非幂律"
                    f"（power_law_n is None，流体={getattr(cement_fluid, 'name', None)!r}）："
                    "幂律间隙修正静默空转（n_rep 回落 1.0、因子恒 1）。"
                    "若本井 LEAD/TAIL 为 Bingham/牛顿，该开关不产生任何效应——"
                    "勿把由此得到的零差异读作「物理中性」。",
                    UserWarning,
                    stacklevel=2,
                )
            from cemdisp.models2d.hb_closure import PowerLawGapClosure
            closure = PowerLawGapClosure(n_rep)
        else:
            closure = None
        # === A-3b（Phase A Task 6）：HB 闭包非线性外迭代（opt-in，默认关）======
        # enable_hb_closure=True 且有水泥相 ⇒ 走 solve_stream_function_nonlinear
        # （B&F25 HB 两层闭包；物理口径参数 + velocity_scale = ŵ = q_half/π，
        # R-T5-1 REVISED 推荐路径；Gb 恒 (0,0)，浮力全部经 b_field——不得两处
        # 同时算同一份浮力）。不收敛 ⇒ 显式告警 + 回退牛顿线性闭包（R-T5-3）。
        # 默认关 ⇒ 下方线性调用逐位 = HEAD。
        if self.enable_hb_closure and cement_fluid is not None:
            psi = self._solve_stream_function_hb(
                geom, c_bar, b_field, float(q_m3s),
                mud_fluid, cement_fluid, eta1, eta2, m_ratio, shear_rate, wall, ny, nz)
        else:
            psi = solve_stream_function(geom, c_bar, eta1, eta2, m_ratio, b_field,
                                        closure=closure, wall=wall, ny=ny, nz=nz)
        w_unit, v_unit = velocity_from_stream_function(psi, geom)
        q_half = float(q_m3s) / 2.0
        w = w_unit * (q_half / np.pi)
        v = v_unit * q_half
        return w, v

    def _hb_closure_for(self, mud_fluid: FluidSpec, cement_fluid: FluidSpec,
                        gamma0: float) -> HBClosure:
        """构造/复用 HB 闭包（**物理口径**参数；B&F25 (2.13)-(2.15) 求值口径）。

        - 被顶替液（泥浆，流体 1/壁面带）：n₁/κ₁ 取 `_phase_power_law_params`
          （Bingham/牛顿 → (1, PV)、POWER_LAW → (n, K)）；**τ_Y1 = 泥浆 spec
          ``yield_stress_pa``**（无则 0）——用户裁定 R-T6-2（2026-09-17，"想想
          办法"）：泥浆 spec YP 在闭包中真实可用，收敛性由接线级工程保障
          （时间步间 warm-start，见 `_solve_stream_function_hb`；探针与分级
          证据见 task-6-report.md 修复轮 1）；
        - 顶替液（水泥，流体 2/中线带）：n₂/κ₂ 取 `_phase_power_law_params`
          （POWER_LAW → (power_law_n, consistency_k)）、τ_Y2 经
          `_hb_closure_cement_tau_y`（R-T6-1(b)：H1 ⇒ 0、H3 ⇒ 与屈服门同源同值）；
        - ``m``/``B`` 由 ``hb_groups`` 在代表性剪切率（与 F²/b_num 同口径）下算出，
          **仅供报告/索引与 R-T1-6 地板参照**，不参与运行时求值（构造后冻结）。
        """
        n1, k1 = self._phase_power_law_params(mud_fluid)
        n2, k2 = self._phase_power_law_params(cement_fluid)
        # τ_Y1 = 泥浆 spec YP（R-T6-2，2026-09-17 用户裁定：泥浆 spec YP 在闭包中
        # 真实可用；收敛性由接线级 warm-start 工程保障，见 _solve_stream_function_hb）。
        tau_y1 = self._fluid_yield_stress(mud_fluid)
        tau_y2 = self._hb_closure_cement_tau_y(cement_fluid)
        self._hb_flush_tau_y_skips()
        key = (n1, n2, k1, k2, tau_y1, tau_y2)
        cached = self._hb_closure_memo.get(key)
        if cached is not None:
            return cached
        # m/B 是报告口径量：gamma0 取代表性剪切率（可能为 0 ⇒ 1e-9 护栏，
        # 仅防 hb_groups 构造校验拒绝；不影响运行时求值）。
        groups = hb_groups(k1, k2, n1, n2, tau_y1, tau_y2, max(float(gamma0), 1e-9))
        closure = HBClosure(n=(n1, n2), kappa=(k1, k2), tau_y=(tau_y1, tau_y2),
                            m=float(groups["m"]), B=float(groups["B"]), ny_gap=201)
        self._hb_closure_memo[key] = closure
        return closure

    def _solve_stream_function_hb(self, geom: Dict[str, Array], c_bar: Array,
                                  b_field: Array, q_m3s: float,
                                  mud_fluid: FluidSpec, cement_fluid: FluidSpec,
                                  eta1: float, eta2: float, m_ratio: float,
                                  shear_rate: float, wall: Array | None,
                                  ny: int, nz: int) -> Array:
        """A-3b：HB 闭包的 (4.22) 非线性外迭代 + R-T5-3 回退（本时间步不静默）。

        - ``velocity_scale = ŵ = q_half/π``（生产换算锚 ``w = w_unit·(q_half/π)``，
          R-T5-1 REVISED 推荐路径：物理速度喂反求 + 物理 (κ, τ_Y) 参数）；
        - ``Gb=(0,0)``：Phase A 浮力全部经 (4.22) 的 ``b_field``（Task 4 警告：
          不得两处同时算同一份浮力）；
        - 外迭代 ``max_outer=400``（接线配置）：悬崖邻域 ū↔G↔I₁ 耦合的 ‖ΔĪ₁‖
          收敛率实测 ~0.85-0.95/轮（ω=0.5，随场规模/悬崖格数变慢），Task 5
          默认 50 轮不够 ⇒ 提到 400 轮（小场冷启动 ~73 轮、快变步 ~102-155 轮、
          呼101 量级场 40×250 实测 ~280 轮）。
          ω 保持 Task 5 默认 0.5；自适应欠松弛（R-T6-2 阶梯 (ii)，机制已实现于
          ``solve_stream_function_nonlinear``）生产接线**不启用**——实测对当前
          失败模式无改善（修复轮 1 证据）；
        - 初始 G（启发式，可能低估）：warm-start 优先（上一步收敛 G 场）；冷启动
          用**牛顿当量两段式**（先解牛顿线性 Ψ、按 ``G = H·ū/I₁_牛顿`` 逐格反推
          注入）——因 ``I₁_HB ≤ I₁_牛顿`` 它系统性低估驱动，仅作外迭代初值；
          稳健求根（反求工具的割线 + 夹逼二分回退，R-T6-2 REVISED）下外迭代
          会自行修正，不依赖其准确；
        - static 分支判据（R-T6-2 REVISED 第 2 条）：仅 |ū| ≤ 1e-12×max|ū| 的格
          取冻结/地板，其余格一律找流动根（修复轮 1 的"救援格滞回"已废——
          它与整批 rtol 阶梯组合曾把本应流动的格扫进冻结分支，实测退化解：
          τ_Y2 扫描 max|w| 逐位相同、生产场 288/288 全 undefined）；
        - 外迭代 ``RuntimeError`` ⇒ **显式 RuntimeWarning + 回退牛顿线性闭包**
          （``closure=None``，即既有线性路径；R-T5-3：不静默，回退由调用方负责）。
        """
        if wall is not None and not self._hb_wall_on_hb_path_warned:
            # B-2（enable_stream_yield_gate）与 HB 同开：非线性入口无 wall 形参，
            # 冻结度不进 HB 算子（已知边界，显式告警不静默）；回退线性路径仍消费 wall。
            self._hb_wall_on_hb_path_warned = True
            warnings.warn(
                "enable_hb_closure=True 路径当前不消费屈服门冻结度 wall"
                "（solve_stream_function_nonlinear 无 wall 入口，Phase A 边界）；"
                "本时间步的冻结度不进 HB 算子，仅回退线性路径消费。",
                UserWarning,
                stacklevel=2,
            )
        closure = self._hb_closure_for(mud_fluid, cement_fluid, shear_rate)
        q_half = float(q_m3s) / 2.0
        velocity_scale = q_half / np.pi
        # R-T6-2 阶梯 (i)：时间步间 warm-start——上一时间步**收敛**的 G 场作首轮
        # 闭包状态（瞬态场连续演化 ⇒ 加速器）；首步/回退步无缓存 ⇒ 冷启动
        # （牛顿当量两段式启发式，见方法 docstring：可能低估 G，稳健求根下
        # 外迭代自行修正）。跨步只传递初值、不冻结任何格的归属（每步重新反求）。
        prev_g = self._hb_prev_G
        initial_g = None
        if prev_g is not None:
            initial_g = prev_g
        else:
            # 冷启动（R-T6-2 阶梯 (i) 的冷启动半段）：**牛顿当量两段式**——先解
            # 一次牛顿线性 Ψ，按 G = H·ū/I₁_牛顿 逐格反推注入（启发式：因
            # I₁_HB ≤ I₁_牛顿 系统性低估驱动，外迭代在稳健求根下自行修正）。
            psi0 = solve_stream_function(geom, c_bar, eta1, eta2, m_ratio, b_field,
                                         closure=None, ny=ny, nz=nz)
            w0, v0 = velocity_from_stream_function(psi0, geom)
            u0 = np.hypot(w0, _VELOCITY_COMPONENT_RATIO * v0) * velocity_scale
            h_arr = np.asarray(geom["H"], dtype=float)
            i1_newton = np.asarray(
                mobility_i1(c_bar, m_ratio, eta1=eta1, eta2=eta2, H=h_arr), dtype=float)
            closure.set_pressure_gradient(h_arr * u0 / np.maximum(i1_newton, 1e-30))
        try:
            psi = solve_stream_function_nonlinear(
                geom, c_bar, closure, b_field, Gb=(0.0, 0.0),
                velocity_scale=velocity_scale, initial_G=initial_g,
                max_outer=400)
        except RuntimeError as exc:
            # 回退步的 G 是迭代中途态，不缓存（下一步仍冷启动，R-T5-3 最后防线）。
            self._hb_prev_G = None
            warnings.warn(
                f"HB 闭包非线性外迭代未收敛，本时间步回退牛顿线性闭包路径"
                f"（R-T5-3，不静默）：{exc}",
                RuntimeWarning,
                stacklevel=2,
            )
            return solve_stream_function(geom, c_bar, eta1, eta2, m_ratio, b_field,
                                         closure=None, wall=wall, ny=ny, nz=nz)
        # 收敛 ⇒ 缓存收敛 G 场，供下一时间步 warm-start（阶梯 (i)）。
        self._hb_prev_G = closure.current_G
        return psi

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

        双路径（Task 9，2026-09-15）：

        - ``enable_stream_function=True``（默认）：速度 (w, v) 由 Z&F22 (4.22)
          流函数椭圆方程解经 (2.2) 换算（`_velocity_stream_function`，含浮力
          双通道推导与 Ψ 物理缩放的度量换算推导）；
        - ``enable_stream_function=False``：旧代数流动度路径——
          1. 计算局部混合流体的表观粘度与密度；
          2. 以幂律缝隙律 ``(b/b̄)^(1+1/n)/η_mix·I₁`` 构造偏心通道主导局部流动度
             （Task 5 起封装在 `_mobility_profile`/`_mobility_base`）；
          3. 根据密度差（顶替液 vs 被顶替液）与浮力数 b 计算浮力修正项
             （方位式 2.5b/4.24 + 轴向式 4.14/4.22）；
          4. 用浮力修正项调整宽边/窄边速度分配；
          5. 由截面排量约束得到轴向速度 ``w``。

        2026-09-07 精简：flusher 形参删除（相降级，见 _compute_props 注）。
        Task 5（2026-09-14）：流动度构造抽为 `_mobility_profile` 纯函数。
        Task 9（2026-09-15）：新路径接入流函数；base 构造 DRY 收敛到
        `_mobility_base`；返回签名（12 元组）不变。

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
        # mu_reg 保留用于 Re 诊断
        c_bar = np.clip(lead + tail, 0.0, 1.0)  # 局部水泥浓度

        # === T9 新路径（enable_stream_function=True，默认）：(4.22) 椭圆方程 ===
        # 浮力（平均密度 ρ·f + 分层 −Δρ·I₂/(H·I₁)·f）完整经 (4.22) 的 b 向量进入
        # （推导与双重计入排查见 `_velocity_stream_function` docstring）；
        # _mobility_profile/K_AXIAL/f_amp 速度乘子均不消费（旧代数近似废止）。
        if self.enable_stream_function:
            w, v = self._velocity_stream_function(
                lead, tail, geom, q_m3s, w_prev, mud_fluid, lead_fluid, tail_fluid,
                wall=(wall if self.enable_stream_yield_gate else None),
            )
            return w, v, mu_reg, rho, mud, Re, mu_turbulent, m_field, tau_y, eta2, n_mix, kappa_mix

        # === 旧代数路径（enable_stream_function=False；逐位复现 76a91c1，可回退）===
        # T1-3b: I₁(c̄,m) 乘子（Zhang 2022 式 4.22，S ∝ 1/(2I₁)）
        # 牛顿极限 m→1 时 I₁=1/3，不改变 base 形状；m≠1 时修正方位分布
        # 2026-09-07 R0 分支删除：m 恒由 _compute_props 自动计算（enable_d2dga_auto_m
        # 恒 True，标量 d2dga_viscosity_ratio 路径为旧论文 R0 状态已移除）。
        m_local = float(np.mean(m_field))
        i1_base = d2dga_dispersion_I1(c_bar, m_local)

        # === 速度场流动度：偏心通道主导 + 浮力修正 ===
        # density_contrast > 0 表示顶替液更重（水泥重 vs 泥浆轻），有助于窄边推进
        # density_contrast < 0 表示顶替液更轻，加剧宽边窜流
        # Task 3: 顶替液密度统一走 buoyancy.displacing_density_kg_m3（全仓唯一口径，
        # 0.67×领浆 + 0.33×尾浆），消除与 summary 段浮力数口径不一致导致的 b 符号翻转。
        rho_disp = buoyancy.displacing_density_kg_m3(lead_fluid, tail_fluid, mud_fluid)

        if self.enable_true_buoyancy and self.enable_d2dga:
            # T1-3b: 体力向量注入流动度（式 2.5b/4.24），替换 (2φ−1) 简化代理
            beta_deg_local = float(np.mean(geom.get("inc_deg", np.zeros(self.nz))))
            # Task 4: F² 按 Z&F22 (2.6) 现算（原先内部硬编码 1.0 → 浮力缺席动力学），
            # 量级 O(10⁻³)（八井 [1.2e-3, 3.1e-2]）；ŵ₀ = q/A、μ̂₁ 取泥浆、d̂ = mean(geom["H"])。
            f2_local = self._froude_squared_at(geom, q_m3s, w_prev, mud_fluid)
            # Task 5: 领先阶轴向浮力数 b（Z&F22 p.8）按 summary 同口径现算（`_buoyancy_number_at`），
            # 经 `_mobility_profile` 接进 (4.14)/(4.22) 动力学——浮力第一次以 b 幅值参与 pref；
            # b>0（重顶替轻）⇒ 轴向浮力项抬高窄边流动度份额。
            b_num_local = self._buoyancy_number_at(geom, q_m3s, w_prev, mud_fluid,
                                                   lead_fluid, tail_fluid)
            rho_displaced = mud_fluid.density_kg_m3 / 1000.0
            delta_rho = (rho - rho_displaced)  # g/cc 局部密度差
            # Task 5: 流动度构造（幂律缝隙律 base + (4.23) 闭包 + (4.14)/(4.22) 浮力形状）
            # 抽为纯函数，便于单测；返回未饱和乘积 base·buoyancy_shape。
            mobility = self._mobility_profile(
                c_bar, b_num_local, geom, i1_base, m_local, beta_deg_local, f2_local,
                eta1=eta1, eta2=eta2, n_mix=n_mix, delta_rho=delta_rho,
            )
        else:
            # R0/R1/R2: 保留 buoyancy_shape 代理（旧论文状态）；base 构造保持原样
            # T9 DRY 收敛：base 构造收敛到 `_mobility_base`（运算逐字保持，逐位等价）
            base = self._mobility_base(c_bar, geom, i1_base, eta1, eta2, n_mix)
            phi = geom["phi"][:, None]
            ebar = geom["e"][None, :]
            density_contrast = (rho_disp - mud_fluid.density_kg_m3) / mud_fluid.density_kg_m3
            stable = float(np.clip(8.0 * density_contrast, -0.35, 0.45))
            buoyancy_shape = 1.0 + stable * ebar * (2.0 * phi - 1.0)
            mobility = base * buoyancy_shape
        pref = np.maximum(mobility, 1.0e-8)
        # T1-5: wall=1 处壁面静止层，流动度归零（式 2.35-2.41）
        if wall is not None:
            pref = pref * (1.0 - wall)

        # 由截面排量约束得到轴向速度 w（2026-09-02 守恒修正：截面权重去掉多余因子2，
        # 使同心极限 w=Q/A、全环空通量 2∫w b dy=Q；旧口径仅输运 Q/2 导致 η_E 系统偏低）。
        dy = np.gradient(geom["y"])[:, None]
        if self.enable_regime_split:
            # M2: 局部流态修正固定点迭代（Maleki & Frigaard 2017 式58-66）
            # 浓度相关量（mobility=base·buoyancy_shape/wall/黏度/密度/b）在迭代外缓存；
            # 黏度场保持 w_prev 一步滞后（既有约定），迭代只重算 Re_p/R/pref/area_weight/w。
            # rho_kg_m3 与 kappa_mix(Pa·s^n)/tau_y(Pa) 同单位系，保证 Re_p/He 无量纲正确。
            from cemdisp.models2d import regime_closure as rc
            wall_factor = np.ones_like(mobility) if wall is None else (1.0 - wall)
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
                # Task 5: base·buoyancy_shape 由 _mobility_profile 的返回值 mobility 承载
                pref_k = np.maximum(mobility * R_new, 1.0e-8) * wall_factor
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
            pref_final = np.maximum(mobility * R, 1.0e-8) * wall_factor
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
        # A-3b：τy 显式跳过留痕按 run 重置（R4 每跑一次性汇总告警，不跨 run 累积）；
        # warm-start 缓存一并重置（G 场形状绑定本 run 的网格，跨 run 不复用）。
        self._active_well_name = well_spec.well_name
        self._hb_tau_y_skips = []
        self._hb_tau_y_skips_reported = False
        self._hb_prev_G = None
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
                if self.enable_stream_function:
                    # T9 新路径：半拉格朗日平流用真实间隙平均速度——Z&F22 的输运
                    # (2.1) 以 (v̄,w̄) 输运浓度；(4.25) 通量形式里的 f_amp 是
                    # (4.28) 的**通量**函数而非速度乘子（把通量放大当速度乘子
                    # = B1 缺陷：模型凭空造水泥 27-35%，2026-09-10 取证）。
                    # f_amp 的速度层消费在新路径废止；f_amp 本身不再计算
                    # （新路径无消费点；q₀/I₃ 通量层不经此变量）。
                    # 旧路径（enable_stream_function=False）保留原行为逐位复现。
                    w_d2dga, v_d2dga = w, v
                else:
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

                # 2026-09-14 Task 7：自创拉普拉斯弥散已删除（Z&F22 p.11
                # "we have no diffusive terms"）。弥散由 q₀（平均平流通量）+ I₃
                # （浮力通量，式 4.25 第二项 / (4.26)）分层通量闭合承载。
                # 原设计"弥散 → 再次执行四相过填修正"的修补步一并删除：
                # 无弥散时平流后的过填修正已保证 sum≤1，第二次修正与第一次之间
                # 无任何浓度场修改，属同一步内的重复执行（no-op）。

                # R2: I3 浮力弥散通量（式 4.25 第二项）—— 仅作用于水泥相(lead+tail)
                if self.enable_d2dga_i3_flux and self.enable_d2dga:
                    cement_for_flux = np.clip(lead + tail, 0.0, 1.0)
                    # 浮力向量 f（用当前井段平均井斜）
                    beta_deg_local = float(np.mean(geom["inc_deg"])) if "inc_deg" in geom else 0.0
                    # Task 4: F² 按 Z&F22 (2.6) 现算，与 _compute_velocity 内同一口径
                    # （排量取最近一次有效泵注 q，速度取本步 w），量级 O(10⁻³)。
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
                    # ⚠️ Task 4 警示：该裁剪是**硬非线性**，**不得成为 I3 项的实际上限**。
                    # 当前量级下远未饱和（div_q ~ 10⁻⁵–10⁻⁴ /s ≪ step_limit ~ 40 /s，
                    # dt_step≈2s、ds≈10m、alpha_cfl=0.5），但若后续再放大浮力向量 f
                    # （或 dt_step 变小），I3 项会被这层 clip 吃掉，使 F² 定标在 R2 上失真；
                    # 届时须复核 div_q 是否触及 step_limit。
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
        # Task 5: 口径实现下沉到 `_buoyancy_number_at`（动力学逐时间步与 summary 共用，
        # 保证 b 的两处口径逐位一致）；半间隙 d̂ = mean(geom["H"])（R20 口径，b=2H 自洽）。
        b_number = self._buoyancy_number_at(
            geom, last_pump_rate_m3s, w_prev, mud_fluid, lead_fluid, tail_fluid,
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
