"""
套管内一维前沿追踪求解器

本模块实现套管内流体的一维前沿追踪计算，用于从地面开泵开始追踪流体在套管内下行的过程。

核心算法：
1. 解析计算各流体前缘到达鞋口的时间（体积推进法）
2. 查询任意时刻鞋口处正在流出的流体类型
3. 累计泵入体积超过管内容积时切换为后续流体

与环空二维模型耦合时，鞋口出流状态就是环空入口边界的物理来源：
在水泥前缘到达鞋口前，出流仍为初始管内钻井液；水泥到达后才切换为水泥。
若环空顶替时间不足，应延长环空二维求解总时长，而不是改变本1D前沿追踪逻辑。

主要类：
- CasingFlowResult: 套管内输运结果，包含：
    * fronts: 各流体前缘的最终位置和时间
    * schedule_steps: 施工步骤记录
    * pipe_cross_section_m2: 管内截面积
    * shoe_md_m: 鞋口深度

- CasingFlowSolver: 前沿追踪求解器，关键方法：
    * run(): 执行1D输运计算
    * pipe_exit_state_at(): 查询任意时刻鞋口出流状态

模型假设：
- 套管内流体以施工程序给定排量注入
- 流体前缘按体积推进法追踪，首版不考虑管内扩散
- 鞋口深度为从地面到鞋口的总测深
- 前缘到达鞋口后，对应流体从鞋口进入环空

物理修正模型（Phase 2 & 3）：
- 重力沉降修正：基于密度差驱动的前缘到达时间修正，支持井斜角投影
  [Ref: Romero & Carter, SPE 55927, 1999; Ekici et al., SPE 166112, 2013]
- 屈服应力修正：考虑流体屈服应力对重力沉降的抑制效应
  [Ref: Shah & Sutton, SPE 18036, 1990; Maleki & Frigaard, JFM 846, 2018]
- 停泵沉降增强：停泵期间考虑凝胶强度发展和屈服应力对沉降的综合影响
  [Ref: Kelessidis et al., JPT 2006; Maglione et al., SPE 56995, 1999]
- 轴向弥散：Taylor-Aris 型弥散模型，适用于层流条件
  [Ref: Taylor, Proc. R. Soc. A 219, 1953; Aris, Proc. R. Soc. A 235, 1956]
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep, PumpingStageEvent
from cemdisp.data.well_spec import WellSpec
from cemdisp.transport1d.interface_tracking import InterfaceFront
from cemdisp.transport1d.pipe_exit_state import PipeExitState
from cemdisp.transport1d.shoe_timeline import ShoeEvent, ShoeEventKind, ShoeTimeline


@dataclass(frozen=True)
class CasingFlowResult:
    """套管内一维输运结果。

    存储各流体前缘的最终位置和到达时间，以及管内几何参数。
    """

    fronts: tuple[InterfaceFront, ...] = field(default_factory=tuple)
    schedule_steps: tuple[PumpingScheduleStep, ...] = field(default_factory=tuple)
    pipe_cross_section_m2: float = 0.0
    shoe_md_m: float = 0.0
    pumping_end_time_s: float = 0.0
    cement_end_time_s: float = 0.0
    shoe_timeline: ShoeTimeline = field(default_factory=lambda: ShoeTimeline(events=[]))
    notes: tuple[str, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class _ScheduledStep:
    """内部使用的施工步骤时间窗。

    将施工步骤转换为带起止时间和累计体积的内部结构，便于时间查询。
    """

    step: PumpingScheduleStep
    start_time_s: float
    end_time_s: float
    cumulative_volume_start_m3: float
    cumulative_volume_end_m3: float


class CasingFlowSolver:
    """套管内一维前沿追踪求解器。

    模型假设：
    - 套管内流体以施工程序给定排量注入；
    - 流体前缘按体积推进法追踪，首版不考虑管内扩散；
    - 鞋口深度为从地面到鞋口的总测深；
    - 前缘到达鞋口后，对应流体从鞋口进入环空。

    物理修正（默认启用）：
    - 重力沉降修正（enable_gravity=True）：根据流体密度差调整前缘到达时间，
      支持井斜角投影修正和屈服应力修正。
    - 轴向弥散（enable_axial_dispersion=True）：Taylor-Aris 型弥散模型，
      将尖锐前缘转换为平滑过渡带。

    停泵沉降增强模型（phase 3）：
    - 停泵期间考虑凝胶强度发展（指数增长模型）和屈服应力对沉降的抑制。
    - 凝胶强度发展参考 Kelessidis et al. (JPT, 2006) 的实验数据拟合。

    Literature references:
    - Romero & Carter, SPE 55927 (1999): Gravity settling in inclined wells
    - Shah & Sutton, SPE 18036 (1990): Yield stress effects on settling
    - Taylor, Proc. R. Soc. A 219 (1953): Laminar dispersion
    - Kelessidis et al., JPT (2006): Gel strength development
    """

    def __init__(
        self,
        *,
        dt: float = 2.0,
        enable_gravity: bool = True,
        g_constant: float = 9.81,
        settling_velocity_factor: float = 0.0015,
        enable_axial_dispersion: bool = True,
        dispersion_alpha: float = 0.25,
        gelation_time_s: float = 600.0,
        gelation_max_factor: float = 0.95,
        enable_buoyancy_physics: bool = True,
        buoyancy_correction_factor: float = 1.0,
        enable_mixing_enhancement: bool = True,
        mixing_enhancement_factor: float = 5.0,
        max_mixing_enhancement: float = 10.0,
        has_plug: bool = False,
        enable_utube: bool = False,
    ) -> None:
        """初始化求解器。

        Args:
            dt: 时间步长（秒），仅用于内部时间查询的容差判断
            enable_gravity: 是否启用套管内重力修正；默认启用以提高物理真实性。
                基于密度差驱动的沉降模型，参考 Romero & Carter (SPE 55927, 1999)。
                可通过 settling_velocity_factor 调节修正幅度。
            g_constant: 重力加速度（m/s²），用于按现场重力条件缩放简化修正项
            settling_velocity_factor: 沉降速度系数，单位为 m/s 每 kg/m³ 密度差。
                默认 0.0015，参考现场数据标定（Ekici et al., SPE 166112, 2013）。
                仅在 enable_buoyancy_physics=False 的旧经验乘子路径中使用。
            enable_axial_dispersion: 是否启用管内轴向弥散；默认启用。
                基于 Taylor-Aris 弥散理论（Taylor 1953, Aris 1956），
                将尖锐流体前缘转换为平滑 S 形过渡带。
            dispersion_alpha: 无量纲弥散系数，默认 0.25。
                对应层流条件下的经验取值范围 [0.1, 0.5]。
            gelation_time_s: 凝胶强度发展特征时间（秒），默认 600（10 分钟）。
                用于停泵沉降增强模型，参考 Kelessidis et al. (JPT, 2006)。
            gelation_max_factor: 凝胶强度最大抑制因子，无量纲，默认 0.95。
                表示凝胶完全发展后对沉降的最大抑制程度（0~1）。
            enable_buoyancy_physics: 是否启用基于 Atwood 数的物理化浮力修正（T1-11）。
                默认启用；False 时回退到以水密度为基准的旧经验乘子（向后兼容）。
                参考 Dai et al. (2024) Petroleum Research。
            buoyancy_correction_factor: 浮力修正标定系数，无量纲，默认 1.0。
                需六井数据反标定（审查修正：0.1 时修正量与旧公式几乎相同）。
            enable_mixing_enhancement: 是否启用界面混浆增强分散（T1-10）。
                在密度不稳定界面（重驱轻 + 高 Re）对 Taylor-Aris 分散系数
                施加增强因子，使过渡带反映物理混浆效应。默认启用。
            mixing_enhancement_factor: 混浆增强系数 k_mix，无量纲，默认 5.0。
                经验值，需六井数据标定。
            max_mixing_enhancement: 混浆增强因子上限，无量纲，默认 10.0。
                防止过渡带不物理地过宽。
            has_plug: 是否有胶塞（尾管固井常配胶塞）；有则胶塞刮拭阻止混浆，
                混浆增强因子恒为 1，默认 False。
            enable_utube: U型管/自由下落修正开关（预留接口），默认关闭。
                当前为空实现钩子 _utube_corrected_arrival_time，未来实现（T1-12）。
        """
        if not math.isfinite(dt) or dt <= 0.0:
            raise ValueError("dt 必须为大于0的有限数值")
        if not math.isfinite(g_constant) or g_constant <= 0.0:
            raise ValueError("g_constant 必须为大于0的有限数值")
        if not math.isfinite(settling_velocity_factor) or settling_velocity_factor < 0.0:
            raise ValueError("settling_velocity_factor 必须为非负有限数值")
        if not math.isfinite(dispersion_alpha) or dispersion_alpha < 0.0:
            raise ValueError("dispersion_alpha 必须为非负有限数值")
        if not math.isfinite(gelation_time_s) or gelation_time_s <= 0.0:
            raise ValueError('gelation_time_s must be a positive finite number')
        if not math.isfinite(gelation_max_factor) or not (0.0 <= gelation_max_factor <= 1.0):
            raise ValueError('gelation_max_factor must be in [0, 1]')
        if not math.isfinite(buoyancy_correction_factor) or buoyancy_correction_factor < 0.0:
            raise ValueError("buoyancy_correction_factor 必须为非负有限数值")
        if not math.isfinite(mixing_enhancement_factor) or mixing_enhancement_factor < 0.0:
            raise ValueError("mixing_enhancement_factor 必须为非负有限数值")
        if not math.isfinite(max_mixing_enhancement) or max_mixing_enhancement < 1.0:
            raise ValueError("max_mixing_enhancement 必须为 >=1 的有限数值")
        self.dt: float = dt
        self.enable_gravity: bool = enable_gravity
        self.g_constant: float = g_constant
        self.settling_velocity_factor: float = settling_velocity_factor
        self.enable_axial_dispersion: bool = enable_axial_dispersion
        self.dispersion_alpha: float = dispersion_alpha
        self.gelation_time_s: float = gelation_time_s
        self.gelation_max_factor: float = gelation_max_factor
        self.enable_buoyancy_physics: bool = enable_buoyancy_physics
        self.buoyancy_correction_factor: float = buoyancy_correction_factor
        self.enable_mixing_enhancement: bool = enable_mixing_enhancement
        self.mixing_enhancement_factor: float = mixing_enhancement_factor
        self.max_mixing_enhancement: float = max_mixing_enhancement
        self.has_plug: bool = has_plug
        self.enable_utube: bool = enable_utube
        self._scheduled_steps_by_result_id: dict[int, tuple[_ScheduledStep, ...]] = {}
        self._initial_fluid_by_result_id: dict[int, str] = {}
        self._fluids_by_result_id: dict[int, tuple[FluidSpec, ...]] = {}

    def run(
        self,
        well_spec: WellSpec,
        fluids: tuple[FluidSpec, ...],
        schedule: PumpingSchedule,
    ) -> CasingFlowResult:
        """运行套管内1D前沿追踪。"""

        pipe_area_m2 = self._pipe_cross_section_area(well_spec)
        shoe_depth_m = well_spec.shoe_md_m
        pipe_volume_m3 = shoe_depth_m * pipe_area_m2
        scheduled_steps = self._build_scheduled_steps(schedule)
        initial_fluid = self._initial_fluid_name(fluids, schedule)
        fluid_by_name = {fluid.name: fluid for fluid in fluids}

        # 为每个注入步骤建立"前缘"：前缘位置由累计泵入体积推动，
        # 而不是只由该流体自身注入体积决定。
        fronts: list[InterfaceFront] = []
        for i, scheduled in enumerate(scheduled_steps):
            arrival_time_s = self._front_arrival_time(scheduled, scheduled_steps, pipe_volume_m3)
            if self.enable_gravity and arrival_time_s is not None:
                arrival_time_s = self._gravity_corrected_arrival_time(
                    arrival_time_s,
                    scheduled.step.fluid_name,
                    self._displaced_fluid_name(scheduled_steps, i, initial_fluid),
                    fluids,
                    well_spec,
                )
            final_distance_m = min(shoe_depth_m, scheduled.cumulative_volume_end_m3 / pipe_area_m2)
            fronts.append(
                InterfaceFront(
                    fluid_name=scheduled.step.fluid_name,
                    distance_m=final_distance_m,
                    time_s=arrival_time_s if arrival_time_s is not None else scheduled.end_time_s,
                )
            )

        # 水泥浆停止时刻定义为"最后一段水泥浆尾缘越过鞋口"的地面累计时间。
        # 这与参考项目的环空终止口径一致：后续替浆刚到环空入口时停止二维顶替评价。
        cement_end_time_s: float | None = None
        for i, scheduled in enumerate(scheduled_steps):
            fluid = fluid_by_name.get(scheduled.step.fluid_name)
            if fluid is None or fluid.role not in {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}:
                continue
            rear_arrival_time_s = self._rear_arrival_time(scheduled, scheduled_steps, pipe_volume_m3)
            if rear_arrival_time_s is not None:
                if self.enable_gravity:
                    rear_arrival_time_s = self._gravity_corrected_arrival_time(
                        rear_arrival_time_s,
                        scheduled.step.fluid_name,
                        self._displaced_fluid_name(scheduled_steps, i, initial_fluid),
                        fluids,
                        well_spec,
                    )
                cement_end_time_s = rear_arrival_time_s

        result = CasingFlowResult(
            fronts=tuple(fronts),
            schedule_steps=schedule.steps,
            pipe_cross_section_m2=pipe_area_m2,
            shoe_md_m=shoe_depth_m,
            pumping_end_time_s=scheduled_steps[-1].end_time_s if scheduled_steps else 0.0,
            cement_end_time_s=cement_end_time_s if cement_end_time_s is not None else (scheduled_steps[-1].end_time_s if scheduled_steps else 0.0),
            shoe_timeline=self._build_shoe_timeline(
                well_spec=well_spec,
                fluids=fluids,
                scheduled_steps=scheduled_steps,
                initial_fluid=initial_fluid,
                legacy_pipe_volume_m3=pipe_volume_m3,
            ),
            notes=(
                "套管内采用理想界面前缘追踪，未加入管内扩散。",
                f"初始管内流体按 {initial_fluid} 处理。",
            ),
        )
        self._scheduled_steps_by_result_id[id(result)] = scheduled_steps
        self._initial_fluid_by_result_id[id(result)] = initial_fluid
        self._fluids_by_result_id[id(result)] = fluids
        return result

    def pipe_exit_state_at(self, result: CasingFlowResult, time_s: float) -> PipeExitState:
        """给定时间，返回鞋口出流状态。"""

        if not math.isfinite(time_s) or time_s < 0.0:
            raise ValueError("time_s 必须为非负有限数值")

        scheduled_steps = self._scheduled_steps_for_result(result)
        pipe_volume_m3 = result.shoe_md_m * result.pipe_cross_section_m2
        initial_fluid = self._initial_fluid_by_result_id.get(id(result), "初始管内流体")
        state = self._pipe_exit_state_from_volume(
            scheduled_steps=scheduled_steps,
            time_s=time_s,
            pipe_volume_m3=pipe_volume_m3,
            initial_fluid=initial_fluid,
        )

        if self.enable_gravity and state.flow_rate_m3_s < 1.0e-9:
            cumulative_at_time = self._cumulative_volume_at(scheduled_steps, time_s)
            fluid_name = self._settled_exit_fluid_name_enhanced(
                result=result,
                scheduled_steps=scheduled_steps,
                time_s=time_s,
                cumulative_at_time=cumulative_at_time,
                pipe_volume_m3=pipe_volume_m3,
                default_fluid_name=state.phase_fractions[0][0],
            )
            return PipeExitState(
                time_s=time_s,
                flow_rate_m3_s=state.flow_rate_m3_s,
                stage_name=state.stage_name,
                phase_fractions=((fluid_name, 1.0),),
            )

        return state

    def _compute_dispersion_coefficient(
        self,
        pipe_radius_m: float,
        fluid: FluidSpec,
        mean_velocity_m_s: float,
    ) -> float:
        """管内层流轴向弥散系数（Taylor-Aris + 屈服应力修正）。

        T1-8: Newtonian → Taylor-Aris; Power-law → Batot et al. (2016) 式(28);
        Bingham → Fan & Wang (1966); HB → 等效 Bingham 近似。

        Args:
            pipe_radius_m: 套管内半径 [m]
            fluid: 管内流体规格
            mean_velocity_m_s: 截面平均速度 [m/s]

        Returns:
            有效轴向弥散系数 D_eff [m²/s]（仅对流弥散部分）
        """
        if mean_velocity_m_s < 1e-9:
            return 0.0

        import math as _math
        from cemdisp.data.fluid_spec import RheologyModel

        U = mean_velocity_m_s
        R = pipe_radius_m
        d_mol = 1.0e-9  # 分子扩散系数 [m²/s]
        Pe = 2.0 * U * R / d_mol  # Péclet 数

        if fluid.rheology_model == RheologyModel.NEWTONIAN:
            return d_mol * Pe**2 / 192.0

        elif fluid.rheology_model == RheologyModel.POWER_LAW:
            n = fluid.power_law_n if fluid.power_law_n is not None else 1.0
            # Batot et al. (2016) 式(28)
            factor = 24.0 * n**2 / ((3.0 * n + 1.0) * (5.0 * n + 1.0))
            return d_mol * Pe**2 / 192.0 * factor

        elif fluid.rheology_model in (RheologyModel.BINGHAM, RheologyModel.HERSCHEL_BULKLEY):
            # T1-8: Fan & Wang (1966) Bingham Taylor 弥散
            tau_y = fluid.yield_stress_pa if fluid.yield_stress_pa is not None else 0.0
            if fluid.rheology_model == RheologyModel.HERSCHEL_BULKLEY:
                k_cons = fluid.consistency_k if fluid.consistency_k is not None else 0.01
                n_hb = fluid.power_law_n if fluid.power_law_n is not None else 1.0
                shear_rate_w = max(8.0 * U / (2.0 * R), 1e-6)
                mu_app = tau_y / max(shear_rate_w, 1e-9) + k_cons * shear_rate_w**(n_hb - 1.0)
            else:
                mu_app = fluid.plastic_viscosity_pa_s if fluid.plastic_viscosity_pa_s is not None else 0.01

            xi0 = tau_y * R / max(4.0 * mu_app * U, 1e-12)
            xi0 = min(xi0, 0.999)

            if xi0 < 1e-9:
                return d_mol * Pe**2 / 192.0

            # Fan & Wang (1966) k(ξ₀) 多项式
            xi = xi0; xi2 = xi**2; xi4 = xi2**2; xi8 = xi4**2
            num = (3.0/8.0 - 44.0/35.0*xi + 16.0/15.0*xi2 + xi4
                   - 28.0/15.0*xi4*xi - 3.0/5.0*xi4*xi2 + 8.0/5.0*xi4*xi2*xi
                   - 29.0/56.0*xi8 + 1.0/5.0*xi8*xi2 - xi8*_math.log(max(xi, 1e-12)))
            den = 2.0 * (3.0 + 2.0*xi + xi2) * (1.0 - xi)**4
            k = max(num / max(den, 1e-12), 0.0)

            return d_mol * Pe**2 / 48.0 * k

        else:
            return self.dispersion_alpha * mean_velocity_m_s * pipe_radius_m

    def _effective_viscosity(
        self,
        fluid: FluidSpec,
        mean_velocity_m_s: float,
        pipe_radius_m: float,
    ) -> float:
        """计算流体在给定剪切率下的表观黏度。

        使用 Dai 2024 eq. A.11: μ_eff = τ₀/(u/D) + k·(u/D)^(n-1)，
        剪切率取壁面剪切率近似 γ_w = 8U/(2R)。
        牛顿流体直接返回塑性黏度；无法求值的边界（零流速/零半径）回退到塑性黏度。

        Args:
            fluid: 流体规格
            mean_velocity_m_s: 截面平均速度 [m/s]
            pipe_radius_m: 管内半径 [m]

        Returns:
            表观黏度 [Pa·s]，恒为正数
        """
        if mean_velocity_m_s < 1e-9 or pipe_radius_m < 1e-9:
            return fluid.plastic_viscosity_pa_s or 0.01

        shear_rate = 8.0 * mean_velocity_m_s / (2.0 * pipe_radius_m)  # 壁面剪切率近似
        shear_rate = max(shear_rate, 1e-6)

        if fluid.rheology_model == RheologyModel.NEWTONIAN:
            return fluid.plastic_viscosity_pa_s or 0.01

        tau_y = fluid.yield_stress_pa or 0.0
        k_cons = fluid.consistency_k or 0.01
        n = fluid.power_law_n if fluid.power_law_n is not None else 1.0

        return tau_y / shear_rate + k_cons * shear_rate ** (n - 1.0)

    def _interface_instability_factor(
        self,
        fluid_next: FluidSpec,
        fluid_prev: FluidSpec,
        pipe_radius_m: float,
        mean_velocity_m_s: float,
    ) -> float:
        """界面失稳增强因子（垂直井适配）。

        返回 >1 表示失稳增强分散，=1 表示稳定（仅 Taylor-Aris）。
        受 Dai 2024 启发的垂直井适配判据。注意：Dai 2024 原判据（eq. A.10/A.14/A.15）
        依赖 cos(β) 和 Fr，专为斜井设计；垂直井（β≈0）时 vt→0、Fr→∞，原判据不可直接用。
        此处简化为 At>0 AND Re>100 的密度不稳定判据。

        物理依据：
        - 重驱轻（At > 0）→ Rayleigh-Taylor 型界面失稳 → 混合增强
        - 黏度差越大 → 界面越不稳定
        - 有胶塞时 → 胶塞刮拭阻止混合 → 增强因子=1

        Args:
            fluid_next: 后继流体（当前注入流体）
            fluid_prev: 前置流体（被顶替流体）
            pipe_radius_m: 管内半径 [m]
            mean_velocity_m_s: 截面平均速度 [m/s]

        Returns:
            界面失稳增强因子 [1, max_mixing_enhancement]
        """
        if self.has_plug:
            return 1.0

        rho_h = max(fluid_next.density_kg_m3, fluid_prev.density_kg_m3)
        rho_l = min(fluid_next.density_kg_m3, fluid_prev.density_kg_m3)
        at = (rho_h - rho_l) / (rho_h + rho_l)  # Atwood number

        if at < 1e-9:
            return 1.0  # 等密度，无失稳

        # 有效黏度（几何平均，Dai 2024 eq. A.12）
        mu_next = self._effective_viscosity(fluid_next, mean_velocity_m_s, pipe_radius_m)
        mu_prev = self._effective_viscosity(fluid_prev, mean_velocity_m_s, pipe_radius_m)
        mu_mean = math.sqrt(max(mu_next * mu_prev, 1e-12))

        # Reynolds number
        rho_avg = (fluid_next.density_kg_m3 + fluid_prev.density_kg_m3) / 2.0
        D = 2.0 * pipe_radius_m
        Re = rho_avg * mean_velocity_m_s * D / max(mu_mean, 1e-12)

        # 垂直井密度不稳定判据：重驱轻 + Re 足够大 → 失稳
        if at > 0 and Re > 100.0:
            enhancement = 1.0 + self.mixing_enhancement_factor * at * math.sqrt(Re / 100.0)
            return min(enhancement, self.max_mixing_enhancement)
        return 1.0

    def _apply_dispersion_to_timeline(
        self,
        events: list[ShoeEvent],
        well_spec: WellSpec,
        scheduled_steps: tuple[_ScheduledStep, ...],
        fluids: tuple[FluidSpec, ...],
    ) -> list[ShoeEvent]:
        """对离散鞋口时间线施加轴向弥散。

        在每个流体前缘 (FRONT_ARRIVAL) 附近，用弥散宽度 σ_t 生成过渡态事件，
        将阶跃切换转换为平滑的 S 形过渡带。
        """
        if not self.enable_axial_dispersion:
            return events

        pipe_radius_m = (well_spec.liner_id_mm or 100.0) / 2000.0
        fluid_by_name = {f.name: f for f in fluids}

        dispersed_events: list[ShoeEvent] = []

        for i, event in enumerate(events):
            if event.kind != ShoeEventKind.FRONT_ARRIVAL:
                dispersed_events.append(event)
                continue

            fluid_name = event.phase_fractions[0][0] if event.phase_fractions else ""
            fluid = fluid_by_name.get(fluid_name)
            if fluid is None:
                dispersed_events.append(event)
                continue

            U = event.flow_rate_m3_s / (math.pi * pipe_radius_m ** 2)
            D_eff = self._compute_dispersion_coefficient(pipe_radius_m, fluid, U)

            # 找到前一个流体（在弥散系数计算前确定，供混浆增强注入使用）
            prev_fluid = ""
            for j in range(i - 1, -1, -1):
                if events[j].phase_fractions:
                    prev_fluid = events[j].phase_fractions[0][0]
                    break

            # T1-10: 混浆增强分散（仅当前流体与前置流体均已知时注入）
            if self.enable_mixing_enhancement:
                prev_fluid_spec = fluid_by_name.get(prev_fluid)
                if prev_fluid_spec is not None:
                    instability = self._interface_instability_factor(
                        fluid, prev_fluid_spec, pipe_radius_m, U
                    )
                    D_eff *= instability

            if D_eff < 1e-12:
                dispersed_events.append(event)
                continue

            # 到达时间（活塞流）
            t_arrival = event.time_s
            # 弥散时间宽度: σ_t = sqrt(2 × D_eff × t_travel) / U
            # t_travel = pipe_volume / Q ≈ shoe_md_m / U
            t_travel = well_spec.shoe_md_m / max(U, 1e-9)
            sigma_t = math.sqrt(2.0 * D_eff * t_travel) / max(U, 1e-9)
            sigma_t = max(sigma_t, self.dt)  # 至少一个时间步

            next_fluid = fluid_name

            # 在 [t_arrival - σ, t_arrival + σ] 范围内生成过渡子事件
            n_sub = 5  # 每个过渡带的子事件数
            for k in range(n_sub):
                t_sub = t_arrival + sigma_t * (2.0 * k / (n_sub - 1) - 1.0)  # [-σ, +σ]
                frac = 0.5 * (1.0 + math.erf(k / (n_sub - 1.0) * 2.0 - 1.0))  # erf 过渡

                dispersed_events.append(ShoeEvent(
                    time_s=t_sub,
                    kind=ShoeEventKind.FRONT_ARRIVAL,
                    flow_rate_m3_s=event.flow_rate_m3_s,
                    stage_name=event.stage_name,
                    phase_fractions=(
                        (next_fluid, frac),
                        (prev_fluid, 1.0 - frac),
                    ),
                ))

        return sorted(dispersed_events, key=lambda e: e.time_s)

    def _build_shoe_timeline(
        self,
        *,
        well_spec: WellSpec,
        fluids: tuple[FluidSpec, ...],
        scheduled_steps: tuple[_ScheduledStep, ...],
        initial_fluid: str,
        legacy_pipe_volume_m3: float,
    ) -> ShoeTimeline:
        """按现有体积追踪数学生成鞋口出流事件时间轴。

        时间轴只记录状态发生变化的关键时刻：施工步骤切换、流体前缘到达、
        流体尾缘离开以及施工结束；不引入新的扩散或混合物理模型。
        """

        pipe_volume_m3 = self._timeline_pipe_volume(well_spec, legacy_pipe_volume_m3)
        event_points: list[tuple[float, ShoeEventKind, tuple[tuple[str, float], ...] | None]] = []
        for i, scheduled in enumerate(scheduled_steps):
            displaced_fluid = self._displaced_fluid_name(scheduled_steps, i, initial_fluid)
            event_points.append((scheduled.start_time_s, self._event_kind_for_step(scheduled.step), None))
            front_time_s = self._front_arrival_time(scheduled, scheduled_steps, pipe_volume_m3)
            if front_time_s is not None:
                if self.enable_gravity:
                    front_time_s = self._gravity_corrected_arrival_time(
                        front_time_s, scheduled.step.fluid_name, displaced_fluid, fluids, well_spec
                    )
                # T1-12 预留：U型管/自由下落修正钩子（重力修正之后调用）
                front_time_s = self._utube_corrected_arrival_time(
                    front_time_s, scheduled.step.fluid_name, fluids, well_spec
                )
                event_points.append((front_time_s, ShoeEventKind.FRONT_ARRIVAL, ((scheduled.step.fluid_name, 1.0),)))
            rear_time_s = self._rear_arrival_time(scheduled, scheduled_steps, pipe_volume_m3)
            if rear_time_s is not None:
                if self.enable_gravity:
                    rear_time_s = self._gravity_corrected_arrival_time(
                        rear_time_s, scheduled.step.fluid_name, displaced_fluid, fluids, well_spec
                    )
                # T1-12 预留：U型管/自由下落修正钩子（重力修正之后调用）
                rear_time_s = self._utube_corrected_arrival_time(
                    rear_time_s, scheduled.step.fluid_name, fluids, well_spec
                )
                event_points.append((rear_time_s, ShoeEventKind.REAR_EXIT, None))
        if scheduled_steps:
            event_points.append((scheduled_steps[-1].end_time_s, ShoeEventKind.END, None))

        events: list[ShoeEvent] = []
        for time_s, kind, phase_override in sorted(event_points, key=lambda item: item[0]):
            state = self._pipe_exit_state_from_volume(
                scheduled_steps=scheduled_steps,
                time_s=time_s,
                pipe_volume_m3=pipe_volume_m3,
                initial_fluid=initial_fluid,
            )
            events.append(
                ShoeEvent(
                    time_s=time_s,
                    kind=kind,
                    flow_rate_m3_s=state.flow_rate_m3_s,
                    stage_name=state.stage_name,
                    phase_fractions=phase_override if phase_override is not None else state.phase_fractions,
                )
            )
        # 应用轴向弥散处理：将尖锐的阶跃前缘转换为平滑的 S 形过渡带
        events = self._apply_dispersion_to_timeline(events, well_spec, scheduled_steps, fluids)
        return ShoeTimeline(events=events)

    @staticmethod
    def _event_kind_for_step(step: PumpingScheduleStep) -> ShoeEventKind:
        """把地面施工事件标签映射为鞋口时间轴事件类型。"""

        if step.event_tag == PumpingStageEvent.SHUTDOWN:
            return ShoeEventKind.SHUTDOWN
        if step.event_tag == PumpingStageEvent.RESTART:
            return ShoeEventKind.RESTART
        return ShoeEventKind.RATE_SWITCH

    @staticmethod
    def _pipe_exit_state_from_volume(
        *,
        scheduled_steps: tuple[_ScheduledStep, ...],
        time_s: float,
        pipe_volume_m3: float,
        initial_fluid: str,
    ) -> PipeExitState:
        """用体积迟到量查询鞋口出流状态，复用旧求解器的核心口径。"""

        active_step = CasingFlowSolver._active_step_at(scheduled_steps, time_s)
        flow_rate_m3_s = 0.0 if active_step is None else active_step.step.rate_m3_min / 60.0
        stage_name = "施工结束后保持" if active_step is None else active_step.step.step_name

        # 鞋口流出流体由"累计泵入体积 - 管内容积"定位到地面注入体积坐标。
        cumulative_at_time = CasingFlowSolver._cumulative_volume_at(scheduled_steps, time_s)
        fluid_name = initial_fluid
        if cumulative_at_time >= pipe_volume_m3:
            delayed_volume_m3 = cumulative_at_time - pipe_volume_m3
            fluid_name = CasingFlowSolver._fluid_by_injected_volume(scheduled_steps, delayed_volume_m3)

        return PipeExitState(
            time_s=time_s,
            flow_rate_m3_s=flow_rate_m3_s,
            stage_name=stage_name,
            phase_fractions=((fluid_name, 1.0),),
        )

    @staticmethod
    def _timeline_pipe_volume(well_spec: WellSpec, legacy_pipe_volume_m3: float) -> float:
        """计算时间轴使用的鞋口迟到体积，兼容双径向井上段内径。"""

        if well_spec.shoe_lag_volume_m3 is not None:
            return well_spec.shoe_lag_volume_m3
        if not well_spec.is_dual_diameter:
            return legacy_pipe_volume_m3
        assert well_spec.upper_section_bottom_md_m is not None
        assert well_spec.upper_liner_id_mm is not None
        assert well_spec.liner_id_mm is not None
        upper_area_m2 = math.pi * (well_spec.upper_liner_id_mm / 1000.0) ** 2 / 4.0
        lower_area_m2 = math.pi * (well_spec.liner_id_mm / 1000.0) ** 2 / 4.0
        upper_length_m = min(well_spec.upper_section_bottom_md_m, well_spec.shoe_md_m)
        lower_length_m = max(well_spec.shoe_md_m - well_spec.upper_section_bottom_md_m, 0.0)
        return upper_length_m * upper_area_m2 + lower_length_m * lower_area_m2

    @staticmethod
    def _pipe_cross_section_area(well_spec: WellSpec) -> float:
        """根据管内径计算等效截面积。

        优先级：shoe_lag_volume_m3 > pipe_id_profile > liner_id_mm。
        shoe_lag_volume_m3 是从地面到鞋口的管内迟到体积（m³），
        除以鞋口深度得到等效截面积。
        """
        if well_spec.shoe_lag_volume_m3 is not None:
            return well_spec.shoe_lag_volume_m3 / well_spec.shoe_md_m

        if well_spec.pipe_id_profile:
            import numpy as np
            depths = np.array([p.depth_md_m for p in well_spec.pipe_id_profile])
            ids_mm = np.array([p.value for p in well_spec.pipe_id_profile])
            mask = (depths >= 0.0) & (depths <= well_spec.shoe_md_m)
            d = depths[mask]
            ids = ids_mm[mask]
            if len(d) < 2:
                d = np.array([0.0, well_spec.shoe_md_m])
                ids = np.array([ids_mm[0], ids_mm[-1]])
            ids_m = ids / 1000.0
            areas = math.pi * ids_m**2 / 4.0
            total_volume = float(np.trapezoid(areas, x=d))
            return total_volume / well_spec.shoe_md_m

        if well_spec.liner_id_mm is None:
            raise ValueError("well_spec.liner_id_mm 不能为空，套管内1D模型需要内径")
        liner_id_m = well_spec.liner_id_mm / 1000.0
        return math.pi * liner_id_m**2 / 4.0

    @staticmethod
    def _build_scheduled_steps(schedule: PumpingSchedule) -> tuple[_ScheduledStep, ...]:
        """把施工步骤转成带起止时间和累计体积的内部结构。"""

        scheduled_steps: list[_ScheduledStep] = []
        elapsed_s = 0.0
        cumulative_volume_m3 = 0.0
        for step in schedule.steps:
            duration_s = 0.0 if step.rate_m3_min <= 0.0 else step.volume_m3 / step.rate_m3_min * 60.0
            start_time_s = step.start_time_s if step.start_time_s is not None else elapsed_s
            end_time_s = step.end_time_s if step.end_time_s is not None else start_time_s + duration_s
            scheduled_steps.append(
                _ScheduledStep(
                    step=step,
                    start_time_s=start_time_s,
                    end_time_s=end_time_s,
                    cumulative_volume_start_m3=cumulative_volume_m3,
                    cumulative_volume_end_m3=cumulative_volume_m3 + step.volume_m3,
                )
            )
            elapsed_s = end_time_s
            cumulative_volume_m3 += step.volume_m3
        return tuple(scheduled_steps)

    @staticmethod
    def _initial_fluid_name(fluids: tuple[FluidSpec, ...], schedule: PumpingSchedule) -> str:
        """确定开泵前套管内默认流体名称。"""

        if not schedule.steps:
            raise ValueError("施工程序不能为空，至少需要一个注入步骤")
        mud = next((fluid for fluid in fluids if fluid.role == FluidRole.MUD), None)
        if mud is not None:
            return mud.name
        return schedule.steps[0].fluid_name

    @staticmethod
    def _front_arrival_time(
        front_step: _ScheduledStep,
        scheduled_steps: tuple[_ScheduledStep, ...],
        pipe_volume_m3: float,
    ) -> float | None:
        """计算某个流体前缘到达鞋口的时间。"""

        target_volume_m3 = front_step.cumulative_volume_start_m3 + pipe_volume_m3
        for scheduled in scheduled_steps:
            if target_volume_m3 <= scheduled.cumulative_volume_end_m3 + 1.0e-12:
                if scheduled.step.rate_m3_min <= 0.0:
                    return scheduled.end_time_s
                volume_into_step_m3 = max(target_volume_m3 - scheduled.cumulative_volume_start_m3, 0.0)
                return scheduled.start_time_s + volume_into_step_m3 / scheduled.step.rate_m3_min * 60.0
        return None

    @staticmethod
    def _rear_arrival_time(
        rear_step: _ScheduledStep,
        scheduled_steps: tuple[_ScheduledStep, ...],
        pipe_volume_m3: float,
    ) -> float | None:
        """计算某个流体尾缘到达鞋口的时间。"""

        target_volume_m3 = rear_step.cumulative_volume_end_m3 + pipe_volume_m3
        for scheduled in scheduled_steps:
            if target_volume_m3 <= scheduled.cumulative_volume_end_m3 + 1.0e-12:
                if scheduled.step.rate_m3_min <= 0.0:
                    return scheduled.end_time_s
                volume_into_step_m3 = max(target_volume_m3 - scheduled.cumulative_volume_start_m3, 0.0)
                return scheduled.start_time_s + volume_into_step_m3 / scheduled.step.rate_m3_min * 60.0
        return None

    def _scheduled_steps_for_result(self, result: CasingFlowResult) -> tuple[_ScheduledStep, ...]:
        """获取结果对应的内部时间窗；必要时从冻结步骤重建。"""

        scheduled_steps = self._scheduled_steps_by_result_id.get(id(result))
        if scheduled_steps is not None:
            return scheduled_steps
        return self._build_scheduled_steps(PumpingSchedule(steps=result.schedule_steps))

    def _gravity_corrected_arrival_time(
        self,
        arrival_time_s: float,
        fluid_name: str,
        displaced_fluid_name: str,
        fluids: tuple[FluidSpec, ...],
        well_spec: WellSpec | None = None,
    ) -> float:
        """按流体密度对前缘到达时间做重力修正，支持井斜角投影和屈服应力修正。

        T1-11 物理化修正流程（enable_buoyancy_physics=True，默认）：
        1. 基础浮力修正：基于 Atwood 数（无量纲密度差，密度差基准为被顶替流体，
           而非旧公式中以水密度为基准），参考 Dai et al. (2024) Petroleum Research；
        2. 井斜角投影修正：当 well_spec.inclination_profile 可用时，
           用平均井斜角的余弦值投影重力分量（垂直分量 = g * cos(theta)）；
        3. 屈服应力修正：对有屈服应力的流体，修正沉降速度以反映屈服应力的抑制效应，
           参考 Maleki & Frigaard, JFM 846 (2018)；
        4. 方向判断：重驱轻（ρ_fluid > ρ_displaced）→ 浮力加速 → 到达时间缩短；
           轻驱重 → 浮力减速 → 到达时间延长。

        enable_buoyancy_physics=False 时回退到旧经验乘子（_legacy_gravity_correction），
        保持向后兼容。

        Literature:
        - Dai et al., Petroleum Research (2024): Atwood 数表征界面密度差
        - Romero & Carter, SPE 55927 (1999): Gravity settling in inclined wells
        - Shah & Sutton, SPE 18036 (1990): Yield stress effects on settling
        - Maleki & Frigaard, JFM 846 (2018): 屈服应力抑制浮力效应

        Args:
            arrival_time_s: 活塞流模型计算的前缘到达时间（秒）
            fluid_name: 流体名称（当前注入流体）
            displaced_fluid_name: 被顶替流体名称（上一步流体，或初始管内流体）
            fluids: 全部流体规格元组
            well_spec: 井筒规格，用于获取井斜剖面和管内径

        Returns:
            修正后的前缘到达时间（秒），保证非负
        """
        if not self.enable_buoyancy_physics:
            # 回退到旧经验乘子（向后兼容）
            return self._legacy_gravity_correction(arrival_time_s, fluid_name, fluids, well_spec)

        rho_fluid = self._get_fluid_density(fluid_name, fluids)
        rho_displaced = self._get_fluid_density(displaced_fluid_name, fluids)

        # Atwood 数（无量纲密度差）
        at = abs(rho_fluid - rho_displaced) / (rho_fluid + rho_displaced)

        # 浮力修正因子
        gravity_scale = self.g_constant / 9.81
        gravity_factor = self.buoyancy_correction_factor * at * gravity_scale

        # 井斜角投影修正：重力沿管轴分量 = g * cos(inclination)
        if well_spec is not None and well_spec.inclination_profile:
            avg_inclination_rad = self._average_inclination_rad(well_spec)
            gravity_factor *= max(math.cos(avg_inclination_rad), 0.0)

        # 屈服应力修正：屈服应力会抑制浮力滑移，减小有效浮力效应
        fluid = next((f for f in fluids if f.name == fluid_name), None)
        if fluid is not None and fluid.yield_stress_pa is not None and fluid.yield_stress_pa > 0.0:
            pipe_radius_m = self._effective_pipe_radius_m(well_spec)
            delta_rho = abs(rho_fluid - rho_displaced)
            tau_critical = delta_rho * self.g_constant * pipe_radius_m
            if tau_critical > 1e-6:
                yield_ratio = min(fluid.yield_stress_pa / tau_critical, 1.0)
                gravity_factor *= (1.0 - 0.8 * yield_ratio)

        # 方向判断：重驱轻 → 浮力加速；轻驱重 → 浮力减速
        if rho_fluid > rho_displaced:
            return max(arrival_time_s * (1.0 - gravity_factor), 0.0)
        else:
            return arrival_time_s * (1.0 + gravity_factor)

    def _legacy_gravity_correction(
        self,
        arrival_time_s: float,
        fluid_name: str,
        fluids: tuple[FluidSpec, ...],
        well_spec: WellSpec | None = None,
    ) -> float:
        """旧经验乘子重力修正（向后兼容路径）。

        enable_buoyancy_physics=False 时使用：以水密度 1000 kg/m³ 为密度差基准的
        经验乘子（settling_velocity_factor × (ρ - 1000) / 1000），保留井斜角投影
        和屈服应力抑制逻辑。

        Args:
            arrival_time_s: 活塞流模型计算的前缘到达时间（秒）
            fluid_name: 流体名称
            fluids: 全部流体规格元组
            well_spec: 井筒规格，用于获取井斜剖面和管内径

        Returns:
            修正后的前缘到达时间（秒），保证非负
        """
        fluid_density = self._get_fluid_density(fluid_name, fluids)
        fluid = next((f for f in fluids if f.name == fluid_name), None)

        # 基础重力因子
        gravity_scale = self.g_constant / 9.81
        gravity_factor = self.settling_velocity_factor * (fluid_density - 1000.0) / 1000.0 * gravity_scale

        # 井斜角投影修正：重力沿管轴分量 = g * cos(inclination)
        if well_spec is not None and well_spec.inclination_profile:
            avg_inclination_rad = self._average_inclination_rad(well_spec)
            cos_inclination = math.cos(avg_inclination_rad)
            gravity_factor *= max(cos_inclination, 0.0)

        # 屈服应力修正：屈服应力会抑制沉降，减小有效沉降速度
        if fluid is not None and fluid.yield_stress_pa is not None and fluid.yield_stress_pa > 0.0:
            pipe_radius_m = self._effective_pipe_radius_m(well_spec)
            delta_rho = abs(fluid_density - 1000.0)
            tau_critical = delta_rho * self.g_constant * pipe_radius_m
            if tau_critical > 1e-6:
                yield_ratio = min(fluid.yield_stress_pa / tau_critical, 1.0)
                yield_suppression = 1.0 - 0.8 * yield_ratio
                gravity_factor *= yield_suppression

        return max(arrival_time_s * (1.0 - gravity_factor), 0.0)

    @staticmethod
    def _displaced_fluid_name(
        scheduled_steps: tuple[_ScheduledStep, ...],
        current_index: int,
        initial_fluid: str,
    ) -> str:
        """返回当前步骤的被顶替流体名（上一步流体，或初始泥浆）。

        Args:
            scheduled_steps: 内部施工步骤时间窗
            current_index: 当前步骤在时间窗中的索引
            initial_fluid: 开泵前套管内默认流体名称（初始泥浆）

        Returns:
            被顶替流体名称；首步（current_index <= 0）回退到初始泥浆。
        """
        if current_index <= 0:
            return initial_fluid
        return scheduled_steps[current_index - 1].step.fluid_name

    def _utube_corrected_arrival_time(
        self,
        arrival_time_s: float,
        fluid_name: str,
        fluids: tuple[FluidSpec, ...],
        well_spec: WellSpec | None = None,
    ) -> float:
        """U型管/自由下落修正钩子。当前为空实现，预留未来扩展。

        未来实现方向（T1-12）：
        1. 计算套管-环空静压平衡
        2. 检测自由下落条件（ΔP_hydro > ΔP_friction）
        3. 修正到达时间（加速效应）

        Args:
            arrival_time_s: 重力修正后的前缘到达时间（秒）
            fluid_name: 流体名称
            fluids: 全部流体规格元组
            well_spec: 井筒规格

        Returns:
            修正后的到达时间（秒）；enable_utube=False 或未实现时原样返回
        """
        if not self.enable_utube:
            return arrival_time_s
        # 未来：填入 U 型管物理计算
        return arrival_time_s

    def _average_inclination_rad(self, well_spec: WellSpec) -> float:
        """计算井段内平均井斜角（弧度）。

        从 well_spec.inclination_profile 提取井斜数据，按深度加权平均。
        若无井斜数据，假设为垂直井（0 弧度）。

        Args:
            well_spec: 井筒规格

        Returns:
            平均井斜角（弧度），范围 [0, pi/2]
        """
        if not well_spec.inclination_profile:
            return 0.0

        points = well_spec.inclination_profile
        if len(points) == 1:
            return math.radians(points[0].value)

        total_weighted = 0.0
        total_length = 0.0
        for i in range(len(points) - 1):
            p0 = points[i]
            p1 = points[i + 1]
            segment_length = p1.depth_md_m - p0.depth_md_m
            if segment_length <= 0.0:
                continue
            avg_inclination_deg = (p0.value + p1.value) / 2.0
            total_weighted += avg_inclination_deg * segment_length
            total_length += segment_length

        if total_length <= 0.0:
            return math.radians(points[0].value)

        return math.radians(total_weighted / total_length)

    def _effective_pipe_radius_m(self, well_spec: WellSpec | None = None) -> float:
        """获取有效管内半径（米），用于屈服应力修正计算。

        优先使用 well_spec.liner_id_mm，若不可用则使用默认值。

        Args:
            well_spec: 井筒规格

        Returns:
            有效管内半径（米）
        """
        if well_spec is not None and well_spec.liner_id_mm is not None:
            return well_spec.liner_id_mm / 2000.0
        return 0.05

    def _settled_exit_fluid_name(
        self,
        *,
        result: CasingFlowResult,
        scheduled_steps: tuple[_ScheduledStep, ...],
        time_s: float,
        cumulative_at_time: float,
        pipe_volume_m3: float,
        default_fluid_name: str,
    ) -> str:
        """停泵期间按密度差估算沉降后鞋口流体。"""

        fluids = self._fluids_by_result_id.get(id(result), ())
        initial_fluid = self._initial_fluid_by_result_id.get(id(result), default_fluid_name)
        current_density = self._get_fluid_density(default_fluid_name, fluids)
        initial_density = self._get_fluid_density(initial_fluid, fluids)
        # 停泵期间，排量为零，体积推进停止；密度差会让重流体下沉、轻流体上浮。
        # v_settle = k * Δρ，其中 k 为经验沉降速度系数，Δρ 为当前流体与初始管内流体密度差。
        delta_rho = current_density - initial_density
        v_settle_m_s = self.settling_velocity_factor * delta_rho
        shutdown_duration_s = self._shutdown_duration_at(scheduled_steps, time_s)
        # 将沉降距离折算成等效体积位移，仅用于查询鞋口相态，不改写原体积追踪前缘。
        settled_volume_m3 = cumulative_at_time + v_settle_m_s * shutdown_duration_s * result.pipe_cross_section_m2
        if settled_volume_m3 < pipe_volume_m3:
            return initial_fluid
        delayed_volume_m3 = max(settled_volume_m3 - pipe_volume_m3, 0.0)
        return self._fluid_by_injected_volume(scheduled_steps, delayed_volume_m3)

    def _settled_exit_fluid_name_enhanced(
        self,
        *,
        result: CasingFlowResult,
        scheduled_steps: tuple[_ScheduledStep, ...],
        time_s: float,
        cumulative_at_time: float,
        pipe_volume_m3: float,
        default_fluid_name: str,
    ) -> str:
        """停泵期间考虑凝胶强度和屈服应力的增强沉降模型。

        在基础密度差沉降模型之上，增加两个物理修正：
        1. 凝胶强度发展：停泵期间水泥浆逐步建立凝胶结构，抑制沉降。
           采用指数增长模型: gel_factor = gelation_max_factor * (1 - exp(-t / tau_gel))
           参考 Kelessidis et al. (JPT, 2006) 的实验数据。
        2. 屈服应力效应：有屈服应力的流体（Bingham/Herschel-Bulkley）在停泵后
           快速建立结构，沉降速度进一步降低。
           参考 Maglione et al. (SPE 56995, 1999) 的流变-沉降耦合模型。

        Args:
            result: CasingFlowResult 实例
            scheduled_steps: 内部施工步骤时间窗
            time_s: 查询时刻（秒）
            cumulative_at_time: 该时刻的累计泵入体积（m3）
            pipe_volume_m3: 管内容积（m3）
            default_fluid_name: 默认流体名称（当无沉降时返回此名称）

        Returns:
            考虑增强沉降效应后的鞋口流体名称
        """
        fluids = self._fluids_by_result_id.get(id(result), ())
        initial_fluid = self._initial_fluid_by_result_id.get(id(result), default_fluid_name)
        current_fluid = next((f for f in fluids if f.name == default_fluid_name), None)
        current_density = self._get_fluid_density(default_fluid_name, fluids)
        initial_density = self._get_fluid_density(initial_fluid, fluids)

        delta_rho = current_density - initial_density
        if delta_rho <= 0.0:
            return default_fluid_name

        shutdown_duration_s = self._shutdown_duration_at(scheduled_steps, time_s)
        if shutdown_duration_s <= 0.0:
            return default_fluid_name

        v_settle_m_s = self.settling_velocity_factor * delta_rho

        # 修正 1: 凝胶强度发展
        gel_factor = self.gelation_max_factor * (1.0 - math.exp(-shutdown_duration_s / self.gelation_time_s))

        # 修正 2: 屈服应力效应
        yield_suppression = 0.0
        if current_fluid is not None and current_fluid.yield_stress_pa is not None and current_fluid.yield_stress_pa > 0.0:
            pipe_radius_m = self._effective_pipe_radius_m()
            tau_critical = delta_rho * self.g_constant * pipe_radius_m
            if tau_critical > 1e-6:
                yield_ratio = min(current_fluid.yield_stress_pa / tau_critical, 1.0)
                yield_suppression = 0.5 * yield_ratio

        total_suppression = min(gel_factor + yield_suppression, 0.99)
        v_effective = v_settle_m_s * (1.0 - total_suppression)

        settled_volume_m3 = cumulative_at_time + v_effective * shutdown_duration_s * result.pipe_cross_section_m2
        if settled_volume_m3 < pipe_volume_m3:
            return initial_fluid
        delayed_volume_m3 = max(settled_volume_m3 - pipe_volume_m3, 0.0)
        return self._fluid_by_injected_volume(scheduled_steps, delayed_volume_m3)

    @staticmethod
    def _shutdown_duration_at(scheduled_steps: tuple[_ScheduledStep, ...], time_s: float) -> float:
        """计算当前停泵段已经持续的时间。"""

        for scheduled in scheduled_steps:
            if scheduled.start_time_s <= time_s < scheduled.end_time_s - 1.0e-12:
                if scheduled.step.rate_m3_min <= 0.0:
                    return max(time_s - scheduled.start_time_s, 0.0)
                return 0.0
        last_flow_end_s = 0.0
        for scheduled in scheduled_steps:
            if scheduled.end_time_s <= time_s and scheduled.step.rate_m3_min > 0.0:
                last_flow_end_s = scheduled.end_time_s
        return max(time_s - last_flow_end_s, 0.0)

    @staticmethod
    def _get_fluid_density(fluid_name: str, fluids: tuple[FluidSpec, ...]) -> float:
        """按流体名称读取密度；未知流体用清水密度保持查询稳健。"""

        for fluid in fluids:
            if fluid.name == fluid_name:
                return fluid.density_kg_m3
        return 1000.0

    @staticmethod
    def _active_step_at(scheduled_steps: tuple[_ScheduledStep, ...], time_s: float) -> _ScheduledStep | None:
        """查找某时刻正在执行的施工步骤。"""

        for scheduled in scheduled_steps:
            if scheduled.start_time_s <= time_s < scheduled.end_time_s - 1.0e-12:
                return scheduled
        return None

    @staticmethod
    def _cumulative_volume_at(scheduled_steps: tuple[_ScheduledStep, ...], time_s: float) -> float:
        """计算某时刻累计泵入体积。"""

        cumulative_m3 = 0.0
        for scheduled in scheduled_steps:
            if time_s >= scheduled.end_time_s:
                cumulative_m3 = scheduled.cumulative_volume_end_m3
                continue
            if time_s < scheduled.start_time_s:
                return cumulative_m3
            elapsed_s = max(time_s - scheduled.start_time_s, 0.0)
            return scheduled.cumulative_volume_start_m3 + scheduled.step.rate_m3_min * elapsed_s / 60.0
        return cumulative_m3

    @staticmethod
    def _fluid_by_injected_volume(scheduled_steps: tuple[_ScheduledStep, ...], volume_m3: float) -> str:
        """按进入管柱顶部的体积坐标查找对应流体。"""

        for scheduled in scheduled_steps:
            if scheduled.cumulative_volume_start_m3 <= volume_m3 < scheduled.cumulative_volume_end_m3 - 1.0e-12:
                return scheduled.step.fluid_name
        return scheduled_steps[-1].step.fluid_name