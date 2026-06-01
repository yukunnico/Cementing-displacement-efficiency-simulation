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
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
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
    """

    def __init__(
        self,
        *,
        dt: float = 2.0,
        enable_gravity: bool = False,
        g_constant: float = 9.81,
        settling_velocity_factor: float = 0.001,
        enable_axial_dispersion: bool = False,
        dispersion_alpha: float = 0.2,
    ) -> None:
        """初始化求解器。

        Args:
            dt: 时间步长（秒），仅用于内部时间查询的容差判断
            enable_gravity: 是否启用套管内重力修正；默认关闭以保持旧模型行为
            g_constant: 重力加速度（m/s²），用于按现场重力条件缩放简化修正项
            settling_velocity_factor: 沉降速度系数，单位为 m/s 每 kg/m³ 密度差
            enable_axial_dispersion: 是否启用管内轴向弥散；默认关闭以保持旧模型行为
            dispersion_alpha: 无量纲弥散系数，默认 0.2
        """
        if not math.isfinite(dt) or dt <= 0.0:
            raise ValueError("dt 必须为大于0的有限数值")
        if not math.isfinite(g_constant) or g_constant <= 0.0:
            raise ValueError("g_constant 必须为大于0的有限数值")
        if not math.isfinite(settling_velocity_factor) or settling_velocity_factor < 0.0:
            raise ValueError("settling_velocity_factor 必须为非负有限数值")
        if not math.isfinite(dispersion_alpha) or dispersion_alpha < 0.0:
            raise ValueError("dispersion_alpha 必须为非负有限数值")
        self.dt: float = dt
        self.enable_gravity: bool = enable_gravity
        self.g_constant: float = g_constant
        self.settling_velocity_factor: float = settling_velocity_factor
        self.enable_axial_dispersion: bool = enable_axial_dispersion
        self.dispersion_alpha: float = dispersion_alpha
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
        for scheduled in scheduled_steps:
            arrival_time_s = self._front_arrival_time(scheduled, scheduled_steps, pipe_volume_m3)
            if self.enable_gravity and arrival_time_s is not None:
                arrival_time_s = self._gravity_corrected_arrival_time(arrival_time_s, scheduled.step.fluid_name, fluids)
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
        for scheduled in scheduled_steps:
            fluid = fluid_by_name.get(scheduled.step.fluid_name)
            if fluid is None or fluid.role not in {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}:
                continue
            rear_arrival_time_s = self._rear_arrival_time(scheduled, scheduled_steps, pipe_volume_m3)
            if rear_arrival_time_s is not None:
                if self.enable_gravity:
                    rear_arrival_time_s = self._gravity_corrected_arrival_time(
                        rear_arrival_time_s,
                        scheduled.step.fluid_name,
                        fluids,
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
            fluid_name = self._settled_exit_fluid_name(
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
        """计算管内层流轴向弥散系数（Taylor-Aris 型）。

        Args:
            pipe_radius_m: 套管内半径 (m)
            fluid: 当前管内流体规格
            mean_velocity_m_s: 截面平均速度 (m/s)

        Returns:
            有效轴向弥散系数 D_eff (m²/s)
        """
        if mean_velocity_m_s < 1e-9:
            return 0.0

        from cemdisp.data.fluid_spec import RheologyModel

        if fluid.rheology_model == RheologyModel.NEWTONIAN:
            # 牛顿流体: D_eff = U²R² / (192 × D_mol)
            # D_mol ~ 1e-9 for typical fluids
            d_mol = 1.0e-9
            return (mean_velocity_m_s ** 2) * (pipe_radius_m ** 2) / (192.0 * d_mol)

        elif fluid.rheology_model == RheologyModel.POWER_LAW:
            n = fluid.power_law_n if fluid.power_law_n is not None else 1.0
            # 幂律修正: κ ≈ 48n + 144 (拟合)
            k_factor = 48.0 * n + 144.0
            d_mol = 1.0e-9
            return (mean_velocity_m_s ** 2) * (pipe_radius_m ** 2) / (k_factor * d_mol)

        elif fluid.rheology_model in (RheologyModel.BINGHAM, RheologyModel.HERSCHEL_BULKLEY):
            # 有屈服应力流体：中心塞流区抑制弥散 → 用更小的 D_eff
            return self.dispersion_alpha * mean_velocity_m_s * pipe_radius_m

        else:
            return self.dispersion_alpha * mean_velocity_m_s * pipe_radius_m

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

            # 找到前一个流体
            prev_fluid = ""
            for j in range(i - 1, -1, -1):
                if events[j].phase_fractions:
                    prev_fluid = events[j].phase_fractions[0][0]
                    break
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
        for scheduled in scheduled_steps:
            event_points.append((scheduled.start_time_s, self._event_kind_for_step(scheduled.step), None))
            front_time_s = self._front_arrival_time(scheduled, scheduled_steps, pipe_volume_m3)
            if front_time_s is not None:
                if self.enable_gravity:
                    front_time_s = self._gravity_corrected_arrival_time(front_time_s, scheduled.step.fluid_name, fluids)
                event_points.append((front_time_s, ShoeEventKind.FRONT_ARRIVAL, ((scheduled.step.fluid_name, 1.0),)))
            rear_time_s = self._rear_arrival_time(scheduled, scheduled_steps, pipe_volume_m3)
            if rear_time_s is not None:
                if self.enable_gravity:
                    rear_time_s = self._gravity_corrected_arrival_time(rear_time_s, scheduled.step.fluid_name, fluids)
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
        fluids: tuple[FluidSpec, ...],
    ) -> float:
        """按流体密度对前缘到达时间做简化重力修正。"""

        fluid_density = self._get_fluid_density(fluid_name, fluids)
        # 在套管内，向下泵注时重力沿管轴方向辅助较重流体前缘下行。
        # 这里不改变原有体积追踪，只把到达鞋口的解析时间按密度做小幅提前。
        # g_constant 用于把现场重力加速度归一到标准重力，便于后续扩展井斜/重力场修正。
        gravity_scale = self.g_constant / 9.81
        gravity_factor = self.settling_velocity_factor * (fluid_density - 1000.0) / 1000.0 * gravity_scale
        return max(arrival_time_s * (1.0 - gravity_factor), 0.0)

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