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
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep
from cemdisp.data.well_spec import WellSpec
from cemdisp.transport1d.interface_tracking import InterfaceFront
from cemdisp.transport1d.pipe_exit_state import PipeExitState


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
    ) -> None:
        """初始化求解器。

        Args:
            dt: 时间步长（秒），仅用于内部时间查询的容差判断
            enable_gravity: 是否启用套管内重力修正；默认关闭以保持旧模型行为
            g_constant: 重力加速度（m/s²），用于按现场重力条件缩放简化修正项
            settling_velocity_factor: 沉降速度系数，单位为 m/s 每 kg/m³ 密度差
        """
        if not math.isfinite(dt) or dt <= 0.0:
            raise ValueError("dt 必须为大于0的有限数值")
        if not math.isfinite(g_constant) or g_constant <= 0.0:
            raise ValueError("g_constant 必须为大于0的有限数值")
        if not math.isfinite(settling_velocity_factor) or settling_velocity_factor < 0.0:
            raise ValueError("settling_velocity_factor 必须为非负有限数值")
        self.dt: float = dt
        self.enable_gravity: bool = enable_gravity
        self.g_constant: float = g_constant
        self.settling_velocity_factor: float = settling_velocity_factor
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

        # 为每个注入步骤建立“前缘”：前缘位置由累计泵入体积推动，
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

        result = CasingFlowResult(
            fronts=tuple(fronts),
            schedule_steps=schedule.steps,
            pipe_cross_section_m2=pipe_area_m2,
            shoe_md_m=shoe_depth_m,
            pumping_end_time_s=scheduled_steps[-1].end_time_s if scheduled_steps else 0.0,
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
        active_step = self._active_step_at(scheduled_steps, time_s)
        flow_rate_m3_s = 0.0 if active_step is None else active_step.step.rate_m3_min / 60.0
        stage_name = "施工结束后保持" if active_step is None else active_step.step.step_name

        # 鞋口流出流体由“当前时间前，累计泵入体积超过管内容积的最近步骤”决定。
        cumulative_at_time = self._cumulative_volume_at(scheduled_steps, time_s)
        fluid_name = self._initial_fluid_by_result_id.get(id(result), "初始管内流体")
        if cumulative_at_time >= pipe_volume_m3:
            delayed_volume_m3 = cumulative_at_time - pipe_volume_m3
            fluid_name = self._fluid_by_injected_volume(scheduled_steps, delayed_volume_m3)

        if self.enable_gravity and flow_rate_m3_s < 1.0e-9:
            fluid_name = self._settled_exit_fluid_name(
                result=result,
                scheduled_steps=scheduled_steps,
                time_s=time_s,
                cumulative_at_time=cumulative_at_time,
                pipe_volume_m3=pipe_volume_m3,
                default_fluid_name=fluid_name,
            )

        return PipeExitState(
            time_s=time_s,
            flow_rate_m3_s=flow_rate_m3_s,
            stage_name=stage_name,
            phase_fractions=((fluid_name, 1.0),),
        )

    @staticmethod
    def _pipe_cross_section_area(well_spec: WellSpec) -> float:
        """根据套管内径计算管内截面积。"""

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
