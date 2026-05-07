"""水泥浆体积/质量守恒诊断。

本模块只做后处理验证：把环空二维结果中的水泥占据量，与施工程序中进入井筒的
水泥浆体积进行对比。这里不引入 CBL 校准、泥饼清除、温度、凝胶、湍流或质量
响应折减等工程修正，避免污染当前 cemdisp 求解主链路。
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from numpy.typing import NDArray

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep
from cemdisp.data.well_spec import DepthValuePoint, WellSpec
from cemdisp.models2d.annulus_d2dga import AnnulusSimulationResult


Array = NDArray[np.float64]
_CEMENT_ROLES = {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}
_CEMENT_NAME_MARKERS = ("水泥", "领浆", "尾浆", "中间浆")


@dataclass(frozen=True)
class CementMassBalanceDiagnostics:
    """水泥浆守恒诊断结果。

    字段均采用英文命名，便于代码调用；字段含义通过中文注释说明。
    """

    well_name: str
    scheduled_cement_volume_m3: float  # 当前模拟时长内送入入口的水泥浆体积
    result_cement_volume_m3: float  # 二维环空结果中水泥浆占据体积
    annular_volume_m3: float  # WellSpec 计算得到的物理环空体积
    volume_error_m3: float  # 结果体积 - 计划水泥浆体积
    volume_relative_error: float  # 相对计划水泥浆体积的误差
    retained_volume_fraction: float  # 结果水泥体积 / 计划水泥体积
    scheduled_cement_mass_kg: float | None = None  # 若能匹配流体密度，则给出当前模拟时长内的入口质量
    result_cement_mass_kg: float | None = None  # 若能匹配流体密度，则给出结果质量
    mass_error_kg: float | None = None
    mass_relative_error: float | None = None
    passed: bool = True  # 是否落在给定体积相对误差容差内
    tolerance: float = 0.25

    def as_dict(self) -> dict[str, float | str | bool | None]:
        """转换为适合报告或测试断言使用的扁平字典。"""
        return {
            "井名": self.well_name,
            "计划水泥浆体积_m3": self.scheduled_cement_volume_m3,
            "结果水泥浆体积_m3": self.result_cement_volume_m3,
            "物理环空体积_m3": self.annular_volume_m3,
            "水泥浆体积误差_m3": self.volume_error_m3,
            "水泥浆体积相对误差": self.volume_relative_error,
            "水泥浆保留体积分数": self.retained_volume_fraction,
            "计划水泥浆质量_kg": self.scheduled_cement_mass_kg,
            "结果水泥浆质量_kg": self.result_cement_mass_kg,
            "水泥浆质量误差_kg": self.mass_error_kg,
            "水泥浆质量相对误差": self.mass_relative_error,
            "是否通过体积守恒容差": self.passed,
            "体积相对误差容差": self.tolerance,
        }


def validate_cement_mass_balance(
    result: AnnulusSimulationResult,
    schedule: PumpingSchedule,
    well_spec: WellSpec,
    fluids: tuple[FluidSpec, ...] = (),
    *,
    tolerance: float = 0.25,
) -> CementMassBalanceDiagnostics:
    """计算水泥浆体积/质量守恒诊断。

    Args:
        result: 环空二维求解结果。
        schedule: 地面施工程序，用于统计计划水泥浆注入体积。
        well_spec: 井筒规格，用于独立计算物理环空体积。
        fluids: 可选流体规格；提供后可按密度计算质量诊断。
        tolerance: 体积相对误差通过阈值。当前用于冒烟测试和后处理提示，不反向影响求解器。

    Returns:
        ``CementMassBalanceDiagnostics``，包含体积、质量和通过性字段。
    """
    if tolerance < 0.0:
        raise ValueError("tolerance 必须为非负数")

    final_time_s = _result_final_time_s(result)
    fluid_by_name = {fluid.name: fluid for fluid in fluids}
    cement_step_names = _cement_step_names(schedule, fluid_by_name)
    scheduled_volume = _scheduled_cement_volume(schedule, cement_step_names, final_time_s)
    scheduled_mass = _scheduled_cement_mass(schedule, fluid_by_name, cement_step_names, final_time_s)

    result_volume = 2.0 * _trapez2d(result.geom["b"] * np.clip(result.cement_field, 0.0, 1.0), result.geom)
    annular_volume = _physical_annular_volume(well_spec)
    result_mass = _result_cement_mass(result, fluid_by_name)

    volume_error = result_volume - scheduled_volume
    volume_relative_error = _relative_error(volume_error, scheduled_volume)
    retained_fraction = result_volume / scheduled_volume if scheduled_volume > 0.0 else 0.0

    mass_error: float | None = None
    mass_relative_error: float | None = None
    if scheduled_mass is not None and result_mass is not None:
        mass_error = result_mass - scheduled_mass
        mass_relative_error = _relative_error(mass_error, scheduled_mass)

    return CementMassBalanceDiagnostics(
        well_name=result.well_name,
        scheduled_cement_volume_m3=scheduled_volume,
        result_cement_volume_m3=result_volume,
        annular_volume_m3=annular_volume,
        volume_error_m3=volume_error,
        volume_relative_error=volume_relative_error,
        retained_volume_fraction=retained_fraction,
        scheduled_cement_mass_kg=scheduled_mass,
        result_cement_mass_kg=result_mass,
        mass_error_kg=mass_error,
        mass_relative_error=mass_relative_error,
        passed=abs(volume_relative_error) <= tolerance,
        tolerance=tolerance,
    )


def _cement_step_names(schedule: PumpingSchedule, fluid_by_name: dict[str, FluidSpec]) -> set[str]:
    """识别施工程序中属于水泥浆的流体名称。"""
    names: set[str] = set()
    for step in schedule.steps:
        fluid = fluid_by_name.get(step.fluid_name)
        if fluid is not None and fluid.role in _CEMENT_ROLES:
            names.add(step.fluid_name)
        elif fluid is None and any(marker in step.fluid_name for marker in _CEMENT_NAME_MARKERS):
            # 未提供 FluidSpec 时保留中文名称兜底，便于轻量测试和早期数据检查。
            names.add(step.fluid_name)
    return names


def _scheduled_cement_mass(
    schedule: PumpingSchedule,
    fluid_by_name: dict[str, FluidSpec],
    cement_step_names: set[str],
    final_time_s: float,
) -> float | None:
    """按当前模拟结束时间累计已送入入口的水泥质量；密度缺失则返回 None。"""
    mass = 0.0
    elapsed_s = 0.0
    for step in schedule.steps:
        if elapsed_s >= final_time_s:
            break
        duration_s = _step_duration_s(step)
        overlap_s = max(min(final_time_s, elapsed_s + duration_s) - elapsed_s, 0.0)
        if overlap_s <= 0.0 or step.fluid_name not in cement_step_names:
            elapsed_s += duration_s
            continue
        fluid = fluid_by_name.get(step.fluid_name)
        if fluid is None:
            return None
        delivered_volume = step.volume_m3 if duration_s <= 0.0 else step.volume_m3 * overlap_s / duration_s
        mass += delivered_volume * fluid.density_kg_m3
        elapsed_s += duration_s
    return mass


def _scheduled_cement_volume(
    schedule: PumpingSchedule,
    cement_step_names: set[str],
    final_time_s: float,
) -> float:
    """按当前模拟结束时间累计已送入入口的水泥体积。"""
    delivered_volume = 0.0
    elapsed_s = 0.0
    for step in schedule.steps:
        if elapsed_s >= final_time_s:
            break
        duration_s = _step_duration_s(step)
        overlap_s = max(min(final_time_s, elapsed_s + duration_s) - elapsed_s, 0.0)
        if overlap_s > 0.0 and step.fluid_name in cement_step_names:
            delivered_volume += step.volume_m3 if duration_s <= 0.0 else step.volume_m3 * overlap_s / duration_s
        elapsed_s += duration_s
    return delivered_volume


def _result_cement_mass(result: AnnulusSimulationResult, fluid_by_name: dict[str, FluidSpec]) -> float | None:
    """根据结果中的领浆/尾浆场和对应密度估算环空内水泥质量。"""
    mass = 0.0
    has_phase = False
    for phase_field, role in ((result.lead_field, FluidRole.LEAD), (result.tail_field, FluidRole.TAIL)):
        if phase_field.size == 0:
            continue
        fluid = next((item for item in fluid_by_name.values() if item.role == role), None)
        if fluid is None:
            return None
        phase_volume = 2.0 * _trapez2d(result.geom["b"] * np.clip(phase_field, 0.0, 1.0), result.geom)
        mass += phase_volume * fluid.density_kg_m3
        has_phase = True
    return mass if has_phase else None


def _profile_to_arrays(points: tuple[DepthValuePoint, ...]) -> tuple[Array, Array]:
    """将井深剖面点转换为数组，供体积积分使用。"""
    depths = np.array([point.depth_md_m for point in points], dtype=float)
    values = np.array([point.value for point in points], dtype=float)
    return depths, values


def _physical_annular_volume(well_spec: WellSpec) -> float:
    """由井径剖面和尾管外径计算物理环空体积。"""
    cal_md, cal_hole = _profile_to_arrays(well_spec.hole_diameter_profile)
    od = np.full_like(cal_md, float(well_spec.liner_od_mm or 0.0), dtype=float)
    area = np.pi * ((cal_hole / 1000.0) ** 2 - (od / 1000.0) ** 2) / 4.0
    return float(np.trapezoid(area, x=cal_md))


def _trapez2d(arr: Array, geom: dict[str, Array]) -> float:
    """对半环空二维网格做梯形积分。"""
    return float(np.trapezoid(np.trapezoid(arr, x=geom["s"], axis=1), x=geom["y"], axis=0))


def _relative_error(error: float, reference: float) -> float:
    """计算相对误差；参考值为零时返回零，避免早期空水泥程序报错。"""
    if reference == 0.0:
        return 0.0
    return error / reference


def _result_final_time_s(result: AnnulusSimulationResult) -> float:
    """读取结果时间轴的最终时刻。"""
    if result.time_points_s:
        return float(result.time_points_s[-1])
    return 0.0


def _step_duration_s(step: PumpingScheduleStep) -> float:
    """按当前步骤体积和排量换算持续时间。"""
    if step.rate_m3_min <= 0.0:
        return 0.0
    return float(step.volume_m3 / step.rate_m3_min * 60.0)
