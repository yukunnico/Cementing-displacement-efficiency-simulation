"""
呼102尾管段标准数据加载器

本模块实现呼102井139.70mm尾管段固井的标准数据加载功能。

主要功能：
- 从井径/井斜CSV文件读取剖面数据
- 构建井筒几何参数（井段范围、套管尺寸、偏心度等）
- 定义钻井液、替浆液、尾管水泥浆，以及可选冲洗液/隔离液的物性参数
- 构建现场记录施工日程（默认两步：尾浆注入+替浆液推进；提取数据实际为五步）
- 提供 legacy 环空入口边界状态提供器（仅用于旧模型对比）

物理参数说明：
- 井段范围: 6823.10m - 7735.00m
- 尾管尺寸: 139.70mm OD, 108.10mm ID (考虑壁厚后)
- 尾管水泥浆: 密度1.90g/cm³；幂律流变 n=0.766, K=1.093（提取数据）
- 替浆液: 74m³, 密度2.02g/cm³, 排量0.8m³/min（提取数据）
- 钻井液（环空初始液）: 密度2.02g/cm³, Bingham PV=80mPa·s, YP=15Pa

数据来源：
- well_spec.csv: 井号、井段范围、套管尺寸、评价窗口（100413.PDF）
- fluid_spec.csv: 流体密度、流变参数（well_spacer_summary_fixed.csv汇总表）
- pumping_schedule.csv: 五步施工程序（well_spacer_summary_fixed.csv汇总表）
- validation_data.csv: CBL合格率66.65%（100413.PDF）

关于施工程序：
- 提取数据显示实际为五步：先导浆(52m³)→驱油隔离液(42m³)→领浆(81m³)→尾浆(139m³)→替浆(74m³)
- 当前loader默认两步（尾浆+替浆）以保持向后兼容
- include_wash_spacer=True时使用五步程序参数（但不包含先导浆/领浆步骤，仅平衡液+隔离液）

legacy 边界模式选项：
- "sustained_tail": 替浆期间环空入口保持尾浆（默认）
- "volume_limited": 管内推进期间环空入口保持尾浆，排量设为0
- "tail_then_mud": 替浆期结束后环空入口切换为替浆液
"""

from __future__ import annotations

import csv
from collections.abc import Callable, Iterable
from pathlib import Path

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep
from cemdisp.data.validation_data import ValidationData
from cemdisp.data.well_spec import DepthValuePoint, EvaluationWindow, WellSpec
from cemdisp.models2d.boundary_bridge import AnnulusInletState


PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_REFERENCE_ROOT = PROJECT_ROOT / "参考文档" / "呼102"
DEFAULT_CALIPER_CSV = PROJECT_ROOT / "hu102model" / "hu102_tail_caliper_inclination.csv"

# 呼102尾管段井段参数（来源: well_spec.csv → 100413.PDF）
# 井段范围: 6823.10m - 7735.00m
# 尾管尺寸: 139.70mm OD, 108.10mm ID (壁厚15.8mm)
HU102_TOP_MD_M = 6823.10        # 井段顶部测深
HU102_BOTTOM_MD_M = 7735.00     # 井段底部测深
HU102_SHOE_MD_M = 7735.00       # 套管鞋深度
HU102_HANGER_MD_M = 6823.10     # 悬挂器位置深度
HU102_CASING_ID_MM = 219.10     # 套管内径（双层套管等效内径）
HU102_LINER_OD_MM = 139.70      # 尾管外径
HU102_LINER_WALL_THICKNESS_MM = 15.80  # 尾管壁厚
HU102_LINER_ID_MM = HU102_LINER_OD_MM - 2.0 * HU102_LINER_WALL_THICKNESS_MM  # 尾管内径

# 呼102尾管段施工参数（来源: pumping_schedule.csv → well_spacer_summary_fixed.csv汇总表）
# 提取数据显示实际为五步程序：先导浆(52m³)→驱油隔离液(42m³)→领浆(81m³)→尾浆(139m³)→替浆(74m³)
# 当前loader默认两步（尾浆+替浆）保持向后兼容；include_wash_spacer=True时使用提取的五步参数
HU102_CEMENT_VOLUME_M3 = 139.0       # 尾管尾浆体积（提取数据）
HU102_CEMENT_DENSITY_KG_M3 = 1900.0  # 尾管尾浆密度（提取数据）
HU102_LEAD_CEMENT_VOLUME_M3 = 81.0   # 尾管领浆体积（提取数据）
HU102_LEAD_CEMENT_DENSITY_KG_M3 = 1900.0  # 尾管领浆密度（提取数据）
HU102_DISPLACEMENT_VOLUME_M3 = 74.0  # 替浆体积（提取数据）
HU102_DISPLACEMENT_DENSITY_KG_M3 = 2020.0  # 替浆液密度（提取数据）
# 注意：旧loader使用35t水泥+2.10g/cm³（由10042.xlsx推算），现更新为提取的139m³+1.90g/cm³
HU102_RATE_M3_MIN = 0.8              # 泵注排量（提取数据）

# 呼102五步程序参数（来源: pumping_schedule.csv）
# 仅在include_wash_spacer=True时使用，与当前默认两步程序共存
HU102_PILOT_VOLUME_M3 = 52.0         # 先导浆体积
HU102_PILOT_DENSITY_KG_M3 = 1600.0   # 先导浆密度
HU102_SPACER_OIL_VOLUME_M3 = 42.0    # 驱油隔离液体积
HU102_SPACER_OIL_DENSITY_KG_M3 = 1850.0  # 驱油隔离液密度

# 呼102流变参数（来源: fluid_spec.csv → well_spacer_summary_fixed.csv汇总表）
# 钻井液/替浆液：Bingham模型
HU102_MUD_PV_PA_S = 0.080            # 环空初始钻井液塑性粘度
HU102_MUD_YP_PA = 15.0              # 环空初始钻井液屈服值
HU102_DISPLACEMENT_PV_PA_S = 0.080   # 替浆液塑性粘度（与钻井液一致）
HU102_DISPLACEMENT_YP_PA = 15.0      # 替浆液屈服值（与钻井液一致）
# 尾管尾浆：幂律模型（提取数据）
HU102_CEMENT_POWER_LAW_N = 0.766     # 水泥浆流性指数（提取数据）
HU102_CEMENT_CONSISTENCY_K = 1.093   # 水泥浆稠度系数（提取数据）
# 尾管领浆：幂律模型（提取数据）
HU102_LEAD_POWER_LAW_N = 0.838       # 领浆流性指数（提取数据）
HU102_LEAD_CONSISTENCY_K = 0.587     # 领浆稠度系数（提取数据）
# 驱油隔离液：幂律模型（提取数据）
HU102_SPACER_POWER_LAW_N = 0.587     # 隔离液流性指数（提取数据）
HU102_SPACER_CONSISTENCY_K = 0.880   # 隔离液稠度系数（提取数据）

# 呼102前置液参数（来源: fluid_spec.csv）
# 平衡液/冲洗液（邻井呼103代理）
HU102_WASH_DENSITY_KG_M3 = 1880.0    # 冲洗液密度（提取数据）
HU102_WASH_PV_PA_S = 0.035           # 冲洗液塑性粘度（邻井代理）
HU102_WASH_YP_PA = 8.0              # 冲洗液屈服值（邻井代理）
HU102_WASH_VOLUME_M3 = 10.0          # 冲洗液体积（邻井代理）

# 注意：呼102尾管主作业（10042.xlsx Row 26）仅记录尾浆+替浆两步，未见先导浆/领浆/隔离液
# 当前loader保留两步默认程序以向后兼容；五步参数来源于well_spacer_summary_fixed.csv汇总表
# 隔离液密度来自fluid_spec.csv（1.85g/cm³），与旧loader（2.05g/cm³）不同


def _read_profile_rows(caliper_csv_path: Path) -> tuple[tuple[float, float, float], ...]:
    with caliper_csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows: list[tuple[float, float, float]] = []
        for row in csv.DictReader(handle):
            rows.append(
                (
                    float(row["depth_md_m"]),
                    float(row["hole_diameter_mm"]),
                    float(row["inclination_deg"]),
                )
            )
    if not rows:
        raise ValueError(f"Hu102 井径/井斜 CSV 为空: {caliper_csv_path}")
    return tuple(sorted(rows))


def _standoff_value(depth_md_m: float, hole_diameter_mm: float, liner_od_mm: float) -> float:
    """计算居中度(standoff)剖面。

    基于邻井呼探1-002的扶正器数据估算：
    - 呼探1-002: 95只整体式弹扶，目的层44m间距，非目的层55m间距
    - 呼102井段较短(911.9m)，采用相似间距策略

    居中度定义：standoff = 1 - 偏心度
    - standoff = 1.0: 完全居中
    - standoff = 0.0: 完全偏心
    """
    # 基于邻井数据的居中度估算
    # 目的层段(7405-7540m): 44m间距 → standoff ≈ 0.70
    # 非目的层段: 55m间距 → standoff ≈ 0.65
    standoff = 0.65  # 默认值（非目的层）
    if 7405.00 <= depth_md_m <= 7540.00:
        standoff = 0.70  # 目的层段（油气水层）
    elif depth_md_m > 7540.00:
        standoff = 0.68  # 底部段（略好于非目的层）
    # 间隙修正因子：间隙越大，居中度越好
    clearance_mm = max(hole_diameter_mm - liner_od_mm, 5.0)
    standoff *= min(max(clearance_mm / 70.0, 0.55), 1.0)
    return min(max(standoff, 0.30), 0.85)


def _build_hole_profile(profile_rows: Iterable[tuple[float, float, float]]) -> tuple[DepthValuePoint, ...]:
    return tuple(DepthValuePoint(depth_md_m=depth, value=hole_diameter) for depth, hole_diameter, _ in profile_rows)


def _build_inclination_profile(profile_rows: Iterable[tuple[float, float, float]]) -> tuple[DepthValuePoint, ...]:
    return tuple(DepthValuePoint(depth_md_m=depth, value=inclination) for depth, _, inclination in profile_rows)


def _build_standoff_profile(profile_rows: Iterable[tuple[float, float, float]]) -> tuple[DepthValuePoint, ...]:
    return tuple(
        DepthValuePoint(
            depth_md_m=depth,
            value=_standoff_value(depth, hole_diameter, HU102_LINER_OD_MM),
        )
        for depth, hole_diameter, _ in profile_rows
    )


def load_hu102_tailpipe(
    *,
    caliper_csv_path: Path | None = None,
    reference_root: Path | None = None,
    include_wash_spacer: bool = False,  # 严格按呼102主作业实录：默认不注入邻井代理前置液/隔离液
) -> tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载呼102尾管段首版模型输入。

    Args:
        caliper_csv_path: 可选井径/井斜 CSV 路径。
        reference_root: 可选参考资料根目录。
        include_wash_spacer: 是否把 0708 邻井代理的冲洗液/隔离液步骤加入泵注程序。
            默认为 False，严格按 10042.xlsx 主作业记录仅保留尾浆+替浆两步。

    Returns:
        井筒参数、流体参数、泵注程序与验证资料路径。
    """

    resolved_reference_root = reference_root or DEFAULT_REFERENCE_ROOT
    resolved_caliper_csv_path = caliper_csv_path or DEFAULT_CALIPER_CSV

    profile_rows = _read_profile_rows(resolved_caliper_csv_path)
    well_spec = WellSpec(
        well_name="呼102",
        top_md_m=HU102_TOP_MD_M,
        bottom_md_m=HU102_BOTTOM_MD_M,
        shoe_md_m=HU102_SHOE_MD_M,
        hanger_md_m=HU102_HANGER_MD_M,
        casing_id_mm=HU102_CASING_ID_MM,
        liner_od_mm=HU102_LINER_OD_MM,
        liner_id_mm=HU102_LINER_ID_MM,
        hole_diameter_profile=_build_hole_profile(profile_rows),
        inclination_profile=_build_inclination_profile(profile_rows),
        standoff_profile=_build_standoff_profile(profile_rows),
        evaluation_windows=(
            EvaluationWindow(name="CBL评价井段", top_md_m=6840.0, bottom_md_m=7665.0, window_type="cbl"),
            EvaluationWindow(name="目标层段一", top_md_m=7405.0, bottom_md_m=7480.0, window_type="target"),
            EvaluationWindow(name="目标层段二", top_md_m=7502.0, bottom_md_m=7540.0, window_type="target"),
        ),
        reference_root=resolved_reference_root,
        notes=(
            "仅对应呼102井 139.70mm 尾管段固井，不含回接固井与其他套管段。",
            "7120–7735m 井径/井斜取自 20215.xlsx Sheet4 派生 CSV；6823.10–7119.80m 双层套管段按等效井径处理。",
        ),
    )

    # 钻井液/替浆液：Bingham模型（提取数据）
    mud_fluid = FluidSpec(
        name="钻井液",
        role=FluidRole.MUD,
        density_kg_m3=HU102_DISPLACEMENT_DENSITY_KG_M3,
        rheology_model=RheologyModel.BINGHAM,
        plastic_viscosity_pa_s=HU102_MUD_PV_PA_S,
        yield_stress_pa=HU102_MUD_YP_PA,
    )
    displacement_fluid = FluidSpec(
        name="替浆液",
        role=FluidRole.DISPLACEMENT,
        density_kg_m3=HU102_DISPLACEMENT_DENSITY_KG_M3,
        rheology_model=RheologyModel.BINGHAM,
        plastic_viscosity_pa_s=HU102_DISPLACEMENT_PV_PA_S,
        yield_stress_pa=HU102_DISPLACEMENT_YP_PA,
    )
    # 冲洗液：Bingham模型（邻井呼103代理）
    wash_fluid = FluidSpec(
        name="冲洗液",
        role=FluidRole.WASH,
        density_kg_m3=HU102_WASH_DENSITY_KG_M3,
        rheology_model=RheologyModel.BINGHAM,
        plastic_viscosity_pa_s=HU102_WASH_PV_PA_S,
        yield_stress_pa=HU102_WASH_YP_PA,
    )
    # 隔离液：幂律模型（提取数据）
    spacer_fluid = FluidSpec(
        name="隔离液",
        role=FluidRole.SPACER,
        density_kg_m3=HU102_SPACER_OIL_DENSITY_KG_M3,
        rheology_model=RheologyModel.POWER_LAW,
        power_law_n=HU102_SPACER_POWER_LAW_N,
        consistency_k=HU102_SPACER_CONSISTENCY_K,
    )
    # 尾管领浆：幂律模型（提取数据）
    lead_fluid = FluidSpec(
        name="尾管领浆",
        role=FluidRole.LEAD,
        density_kg_m3=HU102_LEAD_CEMENT_DENSITY_KG_M3,
        rheology_model=RheologyModel.POWER_LAW,
        power_law_n=HU102_LEAD_POWER_LAW_N,
        consistency_k=HU102_LEAD_CONSISTENCY_K,
    )
    # 尾管尾浆：幂律模型（提取数据）
    cement_fluid = FluidSpec(
        name="尾管水泥浆",
        role=FluidRole.TAIL,
        density_kg_m3=HU102_CEMENT_DENSITY_KG_M3,
        rheology_model=RheologyModel.POWER_LAW,
        power_law_n=HU102_CEMENT_POWER_LAW_N,
        consistency_k=HU102_CEMENT_CONSISTENCY_K,
    )

    # 组装流体列表（顺序：泥浆、替浆液、领浆、尾浆、冲洗液、隔离液）
    fluids = (mud_fluid, displacement_fluid, lead_fluid, cement_fluid, wash_fluid, spacer_fluid)

    # 前置液/隔离液步骤：仅在显式要求敏感性分析时加入。
    # 呼102主作业日报未记录该两类流体，严格现场模式默认不使用邻井代理值。
    # 数据来源：Hu102二次技套(20258.doc)、Hu103回接(20314.doc/20323.doc)、呼探1-002
    optional_front_steps = ()
    if include_wash_spacer:
        optional_front_steps = (
            PumpingScheduleStep(
                step_name="注入平衡液",
                fluid_name="冲洗液",  # 角色映射为WASH，三相模型中归入隔离液相
                volume_m3=HU102_WASH_VOLUME_M3,
                rate_m3_min=HU102_RATE_M3_MIN,
                remarks=f"平衡液/冲洗液 {HU102_WASH_VOLUME_M3}m³，密度{HU102_WASH_DENSITY_KG_M3/1000:.2f}g/cm³（呼103邻井代理）。",
            ),
            PumpingScheduleStep(
                step_name="注入驱油隔离液",
                fluid_name="隔离液",
                volume_m3=HU102_SPACER_OIL_VOLUME_M3,
                rate_m3_min=HU102_RATE_M3_MIN,
                remarks=f"驱油隔离液 {HU102_SPACER_OIL_VOLUME_M3}m³，密度{HU102_SPACER_OIL_DENSITY_KG_M3/1000:.2f}g/cm³（提取数据）。",
            ),
        )

    schedule = PumpingSchedule(
        steps=optional_front_steps + (
            PumpingScheduleStep(
                step_name="注入尾管水泥浆",
                fluid_name="尾管水泥浆",
                volume_m3=HU102_CEMENT_VOLUME_M3,
                rate_m3_min=HU102_RATE_M3_MIN,
                remarks=f"尾管尾浆 {HU102_CEMENT_VOLUME_M3}m³，密度{HU102_CEMENT_DENSITY_KG_M3/1000:.2f}g/cm³（提取数据）。",
            ),
            PumpingScheduleStep(
                step_name="替浆液推进",
                fluid_name="替浆液",
                volume_m3=HU102_DISPLACEMENT_VOLUME_M3,
                rate_m3_min=HU102_RATE_M3_MIN,
                remarks=f"替浆液 {HU102_DISPLACEMENT_VOLUME_M3}m³，密度{HU102_DISPLACEMENT_DENSITY_KG_M3/1000:.2f}g/cm³（提取数据）。",
            ),
        ),
        notes=(
            "提取数据（pumping_schedule.csv）显示实际为五步：先导浆(52m³)→驱油隔离液(42m³)→领浆(81m³)→尾浆(139m³)→替浆(74m³)。",
            "当前loader默认两步（尾浆+替浆）保持向后兼容；include_wash_spacer=True时使用邻井代理的前置液参数。",
            "五步参数来源：well_spacer_summary_fixed.csv汇总表；两步参数来源：10042.xlsx Row 26。",
            "水泥浆密度从旧版2.10g/cm³更新为1.90g/cm³（提取数据），体积从16.67m³更新为139m³。",
        ),
    )

    validation_data = ValidationData(
        cbl_summary_path=resolved_reference_root / "1004" / "10041" / "100413.PDF",
        cbl_pass_rate=0.6665,
        job_report_path=resolved_reference_root / "1004" / "10042.xlsx",
        pump_pressure_series_path=resolved_reference_root / "1004" / "100492.xlsx",
        notes=(
            "100413.PDF 给出 CBL 合格率 66.65%，评价井段 6840–7665m。",
            "10042.xlsx Row 26 提供尾管固井主作业实录：水泥浆35t、密度2.10g/cm³、替浆量74m³（现场记录）。",
            "提取数据（well_spacer_summary_fixed.csv汇总表）显示实际水泥浆密度1.90g/cm³、体积139m³（方案设计）。",
            "流变参数（n、K）已从占位值更新为提取数据值。",
        ),
    )
    return well_spec, fluids, schedule, validation_data


def build_hu102_annulus_inlet_provider(
    schedule: PumpingSchedule,
    fluids: tuple[FluidSpec, ...],
    annulus_boundary_mode: str = "sustained_tail",
) -> Callable[[float], AnnulusInletState]:
    """为 Hu102 尾管段首版模型构建环空入口边界提供器。"""

    fluid_role_by_name: dict[str, FluidRole] = {fluid.name: fluid.role for fluid in fluids}
    cement_fluid_names = {fluid.name for fluid in fluids if fluid.role in {FluidRole.LEAD, FluidRole.TAIL}}
    if not cement_fluid_names:
        raise ValueError("至少需要一个角色为 LEAD 或 TAIL 的水泥浆流体")
    default_cement_name = next(iter(sorted(cement_fluid_names)))

    def _phase_fractions_for_fluid(fluid_name: str) -> tuple[tuple[str, float], ...]:
        role = fluid_role_by_name.get(fluid_name, FluidRole.MUD)
        # 三相映射规则：领浆/尾浆归为水泥相，冲洗液/隔离液归为隔离液相，其余归为泥浆相。
        if role in {FluidRole.LEAD, FluidRole.TAIL}:
            return (("cement", 1.0),)
        if role in {FluidRole.WASH, FluidRole.SPACER}:
            return (("spacer", 1.0),)
        return (("mud", 1.0),)

    steps = schedule.steps
    if len(steps) < 2:
        raise ValueError("Hu102 首版边界提供器至少需要水泥浆步骤和替浆步骤")
    # 多流体程序下替浆步骤不再固定为第 2 步，因此按流体角色定位替浆阶段。
    displacement_step_index = next(
        (
            index
            for index, step in enumerate(steps)
            if fluid_role_by_name[step.fluid_name] == FluidRole.DISPLACEMENT
        ),
        None,
    )
    if displacement_step_index is None:
        raise ValueError("Hu102 边界提供器需要至少一个替浆液步骤")
    displacement_step = steps[displacement_step_index]
    displacement_rate_m3_s = displacement_step.rate_m3_min / 60.0

    def _provider(time_s: float) -> AnnulusInletState:
        elapsed_s = 0.0
        for index, step in enumerate(steps):
            duration_s = 0.0 if step.rate_m3_min <= 0.0 else step.volume_m3 / step.rate_m3_min * 60.0
            if time_s < elapsed_s + duration_s - 1e-12:
                if index == displacement_step_index:
                    if annulus_boundary_mode == "sustained_tail":
                        return AnnulusInletState(
                            time_s=time_s,
                            flow_rate_m3_s=displacement_rate_m3_s,
                            stage_name="替浆推进（尾浆等效入环空）",
                            phase_fractions=_phase_fractions_for_fluid(default_cement_name),
                        )
                    if annulus_boundary_mode == "volume_limited":
                        return AnnulusInletState(
                            time_s=time_s,
                            flow_rate_m3_s=0.0,
                            stage_name="替浆期管内推进（环空入口保持尾浆）",
                            phase_fractions=_phase_fractions_for_fluid(default_cement_name),
                        )
                    if annulus_boundary_mode == "tail_then_mud":
                        return AnnulusInletState(
                            time_s=time_s,
                            flow_rate_m3_s=displacement_rate_m3_s,
                            stage_name="替浆钻井液入环空",
                            phase_fractions=_phase_fractions_for_fluid(step.fluid_name),
                        )
                    raise ValueError(f"Unsupported annulus boundary mode: {annulus_boundary_mode}")

                return AnnulusInletState(
                    time_s=time_s,
                    flow_rate_m3_s=step.rate_m3_min / 60.0,
                    stage_name=step.step_name,
                    phase_fractions=_phase_fractions_for_fluid(step.fluid_name),
                )
            elapsed_s += duration_s

        if annulus_boundary_mode in {"sustained_tail", "volume_limited"}:
            return AnnulusInletState(
                time_s=time_s,
                flow_rate_m3_s=0.0,
                stage_name="替浆结束后保持",
                phase_fractions=_phase_fractions_for_fluid(default_cement_name),
            )
        return AnnulusInletState(
            time_s=time_s,
            flow_rate_m3_s=0.0,
            stage_name="替浆结束后保持",
            phase_fractions=(("mud", 1.0),),
        )

    return _provider
