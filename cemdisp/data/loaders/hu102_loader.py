"""
呼102尾管段标准数据加载器

本模块实现呼102井139.70mm尾管段固井的标准数据加载功能。

主要功能：
- 从井径/井斜CSV文件读取剖面数据
- 构建井筒几何参数（井段范围、套管尺寸、偏心度等）
- 定义钻井液、替浆液、尾管水泥浆的物性参数
- 构建现场记录施工日程（尾浆注入+替浆液推进两步）
- 提供环空入口边界状态提供器（支持多种边界模式）

物理参数说明：
- 井段范围: 6823.10m - 7735.00m
- 尾管尺寸: 139.70mm OD, 108.10mm ID (考虑壁厚后)
- 水泥浆: 35t, 密度2.10g/cm³, 幂律流变 n=0.722, K=0.684
- 替浆液: 74m³, 密度2.02g/cm³, 排量1.30m³/min
- 钻井液（环空初始液）: 密度2.02g/cm³, Bingham PV=80mPa·s, YP=15Pa

现场记录来源（10042.xlsx Row 26, 2022-11-22）：
- 注水泥35.00t, 水泥浆平均密度2.10g/cm³, 替浆液密度2.02g/cm³, 井液74.00m³
- 泵注时间：2022-11-21 17:00–21:00（4小时）
- 现场记录中无冲洗液/隔离液/领浆的注入量（方案A：按现场记录）

可选补充流体（0708邻井代理，暂不强制注入）：
- 冲洗液(WASH): ρ=1880, PV=0.025, YP=1.5 (呼103邻井)
- 隔离液(SPACER): ρ=1850, PV=0.035, YP=8 (呼103邻井)

边界模式选项：
- "sustained_tail": 替浆期间环空入口保持尾浆（默认）
- "volume_limited": 管内推进期间环空入口保持尾浆，排量设为0
- "tail_then_mud": 替浆期结束后环空入口切换为替浆液
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Callable, Dict, Iterable, Tuple

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep
from cemdisp.data.validation_data import ValidationData
from cemdisp.data.well_spec import DepthValuePoint, EvaluationWindow, WellSpec
from cemdisp.models2d.boundary_bridge import AnnulusInletState


PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_REFERENCE_ROOT = PROJECT_ROOT / "参考文档" / "呼102"
DEFAULT_CALIPER_CSV = PROJECT_ROOT / "hu102model" / "hu102_tail_caliper_inclination.csv"

# 呼102尾管段井段参数（单位：m, mm）
HU102_TOP_MD_M = 6823.10        # 井段顶部测深
HU102_BOTTOM_MD_M = 7735.00     # 井段底部测深
HU102_SHOE_MD_M = 7735.00       # 套管鞋深度
HU102_HANGER_MD_M = 6823.10     # 悬挂器位置深度
HU102_CASING_ID_MM = 219.10     # 套管内径
HU102_LINER_OD_MM = 139.70      # 尾管外径
HU102_LINER_WALL_THICKNESS_MM = 15.80  # 尾管壁厚
HU102_LINER_ID_MM = HU102_LINER_OD_MM - 2.0 * HU102_LINER_WALL_THICKNESS_MM  # 尾管内径

# 呼102施工参数
HU102_CEMENT_MASS_T = 35.0              # 水泥浆质量
HU102_CEMENT_DENSITY_KG_M3 = 2100.0     # 水泥浆密度
HU102_DISPLACEMENT_VOLUME_M3 = 74.0    # 替浆体积
HU102_DISPLACEMENT_DENSITY_KG_M3 = 2020.0  # 钻井液密度（替浆用）
HU102_RATE_M3_MIN = 1.30                # 泵注排量

# 呼102流变参数 — 钻井液/替浆液/水泥浆
HU102_MUD_PV_PA_S = 0.080            # 环空初始钻井液塑性粘度（文献暂定）
HU102_MUD_YP_PA = 15.0               # 环空初始钻井液屈服值（文献暂定）
HU102_DISPLACEMENT_PV_PA_S = 0.080   # 替浆液塑性粘度（与钻井液一致）
HU102_DISPLACEMENT_YP_PA = 15.0      # 替浆液屈服值（与钻井液一致）
HU102_CEMENT_POWER_LAW_N = 0.722     # 水泥浆流性指数（幂律，legacy模型沿用）
HU102_CEMENT_CONSISTENCY_K = 0.684   # 水泥浆稠度系数（幂律，legacy模型沿用）

# 呼102可选补充流体参数 — 0708邻井代理（暂不强制注入）
HU102_WASH_DENSITY_KG_M3 = 1880.0    # 冲洗液密度（呼103邻井回接作业）
HU102_WASH_PV_PA_S = 0.025           # 冲洗液塑性粘度（呼103邻井）
HU102_WASH_YP_PA = 1.5               # 冲洗液屈服值（呼103邻井）
HU102_SPACER_DENSITY_KG_M3 = 1850.0  # 驱油隔离液密度（Hu102同井其他作业/呼103邻井）
HU102_SPACER_PV_PA_S = 0.035         # 驱油隔离液塑性粘度（呼103邻井）
HU102_SPACER_YP_PA = 8.0             # 驱油隔离液屈服值（呼103邻井）


def _read_profile_rows(caliper_csv_path: Path) -> Tuple[Tuple[float, float, float], ...]:
    with caliper_csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = []
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
    standoff = 0.58
    if 6840.00 <= depth_md_m < 7119.80:
        standoff = 0.52
    elif 7119.80 <= depth_md_m < 7405.00:
        standoff = 0.62
    elif 7405.00 <= depth_md_m <= 7540.00:
        standoff = 0.56
    elif depth_md_m > 7540.00:
        standoff = 0.60
    clearance_mm = max(hole_diameter_mm - liner_od_mm, 5.0)
    standoff *= min(max(clearance_mm / 70.0, 0.55), 1.0)
    return min(max(standoff, 0.30), 0.82)


def _build_hole_profile(profile_rows: Iterable[Tuple[float, float, float]]) -> Tuple[DepthValuePoint, ...]:
    return tuple(DepthValuePoint(depth_md_m=depth, value=hole_diameter) for depth, hole_diameter, _ in profile_rows)


def _build_inclination_profile(profile_rows: Iterable[Tuple[float, float, float]]) -> Tuple[DepthValuePoint, ...]:
    return tuple(DepthValuePoint(depth_md_m=depth, value=inclination) for depth, _, inclination in profile_rows)


def _build_standoff_profile(profile_rows: Iterable[Tuple[float, float, float]]) -> Tuple[DepthValuePoint, ...]:
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
) -> Tuple[WellSpec, Tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载呼102尾管段首版模型输入。"""

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

    fluids = (
        FluidSpec(
            name="钻井液",
            role=FluidRole.MUD,
            density_kg_m3=HU102_DISPLACEMENT_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU102_MUD_PV_PA_S,
            yield_stress_pa=HU102_MUD_YP_PA,
        ),
        FluidSpec(
            name="替浆液",
            role=FluidRole.DISPLACEMENT,
            density_kg_m3=HU102_DISPLACEMENT_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU102_DISPLACEMENT_PV_PA_S,
            yield_stress_pa=HU102_DISPLACEMENT_YP_PA,
        ),
        FluidSpec(
            name="尾管水泥浆",
            role=FluidRole.TAIL,
            density_kg_m3=HU102_CEMENT_DENSITY_KG_M3,
            rheology_model=RheologyModel.POWER_LAW,
            power_law_n=HU102_CEMENT_POWER_LAW_N,
            consistency_k=HU102_CEMENT_CONSISTENCY_K,
        ),
    )

    schedule = PumpingSchedule(
        steps=(
            PumpingScheduleStep(
                step_name="注入尾管水泥浆",
                fluid_name="尾管水泥浆",
                volume_m3=HU102_CEMENT_MASS_T / (HU102_CEMENT_DENSITY_KG_M3 / 1000.0),
                rate_m3_min=HU102_RATE_M3_MIN,
                remarks="基于 35t 与 2.10g/cm3 换算 ≈ 16.67m3。",
            ),
            PumpingScheduleStep(
                step_name="替浆液推进",
                fluid_name="替浆液",
                volume_m3=HU102_DISPLACEMENT_VOLUME_M3,
                rate_m3_min=HU102_RATE_M3_MIN,
                remarks="主作业直接记录替浆量 74m3，替浆液密度 2.02g/cm3。",
            ),
        ),
        notes=(
            "按现场记录（10042.xlsx Row 26）：仅尾浆+替浆液两步。",
            "前置液/隔离液/领浆程序缺少主作业直接证据，暂不强制注入（方案A）。",
            "可选补充：WASH(ρ=1.88) 和 SPACER(ρ=1.85) 参数已定义在模块常量中，供后续参数化使用。",
        ),
    )

    validation_data = ValidationData(
        cbl_summary_path=resolved_reference_root / "1004" / "10041" / "100413.PDF",
        job_report_path=resolved_reference_root / "1004" / "10042.xlsx",
        pump_pressure_series_path=resolved_reference_root / "1004" / "100492.xlsx",
        notes=(
            "100413.PDF 给出 CBL 合格率 66.65%，评价井段 6840–7665m。",
            "10042.xlsx 提供尾管固井主作业水泥浆 35t、平均密度 2.10g/cm3、替浆量 74m3。",
            "钻井液与水泥浆流变参数仍为首版暂定值，需后续继续用 0708 或文献补强。",
        ),
    )
    return well_spec, fluids, schedule, validation_data


def build_hu102_annulus_inlet_provider(
    schedule: PumpingSchedule,
    fluids: Tuple[FluidSpec, ...],
    annulus_boundary_mode: str = "sustained_tail",
) -> Callable[[float], AnnulusInletState]:
    """为 Hu102 尾管段首版模型构建环空入口边界提供器。"""

    fluid_role_by_name: Dict[str, FluidRole] = {fluid.name: fluid.role for fluid in fluids}
    cement_fluid_names = {fluid.name for fluid in fluids if fluid.role in {FluidRole.LEAD, FluidRole.TAIL}}
    if not cement_fluid_names:
        raise ValueError("至少需要一个角色为 LEAD 或 TAIL 的水泥浆流体")
    default_cement_name = next(iter(sorted(cement_fluid_names)))

    def _phase_fractions_for_fluid(fluid_name: str) -> Tuple[Tuple[str, float], ...]:
        role = fluid_role_by_name[fluid_name]
        if role in {FluidRole.LEAD, FluidRole.TAIL}:
            return (("cement", 1.0),)
        return (("mud", 1.0),)

    steps = schedule.steps
    if len(steps) < 2:
        raise ValueError("Hu102 首版边界提供器至少需要水泥浆步骤和替浆步骤")
    displacement_step = steps[1]
    displacement_rate_m3_s = displacement_step.rate_m3_min / 60.0

    def _provider(time_s: float) -> AnnulusInletState:
        elapsed_s = 0.0
        for index, step in enumerate(steps):
            duration_s = 0.0 if step.rate_m3_min <= 0.0 else step.volume_m3 / step.rate_m3_min * 60.0
            if time_s < elapsed_s + duration_s - 1e-12:
                if index == 1:
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
