"""
呼1-003井（HT1-003）168.3+139.7mm双径尾管段标准数据加载器

本模块把呼1-003井（HT1-003）现场数据包中的 168.3+139.7mm 双径尾管段
资料整理为 cemdisp 标准输入结构。

当前求解目标为下部 139.7mm 尾管段，因此上部 168.3mm 重叠井段用等效井眼
直径保面积近似；鞋口滞后体积按 HT1-003 明确内径分段累加。
"""

from __future__ import annotations

import math
from pathlib import Path

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep
from cemdisp.data.validation_data import ValidationData
from cemdisp.data.well_spec import DepthValuePoint, EvaluationWindow, WellSpec


PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_REFERENCE_ROOT = PROJECT_ROOT / "参考文档" / "呼1-003"

# 呼1-003井段几何参数，来源：参考文档/呼1-003/提取数据/呼1-003井_固井顶替模型数据包.json
HT1_003_WELL_NAME = "呼1-003井（HT1-003）"
HT1_003_DRILLED_DEPTH_MD_M = 7586.0  # 实际完钻井深/尾管鞋深度
HT1_003_HANGER_MD_M = 5277.754  # 尾管悬挂器位置
HT1_003_TOP_MD_M = HT1_003_HANGER_MD_M  # 模型剖面从悬挂器开始，兼容双径尾管等效处理
HT1_003_UPPER_SECTION_BOTTOM_MD_M = 7059.016  # 168.3mm 上段尾管底界/变径位置
HT1_003_BOTTOM_MD_M = HT1_003_DRILLED_DEPTH_MD_M  # 下段 139.7mm 尾管鞋
HT1_003_SHOE_MD_M = HT1_003_DRILLED_DEPTH_MD_M
HT1_003_CASING_ID_MM = 273.1  # 技术套管内径
HT1_003_UPPER_HOLE_NOMINAL_DIAMETER_MM = 241.3  # 上段井眼名义尺寸
HT1_003_LOWER_HOLE_DIAMETER_MM = 215.9  # 下段井眼名义尺寸
HT1_003_BIT_DIAMETER_LOWER_MM = 215.9  # 下段钻头尺寸
HT1_003_UPPER_LINER_OD_MM = 168.3  # 上段尾管外径
HT1_003_LOWER_LINER_OD_MM = 139.7  # 下段尾管外径
HT1_003_LOWER_LINER_WALL_THICKNESS_MM = 15.88  # 139.7mm 管壁厚
HT1_003_LOWER_LINER_ID_MM = 107.94  # 139.7 - 2*15.88
HT1_003_UPPER_LINER_WALL_THICKNESS_MM = 14.7  # 168.3mm 尾管壁厚
HT1_003_UPPER_LINER_ID_MM = HT1_003_UPPER_LINER_OD_MM - 2.0 * HT1_003_UPPER_LINER_WALL_THICKNESS_MM
HT1_003_LOWER_CENTRALIZER_COUNT = 78  # 139.7mm 下段整体式扶正器数量
HT1_003_UPPER_CENTRALIZER_COUNT = 21  # 168.3mm 上段整体式扶正器数量
HT1_003_CENTRALIZER_COUNT = HT1_003_LOWER_CENTRALIZER_COUNT + HT1_003_UPPER_CENTRALIZER_COUNT


def _equivalent_hole_diameter_mm(actual_hole_mm: float, actual_od_mm: float, reference_od_mm: float) -> float:
    """在固定 reference_od 条件下，构造与原始环空面积等价的井眼直径。"""
    area_term = actual_hole_mm**2 - actual_od_mm**2 + reference_od_mm**2
    return math.sqrt(max(area_term, reference_od_mm**2))


HT1_003_UPPER_HOLE_DIAMETER_MM = _equivalent_hole_diameter_mm(
    actual_hole_mm=HT1_003_UPPER_HOLE_NOMINAL_DIAMETER_MM,
    actual_od_mm=HT1_003_UPPER_LINER_OD_MM,
    reference_od_mm=HT1_003_LOWER_LINER_OD_MM,
)


def _pipe_volume_m3(length_m: float, inner_diameter_mm: float) -> float:
    """按内径和长度估算管内容积，用于鞋口滞后体积计算。"""
    radius_m = inner_diameter_mm / 1000.0 / 2.0
    return math.pi * radius_m**2 * length_m


# 鞋口滞后体积估算
HT1_003_SURFACE_TO_HANGER_EFFECTIVE_ID_MM = math.sqrt(4.0 * 52.0 / (math.pi * 7868.0)) * 1000.0
HT1_003_SHOE_LAG_VOLUME_M3 = (
    _pipe_volume_m3(HT1_003_HANGER_MD_M, HT1_003_SURFACE_TO_HANGER_EFFECTIVE_ID_MM)
    + _pipe_volume_m3(
        HT1_003_UPPER_SECTION_BOTTOM_MD_M - HT1_003_HANGER_MD_M,
        HT1_003_UPPER_LINER_ID_MM,
    )
    + _pipe_volume_m3(
        HT1_003_BOTTOM_MD_M - HT1_003_UPPER_SECTION_BOTTOM_MD_M,
        HT1_003_LOWER_LINER_ID_MM,
    )
)
HT1_003_LINER_ID_MM = HT1_003_LOWER_LINER_ID_MM

# 呼1-003现场流体参数
HT1_003_MUD_DENSITY_KG_M3 = 1950.0
HT1_003_BALANCE_DENSITY_KG_M3 = 1750.0
HT1_003_SPACER_DENSITY_KG_M3 = 2000.0
HT1_003_LEAD_DENSITY_KG_M3 = 2050.0
HT1_003_INTERMEDIATE_DENSITY_KG_M3 = 1950.0
HT1_003_TAIL_DENSITY_KG_M3 = 1950.0
HT1_003_DISPLACEMENT_DENSITY_KG_M3 = 1950.0
HT1_003_PLUG_DENSITY_KG_M3 = 1970.0
HT1_003_BASE_FLUID_DENSITY_KG_M3 = 1020.0
HT1_003_WELL_MUD_DENSITY_KG_M3 = 1950.0

# 钻井液流变参数（施工设计提供）
HT1_003_MUD_PV_PA_S = 0.051  # 51 mPa·s
HT1_003_MUD_YP_PA = 10.0
HT1_003_MUD_POWER_LAW_N = 0.82
HT1_003_MUD_CONSISTENCY_K = 0.21

# 其他流体简化参数（untitled.m提供）
HT1_003_BALANCE_PV_PA_S = 0.055
HT1_003_BALANCE_YP_PA = 9.2
HT1_003_SPACER_PV_PA_S = 0.030
HT1_003_SPACER_YP_PA = 9.0
HT1_003_LEAD_PV_PA_S = 0.050
HT1_003_LEAD_YP_PA = 11.0
HT1_003_INTERMEDIATE_PV_PA_S = 0.180
HT1_003_INTERMEDIATE_YP_PA = 14.0
HT1_003_TAIL_PV_PA_S = 0.180
HT1_003_TAIL_YP_PA = 14.0
HT1_003_PLUG_PV_PA_S = 0.030
HT1_003_PLUG_YP_PA = 8.0
HT1_003_DISPLACEMENT_PV_PA_S = 0.030
HT1_003_DISPLACEMENT_YP_PA = 8.5
HT1_003_BASE_FLUID_PV_PA_S = 0.030
HT1_003_BASE_FLUID_YP_PA = 8.0

# 呼1-003现场施工程序参数（按施工设计文档）
HT1_003_BALANCE_VOLUME_M3 = 28.0  # 先导浆
HT1_003_SPACER_VOLUME_M3 = 26.0  # 驱油隔离液
HT1_003_LEAD_VOLUME_M3 = 26.0  # 领浆（设计值）
HT1_003_INTERMEDIATE_VOLUME_M3 = 26.0  # 中间浆（设计值）
HT1_003_TAIL_VOLUME_M3 = 15.0  # 尾浆（设计值）
HT1_003_PLUG_VOLUME_M3 = 2.0  # 压塞液
HT1_003_FAST_MUD_VOLUME_M3 = 25.0  # 替浆首段
HT1_003_BUFFER_VOLUME_M3 = 14.0  # 保护液
HT1_003_BASE_FLUID_VOLUME_M3 = 1.0  # 基液
HT1_003_WELL_MUD_FAST_VOLUME_M3 = 10.0  # 替浆分段1
HT1_003_WELL_MUD_MID1_VOLUME_M3 = 5.0  # 替浆分段2
HT1_003_WELL_MUD_MID2_VOLUME_M3 = 15.0  # 替浆分段3
HT1_003_WELL_MUD_MID3_VOLUME_M3 = 10.0  # 替浆分段4
HT1_003_WELL_MUD_SLOW_VOLUME_M3 = 11.53  # 替浆末段（碰压）

# 排量参数
HT1_003_BALANCE_RATE_M3_MIN = 1.2
HT1_003_SPACER_RATE_M3_MIN = 1.2
HT1_003_CEMENT_RATE_M3_MIN = 1.2
HT1_003_PLUG_RATE_M3_MIN = 1.6
HT1_003_FAST_MUD_RATE_M3_MIN = 1.6
HT1_003_BUFFER_RATE_M3_MIN = 1.4
HT1_003_BASE_FLUID_RATE_M3_MIN = 1.4
HT1_003_WELL_MUD_FAST_RATE_M3_MIN = 1.4
HT1_003_WELL_MUD_MID1_RATE_M3_MIN = 1.2
HT1_003_WELL_MUD_MID2_RATE_M3_MIN = 1.0
HT1_003_WELL_MUD_MID3_RATE_M3_MIN = 0.8
HT1_003_WELL_MUD_SLOW_RATE_M3_MIN = 0.7


def _depth_points(values: tuple[tuple[float, float], ...]) -> tuple[DepthValuePoint, ...]:
    """把测深-数值元组转换为 WellSpec 使用的剖面点。"""
    return tuple(DepthValuePoint(depth_md_m=depth, value=value) for depth, value in values)


def load_ht1_003_tailpipe(
    *,
    reference_root: Path | None = None,
) -> tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载呼1-003井（HT1-003）168.3+139.7mm双径尾管段标准模型输入。"""

    resolved_reference_root = reference_root or DEFAULT_REFERENCE_ROOT
    well_spec = WellSpec(
        well_name=HT1_003_WELL_NAME,
        top_md_m=HT1_003_TOP_MD_M,
        bottom_md_m=HT1_003_BOTTOM_MD_M,
        shoe_md_m=HT1_003_SHOE_MD_M,
        hanger_md_m=HT1_003_HANGER_MD_M,
        casing_id_mm=HT1_003_CASING_ID_MM,
        liner_od_mm=HT1_003_LOWER_LINER_OD_MM,
        liner_id_mm=HT1_003_LINER_ID_MM,
        hole_diameter_profile=_depth_points(
            (
                (HT1_003_TOP_MD_M, HT1_003_UPPER_HOLE_DIAMETER_MM),
                (6000.0, HT1_003_UPPER_HOLE_DIAMETER_MM),
                (6600.0, HT1_003_UPPER_HOLE_DIAMETER_MM),
                (HT1_003_UPPER_SECTION_BOTTOM_MD_M - 1.0, HT1_003_UPPER_HOLE_DIAMETER_MM),
                (
                    HT1_003_UPPER_SECTION_BOTTOM_MD_M,
                    0.5 * (HT1_003_UPPER_HOLE_DIAMETER_MM + HT1_003_LOWER_HOLE_DIAMETER_MM),
                ),
                (7400.0, HT1_003_LOWER_HOLE_DIAMETER_MM),
                (7500.0, HT1_003_LOWER_HOLE_DIAMETER_MM),
                (HT1_003_BOTTOM_MD_M, HT1_003_LOWER_HOLE_DIAMETER_MM),
            )
        ),
        inclination_profile=_depth_points(
            (
                (HT1_003_TOP_MD_M, 2.0),
                (6000.0, 3.0),
                (6600.0, 3.5),
                (HT1_003_UPPER_SECTION_BOTTOM_MD_M - 1.0, 3.8),
                (HT1_003_UPPER_SECTION_BOTTOM_MD_M, 3.9),
                (7400.0, 4.0),
                (7500.0, 4.1),
                (HT1_003_BOTTOM_MD_M, 4.077923),
            )
        ),
        standoff_profile=_depth_points(
            (
                (HT1_003_TOP_MD_M, 0.64),
                (6000.0, 0.62),
                (6600.0, 0.60),
                (HT1_003_UPPER_SECTION_BOTTOM_MD_M - 1.0, 0.58),
                (HT1_003_UPPER_SECTION_BOTTOM_MD_M, 0.68),
                (7400.0, 0.70),
                (7500.0, 0.72),
                (HT1_003_BOTTOM_MD_M, 0.70),
            )
        ),
        evaluation_windows=(
            EvaluationWindow(name="CBL评价井段", top_md_m=5700.0, bottom_md_m=7586.0, window_type="cbl"),
            EvaluationWindow(name="目标层段", top_md_m=7400.0, bottom_md_m=7500.0, window_type="target"),
        ),
        reference_root=resolved_reference_root,
        notes=(
            "呼1-003井（HT1-003）为168.3mm+139.7mm双径复合尾管控压固井，井深7586m。",
            "上段 168.3mm 尾管按等效井眼直径保面积处理，下段 139.7mm 尾管为模型目标井段。",
            "下段井径采用 215.9mm；上段尾管内径采用明确值 138.9mm，下段尾管内径采用明确值 107.94mm。",
            f"鞋口滞后体积估算为 {HT1_003_SHOE_LAG_VOLUME_M3:.2f}m³，WellSpec.liner_id_mm 使用下段内径 {HT1_003_LINER_ID_MM:.2f}mm。",
            "扶正器布置：168.3mm段21个(5286-7057m)，139.7mm段78个(7057-7586m)，间距2根/只。",
            "钻井液流变参数：n=0.82, K=0.21 Pa·s^n (Power Law)，PV=51mPa·s, YP=10Pa。",
            "稠化时间(152℃/145.1MPa)：领浆430-490min/中间浆270-330min/尾浆200-260min。",
            "API失水量均<50mL；24h抗压强度：领浆>3.5MPa/中间浆>14MPa/尾浆>14MPa。",
        ),
    )

    # HT1-003 流体清单：钻井液使用幂律模型（有实验n/K值），其他使用Bingham代理
    fluids = (
        FluidSpec("钻井液", FluidRole.MUD, HT1_003_MUD_DENSITY_KG_M3, RheologyModel.POWER_LAW,
                  power_law_n=HT1_003_MUD_POWER_LAW_N, consistency_k=HT1_003_MUD_CONSISTENCY_K),
        FluidSpec("替浆液", FluidRole.DISPLACEMENT, HT1_003_DISPLACEMENT_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_DISPLACEMENT_PV_PA_S, HT1_003_DISPLACEMENT_YP_PA),
        FluidSpec("平衡液", FluidRole.WASH, HT1_003_BALANCE_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_BALANCE_PV_PA_S, HT1_003_BALANCE_YP_PA),
        FluidSpec("隔离液", FluidRole.SPACER, HT1_003_SPACER_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_SPACER_PV_PA_S, HT1_003_SPACER_YP_PA),
        FluidSpec("领浆", FluidRole.LEAD, HT1_003_LEAD_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_LEAD_PV_PA_S, HT1_003_LEAD_YP_PA),
        FluidSpec("中间浆", FluidRole.INTERMEDIATE, HT1_003_INTERMEDIATE_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_INTERMEDIATE_PV_PA_S, HT1_003_INTERMEDIATE_YP_PA),
        FluidSpec("尾浆", FluidRole.TAIL, HT1_003_TAIL_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_TAIL_PV_PA_S, HT1_003_TAIL_YP_PA),
        FluidSpec("压塞液", FluidRole.OTHER, HT1_003_PLUG_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_PLUG_PV_PA_S, HT1_003_PLUG_YP_PA),
        FluidSpec("替钻井液", FluidRole.DISPLACEMENT, HT1_003_DISPLACEMENT_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_DISPLACEMENT_PV_PA_S, HT1_003_DISPLACEMENT_YP_PA),
        FluidSpec("基液", FluidRole.DISPLACEMENT, HT1_003_BASE_FLUID_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_BASE_FLUID_PV_PA_S, HT1_003_BASE_FLUID_YP_PA),
        FluidSpec("井浆", FluidRole.DISPLACEMENT, HT1_003_WELL_MUD_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_DISPLACEMENT_PV_PA_S, HT1_003_DISPLACEMENT_YP_PA),
    )

    # HT1-003 地面施工程序（按施工设计文档）
    schedule = PumpingSchedule(
        steps=(
            PumpingScheduleStep("注入先导浆", "平衡液", HT1_003_BALANCE_VOLUME_M3, HT1_003_BALANCE_RATE_M3_MIN,
                                remarks="先导浆 28m³@1.2m³/min，密度1.75g/cm³。"),
            PumpingScheduleStep("注入驱油隔离液", "隔离液", HT1_003_SPACER_VOLUME_M3, HT1_003_SPACER_RATE_M3_MIN,
                                remarks="驱油隔离液 26m³@1.2m³/min，密度2.00g/cm³。"),
            PumpingScheduleStep("注入领浆", "领浆", HT1_003_LEAD_VOLUME_M3, HT1_003_CEMENT_RATE_M3_MIN,
                                remarks="领浆 26m³@1.2m³/min，密度2.05g/cm³。"),
            PumpingScheduleStep("注入中间浆", "中间浆", HT1_003_INTERMEDIATE_VOLUME_M3, HT1_003_CEMENT_RATE_M3_MIN,
                                remarks="中间浆 26m³@1.2m³/min，密度1.95g/cm³。"),
            PumpingScheduleStep("注入尾浆", "尾浆", HT1_003_TAIL_VOLUME_M3, HT1_003_CEMENT_RATE_M3_MIN,
                                remarks="尾浆 15m³@1.2m³/min，密度1.95g/cm³。"),
            PumpingScheduleStep("注入压塞液（管内）", "压塞液", HT1_003_PLUG_VOLUME_M3, HT1_003_PLUG_RATE_M3_MIN,
                                remarks="压塞液 2m³@1.6m³/min，仅作为管内占位。"),
            PumpingScheduleStep("替钻井液(快)", "替钻井液", HT1_003_FAST_MUD_VOLUME_M3, HT1_003_FAST_MUD_RATE_M3_MIN,
                                remarks="替钻井液快替 25m³@1.6m³/min。"),
            PumpingScheduleStep("替保护液", "替浆液", HT1_003_BUFFER_VOLUME_M3, HT1_003_BUFFER_RATE_M3_MIN,
                                remarks="保护液 14m³@1.4m³/min，密度1.97g/cm³。"),
            PumpingScheduleStep("替基液", "基液", HT1_003_BASE_FLUID_VOLUME_M3, HT1_003_BASE_FLUID_RATE_M3_MIN,
                                remarks="基液 1m³@1.4m³/min，密度1.02g/cm³。"),
            PumpingScheduleStep("井浆快替1", "井浆", HT1_003_WELL_MUD_FAST_VOLUME_M3, HT1_003_WELL_MUD_FAST_RATE_M3_MIN,
                                remarks="井浆快替1 10m³@1.4m³/min。"),
            PumpingScheduleStep("井浆中替1", "井浆", HT1_003_WELL_MUD_MID1_VOLUME_M3, HT1_003_WELL_MUD_MID1_RATE_M3_MIN,
                                remarks="井浆中替1 5m³@1.2m³/min。"),
            PumpingScheduleStep("井浆中替2", "井浆", HT1_003_WELL_MUD_MID2_VOLUME_M3, HT1_003_WELL_MUD_MID2_RATE_M3_MIN,
                                remarks="井浆中替2 15m³@1.0m³/min。"),
            PumpingScheduleStep("井浆中替3", "井浆", HT1_003_WELL_MUD_MID3_VOLUME_M3, HT1_003_WELL_MUD_MID3_RATE_M3_MIN,
                                remarks="井浆中替3 10m³@0.8m³/min。"),
            PumpingScheduleStep("井浆慢替（碰压）", "井浆", HT1_003_WELL_MUD_SLOW_VOLUME_M3, HT1_003_WELL_MUD_SLOW_RATE_M3_MIN,
                                remarks="井浆慢替 11.53m³@0.7m³/min，末段碰压。"),
        ),
        notes=(
            "施工顺序按 HT1-003 现场数据：先导浆→隔离液→领浆→中间浆→尾浆→压塞液→分段替浆。",
            "替浆总量 93.53m³，排量从 1.6 递减至 0.7m³/min。",
            "压塞液保留在 PumpingSchedule 中用于管内时序占位；环空入口分相映射时按 mud 相处理。",
        ),
    )

    validation_data = ValidationData(
        job_report_path=resolved_reference_root / "提取数据" / "呼1-003井_固井顶替模型数据包.json",
        notes=(
            "呼1-003首版加载器不读取 JSON 文件，仅把数据包中的现场参数固化为模块常量。",
            "CBL 评价井段和目标层段为暂定窗口，后续需结合 CBL 原始数据更新 ValidationData 路径。",
        ),
    )
    return well_spec, fluids, schedule, validation_data
