"""
呼1-003井（HT1-003）168.3+139.7mm双径尾管段标准数据加载器

本模块把呼1-003井（HT1-003）现场数据包中的 168.3+139.7mm 双径尾管段
资料整理为 cemdisp 标准输入结构。

当前求解目标为下部 139.7mm 尾管段，因此上部 168.3mm 重叠井段用等效井眼
直径保面积近似；鞋口滞后体积按 HT1-003 明确内径分段累加。
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep
from cemdisp.data.validation_data import ValidationData
from cemdisp.data.well_spec import DepthValuePoint, EvaluationWindow, WellSpec


PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_REFERENCE_ROOT = PROJECT_ROOT / "参考文档" / "呼1-003" / "新"

# 呼1-003井段几何参数，来源：参考文档/呼1-003/新/呼1-003井身结构.csv + untitled(1).m
HT1_003_WELL_NAME = "呼1-003井（HT1-003）"
HT1_003_DRILLED_DEPTH_MD_M = 7618.0  # 实际完钻井深/尾管鞋深度 (TD)
HT1_003_HANGER_MD_M = 5307.539  # 尾管悬挂器位置
HT1_003_TOP_MD_M = HT1_003_HANGER_MD_M  # 模型剖面从悬挂器开始
HT1_003_CASING_SHOE_MD_M = 5568.0  # 273.1mm套管鞋深度
HT1_003_UPPER_SECTION_BOTTOM_MD_M = 7089.576  # 168.3mm 上段尾管底界/变径位置
HT1_003_LOWER_SECTION_TOP_MD_M = 7096.0  # 241.3mm裸眼→215.9mm裸眼变径位置
HT1_003_BOTTOM_MD_M = HT1_003_DRILLED_DEPTH_MD_M  # 下段 139.7mm 尾管鞋
HT1_003_SHOE_MD_M = HT1_003_DRILLED_DEPTH_MD_M
HT1_003_CASING_ID_MM = 273.1  # 技术套管外径（名义）
HT1_003_CASING_INNER_DIAMETER_MM = 245.42  # 273.1mm套管内径
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


def _read_well_structure_csv(csv_path: Path) -> list[dict[str, float]]:
    """读取井身结构CSV，返回每行数据字典列表。"""
    rows: list[dict[str, float]] = []
    with csv_path.open(encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                "depth": float(row["depth_well_logging_m_"]),
                "length": float(row["length_segment_array_m_"]),
                "annulus_radius_cm": float(row["annulus_radius_array_cm_"]),
                "volume_annulus_L": float(row["volume_annulus_L_"]),
                "inclination": float(row["deg_for_logging_degree_"]),
            })
    return rows


def _build_hole_diameter_profile(
    csv_rows: list[dict[str, float]],
) -> tuple[tuple[float, float], ...]:
    """从CSV环空体积数据反算等效井径剖面。

    CSV中 volume_annulus_L_ 已考虑实际管柱（168.3mm/139.7mm尾管），
    用此体积除以段长得到真实环空截面积，再反算等效井径（基于下段尾管OD=139.7mm）。
    这样2D求解器用 hole_diameter - 139.7mm 算出的间隙和体积与实际一致。
    """
    ref_od_m = HT1_003_LOWER_LINER_OD_MM / 1000.0
    ref_od_sq = ref_od_m ** 2

    points: list[tuple[float, float]] = []
    for row in csv_rows:
        depth = row["depth"]
        if depth < HT1_003_HANGER_MD_M:
            continue
        seg_length = row["length"]
        vol_L = row["volume_annulus_L"]
        if seg_length > 0.0 and vol_L > 0.0:
            area_m2 = vol_L / 1000.0 / seg_length
            equiv_d_m = math.sqrt(4.0 * area_m2 / math.pi + ref_od_sq)
            equiv_d_mm = equiv_d_m * 1000.0
        else:
            equiv_d_mm = row["annulus_radius_cm"] * 20.0
        points.append((depth, equiv_d_mm))

    if not points or points[0][0] > HT1_003_HANGER_MD_M:
        area0 = math.pi * ((HT1_003_CASING_INNER_DIAMETER_MM / 1000.0) ** 2 - (HT1_003_UPPER_LINER_OD_MM / 1000.0) ** 2) / 4.0
        equiv0_mm = math.sqrt(4.0 * area0 / math.pi + ref_od_sq) * 1000.0
        points.insert(0, (HT1_003_HANGER_MD_M, equiv0_mm))
    if points[-1][0] < HT1_003_BOTTOM_MD_M:
        points.append((HT1_003_BOTTOM_MD_M, points[-1][1]))

    return tuple(points)


def _build_liner_od_profile() -> tuple[tuple[float, float], ...]:
    """构建尾管外径剖面（深度, OD mm）。

    段1: 5307.539~7089.576m → 168.3mm (上段尾管)
    段2: 7089.576~7618m     → 139.7mm (下段尾管)
    """
    return (
        (HT1_003_TOP_MD_M, HT1_003_UPPER_LINER_OD_MM),
        (HT1_003_UPPER_SECTION_BOTTOM_MD_M, HT1_003_UPPER_LINER_OD_MM),
        (HT1_003_UPPER_SECTION_BOTTOM_MD_M + 0.001, HT1_003_LOWER_LINER_OD_MM),
        (HT1_003_BOTTOM_MD_M, HT1_003_LOWER_LINER_OD_MM),
    )


def _build_pipe_id_profile() -> tuple[tuple[float, float], ...]:
    """构建管柱内径剖面（深度, ID mm）。

    从地面到TD的完整管柱内径，用于1D前沿追踪计算管内容积。
    段1: 0~3321.682m      → 129.9mm (149.2mm钻杆, 壁厚9.65mm)
    段2: 3321.682~5307.539m → 107.7mm (127mm钻杆, 壁厚9.65mm)
    段3: 5307.539~7089.576m → 138.9mm (168.3mm尾管, 壁厚14.7mm)
    段4: 7089.576~7618m    → 107.94mm (139.7mm尾管, 壁厚15.88mm)
    """
    dp1_id = 149.2 - 2 * 9.65   # 129.9mm
    dp2_id = 127.0 - 2 * 9.65   # 107.7mm
    liner1_id = 168.3 - 2 * 14.7  # 138.9mm
    liner2_id = 139.7 - 2 * 15.88  # 107.94mm
    return (
        (0.001, dp1_id),
        (3321.682, dp1_id),
        (3321.682 + 0.001, dp2_id),
        (5307.539, dp2_id),
        (5307.539 + 0.001, liner1_id),
        (HT1_003_UPPER_SECTION_BOTTOM_MD_M, liner1_id),
        (HT1_003_UPPER_SECTION_BOTTOM_MD_M + 0.001, liner2_id),
        (HT1_003_BOTTOM_MD_M, liner2_id),
    )


def _build_inclination_profile(
    csv_rows: list[dict[str, float]],
) -> tuple[tuple[float, float], ...]:
    """从CSV数据构建井斜剖面。"""
    points: list[tuple[float, float]] = []
    for row in csv_rows:
        depth = row["depth"]
        if depth < HT1_003_HANGER_MD_M:
            continue
        points.append((depth, row["inclination"]))
    if not points or points[0][0] > HT1_003_HANGER_MD_M:
        points.insert(0, (HT1_003_HANGER_MD_M, 0.0))
    if points[-1][0] < HT1_003_BOTTOM_MD_M:
        points.append((HT1_003_BOTTOM_MD_M, points[-1][1]))
    return tuple(points)


def _equivalent_hole_diameter_mm(actual_hole_mm: float, actual_od_mm: float, reference_od_mm: float) -> float:
    """在固定 reference_od 条件下，构造与原始环空面积等价的井眼直径。"""
    area_term = actual_hole_mm**2 - actual_od_mm**2 + reference_od_mm**2
    return math.sqrt(max(area_term, reference_od_mm**2))


HT1_003_UPPER_HOLE_DIAMETER_MM = _equivalent_hole_diameter_mm(
    actual_hole_mm=HT1_003_UPPER_HOLE_NOMINAL_DIAMETER_MM,
    actual_od_mm=HT1_003_UPPER_LINER_OD_MM,
    reference_od_mm=HT1_003_LOWER_LINER_OD_MM,
)
# 套管段等效井径：将168.3mm尾管在273.1mm套管内的环空面积，等效为139.7mm尾管的虚拟井径
HT1_003_CASING_SECTION_EQUIVALENT_HOLE_MM = _equivalent_hole_diameter_mm(
    actual_hole_mm=HT1_003_CASING_INNER_DIAMETER_MM,
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

# 呼1-003现场流体参数（来源：untitled(1).m 8.11施工过程模拟）
HT1_003_MUD_DENSITY_KG_M3 = 1950.0  # 钻井液(井浆) rou0
HT1_003_BALANCE_DENSITY_KG_M3 = 1850.0  # 先导浆 密度1.85g/cm³
HT1_003_SPACER_DENSITY_KG_M3 = 2000.0  # 驱油隔离液 密度2.00g/cm³
HT1_003_LEAD_DENSITY_KG_M3 = 2050.0  # 领浆 rou4
HT1_003_TAIL_DENSITY_KG_M3 = 2050.0  # 尾浆 rou5（施工设计表）
HT1_003_PLUG_DENSITY_KG_M3 = 1950.0  # 压塞液 rou6
HT1_003_DISPLACEMENT_DENSITY_KG_M3 = 1950.0  # 替浆钻井液 rou7
HT1_003_BUFFER_DENSITY_KG_M3 = 1950.0  # 保护液 rou8
HT1_003_BASE_FLUID_DENSITY_KG_M3 = 1020.0  # 基液 密度1.02g/cm³
HT1_003_WELL_MUD_DENSITY_KG_M3 = 1850.0  # 替浆 密度1.85g/cm³

# 钻井液流变参数（施工设计提供）
HT1_003_MUD_PV_PA_S = 0.051  # 51 mPa·s
HT1_003_MUD_YP_PA = 10.0
HT1_003_MUD_POWER_LAW_N = 0.82
HT1_003_MUD_CONSISTENCY_K = 0.21

# 流体流变参数（来源：untitled(1).m miu/tau参数，Bingham代理）
HT1_003_BALANCE_PV_PA_S = 0.051  # 先导浆 PV=51 mPa·s
HT1_003_BALANCE_YP_PA = 10.0  # YP=10 Pa
HT1_003_SPACER_PV_PA_S = 0.060  # 隔离液 PV=60 mPa·s
HT1_003_SPACER_YP_PA = 11.0  # YP=11 Pa
HT1_003_LEAD_PV_PA_S = 0.160  # 领浆 PV=160 mPa·s
HT1_003_LEAD_YP_PA = 13.0  # YP=13 Pa
HT1_003_TAIL_PV_PA_S = 0.180  # 尾浆 miu5=180 mPa·s
HT1_003_TAIL_YP_PA = 14.0  # tau5=14 Pa
HT1_003_PLUG_PV_PA_S = 0.040  # 压塞液 miu6=40 mPa·s
HT1_003_PLUG_YP_PA = 9.0  # tau6=9 Pa
HT1_003_DISPLACEMENT_PV_PA_S = 0.040  # 替浆钻井液 miu7=40 mPa·s
HT1_003_DISPLACEMENT_YP_PA = 9.5  # tau7=9.5 Pa
HT1_003_BUFFER_PV_PA_S = 0.040  # 保护液 miu8=40 mPa·s
HT1_003_BUFFER_YP_PA = 9.2  # tau8=9.2 Pa
HT1_003_BASE_FLUID_PV_PA_S = 0.030  # 基液 PV=30 mPa·s
HT1_003_BASE_FLUID_YP_PA = 9.0  # YP=9 Pa
HT1_003_WELL_MUD_PV_PA_S = 0.030  # 替浆 PV=30 mPa·s
HT1_003_WELL_MUD_YP_PA = 9.3  # YP=9.3 Pa

# 呼1-003现场施工程序参数（来源：untitled(1).m 8.11施工过程模拟）
HT1_003_BALANCE_VOLUME_M3 = 35.0  # 先导浆 35m³
HT1_003_SPACER_VOLUME_M3 = 35.0  # 驱油隔离液 35m³
HT1_003_LEAD_VOLUME_M3 = 38.0  # 领浆 38m³
HT1_003_TAIL_VOLUME_M3 = 28.5  # 尾浆 28.5m³
HT1_003_PLUG_VOLUME_M3 = 2.0  # 压塞液 2m³
HT1_003_FAST_MUD_VOLUME_M3 = 28.0  # 替浆钻井液 28m³
HT1_003_BUFFER_VOLUME_M3 = 12.0  # 保护液 12m³
HT1_003_BASE_FLUID_VOLUME_M3 = 3.0  # 基液 3m³
HT1_003_WELL_MUD_FAST_VOLUME_M3 = 8.0  # 替浆段1 @1.3m³/min
HT1_003_WELL_MUD_MID1_VOLUME_M3 = 14.0  # 替浆段2 @1.1m³/min
HT1_003_WELL_MUD_MID2_VOLUME_M3 = 14.0  # 替浆段3 @0.9m³/min
HT1_003_WELL_MUD_MID3_VOLUME_M3 = 12.0  # 替浆段4 @0.7m³/min
HT1_003_WELL_MUD_SLOW_VOLUME_M3 = 0.0  # 保留段

# 排量参数（来源：untitled(1).m pump_rate参数）
HT1_003_BALANCE_RATE_M3_MIN = 1.2  # 先导浆
HT1_003_SPACER_RATE_M3_MIN = 1.2  # 隔离液
HT1_003_CEMENT_RATE_M3_MIN = 1.2  # 领浆
HT1_003_TAIL_RATE_M3_MIN = 1.2  # 尾浆
HT1_003_PLUG_RATE_M3_MIN = 1.3  # 压塞液
HT1_003_FAST_MUD_RATE_M3_MIN = 1.3  # 替浆钻井液
HT1_003_BUFFER_RATE_M3_MIN = 1.2  # 保护液
HT1_003_BASE_FLUID_RATE_M3_MIN = 1.2  # 基液
HT1_003_WELL_MUD_FAST_RATE_M3_MIN = 1.3  # 替浆段1
HT1_003_WELL_MUD_MID1_RATE_M3_MIN = 1.1  # 替浆段2
HT1_003_WELL_MUD_MID2_RATE_M3_MIN = 0.9  # 替浆段3
HT1_003_WELL_MUD_MID3_RATE_M3_MIN = 0.7  # 替浆段4
HT1_003_WELL_MUD_SLOW_RATE_M3_MIN = 0.7  # 保留段


def _depth_points(values: tuple[tuple[float, float], ...]) -> tuple[DepthValuePoint, ...]:
    """把测深-数值元组转换为 WellSpec 使用的剖面点。"""
    return tuple(DepthValuePoint(depth_md_m=depth, value=value) for depth, value in values)


def load_ht1_003_tailpipe(
    *,
    reference_root: Path | None = None,
) -> tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载呼1-003井（HT1-003）168.3+139.7mm双径尾管段标准模型输入。"""

    resolved_reference_root = reference_root or DEFAULT_REFERENCE_ROOT

    # 读取实测井身结构数据
    csv_path = resolved_reference_root / "呼1-003井身结构.csv"
    csv_rows = _read_well_structure_csv(csv_path)

    # 从CSV构建井径剖面和井斜剖面
    hole_profile = _build_hole_diameter_profile(csv_rows)
    inc_profile = _build_inclination_profile(csv_rows)
    liner_od_profile = _build_liner_od_profile()
    pipe_id_profile = _build_pipe_id_profile()

    well_spec = WellSpec(
        well_name=HT1_003_WELL_NAME,
        top_md_m=HT1_003_TOP_MD_M,
        bottom_md_m=HT1_003_BOTTOM_MD_M,
        shoe_md_m=HT1_003_SHOE_MD_M,
        hanger_md_m=HT1_003_HANGER_MD_M,
        casing_id_mm=HT1_003_CASING_ID_MM,
        liner_od_mm=HT1_003_LOWER_LINER_OD_MM,
        liner_id_mm=HT1_003_LINER_ID_MM,
        shoe_lag_volume_m3=HT1_003_SHOE_LAG_VOLUME_M3,
        hole_diameter_profile=_depth_points(hole_profile),
        liner_od_profile=_depth_points(liner_od_profile),
        pipe_id_profile=_depth_points(pipe_id_profile),
        inclination_profile=_depth_points(inc_profile),
        standoff_profile=_depth_points(
            (
                (HT1_003_TOP_MD_M, 0.83),
                (5500.0, 0.83),
                (HT1_003_CASING_SHOE_MD_M, 0.83),
                (6000.0, 0.83),
                (6600.0, 0.83),
                (HT1_003_UPPER_SECTION_BOTTOM_MD_M, 0.83),
                (7400.0, 0.83),
                (7500.0, 0.83),
                (HT1_003_BOTTOM_MD_M, 0.83),
            )
        ),
        evaluation_windows=(
            EvaluationWindow(name="CBL评价井段", top_md_m=5568.0, bottom_md_m=7618.0, window_type="cbl"),
            EvaluationWindow(name="目标层段", top_md_m=7400.0, bottom_md_m=7618.0, window_type="target"),
        ),
        reference_root=resolved_reference_root,
        notes=(
            "呼1-003井（HT1-003）为168.3mm+139.7mm双径复合尾管控压固井，井深7618m。",
            "井径剖面直接使用CSV实测值（annulus_radius_array_cm_ × 20 = 井径mm）。",
            "liner_od_profile 按深度分段：5307-7089m=168.3mm，7089-7618m=139.7mm。",
            "pipe_id_profile 按4段管柱内径：129.9/107.7/138.9/107.94mm。",
            "流体参数来源：untitled(1).m 8.11施工过程模拟。",
        ),
    )

    # HT1-003 流体清单（来源：untitled(1).m）
    fluids = (
        FluidSpec("钻井液", FluidRole.MUD, HT1_003_MUD_DENSITY_KG_M3, RheologyModel.POWER_LAW,
                  power_law_n=HT1_003_MUD_POWER_LAW_N, consistency_k=HT1_003_MUD_CONSISTENCY_K),
        FluidSpec("保护液", FluidRole.DISPLACEMENT, HT1_003_BUFFER_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_BUFFER_PV_PA_S, HT1_003_BUFFER_YP_PA),
        FluidSpec("平衡液", FluidRole.WASH, HT1_003_BALANCE_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_BALANCE_PV_PA_S, HT1_003_BALANCE_YP_PA),
        FluidSpec("隔离液", FluidRole.SPACER, HT1_003_SPACER_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_SPACER_PV_PA_S, HT1_003_SPACER_YP_PA),
        FluidSpec("领浆", FluidRole.LEAD, HT1_003_LEAD_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_LEAD_PV_PA_S, HT1_003_LEAD_YP_PA),
        FluidSpec("尾浆", FluidRole.TAIL, HT1_003_TAIL_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_TAIL_PV_PA_S, HT1_003_TAIL_YP_PA),
        FluidSpec("压塞液", FluidRole.OTHER, HT1_003_PLUG_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_PLUG_PV_PA_S, HT1_003_PLUG_YP_PA),
        FluidSpec("替钻井液", FluidRole.DISPLACEMENT, HT1_003_DISPLACEMENT_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_DISPLACEMENT_PV_PA_S, HT1_003_DISPLACEMENT_YP_PA),
        FluidSpec("基液", FluidRole.DISPLACEMENT, HT1_003_BASE_FLUID_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_BASE_FLUID_PV_PA_S, HT1_003_BASE_FLUID_YP_PA),
        FluidSpec("井浆", FluidRole.DISPLACEMENT, HT1_003_WELL_MUD_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_WELL_MUD_PV_PA_S, HT1_003_WELL_MUD_YP_PA),
    )

    # HT1-003 地面施工程序（来源：untitled(1).m 8.11施工过程模拟）
    schedule = PumpingSchedule(
        steps=(
            PumpingScheduleStep("注入先导浆", "平衡液", HT1_003_BALANCE_VOLUME_M3, HT1_003_BALANCE_RATE_M3_MIN,
                                remarks="先导浆 35m³@1.2m³/min，密度1.85g/cm³。"),
            PumpingScheduleStep("注入驱油隔离液", "隔离液", HT1_003_SPACER_VOLUME_M3, HT1_003_SPACER_RATE_M3_MIN,
                                remarks="驱油隔离液 35m³@1.2m³/min，密度2.00g/cm³。"),
            PumpingScheduleStep("注入领浆", "领浆", HT1_003_LEAD_VOLUME_M3, HT1_003_CEMENT_RATE_M3_MIN,
                                remarks="领浆 38m³@1.2m³/min，密度2.05g/cm³。"),
            PumpingScheduleStep("注入尾浆", "尾浆", HT1_003_TAIL_VOLUME_M3, HT1_003_TAIL_RATE_M3_MIN,
                                remarks="尾浆 28.5m³@1.2m³/min，密度1.95g/cm³。"),
            PumpingScheduleStep("注入压塞液（管内）", "压塞液", HT1_003_PLUG_VOLUME_M3, HT1_003_PLUG_RATE_M3_MIN,
                                remarks="压塞液 2m³@1.3m³/min，仅作为管内占位。"),
            PumpingScheduleStep("替钻井液(快)", "替钻井液", HT1_003_FAST_MUD_VOLUME_M3, HT1_003_FAST_MUD_RATE_M3_MIN,
                                remarks="替钻井液 28m³@1.3m³/min。"),
            PumpingScheduleStep("替保护液", "保护液", HT1_003_BUFFER_VOLUME_M3, HT1_003_BUFFER_RATE_M3_MIN,
                                remarks="保护液 12m³@1.2m³/min，密度1.80g/cm³。"),
            PumpingScheduleStep("替基液", "基液", HT1_003_BASE_FLUID_VOLUME_M3, HT1_003_BASE_FLUID_RATE_M3_MIN,
                                remarks="基液 3m³@1.2m³/min，密度1.02g/cm³。"),
            PumpingScheduleStep("井浆快替1", "井浆", HT1_003_WELL_MUD_FAST_VOLUME_M3, HT1_003_WELL_MUD_FAST_RATE_M3_MIN,
                                remarks="替浆段1 8m³@1.3m³/min。"),
            PumpingScheduleStep("井浆中替1", "井浆", HT1_003_WELL_MUD_MID1_VOLUME_M3, HT1_003_WELL_MUD_MID1_RATE_M3_MIN,
                                remarks="替浆段2 14m³@1.1m³/min。"),
            PumpingScheduleStep("井浆中替2", "井浆", HT1_003_WELL_MUD_MID2_VOLUME_M3, HT1_003_WELL_MUD_MID2_RATE_M3_MIN,
                                remarks="替浆段3 14m³@0.9m³/min。"),
            PumpingScheduleStep("井浆中替3", "井浆", HT1_003_WELL_MUD_MID3_VOLUME_M3, HT1_003_WELL_MUD_MID3_RATE_M3_MIN,
                                remarks="替浆段4 12m³@0.7m³/min。"),
        ),
        notes=(
            "施工顺序：先导浆→隔离液→领浆→尾浆→压塞液→分段替浆。",
            "替浆总量 48m³，排量从 1.3 递减至 0.7m³/min。",
            "压塞液保留在 PumpingSchedule 中用于管内时序占位。",
        ),
    )

    validation_data = ValidationData(
        job_report_path=resolved_reference_root / "呼1-003井身结构.csv",
        notes=(
            "呼1-003更新版加载器使用新数据包（参考文档/呼1-003/新/）。",
            "井径数据来自呼1-003井身结构CSV实测值。",
            "CBL 评价井段和目标层段覆盖尾管全段(5568-7618m)。",
        ),
    )
    return well_spec, fluids, schedule, validation_data
