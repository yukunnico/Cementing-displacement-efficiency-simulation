"""
呼1-004井（HT1-004）168.3+139.7mm双径尾管段标准数据加载器

本模块把呼1-004井（HT1-004）现场施工设计文档中的 168.3+139.7mm 双径尾管段
资料整理为 cemdisp 标准输入结构。

数据来源：HT1-004井168.3+139.7mm油层尾管控压固井施工设计 (已审批).doc
  - [实测]：直接来自设计文档的实测数据
  - [代理-HT1-003]：复用HT1-003同类型流体的流变参数

当前求解目标为下部 139.7mm 尾管段，因此上部 168.3mm 重叠井段用等效井眼
直径保面积近似；鞋口滞后体积按四段管柱内径分段累加。
"""

from __future__ import annotations

import math
from pathlib import Path

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep
from cemdisp.data.validation_data import ValidationData
from cemdisp.data.well_spec import DepthValuePoint, EvaluationWindow, WellSpec


PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_REFERENCE_ROOT = PROJECT_ROOT / "参考文档" / "呼1-004"

# ===========================================================================
# 呼1-004 井段几何参数
# 来源：设计文档 第1.2节井身结构 + 第7.1.4节顶替量计算
# ===========================================================================
HT1_004_WELL_NAME = "呼1-004井（HT1-004）裸眼段"
HT1_004_DRILLED_DEPTH_MD_M = 7660.0       # 实际完钻井深/尾管鞋深度 (TD) [实测]
HT1_004_HANGER_MD_M = 5243.207            # 尾管悬挂器喇叭口位置 [实测]
HT1_004_TOP_MD_M = 5578.0                 # 模型剖面从套管鞋开始（仅裸眼段，不含套管重叠段）
HT1_004_CASING_SHOE_MD_M = 5578.0         # 273.1mm技术套管鞋深度 [实测]
HT1_004_UPPER_SECTION_BOTTOM_MD_M = 7378.051  # 168.3mm上段尾管底界/变径位置 [实测]
HT1_004_LOWER_HOLE_TOP_MD_M = 7521.0      # 241.3mm→215.9mm裸眼变径位置 [实测]
HT1_004_BOTTOM_MD_M = HT1_004_DRILLED_DEPTH_MD_M  # 下段 139.7mm 尾管鞋
HT1_004_SHOE_MD_M = HT1_004_DRILLED_DEPTH_MD_M

# --- 套管/尾管尺寸 [实测] ---
HT1_004_CASING_OD_MM = 273.1              # 技术套管外径
HT1_004_CASING_INNER_DIAMETER_MM = 245.37  # 273.1mm套管内径
HT1_004_UPPER_HOLE_NOMINAL_DIAMETER_MM = 241.3  # 上段井眼名义尺寸
HT1_004_LOWER_HOLE_NOMINAL_DIAMETER_MM = 215.9  # 下段井眼名义尺寸
HT1_004_UPPER_LINER_OD_MM = 168.3          # 上段尾管外径
HT1_004_LOWER_LINER_OD_MM = 139.7          # 下段尾管外径（基准OD）
HT1_004_UPPER_LINER_WALL_THICKNESS_MM = 14.7    # 168.3mm尾管壁厚 [实测]
HT1_004_LOWER_LINER_WALL_THICKNESS_MM = 15.88   # 139.7mm尾管壁厚 [实测]
HT1_004_UPPER_LINER_ID_MM = HT1_004_UPPER_LINER_OD_MM - 2.0 * HT1_004_UPPER_LINER_WALL_THICKNESS_MM  # 138.9mm
HT1_004_LOWER_LINER_ID_MM = HT1_004_LOWER_LINER_OD_MM - 2.0 * HT1_004_LOWER_LINER_WALL_THICKNESS_MM  # 107.94mm
HT1_004_LINER_ID_MM = HT1_004_LOWER_LINER_ID_MM

# --- 扶正器数据 [实测] 来源：设计文档第6.2节 ---
HT1_004_UPPER_CENTRALIZER_COUNT = 48       # 168.3mm段整体弹扶, 4根/只
HT1_004_MID_CENTRALIZER_COUNT = 14         # 139.7mm×241.3mm段, 1根/只
HT1_004_LOWER_CENTRALIZER_COUNT = 5        # 139.7mm×215.9mm段, 2根/只
HT1_004_CENTRALIZER_COUNT = (
    HT1_004_UPPER_CENTRALIZER_COUNT
    + HT1_004_MID_CENTRALIZER_COUNT
    + HT1_004_LOWER_CENTRALIZER_COUNT
)

# --- 平均居中度 [实测-模拟] 来源：设计文档第6.3节 ---
HT1_004_AVERAGE_STANDOFF = 0.83

# ===========================================================================
# 管柱内径剖面 (pipe_id_profile)
# 来源：设计文档第7.1.4节顶替量计算 [实测]
# ===========================================================================
_DP1_ID_MM = 149.2 - 2.0 * 9.65      # 129.9mm
_DP2_ID_MM = 127.0 - 2.0 * 9.65      # 107.7mm

# 钻具分段深度 [实测]
_DP1_BOTTOM_MD_M = 4025.734           # 149.2mm钻杆底界
_DP2_BOTTOM_MD_M = HT1_004_HANGER_MD_M  # 127mm钻杆底界 = 悬挂器


def _build_pipe_id_profile() -> tuple[tuple[float, float], ...]:
    """构建管柱内径剖面（深度, ID mm）用于1D前沿追踪。

    段1: 0~4025.734m       → 129.9mm (149.2mm钻杆, 壁厚9.65mm)
    段2: 4025.734~5243.207m → 107.7mm (127mm钻杆, 壁厚9.65mm)
    段3: 5243.207~7378.051m → 138.9mm (168.3mm尾管, 壁厚14.7mm)
    段4: 7378.051~7660m     → 107.94mm (139.7mm尾管, 壁厚15.88mm)
    """
    return (
        (0.001, _DP1_ID_MM),
        (_DP1_BOTTOM_MD_M, _DP1_ID_MM),
        (_DP1_BOTTOM_MD_M + 0.001, _DP2_ID_MM),
        (_DP2_BOTTOM_MD_M, _DP2_ID_MM),
        (_DP2_BOTTOM_MD_M + 0.001, HT1_004_UPPER_LINER_ID_MM),
        (HT1_004_UPPER_SECTION_BOTTOM_MD_M, HT1_004_UPPER_LINER_ID_MM),
        (HT1_004_UPPER_SECTION_BOTTOM_MD_M + 0.001, HT1_004_LOWER_LINER_ID_MM),
        (HT1_004_BOTTOM_MD_M, HT1_004_LOWER_LINER_ID_MM),
    )


# ===========================================================================
# 等效井径换算
# 用途：将上段168.3mm尾管环空面积等效为139.7mm基准外径下的虚拟井径
# ===========================================================================

def _equivalent_hole_diameter_mm(
    actual_hole_mm: float, actual_od_mm: float, reference_od_mm: float
) -> float:
    """在固定 reference_od 条件下，构造与原始环空面积等价的井眼直径。

    面积守恒：D_eq² - OD_ref² = D_actual² - OD_actual²
    [估算] 工程等效换算。
    """
    area_term = actual_hole_mm**2 - actual_od_mm**2 + reference_od_mm**2
    return math.sqrt(max(area_term, reference_od_mm**2))


# 套管段(5243.207-5578m)等效井径: 168.3mm尾管在273.1mm套管(ID=245.37mm)内
HT1_004_CASING_SECTION_EQUIVALENT_HOLE_MM = _equivalent_hole_diameter_mm(
    actual_hole_mm=HT1_004_CASING_INNER_DIAMETER_MM,
    actual_od_mm=HT1_004_UPPER_LINER_OD_MM,
    reference_od_mm=HT1_004_LOWER_LINER_OD_MM,
)

# 上裸眼段(5578-7378.051m)等效井径: 168.3mm尾管在241.3mm名义井眼内
HT1_004_UPPER_OPENHOLE_EQUIVALENT_MM = _equivalent_hole_diameter_mm(
    actual_hole_mm=HT1_004_UPPER_HOLE_NOMINAL_DIAMETER_MM,
    actual_od_mm=HT1_004_UPPER_LINER_OD_MM,
    reference_od_mm=HT1_004_LOWER_LINER_OD_MM,
)


# ===========================================================================
# 井径剖面 (hole_diameter_profile)  -- 从悬挂器到鞋底
# 来源：设计文档第1.4.2节电测井径 [实测]
# 裸眼段(5578m以下)使用电测数据，上段使用等效井径换算
# ===========================================================================

def _build_hole_diameter_profile() -> tuple[tuple[float, float], ...]:
    """构建从套管鞋(5578m)到鞋底(7660m)的裸眼段等效井径剖面。

    段1: 5578-7378.051m    → 上裸眼段等效井径 (168.3mm尾管→139.7mm基准)
    段2: 7378.051-7660m    → 下裸眼段直接使用电测井径 (139.7mm基准)
    """
    points: list[tuple[float, float]] = []

    # === 段1: 上裸眼段 (5578-7378.051m) — 168.3mm尾管在241.3mm井眼中 ===
    # 电测井径数据 [实测]，30m间隔，等效换算到139.7mm基准
    _raw_caliper_upper: list[tuple[float, float]] = [
        (5578.0, 250.82), (5630.0, 251.21), (5660.0, 250.95), (5690.0, 248.03),
        (5720.0, 244.48), (5750.0, 247.78), (5780.0, 248.67), (5810.0, 246.89),
        (5840.0, 248.92), (5870.0, 242.00), (5900.0, 246.00), (5930.0, 242.95),
        (5960.0, 249.94), (5990.0, 246.13), (6020.0, 244.09), (6050.0, 245.24),
        (6080.0, 242.70), (6110.0, 243.59), (6140.0, 245.49), (6170.0, 243.21),
        (6200.0, 247.40), (6230.0, 242.32), (6260.0, 246.63), (6290.0, 243.08),
        (6320.0, 242.00), (6350.0, 242.00), (6380.0, 242.00), (6410.0, 243.21),
        (6440.0, 241.68), (6470.0, 243.21), (6500.0, 242.32), (6530.0, 242.00),
        (6560.0, 245.49), (6600.0, 242.00), (6620.0, 242.00), (6650.0, 242.00),
        (6680.0, 242.00), (6710.0, 242.00), (6740.0, 242.00), (6770.0, 242.19),
        (6800.0, 241.93), (6830.0, 242.06), (6860.0, 242.82), (6890.0, 242.57),
        (6920.0, 244.22), (6950.0, 242.19), (6980.0, 242.00), (7010.0, 242.00),
        (7040.0, 242.00), (7070.0, 241.68), (7100.0, 242.00), (7130.0, 242.00),
        (7160.0, 242.00), (7190.0, 242.00), (7220.0, 242.00), (7250.0, 242.00),
        (7280.0, 242.00), (7310.0, 242.00), (7340.0, 242.00),
        (7378.051, 242.00),
    ]
    for depth, caliper_mm in _raw_caliper_upper:
        equiv = _equivalent_hole_diameter_mm(
            actual_hole_mm=caliper_mm,
            actual_od_mm=HT1_004_UPPER_LINER_OD_MM,
            reference_od_mm=HT1_004_LOWER_LINER_OD_MM,
        )
        points.append((depth, equiv))

    # === 段2: 下裸眼段 (7378.051-7660m) — 139.7mm尾管，直接使用电测井径 ===
    _raw_caliper_lower: list[tuple[float, float]] = [
        (7378.051, 242.00), (7380.0, 242.00), (7400.0, 242.00),
        (7430.0, 242.00), (7460.0, 242.00), (7490.0, 242.00),
        (7520.0, 216.00), (7550.0, 216.00), (7580.0, 216.00),
        (7610.0, 216.00), (7660.0, 216.00),
    ]
    for depth, caliper_mm in _raw_caliper_lower:
        points.append((depth, caliper_mm))

    return tuple(points)


# ===========================================================================
# 尾管外径剖面 (liner_od_profile)
# 来源：设计文档第1.2节 [实测]
# ===========================================================================

def _build_liner_od_profile() -> tuple[tuple[float, float], ...]:
    """构建尾管外径剖面（深度, OD mm），仅裸眼段(5578-7660m)。

    段1: 5578~7378.051m → 168.3mm (上段尾管)
    段2: 7378.051~7660m → 139.7mm (下段尾管)
    """
    return (
        (HT1_004_TOP_MD_M, HT1_004_UPPER_LINER_OD_MM),
        (HT1_004_UPPER_SECTION_BOTTOM_MD_M, HT1_004_UPPER_LINER_OD_MM),
        (HT1_004_UPPER_SECTION_BOTTOM_MD_M + 0.001, HT1_004_LOWER_LINER_OD_MM),
        (HT1_004_BOTTOM_MD_M, HT1_004_LOWER_LINER_OD_MM),
    )


# ===========================================================================
# 井斜剖面 (inclination_profile)
# 来源：设计文档第1.4.1节 [实测]
# ===========================================================================

def _build_inclination_profile() -> tuple[tuple[float, float], ...]:
    """构建井斜剖面（深度, 井斜角°）。

    30m间隔，0-7660m。该井基本为直井，最大井斜8.35°(7660m处)。
    """
    # [实测] 从文档1.4.1节井斜/方位表提取
    _raw_inclination: list[tuple[float, float]] = [
        (0.0, 0.0),
        (5600.0, 0.75), (5630.0, 0.75), (5660.0, 0.15), (5690.0, 0.21),
        (5720.0, 0.75), (5750.0, 1.10), (5780.0, 1.25), (5810.0, 1.23),
        (5840.0, 0.94), (5870.0, 0.52), (5900.0, 1.42), (5930.0, 1.03),
        (5960.0, 1.39), (5990.0, 0.88), (6020.0, 1.24), (6050.0, 0.95),
        (6080.0, 1.04), (6110.0, 1.44), (6140.0, 1.01), (6170.0, 1.04),
        (6200.0, 1.11), (6230.0, 1.06), (6260.0, 0.98), (6290.0, 1.24),
        (6320.0, 0.74), (6350.0, 1.28), (6380.0, 1.24), (6410.0, 1.68),
        (6440.0, 0.82), (6470.0, 1.84), (6500.0, 1.29), (6530.0, 1.16),
        (6560.0, 1.21), (6600.0, 1.74), (6620.0, 1.33), (6650.0, 2.02),
        (6680.0, 1.36), (6710.0, 1.82), (6740.0, 1.34), (6770.0, 1.43),
        (6800.0, 1.67), (6830.0, 1.78), (6860.0, 1.83), (6890.0, 1.32),
        (6920.0, 1.51), (6950.0, 1.59), (6980.0, 1.66), (7010.0, 1.35),
        (7040.0, 1.56), (7070.0, 1.62), (7100.0, 2.26), (7130.0, 2.28),
        (7160.0, 2.74), (7190.0, 3.13), (7220.0, 3.20), (7250.0, 3.34),
        (7280.0, 4.11), (7310.0, 4.68), (7340.0, 4.23), (7380.0, 4.63),
        (7400.0, 5.07), (7430.0, 5.48), (7460.0, 6.32), (7490.0, 7.01),
        (7520.0, 6.58), (7550.0, 6.92), (7580.0, 7.54), (7610.0, 7.91),
        (7660.0, 8.35),
    ]

    # 过滤到套管鞋(5578m)以下
    points: list[tuple[float, float]] = []
    for depth, inc in _raw_inclination:
        if depth < HT1_004_TOP_MD_M - 1.0:
            continue
        points.append((depth, inc))
    if not points or points[0][0] > HT1_004_TOP_MD_M:
        points.insert(0, (HT1_004_TOP_MD_M, 1.44))  # 5578m处井斜约1.44°
    if points[-1][0] < HT1_004_BOTTOM_MD_M:
        points.append((HT1_004_BOTTOM_MD_M, points[-1][1]))
    return tuple(points)


# ===========================================================================
# 居中度剖面
# 来源：设计文档第6.3节 — 平均居中度83% [实测-模拟]
# ===========================================================================

def _build_standoff_profile() -> tuple[tuple[float, float], ...]:
    """构建居中度剖面，裸眼段采用固定值0.83（与设计文档6.3节一致）。"""
    return (
        (HT1_004_TOP_MD_M, HT1_004_AVERAGE_STANDOFF),
        (HT1_004_UPPER_SECTION_BOTTOM_MD_M, HT1_004_AVERAGE_STANDOFF),
        (HT1_004_BOTTOM_MD_M, HT1_004_AVERAGE_STANDOFF),
    )


# ===========================================================================
# 鞋口滞后体积
# 来源：四段管柱内容积分段累加 [估算]
# ===========================================================================

def _pipe_volume_m3(length_m: float, inner_diameter_mm: float) -> float:
    """按内径和长度估算管内容积 (m³)。"""
    radius_m = inner_diameter_mm / 1000.0 / 2.0
    return math.pi * radius_m**2 * length_m


HT1_004_SHOE_LAG_VOLUME_M3 = (
    _pipe_volume_m3(_DP1_BOTTOM_MD_M, _DP1_ID_MM)                          # 149.2mm钻杆
    + _pipe_volume_m3(_DP2_BOTTOM_MD_M - _DP1_BOTTOM_MD_M, _DP2_ID_MM)     # 127mm钻杆
    + _pipe_volume_m3(                                                      # 168.3mm尾管
        HT1_004_UPPER_SECTION_BOTTOM_MD_M - HT1_004_HANGER_MD_M,
        HT1_004_UPPER_LINER_ID_MM,
    )
    + _pipe_volume_m3(                                                      # 139.7mm尾管
        HT1_004_BOTTOM_MD_M - HT1_004_UPPER_SECTION_BOTTOM_MD_M,
        HT1_004_LOWER_LINER_ID_MM,
    )
)


# ===========================================================================
# 流体密度参数 [实测] 来源：设计文档第7节浆柱结构设计
# ===========================================================================
HT1_004_MUD_DENSITY_KG_M3 = 1900.0            # 钻井液 密度1.90g/cm³ [实测]
HT1_004_LEAD_MUD_DENSITY_KG_M3 = 1750.0       # 先导浆 密度1.75g/cm³ [实测]
HT1_004_SPACER1_DENSITY_KG_M3 = 1950.0        # 隔离液1 密度1.95g/cm³ [实测]
HT1_004_SPACER2_DENSITY_KG_M3 = 1750.0        # 隔离液2 密度1.75g/cm³ [实测]
HT1_004_LEAD_DENSITY_KG_M3 = 1930.0           # 领浆 密度1.93g/cm³ [实测]
HT1_004_TAIL_DENSITY_KG_M3 = 1900.0           # 尾浆 密度1.90g/cm³ [实测]
HT1_004_PLUG_DENSITY_KG_M3 = 1700.0           # 压塞液 密度1.70g/cm³ [实测]
HT1_004_DISPLACEMENT_DENSITY_KG_M3 = 1900.0   # 替钻井液 密度1.90g/cm³ [实测]
HT1_004_BUFFER_DENSITY_KG_M3 = 1900.0         # 保护液 密度1.90g/cm³ [实测]
HT1_004_BASE_FLUID_DENSITY_KG_M3 = 1020.0     # 基液 密度1.02g/cm³ [实测]
HT1_004_WELL_MUD_DENSITY_KG_M3 = 1900.0       # 井浆(替浆) 密度1.90g/cm³ [实测]

# ===========================================================================
# 流体流变参数
# - 钻井液: [实测] 从设计文档六速读数推算 PV=53mPa·s, YP=8.5Pa
# - 先导浆: [实测] 文档要求 PV≤30mPa·s, YP≤7Pa
# - 水泥浆/隔离液/后置液: [代理-HT1-003] 复用HT1-003同类型流体参数
# ===========================================================================
HT1_004_MUD_PV_PA_S = 0.053          # 钻井液 PV=53mPa·s [实测]
HT1_004_MUD_YP_PA = 8.5              # YP=8.5Pa [实测]

HT1_004_LEAD_MUD_PV_PA_S = 0.030     # 先导浆 PV=30mPa·s [实测-设计要求]
HT1_004_LEAD_MUD_YP_PA = 7.0         # YP=7Pa [实测-设计要求]

HT1_004_SPACER_PV_PA_S = 0.060       # 隔离液 PV=60mPa·s [代理-HT1-003]
HT1_004_SPACER_YP_PA = 11.0          # YP=11Pa [代理-HT1-003]

HT1_004_LEAD_PV_PA_S = 0.160         # 领浆 PV=160mPa·s [代理-HT1-003]
HT1_004_LEAD_YP_PA = 13.0            # YP=13Pa [代理-HT1-003]

HT1_004_TAIL_PV_PA_S = 0.180         # 尾浆 PV=180mPa·s [代理-HT1-003]
HT1_004_TAIL_YP_PA = 14.0            # YP=14Pa [代理-HT1-003]

HT1_004_PLUG_PV_PA_S = 0.040         # 压塞液 PV=40mPa·s [代理-HT1-003]
HT1_004_PLUG_YP_PA = 9.0             # YP=9Pa [代理-HT1-003]

HT1_004_DISPLACEMENT_PV_PA_S = 0.040  # 替钻井液 PV=40mPa·s [代理-HT1-003]
HT1_004_DISPLACEMENT_YP_PA = 9.5      # YP=9.5Pa [代理-HT1-003]

HT1_004_BUFFER_PV_PA_S = 0.040       # 保护液 PV=40mPa·s [代理-HT1-003]
HT1_004_BUFFER_YP_PA = 9.2           # YP=9.2Pa [代理-HT1-003]

HT1_004_BASE_FLUID_PV_PA_S = 0.030   # 基液 PV=30mPa·s [代理-HT1-003]
HT1_004_BASE_FLUID_YP_PA = 9.0       # YP=9Pa [代理-HT1-003]

HT1_004_WELL_MUD_PV_PA_S = 0.030     # 井浆(替浆) PV=30mPa·s [代理-HT1-003]
HT1_004_WELL_MUD_YP_PA = 9.3         # YP=9.3Pa [代理-HT1-003]

# ===========================================================================
# 泵注量/排量参数 [实测] 来源：设计文档第7.2节施工工艺流程
# ===========================================================================
HT1_004_LEAD_MUD_VOLUME_M3 = 25.0       # 先导浆
HT1_004_SPACER1_VOLUME_M3 = 16.0        # 隔离液1
HT1_004_SPACER2_VOLUME_M3 = 10.0        # 隔离液2
HT1_004_LEAD_VOLUME_M3 = 48.0           # 领浆
HT1_004_TAIL_VOLUME_M3 = 28.0           # 尾浆
HT1_004_PLUG_VOLUME_M3 = 2.0            # 压塞液
HT1_004_DISPLACEMENT_FAST_VOLUME_M3 = 29.0   # 替钻井液(快)
HT1_004_BUFFER_VOLUME_M3 = 14.0         # 保护液
HT1_004_BASE_FLUID_VOLUME_M3 = 1.0      # 基液
HT1_004_WELL_MUD_1_VOLUME_M3 = 14.0     # 替浆段2 @1.2m³/min
HT1_004_WELL_MUD_2_VOLUME_M3 = 10.0     # 替浆段3 @1.0m³/min
HT1_004_WELL_MUD_3_VOLUME_M3 = 10.0     # 替浆段4 @0.9m³/min
HT1_004_WELL_MUD_4_VOLUME_M3 = 10.0     # 替浆段5 @0.8m³/min
HT1_004_WELL_MUD_5_VOLUME_M3 = 9.4      # 替浆段6 @0.7m³/min（原7.1+2.3补足，保证尾浆全入环空）

HT1_004_LEAD_MUD_RATE_M3_MIN = 1.4      # 先导浆
HT1_004_SPACER_RATE_M3_MIN = 1.2        # 隔离液
HT1_004_CEMENT_RATE_M3_MIN = 1.2        # 领浆/尾浆
HT1_004_PLUG_RATE_M3_MIN = 1.2          # 压塞液
HT1_004_DISPLACEMENT_FAST_RATE_M3_MIN = 1.5  # 替钻井液(快)
HT1_004_BUFFER_RATE_M3_MIN = 1.4        # 保护液
HT1_004_BASE_FLUID_RATE_M3_MIN = 1.4    # 基液
HT1_004_WELL_MUD_1_RATE_M3_MIN = 1.2    # 替浆段2
HT1_004_WELL_MUD_2_RATE_M3_MIN = 1.0    # 替浆段3
HT1_004_WELL_MUD_3_RATE_M3_MIN = 0.9    # 替浆段4
HT1_004_WELL_MUD_4_RATE_M3_MIN = 0.8    # 替浆段5
HT1_004_WELL_MUD_5_RATE_M3_MIN = 0.7    # 替浆段6


# ===========================================================================
# 辅助函数
# ===========================================================================

def _depth_points(values: tuple[tuple[float, float], ...]) -> tuple[DepthValuePoint, ...]:
    """把测深-数值元组转换为 WellSpec 使用的剖面点。"""
    return tuple(DepthValuePoint(depth_md_m=depth, value=value) for depth, value in values)


# ===========================================================================
# 主加载函数
# ===========================================================================

def load_ht1_004_tailpipe(
    *,
    reference_root: Path | None = None,
) -> tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载呼1-004井（HT1-004）168.3+139.7mm双径尾管段标准模型输入。

    Returns:
        (well_spec, fluids, schedule, validation_data)
    """

    resolved_reference_root = reference_root or DEFAULT_REFERENCE_ROOT

    # --- 构建剖面 ---
    hole_profile = _build_hole_diameter_profile()
    liner_od_profile = _build_liner_od_profile()
    pipe_id_profile = _build_pipe_id_profile()
    inc_profile = _build_inclination_profile()
    standoff_profile = _build_standoff_profile()

    # --- WellSpec ---
    well_spec = WellSpec(
        well_name=HT1_004_WELL_NAME,
        top_md_m=HT1_004_TOP_MD_M,
        bottom_md_m=HT1_004_BOTTOM_MD_M,
        shoe_md_m=HT1_004_SHOE_MD_M,
        hanger_md_m=None,  # 悬挂器在裸眼段之上，缩域后不属于模型域
        casing_id_mm=HT1_004_CASING_OD_MM,
        liner_od_mm=HT1_004_LOWER_LINER_OD_MM,
        liner_id_mm=HT1_004_LINER_ID_MM,
        liner_wall_thickness_mm=HT1_004_LOWER_LINER_WALL_THICKNESS_MM,
        shoe_lag_volume_m3=HT1_004_SHOE_LAG_VOLUME_M3,
        hole_diameter_profile=_depth_points(hole_profile),
        liner_od_profile=_depth_points(liner_od_profile),
        pipe_id_profile=_depth_points(pipe_id_profile),
        inclination_profile=_depth_points(inc_profile),
        standoff_profile=_depth_points(standoff_profile),
        # 双径尾管上段字段
        upper_section_bottom_md_m=HT1_004_UPPER_SECTION_BOTTOM_MD_M,
        upper_liner_od_mm=HT1_004_UPPER_LINER_OD_MM,
        upper_liner_id_mm=HT1_004_UPPER_LINER_ID_MM,
        evaluation_windows=(
            EvaluationWindow(
                name="CBL评价井段",
                top_md_m=HT1_004_CASING_SHOE_MD_M,
                bottom_md_m=HT1_004_BOTTOM_MD_M,
                window_type="cbl",
            ),
            EvaluationWindow(
                name="目标层段",
                top_md_m=7400.0,
                bottom_md_m=HT1_004_BOTTOM_MD_M,
                window_type="target",
            ),
        ),
        reference_root=resolved_reference_root,
        notes=(
            "呼1-004井（HT1-004）168.3+139.7mm双径尾管裸眼段(5578-7660m)控压固井。",
            "模型域仅覆盖裸眼段（套管鞋5578m→鞋底7660m），不含套管重叠段。",
            "井径剖面直接使用设计文档电测井径表实测值，上段等效换算至139.7mm基准OD。",
            "liner_od_profile：5578-7378m=168.3mm，7378-7660m=139.7mm。",
            "pipe_id_profile按4段管柱内径：129.9/107.7/138.9/107.94mm。",
            "钻井液/先导浆流变来自HT1-004文档实测；水泥浆/隔离液流变代理自HT1-003。",
            "居中度采用设计文档模拟值(平均83%)固定剖面。",
        ),
    )

    # --- 流体清单 ---
    fluids = (
        FluidSpec("钻井液", FluidRole.MUD, HT1_004_MUD_DENSITY_KG_M3,
                  RheologyModel.BINGHAM, HT1_004_MUD_PV_PA_S, HT1_004_MUD_YP_PA),
        FluidSpec("保护液", FluidRole.DISPLACEMENT, HT1_004_BUFFER_DENSITY_KG_M3,
                  RheologyModel.BINGHAM, HT1_004_BUFFER_PV_PA_S, HT1_004_BUFFER_YP_PA),
        FluidSpec("先导浆", FluidRole.WASH, HT1_004_LEAD_MUD_DENSITY_KG_M3,
                  RheologyModel.BINGHAM, HT1_004_LEAD_MUD_PV_PA_S, HT1_004_LEAD_MUD_YP_PA),
        FluidSpec("隔离液1", FluidRole.SPACER, HT1_004_SPACER1_DENSITY_KG_M3,
                  RheologyModel.BINGHAM, HT1_004_SPACER_PV_PA_S, HT1_004_SPACER_YP_PA),
        FluidSpec("隔离液2", FluidRole.SPACER, HT1_004_SPACER2_DENSITY_KG_M3,
                  RheologyModel.BINGHAM, HT1_004_SPACER_PV_PA_S, HT1_004_SPACER_YP_PA),
        FluidSpec("领浆", FluidRole.LEAD, HT1_004_LEAD_DENSITY_KG_M3,
                  RheologyModel.BINGHAM, HT1_004_LEAD_PV_PA_S, HT1_004_LEAD_YP_PA),
        FluidSpec("尾浆", FluidRole.TAIL, HT1_004_TAIL_DENSITY_KG_M3,
                  RheologyModel.BINGHAM, HT1_004_TAIL_PV_PA_S, HT1_004_TAIL_YP_PA),
        FluidSpec("压塞液", FluidRole.OTHER, HT1_004_PLUG_DENSITY_KG_M3,
                  RheologyModel.BINGHAM, HT1_004_PLUG_PV_PA_S, HT1_004_PLUG_YP_PA),
        FluidSpec("替钻井液", FluidRole.DISPLACEMENT, HT1_004_DISPLACEMENT_DENSITY_KG_M3,
                  RheologyModel.BINGHAM, HT1_004_DISPLACEMENT_PV_PA_S, HT1_004_DISPLACEMENT_YP_PA),
        FluidSpec("基液", FluidRole.DISPLACEMENT, HT1_004_BASE_FLUID_DENSITY_KG_M3,
                  RheologyModel.BINGHAM, HT1_004_BASE_FLUID_PV_PA_S, HT1_004_BASE_FLUID_YP_PA),
        FluidSpec("井浆", FluidRole.DISPLACEMENT, HT1_004_WELL_MUD_DENSITY_KG_M3,
                  RheologyModel.BINGHAM, HT1_004_WELL_MUD_PV_PA_S, HT1_004_WELL_MUD_YP_PA),
    )

    # --- 施工程序 ---
    # 注意：步骤7（释放胶塞）为停泵操作(pumping=false)，不计入 PumpingSchedule
    # 文档中替浆步骤9-16按分段降排量建模为独立步骤
    schedule = PumpingSchedule(
        steps=(
            PumpingScheduleStep("注入先导浆", "先导浆",
                                HT1_004_LEAD_MUD_VOLUME_M3, HT1_004_LEAD_MUD_RATE_M3_MIN,
                                remarks="先导浆 25m³@1.4m³/min，密度1.75g/cm³。"),
            PumpingScheduleStep("注入隔离液1", "隔离液1",
                                HT1_004_SPACER1_VOLUME_M3, HT1_004_SPACER_RATE_M3_MIN,
                                remarks="低失水驱油隔离液1 16m³@1.2m³/min，密度1.95g/cm³。"),
            PumpingScheduleStep("注入隔离液2", "隔离液2",
                                HT1_004_SPACER2_VOLUME_M3, HT1_004_SPACER_RATE_M3_MIN,
                                remarks="低失水驱油隔离液2 10m³@1.2m³/min，密度1.75g/cm³。"),
            PumpingScheduleStep("注入领浆", "领浆",
                                HT1_004_LEAD_VOLUME_M3, HT1_004_CEMENT_RATE_M3_MIN,
                                remarks="领浆 48m³@1.2m³/min，密度1.93g/cm³，设计占高1659m。"),
            PumpingScheduleStep("注入尾浆", "尾浆",
                                HT1_004_TAIL_VOLUME_M3, HT1_004_CEMENT_RATE_M3_MIN,
                                remarks="尾浆 28m³@1.2m³/min，密度1.90g/cm³，设计占高1060m。"),
            PumpingScheduleStep("替入压塞液", "压塞液",
                                HT1_004_PLUG_VOLUME_M3, HT1_004_PLUG_RATE_M3_MIN,
                                remarks="压塞液 2m³@1.2m³/min，密度1.70g/cm³。"),
            PumpingScheduleStep("替钻井液(快)", "替钻井液",
                                HT1_004_DISPLACEMENT_FAST_VOLUME_M3,
                                HT1_004_DISPLACEMENT_FAST_RATE_M3_MIN,
                                remarks="替钻井液 29m³@1.5m³/min，密度1.90g/cm³。"),
            PumpingScheduleStep("替保护液", "保护液",
                                HT1_004_BUFFER_VOLUME_M3, HT1_004_BUFFER_RATE_M3_MIN,
                                remarks="保护液 14m³@1.4m³/min，密度1.90g/cm³。"),
            PumpingScheduleStep("替基液", "基液",
                                HT1_004_BASE_FLUID_VOLUME_M3, HT1_004_BASE_FLUID_RATE_M3_MIN,
                                remarks="基液 1m³@1.4m³/min，密度1.02g/cm³。"),
            PumpingScheduleStep("井浆替入1", "井浆",
                                HT1_004_WELL_MUD_1_VOLUME_M3,
                                HT1_004_WELL_MUD_1_RATE_M3_MIN,
                                remarks="替浆段2 14m³@1.2m³/min，密度1.90g/cm³。"),
            PumpingScheduleStep("井浆替入2", "井浆",
                                HT1_004_WELL_MUD_2_VOLUME_M3,
                                HT1_004_WELL_MUD_2_RATE_M3_MIN,
                                remarks="替浆段3 10m³@1.0m³/min。"),
            PumpingScheduleStep("井浆替入3", "井浆",
                                HT1_004_WELL_MUD_3_VOLUME_M3,
                                HT1_004_WELL_MUD_3_RATE_M3_MIN,
                                remarks="替浆段4 10m³@0.9m³/min。"),
            PumpingScheduleStep("井浆替入4", "井浆",
                                HT1_004_WELL_MUD_4_VOLUME_M3,
                                HT1_004_WELL_MUD_4_RATE_M3_MIN,
                                remarks="替浆段5 10m³@0.8m³/min。"),
            PumpingScheduleStep("井浆替入5", "井浆",
                                HT1_004_WELL_MUD_5_VOLUME_M3,
                                HT1_004_WELL_MUD_5_RATE_M3_MIN,
                                remarks="替浆段6 7.1m³@0.7m³/min。"),
        ),
        notes=(
            "施工顺序：先导浆→隔离液1→隔离液2→领浆→尾浆→释放胶塞(停泵)"
            "→压塞液→替钻井液→保护液→基液→分段降排量替浆。",
            "替浆总量(不含压塞液): 29+14+1+14+10+10+10+9.4=97.4m³（末段+2.3补足尾浆全入环空）。",
            "释放胶塞步骤不计入PumpingSchedule（停泵操作）。",
            "分段降排量: 1.5→1.4→1.4→1.2→1.0→0.9→0.8→0.7 m³/min。",
        ),
    )

    validation_data = ValidationData(
        job_report_path=resolved_reference_root / "HT1-004井168.3+139.7mm油层尾管控压固井施工设计 (已审批) .doc",
        notes=(
            "呼1-004数据来自施工设计文档（已审批版），仅模拟裸眼段(5578-7660m)。",
            "井径数据来自设计文档1.4.2节电测井径实测值。",
            "井斜数据来自设计文档1.4.1节实测值。",
            "钻井液/先导浆流变来自HT1-004文档实测；水泥浆/隔离液流变代理自HT1-003。",
            "CBL评价井段和目标层段覆盖裸眼全段(5578-7660m)。",
        ),
    )
    return well_spec, fluids, schedule, validation_data
