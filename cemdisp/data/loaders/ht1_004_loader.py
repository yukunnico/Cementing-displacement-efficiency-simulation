"""
呼1-004井（HT1-004）168.3+139.7mm双径尾管段标准数据加载器

本模块读取呼1-004井身结构CSV实测数据，整理为 cemdisp 标准输入结构。

数据来源：
  - 呼1-004井身结构.csv（256行，30m间距实测数据）[实测]
  - HT1-004井施工设计文档（已审批）.doc [实测]
  - MATLAB压力计算脚本（流体流变/泵注参数）[参照]
  - HT1-003同类型流体参数 [代理-HT1-003]

当前求解目标为下部 139.7mm 尾管段，因此上部 168.3mm 重叠井段用等效井眼
直径保面积近似；鞋口滞后体积按四段管柱内径分段累加。
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
DEFAULT_REFERENCE_ROOT = PROJECT_ROOT / "参考文档" / "呼1-004"

# ===========================================================================
# 井段几何参数 [实测] 来源：设计文档 + 呼1-004井身结构.csv
# ===========================================================================
HT1_004_WELL_NAME = "呼1-004井（HT1-004）"
HT1_004_DRILLED_DEPTH_MD_M = 7660.0       # 实际完钻井深/TD
HT1_004_HANGER_MD_M = 5243.207            # 尾管悬挂器喇叭口
HT1_004_TOP_MD_M = HT1_004_HANGER_MD_M    # 模型域从悬挂器开始
HT1_004_CASING_SHOE_MD_M = 5578.0         # 273.1mm套管鞋
HT1_004_UPPER_SECTION_BOTTOM_MD_M = 7378.051  # 168.3→139.7mm尾管变径
HT1_004_LOWER_HOLE_TOP_MD_M = 7521.0      # 241.3→215.9mm井眼变径
HT1_004_BOTTOM_MD_M = HT1_004_DRILLED_DEPTH_MD_M
HT1_004_SHOE_MD_M = HT1_004_DRILLED_DEPTH_MD_M

# --- 套管/管柱尺寸 [实测] ---
HT1_004_CASING_OD_MM = 273.1
HT1_004_CASING_INNER_DIAMETER_MM = 245.37
HT1_004_UPPER_LINER_OD_MM = 168.3
HT1_004_LOWER_LINER_OD_MM = 139.7
HT1_004_UPPER_LINER_WALL_THICKNESS_MM = 15.88  # 来源：MATLAB脚本基础参数代码.docx
HT1_004_LOWER_LINER_WALL_THICKNESS_MM = 15.88
HT1_004_UPPER_LINER_ID_MM = HT1_004_UPPER_LINER_OD_MM - 2.0 * HT1_004_UPPER_LINER_WALL_THICKNESS_MM
HT1_004_LOWER_LINER_ID_MM = HT1_004_LOWER_LINER_OD_MM - 2.0 * HT1_004_LOWER_LINER_WALL_THICKNESS_MM
HT1_004_LINER_ID_MM = HT1_004_LOWER_LINER_ID_MM

# --- 平均居中度 [实测-模拟] 来源：设计文档第6.3节 ---
HT1_004_AVERAGE_STANDOFF = 0.83

# ===========================================================================
# CSV 剖面数据读取
# ===========================================================================

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


# ===========================================================================
# 等效井径换算
# ===========================================================================

def _equivalent_hole_diameter_mm(
    actual_hole_mm: float, actual_od_mm: float, reference_od_mm: float
) -> float:
    """在固定 reference_od 条件下，构造与原始环空面积等价的井眼直径。"""
    area_term = actual_hole_mm**2 - actual_od_mm**2 + reference_od_mm**2
    return math.sqrt(max(area_term, reference_od_mm**2))


# ===========================================================================
# 剖面构建（基于CSV数据）
# ===========================================================================

def _build_hole_diameter_profile(
    csv_rows: list[dict[str, float]],
) -> tuple[tuple[float, float], ...]:
    """从CSV环空体积数据反算等效井径剖面。

    CSV中 volume_annulus_L_ 已考虑实际管柱外径（149.2/127mm钻杆、168.3/139.7mm尾管），
    用此体积除以段长得到真实环空截面积，再反算等效井径（基于下段尾管OD=139.7mm）。
    这样2D求解器用 hole_diameter - 139.7mm 算出的间隙和体积与实际一致。
    """
    ref_od_m = HT1_004_LOWER_LINER_OD_MM / 1000.0
    ref_od_sq = ref_od_m ** 2

    points: list[tuple[float, float]] = []
    for row in csv_rows:
        depth = row["depth"]
        if depth < HT1_004_HANGER_MD_M - 0.1:
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

    if not points or points[0][0] > HT1_004_HANGER_MD_M:
        # 悬挂器处：273.1mm套管(ID=245.37) + 168.3mm尾管 → 139.7mm基准
        equiv0 = _equivalent_hole_diameter_mm(
            HT1_004_CASING_INNER_DIAMETER_MM,
            HT1_004_UPPER_LINER_OD_MM,
            HT1_004_LOWER_LINER_OD_MM,
        )
        points.insert(0, (HT1_004_HANGER_MD_M, equiv0))
    if points[-1][0] < HT1_004_BOTTOM_MD_M:
        points.append((HT1_004_BOTTOM_MD_M, points[-1][1]))

    return tuple(points)


def _build_inclination_profile(
    csv_rows: list[dict[str, float]],
) -> tuple[tuple[float, float], ...]:
    """从CSV数据构建井斜剖面。"""
    points: list[tuple[float, float]] = []
    for row in csv_rows:
        depth = row["depth"]
        if depth < HT1_004_HANGER_MD_M - 0.1:
            continue
        points.append((depth, row["inclination"]))
    if not points or points[0][0] > HT1_004_HANGER_MD_M:
        points.insert(0, (HT1_004_HANGER_MD_M, 0.75))
    if points[-1][0] < HT1_004_BOTTOM_MD_M:
        points.append((HT1_004_BOTTOM_MD_M, points[-1][1]))
    return tuple(points)


def _build_standoff_profile(
    csv_rows: list[dict[str, float]],
) -> tuple[tuple[float, float], ...]:
    """构建居中度剖面，全井段固定83%（设计文档6.3节模拟值）。"""
    return (
        (HT1_004_TOP_MD_M, HT1_004_AVERAGE_STANDOFF),
        (HT1_004_CASING_SHOE_MD_M, HT1_004_AVERAGE_STANDOFF),
        (HT1_004_UPPER_SECTION_BOTTOM_MD_M, HT1_004_AVERAGE_STANDOFF),
        (HT1_004_BOTTOM_MD_M, HT1_004_AVERAGE_STANDOFF),
    )


# ===========================================================================
# 管柱内径剖面 (pipe_id_profile) — 用于1D前沿追踪
# ===========================================================================
_DP1_ID_MM = 149.2 - 2.0 * 9.65      # 129.9mm
_DP2_ID_MM = 127.0 - 2.0 * 9.65      # 107.7mm
_DP1_BOTTOM_MD_M = 4025.734
_DP2_BOTTOM_MD_M = HT1_004_HANGER_MD_M


def _build_pipe_id_profile() -> tuple[tuple[float, float], ...]:
    """构建管柱内径剖面（深度, ID mm）。

    段1: 0~4025.734m        → 129.9mm (149.2mm钻杆)
    段2: 4025.734~5243.207m  → 107.7mm (127mm钻杆)
    段3: 5243.207~7378.051m  → 136.54mm (168.3mm尾管)
    段4: 7378.051~7660m      → 107.94mm (139.7mm尾管)
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


def _build_liner_od_profile() -> tuple[tuple[float, float], ...]:
    """构建尾管外径剖面（深度, OD mm）。

    段1: 5243.207~7378.051m → 168.3mm (上段尾管)
    段2: 7378.051~7660m     → 139.7mm (下段尾管)
    """
    return (
        (HT1_004_TOP_MD_M, HT1_004_UPPER_LINER_OD_MM),
        (HT1_004_UPPER_SECTION_BOTTOM_MD_M, HT1_004_UPPER_LINER_OD_MM),
        (HT1_004_UPPER_SECTION_BOTTOM_MD_M + 0.001, HT1_004_LOWER_LINER_OD_MM),
        (HT1_004_BOTTOM_MD_M, HT1_004_LOWER_LINER_OD_MM),
    )


# ===========================================================================
# 鞋口滞后体积
# ===========================================================================

def _pipe_volume_m3(length_m: float, inner_diameter_mm: float) -> float:
    """按内径和长度估算管内容积 (m³)。"""
    radius_m = inner_diameter_mm / 1000.0 / 2.0
    return math.pi * radius_m**2 * length_m


HT1_004_SHOE_LAG_VOLUME_M3 = (
    _pipe_volume_m3(_DP1_BOTTOM_MD_M, _DP1_ID_MM)
    + _pipe_volume_m3(_DP2_BOTTOM_MD_M - _DP1_BOTTOM_MD_M, _DP2_ID_MM)
    + _pipe_volume_m3(
        HT1_004_UPPER_SECTION_BOTTOM_MD_M - HT1_004_HANGER_MD_M,
        HT1_004_UPPER_LINER_ID_MM,
    )
    + _pipe_volume_m3(
        HT1_004_BOTTOM_MD_M - HT1_004_UPPER_SECTION_BOTTOM_MD_M,
        HT1_004_LOWER_LINER_ID_MM,
    )
)

# ===========================================================================
# 流体密度 [优化参数] 来源：优化参数.docx（2026-06-11）
# ===========================================================================
HT1_004_MUD_DENSITY_KG_M3 = 1900.0            # 钻井液 1.90（优化参数）
HT1_004_LEAD_MUD_DENSITY_KG_M3 = 1750.0       # 先导浆 1.75（优化参数）
HT1_004_SPACER1_DENSITY_KG_M3 = 1950.0        # 隔离液1 1.95（优化参数）
HT1_004_SPACER2_DENSITY_KG_M3 = 1750.0        # 隔离液2 1.75（优化参数）
HT1_004_LEAD_DENSITY_KG_M3 = 1930.0           # 领浆 1.93（优化参数）
HT1_004_TAIL_DENSITY_KG_M3 = 1900.0           # 尾浆 1.90（优化参数）
HT1_004_PLUG_DENSITY_KG_M3 = 1700.0           # 压塞液 1.70（优化参数）
HT1_004_DISPLACEMENT_DENSITY_KG_M3 = 1900.0   # 钻井液 1.90（优化参数）
HT1_004_BUFFER_DENSITY_KG_M3 = 1900.0         # 保护液 1.90（优化参数）
HT1_004_BASE_FLUID_DENSITY_KG_M3 = 1020.0     # 基液 1.02（优化参数）
HT1_004_WELL_MUD_DENSITY_KG_M3 = 1900.0       # 井浆(替浆段) 1.90（优化参数）

# ===========================================================================
# 流体流变参数 [优化参数] 来源：优化参数.docx + MATLAB基础参数代码.docx
# PV=mPa·s/1000→Pa·s, YP=Pa
# ===========================================================================
HT1_004_MUD_PV_PA_S = 0.053             # 井浆 PV=53mPa·s (MATLAB miu0)
HT1_004_MUD_YP_PA = 8.5                 # YP=8.5Pa (MATLAB tau0)

HT1_004_LEAD_MUD_PV_PA_S = 0.058        # 先导浆 PV=58mPa·s (MATLAB miu1)
HT1_004_LEAD_MUD_YP_PA = 9.8            # YP=9.8Pa (MATLAB tau1)

HT1_004_SPACER1_PV_PA_S = 0.058         # 隔离液1 PV=58mPa·s (MATLAB miu2)
HT1_004_SPACER1_YP_PA = 9.8             # YP=9.8Pa (MATLAB tau2)

HT1_004_SPACER2_PV_PA_S = 0.065         # 隔离液2 PV=65mPa·s (MATLAB miu3)
HT1_004_SPACER2_YP_PA = 10.0            # YP=10Pa (MATLAB tau3)

HT1_004_LEAD_PV_PA_S = 0.170            # 领浆 PV=170mPa·s (优化参数)
HT1_004_LEAD_YP_PA = 13.0               # YP=13Pa (优化参数)

HT1_004_TAIL_PV_PA_S = 0.180            # 尾浆 PV=180mPa·s (MATLAB miu5)
HT1_004_TAIL_YP_PA = 14.0               # YP=14Pa (MATLAB tau5)

HT1_004_PLUG_PV_PA_S = 0.050            # 压塞液 PV=50mPa·s (MATLAB miu6)
HT1_004_PLUG_YP_PA = 9.0                # YP=9Pa (MATLAB tau6)

HT1_004_DISPLACEMENT_PV_PA_S = 0.050    # 替钻井液 PV=50mPa·s (MATLAB miu7)
HT1_004_DISPLACEMENT_YP_PA = 9.5        # YP=9.5Pa (MATLAB tau7)

HT1_004_BUFFER_PV_PA_S = 0.050          # 保护液 PV=50mPa·s (MATLAB miu8)
HT1_004_BUFFER_YP_PA = 9.2              # YP=9.2Pa (MATLAB tau8)

HT1_004_BASE_FLUID_PV_PA_S = 0.050      # 基液 PV=50mPa·s (MATLAB miu9)
HT1_004_BASE_FLUID_YP_PA = 9.0          # YP=9Pa (MATLAB tau9)

HT1_004_WELL_MUD_PV_PA_S = 0.055        # 井浆(替浆段) PV=55mPa·s (MATLAB miu91-95)
HT1_004_WELL_MUD_YP_PA = 9.5            # YP=9.5Pa (MATLAB tau91-95)

# 合成 FLUSHER（冲洗液）参数：本井设计为先导浆+隔离液+领/尾浆，无 cement 前独立冲洗液，
# 故加入合成 FLUSHER 以验证 mud-spacer-flusher-cement 序列可表达。
HT1_004_FLUSHER_DENSITY_KG_M3 = 1880.0
HT1_004_FLUSHER_PV_PA_S = 0.050
HT1_004_FLUSHER_YP_PA = 9.0
HT1_004_FLUSHER_VOLUME_M3 = 5.0
HT1_004_FLUSHER_RATE_M3_MIN = 1.2

# ===========================================================================
# 泵注体积/排量 [优化参数] 来源：优化参数.docx（2026-06-11）
# 替浆分5级降排量：1.15→1.05→0.95→0.85→0.75m³/min，总量51.1m³
# 分段体积14+10+10+10+7.1=51.1m³，末段7.1+2.3补足=9.4m³
# ===========================================================================
HT1_004_LEAD_MUD_VOLUME_M3 = 25.0       # 优化：先导液 25m³
HT1_004_SPACER1_VOLUME_M3 = 16.0        # 优化：隔离液1 16m³
HT1_004_SPACER2_VOLUME_M3 = 10.0        # 优化：隔离液2 10m³
HT1_004_LEAD_VOLUME_M3 = 48.0           # 优化：领浆 48m³
HT1_004_TAIL_VOLUME_M3 = 28.0           # 优化：尾浆 28m³
HT1_004_PLUG_VOLUME_M3 = 2.0            # 优化：压塞液 2m³
HT1_004_DISPLACEMENT_FAST_VOLUME_M3 = 29.0   # 优化：钻井液 29m³
HT1_004_BUFFER_VOLUME_M3 = 14.0         # 优化：保护液 14m³
HT1_004_BASE_FLUID_VOLUME_M3 = 1.0      # 优化：基浆 1m³
HT1_004_WELL_MUD_1_VOLUME_M3 = 14.0     # 替浆段1 14m³ @1.15m³/min
HT1_004_WELL_MUD_2_VOLUME_M3 = 10.0     # 替浆段2 10m³ @1.05m³/min
HT1_004_WELL_MUD_3_VOLUME_M3 = 10.0     # 替浆段3 10m³ @0.95m³/min
HT1_004_WELL_MUD_4_VOLUME_M3 = 10.0     # 替浆段4 10m³ @0.85m³/min
HT1_004_WELL_MUD_5_VOLUME_M3 = 9.4      # 替浆段5 7.1+2.3补足=9.4m³ @0.75m³/min

HT1_004_LEAD_MUD_RATE_M3_MIN = 1.4      # 优化：先导液 1.4m³/min
HT1_004_SPACER_RATE_M3_MIN = 1.2        # 优化：隔离液 1.2m³/min
HT1_004_CEMENT_RATE_M3_MIN = 1.2        # 优化：领浆 1.2m³/min
HT1_004_TAIL_RATE_M3_MIN = 1.25         # 优化：尾浆 1.25m³/min
HT1_004_PLUG_RATE_M3_MIN = 1.2          # 优化：压塞液 1.2m³/min
HT1_004_DISPLACEMENT_FAST_RATE_M3_MIN = 1.5  # 优化：钻井液 1.5m³/min
HT1_004_BUFFER_RATE_M3_MIN = 1.4        # 优化：保护液 1.4m³/min
HT1_004_BASE_FLUID_RATE_M3_MIN = 1.4    # 优化：基浆 1.4m³/min
HT1_004_WELL_MUD_1_RATE_M3_MIN = 1.15   # 优化替浆第1段 1.15m³/min
HT1_004_WELL_MUD_2_RATE_M3_MIN = 1.05   # 优化替浆第2段 1.05m³/min
HT1_004_WELL_MUD_3_RATE_M3_MIN = 0.95   # 优化替浆第3段 0.95m³/min
HT1_004_WELL_MUD_4_RATE_M3_MIN = 0.85   # 优化替浆第4段 0.85m³/min
HT1_004_WELL_MUD_5_RATE_M3_MIN = 0.75   # 优化替浆第5段 0.75m³/min


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
    include_wash_spacer: bool = False,
) -> tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载呼1-004井（HT1-004）168.3+139.7mm双径尾管段标准模型输入。

    Args:
        reference_root: 可选参考资料根目录。
        include_wash_spacer: 是否注入合成 FLUSHER（冲洗液）步骤以验证
            mud-spacer-flusher-cement 序列可表达。默认为 False（严格现场模式）。

    Returns:
        (well_spec, fluids, schedule, validation_data)
    """

    resolved_reference_root = reference_root or DEFAULT_REFERENCE_ROOT

    # --- 读取实测井身结构CSV ---
    csv_path = resolved_reference_root / "呼1-004井身结构.csv"
    csv_rows = _read_well_structure_csv(csv_path)

    # --- 构建剖面 ---
    hole_profile = _build_hole_diameter_profile(csv_rows)
    inc_profile = _build_inclination_profile(csv_rows)
    standoff_profile = _build_standoff_profile(csv_rows)
    liner_od_profile = _build_liner_od_profile()
    pipe_id_profile = _build_pipe_id_profile()

    # --- WellSpec ---
    well_spec = WellSpec(
        well_name=HT1_004_WELL_NAME,
        top_md_m=HT1_004_TOP_MD_M,
        bottom_md_m=HT1_004_BOTTOM_MD_M,
        shoe_md_m=HT1_004_SHOE_MD_M,
        hanger_md_m=HT1_004_HANGER_MD_M,
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
            "呼1-004井（HT1-004）168.3+139.7mm双径尾管控压固井，井深7660m。",
            "井径/井斜剖面来自呼1-004井身结构.csv（256行30m间距实测）。",
            "等效井径按环空体积反算至139.7mm基准OD。",
            "pipe_id_profile按4段管柱内径：129.9/107.7/136.54/107.94mm。",
            "流体流变参数参照优化参数.docx + MATLAB基础参数代码.docx。",
            "泵注程序参照优化参数.docx（2026-06-11）；末段替浆+2.3m³补足。",
            "居中度采用设计文档6.3节模拟值(平均83%)固定剖面。",
        ),
    )

    # --- 流体清单（流变参照MATLAB脚本） ---
    fluids = (
        FluidSpec("钻井液", FluidRole.MUD, HT1_004_MUD_DENSITY_KG_M3,
                  RheologyModel.BINGHAM, HT1_004_MUD_PV_PA_S, HT1_004_MUD_YP_PA),
        FluidSpec("保护液", FluidRole.DISPLACEMENT, HT1_004_BUFFER_DENSITY_KG_M3,
                  RheologyModel.BINGHAM, HT1_004_BUFFER_PV_PA_S, HT1_004_BUFFER_YP_PA),
        FluidSpec("先导浆", FluidRole.WASH, HT1_004_LEAD_MUD_DENSITY_KG_M3,
                  RheologyModel.BINGHAM, HT1_004_LEAD_MUD_PV_PA_S, HT1_004_LEAD_MUD_YP_PA),
        FluidSpec("隔离液1", FluidRole.SPACER, HT1_004_SPACER1_DENSITY_KG_M3,
                  RheologyModel.BINGHAM, HT1_004_SPACER1_PV_PA_S, HT1_004_SPACER1_YP_PA),
        FluidSpec("隔离液2", FluidRole.SPACER, HT1_004_SPACER2_DENSITY_KG_M3,
                  RheologyModel.BINGHAM, HT1_004_SPACER2_PV_PA_S, HT1_004_SPACER2_YP_PA),
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
        FluidSpec("冲洗液（FLUSHER）", FluidRole.FLUSHER, HT1_004_FLUSHER_DENSITY_KG_M3,
                  RheologyModel.BINGHAM, HT1_004_FLUSHER_PV_PA_S, HT1_004_FLUSHER_YP_PA),
    )

    # --- 施工程序（参照MATLAB脚本 + 设计文档） ---
    # 合成 FLUSHER 步骤：仅在 include_wash_spacer=True 时注入，用于验证 mud-spacer-flusher-cement 序列可表达。
    flusher_step: tuple[PumpingScheduleStep, ...] = ()
    if include_wash_spacer:
        flusher_step = (
            PumpingScheduleStep("注入冲洗液（FLUSHER）", "冲洗液（FLUSHER）",
                                HT1_004_FLUSHER_VOLUME_M3, HT1_004_FLUSHER_RATE_M3_MIN,
                                remarks="合成冲洗液（FLUSHER）5m³@1.2m³/min，密度1.88g/cm³（验证序列可表达）。"),
        )

    schedule = PumpingSchedule(
        steps=(
            PumpingScheduleStep("注入先导浆", "先导浆",
                                HT1_004_LEAD_MUD_VOLUME_M3, HT1_004_LEAD_MUD_RATE_M3_MIN,
                                remarks="先导浆 25m³@1.4m³/min，密度1.75g/cm³（优化参数）。"),
            PumpingScheduleStep("注入隔离液1", "隔离液1",
                                HT1_004_SPACER1_VOLUME_M3, HT1_004_SPACER_RATE_M3_MIN,
                                remarks="隔离液1 16m³@1.2m³/min，密度1.95g/cm³（优化参数）。"),
            PumpingScheduleStep("注入隔离液2", "隔离液2",
                                HT1_004_SPACER2_VOLUME_M3, HT1_004_SPACER_RATE_M3_MIN,
                                remarks="隔离液2 10m³@1.2m³/min，密度1.75g/cm³（优化参数）。"),
            *flusher_step,
            PumpingScheduleStep("注入领浆", "领浆",
                                HT1_004_LEAD_VOLUME_M3, HT1_004_CEMENT_RATE_M3_MIN,
                                remarks="领浆 48m³@1.2m³/min，密度1.93g/cm³，PV=170mPa·s/YP=13Pa（优化参数）。"),
            PumpingScheduleStep("注入尾浆", "尾浆",
                                HT1_004_TAIL_VOLUME_M3, HT1_004_TAIL_RATE_M3_MIN,
                                remarks="尾浆 28m³@1.25m³/min，密度1.90g/cm³（优化参数）。"),
            PumpingScheduleStep("替入压塞液", "压塞液",
                                HT1_004_PLUG_VOLUME_M3, HT1_004_PLUG_RATE_M3_MIN,
                                remarks="压塞液 2m³@1.2m³/min，密度1.70g/cm³（优化参数）。"),
            PumpingScheduleStep("替钻井液(快)", "替钻井液",
                                HT1_004_DISPLACEMENT_FAST_VOLUME_M3,
                                HT1_004_DISPLACEMENT_FAST_RATE_M3_MIN,
                                remarks="钻井液 29m³@1.5m³/min，密度1.90g/cm³（优化参数）。"),
            PumpingScheduleStep("替保护液", "保护液",
                                HT1_004_BUFFER_VOLUME_M3, HT1_004_BUFFER_RATE_M3_MIN,
                                remarks="保护液 14m³@1.4m³/min，密度1.90g/cm³（优化参数）。"),
            PumpingScheduleStep("替基液", "基液",
                                HT1_004_BASE_FLUID_VOLUME_M3, HT1_004_BASE_FLUID_RATE_M3_MIN,
                                remarks="基液 1m³@1.4m³/min，密度1.02g/cm³（优化参数）。"),
            PumpingScheduleStep("井浆替入1", "井浆",
                                HT1_004_WELL_MUD_1_VOLUME_M3,
                                HT1_004_WELL_MUD_1_RATE_M3_MIN,
                                remarks="替浆段1 14m³@1.15m³/min（优化参数）。"),
            PumpingScheduleStep("井浆替入2", "井浆",
                                HT1_004_WELL_MUD_2_VOLUME_M3,
                                HT1_004_WELL_MUD_2_RATE_M3_MIN,
                                remarks="替浆段2 10m³@1.05m³/min（优化参数）。"),
            PumpingScheduleStep("井浆替入3", "井浆",
                                HT1_004_WELL_MUD_3_VOLUME_M3,
                                HT1_004_WELL_MUD_3_RATE_M3_MIN,
                                remarks="替浆段3 10m³@0.95m³/min（优化参数）。"),
            PumpingScheduleStep("井浆替入4", "井浆",
                                HT1_004_WELL_MUD_4_VOLUME_M3,
                                HT1_004_WELL_MUD_4_RATE_M3_MIN,
                                remarks="替浆段4 10m³@0.85m³/min（优化参数）。"),
            PumpingScheduleStep("井浆替入5", "井浆",
                                HT1_004_WELL_MUD_5_VOLUME_M3,
                                HT1_004_WELL_MUD_5_RATE_M3_MIN,
                                remarks="替浆段5 9.4m³@0.75m³/min（优化7.1+2.3m³补足尾浆全入环空）。"),
        ),
        notes=(
            "施工顺序：先导浆→隔离液1→隔离液2"
            + ("→冲洗液（FLUSHER）" if include_wash_spacer else "")
            + "→领浆→尾浆→压塞液→钻井液→保护液→基液→5级降排量替浆。",
            "替浆总量(不含压塞液): 29+14+1+14+10+10+10+9.4=97.4m³。",
            "泵注参数参照优化参数.docx（2026-06-11），末段补足2.3m³保证尾浆全入环空。",
            "排量: 1.4→1.2→1.2→1.2→1.25→1.2→1.5→1.4→1.4→1.15→1.05→0.95→0.85→0.75 m³/min。",
            "合成FLUSHER步骤仅在 include_wash_spacer=True 时注入，用于验证序列可表达。",
        ),
    )

    validation_data = ValidationData(
        job_report_path=resolved_reference_root / "呼1-004井身结构.csv",
        notes=(
            "呼1-004井身结构数据来自高精度CSV实测（256行30m间距）。",
            "井径/井斜剖面直接从CSV volume_annulus_L 反算等效井径。",
            "流体流变参数参照优化参数.docx + MATLAB基础参数代码.docx。",
            "泵注程序参照优化参数.docx五级降排量方案；末段补足2.3m³。",
            "CBL评价井段和目标层段覆盖裸眼全段(5578-7660m)。",
        ),
    )
    return well_spec, fluids, schedule, validation_data
