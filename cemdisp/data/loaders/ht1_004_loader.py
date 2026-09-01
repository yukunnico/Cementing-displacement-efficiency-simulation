"""
呼1-004井（HT1-004）168.3+139.7mm双径尾管段标准数据加载器

本模块读取呼1-004井现场提取包资料，整理为 cemdisp 标准输入结构。

数据来源（2026-08-16 核对报告 00_执行记录/ht1_003_004_loader核对_2026-08-16.md 落地）：
  - 井径/井斜：现场提取包 caliper_profile.csv / inclination_profile.csv（与 hu102/hu103 一致）
  - 井身结构：参考文档/呼1-004/呼1-004井身结构.csv + 施工设计文档（已审批）.doc
  - 流体密度/泵注体积：现场提取包 fluid_properties.csv / pumping_schedule.csv（field_measured）
  - 流变/排量：**优化参数.docx（2026-06-11）——非现场实测**，本模块默认输入即"优化参数化输入"，
    不可写成现场实测验证；论文使用应定位为"优化/应用案例"
  - 实际施工版（load_ht1_004_tailpipe_actual）：按施工记录表/作业史实际排量与替浆 97.1m³ 分 5 段重建

2026-08-16 更新要点：
- 168.3mm 上段壁厚 15.88→14.7（ID 136.54→138.9，现场设计 6.1）；鞋口滞后体积随之复核（~99.4m³）
- 变径深度 7378.051→7376.656（变扣 7376.656m，变径 BGT2*BG-FJU@7378.051m 双口径注明）
- CBL：补 cbl_pass_rate=0.003（尾管段 5245-7581m 数字化口径）；官方全井 0.2999 记入 notes；
  CBL 窗 5578-7660→5245-7581；目标窗 7400-7660→目的井段 7495-7526+7531-7550（两段）+ 红线 7482-7560
- 数据标签落地：密度 field_measured、流变/排量 optimized_input、合成 FLUSHER model_assumption（现场无此流体）
- 井径/井斜剖面来源切到现场提取包（原呼1-004井身结构.csv 反算口径以 LEGACY 保留）

当前求解目标为下部 139.7mm 尾管段，因此上部 168.3mm 重叠井段用等效井眼
直径保面积近似；鞋口滞后体积按四段管柱内径分段累加。
"""

from __future__ import annotations

import csv
import math
import warnings
from pathlib import Path

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep
from cemdisp.data.validation_data import ValidationData
from cemdisp.data.well_spec import DepthValuePoint, EvaluationWindow, WellSpec


PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_REFERENCE_ROOT = PROJECT_ROOT / "参考文档" / "呼1-004"
# 现场提取包 caliper_profile.csv：设计一.1.4.2 电测井径表 72 段（含套管段 245.37mm 与裸眼段）。
# LEGACY(2026-08-16 前): 参考文档/呼1-004/呼1-004井身结构.csv（环空体积反算等效井径，非现场提取包）。
DEFAULT_CALIPER_CSV = PROJECT_ROOT / "参考文档" / "现场资料提取" / "ht1_004_呼1-004" / "caliper_profile.csv"
# 现场提取包 inclination_profile.csv：设计一.1.4.1 井斜方位表 70 点（0-7660m，最大 8.35°@7660m）。
DEFAULT_INCLINATION_CSV = PROJECT_ROOT / "参考文档" / "现场资料提取" / "ht1_004_呼1-004" / "inclination_profile.csv"
# CBL 数字化产物：尾管评价段 5245-7581m 一界面中等及以上占比 0.3%（interpreted）。
DEFAULT_CBL_DIGITIZATION_CSV = PROJECT_ROOT / "参考文档" / "现场资料提取" / "ht1_004_呼1-004" / "cbl_digitization" / "cbl_digitization.csv"

# ===========================================================================
# 井段几何参数 [实测] 来源：设计文档 + 呼1-004井身结构.csv / 现场提取包 well_geometry.csv
# ===========================================================================
HT1_004_WELL_NAME = "呼1-004井（HT1-004）"
HT1_004_DRILLED_DEPTH_MD_M = 7660.0       # 实际完钻井深/TD
HT1_004_HANGER_MD_M = 5243.207            # 尾管悬挂器喇叭口（本体下深 5251.11m；CBL 头 5245.010-5251.110m）
HT1_004_TOP_MD_M = HT1_004_HANGER_MD_M    # 模型域从悬挂器开始（设计 4.4.2）
HT1_004_CASING_SHOE_MD_M = 5578.0         # 273.1mm套管鞋
HT1_004_UPPER_SECTION_BOTTOM_MD_M = 7376.656  # 168.3→139.7mm尾管变径/变扣（139.7 段自 7376.656m 起）；
                                              # 变径 BGT2*BG-FJU@7378.051m 双口径，任务锚点取 7376.66-7660。
                                              # LEGACY(2026-08-16 前): 7378.051（误用变径变扣深度）
HT1_004_LOWER_HOLE_TOP_MD_M = 7521.0      # 241.3→215.9mm井眼变径（钻头程序；电测井径表 7490m 以下 216mm，差 31m 已注明）
HT1_004_BOTTOM_MD_M = HT1_004_DRILLED_DEPTH_MD_M
HT1_004_SHOE_MD_M = HT1_004_DRILLED_DEPTH_MD_M

# --- 套管/管柱尺寸 [实测] ---
# 2026-08-29 语义统一：casing_id_mm 按 PACKAGE_REFERENCE 文档语义存"外层套管内径"——
# 273.05mm 技套（名义 273.1）真实 ID=245.37（设计电测井径表口径；计算值 245.42）；OD 公称 273.1 存档。
# 字段不被求解器消费（环空计算由 casing_inner_diameter 与 liner_od_profile/hole_diameter_profile 完成）。
# LEGACY(2026-08-29 前): 本字段传 OD 273.1（名义口径，命名历史遗留）。
HT1_004_CASING_OD_MM = 273.1
HT1_004_CASING_INNER_DIAMETER_MM = 245.37  # 设计电测井径表口径（计算值 245.42）
HT1_004_UPPER_LINER_OD_MM = 168.3
HT1_004_LOWER_LINER_OD_MM = 139.7
# 上段 168.3mm 尾管壁厚 14.7mm（现场设计 6.1/施工记录表：168.3-2*14.7=138.9mm）。
# LEGACY(2026-08-16 前): 15.88（误用 139.7mm 段壁厚，导致上段 ID 136.54mm 偏小）。
HT1_004_UPPER_LINER_WALL_THICKNESS_MM = 14.7
HT1_004_LOWER_LINER_WALL_THICKNESS_MM = 15.88
HT1_004_UPPER_LINER_ID_MM = HT1_004_UPPER_LINER_OD_MM - 2.0 * HT1_004_UPPER_LINER_WALL_THICKNESS_MM  # 138.9mm
HT1_004_LOWER_LINER_ID_MM = HT1_004_LOWER_LINER_OD_MM - 2.0 * HT1_004_LOWER_LINER_WALL_THICKNESS_MM  # 107.94mm
HT1_004_LINER_ID_MM = HT1_004_LOWER_LINER_ID_MM

# --- 平均居中度 [实测-模拟] 来源：设计文档第6.3节（model_assumption，无连续实测） ---
HT1_004_AVERAGE_STANDOFF = 0.83

# ===========================================================================
# CSV 剖面数据读取（现场提取包）
# ===========================================================================

def _read_caliper_rows(caliper_csv_path: Path) -> tuple[tuple[float, float], ...]:
    """读取现场提取包井径剖面 CSV（md_m / caliper_mm），按 md_m 升序返回。

    含套管内/重叠段环空行（caliper=技套内径 245.37mm，非裸眼井径）与裸眼电测行。
    """
    rows: list[tuple[float, float]] = []
    with caliper_csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            md = row.get("md_m")
            cal = row.get("caliper_mm")
            if md is None or cal is None:
                continue
            rows.append((float(md), float(cal)))
    if not rows:
        raise ValueError(f"井径 CSV 为空: {caliper_csv_path}")
    # 轻量 sanity check（2026-08-29 校准）：汇总统计行曾混入 caliper_profile.csv（本次已删 4 行），
    # 这里按"重复深度必错"直接拒绝，并对异常密集采样告警，防止汇总统计行再次混入。
    depths = [md for md, _ in rows]
    duplicates = sorted({d for d in depths if depths.count(d) > 1})
    if duplicates:
        raise ValueError(f"井径 CSV 存在重复深度（疑汇总统计行混入）: {duplicates}: {caliper_csv_path}")
    sorted_depths = sorted(depths)
    min_gap = min((b - a) for a, b in zip(sorted_depths, sorted_depths[1:]))
    if min_gap < 1.0:
        warnings.warn(
            f"井径 CSV 相邻深度最小间隔仅 {min_gap:.3f}m，疑含非测点行/汇总行，请人工复核: {caliper_csv_path}"
        )
    return tuple(sorted(rows))


def _read_inclination_rows(inclination_csv_path: Path) -> tuple[tuple[float, float], ...]:
    """读取现场提取包井斜剖面 CSV（md_m / inclination_deg），按 md_m 升序返回。"""
    rows: list[tuple[float, float]] = []
    with inclination_csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            md = row.get("md_m")
            inc = row.get("inclination_deg")
            if md is None or inc is None:
                continue
            rows.append((float(md), float(inc)))
    if not rows:
        raise ValueError(f"井斜 CSV 为空: {inclination_csv_path}")
    return tuple(sorted(rows))


# ===========================================================================
# 等效井径换算（保留供旧口径追溯；当前剖面直接用现场提取包裸眼实测 caliper_mm）
# ===========================================================================

def _equivalent_hole_diameter_mm(
    actual_hole_mm: float, actual_od_mm: float, reference_od_mm: float
) -> float:
    """在固定 reference_od 条件下，构造与原始环空面积等价的井眼直径。"""
    area_term = actual_hole_mm**2 - actual_od_mm**2 + reference_od_mm**2
    return math.sqrt(max(area_term, reference_od_mm**2))


# ===========================================================================
# 剖面构建（基于现场提取包 CSV）
# ===========================================================================

def _build_hole_profile(caliper_rows: tuple[tuple[float, float], ...]) -> tuple[tuple[float, float], ...]:
    """从现场提取包井径行构建井径剖面（裸眼实测 caliper_mm）。

    5243.207-5578m 为 273.1mm 套管内重叠段，取技套内径 245.37mm 作为等效井径
    （与 168.3mm 尾管环空间隙一致，由 liner_od_profile 表达）。
    """
    points: list[tuple[float, float]] = [(md, cal) for md, cal in caliper_rows if md >= HT1_004_HANGER_MD_M - 0.1]
    if not points or points[0][0] > HT1_004_HANGER_MD_M:
        points.insert(0, (HT1_004_HANGER_MD_M, HT1_004_CASING_INNER_DIAMETER_MM))
    if points[-1][0] < HT1_004_BOTTOM_MD_M:
        points.append((HT1_004_BOTTOM_MD_M, points[-1][1]))
    return tuple(points)


def _build_inclination_profile(incl_rows: tuple[tuple[float, float], ...]) -> tuple[tuple[float, float], ...]:
    """从现场提取包井斜行构建井斜剖面（悬挂器-5600m 之间取首测点 0.75° 外推）。"""
    points: list[tuple[float, float]] = [(md, inc) for md, inc in incl_rows if md >= HT1_004_HANGER_MD_M - 0.1]
    if not points or points[0][0] > HT1_004_HANGER_MD_M:
        points.insert(0, (HT1_004_HANGER_MD_M, 0.75))
    if points[-1][0] < HT1_004_BOTTOM_MD_M:
        points.append((HT1_004_BOTTOM_MD_M, points[-1][1]))
    return tuple(points)


def _build_standoff_profile() -> tuple[tuple[float, float], ...]:
    """构建居中度剖面，全井段固定83%（设计文档6.3节模拟值，model_assumption）。"""
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
    段3: 5243.207~7376.656m  → 138.9mm (168.3mm尾管, 壁厚14.7mm)
    段4: 7376.656~7660m      → 107.94mm (139.7mm尾管, 壁厚15.88mm)
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

    段1: 5243.207~7376.656m → 168.3mm (上段尾管)
    段2: 7376.656~7660m     → 139.7mm (下段尾管)
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


# 4 段管柱 ID 累加（上段 ID 已按 14.7mm 壁厚修正为 138.9mm，结果 99.36m³，名义内径算至 7660）。
# 三口径并存（2026-08-29 校准）：99.36（名义内径算至 7660）≈ 作业史 99.31；
# 设计/碰压口径 97.1（149.2 段实测 12.9 L/m、算至球座 7560.286、含压塞 2）。loader 取 4 段管柱 ID 累加值。
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
# 流体密度 [field_measured] 来源：现场提取包 fluid_properties.csv（密度与现场全部一致）
# ===========================================================================
HT1_004_MUD_DENSITY_KG_M3 = 1900.0            # 钻井液 1.90（固井前，现场）
HT1_004_LEAD_MUD_DENSITY_KG_M3 = 1750.0       # 先导浆 1.75（现场）
HT1_004_SPACER1_DENSITY_KG_M3 = 1950.0        # 隔离液1 1.95（现场）
HT1_004_SPACER2_DENSITY_KG_M3 = 1750.0        # 隔离液2 1.75（现场）
HT1_004_LEAD_DENSITY_KG_M3 = 1930.0           # 领浆 1.93（现场）
HT1_004_TAIL_DENSITY_KG_M3 = 1900.0           # 尾浆 1.90（现场）
HT1_004_PLUG_DENSITY_KG_M3 = 1700.0           # 压塞液 1.70（实际后置液，现场）
HT1_004_DISPLACEMENT_DENSITY_KG_M3 = 1900.0   # 替钻井液 1.90（现场）
HT1_004_BUFFER_DENSITY_KG_M3 = 1900.0         # 保护液 1.90（现场）
HT1_004_BASE_FLUID_DENSITY_KG_M3 = 1020.0     # 基液 1.02（设计）
HT1_004_WELL_MUD_DENSITY_KG_M3 = 1900.0       # 井浆(替浆段) 1.90（现场）

# ===========================================================================
# 流体流变参数 [optimized_input] 来源：优化参数.docx（2026-06-11）+ MATLAB基础参数代码.docx
# —— 非现场实测，仅作优化/应用案例输入；化验报告幂律实测见 _POWER_LAW 常量（实际版用）
# PV=mPa·s/1000→Pa·s, YP=Pa
# ===========================================================================
HT1_004_MUD_PV_PA_S = 0.053             # 井浆 PV=53mPa·s (MATLAB miu0)
HT1_004_MUD_YP_PA = 8.5                 # YP=8.5Pa (MATLAB tau0)

HT1_004_LEAD_MUD_PV_PA_S = 0.058        # 先导浆 PV=58mPa·s (MATLAB miu1) [optimized]
HT1_004_LEAD_MUD_YP_PA = 9.8            # YP=9.8Pa (MATLAB tau1) [optimized]

HT1_004_SPACER1_PV_PA_S = 0.058         # 隔离液1 PV=58mPa·s (MATLAB miu2) [optimized]
HT1_004_SPACER1_YP_PA = 9.8             # YP=9.8Pa (MATLAB tau2) [optimized]

HT1_004_SPACER2_PV_PA_S = 0.065         # 隔离液2 PV=65mPa·s (MATLAB miu3) [optimized]
HT1_004_SPACER2_YP_PA = 10.0            # YP=10Pa (MATLAB tau3) [optimized]

HT1_004_LEAD_PV_PA_S = 0.170            # 领浆 PV=170mPa·s [optimized；化验幂律 n=0.853/K=0.746]
HT1_004_LEAD_YP_PA = 13.0               # YP=13Pa [optimized]

HT1_004_TAIL_PV_PA_S = 0.180            # 尾浆 PV=180mPa·s (MATLAB miu5) [optimized；化验幂律 n=0.869/K=0.669]
HT1_004_TAIL_YP_PA = 14.0               # YP=14Pa (MATLAB tau5) [optimized]

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

# 化验报告幂律实测（field_measured，实际版 load_ht1_004_tailpipe_actual 用）
# 化验报告（132C↓93C）：钻井液 123/70/51/31/7/5 → n=0.635/K=0.738；
# 领浆 → n=0.853/K=0.746；尾浆 → n=0.869/K=0.669；隔离液1/2 同流变 → n=0.611/K=0.736。
HT1_004_MUD_POWER_LAW_N = 0.635
HT1_004_MUD_CONSISTENCY_K = 0.738
HT1_004_LEAD_POWER_LAW_N = 0.853
HT1_004_LEAD_CONSISTENCY_K = 0.746
HT1_004_TAIL_POWER_LAW_N = 0.869
HT1_004_TAIL_CONSISTENCY_K = 0.669
HT1_004_SPACER_POWER_LAW_N = 0.611
HT1_004_SPACER_CONSISTENCY_K = 0.736

# 合成 FLUSHER（冲洗液）参数 [model_assumption]：本井设计为先导浆+隔离液+领/尾浆，
# 现场无 cement 前独立冲洗液；保留合成 FLUSHER 仅用于验证 mud-spacer-flusher-cement 序列可表达
# （include_wash_spacer=True 时才注入，默认不注入行为保持）。
HT1_004_FLUSHER_DENSITY_KG_M3 = 1880.0
HT1_004_FLUSHER_PV_PA_S = 0.050
HT1_004_FLUSHER_YP_PA = 9.0
HT1_004_FLUSHER_VOLUME_M3 = 5.0
HT1_004_FLUSHER_RATE_M3_MIN = 1.2

# ===========================================================================
# 泵注体积 [field_measured，与现场一致] / 排量 [optimized_input] 来源：优化参数.docx（2026-06-11）
# 默认输入为优化参数化：替浆分5级降排量 1.15→1.05→0.95→0.85→0.75m³/min，
# 分段体积14+10+10+10+7.1=51.1m³，末段7.1+2.3补足=9.4m³；总替浆(含前置)97.4m³（非现场 97.1）。
# 实际施工排量见 HT1_004_ACTUAL_* 常量与 load_ht1_004_tailpipe_actual。
# ===========================================================================
HT1_004_LEAD_MUD_VOLUME_M3 = 25.0       # 先导液 25m³（现场）
HT1_004_SPACER1_VOLUME_M3 = 16.0        # 隔离液1 16m³（现场）
HT1_004_SPACER2_VOLUME_M3 = 10.0        # 隔离液2 10m³（现场）
HT1_004_LEAD_VOLUME_M3 = 48.0           # 领浆 48m³（现场）
HT1_004_TAIL_VOLUME_M3 = 28.0           # 尾浆 28m³（现场）
HT1_004_PLUG_VOLUME_M3 = 2.0            # 压塞液 2m³（现场）
HT1_004_DISPLACEMENT_FAST_VOLUME_M3 = 29.0   # 钻井液 29m³（现场）
HT1_004_BUFFER_VOLUME_M3 = 14.0         # 保护液 14m³（现场）
HT1_004_BASE_FLUID_VOLUME_M3 = 1.0      # 基液 1m³（设计）
HT1_004_WELL_MUD_1_VOLUME_M3 = 14.0     # 替浆段1 14m³ @1.15m³/min [optimized]
HT1_004_WELL_MUD_2_VOLUME_M3 = 10.0     # 替浆段2 10m³ @1.05m³/min [optimized]
HT1_004_WELL_MUD_3_VOLUME_M3 = 10.0     # 替浆段3 10m³ @0.95m³/min [optimized]
HT1_004_WELL_MUD_4_VOLUME_M3 = 10.0     # 替浆段4 10m³ @0.85m³/min [optimized]
HT1_004_WELL_MUD_5_VOLUME_M3 = 9.4      # 替浆段5 7.1+2.3补足=9.4m³ @0.75m³/min [optimized；非现场]

HT1_004_LEAD_MUD_RATE_M3_MIN = 1.4      # 先导液 1.4m³/min [optimized；设计1.4/实际1.2]
HT1_004_SPACER_RATE_M3_MIN = 1.2        # 隔离液 1.2m³/min [optimized；实际1.0]
HT1_004_CEMENT_RATE_M3_MIN = 1.2        # 领浆 1.2m³/min（实际/设计）
HT1_004_TAIL_RATE_M3_MIN = 1.25         # 尾浆 1.25m³/min [optimized；实际1.0/设计1.2]
HT1_004_PLUG_RATE_M3_MIN = 1.2          # 压塞液 1.2m³/min（实际 7a 亦 1.2）
HT1_004_DISPLACEMENT_FAST_RATE_M3_MIN = 1.5  # 钻井液 1.5m³/min（实际 7b 亦 1.5）
HT1_004_BUFFER_RATE_M3_MIN = 1.4        # 保护液 1.4m³/min [optimized；作业史1.5]
HT1_004_BASE_FLUID_RATE_M3_MIN = 1.4    # 基液 1.4m³/min [optimized]
HT1_004_WELL_MUD_1_RATE_M3_MIN = 1.15   # 替浆段1 1.15m³/min [optimized；设计1.2]
HT1_004_WELL_MUD_2_RATE_M3_MIN = 1.05   # 替浆段2 1.05m³/min [optimized；设计1.0]
HT1_004_WELL_MUD_3_RATE_M3_MIN = 0.95   # 替浆段3 0.95m³/min [optimized；设计0.9]
HT1_004_WELL_MUD_4_RATE_M3_MIN = 0.85   # 替浆段4 0.85m³/min [optimized；设计0.8]
HT1_004_WELL_MUD_5_RATE_M3_MIN = 0.75   # 替浆段5 0.75m³/min [optimized；设计0.7]

# ===========================================================================
# 实际施工版参数（load_ht1_004_tailpipe_actual，施工记录表/作业史 field_measured）
# 实际替浆 97.1m³ 分 5 段（作业史 7a-7e）：压塞2 + 钻井液29 + 保护14 + 井浆30 + 井浆末段22。
# ===========================================================================
HT1_004_ACTUAL_LEAD_MUD_RATE_M3_MIN = 1.2   # 先导浆 实际 1.2（设计 1.4，作业史 1.5）
HT1_004_ACTUAL_SPACER_RATE_M3_MIN = 1.0     # 隔离液 实际 1.0（设计 1.2）
HT1_004_ACTUAL_TAIL_RATE_M3_MIN = 1.0       # 尾浆 实际 1.0（施工记录表；作业史 1.2）
HT1_004_ACTUAL_PLUG_RATE_M3_MIN = 1.2       # 压塞液 实际 1.2（作业史 7a，泵压13MPa）
HT1_004_ACTUAL_DISPLACEMENT_FAST_RATE_M3_MIN = 1.5  # 钻井液 实际 1.5（作业史 7b，泵压16MPa）
HT1_004_ACTUAL_BUFFER_RATE_M3_MIN = 1.5     # 保护液 实际 1.5（作业史 7c，泵压16MPa）
HT1_004_ACTUAL_WELL_MUD_SEQ1_VOLUME_M3 = 30.0  # 井浆 30m³（作业史 7d，排量1.2-1.0，泵压16MPa）
HT1_004_ACTUAL_WELL_MUD_SEQ1_RATE_M3_MIN = 1.2
HT1_004_ACTUAL_WELL_MUD_SEQ2_VOLUME_M3 = 22.0  # 井浆末段 22m³（作业史 7e，排量1.0-0.8-0.6，泵压11-9MPa）
HT1_004_ACTUAL_WELL_MUD_SEQ2_RATE_M3_MIN = 1.0


# ===========================================================================
# 辅助函数
# ===========================================================================

def _depth_points(values: tuple[tuple[float, float], ...]) -> tuple[DepthValuePoint, ...]:
    """把测深-数值元组转换为 WellSpec 使用的剖面点。"""
    return tuple(DepthValuePoint(depth_md_m=depth, value=value) for depth, value in values)


def _build_well_spec(
    resolved_reference_root: Path,
    caliper_csv_path: Path,
    inclination_csv_path: Path,
) -> WellSpec:
    """构建呼1-004 井筒规格（几何/剖面/评价窗，默认优化版与实际版共用）。"""
    caliper_rows = _read_caliper_rows(caliper_csv_path)
    incl_rows = _read_inclination_rows(inclination_csv_path)
    hole_profile = _build_hole_profile(caliper_rows)
    inc_profile = _build_inclination_profile(incl_rows)
    standoff_profile = _build_standoff_profile()
    liner_od_profile = _build_liner_od_profile()
    pipe_id_profile = _build_pipe_id_profile()

    return WellSpec(
        well_name=HT1_004_WELL_NAME,
        top_md_m=HT1_004_TOP_MD_M,
        bottom_md_m=HT1_004_BOTTOM_MD_M,
        shoe_md_m=HT1_004_SHOE_MD_M,
        hanger_md_m=HT1_004_HANGER_MD_M,
        casing_id_mm=HT1_004_CASING_INNER_DIAMETER_MM,  # 2026-08-29 语义统一：存真实内径 245.37（原传 OD 273.1）
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
            # 验证窗口口径（汇总 §10.2）：CBL 评价窗取尾管段 5245-7581m；目的井段 7495-7526+7531-7550 两段；
            # 油气水层段 7482-7560m 为官方红线关注段。
            EvaluationWindow(name="CBL评价井段(尾管段)", top_md_m=5245.0, bottom_md_m=7581.0, window_type="cbl"),
            EvaluationWindow(name="目的井段1(地层目标)", top_md_m=7495.0, bottom_md_m=7526.0, window_type="formation_target"),
            EvaluationWindow(name="目的井段2(地层目标)", top_md_m=7531.0, bottom_md_m=7550.0, window_type="formation_target"),
            EvaluationWindow(name="油气水层段(红线)", top_md_m=7482.0, bottom_md_m=7560.0, window_type="oil_gas_show"),
        ),
        reference_root=resolved_reference_root,
        notes=(
            "呼1-004井（HT1-004）168.3+139.7mm双径尾管控压固井，井深7660m。",
            "井径/井斜剖面来自现场提取包 caliper_profile.csv / inclination_profile.csv（设计一.1.4.2/1.4.1）。",
            "悬挂器取喇叭口 5243.207m（本体下深 5251.11m；CBL 头 5245.010-5251.110m，差 1.8m 已注明）；"
            "变径取 7376.656m（变径变扣 7378.051m 双口径）。",
            "pipe_id_profile按4段管柱内径：129.9/107.7/138.9/107.94mm；上段壁厚 14.7mm（LEGACY 误用 15.88）。",
            "默认输入为'优化参数化'：流体流变（领浆PV170/YP13、尾浆180/14等）与替浆排量（1.15-0.75、末段9.4补足）"
            "来自优化参数.docx（2026-06-11），非现场实测；密度与泵注体积为现场值。",
            "居中度采用设计文档6.3节模拟值(平均83%)固定剖面（model_assumption）。",
            "2026-08-29 校准补记：鞋口滞后三口径 99.36（名义内径算至 7660）≈ 作业史 99.31 vs 设计/碰压口径 97.1"
            "（149.2 段实测 12.9 L/m、算至球座 7560.286、含压塞 2），loader 取管柱 ID 累加值；"
            "替浆介质密度双口径 1.90（设计/总结）vs 1.91（记录表）；"
            "'灰量 161t'为图头浆体质量（领 100+尾 61）非干灰（干灰 72.7t）；"
            "技套 273.05（图头/名义）vs 273.1 双口径；"
            "MATLAB 基础参数代码.docx 的 168.3-2×15.88 为上游错误源，勿回滚 2026-08-16 壁厚修正。",
        ),
    )


def _build_fluids() -> tuple[FluidSpec, ...]:
    """呼1-004 流体清单（默认优化版；密度 field_measured、流变 optimized_input）。"""
    return (
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
        # 合成 FLUSHER（model_assumption，现场无此流体）：保留定义以维持 mud-spacer-flusher-cement
        # 序列可表达及测试兼容；include_wash_spacer=False 默认不进 schedule。
        FluidSpec("冲洗液（FLUSHER）", FluidRole.FLUSHER, HT1_004_FLUSHER_DENSITY_KG_M3,
                  RheologyModel.BINGHAM, HT1_004_FLUSHER_PV_PA_S, HT1_004_FLUSHER_YP_PA),
    )


def _build_fluids_actual() -> tuple[FluidSpec, ...]:
    """呼1-004 流体清单（实际施工版）：领浆/尾浆改用化验报告幂律实测（field_measured）。

    其余流体与默认版一致（密度现场、流变 optimized/proxy）。
    """
    return (
        FluidSpec("钻井液", FluidRole.MUD, HT1_004_MUD_DENSITY_KG_M3,
                  RheologyModel.POWER_LAW, power_law_n=HT1_004_MUD_POWER_LAW_N,
                  consistency_k=HT1_004_MUD_CONSISTENCY_K),
        FluidSpec("保护液", FluidRole.DISPLACEMENT, HT1_004_BUFFER_DENSITY_KG_M3,
                  RheologyModel.BINGHAM, HT1_004_BUFFER_PV_PA_S, HT1_004_BUFFER_YP_PA),
        FluidSpec("先导浆", FluidRole.WASH, HT1_004_LEAD_MUD_DENSITY_KG_M3,
                  RheologyModel.BINGHAM, HT1_004_LEAD_MUD_PV_PA_S, HT1_004_LEAD_MUD_YP_PA),
        FluidSpec("隔离液1", FluidRole.SPACER, HT1_004_SPACER1_DENSITY_KG_M3,
                  RheologyModel.POWER_LAW, power_law_n=HT1_004_SPACER_POWER_LAW_N,
                  consistency_k=HT1_004_SPACER_CONSISTENCY_K),
        FluidSpec("隔离液2", FluidRole.SPACER, HT1_004_SPACER2_DENSITY_KG_M3,
                  RheologyModel.POWER_LAW, power_law_n=HT1_004_SPACER_POWER_LAW_N,
                  consistency_k=HT1_004_SPACER_CONSISTENCY_K),
        FluidSpec("领浆", FluidRole.LEAD, HT1_004_LEAD_DENSITY_KG_M3,
                  RheologyModel.POWER_LAW, power_law_n=HT1_004_LEAD_POWER_LAW_N,
                  consistency_k=HT1_004_LEAD_CONSISTENCY_K),
        FluidSpec("尾浆", FluidRole.TAIL, HT1_004_TAIL_DENSITY_KG_M3,
                  RheologyModel.POWER_LAW, power_law_n=HT1_004_TAIL_POWER_LAW_N,
                  consistency_k=HT1_004_TAIL_CONSISTENCY_K),
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


# ===========================================================================
# 主加载函数
# ===========================================================================

def load_ht1_004_tailpipe(
    *,
    reference_root: Path | None = None,
    include_wash_spacer: bool = False,
    caliper_csv_path: Path | None = None,
    inclination_csv_path: Path | None = None,
) -> tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载呼1-004井（HT1-004）168.3+139.7mm双径尾管段标准模型输入（**优化参数化**）。

    ⚠️ 数据标签声明：本函数为**优化参数化输入，非现场实测**——流体流变（领浆PV170/YP13、
    尾浆180/14、先导58/9.8、隔离液58/9.8、65/10）与替浆排量（1.15/1.05/0.95/0.85/0.75、末段
    9.4=7.1+2.3 补足、总替浆 97.4m³）来自 优化参数.docx（2026-06-11），不能写成现场实测验证；
    论文使用应定位为"优化/应用案例"。密度与泵注体积为现场值（field_measured）。

    实际施工版（施工记录表/作业史，替浆 97.1m³ 分 5 段、排量 1.2/1.0）见
    load_ht1_004_tailpipe_actual。

    Args:
        reference_root: 可选参考资料根目录。
        include_wash_spacer: 是否注入合成 FLUSHER（冲洗液）步骤以验证
            mud-spacer-flusher-cement 序列可表达。默认为 False（严格现场模式，
            现场无 FLUSHER 流体，model_assumption）。
        caliper_csv_path: 可选现场提取包井径 CSV 路径。
        inclination_csv_path: 可选现场提取包井斜 CSV 路径。

    Returns:
        (well_spec, fluids, schedule, validation_data)
    """

    resolved_reference_root = reference_root or DEFAULT_REFERENCE_ROOT
    resolved_caliper_csv_path = caliper_csv_path or DEFAULT_CALIPER_CSV
    resolved_inclination_csv_path = inclination_csv_path or DEFAULT_INCLINATION_CSV

    well_spec = _build_well_spec(resolved_reference_root, resolved_caliper_csv_path, resolved_inclination_csv_path)
    fluids = _build_fluids()

    # 合成 FLUSHER 步骤：仅在 include_wash_spacer=True 时注入（model_assumption）。
    flusher_step: tuple[PumpingScheduleStep, ...] = ()
    if include_wash_spacer:
        flusher_step = (
            PumpingScheduleStep("注入冲洗液（FLUSHER）", "冲洗液（FLUSHER）",
                                HT1_004_FLUSHER_VOLUME_M3, HT1_004_FLUSHER_RATE_M3_MIN,
                                remarks="合成冲洗液（FLUSHER）5m³@1.2m³/min，密度1.88g/cm³（model_assumption，验证序列可表达）。"),
        )

    # 默认优化版泵注程序（优化参数.docx 2026-06-11）
    schedule = PumpingSchedule(
        steps=(
            PumpingScheduleStep("注入先导浆", "先导浆",
                                HT1_004_LEAD_MUD_VOLUME_M3, HT1_004_LEAD_MUD_RATE_M3_MIN,
                                remarks="先导浆 25m³@1.4m³/min，密度1.75g/cm³（排量optimized；实际1.2）。"),
            PumpingScheduleStep("注入隔离液1", "隔离液1",
                                HT1_004_SPACER1_VOLUME_M3, HT1_004_SPACER_RATE_M3_MIN,
                                remarks="隔离液1 16m³@1.2m³/min，密度1.95g/cm³（排量optimized；实际1.0）。"),
            PumpingScheduleStep("注入隔离液2", "隔离液2",
                                HT1_004_SPACER2_VOLUME_M3, HT1_004_SPACER_RATE_M3_MIN,
                                remarks="隔离液2 10m³@1.2m³/min，密度1.75g/cm³（排量optimized；实际1.0）。"),
            *flusher_step,
            PumpingScheduleStep("注入领浆", "领浆",
                                HT1_004_LEAD_VOLUME_M3, HT1_004_CEMENT_RATE_M3_MIN,
                                remarks="领浆 48m³@1.2m³/min，密度1.93g/cm³，PV=170mPa·s/YP=13Pa（流变optimized）。"),
            PumpingScheduleStep("注入尾浆", "尾浆",
                                HT1_004_TAIL_VOLUME_M3, HT1_004_TAIL_RATE_M3_MIN,
                                remarks="尾浆 28m³@1.25m³/min，密度1.90g/cm³（排量optimized；实际1.0）。"),
            PumpingScheduleStep("替入压塞液", "压塞液",
                                HT1_004_PLUG_VOLUME_M3, HT1_004_PLUG_RATE_M3_MIN,
                                remarks="压塞液 2m³@1.2m³/min，密度1.70g/cm³（实际后置液）。"),
            PumpingScheduleStep("替钻井液(快)", "替钻井液",
                                HT1_004_DISPLACEMENT_FAST_VOLUME_M3,
                                HT1_004_DISPLACEMENT_FAST_RATE_M3_MIN,
                                remarks="钻井液 29m³@1.5m³/min，密度1.90g/cm³（实际 7b 亦 1.5）。"),
            PumpingScheduleStep("替保护液", "保护液",
                                HT1_004_BUFFER_VOLUME_M3, HT1_004_BUFFER_RATE_M3_MIN,
                                remarks="保护液 14m³@1.4m³/min，密度1.90g/cm³（排量optimized；作业史1.5）。"),
            PumpingScheduleStep("替基液", "基液",
                                HT1_004_BASE_FLUID_VOLUME_M3, HT1_004_BASE_FLUID_RATE_M3_MIN,
                                remarks="基液 1m³@1.4m³/min，密度1.02g/cm³（排量optimized）。"),
            PumpingScheduleStep("井浆替入1", "井浆",
                                HT1_004_WELL_MUD_1_VOLUME_M3,
                                HT1_004_WELL_MUD_1_RATE_M3_MIN,
                                remarks="替浆段1 14m³@1.15m³/min（排量optimized；设计1.2）。"),
            PumpingScheduleStep("井浆替入2", "井浆",
                                HT1_004_WELL_MUD_2_VOLUME_M3,
                                HT1_004_WELL_MUD_2_RATE_M3_MIN,
                                remarks="替浆段2 10m³@1.05m³/min（排量optimized；设计1.0）。"),
            PumpingScheduleStep("井浆替入3", "井浆",
                                HT1_004_WELL_MUD_3_VOLUME_M3,
                                HT1_004_WELL_MUD_3_RATE_M3_MIN,
                                remarks="替浆段3 10m³@0.95m³/min（排量optimized；设计0.9）。"),
            PumpingScheduleStep("井浆替入4", "井浆",
                                HT1_004_WELL_MUD_4_VOLUME_M3,
                                HT1_004_WELL_MUD_4_RATE_M3_MIN,
                                remarks="替浆段4 10m³@0.85m³/min（排量optimized；设计0.8）。"),
            PumpingScheduleStep("井浆替入5", "井浆",
                                HT1_004_WELL_MUD_5_VOLUME_M3,
                                HT1_004_WELL_MUD_5_RATE_M3_MIN,
                                remarks="替浆段5 9.4m³@0.75m³/min（optimized：7.1+2.3m³补足尾浆全入环空，非现场）。"),
        ),
        notes=(
            "⚠️ 默认输入为优化参数化（优化参数.docx 2026-06-11），非现场实测；"
            "施工顺序：先导浆→隔离液1→隔离液2"
            + ("→冲洗液（FLUSHER）" if include_wash_spacer else "")
            + "→领浆→尾浆→压塞液→钻井液→保护液→基液→5级降排量替浆。",
            "替浆总量(不含压塞液): 29+14+1+14+10+10+10+9.4=97.4m³（optimized，含末段补足2.3m³；现场97.1）。",
            "实际施工版排量/替浆分段见 load_ht1_004_tailpipe_actual（作业史：先导1.2/隔离1.0/尾浆1.0/替浆97.1）。",
            "合成FLUSHER步骤仅在 include_wash_spacer=True 时注入（model_assumption）。",
        ),
    )

    validation_data = ValidationData(
        cbl_summary_path=DEFAULT_CBL_DIGITIZATION_CSV,
        cbl_pass_rate=0.003,
        job_report_path=resolved_reference_root / "呼1-004井身结构.csv",
        notes=(
            "cbl_pass_rate=0.003 为尾管评价段 5245-7581m 一界面中等及以上占比（数字化 0.3%，interpreted，来源 cbl_digitization.csv）。",
            "官方 CBL 结论：不合格，一、二界面综合合格率 29.99%（全井 11-7581m，未达 70% 红线；CPLog/LEAD5.0 2026-07-23 解释）。",
            "RCD 水泥充填成像判'合格'与 CBL 声学胶结判'不合格'不矛盾（不同物理量）；模型验证以 CBL 一界面为准。",
            "默认输入为优化参数化（流变/排量来自优化参数.docx），非现场实测；密度与泵注体积为现场值。",
        ),
    )
    return well_spec, fluids, schedule, validation_data


def load_ht1_004_tailpipe_actual(
    *,
    reference_root: Path | None = None,
    include_wash_spacer: bool = False,
    caliper_csv_path: Path | None = None,
    inclination_csv_path: Path | None = None,
) -> tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载呼1-004井（HT1-004）168.3+139.7mm双径尾管段"实际施工版"模型输入。

    与默认优化版（load_ht1_004_tailpipe）的差异：
    - 排量改用施工记录表/作业史实际值：先导浆 1.2、隔离液 1.0、尾浆 1.0（作业史记 1.2）；
    - 替浆按作业史 7a-7e 分 5 段（压塞2 + 钻井液29 + 保护14 + 井浆30 + 井浆末段22），总 97.1m³；
    - 领浆/尾浆流变改用化验报告幂律实测（领浆 n=0.853/K=0.746、尾浆 n=0.869/K=0.669）。

    Returns:
        (well_spec, fluids, actual_schedule, validation_data)
    """
    resolved_reference_root = reference_root or DEFAULT_REFERENCE_ROOT
    resolved_caliper_csv_path = caliper_csv_path or DEFAULT_CALIPER_CSV
    resolved_inclination_csv_path = inclination_csv_path or DEFAULT_INCLINATION_CSV

    well_spec = _build_well_spec(resolved_reference_root, resolved_caliper_csv_path, resolved_inclination_csv_path)
    fluids = _build_fluids_actual()

    flusher_step: tuple[PumpingScheduleStep, ...] = ()
    if include_wash_spacer:
        flusher_step = (
            PumpingScheduleStep("注入冲洗液（FLUSHER）", "冲洗液（FLUSHER）",
                                HT1_004_FLUSHER_VOLUME_M3, HT1_004_FLUSHER_RATE_M3_MIN,
                                remarks="合成冲洗液（FLUSHER）5m³@1.2m³/min，密度1.88g/cm³（model_assumption）。"),
        )

    # 实际施工泵注程序（施工记录表/作业史）
    schedule = PumpingSchedule(
        steps=(
            PumpingScheduleStep("注入先导浆(实际)", "先导浆",
                                HT1_004_LEAD_MUD_VOLUME_M3, HT1_004_ACTUAL_LEAD_MUD_RATE_M3_MIN,
                                remarks="实际: 先导浆 25m³@1.2m³/min，密度1.75g/cm³，泵压21-24.6MPa。"),
            PumpingScheduleStep("注入隔离液1(实际)", "隔离液1",
                                HT1_004_SPACER1_VOLUME_M3, HT1_004_ACTUAL_SPACER_RATE_M3_MIN,
                                remarks="实际: 隔离液1 16m³@1.0m³/min，密度1.95g/cm³，泵压19MPa。"),
            PumpingScheduleStep("注入隔离液2(实际)", "隔离液2",
                                HT1_004_SPACER2_VOLUME_M3, HT1_004_ACTUAL_SPACER_RATE_M3_MIN,
                                remarks="实际: 隔离液2 10m³@1.0m³/min，密度1.75g/cm³，泵压18-19MPa。"),
            *flusher_step,
            PumpingScheduleStep("注入领浆(实际)", "领浆",
                                HT1_004_LEAD_VOLUME_M3, HT1_004_CEMENT_RATE_M3_MIN,
                                remarks="实际: 领浆 48m³@1.2m³/min，密度1.93g/cm³，泵压21-22MPa；流变按化验 n=0.853/K=0.746。"),
            PumpingScheduleStep("注入尾浆(实际)", "尾浆",
                                HT1_004_TAIL_VOLUME_M3, HT1_004_ACTUAL_TAIL_RATE_M3_MIN,
                                remarks="实际: 尾浆 28m³@1.0m³/min，密度1.90g/cm³，泵压22MPa（作业史记1.2）；流变按化验 n=0.869/K=0.669。"),
            PumpingScheduleStep("注入压塞液(实际)", "压塞液",
                                HT1_004_PLUG_VOLUME_M3, HT1_004_ACTUAL_PLUG_RATE_M3_MIN,
                                remarks="实际: 压塞液 2m³@1.2m³/min，密度1.70g/cm³，泵压13MPa（作业史 7a；记录表 15MPa 变体，2026-08-29 补记）。"),
            PumpingScheduleStep("替钻井液(实际)", "替钻井液",
                                HT1_004_DISPLACEMENT_FAST_VOLUME_M3,
                                HT1_004_ACTUAL_DISPLACEMENT_FAST_RATE_M3_MIN,
                                remarks="实际: 钻井液 29m³@1.5m³/min，密度1.90g/cm³，泵压16MPa（作业史 7b）。"),
            PumpingScheduleStep("替保护液(实际)", "保护液",
                                HT1_004_BUFFER_VOLUME_M3, HT1_004_ACTUAL_BUFFER_RATE_M3_MIN,
                                remarks="实际: 保护液 14m³@1.5m³/min，密度1.90g/cm³，泵压16MPa（作业史 7c）。"),
            PumpingScheduleStep("替井浆(实际)", "井浆",
                                HT1_004_ACTUAL_WELL_MUD_SEQ1_VOLUME_M3,
                                HT1_004_ACTUAL_WELL_MUD_SEQ1_RATE_M3_MIN,
                                remarks="实际: 井浆 30m³@1.2-1.0m³/min，密度1.90g/cm³，泵压16MPa（作业史 7d）。"),
            PumpingScheduleStep("替井浆末段(实际)", "井浆",
                                HT1_004_ACTUAL_WELL_MUD_SEQ2_VOLUME_M3,
                                HT1_004_ACTUAL_WELL_MUD_SEQ2_RATE_M3_MIN,
                                remarks="实际: 井浆末段 22m³@1.0-0.8-0.6m³/min，密度1.90g/cm³，泵压11-9MPa；替浆到量碰压9-15MPa（作业史 7e）。"),
        ),
        notes=(
            "实际施工顺序（施工记录表/作业史）：先导浆→隔离液1→隔离液2→领浆→尾浆→压塞液→钻井液→保护液→井浆→井浆末段。",
            "实际替浆分 5 段（作业史 7a-7e）：压塞液2 + 钻井液29 + 保护液14 + 井浆30 + 井浆末段22 = 97.0m³；"
            "固井总结/设计 7.2 口径 97.1m³（含基液1、钻井液 51.1），差 0.1m³ 为文档间取整差异，未编造。",
            "实际排量：先导浆 1.2、隔离液 1.0、尾浆 1.0（施工记录表；作业史记 1.2）、替浆 1.5→1.2→1.0 递减。",
            "合成FLUSHER步骤仅在 include_wash_spacer=True 时注入（model_assumption）。",
        ),
    )

    validation_data = ValidationData(
        cbl_summary_path=DEFAULT_CBL_DIGITIZATION_CSV,
        cbl_pass_rate=0.003,
        job_report_path=resolved_reference_root / "呼1-004井身结构.csv",
        notes=(
            "cbl_pass_rate=0.003 为尾管评价段 5245-7581m 一界面中等及以上占比（数字化 0.3%，interpreted）。",
            "官方 CBL 结论：不合格，一、二界面综合合格率 29.99%（全井 11-7581m，未达 70% 红线）。",
            "实际施工版与默认优化版共用几何；差异在泵注排量/替浆分段与领尾浆流变（化验幂律）。",
        ),
    )
    return well_spec, fluids, schedule, validation_data
