"""
呼1-003井（HT1-003）168.3+139.7mm双径尾管段标准数据加载器

本模块把呼1-003井（HT1-003）现场数据包中的 168.3+139.7mm 双径尾管段
资料整理为 cemdisp 标准输入结构。

当前求解目标为下部 139.7mm 尾管段，因此上部 168.3mm 重叠井段用等效井眼
直径保面积近似；鞋口滞后体积按 HT1-003 明确内径分段累加。

2026-08-16 按现场提取包核对（00_执行记录/ht1_003_004_loader核对_2026-08-16.md）更新：
- 井径/井斜剖面改用现场提取包 caliper_profile.csv / inclination_profile.csv（与 hu102/hu103 一致）；
  原 `参考文档/呼1-003/新/呼1-003井身结构.csv` 环空体积反算等效井径以 LEGACY 保留
- 尾浆密度修正 2.05→1.95（现场化验实测）；先导浆密度 1.85→1.75、体积 35→28；
  隔离液拆为隔离液1(16m³/ρ2.05)+隔离液2(10m³/ρ1.95) 两段（删单一隔离液 35m³/ρ2.00）
- 领浆 38→39、尾浆 28.5→28、井浆(替浆段)密度 1.85→1.95
- 扶正器 21+78=99 → 现场 75+22=97（按分段，设计 6.6.2/作业史）
- 鞋口滞后体积修复 7868 常数 bug（复制 hu101 TD 的笔误）→ 4 段管柱 ID 累加
- 流变改用化验报告幂律实测（领浆 n=0.597/K=1.622、尾浆 n=0.585/K=1.673、
  隔离液 n=0.668/K=1.245、钻井液 n=0.631/K=0.751）；旧 untitled(1).m Bingham 代理以 LEGACY 保留
- 泵注重建为"实际版"（默认，施工记录表/作业史实际）+ "设计版"（设计 7.2）双版本，参照 hu103 双版本做法；
  删"基液 3"步骤（7.2 主流程无基液），修正 notes"替浆总量 48m³"笔误
- CBL：结构化补入数字化结果 cbl_pass_rate=0.787（尾管评价段 5307.54-7514.21m 口径，
  data_type=interpreted）；官方结论为"不合格"，判据是 25m 连续胶结红线而非占比口径
- 目标窗 7400-7618 → 目的井段 7442-7618（区分 CBL 评价窗与地层目标窗）
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
# 原始参考目录保留用于 job_report_path 等原始资料追溯；井径/井斜剖面已切到现场提取包。
DEFAULT_REFERENCE_ROOT = PROJECT_ROOT / "参考文档" / "呼1-003" / "新"
# 现场提取包 caliper_profile.csv：设计一.1.4.2 电测井径 66 点（5570–7618m，30m 间隔）。
# LEGACY(2026-08-16 前): 参考文档/呼1-003/新/呼1-003井身结构.csv（环空体积反算等效井径，非现场提取包）。
DEFAULT_CALIPER_CSV = PROJECT_ROOT / "参考文档" / "现场资料提取" / "ht1_003_呼1-003" / "caliper_profile.csv"
# 现场提取包 inclination_profile.csv：设计一.1.4.1 井斜方位温度 67 点（5560–7618m）。
DEFAULT_INCLINATION_CSV = PROJECT_ROOT / "参考文档" / "现场资料提取" / "ht1_003_呼1-003" / "inclination_profile.csv"
# CBL 数字化产物：尾管评价段 5307.54–7514.21m 一界面中等及以上占比 78.7%（interpreted）。
DEFAULT_CBL_DIGITIZATION_CSV = PROJECT_ROOT / "参考文档" / "现场资料提取" / "ht1_003_呼1-003" / "cbl_digitization.csv"

# 呼1-003井段几何参数（单位：m, mm）—— 与核对报告 §1.1 确认一致
HT1_003_WELL_NAME = "呼1-003井（HT1-003）"
HT1_003_DRILLED_DEPTH_MD_M = 7618.0  # 实际完钻井深/尾管鞋深度 (TD)
HT1_003_HANGER_MD_M = 5307.539  # 尾管悬挂器喇叭口（悬挂器本体 5315.439m，作业史；CBL 头 5307.540-5315.440m）
HT1_003_TOP_MD_M = HT1_003_HANGER_MD_M  # 模型剖面从悬挂器开始（设计 4.4.2 建议模拟段 5307.539-7618m）
HT1_003_CASING_SHOE_MD_M = 5568.0  # 273.1mm套管鞋深度（裸眼自 5568m）
HT1_003_UPPER_SECTION_BOTTOM_MD_M = 7089.576  # 168.3mm 上段尾管底界/变径位置（变径变扣 7091.016m）
HT1_003_LOWER_SECTION_TOP_MD_M = 7096.0  # 241.3mm裸眼→215.9mm裸眼变径位置
HT1_003_BOTTOM_MD_M = HT1_003_DRILLED_DEPTH_MD_M  # 下段 139.7mm 尾管鞋
HT1_003_SHOE_MD_M = HT1_003_DRILLED_DEPTH_MD_M
# 2026-08-29 语义统一：casing_id_mm 按 PACKAGE_REFERENCE 文档语义存"外层套管内径"——
# 273.05mm 技套（名义 273.1）真实 ID=245.37（设计五.5.1 原文+CBL 图头 273.050-2×13.84）；
# OD 公称 273.1 存档于 CASING_OD_MM。字段不被求解器消费（环空计算由 casing_inner_diameter 与
# liner_od_profile/hole_diameter_profile 完成），纯元数据口径。
# LEGACY(2026-08-29 前): 本字段存 OD 273.1（名义口径，命名历史遗留）。
HT1_003_CASING_ID_MM = 245.37  # 技术套管真实内径（语义统一后口径）
HT1_003_CASING_OD_MM = 273.1  # 技术套管外径（OD 公称，存档）
HT1_003_CASING_INNER_DIAMETER_MM = 245.37  # 273.1mm套管内径（设计五.5.1 原文 + CBL 图头 273.050-2×13.84，2026-08-29 校准）；LEGACY(2026-08-29 前): 245.42（计算值）
HT1_003_UPPER_HOLE_NOMINAL_DIAMETER_MM = 241.3  # 上段井眼名义尺寸（钻头程序 241.3*7096）
HT1_003_LOWER_HOLE_DIAMETER_MM = 215.9  # 下段井眼名义尺寸
HT1_003_BIT_DIAMETER_LOWER_MM = 215.9  # 下段钻头尺寸
HT1_003_UPPER_LINER_OD_MM = 168.3  # 上段尾管外径（实际下入 168.28mm）
HT1_003_LOWER_LINER_OD_MM = 139.7  # 下段尾管外径
HT1_003_LOWER_LINER_WALL_THICKNESS_MM = 15.88  # 139.7mm 管壁厚
HT1_003_LOWER_LINER_ID_MM = 107.94  # 139.7 - 2*15.88
HT1_003_UPPER_LINER_WALL_THICKNESS_MM = 14.7  # 168.3mm 尾管壁厚（设计 6.1/作业史）
HT1_003_UPPER_LINER_ID_MM = HT1_003_UPPER_LINER_OD_MM - 2.0 * HT1_003_UPPER_LINER_WALL_THICKNESS_MM
# 扶正器数量（现场，设计 6.6.2/作业史）：168.3mm 段 75 只 + 139.7mm 段 22 只 = 97 只。
# LEGACY(2026-08-16 前): UPPER=21、LOWER=78、总 99（数量与分段均与现场不符）。
HT1_003_LOWER_CENTRALIZER_COUNT = 22  # 139.7mm 下段整体式扶正器（间距2根/只，22.8m）
HT1_003_UPPER_CENTRALIZER_COUNT = 75  # 168.3mm 上段整体式扶正器（间距2根/只，22.8m）
HT1_003_CENTRALIZER_COUNT = HT1_003_LOWER_CENTRALIZER_COUNT + HT1_003_UPPER_CENTRALIZER_COUNT


def _read_caliper_rows(caliper_csv_path: Path) -> tuple[tuple[float, float], ...]:
    """读取现场提取包井径剖面 CSV（md_m / caliper_mm），按 md_m 升序返回。

    列名与 hu102/hu103 现场提取包 caliper_profile.csv 一致（裸眼实测井径）。
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


def _build_hole_profile(caliper_rows: tuple[tuple[float, float], ...]) -> tuple[tuple[float, float], ...]:
    """从现场提取包井径行构建井径剖面（裸眼实测 caliper_mm）。

    悬挂器-技套鞋段（5307.539-5568m）为 273.1mm 套管内重叠段，无裸眼井径测点，
    取技套内径 245.42mm 作为等效井径（与 168.3mm 尾管环空间隙一致，由 liner_od_profile 表达）。
    """
    points: list[tuple[float, float]] = [(md, cal) for md, cal in caliper_rows if md >= HT1_003_HANGER_MD_M]
    if not points or points[0][0] > HT1_003_HANGER_MD_M:
        points.insert(0, (HT1_003_HANGER_MD_M, HT1_003_CASING_INNER_DIAMETER_MM))
    if points[-1][0] < HT1_003_BOTTOM_MD_M:
        points.append((HT1_003_BOTTOM_MD_M, points[-1][1]))
    return tuple(points)


def _build_inclination_profile(incl_rows: tuple[tuple[float, float], ...]) -> tuple[tuple[float, float], ...]:
    """从现场提取包井斜行构建井斜剖面（悬挂器以上取 0° 外推）。"""
    points: list[tuple[float, float]] = [(md, inc) for md, inc in incl_rows if md >= HT1_003_HANGER_MD_M]
    if not points or points[0][0] > HT1_003_HANGER_MD_M:
        points.insert(0, (HT1_003_HANGER_MD_M, 0.0))
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


# 管柱内径剖面常量（与 _build_pipe_id_profile / 鞋口滞后体积共用）
_DP1_ID_MM = 149.2 - 2 * 9.65   # 129.9mm（149.2mm 钻杆，壁厚 9.65mm）
_DP2_ID_MM = 127.0 - 2 * 9.65   # 107.7mm（127mm 钻杆，壁厚 9.65mm）
_DP1_BOTTOM_MD_M = 3321.682


def _build_pipe_id_profile() -> tuple[tuple[float, float], ...]:
    """构建管柱内径剖面（深度, ID mm）。

    从地面到TD的完整管柱内径，用于1D前沿追踪计算管内容积。
    段1: 0~3321.682m      → 129.9mm (149.2mm钻杆, 壁厚9.65mm)
    段2: 3321.682~5307.539m → 107.7mm (127mm钻杆, 壁厚9.65mm)
    段3: 5307.539~7089.576m → 138.9mm (168.3mm尾管, 壁厚14.7mm)
    段4: 7089.576~7618m    → 107.94mm (139.7mm尾管, 壁厚15.88mm)
    """
    return (
        (0.001, _DP1_ID_MM),
        (_DP1_BOTTOM_MD_M, _DP1_ID_MM),
        (_DP1_BOTTOM_MD_M + 0.001, _DP2_ID_MM),
        (HT1_003_HANGER_MD_M, _DP2_ID_MM),
        (HT1_003_HANGER_MD_M + 0.001, HT1_003_UPPER_LINER_ID_MM),
        (HT1_003_UPPER_SECTION_BOTTOM_MD_M, HT1_003_UPPER_LINER_ID_MM),
        (HT1_003_UPPER_SECTION_BOTTOM_MD_M + 0.001, HT1_003_LOWER_LINER_ID_MM),
        (HT1_003_BOTTOM_MD_M, HT1_003_LOWER_LINER_ID_MM),
    )


def _equivalent_hole_diameter_mm(actual_hole_mm: float, actual_od_mm: float, reference_od_mm: float) -> float:
    """在固定 reference_od 条件下，构造与原始环空面积等价的井眼直径。"""
    area_term = actual_hole_mm**2 - actual_od_mm**2 + reference_od_mm**2
    return math.sqrt(max(area_term, reference_od_mm**2))


# LEGACY(2026-08-16 前): 等效井径常量（环空体积反算口径，剖面已切到现场提取包裸眼实测井径，
# 不再用于 _build_hole_profile；保留常量供旧口径追溯）。
HT1_003_UPPER_HOLE_DIAMETER_MM = _equivalent_hole_diameter_mm(
    actual_hole_mm=HT1_003_UPPER_HOLE_NOMINAL_DIAMETER_MM,
    actual_od_mm=HT1_003_UPPER_LINER_OD_MM,
    reference_od_mm=HT1_003_LOWER_LINER_OD_MM,
)
# LEGACY(2026-08-16 前): 套管段等效井径：将168.3mm尾管在273.1mm套管内的环空面积，
# 等效为139.7mm尾管的虚拟井径（旧 _build_hole_diameter_profile 用）。
HT1_003_CASING_SECTION_EQUIVALENT_HOLE_MM = _equivalent_hole_diameter_mm(
    actual_hole_mm=HT1_003_CASING_INNER_DIAMETER_MM,
    actual_od_mm=HT1_003_UPPER_LINER_OD_MM,
    reference_od_mm=HT1_003_LOWER_LINER_OD_MM,
)


def _pipe_volume_m3(length_m: float, inner_diameter_mm: float) -> float:
    """按内径和长度估算管内容积，用于鞋口滞后体积计算。"""
    radius_m = inner_diameter_mm / 1000.0 / 2.0
    return math.pi * radius_m**2 * length_m


# 鞋口滞后体积（2026-08-29 校准）：设计七.1.4 顶替量计算表 42.9+18.1+27+3.9=91.9m³
# （胶塞试通过 61m³ 实测验证；149.2 段用 12.91 L/m 折算口径；数采计量 91.009m³ 与理论差 1%）。
# LEGACY(2026-08-29 前): 4 段管柱 ID 直接累加 ≈93.95m³（149.2 段用 12.91 L/m 折算与水力 ID 129.9 双口径差异）；
# 更早 LEGACY(2026-08-16 前): 7868 常数 bug（复制 hu101 TD 的笔误，鞋口滞后 66.92m³ 偏低约 25m³）已修复。
_HT1_003_SHOE_LAG_PIPE_SUM_LEGACY_M3 = (
    _pipe_volume_m3(_DP1_BOTTOM_MD_M, _DP1_ID_MM)
    + _pipe_volume_m3(HT1_003_HANGER_MD_M - _DP1_BOTTOM_MD_M, _DP2_ID_MM)
    + _pipe_volume_m3(
        HT1_003_UPPER_SECTION_BOTTOM_MD_M - HT1_003_HANGER_MD_M,
        HT1_003_UPPER_LINER_ID_MM,
    )
    + _pipe_volume_m3(
        HT1_003_BOTTOM_MD_M - HT1_003_UPPER_SECTION_BOTTOM_MD_M,
        HT1_003_LOWER_LINER_ID_MM,
    )
)  # ≈93.95，LEGACY(2026-08-29 前) 保留供口径追溯
HT1_003_SURFACE_TO_HANGER_EFFECTIVE_ID_MM = math.sqrt(4.0 * 52.0 / (math.pi * 7868.0)) * 1000.0  # LEGACY
HT1_003_SHOE_LAG_VOLUME_M3 = 91.9  # 设计七.1.4 顶替量计算表（42.9+18.1+27+3.9），见上注释与 notes。
HT1_003_LINER_ID_MM = HT1_003_LOWER_LINER_ID_MM

# ===========================================================================
# 呼1-003流体参数（2026-08-16 按现场提取包/化验报告更新；旧 untitled(1).m 代理值以 LEGACY 保留）
# 密度：现场值均 field_measured；流变：钻井液/领浆/尾浆/隔离液为化验报告幂律实测（field），
# 先导浆/压塞液/保护液/基液为设计值（design_value），替浆钻井液为施工记录表实测塑粘（field）。
# ===========================================================================
HT1_003_MUD_DENSITY_KG_M3 = 1950.0  # 钻井液(井浆) 1.95（现场）
HT1_003_BALANCE_DENSITY_KG_M3 = 1750.0  # 先导浆 1.75（现场）；LEGACY(2026-08-16 前): 1850.0
HT1_003_SPACER1_DENSITY_KG_M3 = 2050.0  # 隔离液1 2.05（现场化验）
HT1_003_SPACER2_DENSITY_KG_M3 = 1950.0  # 隔离液2 1.95（现场化验）
HT1_003_SPACER_DENSITY_KG_M3 = 2000.0  # LEGACY(2026-08-16 前): 单一隔离液 2.00，已拆分为 SPACER1/SPACER2
HT1_003_LEAD_DENSITY_KG_M3 = 2050.0  # 领浆 2.05（现场）
HT1_003_TAIL_DENSITY_KG_M3 = 1950.0  # 尾浆 1.95（现场化验实测）；LEGACY(2026-08-16 前): 2050.0（误用 2.05）
HT1_003_PLUG_DENSITY_KG_M3 = 1950.0  # 压塞液 1.95（设计 7.2）
HT1_003_DISPLACEMENT_DENSITY_KG_M3 = 1950.0  # 替浆钻井液 1.95（现场）
HT1_003_BUFFER_DENSITY_KG_M3 = 1950.0  # 保护液 1.95（设计 7.2）
HT1_003_BASE_FLUID_DENSITY_KG_M3 = 1020.0  # 基液 1.02（设计 7.2）
HT1_003_WELL_MUD_DENSITY_KG_M3 = 1950.0  # 井浆(替浆段) 1.95（现场）；LEGACY(2026-08-16 前): 1850.0

# 钻井液/水泥浆/隔离液流变（2026-08-16 改用化验报告幂律实测，field_measured）
# 化验报告（129C↓93C）：钻井液六速 122/71/52/21/7/6 → n=0.631/K=0.751；
# 领浆 207/121/93/60/15/9 → n=0.597/K=1.622；尾浆 196/118/89/56/14/10 → n=0.585/K=1.673；
# 隔离液 239/177/123/89/10/6 → n=0.668/K=1.245（隔离液1/2 同流变）。
HT1_003_MUD_POWER_LAW_N = 0.631  # 钻井液 流性指数 n（化验报告）；LEGACY(2026-08-16 前): 0.82
HT1_003_MUD_CONSISTENCY_K = 0.751  # 钻井液 稠度系数 K Pa·s^n（化验报告）；LEGACY: 0.21
HT1_003_LEAD_POWER_LAW_N = 0.597  # 领浆 流性指数 n（化验报告）
HT1_003_LEAD_CONSISTENCY_K = 1.622  # 领浆 稠度系数 K（化验报告）
HT1_003_TAIL_POWER_LAW_N = 0.585  # 尾浆 流性指数 n（化验报告）
HT1_003_TAIL_CONSISTENCY_K = 1.673  # 尾浆 稠度系数 K（化验报告）
HT1_003_SPACER_POWER_LAW_N = 0.668  # 隔离液 流性指数 n（化验报告）
HT1_003_SPACER_CONSISTENCY_K = 1.245  # 隔离液 稠度系数 K（化验报告）
# LEGACY(2026-08-16 前): 钻井液/领浆/尾浆/隔离液 Bingham 代理（untitled(1).m miu/tau 参数）
HT1_003_MUD_PV_PA_S = 0.051  # LEGACY 钻井液 PV=51 mPa·s（设计 Bingham）
HT1_003_MUD_YP_PA = 10.0  # LEGACY 钻井液 YP=10 Pa
HT1_003_LEAD_PV_PA_S = 0.160  # LEGACY 领浆 PV=160 mPa·s
HT1_003_LEAD_YP_PA = 13.0  # LEGACY 领浆 YP=13 Pa
HT1_003_TAIL_PV_PA_S = 0.180  # LEGACY 尾浆 PV=180 mPa·s
HT1_003_TAIL_YP_PA = 14.0  # LEGACY 尾浆 YP=14 Pa
HT1_003_SPACER_PV_PA_S = 0.060  # LEGACY 隔离液 PV=60 mPa·s
HT1_003_SPACER_YP_PA = 11.0  # LEGACY 隔离液 YP=11 Pa

# 先导浆流变：设计初稿 Bingham PV55/YP9.2（design_value，现场无流变实测）
HT1_003_BALANCE_PV_PA_S = 0.055  # 先导浆 PV=55 mPa·s（设计初稿）；LEGACY: 0.051
HT1_003_BALANCE_YP_PA = 9.2  # 先导浆 YP=9.2 Pa（设计初稿）；LEGACY: 10.0

# 压塞液/保护液/基液流变：设计值（design_value，现场无流变实测）
HT1_003_PLUG_PV_PA_S = 0.030  # 压塞液 PV=30 mPa·s（设计）；LEGACY: 0.040
HT1_003_PLUG_YP_PA = 8.0  # 压塞液 YP=8 Pa（设计）；LEGACY: 9.0
HT1_003_BUFFER_PV_PA_S = 0.030  # 保护液 PV=30 mPa·s（设计）；LEGACY: 0.040
HT1_003_BUFFER_YP_PA = 8.2  # 保护液 YP=8.2 Pa（设计）；LEGACY: 9.2
HT1_003_BASE_FLUID_PV_PA_S = 0.030  # 基液 PV=30 mPa·s（设计）；LEGACY: 0.030
HT1_003_BASE_FLUID_YP_PA = 8.0  # 基液 YP=8 Pa（设计）；LEGACY: 9.0

# 替浆钻井液：塑粘 62 mPa·s（施工记录表 替浆 field_measured），YP 取洗井屈服值 10 Pa（现场）
HT1_003_DISPLACEMENT_PV_PA_S = 0.062  # 替浆钻井液 PV=62 mPa·s（记录表）；LEGACY: 0.040
HT1_003_DISPLACEMENT_YP_PA = 10.0  # 替浆钻井液 YP=10 Pa（洗井屈服值）；LEGACY: 9.5
# 井浆(替浆段)流变：沿用旧代理 PV30/YP9.3（LEGACY，无实测）
HT1_003_WELL_MUD_PV_PA_S = 0.030
HT1_003_WELL_MUD_YP_PA = 9.3

# ===========================================================================
# 呼1-003施工程序参数（2026-08-16 按现场提取包 pumping_schedule.csv 重建）
# 实际版（默认，施工记录表/作业史）：先导浆28/隔离液16+10两段/领浆39/尾浆28/压塞2/替浆91.9。
# 设计版（load_ht1_003_tailpipe_design，设计 7.2）：替浆分 7 步（压塞2/钻井液25/保护14/钻井液8/14/16/12.9）。
# LEGACY(2026-08-16 前): untitled(1).m 代理体积（35/35/38/28.5/…）与"基液 3"步骤已废弃。
# ===========================================================================
HT1_003_BALANCE_VOLUME_M3 = 28.0  # 先导浆 28m³（现场/设计 7.2）；LEGACY: 35.0
HT1_003_SPACER1_VOLUME_M3 = 16.0  # 隔离液1 16m³（现场）
HT1_003_SPACER2_VOLUME_M3 = 10.0  # 隔离液2 10m³（现场）
HT1_003_SPACER_VOLUME_M3 = 35.0  # LEGACY(2026-08-16 前): 单一隔离液 35m³，已拆分为 SPACER1/SPACER2
HT1_003_LEAD_VOLUME_M3 = 39.0  # 领浆 39m³（现场/设计 7.2）；LEGACY: 38.0
HT1_003_TAIL_VOLUME_M3 = 28.0  # 尾浆 28m³（现场/设计 7.2）；LEGACY: 28.5
HT1_003_PLUG_VOLUME_M3 = 2.0  # 压塞液 2m³（现场/设计 7.2）
HT1_003_DISPLACEMENT_VOLUME_M3 = 91.9  # 替浆总 91.9m³（设计 7.2 理论顶替量含压塞2/中置14；数采 91.01）
# 设计 7.2 替浆分步体积（design_value，见 pumping_schedule.csv step 106-111）
HT1_003_DESIGN_DISP_PLUG_VOLUME_M3 = 2.0  # 压塞液 2m³
HT1_003_DESIGN_DISP_MUD1_VOLUME_M3 = 25.0  # 钻井液 25m³（固井计量罐）
HT1_003_DESIGN_DISP_BUFFER_VOLUME_M3 = 14.0  # 保护液 14m³（批混车混配）
HT1_003_DESIGN_DISP_MUD2A_VOLUME_M3 = 8.0  # 钻井液 8m³
HT1_003_DESIGN_DISP_MUD2B_VOLUME_M3 = 14.0  # 钻井液 14m³
HT1_003_DESIGN_DISP_MUD2C_VOLUME_M3 = 16.0  # 钻井液 16m³
HT1_003_DESIGN_DISP_MUD2D_VOLUME_M3 = 12.9  # 钻井液 12.9m³（末段；理论顶替总量 91.9m³）

# 排量（2026-08-16 重建；实际版按施工记录表/作业史实际排量，设计版按设计 7.2）
HT1_003_BALANCE_RATE_M3_MIN = 1.2  # 先导浆 实际 1.2（设计亦 1.2）
HT1_003_SPACER_RATE_M3_MIN = 1.0  # 隔离液 实际 1.0（设计 1.2）；LEGACY: 1.2
HT1_003_CEMENT_RATE_M3_MIN = 1.2  # 领浆 1.2（实际/设计）
HT1_003_TAIL_RATE_M3_MIN = 1.2  # 尾浆 实际 1.2（设计 1.4）
HT1_003_PLUG_RATE_M3_MIN = 1.3  # 压塞液（管内占位）
HT1_003_DISPLACEMENT_RATE_M3_MIN = 1.4  # 替浆单段 实际 1.4（实际 1.4-0.8 递减）
# 设计 7.2 排量（design_value）
HT1_003_DESIGN_SPACER_RATE_M3_MIN = 1.2  # 隔离液(设计) 1.2
HT1_003_DESIGN_TAIL_RATE_M3_MIN = 1.4  # 尾浆(设计) 1.4
HT1_003_DESIGN_DISP_PLUG_RATE_M3_MIN = 1.0  # 压塞(设计) 1.0-1.4
HT1_003_DESIGN_DISP_MUD1_RATE_M3_MIN = 1.6  # 钻井液(设计) 1.6
HT1_003_DESIGN_DISP_BUFFER_RATE_M3_MIN = 1.4  # 保护液(设计) 1.4
HT1_003_DESIGN_DISP_MUD2A_RATE_M3_MIN = 1.2  # 钻井液(设计) 1.2
HT1_003_DESIGN_DISP_MUD2B_RATE_M3_MIN = 1.0  # 钻井液(设计) 1.0
HT1_003_DESIGN_DISP_MUD2C_RATE_M3_MIN = 0.8  # 钻井液(设计) 0.8
HT1_003_DESIGN_DISP_MUD2D_RATE_M3_MIN = 0.7  # 钻井液(设计) 0.7（末段）
# LEGACY(2026-08-16 前): untitled(1).m 替浆分段（替钻井液28/保护12/基液3/井浆8+14+14+12=91，基液步骤已删）
HT1_003_FAST_MUD_VOLUME_M3 = 28.0
HT1_003_BUFFER_VOLUME_M3 = 12.0
HT1_003_BASE_FLUID_VOLUME_M3 = 3.0
HT1_003_WELL_MUD_FAST_VOLUME_M3 = 8.0
HT1_003_WELL_MUD_MID1_VOLUME_M3 = 14.0
HT1_003_WELL_MUD_MID2_VOLUME_M3 = 14.0
HT1_003_WELL_MUD_MID3_VOLUME_M3 = 12.0
HT1_003_WELL_MUD_SLOW_VOLUME_M3 = 0.0
HT1_003_FAST_MUD_RATE_M3_MIN = 1.3  # LEGACY 替钻井液 1.3
HT1_003_BUFFER_RATE_M3_MIN = 1.2  # LEGACY 保护液 1.2
HT1_003_BASE_FLUID_RATE_M3_MIN = 1.2  # LEGACY 基液 1.2
HT1_003_WELL_MUD_FAST_RATE_M3_MIN = 1.3  # LEGACY 井浆快替 1.3
HT1_003_WELL_MUD_MID1_RATE_M3_MIN = 1.1  # LEGACY 井浆中替1 1.1
HT1_003_WELL_MUD_MID2_RATE_M3_MIN = 0.9  # LEGACY 井浆中替2 0.9
HT1_003_WELL_MUD_MID3_RATE_M3_MIN = 0.7  # LEGACY 井浆中替3 0.7
HT1_003_WELL_MUD_SLOW_RATE_M3_MIN = 0.7  # LEGACY 保留段 0.7


def _depth_points(values: tuple[tuple[float, float], ...]) -> tuple[DepthValuePoint, ...]:
    """把测深-数值元组转换为 WellSpec 使用的剖面点。"""
    return tuple(DepthValuePoint(depth_md_m=depth, value=value) for depth, value in values)


def _build_well_spec(
    resolved_reference_root: Path,
    caliper_csv_path: Path,
    inclination_csv_path: Path,
) -> WellSpec:
    """构建呼1-003 井筒规格（几何/剖面/评价窗，实际版与设计版共用）。"""
    caliper_rows = _read_caliper_rows(caliper_csv_path)
    incl_rows = _read_inclination_rows(inclination_csv_path)
    hole_profile = _build_hole_profile(caliper_rows)
    inc_profile = _build_inclination_profile(incl_rows)
    liner_od_profile = _build_liner_od_profile()
    pipe_id_profile = _build_pipe_id_profile()

    return WellSpec(
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
            EvaluationWindow(name="CBL评价井段(尾管段)", top_md_m=5568.0, bottom_md_m=7618.0, window_type="cbl"),
            # 数字化窗（cbl_pass_rate=0.787 对应段：水泥面-人工井底，interpreted，cbl_digitization.csv）；
            # 与官方"25m 连续判据不合格"口径不同，论文引用前须裁定（见 ValidationData notes）。
            EvaluationWindow(name="CBL数字化窗(78.7%对应)", top_md_m=5307.54, bottom_md_m=7514.21, window_type="cbl_digitization"),
            # 地层目标（设计 二.2.3 油气显示，field_measured）。
            EvaluationWindow(name="气层-K1q主力油层段(地层目标)", top_md_m=7418.0, bottom_md_m=7440.0, window_type="formation_target"),
            EvaluationWindow(name="目的井段(地层目标)", top_md_m=7442.0, bottom_md_m=7618.0, window_type="formation_target"),
        ),
        reference_root=resolved_reference_root,
        notes=(
            "呼1-003井（HT1-003）为168.3mm+139.7mm双径复合尾管控压固井，井深7618m。",
            "井径/井斜剖面使用现场提取包 caliper_profile.csv（设计一.1.4.2 电测井径 66 点）与 inclination_profile.csv（设计一.1.4.1 井斜 67 点）。",
            "liner_od_profile 按深度分段：5307-7089m=168.3mm，7089-7618m=139.7mm；pipe_id_profile 按4段管柱内径：129.9/107.7/138.9/107.94mm。",
            "悬挂器取喇叭口 5307.539m（本体 5315.439m、CBL 头 5307.540-5315.440m）；变径 7089.576m（变径变扣 7091.016m）。",
            "鞋口滞后体积 91.9m³（2026-08-29 校准：设计七.1.4 顶替量计算表 42.9+18.1+27+3.9；胶塞试通过 61m³ 实测验证；"
            "数采计量 91.009m³ 与理论差 1%）；LEGACY(2026-08-29 前): 4 段管柱 ID 累加 ≈93.95m³"
            "（149.2 段体积口径 12.91 L/m 折算 vs 水力 ID 129.9 双口径并存）；更早 LEGACY 7868 常数 bug 已于 2026-08-16 修复。",
            "压塞液密度保留 1.95（设计 7.2 口径）；记录表'后置液 密度1.02 基液'为另一口径"
            "（数采密度计 0.998-1.000 存疑），需甲方确认（2026-08-29 注记）。",
            "2026-08-29 校准补记：井径均值 CSV 重算 247.00/220.27 vs 设计自称 247.3/218.9 为统计口径差；"
            "大肚子 5573-5577 平均 370.84/最大 387.43（0.5m 曲线极值，30m 表未含）；"
            "悬挂器本体最小内径 138.22（5307.539-5316.036 段容差 0.05m³ 可忽略）。",
            "流体流变改用化验报告幂律实测（领/尾浆/隔离液/钻井液），先导浆/压塞/保护/基液为设计值。",
            "居中度取设计模拟值 0.83（设计 6.3，裸眼段无连续实测；model_assumption）。",
        ),
    )


def _build_fluids() -> tuple[FluidSpec, ...]:
    """呼1-003 流体清单（实际版与设计版共用；2026-08-16 化验幂律重建）。"""
    return (
        FluidSpec("钻井液", FluidRole.MUD, HT1_003_MUD_DENSITY_KG_M3, RheologyModel.POWER_LAW,
                  power_law_n=HT1_003_MUD_POWER_LAW_N, consistency_k=HT1_003_MUD_CONSISTENCY_K),
        FluidSpec("平衡液", FluidRole.WASH, HT1_003_BALANCE_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_BALANCE_PV_PA_S, HT1_003_BALANCE_YP_PA),
        FluidSpec("隔离液1", FluidRole.SPACER, HT1_003_SPACER1_DENSITY_KG_M3, RheologyModel.POWER_LAW,
                  power_law_n=HT1_003_SPACER_POWER_LAW_N, consistency_k=HT1_003_SPACER_CONSISTENCY_K),
        FluidSpec("隔离液2", FluidRole.SPACER, HT1_003_SPACER2_DENSITY_KG_M3, RheologyModel.POWER_LAW,
                  power_law_n=HT1_003_SPACER_POWER_LAW_N, consistency_k=HT1_003_SPACER_CONSISTENCY_K),
        FluidSpec("领浆", FluidRole.LEAD, HT1_003_LEAD_DENSITY_KG_M3, RheologyModel.POWER_LAW,
                  power_law_n=HT1_003_LEAD_POWER_LAW_N, consistency_k=HT1_003_LEAD_CONSISTENCY_K),
        FluidSpec("尾浆", FluidRole.TAIL, HT1_003_TAIL_DENSITY_KG_M3, RheologyModel.POWER_LAW,
                  power_law_n=HT1_003_TAIL_POWER_LAW_N, consistency_k=HT1_003_TAIL_CONSISTENCY_K),
        FluidSpec("压塞液", FluidRole.OTHER, HT1_003_PLUG_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_PLUG_PV_PA_S, HT1_003_PLUG_YP_PA),
        FluidSpec("保护液", FluidRole.DISPLACEMENT, HT1_003_BUFFER_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_BUFFER_PV_PA_S, HT1_003_BUFFER_YP_PA),
        FluidSpec("基液", FluidRole.DISPLACEMENT, HT1_003_BASE_FLUID_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_BASE_FLUID_PV_PA_S, HT1_003_BASE_FLUID_YP_PA),
        FluidSpec("替钻井液", FluidRole.DISPLACEMENT, HT1_003_DISPLACEMENT_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_DISPLACEMENT_PV_PA_S, HT1_003_DISPLACEMENT_YP_PA),
        FluidSpec("井浆", FluidRole.DISPLACEMENT, HT1_003_WELL_MUD_DENSITY_KG_M3, RheologyModel.BINGHAM,
                  HT1_003_WELL_MUD_PV_PA_S, HT1_003_WELL_MUD_YP_PA),
    )


def load_ht1_003_tailpipe(
    *,
    reference_root: Path | None = None,
    caliper_csv_path: Path | None = None,
    inclination_csv_path: Path | None = None,
) -> tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载呼1-003井（HT1-003）168.3+139.7mm双径尾管段"实际施工版"模型输入（默认）。

    实际版泵注程序（施工记录表/作业史）：先导浆28@1.2 → 隔离液1 16@1.0 → 隔离液2 10@1.0 →
    领浆39@1.2 → 尾浆28@1.2 → 压塞液2（管内占位）→ 替浆91.9@1.4（实际单段 1.4-0.8 递减）。

    设计版（设计 7.2）见 load_ht1_003_tailpipe_design。

    Args:
        reference_root: 可选参考资料根目录（默认参考文档/呼1-003/新）。
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

    # 呼1-003 实际施工泵注程序（2026-08-16 重建，来源：施工记录表/作业史/固井总结）
    schedule = PumpingSchedule(
        steps=(
            PumpingScheduleStep("注入先导浆（平衡液）", "平衡液", HT1_003_BALANCE_VOLUME_M3, HT1_003_BALANCE_RATE_M3_MIN,
                                remarks="实际: 先导浆 28m³@1.2m³/min，密度1.75g/cm³，泵压14-17.8MPa。"),
            PumpingScheduleStep("注入隔离液1", "隔离液1", HT1_003_SPACER1_VOLUME_M3, HT1_003_SPACER_RATE_M3_MIN,
                                remarks="实际: 隔离液1 16m³@1.0m³/min，密度2.05g/cm³，泵压19MPa（设计排量1.2）。"),
            PumpingScheduleStep("注入隔离液2", "隔离液2", HT1_003_SPACER2_VOLUME_M3, HT1_003_SPACER_RATE_M3_MIN,
                                remarks="实际: 隔离液2 10m³@1.0m³/min，密度1.95g/cm³，泵压18-19MPa。"),
            PumpingScheduleStep("注入领浆", "领浆", HT1_003_LEAD_VOLUME_M3, HT1_003_CEMENT_RATE_M3_MIN,
                                remarks="实际: 领浆 39m³@1.2m³/min，密度2.05g/cm³，泵压22MPa；灰量157t。"),
            PumpingScheduleStep("注入尾浆", "尾浆", HT1_003_TAIL_VOLUME_M3, HT1_003_TAIL_RATE_M3_MIN,
                                remarks="实际: 尾浆 28m³@1.2m³/min，密度1.95g/cm³，泵压22MPa。"),
            PumpingScheduleStep("注入压塞液（管内）", "压塞液", HT1_003_PLUG_VOLUME_M3, HT1_003_PLUG_RATE_M3_MIN,
                                remarks="压塞液 2m³@1.3m³/min，仅作为管内占位（碰压19-25MPa）。"),
            PumpingScheduleStep("替浆（单段，实际）", "替钻井液", HT1_003_DISPLACEMENT_VOLUME_M3,
                                HT1_003_DISPLACEMENT_RATE_M3_MIN,
                                remarks="实际: 替浆 91.9m³@1.4m³/min（单段，实际1.4-0.8递减，约80min），泵压15-19-26-19MPa；"
                                        "现场理论顶替量(含压塞2+中置14+钻井液75.9)；数采总泵量91.01m³。"),
        ),
        notes=(
            "实际施工顺序：先导浆→隔离液1→隔离液2→领浆→尾浆→压塞液（管内占位）→替浆（单段91.9m³）。",
            "替浆 91.9m³ 为现场理论顶替量（设计7.2含压塞2+中置/保护14+钻井液75.9）；数采佐证 91.01m³。",
            "LEGACY(2026-08-16 前): untitled(1).m 代理替浆（28+12+3+8+14+14+12=91，含'基液3'步骤）已重建；"
            "旧 notes'替浆总量 48m³'为笔误，已修正为 91.9m³。",
        ),
    )

    validation_data = ValidationData(
        cbl_summary_path=DEFAULT_CBL_DIGITIZATION_CSV,
        cbl_pass_rate=0.787,
        job_report_path=resolved_reference_root / "呼1-003井身结构.csv",
        notes=(
            "cbl_pass_rate=0.787 为尾管评价段 5307.54-7514.21m（水泥面-人工井底）一界面中等及以上占比"
            "（数字化口径 78.7%，interpreted，来源 cbl_digitization.csv）。",
            "官方 CBL 结论为'不合格'，判据是七条红线（上层套管鞋/尾管重合段及以上 25m 连续胶结中等及以上、"
            "油气水层段一、二界面连续胶结），与总占比口径不同；官方未给数值合格率，占比较值仅作参考。",
            "现场提取包井径/井斜/泵注/流体按 2026-08-16 核对报告落地；几何与流变以 field_measured 为准。",
        ),
    )
    return well_spec, fluids, schedule, validation_data


def load_ht1_003_tailpipe_design(
    *,
    reference_root: Path | None = None,
    caliper_csv_path: Path | None = None,
    inclination_csv_path: Path | None = None,
) -> tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载呼1-003井（HT1-003）168.3+139.7mm双径尾管段"设计版"模型输入（设计 7.2）。

    设计版泵注程序（设计 7.2）：先导浆28@1.2 → 隔离液1 16@1.2 → 隔离液2 10@1.2 → 领浆39@1.2 →
    尾浆28@1.4 → 压塞2@1.0 → 替钻井液25@1.6 → 保护14@1.4 → 替钻井液8@1.2/14@1.0/16@0.8/12.9@0.7
    （理论顶替总量 91.9m³）。

    与默认实际版（load_ht1_003_tailpipe）的差异仅在 PumpingSchedule 排量与分段；几何/流体/校验资料共用。

    Returns:
        (well_spec, fluids, design_schedule, validation_data)
    """
    resolved_reference_root = reference_root or DEFAULT_REFERENCE_ROOT
    resolved_caliper_csv_path = caliper_csv_path or DEFAULT_CALIPER_CSV
    resolved_inclination_csv_path = inclination_csv_path or DEFAULT_INCLINATION_CSV

    well_spec = _build_well_spec(resolved_reference_root, resolved_caliper_csv_path, resolved_inclination_csv_path)
    fluids = _build_fluids()

    # 呼1-003 设计版泵注程序（设计 7.2，design_value）
    schedule = PumpingSchedule(
        steps=(
            PumpingScheduleStep("注入先导浆（设计）", "平衡液", HT1_003_BALANCE_VOLUME_M3, HT1_003_BALANCE_RATE_M3_MIN,
                                remarks="设计: 先导浆 28m³@1.2m³/min，密度1.75g/cm³，泵压12.7-16.1MPa。"),
            PumpingScheduleStep("注入隔离液1（设计）", "隔离液1", HT1_003_SPACER1_VOLUME_M3, HT1_003_DESIGN_SPACER_RATE_M3_MIN,
                                remarks="设计: 隔离液1 16m³@1.2m³/min，密度2.05g/cm³（实际排量1.0）。"),
            PumpingScheduleStep("注入隔离液2（设计）", "隔离液2", HT1_003_SPACER2_VOLUME_M3, HT1_003_DESIGN_SPACER_RATE_M3_MIN,
                                remarks="设计: 隔离液2 10m³@1.2m³/min，密度1.95g/cm³（实际排量1.0）。"),
            PumpingScheduleStep("注入领浆（设计）", "领浆", HT1_003_LEAD_VOLUME_M3, HT1_003_CEMENT_RATE_M3_MIN,
                                remarks="设计: 领浆 39m³@1.2m³/min，密度2.05g/cm³，40min，泵压17-15MPa。"),
            PumpingScheduleStep("注入尾浆（设计）", "尾浆", HT1_003_TAIL_VOLUME_M3, HT1_003_DESIGN_TAIL_RATE_M3_MIN,
                                remarks="设计: 尾浆 28m³@1.4m³/min（1.4-1.2），密度1.95g/cm³，30min。"),
            PumpingScheduleStep("注入压塞液（设计）", "压塞液", HT1_003_DESIGN_DISP_PLUG_VOLUME_M3, HT1_003_DESIGN_DISP_PLUG_RATE_M3_MIN,
                                remarks="设计: 压塞液 2m³@1.0m³/min（1.0-1.4），密度1.95g/cm³。"),
            PumpingScheduleStep("替钻井液(设计一)", "替钻井液", HT1_003_DESIGN_DISP_MUD1_VOLUME_M3, HT1_003_DESIGN_DISP_MUD1_RATE_M3_MIN,
                                remarks="设计: 钻井液 25m³@1.6m³/min（固井计量罐）。"),
            PumpingScheduleStep("替保护液(设计)", "保护液", HT1_003_DESIGN_DISP_BUFFER_VOLUME_M3, HT1_003_DESIGN_DISP_BUFFER_RATE_M3_MIN,
                                remarks="设计: 保护液 14m³@1.4m³/min（批混车混配）。"),
            PumpingScheduleStep("替钻井液(设计二)", "替钻井液", HT1_003_DESIGN_DISP_MUD2A_VOLUME_M3, HT1_003_DESIGN_DISP_MUD2A_RATE_M3_MIN,
                                remarks="设计: 钻井液 8m³@1.2m³/min。"),
            PumpingScheduleStep("替钻井液(设计三)", "替钻井液", HT1_003_DESIGN_DISP_MUD2B_VOLUME_M3, HT1_003_DESIGN_DISP_MUD2B_RATE_M3_MIN,
                                remarks="设计: 钻井液 14m³@1.0m³/min。"),
            PumpingScheduleStep("替钻井液(设计四)", "替钻井液", HT1_003_DESIGN_DISP_MUD2C_VOLUME_M3, HT1_003_DESIGN_DISP_MUD2C_RATE_M3_MIN,
                                remarks="设计: 钻井液 16m³@0.8m³/min。"),
            PumpingScheduleStep("替钻井液(设计末段)", "替钻井液", HT1_003_DESIGN_DISP_MUD2D_VOLUME_M3, HT1_003_DESIGN_DISP_MUD2D_RATE_M3_MIN,
                                remarks="设计: 钻井液 12.9m³@0.7m³/min（末段；理论顶替总量含后置液 91.9m³）。"),
        ),
        notes=(
            "设计版泵注程序（设计 7.2）：先导浆→隔离液1→隔离液2→领浆→尾浆→压塞液→替钻井液→保护液→替钻井液分段递减。",
            "设计排量：前置 1.2（隔离液）/ 尾浆 1.4 / 替浆 1.6→0.7 递减；理论顶替总量 91.9m³。",
            "设计 7.2 主流程无基液步骤（LEGACY untitled(1).m 代理含'基液3m³'步骤，已删）。",
        ),
    )

    validation_data = ValidationData(
        cbl_summary_path=DEFAULT_CBL_DIGITIZATION_CSV,
        cbl_pass_rate=0.787,
        job_report_path=resolved_reference_root / "呼1-003井身结构.csv",
        notes=(
            "cbl_pass_rate=0.787 为尾管评价段 5307.54-7514.21m 一界面中等及以上占比（数字化 78.7%，interpreted）。",
            "官方 CBL 结论'不合格'按 25m 连续胶结红线判据；占比较值仅作参考。",
            "设计版与默认实际版共用几何/流体/校验资料，差异仅在 PumpingSchedule。",
        ),
    )
    return well_spec, fluids, schedule, validation_data
