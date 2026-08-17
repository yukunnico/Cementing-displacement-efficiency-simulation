"""
呼探1-002井（HT1-002）139.7mm完井尾管段标准数据加载器

本模块把呼探1-002井现场数据包中的 139.7mm 完井尾管资料整理为 cemdisp
标准输入结构。与呼探1井不同，本井目标段为单一 139.7mm 尾管，不构造
双径尾管等效几何；上部地面至悬挂器的管内容积只用于鞋口滞后体积估算。

2026-08-16 更新（核对报告：其他井loader核对_2026-08-16.md，hu2/HT1-002 节）：
- 井径/井斜剖面从合成代理切换为现场提取包实测（64 点，30m 间隔，field_measured）；
  5292.5–5631m 为 219.1mm 技套内重叠段（技套鞋 5630m），取技套内径 193.7mm 等效。
- 隔离液密度 2.10 → 现场 2.05 g/cm³（提取包统一 2.05，field_measured；旧"实际 2.10"为历史口径）。
- 钻井液 YP 5 → 8.5 Pa（化验栏）；隔离液/压塞液/中置液流变按化验 PV30/YP8。
- 压塞液/中置液排量 0.6 → 0.8 m³/min（pumping_schedule，field_measured）。
- 居中度改照设计值 78%（斯伦贝谢模拟，无实测，model_assumption）。
- 目标层段 7400–7500 → 油气显示层 7402–7500（施工记录表 R3）。
- 旧代理值一律以 LEGACY(2026-08-16 前) 注释保留。
"""

from __future__ import annotations

import warnings

import csv
from collections.abc import Callable
import math
from pathlib import Path

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep
from cemdisp.data.validation_data import ValidationData
from cemdisp.data.well_spec import DepthValuePoint, EvaluationWindow, WellSpec
from cemdisp.models2d.boundary_bridge import AnnulusInletState
from cemdisp.data.provenance import WELL_PROVENANCE
from cemdisp.models2d.boundary_bridge import build_sync_card
from cemdisp.transport1d.casing_flow import CasingFlowSolver



# 项目根目录用于定位参考文档，避免在求解器中硬编码单井文件路径。
PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_REFERENCE_ROOT = PROJECT_ROOT / "参考文档" / "呼探1-002"
# 现场提取包井径/井斜 CSV：施工设计 Table7/Table8 电测数据（30m 间隔，64 点，field_measured）。
# LEGACY(2026-08-16 前): 合成剖面（平均井径 + 两段大肚子常量 + 井斜代理）。
DEFAULT_CALIPER_CSV = PROJECT_ROOT / "参考文档" / "现场资料提取" / "ht1_002_呼探1-002" / "caliper_profile.csv"
DEFAULT_INCLINATION_CSV = PROJECT_ROOT / "参考文档" / "现场资料提取" / "ht1_002_呼探1-002" / "inclination_profile.csv"

# 呼探1-002井段参数，来源：参考文档/呼探1-002/提取数据/呼探1-002井_固井顶替模型数据包.json +
# 现场提取包 ht1_002_呼探1-002/（2026-08-16 已核对，well_geometry.csv / casing_liner_string.csv）。
HU2_WELL_NAME = "呼探1-002井（HT1-002）"
HU2_WELL_DEPTH_MD_M = 7559.0  # 实际完钻井深，鞋以下留有约 5m 口袋段（well_basic_info，field_measured）。
HU2_HANGER_MD_M = 5292.5  # 139.7mm 完井尾管顶部/悬挂器位置（设计值，well_geometry，field_measured）。
HU2_TOP_MD_M = HU2_HANGER_MD_M  # 模型剖面从尾管悬挂器开始（simulation_top，model_assumption）。
HU2_BOTTOM_MD_M = 7554.0  # 139.7mm 尾管下深/固井井段底部（liner_shoe，field_measured）。
HU2_SHOE_MD_M = HU2_BOTTOM_MD_M  # 环空入口对应尾管鞋深度。
HU2_BIT_DIAMETER_MM = 190.5  # 五开钻头尺寸/名义井径（field_measured，仅作资料留痕）。
HU2_AVERAGE_HOLE_DIAMETER_MM = 193.05  # 目标尾管段平均井径（设计 1.4.2 口径；提取包 CSV 全段均值 194.71，field_measured）。
# 上部 219.1mm 技术套管（四开，3753–5630m）内径：219.1 − 2×12.7 = 193.7mm（casing_liner_string，field_measured）。
# 该值兼作悬挂器以上管内容积计算的内径代理，以及 5292.5–5631m 技套内重叠段环空外边界等效井径（model_assumption）。
# LEGACY(2026-08-16 前): 195.0（旧代理值，与现场 193.7 差 1.3mm）。
HU2_CASING_ID_MM = 193.7
HU2_LINER_OD_MM = 139.7  # 完井尾管外径，本井为单一口径（field_measured）。
HU2_LINER_WALL_THICKNESS_MM = 15.88  # 139.7mm BG140V/BG-T2 套管壁厚（施工记录表 R12，field_measured）。
HU2_LINER_ID_ACTUAL_MM = HU2_LINER_OD_MM - 2.0 * HU2_LINER_WALL_THICKNESS_MM  # 107.94mm，field_measured
HU2_CENTRALIZER_COUNT = 95  # φ139.7*190.5mm 整体式弹扶数量（field_measured，casing_liner_string）。
HU2_TARGET_CENTRALIZER_SPACING_M = 44.0  # 目的层扶正器间距。
HU2_NON_TARGET_CENTRALIZER_SPACING_M = 55.0  # 非目的层扶正器间距。
# 大肚子 1（5631–5940m、设计平均 206.65/最大 268.48）：来源=施工设计 1.4.2 节（设计口径，model_assumption）。
# 30m 采样提取包 CSV 段内均值 211.35，与设计值差 4.7mm；常量保留供口径追溯。
HU2_LARGE_HOLE_1_TOP_MD_M = 5631.0  # 第一段大肚子井眼起点。
HU2_LARGE_HOLE_1_BOTTOM_MD_M = 5940.0  # 第一段大肚子井眼终点。
HU2_LARGE_HOLE_1_AVERAGE_DIAMETER_MM = 206.65  # 第一段大肚子平均井径（设计 1.4.2）。
HU2_LARGE_HOLE_1_MAX_DIAMETER_MM = 268.48  # 第一段大肚子最大井径，作为资料备注。
# 大肚子 2（7355–7360m、设计平均 226.59/最大 249.90）：来源=施工设计 1.4.2 节（设计口径，model_assumption）。
# 30m 采样提取包 CSV 在 7355–7360 无测点（全段 190.5），无法复现该 5m 段，常量仅作设计备注。
HU2_LARGE_HOLE_2_TOP_MD_M = 7355.0  # 第二段大肚子井眼起点。
HU2_LARGE_HOLE_2_BOTTOM_MD_M = 7360.0  # 第二段大肚子井眼终点。
HU2_LARGE_HOLE_2_AVERAGE_DIAMETER_MM = 226.59  # 第二段大肚子平均井径（设计 1.4.2）。
HU2_LARGE_HOLE_2_MAX_DIAMETER_MM = 249.90  # 第二段大肚子最大井径，作为资料备注。
# 井斜极值：设计资料记为最大 3.9°@7443m（HU2_MAX_INCLINATION_* 为设计口径）；提取包 30m CSV
# 实测最大 3.5°@7440m（无 7443 测点），两处并存，仅作资料备注，剖面已改用逐点实测。
HU2_MAX_INCLINATION_DEG = 3.9  # 最大井斜（设计资料，model_assumption）。
HU2_MAX_INCLINATION_MD_M = 7443.0  # 最大井斜所在测深（设计资料）。


def _pipe_volume_m3(length_m: float, inner_diameter_mm: float) -> float:
    """按内径和长度估算管内容积，用于鞋口滞后体积计算。"""

    # 由毫米内径换算为米半径后计算圆管容积。
    radius_m = inner_diameter_mm / 1000.0 / 2.0
    return math.pi * radius_m**2 * length_m


# 鞋口滞后体积估算：地面至悬挂器按 219.1mm 技术套管内径代理（193.7），尾管段按 107.94mm 实际内径累加。
# 该体积只用于 field_order_realistic 边界的到鞋延迟，不代表环空求解器几何直径。
# 尾管段 5308.568–7052.26m 实为 14.27mm 薄壁段（ID 111.16），此处按厚壁段 ID 107.94 统一估算（简化，model_assumption）。
HU2_SURFACE_TO_HANGER_EFFECTIVE_ID_MM = HU2_CASING_ID_MM
HU2_SHOE_LAG_VOLUME_M3 = _pipe_volume_m3(
    HU2_HANGER_MD_M,
    HU2_SURFACE_TO_HANGER_EFFECTIVE_ID_MM,
) + _pipe_volume_m3(
    HU2_SHOE_MD_M - HU2_HANGER_MD_M,
    HU2_LINER_ID_ACTUAL_MM,
)
HU2_LINER_ID_MM = HU2_LINER_ID_ACTUAL_MM

# 呼探1-002现场流体参数（2026-08-16 按提取包 fluid_properties.csv/化验报告核对，field_measured 为主）；
# Bingham 参数按实验值或代理值，水泥浆按实验幂律参数。
HU2_MUD_DENSITY_KG_M3 = 2060.0  # 油基钻井液密度 2.06（field_measured）。
HU2_DISPLACEMENT_DENSITY_KG_M3 = 2060.0  # 替浆钻井液密度 2.06（替浆-1..5，pumping_schedule，field_measured）。
HU2_BALANCE_DENSITY_KG_M3 = 1850.0  # 先导浆密度 1.85（field_measured）。
HU2_SPACER_DENSITY_KG_M3 = 2050.0  # 驱油隔离液密度 2.05（field_measured，pumping_schedule/fluid_properties）。
# LEGACY(2026-08-16 前): 2100.0（旧 loader 记"实际 2.10"，与提取包 2.05 冲突，以提取包为准；
# 现场是否存在 2.10 口径需复核施工记录表 R20，本版不采信）。
HU2_LEAD_DENSITY_KG_M3 = 2100.0  # 领浆密度 2.10（field_measured，化验）。
HU2_INTERMEDIATE_DENSITY_KG_M3 = 1950.0  # 中间浆密度 1.95（field_measured，化验）。
HU2_TAIL_DENSITY_KG_M3 = 1950.0  # 尾浆密度 1.95（field_measured，化验）。
HU2_MUD_PV_PA_S = 0.058  # 钻井液 PV：58mPa·s（field_measured，化验 133→93℃）。
HU2_MUD_YP_PA = 8.5  # 钻井液 YP=8.5Pa（field_measured，化验；老化后 10.5、循环后 6）。
# LEGACY(2026-08-16 前): 5.0（终切力 5Pa 代理）
HU2_DISPLACEMENT_PV_PA_S = 0.058  # 替浆液按同密度钻井液流变处理（field_measured 同钻井液）。
HU2_DISPLACEMENT_YP_PA = 8.5  # 替浆液 YP=8.5Pa（与钻井液一致，field_measured）。
HU2_BALANCE_PV_PA_S = 0.030  # 平衡液缺流变实测，使用代理 PV（proxy）。
HU2_BALANCE_YP_PA = 3.0  # 平衡液缺流变实测，使用代理 YP（proxy）。
HU2_SPACER_PV_PA_S = 0.030  # 隔离液塑粘 30mPa·s（field_measured，化验 Table7）。
HU2_SPACER_YP_PA = 8.0  # 隔离液 YP=8Pa（field_measured，化验）；LEGACY(2026-08-16 前): 5.0
HU2_LEAD_POWER_LAW_N = 0.811  # 领浆 n（field_measured，化验 133℃）。
HU2_LEAD_CONSISTENCY_K = 0.876  # 领浆 K Pa·s^n（field_measured）。
HU2_INTERMEDIATE_POWER_LAW_N = 0.871  # 中间浆 n（field_measured）。
HU2_INTERMEDIATE_CONSISTENCY_K = 0.504  # 中间浆 K（field_measured）。
HU2_TAIL_POWER_LAW_N = 0.886  # 尾浆 n（field_measured）。
HU2_TAIL_CONSISTENCY_K = 0.453  # 尾浆 K（field_measured）。
HU2_PLUG_DENSITY_KG_M3 = 1900.0  # 压塞液密度 1.90（field_measured，施工记录表 R30）。
HU2_PLUG_PV_PA_S = 0.030  # 压塞液塑粘 30mPa·s（隔离液类型，field_measured）。
HU2_PLUG_YP_PA = 8.0  # 压塞液 YP=8Pa（隔离液类型，field_measured）。
HU2_MIDDLE_FLUID_DENSITY_KG_M3 = 1900.0  # 中置液密度 1.90（field_measured，施工记录表 R30）。
HU2_MIDDLE_FLUID_PV_PA_S = 0.030  # 中置液塑粘 30mPa·s（隔离液类型，field_measured）。
HU2_MIDDLE_FLUID_YP_PA = 8.0  # 中置液 YP=8Pa（隔离液类型，field_measured）。

# 呼探1-002现场施工程序参数，按地面注入顺序排列（pumping_schedule.csv，field_measured）。
HU2_BALANCE_VOLUME_M3 = 20.0  # 先导浆 20m³@0.8。
HU2_SPACER_VOLUME_M3 = 15.0  # 驱油隔离液 15m³@0.8。
HU2_LEAD_VOLUME_M3 = 12.0  # 领浆 12m³@0.8。
HU2_INTERMEDIATE_VOLUME_M3 = 14.0  # 中间浆 14m³@0.8。
HU2_TAIL_VOLUME_M3 = 12.0  # 尾浆 12m³@0.8。
HU2_PLUG_VOLUME_M3 = 5.0  # 压塞液 5m³（Table37 排量 0.8）。
HU2_FIRST_DISPLACEMENT_VOLUME_M3 = 12.0  # 替浆-1 钻井液 12m³@0.8。
HU2_MIDDLE_DISPLACEMENT_VOLUME_M3 = 15.0  # 注中置液 15m³@0.8。
# 重泥浆替浆 47m³ 现场四段：替浆-2 15@0.8 + 替浆-3 5@0.5 + 替浆-4 18@0.3 + 替浆-5 9@0.4；
# 模型合并为快替 30@0.8 + 慢替 17@0.3（总 47 一致，分段为模型代理，LEGACY(2026-08-16 前) 同）。
HU2_FAST_DISPLACEMENT_VOLUME_M3 = 30.0  # 现场 15+5 合并（0.5 段并入快替，模型代理）。
HU2_SLOW_DISPLACEMENT_VOLUME_M3 = 17.0  # 现场 18+9 合并（0.4 段并入慢替，模型代理）。
HU2_MAIN_RATE_M3_MIN = 0.8  # 主排量（先导浆/隔离液/水泥浆/替浆-1，field_measured）。
HU2_PLUG_RATE_M3_MIN = 0.8  # 压塞液排量 0.8（pumping_schedule Table37，field_measured）。
# LEGACY(2026-08-16 前): 0.6（旧值，与现场 0.8 不符）
HU2_MIDDLE_DISPLACEMENT_RATE_M3_MIN = 0.8  # 中置液排量 0.8（pumping_schedule step10 @0.8，field_measured）。
# LEGACY(2026-08-16 前): 0.6（旧值，loader notes 曾自注"现场 0.8→0.5，模型用 0.6 代理"）
HU2_FAST_DISPLACEMENT_RATE_M3_MIN = 0.8  # 快替 0.8（模型代理，现场 0.8/0.5 合并段）。
HU2_SLOW_DISPLACEMENT_RATE_M3_MIN = 0.3  # 慢替 0.3（模型代理，现场 0.3/0.4 合并段）。


def _depth_points(values: tuple[tuple[float, float], ...]) -> tuple[DepthValuePoint, ...]:
    """把测深-数值元组转换为 WellSpec 使用的剖面点。"""

    # WellSpec 使用不可变剖面点，便于后续求解器按深度插值。
    return tuple(DepthValuePoint(depth_md_m=depth, value=value) for depth, value in values)


def _read_caliper_rows(caliper_csv_path: Path) -> tuple[tuple[float, float], ...]:
    """读取现场提取包井径剖面 CSV（md_m / caliper_mm），按 md_m 升序返回。"""
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
    """从现场提取包井径行构建井径剖面（裸眼实测 caliper_mm，field_measured）。

    模型段顶 5292.5–5631m 为 219.1mm 技套内重叠段（技套鞋 5630m），无裸眼井径测点，
    取技套内径 193.7mm 作为环空外边界等效井径（model_assumption，与旧平均井径 193.05 接近）。
    """
    points: list[tuple[float, float]] = [(md, cal) for md, cal in caliper_rows if md >= HU2_TOP_MD_M]
    if not points or points[0][0] > HU2_TOP_MD_M:
        points.insert(0, (HU2_TOP_MD_M, HU2_CASING_ID_MM))
    if points[-1][0] < HU2_BOTTOM_MD_M:
        points.append((HU2_BOTTOM_MD_M, points[-1][1]))
    return tuple(points)


def _build_inclination_profile(incl_rows: tuple[tuple[float, float], ...]) -> tuple[tuple[float, float], ...]:
    """从现场提取包井斜行构建井斜剖面（裸眼实测 inclination_deg，field_measured）。

    模型段顶 5292.5–5631m 无测斜点，取首测点值（5631m，0.1°）近垂直外推（proxy）。
    """
    points: list[tuple[float, float]] = [(md, inc) for md, inc in incl_rows if md >= HU2_TOP_MD_M]
    if not points or points[0][0] > HU2_TOP_MD_M:
        points.insert(0, (HU2_TOP_MD_M, points[0][1]))
    if points[-1][0] < HU2_BOTTOM_MD_M:
        points.append((HU2_BOTTOM_MD_M, points[-1][1]))
    return tuple(points)


def _phase_fractions_for_role(role: FluidRole, *, split_cement_phases: bool) -> tuple[tuple[str, float], ...]:
    """把标准流体角色映射为环空二维模型相名称。"""

    # 分相口径下：领浆和中间浆合并为 lead，尾浆单独为 tail。
    if split_cement_phases and role in {FluidRole.LEAD, FluidRole.INTERMEDIATE}:
        return (("lead", 1.0),)
    if split_cement_phases and role == FluidRole.TAIL:
        return (("tail", 1.0),)
    # 默认三相口径下：所有水泥浆统一归入 cement。
    if role in {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}:
        return (("cement", 1.0),)
    if role in {FluidRole.WASH, FluidRole.SPACER}:
        return (("spacer", 1.0),)
    return (("mud", 1.0),)


def _phase_fractions_for_fluid(
    fluid_name: str,
    role_by_name: dict[str, FluidRole],
    *,
    split_cement_phases: bool,
) -> tuple[tuple[str, float], ...]:
    """按流体名称映射入口相分数，支持 lead/tail 分相。"""

    # 未识别流体保守映射为 mud，避免把管内压塞液误认为水泥相。
    role = role_by_name.get(fluid_name, FluidRole.MUD)
    return _phase_fractions_for_role(role, split_cement_phases=split_cement_phases)


def load_hu2_tailpipe(
    *,
    reference_root: Path | None = None,
) -> tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载呼探1-002井（HT1-002）139.7mm完井尾管段标准模型输入。"""

    # 允许调用方传入替代资料根目录，默认使用项目内参考文档路径。
    resolved_reference_root = reference_root or DEFAULT_REFERENCE_ROOT
    caliper_rows = _read_caliper_rows(DEFAULT_CALIPER_CSV)
    incl_rows = _read_inclination_rows(DEFAULT_INCLINATION_CSV)
    well_spec = WellSpec(
        well_name=HU2_WELL_NAME,
        top_md_m=HU2_TOP_MD_M,
        bottom_md_m=HU2_BOTTOM_MD_M,
        shoe_md_m=HU2_SHOE_MD_M,
        hanger_md_m=HU2_HANGER_MD_M,
        casing_id_mm=HU2_CASING_ID_MM,
        liner_od_mm=HU2_LINER_OD_MM,
        liner_id_mm=HU2_LINER_ID_MM,
        hole_diameter_profile=_depth_points(_build_hole_profile(caliper_rows)),
        inclination_profile=_depth_points(_build_inclination_profile(incl_rows)),
        standoff_profile=_depth_points(
            (
                # 居中度按设计值 78%（斯伦贝谢模拟，centralization_profile.csv 仅记两段 78%，
                # 无实测连续分布，model_assumption）；鞋底略降。
                # LEGACY(2026-08-16 前): 0.58–0.75 扶正器代理剖面（低于设计 78%）。
                (HU2_TOP_MD_M, 0.78),
                (5900.0, 0.78),
                (6500.0, 0.78),
                (7000.0, 0.78),
                (7350.0, 0.78),
                (7400.0, 0.78),
                (7500.0, 0.78),
                (HU2_SHOE_MD_M, 0.76),
            )
        ),
        evaluation_windows=(
            # CBL评价井段覆盖完整 139.7mm 完井尾管固井井段（与施工设计一致）。
            EvaluationWindow(name="CBL评价井段", top_md_m=HU2_TOP_MD_M, bottom_md_m=HU2_BOTTOM_MD_M, window_type="cbl"),
            # 目标层段按施工记录表 R3 油气显示层 7402–7500m；LEGACY(2026-08-16 前): 7400–7500、type="target"。
            EvaluationWindow(name="目标层段（油气显示层）", top_md_m=7402.0, bottom_md_m=7500.0, window_type="formation_target"),
            # 地层目标（target_intervals.csv，施工设计 Table0/5.1 节，field_measured）。
            EvaluationWindow(name="目的层-K1q清水河组(地层目标)", top_md_m=5702.0, bottom_md_m=6826.9, window_type="formation_target"),
            EvaluationWindow(name="目的层-J3k2喀拉扎组(地层目标)", top_md_m=6826.9, bottom_md_m=7554.0, window_type="formation_target"),
            EvaluationWindow(name="高压水层(地层目标)", top_md_m=5998.0, bottom_md_m=6002.0, window_type="formation_target"),
        ),
        reference_root=resolved_reference_root,
        notes=(
            "呼探1-002目标段为单一 139.7mm 完井尾管，不采用呼探1双径尾管等效几何。",
            "井径/井斜剖面改用现场提取包实测（各 64 点，30m 间隔，field_measured）；5292.5–5631m 为 219.1mm 技套内重叠段，取技套内径 193.7mm 等效（model_assumption）。",
            "大肚子 1（5631–5940m，设计平均 206.65）与大肚子 2（7355–7360m，设计平均 226.59）均来自施工设计 1.4.2 节（设计口径，model_assumption）；30m 采样 CSV 无法复现大肚子 2（7355–7360 无测点）。",
            "居中度按设计值 78%（斯伦贝谢模拟，无实测）构造剖面，model_assumption。",
            f"鞋口滞后体积估算为 {HU2_SHOE_LAG_VOLUME_M3:.2f}m³：地面至悬挂器按 193.7mm（219.1 技套内径）代理内径，尾管段按 {HU2_LINER_ID_ACTUAL_MM:.2f}mm 内径。",
            "CBL评价井段为 5292.5–7554m 完整尾管段；质量标签仅为'合格'弱监督代理（现场无数值 CBL 合格率），不等同顶替效率真值。",
        ),
    )

    # 九类环空模型核心流体；压塞液/中置液流变按隔离液类型化验值（PV30/YP8，field_measured）。
    fluids = (
        FluidSpec("钻井液", FluidRole.MUD, HU2_MUD_DENSITY_KG_M3, RheologyModel.BINGHAM, HU2_MUD_PV_PA_S, HU2_MUD_YP_PA),
        FluidSpec("替浆液", FluidRole.DISPLACEMENT, HU2_DISPLACEMENT_DENSITY_KG_M3, RheologyModel.BINGHAM, HU2_DISPLACEMENT_PV_PA_S, HU2_DISPLACEMENT_YP_PA),
        FluidSpec("平衡液", FluidRole.WASH, HU2_BALANCE_DENSITY_KG_M3, RheologyModel.BINGHAM, HU2_BALANCE_PV_PA_S, HU2_BALANCE_YP_PA),
        FluidSpec("隔离液", FluidRole.SPACER, HU2_SPACER_DENSITY_KG_M3, RheologyModel.BINGHAM, HU2_SPACER_PV_PA_S, HU2_SPACER_YP_PA),
        FluidSpec("领浆", FluidRole.LEAD, HU2_LEAD_DENSITY_KG_M3, RheologyModel.POWER_LAW, power_law_n=HU2_LEAD_POWER_LAW_N, consistency_k=HU2_LEAD_CONSISTENCY_K),
        FluidSpec("中间浆", FluidRole.INTERMEDIATE, HU2_INTERMEDIATE_DENSITY_KG_M3, RheologyModel.POWER_LAW, power_law_n=HU2_INTERMEDIATE_POWER_LAW_N, consistency_k=HU2_INTERMEDIATE_CONSISTENCY_K),
        FluidSpec("尾浆", FluidRole.TAIL, HU2_TAIL_DENSITY_KG_M3, RheologyModel.POWER_LAW, power_law_n=HU2_TAIL_POWER_LAW_N, consistency_k=HU2_TAIL_CONSISTENCY_K),
        FluidSpec("压塞液", FluidRole.OTHER, HU2_PLUG_DENSITY_KG_M3, RheologyModel.BINGHAM, HU2_PLUG_PV_PA_S, HU2_PLUG_YP_PA),
        FluidSpec("中置液", FluidRole.DISPLACEMENT, HU2_MIDDLE_FLUID_DENSITY_KG_M3, RheologyModel.BINGHAM, HU2_MIDDLE_FLUID_PV_PA_S, HU2_MIDDLE_FLUID_YP_PA),
        FluidSpec("井浆", FluidRole.DISPLACEMENT, HU2_DISPLACEMENT_DENSITY_KG_M3, RheologyModel.BINGHAM, HU2_DISPLACEMENT_PV_PA_S, HU2_DISPLACEMENT_YP_PA),
    )

    # 地面施工程序按现场注入顺序排列（pumping_schedule.csv，field_measured）；
    # 47m³ 重泥浆替浆现场四段（15@0.8 + 5@0.5 + 18@0.3 + 9@0.4）拆分为快替 30m³ 与慢替 17m³（模型代理，总量 47 一致）。
    schedule = PumpingSchedule(
        steps=(
            PumpingScheduleStep("注入平衡液", "平衡液", HU2_BALANCE_VOLUME_M3, HU2_MAIN_RATE_M3_MIN, remarks="先导浆 20m³@0.8m³/min，角色 WASH（field_measured）。"),
            PumpingScheduleStep("注入隔离液", "隔离液", HU2_SPACER_VOLUME_M3, HU2_MAIN_RATE_M3_MIN, remarks="驱油隔离液 15m³@0.8m³/min，密度 2.05g/cm³（提取包 field_measured；旧记'实际 2.10'为历史口径）。"),
            PumpingScheduleStep("注入领浆", "领浆", HU2_LEAD_VOLUME_M3, HU2_MAIN_RATE_M3_MIN, remarks="领浆 12m³@0.8m³/min，ρ2.10（field_measured）。"),
            PumpingScheduleStep("注入中间浆", "中间浆", HU2_INTERMEDIATE_VOLUME_M3, HU2_MAIN_RATE_M3_MIN, remarks="中间浆 14m³@0.8m³/min，ρ1.95，角色 INTERMEDIATE（field_measured）。"),
            PumpingScheduleStep("注入尾浆", "尾浆", HU2_TAIL_VOLUME_M3, HU2_MAIN_RATE_M3_MIN, remarks="尾浆 12m³@0.8m³/min，ρ1.95（field_measured）。"),
            PumpingScheduleStep("注入压塞液（管内）", "压塞液", HU2_PLUG_VOLUME_M3, HU2_PLUG_RATE_M3_MIN, remarks="压塞液 5m³@0.8m³/min（现场 Table37），仅作为管内占位（field_measured）。"),
            PumpingScheduleStep("替浆泥浆(快)", "井浆", HU2_FIRST_DISPLACEMENT_VOLUME_M3, HU2_MAIN_RATE_M3_MIN, remarks="替浆-1 钻井液 12m³@0.8m³/min，角色 DISPLACEMENT（field_measured）。"),
            PumpingScheduleStep("替浆中置液", "中置液", HU2_MIDDLE_DISPLACEMENT_VOLUME_M3, HU2_MIDDLE_DISPLACEMENT_RATE_M3_MIN, remarks="注中置液 15m³@0.8m³/min，ρ1.90（field_measured）。"),
            PumpingScheduleStep("井浆快替", "井浆", HU2_FAST_DISPLACEMENT_VOLUME_M3, HU2_FAST_DISPLACEMENT_RATE_M3_MIN, remarks="重泥浆 47m³ 中快替 30m³@0.8m³/min（现场 15@0.8+5@0.5 合并，模型代理）。"),
            PumpingScheduleStep("井浆慢替", "井浆", HU2_SLOW_DISPLACEMENT_VOLUME_M3, HU2_SLOW_DISPLACEMENT_RATE_M3_MIN, remarks="重泥浆 47m³ 中慢替 17m³@0.3m³/min（现场 18@0.3+9@0.4 合并，模型代理）。"),
        ),
        notes=(
            "施工顺序按现场记录：先导浆→驱油隔离液→领浆→中间浆→尾浆→压塞液→替浆-1→中置液→替浆-2/3/4/5（现场四段合并为快/慢替）。",
            "替浆总量 79m³ = 5 + 12 + 15 + 47m³，其中压塞液用于管内时序占位（field_measured）。",
            "隔离液密度以提取包 2.05g/cm³ 为准（field_measured）；旧 loader 记'实际 2.10'与提取包冲突，属历史口径，不再混用。",
            "替浆-2/3/4/5 现场四段为 15@0.8 + 5@0.5 + 18@0.3 + 9@0.4（pumping_schedule，field_measured）；模型合并为 30@0.8 + 17@0.3，总量 47 一致。",
        ),
    )

    # 校验资料路径指向数据包 JSON；本加载器不在运行时解析 JSON，只固化标准输入常量。
    validation_data = ValidationData(
        job_report_path=resolved_reference_root / "提取数据" / "呼探1-002井_固井顶替模型数据包.json",
        notes=(
            "呼探1-002首版加载器不读取 JSON 文件，仅把数据包中的现场参数固化为模块常量。",
            "CBL 定量合格率缺失：三份'报告'PDF 经视觉核查均为水泥浆实验报告而非 CBL 评价图，仅定性'合格'，cbl_pass_rate 保持 None。",
            "数据标签口径：井径/井斜/密度/流变/泵注为 field_measured（提取包）；等效内径、大肚子设计口径、居中度 78%、技套内重叠段等效井径为 model_assumption。",
        ),
    )
    return well_spec, fluids, schedule, validation_data


def export_hu2_sync_card_markdown(output_dir: Path) -> Path:
    """导出 呼探1-002井（HT1-002） 同步画像卡为 Markdown 文件。"""

    well_spec, fluids, schedule, _ = load_hu2_tailpipe()
    solver = CasingFlowSolver()
    result = solver.run(well_spec, fluids, schedule)
    provenance = WELL_PROVENANCE.get("呼探1-002井（HT1-002）")
    if provenance is None:
        from cemdisp.data.provenance import _WELL_NAME_ALIASES
        canonical = _WELL_NAME_ALIASES.get("呼探1-002井（HT1-002）", "呼探1-002井（HT1-002）")
        provenance = WELL_PROVENANCE[canonical]
    sync_card = build_sync_card("呼探1-002井（HT1-002）", result.shoe_timeline, provenance)

    output_path = output_dir / ("呼探1-002井（HT1-002）_同步画像卡.md")
    lines = [
        "# 呼探1-002井（HT1-002） 同步画像卡",
        "",
        "- 井名：" + str(sync_card["井名"]),
        "- 鞋口同步事件数：" + str(sync_card["鞋口同步口径"]["事件数"]),
    ]
    first_time = sync_card["鞋口同步口径"]["首事件时间_s"]
    last_time = sync_card["鞋口同步口径"]["末事件时间_s"]
    if first_time is not None:
        lines.append("- 首事件时间：{:.1f} s".format(first_time))
    if last_time is not None:
        lines.append("- 末事件时间：{:.1f} s".format(last_time))
    proxy_note = sync_card["代理提醒"]
    if proxy_note:
        lines.append("- 代理提醒：" + proxy_note)
    output_path.write_text("\n".join(lines), encoding="utf-8")
    return output_path