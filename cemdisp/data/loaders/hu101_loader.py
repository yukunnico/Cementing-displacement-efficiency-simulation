"""
呼101尾管段标准数据加载器

本模块把 legacy 呼101脚本中的本井资料整理为 cemdisp 标准输入结构，
用于 Hu101 runner 复用通用 2D 环空求解、导出和 1D-2D 边界桥接流程。

2026-08-16 更新（核对报告：其他井loader核对_2026-08-16.md，hu101 节）：
- 井径/井斜剖面从合成剖面切换为现场提取包实测（108 点，20m 间隔，field_measured）；
  上部 168.3mm 段（5400–6796.329m）实测井径按 139.7mm 参考外径做保面积等效（model_assumption），
  以适配求解器单一 139.7mm 外径几何；下部 139.7mm 段直接用实测值。
- liner_id_mm=91.73（52m³ 鞋口滞后体积反推）显式标 model_assumption：非实测单一内径，
  为复合尾管（168.3 上段 + 139.7 下段，厚壁 ID 108.1/薄壁 ID 111.16）的等效几何口径。
- standoff 0.38–0.48 剖面显式标 model_assumption：现场仅有扶正器布置（132 只）无实测居中度，
  且悬挂器坐挂失败、最终座底固井，居中度实际更低。
- cbl_pass_rate=0.6277 结构化进 ValidationData（参照 hu102=0.6665 / hu103=0.1206 写法），
  注明正式解释测量段 5390–7810m 口径；CBL 评价窗保持 legacy 的 5700–7810m（裸眼段顶+测量段底组合）。
- 中置液排量 1.2 → 1.5 m³/min（现场 pumping_schedule，field_measured）。
- casing_id_mm 命名语义问题（实存 OD 273.1）仅保留/强化注释，不改名（全仓统一需动 API，超出本任务）。
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



PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_REFERENCE_ROOT = PROJECT_ROOT / "参考文档" / "呼101"
# 现场提取包井径/井斜 CSV：固井作业史 2011113.doc 电测数据（108 点，20m 间隔，field_measured）。
# LEGACY(2026-08-16 前): 合成剖面（上段等效 242.84 + 下段 215.9 常量；井斜 0.6→1.9°）。
DEFAULT_CALIPER_CSV = PROJECT_ROOT / "参考文档" / "现场资料提取" / "hu101_呼101" / "caliper_profile.csv"
DEFAULT_INCLINATION_CSV = PROJECT_ROOT / "参考文档" / "现场资料提取" / "hu101_呼101" / "inclination_profile.csv"

# 呼101尾管段井身结构参数（2026-08-16 按现场提取包 well_geometry/casing_liner_string/100312.PDF 核对）。
HU101_TOP_MD_M = 5400.0  # 模型段顶：水泥返高 5400m 口径（高于尾管串顶 5402.885m 约 2.9m，notes 注明）。
HU101_BOTTOM_MD_M = 7868.0  # 完钻/尾管鞋深度（浮鞋下深，field_measured）。
HU101_SHOE_MD_M = 7868.0
HU101_HANGER_MD_M = 5407.46  # 悬挂器本体下深（field_measured）；悬挂器跨度 5402.85–5407.46m。
# 字段名 casing_id_mm 实际上存 273.1mm 技术套管"外径"（well_geometry；技套鞋 5699.8m）。
# 命名语义历史遗留（全仓 hu101/hu102/hu103/ht1_003 均同），统一改名需动 API，超本任务范围，仅注明不改名。
HU101_TECH_CASING_EQUIV_ID_MM = 273.10  # 273.1mm 技术套管外径（OD，非 ID）。
HU101_UPPER_SECTION_BOTTOM_MD_M = 6796.0  # 168.3mm 段底（变径变扣 6796.329m，差 0.329m 取整）；LEGACY 同。
HU101_UPPER_ACTUAL_HOLE_DIAMETER_MM = 260.35  # 上段实际井径（设计口径；实测 5700–6796 段均值 259.46）。
HU101_UPPER_ACTUAL_LINER_OD_MM = 168.30  # 上段尾管外径（5402.885–6796.136m，field_measured）。
HU101_LOWER_HOLE_DIAMETER_MM = 215.90  # 下段井眼名义/实测井径（139.7mm 段，7048–7868 段均值 216.24）。
HU101_LOWER_LINER_OD_MM = 139.70  # 下段尾管外径（6796.329–7868m，field_measured）。
HU101_LINER_WALL_THICKNESS_MM = 15.80  # 139.7mm 段厚壁段壁厚（ID 108.1）；上部薄壁段 14.27mm（ID 111.16）存在。
HU101_SHOE_LAG_VOLUME_M3 = 52.0  # 鞋口滞后体积 52m³（现场记录，model_assumption 口径）。


def _equivalent_hole_diameter_mm(actual_hole_mm: float, actual_od_mm: float, reference_od_mm: float) -> float:
    """在固定 reference_od 条件下，构造与原始环空面积等价的井眼直径。"""

    area_term = actual_hole_mm**2 - actual_od_mm**2 + reference_od_mm**2
    return math.sqrt(max(area_term, reference_od_mm**2))


HU101_UPPER_HOLE_DIAMETER_MM = _equivalent_hole_diameter_mm(
    actual_hole_mm=HU101_UPPER_ACTUAL_HOLE_DIAMETER_MM,
    actual_od_mm=HU101_UPPER_ACTUAL_LINER_OD_MM,
    reference_od_mm=HU101_LOWER_LINER_OD_MM,
)
# HU101_LINER_ID_MM：由 52m³ 鞋口滞后体积 / 7868m 反推的等效内径（约 91.73mm）。
# **model_assumption**：非实测单一内径——139.7mm 段真实内径为厚壁 108.1mm / 薄壁 111.16mm，
# 该等效值仅服务 1D 到鞋时间口径（管内容积），论文不可写成现场内径。
HU101_LINER_ID_MM = math.sqrt(4.0 * HU101_SHOE_LAG_VOLUME_M3 / (math.pi * HU101_SHOE_MD_M)) * 1000.0

# 呼101现场施工与流体参数（2026-08-16 按提取包 fluid_properties.csv/pumping_schedule.csv 核对，
# 密度/体积/流变均 field_measured；此处不逐项改变）。
HU101_LEAD_VOLUME_M3 = 47.0  # 领浆 47m³@1.2，ρ2.10（field_measured，实际灰量 121t）。
HU101_TAIL_VOLUME_M3 = 23.0  # 尾浆 23m³@1.2，ρ1.90（field_measured）。
HU101_MUD_DENSITY_KG_M3 = 1960.0  # 油基钻井液 ρ1.96（field_measured，65℃）。
HU101_BALANCE_DENSITY_KG_M3 = 1850.0  # 平衡液 ρ1.85（field_measured）。
HU101_SPACER_DENSITY_KG_M3 = 2000.0  # 驱油隔离液 ρ2.00（field_measured）。
HU101_LEAD_DENSITY_KG_M3 = 2100.0  # 领浆 ρ2.10（field_measured，化验）。
HU101_TAIL_DENSITY_KG_M3 = 1900.0  # 尾浆 ρ1.90（field_measured，化验）。
HU101_MUD_PV_PA_S = 0.058  # 钻井液 65℃ PV=58mPa·s（field_measured）。
HU101_MUD_YP_PA = 9.2  # 钻井液 65℃ YP=9.2 Pa（field_measured，化验；施工前洗井 YP=5 为施工前口径）。
HU101_BALANCE_PV_PA_S = 0.030  # 平衡液缺流变实测，代理 PV（proxy）。
HU101_BALANCE_YP_PA = 3.0  # 平衡液缺流变实测，代理 YP（proxy）。
HU101_SPACER_PV_PA_S = 0.030  # 驱油隔离液塑粘 30mPa·s（field_measured）；YP=5 为代理（proxy）。
HU101_SPACER_YP_PA = 5.0
HU101_LEAD_POWER_LAW_N = 0.719  # 领浆 n（field_measured，化验 93℃）。
HU101_LEAD_CONSISTENCY_K = 0.815  # 领浆 K Pa·s^n（field_measured）。
HU101_TAIL_POWER_LAW_N = 0.722  # 尾浆 n（field_measured，化验 93℃）。
HU101_TAIL_CONSISTENCY_K = 0.684  # 尾浆 K Pa·s^n（field_measured）。


def _depth_points(values: tuple[tuple[float, float], ...]) -> tuple[DepthValuePoint, ...]:
    """把测深-数值元组转换为 WellSpec 使用的剖面点。"""

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

    求解器为单一 139.7mm 外径几何：上部 168.3mm 段（<6796.329m）实测井径按
    139.7mm 参考外径转换为保面积等效（model_assumption）；下部 139.7mm 段直接用实测值。
    模型段顶 5400–5720m 为 273.1mm 技套内重叠段（技套鞋 5699.8m），无裸眼测点，
    取旧等效上部井径 242.84mm 外推（model_assumption，LEGACY 常量）。
    """
    converted: list[tuple[float, float]] = []
    for md, cal in caliper_rows:
        if md < HU101_UPPER_SECTION_BOTTOM_MD_M:
            # 上部 168.3mm 段：保面积等效到 139.7mm 参考外径。
            converted.append((md, _equivalent_hole_diameter_mm(cal, HU101_UPPER_ACTUAL_LINER_OD_MM, HU101_LOWER_LINER_OD_MM)))
        else:
            converted.append((md, cal))
    points = [p for p in converted if p[0] >= HU101_TOP_MD_M]
    if not points or points[0][0] > HU101_TOP_MD_M:
        points.insert(0, (HU101_TOP_MD_M, HU101_UPPER_HOLE_DIAMETER_MM))
    if points[-1][0] < HU101_BOTTOM_MD_M:
        points.append((HU101_BOTTOM_MD_M, points[-1][1]))
    return tuple(points)


def _build_inclination_profile(incl_rows: tuple[tuple[float, float], ...]) -> tuple[tuple[float, float], ...]:
    """从现场提取包井斜行构建井斜剖面（裸眼实测 inclination_deg，field_measured）。

    模型段顶 5400–5720m 无测斜点，取首测点值（5720m，0.64°）外推（proxy）。
    """
    points: list[tuple[float, float]] = [(md, inc) for md, inc in incl_rows if md >= HU101_TOP_MD_M]
    if not points or points[0][0] > HU101_TOP_MD_M:
        points.insert(0, (HU101_TOP_MD_M, points[0][1]))
    if points[-1][0] < HU101_BOTTOM_MD_M:
        points.append((HU101_BOTTOM_MD_M, points[-1][1]))
    return tuple(points)


def _phase_fractions_for_role(role: FluidRole) -> tuple[tuple[str, float], ...]:
    """把标准流体角色映射为环空二维模型的三相名称。"""

    # 领浆/尾浆同属水泥相；平衡液/隔离液同属隔离液相；其余视为泥浆相。
    if role in {FluidRole.LEAD, FluidRole.TAIL}:
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
    """按流体名称映射入口相分数，支持 Hu101 的 lead/tail 分相。"""

    role = role_by_name.get(fluid_name, FluidRole.MUD)
    if split_cement_phases and role == FluidRole.LEAD:
        return (("lead", 1.0),)
    if split_cement_phases and role == FluidRole.TAIL:
        return (("tail", 1.0),)
    return _phase_fractions_for_role(role)


# 实测居中度剖面（从呼101尾管居中度检测图 Pipe Standoff 读取）。
# 扶正器之间偏下限、扶正器处偏上限；来源与 scripts/hu101_standoff_measured_vs_assumed.py 一致。
_MEASURED_STANDOFF_BETWEEN = (
    (HU101_TOP_MD_M, 0.78), (5700.0, 0.72), (6000.0, 0.65), (6300.0, 0.58),
    (6600.0, 0.52), (6900.0, 0.48), (7200.0, 0.42), (7500.0, 0.32),
    (7700.0, 0.25), (HU101_BOTTOM_MD_M, 0.22),
)
_MEASURED_STANDOFF_AT = (
    (HU101_TOP_MD_M, 0.88), (5700.0, 0.85), (6000.0, 0.80), (6300.0, 0.76),
    (6600.0, 0.72), (6900.0, 0.70), (7200.0, 0.68), (7500.0, 0.65),
    (7700.0, 0.62), (HU101_BOTTOM_MD_M, 0.60),
)
_ASSUMED_STANDOFF = (
    (HU101_TOP_MD_M, 0.45), (6100.0, 0.38), (6796.0, 0.44),
    (7200.0, 0.48), (7600.0, 0.42), (HU101_BOTTOM_MD_M, 0.46),
)


def _resolve_standoff_profile(measured_standoff):
    """返回 (md, standoff) 点序列：None=名义剖面，'between_centralizers'/'at_centralizers'=实测。"""
    if measured_standoff is None:
        return _ASSUMED_STANDOFF
    if measured_standoff == "between_centralizers":
        return _MEASURED_STANDOFF_BETWEEN
    if measured_standoff == "at_centralizers":
        return _MEASURED_STANDOFF_AT
    raise ValueError(
        f"measured_standoff 必须为 None/between_centralizers/at_centralizers，得到 {measured_standoff!r}")


def load_hu101_tailpipe(
    *,
    reference_root: Path | None = None,
    measured_standoff: str | None = None,
) -> tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载呼101尾管段标准模型输入。

    measured_standoff: 若为 "between_centralizers" 或 "at_centralizers"，用从
        居中度检测图读取的实测剖面（扶正器间/扶正器处）替换 model_assumption 的
        0.38–0.48 名义剖面；None（默认）保持名义剖面。
    """

    resolved_reference_root = reference_root or DEFAULT_REFERENCE_ROOT
    caliper_rows = _read_caliper_rows(DEFAULT_CALIPER_CSV)
    incl_rows = _read_inclination_rows(DEFAULT_INCLINATION_CSV)
    standoff_points = _resolve_standoff_profile(measured_standoff)
    well_spec = WellSpec(
        well_name="呼101",
        top_md_m=HU101_TOP_MD_M,
        bottom_md_m=HU101_BOTTOM_MD_M,
        shoe_md_m=HU101_SHOE_MD_M,
        hanger_md_m=HU101_HANGER_MD_M,
        casing_id_mm=HU101_TECH_CASING_EQUIV_ID_MM,
        liner_od_mm=HU101_LOWER_LINER_OD_MM,
        liner_id_mm=HU101_LINER_ID_MM,
        hole_diameter_profile=_depth_points(_build_hole_profile(caliper_rows)),
        inclination_profile=_depth_points(_build_inclination_profile(incl_rows)),
        standoff_profile=_depth_points(standoff_points),
        evaluation_windows=(
            # 正式 CBL 解释测量段 5390–7810m（100312.PDF，cbl_pass_rate=0.6277 对应整测量段口径，field_measured）；
            # 其中 5390–5699.8m 为双层套管段不评价（且在模型域水泥返高 5400m 之上），可评价段为 5699.8–7810m，
            # 7810–7868m 悬空段（尾管鞋未测 58m）无 CBL 数据。
            EvaluationWindow(name="CBL评价井段(单层套管可评价段)", top_md_m=5699.8, bottom_md_m=7810.0, window_type="cbl"),
            # 初评表：7537–7674m 大面积连续差段（J3k 喀拉扎组+漏层段，interpreted）。
            EvaluationWindow(name="CBL质量段(连续差段)", top_md_m=7537.0, bottom_md_m=7674.0, window_type="cbl_quality"),
            # 地层目标（target_intervals.csv，井史 2011111.doc 2.2 录井解释，field_measured）。
            EvaluationWindow(name="主要气层-K1q(地层目标)", top_md_m=7492.0, bottom_md_m=7536.0, window_type="formation_target"),
            EvaluationWindow(name="气层段(CBL评价重点)", top_md_m=7492.0, bottom_md_m=7735.0, window_type="formation_target"),
            EvaluationWindow(name="目的层-喀拉扎组J3k(地层目标)", top_md_m=7537.0, bottom_md_m=7868.0, window_type="formation_target"),
            EvaluationWindow(name="水层-K1s(地层目标)", top_md_m=6152.0, bottom_md_m=6156.0, window_type="formation_target"),
        ),
        reference_root=resolved_reference_root,
        notes=(
            "呼101上部168.3mm+下部139.7mm复合尾管（变径变扣 6796.329m）；求解器为单一 139.7mm 外径几何，"
            "上部实测井径按保面积转换为等效（model_assumption），下部直接使用实测值。",
            "liner_id_mm 使用 52m³ 鞋口滞后反推的等效内径（约 91.73mm），**model_assumption**：非实测单一内径（"
            "139.7mm 段厚壁 ID 108.1 / 薄壁 ID 111.16），仅服务 1D 到鞋时间口径，论文不可写成现场内径。",
            "井径/井斜剖面为现场提取包实测（108 点，20m 间隔，field_measured）：井径 204–265mm，"
            "井斜 0.19–8.19°（均值 2.81°，7420m 处最大 8.19°）；5400–5720m 技套内重叠段取首测点外推。",
            "居中度 standoff（0.38–0.48）为 **model_assumption**：现场仅扶正器布置（132 只，11–22m 间距）"
            "无实测居中度，且悬挂器坐挂失败、座底固井，实际居中度更低。",
            "casing_id_mm 实存 273.1mm 技术套管外径（命名语义历史遗留，全仓未统一；本任务不改名）。",
            "回接段（0–5399.78m，193.7mm 套管）数据明确排除：模型域 5400–7868m 不含回接段，防资料混入。",
        ),
    )

    fluids = (
        FluidSpec("钻井液", FluidRole.MUD, HU101_MUD_DENSITY_KG_M3, RheologyModel.BINGHAM, HU101_MUD_PV_PA_S, HU101_MUD_YP_PA),
        FluidSpec("平衡液", FluidRole.WASH, HU101_BALANCE_DENSITY_KG_M3, RheologyModel.BINGHAM, HU101_BALANCE_PV_PA_S, HU101_BALANCE_YP_PA),
        FluidSpec("驱油隔离液", FluidRole.SPACER, HU101_SPACER_DENSITY_KG_M3, RheologyModel.BINGHAM, HU101_SPACER_PV_PA_S, HU101_SPACER_YP_PA),
        FluidSpec("领浆", FluidRole.LEAD, HU101_LEAD_DENSITY_KG_M3, RheologyModel.POWER_LAW, power_law_n=HU101_LEAD_POWER_LAW_N, consistency_k=HU101_LEAD_CONSISTENCY_K),
        FluidSpec("尾浆", FluidRole.TAIL, HU101_TAIL_DENSITY_KG_M3, RheologyModel.POWER_LAW, power_law_n=HU101_TAIL_POWER_LAW_N, consistency_k=HU101_TAIL_CONSISTENCY_K),
        FluidSpec("轻泥浆", FluidRole.DISPLACEMENT, HU101_MUD_DENSITY_KG_M3, RheologyModel.BINGHAM, HU101_MUD_PV_PA_S, HU101_MUD_YP_PA),
        FluidSpec("中置液", FluidRole.DISPLACEMENT, HU101_MUD_DENSITY_KG_M3, RheologyModel.BINGHAM, HU101_MUD_PV_PA_S, HU101_MUD_YP_PA),
        FluidSpec("井浆", FluidRole.DISPLACEMENT, HU101_MUD_DENSITY_KG_M3, RheologyModel.BINGHAM, HU101_MUD_PV_PA_S, HU101_MUD_YP_PA),
    )

    schedule = PumpingSchedule(
        steps=(
            PumpingScheduleStep("注平衡液", "平衡液", 25.0, 1.2, remarks="呼101现场抽取：25m³@1.2，ρ1.85（field_measured）。"),
            PumpingScheduleStep("注驱油隔离液", "驱油隔离液", 25.0, 1.2, remarks="呼101现场抽取：25m³@1.2，ρ2.00（field_measured）。"),
            PumpingScheduleStep("注领浆", "领浆", HU101_LEAD_VOLUME_M3, 1.2, remarks="呼101现场抽取：47m³@1.2，ρ2.10，实际灰量 121t（field_measured）。"),
            PumpingScheduleStep("注尾浆", "尾浆", HU101_TAIL_VOLUME_M3, 1.2, remarks="呼101现场抽取：23m³@1.2，ρ1.90（field_measured）。"),
            PumpingScheduleStep("注后置液(管内)", "井浆", 2.0, 0.6, remarks="中置液/压塞液 2m³@0.6，仅作为管内压塞/占位流体（field_measured）。"),
            PumpingScheduleStep("注轻泥浆", "轻泥浆", 26.0, 1.5, remarks="按现场轻泥浆 26m³@1.5，ρ1.85（field_measured）。"),
            PumpingScheduleStep("注中置液", "中置液", 10.0, 1.5, remarks="中置液 10m³@1.5，ρ2.00（现场 pumping_schedule，field_measured）；LEGACY(2026-08-16 前): 1.2。"),
            PumpingScheduleStep("井浆快替", "井浆", 40.0, 1.0, remarks="设计表 40m³@1.0（现场替浆排量 1.5→0.55 变排量，模型按分段代理）。"),
            PumpingScheduleStep("井浆慢替", "井浆", 23.4, 0.55, remarks="补足现场总替量 101.4m³；理论碰压量 102.2m³，CSV 记替浆 52m³（口径差异见 notes）。"),
        ),
        notes=(
            "现场顺序：平衡液→驱油隔离液→领浆→尾浆→后置液(2m³)→轻泥浆→中置液→替浆（快/慢）。",
            "替浆合计 101.4m³（快替 40 + 慢替 23.4），与理论碰压量 102.2m³ 一致（差 0.8m³ 含压塞口径差异）；"
            "提取包 pumping_schedule 记'替浆钻井液 52m³ @1.5-0.55'与 109min 时长不符、字段疑似不完整，未直接采用，待人工复核。",
            "循环排混浆 85m³@2.4 与碰压/候凝为非顶替主输入步骤，不进入模型泵注序列。",
        ),
    )

    validation_data = ValidationData(
        cbl_summary_path=resolved_reference_root / "1003" / "100312.PDF",
        cbl_pass_rate=0.6277,
        job_report_path=resolved_reference_root / "提取数据" / "h101_data_extraction_report.md",
        notes=(
            "100312.PDF 正式 CBL/VDL 评价图：测量井段 5390–7810m，水泥胶结质量合格率 62.77%（不合格，未达 70% 红线），"
            "结构化 cbl_pass_rate=0.6277（interpreted，Vision+图头 OCR）。",
            "CBL 评价窗 5700–7810m 为'裸眼段顶 5700 + 测量段底 7810'组合口径（继承 legacy）；"
            "5390–5699.8m 双层套管不评价、7810–7868m 悬空段（尾管鞋 7868m 未测 58m）无 CBL 数据。",
            "5390–5955m 段为 Vision 逐段分析（high 置信）；5955–7810m 为自动化像素推断（low-medium 置信），论文定量验证只引用正式 62.77%。",
            "cbl_summary_path 指向参考文档/呼101/1003/100312.PDF；若原始文件路径变动仅需更新此处。",
        ),
    )
    return well_spec, fluids, schedule, validation_data


def export_hu101_sync_card_markdown(output_dir: Path) -> Path:
    """导出 呼101 同步画像卡为 Markdown 文件。"""

    well_spec, fluids, schedule, _ = load_hu101_tailpipe()
    solver = CasingFlowSolver()
    result = solver.run(well_spec, fluids, schedule)
    provenance = WELL_PROVENANCE.get("呼101")
    if provenance is None:
        from cemdisp.data.provenance import _WELL_NAME_ALIASES
        canonical = _WELL_NAME_ALIASES.get("呼101", "呼101")
        provenance = WELL_PROVENANCE[canonical]
    sync_card = build_sync_card("呼101", result.shoe_timeline, provenance)

    output_path = output_dir / ("呼101_同步画像卡.md")
    lines = [
        "# 呼101 同步画像卡",
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