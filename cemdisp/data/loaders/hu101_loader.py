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
  【2026-09-02 退役】91.73 反推与 52m³ 滞后均已退役，见下方 2026-09-02 校准更新。
- standoff 0.38–0.48 剖面显式标 model_assumption：现场仅有扶正器布置（132 只）无实测居中度，
  且悬挂器坐挂失败、最终座底固井，居中度实际更低。
- cbl_pass_rate=0.6277 结构化进 ValidationData（参照 hu102=0.6665 / hu103=0.1206 写法），
  注明正式解释测量段 5390–7810m 口径；CBL 评价窗保持 legacy 的 5700–7810m（裸眼段顶+测量段底组合）。
- 中置液排量 1.2 → 1.5 m³/min（现场 pumping_schedule，field_measured）。
- casing_id_mm 命名语义问题（实存 OD）仅保留/强化注释，不改名（全仓统一需动 API，超出本任务）。

2026-08-29 校准更新（0708 原件三方核对：重新提取校准_2026-08-29/hu101_校准核对.md §8）：
- 技套等效外径 273.10 → 273.05mm（0708 实际施工记录 20133/20134.doc+100312 图头；273.1 为 10033 钻井设计公称）；
- 井径/井斜剖面 108 点 → 109 点（5700–7868m，2011116.xls Sheet2 电子版重建，2011113.doc 逐点互证）；
  井斜最大 8.19° 在 7440m（旧注 7420m 有误）；
- 泵注：隔离液/领浆/尾浆排量 1.2→1.0（现场简述时序值；设计 1.2）；中置液 1.5→1.0（现场名保护液，
  1.5 无现场依据）；替浆现场为井浆 63.4m³@1.0-0.55 连续降排量，40+23.4 为模型分段代理
  （旧注"CSV 记替浆 52m³ 待人工复核"已复核：52 为设计值，2011114 浆柱）；
- 领浆 n/K 0.719/0.815 → 0.844/0.381、尾浆 0.722/0.684 → 0.830/0.352（主检 2011122.pdf W301-22094，89℃；
  旧复检值 2011121.doc 93℃ 为另一份独立报告，LEGACY 保留；建议论文引用前做双口径敏感性重算）；
- 平衡液流变 proxy（0.030/3.0）→ 0708 设计六速口径 PV41/YP9.2（2011111 §1.4.2，60℃ 六速 100/59/43/26/6/5）；
- 鞋口滞后 52m³ 标 missing（0708 七份文本未找到出处，2026-08-29 查证）；liner_id 91.73 反推依赖它，
  91.73/52 为 legacy 等效口径——0708 分段真实内径链（149.2 钻杆 ID129.9 0–5397.21 + 168.3 ID138.9 +
  139.7 ID111.16/108.1）管内容积≈102.6m³，与理论碰压 102.2m³（2011114）呼应。

2026-09-02 管容链修复（0708 原件核实，tests/test_pipe_capacity_chain_fix.py 锁定）：
- shoe_lag_volume_m3=102.6m³：0708 分段真实内径链（149.2 送入钻杆 ID129.9 0–5397.21m +
  尾管两段）管内容积，与 2011114.txt 行157"理论需102.2方碰压"（实泵 101.8，行158）呼应；
  直接驱动 casing_flow._timeline_pipe_volume 与 _pipe_cross_section_area 双链。
- liner_id=91.73（52m³ 无出处滞后反推）退役：改用 139.7mm 段厚壁真实内径 108.10mm
  （=139.7−2×15.8，与 HU101_LINER_WALL_THICKNESS_MM 自洽）；该字段仅剩弥散半径/
  屈服半径两个次口径消费，管容主口径由 shoe_lag 驱动。

2026-09-15 R30 裁定（Task 10，tests/history/test_hu101_loader_standoff.py 锁定）：
- standoff 名义剖面回退 LEGACY 0.38–0.48（均值 0.429）：撤销 2026-09-11 临时改动
  0.80（反推情景，代码自认不得作验证数字，见 _ASSUMED_STANDOFF 处三代注释链）；
  上方 2026-08-16 块的"0.38–0.48"表述回退后重新自洽。
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
# 现场提取包井径/井斜 CSV：固井作业史 2011113.doc 电测数据（109 点，5700–7868m，
# 2011116.xls Sheet2 电子版重建，2011113.doc 逐点互证，field_measured）。
# LEGACY(2026-08-29 前): 108 点（20m 间隔）；LEGACY(2026-08-16 前): 合成剖面（上段等效 242.84 + 下段 215.9 常量；井斜 0.6→1.9°）。
DEFAULT_CALIPER_CSV = PROJECT_ROOT / "参考文档" / "现场资料提取" / "hu101_呼101" / "caliper_profile.csv"
DEFAULT_INCLINATION_CSV = PROJECT_ROOT / "参考文档" / "现场资料提取" / "hu101_呼101" / "inclination_profile.csv"

# 呼101尾管段井身结构参数（2026-08-16 按现场提取包 well_geometry/casing_liner_string/100312.PDF 核对）。
HU101_TOP_MD_M = 5400.0  # 模型段顶：水泥返高 5400m 口径（高于尾管串顶 5402.885m 约 2.9m，notes 注明）。
HU101_BOTTOM_MD_M = 7868.0  # 完钻/尾管鞋深度（浮鞋下深，field_measured）。
HU101_SHOE_MD_M = 7868.0
HU101_HANGER_MD_M = 5407.46  # 悬挂器本体下深（field_measured）；悬挂器跨度 5402.85–5407.46m。
# 2026-08-29 语义统一：casing_id_mm 按 PACKAGE_REFERENCE 文档语义存"外层套管内径"——
# 273.05mm 技套真实 ID=245.37mm（壁厚 13.84；0708 实际施工记录 20133/20134.doc+100312 图头；
# 273.05 为 OD 公称，见 HU101_TECH_CASING_OD_MM；10033 设计写 273.1）。技套鞋 5699.8m。
# 该字段不被求解器消费（环空几何由 hole_diameter_profile 表达），纯元数据口径。
# LEGACY(2026-08-29 前): 本字段存 OD 273.10（10033 设计公称口径，命名历史遗留）。
HU101_TECH_CASING_OD_MM = 273.05  # 273.05mm 技术套管外径（OD 公称；20133/20134.doc 实际施工记录）。
HU101_TECH_CASING_EQUIV_ID_MM = 245.37  # 273.05mm 技术套管真实内径（casing_id_mm 语义统一后口径）。
HU101_UPPER_SECTION_BOTTOM_MD_M = 6796.0  # 168.3mm 段底（变径变扣 6796.329m，差 0.329m 取整）；LEGACY 同。
HU101_UPPER_ACTUAL_HOLE_DIAMETER_MM = 260.35  # 上段实际井径（设计口径；实测 5700–6796 段均值 259.46）。
HU101_UPPER_ACTUAL_LINER_OD_MM = 168.30  # 上段尾管外径（5402.885–6796.136m，field_measured）。
HU101_LOWER_HOLE_DIAMETER_MM = 215.90  # 下段井眼名义/实测井径（139.7mm 段，7048–7868 段均值 216.24）。
HU101_LOWER_LINER_OD_MM = 139.70  # 下段尾管外径（6796.329–7868m，field_measured）。
HU101_LINER_WALL_THICKNESS_MM = 15.80  # 139.7mm 段厚壁段壁厚（ID 108.1）；上部薄壁段 14.27mm（ID 111.16）存在。
# 鞋口滞后体积（2026-09-02 现场核实修复）：102.6m³ = 0708 分段真实内径链管内容积
# （149.2 送入钻杆 ID129.9 0–5397.21m + 尾管两段），与 2011114.txt 行157"理论需102.2方碰压"
# （实泵 101.8，行158）呼应；直接驱动 casing_flow 管容双链（_timeline_pipe_volume /
# _pipe_cross_section_area）。
# LEGACY(2026-09-02 前): 52.0——0708 七份文本未找到出处（missing，2026-08-29 查证），已退役。
HU101_SHOE_LAG_VOLUME_M3 = 102.6


def _equivalent_hole_diameter_mm(actual_hole_mm: float, actual_od_mm: float, reference_od_mm: float) -> float:
    """在固定 reference_od 条件下，构造与原始环空面积等价的井眼直径。"""

    area_term = actual_hole_mm**2 - actual_od_mm**2 + reference_od_mm**2
    return math.sqrt(max(area_term, reference_od_mm**2))


HU101_UPPER_HOLE_DIAMETER_MM = _equivalent_hole_diameter_mm(
    actual_hole_mm=HU101_UPPER_ACTUAL_HOLE_DIAMETER_MM,
    actual_od_mm=HU101_UPPER_ACTUAL_LINER_OD_MM,
    reference_od_mm=HU101_LOWER_LINER_OD_MM,
)
# HU101_LINER_ID_MM：139.7mm 段厚壁真实内径 108.10mm（=139.70−2×15.80，与壁厚常量自洽）。
# 【2026-09-02 退役】旧公式 sqrt(4×shoe_lag/(π×shoe_md))×1000 由 52m³ 无出处滞后反推出
# 91.73mm 等效内径，随 52m³ 一并退役。1D 管容主口径现由 shoe_lag_volume_m3=102.6 直接驱动；
# 本字段仅剩弥散半径（casing_flow:532）与屈服半径（casing_flow:1020）两个次口径消费。
# **field_derived**：厚壁 ID 108.1（139.7mm 段）；上部薄壁段 14.27mm 壁厚对应 ID 111.16 存在。
HU101_LINER_ID_MM = HU101_LOWER_LINER_OD_MM - 2.0 * HU101_LINER_WALL_THICKNESS_MM

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
# 平衡液流变（2026-08-29 校准）：0708 设计六速 2011111 §1.4.2（60℃ 六速 100/59/43/26/6/5）→ PV41/YP9.2。
# LEGACY(2026-08-29 前): 0.030/3.0（proxy 代理）。
HU101_BALANCE_PV_PA_S = 0.041  # 平衡液 PV=41mPa·s（0708 设计六速口径，2011111 §1.4.2）。
HU101_BALANCE_YP_PA = 9.2  # 平衡液 YP=9.2Pa（同上）。
HU101_SPACER_PV_PA_S = 0.030  # 驱油隔离液塑粘 30mPa·s（field_measured）；YP=5 为代理（proxy）。
HU101_SPACER_YP_PA = 5.0
# 领/尾浆流变（2026-08-29 校准）：主检报告 2011122.pdf（W301-22094，89℃）实测幂律；
# 另一份独立复检报告 2011121.doc（93℃）记 0.719/0.815、0.722/0.684，为双报告并存，LEGACY 保留。
# 建议论文引用前做双口径敏感性重算。
HU101_LEAD_POWER_LAW_N = 0.844  # 领浆 n（主检 2011122.pdf W301-22094，89℃）；LEGACY(2026-08-29 前): 0.719（2011121.doc 93℃ 复检）。
HU101_LEAD_CONSISTENCY_K = 0.381  # 领浆 K Pa·s^n（同上）；LEGACY(2026-08-29 前): 0.815。
HU101_TAIL_POWER_LAW_N = 0.830  # 尾浆 n（主检 2011122.pdf W301-22094，89℃）；LEGACY(2026-08-29 前): 0.722。
HU101_TAIL_CONSISTENCY_K = 0.352  # 尾浆 K Pa·s^n（同上）；LEGACY(2026-08-29 前): 0.684。
# 复检口径常量（2011121.doc 固井公司化验中心，93℃）：供 rheology_source="recheck" 双口径敏感性使用。
HU101_LEAD_RECHECK_POWER_LAW_N = 0.719
HU101_LEAD_RECHECK_CONSISTENCY_K = 0.815
HU101_TAIL_RECHECK_POWER_LAW_N = 0.722
HU101_TAIL_RECHECK_CONSISTENCY_K = 0.684


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
# 居中度剖面三代记录（勿删历史，逐代留痕）：
# ① LEGACY（2026-09-11 前，现行口径）：名义剖面 0.38–0.48，均值 0.429，model_assumption
#    （现场仅扶正器布置 132 只，无实测居中度）。
# ② 2026-09-11 用户指令临时改动：全井常数 0.80（反推情景——与现场记录"悬挂器坐挂失败、
#    座底固井，居中度下降"（2011114.doc 七、固井质量分析 / cbl_evaluation.csv）方向相反，
#    代码自认**不得作为模型验证数字引用**；基线结果备份在
#    results/_呼101_基线备份_居中度0.43_2026-09-11/）。
# ③ 2026-09-15 R30 裁定回退 LEGACY：现场坐挂失败居中度下降 + ②口径代码自认不得作
#    验证数字 + 用户 2026-09-14 贴现场指令（坐挂失败井居中度更低）。端点结构沿用
#    HU101_TOP_MD_M→5400 起、HU101_BOTTOM_MD_M→7868 止。
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
    rheology_source: str = "primary",
) -> tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载呼101尾管段标准模型输入。

    measured_standoff: 若为 "between_centralizers" 或 "at_centralizers"，用从
        居中度检测图读取的实测剖面（扶正器间/扶正器处）替换 model_assumption 的
        0.38–0.48 名义剖面；None（默认）保持名义剖面。
    rheology_source: "primary"（默认）=主检报告 2011122.pdf（W301-22094，89℃）领/尾浆幂律；
        "recheck"=复检报告 2011121.doc（93℃）口径，供双化验口径敏感性重算。
    """

    if rheology_source not in ("primary", "recheck"):
        raise ValueError(
            f"rheology_source 必须为 'primary'/'recheck'，得到 {rheology_source!r}")
    lead_n = HU101_LEAD_POWER_LAW_N if rheology_source == "primary" else HU101_LEAD_RECHECK_POWER_LAW_N
    lead_k = HU101_LEAD_CONSISTENCY_K if rheology_source == "primary" else HU101_LEAD_RECHECK_CONSISTENCY_K
    tail_n = HU101_TAIL_POWER_LAW_N if rheology_source == "primary" else HU101_TAIL_RECHECK_POWER_LAW_N
    tail_k = HU101_TAIL_CONSISTENCY_K if rheology_source == "primary" else HU101_TAIL_RECHECK_CONSISTENCY_K
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
        # 全井管容链现场核实值（2026-09-02）：102.6m³ = 0708 分段真实内径链管内容积，
        # 驱动 casing_flow._timeline_pipe_volume 与 _pipe_cross_section_area 双链。
        shoe_lag_volume_m3=HU101_SHOE_LAG_VOLUME_M3,
        hole_diameter_profile=_depth_points(_build_hole_profile(caliper_rows)),
        inclination_profile=_depth_points(_build_inclination_profile(incl_rows)),
        standoff_profile=_depth_points(standoff_points),
        # e_clip 裁定（2026-09-06）：仅当使用实测居中度剖面时标记 standoff_measured=True，
        # 2D 求解器据此放开 e 截断上限（实测输入才允许放大模型响应）。
        # 名义剖面（measured_standoff=None）保持 model_assumption 属性。
        standoff_measured=(measured_standoff is not None),
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
            "管容链现场核实（2026-09-02）：shoe_lag_volume_m3=102.6m³ 为 0708 分段真实内径链管内容积"
            "（149.2 送入钻杆 ID129.9 0–5397.21m + 168.3 ID138.9 + 139.7 ID111.16/108.1 尾管两段），"
            "与 2011114.txt 行157'理论需102.2方碰压'、行158'累计泵冲到量碰压…=101.8方（96%上水效率）'呼应；"
            "liner_id=91.73（由 52m³ 无出处鞋口滞后反推，missing）已退役，改用 139.7mm 段厚壁真实内径 "
            "108.10mm；1D 管容主口径由 shoe_lag 驱动，liner_id 仅剩弥散/屈服半径次口径消费。",
            "井径/井斜剖面为现场提取包实测（109 点，5700–7868m，2011116.xls Sheet2 电子版重建，"
            "2011113.doc 逐点互证，field_measured）：井径 204–265mm，"
            "井斜 0.19–8.19°（均值 2.81°，7440m 处最大 8.19°）；5400–5720m 技套内重叠段取首测点外推。",
            "层位/井况补记（2026-08-29 校准，不改 EvaluationWindow）：实测水泥面 5402.85m（100312 图头）、"
            "人工井底 7809.4m、目的井段 6153–7741m。",
            "居中度 standoff 已被临时改为**全井常数 0.80**（2026-09-11 用户指令的反推情景；"
            "LEGACY: 0.38–0.48 model_assumption）。⚠️ 现场记录为'悬挂器坐挂失败、座底固井，居中度下降'"
            "（132 只扶正器布置，11–22m 间距，无实测居中度），本改动与该记录方向相反，"
            "**不得作为模型验证数字引用**。",
            "casing_id_mm 已按 2026-08-29 语义统一存 273.05mm 技术套管真实内径 245.37（OD 273.05 见 "
            "HU101_TECH_CASING_OD_MM；0708 实际施工记录 20133/20134.doc+100312 图头；字段不被求解器消费）。",
            "回接段（0–5399.78m，193.7mm 套管）数据明确排除：模型域 5400–7868m 不含回接段，防资料混入。",
        ),
    )

    fluids = (
        FluidSpec("钻井液", FluidRole.MUD, HU101_MUD_DENSITY_KG_M3, RheologyModel.BINGHAM, HU101_MUD_PV_PA_S, HU101_MUD_YP_PA),
        FluidSpec("平衡液", FluidRole.WASH, HU101_BALANCE_DENSITY_KG_M3, RheologyModel.BINGHAM, HU101_BALANCE_PV_PA_S, HU101_BALANCE_YP_PA),
        FluidSpec("驱油隔离液", FluidRole.SPACER, HU101_SPACER_DENSITY_KG_M3, RheologyModel.BINGHAM, HU101_SPACER_PV_PA_S, HU101_SPACER_YP_PA),
        FluidSpec("领浆", FluidRole.LEAD, HU101_LEAD_DENSITY_KG_M3, RheologyModel.POWER_LAW, power_law_n=lead_n, consistency_k=lead_k),
        FluidSpec("尾浆", FluidRole.TAIL, HU101_TAIL_DENSITY_KG_M3, RheologyModel.POWER_LAW, power_law_n=tail_n, consistency_k=tail_k),
        FluidSpec("轻泥浆", FluidRole.DISPLACEMENT, HU101_MUD_DENSITY_KG_M3, RheologyModel.BINGHAM, HU101_MUD_PV_PA_S, HU101_MUD_YP_PA),
        FluidSpec("中置液", FluidRole.DISPLACEMENT, HU101_MUD_DENSITY_KG_M3, RheologyModel.BINGHAM, HU101_MUD_PV_PA_S, HU101_MUD_YP_PA),
        FluidSpec("井浆", FluidRole.DISPLACEMENT, HU101_MUD_DENSITY_KG_M3, RheologyModel.BINGHAM, HU101_MUD_PV_PA_S, HU101_MUD_YP_PA),
    )

    schedule = PumpingSchedule(
        steps=(
            PumpingScheduleStep("注平衡液", "平衡液", 25.0, 1.2, remarks="呼101现场抽取：25m³@1.2，ρ1.85（field_measured）。"),
            PumpingScheduleStep("注驱油隔离液", "驱油隔离液", 25.0, 1.0, remarks="呼101现场抽取：25m³@1.0，ρ2.00（field_measured）；现场简述时序值，设计 1.2（2026-08-29 校准）。"),
            PumpingScheduleStep("注领浆", "领浆", HU101_LEAD_VOLUME_M3, 1.0, remarks="呼101现场抽取：47m³@1.0，ρ2.10，实际灰量 121t（field_measured）；现场简述时序值，设计 1.2（2026-08-29 校准）。"),
            PumpingScheduleStep("注尾浆", "尾浆", HU101_TAIL_VOLUME_M3, 1.0, remarks="呼101现场抽取：23m³@1.0，ρ1.90（field_measured）；现场简述时序值，设计 1.2（2026-08-29 校准）。"),
            PumpingScheduleStep("注后置液(管内)", "井浆", 2.0, 0.6, remarks="中置液/压塞液 2m³@0.6，仅作为管内压塞/占位流体（field_measured）。"),
            PumpingScheduleStep("注轻泥浆", "轻泥浆", 26.0, 1.5, remarks="按现场轻泥浆 26m³@1.5，ρ1.85（field_measured）。"),
            PumpingScheduleStep("注中置液", "中置液", 10.0, 1.0, remarks="中置液 10m³@1.0（现场名保护液，排量 1.0；1.5 无现场依据，2026-08-29 校准）；LEGACY(2026-08-29 前): 1.5；LEGACY(2026-08-16 前): 1.2。"),
            PumpingScheduleStep("井浆快替", "井浆", 40.0, 1.0, remarks="现场为井浆 63.4m³@1.0-0.55 连续降排量，40+23.4 为模型分段代理（本段 40m³@1.0）；旧注'CSV 记替浆 52m³ 待人工复核'已复核：52 为设计值（2011114 浆柱）。"),
            PumpingScheduleStep("井浆慢替", "井浆", 23.4, 0.55, remarks="现场为井浆 63.4m³@1.0-0.55 连续降排量，40+23.4 为模型分段代理（本段 23.4m³@0.55，补足总替量 101.4m³）。"),
        ),
        notes=(
            "现场顺序：平衡液→驱油隔离液→领浆→尾浆→后置液(2m³)→轻泥浆→中置液→替浆（快/慢）。",
            "替浆合计 101.4m³（快替 40 + 慢替 23.4），与理论碰压量 102.2m³ 一致（差 0.8m³ 含压塞口径差异）；"
            "现场实际为井浆 63.4m³@1.0-0.55 连续降排量，40/23.4 为模型分段代理；"
            "旧注'提取包 pumping_schedule 记替浆 52m³ 待人工复核'已复核：52 为设计值（2011114 浆柱，2026-08-29 校准）。",
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