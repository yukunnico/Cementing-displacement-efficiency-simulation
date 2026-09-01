"""
呼探1-001井（HT1-001）139.7+168.3mm双径尾管段标准数据加载器

本模块把呼探1-001井（HT1-001）现场数据包中的 139.7+168.3mm 双径尾管段
资料整理为 cemdisp 标准输入结构。呼探1-001井（HT1-001）不是呼探1井，
两口井虽然具有相同的双径尾管结构，但井眼、流体流变和施工程序参数不同。

当前求解目标为下部 139.7mm 尾管段，因此上部 168.3mm 重叠井段用等效井眼
直径保面积近似；鞋口滞后体积按 HT1-001 明确内径分段累加。

2026-08-16 核对修正（参考文档/现场资料提取/00_执行记录/其他井loader核对_2026-08-16.md §4）：
- 悬挂器 5469.711 -> 5460.159m（施工设计 simulation_top/liner_hanger_md）；
- 井径剖面由合成（上段等效 222.3mm）切换为现场提取包 69 点实测 caliper_profile.csv
  （上段 5670–7441 实测均值 247.7mm，环空体积较旧等效提高约 9%）；

2026-08-29 校准更新（0708 原件三方核对：重新提取校准_2026-08-29/ht1_001_校准核对.md §8）：
- 悬挂器改回实测 5469.711m（0708 三源：作业史"设计返高:5469.711m"+套管数据表三方核对版逐根累计闭合）；
  2026-08-16 的"5469.711→5460.159"修正方向反转登记——5460.159 为施工设计施工前预估（差 9.552m）；
  重叠段随之缩短为 200.289m（5469.711–5670m）；
- 鞋口滞后体积重估：83.2m³（52+26.0+5.2，无 0708 依据）→ 94.5m³（施工设计 6.4.1 分段内容积链
  45.8+17.4+26.0+4.3，理论总量含后置液；小结 5.3 全管内容积 95m³ 作上界；实际泵入 93.7m³ 含压塞/纯替浆 91.7m³ 各口径注记）；
- n/K 备注澄清：docx 复检报告与 PDF 委托检测报告（W301-2025045）为两份独立报告（非 OCR 口径差异），loader 取 docx 口径；
- 20512/20514 补规格：365.13mm（二开）浮箍浮鞋出厂检验报告，非尾管段附件、非 CBL。
- 井斜剖面由合成 2.0–7.8° 切换为 69 点实测 inclination_profile.csv（0.15–1.60°，
  最大 1.7665°@7665m，近似直井）——旧井斜量级错误，直接影响模型重力项；
- 下段尾管内径 108.04 -> 107.94mm（139.7−2×15.88）；
- 居中度按设计模拟 80.4% 取均一 standoff=0.804（无实测）；
- CBL 定量确证缺失（PDF 实为化验报告），评价窗分口径标注 cbl / formation_target / model_focus。
"""

from __future__ import annotations

import warnings

import csv
from collections.abc import Callable, Iterable
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
# LEGACY(2026-08-16 前): DEFAULT_REFERENCE_ROOT=参考文档/呼探1-001（旧数据包目录）
DEFAULT_REFERENCE_ROOT = PROJECT_ROOT / "参考文档" / "现场资料提取" / "ht1_001_呼探1-001"
# 现场提取包剖面 CSV（hu102/hu103 同款直读方式；均 69 点 30m 间隔电测/测斜）
DEFAULT_CALIPER_CSV = DEFAULT_REFERENCE_ROOT / "caliper_profile.csv"
DEFAULT_INCLINATION_CSV = DEFAULT_REFERENCE_ROOT / "inclination_profile.csv"

# 呼探1-001井段几何参数，来源：参考文档/现场资料提取/ht1_001_呼探1-001/（well_geometry.csv / casing_liner_string.csv）。
HT1_001_WELL_NAME = "呼探1-001井（HT1-001）"
HT1_001_DRILLED_DEPTH_MD_M = 7746.0  # 实际完钻井深/尾管鞋深度（施工小结 section 9）。
# 2026-08-29 复核改回实测 5469.711m：0708 三源（作业史"设计返高:5469.711m"+套管数据表三方核对版
# 逐根累计闭合，套长 2276.289m）为实际下入口径。
# LEGACY(2026-08-29 前): 5460.159（施工设计 section 4.2/6.1 liner_hanger_md，施工前预估，差 9.552m；
# 2026-08-16 曾按该口径修正，方向已反转）。
HT1_001_HANGER_MD_M = 5469.711  # 尾管悬挂器/回接筒顶（实测口径）。
HT1_001_TOP_MD_M = HT1_001_HANGER_MD_M  # 模型剖面从悬挂器开始，兼容双径尾管等效处理。
HT1_001_UPPER_SECTION_BOTTOM_MD_M = 7174.938  # 168.3mm 上段尾管底界/变径位置（现场变扣 7174.941m，差 0.003m）。
HT1_001_BOTTOM_MD_M = HT1_001_DRILLED_DEPTH_MD_M  # 下段 139.7mm 尾管鞋。
HT1_001_SHOE_MD_M = HT1_001_DRILLED_DEPTH_MD_M
# 2026-08-29 语义统一：casing_id_mm 按 PACKAGE_REFERENCE 文档语义存"外层套管内径"——
# 273.1mm 技套真实 ID=245.42mm（273.1−2×13.84，well_geometry casing_273_id）；OD 273.1 存档于 CASING_OD_MM。
# 字段不被求解器消费（环空几何由 hole_diameter_profile/重叠段等效井径表达），纯元数据口径。
# LEGACY(2026-08-29 前): 本字段存 OD 273.1（沿旧 hu102 口径，命名历史遗留）。
HT1_001_CASING_ID_MM = 245.42
HT1_001_CASING_OD_MM = 273.1  # 273.1mm 技术套管外径（OD 公称，存档）。
HT1_001_CASING_ID_TRUE_MM = 245.42  # 273.1mm 技套内径（与 CASING_ID_MM 同口径，供重叠段等效井径使用）。
HT1_001_UPPER_HOLE_NOMINAL_DIAMETER_MM = 241.3  # 上段井眼名义尺寸/钻头（5670–7441m）。
# LEGACY(2026-08-16 前): HT1_001_LOWER_HOLE_DIAMETER_MM=229.46 为 7441–7746 设计平均；
# 实测 caliper_profile.csv 该段均值 230.65mm，现由实测剖面直接表达，常量仅存档。
HT1_001_LOWER_HOLE_DIAMETER_MM = 229.46
HT1_001_BIT_DIAMETER_LOWER_MM = 215.9  # 下段钻头尺寸（7441–7746m 名义 215.9mm）。
HT1_001_UPPER_LINER_OD_MM = 168.3  # 上段 168.3mm 尾管外径（5469.711–7174.941m）。
HT1_001_LOWER_LINER_OD_MM = 139.7  # 下段 139.7mm 尾管外径，作为通用求解器参考外径。
HT1_001_LOWER_LINER_WALL_THICKNESS_MM = 15.88  # 139.7mm 管壁厚（施工设计/套管数据表）。
# LEGACY(2026-08-16 前): HT1_001_LOWER_LINER_ID_MM=108.04（数据包给定口径）；
# 核对后按 139.7−2×15.88=107.94mm（well_geometry liner_inner_diameter_139 / 套管数据表）。
HT1_001_LOWER_LINER_ID_MM = 107.94
# 上段 168.3mm 壁厚：三方核对套管表标注 14.27mm，但 168.3−2×14.27=139.76≠ID 138.76，
# 提取内部不一致；loader 保留 138.9mm（≈壁厚 14.7mm）口径，差异见核对报告 §4.1 与 loader notes。
HT1_001_UPPER_LINER_WALL_THICKNESS_MM = 14.7
HT1_001_UPPER_LINER_ID_MM = HT1_001_UPPER_LINER_OD_MM - 2.0 * HT1_001_UPPER_LINER_WALL_THICKNESS_MM  # 138.9
HT1_001_LOWER_CENTRALIZER_COUNT = 24  # 139.7mm 下段整体式扶正器数量（每 2 根 1 只）。
HT1_001_UPPER_CENTRALIZER_COUNT = 77  # 168.3mm 上段整体式扶正器数量。
HT1_001_CENTRALIZER_COUNT = HT1_001_LOWER_CENTRALIZER_COUNT + HT1_001_UPPER_CENTRALIZER_COUNT

# 重叠段 5469.711–5670m（200.289m）：168.3mm 尾管在 273.1mm 技套内（双层套管段），无裸眼井径测点，
# 等效井径取技套 ID 245.42mm（model_assumption，非实测）。
HT1_001_OVERLAP_EQUIVALENT_HOLE_DIAMETER_MM = HT1_001_CASING_ID_TRUE_MM
# 顶段井斜代理：取 5660m 首测点值 0.32°（重叠段无连续测斜，近似直井）。
HT1_001_TOP_INCLINATION_DEG_PROXY = 0.32
# 设计居中度 80.4%（施工设计 section 6.3 软件模拟，无实测；centralization_profile.csv 标 model_assumption）。
HT1_001_STANDOFF_DESIGN = 0.804


def _equivalent_hole_diameter_mm(actual_hole_mm: float, actual_od_mm: float, reference_od_mm: float) -> float:
    """在固定 reference_od 条件下，构造与原始环空面积等价的井眼直径。"""

    # 面积守恒：D_eq² - OD_ref² = D_actual² - OD_actual²。
    area_term = actual_hole_mm**2 - actual_od_mm**2 + reference_od_mm**2
    return math.sqrt(max(area_term, reference_od_mm**2))


# LEGACY(2026-08-16 前): 上段 168.3mm 段按保面积等效孔径计算，值约 222.3mm，
# 显著低于实测 caliper_profile.csv（5670–7441 实测均值 247.7mm，环空体积低估约 9%）。
# 2026-08-16 起井径剖面改读实测剖面，该等效值不再用于建模，仅保留常量供追溯。
HT1_001_UPPER_HOLE_DIAMETER_MM = _equivalent_hole_diameter_mm(
    actual_hole_mm=HT1_001_UPPER_HOLE_NOMINAL_DIAMETER_MM,
    actual_od_mm=HT1_001_UPPER_LINER_OD_MM,
    reference_od_mm=HT1_001_LOWER_LINER_OD_MM,
)


def _pipe_volume_m3(length_m: float, inner_diameter_mm: float) -> float:
    """按内径和长度估算管内容积，用于鞋口滞后体积计算。"""

    # 管内容积 = πr²L，内径单位从 mm 转换为 m。
    radius_m = inner_diameter_mm / 1000.0 / 2.0
    return math.pi * radius_m**2 * length_m


# 鞋口滞后体积锚点（2026-08-29 重估）：施工设计 6.4.1 分段内容积链 149.2 钻杆 45.8 + 127 钻杆 17.4
# + 168.3 段 26.0 + 139.7 段（阻位以上）4.3 ≈ 理论总量 94.5m³（含后置液）；
# 小结 5.3 全管内容积 95m³ 作上界；实际泵入替浆段 93.7m³（含压塞液 2）、纯替浆 91.7m³，各口径见 notes。
# LEGACY(2026-08-29 前): 83.2m³（52+26.0+5.2，52 系 52/π/7868 反推等效内径口径，无 0708 依据）。
HT1_001_SURFACE_TO_HANGER_EFFECTIVE_ID_MM = math.sqrt(4.0 * 52.0 / (math.pi * 7868.0)) * 1000.0  # LEGACY(2026-08-29 前) 派生代理，已不参与鞋口滞后计算
HT1_001_SHOE_LAG_VOLUME_M3 = 94.5  # 设计 6.4.1 理论总量（含后置液），见上注释与 notes。
HT1_001_LINER_ID_MM = HT1_001_LOWER_LINER_ID_MM

# ---- 数据标签（核对口径 2026-08-16；对齐 00_执行记录/现场提取数据_vs_cemdisp_loader核对汇总 §10.1）----
# 说明：field_measured=现场实测/设计实录；interpreted=由实测/设计推断；model_assumption=模型假设或等效/代理；
# loader_legacy=2026-08-16 前 loader 遗留值。标签仅作元数据登记，不影响数值。
HT1_001_DATA_LABELS: dict[str, str] = {
    "几何-完钻/鞋深": "field_measured（施工小结实际下深 7746m）",
    "几何-悬挂器": "field_measured（0708 实测 5469.711m：作业史'设计返高:5469.711m'+套管数据表三方核对版逐根累计闭合；"
                 "LEGACY(2026-08-29 前) 5460.159 系施工设计施工前预估，差 9.552m）",
    "几何-变径(168.3→139.7)": "field_measured（施工设计变扣 7174.941m；loader 取 7174.938m，差 0.003m）",
    "几何-casing_id_mm": "field_measured（2026-08-29 语义统一：存 273.1mm 技套真实内径 245.42mm；"
                 "OD 273.1 为 well_geometry casing_273_od；字段不被求解器消费，纯元数据）",
    "几何-上段井径剖面": "field_measured（69 点电测 caliper_profile.csv，5670–7441 均值 247.7mm）",
    "几何-下段井径剖面": "field_measured（69 点电测，7441–7746 均值 230.65mm）；设计平均 229.46mm 仅存档(loader_legacy)",
    "几何-重叠段(5469.711–5670,200.289m)": "model_assumption（168.3mm 尾管在 273.1mm 技套内，按技套 ID 245.42mm 等效，非裸眼实测）",
    "几何-下段尾管ID": "field_measured/interpreted（139.7−2×15.88=107.94mm）；LOADER_LEGACY 旧值 108.04（数据包口径）",
    "几何-上段尾管ID": "field_measured/interpreted（取 138.9mm；现场壁厚标注 14.27 与 ID 138.76 内部不一致）",
    "几何-井斜剖面": "field_measured（69 点电测，0.15–1.60°，最大 1.7665°@7665；近似直井）；LEGACY 合成 2.0–7.8° 废止",
    "几何-居中度 standoff": "model_assumption（设计模拟 80.4% 无实测，按均一 0.804）",
    "泵注程序": "field_measured（施工小结 12 步逐项精确匹配 pumping_schedule.csv，含替浆拆分口径）",
    "流变-钻井液": "field_measured（化验六速 PV51/YP6 @65C）",
    "流变-水泥浆 n/K": "field_measured（docx 复检报告口径，128→93℃降温流变：0.811/0.876、0.871/0.504、0.886/0.453）；"
                 "PDF《钻井流体分析实验中心》委托检测报告（W301-2025045，89℃幂律：0.825/0.842、0.816/0.512、0.807/0.527）"
                 "为另一份独立报告（2026-08-29 澄清：两份独立报告，非 OCR/来源口径差异）；loader 取 docx 口径",
    "流变-平衡液": "proxy（密度 field_measured；无实测流变，PV0.030/YP3.0 为代理）",
    "CBL 定量": "缺失（施工小结仅定性\"合格\"；油层尾管报告.pdf 为化验报告非 CBL；窗为 model_assumption）",
    "目标层段": "field_measured（录井解释显示层 target_intervals.csv）+ model_focus（7400–7500 模型关注窗）",
}

# 呼探1-001现场流体参数；HT1-001 与呼探1井不同，缺项仅按题设代理值补齐。
HT1_001_MUD_DENSITY_KG_M3 = 1920.0
HT1_001_BALANCE_DENSITY_KG_M3 = 1750.0
HT1_001_SPACER_DENSITY_KG_M3 = 1980.0
HT1_001_LEAD_DENSITY_KG_M3 = 2050.0
HT1_001_INTERMEDIATE_DENSITY_KG_M3 = 1900.0
HT1_001_TAIL_DENSITY_KG_M3 = 1900.0
HT1_001_DISPLACEMENT_DENSITY_KG_M3 = 1920.0
HT1_001_PLUG_DENSITY_KG_M3 = 1980.0
HT1_001_BUFFER_DENSITY_KG_M3 = 1980.0
HT1_001_BASE_FLUID_DENSITY_KG_M3 = 1020.0
HT1_001_WELL_MUD_DENSITY_KG_M3 = 1920.0
HT1_001_MUD_PV_PA_S = 0.051
HT1_001_MUD_YP_PA = 6.0
HT1_001_BALANCE_PV_PA_S = 0.030
HT1_001_BALANCE_YP_PA = 3.0
HT1_001_DISPLACEMENT_PV_PA_S = 0.051
HT1_001_DISPLACEMENT_YP_PA = 6.0
HT1_001_BASE_FLUID_PV_PA_S = 0.030
HT1_001_BASE_FLUID_YP_PA = 3.0
HT1_001_SPACER_POWER_LAW_N = 0.545
HT1_001_SPACER_CONSISTENCY_K = 1.338
HT1_001_LEAD_POWER_LAW_N = 0.811
HT1_001_LEAD_CONSISTENCY_K = 0.876
HT1_001_INTERMEDIATE_POWER_LAW_N = 0.871
HT1_001_INTERMEDIATE_CONSISTENCY_K = 0.504
HT1_001_TAIL_POWER_LAW_N = 0.886
HT1_001_TAIL_CONSISTENCY_K = 0.453

# 呼探1-001现场施工程序参数，按地面注入顺序排列。
HT1_001_BALANCE_VOLUME_M3 = 40.0
HT1_001_SPACER_VOLUME_M3 = 20.0
HT1_001_LEAD_VOLUME_M3 = 20.6
HT1_001_INTERMEDIATE_VOLUME_M3 = 28.7
HT1_001_TAIL_VOLUME_M3 = 22.1
HT1_001_PLUG_VOLUME_M3 = 2.0
HT1_001_FAST_MUD_VOLUME_M3 = 25.0
HT1_001_BUFFER_VOLUME_M3 = 10.0
HT1_001_BASE_FLUID_VOLUME_M3 = 3.0
HT1_001_WELL_MUD_FAST_VOLUME_M3 = 35.0
HT1_001_WELL_MUD_SLOW_VOLUME_M3 = 18.7
HT1_001_BALANCE_RATE_M3_MIN = 1.4
HT1_001_SPACER_RATE_M3_MIN = 1.0
HT1_001_CEMENT_RATE_M3_MIN = 1.0
HT1_001_PLUG_RATE_M3_MIN = 0.6
HT1_001_FAST_MUD_RATE_M3_MIN = 1.4
HT1_001_BUFFER_RATE_M3_MIN = 1.0
HT1_001_BASE_FLUID_RATE_M3_MIN = 1.0
HT1_001_WELL_MUD_FAST_RATE_M3_MIN = 1.0
HT1_001_WELL_MUD_SLOW_RATE_M3_MIN = 0.6


def _depth_points(values: tuple[tuple[float, float], ...]) -> tuple[DepthValuePoint, ...]:
    """把测深-数值元组转换为 WellSpec 使用的剖面点。"""

    # 统一把轻量元组转为冻结数据类，便于 WellSpec 校验和后续插值。
    return tuple(DepthValuePoint(depth_md_m=depth, value=value) for depth, value in values)


def _read_profile_rows(caliper_csv_path: Path, inclination_csv_path: Path) -> tuple[tuple[float, float, float], ...]:
    """读取现场提取包井径/井斜 CSV 并按测深合并为 (depth, hole_mm, inclination_deg) 元组。

    2026-08-16 起直接读取 caliper_profile.csv + inclination_profile.csv（均 69 点、30m 间隔）；
    两文件同一测深采样，逐点合并。重叠段 5460.159–5660m 无电测点，由调用方前插等效点。

    Raises:
        ValueError: 任一 CSV 为空，或井径深度在井斜 CSV 中缺对应测点。
    """

    caliper_by_depth: dict[float, float] = {}
    with caliper_csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            depth = row.get("md_m", row.get("depth_md_m"))
            hole = row.get("caliper_mm", row.get("hole_diameter_mm"))
            if depth is None or hole is None:
                continue
            caliper_by_depth[float(depth)] = float(hole)
    inclination_by_depth: dict[float, float] = {}
    with inclination_csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            depth = row.get("md_m", row.get("depth_md_m"))
            inclination = row.get("inclination_deg")
            if depth is None or inclination is None:
                continue
            inclination_by_depth[float(depth)] = float(inclination)
    if not caliper_by_depth:
        raise ValueError(f"呼探1-001 井径 CSV 为空: {caliper_csv_path}")
    if not inclination_by_depth:
        raise ValueError(f"呼探1-001 井斜 CSV 为空: {inclination_csv_path}")
    rows: list[tuple[float, float, float]] = []
    for depth in sorted(caliper_by_depth):
        inclination = inclination_by_depth.get(depth)
        if inclination is None:
            raise ValueError(f"呼探1-001 井径深度 {depth}m 在井斜 CSV 中无对应测点")
        rows.append((depth, caliper_by_depth[depth], inclination))
    return tuple(rows)


def _build_hole_profile(profile_rows: Iterable[tuple[float, float, float]]) -> tuple[DepthValuePoint, ...]:
    """从合并剖面行提取井径剖面。"""

    return tuple(DepthValuePoint(depth_md_m=depth, value=hole_diameter) for depth, hole_diameter, _ in profile_rows)


def _build_inclination_profile(profile_rows: Iterable[tuple[float, float, float]]) -> tuple[DepthValuePoint, ...]:
    """从合并剖面行提取井斜剖面（度）。"""

    return tuple(DepthValuePoint(depth_md_m=depth, value=inclination) for depth, _, inclination in profile_rows)


def _build_standoff_profile(profile_rows: Iterable[tuple[float, float, float]]) -> tuple[DepthValuePoint, ...]:
    """从合并剖面行构建居中度剖面。

    现场无实测居中度，按设计模拟 80.4%（centralization_profile.csv 标 model_assumption）
    取均一 standoff=0.804，不再沿用旧 0.55–0.75 合成剖面（低于设计 80.4%）。
    """

    return tuple(DepthValuePoint(depth_md_m=depth, value=HT1_001_STANDOFF_DESIGN) for depth, _, _ in profile_rows)


def _ht1_001_evaluation_windows() -> tuple[EvaluationWindow, ...]:
    """呼探1-001 验证窗口（分口径维护，核对汇总 §10.2 口径）。

    - cbl_eval_window：现场 CBL 定量完全缺失（施工小结仅定性\"合格\"），窗为模型假设（model_assumption），
      切勿写成现场窗；待独立 CBL/VDL 报告到位后可替换。
    - formation_target：现场录井解释显示/目的层（target_intervals.csv），为 field_measured。
    - model_focus_window：7400–7500m 为模型关注窗（历史首版暂定，现场无此窗；覆盖部分油气显示层）。
    """

    return (
        EvaluationWindow(
            name="CBL对比窗(模型假设,现场无CBL定量)",
            top_md_m=HT1_001_TOP_MD_M,
            bottom_md_m=HT1_001_BOTTOM_MD_M,
            window_type="model_focus",
        ),
        EvaluationWindow(name="K1q差气层1(录井显示)", top_md_m=7409.0, bottom_md_m=7423.0, window_type="formation_target"),
        EvaluationWindow(name="K1q差气层2(录井显示高峰)", top_md_m=7460.0, bottom_md_m=7462.0, window_type="formation_target"),
        EvaluationWindow(name="K1q油气同层(录井显示)", top_md_m=7504.0, bottom_md_m=7566.0, window_type="formation_target"),
        EvaluationWindow(name="J3k3气层(录井显示)", top_md_m=7642.0, bottom_md_m=7661.0, window_type="formation_target"),
        EvaluationWindow(name="K1s水层1(录井解释)", top_md_m=6060.0, bottom_md_m=6064.0, window_type="formation_target"),
        EvaluationWindow(name="K1s水层2(录井解释)", top_md_m=6092.0, bottom_md_m=6098.0, window_type="formation_target"),
        EvaluationWindow(name="K1q+J3k2主要目的层", top_md_m=6828.0, bottom_md_m=7746.0, window_type="formation_target"),
        EvaluationWindow(name="模型关注窗(替浆/掺混分析)", top_md_m=7400.0, bottom_md_m=7500.0, window_type="model_focus"),
    )


def _phase_fractions_for_role(role: FluidRole, *, split_cement_phases: bool) -> tuple[tuple[str, float], ...]:
    """把标准流体角色映射为环空二维模型相名称。"""

    # 分相口径下：领浆和中间浆合并为 lead，相当于前置水泥相；尾浆单独为 tail。
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

    # 未识别流体保守归入 mud 相，避免施工末端替浆流体误作水泥相。
    role = role_by_name.get(fluid_name, FluidRole.MUD)
    return _phase_fractions_for_role(role, split_cement_phases=split_cement_phases)


def load_ht1_001_tailpipe(
    *,
    reference_root: Path | None = None,
) -> tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载呼探1-001井（HT1-001）139.7+168.3mm双径尾管段标准模型输入。"""

    # 允许调用方覆盖参考资料根目录，默认指向项目内呼探1-001资料包。
    resolved_reference_root = reference_root or DEFAULT_REFERENCE_ROOT
    # 2026-08-16 起从现场提取包直读 69 点实测井径/井斜（hu102/hu103 同款方式）。
    measured_rows = _read_profile_rows(DEFAULT_CALIPER_CSV, DEFAULT_INCLINATION_CSV)
    # 重叠段 5460.159–5660m 无电测点：井径按 273.1mm 技套 ID 245.42mm 等效（model_assumption），
    # 井斜取 5660m 首测点 0.32° 作顶段代理（近似直井）。
    profile_rows = (
        (
            HT1_001_TOP_MD_M,
            HT1_001_OVERLAP_EQUIVALENT_HOLE_DIAMETER_MM,
            HT1_001_TOP_INCLINATION_DEG_PROXY,
        ),
    ) + measured_rows
    well_spec = WellSpec(
        well_name=HT1_001_WELL_NAME,
        top_md_m=HT1_001_TOP_MD_M,
        bottom_md_m=HT1_001_BOTTOM_MD_M,
        shoe_md_m=HT1_001_SHOE_MD_M,
        hanger_md_m=HT1_001_HANGER_MD_M,
        casing_id_mm=HT1_001_CASING_ID_MM,
        liner_od_mm=HT1_001_LOWER_LINER_OD_MM,
        liner_id_mm=HT1_001_LINER_ID_MM,
        hole_diameter_profile=_build_hole_profile(profile_rows),
        inclination_profile=_build_inclination_profile(profile_rows),
        standoff_profile=_build_standoff_profile(profile_rows),
        evaluation_windows=_ht1_001_evaluation_windows(),
        reference_root=resolved_reference_root,
        notes=(
            "呼探1-001井（HT1-001）不是呼探1井，本加载器使用 HT1-001 专属井眼、流体和施工程序参数。",
            "2026-08-29 复核：悬挂器改回实测 5469.711m（0708 三源；2026-08-16 曾按施工设计施工前预估改 5460.159m，方向反转登记）；"
            "井径/井斜为现场提取包 69 点实测剖面（caliper_profile.csv / inclination_profile.csv）。",
            "井径：上段 5670–7441 实测均值 247.7mm（名义 241.3）、下段 7441–7746 实测均值 230.65mm（名义 215.9）；旧上段等效 222.3mm 废除（环空体积低估约 9%）；"
            "小结 1.4.2 官方 1m 曲线统计值 247.83/229.46 与 CSV 30m 抽样 247.69/230.65 之差为统计口径，非数据错误（2026-08-29 补记）。",
            "井斜：实测近似直井 0.15–1.60°（最大 1.7665°@7665m），替换旧合成 2.0–7.8°（量级错误，原沿用呼探1井代理，直接影响重力项）。",
            "重叠段 5469.711–5670m（200.289m）为 168.3mm 尾管在 273.1mm 技套内的双层套管段，无裸眼电测点，按技套 ID 245.42mm 等效（model_assumption）。",
            "下段尾管内径 107.94mm（139.7−2×15.88）；上段尾管内径 138.9mm（现场壁厚标注 14.27mm 与 ID 138.76mm 内部不一致，保留口径见核对报告 §4.1）。",
            f"鞋口滞后体积 {HT1_001_SHOE_LAG_VOLUME_M3:.1f}m³：按施工设计 6.4.1 分段内容积链（45.8+17.4+26.0+4.3≈94.5，理论总量含后置液）；"
            f"小结 5.3 全管内容积 95m³ 作上界、实际泵入替浆段 93.7m³（含压塞液 2）/纯替浆 91.7m³ 各口径并存；"
            f"LEGACY(2026-08-29 前): 83.2m³（52+26.0+5.2）无 0708 依据。WellSpec.liner_id_mm 使用下段内径 {HT1_001_LINER_ID_MM:.2f}mm。",
            "居中度无实测，按设计模拟 80.4%（施工设计 section 6.3）取均一 standoff=0.804（model_assumption）。",
            "CBL 定量数据确认完全缺失：施工小结仅定性\"合格\"；油层尾管报告.pdf 为水泥浆化验报告（非 CBL）；"
            "20512/20514 为 365.13mm（二开）浮箍浮鞋出厂检验报告，非尾管段附件、非 CBL；20511 为浅层测井。CBL 评价窗为模型假设，勿写成现场窗。",
            "验证窗口分口径：cbl（模型假设）/ formation_target（录井显示层，field_measured）/ model_focus（7400–7500 模型关注窗）；油气显示 K1q 7409–7423/7460–7462/7504–7566、J3k3 气层 7642–7661、主要目的层 K1q+J3k2 6828–7746。",
            "数据标签口径见模块常量 HT1_001_DATA_LABELS（field_measured / interpreted / model_assumption / loader_legacy）。",
        ),
    )

    # HT1-001 流体清单：隔离液、领浆、中间浆、尾浆使用幂律模型，替浆类流体使用 Bingham 代理。
    fluids = (
        FluidSpec("钻井液", FluidRole.MUD, HT1_001_MUD_DENSITY_KG_M3, RheologyModel.BINGHAM, HT1_001_MUD_PV_PA_S, HT1_001_MUD_YP_PA),
        FluidSpec("替浆液", FluidRole.DISPLACEMENT, HT1_001_DISPLACEMENT_DENSITY_KG_M3, RheologyModel.BINGHAM, HT1_001_DISPLACEMENT_PV_PA_S, HT1_001_DISPLACEMENT_YP_PA),
        FluidSpec("平衡液", FluidRole.WASH, HT1_001_BALANCE_DENSITY_KG_M3, RheologyModel.BINGHAM, HT1_001_BALANCE_PV_PA_S, HT1_001_BALANCE_YP_PA),
        FluidSpec("隔离液", FluidRole.SPACER, HT1_001_SPACER_DENSITY_KG_M3, RheologyModel.POWER_LAW, power_law_n=HT1_001_SPACER_POWER_LAW_N, consistency_k=HT1_001_SPACER_CONSISTENCY_K),
        FluidSpec("领浆", FluidRole.LEAD, HT1_001_LEAD_DENSITY_KG_M3, RheologyModel.POWER_LAW, power_law_n=HT1_001_LEAD_POWER_LAW_N, consistency_k=HT1_001_LEAD_CONSISTENCY_K),
        FluidSpec("中间浆", FluidRole.INTERMEDIATE, HT1_001_INTERMEDIATE_DENSITY_KG_M3, RheologyModel.POWER_LAW, power_law_n=HT1_001_INTERMEDIATE_POWER_LAW_N, consistency_k=HT1_001_INTERMEDIATE_CONSISTENCY_K),
        FluidSpec("尾浆", FluidRole.TAIL, HT1_001_TAIL_DENSITY_KG_M3, RheologyModel.POWER_LAW, power_law_n=HT1_001_TAIL_POWER_LAW_N, consistency_k=HT1_001_TAIL_CONSISTENCY_K),
        FluidSpec("压塞液", FluidRole.OTHER, HT1_001_PLUG_DENSITY_KG_M3, RheologyModel.BINGHAM, HT1_001_DISPLACEMENT_PV_PA_S, HT1_001_DISPLACEMENT_YP_PA),
        FluidSpec("替钻井液", FluidRole.DISPLACEMENT, HT1_001_DISPLACEMENT_DENSITY_KG_M3, RheologyModel.BINGHAM, HT1_001_DISPLACEMENT_PV_PA_S, HT1_001_DISPLACEMENT_YP_PA),
        FluidSpec("中置液", FluidRole.DISPLACEMENT, HT1_001_BUFFER_DENSITY_KG_M3, RheologyModel.BINGHAM, HT1_001_DISPLACEMENT_PV_PA_S, HT1_001_DISPLACEMENT_YP_PA),
        FluidSpec("基液", FluidRole.DISPLACEMENT, HT1_001_BASE_FLUID_DENSITY_KG_M3, RheologyModel.BINGHAM, HT1_001_BASE_FLUID_PV_PA_S, HT1_001_BASE_FLUID_YP_PA),
        FluidSpec("井浆", FluidRole.DISPLACEMENT, HT1_001_WELL_MUD_DENSITY_KG_M3, RheologyModel.BINGHAM, HT1_001_DISPLACEMENT_PV_PA_S, HT1_001_DISPLACEMENT_YP_PA),
    )

    # HT1-001 地面施工程序：与呼探1井不同，含 3m³ 基液步骤，井浆替浆总量拆分为 35+18.7m³。
    schedule = PumpingSchedule(
        steps=(
            PumpingScheduleStep("注入平衡液", "平衡液", HT1_001_BALANCE_VOLUME_M3, HT1_001_BALANCE_RATE_M3_MIN, remarks="平衡液 40m³@1.4m³/min，角色 WASH。"),
            PumpingScheduleStep("注入隔离液", "隔离液", HT1_001_SPACER_VOLUME_M3, HT1_001_SPACER_RATE_M3_MIN, remarks="隔离液 20m³@1.0m³/min，角色 SPACER，使用 HT1-001 幂律流变。"),
            PumpingScheduleStep("注入领浆", "领浆", HT1_001_LEAD_VOLUME_M3, HT1_001_CEMENT_RATE_M3_MIN, remarks="领浆 20.6m³@1.0m³/min。"),
            PumpingScheduleStep("注入中间浆", "中间浆", HT1_001_INTERMEDIATE_VOLUME_M3, HT1_001_CEMENT_RATE_M3_MIN, remarks="中间浆 28.7m³@1.0m³/min，角色 INTERMEDIATE。"),
            PumpingScheduleStep("注入尾浆", "尾浆", HT1_001_TAIL_VOLUME_M3, HT1_001_CEMENT_RATE_M3_MIN, remarks="尾浆 22.1m³@1.0m³/min，使用现场实际体积。"),
            PumpingScheduleStep("注入压塞液（管内）", "压塞液", HT1_001_PLUG_VOLUME_M3, HT1_001_PLUG_RATE_M3_MIN, remarks="压塞液 2m³@0.6m³/min，仅作为管内占位，不作为水泥入环空体积。"),
            PumpingScheduleStep("替钻井液(快)", "替钻井液", HT1_001_FAST_MUD_VOLUME_M3, HT1_001_FAST_MUD_RATE_M3_MIN, remarks="替钻井液快替 25m³@1.4m³/min。"),
            PumpingScheduleStep("替保护液/中置液", "中置液", HT1_001_BUFFER_VOLUME_M3, HT1_001_BUFFER_RATE_M3_MIN, remarks="替保护液/中置液 10m³@1.0m³/min。"),
            PumpingScheduleStep("替基液", "基液", HT1_001_BASE_FLUID_VOLUME_M3, HT1_001_BASE_FLUID_RATE_M3_MIN, remarks="基液 3m³@1.0m³/min，密度 1020kg/m³。"),
            PumpingScheduleStep("井浆快替", "井浆", HT1_001_WELL_MUD_FAST_VOLUME_M3, HT1_001_WELL_MUD_FAST_RATE_M3_MIN, remarks="井浆快替 35m³@1.0m³/min（现场\"替钻井液2 53.7m³@1.2–0.6\"拆分的前段）。"),
            PumpingScheduleStep("井浆慢替", "井浆", HT1_001_WELL_MUD_SLOW_VOLUME_M3, HT1_001_WELL_MUD_SLOW_RATE_M3_MIN, remarks="井浆慢替 18.7m³@0.6m³/min（现场\"替钻井液2 53.7m³@1.2–0.6\"拆分的后段，合计 53.7m³）。"),
        ),
        notes=(
            "施工顺序按 HT1-001 现场数据：平衡液→隔离液→领浆→中间浆→尾浆→压塞液→四段替浆并含基液步骤。",
            "替浆总量 91.7m³ = 25 + 10 + 3 + 35 + 18.7m³，排量从 1.4 递减至 0.6m³/min。",
            "现场\"替钻井液2\"为单段 53.7m³@1.2–0.6（施工小结）；模型拆为井浆快替 35 + 井浆慢替 18.7 两段，总量一致（pumping_schedule.csv 口径）。",
            "泵注 12 步体积/密度/排量与现场提取包 pumping_schedule.csv 100% 吻合（field_measured）。",
            "压塞液保留在 PumpingSchedule 中用于管内时序占位；环空入口分相映射时按 mud 相处理。",
        ),
    )

    # 验证资料路径指向 HT1-001 现场提取包；本加载器不在运行时解析原始文件。
    validation_data = ValidationData(
        job_report_path=resolved_reference_root / "construction_events.csv",
        notes=(
            "呼探1-001加载器把现场提取包参数固化为模块常量（2026-08-16 起参考根目录指向 参考文档/现场资料提取/ht1_001_呼探1-001）。",
            "CBL 定量数据确认完全缺失：施工小结仅定性\"合格\"；油层尾管报告.pdf 为水泥浆化验报告（非 CBL）；"
            "20512/20514 为 365.13mm（二开）浮箍浮鞋出厂检验报告，非尾管段附件、非 CBL；20511 为浅层测井（1103–3940m 非尾管段）。",
            "cbl_pass_rate 保持 None（现场无数值）；评价窗为模型假设（model_assumption），待外部获取独立 CBL/VDL 报告后可升级。",
            "泵注（12 步）与流体密度/流变（钻井液 PV51/YP6、隔离液 n0.545/K1.338、领浆 n0.811/K0.876、中间浆 n0.871/K0.504、尾浆 n0.886/K0.453）已与现场提取包 100% 吻合（field_measured）。",
        ),
    )
    return well_spec, fluids, schedule, validation_data


def export_ht1_001_sync_card_markdown(output_dir: Path) -> Path:
    """导出 呼探1-001井（HT1-001） 同步画像卡为 Markdown 文件。"""

    well_spec, fluids, schedule, _ = load_ht1_001_tailpipe()
    solver = CasingFlowSolver()
    result = solver.run(well_spec, fluids, schedule)
    provenance = WELL_PROVENANCE.get("呼探1-001井（HT1-001）")
    if provenance is None:
        from cemdisp.data.provenance import _WELL_NAME_ALIASES
        canonical = _WELL_NAME_ALIASES.get("呼探1-001井（HT1-001）", "呼探1-001井（HT1-001）")
        provenance = WELL_PROVENANCE[canonical]
    sync_card = build_sync_card("呼探1-001井（HT1-001）", result.shoe_timeline, provenance)

    output_path = output_dir / ("呼探1-001井（HT1-001）_同步画像卡.md")
    lines = [
        "# 呼探1-001井（HT1-001） 同步画像卡",
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
