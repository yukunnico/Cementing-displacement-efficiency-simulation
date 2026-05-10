"""
呼探1井尾管段标准数据加载器

本模块把呼探1井现场数据包中的尾管段资料整理为 cemdisp 标准输入结构。
呼探1井包含上部 168.3mm 与下部 139.7mm 双径尾管；当前求解目标为下部
139.7mm 尾管段，因此上部重叠井段用等效井眼直径保面积近似。

注意：呼探1井与 HT1-001井是两口不同的井，本模块专用于呼探1井。
"""

from __future__ import annotations

from collections.abc import Callable
import math
from pathlib import Path

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep
from cemdisp.data.validation_data import ValidationData
from cemdisp.data.well_spec import DepthValuePoint, EvaluationWindow, WellSpec
from cemdisp.models2d.boundary_bridge import AnnulusInletState


PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_REFERENCE_ROOT = PROJECT_ROOT / "参考文档" / "呼探1"

# 呼探1井段参数，数据来源：参考文档/呼探1/提取数据/well_spec.csv（证据等级A）。
HU1_WELL_NAME = "呼探1井"  # 呼探1井（HT1-001为另一口井）。
HU1_DRILLED_DEPTH_MD_M = 7746.0  # 实际完钻井深/尾管鞋深度。
HU1_HANGER_MD_M = 5469.711  # 尾管悬挂器位置。
HU1_TOP_MD_M = HU1_HANGER_MD_M  # 模型剖面从悬挂器开始，兼容双径尾管等效处理。
HU1_UPPER_SECTION_BOTTOM_MD_M = 7174.938  # 168.3mm 上段尾管底界/变径位置。
HU1_BOTTOM_MD_M = HU1_DRILLED_DEPTH_MD_M  # 下段 139.7mm 尾管鞋。
HU1_SHOE_MD_M = HU1_DRILLED_DEPTH_MD_M
# 呼探1井段参数，数据来源：参考文档/呼探1/提取数据/well_spec.csv（作业史 + 套管数据表，证据等级A）。
HU1_CASING_ID_MM = 273.1  # 技术套管内径暂用值，参考呼101类似深井（well_spec.csv未提供此值）。
HU1_UPPER_HOLE_NOMINAL_DIAMETER_MM = 241.3  # 上段井眼名义尺寸（well_spec.csv未提供）。
HU1_LOWER_HOLE_DIAMETER_MM = 215.9  # 下段五开钻头尺寸（作业史，证据等级A）。
HU1_UPPER_LINER_OD_MM = 168.3  # 168.3mm上段尾管外径（作业史）。
HU1_LOWER_LINER_OD_MM = 139.7  # 139.7mm下段尾管外径（作业史）。
HU1_LOWER_LINER_WALL_THICKNESS_MM = 15.88  # 139.7mm BG140V/BGT2套管壁厚（套管数据表，证据等级A）。
HU1_LOWER_LINER_ID_MM = HU1_LOWER_LINER_OD_MM - 2.0 * HU1_LOWER_LINER_WALL_THICKNESS_MM
HU1_UPPER_LINER_ID_PROXY_MM = 138.90  # 168.3mm上段尾管内径（套管数据表壁厚14.7mm，证据等级A）。
HU1_CENTRALIZER_COUNT = 101  # 现场整体式套管扶正器数量。


def _equivalent_hole_diameter_mm(actual_hole_mm: float, actual_od_mm: float, reference_od_mm: float) -> float:
    """在固定 reference_od 条件下，构造与原始环空面积等价的井眼直径。"""

    # 面积守恒：D_eq² - OD_ref² = D_actual² - OD_actual²。
    area_term = actual_hole_mm**2 - actual_od_mm**2 + reference_od_mm**2
    return math.sqrt(max(area_term, reference_od_mm**2))


HU1_UPPER_HOLE_DIAMETER_MM = _equivalent_hole_diameter_mm(
    actual_hole_mm=HU1_UPPER_HOLE_NOMINAL_DIAMETER_MM,
    actual_od_mm=HU1_UPPER_LINER_OD_MM,
    reference_od_mm=HU1_LOWER_LINER_OD_MM,
)


def _pipe_volume_m3(length_m: float, inner_diameter_mm: float) -> float:
    """按内径和长度估算管内容积，用于鞋口滞后体积计算。"""

    radius_m = inner_diameter_mm / 1000.0 / 2.0
    return math.pi * radius_m**2 * length_m


# 鞋口滞后体积估算：上部入井管柱按呼101深井等效内径代理，尾管段按双径尾管内径分段累加。
# 该体积只用于 field_order_realistic 边界的到鞋延迟，不代表真实单一管柱内径。
HU1_SURFACE_TO_HANGER_EFFECTIVE_ID_MM = math.sqrt(4.0 * 52.0 / (math.pi * 7868.0)) * 1000.0
HU1_SHOE_LAG_VOLUME_M3 = (
    _pipe_volume_m3(HU1_HANGER_MD_M, HU1_SURFACE_TO_HANGER_EFFECTIVE_ID_MM)
    + _pipe_volume_m3(HU1_UPPER_SECTION_BOTTOM_MD_M - HU1_HANGER_MD_M, HU1_UPPER_LINER_ID_PROXY_MM)
    + _pipe_volume_m3(HU1_BOTTOM_MD_M - HU1_UPPER_SECTION_BOTTOM_MD_M, HU1_LOWER_LINER_ID_MM)
)
HU1_LINER_ID_MM = math.sqrt(4.0 * HU1_SHOE_LAG_VOLUME_M3 / (math.pi * HU1_SHOE_MD_M)) * 1000.0

# 呼探1现场流体参数；密度由 g/cm³ 换算为 kg/m³，流变缺项按题设代理值补齐。
# 呼探1现场流体参数；密度由 g/cm³ 换算为 kg/m³。
# 数据来源：参考文档/呼探1/提取数据/fluid_spec.csv（化验报告 + 技术总结，证据等级A）。
# 驱油隔离液：化验报告密度1.98 g/cm³，n=0.545，K=1.338 Pa·s^n。
# 水泥浆（领/中/尾浆）：化验报告密度与流变参数。
# 压塞液/平衡液/替浆泥浆：同一油基钻井液体系，PV=51 mPa·s，YP=6 Pa（化验报告）。
HU1_MUD_DENSITY_KG_M3 = 1930.0  # 替浆泥浆密度1.92 g/cm³（技术总结）。
HU1_DISPLACEMENT_DENSITY_KG_M3 = 1920.0  # 替浆泥浆密度1.92 g/cm³（技术总结）。
HU1_BALANCE_DENSITY_KG_M3 = 1750.0  # 平衡液密度1.75 g/cm³（技术总结）。
HU1_SPACER_DENSITY_KG_M3 = 1980.0  # 驱油隔离液现场实际密度1.98 g/cm³（化验报告）。
HU1_LEAD_DENSITY_KG_M3 = 2050.0  # 领浆密度2.05 g/cm³（化验报告）。
HU1_INTERMEDIATE_DENSITY_KG_M3 = 1900.0  # 中间浆密度1.90 g/cm³（化验报告）。
HU1_TAIL_DENSITY_KG_M3 = 1900.0  # 尾浆密度1.90 g/cm³（化验报告）。
HU1_PLUG_DENSITY_KG_M3 = 1600.0  # 压塞液密度1.60 g/cm³（作业史）。
HU1_MUD_PV_PA_S = 0.051  # 替浆泥浆 PV=51 mPa·s（化验报告油基体系）。
HU1_MUD_YP_PA = 6.0  # 替浆泥浆 YP=6 Pa（化验报告油基体系）。
HU1_DISPLACEMENT_PV_PA_S = 0.051  # 替浆泥浆 PV=51 mPa·s（同上）。
HU1_DISPLACEMENT_YP_PA = 6.0  # 替浆泥浆 YP=6 Pa（同上）。
HU1_BALANCE_PV_PA_S = 0.051  # 平衡液 PV=51 mPa·s（化验报告油基体系）。
HU1_BALANCE_YP_PA = 6.0  # 平衡液 YP=6 Pa（化验报告油基体系）。
HU1_SPACER_PV_PA_S = 0.051  # 驱油隔离液 PV 代理（同油基体系实测PV=51 mPa·s）。
HU1_SPACER_YP_PA = 6.0  # 驱油隔离液 YP 代理（同油基体系实测YP=6 Pa）。
HU1_SPACER_POWER_LAW_N = 0.545  # 驱油隔离液幂律指数（化验报告）。
HU1_SPACER_CONSISTENCY_K = 1.338  # 驱油隔离液稠度系数 Pa·s^n（化验报告）。
HU1_LEAD_POWER_LAW_N = 0.825  # 领浆幂律指数（化验报告）。
HU1_LEAD_CONSISTENCY_K = 0.842  # 领浆稠度系数 Pa·s^n（化验报告）。
HU1_INTERMEDIATE_POWER_LAW_N = 0.816  # 中间浆幂律指数（化验报告）。
HU1_INTERMEDIATE_CONSISTENCY_K = 0.512  # 中间浆稠度系数 Pa·s^n（化验报告）。
HU1_TAIL_POWER_LAW_N = 0.807  # 尾浆幂律指数（化验报告）。
HU1_TAIL_CONSISTENCY_K = 0.527  # 尾浆稠度系数 Pa·s^n（化验报告）。

# 呼探1现场施工程序参数，按地面注入顺序排列。
# 数据来源：参考文档/呼探1/提取数据/pumping_schedule.csv（作业史，证据等级A）。
# 体积为实际施工值：平衡液40m³、隔离液20m³（作业史"注先导浆40m³（含平衡液）"），
# 领浆20.6m³、中间浆28.7m³、尾浆22.1m³、压塞液2m³（作业史），
# 首段替泥浆25m³（00:03-00:22替钻井液25m³）、
# 中置液10m³（00:22-00:33替保护液/中置液10m³）、
# 末段替泥浆53.7m³（00:35-01:34替钻井液53.7m³，合计替浆88.8m³）。
HU1_BALANCE_VOLUME_M3 = 40.0  # 平衡液实际40m³（作业史）。
HU1_SPACER_VOLUME_M3 = 20.0  # 驱油隔离液实际20m³（作业史）。
HU1_LEAD_VOLUME_M3 = 20.6  # 领浆实际20.6m³（作业史）。
HU1_INTERMEDIATE_VOLUME_M3 = 28.7  # 中间浆实际28.7m³（作业史）。
HU1_TAIL_VOLUME_M3 = 22.1  # 尾浆实际22.1m³（作业史）。
HU1_PLUG_VOLUME_M3 = 2.0  # 压塞液实际2m³（作业史，设计5m³）。
HU1_LIGHT_MUD_VOLUME_M3 = 25.0  # 首段替泥浆25m³@1.4m³/min（作业史）。
HU1_BUFFER_VOLUME_M3 = 10.0  # 中置液10m³@1.0m³/min（作业史替保护液/中置液）。
HU1_FAST_DISPLACEMENT_VOLUME_M3 = 25.0  # 快替段25m³（归入首段替泥浆，总量25m³）。
HU1_SLOW_DISPLACEMENT_VOLUME_M3 = 53.7  # 末段替泥浆53.7m³@1.2-0.6m³/min（作业史）。
HU1_FRONT_RATE_M3_MIN = 1.4  # 平衡液/隔离液/首段替浆排量1.4m³/min。
HU1_CEMENT_RATE_M3_MIN = 1.0  # 水泥浆排量1.0m³/min。
HU1_PLUG_RATE_M3_MIN = 0.6  # 压塞液排量0.6m³/min。
HU1_LIGHT_MUD_RATE_M3_MIN = 1.4  # 首段替泥浆排量1.4m³/min。
HU1_BUFFER_RATE_M3_MIN = 1.0  # 中置液排量1.0m³/min。
HU1_FAST_DISPLACEMENT_RATE_M3_MIN = 1.0  # 快替排量1.0m³/min。
HU1_SLOW_DISPLACEMENT_RATE_M3_MIN = 0.6  # 末段替泥浆排量0.6m³/min。


def _depth_points(values: tuple[tuple[float, float], ...]) -> tuple[DepthValuePoint, ...]:
    """把测深-数值元组转换为 WellSpec 使用的剖面点。"""

    return tuple(DepthValuePoint(depth_md_m=depth, value=value) for depth, value in values)


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

    role = role_by_name.get(fluid_name, FluidRole.MUD)
    return _phase_fractions_for_role(role, split_cement_phases=split_cement_phases)


def load_hu1_tailpipe(
    *,
    reference_root: Path | None = None,
) -> tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载呼探1井尾管段标准模型输入。"""

    resolved_reference_root = reference_root or DEFAULT_REFERENCE_ROOT
    well_spec = WellSpec(
        well_name=HU1_WELL_NAME,
        top_md_m=HU1_TOP_MD_M,
        bottom_md_m=HU1_BOTTOM_MD_M,
        shoe_md_m=HU1_SHOE_MD_M,
        hanger_md_m=HU1_HANGER_MD_M,
        casing_id_mm=HU1_CASING_ID_MM,
        liner_od_mm=HU1_LOWER_LINER_OD_MM,
        liner_id_mm=HU1_LINER_ID_MM,
        hole_diameter_profile=_depth_points(
            (
                # 呼探1缺实测井径 CSV，上段按等效井径保留双径环空面积。
                (HU1_TOP_MD_M, HU1_UPPER_HOLE_DIAMETER_MM),
                (6000.0, HU1_UPPER_HOLE_DIAMETER_MM),
                (6600.0, HU1_UPPER_HOLE_DIAMETER_MM),
                (HU1_UPPER_SECTION_BOTTOM_MD_M - 1.0, HU1_UPPER_HOLE_DIAMETER_MM),
                # 7174.938m 变径位置做平滑过渡，避免剖面突跳过硬。
                (
                    HU1_UPPER_SECTION_BOTTOM_MD_M,
                    0.5 * (HU1_UPPER_HOLE_DIAMETER_MM + HU1_LOWER_HOLE_DIAMETER_MM),
                ),
                (7400.0, HU1_LOWER_HOLE_DIAMETER_MM),
                (7500.0, HU1_LOWER_HOLE_DIAMETER_MM),
                (7710.0, HU1_LOWER_HOLE_DIAMETER_MM),
                (HU1_BOTTOM_MD_M, HU1_LOWER_HOLE_DIAMETER_MM),
            )
        ),
        inclination_profile=_depth_points(
            (
                # 井斜剖面无实测 CSV，参考呼101深井类似工况由 2° 缓增至约 7–8°。
                (HU1_TOP_MD_M, 2.0),
                (6000.0, 3.4),
                (6600.0, 5.2),
                (HU1_UPPER_SECTION_BOTTOM_MD_M - 1.0, 6.6),
                (HU1_UPPER_SECTION_BOTTOM_MD_M, 6.8),
                (7400.0, 7.6),
                (7500.0, 7.8),
                (7710.0, 7.3),
                (HU1_BOTTOM_MD_M, 7.0),
            )
        ),
        standoff_profile=_depth_points(
            (
                # 101只整体式扶正器代理居中度：总体控制在 0.55–0.75 范围。
                (HU1_TOP_MD_M, 0.64),
                (6000.0, 0.66),
                (6600.0, 0.62),
                (HU1_UPPER_SECTION_BOTTOM_MD_M - 1.0, 0.60),
                (HU1_UPPER_SECTION_BOTTOM_MD_M, 0.58),
                (7400.0, 0.70),
                (7500.0, 0.72),
                (7710.0, 0.62),
                (HU1_BOTTOM_MD_M, 0.58),
            )
        ),
        evaluation_windows=(
            EvaluationWindow(name="CBL评价井段", top_md_m=5700.0, bottom_md_m=7710.0, window_type="cbl"),
            EvaluationWindow(name="目标层段", top_md_m=7400.0, bottom_md_m=7500.0, window_type="target"),
        ),
        reference_root=resolved_reference_root,
        notes=(
            "呼探1缺实测井径/井斜 CSV，井眼直径按名义钻头尺寸 215.9mm（作业史），井斜按呼101类似深井代理。",
            "上段 168.3mm 尾管按等效井眼直径保面积处理，下段 139.7mm 尾管为模型目标井段。",
            f"鞋口滞后体积估算为 {HU1_SHOE_LAG_VOLUME_M3:.2f}m³，并反推等效 liner_id_mm={HU1_LINER_ID_MM:.2f}mm。",
            "CBL评价井段 5700–7710m 与目标层段 7400–7500m 为首版暂定，后续应随 CBL 数据确认调整。",
        ),
    )

    fluids = (
        FluidSpec("钻井液", FluidRole.MUD, HU1_MUD_DENSITY_KG_M3, RheologyModel.BINGHAM, HU1_MUD_PV_PA_S, HU1_MUD_YP_PA),
        FluidSpec("替浆液", FluidRole.DISPLACEMENT, HU1_DISPLACEMENT_DENSITY_KG_M3, RheologyModel.BINGHAM, HU1_DISPLACEMENT_PV_PA_S, HU1_DISPLACEMENT_YP_PA),
        FluidSpec("平衡液", FluidRole.WASH, HU1_BALANCE_DENSITY_KG_M3, RheologyModel.BINGHAM, HU1_BALANCE_PV_PA_S, HU1_BALANCE_YP_PA),
        FluidSpec("隔离液", FluidRole.SPACER, HU1_SPACER_DENSITY_KG_M3, RheologyModel.BINGHAM, HU1_SPACER_PV_PA_S, HU1_SPACER_YP_PA),
        FluidSpec("领浆", FluidRole.LEAD, HU1_LEAD_DENSITY_KG_M3, RheologyModel.POWER_LAW, power_law_n=HU1_LEAD_POWER_LAW_N, consistency_k=HU1_LEAD_CONSISTENCY_K),
        FluidSpec("中间浆", FluidRole.INTERMEDIATE, HU1_INTERMEDIATE_DENSITY_KG_M3, RheologyModel.POWER_LAW, power_law_n=HU1_INTERMEDIATE_POWER_LAW_N, consistency_k=HU1_INTERMEDIATE_CONSISTENCY_K),
        FluidSpec("尾浆", FluidRole.TAIL, HU1_TAIL_DENSITY_KG_M3, RheologyModel.POWER_LAW, power_law_n=HU1_TAIL_POWER_LAW_N, consistency_k=HU1_TAIL_CONSISTENCY_K),
        FluidSpec("压塞液", FluidRole.OTHER, HU1_PLUG_DENSITY_KG_M3, RheologyModel.POWER_LAW, power_law_n=0.689, consistency_k=0.557),
        FluidSpec("轻泥浆", FluidRole.DISPLACEMENT, HU1_DISPLACEMENT_DENSITY_KG_M3, RheologyModel.BINGHAM, HU1_DISPLACEMENT_PV_PA_S, HU1_DISPLACEMENT_YP_PA),
        FluidSpec("中置液", FluidRole.DISPLACEMENT, HU1_DISPLACEMENT_DENSITY_KG_M3, RheologyModel.BINGHAM, HU1_DISPLACEMENT_PV_PA_S, HU1_DISPLACEMENT_YP_PA),
        FluidSpec("井浆", FluidRole.DISPLACEMENT, HU1_DISPLACEMENT_DENSITY_KG_M3, RheologyModel.BINGHAM, HU1_DISPLACEMENT_PV_PA_S, HU1_DISPLACEMENT_YP_PA),
    )

    schedule = PumpingSchedule(
        steps=(
            PumpingScheduleStep("注入平衡液", "平衡液", HU1_BALANCE_VOLUME_M3, HU1_FRONT_RATE_M3_MIN, remarks="平衡液 40m³@1.4m³/min，密度1.75g/cm³（作业史实际值）。"),
            PumpingScheduleStep("注入隔离液", "隔离液", HU1_SPACER_VOLUME_M3, HU1_FRONT_RATE_M3_MIN, remarks="驱油隔离液 20m³@1.4m³/min，密度1.98g/cm³（作业史实际值）。"),
            PumpingScheduleStep("注入领浆", "领浆", HU1_LEAD_VOLUME_M3, HU1_CEMENT_RATE_M3_MIN, remarks="领浆 20.6m³@1.0m³/min，密度2.05g/cm³（作业史实际值）。"),
            PumpingScheduleStep("注入中间浆", "中间浆", HU1_INTERMEDIATE_VOLUME_M3, HU1_CEMENT_RATE_M3_MIN, remarks="中间浆 28.7m³@1.0m³/min，密度1.90g/cm³（作业史实际值）。"),
            PumpingScheduleStep("注入尾浆", "尾浆", HU1_TAIL_VOLUME_M3, HU1_CEMENT_RATE_M3_MIN, remarks="尾浆 22.1m³@1.0m³/min，密度1.90g/cm³（作业史实际值）。"),
            PumpingScheduleStep("注入压塞液（管内）", "压塞液", HU1_PLUG_VOLUME_M3, HU1_PLUG_RATE_M3_MIN, remarks="压塞液 2m³@0.6m³/min（设计5m³），密度1.60g/cm³，仅管内占位不进入环空。"),
            PumpingScheduleStep("替浆首段", "轻泥浆", HU1_LIGHT_MUD_VOLUME_M3, HU1_LIGHT_MUD_RATE_M3_MIN, remarks="首段替泥浆 25m³@1.4m³/min（作业史00:03-00:22替钻井液25m³）。"),
            PumpingScheduleStep("替浆中置液", "中置液", HU1_BUFFER_VOLUME_M3, HU1_BUFFER_RATE_M3_MIN, remarks="中置液 10m³@1.0m³/min（作业史00:22-00:33替保护液/中置液10m³）。"),
            PumpingScheduleStep("井浆快替", "井浆", HU1_FAST_DISPLACEMENT_VOLUME_M3, HU1_FAST_DISPLACEMENT_RATE_M3_MIN, remarks="快替 25m³@1.0m³/min（归入首段替浆，总量25m³）。"),
            PumpingScheduleStep("井浆慢替", "井浆", HU1_SLOW_DISPLACEMENT_VOLUME_M3, HU1_SLOW_DISPLACEMENT_RATE_M3_MIN, remarks="末段替泥浆 53.7m³@1.2→0.6m³/min（作业史00:35-01:34，合计替浆88.7m³）。"),
        ),
        notes=(
            "施工顺序按现场数据（作业史）：平衡液40m³→隔离液20m³→领浆20.6m³→中间浆28.7m³→尾浆22.1m³→压塞液2m³→四段替浆。",
            "替浆总量 88.7m³ = 首段替泥浆25 + 中置液10 + 末段替泥浆53.7m³（作业史合计替钻井液88.8m³）。",
            "压塞液保留在 PumpingSchedule 中用于管内时序占位；环空入口分相映射时按 mud 相处理。",
        ),
    )

    validation_data = ValidationData(
        job_report_path=resolved_reference_root / "提取数据" / "呼探1井_固井顶替模型数据包.json",
        notes=(
            "呼探1首版加载器不读取 JSON 文件，仅把数据包中的现场参数固化为模块常量。",
            "CBL 评价井段和目标层段为暂定窗口，后续需结合 CBL 原始数据更新 ValidationData 路径。",
        ),
    )
    return well_spec, fluids, schedule, validation_data


def build_hu1_annulus_inlet_provider(
    schedule: PumpingSchedule,
    fluids: tuple[FluidSpec, ...],
    annulus_boundary_mode: str = "field_order_realistic",
    *,
    split_cement_phases: bool = False,
) -> Callable[[float], AnnulusInletState]:
    """为呼探1尾管段构建环空入口边界提供器。

    支持两类口径：
    1. field_order_realistic：按估算鞋口滞后体积修正地面施工流体到达鞋口的顺序；
    2. sustained_tail / volume_limited / tail_then_mud：保留 legacy 对比模式。
    """

    role_by_name: dict[str, FluidRole] = {fluid.name: fluid.role for fluid in fluids}
    if not any(fluid.role == FluidRole.TAIL for fluid in fluids):
        raise ValueError("Hu1 边界提供器需要尾浆流体")

    # legacy 模式只跟踪前置液和水泥浆入环空，替浆期按固定口径处理。
    annulus_entry_steps = tuple(
        step
        for step in schedule.steps
        if role_by_name.get(step.fluid_name, FluidRole.MUD)
        in {FluidRole.WASH, FluidRole.SPACER, FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}
    )

    def _surface_state(time_s: float) -> tuple[float, float, str]:
        """返回地面累计注入体积、当前排量和当前施工阶段名。"""

        elapsed_s = 0.0
        injected_volume_m3 = 0.0
        for step in schedule.steps:
            duration_s = 0.0 if step.rate_m3_min <= 0.0 else step.volume_m3 / step.rate_m3_min * 60.0
            if time_s < elapsed_s + duration_s - 1e-12:
                active_volume = max(time_s - elapsed_s, 0.0) / 60.0 * step.rate_m3_min
                return injected_volume_m3 + active_volume, step.rate_m3_min / 60.0, step.step_name
            injected_volume_m3 += step.volume_m3
            elapsed_s += duration_s
        return injected_volume_m3, 0.0, "施工结束后保持"

    def _annulus_step_by_arrival_volume(arrival_volume_m3: float) -> PumpingScheduleStep | None:
        """按鞋口滞后后的到达体积定位 legacy 口径的前置液/水泥阶段。"""

        cumulative_m3 = 0.0
        for step in annulus_entry_steps:
            cumulative_m3 += step.volume_m3
            if arrival_volume_m3 < cumulative_m3 - 1e-12:
                return step
        return None

    def _surface_step_by_arrival_volume(arrival_volume_m3: float) -> PumpingScheduleStep | None:
        """按鞋口滞后后的到达体积定位所有地面泵注流体阶段。"""

        cumulative_m3 = 0.0
        for step in schedule.steps:
            cumulative_m3 += step.volume_m3
            if arrival_volume_m3 < cumulative_m3 - 1e-12:
                return step
        return None

    def _provider(time_s: float) -> AnnulusInletState:
        surface_volume_m3, flow_rate_m3_s, surface_stage_name = _surface_state(time_s)
        arrival_volume_m3 = max(surface_volume_m3 - HU1_SHOE_LAG_VOLUME_M3, 0.0)

        if annulus_boundary_mode == "field_order_realistic":
            if arrival_volume_m3 <= 0.0:
                return AnnulusInletState(
                    time_s,
                    flow_rate_m3_s,
                    f"{surface_stage_name}（鞋口前仍为钻井液）",
                    (("mud", 1.0),),
                )
            arrival_step = _surface_step_by_arrival_volume(arrival_volume_m3)
            if arrival_step is None:
                return AnnulusInletState(time_s, 0.0, "施工结束后保持（环空末端为替浆泥浆）", (("mud", 1.0),))
            return AnnulusInletState(
                time_s,
                flow_rate_m3_s,
                f"{arrival_step.step_name}（{HU1_SHOE_LAG_VOLUME_M3:.1f}m³鞋口滞后修正）",
                _phase_fractions_for_fluid(
                    arrival_step.fluid_name,
                    role_by_name,
                    split_cement_phases=split_cement_phases,
                ),
            )

        annulus_step = _annulus_step_by_arrival_volume(arrival_volume_m3)
        if arrival_volume_m3 <= 0.0 or annulus_step is None:
            if annulus_boundary_mode == "sustained_tail":
                return AnnulusInletState(
                    time_s,
                    flow_rate_m3_s,
                    f"{surface_stage_name}（环空保持尾浆）",
                    _phase_fractions_for_role(FluidRole.TAIL, split_cement_phases=split_cement_phases),
                )
            if annulus_boundary_mode == "volume_limited":
                return AnnulusInletState(
                    time_s,
                    0.0,
                    f"{surface_stage_name}（环空入口保持尾浆）",
                    _phase_fractions_for_role(FluidRole.TAIL, split_cement_phases=split_cement_phases),
                )
            if annulus_boundary_mode == "tail_then_mud":
                return AnnulusInletState(time_s, flow_rate_m3_s, f"{surface_stage_name}入环空", (("mud", 1.0),))
            raise ValueError(f"Unsupported annulus boundary mode: {annulus_boundary_mode}")

        return AnnulusInletState(
            time_s,
            flow_rate_m3_s,
            f"{annulus_step.step_name}（{HU1_SHOE_LAG_VOLUME_M3:.1f}m³鞋口滞后修正）",
            _phase_fractions_for_fluid(
                annulus_step.fluid_name,
                role_by_name,
                split_cement_phases=split_cement_phases,
            ),
        )

    return _provider
