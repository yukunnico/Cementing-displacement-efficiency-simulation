"""
呼探1-002井（HT1-002）139.7mm完井尾管段标准数据加载器

本模块把呼探1-002井现场数据包中的 139.7mm 完井尾管资料整理为 cemdisp
标准输入结构。与呼探1井不同，本井目标段为单一 139.7mm 尾管，不构造
双径尾管等效几何；上部地面至悬挂器的管内容积只用于鞋口滞后体积估算。
"""

from __future__ import annotations

import warnings

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

# 呼探1-002井段参数，来源：参考文档/呼探1-002/提取数据/呼探1-002井_固井顶替模型数据包.json。
HU2_WELL_NAME = "呼探1-002井（HT1-002）"
HU2_WELL_DEPTH_MD_M = 7559.0  # 实际完钻井深，鞋以下留有约 5m 口袋段。
HU2_HANGER_MD_M = 5292.5  # 139.7mm 完井尾管顶部/悬挂器位置。
HU2_TOP_MD_M = HU2_HANGER_MD_M  # 模型剖面从尾管悬挂器开始。
HU2_BOTTOM_MD_M = 7554.0  # 139.7mm 尾管下深/固井井段底部。
HU2_SHOE_MD_M = HU2_BOTTOM_MD_M  # 环空入口对应尾管鞋深度。
HU2_BIT_DIAMETER_MM = 190.5  # 五开钻头尺寸，仅作资料留痕。
HU2_AVERAGE_HOLE_DIAMETER_MM = 193.05  # 目标尾管段平均井径。
HU2_CASING_ID_MM = 195.0  # 上部 219.1mm 技术套管内径代理值。
HU2_LINER_OD_MM = 139.7  # 完井尾管外径，本井为单一口径。
HU2_LINER_WALL_THICKNESS_MM = 15.88  # 139.7mm BG140V/BG-T2 套管壁厚。
HU2_LINER_ID_ACTUAL_MM = HU2_LINER_OD_MM - 2.0 * HU2_LINER_WALL_THICKNESS_MM
HU2_CENTRALIZER_COUNT = 95  # φ139.7*190.5mm 整体式弹扶数量。
HU2_TARGET_CENTRALIZER_SPACING_M = 44.0  # 目的层扶正器间距。
HU2_NON_TARGET_CENTRALIZER_SPACING_M = 55.0  # 非目的层扶正器间距。
HU2_LARGE_HOLE_1_TOP_MD_M = 5631.0  # 第一段大肚子井眼起点。
HU2_LARGE_HOLE_1_BOTTOM_MD_M = 5940.0  # 第一段大肚子井眼终点。
HU2_LARGE_HOLE_1_AVERAGE_DIAMETER_MM = 206.65  # 第一段大肚子平均井径。
HU2_LARGE_HOLE_1_MAX_DIAMETER_MM = 268.48  # 第一段大肚子最大井径，作为资料备注。
HU2_LARGE_HOLE_2_TOP_MD_M = 7355.0  # 第二段大肚子井眼起点。
HU2_LARGE_HOLE_2_BOTTOM_MD_M = 7360.0  # 第二段大肚子井眼终点。
HU2_LARGE_HOLE_2_AVERAGE_DIAMETER_MM = 226.59  # 第二段大肚子平均井径。
HU2_LARGE_HOLE_2_MAX_DIAMETER_MM = 249.90  # 第二段大肚子最大井径，作为资料备注。
HU2_MAX_INCLINATION_DEG = 3.9  # 最大井斜，表明本井近垂直。
HU2_MAX_INCLINATION_MD_M = 7443.0  # 最大井斜所在测深。


def _pipe_volume_m3(length_m: float, inner_diameter_mm: float) -> float:
    """按内径和长度估算管内容积，用于鞋口滞后体积计算。"""

    # 由毫米内径换算为米半径后计算圆管容积。
    radius_m = inner_diameter_mm / 1000.0 / 2.0
    return math.pi * radius_m**2 * length_m


# 鞋口滞后体积估算：地面至悬挂器按 219.1mm 技术套管内径代理，尾管段按 108.04mm 实际内径累加。
# 该体积只用于 field_order_realistic 边界的到鞋延迟，不代表环空求解器几何直径。
HU2_SURFACE_TO_HANGER_EFFECTIVE_ID_MM = HU2_CASING_ID_MM
HU2_SHOE_LAG_VOLUME_M3 = _pipe_volume_m3(
    HU2_HANGER_MD_M,
    HU2_SURFACE_TO_HANGER_EFFECTIVE_ID_MM,
) + _pipe_volume_m3(
    HU2_SHOE_MD_M - HU2_HANGER_MD_M,
    HU2_LINER_ID_ACTUAL_MM,
)
HU2_LINER_ID_MM = HU2_LINER_ID_ACTUAL_MM

# 呼探1-002现场流体参数；Bingham 参数按实验值或代理值，水泥浆按实验幂律参数。
HU2_MUD_DENSITY_KG_M3 = 2060.0
HU2_DISPLACEMENT_DENSITY_KG_M3 = 2060.0
HU2_BALANCE_DENSITY_KG_M3 = 1850.0
HU2_SPACER_DENSITY_KG_M3 = 2100.0
HU2_LEAD_DENSITY_KG_M3 = 2100.0
HU2_INTERMEDIATE_DENSITY_KG_M3 = 1950.0
HU2_TAIL_DENSITY_KG_M3 = 1950.0
HU2_MUD_PV_PA_S = 0.058  # 钻井液 PV：58mPa·s。
HU2_MUD_YP_PA = 5.0  # 钻井液终切力 5Pa，作为 Bingham 屈服值。
HU2_DISPLACEMENT_PV_PA_S = 0.058  # 替浆液按同密度钻井液流变处理。
HU2_DISPLACEMENT_YP_PA = 5.0
HU2_BALANCE_PV_PA_S = 0.030  # 平衡液缺流变实测，使用代理 PV。
HU2_BALANCE_YP_PA = 3.0  # 平衡液缺流变实测，使用代理 YP。
HU2_SPACER_PV_PA_S = 0.030  # 隔离液缺流变实测，使用代理 PV。
HU2_SPACER_YP_PA = 5.0  # 隔离液缺流变实测，使用代理 YP。
HU2_LEAD_POWER_LAW_N = 0.811
HU2_LEAD_CONSISTENCY_K = 0.876
HU2_INTERMEDIATE_POWER_LAW_N = 0.871
HU2_INTERMEDIATE_CONSISTENCY_K = 0.504
HU2_TAIL_POWER_LAW_N = 0.886
HU2_TAIL_CONSISTENCY_K = 0.453
HU2_PLUG_DENSITY_KG_M3 = 1900.0  # 作业史记录压塞液实际密度 1.90g/cm³。
HU2_MIDDLE_FLUID_DENSITY_KG_M3 = 1900.0  # 作业史记录中置液/中间保护液实际密度 1.90g/cm³。

# 呼探1-002现场施工程序参数，按地面注入顺序排列。
HU2_BALANCE_VOLUME_M3 = 20.0
HU2_SPACER_VOLUME_M3 = 15.0
HU2_LEAD_VOLUME_M3 = 12.0
HU2_INTERMEDIATE_VOLUME_M3 = 14.0
HU2_TAIL_VOLUME_M3 = 12.0
HU2_PLUG_VOLUME_M3 = 5.0
HU2_FIRST_DISPLACEMENT_VOLUME_M3 = 12.0
HU2_MIDDLE_DISPLACEMENT_VOLUME_M3 = 15.0
HU2_FAST_DISPLACEMENT_VOLUME_M3 = 30.0
HU2_SLOW_DISPLACEMENT_VOLUME_M3 = 17.0
HU2_MAIN_RATE_M3_MIN = 0.8
HU2_PLUG_RATE_M3_MIN = 0.6
HU2_MIDDLE_DISPLACEMENT_RATE_M3_MIN = 0.6
HU2_FAST_DISPLACEMENT_RATE_M3_MIN = 0.8
HU2_SLOW_DISPLACEMENT_RATE_M3_MIN = 0.3


def _depth_points(values: tuple[tuple[float, float], ...]) -> tuple[DepthValuePoint, ...]:
    """把测深-数值元组转换为 WellSpec 使用的剖面点。"""

    # WellSpec 使用不可变剖面点，便于后续求解器按深度插值。
    return tuple(DepthValuePoint(depth_md_m=depth, value=value) for depth, value in values)


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
    well_spec = WellSpec(
        well_name=HU2_WELL_NAME,
        top_md_m=HU2_TOP_MD_M,
        bottom_md_m=HU2_BOTTOM_MD_M,
        shoe_md_m=HU2_SHOE_MD_M,
        hanger_md_m=HU2_HANGER_MD_M,
        casing_id_mm=HU2_CASING_ID_MM,
        liner_od_mm=HU2_LINER_OD_MM,
        liner_id_mm=HU2_LINER_ID_MM,
        hole_diameter_profile=_depth_points(
            (
                # 井径剖面由平均井径与两段大肚子井段构造，不引入双径尾管几何。
                (HU2_TOP_MD_M, HU2_AVERAGE_HOLE_DIAMETER_MM),
                (HU2_LARGE_HOLE_1_TOP_MD_M, HU2_LARGE_HOLE_1_AVERAGE_DIAMETER_MM),
                (HU2_LARGE_HOLE_1_BOTTOM_MD_M, HU2_AVERAGE_HOLE_DIAMETER_MM),
                (7000.0, HU2_AVERAGE_HOLE_DIAMETER_MM),
                (HU2_LARGE_HOLE_2_TOP_MD_M, HU2_LARGE_HOLE_2_AVERAGE_DIAMETER_MM),
                (HU2_LARGE_HOLE_2_BOTTOM_MD_M, HU2_AVERAGE_HOLE_DIAMETER_MM),
                (HU2_SHOE_MD_M, HU2_AVERAGE_HOLE_DIAMETER_MM),
            )
        ),
        inclination_profile=_depth_points(
            (
                # 井斜按近垂直井代理剖面构造，最大井斜 3.9° 位于 7443m。
                (HU2_TOP_MD_M, 1.0),
                (5800.0, 1.5),
                (6500.0, 2.5),
                (7000.0, 3.0),
                (HU2_MAX_INCLINATION_MD_M, HU2_MAX_INCLINATION_DEG),
                (HU2_SHOE_MD_M, 3.5),
            )
        ),
        standoff_profile=_depth_points(
            (
                # 95 只扶正器代理居中度：非目的层偏低，目的层近底部提高。
                (HU2_TOP_MD_M, 0.60),
                (5800.0, 0.62),
                (6500.0, 0.58),
                (7000.0, 0.65),
                (7350.0, 0.68),
                (7400.0, 0.72),
                (7500.0, 0.75),
                (HU2_SHOE_MD_M, 0.70),
            )
        ),
        evaluation_windows=(
            # CBL评价井段覆盖完整 139.7mm 完井尾管固井井段。
            EvaluationWindow(name="CBL评价井段", top_md_m=HU2_TOP_MD_M, bottom_md_m=HU2_BOTTOM_MD_M, window_type="cbl"),
            # 目标层段按资料给出的 7400–7500m 单独评价。
            EvaluationWindow(name="目标层段", top_md_m=7400.0, bottom_md_m=7500.0, window_type="target"),
        ),
        reference_root=resolved_reference_root,
        notes=(
            "呼探1-002目标段为单一 139.7mm 完井尾管，不采用呼探1双径尾管等效几何。",
            "井径剖面由平均井径 193.05mm 与 5631–5940m、7355–7360m 两段大肚子井段构造。",
            f"鞋口滞后体积估算为 {HU2_SHOE_LAG_VOLUME_M3:.2f}m³：地面至悬挂器按 195mm 代理内径，尾管段按 {HU2_LINER_ID_ACTUAL_MM:.2f}mm 内径。",
            "CBL评价井段为 5292.5–7554m 完整尾管段；质量标签仅为‘合格’弱监督代理，不等同顶替效率真值。",
        ),
    )

    # 七类环空模型核心流体，加上压塞液/中置液/井浆名称以覆盖现场 PumpingSchedule。
    fluids = (
        FluidSpec("钻井液", FluidRole.MUD, HU2_MUD_DENSITY_KG_M3, RheologyModel.BINGHAM, HU2_MUD_PV_PA_S, HU2_MUD_YP_PA),
        FluidSpec("替浆液", FluidRole.DISPLACEMENT, HU2_DISPLACEMENT_DENSITY_KG_M3, RheologyModel.BINGHAM, HU2_DISPLACEMENT_PV_PA_S, HU2_DISPLACEMENT_YP_PA),
        FluidSpec("平衡液", FluidRole.WASH, HU2_BALANCE_DENSITY_KG_M3, RheologyModel.BINGHAM, HU2_BALANCE_PV_PA_S, HU2_BALANCE_YP_PA),
        FluidSpec("隔离液", FluidRole.SPACER, HU2_SPACER_DENSITY_KG_M3, RheologyModel.BINGHAM, HU2_SPACER_PV_PA_S, HU2_SPACER_YP_PA),
        FluidSpec("领浆", FluidRole.LEAD, HU2_LEAD_DENSITY_KG_M3, RheologyModel.POWER_LAW, power_law_n=HU2_LEAD_POWER_LAW_N, consistency_k=HU2_LEAD_CONSISTENCY_K),
        FluidSpec("中间浆", FluidRole.INTERMEDIATE, HU2_INTERMEDIATE_DENSITY_KG_M3, RheologyModel.POWER_LAW, power_law_n=HU2_INTERMEDIATE_POWER_LAW_N, consistency_k=HU2_INTERMEDIATE_CONSISTENCY_K),
        FluidSpec("尾浆", FluidRole.TAIL, HU2_TAIL_DENSITY_KG_M3, RheologyModel.POWER_LAW, power_law_n=HU2_TAIL_POWER_LAW_N, consistency_k=HU2_TAIL_CONSISTENCY_K),
        FluidSpec("压塞液", FluidRole.OTHER, HU2_PLUG_DENSITY_KG_M3, RheologyModel.BINGHAM, HU2_DISPLACEMENT_PV_PA_S, HU2_DISPLACEMENT_YP_PA),
        FluidSpec("中置液", FluidRole.DISPLACEMENT, HU2_MIDDLE_FLUID_DENSITY_KG_M3, RheologyModel.BINGHAM, HU2_DISPLACEMENT_PV_PA_S, HU2_DISPLACEMENT_YP_PA),
        FluidSpec("井浆", FluidRole.DISPLACEMENT, HU2_DISPLACEMENT_DENSITY_KG_M3, RheologyModel.BINGHAM, HU2_DISPLACEMENT_PV_PA_S, HU2_DISPLACEMENT_YP_PA),
    )

    # 地面施工程序按现场注入顺序排列，47m³ 重泥浆替浆拆分为快替 30m³ 与慢替 17m³。
    schedule = PumpingSchedule(
        steps=(
            PumpingScheduleStep("注入平衡液", "平衡液", HU2_BALANCE_VOLUME_M3, HU2_MAIN_RATE_M3_MIN, remarks="平衡液 20m³@0.8m³/min，角色 WASH。"),
            PumpingScheduleStep("注入隔离液", "隔离液", HU2_SPACER_VOLUME_M3, HU2_MAIN_RATE_M3_MIN, remarks="隔离液实际密度 2.10g/cm³，15m³@0.8m³/min，角色 SPACER。"),
            PumpingScheduleStep("注入领浆", "领浆", HU2_LEAD_VOLUME_M3, HU2_MAIN_RATE_M3_MIN, remarks="领浆 12m³@0.8m³/min。"),
            PumpingScheduleStep("注入中间浆", "中间浆", HU2_INTERMEDIATE_VOLUME_M3, HU2_MAIN_RATE_M3_MIN, remarks="中间浆 14m³@0.8m³/min，角色 INTERMEDIATE。"),
            PumpingScheduleStep("注入尾浆", "尾浆", HU2_TAIL_VOLUME_M3, HU2_MAIN_RATE_M3_MIN, remarks="尾浆 12m³@0.8m³/min。"),
            PumpingScheduleStep("注入压塞液（管内）", "压塞液", HU2_PLUG_VOLUME_M3, HU2_PLUG_RATE_M3_MIN, remarks="压塞液 5m³@0.6m³/min，仅作为管内占位。"),
            PumpingScheduleStep("替浆泥浆(快)", "井浆", HU2_FIRST_DISPLACEMENT_VOLUME_M3, HU2_MAIN_RATE_M3_MIN, remarks="替浆泥浆 12m³@0.8m³/min，角色 DISPLACEMENT。"),
            PumpingScheduleStep("替浆中置液", "中置液", HU2_MIDDLE_DISPLACEMENT_VOLUME_M3, HU2_MIDDLE_DISPLACEMENT_RATE_M3_MIN, remarks="中置液 15m³，现场 0.8→0.5m³/min，模型用 0.6m³/min 代理。"),
            PumpingScheduleStep("井浆快替", "井浆", HU2_FAST_DISPLACEMENT_VOLUME_M3, HU2_FAST_DISPLACEMENT_RATE_M3_MIN, remarks="重泥浆 47m³ 中快替 30m³@0.8m³/min。"),
            PumpingScheduleStep("井浆慢替", "井浆", HU2_SLOW_DISPLACEMENT_VOLUME_M3, HU2_SLOW_DISPLACEMENT_RATE_M3_MIN, remarks="重泥浆 47m³ 中慢替 17m³@0.3m³/min。"),
        ),
        notes=(
            "施工顺序按现场记录：平衡液→隔离液→领浆→中间浆→尾浆→压塞液→三段替浆泥浆/中置液/重泥浆。",
            "替浆总量 79m³ = 5 + 12 + 15 + 47m³，其中压塞液用于管内时序占位。",
            "隔离液使用现场实际密度 2.10g/cm³，不与设计值 2.05g/cm³ 混用。",
        ),
    )

    # 校验资料路径指向数据包 JSON；本加载器不在运行时解析 JSON，只固化标准输入常量。
    validation_data = ValidationData(
        job_report_path=resolved_reference_root / "提取数据" / "呼探1-002井_固井顶替模型数据包.json",
        notes=(
            "呼探1-002首版加载器不读取 JSON 文件，仅把数据包中的现场参数固化为模块常量。",
            "质量标签为‘合格’定性结论，未提供数值型 CBL 合格率，不能作为顶替效率真值。",
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
