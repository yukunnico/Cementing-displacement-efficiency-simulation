"""
呼探1井尾管段标准数据加载器

本模块把呼探1井现场数据包中的尾管段资料整理为 cemdisp 标准输入结构。
呼探1井包含上部 168.3mm 与下部 139.7mm 双径尾管；当前求解目标为下部
139.7mm 尾管段，因此上部重叠井段用等效井眼直径保面积近似。

注意：呼探1井与 HT1-001井是两口不同的井，本模块专用于呼探1井。
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



PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_REFERENCE_ROOT = PROJECT_ROOT / "参考文档" / "呼探1"

# 呼探1井段参数，来源：参考文档/呼探1/提取数据/呼探1井_固井顶替模型数据包.json。
HU1_WELL_NAME = "呼探1井"
HU1_DRILLED_DEPTH_MD_M = 7746.0  # 实际完钻井深/尾管鞋深度。
HU1_HANGER_MD_M = 5469.711  # 尾管悬挂器位置。
HU1_TOP_MD_M = HU1_HANGER_MD_M  # 模型剖面从悬挂器开始，兼容双径尾管等效处理。
HU1_UPPER_SECTION_BOTTOM_MD_M = 7174.938  # 168.3mm 上段尾管底界/变径位置。
HU1_BOTTOM_MD_M = HU1_DRILLED_DEPTH_MD_M  # 下段 139.7mm 尾管鞋。
HU1_SHOE_MD_M = HU1_DRILLED_DEPTH_MD_M
HU1_CASING_ID_MM = 273.1  # 技术套管内径暂用值，参考呼101类似深井。
HU1_UPPER_HOLE_NOMINAL_DIAMETER_MM = 241.3  # 上段井眼名义尺寸。
HU1_LOWER_HOLE_DIAMETER_MM = 215.9  # 下段目标井段名义井眼尺寸。
HU1_UPPER_LINER_OD_MM = 168.3  # 上段尾管外径。
HU1_LOWER_LINER_OD_MM = 139.7  # 下段尾管外径，作为通用求解器参考外径。
HU1_LOWER_LINER_WALL_THICKNESS_MM = 15.8  # 139.7mm 管壁厚，参考呼102同口径。
HU1_LOWER_LINER_ID_MM = HU1_LOWER_LINER_OD_MM - 2.0 * HU1_LOWER_LINER_WALL_THICKNESS_MM
HU1_UPPER_LINER_ID_PROXY_MM = 150.42  # 168.3mm 上段尾管内径代理，沿用呼103上段估算口径。
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
HU1_MUD_DENSITY_KG_M3 = 1930.0
HU1_DISPLACEMENT_DENSITY_KG_M3 = 1920.0
HU1_BALANCE_DENSITY_KG_M3 = 1750.0
HU1_SPACER_DENSITY_KG_M3 = 1980.0
HU1_LEAD_DENSITY_KG_M3 = 2050.0
HU1_INTERMEDIATE_DENSITY_KG_M3 = 1900.0
HU1_TAIL_DENSITY_KG_M3 = 1900.0
HU1_MUD_PV_PA_S = 0.058  # 钻井液 PV 暂缺，参考呼101类似工况。
HU1_MUD_YP_PA = 5.0  # 钻井液 YP 暂缺，参考呼101类似工况。
HU1_DISPLACEMENT_PV_PA_S = 0.058
HU1_DISPLACEMENT_YP_PA = 5.0
HU1_BALANCE_PV_PA_S = 0.030
HU1_BALANCE_YP_PA = 3.0
HU1_SPACER_PV_PA_S = 0.030
HU1_SPACER_YP_PA = 5.0
HU1_LEAD_POWER_LAW_N = 0.825
HU1_LEAD_CONSISTENCY_K = 0.842
HU1_INTERMEDIATE_POWER_LAW_N = 0.816
HU1_INTERMEDIATE_CONSISTENCY_K = 0.512
HU1_TAIL_POWER_LAW_N = 0.807
HU1_TAIL_CONSISTENCY_K = 0.527

# 呼探1现场施工程序参数，按地面注入顺序排列。
HU1_BALANCE_VOLUME_M3 = 40.0
HU1_SPACER_VOLUME_M3 = 20.0
HU1_LEAD_VOLUME_M3 = 20.6
HU1_INTERMEDIATE_VOLUME_M3 = 28.7
HU1_TAIL_VOLUME_M3 = 22.1
HU1_PLUG_VOLUME_M3 = 2.0
HU1_LIGHT_MUD_VOLUME_M3 = 26.0
HU1_BUFFER_VOLUME_M3 = 10.0
HU1_FAST_DISPLACEMENT_VOLUME_M3 = 40.0
HU1_SLOW_DISPLACEMENT_VOLUME_M3 = 17.7
HU1_FRONT_RATE_M3_MIN = 1.4
HU1_CEMENT_RATE_M3_MIN = 1.0
HU1_PLUG_RATE_M3_MIN = 0.6
HU1_LIGHT_MUD_RATE_M3_MIN = 1.4
HU1_BUFFER_RATE_M3_MIN = 1.2
HU1_FAST_DISPLACEMENT_RATE_M3_MIN = 1.0
HU1_SLOW_DISPLACEMENT_RATE_M3_MIN = 0.6


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
            "呼探1缺实测井径/井斜 CSV，当前按名义井径、呼101类似井斜和 101 只扶正器代理构造剖面。",
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
        FluidSpec("压塞液", FluidRole.OTHER, HU1_DISPLACEMENT_DENSITY_KG_M3, RheologyModel.BINGHAM, HU1_DISPLACEMENT_PV_PA_S, HU1_DISPLACEMENT_YP_PA),
        FluidSpec("轻泥浆", FluidRole.DISPLACEMENT, HU1_DISPLACEMENT_DENSITY_KG_M3, RheologyModel.BINGHAM, HU1_DISPLACEMENT_PV_PA_S, HU1_DISPLACEMENT_YP_PA),
        FluidSpec("中置液", FluidRole.DISPLACEMENT, HU1_DISPLACEMENT_DENSITY_KG_M3, RheologyModel.BINGHAM, HU1_DISPLACEMENT_PV_PA_S, HU1_DISPLACEMENT_YP_PA),
        FluidSpec("井浆", FluidRole.DISPLACEMENT, HU1_DISPLACEMENT_DENSITY_KG_M3, RheologyModel.BINGHAM, HU1_DISPLACEMENT_PV_PA_S, HU1_DISPLACEMENT_YP_PA),
    )

    schedule = PumpingSchedule(
        steps=(
            PumpingScheduleStep("注入平衡液", "平衡液", HU1_BALANCE_VOLUME_M3, HU1_FRONT_RATE_M3_MIN, remarks="平衡液 40m³@1.4m³/min，角色 WASH。"),
            PumpingScheduleStep("注入隔离液", "隔离液", HU1_SPACER_VOLUME_M3, HU1_FRONT_RATE_M3_MIN, remarks="隔离液 20m³@1.4m³/min，角色 SPACER。"),
            PumpingScheduleStep("注入领浆", "领浆", HU1_LEAD_VOLUME_M3, HU1_CEMENT_RATE_M3_MIN, remarks="领浆 20.6m³@1.0m³/min。"),
            PumpingScheduleStep("注入中间浆", "中间浆", HU1_INTERMEDIATE_VOLUME_M3, HU1_CEMENT_RATE_M3_MIN, remarks="中间浆 28.7m³@1.0m³/min，角色 INTERMEDIATE。"),
            PumpingScheduleStep("注入尾浆", "尾浆", HU1_TAIL_VOLUME_M3, HU1_CEMENT_RATE_M3_MIN, remarks="尾浆 22.1m³@1.0m³/min。"),
            PumpingScheduleStep("注入压塞液（管内）", "压塞液", HU1_PLUG_VOLUME_M3, HU1_PLUG_RATE_M3_MIN, remarks="压塞液 2m³@0.6m³/min，仅作为管内占位，不作为水泥入环空体积。"),
            PumpingScheduleStep("替浆轻泥浆", "轻泥浆", HU1_LIGHT_MUD_VOLUME_M3, HU1_LIGHT_MUD_RATE_M3_MIN, remarks="替浆方案：轻泥浆 26m³@1.4m³/min。"),
            PumpingScheduleStep("替浆中置液", "中置液", HU1_BUFFER_VOLUME_M3, HU1_BUFFER_RATE_M3_MIN, remarks="替浆方案：中置液 10m³@1.2m³/min。"),
            PumpingScheduleStep("井浆快替", "井浆", HU1_FAST_DISPLACEMENT_VOLUME_M3, HU1_FAST_DISPLACEMENT_RATE_M3_MIN, remarks="替浆方案：井浆快替 40m³@1.0m³/min。"),
            PumpingScheduleStep("井浆慢替", "井浆", HU1_SLOW_DISPLACEMENT_VOLUME_M3, HU1_SLOW_DISPLACEMENT_RATE_M3_MIN, remarks="替浆方案：井浆慢替 17.7m³@0.6m³/min。"),
        ),
        notes=(
            "施工顺序按现场数据：平衡液→隔离液→领浆→中间浆→尾浆→压塞液→四段替浆。",
            "替浆总量 93.7m³ = 26 + 10 + 40 + 17.7m³，排量从 1.4 递减至 0.6m³/min。",
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


def export_hu1_sync_card_markdown(output_dir: Path) -> Path:
    """导出 呼探1井 同步画像卡为 Markdown 文件。"""

    well_spec, fluids, schedule, _ = load_hu1_tailpipe()
    solver = CasingFlowSolver()
    result = solver.run(well_spec, fluids, schedule)
    provenance = WELL_PROVENANCE.get("呼探1井")
    if provenance is None:
        from cemdisp.data.provenance import _WELL_NAME_ALIASES
        canonical = _WELL_NAME_ALIASES.get("呼探1井", "呼探1井")
        provenance = WELL_PROVENANCE[canonical]
    sync_card = build_sync_card("呼探1井", result.shoe_timeline, provenance)

    output_path = output_dir / ("呼探1井_同步画像卡.md")
    lines = [
        "# 呼探1井 同步画像卡",
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
