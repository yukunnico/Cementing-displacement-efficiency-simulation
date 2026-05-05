"""
呼探1-001井（HT1-001）139.7+168.3mm双径尾管段标准数据加载器

本模块把呼探1-001井（HT1-001）现场数据包中的 139.7+168.3mm 双径尾管段
资料整理为 cemdisp 标准输入结构。呼探1-001井（HT1-001）不是呼探1井，
两口井虽然具有相同的双径尾管结构，但井眼、流体流变和施工程序参数不同。

当前求解目标为下部 139.7mm 尾管段，因此上部 168.3mm 重叠井段用等效井眼
直径保面积近似；鞋口滞后体积按 HT1-001 明确内径分段累加。
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
DEFAULT_REFERENCE_ROOT = PROJECT_ROOT / "参考文档" / "呼探1-001"

# 呼探1-001井段几何参数，来源：参考文档/呼探1-001/提取数据/呼探1-001井_固井顶替模型数据包.json。
HT1_001_WELL_NAME = "呼探1-001井（HT1-001）"
HT1_001_DRILLED_DEPTH_MD_M = 7746.0  # 实际完钻井深/尾管鞋深度。
HT1_001_HANGER_MD_M = 5469.711  # 尾管悬挂器位置。
HT1_001_TOP_MD_M = HT1_001_HANGER_MD_M  # 模型剖面从悬挂器开始，兼容双径尾管等效处理。
HT1_001_UPPER_SECTION_BOTTOM_MD_M = 7174.938  # 168.3mm 上段尾管底界/变径位置。
HT1_001_BOTTOM_MD_M = HT1_001_DRILLED_DEPTH_MD_M  # 下段 139.7mm 尾管鞋。
HT1_001_SHOE_MD_M = HT1_001_DRILLED_DEPTH_MD_M
HT1_001_CASING_ID_MM = 273.1  # 技术套管内径暂用与呼探1相同代理值。
HT1_001_UPPER_HOLE_NOMINAL_DIAMETER_MM = 241.3  # 上段井眼名义尺寸。
HT1_001_LOWER_HOLE_DIAMETER_MM = 229.46  # HT1-001 下段平均井径，不使用呼探1井 215.9mm。
HT1_001_BIT_DIAMETER_LOWER_MM = 215.9  # 下段钻头尺寸，仅作为现场设计口径记录。
HT1_001_UPPER_LINER_OD_MM = 168.3  # 上段尾管外径。
HT1_001_LOWER_LINER_OD_MM = 139.7  # 下段尾管外径，作为通用求解器参考外径。
HT1_001_LOWER_LINER_WALL_THICKNESS_MM = 15.88  # HT1-001 明确 139.7mm 管壁厚。
# 按 HT1-001 数据包给定口径使用 108.04mm，避免后续鞋口滞后体积与资料口径不一致。
HT1_001_LOWER_LINER_ID_MM = 108.04
HT1_001_UPPER_LINER_WALL_THICKNESS_MM = 14.7  # HT1-001 明确 168.3mm 管壁厚。
HT1_001_UPPER_LINER_ID_MM = HT1_001_UPPER_LINER_OD_MM - 2.0 * HT1_001_UPPER_LINER_WALL_THICKNESS_MM
HT1_001_LOWER_CENTRALIZER_COUNT = 24  # 139.7mm 下段整体式扶正器数量。
HT1_001_UPPER_CENTRALIZER_COUNT = 77  # 168.3mm 上段整体式扶正器数量。
HT1_001_CENTRALIZER_COUNT = HT1_001_LOWER_CENTRALIZER_COUNT + HT1_001_UPPER_CENTRALIZER_COUNT


def _equivalent_hole_diameter_mm(actual_hole_mm: float, actual_od_mm: float, reference_od_mm: float) -> float:
    """在固定 reference_od 条件下，构造与原始环空面积等价的井眼直径。"""

    # 面积守恒：D_eq² - OD_ref² = D_actual² - OD_actual²。
    area_term = actual_hole_mm**2 - actual_od_mm**2 + reference_od_mm**2
    return math.sqrt(max(area_term, reference_od_mm**2))


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


# 鞋口滞后体积估算：地面至悬挂器沿用同上部套管代理，尾管段使用 HT1-001 明确双径内径。
# 该体积只用于 field_order_realistic 边界的到鞋延迟，不代表真实单一管柱内径。
HT1_001_SURFACE_TO_HANGER_EFFECTIVE_ID_MM = math.sqrt(4.0 * 52.0 / (math.pi * 7868.0)) * 1000.0
HT1_001_SHOE_LAG_VOLUME_M3 = (
    _pipe_volume_m3(HT1_001_HANGER_MD_M, HT1_001_SURFACE_TO_HANGER_EFFECTIVE_ID_MM)
    + _pipe_volume_m3(
        HT1_001_UPPER_SECTION_BOTTOM_MD_M - HT1_001_HANGER_MD_M,
        HT1_001_UPPER_LINER_ID_MM,
    )
    + _pipe_volume_m3(
        HT1_001_BOTTOM_MD_M - HT1_001_UPPER_SECTION_BOTTOM_MD_M,
        HT1_001_LOWER_LINER_ID_MM,
    )
)
HT1_001_LINER_ID_MM = HT1_001_LOWER_LINER_ID_MM

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
    well_spec = WellSpec(
        well_name=HT1_001_WELL_NAME,
        top_md_m=HT1_001_TOP_MD_M,
        bottom_md_m=HT1_001_BOTTOM_MD_M,
        shoe_md_m=HT1_001_SHOE_MD_M,
        hanger_md_m=HT1_001_HANGER_MD_M,
        casing_id_mm=HT1_001_CASING_ID_MM,
        liner_od_mm=HT1_001_LOWER_LINER_OD_MM,
        liner_id_mm=HT1_001_LINER_ID_MM,
        hole_diameter_profile=_depth_points(
            (
                # HT1-001 上段按等效井径保留 168.3mm 尾管环空面积。
                (HT1_001_TOP_MD_M, HT1_001_UPPER_HOLE_DIAMETER_MM),
                (6000.0, HT1_001_UPPER_HOLE_DIAMETER_MM),
                (6600.0, HT1_001_UPPER_HOLE_DIAMETER_MM),
                (HT1_001_UPPER_SECTION_BOTTOM_MD_M - 1.0, HT1_001_UPPER_HOLE_DIAMETER_MM),
                # 7174.938m 变径位置做平滑过渡，过渡到 HT1-001 下段 229.46mm 平均井径。
                (
                    HT1_001_UPPER_SECTION_BOTTOM_MD_M,
                    0.5 * (HT1_001_UPPER_HOLE_DIAMETER_MM + HT1_001_LOWER_HOLE_DIAMETER_MM),
                ),
                (7400.0, HT1_001_LOWER_HOLE_DIAMETER_MM),
                (7500.0, HT1_001_LOWER_HOLE_DIAMETER_MM),
                (HT1_001_BOTTOM_MD_M, HT1_001_LOWER_HOLE_DIAMETER_MM),
            )
        ),
        inclination_profile=_depth_points(
            (
                # HT1-001 暂无专属井斜剖面，沿用呼探1井代理模式由 2° 缓增至约 7–8°。
                (HT1_001_TOP_MD_M, 2.0),
                (6000.0, 3.4),
                (6600.0, 5.2),
                (HT1_001_UPPER_SECTION_BOTTOM_MD_M - 1.0, 6.6),
                (HT1_001_UPPER_SECTION_BOTTOM_MD_M, 6.8),
                (7400.0, 7.6),
                (7500.0, 7.8),
                (7710.0, 7.3),
                (HT1_001_BOTTOM_MD_M, 7.0),
            )
        ),
        standoff_profile=_depth_points(
            (
                # 101只整体式扶正器代理居中度：总体控制在 0.55–0.75 范围。
                (HT1_001_TOP_MD_M, 0.64),
                (6000.0, 0.66),
                (6600.0, 0.62),
                (HT1_001_UPPER_SECTION_BOTTOM_MD_M - 1.0, 0.60),
                (HT1_001_UPPER_SECTION_BOTTOM_MD_M, 0.58),
                (7400.0, 0.70),
                (7500.0, 0.72),
                (7710.0, 0.62),
                (HT1_001_BOTTOM_MD_M, 0.58),
            )
        ),
        evaluation_windows=(
            EvaluationWindow(name="CBL评价井段", top_md_m=5700.0, bottom_md_m=7710.0, window_type="cbl"),
            EvaluationWindow(name="目标层段", top_md_m=7400.0, bottom_md_m=7500.0, window_type="target"),
        ),
        reference_root=resolved_reference_root,
        notes=(
            "呼探1-001井（HT1-001）不是呼探1井，本加载器使用 HT1-001 专属井眼、流体和施工程序参数。",
            "上段 168.3mm 尾管按等效井眼直径保面积处理，下段 139.7mm 尾管为模型目标井段。",
            "下段平均井径采用 229.46mm；上段尾管内径采用明确值 138.9mm，下段尾管内径采用明确值 108.04mm。",
            f"鞋口滞后体积估算为 {HT1_001_SHOE_LAG_VOLUME_M3:.2f}m³，WellSpec.liner_id_mm 使用下段内径 {HT1_001_LINER_ID_MM:.2f}mm。",
            "CBL评价井段 5700–7710m 与目标层段 7400–7500m 为首版暂定，后续应随 CBL 数据确认调整。",
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
            PumpingScheduleStep("井浆快替", "井浆", HT1_001_WELL_MUD_FAST_VOLUME_M3, HT1_001_WELL_MUD_FAST_RATE_M3_MIN, remarks="井浆快替 35m³@1.0m³/min。"),
            PumpingScheduleStep("井浆慢替", "井浆", HT1_001_WELL_MUD_SLOW_VOLUME_M3, HT1_001_WELL_MUD_SLOW_RATE_M3_MIN, remarks="井浆慢替 18.7m³@0.6m³/min。"),
        ),
        notes=(
            "施工顺序按 HT1-001 现场数据：平衡液→隔离液→领浆→中间浆→尾浆→压塞液→四段替浆并含基液步骤。",
            "替浆总量 91.7m³ = 25 + 10 + 3 + 35 + 18.7m³，排量从 1.4 递减至 0.6m³/min。",
            "压塞液保留在 PumpingSchedule 中用于管内时序占位；环空入口分相映射时按 mud 相处理。",
        ),
    )

    # 验证资料路径指向 HT1-001 数据包；本加载器首版不在运行时解析 JSON。
    validation_data = ValidationData(
        job_report_path=resolved_reference_root / "提取数据" / "呼探1-001井_固井顶替模型数据包.json",
        notes=(
            "呼探1-001首版加载器不读取 JSON 文件，仅把数据包中的现场参数固化为模块常量。",
            "CBL 评价井段和目标层段为暂定窗口，后续需结合 CBL 原始数据更新 ValidationData 路径。",
        ),
    )
    return well_spec, fluids, schedule, validation_data


def build_ht1_001_annulus_inlet_provider(
    schedule: PumpingSchedule,
    fluids: tuple[FluidSpec, ...],
    annulus_boundary_mode: str = "field_order_realistic",
    *,
    split_cement_phases: bool = False,
) -> Callable[[float], AnnulusInletState]:
    """为呼探1-001井（HT1-001）尾管段构建环空入口边界提供器。

    支持两类口径：
    1. field_order_realistic：按估算鞋口滞后体积修正地面施工流体到达鞋口的顺序；
    2. sustained_tail / volume_limited / tail_then_mud：保留 legacy 对比模式。
    """

    # 建立流体名称到角色的索引，供环空入口相分数映射使用。
    role_by_name: dict[str, FluidRole] = {fluid.name: fluid.role for fluid in fluids}
    if not any(fluid.role == FluidRole.TAIL for fluid in fluids):
        raise ValueError("HT1-001 边界提供器需要尾浆流体")

    # legacy 模式只跟踪前置液和水泥浆入环空，替浆期按固定口径处理。
    annulus_entry_steps = tuple(
        step
        for step in schedule.steps
        if role_by_name.get(step.fluid_name, FluidRole.MUD)
        in {FluidRole.WASH, FluidRole.SPACER, FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}
    )

    def _surface_state(time_s: float) -> tuple[float, float, str]:
        """返回地面累计注入体积、当前排量和当前施工阶段名。"""

        # 按地面施工顺序累加体积，定位给定时刻所在步骤。
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

        # 只在前置液和水泥浆体积范围内寻找环空入口阶段。
        cumulative_m3 = 0.0
        for step in annulus_entry_steps:
            cumulative_m3 += step.volume_m3
            if arrival_volume_m3 < cumulative_m3 - 1e-12:
                return step
        return None

    def _surface_step_by_arrival_volume(arrival_volume_m3: float) -> PumpingScheduleStep | None:
        """按鞋口滞后后的到达体积定位所有地面泵注流体阶段。"""

        # field_order_realistic 口径跟踪所有地面注入流体到达鞋口后的实际顺序。
        cumulative_m3 = 0.0
        for step in schedule.steps:
            cumulative_m3 += step.volume_m3
            if arrival_volume_m3 < cumulative_m3 - 1e-12:
                return step
        return None

    def _provider(time_s: float) -> AnnulusInletState:
        # 先计算地面状态，再扣除 HT1-001 鞋口滞后体积获得到鞋体积。
        surface_volume_m3, flow_rate_m3_s, surface_stage_name = _surface_state(time_s)
        arrival_volume_m3 = max(surface_volume_m3 - HT1_001_SHOE_LAG_VOLUME_M3, 0.0)

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
                f"{arrival_step.step_name}（{HT1_001_SHOE_LAG_VOLUME_M3:.1f}m³鞋口滞后修正）",
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
            f"{annulus_step.step_name}（{HT1_001_SHOE_LAG_VOLUME_M3:.1f}m³鞋口滞后修正）",
            _phase_fractions_for_fluid(
                annulus_step.fluid_name,
                role_by_name,
                split_cement_phases=split_cement_phases,
            ),
        )

    return _provider
