"""
呼探1-002井（HT1-002）139.7mm完井尾管段标准数据加载器

本模块把呼探1-002井现场数据包中的 139.7mm 完井尾管资料整理为 cemdisp
标准输入结构。与呼探1井不同，本井目标段为单一 139.7mm 尾管，不构造
双径尾管等效几何；上部地面至悬挂器的管内容积只用于鞋口滞后体积估算。
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


def build_hu2_annulus_inlet_provider(
    schedule: PumpingSchedule,
    fluids: tuple[FluidSpec, ...],
    annulus_boundary_mode: str = "field_order_realistic",
    *,
    split_cement_phases: bool = False,
) -> Callable[[float], AnnulusInletState]:
    """为呼探1-002井 139.7mm 尾管段构建环空入口边界提供器。

    支持两类口径：
    1. field_order_realistic：按估算鞋口滞后体积修正地面施工流体到达鞋口的顺序；
    2. sustained_tail / volume_limited / tail_then_mud：保留 legacy 对比模式。
    """

    # 建立流体名称到角色的映射，用于把现场阶段转换为二维环空相分数。
    role_by_name: dict[str, FluidRole] = {fluid.name: fluid.role for fluid in fluids}
    if not any(fluid.role == FluidRole.TAIL for fluid in fluids):
        raise ValueError("Hu2 边界提供器需要尾浆流体")

    # legacy 模式只跟踪前置液和水泥浆入环空，替浆期按固定口径处理。
    annulus_entry_steps = tuple(
        step
        for step in schedule.steps
        if role_by_name.get(step.fluid_name, FluidRole.MUD)
        in {FluidRole.WASH, FluidRole.SPACER, FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}
    )

    def _surface_state(time_s: float) -> tuple[float, float, str]:
        """返回地面累计注入体积、当前排量和当前施工阶段名。"""

        # 按 PumpingSchedule 顺序积分地面注入体积，用于鞋口滞后体积修正。
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

        # 仅在前置液和水泥浆累计体积中查找，用于旧边界模式对比。
        cumulative_m3 = 0.0
        for step in annulus_entry_steps:
            cumulative_m3 += step.volume_m3
            if arrival_volume_m3 < cumulative_m3 - 1e-12:
                return step
        return None

    def _surface_step_by_arrival_volume(arrival_volume_m3: float) -> PumpingScheduleStep | None:
        """按鞋口滞后后的到达体积定位所有地面泵注流体阶段。"""

        # field_order_realistic 使用完整施工顺序，压塞液和替浆液也参与到鞋口到达状态判断。
        cumulative_m3 = 0.0
        for step in schedule.steps:
            cumulative_m3 += step.volume_m3
            if arrival_volume_m3 < cumulative_m3 - 1e-12:
                return step
        return None

    def _provider(time_s: float) -> AnnulusInletState:
        """按给定时间返回环空入口边界状态。"""

        # 先计算地面累计体积，再扣除估算管内容积得到鞋口到达体积。
        surface_volume_m3, flow_rate_m3_s, surface_stage_name = _surface_state(time_s)
        arrival_volume_m3 = max(surface_volume_m3 - HU2_SHOE_LAG_VOLUME_M3, 0.0)

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
                f"{arrival_step.step_name}（{HU2_SHOE_LAG_VOLUME_M3:.1f}m³鞋口滞后修正）",
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
            f"{annulus_step.step_name}（{HU2_SHOE_LAG_VOLUME_M3:.1f}m³鞋口滞后修正）",
            _phase_fractions_for_fluid(
                annulus_step.fluid_name,
                role_by_name,
                split_cement_phases=split_cement_phases,
            ),
        )

    return _provider
