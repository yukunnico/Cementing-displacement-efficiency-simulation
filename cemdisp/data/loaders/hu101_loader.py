"""
呼101尾管段标准数据加载器

本模块把 legacy 呼101脚本中的本井资料整理为 cemdisp 标准输入结构，
用于 Hu101 runner 复用通用 2D 环空求解、导出和 1D-2D 边界桥接流程。
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
DEFAULT_REFERENCE_ROOT = PROJECT_ROOT / "参考文档" / "呼101"

# 呼101尾管段井身结构参数，来自 hu101 legacy 模型说明。
HU101_TOP_MD_M = 5400.0
HU101_BOTTOM_MD_M = 7868.0
HU101_SHOE_MD_M = 7868.0
HU101_HANGER_MD_M = 5407.46
HU101_TECH_CASING_EQUIV_ID_MM = 273.10
HU101_UPPER_SECTION_BOTTOM_MD_M = 6796.0
HU101_UPPER_ACTUAL_HOLE_DIAMETER_MM = 260.35
HU101_UPPER_ACTUAL_LINER_OD_MM = 168.30
HU101_LOWER_HOLE_DIAMETER_MM = 215.90
HU101_LOWER_LINER_OD_MM = 139.70
HU101_LINER_WALL_THICKNESS_MM = 15.80
HU101_SHOE_LAG_VOLUME_M3 = 52.0


def _equivalent_hole_diameter_mm(actual_hole_mm: float, actual_od_mm: float, reference_od_mm: float) -> float:
    """在固定 reference_od 条件下，构造与原始环空面积等价的井眼直径。"""

    area_term = actual_hole_mm**2 - actual_od_mm**2 + reference_od_mm**2
    return math.sqrt(max(area_term, reference_od_mm**2))


HU101_UPPER_HOLE_DIAMETER_MM = _equivalent_hole_diameter_mm(
    actual_hole_mm=HU101_UPPER_ACTUAL_HOLE_DIAMETER_MM,
    actual_od_mm=HU101_UPPER_ACTUAL_LINER_OD_MM,
    reference_od_mm=HU101_LOWER_LINER_OD_MM,
)
HU101_LINER_ID_MM = math.sqrt(4.0 * HU101_SHOE_LAG_VOLUME_M3 / (math.pi * HU101_SHOE_MD_M)) * 1000.0

# 呼101现场施工与流体参数，来自 h101 抽取报告和 legacy 现场顺序脚本。
HU101_LEAD_VOLUME_M3 = 47.0
HU101_TAIL_VOLUME_M3 = 23.0
HU101_MUD_DENSITY_KG_M3 = 1960.0
HU101_BALANCE_DENSITY_KG_M3 = 1850.0
HU101_SPACER_DENSITY_KG_M3 = 2000.0
HU101_LEAD_DENSITY_KG_M3 = 2100.0
HU101_TAIL_DENSITY_KG_M3 = 1900.0
HU101_MUD_PV_PA_S = 0.058
HU101_MUD_YP_PA = 9.2  # 现场数据：65℃下YP=9.2 Pa（来自尾管设计）
HU101_BALANCE_PV_PA_S = 0.030
HU101_BALANCE_YP_PA = 3.0
HU101_SPACER_PV_PA_S = 0.030
HU101_SPACER_YP_PA = 5.0
HU101_LEAD_POWER_LAW_N = 0.719
HU101_LEAD_CONSISTENCY_K = 0.815
HU101_TAIL_POWER_LAW_N = 0.722
HU101_TAIL_CONSISTENCY_K = 0.684


def _depth_points(values: tuple[tuple[float, float], ...]) -> tuple[DepthValuePoint, ...]:
    """把测深-数值元组转换为 WellSpec 使用的剖面点。"""

    return tuple(DepthValuePoint(depth_md_m=depth, value=value) for depth, value in values)


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


def load_hu101_tailpipe(
    *,
    reference_root: Path | None = None,
) -> tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载呼101尾管段标准模型输入。"""

    resolved_reference_root = reference_root or DEFAULT_REFERENCE_ROOT
    well_spec = WellSpec(
        well_name="呼101",
        top_md_m=HU101_TOP_MD_M,
        bottom_md_m=HU101_BOTTOM_MD_M,
        shoe_md_m=HU101_SHOE_MD_M,
        hanger_md_m=HU101_HANGER_MD_M,
        casing_id_mm=HU101_TECH_CASING_EQUIV_ID_MM,
        liner_od_mm=HU101_LOWER_LINER_OD_MM,
        liner_id_mm=HU101_LINER_ID_MM,
        hole_diameter_profile=_depth_points(
            (
                (HU101_TOP_MD_M, HU101_UPPER_HOLE_DIAMETER_MM),
                (5700.0, HU101_UPPER_HOLE_DIAMETER_MM),
                (6050.0, HU101_UPPER_HOLE_DIAMETER_MM),
                (6210.0, HU101_UPPER_HOLE_DIAMETER_MM),
                (6400.0, HU101_UPPER_HOLE_DIAMETER_MM),
                (HU101_UPPER_SECTION_BOTTOM_MD_M - 1.0, HU101_UPPER_HOLE_DIAMETER_MM),
                (HU101_UPPER_SECTION_BOTTOM_MD_M, 0.5 * (HU101_UPPER_HOLE_DIAMETER_MM + HU101_LOWER_HOLE_DIAMETER_MM)),
                (7400.0, HU101_LOWER_HOLE_DIAMETER_MM),
                (7735.0, HU101_LOWER_HOLE_DIAMETER_MM),
                (HU101_BOTTOM_MD_M, HU101_LOWER_HOLE_DIAMETER_MM),
            )
        ),
        inclination_profile=_depth_points(
            (
                (HU101_TOP_MD_M, 0.6),
                (6100.0, 1.1),
                (6796.0, 1.4),
                (7400.0, 1.7),
                (HU101_BOTTOM_MD_M, 1.9),
            )
        ),
        standoff_profile=_depth_points(
            (
                (HU101_TOP_MD_M, 0.45),
                (6100.0, 0.38),
                (6796.0, 0.44),
                (7200.0, 0.48),
                (7600.0, 0.42),
                (HU101_BOTTOM_MD_M, 0.46),
            )
        ),
        evaluation_windows=(
            EvaluationWindow(name="CBL评价井段", top_md_m=5700.0, bottom_md_m=7810.0, window_type="cbl"),
            EvaluationWindow(name="目标层段", top_md_m=7492.0, bottom_md_m=7735.0, window_type="target"),
        ),
        reference_root=resolved_reference_root,
        notes=(
            "呼101上部168.3mm+下部139.7mm变径井段，当前用等效上部井眼直径保面积近似，以兼容通用单外径求解器。",
            "liner_id_mm 使用52m³鞋口滞后反推的等效内径，仅服务于1D到鞋时间口径，不代表真实单一尾管内径。",
            "井径、井斜和居中度剖面来自 hu101 legacy D2DGA 模型及改进模型的名义/退化剖面。",
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
            PumpingScheduleStep("注平衡液", "平衡液", 25.0, 1.2, remarks="呼101现场抽取：25m³。"),
            PumpingScheduleStep("注驱油隔离液", "驱油隔离液", 25.0, 1.2, remarks="呼101现场抽取：25m³。"),
            PumpingScheduleStep("注领浆", "领浆", HU101_LEAD_VOLUME_M3, 1.2, remarks="呼101现场抽取：47m³。"),
            PumpingScheduleStep("注尾浆", "尾浆", HU101_TAIL_VOLUME_M3, 1.2, remarks="呼101现场抽取：23m³。"),
            PumpingScheduleStep("注后置液(管内)", "井浆", 2.0, 0.6, remarks="仅作为管内压塞/占位流体。"),
            PumpingScheduleStep("注轻泥浆", "轻泥浆", 26.0, 1.5, remarks="按现场轻泥浆排量建模。"),
            PumpingScheduleStep("注中置液", "中置液", 10.0, 1.2, remarks="按主替浆排量建模。"),
            PumpingScheduleStep("井浆快替", "井浆", 40.0, 1.0, remarks="设计表40m³@1.0m³/min。"),
            PumpingScheduleStep("井浆慢替", "井浆", 23.4, 0.55, remarks="补足现场总替量101.4m³。"),
        ),
        notes=(
            "现场顺序：平衡液→驱油隔离液→领浆→尾浆→后置液→轻泥浆→中置液→井浆快替/慢替。",
            "替浆合计101.4m³，其中井浆慢替用于补足现场总替量。",
        ),
    )

    validation_data = ValidationData(
        cbl_summary_path=resolved_reference_root / "1003" / "100312.PDF",
        job_report_path=resolved_reference_root / "提取数据" / "h101_data_extraction_report.md",
        notes=(
            "100312.PDF CBL测井固井质量合格率62.77%。",
            "CBL评价口径按 legacy 模型的 5700–7810m 区间设置。",
        ),
    )
    return well_spec, fluids, schedule, validation_data


def build_hu101_annulus_inlet_provider(
    schedule: PumpingSchedule,
    fluids: tuple[FluidSpec, ...],
    annulus_boundary_mode: str = "field_order_realistic",
    *,
    split_cement_phases: bool = False,
) -> Callable[[float], AnnulusInletState]:
    """为呼101尾管段构建硬编码环空入口边界提供器。

    支持两类口径：
    1. sustained_tail / volume_limited / tail_then_mud：保持 legacy 对比口径；
    2. field_order_realistic：按 52m³ 鞋口滞后，把后续替浆到鞋后的 mud invasion 单独送入环空。
    """

    warnings.warn(
        "build_hu101_annulus_inlet_provider 已废弃，请改用 CasingFlowSolver + build_coupled_annulus_inlet_provider",
        DeprecationWarning,
        stacklevel=2,
    )
    role_by_name = {fluid.name: fluid.role for fluid in fluids}
    if not any(fluid.role == FluidRole.TAIL for fluid in fluids):
        raise ValueError("Hu101 边界提供器需要尾浆流体")
    shoe_lag_volume_m3 = HU101_SHOE_LAG_VOLUME_M3
    annulus_entry_steps = tuple(
        step
        for step in schedule.steps
        if role_by_name.get(step.fluid_name, FluidRole.MUD) in {FluidRole.WASH, FluidRole.SPACER, FluidRole.LEAD, FluidRole.TAIL}
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
        """按52m³鞋口滞后后的到达体积定位实际进入环空的前置液/水泥阶段。"""

        cumulative_m3 = 0.0
        for step in annulus_entry_steps:
            cumulative_m3 += step.volume_m3
            if arrival_volume_m3 < cumulative_m3 - 1e-12:
                return step
        return None

    def _surface_step_by_arrival_volume(arrival_volume_m3: float) -> PumpingScheduleStep | None:
        """按52m³鞋口滞后后的到达体积定位所有到鞋流体阶段。"""

        cumulative_m3 = 0.0
        for step in schedule.steps:
            cumulative_m3 += step.volume_m3
            if arrival_volume_m3 < cumulative_m3 - 1e-12:
                return step
        return None

    def _provider(time_s: float) -> AnnulusInletState:
        surface_volume_m3, flow_rate_m3_s, surface_stage_name = _surface_state(time_s)
        arrival_volume_m3 = max(surface_volume_m3 - shoe_lag_volume_m3, 0.0)

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
                f"{arrival_step.step_name}（52m³鞋口滞后修正）",
                _phase_fractions_for_fluid(
                    arrival_step.fluid_name,
                    role_by_name,
                    split_cement_phases=split_cement_phases,
                ),
            )

        annulus_step = _annulus_step_by_arrival_volume(arrival_volume_m3)
        if arrival_volume_m3 <= 0.0 or annulus_step is None:
            if annulus_boundary_mode == "sustained_tail":
                return AnnulusInletState(time_s, flow_rate_m3_s, f"{surface_stage_name}（环空保持尾浆）", _phase_fractions_for_role(FluidRole.TAIL))
            if annulus_boundary_mode == "volume_limited":
                return AnnulusInletState(time_s, 0.0, f"{surface_stage_name}（环空入口保持尾浆）", _phase_fractions_for_role(FluidRole.TAIL))
            if annulus_boundary_mode == "tail_then_mud":
                return AnnulusInletState(time_s, flow_rate_m3_s, f"{surface_stage_name}入环空", _phase_fractions_for_role(FluidRole.MUD))
            raise ValueError(f"Unsupported annulus boundary mode: {annulus_boundary_mode}")

        return AnnulusInletState(
            time_s,
            flow_rate_m3_s,
            f"{annulus_step.step_name}（52m³鞋口滞后修正）",
            _phase_fractions_for_fluid(
                annulus_step.fluid_name,
                role_by_name,
                split_cement_phases=split_cement_phases,
            ),
        )

    return _provider


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
