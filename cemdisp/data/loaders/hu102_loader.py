"""
呼102尾管段标准数据加载器

本模块实现呼102井139.70mm尾管段固井的标准数据加载功能。

主要功能：
- 从井径/井斜CSV文件读取剖面数据
- 构建井筒几何参数（井段范围、套管尺寸、偏心度等）
- 定义钻井液、替浆液、尾管水泥浆，以及可选冲洗液/隔离液的物性参数
- 构建现场记录施工日程（默认严格按现场：尾浆注入+替浆液推进两步）
- 提供 legacy 环空入口边界状态提供器（仅用于旧模型对比）

物理参数说明：
- 井段范围: 6823.10m - 7735.00m
- 尾管尺寸: 139.70mm OD, 108.10mm ID (考虑壁厚后)
- 水泥浆: 35t, 密度2.10g/cm³；幂律流变 n=0.722, K=0.684 当前仍为 legacy 占位
- 替浆液: 74m³, 密度2.02g/cm³, 排量0.378m³/min
- 钻井液（环空初始液）: 密度2.02g/cm³, Bingham PV=80mPa·s, YP=15Pa

现场记录来源（10042.xlsx Row 26, 2022-11-22）：
- 注水泥35.00t, 水泥浆平均密度2.10g/cm³, 替浆液密度2.02g/cm³, 井液74.00m³
- 泵注时间：2022-11-21 17:00–21:00（4小时）
- 现场记录中无冲洗液/隔离液/领浆的注入量（方案A：按现场记录）

可选补充流体（0708邻井代理，严格现场模式下不注入）：
- 冲洗液(WASH): ρ=1880, PV=0.025, YP=1.5 (呼103邻井)
- 隔离液(SPACER): ρ=1850, PV=0.035, YP=8 (呼103邻井)

legacy 边界模式选项：
- "sustained_tail": 替浆期间环空入口保持尾浆（默认）
- "volume_limited": 管内推进期间环空入口保持尾浆，排量设为0
- "tail_then_mud": 替浆期结束后环空入口切换为替浆液
"""

from __future__ import annotations

import warnings

import csv
from collections.abc import Callable, Iterable
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
DEFAULT_REFERENCE_ROOT = PROJECT_ROOT / "参考文档" / "呼102"
DEFAULT_CALIPER_CSV = PROJECT_ROOT / "hu102model" / "hu102_tail_caliper_inclination.csv"

# 呼102尾管段井段参数（单位：m, mm）
HU102_TOP_MD_M = 6823.10        # 井段顶部测深
HU102_BOTTOM_MD_M = 7735.00     # 井段底部测深
HU102_SHOE_MD_M = 7735.00       # 套管鞋深度
HU102_HANGER_MD_M = 6823.10     # 悬挂器位置深度
HU102_CASING_ID_MM = 219.10     # 套管内径
HU102_LINER_OD_MM = 139.70      # 尾管外径
HU102_LINER_WALL_THICKNESS_MM = 15.80  # 尾管壁厚
HU102_LINER_ID_MM = HU102_LINER_OD_MM - 2.0 * HU102_LINER_WALL_THICKNESS_MM  # 尾管内径

# 呼102施工参数
HU102_CEMENT_MASS_T = 35.0              # 水泥浆质量
HU102_CEMENT_DENSITY_KG_M3 = 2100.0     # 水泥浆密度
HU102_DISPLACEMENT_VOLUME_M3 = 74.0    # 替浆体积
HU102_DISPLACEMENT_DENSITY_KG_M3 = 2020.0  # 钻井液密度（替浆用）
# 按现场 17:00–21:00 累计 90.67m³ / 240min 校正平均排量，避免原1.30m³/min压缩顶替过程。
HU102_RATE_M3_MIN = 0.378               # 泵注排量

# 呼102流变参数 — 钻井液/替浆液/水泥浆
HU102_MUD_PV_PA_S = 0.080            # 环空初始钻井液塑性粘度（文献暂定）
HU102_MUD_YP_PA = 15.0               # 环空初始钻井液屈服值（文献暂定）
HU102_DISPLACEMENT_PV_PA_S = 0.080   # 替浆液塑性粘度（与钻井液一致）
HU102_DISPLACEMENT_YP_PA = 15.0      # 替浆液屈服值（与钻井液一致）
HU102_CEMENT_POWER_LAW_N = 0.722     # 水泥浆流性指数（幂律，占位值；本井主作业未见实测 n）
HU102_CEMENT_CONSISTENCY_K = 0.684   # 水泥浆稠度系数（幂律，占位值；本井主作业未见实测 K）

# 呼102前置液/隔离液参数 — 基于呼探1-002邻井同口径139.7mm尾管数据
# 数据来源：
#   - 呼探1-002 139.7mm尾管：隔离液 2.05g/cm³(设计)/2.10g/cm³(现场)，15m³，冲洗效率97.7%
#   - 呼探1-002 数据抽取报告：化验报告+技术总结+作业史
#   - 注：Hu102尾管主作业日报(10042.xlsx)未找到隔离液记录，此处使用邻井代理
HU102_WASH_DENSITY_KG_M3 = 2050.0    # 平衡液/冲洗液密度（呼探1-002邻井代理，与隔离液同体系）
HU102_WASH_PV_PA_S = 0.035           # 平衡液/冲洗液塑性粘度（呼探1-002邻井代理）
HU102_WASH_YP_PA = 8.0               # 平衡液/冲洗液屈服值（呼探1-002邻井代理）
HU102_WASH_VOLUME_M3 = 10.0          # 平衡液/冲洗液设计体积（呼探1-002邻井代理）
# 合成 FLUSHER（冲洗液）参数：用于验证 mud-spacer-flusher-cement 序列可表达，
# 本井现场记录无真实冲洗液，故采用邻井呼103代理值（ρ=1880, PV=0.025, YP=1.5）。
HU102_FLUSHER_DENSITY_KG_M3 = 1880.0
HU102_FLUSHER_PV_PA_S = 0.025
HU102_FLUSHER_YP_PA = 1.5
HU102_FLUSHER_VOLUME_M3 = 5.0
HU102_SPACER_DENSITY_KG_M3 = 2050.0  # 驱油隔离液密度（呼探1-002设计值2.05g/cm³）
HU102_SPACER_PV_PA_S = 0.035         # 驱油隔离液塑性粘度（呼探1-002邻井代理）
HU102_SPACER_YP_PA = 8.0             # 驱油隔离液屈服值（呼探1-002邻井代理）
HU102_SPACER_VOLUME_M3 = 15.0        # 驱油隔离液体积（呼探1-002现场记录15m³）


def _read_profile_rows(caliper_csv_path: Path) -> tuple[tuple[float, float, float], ...]:
    with caliper_csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows: list[tuple[float, float, float]] = []
        for row in csv.DictReader(handle):
            rows.append(
                (
                    float(row["depth_md_m"]),
                    float(row["hole_diameter_mm"]),
                    float(row["inclination_deg"]),
                )
            )
    if not rows:
        raise ValueError(f"Hu102 井径/井斜 CSV 为空: {caliper_csv_path}")
    return tuple(sorted(rows))


def _standoff_value(depth_md_m: float, hole_diameter_mm: float, liner_od_mm: float) -> float:
    """计算居中度(standoff)剖面。

    基于邻井呼探1-002的扶正器数据估算：
    - 呼探1-002: 95只整体式弹扶，目的层44m间距，非目的层55m间距
    - 呼102井段较短(911.9m)，采用相似间距策略

    居中度定义：standoff = 1 - 偏心度
    - standoff = 1.0: 完全居中
    - standoff = 0.0: 完全偏心
    """
    # 基于邻井数据的居中度估算
    # 目的层段(7405-7540m): 44m间距 → standoff ≈ 0.70
    # 非目的层段: 55m间距 → standoff ≈ 0.65
    standoff = 0.65  # 默认值（非目的层）
    if 7405.00 <= depth_md_m <= 7540.00:
        standoff = 0.70  # 目的层段（油气水层）
    elif depth_md_m > 7540.00:
        standoff = 0.68  # 底部段（略好于非目的层）
    # 间隙修正因子：间隙越大，居中度越好
    clearance_mm = max(hole_diameter_mm - liner_od_mm, 5.0)
    standoff *= min(max(clearance_mm / 70.0, 0.55), 1.0)
    return min(max(standoff, 0.30), 0.85)


def _build_hole_profile(profile_rows: Iterable[tuple[float, float, float]]) -> tuple[DepthValuePoint, ...]:
    return tuple(DepthValuePoint(depth_md_m=depth, value=hole_diameter) for depth, hole_diameter, _ in profile_rows)


def _build_inclination_profile(profile_rows: Iterable[tuple[float, float, float]]) -> tuple[DepthValuePoint, ...]:
    return tuple(DepthValuePoint(depth_md_m=depth, value=inclination) for depth, _, inclination in profile_rows)


def _build_standoff_profile(profile_rows: Iterable[tuple[float, float, float]]) -> tuple[DepthValuePoint, ...]:
    return tuple(
        DepthValuePoint(
            depth_md_m=depth,
            value=_standoff_value(depth, hole_diameter, HU102_LINER_OD_MM),
        )
        for depth, hole_diameter, _ in profile_rows
    )


def load_hu102_tailpipe(
    *,
    caliper_csv_path: Path | None = None,
    reference_root: Path | None = None,
    include_wash_spacer: bool = False,  # 严格按呼102主作业实录：默认不注入邻井代理前置液/隔离液
) -> tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载呼102尾管段首版模型输入。

    Args:
        caliper_csv_path: 可选井径/井斜 CSV 路径。
        reference_root: 可选参考资料根目录。
        include_wash_spacer: 是否把 0708 邻井代理的冲洗液/隔离液步骤加入泵注程序。
            默认为 False，严格按 10042.xlsx 主作业记录仅保留尾浆+替浆两步。

    Returns:
        井筒参数、流体参数、泵注程序与验证资料路径。
    """

    resolved_reference_root = reference_root or DEFAULT_REFERENCE_ROOT
    resolved_caliper_csv_path = caliper_csv_path or DEFAULT_CALIPER_CSV

    profile_rows = _read_profile_rows(resolved_caliper_csv_path)
    well_spec = WellSpec(
        well_name="呼102",
        top_md_m=HU102_TOP_MD_M,
        bottom_md_m=HU102_BOTTOM_MD_M,
        shoe_md_m=HU102_SHOE_MD_M,
        hanger_md_m=HU102_HANGER_MD_M,
        casing_id_mm=HU102_CASING_ID_MM,
        liner_od_mm=HU102_LINER_OD_MM,
        liner_id_mm=HU102_LINER_ID_MM,
        hole_diameter_profile=_build_hole_profile(profile_rows),
        inclination_profile=_build_inclination_profile(profile_rows),
        standoff_profile=_build_standoff_profile(profile_rows),
        evaluation_windows=(
            EvaluationWindow(name="CBL评价井段", top_md_m=6840.0, bottom_md_m=7665.0, window_type="cbl"),
            EvaluationWindow(name="目标层段一", top_md_m=7405.0, bottom_md_m=7480.0, window_type="target"),
            EvaluationWindow(name="目标层段二", top_md_m=7502.0, bottom_md_m=7540.0, window_type="target"),
        ),
        reference_root=resolved_reference_root,
        notes=(
            "仅对应呼102井 139.70mm 尾管段固井，不含回接固井与其他套管段。",
            "7120–7735m 井径/井斜取自 20215.xlsx Sheet4 派生 CSV；6823.10–7119.80m 双层套管段按等效井径处理。",
        ),
    )

    fluids = (
        FluidSpec(
            name="钻井液",
            role=FluidRole.MUD,
            density_kg_m3=HU102_DISPLACEMENT_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU102_MUD_PV_PA_S,
            yield_stress_pa=HU102_MUD_YP_PA,
        ),
        FluidSpec(
            name="替浆液",
            role=FluidRole.DISPLACEMENT,
            density_kg_m3=HU102_DISPLACEMENT_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU102_DISPLACEMENT_PV_PA_S,
            yield_stress_pa=HU102_DISPLACEMENT_YP_PA,
        ),
        FluidSpec(
            name="尾管水泥浆",
            role=FluidRole.TAIL,
            density_kg_m3=HU102_CEMENT_DENSITY_KG_M3,
            rheology_model=RheologyModel.POWER_LAW,
            power_law_n=HU102_CEMENT_POWER_LAW_N,
            consistency_k=HU102_CEMENT_CONSISTENCY_K,
        ),
        FluidSpec(
            name="冲洗液",
            role=FluidRole.WASH,
            density_kg_m3=HU102_WASH_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU102_WASH_PV_PA_S,
            yield_stress_pa=HU102_WASH_YP_PA,
        ),
        FluidSpec(
            name="隔离液",
            role=FluidRole.SPACER,
            density_kg_m3=HU102_SPACER_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU102_SPACER_PV_PA_S,
            yield_stress_pa=HU102_SPACER_YP_PA,
        ),
        FluidSpec(
            name="冲洗液（FLUSHER）",
            role=FluidRole.FLUSHER,
            density_kg_m3=HU102_FLUSHER_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU102_FLUSHER_PV_PA_S,
            yield_stress_pa=HU102_FLUSHER_YP_PA,
        ),
    )

    # 前置液/隔离液步骤：仅在显式要求敏感性分析时加入。
    # 呼102主作业日报未记录该两类流体，严格现场模式默认不使用邻井代理值。
    # 数据来源：Hu102二次技套(20258.doc)、Hu103回接(20314.doc/20323.doc)、呼探1-002
    optional_front_steps = ()
    if include_wash_spacer:
        optional_front_steps = (
            PumpingScheduleStep(
                step_name="注入平衡液",
                fluid_name="冲洗液",  # 角色映射为WASH，三相模型中归入隔离液相
                volume_m3=HU102_WASH_VOLUME_M3,
                rate_m3_min=HU102_RATE_M3_MIN,
                remarks=f"平衡液/冲洗液 {HU102_WASH_VOLUME_M3}m³，密度{HU102_WASH_DENSITY_KG_M3/1000:.2f}g/cm³（呼103邻井代理）。",
            ),
            PumpingScheduleStep(
                step_name="注入驱油隔离液",
                fluid_name="隔离液",
                volume_m3=HU102_SPACER_VOLUME_M3,
                rate_m3_min=HU102_RATE_M3_MIN,
                remarks=f"驱油隔离液 {HU102_SPACER_VOLUME_M3}m³，密度{HU102_SPACER_DENSITY_KG_M3/1000:.2f}g/cm³（Hu102二次技套/呼103邻井代理）。",
            ),
        )

    # 合成 FLUSHER 步骤：与前置液/隔离液同属邻井代理敏感性输入，仅 include_wash_spacer=True 时注入，
    # 位于隔离液之后、水泥浆之前，验证 mud-spacer-flusher-cement 序列可表达。
    flusher_step = PumpingScheduleStep(
        step_name="注入冲洗液（FLUSHER）",
        fluid_name="冲洗液（FLUSHER）",
        volume_m3=HU102_FLUSHER_VOLUME_M3,
        rate_m3_min=HU102_RATE_M3_MIN,
        remarks=f"合成冲洗液（FLUSHER）{HU102_FLUSHER_VOLUME_M3}m³，密度{HU102_FLUSHER_DENSITY_KG_M3/1000:.2f}g/cm³（呼103邻井代理，验证序列可表达）。",
    )

    front_steps = optional_front_steps
    if include_wash_spacer:
        front_steps = front_steps + (flusher_step,)

    schedule = PumpingSchedule(
        steps=front_steps + (
            PumpingScheduleStep(
                step_name="注入尾管水泥浆",
                fluid_name="尾管水泥浆",
                volume_m3=HU102_CEMENT_MASS_T / (HU102_CEMENT_DENSITY_KG_M3 / 1000.0),
                rate_m3_min=HU102_RATE_M3_MIN,
                remarks="基于 35t 与 2.10g/cm3 换算 ≈ 16.67m3。",
            ),
            PumpingScheduleStep(
                step_name="替浆液推进",
                fluid_name="替浆液",
                volume_m3=HU102_DISPLACEMENT_VOLUME_M3,
                rate_m3_min=HU102_RATE_M3_MIN,
                remarks="主作业直接记录替浆量 74m3，替浆液密度 2.02g/cm3。",
            ),
        ),
        notes=(
            "按现场记录（10042.xlsx Row 26）：尾浆+替浆液两步为主程序。",
            "严格现场模式（include_wash_spacer=False）默认不注入前置液/隔离液/合成FLUSHER。",
            "10042.xlsx 主作业记录未见冲洗液、隔离液、领浆或独立FLUSHER注入量。",
            "平衡液/隔离液/合成FLUSHER参数仅保留为 include_wash_spacer=True 时的邻井代理敏感性输入，不作为呼102现场实录。",
            "include_wash_spacer=True 时才加入邻井代理的平衡液(10m³)、驱油隔离液(15m³)及合成FLUSHER(5m³)。",
        ),
    )

    validation_data = ValidationData(
        cbl_summary_path=resolved_reference_root / "1004" / "10041" / "100413.PDF",
        cbl_pass_rate=0.6665,
        job_report_path=resolved_reference_root / "1004" / "10042.xlsx",
        pump_pressure_series_path=resolved_reference_root / "1004" / "100492.xlsx",
        notes=(
            "100413.PDF 给出 CBL 合格率 66.65%，评价井段 6840–7665m。",
            "10042.xlsx 提供尾管固井主作业水泥浆 35t、平均密度 2.10g/cm3、替浆量 74m3。",
            "钻井液与水泥浆流变参数仍为首版暂定值，需后续继续用 0708 或文献补强。",
        ),
    )
    return well_spec, fluids, schedule, validation_data


def export_hu102_sync_card_markdown(output_dir: Path) -> Path:
    """导出 呼102 同步画像卡为 Markdown 文件。"""

    well_spec, fluids, schedule, _ = load_hu102_tailpipe()
    solver = CasingFlowSolver()
    result = solver.run(well_spec, fluids, schedule)
    provenance = WELL_PROVENANCE.get("呼102")
    if provenance is None:
        from cemdisp.data.provenance import _WELL_NAME_ALIASES
        canonical = _WELL_NAME_ALIASES.get("呼102", "呼102")
        provenance = WELL_PROVENANCE[canonical]
    sync_card = build_sync_card("呼102", result.shoe_timeline, provenance)

    output_path = output_dir / ("呼102_同步画像卡.md")
    lines = [
        "# 呼102 同步画像卡",
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
