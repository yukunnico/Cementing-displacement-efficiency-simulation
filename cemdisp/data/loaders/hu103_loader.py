"""
呼103尾管段标准数据加载器

本模块实现呼103井 139.70mm 完井尾管段固井的标准数据加载功能。

主要功能：
- 从井径/井斜 CSV 文件读取剖面数据
- 构建井筒几何参数（井段范围、套管尺寸、偏心度等）
- 定义钻井液、平衡液、冲洗液、隔离液、领浆、中间浆、尾浆与替浆液物性参数
- 构建现场设计施工日程（平衡液+冲洗液+隔离液+领浆+中间浆+尾浆+压塞液+四段替浆）
- 提供 legacy 环空入口边界状态提供器（仅用于旧模型对比）

物理参数说明：
- 井段范围: 7330.60m - 7770.00m
- 尾管尺寸: 139.70mm OD, 121.36mm ID（缺实测壁厚，首版暂用9.17mm）
- 上段尾管内径代理: 150.42mm（φ168.30mm 上段尾管 ID 代理值，假设壁厚8.94mm）
- 井眼尺寸: 215.90mm（下段尾管井眼名义尺寸）
- 钻井液（环空初始液）: 密度1.98g/cm³, Bingham PV=54mPa·s, YP=12.5Pa
- 领浆: 91m³, 密度2.05g/cm³, 幂律流变 n=0.82, K=0.67
- 中间浆: 10m³代理, 密度2.05g/cm³, 幂律流变 n=0.76, K=1.11
- 尾浆: 19m³, 密度2.05g/cm³, 幂律流变 n=0.76, K=1.14
- 替浆液: 110m³, 密度1.50g/cm³, 四段排量 1.5/1.0/0.8/0.6m³/min

现场资料来源（呼103井数据包第二版）：
- φ139.70mm 完井尾管顶部/变扣位置 7330.6m，尾管鞋/完钻井深 7770.0m
- CBL 评价井段 7338.0–7712.0m，CBL 合格率 12.06%
- 流体物性来自 model_ready_candidate 与 recommended_inputs
- 中间浆体积与替浆分段体积为首版代理值，已在施工步骤备注中标注

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
DEFAULT_REFERENCE_ROOT = PROJECT_ROOT / "参考文档" / "呼103"
DEFAULT_CALIPER_CSV = PROJECT_ROOT / "cemdisp" / "data" / "loaders" / "hu103_tail_caliper_inclination.csv"

# 呼103尾管段井段参数（单位：m, mm）
HU103_TOP_MD_M = 7330.6          # φ139.70mm 完井尾管顶部/变扣位置
HU103_BOTTOM_MD_M = 7770.0       # 完钻井深/尾管鞋
HU103_SHOE_MD_M = 7770.0         # 套管鞋深度
HU103_HANGER_MD_M = 7330.6       # 悬挂器位置深度
HU103_CASING_ID_MM = 150.42      # φ168.30mm 上段尾管 ID 代理值，假设壁厚8.94mm
HU103_LINER_OD_MM = 139.7        # 完井尾管外径
HU103_LINER_WALL_THICKNESS_MM = 9.17  # 缺实测壁厚，首版暂用9.17mm
HU103_LINER_ID_MM = HU103_LINER_OD_MM - 2.0 * HU103_LINER_WALL_THICKNESS_MM  # 完井尾管内径
HU103_BIT_DIAMETER_MM = 215.9    # 下段尾管井眼名义尺寸
HU103_STANDOFF_PROXY_PCT = 86.5  # 设计居中度代理值

# 呼103施工参数
HU103_PUMP_RATE_INJECTION_M3_MIN = 1.3       # 注前置液及水泥浆主排量
HU103_PUMP_RATE_DISP_HIGH_M3_MIN = 1.5       # 替浆第一阶段排量
HU103_PUMP_RATE_DISP_MED_M3_MIN = 1.0        # 替浆第二阶段排量
HU103_PUMP_RATE_DISP_LOW_M3_MIN = 0.8        # 替浆第三阶段排量
HU103_PUMP_RATE_DISP_FINAL_M3_MIN = 0.6      # 替浆第四阶段排量
HU103_BALANCE_VOLUME_M3 = 18.0               # 平衡液体积
HU103_FLUSH_VOLUME_M3 = 14.0                 # 驱油冲洗液体积
HU103_SPACER_VOLUME_M3 = 18.0                # 驱油隔离液体积
HU103_LEAD_VOLUME_M3 = 91.0                  # 领浆设计量，含附加10m³
HU103_INTERMEDIATE_VOLUME_M3 = 10.0          # 中间浆体积代理值，数据包未单独给出
HU103_TAIL_VOLUME_M3 = 19.0                  # 尾浆设计量，含下塞1.18m³
HU103_PLUG_VOLUME_M3 = 5.0                   # 压塞液体积
HU103_DISPLACEMENT_VOLUME_M3 = 110.0         # 顶替轻泥浆体积

# 呼103流变参数 — 钻井液/环空初始液
HU103_MUD_DENSITY_KG_M3 = 1980.0             # 固井时井液密度
HU103_MUD_PV_PA_S = 0.054                    # 54mPa·s 换算为 Pa·s
HU103_MUD_YP_PA = 12.5                       # 10min 凝胶值作为屈服值代理

# 呼103流变参数 — 平衡液/冲洗液/隔离液
HU103_BALANCE_DENSITY_KG_M3 = 1750.0         # 平衡液密度 1.75g/cm³
HU103_BALANCE_PV_PA_S = 0.025                # 缺实测，暂用冲洗液体系代理
HU103_BALANCE_YP_PA = 1.5                    # 缺实测，暂用冲洗液体系代理
HU103_FLUSH_DENSITY_KG_M3 = 1020.0           # 冲洗液密度 1.02g/cm³
HU103_FLUSH_PV_PA_S = 0.025                  # 冲洗液塑性粘度代理
HU103_FLUSH_YP_PA = 1.5                      # 冲洗液屈服值代理
HU103_SPACER_DENSITY_KG_M3 = 1800.0          # 隔离液设计密度 1.80g/cm³
HU103_SPACER_PV_PA_S = 0.035                 # 35mPa·s 换算为 Pa·s
HU103_SPACER_YP_PA = 8.0                     # 隔离液现场屈服值

# 呼103流变参数 — 三段水泥浆
HU103_LEAD_DENSITY_KG_M3 = 2050.0            # 领浆密度 2.05g/cm³
HU103_LEAD_N = 0.82                          # 领浆幂律流性指数
HU103_LEAD_K_PA_S_N = 0.67                   # 领浆幂律稠度系数
HU103_INTERMEDIATE_DENSITY_KG_M3 = 2050.0    # 中间浆密度 2.05g/cm³
HU103_INTERMEDIATE_N = 0.76                  # 中间浆幂律流性指数
HU103_INTERMEDIATE_K_PA_S_N = 1.11           # 中间浆幂律稠度系数
HU103_TAIL_DENSITY_KG_M3 = 2050.0            # 尾浆密度 2.05g/cm³
HU103_TAIL_N = 0.76                          # 尾浆幂律流性指数
HU103_TAIL_K_PA_S_N = 1.14                   # 尾浆幂律稠度系数

# 呼103流变参数 — 压塞液/替浆液
HU103_PLUG_DENSITY_KG_M3 = 1500.0            # 压塞液密度 1.50g/cm³
HU103_DISPLACEMENT_DENSITY_KG_M3 = 1500.0    # 替浆轻泥浆密度 1.50g/cm³
HU103_DISPLACEMENT_PV_PA_S = 0.054           # 替浆液塑性粘度代理，与固井时钻井液一致
HU103_DISPLACEMENT_YP_PA = 12.5              # 替浆液屈服值代理，与固井时钻井液一致


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
        raise ValueError(f"Hu103 井径/井斜 CSV 为空: {caliper_csv_path}")
    return tuple(sorted(rows))


def _standoff_value(depth_md_m: float, hole_diameter_mm: float, liner_od_mm: float) -> float:
    """计算居中度(standoff)剖面。

    基于呼103设计居中度代理值估算：
    - 数据包给出设计居中度代理值 86.5%
    - 首版加载器不引入扶正器明细，只做轻微分段修正

    居中度定义：standoff = 1 - 偏心度
    - standoff = 1.0: 完全居中
    - standoff = 0.0: 完全偏心
    """
    # 以设计居中度代理值为基准，底部尾管鞋附近略降，CBL评价主段保持设计值。
    standoff = HU103_STANDOFF_PROXY_PCT / 100.0
    if depth_md_m > 7712.0:
        standoff -= 0.02
    elif depth_md_m < 7338.0:
        standoff -= 0.01
    # 间隙修正因子：井眼明显扩径时略降低居中效果，避免把名义设计值外推过满。
    clearance_mm = max(hole_diameter_mm - liner_od_mm, 5.0)
    nominal_clearance_mm = max(HU103_BIT_DIAMETER_MM - liner_od_mm, 5.0)
    standoff *= min(max(nominal_clearance_mm / clearance_mm, 0.92), 1.03)
    return min(max(standoff, 0.30), 0.90)


def _build_hole_profile(profile_rows: Iterable[tuple[float, float, float]]) -> tuple[DepthValuePoint, ...]:
    return tuple(DepthValuePoint(depth_md_m=depth, value=hole_diameter) for depth, hole_diameter, _ in profile_rows)


def _build_inclination_profile(profile_rows: Iterable[tuple[float, float, float]]) -> tuple[DepthValuePoint, ...]:
    return tuple(DepthValuePoint(depth_md_m=depth, value=inclination) for depth, _, inclination in profile_rows)


def _build_standoff_profile(profile_rows: Iterable[tuple[float, float, float]]) -> tuple[DepthValuePoint, ...]:
    return tuple(
        DepthValuePoint(
            depth_md_m=depth,
            value=_standoff_value(depth, hole_diameter, HU103_LINER_OD_MM),
        )
        for depth, hole_diameter, _ in profile_rows
    )


def load_hu103_tailpipe(
    *,
    caliper_csv_path: Path | None = None,
    reference_root: Path | None = None,
) -> tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载呼103完井尾管段首版模型输入。

    Args:
        caliper_csv_path: 可选井径/井斜 CSV 路径。
        reference_root: 可选参考资料根目录。

    Returns:
        井筒参数、流体参数、泵注程序与验证资料路径。
    """

    resolved_reference_root = reference_root or DEFAULT_REFERENCE_ROOT
    resolved_caliper_csv_path = caliper_csv_path or DEFAULT_CALIPER_CSV

    profile_rows = _read_profile_rows(resolved_caliper_csv_path)
    well_spec = WellSpec(
        well_name="呼103",
        top_md_m=HU103_TOP_MD_M,
        bottom_md_m=HU103_BOTTOM_MD_M,
        shoe_md_m=HU103_SHOE_MD_M,
        hanger_md_m=HU103_HANGER_MD_M,
        casing_id_mm=HU103_CASING_ID_MM,
        liner_od_mm=HU103_LINER_OD_MM,
        liner_id_mm=HU103_LINER_ID_MM,
        hole_diameter_profile=_build_hole_profile(profile_rows),
        inclination_profile=_build_inclination_profile(profile_rows),
        standoff_profile=_build_standoff_profile(profile_rows),
        evaluation_windows=(
            EvaluationWindow(name="CBL评价井段", top_md_m=7338.0, bottom_md_m=7712.0, window_type="cbl"),
            EvaluationWindow(name="全井段", top_md_m=7330.6, bottom_md_m=7770.0, window_type="custom"),
        ),
        reference_root=resolved_reference_root,
        notes=(
            "仅对应呼103井 139.70mm 完井尾管段固井，不含回接固井与其他套管段。",
            "上段 φ168.30mm 尾管内径 150.42mm 为代理值，按壁厚8.94mm估算。",
            "φ139.70mm 完井尾管缺实测壁厚，首版暂用9.17mm计算内径。",
            "居中度剖面以设计居中度 86.5% 为代理，并按井段位置和井径间隙做轻微修正。",
        ),
    )

    fluids = (
        FluidSpec(
            name="钻井液",
            role=FluidRole.MUD,
            density_kg_m3=HU103_MUD_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU103_MUD_PV_PA_S,
            yield_stress_pa=HU103_MUD_YP_PA,
        ),
        FluidSpec(
            name="平衡液",
            role=FluidRole.SPACER,
            density_kg_m3=HU103_BALANCE_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU103_BALANCE_PV_PA_S,
            yield_stress_pa=HU103_BALANCE_YP_PA,
        ),
        FluidSpec(
            name="驱油冲洗液",
            role=FluidRole.WASH,
            density_kg_m3=HU103_FLUSH_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU103_FLUSH_PV_PA_S,
            yield_stress_pa=HU103_FLUSH_YP_PA,
        ),
        FluidSpec(
            name="隔离液",
            role=FluidRole.SPACER,
            density_kg_m3=HU103_SPACER_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU103_SPACER_PV_PA_S,
            yield_stress_pa=HU103_SPACER_YP_PA,
        ),
        FluidSpec(
            name="领浆",
            role=FluidRole.LEAD,
            density_kg_m3=HU103_LEAD_DENSITY_KG_M3,
            rheology_model=RheologyModel.POWER_LAW,
            power_law_n=HU103_LEAD_N,
            consistency_k=HU103_LEAD_K_PA_S_N,
        ),
        FluidSpec(
            name="中间浆",
            role=FluidRole.INTERMEDIATE,
            density_kg_m3=HU103_INTERMEDIATE_DENSITY_KG_M3,
            rheology_model=RheologyModel.POWER_LAW,
            power_law_n=HU103_INTERMEDIATE_N,
            consistency_k=HU103_INTERMEDIATE_K_PA_S_N,
        ),
        FluidSpec(
            name="尾浆",
            role=FluidRole.TAIL,
            density_kg_m3=HU103_TAIL_DENSITY_KG_M3,
            rheology_model=RheologyModel.POWER_LAW,
            power_law_n=HU103_TAIL_N,
            consistency_k=HU103_TAIL_K_PA_S_N,
        ),
        FluidSpec(
            name="压塞液",
            role=FluidRole.OTHER,
            density_kg_m3=HU103_PLUG_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU103_DISPLACEMENT_PV_PA_S,
            yield_stress_pa=HU103_DISPLACEMENT_YP_PA,
        ),
        FluidSpec(
            name="替浆液",
            role=FluidRole.DISPLACEMENT,
            density_kg_m3=HU103_DISPLACEMENT_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU103_DISPLACEMENT_PV_PA_S,
            yield_stress_pa=HU103_DISPLACEMENT_YP_PA,
        ),
    )

    schedule = PumpingSchedule(
        steps=(
            PumpingScheduleStep(
                step_name="注入平衡液",
                fluid_name="平衡液",
                volume_m3=HU103_BALANCE_VOLUME_M3,
                rate_m3_min=HU103_PUMP_RATE_INJECTION_M3_MIN,
                remarks="平衡液 18m3，密度1.75g/cm3；角色按 SPACER 处理，在二维模型中归入隔离液相。",
            ),
            PumpingScheduleStep(
                step_name="注入驱油冲洗液",
                fluid_name="驱油冲洗液",
                volume_m3=HU103_FLUSH_VOLUME_M3,
                rate_m3_min=HU103_PUMP_RATE_INJECTION_M3_MIN,
                remarks="驱油冲洗液 14m3，密度1.02g/cm3；角色为 WASH，在二维模型中归入隔离液相。",
            ),
            PumpingScheduleStep(
                step_name="注入驱油隔离液",
                fluid_name="隔离液",
                volume_m3=HU103_SPACER_VOLUME_M3,
                rate_m3_min=HU103_PUMP_RATE_INJECTION_M3_MIN,
                remarks="驱油隔离液 18m3，密度1.80g/cm3。",
            ),
            PumpingScheduleStep(
                step_name="注入领浆",
                fluid_name="领浆",
                volume_m3=HU103_LEAD_VOLUME_M3,
                rate_m3_min=HU103_PUMP_RATE_INJECTION_M3_MIN,
                remarks="领浆设计量 91m3，含附加10m3。",
            ),
            PumpingScheduleStep(
                step_name="注入中间浆",
                fluid_name="中间浆",
                volume_m3=HU103_INTERMEDIATE_VOLUME_M3,
                rate_m3_min=HU103_PUMP_RATE_INJECTION_M3_MIN,
                remarks="中间浆体积为 10m3 代理值；数据包未单独给出中间浆体积。",
            ),
            PumpingScheduleStep(
                step_name="注入尾浆",
                fluid_name="尾浆",
                volume_m3=HU103_TAIL_VOLUME_M3,
                rate_m3_min=HU103_PUMP_RATE_INJECTION_M3_MIN,
                remarks="尾浆设计量 19m3，含下塞1.18m3。",
            ),
            PumpingScheduleStep(
                step_name="压塞液",
                fluid_name="压塞液",
                volume_m3=HU103_PLUG_VOLUME_M3,
                rate_m3_min=HU103_PUMP_RATE_INJECTION_M3_MIN,
                remarks="压塞液 5m3，密度1.50g/cm3。",
            ),
            PumpingScheduleStep(
                step_name="替浆（第一阶段高速）",
                fluid_name="替浆液",
                volume_m3=33.0,
                rate_m3_min=HU103_PUMP_RATE_DISP_HIGH_M3_MIN,
                remarks="替浆总量 110m3 的约30%代理分配，高速阶段。",
            ),
            PumpingScheduleStep(
                step_name="替浆（第二阶段中速）",
                fluid_name="替浆液",
                volume_m3=33.0,
                rate_m3_min=HU103_PUMP_RATE_DISP_MED_M3_MIN,
                remarks="替浆总量 110m3 的约30%代理分配，中速阶段。",
            ),
            PumpingScheduleStep(
                step_name="替浆（第三阶段低速）",
                fluid_name="替浆液",
                volume_m3=22.0,
                rate_m3_min=HU103_PUMP_RATE_DISP_LOW_M3_MIN,
                remarks="替浆总量 110m3 的约20%代理分配，低速阶段。",
            ),
            PumpingScheduleStep(
                step_name="替浆（第四阶段终速）",
                fluid_name="替浆液",
                volume_m3=22.0,
                rate_m3_min=HU103_PUMP_RATE_DISP_FINAL_M3_MIN,
                remarks="替浆总量 110m3 的约20%代理分配，终速阶段。",
            ),
        ),
        notes=(
            "按呼103井数据包第二版施工顺序构建：平衡液、冲洗液、隔离液、领浆、中间浆、尾浆、压塞液、四段替浆。",
            "中间浆体积为 10m3 首版代理值，因数据包未单独给出中间浆体积。",
            "替浆四段体积按 30%/30%/20%/20% 代理分配，总量保持 110m3。",
            "平衡液密度较高，角色按 SPACER 处理；冲洗液角色为 WASH，二者在二维三相模型中均归入隔离液相。",
        ),
    )

    validation_data = ValidationData(
        cbl_summary_path=resolved_reference_root / "呼103井_CBL评价报告.pdf",
        notes=(
            "呼103井数据包第二版给出 CBL 合格率 12.06%，评价井段 7338.0–7712.0m。",
            "cbl_summary_path 指向参考文档/呼103下的 CBL PDF；若原始文件名不同，后续仅需更新路径。",
            "流体物性来自 model_ready_candidate 与 recommended_inputs；代理值已在对应常量和施工步骤备注中标注。",
        ),
    )
    return well_spec, fluids, schedule, validation_data


def build_hu103_annulus_inlet_provider(
    schedule: PumpingSchedule,
    fluids: tuple[FluidSpec, ...],
    annulus_boundary_mode: str = "sustained_tail",
) -> Callable[[float], AnnulusInletState]:
    """为 Hu103 完井尾管段首版模型构建环空入口边界提供器。"""

    warnings.warn(
        "build_hu103_annulus_inlet_provider 已废弃，请改用 CasingFlowSolver + build_coupled_annulus_inlet_provider",
        DeprecationWarning,
        stacklevel=2,
    )
    fluid_role_by_name: dict[str, FluidRole] = {fluid.name: fluid.role for fluid in fluids}
    cement_fluid_names = {
        fluid.name for fluid in fluids if fluid.role in {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}
    }
    if not cement_fluid_names:
        raise ValueError("至少需要一个角色为 LEAD、INTERMEDIATE 或 TAIL 的水泥浆流体")
    default_cement_name = "尾浆" if "尾浆" in cement_fluid_names else next(iter(sorted(cement_fluid_names)))

    def _phase_fractions_for_fluid(fluid_name: str) -> tuple[tuple[str, float], ...]:
        role = fluid_role_by_name.get(fluid_name, FluidRole.MUD)
        # 三相映射规则：领浆/中间浆/尾浆归为水泥相，冲洗液/隔离液归为隔离液相，其余归为泥浆相。
        if role in {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}:
            return (("cement", 1.0),)
        if role in {FluidRole.WASH, FluidRole.SPACER}:
            return (("spacer", 1.0),)
        return (("mud", 1.0),)

    steps = schedule.steps
    if len(steps) < 2:
        raise ValueError("Hu103 首版边界提供器至少需要水泥浆步骤和替浆步骤")
    # 多流体程序下替浆步骤不再固定为单一步骤，因此按流体角色定位第一个替浆阶段。
    displacement_step_index = next(
        (
            index
            for index, step in enumerate(steps)
            if fluid_role_by_name[step.fluid_name] == FluidRole.DISPLACEMENT
        ),
        None,
    )
    if displacement_step_index is None:
        raise ValueError("Hu103 边界提供器需要至少一个替浆液步骤")
    def _provider(time_s: float) -> AnnulusInletState:
        elapsed_s = 0.0
        for index, step in enumerate(steps):
            duration_s = 0.0 if step.rate_m3_min <= 0.0 else step.volume_m3 / step.rate_m3_min * 60.0
            if time_s < elapsed_s + duration_s - 1e-12:
                if index >= displacement_step_index and fluid_role_by_name[step.fluid_name] == FluidRole.DISPLACEMENT:
                    if annulus_boundary_mode == "sustained_tail":
                        return AnnulusInletState(
                            time_s=time_s,
                            flow_rate_m3_s=step.rate_m3_min / 60.0,
                            stage_name="替浆推进（尾浆等效入环空）",
                            phase_fractions=_phase_fractions_for_fluid(default_cement_name),
                        )
                    if annulus_boundary_mode == "volume_limited":
                        return AnnulusInletState(
                            time_s=time_s,
                            flow_rate_m3_s=0.0,
                            stage_name="替浆期管内推进（环空入口保持尾浆）",
                            phase_fractions=_phase_fractions_for_fluid(default_cement_name),
                        )
                    if annulus_boundary_mode == "tail_then_mud":
                        return AnnulusInletState(
                            time_s=time_s,
                            flow_rate_m3_s=step.rate_m3_min / 60.0,
                            stage_name="替浆钻井液入环空",
                            phase_fractions=_phase_fractions_for_fluid(step.fluid_name),
                        )
                    raise ValueError(f"不支持的环空边界模式: {annulus_boundary_mode}")

                return AnnulusInletState(
                    time_s=time_s,
                    flow_rate_m3_s=step.rate_m3_min / 60.0,
                    stage_name=step.step_name,
                    phase_fractions=_phase_fractions_for_fluid(step.fluid_name),
                )
            elapsed_s += duration_s

        if annulus_boundary_mode in {"sustained_tail", "volume_limited"}:
            return AnnulusInletState(
                time_s=time_s,
                flow_rate_m3_s=0.0,
                stage_name="替浆结束后保持",
                phase_fractions=_phase_fractions_for_fluid(default_cement_name),
            )
        return AnnulusInletState(
            time_s=time_s,
            flow_rate_m3_s=0.0,
            stage_name="替浆结束后保持",
            phase_fractions=(("mud", 1.0),),
        )

    return _provider


def export_hu103_sync_card_markdown(output_dir: Path) -> Path:
    """导出 呼103 同步画像卡为 Markdown 文件。"""

    well_spec, fluids, schedule, _ = load_hu103_tailpipe()
    solver = CasingFlowSolver()
    result = solver.run(well_spec, fluids, schedule)
    provenance = WELL_PROVENANCE.get("呼103")
    if provenance is None:
        from cemdisp.data.provenance import _WELL_NAME_ALIASES
        canonical = _WELL_NAME_ALIASES.get("呼103", "呼103")
        provenance = WELL_PROVENANCE[canonical]
    sync_card = build_sync_card("呼103", result.shoe_timeline, provenance)

    output_path = output_dir / ("呼103_同步画像卡.md")
    lines = [
        "# 呼103 同步画像卡",
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
