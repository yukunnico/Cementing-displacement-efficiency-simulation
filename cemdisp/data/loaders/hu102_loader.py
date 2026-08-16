"""
呼102尾管段标准数据加载器

本模块实现呼102井139.70mm尾管段固井的标准数据加载功能。

主要功能：
- 从现场提取包井径/井斜CSV读取实测剖面（7120–7735m 61 点实测；顶段双层套管按等效井径）
- 构建井筒几何参数（井段范围、套管尺寸、偏心度等）
- 定义钻井液、替浆液、领浆、尾管水泥浆、平衡液、隔离液、压塞液、后置液等物性参数
- 构建 2026-08 重建的完整现场泵注程序（前置液→试隔离液→隔离液→领浆→尾浆→压塞液→后置液→替浆→循环排混浆）
- 提供 legacy 环空入口边界状态提供器（仅用于旧模型对比）

物理参数说明（2026-08 已按现场提取包更新，旧代理值以 LEGACY 注释保留）：
- 井段范围: 6823.10m - 7735.00m（6823.10–7119.80m 为双层套管段，按等效井径处理）
- 尾管尺寸: 139.70mm OD, 108.10mm ID（底部15.8mm壁厚段；上部14.27mm壁厚段ID 111.16mm）
- 领浆: 10m³, 施工密度2.10g/cm³, 幂律 n=0.737, K=0.947（现场实测 20234.doc）
- 尾浆: 7m³, 施工密度2.10g/cm³, 幂律 n=0.737, K=0.947（现场实测，与领浆同体系）
- 替浆液: 74m³(汇总口径;泵冲记录72m³), 密度2.02g/cm³, 平均排量≈0.84m³/min（现场88min）
- 钻井液（环空初始液）: 密度2.02g/cm³, Bingham PV=66mPa·s, YP=10Pa（现场实测 20211.doc s1.2）

现场记录来源（2026-08 更新；主作业 2022-11-21，20211.doc 施工小结 + 10042.xlsx 日报）：
- 注水泥35.00t，现场拆分为领浆10m³ + 尾浆7m³（均 2.10g/cm³）；替浆液 2.02g/cm³、74m³
- 泵注时间：2022-11-21 17:00–21:00（含停泵/压塞/下胶塞等工序）
- 施工存在下部单流阀失效导致的方案临场调整：尾管内全部替入2.10g/cm³隔离液、不留上塞

LEGACY(2026-08 前)可选补充流体（邻井代理，仅 include_wash_spacer=True 时注入合成 FLUSHER）：
- 冲洗液(WASH): ρ=2050, PV=0.035, YP=8（呼探1-002/呼103邻井代理）
- 合成 FLUSHER: ρ=1880, PV=0.025, YP=1.5（呼103邻井代理，验证序列可表达）

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
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep, PumpingStageEvent
from cemdisp.data.validation_data import ValidationData
from cemdisp.data.well_spec import DepthValuePoint, EvaluationWindow, WellSpec
from cemdisp.models2d.boundary_bridge import AnnulusInletState
from cemdisp.data.provenance import WELL_PROVENANCE
from cemdisp.models2d.boundary_bridge import build_sync_card
from cemdisp.transport1d.casing_flow import CasingFlowSolver



PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_REFERENCE_ROOT = PROJECT_ROOT / "参考文档" / "呼102"
# 现场提取包 caliper_profile.csv：20211.doc §1.6 电测井径+井斜，7120–7735m，10m 间隔，61 点实测。
# LEGACY(2026-08 前): PROJECT_ROOT / "hu102model" / "hu102_tail_caliper_inclination.csv"（固定 215.9mm 合成剖面）。
DEFAULT_CALIPER_CSV = PROJECT_ROOT / "参考文档" / "现场资料提取" / "hu102_呼102" / "caliper_profile.csv"

# 呼102尾管段井段参数（单位：m, mm）—— 与核对汇总 §3.2 确认一致
HU102_TOP_MD_M = 6823.10        # 井段顶部测深（尾管顶/悬挂器顶）
HU102_BOTTOM_MD_M = 7735.00     # 井段底部测深
HU102_SHOE_MD_M = 7735.00       # 套管鞋深度
HU102_HANGER_MD_M = 6823.10     # 悬挂器位置深度（FHDS550 坐挂顶 6821.895，坐挂位 7665.516）
HU102_CASING_ID_MM = 219.10     # 上层技术套管外径（悬挂器以上；well_geometry 口径）
HU102_LINER_OD_MM = 139.70      # 尾管外径
HU102_LINER_WALL_THICKNESS_MM = 15.80  # 尾管壁厚（底部厚壁段 15.8mm，对应 ID 108.10mm）
HU102_LINER_ID_MM = HU102_LINER_OD_MM - 2.0 * HU102_LINER_WALL_THICKNESS_MM  # 尾管内径（108.10mm）
# 6823.10–7119.80m 双层套管段（139.7mm 尾管在 219.1mm 技术套管内）等效井径。
# 现场该段无裸眼井径测点（CBL 6840–7119.8m 标注"双层套管不评价"），保留旧 loader 等效值并显式标注。
HU102_DOUBLE_CASING_EQUIVALENT_DIAMETER_MM = 215.9   # LEGACY(2026-08 前)等效值；非实测

# 呼102施工参数（2026-08 按现场提取包 pumping_schedule.csv 重建）
HU102_CEMENT_MASS_T = 35.0              # 水泥浆总质量（LEGACY 汇总口径，35t≈领浆10+尾浆7=17m³）
HU102_CEMENT_DENSITY_KG_M3 = 2100.0     # 水泥浆施工密度（现场 20211.doc，领浆/尾浆均 2.10）
HU102_LEAD_DENSITY_KG_M3 = 2100.0       # 领浆施工密度 2.10（现场）
HU102_TAIL_DENSITY_KG_M3 = 2100.0       # 尾浆施工密度 2.10（现场）
HU102_LEAD_VOLUME_M3 = 10.0             # 领浆体积（现场，封固 6720–7300m）
HU102_TAIL_VOLUME_M3 = 7.0              # 尾浆体积（现场，封固 7300–7735m）
HU102_DISPLACEMENT_VOLUME_M3 = 74.0     # 替浆体积（汇总口径；现场泵冲记录到量 72m³）
HU102_DISPLACEMENT_DENSITY_KG_M3 = 2020.0  # 替浆钻井液密度（现场 2.02，油基）
HU102_BALANCE_VOLUME_M3 = 15.0          # 前置液（平衡液）体积（现场）
HU102_BALANCE_DENSITY_KG_M3 = 1900.0    # 平衡液密度（现场 1.90）
HU102_SPACER_TEST_VOLUME_M3 = 5.0       # 试隔离液体积（现场）
HU102_SPACER_VOLUME_M3 = 15.0           # 隔离液体积（现场 15m³，密度 2.05）
HU102_PLUG_VOLUME_M3 = 7.0              # 压塞液体积（现场 7m³，密度 2.10）
HU102_PLUG_DENSITY_KG_M3 = 2100.0       # 压塞液密度（现场 2.10）
HU102_AFTERFLUID_VOLUME_M3 = 6.0        # 后置液（保护液）体积（现场 6m³，密度 1.90）
HU102_AFTERFLUID_DENSITY_KG_M3 = 1900.0 # 后置液（保护液）密度（现场 1.90）
HU102_CIRCULATION_VOLUME_M3 = 41.0      # 循环排混浆体积（现场 41m³，次日后处理）

# 泵注排量（2026-08 重建；每步取现场记录排量段的代表性值，取值口径见对应步骤备注）
HU102_RATE_BALANCE_M3_MIN = 0.45        # 平衡液 0.3–0.6 → 取 0.45
HU102_RATE_SPACER_TEST_M3_MIN = 0.6     # 试隔离液 0.6
HU102_RATE_SPACER_M3_MIN = 0.7          # 隔离液 0.6–0.8 → 取 0.7
HU102_RATE_CEMENT_M3_MIN = 0.7          # 领浆/尾浆 0.6–0.8 → 取 0.7
HU102_RATE_PLUG_M3_MIN = 0.6            # 压塞液 0.3–0.8 → 取 0.6
HU102_RATE_AFTERFLUID_M3_MIN = 0.7      # 后置液 0.6–0.8 → 取 0.7
# 替浆：现场 18:57–20:25 共 88min，排量逐步 0.9→0.46；平均 74/88 ≈ 0.841 m³/min（源自现场时长）。
HU102_RATE_DISPLACEMENT_M3_MIN = 74.0 / 88.0
HU102_RATE_CIRCULATION_M3_MIN = 1.3     # 循环排混浆（现场 1.3）
# LEGACY(2026-08 前)：按 17:00–21:00 累计 90.67m³/240min 校正的旧平均排量，已废弃。
HU102_RATE_M3_MIN = 0.378

# 呼102流变参数（2026-08 现场实测替换，旧代理值保留注释）
HU102_MUD_PV_PA_S = 0.066            # 环空初始钻井液塑性粘度 PV=66mPa·s（现场实测 20211.doc s1.2）；LEGACY: 0.080
HU102_MUD_YP_PA = 10.0               # 环空初始钻井液屈服值 YP=10Pa（现场实测）；LEGACY: 15.0
HU102_DISPLACEMENT_PV_PA_S = 0.066   # 替浆液塑性粘度（现场替浆用油基钻井液，与环空初始液一致）；LEGACY: 0.080
HU102_DISPLACEMENT_YP_PA = 10.0      # 替浆液屈服值（与钻井液一致）；LEGACY: 15.0
HU102_CEMENT_POWER_LAW_N = 0.737     # 尾浆流性指数 n（现场实测 20234.doc 93C）；LEGACY: 0.722
HU102_CEMENT_CONSISTENCY_K = 0.947   # 尾浆稠度系数 K Pa·s^n（现场实测）；LEGACY: 0.684
HU102_LEAD_POWER_LAW_N = 0.737       # 领浆流性指数 n（现场实测，与尾浆同体系 20234.doc）
HU102_LEAD_CONSISTENCY_K = 0.947     # 领浆稠度系数 K Pa·s^n（现场实测）

# 呼102前置液/隔离液/压塞液/后置液参数
# 现场 20211.doc 只给出各前置流体密度，无流变实测；流变沿用旧 loader 邻井代理值并标注 proxy。
# 平衡液（前置液）：密度 1.90（现场）；流变 PV/YP 为代理（LEGACY: 呼探1-002/呼103 邻井 WASH 体系）。
HU102_BALANCE_PV_PA_S = 0.025
HU102_BALANCE_YP_PA = 1.5
# 隔离液：密度 2.05（现场 20211.doc）；流变 PV/YP 为代理（呼探1-002 邻井）。
HU102_SPACER_DENSITY_KG_M3 = 2050.0
HU102_SPACER_PV_PA_S = 0.035
HU102_SPACER_YP_PA = 8.0
# 压塞液：密度 2.10（现场）；流变复用钻井液（PV=0.066/YP=10）。
# 后置液（保护液）：密度 1.90（现场）；流变复用钻井液（proxy）。
# 压塞液/后置液流变参数（proxy，复用钻井液）
HU102_PLUG_PV_PA_S = 0.066
HU102_PLUG_YP_PA = 10.0
HU102_AFTERFLUID_PV_PA_S = 0.066
HU102_AFTERFLUID_YP_PA = 10.0

# LEGACY(2026-08 前) 邻井代理流体参数：仅保留用于 include_wash_spacer=True 的敏感性 FLUSHER 步骤，
# 以及"冲洗液"（WASH）流体定义（测试兼容）。不再作为呼102现场实录。
HU102_WASH_DENSITY_KG_M3 = 2050.0    # LEGACY 平衡液/冲洗液密度（呼探1-002邻井代理）
HU102_WASH_PV_PA_S = 0.035           # LEGACY
HU102_WASH_YP_PA = 8.0               # LEGACY
HU102_WASH_VOLUME_M3 = 10.0          # LEGACY 平衡液/冲洗液设计体积（呼探1-002邻井代理）
HU102_FLUSHER_DENSITY_KG_M3 = 1880.0 # LEGACY 合成 FLUSHER（呼103邻井代理 ρ=1880, PV=0.025, YP=1.5）
HU102_FLUSHER_PV_PA_S = 0.025
HU102_FLUSHER_YP_PA = 1.5
HU102_FLUSHER_VOLUME_M3 = 5.0


def _read_profile_rows(caliper_csv_path: Path) -> tuple[tuple[float, float, float], ...]:
    """读取井径/井斜剖面 CSV。

    兼容两套列名（2026-08 起现场提取包 caliper_profile.csv 为默认输入）：
    - 现场提取包: md_m / caliper_mm / inclination_deg
    - 旧模型 CSV: depth_md_m / hole_diameter_mm / inclination_deg
    """
    with caliper_csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows: list[tuple[float, float, float]] = []
        for row in csv.DictReader(handle):
            depth = row.get("depth_md_m", row.get("md_m"))
            hole = row.get("hole_diameter_mm", row.get("caliper_mm"))
            inclination = row.get("inclination_deg")
            if depth is None or hole is None or inclination is None:
                continue
            rows.append((float(depth), float(hole), float(inclination)))
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

    measured_rows = _read_profile_rows(resolved_caliper_csv_path)
    # 6823.10–7119.80m 双层套管段无裸眼井径测点，前插一个等效顶段点；
    # 井斜取 7120m 首测点值 4.82° 作为顶段代理（连续测斜表缺，低置信）。
    top_rows = (
        (HU102_TOP_MD_M, HU102_DOUBLE_CASING_EQUIVALENT_DIAMETER_MM, 4.82),
    )
    profile_rows = top_rows + measured_rows
    well_spec = WellSpec(
        well_name="呼102",
        top_md_m=HU102_TOP_MD_M,
        bottom_md_m=HU102_BOTTOM_MD_M,
        shoe_md_m=HU102_SHOE_MD_M,
        hanger_md_m=HU102_HANGER_MD_M,
        casing_id_mm=HU102_CASING_ID_MM,
        liner_od_mm=HU102_LINER_OD_MM,
        liner_id_mm=HU102_LINER_ID_MM,
        liner_wall_thickness_mm=HU102_LINER_WALL_THICKNESS_MM,
        hole_diameter_profile=_build_hole_profile(profile_rows),
        inclination_profile=_build_inclination_profile(profile_rows),
        standoff_profile=_build_standoff_profile(profile_rows),
        evaluation_windows=(
            # 拆分 CBL 质量窗与地层目标窗（核对汇总 §10.2）：
            # CBL 评价窗/质量窗来自 100413.PDF 分段解释；地层目标窗来自 20212.doc 基础信息。
            EvaluationWindow(name="CBL评价井段", top_md_m=6840.0, bottom_md_m=7665.0, window_type="cbl"),
            EvaluationWindow(name="CBL质量段一(胶结差/不合格)", top_md_m=7405.0, bottom_md_m=7480.0, window_type="cbl_quality"),
            EvaluationWindow(name="CBL质量段二(胶结差/不合格)", top_md_m=7502.0, bottom_md_m=7540.0, window_type="cbl_quality"),
            EvaluationWindow(name="J3k目的层(地层目标)", top_md_m=7422.0, bottom_md_m=7640.0, window_type="formation_target"),
            EvaluationWindow(name="油气显示层1", top_md_m=7432.914, bottom_md_m=7434.268, window_type="oil_gas_show"),
            EvaluationWindow(name="油气显示层2", top_md_m=7638.05, bottom_md_m=7639.403, window_type="oil_gas_show"),
        ),
        reference_root=resolved_reference_root,
        notes=(
            "仅对应呼102井 139.70mm 尾管段固井，不含回接固井与其他套管段。",
            "7120–7735m 井径/井斜为 61 点现场实测（20211.doc §1.6，10m 间隔，184.85–200.58mm）；6823.10–7119.80m 双层套管段按等效井径 215.9mm 处理（等效，非实测）。",
            "最大井斜 11.49°：现场摘要一处记@7450m（20211.doc 井径表）、一处记@7490m（测斜表），两处并存、峰值一致，模型按逐点实测剖面取值。",
        ),
    )

    # 2026-08 重建：领浆/尾浆独立建模；前置液/隔离液/压塞液/后置液按现场密度重建。
    # 流变标注：钻井液/水泥浆为现场实测（20211.doc / 20234.doc）；平衡液、隔离液仅密度实测，
    # 流变沿用邻井代理；压塞液、后置液仅密度实测，流变复用钻井液（proxy，非实测）。
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
            name="平衡液",
            role=FluidRole.WASH,
            density_kg_m3=HU102_BALANCE_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU102_BALANCE_PV_PA_S,
            yield_stress_pa=HU102_BALANCE_YP_PA,
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
            name="领浆",
            role=FluidRole.LEAD,
            density_kg_m3=HU102_LEAD_DENSITY_KG_M3,
            rheology_model=RheologyModel.POWER_LAW,
            power_law_n=HU102_LEAD_POWER_LAW_N,
            consistency_k=HU102_LEAD_CONSISTENCY_K,
        ),
        FluidSpec(
            name="尾管水泥浆",
            role=FluidRole.TAIL,
            density_kg_m3=HU102_TAIL_DENSITY_KG_M3,
            rheology_model=RheologyModel.POWER_LAW,
            power_law_n=HU102_CEMENT_POWER_LAW_N,
            consistency_k=HU102_CEMENT_CONSISTENCY_K,
        ),
        FluidSpec(
            name="压塞液",
            role=FluidRole.OTHER,
            density_kg_m3=HU102_PLUG_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU102_PLUG_PV_PA_S,
            yield_stress_pa=HU102_PLUG_YP_PA,
        ),
        FluidSpec(
            name="后置液",
            role=FluidRole.OTHER,
            density_kg_m3=HU102_AFTERFLUID_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU102_AFTERFLUID_PV_PA_S,
            yield_stress_pa=HU102_AFTERFLUID_YP_PA,
        ),
        # LEGACY 邻井代理流体：保留定义以维持 mud-spacer-flusher-cement 序列可表达及测试兼容，
        # 不作为呼102现场实录（include_wash_spacer=False 默认不进 schedule）。
        FluidSpec(
            name="冲洗液",
            role=FluidRole.WASH,
            density_kg_m3=HU102_WASH_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU102_WASH_PV_PA_S,
            yield_stress_pa=HU102_WASH_YP_PA,
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

    # 2026-08 重建完整现场泵注程序（20211.doc 施工记录泵注部分）：
    # 前置液(平衡液)→试隔离液→隔离液→领浆→尾浆→压塞液→后置液(保护液)→替浆→循环排混浆。
    # 说明：施工存在下部单流阀失效导致的方案临场调整（尾管内全部替入2.10隔离液、不留上塞）；
    # 循环排混浆为次日后处理（00:30），与水泥顶替过程解耦。
    field_schedule_steps = (
        PumpingScheduleStep(
            step_name="注入平衡液（前置液）",
            fluid_name="平衡液",
            volume_m3=HU102_BALANCE_VOLUME_M3,
            rate_m3_min=HU102_RATE_BALANCE_M3_MIN,
            remarks="现场20211.doc: 前置液15m³，密度1.90，排量0.3-0.6（取0.45）。",
        ),
        PumpingScheduleStep(
            step_name="注入试隔离液",
            fluid_name="隔离液",
            volume_m3=HU102_SPACER_TEST_VOLUME_M3,
            rate_m3_min=HU102_RATE_SPACER_TEST_M3_MIN,
            remarks="现场: 试隔离液5m³，密度2.05，排量0.6，泵压15.5MPa。",
        ),
        PumpingScheduleStep(
            step_name="注入隔离液",
            fluid_name="隔离液",
            volume_m3=HU102_SPACER_VOLUME_M3,
            rate_m3_min=HU102_RATE_SPACER_M3_MIN,
            remarks="现场: 隔离液15m³，密度2.05，排量0.6-0.8（取0.7），泵压13-15.7MPa。",
        ),
        PumpingScheduleStep(
            step_name="注入领浆",
            fluid_name="领浆",
            volume_m3=HU102_LEAD_VOLUME_M3,
            rate_m3_min=HU102_RATE_CEMENT_M3_MIN,
            remarks="现场: 领浆10m³，密度2.10，排量0.6-0.8（取0.7），封固6720-7300m。",
        ),
        PumpingScheduleStep(
            step_name="注入尾浆",
            fluid_name="尾管水泥浆",
            volume_m3=HU102_TAIL_VOLUME_M3,
            rate_m3_min=HU102_RATE_CEMENT_M3_MIN,
            remarks="现场: 尾浆7m³，密度2.10，排量0.6-0.8（取0.7），封固7300-7735m。",
        ),
        PumpingScheduleStep(
            step_name="注入压塞液",
            fluid_name="压塞液",
            volume_m3=HU102_PLUG_VOLUME_M3,
            rate_m3_min=HU102_RATE_PLUG_M3_MIN,
            remarks="现场: 压塞液7m³，密度2.10，排量0.3-0.8（取0.6），泵压6-13MPa。",
        ),
        PumpingScheduleStep(
            step_name="注入后置液（保护液）",
            fluid_name="后置液",
            volume_m3=HU102_AFTERFLUID_VOLUME_M3,
            rate_m3_min=HU102_RATE_AFTERFLUID_M3_MIN,
            remarks="现场: 保护液6m³，密度1.90，排量0.6-0.8（取0.7）。",
        ),
        PumpingScheduleStep(
            step_name="替浆（钻井液）",
            fluid_name="替浆液",
            volume_m3=HU102_DISPLACEMENT_VOLUME_M3,
            rate_m3_min=HU102_RATE_DISPLACEMENT_M3_MIN,
            remarks="现场: 替浆74m³(汇总口径;泵冲到量72m³)，密度2.02，排量0.9→0.46，平均≈74/88min=0.841。",
        ),
        PumpingScheduleStep(
            step_name="循环排混浆",
            fluid_name="钻井液",
            volume_m3=HU102_CIRCULATION_VOLUME_M3,
            rate_m3_min=HU102_RATE_CIRCULATION_M3_MIN,
            event_tag=PumpingStageEvent.RESTART,
            remarks="现场(次日后处理): 敞压循环排混浆41m³，排量1.3；属候凝前清理，与水泥顶替过程解耦。",
        ),
    )

    # include_wash_spacer=True 时插入合成 FLUSHER 步骤（邻井代理敏感性输入，验证
    # mud-spacer-flusher-cement 序列可表达），位于隔离液之后、领浆之前。
    # 注：2026-08 起现场平衡液/隔离液已作为实录进入默认 schedule，该参数仅额外注入 FLUSHER。
    flusher_step = PumpingScheduleStep(
        step_name="注入冲洗液（FLUSHER）",
        fluid_name="冲洗液（FLUSHER）",
        volume_m3=HU102_FLUSHER_VOLUME_M3,
        rate_m3_min=HU102_RATE_M3_MIN,
        remarks=f"合成冲洗液（FLUSHER）{HU102_FLUSHER_VOLUME_M3}m³，密度{HU102_FLUSHER_DENSITY_KG_M3/1000:.2f}g/cm³（呼103邻井代理，验证序列可表达）。",
    )

    if include_wash_spacer:
        # 插入到"注入隔离液"（index 2）之后、"注入领浆"（index 3）之前
        steps = (
            field_schedule_steps[:3]
            + (flusher_step,)
            + field_schedule_steps[3:]
        )
    else:
        steps = field_schedule_steps

    schedule = PumpingSchedule(
        steps=steps,
        notes=(
            "2026-08 重建完整现场泵注程序（20211.doc 施工记录）：前置液→试隔离液→隔离液→领浆→尾浆→压塞液→后置液→替浆→循环排混浆。",
            "施工存在下部单流阀失效导致的方案临场调整（尾管内全部替入2.10隔离液、不留上塞）；替浆与循环之间有关井/拔中心管等非泵注工序。",
            "循环排混浆(41m³)为次日后处理，仅作现场记录，与水泥顶替过程解耦，可由消费方按需剔除。",
            "替换 LEGACY(2026-08 前) 两步简化程序（尾管水泥浆+替浆液推进）。",
            "include_wash_spacer=True 时额外注入合成 FLUSHER(5m³) 邻井代理敏感性步骤，验证 mud-spacer-flusher-cement 序列可表达。",
        ),
    )

    validation_data = ValidationData(
        cbl_summary_path=resolved_reference_root / "1004" / "10041" / "100413.PDF",
        cbl_pass_rate=0.6665,
        job_report_path=resolved_reference_root / "1004" / "10042.xlsx",
        pump_pressure_series_path=resolved_reference_root / "1004" / "100492.xlsx",
        notes=(
            "100413.PDF 给出 CBL 合格率 66.65%，评价井段 6840–7665m（与采用值对照表一致）。",
            "CBL 分段解释（100413.PDF）：7119.8–7405 良好-中等；7405–7480 胶结差(不合格)；7480–7502 中等；7502–7540 胶结差(不合格)；7540–7665 良好-中等。",
            "钻井液 PV=66mPa·s/YP=10Pa、领/尾浆 n=0.737/K=0.947 已更新为现场实测（2026-08）。",
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
