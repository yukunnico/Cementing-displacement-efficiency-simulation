"""
呼103尾管段标准数据加载器

本模块实现呼103井 139.70mm 完井尾管段固井的标准数据加载功能。

主要功能：
- 从井径/井斜 CSV 读取现场实测分段曲线（5755–7770m 电测井径/井斜；5536–5755m 重叠段等效外推）
- 构建复合尾管井筒几何（上段 168.3mm / 下段 139.7mm，dual-diameter 字段 + liner_od_profile）
- 定义钻井液、平衡液、隔离液1/2、领浆、中间浆、尾浆、压塞液、中置液与替浆液物性参数
- 构建设计版（load_hu103_tailpipe，20313.doc）与实际版（load_hu103_tailpipe_actual，20314.doc）两套泵注程序
- 提供 legacy 环空入口边界状态提供器（仅用于旧模型对比）

物理参数说明（2026-08 已按现场提取包更新，旧代理值以 LEGACY 注释保留）：
- 井段范围: 5536.662m（回接筒顶/尾管串顶） - 7770.00m（尾管鞋）
- 悬挂器: 5536.662–5545.972m（NSSX-CYFMD 封隔一体式，273.1x193.7 转 168.3mm）
- 复合尾管: 上段 168.30mm OD / 138.90mm ID（5546.022–7330.694m），下段 139.70mm OD / 107.94mm ID（7330.694–7768.51m）
- 井眼尺寸: 260.0mm(5750–7382 扩眼段) / 241.3mm(7382–7563) / 215.9mm(7563–7770)，井径/井斜用实测分段曲线
- 钻井液（环空初始液）: 密度1.98g/cm³, Bingham PV=54mPa·s, YP=10Pa（现场 65C, 20313.doc 1.3.1）
- 设计版（load_hu103_tailpipe 默认）: 平衡液28 + 隔离液17.5×2 + 领浆22 + 中间浆35 + 尾浆26 + 压塞液7 + 替浆83.2（20313.doc 5.2.4；2026-08-29 校准）
- 实际版（load_hu103_tailpipe_actual）: 平衡液48 + 隔离液35 + 水泥浆(三凝,分相按设计) + 压塞液3.0(ρ2.20) + 替浆90.9，全程 0.6m³/min（20314.doc 施工记录表）
- CBL: 139.7mm 段 7338–7712m 合格率 12.06%（结构化 cbl_pass_rate=0.1206）

现场资料来源（2026-08 更新，现场提取包 hu103_呼103/）：
- 井身结构/悬挂器/尾管尺寸：20313.doc（施工设计）、20315.doc（施工总结）、20318.doc（作业史）
- 井径/井斜：20313.doc §1.5.1/1.5.2 电测分段数据（5755–7770m）
- 泵注设计版：20313.doc 5.2.4；实际版：20314.doc 施工记录表
- 水泥浆流变：203111.docx 实验报告（领浆 n=0.82/K=0.67、中间浆 n=0.76/K=1.11、尾浆 n=0.76/K=1.14）；
  隔离液流变实测（203111 表7，140↘93℃）：n=0.54/K=2.12（2026-08-29 校准，旧 Bingham 代理 PV35/YP8 废弃）
- CBL：100516.pdf 139.7mm 段 7338.0–7712.0m 合格率 12.06%（结构化 cbl_pass_rate=0.1206）

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

# 呼103尾管段井段参数（单位：m, mm）—— 2026-08 按现场提取包更新
HU103_TOP_MD_M = 5536.662        # 井段顶部（实际回接筒顶 5536.66m，20318.doc）
HU103_BOTTOM_MD_M = 7770.0       # 完钻井深/尾管鞋
HU103_SHOE_MD_M = 7770.0         # 套管鞋深度
HU103_HANGER_MD_M = 5545.972     # 悬挂器底/坐挂点（悬挂器实际跨度 5536.662–5545.972m）
# LEGACY(2026-08 前): HU103_TOP_MD_M/HANGER 误用 7330.6（实为 139.7mm 段变扣位置）
HU103_LINER_OD_MM = 139.7        # 下段完井尾管外径（139.7mm 段 7330.694–7768.51m）
HU103_LINER_WALL_THICKNESS_MM = 15.88  # 139.7mm 段壁厚（BG140V 15.88mm，20313.doc 套管数据）；LEGACY: 9.17
HU103_LINER_ID_MM = 107.94       # 139.7mm 段内径（实测/设计口径 20313.doc）；LEGACY: 121.36
HU103_UPPER_LINER_OD_MM = 168.3  # 上段 168.3mm 尾管外径（5546.022–7330.694m，20315.doc 施工总结）
HU103_UPPER_LINER_ID_MM = 138.9  # 上段 168.3mm 尾管内径（BG140V 14.7mm 壁厚）；LEGACY: 150.42
HU103_UPPER_SECTION_BOTTOM_MD_M = 7330.694  # 上段（168.3mm）底/变扣位置
HU103_CASING_ID_MM = 138.9       # 上层尾管内径口径（= 168.3mm 段 ID；dual-diameter 上层用）
HU103_BIT_DIAMETER_MM = 215.9    # 下段 139.7mm 尾管井眼名义尺寸（7563–7770m）
HU103_STANDOFF_PROXY_PCT = 77.8  # 设计居中度代理值（20313.doc 6.3 软件模拟 77.8%，无实测）；LEGACY: 86.5

# 呼103施工参数 —— 设计版（20313.doc 5.2.4 注入/替浆浆柱结构）
HU103_DESIGN_RATE_FRONT_M3_MIN = 1.6        # 平衡液/隔离液设计排量 1.6m³/min
HU103_DESIGN_RATE_CEMENT_M3_MIN = 1.8       # 领浆/中间浆/尾浆/压塞液/替浆设计排量 1.8m³/min
HU103_BALANCE_VOLUME_M3 = 28.0              # 平衡液(轻泥浆)设计 28m³，ρ1.88（800m 占高）
HU103_SPACER1_VOLUME_M3 = 17.5              # 驱油隔离液1 设计 17.5m³，ρ2.00（500m）
HU103_SPACER2_VOLUME_M3 = 17.5              # 驱油隔离液2 设计 17.5m³，ρ1.95（500m）
HU103_LEAD_VOLUME_M3 = 22.0                 # 领浆设计 22m³（封固 5350–6000m）
HU103_INTERMEDIATE_VOLUME_M3 = 35.0         # 中间浆设计 35m³（封固 6000–7000m；0708 多数口径 4 处：7.1.1 用量表总量 83=22+35+26、7.2 流程、实际 35）；LEGACY(2026-08-29 前): 36（浆柱表孤值）；LEGACY 代理: 10
HU103_TAIL_VOLUME_M3 = 26.0                 # 尾浆设计 26m³（封固 7000–7770m）；LEGACY 代理: 19（采用值对照表记 19.0，与设计 26 存在口径差异，见 notes）
HU103_PLUG_VOLUME_M3 = 7.0                  # 压塞液设计 7m³，ρ2.10（管内占高 582m）
HU103_DISP_MUD1_VOLUME_M3 = 23.0            # 替钻井液 23m³（管内占高 1520m）
HU103_DISP_MIDDLE_VOLUME_M3 = 15.0          # 替中置液 15m³，ρ1.95（占高 1600m）
HU103_DISP_MUD2_VOLUME_M3 = 45.2  # 替钻井液 45.2m³（设计 1.8–0.8 变排量段）；LEGACY(2026-08-29 前): 46.2（排量表笔误）；20313 七.1.4 总量 90.2 自洽（23+15+45.2+压塞7=90.2）。
HU103_DISPLACEMENT_VOLUME_M3 = (
    HU103_DISP_MUD1_VOLUME_M3 + HU103_DISP_MIDDLE_VOLUME_M3 + HU103_DISP_MUD2_VOLUME_M3
)  # 设计替浆总量 83.2m³（23+15+45.2；+压塞液7=90.2 与 20313 七.1.4 自洽）；LEGACY(2026-08-29 前): 84.2（46.2 排量表笔误）。

# 呼103施工参数 —— 实际版（20314.doc 施工记录表；环空憋堵迫使降排量）
HU103_ACTUAL_RATE_M3_MIN = 0.6              # 实际全程排量 0.6m³/min（设计 1.8 的 1/3）
HU103_ACTUAL_BALANCE_VOLUME_M3 = 48.0       # 实际平衡液 48m³（泵压 8–12MPa）
HU103_ACTUAL_SPACER_VOLUME_M3 = 35.0        # 实际隔离液 35m³，ρ1.92（泵压 15MPa）
HU103_ACTUAL_DISPLACEMENT_VOLUME_M3 = 90.9  # 实际替浆 90.9m³（泵压 13–9MPa）
# 实际水泥浆为"三凝"体系且无分相体积记录；分相体积保守沿用设计版 22/36/26（见 load_hu103_tailpipe_actual）

# 呼103流变参数 — 钻井液/环空初始液（现场 65C，20313.doc 1.3.1）
HU103_MUD_DENSITY_KG_M3 = 1980.0             # 固井时井液密度（白油基）
HU103_MUD_PV_PA_S = 0.054                    # PV=54mPa·s 换算为 Pa·s
HU103_MUD_YP_PA = 10.0                       # YP=10Pa（现场 65C）；LEGACY: 12.5（原用终切凝胶值）

# 呼103流变参数 — 平衡液/隔离液（密度现场，流变无实测时用代理并标注）
HU103_BALANCE_DENSITY_KG_M3 = 1880.0         # 平衡液密度 1.88（现场实际 20314.doc）；LEGACY: 1750
HU103_BALANCE_PV_PA_S = 0.025                # 平衡液 PV=25mPa·s（现场实际）；LEGACY 代理口径
HU103_BALANCE_YP_PA = 1.5                    # 平衡液 YP=1.5Pa（现场实际）
HU103_SPACER1_DENSITY_KG_M3 = 2000.0         # 驱油隔离液1 设计密度 2.00（20313.doc 7.1.2）
HU103_SPACER2_DENSITY_KG_M3 = 1950.0         # 驱油隔离液2 设计密度 1.95
HU103_SPACER_ACTUAL_DENSITY_KG_M3 = 1920.0   # 实际隔离液密度 1.92（20314.doc）
# 隔离液流变（2026-08-29 校准）：203111 表7 实测（140↘93℃），六速 128/74/54/33/6/5 → POWER_LAW n=0.54/K=2.12。
# LEGACY(2026-08-29 前): Bingham 代理 PV35/YP8（无实测）。
HU103_SPACER_POWER_LAW_N = 0.54              # 隔离液流性指数 n（203111 表7 实测）
HU103_SPACER_CONSISTENCY_K = 2.12            # 隔离液稠度系数 K Pa·s^n（同上）
HU103_SPACER_PV_PA_S = 0.035                 # LEGACY(2026-08-29 前) 隔离液塑性粘度代理（无实测），已由幂律实测替换
HU103_SPACER_YP_PA = 8.0                     # LEGACY(2026-08-29 前) 隔离液屈服值代理（无实测），已由幂律实测替换
# LEGACY(2026-08 前): FLUSH(冲洗液) ρ1.02、SPACER ρ1.80 等代理已按现场设计/实际密度重建。
HU103_FLUSH_DENSITY_KG_M3 = 1020.0           # LEGACY 冲洗液密度（现场设计无独立冲洗液段，保留常量）
HU103_FLUSH_PV_PA_S = 0.025
HU103_FLUSH_YP_PA = 1.5

# 呼103流变参数 — 三段水泥浆（203111.docx Table7 实验报告，与核对汇总 §4.2 确认一致）
HU103_LEAD_DENSITY_KG_M3 = 2050.0            # 领浆密度 2.05g/cm³
HU103_LEAD_N = 0.82                          # 领浆幂律流性指数
HU103_LEAD_K_PA_S_N = 0.67                   # 领浆幂律稠度系数
HU103_INTERMEDIATE_DENSITY_KG_M3 = 2050.0    # 中间浆密度 2.05g/cm³
HU103_INTERMEDIATE_N = 0.76                  # 中间浆幂律流性指数
HU103_INTERMEDIATE_K_PA_S_N = 1.11           # 中间浆幂律稠度系数
HU103_TAIL_DENSITY_KG_M3 = 2050.0            # 尾浆密度 2.05g/cm³
HU103_TAIL_N = 0.76                          # 尾浆幂律流性指数
HU103_TAIL_K_PA_S_N = 1.14                   # 尾浆幂律稠度系数

# 呼103流变参数 — 压塞液/替浆液（2026-08 按现场密度重建；流变复用钻井液 proxy）
HU103_PLUG_DENSITY_KG_M3 = 2100.0            # 压塞液密度 2.10（设计 20313.doc 7.1.3）；LEGACY: 1500
HU103_MIDDLE_FLUID_DENSITY_KG_M3 = 1950.0    # 中置液密度 1.95（设计 7.1.3）；0708 分歧 1.95（5.2.4 表）/2.0（7.1.3/9.2/203111 表3 三处）；实际 1.90（20315），2026-08-29 校准扩注
HU103_DISPLACEMENT_DENSITY_KG_M3 = 1980.0    # 替浆钻井液密度 1.98（设计/实际均为钻井液）；LEGACY: 1500
HU103_DISPLACEMENT_PV_PA_S = 0.054           # 替浆液塑性粘度代理，与固井时钻井液一致
HU103_DISPLACEMENT_YP_PA = 10.0              # 替浆液屈服值代理，与钻井液一致（YP=10）


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


def _nominal_bit_diameter_mm(depth_md_m: float) -> float:
    """按深度返回名义井径（20313.doc 井身质量情况：5750–7382 扩眼 260 / 7382–7563 241.3 / 7563–7770 215.9）。"""
    if depth_md_m < 7382.0:
        return 260.0
    if depth_md_m < 7563.0:
        return 241.3
    return 215.9


def _liner_od_at(depth_md_m: float) -> float:
    """返回复合尾管在给定深度的外径（上段 168.3mm / 下段 139.7mm）。"""
    if depth_md_m < HU103_UPPER_SECTION_BOTTOM_MD_M:
        return HU103_UPPER_LINER_OD_MM
    return HU103_LINER_OD_MM


def _standoff_value(depth_md_m: float, hole_diameter_mm: float, liner_od_mm: float, nominal_bit_mm: float) -> float:
    """计算居中度(standoff)剖面。

    基于呼103设计居中度代理值估算：
    - 20313.doc 6.3 软件模拟设计居中度 77.8%（无实测居中度曲线）
    - 不引入扶正器明细，只做轻微分段修正

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
    nominal_clearance_mm = max(nominal_bit_mm - liner_od_mm, 5.0)
    standoff *= min(max(nominal_clearance_mm / clearance_mm, 0.92), 1.03)
    return min(max(standoff, 0.30), 0.90)


def _build_hole_profile(profile_rows: Iterable[tuple[float, float, float]]) -> tuple[DepthValuePoint, ...]:
    return tuple(DepthValuePoint(depth_md_m=depth, value=hole_diameter) for depth, hole_diameter, _ in profile_rows)


def _build_inclination_profile(profile_rows: Iterable[tuple[float, float, float]]) -> tuple[DepthValuePoint, ...]:
    return tuple(DepthValuePoint(depth_md_m=depth, value=inclination) for depth, _, inclination in profile_rows)


def _build_liner_od_profile(profile_rows: Iterable[tuple[float, float, float]]) -> tuple[DepthValuePoint, ...]:
    """构建尾管外径剖面（复合尾管：上段 168.3mm / 下段 139.7mm，变扣 7330.694m）。"""
    return tuple(
        DepthValuePoint(depth_md_m=depth, value=_liner_od_at(depth))
        for depth, hole_diameter, _ in profile_rows
    )


def _build_standoff_profile(profile_rows: Iterable[tuple[float, float, float]]) -> tuple[DepthValuePoint, ...]:
    return tuple(
        DepthValuePoint(
            depth_md_m=depth,
            value=_standoff_value(depth, hole_diameter, _liner_od_at(depth), _nominal_bit_diameter_mm(depth)),
        )
        for depth, hole_diameter, _ in profile_rows
    )


def load_hu103_tailpipe(
    *,
    caliper_csv_path: Path | None = None,
    reference_root: Path | None = None,
) -> tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载呼103完井尾管段"设计版"模型输入（默认，20313.doc 5.2.4）。

    Args:
        caliper_csv_path: 可选井径/井斜 CSV 路径（默认指向现场实测分段曲线合并 CSV）。
        reference_root: 可选参考资料根目录。

    Returns:
        井筒参数（复合尾管 dual-diameter）、流体参数、设计版泵注程序与验证资料路径。

    Note:
        实际施工版（20314.doc 施工记录表，全程 0.6m³/min）见 load_hu103_tailpipe_actual。
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
        liner_wall_thickness_mm=HU103_LINER_WALL_THICKNESS_MM,
        hole_diameter_profile=_build_hole_profile(profile_rows),
        inclination_profile=_build_inclination_profile(profile_rows),
        standoff_profile=_build_standoff_profile(profile_rows),
        liner_od_profile=_build_liner_od_profile(profile_rows),
        evaluation_windows=(
            EvaluationWindow(name="CBL评价井段(139.7mm段)", top_md_m=7338.0, bottom_md_m=7712.0, window_type="cbl"),
            EvaluationWindow(name="CBL评价井段(168.3mm段)", top_md_m=5540.0, bottom_md_m=7330.6, window_type="cbl"),
            EvaluationWindow(name="全尾管段", top_md_m=5536.662, bottom_md_m=7770.0, window_type="custom"),
            # 地层目标（target_intervals.csv，20313.doc 井史 2.3 油气显示，design/field_measured）。
            EvaluationWindow(name="气层-K1q(地层目标)", top_md_m=7507.0, bottom_md_m=7557.0, window_type="formation_target"),
            EvaluationWindow(name="差气层-K1q(地层目标)", top_md_m=7476.0, bottom_md_m=7479.0, window_type="formation_target"),
            EvaluationWindow(name="差气层-J3k(地层目标)", top_md_m=7563.0, bottom_md_m=7619.0, window_type="formation_target"),
            EvaluationWindow(name="高压水层-K1s(地层目标)", top_md_m=6120.0, bottom_md_m=6122.0, window_type="formation_target"),
        ),
        reference_root=resolved_reference_root,
        # 复合尾管：上层 168.3mm 段字段（dual-diameter），与 liner_od_profile 共同表达复合结构
        upper_section_bottom_md_m=HU103_UPPER_SECTION_BOTTOM_MD_M,
        upper_liner_od_mm=HU103_UPPER_LINER_OD_MM,
        upper_liner_id_mm=HU103_UPPER_LINER_ID_MM,
        notes=(
            "呼103井 168.3+139.7mm 复合油层尾管固井段（5536.662–7770m），不含回接固井与其他套管段。",
            "悬挂器 5536.662–5545.972m（NSSX-CYFMD 封隔一体式）；上段 168.3mm(5546.022–7330.694m)、下段 139.7mm(7330.694–7768.51m)。",
            "139.7mm 段 ID 107.94mm（壁厚15.88mm）、168.3mm 段 ID 138.9mm（壁厚14.7mm）为实测/设计口径（20313.doc/20315.doc）。",
            "井径/井斜用现场实测分段曲线（20313.doc 1.5.1/1.5.2）；5536.662–5755m 重叠段无实测，取首测点等效外推。",
            "居中度仅设计模拟值 77.8%（20313.doc 6.3），无实测，模型按 proxy 处理。",
            "复合尾管由 dual-diameter 上层字段 + liner_od_profile 表达；LEGACY(2026-08 前) 为 139.7mm 单外径代理模型。",
            "2026-08-29 校准补记：168.3mm 段实际下深 5545.972–7330.604m（20318，与设计口径 5546.022–7330.694m 并存）；"
            "回接筒 5536.662–5539.762m（OD 226/ID 202）。",
        ),
    )

    # 2026-08 重建非水泥浆流体：密度/角色按现场值（平衡液=mud/ρ1.88、隔离液1/2=spacer/2.00/1.95、
    # 压塞液=2.10、中置液=1.95、替浆=钻井液1.98）；流变无实测者复用钻井液或代理并标注。
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
            role=FluidRole.MUD,  # 现场角色=mud（轻泥浆，20313.doc/20314.doc）
            density_kg_m3=HU103_BALANCE_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU103_BALANCE_PV_PA_S,
            yield_stress_pa=HU103_BALANCE_YP_PA,
        ),
        FluidSpec(
            name="隔离液1",
            role=FluidRole.SPACER,
            density_kg_m3=HU103_SPACER1_DENSITY_KG_M3,
            rheology_model=RheologyModel.POWER_LAW,
            power_law_n=HU103_SPACER_POWER_LAW_N,
            consistency_k=HU103_SPACER_CONSISTENCY_K,
        ),
        FluidSpec(
            name="隔离液2",
            role=FluidRole.SPACER,
            density_kg_m3=HU103_SPACER2_DENSITY_KG_M3,
            rheology_model=RheologyModel.POWER_LAW,
            power_law_n=HU103_SPACER_POWER_LAW_N,
            consistency_k=HU103_SPACER_CONSISTENCY_K,
        ),
        FluidSpec(
            name="隔离液(实际)",
            role=FluidRole.SPACER,
            density_kg_m3=HU103_SPACER_ACTUAL_DENSITY_KG_M3,
            rheology_model=RheologyModel.POWER_LAW,
            power_law_n=HU103_SPACER_POWER_LAW_N,
            consistency_k=HU103_SPACER_CONSISTENCY_K,
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
            name="中置液",
            role=FluidRole.DISPLACEMENT,
            density_kg_m3=HU103_MIDDLE_FLUID_DENSITY_KG_M3,
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

    # 设计版泵注程序（20313.doc 5.2.4 注入浆柱结构/替浆浆柱结构）
    schedule = PumpingSchedule(
        steps=(
            PumpingScheduleStep(
                step_name="注入平衡液(轻泥浆)",
                fluid_name="平衡液",
                volume_m3=HU103_BALANCE_VOLUME_M3,
                rate_m3_min=HU103_DESIGN_RATE_FRONT_M3_MIN,
                remarks="设计: 平衡液28m³，ρ1.88，环空占高800m(260mm井眼)。",
            ),
            PumpingScheduleStep(
                step_name="注入驱油隔离液1",
                fluid_name="隔离液1",
                volume_m3=HU103_SPACER1_VOLUME_M3,
                rate_m3_min=HU103_DESIGN_RATE_FRONT_M3_MIN,
                remarks="设计: 隔离液1 17.5m³，ρ2.00，占高500m。",
            ),
            PumpingScheduleStep(
                step_name="注入驱油隔离液2",
                fluid_name="隔离液2",
                volume_m3=HU103_SPACER2_VOLUME_M3,
                rate_m3_min=HU103_DESIGN_RATE_FRONT_M3_MIN,
                remarks="设计: 隔离液2 17.5m³，ρ1.95，占高500m。",
            ),
            PumpingScheduleStep(
                step_name="注入领浆",
                fluid_name="领浆",
                volume_m3=HU103_LEAD_VOLUME_M3,
                rate_m3_min=HU103_DESIGN_RATE_CEMENT_M3_MIN,
                remarks="设计: 领浆22m³，ρ2.05，封固5350-6000m。",
            ),
            PumpingScheduleStep(
                step_name="注入中间浆",
                fluid_name="中间浆",
                volume_m3=HU103_INTERMEDIATE_VOLUME_M3,
                rate_m3_min=HU103_DESIGN_RATE_CEMENT_M3_MIN,
                remarks="设计: 中间浆35m³，ρ2.05，封固6000-7000m（2026-08-29 校准，0708 多数口径 35；LEGACY 36 系浆柱表孤值）。",
            ),
            PumpingScheduleStep(
                step_name="注入尾浆",
                fluid_name="尾浆",
                volume_m3=HU103_TAIL_VOLUME_M3,
                rate_m3_min=HU103_DESIGN_RATE_CEMENT_M3_MIN,
                remarks="设计: 尾浆26m³，ρ2.05，封固7000-7770m；旧记'尾浆 19.0（含下塞 1.18）'系 legacy 残留、0708 无出处（实际尾浆 26、下塞 0.67，2026-08-29 校准）。",
            ),
            PumpingScheduleStep(
                step_name="注入压塞液",
                fluid_name="压塞液",
                volume_m3=HU103_PLUG_VOLUME_M3,
                rate_m3_min=HU103_DESIGN_RATE_CEMENT_M3_MIN,
                remarks="设计: 压塞液7m³，ρ2.10，管内占高582m。",
            ),
            PumpingScheduleStep(
                step_name="替钻井液(一段)",
                fluid_name="替浆液",
                volume_m3=HU103_DISP_MUD1_VOLUME_M3,
                rate_m3_min=HU103_DESIGN_RATE_CEMENT_M3_MIN,
                remarks="设计: 替钻井液23m³，管内占高1520m。",
            ),
            PumpingScheduleStep(
                step_name="替中置液",
                fluid_name="中置液",
                volume_m3=HU103_DISP_MIDDLE_VOLUME_M3,
                rate_m3_min=HU103_DESIGN_RATE_CEMENT_M3_MIN,
                remarks="设计: 中置液15m³，ρ1.95，占高1600m。",
            ),
            PumpingScheduleStep(
                step_name="替钻井液(二段)",
                fluid_name="替浆液",
                volume_m3=HU103_DISP_MUD2_VOLUME_M3,
                rate_m3_min=HU103_DESIGN_RATE_CEMENT_M3_MIN,
                remarks="设计: 替钻井液46.2m³，设计排量1.8-0.8变排量（此处取名义1.8）。",
            ),
        ),
        notes=(
            "设计版泵注程序（20313.doc 5.2.4）：平衡液→隔离液1→隔离液2→领浆→中间浆→尾浆→压塞液→替钻井液→替中置液→替钻井液。",
            "设计排量 1.6(前置)/1.8(水泥及替浆)m³/min；实际施工因环空憋堵降至 0.6m³/min，见 load_hu103_tailpipe_actual。",
            "尾浆设计 26m³；旧注'采用值对照表记 tail=19.0（含下塞1.18m³）'系 legacy 残留、0708 无出处（实际尾浆 26、下塞 0.67，2026-08-29 校准标注）。",
            "设计替浆总量 83.2m³（23+15+45.2；+压塞液7=90.2 与 20313 七.1.4 总量 90.2 自洽）；LEGACY(2026-08-29 前): 84.2（46.2 排量表笔误）。",
            "中间浆设计 35m³（0708 多数口径 4 处：7.1.1 用量表总量 83=22+35+26、7.2 流程、实际 35）；LEGACY(2026-08-29 前): 36 系浆柱表孤值。",
            "LEGACY(2026-08 前): 平衡液18/冲洗液14/隔离液18/领浆91/中间浆10/尾浆19/压塞液5/替浆110 代理程序已重建。",
        ),
    )

    validation_data = _hu103_validation_data(resolved_reference_root)
    return well_spec, fluids, schedule, validation_data


def _hu103_validation_data(reference_root: Path) -> ValidationData:
    """呼103 校验资料集合（设计版与实际版共用）。"""
    return ValidationData(
        cbl_summary_path=reference_root / "呼103井_CBL评价报告.pdf",
        cbl_pass_rate=0.1206,
        notes=(
            "139.7mm 尾管段（7338.0–7712.0m）CBL 合格率 12.06%（100516.pdf），不合格；结构化 cbl_pass_rate=0.1206。",
            "不得混用 168.3mm 段 0.04%（100515.pdf）；旧注'整段加权综合 6.05%'系简单平均误标，"
            "正确长度加权综合为 2.12%（(0.04×1790.6+12.06×374)/2164.6，2026-08-29 校准）。",
            "cbl_summary_path 指向参考文档/呼103下的 CBL PDF；若原始文件名不同，后续仅需更新路径。",
            "2026-08 已按现场提取包重建几何/剖面/泵注/非水泥浆流体；水泥浆流变沿用实验报告（203111.docx）。",
        ),
    )


def load_hu103_tailpipe_actual(
    *,
    caliper_csv_path: Path | None = None,
    reference_root: Path | None = None,
) -> tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载呼103完井尾管段"实际施工"版输入（20314.doc 施工记录表）。

    与设计版（load_hu103_tailpipe）的差异仅在 PumpingSchedule：
    - 实际全程排量 0.6m³/min（环空憋堵迫使降排量，约为设计 1.8m³/min 的 1/3）；
    - 实际平衡液 48m³、隔离液 35m³、压塞液 3.0m³（ρ2.20）、替浆 90.9m³；
    - 实际水泥浆为"三凝"体系且未记录分相体积，分相 22/35/26 保守沿用设计版。

    Returns:
        井筒参数（同设计版）、流体参数（同设计版）、实际施工泵注程序、校验资料。
    """

    well_spec, fluids, _, _ = load_hu103_tailpipe(
        caliper_csv_path=caliper_csv_path,
        reference_root=reference_root,
    )
    resolved_reference_root = reference_root or DEFAULT_REFERENCE_ROOT

    # 实际版泵注程序（20314.doc 施工记录表）
    schedule = PumpingSchedule(
        steps=(
            PumpingScheduleStep(
                step_name="注入平衡液(实际)",
                fluid_name="平衡液",
                volume_m3=HU103_ACTUAL_BALANCE_VOLUME_M3,
                rate_m3_min=HU103_ACTUAL_RATE_M3_MIN,
                remarks="实际: 平衡液48m³@0.6，ρ1.88，泵压8-12MPa。",
            ),
            PumpingScheduleStep(
                step_name="注入隔离液(实际)",
                fluid_name="隔离液(实际)",
                volume_m3=HU103_ACTUAL_SPACER_VOLUME_M3,
                rate_m3_min=HU103_ACTUAL_RATE_M3_MIN,
                remarks="实际: 隔离液35m³@0.6，ρ1.92，泵压15MPa。",
            ),
            PumpingScheduleStep(
                step_name="注入领浆(三凝)",
                fluid_name="领浆",
                volume_m3=HU103_LEAD_VOLUME_M3,
                rate_m3_min=HU103_ACTUAL_RATE_M3_MIN,
                remarks="实际水泥浆为三凝体系，未记录分相体积，分相 22m³ 沿用设计版。",
            ),
            PumpingScheduleStep(
                step_name="注入中间浆(三凝)",
                fluid_name="中间浆",
                volume_m3=HU103_INTERMEDIATE_VOLUME_M3,
                rate_m3_min=HU103_ACTUAL_RATE_M3_MIN,
                remarks="实际中间浆分相 36m³ 沿用设计版（三凝未分相记录）。",
            ),
            PumpingScheduleStep(
                step_name="注入尾浆(三凝)",
                fluid_name="尾浆",
                volume_m3=HU103_TAIL_VOLUME_M3,
                rate_m3_min=HU103_ACTUAL_RATE_M3_MIN,
                remarks="实际尾浆分相 26m³ 沿用设计版（三凝未分相记录）。",
            ),
            PumpingScheduleStep(
                step_name="注入压塞液(实际)",
                fluid_name="压塞液",
                volume_m3=3.0,
                rate_m3_min=HU103_ACTUAL_RATE_M3_MIN,
                remarks="实际: 压塞液 3.0m³@0.6，平均密度 2.20（20315'压塞3.0m3 平均密度2.20'；20314，2026-08-29 校准补入）。",
            ),
            PumpingScheduleStep(
                step_name="替浆(实际)",
                fluid_name="替浆液",
                volume_m3=HU103_ACTUAL_DISPLACEMENT_VOLUME_M3,
                rate_m3_min=HU103_ACTUAL_RATE_M3_MIN,
                remarks="实际: 替浆90.9m³@0.6，泵压13-9MPa；到量碰压、放回水断流。",
            ),
        ),
        notes=(
            "实际版泵注程序（20314.doc 施工记录表）：平衡液48 + 隔离液35 + 水泥浆(三凝) + 压塞液3.0(ρ2.20) + 替浆90.9，全程 0.6m³/min。",
            "实际排量仅为设计 1.8m³/min 的 1/3，直接原因是下尾管后环空憋堵（20315.doc（三）井循环过程）。",
            "实际水泥浆为三凝体系且未记录分相体积，分相 22/35/26 保守沿用设计版，未编造。",
            "压塞液实际 3.0m³@0.6、平均密度 2.20（20315/20314，2026-08-29 校准补入）。",
            "环空初始液保持 1.98（设计口径）；实际固井前 2.02/PV40/YP7（20314），若按实际口径可换（2026-08-29 注记）。",
            "候凝：环空加压 5.3MPa 候凝 48h；碰压后放回水断流（20314.doc）。",
        ),
    )
    validation_data = _hu103_validation_data(resolved_reference_root)
    return well_spec, fluids, schedule, validation_data


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
