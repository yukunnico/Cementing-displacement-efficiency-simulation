"""
呼探1井（hu1）139.7mm 油层尾管段标准数据加载器

本模块把呼探1井现场提取包（参考文档/现场资料提取/hu1_呼探1/）中的尾管段资料
整理为 cemdisp 标准输入结构。

================================================================================
2026-08-16 整体重建声明
================================================================================
旧版本文件为 HT1-001（呼探1-001）误复制件：TD=7746、hanger=5469.711、
168.3+139.7 双径尾管、五段水泥泵注（40/20/20.6/28.7/22.1）等全部是 HT1-001
井的值，与呼探1现场（TD 7601 / 单 139.7mm 尾管 / 泵注冲洗液15+隔离液15+
领浆24+尾浆36+替浆87）完全不符。**旧版驱动的所有模拟结果作废。**

本次整体重写依据：`参考文档/现场资料提取/hu1_呼探1/`（16 CSV，唯一现场值来源）
与核对文档 `00_执行记录/其他井loader核对_2026-08-16.md` §2（hu1 节）。

井号身份注记：呼探1 与 呼探1-001（HT1-001）是否同井存疑，待现场裁定。
本模块按 hu1_呼探1 提取包独立建井，与 ht1_001_loader 分列，不互相引用。

物理参数要点（2026-08-16 重建）：
- 井段范围: 3523.27m（悬挂器设计位置/尾管顶） - 7601.00m（完钻/浮鞋）
- 悬挂器: 设计 3523.27m / 实际 3532.21m（双口径，模型取设计值）
- 单一直径 139.7mm 尾管（无 168.3mm 段）：
  上部薄壁段 12.09mm（3523-5300m，ID 115.52）+ 下部厚壁段 14.27mm（5300-7601m，ID 111.16）
- 五开钻头 190.5mm；裸眼井径 96 点现场实测（5694-7601m，均值 194.33mm，最大 233.34@5780m）
- 井斜 90 点现场实测（5694-7601m，3.04°→15.09°，均值 7.34°）
- 泵注重建（204151.doc 施工总结）：冲洗液15 + 加重隔离液15 + 领浆24 + 尾浆36 +
  后置液2 + 替浆87（管串36 + 环空49），全程排量 0.65→0.36 m³/min
- 流变标注：钻井液 PV/YP 为现场实测（204111.doc）；水泥浆 n/K 无实测记录，
  以同区块同密度体系实验值作 proxy；前置液/后置液/替浆液流变亦为 proxy
- 居中度：现场无实测（扶正器 136 只、间距约 33m），standoff 为 model_assumption
- CBL：仅定性分段（3 张 CBL/VDL 图，7474-7514m 不合格 40m，整体"合格"），
  无数值合格率 → cbl_pass_rate=None，不编造
"""

from __future__ import annotations

import csv
from pathlib import Path

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep
from cemdisp.data.validation_data import ValidationData
from cemdisp.data.well_spec import DepthValuePoint, EvaluationWindow, WellSpec
from cemdisp.data.provenance import WELL_PROVENANCE
from cemdisp.models2d.boundary_bridge import build_sync_card
from cemdisp.transport1d.casing_flow import CasingFlowSolver


PROJECT_ROOT = Path(__file__).resolve().parents[3]
# 原始参考目录保留用于 job_report_path 等原始资料追溯（旧数据包，2026 年前构建）。
DEFAULT_REFERENCE_ROOT = PROJECT_ROOT / "参考文档" / "呼探1"
# 现场提取包（唯一现场值来源）：
# - caliper_profile.csv: 204111.doc 电测井径 96 点（5694–7601m，20m 间隔）
# - inclination_profile.csv: 204111.doc 连续测斜 90 点（5694–7601m）
DEFAULT_EXTRACTION_ROOT = PROJECT_ROOT / "参考文档" / "现场资料提取" / "hu1_呼探1"
DEFAULT_CALIPER_CSV = DEFAULT_EXTRACTION_ROOT / "caliper_profile.csv"
DEFAULT_INCLINATION_CSV = DEFAULT_EXTRACTION_ROOT / "inclination_profile.csv"
DEFAULT_CBL_EVALUATION_CSV = DEFAULT_EXTRACTION_ROOT / "cbl_evaluation.csv"

# ===========================================================================
# 呼探1井段几何参数（单位：m, mm）—— 与核对报告 §2.1 现场值一致
# 来源：well_basic_info.csv / well_geometry.csv / casing_liner_string.csv
# ===========================================================================
HU1_WELL_NAME = "呼探1井"
HU1_DRILLED_DEPTH_MD_M = 7601.0        # 实际完钻井深（well_basic_info：2020-09-26 钻至 7601.00m 完钻）
HU1_BOTTOM_MD_M = HU1_DRILLED_DEPTH_MD_M
HU1_SHOE_MD_M = HU1_DRILLED_DEPTH_MD_M  # 浮鞋 7601m（casing_liner_string；139.7mm 尾管下深 7600.50m）
# 悬挂器：设计 3523.27m（204111.doc 设计悬挂器位置）/ 实际 3532.21m（10068.doc 实际悬挂器位置）。
# 模型取设计值作剖面顶，实际值在 notes 双口径标注。
HU1_HANGER_MD_M = 3523.27
HU1_TOP_MD_M = HU1_HANGER_MD_M
# 上层套管结构（等效井径用）：三开 273.1mm 技套 OD，壁厚 13.84 → ID 245.42；
# 四开 219.1mm 技术尾管 ID 190.5（well_geometry casing_id_mm=190.5），其鞋 5693.89m。
# casing_id_mm 字段语义在本仓库为上层套管"外径"（同上 hu102/hu103/ht1_003，命名历史遗留）。
HU1_CASING_ID_MM = 273.1                # 三开技术套管外径（名义）
HU1_SURFACE_CASING_ID_MM = 245.42       # 273.1mm 技套内径 = 273.1 - 2×13.84
HU1_INTERMEDIATE_LINER_TOP_MD_M = 3623.27  # 四开 219.1mm 技术尾管顶
HU1_INTERMEDIATE_LINER_ID_MM = 190.5    # 219.1mm 技术尾管内径（= 五开钻头尺寸）
HU1_OPEN_HOLE_TOP_MD_M = 5694.1         # 裸眼段顶（219.1 尾管鞋 5693.89m；井径首测点 5694.1m）
# 单一直径 139.7mm 尾管：上部薄壁段 12.09mm（3523–5300m）/ 下部厚壁段 14.27mm（5300–7601m）
# （204111.doc 套管规格；casing_liner_string 尾部套管 12.09/14.27 两段）。
HU1_LINER_OD_MM = 139.7
HU1_LINER_WALL_THICKNESS_MM = 14.27     # 厚壁段壁厚（下段；作参考口径）
HU1_UPPER_LINER_WALL_THICKNESS_MM = 12.09  # 薄壁段壁厚（上段 3523–5300m）
HU1_LOWER_LINER_ID_MM = HU1_LINER_OD_MM - 2.0 * HU1_LINER_WALL_THICKNESS_MM   # 111.16
HU1_UPPER_LINER_ID_MM = HU1_LINER_OD_MM - 2.0 * HU1_UPPER_LINER_WALL_THICKNESS_MM  # 115.52
# liner_id 单值口径：取厚壁段 ID 111.16（与 hu102/hu103 一致的简化）。
# 注：well_geometry 标 casing_id_mm=107.98（疑为 15.88 壁厚复制误差），现场内部不一致，
# 按 204111.doc 壁厚计算口径 111.16 采用（厚壁段）。
HU1_LINER_ID_MM = HU1_LOWER_LINER_ID_MM
HU1_BIT_DIAMETER_MM = 190.5             # 五开钻头尺寸（裸眼名义井径）
HU1_CENTRALIZER_COUNT = 136             # 现场整体弹性扶正器数量（extraction_notes，间距约 33m）

# ===========================================================================
# 呼探1流体参数（单位密度 kg/m³；来源 204111.doc 固井设计 / 204151.doc 施工总结）
# 密度全为 field_measured；流变标注见各常量注释。
# ===========================================================================
HU1_MUD_DENSITY_KG_M3 = 2120.0          # 固井前钻井液（油基）2.12（204111.doc，field_measured）
HU1_MUD_PV_PA_S = 0.066                 # 塑性粘度 66mPa·s（204111.doc，field_measured）
HU1_MUD_YP_PA = 10.0                    # 动切力 10Pa（204111.doc，field_measured）
HU1_WASH_DENSITY_KG_M3 = 2000.0         # 冲洗型隔离液（冲洗液）2.00（204151.doc，field_measured）
HU1_SPACER_DENSITY_KG_M3 = 2050.0       # 加重隔离液 2.05（204151.doc，field_measured）
# 冲洗液/加重隔离液流变：现场无实测六速数据 → proxy（沿用同区块呼102 隔离液代理 PV35/YP8）。
HU1_WASH_PV_PA_S = 0.035                # proxy（无实测）
HU1_WASH_YP_PA = 8.0                    # proxy（无实测）
HU1_SPACER_PV_PA_S = 0.035              # proxy（无实测）
HU1_SPACER_YP_PA = 8.0                  # proxy（无实测）
HU1_LEAD_DENSITY_KG_M3 = 2100.0         # 领浆 2.10（204151.doc/204111.doc，field_measured；3523–5300m 段）
HU1_TAIL_DENSITY_KG_M3 = 1900.0         # 尾浆 1.90（204151.doc/204111.doc，field_measured；5300–7601m 段）
# 领/尾浆流变：20441.doc 稠化/强度/污染实验完整，但无六速幂律 n/K 实测记录
# （rheometer_readings 均为污染/稠化实验，无 n/K）→ proxy（同区块同密度体系实验值）：
#   领浆 ρ2.10 → 呼102 领浆化验 n=0.737/K=0.947（20234.doc，同 G 级 D 系列高密度体系）；
#   尾浆 ρ1.90 → 呼1-001 尾浆化验 n=0.886/K=0.453（fluid_properties，同密度 1.90）。
HU1_LEAD_POWER_LAW_N = 0.737            # proxy（呼102 领浆化验值，无呼探1 实测）
HU1_LEAD_CONSISTENCY_K = 0.947          # proxy（同上）
HU1_TAIL_POWER_LAW_N = 0.886            # proxy（呼1-001 尾浆化验值，同密度）
HU1_TAIL_CONSISTENCY_K = 0.453          # proxy（同上）
HU1_PLUG_DENSITY_KG_M3 = 2050.0         # 后置液 2.05（204151.doc 后置液 2m³，field_measured）
# 后置液流变：无实测 → proxy（复用隔离液 PV35/YP8）。
HU1_PLUG_PV_PA_S = 0.035                # proxy（无实测）
HU1_PLUG_YP_PA = 8.0                    # proxy（无实测）
HU1_DISPLACEMENT1_DENSITY_KG_M3 = 1950.0  # 替浆液1（管串）1.95（204151.doc，field_measured，36m³）
HU1_DISPLACEMENT2_DENSITY_KG_M3 = 2080.0  # 替浆液2（环空）2.08（204151.doc，field_measured，49m³）
# 替浆液流变：无实测 → proxy（复用固井前钻井液实测 PV66/YP10，替浆液为钻井液型）。
HU1_DISPLACEMENT_PV_PA_S = 0.066        # proxy（复用钻井液）
HU1_DISPLACEMENT_YP_PA = 10.0           # proxy（复用钻井液）

# ===========================================================================
# 呼探1施工程序参数（204151.doc 施工总结 / pumping_schedule.csv / construction_events.csv）
# 现场实际七步 + 替浆三分段；无独立"设计版"，默认加载即为实际施工版。
# ===========================================================================
HU1_WASH_VOLUME_M3 = 15.0               # 冲洗液 15m³（实际）
HU1_SPACER_VOLUME_M3 = 15.0             # 加重隔离液 15m³（实际）
HU1_LEAD_VOLUME_M3 = 24.0               # 领浆 24m³（实际，3523–5300m 段）
HU1_TAIL_VOLUME_M3 = 36.0               # 尾浆 36m³（实际，14+22，5300–7601m 段）
HU1_PLUG_VOLUME_M3 = 2.0                # 后置液 2m³（替浆段首段，ρ2.05）
HU1_DISPLACEMENT1_VOLUME_M3 = 36.0      # 替浆液1（管串）36m³（ρ1.95）
HU1_DISPLACEMENT2_VOLUME_M3 = 49.0      # 替浆液2（环空）49m³（ρ2.08）
HU1_DISPLACEMENT_TOTAL_VOLUME_M3 = 87.0  # 替浆总量 87m³ = 2 + 36 + 49（现场口径）
# 排量（204151.doc）：前置/水泥 0.65；替浆段 0.65→0.5→0.36→0.42 递减（240min）。
HU1_FRONT_RATE_M3_MIN = 0.65            # 冲洗液/隔离液/领浆/尾浆排量
HU1_PLUG_RATE_M3_MIN = 0.5              # 后置液排量（替浆段内代表性值）
HU1_DISPLACEMENT1_RATE_M3_MIN = 0.5     # 替浆液1（管串）代表性排量
HU1_DISPLACEMENT2_RATE_M3_MIN = 0.36    # 替浆液2（环空）代表性排量（末段 0.42 为碰压前）

# 替浆时间校验：2/0.5 + 36/0.5 + 49/0.36 = 4 + 72 + 136 ≈ 212min，与现场 240min 同量级（现场含降排量/停泵）。


def _read_caliper_rows(caliper_csv_path: Path) -> tuple[tuple[float, float], ...]:
    """读取现场提取包井径剖面 CSV，按 md 升序返回 (md_m, caliper_mm)。

    兼容两套列名：
    - 呼探1提取包: measured_depth_m / avg_caliper_mm
    - ht1_003/hu102 提取包: md_m / caliper_mm
    """
    rows: list[tuple[float, float]] = []
    with caliper_csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            md = row.get("md_m", row.get("measured_depth_m"))
            cal = row.get("caliper_mm", row.get("avg_caliper_mm"))
            if md is None or cal is None:
                continue
            rows.append((float(md), float(cal)))
    if not rows:
        raise ValueError(f"呼探1 井径 CSV 为空: {caliper_csv_path}")
    return tuple(sorted(rows))


def _read_inclination_rows(inclination_csv_path: Path) -> tuple[tuple[float, float], ...]:
    """读取现场提取包井斜剖面 CSV，按 md 升序返回 (md_m, inclination_deg)。"""
    rows: list[tuple[float, float]] = []
    with inclination_csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            md = row.get("md_m", row.get("measured_depth_m"))
            inc = row.get("inclination_deg")
            if md is None or inc is None:
                continue
            rows.append((float(md), float(inc)))
    if not rows:
        raise ValueError(f"呼探1 井斜 CSV 为空: {inclination_csv_path}")
    return tuple(sorted(rows))


def _depth_points(values: tuple[tuple[float, float], ...]) -> tuple[DepthValuePoint, ...]:
    """把测深-数值元组转换为 WellSpec 使用的剖面点。"""
    return tuple(DepthValuePoint(depth_md_m=depth, value=value) for depth, value in values)


def _build_hole_profile(caliper_rows: tuple[tuple[float, float], ...]) -> tuple[DepthValuePoint, ...]:
    """构建井径剖面。

    3523.27–5694m 段为上层套管柱内重叠段（139.7mm 尾管在套管内），无裸眼井径测点，
    按上层套管内径等效（273.1mm 技套段 ID 245.42 → 219.1mm 尾管段 ID 190.5）；
    5694–7601m 段为裸眼，96 点现场实测。
    """
    points: list[tuple[float, float]] = [
        (HU1_TOP_MD_M, HU1_SURFACE_CASING_ID_MM),
        (HU1_INTERMEDIATE_LINER_TOP_MD_M, HU1_SURFACE_CASING_ID_MM),
        (HU1_INTERMEDIATE_LINER_TOP_MD_M + 0.001, HU1_INTERMEDIATE_LINER_ID_MM),
        (HU1_OPEN_HOLE_TOP_MD_M - 0.001, HU1_INTERMEDIATE_LINER_ID_MM),
    ]
    for md, cal in caliper_rows:
        points.append((md, cal))
    if points[-1][0] < HU1_BOTTOM_MD_M:
        points.append((HU1_BOTTOM_MD_M, points[-1][1]))
    return _depth_points(tuple(points))


def _build_inclination_profile(incl_rows: tuple[tuple[float, float], ...]) -> tuple[DepthValuePoint, ...]:
    """构建井斜剖面。

    3523.27–5694m 段无连续测斜（悬挂器以上套管段），以 5694m 首测点值 3.04° 近似
    （model_assumption）；5694–7601m 为 90 点现场实测（3.04°→15.09°）。
    """
    first_inc = incl_rows[0][1] if incl_rows else 0.0
    points: list[tuple[float, float]] = [
        (HU1_TOP_MD_M, first_inc),
        (HU1_OPEN_HOLE_TOP_MD_M - 0.001, first_inc),
    ]
    for md, inc in incl_rows:
        points.append((md, inc))
    if points[-1][0] < HU1_BOTTOM_MD_M:
        points.append((HU1_BOTTOM_MD_M, points[-1][1]))
    return _depth_points(tuple(points))


def _build_standoff_profile() -> tuple[DepthValuePoint, ...]:
    """构建居中度(standoff)剖面 —— model_assumption 均一值。

    现场无居中度实测（centralization_profile.csv data_type=missing；扶正器 136 只、
    间距约 33m 为 extraction_notes 记录）。取均一 0.65（保守代理，标注 model_assumption）；
    参照邻井 hu102（0.65–0.70）、呼1-003（0.83 设计值）量级。
    """
    return _depth_points(
        (
            (HU1_TOP_MD_M, 0.65),
            (HU1_BOTTOM_MD_M, 0.65),
        )
    )


def _build_well_spec(
    resolved_reference_root: Path,
    caliper_csv_path: Path,
    inclination_csv_path: Path,
) -> WellSpec:
    """构建呼探1井井筒规格（几何/剖面/评价窗）。"""
    caliper_rows = _read_caliper_rows(caliper_csv_path)
    incl_rows = _read_inclination_rows(inclination_csv_path)

    return WellSpec(
        well_name=HU1_WELL_NAME,
        top_md_m=HU1_TOP_MD_M,
        bottom_md_m=HU1_BOTTOM_MD_M,
        shoe_md_m=HU1_SHOE_MD_M,
        hanger_md_m=HU1_HANGER_MD_M,
        casing_id_mm=HU1_CASING_ID_MM,
        liner_od_mm=HU1_LINER_OD_MM,
        liner_id_mm=HU1_LINER_ID_MM,
        liner_wall_thickness_mm=HU1_LINER_WALL_THICKNESS_MM,
        hole_diameter_profile=_build_hole_profile(caliper_rows),
        inclination_profile=_build_inclination_profile(incl_rows),
        standoff_profile=_build_standoff_profile(),
        evaluation_windows=(
            # CBL 评价窗：电测测量段 3530–7601m（construction_events.csv；CBL/VDL 图覆盖 5630–7607m）。
            EvaluationWindow(name="CBL评价井段", top_md_m=3530.0, bottom_md_m=7601.0, window_type="cbl"),
            # CBL 分段质量（cbl_evaluation.csv 定性 16 段中取最差两段，interpreted；完整分段见提取包 CSV）。
            EvaluationWindow(name="CBL质量段(不合格40m)", top_md_m=7474.0, bottom_md_m=7514.0, window_type="cbl_quality"),
            EvaluationWindow(name="CBL质量段(差)", top_md_m=7415.0, bottom_md_m=7430.0, window_type="cbl_quality"),
            # 目标窗：target_intervals.csv 明确"139.7mm尾管目标段 missing"（现场无地层目标窗），
            # 仅给出模型聚焦窗（油层尾管封固段 3523.27–7601m），勿写成现场目标。
            EvaluationWindow(name="尾管封固段(模型聚焦窗)", top_md_m=3523.27, bottom_md_m=7601.0, window_type="model_focus"),
        ),
        reference_root=resolved_reference_root,
        notes=(
            "呼探1井（hu1）为单 139.7mm 油层尾管固井，无 168.3mm 段：上部薄壁段 12.09mm（3523–5300m，ID 115.52）/ "
            "下部厚壁段 14.27mm（5300–7601m，ID 111.16），liner_id_mm 取厚壁段 111.16 单值口径。",
            "悬挂器设计 3523.27m / 实际 3532.21m（双口径，模型取设计值作剖面顶）；尾管下深 7600.50m。",
            "井径剖面：3523–5694m 上层套管柱重叠段按等效内径（273.1 技套 ID 245.42 → 219.1 尾管 ID 190.5），"
            "5694–7601m 为 96 点现场实测（均值 194.33mm，最大 233.34@5780m，最小 185.27@7540m，field_measured）。",
            "井斜剖面：5694–7601m 为 90 点现场实测（3.04°→15.09°，均值 7.34°）；3523–5694m 无连续测斜，"
            "以 5694m 首测点 3.04° 近似（model_assumption）。",
            "居中度 standoff 无现场实测（扶正器 136 只、间距约 33m 为 extraction_notes 记录；"
            "centralization_profile.csv 标 missing），取均一 0.65，model_assumption。",
            "井号身份注记：呼探1 与 呼探1-001（HT1-001）是否同井存疑，待现场裁定；"
            "本 loader 按 hu1_呼探1 提取包独立建井，与 ht1_001_loader 分列。",
            "CBL 评价窗 3530–7601m 为电测测量段（construction_events.csv）；target_intervals.csv 明确"
            "139.7mm 尾管目标段缺失，目标窗仅为模型聚焦窗（model_focus）。",
        ),
    )


def _build_fluids() -> tuple[FluidSpec, ...]:
    """呼探1 流体清单（2026-08-16 重建）。

    流变标注：钻井液 PV/YP 为现场实测；领/尾浆 n/K 为同区块同密度体系 proxy；
    冲洗液/隔离液/后置液/替浆液流变为 proxy（现场均无六速实测）。
    """
    return (
        FluidSpec(
            name="钻井液",
            role=FluidRole.MUD,
            density_kg_m3=HU1_MUD_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU1_MUD_PV_PA_S,
            yield_stress_pa=HU1_MUD_YP_PA,
        ),
        FluidSpec(
            name="冲洗液",
            role=FluidRole.WASH,
            density_kg_m3=HU1_WASH_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU1_WASH_PV_PA_S,
            yield_stress_pa=HU1_WASH_YP_PA,
        ),
        FluidSpec(
            name="加重隔离液",
            role=FluidRole.SPACER,
            density_kg_m3=HU1_SPACER_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU1_SPACER_PV_PA_S,
            yield_stress_pa=HU1_SPACER_YP_PA,
        ),
        FluidSpec(
            name="领浆",
            role=FluidRole.LEAD,
            density_kg_m3=HU1_LEAD_DENSITY_KG_M3,
            rheology_model=RheologyModel.POWER_LAW,
            power_law_n=HU1_LEAD_POWER_LAW_N,
            consistency_k=HU1_LEAD_CONSISTENCY_K,
        ),
        FluidSpec(
            name="尾浆",
            role=FluidRole.TAIL,
            density_kg_m3=HU1_TAIL_DENSITY_KG_M3,
            rheology_model=RheologyModel.POWER_LAW,
            power_law_n=HU1_TAIL_POWER_LAW_N,
            consistency_k=HU1_TAIL_CONSISTENCY_K,
        ),
        FluidSpec(
            name="后置液",
            role=FluidRole.OTHER,
            density_kg_m3=HU1_PLUG_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU1_PLUG_PV_PA_S,
            yield_stress_pa=HU1_PLUG_YP_PA,
        ),
        FluidSpec(
            name="替浆液1",
            role=FluidRole.DISPLACEMENT,
            density_kg_m3=HU1_DISPLACEMENT1_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU1_DISPLACEMENT_PV_PA_S,
            yield_stress_pa=HU1_DISPLACEMENT_YP_PA,
        ),
        FluidSpec(
            name="替浆液2",
            role=FluidRole.DISPLACEMENT,
            density_kg_m3=HU1_DISPLACEMENT2_DENSITY_KG_M3,
            rheology_model=RheologyModel.BINGHAM,
            plastic_viscosity_pa_s=HU1_DISPLACEMENT_PV_PA_S,
            yield_stress_pa=HU1_DISPLACEMENT_YP_PA,
        ),
    )


def _build_schedule() -> PumpingSchedule:
    """呼探1 现场施工泵注程序（2026-08-16 重建，来源 204151.doc 施工总结）。

    施工顺序：冲洗液15 → 加重隔离液15 → 领浆24 → 尾浆36 → 后置液2（替浆首段）→
    替浆液1（管串）36 → 替浆液2（环空）49；替浆总量 87m³。
    """
    return PumpingSchedule(
        steps=(
            PumpingScheduleStep(
                step_name="注入冲洗液",
                fluid_name="冲洗液",
                volume_m3=HU1_WASH_VOLUME_M3,
                rate_m3_min=HU1_FRONT_RATE_M3_MIN,
                remarks="现场 204151.doc: 冲洗型隔离液 15m³，密度2.00g/cm³，排量0.65m³/min，泵压6.8-7.6MPa。",
            ),
            PumpingScheduleStep(
                step_name="注入加重隔离液",
                fluid_name="加重隔离液",
                volume_m3=HU1_SPACER_VOLUME_M3,
                rate_m3_min=HU1_FRONT_RATE_M3_MIN,
                remarks="现场 204151.doc: 加重隔离液 15m³，密度2.05g/cm³，排量0.65m³/min，泵压6.5-8.0MPa。",
            ),
            PumpingScheduleStep(
                step_name="注入领浆",
                fluid_name="领浆",
                volume_m3=HU1_LEAD_VOLUME_M3,
                rate_m3_min=HU1_FRONT_RATE_M3_MIN,
                remarks="现场 204151.doc: 领浆 24m³（实际），密度2.10g/cm³，排量0.65m³/min，泵压8-8.5MPa，封固3523-5300m。",
            ),
            PumpingScheduleStep(
                step_name="注入尾浆",
                fluid_name="尾浆",
                volume_m3=HU1_TAIL_VOLUME_M3,
                rate_m3_min=HU1_FRONT_RATE_M3_MIN,
                remarks="现场 204151.doc: 尾浆 36m³（14+22），密度1.90g/cm³，排量0.65m³/min，泵压9-13MPa，封固5300-7601m。",
            ),
            PumpingScheduleStep(
                step_name="注入后置液（替浆首段）",
                fluid_name="后置液",
                volume_m3=HU1_PLUG_VOLUME_M3,
                rate_m3_min=HU1_PLUG_RATE_M3_MIN,
                remarks="现场 204151.doc: 后置液 2m³，密度2.05g/cm³（含于替浆87m³，管内占位，碰压前段）。",
            ),
            PumpingScheduleStep(
                step_name="替浆液1（管串）",
                fluid_name="替浆液1",
                volume_m3=HU1_DISPLACEMENT1_VOLUME_M3,
                rate_m3_min=HU1_DISPLACEMENT1_RATE_M3_MIN,
                remarks="现场 204151.doc: 替浆液1（管串）36m³，密度1.95g/cm³，排量0.65→0.5（取0.5）。",
            ),
            PumpingScheduleStep(
                step_name="替浆液2（环空）",
                fluid_name="替浆液2",
                volume_m3=HU1_DISPLACEMENT2_VOLUME_M3,
                rate_m3_min=HU1_DISPLACEMENT2_RATE_M3_MIN,
                remarks="现场 204151.doc: 替浆液2（环空）49m³，密度2.08g/cm³，排量0.36-0.42（取0.36），碰压6MPa。",
            ),
        ),
        notes=(
            "现场施工顺序（204151.doc 2020-10-25）：冲洗液→加重隔离液→领浆→尾浆→投钻杆胶塞→替浆87m³"
            "（后置2 + 管串36 + 环空49）。",
            "替浆总量 87m³ 为现场口径（pumping_schedule.csv）；排量全程 0.65→0.5→0.36→0.42 m³/min 递减"
            "（约 240min），替浆液1/2 排量取代表性值。",
            "2026-08-16 重建：旧 HT1-001 误复制泵注（平衡液40/隔离液20/领浆20.6/中间浆28.7/尾浆22.1/替浆93.7）全部废弃。",
            "管内容积口径：尾管单值 ID 111.16mm 估算；提取包缺悬挂器以上钻杆尺寸，未分段 shoe_lag_volume，"
            "1D 到鞋时间按 liner_id_mm 单值口径近似（与 hu102 一致）。",
        ),
    )


def load_hu1_tailpipe(
    *,
    reference_root: Path | None = None,
    caliper_csv_path: Path | None = None,
    inclination_csv_path: Path | None = None,
) -> tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载呼探1井（hu1）139.7mm 尾管段"现场实际施工版"模型输入（默认）。

    2026-08-16 整体重建：旧版为 HT1-001 误复制件（TD/hanger/双径尾管/泵注均为错井值），
    旧版驱动的所有模拟结果作废。本函数按 hu1_呼探1 提取包重建。

    Args:
        reference_root: 可选参考资料根目录（默认参考文档/呼探1，旧数据包追溯）。
        caliper_csv_path: 可选现场提取包井径 CSV 路径（默认 hu1_呼探1/caliper_profile.csv）。
        inclination_csv_path: 可选现场提取包井斜 CSV 路径（默认 hu1_呼探1/inclination_profile.csv）。

    Returns:
        (well_spec, fluids, schedule, validation_data)
    """
    resolved_reference_root = reference_root or DEFAULT_REFERENCE_ROOT
    resolved_caliper_csv_path = caliper_csv_path or DEFAULT_CALIPER_CSV
    resolved_inclination_csv_path = inclination_csv_path or DEFAULT_INCLINATION_CSV

    well_spec = _build_well_spec(resolved_reference_root, resolved_caliper_csv_path, resolved_inclination_csv_path)
    fluids = _build_fluids()
    schedule = _build_schedule()

    validation_data = ValidationData(
        cbl_summary_path=DEFAULT_CBL_EVALUATION_CSV,
        cbl_pass_rate=None,
        job_report_path=resolved_reference_root / "提取数据" / "呼探1井_固井顶替模型数据包.json",
        notes=(
            "CBL 仅定性（3 张 CBL/VDL 图，cbl_evaluation.csv）：100611.jpg 7210-7560m 良好 + 7415-7430m 差；"
            "100612.jpg 7333-7607m 中 7474-7514m（40m）不合格、7514-7586m 中等；100613.jpg 5630-7425m 90% 以上优良。"
            "整体'合格'（完井报告 10068.doc）。现场无数值合格率，cbl_pass_rate=None，不编造。",
            "CBL 电测测量段 3530-7601m（construction_events.csv）；目标层段 target_intervals.csv 标 missing。",
            "2026-08-16 重建：几何/泵注/流体全部取自 hu1_呼探1 提取包；密度 field_measured、"
            "水泥浆 n/K 与前置液/替浆液流变为 proxy（现场无六速实测）。",
        ),
    )
    return well_spec, fluids, schedule, validation_data


def load_hu1_tailpipe_actual(
    *,
    reference_root: Path | None = None,
    caliper_csv_path: Path | None = None,
    inclination_csv_path: Path | None = None,
) -> tuple[WellSpec, tuple[FluidSpec, ...], PumpingSchedule, ValidationData]:
    """加载呼探1井尾管段"现场实际施工版"模型输入（与 load_hu1_tailpipe 同义）。

    呼探1现场仅有一套实际施工记录（204151.doc），无独立设计版；本函数为命名对称提供
    的统一实际版入口（参照 hu103/ht1_004 的 load_<well>_tailpipe_actual 模式）。
    """
    return load_hu1_tailpipe(
        reference_root=reference_root,
        caliper_csv_path=caliper_csv_path,
        inclination_csv_path=inclination_csv_path,
    )


def export_hu1_sync_card_markdown(output_dir: Path) -> Path:
    """导出 呼探1井 同步画像卡为 Markdown 文件。"""

    well_spec, fluids, schedule, _ = load_hu1_tailpipe()
    solver = CasingFlowSolver()
    result = solver.run(well_spec, fluids, schedule)
    provenance = WELL_PROVENANCE.get(HU1_WELL_NAME)
    if provenance is None:
        from cemdisp.data.provenance import _WELL_NAME_ALIASES
        canonical = _WELL_NAME_ALIASES.get(HU1_WELL_NAME, HU1_WELL_NAME)
        provenance = WELL_PROVENANCE[canonical]
    sync_card = build_sync_card(HU1_WELL_NAME, result.shoe_timeline, provenance)

    output_path = output_dir / "呼探1井_同步画像卡.md"
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