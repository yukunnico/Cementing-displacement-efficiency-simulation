"""
逐井数据加载器公共基础模块

本模块提取各井 loader 中重复出现的公共逻辑，包括：
- CSV 剖面文件读取（井径/井斜）
- 等效井径换算（双径尾管面积守恒）
- 管内容积估算（鞋口滞后体积）
- 深度-数值剖面点构建
- WellSpec / FluidSpec / PumpingSchedule 统一组装

数据来源标注约定：
- [实测]：直接来自现场实测数据或正式报告
- [代理]：来自邻井或类似工况的近似值
- [估算]：基于工程经验或公式推算的值
"""

from __future__ import annotations

import csv
import math
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep
from cemdisp.data.validation_data import ValidationData
from cemdisp.data.well_spec import DepthValuePoint, EvaluationWindow, WellSpec


# ---------------------------------------------------------------------------
# CSV 剖面读取
# ---------------------------------------------------------------------------


def read_profile_csv(
    csv_path: Path,
    encoding: str = "utf-8-sig",
) -> tuple[tuple[float, float, float], ...]:
    """读取井径/井斜 CSV 文件。

    CSV 文件需包含列：depth_md_m, hole_diameter_mm, inclination_deg。
    返回按测深升序排列的 (depth_md_m, hole_diameter_mm, inclination_deg) 元组。

    Args:
        csv_path: CSV 文件路径。
        encoding: 文件编码，默认 utf-8-sig（兼容 Excel 导出的 BOM）。

    Returns:
        按测深排序的三元组元组。

    Raises:
        FileNotFoundError: 文件不存在。
        ValueError: CSV 为空或缺少必要列。
    """

    with csv_path.open("r", encoding=encoding, newline="") as handle:
        rows: list[tuple[float, float, float]] = []
        for row in csv.DictReader(handle):
            rows.append((
                float(row["depth_md_m"]),
                float(row["hole_diameter_mm"]),
                float(row["inclination_deg"]),
            ))
    if not rows:
        raise ValueError(f"井径/井斜 CSV 为空: {csv_path}")
    return tuple(sorted(rows))


# ---------------------------------------------------------------------------
# 几何计算工具
# ---------------------------------------------------------------------------


def equivalent_hole_diameter_mm(
    actual_hole_mm: float,
    actual_od_mm: float,
    reference_od_mm: float,
) -> float:
    """在固定 reference_od 条件下，构造与原始环空面积等价的井眼直径。

    面积守恒：D_eq^2 - OD_ref^2 = D_actual^2 - OD_actual^2。
    用于双径尾管上段井眼等效处理，使 2D 求解器用单一外径即可保持
    环空截面积一致。

    [估算] 此为工程等效换算，不来自实测。

    Args:
        actual_hole_mm: 实际井眼直径 (mm)。
        actual_od_mm: 实际管柱外径 (mm)。
        reference_od_mm: 参考管柱外径 (mm)，即求解器使用的统一外径。

    Returns:
        等效井眼直径 (mm)。
    """

    area_term = actual_hole_mm**2 - actual_od_mm**2 + reference_od_mm**2
    return math.sqrt(max(area_term, reference_od_mm**2))


def pipe_volume_m3(length_m: float, inner_diameter_mm: float) -> float:
    """按内径和长度估算管内容积 (m^3)。

    [估算] 用于鞋口滞后体积计算，不代表真实管柱内径。

    Args:
        length_m: 管段长度 (m)。
        inner_diameter_mm: 管内径 (mm)。

    Returns:
        管内容积 (m^3)。
    """

    radius_m = inner_diameter_mm / 1000.0 / 2.0
    return math.pi * radius_m**2 * length_m


# ---------------------------------------------------------------------------
# 剖面点构建
# ---------------------------------------------------------------------------


def build_depth_points(
    values: tuple[tuple[float, float], ...],
) -> tuple[DepthValuePoint, ...]:
    """把 (测深, 数值) 元组序列转换为 WellSpec 使用的剖面点。"""

    return tuple(DepthValuePoint(depth_md_m=depth, value=value) for depth, value in values)


def build_hole_profile(
    profile_rows: Iterable[tuple[float, float, float]],
) -> tuple[DepthValuePoint, ...]:
    """从 CSV 剖面行构建井径剖面点。"""

    return tuple(
        DepthValuePoint(depth_md_m=depth, value=hole_diameter)
        for depth, hole_diameter, _ in profile_rows
    )


def build_inclination_profile(
    profile_rows: Iterable[tuple[float, float, float]],
) -> tuple[DepthValuePoint, ...]:
    """从 CSV 剖面行构建井斜剖面点。"""

    return tuple(
        DepthValuePoint(depth_md_m=depth, value=inclination)
        for depth, _, inclination in profile_rows
    )


def build_standoff_profile(
    profile_rows: Iterable[tuple[float, float, float]],
    standoff_func: Any,
    liner_od_mm: float,
) -> tuple[DepthValuePoint, ...]:
    """从 CSV 剖面行和居中度计算函数构建居中度剖面点。

    Args:
        profile_rows: CSV 剖面行，每行为 (depth, hole_diameter, inclination)。
        standoff_func: 居中度计算回调，签名 (depth, hole_diameter, liner_od) -> float。
        liner_od_mm: 尾管外径 (mm)。

    Returns:
        居中度剖面点元组。
    """

    return tuple(
        DepthValuePoint(
            depth_md_m=depth,
            value=standoff_func(depth, hole_diameter, liner_od_mm),
        )
        for depth, hole_diameter, _ in profile_rows
    )


# ---------------------------------------------------------------------------
# 流体与施工程序组装
# ---------------------------------------------------------------------------


def build_fluids(
    fluid_definitions: tuple[dict[str, Any], ...],
) -> tuple[FluidSpec, ...]:
    """根据字典定义批量构建 FluidSpec 元组。

    每个字典需包含以下键：
    - name (str): 流体名称
    - role (FluidRole): 流体角色
    - density_kg_m3 (float): 密度
    - rheology_model (RheologyModel): 流变模型

    Bingham 模型额外需要：
    - plastic_viscosity_pa_s (float)
    - yield_stress_pa (float)

    幂律模型额外需要：
    - power_law_n (float)
    - consistency_k (float)

    Args:
        fluid_definitions: 流体定义字典元组。

    Returns:
        FluidSpec 元组。
    """

    return tuple(FluidSpec(**defn) for defn in fluid_definitions)


def build_pumping_schedule(
    steps_definitions: tuple[dict[str, Any], ...],
    notes: tuple[str, ...] = (),
) -> PumpingSchedule:
    """根据字典定义构建 PumpingSchedule。

    每个步骤字典需包含以下键：
    - step_name (str): 步骤名称
    - fluid_name (str): 流体名称
    - volume_m3 (float): 注入体积
    - rate_m3_min (float): 泵注排量

    可选键：
    - remarks (str): 备注，默认空字符串
    - event_tag (PumpingStageEvent): 作业阶段标签，默认 None

    Args:
        steps_definitions: 步骤定义字典元组。
        notes: 施工程序整体备注。

    Returns:
        PumpingSchedule 实例。
    """

    steps = tuple(PumpingScheduleStep(**defn) for defn in steps_definitions)
    return PumpingSchedule(steps=steps, notes=notes)
