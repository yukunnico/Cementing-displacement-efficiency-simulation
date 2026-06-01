"""
井筒与评价井段数据结构

本模块定义了描述单口井井筒几何参数和评价井段的核心数据结构。

主要类：
- DepthValuePoint: 测深-数值剖面上的单个数据点（井径、井斜、偏心度等）
- EvaluationWindow: 需要单独评价的井段窗口（如CBL评价井段、目标层段等）
- WellSpec: 单口井的完整输入规格，包含井段范围、套管参数、剖面数据等

数据约束：
- 所有深度值必须为正数且有限
- 井段顶部深度必须小于底部深度
- 鞋口深度必须位于井段范围内
- 评价窗口不得超出井段范围
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Tuple
import math


def _require_finite(name: str, value: float) -> None:
    """校验数值是否为有限值的内部辅助函数。"""
    if not math.isfinite(value):
        raise ValueError(f"{name}必须为有限数值")


def _require_positive(name: str, value: float) -> None:
    """校验数值是否大于零的内部辅助函数。"""
    _require_finite(name, value)
    if value <= 0.0:
        raise ValueError(f"{name}必须大于0")


def _check_dual_diameter_fields(spec: "WellSpec") -> None:
    """校验双径向井上段字段是否满足全有或全无约束。"""
    any_present = (
        spec.upper_section_bottom_md_m is not None
        or spec.upper_liner_od_mm is not None
        or spec.upper_liner_id_mm is not None
    )
    all_present = (
        spec.upper_section_bottom_md_m is not None
        and spec.upper_liner_od_mm is not None
        and spec.upper_liner_id_mm is not None
    )
    if any_present and not all_present:
        raise ValueError(
            "upper_section_bottom_md_m、upper_liner_od_mm、upper_liner_id_mm "
            "须全部提供或全部省略（全有或全无约束）"
        )


@dataclass(frozen=True)
class DepthValuePoint:
    """测深-数值剖面点。"""

    depth_md_m: float
    value: float

    def __post_init__(self) -> None:
        _require_positive("depth_md_m", self.depth_md_m)
        _require_finite("value", self.value)


@dataclass(frozen=True)
class EvaluationWindow:
    """需要单独评价的井段窗口。"""

    name: str
    top_md_m: float
    bottom_md_m: float
    window_type: str = "custom"

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("评价窗口名称不能为空")
        _require_positive("top_md_m", self.top_md_m)
        _require_positive("bottom_md_m", self.bottom_md_m)
        if self.top_md_m >= self.bottom_md_m:
            raise ValueError("评价窗口 top_md_m 必须小于 bottom_md_m")


@dataclass(frozen=True)
class WellSpec:
    """单口井的标准输入规格。

    支持双径向井（dual-diameter）可选字段，用于描述上段套管/衬管参数。
    上段字段（upper_section_bottom_md_m、upper_liner_od_mm、upper_liner_id_mm）
    须同时提供或同时省略，违反此约束将触发 ValueError。

    属性：
        is_dual_diameter: 上段字段全部提供时返回 True，否则返回 False
    """

    well_name: str
    top_md_m: float
    bottom_md_m: float
    shoe_md_m: float
    hanger_md_m: Optional[float] = None
    casing_id_mm: Optional[float] = None
    liner_od_mm: Optional[float] = None
    liner_id_mm: Optional[float] = None
    # --- 剖面数据 ---
    hole_diameter_profile: Tuple[DepthValuePoint, ...] = field(default_factory=tuple)
    liner_od_profile: Tuple[DepthValuePoint, ...] = field(default_factory=tuple)
    pipe_id_profile: Tuple[DepthValuePoint, ...] = field(default_factory=tuple)
    inclination_profile: Tuple[DepthValuePoint, ...] = field(default_factory=tuple)
    standoff_profile: Tuple[DepthValuePoint, ...] = field(default_factory=tuple)
    evaluation_windows: Tuple[EvaluationWindow, ...] = field(default_factory=tuple)
    reference_root: Optional[Path] = None
    notes: Tuple[str, ...] = field(default_factory=tuple)
    # --- 双径向井上段标识字段（全有或全无） ---
    upper_section_bottom_md_m: Optional[float] = None  # 上段底部深度（测深）
    upper_liner_od_mm: Optional[float] = None         # 上段衬管外径
    upper_liner_id_mm: Optional[float] = None         # 上段衬管内径
    # --- 鞋口延迟体积 ---
    shoe_lag_volume_m3: Optional[float] = None        # 鞋口处迟到体积（m³）
    # --- 衬管壁厚 ---
    liner_wall_thickness_mm: Optional[float] = None   # 衬管壁厚（mm）
    # --- 管内径剖面（1D前沿追踪用，深度→内径mm） ---
    pipe_id_profile: Tuple[DepthValuePoint, ...] = field(default_factory=tuple)

    @property
    def is_dual_diameter(self) -> bool:
        """上段字段（upper_section_bottom_md_m / upper_liner_od_mm / upper_liner_id_mm）全部提供时返回 True。"""
        return (
            self.upper_section_bottom_md_m is not None
            and self.upper_liner_od_mm is not None
            and self.upper_liner_id_mm is not None
        )

    def __post_init__(self) -> None:
        if not self.well_name.strip():
            raise ValueError("well_name不能为空")
        _require_positive("top_md_m", self.top_md_m)
        _require_positive("bottom_md_m", self.bottom_md_m)
        _require_positive("shoe_md_m", self.shoe_md_m)
        if self.top_md_m >= self.bottom_md_m:
            raise ValueError("top_md_m 必须小于 bottom_md_m")
        if not (self.top_md_m <= self.shoe_md_m <= self.bottom_md_m):
            raise ValueError("shoe_md_m 必须位于井段范围内")
        if self.hanger_md_m is not None:
            _require_positive("hanger_md_m", self.hanger_md_m)
            if not (self.top_md_m <= self.hanger_md_m <= self.bottom_md_m):
                raise ValueError("hanger_md_m 必须位于井段范围内")
        for name, value in (
            ("casing_id_mm", self.casing_id_mm),
            ("liner_od_mm", self.liner_od_mm),
            ("liner_id_mm", self.liner_id_mm),
        ):
            if value is not None:
                _require_positive(name, value)
        if self.reference_root is not None and not isinstance(self.reference_root, Path):
            raise TypeError("reference_root 必须为 pathlib.Path 或 None")
        for window in self.evaluation_windows:
            if window.top_md_m < self.top_md_m or window.bottom_md_m > self.bottom_md_m:
                raise ValueError(f"评价窗口 {window.name} 超出井段范围")
        # --- 双径向井上段字段：全有或全无约束 ---
        _check_dual_diameter_fields(self)
        # --- 上段底部深度须位于井段范围内 ---
        if self.upper_section_bottom_md_m is not None:
            _require_positive("upper_section_bottom_md_m", self.upper_section_bottom_md_m)
            if not (self.top_md_m < self.upper_section_bottom_md_m < self.bottom_md_m):
                raise ValueError("upper_section_bottom_md_m 必须在井段范围内")
        # --- 上段衬管尺寸须为正数 ---
        for name, value in (
            ("upper_liner_od_mm", self.upper_liner_od_mm),
            ("upper_liner_id_mm", self.upper_liner_id_mm),
        ):
            if value is not None:
                _require_positive(name, value)
        # --- shoe_lag_volume_m3 须为正有限值 ---
        if self.shoe_lag_volume_m3 is not None:
            _require_positive("shoe_lag_volume_m3", self.shoe_lag_volume_m3)
        # --- liner_wall_thickness_mm 须为正数 ---
        if self.liner_wall_thickness_mm is not None:
            _require_positive("liner_wall_thickness_mm", self.liner_wall_thickness_mm)