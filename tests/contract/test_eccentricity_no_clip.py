# -*- coding: utf-8 -*-
"""Task 10 契约测试：e_clip 截断移除，e = 1−standoff 文献口径直取（2026-09-15）。

文献锚点：偏心度 e ∈ [0, 1)（Pelipenko04 (2.1)）。
旧口径 ``e = clip(1−standoff, 0.05, 0.55)`` 把 standoff<0.45 的全部井压进
e=0.55 死区（standoff 0.35 与 0.45 逐位同几何），并人为封顶强偏心响应；
本契约锁定新口径：

1. standoff=0.35 与 0.45 给出**不同**的 e（0.65 vs 0.55，死区消除）；
2. 端点行为：standoff→0 ⇒ e→1−1e-6（<1，开区间上界）；standoff→1 ⇒ e→1e-6；
3. ``e_clip_max``/``e_clip_measured_max``/``enable_e_clip_ruling`` 三形参保留但
   失效——偏离 legacy 默认（0.55/0.90/True）的传值触发一次性 DeprecationWarning
   （中文消息，含 Pelipenko04 (2.1) 引文），legacy 等值传参与默认传参静默；
4. 基准算例 runner 不再传 ``e_clip_max``（Task 12 起完成弃用形参退役迁移）：
   构造零弃用警告，e=0.8 算例几何按论文原值直取。
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from cemdisp.data.well_spec import DepthValuePoint, WellSpec
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
from cemdisp.runners.zhang2022_benchmark import (
    ZHANG2022_CASE_BY_ID,
    build_case_solver,
    build_case_well_spec,
)


def _make_solver(**kw) -> AnnulusD2DGASolver:
    return AnnulusD2DGASolver(dt=4.0, nz=20, ny=10, total_t=40.0, **kw)


def _toy_well(standoff: float, standoff_measured: bool = False) -> WellSpec:
    pts = lambda d, v: DepthValuePoint(depth_md_m=d, value=v)
    return WellSpec(
        well_name="toy",
        top_md_m=1000.0, bottom_md_m=1100.0, shoe_md_m=1100.0, hanger_md_m=1000.0,
        casing_id_mm=200.0, liner_od_mm=139.7, liner_id_mm=108.0,
        hole_diameter_profile=[pts(1000.0, 215.9), pts(1100.0, 215.9)],
        inclination_profile=[pts(1000.0, 5.0), pts(1100.0, 5.0)],
        standoff_profile=[pts(1000.0, standoff), pts(1100.0, standoff)],
        standoff_measured=standoff_measured,
        evaluation_windows=[],
    )


def _geom_e(standoff: float, **solver_kw) -> np.ndarray:
    """构造给定常数居中度的玩具井几何，返回 e 数组（均匀剖面 ⇒ e 逐深一致）。"""
    s = _make_solver(**solver_kw)
    return s._build_geom(_toy_well(standoff))["e"]


# ---------------------------------------------------------------------------
# 1. 死区消除：0.35 与 0.45 必须不同
# ---------------------------------------------------------------------------
class TestDeadZoneRemoved:
    def test_standoff_035_vs_045_give_different_e(self):
        """standoff=0.35 → e=0.65，standoff=0.45 → e=0.55（旧口径两者同为 0.55）。"""
        e_low_so = _geom_e(0.35)
        e_high_so = _geom_e(0.45)
        # 旧口径 clip(1−SO, 0.05, 0.55) 使两者同为 0.55 死区
        assert not np.allclose(e_low_so, e_high_so), (
            f"standoff 0.35 与 0.45 的 e 仍相同（死区未消除）：{e_low_so[0]} vs {e_high_so[0]}"
        )
        assert np.allclose(e_low_so, 0.65, atol=1e-12), f"e(0.35) 应为 0.65，得到 {e_low_so[0]}"
        assert np.allclose(e_high_so, 0.55, atol=1e-12), f"e(0.45) 应为 0.55，得到 {e_high_so[0]}"

    def test_strong_eccentricity_not_capped(self):
        """强偏心不再封顶：standoff=0.05 → e=0.95（旧口径被压到 0.55/0.90）。"""
        assert np.allclose(_geom_e(0.05), 0.95, atol=1e-12)


# ---------------------------------------------------------------------------
# 2. 端点行为：e ∈ [0, 1)（Pelipenko04 (2.1) 开区间上界）
# ---------------------------------------------------------------------------
class TestEndpointOpenInterval:
    def test_standoff_zero_gives_e_below_one(self):
        """standoff→0 ⇒ e→1−1e-6（护栏夹取，e 严格小于 1）。"""
        e = _geom_e(0.0)
        assert np.allclose(e, 1.0 - 1.0e-6, atol=1e-15), f"e 应为 1−1e-6，得到 {e[0]}"
        assert e.max() < 1.0, "e 必须严格 < 1（Pelipenko04 (2.1) 开区间）"

    def test_standoff_one_gives_e_above_zero(self):
        """standoff→1 ⇒ e→1e-6（同心极限护栏，避免零间隙除零）。"""
        e = _geom_e(1.0)
        assert np.allclose(e, 1.0e-6, atol=1e-15), f"e 应为 1e-6，得到 {e[0]}"
        assert e.min() > 0.0, "e 必须严格 > 0"

    def test_no_lower_dead_floor_at_005(self):
        """旧口径下限 0.05 已移除：standoff=0.97 → e=0.03（不再被抬到 0.05）。"""
        assert np.allclose(_geom_e(0.97), 0.03, atol=1e-12)


# ---------------------------------------------------------------------------
# 3. 弃用形参：偏离 legacy 默认才警告，等值/默认静默
# ---------------------------------------------------------------------------
class TestDeprecatedEClipParams:
    @pytest.mark.parametrize("kw", [
        {"e_clip_max": 0.90},
        {"e_clip_max": 1.0},
        {"e_clip_measured_max": 0.95},
        {"enable_e_clip_ruling": False},
        {"e_clip_max": 1.0, "e_clip_measured_max": 0.95, "enable_e_clip_ruling": False},
    ])
    def test_nondefault_params_warn_once(self, kw):
        """任一形参偏离 legacy 默认（0.55/0.90/True）→ 恰一次 DeprecationWarning。"""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            _make_solver(**kw)
        dep = [w for w in caught if issubclass(w.category, DeprecationWarning)]
        assert len(dep) == 1, f"应恰一次 DeprecationWarning，得到 {len(dep)}：{[str(w.message) for w in dep]}"
        message = str(dep[0].message)
        assert "Pelipenko04" in message and "2.1" in message, (
            f"弃用消息应含 Pelipenko04 (2.1) 文献引文：{message}"
        )

    @pytest.mark.parametrize("kw", [
        {},
        {"e_clip_max": 0.55},
        {"e_clip_measured_max": 0.90},
        {"enable_e_clip_ruling": True},
        {"e_clip_max": 0.55, "e_clip_measured_max": 0.90, "enable_e_clip_ruling": True},
    ])
    def test_default_and_legacy_equal_params_silent(self, kw):
        """默认与 legacy 等值传参（0.55/0.90/True）完全静默（调用方无感迁移）。"""
        with warnings.catch_warnings():
            warnings.simplefilter("error")  # 任何警告都视为异常
            _make_solver(**kw)

    def test_params_kept_but_inert(self):
        """三形参属性仍保留（旧脚本读值兼容），但不再被 _build_geom 消费。"""
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            s = _make_solver(e_clip_max=0.90, e_clip_measured_max=0.80, enable_e_clip_ruling=False)
        assert s.e_clip_max == 0.90
        assert s.e_clip_measured_max == 0.80
        assert s.enable_e_clip_ruling is False
        # 形参失效：standoff=0.05 的井 e=0.95，不被 0.90/0.80 上限截断
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            e_uncapped = _geom_e(0.05, e_clip_max=0.90)
        assert np.allclose(e_uncapped, 0.95, atol=1e-12)
        # enable_e_clip_ruling + standoff_measured 的旧裁定组合同样失效
        s2 = _make_solver()
        geom = s2._build_geom(_toy_well(0.05, standoff_measured=True))
        assert np.allclose(geom["e"], 0.95, atol=1e-12), (
            f"实测井不应再走 0.90 裁定上限，e 应为 0.95，得到 {geom['e'][0]}"
        )


# ---------------------------------------------------------------------------
# 4. 基准算例 runner 路径：不再传弃用形参 e_clip_max（Task 12 迁移完成）
# ---------------------------------------------------------------------------
class TestBenchmarkRunnerPath:
    def test_runner_no_longer_passes_e_clip_max(self):
        """runner 停止传 e_clip_max：零弃用警告且 e=0.8 几何正确。

        （Task 10 迁移期契约"显式传 1.0 收警告仍工作"已随 Task 12 runner 停传
        而升级为"零警告"口径；几何不变性断言保持。）
        """
        case = ZHANG2022_CASE_BY_ID[1]  # e=0.8 算例
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            solver = build_case_solver(case, nz=20)
        dep = [w for w in caught if issubclass(w.category, DeprecationWarning)]
        assert not dep, f"runner 不应再触发 e_clip 弃用警告：{[str(w.message) for w in dep]}"
        ws = build_case_well_spec(case)
        e = solver._build_geom(ws)["e"][0]
        assert np.isclose(e, 0.8, atol=1e-12), f"e=0.8 算例几何应保持 0.8，得到 {e}"
