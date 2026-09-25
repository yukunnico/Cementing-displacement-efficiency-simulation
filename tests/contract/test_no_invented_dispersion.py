"""契约测试:环空 2D 求解器不得含自创拉普拉斯弥散。

文献依据:Zhang & Frigaard (2022, JFM 947:A32) p.11
``we have no diffusive terms``——原模型无任何人工扩散项,
弥散由 q₀(平均平流通量)+ I₃(浮力通量,Z&F22 式 4.25 第二项 / (4.26) I₃)
分层通量闭合承载。

2026-09-14 源模型口径重构 Task 7:删除 ``_smooth_dispersion`` 及其
run 循环中的三处调用(lead/tail/spacer),并连带删除原设计中
"弥散 → 再次执行四相过填修正"的修补步(弥散删除后该修正与平流后的
过填修正之间无任何浓度场修改,成为同一步内的重复执行)。
"""
from __future__ import annotations

import inspect
import warnings

import pytest

import cemdisp.models2d.annulus_d2dga as m


def test_no_laplacian_dispersion_in_concentration_update():
    """run() 浓度更新中不得出现自创拉普拉斯平滑(Z&F22 p.11: no diffusive terms)。"""
    src = inspect.getsource(m.AnnulusD2DGASolver.run)
    assert "_smooth_dispersion" not in src


def test_smooth_dispersion_removed_from_module():
    """``_smooth_dispersion`` 定义与全部调用点已从模块中删除。"""
    src = inspect.getsource(m)
    assert "_smooth_dispersion" not in src


def test_dispersion_params_fully_removed():
    """⚠️ 2026-09-18 收紧：4 个 dispersion_* 形参已从求解器**彻底删除**。

    原契约是"形参保留 + 传值触发 DeprecationWarning"（2026-09-14 Task 7 的过渡态）。
    自创拉普拉斯弥散既已删除，保留形参只是把死开关留在 API 面上；现改为直接删除——
    任何传值都是 ``TypeError``，比"静默忽略 + 告警"更难误用。
    """
    import inspect
    params = inspect.signature(m.AnnulusD2DGASolver.__init__).parameters
    for name in ("dispersion_axial", "dispersion_azimuthal",
                 "dispersion_dt_ref", "dispersion_dt_scale"):
        assert name not in params, f"{name} 应已从 __init__ 移除"
    # 属性面同样不得残留（旧版保留 self.dispersion_* 供脚本读值）
    solver = m.AnnulusD2DGASolver(nz=10, ny=5, total_t=20.0)
    for name in ("dispersion_axial", "dispersion_azimuthal",
                 "dispersion_dt_ref", "dispersion_dt_scale"):
        assert not hasattr(solver, name), f"self.{name} 应已移除"


def test_removed_dispersion_params_raise_type_error():
    """传值必须显式失败（不得再被静默吞掉）。"""
    with pytest.raises(TypeError):
        m.AnnulusD2DGASolver(nz=10, ny=5, total_t=20.0, dispersion_dt_scale=1.0)


def test_no_deprecation_warning_on_default_construction():
    """默认构造不触发任何 DeprecationWarning（8 个 runner 无感）。"""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        m.AnnulusD2DGASolver(nz=10, ny=5, total_t=20.0)
    dep = [w for w in caught if issubclass(w.category, DeprecationWarning)]
    assert not dep, f"不应触发弃用警告: {[str(w.message) for w in dep]}"
