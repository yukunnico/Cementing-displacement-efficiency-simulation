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


def test_deprecation_warning_on_non_default_dispersion_args():
    """显式传入任一 dispersion_* 形参时触发 DeprecationWarning(中文消息)。"""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        m.AnnulusD2DGASolver(nz=10, ny=5, total_t=20.0, dispersion_dt_scale=1.0)
    dep = [w for w in caught if issubclass(w.category, DeprecationWarning)]
    assert dep, "显式传 dispersion_* 应触发 DeprecationWarning"
    assert any("弥散" in str(w.message) for w in dep), "弃用消息应为中文且说明弥散已删除"


@pytest.mark.parametrize("kwargs", [
    {},  # 全默认(形参默认 None)
    {"dispersion_axial": None, "dispersion_azimuthal": None,
     "dispersion_dt_ref": None, "dispersion_dt_scale": None},  # 显式 None 视同默认
])
def test_no_deprecation_warning_when_args_absent(kwargs):
    """默认构造或显式传 None:不触发 DeprecationWarning(8 个 runner 无感)。"""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        m.AnnulusD2DGASolver(nz=10, ny=5, total_t=20.0, **kwargs)
    dep = [w for w in caught if issubclass(w.category, DeprecationWarning)]
    assert not dep, f"不应触发弃用警告: {[str(w.message) for w in dep]}"
