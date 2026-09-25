# -*- coding: utf-8 -*-
"""调用签名闸门：形参被删却仍被传 ⇒ 运行即 `TypeError`，编译期无人察觉。

本闸门用 AST 全仓扫描对 `AnnulusD2DGASolver` / `CasingFlowSolver` 的调用关键字，
与 `inspect.signature` 的真实形参比对，把"静默腐烂"提前到测试期（Task 5）。
"""
from pathlib import Path

from scripts.entrypoints.check_call_signatures import find_bad_kwargs

ROOT = Path(__file__).resolve().parents[2]


def test_no_bad_keyword_arguments_anywhere():
    """全仓生产调用点冻结：不存在任何"未知形参"调用。"""
    bad = find_bad_kwargs(ROOT)
    assert bad == [], "存在未知形参调用点：" + "; ".join(
        f"{f}:{n}:{c}.{k}" for f, n, c, k in bad)


# --------------------------------------------------------------------------- #
# 有牙证明：样例文件里的坏调用必须被抓到（防"扫描器恒返回空"的假绿）
# --------------------------------------------------------------------------- #

def test_scanner_detects_a_known_bad_call(tmp_path):
    p = tmp_path / "sample.py"
    p.write_text("from cemdisp.models2d import AnnulusD2DGASolver\n"
                 "s = AnnulusD2DGASolver(enable_d2dga_auto_m=True)\n", encoding="utf-8")
    assert "enable_d2dga_auto_m" in {k for *_, k in find_bad_kwargs(tmp_path)}


def test_scanner_detects_bad_call_on_casing_solver(tmp_path):
    """两个目标类都要扫——`CasingFlowSolver` 同样不能漏。"""
    p = tmp_path / "casing.py"
    p.write_text("from cemdisp.transport1d import CasingFlowSolver\n"
                 "s = CasingFlowSolver(enable_graviy=True)\n", encoding="utf-8")
    assert "enable_graviy" in {k for *_, k in find_bad_kwargs(tmp_path)}


def test_scanner_accepts_real_parameters(tmp_path):
    """不误报：真实存在的形参（含布尔开关、显式偏离默认的组合）不得被标记。"""
    p = tmp_path / "ok.py"
    p.write_text(
        "from cemdisp.models2d import AnnulusD2DGASolver\n"
        "from cemdisp.transport1d import CasingFlowSolver\n"
        "a = AnnulusD2DGASolver(nz=10, ny=5, total_t=20.0, enable_d2dga=True,\n"
        "                      enable_d2dga_i3_flux=True, enable_true_buoyancy=False,\n"
        "                      open_outlet=True, enable_yield_gate=False)\n"
        "b = CasingFlowSolver(enable_gravity=True)\n", encoding="utf-8")
    assert find_bad_kwargs(tmp_path) == []


# --------------------------------------------------------------------------- #
# 豁免面：`with pytest.raises(...)` 内的**故意反例**不算坏调用
# （`tests/contract/test_no_invented_dispersion.py` 用它对已删形参做
# "传值必须显式失败"契约；豁免面必须窄到只剩这一种形态）
# --------------------------------------------------------------------------- #

def test_scanner_skips_deliberate_raises_negative(tmp_path):
    p = tmp_path / "negative.py"
    p.write_text("import pytest\n"
                 "from cemdisp.models2d import AnnulusD2DGASolver\n"
                 "def test_x():\n"
                 "    with pytest.raises(TypeError):\n"
                 "        AnnulusD2DGASolver(dispersion_dt_scale=1.0)\n", encoding="utf-8")
    assert find_bad_kwargs(tmp_path) == []


def test_scanner_still_flags_same_file_outside_raises(tmp_path):
    """豁免面窄：同一文件里 `raises` 块**之外**的坏调用照样被抓。"""
    p = tmp_path / "mixed.py"
    p.write_text("import pytest\n"
                 "from cemdisp.models2d import AnnulusD2DGASolver\n"
                 "def test_x():\n"
                 "    with pytest.raises(TypeError):\n"
                 "        AnnulusD2DGASolver(dispersion_dt_scale=1.0)\n"
                 "def test_y():\n"
                 "    AnnulusD2DGASolver(enable_d2dga_auto_m=True)\n", encoding="utf-8")
    assert {k for *_, k in find_bad_kwargs(tmp_path)} == {"enable_d2dga_auto_m"}
