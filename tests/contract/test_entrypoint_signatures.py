# -*- coding: utf-8 -*-
"""调用签名闸门：形参被删却仍被传 ⇒ 运行即 `TypeError`，编译期无人察觉。

本闸门用 AST 全仓扫描对 `AnnulusD2DGASolver` / `CasingFlowSolver` 的调用关键字，
与 `inspect.signature` 的真实形参比对，把"静默腐烂"提前到测试期（Task 5）。

Task 14 补两件事：
  1. **dict→splat**（同文件 / 跨文件形参中转）的键集校验 —— 否则 `dict(...)+**kw`
     形态会让闸门报绿而脚本运行即崩；
  2. **豁免面收紧** —— 必须"`raises` 上下文 ∧ 首个位置参数是 `TypeError` ∧ 文件在
     `tests/` 下"三条同时成立，非测试文件自定义的 `def raises(e)` 不豁免。
"""
from pathlib import Path

from scripts.entrypoints import check_call_signatures as gate
from scripts.entrypoints.check_call_signatures import find_bad_kwargs

ROOT = Path(__file__).resolve().parents[2]


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def test_no_bad_keyword_arguments_anywhere():
    """全仓生产调用点冻结：不存在任何"未知形参"调用（含 dict→splat 形态）。"""
    bad = find_bad_kwargs(ROOT)
    assert bad == [], "存在未知形参调用点：" + "; ".join(
        f"{f}:{n}:{c}.{k}" for f, n, c, k in bad)


# --------------------------------------------------------------------------- #
# 有牙证明：样例文件里的坏调用必须被抓到（防"扫描器恒返回空"的假绿）
# --------------------------------------------------------------------------- #

def test_scanner_detects_a_known_bad_call(tmp_path):
    _write(tmp_path / "sample.py",
           "from cemdisp.models2d import AnnulusD2DGASolver\n"
           "s = AnnulusD2DGASolver(enable_d2dga_auto_m=True)\n")
    assert "enable_d2dga_auto_m" in {k for *_, k in find_bad_kwargs(tmp_path)}


def test_scanner_detects_bad_call_on_casing_solver(tmp_path):
    """两个目标类都要扫——`CasingFlowSolver` 同样不能漏。"""
    _write(tmp_path / "casing.py",
           "from cemdisp.transport1d import CasingFlowSolver\n"
           "s = CasingFlowSolver(enable_graviy=True)\n")
    assert "enable_graviy" in {k for *_, k in find_bad_kwargs(tmp_path)}


def test_scanner_accepts_real_parameters(tmp_path):
    """不误报：真实存在的形参（含布尔开关、显式偏离默认的组合）不得被标记。"""
    _write(tmp_path / "ok.py",
           "from cemdisp.models2d import AnnulusD2DGASolver\n"
           "from cemdisp.transport1d import CasingFlowSolver\n"
           "a = AnnulusD2DGASolver(nz=10, ny=5, total_t=20.0, enable_d2dga=True,\n"
           "                      enable_d2dga_i3_flux=True, enable_true_buoyancy=False,\n"
           "                      open_outlet=True, enable_yield_gate=False)\n"
           "b = CasingFlowSolver(enable_gravity=True)\n")
    assert find_bad_kwargs(tmp_path) == []


# --------------------------------------------------------------------------- #
# Task 14 Step 3/4：dict→splat 形态（原闸门的盲区，"绿而无用"的根源）
# --------------------------------------------------------------------------- #

def test_scanner_detects_bad_key_behind_dict_splat(tmp_path):
    """`kw = dict(...)` + `**kw`：坏键藏在展开里也必须被抓。"""
    _write(tmp_path / "sample.py",
           "from cemdisp.models2d import AnnulusD2DGASolver\n"
           "kw = dict(nz=10, ny=5, total_t=20.0, dispersion_axial=0.0)\n"
           "s = AnnulusD2DGASolver(**kw)\n")
    assert {k for *_, k in find_bad_kwargs(tmp_path)} == {"dispersion_axial"}


def test_scanner_detects_bad_key_after_dict_update(tmp_path):
    """`kw.update(...)` 追加的坏键同样要被抓（并集语义）。"""
    _write(tmp_path / "sample.py",
           "from cemdisp.models2d import AnnulusD2DGASolver\n"
           "kw = dict(nz=10, ny=5, total_t=20.0)\n"
           "kw.update(c_min=0.0)\n"
           "s = AnnulusD2DGASolver(**kw)\n")
    assert {k for *_, k in find_bad_kwargs(tmp_path)} == {"c_min"}


def test_scanner_detects_bad_key_in_cross_file_dict_splat(tmp_path):
    """dict 在**另一个文件**里构造、经函数形参中转后再展开 —— 也必须被抓（R113）。

    同文件 dict 字面量分析抓不到这一形态，这也是"至少 4 例 splat 落在另一文件"的修法。
    """
    _write(tmp_path / "runner.py",
           "from cemdisp.models2d import AnnulusD2DGASolver\n"
           "def run_case(well, *, solver_kw=None):\n"
           "    kw = dict(total_t=1.0, nz=10)\n"
           "    if solver_kw:\n"
           "        kw.update(solver_kw)\n"
           "    return AnnulusD2DGASolver(**kw)\n")
    _write(tmp_path / "main.py",
           "from runner import run_case\n"
           "run_case('w', solver_kw={'dispersion_azimuthal': 0.015})\n")
    assert {k for *_, k in find_bad_kwargs(tmp_path)} == {"dispersion_azimuthal"}


def test_scanner_detects_bad_key_in_variant_table_splat(tmp_path):
    """模块级变体表 + `for name, kw in VARIANTS.items(): Solver(**kw)` 也要覆盖。"""
    _write(tmp_path / "variants.py",
           "from cemdisp.models2d import AnnulusD2DGASolver\n"
           "VARIANTS = {'base': {}, 'zero': {'dispersion_axial': 0.0}}\n"
           "for name, kw in VARIANTS.items():\n"
           "    AnnulusD2DGASolver(nz=10, ny=5, total_t=20.0, **kw)\n")
    assert {k for *_, k in find_bad_kwargs(tmp_path)} == {"dispersion_axial"}


def test_scanner_accepts_clean_dict_splat(tmp_path):
    """不误报：dict 里全是真实形参时不得报错（防"恒返回非空"的假红）。"""
    _write(tmp_path / "ok.py",
           "from cemdisp.models2d import AnnulusD2DGASolver\n"
           "kw = dict(nz=10, ny=5, total_t=20.0)\n"
           "kw.update(enable_yield_gate=False)\n"
           "s = AnnulusD2DGASolver(**kw)\n")
    assert find_bad_kwargs(tmp_path) == []


def test_scanner_does_not_confuse_same_named_functions_across_files(tmp_path):
    """同名函数跨文件不得串味：A 文件的坏键不能算到 B 文件的调用点上（防假报）。"""
    _write(tmp_path / "helper_a.py",
           "from cemdisp.models2d import AnnulusD2DGASolver\n"
           "def run_case(*, solver_kw=None):\n"
           "    return AnnulusD2DGASolver(nz=10, ny=5, total_t=20.0, **(solver_kw or {}))\n")
    _write(tmp_path / "user_a.py",
           "from helper_a import run_case\n"
           "run_case(solver_kw={'dispersion_axial': 0.0})\n")
    _write(tmp_path / "helper_b.py",
           "from cemdisp.models2d import AnnulusD2DGASolver\n"
           "def run_case(*, solver_kw=None):\n"
           "    return AnnulusD2DGASolver(nz=10, ny=5, total_t=20.0, **(solver_kw or {}))\n")
    _write(tmp_path / "user_b.py",
           "from helper_b import run_case\n"
           "run_case(solver_kw={'enable_yield_gate': False})\n")
    hits = find_bad_kwargs(tmp_path)
    assert {(Path(f).name, k) for f, _, _, k in hits} == {("helper_a.py", "dispersion_axial")}


# --------------------------------------------------------------------------- #
# 豁免面：三条**同时**成立才算豁免（`raises` 上下文 ∧ 首个位置参数 TypeError ∧ tests/）
# --------------------------------------------------------------------------- #

def test_scanner_exempts_deliberate_negative_under_tests(tmp_path):
    """正例：`tests/` 下 `pytest.raises(TypeError)` 内的坏调用 = 故意反例，豁免。"""
    _write(tmp_path / "tests" / "negative.py",
           "import pytest\n"
           "from cemdisp.models2d import AnnulusD2DGASolver\n"
           "def test_x():\n"
           "    with pytest.raises(TypeError):\n"
           "        AnnulusD2DGASolver(dispersion_dt_scale=1.0)\n")
    assert find_bad_kwargs(tmp_path) == []


def test_scanner_flags_custom_raises_context_manager_outside_tests(tmp_path):
    """反例（R112）：非 `tests/` 文件里自定义 `def raises(e)` 上下文管理器**不豁免**。"""
    _write(tmp_path / "lib.py",
           "import contextlib\n"
           "from cemdisp.models2d import AnnulusD2DGASolver\n"
           "@contextlib.contextmanager\n"
           "def raises(e):\n"
           "    yield e\n"
           "def f():\n"
           "    with raises(TypeError):\n"
           "        AnnulusD2DGASolver(enable_d2dga_auto_m=True)\n")
    assert {k for *_, k in find_bad_kwargs(tmp_path)} == {"enable_d2dga_auto_m"}


def test_scanner_does_not_exempt_tests_file_without_typeerror(tmp_path):
    """第二条要件：`tests/` 下 `pytest.raises(ValueError)` 里传坏值**不算**反例。"""
    _write(tmp_path / "tests" / "negative.py",
           "import pytest\n"
           "from cemdisp.models2d import AnnulusD2DGASolver\n"
           "def test_x():\n"
           "    with pytest.raises(ValueError):\n"
           "        AnnulusD2DGASolver(dispersion_dt_scale=1.0)\n")
    assert {k for *_, k in find_bad_kwargs(tmp_path)} == {"dispersion_dt_scale"}


def test_scanner_does_not_exempt_calls_under_a_repo_without_tests_dir(tmp_path):
    """第三条要件：同一段代码放在没有 `tests/` 目录的树里 ⇒ 不豁免。"""
    _write(tmp_path / "negative.py",
           "import pytest\n"
           "from cemdisp.models2d import AnnulusD2DGASolver\n"
           "def test_x():\n"
           "    with pytest.raises(TypeError):\n"
           "        AnnulusD2DGASolver(dispersion_dt_scale=1.0)\n")
    assert {k for *_, k in find_bad_kwargs(tmp_path)} == {"dispersion_dt_scale"}


def test_scanner_still_flags_same_file_outside_raises(tmp_path):
    """豁免面窄：同一文件里 `raises` 块**之外**的坏调用照样被抓。"""
    _write(tmp_path / "tests" / "mixed.py",
           "import pytest\n"
           "from cemdisp.models2d import AnnulusD2DGASolver\n"
           "def test_x():\n"
           "    with pytest.raises(TypeError):\n"
           "        AnnulusD2DGASolver(dispersion_dt_scale=1.0)\n"
           "def test_y():\n"
           "    AnnulusD2DGASolver(enable_d2dga_auto_m=True)\n")
    assert {k for *_, k in find_bad_kwargs(tmp_path)} == {"enable_d2dga_auto_m"}


def test_scanner_supports_async_with_exemption(tmp_path):
    """`async with` 与 `with` 一视同仁（R117②）：异步反例同样豁免。"""
    _write(tmp_path / "tests" / "negative_async.py",
           "import pytest\n"
           "from cemdisp.models2d import AnnulusD2DGASolver\n"
           "async def test_x():\n"
           "    async with pytest.raises(TypeError):\n"
           "        AnnulusD2DGASolver(dispersion_dt_scale=1.0)\n")
    assert find_bad_kwargs(tmp_path) == []


def test_exemption_set_is_frozen_to_the_one_deliberate_negative():
    """全仓豁免实例集合冻结：只有 `test_no_invented_dispersion.py:57` 这一处。

    同时证明 Task 5 修好的那几处**不在**豁免集合里——它们是"改好了"，不是"被豁免"；
    第 4 处（`dispersion_scale_sensitivity_scan.py`）已归档 ⇒ 退出扫描范围。
    """
    result = gate._scan(ROOT)
    assert {(Path(h.path).name, h.lineno) for h in result.exempt} == {
        ("test_no_invented_dispersion.py", 57)}
    repaired = {"ht1_004_ablation.py", "closure_contribution_scan.py",
                "density_contrast_sensitivity_scan.py"}
    for hit in result.bad + result.exempt:
        assert Path(hit.path).name not in repaired, f"{hit} 应已修好，不该出现在任何一类"
    archived = (ROOT / "archive" / "scripts_probes_obsolete_2026-09-26"
                / "dispersion_scale_sensitivity_scan.py")
    assert archived.is_file(), "第 4 处修复点应已归档（保留历史、退出扫描面）"


# --------------------------------------------------------------------------- #
# 扫描器健壮性（R117）：盲区必须显式化，不得静默
# --------------------------------------------------------------------------- #

def test_scanner_reports_unparsable_files_instead_of_silence(tmp_path):
    """语法错误文件不得静默跳过：必须在 `skipped` 里报出文件名 + 异常类型（R117①）。"""
    _write(tmp_path / "broken.py", "def (:\n")
    _write(tmp_path / "fine.py", "x = 1\n")
    result = gate._scan(tmp_path)
    assert [(Path(p).name, kind) for p, kind in result.skipped] == [("broken.py", "SyntaxError")]


def test_target_classes_accept_no_var_keyword():
    """目标类当前**无** `**kwargs` ⇒ 关键字校验全开（R117③）。

    若某天构造签名新增 `**kwargs`，本条会红——提示闸门已自动降级为"该类停用 + 告警"
    （见下一条用例），而不是把全部关键字误报成"形参不存在"。
    """
    assert gate._var_keyword_targets() == set()


def test_scanner_disables_check_for_target_taking_var_keyword(tmp_path, monkeypatch):
    """目标类若新增 `**kwargs`：显式停用该类校验并记录告警，不产生全量误报（R117③）。"""
    _write(tmp_path / "sample.py",
           "from cemdisp.models2d import AnnulusD2DGASolver\n"
           "s = AnnulusD2DGASolver(anything_goes=True)\n")
    monkeypatch.setattr(gate, "_var_keyword_targets", lambda: {"AnnulusD2DGASolver"})
    result = gate._scan(tmp_path)
    assert result.bad == []
    assert result.disabled_targets == ["AnnulusD2DGASolver"]
