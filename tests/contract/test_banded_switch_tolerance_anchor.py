# -*- coding: utf-8 -*-
"""`enable_banded_solve` 两口径的容差锚（Task 8 Step 5 / 协调者裁定）。

存在理由
--------
`tests/contract/_default_path_anchor_hu101.json` 是**逐位**锚，只能锚定**一个**
口径（现为 banded=True 的新默认）。若没有本文件，"关掉带状求解器回到 spsolve"
这条回退路径就没有任何端到端守护——它照样能让全部测试变绿，回归无人发现。

口径
----
- `banded=True`（新默认）与 `banded=False`（历史 spsolve 口径）在**同一算例**上
  跑完整环空求解，比较 η_E / η_N；
- 容差取 **1e-2 pp**（即 η 的 1e-4）：依据 = 本任务实测同一微网格锚算例上
  |Δη_E| = 6.0e-4 pp、|Δη_N| = 1.5e-3 pp（见 `results/内部自洽加固_2026-09-25/
  位移台账.csv` 的"线性求解换带状Cholesky"行），取 1e-2 pp 相对实测最大位移
  （1.525e-3 pp）留约 **6.6×** 余量，既吸收"舍入被混沌放大"的规模效应，
  又远严于任何物理量级（1e-2 pp ≪ 判据阈值）；
- ⚠️ 本测试**不是**"零数值影响"的证明（R10）：换求解器**必然**让端到端逐位数字变化，
  证明路径是 `test_banded_solve_equivalence.py` 的离散系统等价性 + 本容差锚 +
  逐位锚重锚（R19）。
- 文件末尾另有 **solver 级不变量守卫**（修复轮 2，OPEN FINDING 2）：
  `test_solver_level_switch_reaches_every_solve_call_site` 直接在 solver 级断言
  `enable_banded_solve=False ⇒ 带状分解一次都不调用`（**直接覆盖 3 处流函数线性调用点**：
  线性路径、HB 冷启动启发式、HB 非线性入口），是**接线级的主守卫**；
  `test_banded_solve_equivalence.py` 的两条 spy 测试是更快的单元级守卫
  （分别钉住开关分支选择与非线性入口的形参透传）。
  ⚠️ **第 4 处调用点未被任何测试直接覆盖**：`_solve_stream_function_hb` 的 HB **回退分支**
  （`RuntimeError` ⇒ 牛顿线性闭包，`cemdisp/models2d/annulus_d2dga.py` 内
  `return solve_stream_function(..., banded=self.enable_banded_solve)`，即
  `annulus_d2dga.py:1857-1859`）在本窗口内未触发；它是**独立**调用点、带**自己的**
  `banded=` 透传，删掉它**不会被任何一个测试抓到**（本守卫对它是**间接**的：仅当回退
  分支被触发时才可能连带变红）。用户裁定 I-1：**不得**声称"覆盖全部 4 处"。

网格与时长（nz=30/ny=12/total_t=14000s）与**逐位锚 `_default_path_anchor_hu101.json`
完全同口径**（同 loader、同 T1 三开关、同默认环空开关），只多一个显式的
`enable_banded_solve`：这样两口径的差 **就是** 位移台账那一行所记录的量（η_E/η_N），
判据与留痕指向同一个数，不存在"另跑一个小算例自证"的口径漂移。
计算量 = 锚算例 × 2（各约 10 s），可接受。
"""
import warnings

import numpy as np
import pytest

from cemdisp.data.loaders import hu101_loader
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver

# 实测依据（本锚算例，见 results/内部自洽加固_2026-09-25/位移台账.csv）：
# |Δη_E| = 6.0e-4 pp、|Δη_N| = 1.5e-3 pp ⇒ 本容差相对实测最大位移（1.525e-3 pp）
# 留约 6.6× 余量（= 1e-2 / 1.525e-3）。
# 单位：pp（百分点）。
TOL_PP = 1.0e-2

# 与逐位锚完全同口径（微网格，非生产数字）
_NZ, _NY, _TOTAL_T = 30, 12, 14000.0


def _run(*, banded: bool) -> dict:
    well, fluids, schedule, _ = hu101_loader.load_hu101_tailpipe()
    cr = CasingFlowSolver(enable_gravity=True, mixing_contact_time=True,
                          plug_face_zero_mixing=True, has_plug=True).run(
        well, fluids, schedule)
    inlet = build_coupled_annulus_inlet_provider(
        cr, CasingFlowSolver(enable_gravity=True), fluids, split_cement_phases=True)
    solver = AnnulusD2DGASolver(total_t=_TOTAL_T, nz=_NZ, ny=_NY,
                                enable_banded_solve=banded)
    res = solver.run(well, fluids, inlet, schedule=schedule)
    fr = res.summary["最终结果"]
    return {
        "eta_E": float(fr["全井段最终有效顶替效率"]),
        "eta_N": float(fr["窄四分位效率"]),
    }


def test_banded_switch_round_trip_within_tolerance():
    """banded=True 与 False 必须回到同一物理答案（容差 1e-2 pp，见模块 docstring）。"""
    new = _run(banded=True)
    old = _run(banded=False)
    d_e = abs(new["eta_E"] - old["eta_E"]) * 100.0
    d_n = abs(new["eta_N"] - old["eta_N"]) * 100.0
    print(f"\n[两口径] eta_E: banded={new['eta_E']!r} spsolve={old['eta_E']!r} "
          f"Δ={d_e:.2e} pp")
    print(f"[两口径] eta_N: banded={new['eta_N']!r} spsolve={old['eta_N']!r} "
          f"Δ={d_n:.2e} pp")
    assert d_e <= TOL_PP, (
        f"η_E 两求解口径差 {d_e:.4g} pp 超容差 {TOL_PP} pp"
        "——两条路径未回到同一物理答案（不是舍入级差异，须停下上报）")
    assert d_n <= TOL_PP, (
        f"η_N 两求解口径差 {d_n:.4g} pp 超容差 {TOL_PP} pp"
        "——两条路径未回到同一物理答案（不是舍入级差异，须停下上报）")
    assert np.isfinite(new["eta_E"]) and np.isfinite(new["eta_N"])
    # 判别力下限：若两口径跑出**同一个** η_E，说明本测试没有真正驱动到 banded
    # 分支（例如开关没透传、或算例本身五场恒零），此时应当怀疑测试失效而非庆祝零差。
    assert d_e > 0.0 or d_n > 0.0, (
        "两口径 η 完全逐位相同——banded 开关可能未真正透传（本测试失去判别力），"
        "请先确认 `enable_banded_solve` 到达了 solve_stream_function")


# --------------------------------------------------------------------------- #
# solver 级**不变量**守卫（修复轮 2，OPEN FINDING 2）
#
# 为什么不做"单跳接线守卫"：``enable_banded_solve`` 在环空求解里共 **4 处**流函数线性
# 调用点（线性路径；HB 非线性路径的冷启动启发式 / 非线性入口 / 回退），单跳测试只能
# 覆盖其中 1 处、并会暗示并不存在的对称性覆盖（其余 3 处同样无跳级测试）。故改为在
# **solver 级**直接断言不变量：``enable_banded_solve=False`` ⇒ 整个 run 内带状分解
# **一次都不调用**（且 spsolve 确实被调用，防"空跑也能过"）；``=True`` ⇒ 至少调用一次。
# 两个参数（HB 关 / HB 开）合起来**直接覆盖 4 处中的 3 处**；第 4 处
# （``_solve_stream_function_hb`` 的 HB 回退分支，``annulus_d2dga.py:1857-1859``）
# 本窗口不触发，本守卫对它只是**间接**的（见本文件 docstring 的 ⚠️ 段与 :162-163）。
# --------------------------------------------------------------------------- #

# 小窗（开关**每一步**都被消费 ⇒ 无需跑到前缘）；实测线性 ~0.04 s / HB ~0.3 s。
_NZ_SMALL, _NY_SMALL, _TOTAL_T_SMALL = 12, 8, 60.0


@pytest.fixture(scope="module")
def _pipeline():
    """hu101 生产 loader + 管内段 1D + 入口桥（模块级建一次）。

    入口桥跨 run 复用与 ``.tmp_research/t8_time.py`` 的 A/B/A/B 口径一致（同一 provider
    连跑 4 次），故可安全缓存供两个参数化用例共享。
    """
    well, fluids, schedule, _ = hu101_loader.load_hu101_tailpipe()
    cr = CasingFlowSolver(enable_gravity=True, mixing_contact_time=True,
                          plug_face_zero_mixing=True, has_plug=True).run(
        well, fluids, schedule)
    inlet = build_coupled_annulus_inlet_provider(
        cr, CasingFlowSolver(enable_gravity=True), fluids, split_cement_phases=True)
    return well, fluids, schedule, inlet


def _spy_factorizations(monkeypatch):
    """给两条求解路径各挂计数器 spy（返回 ``{"banded": n, "sparse": n}``）。

    - ``cholesky_banded``：``solve_stream_function`` 在函数体内
      ``from scipy.linalg import cholesky_banded`` ⇒ 打 ``scipy.linalg`` 上的名字即命中；
    - ``spsolve``：模块级引用 ``stream_function.spla.spsolve``（历史路径逐字消费）。
    """
    import scipy.linalg as sla

    import cemdisp.models2d.stream_function as sf

    calls = {"banded": 0, "sparse": 0}
    real_cb = sla.cholesky_banded
    real_ss = sf.spla.spsolve

    def _cb(*args, **kwargs):
        calls["banded"] += 1
        return real_cb(*args, **kwargs)

    def _ss(*args, **kwargs):
        calls["sparse"] += 1
        return real_ss(*args, **kwargs)

    monkeypatch.setattr(sla, "cholesky_banded", _cb)
    monkeypatch.setattr(sf.spla, "spsolve", _ss)
    return calls


@pytest.mark.parametrize("hb", [False, True], ids=["线性路径", "HB非线性路径"])
def test_solver_level_switch_reaches_every_solve_call_site(monkeypatch, _pipeline, hb):
    """``enable_banded_solve`` 必须到达**每一处**流函数线性解（solver 级不变量）。

    覆盖情况（如实标注）：
      * ``hb=False`` 参数 = **线性路径**调用点；
      * ``hb=True``  参数 = **HB 非线性路径**（冷启动启发式 + 非线性入口，含首轮/
        欠松弛重解/收敛终解）；
      * **未**直接覆盖：HB 的**回退分支**（``RuntimeError`` ⇒ 牛顿线性闭包）——本窗口
        内未触发。它是**独立**调用点（``annulus_d2dga.py:1857-1859``，
        自带 ``banded=self.enable_banded_solve``），**不**与冷启动启发式共用透传
        ⇒ 删掉它的 ``banded=`` 不会被任何测试抓到；本守卫对它是**间接**的。
    """
    well, fluids, schedule, inlet = _pipeline
    calls = _spy_factorizations(monkeypatch)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")          # HB/dead-switch 告警与本断言无关
        AnnulusD2DGASolver(total_t=_TOTAL_T_SMALL, nz=_NZ_SMALL, ny=_NY_SMALL,
                           enable_banded_solve=False, enable_hb_closure=hb).run(
            well, fluids, inlet, schedule=schedule)
    assert calls["sparse"] >= 1, (
        f"enable_banded_solve=False 且 hb={hb} 时**一次线性解都没发生**——"
        "本用例空跑（窗口过短或路径未启用），下面的 0 次断言不具判别力")
    assert calls["banded"] == 0, (
        f"enable_banded_solve=False 仍调用带状分解 {calls['banded']} 次（hb={hb}）——"
        "开关在某处 solve_stream_function 调用点静默失效（漏透传）")

    calls["banded"] = calls["sparse"] = 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        AnnulusD2DGASolver(total_t=_TOTAL_T_SMALL, nz=_NZ_SMALL, ny=_NY_SMALL,
                           enable_banded_solve=True, enable_hb_closure=hb).run(
            well, fluids, inlet, schedule=schedule)
    assert calls["banded"] >= 1, (
        f"enable_banded_solve=True 一次带状分解都没调用（hb={hb}）——开关方向反了"
        "或透传丢失")
    assert calls["sparse"] == 0, (
        f"enable_banded_solve=True 仍调用 spsolve {calls['sparse']} 次（hb={hb}）——"
        "分支选择错误")
