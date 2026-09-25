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

网格与时长（nz=30/ny=12/total_t=14000s）与**逐位锚 `_default_path_anchor_hu101.json`
完全同口径**（同 loader、同 T1 三开关、同默认环空开关），只多一个显式的
`enable_banded_solve`：这样两口径的差 **就是** 位移台账那一行所记录的量（η_E/η_N），
判据与留痕指向同一个数，不存在"另跑一个小算例自证"的口径漂移。
计算量 = 锚算例 × 2（各约 10 s），可接受。
"""
import numpy as np

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
