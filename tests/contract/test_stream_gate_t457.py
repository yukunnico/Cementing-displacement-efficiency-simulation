# -*- coding: utf-8 -*-
"""Phase 2（2d）专项测试：T4 解析解+地板残余 / T5 连续门自解析极限 / T7 网格收敛趋势。

判据来源 = 续作计划 §3.3 改写版 + Phase 2 spec §5（**原判据不可用，勿照抄总纲**）：

- **T4 拆两条**：①Bingham 槽流解析解对照（复用 ``tests/contract/test_gap_solver.py:55``
  同构：``solve_fixed_G`` + 单流体 + ``u_exact`` 逐位）；②严格 ``wall≡1`` 单元
  （**排除过渡界面元**——界面元实测保留 26–36% 基线速度，不受地板约束）的地板残余
  ``I₁_eff/I₁ ≤ 1e-5``：``stream_function.py:435``
  ``conductance = np.maximum(1.0 - w_arr, _WALL_CONDUCTANCE_FLOOR)``，地板 = 1e-6
  ⇒ 严格冻结元比值恒 = 1e-6。**不得**用 "<0.1%/步" 作单测判据（仓内塞流解析解在
  gap_solver 路径，与 wall 地板不同路径）。
- **T5 连续门自解析极限**（单列构造，仿 ``tests/contract/test_yield_gate_continuous.py:46``）：
  ``wall_i = clip(1 − τw_extrap_i/(f·τy), 0, 1)`` 逐位自洽 + 两侧极限 + 参考元恒 0。
  ``_yield_gate_wall`` 自 2026-09-15 起为**连续**口径（无二值分支）；Balmforth 坡面式
  有自由面、本模型无 ⇒ **不得作解析真值**（§3.3）。
- **T7 网格收敛（判据重定义）**：pin 网格阶梯 {(20,80),(30,120),(40,160)}，
  gate-on/off **分别**做；"冻结过渡带" = ``res.wall_field`` 中 ``0<wall<1−δ`` 单元计数
  （如实记录）；判据 = **两级 |Δη_N| 单调递减**（趋势断言，**不写固定 <1%**——
  实测同一实现下随网格对可摆动 +24% ↔ +0.22%，固定阈值无区分力）。
  **Ruling R-2e-6（2026-10-07）**：合成算例取 ``total_t=600s``——240s 瞬态期两级 Δ
  非单调（噪声级 ≤0.35pp，gate 双档同象），600s 下 gate-on/off 双档单调成立
  （实测 0.66pp→0.43pp）；若错 = 趋势断言在瞬态期失去判别力，须按 total_t 复查。

测试红线（spec §5 / 续作 §3.3）：本文件**不写**求解器构造的 ``pytest.raises(TypeError)``
坏键反例（冻结豁免集不扩容）。
"""
from __future__ import annotations

import numpy as np
import pytest

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.well_spec import DepthValuePoint, WellSpec
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import AnnulusInletState


# ---------------------------------------------------------------------------
# T4-① Bingham 塞流解析解对照（gap_solver 路径；同构 test_gap_solver.py:55）
# ---------------------------------------------------------------------------
class TestT4AnalyticPlug:
    @pytest.mark.parametrize("G,tauY,kappa", [(10.0, 3.0, 1.0), (8.0, 2.0, 1.5)])
    def test_bingham_plug_and_yield_band_match_slot_solution(self, G, tauY, kappa):
        """塞流区 du/dy=0；屈服带与 Bingham 槽流解析解逐位一致（rtol 1e-6）。

        u(ỹ) = (G/2κ)(1−ỹ²) − (τY/κ)(1−ỹ)，塞流起于 y=τY/|G|（B&F25 (A1) 口径）。
        """
        from cemdisp.models2d.gap_solver import solve_fixed_G

        sol = solve_fixed_G(c_bar=0.0, n=(1.0, 1.0), kappa=(kappa, kappa),
                            tau_y=(tauY, tauY), G=(G, 0.0))
        assert sol.converged
        y_plug = tauY / G
        du = np.gradient(sol.u, sol.y)
        mask_plug = sol.y < y_plug - 1e-6
        assert np.max(np.abs(du[mask_plug])) < 1e-6, "塞流区存在残余剪切"
        mask_yield = sol.y > y_plug + 1e-6
        u_exact = (G / (2.0 * kappa)) * (1.0 - sol.y**2) - (tauY / kappa) * (1.0 - sol.y)
        assert np.allclose(sol.u[mask_yield], u_exact[mask_yield], rtol=1e-6, atol=1e-8)


# ---------------------------------------------------------------------------
# T4-② 严格 wall≡1 单元的地板残余（stream_function 消费路径）
# ---------------------------------------------------------------------------
class TestT4WallFloorResidual:
    def test_strict_frozen_cells_conductance_at_floor(self, monkeypatch):
        """严格 wall≡1 单元 I₁_eff/I₁ = 地板 1e-6 ≤ 1e-5；界面元/活跃元分列断言。

        捕获方式：np.maximum 纯透传 spy（第二实参 == _WALL_CONDUCTANCE_FLOOR 时记录
        输出 = conductance 场）。spy 只读不改——数值路径零扰动。
        """
        from cemdisp.models2d import stream_function as sf

        ny, nz = 8, 6
        phi = np.linspace(0.0, np.pi, ny)
        s = np.linspace(0.0, 300.0, nz)
        geom = {"phi": phi, "H": np.full((ny, nz), 0.045), "s": s,
                "hole_mm": 250.0, "od_mm": 168.3}
        c_bar = np.full((ny, nz), 0.5)
        b_field = np.zeros((2, ny, nz))
        # 均匀幅值输入源项为零（R29 口径）⇒ b_φ 须携带 φ-梯度才有非平凡源
        b_field[0] = 0.1 * np.cos(phi)[:, None]
        b_field[1] = 0.05 * np.sin(phi)[:, None]

        wall = np.zeros((ny, nz))
        wall[3:5, :] = 1.0     # 严格冻结带（判据对象）
        wall[2, :] = 0.5       # 过渡界面元（排除——实测保留 26–36% 基线速度）

        captured: list[np.ndarray] = []
        real_max = np.maximum

        def spy(a, b, *args, **kw):
            out = real_max(a, b, *args, **kw)
            try:
                if np.isscalar(b) and abs(float(b) - sf._WALL_CONDUCTANCE_FLOOR) < 1e-15:
                    captured.append(np.array(out, dtype=float, copy=True))
            except Exception:
                pass
            return out

        monkeypatch.setattr(np, "maximum", spy)
        psi = sf.solve_stream_function(geom, c_bar, 0.05, 0.3, 0.05 / 0.3,
                                       b_field, wall=wall)

        assert captured, ("未捕获地板调用——wall 消费位点已迁移，本测试须随迁"
                          "（不得静默失效）")
        ratio = captured[0]                      # I₁_eff/I₁ = conductance
        strict = wall == 1.0
        trans = (wall > 0.0) & (wall < 1.0)
        active = wall == 0.0
        assert np.all(ratio[strict] <= 1e-5), \
            f"严格冻结元地板残余 > 1e-5: max={ratio[strict].max():.3e}"
        assert np.allclose(ratio[strict], sf._WALL_CONDUCTANCE_FLOOR), \
            "严格冻结元比值应恰为地板值"
        assert np.all(ratio[trans] > 1e-5), "过渡界面元被误入地板判据（应排除）"
        assert np.allclose(ratio[active], 1.0), "活跃元必须零扰动（比值=1）"
        # 单位通量 BC 与有限性（椭圆结构不变的旁证）
        assert np.all(np.isfinite(psi))
        assert np.allclose(psi[0], 0.0) and np.allclose(psi[-1], 1.0)

    def test_wall_none_path_bitwise_untouched(self):
        """wall=None ⇒ I₁_eff = I₁（逐位无扰，关2 语义在算子层的旁证）。"""
        from cemdisp.models2d import stream_function as sf

        ny, nz = 6, 4
        phi = np.linspace(0.0, np.pi, ny)
        geom = {"phi": phi, "H": np.full((ny, nz), 0.04), "s": np.linspace(0, 200, nz),
                "hole_mm": 250.0, "od_mm": 168.3}
        c_bar = np.full((ny, nz), 0.4)
        b_field = np.zeros((2, ny, nz))
        b_field[0] = 0.1 * np.cos(phi)[:, None]
        psi_none = sf.solve_stream_function(geom, c_bar, 0.05, 0.3, 0.05 / 0.3,
                                            b_field, wall=None)
        psi_zero = sf.solve_stream_function(geom, c_bar, 0.05, 0.3, 0.05 / 0.3,
                                            b_field, wall=np.zeros((ny, nz)))
        # wall=0 场 ⇒ conductance=max(1,1e-6)=1 ⇒ 与 None 路径同一算术结果
        assert np.array_equal(psi_none, psi_zero)


# ---------------------------------------------------------------------------
# T5 连续屈服门自解析极限（单列构造，仿 test_yield_gate_continuous.py:46）
# ---------------------------------------------------------------------------
_F_SAFETY = 1.15
_TAU_Y = 8.0


def _single_column(tau_w_targets_tail, b_ref=0.04, mu=0.1, w_ref=1.0):
    """row0 = 参考元（唯一 w>0）；第 i≥1 行外推壁剪精确钉为 targets[i-1]。

    外推链 τw_extrap_i = τw_ref·(b_i/b_ref)，τw_ref = mu·6·w_ref/b_ref
    ⇒ 反解 b_i = b_ref·τw_i/τw_ref（与既有连续门测试同构）。
    """
    tau_w_ref = mu * 6.0 * w_ref / b_ref
    tail = np.asarray(tau_w_targets_tail, dtype=float)
    b = np.concatenate(([b_ref], b_ref * tail / tau_w_ref))[:, None]
    ny = b.shape[0]
    w = np.zeros((ny, 1))
    w[0, 0] = w_ref
    return w, b, np.full((ny, 1), mu), tau_w_ref


class TestT5ContinuousGateAnalyticLimit:
    def test_wall_matches_own_analytic_formula_bitwise(self):
        """wall_i = clip(1 − τw_extrap_i/(f·τy), 0, 1) **逐位**自洽（atol=0）。

        自解析极限：连续门的解析真值 = 其自身闭式（Balmforth 坡面式有自由面，
        本模型无，不得作真值——§3.3）。覆盖阈值内/恰等/阈值外/深外四段。
        """
        thr = _F_SAFETY * _TAU_Y                      # 9.2
        targets = np.array([0.25 * thr, 0.5 * thr, 0.9 * thr,
                            thr, 1.5 * thr, 100.0 * thr])
        w, b, mu_reg, _ = _single_column(targets)
        ny = b.shape[0]
        cement_ever = np.ones((ny, 1))
        cement_local = np.full((ny, 1), 0.9)
        tau_y = np.full((ny, 1), _TAU_Y)
        wall = AnnulusD2DGASolver._yield_gate_wall(
            w, b, mu_reg, tau_y, cement_ever, cement_local, _F_SAFETY)
        expected = np.clip(1.0 - targets / thr, 0.0, 1.0)
        assert np.allclose(wall[1:, 0], expected, rtol=0.0, atol=1e-12), \
            f"连续门偏离自解析式: {wall[1:,0]} vs {expected}"
        # 两侧极限与边界语义：τw<f·τy ⇒ wall>0（阈值内冻结度）；τw≥f·τy ⇒ wall=0
        assert wall[1, 0] == pytest.approx(0.75, abs=1e-12)   # 0.25·thr
        assert wall[4, 0] == 0.0                              # 恰等归可流动侧
        assert wall[-1, 0] == 0.0                             # 深外
        assert wall[0, 0] == 0.0                              # 参考元恒 0（R3 不变量①）

    def test_zero_tauy_gives_no_freeze_and_no_nan(self):
        """τy=0 ⇒ f·τy=0 ⇒ 无屈服应力 ⇒ wall≡0 且无 NaN/inf（0/0 防护语义）。"""
        w, b, mu_reg, _ = _single_column([1.0, 5.0, 20.0])
        ny = b.shape[0]
        wall = AnnulusD2DGASolver._yield_gate_wall(
            w, b, mu_reg, np.zeros((ny, 1)), np.ones((ny, 1)),
            np.full((ny, 1), 0.9), _F_SAFETY)
        assert np.all(np.isfinite(wall))
        assert np.all(wall == 0.0)


# ---------------------------------------------------------------------------
# T7 网格收敛趋势（pin 阶梯，gate-on/off 分别做）
# ---------------------------------------------------------------------------
_Q_M3S = 1.0 / 60.0          # 生产级排量（低流速会把 Bingham 泥浆格点推上屈服悬崖）
_LADDER = ((20, 80), (30, 120), (40, 160))
_TOTAL_T = 600.0             # R-2e-6：240s 瞬态期趋势为噪声级非单调，600s 双档单调


def _t7_well() -> WellSpec:
    top, bottom = 100.0, 400.0
    return WellSpec(
        well_name="t7_mesh_trend",
        top_md_m=top, bottom_md_m=bottom, shoe_md_m=bottom,
        hole_diameter_profile=(DepthValuePoint(top, 260.0), DepthValuePoint(bottom, 250.0)),
        liner_od_profile=(DepthValuePoint(top, 168.3), DepthValuePoint(bottom, 168.3)),
        inclination_profile=(DepthValuePoint(top, 2.0), DepthValuePoint(bottom, 5.0)),
        standoff_profile=(DepthValuePoint(top, 0.75), DepthValuePoint(bottom, 0.60)),
    )


def _t7_fluids() -> tuple:
    mud = FluidSpec("mud", FluidRole.MUD, 1050.0, RheologyModel.BINGHAM,
                    plastic_viscosity_pa_s=0.022, yield_stress_pa=6.0)
    spacer = FluidSpec("spacer", FluidRole.SPACER, 1120.0, RheologyModel.BINGHAM,
                       plastic_viscosity_pa_s=0.030, yield_stress_pa=9.0)
    lead = FluidSpec("lead", FluidRole.LEAD, 1890.0, RheologyModel.POWER_LAW,
                     power_law_n=0.75, consistency_k=0.55)
    tail = FluidSpec("tail", FluidRole.TAIL, 1920.0, RheologyModel.POWER_LAW,
                     power_law_n=0.70, consistency_k=0.90)
    return mud, spacer, lead, tail


def _t7_provider(t: float) -> AnnulusInletState:
    if t < 40.0:
        return AnnulusInletState(t, _Q_M3S, "spacer", (("spacer", 1.0),))
    if t < 420.0:
        return AnnulusInletState(t, _Q_M3S, "lead", (("lead", 1.0),))
    return AnnulusInletState(t, _Q_M3S, "tail", (("tail", 1.0),))


@pytest.mark.parametrize("gate", [False, True], ids=["gate_off", "gate_on"])
def test_t7_two_level_eta_n_difference_monotonically_decreases(gate):
    """两级 |Δη_N| 单调递减（趋势断言）；冻结过渡带计数如实记录在失败信息。

    不写固定 <1% 阈值（§3.3：无区分力）；Δ 全零 = 断言失去判别力 ⇒ 响亮失败。
    """
    etas, bands = [], []
    for ny, nz in _LADDER:
        solver = AnnulusD2DGASolver(
            nz=nz, ny=ny, dt=4.0, total_t=_TOTAL_T,
            enable_cfl_adaptive=True, enable_stream_yield_gate=gate,
            open_outlet=True,
        )
        res = solver.run(_t7_well(), _t7_fluids(), _t7_provider)
        etas.append(float(res.summary["最终结果"]["窄四分位效率"]))
        wf = np.asarray(res.wall_field, dtype=float)
        bands.append(int(np.sum((wf > 0.0) & (wf < 1.0 - 1e-6))))
    d1, d2 = abs(etas[1] - etas[0]), abs(etas[2] - etas[1])
    assert d1 > 0.0 and d2 > 0.0, (
        f"两级 Δ 出现零值（η_N={etas}）——阶梯失去判别力，检查算例构造")
    assert d2 < d1, (
        f"|Δη_N| 未单调递减: {d1:.6f} → {d2:.6f}；η_N={etas}；"
        f"冻结过渡带(0<wall<1−δ)计数={bands}")
