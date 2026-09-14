# -*- coding: utf-8 -*-
"""流函数速度场路径接入求解器的契约测试（Task 9，2026-09-15）。

三组契约：

1. **R7 冻结锚**（controller 裁定第 1 条）：``enable_stream_function=False``
   的输出与改动前 HEAD 76a91c1 的冻结锚**逐位一致**（float repr 往返精确 +
   float64 位级 sha256 摘要）。锚值出处：``.tmp_research/task9_probe/
   anchor_head76a91c1.json``（HEAD 76a91c1，refactor/d2dga-source-fidelity，
   Python 3.13.5 / numpy 2.1.3，2026-09-15 采集——改动前、任何 T9 编辑之前）。
2. **新路径接线与守恒**：默认 enable_stream_function=True；(4.22) 椭圆解在
   动力学时间步被消费；每列弧长通量 = q_half（T9 度量换算 D2）。
3. **新路径浮力判别与性能**：b 符号判别（论文 Table 3 机制，Z&F22 A32-22
   页竖直井机制）；ny=40, nz=250 单次椭圆解 < 0.15 s（controller 裁定第 7 条，
   Task 8 实测 ~24 ms）。
"""
from __future__ import annotations

import hashlib
import time

import numpy as np

import cemdisp.models2d.annulus_d2dga as _ann_mod
from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.well_spec import DepthValuePoint, WellSpec
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import AnnulusInletState
from cemdisp.models2d.stream_function import solve_stream_function
from cemdisp.runners.zhang2022_benchmark import (
    ZHANG2022_CASE_BY_ID,
    case_mean_velocity_m_s,
    paper_eta_e,
    paper_t_br,
    run_case,
)

# ---------------------------------------------------------------------------
# R7 冻结锚（出处见模块 docstring；逐位 = float64 位级）
# ---------------------------------------------------------------------------
_ANCHOR_PROVENANCE = (
    "HEAD 76a91c1 (refactor/d2dga-source-fidelity), 2026-09-15 采集于任何 T9 编辑前；"
    "Python 3.13.5 / numpy 2.1.3；脚本 .tmp_research/task9_probe/collect_anchor.py"
)
# field_like_nz48（合成现场级：井斜 6°、standoff 0.75→0.55、Bingham 四相、屈服门开）
_A_FIELD = {
    "eta_E": 0.10383284754023875,
    "eta_N": 0.014260352660175562,
    "cement_occ": 0.10383284754023875,
    "channeling": 0.9999999999966667,
    "mixing": 0.12922583337920518,
    "b_number": 10280.49171010473,
    "front_wide": 300.0,
    "front_narrow": 0.0,
    "front_mid": 6.382978723404255,
    "mean_wall": 0.19791666666666666,
    "cement_digest": "040eeffb4babc4d0",
    "spacer_digest": "c92ae788e3118928",
    "wall_digest": "93ca9c3502653840",
}
# zhang_case10_nz140（runner 标准口径）
_A_CASE10 = {
    "eta_E_paper_1p2L": 0.9998828053018046,
    "t_br_s": 4.274302414393974,
    "mass_err": 0.33544637519882115,
    "peak_ratio": 1.3175,
    "b_summary": 100.00045026354375,
    "cement_digest": "b7b659128b129598",
}


def _digest(arr: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(arr, dtype=np.float64).tobytes()).hexdigest()[:16]


def _field_like_well_spec() -> WellSpec:
    top, bottom = 100.0, 100.0 + 300.0
    return WellSpec(
        well_name="T9_anchor_field_like",
        top_md_m=top,
        bottom_md_m=bottom,
        shoe_md_m=bottom,
        hole_diameter_profile=(DepthValuePoint(top, 260.0), DepthValuePoint(bottom, 250.0)),
        liner_od_profile=(DepthValuePoint(top, 168.3), DepthValuePoint(bottom, 168.3)),
        inclination_profile=(DepthValuePoint(top, 2.0), DepthValuePoint(bottom, 6.0)),
        standoff_profile=(DepthValuePoint(top, 0.75), DepthValuePoint(bottom, 0.55)),
    )


def _field_like_fluids() -> tuple:
    mud = FluidSpec("a_mud", FluidRole.MUD, 1050.0, RheologyModel.BINGHAM,
                    plastic_viscosity_pa_s=0.022, yield_stress_pa=6.0)
    spacer = FluidSpec("a_spacer", FluidRole.SPACER, 1120.0, RheologyModel.BINGHAM,
                       plastic_viscosity_pa_s=0.030, yield_stress_pa=9.0)
    lead = FluidSpec("a_lead", FluidRole.LEAD, 1890.0, RheologyModel.BINGHAM,
                     plastic_viscosity_pa_s=0.085, yield_stress_pa=12.0)
    tail = FluidSpec("a_tail", FluidRole.TAIL, 1920.0, RheologyModel.BINGHAM,
                     plastic_viscosity_pa_s=0.120, yield_stress_pa=18.0)
    return mud, spacer, lead, tail


def _field_like_provider():
    q = 0.030 / 60.0  # m³/s

    def provider(t: float) -> AnnulusInletState:
        if t < 480.0:
            return AnnulusInletState(t, q, "spacer", (("spacer", 1.0),))
        if t < 1680.0:
            return AnnulusInletState(t, q, "lead", (("lead", 1.0),))
        return AnnulusInletState(t, q, "tail", (("tail", 1.0),))

    return provider


def _run_field_like(enable_stream_function: bool):
    solver = AnnulusD2DGASolver(
        nz=48, ny=24, dt=2.0, total_t=3000.0, enable_cfl_adaptive=False,
        open_outlet=True, e_clip_max=0.55,
        enable_stream_function=enable_stream_function,
    )
    return solver.run(_field_like_well_spec(), _field_like_fluids(), _field_like_provider())


class TestR7FrozenAnchorOldPathBitwise:
    """R7 护栏：enable_stream_function=False 逐位复现 HEAD 76a91c1 冻结锚。"""

    def test_field_like_bitwise(self):
        res = _run_field_like(enable_stream_function=False)
        s = res.summary["最终结果"]
        assert float(s["全井段最终有效顶替效率"]) == _A_FIELD["eta_E"]
        assert float(s["窄四分位效率"]) == _A_FIELD["eta_N"]
        assert float(s["最终水泥浆占据率"]) == _A_FIELD["cement_occ"]
        assert float(s["最终窜槽指数"]) == _A_FIELD["channeling"]
        assert float(s["最终混浆指数"]) == _A_FIELD["mixing"]
        assert float(s["浮力数_b"]) == _A_FIELD["b_number"]
        assert float(res.metrics["front_wide_m"].iloc[-1]) == _A_FIELD["front_wide"]
        assert float(res.metrics["front_narrow_m"].iloc[-1]) == _A_FIELD["front_narrow"]
        assert float(res.metrics["front_mid_m"].iloc[-1]) == _A_FIELD["front_mid"]
        assert float(res.metrics["mean_wall_mud"].iloc[-1]) == _A_FIELD["mean_wall"]
        assert _digest(res.cement_field) == _A_FIELD["cement_digest"]
        assert _digest(res.spacer_field) == _A_FIELD["spacer_digest"]
        assert _digest(res.wall_field) == _A_FIELD["wall_digest"]

    def test_zhang_case10_bitwise(self):
        case = ZHANG2022_CASE_BY_ID[10]
        row = run_case(case, nz=140, ny=40, keep_result=True,
                       solver_overrides={"enable_stream_function": False})
        res, ws = row["_result"], row["_well_spec"]
        w0 = case_mean_velocity_m_s(case)
        assert paper_eta_e(res, ws, w0) == _A_CASE10["eta_E_paper_1p2L"]
        t_br, _ = paper_t_br(res, ws, w0)
        assert t_br == _A_CASE10["t_br_s"]
        assert float(row["mass_conservation_error"]) == _A_CASE10["mass_err"]
        assert float(row["体积创造峰值"]) == _A_CASE10["peak_ratio"]
        assert float(res.summary["buoyancy_number"]) == _A_CASE10["b_summary"]
        assert _digest(res.cement_field) == _A_CASE10["cement_digest"]


class TestStreamFunctionPathWiring:
    """新路径接线：默认 True、椭圆解在动力学时间步被消费、列通量守恒。"""

    def test_default_flag_is_true(self):
        assert AnnulusD2DGASolver().enable_stream_function is True

    def test_dynamics_consumes_stream_function(self, monkeypatch):
        """动力学时间步（_compute_velocity 栈）必须消费 (4.22) 椭圆解。"""
        calls: list[str] = []
        real = _ann_mod.solve_stream_function

        def _spy(*args, **kwargs):
            caller = __import__("inspect").currentframe().f_back
            calls.append(caller.f_code.co_name)
            return real(*args, **kwargs)

        monkeypatch.setattr(_ann_mod, "solve_stream_function", _spy)
        _run_field_like(enable_stream_function=True)
        assert calls, "run() 全程未调用 solve_stream_function"
        assert all(name == "_velocity_stream_function" for name in calls)

    def test_new_path_column_flux_conservation(self):
        """新路径每列弧长通量 ∫w·b·dy = q_half（T9 推导 D2，梯形求积）。

        求积口径：模型体积核算（_trapez2d/half_volume/bulk_fill）全用
        ``np.trapezoid``——差分-梯形恒等式使其与模块 φ-度量不变量精确衔接
        （trap_y(w·b) = π·trap_φ(2r_aH·w_unit) = π ⇒ 缩放 q_half/π 后逐列
        精确 = q_half，任意浮力场）。
        """
        solver = AnnulusD2DGASolver(
            nz=24, ny=16, dt=4.0, total_t=400.0, enable_cfl_adaptive=False,
            enable_stream_function=True,
        )
        ws = WellSpec(
            well_name="flux_anchor", top_md_m=100.0, bottom_md_m=300.0, shoe_md_m=300.0,
            hole_diameter_profile=(DepthValuePoint(100.0, 260.0), DepthValuePoint(300.0, 260.0)),
            liner_od_profile=(DepthValuePoint(100.0, 168.3), DepthValuePoint(300.0, 168.3)),
            inclination_profile=(DepthValuePoint(100.0, 3.0), DepthValuePoint(300.0, 3.0)),
            standoff_profile=(DepthValuePoint(100.0, 0.7), DepthValuePoint(300.0, 0.7)),
        )
        mud = FluidSpec("m", FluidRole.MUD, 1050.0, RheologyModel.BINGHAM,
                        plastic_viscosity_pa_s=0.022, yield_stress_pa=6.0)
        tail = FluidSpec("t", FluidRole.TAIL, 1900.0, RheologyModel.BINGHAM,
                         plastic_viscosity_pa_s=0.12, yield_stress_pa=18.0)
        geom = solver._build_geom(ws)
        lead = np.full((solver.ny, solver.nz), 0.5)
        w, v, *_ = solver._compute_velocity(
            lead, np.zeros_like(lead), np.zeros_like(lead), geom, 0.03,
            np.full((solver.ny, solver.nz), 0.4), mud, tail, None, None, wall=None,
        )
        col_flux = np.trapezoid(w * geom["b"], x=geom["y"], axis=0)
        np.testing.assert_allclose(col_flux, 0.015, rtol=2e-3)


class TestStreamFunctionPathPhysics:
    """新路径浮力判别（论文 Table 3 机制方向）。"""

    def test_buoyancy_sign_discriminates_breakthrough(self):
        """重水泥（case2, b=+100）须比轻水泥（case1, b=−50）更晚突破、η_E 更高。

        论文 Table 3：case1 t_br=0.44/η_E=0.66 vs case2 t_br=0.95/η_E=0.95
        （Z&F22 竖直井机制：b 的 φ-梯度经 (4.22) 驱动次流，A32-22 页）。
        T9 首跑实测（nz=140）：t_br 0.119 vs 0.978——机制方向与量级由本测试锁定。
        """
        r1 = run_case(ZHANG2022_CASE_BY_ID[1], nz=60, ny=40)
        r2 = run_case(ZHANG2022_CASE_BY_ID[2], nz=60, ny=40)
        assert r2["t_br_模型"] - r1["t_br_模型"] > 0.10, (r1["t_br_模型"], r2["t_br_模型"])
        assert r2["eta_E_模型"] > r1["eta_E_模型"], (r1["eta_E_模型"], r2["eta_E_模型"])


class TestStreamFunctionSolvePerformance:
    """性能冒烟（controller 裁定第 7 条）：ny=40, nz=250 单次求解 < 0.15 s。"""

    def test_single_solve_under_150ms(self):
        y = np.linspace(0.0, np.pi * 0.1071, 40)
        phi = y / y[-1]
        nz = 250
        H = np.broadcast_to(0.0458 * (1 + 0.3 * np.cos(np.pi * phi))[:, None], (40, nz)).copy()
        geom = {"y": y, "phi": phi, "H": H, "b": 2 * H, "s": np.linspace(0, 2400.0, nz),
                "hole_mm": np.full((1, nz), 260.0), "od_mm": np.full((1, nz), 168.3)}
        c = 0.5 * (1.0 + np.tanh((geom["s"][None, :] - geom["s"].mean()) / 200.0))
        c = np.clip(np.broadcast_to(c * (0.7 + 0.3 * phi)[:, None], (40, nz)), 0.0, 1.0)
        b_field = np.zeros((2, 40, nz))
        b_field[0] = 8.0e4 * np.sin(np.pi * phi)[:, None] * c
        b_field[1] = 5.0e3 * np.tanh(np.linspace(-1, 1, nz))[None, :] * c
        t0 = time.perf_counter()
        psi = solve_stream_function(geom, c, np.ones_like(H), np.ones_like(H), 2.0, b_field)
        elapsed = time.perf_counter() - t0
        assert psi.shape == (40, nz)
        assert elapsed < 0.15, f"单次流函数求解 {elapsed*1000:.1f} ms 超 150 ms 预算"
