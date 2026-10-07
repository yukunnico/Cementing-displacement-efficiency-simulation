# -*- coding: utf-8 -*-
"""Phase 4b 契约测试：2D 环空逐列温度（`enable_depthwise_temperature`）（2026-10-07）。

权威 = docs/superpowers/specs/2026-10-07-phase4-depthwise-temperature-design.md
§0 F-1/F-2/F-3/F-4 裁定 + §1 项 2/3/4 + §2「4b 验收」；行号 = 同日期测绘调研
项 5/6/7/8/9/10/11；四站点同传红线 = 2026-10-06 Phase 2 spec §2.2。

五段关2 + 关3 + 关5 + 四站点同传：

1. **关2① 不传场/不开逐列 ⇒ 短算例逐位 = HEAD**（4a 同型锚值法：改前 HEAD
   捕获的 summary 浮点 + 场字节哈希 + `_compute_props` 直调哈希 + 派生参数）。
2. **关2② Constant 场 + 开逐列 ⇒ 快捷路径命中（计数>0）且输出与不开逐列逐位
   相等**（`_representative_temperature` :1226-1228 短路保留为标量）。
3. **关2③ two_layer 标量分支逐位**（math.sqrt→np.sqrt 替换后标量入标量出逐位；
   测试内嵌改前参照实现）。
4. **关2④ `_scalar_viscosity` 标量入逐位 + 广播形状放宽**（(1,nz)/(ny,nz) 列场
   透传；其他非均匀仍 raise；全场一致数组仍取首值）。
5. **关2⑤ T6 同型重钉**：λ_op·(1/F²) 乘积不变量在「列场退化标量」（Constant
   场+开逐列）下逐位 = 现状（F-3 代价条款；λ_op 与 F² 站点记录器逐位对照）。
6. **关3**：静温梯度场逐列 vs 域均标量对照——深段（高温）τ₀/μp 上翘显式化
   （组 B 水泥 T>100 一次式上行段）；列参数 == 手写 `fluid_at` 逐列（精确）；
   动力学响应（ON≠OFF）如实。
7. **关5**：memo 尺寸上界（静温档逐列派生每 run 每相一次，不随步数增长）；
   审计环聚合计数正确性 + 环不爆（<定长上限）；`np.maximum` 列位点抽查
   （`_compute_props` :1453 型位点在列参下不抛 python-max TypeError）。
8. **四站点同传**：逐列 + `include_yield_term=True` ⇒ η 列 == 逐列
   `fluid_apparent_viscosity(include=True)`（拆分旗标贯穿列站点），λ_op/F²/b_num
   仍域均标量（F-3 边界），组合下负对照不破坏（ON-split 与 OFF-split 动力学可分）。

探针 A（呼101 列化关 hard 锚）由 results/_probe_4b逐列_2026-10-07/ 承担，不在本文件。
"""
from __future__ import annotations

import hashlib
import math

import numpy as np
import pytest

import cemdisp.models2d.stream_function as sf_module
from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.rheology_vs_temperature import fluid_at
from cemdisp.data.temperature_field import ConstantTemperatureField
from cemdisp.data.well_spec import DepthValuePoint, WellSpec
from cemdisp.models2d import buoyancy
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import AnnulusInletState
from cemdisp.models2d.two_layer import mobility_i1, mobility_i2

_Q_M3S = 1.0 / 60.0

# --------------------------------------------------------------------------- #
# HEAD 逐位锚（2026-10-07 改前捕获 = results/_probe_4b逐列_2026-10-07/head_anchors.json；
# 与 capture_head_anchors_20261007.py 同装配；禁手改——重捕获脚本可复证）
# --------------------------------------------------------------------------- #
_ANCHOR_TOFF = {
    "eta_E": 0.3946295373348512, "eta_N": 0.36065145369450724,
    "b": 308.41786014961826,
    "cement_sha256": "cdff15e5fe098888fa105854609b97429a1e72a15f3f4c2b299691afca5a0a3b",
}
_ANCHOR_CONSTANT = {
    "eta_E": 0.40503766614396924, "eta_N": 0.3466329656571022,
    "b": 91.33759550293782,
    "cement_sha256": "a5a0ac156cb46f5f068fb0ab1cdd7a16938175cf8790b350fe65d5c5e66219c7",
    "props_sha256": [
        "02d7d61e3608a9bec983b226a07b8ce88444a5fc0ec947a00ef58ba1ba54b331",
        "9c1230c26a0f12bfa6c323d5d6e6baa95ec6aaa61f94fa14edf242618057c70e",
        "07c6736dd542b33c6736aca402427957781e2bd5b71da47eac73c6cec9fe69a9",
        "cb1a222f4cdc618fe4715d3f36f6748480f3dc8fe73d16a3843cf378ff18b04f",
        "078adf21b908a4e5935d7c7c6bf7865c27542dc87470f4292c67a3adc519392e",
        "5f48f5e565357e7364d8366da82af95f246bf04e094e6748c27b70be5781c8ab",
        "a577809c8d12ed6f97074b3c58b469fc8aaabce92cd4d71fd50d2eb1aa70d21e",
        "825fae38cced1bd686db1cf9e5f1ed4a53d370281f30467164fe15b0834b66f0",
        "a43548966ad56af3db2fad43f0b4daf0d1619f3084e01d3fe5e72f42c7850898",
    ],
    "phase_props": {
        "mud": (0.07428696678437698, 12.968444096844806),
        "lead": (0.09151000000000001, 2.85208),
        "tail": (0.09151000000000001, 2.85208),
        "spacer": (0.0580395677469, 6.702165213137),
    },
}
_ANCHOR_RAMP = {
    "eta_E": 0.40497080139046954, "eta_N": 0.33133718679127844,
    "b": 57.538683764440414,
    "cement_sha256": "b3b38f021162f837c8952e98ec4ecb4954dbb774ca1b4114b53ef5e07924795b",
    "rep_T": 12.740000000000002,
    "props_sha256": [
        "723cfb9333d31cfd07fef9bd48a0b9efbd9da53b29dde497d7f30162ac2d4269",
        "9c1230c26a0f12bfa6c323d5d6e6baa95ec6aaa61f94fa14edf242618057c70e",
        "07c6736dd542b33c6736aca402427957781e2bd5b71da47eac73c6cec9fe69a9",
        "007c2c16fc40ce9eeea2db6ceb4e83f0a9e778b2bb78b87af8bc91cc78ee1b4b",
        "52ad7108412b035f939959d699fdf506fa56752b9d365893b3cba2370ae55d96",
        "28ec58b4c2360969f813735789c648d6ad2a4729a35ede2251790a2f55b6ec4d",
        "ae488801a62a70a35115152c615fea21f6d803cc1e354daa08fc9d3d92f85252",
        "825fae38cced1bd686db1cf9e5f1ed4a53d370281f30467164fe15b0834b66f0",
        "2df55a0f540b92e4a9a5bd04150e6ea19f7500c4efe478cddbbe35fd5a32dbee",
    ],
}


# --------------------------------------------------------------------------- #
# 短算例装配（与 capture_head_anchors_20261007.py / test_temperature_phase_props 同源）
# --------------------------------------------------------------------------- #

def _well_spec() -> WellSpec:
    top, bottom = 100.0, 400.0
    return WellSpec(
        well_name="P4b_columnwise_2d",
        top_md_m=top, bottom_md_m=bottom, shoe_md_m=bottom,
        hole_diameter_profile=(DepthValuePoint(top, 260.0), DepthValuePoint(bottom, 250.0)),
        liner_od_profile=(DepthValuePoint(top, 168.3), DepthValuePoint(bottom, 168.3)),
        inclination_profile=(DepthValuePoint(top, 2.0), DepthValuePoint(bottom, 5.0)),
        standoff_profile=(DepthValuePoint(top, 0.75), DepthValuePoint(bottom, 0.60)),
    )


def _fluids() -> tuple:
    mud = FluidSpec("mud", FluidRole.MUD, 1050.0, RheologyModel.BINGHAM,
                    plastic_viscosity_pa_s=0.022, yield_stress_pa=6.0)
    spacer = FluidSpec("spacer", FluidRole.SPACER, 1120.0, RheologyModel.BINGHAM,
                       plastic_viscosity_pa_s=0.030, yield_stress_pa=9.0)
    lead = FluidSpec("lead", FluidRole.LEAD, 1890.0, RheologyModel.POWER_LAW,
                     power_law_n=0.75, consistency_k=0.55)
    tail = FluidSpec("tail", FluidRole.TAIL, 1920.0, RheologyModel.POWER_LAW,
                     power_law_n=0.70, consistency_k=0.90)
    return mud, spacer, lead, tail


def _provider(t: float) -> AnnulusInletState:
    if t < 40.0:
        return AnnulusInletState(t, _Q_M3S, "spacer", (("spacer", 1.0),))
    if t < 420.0:
        return AnnulusInletState(t, _Q_M3S, "lead", (("lead", 1.0),))
    return AnnulusInletState(t, _Q_M3S, "tail", (("tail", 1.0),))


def _solver(**overrides) -> AnnulusD2DGASolver:
    params = dict(nz=24, ny=12, dt=2.0, total_t=240.0, enable_cfl_adaptive=False,
                  open_outlet=True)
    params.update(overrides)
    return AnnulusD2DGASolver(**params)


class _RampField:
    """斜坡温度场（t 相关 ⇒ 表格场同型；列每步变）：10 + 0.01·md + 0.001·t。"""

    def T(self, md_m: float, t_s: float) -> float:
        return 10.0 + 0.01 * float(md_m) + 0.001 * float(t_s)


class _StaticGradientField:
    """静温梯度场（t 无关 ⇒ 4b 主验收档）：T = 120 + 0.02·md，全程 T>100°C。

    取 >100°C 段 = 组 B 水泥 τ₀(T) 一次式**上翘**区间（0.03206·T+1.979，Phase 0
    公式层），关3 方向断言的数据前提。
    """

    def T(self, md_m: float, t_s: float) -> float:
        return 120.0 + 0.02 * float(md_m)


def _sha(arr) -> str:
    return hashlib.sha256(np.asarray(arr, dtype=float).tobytes()).hexdigest()


def _run_summary(res) -> tuple:
    f = res.summary["最终结果"]
    return (float(f["全井段最终有效顶替效率"]), float(f["窄四分位效率"]),
            float(f["浮力数_b"]))


def _direct_props(solver: AnnulusD2DGASolver, res, t: float = 240.0):
    """`_compute_props` 末步 geom 直调（固定输入 ⇒ 场公式列参位点抽查）。"""
    geom = res.geom
    mud, spacer, lead, tail = _fluids()
    w_prev = np.full((solver.ny, solver.nz), 0.45)
    lead0 = np.zeros((solver.ny, solver.nz)); lead0[:, 0] = 0.5
    zero = np.zeros_like(lead0)
    return solver._compute_props(lead0, zero, zero, w_prev, geom,
                                 mud, lead, tail, spacer, t=t)


# --------------------------------------------------------------------------- #
# 关2①：不开逐列（且不传场 / T-on）⇒ 逐位 = HEAD 锚
# --------------------------------------------------------------------------- #
class TestOffBitwiseHeadAnchor:
    def test_toff_short_run_bitwise(self) -> None:
        solver = _solver()
        res = solver.run(_well_spec(), _fluids(), _provider)
        assert _run_summary(res) == (_ANCHOR_TOFF["eta_E"], _ANCHOR_TOFF["eta_N"],
                                     _ANCHOR_TOFF["b"])
        assert _sha(res.cement_field) == _ANCHOR_TOFF["cement_sha256"]
        assert solver._col_memo == {} and solver._col_batches == 0
        assert solver._col_audit_counts == {}

    def test_ton_constant_bitwise(self) -> None:
        solver = _solver(enable_temperature_rheology=True)
        res = solver.run(_well_spec(), _fluids(), _provider,
                         temperature_field=ConstantTemperatureField(60.0))
        assert _run_summary(res) == (_ANCHOR_CONSTANT["eta_E"],
                                     _ANCHOR_CONSTANT["eta_N"],
                                     _ANCHOR_CONSTANT["b"])
        assert _sha(res.cement_field) == _ANCHOR_CONSTANT["cement_sha256"]
        assert [_sha(a) for a in _direct_props(solver, res)] == \
            _ANCHOR_CONSTANT["props_sha256"]
        # 派生四相参数逐位（`_phase_props` 标量路径未被列版改动波及）
        solver._phase_memo.clear()
        geom = res.geom
        der = {n: solver._phase_props(f, geom, 240.0) for n, f in
               zip(("mud", "spacer", "lead", "tail"), _fluids())}
        for n, (pv, ty) in _ANCHOR_CONSTANT["phase_props"].items():
            assert der[n].plastic_viscosity_pa_s == pv
            assert der[n].yield_stress_pa == ty

    def test_ton_ramp_columnwise_off_bitwise(self) -> None:
        """T-on 非均匀代表标量（域均）链逐位 = HEAD：列化关时均值分支不动。"""
        solver = _solver(enable_temperature_rheology=True)  # depthwise 默认关
        res = solver.run(_well_spec(), _fluids(), _provider, temperature_field=_RampField())
        assert _run_summary(res) == (_ANCHOR_RAMP["eta_E"], _ANCHOR_RAMP["eta_N"],
                                     _ANCHOR_RAMP["b"])
        assert _sha(res.cement_field) == _ANCHOR_RAMP["cement_sha256"]
        assert solver._representative_temperature(res.geom, 240.0) == _ANCHOR_RAMP["rep_T"]
        assert [_sha(a) for a in _direct_props(solver, res)] == _ANCHOR_RAMP["props_sha256"]


# --------------------------------------------------------------------------- #
# 关2②：Constant 场 + 开逐列 ⇒ 快捷路径命中且与标量路径逐位
# --------------------------------------------------------------------------- #
class TestConstantDegenerateBitwise:
    def test_uniform_shortcut_hits_and_bitwise(self) -> None:
        off = _solver(enable_temperature_rheology=True)
        on = _solver(enable_temperature_rheology=True, enable_depthwise_temperature=True)
        r_off = off.run(_well_spec(), _fluids(), _provider,
                        temperature_field=ConstantTemperatureField(60.0))
        r_on = on.run(_well_spec(), _fluids(), _provider,
                      temperature_field=ConstantTemperatureField(60.0))
        assert on._uniform_field_hits > 0           # 快捷路径命中计数（关2 证据）
        assert on._col_batches == 0                  # 均匀 ⇒ 零列派生批
        assert on._col_memo == {}
        assert _run_summary(r_off) == _run_summary(r_on)
        assert r_off.cement_field.tobytes() == r_on.cement_field.tobytes()
        assert _sha(_direct_props(on, r_on, t=240.0)[0]) == \
            _sha(_direct_props(off, r_off, t=240.0)[0])
        # λ_op/γ̇_rep 观测位逐位（均纯标量）
        assert on._lambda_op_last == off._lambda_op_last
        assert on._shear_rate_rep_last == off._shear_rate_rep_last


# --------------------------------------------------------------------------- #
# 关2③：two_layer np.sqrt 替换后标量分支逐位
# --------------------------------------------------------------------------- #
def _ref_mobility_i1(c_bar, m, eta1=1.0, eta2=1.0, H=1.0):
    """改前 HEAD 参照（math.sqrt 版，逐字）。"""
    c = np.asarray(c_bar, dtype=float)
    sq_m = math.sqrt(m)
    out = (sq_m * c**3 + (1.0 - c**3) / sq_m) / 3.0
    out = out * (H**3 / math.sqrt(eta1 * eta2))
    return float(out) if np.isscalar(c_bar) else out.astype(float, copy=False)


def _ref_mobility_i2(c_bar, m, eta1=1.0, eta2=1.0, H=1.0):
    c = np.asarray(c_bar, dtype=float)
    sq_m = math.sqrt(m)
    out = (2.0 * sq_m * c**3 * (1.0 - c) + c * (1.0 - c) ** 2 * (1.0 + 2.0 * c) / sq_m) / 6.0
    out = out * (H**4 / math.sqrt(eta1 * eta2))
    return float(out) if np.isscalar(c_bar) else out.astype(float, copy=False)


class TestTwoLayerScalarBitwise:
    def test_scalar_grid_bitwise(self) -> None:
        rng = np.random.default_rng(20261007)
        cs = list(np.linspace(0.0, 1.0, 21)) + [0.003, 0.997, float(rng.random())]
        ms = [1e-3, 0.1, 0.356, 1.0, 2.7, 10.0, float(rng.uniform(0.01, 5))]
        etas = [0.02, 0.3, 1.0, 3.0, float(rng.uniform(0.05, 2))]
        for c in cs:
            for m in ms:
                for e1, e2 in zip(etas, etas[::-1]):
                    for f_new, f_ref in ((mobility_i1, _ref_mobility_i1),
                                         (mobility_i2, _ref_mobility_i2)):
                        a = f_new(c, m, eta1=e1, eta2=e2, H=0.05)
                        b = f_ref(c, m, eta1=e1, eta2=e2, H=0.05)
                        assert isinstance(a, float) and a == b
        # 数组入（c_bar 场）同样逐位
        c_arr = np.linspace(0.0, 1.0, 37)
        for f_new, f_ref in ((mobility_i1, _ref_mobility_i1), (mobility_i2, _ref_mobility_i2)):
            assert np.array_equal(f_new(c_arr, 0.5, eta1=0.2, eta2=0.6, H=0.04),
                                  f_ref(c_arr, 0.5, eta1=0.2, eta2=0.6, H=0.04))

    def test_column_m_broadcast_field(self) -> None:
        """m/η 列场（(1,nz)）⇒ 逐列独立 (4.21)：与逐列标量调用等值。"""
        nz = 6
        m_col = np.linspace(0.2, 3.0, nz).reshape(1, nz)
        e1_col = np.linspace(0.05, 0.4, nz).reshape(1, nz)
        e2_col = np.linspace(0.1, 0.9, nz).reshape(1, nz)
        c = np.full((3, nz), 0.4)
        H = np.full((3, nz), 0.05)
        got = mobility_i1(c, m_col, eta1=e1_col, eta2=e2_col, H=H)
        for j in range(nz):
            per_col = _ref_mobility_i1(c[:, j], float(m_col[0, j]),
                                       eta1=float(e1_col[0, j]),
                                       eta2=float(e2_col[0, j]), H=H[:, j])
            np.testing.assert_array_equal(got[:, j], np.asarray(per_col))
        got2 = mobility_i2(c, m_col, eta1=e1_col, eta2=e2_col, H=H)
        for j in range(nz):
            per_col = _ref_mobility_i2(c[:, j], float(m_col[0, j]),
                                       eta1=float(e1_col[0, j]),
                                       eta2=float(e2_col[0, j]), H=H[:, j])
            np.testing.assert_array_equal(got2[:, j], np.asarray(per_col))

    def test_scalar_in_field_c_array_out(self) -> None:
        """c_bar 标量 + m 数组 ⇒ 场出（不炸 float(out)；测绘项10 返回分支重排）。"""
        out = mobility_i1(0.5, np.array([[0.3, 0.7, 1.2]]))
        assert isinstance(out, np.ndarray) and out.shape == (1, 3)


# --------------------------------------------------------------------------- #
# 关2④：`_scalar_viscosity` 标量入逐位 + 广播形状放宽
# --------------------------------------------------------------------------- #
class TestScalarViscosityRelaxation:
    def test_scalar_and_uniform_bitwise(self) -> None:
        assert sf_module._scalar_viscosity(0.42, "e") == 0.42
        assert sf_module._scalar_viscosity(np.float64(0.42), "e") == 0.42
        # 全场一致数组 ⇒ 仍取首值（旧捷径不动）
        arr = np.full((1, 5), 0.31)
        assert sf_module._scalar_viscosity(arr, "e", shape=(3, 5)) == 0.31
        assert sf_module._scalar_viscosity(arr, "e") == 0.31

    def test_column_shape_passthrough_and_reject(self) -> None:
        col = np.array([[0.1, 0.2, 0.5, 0.9]])
        got = sf_module._scalar_viscosity(col, "eta1", shape=(2, 4))
        assert isinstance(got, np.ndarray) and np.array_equal(got, col)
        field = np.arange(8.0).reshape(2, 4) + 1.0
        assert np.array_equal(sf_module._scalar_viscosity(field, "eta2",
                                                          shape=(2, 4)), field)
        # 无 shape（旧调用方）⇒ 非均匀一律 raise（旧语义保留）
        with pytest.raises(ValueError):
            sf_module._scalar_viscosity(col, "eta1")
        # 给 shape 但广播不上 ⇒ 仍 raise
        with pytest.raises(ValueError):
            sf_module._scalar_viscosity(np.array([[1.0, 2.0]]), "eta1", shape=(2, 4))

    def test_solve_stream_array_m_positive_guard(self) -> None:
        geom = _tiny_geom(nz=4, ny=3)
        c = np.full((3, 4), 0.5)
        b = np.zeros((2, 3, 4))
        m_bad = np.array([[-0.5, 0.5, 0.5, 0.5]])
        with pytest.raises(ValueError, match="必须为正"):
            sf_module.solve_stream_function(geom, c, np.array([[0.2, 0.2, 0.2, 0.2]]),
                                            1.0, m_bad, b, ny=3, nz=4)


def _tiny_geom(*, nz: int, ny: int) -> dict:
    s = np.linspace(0.0, 1.0, nz)
    y = np.linspace(0.0, 1.0, ny)
    return {
        "s": s, "md": s, "y": y, "phi": y,
        "H": np.full((ny, nz), 0.05), "b": np.full((ny, nz), 0.1),
        "hole_mm": np.full((1, nz), 300.0), "od_mm": np.full((1, nz), 200.0),
    }


# --------------------------------------------------------------------------- #
# 关2⑤（F-3 代价条款）：λ_op·(1/F²) 乘积不变量在列场退化标量下逐位重钉
# --------------------------------------------------------------------------- #
class TestT6RepinDegenerateProduct:
    def _run_capture(self, *, depthwise: bool):
        solver = _solver(enable_temperature_rheology=True,
                         enable_depthwise_temperature=depthwise)
        orig = solver._froude_squared_at
        f2_log: list = []

        def spy(*args, **kwargs):
            out = orig(*args, **kwargs)
            f2_log.append(out)
            return out

        solver._froude_squared_at = spy
        res = solver.run(_well_spec(), _fluids(), _provider,
                         temperature_field=ConstantTemperatureField(60.0))
        return solver, res, f2_log

    def test_lambda_and_f2_bitwise_under_degenerate_column(self) -> None:
        s_off, r_off, f2_off = self._run_capture(depthwise=False)
        s_on, r_on, f2_on = self._run_capture(depthwise=True)
        # λ_op（域均标量口径）与逐步 F² 记录序列逐位 ⇒ 乘积不变量同位重钉
        assert s_on._lambda_op_last == s_off._lambda_op_last
        assert len(f2_on) == len(f2_off)
        assert all(a == b for a, b in zip(f2_on, f2_off))
        assert isinstance(f2_on[0], float) and isinstance(s_on._lambda_op_last, float)
        assert _run_summary(r_off) == _run_summary(r_on)


# --------------------------------------------------------------------------- #
# 关3：静温梯度列 vs 域均标量对照（深段高温 τ₀/μp 上翘显式化）
# --------------------------------------------------------------------------- #
class TestGradientDirectionAndResponse:
    def test_column_params_equal_manual_fluid_at(self) -> None:
        solver = _solver(enable_temperature_rheology=True,
                         enable_depthwise_temperature=True)
        res = solver.run(_well_spec(), _fluids(), _provider,
                         temperature_field=_StaticGradientField())
        geom = res.geom
        lead = _fluids()[2]
        b = solver._phase_props(lead, geom, 240.0)
        assert AnnulusD2DGASolver._is_col_bundle(b)
        T_col = np.asarray(b["T_col"], dtype=float)
        assert T_col.shape == geom["md"].shape
        cols = np.asarray(b["yield_stress_pa"], dtype=float).reshape(-1)
        pvs = np.asarray(b["plastic_viscosity_pa_s"], dtype=float).reshape(-1)
        for j, T in enumerate(T_col):
            ref = fluid_at(lead, float(T), None)   # P：本装配无压力场 ⇒ None（p_default 口径）
            assert cols[j] == ref.yield_stress_pa
            assert pvs[j] == ref.plastic_viscosity_pa_s

    def test_deep_hot_segment_tau0_mu_uplift(self) -> None:
        """方向断言：T>100°C 组 B 一次式上翘 ⇒ 深列（高温）τ₀/μp > 域均代表标量。"""
        solver = _solver(enable_temperature_rheology=True,
                         enable_depthwise_temperature=True)
        res = solver.run(_well_spec(), _fluids(), _provider,
                         temperature_field=_StaticGradientField())
        geom = res.geom
        lead = _fluids()[2]
        b = solver._phase_props(lead, geom, 240.0)
        T_col = np.asarray(b["T_col"], dtype=float)
        deep = int(np.argmax(geom["md"]))       # geom md 降序 ⇒ index 0 = 最深
        assert geom["md"][deep] == float(np.max(geom["md"]))
        cols = np.asarray(b["yield_stress_pa"], dtype=float).reshape(-1)
        pvs = np.asarray(b["plastic_viscosity_pa_s"], dtype=float).reshape(-1)
        rep = b["scalar"]                        # 域均代表标量（= 标量路径本步结果）
        assert T_col[deep] > float(np.mean(T_col)) > 100.0   # 断点以上区间（上翘前提）
        assert cols[deep] > rep.yield_stress_pa              # τ₀ 上翘
        assert pvs[deep] > rep.plastic_viscosity_pa_s        # μp 上翘（组 B 二次式 T>24.3 升段）

    def test_field_uses_column_props_and_dynamics_differ(self) -> None:
        """列化 ⇒ 深段流体属性 ≠ 域均属性（场公式吃列参）+ 整机响应可分。"""
        off = _solver(enable_temperature_rheology=True)
        on = _solver(enable_temperature_rheology=True, enable_depthwise_temperature=True)
        r_off = off.run(_well_spec(), _fluids(), _provider,
                        temperature_field=_StaticGradientField())
        r_on = on.run(_well_spec(), _fluids(), _provider,
                        temperature_field=_StaticGradientField())
        assert _run_summary(r_off) != _run_summary(r_on)   # Δη 如实（小量非零）
        # 场 τy：深列（j=0，md 降序）≠ 域均标量场同列
        p_off = _direct_props(off, r_off, t=240.0)
        p_on = _direct_props(on, r_on, t=240.0)
        tau_off, tau_on = p_off[3], p_on[3]
        assert np.isfinite(tau_on).all()
        assert not np.array_equal(tau_off[:, 0], tau_on[:, 0])
        # κ_mix 列化（np.maximum 位点）：列参下不抛 TypeError（隐式验证 = 本调用成功）
        assert tau_on.shape == (on.ny, on.nz)


# --------------------------------------------------------------------------- #
# 关5：memo 尺寸上界 + 审计聚合
# --------------------------------------------------------------------------- #
class TestMemoAndAuditObservations:
    def test_static_field_column_memo_bounded_per_run(self) -> None:
        """静温档（t 无关）：列派生每 run 每相一次 ⇒ 批次/尺寸不随步数增长。"""
        base = dict(enable_temperature_rheology=True,
                    enable_depthwise_temperature=True)
        s1 = _solver(**base)
        s1.run(_well_spec(), _fluids(), _provider, temperature_field=_StaticGradientField())
        s2 = _solver(**base, total_t=480.0)   # 步数翻倍
        s2.run(_well_spec(), _fluids(), _provider, temperature_field=_StaticGradientField())
        assert s2._col_batches == s1._col_batches          # 与步数无关
        assert len(s2._col_memo) == len(s1._col_memo)
        # 上界：≤ 相数×常数（链上派生对象 + 末尾诊断原始相 ⇒ 每相两键以内）
        assert len(s1._col_memo) <= 4 * 2 * 2
        assert s1._uniform_field_hits == 0                 # 梯度场不走均匀短路

    def test_columnwise_audit_ring_aggregated_not_overflow(self) -> None:
        """列批 fluid_at 事件排空聚合：环不爆（<定长 10000）+ 计数非零合理。"""
        import cemdisp.data.rheology_vs_temperature as rvt
        solver = _solver(enable_temperature_rheology=True,
                         enable_depthwise_temperature=True)
        rvt.reset_audit()
        res = solver.run(_well_spec(), _fluids(), _provider, temperature_field=_RampField())
        assert len(rvt.get_audit()) <= rvt._AUDIT_MAX      # 列×步查询打爆环 ⇒ 排空后必然有界
        agg = solver._col_audit_counts
        nz = solver.nz
        assert solver._col_batches > 0
        # 聚合语义 = 现环汇总：mud 全程出 [40,80] 域 ⇒ 每列批次 nz+1 条 clamp 起步
        assert agg.get("clamp", 0) >= nz
        assert sum(agg.values()) > 0
        assert res.summary["temperature_rheology_audit"] == solver._temp_rheo_audit_counts

    def test_uniform_field_no_column_batches(self) -> None:
        solver = _solver(enable_temperature_rheology=True,
                         enable_depthwise_temperature=True)
        solver.run(_well_spec(), _fluids(), _provider,
                   temperature_field=ConstantTemperatureField(60.0))
        assert solver._col_batches == 0
        assert solver._col_audit_counts == {}
        assert solver._uniform_field_hits > 0


# --------------------------------------------------------------------------- #
# 四站点同传：逐列 + include_yield_term=True（拆分旗标贯穿列站点，负对照不破）
# --------------------------------------------------------------------------- #
class TestColumnwiseWithYieldSplit:
    def test_eta_column_matches_site_wise_scalar_formula(self) -> None:
        """η 列 == 逐列 `fluid_apparent_viscosity(include_yield_term=True)`（同旗标）。"""
        solver = _solver(enable_temperature_rheology=True,
                         enable_depthwise_temperature=True, include_yield_term=True)
        res = solver.run(_well_spec(), _fluids(), _provider,
                         temperature_field=_StaticGradientField())
        geom = res.geom
        lead = _fluids()[2]
        b = solver._phase_props(lead, geom, 240.0)
        g = 6.0 * 0.45 / max(float(np.mean(geom["b"])), 1e-12)
        eta_col = np.asarray(solver._eta_column(b, g), dtype=float).reshape(-1)
        T_col = np.asarray(b["T_col"], dtype=float)
        for j, T in enumerate(T_col):
            ref = buoyancy.fluid_apparent_viscosity(
                fluid_at(lead, float(T), None), g, include_yield_term=True)
            assert eta_col[j] == pytest.approx(ref, rel=0.0, abs=0.0)  # 逐位
            assert 1e-5 <= eta_col[j] <= 3.0 + 0.0  # clip 结构常数域（场/标量口径同构）

    def test_scalar_sites_stay_representative(self) -> None:
        """F-3 边界：λ_op/γ̇_rep/F²/b_num 在列+拆分组合下仍为域均标量 float。"""
        solver = _solver(enable_temperature_rheology=True,
                         enable_depthwise_temperature=True, include_yield_term=True)
        orig = solver._froude_squared_at
        f2_log: list = []

        def spy(*a, **k):
            out = orig(*a, **k)
            f2_log.append(out)
            return out

        solver._froude_squared_at = spy
        solver.run(_well_spec(), _fluids(), _provider,
                   temperature_field=_StaticGradientField())
        assert all(isinstance(v, float) for v in f2_log)
        assert isinstance(solver._lambda_op_last, float)
        assert isinstance(solver._shear_rate_rep_last, float)

    def test_negative_control_split_on_off_resolvable(self) -> None:
        """负对照（复用 test_yield_term_split T6 思路）：列模式下拆分旗标可分辨
        （只传部分站点会破坏的"同传"语义 ⇒ 四站点消费同一旗标 ⇒ ON≠OFF 动力学）。"""
        kw = dict(enable_temperature_rheology=True, enable_depthwise_temperature=True)
        on = _solver(**kw, include_yield_term=True)
        off = _solver(**kw, include_yield_term=False)
        r_on = on.run(_well_spec(), _fluids(), _provider,
                      temperature_field=_StaticGradientField())
        r_off = off.run(_well_spec(), _fluids(), _provider,
                        temperature_field=_StaticGradientField())
        assert _run_summary(r_on) != _run_summary(r_off)
        # 且两档的 η 列与标量代表按各自旗标走（同一 self.include_yield_term 源）
        g = 0.5
        b = on._phase_props(_fluids()[2], r_on.geom, 200.0)
        e_on = np.asarray(on._eta_column(b, g), dtype=float).reshape(-1)
        b_off = off._phase_props(_fluids()[2], r_off.geom, 200.0)
        e_off = np.asarray(off._eta_column(b_off, g), dtype=float).reshape(-1)
        assert np.all(e_on >= e_off)      # 拆分只加不减（τy/γ̇ 项 + clip 上限内）
