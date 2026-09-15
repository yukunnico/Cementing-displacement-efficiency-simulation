# -*- coding: utf-8 -*-
"""R32 迭代 1 契约测试：新路径 (4.25) 守恒通量形式输运（2026-09-15）。

文献锚：(2.1)/(4.9) 守恒输运 + (4.25) 第一项 = 速度·q₀(c̄,m)（(4.27)/(4.28)，
Yang & Yortsos 1997 TFE——q₀'(0)=1.5 前缘稀疏波、q₀'(1)=0）。donor-cell 上风
离散，(y,s) 度量守恒律两支对称带 b。仅新路径消费（旧路径半拉格朗日由 R7
冻结锚锁定，本文件不涉及）。
"""
from __future__ import annotations

import numpy as np
import pytest

from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
from cemdisp.models2d.two_layer import isotropic_flux_q0


def _geom_ecc(e: float = 0.6, ny: int = 16, nz: int = 12, half_gap: float = 0.0458):
    phi = np.linspace(0.0, 1.0, ny)
    H = half_gap * (1.0 + e * np.cos(np.pi * phi))[:, None] * np.ones((1, nz))
    y = np.linspace(0.0, np.pi * 0.1071, ny)
    s = np.linspace(0.0, 10.0, nz)
    return {"y": y, "phi": phi, "H": H, "b": 2.0 * H, "s": s,
            "hole_mm": np.full((1, nz), 260.0), "od_mm": np.full((1, nz), 168.3)}


def _solver():
    return AnnulusD2DGASolver(nz=10, ny=8, enable_q0_flux=True)


def _trap_weights(g):
    """y 向 trapezoid 求积权重（= 方案顶点中心测度 h_i，= bulk_fill 的 np.trapezoid）。"""
    ny = g["H"].shape[0]
    wq_ = np.full(ny, g["y"][1] - g["y"][0])
    wq_[0] = wq_[-1] = wq_[0] / 2.0
    return wq_


class TestFluxFormTransport:
    """新路径 (4.25) 通量形式输运的守恒/边界/特征速度契约。"""

    def test_pure_single_fluid_states_are_steady(self):
        """c̄=1（纯水泥，入口亦水泥）与 c̄=0（纯泥浆，零入流）态输运不动。

        q₀(1)=1、q₀(0)=0：均匀单一流体态的通量场无散度（w 均匀、b 常数），
        ∂t c̄ ≡ 0——(4.25) 的端点精确性（isotropic_flux_q0 契约）在输运层保持。
        """
        s = _solver()
        g = _geom_ecc(e=0.0, ny=8, nz=10)
        w = np.full((8, 10), 0.05)
        v = np.zeros((8, 10))
        for cement_in, lead0 in ((1.0, 1.0), (0.0, 0.0)):
            lead = np.full((8, 10), lead0)
            tail = np.zeros((8, 10))
            spacer = np.zeros((8, 10))
            f_in = (lead0, 0.0, 0.0)   # 入口 declares 与场一致的单相态
            lead2, tail2, spacer2 = s._transport_flux_form(
                lead, tail, spacer, w, v, g, 0.5, 1.0, *f_in)
            assert np.array_equal(lead2, lead)
            assert np.array_equal(tail2, tail)
            assert np.array_equal(spacer2, spacer)

    def test_domain_integral_conserved_closed_boundaries(self):
        """零入流闭合域（出入口行置纯态 ⇒ 端面通量零）：Σb·(lead+tail+spacer)
        逐位守恒——donor-cell 散度形式的构造性守恒（持续造体积的消除保证）。"""
        s = _solver()
        g = _geom_ecc(e=0.6, ny=16, nz=12)
        ny, nz = g["H"].shape
        rng = np.random.default_rng(7)
        c = np.clip(0.5 + 0.4 * rng.standard_normal((ny, nz)) *
                    np.exp(-((g["phi"][:, None] - 0.5) ** 2) / 0.04), 0.0, 1.0)
        lead, tail = 0.6 * c, 0.3 * c
        spacer = 0.1 * c
        w = 0.05 * np.cos(np.pi * g["phi"])[:, None] * np.ones((1, nz))
        v = 0.02 * np.sin(np.pi * g["phi"])[:, None] * np.ones((1, nz))
        lead[:, 0] = 0.0; tail[:, 0] = 0.0; spacer[:, 0] = 0.0   # 出入口行置纯态
        lead[:, -1] = 0.0; tail[:, -1] = 0.0; spacer[:, -1] = 0.0
        zero_in = (0.0, 0.0, 0.0)
        lead2, tail2, spacer2 = s._transport_flux_form(
            lead.copy(), tail.copy(), spacer.copy(), w, v, g, 0.7, 1.5, *zero_in)
        wt = _trap_weights(g)[:, None]
        before = float(np.sum(wt * g["b"] * (lead + tail + spacer)))
        after = float(np.sum(wt * g["b"] * (lead2 + tail2 + spacer2)))
        assert after == pytest.approx(before, rel=1e-12)

    def test_inlet_flux_volume_matches_q_half(self):
        """纯水泥入口：每步流入体积 = Σb·w_inlet·dy·dt = q_half·dt（守恒律边界项）。"""
        s = _solver()
        g = _geom_ecc(e=0.6, ny=16, nz=12)
        ny, nz = g["H"].shape
        lead = np.zeros((ny, nz))
        tail = np.zeros((ny, nz))
        spacer = np.zeros((ny, nz))
        w = np.full((ny, nz), 0.04)
        v = np.zeros((ny, nz))
        dt = 0.5
        lead2, tail2, _ = s._transport_flux_form(
            lead, tail, spacer, w, v, g, dt, 1.0, 0.0, 1.0, 0.0)
        q_half = float(np.trapezoid(g["b"][:, 0] * w[:, 0], x=g["y"]))
        ds = g["s"][1] - g["s"][0]
        wt = _trap_weights(g)[:, None]
        gained_vol = float(np.sum(wt * g["b"] * (tail2 - tail))) * ds
        assert gained_vol == pytest.approx(q_half * dt, rel=2e-3)

    def test_characteristic_speed_is_u_times_q0prime(self):
        """minmod 限制器 φ∈[0,1] ⇒ 每面通量介于 donor（上风）与中心之间：单步
        更新逐点夹在上风差分与中心差分两参照之间——准线性速度 ∈ [U·q₀'的上风/
        中心实现]，光滑单调区取中心（2 阶）。"""
        s = _solver()
        g = _geom_ecc(e=0.0, ny=8, nz=200)
        ny, nz = g["H"].shape
        U = 0.04
        w = np.full((ny, nz), U)
        v = np.zeros((ny, nz))
        eps = 1e-4
        s_arr = g["s"]
        L = s_arr[-1]
        c = (0.4 + eps * np.sin(2.0 * np.pi * s_arr / L))[None, :] * np.ones((ny, 1))
        lead, tail = c.copy(), np.zeros_like(c)
        dt = 0.05
        lead2, _, _ = s._transport_flux_form(
            lead, tail, np.zeros_like(c), w, v, g, dt, 1.0, 0.0, 1.0, 0.0)
        q0 = np.asarray(isotropic_flux_q0(c, 1.0), dtype=float)
        ds = s_arr[1] - s_arr[0]
        dc_dt_num = (lead2 - lead)[1, 2:-2] / dt
        ref_center = -U * (q0[1, 3:-1] - q0[1, 1:-3]) / (2.0 * ds)   # 中心差分（2 阶参照）
        # minmod 光滑区 φ=1+O(h) ⇒ 更新 = 中心差分 + O(h²)——与中心差分偏差
        # 应为百分位以下（一阶 donor 同网格偏差同量级处）但方向由限制器定；
        # 此处断言光滑单调区近中心（2 阶行为；面通量逐面有界已由
        # test_gaussian_peak_retained_better_than_donor 判别）。
        mid = slice(70, 90)
        assert np.allclose(dc_dt_num[mid], ref_center[mid], rtol=5e-2), (
            np.max(np.abs(dc_dt_num[mid] - ref_center[mid]) /
                   np.max(np.abs(ref_center[mid]))))

    def test_q0_switch_flag(self):
        """enable_q0_flux 默认 True、可关闭（回退开关；False 时 run() 不消费）。"""
        assert AnnulusD2DGASolver().enable_q0_flux is True
        assert AnnulusD2DGASolver(enable_q0_flux=False).enable_q0_flux is False


# ---------------------------------------------------------------------------
# R33 迭代 2：FCT 反扩散（minmod TVD 通量限制器）与入口度量归一
# ---------------------------------------------------------------------------
class TestFCTAntiDiffusion:
    """minmod 通量限制器（Boris-Book/Zalesak FCT 谱系，Sweby 1984）契约。"""

    @staticmethod
    def _solver():
        return AnnulusD2DGASolver(nz=10, ny=8, enable_q0_flux=True)

    def test_step_transport_no_new_extrema(self):
        """TVD 单调保持：阶跃输运不产生新极值（minmod 的单调性保证）。"""
        s = self._solver()
        g = _geom_ecc(e=0.0, ny=8, nz=200)
        ny, nz = g["H"].shape
        w = np.full((ny, nz), 0.04)
        v = np.zeros((ny, nz))
        lead = (g["s"][None, :] > 5.0) * 1.0 * np.ones((ny, 1)) * 0.9
        tail = np.zeros_like(lead)
        c_min, c_max = float(lead.min()), float(lead.max())
        for _ in range(120):
            lead, tail, _ = s._transport_flux_form(
                lead, tail, np.zeros_like(lead), w, v, g, 0.05, 1.0, 0.0, 1.0, 0.0)
            lead = np.clip(lead, 0.0, 1.0)
        assert float(lead.min()) >= c_min - 1e-12
        assert float(lead.max()) <= c_max + 1e-12

    def test_gaussian_peak_retained_better_than_donor(self):
        """FCT 反扩散判别（高斯脉冲峰值保持）：同参数平流 N 步后 FCT 的峰值
        显著高于一阶 donor-cell 参照（内联独立复算），且接近初值。"""
        g = _geom_ecc(e=0.0, ny=8, nz=400)
        ny, nz = g["H"].shape
        b = g["b"]
        w = np.full((ny, nz), 0.04)
        v = np.zeros((ny, nz))
        s_arr = g["s"]
        ds = s_arr[1] - s_arr[0]
        dt = 0.05
        n_steps = 100
        lead0 = (0.9 * np.exp(-0.5 * ((s_arr[None, :] - 2.0) / 0.15) ** 2)
                 * np.ones((ny, 1)))

        # 纯 donor-cell 参照（独立复算，一阶上风；定性 = 迭代 1 的输运执行级）
        c_donor = lead0.copy()
        for _ in range(n_steps):
            w_face = np.zeros((ny, nz + 1))
            w_face[:, 1:nz] = 0.5 * (w[:, :-1] + w[:, 1:])
            w_face[:, 0] = w[:, 0]
            w_face[:, nz] = w[:, -1]
            Gs = np.zeros((ny, nz + 1))
            Gs[:, 1:nz] = np.where(w_face[:, 1:nz] >= 0.0, c_donor[:, :-1], c_donor[:, 1:])
            Gs[:, 0] = c_donor[:, 0]
            Gs[:, nz] = c_donor[:, -1]
            dc = -dt * ((w_face[:, 1:] * Gs[:, 1:]) - (w_face[:, :-1] * Gs[:, :-1])) / (ds * b)
            c_donor = np.clip(c_donor + dc, 0.0, 1.0)

        s = _solver()
        lead = lead0.copy()
        for _ in range(n_steps):
            lead, _, _ = s._transport_flux_form(
                lead, np.zeros_like(lead), np.zeros_like(lead), w, v, g, dt, 1.0,
                0.0, 1.0, 0.0)
            lead = np.clip(lead, 0.0, 1.0)

        peak_fct = float(lead.max())
        peak_donor = float(c_donor.max())
        peak_init = float(lead0.max())
        assert peak_fct > peak_donor + 0.05, (peak_fct, peak_donor)
        assert peak_fct > 0.95 * peak_init, (peak_fct, peak_init)


class TestInletMetricNormalization:
    """入口通量度量归一（R33 迭代 2，case1 缺陷修复）契约。"""

    def test_spiky_w_inlet_volume_still_q_half(self):
        """欠解析 w 尖峰（case1 缺陷回归）：入口流入体积仍 = q_half·dt。"""
        s = AnnulusD2DGASolver(nz=10, ny=8, enable_q0_flux=True)
        g = _geom_ecc(e=0.6, ny=16, nz=12)
        ny, nz = g["H"].shape
        # 宽边（φ=0 第一行）放置 w 尖峰：矩形采样会被半格主导
        w = 0.005 * np.ones((ny, nz))
        w[0, :] = 0.5   # 尖峰：单行 = 其余行 100 倍
        v = np.zeros((ny, nz))
        lead = np.zeros((ny, nz))
        tail = np.zeros((ny, nz))
        spacer = np.zeros((ny, nz))
        dt = 0.1
        lead2, tail2, _ = s._transport_flux_form(
            lead, tail, spacer, w, v, g, dt, 1.0, 0.0, 1.0, 0.0)
        q_half = float(np.trapezoid(g["b"][:, 0] * w[:, 0], x=g["y"]))
        ds = g["s"][1] - g["s"][0]
        wt = _trap_weights(g)[:, None]
        gained_vol = float(np.sum(wt * g["b"] * (tail2 - tail))) * ds
        assert gained_vol == pytest.approx(q_half * dt, rel=2e-3)
