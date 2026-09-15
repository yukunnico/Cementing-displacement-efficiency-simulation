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
        before = float(np.sum(g["b"] * (lead + tail + spacer)))
        after = float(np.sum(g["b"] * (lead2 + tail2 + spacer2)))
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
        q_half = float(np.sum(g["b"][:, 0] * w[0, 0])) * (g["y"][1] - g["y"][0])
        dy, ds = g["y"][1] - g["y"][0], g["s"][1] - g["s"][0]
        gained_vol = float(np.sum(g["b"] * (tail2 - tail))) * dy * ds
        assert gained_vol == pytest.approx(q_half * dt, rel=2e-3)

    def test_characteristic_speed_is_u_times_q0prime(self):
        """均匀 b + 均匀 w=U 时单步更新 = donor-cell 上风差分 −U·Δs⁻¹Δq₀(c̄)
        ——准线性速度 U·q₀'(c̄) 进入动力学（前缘 q₀'(0)=1.5 稀疏波）。"""
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
        # donor-cell（w>0 上风取左胞）：dc/dt = −U·(q₀_j − q₀_{j−1})/ds
        q0 = np.asarray(isotropic_flux_q0(c, 1.0), dtype=float)
        dc_dt_ref = -U * (q0[1, 2:-2] - q0[1, 1:-3]) / (s_arr[1] - s_arr[0])
        dc_dt_num = (lead2 - lead)[1, 2:-2] / dt
        assert np.allclose(dc_dt_num, dc_dt_ref, rtol=1e-10, atol=1e-15)

    def test_q0_switch_flag(self):
        """enable_q0_flux 默认 True、可关闭（回退开关；False 时 run() 不消费）。"""
        assert AnnulusD2DGASolver().enable_q0_flux is True
        assert AnnulusD2DGASolver(enable_q0_flux=False).enable_q0_flux is False
