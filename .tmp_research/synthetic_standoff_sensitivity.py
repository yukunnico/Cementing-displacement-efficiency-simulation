# -*- coding: utf-8 -*-
"""环空 2D 模型居中度敏感性合成算例（只读取证，不改包代码）。

设计：同一 toy 井（域长 100m，Bingham 泥浆+水泥，恒定排量），只改 standoff：
  A: standoff=0.80 (e=0.20，接近 ht1_003/ht1_004)
  B: standoff=0.65 (e=0.35，hu1)
  C: standoff=0.45 (e=0.55，hu101 上限)
  D: standoff=0.45 + 关屈服门（enable_yield_gate=False，隔离 B2 屈服门贡献）
库存比设计为 <1（注入体积 ≈ 0.8×物理环空体积），保证前缘不出域（或仅宽边微出），
使 η_E/η_N 的窄边亏空能直接体现居中度的影响。
"""
from __future__ import annotations

import io
import sys
import contextlib

import numpy as np

from cemdisp.data.fluid_spec import FluidSpec, FluidRole, RheologyModel
from cemdisp.data.well_spec import DepthValuePoint, EvaluationWindow, WellSpec
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver, _trapez2d
from cemdisp.models2d.boundary_bridge import AnnulusInletState


def make_well(standoff: float) -> WellSpec:
    pts = lambda d, v: DepthValuePoint(depth_md_m=d, value=v)
    return WellSpec(
        well_name=f"synthetic_so{standoff:.2f}",
        top_md_m=1000.0, bottom_md_m=1100.0, shoe_md_m=1100.0, hanger_md_m=1000.0,
        casing_id_mm=200.0, liner_od_mm=139.7, liner_id_mm=108.0,
        hole_diameter_profile=[pts(1000.0, 215.9), pts(1100.0, 215.9)],
        inclination_profile=[pts(1000.0, 3.0), pts(1100.0, 3.0)],
        standoff_profile=[pts(1000.0, standoff), pts(1100.0, standoff)],
        evaluation_windows=[EvaluationWindow(name="w", top_md_m=1000.0, bottom_md_m=1100.0, window_type="full")],
    )


MUD = FluidSpec(name="mud", role=FluidRole.MUD, density_kg_m3=1900.0,
                rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.053, yield_stress_pa=8.5)
LEAD = FluidSpec(name="lead", role=FluidRole.LEAD, density_kg_m3=1930.0,
                 rheology_model=RheologyModel.BINGHAM, plastic_viscosity_pa_s=0.180, yield_stress_pa=14.0)
FLUIDS = (MUD, LEAD)

Q = 0.010  # m3/s = 600 L/min


def inlet(t: float) -> AnnulusInletState:
    return AnnulusInletState(
        time_s=t, flow_rate_m3_s=Q, stage_name="cement",
        phase_fractions=(("cement", 1.0), ("lead", 1.0)),
    )


def run_case(tag: str, standoff: float, **solver_kw) -> dict:
    well = make_well(standoff)
    tmp = AnnulusD2DGASolver()
    v_ann = tmp._physical_annular_volume(well)
    # 库存比目标 ~0.80：注入体积 = 0.8 * v_ann，时间 = 体积 / Q
    total_t = 0.80 * v_ann / Q
    # 生产口径：默认参数（B2 屈服门开、CFL 自适应开、弥散 0.018/0.015、e_clip_max=0.55）
    solver = AnnulusD2DGASolver(total_t=total_t, nz=100, ny=40, save_interval=10_000, **solver_kw)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        res = solver.run(well, FLUIDS, inlet)
    g = res.geom
    fin = res.metrics.iloc[-1]
    cement = res.cement_field
    s_max = float(g["s"][-1])
    v_dom = 2.0 * _trapez2d(g["b"] * cement, g)
    # 窄边（最后 ny//4 行）b 加权平均水泥浓度
    n_q = max(1, solver.ny // 4)
    b_q, c_q, y_q = g["b"][-n_q:, :], cement[-n_q:, :], g["y"][-n_q:]
    geo_q = {"y": y_q, "s": g["s"]}
    narrow_fill = _trapez2d(b_q * c_q, geo_q) / max(_trapez2d(b_q, geo_q), 1e-12)
    # 域内水泥库存比（相对物理环空体积）
    out = {
        "case": tag,
        "standoff": standoff,
        "e_clip_max": solver.e_clip_max,
        "库存比(注入/环空)": round(0.80, 3),
        "total_t_s": round(total_t, 1),
        "eta_E": round(float(fin["effective_efficiency"]), 4),
        "eta_N": round(float(res.summary["eta_narrow"]), 4),
        "窄边b加权填充": round(float(narrow_fill), 4),
        "宽边前缘m": round(float(fin["front_wide_m"]), 1),
        "窄边前缘m": round(float(fin["front_narrow_m"]), 1),
        "域长m": round(s_max, 1),
        "窜槽指数": round(float(fin["channeling_index"]), 4),
        "壁面冻结占比": round(float(fin["mean_wall_mud"]), 4),
        "混浆指数": round(float(fin["mixing_index"]), 4),
        "域内水泥/环空": round(v_dom / v_ann, 4),
    }
    return out


def main() -> None:
    cases = [
        ("A_e0.20", 0.80, {}),
        ("B_e0.35", 0.65, {}),
        ("C_e0.55", 0.45, {}),
        ("D_e0.55_无屈服门", 0.45, dict(enable_yield_gate=False)),
        ("E_e0.55_无弥散", 0.45, dict(dispersion_axial=0.0, dispersion_azimuthal=0.0)),
    ]
    rows = []
    for tag, so, kw in cases:
        row = run_case(tag, so, **kw)
        rows.append(row)
        print(row, flush=True)
    # 汇总表
    keys = list(rows[0].keys())
    print("\n" + "\t".join(keys))
    for r in rows:
        print("\t".join(str(r[k]) for k in keys))


if __name__ == "__main__":
    main()
