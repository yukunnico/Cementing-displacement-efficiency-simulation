"""4b 改前 HEAD 锚捕获（2026-10-07，短算例逐位锚=测试关2① 数据源）。

HEAD（Phase 4b 改动前）跑三类短算例并把关键浮点/字节哈希写入 head_anchors.json，
tests/contract/test_annulus_columnwise_temperature.py 以字面量内嵌复跑断言逐位。
运行（改前一次、改后校验一次）：
  PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python capture_head_anchors_20261007.py [--out <json>]
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

BRANCH_ROOT = Path(__file__).resolve().parents[2]
import sys
if str(BRANCH_ROOT) not in sys.path:
    sys.path.insert(0, str(BRANCH_ROOT))
import cemdisp  # noqa: E402
assert Path(cemdisp.__file__).resolve().is_relative_to(BRANCH_ROOT), \
    f"cemdisp 指向错误树: {cemdisp.__file__}"

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel  # noqa: E402
from cemdisp.data.rheology_vs_temperature import fluid_at  # noqa: E402
from cemdisp.data.temperature_field import ConstantTemperatureField  # noqa: E402
from cemdisp.data.well_spec import DepthValuePoint, WellSpec  # noqa: E402
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver  # noqa: E402
from cemdisp.models2d.boundary_bridge import AnnulusInletState  # noqa: E402

_Q_M3S = 1.0 / 60.0  # 生产级排量（同 test_temperature_phase_props）


def _well_spec() -> WellSpec:
    top, bottom = 100.0, 400.0
    return WellSpec(
        well_name="P4b_columnwise_2d",
        top_md_m=top,
        bottom_md_m=bottom,
        shoe_md_m=bottom,
        hole_diameter_profile=(DepthValuePoint(top, 260.0), DepthValuePoint(bottom, 250.0)),
        liner_od_profile=(DepthValuePoint(top, 168.3), DepthValuePoint(bottom, 168.3)),
        inclination_profile=(DepthValuePoint(top, 2.0), DepthValuePoint(bottom, 5.0)),
        standoff_profile=(DepthValuePoint(top, 0.75), DepthValuePoint(bottom, 0.60)),
    )


def _fluids() -> tuple:
    """Bingham 泥浆/隔离液 + POWER_LAW 水泥（同 T1-2 契约测试装配）。"""
    mud = FluidSpec("mud", FluidRole.MUD, 1050.0, RheologyModel.BINGHAM,
                    plastic_viscosity_pa_s=0.022, yield_stress_pa=6.0)
    spacer = FluidSpec("spacer", FluidRole.SPACER, 1120.0, RheologyModel.BINGHAM,
                       plastic_viscosity_pa_s=0.030, yield_stress_pa=9.0)
    lead = FluidSpec("lead", FluidRole.LEAD, 1890.0, RheologyModel.POWER_LAW,
                     power_law_n=0.75, consistency_k=0.55)
    tail = FluidSpec("tail", FluidRole.TAIL, 1920.0, RheologyModel.POWER_LAW,
                     power_law_n=0.70, consistency_k=0.90)
    return mud, spacer, lead, tail


class _RampField:
    """斜坡温度场（非均匀 ⇒ 关2 锚含"代表温度=深度均值"分支）：10 + 0.01·md + 0.001·t。"""

    def T(self, md_m: float, t_s: float) -> float:
        return 10.0 + 0.01 * float(md_m) + 0.001 * float(t_s)


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


def _digest(arr) -> str:
    return hashlib.sha256(np.asarray(arr, dtype=float).tobytes()).hexdigest()


def _run_case(name: str, *, temperature_rheology: bool, field) -> dict:
    solver = _solver(enable_temperature_rheology=temperature_rheology)
    res = solver.run(_well_spec(), _fluids(), _provider, temperature_field=field)
    final = res.summary["最终结果"]
    out = {
        "eta_E": float(final["全井段最终有效顶替效率"]),
        "eta_N": float(final["窄四分位效率"]),
        "buoyancy_b": float(final["浮力数_b"]),
        "cement_sha256": _digest(res.cement_field),
        "mu_sha256": _digest(np.asarray(res.geom["md"])),  # geom md 序锚（列对齐校验用）
        "step_T_c": None if final is None or not temperature_rheology
        else float(solver._step_T_c),
    }
    # 代表温度 + _compute_props 末步直调锚（场公式吃标量参=HEAD 语义）
    if temperature_rheology:
        geom = res.geom
        rep = solver._representative_temperature(geom, 240.0)
        out["rep_T"] = float(rep)
        mud, spacer, lead, tail = _fluids()
        w_prev = np.full((solver.ny, solver.nz), 0.45)
        lead0 = np.zeros((solver.ny, solver.nz))
        lead0[:, 0] = 0.5
        tail0 = np.zeros_like(lead0)
        sp0 = np.zeros_like(lead0)
        props = solver._compute_props(lead0, tail0, sp0, w_prev, geom,
                                      mud, lead, tail, spacer, t=240.0)
        out["compute_props_sha256"] = [_digest(a) for a in props]
        # 派生四相参数锚（_phase_props 直调，memo 清后）
        solver._phase_memo.clear()
        der = {n: solver._phase_props(f, geom, 240.0) for n, f in
               (("mud", mud), ("lead", lead), ("tail", tail), ("spacer", spacer))}
        out["phase_props"] = {
            n: {"PV": None if d.plastic_viscosity_pa_s is None else float(d.plastic_viscosity_pa_s),
                "TY": None if d.yield_stress_pa is None else float(d.yield_stress_pa),
                "N": None if d.power_law_n is None else float(d.power_law_n),
                "K": None if d.consistency_k is None else float(d.consistency_k)}
            for n, d in der.items()}
    return {"case": name, **out}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent / "head_anchors.json"))
    args = ap.parse_args()
    cases = [
        _run_case("toff_no_field", temperature_rheology=False, field=None),
        _run_case("ton_constant60", temperature_rheology=True,
                  field=ConstantTemperatureField(60.0)),
        _run_case("ton_ramp", temperature_rheology=True, field=_RampField()),
    ]
    Path(args.out).write_text(json.dumps(cases, indent=1), encoding="utf-8")
    print(f"anchors → {args.out}")
    for c in cases:
        print(f"{c['case']}: eta_E={c['eta_E']!r} eta_N={c['eta_N']!r} b={c['buoyancy_b']!r}")


if __name__ == "__main__":
    main()
