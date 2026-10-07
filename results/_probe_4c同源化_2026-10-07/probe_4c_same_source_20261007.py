# -*- coding: utf-8 -*-
"""Phase 4c 验收 probe：呼101 静温档 1D/2D 鞋口同源温差对账（2026-10-07）。

操作化口径（spec §2「4a+4c 验收」/ 测绘项12）：

A. **同源差恒 0**：4a 后 1D 鞋口消费点（重力入口 = 现成鞋深+到达时刻）实际吃到的
   温度（`fluid_at` spy 捕获）vs **同一场对象**在 (shoe_md, t) 的求值 ⇒ 差恒 0；
   三方对照 `AnnulusD2DGASolver._temperature_array` 在 (shoe_md, t) 的 2D 列值
   （同场对象 ⇒ 三方全等）。
B. **旧口径对照**（证明断裂语义已移位、非被"归一化"掉）：
   1D 4a 前 = 域顶 t=0 单点标量（T(top,0)≈111.04°C）；
   2D 逐步代表 = 深度均值（F-3 维持原样，域均口径不动）≈132.75°C；
   差 = 0.017598×Δmd ≈ 21.72°C。4a 后 1D 主消费点集合不含域顶 ⇒ 断裂不再发生；
   域顶值仅存在于回退标量 `_step_T_c`（无 md 直调路径）。
C. **关5 落盘**：弥散入口 (md,t) 体积链反推记录与闭合残差（管容/等效截面 − 鞋深）。

运行：
  PYTHONIOENCODING=utf-8 PYTHONUTF8=1 "D:/apps/Anaconda/envs/cementT/python.exe" \
      probe_4c_same_source_20261007.py
输出：probe_4c_same_source_20261007.json（同目录）
"""
from __future__ import annotations

import json
from pathlib import Path

from cemdisp.data.loaders import load_hu101_tailpipe
from cemdisp.data.temperature_field import (
    GEO_GRAD_C_PER_M,
    GeothermalTemperatureField,
)
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver
from cemdisp.transport1d import casing_flow as casing_module
from cemdisp.transport1d.casing_flow import CasingFlowSolver

HERE = Path(__file__).resolve().parent
OUT = HERE / "probe_4c_same_source_20261007.json"

# --------------------------------------------------------------------------- #
# 装配：呼101 静温档（无瞬态表井 ⇒ Geothermal 线，09-30 §1 统一式）
# --------------------------------------------------------------------------- #
well_spec, fluids, schedule, _validation = load_hu101_tailpipe()
field = GeothermalTemperatureField()

# casing 生产口径与 runners/hu101_tailpipe.py 一致（T1 裁定双开关全开 + has_plug），
# 叠加 enable_temperature_rheology=True（4a/4c 验收为 T-on 语义；runner 本身不动）。
solver = CasingFlowSolver(
    enable_gravity=True,
    mixing_contact_time=True,
    plug_face_zero_mixing=True,
    has_plug=True,
    enable_temperature_rheology=True,
)

calls: list[tuple[str, float | None, float | None]] = []
_orig_phase = solver._phase_props


def _phase_spy(fluid, md_m=None, t_s=None):
    calls.append((fluid.name, md_m, t_s))
    return _orig_phase(fluid, md_m, t_s)


solver._phase_props = _phase_spy
result = solver.run(well_spec, fluids, schedule, temperature_field=field)

shoe_m = float(well_spec.shoe_md_m)
top_m = float(well_spec.top_md_m)
bottom_m = float(well_spec.bottom_md_m)

# --------------------------------------------------------------------------- #
# A. 鞋口同源差：1D 实际消费温度 vs 同一场对象在 (shoe, t) 求值（差恒 0）
# --------------------------------------------------------------------------- #
tail = next(f for f in fluids if f.name == "尾浆")
shoe_times = sorted({t for _n, md, t in calls
                     if md is not None and md == shoe_m and t is not None})
_orig_fa = casing_module.fluid_at


def _capture_T(md: float, t: float) -> float:
    """清 memo 后单点求值，捕获 `fluid_at` 实际吃到的温度。"""
    solver._phase_memo.clear()
    log: list[tuple] = []

    def _fa(fluid, temperature_c, pressure_mpa=None, **kw):
        log.append((fluid.name, temperature_c, pressure_mpa))
        return _orig_fa(fluid, temperature_c, pressure_mpa, **kw)

    casing_module.fluid_at = _fa
    try:
        solver._phase_props(tail, md, t)
    finally:
        casing_module.fluid_at = _orig_fa
    assert log, "单点求值必有 fluid_at 消费"
    return log[0][1]


ann = AnnulusD2DGASolver(dt=2.0, nz=8, ny=4, total_t=10.0)
ann._temperature_field = field

a_rows = []
for t_q in shoe_times[:8]:
    t_1d = _capture_T(shoe_m, t_q)
    t_field = field.T(shoe_m, t_q)
    t_2d_col = float(ann._temperature_array([shoe_m], t_q)[0])  # 2D 同 md 列值
    a_rows.append({
        "t_s": t_q,
        "T_1d_shoe_c": t_1d,
        "T_field_shoe_c": t_field,
        "T_2d_column_shoe_c": t_2d_col,
        "same_source_diff_1d_vs_field_c": t_1d - t_field,
        "same_source_diff_1d_vs_2d_column_c": t_1d - t_2d_col,
    })

all_zero_vs_field = bool(a_rows) and all(
    r["same_source_diff_1d_vs_field_c"] == 0.0 for r in a_rows)
all_zero_vs_2d = bool(a_rows) and all(
    r["same_source_diff_1d_vs_2d_column_c"] == 0.0 for r in a_rows)

# --------------------------------------------------------------------------- #
# B. 旧口径对照（21.7°C 断裂算术，测绘项12 复现）+ 移位证明
# --------------------------------------------------------------------------- #
T_top_c = field.T(top_m, 0.0)                 # 1D 4a 前口径：域顶 t=0 单点
T_mean_c = field.T(0.5 * (top_m + bottom_m), 0.0)  # 2D 逐步口径：深度均值（F-3 维持）
old_fracture_c = T_mean_c - T_top_c           # ≈21.72
T_shoe_c = field.T(shoe_m, 0.0)               # 4a 后 1D 鞋口消费口径
used_md = {md for _n, md, _t in calls if md is not None}
domaintop_consumed = top_m in used_md         # 必须 False ⇒ 断裂语义已移位
fallback_scalar = solver._step_T_c            # 域顶值仍在（仅回退路径）
assert fallback_scalar == T_top_c             # 回退语义 = 改前，逐位保留

# --------------------------------------------------------------------------- #
# C. 关5：弥散入口 (md,t) 体积链反推记录 + 闭合残差落盘
# --------------------------------------------------------------------------- #
recs = solver._dispersion_depthwise_records
closure = {
    "n_dispersion_events": len(recs),
    "closure_residual_m_min": min((r[3] for r in recs), default=None),
    "closure_residual_m_max": max((r[3] for r in recs), default=None),
    "sample": [
        {"t_event_s": r[0], "md_pathwise_m": r[1], "t_pathwise_s": r[2],
         "closure_residual_m": r[3]}
        for r in recs[:12]
    ],
}

payload = {
    "well": "hu101（呼101）静温档 GeothermalTemperatureField",
    "domain": {"top_md_m": top_m, "bottom_md_m": bottom_m, "shoe_md_m": shoe_m},
    "geo_grad_c_per_m": GEO_GRAD_C_PER_M,
    "A_same_source_shoe": {
        "rows": a_rows,
        "all_zero_vs_field": all_zero_vs_field,
        "all_zero_vs_2d_column": all_zero_vs_2d,
    },
    "B_old_caliber_contrast": {
        "T_1d_domaintop_fallback_c": T_top_c,
        "T_2d_domainmean_representative_c": T_mean_c,
        "old_fracture_c": old_fracture_c,
        "T_1d_shoe_depthwise_c": T_shoe_c,
        "new_gap_shoe_minus_domaintop_c": T_shoe_c - T_top_c,
        "domaintop_in_primary_consumption": domaintop_consumed,
        "note": "2D 域均代表标量按 F-3 原样保留（不『归一化』）；"
                "同源=同场同点，鞋口三方差恒 0；域顶仅回退路径。",
    },
    "C_dispersion_volume_chain": closure,
    "summary": {
        "shoe_same_source_diff_zero": all_zero_vs_field and all_zero_vs_2d,
        "old_fracture_reproduced": abs(old_fracture_c - 0.017598 * (bottom_m - top_m) / 2.0) < 1e-9,
        "fracture_semantics_shifted": not domaintop_consumed,
    },
}

OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

print("=== 4c 同源 probe（呼101 静温档）===")
print(f"域: top={top_m} / shoe={shoe_m} / bottom={bottom_m}")
print(f"[A] 鞋口同源差: vs 场求值 max|diff| = "
      f"{max(abs(r['same_source_diff_1d_vs_field_c']) for r in a_rows)} "
      f"(n={len(a_rows)}, 恒0={all_zero_vs_field}); "
      f"vs 2D列值 恒0={all_zero_vs_2d}")
for r in a_rows[:4]:
    print(f"    t={r['t_s']:10.3f}s  T_1d={r['T_1d_shoe_c']:.6f}°C  "
          f"T_field={r['T_field_shoe_c']:.6f}°C  Δ={r['same_source_diff_1d_vs_field_c']}")
print(f"[B] 旧口径对照: 域顶 {T_top_c:.4f}°C vs 域均 {T_mean_c:.4f}°C "
      f"⇒ 断裂 {old_fracture_c:.4f}°C（测绘项12 ≈21.72）")
print(f"    新口径 1D 鞋口 = {T_shoe_c:.4f}°C；域顶进入主消费集合？{domaintop_consumed} "
      f"⇒ 断裂语义移位（域顶仅回退 _step_T_c={fallback_scalar:.4f}）")
print(f"[C] 弥散沿程记录 n={closure['n_dispersion_events']} "
      f"闭合残差范围 [{closure['closure_residual_m_min']}, {closure['closure_residual_m_max']}] m")
ok = (payload["summary"]["shoe_same_source_diff_zero"]
      and payload["summary"]["old_fracture_reproduced"]
      and payload["summary"]["fracture_semantics_shifted"])
print(f"结论: {'PASS' if ok else 'FAIL'}（JSON → {OUT.name}）")
raise SystemExit(0 if ok else 1)
