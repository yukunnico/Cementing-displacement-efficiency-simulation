"""2026-09-08 路线 B（管内浓度剖面真解）Task 3：8 井双开关验证 + alpha 消融前置。

复用 scripts/reruns/rerun8_three_fixes_20260906.py 的**纯默认 solver 口径**
（CasingFlowSolver(enable_gravity=True) + split_cement_phases provider +
AnnulusD2DGASolver(total_t=stop, nz=250, ny=40)），勿用 _sensitivity_common 的
CORRECTED_KW 骨架（口径不同）。变体只在 CasingFlowSolver 构造参数上开/关。

变体矩阵：8 井 × 4 变体 {基线双关, contact_time开, plug_face开, 双开}。
hu101 只跑名义口径（不跑 0906 脚本的实测 standoff 变体）。

alpha 消融前置：hu103/hu102 × dispersion_alpha∈{0.04,0.12,0.25,0.50}
（基线开关口径），量化"混浆带宽→η"灵敏度上界——alpha 减半 η 都不动
<0.2pp，则路线 B 的 η 效应预期同样微小，论文叙事按"混浆带物理化修正"写。

验收红线：双开 vs 基线——η_E 全井 |Δ|<0.5pp、stop 逐位不变；
违反则脚本打印 [BLOCKED] 并以退出码 3 结束（数字已全部落盘）。

注意：8 井 loader 生产口径 has_plug=False，plug_face_zero_mixing 判据第二项
（has_plug 门）恒假 → plug_face开 变体预期与基线逐位一致（test_plug_face_off_bitwise
已证）；本矩阵如实量化该口径下的效应。

输出：results/浓度剖面真解路线B_2026-09-08/
  - 汇总.csv                    （32 行矩阵）
  - {井}_变体对照.csv           （每井 4 变体行，增量覆盖写）
  - alpha消融_hu103_hu102.csv   （8 行，增量覆盖写）
  - 红线核对.json               （双开 vs 基线逐井判定 + alpha=0.25 确定性锚）
幂等：输出目录已存在时覆盖重写。

用法：
  python scripts/reruns/rerun8_mixing_contact_time_20260908.py smoke   # hu2 单变体冒烟
  python scripts/reruns/rerun8_mixing_contact_time_20260908.py full    # 全量（默认）
"""
from __future__ import annotations

import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver, _trapez2d  # noqa: E402
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider  # noqa: E402
from cemdisp.transport1d import CasingFlowSolver  # noqa: E402
from scripts.lib.mass_balance_diag import WELLS, integrate_injection  # noqa: E402

OUT = PROJECT_ROOT / "results" / "浓度剖面真解路线B_2026-09-08"
OUT.mkdir(parents=True, exist_ok=True)
NZ = 250
NY = 40

# 变体矩阵：只在路线 B 两个开关上开/关，其余构造参数与 0906 run_one 完全一致
VARIANTS: list[tuple[str, dict]] = [
    ("基线双关", {}),
    ("contact_time开", {"mixing_contact_time": True}),
    ("plug_face开", {"plug_face_zero_mixing": True}),
    ("双开", {"mixing_contact_time": True, "plug_face_zero_mixing": True}),
]

# alpha 消融（基线开关口径）：0.25 为生产默认，兼作确定性锚
ALPHAS = (0.04, 0.12, 0.25, 0.50)
ABLATION_WELLS = ("hu103", "hu102")


def run_one(name: str, variant: str, casing_kw: dict) -> tuple[dict, float]:
    """跑单井单变体，行字段与 0906 run_one 对齐（加"变体"列）；返回 (行, 原值stop)。"""
    t0 = time.time()
    well_spec, fluids, schedule, _ = WELLS[name]()
    casing_solver = CasingFlowSolver(enable_gravity=True, **casing_kw)
    casing_result = casing_solver.run(well_spec, fluids, schedule)
    provider = build_coupled_annulus_inlet_provider(
        casing_result, casing_solver, fluids, split_cement_phases=True)
    stop = float(casing_result.cement_end_time_s)
    inj = integrate_injection(provider, stop)

    solver = AnnulusD2DGASolver(total_t=stop, nz=NZ, ny=NY)
    res = solver.run(well_spec, fluids, provider, schedule=schedule)
    fin = res.metrics.iloc[-1]
    g = res.geom
    v_dom = 2.0 * _trapez2d(g["b"] * res.cement_field, g)
    v_full = 2.0 * _trapez2d(g["b"], g)
    s_max = float(g["s"][-1])
    sp_f = solver._pick_fluids(fluids)[3]
    row = {
        "井": name,
        "变体": variant,
        "η_E": round(float(fin["effective_efficiency"]), 4),
        "η_N": round(float(res.summary.get("eta_narrow", np.nan)), 4),
        "stop_s": round(stop, 1),
        "入环空水泥_m3": round(inj["cement"], 2),
        "域内水泥_m3": round(v_dom, 2),
        "域满体积_m3": round(v_full, 2),
        "守恒率": round(v_dom / max(inj["cement"], 1e-9), 4),
        "宽边前缘m": round(float(fin["front_wide_m"]), 1),
        "窄边前缘m": round(float(fin["front_narrow_m"]), 1),
        "窄边到位率": round(float(fin["front_narrow_m"]) / s_max, 4),
        "混浆指数": round(float(fin["mixing_index"]), 4),
        "壁面冻结占比": round(float(fin["mean_wall_mud"]), 4),
        "等效隔离液密度": (sp_f.density_kg_m3 if sp_f is not None else None),
        "等效隔离液名": (sp_f.name if sp_f is not None else None),
        "standoff实测标记": bool(getattr(well_spec, "standoff_measured", False)),
        "耗时s": round(time.time() - t0, 1),
    }
    print(json.dumps(row, ensure_ascii=False), flush=True)
    return row, stop


def _error_row(name: str, variant: str, exc: Exception) -> dict:
    return {"井": name, "变体": variant, "状态": "error",
            "错误": f"{type(exc).__name__}: {exc}"}


def run_ablation() -> list[dict]:
    """hu103/hu102 × alpha∈{0.04,0.12,0.25,0.50}（基线开关口径），增量落盘。"""
    rows: list[dict] = []
    for name in ABLATION_WELLS:
        for alpha in ALPHAS:
            variant = f"alpha={alpha}"
            try:
                row, _ = run_one(name, variant, {"dispersion_alpha": alpha})
                row["dispersion_alpha"] = alpha
            except Exception as exc:
                traceback.print_exc()
                row = _error_row(name, variant, exc)
                row["dispersion_alpha"] = alpha
            rows.append(row)
            pd.DataFrame(rows).to_csv(
                OUT / "alpha消融_hu103_hu102.csv", index=False, encoding="utf-8-sig")
    return rows


def check_red_line(per_well: dict[str, list[dict]],
                   stop_raw: dict[tuple[str, str], float]) -> bool:
    """双开 vs 基线：η_E 全井 |Δ|<0.5pp、stop 逐位不变。返回是否全部通过。"""
    print("\n===== 验收红线核对：双开 vs 基线 =====")
    verdict: dict = {}
    all_pass = True
    for name, rows in per_well.items():
        base = next((r for r in rows if r.get("变体") == "基线双关" and "η_E" in r), None)
        both = next((r for r in rows if r.get("变体") == "双开" and "η_E" in r), None)
        if base is None or both is None:
            verdict[name] = {"判定": "FAIL", "原因": "基线或双开缺行（该井存在 error）"}
            all_pass = False
            print(f"{name}: 缺行 → FAIL")
            continue
        d_eta_pp = (both["η_E"] - base["η_E"]) * 100.0
        s_base = stop_raw.get((name, "基线双关"))
        s_both = stop_raw.get((name, "双开"))
        stop_same = (s_base is not None and s_base == s_both)
        ok = (abs(d_eta_pp) < 0.5) and stop_same
        verdict[name] = {
            "η_E_基线": base["η_E"], "η_E_双开": both["η_E"],
            "Δη_E_pp": round(d_eta_pp, 4),
            "stop_基线_s": s_base, "stop_双开_s": s_both,
            "stop逐位一致": stop_same,
            "判定": "PASS" if ok else "FAIL",
        }
        if not ok:
            all_pass = False
        print(f"{name}: Δη_E={d_eta_pp:+.4f}pp, stop逐位一致={stop_same} "
              f"→ {'PASS' if ok else 'FAIL'}")

    # alpha=0.25 确定性锚：消融表 0.25 行应与矩阵基线行逐位一致（同口径同默认）
    abl_path = OUT / "alpha消融_hu103_hu102.csv"
    if abl_path.exists():
        abl = pd.read_csv(abl_path)
        anchors = {}
        for name in ABLATION_WELLS:
            r025 = abl[(abl["井"] == name) & (abl["变体"] == "alpha=0.25")]
            base = per_well.get(name, [])
            b = next((r for r in base if r.get("变体") == "基线双关" and "η_E" in r), None)
            if not r025.empty and b is not None and "η_E" in r025.iloc[0]:
                same = bool(
                    r025.iloc[0]["η_E"] == b["η_E"]
                    and abs(float(r025.iloc[0]["stop_s"]) - float(b["stop_s"])) < 1e-6
                )
                anchors[name] = {"alpha0.25==基线": same,
                                 "判定": "PASS" if same else "FAIL"}
                if not same:
                    all_pass = False
        verdict["确定性锚_alpha0.25vs基线"] = anchors
        print(f"确定性锚 alpha=0.25 vs 基线: {anchors}")

    (OUT / "红线核对.json").write_text(
        json.dumps(verdict, ensure_ascii=False, indent=2), encoding="utf-8")
    if not all_pass:
        print("[BLOCKED] 验收红线未通过——报告标 BLOCKED（数字已落盘）")
        return False
    print("[PASS] 验收红线全部通过")
    return True


def main_full() -> int:
    print("===== 阶段A：alpha 消融（hu103/hu102 × 4 值） =====", flush=True)
    run_ablation()
    print("\n===== 阶段B：8 井双开关矩阵 =====", flush=True)
    per_well, stop_raw = run_matrix_named()
    ok = check_red_line(per_well, stop_raw)
    print("\n[done] 全部落盘完成" if ok else "\n[done] 全部落盘完成（红线 FAIL）")
    return 0 if ok else 3


def run_matrix_named() -> tuple[dict[str, list[dict]], dict[tuple[str, str], float]]:
    """同 run_matrix，但以 {井: 行列表} 返回，供红线核对。"""
    per_well: dict[str, list[dict]] = {}
    stop_raw: dict[tuple[str, str], float] = {}
    rows: list[dict] = []
    for name in WELLS:
        per_well.setdefault(name, [])
        for variant, kw in VARIANTS:
            try:
                row, stop = run_one(name, variant, kw)
                stop_raw[(name, variant)] = stop
            except Exception as exc:
                traceback.print_exc()
                row = _error_row(name, variant, exc)
            rows.append(row)
            per_well[name].append(row)
            pd.DataFrame(rows).to_csv(OUT / "汇总.csv", index=False, encoding="utf-8-sig")
            pd.DataFrame(per_well[name]).to_csv(
                OUT / f"{name}_变体对照.csv", index=False, encoding="utf-8-sig")
    return per_well, stop_raw


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "full"
    if mode == "smoke":
        row, stop = run_one("hu2", "基线双关", {})
        print(f"[smoke] ok, stop={stop!r}")
    elif mode == "full":
        raise SystemExit(main_full())
    else:
        raise SystemExit(f"未知模式: {mode}（可用 smoke/full）")
