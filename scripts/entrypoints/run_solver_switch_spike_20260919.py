"""求解器开关 × 输入变体 的通用探针入口（2026-09-19 起）。

用途
----
现有脚本各自绑定一种求解器口径（``run_sensitivity_current_20260916.py`` 硬编码
``CORRECTED_KW``、``run_p1_yield_gate_20260919.py`` 固定叠 ``enable_stream_yield_gate``）。
本脚本把"求解器开关"变成命令行参数，避免为每个新开关再增生一个脚本。

当前批次 = **B-3 spike**：``enable_power_law_gap_correction=True``（幂律间隙律
``I₁·(H/H̄)^{1/n−1}``）下，水泥浆 n 杠杆是否复活。

判据（对照 ``results/敏感性变体_当前口径_2026-09-16/`` 的 B-3 关同变体）：
- 09-16 实测（B-3 关）：n±0.1 仅 ±1.75 pp（η_E）/ ±4~5 pp（η_N）；
- 若本批次跃升到几十 pp 量级 ⇒ 缝隙律确是被塌缩掉的通道；
- 若仍是几 pp ⇒ 通道判断错误，需回到逐格剪切率黏度那条路。

口径纪律
--------
- 基座 = ``rerun_all_wells_corrected.CORRECTED_KW``（与权威八井同口径），开关由
  ``--set`` 叠加；不复制常量。
- 输出目录由 ``--out`` 指定，**默认绝不写权威目录**。
- 断点续跑：已存在的 ``*_结果摘要.json`` 复用。

用法
----
    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -u \
        scripts/entrypoints/run_solver_switch_spike_20260919.py \
        --set enable_power_law_gap_correction=True \
        --wells 呼101,呼103 --variants baseline,cement_n_p0.1,cement_n_m0.1 \
        --out results/spike_B3缝隙律_2026-09-19
"""
from __future__ import annotations

import argparse
import csv
import importlib
import json
import sys
import time
from pathlib import Path

from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver

_SCRIPTS_DIR = Path(__file__).resolve().parents[1]
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))
from entrypoints.rerun_all_wells_corrected import (  # noqa: E402
    CORRECTED_KW,
    NZ,
    _stop_t,
    _total_t,
)
from entrypoints.run_sensitivity_current_20260916 import shift_cement_n  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[2]
G3_DIR = PROJECT_ROOT / "results" / "源模型口径重跑_2026-09-14" / "单井结果"      # B-3 关基线
SENS_DIR = PROJECT_ROOT / "results" / "敏感性变体_当前口径_2026-09-16"            # B-3 关变体

# 井名 → (loader 模块, loader 函数, G3 基线文件名, 09-16 目录中的井前缀)
WELLS: dict[str, tuple[str, str, str, str]] = {
    "呼101": ("cemdisp.data.loaders.hu101_loader", "load_hu101_tailpipe",
              "hu101_corrected_on.json", "呼101"),
    "呼103": ("cemdisp.data.loaders.hu103_loader", "load_hu103_tailpipe",
              "hu103_corrected_on.json", "呼103"),
}

# 变体名 → (fluids 变换, 09-16 目录中的同变体名；None = 基线档，不查对照)
VARIANT_FNS = {
    "baseline": (lambda fs: fs, None),
    "cement_n_p0.1": (lambda fs: shift_cement_n(fs, +0.1), "cement_n_p0.1"),
    "cement_n_m0.1": (lambda fs: shift_cement_n(fs, -0.1), "cement_n_m0.1"),
}


def _coerce(v: str):
    low = v.strip().lower()
    if low in ("true", "false"):
        return low == "true"
    try:
        return int(v)
    except ValueError:
        pass
    try:
        return float(v)
    except ValueError:
        return v


def _run_one(loader, fluid_fn, kw) -> tuple[dict, dict, float]:
    well, fluids, schedule, _ = loader()
    fluids = fluid_fn(fluids)

    casing = CasingFlowSolver(
        enable_gravity=True, mixing_contact_time=True,
        plug_face_zero_mixing=True, has_plug=True,
    )
    cr = casing.run(well, fluids, schedule)
    inlet = build_coupled_annulus_inlet_provider(
        cr, casing, fluids, split_cement_phases=True,
    )
    tt = min(_total_t(schedule) + 1200.0, _stop_t(cr, fluids))
    solver = AnnulusD2DGASolver(total_t=tt, nz=NZ, enable_cfl_adaptive=True, **kw)
    t0 = time.perf_counter()
    res = solver.run(well, fluids, inlet, schedule=schedule)
    return res.summary["最终结果"], res.summary, round(time.perf_counter() - t0, 1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    help="叠加到 CORRECTED_KW 的求解器开关，可重复")
    ap.add_argument("--wells", default="呼101,呼103")
    ap.add_argument("--variants", default="baseline,cement_n_p0.1,cement_n_m0.1")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    overrides = {}
    for kv in args.set:
        k, _, v = kv.partition("=")
        overrides[k.strip()] = _coerce(v)
    kw = dict(CORRECTED_KW, **overrides)
    out_dir = (PROJECT_ROOT / args.out) if not Path(args.out).is_absolute() else Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    wells = [w.strip() for w in args.wells.split(",") if w.strip()]
    variants = [v.strip() for v in args.variants.split(",") if v.strip()]
    print(f"开关叠加 = {overrides}\n最终 kw = {kw}\n输出 = {out_dir}", flush=True)

    rows: list[dict] = []
    for wk in wells:
        mod, fn, base_file, sens_prefix = WELLS[wk]
        loader = getattr(importlib.import_module(mod), fn)
        base = json.loads((G3_DIR / base_file).read_text(encoding="utf-8"))
        be, bn = float(base["eta_E"]), float(base["eta_N"])
        print(f"\n=== {wk}  G3(B-3 关)基线 η_E={be:.4f} η_N={bn:.4f} ===", flush=True)

        # 本开关口径下的基线档（用于算杠杆）
        bl = {"η_E": None, "η_N": None}
        for var in variants:
            fluid_fn, sens_name = VARIANT_FNS[var]
            cj = out_dir / f"{wk}_{var}_结果摘要.json"
            if cj.exists():
                summary = json.loads(cj.read_text(encoding="utf-8"))
                final, elapsed = summary["最终结果"], None
            else:
                final, summary, elapsed = _run_one(loader, fluid_fn, kw)
                cj.write_text(json.dumps(summary, ensure_ascii=False, indent=2),
                              encoding="utf-8")
            eta_e = float(final["全井段最终有效顶替效率"])
            eta_n = float(final["窄四分位效率"])
            if var == "baseline":
                bl = {"η_E": eta_e, "η_N": eta_n}

            # 对照：同变体在 09-16（B-3 关）下的值
            ref_e = ref_n = None
            if sens_name is not None:
                rp = SENS_DIR / f"{sens_prefix}_{sens_name}_结果摘要.json"
                if rp.exists():
                    rs = json.loads(rp.read_text(encoding="utf-8"))
                    rf = rs["最终结果"] if "最终结果" in rs else rs
                    ref_e = float(rf["全井段最终有效顶替效率"])
                    ref_n = float(rf["窄四分位效率"])

            rows.append({
                "井名": wk, "变体": var, "开关口径": json.dumps(overrides, ensure_ascii=False),
                "η_E": eta_e, "η_N": eta_n,
                "基线η_E": be, "基线η_N": bn,
                "本口径Δη_E_pp": (eta_e - be) * 100.0,
                "本口径Δη_N_pp": (eta_n - bn) * 100.0,
                "对照(B3关)η_E": ref_e, "对照(B3关)η_N": ref_n,
                "对照Δη_E_pp": None if ref_e is None else (ref_e - be) * 100.0,
                "对照Δη_N_pp": None if ref_n is None else (ref_n - bn) * 100.0,
                "耗时_s": elapsed if elapsed is not None else "",
            })
            print(f"  [{'复用' if elapsed is None else '计算'}] {wk}×{var}: "
                  f"η_E={eta_e:.4f} η_N={eta_n:.4f}"
                  + (f"  (对照B3关 η_E={ref_e:.4f} η_N={ref_n:.4f})" if ref_e else "")
                  + (f" ({elapsed}s)" if elapsed else ""), flush=True)

        # 杠杆量级对比（本口径 vs B3关）
        if bl["η_E"] is not None:
            print(f"  --- {wk} n 杠杆（本口径）---", flush=True)
            for r in [x for x in rows if x["井名"] == wk and x["变体"] != "baseline"]:
                print(f"    {r['变体']}: Δη_E={r['本口径Δη_E_pp']:+.2f}pp "
                      f"Δη_N={r['本口径Δη_N_pp']:+.2f}pp   |   "
                      f"对照(B3关): Δη_E={r['对照Δη_E_pp']:+.2f}pp "
                      f"Δη_N={r['对照Δη_N_pp']:+.2f}pp", flush=True)

    csv_path = out_dir / "汇总表.csv"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\n完成：{len(rows)} 行 → {csv_path}", flush=True)


if __name__ == "__main__":
    main()
