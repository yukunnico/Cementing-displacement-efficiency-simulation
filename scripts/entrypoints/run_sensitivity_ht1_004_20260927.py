"""呼1-004 敏感性补充批（2026-09-27）——薄驱动，复用 09-16 矩阵的变体定义与流水线。

背景与裁定
----------
用户裁定（2026-09-27）：论文敏感性章节用三口井——呼101/呼103 复用
`results/敏感性变体_当前口径_2026-09-16/`（已实证与当前 HEAD 漂移仅 0.002pp），
呼1-004（ht1_004，资料最全 + 09-26 新口径重跑 + 补居中度良好端）由本批新跑。
主指标 = 全段 η_N，CBL 评价窗 η_N 与窜槽指数作辅。

口径（与 run_sensitivity_current_20260916.py 严格同源）
------------------------------------------------------
- 直接 import 其 build_variants()/run_variant()：19 变体（12 杠杆 + standoff 补档，
  与 ±0.1 合成 9 档）；1D CasingFlowSolver T1 生产开关；2D AnnulusD2DGASolver
  (total_t, nz=NZ, enable_cfl_adaptive=True, **CORRECTED_KW)。

与 09-16 批的两处刻意差异（论文口径脚注需披露）
----------------------------------------------
1. 基线不读 `源模型口径重跑_2026-09-14/单井结果/`：呼1-004 经前缘保序修复
   （39a2303，η_E +0.73pp / η_N +1.74pp），09-14 基线已过期 ⇒ 本批先跑 zero
   变体（恒等变换），以其结果为该井自身基线，批内 Δ 全部相对它计算。
2. 输出到新目录 results/敏感性变体_呼1-004_2026-09-27/，不触碰 09-16 既有产物
   （硬防护：OUT_DIR 若落在旧目录内直接 RuntimeError）。

用法
----
    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 conda run --no-capture-output -n shenjingwangluo \
        python scripts/entrypoints/run_sensitivity_ht1_004_20260927.py

断点续跑：同名变体 JSON 已存在则复用（与 09-16 脚本同语义，新目录内安全）。
"""
from __future__ import annotations

import csv
import importlib
import json
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------- 路径与复用
PROJECT_ROOT = Path(__file__).resolve().parents[2]
_SCRIPTS_DIR = PROJECT_ROOT / "scripts"
for p in (str(PROJECT_ROOT), str(_SCRIPTS_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from entrypoints.run_sensitivity_current_20260916 import (  # noqa: E402
    _identity,
    build_variants,
    run_variant,
)

WELL_KEY = "呼1-004"
LOADER_MOD = "cemdisp.data.loaders.ht1_004_loader"
LOADER_FN = "load_ht1_004_tailpipe"

OUT_DIR = PROJECT_ROOT / "results" / "敏感性变体_呼1-004_2026-09-27"
_OLD_DIR = PROJECT_ROOT / "results" / "敏感性变体_当前口径_2026-09-16"
if _OLD_DIR in OUT_DIR.parents or OUT_DIR == _OLD_DIR:
    raise RuntimeError(f"输出目录防护触发：{OUT_DIR} 不得落在 09-16 既有产物目录内")


def _run_one(loader, tag: str, well_fn, fluid_fn, sched_fn, baseline: dict | None,
             run_opts: dict | None = None) -> dict:
    """跑（或复用）一个变体，返回汇总 CSV 行。baseline=None 表示 zero 自身。

    run_opts 由 09-16 共享装配层定义（2026-10-01 5 元组扩展）；本批变体
    run_opts 恒为 None ⇒ 纯 T-off，行为与扩展前逐位一致。
    """
    case_json = OUT_DIR / f"{WELL_KEY}_{tag}_结果摘要.json"
    if case_json.exists():                       # 断点续跑（新目录内安全）
        summary = json.loads(case_json.read_text(encoding="utf-8"))
        final, elapsed = summary["最终结果"], None
        print(f"  [复用] {WELL_KEY} × {tag}", flush=True)
    else:
        t0 = time.perf_counter()
        summary, _extra = run_variant(
            loader, well_fn, fluid_fn, sched_fn, run_opts, well_key=WELL_KEY,
        )
        final = summary["最终结果"]
        elapsed = round(time.perf_counter() - t0, 1)
        case_json.write_text(
            json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8",
        )
        print(f"  [计算] {WELL_KEY} × {tag}: "
              f"η_E={float(final['全井段最终有效顶替效率']):.4f} "
              f"η_N={float(final['窄四分位效率']):.4f} ({elapsed}s)", flush=True)

    eta_e = float(final["全井段最终有效顶替效率"])
    eta_n = float(final["窄四分位效率"])
    row = {
        "井名": WELL_KEY,
        "变体": tag,
        "η_E": eta_e,
        "η_N": eta_n,
        "窜槽": float(final["最终窜槽指数"]),
        "混浆": float(final["最终混浆指数"]),
        "失稳": float(final["最终失稳指数"]),
        "Δη_E_pp": "" if baseline is None else (eta_e - baseline["eta_E"]) * 100.0,
        "Δη_N_pp": "" if baseline is None else (eta_n - baseline["eta_N"]) * 100.0,
        "基线η_E": "" if baseline is None else baseline["eta_E"],
        "基线η_N": "" if baseline is None else baseline["eta_N"],
        "耗时_s": elapsed if elapsed is not None else "",
    }
    return row


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    loader = getattr(importlib.import_module(LOADER_MOD), LOADER_FN)

    # 前置校验：standoff 剖面必须存在（shift_standoff 的硬前提）
    well0, _, _, _ = loader()
    if not well0.standoff_profile:
        raise RuntimeError(f"{WELL_KEY}: loader 未给出 standoff_profile，无法做 standoff 档")

    # ① zero 变体 = 本井自身基线（前缘保序修复后，09-14 基线已过期）
    print(f"=== {WELL_KEY} zero 变体（自身基线，当前代码+09-16 同源口径）===", flush=True)
    zero_row = _run_one(loader, "zero", _identity, _identity, _identity, None)
    baseline = {"eta_E": zero_row["η_E"], "eta_N": zero_row["η_N"]}
    print(f"基线：η_E={baseline['eta_E']:.4f}  η_N={baseline['eta_N']:.4f}", flush=True)

    # ② 19 变体（与 09-16 矩阵逐个同名，三井横向可比；run_opts 恒 None=T-off）
    rows = [zero_row]
    for var_name, well_fn, fluid_fn, sched_fn, run_opts in build_variants():
        rows.append(_run_one(loader, var_name, well_fn, fluid_fn, sched_fn, baseline,
                             run_opts))

    # ③ 汇总 CSV/MD（列结构与 09-16 汇总表一致，同一套分析代码可解析）
    csv_path = OUT_DIR / "汇总表_敏感性变体_呼1-004.csv"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    md = [
        "# 敏感性变体汇总（呼1-004 补充批，2026-09-27）",
        "",
        "口径：与 `results/敏感性变体_当前口径_2026-09-16/` 严格同源（19 变体同名、",
        "nz/CFL/CORRECTED_KW/T1 套管生产开关一致）；**基线 = 本批 zero 变体**（当前代码，",
        "含前缘保序修复 39a2303），不读 09-14 corrected_on（对呼1-004 已过期）。",
        "",
        "| 井名 | 变体 | η_E | η_N | Δη_E/pp | Δη_N/pp | 耗时/s |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        de = f"{r['Δη_E_pp']:+.2f}" if r["Δη_E_pp"] != "" else "—"
        dn = f"{r['Δη_N_pp']:+.2f}" if r["Δη_N_pp"] != "" else "—"
        md.append(f"| {r['井名']} | {r['变体']} | {r['η_E']:.4f} | {r['η_N']:.4f} "
                  f"| {de} | {dn} | {r['耗时_s']} |")
    (OUT_DIR / "汇总表_敏感性变体_呼1-004.md").write_text("\n".join(md), encoding="utf-8")

    print(f"\n完成：{len(rows)} 个变体（含 zero）→ {OUT_DIR}", flush=True)


if __name__ == "__main__":
    main()
