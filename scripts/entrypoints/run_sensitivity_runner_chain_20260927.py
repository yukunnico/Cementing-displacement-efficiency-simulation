"""Runner 链敏感性三井批（2026-09-27）——§4/§5 单链统一口径。

背景与裁定
----------
§5 敏感性章原拟用 09-16 脚本链矩阵（CORRECTED_KW + tt=min(泵总+1200, stop)），
但脚本链与 §4 权威 runner 链的基线差呈井依赖：呼101 η_N 0.12pp / 呼103 4.88pp /
呼1-004 0.00pp。用户裁定（2026-09-27）：三井全部改 runner 链重跑，论文单链统一，
消除同文双基线。呼1-004 一并重跑（成本 ~13min）以保持矩阵内部单链纯净；其
09-27 脚本链批（results/敏感性变体_呼1-004_2026-09-27/）保留作跨链对照证据。

口径（= cemdisp/runners/<well>.py 的忠实镜像）
----------------------------------------------
- 1D：CasingFlowSolver(enable_gravity=True, mixing_contact_time=True,
  plug_face_zero_mixing=True, has_plug=True)（T1 生产口径，与 runner 同）
- 入口：build_coupled_annulus_inlet_provider(cr, casing, fluids2, split_cement_phases=True)
- total_t：各井 runner 的 annulus_stop_time_s（F2 口径，cement_end_time_s 优先）
- 2D：AnnulusD2DGASolver(total_t=tt, nz=250)——纯默认，无 CORRECTED_KW
- 变体：复用 09-16 脚本 build_variants()（19 个同名变体，三井横向可比）

验收标准（zero 变体 vs §4 权威摘要「最终结果」块）
--------------------------------------------------
- 呼101 / 呼1-004（权威=当前 HEAD 09-26 跑）：|Δ| ≤ 0.05pp，否则判"镜像失真"停跑排查；
- 呼103（权威=09-15 旧代码产物）：差异仅归因报告（前缘保序 39a2303 / 停泵沉降
  f4f5ca0 / 带状 Cholesky 58e9b12 之后的代码差），不判失败；此时以本批 zero 为
  §5 基线，并提示 §4 权威目录需同步重跑（D3）以保持单链。

两段式（先报后动）
------------------
    python scripts/entrypoints/run_sensitivity_runner_chain_20260927.py --phase zero
    python ... --phase matrix --well 呼101    # 每井一个后台进程，最多 3 路并行

输出：results/敏感性变体_runner链_2026-09-27/
    {井}_{变体}_结果摘要.json（res.summary 全量，schema 与 09-16 批一致）
    验收_zero对照.csv/.md · 汇总表_敏感性变体_runner链.csv/.md（matrix 阶段末重建）
"""
from __future__ import annotations

import argparse
import csv
import importlib
import json
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
_SCRIPTS_DIR = PROJECT_ROOT / "scripts"
for p in (str(PROJECT_ROOT), str(_SCRIPTS_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

# 变体定义与 09-16 批同源（同名、同变换闭包），保证三井横向可比
from entrypoints.run_sensitivity_current_20260916 import (  # noqa: E402
    _identity,
    annulus_kwargs_from_opts,
    build_temperature_fields,
    build_variants,
    casing_kwargs_from_opts,
    normalize_run_opts,
)

from cemdisp.models2d import AnnulusD2DGASolver  # noqa: E402
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider  # noqa: E402
from cemdisp.transport1d import CasingFlowSolver  # noqa: E402

OUT_DIR = PROJECT_ROOT / "results" / "敏感性变体_runner链_2026-09-27"

WELLS: dict[str, dict] = {
    "呼101": {
        "loader": ("cemdisp.data.loaders.hu101_loader", "load_hu101_tailpipe"),
        "stop": ("cemdisp.runners.hu101_tailpipe", "annulus_stop_time_s"),
        "auth": PROJECT_ROOT / "results" / "呼101尾管_1D2D耦合模型" / "呼101尾管_1D2D耦合模型_结果摘要.json",
        "auth_current_head": True,
    },
    "呼103": {
        "loader": ("cemdisp.data.loaders.hu103_loader", "load_hu103_tailpipe"),
        "stop": ("cemdisp.runners.hu103_tailpipe", "annulus_stop_time_s"),
        "auth": PROJECT_ROOT / "results" / "呼103尾管_1D2D耦合模型" / "呼103尾管_1D2D耦合模型_结果摘要.json",
        "auth_current_head": False,   # 09-15 旧代码产物，差异仅归因报告
    },
    "呼1-004": {
        "loader": ("cemdisp.data.loaders.ht1_004_loader", "load_ht1_004_tailpipe"),
        "stop": ("cemdisp.runners.ht1_004_tailpipe", "annulus_stop_time_s"),
        "auth": PROJECT_ROOT / "results" / "呼1-004_1D2D耦合模型" / "呼1-004_1D2D耦合模型_结果摘要.json",
        "auth_current_head": True,
    },
}


def _resolve(pair: tuple[str, str]):
    return getattr(importlib.import_module(pair[0]), pair[1])


def run_variant_runner_chain(loader, stop_fn, well_fn, fluid_fn, sched_fn,
                             run_opts: dict | None = None, well_key: str = ""):
    """单变体 runner 链流水线——忠实镜像 cemdisp/runners/<well>.py 的构造。

    run_opts 为 2026-10-01 5 元组扩展位（09-16 共享装配层定义）；本批变体恒
    None ⇒ 纯 T-off，与扩展前逐位一致。给 run_opts 时按同源 helper 消费
    （温度总开关/温度场/屈服门覆盖），基线 kw 仍取 runner 纯默认口径。
    """
    opts = normalize_run_opts(run_opts)
    well, fluids, schedule, _ = loader()
    well2, fluids2, schedule2 = well_fn(well), fluid_fn(fluids), sched_fn(schedule)

    field_1d, field_2d, _note = build_temperature_fields(
        well_key, opts["temperature_mode"]
    )
    casing = CasingFlowSolver(
        enable_gravity=True,
        mixing_contact_time=True,
        plug_face_zero_mixing=True,
        has_plug=True,
        **casing_kwargs_from_opts(opts),
    )
    if opts["enable_temperature_rheology"] and field_1d is not None:
        cr = casing.run(well2, fluids2, schedule2, temperature_field=field_1d)
    else:
        cr = casing.run(well2, fluids2, schedule2)
    inlet = build_coupled_annulus_inlet_provider(
        cr, casing, fluids2, split_cement_phases=True,
    )
    tt = stop_fn(casing_result=cr, fluids=fluids2)
    # 纯默认 + run_opts 覆盖（runner 链无 CORRECTED_KW；enable_yield_gate 默认即 True）
    solver = AnnulusD2DGASolver(
        total_t=tt, nz=250, **annulus_kwargs_from_opts(opts),
    )
    if opts["enable_temperature_rheology"] and field_2d is not None:
        res = solver.run(
            well2, fluids2, inlet, schedule=schedule2, temperature_field=field_2d,
        )
    else:
        res = solver.run(well2, fluids2, inlet, schedule=schedule2)
    return res.summary["最终结果"], res.summary


def _final_of(summary_json: Path) -> tuple[float, float]:
    data = json.loads(summary_json.read_text(encoding="utf-8"))
    fr = data["最终结果"]
    return float(fr["全井段最终有效顶替效率"]), float(fr["窄四分位效率"])


def _run_one(well_key: str, loader, stop_fn, tag: str,
             well_fn, fluid_fn, sched_fn, baseline: dict | None,
             run_opts: dict | None = None) -> dict:
    case_json = OUT_DIR / f"{well_key}_{tag}_结果摘要.json"
    if case_json.exists():                       # 断点续跑（新目录内安全）
        summary = json.loads(case_json.read_text(encoding="utf-8"))
        final, elapsed = summary["最终结果"], None
        print(f"  [复用] {well_key} × {tag}", flush=True)
    else:
        t0 = time.perf_counter()
        final, summary = run_variant_runner_chain(
            loader, stop_fn, well_fn, fluid_fn, sched_fn, run_opts,
            well_key=well_key,
        )
        elapsed = round(time.perf_counter() - t0, 1)
        case_json.write_text(
            json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8",
        )
        print(f"  [计算] {well_key} × {tag}: "
              f"η_E={float(final['全井段最终有效顶替效率']):.4f} "
              f"η_N={float(final['窄四分位效率']):.4f} ({elapsed}s)", flush=True)

    eta_e = float(final["全井段最终有效顶替效率"])
    eta_n = float(final["窄四分位效率"])
    return {
        "井名": well_key,
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


def phase_zero(well_keys: list[str]) -> None:
    """zero 变体 + 权威验收对照（呼101/呼1-004 须 ≤0.05pp；呼103 仅归因报告）。"""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for wk in well_keys:
        cfg = WELLS[wk]
        loader, stop_fn = _resolve(cfg["loader"]), _resolve(cfg["stop"])
        print(f"=== {wk} zero（runner 链镜像）===", flush=True)
        row = _run_one(wk, loader, stop_fn, "zero", _identity, _identity, _identity, None)
        if cfg["auth"].exists():
            a_e, a_n = _final_of(cfg["auth"])
            d_e = (row["η_E"] - a_e) * 100.0
            d_n = (row["η_N"] - a_n) * 100.0
            ok = abs(d_e) <= 0.05 and abs(d_n) <= 0.05
            if ok:
                verdict = "PASS(|Δ|≤0.05pp)"
            elif cfg["auth_current_head"]:
                verdict = "FAIL(权威=当前HEAD仍超差→镜像失真,停跑排查)"
            else:
                verdict = "REPORT(权威=09-15旧代码,差异归因后以本zero为§5基线)"
            rows.append({
                "井": wk, "zeroη_E": row["η_E"], "zeroη_N": row["η_N"],
                "权威η_E": a_e, "权威η_N": a_n,
                "Δη_E_pp": d_e, "Δη_N_pp": d_n,
                "权威是否当前HEAD": "是" if cfg["auth_current_head"] else "否(09-15)",
                "判定": verdict,
            })
        else:
            print(f"  ⚠ 权威摘要缺失：{cfg['auth']}", flush=True)

    if rows:
        with (OUT_DIR / "验收_zero对照.csv").open("w", encoding="utf-8-sig", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        md = ["# zero 变体验收对照（runner 链镜像 vs §4 权威摘要）", "",
              "| 井 | zero η_E | zero η_N | 权威 η_E | 权威 η_N | Δη_E/pp | Δη_N/pp | 判定 |",
              "|---|---|---|---|---|---|---|---|"]
        for r in rows:
            md.append(f"| {r['井']} | {r['zeroη_E']:.4f} | {r['zeroη_N']:.4f} "
                      f"| {r['权威η_E']:.4f} | {r['权威η_N']:.4f} "
                      f"| {r['Δη_E_pp']:+.3f} | {r['Δη_N_pp']:+.3f} | {r['判定']} |")
        (OUT_DIR / "验收_zero对照.md").write_text("\n".join(md), encoding="utf-8")
        print("\n".join(md), flush=True)


def phase_matrix(well_keys: list[str]) -> None:
    """19 变体（依赖 zero 已存在）；结束后从全部落盘 JSON 幂等重建合并汇总表。"""
    variants = build_variants()
    for wk in well_keys:
        cfg = WELLS[wk]
        zero_json = OUT_DIR / f"{wk}_zero_结果摘要.json"
        if not zero_json.exists():
            raise RuntimeError(f"{wk}: 缺 zero 基线，请先跑 --phase zero")
        z_e, z_n = _final_of(zero_json)
        baseline = {"eta_E": z_e, "eta_N": z_n}
        loader, stop_fn = _resolve(cfg["loader"]), _resolve(cfg["stop"])
        print(f"\n=== {wk} 19 变体（基线 η_E={z_e:.4f} η_N={z_n:.4f}）===", flush=True)
        for var_name, well_fn, fluid_fn, sched_fn, run_opts in variants:
            _run_one(wk, loader, stop_fn, var_name, well_fn, fluid_fn, sched_fn,
                     baseline, run_opts)
    _rebuild_combined_table()


def _rebuild_combined_table() -> None:
    """从落盘 JSON 幂等重建三井合并汇总（缺档跳过）。"""
    variant_names = ["zero"] + [v[0] for v in build_variants()]
    rows = []
    for wk in WELLS:
        zero_json = OUT_DIR / f"{wk}_zero_结果摘要.json"
        if not zero_json.exists():
            continue
        z_e, z_n = _final_of(zero_json)
        for tag in variant_names:
            cj = OUT_DIR / f"{wk}_{tag}_结果摘要.json"
            if not cj.exists():
                continue
            fr = json.loads(cj.read_text(encoding="utf-8"))["最终结果"]
            e, n = float(fr["全井段最终有效顶替效率"]), float(fr["窄四分位效率"])
            rows.append({
                "井名": wk, "变体": tag, "η_E": e, "η_N": n,
                "窜槽": float(fr["最终窜槽指数"]), "混浆": float(fr["最终混浆指数"]),
                "失稳": float(fr["最终失稳指数"]),
                "Δη_E_pp": "" if tag == "zero" else (e - z_e) * 100.0,
                "Δη_N_pp": "" if tag == "zero" else (n - z_n) * 100.0,
                "基线η_E": z_e, "基线η_N": z_n, "耗时_s": "",
            })
    if not rows:
        return
    with (OUT_DIR / "汇总表_敏感性变体_runner链.csv").open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    md = ["# 敏感性变体汇总（runner 链三井批，2026-09-27）", "",
          "口径：cemdisp/runners/<井>.py 构造的忠实镜像（1D T1 生产开关、F2 停算、",
          "2D 纯默认 nz=250）；基线=各井本批 zero 变体；与 §4 权威摘要的验收对照见",
          "验收_zero对照.md。变体名与 09-16 脚本链矩阵逐个同名，可跨链对照。", "",
          "| 井名 | 变体 | η_E | η_N | Δη_E/pp | Δη_N/pp |",
          "|---|---|---|---|---|---|"]
    for r in rows:
        de = f"{r['Δη_E_pp']:+.2f}" if r["Δη_E_pp"] != "" else "—"
        dn = f"{r['Δη_N_pp']:+.2f}" if r["Δη_N_pp"] != "" else "—"
        md.append(f"| {r['井名']} | {r['变体']} | {r['η_E']:.4f} | {r['η_N']:.4f} | {de} | {dn} |")
    (OUT_DIR / "汇总表_敏感性变体_runner链.md").write_text("\n".join(md), encoding="utf-8")
    print(f"\n合并汇总重建完成：{len(rows)} 行 → {OUT_DIR}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description="runner 链敏感性三井批（两段式）")
    ap.add_argument("--phase", choices=["zero", "matrix"], required=True)
    ap.add_argument("--well", choices=[*WELLS.keys(), "all"], default="all")
    args = ap.parse_args()
    well_keys = list(WELLS) if args.well == "all" else [args.well]
    if args.phase == "zero":
        phase_zero(well_keys)
    else:
        phase_matrix(well_keys)


if __name__ == "__main__":
    main()
