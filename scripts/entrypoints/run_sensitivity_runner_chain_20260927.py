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
import os
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

# 本次运行输出目录覆盖（--out-dir；None ⇒ 用 OUT_DIR，行为与 2026-09-27 原版一致）。
# 2026-10-06 执行窗口裁定 A 案：0.6 试点必须写**新日期目录**，而本脚本原 OUT_DIR 已冻结
# （含 60 个 *_结果摘要.json，`_assert_not_frozen` 会拒绝重跑）。
_ACTIVE_OUT_DIR: Path | None = None


def _out_dir() -> Path:
    """本次运行输出目录：--out-dir 覆盖优先，否则 = OUT_DIR（默认行为不变）。"""
    return _ACTIVE_OUT_DIR if _ACTIVE_OUT_DIR is not None else OUT_DIR

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
    # 2026-10-06 新增重点井。⚠️ 计划 §3 表格写的 auth 路径
    # "results/呼1-003尾管_1D2D耦合模型/…" 在磁盘上**不存在**（ht 井目录名不带“尾管”）；
    # 依计划「以磁盘实际文件名为准」修正为下列真实路径。
    "呼1-003": {
        "loader": ("cemdisp.data.loaders.ht1_003_loader", "load_ht1_003_tailpipe"),
        "stop": ("cemdisp.runners.ht1_003_tailpipe", "annulus_stop_time_s"),
        "auth": PROJECT_ROOT / "results" / "呼1-003_1D2D耦合模型" / "呼1-003_1D2D耦合模型_结果摘要.json",
        # 该摘要生成于 2026-09-30（早于本轮 HEAD ca8bba7 的 T2 系列提交）⇒
        # 不作「当前 HEAD」判据，差异仅归因报告（计划 §3 0.6「呼1-003 首建归因报告」）。
        "auth_current_head": False,
    },
}


def _resolve(pair: tuple[str, str]):
    return getattr(importlib.import_module(pair[0]), pair[1])


# 冻结目录（红线：零写入）——**显式 --out-dir 亦不得指向**（2026-10-06 加）
FROZEN_DIRS: tuple[Path, ...] = (OUT_DIR,)


def _assert_not_frozen(out_dir: Path) -> None:
    """冻结守卫：目录内已有 ``*_结果摘要.json`` 即视为已跑批冻结、拒绝写入。

    在 main()/phase_zero()/phase_matrix()/phase_pilot() 等一切可被外部调用的
    入口**顶部、任何 mkdir/写文件之前**调用。确需重生成请显式设环境变量
    ``SENS_ALLOW_FROZEN_REGEN=1``（取值恰为 "1" 才放行；未设或其它值一律拒绝）。

    **两层语义（2026-10-06 执行窗口修订）**：
    ① **denylist 无条件生效**——out_dir 命中 :data:`FROZEN_DIRS`（本批冻结目录
       2026-09-27）一律拒绝，显式 ``--out-dir`` 也不能绕过；
    ② **内容守卫只对默认目录生效**——不传 ``--out-dir`` 时行为与 2026-09-27 原版
       完全一致；传了 ``--out-dir`` 即视为「新批次目录」，不设内容守卫。
       理由：新目录一旦跑过第一口井，其自身产物会触发旧口径的守卫，把**逐井
       追加 / 断点续跑**整条路堵死（本轮实测：第二口井起全部被拒）。
    """
    out_dir = Path(out_dir)
    if out_dir.resolve() in {d.resolve() for d in FROZEN_DIRS}:
        raise RuntimeError(
            "目录为冻结目录（红线：零写入）：" + str(out_dir) + chr(10)
            + "新批次请用 --out-dir 指向新日期目录"
        )
    if _ACTIVE_OUT_DIR is not None:
        return                      # 显式 --out-dir：新批次目录，不做内容守卫
    if os.environ.get("SENS_ALLOW_FROZEN_REGEN") == "1":
        return
    if out_dir.exists() and any(out_dir.glob("*_结果摘要.json")):
        raise RuntimeError(
            "目录已冻结（已存在 *_结果摘要.json）：" + str(out_dir) + chr(10)
            + "确需重生成请显式设环境变量 SENS_ALLOW_FROZEN_REGEN=1"
        )


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
    case_json = _out_dir() / f"{well_key}_{tag}_结果摘要.json"
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
    _assert_not_frozen(_out_dir())   # 入口顶部，早于任何 mkdir/写文件
    _out_dir().mkdir(parents=True, exist_ok=True)
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
        with (_out_dir() / "验收_zero对照.csv").open("w", encoding="utf-8-sig", newline="") as fh:
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
        (_out_dir() / "验收_zero对照.md").write_text("\n".join(md), encoding="utf-8")
        print("\n".join(md), flush=True)


def phase_matrix(well_keys: list[str]) -> None:
    """19 变体（依赖 zero 已存在）；结束后从全部落盘 JSON 幂等重建合并汇总表。"""
    _assert_not_frozen(_out_dir())   # 入口顶部，早于任何落盘 JSON 复用/重算与汇总写盘
    variants = build_variants()
    for wk in well_keys:
        cfg = WELLS[wk]
        zero_json = _out_dir() / f"{wk}_zero_结果摘要.json"
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


# ---------------------------------------------------------------- D8-A 试点
# 2026-10-06 执行窗口裁定（计划 §3 0.6 / §1 Q12-D8-A）：只做**镜像脚本试点**；
# 产品化（改 8 井 runner 构造）等试点结果另裁（§6-8）。本段只跑每井
# (zero, Ton_static) 一对：zero = T-off 恒等（镜像基线），Ton_static = 温变开启静温档。
PILOT_TAG_TON = "Ton_static"


def _pilot_pairs() -> list[tuple]:
    """(tag, well_fn, fluid_fn, sched_fn, run_opts) 试点对。"""
    return [
        ("zero", _identity, _identity, _identity, None),
        (PILOT_TAG_TON, _identity, _identity, _identity,
         {"enable_temperature_rheology": True, "temperature_mode": "static",
          "enable_yield_gate": None}),
    ]


def phase_pilot(well_keys: list[str]) -> None:
    """D8-A runner 链 T-on 试点：每井 zero（对照）+ Ton_static（温变）。

    验收（计划 §3 Phase 0.6）：zero 对照呼101 / 呼1-004 权威 ≤0.05pp；
    呼1-003 首建（权威=2026-09-30 旧产物）⇒ 差异仅归因报告，不判失败。
    """
    _assert_not_frozen(_out_dir())
    _out_dir().mkdir(parents=True, exist_ok=True)
    for wk in well_keys:
        cfg = WELLS[wk]
        loader, stop_fn = _resolve(cfg["loader"]), _resolve(cfg["stop"])
        print(f"=== {wk} pilot（zero + {PILOT_TAG_TON}）===", flush=True)
        for tag, well_fn, fluid_fn, sched_fn, run_opts in _pilot_pairs():
            _run_one(wk, loader, stop_fn, tag, well_fn, fluid_fn,
                     sched_fn, None, run_opts)
    # 幂等重建：表格覆盖目录内**已跑过**的全部井（支持逐井追加 / 断点续跑），
    # 而不是只输出本次 --well 的那一口井。
    rows = []
    for wk in [w for w in WELLS
               if (_out_dir() / f"{w}_zero_结果摘要.json").exists()
               and (_out_dir() / f"{w}_{PILOT_TAG_TON}_结果摘要.json").exists()]:
        cfg = WELLS[wk]
        _ze, _zn = _final_of(_out_dir() / f"{wk}_zero_结果摘要.json")
        z = {"η_E": _ze, "η_N": _zn}
        _te, _tn = _final_of(_out_dir() / f"{wk}_{PILOT_TAG_TON}_结果摘要.json")
        t = {"η_E": _te, "η_N": _tn}
        auth_e = auth_n = None
        verdict = "无权威摘要（首建）"
        if cfg["auth"].exists():
            auth_e, auth_n = _final_of(cfg["auth"])
            de, dn = (z["η_E"] - auth_e) * 100.0, (z["η_N"] - auth_n) * 100.0
            if abs(de) <= 0.05 and abs(dn) <= 0.05:
                verdict = "PASS(|Δ|≤0.05pp)"
            elif cfg["auth_current_head"]:
                verdict = "FAIL(权威=当前HEAD仍超差→镜像失真,停跑排查)"
            else:
                verdict = "REPORT(权威非当前HEAD,差异归因)"
        rows.append({
            "井": wk,
            "zeroη_E": z["η_E"], "zeroη_N": z["η_N"],
            "Ton_staticη_E": t["η_E"], "Ton_staticη_N": t["η_N"],
            "ΔTon−zero_ηE_pp": (t["η_E"] - z["η_E"]) * 100.0,
            "ΔTon−zero_ηN_pp": (t["η_N"] - z["η_N"]) * 100.0,
            "权威η_E": auth_e if auth_e is not None else "",
            "权威η_N": auth_n if auth_n is not None else "",
            "zero−权威_ηE_pp": (z["η_E"] - auth_e) * 100.0 if auth_e is not None else "",
            "zero−权威_ηN_pp": (z["η_N"] - auth_n) * 100.0 if auth_n is not None else "",
            "判定": verdict,
        })
    with (_out_dir() / "试点_runner链_Ton.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    md = ["# runner 链 T-on 试点（D8-A，2026-10-06）", "",
          "口径：cemdisp/runners/<井>.py 构造的忠实镜像（1D T1 生产开关、F2 停算、"
          "2D 纯默认 nz=250）；Ton_static = enable_temperature_rheology=True + "
          "temperature_mode=static（无瞬态表井取静温剖面）。",
          "⚠️ 本试点按 Phase 0.0 密度「就近取」口径（Q16）执行。", "",
          "| 井 | zero η_E | zero η_N | Ton_static η_E | Ton_static η_N | "
          "ΔTon−zero η_N/pp | 权威 η_N | zero−权威 η_N/pp | 判定 |",
          "|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        ae = f"{r['权威η_N']:.4f}" if r["权威η_N"] != "" else "—"
        dn = f"{r['zero−权威_ηN_pp']:+.3f}" if r["zero−权威_ηN_pp"] != "" else "—"
        md.append(
            f"| {r['井']} | {r['zeroη_E']:.4f} | {r['zeroη_N']:.4f} | "
            f"{r['Ton_staticη_E']:.4f} | {r['Ton_staticη_N']:.4f} | "
            f"{r['ΔTon−zero_ηN_pp']:+.3f} | {ae} | {dn} | {r['判定']} |"
        )
    (_out_dir() / "试点_runner链_Ton.md").write_text("\n".join(md), encoding="utf-8")
    print("\n".join(md), flush=True)


def _rebuild_combined_table() -> None:
    """从落盘 JSON 幂等重建三井合并汇总（缺档跳过）。"""
    variant_names = ["zero"] + [v[0] for v in build_variants()]
    rows = []
    for wk in WELLS:
        zero_json = _out_dir() / f"{wk}_zero_结果摘要.json"
        if not zero_json.exists():
            continue
        z_e, z_n = _final_of(zero_json)
        for tag in variant_names:
            cj = _out_dir() / f"{wk}_{tag}_结果摘要.json"
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
    with (_out_dir() / "汇总表_敏感性变体_runner链.csv").open("w", encoding="utf-8-sig", newline="") as fh:
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
    (_out_dir() / "汇总表_敏感性变体_runner链.md").write_text("\n".join(md), encoding="utf-8")
    print(f"\n合并汇总重建完成：{len(rows)} 行 → {_out_dir()}", flush=True)


def main() -> None:
    global _ACTIVE_OUT_DIR
    ap = argparse.ArgumentParser(
        description="runner 链敏感性批（zero / matrix / pilot 三段式）"
    )
    ap.add_argument("--phase", choices=["zero", "matrix", "pilot"], required=True)
    ap.add_argument("--well", choices=[*WELLS.keys(), "all"], default="all")
    ap.add_argument("--out-dir", default=None,
                    help="覆盖输出目录（默认 = OUT_DIR，行为与 2026-09-27 原版"
                         "一致）。新批次必须指向新日期目录"
                         "（红线：补跑新目录带日期后缀、禁止覆盖既有产物）")
    args = ap.parse_args()
    if args.out_dir:
        _ACTIVE_OUT_DIR = Path(args.out_dir).resolve()
    _assert_not_frozen(_out_dir())   # main() 顶部，早于任何 mkdir/写文件
    well_keys = list(WELLS) if args.well == "all" else [args.well]
    if args.phase == "zero":
        phase_zero(well_keys)
    elif args.phase == "matrix":
        phase_matrix(well_keys)
    else:
        phase_pilot(well_keys)


if __name__ == "__main__":
    main()
