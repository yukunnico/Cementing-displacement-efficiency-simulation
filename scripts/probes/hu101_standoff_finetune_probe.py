"""呼101井居中度微调探针（spike，2026-09-10）。

目的：回答"微调呼101居中度输入，能不能把顶替效率推到现场 CBL 合格率 62.77%"。

标靶：呼101 CBL 评价井段(单层套管可评价段) 5699.8–7810m 的 η_E
      当前生产口径 = 0.50346（runners/hu101_tailpipe.py 结果摘要 JSON），
      现场 CBL 合格率 = 62.77%，偏差 −12.42pp（8 井唯一负偏差）。

三段输出：
  1. 现有三套剖面的生产口径重跑：
       legacy 假设     0.38–0.48（loader 默认，standoff_measured=False → e 截断 0.55）
       检测图-扶正器间 0.22–0.78（measured_standoff="between_centralizers" → e 截断放开 0.90）
       检测图-扶正器处 0.60–0.88（measured_standoff="at_centralizers"     → e 截断放开 0.90）
     ⚠️ results/呼101_居中度实测对比/ 是旧脚本（无 T1 开关、无 schedule、手改剖面绕过
        standoff_measured 裁定）产物，口径与生产不符，本脚本重跑取代其结论。
  2. 均匀居中度扫描 SO ∈ {0.50,0.60,0.70,0.80,0.90,0.95}（全井常数）→ η_E 响应曲线。
  3. 反解"达到 62.77% 所需的均居中度"，并与现场证据（悬挂器坐挂失败/座底固井 → 居中度应更低）
     做冲突判定。

口径：与 cemdisp/runners/hu101_tailpipe.py 逐位一致
      （CasingFlowSolver 三开关 + annulus_stop_time_s + split_cement_phases
        + AnnulusD2DGASolver(total_t=stop, nz=250) + run(..., schedule=schedule)）。
      1D 套管流与入口桥接均不消费居中度，故 1D 结果全 case 复用。
      均匀扫描 case 设 standoff_measured=False（保守路径）；因 SO≥0.50 → e=1−SO≤0.50<0.55，
      扫描区间内不发生截断，两条裁定路径在此区间数值等价。

本脚本为一次性探针，不改 loader 默认值。
"""

from __future__ import annotations

import csv
import json
import time
from dataclasses import replace
from pathlib import Path
from typing import Any, cast

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from cemdisp.data.loaders import load_hu101_tailpipe  # noqa: E402
from cemdisp.data.well_spec import DepthValuePoint  # noqa: E402
from cemdisp.models2d import AnnulusD2DGASolver  # noqa: E402
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider  # noqa: E402
from cemdisp.runners.hu101_tailpipe import annulus_stop_time_s  # noqa: E402
from cemdisp.transport1d import CasingFlowSolver  # noqa: E402

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = PROJECT_ROOT / "results" / "呼101_居中度微调探针_2026-09-10"

NZ = 250  # hu101 legacy 口径
CBL_WINDOW = "CBL评价井段(单层套管可评价段)"
TARGET = 0.6277        # 现场 CBL 合格率（呼101，5390–7810m 图头读数）
CBL_WINDOW_DEPTH = (5699.8, 7810.0)

TOP_MD, BOTTOM_MD = 5400.0, 7868.0
SO_SWEEP = (0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

# 与 cemdisp/runners/hu101_tailpipe.py 完全一致的生产口径开关
CASING_KW = dict(
    enable_gravity=True,
    mixing_contact_time=True,
    plug_face_zero_mixing=True,
    has_plug=True,
)


def uniform_standoff_profile(so: float) -> tuple[DepthValuePoint, ...]:
    """全井常数居中度剖面。"""
    return (
        DepthValuePoint(TOP_MD, so),
        DepthValuePoint(BOTTOM_MD, so),
    )


def mean_standoff_of(profile: tuple[DepthValuePoint, ...]) -> float:
    """按 200 点等距插值求剖面均值（与旧对比脚本同口径）。"""
    mds = np.array([p.depth_md_m for p in profile], dtype=float)
    vals = np.array([p.value for p in profile], dtype=float)
    return float(np.mean(np.interp(np.linspace(TOP_MD, BOTTOM_MD, 200), mds, vals)))


def run_case(
    *,
    label: str,
    kind: str,
    well_spec: Any,
    fluids: Any,
    schedule: Any,
    casing_result: Any,
    casing_solver: Any,
    stop_time_s: float,
) -> dict[str, Any]:
    """按生产口径跑单个居中度 case（1D 结果复用）。"""
    t0 = time.perf_counter()
    inlet_provider = build_coupled_annulus_inlet_provider(
        casing_result, casing_solver, fluids, split_cement_phases=True,
    )
    solver = AnnulusD2DGASolver(total_t=stop_time_s, nz=NZ)
    result = solver.run(well_spec, fluids, inlet_provider, schedule=schedule)
    elapsed = time.perf_counter() - t0

    final = cast(dict, result.summary["最终结果"])
    windows = cast(dict, result.summary["评价窗效率"])
    cbl = windows.get(CBL_WINDOW)
    if cbl is None:
        raise KeyError(f"结果摘要缺窗口 {CBL_WINDOW!r}，实有 {list(windows)}")

    profile = well_spec.standoff_profile
    metrics: dict[str, Any] = {
        "case": label,
        "kind": kind,
        "mean_standoff": round(mean_standoff_of(profile), 4),
        "eta_E_well": float(final["全井段最终有效顶替效率"]),
        "eta_N_well": float(final["窄四分位效率"]),
        "eta_E_cbl": float(cbl["eta_E"]),
        "eta_N_cbl": float(cbl["eta_N"]),
        "channeling_index": float(final["最终窜槽指数"]),
        "mixing_index": float(final["最终混浆指数"]),
        "instability_index": float(final["最终失稳指数"]),
        "elapsed_s": round(elapsed, 1),
    }
    print(
        f"  [{label:<22s}] 均SO={metrics['mean_standoff']:.4f}  "
        f"ηE全井={metrics['eta_E_well']:.4f}  ηE_CBL窗={metrics['eta_E_cbl']:.4f}  "
        f"ηN_CBL窗={metrics['eta_N_cbl']:.4f}  ({elapsed:.1f}s)"
    )
    return metrics


def interp_crossing(so: list[float], eta: list[float], target: float) -> float | None:
    """在 (SO, η) 曲线上线性反解 η=target 所需的 SO；不可达返回 None。"""
    so_arr, eta_arr = np.asarray(so, dtype=float), np.asarray(eta, dtype=float)
    order = np.argsort(so_arr)
    so_arr, eta_arr = so_arr[order], eta_arr[order]
    if eta_arr.max() < target:
        return None
    for i in range(len(so_arr) - 1):
        lo, hi = eta_arr[i], eta_arr[i + 1]
        if (lo - target) * (hi - target) <= 0 and hi != lo:
            frac = (target - lo) / (hi - lo)
            return float(so_arr[i] + frac * (so_arr[i + 1] - so_arr[i]))
    return None


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("  呼101井 居中度微调探针（标靶 = 现场 CBL 合格率 62.77%）")
    print("=" * 78)

    # ---- 1D 套管流（不消费居中度，全 case 复用）----
    base_well, fluids, schedule, _ = load_hu101_tailpipe()
    casing_solver = CasingFlowSolver(**CASING_KW)  # type: ignore[arg-type]
    casing_result = casing_solver.run(base_well, fluids, schedule)
    stop_time_s = annulus_stop_time_s(casing_result=casing_result, fluids=fluids)
    print(f"  1D 碰压停泵时刻 = {stop_time_s:.1f}s ；2D 求解 total_t 取同值，nz={NZ}")

    cases: list[dict[str, Any]] = []

    # ---- 第一段：现有三套剖面（生产口径重跑）----
    print("\n[1/3] 现有三套剖面（生产口径重跑）")
    for label, mso in (
        ("legacy假设(0.38-0.48)", None),
        ("检测图-扶正器间", "between_centralizers"),
        ("检测图-扶正器处", "at_centralizers"),
    ):
        well, fl, sch, _ = load_hu101_tailpipe(measured_standoff=mso)
        cases.append(run_case(
            label=label, kind="existing", well_spec=well, fluids=fl, schedule=sch,
            casing_result=casing_result, casing_solver=casing_solver, stop_time_s=stop_time_s,
        ))

    # ---- 第二段：均匀居中度扫描 ----
    print(f"\n[2/3] 均匀居中度扫描 SO ∈ {SO_SWEEP}")
    for so in SO_SWEEP:
        well = replace(
            base_well,
            standoff_profile=uniform_standoff_profile(so),
            standoff_measured=False,  # 保守路径；SO≥0.50 时无截断，与放开路径数值等价
        )
        cases.append(run_case(
            label=f"均匀SO={so:.2f}", kind="sweep", well_spec=well, fluids=fluids,
            schedule=schedule, casing_result=casing_result, casing_solver=casing_solver,
            stop_time_s=stop_time_s,
        ))

    # ---- 汇总 ----
    sweep = [c for c in cases if c["kind"] == "sweep"]
    so_list = [c["mean_standoff"] for c in sweep]
    eta_well = [c["eta_E_well"] for c in sweep]
    eta_cbl = [c["eta_E_cbl"] for c in sweep]
    so_star_well = interp_crossing(so_list, eta_well, TARGET)
    so_star_cbl = interp_crossing(so_list, eta_cbl, TARGET)

    base = cases[0]
    print("\n" + "=" * 78)
    print("  汇总")
    print("=" * 78)
    print(f"  {'case':<22s}{'均SO':>8s}{'ηE全井':>10s}{'ηE_CBL窗':>11s}{'ηN_CBL窗':>11s}")
    for c in cases:
        print(f"  {c['case']:<22s}{c['mean_standoff']:>8.4f}{c['eta_E_well']:>10.4f}"
              f"{c['eta_E_cbl']:>11.4f}{c['eta_N_cbl']:>11.4f}")
    print(f"\n  基线（legacy假设）：ηE全井 {base['eta_E_well']:.4f} / ηE_CBL窗 {base['eta_E_cbl']:.4f}")
    print(f"  标靶（现场CBL 62.77%）：差 {TARGET - base['eta_E_cbl']:+.4f}（{100*(TARGET-base['eta_E_cbl']):+.2f} pp）")
    print(f"  反解所需均居中度：CBL窗口径 SO*={so_star_cbl if so_star_cbl else '>0.95（不可达）'}"
          f" ；全井口径 SO*={so_star_well if so_star_well else '>0.95（不可达）'}")

    # ---- 落盘 ----
    fields = ["case", "kind", "mean_standoff", "eta_E_well", "eta_N_well",
              "eta_E_cbl", "eta_N_cbl", "channeling_index", "mixing_index",
              "instability_index", "elapsed_s"]
    csv_path = OUTPUT_DIR / "居中度微调对比表.csv"
    with csv_path.open("w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for c in cases:
            writer.writerow(c)

    payload = {
        "description": "呼101井居中度微调探针（标靶=现场CBL合格率62.77%）",
        "caliber": "cemdisp/runners/hu101_tailpipe.py 生产口径（T1三开关+annulus_stop_time+nz=250+schedule）",
        "cbl_window": {"name": CBL_WINDOW, "depth_m": CBL_WINDOW_DEPTH, "cbl_pass_rate": TARGET},
        "baseline_reproduction_check": {
            "eta_E_well": base["eta_E_well"], "eta_E_cbl": base["eta_E_cbl"],
            "expected_from_runner": {"eta_E_well": 0.49194279093047893, "eta_E_cbl": 0.5034653062957734},
        },
        "stop_time_s": stop_time_s,
        "cases": cases,
        "required_mean_standoff_for_target": {
            "by_cbl_window_eta_E": so_star_cbl, "by_well_eta_E": so_star_well,
        },
    }
    json_path = OUTPUT_DIR / "居中度微调摘要.json"
    json_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    # ---- 响应曲线图 ----
    fig, ax = plt.subplots(figsize=(9, 6), dpi=160)
    ax.plot(so_list, eta_cbl, "o-", color="#c0392b", lw=2, label="均匀扫描 — CBL窗 η$_E$")
    ax.plot(so_list, eta_well, "s--", color="#2980b9", lw=1.6, label="均匀扫描 — 全井 η$_E$")
    markers = {"legacy假设(0.38-0.48)": "^", "检测图-扶正器间": "v", "检测图-扶正器处": "D"}
    for c in cases[:3]:
        ax.scatter([c["mean_standoff"]], [c["eta_E_cbl"]], marker=markers[c["case"]], s=90,
                   facecolors="none", edgecolors="#27ae60", linewidths=1.8, zorder=5)
        ax.annotate(c["case"], (c["mean_standoff"], c["eta_E_cbl"]),
                    textcoords="offset points", xytext=(6, -12), fontsize=9, color="#27ae60")
    ax.axhline(TARGET, color="#7f8c8d", ls=":", lw=1.5)
    ax.annotate(f"现场 CBL 合格率 {TARGET:.2%}", (0.51, TARGET), textcoords="offset points",
                xytext=(0, 6), fontsize=10, color="#7f8c8d")
    if so_star_cbl is not None:
        ax.axvline(so_star_cbl, color="#8e44ad", ls="-.", lw=1.4)
        ax.annotate(f"所需均SO*={so_star_cbl:.3f}", (so_star_cbl, 0.55),
                    textcoords="offset points", xytext=(6, 0), fontsize=10, color="#8e44ad")
    ax.set_xlabel("全井均居中度 SO")
    ax.set_ylabel("顶替效率 η$_E$")
    ax.set_title("呼101井 居中度 → 顶替效率响应（生产口径，nz=250）")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "居中度-效率响应曲线.png")
    plt.close(fig)

    print(f"\n结果已落盘：{OUTPUT_DIR}")
    print("  - 居中度微调对比表.csv / 居中度微调摘要.json / 居中度-效率响应曲线.png")


if __name__ == "__main__":
    main()
