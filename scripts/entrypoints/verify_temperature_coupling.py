# -*- coding: utf-8 -*-
"""T1-5 三级验证脚本 + 口径差对照表（温压耦合 Task 10，用户检查点）。

三关口径（task-10-brief 裁定定稿，覆盖计划 §7 L1 原文）：

1. **关1 = T-off 基线逐位**：T-off 生产口径跑呼1-004，结果摘要与
   ``results/_baseline_T_off/基线摘要_T_off.json``（Task 1 锚）树 diff + 文件字节
   两级对照，判据 diff=0。
2. **关2 = T≡60°C 自洽逐位**：T-on + ``ConstantTemperatureField(60)`` 对照
   「手写同参参照」——参照 = T-off 求解器吃测试侧**预先** ``fluid_at(f, 60)``
   派生的静态 FluidSpec（求解器内部不再进 fluid_at），逐位 diff=0。
   豁免键：``temperature_rheology_audit``（仅 T-on 追加的审计计数，无物理数值）。
3. **关3 = 全开对照 + L3 物理方向**：T-on + 交付瞬态表（管内 T_in / 环空 T_out）
   对照 T-off 出摘要差异表；L3 方向检查（公式行为 / 深段 η / 宽窄边），
   列预期符号与实测符号。

附：**口径差对照表**（计划 §3.7 全网格温点上「公式值 vs 原 loader 静态参数」
逐相全表 + 分派表 9 行逐相标注）→ ``<out-dir>/口径差对照表.md``。

生产 runner 口径声明
--------------------
四次跑批与 ``cemdisp.runners.ht1_004_tailpipe.run_ht1_004_tailpipe_initial``
**逐调用等价**：同 loader（load_ht1_004_tailpipe）、同 casing 生产开关
（enable_gravity + T1 三开关 mixing_contact_time/plug_face_zero_mixing/has_plug）、
同 ``build_coupled_annulus_inlet_provider(..., split_cement_phases=True)``、
同 ``annulus_stop_time_s``、同 ``AnnulusD2DGASolver(total_t, nz=250)``（其余默认）、
同 summary 序列化（dict(result.summary) + 注入流体现场符合性检查 +
json.dumps ensure_ascii=False, indent=2）。
**唯一差别 = 不跑 run_and_export 的 CSV/PNG/NPZ/GIF 导出**（对照对象是摘要；
也避免刷 ``results/呼1-004_1D2D耦合模型/`` 的 tracked 文件破坏工作树纪律）。

参照构造法（关2，报告同步写清）
------------------------------
* 1D 参照流体：``casing_ref = tuple(fluid_at(f, 60.0) for f in raw)``——
  12 相逐个在测试侧派生（casing 从不合成隔离液，按名逐相消费）。
* 2D 参照流体：wash/spacer 三相（先导浆/隔离液1/隔离液2）**先按 solver 同规则**
  体积加权合成（``_wash_spacer_volume_weights`` + ``_composite_spacer_fluid``
  ——静态方法，与 solver 内 run() 合成路径同一实现），再对合成体
  ``fluid_at(comp, 60.0)`` 派生，以**单相**替换进参照元组（solver 见 len==1
  不再合成 ⇒ 与 T-on「先合成后派生」逐位对齐）；其余 9 相逐个 fluid_at 派生。
* provider/取 provenance 用原始 raw 元组（只消费名字/角色/密度，两侧同源）。
* 参照跑 = 全 T-off（enable_temperature_rheology=False，不注入温度场）。

产物（默认 ``results/温度耦合验证/``，可用 ``--out-dir`` 覆盖）：三级验证报告.md、
口径差对照表.md、4 份摘要 JSON。

用法::

    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/entrypoints/verify_temperature_coupling.py
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

_ROOT = Path(__file__).resolve().parents[2]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from cemdisp.data.fluid_provenance import build_injected_fluid_provenance_summary
from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.loaders.ht1_004_loader import load_ht1_004_tailpipe
from cemdisp.data.rheology_vs_temperature import _route, fluid_at, reset_audit
from cemdisp.data.temperature_field import (
    ConstantTemperatureField,
    load_delivered_pair,
)
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.runners.ht1_004_tailpipe import annulus_stop_time_s
from cemdisp.transport1d import CasingFlowSolver

DEFAULT_OUT_DIR = _ROOT / "results" / "温度耦合验证"
# ⚠️ `results/温度耦合验证/` 是 **tracked 冻结产物**（v1，口径差表 P=0.1 MPa 缺省口径）——
# P-1 落地后 v1 作废，重跑**必须**换新目录（红线：禁止覆盖既有产物）：
#     python scripts/entrypoints/verify_temperature_coupling.py --out-dir results/温度耦合验证_P1_<日期>
OUT_DIR = DEFAULT_OUT_DIR
BASELINE_JSON = _ROOT / "results" / "_baseline_T_off" / "基线摘要_T_off.json"
AUDIT_KEY = "temperature_rheology_audit"
WS_ROLES = (FluidRole.WASH, FluidRole.SPACER)
T_REF = 60.0  # 关2 参照派生温度（= ConstantTemperatureField(60)）


def _out_rel() -> str:
    """`OUT_DIR` 相对仓库根的显示路径——正文/日志一律经此派生，切目录后不得残留旧路径。"""
    try:
        return OUT_DIR.relative_to(_ROOT).as_posix()
    except ValueError:  # 目录在仓外
        return str(OUT_DIR)

# 计划 §3.7 回归锚点网格（100+ 取 100.1，与 tests/contract/test_rheology_vs_temperature.py 同口径）
CEM_GRID: list[tuple[str, float]] = [
    ("20", 20.0), ("60", 60.0), ("100−", 100.0), ("100+", 100.1),
    ("120", 120.0), ("155", 155.0), ("170", 170.0),
]
MUD_GRID: list[tuple[str, float]] = [
    ("40", 40.0), ("50", 50.0), ("60", 60.0), ("70", 70.0), ("80", 80.0),
]


# --------------------------------------------------------------------------- #
# 跑批（生产 runner 逐调用等价，见模块 docstring）
# --------------------------------------------------------------------------- #
def _casing_solver(enable_t: bool) -> CasingFlowSolver:
    return CasingFlowSolver(
        enable_gravity=True,
        mixing_contact_time=True,
        plug_face_zero_mixing=True,
        has_plug=True,
        enable_temperature_rheology=enable_t,
    )


def run_variant(
    tag: str,
    *,
    enable_t: bool,
    field_casing: Any = None,
    field_annulus: Any = None,
    casing_fluids: tuple[FluidSpec, ...] | None = None,
    annulus_fluids: tuple[FluidSpec, ...] | None = None,
) -> dict:
    """跑一次呼1-004 全链（1D casing → 鞋口 → 2D 环空），落摘要 JSON。

    字段 *_fluids 为 None ⇒ 用 loader 原始元组（T-on 两侧均如此——派生在
    solver 内部经 enable_temperature_rheology 自行完成）。
    """
    well, fluids_raw, schedule, _ = load_ht1_004_tailpipe()
    casing_fluids = fluids_raw if casing_fluids is None else casing_fluids
    annulus_fluids = fluids_raw if annulus_fluids is None else annulus_fluids

    t0 = time.perf_counter()
    casing = _casing_solver(enable_t)
    casing_result = casing.run(well, casing_fluids, schedule,
                               temperature_field=field_casing)
    provider = build_coupled_annulus_inlet_provider(
        casing_result, casing, fluids_raw, split_cement_phases=True)
    stop_t = annulus_stop_time_s(casing_result=casing_result, fluids=fluids_raw)
    annulus = AnnulusD2DGASolver(
        total_t=stop_t, nz=250, enable_temperature_rheology=enable_t)
    result = annulus.run(
        well, annulus_fluids, provider,
        schedule=schedule, temperature_field=field_annulus,
    )

    payload = dict(result.summary)
    payload["注入流体现场符合性检查"] = build_injected_fluid_provenance_summary(
        well.well_name, schedule, fluids_raw)
    elapsed = time.perf_counter() - t0
    text = json.dumps(payload, ensure_ascii=False, indent=2)
    (OUT_DIR / f"{tag}_摘要.json").write_text(text, encoding="utf-8")
    print(
        f"[{tag}] 完成 {elapsed:.1f}s stop_t={stop_t:.3f}s "
        f"η_E={payload.get('effective_efficiency')!r} "
        f"η_N={payload.get('eta_narrow')!r} 键数={len(payload)}",
        flush=True,
    )
    return {
        "tag": tag,
        "payload": payload,
        "json_text": text,
        "stop_t": stop_t,
        "elapsed": elapsed,
        "result": result,
    }


def build_reference_fluids() -> tuple[tuple[FluidSpec, ...], tuple[FluidSpec, ...],
                                      FluidSpec, FluidSpec]:
    """关2 手写同参参照流体（构造法见模块 docstring「参照构造法」节）。

    返回 (casing_ref, annulus_ref, comp_raw, comp60)。
    """
    _, fluids, schedule, _ = load_ht1_004_tailpipe()
    casing_ref = tuple(fluid_at(f, T_REF) for f in fluids)

    ws = [f for f in fluids if f.role in WS_ROLES]
    weights = AnnulusD2DGASolver._wash_spacer_volume_weights(ws, schedule)
    comp_raw = AnnulusD2DGASolver._composite_spacer_fluid(ws, weights)
    comp60 = fluid_at(comp_raw, T_REF)

    ann_ref: list[FluidSpec] = []
    inserted = False
    for f in fluids:
        if f.role in WS_ROLES:
            if not inserted:
                ann_ref.append(comp60)
                inserted = True
            continue
        ann_ref.append(fluid_at(f, T_REF))
    return casing_ref, tuple(ann_ref), comp_raw, comp60


# --------------------------------------------------------------------------- #
# 逐位对照
# --------------------------------------------------------------------------- #
def deep_diff(a: Any, b: Any, path: str = "$", out: list | None = None) -> list:
    """JSON 树逐叶子对照：float 用 ==（IEEE 逐位），NaN≡NaN 记等。"""
    if out is None:
        out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            if k not in a:
                out.append((f"{path}.{k}", "<缺键>", b[k]))
            elif k not in b:
                out.append((f"{path}.{k}", a[k], "<缺键>"))
            else:
                deep_diff(a[k], b[k], f"{path}.{k}", out)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append((f"{path}[长度]", len(a), len(b)))
        for i, (x, y) in enumerate(zip(a, b)):
            deep_diff(x, y, f"{path}[{i}]", out)
    else:
        same = a == b or (
            isinstance(a, float) and isinstance(b, float)
            and math.isnan(a) and math.isnan(b)
        )
        if not same:
            out.append((path, a, b))
    return out


def fmt_diffs(diffs: list, limit: int = 20) -> str:
    if not diffs:
        return "无（diff=0）"
    lines = [f"共 {len(diffs)} 处不等，前 {min(limit, len(diffs))} 条："]
    for p, x, y in diffs[:limit]:
        lines.append(f"- `{p}`：T侧={x!r} vs 对照={y!r}")
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# 口径差对照表（§3.7 全网格：公式值 vs loader 静态参数）
# --------------------------------------------------------------------------- #
# 分派表 9 行（task-3-brief 裁定定稿）与呼1-004 逐相命中
DISPATCH_ROWS = [
    ("行1", "钻井液 → 钻井液式（八井通用）", ["钻井液"]),
    ("行2", "先导浆、平衡液（仅非 SPACER role）→ 钻井液式（model_assumption）；"
     "role=SPACER 一律优先走行4/5（C1，含 2D 合成相）", ["先导浆"]),
    ("行3", "压塞液/替钻井液/井浆/基液/保护液 → 不替换（替浆链常数现状）",
     ["压塞液", "替钻井液", "基液", "井浆", "保护液"]),
    ("行4", "隔离液 ρ∈[1.95,2.05] → 六系数密度插值（Q2b）", ["隔离液1"]),
    ("行5", "隔离液 ρ 域外 → clamp 就近端（model_assumption）", ["隔离液2"]),
    ("行6", "领/尾/中间浆 ρ∈[1.88,1.92] → cement_B", ["尾浆"]),
    ("行7", "领/尾/中间浆 ρ∈[2.08,2.12] → cement_A", []),
    ("行8", "领/尾/中间浆 ρ∈(1.92,2.08) → cement_interp", ["领浆"]),
    ("行9", "冲洗液（FLUSHER） → 不替换", ["冲洗液（FLUSHER）"]),
]
# 代表温点：泥浆族=井内 clamp 端 80°C；水泥/隔离液=井内中段 120°C（P=0.1 常压缺省）
REP_T = {"mud": 80.0, "spacer": 120.0, "cement": 120.0}


def _pct(formula: float, loader: float) -> str:
    if loader == 0:
        return "—"
    return f"{(formula - loader) / loader * 100.0:+.1f}%"


def build_caliber_table(
    p_mpa: float | None = None,
    p_source: str = "",
) -> tuple[str, list[dict]]:
    """生成口径差对照表 markdown + 代表温点关键数字列表。

    ``p_mpa``（P-1，2026-10-06）：隔离液族的**真实压力** [MPa]。``None`` ⇒ 常压 0.1 MPa
    缺省（v1 口径，记 ``p_default`` 审计）；给值 ⇒ 隔离液族按该压力求值（v2 口径），
    并在表中增列「@P」对照。**效应面只有隔离液族**——水泥/钻井液公式纯温度，不吃 P。
    """
    _, fluids, _, _ = load_ht1_004_tailpipe()
    row_of: dict[str, tuple[str, str]] = {}
    for row_id, desc, names in DISPATCH_ROWS:
        for n in names:
            row_of[n] = (row_id, desc)

    rep_rows: list[dict] = []
    md: list[str] = []
    md.append("# 口径差对照表：公式值 vs 原 loader 静态参数（呼1-004）\n")
    md.append(
        "- 公式侧 = `cemdisp.data.rheology_vs_temperature.fluid_at`（T0-2 逐字录入公式，"
        "Task 3 实现）；loader 侧 = `ht1_004_loader._build_fluids()` 静态常量——两侧来源独立。\n"
        "- 网格 = 计划 §3.7 全网格温点：水泥/隔离液族 20/60/100−(=100.0)/100+(=100.1)/120/155/170 °C；"
        "泥浆族 40/50/60/70/80 °C。\n"
        + (
            ("- 隔离液压力 P：**P-1 已接线** —— 本表按 `HydrostaticPressureField` 的鞋深单点"
             f" **P={p_mpa:.4f} MPa**（{p_source}）求值；表中「@P」列 = 该压力下的公式值，"
             "未标 P 的列 = v1 常压缺省 0.1 MPa 口径（**v1 已作废，仅留对照**）。\n")
            if p_mpa is not None else
            "- 隔离液压力 P：与生产 T-on 路径同口径 **P 未传入 → 常压 0.1 MPa 缺省**"
            "（压力耦合属 Phase P，未接线）。\n"
        ) +
        "- μp 单位 Pa·s、τy 单位 Pa；Δ = 公式 − loader；Δ% = (公式−loader)/loader。\n"
        "- 「不替换」行（分派表行3/行9）：T-on 对该相原样返回 loader 常数 ⇒ Δ≡0。\n"
    )
    md.append("\n## 分派表 9 行 × 呼1-004 命中相\n")
    md.append("| 分派行 | 规则（task-3-brief 裁定） | 本井命中相 |")
    md.append("|---|---|---|")
    for row_id, desc, names in DISPATCH_ROWS:
        md.append(f"| {row_id} | {desc} | {'、'.join(names) if names else '—'} |")

    for f in fluids:
        family, note = _route(f)
        row_id, row_desc = row_of.get(f.name, ("—", "未映射"))
        loader_tauy = f.yield_stress_pa
        loader_mup = f.plastic_viscosity_pa_s
        d_gcc = f.density_kg_m3 / 1000.0
        md.append(f"\n## {f.name}（role={f.role.value}，ρ={d_gcc:.3f} g/cm³）\n")
        md.append(f"- 路由：family=`{family}`"
                  + (f"，note=`{note}`" if note else "")
                  + f"；{row_id}：{row_desc}")
        if family == "no_replace":
            md.append(f"- **不替换**：T-on 沿用 loader 常数（τy={loader_tauy} Pa，"
                      f"μp={loader_mup} Pa·s），Δ≡0。")
            rep_rows.append({
                "phase": f.name, "family": family, "row": row_id,
                "T": None, "tauy_f": loader_tauy, "tauy_l": loader_tauy,
                "mup_f": loader_mup, "mup_l": loader_mup,
            })
            continue
        grid = MUD_GRID if family == "mud" else CEM_GRID
        show_p = (family == "spacer" and p_mpa is not None)
        cols = ("| 温点 | 公式τy | loader YP | Δτy | Δτy% | 公式μp | loader PV | Δμp | Δμp% |"
                + (" 公式τy@P | Δτy%@P | 公式μp@P | Δμp%@P |" if show_p else ""))
        md.append("\n" + cols)
        md.append("|---" * (9 + (4 if show_p else 0)) + "|")
        reset_audit()
        rep_t = REP_T[family]
        for label, t in grid:
            out = fluid_at(f, t)
            ft = float(out.yield_stress_pa)
            fm = float(out.plastic_viscosity_pa_s)
            ftp = fmp = None
            extra = ""
            if show_p:
                outp = fluid_at(f, t, p_mpa)
                ftp = float(outp.yield_stress_pa)
                fmp = float(outp.plastic_viscosity_pa_s)
                extra = (f" {ftp:.4f} | {_pct(ftp, loader_tauy)} |"
                         f" {fmp:.6f} | {_pct(fmp, loader_mup)} |")
            md.append(
                f"| {label} °C | {ft:.4f} | {loader_tauy} | "
                f"{ft - loader_tauy:+.4f} | {_pct(ft, loader_tauy)} | "
                f"{fm:.6f} | {loader_mup} | {fm - loader_mup:+.6f} | "
                f"{_pct(fm, loader_mup)} |" + extra
            )
            if abs(t - rep_t) < 1e-9:
                rep_rows.append({
                    "phase": f.name, "family": family, "row": row_id,
                    "T": t, "tauy_f": ft, "tauy_l": loader_tauy,
                    "mup_f": fm, "mup_l": loader_mup,
                    "tauy_f_p": ftp, "mup_f_p": fmp,
                    "p_mpa": p_mpa if show_p else None,
                })
        md.append("")

    md.append("\n## 代表温点汇总（泥浆族=80 °C 井内 clamp 端；水泥/隔离液=120 °C，P=0.1 MPa）\n")
    md.append("| 相 | 分派行 | 代表T | τy 公式 vs loader | τy Δ% | μp 公式 vs loader | μp Δ% |")
    md.append("|---|---|---|---|---|---|---|")
    for r in rep_rows:
        t_label = f"{r['T']:.0f} °C" if r["T"] is not None else "不替换"
        md.append(
            f"| {r['phase']} | {r['row']} | {t_label} | "
            f"{r['tauy_f']:.4f} vs {r['tauy_l']} | {_pct(r['tauy_f'], r['tauy_l'])} | "
            f"{r['mup_f']:.6f} vs {r['mup_l']} | {_pct(r['mup_f'], r['mup_l'])} |"
        )
    return "\n".join(md) + "\n", rep_rows


# --------------------------------------------------------------------------- #
# 关3：L3 物理方向检查
# --------------------------------------------------------------------------- #
def _band_mean(df, top: float, bot: float, col: str) -> float:
    m = (df["井深_m"].to_numpy(float) >= top) & (df["井深_m"].to_numpy(float) <= bot)
    w = df.loc[m, "环空间隙_m"].to_numpy(float)
    v = df.loc[m, col].to_numpy(float)
    return float(np.average(v, weights=w))


def _sign(x: float) -> str:
    if x > 0:
        return "+"
    if x < 0:
        return "−"
    return "0"


def l3_checks(off: dict, tab: dict) -> tuple[list[dict], dict]:
    """L3 物理方向检查：返回 (检查行列表, 附带数值 dict)。

    预期符号的机制依据（写进报告）：
    * L3-a 水泥 τy(T)：二次式开口向上（组A顶点≈70.1 °C、组B顶点≈88.5 °C）+
      100 °C 后升斜段 ⇒ 全网格先降后升；井内温区 [98.7,150.8] °C 已过顶点 ⇒ 深(热)端升。
    * L3-b 水泥 μp(T)：判定式只认**本井相的井内方向**——尾浆=纯组B（网格应单调升、
      井内升）、领浆=组B主导插值（井内升）。组A 指数降是**公式族脚注**（顶点≈24.2 °C
      以上组B 升），不进本井判定式；按 task-3 裁定**公式逐字为准**，预期格只写本井
      实际判定式（组B/插值=井内升），组A 的「降」只留机制脚注。
    * L3-c 深段 η 相对下拉（计划 §7 L3 方向检查）：**实现无深度梯度**——物性吃
      `_representative_temperature` 的**域均标量代表温**（annulus_d2dga.py:1144-1166），
      温度不分深。实际机制 = 温度效应经**时间次序**传递：早期均温低 ⇒ 水泥 τy(T) 低 ⇒
      屈服门弱 ⇒ 先到位相受益（方向与公式行为相容）。故 L3-c 只作**方向记录**
      （(Δ深 − Δ浅) 期望 < 0），**不得**归因为"浅端(低温)门减弱更多"的深度梯度叙事；
      实测方向不符即如实记「否」，不改判据。
    * L3-d 宽窄边：冻结集中在窄边（τw 低），水泥 τy↓ ⇒ 门减弱使窄边受益最大
      ⇒ 预期 Δη_N − Δη_E > 0（机制=屈服门，计划 L3 主张）。
    """
    well, fluids, _, _ = load_ht1_004_tailpipe()
    lead = next(f for f in fluids if f.name == "领浆")
    tail = next(f for f in fluids if f.name == "尾浆")
    _, t_out = load_delivered_pair()

    # 井内温区（环空表 T_out，域 5243.207-7660 m × 全时间轴）
    z = t_out.depth_m
    dom = (z >= well.top_md_m) & (z <= well.bottom_md_m)
    t_lo = float(t_out.table[dom].min())
    t_hi = float(t_out.table[dom].max())

    # 公式层：全网格序列 + 井内端点
    grid_t = [t for _, t in CEM_GRID]
    seq = {}
    for name, f in (("领浆", lead), ("尾浆", tail)):
        seq[name] = {
            "tauy": [float(fluid_at(f, t).yield_stress_pa) for t in grid_t],
            "mup": [float(fluid_at(f, t).plastic_viscosity_pa_s) for t in grid_t],
            "tauy_lo": float(fluid_at(f, t_lo).yield_stress_pa),
            "tauy_hi": float(fluid_at(f, t_hi).yield_stress_pa),
            "mup_lo": float(fluid_at(f, t_lo).plastic_viscosity_pa_s),
            "mup_hi": float(fluid_at(f, t_hi).plastic_viscosity_pa_s),
        }

    checks: list[dict] = []

    # L3-a：τy 先降后升（全网格）+ 井内升
    a_ok = True
    a_detail = []
    for name in ("领浆", "尾浆"):
        s = seq[name]["tauy"]
        dips = min(s) < s[0] and min(s) < s[-1]
        rises = seq[name]["tauy_hi"] > seq[name]["tauy_lo"]
        a_ok = a_ok and dips and rises
        a_detail.append(
            f"{name}: 网格[{', '.join(f'{v:.2f}' for v in s)}] 先降后升={'是' if dips else '否'}；"
            f"井内 {t_lo:.1f}→{t_hi:.1f} °C: "
            f"{seq[name]['tauy_lo']:.3f}→{seq[name]['tauy_hi']:.3f} Pa "
            f"({_sign(seq[name]['tauy_hi'] - seq[name]['tauy_lo'])})"
        )
    checks.append({
        "id": "L3-a", "对象": "水泥 τy(T) 公式行为",
        "预期符号": "全网格先降后升；井内(过顶点) 深端 +（升）",
        "实测": "；".join(a_detail),
        "一致": "是" if a_ok else "否",
    })

    # L3-b：μp 方向（计划散文 vs 分组公式；判定=实测序列与分组公式预期相符）
    b_detail = []
    b_ok = True
    for name in ("领浆", "尾浆"):
        s = seq[name]["mup"]
        grid_dir = "升" if all(s[i + 1] >= s[i] for i in range(len(s) - 1)) else "非单调"
        inw = seq[name]["mup_hi"] - seq[name]["mup_lo"]
        b_detail.append(
            f"{name}: 20→170 °C 网格 {grid_dir} "
            f"[{', '.join(f'{v * 1000:.1f}' for v in s)} mPa·s]；"
            f"井内 {seq[name]['mup_lo'] * 1000:.1f}→{seq[name]['mup_hi'] * 1000:.1f} mPa·s "
            f"({_sign(inw)})"
        )
        # 分组公式预期：尾浆=纯组B（网格应单调升、井内升）；
        # 领浆=插值(w=0.15, B 主导)：井内应升
        if name == "尾浆":
            b_ok = b_ok and grid_dir == "升" and inw > 0
        else:
            b_ok = b_ok and inw > 0
    checks.append({
        "id": "L3-b", "对象": "水泥 μp(T) 公式行为",
        "预期符号": "本井判定式（组B/插值）：井内 μp 升（+）；尾浆网格单调升",
        "实测": "；".join(b_detail),
        "一致": "是" if b_ok else "否",
    })

    # L3-c：深段 η 下拉
    df_off = off["result"].depth_profiles
    df_on = tab["result"].depth_profiles
    top, bot = float(off["payload"]["井段_m"][0]), float(off["payload"]["井段_m"][1])
    mid = 0.5 * (top + bot)
    col = "平均有效顶替效率"
    off_sh, on_sh = _band_mean(df_off, top, mid, col), _band_mean(df_on, top, mid, col)
    off_dp, on_dp = _band_mean(df_off, mid, bot, col), _band_mean(df_on, mid, bot, col)
    d_sh, d_dp = on_sh - off_sh, on_dp - off_dp
    rel = d_dp - d_sh
    checks.append({
        "id": "L3-c", "对象": "深段 η 相对下拉（域均标量代表温，方向记录）",
        "预期符号": "(Δ深 − Δ浅) < 0（深段相对下拉）",
        "实测": (
            f"浅段[{top:.0f},{mid:.0f}] {off_sh:.6f}→{on_sh:.6f} (Δ={d_sh:+.6f})；"
            f"深段[{mid:.0f},{bot:.0f}] {off_dp:.6f}→{on_dp:.6f} (Δ={d_dp:+.6f})；"
            f"Δ深−Δ浅={rel:+.6f} ({_sign(rel)})"
        ),
        "一致": "是" if rel < 0 else "否",
    })

    # L3-d：宽窄边
    e_off = float(off["payload"]["effective_efficiency"])
    e_on = float(tab["payload"]["effective_efficiency"])
    n_off = float(off["payload"]["eta_narrow"])
    n_on = float(tab["payload"]["eta_narrow"])
    d_e, d_n = e_on - e_off, n_on - n_off
    gap = d_n - d_e
    wide_off = _band_mean(df_off, top, bot, "宽边有效效率")
    wide_on = _band_mean(df_on, top, bot, "宽边有效效率")
    edge_off = _band_mean(df_off, top, bot, "窄边有效效率")
    edge_on = _band_mean(df_on, top, bot, "窄边有效效率")
    checks.append({
        "id": "L3-d", "对象": "宽窄边（门减弱窄边受益最大）",
        "预期符号": "Δη_N − Δη_E > 0",
        "实测": (
            f"η_E {e_off:.6f}→{e_on:.6f} (Δ={d_e:+.6f})；"
            f"η_N {n_off:.6f}→{n_on:.6f} (Δ={d_n:+.6f})；"
            f"ΔN−ΔE={gap:+.6f} ({_sign(gap)})；"
            f"剖面宽边均值 {wide_off:.6f}→{wide_on:.6f} (Δ={wide_on - wide_off:+.6f})、"
            f"窄边 {edge_off:.6f}→{edge_on:.6f} (Δ={edge_on - edge_off:+.6f})"
        ),
        "一致": "是" if gap > 0 else "否",
    })

    aux = {
        "t_lo": t_lo, "t_hi": t_hi, "seq": seq,
        "d_sh": d_sh, "d_dp": d_dp, "d_e": d_e, "d_n": d_n,
    }
    return checks, aux


def summary_diff_rows(off: dict, tab: dict) -> list[tuple[str, Any, Any, str]]:
    """关3 摘要差异表（关键指标，T-off vs T-table）。"""
    rows: list[tuple[str, Any, Any, str]] = []

    def add(key: str, a: Any, b: Any) -> None:
        if isinstance(a, (int, float)) and isinstance(b, (int, float)):
            d = f"{b - a:+.6g}" if not (isinstance(a, float) and math.isnan(a)) else "—"
        else:
            d = "≠" if a != b else "="
        rows.append((key, a, b, d))

    for k in ("effective_efficiency", "eta_narrow", "channeling_index",
              "mixing_index", "buoyancy_number"):
        add(k, off["payload"].get(k), tab["payload"].get(k))
    fr_off = off["payload"]["最终结果"]
    fr_on = tab["payload"]["最终结果"]
    for k in fr_off:
        add(f"最终结果/{k}", fr_off[k], fr_on.get(k))
    w_off = off["payload"]["评价窗效率"]
    w_on = tab["payload"]["评价窗效率"]
    for name in w_off:
        for sub in ("eta_E", "eta_N"):
            add(f"评价窗/{name}/{sub}", w_off[name][sub],
                w_on.get(name, {}).get(sub))
    add("stop_t_s", off["stop_t"], tab["stop_t"])
    return rows


# --------------------------------------------------------------------------- #
# 报告
# --------------------------------------------------------------------------- #
def _git_head() -> str:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=_ROOT, capture_output=True, text=True, timeout=20,
        )
        return out.stdout.strip() or "?"
    except Exception:
        return "?"


def write_report(ctx: dict) -> str:
    g1, g2, g3 = ctx["gate1"], ctx["gate2"], ctx["gate3"]
    lines: list[str] = []
    add = lines.append
    add("# 温压耦合 T1-5 三级验证报告（Task 10 用户检查点）\n")
    add(f"- 日期：{ctx['timestamp']}")
    add(f"- 工作目录：`{_ROOT}`")
    add(f"- 分支 HEAD：`{ctx['head']}`（feat/temperature-coupling）")
    add("- 环境：`D:\\apps\\Anaconda\\envs\\cementT\\python.exe`；"
        "`PYTHONIOENCODING=utf-8 PYTHONUTF8=1`")
    add("- 脚本：`scripts/entrypoints/verify_temperature_coupling.py`"
        "（本报告由该脚本一次性生成）")
    add(f"- 产物：`{_out_rel()}/`（本报告、口径差对照表.md、4 份摘要 JSON）\n")

    add("## 三关结论\n")
    add("| 关 | 判据 | 结果 | 证据 |")
    add("|---|---|---|---|")
    add(f"| 关1 T-off=基线逐位 | 摘要树 diff=0 + 文件字节一致 | "
        f"**{g1['verdict']}** | diff 叶子 {g1['n_diff']}；"
        f"字节级 {'一致' if g1['byte_same'] else '不一致'} |")
    add(f"| 关2 T≡60 °C 自洽逐位 | T-on vs 手写参照 树 diff=0（豁免 `{AUDIT_KEY}`） | "
        f"**{g2['verdict']}** | diff 叶子 {g2['n_diff']} |")
    add(f"| 关3 全开对照+L3 方向 | L3 方向检查全部一致 | "
        f"**{g3['verdict']}** | L3 {g3['n_pass']}/{g3['n_total']} 项一致"
        + (f"（不一致：{g3['fail_ids']}）" if g3["fail_ids"] else "") + " |")
    add("")
    add(f"**总结论：{ctx['overall']}**\n")

    add("## 跑批清单与口径声明\n")
    add("| 变体 | 开关 | 温度场 | stop_t (s) | 耗时 (s) | η_E | η_N |")
    add("|---|---|---|---|---|---|---|")
    for r in ctx["runs"]:
        add(f"| {r['tag']} | {'T-on' if r['enable_t'] else 'T-off'} | {r['field']} | "
            f"{r['stop_t']:.3f} | {r['elapsed']:.1f} | "
            f"{r['payload'].get('effective_efficiency')!r} | "
            f"{r['payload'].get('eta_narrow')!r} |")
    add("")
    add("- **生产 runner 口径**：与 `run_ht1_004_tailpipe_initial` 逐调用等价"
        "（loader / casing 三开关+enable_gravity / provider split_cement_phases / "
        "annulus_stop_time_s / AnnulusD2DGASolver(total_t, nz=250) / 摘要序列化同参），"
        "仅不跑 CSV/PNG/GIF 导出（对照对象=摘要；避免刷 `results/呼1-004_1D2D耦合模型/` "
        "tracked 文件）。关1 diff=0 即证该等价复现与 Task 1 生产跑批同码同值。\n")

    # ---- 关1 ----
    add("## 关1：T-off 生产口径 = Task 1 基线（逐位）\n")
    add(f"- 对照对象：`{BASELINE_JSON.relative_to(_ROOT)}`（Task 1 归档锚）")
    add(f"- 方法：双方 `json.loads` 后整树逐叶子 `==`（float=IEEE 逐位，NaN≡NaN）；"
        f"另做**文件字节级**对照（同一 dumps 参数 ensure_ascii=False, indent=2）")
    add(f"- 结果：**{g1['verdict']}** —— diff 叶子数 **{g1['n_diff']}**；"
        f"文件字节{'一致' if g1['byte_same'] else '**不一致**'}")
    if g1["diffs"]:
        add("\n" + fmt_diffs(g1["diffs"]))
    add("\n关键字段对照（baseline / T-off 实跑）：\n")
    add("| 字段 | 基线 | T-off 实跑 |")
    add("|---|---|---|")
    for k, v in g1["key_rows"]:
        add(f"| {k} | {v[0]!r} | {v[1]!r} |")
    add("")

    # ---- 关2 ----
    add("## 关2：T≡60 °C 常数自洽（逐位）\n")
    add("### 参照构造法（手写同参，不经 solver 内 fluid_at）\n")
    add("1. **1D 参照流体**：`tuple(fluid_at(f, 60) for f in raw)`——呼1-004 的 12 相"
        "在测试侧逐相派生为静态 FluidSpec（casing 按名逐相消费，从不合成隔离液）；")
    add("2. **2D 参照流体**：wash/spacer 三相（先导浆/隔离液1/隔离液2）先按 **solver 同规则**"
        "体积加权合成——`AnnulusD2DGASolver._wash_spacer_volume_weights(ws, schedule)` + "
        "`_composite_spacer_fluid(ws, weights)`（run() 内同一对静态方法），"
        "再 `fluid_at(comp, 60)` 派生，以**单相**进参照元组 ⇒ solver 见 len==1 "
        "不再合成，与 T-on「先合成、后派生」顺序逐位对齐；其余 9 相逐相 fluid_at 派生；")
    add("3. provider / provenance 用**原始 raw 元组**（只消费名字/角色/密度，两侧同源）；")
    add("4. 参照跑 = 全 T-off（`enable_temperature_rheology=False`、不注入温度场）；"
        "T-on 跑 = 双求解器开关开 + `ConstantTemperatureField(60)` 注入、"
        "**raw 元组**（派生在 solver 内完成）。\n")
    add(f"结果：**{g2['verdict']}** —— diff 叶子数 **{g2['n_diff']}**"
        f"（豁免键 `{AUDIT_KEY}`：T-on={g2['audit_on']!r}，参照侧无此键——"
        "仅 T-on 追加的审计计数，无物理数值）")
    if g2["diffs"]:
        add("\n" + fmt_diffs(g2["diffs"]))
    add("\n关键字段对照（T-on Const60 / 参照）：\n")
    add("| 字段 | T-on | 参照 |")
    add("|---|---|---|")
    for k, v in g2["key_rows"]:
        add(f"| {k} | {v[0]!r} | {v[1]!r} |")
    add("")

    # ---- 关3 ----
    add("## 关3：全开对照（T_in/T_out 瞬态表）+ L3 物理方向\n")
    add("### 温度表统计（交付表，环空域 5243.207–7660 m × 全时间轴）\n")
    add(f"- T_out 域内：[{ctx['t_lo']:.2f}, {ctx['t_hi']:.2f}] °C"
        "（管内 T_in 域内 [79.67, 150.76] °C；时间轴 0–11940 s 覆盖 stop_t）")
    add(f"- T-on 审计（表格场 run）：`{AUDIT_KEY}` = {ctx['tab_audit']!r}；"
        f"温度场越界次数：T_in oob={ctx['oob_in']}，T_out oob={ctx['oob_out']}"
        "（0=查询全在表域内）\n")
    add("### 摘要差异表（T-off → T-on 表格场）\n")
    add("| 指标 | T-off | T-on 表格 | Δ |")
    add("|---|---|---|---|")
    for k, a, b, d in ctx["diff_rows"]:
        add(f"| {k} | {a!r} | {b!r} | {d} |")
    add("")
    add("### L3 物理方向检查（预期符号 vs 实测符号）\n")
    add("| 检查 | 对象 | 预期符号 | 实测 | 一致 |")
    add("|---|---|---|---|---|")
    for c in ctx["l3"]:
        add(f"| {c['id']} | {c['对象']} | {c['预期符号']} | {c['实测']} | {c['一致']} |")
    add("")
    add("机制说明（预期符号依据）：\n"
        "\n"
        "- **L3-c 不作深度梯度归因**——2D 物性吃 `_representative_temperature` 的"
        "**域均标量代表温**（不分深），故「浅端(低温)门减弱更多」的逐深 τy(T) 叙事"
        "与实现不符、已撤回。实际机制：温度效应经**时间次序**传递（早期均温低 → "
        "水泥 τy(T) 低 → 屈服门弱 → 先到位相受益），方向与公式行为相容；"
        "本行只作 (Δ深 − Δ浅) 的**方向记录**，实测不符即如实记「否」并入总结论。\n"
        "\n"
        "- **L3-d** 预期挂在**屈服门**（计划 §7 L3 主张）：T-off 门吃 loader 常数 "
        "τy(领13/尾14 Pa)，T-on 吃公式 τy(T)≈5~7 Pa ⇒ 门全域减弱；"
        "泥浆族井内 T≥98.7 °C 恒 clamp@80 °C ⇒ 泥浆 τy 无深度差异；"
        "冻结集中窄边 ⇒ 门减弱窄边受益最大 ⇒ Δη_N − Δη_E > 0。\n"
        "\n"
        "- **L3-b 脚注（组A 的「降」不进本井判定式）**：计划散文「μp 降」仅组A"
        "（指数降）成立；组B 二次式顶点≈24.2 °C，井内温区 [98.7,150.8] °C 已过顶点 ⇒ "
        "升。本井领浆=组B 主导插值、尾浆=纯组B ⇒ 判定式只认**井内升（+）**"
        "（尾浆另加网格单调升）；按 task-3 裁定公式逐字为准，不判关失败。\n")

    # ---- 口径差 ----
    add("## 口径差对照表（§3.7 全网格，公式 vs loader）\n")
    add(f"- 全表见 `{_out_rel()}/口径差对照表.md`。")
    add("- 代表温点关键数字（泥浆族=80 °C 井内 clamp 端；水泥/隔离液=120 °C，P=0.1 MPa）：\n")
    add("| 相 | 分派行 | τy 公式 vs loader | Δ% | μp 公式 vs loader | Δ% |")
    add("|---|---|---|---|---|---|")
    for r in ctx["rep_rows"]:
        if r["T"] is None:
            add(f"| {r['phase']} | {r['row']} | 不替换（沿用 {r['tauy_l']}） | 0% | "
                f"不替换（沿用 {r['mup_l']}） | 0% |")
        else:
            add(f"| {r['phase']} | {r['row']} | {r['tauy_f']:.4f} vs {r['tauy_l']} | "
                f"{_pct(r['tauy_f'], r['tauy_l'])} | {r['mup_f']:.6f} vs {r['mup_l']} | "
                f"{_pct(r['mup_f'], r['mup_l'])} |")
    add("")

    add("## 自审与疑虑\n")
    for note in ctx["notes"]:
        add(f"- {note}")
    add("")
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main(argv: list[str] | None = None) -> int:
    global OUT_DIR
    import argparse
    ap = argparse.ArgumentParser(description="T1-5 三级验证 + 口径差对照表")
    ap.add_argument(
        "--out-dir", default=None,
        help="输出目录（默认 results/温度耦合验证，**该目录为冻结 v1 产物**；"
             "P-1 后重跑请给新目录，如 results/温度耦合验证_P1_20261006）",
    )
    ap.add_argument(
        "--no-pressure", action="store_true",
        help="口径差表退回 v1 的常压 0.1 MPa 缺省口径（默认按 P-1 静压场求值）",
    )
    args = ap.parse_args(argv)
    if args.out_dir:
        OUT_DIR = Path(args.out_dir)
        if not OUT_DIR.is_absolute():
            OUT_DIR = _ROOT / OUT_DIR
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().isoformat(timespec="seconds")
    print(f"== T1-5 三级验证开始 {ts} ==", flush=True)

    # 关2 参照流体（构造一次）
    t_build = time.perf_counter()
    casing_ref, ann_ref, comp_raw, comp60 = build_reference_fluids()
    print(f"[参照] 构造完成 {time.perf_counter() - t_build:.2f}s "
          f"合成隔构={comp_raw.name!r} ρ={comp_raw.density_kg_m3:.1f}", flush=True)

    # 1) 关1 + 关3 对照底：T-off
    run_off = run_variant("T_off", enable_t=False)
    # 2) 关2 T-on Const60
    const_field = ConstantTemperatureField(T_REF)
    run_on60 = run_variant("T_const60", enable_t=True,
                           field_casing=const_field, field_annulus=const_field)
    # 3) 关2 参照（全 T-off + 预派生流体）
    run_ref = run_variant("T_ref_const60", enable_t=False,
                          casing_fluids=casing_ref, annulus_fluids=ann_ref)
    # 4) 关3 全开（管内 T_in / 环空 T_out 瞬态表）
    t_in, t_out = load_delivered_pair()
    run_tab = run_variant("T_table", enable_t=True,
                          field_casing=t_in, field_annulus=t_out)

    # ---------- 关1 ----------
    baseline_text = BASELINE_JSON.read_text(encoding="utf-8")
    baseline = json.loads(baseline_text)
    mine = json.loads(run_off["json_text"])
    d1 = deep_diff(mine, baseline)
    byte_same = (OUT_DIR / "T_off_摘要.json").read_bytes() == BASELINE_JSON.read_bytes()
    key_names = ["effective_efficiency", "eta_narrow", "mixing_index",
                 "buoyancy_number"]
    key_rows = [(k, (baseline.get(k), mine.get(k))) for k in key_names]
    key_rows.append((
        "评价窗效率/CBL评价井段(尾管段)/eta_E",
        (baseline["评价窗效率"]["CBL评价井段(尾管段)"]["eta_E"],
         mine["评价窗效率"]["CBL评价井段(尾管段)"]["eta_E"]),
    ))
    key_rows.append((
        "评价窗效率/CBL评价井段(尾管段)/eta_N",
        (baseline["评价窗效率"]["CBL评价井段(尾管段)"]["eta_N"],
         mine["评价窗效率"]["CBL评价井段(尾管段)"]["eta_N"]),
    ))
    gate1 = {
        "verdict": "过" if (not d1 and byte_same) else "未过",
        "diffs": d1, "n_diff": len(d1), "byte_same": byte_same,
        "key_rows": key_rows,
    }

    # ---------- 关2 ----------
    s_on = dict(run_on60["payload"])
    s_ref = dict(run_ref["payload"])
    audit_on = s_on.pop(AUDIT_KEY, None)
    s_ref.pop(AUDIT_KEY, None)
    d2 = deep_diff(s_on, s_ref)
    key2 = []
    for k in key_names:
        key2.append((k, (run_on60["payload"].get(k), run_ref["payload"].get(k))))
    key2.append(("stop_t_s", (run_on60["stop_t"], run_ref["stop_t"])))
    gate2 = {
        "verdict": "过" if not d2 else "未过",
        "diffs": d2, "n_diff": len(d2), "audit_on": audit_on,
        "key_rows": key2,
    }

    # ---------- 关3 ----------
    tab_payload = dict(run_tab["payload"])
    tab_audit = tab_payload.get(AUDIT_KEY)
    diff_rows = summary_diff_rows(run_off, run_tab)
    l3, aux = l3_checks(run_off, run_tab)
    l3_fail = [c["id"] for c in l3
               if not str(c["一致"]).startswith("是")]
    gate3 = {
        "verdict": "过" if not l3_fail else "未过",
        "n_pass": sum(1 for c in l3 if str(c["一致"]).startswith("是")),
        "n_total": len(l3),
        "fail_ids": "、".join(l3_fail) if l3_fail else "无",
    }

    overall = "全过" if (gate1["verdict"] == "过" and gate2["verdict"] == "过"
                        and gate3["verdict"] == "过") else "未过（原样列数字，不修不调）"

    # ---------- 口径差对照表（v2：P-1 真实静压口径）----------
    p_mpa: float | None = None
    p_source = "压力场未接线（v1 口径）"
    if not args.no_pressure:
        from cemdisp.data.pressure_field import (
            HydrostaticPressureField, insitu_column_density,
        )
        _well_p, _fluids_p, _sched_p, _ = load_ht1_004_tailpipe()
        _rho_bar = insitu_column_density(_fluids_p, _sched_p)
        _pf = HydrostaticPressureField.from_well(_well_p, _rho_bar)
        p_mpa = float(_pf.P(float(_well_p.shoe_md_m), 0.0))
        p_source = (f"HydrostaticPressureField 鞋深单点，ρ̄={_rho_bar:.2f} kg/m³"
                    "（在场相密度按设计泵注体积加权）")
        print(f"[P-1] 呼1-004 鞋深静压 P={p_mpa:.4f} MPa（ρ̄={_rho_bar:.2f}）", flush=True)
    caliber_md, rep_rows = build_caliber_table(p_mpa, p_source)
    (OUT_DIR / "口径差对照表.md").write_text(caliber_md, encoding="utf-8")

    # ---------- 报告 ----------
    notes = [
        "C1 路由修复（2026-10-01 终审）：`_route` 改为 **role==SPACER 优先于名字子串**——"
        "2D 合成等效隔离液（名=「先导浆+隔离液1+隔离液2」、role=SPACER、ρ≈1.82）现走"
        "隔离液族（域外就近借 1.95 式 + borrow/model_assumption 审计）；修复前名字子串"
        "先判 ⇒ 误落泥浆式（τy≈10.98 vs 正确 ≈4.89，T 被误 clamp [40,80]）。"
        "单独「先导浆」(role=WASH)、「平衡液」(role=WASH/MUD) 照旧泥浆式，1D 单相派生不受影响。",
        "关2 豁免键 `" + AUDIT_KEY + "`：T-on 专属审计计数（clamp/borrow/p_default 次数），"
        "不含物理数值；引用两侧其余全部键逐叶子对照。",
        "L3-b 为公式层记录项：计划散文「μp 降」仅组A（顶点≈24 °C 以上组B 升）成立；"
        "本井领/尾浆为组B/插值主导、井内温区 μp 实测升——按 task-3 裁定公式逐字为准，"
        "不判关失败，只如实列预期与实测。",
        "关3 stop_t 两侧可能不同（T-on 改变 casing 前缘时刻 ⇒ annulus_stop_time_s 变）——"
        "为生产语义（各跑自算停算时刻），差异表单列 stop_t_s 行披露。",
        "代表温点口径：泥浆族取 80 °C（井内 T≥98.7 °C 恒 clamp 到域端 80）；"
        "水泥/隔离液取 120 °C（井内温区 [98.7, 150.8] 中段），隔离液 P=0.1 MPa 常压缺省。",
        f"本脚本只新增自身与 `{_out_rel()}/**`；求解器/既有 results 零改动。",
    ]
    ctx = {
        "timestamp": ts,
        "head": _git_head(),
        "gate1": gate1, "gate2": gate2, "gate3": gate3,
        "overall": overall,
        "runs": [
            {"tag": "T_off", "enable_t": False, "field": "无（T-off）",
             "stop_t": run_off["stop_t"], "elapsed": run_off["elapsed"],
             "payload": run_off["payload"]},
            {"tag": "T_const60", "enable_t": True, "field": "Constant(60)",
             "stop_t": run_on60["stop_t"], "elapsed": run_on60["elapsed"],
             "payload": run_on60["payload"]},
            {"tag": "T_ref_const60", "enable_t": False,
             "field": "无（参照=预派生流体）",
             "stop_t": run_ref["stop_t"], "elapsed": run_ref["elapsed"],
             "payload": run_ref["payload"]},
            {"tag": "T_table", "enable_t": True,
             "field": "T_in(管内)/T_out(环空) 瞬态表",
             "stop_t": run_tab["stop_t"], "elapsed": run_tab["elapsed"],
             "payload": run_tab["payload"]},
        ],
        "diff_rows": diff_rows,
        "l3": l3,
        "tab_audit": tab_audit,
        "oob_in": t_in.oob_count,
        "oob_out": t_out.oob_count,
        "t_lo": aux["t_lo"], "t_hi": aux["t_hi"],
        "rep_rows": rep_rows,
        "notes": notes,
    }
    report = write_report(ctx)
    (OUT_DIR / "三级验证报告.md").write_text(report, encoding="utf-8")

    # ---------- 控制台结论（用户检查点） ----------
    print("\n=== 三关结论 ===")
    print(f"关1（T-off=基线逐位）：{gate1['verdict']} —— diff 叶子 {gate1['n_diff']}，"
          f"文件字节{'一致' if gate1['byte_same'] else '不一致'}")
    print(f"关2（T≡60°C 自洽逐位）：{gate2['verdict']} —— diff 叶子 {gate2['n_diff']}"
          f"（豁免审计键）")
    print(f"关3（全开对照+L3 方向）：{gate3['verdict']} —— L3 {gate3['n_pass']}/{gate3['n_total']} 一致"
          + (f"（不一致：{gate3['fail_ids']}）" if gate3["fail_ids"] != "无" else ""))
    print(f"总结论：{overall}")
    print("\n=== 口径差对照表代表温点（公式 vs loader）===")
    for r in rep_rows:
        if r["T"] is None:
            print(f"{r['phase']}: 不替换 Δ=0%（沿用 τy={r['tauy_l']} / μp={r['mup_l']}）")
        else:
            print(f"{r['phase']} @{r['T']:.0f}°C: "
                  f"τy {r['tauy_f']:.4f} vs {r['tauy_l']} ({_pct(r['tauy_f'], r['tauy_l'])}) / "
                  f"μp {r['mup_f']:.6f} vs {r['mup_l']} ({_pct(r['mup_f'], r['mup_l'])})")
            if r.get("tauy_f_p") is not None:
                print(f"    @P={r['p_mpa']:.4f} MPa: "
                      f"τy {r['tauy_f_p']:.4f} ({_pct(r['tauy_f_p'], r['tauy_l'])}) / "
                      f"μp {r['mup_f_p']:.6f} ({_pct(r['mup_f_p'], r['mup_l'])})")
    print(f"\n产物：{OUT_DIR}")
    return 0 if overall == "全过" else 1


if __name__ == "__main__":
    raise SystemExit(main())
