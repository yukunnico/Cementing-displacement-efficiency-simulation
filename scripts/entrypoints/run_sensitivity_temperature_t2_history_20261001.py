# -*- coding: utf-8 -*-
"""T2 温压敏感性·过程观察层批（2026-10-01）——时程 + 事件时刻 + 等库存比截面。

用户裁定（2026-10-01，方案 1+2 组合，本脚本只照单执行）
------------------------------------------------------
1. **方案 ①**：每变体出**时程表**（快照网格，降采样 ≤300 行）+ **事件时刻表**
   （隔离液出鞋 / 尾浆出鞋 / 泵注 50% / 碰压前，每事件一行含全部指标）。
2. **方案 ②**：每变体出**等库存比截面表**（cement_occ 首达 0.5 / 0.8 / 1.0
   的时刻与指标；未达档如实标「未达到」并给末态对照）。
3. **只加观察层，不碰 solver**：时程量全部来自既有
   ``res.metrics``（逐行）与 ``res.*_snapshots``（每 60 步 + 末步快照）的
   **跑后后处理**；共享装配层仅把 ``run_variant`` 内核抽成
   ``run_variant_res``（逐行搬移、数值路径零改动），``run_variant`` 保持薄壳。
4. **T-off 逐位红线**：solver summary 不加任何键；本批每个变体另存
   ``*_结果摘要.json``（= res.summary 原样），验收时与
   ``敏感性变体_温压T2_2026-10-01/`` 同名文件数值全等对照。
5. **输出新目录** ``results/敏感性变体_温压T2时程_2026-10-01/``；
   三旧敏感性目录 + T2 45 变体目录**只读冻结**（本脚本有目录守卫）。
6. 配置（run_opts / 温度档 / 变体名）与 T2 批**同名同款**，只加观察层。

变体矩阵（13，注册表直接复用 T2 的 ``build_temperature_variants()``）
---------------------------------------------------------------------
- 呼1-004（5）：``Toff_zero``、``Ton_static_rate_x1.0``、
  ``Ton_static_rate_x1.4``、``Ton_const60_rate_x1.0``、``Ton_table_rate_x1.0``
- 呼101 / 呼103 各（4）：``Toff_zero``、``Ton_static_rate_x1.0``、
  ``Ton_static_rate_x1.4``、``Ton_const60_rate_x1.0``

时程 10 指标（列序即写盘序）
----------------------------
t_s, η_E, η_N, 饥饿份额, cement_occ, front_narrow_m, front_wide_m,
宽窄边差_m, 屈服门活化率_b加权, 屈服门_wall占比

- η_E / cement_occ / 前缘：来自 ``res.metrics``（逐时间步，恒等式
  η_E ≡ effective_efficiency ≡ bulk_cement_fill = cement_occ）；
- η_N / 饥饿份额 / 屈服门活化率双口径：**快照后处理**
  （η_N = ``_narrow_quarter_efficiency``，饥饿 = ``starved_volume_fraction``，
  门 = b 加权 mean(wall) 与 wall>0 占比，口径同 T2 末态判别量）。
- 事件时刻取值：metrics 列在全逐行网格上线性插值、快照列在快照网格上
  线性插值；落在仿真区间外则端点钳位并在「采样」列标注。

输出
----
results/敏感性变体_温压T2时程_2026-10-01/
    {井}_{变体}_结果摘要.json     res.summary 原样（与 T2 同 schema）
    {井}_{变体}_时程.csv          过程观察时程（≤300 行，UTF-8-sig）
    {井}_{变体}_事件时刻表.csv     4 事件 × 全指标
    {井}_{变体}_等库存比截面.csv    3 档库存比 × 全指标
    事件时刻汇总_温压T2时程.csv/.md
    等库存比截面汇总_温压T2时程.csv/.md
    （图：scripts/plots/plot_sensitivity_temperature_history_20261001.py）

用法
----
    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 conda run -n cementT \
        python scripts/entrypoints/run_sensitivity_temperature_t2_history_20261001.py
    # 子集：--well 呼1-004 --variant Toff_zero；强制重算：--force
"""
from __future__ import annotations

import argparse
import csv
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
_SCRIPTS_DIR = PROJECT_ROOT / "scripts"
for p in (str(PROJECT_ROOT), str(_SCRIPTS_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from entrypoints.rerun_all_wells_corrected import _stop_t, _total_t  # noqa: E402
from entrypoints.run_sensitivity_current_20260916 import run_variant_res  # noqa: E402
from entrypoints.run_sensitivity_temperature_t2_20261001 import (  # noqa: E402
    WELLS,
    _jsonable,
    build_temperature_variants,
)

OUT_DIR = PROJECT_ROOT / "results" / "敏感性变体_温压T2时程_2026-10-01"

# 只读冻结目录（红线：零写入；输出必须落在它们之外）
FROZEN_DIRS = (
    PROJECT_ROOT / "results" / "敏感性变体_当前口径_2026-09-16",
    PROJECT_ROOT / "results" / "敏感性变体_呼1-004_2026-09-27",
    PROJECT_ROOT / "results" / "敏感性变体_runner链_2026-09-27",
    PROJECT_ROOT / "results" / "敏感性变体_温压T2_2026-10-01",
)
for _d in FROZEN_DIRS:
    if OUT_DIR == _d or _d in OUT_DIR.parents:
        raise RuntimeError(f"输出目录防护触发：{OUT_DIR} 落在冻结目录 {_d} 内")

# 变体矩阵（13，名字与 T2 注册表同名同款）
REQUIRED_VARIANTS: dict[str, list[str]] = {
    "呼1-004": [
        "Toff_zero",
        "Ton_static_rate_x1.0",
        "Ton_static_rate_x1.4",
        "Ton_const60_rate_x1.0",
        "Ton_table_rate_x1.0",
    ],
    "呼101": [
        "Toff_zero",
        "Ton_static_rate_x1.0",
        "Ton_static_rate_x1.4",
        "Ton_const60_rate_x1.0",
    ],
    "呼103": [
        "Toff_zero",
        "Ton_static_rate_x1.0",
        "Ton_static_rate_x1.4",
        "Ton_const60_rate_x1.0",
    ],
}
TOTAL_VARIANTS = sum(len(v) for v in REQUIRED_VARIANTS.values())  # 13

# 时程 10 指标（列序即写盘序）
METRIC_COLS = [
    "t_s", "η_E", "η_N", "饥饿份额", "cement_occ",
    "front_narrow_m", "front_wide_m", "宽窄边差_m",
    "屈服门活化率_b加权", "屈服门_wall占比",
]
# 汇总表在指标前多出的识别列
ID_COLS = ["井名", "变体", "状态"]

MAX_HISTORY_ROWS = 300
OCC_LEVELS = (0.5, 0.8, 1.0)
# 事件时刻表行序
EVENT_ORDER = ("隔离液出鞋", "尾浆出鞋", "泵注50%", "碰压前")

# 1D 鞋口相名兜底分类（FluidSpec.role 优先，名字子串回退；与 casing_depth_profile 同源）
_SPACER_NAME_HINTS = ("隔离", "平衡", "先导", "冲洗")
_TAIL_NAME_HINTS = ("尾浆", "尾管水泥")


def _trapez2d(field, geom):
    from cemdisp.models2d.annulus_d2dga import _trapez2d as _t
    return float(_t(np.asarray(field, dtype=float), geom))


# ---------------------------------------------------------------- 观察层核心
def build_history(res) -> tuple[list[dict], dict]:
    """快照网格上的 10 列时程 + QC 字典（**纯后处理，不进 solver**）。

    网格 = ``res.snapshot_times_s``（每 save_interval=60 步 + 末步）。
    metrics 列在该时刻查表取值（快照与指标同轮同 record_time ⇒ 精确对齐）；
    快照列逐快照后处理。QC：快照口径 cement_occ vs metrics bulk_cement_fill、
    末态 η_N vs summary 窄四分位效率。
    """
    from cemdisp.diagnostics.internal_consistency import starved_volume_fraction
    from cemdisp.diagnostics.displacement_metrics import _narrow_quarter_efficiency

    geom = res.geom
    b = np.asarray(geom["b"], dtype=float)
    denom = max(_trapez2d(b, geom), 1e-12)

    times = np.asarray(res.snapshot_times_s, dtype=float)
    if times.size == 0:
        raise RuntimeError("res.snapshot_times_s 为空——快照观察层无源数据")
    m = res.metrics
    mt = np.asarray(m["time_s"], dtype=float)
    pos = np.searchsorted(mt, times)
    if not np.allclose(mt[np.clip(pos, 0, len(mt) - 1)], times, atol=1e-6):
        raise RuntimeError("快照时刻未能在 metrics 中精确对齐（record_time 应同轮）")

    col = {name: np.asarray(m[name], dtype=float) for name in (
        "effective_efficiency", "bulk_cement_fill", "front_narrow_m", "front_wide_m",
    )}

    rows: list[dict] = []
    occ_snap = np.empty(len(times), dtype=float)
    eta_n = np.empty(len(times), dtype=float)
    for i, t in enumerate(times):
        cement = np.asarray(res.cement_snapshots[i], dtype=float)
        wall = np.asarray(res.wall_snapshots[i], dtype=float)
        eta_e = float(col["effective_efficiency"][pos[i]])
        occ = float(col["bulk_cement_fill"][pos[i]])
        fn = float(col["front_narrow_m"][pos[i]])
        fw = float(col["front_wide_m"][pos[i]])
        eta_n_i = _narrow_quarter_efficiency(cement, geom)
        starved = starved_volume_fraction(cement, geom)
        occ_i = _trapez2d(b * cement, geom) / denom
        gate_b = _trapez2d(b * wall, geom) / denom
        gate_frac = float(np.mean(wall > 0.0))
        occ_snap[i] = occ_i
        eta_n[i] = eta_n_i
        rows.append({
            "t_s": float(t),
            "η_E": eta_e,
            "η_N": eta_n_i,
            "饥饿份额": starved,
            "cement_occ": occ,
            "front_narrow_m": fn,
            "front_wide_m": fw,
            "宽窄边差_m": fw - fn,
            "屈服门活化率_b加权": gate_b,
            "屈服门_wall占比": gate_frac,
        })

    final = res.summary["最终结果"]
    qc = {
        "快照数": int(len(times)),
        "metrics行数": int(len(m)),
        "occ快照vs逐行_maxdiff": float(np.max(np.abs(occ_snap - col["bulk_cement_fill"][pos]))),
        "末态η_N_快照vs summary": abs(
            eta_n[-1] - float(final["窄四分位效率"])
        ),
        "末态η_E_快照vs summary": abs(
            float(col["effective_efficiency"][pos[-1]])
            - float(final["全井段最终有效顶替效率"])
        ),
        "stop_t_s": float(mt[-1]),
    }
    return rows, qc


def _interp(t: float, xs: np.ndarray, ys: np.ndarray) -> tuple[float, str]:
    """线性插值取值；t 超出网格范围 ⇒ 端点钳位并标注。"""
    if len(xs) == 0:
        return float("nan"), "空网格"
    if t < xs[0]:
        return float(ys[0]), "端点钳位(前)"
    if t > xs[-1]:
        return float(ys[-1]), "端点钳位(后)"
    return float(np.interp(t, xs, ys)), "线性插值"


def build_series_sources(res, hist_rows: list[dict]) -> dict:
    """两套插值源：metrics 逐行网格（4 列）+ 快照网格（4 列，取自已建时程）。"""
    m = res.metrics
    src: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for key, col in (
        ("η_E", "effective_efficiency"),
        ("cement_occ", "bulk_cement_fill"),
        ("front_narrow_m", "front_narrow_m"),
        ("front_wide_m", "front_wide_m"),
    ):
        src[key] = (np.asarray(m["time_s"], dtype=float), np.asarray(m[col], dtype=float))
    src["_history"] = (np.asarray([r["t_s"] for r in hist_rows]), hist_rows)
    return src


def sample_at(t: float, src: dict) -> tuple[dict, str]:
    """在时刻 t 取全部 10 指标（metrics 列 / 快照列分别在各自网格插值）。"""
    out: dict = {"t_s": float(t)}
    flags: list[str] = []
    for key in ("η_E", "cement_occ", "front_narrow_m", "front_wide_m"):
        xs, ys = src[key]
        v, flag = _interp(t, xs, ys)
        out[key] = v
        flags.append(flag)
    out["宽窄边差_m"] = out["front_wide_m"] - out["front_narrow_m"]
    h_t, h_rows = src["_history"]
    for key in ("η_N", "饥饿份额", "屈服门活化率_b加权", "屈服门_wall占比"):
        ys = np.asarray([r[key] for r in h_rows], dtype=float)
        v, flag = _interp(t, h_t, ys)
        out[key] = v
        flags.append(flag)
    status = "网格" if all(f == "线性插值" for f in flags) else (
        "端点钳位" if any("钳位" in f for f in flags) else "混合"
    )
    return out, status


# ---------------------------------------------------------------- 事件/截面
def event_times(cr, schedule2, fluids2, metrics_t: np.ndarray,
                stop_t: float) -> dict[str, tuple[float | None, str]]:
    """4 事件时刻：隔离液出鞋 / 尾浆出鞋 / 泵注 50% / 碰压前。

    - 出鞋：``cr.shoe_timeline.events`` 中相分数 > 0 的首个对应相事件；
      FluidSpec.role 判相，名字子串兜底（同 casing_depth_profile 分类）。
    - 泵注 50%：``0.5 × _total_t(本变体泵序)``（排量缩放变体时长随 1/f 变）。
    - 碰压前：metrics 中 ≤ stop_t 的最后一行采样时刻
      （stop_t = _stop_t = 尾浆全部入库 ≡ 碰压断面；若 stop_t 超出仿真区间
      则取仿真末行）。
    """
    from cemdisp.data.fluid_spec import FluidRole

    roles = {f.name: f.role for f in fluids2}

    def _first_shoe(roles_wanted, name_hints) -> tuple[float | None, str]:
        for ev in cr.shoe_timeline.events:
            for nm, frac in ev.phase_fractions:
                if frac <= 0.0:
                    continue
                role = roles.get(nm)
                hit = (role in roles_wanted) if role is not None else False
                if not hit and role is None:
                    hit = any(h in str(nm) for h in name_hints)
                if hit:
                    return float(ev.time_s), "鞋口事件"
        return None, "未达到（鞋口无该相事件）"

    t_spacer, s_spacer = _first_shoe({FluidRole.SPACER}, _SPACER_NAME_HINTS)
    t_tail, s_tail = _first_shoe({FluidRole.TAIL}, _TAIL_NAME_HINTS)
    t_pump50 = 0.5 * float(_total_t(schedule2))
    bump_idx = np.where(metrics_t <= stop_t + 1e-9)[0]
    if bump_idx.size:
        t_bump = float(metrics_t[bump_idx.max()])
        s_bump = "仿真内最后采样"
    else:
        t_bump = float(metrics_t[0])
        s_bump = "stop_t 早于首个采样（取首行）"

    return {
        "隔离液出鞋": (t_spacer, s_spacer),
        "尾浆出鞋": (t_tail, s_tail),
        "泵注50%": (t_pump50, "泵序解析"),
        "碰压前": (t_bump, s_bump),
    }


def occ_crossings(res, levels=OCC_LEVELS) -> dict[float, tuple[float | None, str]]:
    """cement_occ(t) 首达各档的时刻（逐行网格两采样点间线性插值）。"""
    m = res.metrics
    mt = np.asarray(m["time_s"], dtype=float)
    occ = np.asarray(m["bulk_cement_fill"], dtype=float)
    out: dict[float, tuple[float | None, str]] = {}
    for lv in levels:
        hit = np.where(occ >= lv)[0]
        if hit.size == 0:
            out[lv] = (None, "未达到")
            continue
        i = int(hit[0])
        if i == 0:
            out[lv] = (float(mt[0]), "首采样已达")
        else:
            o0, o1 = float(occ[i - 1]), float(occ[i])
            t0, t1 = float(mt[i - 1]), float(mt[i])
            if o1 - o0 <= 1e-15:
                out[lv] = (t1, "线性插值(退化)")
            else:
                frac = (lv - o0) / (o1 - o0)
                out[lv] = (t0 + frac * (t1 - t0), "线性插值")
    return out


# ---------------------------------------------------------------- 降采样
def downsample(rows: list[dict], force_times: list[float],
               max_rows: int = MAX_HISTORY_ROWS) -> list[dict]:
    """等时距降采样 + 事件/穿越时刻强制保留（线性插值补行），按 t_s 排序去重。"""
    n = len(rows)
    if n > max_rows:
        idx = np.unique(np.linspace(0, n - 1, max_rows).round().astype(int))
        base = [rows[i] for i in idx]
    else:
        base = list(rows)
    t_grid = np.asarray([r["t_s"] for r in rows], dtype=float)
    extra: list[dict] = []
    for ft in force_times:
        if ft is None or not np.isfinite(ft):
            continue
        if np.any(np.isclose(t_grid, ft, atol=1e-9)):
            continue  # 已在网格上
        r = {"t_s": float(ft)}
        for key in METRIC_COLS[1:]:
            r[key] = float(np.interp(ft, t_grid, np.asarray(
                [x[key] for x in rows], dtype=float)))
        extra.append(r)
    merged = base + extra
    merged.sort(key=lambda r: r["t_s"])
    out: list[dict] = []
    for r in merged:
        if out and abs(out[-1]["t_s"] - r["t_s"]) < 1e-9:
            continue
        out.append(r)
    return out


# ---------------------------------------------------------------- 落盘
def _write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def _fmt(v, nd: int = 4) -> str:
    if v is None or v == "" or (isinstance(v, float) and not np.isfinite(v)):
        return "—"
    if isinstance(v, str):
        try:
            v = float(v)          # 落盘 CSV 重读的数值串同样按 nd 位出
        except ValueError:
            return v
    if not np.isfinite(float(v)):
        return "—"
    return f"{float(v):.{nd}f}"


def _md_table(rows: list[dict], cols: list[str], numeric: set[str]) -> list[str]:
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in rows:
        cells = []
        for c in cols:
            v = r.get(c, "")
            cells.append(_fmt(v, 4) if c in numeric else str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return lines


def _write_history_md(well: str, name: str, rows: list[dict]) -> None:
    """时程表的 MD 镜像（与 CSV 同行同值，供 obsidian 直读）。"""
    numeric = set(METRIC_COLS)
    md = [
        f"# {well} × {name} 过程观察时程（T2 时程批 2026-10-01）",
        "",
        "口径：09-16 脚本链；行 = 快照网格（每 60 步 + 末步）"
        "∪ 强制采样行（4 事件时刻 + 库存比穿越时刻，线性插值）；"
        "降采样上限 300 行。η_E/cement_occ/前缘来自逐行 metrics，"
        "η_N/饥饿/屈服门来自快照后处理。",
        "",
        f"共 {len(rows)} 行。",
        "",
    ]
    md += _md_table(rows, METRIC_COLS, numeric)
    (OUT_DIR / f"{well}_{name}_时程.md").write_text("\n".join(md), encoding="utf-8")


# ---------------------------------------------------------------- 主批
def _paths(well: str, variant: str) -> dict[str, Path]:
    stem = OUT_DIR / f"{well}_{variant}"
    return {
        "summary": Path(f"{stem}_结果摘要.json"),
        "history": Path(f"{stem}_时程.csv"),
        "events": Path(f"{stem}_事件时刻表.csv"),
        "occ": Path(f"{stem}_等库存比截面.csv"),
    }


def _complete(well: str, variant: str) -> bool:
    return all(p.exists() for p in _paths(well, variant).values())


def run_one(well: str, loader, name: str, transforms, run_opts: dict) -> dict:
    """跑一个变体并产出时程/事件/截面三表 + summary JSON。返回事件/截面行。"""
    t0 = time.perf_counter()
    well_fn, fluid_fn, sched_fn = transforms
    res, cr, schedule2, extra = run_variant_res(
        loader, well_fn, fluid_fn, sched_fn, run_opts, well_key=well,
    )
    elapsed = round(time.perf_counter() - t0, 1)

    hist_rows, qc = build_history(res)
    src = build_series_sources(res, hist_rows)
    mt = np.asarray(res.metrics["time_s"], dtype=float)
    stop_t = float(extra["stop_t_s"])
    # fluids2 需从 loader 现取（与 run_variant_res 同一变换链，纯读变换无副作用）
    _well0, fluids0, _sched0, _ = loader()
    fluids2 = fluid_fn(fluids0)
    ev_t = event_times(cr, schedule2, fluids2, mt, stop_t)
    crossings = occ_crossings(res)

    # ---- 事件时刻表 ----
    event_rows: list[dict] = []
    for label in EVENT_ORDER:
        t, ev_status = ev_t[label]
        if t is None:
            metrics_last = hist_rows[-1]
            row = {"井名": well, "变体": name, "事件": label, "t_s": "",
                   "t_min": "", "状态": ev_status}
            row.update({k: metrics_last[k] for k in METRIC_COLS[1:]})
        else:
            vals, sample_flag = sample_at(float(t), src)
            row = {"井名": well, "变体": name, "事件": label,
                   "t_s": float(t), "t_min": float(t) / 60.0,
                   "状态": ev_status if sample_flag == "网格"
                   else f"{ev_status}；采样={sample_flag}"}
            row.update({k: vals[k] for k in METRIC_COLS[1:]})
        event_rows.append(row)

    # ---- 等库存比截面表 ----
    occ_rows: list[dict] = []
    for lv in OCC_LEVELS:
        t, status = crossings[lv]
        if t is None:
            last = hist_rows[-1]
            row = {"井名": well, "变体": name, "库存比档": lv, "t_s": "",
                   "t_min": "", "状态": "未达到（末态对照）"}
            row.update({k: last[k] for k in METRIC_COLS[1:]})
            row["cement_occ"] = last["cement_occ"]  # 末态实际库存比（< 档位）
        else:
            vals, _sample_flag = sample_at(float(t), src)
            row = {"井名": well, "变体": name, "库存比档": lv,
                   "t_s": float(t), "t_min": float(t) / 60.0, "状态": status}
            row.update({k: vals[k] for k in METRIC_COLS[1:]})
        occ_rows.append(row)

    # ---- 时程（降采样 + 强制保留事件/穿越时刻）----
    force_times = [t for t, _ in ev_t.values() if t is not None]
    force_times += [t for t, st in crossings.values()
                    if t is not None and st != "未达到"]
    hist_final = downsample(hist_rows, force_times)

    p = _paths(well, name)
    p["summary"].write_text(
        json.dumps(_jsonable(res.summary), ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    _write_csv(p["history"], hist_final, METRIC_COLS)
    ev_fields = ["井名", "变体", "事件", "t_s", "t_min", "状态"] + METRIC_COLS[1:]
    _write_csv(p["events"], event_rows, ev_fields)
    occ_fields = ["井名", "变体", "库存比档", "t_s", "t_min", "状态"] + METRIC_COLS[1:]
    _write_csv(p["occ"], occ_rows, occ_fields)

    print(f"  [完成] {well} × {name}: 快照{qc['快照数']}/"
          f"metrics{qc['metrics行数']} 行 → 时程{len(hist_final)} 行 "
          f"(occ口径差 {qc['occ快照vs逐行_maxdiff']:.2e}, "
          f"末态η_N差 {qc['末态η_N_快照vs summary']:.2e}, {elapsed}s)", flush=True)
    return {"well": well, "variant": name, "events": event_rows,
            "occ": occ_rows, "qc": qc, "elapsed": elapsed}


def rebuild_summaries(results: list[dict]) -> None:
    """从本批已落盘行幂等重建两张汇总表（支持断点/子集）。"""
    ev_all: list[dict] = []
    occ_all: list[dict] = []
    for r in results:
        ev_all.extend(r["events"])
        occ_all.extend(r["occ"])
    if not ev_all:
        return
    ev_fields = ["井名", "变体", "事件", "t_s", "t_min", "状态"] + METRIC_COLS[1:]
    occ_fields = ["井名", "变体", "库存比档", "t_s", "t_min", "状态"] + METRIC_COLS[1:]
    _write_csv(OUT_DIR / "事件时刻汇总_温压T2时程.csv", ev_all, ev_fields)
    _write_csv(OUT_DIR / "等库存比截面汇总_温压T2时程.csv", occ_all, occ_fields)

    numeric = set(METRIC_COLS[1:]) | {"t_s", "t_min"}
    md = [
        "# T2 时程观察·事件时刻汇总（2026-10-01）",
        "",
        "口径：09-16 脚本链（CORRECTED_KW + CFL、tt=min(泵总+1200, stop_t)）；"
        "只加观察层（solver 零改动、summary 零新增键）；**基线 = 本批 `Toff_zero`**。",
        "事件定义：隔离液/尾浆出鞋 = 1D shoe_timeline 首个该相出流事件；"
        "泵注50% = 0.5×本变体泵序总时长；碰压前 = metrics 中 ≤ stop_t 的最后采样。",
        "",
    ]
    md += _md_table(ev_all, ev_fields, numeric)
    (OUT_DIR / "事件时刻汇总_温压T2时程.md").write_text("\n".join(md), encoding="utf-8")

    md2 = [
        "# T2 时程观察·等库存比截面汇总（2026-10-01）",
        "",
        "cement_occ(t)（= bulk_cement_fill，与 η_E 恒等）首达 0.5 / 0.8 / 1.0 "
        "的时刻与该时刻指标；两采样点间线性插值。未达档标「未达到（末态对照）」"
        "并给末态全指标。",
        "",
    ]
    md2 += _md_table(occ_all, occ_fields, numeric)
    (OUT_DIR / "等库存比截面汇总_温压T2时程.md").write_text(
        "\n".join(md2), encoding="utf-8",
    )
    print(f"重建汇总：事件 {len(ev_all)} 行、截面 {len(occ_all)} 行 → {OUT_DIR}",
          flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description="T2 时程观察批（13 变体，断点续跑）")
    ap.add_argument("--well", choices=[*REQUIRED_VARIANTS.keys(), "all"], default="all")
    ap.add_argument("--variant", default=None, help="只跑/只列该变体名")
    ap.add_argument("--force", action="store_true", help="忽略既有产物全部重算")
    args = ap.parse_args()

    registry = build_temperature_variants()
    # 注册表自检：13 变体必须在 T2 注册表中同名存在（配置同款）
    total = 0
    for well, names in REQUIRED_VARIANTS.items():
        have = {v[0] for v in registry[well]}
        missing = [n for n in names if n not in have]
        if missing:
            raise RuntimeError(f"{well}: T2 注册表缺变体 {missing}")
        total += len(names)
    if total != TOTAL_VARIANTS or total != 13:
        raise RuntimeError(f"变体矩阵应为 13，实际 {total}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    well_keys = list(REQUIRED_VARIANTS) if args.well == "all" else [args.well]

    for well in well_keys:
        mod, fn = WELLS[well]
        loader = getattr(importlib.import_module(mod), fn)
        reg = {v[0]: v for v in registry[well]}
        names = REQUIRED_VARIANTS[well]
        if args.variant is not None:
            names = [n for n in names if n == args.variant]
            if not names:
                raise SystemExit(f"{well} 无变体 {args.variant!r}")
        print(f"\n=== {well}  {len(names)} 变体 ===", flush=True)
        for name in names:
            if not args.force and _complete(well, name):
                print(f"  [复用] {well} × {name}", flush=True)
                continue
            _n, well_fn, fluid_fn, sched_fn, run_opts = reg[name]
            assert _n == name
            run_one(well, loader, name, (well_fn, fluid_fn, sched_fn), run_opts)

    # 汇总：一律从落盘 CSV 重读（子集/断点跑也能重建全量，格式唯一）
    all_results = []
    for well in REQUIRED_VARIANTS:
        for name in REQUIRED_VARIANTS[well]:
            if not _complete(well, name):
                continue
            # 断点复用的变体补写时程 MD（跑过的已在 run_one 写过，幂等覆盖）
            hist = _read_history(well, name)
            if hist:
                _write_history_md(well, name, hist)
            all_results.append(_reload_result(well, name))
    rebuild_summaries(all_results)
    return 0


def _read_history(well: str, variant: str) -> list[dict]:
    p = _paths(well, variant)["history"]
    if not p.exists():
        return []
    with p.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def _reload_result(well: str, name: str) -> dict:
    """断点复用：从落盘 CSV 重读事件/截面行（汇总幂等重建用）。"""
    p = _paths(well, name)
    with p["events"].open(encoding="utf-8-sig", newline="") as fh:
        events = list(csv.DictReader(fh))
    with p["occ"].open(encoding="utf-8-sig", newline="") as fh:
        occ = list(csv.DictReader(fh))
    return {"well": well, "variant": name, "events": events, "occ": occ,
            "qc": {}, "elapsed": None}


if __name__ == "__main__":
    raise SystemExit(main())
