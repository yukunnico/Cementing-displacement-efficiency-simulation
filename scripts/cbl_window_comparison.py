"""B1: 模型评价窗 × 现场 CBL 窗口对照模块（论文素材，8 井）。

背景（模型 vs CBL 差异归因，2026-08-20 调研结论）：模型 η_E 为全井段体积分数均值，
CBL 合格率为评价段胶结等级统计，二者存在三重结构性错配（全井均值 vs 评价段、
泵停时刻 vs 候凝 24h+、体积分数 vs 胶结等级）。本脚本把对照下沉到窗口粒度：
对 WellSpec.evaluation_windows 中每个窗口取模型 per-window η_E/η_N
（M0 评价窗效率，res.summary["评价窗效率"]），与 cbl_evaluation.csv 中该窗口
重叠的 CBL 数字段/定性段配对，产出"窗口对照表"。

可比性三档判定规则（数字段=有数值合格率的 CBL 行；定性段=仅有质量等级文本的行）：
- 可比：窗口与某数字段高度重合（窗口覆盖率≥0.8 且 CBL 段覆盖率≥0.8），
  且窗口内其它数字段贡献可忽略（<0.1 窗口长度）——单一口径直接对照；
- 部分可比：有数字段重叠但不满足"可比"（部分重叠/窗口跨多个数字段口径），
  或仅有定性段且窗口覆盖率≥0.5（定性对照）；
- 不可比：无重叠 CBL 段、或重叠过小（全部窗口覆盖率<0.3）、或该井无 CBL 数据。

已知错配特例（备注列自动标注）：
- hu101：CBL 62.77% 数字段 5390–7810 含 5390–5699.8 双层套管段，模型窗
  "CBL评价井段(单层套管可评价段)" 5699.8–7810 已剔除之（交集口径）；
- ht1_004：全域 29.99%(11–7581) 与尾管段 0.3%(5245–7581) 双口径，模型域
  5245–7660 与全域口径大面积错配，窗口对照同时输出两行。

用法：
  python scripts/cbl_window_comparison.py --well hu101            # 单井：重跑窗口指标+出表
  python scripts/cbl_window_comparison.py --all                   # 8 井全跑（默认目标）
  python scripts/cbl_window_comparison.py --tables-only           # 只用已落盘窗口指标重出对照表
  python scripts/cbl_window_comparison.py --well hu101 --force    # 强制重算该井窗口指标

输出：
  results/最终基线_2026-08-29/cfl_on/单井结果/<well>_windows.json   # 窗口指标（断点续跑）
  results/最终基线_2026-08-29/CBL窗口对照/<well>_窗口对照.csv        # 每井窗口对照表
  results/最终基线_2026-08-29/CBL窗口对照总表.csv                    # 8 井汇总
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
import time
from dataclasses import replace
from pathlib import Path
from typing import cast

import numpy as np

_SCRIPTS_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SCRIPTS_DIR.parent
for _p in (str(_SCRIPTS_DIR), str(_PROJECT_ROOT)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from cemdisp.data.well_spec import DepthValuePoint, WellSpec  # noqa: E402
from cemdisp.models2d import AnnulusD2DGASolver  # noqa: E402
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider  # noqa: E402
from cemdisp.transport1d import CasingFlowSolver  # noqa: E402
from rerun_all_wells_corrected import CORRECTED_KW, NZ, WELLS, _stop_t, _total_t  # noqa: E402

BASELINE_DIR = _PROJECT_ROOT / "results" / "最终基线_2026-08-29"
CFL_ON_DIR = BASELINE_DIR / "cfl_on"
RESULTS_DIR = CFL_ON_DIR / "单井结果"
COMP_DIR = BASELINE_DIR / "CBL窗口对照"
FIELD_DIR = _PROJECT_ROOT / "参考文档" / "现场资料提取"
TOTAL_CSV = BASELINE_DIR / "CBL窗口对照总表.csv"

WELL_NAME_CN = {
    "hu101": "呼101", "hu102": "呼102", "hu103": "呼103", "hu1": "呼探1",
    "hu2": "呼探1-002", "ht1_001": "呼探1-001", "ht1_003": "呼1-003", "ht1_004": "呼1-004",
}


# ---------------------------------------------------------------- 窗口指标重跑

def run_well_windows(well_id: str) -> dict:
    """以最终基线同款配置（CORRECTED_KW, nz=250, CFL on）重跑单井，提取评价窗效率。

    与 rerun_all_wells_corrected.run_one 唯一区别：保留 res.summary["评价窗效率"]
    （M0 per-window b 加权 2D 积分），并记录模型域与一致性校验。
    """
    loader = next(l for w, l in WELLS if w == well_id)
    t0 = time.perf_counter()
    well, fluids, schedule, _ = loader()
    cr = CasingFlowSolver(enable_gravity=True).run(well, fluids, schedule)
    inlet = build_coupled_annulus_inlet_provider(
        cr, CasingFlowSolver(enable_gravity=True), fluids, split_cement_phases=True)
    tt = min(_total_t(schedule) + 1200.0, _stop_t(cr, fluids) + 600.0)
    res = AnnulusD2DGASolver(
        total_t=tt, nz=NZ, enable_cfl_adaptive=True, **CORRECTED_KW
    ).run(well, fluids, inlet)
    fr = cast(dict, res.summary["最终结果"])
    eta_e = float(fr["全井段最终有效顶替效率"])

    # 一致性校验：同配置同数据应为确定性重放，与 A3 基线 JSON 逐位一致
    base_path = RESULTS_DIR / f"{well_id}_corrected_on.json"
    matches, base_eta_e = None, None
    if base_path.exists():
        base = json.loads(base_path.read_text(encoding="utf-8"))
        base_eta_e = float(base["eta_E"])
        matches = abs(base_eta_e - eta_e) < 1e-9

    windows = {}
    for name, w in res.summary["评价窗效率"].items():
        windows[name] = {
            "window_type": w["window_type"],
            "eta_E": float(w["eta_E"]),
            "eta_N": float(w["eta_N"]),
        }
    # 窗口深度区间从 WellSpec 补齐（summary 里没有，对照表需要）
    for w in well.evaluation_windows:
        if w.name in windows:
            windows[w.name]["top_md_m"] = float(w.top_md_m)
            windows[w.name]["bottom_md_m"] = float(w.bottom_md_m)

    return {
        "well": well_id,
        "well_name_cn": WELL_NAME_CN.get(well_id, well_id),
        "eta_E": eta_e,
        "eta_N": float(fr.get("窄四分位效率", float("nan"))),
        "model_domain_md": [float(well.top_md_m), float(well.bottom_md_m)],
        "windows": windows,
        "matches_baseline": matches,
        "baseline_eta_E": base_eta_e,
        "config": {"corrected": True, "cfl_mode": "on", "nz": NZ, **CORRECTED_KW},
        "git_commit": _git_commit(),
        "elapsed_s": round(time.perf_counter() - t0, 1),
    }


def _git_commit() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=_PROJECT_ROOT,
            capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return "unknown"


def _load_or_compute_windows(well_id: str, *, force: bool = False) -> dict:
    """窗口指标落盘/断点续跑：已存在 JSON 直接复用，除非 --force。"""
    path = RESULTS_DIR / f"{well_id}_windows.json"
    if not force and path.exists():
        data = json.loads(path.read_text(encoding="utf-8"))
        print(f"  [复用] {path.name}（eta_E={data['eta_E']:.4f}，"
              f"matches_baseline={data.get('matches_baseline')}）", flush=True)
        return data
    print(f"  [重跑] {well_id}（CORRECTED_KW, nz={NZ}, CFL on）...", flush=True)
    data = run_well_windows(well_id)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    flag = "" if data["matches_baseline"] else "  ⚠️ 与基线 JSON 不一致！"
    print(f"  [落盘] {path.name} eta_E={data['eta_E']:.4f} ({data['elapsed_s']}s){flag}", flush=True)
    return data


# ---------------------------------------------------------------- CBL csv 读取

# 各井 cbl_evaluation.csv 表头不一（4 种 schema），按字段名族映射
_COL_TOP = ("cbl_top_md_m", "md_top_m", "depth_top_m")
_COL_BOT = ("cbl_bottom_md_m", "md_bottom_m", "depth_bottom_m")
_COL_PASS = ("cbl_pass_rate",)
_COL_CLASS = ("cbl_quality_class", "bond_quality", "qualified", "quality_grade")
_COL_INCLUDE = ("include_in_validation",)
_COL_TYPE = ("data_type",)
_COL_NOTE = ("notes", "interpretation_summary", "vd1_description")


def _pick(row: dict, names: tuple[str, ...]) -> str:
    for n in names:
        v = row.get(n)
        if v not in (None, ""):
            return str(v).strip()
    return ""


def _parse_pass_rate(row: dict) -> float | None:
    """数值合格率解析：cbl_pass_rate 直接解析；hu102 综合行合格率藏在
    cbl_amplitude_pct（'66.65%合格率'），仅对该列做正则兜底（不扫全行，
    避免把备注里引用的'62.77%合格率'等叙述性文本误判为该段数字）。"""
    raw = _pick(row, _COL_PASS)
    if raw:
        try:
            return float(raw)
        except ValueError:
            pass
    m = re.search(r"(\d+(?:\.\d+)?)\s*%合格率", row.get("cbl_amplitude_pct", "") or "")
    return float(m.group(1)) if m else None


def read_cbl_rows(well_id: str) -> list[dict]:
    """读取该井 cbl_evaluation.csv，归一化为统一行结构。

    返回行字段：top_md_m / bottom_md_m / pass_rate(数字或 None) / quality_text /
    include_in_validation / data_type。well_id 目录不存在或无 csv 返回 []。
    """
    dirs = sorted(FIELD_DIR.glob(f"{well_id}_*"))
    if not dirs:
        return []
    path = dirs[0] / "cbl_evaluation.csv"
    if not path.exists():
        return []
    rows: list[dict] = []
    with open(path, encoding="utf-8-sig") as fh:
        for raw in csv.DictReader(fh):
            top = _pick(raw, _COL_TOP)
            bot = _pick(raw, _COL_BOT)
            if not top or not bot:
                continue  # 表头行以外的残缺行
            try:
                t, b = float(top), float(bot)
            except ValueError:
                continue
            if b <= t:
                continue
            inc_raw = _pick(raw, _COL_INCLUDE)
            rows.append({
                "top_md_m": t,
                "bottom_md_m": b,
                "pass_rate": _parse_pass_rate(raw),
                "quality_text": _pick(raw, _COL_CLASS) or _pick(raw, _COL_NOTE)[:60],
                "include_in_validation": inc_raw if inc_raw else "1",
                "data_type": _pick(raw, _COL_TYPE),
            })
    return rows


# ---------------------------------------------------------------- 窗口对照

def _overlap(a_top: float, a_bot: float, b_top: float, b_bot: float) -> float:
    """区间交集长度；无交集返回 0。"""
    return max(0.0, min(a_bot, b_bot) - max(a_top, b_top))


def compare_window(win: dict, cbl_rows: list[dict]) -> list[dict]:
    """单窗口 × 全部 CBL 行配对，产出对照行（可比性三档）。

    数字段每条重叠都输出一行（保双口径，如 ht1_004 全域 vs 尾管段）；
    定性段仅取交集最长的一条（明细并入备注）。窗口无任何重叠时输出一行空对照。
    """
    w_top, w_bot = win["top_md_m"], win["bottom_md_m"]
    w_len = w_bot - w_top
    numeric, qualitative = [], []
    for r in cbl_rows:
        inter = _overlap(w_top, w_bot, r["top_md_m"], r["bottom_md_m"])
        if inter <= 0.0:
            continue
        entry = {
            **r,
            "inter_m": inter,
            "win_cov": inter / w_len if w_len > 0 else 0.0,
            "cbl_cov": inter / (r["bottom_md_m"] - r["top_md_m"]),
        }
        (numeric if r["pass_rate"] is not None else qualitative).append(entry)

    rows: list[dict] = []
    if numeric:
        numeric.sort(key=lambda e: e["cbl_cov"], reverse=True)
        # 主口径段：唯一覆盖窗口≥80% 且被窗口覆盖≥80% 的数字段（至多一条）。
        # 存在主口径时直接判"可比"，其余数字段降为"次要口径"行（不再触发跨段降级——
        # 局部小段口径并存是常态，见 hu101 62.77% 总口径 vs 81.1/83.6 局部段）。
        primary = next((e for e in numeric if e["win_cov"] >= 0.8 and e["cbl_cov"] >= 0.8), None)
        for i, e in enumerate(numeric):
            if e is primary:
                others = len(numeric) - 1
                tag = "可比"
                note = f"另有{others}个局部数字段口径并存" if others else ""
            elif max(x["win_cov"] for x in numeric) < 0.3:
                tag, note = "不可比", "重叠过小"
            else:
                tag = "部分可比"
                if primary is not None:
                    note = "次要口径"
                elif len(numeric) > 1:
                    note = "窗口跨多个CBL数字段口径（无单一总口径）"
                else:
                    note = "部分重叠"
            rows.append(_mk_row(win, e["top_md_m"], e["bottom_md_m"], f"{e['pass_rate']:g}",
                                e["quality_text"], e["win_cov"], e["cbl_cov"],
                                e["include_in_validation"], tag, note))
    elif qualitative:
        qualitative.sort(key=lambda e: e["inter_m"], reverse=True)
        e = qualitative[0]
        if e["win_cov"] >= 0.5:
            detail = "；".join(f"{x['quality_text']}:{x['top_md_m']:g}-{x['bottom_md_m']:g}"
                               for x in qualitative[:6])
            rows.append(_mk_row(win, e["top_md_m"], e["bottom_md_m"], "",
                                e["quality_text"], e["win_cov"], e["cbl_cov"],
                                e["include_in_validation"], "部分可比(定性)",
                                f"无定量合格率；窗口内定性段：{detail}"))
        else:
            rows.append(_mk_row(win, e["top_md_m"], e["bottom_md_m"], "",
                                e["quality_text"], e["win_cov"], e["cbl_cov"],
                                e["include_in_validation"], "不可比",
                                "仅定性段且重叠过小"))
    else:
        rows.append(_mk_row(win, None, None, "", "无重叠CBL段", None, None, "",
                            "不可比", "窗口内无任何 CBL 评价段"))
    return rows


def _mk_row(win: dict, c_top, c_bot, pass_rate, quality, win_cov, cbl_cov, inc, tag, note):
    return {
        "窗口名": win["name"],
        "窗口类型": win["window_type"],
        "窗口顶md": f"{win['top_md_m']:g}",
        "窗口底md": f"{win['bottom_md_m']:g}",
        "模型窗口η_E": _fmt(win.get("eta_E")),
        "模型窗口η_N": _fmt(win.get("eta_N")),
        "CBL段顶md": "" if c_top is None else f"{c_top:g}",
        "CBL段底md": "" if c_bot is None else f"{c_bot:g}",
        "CBL合格率%": pass_rate,
        "CBL质量口径": quality,
        "窗口覆盖率": _fmt(win_cov),
        "CBL段覆盖率": _fmt(cbl_cov),
        "纳入现场验证": inc,
        "可比性": tag,
        "备注": note,
    }


def _fmt(v) -> str:
    return "" if v is None else f"{v:.4f}"


# ---------------------------------------------------------------- 主流程

CSV_FIELDS = ["井号", "井名",
    "窗口名", "窗口类型", "窗口顶md", "窗口底md", "模型窗口η_E", "模型窗口η_N",
    "CBL段顶md", "CBL段底md", "CBL合格率%", "CBL质量口径", "窗口覆盖率",
    "CBL段覆盖率", "纳入现场验证", "可比性", "备注"]


def build_well_table(well_id: str, data: dict, cbl_rows: list[dict]) -> list[dict]:
    """单井窗口对照表：每窗口至少一行。"""
    out: list[dict] = []
    for name, win in data["windows"].items():
        win = {**win, "name": name}
        for row in compare_window(win, cbl_rows):
            out.append({"井号": well_id, "井名": data["well_name_cn"], **row})
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="B1: 模型评价窗 × CBL 窗口对照（8 井）")
    ap.add_argument("--well", help="单井 well_id（如 hu101）")
    ap.add_argument("--all", action="store_true", help="8 井全跑")
    ap.add_argument("--tables-only", action="store_true",
                    help="只用已落盘 <well>_windows.json 重出对照表（不重跑）")
    ap.add_argument("--force", action="store_true", help="强制重算窗口指标")
    args = ap.parse_args()

    wells = [args.well] if args.well else ([w for w, _ in WELLS] if args.all or args.tables_only
                                           else None)
    if not wells:
        ap.error("需指定 --well <id> / --all / --tables-only 之一")

    all_rows: list[dict] = []
    for well_id in wells:
        print(f"\n=== {well_id} ===", flush=True)
        cbl_rows = read_cbl_rows(well_id)
        if not cbl_rows:
            print("  ⚠️ 无 cbl_evaluation.csv（或无有效行），全部窗口标不可比", flush=True)
        if args.tables_only:
            path = RESULTS_DIR / f"{well_id}_windows.json"
            data = json.loads(path.read_text(encoding="utf-8"))
        else:
            data = _load_or_compute_windows(well_id, force=args.force)
        rows = build_well_table(well_id, data, cbl_rows)
        COMP_DIR.mkdir(parents=True, exist_ok=True)
        per_well = COMP_DIR / f"{well_id}_窗口对照.csv"
        _write_csv(per_well, rows)
        n_cmp = sum(1 for r in rows if r["可比性"] == "可比")
        print(f"  对照表: {per_well.relative_to(_PROJECT_ROOT)} "
              f"（{len(rows)} 行，其中可比 {n_cmp}）", flush=True)
        all_rows.extend(rows)

    _write_csv(TOTAL_CSV, all_rows)
    print(f"\n总表: {TOTAL_CSV.relative_to(_PROJECT_ROOT)}（{len(all_rows)} 行）", flush=True)


def _write_csv(path: Path, rows: list[dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
