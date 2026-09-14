# -*- coding: utf-8 -*-
"""8-well runner batch driver (stop-flag-fixed rerun, 2026-09-02).

Subprocess-runs each well runner, times it, validates summary JSON,
keeps resumable status.json, then builds comparison table + report.

Usage:
  python run_all.py            # batch run (skip done wells)
  python run_all.py --force    # rerun all
"""
import csv
import json
import os
import subprocess
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]

WELLS = ["hu101", "hu102", "hu103", "hu1", "hu2", "ht1_001", "ht1_003", "ht1_004"]
CN = {
    "hu101": "呼101", "hu102": "呼102", "hu103": "呼103", "hu1": "呼探1",
    "hu2": "呼探1-002", "ht1_001": "呼探1-001", "ht1_003": "呼1-003", "ht1_004": "呼1-004",
}

SUMMARY_REL = {
    "hu101":   "results/呼101尾管_1D2D耦合模型/呼101尾管_1D2D耦合模型_结果摘要.json",
    "hu102":   "results/呼102尾管_1D2D耦合模型/呼102尾管_1D2D耦合模型_结果摘要.json",
    "hu103":   "results/呼103尾管_1D2D耦合模型/呼103尾管_1D2D耦合模型_结果摘要.json",
    "hu1":     "results/呼探1尾管_1D2D耦合模型/呼探1尾管_1D2D耦合模型_结果摘要.json",
    "hu2":     "results/呼探1-002尾管_1D2D耦合模型/呼探1-002尾管_1D2D耦合模型_结果摘要.json",
    "ht1_001": "results/呼探1-001尾管_1D2D耦合模型/呼探1-001尾管_1D2D耦合模型_结果摘要.json",
    "ht1_003": "results/呼1-003_1D2D耦合模型/呼1-003_1D2D耦合模型_结果摘要.json",
    "ht1_004": "results/呼1-004_1D2D耦合模型/呼1-004_1D2D耦合模型_结果摘要.json",
}

KEY_E = "全井段最终有效顶替效率"
KEY_OCC = "最终水泥浆占据率"
KEY_CH = "最终窜槽指数"
KEY_MX = "最终混浆指数"
KEY_UF = "最终失稳指数"
KEY_CBL = "CBL评价井段(单层套管可评价段)"
STOP_NOTE = "runner基线停止口径=尾浆全部进入环空(cement_end_time_s)，封顶泵注结束"


def read_summary(well):
    p = REPO / SUMMARY_REL[well]
    try:
        with open(p, encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return None


def read_metrics(well):
    s = read_summary(well)
    if s is None:
        return None
    fin = s.get("最终结果", {}) or {}
    cbl = (s.get("评价窗效率", {}) or {}).get(KEY_CBL, {}) or {}
    return {
        "eta_E_全井段": fin.get(KEY_E),
        "占据率": fin.get(KEY_OCC),
        "窜槽指数": fin.get(KEY_CH),
        "混浆指数": fin.get(KEY_MX),
        "失稳指数": fin.get(KEY_UF),
        "eta_E_CBL评价窗": cbl.get("eta_E"),
    }


def fmt(x, nd=4):
    return "-" if x is None else f"{x:.{nd}f}"


def load_status():
    p = HERE / "status.json"
    if p.exists():
        try:
            return json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            return {}
    return {}


def save_status(st):
    (HERE / "status.json").write_text(
        json.dumps(st, ensure_ascii=False, indent=2), encoding="utf-8"
    )


def refresh_new_snapshot():
    """Re-extract new metrics for all finished wells into 新值快照.json."""
    snap_path = HERE / "新值快照.json"
    snap = {}
    if snap_path.exists():
        try:
            snap = json.loads(snap_path.read_text(encoding="utf-8"))
        except Exception:
            snap = {}
    n_ok = 0
    for w in WELLS:
        m = read_metrics(w)
        if m is not None:
            base = snap.get(w, {})
            snap[w] = {"井名": CN[w], "路径": SUMMARY_REL[w], **base, **m}
            n_ok += 1
    snap_path.write_text(json.dumps(snap, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[snap] 新值快照.json refreshed ({n_ok}/8)", flush=True)


def build_compare():
    """Build comparison CSV + markdown report from snapshots and status."""
    old = {}
    try:
        old = json.loads((HERE / "旧值快照.json").read_text(encoding="utf-8"))
    except Exception:
        pass
    new = {}
    try:
        new = json.loads((HERE / "新值快照.json").read_text(encoding="utf-8"))
    except Exception:
        pass
    st = load_status()

    rows = []
    for w in WELLS:
        o = old.get(w, {})
        n = new.get(w, {})
        s = st.get(w, {})
        delta = None
        if o.get("eta_E_全井段") is not None and n.get("eta_E_全井段") is not None:
            delta = n["eta_E_全井段"] - o["eta_E_全井段"]
        rows.append({
            "well": w,
            "旧eta_E": fmt(o.get("eta_E_全井段")),
            "新eta_E": fmt(n.get("eta_E_全井段")),
            "Δ": fmt(delta),
            "旧占据率": fmt(o.get("占据率")),
            "新占据率": fmt(n.get("占据率")),
            "旧窜槽": fmt(o.get("窜槽指数")),
            "新窜槽": fmt(n.get("窜槽指数")),
            "旧混浆": fmt(o.get("混浆指数")),
            "新混浆": fmt(n.get("混浆指数")),
            "旧失稳": fmt(o.get("失稳指数")),
            "新失稳": fmt(n.get("失稳指数")),
            "旧CBL评价窗eta_E": fmt(o.get("eta_E_CBL评价窗")),
            "新CBL评价窗eta_E": fmt(n.get("eta_E_CBL评价窗")),
            "stop口径说明": STOP_NOTE,
            "elapsed_s": s.get("elapsed_s", ""),
            "status": s.get("status", "not_run"),
        })
    csv_path = HERE / "对照表.csv"
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print("[cmp] written duibibiao CSV", flush=True)

    # ---- markdown report, built line by line (short strings only) ----
    md = []
    md.append("# 全井 runner 重跑汇总（停止标志修复后）2026-09-02")
    md.append("")
    md.append("## 口径说明")
    md.append("- 本次重跑为 runner 基线配置：python -m cemdisp.runners.<well>_tailpipe，")
    md.append("  默认参数 nz=250，无 CORRECTED M 系列 kwargs；与生产 RR(corrected) 口径不同。")
    md.append("- 停止标志为修复后口径：annulus_stop_time_s 优先取 cement_end_time_s")
    md.append("  （尾浆全部进入环空时刻），casing_flow.py 尾缘配对 + max 约束 + 泵注结束封顶。")
    md.append("- runner 覆盖 results/呼X*/ 旧目录；6 月旧版在 git HEAD 可恢复（预期）。")
    md.append("- 旧值留档 旧值快照.json；新值 新值快照.json；逐井日志 logs/<well>.log；")
    md.append("  断点状态 status.json（中断后 python run_all.py 续跑）。")
    md.append("- CBL评价窗 eta_E 仅 hu101 旧摘要含此字段，其余井旧摘要无（快照中为 null）。")
    md.append("")
    md.append("## 新旧对照表")
    md.append("")
    md.append("| 井 | 旧eta_E | 新eta_E | Δ | 旧占据率 | 新占据率 | 窜槽(旧→新) | 混浆(旧→新) | 失稳(旧→新) | CBL评价窗eta_E(旧→新) | 耗时s | 状态 |")
    md.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        md.append(
            "| " + CN[r["well"]] + " (" + r["well"] + ") "
            + "| " + r["旧eta_E"] + " | " + r["新eta_E"] + " | " + r["Δ"]
            + " | " + r["旧占据率"] + " | " + r["新占据率"]
            + " | " + r["旧窜槽"] + "→" + r["新窜槽"]
            + " | " + r["旧混浆"] + "→" + r["新混浆"]
            + " | " + r["旧失稳"] + "→" + r["新失稳"]
            + " | " + r["旧CBL评价窗eta_E"] + "→" + r["新CBL评价窗eta_E"]
            + " | " + str(r["elapsed_s"]) + " | " + r["status"] + " |"
        )
    md.append("")
    md.append("## 运行记录")
    md.append("")
    for w in WELLS:
        s = st.get(w, {})
        if s.get("status") == "failed":
            md.append("- " + CN[w] + " (" + w + "): **FAILED**, elapsed=" + str(s.get("elapsed_s")) + "s")
            md.append("  ```text")
            tail = (s.get("error_tail") or "").strip().splitlines()[-8:]
            for ln in tail:
                md.append("  " + ln)
            md.append("  ```")
        elif s.get("status") == "done":
            md.append("- " + CN[w] + " (" + w + "): 成功, 耗时 " + str(s.get("elapsed_min")) + " 分钟, 摘要校验通过")
        else:
            md.append("- " + CN[w] + " (" + w + "): 未运行")
    md.append("")
    done_n = sum(1 for w in WELLS if st.get(w, {}).get("status") == "done")
    fail_n = sum(1 for w in WELLS if st.get(w, {}).get("status") == "failed")
    if done_n == 8:
        md.append("**完成 8/8 井。**")
    else:
        md.append("**完成 " + str(done_n) + "/8 井，失败 " + str(fail_n) + " 井（中断可续跑）。**")
    md_text = "\n".join(md) + "\n"
    (HERE / "汇总.md").write_text(md_text, encoding="utf-8")
    print("[cmp] written huizong MD", flush=True)


def run_one(well):
    st = load_status()
    if st.get(well, {}).get("status") == "done" and "--force" not in sys.argv:
        print(f"[skip] {well} already done", flush=True)
        return
    mod = f"cemdisp.runners.{well}_tailpipe"
    log_path = HERE / "logs" / f"{well}.log"
    log_path.parent.mkdir(exist_ok=True)
    print(f"[run ] {well} start", flush=True)
    t0 = time.time()
    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "utf-8"
    env["PYTHONUTF8"] = "1"
    try:
        proc = subprocess.run(
            [sys.executable, "-m", mod],
            cwd=str(REPO),
            env=env,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=3600,
        )
        elapsed = time.time() - t0
        log_path.write_text(
            (proc.stdout or "") + "\n===== STDERR =====\n" + (proc.stderr or ""),
            encoding="utf-8",
        )
        rc = proc.returncode
        ok = rc == 0
        err_tail = (proc.stderr or "")[-2000:]
        if ok:
            s = read_summary(well)
            if s is None:
                ok = False
                err_tail = "returncode=0 but summary json missing/unparseable"
    except subprocess.TimeoutExpired:
        elapsed = time.time() - t0
        ok = False
        err_tail = f"timeout after {elapsed:.0f}s"
        log_path.write_text(f"TIMEOUT {elapsed:.0f}s", encoding="utf-8")
    except Exception:
        elapsed = time.time() - t0
        ok = False
        err_tail = traceback.format_exc()[-2000:]
        log_path.write_text(err_tail, encoding="utf-8")

    st = load_status()
    st[well] = {
        "status": "done" if ok else "failed",
        "elapsed_s": round(elapsed, 1),
        "elapsed_min": round(elapsed / 60, 1),
        "error_tail": None if ok else err_tail[-2000:],
        "log": str(log_path),
        "finished_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    save_status(st)
    print(f"[{'done' if ok else 'FAIL'}] {well} elapsed={elapsed/60:.1f}min", flush=True)


def main():
    t0 = time.time()
    for w in WELLS:
        run_one(w)
        refresh_new_snapshot()
        build_compare()
    print(f"[all ] total {(time.time()-t0)/60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
