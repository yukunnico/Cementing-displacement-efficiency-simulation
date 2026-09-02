# -*- coding: utf-8 -*-
"""提取 8 井结果摘要关键字段快照（跑前旧值 / 跑后新值均可用）。

只读 JSON 并提取标量字段，丢弃 tier0 大数组，避免上下文膨胀。
用法: python extract_snapshot.py <tag>   # tag = old / new
输出: 旧值快照.json 或 新值快照.json（UTF-8）
"""
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
RESULTS = REPO / "results"
HERE = Path(__file__).resolve().parent

WELLS = {
    "hu101":   ("results/呼101尾管_1D2D耦合模型/呼101尾管_1D2D耦合模型_结果摘要.json", "呼101"),
    "hu102":   ("results/呼102尾管_1D2D耦合模型/呼102尾管_1D2D耦合模型_结果摘要.json", "呼102"),
    "hu103":   ("results/呼103尾管_1D2D耦合模型/呼103尾管_1D2D耦合模型_结果摘要.json", "呼103"),
    "hu1":     ("results/呼探1尾管_1D2D耦合模型/呼探1尾管_1D2D耦合模型_结果摘要.json", "呼探1"),
    "hu2":     ("results/呼探1-002尾管_1D2D耦合模型/呼探1-002尾管_1D2D耦合模型_结果摘要.json", "呼探1-002"),
    "ht1_001": ("results/呼探1-001尾管_1D2D耦合模型/呼探1-001尾管_1D2D耦合模型_结果摘要.json", "呼探1-001"),
    "ht1_003": ("results/呼1-003_1D2D耦合模型/呼1-003_1D2D耦合模型_结果摘要.json", "呼1-003"),
    "ht1_004": ("results/呼1-004_1D2D耦合模型/呼1-004_1D2D耦合模型_结果摘要.json", "呼1-004"),
}


def pick(d, *keys):
    cur = d
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return None
        cur = cur[k]
    return cur


def extract(summary: dict) -> dict:
    fin = summary.get("最终结果", {}) or {}
    cbl = (summary.get("评价窗效率", {}) or {}).get("CBL评价井段(单层套管可评价段)", {}) or {}
    return {
        "eta_E_全井段": pick(fin, "全井段最终有效顶替效率"),
        "占据率": pick(fin, "最终水泥浆占据率"),
        "窜槽指数": pick(fin, "最终窜槽指数"),
        "混浆指数": pick(fin, "最终混浆指数"),
        "失稳指数": pick(fin, "最终失稳指数"),
        "窄四分位效率": pick(fin, "窄四分位效率"),
        "eta_E_CBL评价窗": pick(cbl, "eta_E"),
        "eta_N_CBL评价窗": pick(cbl, "eta_N"),
    }


def main() -> None:
    tag = sys.argv[1] if len(sys.argv) > 1 else "old"
    out = {}
    for well, (rel, cn) in WELLS.items():
        p = REPO / rel
        try:
            with open(p, encoding="utf-8") as f:
                summary = json.load(f)
            out[well] = {"井名": cn, "路径": rel, **extract(summary)}
            print(f"[{tag}] {well}: OK")
        except Exception as exc:  # noqa: BLE001
            out[well] = {"井名": cn, "路径": rel, "error": repr(exc)}
            print(f"[{tag}] {well}: ERROR {exc!r}")
    dest = HERE / ("旧值快照.json" if tag == "old" else "新值快照.json")
    dest.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[{tag}] written -> {dest}")


if __name__ == "__main__":
    main()
