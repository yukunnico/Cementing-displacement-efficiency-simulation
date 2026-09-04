# -*- coding: utf-8 -*-
"""对比三口井（hu101/ht1_003/ht1_004）新旧两轮结果的最终环空占位分布。
新 = results/胶塞语义修复后终跑_2026-09-03（B1+管容链+胶塞语义后）
旧 = results/呼101尾管_1D2D耦合模型 等 runner 目录（09-02 中午，B1 修复前）
只读，不写任何文件。
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(r"D:\users\desktop\research\控压固井项目\cement model\results")
NEW = ROOT / "胶塞语义修复后终跑_2026-09-03"
OLD = {
    "hu101": ROOT / "呼101尾管_1D2D耦合模型" / "呼101尾管_1D2D耦合模型_深度剖面.csv",
    "ht1_003": ROOT / "呼1-003_1D2D耦合模型" / "呼1-003_1D2D耦合模型_深度剖面.csv",
    "ht1_004": ROOT / "呼1-004_1D2D耦合模型" / "呼1-004_1D2D耦合模型_深度剖面.csv",
}


def find_col(cols, *keys):
    for c in cols:
        if all(k in c for k in keys):
            return c
    return None


def summarize_profile(path, tag):
    df = pd.read_csv(path)
    cols = list(df.columns)
    print(f"  [{tag}] 行数={len(df)}  列={cols}")
    dcol = find_col(cols, "深度") or find_col(cols, "md") or cols[0]
    ccol = find_col(cols, "水泥", "浓度")
    mcol = find_col(cols, "钻井液", "浓度")
    scol = find_col(cols, "前置", "浓度") or find_col(cols, "隔离", "浓度")
    fcol = find_col(cols, "冲洗", "浓度")
    d = df[dcol].to_numpy(float)
    print(f"    深度列={dcol}: {d.min():.1f} ~ {d.max():.1f} m")
    for name, col in (("水泥", ccol), ("钻井液", mcol), ("前置液", scol), ("冲洗液", fcol)):
        if col is None:
            print(f"    {name}: 列缺失")
            continue
        v = df[col].to_numpy(float)
        n = len(v)
        print(
            f"    {name}({col}): 全段均值={v.mean():.4f}  "
            f"底1/3={v[: n // 3].mean():.3f}  中1/3={v[n // 3 : 2 * n // 3].mean():.3f}  "
            f"顶1/3={v[2 * n // 3 :].mean():.3f}  行>0.5={int((v > 0.5).sum())}"
        )
        band = np.where(v >= 0.1)[0]
        if band.size:
            print(
                f"        >=0.1 带位置: 深度 {d[band].min():.1f} ~ {d[band].max():.1f} m"
                f"（行数 {band.size}/{n}）"
            )
    if ccol is not None:
        v = df[ccol].to_numpy(float)
        hit = np.where(v >= 0.5)[0]
        if hit.size:
            print(f"    水泥顶界(最浅 c>=0.5): {d[hit].min():.1f} m（域顶={d.min() if d[0] < d[-1] else d.max():.1f}）")


def summarize_json(well):
    sm = json.loads((NEW / f"{well}_摘要.json").read_text(encoding="utf-8"))
    print("  [新·摘要] 最终结果:", json.dumps(sm.get("最终结果", {}), ensure_ascii=False))
    print("  [新·摘要] 低尾指标:", json.dumps(sm.get("低尾指标", {}), ensure_ascii=False))
    for k, v in sm.get("评价窗效率", {}).items():
        print(f"    窗口[{k}] ({v.get('window_type')}): eta_E={v.get('eta_E')}  eta_N={v.get('eta_N')}")
    t0 = sm.get("tier0_diagnostics", {})
    mr = dict(t0.get("muskat_regime", {}))
    mr.pop("c_bar_grid", None)
    print("  [新·摘要] tier0.muskat:", json.dumps(mr, ensure_ascii=False))


for w in ("hu101", "ht1_003", "ht1_004"):
    print("=" * 28, w, "=" * 28)
    summarize_profile(NEW / f"{w}_深度剖面.csv", "新·终跑09-03")
    summarize_json(w)
    if OLD[w].exists():
        summarize_profile(OLD[w], "旧·runner目录09-02午")
    else:
        print("  [旧] 文件不存在:", OLD[w])
