# -*- coding: utf-8 -*-
"""魔鬼代言人独立核查（只读 results/，不改任何代码）。
L2: Spearman 重算 + 混杂因子（库存比/域长/入库完成度）偏相关
L6: hu102/hu101 深度剖面 偏心度指标 分布
L7: 8 井评价窗 eta_E/eta_N 全量扫描（找 CBL 窗与 1e-13 类值）
L9: 摘要.json 顶层字段清单（指标并列口径现状）
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(r"D:\users\desktop\research\控压固井项目\cement model")
RES = ROOT / "results" / "胶塞语义修复后终跑_2026-09-03"

# ---------- L2: Spearman + 偏相关 ----------
df = pd.read_csv(RES / "汇总.csv")
print("=" * 70)
print("[L2] 汇总列：", list(df.columns))
cols = ["本轮η_E", "η_N", "库存比", "域长m", "入环空水泥_m3", "域内水泥_m3",
        "物理环空体积_m3", "壁面冻结占比", "窄边到位率", "混浆指数"]
e = df["e域均"] if "e域均" in df.columns else None
# e 不在汇总里：从摘要 tier0 读
e_vals, well_order = [], list(df["井"])
for w in well_order:
    sm = json.loads((RES / f"{w}_摘要.json").read_text(encoding="utf-8"))
    e_vals.append(float(sm["tier0_diagnostics"]["muskat_regime"]["eccentricity"]))
df["e"] = e_vals
print(df[["井", "e"] + cols].to_string(index=False))

def sp(x, y):
    r, p = spearmanr(x, y)
    return round(r, 4), round(p, 4)

for label, sub in [("8井", df), ("7井(除hu101)", df[df["井"] != "hu101"]),
                   ("6井(e0.17-0.35)", df[~df["井"].isin(["hu101", "hu102"])])]:
    print(f"\n-- {label} (n={len(sub)}) --")
    for tgt in ["本轮η_E", "η_N", "窄边到位率", "混浆指数"]:
        print(f"  e vs {tgt}: rho,p = {sp(sub['e'], sub[tgt])}")

# 偏相关：控制混杂后 e vs eta_N / eta_E（秩残差法）
def rank_residual(y, X):
    """y 对 X 秩回归后的残差（X 为 DataFrame 多列）。"""
    yr = pd.Series(y).rank()
    Xr = X.rank()
    Xd = np.column_stack([np.ones(len(Xr))] + [Xr[c].values for c in Xr.columns])
    beta, *_ = np.linalg.lstsq(Xd, yr.values, rcond=None)
    return yr.values - Xd @ beta

print("\n[L2-混杂] e/库存比/域长 互相的秩相关：")
for a, b in [("e", "库存比"), ("e", "域长m"), ("e", "入环空水泥_m3"),
             ("库存比", "域长m")]:
    print(f"  {a} vs {b}: {sp(df[a], df[b])}")

for ctrl in [["库存比"], ["域长m"], ["库存比", "域长m"]]:
    re_ = rank_residual(df["e"], df[ctrl])
    for tgt in ["本轮η_E", "η_N"]:
        rt = rank_residual(df[tgt], df[ctrl])
        r, p = spearmanr(re_, rt)
        print(f"  偏相关 e vs {tgt} | 控制{ctrl}: rho={r:.4f} p={p:.4f}")

# ---------- L6: 深度剖面偏心度 ----------
print("\n" + "=" * 70)
print("[L6] 深度剖面 偏心度指标 分布：")
for w in ["hu102", "hu101", "hu1"]:
    dp = pd.read_csv(RES / f"{w}_深度剖面.csv", encoding="utf-8-sig")
    ec = pd.to_numeric(dp["偏心度指标"], errors="coerce")
    print(f"  {w}: min={ec.min():.4f} max={ec.max():.4f} mean={ec.mean():.4f} "
          f"n={ec.notna().sum()}  >0.5占比={(ec > 0.5).mean():.3f} "
          f"==0.55占比={(np.isclose(ec, 0.55, atol=1e-6)).mean():.3f}")
    so = pd.to_numeric(dp.get("居中度"), errors="coerce")
    if so is not None:
        print(f"      居中度: min={so.min():.3f} max={so.max():.3f}")

# ---------- L7 + L9: 评价窗与摘要字段 ----------
print("\n" + "=" * 70)
print("[L7/L9] 8 井评价窗效率扫描：")
tiny_hits = []
for w in well_order:
    sm = json.loads((RES / f"{w}_摘要.json").read_text(encoding="utf-8"))
    wins = sm.get("评价窗效率", {})
    cbl = {k: v for k, v in wins.items() if "CBL" in k}
    others = {k: v for k, v in wins.items() if "CBL" not in k}
    print(f"\n-- {w} --  窗数={len(wins)}")
    for k, v in cbl.items():
        print(f"   [CBL] {k}: eta_E={v.get('eta_E')} eta_N={v.get('eta_N')} type={v.get('window_type')}")
    # 非 CBL 窗只列 eta_E 异常小的
    for k, v in others.items():
        ee = v.get("eta_E")
        if isinstance(ee, (int, float)) and 0.0 < abs(ee) < 1e-6:
            tiny_hits.append((w, k, ee))
        if isinstance(ee, (int, float)) and ee == 0.0:
            print(f"   [ZERO] {k}: eta_E=0 eta_N={v.get('eta_N')}")
    if not cbl:
        print("   (无 CBL 窗)")
print("\n[L7] (0,1e-6) 区间 eta_E 命中：", tiny_hits if tiny_hits else "无")

print("\n[L9] hu1 摘要顶层字段：")
sm1 = json.loads((RES / "hu1_摘要.json").read_text(encoding="utf-8"))
print("  ", list(sm1.keys()))
print("  最终结果：", list(sm1["最终结果"].keys()))
print("  低尾指标：", list(sm1["低尾指标"].keys()))
print("  顶层是否含 eta_narrow/channeling/mixing：",
      [k for k in ("eta_narrow", "channeling_index", "mixing_index") if k in sm1])
print("  顶层是否含 前缘到位率 类字段：",
      [k for k in sm1 if "前缘" in k or "front" in k or "到位" in k])
print("  汇总.csv 是否含到位率列：", [c for c in df.columns if "到位" in c])

# η_E 恒等式复核
df["eta_calc"] = df["域内水泥_m3"] / df["物理环空体积_m3"]
print("\n[复核] η_E vs 域内/环空 最大偏差：",
      float((df["本轮η_E"] - df["eta_calc"]).abs().max()))
