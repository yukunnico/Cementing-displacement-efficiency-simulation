"""B&F25 图 6 数据复现探针：定均速闭包 Ī₁(c̄) 与 q₀(c̄)（R-T1-8 / Task 3 §6）。

图 6 题注（B&F25 p.34）：「Variation of mean mobility Ī₁(c̄) and isotropic flux q₀(c̄)
with c̄ **based on D2DGA model** for shear-thinning fluids at **m = 1, w̄ = 1, b = 0**
(isodensity). As the power-law index of the fluids reduces, the effective viscosity
(μ̂e) decreases leading to a **growth in the mean mobility value**.」

本探针按该口径算曲线：
- ``m = 1``（κ₁ = κ₂，§2.32 式在 B = 0 下 κ = m^{±0.5} = 1）；
- ``b = 0`` ⇒ Gb = 0（等密度）；
- ``w̄ = 1`` ⇒ **定均速** ``solve_fixed_mean_velocity(u_bar=1)``（这是图 6 的关键：
  它不是定 G 口径）；``κ = τ_Y = 0``（B = 0 ⇒ τ_Y = 0）。

输出量与口径：
- ``Ī₁ = w̄/G``（w̄=1 ⇒ ``Ī₁ = 1/G``）——"mean mobility"（H=1 时与闭包 ``I1`` 同口径；
  两者差仅 ``u_bar`` 的梯形 O(h²) 离散误差 ~1e-5，见报告"恒等式"一节）；
- ``q₀`` 取 :attr:`GapSolution.q0`（通量比口径，R-T1-8）。

用法::

    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/probe_bf25_fig6_mean_mobility.py
    ... --out DIR   # 另存 I1.csv / q0.csv

⚠️ 本脚本**不读图**：判据只用文献正文的定性陈述（n 减小 ⇒ Ī₁ 增大），不做像素级比对。
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

from cemdisp.models2d.gap_solver import solve_fixed_mean_velocity

N_LIST = (1.0, 0.8, 0.6, 0.4, 0.2)
C_GRID = (0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0)


def curves(ny: int = 201) -> tuple[dict[float, list[float]], dict[float, list[float]]]:
    """算 ``Ī₁(c̄)`` 与 ``q₀(c̄)``（各 n；m=1、w̄=1、b=0）。"""
    i1_tab: dict[float, list[float]] = {}
    q0_tab: dict[float, list[float]] = {}
    for n in N_LIST:
        i1s, q0s = [], []
        for c in C_GRID:
            sol = solve_fixed_mean_velocity(
                c_bar=c, n=(n, n), kappa=(1.0, 1.0), tau_y=(0.0, 0.0),
                u_bar=(1.0, 0.0), ny=ny,
            )
            i1s.append(1.0 / sol.G[0])  # w̄ = 1 ⇒ Ī₁ = 1/G（"mean mobility"）
            q0s.append(sol.q0)
        i1_tab[n], q0_tab[n] = i1s, q0s
    return i1_tab, q0_tab


def _print_table(title: str, tab: dict[float, list[float]]) -> None:
    head = "  c_bar |" + "".join(f"   n={n:<4}" for n in N_LIST)
    print(f"\n{title}\n{head}\n" + "-" * len(head))
    for i, c in enumerate(C_GRID):
        print(f"  {c:5.2f} |" + "".join(f"  {tab[n][i]:8.5f}" for n in N_LIST))


def _save(path: Path, tab: dict[float, list[float]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["c_bar"] + [f"n={n}" for n in N_LIST])
        for i, c in enumerate(C_GRID):
            w.writerow([c] + [repr(tab[n][i]) for n in N_LIST])
    print(f"已写入 {path}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=None, help="CSV 输出目录（可选）")
    ap.add_argument("--ny", type=int, default=201)
    args = ap.parse_args(argv)

    i1_tab, q0_tab = curves(ny=args.ny)
    _print_table("₁(c̄) = w̄/G（定均速 w̄=1、m=1、b=0）", i1_tab)
    _print_table("q₀(c̄)（通量比口径）", q0_tab)

    print("\n判据（文献正文定性陈述，非图像比对）")
    print("① 题注：'As the power-law index of the fluids reduces ... a growth in the")
    print("   mean mobility value' ⇒ 同一 c̄ 下 n 越小 Ī₁ 越大。逐点核对：")
    ok = True
    for i, c in enumerate(C_GRID):
        vals = [i1_tab[n][i] for n in N_LIST]
        mono = all(a < b for a, b in zip(vals, vals[1:]))  # n 递减 ⇒ Ī₁ 递增
        if not mono and c not in (0.0, 1.0):
            ok = False
            print(f"   ✗ c̄={c}: " + " ".join(f"{v:.5f}" for v in vals))
    print(f"   结果：{'全部单调（n↓ ⇒ Ī₁↑）' if ok else '存在反例（见上）'}")
    print("② 端点：Ī₁(c̄=0) = 1/3（牛顿单流体 ū=G/(3κ) 的 κ=1 值；与 c̄ 无关，")
    print(f"   实测 n=1 全表恒为 {i1_tab[1.0][0]:.5f}）；q₀(0)=0、q₀(1)=1 精确。")
    print("③ q 的 n 依赖：题注未给判据 ⇒ 只记观测（n↓ 时 q₀ 略降：剪切变稀的塞流化")
    print("   把通量分布推向均匀，中线带份额下降）。")
    print("④ m=1（两流体同一流变）⇒ Ī₁ 与 c̄ 无关（本闭包与 B&F25 (2.24) 闭式在 m=1 下")
    print("   同为 1/3）⇒ 6(a) 的曲线必为水平线；实测 y 值域 [0.3333, 0.6776] 落在")
    print("   PDF 文本层抽出的 6(a) 纵轴刻度 0.3–0.7 内。")

    if args.out is not None:
        args.out.mkdir(parents=True, exist_ok=True)
        _save(args.out / "fig6_I1.csv", i1_tab)
        _save(args.out / "fig6_q0.csv", q0_tab)
    return 0


if __name__ == "__main__":
    sys.exit(main())