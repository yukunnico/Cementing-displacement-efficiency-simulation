# -*- coding: utf-8 -*-
"""T2 时程观察层出图（2026-10-01）——每井时程多子图 + 等库存比截面对比小图。

数据源：``results/敏感性变体_温压T2时程_2026-10-01/``（批脚本落盘 CSV，只读）。

图
--
- ``{井}_时程观察.png``：四行（η_N / 饥饿份额 / cement_occ / 窄边前缘）
  × 各变体曲线；竖线 = 基线 ``Toff_zero`` 的 4 个事件时刻
  （隔离液出鞋 / 尾浆出鞋 / 泵注50% / 碰压前）。
- ``{井}_等库存比截面.png``：左 η_N、右 饥饿份额；横轴 = 库存比 0.5/0.8/1.0；
  分组柱 = 各变体（T-off 灰阶、T-on 彩色）；未达档标「未达到」不画柱。

用法
----
    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 conda run -n cementT \
        python scripts/plots/plot_sensitivity_temperature_history_20261001.py
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = PROJECT_ROOT / "results" / "敏感性变体_温压T2时程_2026-10-01"

# 中文回退：matplotlib 3.10 在 family='sans-serif' 下不沿 font.sans-serif 逐字回退，
# 把 font.family 设成同一列表，保证中文回落到 SimHei/微软雅黑（避免豆腐块）
plt.rcParams["font.sans-serif"] = [
    "Times New Roman", "SimHei", "Microsoft YaHei", "DejaVu Sans",
]
plt.rcParams["font.family"] = list(plt.rcParams["font.sans-serif"])
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 150
plt.rcParams["savefig.dpi"] = 200

# 变体 → (显示名, 线色, 线型)；T-off 灰阶、T-on 彩色
VARIANT_STYLE: dict[str, tuple[str, str, str]] = {
    "Toff_zero": ("T-off 基线", "#444444", "-"),
    "Ton_const60_rate_x1.0": ("T-on 60°C r1.0", "#ff7f0e", "-"),
    "Ton_static_rate_x1.0": ("T-on 静温 r1.0", "#1f77b4", "-"),
    "Ton_static_rate_x1.4": ("T-on 静温 r1.4", "#2ca02c", "--"),
    "Ton_table_rate_x1.0": ("T-on 表格 r1.0", "#d62728", ":"),
}
WELLS = ("呼1-004", "呼101", "呼103")
EVENTS = ("隔离液出鞋", "尾浆出鞋", "泵注50%", "碰压前")
OCC_LEVELS = ("0.5", "0.8", "1.0")
TIME_SERIES = (
    ("η_N", "窄边效率 η_N"),
    ("饥饿份额", "饥饿份额"),
    ("cement_occ", "库存比 cement_occ(=η_E)"),
    ("front_narrow_m", "窄边前缘 / m"),
)


def _read(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def _variants_of(well: str) -> list[str]:
    """该井在盘的变体（按 VARIANT_STYLE 序，未列名的殿后）。"""
    found = {p.name[len(well) + 1:-len("_时程.csv")]
             for p in OUT_DIR.glob(f"{well}_*_时程.csv")}
    ordered = [v for v in VARIANT_STYLE if v in found]
    ordered += sorted(found - set(VARIANT_STYLE))
    return ordered


def plot_history(well: str, variants: list[str]) -> Path:
    fig, axes = plt.subplots(len(TIME_SERIES), 1, figsize=(9.2, 10.5),
                             sharex=True)
    t0 = _read(OUT_DIR / f"{well}_Toff_zero_事件时刻表.csv")
    t0_map = {r["事件"]: r for r in t0}

    for ax, (col, title) in zip(axes, TIME_SERIES):
        for v in variants:
            rows = _read(OUT_DIR / f"{well}_{v}_时程.csv")
            if not rows:
                continue
            t_min = np.array([float(r["t_s"]) / 60.0 for r in rows])
            y = np.array([float(r[col]) for r in rows])
            label, color, ls = VARIANT_STYLE.get(v, (v, None, "-"))
            ax.plot(t_min, y, color=color, ls=ls, lw=1.6, label=label,
                    marker="o" if len(rows) <= 60 else None, ms=2.6)
        ax.set_ylabel(title, fontsize=10)
        ax.grid(alpha=0.3, lw=0.5)
        ax.set_ylim(bottom=0)
        # 4 条事件竖线（基线 Toff_zero 的时刻）；标签只标在首图内、
        # 奇偶交错高度——「隔离液出鞋≈泵注50%」两线只差几分钟，同高必叠
        for i, ev in enumerate(EVENTS):
            r = t0_map.get(ev)
            if not r or r.get("t_s", "") in ("", None):
                continue
            x = float(r["t_s"]) / 60.0
            ax.axvline(x, color="#888888", ls=":", lw=1.0, alpha=0.85)
            if ax is axes[0]:
                x0, x1 = ax.get_xlim()
                if x >= x0 + 0.85 * (x1 - x0):
                    ha = "right"
                elif x <= x0 + 0.15 * (x1 - x0):
                    ha = "left"
                else:
                    ha = "center"
                ax.annotate(
                    ev, xy=(x, 0.97 - i * 0.085),
                    xycoords=("data", "axes fraction"), ha=ha, va="top",
                    fontsize=7.5, color="#555555",
                )
    axes[0].legend(fontsize=9, loc="lower right", framealpha=0.9)
    axes[0].set_title(f"{well} 过程观察时程（竖线=基线 4 事件时刻）",
                      fontsize=12, fontweight="bold")
    axes[-1].set_xlabel("地面累计时间 / min", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    out = OUT_DIR / f"{well}_时程观察.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def plot_occ_sections(well: str, variants: list[str]) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6))
    metrics = (("η_N", "窄边效率 η_N"), ("饥饿份额", "饥饿份额"))
    x = np.arange(len(OCC_LEVELS), dtype=float)
    width = 0.8 / max(len(variants), 1)
    missing: list[str] = []

    for ax, (col, title) in zip(axes, metrics):
        for i, v in enumerate(variants):
            rows = {r["库存比档"]: r
                    for r in _read(OUT_DIR / f"{well}_{v}_等库存比截面.csv")}
            label, color, _ls = VARIANT_STYLE.get(v, (v, None, "-"))
            vals, reached = [], []
            for lv in OCC_LEVELS:
                r = rows.get(lv)
                if r is None:
                    vals.append(0.0)
                    reached.append(False)
                    continue
                ok = not str(r.get("状态", "")).startswith("未达到")
                vals.append(float(r[col]) if ok else 0.0)
                reached.append(ok)
                if not ok and col == metrics[0][0]:
                    missing.append(f"{v}@occ{lv}")
            left = x + (i - (len(variants) - 1) / 2.0) * width
            ax.bar(left, vals, width=width * 0.95, color=color,
                   edgecolor="white", lw=0.4, label=label)
            for xi, (val, ok) in enumerate(zip(vals, reached)):
                if not ok:
                    ax.annotate("未达到", xy=(left[xi], 0.02), ha="center",
                                va="bottom", fontsize=7.5, rotation=90,
                                color="#aa0000")
                elif col == "η_N":
                    ax.annotate(f"{val:.3f}", xy=(left[xi], val), ha="center",
                                va="bottom", fontsize=7, color="#333333")
        ax.set_xticks(x)
        ax.set_xticklabels([f"occ={lv}" for lv in OCC_LEVELS])
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_ylabel("指标值", fontsize=10)
        ax.grid(alpha=0.3, axis="y", lw=0.5)
        ax.set_ylim(0, 1.12)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=min(len(labels), 5),
               fontsize=9, frameon=False, bbox_to_anchor=(0.5, 0.945))
    fig.suptitle(f"{well} 等库存比截面对比（T-on vs T-off）",
                 fontsize=12, fontweight="bold", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.87))
    out = OUT_DIR / f"{well}_等库存比截面.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def main() -> int:
    if not OUT_DIR.is_dir():
        raise SystemExit(f"缺输出目录：{OUT_DIR}（先跑批脚本）")
    made: list[Path] = []
    for well in WELLS:
        variants = _variants_of(well)
        if not variants:
            print(f"[跳过] {well}: 无时程 CSV")
            continue
        print(f"{well}: {len(variants)} 变体 → {', '.join(variants)}")
        made.append(plot_history(well, variants))
        made.append(plot_occ_sections(well, variants))
    for p in made:
        print(f"  [图] {p.relative_to(PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
