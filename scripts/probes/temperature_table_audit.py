"""
T0-3 温度表只读审计 probe（scripts/probes/temperature_table_audit.py）

用 ``TableTemperatureField`` / ``load_delivered_pair`` 读交付的呼1-004 管内/环空
温度表（T_in.xlsx / T_out.xlsx，333 深度 × 200 分钟），出三联图并打印审计摘要：

- 子图1：T_in 热图（x=时间 min，y=深度 m，y 轴翻转朝下）
- 子图2：T_out 热图
- 子图3：初始时刻（col0）两表剖面 + 地温静温线 T(z)=16.006+1.7598e-2·z 对照
  + 本井电测静温锚点（temperature_pressure_profile.csv，可选读入）
- 图注：表值比电测静温偏冷 4–16°C（工况口径差异）

锚点口径：地温线即温度表 col0（最大残差 ~0.05°C），"偏冷 4–16°C" 是**表 vs 电测
静温**的差（5241 m：124−108.2≈15.8°C；7660 m：155−150.8≈4.2°C），不是表 vs 地温线。

只读约束：不改任何输入数据与库代码，唯一落盘产物是审计图
``results/温度表审计/T表审计.png``（加载用 ``use_cache=False``，不写 npz 缓存）。

用法（cementT 环境，工作目录=worktree 根）::

    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/probes/temperature_table_audit.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from cemdisp.data.temperature_field import (  # noqa: E402  (需先插 sys.path)
    load_delivered_pair,
)

# matplotlib 3.10 在 family='sans-serif' 下不沿 font.sans-serif 逐字回退，
# 把 font.family 设为同一列表，保证中文回落到微软雅黑/黑体（避免豆腐块）
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["font.family"] = list(plt.rcParams["font.sans-serif"])
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 11

OUT_DIR = _PROJECT_ROOT / "results" / "温度表审计"
OUT_PNG = OUT_DIR / "T表审计.png"

# 地温静温拟合线（温压耦合改进计划_2026-09-30 §1 裁定；对温度表 col0 拟合残差 ~0.04°C）
GEO_A = 16.006          # °C，地表截距
GEO_B = 1.7598e-2       # °C/m，地温梯度

# 本井电测静温锚点（可选读入；"偏冷 4–16°C" 口径即表 vs 这些锚点）
ANCHOR_CSV = (
    _PROJECT_ROOT / "参考文档" / "现场资料提取" / "ht1_004_呼1-004"
    / "temperature_pressure_profile.csv"
)

# 节点命中抽查索引（(行, 列) = (深度点, 分钟)）
SPOT_CHECKS = ((0, 0), (166, 100), (332, 199))


def _fmt_range(x: np.ndarray) -> str:
    return f"[{float(np.min(x)):.2f}, {float(np.max(x)):.2f}] °C"


def build_figure(t_in, t_out) -> tuple[plt.Figure, dict]:
    """画三联审计图；返回 figure 与计算出的对照量（供摘要打印）。"""
    z = t_in.depth_m                       # (333,) m，30 → 7660
    t_min = t_in.time_s / 60.0             # (200,) min，0 → 199
    z_top, z_bot = float(z[0]), float(z[-1])
    t_left, t_right = float(t_min[0]), float(t_min[-1])

    fig = plt.figure(figsize=(15.5, 9.0))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.15, 1.0], hspace=0.28, wspace=0.18)
    ax_in = fig.add_subplot(gs[0, 0])
    ax_out = fig.add_subplot(gs[0, 1])
    ax_prof = fig.add_subplot(gs[1, :])

    # 两张热图共用色标范围，便于管内/环空横向比较
    vmin = float(min(t_in.table.min(), t_out.table.min()))
    vmax = float(max(t_in.table.max(), t_out.table.max()))

    for ax, tab, title in (
        (ax_in, t_in.table, "管内 T_in(z, t) 热图"),
        (ax_out, t_out.table, "环空 T_out(z, t) 热图"),
    ):
        im = ax.imshow(
            tab,
            aspect="auto",
            origin="upper",                    # 行0=最浅 → 深度向下增大
            extent=[t_left, t_right, z_bot, z_top],  # (左,右,下,上)：上=浅、下=深
            cmap="inferno",
            vmin=vmin,
            vmax=vmax,
            interpolation="nearest",
        )
        ax.set_title(title, fontsize=12)
        ax.set_xlabel("时间 (min)")
        fig.colorbar(im, ax=ax, pad=0.015, label="温度 (°C)")
    ax_in.set_ylabel("深度 (m)")
    ax_out.set_ylabel("深度 (m)")
    ax_in.grid(False)
    ax_out.grid(False)

    # ---- 子图3：初始时刻（col0）两表剖面 + 地温静温线 + 电测静温锚点 ----
    T_geo = GEO_A + GEO_B * z
    T0_in = t_in.table[:, 0]
    T0_out = t_out.table[:, 0]
    geo_resid = float(np.max(np.abs(T_geo - T0_in)))   # 地温线 vs 表首列

    anchors = _load_electric_log_anchors(z)
    # 锚点处表值（col0 沿深度线性插值）与偏冷量
    a_md = anchors["md_m"].to_numpy(dtype=float)
    a_T = anchors["temperature_c"].to_numpy(dtype=float)
    a_T_table = np.interp(a_md, z, T0_in)
    a_delta = a_T - a_T_table                          # 正 = 表比电测偏冷

    ax_prof.plot(T_geo, z, color="k", ls="-.", lw=2.0,
                 label=r"地温静温线 $T(z)=16.006+1.7598\times10^{-2}\,z$（≈col0）")
    ax_prof.plot(T0_in, z, color="#1f77b4", lw=1.8, label="T_in 初始时刻 (col0)")
    ax_prof.plot(T0_out, z, color="#d62728", lw=1.8, ls="--",
                 label="T_out 初始时刻 (col0)（与 T_in 重合）")
    if len(a_md):
        ax_prof.scatter(
            a_T, a_md, s=70, marker="D", c="#2ca02c", zorder=5,
            edgecolors="k", linewidths=0.6,
            label="电测静温锚点（现场资料）",
        )
        for k, (md, Tft, dlt) in enumerate(zip(a_md, a_T, a_delta)):
            deep = k == len(a_md) - 1       # 最深锚点：文字放点右侧空白，避让剖面线
            ax_prof.annotate(
                f"{md:.0f} m：{Tft:.0f} °C（偏冷 {dlt:.1f} °C）",
                xy=(Tft, md),
                xytext=(12, 4) if deep else (12, -12),
                textcoords="offset points",
                ha="left",
                va="bottom" if deep else "top",
                fontsize=9,
                color="#1a6b1a",
            )
    # 右侧留白（容纳最深锚点标注），上下留白（锚点标记不贴边）
    ax_prof.set_xlim(
        float(min(T0_in.min(), T_geo.min())) - 3.0,
        float(max(T0_in.max(), T_geo.max(), a_T.max() if len(a_md) else 0.0)) + 42.0,
    )
    ax_prof.set_ylim(z_bot + 90.0, z_top - 60.0)   # y 轴翻转：浅在上
    ax_prof.set_xlabel("温度 (°C)")
    ax_prof.set_ylabel("深度 (m)")
    ax_prof.set_title("初始时刻 (col0) 两表剖面 vs 地温静温线 / 电测静温锚点（呼1-004）",
                      fontsize=12)
    ax_prof.grid(alpha=0.3)
    ax_prof.legend(loc="upper right", fontsize=10, framealpha=0.9)

    if len(a_md):
        d_txt = "，".join(
            f"{md:.0f} m {d:.1f}" for md, d in zip(a_md, a_delta)
        )
        offset_txt = f"实测偏冷量：{d_txt} °C（≈4–16 °C）"
    else:
        offset_txt = "锚点文件缺失，未绘制电测静温对照"
    ax_prof.text(
        0.01,
        0.02,
        "表值比电测静温偏冷 4–16 °C（工况口径差异）\n" + offset_txt,
        transform=ax_prof.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.4", fc="#fff6d5", ec="#c8a400", alpha=0.95),
    )

    fig.suptitle(
        "T0-3 温度表审计 · 呼1-004（HT1-004_T.m 交付件：333 深度 × 200 min）",
        fontsize=14,
        y=0.985,
    )
    fig.text(
        0.5,
        0.015,
        "注：交付温度表比电测静温偏冷 4–16 °C（工况口径差异，结果标注即可，不做锚点校正）；"
        "深度轴取呼1-004井身结构表第 2 列，时间轴每列 1 min，col0=初始时刻。",
        ha="center",
        fontsize=10,
        color="#444444",
    )

    stats = {
        "geo_resid": geo_resid,
        "anchors_md": a_md,
        "anchors_T": a_T,
        "anchors_T_table": a_T_table,
        "anchors_delta": a_delta,
        "T0_in": T0_in,
        "T0_out": T0_out,
        "T_geo": T_geo,
        "t_min": t_min,
    }
    return fig, stats


def _load_electric_log_anchors(z: np.ndarray) -> pd.DataFrame:
    """读本井电测静温锚点（field_measured + notes 含"电测" + 落在深度域内）。

    锚点文件缺失/列不符时返回空表（图与摘要退化为无锚点对照，不报错）。
    """
    if not ANCHOR_CSV.exists():
        print(f"[警告] 锚点 CSV 缺失，跳过电测静温对照: {ANCHOR_CSV}")
        return pd.DataFrame(columns=["md_m", "temperature_c"])
    try:
        df = pd.read_csv(ANCHOR_CSV)
        mask = (
            df["temperature_c"].notna()
            & (df["data_type"] == "field_measured")
            & df["notes"].astype(str).str.contains("电测", na=False)
            & (df["md_m"] >= float(z[0]))
            & (df["md_m"] <= float(z[-1]))
        )
        return df.loc[mask, ["md_m", "temperature_c"]].reset_index(drop=True)
    except Exception as exc:  # noqa: BLE001 —— 只读 probe 不因锚点文件失败而中断
        print(f"[警告] 锚点 CSV 读取失败，跳过电测静温对照: {exc}")
        return pd.DataFrame(columns=["md_m", "temperature_c"])


def _anchor_line(stats: dict) -> str:
    """电测静温锚点偏冷量摘要行（锚点缺失时给退化文案）。"""
    md = stats["anchors_md"]
    if len(md) == 0:
        return "电测锚点对照  : 锚点文件缺失/无匹配行（见上方 [警告]）"
    parts = [
        f"{m:.0f} m 实测 {Tf:.0f} vs 表 {Tb:.1f} → 偏冷 {d:.1f}"
        for m, Tf, Tb, d in zip(
            md,
            stats["anchors_T"],
            stats["anchors_T_table"],
            stats["anchors_delta"],
        )
    ]
    dmin = float(np.min(stats["anchors_delta"]))
    dmax = float(np.max(stats["anchors_delta"]))
    ok = "✓" if (3.5 <= dmin <= 4.5 and 15.0 <= dmax <= 16.5) else "!"
    return (
        f"电测锚点对照  : {' ; '.join(parts)} °C → 范围 "
        f"{dmin:.1f}~{dmax:.1f} °C（口径声明 4–16 °C {ok}）"
    )


def audit_summary(t_in, t_out, stats) -> None:
    """打印只读审计摘要（形状/域/井底两表差/clamp 计数等）。"""
    shape_in, shape_out = t_in.table.shape, t_out.table.shape
    z = t_in.depth_m
    t_s = t_in.time_s
    dz = np.diff(z)
    bottom_diff = float(np.max(np.abs(t_in.table[-1] - t_out.table[-1])))

    lines = [
        "=" * 78,
        "T0-3 温度表审计摘要（只读 probe）— 呼1-004 / T_in.xlsx + T_out.xlsx",
        "=" * 78,
        f"表形状        : T_in {shape_in} , T_out {shape_out}（深度×时间；"
        f"期望 {(333, 200)} {'✓' if shape_in == shape_out == (333, 200) else '✗'}）",
        f"深度域        : {z[0]:.1f} → {z[-1]:.1f} m（{len(z)} 点，步长 "
        f"{dz.min():.2f}~{dz.max():.2f} m，中位 {float(np.median(dz)):.2f} m）",
        f"时间域        : {t_s[0] / 60:.0f} → {t_s[-1] / 60:.0f} min"
        f"（= {t_s[0]:.0f} → {t_s[-1]:.0f} s；每列 1 min，col0=初始时刻）",
        f"全表温度范围  : T_in {_fmt_range(t_in.table)} , T_out {_fmt_range(t_out.table)}",
        f"col0 两表一致性: max|T_in−T_out| 首列 = "
        f"{float(np.max(np.abs(t_in.table[:, 0] - t_out.table[:, 0]))):.3e} °C"
        f"{'（逐位相同 ✓）' if np.array_equal(t_in.table[:, 0], t_out.table[:, 0]) else ''}",
        f"井底两表差    : max|T_in−T_out| 最深行 = {bottom_diff:.3e} °C "
        f"{'✓ 共享同值' if bottom_diff <= 1e-9 else '✗ 不一致'}",
        f"地温线 vs col0: T(z)=16.006+1.7598e-2·z 与表首列 max|残差| = "
        f"{stats['geo_resid']:.3f} °C（表首列即线性地温）",
        _anchor_line(stats),
        f"clamp 计数    : T_in oob_count={t_in.oob_count} , "
        f"T_out oob_count={t_out.oob_count}（probe 仅域内查询）",
    ]

    # 节点命中抽查：T(z_i, t_j) 须逐位等于表内值
    hit = 0
    for i, j in SPOT_CHECKS:
        got = t_in.T(float(t_in.depth_m[i]), float(t_in.time_s[j]))
        ok = got == t_in.table[i, j]
        hit += int(ok)
    lines.append(
        f"节点命中抽查  : {hit}/{len(SPOT_CHECKS)} 与表内值逐位相等"
        f"（(行,列)={SPOT_CHECKS}）"
    )

    # clamp 自检：故意越界一次，验证审计计数与落点，随后复位
    n_before = t_in.oob_count
    t_in.T(float(z[0] - 100.0), float(t_s[-1] + 600.0))
    ev = t_in.oob_events[-1]
    lines.append(
        f"clamp 自检    : 越界 ({ev.md_m:.0f} m, {ev.t_s / 60:.0f} min) → clamp 至 "
        f"({ev.md_clamped_m:.1f} m, {ev.t_clamped_s / 60:.0f} min)；"
        f"T_in oob_count {n_before}→{t_in.oob_count} ✓"
    )
    t_in.reset_audit()
    t_out.reset_audit()
    lines += [
        f"clamp 计数(复位后): T_in={t_in.oob_count} , T_out={t_out.oob_count}",
        f"图件          : {OUT_PNG.relative_to(_PROJECT_ROOT)}",
        "=" * 78,
    ]
    print("\n".join(lines))


def main() -> int:
    # use_cache=False：不写 npz（只读约束），每次直接解析交付 xlsx + 井身结构 CSV
    t_in, t_out = load_delivered_pair(use_cache=False)

    fig, stats = build_figure(t_in, t_out)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
    plt.close(fig)

    audit_summary(t_in, t_out, stats)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
