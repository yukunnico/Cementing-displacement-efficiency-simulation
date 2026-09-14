"""呼101 尾管井：模型计算顶替情况 vs 现场固井质量 —— 逐段对照图

修复记录（2026-09-11）：
  旧版依赖 `参考文档/hu101_1d2d_paper_project/…/呼101尾管_各段固井质量与模型顶替效率.csv`，
  该目录已整个不存在（孤儿脚本，产出的 PNG 停在 08-01）。且旧版把"模型顶替效率_%"当外部输入，
  模型口径一变（如 2026-09-11 居中度改 0.80）该列即为过期值。
  本版：
    - 现场侧改用 `cbl_evaluation.csv`（人工整理的定性分段表，带溯源与置信度）
    - 模型侧改为**从当前深度剖面 CSV 实时计算**逐段均值
    - 新增逐段对照表 CSV
    - 图内标注模型口径（e 域均值），防止把临时反推情景误当验证口径

⚠️ 口径声明：现场侧是**声幅类胶结定性评价**，模型侧是**体积浓度分数**，两者不是同一物理量。
   本图是**并置对照**，不是一一对应的验证图。现场定量合格率只有整段一个数：62.77%。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[2]

MODEL_DIR = PROJECT_ROOT / "results" / "呼101尾管_1D2D耦合模型"
MODEL_PROFILE = MODEL_DIR / "呼101尾管_1D2D耦合模型_深度剖面.csv"
MODEL_SUMMARY = MODEL_DIR / "呼101尾管_1D2D耦合模型_结果摘要.json"
FIELD_SEGMENTS = PROJECT_ROOT / "参考文档" / "现场资料提取" / "hu101_呼101" / "cbl_evaluation.csv"

OUTPUT_PNG = MODEL_DIR / "呼101_固井质量与模型效率沿深度分布图.png"
OUTPUT_CSV = MODEL_DIR / "呼101_逐段模型vs现场对照表.csv"

DEPTH_TOP, DEPTH_BOTTOM = 5380.0, 7870.0   # 取 5380 起：现场分段自 5390m 起，避免首段标签被轴裁掉
WHOLE_WELL_LENGTH_M = 1000.0          # 超过此长度的行视为整井汇总，不作为逐段
FIELD_OFFICIAL_PASS_RATE = 62.77      # 官方 CBL 合格率（5390–7810m，100312.PDF 图头）

CBL_QUALITY_COLOR = {
    "良好": "#2E7D32",
    "良好至中等": "#9CCC65",
    "中等": "#F9A825",
    "差至中等": "#EF6C00",
    "差": "#C62828",
    "未分类": "#9E9E9E",
}
MODEL_TIER_COLOR = {"高效率(≥80%)": "#4CAF50", "中等效率(50–80%)": "#DAA520", "低效率(<50%)": "#D32F2F"}


def _pick_cjk_font() -> str | None:
    """从 matplotlib 实际可用字体里挑一个中文字体（旧的 try/except 对缺字体无效）。"""
    available = {f.name for f in fm.fontManager.ttflist}
    for cand in ("Microsoft YaHei", "SimHei", "SimSun", "Arial Unicode MS",
                 "Noto Sans CJK SC", "WenQuanYi Micro Hei"):
        if cand in available:
            return cand
    return None


def grade_bucket(text: object) -> str:
    """把现场质量描述归到可配色的一档。"""
    t = str(text).strip()
    if "交替" in t or t.startswith("差至中等"):
        return "差至中等"
    if "至中等" in t or "至良好" in t:
        return "良好至中等"
    if t.startswith("不合格") or "差为主" in t or t.startswith("差"):
        return "差"
    if t.startswith("中等"):
        return "中等"
    if t.startswith("优良") or t.startswith("良好"):
        return "良好"
    return "未分类"


def load_field_segments(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"现场分段表不存在：{path}")
    df = pd.read_csv(path, encoding="utf-8-sig")
    df = df.dropna(subset=["cbl_top_md_m", "cbl_bottom_md_m"]).copy()
    df["段长_m"] = df["cbl_bottom_md_m"] - df["cbl_top_md_m"]
    df = df[(df["段长_m"] < WHOLE_WELL_LENGTH_M) & (df["cbl_bottom_md_m"] > DEPTH_TOP)]
    df["等级"] = df["cbl_quality_class"].map(grade_bucket)
    return df.sort_values("cbl_top_md_m").reset_index(drop=True)


def load_model_profile(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"模型深度剖面不存在：{path}")
    df = pd.read_csv(path, encoding="utf-8-sig")
    df = df[(df["井深_m"] >= DEPTH_TOP) & (df["井深_m"] <= DEPTH_BOTTOM)].copy()
    return df.sort_values("井深_m").reset_index(drop=True)


def model_caliber() -> dict[str, float]:
    """读模型摘要里的口径信息（e 域均值），用于图内标注。"""
    if not MODEL_SUMMARY.exists():
        return {}
    data = json.loads(MODEL_SUMMARY.read_text(encoding="utf-8"))
    e = float(data.get("tier0_diagnostics", {}).get("muskat_regime", {})
              .get("eccentricity", float("nan")))
    return {
        "eta_E_well": float(data["最终结果"]["全井段最终有效顶替效率"]),
        "e_mean": e,
        "standoff_mean": 1.0 - e if np.isfinite(e) else float("nan"),
        "eta_E_cbl": float(data.get("评价窗效率", {})
                           .get("CBL评价井段(单层套管可评价段)", {}).get("eta_E", float("nan"))),
    }


def segment_stats(field_df: pd.DataFrame, model_df: pd.DataFrame) -> pd.DataFrame:
    """逐段：模型效率均值 vs 现场质量等级。"""
    rows = []
    for i, r in field_df.iterrows():
        top, bottom = float(r["cbl_top_md_m"]), float(r["cbl_bottom_md_m"])
        m = model_df[(model_df["井深_m"] >= top) & (model_df["井深_m"] < bottom)]
        rows.append({
            "段号": i + 1,
            "顶界_m": round(top, 1),
            "底界_m": round(bottom, 1),
            "段长_m": round(bottom - top, 1),
            "现场质量等级": r["等级"],
            "现场质量描述": r["cbl_quality_class"],
            "现场合格率_pct": r.get("cbl_pass_rate", np.nan),
            "现场数据来源": r.get("data_type", ""),
            "现场置信度": r.get("confidence", ""),
            "模型ηE_均值_pct": round(float(m["平均有效顶替效率"].mean()) * 100, 1) if len(m) else np.nan,
            "模型ηN_均值_pct": round(float(m["窄边有效效率"].mean()) * 100, 1) if len(m) else np.nan,
            "模型宽边_均值_pct": round(float(m["宽边有效效率"].mean()) * 100, 1) if len(m) else np.nan,
            "段内模型点数": len(m),
        })
    return pd.DataFrame(rows)


def plot_overview(field_df: pd.DataFrame, model_df: pd.DataFrame, seg_df: pd.DataFrame,
                  caliber: dict[str, float], output_path: Path) -> None:
    fig = plt.figure(figsize=(19, 12))
    gs = fig.add_gridspec(1, 4, width_ratios=[0.9, 1.3, 3.0, 2.6], wspace=0.08)

    def setup_depth_axis(ax, title: str) -> None:
        ax.set_ylim(DEPTH_TOP, DEPTH_BOTTOM)
        ax.invert_yaxis()
        ax.set_title(title, fontsize=11)
        ax.set_xticks([])
        ax.set_yticks([])
        for side in ("top", "right", "bottom"):
            ax.spines[side].set_visible(False)

    # ---------- 面板 1：深度段标签 ----------
    ax1 = fig.add_subplot(gs[0, 0])
    setup_depth_axis(ax1, "井深分段")
    for _, r in seg_df.iterrows():
        ax1.axhline(r["顶界_m"], color="black", lw=0.6)
        ax1.axhline(r["底界_m"], color="black", lw=0.6)
        ax1.text(0.5, (r["顶界_m"] + r["底界_m"]) / 2,
                 f"{r['顶界_m']:.0f}\n–\n{r['底界_m']:.0f}",
                 ha="center", va="center", fontsize=8)
    ax1.set_xlim(0, 1)

    # ---------- 面板 2：现场 CBL 定性质量色带 ----------
    ax2 = fig.add_subplot(gs[0, 1])
    setup_depth_axis(ax2, "现场 CBL 胶结质量\n（定性分段表 · 声幅类）")
    ax2.set_xlim(0, 1)
    for _, r in seg_df.iterrows():
        rect = Rectangle((0.05, r["顶界_m"]), 0.9, r["底界_m"] - r["顶界_m"],
                         facecolor=CBL_QUALITY_COLOR[r["现场质量等级"]],
                         edgecolor="black", lw=0.5)
        ax2.add_patch(rect)
        q = r["现场质量等级"]
        ax2.text(0.5, (r["顶界_m"] + r["底界_m"]) / 2, q, ha="center", va="center",
                 fontsize=9, fontweight="bold",
                 color="white" if q in ("差", "良好", "差至中等") else "black")

    # ---------- 面板 3：模型逐段 η_E 横向条形 ----------
    ax3 = fig.add_subplot(gs[0, 2])
    for _, r in seg_df.iterrows():
        eff = r["模型ηE_均值_pct"]
        mid = (r["顶界_m"] + r["底界_m"]) / 2
        if not np.isfinite(eff):
            continue
        tier = "高效率(≥80%)" if eff >= 80 else ("中等效率(50–80%)" if eff >= 50 else "低效率(<50%)")
        ax3.barh(mid, eff, height=r["底界_m"] - r["顶界_m"],
                 color=MODEL_TIER_COLOR[tier], edgecolor="black", lw=0.5, alpha=0.88)
        ax3.text(eff + 1.5, mid, f"{eff:.1f}%", ha="left", va="center", fontsize=9, fontweight="bold")
        ax3.text(2, mid, f"现场：{r['现场质量等级']}", ha="left", va="center",
                 fontsize=8.5, color="#37474F")
    ax3.set_ylim(DEPTH_TOP, DEPTH_BOTTOM)
    ax3.invert_yaxis()
    ax3.set_xlim(0, 112)
    ax3.set_xlabel("模型计算顶替效率 / %（段内均值）", fontsize=11)
    ax3.set_yticks([])
    ax3.grid(axis="x", linestyle="--", alpha=0.4)
    ax3.legend(handles=[Patch(facecolor=c, edgecolor="black", label=k)
                        for k, c in MODEL_TIER_COLOR.items()],
               loc="upper right", fontsize=9, title="模型顶替效率")

    # ---------- 面板 4：模型深度剖面 ----------
    ax4 = fig.add_subplot(gs[0, 3])
    depth = model_df["井深_m"].to_numpy()
    for col, color, label, lw in (
        ("平均有效顶替效率", "#1565C0", "周向平均", 2.0),
        ("宽边有效效率", "#2E7D32", "宽边", 1.4),
        ("中线有效效率", "#F9A825", "中线", 1.4),
        ("窄边有效效率", "#C62828", "窄边", 1.4),
    ):
        ax4.plot(model_df[col].to_numpy() * 100, depth, color=color, lw=lw, label=label, alpha=0.9)
    ax4.set_ylim(DEPTH_TOP, DEPTH_BOTTOM)
    ax4.invert_yaxis()
    ax4.set_xlim(0, 105)
    ax4.set_xlabel("模型沿深度有效率 / %", fontsize=11)
    ax4.set_yticks([])
    ax4.grid(axis="x", linestyle="--", alpha=0.4)
    ax4.legend(loc="lower left", fontsize=9)
    ax4_right = ax4.twinx()
    ax4_right.set_ylim(DEPTH_TOP, DEPTH_BOTTOM)
    ax4_right.invert_yaxis()
    ax4_right.set_ylabel("井深 / m", fontsize=11)

    # ---------- 标题与口径标注 ----------
    fig.suptitle("呼101 井尾管段：现场 CBL 胶结质量与模型计算顶替效率 沿深度逐段对照",
                 fontsize=16, fontweight="bold", y=0.985)
    caliber_txt = (f"模型口径：e 域均值 {caliber.get('e_mean', float('nan')):.3f}"
                   f"（居中度 {caliber.get('standoff_mean', float('nan')):.3f}）"
                   f"；η_E 全井 {caliber.get('eta_E_well', float('nan')):.4f}"
                   f"，CBL 窗 {caliber.get('eta_E_cbl', float('nan')):.4f}")
    fig.text(0.5, 0.945, caliber_txt, ha="center", fontsize=10.5, color="#37474F")
    fig.text(0.5, 0.918,
             "※ 口径声明：现场侧为声幅类胶结定性评价，模型侧为体积浓度分数——二者不是同一物理量，"
             "本图为并置对照，非一一对应验证。",
             ha="center", fontsize=10, style="italic", color="#B71C1C")
    fig.text(0.5, 0.012,
             f"现场分段源：参考文档/现场资料提取/hu101_呼101/cbl_evaluation.csv"
             f"（5390–5955m 为 Vision 逐段，5955–7810m 为像素推断，置信度低–中）；"
             f"官方整段 CBL 合格率 {FIELD_OFFICIAL_PASS_RATE:.2f}%（5390–7810m）。"
             f"模型源：results/呼101尾管_1D2D耦合模型/（口径随 runner 重跑自动更新）。",
             ha="center", fontsize=9, color="#546E7A")

    # 手动布局（不用 tight_layout：twinx + gridspec 组合会触发不兼容告警）
    fig.subplots_adjust(left=0.03, right=0.965, top=0.89, bottom=0.055, wspace=0.08)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[信息] 图表已保存：{output_path}")


def main() -> None:
    font = _pick_cjk_font()
    if font:
        plt.rcParams["font.family"] = font
        plt.rcParams["axes.unicode_minus"] = False
        print(f"[信息] 中文字体：{font}")
    else:
        print("[警告] 未找到中文字体，图内中文可能显示为方框。", file=sys.stderr)

    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    print("[信息] 读取现场定性分段表 …")
    field_df = load_field_segments(FIELD_SEGMENTS)
    print(f"       可用分段数：{len(field_df)}"
          f"（{field_df['cbl_top_md_m'].min():.0f}–{field_df['cbl_bottom_md_m'].max():.0f}m）")

    print("[信息] 读取模型深度剖面 …")
    model_df = load_model_profile(MODEL_PROFILE)
    print(f"       剖面点数：{len(model_df)}"
          f"（{model_df['井深_m'].min():.1f}–{model_df['井深_m'].max():.1f}m）")

    caliber = model_caliber()
    print(f"[信息] 模型口径：e={caliber.get('e_mean'):.3f}  η_E全井={caliber.get('eta_E_well'):.4f}")

    seg_df = segment_stats(field_df, model_df)
    seg_df.to_csv(OUTPUT_CSV, index=False, encoding="utf-8-sig")
    print(f"[信息] 逐段对照表已保存：{OUTPUT_CSV}")
    print(seg_df[["段号", "顶界_m", "底界_m", "现场质量等级", "模型ηE_均值_pct", "模型ηN_均值_pct"]]
          .to_string(index=False))

    plot_overview(field_df, model_df, seg_df, caliber, OUTPUT_PNG)


if __name__ == "__main__":
    main()
