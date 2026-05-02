"""
呼101井尾管 D2DGA 改进驱动脚本
================================

本脚本不覆盖原始 `hu101_d2dga_model.py`，而是在其基础上做可审计的外层改进：
1. 将悬挂器坐挂失败从全局 α 解释项中前移到几何/居中度代理修正；
2. 对 6050-6210m CBL 强异常段、6100-6400m 局部异常段施加局部 standoff 退化；
3. 对比 sustained_tail、volume_limited、tail_then_mud 三种环空入口边界模式；
4. 将模拟深度剖面与 100312 CBL 剖面代理值进行同深度对比，输出剖面误差指标。

注意：CBL 合格率/剖面是固井质量代理标签，不等同于实测水力顶替效率。
"""

from __future__ import annotations

import importlib
import json
from pathlib import Path
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import hu101_d2dga_model as base


OUT_DIR = Path(__file__).resolve().parent / "改进模型输出"
CBL_PROFILE_PATH = Path("D:/users/desktop/research/控压固井项目/参考文档/呼101/提取数据/100312_CBL剖面_Excel版.csv")
BOUNDARY_MODES = ["sustained_tail", "volume_limited", "tail_then_mud"]
PENALTY_WEIGHTS_SUMMARY = (
    base.CHANNELING_PENALTY_WEIGHT,
    base.MIXING_PENALTY_WEIGHT,
    base.INSTABILITY_PENALTY_WEIGHT,
)


def load_cbl_profile() -> pd.DataFrame:
    """读取 100312 CBL 深度剖面，保留数值字段。"""
    columns = [
        "row_id",
        "depth_md_m",
        "cbl_amplitude_pct",
        "quality_proxy",
        "quality_proxy_pct",
        "quality_grade_cn",
        "segment_type_cn",
        "segment_note_cn",
        "is_target_interval_cn",
        "is_double_liner_excluded_cn",
    ]
    last_error: UnicodeDecodeError | None = None
    for encoding in ("utf-8-sig", "gbk", "gb18030"):
        try:
            df = pd.read_csv(CBL_PROFILE_PATH, skiprows=2, names=columns, encoding=encoding)
            break
        except UnicodeDecodeError as exc:
            last_error = exc
    else:
        raise RuntimeError(f"无法读取CBL剖面CSV编码: {CBL_PROFILE_PATH}") from last_error
    for col in ["depth_md_m", "cbl_amplitude_pct", "quality_proxy", "quality_proxy_pct"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df.dropna(subset=["depth_md_m", "quality_proxy"]).reset_index(drop=True)


def make_hanger_failure_standoff(original: Callable[[base.Array, base.Array, base.Array], base.Array]) -> Callable[[base.Array, base.Array, base.Array], base.Array]:
    """生成带悬挂器失败/局部 CBL 异常约束的 standoff 剖面函数。"""

    def improved_standoff(md: base.Array, hole: base.Array, od: base.Array) -> base.Array:
        standoff = original(md, hole, od).copy()

        # 悬挂器坐挂失败的可审计代理：从悬挂器以下全段引入温和偏心退化。
        standoff *= np.where(md >= 5700.0, 0.94, 1.0)

        # CBL 剖面强异常：6050-6210m 质量代理约 0.07，按局部严重偏心/窄边滞留处理。
        severe = (md >= 6050.0) & (md <= 6210.0)
        standoff[severe] = np.minimum(standoff[severe], 0.34)

        # 已识别 6100-6400m 局部异常段：6210m 后恢复较快，但仍保留中等退化。
        moderate = (md > 6210.0) & (md <= 6400.0)
        standoff[moderate] = np.minimum(standoff[moderate], 0.48)

        # 7400m 以下近最大井斜段追加轻微偏心退化，避免油气水敏感层过度理想化。
        lower = (md >= 7400.0) & (md <= 7735.0)
        standoff[lower] = np.minimum(standoff[lower], original(md, hole, od)[lower] * 0.92)

        return np.clip(standoff, 0.28, 0.82)

    return improved_standoff


def compute_profile_diagnostics(run_dir: Path, cbl_profile: pd.DataFrame) -> dict[str, float]:
    """计算模拟剖面与 CBL 质量代理剖面的同深度诊断指标。"""
    profile_path = run_dir / "呼101尾管_D2DGA_深度剖面.csv"
    profile = pd.read_csv(profile_path, encoding="utf-8-sig").sort_values("井深_m")

    simulated = np.interp(
        cbl_profile["depth_md_m"].to_numpy(dtype=float),
        profile["井深_m"].to_numpy(dtype=float),
        profile["平均有效顶替效率"].to_numpy(dtype=float),
    )
    observed = cbl_profile["quality_proxy"].to_numpy(dtype=float)
    diff = simulated - observed

    severe_mask = (cbl_profile["depth_md_m"] >= 6050.0) & (cbl_profile["depth_md_m"] <= 6210.0)
    anomaly_mask = (cbl_profile["depth_md_m"] >= 6100.0) & (cbl_profile["depth_md_m"] <= 6400.0)
    target_mask = (cbl_profile["depth_md_m"] >= 6153.0) & (cbl_profile["depth_md_m"] <= 7741.0)

    comparison = cbl_profile[["depth_md_m", "quality_proxy", "quality_proxy_pct"]].copy()
    comparison["simulated_effective_efficiency"] = simulated
    comparison["simulated_minus_cbl_proxy"] = diff
    comparison.to_csv(run_dir / "呼101_模拟剖面_vs_CBL代理剖面.csv", index=False, encoding="utf-8-sig")

    def mean_for(mask: pd.Series | np.ndarray, arr: np.ndarray) -> float:
        mask_arr = np.asarray(mask, dtype=bool)
        return float(np.mean(arr[mask_arr])) if np.any(mask_arr) else float("nan")

    return {
        "CBL剖面对比点数": float(len(cbl_profile)),
        "模拟_minus_CBL代理_平均差": float(np.mean(diff)),
        "模拟_minus_CBL代理_MAE": float(np.mean(np.abs(diff))),
        "模拟_minus_CBL代理_RMSE": float(np.sqrt(np.mean(diff**2))),
        "CBL代理_全剖面均值": float(np.mean(observed)),
        "模拟剖面_全剖面均值": float(np.mean(simulated)),
        "6050_6210m_CBL代理均值": mean_for(severe_mask, observed),
        "6050_6210m_模拟均值": mean_for(severe_mask, simulated),
        "6100_6400m_CBL代理均值": mean_for(anomaly_mask, observed),
        "6100_6400m_模拟均值": mean_for(anomaly_mask, simulated),
        "6153_7741m_CBL代理均值": mean_for(target_mask, observed),
        "6153_7741m_模拟均值": mean_for(target_mask, simulated),
    }


def recalibrate_residual_alpha(summary: dict[str, object]) -> dict[str, float]:
    """在显式几何退化后，重标定残余质量惩罚α，避免与原α双重计入悬挂器失败。"""
    final = summary["最终结果"]
    hydraulic_efficiency = float(final["CBL评价井段模拟有效顶替效率"])
    field_reference = float(final["资料CBL合格率_代理顶替效率"])
    channeling = float(final["最终窜槽指数"])
    mixing = float(final["最终混浆指数"])
    instability = float(final["最终失稳指数"])
    raw_penalty = (
        base.CHANNELING_PENALTY_WEIGHT * channeling
        + base.MIXING_PENALTY_WEIGHT * mixing
        + base.INSTABILITY_PENALTY_WEIGHT * instability
    )
    target_quality_factor = float(np.clip(field_reference / max(hydraulic_efficiency, 1e-12), 0.0, 1.0))
    residual_alpha = float(np.clip((1.0 - target_quality_factor) / max(raw_penalty, 1e-12), 0.0, 1.0))
    recalibrated_quality_response = hydraulic_efficiency * np.clip(1.0 - residual_alpha * raw_penalty, 0.0, 1.0)
    return {
        "原始α": float(base.QUALITY_PENALTY_SCALE),
        "残余重标定α": residual_alpha,
        "原始α质量响应效率": float(final["CBL评价井段质量响应效率"]),
        "残余α质量响应效率": float(recalibrated_quality_response),
        "残余α质量响应_minus_资料": float(recalibrated_quality_response - field_reference),
        "目标质量因子": target_quality_factor,
        "惩罚总和S": raw_penalty,
        "说明": "显式加入悬挂器失败与局部standoff退化后，原α=0.671会重复吸收同一异常；残余α按总体CBL合格率后验归一化，仅代表未被几何项解释的剩余质量惩罚。该数值由CBL目标反算得到，不是独立预测或现场验证。",
    }


def rewrite_requested_charts(
    run_dir: Path,
    geom: dict[str, base.Array],
    x: base.Array,
    wall: base.Array,
    metrics: pd.DataFrame,
) -> None:
    """按当前需求覆盖改进模型输出中的三张图。"""
    charts = run_dir / "图表"
    charts.mkdir(parents=True, exist_ok=True)

    depth_profiles = pd.read_csv(run_dir / "呼101尾管_D2DGA_深度剖面.csv", encoding="utf-8-sig")
    final = metrics.iloc[-1]

    plt.figure(figsize=(7, 5))
    names = ["顶替效率", "资料CBL合格率"]
    vals = [float(final["cbl_eval_interval_efficiency"]), base.FIELD_REFERENCE_EFFICIENCY]
    bars = plt.bar(names, vals, color=["#3B82F6", "#F97316"])
    for bar, val in zip(bars, vals):
        plt.text(bar.get_x() + bar.get_width() / 2, val + 0.02, f"{val:.3f}", ha="center")
    plt.ylim(0, 1.1)
    plt.ylabel("效率/合格率")
    plt.title("呼101尾管模拟顶替效率与资料CBL合格率对比")
    plt.tight_layout()
    plt.savefig(charts / "呼101尾管_模拟与资料效率对比.png", dpi=220)
    plt.close()

    plt.figure(figsize=(7, 6))
    plt.plot(depth_profiles["平均有效顶替效率"], depth_profiles["井深_m"], label="有效顶替效率")
    plt.axvspan(
        0,
        base.FIELD_REFERENCE_EFFICIENCY,
        color="#F97316",
        alpha=0.12,
        label=f"资料CBL合格率{base.FIELD_REFERENCE_EFFICIENCY * 100:.2f}%",
    )
    plt.axhspan(7492.0, 7735.0, color="#EF4444", alpha=0.10)
    plt.gca().invert_yaxis()
    plt.xlabel("效率")
    plt.ylabel("井深 / m")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(charts / "呼101尾管_深度效率剖面对比.png", dpi=220)
    plt.close()

    cement = x[base.TRACKED.index("lead")] + x[base.TRACKED.index("tail")]
    eff = cement * (1.0 - wall)
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    field_specs = [
        (cement, "水泥浓度", "viridis"),
        (eff, "有效顶替效率", "YlGnBu"),
    ]
    extent = (float(geom["md"].min()), float(geom["md"].max()), 1.0, 0.0)
    for ax, (field, title, cmap) in zip(axes, field_specs):
        image = ax.imshow(field, aspect="auto", extent=extent, vmin=0.0, vmax=1.0, cmap=cmap)
        ax.set_ylabel("方位：宽边0 → 窄边1")
        ax.set_title(title)
        fig.colorbar(image, ax=ax, label="比例")
    axes[-1].set_xlabel("井深 / m")
    fig.suptitle("呼101尾管最终二维方位-轴向场分布")
    fig.tight_layout()
    fig.savefig(charts / "呼101尾管_最终二维场_浓度效率残余泥浆.png", dpi=220)
    plt.close(fig)


def run_mode(mode: str, cbl_profile: pd.DataFrame, original_standoff: Callable[[base.Array, base.Array, base.Array], base.Array]) -> dict[str, object]:
    """运行单个入口边界模式并写入增强摘要。"""
    importlib.reload(base)
    base.ANNULUS_BOUNDARY_MODE = mode
    base.OUT_DIR = OUT_DIR / mode
    base.OUT_DIR.mkdir(parents=True, exist_ok=True)
    base.standoff_profile = make_hanger_failure_standoff(original_standoff)

    geom, x, wall, metrics = base.simulate()
    summary = base.save_outputs(geom, x, wall, metrics)
    rewrite_requested_charts(base.OUT_DIR, geom, x, wall, metrics)
    diagnostics = compute_profile_diagnostics(base.OUT_DIR, cbl_profile)

    summary["改进项"] = [
        "新增悬挂器坐挂失败代理：5700m以下standoff整体温和退化。",
        "新增6050-6210m CBL强异常段局部严重standoff退化。",
        "新增6100-6400m局部异常段中等standoff退化。",
        "保留三种入口边界模式对比，避免单一边界模式误判。",
        "新增模拟深度剖面与CBL质量代理剖面同深度误差诊断。",
    ]
    summary["CBL剖面对比诊断"] = diagnostics
    summary["残余质量惩罚重标定"] = recalibrate_residual_alpha(summary)
    summary["假设说明"] = [
        item.replace(
            "悬挂器坐挂失败(技术总结记载)未在几何中建模，其影响通过α标定系数吸收。",
            "悬挂器坐挂失败已在本改进脚本中通过standoff代理退化显式进入几何；剩余差异再由残余α后验归一化。",
        )
        for item in summary.get("假设说明", [])
    ]
    summary["重要限制"] = [
        "CBL剖面为胶结/固井质量代理标签，不是水力顶替效率实测值。",
        "局部standoff退化为基于悬挂器失败与CBL异常段的代理修正，仍需实测caliper/偏心数据约束。",
        "残余α是使用总体CBL合格率反算的后验归一化参数；残余α质量响应效率不应作为独立预测值表述。",
        "6050-6210m强异常段的CBL代理远低于模拟水动力效率，提示该段可能包含胶结、气窜或测井响应等D2DGA水动力模型外机制。",
        "窜槽指数与失稳指数并非完全独立，当前线性加权惩罚可能存在相关项重复计入，应作为经验质量响应模型使用。",
        "本脚本不改动原始模型，所有改进输出均写入改进模型输出目录。",
    ]

    enhanced_summary_path = base.OUT_DIR / "呼101尾管_D2DGA_改进结果摘要.json"
    enhanced_summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def write_overall_report(results: list[dict[str, object]]) -> None:
    """写入跨边界模式汇总报告和总表。"""
    rows = []
    for summary in results:
        final = summary["最终结果"]
        diag = summary["CBL剖面对比诊断"]
        rows.append(
            {
                "边界模式": summary["环空入口边界模式"],
                "全井段最终有效顶替效率": final["全井段最终有效顶替效率"],
                "CBL评价井段模拟有效顶替效率": final["CBL评价井段模拟有效顶替效率"],
                "CBL评价井段质量响应效率": final["CBL评价井段质量响应效率"],
                "残余α质量响应效率": summary["残余质量惩罚重标定"]["残余α质量响应效率"],
                "残余重标定α": summary["残余质量惩罚重标定"]["残余重标定α"],
                "资料CBL合格率代理值": final["资料CBL合格率_代理顶替效率"],
                "最终窜槽指数": final["最终窜槽指数"],
                "最终混浆指数": final["最终混浆指数"],
                "最终失稳指数": final["最终失稳指数"],
                "剖面MAE": diag["模拟_minus_CBL代理_MAE"],
                "6050_6210m_CBL代理均值": diag["6050_6210m_CBL代理均值"],
                "6050_6210m_模拟均值": diag["6050_6210m_模拟均值"],
                "6100_6400m_CBL代理均值": diag["6100_6400m_CBL代理均值"],
                "6100_6400m_模拟均值": diag["6100_6400m_模拟均值"],
            }
        )

    table = pd.DataFrame(rows)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    table_path = OUT_DIR / "呼101_改进模型_边界模式对比总表.csv"
    table.to_csv(table_path, index=False, encoding="utf-8-sig")

    primary = table.loc[table["边界模式"] == "sustained_tail"].iloc[0]
    report_lines = [
        "# 呼101尾管 D2DGA 模型改进报告",
        "",
        "## 1. 改进目标",
        "本次改进在不覆盖原始 `hu101_d2dga_model.py` 的前提下，将原先主要由全局 α 吸收的现场异常前移到可审计的几何代理项：悬挂器坐挂失败、6050-6210m CBL强异常段、6100-6400m局部异常段，并补充三种环空入口边界模式对比。",
        "",
        "## 2. 方法与约束",
        "- 悬挂器坐挂失败：5700m以下 standoff 乘以0.94，代表全段温和偏心退化。",
        "- 6050-6210m强异常段：standoff 上限压到0.34，代表局部严重偏心/窄边滞留代理。",
        "- 6100-6400m异常段：standoff 上限压到0.48，代表异常段后续中等退化。",
        "- 7400-7735m敏感层：追加轻微偏心退化，避免最大井斜附近过度理想化。",
        "- CBL剖面：只作为固井质量代理标签参与剖面对比，不等同于实测水力顶替效率。",
        "",
        "## 3. 核心结果（基准对比口径 sustained_tail）",
        f"- 全井段最终有效顶替效率：{primary['全井段最终有效顶替效率']:.4f}",
        f"- CBL评价井段模拟有效顶替效率：{primary['CBL评价井段模拟有效顶替效率']:.4f}",
        f"- CBL评价井段质量响应效率（沿用原α=0.671）：{primary['CBL评价井段质量响应效率']:.4f}",
        f"- CBL评价井段质量响应效率（显式几何退化后残余α重标定）：{primary['残余α质量响应效率']:.4f}",
        f"- 残余重标定α：{primary['残余重标定α']:.4f}",
        f"- 资料CBL合格率代理值：{primary['资料CBL合格率代理值']:.4f}",
        f"- 剖面MAE（模拟有效效率 vs CBL代理）：{primary['剖面MAE']:.4f}",
        f"- 6050-6210m：CBL代理均值 {primary['6050_6210m_CBL代理均值']:.4f}，模拟均值 {primary['6050_6210m_模拟均值']:.4f}",
        "",
        "## 4. 三种入口边界模式对比",
        table.to_markdown(index=False),
        "",
        "## 5. 结论",
        "改进后模型把现场已知异常从单一全局补偿项拆解为局部、可解释的几何退化和边界模式敏感性。水动力预测应优先看 CBL评价井段模拟有效顶替效率；残余α质量响应效率是按总体 CBL 合格率反算的后验归一化值，不是独立预测或现场验证。6050-6210m 强异常段仍存在明显机制缺口，CBL代理均值仅0.0738，而 sustained_tail 模式水动力模拟均值仍为0.6923，说明该段可能受胶结质量、气窜、测井响应或其他非水动力因素控制。下一步若补充实测井径、扶正器实际下入/坐挂状态或声波/超声胶结资料，可将本次 standoff 代理项替换为实测几何输入，并单独建立胶结质量响应模型。",
        "",
        "## 6. 解释边界与复核意见",
        "- CBL合格率/剖面只作为固井质量代理标签，不等同于实测水力顶替效率。",
        "- 残余α由 `α=(1-CBL_ref/η_hyd)/S` 反算；当 `η_hyd >= CBL_ref` 时，质量响应等于CBL参考值是公式构造结果。",
        "- 当前惩罚模型中失稳指数包含窜槽因素，和显式窜槽项存在相关性；线性叠加结果应作为经验诊断，不宜过度物理解释。",
        "- 三种入口边界模式给出显著不确定性范围：CBL井段水动力效率约0.5884-0.8398，应比单一标定值更能代表模型不确定性。",
    ]
    (OUT_DIR / "呼101尾管_D2DGA_模型改进报告.md").write_text("\n".join(report_lines), encoding="utf-8")


def main() -> None:
    cbl_profile = load_cbl_profile()
    original_standoff = base.standoff_profile
    results = [run_mode(mode, cbl_profile, original_standoff) for mode in BOUNDARY_MODES]
    write_overall_report(results)
    print(json.dumps({"输出目录": str(OUT_DIR), "边界模式": BOUNDARY_MODES}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
