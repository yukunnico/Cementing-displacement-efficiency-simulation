"""
呼102井尾管 D2DGA 三因素敏感性分析
=====================================

本脚本对呼102井尾管固井顶替模型进行三因素全因子敏感性分析，
评估排量、塑性黏度(PV)和屈服值(YP)对顶替效率的影响程度。

    分析流程：
    1. 将原始模型的幂律流变参数(n, k)等效转换为Bingham参数(PV, YP)
    2. 定义三因素三水平（3×3×3=27个场景）的全因子实验设计
    3. 对每个场景：将Bingham参数反算为幂律参数→修改模型流体物性→运行模拟
    4. 分析各因素对CBL评价井段水动力效率(cbl_eval_interval_efficiency)的主效应和重要度排序
    5. 从高性能场景中提取建议参数区间

因素水平设置：
- 排量：1.00, 1.30, 1.60 m³/min
- 塑性黏度：60, 80, 100 mPa·s
- 屈服值：15, 22, 30 Pa

输出文件：
1. 呼102尾管_D2DGA_三因素敏感性_场景结果.csv — 所有场景的模拟结果
2. 呼102尾管_D2DGA_三因素敏感性_结果摘要.json — 影响排序、最优场景和建议区间
3. 呼102尾管_D2DGA_三因素敏感性_结果摘要.md — 可读性更好的Markdown格式摘要

注意：本分析仅评估顶替效率响应，不包含泵压窗口、ECD上限与漏失风险约束。
"""

from __future__ import annotations

import json
from collections.abc import Sequence
import importlib
from itertools import product
from pathlib import Path
from typing import Callable

import numpy as np
from numpy.typing import NDArray
import pandas as pd

if __package__:
    from . import hu102_tail_d2dga_model as base_model  # 导入基础D2DGA模型
else:
    base_model = importlib.import_module("hu102_tail_d2dga_model")  # 导入基础D2DGA模型

BoundaryState = Callable[[float], tuple[NDArray[np.float64], float, str]]

# ========== 路径配置 ==========
OUT_DIR = Path(__file__).resolve().parent
SCENARIO_CSV = OUT_DIR / "呼102尾管_D2DGA_三因素敏感性_场景结果.csv"
SUMMARY_JSON = OUT_DIR / "呼102尾管_D2DGA_三因素敏感性_结果摘要.json"
SUMMARY_MD = OUT_DIR / "呼102尾管_D2DGA_三因素敏感性_结果摘要.md"

# ========== Fann流变仪等效剪切速率 ==========
# 600rpm和300rpm对应的等效剪切速率，用于幂律↔Bingham参数转换
GAMMA_600 = 1022.0  # 600 rpm 等效剪切速率, 1/s
GAMMA_300 = 511.0   # 300 rpm 等效剪切速率, 1/s

# ========== 三因素水平设置 ==========
RATE_LEVELS_M3_MIN = [1.00, 1.30, 1.60]   # 排量水平，单位 m³/min
PV_LEVELS_MPA_S = [60.0, 80.0, 100.0]     # 塑性黏度水平，单位 mPa·s
YP_LEVELS_PA = [15.0, 22.0, 30.0]          # 屈服值水平，单位 Pa

# ========== 目标指标 ==========
# 用于评估因素影响和选择最优场景的指标
TARGET_METRIC = "cbl_eval_interval_efficiency"  # CBL评价井段水动力有效顶替效率
FULL_WELL_TARGET_METRIC = "effective_efficiency"  # 全井段有效顶替效率


def powerlaw_to_bingham(k_pa_sn: float, n: float) -> tuple[float, float]:
    """
    将幂律参数在300/600rpm点等效转换为Bingham PV/YP。

    转换方法：在Fann六速流变仪的300rpm和600rpm两个剪切速率点，
    使幂律模型和Bingham模型的切应力相等，从而求解PV和YP。

    计算步骤：
    1. 由幂律模型计算300rpm和600rpm处的切应力：
       τ_600 = k × γ_600^n
       τ_300 = k × γ_300^n
    2. 由Bingham模型反算PV和YP：
       PV = (τ_600 - τ_300) / (γ_600 - γ_300)
       YP = τ_300 - PV × γ_300

    Args:
        k_pa_sn: 幂律稠度系数，单位 Pa·s^n
        n: 幂律指数，无量纲

    Returns:
        (pv_mpa_s, yp_pa):
        - pv_mpa_s: 塑性黏度，单位 mPa·s（= PV_Pa·s × 1000）
        - yp_pa: 屈服值，单位 Pa
    """
    tau_600 = k_pa_sn * (GAMMA_600 ** n)
    tau_300 = k_pa_sn * (GAMMA_300 ** n)
    pv_pa_s = (tau_600 - tau_300) / (GAMMA_600 - GAMMA_300)
    yp_pa = tau_300 - pv_pa_s * GAMMA_300
    return pv_pa_s * 1000.0, yp_pa


def bingham_to_powerlaw(pv_mpa_s: float, yp_pa: float) -> tuple[float, float]:
    """
    将Bingham PV/YP在300/600rpm点拟合为幂律参数n,k。

    反向转换：给定Bingham模型的PV和YP，
    在300rpm和600rpm两个点计算切应力，
    然后拟合幂律参数使两点的切应力匹配。

    计算步骤：
    1. 由Bingham模型计算切应力：
       τ_600 = YP + PV × γ_600
       τ_300 = YP + PV × γ_300
    2. 由两点切应力拟合幂律参数：
       n = ln(τ_600/τ_300) / ln(γ_600/γ_300)
       k = τ_600 / γ_600^n

    Args:
        pv_mpa_s: 塑性黏度，单位 mPa·s
        yp_pa: 屈服值，单位 Pa

    Returns:
        (n, k):
        - n: 幂律指数
        - k: 稠度系数，单位 Pa·s^n
    """
    pv_pa_s = pv_mpa_s / 1000.0
    tau_600 = max(yp_pa + pv_pa_s * GAMMA_600, 1e-6)
    tau_300 = max(yp_pa + pv_pa_s * GAMMA_300, 1e-6)
    n = float(np.log(tau_600 / tau_300) / np.log(GAMMA_600 / GAMMA_300))
    k = float(tau_600 / (GAMMA_600 ** n))
    return n, k


def make_boundary_state(rate_m3_min: float) -> BoundaryState:
    """
    创建指定排量下的注入边界时序函数。

    由于敏感性分析需要改变排量，而基础模型的boundary_state函数
    使用硬编码的排量值，因此需要生成一个新的边界函数来替代。

    生成的函数与base_model.boundary_state接口一致，
    但使用自定义的排量参数。

    Args:
        rate_m3_min: 注入排量，单位 m³/min

    Returns:
        边界时序函数，签名为 (t: float) -> (vec, q_m3s, stage)
    """

    def boundary_state_custom(t: float) -> tuple[NDArray[np.float64], float, str]:
        """自定义排量的注入边界时序函数。"""
        tail_slurry_volume_m3 = 35.0 / 2.10  # 尾浆体积
        displacement_volume_m3 = 74.0         # 替浆量
        schedule = [
            ("balance", 0.0, rate_m3_min),
            ("spacer", 0.0, rate_m3_min),
            ("lead", 0.0, rate_m3_min),
            ("tail", tail_slurry_volume_m3, rate_m3_min),
        ]
        t0 = 0.0
        vec = np.zeros(len(base_model.TRACKED))
        for name, volume, rate in schedule:
            duration = volume / rate * 60.0
            if t < t0 + duration - 1e-12:
                vec[base_model.TRACKED.index(name)] = 1.0
                return vec, rate / 60.0, "注入尾管水泥浆"
            t0 += duration

        push_duration_s = displacement_volume_m3 / rate_m3_min * 60.0
        if t < t0 + push_duration_s:
            vec[base_model.TRACKED.index("tail")] = 1.0
            return vec, rate_m3_min / 60.0, "替浆推进"

        vec[base_model.TRACKED.index("tail")] = 1.0
        return vec, 0.0, "停泵保持"

    return boundary_state_custom


def run_case(rate_m3_min: float, pv_mpa_s: float, yp_pa: float) -> dict[str, float | str]:
    """
    运行单个敏感性分析场景并返回最终关键指标。

    执行步骤：
    1. 将Bingham参数(PV, YP)转换为幂律参数(n, k)
    2. 临时替换基础模型的流体物性和边界条件
    3. 运行完整模拟
    4. 恢复原始模型参数（使用try/finally确保恢复）
    5. 返回关键指标

    注意：通过直接修改模块级变量(base_model.FLUIDS, base_model.boundary_state)
    来实现参数替换，这是为了复用基础模型的simulate()函数而采用的技巧。
    使用try/finally确保即使模拟出错也能恢复原始参数。

    Args:
        rate_m3_min: 排量，单位 m³/min
        pv_mpa_s: 塑性黏度，单位 mPa·s
        yp_pa: 屈服值，单位 Pa

    Returns:
        包含以下字段的指标字典：
        - rate_m3_min, pv_mpa_s, yp_pa: 输入参数
        - n_fit, k_fit: 拟合得到的幂律参数
        - effective_efficiency: 全井段有效顶替效率
        - cbl_eval_interval_efficiency: CBL评价井段水动力效率
        - cbl_quality_proxy: CBL质量响应效率（目标指标）
        - target_interval_efficiency: 油气水层段效率
        - channeling_index: 窜槽指数
        - mixing_index: 混浆指数
        - instability_index: 失稳指数
    """
    # 保存原始模型参数，用于恢复
    original_fluids = base_model.FLUIDS
    original_boundary_state = base_model.boundary_state

    try:
        # 步骤1：Bingham → 幂律参数转换
        n_fit, k_fit = bingham_to_powerlaw(pv_mpa_s=pv_mpa_s, yp_pa=yp_pa)

        # 步骤2：替换流体物性（仅修改领浆和尾浆的幂律参数，密度不变）
        fluids = dict(original_fluids)
        lead_density = original_fluids["lead"].density_gcc
        tail_density = original_fluids["tail"].density_gcc
        fluids["lead"] = base_model.Fluid("领浆", lead_density, "power_law", n=n_fit, k=k_fit)
        fluids["tail"] = base_model.Fluid("尾管水泥浆", tail_density, "power_law", n=n_fit, k=k_fit)

        # 步骤3：替换模块级变量
        setattr(base_model, "FLUIDS", fluids)
        setattr(base_model, "boundary_state", make_boundary_state(rate_m3_min))

        # 步骤4：运行模拟
        _, _, _, metrics = base_model.simulate()
        final = metrics.iloc[-1]

        # 步骤5：收集结果
        return {
            "rate_m3_min": rate_m3_min,
            "pv_mpa_s": pv_mpa_s,
            "yp_pa": yp_pa,
            "n_fit": n_fit,
            "k_fit": k_fit,
            "effective_efficiency": float(final["effective_efficiency"]),
            "cbl_eval_interval_efficiency": float(final["cbl_eval_interval_efficiency"]),
            "cbl_quality_proxy": float(final["cbl_quality_proxy"]),
            "target_interval_efficiency": float(final["target_interval_efficiency"]),
            "channeling_index": float(final["channeling_index"]),
            "mixing_index": float(final["mixing_index"]),
            "instability_index": float(final["instability_index"]),
        }
    finally:
        # 确保恢复原始模型参数
        setattr(base_model, "FLUIDS", original_fluids)
        setattr(base_model, "boundary_state", original_boundary_state)


def validate_levels(values: Sequence[float], name: str) -> list[float]:
    """
    校验并规范化敏感性分析因素水平。

    Args:
        values: 待校验的因素水平序列。
        name: 因素名称，用于错误提示。

    Returns:
        去重并按升序排列后的正值列表。

    Raises:
        ValueError: 当输入为空、非有限数或非正数时抛出。
    """
    if not values:
        raise ValueError(f"{name}至少需要输入一个数值")

    clean_values: list[float] = []
    for value in values:
        number = float(value)
        if not np.isfinite(number):
            raise ValueError(f"{name}包含非有限数值：{value}")
        if number <= 0.0:
            raise ValueError(f"{name}必须为正数：{value}")
        clean_values.append(number)
    return sorted(set(clean_values))


def baseline_result() -> dict[str, float | str]:
    """运行原始模型基线场景并返回关键指标。"""
    base_n = base_model.FLUIDS["tail"].n
    base_k = base_model.FLUIDS["tail"].k
    base_pv_mpa_s, base_yp_pa = powerlaw_to_bingham(base_k, base_n)
    _, _, _, baseline_metrics = base_model.simulate()
    baseline = baseline_metrics.iloc[-1]
    return {
        "rate_m3_min": 1.30,
        "pv_mpa_s": float(base_pv_mpa_s),
        "yp_pa": float(base_yp_pa),
        "n_fit": float(base_n),
        "k_fit": float(base_k),
        "effective_efficiency": float(baseline["effective_efficiency"]),
        "cbl_eval_interval_efficiency": float(baseline["cbl_eval_interval_efficiency"]),
        "target_interval_efficiency": float(baseline["target_interval_efficiency"]),
        "channeling_index": float(baseline["channeling_index"]),
        "mixing_index": float(baseline["mixing_index"]),
        "instability_index": float(baseline["instability_index"]),
        "scenario_type": "baseline_original",
    }


def run_factorial_cases(
    rate_levels_m3_min: Sequence[float],
    pv_levels_mpa_s: Sequence[float],
    yp_levels_pa: Sequence[float],
    *,
    include_baseline: bool = False,
    scenario_type: str = "custom_factorial",
    progress_callback: Callable[[int, int, float, float, float], None] | None = None,
) -> pd.DataFrame:
    """
    按输入的排量、塑性黏度和屈服值水平运行全因子场景。

    Args:
        rate_levels_m3_min: 排量水平，单位 m³/min。
        pv_levels_mpa_s: 塑性黏度水平，单位 mPa·s。
        yp_levels_pa: 屈服值水平，单位 Pa。
        include_baseline: 是否在结果中追加原始模型基线场景。
        scenario_type: 自定义场景类型标签。
        progress_callback: 场景完成后的回调函数，参数依次为
            (已完成场景数, 总场景数, 排量, 塑性黏度, 屈服值)。

    Returns:
        场景结果表，每行包含输入参数和最终全井段有效顶替效率等指标。
    """
    rates = validate_levels(rate_levels_m3_min, "排量")
    pvs = validate_levels(pv_levels_mpa_s, "塑性黏度")
    yps = validate_levels(yp_levels_pa, "屈服值")

    rows: list[dict[str, float | str]] = []
    if include_baseline:
        rows.append(baseline_result())

    total_cases = len(rates) * len(pvs) * len(yps)
    completed_cases = 0
    for rate_m3_min, pv_mpa_s, yp_pa in product(rates, pvs, yps):
        case = run_case(rate_m3_min=rate_m3_min, pv_mpa_s=pv_mpa_s, yp_pa=yp_pa)
        case["scenario_type"] = scenario_type
        rows.append(case)
        completed_cases += 1
        if progress_callback is not None:
            progress_callback(completed_cases, total_cases, rate_m3_min, pv_mpa_s, yp_pa)

    return pd.DataFrame(rows)


def summarize_factorial_results(
    result_df: pd.DataFrame,
    *,
    target_metric: str = FULL_WELL_TARGET_METRIC,
    top_ratio: float = 0.20,
) -> dict[str, object]:
    """
    汇总场景计算结果，给出影响排序、最优场景和建议参数区间。

    Args:
        result_df: run_factorial_cases生成的场景结果表。
        target_metric: 用于排序和推荐的目标指标。
        top_ratio: 提取建议区间时选取的高性能场景比例。

    Returns:
        可序列化的结果摘要字典。
    """
    if target_metric not in result_df.columns:
        raise ValueError(f"结果表缺少目标指标列：{target_metric}")

    scenario_df = pd.DataFrame(result_df.loc[result_df["scenario_type"] != "baseline_original"]).copy()
    if scenario_df.empty:
        raise ValueError("至少需要一个非基线场景用于结果汇总")

    ranking_df = factor_importance(scenario_df, target_metric)
    recommendation = suggest_intervals(scenario_df, target_metric, top_ratio=top_ratio)
    best_row = scenario_df.sort_values(by=target_metric, ascending=False).iloc[0]

    return {
        "target_metric": target_metric,
        "scenario_count": int(len(scenario_df)),
        "factors": {
            "rate_m3_min": sorted(float(value) for value in scenario_df["rate_m3_min"].drop_duplicates().tolist()),
            "pv_mpa_s": sorted(float(value) for value in scenario_df["pv_mpa_s"].drop_duplicates().tolist()),
            "yp_pa": sorted(float(value) for value in scenario_df["yp_pa"].drop_duplicates().tolist()),
        },
        "impact_ranking": ranking_df.to_dict(orient="records"),
        "best_case": {
            "rate_m3_min": float(best_row["rate_m3_min"]),
            "pv_mpa_s": float(best_row["pv_mpa_s"]),
            "yp_pa": float(best_row["yp_pa"]),
            "effective_efficiency": float(best_row["effective_efficiency"]),
            "cbl_eval_interval_efficiency": float(best_row["cbl_eval_interval_efficiency"]),
            "target_interval_efficiency": float(best_row["target_interval_efficiency"]),
        },
        "recommended_ranges": recommendation,
    }


def factor_importance(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    """
    按主效应均值跨度评估因素重要度。

    方法：对每个因素，计算各水平下目标指标的均值，
    取最大均值与最小均值之差（主效应跨度）作为该因素的影响程度指标。
    跨度越大，说明该因素对目标指标的影响越显著。

    同时计算整体斜率（最高水平均值 - 最低水平均值）/ (水平间距)，
    反映因素对指标的影响方向（正/负）。

    Args:
        df: 全因子场景结果DataFrame
        metric: 目标指标列名

    Returns:
        因素重要度排序DataFrame，包含：
        - factor: 因素列名
        - factor_name: 因素中文名称
        - mean_effect_span: 主效应跨度（各水平均值的极差）
        - overall_slope: 整体斜率
        - importance_share: 贡献占比（该因素跨度占总跨度的比例）
        - rank: 重要度排名
    """
    rows: list[dict[str, float | str]] = []
    factors = {
        "rate_m3_min": "排量(m3/min)",
        "pv_mpa_s": "塑性黏度(mPa·s)",
        "yp_pa": "屈服值(Pa)",
    }
    for factor_col, factor_name in factors.items():
        # 按因素水平分组计算目标指标均值
        factor_values = sorted(float(value) for value in df[factor_col].drop_duplicates().tolist())
        metric_means = [float(df.loc[df[factor_col] == value, metric].mean()) for value in factor_values]
        # 主效应跨度 = 最大水平均值 - 最小水平均值
        span = max(metric_means) - min(metric_means)
        # 整体斜率 = (最高水平均值 - 最低水平均值) / (水平间距)
        factor_delta = factor_values[-1] - factor_values[0]
        slope = 0.0 if factor_delta == 0.0 else (metric_means[-1] - metric_means[0]) / factor_delta
        rows.append({
            "factor": factor_col,
            "factor_name": factor_name,
            "mean_effect_span": span,
            "overall_slope": slope,
            "level_count": float(len(factor_values)),
            "note": "仅一个水平，无法评估影响" if len(factor_values) == 1 else "",
        })
    # 按主效应跨度降序排列
    out = pd.DataFrame(rows).sort_values("mean_effect_span", ascending=False).reset_index(drop=True)
    # 计算贡献占比
    total_span = float(out["mean_effect_span"].sum())
    out["importance_share"] = out["mean_effect_span"] / total_span if total_span > 0 else 0.0
    out["rank"] = np.arange(1, len(out) + 1)
    return out


def suggest_intervals(df: pd.DataFrame, metric: str, top_ratio: float = 0.20) -> dict[str, dict[str, float | str]]:
    """
    从高性能场景提取建议参数区间。

    方法：选取目标指标排名前top_ratio的场景，
    提取这些场景中各因素参数的最小值和最大值作为建议区间。

    Args:
        df: 全因子场景结果DataFrame
        metric: 目标指标列名
        top_ratio: 选取比例（默认20%）

    Returns:
        建议参数区间字典，包含：
        - selection_rule: 选取规则（场景数、最大指标值、指标名）
        - rate_m3_min: 排量建议区间 {min, max}
        - pv_mpa_s: 塑性黏度建议区间 {min, max}
        - yp_pa: 屈服值建议区间 {min, max}
    """
    top_n = min(len(df), max(1, int(np.ceil(len(df) * top_ratio))))
    top_df = df.nlargest(top_n, metric)
    max_metric = float(df[metric].max())

    return {
        "selection_rule": {
            "top_n": float(top_n),
            "max_metric": max_metric,
            "metric_name": metric,
        },
        "rate_m3_min": {
            "min": float(top_df["rate_m3_min"].min()),
            "max": float(top_df["rate_m3_min"].max()),
        },
        "pv_mpa_s": {
            "min": float(top_df["pv_mpa_s"].min()),
            "max": float(top_df["pv_mpa_s"].max()),
        },
        "yp_pa": {
            "min": float(top_df["yp_pa"].min()),
            "max": float(top_df["yp_pa"].max()),
        },
    }


def main() -> None:
    """
    主入口：执行三因素敏感性分析。

    流程：
    1. 将原始幂律参数等效转换为Bingham参数，作为基线参考
    2. 运行原始模型基线场景（不改流变、不改排量）
    3. 运行3×3×3=27个全因子场景
    4. 分析因素重要度排序
    5. 提取建议参数区间和最优场景
    6. 保存结果到CSV、JSON和Markdown文件
    """
    # 步骤1-4：运行原始基线与三因素全因子场景（3×3×3=27个）
    result_df = run_factorial_cases(
        RATE_LEVELS_M3_MIN,
        PV_LEVELS_MPA_S,
        YP_LEVELS_PA,
        include_baseline=True,
        scenario_type="factorial_3x3x3",
    )
    factorial_df = pd.DataFrame(result_df.loc[result_df["scenario_type"] == "factorial_3x3x3"]).copy()
    baseline_row = result_df[result_df["scenario_type"] == "baseline_original"].iloc[0]

    # 步骤5：因素重要度分析
    ranking_df = factor_importance(factorial_df, TARGET_METRIC)
    recommendation = suggest_intervals(factorial_df, TARGET_METRIC, top_ratio=0.20)

    # 最优场景
    best_row = factorial_df.sort_values(by=TARGET_METRIC, ascending=False).iloc[0]

    # 保存场景结果CSV
    result_df.to_csv(SCENARIO_CSV, index=False, encoding="utf-8-sig")

    # 步骤6：生成结果摘要
    summary = {
        "analysis_name": "呼102尾管_D2DGA 三因素敏感性分析",
        "factors": {
            "rate_m3_min": RATE_LEVELS_M3_MIN,
            "pv_mpa_s": PV_LEVELS_MPA_S,
            "yp_pa": YP_LEVELS_PA,
        },
        "target_metric": TARGET_METRIC,
        "baseline_equivalent_rheology": {
            "tail_power_law_n": float(baseline_row["n_fit"]),
            "tail_power_law_k": float(baseline_row["k_fit"]),
            "equivalent_pv_mpa_s": float(baseline_row["pv_mpa_s"]),
            "equivalent_yp_pa": float(baseline_row["yp_pa"]),
            "note": "由原始幂律参数在300/600rpm点等效反算。",
        },
        "impact_ranking": ranking_df.to_dict(orient="records"),
        "best_case": {
            "rate_m3_min": float(best_row["rate_m3_min"]),
            "pv_mpa_s": float(best_row["pv_mpa_s"]),
            "yp_pa": float(best_row["yp_pa"]),
            "cbl_quality_proxy": float(best_row["cbl_quality_proxy"]),
            "cbl_eval_interval_efficiency": float(best_row["cbl_eval_interval_efficiency"]),
            "effective_efficiency": float(best_row["effective_efficiency"]),
        },
        "recommended_ranges": recommendation,
        "outputs": {
            "scenario_csv": str(SCENARIO_CSV),
            "summary_json": str(SUMMARY_JSON),
            "summary_md": str(SUMMARY_MD),
        },
    }

    SUMMARY_JSON.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    # 生成Markdown格式的影响排序行
    rank_lines = [
        f"{int(row['rank'])}. {row['factor_name']} | 主效应跨度={row['mean_effect_span']:.5f} | 贡献占比={row['importance_share']:.2%}"
        + (f" | {row['note']}" if row.get("note") else "")
        for _, row in ranking_df.iterrows()
    ]

    # 生成Markdown摘要
    md_lines = [
        "# 呼102尾管 D2DGA 三因素敏感性分析",
        "",
        "## 分析设置",
        f"- 场景规模：3×3×3 全因子，共 {len(factorial_df)} 个场景",
        f"- 目标指标：{TARGET_METRIC}（CBL评价井段质量响应效率，无量纲）",
        f"- 排量水平（m3/min）：{RATE_LEVELS_M3_MIN}",
        f"- 塑性黏度水平（mPa·s）：{PV_LEVELS_MPA_S}",
        f"- 屈服值水平（Pa）：{YP_LEVELS_PA}",
        "",
        "## 影响排序（按主效应跨度）",
        *[f"- {line}" for line in rank_lines],
        "",
        "## 最优场景（以目标指标最大为准）",
        f"- 排量：{best_row['rate_m3_min']:.2f} m3/min",
        f"- 塑性黏度：{best_row['pv_mpa_s']:.1f} mPa·s",
        f"- 屈服值：{best_row['yp_pa']:.1f} Pa",
        f"- CBL质量响应效率：{best_row['cbl_quality_proxy']:.4f}",
        f"- CBL评价井段水动力效率：{best_row['cbl_eval_interval_efficiency']:.4f}",
        "",
        "## 建议参数区间（取目标指标前20%场景）",
        f"- 排量：{recommendation['rate_m3_min']['min']:.2f} - {recommendation['rate_m3_min']['max']:.2f} m3/min",
        f"- 塑性黏度：{recommendation['pv_mpa_s']['min']:.1f} - {recommendation['pv_mpa_s']['max']:.1f} mPa·s",
        f"- 屈服值：{recommendation['yp_pa']['min']:.1f} - {recommendation['yp_pa']['max']:.1f} Pa",
        "",
        "## 说明",
        "- 塑性黏度/屈服值通过300/600rpm点等效拟合为幂律参数后接入原模型。",
        "- 本分析仅评估顶替效率响应，不包含泵压窗口、ECD上限与漏失风险约束。",
    ]
    SUMMARY_MD.write_text("\n".join(md_lines), encoding="utf-8")

    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
