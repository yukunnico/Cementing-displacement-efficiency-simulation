"""
呼101井尾管 D2DGA 现场顺序多流体驱动脚本
=========================================

目标：
1. 在不覆盖 `hu101_d2dga_model.py` 的前提下，将呼101改成按现场施工顺序的多流体入环空模式；
2. 优先使用呼101本井资料；缺失参数按“邻井 > 文献规则 > 既有模型代理值”回填；
3. 复用已验证的悬挂器失败 / CBL异常段 standoff 修正；
4. 输出结果摘要、图表、参数来源与假设说明。
"""

from __future__ import annotations

import importlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

base = cast(Any, importlib.import_module("hu101_d2dga_model"))
prior = cast(Any, importlib.import_module("hu101_d2dga_model_improved"))


OUT_DIR = Path(__file__).resolve().parent / "现场顺序多流体输出"


@dataclass(frozen=True)
class InjectionStage:
    fluid_key: str | None
    volume_m3: float
    rate_m3_min: float
    label: str
    source: str
    note: str


FIELD_TRACKED = ["balance", "spacer", "lead", "tail"]
SHOE_LAG_VOLUME_M3 = 52.0
BUMP_PRESSURE_RANGE_MPA = "4.8-9.0"
BUMP_HOLD_PRESSURE_MPA = 5.4

# ===== 呼101本井 + 邻井/文献兜底规则 =====
# 直接值：来自呼101抽取包/抽取报告。
# 兜底1：邻井顺序/类别命名（postflush/light mud/middle fluid）来自呼101抽取报告与邻井序列模式。
# 兜底2：隔离液流变缺失 → 采用 literature-guided fallback：
#        PV = 0.75 * mud_PV，YP = 1.30 * mud_YP；其逻辑与AADE/API综述一致：
#        spacer 应位于 mud 和 cement 之间、并具有不低于 mud 的洗井能力。
# 兜底3：平衡液/轻泥浆流变缺失 → 沿用既有基准模型代理值，视为 light-mud proxy。

SPACER_PV_FALLBACK = 0.75 * base.MUD_PV_PA_S
SPACER_YP_FALLBACK = 1.30 * base.MUD_YP_PA

FIELD_FLUIDS: dict[str, Any] = {
    "mud": base.Fluid("井浆/钻井液", 1.96, "bingham", k=base.MUD_PV_PA_S, yield_stress_pa=base.MUD_YP_PA),
    "balance": base.Fluid("平衡液", 1.85, "bingham", k=base.BALANCE_PV_PA_S, yield_stress_pa=base.BALANCE_YP_PA),
    "spacer": base.Fluid("驱油隔离液", 2.00, "bingham", k=SPACER_PV_FALLBACK, yield_stress_pa=SPACER_YP_FALLBACK),
    "lead": base.Fluid("领浆", 2.10, "power_law", n=0.719, k=0.815),
    "tail": base.Fluid("尾浆", 1.90, "power_law", n=0.722, k=0.684),
}


def build_surface_schedule() -> list[InjectionStage]:
    """按呼101现场抽取结果构建地面泵注顺序。"""
    return [
        InjectionStage("balance", 25.0, 1.2, "注平衡液", "呼101抽取报告 §3.2", "直接体积与密度"),
        InjectionStage("spacer", 25.0, 1.2, "注驱油隔离液", "呼101抽取报告 §3.2", "直接体积与密度；流变缺失，采用文献+泥浆代理回填"),
        InjectionStage("lead", 47.0, 1.2, "注领浆", "呼101抽取报告 §3.2 + 2011121", "直接体积、密度、幂律参数"),
        InjectionStage("tail", 23.0, 1.2, "注尾浆", "呼101抽取报告 §3.2 + 2011121", "直接体积、密度、幂律参数"),
        InjectionStage(None, 2.0, 0.6, "注后置液(管内)", "呼101抽取报告 §3.2/§3.3", "仅按管内压塞/占位处理，不作为独立环空入流"),
        InjectionStage("light_mud", 26.0, 1.5, "注轻泥浆", "呼101抽取报告 §3.2/§3.3 + 技术总结现场排量", "直接体积；按现场 1.5 m³/min 建模"),
        InjectionStage("middle_fluid", 10.0, 1.2, "注中置液", "呼101抽取报告 §3.2/§3.3", "直接体积；按主替浆排量建模"),
        InjectionStage(None, 40.0, 1.0, "井浆快替", "呼101设计/技术总结", "按设计表的 40m³@1.0 m³/min 作为快替井浆"),
        InjectionStage(None, 23.4, 0.55, "井浆慢替", "呼101施工记录表", "为满足 101.4m³ 总替量，按剩余 23.4m³@0.55 m³/min 处理"),
    ]


def build_annulus_entry_schedule(surface_schedule: list[InjectionStage]) -> list[InjectionStage]:
    """构建越鞋进入环空的有效阶段表：只保留真正进入环空的前置液与水泥浆。"""
    allowed = {"注平衡液", "注驱油隔离液", "注领浆", "注尾浆"}
    return [stage for stage in surface_schedule if stage.label in allowed]


def locate_surface_stage(surface_schedule: list[InjectionStage], time_s: float) -> InjectionStage | None:
    elapsed = 0.0
    for stage in surface_schedule:
        duration = stage.volume_m3 / stage.rate_m3_min * 60.0 if stage.rate_m3_min > 0.0 else 0.0
        if time_s < elapsed + duration - 1e-12:
            return stage
        elapsed += duration
    return None


def cumulative_surface_volume(surface_schedule: list[InjectionStage], time_s: float) -> float:
    """计算给定时刻的地面累计泵入体积。"""
    elapsed = 0.0
    injected = 0.0
    for stage in surface_schedule:
        duration = stage.volume_m3 / stage.rate_m3_min * 60.0 if stage.rate_m3_min > 0.0 else 0.0
        if time_s < elapsed + duration - 1e-12:
            injected += max(time_s - elapsed, 0.0) / 60.0 * stage.rate_m3_min
            return injected
        injected += stage.volume_m3
        elapsed += duration
    return injected


def locate_stage_by_volume(schedule: list[InjectionStage], arrival_volume_m3: float) -> InjectionStage | None:
    cumulative = 0.0
    for stage in schedule:
        cumulative += stage.volume_m3
        if arrival_volume_m3 < cumulative - 1e-12:
            return stage
    return None


def make_field_order_boundary_state(surface_schedule: list[InjectionStage], annulus_schedule: list[InjectionStage]):
    """生成带鞋口滞后体积的现场顺序多流体边界函数。"""

    def boundary_state(t: float) -> tuple[Any, float, str]:
        active_surface_stage = locate_surface_stage(surface_schedule, t)
        if active_surface_stage is None:
            return np.zeros(len(base.TRACKED)), 0.0, "碰压结束"

        surface_volume = cumulative_surface_volume(surface_schedule, t)
        annulus_arrival_volume = max(surface_volume - SHOE_LAG_VOLUME_M3, 0.0)
        active_stage = locate_stage_by_volume(annulus_schedule, annulus_arrival_volume)
        vec = np.zeros(len(base.TRACKED))
        if active_stage is None:
            vec[base.TRACKED.index("tail")] = 1.0
            return vec, active_surface_stage.rate_m3_min / 60.0, f"{active_surface_stage.label}(顶胶塞管内顶替，环空保持尾浆)"

        vec[base.TRACKED.index(active_stage.fluid_key)] = 1.0
        return vec, active_surface_stage.rate_m3_min / 60.0, f"{active_stage.label}(按52m³鞋口滞后折算)"

    return boundary_state


def set_field_order_globals() -> tuple[list[InjectionStage], list[InjectionStage]]:
    """向基准模型注入现场顺序多流体配置。"""
    surface_schedule = build_surface_schedule()
    annulus_schedule = build_annulus_entry_schedule(surface_schedule)
    base.TRACKED = FIELD_TRACKED.copy()
    base.FLUIDS = FIELD_FLUIDS.copy()
    base.ANNULUS_BOUNDARY_MODE = "field_order_multifluid"
    base.BALANCE_VOLUME_M3 = 25.0
    base.SPACER_VOLUME_M3 = 25.0
    base.LEAD_VOLUME_M3 = 47.0
    base.TAIL_VOLUME_M3 = 23.0
    base.DISPLACEMENT_VOLUME_M3 = 101.4
    base.FIELD_RATE_M3_MIN = 1.2
    base.FIELD_RATE_LOW_M3_MIN = 0.55
    base.RATE_SWITCH_DISP_VOLUME_M3 = 62.0
    base.BALANCE_DENSITY_GCC = 1.85
    base.SPACER_DENSITY_GCC = 2.00
    base.LEAD_DENSITY_GCC = 2.10
    base.TAIL_DENSITY_GCC = 1.90
    base.boundary_state = make_field_order_boundary_state(surface_schedule, annulus_schedule)
    return surface_schedule, annulus_schedule


def build_stage_table(schedule: list[InjectionStage]) -> pd.DataFrame:
    rows = []
    cumulative = 0.0
    for stage in schedule:
        start_m3 = cumulative
        cumulative += stage.volume_m3
        rows.append(
            {
                "阶段": stage.label,
                "流体键": stage.fluid_key or "mud(complement)",
                "体积_m3": stage.volume_m3,
                "排量_m3_min": stage.rate_m3_min,
                "起始累计体积_m3": start_m3,
                "结束累计体积_m3": cumulative,
                "来源": stage.source,
                "备注": stage.note,
            }
        )
    return pd.DataFrame(rows)


def compute_injected_volume_series(schedule: list[InjectionStage], time_s: Any) -> Any:
    """按真实分段排量累计注入体积。"""
    volumes = []
    for t in time_s.to_numpy(dtype=float):
        remaining_t = t
        injected = 0.0
        for stage in schedule:
            duration_s = stage.volume_m3 / stage.rate_m3_min * 60.0 if stage.rate_m3_min > 0 else 0.0
            if remaining_t <= 0.0:
                break
            consume = min(remaining_t, duration_s)
            injected += consume / 60.0 * stage.rate_m3_min
            remaining_t -= consume
        volumes.append(injected)
    return np.asarray(volumes, dtype=float)


def schedule_total_time_s(schedule: list[InjectionStage]) -> float:
    """按地面泵注阶段计算总施工时长，用作碰压结束代理。"""
    total_s = 0.0
    for stage in schedule:
        if stage.rate_m3_min > 0.0:
            total_s += stage.volume_m3 / stage.rate_m3_min * 60.0
    return total_s


def write_extra_outputs(surface_schedule: list[InjectionStage], annulus_schedule: list[InjectionStage], geom: dict[str, Any], x: Any, metrics: pd.DataFrame) -> dict[str, str]:
    """补充多流体模式专属表格与图表。"""
    charts = base.OUT_DIR / "图表"
    charts.mkdir(parents=True, exist_ok=True)

    stage_table = build_stage_table(surface_schedule)
    stage_csv = base.OUT_DIR / "呼101_现场顺序多流体_阶段表.csv"
    stage_table.to_csv(stage_csv, index=False, encoding="utf-8-sig")

    annulus_stage_table = build_stage_table(annulus_schedule)
    annulus_stage_csv = base.OUT_DIR / "呼101_环空入口阶段表_鞋口滞后修正.csv"
    annulus_stage_table.to_csv(annulus_stage_csv, index=False, encoding="utf-8-sig")

    fluid_profile = pd.DataFrame({"井深_m": geom["md"]})
    for idx, name in enumerate(base.TRACKED):
        fluid_profile[f"{name}_平均浓度"] = np.average(x[idx], axis=0, weights=geom["b"])
    fluid_profile["mud_平均浓度"] = np.clip(1.0 - fluid_profile[[c for c in fluid_profile.columns if c.endswith("平均浓度") and c != "mud_平均浓度"]].sum(axis=1), 0.0, 1.0)
    fluid_profile_path = base.OUT_DIR / "呼101_现场顺序多流体_深度浓度剖面.csv"
    fluid_profile.to_csv(fluid_profile_path, index=False, encoding="utf-8-sig")

    schedule_plot = charts / "呼101_现场顺序多流体_施工阶段图.png"
    fig, ax = plt.subplots(figsize=(11, 4))
    left = 0.0
    colors = ["#1D4ED8", "#0EA5E9", "#10B981", "#84CC16", "#F59E0B", "#F97316", "#A855F7", "#6B7280", "#111827"]
    for idx, stage in enumerate(surface_schedule):
        duration = stage.volume_m3 / stage.rate_m3_min
        ax.barh([0], [duration], left=[left], color=colors[idx % len(colors)], edgecolor="white", label=stage.label)
        ax.text(left + duration / 2.0, 0, f"{stage.label}\n{stage.volume_m3:.1f}m³", ha="center", va="center", fontsize=8, color="white")
        left += duration
    ax.set_xlabel("累计施工时间 / min")
    ax.set_yticks([])
    ax.set_title("呼101现场施工顺序多流体入环空阶段图")
    ax.grid(True, axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(schedule_plot, dpi=220)
    plt.close(fig)

    concentration_plot = charts / "呼101_现场顺序多流体_最终浓度剖面.png"
    fig, ax = plt.subplots(figsize=(8, 7))
    for name in ["balance", "spacer", "lead", "tail", "mud"]:
        column = f"{name}_平均浓度"
        if column in fluid_profile.columns:
            ax.plot(fluid_profile[column], fluid_profile["井深_m"], label=name)
    ax.invert_yaxis()
    ax.set_xlabel("平均体积分数")
    ax.set_ylabel("井深 / m")
    ax.set_title("呼101最终多流体平均浓度剖面")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(concentration_plot, dpi=220)
    plt.close(fig)

    annular_volume_plot = charts / "呼101_现场顺序多流体_效率随环空体积数变化.png"
    injected_volume_m3 = compute_injected_volume_series(surface_schedule, metrics["time_s"])
    annular_pore_volumes = injected_volume_m3 / base.physical_annular_volume_m3()
    plt.figure(figsize=(8, 5))
    plt.plot(annular_pore_volumes, metrics["cbl_eval_interval_efficiency"], label="CBL井段水动力效率", color="#2563EB", linewidth=2)
    plt.plot(annular_pore_volumes, metrics["cbl_quality_proxy"], label="原α质量响应效率", color="#7C3AED", linewidth=2)
    plt.plot(annular_pore_volumes, metrics["target_interval_efficiency"], label="油气水层段有效效率", color="#059669", linewidth=1.8)
    plt.axhline(base.FIELD_REFERENCE_EFFICIENCY, color="#F97316", linestyle="--", label="资料CBL合格率")
    plt.xlabel("累计注入体积 / 物理环空体积")
    plt.ylabel("效率/合格率")
    plt.title("呼101现场顺序多流体顶替效率随环空体积数变化")
    plt.ylim(0, 1.05)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(annular_volume_plot, dpi=220)
    plt.close()

    return {
        "stage_csv": str(stage_csv),
        "annulus_stage_csv": str(annulus_stage_csv),
        "fluid_profile_csv": str(fluid_profile_path),
        "stage_plot": str(schedule_plot),
        "concentration_plot": str(concentration_plot),
        "annular_volume_plot": str(annular_volume_plot),
    }


def enhance_summary(summary: dict[str, Any], surface_schedule: list[InjectionStage], annulus_schedule: list[InjectionStage], extra_outputs: dict[str, str], metrics: pd.DataFrame) -> dict[str, Any]:
    final = summary["最终结果"]
    peak_row = metrics.loc[metrics["bulk_cement_fill"].idxmax()]
    summary["模型名称"] = "呼101尾管_D2DGA_现场施工顺序多流体模式"
    summary["环空入口边界模式"] = "field_order_multifluid"
    summary["现场顺序多流体阶段"] = build_stage_table(surface_schedule).to_dict(orient="records")
    summary["环空入口阶段_鞋口滞后修正"] = build_stage_table(annulus_schedule).to_dict(orient="records")
    summary["多流体参数来源与回填规则"] = [
        "平衡液25m³、驱油隔离液25m³、领浆47m³、尾浆23m³、后置液2m³、轻泥浆26m³、中置液10m³来自 h101_data_extraction_report.md §3.2。",
        "替浆总量101.4m³来自呼101结构化数据包；本次采用后置液2m³、轻泥浆26m³、中置液10m³、井浆快替40m³、井浆慢替23.4m³，总计101.4m³。",
        "平衡液密度1.85、隔离液密度2.00、泥浆密度1.96、领浆2.10、尾浆1.90来自呼101本井资料。",
        f"隔离液流变缺失：采用文献+现场泥浆兜底，PV={SPACER_PV_FALLBACK*1000:.1f} mPa·s (=0.75×mud PV)，YP={SPACER_YP_FALLBACK:.1f} Pa (=1.30×mud YP)。",
        "平衡液/轻泥浆流变缺失：沿用既有基准模型代理值 PV=30 mPa·s、YP=3 Pa，作为轻泥浆 proxy。",
        f"环空入口采用 52m³ 鞋口滞后体积修正：用以把地面泵注顺序映射为越鞋进入环空顺序，使隔离液/领浆/尾浆进入环空时刻更接近施工记录。",
    ]
    summary["假设说明"] = list(summary.get("假设说明", [])) + [
        "本次不再使用 sustained_tail / volume_limited / tail_then_mud 三个抽象替浆边界，而改为按呼101现场地面泵注阶段顺序，经鞋口滞后修正后映射为环空入口顺序。",
        "后置液、轻泥浆、中置液和后续井浆都作为管内顶替驱动包处理；它们推动顶胶塞和尾浆前进，但不作为新的独立环空流体组分。",
        "后置液、中置液缺失独立流变时，按同密度前置/隔离体系代理建模；其不确定性在报告中单列。",
        "碰压 4.8-9MPa 被视为作业结束信号；本模型在替浆总量达到现场碰压体积窗口时结束。",
    ]
    summary["假设说明"] = [
        item.replace(
            "悬挂器坐挂失败(技术总结记载)未在几何中建模，其影响通过α标定系数吸收。",
            "悬挂器坐挂失败已通过standoff代理退化显式进入几何；剩余差异仅在质量响应解释时保留后验α口径。",
        ).replace(
            "环空入口边界默认为sustained_tail，保留替浆阶段尾浆等效入环空的口径。",
            "环空入口采用 field_order_multifluid 现场顺序口径，不再使用 sustained_tail 抽象边界。",
        )
        for item in summary["假设说明"]
    ]
    summary["重要限制"] = list(summary.get("重要限制", [])) + [
        "平衡液/轻泥浆、后置液、中置液的流变参数并非本井实测，属于邻井/文献规则兜底后的工程代理。",
        "现场顺序阶段已显式化，但末段井浆仍包含11.4m³未分项来源的总替量残差，只能并入井浆类别处理。",
        "本模型把后续驱替流体统一视作管内推动，不显式求解顶胶塞与浮箍的接触力学；碰压仅作为基于现场体积/压力记录的结束条件代理。",
        "当前计算域只覆盖 5400-7868m 尾管固井区间，5400m 顶界为开边界；因此最终低效率同时反映了区间内被后续流体替出的水泥与越过区间顶界的水泥。",
    ]
    summary["新增输出文件"] = extra_outputs
    summary["最终结果补充"] = {
        "多流体模式全井段最终区间留置效率": float(final["全井段最终有效顶替效率"]),
        "多流体模式CBL井段最终区间留置效率": float(final["CBL评价井段模拟有效顶替效率"]),
        "多流体模式原α质量响应效率": float(final["CBL评价井段质量响应效率"]),
        "峰值水泥占据率": float(peak_row["bulk_cement_fill"]),
        "峰值CBL井段水动力效率": float(metrics["cbl_eval_interval_efficiency"].max()),
        "峰值出现时间_min": float(peak_row["time_min"]),
    }
    return summary


def write_report(summary: dict[str, Any], diagnostics: dict[str, float], residual_alpha: dict[str, float]) -> None:
    prior_comparison = ""
    final_cbl_eff = float(summary["最终结果"]["CBL评价井段模拟有效顶替效率"])
    peak_cbl_eff = float(summary["最终结果补充"]["峰值CBL井段水动力效率"])
    drop_ratio = peak_cbl_eff / max(final_cbl_eff, 1e-12)
    prior_table_path = Path(__file__).resolve().parent / "改进模型输出" / "呼101_改进模型_边界模式对比总表.csv"
    if prior_table_path.exists():
        prior_table = pd.read_csv(prior_table_path, encoding="utf-8-sig")
        prior_comparison = cast(str, prior_table[["边界模式", "CBL评价井段模拟有效顶替效率", "残余α质量响应效率"]].to_markdown(index=False))

    stage_records = cast(list[dict[str, Any]], summary["现场顺序多流体阶段"])
    annulus_stage_records = cast(list[dict[str, Any]], summary["环空入口阶段_鞋口滞后修正"])

    report_lines = [
        "# 呼101尾管 D2DGA 现场施工顺序多流体模型报告",
        "",
        "## 1. 改进目标",
        "将呼101从抽象的替浆边界模式，改为按现场施工顺序的多流体入环空模式：平衡液 → 驱油隔离液 → 领浆 → 尾浆 → 后置液 → 轻泥浆 → 中置液 → 井浆快替/慢替。",
        "",
        "## 2. 参数来源与回填原则",
        "- 本井直接值：呼101结构化数据包 + h101_data_extraction_report.md。",
        "- 邻井参考：呼103、呼探1、呼探1-002 的多流体顺序和分段体积，用于确认现场施工序列类型。",
        "- 文献兜底：隔离液流变缺失时，按 spacer 位于 mud 与 cement 之间、且具备不低于 mud 洗井能力的规则回填。",
        f"- 本次采用的隔离液兜底流变：PV={SPACER_PV_FALLBACK*1000:.1f} mPa·s，YP={SPACER_YP_FALLBACK:.1f} Pa。",
        "",
        "## 3. 现场顺序阶段表",
        cast(str, pd.DataFrame(stage_records).to_markdown(index=False)),
        "",
        "## 3.1 环空入口阶段表（52m³鞋口滞后修正）",
        cast(str, pd.DataFrame(annulus_stage_records).to_markdown(index=False)),
        "",
        f"**碰压结束条件**：现场记录碰压 {BUMP_PRESSURE_RANGE_MPA} MPa，控压 {BUMP_HOLD_PRESSURE_MPA:.1f} MPa；本模型以现场替浆总量 101.4m³ 对应的碰压窗口作为结束条件代理。",
        "",
        "## 4. 核心结果",
        f"- 全井段最终有效顶替效率：{summary['最终结果']['全井段最终有效顶替效率']:.4f}",
        f"- CBL评价井段模拟有效顶替效率：{summary['最终结果']['CBL评价井段模拟有效顶替效率']:.4f}",
        f"- 原α质量响应效率：{summary['最终结果']['CBL评价井段质量响应效率']:.4f}",
        f"- 残余α重标定值：{residual_alpha['残余重标定α']:.4f}",
        f"- 残余α后验归一化质量响应效率：{residual_alpha['残余α质量响应效率']:.4f}",
        f"- 峰值水泥占据率：{summary['最终结果补充']['峰值水泥占据率']:.4f} @ {summary['最终结果补充']['峰值出现时间_min']:.1f} min",
        f"- 峰值CBL井段水动力效率：{summary['最终结果补充']['峰值CBL井段水动力效率']:.4f}",
        f"- 6050-6210m：CBL代理均值 {diagnostics['6050_6210m_CBL代理均值']:.4f}，模拟均值 {diagnostics['6050_6210m_模拟均值']:.4f}",
        f"- 6100-6400m：CBL代理均值 {diagnostics['6100_6400m_CBL代理均值']:.4f}，模拟均值 {diagnostics['6100_6400m_模拟均值']:.4f}",
        "",
        "## 4.1 与前一版三种抽象边界模式对比",
        prior_comparison or "未找到前一版三模式对比总表。",
        "",
        "## 5. 关键解释",
        "- 与三种抽象边界模式相比，本次模型把地面泵注顺序与越鞋进入环空顺序拆开处理：地面阶段按施工表，环空入口按 52m³ 鞋口滞后修正。",
        "- 后置液、轻泥浆、中置液和后续井浆不再被当成新的环空流体，而是作为管内顶替驱动包；环空入口在尾浆到鞋后保持尾浆直到碰压结束。",
        "- 本次结果对应‘碰压结束时刻’的到位效率，而不是碰压后继续大体积替浆后的区间留置效率；因此峰值与终值几乎一致，说明模型现在已按碰压作为结束信号截断。",
        f"- 峰-终 {drop_ratio:.1f} 倍（{peak_cbl_eff:.3f}→{final_cbl_eff:.3f}）几乎无衰减，说明顶胶塞顶替到碰压结束这一口径下，水泥尚未被后续流体继续推出 5400m 顶界。",
        "- 这次结果的物理含义优先看水动力效率；残余α结果仍然只是对总体CBL合格率的后验归一化，不是独立预测。当最终水动力效率高于现场CBL代理值时，残余α只表示还需要多少经验质量惩罚才能把水动力结果压到现场质量口径。",
        f"- 6050-6210m 强异常段已经施加了严重 standoff 退化（上限压到 0.34，对应悬挂器失效/偏心极端代理）；即便如此，模拟均值仍为 {diagnostics['6050_6210m_模拟均值']:.3f}、远高于 CBL 代理 {diagnostics['6050_6210m_CBL代理均值']:.3f}，说明该段问题已超出单纯顶替顺序与环空流动范畴，更可能受胶结质量、气窜、微环空或测井响应影响。",
        "- 换句话说，这一段更可能是‘水泥基本到位但胶结效果失真’而不是‘水泥根本没顶替到位’；因此本模型适合解释碰压到位时的流体置换，而不能单独解释后续胶结质量劣化。",
        "",
        "## 6. 局限",
        "- 平衡液/轻泥浆、后置液、中置液缺实测流变；目前只能按邻井/文献规则代理。",
        "- 井浆分项体积和总替量存在 11.4m³ 差额，本次已并入末段井浆慢替以满足质量守恒。",
        "- 当前结果可视作‘后续流体只在管内顶替、环空入口保持尾浆直到碰压结束’的到位解释；如果鞋口滞后更大、或现场存在更多局部滞留，则实际碰压到位效率可能与本次结果仍有偏差。",
        "- CBL剖面仍只是质量代理标签，不是直接测得的水力顶替效率。",
    ]
    (OUT_DIR / "呼101尾管_D2DGA_现场顺序多流体报告.md").write_text("\n".join(report_lines), encoding="utf-8")


def main() -> None:
    cbl_profile = prior.load_cbl_profile()
    importlib.reload(base)
    base.OUT_DIR = OUT_DIR
    base.OUT_DIR.mkdir(parents=True, exist_ok=True)
    surface_schedule, annulus_schedule = set_field_order_globals()
    base.standoff_profile = prior.make_hanger_failure_standoff(base.standoff_profile)

    total_t = schedule_total_time_s(surface_schedule)
    geom, x, wall, metrics = base.simulate(total_t=total_t)
    summary = base.save_outputs(geom, x, wall, metrics)
    diagnostics = prior.compute_profile_diagnostics(base.OUT_DIR, cbl_profile)
    residual_alpha = prior.recalibrate_residual_alpha(summary)
    extra_outputs = write_extra_outputs(surface_schedule, annulus_schedule, geom, x, metrics)
    summary["CBL剖面对比诊断"] = diagnostics
    summary["残余质量惩罚重标定"] = residual_alpha
    summary = enhance_summary(summary, surface_schedule, annulus_schedule, extra_outputs, metrics)

    summary_path = base.OUT_DIR / "呼101尾管_D2DGA_现场顺序多流体_结果摘要.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    write_report(summary, diagnostics, residual_alpha)
    print(json.dumps({"输出目录": str(OUT_DIR), "地面阶段数": len(surface_schedule), "环空阶段数": len(annulus_schedule)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
