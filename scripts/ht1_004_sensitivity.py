"""
呼1-004井尾浆参数敏感性分析脚本

对尾浆的注入排量、屈服应力、塑性黏度进行全因子敏感性分析，
评估各参数对全井有效顶替效率的影响。

参数空间：
- 排量: 0.30~1.60 m^3/min, 步长0.10, 共14个水平
- 屈服应力: 5~25 Pa, 步长1, 共21个水平
- 塑性黏度: 50~350 mPa·s, 步长20, 共16个水平
- 全因子组合: 14 × 21 × 16 = 4704次模拟

输出：
- CSV汇总表: results/呼1-004_敏感性分析/呼1-004_敏感性分析结果.csv
- 单因素曲线图 ×3
- 两两交互热力图 ×3

使用方式：
    python scripts/ht1_004_sensitivity.py
    python scripts/ht1_004_sensitivity.py --resume  # 断点续跑
    python scripts/ht1_004_sensitivity.py --plot-only  # 仅生成图表
"""

from __future__ import annotations

import argparse
import dataclasses
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # 非交互式后端，适合服务器/批量运行
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

# 确保项目根目录在 PYTHONPATH 中
_PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.loaders.ht1_004_loader import load_ht1_004_tailpipe
from cemdisp.data.pumping_schedule import PumpingSchedule
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver


# ============================================================
# 中文字体配置
# ============================================================
def _setup_chinese_font() -> None:
    """配置 matplotlib 中文字体，优先使用 SimHei。"""
    plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False  # 负号正常显示


# ============================================================
# 参数空间定义
# ============================================================
# 排量范围 (m^3/min)
RATE_MIN = 0.30
RATE_MAX = 1.60
RATE_STEP = 0.10
RATE_BASELINE = 1.40

# 屈服应力范围 (Pa)
YP_MIN = 5
YP_MAX = 25
YP_STEP = 2
YP_BASELINE = 14.0

# 塑性黏度范围 (mPa·s)
PV_MIN = 50
PV_MAX = 350
PV_STEP = 20
PV_BASELINE = 180.0

# 生成参数值数组
RATE_VALUES = np.arange(RATE_MIN, RATE_MAX + RATE_STEP / 2, RATE_STEP)  # 14个
YP_VALUES = np.arange(YP_MIN, YP_MAX + YP_STEP / 2, YP_STEP)  # 21个
PV_VALUES = np.arange(PV_MIN, PV_MAX + PV_STEP / 2, PV_STEP)  # 16个

# 输出目录
OUTPUT_DIR = _PROJECT_ROOT / "results" / "呼1-004_敏感性分析"
CSV_PATH = OUTPUT_DIR / "呼1-004_敏感性分析结果.csv"
CHECKPOINT_PATH = OUTPUT_DIR / "呼1-004_敏感性分析_中间结果.csv"


# ============================================================
# 数据修改辅助函数
# ============================================================
def _find_tail_fluid_index(fluids: tuple[FluidSpec, ...]) -> int:
    """找到尾浆在 fluids 元组中的索引。"""
    for i, fluid in enumerate(fluids):
        if fluid.name == "尾浆":
            return i
    raise ValueError("未找到名为'尾浆'的流体")


def _find_tail_step_index(schedule: PumpingSchedule) -> int:
    """找到尾浆注入步骤在 PumpingSchedule 中的索引。"""
    for i, step in enumerate(schedule.steps):
        if step.step_name == "注入尾浆":
            return i
    raise ValueError("未找到名为'注入尾浆'的泵送步骤")


def _create_modified_data(
    base_fluids: tuple[FluidSpec, ...],
    base_schedule: PumpingSchedule,
    tail_yp_pa: float,
    tail_pv_pa_s: float,
    tail_rate_m3_min: float,
) -> tuple[tuple[FluidSpec, ...], PumpingSchedule]:
    """基于基线数据创建修改后的流体元组和泵送计划。

    参数：
        base_fluids: 基线流体元组
        base_schedule: 基线泵送计划
        tail_yp_pa: 尾浆屈服应力 (Pa)
        tail_pv_pa_s: 尾浆塑性黏度 (Pa·s)
        tail_rate_m3_min: 尾浆排量 (m^3/min)

    返回：
        修改后的 (fluids, schedule) 元组
    """
    # 修改尾浆流体参数
    tail_idx = _find_tail_fluid_index(base_fluids)
    new_tail = dataclasses.replace(
        base_fluids[tail_idx],
        plastic_viscosity_pa_s=tail_pv_pa_s,
        yield_stress_pa=tail_yp_pa,
    )
    new_fluids = tuple(new_tail if i == tail_idx else f for i, f in enumerate(base_fluids))

    # 修改尾浆排量
    tail_step_idx = _find_tail_step_index(base_schedule)
    new_step = dataclasses.replace(
        base_schedule.steps[tail_step_idx],
        rate_m3_min=tail_rate_m3_min,
    )
    new_steps = tuple(
        new_step if i == tail_step_idx else s
        for i, s in enumerate(base_schedule.steps)
    )
    new_schedule = PumpingSchedule(steps=new_steps)

    return new_fluids, new_schedule


# ============================================================
# 单次模拟运行
# ============================================================
def _annulus_stop_time_s(casing_result, fluids: tuple[FluidSpec, ...]) -> float:
    """返回呼1-004环空二维顶替应停止的地面累计时间。

    对齐 cemdisp.runners.ht1_004_tailpipe.annulus_stop_time_s：
    当水泥浆之后的首个非水泥流体到达鞋口时停止；若未找到，
    回退到 cement_end_time_s。
    """
    cement_roles = {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}
    fluid_by_name = {fluid.name: fluid for fluid in fluids}
    found_cement = False
    for front in casing_result.fronts:
        fluid = fluid_by_name.get(front.fluid_name)
        if fluid is None:
            continue
        if fluid.role in cement_roles:
            found_cement = True
            continue
        if found_cement:
            return float(front.time_s)
    return float(casing_result.cement_end_time_s)


def _run_single(
    well_spec,
    fluids: tuple[FluidSpec, ...],
    schedule: PumpingSchedule,
) -> tuple[float, float]:
    """运行一次完整 1D-2D 耦合模拟，返回 (全井有效顶替效率, 耗时秒数)。

    参数：
        well_spec: 井身结构规格
        fluids: 流体元组（已修改）
        schedule: 泵送计划（已修改）

    返回：
        (effective_efficiency, elapsed_s)
    """
    t0 = time.perf_counter()

    # 1D 套管内流动
    casing_solver = CasingFlowSolver(enable_gravity=True)
    casing_result = casing_solver.run(well_spec, fluids, schedule)

    # 鞋口出流 → 环空入口桥接
    coupled_provider = build_coupled_annulus_inlet_provider(
        casing_result, casing_solver, fluids, split_cement_phases=False,
    )

    # 环空二维顶替：停止时间对齐呼1-004 runner 的严格现场口径，
    # 即水泥之后的首个非水泥流体到达鞋口时停止；批量敏感性仍用 nz=140 加速。
    total_t_s = _annulus_stop_time_s(casing_result, fluids)
    annulus_solver = AnnulusD2DGASolver(total_t=total_t_s, nz=140)
    result = annulus_solver.run(well_spec, fluids, coupled_provider)

    # 提取全井有效顶替效率
    final_metrics = result.metrics.iloc[-1]
    efficiency = float(final_metrics["effective_efficiency"])
    elapsed_s = time.perf_counter() - t0

    return efficiency, elapsed_s


# ============================================================
# 敏感性分析主循环（含断点续跑）
# ============================================================
def run_sensitivity_analysis(resume: bool = False, limit: int = 0) -> pd.DataFrame:
    """执行全因子敏感性分析，返回结果 DataFrame。

    参数：
        resume: 是否从中间结果续跑
        limit: 限制运行组数，0=不限制（运行全部）

    返回：
        包含所有参数组合和对应效率的 DataFrame
    """
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # 加载基线数据
    well_spec, base_fluids, base_schedule, _ = load_ht1_004_tailpipe()

    # 构建全因子参数网格
    all_combos = []
    for rate in RATE_VALUES:
        for yp in YP_VALUES:
            for pv in PV_VALUES:
                all_combos.append((float(rate), float(yp), float(pv)))

    total_runs = len(all_combos)
    print(f"全因子组合数: {total_runs}")
    print(f"  排量: {len(RATE_VALUES)} 个水平 ({RATE_MIN}~{RATE_MAX}, 步长{RATE_STEP})")
    print(f"  屈服应力: {len(YP_VALUES)} 个水平 ({YP_MIN}~{YP_MAX}, 步长{YP_STEP})")
    print(f"  塑性黏度: {len(PV_VALUES)} 个水平 ({PV_MIN}~{PV_MAX}, 步长{PV_STEP})")

    # 断点续跑：加载已完成的组合
    completed_set: set[tuple[float, float, float]] = set()
    existing_rows: list[dict] = []

    if resume and CHECKPOINT_PATH.exists():
        df_existing = pd.read_csv(CHECKPOINT_PATH)
        for _, row in df_existing.iterrows():
            key = (round(row["排量_m3min"], 2), round(row["屈服应力_Pa"], 1), round(row["黏度_mPas"], 0))
            completed_set.add(key)
            existing_rows.append(row.to_dict())
        print(f"断点续跑: 已加载 {len(completed_set)} 个已完成组合")

    # 筛选待运行的组合
    remaining = [(r, y, p) for r, y, p in all_combos
                 if (round(r, 2), round(y, 1), round(p, 0)) not in completed_set]

    # 限制运行组数时，均匀采样覆盖整个参数空间（而非只取前N个）
    if limit > 0 and len(remaining) > limit:
        step = len(remaining) / limit
        remaining = [remaining[int(i * step)] for i in range(limit)]
        # 更新进度条的总数为实际要运行的数量
        total_runs = len(existing_rows) + len(remaining)
        print(f"限制模式: 均匀采样 {limit} 组覆盖整个参数空间")

    if not remaining:
        print("所有组合已完成，直接生成图表。")
        return pd.DataFrame(existing_rows)

    print(f"待运行: {len(remaining)} / {total_runs}")
    print("=" * 60)

    # 运行模拟（带进度条）
    results: list[dict] = list(existing_rows)
    completed_before = len(existing_rows)

    # 进度条：显示已完成/总数、效率、速度、预计剩余时间
    pbar = tqdm(
        remaining,
        desc="敏感性分析",
        unit="次",
        initial=completed_before,
        total=total_runs,
        ncols=100,
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}] {postfix}",
    )

    for run_idx, (rate, yp, pv) in enumerate(pbar):
        # 将 mPa·s 转换为 Pa·s
        pv_pa_s = pv / 1000.0

        # 创建修改后的数据
        mod_fluids, mod_schedule = _create_modified_data(
            base_fluids, base_schedule,
            tail_yp_pa=yp, tail_pv_pa_s=pv_pa_s, tail_rate_m3_min=rate,
        )

        # 运行模拟
        efficiency, elapsed_s = _run_single(well_spec, mod_fluids, mod_schedule)

        row = {
            "排量_m3min": round(rate, 2),
            "屈服应力_Pa": round(yp, 1),
            "黏度_mPas": round(pv, 0),
            "全井有效顶替效率": round(efficiency, 6),
            "耗时_s": round(elapsed_s, 1),
        }
        results.append(row)

        # 更新进度条后缀：显示当前参数和效率
        pbar.set_postfix_str(
            f"效率={efficiency:.4f} | {rate:.2f}m^3/min, {yp:.0f}Pa, {pv:.0f}mPa·s"
        )

        # 中间保存（每100次）
        global_idx = len(results)
        if global_idx % 100 == 0:
            df_checkpoint = pd.DataFrame(results)
            df_checkpoint.to_csv(CHECKPOINT_PATH, index=False, encoding="utf-8-sig")

    pbar.close()

    # 最终保存完整CSV
    df_result = pd.DataFrame(results)
    df_result.to_csv(CSV_PATH, index=False, encoding="utf-8-sig")
    print(f"\n完整结果已保存: {CSV_PATH}")

    # 清理中间文件
    if CHECKPOINT_PATH.exists():
        CHECKPOINT_PATH.unlink()

    return df_result


# ============================================================
# 图表生成
# ============================================================
def _extract_subset(
    df: pd.DataFrame,
    fixed_params: dict[str, float],
    x_col: str,
) -> pd.DataFrame:
    """从全量结果中提取固定某些参数后的子集。

    参数：
        df: 全量结果 DataFrame
        fixed_params: 固定参数 {列名: 值}
        x_col: X轴列名（不在fixed_params中）

    返回：
        过滤后的 DataFrame，按 x_col 排序
    """
    mask = pd.Series(True, index=df.index)
    for col, val in fixed_params.items():
        mask = mask & (df[col].round(4) == round(val, 4))
    subset = df[mask].copy()
    return subset.sort_values(x_col)


def plot_single_factor_curves(df: pd.DataFrame) -> None:
    """生成3张单因素敏感性曲线图。

    每个图横轴为扫描参数，纵轴为全井有效顶替效率的均值±标准差。
    对横轴每个取值，汇总该参数的所有组合（其他参数取不同值）计算统计量，
    从而在数据稀疏时也能给出有意义的趋势曲线。
    """
    _setup_chinese_font()

    configs = [
        {
            "x_col": "排量_m3min",
            "x_label": "尾浆排量 (m^3/min)",
            "title": "排量对全井有效顶替效率的影响",
            "filename": "呼1-004_排量敏感性.png",
        },
        {
            "x_col": "屈服应力_Pa",
            "x_label": "尾浆屈服应力 (Pa)",
            "title": "屈服应力对全井有效顶替效率的影响",
            "filename": "呼1-004_屈服应力敏感性.png",
        },
        {
            "x_col": "黏度_mPas",
            "x_label": "尾浆塑性黏度 (mPa·s)",
            "title": "塑性黏度对全井有效顶替效率的影响",
            "filename": "呼1-004_黏度敏感性.png",
        },
    ]

    for cfg in configs:
        # 按横轴参数分组，计算均值和标准差（汇总所有其他参数的组合）
        grouped = df.groupby(cfg["x_col"])["全井有效顶替效率"]
        means = grouped.mean()
        stds = grouped.std()

        if means.empty:
            print(f"警告: {cfg['filename']} 数据为空，跳过")
            continue

        fig, ax = plt.subplots(figsize=(8, 5))
        x_vals = means.index.values
        y_vals = means.values

        # 均值曲线 + 误差带（±1σ）
        ax.plot(x_vals, y_vals, "o-", color="#2196F3", linewidth=2, markersize=5, label="均值")
        ax.fill_between(
            x_vals,
            y_vals - stds.values,
            y_vals + stds.values,
            alpha=0.2, color="#2196F3", label="±1 标准差",
        )

        ax.set_xlabel(cfg["x_label"], fontsize=12)
        ax.set_ylabel("全井有效顶替效率", fontsize=12)
        ax.set_title(cfg["title"], fontsize=14)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 1.05)

        fig.tight_layout()
        save_path = OUTPUT_DIR / cfg["filename"]
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"已保存: {save_path}")


def _format_fixed_param(col: str, val: float) -> str:
    if col == "排量_m3min":
        return f"排量={val:.2f} m^3/min"
    if col == "屈服应力_Pa":
        return f"屈服应力={val:.0f} Pa"
    if col == "黏度_mPas":
        return f"黏度={val:.0f} mPa·s"
    return f"{col}={val:g}"


def _format_axis_tick(col: str, val: float) -> str:
    if col == "排量_m3min":
        return f"{val:.2f}"
    return f"{val:.0f}"


def plot_heatmaps(df: pd.DataFrame) -> None:
    """生成3张两两交互热力图。"""
    _setup_chinese_font()

    configs = [
        {
            "x_col": "排量_m3min",
            "y_col": "屈服应力_Pa",
            "x_label": "尾浆排量 (m^3/min)",
            "y_label": "尾浆屈服应力 (Pa)",
            "fixed": {"黏度_mPas": PV_BASELINE},
            "title": f"排量-屈服应力 交互效应 (黏度={PV_BASELINE:.0f} mPa·s)",
            "filename": "呼1-004_排量-屈服应力_热力图.png",
        },
        {
            "x_col": "排量_m3min",
            "y_col": "黏度_mPas",
            "x_label": "尾浆排量 (m^3/min)",
            "y_label": "尾浆塑性黏度 (mPa·s)",
            "fixed": {"屈服应力_Pa": YP_BASELINE},
            "title": f"排量-黏度 交互效应 (屈服应力={YP_BASELINE:.0f} Pa)",
            "filename": "呼1-004_排量-黏度_热力图.png",
        },
        {
            "x_col": "屈服应力_Pa",
            "y_col": "黏度_mPas",
            "x_label": "尾浆屈服应力 (Pa)",
            "y_label": "尾浆塑性黏度 (mPa·s)",
            "fixed": {"排量_m3min": RATE_BASELINE},
            "title": f"屈服应力-黏度 交互效应 (排量={RATE_BASELINE:.1f} m^3/min)",
            "filename": "呼1-004_屈服应力-黏度_热力图.png",
        },
    ]

    for cfg in configs:
        # 过滤固定参数；若精确基线值不存在则取最近邻
        subset = df.copy()
        actual_fixed = {}
        for col, val in cfg["fixed"].items():
            unique_vals = sorted(df[col].unique())
            if val in unique_vals:
                used_val = val
            else:
                used_val = min(unique_vals, key=lambda v: abs(v - val))
                print(f"  注意: {cfg['filename']} 固定参数 {col} 基线={val} 不存在，使用最近值={used_val}")
            actual_fixed[col] = used_val
            subset = subset[subset[col].round(4) == round(used_val, 4)]
        # 更新标题中的固定参数值为实际使用的值
        title_fixed_str = ", ".join(
            _format_fixed_param(col, actual_fixed[col])
            for col in cfg["fixed"]
        )
        cfg_title = cfg["title"].split("(")[0].strip() + f" ({title_fixed_str})"

        if subset.empty:
            print(f"警告: {cfg['filename']} 数据为空，跳过")
            continue

        # 构建透视表
        pivot = subset.pivot_table(
            index=cfg["y_col"],
            columns=cfg["x_col"],
            values="全井有效顶替效率",
            aggfunc="mean",
        )

        fig, ax = plt.subplots(figsize=(10, 7))
        im = ax.imshow(
            pivot.values,
            aspect="auto",
            cmap="RdYlGn",
            vmin=0,
            vmax=1,
            origin="lower",
        )

        # 设置刻度
        ax.set_xticks(range(len(pivot.columns)))
        ax.set_xticklabels([_format_axis_tick(cfg["x_col"], v) for v in pivot.columns], fontsize=9)
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels([_format_axis_tick(cfg["y_col"], v) for v in pivot.index], fontsize=9)

        ax.set_xlabel(cfg["x_label"], fontsize=12)
        ax.set_ylabel(cfg["y_label"], fontsize=12)
        ax.set_title(cfg_title, fontsize=14)

        # 颜色条
        cbar = fig.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label("全井有效顶替效率", fontsize=11)

        # 在每个格子中标注数值
        for i in range(len(pivot.index)):
            for j in range(len(pivot.columns)):
                val = pivot.values[i, j]
                if not np.isnan(val):
                    # 深色格子用白色字，浅色格子用黑色字
                    text_color = "white" if val < 0.3 or val > 0.85 else "black"
                    ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                            fontsize=7, color=text_color)

        fig.tight_layout()
        save_path = OUTPUT_DIR / cfg["filename"]
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"已保存: {save_path}")


# ============================================================
# 主入口
# ============================================================
def main() -> None:
    parser = argparse.ArgumentParser(description="呼1-004井尾浆参数敏感性分析")
    parser.add_argument(
        "--resume", action="store_true",
        help="断点续跑：从中间结果CSV继续",
    )
    parser.add_argument(
        "--plot-only", action="store_true",
        help="仅生成图表（需要已有完整CSV结果）",
    )
    parser.add_argument(
        "--limit", type=int, default=0,
        help="限制运行组数（调试用），0=不限制",
    )
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if args.plot_only:
        # 仅绘图模式
        if not CSV_PATH.exists():
            print(f"错误: 未找到结果文件 {CSV_PATH}")
            sys.exit(1)
        df = pd.read_csv(CSV_PATH)
        print(f"已加载 {len(df)} 行结果")
    else:
        # 运行敏感性分析
        print("=" * 60)
        print("呼1-004井尾浆参数敏感性分析")
        print("=" * 60)
        t_start = time.perf_counter()
        df = run_sensitivity_analysis(resume=args.resume, limit=args.limit)
        t_total = time.perf_counter() - t_start
        print(f"\n总耗时: {t_total / 3600:.1f} 小时")

    # 生成图表
    print("\n生成图表...")
    plot_single_factor_curves(df)
    plot_heatmaps(df)

    # 打印统计摘要
    print("\n" + "=" * 60)
    print("统计摘要")
    print("=" * 60)
    print(f"总组合数: {len(df)}")
    print(f"效率范围: {df['全井有效顶替效率'].min():.4f} ~ {df['全井有效顶替效率'].max():.4f}")
    print(f"效率均值: {df['全井有效顶替效率'].mean():.4f}")
    print(f"效率标准差: {df['全井有效顶替效率'].std():.4f}")

    # 基线效率
    baseline_row = df[
        (df["排量_m3min"].round(2) == RATE_BASELINE) &
        (df["屈服应力_Pa"].round(1) == YP_BASELINE) &
        (df["黏度_mPas"].round(0) == PV_BASELINE)
    ]
    if not baseline_row.empty:
        baseline_eff = baseline_row["全井有效顶替效率"].iloc[0]
        print(f"基线效率 (排量={RATE_BASELINE}, YP={YP_BASELINE}, PV={PV_BASELINE}): {baseline_eff:.4f}")

    print(f"\n结果已保存到: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
