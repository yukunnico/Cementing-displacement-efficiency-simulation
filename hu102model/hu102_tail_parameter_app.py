"""
呼102尾管 D2DGA 参数输入页面
=============================

提供一个本地 Tkinter 图形界面，用于输入屈服值、塑性黏度和排量，
并调用三因素敏感性模型计算不同影响因素组合下的全井段有效顶替效率。

运行方式：
    python hu102_tail_parameter_app.py
"""

from __future__ import annotations

from collections.abc import Mapping
import importlib
import json
from numbers import Real
import queue
import threading
from pathlib import Path
import tkinter as tk
from tkinter import messagebox, ttk

import pandas as pd

# 按模块名动态导入三因素敏感性脚本，便于与主模型脚本解耦。
sensitivity = importlib.import_module("hu102_tail_three_factor_sensitivity")


# 统一用户输入场景的输出文件位置，确保每次运行产物路径固定。
OUT_DIR = Path(__file__).resolve().parent
USER_SCENARIO_CSV = OUT_DIR / "呼102尾管_D2DGA_用户输入参数_场景结果.csv"
USER_SUMMARY_JSON = OUT_DIR / "呼102尾管_D2DGA_用户输入参数_结果摘要.json"
USER_SUMMARY_MD = OUT_DIR / "呼102尾管_D2DGA_用户输入参数_结果摘要.md"
MAX_SCENARIO_COUNT = 64


def require_mapping(value: object, name: str) -> Mapping[str, object]:
    """确保摘要字段为映射类型。"""
    if not isinstance(value, Mapping):
        raise TypeError(f"{name}格式错误：需要字典类型")
    return value


def require_ranking(value: object) -> list[Mapping[str, object]]:
    """确保影响排序字段为字典列表。"""
    if not isinstance(value, list):
        raise TypeError("impact_ranking格式错误：需要列表类型")

    rows: list[Mapping[str, object]] = []
    for item in value:
        if not isinstance(item, Mapping):
            raise TypeError("impact_ranking格式错误：列表元素需要字典类型")
        rows.append(item)
    return rows


def number_field(mapping: Mapping[str, object], key: str) -> float:
    """从摘要映射中读取数值字段。"""
    value = mapping[key]
    if isinstance(value, Real | str):
        return float(value)
    raise TypeError(f"{key}格式错误：需要数值")


def parse_number_list(text: str, name: str) -> list[float]:
    """
    解析用户输入的数值列表。

    支持英文逗号、中文逗号、分号、空格和换行分隔，例如：
    ``1.0, 1.3, 1.6`` 或 ``60 80 100``。
    """
    # 将常见分隔符统一成英文逗号，后续只需一套拆分逻辑。
    normalized = text.replace("，", ",").replace("；", ",").replace(";", ",")
    normalized = normalized.replace("\n", ",").replace("\t", ",").replace(" ", ",")
    parts = [part.strip() for part in normalized.split(",") if part.strip()]
    if not parts:
        raise ValueError(f"{name}至少需要输入一个数值")

    values: list[float] = []
    for part in parts:
        try:
            values.append(float(part))
        except ValueError as exc:
            raise ValueError(f"{name}包含无法识别的数值：{part}") from exc
    return values


def save_user_outputs(result_df: pd.DataFrame, summary: dict[str, object]) -> dict[str, str]:
    """保存用户输入场景的CSV、JSON和Markdown摘要。"""
    result_df.to_csv(USER_SCENARIO_CSV, index=False, encoding="utf-8-sig")

    outputs = {
        "scenario_csv": str(USER_SCENARIO_CSV),
        "summary_json": str(USER_SUMMARY_JSON),
        "summary_md": str(USER_SUMMARY_MD),
    }
    # 在摘要中写入本次分析名称与输出文件清单，便于归档追溯。
    summary_with_outputs = dict(summary)
    summary_with_outputs["analysis_name"] = "呼102尾管_D2DGA 用户输入参数全井段顶替效率计算"
    summary_with_outputs["outputs"] = outputs
    USER_SUMMARY_JSON.write_text(json.dumps(summary_with_outputs, ensure_ascii=False, indent=2), encoding="utf-8")

    best_case = require_mapping(summary["best_case"], "best_case")
    recommended_ranges = require_mapping(summary["recommended_ranges"], "recommended_ranges")
    rate_range = require_mapping(recommended_ranges["rate_m3_min"], "rate_m3_min")
    pv_range = require_mapping(recommended_ranges["pv_mpa_s"], "pv_mpa_s")
    yp_range = require_mapping(recommended_ranges["yp_pa"], "yp_pa")
    impact_ranking = require_ranking(summary["impact_ranking"])
    # 影响排序条目单独格式化，便于在 Markdown 中直接展示。
    rank_lines = [
        f"- {int(number_field(row, 'rank'))}. {row['factor_name']}：主效应跨度={number_field(row, 'mean_effect_span'):.5f}，贡献占比={number_field(row, 'importance_share'):.2%}"
        + (f"；{row['note']}" if row.get("note") else "")
        for row in impact_ranking
    ]
    selection_rule = require_mapping(recommended_ranges["selection_rule"], "selection_rule")
    top_n = int(number_field(selection_rule, "top_n"))
    # 生成可读性更好的汇报模板，供现场参数讨论直接引用。
    md_lines = [
        "# 呼102尾管 D2DGA 用户输入参数结果摘要",
        "",
        f"- 场景数量：{summary['scenario_count']}",
        "- 目标指标：effective_efficiency（全井段有效顶替效率）",
        "",
        "## 最优场景（按全井段有效顶替效率）",
        f"- 排量：{number_field(best_case, 'rate_m3_min'):.3f} m3/min",
        f"- 塑性黏度：{number_field(best_case, 'pv_mpa_s'):.3f} mPa·s",
        f"- 屈服值：{number_field(best_case, 'yp_pa'):.3f} Pa",
        f"- 全井段有效顶替效率：{number_field(best_case, 'effective_efficiency'):.4f}",
        f"- CBL评价井段水动力效率：{number_field(best_case, 'cbl_eval_interval_efficiency'):.4f}",
        f"- CBL评价井段质量响应效率：{number_field(best_case, 'cbl_quality_proxy'):.4f}",
        "",
        "## 影响因素排序",
        *rank_lines,
        "",
        f"## 建议参数区间（目标指标前20%场景，实际选取{top_n}个场景）",
        f"- 排量：{number_field(rate_range, 'min'):.3f} - {number_field(rate_range, 'max'):.3f} m3/min",
        f"- 塑性黏度：{number_field(pv_range, 'min'):.3f} - {number_field(pv_range, 'max'):.3f} mPa·s",
        f"- 屈服值：{number_field(yp_range, 'min'):.3f} - {number_field(yp_range, 'max'):.3f} Pa",
        "",
        "## 输出文件",
        f"- 场景结果CSV：`{USER_SCENARIO_CSV}`",
        f"- 结果摘要JSON：`{USER_SUMMARY_JSON}`",
        f"- 结果摘要Markdown：`{USER_SUMMARY_MD}`",
    ]
    USER_SUMMARY_MD.write_text("\n".join(md_lines), encoding="utf-8")
    return outputs


class ParameterApp:
    """参数输入页面。"""

    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("呼102尾管 D2DGA 顶替效率参数计算")
        self.root.geometry("1080x720")
        self.result_queue: queue.Queue[tuple[str, object]] = queue.Queue()

        self.rate_var = tk.StringVar(value="1.00, 1.30, 1.60")
        self.pv_var = tk.StringVar(value="60, 80, 100")
        self.yp_var = tk.StringVar(value="15, 22, 30")
        self.status_var = tk.StringVar(value="请输入参数后点击“开始计算”。多个值会自动组合成不同影响因素场景。")
        self.progress_var = tk.DoubleVar(value=0.0)

        self._build_layout()

    def _build_layout(self) -> None:
        # 主容器：包含输入区、进度区、结果表格和摘要文本。
        container = ttk.Frame(self.root, padding=16)
        container.pack(fill=tk.BOTH, expand=True)

        title = ttk.Label(container, text="呼102尾管 D2DGA 全井段顶替效率计算", font=("Microsoft YaHei", 16, "bold"))
        title.pack(anchor=tk.W, pady=(0, 10))

        # 输入区支持多值，程序会做笛卡尔组合形成多场景计算。
        input_frame = ttk.LabelFrame(container, text="输入参数（可输入单值或多个值）", padding=12)
        input_frame.pack(fill=tk.X, pady=(0, 10))

        self._add_input_row(input_frame, 0, "排量 / m³·min⁻¹", self.rate_var, "例如：1.0, 1.3, 1.6")
        self._add_input_row(input_frame, 1, "塑性黏度 PV / mPa·s", self.pv_var, "例如：60, 80, 100")
        self._add_input_row(input_frame, 2, "屈服值 YP / Pa", self.yp_var, "例如：15, 22, 30")

        button_frame = ttk.Frame(input_frame)
        button_frame.grid(row=3, column=1, sticky=tk.W, pady=(10, 0))
        self.run_button = ttk.Button(button_frame, text="开始计算", command=self.start_calculation)
        self.run_button.pack(side=tk.LEFT)
        ttk.Button(button_frame, text="恢复默认三水平", command=self.reset_defaults).pack(side=tk.LEFT, padx=(8, 0))

        ttk.Label(container, textvariable=self.status_var, foreground="#1D4ED8").pack(anchor=tk.W, pady=(0, 8))

        self.progress_bar = ttk.Progressbar(container, variable=self.progress_var, maximum=100.0)
        self.progress_bar.pack(fill=tk.X, pady=(0, 8))

        # 表格列与输出 CSV 保持一致，方便人工核对与二次分析。
        columns = (
            "rate_m3_min",
            "pv_mpa_s",
            "yp_pa",
            "effective_efficiency",
            "cbl_eval_interval_efficiency",
            "cbl_quality_proxy",
            "target_interval_efficiency",
            "channeling_index",
        )
        self.table = ttk.Treeview(container, columns=columns, show="headings", height=14)
        headers = {
            "rate_m3_min": "排量",
            "pv_mpa_s": "PV",
            "yp_pa": "YP",
            "effective_efficiency": "全井段效率",
            "cbl_eval_interval_efficiency": "CBL井段效率",
            "cbl_quality_proxy": "质量响应效率",
            "target_interval_efficiency": "油气水层效率",
            "channeling_index": "窜槽指数",
        }
        for column in columns:
            self.table.heading(column, text=headers[column])
            self.table.column(column, width=120, anchor=tk.CENTER)
        self.table.pack(fill=tk.BOTH, expand=True)

        self.summary_text = tk.Text(container, height=10, wrap=tk.WORD)
        self.summary_text.pack(fill=tk.X, pady=(10, 0))
        self.summary_text.insert(tk.END, "结果摘要会显示在这里。\n")
        self.summary_text.configure(state=tk.DISABLED)

    def _add_input_row(self, parent: tk.Misc, row: int, label: str, variable: tk.StringVar, hint: str) -> None:
        ttk.Label(parent, text=label, width=22).grid(row=row, column=0, sticky=tk.W, pady=5)
        entry = ttk.Entry(parent, textvariable=variable, width=42)
        entry.grid(row=row, column=1, sticky=tk.W, pady=5)
        ttk.Label(parent, text=hint, foreground="#64748B").grid(row=row, column=2, sticky=tk.W, padx=(10, 0), pady=5)

    def reset_defaults(self) -> None:
        """恢复脚本内置的三因素三水平默认值。"""
        self.rate_var.set("1.00, 1.30, 1.60")
        self.pv_var.set("60, 80, 100")
        self.yp_var.set("15, 22, 30")

    def start_calculation(self) -> None:
        """读取输入并启动后台计算。"""
        try:
            rates = parse_number_list(self.rate_var.get(), "排量")
            pvs = parse_number_list(self.pv_var.get(), "塑性黏度")
            yps = parse_number_list(self.yp_var.get(), "屈服值")
            # 使用去重后的水平数计算总场景数，避免重复输入放大计算量。
            scenario_count = len(set(rates)) * len(set(pvs)) * len(set(yps))
            if scenario_count > MAX_SCENARIO_COUNT:
                raise ValueError(f"当前会生成{scenario_count}个场景，超过上限{MAX_SCENARIO_COUNT}；请减少输入水平数量。")
        except ValueError as exc:
            messagebox.showerror("输入错误", str(exc))
            return

        self.run_button.configure(state=tk.DISABLED)
        self.progress_var.set(0.0)
        self.status_var.set(f"正在计算 {scenario_count} 个场景，请稍候；模型会逐场景运行完整D2DGA模拟。")
        self._clear_results()

        # 数值模拟放到后台线程，避免 Tk 主线程卡顿。
        worker = threading.Thread(target=self._calculate_in_worker, args=(rates, pvs, yps), daemon=True)
        worker.start()
        self.root.after(250, self._poll_result_queue)

    def _calculate_in_worker(self, rates: list[float], pvs: list[float], yps: list[float]) -> None:
        """后台线程执行耗时的数值模拟。"""
        try:
            result_df = sensitivity.run_factorial_cases(
                rates,
                pvs,
                yps,
                include_baseline=False,
                scenario_type="user_input_factorial",
                progress_callback=self._report_progress,
            )
            summary = sensitivity.summarize_factorial_results(
                result_df,
                target_metric=sensitivity.FULL_WELL_TARGET_METRIC,
            )
            outputs = save_user_outputs(result_df, summary)
            self.result_queue.put(("success", {"result_df": result_df, "summary": summary, "outputs": outputs}))
        except Exception as exc:
            self.result_queue.put(("error", str(exc)))

    def _report_progress(self, completed: int, total: int, rate: float, pv: float, yp: float) -> None:
        """从后台线程回传当前计算进度。"""
        self.result_queue.put(("progress", {
            "completed": completed,
            "total": total,
            "rate": rate,
            "pv": pv,
            "yp": yp,
        }))

    def _poll_result_queue(self) -> None:
        """轮询后台线程结果并更新界面。"""
        try:
            status, payload = self.result_queue.get_nowait()
        except queue.Empty:
            # 队列暂无数据时继续轮询，保持界面响应。
            self.root.after(250, self._poll_result_queue)
            return

        self.run_button.configure(state=tk.NORMAL)
        if status == "progress":
            self.run_button.configure(state=tk.DISABLED)
            if isinstance(payload, dict):
                completed = int(number_field(payload, "completed"))
                total = int(number_field(payload, "total"))
                rate = number_field(payload, "rate")
                pv = number_field(payload, "pv")
                yp = number_field(payload, "yp")
                self.progress_var.set(completed / max(total, 1) * 100.0)
                self.status_var.set(
                    f"已完成 {completed}/{total} 个场景：排量={rate:.3f} m3/min，PV={pv:.3f} mPa·s，YP={yp:.3f} Pa"
                )
            # 继续轮询后续进度或最终结果。
            self.root.after(250, self._poll_result_queue)
            return

        self.run_button.configure(state=tk.NORMAL)
        if status == "error":
            self.status_var.set("计算失败，请检查输入或模型文件。")
            messagebox.showerror("计算失败", str(payload))
            return

        if isinstance(payload, dict):
            result_df = payload["result_df"]
            summary = payload["summary"]
            outputs = payload["outputs"]
            if isinstance(result_df, pd.DataFrame) and isinstance(summary, dict) and isinstance(outputs, dict):
                self.progress_var.set(100.0)
                self._show_results(result_df, summary, outputs)

    def _clear_results(self) -> None:
        for item in self.table.get_children():
            self.table.delete(item)
        self.summary_text.configure(state=tk.NORMAL)
        self.summary_text.delete("1.0", tk.END)
        self.summary_text.configure(state=tk.DISABLED)

    def _show_results(self, result_df: pd.DataFrame, summary: dict[str, object], outputs: dict[str, str]) -> None:
        """将计算结果展示到表格和摘要框。"""
        # 优先展示最优场景，便于第一时间定位推荐参数组合。
        display_df = result_df.sort_values("effective_efficiency", ascending=False).reset_index(drop=True)
        columns = list(self.table["columns"])
        for _, row in display_df.iterrows():
            # 指标列保留4位小数，输入参数列保留3位小数。
            values = [f"{float(row[column]):.4f}" if column.endswith("efficiency") or column.endswith("proxy") or column.endswith("index") else f"{float(row[column]):.3f}" for column in columns]
            self.table.insert("", tk.END, values=values)

        best_case = require_mapping(summary["best_case"], "best_case")
        recommended_ranges = require_mapping(summary["recommended_ranges"], "recommended_ranges")
        rate_range = require_mapping(recommended_ranges["rate_m3_min"], "rate_m3_min")
        pv_range = require_mapping(recommended_ranges["pv_mpa_s"], "pv_mpa_s")
        yp_range = require_mapping(recommended_ranges["yp_pa"], "yp_pa")
        selection_rule = require_mapping(recommended_ranges["selection_rule"], "selection_rule")
        top_n = int(number_field(selection_rule, "top_n"))
        lines = [
            "计算完成。",
            f"场景数量：{summary['scenario_count']}",
            f"最优排量：{number_field(best_case, 'rate_m3_min'):.3f} m3/min",
            f"最优塑性黏度：{number_field(best_case, 'pv_mpa_s'):.3f} mPa·s",
            f"最优屈服值：{number_field(best_case, 'yp_pa'):.3f} Pa",
            f"最高全井段有效顶替效率：{number_field(best_case, 'effective_efficiency'):.4f}",
            f"建议区间依据：目标指标前20%场景，实际选取{top_n}个场景",
            f"建议排量区间：{number_field(rate_range, 'min'):.3f} - {number_field(rate_range, 'max'):.3f} m3/min",
            f"建议PV区间：{number_field(pv_range, 'min'):.3f} - {number_field(pv_range, 'max'):.3f} mPa·s",
            f"建议YP区间：{number_field(yp_range, 'min'):.3f} - {number_field(yp_range, 'max'):.3f} Pa",
            "",
            "输出文件：",
            f"- {outputs['scenario_csv']}",
            f"- {outputs['summary_json']}",
            f"- {outputs['summary_md']}",
        ]
        self.summary_text.configure(state=tk.NORMAL)
        self.summary_text.delete("1.0", tk.END)
        self.summary_text.insert(tk.END, "\n".join(lines))
        self.summary_text.configure(state=tk.DISABLED)
        self.status_var.set("计算完成：已按全井段有效顶替效率排序，并保存CSV/JSON/Markdown结果。")


def main() -> None:
    """启动参数输入页面。"""
    root = tk.Tk()
    app = ParameterApp(root)
    app.root.mainloop()


if __name__ == "__main__":
    main()
