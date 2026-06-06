# -*- coding: utf-8 -*-
"""cemdisp模型全井测试脚本 - 测试所有6口井的尾管段模型能否正常运行并生成结果汇总表."""

from __future__ import annotations

import json
import sys
import time
import traceback
from collections import OrderedDict
from dataclasses import dataclass, field
from io import StringIO
from pathlib import Path
from typing import Any

# 确保项目根目录在 sys.path 中
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


@dataclass
class WellResult:
    """单口井的运行结果记录。"""
    well_name: str
    runner_name: str
    success: bool = False
    error_message: str = ""
    run_time_s: float = 0.0
    summary: dict[str, Any] = field(default_factory=dict)
    final_result: dict[str, float] = field(default_factory=dict)


def run_well(
    well_name: str,
    runner_name: str,
    runner_func,
) -> WellResult:
    """运行单口井的模型并捕获结果。"""
    result = WellResult(well_name=well_name, runner_name=runner_name)
    print(f"\n{'='*60}")
    print(f"开始运行: {well_name} ({runner_name})")
    print(f"{'='*60}")

    start_time = time.time()
    try:
        # 捕获 stdout 以获取 summary 输出
        old_stdout = sys.stdout
        sys.stdout = captured = StringIO()
        try:
            runner_func()
        finally:
            sys.stdout = old_stdout

        result.run_time_s = time.time() - start_time
        result.success = True

        # 从捕获的输出中解析 JSON 摘要
        output = captured.getvalue()
        # 查找最后一个 JSON 块（在 === 标记之后）
        json_blocks = []
        current_block = []
        in_json = False
        for line in output.split('\n'):
            if line.strip().startswith('{'):
                in_json = True
                current_block = [line]
            elif in_json:
                current_block.append(line)
                if line.strip().endswith('}'):
                    try:
                        block_text = '\n'.join(current_block)
                        parsed = json.loads(block_text)
                        json_blocks.append(parsed)
                    except json.JSONDecodeError:
                        pass
                if line.strip() == '}':
                    in_json = False
                    current_block = []

        if json_blocks:
            result.summary = json_blocks[-1]
            final = result.summary.get("最终结果", {})
            result.final_result = final

        print(f"[成功] {well_name} 运行完成，耗时 {result.run_time_s:.1f} 秒")

    except Exception as e:
        result.run_time_s = time.time() - start_time
        result.error_message = f"{type(e).__name__}: {e}"
        print(f"[失败] {well_name} 运行失败: {result.error_message}")
        traceback.print_exc()

    return result


def format_number(value: float, decimals: int = 4) -> str:
    """格式化数字显示。"""
    if value == 0.0 and decimals == 0:
        return "N/A"
    return f"{value:.{decimals}f}"


def print_summary_table(results: list[WellResult]) -> None:
    """打印结果汇总表。"""
    print("\n")
    print("=" * 120)
    print("cemdisp 模型全井测试结果汇总")
    print("=" * 120)

    # 表头
    header = f"{'井号':<12} {'状态':<8} {'耗时(s)':<10} {'有效顶替效率':<14} {'水泥浆占据率':<14} {'窜槽指数':<12} {'混浆指数':<12} {'失稳指数':<12}"
    print(header)
    print("-" * 120)

    success_count = 0
    fail_count = 0

    for r in results:
        status = "成功" if r.success else "失败"
        if r.success:
            success_count += 1
        else:
            fail_count += 1

        if r.success and r.final_result:
            fr = r.final_result
            row = (
                f"{r.well_name:<12} "
                f"{status:<8} "
                f"{r.run_time_s:<10.1f} "
                f"{format_number(fr.get('全井段最终有效顶替效率', 0)):<14} "
                f"{format_number(fr.get('最终水泥浆占据率', 0)):<14} "
                f"{format_number(fr.get('最终窜槽指数', 0)):<12} "
                f"{format_number(fr.get('最终混浆指数', 0)):<12} "
                f"{format_number(fr.get('最终失稳指数', 0)):<12}"
            )
        else:
            row = (
                f"{r.well_name:<12} "
                f"{status:<8} "
                f"{r.run_time_s:<10.1f} "
                f"{'N/A':<14} "
                f"{'N/A':<14} "
                f"{'N/A':<12} "
                f"{'N/A':<12} "
                f"{'N/A':<12}"
            )
        print(row)

    print("-" * 120)
    print(f"总计: {len(results)} 口井，成功 {success_count} 口，失败 {fail_count} 口")
    print("=" * 120)


def print_detailed_results(results: list[WellResult]) -> None:
    """打印每口井的详细结果。"""
    print("\n")
    print("=" * 80)
    print("详细结果")
    print("=" * 80)

    for r in results:
        print(f"\n--- {r.well_name} ---")
        print(f"运行器: {r.runner_name}")
        print(f"状态: {'成功' if r.success else '失败'}")
        print(f"耗时: {r.run_time_s:.1f} 秒")

        if r.error_message:
            print(f"错误信息: {r.error_message}")

        if r.success and r.final_result:
            fr = r.final_result
            print(f"全井段最终有效顶替效率: {fr.get('全井段最终有效顶替效率', 'N/A')}")
            print(f"CBL评价井段水力有效顶替效率: {fr.get('CBL评价井段水力有效顶替效率', 'N/A')}")
            print(f"目标层段水力有效顶替效率: {fr.get('目标层段水力有效顶替效率', 'N/A')}")
            print(f"最终水泥浆占据率: {fr.get('最终水泥浆占据率', 'N/A')}")
            print(f"最终窜槽指数: {fr.get('最终窜槽指数', 'N/A')}")
            print(f"最终混浆指数: {fr.get('最终混浆指数', 'N/A')}")
            print(f"最终失稳指数: {fr.get('最终失稳指数', 'N/A')}")


def main() -> None:
    """主测试函数。"""
    print("cemdisp 模型全井测试")
    print(f"项目根目录: {PROJECT_ROOT}")
    print(f"测试时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    # 定义要测试的井（7口）
    wells = [
        ("呼101", "run_hu101_tailpipe_initial"),
        ("呼102", "run_hu102_tailpipe_initial"),
        ("呼103", "run_hu103_tailpipe_initial"),
        ("呼探1", "run_hu1_tailpipe_initial"),
        ("呼探1-002", "run_hu2_tailpipe_initial"),
        ("呼探1-001", "run_ht1_001_tailpipe_initial"),
        ("呼探1-003", "run_ht1_003_tailpipe_initial"),
    ]

    # 导入运行器模块
    try:
        from cemdisp.runners import (
            run_hu101_tailpipe_initial,
            run_hu102_tailpipe_initial,
            run_hu103_tailpipe_initial,
            run_hu1_tailpipe_initial,
            run_hu2_tailpipe_initial,
            run_ht1_001_tailpipe_initial,
        )
        from cemdisp.runners.ht1_003_tailpipe import run_ht1_003_tailpipe_initial
    except ImportError as e:
        print(f"错误: 无法导入运行器模块 - {e}")
        sys.exit(1)

    # 运行器函数映射
    runner_map = {
        "run_hu101_tailpipe_initial": run_hu101_tailpipe_initial,
        "run_hu102_tailpipe_initial": run_hu102_tailpipe_initial,
        "run_hu103_tailpipe_initial": run_hu103_tailpipe_initial,
        "run_hu1_tailpipe_initial": run_hu1_tailpipe_initial,
        "run_hu2_tailpipe_initial": run_hu2_tailpipe_initial,
        "run_ht1_001_tailpipe_initial": run_ht1_001_tailpipe_initial,
        "run_ht1_003_tailpipe_initial": run_ht1_003_tailpipe_initial,
    }

    # 运行所有井
    results: list[WellResult] = []
    total_start_time = time.time()

    for well_name, runner_name in wells:
        runner_func = runner_map.get(runner_name)
        if runner_func is None:
            print(f"警告: 找不到运行器 {runner_name}")
            continue

        result = run_well(well_name, runner_name, runner_func)
        results.append(result)

    total_time = time.time() - total_start_time

    # 打印结果汇总
    print_summary_table(results)
    print_detailed_results(results)

    # 保存结果到 JSON 文件
    output_data = {
        "test_time": time.strftime("%Y-%m-%d %H:%M:%S"),
        "total_time_s": total_time,
        "results": []
    }

    for r in results:
        output_data["results"].append({
            "well_name": r.well_name,
            "runner_name": r.runner_name,
            "success": r.success,
            "run_time_s": r.run_time_s,
            "error_message": r.error_message,
            "final_result": r.final_result,
        })

    output_path = PROJECT_ROOT / "results" / "all_wells_test_result.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(output_data, f, ensure_ascii=False, indent=2)

    print(f"\n结果已保存到: {output_path}")
    print(f"总耗时: {total_time:.1f} 秒")

    # 返回退出码
    fail_count = sum(1 for r in results if not r.success)
    sys.exit(fail_count)


if __name__ == "__main__":
    main()
