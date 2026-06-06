# -*- coding: utf-8 -*-
"""cemdisp模型全井测试结果汇总报告生成脚本."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def load_json(path: Path) -> dict[str, Any]:
    """加载JSON文件."""
    with open(path, "r", encoding="utf-8-sig") as f:
        return json.load(f)


def main() -> None:
    """主函数."""
    results_dir = PROJECT_ROOT / "results"

    # 定义各井结果文件路径
    wells = [
        {
            "name": "呼101",
            "dir": "呼101尾管_1D2D耦合模型",
            "file": "呼101尾管_1D2D耦合模型_结果摘要.json",
            "status": "成功",
            "time_s": 32.3,
        },
        {
            "name": "呼102",
            "dir": "呼102尾管_1D2D耦合模型",
            "file": "呼102尾管_1D2D耦合模型_结果摘要.json",
            "status": "失败",
            "time_s": 0.0,
            "error": "缺少数据文件: hu102model/hu102_tail_caliper_inclination.csv",
        },
        {
            "name": "呼103",
            "dir": "呼103尾管_1D2D耦合模型",
            "file": "呼103尾管_1D2D耦合模型_结果摘要.json",
            "status": "成功",
            "time_s": 25.1,
        },
        {
            "name": "呼探1",
            "dir": "呼探1尾管_1D2D耦合模型",
            "file": "呼探1尾管_1D2D耦合模型_结果摘要.json",
            "status": "成功",
            "time_s": 20.5,
        },
        {
            "name": "呼探1-002",
            "dir": "呼探1-002尾管_1D2D耦合模型",
            "file": "呼探1-002尾管_1D2D耦合模型_结果摘要.json",
            "status": "成功",
            "time_s": 22.5,
        },
        {
            "name": "呼探1-001",
            "dir": "呼探1-001尾管_1D2D耦合模型",
            "file": "呼探1-001尾管_1D2D耦合模型_结果摘要.json",
            "status": "成功",
            "time_s": 21.7,
        },
    ]

    # 加载各井结果
    for well in wells:
        if well["status"] == "成功":
            json_path = results_dir / well["dir"] / well["file"]
            try:
                data = load_json(json_path)
                well["summary"] = data
                well["final_result"] = data.get("最终结果", {})
            except Exception as e:
                well["status"] = "失败"
                well["error"] = f"读取结果文件失败: {e}"

    # 打印报告
    print("=" * 100)
    print("cemdisp 模型全井测试结果汇总报告")
    print("=" * 100)
    print()

    # 1. 运行状态汇总
    print("一、运行状态汇总")
    print("-" * 100)
    print(f"{'井号':<12} {'状态':<8} {'耗时(秒)':<12} {'备注':<40}")
    print("-" * 100)
    for well in wells:
        status = well["status"]
        time_s = well.get("time_s", 0)
        error = well.get("error", "")
        print(f"{well['name']:<12} {status:<8} {time_s:<12.1f} {error:<40}")
    print("-" * 100)
    success_count = sum(1 for w in wells if w["status"] == "成功")
    print(f"总计: {len(wells)} 口井，成功 {success_count} 口，失败 {len(wells) - success_count} 口")
    print()

    # 2. 各井关键指标对比
    print("二、各井关键指标对比")
    print("-" * 100)
    header = f"{'井号':<12} {'有效顶替效率':<16} {'水泥浆占据率':<16} {'窜槽指数':<12} {'混浆指数':<12} {'失稳指数':<12}"
    print(header)
    print("-" * 100)
    for well in wells:
        if well["status"] == "成功" and "final_result" in well:
            fr = well["final_result"]
            print(
                f"{well['name']:<12} "
                f"{fr.get('全井段最终有效顶替效率', 0):<16.4f} "
                f"{fr.get('最终水泥浆占据率', 0):<16.4f} "
                f"{fr.get('最终窜槽指数', 0):<12.4f} "
                f"{fr.get('最终混浆指数', 0):<12.4f} "
                f"{fr.get('最终失稳指数', 0):<12.4f}"
            )
        else:
            print(f"{well['name']:<12} {'N/A':<16} {'N/A':<16} {'N/A':<12} {'N/A':<12} {'N/A':<12}")
    print("-" * 100)
    print()

    # 3. CBL评价井段效率
    print("三、CBL评价井段水力有效顶替效率")
    print("-" * 100)
    print(f"{'井号':<12} {'CBL评价井段效率':<20} {'目标层段效率':<20}")
    print("-" * 100)
    for well in wells:
        if well["status"] == "成功" and "final_result" in well:
            fr = well["final_result"]
            print(
                f"{well['name']:<12} "
                f"{fr.get('CBL评价井段水力有效顶替效率', 0):<20.4f} "
                f"{fr.get('目标层段水力有效顶替效率', 0):<20.4f}"
            )
        else:
            print(f"{well['name']:<12} {'N/A':<20} {'N/A':<20}")
    print("-" * 100)
    print()

    # 4. 井段信息
    print("四、井段信息")
    print("-" * 100)
    print(f"{'井号':<12} {'井段范围(m)':<30} {'环空体积(m3)':<16}")
    print("-" * 100)
    for well in wells:
        if well["status"] == "成功" and "summary" in well:
            s = well["summary"]
            interval = s.get("井段_m", [0, 0])
            volume = s.get("物理环空体积_m3", 0)
            print(f"{well['name']:<12} {interval[0]:.2f} - {interval[1]:.2f}{'':<10} {volume:<16.2f}")
        else:
            print(f"{well['name']:<12} {'N/A':<30} {'N/A':<16}")
    print("-" * 100)
    print()

    # 5. 错误说明
    print("五、问题说明")
    print("-" * 100)
    print("1. 呼102井运行失败:")
    print("   - 错误原因: 缺少数据文件 hu102model/hu102_tail_caliper_inclination.csv")
    print("   - 解决方案: 需要补充呼102井的井径/井斜数据文件")
    print()
    print("2. 呼103井目标层段效率为0:")
    print("   - 原因: 该井未设置目标层段评价窗口")
    print()
    print("3. 呼101井失稳指数为1.0:")
    print("   - 原因: 该井段较长(5400-7868m)，可能存在重力分异效应")
    print("-" * 100)
    print()

    # 6. 结果文件位置
    print("六、结果文件位置")
    print("-" * 100)
    for well in wells:
        if well["status"] == "成功":
            print(f"  {well['name']}: {results_dir / well['dir']}")
    print("-" * 100)
    print()

    print("=" * 100)
    print("测试完成")
    print("=" * 100)


if __name__ == "__main__":
    main()
