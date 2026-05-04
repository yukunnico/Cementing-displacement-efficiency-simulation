"""
运行器子包 — 提供各井段顶替效率模型的执行与导出功能。

运行器负责：
1. 加载井段输入数据（WellSpec / FluidSpec / PumpingSchedule）
2. 配置求解器与边界条件
3. 执行模拟并导出结果（CSV / JSON / Markdown / PNG）

每个运行器模块对应一口井的一个固井段，
输出文件名和图表标题均使用中文。
"""

from cemdisp.runners.hu101_tailpipe import run_hu101_tailpipe_initial
from cemdisp.runners.hu102_tailpipe import run_hu102_tailpipe_initial


__all__ = ["run_hu101_tailpipe_initial", "run_hu102_tailpipe_initial"]
