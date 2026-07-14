# HT1-004 尾浆敏感性脚本设计

**日期：** 2026-06-11

## 目标

新增 `scripts/ht1_004_sensitivity.py`，用于呼1-004井尾浆参数敏感性分析。脚本以 `scripts/ht1_003_sensitivity.py` 为模板，保持主要功能和尾浆三参数扫描范围不变，并按用户确认的方案 B：对齐现有 HT1-004 runner 的环空停止时间口径，同时将二维环空网格保持为 `nz=140` 以适合批量计算。

## 范围

脚本需要实现：

- 全因子扫描尾浆排量、屈服应力、塑性黏度。
- 每组参数运行一次 1D-2D 耦合顶替模拟。
- 输出完整 CSV 和中间 checkpoint CSV。
- 支持 `--resume`、`--plot-only`、`--limit`。
- 输出三张单因素曲线和三张两两交互热力图。
- 打印统计摘要和基线效率。

## 参数空间

保持与呼1-003敏感性脚本一致：

- 排量：`0.30 ~ 1.60 m^3/min`，步长 `0.10`。
- 屈服应力：`5 ~ 25 Pa`，步长 `1`。
- 塑性黏度：`50 ~ 350 mPa·s`，步长 `20`。
- 全因子组合数：`14 × 21 × 16 = 4704`。

基线值继续使用：

- 排量：`1.20 m^3/min`。
- 屈服应力：`14.0 Pa`。
- 塑性黏度：`180.0 mPa·s`。

## 数据来源

新脚本使用：

```python
from cemdisp.data.loaders.ht1_004_loader import load_ht1_004_tailpipe
```

并通过：

```python
well_spec, base_fluids, base_schedule, _ = load_ht1_004_tailpipe()
```

加载呼1-004井基线模型。

尾浆定位规则保持为：

- 流体名称：`"尾浆"`。
- 泵注步骤名称：`"注入尾浆"`。

## 模拟流程

每个参数组合执行以下流程：

1. 将扫描得到的塑性黏度从 `mPa·s` 转换为 `Pa·s`。
2. 使用 `dataclasses.replace()` 修改 `FluidSpec("尾浆")` 的：
   - `plastic_viscosity_pa_s`
   - `yield_stress_pa`
3. 使用 `dataclasses.replace()` 修改 `PumpingSchedule` 中 `"注入尾浆"` 步骤的：
   - `rate_m3_min`
4. 使用 `CasingFlowSolver(enable_gravity=True)` 运行套管内 1D 前缘追踪。
5. 使用 `build_coupled_annulus_inlet_provider()` 将鞋口出流桥接为环空入口边界。
6. 计算环空二维停止时间。
7. 使用 `AnnulusD2DGASolver(total_t=stop_time, nz=140)` 运行环空 2D D2DGA 顶替模拟。
8. 提取最终 `result.metrics.iloc[-1]["effective_efficiency"]`。

## 环空停止时间口径

按用户确认，新脚本对齐现有 `cemdisp/runners/ht1_004_tailpipe.py` 的现场停止时间逻辑：

- 从 `casing_result.fronts` 中按顺序查找鞋口前缘。
- 先确认已经出现水泥类流体，水泥类角色包括：
  - `FluidRole.LEAD`
  - `FluidRole.INTERMEDIATE`
  - `FluidRole.TAIL`
- 一旦水泥之后的首个非水泥流体到达鞋口，就以该前缘时间作为环空二维停止时间。
- 若未找到水泥之后的非水泥流体，则回退到 `casing_result.cement_end_time_s`。

该逻辑会在新脚本内实现为 `_annulus_stop_time_s()`，避免直接依赖 runner 内部函数，同时保持行为一致。

## 输出

输出目录：

```text
results/呼1-004_敏感性分析/
```

主要文件：

```text
呼1-004_敏感性分析结果.csv
呼1-004_敏感性分析_中间结果.csv
呼1-004_排量敏感性.png
呼1-004_屈服应力敏感性.png
呼1-004_黏度敏感性.png
呼1-004_排量-屈服应力_热力图.png
呼1-004_排量-黏度_热力图.png
呼1-004_屈服应力-黏度_热力图.png
```

CSV 字段保持与 003 敏感性脚本一致：

- `排量_m3min`
- `屈服应力_Pa`
- `黏度_mPas`
- `全井有效顶替效率`
- `耗时_s`

## 命令行接口

保留以下用法：

```powershell
python scripts/ht1_004_sensitivity.py
python scripts/ht1_004_sensitivity.py --resume
python scripts/ht1_004_sensitivity.py --plot-only
python scripts/ht1_004_sensitivity.py --limit 100
```

## 验证

实现后执行：

```powershell
python -m py_compile scripts/ht1_004_sensitivity.py
python scripts/ht1_004_sensitivity.py --limit 1
```

成功标准：

- 语法检查通过。
- `--limit 1` 能完成一次 1D-2D 模拟。
- 生成 CSV 和图表。
- 统计摘要能正常打印。

## 非目标

本次不做以下改动：

- 不重构 003 和 004 为通用敏感性框架。
- 不改变三参数范围和步长。
- 不并行化 4704 次模拟。
- 不改变物理模型实现。
- 不修改 HT1-004 loader 或 runner。

## 自检

- 没有未定需求或占位符。
- 方案与用户确认一致：对齐 004 runner 停止时间，`nz=140`。
- 输出命名与现有中文结果目录风格一致。
- 实现范围集中在新增脚本，避免影响已有 003 脚本和 004 runner。
