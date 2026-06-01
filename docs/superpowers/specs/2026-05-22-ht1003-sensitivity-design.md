# 呼探1-003井尾浆参数敏感性分析设计文档

**日期**: 2026-05-22
**目标**: 对呼探1-003井（HT1-003）尾浆的注入排量、屈服应力、塑性黏度进行全因子敏感性分析，评估各参数对全井有效顶替效率的影响。

## 1. 背景

呼探1-003井的1D-2D耦合模拟已完成（commit 3d1ab46），基线结果：
- 全井有效顶替效率: 0.803
- CBL评价井段效率: 0.916
- 尾浆基线参数: 排量=1.2 m³/min, 屈服应力=14.0 Pa, 塑性黏度=180 mPa·s

现有 `p3_p4_integration.py` 只扫描模型内部标定参数（M, alpha），无物理流体参数扫描能力。

## 2. 参数空间

| 参数 | 范围 | 步长 | 水平数 | 基线值 |
|------|------|------|--------|--------|
| 尾浆排量 (m³/min) | 0.30~1.50 | 0.10 | 13 | 1.20 |
| 尾浆屈服应力 (Pa) | 4~24 | 2 | 11 | 14 |
| 尾浆塑性黏度 (mPa·s) | 20~400 | 20 | 20 | 180 |

全因子组合: 13 × 11 × 20 = **2,860 次模拟**

## 3. 实现方案

### 3.1 文件位置

新建独立脚本: `scripts/ht1_003_sensitivity.py`

### 3.2 数据修改策略

使用 `dataclasses.replace()` 创建修改后的副本，不改变原始 loader 数据：

- **尾浆流体参数**: 找到 fluids 元组中 `name="尾浆"` 的 FluidSpec，用 `replace(plastic_viscosity_pa_s=new_pv, yield_stress_pa=new_yp)` 创建新副本
- **尾浆排量**: 找到 PumpingSchedule 中 `step_name="注入尾浆"` 的步骤，用 `replace(rate_m3_min=new_rate)` 创建新副本
- 重建 fluids 元组和 PumpingSchedule 后喂给模拟管线

### 3.3 顶替结束时间的自洽性

改变尾浆排量会影响 `cement_end_time_s`（水泥浆尾缘越过鞋口的时间），但模型自动处理：
- `CasingFlowSolver.run()` 基于修改后的泵送计划重新计算 `cement_end_time_s`
- `_rear_arrival_time()` = (该步骤累计泵入结束体积 + 管内容积) 对应的地面时间
- 排量降低 → 尾浆步骤持续时间增加 → `cement_end_time_s` 延后
- `AnnulusD2DGASolver` 使用新的 `cement_end_time_s` 作为总时长
- 每次参数组合都重新运行完整 1D 管线，结束时间判断自洽

### 3.4 模拟流程

每个参数组合执行完整 1D-2D 耦合管线：

```
load_ht1_003_tailpipe()  →  基线数据
    ↓
dataclasses.replace()    →  修改尾浆参数/排量
    ↓
CasingFlowSolver(enable_gravity=True).run()  →  1D 套管内前沿
    ↓
build_coupled_annulus_inlet_provider(split_cement_phases=True)  →  桥接
    ↓
AnnulusD2DGASolver(total_t=cement_end_time_s, nz=500).run()  →  2D 环空
    ↓
提取 result.metrics.iloc[-1]["effective_efficiency"]  →  全井有效顶替效率
```

### 3.5 输出设计

**输出目录**: `results/呼1-003_敏感性分析/`

**CSV 汇总表**: `呼1-003_敏感性分析结果.csv`

列: 排量_m3min, 屈服应力_Pa, 黏度_mPas, 全井有效顶替效率, 耗时_s

**图表** (中文标签, matplotlib + SimHei):

1. **单因素曲线图** ×3:
   - `呼1-003_排量敏感性.png`: X=排量, Y=效率 (屈服应力=14, 黏度=180 的子集)
   - `呼1-003_屈服应力敏感性.png`: X=屈服应力, Y=效率 (排量=1.2, 黏度=180 的子集)
   - `呼1-003_黏度敏感性.png`: X=黏度, Y=效率 (排量=1.2, 屈服应力=14 的子集)

2. **两两交互热力图** ×3:
   - `呼1-003_排量-屈服应力_热力图.png`: X=排量, Y=屈服应力, 颜色=效率 (黏度固定180)
   - `呼1-003_排量-黏度_热力图.png`: X=排量, Y=黏度, 颜色=效率 (屈服应力固定14)
   - `呼1-003_屈服应力-黏度_热力图.png`: X=屈服应力, Y=黏度, 颜色=效率 (排量固定1.2)

### 3.6 执行策略

- **进度显示**: 每10次模拟打印进度
- **中间保存**: 每100次模拟保存中间CSV
- **断点续跑**: 检测已有中间CSV，跳过已完成的组合
- **最终保存**: 完成后保存完整CSV + 所有图表

## 4. 关键接口

| 接口 | 来源 |
|------|------|
| `load_ht1_003_tailpipe()` | `cemdisp.data.loaders.ht1_003_loader` |
| `dataclasses.replace()` | Python 标准库 |
| `CasingFlowSolver` | `cemdisp.transport1d.casing_flow` |
| `build_coupled_annulus_inlet_provider()` | `cemdisp.models2d.boundary_bridge` |
| `AnnulusD2DGASolver` | `cemdisp.models2d.annulus_d2dga` |

## 5. 预计计算量

- 2,860 次模拟 × ~15-60s/次（nz=500）
- 预计总时间: 12-48 小时
