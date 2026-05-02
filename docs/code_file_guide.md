# cement model 代码文件功能说明

## 1. 文档目的

这份文档专门解释 `cement model` 项目里各个代码文件的作用、相互关系，以及推荐的阅读顺序。它的目标不是替代代码本身，而是帮助你在读代码前先建立一张“地图”，知道：

- 每个文件负责什么
- 数据是如何在模块之间流动的
- 应该先读哪里，后读哪里
- 哪些文件是核心建模逻辑，哪些文件只是输出和辅助层

---

## 2. 项目整体结构

当前项目的主目录如下：

- `cemdisp/`：核心模型包
- `scripts/`：主运行脚本
- `results/`：模型输出结果
- `docs/`：说明文档
- `tests/`：单元测试
- `cankao/`：参考文献原文与资料

你真正需要重点关注的是：

1. `cemdisp/`
2. `scripts/simulate_10033_tailpipe.py`
3. `results/` 里的输出文件

---

## 3. 当前模型链路总览

当前项目不是一个全二维 CFD 求解器，而是一个**筛选级、降阶化**的尾管固井顶替效率分析工具。它的主链路可以概括为：

```text
工程案例数据(cases.py)
    → 基础几何/流体对象(geometry.py, fluid.py)
    → Yang2019 经验回归基线(yang2019.py)
    → 绝对基线修正规则(screening.py)
    → 偏心修正终值(eccentricity_corrected.py)
    → Sun 2020 工程判别(sun2020_screening.py)
    → 时间曲线(transient.py)
    → GIF 动态顶替示意图(reporting.py)
    → Markdown / CSV / JSON / PNG / GIF 输出(reporting.py)
```

所以如果你要理解整个模型是怎么工作的，真正最关键的是：

- `cases.py`
- `yang2019.py`
- `eccentricity_corrected.py`
- `sun2020_screening.py`
- `transient.py`
- `screening.py`
- `reporting.py`

---

## 4. 推荐阅读顺序

建议不要一上来就读 `screening.py`，因为那里面是总装逻辑，变量和对象比较多。更好的方式是按下面顺序读。

### 第一阶段：先理解基础对象

#### 4.1 `cemdisp/geometry.py`

这是**几何基础文件**，定义了 `WellboreGeometry`。

你在这里可以理解：

- 井眼直径 `hole_diameter`
- 套管外径 `casing_outer_diameter`
- 偏心度 `eccentricity`
- 环空面积 `annular_area`
- 宽边/窄边间隙 `max_gap` / `min_gap`
- `standoff_ratio`

这个文件回答的问题是：

> “环空在几何上长什么样？”

如果你不先理解它，后面所有偏心效应都很难看懂。

#### 4.2 `cemdisp/fluid.py`

这是**流体与流变基础文件**。

它定义了：

- `RheologyModel`：支持哪几种流变模型
- `FluidRheology`：单个流体对象
- `FluidPair`：泥浆 / 隔离液 / 水泥浆三相组合

这个文件回答的问题是：

> “参与顶替的每种流体怎么被描述？”

---

### 第二阶段：理解案例数据怎么进入模型

#### 4.3 `cemdisp/cases.py`

这是**工程案例封装文件**，非常重要。

这里把：

- `10033.pdf` 中的案例数据
- 参考文档中的参数包

整理成了代码可直接调用的对象。

关键类型：

- `ReferenceOperatingEnvelope`
- `PublishedLinerCase`

关键常量：

- `H101_TECHNICAL_LINER_219`
- `H101_OIL_LINER_139`

你在这里会看到：

- 固井段长度
- 环空容积
- 泥浆/隔离液/水泥浆密度
- 隔离液体积
- 参考流速与参考黏度

这个文件回答的问题是：

> “现场资料最终是怎么变成模型输入的？”

---

### 第三阶段：理解经验模型和修正层

#### 4.4 `cemdisp/correlations/base.py`

这是**经验模型抽象基类**。

你不需要花太多时间，但建议先看一遍，因为它规定了：

- 什么叫一个“相关式模型”
- 统一用什么接口 `evaluate(...)`
- 参数范围怎么检查

它回答的问题是：

> “所有经验模型应该长成什么统一形式？”

#### 4.5 `cemdisp/correlations/yang2019.py`

这是当前项目的**经验基线模型**。

作用：

- 根据密度比、黏度比、流速，计算 Yang2019 的原始回归输出

你要特别注意：

- 它输出的不是无条件可信的“绝对顶替效率”
- 如果结果超出 `[0,1]`，项目现在只把它当作**诊断量**
- 当前实现只接受**牛顿流体**

这个文件回答的问题是：

> “项目里的经验基线是怎么来的？”

#### 4.6 `cemdisp/correlations/eccentricity_corrected.py`

这是当前项目里最核心的一个文件之一。

作用：

- 基于偏心环空的宽边/窄边速度差异，对终值做修正
- 输出偏心修正后的终值
- 输出 GIF 动态示意图所需的宽边/中线/窄边前沿位置

关键类型：

- `EccentricityEndpointResult`
- `FrontProgressSample`
- `EccentricityCorrectedKinematicModel`

这个文件回答的问题是：

> “偏心为什么会降低顶替效果？代码里是怎么表达的？”

#### 4.7 `cemdisp/correlations/sun2020_screening.py`

这是当前项目中最接近“文献判别层”的文件。

作用：

- 直接计算 `Sun 2020` 中能落地的工程判别指标：
  - `F*`
  - `m`
  - `e*`
  - 若条件允许则加 `B1/B2`

当前它不是完整二维数值层，而是：

> “用文献中的无量纲指标来判定当前工况好不好。”

这个文件回答的问题是：

> “为什么图里会出现 `e*=0.33` 阈值？这个阈值从哪里来？”

#### 4.8 `cemdisp/correlations/lockyear1984.py`

这是一个**占位相关式文件**。

当前它已经纳入项目结构，也已经有统一的接口骨架，但还没有把文献公式正式接进来。

你可以把它理解为：

> “项目为后续继续补充老牌经验模型预留的位置。”

#### 4.9 `cemdisp/correlations/mclean1967.py`

这也是一个**占位相关式文件**。

和 `lockyear1984.py` 一样，它目前：

- 文件存在
- 接口已经准备好
- 但公式还没有真正实现

它的作用同样是：

> “提醒当前项目未来还能继续补充历史经验模型。”

---

### 第四阶段：理解时间曲线和总装流程

#### 4.10 `cemdisp/transient.py`

这是**时间演化层**。

作用：

- 从偏心修正终值出发，生成随时间逼近终值的曲线

注意：

- 它不是完整界面 PDE 的数值解
- 它是一个**降阶时间表达层**

这个文件回答的问题是：

> “时间曲线是怎么从终值推出来的？”

#### 4.11 `cemdisp/screening.py`

这是整个项目的**总装层**，最重要。

如果前面的文件你都看过，这个文件就会变得很好理解；如果你前面没看，直接看这里会很绕。

它负责：

1. 组装案例与情景
2. 调 Yang2019 跑基线
3. 判断 Yang 输出是否超界
4. 如果超界，回退到隔离液覆盖率上限做绝对基线
5. 计算偏心修正终值
6. 计算 Sun 2020 判别指标
7. 生成时间曲线
8. 生成动态示意图前沿采样
9. 生成流速-效率关系点

这个文件回答的问题是：

> “整个项目到底是怎么从输入数据一路算到最终结果的？”

---

### 第五阶段：理解输出层

#### 4.10 `cemdisp/reporting.py`

这是**结果导出层**。

它把 `screening.py` 生成的标准结果包转换成：

- 简版 Markdown 报告
- 详细说明报告
- CSV
- JSON
- 4 张 PNG
- 2 个 GIF

这个文件回答的问题是：

> “结果文件为什么会长成这样？”

#### 4.11 `scripts/simulate_10033_tailpipe.py`

这是**主运行脚本**。

如果你只想复现当前结果，不想进入包内部，这个文件最有用。

它做的事情很简单：

1. 读取两个尾管案例
2. 调 `screen_tailpipe_case(...)`
3. 打印摘要
4. 调用 `reporting.py` 导出所有结果

这个文件回答的问题是：

> “我怎么一键跑出当前结果？”

---

## 5. 其他文件的作用

### 5.1 `cemdisp/utils/validators.py`

参数范围校验工具。主要作用是：

- 如果输入超出文献标定范围，就发警告

### 5.2 `cemdisp/utils/dimensionless.py`

无量纲数和辅助计算函数，目前不是主链路核心，但以后扩展会更重要。

### 5.3 `cemdisp/models2d/base.py`

这是未来完整二维模型的预留接口，目前只是骨架，还没有真正求解器。

### 5.4 `cemdisp/__init__.py`

包导出文件。把常用对象统一导出，方便外部调用。

---

## 6. 推荐的实际阅读方式

如果你要真正读懂现在的项目，我建议按下面顺序：

1. `docs/displacement_efficiency_review.md`
2. `cemdisp/geometry.py`
3. `cemdisp/fluid.py`
4. `cemdisp/cases.py`
5. `cemdisp/correlations/yang2019.py`
6. `cemdisp/correlations/eccentricity_corrected.py`
7. `cemdisp/correlations/sun2020_screening.py`
8. `cemdisp/correlations/lockyear1984.py`（快速扫一眼即可）
9. `cemdisp/correlations/mclean1967.py`（快速扫一眼即可）
10. `cemdisp/transient.py`
11. `cemdisp/screening.py`
12. `cemdisp/reporting.py`
13. `scripts/simulate_10033_tailpipe.py`

一句话理解这个顺序：

> **先看“对象是什么”，再看“模型怎么算”，最后看“结果怎么导出”。**

---

## 7. 测试文件怎么用

如果你读完代码还不放心，可以看 `tests/` 目录。它的作用不是讲原理，而是帮助你理解“代码认为自己应该满足什么行为”。

比较值得看的是：

- `tests/test_screening.py`
- `tests/test_eccentricity_corrected.py`
- `tests/test_transient.py`
- `tests/test_reporting.py`

---

## 8. 当前代码的边界

虽然现在这套代码已经能出图、出表、出 GIF，但你读代码时一定要记住它的边界：

- 这不是完整二维 CFD 模型
- 也不是完整的 Sun 2020 数值界面推进实现
- 当前是**文献驱动的降阶筛选模型**
- 它更适合做：
  - 参数趋势比较
  - 偏心惩罚分析
  - 时程示意表达
- `screening.py` 里的“覆盖率回退规则”只用于**偏心层和时间层的绝对终值基线构造**，并不意味着所有 Yang 超界情景都会被直接当作有效顶替效率。
- 不适合直接替代：
  - 现场实测评价
  - 全二维数值求解
  - 实验标定结果

---

## 9. 你下一步怎么读最省力

如果你不想一下子看太多，我建议你先只看这 3 个文件：

1. `cemdisp/cases.py`
2. `cemdisp/screening.py`
3. `cemdisp/reporting.py`

因为这 3 个文件已经能让你快速回答：

- 输入是什么
- 中间怎么串起来
- 最后输出了什么

如果你想再深入算法细节，再去看：

- `cemdisp/correlations/yang2019.py`
- `cemdisp/correlations/eccentricity_corrected.py`
- `cemdisp/correlations/sun2020_screening.py`

---

## 10. 文档位置

这份文件已保存到：

`D:\users\desktop\research\控压固井项目\cement model\docs\code_file_guide.md`
