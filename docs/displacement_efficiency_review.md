# 固井水泥浆顶替效率模型综述

> 目的：为“固井水泥浆顶替效率”后续 Python 建模提供统一的文献脉络、变量命名、模型边界与实现优先级。
>
> 说明：本综述基于已检索到的公开论文页面、摘要、会议文献和中文行业论文题录整理。凡涉及具体公式常数、适用参数范围和回归系数者，正式编码前仍应再对原文全文逐条复核。

---

## 1. 问题定义与建模目标

固井顶替效率研究的核心，是描述在环空中“水泥浆/隔离液对钻井液的清除程度”。工程上它直接关系到：

1. 壁面残留钻井液厚度；
2. 窜槽与沟流风险；
3. 水泥环胶结质量与层间封隔；
4. 偏心、斜井、窄窗口等复杂井况下的参数优化。

在文献中，“顶替效率”并非只有一种定义，常见表述包括：

- **体积分数型效率**：环空内被水泥浆占据的体积分数；
- **泥浆清除率**：钻井液被驱替掉的比例；
- **壁面覆盖质量**：重点关注窄边和低侧残留；
- **界面稳定性指标**：通过界面形状、指进、混浆区长度间接表征。

因此，后续代码实现应避免只保留单一 `efficiency` 标量，而应保留“定义来源”和“计算口径”。

---

## 2. 研究对象与关键变量

### 2.1 几何与井眼变量

- 井眼直径 `hole_diameter`
- 套管外径 `casing_outer_diameter`
- 环空间隙 `annular_gap`
- 偏心度 `eccentricity`
- 最窄/最宽间隙
- 井斜角 `inclination_deg`
- 环空长度 `length`

### 2.2 流体变量

- 钻井液、隔离液、水泥浆密度 `density`
- 塑性黏度 `plastic_viscosity`
- 屈服值/屈服应力 `yield_stress`
- 稠度系数 `consistency_index`
- 流性指数 `flow_behavior_index`
- 流变模型：Newtonian / Bingham / Power-law / Herschel-Bulkley

### 2.3 工艺变量

- 排量 / 流速 `flow_rate`
- 套管旋转与上下活动
- 前置液体系与体积倍数
- 注入顺序（泥浆 → 隔离液 → 水泥浆）

### 2.4 常见无量纲量

- 密度比 / 密度差
- 黏度比
- Reynolds 数
- Hedstrom 数
- 浮力-黏性竞争参数（不同文献定义略有差异）

---

## 3. 顶替效率模型分类

文献中的顶替模型大体可以分为三类。

### 3.1 经验/实验关联式模型

这类模型以实验、现场或 CFD 数据为基础，通过回归或查图给出顶替效率与若干主控变量之间的关系。优点是实现快、工程可用性强；缺点是外推性弱、机理透明度低。

#### 代表工作与特点

1. **McLean, Manry, Whitaker (1967)**
   - 早期经典研究，讨论了固井顶替的力学过程；
   - 常被作为“工程经验模型”起点引用；
   - 对后续低排量、密度差、黏性控制等研究影响较大。

2. **Lockyear & Hibbert (1989)**
   - 强调初次固井成功的综合主控因素；
   - 关注密度、流变和作业设计对泥浆清除的影响；
   - 更偏工程设计准则而非统一显式公式。

3. **Savery & Darbe (AADE-07-NTCE-49)**
   - 基于 CFD 生成体积分数-时间曲线，并转化为工程上易用的参数图；
   - 适合快速比选密度、黏度和偏心环空条件下的泥浆清除趋势；
   - 本质上属于“CFD 支撑的经验关联式”。

4. **杨谋等（2019，《天然气工业》）**
   - 中文文献中对工程实现最友好的一类；
   - 以数值模拟结果为样本，建立顶替效率评估新模型；
   - 适合做 Python 第一阶段原型。

    其文中式（7）可直接写为：

    \[
    \varphi_c = 0.9350\left(\frac{\rho_m}{\rho_s}\right)^{0.093}
    \left(\frac{\rho_s}{\rho_c}\right)^{1.094}
    \left(\frac{\mu_m}{\mu_s}\right)^{0.104}
    \left(\frac{\mu_s}{\mu_c}\right)^{-0.267}
    v^{0.058}
    \]

    其中：\(\rho\) 单位为 g/cm³，\(\mu\) 单位为 Pa·s，\(v\) 单位为 m/s，\(\varphi_c\) 为水泥浆体积分数。作者给出的验证算例中，CFD 值为 0.9226，拟合值为 0.8893，相对误差为 4.6%。该模型依托特定物理模型条件（井眼直径 317.5 mm、套管外径 222.63 mm、模型长度 15 m）拟合得到，外推到其他几何时应视为经验估计。

#### 建模启发

- 适合先实现成 `CorrelationModel`；
- 输入一般为：密度差、流变参数、排量、偏心度、隔离液体积倍数；
- 输出是单一效率值或某一位置的体积分数；
- 最适合作为“快速扫参”和“参数敏感性分析”入口。

### 3.2 简化机理模型：Hele-Shaw / 2DGA

这是学术界最系统、最适合向“研究型代码”推进的一条路线。其思路是把窄环空中的三维流动做跨间隙平均，化为二维问题，在保持主要物理机制的同时显著降低求解成本。

#### 代表工作与特点

1. **Bittleston, Ferguson & Frigaard (2002)**
   - 经典奠基工作；
   - 面向偏心环空、层流、非牛顿流体位移；
   - 给出了 Hele-Shaw/2DGA 类问题的基本框架。

2. **Pelipenko & Frigaard (2004)**
   - 进一步研究稳态和瞬态位移；
   - 对界面推进、稳定性与残留形成机理进行了更系统处理。

3. **Carrasco-Teja et al. (2008)**
   - 将问题推广到水平井/高斜井；
   - 研究偏心环空中的黏塑性流体顶替。

4. **Bercovier, Engelman, Tardy (2015)**
   - 考虑套管运动（旋转/往复）的环空顶替模型；
   - 更接近现场工艺条件。

5. **Maleki & Frigaard (2017, 2018, 2019)**
   - 从层流进一步扩展到湍流和混合流态；
   - 说明 2DGA 路线不是只适用于极简层流场景。

#### 核心思想

- 使用窄间隙假设，把速度场、压力场和界面演化写成二维问题；
- 用流函数、体积分数或浓度输运方程描述界面推进；
- 流变方程常选 Bingham 或 Herschel-Bulkley；
- 偏心度、井斜角、密度差与流变耦合决定界面是否稳定、窄边是否残留。

#### 建模启发

- 这是后续最值得作为“项目主体模型”的路线；
- 适合放到 `models2d/` 下单独发展；
- 第一版应先实现接口与数据结构，不急于直接求 PDE。

### 3.3 高保真数值模拟：CFD / VOF / 多相流

这类方法直接求解多相流问题，通常采用 FLUENT、OpenFOAM 或其他 CFD 平台。

#### 代表工作与特点

1. **Sun et al. (2020)**
   - 建立偏心环空注水泥界面运动数值模型；
   - 研究密度差、偏心度、流变对界面平稳性的影响。

2. **张多源等相关工作**
   - 关注水泥浆和钻井液流变参数对顶替效率的影响；
   - 多用 VOF、多相流和紊流模型。

3. **Saasen et al. (2020)**
   - 大尺寸偏心环空实验与顶替效率研究；
   - 为数值模型校核提供了较有价值的实验背景。

#### 建模启发

- 不适合作为当前项目第一阶段纯 Python 自研目标；
- 更适合后续作为校核工具、数据来源或参数拟合样本来源；
- 未来若要深入，可采用“Python 驱动 CFD 后处理”的模式，而不是从零写 3D 求解器。

---

## 4. 模型比较与适用边界

| 模型类型 | 代表文献 | 主要输入 | 主要输出 | 优点 | 局限 | Python 优先级 |
|---|---|---|---|---|---|---|
| 经验/拟合模型 | McLean 1967；杨谋等 2019；Savery & Darbe 2007 | 密度差、黏度、偏心度、排量、隔离液参数 | 顶替效率或清除率 | 快，适合扫参 | 泛化弱 | 高 |
| Hele-Shaw / 2DGA | Bittleston 2002；Pelipenko 2004；Carrasco-Teja 2008 | 几何、偏心度、流变、密度差、井斜 | 界面演化、残留、效率 | 机理清楚、学术价值高 | 数值实现更复杂 | 高 |
| CFD / VOF | Sun 2020；Saasen 2020 等 | 全部几何与工艺参数 | 全场流动与界面细节 | 最细致 | 计算成本高 | 低 |

---

## 5. 对当前项目的实现建议

### 5.1 建议的三阶段路线

#### 第一阶段：经验/半经验模型

目标：尽快形成一个能跑、能扫参、能画图的原型。

建议优先实现：

1. 通用几何与流体数据结构；
2. 密度比、黏度比、环空间隙、偏心相关量；
3. 一到两个文献可追溯的经验/拟合模型；
4. 结果图：效率-排量、效率-偏心度、效率-密度差。

#### 第二阶段：2DGA 接口与数值框架

目标：为真正的机理模型留接口，而不在项目骨架阶段过度求解。

建议优先做：

1. `Model2D` 抽象基类；
2. 网格规范 `GridSpec`；
3. 结果对象 `DisplacementField`；
4. 无量纲量与校验模块复用。

#### 第三阶段：机理求解与校核

目标：接入实际 2DGA 数值求解和文献算例。

建议届时再增加：

1. 稀疏线性系统求解；
2. 体积分数/浓度输运；
3. 文献算例回归；
4. 与旧脚本或外部 CFD 结果的对照。

### 5.2 与当前工作区旧脚本的关系

当前目录中的 `环空六浆柱注入阶段模型(1)` 更像一套历史原型，特点是：

- 单脚本主程序驱动；
- 数据文件与代码耦合较紧；
- 已经包含偏心摩阻修正与环空流动计算；
- 更接近“工程算例脚本”，不太像可复用包。

因此，建议：

- **保留旧脚本原样**，不直接改；
- 在新包中建立统一数据结构与接口；
- 未来把旧脚本里的局部公式逐步迁入 `geometry.py`、`dimensionless.py` 或 `legacy/` 适配层。

---

## 6. 推荐的统一变量命名

为了避免文献符号和代码变量漂移，建议在代码中固定使用下表。

| 文献常见符号 | 建议代码名 | 含义 |
|---|---|---|
| \(\rho\) | `density` | 流体密度 |
| \(\tau_y\) | `yield_stress` | 屈服应力 |
| \(K\) | `consistency_index` | 稠度系数 |
| \(n\) | `flow_behavior_index` | 流性指数 |
| \(D_h\) | `hydraulic_diameter` | 水力直径 |
| \(e\) | `eccentricity` | 偏心度 |
| \(Q\) | `flow_rate` | 排量 |
| \(\Delta\rho\) | `density_difference` | 密度差 |
| `standoff` | `standoff_ratio` | 最窄/最宽间隙比 |

---

## 7. 研究缺口与注意事项

1. **顶替效率定义不统一**：不同论文中的 `efficiency` 不应直接横向比较。
2. **经验模型高度依赖样本范围**：对深井、高温高压、控压固井场景外推时要加警示。
3. **控压固井专用验证仍偏少**：多数顶替关联式来自常规一次固井场景。
4. **旧脚本与文献符号不统一**：后续迁移时要优先解决变量含义和单位问题。
5. **流变参数随温压变化**：若后续上升到真实井筒模拟，必须考虑温压耦合。

---

## 8. 参考文献线索

以下是当前阶段最值得继续核对全文的文献线索。

1. Bittleston, S. H., Ferguson, J., & Frigaard, I. A. (2002). *Mud removal and cement placement during primary cementing of an oil well – Laminar non-Newtonian displacements in an eccentric annular Hele-Shaw cell*. Journal of Engineering Mathematics.
2. Pelipenko, S., & Frigaard, I. A. (2004). *Mud removal and cement placement during primary cementing of an oil well – steady-state displacements*. Journal of Engineering Mathematics.
3. Carrasco-Teja, M., Frigaard, I. A., Seymour, B. R., & Storey, S. P. (2008). *Viscoplastic fluid displacements in horizontal narrow eccentric annuli*. Journal of Fluid Mechanics.
4. Bercovier, M., Engelman, M., & Tardy, P. (2015). *A model for annular displacements of wellbore completion fluids involving casing movement*. Journal of Petroleum Science and Engineering.
5. Maleki, A., & Frigaard, I. A. (2017/2018/2019). 关于湍流与混合流态下固井顶替的一组研究。
6. Sun, J., Li, Z., Luo, P., & Huang, S. (2020). *Numerical modeling of motion of displacement interface in eccentric annulus during primary cementing*. Energy Science & Engineering.
7. Savery, M., & Darbe, R. (2007). *Predicting Mud Removal During Cementing – A New and Simple Approach*. AADE-07-NTCE-49.
8. 杨谋，唐大千等（2019）. **顶替效率评估的新模型**. 《天然气工业》.
9. 侯婷等（2015）. **注水泥固井顶替的极限偏心度计算**. 《石油学报》.

---

## 9. 本综述对代码结构的直接约束

本综述不是独立文档，它将直接约束后续代码设计：

- `cemdisp/fluid.py`：承接流变参数与流体类型；
- `cemdisp/geometry.py`：承接环空几何、偏心度与 standoff；
- `cemdisp/correlations/`：放经验/拟合模型；
- `cemdisp/models2d/`：放 2DGA 类机理模型接口；
- `cemdisp/utils/dimensionless.py`：放公共无量纲量与几何辅助函数。

换句话说，**文献综述负责定义变量、边界和优先级；代码骨架负责把这些约束固化为接口。**
