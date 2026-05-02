# cement model 建模思路、模块来源与公式说明

## 1. 文档目的

这份文档专门解释当前 `cement model` 项目的建模思路，回答以下几个问题：

1. 整个模型链路是怎么搭起来的；
2. 每个模块分别采用了哪篇论文中的哪一类模型；
3. 当前代码里真正实现了哪些公式；
4. 哪些地方只是“受论文启发的降阶表达”，还不是原文完整数值模型；
5. 项目中引用到的论文原文在本地哪里。

这份文档的定位是“模型说明书”，不是论文综述，也不是代码注释替代品。

---

## 2. 当前项目的总体建模思路

当前 `cement model` 不是完整 CFD 模型，也不是完整二维 Hele–Shaw / 2DGA 数值求解器，而是一个**文献驱动的降阶筛选模型**。它的总体思路是：

```text
现场案例数据
    ↓
几何与流体对象封装
    ↓
Yang2019 经验回归基线
    ↓
绝对基线有效性判断
    ↓
偏心环空降阶修正层
    ↓
Sun 2020 工程判别层
    ↓
时间演化层
    ↓
静态图 / 表格 / 报告 / GIF 动态示意图
```

它的核心思想不是“尽可能求得最精确流场”，而是：

- 用一个可实现、可解释、可快速复算的模型链路；
- 在数据不完备的工程场景下，优先得到**趋势可信**的筛选结果；
- 明确区分：
  - 什么是**绝对可报告量**；
  - 什么只是**诊断量或相对修正量**。

---

## 3. 各模块采用的模型与论文来源

### 3.1 `cemdisp/cases.py`

**模块职责**：把 `10033.pdf` 和参考工程资料里的案例数据组织成可计算对象。  
**模型来源**：不是论文模型，而是工程数据封装层。  
**对应对象**：

- `ReferenceOperatingEnvelope`
- `PublishedLinerCase`

**作用**：

- 统一尾管案例输入；
- 提供参考流速、参考黏度、相容性说明；
- 将复杂组合环空折算为当前筛选模型可接受的“等效环空”。

---

### 3.2 `cemdisp/geometry.py`

**模块职责**：定义环空几何。  
**模型来源**：几何基础层，不对应某一篇论文。  
**核心量**：

- 环空面积 `A_ann`
- 窄边间隙 `h_min`
- 宽边间隙 `h_max`
- 偏心度 `e`
- `standoff_ratio = h_min / h_max`

**作用**：

- 为经验模型、偏心修正层、Sun 2020 判据层提供统一几何输入。

---

### 3.3 `cemdisp/fluid.py`

**模块职责**：定义泥浆 / 隔离液 / 水泥浆的流变与密度。  
**模型来源**：流变基础层，不对应某一篇论文。  
**核心对象**：

- `FluidRheology`
- `FluidPair`
- `RheologyModel`

**当前作用**：

- Yang2019 当前仅接受牛顿流体输入；
- Sun2020 判据层直接使用塑性黏度和密度；
- 未来若接完整 2D 模型，这一层将继续扩展。

---

### 3.4 `cemdisp/correlations/yang2019.py`

**论文来源**：

- 杨谋，唐大千，等，《顶替效率评估的新模型》，《天然气工业》，2019

**模块职责**：提供当前项目的经验基线模型。  
**当前角色**：**参考基线层**，不再作为唯一的最终绝对顶替效率来源。

#### 当前实现的公式

代码中实现的是文中式（7）：

\[
\varphi_c = 0.9350\left(\frac{\rho_m}{\rho_s}\right)^{0.093}
\left(\frac{\rho_s}{\rho_c}\right)^{1.094}
\left(\frac{\mu_m}{\mu_s}\right)^{0.104}
\left(\frac{\mu_s}{\mu_c}\right)^{-0.267}
v^{0.058}
\]

其中：

- \(\rho_m\)：泥浆密度
- \(\rho_s\)：隔离液密度
- \(\rho_c\)：水泥浆密度
- \(\mu_m\)：泥浆黏度
- \(\mu_s\)：隔离液黏度
- \(\mu_c\)：水泥浆黏度
- \(v\)：环空流速

文献单位：

- 密度：g/cm³
- 黏度：Pa·s
- 流速：m/s

当前代码内部会把项目中的 SI 密度（kg/m³）自动换算到 g/cm³。

#### 当前模块的边界

- 仅支持牛顿流体；
- 对标定几何以外工况只作为经验外推；
- 若输出超出 `[0,1]`，当前项目只将其视为**杨模型原始指数**，不再直接当作绝对顶替效率报告。

---

### 3.5 `cemdisp/correlations/eccentricity_corrected.py`

**论文来源**：

- 思想上主要受以下路线启发：
  - Bittleston, Ferguson & Frigaard (2002)
  - Pelipenko & Frigaard (2004)
  - Moyers-Gonzalez & Frigaard (2009)

**模块职责**：把偏心环空的宽边/窄边推进差异压缩成一个可计算的降阶修正层。  
**当前角色**：**偏心修正终值层 + GIF 动态前沿采样层**。  
**重要说明**：它不是上述论文的完整 PDE 解，而是**受其物理图景启发的工程化降阶实现**。

#### 当前实现的核心思路

偏心后，周向不同位置的局部间隙不同。代码里用一个相对间隙尺度表示：

\[
g(\theta) = 1 - e \cos\theta
\]

然后假设局部流动能力近似与间隙平方成正比：

\[
M(\theta) \propto g(\theta)^2
\]

再通过周向平均把它归一化为局部速度比：

\[
r(\theta) = \frac{g(\theta)^2}{\langle g(\theta)^2 \rangle}
\]

其中：

- \(e\)：偏心度
- \(\theta=0\)：窄边方向
- \(\theta=\pi\)：宽边方向

#### 当前实现的终值修正公式

先定义一个环空体积周转时间：

\[
t_{turn} = \frac{A_{ann} L}{Q}
\]

其中：

- \(A_{ann}\)：环空面积
- \(L\)：固井段长度
- \(Q\)：体积流量

在任意时刻，局部替代分数写成：

\[
f_i(t) = \mathrm{clip}\left(\frac{t}{t_{turn}} r_i,\ 0,\ 1\right)
\]

再对全周取平均，得到当前整体修正因子：

\[
C(t) = \frac{1}{N} \sum_{i=1}^{N} f_i(t)
\]

最后，偏心修正终值取一个周转时刻：

\[
\eta_{corr} = \eta_{ref} \cdot C(t_{turn})
\]

其中：

- \(\eta_{ref}\)：绝对基线效率
- \(\eta_{corr}\)：偏心修正后的终值

#### 当前模块还能输出什么

它还提供 GIF 所需的三条前沿：

- 窄边前沿
- 平均前沿
- 宽边前沿

用于表达偏心导致的“宽边领先、窄边滞后”。

---

### 3.6 `cemdisp/correlations/sun2020_screening.py`

**论文来源**：

- Sun et al. (2020), *Numerical modeling of motion of displacement interface in eccentric annulus*

**模块职责**：接入当前工程数据下可直接计算的 Sun 2020 工程判别指标。  
**当前角色**：**判别层 / 解释层**，不是完整离散界面求解器。

#### 当前实现的公式

1. **修正密度差指标**

\[
F^* = \frac{(\rho_c - \rho_m) g h_{wide}^2}{\mu_c v}
\]

其中：

- \(\rho_c\)：水泥浆密度
- \(\rho_m\)：泥浆密度
- \(g\)：重力加速度
- \(h_{wide}\)：宽边最大半间隙
- \(\mu_c\)：水泥浆塑性黏度
- \(v\)：环空流速

2. **塑性黏度比**

\[
m = \frac{\mu_m}{\mu_c}
\]

3. **偏心度指标**

\[
e^* = e
\]

4. **若屈服值可得时的 Bingham 数**

\[
B_1 = \frac{\tau_{y,c} h_{wide}}{\mu_c v}
\]

\[
B_2 = \frac{\tau_{y,m} h_{wide}}{\mu_c v}
\]

其中：

- \(\tau_{y,c}\)：水泥浆屈服值
- \(\tau_{y,m}\)：泥浆屈服值

#### 当前判据含义

代码中已经接入的文献判别逻辑包括：

- `e* < 0.33`：参考文献建议阈值
- `m < 1`：黏度比小于 1 更有利
- 若 `B2` 可算，则进一步检查：

\[
F^* > 0 \quad \text{且} \quad B_2 < 5
\]

#### 当前模块的边界

- 当前工程数据中通常缺少稳定可靠的屈服值输入；
- 因此现在 `B1/B2` 常显示为“未知”；
- 该模块目前主要用于**判别与解释结果**，而不是替代完整界面推进求解。

---

### 3.7 `cemdisp/transient.py`

**论文来源**：

- 不是某一篇论文公式的直接照搬；
- 它建立在 `eccentricity_corrected.py` 的修正因子之上，是当前项目的降阶时间表达层。

**模块职责**：从偏心修正终值生成时间曲线。  
**当前角色**：**时间演化层**。

#### 当前实现的公式

先定义归一化周转倍数：

\[
\lambda = \frac{t}{t_{turn}}
\]

然后在一系列时间步上计算：

\[
\eta(t) = \eta_{ref} \cdot C(t)
\]

其中：

- \(\eta_{ref}\)：参考基线
- \(C(t)\)：偏心修正层在该时刻的修正因子

所以它的逻辑不是“单独发明一条时间曲线”，而是：

> “把偏心修正终值的计算规则延伸到不同时间点上。”

---

### 3.8 `cemdisp/screening.py`

**模块职责**：总装层，负责把所有子模块串起来。  
**模型来源**：不是单独论文，而是当前项目的**统一编排层**。  
**当前角色**：主工作流入口。

#### 当前工作流

1. 读取案例数据；
2. 构造情景参数；
3. 计算 Yang2019 基线；
4. 判断 Yang 输出是否在 `[0,1]`；
5. 若有效，则用 Yang 结果作绝对基线；
6. 若超界，则采用：

\[
\eta_{ref} = \min(1.0, \text{spacer coverage ratio})
\]

7. 对多个偏心度计算偏心修正终值；
8. 同步生成 Sun 2020 判据；
9. 生成时间曲线；
10. 生成 GIF 动态示意图采样；
11. 生成固定偏心、固定时刻下的流速-效率关系点。

#### 这里最重要的一点

`screening.py` 里的“覆盖率回退规则”只用于：

- 偏心层的绝对基线
- 时间层的绝对基线

它**不意味着所有 Yang 超界情景都被当成有效顶替效率**。

---

### 3.9 `cemdisp/reporting.py`

**模块职责**：把结果对象转换成报告、图表、表格和 GIF。  
**模型来源**：输出层，不对应单独论文模型。  
**当前角色**：交付层。

#### 当前输出物

- 简版 Markdown 报告
- 详细说明报告
- 情景 CSV
- JSON 结果
- 时间曲线 CSV
- 4 张 PNG
- 2 个 GIF

这个模块不改变模型本身，只负责把结果解释清楚、展示出来。

---

### 3.10 `scripts/simulate_10033_tailpipe.py`

**模块职责**：主运行入口。  
**模型来源**：脚本层。  
**当前角色**：一键复现所有结果。

运行它时，会依次：

1. 读取两个尾管案例；
2. 调 `screen_tailpipe_case(...)`；
3. 打印终端摘要；
4. 导出全部报告、图表和 GIF。

---

## 4. 关键公式在当前项目中的角色总结

### 4.1 当前真正直接实现的核心公式

1. **Yang2019 式（7）**
2. **偏心修正层的局部速度比、修正因子和终值公式**
3. **Sun 2020 指标层中的 `F*`、`m`、`e*` 以及条件可得时的 `B1/B2`**
4. **时间曲线公式：`η(t)=η_ref·C(t)`**

### 4.2 当前没有完整实现的部分

以下内容目前还没有作为完整方程求解器落地：

- Sun 2020 的完整二维界面离散推进
- Bittleston / Pelipenko 的完整 Hele–Shaw PDE 求解
- Moyers-Gonzalez 的完整静止泥浆通道解析边界
- Lockyear1984 / McLean1967 的完整相关式实现

这些目前在项目里的角色更多是：

- 物理启发
- 判据解释
- 后续扩展方向

---

## 5. 本地论文链接

下面列出当前项目中最关键的本地论文文件路径，便于你直接回看原文。

### 已在项目中直接引用的文献

- Yang2019 原文（本地参考来源在文档中已有题录，对应模型已直接实现）
  - `D:\users\desktop\research\控压固井项目\cement model\docs\displacement_efficiency_review.md`

- Sun 2020 原文 PDF
  - `D:\users\desktop\research\控压固井项目\cement model\cankao\Energy Science   Engineering - 2020 - Sun - Numerical modeling of motion of displacement interface in eccentric annulus.pdf`

- Sun 2020 提取文本
  - `D:\users\desktop\research\控压固井项目\cement model\cankao\sun2020_full.txt`

### 作为偏心修正层物理启发的重要原文

- Bittleston 2002 原文 PDF
  - `D:\users\desktop\research\控压固井项目\cement model\cankao\A_1020370417367.pdf`

- Bittleston 2002 提取文本
  - `D:\users\desktop\research\控压固井项目\cement model\cankao\A_1020370417367_full.txt`

- Pelipenko 2004 原文 PDF
  - `D:\users\desktop\research\控压固井项目\cement model\cankao\B_ENGI.0000009499.63859.f0.pdf`

- Pelipenko 2004 提取文本
  - `D:\users\desktop\research\控压固井项目\cement model\cankao\B_ENGI.0000009499.63859.f0_full.txt`

- Moyers-Gonzalez 2009 原文 PDF
  - `D:\users\desktop\research\控压固井项目\cement model\cankao\s10665-008-9260-0.pdf`

- Moyers-Gonzalez 2009 提取文本
  - `D:\users\desktop\research\控压固井项目\cement model\cankao\s10665-008-9260-0_full.txt`

### 其他重要参考

- Tardy & Bittleston 2015 PDF
  - `D:\users\desktop\research\控压固井项目\cement model\cankao\1-s2.0-S0920410514004306-main.pdf`

- Tehrani 1993 PDF
  - `D:\users\desktop\research\控压固井项目\cement model\cankao\BF00194015.pdf`

- Foroushan et al. SPE-199553 PDF
  - `D:\users\desktop\research\控压固井项目\cement model\cankao\spe-199553-ms.pdf`

- LAPSE 2023 PDF
  - `D:\users\desktop\research\控压固井项目\cement model\cankao\LAPSE-2023.6275-1v1.pdf`

---

## 6. 当前这份文档的边界

这份文档解释的是：

- 当前代码里**已经实现**了什么；
- 当前代码里**只是受到论文启发**的是什么；
- 哪些公式是直接可回指到代码的；
- 哪些仍属于未来可扩展方向。

所以它比“代码文件功能说明”更偏模型说明，但依然不是论文复现报告。

如果你愿意，我下一步可以继续给你做一份配套文档：

**“模型变量表与符号对照表（代码变量名 ↔ 文献符号 ↔ 工程含义）”**。
