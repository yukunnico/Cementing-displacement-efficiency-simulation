# 10033 关键结果文件建模说明

## 1. 文档目的

这份文档不是泛泛解释整个项目，而是只针对下面 4 个结果文件，说明：

1. 它们分别表达什么；
2. 它们在代码里是怎么生成的；
3. 它们背后用到了哪一层模型；
4. 哪些公式是直接来自论文，哪些只是受论文启发的降阶表达；
5. 这些结果应该如何正确解读。

本次重点解释的 4 个结果文件是：

- `results/10033_tailpipe_139_displacement.gif`
- `results/10033_tailpipe_219_displacement.gif`
- `results/10033_tailpipe_efficiency.png`
- `results/10033_tailpipe_efficiency_time.png`

---

## 2. 先说结论：这 4 个文件不是来自同一层模型

这 4 个结果文件虽然都来自同一个项目，但它们在模型层次上并不一样：

| 文件 | 本质 | 主要模型层 |
|---|---|---|
| `10033_tailpipe_139_displacement.gif` | 动态顶替示意图 | 偏心修正层 + 动态前沿采样 |
| `10033_tailpipe_219_displacement.gif` | 动态顶替示意图 | 偏心修正层 + 动态前沿采样 |
| `10033_tailpipe_efficiency.png` | 主结果图 | Yang 基线 + 覆盖率回退 + 偏心修正终值 |
| `10033_tailpipe_efficiency_time.png` | 时间演化图 | 偏心修正终值 + 降阶时间曲线 |

所以一定要注意：

- GIF 不是“二维 CFD 动画”；
- `efficiency.png` 不是直接画 Yang 原始输出；
- `efficiency_time.png` 也不是直接来自 Sun 2020 的完整界面数值解。

---

## 3. 这 4 个结果文件共享的总工作流

这 4 个结果文件共享同一个总工作流背景：

```text
cases.py
  ↓
geometry.py + fluid.py
  ↓
screening.py
  ├─ yang2019.py
  ├─ eccentricity_corrected.py
  ├─ sun2020_screening.py
  └─ transient.py
  ↓
reporting.py
  ↓
results/*.png / *.gif
```

对应主入口是：

- `scripts/simulate_10033_tailpipe.py`

但要特别注意：

- **共享工作流** 不等于 **每个结果文件都直接数值依赖同样的子模块**；
- GIF、主结果图、时间曲线图虽然都经过 `screening.py` 总装层，但它们真正直接消费的中间结果并不完全相同；
- 下面各小节会分别说明每个结果文件的**直接生成链路**。

---

## 4. 两个 GIF 的建模说明

### 4.1 对应文件

- `results/10033_tailpipe_139_displacement.gif`
- `results/10033_tailpipe_219_displacement.gif`

这两个 GIF 使用的是同一套建模逻辑，只是案例参数不同。

### 4.2 代码直接生成链路

代码入口是：

- `scripts/simulate_10033_tailpipe.py`

最终导出函数是：

- `cemdisp/reporting.py` 中的 `write_displacement_gifs()`

它所消费的中间结果来自：

- `cemdisp/screening.py` 中 `screen_tailpipe_case()` 生成的 `schematic_series`

而 `schematic_series` 又来自：

- `cemdisp/correlations/eccentricity_corrected.py` 中 `sample_front_progress()`

因此，这两个 GIF 的**直接数值链路**可以写成：

```text
simulate_10033_tailpipe.py
    → screening.py / screen_tailpipe_case()
    → eccentricity_corrected.py / sample_front_progress()
    → reporting.py / write_displacement_gifs()
```

需要特别说明：

- GIF 的动画帧本身**不直接使用** `transient.py` 的时间曲线点；
- GIF 的动画帧本身也**不直接使用** `yang2019.py` 的原始回归输出；
- 它们真正直接使用的是偏心修正层生成的 `schematic_series`。

### 4.3 这些 GIF 采用的模型思路

这里不是完整二维界面推进求解，而是**偏心环空的降阶运动学示意模型**。

它的核心思想是：

1. 偏心以后，宽边间隙大、窄边间隙小；
2. 间隙不同会导致局部流动能力不同；
3. 同样的推进时间下，宽边前沿更快，窄边前沿更慢；
4. 因此可以采样出：
   - 窄边前沿
   - 平均前沿
   - 宽边前沿
5. 再把这三条前沿画成动画，得到动态顶替示意图。

### 4.4 用到的公式

#### （1）相对间隙尺度

当前代码中采用：

\[
g(\theta) = 1 - e\cos\theta
\]

其中：

- \(e\)：偏心度
- \(\theta\)：周向角度

#### （2）局部流动能力近似

当前代码假定局部流动能力与间隙平方近似成正比：

\[
M(\theta) \propto g(\theta)^2
\]

#### （3）局部速度比

\[
r(\theta) = \frac{g(\theta)^2}{\langle g(\theta)^2 \rangle}
\]

这个式子不是直接照抄单篇论文的方程编号，而是**受 Bittleston 2002 / Pelipenko 2004 / Moyers-Gonzalez 2009 的偏心环空物理图景启发**，在当前项目中实现成一个降阶速度修正层。

#### （4）周转时间

\[
t_{turn} = \frac{A_{ann}L}{Q}
\]

其中：

- \(A_{ann}\)：环空面积
- \(L\)：固井段长度
- \(Q\)：体积流量

#### （5）前沿推进分数

在任意时刻：

\[
f_i(t)=\mathrm{clip}\left(\frac{t}{t_{turn}}r_i,0,1\right)
\]

其中：

- 宽边前沿：取 \(\theta=\pi\)
- 窄边前沿：取 \(\theta=0\)
- 平均前沿：用周向平均修正因子代替

### 4.5 这些 GIF 直接用到了哪篇论文

严格说：

- **没有哪一篇论文被完整、逐式实现成这两个 GIF**；
- 它们采用的是**受文献启发的降阶运动学表达**。

最准确的表述应该是：

> 这两个 GIF 的核心思路受 Bittleston 2002、Pelipenko 2004、Moyers-Gonzalez 2009 对偏心环空位移不均匀性的物理图景启发，但当前代码中的 `gap² → velocity ratio → front fraction` 是本项目自定义的工程化降阶表达，不是这些论文 PDE 的完整数值求解结果。

### 4.6 这两个 GIF 应该怎么读

它们表达的是：

- 偏心越大，宽边越快、窄边越慢；
- 前沿推进不再是“齐头并进”，而是明显分叉；
- 蓝色区域表示“窄边前沿以内已经完成顶替的保守区域”。

所以 GIF 的重点不是“精确界面形状”，而是：

> **用动态形式展示偏心导致的顶替不同步。**

---

## 5. `10033_tailpipe_efficiency.png` 的建模说明

### 5.1 这个图表达什么

这个图表达的是：

> 在当前案例下，偏心度增加时，最终的相对修正顶替结果如何下降。

横轴：

- 套管偏心度 `e*`

纵轴：

- 相对修正顶替效率

图中红色虚线：

- `参考文献建议阈值 e*=0.33`

### 5.2 它在代码里怎么生成

代码路径是：

- `screening.py` 生成 `eccentricity_results`
- `reporting.py` 的 `write_plots()` 生成 `10033_tailpipe_efficiency.png`

### 5.3 它背后用了哪些模型层

这个图不是单一模型输出，而是三层叠加：

1. **Yang2019 基线层**
2. **覆盖率回退规则**
3. **偏心修正终值层**

### 5.4 Yang2019 直接实现的公式

代码里直接实现的是：

\[
\varphi_c = 0.9350\left(\frac{\rho_m}{\rho_s}\right)^{0.093}
\left(\frac{\rho_s}{\rho_c}\right)^{1.094}
\left(\frac{\mu_m}{\mu_s}\right)^{0.104}
\left(\frac{\mu_s}{\mu_c}\right)^{-0.267}
v^{0.058}
\]

这是项目里**直接实现**的论文公式。

但是要注意：

- 这一步先得到的是 `yang_raw_index`
- 只有当这个值落在 `[0,1]` 内时，才可以当作绝对顶替效率参考

### 5.5 覆盖率回退规则

如果 Yang 超界，则当前项目不会再把它直接当效率，而是回退到：

\[
\eta_{ref} = \min(1.0, \text{spacer coverage ratio})
\]

这不是论文原式，而是项目里的**工程保守规则**。

### 5.6 偏心终值修正公式

再在绝对基线上叠加偏心修正：

\[
\eta_{corr} = \eta_{ref}\cdot C(t_{turn})
\]

其中：

\[
C(t)=\frac{1}{N}\sum_{i=1}^{N}\mathrm{clip}\left(\frac{t}{t_{turn}}r_i,0,1\right)
\]

而终值点取：

\[
t=t_{turn}
\]

### 5.7 图中的阈值线是什么意思

图里的：

- `参考文献建议阈值 e*=0.33`

来自：

- `Sun et al. (2020)`

但这里要把表述说得更精确：

- 当前代码**没有**根据 `sun2020_indicators` 现场计算出一条“案例专属阈值线”；
- 这条红色虚线是在 `reporting.py` 中按 `e*=0.33` **固定画出** 的参考线；
- 它的作用是**文献对照和工程解释**，而不是当前项目数值求解的一部分。

它的含义不是“硬性物理分界”，而是：

> 当偏心度超过这个值后，宽边与窄边的界面推进差异明显增大，窄边残余和窜流风险显著上升，因此可作为工程警戒阈值。

所以图里画它，是为了帮助你看：

- 当前案例的偏心度如果继续上升，什么时候开始进入更不利区间。

### 5.8 这张图最准确的模型归属表述

> `10033_tailpipe_efficiency.png` 的横轴—纵轴关系不是直接来自 Sun 2020 或 Yang 2019 的单一原始图，而是“Yang2019 基线 + 覆盖率回退 + 偏心环空降阶修正”的组合结果；图中的 `e*=0.33` 红线是按 Sun 2020 工程建议值固定标注的参考线，用于解释而不是参与当前数值计算。

---

## 6. `10033_tailpipe_efficiency_time.png` 的建模说明

### 6.1 这个图表达什么

这个图表达的是：

> 在基准环空流速下，不同偏心工况的相对修正顶替效率是如何随时间逼近稳定终值的。

### 6.2 图的条件是什么

这个图不是“任意流速下”的结果，而是：

- 采用 `参考-基准` 情景的环空流速
- 技术尾管：`0.400 m/s`
- 油层尾管：`0.620 m/s`

同时比较：

- `e=0.00`
- `e=0.30`
- `e=0.60`

### 6.3 它在代码里怎么生成

代码链路是：

- `screening.py` 中调用 `build_efficiency_time_series(...)`
- `transient.py` 中逐步生成时间曲线
- `reporting.py` 中 `write_plots()` 导出 `10033_tailpipe_efficiency_time.png`

### 6.4 它背后的公式

定义归一化周转倍数：

\[
\lambda = \frac{t}{t_{turn}}
\]

再在每个时间点计算：

\[
\eta(t) = \eta_{ref}\cdot C(t)
\]

其中：

\[
C(t)=\frac{1}{N}\sum_{i=1}^{N}\mathrm{clip}(\lambda r_i,0,1)
\]

这和偏心终值层用的是同一套修正思想，只不过终值图取的是：

\[
t=t_{turn}
\]

而时间图展示的是：

\[
0 \le t \le t_{turn}
\]

### 6.5 它直接来自哪篇论文

严格说：

- 它**不是**某一篇论文里逐式可对应的完整时间推进公式；
- 它是当前项目在偏心修正层上的**时间扩展表达**。

最准确的表述应该是：

> `10033_tailpipe_efficiency_time.png` 不是 Sun 2020 完整二维离散界面推进的直接复现，而是建立在偏心修正终值模型之上的降阶时间演化表达，主要用于表示不同偏心工况达到稳定终值所需的相对时程差异。

再进一步说清楚：

- 这张图的纵轴数据直接来自 `transient.py`；
- `transient.py` 本身又直接调用 `eccentricity_corrected.py` 的修正因子；
- 因此它和 `Sun 2020` 的关系主要体现在**结果解释层**，而不是“逐式照搬求解”。

---

## 7. 四个结果文件与论文/模型的对应关系总表

| 结果文件 | 主要模型 | 公式来源 | 实现性质 |
|---|---|---|---|
| `10033_tailpipe_139_displacement.gif` | 偏心前沿推进采样 | `gap^2 → velocity ratio → front fraction` | 受 Bittleston/Pelipenko/Moyers 启发的降阶表达 |
| `10033_tailpipe_219_displacement.gif` | 偏心前沿推进采样 | 同上 | 同上 |
| `10033_tailpipe_efficiency.png` | Yang 基线 + 覆盖率回退 + 偏心终值修正 | Yang2019 式(7) + 本项目规则 + 降阶偏心修正 + Sun2020 阈值注释 | 混合型：论文原式 + 本地工程化组合 + 文献参考线 |
| `10033_tailpipe_efficiency_time.png` | 偏心终值的时间扩展 | `η(t)=η_ref·C(t)` | 本项目降阶时间表达 |

---

## 8. 本地论文链接

### Yang2019

- 文献综述中已整理公式与来源：
  - `D:\users\desktop\research\控压固井项目\cement model\docs\displacement_efficiency_review.md`

### Sun 2020

- PDF：
  - `D:\users\desktop\research\控压固井项目\cement model\cankao\Energy Science   Engineering - 2020 - Sun - Numerical modeling of motion of displacement interface in eccentric annulus.pdf`
- 文本提取：
  - `D:\users\desktop\research\控压固井项目\cement model\cankao\sun2020_full.txt`

### Bittleston 2002

- PDF：
  - `D:\users\desktop\research\控压固井项目\cement model\cankao\A_1020370417367.pdf`
- 文本提取：
  - `D:\users\desktop\research\控压固井项目\cement model\cankao\A_1020370417367_full.txt`

### Pelipenko 2004

- PDF：
  - `D:\users\desktop\research\控压固井项目\cement model\cankao\B_ENGI.0000009499.63859.f0.pdf`
- 文本提取：
  - `D:\users\desktop\research\控压固井项目\cement model\cankao\B_ENGI.0000009499.63859.f0_full.txt`

### Moyers-Gonzalez 2009

- PDF：
  - `D:\users\desktop\research\控压固井项目\cement model\cankao\s10665-008-9260-0.pdf`
- 文本提取：
  - `D:\users\desktop\research\控压固井项目\cement model\cankao\s10665-008-9260-0_full.txt`

---

## 9. 这份文档的边界

这份说明文档已经尽量把四个结果文件背后的模型来源讲清楚，但你要注意：

- GIF 和时间曲线不是完整二维界面数值解；
- 它们是当前项目的降阶表达；
- Yang2019 是直接实现的论文公式；
- Sun2020 当前是判别层，不是完整界面推进求解器。

所以最准确的理解应该是：

> 当前项目已经具备“文献驱动 + 工程可落地”的筛选能力，但还没有达到“完整论文数值复现”的层次。
