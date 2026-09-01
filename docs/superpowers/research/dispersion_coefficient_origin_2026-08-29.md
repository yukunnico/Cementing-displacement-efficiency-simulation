# 弥散系数溯源报告：0.018/0.015 能否在 Zhang & Frigaard 2022 中找到依据

> 任务：P0 计划 Task C1（`docs/superpowers/plans/2026-08-29-model-foundation-closure-p0-and-trackb-roadmap.md` Step 2）
> 日期：2026-08-29
> 性质：只读溯源核对，未改任何代码
> 输入：Zhang & Frigaard 2022（JFM 947 A32，D2DGA 原始论文）提取文本、Zhang & Frigaard 2023（实验篇）、Bararpour & Frigaard 2025、当前代码 `cemdisp/models2d/annulus_d2dga.py` 与 `cemdisp/models2d/d2dga_flux.py`、git 历史

---

## 0. 结论速览

| 问题 | 结论 |
|---|---|
| Q1 Zhang22 弥散项理论形式 | 弥散全部编码在**双曲通量函数**（式 4.25 第二项 I₃、式 4.28 放大因子）中；浓度方程**无任何扩散/弥散项**（Pe≈10⁷–10⁹，分子扩散被显式忽略）。**全文未给出任何弥散系数数值**，"dispersion coefficient/D_s" 一词在全文不存在 |
| Q2 0.018/0.015 有无依据 | **无据**。该系数是 2026-05-07 提交 `fa44ace`（"重修版本1"）引入的显式拉普拉斯平滑常数，docstring 中"通常 0.012–0.020"区间无文献来源；与原文任何公式均无对应。属纯标定/人为值 |
| Q3 层流不落 Taylor 弥散区 | **确认存在**。Zhang22 §4（p. A32-13）原文明确："the laminar flow is far from the Taylor dispersion regime"；Zhang23（p. 615-617）与 Bararpour25 均重申。因此"套用 Taylor 弥散公式给系数"这条路在物理上本身不成立 |
| Q4 最小处理方案 | 推荐 **方案 A（声明标定参数 + 敏感性附表）** 为主，**方案 C（理论量级参照）** 作辅助论述；方案 B（κ→0 纯物理通量）作为 Track B 方向，不在本轮 P0 内 |

**核心判断**：模型结果被 0.018/0.015 主导（HT1-004 冒烟网格扫描 scale 0→eff 1.74%、scale 1→66.5%， archived CSV 见 `results/m1_dispersion_scale/m1_dispersion_scale_scan.csv`），而该系数无文献出处、且在物理量纲上等效于一个**随网格分辨率平方变化的伪扩散率**（NZ=250 时轴向 ≈0.33 m²/s，网格加密一即衰减 4 倍）——论文若不加声明直接写"显式扩散"，是可被审稿人击穿的数据完整性风险点。

---

## 1. Q1：Zhang22 中弥散项的理论形式

### 1.1 浓度方程本身没有扩散项

Zhang22 §3（p. A32-11）在分析 2DGA 与 3-D 模型差异时明确：

> "Regarding the actual dispersion, note that in both models we have **no diffusive terms** in the concentration equations. For the ranges of typical parameters in these displacement flows, the **Péclet numbers would typically be in the range 10⁷–10⁹**, which motivates the neglect of diffusion terms, although they can also be included (Maleki & Frigaard 2017)."

即：分子扩散 D_m 被完全忽略（Pe 太高），弥散不是 D_m 的函数。§2（p. A32-13，line 270 附近）亦重申 (2.1) 右端不含扩散项（湍流情形除外，由 Maleki & Frigaard 2017 单独处理）。

### 1.2 弥散的实际载体：间隙分层的通量修正（§4）

§4（"Accounting for gap-scale dispersion"，p. A32-13 起）的推导路径：

1. **层流高 Pe 位移流不属于 Taylor 弥散区**，而属于 Yang & Yortsos (1997) 的"大 Pe、小纵横比"间隙尺度位移流态（TFE 横向平衡模型），其一维封闭通量即式 (4.1)–(4.2)（含浮力推广源自 Lajeunesse et al. 1999）；
2. **间隙分层假设**（§4.2，式 4.9–4.10）：顶替液占据间隙中心（y∈[−y_i, y_i]），被替液贴壁，间隙平均浓度 ḉ = y_i/H；
3. **两层牛顿剪切流解**（式 4.14–4.21）：流动度 I₁、浮力流动度 I₂；
4. **D2DGA 浓度通量**（§4.2.2，式 4.25）：

   $$(q_\varphi, q_\xi) = \frac{r_a}{2}\Bigl(-\partial_\xi\Psi,\ \frac{1}{r_a}\partial_\varphi\Psi\Bigr)\,\bar c\,\underbrace{\frac{m\bar c^2 + 1.5(1-\bar c^2)}{m\bar c^3 + 1-\bar c^3}}_{\text{式(4.28) 放大因子}} + \frac{\Delta\rho H^3}{6\eta_2} I_3(\bar c, m)\,[-f_\xi, f_\varphi]$$

   其中浮力弥散函数（式 4.26）：
   $$I_3(\bar c, m) = \frac{\bar c^2(1-\bar c)^3[4m\bar c + 3(1-\bar c)]}{2m[m\bar c^3 + 1-\bar c^3]}$$

5. 原文对式 (4.28) 的定位（p. A32-20）："On inspection, we see that this amplification factor is **the first term in (4.2), as derived by Yang & Yortsos (1997)**" —— 即弥散来自间隙尺度速度剖面的**对流性展宽**，是通量函数对 ḉ 的非线性依赖，不是 Fickian 扩散项。

**要点**：Zhang22 的"弥散"= 修正后的**双曲（对流型）通量函数**，弥散强度由 ḉ、m、Δρ、H、η₂ 无量纲组合内生决定，**没有任何可调弥散系数**，更没有 0.018/0.015 这类常数。全文（含附录）检索 "coefficient / D_s / diffusiv" 无一处给出弥散系数数值。

### 1.3 对照：当前代码的物理弥散通道已实现

| 原文公式 | 代码位置 | 状态 |
|---|---|---|
| 式 (4.28) 通量放大因子 | `d2dga_flux.py:45-83` `d2dga_flux_amplification`，接线 `annulus_d2dga.py:1039`（速度场半拉格朗日回溯） | 已实现（R1） |
| 式 (4.26) I₃ + 式 (4.25) 第二项浮力通量 | `d2dga_flux.py:86-142`，接线 `annulus_d2dga.py:1091-1117` | 已实现（R2，默认开关 `enable_d2dga_i3_flux`） |
| 式 (4.22)-(4.23) 两层黏度闭包 S | `annulus_d2dga.py:764-777`（I₁ 乘子，I₁ 取 Bararpour25 式 2.24） | 已实现（R3 速度场闭包） |
| 式 (4.13) 浮力力向量 f | `annulus_d2dga.py:794-803` | 已实现（R3） |

即：**原文物理弥散机制（通量修正）已在代码中，0.018/0.015 是叠加在其上的第三条人工平滑通道**（`annulus_d2dga.py:1078-1081`，对 lead/tail 用 0.018/0.015，spacer/flusher 用 0.012/0.012），两者并行存在。原文中不存在与后者对应的对象。

---

## 2. Q2：0.018/0.015 的来源考证

### 2.1 git 考古

- 引入提交：`fa44ace`（2026-05-07，"重修版本1"），`_smooth_dispersion` 与常数 0.018/0.015/0.012 一次性写入，同时替换了更早的量纲化扩散模型 `dax/day`（后者有随速度、间隙、偏心度、浮力数变化的公式，见 `a8878f2` 2026-05-03 "改进弥散模型" 时的形态）。
- 提交信息与 diff 中**无任何文献引用或推导**；docstring 中"论文版采用显式二阶差分……轴向弥散系数通常 0.012–0.020"指的是本队论文稿的做法描述，"通常"区间没有出处。
- 后续仅有的两个相关提交均为工程性修补：`62fa526`/`64512d6`（M1：系数按 dt 归一 + spacer/flusher 字面量 0.012 逐位复现），不涉及系数取值依据。

### 2.2 与原文公式无对应

§4 全部公式（4.1–4.29）中不存在任何 per-step 常数弥散系数；式 (4.25) 的两个通量项均由几何/物性参数（H、Δρ、η₂、m、f 向量）无量纲化内生给出。0.018/0.015 无法映射到其中任何一项。

### 2.3 物理量纲审查：等效伪扩散率（本报告估算）

`_smooth_dispersion` 的更新式 `f_i += κ(f_{i+1}−2f_i+f_{i−1})` 与显式扩散格式 `f_i += (D·dt/Δx²)(f_{i+1}−2f_i+f_{i−1})` 对比，得等效扩散率：

$$D_\text{eff} = \frac{\kappa}{\Delta t}\,\Delta x^2 = \frac{\kappa_0}{\Delta t_\text{ref}}\,\Delta x^2 \quad (\text{M1 归一后：}\kappa=\kappa_0\,\Delta t/\Delta t_\text{ref})$$

以 HT1-004 尾管段（悬挂器 5243.2 m → 井底 7376.7 m，L≈2133 m）为例（轴向坐标 Δs = L/(NZ−1)）：

| 网格 | Δs (m) | D_eff,axial (m²/s) |
|---|---|---|
| NZ=30（冒烟网格，archived 扫描即此口径） | ≈73.6 | ≈24.4 |
| NZ=250（八井生产基线） | ≈8.57 | **≈0.33** |
| NZ=500（扫描脚本常量） | ≈4.28 | ≈0.082 |

参照量级（本报告估计，非文献值）：间隙尺度对流弥散的等效扩散率 D ~ ΔU·b，取环空平均流速 U~1–2 m/s、半间隙 b~0.012 m，得 D ~ 0.01–0.02 m²/s。即 **NZ=250 生产网格下人工弥散约为物理机制的 15–30 倍，且 D_eff ∝ Δs²——网格加密一倍弥散衰减 4 倍，网格收敛性不成立**。M1 的 dt 归一只修好了时间维，空间维（Δx²）仍是网格依赖的。

方位角方向：y 坐标为弧长 `y = π·mean_radius`（`annulus_d2dga.py:401`），HT1-004 下 y 跨度≈0.30 m、NY=40 → Δy≈7.7 mm，D_eff,azimuthal ≈ (0.015/4)×(7.7e-3)² ≈ **2.2×10⁻⁷ m²/s**，与轴向相差 6 个数量级。所谓"方位角弥散 0.015"在物理单位下几乎不可见，其实际作用是指标空间的逐步前缘平滑——两个系数都不构成有物理意义的"方位角弥散系数"。

### 2.4 判定

**0.018/0.015 为纯标定/人为值**：无文献出处、无推导、无原文对应公式；其量纲身份是"每参考步长的指标空间平滑幅值"，而非物理弥散系数。文档记忆中的"弥散钉死"（HT1-004 scale 0→1.74%、scale 1→66.5%，mixing 0.53、前缘长度 81 m→2099 m）与该判定一致：结果形态由该常数主导，而非由 D2DGA 物理通量主导。

> 口径注意：archived CSV（`results/m1_dispersion_scale/m1_dispersion_scale_scan.csv`，提交 `1b91789` 2026-08-29）只有 scale∈{0,1}×CFL∈{on,off} 四行、单 case 耗时 ≈4 s，与脚本冒烟模式（NZ=30, NY=10）吻合，**是冒烟网格数字**。定性结论（κ 主导结果）在冒烟网格下成立且方向明确，但论文敏感性附表必须在生产网格（NZ=250/500）重跑（见 §4 方案 A）。

---

## 3. Q3：原文是否有"层流不落在 Taylor 弥散区"的论述

**有，且是三篇论文一贯的立论前提**：

1. **Zhang22 §4（p. A32-13）**（弥散建模动机段）：
   > "In Maleki & Frigaard (2017) this idea is made more explicit in modelling the Taylor dispersion in **turbulent** annular flows. However, **here the laminar flow is far from the Taylor dispersion regime**. Instead, the flows are characterised as large Péclet number, small aspect ratio miscible displacement flows on the annulus gap-scale."

2. **Zhang23（实验篇，p. 615-617）**：
   > "In turbulent flow regimes, the fluids mix rapidly across the annular gap, leading to spreading of the displacement front via Taylor dispersion; see Maleki & Frigaard (2017). However, **laminar displacements are at high Péclet number, outside of the laminar Taylor dispersion regime. Thus, advective dispersion over shorter time scales must be accounted for** in any realistic model: hence, the D2DGA model of Zhang & Frigaard (2022)."

   Zhang23 §2（p. 322-324 行）另给出实验口径 Pe ⩾ 10⁵、"To resolve diffusive effects with Pe in this range requires very fine meshes"。

3. **Bararpour & Frigaard 2025（HB 弥散篇）**：line 126-139 重申 "Certainly we are far from the classical Taylor dispersion regime where… This results in **advective dispersion**…"；line 17："Taylor dispersion and laminar advective dispersion occur, depending on flow regime."

**推论**：Taylor 弥散（D_T ~ b²U²/D_m 类公式）只在湍流/跨间隙扩散平衡建立时成立；本模型对象是层流高 Pe 顶替，**任何"按 Taylor 公式反推一个弥散系数"的做法都与原文立论相抵触**。这也正是综合报告（obsidian 模型与文献差距分析 2026-08-23）"幂律层流指数 > b³ 陷阱/层流不落 Taylor 弥散区"论断的原文出处。

---

## 4. Q4：若系数无据——最小处理方案

### 方案 A（推荐，P0 内可完成）：公开声明为标定参数 + 敏感性附表

- 论文 §1.5"数值格式"处不写"显式扩散"一笔带过，改为明确一段："浓度方程在 D2DGA 通量闭包（式 4.25/4.26/4.28）之外，附加显式 Laplacian 数值平滑以稳定半拉格朗日前缘，轴向/方位系数 κ_s=0.018/0.015（前置液 0.012），**属数值正则化/标定参数，非物理弥散系数，原文无对应项**；按 Δt 归一（M1）"。
- 附敏感性表：在生产网格（NZ=250）重跑 scale ∈ {0, 0.25, 0.5, 1.0}（如时间允许加 1.5/2.0），报告 eff / mixing / 窜槽指数 / 前缘长度；同时报告每档对应的等效 D_eff（m²/s，见 §2.3 换算式），把"网格依赖"摊开在明面上。
- 现成脚本：`scripts/dispersion_scale_sensitivity_scan.py`（改 NZ 为 250 重跑即可；archived 结果是冒烟网格，不能直接入论文）。
- 成本：半天～1 天（8 井中选 2–3 口代表井即可，不必 8 井全扫）。

### 方案 B（Track B 方向，本轮不做）：κ→0，回归纯物理通量

原理上 D2DGA 通量修正本身就是弥散的物理载体，人工平滑应可退场。但 archived 扫描显示 scale→0 时 HT1-004 eff 仅 1.74%、前缘 81 m——说明当前半拉格朗日+双线性插值格式的数值扩散/其他失真在 κ=0 下另有主导问题（与 2026-08-23 环空失真修正调研的 M0–M4/I3 结论一致：e 饱和、b³ 分流、I3 弱效等），**直接删 κ 会把问题暴露成"模型整体失效"而非"弥散归位"**。应在 Track B（TVD/守恒形式重构、网格收敛性验证）中处理，届时 κ 以网格收敛极限的形式消失。

### 方案 C（辅助论述）：给出理论量级参照，说明 κ 的量级身份

如审稿人问"为何取 0.018"，可答：按 §2.3 换算，NZ=250 下等效轴向扩散率 ≈0.33 m²/s，与间隙尺度对流弥散量级估计 D~ΔU·b ≈ 0.01–0.02 m²/s 相比偏大约一个数量级、但同属"间隙尺度展宽远大于分子扩散（Pe≥10⁵）"的区制；取值经 HT1-004 等井的前缘形态/效率标定，敏感性见附表。**注意此估算只作辩护材料，不写进论文正文当"依据"**。

### 不建议

- 套 Taylor-Aris / Taylor 弥散公式给 κ 换算依据（与原文层流高 Pe 立论冲突，见 §3）；
- 在论文中把 0.018/0.015 写成"依据 Zhang & Frigaard (2022) 的间隙尺度弥散"——原文无此系数，属引用失实，审稿人核对 (4.25)-(4.28) 即穿帮。

---

## 5. 逐问结论

1. **Zhang22 弥散项理论形式**：间隙分层假设下的双曲通量修正——式 (4.25) 对流通量 × 放大因子 (4.28)（= Yang & Yortsos (1997) 式 (4.2) 首项）+ 浮力弥散通量（式 4.26 I₃）；浓度方程无扩散项，Pe≈10⁷–10⁹（实验口径 Pe⩾10⁵）；**无任何数值弥散系数**。
2. **0.018/0.015 对应关系**：无对应。2026-05-07 `fa44ace` 引入的数值平滑常数；docstring 区间无出处；等效物理扩散率网格依赖（D∝Δs²，NZ=250 时 ≈0.33 m²/s，方位角方向 ≈2×10⁻⁷ m²/s）。**纯标定/人为值**。
3. **层流不落 Taylor 弥散区**：原文明确存在（Zhang22 §4 p.A32-13；Zhang23 p.615-617；Bararpour25 line 126-139），是三篇论文共同立论；因此弥散必须走对流型通量修正（已实现）而非扩散系数。
4. **最小处理方案**：方案 A（论文声明标定参数 + 生产网格敏感性附表，0.5–1 天）为 P0 主选；方案 B（κ→0）归 Track B；方案 C 作审稿答辩储备。论文中**不得**把该系数表述为源自 Zhang22。

---

## 附：证据位置索引

| 证据 | 位置 |
|---|---|
| Pe 10⁷–10⁹ / 无扩散项 | Zhang22 提取文本 line 790-793（p. A32-11） |
| 层流远离 Taylor 弥散区 | Zhang22 line 918-922（p. A32-13，§4 开篇） |
| 通量式 (4.25)、I₃ (4.26)、放大因子 (4.28)、与 (4.2) 关系 | Zhang22 line 1509-1553（p. A32-20） |
| Zhang23 层流高 Pe / 对流弥散 | Zhang23 line 613-618、line 322-324、line 1126 |
| Bararpour25 远离经典 Taylor 区 | Bararpour25 line 17、line 126-139 |
| 0.018/0.015 引入提交 | `fa44ace`（2026-05-07）diff：`_smooth_dispersion` 全文 |
| 现行接线（κ/dt 归一 + 0.012） | `cemdisp/models2d/annulus_d2dga.py:1063-1081`，`_smooth_dispersion` 定义 554-584 |
| 物理通量通道已实现 | `d2dga_flux.py`（4.28/4.26/4.25 第二项/2.24/2.25）；接线 `annulus_d2dga.py:1039, 775, 800, 1091-1117` |
| y 坐标弧长定义 | `annulus_d2dga.py:399-401` |
| HT1-004 尾管段深度 | `cemdisp/data/loaders/ht1_004_loader.py:54-56, 194-195`（5243.2–7376.7 m） |
| 敏感性 archived（冒烟网格） | `results/m1_dispersion_scale/m1_dispersion_scale_scan.csv`（提交 `1b91789`） |
| 扫描脚本（需改 NZ 重跑） | `scripts/dispersion_scale_sensitivity_scan.py` |
