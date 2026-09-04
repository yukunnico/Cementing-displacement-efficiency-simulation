# 文献 B：环空顶替数值模拟方法学 —— 注释书目与综合（联网核实版）

> 调研日期：2026-09-03。全部条目经 tavily 联网核实（出版方页面 / Crossref API / Semantic Scholar / 原文 PDF / 权威综述参考文献表）。
> 核实方式逐条标注。凡未能核实作者/年份/卷号的条目一律不进入正式书目，见文末"待核/弃用清单"。
> 本调研为对照裁定性质：已知模型弱点（b³ 闭合、人工弥散、正则化、停泵冻结口径）与文献逐一对照，文献的批评如实报告。

---

## 一、注释书目

### A. 2D 缝隙平均（Hele-Shaw/2DGA）谱系源头 —— 对应 Q1、Q2

**[A1] Bittleston, S.H., Ferguson, J., Frigaard, I.A. (2002).** Mud removal and cement placement during primary cementing of an oil well — Laminar non-Newtonian displacements in an eccentric annular Hele-Shaw cell. *Journal of Engineering Mathematics*, 43(2), 229–253.
- 核实：LAPSE-2023.29256 综述参考文献表（psecommunity.org PDF）+ arXiv:2303.11043 参考文献表双重确认。
- 要点：**本谱系的源头**。把窄偏心环空展开为变间隙 Hele-Shaw cell，径向缝隙平均降维为（周向×轴向）2D 模型：流函数满足准线性 Poisson 型方程 + 浓度对流方程组 + 缝隙平均浓度闭合。演示了钻井液（n=0.7, τy=4.79 Pa）被隔离液顶替时窄边滞留通道的预测。**我们现在用的缝隙平均 2D 环空正是这一谱系。**

**[A2] Pelipenko, S., Frigaard, I.A. (2004a).** On steady state displacements in primary cementing of an oil well. *Journal of Engineering Mathematics*, 48, 1–26.（另见 pelipenko.co.uk/files/JEM2004.pdf 原文）
- 核实：原文 PDF 可访问（pelipenko.co.uk）+ arXiv:2303.11043 参考文献表。**注意：卷号存在 46 vs 48 两种记载（LAPSE 综述写 46，arXiv 写 48），引用时建议以原文 PDF 为准。**
- 要点：稳态顶替位移分析，给出屈服流体顶替的稳态判据（何种条件下界面停止推进、泥浆滞留）。

**[A3] Pelipenko, Frigaard (2004b).** Two-dimensional computational simulation of eccentric annular cementing displacements. *IMA Journal of Applied Mathematics*, 69, 557–583.
- 核实：Frigaard 官方出版列表（blogs.ubc.ca/frigaard/publications）第 128 条。
- 要点：**2D 缝隙平均模型的数值实现与收敛性**（本谱系第一代数值解法、非牛顿闭合律保持正性的构造）。对应 Q1 网格/收敛证据与 Q2 闭合。

**[A4] Pelipenko, Frigaard (2004c).** Visco-plastic fluid displacements in near-vertical eccentric annuli: lubrication modelling. *Journal of Fluid Mechanics*, 520, 343–377.
- 核实：Frigaard 官方出版列表第 129 条。
- 要点：**对 HB 流体在偏心环空推导润滑理论闭合**（压力梯度-流量非线性关系），替代简单牛顿 b³ 律。是"缝隙平均闭合修正"的标准引文（Q2）。

**[A5] Maleki, A., Frigaard, I.A. (2017).** Primary cementing of oil & gas wells in turbulent & mixed regimes. *Journal of Engineering Mathematics*, 107, 201–230.
- 核实：Frigaard 出版列表第 45 条 + arXiv:2303.11043 参考文献表。
- 要点：2DGA 扩展到湍流/混合流态；湍流 Taylor 弥散项进入 2DGA 框架（Q5）。

**[A6] Maleki & Frigaard (2019).** Comparing laminar and turbulent primary cementing flows. *J. Petroleum Science and Engineering*, 177, 808–821.
- 核实：Frigaard 出版列表第 29 条。
- 要点：层流/湍流固井顶替对比：层流顶替由浮力+流变分层主导，湍流由弥散主导——**流态决定弥散机制**（Q5）。

### B. D2DGA（弥散型缝隙平均）与 3D 对照 —— 对应 Q1、Q2、Q5（与我们框架最直接同源）

**[B1] Zhang, R., Frigaard, I.A. (2022).** Primary cementing of vertical wells: displacement and dispersion effects in narrow eccentric annuli. *Journal of Fluid Mechanics*, 947, A32. DOI: 10.1017/jfm.2022.626
- 核实：Cambridge Core 页面（在线 2022-08-25）+ Crossref API（作者 R. Zhang, I.A. Frigaard，卷 947）+ OMAE 引用格式（947: A32）。
- 要点：**"D2DGA"一词的出处**（dispersive 2D gap-averaged）：指出经典 2DGA 假设整个间隙内单一流体完全混合，在层流域失效；改用**跨间隙两层流**（顶替液沿间隙中心对称驱替）重推导闭合，把"推进型弥散"显式纳入。§2.4 专门做两模型**质量守恒验证**；给出代表性算例"显示 Hele-Shaw 方法的成功与失败"（successes and failures），并指出周向二次流可以很显著。原文（Bararpour 2025 转述）：经典 2DGA **"overestimating the volumetric efficiency"（高估体积效率）**。**这是与我们模型同谱系、但物理弥散处理更严格的版本。**

**[B2] Zhang & Frigaard (2023).** Primary cementing of vertical wells: displacement and dispersion effects in narrow eccentric annuli. Part 2: Flow behaviour and classification. *Journal of Fluid Mechanics*. DOI: 10.1017/jfm.2023.697
- 核实：Cambridge Core 页面（在线 2023-10-06，DOI 确认）+ ResearchGate 条目（题名含 Part 2: Flow behaviour and classification）。
- 要点：垂直偏心环空顶替实验（2.4–3.6 m 段等）+ 流态分类；实验中直接观察到静态层、通道、弥散尖峰（static layers and channels, dispersive spikes/fronts）。

**[B3] Bararpour, F., Frigaard, I.A. (2025).** Capturing dispersion of Herschel–Bulkley fluids in miscible primary cementing displacement flows. *Journal of Fluid Mechanics*, 1022, A15. DOI: 10.1017/jfm.2025.10773
- 核实：Cambridge Core 页面 + Crossref API（作者 Fatemeh Bararpour & Ian A. Frigaard，卷 1022，print 2025-11-10）。
- 要点：D2DGA 扩展到 HB 流体：压力梯度-流量闭合按"跨间隙分层"构造，用 TVD/minmod 通量限制器求解浓度方程以**消除陡浓度梯度处的数值振荡**。三条对方法学的直接裁定：
  1. "单流体全混跨间隙"假设**只对湍流剪切流有效**（转引 Maleki & Frigaard 2017）；层流域间隙内存在推进型界面，经典 2DGA 因此**高估体积效率**；
  2. D2DGA 的弥散项是**从降维前的剪切流方程重新推导的**（物理闭合），不是自由参数；
  3. 作者自认："D2DGA 对不同广义牛顿流变的外部验证（3D 计算/实验）仍是缺失的一步（external verification is an important missing step）"。

**[B4] Renteria, A., Frigaard, I.A. (2020).** Primary cementing of horizontal wells. Displacement flows in eccentric horizontal annuli. Part 1: experiments. *JFM*, 905, A7.
- 核实：Frigaard 出版列表第 14 条 + Cambridge 文章页。
- 要点：水平偏心环空顶替实验基准。结论（转述自 [B5]/[B6]）：2DGA 能定性预测许多现象，**"弥散与惯性效应超出 2DGA 的预测能力"**。

**[B5] Sarmadi, Renteria, Frigaard (2021).** Primary cementing of horizontal wells. … Part 2: computations. *JFM*, 915, A83.
- 核实：Frigaard 出版列表第 9 条。
- 要点：与 [B4] 实验对应的 3D 计算；**实验与 3D 计算都显示水平顶替流的弥散显著大于 Bittleston 2002 / Maleki & Frigaard 2017 模型**（转述自 [B6] 引言）。

**[B6] Renteria, A., Frigaard, I.A. (2023).** Horizontal cementing displacement flows of shear-thinning fluids, with and without casing rotation. *Geoenergy Science and Engineering*, 226, 211747.（PII S2949891023003342）
- 核实：ScienceDirect PII 页 + OMAE 引用条目。
- 要点：剪切稀化流体水平顶替（含转动）；引言明说 2DGA 谱系模型**弥散偏小**。

### C. 1D-2D 耦合与商业/工业模拟器 —— 对应 Q1、Q4

**[C1] Tardy, P.M.J., Bittleston, S.H. (2015).** A model for annular displacements of wellbore completion fluids involving casing movement. *J. Petroleum Science and Engineering*, 126, 105–123. DOI: 10.1016/j.petrol.2014.12.018
- 核实：Semantic Scholar（含完整 BibTeX：JPSE 126, 105-123, DOI）。
- 要点：**工业级 1D-3D 耦合的范本**（Schlumberger 系）：套管内 1D + 环空三维润滑（narrow-gap）+ 窄缝压力椭圆方程，含转动/往复；质量/动量通过井底/浮鞋耦合。同时给出诚实警告：**"此类模型可预测失稳起点，但失稳解析后的流动往往违反其建立的尺度假设（如忽略惯性）"**——缝隙平均模型预测失稳后的细节不可尽信。
- 备注：另有 Tardy (2018) JPSE 162:114–136 "A 3D model for annular displacements … with casing movement"（核实：Semantic Scholar + Foroushan 2021 综述引文表），把环空升级为完整 3D 润滑速度场+2D 椭圆压力方程——比 2DGA 更全面，可分辨周向二次流。

**[C2] Chiney, A., Yerubandi, K.B. (2013).** 3D Displacement Simulator Realistically Predicts Free Fall during Cementing. SPE/IADC Middle East Drilling Technology Conference & Exhibition, Dubai, 2013-10-07~09.（SPE-166867-MS）
- 核实：Wang et al. 2024 (Energies 17(5):1226) 参考文献表第 [11] 条（题名+会议+时间完整）。
- 要点：商业 3D 顶替模拟器对**套管内自由落体（free fall）与胶塞/dart** 的显式建模——Q4 的直接证据：工业模拟器把 dart/胶塞语义、U 型管、自由落体作为一等公民处理。

**[C3] Wang, N., Lamb, C., Ashok, P., van Oort, E., Granier, G., Gobert, T. (2024).** Advanced Mud Displacement Modeling for Slim Hole Cementing Operations. *Energies*, 17(5), 1226. DOI: 10.3390/en17051226
- 核实：Crossref API（六位作者、卷 17、题名）+ MDPI 全文页。
- 要点：**"套管内自由落体 + 环空单相"子模型解耦**（"free fall of the cement in the casing while the flow in the annulus is single-phase"），再整体校核；经非常规页岩井细通道固井现场案例验证。**与我们的"套管内 1D + 环空 2D"解耦哲学同构**，且明确：环空只见到水泥尾浆（替浆不进环空）。
- 备注：同一参考文献表还有 Chen, Z.; Chaudhary, S.; Shine, J. "Intermixing of cementing fluids: Understanding mud displacement and cement placement"（IADC/SPE Drilling, Fort Worth；**年份未能独立核实**，见待核清单）。

**[C4] Dai, H., Liu, G. (2018).** An Overview of Annular Displacement Efficiency in Cementing Jobs Using an Efficient Numerical Model. AADE-18-FTCE-094, 2018 AADE Fluids Technical Conference, Houston, 2018-04-10~11.
- 核实：AADE 官方 PDF（aade.org/download_file/1335/394）全文提取，作者 "Hu Dai and Gefei Liu"（Pegasus Vertex Inc.）出现在正文。
- 要点：PVI 商业模拟器（CEMPRO+ 系）的顶替效率总览：**"Snapshots of end-of-job mud/cement concentration in the annulus"** ——顶替效率的工业惯例口径是**停泵/作业结束时刻的环空浓度快照**（Q4）。验证方式：与另一主流量模拟器对比 + 文献数值对比（Dai & Liu 2017; Wang & Dai 2018）。

**[C5] Jung, Frigaard (2021).** Evaluation of common cementing practices affecting primary cementing quality. *J. Petroleum Science and Engineering*, 208, 109622. DOI: 10.1016/j.petrol.2021.109622
- 核实：Frigaard 出版列表第 3 条。
- 对应 Q4：常规固井工艺参数（排量/流变/偏心等）对顶替质量的系统评价。

### D. 3D CFD 与实验基准 —— 对应 Q1

**[D1] Skadsem, H., Kragset, S., Lund, B., Ytrehus, J.D., Taghipour, A. (2019).** Annular displacement in a highly inclined irregular wellbore: Experimental and three-dimensional numerical simulations. *J. Petroleum Science and Engineering*, 172, 998–1013.
- 核实：arXiv:2303.11043 参考文献表 + Semantic Scholar（Corpus 相关文献列表）。
- 要点：**大斜度不规则井眼 3D CFD vs 实验基准**——全 3D CFD 路线的代表验证基准（Q1）。

**[D2] Malekmohammadi, S., Carrasco-Teja, M., Storey, S., Frigaard, I.A., Martinez, D.M. (2010).** An experimental study of displacement flow phenomena in narrow vertical eccentric annuli. *JFM*, 649, 371–398.
- 核实：Frigaard 出版列表第 96 条（armarocks.net 引文表作 "laminar displacement flows"，题名两源略有差异，以 Frigaard 官方列表为准）。
- 要点：**窄垂直偏心环空顶替实验的经典基准**（LDA/PIV 级别的界面与浓度观察，通道/静态层直接可视化）。

**[D3] Carrasco-Teja, Frigaard, Seymour, Storey (2008).** Visco-plastic fluid displacements in horizontal narrow eccentric annuli. *JFM*, 605, 293–327.
- 核实：Frigaard 出版列表第 108 条。
- 要点：水平窄偏心环空屈服流体顶替实验+建模（Q1/Q3）。

**[D4] Carrasco-Teja & Frigaard (2010).** Non-Newtonian fluid displacements in horizontal narrow eccentric annuli: Effects of motion of the inner cylinder. *JFM*, 653, 137–173.
- 核实：Frigaard 出版列表第 93 条。
- 要点：内管转动/往复对水平偏心环空顶替的作用——静态滞留层可被转动剪切动员（Q3 的"静态区能否被动员"提供实验证据）。

**[D5] Singh, R., Ahmed, R., Karami, H., Nasser, M., Hussein, I. (2021).** CFD Analysis of Turbulent Flow of Power-Law Fluid in a Partially Blocked Eccentric Annulus. *Energies*, 14(3), 731. DOI: 10.3390/en14030731
- 核实：Crossref API（五位作者）。
- **对应 Q1/Q3**：部分堵塞（泥浆通道/桥堵）偏心环空的幂律湍流 CFD——堵塞段分流重分配，CFD 路线处理窄边滞留的代表。

**[D6] Eslami, Frigaard, Taghavi (2017).** Viscoplastic fluid displacement flows in horizontal channels: Numerical simulations. *J. Non-Newtonian Fluid Mech.*, 249, 79–96.
- 核实：Frigaard 出版列表第 43 条。
- 要点：屈服流体顶替的 2D 数值模拟 vs 实验基准（Q1/Q3）。

**[D7] Etrati, Roustaei, Frigaard (2020).** Strategies for mud-removal from washouts during cementing of surface casing. *J. Petroleum Science and Engineering*, 195, 107454.
- 核实：Frigaard 扩大出版列表第 17 条。
- 要点：扩径段（washout）泥浆清除策略——不规则井眼几何滞留（Q1/Q3 延伸）。

### E. 缝隙平均闭合的定量误差研究 —— 对应 Q2

**[E1] Krishna, S., Weiskirchner, S., Ravi, K., Prohaska-Marchried, M. (2025).** Semi-analytical modeling and physical insights into eccentric annular flow of Herschel–Bulkley fluids. *European Physical Journal Plus*. DOI: 10.1140/epjp/s13360-025-07054-w
- 核实：Crossref API（四位作者）+ Springer 页面（数值指标：与实验平均百分比误差 0.8–13.2%，MAPE 7.64%，R²=0.9831；与两个独立数据集误差 10.28%/7.83%）。
- Krishna 等 (2025) 半解析 HB 偏心环空闭合模型：**对"近似法"类闭合（即 b³ 类缝隙律修正）给出量化误差**——半解析闭合相对实验 MAPE 7.64%，而简单近似法误差可至 34.4%。

**[E2] Song, X. 等 (2026).** A refined semi-analytical model for hydraulic analysis of laminar power-law fluid flow in eccentric annuli. *Results in Engineering*. DOI: 10.1016/j.rineng.2026.109096
- 核实：Crossref 检索（DOI + 第一作者 Xuncheng Song）+ ScienceDirect PII S2590123026001398。
- 要点：幂律偏心环空半解析闭合，**明确以"克服以往近似法的不准确性"为目标**；给出宽侧峰值速度/剪切、非对称速度剖面、压降-流量直接换算。**b³ 类闭合的当代修正方向（幂律指数进入缝隙律）。**

**[E3] Taylor, G.I. (1953).** Dispersion of soluble matter in solvent flowing slowly through a tube. *Proceedings of the Royal Society A*, 219, 186–203.
- 核实：ADS 条目 1953RSPSA.219..186T + Semantic Scholar。
- 对应 Q5：轴向弥散理论源头。**[E4] Aris, R. (1956).** On the dispersion of a solute in a fluid flowing through a tube. *Proc. R. Soc. A*, 235, 67–77（核实：Semantic Scholar 条目存在；页码为通行记录）。

### F. 屈服流体静态区与正则化陷阱 —— 对应 Q3

**[F1] Papanastasiou, T.C. (1987).** Flows of materials with yield. *Journal of Rheology*, 31(5), 385–404.
- 核实：Springer Computational Geosciences 2017 论文参考文献表（完整题录）。
- 要点：指数正则化（Bingham-Papanastasiou）原始文献。**正则化把屈服流体换成"高粘流体"**，是伪屈服/伪流动问题的根源性设定。

**[F2] Putz, A., Frigaard, I.A., Martinez, D.M. (2009).** The lubrication paradox & use of regularisation methods for lubrication flows. *J. Non-Newtonian Fluid Mech.*, 163, 62–77.
- 核实：Frigaard 出版列表第 101 条。
- 要点：**"润滑悖论"**：正则化方法在窄缝/润滑极限流中会**定性改变解**——近静态流体层被正则化高粘度人为"润滑化"而错误流动。**这正是"窄边静止泥浆在正则化+粗网格下被错误流动起来"的机理文献**（Q3 核心证据）。

**[F3] Treskatis, T., Roustaei, A., Frigaard, I.A., Wachs, A. (2018).** Practical guidelines for fast, efficient and robust simulations of yield-stress flows without regularisation: a study of accelerated proximal gradient and augmented Lagrangian methods. *J. Non-Newtonian Fluid Mech.*, 262, 149–164.
- 核实：Frigaard 出版列表第 33 条。
- 要点：**无正则化数值方法（增广拉格朗日/近似点法）实用指南**——屈服流体模拟的"正确做法"文献标准（Q3）。

**[F4] Balmforth, Frigaard, Ovarlez (2014).** Yielding to stress: Recent developments in viscoplastic fluid mechanics. *Annual Review of Fluid Mechanics*, 46, 121–146.
- 核实：Frigaard 出版列表第 67 条。
- 要点：权威综述，明确正则化=近似、非正则化方法=精确屈服面的两分法。

**[F5] Frigaard, I.A., Paso, K.G., de Souza Mendes, P.R. (2017).** Bingham's model in the oil and gas industry. *Rheologica Acta*, 56(3), 259–282.
- 核实：Frigaard 出版列表第 51 条。
- 要点：油气领域流变模型选择与数值处理（正则化参数与物性解耦讨论）。

**[F6] Roustaei & Frigaard (2015).** Residual drilling mud during conditioning of uneven boreholes in primary cementing. Part 1: Rheology and geometry effects in non-inertial flows. *JNNFM*, 220, 87–98；Part 2: Steady laminar inertial flows. *JNNFM*, 226, 1–15.
- 核实：Frigaard 出版列表第 63、61 条。
- **Q3 关键证据**：不规则井眼中残余泥浆的"结垢层"（fouling layers）**在充分发展的层流/惯性流中持续存在**——静态滞留是鲁棒的物理现象，不是数值伪影；反之，如果数值方法把窄边静止泥浆算流动了，那是数值误差而非物理。

**[F7] Wielage-Burchard & Frigaard (2011).** Static wall layers in plane channel displacement flows. *JNNFM*, 166(5-6), 245–261.
- 核实：Frigaard 出版列表第 84 条。
- 要点：静态壁面层的解析/数值证据（Q3）。

### G. 弥散取值与工艺评价口径 —— 对应 Q4、Q5

**[G1] Maleki & Frigaard (2016).** Axial dispersion in weakly turbulent flows of yield stress fluids. *JNNFM*, 235, 1–19.
- 核实：Frigaard 出版列表第 59 条。
- 要点：**屈服应力流体湍流轴向弥散系数的定量来源**（Q5）。

**[G2] Foroushan, H.K., Lund, B., Ytrehus, J.D., Saasen, A. (2021).** Cement Placement: An Overview of Fluid Displacement Techniques and Modelling. *Energies*, 14(3), 573. DOI: 10.3390/en14030573
- 核实：Crossref API（四位作者）+ MDPI 全文。
- 要点：顶替技术+建模总综述（VOF/CFD/2DGA/润滑模型谱系表）；顶替效率影响因素矩阵（密度/流变/偏心/转速/排量）。
- 备注：LAPSE-2023.29256（psecommunity.org）为同题扩展版综述 PDF，含 2DGA 推导细节。

**[G3] 杨建波, 邓建民(西南石油大学), 冯予淇, 赵海艳, 王丹丹 (2008).** 低速注水泥时密度差对顶替效率影响规律的数值模拟研究. **石油钻探技术**, 36(5), 62–65.
- 核实：石油钻探技术官网 PDF 预览页（syzt.com.cn/cn/article/pdf/preview/109.pdf）：刊名、卷期页码、五位作者全部核实。
- 要点：FLUENT 三维 VOF 类模拟：正密度差顶替效果最好、负密度差最差；偏心+负密度差最易失稳——中文 CFD 验证文献（Q1）。

**[G4] 方春飞 等 (2016).** 井径不规则性对固井顶替效率影响规律研究. **石油机械**, 2016 年第 10 期（rhhz.net/syjxzz/html/20161001.htm）.
- 核实：rhhz.net 全文页（作者简介：方春飞，1983 年生，中国石油大学（北京）背景；期刊页眉含石油机械）。**仅第一作者+题名+期刊可核，合著者与页码未核。**
- 要点：井径不规则（扩径/缩径）对顶替效率的影响——中文井眼不规则滞留研究（Q1/Q3 延伸）。

**[G5] 《幂律流体偏心环状管流的数值模拟》. 油气储运, 2005 年第 9 期. DOI: 10.6047/j.issn.1000-8241.2005.09.005**
- **作者名未能核实**（三次检索未获作者列表），移入文末待核清单。

---

## 二、按 Q1–Q5 综合

### Q1 数值模拟谱系与验证基准

**谱系四条主线**（按"降维程度"排序）：

1. **全 3D CFD（VOF/界面捕捉）**：Skadsem 2019（不规则井眼）、Singh 2021（部分堵塞环空）、Wang 2024（细通道，CFD 子模型）[D1,D5,C3]、杨建波 2008（FLUENT）[G3]。验证基准实验：Malekmohammadi 2010 [D2]、Renteria & Frigaard 2020 [B4]、Zhang & Frigaard 2023 Part 2 [B2]、Carrasco-Teja 2008/2010 [D3,D4]。
2. **3D 润滑/窄缝模型**（保留 3D 速度场、窄缝近似动量方程）：Tardy 2018 [C1 备注]。
3. **2D 缝隙平均（Hele-Shaw/2DGA）**：Bittleston 2002 [A1] → Pelipenko 2004 [A2-A4] → Maleki & Frigaard 2017 [A5] → Zhang & Frigaard 2022 D2DGA [B1] → Bararpour & Frigaard 2025 [B3]。
4. **1D-2D/1D-3D 耦合（套管内 1D + 环空）**：Tardy & Bittleston 2015 [C1]（工业级范本）、Wang 2024 [C3]、PVI CEMPRO+ 系（Dai & Liu 2018）[C4]。

**验证惯例**（文献实际做法）：
- **质量守恒验证**：Zhang & Frigaard 2022 设专门小节验证 2DGA 与 3D 模型的质量守恒特性 [B1]；
- **网格/时间步收敛**：Tardy 2018 有限体积离散收敛性；Skadsem 2019 3D CFD 网格；Bararpour 2025 用 TVD/minmod 限制器抑制陡浓度锋的数值振荡 [B3]；
- **实验基准分级**：垂直（Malekmohammadi 2010）、水平（Renteria & Frigaard 2020）、不规则井眼（Skadsem 2019）、扩径（Etrati 2020）各有一套实验基准，且**同一团队同时给出 2D 简化模型、3D 计算、实验三者互验**[B1,B4,B5]。

**数值扩散对界面/残留通道预测的影响**（文献明确结论）：
- 经典 2DGA 因"单流体全混跨间隙"假设，**高估体积效率**（Bararpour 2025 转述 Zhang & Frigaard 2022）[B3]；
- 2DGA 定性预测好，但**弥散/惯性效应超出其能力**（Renteria & Frigaard 2020 自评）[B4]；实验与 3D 计算的弥散**显著大于** 2DGA 模型（Sarmadi 2021；Renteria 2023）[B5,B6]；
- Tardy & Bittleston 2015：缝隙平均模型可预测失稳起点，但失稳后的解析流动违反尺度假设（惯性被忽略）——**失稳后细节不可尽信**[C1]。

### Q2 缝隙平均假设（b³ 律）的有效性

**b³ 律的地位**：q ∝ b³∇p 是牛顿流体缝隙流的精确律（Hele-Shaw），2DGA 谱系对非牛顿流体一律改用**非线性压力梯度-流量闭合**，不是直接用 b³：Pelipenko & Frigaard 2004c 为 HB 推导润滑闭合 [A4]；Bittleston 2002 用混合律闭合 [A1]；D2DGA 用跨间隙两层闭合 [B1,B3]。

**幂律缝隙律与 b³ 的偏差方向**（由幂律缝隙流解析解：q ∝ b^(2+1/n)·G^(1/n)）：
- n<1（剪切稀化，水泥浆典型 n≈0.2–0.9）：指数 2+1/n>3，**真实分流比比 b³ 更极端**（宽侧更宽、窄侧更窄）→ **用 b³ 低估流路集中度/通道化程度 → 效率偏高**；
- n>1（剪切增稠）：b³ 反向偏差。
- 文献佐证：Song 2026 [E2]（幂律半解析闭合，显式含 n）、Krishna 2025 [E1]（HB 半解析，量化近似法误差至 34.4% 量级）。

**Hele-Shaw/2DGA 失效条件**（文献明示）：
- **周向二次流显著时**（Zhang & Frigaard 2022 明说）[B1]；
- 层流域（间隙内推进型界面 vs 全混假设）[B3]；
- 弥散/惯性显著时（Renteria 2020/2023）[B4,B6]；
- 大斜度/浮力主导（Pelipenko 2004c 限"近垂直"；Tardy & Bittleston 2015 惯性警告）[A4,C1]。

### Q3 屈服流体静态区数值处理

- **正则化=近似**：Papanastasiou 指数正则化把屈服流体换成低剪高粘流体 [F1]；在窄缝/润滑极限下会触发**"润滑悖论"**——正则化可定性改变解、近静态层被错误动员 [F2]。**这正是"窄边静止泥浆被错误流动起来"的机理。** 偏差方向：正则化+粗网格 → 窄边泥浆伪流动 → 通道被抹平 → **效率偏高**。
- **正确做法**：增广拉格朗日/近似点法（无正则化，屈服面精确）[F3]，或至少网格收敛性研究 + 正则化参数敏感性报告 [F4,F5]。
- **静态层是鲁棒物理**：残余泥浆结垢层在层流/惯性流中持续存在 [F6]；静态壁面层解析证据 [F7]。**停泵冻结窄边静态泥浆在物理上有依据**（但仅限已静止的部分）。
- 内管转动/往复可动员静态滞留层（实验证据 [D4]）——如模型不含转动，静止假设反而更合理。

### Q4 停泵/胶塞边界语义

- **工业模拟器惯例**：胶塞/dart 是一等公民（自由落体、U 型管、碰压显式建模）[C2]；**套管内自由落体+环空单相解耦**是成熟做法 [C3]；胶塞扫过套管内混浆、替浆不进环空是标准工艺语义 [C2,C3]。
- **顶替效率评价时刻**：文献/工业惯例是**停泵/作业结束时刻（end-of-job/EOP）的环空浓度快照**[C4]；Tardy & Bittleston 2015 的模型即终止于替浆完成 [C1]。**候凝后过程（凝胶强度恢复、浮力倒换、微环空）不属于顶替效率模拟口径**，属于环空窜流/气窜文献的范畴（静态凝胶强度、失重等另立口径）。
- **停泵冻结全场**：对"已静止的窄边泥浆"有物理依据（静态层鲁棒 [F6,F7]）；对"正在流动的宽侧"冻结会截断尾缘演化——文献口径同样在 EOP 截断，这一点我们与文献一致。

### Q5 弥散/混合系数

- **理论源头**：Taylor 1953 / Aris 1956 [E3,E4]（层流管流轴向弥散 D_L 理论）。
- **屈服流体湍流弥散**：Maleki & Frigaard 2016 [G1]（屈服应力流体弱湍流轴向弥散系数）。
- **环空层流的"弥散"主要是推进型（advective）而非 Taylor 型**：Zhang & Frigaard 2022 明确"远未到 Taylor 弥散极限，advective dispersion dominates"且经典模型未计入 [B1]；D2DGA 的弥散项是**从剪切流方程重推导的物理闭合（两层流+浮力）**，不是自由参数 [B1,B3]。
- **取值依据与偏差方向**：
  - 人工弥散取**大** → 界面糊化、窄边通道被抹平 → **效率系统性偏高**（与经典 2DGA 全混假设同向）；
  - 人工弥散取**小** → 界面过陡、通道被夸大 → 效率偏低；
  - 文献修正方向：**用 D2DGA 两层弥散闭合替换自由人工弥散系数**，或至少做弥散系数敏感性分析 + 与实验基准（Malekmohammadi 2010 垂直、Renteria 2020 水平）标定 [B1,B3,D2,B4]。

---

## 三、对我们方法学的文献对照裁定

### 3.1 有文献依据的做法

| 我们的做法 | 文献依据 | 强度 |
|---|---|---|
| 1D-2D 耦合（套管内 1D + 环空缝隙平均 2D） | Tardy & Bittleston 2015 [C1]、Wang 2024 [C3] 同构解耦哲学 | 强 |
| 环空 2D 缝隙平均（周向×轴向网格） | Bittleston 2002 → Pelipenko 2004 → Maleki 2017 → Zhang & Frigaard 2022 谱系 [A1,B1] | 强（但见 3.2 缺陷） |
| 屈服门槛（物理屈服门） | Pelipenko 2004c HB 润滑闭合 [A4]；无正则化方法 [F3] | 中-强 |
| 停泵时刻冻结 + EOP 口径 | 工业惯例 end-of-job snapshot [C4]；Tardy & Bittleston 2015 终止于替浆完成 [C1] | 强 |
| 静态泥浆保持静止（冻结） | 静态层鲁棒性 [F6,F7] | 强（仅限已静止部分） |
| 替浆不进环空（胶塞语义） | 工业模拟器胶塞/dart 显式建模 [C2,C3] | 强 |
| CFL 自适应步长 | 通用数值要求（界面捕捉类计算需限制 dt）[C4,B3] | 常规 |

### 3.2 已知缺陷（文献明确批评，如实报告）

1. **b³ 牛顿缝隙律用于剪切稀化水泥浆 = 已知闭合缺陷**。幂律缝隙律 q ∝ b^(2+1/n)·G^(1/n)，n<1 时真实分流比比 b³ 更极端 → **b³ 低估通道化 → 效率偏高**。[A4,E1,E2]（对照项目记忆：b³ 分流高偏放大问题，方向一致）
2. **人工弥散项是自由参数，无物理闭合**。文献谱系已从"经典 2DGA 无弥散（高估效率）"进化到 D2DGA 物理弥散闭合（两层流+浮力重推导）[B1,B3]；我们的人工弥散与经典 2DGA 的"全混假设"一样会使 η_E 系统性偏高，且**人工弥散系数取大 → 效率偏高**是文献一致的偏差方向 [B3,B4,B5]。
3. **单流体全混跨间隙假设在层流失效 → 效率高估**：Bararpour 2025 明说经典 2DGA "overestimating the volumetric efficiency" [B3]。我们若仍按全混 gap 平均，正中此批评。
4. **屈服流体正则化陷阱**：窄边静止泥浆在正则化+粗网格下易被伪流动（润滑悖论 [F2]）；文献标准做法=增广拉格朗日/近似点法 [F3] 或网格收敛+正则化参数敏感性报告 [F4]。若我们未报告网格收敛性与正则化参数敏感性，属于验证证据缺口。
5. **缝隙平均模型在失稳后细节、周向二次流、惯性上不可尽信**[C1,B1]——预测失稳起点可，解析失稳后场不可。
6. **"全域水泥库存比"η_E 口径在文献中无对应**：文献口径是"环空体积效率"（评价段环空内水泥占位率）或评价段效率 [C4,B3]；全域库存比含套管内存量，会把套管内水泥计入分子分母，**与 CBL 评价段不可直接比**。且库存比≥1 时 η_E≈100% 是**指标饱和**（分子≥分母按定义必然），不是顶替质量结论——文献不会接受这个口径作为"顶替好"的证据。

### 3.3 修正方向（按文献证据强度排序）

1. **闭合升级**：b³ → 幂律缝隙律 q ∝ b^(2+1/n)·G^(1/n)（半解析系数可用 [E1,E2]），或直接采用 D2DGA 两层闭合 [B1,B3]。这是文献谱系自己走过的路。
2. **弥散项物理化**：人工弥散系数 → D2DGA 两层弥散闭合（从剪切流重推导，含浮力）[B1,B3]；过渡方案：保留人工弥散但做敏感性分析并与 [D2,B4] 实验基准标定。
3. **屈服数值处理**：报告网格收敛性与正则化参数敏感性；有条件时对照无正则化方法（增广拉格朗日 [F3]）。
4. **指标口径**：并列报告"环空体积效率（评价段）"与"库存比"，η_E 与 CBL 评价段对齐；库存比≥1 时不得用 η_E≈100% 作为顶替好的证据（指标饱和），改为报告评价段效率+通道指标（窄边滞留分数）。
5. **方法学自限声明**：预测失稳起点可、失稳后细节不可信（[C1] 警告）；周向二次流/惯性/弥散超出 2DGA 能力（[B1,B4,B6]）。

---

## 四、待核/弃用清单

| 条目 | 状态 | 原因 |
|---|---|---|
| Pelipenko & Frigaard 2004a 卷号（46 vs 48） | 存疑待核 | LAPSE 综述写 46，arXiv 参考文献表写 48；原文 PDF（pelipenko.co.uk/files/JEM2004.pdf）可访问但未抓到卷页页眉，引用前以原文为准 |
| 《幂律流体偏心环状管流的数值模拟》油气储运 2005(9), DOI 10.6047/j.issn.1000-8241.2005.09.005 | 弃用→待核 | 题名+DOI 可核（pipechina 官网）但**作者名三次检索未获**；如需引用须先补作者 |
| Chen, Z., Chaudhary, S., Shine, J. "Intermixing of cementing fluids…" IADC/SPE Drilling, Fort Worth | 待核 | 存在于 Wang 2024 参考文献表（存在性可核），但年份/SPE 号未独立核实（疑 2016） |
| 方春飞等 2016 合著者与页码 | 部分核实 | 仅第一作者+题名+期刊+年份可核 |
| 方春飞等 2016 期刊名（石油机械 vs 其他） | 基本核实 | rhhz.net/syjxzz 域名对应《石油机械》，页眉见石油机械字样 |
| Malekmohammadi 2010 题名两源差异 | 已裁定 | Frigaard 官方列表："An experimental study of displacement flow phenomena in narrow vertical eccentric annuli"（armarocks 引文表作 "laminar displacement flows"，以官方列表为准） |
| JPT "Advanced Simulation Tool Developed for Deepwater Well Cementing"（2025?） | 未采信 | 会员限定内容，作者/年份/出处不可核 |
| Dwin Kurniawan UTP 学位论文（slot flow approximation） | 未采信 | 学位论文，检索页未确认作者题名一致性 |
| Zhang & Frigaard Part 2 卷号/文章号 | 部分核实 | DOI 10.1017/jfm.2023.697 与在线日期已核，卷号未抓到（Cambridge 卷页为 hash URL） |

---

## 附：核实渠道一览

- Crossref API（api.crossref.org/works/DOI）：用于 Bararpour 2025、Zhang & Frigaard 2022、Tardy 2015、Krishna 2025、Song 2026、Singh 2021、Wang 2024、Foroushan 2021、杨建波 2008（syzt 官网 PDF）等。
- Frigaard 官方出版列表（blogs.ubc.ca/frigaard/publications）：一次性核实 20+ 条题录（Pelipenko 2004×2、Maleki 2016/2017/2019/2018、Treskatis 2018、Putz 2009、Balmforth 2014、Frigaard 2017 Rheol Acta、Roustaei 2015×2、Wielage-Burchard 2011、Renteria 2020、Sarmadi 2021、Carrasco-Teja 2008/2010、Malekmohammadi 2010、Eslami 2017、Etrati 2020、Jung 2021）。
- Cambridge Core 页面：Zhang & Frigaard 2022/2023、Bararpour & Frigaard 2025（含可引用原文表述）。
- Semantic Scholar：Tardy & Bittleston 2015（完整 BibTeX）。
- 原文 PDF：AADE-18-FTCE-094（作者 Dai & Liu）、LAPSE-2023.29256 综述、pelipenko.co.uk JEM2004。
- 中文：石油钻探技术官网 PDF（杨建波 2008）、rhhz.net（方春飞 2016）。
