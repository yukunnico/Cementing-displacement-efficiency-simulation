# 文献 A：居中度/偏心度效应量级与机制（注释书目+综合裁定）

> 深度调研团队·文献 A 调研报告。检索工具：tavily CLI（search/extract），检索日期 2026-09-03。
> 铁律执行说明：下列条目均在本次检索中经 OnePetro / Cambridge Core / Elsevier ScienceDirect / AIP / Springer / MDPI / 中文期刊官网等一手或权威二手页面核实存在；无法核实者一律列入文末"弃用清单"。
> 口径约定：偏心度 e（0=完全居中，1=完全贴边）；居中度/standoff = 1 − e。我们的 η_E 为全井体积效率口径。

---

## ① 注释书目

### A 组：奠基经典（1960s–1990s）

**[A1] McLean, R.H., Manry, C.W., Whitaker, W.W. (1967). Displacement Mechanics in Primary Cementing. *Journal of Petroleum Technology*, 19(2), 251–260.**
- 核实：OnePetro JPT 卷期页完整（onepetro.org/JPT/article/19/02/251/162930）。
- 要点：固井顶替机理奠基论文。首次系统提出偏心环空中水泥沿宽边先窜、窄边泥浆滞留成窜槽的图像；核心思想是"壁面剪切应力是否超过泥浆凝胶强度/屈服应力"的判据；并观察到把水泥变稀追求紊流反而加重窜槽（紊流水泥走宽边、厚泥浆层流滞窄边）。
- 回答：Q3（窄边窜槽机理奠基）、Q2（正方谱系源头）、Q4（工业判据思想源头）。

**[A2] Walton, I.C., Bittleston, S.H. (1991). The axial flow of a Bingham plastic in a narrow eccentric annulus. *Journal of Fluid Mechanics*, 222, 39–60.**
- 核实：Cambridge Core 页面（摘要+结论段落）。
- 要点：Bingham 塑性体窄偏心环空轴流的 δ 展开（缝隙比小参数）解析+数值解。流动由局部槽流解构成：大部分环空呈"伪塞"结构，最宽（有时最窄）处出现真塞。**给出流体在窄边停止流动的简单判据**，并给出静止区（motionless region）范围的上下界。窄边停流的物理根源：局部壁面剪切应力 τ_w = (dp/dz)·b_narrow/2 低于屈服应力 τy。
- 回答：Q3（静止区理论判据）、Q5（局部槽流解而非全局 b³ 的理论出处）。

**[A3] Tehrani, M.A., Bittleston, S.H., Long, P.J.G. (1993). Flow instabilities during annular displacement of one non-Newtonian fluid by another. *Experiments in Fluids*, 14, 246–256.**
- 核实：经 AADE-18-FTCE-094 参考文献 + Energies 2021 综述两处独立转引核实。
- 要点：首次实验报告垂直偏心环空的周向二次流（azimuthal secondary flow）——周向输送与界面失稳的来源；建议层流顶替驱替/被顶替流体间 10–15% 正密度差。
- 回答：Q1（密度差补偿机制）、Q3。

**[A4] Haut, R.C., Crook, R.J. (1979). Primary Cementing: The Mud Displacement Process. SPE-8253-MS, SPE Annual Technical Conference and Exhibition, Las Vegas, Sept. 1979.**
- 核实：OnePetro 页面（cited by 181），DOI 10.2118/8253-MS。
- 要点：经典实验系列：泥浆条件（凝胶强度、失水造壁）直接决定顶替结果；泥饼与静止凝胶泥浆是残留主因，隔离液/冲洗液的功能是破坏凝胶与稀释泥饼。
- 回答：Q3（凝胶强度条件）、Q2 正方谱系。

**[A5] Keller, S., Crook, R., Haut, R., Kulakofsky, D. (1987). Deviated-wellbore cementing: Part 1—problems. *JPT*, 39, 955–960.**
- 核实：经 Physics of Fluids 34, 053610 (2022) 参考文献列表核实。
- 要点：斜井固井顶替问题谱：低边泥浆滞留、固相沉积、水窜；偏心与重力分异在斜井叠加。
- 回答：Q3、Q4（67% 阈值在斜井不充分的依据之一）。

**[A6] Sabins, F.L. (1990). Problems in cementing horizontal wells. *JPT*, 42, 398–400.**
- 核实：经 Physics of Fluids 34, 053610 (2022) 参考文献列表核实。
- 要点：水平井固井顶替问题的早期现场综述：低边窜槽=偏心+凝胶强度+失水造壁的组合。
- 回答：Q3。

**[A7] Lockyear, C.F., Hibbert, A.P. (1989). Integrated Primary Cementing Study Defines Key Factors for Field Success. *JPT*, 41, 1320–1325.**
- 核实：Energies 2021 综述参考文献列表（含卷期页）。
- 要点：现场综合研究：提高居中度对泥浆顶替与水泥放置质量有"显著（dramatic）"影响。
- 回答：Q1（现场方向性证据）、Q4。

**[A8] Lockyear, C.F., Ryan, D.F., Gunningham, M.M. (1990). Cement Channeling: How to Predict and Prevent. *SPE Drilling Engineering*, 5, 201–208.**（会议版 SPE-19865-MS，SPE ATCE San Antonio, 1989）
- 核实：OnePetro SPE Drill Eng 5(3) 201 页面 + OSTI 6405392（cited by 99）。
- 要点：45 ft 长 7-in 管/名义 8½-in 井眼大尺度流动环路实验。结论：单一流体速度剖面不能预测两相界面剖面；防窜槽需控制窄边界面速度（提密度差、紊流或周向交换）；并提出"水泥屈服应力过高时剪切应力可能永远无法超过屈服应力，水泥自身也可能无法进入窄边"——临界压力梯度思想的现场版。
- 回答：Q3（临界压力梯度）、Q2（正方：窄边剪切不足与泵量弱相关）。

### B 组：缝隙平均（2DGA/Hele-Shaw）模型谱系（Q5 主线）

**[B1] Bittleston, S.H., Ferguson, J., Frigaard, I.A. (2002). Mud removal and cement placement during primary cementing of an oil well — Laminar non-Newtonian displacements in an eccentric annular Hele-Shaw cell. *Journal of Engineering Mathematics*, 43(2–4), 229–253.**（AADE-18 参考文献引为 223–253，两处转引卷 43 一致）
- 核实：Springer 页面（link.springer.com/article/10.1023/A:1020370417367）摘要核实。
- 要点：**与本项目 1D-2D 耦合模型同源的方法奠基论文**：将窄环空 N-S 方程按缝隙平均化为二维浓度对流方程+准线性 Poisson 型流函数方程，把偏心环空映射为 Hele-Shaw cell；演示了窄边泥浆通道的再现、稳态位移解与真实工况复杂位移。速度场由局部槽流解组装（非全局 b³ 律）。
- 回答：Q5（缝隙平均方法奠基）、Q3（模型再现窄边通道）。

**[B2] Pelipenko, S., Frigaard, I.A. (2004) 三部曲：**
- [B2a] "Visco-plastic fluid displacements in near-vertical narrow eccentric annuli: prediction of travelling-wave solutions and interfacial instability." *J. Fluid Mech.*, 520, 343–377.
- [B2b] "Two-dimensional computational simulation of eccentric annular cementing displacements." *IMA J. Appl. Math.*, 69, 557–583.
- [B2c] "On steady state displacements in primary cementing of an oil well." *J. Eng. Math.*, 46, 1–26（全文 PDF：pelipenko.co.uk/files/JEM2004.pdf，本次直接打开核实）。
- 核实：Energies 2021 参考文献列表 + JFM 2023 正文转引 + B2c 全文 PDF 直接打开，三重核实。
- 要点：2DGA 框架数学化：行波解与界面失稳判据；[B2b] **确定性计算出偏心环空窄边静态泥浆通道**（JFM 2023 评价："definitively computed static channels of mud on the narrow side of eccentric annuli, as postulated many years earlier (McLean et al. 1967)"）；[B2c] 系统分析 2DGA 与早期水力学类模型的联系（JFM 2023 转述为"showing their conservative nature"；**该"保守性"指代方向未逐字核实，引用需回读原文**）。
- 回答：Q5（2DGA 数学基础）、Q3。

**[B3] Zhao, Y., Wang, Z., Zeng, Q., et al. (2016). Lattice Boltzmann simulation for steady displacement interface in cementing horizontal wells with eccentric annuli. *Journal of Petroleum Science & Engineering*, 145, 213–221.**
- 核实：ScienceDirect 摘要页 + PoF 2022 参考文献列表双重核实。
- 要点：中文团队基于 LBM+Hele-Shaw 方法模拟水平井偏心环空稳态位移界面——缝隙平均思路在中文固井界同样被视为正式手段。
- 回答：Q5。

**[B3b]《海上大位移井固井顶替数值模拟研究》. 《海洋工程装备与技术》, 2014, 1(1): 14–（起止页/作者未核实到一手页面，见弃用清单注记）.**
- 核实：期刊官网摘要页（qk.sjtu.edu.cn/oeet/article/2014/2095-7297/2095-7297-1-1-14.shtml）核实结论与期刊卷期；作者未核实。
- 要点：FLUENT 三维两相模拟（南海大位移井）：偏心度越大速度差异越明显、界面越不稳定、顶替效率越差，**建议将偏心度控制在 0.5 以内**；套管正弦屈曲时顶替效率 90.16%、螺旋屈曲 89.62%；60° 螺旋角旋流扶正器可提效。
- 回答：Q1（中文 CFD 定量）、Q4（中文 0.5 偏心度上限建议）。

### C 组：CFD 参数化与实验/现场定量（Q1/Q2 主线）

**[C1] Dai, H., Liu, G. (2018). An Overview of Annular Displacement Efficiency in Cementing Jobs Using an Efficient Numerical Model. AADE-18-FTCE-094, 2018 AADE Fluids Technical Conference and Exhibition, Houston, April 10–11, 2018.**
- 核实：AADE 官网 PDF 全文提取（aade.org/download_file/1335/394），本次直接下载全文。
- 要点：**与本项目架构最接近的工业级模型**：准静态伪 3D 缝隙平均+VOF（轴向速度用平行板 H-B 局部槽流解+压力迭代，周向速度由连续性恢复，压力假设横截面均一）。验证：与主流商业模拟器泥浆风险图一致；与 ANSYS 在 e=0.5、Re=176、τ0=1.22 Pa 流体下旋转参数化对比一致。
- 定量结论（均出自该文正文与表 4–7）：①0% standoff 无旋转→显著窄边泥浆通道；②30 rpm 旋转可完全消除 0% standoff 窜槽；③**70% standoff + 30 rpm 反而劣于低转速**（旋转把窄边低黏流体带往宽边→宽边壁面剪切下降→泥膜残留的反机制）；④**"完美居中也不给出 100% 顶替"**（原文：a perfectly centered casing will not give a 100% displacement——无滑移壁面泥膜）；⑤排量 4→12 bpm 在各 standoff 下效率均略升（所研究工况多属紊流区）；⑥YP 递增序列规则被其表 4 部分质疑（同流变 60/60 对顶最优）。
- 回答：Q1（CFD 参数化矩阵）、Q2（旋转/排量条件）、Q3（0% standoff+YP>0→通道的条件矩阵）、Q5（伪 3D 缝隙平均工业可行性）。

**[C2] Enayatpour, S., van Oort, E. (2017). Advanced Modeling of Cement Displacement Complexities. SPE-184702-MS, SPE/IADC Drilling Conference and Exhibition, The Hague, 14–16 March 2017. DOI 10.2118/184702-MS.**
- 核实：OnePetro proceedings 页面 + van Oort UT Austin 官方 CV 双重核实。
- 要点：3D 有限体积 VOF 模型系统参数化：**旋转+高偏心降低环空摩擦压降；旋转+低偏心改善水泥顶替效率（CDE）；旋转可部分抵消偏心的不利影响**；偏心使环空压力梯度下降（宽边短路）。定量数字在付费墙内，经 MDPI Energies 17:1226 转述核实其定性结论。
- 回答：Q1、Q2（旋转补偿条件）。

**[C2b] Guzman, J., Mavares, F. (2018). Casing Centralization and Pipe Movement in Cementing Operations for Improved Displacement Efficiency. SPE-191255-MS, SPE Trinidad and Tobago Section Energy Resources Conference, June 2018.**
- 核实：OnePetro 18TTCE 页面（cited by 21）+ ResearchGate 条目（DOI 10.2118/191255-MS）双重核实。
- 要点：现场+模拟："高 standoff 给出高顶替效率"；旋转 5–10 rpm、往复 2 ft/min 即可有效改善。
- 回答：Q1（现场方向性）、Q4。

**[C3] Yang, M., Yang, L., Li, Y., Zhang, L., Wang, C., Zhou, L., Yang, M. (2021). Improving displacement efficiency by optimizing pad fluid injection sequence during primary cementing of eccentric annulus in shale gas horizontal wells. *Journal of Petroleum Science & Engineering*, 204, 108691.**
- 核实：ScienceDirect S092041052100351X 摘要 + MDPI Energies 16:3650 参考文献列表（完整作者与卷页）双重核实。
- 要点：页岩气水平井偏心环空（e=0.4）CFD：三相顶替效率 85.88%；四相（泥浆-隔离液-冲洗液-水泥）89.67%；**把水泥浆用量提高到 1.5 倍环空体积后效率升至 95.16%**。密度差与浆柱结构优化显著减少水平段泥浆滞留。
- 回答：Q1（e=0.4 定量点）、Q2 反方（体积充足可补救）。

**[C3b] 孙劲飞, 李早元, 罗平亚, 张刚刚, 焦少卿 (2019). 水平井偏心环空低速顶替运移机制研究. 《西南石油大学学报(自然科学版)》, 41(1): 111–118.**
- 核实：期刊官网 HTML 全文页核实（含引用格式、作者单位、基金号；SWPU 油气藏地质及开发工程国家重点实验室）。
- 要点：SWPU 团队 CFD+VOF，水平井偏心环空低速顶替：**低偏心度下一次隔离液用量（1 倍环容）即可实现 90% 顶替效率；偏心度 >0.5 时效率相对较低，即使增加隔离液用量也无法进一步提高**（体积饱和现象）；流速 1.0→0.2 m/s 提效幅度：e=0.2 时仅 +1.8%，e=0.5 时 +6.8%，e=0.8 时 +4.2%（拐点在 e≈0.5）。
- 回答：Q1（中文 η-e-流速曲面）、Q2（低偏心 1 倍体积→90%：反方证据；e>0.5 加体积无效：正方/体积饱和证据——一文两侧）。

**[C4] 《低速注水泥时密度差对顶替效率影响规律的数值模拟研究》. 《石油钻采工艺》, 2008, 30(5): 62–65.**（作者未完全核实——PDF 为扫描件，见弃用清单注记；文章编号 1001-0890(2008)05-0062-04 核实）
- 要点：FLUENT 模拟+试验对比：低速顶替时正密度差效果最好、负密度差最差；偏心 0.5、0.2 m/s 工况下较大正密度差能克服偏心影响提高效率；偏心 0.5 时效率降至约 50% 量级。
- 回答：Q1（e=0.5 低返速定量点，采信度中等）。

### D 组：窄边窜槽机理与静止区（Q3 主线）

**[D1] Roustaei, A., Gosselin, A., Frigaard, I.A. (2015). Residual drilling mud during conditioning of uneven boreholes in primary cementing. Part 1: Rheology and geometry effects in non-inertial flows. *Journal of Non-Newtonian Fluid Mechanics*, 220, 87–98.**
- 核实：Frigaard 组托管全文 PDF（nrs.blob.core.windows.net）直接打开核实。
- 要点：不规则井眼+偏心环空非惯性位移模拟。**原句核实（本文档 Q3 的核心证据）**："the mud's yield stress determines the quality of the displacement on the narrow side to a large degree, e.g. without a yield stress the static mud channel cannot form"（泥浆屈服应力在很大程度上决定窄边顶替质量；没有屈服应力就不会形成窄边静态泥浆通道）；τy=10 Pa、e=0.3 时顶替"largely effective"；**"in presence of a sufficiently high yield stress, neither the use of centralizers nor moderate density difference can prevent the development of a mud channel"**（屈服应力足够高时，扶正器与适度密度差都无法阻止泥浆通道形成）；通道宽度跟随井眼不规则性变化，扶正器位置周期性打断通道。
- 回答：Q3（屈服应力是窄边通道的开关）、Q2 正方（高屈服下居中度与密度差都救不了）。

**[D1b] Skadsem, H.J., et al. (2019). Annular displacement in a highly inclined irregular wellbore. *Journal of Petroleum Science & Engineering*（ScienceDirect S0920410518307630，卷页经摘要页核实）.**
- 要点：实验+模拟：偏心促进宽边流动、界面轴向拉长；扩径段残留位置随偏心与井斜而变；向水平的偏转改善扩径段顶替。
- 回答：Q3（扩径+偏心组合）。

**[D2] MalekMohammadi, S., Carrasco-Teja, M., Storey, S., Frigaard, I.A. (2010). An experimental study of laminar displacement flows in narrow vertical eccentric annuli. *Journal of Fluid Mechanics*, 649, 371–398.**
- 核实：Cambridge Core 页面（cited by 62）。
- 要点：窄垂直偏心环空（强偏心）实验。**Q2 的实验裁定句**："即使强偏心下，合适的黏度比、密度比与流速组合仍能得到沿环空全长的稳态行波位移"；小偏心、高黏度比、高密度比、低流速有利于稳态位移（牛顿流体）；**"With a strong enough yield stress and with a large enough eccentricity, unyielded fluid remains behind on the narrow side of the annulus"（屈服应力足够强+偏心足够大→未屈服流体残留在窄边）**；实验与 Hele-Shaw 型模型（Bittleston 2002）预测定性一致。
- 回答：Q2（两派裁定的实验句）、Q3、Q5（Hele-Shaw 模型的实验支持）。

**[D2b] Carrasco-Teja, M., Frigaard, I.A. (2009/2010).** ①"Displacement flows in horizontal, narrow, eccentric annuli with a moving inner cylinder." *Physics of Fluids*, 21, 073102 (2009)；②"Non-Newtonian fluid displacements in horizontal narrow eccentric annuli: effects of slow motion of the inner cylinder." *J. Fluid Mech.*, 653, 137–173 (2010).
- 核实：均经 AIP PoF 34, 053610 (2022) 参考文献列表核实。
- 要点：慢旋转内筒的非牛顿位移：周向流把驱替液带入窄边，旋转缩短窄边残留层长度。
- 回答：Q2（旋转补偿机制）。

**[D2c] （Frigaard 组）(2022). Displacement flows in eccentric annuli with a rotating inner cylinder. *Physics of Fluids*, 34(5), 053610, DOI 10.1063/5.0092026.**（第一作者姓名未核实，作者含 Frigaard/UBC 组——引用时以"（Frigaard 组，2022）"表述，见弃用清单注记）
- 核实：AIP 页面（标题/卷期页/DOI/摘要）核实。
- 要点：牛顿流体水平偏心环空+内筒旋转实验：高偏心时发生"流动分离"（驱替液被限制在宽边，有效顶替延迟）；**"rotation improves the displacement (volumetric efficiency)...and can prevent a narrow mud channel from forming if the excess fluid volume is used"（若使用过量流体体积，旋转可阻止窄边泥浆通道形成）**——Q2 反方的最强证据（体积充足+旋转→窄边可清）。
- 回答：Q2 反方核心、Q1。

**[D3] Zhang, R., Frigaard, I.A. (2023). Primary cementing of vertical wells: displacement and dispersion effects in narrow eccentric annuli. Part 2. Flow behaviour and classification. *Journal of Fluid Mechanics*, 972, A38.**
- 核实：Cambridge Core 页面 + 独立论文引用列表双重核实（卷 972, A38）。
- 要点：牛顿流体窄垂直偏心环空实验+**2DGA 与全 3D 计算正面对比**：两模型预测出相同的稳态/非稳态位移与残余流体位置，窄边残余两者都有；"大偏心对完全高效顶替不利"；周向二次流与弥散决定周向分布均匀性；并转述"即使紊流中冲洗液仍趋于沿宽边窜流（Guillot et al. 2007; Maleki & Frigaard 2018）"。正文亦确认："The current popular way to model this process is via a Hele-Shaw approach"——Hele-Shaw 缝隙平均是当前主流建模方式。
- 回答：Q5（2DGA vs 全 3D 一致性与差异）、Q1。

**[D3b] Zhang, R., Ghorbani, M., Wong, S., Frigaard, I.A. (2023). Vertical cementing displacement flows of shear-thinning fluids. *Physics of Fluids*, 35（卷号经 Semantic Scholar 页面核实，页码未核）.**
- 要点：幂律（剪切变稀、无屈服）流体垂直偏心环空实验+计算。**与 Bingham 的分野：无屈服→窄边是"慢流"而非"停流"，通道被慢流主导**——这正是 b³/b^(n+2) 缝隙律（无静止区）与屈服流体局部槽流解（有停流判据）的分水岭。
- 回答：Q5（b³ 与屈服流体静态区的分水岭）、Q3。

### E 组：湍流与流态证据（Q2/Q3）

**[E1] Maleki, A., Frigaard, I.A. (2018/2019).** ①"Turbulent displacement flows in primary cementing of oil and gas wells." *Physics of Fluids*, 30(12), 123101 (2018)；②"Comparing laminar and turbulent primary cementing flows." *Journal of Petroleum Science & Engineering*, 177, 808–821 (2019)（cited by 38，ScienceDirect 摘要页核实）。
- 要点：②的 highlight 原句核实：**"Provides concrete evidence suggesting that turbulent displacement flows do not necessarily outperform laminar displacement flows"**（给出具体证据表明湍流顶替流不必然优于层流）；关键在窄边速度与流态匹配；①+转述：低黏冲洗液即使紊流也趋于沿宽边窜流（与 Guillot et al. 2007 一致）。
- 回答：Q2/Q3（湍流信条的条件化）、Q5（2DGA 湍流闭合已存在）。

**[E1b] Guillot, M., Desroches, J., Frigaard, I.A. (2007).**（经 JFM 2023 [D3] 正文转引核实，一手未核，见弃用清单注记）冲洗液紊流仍沿宽边窜流。
- 回答：Q2/Q3。

**[E2] Ytrehus, J.D., Lund, B., Taghipour, A., Divyankar, S., Saasen, A. (2017). Experimental Investigation of Wellbore Fluid Displacement in Concentric and Eccentric Annulus. Proc. ASME OMAE 2017.**（经 Energies 2021 综述参考文献列表核实）——同心/偏心环空顶替实验谱系；相关：Lund, B., Taghipour, A., Ytrehus, J.D., Saasen, A. (2020). *Energies*, 13(19)（实验方法）。
- 回答：Q1（实验谱系）。

### F 组：现场/工业阈值与标准（Q4 主线）

**[F1] API RP 10D-2, 2nd ed. (2023). Centralizer Placement and Stop-Collar Testing. API, 58 页；ISO 对应标准 ISO 10427-2:2004（Centralizer placement and stop-collar testing）.**
- 核实：accuristech + globalspec + iTeh-standards 三处商店/目录页核实（卷期、页数、年份）；ISO 10427-2 摘要含 standoff ratio 定义 R = (standoff/annular clearance)×100。
- 要点：给出斜井/狗腿度井眼扶正器间距计算方法（基于扶正器性能与目标 standoff）与停环试验程序。**是方法标准，非数字阈值出处；本次未在其公开摘要层面核实到 67%/70%/75% 的数字。**
- 回答：Q4（standoff 计算与扶正器布置的方法标准出处）。

**[F1b] API Spec 10D（扶正器规格标准，老版）——经 OGJ 文章转述核实**："The old standard of 67% standoff, as formerly recommended by the American Petroleum Institute (API), was considered adequate for most wells. However, the current API specification 10D states that '67% may or may not be sufficient for a good cementation.'"（老 API 推荐 67% standoff；现行 API Spec 10D 明确指出 67% 对良好固井"可能够也可能不够"。）**67% 数字本身未在标准文本层面核实，属工业转述层。**
- 回答：Q4（67% 的"老 API 推荐"地位）。

**[F2] Oil & Gas Journal. "HORIZONTAL WELLS-7: Cementing horizontal holes becoming more common."**
- 核实：ogj.com 页面核实。
- 要点（Q4+Q2 双重）：**"To reduce low-side bypass and top-side channeling, standoff should be 67% or greater"**（减少低边绕流与顶部窜槽，standoff 应≥67%）；**"Without proper standoff, even very high annular flow rate may not develop enough wall shear stress to overcome mud gel strength on the low side of poorly centralized casing"**（无适当居中时，即使很高的环空排量也可能产生不足以克服低边泥浆凝胶强度的壁面剪切应力）——Q2 正方的工业口径。
- 回答：Q4（67% 工业惯例出处）、Q2 正方。

**[F2b] Oil & Gas Journal. "Proper centralizers can improve horizontal well cementing."**
- 核实：ogj.com 页面核实。
- 要点：实验室+现场经验：**standoff 80–90% 被认为足够（adequate），即使对水平井**；实际井下居中度受许多地区性因素影响，数学模型只描述理想情形。
- 回答：Q4（80–90% adequate 口径）。

**[F3] Drilling Contractor 文章（转述 Moroni, N., Ravi, K., Hemphill, T., Sairam, P. 2009, SPE-124726-MS "Pipe Rotation Improves Hole Cleaning and Cement-Slurry Placement: Mathematical Modeling and Field Validation", SPE Offshore Europe, Aberdeen）.**
- 核实：drillingcontractor.org 页面 + AADE-18 参考文献列表双重核实。
- 要点：偏心套管使隔离液/水泥沿阻力最小路径走（宽边短路），低边留下大段泥浆段；旋转提供周向输运清除窄边。
- 回答：Q4（工业口径）、Q2（旋转补偿）。

**[F4] 张晋凯, 周仕明, 陶谦, 方春飞, 张林海, 高元 (2016). 套管低偏心度下的水泥浆顶替界面特性研究. 《石油机械》, 44(7): 1–6. DOI 10.16082/j.cnki.issn.1001-4578.2016.07.001.**
- 核实：期刊官网 HTML 全文+DOI 核实（rhhz.net/syjxzz/20160701.htm；中石化石油工程技术研究院+中国石油大学(北京)）。
- 要点（本次调研最有用的中文定量矩阵）：1000 m 水平段、600 万网格 CFD，偏心度 e=0/0.1/0.2/0.3 × 正密度差 Δρ=100/300 kg/m³，**1000 s 顶替效率：Δρ=300 时 0.86/0.88/0.84/0.83；Δρ=100 时 0.84/0.85/0.80/0.81**。三个关键结论：①**低偏心段（e≤0.3）效率差仅 2–8 个百分点**；②e=0.1 的效率反而高于完全居中（0.88 vs 0.86）——存在"**临界最优偏心度**"（0.1–0.2 之间），界面长度先减后增、效率先增后减（浮力-偏心耦合平衡）；③偏心效应的影响大于密度差；e 超过临界最优偏心度后影响程度明显更大。
- 口径警示：该表效率为 1000 s 顶替进行中的瞬时口径，非终态效率；与全井 η_E 不可直接对表，但"低偏心段敏感性量级（2–8 pp）"的结论可用。
- 回答：Q1（低偏心段定量矩阵）、Q2 反方（低偏心下居中度弱敏感甚至非单调）、Q4（中文"临界最优偏心度"口径）。

**[F4b] "Research on Key Technologies to Improve Cementing Displacement Efficiency." ACS Omega/PMC（PMC9607674，2022；期刊名/DOI 待核，见弃用清单注记）.**
- 要点：中文团队滞留层模型：偏心造成宽窄间隙→窄边流动阻力大→延迟流动甚至整体不流；其模型口径下顶替效率随偏心度增大**单调下降**（与张晋凯非单调结论矛盾，条件不同：本文为简化模型+油基泥浆口径）；偏心 10% 时动切力越高效率越低（YP 反直觉效应）；密度差 0.2 g/cm³ 时效率可达 90%。
- 回答：Q1、Q2（单调性矛盾见对照表）。

### G 组：综述

**[G1] Energies 2021, 14, 573. "Cement Placement: An Overview of Fluid Displacement Techniques and Modelling."**（33 页综述；MDPI 页面+LAPSE-2023.29256 全文 PDF 双重核实）
- 要点：Hele-Shaw/2DGA 谱系地图（Bittleston 2002 → Pelipenko-Frigaard 2004 三部曲）；Tehrani 10–15% 密度差建议；Couturier et al. "只有处处湍流才有效"；McLean"水泥变稀追求紊流反而加重窜槽"的警告；Macondo 窄边泥浆通道实例（Chief Counsel 报告图 1）；**扶正器间套管下垂（casing sag between centralizers）使实际 standoff 远低于设计值**。
- 回答：全部问题的谱系地图与"实际 standoff≠设计 standoff"警示。

**[G2] "Cement Placement Modeling—A Review." *SPE Drilling & Completion*, 38(02), 297（2025）.**（OnePetro DC 38(02) 297 页面核实；作者未核，见弃用清单注记）——最新综述，转引 Maleki Zamenjani (2018) 层流/紊流/混合流态模拟。
- 回答：最新综述背书。

---

## ② Q1–Q5 综合回答

### Q1：standoff/居中度对顶替效率的典型影响量级

**结论：量级高度条件依赖，呈"低偏心不敏感—高偏心断崖"的阈值型形态；没有普适 η-standoff 曲线本身就是文献共识**（量级由窄边屈服门决定，standoff 只是 b_narrow 的输入之一）。分区间汇总（偏心度 e，居中度=1−e）：

| 偏心度区间 | 文献定量结果 | 出处 |
|---|---|---|
| e ≤ 0.2（居中度≥80%） | 效率差仅 **2–8 pp**；e=0.1 甚至优于完全居中（0.88 vs 0.86，Δρ=300 kg/m³） | [F4] 张晋凯 2016 表 2/表 3（瞬时口径） |
| 0.2 < e ≤ 0.5 | 明显恶化但仍可接受：e=0.4 → 85.88%（三相）/89.67%（四相）；e=0.5 低返速时约 50% 量级（待核条目 [C4]）；正弦屈曲 90.16%/螺旋屈曲 89.62% | [C3] Yang 2021；[C4]；[B3b] |
| e > 0.5 | 断崖：效率显著下降且**加隔离液量也无效（体积饱和）** | [C3b] 孙劲飞 2019；[C4] 元坝口径"近似指数衰减"（注：元坝条目已并入 [C4] 同源扫描件组，采信度中等） |
| 旋转补偿 | 0% standoff + 30 rpm 可完全消除窜槽；e 高时旋转可部分抵消偏心不利；但 70% standoff + 30 rpm 反而劣于低转速（反机制） | [C1] Dai & Liu 2018 表 5；[C2] Enayatpour 2017；[D2c] PoF 2022 |
| 居中度 100% | 也不保证 100% 效率（无滑移壁面泥膜） | [C1] Dai & Liu 2018 |

有无 η-standoff 定量曲线：**有，但全部是条件化的**——张晋凯 2016 表 2（e 0–0.3 × Δρ 矩阵）、孙劲飞 2019 图 5（e 0.2/0.5/0.8 × 流速）、Yang 2021（e=0.4 单点+体积扫描）、Dai & Liu 2018 表 5–7（standoff 0/40/70/100% × 转速/往复/排量矩阵）。没有文献声称给出普适曲线。

口径警示：文献效率口径三分——①全井体积效率（我们的 η_E 口径）；②瞬时端面/截面效率（张晋凯 1000 s、元坝口径）；③窄边残留层厚度/通道宽度（Frigaard 组口径）。三类数字不可直接对表。

### Q2：体积充足（≥环空体积）时偏心是否仍导致显著效率损失——两派对照

**两派都有代表文献；分歧根源不在"体积"而在"窄边局部压力梯度 vs 屈服应力"。**

| 立场 | 代表文献 | 核心证据 |
|---|---|---|
| **正方：屈服卡死窄边→体积救不了** | [A1] McLean 1967；[A2] Walton-Bittleston 1991 停流判据；[D1] Roustaei 2015（"屈服应力足够高时，扶正器+适度密度差都无法阻止窄边通道"）；[A8] Lockyear 1990（水泥屈服过高→水泥自身也进不了窄边）；[A4] Haut-Crook 1979（凝胶强度）；[F2] OGJ（"无适当居中时高排量也不足以克服低边凝胶强度"）；[C3b] 孙劲飞 2019（e>0.5 加隔离液无效） | 屈服门关闭时窄边流量→0；静态通道是**稳态解**（[B2b] Pelipenko-Frigaard 确定性算出），多泵水泥只延长通道而不清除它。判据：τ_w=(dp/dz)·b_narrow/2 < τy → 窄边停流，**与累计体积无关** |
| **反方：体积充足+正确流变/旋转/密度差→窄边可清** | [D2] MalekMohammadi 2010（强偏心下稳态行波位移可达）；[C3] Yang 2021（1.5 倍环容→95.16%）；[D2c] PoF 2022（旋转+过量体积→窄边通道可避免）；[C2] Enayatpour 2017（旋转部分抵消偏心）；[C1] Dai & Liu 2018（30 rpm 消除 0% standoff 窜槽）；[C3b] 孙劲飞 2019（低偏心 1 倍隔离液→90%） | 窄边压力梯度>屈服门时窄边仍有慢流，多泵的每一桶都在前推窄边残留；旋转/密度差抬高窄边有效驱替能力；稳态行波位移下窄边被持续外输 |
| **统一框架** | [A2] 判据 + [D1] + [D2] | "体积充足"是**必要非充分**条件。真正控制变量：τy（泥浆屈服/凝胶强度）、b_narrow、dp/dz（受排量与宽边短路程度控制）、Δρ、旋转。**窄边停流 ⟺ (dp/dz)·b_narrow/2 < τy** |

**裁定**："体积充足→效率趋 100%"只在窄边屈服门开启时成立；屈服门关闭时多泵水泥不改变窄边停流态（稳态通道解），效率封顶。这同时解释了文献中"低偏心 1 倍体积→90%"与"高偏心加体积无效"两个看似矛盾的现象。

### Q3：窄边窜槽机理——谁定量描述、典型条件

1. **理论判据**：[A2] Walton-Bittleston 1991——窄边停流判据与静止区范围界：**τ_w = (dp/dz)·b/2，停流 ⟺ (dp/dz)·b_narrow/2 < τy**。这是 Q3 的理论基石。
2. **模型确定性再现**：[B2b] Pelipenko-Frigaard 2004（IMA JAM 69）确定性算出窄边静态通道；[B1] Bittleston 2002 演示通道形成。
3. **实验证实**：[D2] MalekMohammadi 2010："屈服应力足够强+偏心足够大→未屈服流体残留在窄边"；且强偏心+弱屈服+合适黏度/密度比→稳态行波位移（不残留）。
4. **开关变量**：[D1] Roustaei 2015——无屈服应力→通道根本不能形成；屈服应力是窄边通道的开关。
5. **典型条件量级**：
   - τy=10 Pa、e=0.3、非惯性流 → 窄边顶替"基本有效"（[D1] 图 3a,b）；"足够高屈服应力"→静态通道（[D1] 未给单一临界值）。
   - [C1] Dai & Liu 2018 条件矩阵：0% standoff + YP>0 → 通道；100% 居中+流变失配 → 泥膜/泥纹（非通道）；0% standoff+30 rpm → 通道消除。
   - [A4] Haut-Crook 1979：凝胶强度（静态）>动切力（动态）；停泵后凝胶强度接管，静止泥浆更难动。
   - [A8] Lockyear 1990：水泥屈服过高时水泥自身也进不了窄边（驱替液侧也有屈服门）。
6. **补偿机制及其边界**：密度差 10–15%（[A3]，层流建议）；"处处湍流"才有效（Couturier，经 [G1] 转述）；旋转周向输运（[D2b][D2c][F3]），但旋转存在反机制（[C1] 70% standoff+30 rpm 更差）。
7. **体积饱和**：e>0.5 时加隔离液无效（[C3b]）——稳态通道的必然推论：通道一旦形成，多泵的每一桶都在通道两侧流动，通道内静止泥浆永不被推走。

### Q4：工业阈值 67%/70%/75% 的出处与适用前提

**核实状态：67% 是"老 API 推荐"，其标准文本层的原始出处本次未能核实到一手文本（付费墙/未数字化），只能核实到工业转述层**：
- [F1b] OGJ 转述：老 API 推荐 67% standoff；现行 API Spec 10D 指出"67% may or may not be sufficient for a good cementation"。
- [F2] OGJ：为减少低边绕流与顶部窜槽，standoff 应≥67%；大多数服务公司软件可优化扶正器设计。
- [F2b] OGJ：实验室+现场经验，standoff 80–90% 被认为足够，即使水平井。
- [F1] API RP 10D-2（2023 第 2 版）/ISO 10427-2:2004：standoff ratio 定义与扶正器布置计算的方法标准（非数字阈值）。
- 中文口径：偏心度≤0.5 建议（[B3b] 待核条目）；"临界最优偏心度 0.1–0.2"（[F4]）。

换算与自洽性检查：standoff 67% ⟺ e≤0.33（线性换算）。与文献定量结果不矛盾：e≤0.3 处于"低敏感区"（效率差 2–8 pp，[F4]），e>0.5 进入断崖区（[C3b]）——**67% 是"充分不必要"的设计下限**：位于低敏感区上缘之外但仍在断崖之前。

适用前提（文献一致的适用前提）：①窄边压力梯度>屈服门（排除静态通道态）——无此前提则 67% 也救不了（[D1][F2]）；②密度差为正、流变匹配（或至少不倒挂，[C1] 表 4）；③紊流设计时须"处处湍流"（Couturier，经 [G1]）；④斜井/水平井另有低边重力滞留机制叠加，67% 不再充分（[A5][A6]）；⑤实际 standoff 受扶正器间套管下垂影响远低于设计值（[G1]），按设计值套用 67% 规则会高估井筒状态。

### Q5：缝隙平均模型能否表现居中度效应；b³ 的适用性与批评

**缝隙平均（2DGA/Hele-Shaw）方法是本领域公认的正式方法**，且有正面背书：
- 奠基与数学化：[B1] Bittleston 2002；[B2] Pelipenko-Frigaard 2004 三部曲（行波解、静态通道确定性计算、与早期水力学模型的联系分析）。
- 实验/全 3D 背书：[D2] MalekMohammadi 2010（实验与 Bittleston 2002 Hele-Shaw 预测定性一致）；[D3] Zhang-Frigaard 2023 JFM 972:A38（2DGA 与全 3D 预测相同的稳态/非稳态位移与残余位置）；[B3] Zhao 2016（LBM+Hele-Shaw）；[C1] Dai & Liu 2018（工业级伪 3D 缝隙平均与商业模拟器/ANSYS 对标一致）。JFM 2023 原文确认"Hele-Shaw approach 是当前主流建模方式"。
- 2DGA 已确定性算出窄边静态通道（[B2b]）——**"缝隙平均模型能否表现居中度效应"的答案是：能，且是文献主流手段**。

批评/局限谱系（文献中的边界，非"缺陷"）：
1. **纯 b³（牛顿缝隙律）无法产生静止区**：牛顿/幂律缝隙律（b³、b^(n+2)/(n+2)）下窄边速度慢但非零，效率损失平缓（[D3b]）；**屈服流体的局部槽流解才有停流判据**（[A2]）。因此"模型能否表现居中度效应"的关键**不在缝隙平均本身，而在流变模型是否含屈服项、局部槽流解是否带屈服门**。
2. 缝隙平均的适用前提：b/H≪1（环空窄缝）、缓变截面、层流。**强偏心 e→1（贴边）时窄边 b→0，lubrication 极限本身发散**——贴边情形缝隙平均失效（理论推断，文献未逐字核实此句，标注为理论边界而非引文）。
3. 周向二次流与弥散：2DGA 中由连续性诱导的周向速度分量承载，弥散需闭合参数化；[D2] 实验显示周向逆向流与弥散尖峰是真实存在的物理，2DGA 对其的再现是近似层面（[D3] 显示 2DGA 与 3D 在"何处残留"上仍一致）。
4. 湍流：2DGA 湍流闭合已存在（[E1] Maleki-Frigaard 2018/2019），不缺方法。
5. "conservative nature"句（JFM 2023 转述 [B2c]）：2DGA 与早期水力学模型联系"显示其保守性"——**指代方向未逐字核实**，引用需回读 [B2c] 原文。

**对我们项目内部口径的修正**：项目记忆中"b³ 分流高偏放大"的表述，按文献应修正为——局部槽流解（含屈服门）是正解；全局 b³ 律只是牛顿情形的推论而非普适律；若模型以纯 b³ 权重做周向分流且不含屈服门，才会抹平静止区/错误分配窄边流量。牛顿情形下局部槽流解对偏心几何逐点成立，b³ 流量重分配本身不是误差源；误差源是屈服门缺失或被正则化抹平。

---

## ③ 对"我们模型 η_E≈100% + 居中度不敏感"的文献裁定

### 裁定主结论：**取决于条件；在多数文献条件下属正常物理，但必须通过三个检验**

文献共识：居中度对顶替效率的影响由**窄边屈服门（narrow-gap yield gate）**控制——窄边局部壁面剪切应力 vs 泥浆屈服应力/凝胶强度，即 τ_w=(dp/dz)·b_narrow/2 与 τy 的比较。τy、b_narrow、dp/dz、Δρ、旋转共同决定门的状态；standoff 只是 b_narrow 的一个输入。文献中"效率高+居中度弱敏感"与"窄边窜槽"是**同一判据的两个分支**，可以在同一模型内共存（[C1] Dai & Liu 2018 一文就同时给出：100% 居中不保证 100% 效率 + 0% standoff 旋转可清通道）。

### 三个必做检验（把"正常物理"与"模型缺陷"分开）

**检验 1：窄边屈服门检验（决定性）**
- 判据：窄边单元速度归零 ⟺ (dp/dz)·b_narrow/2 < τy。
- 对 8 井逐井、逐段计算窄边屈服门状态（输入：泥浆动切力 τy、窄边缝宽 b_narrow、顶替期压力梯度 dp/dz）。
- 7 井门开（窄边流）→ η_E≈100% 属正常物理；唯一 e=0.547 井门关（窄边停流）→ 窄边窜槽是模型对判据越过阈值的正确表现。
- τy 取值：顶替期用动切力（动 YP）；停泵/候凝期用凝胶强度（我们工艺=胶塞驱动+连续顶替→用动 YP）。**若模型只用塑性黏度、无屈服项→居中度不敏感就是模型缺陷（进入分支 2）**。
- 附加失真源：若模型有屈服门但被正则化/弥散闭包抹平（正则化 τ→0 极限、c_min 冻结问题），窄边慢流会被抹平——文献判定为模型侧失真，需检查（对应我们项目 A4 报告的"b³ 非层流修正+c_min 解冻"方向，文献一致支持屈服门+局部槽流解是正解路径）。

**检验 2：窄边分流权重检验（b³ 权重 vs 局部槽流解）**
- 若模型用纯 b³ 权重分配环空周向分流且不含屈服门，会高估窄边分流（窄边更快被扫过→居中度不敏感+高效率）——即"模型缺陷分支"的一种实现。
- 文献侧修复路径：改用 Walton-Bittleston 式局部槽流解（含屈服门）作为分流权重。

**检验 3：口径检验**
- η_E 是体积效率口径，对窄边通道**部分稀释**：若窄边残留占环空体积 1–2%，η_E 仍可 98%+。**"η_E≈100%"与"窄边窜槽"不矛盾**（Dai & Liu 2018 表 5 中 0% standoff 无旋转时通道显著、效率显著<100%；较高 standoff 时效率高）。
- 我们 8 井：7 井居中度较高+屈服门开→效率≈100% 正常；hu103 e=0.547 屈服门关→窄边残留、η_E<100%——**恰好落进文献的"e>0.5 断崖区"**。
- 风险：体积效率口径会把小体积窄边窜槽稀释进"≈100%"里（与项目记忆 F4"92% 入库稀释"同构）。建议并行报告窄边局部效率或残留层厚度口径（Frigaard 组口径）。

### 裁定三分支判定树

- **分支 1（最可能）：正常物理。** 若检验 1 通过（7 井门开、hu103 门关），文献裁定：η_E≈100% + 居中度不敏感属正常物理，与文献"低偏心段效率差仅 2–8 pp（[F4]）""湍流不必然优于层流（[E1]）""低偏心 1 倍体积→90%（[C3b]）"完全一致；唯一例外井窜槽=判据越阈的正确表现，与"e>0.5 断崖（[C3b]）"一致。支持文献：[F4][C3b][E1][D2][C1]。
- **分支 2：模型缺陷。** 若检验 1/2 发现模型无屈服门或屈服门被正则化抹平，居中度不敏感=模型缺陷。修复路径：局部槽流解+屈服门（[A2][B1]）；c_min 解冻；弥散闭包量级校准。（注：B1 通量修复后 7/8 井效率 0.85–0.999 随居中度单调——井间单调与张晋凯"参数扫描下临界最优偏心度"非单调不矛盾：井间差异来自多参数组合，不是单变量扫描。）
- **分支 3：口径问题（不改模型，改报告）。** 并行报告三种口径：①全井体积 η_E；②窄边局部效率/残留层厚度；③standoff×YP 条件矩阵（仿 [C1] 表 5）。

### 给模型团队的行动建议（文献侧）

1. **逐井窄边屈服门状态表**：输入动 YP、b_narrow、dp/dz，判定每井每段 τ_w=(dp/dz)·b_narrow/2 是否低于 τy；hu103 重点复核。
2. **并行口径报告**（分支 3）：全井 η_E + 窄边局部效率/残留层厚度 + standoff×YP 矩阵。
3. **局部槽流解权重检查**（分支 2 排查）：分流权重是否含屈服门；若纯 b³ 权重无屈服门→按文献为模型缺陷，改局部槽流解权重。
4. **e>0.5 井（hu103）单独报告窄边通道宽度/残留层厚度**，与 [C1] 表 5 的 0% standoff 通道形态对照。
5. **居中度口径警示**：扶正器间套管下垂（sag between centralizers，[G1]）使实际 standoff 远低于设计值；若 8 井居中度是设计/平均口径，真实最小缝隙段居中度更低，屈服门状态会恶化。**平均居中度会低估窄边风险**——建议以"最小缝隙段居中度"复核 hu103。
6. **τy 数据质量**：文献反复强调凝胶强度（静态）>动切力（动态）（[A4][A1]）；顶替期用动 YP，但若井上有停泵等待（我们工艺为胶塞连续顶替，不适用），须用凝胶强度复核。

### 本报告局限声明

- OnePetro/SPE 付费墙：Enayatpour 2017、Guzman 2018 仅核到摘要/转引层（卷期页+DOI+被引数核实；定量结论经 Energies 17:1226 转述核实）；Couturier et al. 仅转引核实。
- 中文条目核实程度不一：张晋凯 2016、孙劲飞 2019 全文核实；2008 石油钻采工艺与元坝组条目为扫描件 OCR 受损，作者/卷期未完全核实（采信度中等）；PMC9607674 期刊名待核；[B3b] 作者待核。
- 张晋凯表 2/3 效率为 1000 s 瞬时口径，与全井终态 η_E 口径不同。
- "conservative nature"句方向未逐字核实；PoF 2022 第一作者未核实到姓名；[D3b] 卷号页码未全核。
- 未覆盖维度：温度、失水/泥饼、防气窜、U 型管效应（超出"文献 A：居中度维度"范围）。

### 弃用清单（无法核实/仅转引，引用前须回读一手）

| 条目 | 弃用/降级原因 |
|---|---|
| "Displacement Efficiency in Eccentric Annuli"（ResearchGate 347725348，疑似 SPE 论文） | 检索未能核实到作者/会议/卷期页（检索词命中词表文件污染），无法核实 |
| 2008《石油钻采工艺》低速注水泥密度差文章（syzt 109.pdf） | 扫描件 OCR 乱码，作者名未核实；结论已标注"采信度中等" |
| 元坝地区尾管固井数值模拟（《石油钻探技术》2010(4)，DOI 10.3969/j.issn.1001-0890.2010.04.011） | PDF OCR 受损，作者名未核实；"近似指数衰减"结论采信度中等 |
| PMC9607674（"Research on Key Technologies..."） | 期刊名/DOI 未在 PMC 页面直接核实，标注"待核"使用 |
| 《海上大位移井固井顶替数值模拟研究》（海洋工程装备与技术 2014, 1(1)） | 期刊官网核实到卷期与结论，作者名未核实 |
| Couturier, Guillot, Hendriks, Ropel（SPE 26582?） | 仅经 Energies 2021 与 JFM 2023 转引核实，一手未核 |
| Guillot, Desroches, Frigaard (2007) | 仅经 JFM 2023 转引核实 |
| armarrorocks 859.pdf | 身份不明的论文 PDF，仅作引用索引，不作为独立证据条目 |
| API Spec 10D / API RP 65 中"67%"数字 | 数字本身未在标准文本层面核实（只核到 OGJ/服务公司转述层），引用必须注明"工业转述层" |
| PoF 34(5), 053610 (2022) 第一作者姓名 | 仅核实到作者含 Frigaard（UBC 组）；引用以"（Frigaard 组，2022）"表述 |
| Parker et al. 1965 SPE-1234 | 仅经第三方引用列表转引核实 |

---

*报告完。所有条目经 tavily search/extract 在线核实（2026-09-03）。*
