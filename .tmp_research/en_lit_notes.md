# 英文文献调研笔记：固井"插旗杆"与"灌香肠"对应术语与文献

> 调研日期：2026-09-02。方法：Tavily 网络搜索（26 轮关键词检索 + 3 次页面全文提取），来源以 SPE/OnePetro、ScienceDirect、Springer、MDPI、OSTI、PMC 及行业技术站点为主。**所有条目均来自实际搜索/提取结果，未编造**；未能核实的信息已标注。撰写语言：中文，引用保留英文原文。

---

## 1. 术语对照表（中文俗语 ↔ 英文正式术语 ↔ 定义）

| 中文俗语 | 英文正式术语 | 英文文献中的定义/表述 | 备注 |
|---|---|---|---|
| **插旗杆**（水泥浆在管柱内/环空提前凝固把管柱"焊"死） | **stuck pipe due to cement / cemented-up drill pipe / green cement sticking / flash set** | "A more serious situation is that of cement setting around the drill pipe. This can be caused by: The flash setting of the cement when setting a cement plug. A leak in a string while setting a plug or cementing a string of casing through a drill pipe and a stabbed-in stinger. Human error, such as miscalculation of displacement time or failure to add CMT retarders when setting a plug."（Drilling Manual, Stuck Pipe Course） | 英文无"插旗杆"直接对应词；工程语境按成因拆为 flash set（闪凝）、premature setting（提前凝固）、green cement（未凝软水泥）粘卡、cement blocks（水泥块卡钻） |
| ——（闪凝） | **flash set** | "Flash Set is the state of cement slurry when became too gelatinous near to stiffening point during mixing and/or pumping."（M. Mohammed, LinkedIn, Oil and Gas Well Cementing Part 3 – Flash Set and False Set）；注意与 **false set**（假凝，可搅拌恢复）区分 | 高温深井主要风险；>120°C 时高异常胶凝风险显著（Springer 2025 综述） |
| ——（未凝软水泥粘卡） | **green cement stuck / soft cement** | "the clean-out or drilling assembly is run after spotting a cement plug; the drill string may encounter soft (green) cement, becoming trapped and stuck."（Drillopedia, Mechanical sticking） | 处理：上提+震击；凝固后需倒扣/套铣 |
| ——（水泥块卡钻） | **stuck in cement blocks** | "The cement around the shoe of a string of casing pipe can break up due to the impact of the drill pipe against the casing during subsequent drilling. It can then fall into the hole in blocks and jam the bit or on top of a stabilizer."（Drilling Manual） | 下回次钻穿未凝水泥鞋时 BHA 激动压力可诱发 flash set |
| **灌香肠**（水泥浆未顶出管内/管内留水泥柱） | **cement left in casing / under-displacement / casing full of cement** | "Brainstorm time: Cement left in Pipe"（r/oilandgasworkers, 2024，公司级事故调查案例）；"contaminated cement to be trapped within casing"（Innovex 浮箍产品说明） | 英文口语也有 "cement in the shoe track" 表述；正式术语多用 under-displacement |
| ——（浮箍失效回流留塞） | **cement fallback / reverse flow / float valve failure** | "Failure of the float to operate properly will allow the cement to flow back into the casing. If this happens it can be countered by holding back[-pressure]."（Drilling Manual, Float Shoe And Float Collar） | 直接对应"灌香肠"最常见成因之一 |
| ——（湿鞋跟） | **wet shoe track** | "A 'Wet' shoe track is defined as the occurrence of unset, contaminated or no cement in the casing section between float collar and shoe after a primary cement [job]."（Hibbeler, SPE-62751-MS, 2000） | 与 dry shoe track（水泥已凝固的鞋跟）相对 |
| ——（管内留塞打塞时） | **cement plug stuck / stinger stuck / plug back-off** | "There is a much higher likelihood of fluid U-tubing out of the drillpipe, resulting in spills on the rig floor... the stinger is abandoned in the well."（JPT, Stinger or Tailpipe Placement of Cement Plugs, 对 SPE-168005 的评述） | 常见后果：塞面以上钻杆被凝、stinger 报废井内 |
| ——（环空水泥返高不连续/窜槽） | **mud channel / cement channeling / discontinuous cement column** | "This early acceleration causes a 'partially empty' or 'discontinuous' zone or gap to form between the free-falling column of fluids and the wellhead."（Beirute, SPE-13045-MS, 1984——注意此处"discontinuous"指管内自由下落段，非环空）；环空窜槽表述为 "mud channel at the narrow side"、"incomplete mud removal"（AADE-18-FTCE-094） | 环空"分段水泥柱"对应文献多用 mud channel / gas channeling / free water channel |
| ——（自由下落/U型管） | **free fall / U-tubing / well "goes on a vacuum"** | "in many jobs the well 'goes on a vacuum', 'U tubes' or 'free-falls' while the heavy fluids (cement slurries, spacers) are being pumped down the casing."（Beirute, SPE-13045-MS） | 灌香肠（管内超充/欠顶替）与插旗杆（迟到凝固）的共同水力学根源 |

> 术语验证结论：英文文献中**不存在**与"插旗杆/灌香肠"一一对应的单一术语；两者分别对应一组成因-结果链条（flash set→stuck；under-displacement/fallback→cement left in casing & wet shoe）。写作引用时建议用上述正式术语组合。

---

## 2. 按主题分节文献

### 2.1 Flash set / premature thickening（提前稠化/闪凝机理）

| 文献 | 出处 | 核心发现 | 与插旗杆/灌香肠的关联 |
|---|---|---|---|
| **"Advances in oilwell cement retarders: a bibliometric and systematic review of mechanisms, challenges, emerging trends, and future directions"** (2025) | Springer, DOI 10.1007/s44416-025-00008-6 | 系统综述缓凝剂机理；指出 BHCT **超过 120°C** 时出现 consistency 曲线"鼓包"和温度曲线"波动"的 flash setting 现象；聚羧酸类缓凝剂在室内实验中出现"胶凝层粘附井壁"的闪凝案例；引用文献 155：Umeokafor & Joel SPE-136973-MS；引用文献 156：Zhang H. et al. Constr Build Mater 2021;274:121994 | 插旗杆第一成因（高温闪凝）的近年权威综述，参考文献表可直接用于溯源 |
| **"Occurrence Mechanism of the Abnormal Gelation Phenomenon of High Temperature Cementing Slurry Induced by a Polycarboxylic Retarder"** | PMC10905580（ACS 期刊全文，2024） | 120°C/60 MPa、缓凝剂掺量**≤1.5% bwoc** 时尾浆易发 abnormal gelation (AGP)；厚化曲线异常（rapid thickening, AGP, thickening time inversion）；依据 API RP 10B-2 Section 9 井模拟稠化试验 + SY/T 5504.1-2013，Chandler 8040D 高温高压稠化仪 | 给出"低掺量反而闪凝"的非单调掺量窗口——插旗杆事故复现的实验室路径 |
| **"Oil well casing cement flash setting problem: causes and identification strategy based on cheese model"** | ResearchGate 340538814 | 用 PABM（phenomenon analysis-based method）对闪凝做定义、解释与定位，提出"奶酪模型"识别策略 | 闪凝事故诊断框架 |
| **"Oil and Gas Well Cementing Part 3 – Flash Set and False Set"** (M. Mohammed) | LinkedIn（行业文章） | flash set 定义（混合/泵注中接近初凝的胶凝化）与 false set 区分 | 术语与现场识别 |
| **Umeokafor C., Joel O. (2010) "Modeling of cement thickening time at high temperatures with different retarder concentrations", SPE-136973-MS** | SPE Nigeria Annual International Conference and Exhibition | 不同缓凝剂掺量下高温稠化时间建模 | 高温稠化时间预测（经 Springer 2025 综述引用核实） |
| **Zhang H. et al. (2021) "Inhibitory effects of functionalized polycarboxylate retarder on aberrant thickening phenomena of oil well cement at high temperature"** | Constr Build Mater 274:121994 | 功能化聚羧酸缓凝剂抑制高温异常稠化 | 闪凝防治化学路线（经 Springer 2025 综述引用核实） |
| **Li Z. et al. (2016) "Contamination of cement slurries with oil based mud and [drilling fluid additives]"** | J. Natural Gas Science & Engineering（ScienceDirect S1875510016300038，被引 56） | 油基泥浆与水泥浆混浆不相容导致稠化/强度性能劣化 | 混浆污染→提前凝固或超缓凝双向风险（插旗杆+候凝缺陷共同根源） |
| **Arbad N. et al. (2020) "A Review of Recent Research on Contamination of Oil Well Cement"** | ChemEngineering 4(2):28（MDPI，被引 40） | 综述钻井液（尤其柴油/油基）对油井水泥污染的研究进展 | 污染闪凝/缓凝失效机理谱系 |
| **"Contamination of Oil-Well Cement with Conventional and Microemulsion Spacers"** | SPE Journal 25(06):3002 (2020, OnePetro) | 微乳液隔离液混入水泥浆会造成长期耐久性问题 | 隔离液-水泥相容性 |
| **"A study on the performance effects and mechanism of a [drilling fluid treatment agent contaminating cement slurry]"** | PMC13101578 (2025) | 逐一考察单个泥浆处理剂对水泥浆性能与机理的影响 | 污染源分解实验设计 |
| **"Performance experiment of ultra high temperature cementing slurry system"** | Frontiers in Materials 11:1383286 (2024) | 160–210°C DRH-2L 缓凝剂掺量-稠化时间数据表（如 180°C 掺量 3.0%→TT 92 min；8.5%→70 min，非线性） | 超高温缓凝掺量非单调性数据 |
| **"Cement Slurry Laboratory Testing"** | better-cementing-for-all.org（行业博客，专业性强） | 稠化时间试验设计准则：用**最高循环温度**（非 BHCT）；泵注时间+**安全系数不少于 2 小时**；尾管/打塞/分级固井需在稠化试验中加入**静止段**（识别 POOH 前胶凝倾向）；CT 注水泥自滤失"porous media"堵塞"with catastrophic results" | 稠化时间安全系数设计原则的直接来源（第 2.7 节标准部分亦引） |
| **Khetani (2021) "A Comprehensive Study of Cementing Operation for HPHT [wells]"** | Stanford Geothermal Workshop (SGW/2021/Khetani.pdf) | 套管固井稠化时间"normally 2 to 3 hrs. But due to safety constraints an hour is additionally add[ed]"；稠化试验用 BHCT（低于 BHST） | HPHT 稠化时间裕量现场口径 |

### 2.2 Stuck pipe/casing during cementing（固井相关卡钻/卡套管：统计、处理、预防）

| 文献/来源 | 出处 | 核心发现 | 关联 |
|---|---|---|---|
| **Drilling Manual, "Troubleshooting Pipe Sticking in Oil and Gas Rigs" (Stuck Pipe Course)** | drillingmanual.com/stuck-pipe-course | 水泥卡钻两大类：①鞋跟水泥块崩落卡 BHA；②水泥在钻杆周围凝固——成因清单：打水泥塞时闪凝、管串（钻杆+插入式 stinger 固井）泄漏、人为失误（顶替时间算错/打塞忘加缓凝剂）、下钻嵌入未凝水泥；以及井控/堵漏应急注水泥造成的 cemented sticking | "插旗杆"成因分类的最完整工程表述 |
| **Drillopedia, "Mechanical pipe sticking – mechanism, indicators, prevention, recovery"** | drillopedia.com/mechanical-sticking | 绿水泥卡钻征兆（蹩扭矩/上提超拉）、处理（最大载荷内上提+震击，"Since the cement is still not fully set, this should be able to pull the string out of the green cement"）、预防（good cementing practices） | 早期处置规程 |
| **Drilling Manual LinkedIn 帖（2025）"Causes and Consequences of Stuck Pipe in Cement"** | LinkedIn | 水泥块卡钻机理：BHA 激动压力使未凝水泥 flash set；清扫新凝水泥时高机械钻速同样诱发 | 激动压力-闪凝耦合机制 |
| **Muqeem M.A., Weekse A.E., Al-Hajji A.A. (2012) "Stuck Pipe Best Practices - A Challenging Approach to Reducing Stuck Pipe Costs", SPE-160845-MS** | SPE Saudi Arabia Section Technical Symposium（被引 93） | "Various industry estimates claim that stuck pipe costs may exceed several hundred million US dollars per year"；沙特阿美卡钻最佳实践体系 | 卡钻经济量级（固井相关为其中一类成因） |
| **Salminen K. et al. (2016) "Stuck Pipe Prediction Using Automated Real-Time Modeling and Data Analysis", IADC/SPE-178888-MS** | IADC/SPE Drilling Conference | 卡钻"cost the petroleum industry hundreds of millions of dollars annually"；实时预测 | 卡钻统计与预测 |
| **Pflügl J.K. et al. (2019) "Mechanical Stuck Pipe Events - Development of Digital [catalogue/classification]"** | Montanuniversität Leoben (pure.unileoben.ac.at AC17211200) | 卡钻事件平均持续 **6 天**、占 NPT **12%**（样本统计） | 卡钻 NPT 量级 |
| **"A Comprehensive Review of the Pipe Sticking Mechanism in Oil Well Drilling Operations" (2024)** | ResearchGate 385420942 | "Historically, stuck-pipe events have been shown to cost the industry several hundred million dollars annually and over 25% of non-productive [time]"；含水泥块/绿水泥/压差等成因章节 | 卡钻成因综述 |
| **Suranta B.Y. et al. (2026) "Identification and Mitigation of Stuck Pipe Problem During [cement drill-out operations]"** | SCOG (Lemigas 期刊) | 建立 pack-off 型卡钻（固井后钻水泥塞阶段）综合诊断框架 | 固井后钻塞阶段卡钻专题 |
| **Rogers H. (2006) "Drillable Tailpipe Disconnect: Used Successfully in More Than 120 [wells]", SPE-102534-MS** | SPE ATCE, San Antonio | 尾管（tailpipe）断开工具：打塞 POOH 时尾管被水泥凝固/卡埋的预防手段，>120 口井案例 | 打塞场景"插旗杆"的工程预防 |
| **Cocking (1991) "Task force approach to reducing stuck pipe costs", SPE/IADC Drilling Conference, Amsterdam** | 经 Politecnico di Milano 论文（Matteucci et al.）参考文献核实 | 卡钻成本攻关任务组（早期经典统计文献） | 卡钻统计源头之一 |

**小结**：英文统计文献中卡钻整体占 NPT 12–25%+，年损失数亿美元；**专门统计"固井作业中卡钻占比"的公开论文未检索到**（Muqeem/Salminen 等按差压/井眼清洁/键槽分类，水泥类成因在 Drilling Manual/Drillopedia 教材型资料中有系统定性描述但无占比数字）——这是可引用的表述边界。

### 2.3 U-tube effect / free fall during primary cementing（自由下落/U型管）

| 文献 | 出处 | 核心发现 | 关联 |
|---|---|---|---|
| **Beirute R.M. (1984) "The Phenomenon of Free Fall During Primary Cementing", SPE-13045-MS**（已提取全文摘要） | SPE ATCE 1984, DOI 10.2118/13045-MS | ①自由下落期间**井口返速≠泵速**；②早期加速段在自由下落液柱与井口间形成 "partially empty / discontinuous zone or gap"（管内不连续/脱空段）；③早期加速可能破坏塞流设计，末期减速可跌破紊流顶替最低排量；④建立自由下落数学模型并入固井模拟器，可预测返出排量、井口与井底压力；⑤与多个现场案例对比吻合 | "灌香肠"的核心水力学文献；管内不连续液柱=欠顶替/超顶替判别基础 |
| **Campos W., Lage A.C.V.M., Poggio A. Jr (1993) "Free-fall-effect calculation ensures better cement-operation design", SPE-21107-PA** | SPE Drilling Engineering 8(3), DOI 10.2118/21107-PA（OSTI 5962210 已提取） | 宏观质量+动量守恒把全场方程降为 **1D 模型**，Runge-Kutta 求解初值问题；**专门处理自由下落中两种流体界面的控制**；微机程序支持复杂井身结构 | 与本项目 1D-2D 耦合框架最直接对标的经典（1D 自由下落+界面追踪） |
| **Kelessidis V.C., Rafferty R., Merlo A., Maglione R. (1994) "Simulator models U-tubing to improve primary cementing"** | Oil & Gas Journal 92(10)（OSTI 5180509 已提取） | U-tubing 起因=casing/环空流体密度差；用 **6 口海上井**数据验证：准确预测 U-tubing 起止时刻、环空返速与井口压力；"Proper modeling of U-tubing can help identify potential downhole problems such as formation fracturing leading to lost circulation or inadequate mud removal leading to a poor cement bond" | U-tubing 模型现场验证标杆 |
| **Chiney A., Yerubandi K.B. (2013) "3D Displacement Simulator Realistically Predicts Free Fall during Cementing"** | SPE/IADC Middle East Drilling Technology Conference, Dubai | 3D 顶替模拟器嵌入自由下落预测（经 MDPI Energies 17(5):1226 参考文献核实） | 1D→3D 自由下落耦合的现代实现 |
| **Bogaerts M. et al. (2019) "Novel 3D Fluid Displacement Simulations Improve Cement Job Design and Planning in the Gulf of Mexico"** | SPE ATCE, Calgary | 3D 流体顶替模拟改进固井设计（经 MDPI 17(5):1226 参考文献） | 商业 3D 模拟器现状 |
| **Ryan D.F., Kellingray D.S., Lockyear C.F. "Improved cement placement on North Sea wells using a cement placement simulator"** | SPE（编号待核，见 Bittleston 2002 参考文献） | 北海井用水泥顶替模拟器改善注水泥就位 | 早期商业模拟器应用 |
| **Enayatpour S., van Oort E. (2017) "Advanced modeling of cement displacement complexities"** | SPE/IADC Drilling Conference, The Hague | 顶替复杂性高级建模（经 MDPI 17(5):1226 参考文献） | 顶替建模进展 |
| **Marquairi R., Brisac J. (1966) "Primary Cementing by Reverse Circulation Solves Critical Problem in the North Hassi-Messaoud Field, Algeria", JPT 18(02):146–150** | JPT | 反循环注水泥解决严重漏失问题（反向 U-tube 的工程利用） | 反循环固井源头文献 |
| **Davies J. et al. (2004) "Reverse Circulation of Primary Cementing Jobs—Evaluation and Case History", SPE-87197-MS** | SPE/IADC Drilling Conference | 反循环固井评价+案例 | 反循环 vs 正循环 |
| **Kuru E., Seatter S. (2005) "Reverse Circulation Placement Technique Versus Conventional Placement Technique: A Comparative Study of Cement Job Hydraulics Design", JCPT 44(7):16–19** | JCPT | 反循环/正循环水力学设计对比 | 同上 |
| **Macfarlan K.H. et al. (2017) "A Comparative Hydraulic Analysis of Conventional- and Reverse-Circulation Primary Cementing in Offshore Wells", SPE Drill. Complet. 32(1):59–68** | SPE DC | 海上井正/反循环固井水力学对比 | 同上 |
| **"A Numerical Study of Density-Unstable Reverse [Circulation Cementing]"** | semanticscholar PDF (a3fe/f24eee…) | 反循环末期窄偏心环空充满高粘水泥浆、摩擦压耗大，管内为低粘流体——密度失稳数值研究 | 反循环顶替界面失稳 |

> 备注：任务书中猜测的 "SPE 81635" 未检索到与自由下落/U-tube 的对应关系（搜索仅返回无关结果）；自由下落的正主是 **SPE-13045（Beirute 1984）**。

### 2.4 Wet shoe / dry shoe、cement fallback（湿鞋跟/回流留塞）

| 文献 | 出处 | 核心发现 | 关联 |
|---|---|---|---|
| **Hibbeler J. (2000) "New Float Collar Design to Eliminate Wet Shoe Tracks", SPE-62751-MS** | IADC/SPE Asia Pacific Drilling Technology, Kuala Lumpur, Sept 2000（OnePetro 131867） | 给出 wet shoe track 正式定义（"unset, contaminated or no cement in the casing section between float collar and shoe"）；新浮箍设计+越南成功案例 | "灌香肠"鞋跟形态的正式术语文献 |
| **Drilling Manual, "Float Shoe And Float Collar In Drilling Operations"** | drillingmanual.com | "Failure of the float to operate properly will allow the cement to flow back into the casing. If this happens it can be countered by holding back [pressure]" | 浮阀失效→回流留塞（fallback）机理与对策 |
| **SLB Energy Glossary: "float shoe" / "shoe track"** | glossary.slb.com | 浮鞋单向阀"prevents reverse flow, or U-tubing, of cement slurry from the annulus into the casing"；shoe track "usually left full of cement on the inside" | 术语定义 |
| **Pegasus Vertex, "Common Well Cementing Problems and Solutions"（白皮书）** | linqx.io/white-paper/…pdf | 系统列举常见固井问题：gas flow after placement、zonal communication、poor displacement efficiency、cement failure、fluid influx during pumping、lost circulation、**poor pumpability**、**wet shoe track**、lifted casing 等，并给 remediation | 问题-对策对照表（灌香肠类问题的工程清单） |
| **ScienceDirect Topics: "Casing Shoe"（教材汇编）** | sciencedirect.com/topics/engineering/casing-shoe | TOC 低于井口形成"closed annulus"；"A fluid can even be trapped in the annulus due to **cement drop-off from wellhead level at the end of the cementing** even though a full cement column to wellhead level is originally planned"——返高回落、圈闭流体、APB | 返高不达/回落（灌香肠环空形态）与 APB 风险 |
| **Purvis D.L. et al. (2009) "Wellbore Re-entries and Repairs: Practical Guidelines for Cementing New Casing Inside Existing Casing", SPE-124381-MS**（Drilling Contractor 转载） | SPE ATCE 2009 | 修复井固井水泥"not quite making it to surface, again leaving some exposed pipe surface"；TOC 实测低于泄漏点（温度测井 1450 ft）导致腐蚀复发 | 欠顶替/返高不足的长期后果案例集 |
| **Innovex / Forum Energy Technologies（Davis-Lynch PDQ、M.O.A.S. 等）产品技术页** | 厂商资料 | 现代 wet shoe system/一体化鞋跟："makes shoe track cement unnecessary as a formation barrier"——取消鞋跟水泥塞的工程趋势 | 防湿鞋跟装备路线 |

### 2.5 Cement channeling / segmented-discontinuous cement column / gas migration（窜槽、返高不连续、气窜）

| 文献 | 出处 | 核心发现 | 关联 |
|---|---|---|---|
| **"An Overview of Annular Displacement Efficiency in Cementing Jobs Using an Efficient Numerical Model", AADE-18-FTCE-094** | AADE Fluids Technical Conference | 环空体积顶替效率定义 E(t)=1−V_mud(t)/V_annulus，按 TOC 以下积分；算例：9-5/8" 尾管居中度差→窄边显著 mud channel；13-5/8" 表层 1830 ft TOC 窄边窜槽；自研 FVM 模型与国际主流模拟器对比 | 窜槽=环空"未充满水泥"的量化框架（灌香肠的环空镜像） |
| **"Cement Placement: An Overview of Fluid Displacement Techniques and Modelling", Energies 14(3):573 (2021)** | MDPI | 顶替技术+建模大综述；gas migration 是"most dangerous and the most frequent kind of annular fluid migration"；引 Martin/Latil/Vetter 1977、Haut & Crook 1979、Jakobsen 1991 等 | 窜槽机理谱系入口文献 |
| **"A Brief Review of Gas Migration in Oilwell Cement Slurries", Energies 14(9):2369 (2021)** | MDPI | 过渡时间多种定义（right-angle set、WOC、SGS→500 lbf/100ft²）；SGS 失去静液柱传递能力→气窜；引 Sabins & Tinsley（522 lbf/100ft² 终点） | 停泵后环空压力亏空与窜槽的时间窗口 |
| **Tinsley J.M., Miller E.C., Sabins F.L., Sutton D.L. (1980) "Study of factors causing annular gas flow following primary cementing", JPT 32:1427–1437** | JPT | 环空气窜成因因素研究（经典） | 停泵后窜槽经典源头 |
| **Sabins F.L. et al. (1979/1980) [transition state concept], SPE 9285；Sabins & Sutton（SPEPE）"The Relationship of Thickening Time, Gel Strength, and Compressive Strength of Oil Well Cements"** | SPE（经 Wellcem/M GPI 综述核实） | "transition state"：既非流体又非固体的中间期，浆体失去静液压力传递能力；起于首个可测胶凝强度（约 21 lb/100ft²），止于气体不能渗透（~500 lb/100ft²） | 停泵后到凝固前的时间窗=窜槽/气窜风险窗 |
| **OFITE "Mitigating Gas Migration by Measuring Static Gel Strength"** | ofite.com | "Gas intrusion begins when cement slurries develop 100 lbf/100 ft² (48 Pa). Gas intrusion ends at 500 lbf/100 ft² (240 Pa)" | SGS 过渡时间量化 |
| **Drilling Contractor "New cementing method uses pipe movement to maximize displacement"** | drillingcontractor.org | standoff 与顶替效率关系；"if gas migration is likely, then a full column of cement in the entire length of the annulus is required. Although well intended, this practice has not always proven effective. The additional cement length increases the hydrostatic pressure initially. Then, as it moves to transition, its sudden loss creates more pressure loss and exacerbates the gas migration tendency."——全环空水泥柱策略反而加剧气窜 | "加大水泥返高"对策的证伪性证据 |
| **Wellcem "Free water in cement slurries for oil and gas wells: Big trouble?"** | wellcem.com | 游离水在长水泥柱顶部/大斜度井窝集："large amounts of free water in a highly deviated well might lead to a communication channel on the high side"；长稠化时间加剧；后期 SCA/套管腐蚀/挤毁 | 游离水通道=环空分段（水柱-水泥柱交替）的直接成因 |
| **"Cementing Techniques for Horizontal Wells"** | Scribd（教材型文档） | "Free water is one cause of unwanted communication channeling during deviated well completions" | 同上 |
| **Bittleston S.H., Ferguson J., Frigaard I.A. (2002) "Mud removal and cement placement during primary cementing of an oil well – Laminar non-Newtonian displacements in an eccentric annular Hele-Shaw cell", J. Eng. Math. 43:229–253** | Springer（已核实卷页码与作者单位：SCR + UBC） | 窄偏心环空 Hele-Shaw 2D 顶替模型：流函数泊松型方程+浓度对流方程组+浓度依赖物性闭合；偏心→窄边泥浆滞留 | 环空窜槽 2D 建模基石（与本项目 D2DGA 同族框架） |
| **Tehrani A., Ferguson J., Bittleston S.H. (1992) "Laminar displacement in annuli: A combined experimental and theoretical study", SPE-24569**；**Tehrani, Bittleston, Long (1993) Flow instabilities during annular displacement of one non-Newtonian fluid by another, Exp. Fluids 14:246–256** | SPE / Exp Fluids | 环空层流顶替实验+理论；非牛顿流体顶替失稳 | 窜槽失稳实验经典 |
| **Pelipenko S., Frigaard I.A. (2004) "On steady state displacements in primary cementing of an oil well", J. Eng. Math. 46:1–26**（及后续两篇） | Springer | 稳态顶替的有效顶替判据（有效顶替数） | 偏心环空顶替设计准则 |
| **Renteria A., Frigaard I.A. (2020) "Primary cementing of horizontal wells. Displacement flows in eccentric horizontal annuli. Part 1. Experiments", J. Fluid Mech. 905:A7** | JFM | 水平井偏心环空顶替实验（系列第一部分） | 窜槽流变-几何耦合最新实验 |
| **Maleki A., Frigaard I.A. [3D primary cementing model]**（经 MDPI 17(5):1226 核实存在） | （期刊待核） | "Maleki and Frigaard also built a 3D model based on the narrow gap assumption and used it to study primary cementing" | 3D 窄隙模型 |

### 2.6 Cement plug placement（打水泥塞：平衡塞、粘卡、污染）

| 文献 | 出处 | 核心发现 | 关联 |
|---|---|---|---|
| **"Current Industry Practices: Placing Cement Plugs in a Wellbore Using a Stinger or Tail-Pipe", SPE-168005-MS (2014)** | SLB 技术论文页（已提取要点） | 平衡塞法配小直径 stinger 时，POOH 过程管内外驱替体积不等→系统失平衡："using the balanced plug method will result in the contamination and potential failure of the plug"；JPT 评述：塞顶部污染 **40–50% 水泥体积**，且"much higher likelihood of fluid U-tubing out of the drillpipe, resulting in spills on the rig floor" | 打塞场景 U-tube→塞内污染+管内留水泥（灌香肠的塞版本） |
| **JPT (2016?) "Stinger or Tailpipe Placement of Cement Plugs"**（SPE-168005 评述） | jpt.spe.org | POOH 后界面位置逐步计算（7307.6 ft 处 spacer/cement 界面下移等定量示例）；stinger 段留水泥 176.51 ft | 打塞管内留塞的定量案例 |
| **Rogers H. (2006) SPE-102534-MS** | SPE ATCE | 尾管可钻断开工具（>120 井案例），防 POOH 时尾管被凝/卡 | 插旗杆预防（打塞） |
| **Bois A.-P. et al. (2019) "Cement Plug Hydraulic Integrity—The Ultimate Objective of Cement Plug Integrity", SPE 191335**（经 JPT Cementing and Zonal Isolation-2019 栏目核实） | SPE ATCE 2019 | 水泥塞水力完整性目标框架 | 塞质量评价 |
| **Aas B. et al. (2012) "Cement placement with tubing left in hole during plug and abandonment operations"** | IADC/SPE Drilling Conference, San Diego（经 LAPSE-2023.29256 参考文献 42 核实） | P&A 中留管注水泥就位 | 留管塞工艺 |
| **Drillopedia "Well P&A Cementing Plug Placement Design"** | drillopedia.com/pa-cement-plug | 塞设计标准清单：API RP 65-2 (2010)、API RP 65-3、API RP 10B 系列、NORSOK D-010 (2021)、OGUK (2023)、Schlumberger Cementing Handbook 2022、Halliburton Cementing Design Guide 2021 | 塞设计标准地图 |
| **SPE/IADC-221446-MS (2024) "Managed Pressure Cementing in [drilling operations]"** | IADC（PDF 可达） | MPD 环境下打 5 个水泥塞的泵注/回压程序（Figure 8/Table 5） | 控压固井+打塞结合（与本项目控压背景直接相关） |

### 2.7 相关标准（API/ISO/NORSOK）

| 标准 | 内容要点（来自搜索结果） |
|---|---|
| **API RP 10B-2 / ISO 10426-2** | 油井水泥浆实验室试验推荐做法；ISO 10426-2:2003 定义 3.1.49 thickening time = "time required for a cement slurry to develop a selected Bearden consistency value...provide an indication of the length of time a cement slurry can remain pumpable under the test conditions"；Section 9 = Well-simulation Thickening Time Tests（高温高压稠化仪 8040D 执行，见 PMC10905580）；Section 15 特殊工况试验 |
| **API 10A** | 油井水泥规范：Class A–H 八类+抗硫酸盐级别 O/MSR/HSR（Heidelberg Materials 页） |
| **API RP 10B-3/4/5/6 (ISO 10426-3/4/5/6)** | 深水配方试验 / 泡沫水泥 / 收缩膨胀 / 静胶凝强度（Wellcem 标准页；ISO/FDIS 10426-4:2019 文件） |
| **API RP 65 (2002) / RP 65-2 (2010) "Isolating Potential Flow Zones During Well Construction"** | 防环空流动分区隔离工程过程；JPT 58(01):53 (2006) "New API Practices for Isolating Potential Flow Zones..."（OnePetro）；美国 OCS/BSEE 法规引用（Federal Register E6-7792）；俄亥俄州规则要求按 65-2 设计防气窜水泥浆且钻塞前抗压 ≥500 psi |
| **API RP 65-3** | Wellbore Plugging and Abandonment（P&A 塞设计验证，Drillopedia） |
| **NORSOK D-010 (Rev.5, 2021)** | 永久屏障要求（P&A） |
| **稠化时间安全系数设计原则** | better-cementing-for-all：泵注作业总时间（混配+泵注+顶替+静止事件如投塞/起塞管/循环多余水泥）+ **≥2 小时**安全系数；用最高循环温度做试验；尾管/打塞/分级固井试验中插入静止段识别胶凝倾向。Khetani (2021)：套管固井 TT 常按 2–3 h+1 h 裕量 |

### 2.8 数值模拟方向（固井水力学 1D/2D、界面追踪、停泵后沉降/凝结）

| 文献 | 出处 | 核心发现 | 关联 |
|---|---|---|---|
| **Beirute 1984 / Campos 1993 / Kelessidis 1994** | 见 2.3 | 1D 自由下落+界面控制模型谱系（1984 提出→1993 Runge-Kutta 1D→1994 六井现场验证） | 与本项目 1D-2D 耦合直接对标 |
| **Bois A.-P. et al. (2023) "Cement Placement Modeling—A Review", SPE-214331-PA, SPE Drill. Complet. 38(02):342–…** | OnePetro（OnePetro DC 38(02)） | "presents an extensive analysis of all the available cement placement computerized models, highlighting their advantages and disadvantages" | 水泥就位建模最新权威综述（选型对照必引） |
| **"Cement Placement: An Overview of Fluid Displacement Techniques and Modelling"（两个版本）** | Energies 14(3):573 (2021) + LAPSE 2023.29256 | Hele-Shaw 降维谱系（Bittleston 2002→Pelipenko 2004）、界面方程、量纲化约；LAPSE 版含 50+ 参考文献表 | 2D 顶替建模综述（含本项目同类模型定位） |
| **"Advanced Mud Displacement Modeling for Slim Hole Cementing Operations", Energies 17(5):1226 (2024)** | MDPI | 小井眼顶替建模，参考文献表含 McLean 1967、Beirute & Flumerfelt 1977、Chiney 2013、Bogaerts 2019、Enayatpour 2017、Li & Novotny 2006（格子-Boltzmann）、Zhang 2023（SPE J 28:509–521 螺旋顶替）、Maleki & Frigaard 3D | 参考文献表可整体复用 |
| **"Numerical Simulation on the Safety and Quality of Cementing by Using Pad Fluid in Horizontal Wells", Energies 16(9):3650 (2023)** | MDPI | N-S 方程+VOF 界面捕捉模拟前置液顶替；顶替效率与 ECD 双目标（安全密度窗口） | ECD-效率耦合模拟 |
| **"Numerical Analysis of Cement Placement into Drilling Fluid in [eccentric annulus]"** | PMC12250809 (2025) | VOF 多相（泥浆/隔离液/水泥/地层流体）界面模拟，非牛顿水泥本构讨论 | CFD 顶替 |
| **"Numerical modeling of motion of displacement interface in eccentric annulus during primary cementing"** | ResearchGate 338826998 | 偏心环空顶替界面运动数值模拟 | 界面追踪 |
| **"Theoretical and Experimental Study on Cementing [Displacement Interface]"** | LAPSE 2023.6275 | 偏心环空两相界面连续条件数学模型；无量纲界面长度/环空长度比 7.4 对比 | 界面模型-实验对照 |
| **Foroushan H.K., Ozbayoglu E., Gomes P.J. (2020) "How Realistic is the Calculated Cementing Displacement Efficiency?", SPE-199553** | IADC/SPE Drilling Conference, Galveston | 顶替效率计算的"现实性"拷问（口径/假设敏感性） | 顶替效率指标审计（与本项目 η 口径议题呼应） |
| **Foroushan H.K. (2018) "Displacement of Fluids in Annuli" (PhD thesis, U. of Tulsa)** | TUDRP | 环空顶替流体力学 | 系统性学位论文 |
| **AADE-18-FTCE-094** | AADE | FVM 固井模拟器 vs 商业模拟器（9-5/8" 尾管与 13-5/8" 表层算例 mud channel 对比） | 模拟器对标 |
| **"A Numerical Study of Density-Unstable Reverse [Circulation]"** | semanticscholar PDF | 反循环固井密度失稳数值研究（含 Marquairi 1966、Davies 2004 等参考文献链） | 反循环模拟 |
| **停泵后水泥沉降/凝结模拟** | — | **未检索到水泥浆停泵后沉降-胶凝全过程 CFD 模拟的专门英文文献**；相关碎片：颗粒沉降 CFD-DEM 方法在其他领域成熟（arXiv 1711.01524：非粘性颗粒沉降 CFD-DEM；MDPI JMSE 11(9):1685：多面体颗粒自由沉降）；水泥侧仅有浆体稳定性评价指标（free water/sedimentation，如 PMC11547183 Zhang 2024 优化浆体 sedimentation<…、API free water 试验）与 SGS 过渡时间框架（2.5 节）。**这是一个明确的研究空白** | 论文创新点候选 |

---

## 3. 机理总结：提前凝固/闪凝成因分类（英文文献归纳）

1. **温度因素（HTHP 闪凝主因）**
   - BHCT > 120°C 时缓凝剂失效窗口收缩，稠化曲线出现"鼓包"（bulges）与温度波动（Springer 2025 综述）；
   - 超高温（160–210°C）缓凝剂掺量-稠化时间呈**非单调**关系（Frontiers 2024：180°C 下掺量 3.0%→TT 92 min，8.5%→70 min）；
   - 低掺量聚羧酸缓凝剂（≤1.5% bwoc）在 120°C/60 MPa 反而诱发 AGP（PMC10905580）。
2. **化学污染/相容性**
   - 油基泥浆及处理剂混入（Li 2016；Arbad 2020；PMC13101578）；
   - 隔离液（含微乳液）混入（SPE J 2020）。
3. **水力-作业因素（插旗杆直接诱因，Drilling Manual 清单）**
   - 打塞时水泥闪凝；管串泄漏（stinger 固井漏失使水泥滞留管内）；
   - 顶替时间误算、忘加缓凝剂（人为失误）；
   - 下钻嵌入未凝水泥；BHA 激动压力+高 ROP 钻扫未凝水泥诱发 flash set（Drilling Manual LinkedIn 帖）；
   - 停泵后 U-tube 持续流动使局部水泥在低速率/静止状态滞留管内至稠化（Campos 1993 / Kelessidis 1994 框架推论——注意：此条为文献机理组合推论，非单一文献原文）。
4. **失水/胶凝连锁**（气窜窗口侧）：
   - 失水形成滤饼桥堵（bridge-off）→静液压力传递损失→过渡期（Sabins 转变态）；
   - 游离水分离窝集（Wellcem）→高边水通道。

## 4. 防治措施汇总（英文文献口径）

| 方向 | 措施 | 来源 |
|---|---|---|
| 配方 | 高温缓凝体系（功能化聚羧酸、DRH-2L 等）、RAS（right-angle set）设计缩短过渡时间、零游离水/沉降稳定（sedimentation）控制、降失水 ~50 mL/30min | Springer 2025；Frontiers 2024；Wellcem；PMC11547183 |
| 实验设计 | 稠化试验用最高循环温度；总作业时间+≥2 h 安全裕量；尾管/打塞试验插入静止段复现 POOH 胶凝 | better-cementing-for-all；API RP 10B-2 §9/§15 |
| 水力学设计 | 自由下落/U-tube 模拟器预测返速-ECD 窗口（避免压破地层与欠顶替）；正/反循环比选 | Beirute 1984；Campos 1993；Kelessidis 1994；Chiney 2013；Kuru & Seatter 2005 |
| 装备 | 双阀浮箍/浮鞋防回流；新浮箍消除湿鞋跟（Hibbeler SPE-62751）；可钻尾管断开工具（SPE-102534）；一体化 wet shoe system | Drilling Manual；Hibbeler 2000；Rogers 2006；Innovex/FET |
| 作业程序 | POOH 前塞顶充分循环；stinger 替代平衡塞或专用计算防污染（40–50% 塞顶污染）；固井后控制钻扫未凝水泥的激动压力与 ROP | SPE-168005/JPT；Drillopedia；Drilling Manual |
| 环空防窜 | 合理返高（非"全环空水泥柱"——Drilling Contractor 证伪）；管柱活动改善窄边顶替；SGS 过渡时间监测 | Drilling Contractor；AADE-18-FTCE-094；OFITE |

## 5. 研究空白（可写入论文的 gap）

1. **固井作业中卡钻（插旗杆）的定量统计缺失**：卡钻总体统计充分（NPT 12–25%、数亿美元/年），但"水泥类成因占卡钻比例"无公开统计——英文文献均为定性教材级描述。
2. **停泵后管内/环空水泥"沉降-胶凝-失返"全过程模拟空白**：1D U-tube 模拟止于泵注结束（Campos/Kelessidis），2D/3D 顶替模拟止于顶替结束（Bittleston/Renteria/Bois 综述均指出），停泵后自由下落残余流动+沉降+SGS 发展的耦合模拟未检索到——与本项目"停泵后顶替"卖点正交衔接（cf. SWPU 空档记忆：停泵后顶替为 SWPU 调研确认的空白方向）。
3. **自由下落 1D 模型与环空 2D 顶替模型的在线耦合**在公开文献中只有 3D 全井模拟器（Chiney 2013、Bogaerts 2019、Enayatpour 2017）与独立 1D/2D 谱系，缺少"1D 管内-环空水力学 + 2D 环空截面顶替"轻量化耦合框架的公开实现——D2DGA 类框架可对位引用 Bois 2023 综述作坐标。
4. **湿鞋跟/回流留塞的概率风险评估**未见系统研究（仅装备改进与个案）。
5. **高温闪凝的井筒实际温-压耦合稠化预测**（井场温度场 vs API §9 试验条件差异）仍有争议空间（bcf4all 指出 BHCT vs 最高循环温度口径）。

## 6. 完整参考文献列表（按主题，均来自检索结果）

**闪凝/提前稠化**
1. Advances in oilwell cement retarders: a bibliometric and systematic review… — Springer, 2025. https://link.springer.com/article/10.1007/s44416-025-00008-6
2. Occurrence Mechanism of the Abnormal Gelation Phenomenon of High Temperature Cementing Slurry Induced by a Polycarboxylic Retarder — PMC10905580 (2024). https://pmc.ncbi.nlm.nih.gov/articles/PMC10905580
3. Umeokafor C., Joel O. Modeling of cement thickening time at high temperatures with different retarder concentrations. SPE-136973-MS, SPE NAICE 2010.
4. Zhang H., et al. Inhibitory effects of functionalized polycarboxylate retarder on aberrant thickening phenomena of oil well cement at high temperature. Constr Build Mater 2021;274:121994.
5. Li Z., et al. Contamination of cement slurries with oil based mud and [drilling fluids]. J Nat Gas Sci Eng 2016. https://www.sciencedirect.com/science/article/abs/pii/S1875510016300038
6. Arbad N., et al. A Review of Recent Research on Contamination of Oil Well Cement. ChemEngineering 2020;4(2):28. https://www.mdpi.com/2305-7084/4/2/28
7. Contamination of Oil-Well Cement with Conventional and Microemulsion Spacers. SPE Journal 2020;25(06):3002. https://onepetro.org/SJ/article/25/06/3002/453926
8. A study on the performance effects and mechanism of a [cementing slurry contaminant] — PMC13101578 (2025). https://pmc.ncbi.nlm.nih.gov/articles/PMC13101578
9. Performance experiment of ultra high temperature cementing slurry system. Frontiers in Materials 2024;11:1383286. https://www.frontiersin.org/journals/materials/articles/10.3389/fmats.2024.1383286/full
10. Oil well casing cement flash setting problem: causes and identification strategy based on cheese model. ResearchGate 340538814. https://www.researchgate.net/publication/340538814
11. Khetani A. A Comprehensive Study of Cementing Operation for HPHT [Geothermal Wells]. Stanford Geothermal Workshop 2021. https://pangea.stanford.edu/ERE/db/GeoConf/papers/SGW/2021/Khetani.pdf
12. Cement Slurry Laboratory Testing. Better Cementing for All. https://better-cementing-for-all.org/cement-slurry-laboratory-testing

**卡钻/固井卡套管**
13. Troubleshooting Pipe Sticking in Oil and Gas Rigs. Drilling Manual. https://www.drillingmanual.com/stuck-pipe-course
14. Mechanical pipe sticking – mechanism, indicators, prevention, recovery. Drillopedia. https://www.drillopedia.com/mechanical-sticking
15. Muqeem M.A., Weekse A.E., Al-Hajji A.A. Stuck Pipe Best Practices – A Challenging Approach to Reducing Stuck Pipe Costs. SPE-160845-MS, 2012. https://onepetro.org/SPESATS/proceedings/12SATS/12SATS/SPE-160845-MS/159022
16. Salminen K., Cheatham C., Smith M., Valiulin K. Stuck Pipe Prediction Using Automated Real-Time Modeling and Data Analysis. IADC/SPE-178888-MS, 2016.
17. Pflügl J.K., et al. Mechanical Stuck Pipe Events – Development of Digital [Catalogue]. Montanuniversität Leoben, 2019. https://pure.unileoben.ac.at/files/27792991/AC17211200.pdf
18. A Comprehensive Review of the Pipe Sticking Mechanism in Oil Well Drilling Operations (2024). https://www.researchgate.net/publication/385420942
19. Suranta B.Y., et al. Identification and Mitigation of Stuck Pipe Problem During [Cement Drill-Out]. SCOG 2026. https://journal.lemigas.esdm.go.id/index.php/SCOG/article/download/1905/1639/7278
20. Rogers H. Drillable Tailpipe Disconnect: Used Successfully in More Than 120 [Wells]. SPE-102534-MS, SPE ATCE 2006. https://onepetro.org/SPEATCE/proceedings-abstract/06ATCE/All-06ATCE/139919
21. Cocking G.A. Task force approach to reducing stuck pipe costs. SPE/IADC Drilling Conf., Amsterdam, 1991.（经 Matteucci et al. 参考文献表核实）

**U-tube / 自由下落**
22. Beirute R.M. The Phenomenon of Free Fall During Primary Cementing. SPE-13045-MS, SPE ATCE 1984. DOI 10.2118/13045-MS. https://colab.ws/articles/10.2118%2F13045-MS
23. Campos W., Lage A.C.V.M., Poggio A. Jr. Free-fall-effect calculation ensures better cement-operation design. SPE Drilling Engineering 1993;8(3). DOI 10.2118/21107-PA. https://www.osti.gov/biblio/5962210
24. Kelessidis V.C., Rafferty R., Merlo A., Maglione R. Simulator models U-tubing to improve primary cementing. Oil & Gas J 1994;92(10). https://www.osti.gov/biblio/5180509
25. Chiney A., Yerubandi K.B. 3D Displacement Simulator Realistically Predicts Free Fall during Cementing. SPE/IADC MEDTE, Dubai, 2013.（经 Energies 17(5):1226 文献表核实）
26. Bogaerts M., et al. Novel 3D Fluid Displacement Simulations Improve Cement Job Design and Planning in the Gulf of Mexico. SPE ATCE, Calgary, 2019.
27. Enayatpour S., van Oort E. Advanced modeling of cement displacement complexities. SPE/IADC Drilling Conf., The Hague, 2017.
28. Marquairi R., Brisac J. Primary Cementing by Reverse Circulation Solves Critical Problem in the North Hassi-Messaoud Field, Algeria. JPT 1966;18(02):146–150.
29. Davies J., et al. Reverse Circulation of Primary Cementing Jobs—Evaluation and Case History. SPE-87197-MS, 2004.
30. Kuru E., Seatter S. Reverse Circulation Placement Technique Versus Conventional Placement Technique. JCPT 2005;44(7):16–19.
31. Macfarlan K.H., et al. A Comparative Hydraulic Analysis of Conventional- and Reverse-Circulation Primary Cementing in Offshore Wells. SPE Drill Complet 2017;32(1):59–68.
32. A Numerical Study of Density-Unstable Reverse [Circulation Cementing]. https://pdfs.semanticscholar.org/a3fe/f24eee0d9d2681c535df2bcf4a7e80e01149.pdf

**湿鞋跟/回流/欠顶替**
33. Hibbeler J. New Float Collar Design to Eliminate Wet Shoe Tracks. SPE-62751-MS, IADC/SPE Asia Pacific Drilling Tech, Kuala Lumpur, 2000. https://onepetro.org/SPEAPDT/proceedings-abstract/00APDT/00APDT/131867
34. Float Shoe And Float Collar In Drilling Operations. Drilling Manual. https://www.drillingmanual.com/float-shoe-collar-equipment-cementing-casing
35. float shoe / shoe track. SLB Energy Glossary. https://glossary.slb.com/terms/f/float_shoe ; https://glossary.slb.com/terms/s/shoe_track.aspx
36. Common Well Cementing Problems and Solutions. Pegasus Vertex white paper. https://linqx.io/white-paper/Common-Well-Cementing-Problems-and-Solutions.pdf
37. Casing Shoe — overview. ScienceDirect Topics. https://www.sciencedirect.com/topics/engineering/casing-shoe
38. Purvis D.L., et al. Wellbore Re-entries and Repairs: Practical Guidelines for Cementing New Casing Inside Existing Casing. SPE-124381-MS, SPE ATCE 2009（Drilling Contractor 转载）

**窜槽/气窜/顶替**
39. An Overview of Annular Displacement Efficiency in Cementing Jobs Using an Efficient Numerical Model. AADE-18-FTCE-094. https://www.aade.org/download_file/1335/394
40. Cement Placement: An Overview of Fluid Displacement Techniques and Modelling. Energies 2021;14(3):573. https://www.mdpi.com/1996-1073/14/3/573 ；LAPSE 版：https://psecommunity.org/.../LAPSE-2023.29256-1v1.pdf
41. A Brief Review of Gas Migration in Oilwell Cement Slurries. Energies 2021;14(9):2369. https://www.mdpi.com/1996-1073/14/9/2369
42. Tinsley J.M., Miller E.C., Sabins F.L., Sutton D.L. Study of factors causing annular gas flow following primary cementing. JPT 1980;32:1427–1437.
43. Sabins F.L., et al. [Transition state of cement slurries]. SPE 9285, 1979/1980；Sabins F.L., Sutton D.L. The Relationship of Thickening Time, Gel Strength, and Compressive Strength of Oil Well Cements. SPEPE.
44. New cementing method uses pipe movement to maximize displacement. Drilling Contractor. https://drillingcontractor.org/new-cementing-method-uses-pipe-movement-to-maximize-displacement-21248
45. Mitigating Gas Migration by Measuring Static Gel Strength of Well Cement. OFITE. https://www.ofite.com/news/mitigating-gas-migration-by-measuring-static-gel-strength-of-well-cement
46. Free water in cement slurries for oil and gas wells: Big trouble? Wellcem. https://wellcem.com/news-and-articles/free-water-in-cement-slurries-for-oil-and-gas-wells-big-trouble
47. Bittleston S.H., Ferguson J., Frigaard I.A. Mud removal and cement placement during primary cementing of an oil well – Laminar non-Newtonian displacements in an eccentric annular Hele-Shaw cell. J Eng Math 2002;43:229–253. DOI 10.1023/A:1020370417367. https://link.springer.com/article/10.1023/A:1020370417367
48. Tehrani A., Ferguson J., Bittleston S.H. Laminar displacement in annuli: a combined experimental and theoretical study. SPE-24569, 1992.
49. Tehrani A., Bittleston S.H., Long P.J.G. Flow instabilities during annular displacement of one non-Newtonian fluid by another. Exp Fluids 1993;14:246–256.
50. Pelipenko S., Frigaard I.A. On steady state displacements in primary cementing of an oil well. J Eng Math 2004;46:1–26.
51. Renteria A., Frigaard I.A. Primary cementing of horizontal wells. Displacement flows in eccentric horizontal annuli. Part 1. Experiments. J Fluid Mech 2020;905:A7.
52. Research on Key Technologies to Improve Cementing Displacement Efficiency. PMC9607674. https://pmc.ncbi.nlm.nih.gov/articles/PMC9607674

**打水泥塞**
53. Current Industry Practices: Placing Cement Plugs in a Wellbore Using a Stinger or Tail-Pipe. SPE-168005-MS, 2014. https://www.slb.com/resource-library/technical-paper/ce/spe-168005
54. Stinger or Tailpipe Placement of Cement Plugs. JPT. https://jpt.spe.org/stinger-or-tailpipe-placement-cement-plugs
55. Balanced Cement Plug Guide In Oil & Gas Wells. Drilling Manual. https://www.drillingmanual.com/cement-plug-procedure-calculation-excel
56. Bois A.-P., et al. Cement Plug Hydraulic Integrity—The Ultimate Objective of Cement Plug Integrity. SPE 191335, 2019（经 JPT Cementing and Zonal Isolation-2019 核实）
57. Aas B., et al. Cement placement with tubing left in hole during plug and abandonment operations. IADC/SPE Drilling Conf., San Diego, 2012.
58. SPE/IADC-221446-MS. Managed Pressure Cementing in [Middle East operations]. 2024. https://iadc.org/wp-content/uploads/2024/09/SPE-221446-MS.pdf
59. Well P&A Cementing Plug Placement Design. Drillopedia. https://www.drillopedia.com/pa-cement-plug

**标准**
60. API RP 10B-2 / ISO 10426-2:2003 (sample PDF). https://cdn.standards.iteh.ai/samples/37866/d83023913ff241d18eb58fa1cae76d15/ISO-10426-2-2003.pdf
61. Standards you need to know when designing a cement slurry. Wellcem. https://wellcem.com/news-and-articles/standards-you-need-to-know-when-designing-a-cement-slurry-for-oil-well-applications
62. New API Practices for Isolating Potential Flow Zones [During Well Construction]. JPT 2006;58(01):53. https://onepetro.org/JPT/article/58/01/53/196207
63. API RP 65 (shallow water flow) incorporation, BSEE/US MMS. Federal Register E6-7792 (2006). https://www.federalregister.gov/documents/2006/05/22/E6-7792/...
64. ISO/FDIS 10426-4:2019 (foamed cement). https://www.tksneftegaz.ru/fileadmin/f/subcommittees/sc6/consideration/12-10_ISO_FDIS_10426-4.pdf
65. OFI Testing Equipment: Testing Equipment for API RP 10B-2. https://www.ofite.com/knowledgebase/api-rp-10b-2

**数值模拟**
66. Bois A.-P., et al. Cement Placement Modeling—A Review. SPE-214331-PA, SPE Drill Complet 2023;38(02):342. https://onepetro.org/DC/article-abstract/doi/10.2118/214331-PA/518528
67. Advanced Mud Displacement Modeling for Slim Hole Cementing Operations. Energies 2024;17(5):1226. https://www.mdpi.com/1996-1073/17/5/1226
68. Numerical Simulation on the Safety and Quality of Cementing by Using Pad Fluid in Horizontal Wells. Energies 2023;16(9):3650. https://www.mdpi.com/1996-1073/16/9/3650
69. Numerical Analysis of Cement Placement into Drilling Fluid in [annuli] (VOF). PMC12250809 (2025). https://pmc.ncbi.nlm.nih.gov/articles/PMC12250809
70. Numerical modeling of motion of displacement interface in eccentric annulus during primary cementing. ResearchGate 338826998.
71. Theoretical and Experimental Study on Cementing [Displacement Interface]. LAPSE-2023.6275. https://psecommunity.org/wp-content/plugins/wpor/includes/file/2302/LAPSE-2023.6275-1v1.pdf
72. Foroushan H.K., Ozbayoglu E., Gomes P.J. How Realistic is the Calculated Cementing Displacement Efficiency? SPE-199553, IADC/SPE Drilling Conf., Galveston, 2020（经 TUDRP 出版列表核实）
73. Foroushan H.K. Displacement of Fluids in Annuli. PhD thesis, University of Tulsa, 2018.
74. Li X., Novotny R.J. Study on cement displacement by lattice-Boltzmann method. SPE ATCE, San Antonio, 2006.
75. Zhang Z., et al. Characteristic of Spiral Displacement Process in Primary Cementing of Vertical Well Washout. SPE J 2023;28:509–521.
76. Study of sedimentation of non-cohesive particles via CFD-DEM. arXiv:1711.01524. https://arxiv.org/pdf/1711.01524
77. Nelson E.B., Guillot D. (Eds.) Well Cementing, 2nd ed. Schlumberger, 2006.（各综述一致引用的工具书）

---

### 附：检索执行记录（可复核）
- 共 26 轮 Tavily 搜索（每主题 ≥2 轮不同关键词）+ 3 次全文提取（colab.ws SPE-13045 页、OSTI 5180509、OSTI 5962210）。
- 检索未果项（如实记录）：① "SPE 81635" 与自由下落的关联未证实；② 固井作业卡钻占比的公开统计；③ 停泵后水泥沉降-凝结 CFD 模拟专文；④ Maleki & Frigaard 3D 模型的具体期刊出处（仅确认存在）。
