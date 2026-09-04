# 既往报告挖掘摘要：居中度敏感性 / 1D-2D 耦合模型可靠性 / 数值问题（供主编综合）

> 挖掘日期：2026-09-03 ｜ 来源目录：`D:\obsidian\obsidian-storage\1.科研\新疆油田控压固井项目\固井顶替效率改进\` ｜ 只读挖掘，未改动任何 obsidian 文件。
> **全局事实基线（所有 2026-09-02 状态声明共同指向）**：cement model commit `6d58fa0` 落地 B1（环空速度截面归一化多除因子 2 守恒 bug，修复前塞流守恒率 0.493→修复后 0.980）+ B2（壁面残泥判据改可逆物理屈服门 τw≤1.15τy）+ B3（管容双链统一）；后续 0708 管容核实（`3d88385`）、F4 erf 收尾、胶塞语义修复（`21d4490`）后，**8 井 η_E 最终基线 0.5859–0.9985（均值 0.9375），且随居中度单调上升（hu101 0.586 → ht1_004 0.9985）**——"居中度体现不出来"的现象在当前口径下已不复存在。凡引用历史报告必须区分"修复前口径"与"仍有效的机理/文献/口径分析"。

---

## A. 各文件核心结论速览

### A1. 居中度敏感性微弱根因诊断_2026-08-23.md（B 类，部分过时）
- **原文状态声明（2026-09-02，照抄要点）**："本报告**部分结论已过时**。'居中度弱响应'现象本身是 B1 守恒 bug 口径下的产物——修复后效率随居中度显著单调；三机制中 e 截断已被 A4 证伪、弥散钉死被 B2 证伪为无影响，仅'η_E 间隙加权稀释窄边信号'的指标层分析仍有效。"
- 初诊（hu101 三套居中度剖面实跑）：实测扶正器间均值 0.52 替换假设 0.43 后 η_E 仅 +0.32%；判三根因＝①e 截断 [0.05,0.55] 饱和（`annulus_d2dga.py:303`）②η_E=∫b·c/∫b 间隙加权稀释窄边（贯穿窜槽时 η_E 仍 78.5%）③弥散 mixing≈0.59 钉死基线。
- 伴生发现：失稳指数饱和恒 1.0（hu101/hu102）；c_min=0.05 壁面永久冻结形成"窄边饿死→冻结→永远饿死"正反馈（窄边水泥终态仅 0.008）；η_N（窄四分位效率，hu101 0.09 vs η_E 0.52）被 try/except 吞掉不进 summary。
- 已排除误判：环空入口均匀注入**不是 bug**（鞋口出流单股均质+偏心速度场自然分流，物理等价）。
- 改进建议 P0-1 去截断（后被 A4 证伪）/P0-2 η_N 并列（已由 M0 实施）/P1-1 弥散重标定（已实施后被 B2 证伪非承重）/P2-1 真实 Hele-Shaw 横流（**未做**）/P2-2 c_min 解冻（已由 B2 取代解决）/P2-3 窄缝穿透二值门槛（未做）。

### A2. 居中度敏感性根因再诊断_A4证伪与b3分流放大_2026-08-23.md（B 类，部分过时）
- **原文状态声明（2026-09-02）**："'η_E 对居中度弱响应是 Hele-Shaw 深度平均本质、非 bug'的重新定性被修复推翻——修复后居中度显著单调；b³ 分流机制分析仍有参考价值，但其在修复后口径下的权重需重验。"
- A4 数值试验（e 上限 0.55→0.999，已回退）：e 去截断对 at_cent（高居中）零作用、对 between_cent 反转为 −1.41%（窜槽 0.992，front_narrow=0）→ **"e 截断是死区主因"被证伪**；均居中度 0.52 掩盖了 standoff 低尾（下部 0.22）。
- 证据链 D5：between (b_wide/b_narrow)³=529.7（baseline 77.2 的 6.9 倍），b_narrow=8.4mm 非除零伪影 → 定性为**层流 b³ 本构在高偏处外推失真**（真实有湍流/屈服死区/横向摩擦缓和）。
- η_N 是唯一随居中度单调的指标（0.056/0.090/0.180）；c_min 死锁是窜槽的二级放大非 first cause。
- 建议：e 上限护栏（已做 M4=0.90）、b³ 非层流修正（已做 M2 但对 8 井中性）、c_min 解冻（已做 B2）、真实横流救回（未做）、低尾辅助指标（已做 M0）。

### A3. 环空顶替失真修正与居中度敏感性_调研报告_2026-08-23.md（B 类，文献与机制栈仍是本主题最全的一份）
- **原文状态声明（2026-09-02）**："M1'弥散钉死是主要数值失真'的权重判断被 B2 证伪；M0 指标层方向有效且已实施；M3 屈服门槛方向有效，但已被 B2 物理屈服门取代升级。"
- 核心产出一：**六层机制栈与主次重排**——屈服门槛（M3"开关层"）> b³ 流态分流（M2"记账层"）> 浮力横流/I3（M5"救援层"）> 弥散（M1）> 指标（M0）> 输入数据（M6 standoff/流变）。
- 核心产出二：**流态判定**——hu101 宽边 Re 2400–20000（过渡-湍流）、窄边 Re<600（深度层流），纯层流 b³ 两侧同律失真；幂律层流指数 2+1/n≈3.39 **大于** 3（"b³ 过陡"只在宽边湍流缓和方向成立，禁全局降指数）；宽窄流量比闭合经验式**文献不存在**，需自研。
- 核心产出三：**Lockyear 屈服门槛闭环**——呼101 无浮力项门槛 standoff≈0.72、含浮力（Δρ=140）≈0.24、YP 温度减半 0.12–0.36；恰好解释 A3/A4 三 case（0.43/0.52 未跨门槛→效率本不该大涨；0.73 跨过→η_N=0.18）。
- 核心产出四：修复后居中度响应预期为"曲棍球杆"阈值型（≳0.67 边际递减是正确物理；0.4–0.67 过渡；≲0.4 尖锐塌陷），与 API 10D 67%/SY/T 5334 66.7%/Foolad 工业准则 e≤0.5–0.6（尾管 0.4）/Z23 e≥0.8 死通道四方互证；**主指标应切换 η_N+评价窗效率+窜槽三件套**。
- 文献栈最全（§3 + 参考文献节），含"湍流救不了窄边"三方反证（Couturier/M-F19/Foolad）与检索盲区 7 条（诚实边界）。

### A4. 环空失真修正总结.md（B 类，实施记录真实）
- **原文状态声明（2026-09-02）**："M1 弥散 dt 归一经 B2 证伪为非承重项；旧 c_min/M3 冻结叙事已被 B2 物理屈服门取代；文中 8 井重跑结果数字为 B1 修复前口径，已被取代。修复实施记录本身仍为真实历史。"
- 分支 `annulus-distortion-fix` 合并 `8e0e5ae`（2026-08-24）：M0（η_N/失稳去饱和/评价窗效率/低尾指标进 summary）+M1（弥散×dt/dt_ref）+M2（regime_closure.py，M-F17 式 58–66）+M3（可重启屈服门 ref_mask 不变量）+M4（e_clip_max 参数化）+I3 局部化（np.maximum 修数组崩溃）全部落地，默认关、门关路径逐位复现。
- 评审抓出 6 个真 bug（含 M3 参考元自冻结→全域 wall=100% Critical bug）；332 测试过。
- 8 井"修正后"数字（0.38–0.66 区间）是 B1 修复前口径，**全部作废**；hu103/hu1 当时 η_E≈0 是 c_min 全域冻结病态值。

### A5. 文献下载清单_环空失真修正调研_2026-08-23.md（C 类，文献清单工具页，仍有效）
- 29+ 项分 8 区：A 本地已核 18 篇（两个库存）；B P1 直接支撑 9 篇（标题/DOI/SPE 号已核）；C 经典谱系 11 篇（SPE 号已核）；D 高价值补充 8 篇；E 开放获取 10 条（含免费直链）；F 中文 9 篇；G 标准 3 项；H 追加 7 篇（Lee 1986 等 67% 谱系）。
- 已更正两处引用错误：Jung & Frigaard 2022 DOI 是 **109622**（109581 是 Renteria 2022）；Sauer 1987 是 **SPE-14197**（非 14797）。
- 明示"A2 三篇此前未被分析使用"：Tehrani 1993 实验原文、Foroushan 2020 IADC/SPE-199553 批判论文、Tardy & Bittleston 2015 活动套管模型——均在 `cankao/` 本地。

### A6. 顶替效率模型参考文献调研_改进方向.md（C 类，6 篇 Frigaard 组精读）
- 时态声明（2026-08-19 审计）：Tier 0 七项诊断、Tier 1 八项物理 **已全部完成**（T1-1~T1-8，含 T1-10 混浆、T1-11 Atwood 浮力）；**Tier 2 五项重构仍待做**（HB 闭包查表 T2-1 / augmented Lagrangian T2-2 / TVD T2-3 / 瞬态 T2-4 / 序列优化器 T2-5）。
- 六篇谱系：Pelipenko04（行波+Muskat 判据，纯后处理低成本）、Moyers-González07（瞬态+脉动+停泵有限时间衰减）、Yang21（序列优化+偏心-浮力平衡 e=0.1 最优、e>0.3 滞留）、Zhang22（D2DGA 奠基）、Zhang23（实验+三分类判据+浮力数阈值 b>80~100 non-dispersive、e≥0.8 窄边显著下降）、Bararpour25（HB 推广+augmented Lagrangian+static wall layer）。
- 待核实项：**I3 公式分母差因子 6/m**（代码 `2m[...]` vs B25 式 `12[...]`，三份文件反复标注高优先级）；m 定义方向需确认一致性。
- 注意：其中"浮力体力向量死代码"表述已过时（T1-3b 注入 R2 I3 通量，非死代码）。

### A7. 尾管固井顶替效率模型合理性评估与改进方向_独立文献调研.md（C 类）
- 独立网络调研（2026-07-27），结论：**框架层面合理**（D2DGA 是国际主流，核心公式均有严格文献推导）；最需改进 3 项＝HB 闭包（最大物理差距）、流动分类判据（纯后处理高回报）、I3 分母核实（可能偏差数倍）。
- 与内部调研一致性 >90%；补充维度：井径不规则性（方春飞 2016：扩大率>1.1 效率显著下降）、套管旋转（AADE-18-FTCE-094/Foroushan21）、水泥污染强度（Khan 2024：30% 污染→抗压强度 −72.6%）、MPC 控压固井。
- "拉普拉斯平滑在 I3 启用后应降级为纯稳定器"的判断被 B2 证实（清理清单明示）。
- 08-19 审计：比对表中"部分实现"项多数已由 Tier 0/1 完成；Tier 2 截至 08-19 仍待做。

### A8. 西南石油大学固井顶替科研成果调研汇总_2026-08-31.md（C 类）
- **状态声明**：网络检索口径（tavily，60+ 轮搜索），每条带 URL 佐证，**未经人工逐篇全文复核**，引用前需核对原文版面；含同名陷阱警示（中石化西南石油工程≠SWPU；陈家琅/王瑞和等均非 SWPU）。
- 三代谱系：刘崇建/郭小阳→李早元/徐璧华/邓建民→孙劲飞/杨谋/王敬朋等；四条主线（李早元 CFD+VOF、杨谋评估模型+ML、邓建民密度差、徐璧华控压固井）。
- 与本课题最同源：王敬朋 2022《呼探1井 Φ139.7mm 尾管精细动态控压固井技术》（同区块同井，必引对标）；★3 五力平衡滞留层模型（ACS Omega）与本项目 b³ 分流/屈服门槛修正直接相关。
- SWPU 三处空档（论文差异化定位依据）：停泵后顶替/静置失重/U 型管、控压固井全过程温压场瞬态建模（石大北京+塔里木主导）、紊流临界机理；ML/PINN 线已被 SWPU 抢先（P1/P2 专利），数据驱动不宜作主卖点。
- 37 条中偏心/顶替数值模拟相关条目见 B6 节收割。

### A9. 结果偏低根因独立调研报告_2026-09-02.md（B 类）
- **原文状态声明（2026-09-02）**："其独立发现（hu103 管容双链分裂、erf 92% 封顶、过填归一化不收敛）真实有效并已由 B1/B2/B3 修复落地；但'修复后预期 0.62–0.68、hu103 回到 0.40–0.60'的预期与'F4/D 项为主要残留缺陷'的权重判断已被 B1 修正——漏判了截面因子 2 守恒 bug 这一主因。"
- 决定性定性：**η_E ≡ 水泥库存比**（`eff=cement`，与 bulk_cement_fill 逐位相等，8 井 15 位有效数字验证）＝"评价时刻已入环空水泥体积/环空域体积"，不是文献意义的"顶替效率"；上限=水泥总量/环空体积；2D 位移物理只能往下压它。
- hu103 管容双链分裂实锤：fronts 链 71.1 vs timeline 链 115.1 m³（差 61.9%），停止早于尾浆到达 699.7s，尾浆供给 0 m³；工艺自洽核算真实管容 ≈82.4（0708 后定值 88.55）。
- 质量守恒三通道：F4 erf 92.07% 封顶（无 frac→1.0 收尾）、F5 顶部反射边界+过填归一化销毁、D 项 nz 加密 η_E 单调下降（0.290/0.231/0.194）=不收敛指纹。
- 论文口径建议：η_E（库存）/入库完成度/η_N（窄边）三行并列声明；禁止写模型-现场相关系数。

### A10. 结果偏低复核_管容双链分裂根因_2026-09-02.md（B 类）
- **原文状态声明（2026-09-02）**："管容双链分裂的发现仍有效（已由 B3 统一双内径感知管容修复落地）；但'修复后预期 0.62–0.68'与'F4/D 项为主要残留缺陷'的权重判断已被 B1 发现修正。"
- 关键增量：分裂是 **08-29 loader 校准后显形**的（timeline 链管容被改大，fronts 链未动）；08-16 的 0.60 是分裂未显形时的旧态。
- 逐井"下降"核查：hu102 ③④差 31.7pp 是 runner 基线口径 c_min 冻结偏重极例；ht1_004 0.1413 是 +600s 尾窗失稳病态值；修复后 CBL 排序一致（78.7>66.65>62.77>2.12>0.3）。
- 次级发现：timeline 链尾缘配对未同步 F2 修复（ht1_001 −117.8s，09-01 已修）。

### A11. 守恒修复与GIF前缘定论_模型状态及论文数字路线_2026-09-02.md（**C 类，本轮终局权威报告，无作废声明**）
- GIF"水泥爬不到顶"定论：三层根因 B1（主因，半速前缘自 2026-05-07 fa44ace 起存在 4 个月）+B3+B2 全部已修；修复后 8 井宽边到位率全 1.000，窄边 hu101 仅 2.0%（=真实窄边窜槽，非伪影）、hu103 77.9%、hu1 93.6%。
- 八井权威表（B1/B2/B3 后）：0.5348/0.9687/0.8529/0.9494/0.9843/0.9888/0.9976/0.9985；0708 管容核实后终跑 0.5859/0.9846/0.9926/0.9659/0.9846/0.9898/0.9984/0.9985；胶塞语义轮（09-03，论文最终基线）仅 hu102 微调 0.9836（−1.88 m³ 滞留尾浆精确闭合）。
- **效率随居中度单调升**（hu101 0.438→ht1_003/004 0.83 区间 η_N），与文献量级吻合；hu101 是唯一低值=真实窄边窜槽（η_N=0.03、29.7% 未屈服冻结层）。
- 论文叙事硬要求：η_E 修复后无区分度，**验证叙事必须换轴**——顶替充分性轴（η_E/η_N）+候凝风险轴（失稳/窜槽），两锚点井 ht1_003（好）/ht1_004（差）两轴自洽；η_E 定义声明必须写；历史表 4/08-16/08-24/08-29/09-02 各轮全不得进正文。
- 遗留 7 条（见 C 节）。

### A12. 固井顶替模型因果分析与偏心度敏感性.md（C 类，数字已挂过时声明）
- 08-19 审计声明：偏心度敏感性数值基于 2026-07-14 旧基线（CFL 前）"仅作定性参考"；`_buoyancy_force_vector` 非死代码的判断仍正确。
- HT1-004 七点 standoff 扫描（nz=500，dt=4s 旧口径）：η_E −0.27/e、失稳 +2.80/e、窜槽 +2.13/e；**敏感性排序：失稳>窜槽>>混浆>>效率**；e=0.55 硬饱和（standoff 0.40/0.30 结果逐位相同）。
- 定性结论"居中度是头号可控因素、失稳指数最早亮红灯"被修复后单调性证实（清理清单 C 类理由）。
- 套管 1D 定位：解析体积推进（活塞流前缘追踪）非 PDE；对环空的影响走五通道（排量时序/相分数时序/总时长/管容延迟/重力弥散）。

### A13. 模型现场差异归因与效率影响因素调研_2026-08-20.md（B 类）
- **原文状态声明（2026-09-02）**："'弥散/混浆未接线'已由 08-22 接线修复、'CFL 开关差 30pp'已裁定收口，弥散相关权重被 B2 证伪；**模型-CBL 结构性口径差异与'standoff 为最大输入不确定性'的结论仍有效**（修复后居中度/standoff 通道反而更加重要）。"
- 结构性必然：模型=全求解域×泵停时刻×体积分数均值 vs CBL=评价井段×候凝 24h+×胶结质量——深度段/时间/物理量三重错配，不能当"误差"读。
- standoff 是最大输入不确定：**8 井全部 model_assumption 无一实测**；实测敏感性 HT1-001 standoff 0.63→0.83 时 η_N 0.379→0.593（+21pp）。
- CBL 对照链路已断（validation/ 07-29 删除）；候选因素五层清单（口径/输入/物理简化/数值/现场未覆盖物理）仍然有效。

### A14. D2DGA三闭包公式速查.md（C 类，公式参考不受修复影响）
- 三闭包 LaTeX 全式（R1 auto-m 黏度比场、R2 I3 浮力弥散通量、R3 真浮力体力+式 4.24 方位修正）+ 通量放大 f(c,m)、I1/I2、两层黏度闭包，符号表与代码位置（`d2dga_flux.py`）。**供主编对照代码用**；内嵌警示：I3 分母 6/m 差异待与原文核对；m 定义 Yang & Yortsos 1997。
- 公式细节我不照抄进结论——具体每式与代码逐项对照**待主编对照代码**（本摘要只提供出处定位）。

### A15. 停止标志修复与环空域口径实验_2026-09-01.md（B 类）
- **原文状态声明（2026-09-02）**："停止标志语义修复（stop≡碰压）与环空域口径结论仍有效；但 §1.3 尾窗效应、§2.2 裸眼域对照与 §4 追记的 η_E 数字均为 B1 修复前口径，绝对值已失效。"
- stop≡碰压：8 井新 stop 与旧扫描逐位一致（Δ=0），语义从"碰巧对"变"按构造对"；尾窗去除效应符号逐井不一（hu103 反向 −12pp），"尾窗一律有害"证伪。
- 环空域口径定论：**8 井域全部含技套内重合段**（hu1 域 2/3 在技套内）；重合段井径口径偏差 13–28%（hu103/ht1_001/hu101）；建议**修井径不裁域**；CBL 评价窗全在裸眼段→域口径对模型-CBL 对比影响有限。
- 追记（09-02）：wall 全冻结死锁缺陷实锤（弥散光晕 1e-20 级也判"已到达"+pref 下限顺序错误→注入死锁，hu103 80s 内 171/250 列全冻 η=0.0013）；C+A 修复落地（wall_seed_c_min=0.005 + 全冻列保流通道）。

### A16. 2026-09-02报告归档清理清单.md + 归档目录（清单本身 C 类）
- A 类归档 4 份（`归档_已被2026-09-02守恒修复取代/`）：校准后八井重跑分析_08-29、尾管顶替效率模型现状与改进空间_综合调研_08-29、模型与参考文献差距分析及优化路线_08-23、**结果下降与居中度敏感性_综合修复调研_08-31**——核心结论均被 B1/B2/B3 推翻。
- B 类原位加声明 11 份；C 类不动 20 份。判定标准="核心结论/修复建议是否被 B1/B2/B3 推翻、取代或证明为次要"。
- 归档内《结果下降与居中度敏感性_综合修复调研_08-31》虽整体归档，但其 Sources 21 条（[S1]–[S21]）与"不做/已排除清单"（❌e 去截断当解药、❌b_narrow≥5mm 下限、❌"湍流救窄边"普适化、❌全局降幂律指数、❌刮泥器建模、❌旧表 4 修回）仍有文献价值，已并入本摘要 B 节。

---

## B. 已核实文献清单收割（按主题分组，逐条：作者(年份). 题名. 出处. + 支撑论点 + 来源文件）

> 通用说明：标注"本地已核"=两个本地库存（`参考文档/顶替效率模型参考文章/`、`cankao/`）已核对；标注"未全文核读"=仅 tavily 定位；标注"⚠️出处冲突"=不同 obsidian 文件给出的卷期不一致，引用前需核对原文。中文文献多为"引用页佐证"，SWPU 汇总明示未经人工逐篇复核。

### B1. Hele-Shaw / D2DGA 框架与闭包（回答"1D+2D 模型有无问题"的理论基准）
1. Bittleston, S.H., Ferguson, J. & Frigaard, I.A. (2002). Mud removal and cement placement during primary cementing of an oil well — Laminar non-Newtonian displacements in an eccentric annular Hele-Shaw cell. J. Eng. Math. 43:229–253.（偏心 Hele-Shaw 顶替开山之作、D2DGA 前身、商业软件基准）｜本地 cankao/｜出自：文献下载清单 A2。
2. Zhang, R. & Frigaard, I.A. (2022). Primary cementing of vertical wells: displacement and dispersion effects in narrow eccentric annuli (Part 1). JFM 947:A32.（D2DGA 奠基；通量放大 f、I3、两层黏度闭包、浮力体力式 4.13/4.24 的直接来源）｜本地已核｜出自：三闭包公式速查、独立文献调研、下载清单 A1。
3. Zhang, R. & Frigaard, I.A. (2023). Part 2. Flow behaviour and classification. JFM（⚠️出处冲突：下载清单 A1 写"JFM 972 A38"，独立文献调研写"JFM 976 A24"——需主编核对原文卷号）.（~120 次实验+3D+D2DGA，三分类判据、浮力数阈值 b>80~100、e≥0.8 窄边效率显著下降、e 不改变分类但增大残量）｜本地已核｜出自：顶替效率模型参考文献调研_改进方向 §1.2。
4. Bararpour, S. & Frigaard, I.A. (2025). Capturing dispersion of Herschel–Bulkley fluids in miscible primary cementing displacement flows.（⚠️出处冲突：下载清单/守恒修复报告写"JFM 1022 A15（开放获取全文）"，独立调研附录取舍"JFM"，三闭包速查写"J. Non-Newtonian Fluid Mech. 2025"——需核对）.（D2DGA 推广到 HB；augmented Lagrangian 真屈服；static wall layer c̄min 判据=替代 c_min=0.05 的文献依据；Muskat 三 regime）｜本地已核｜出自：下载清单 A1、独立文献调研 §1.1。
5. Pelipenko, S. & Frigaard, I.A. (2004). On the displacement of Herschel–Bulkley fluids in narrow eccentric annuli / 行波解+Muskat 稳定性。JFM 520:343–377（归档[S8]口径）；另两篇 Part 系列 J. Eng. Math.（下载清单 A2：Part 2 steady-state J. Eng. Math. 48:1–26；独立调研附录写 J. Eng. Math. 50:125-158 ⚠️卷期多口径需核对）.（Muskat 宽/窄边指进判据纯后处理成本极低；行波三态=稳态/非稳态/静态泥浆通道；ELF 工业规则保守 3–10%）｜本地已核｜出自：改进方向 §1.4、归档综合修复调研[S8][S14]。
6. Moyers-González, M., Frigaard, I.A. et al. (2007). Transient displacement flows in primary cementing. JFM 576:31–59.（流函数三阶非线性 PDE 瞬态演化；停泵有限时间衰减判据式 3.35/3.40；**准静态模型无法捕捉脉动**——当前模型正是准静态）｜本地已核｜出自：改进方向 §1.5。
7. Moyers-González, M. & Frigaard (2009). Kinematic instabilities in two-layer eccentric annular flows, part 2: shear-thinning and yield-stress effects. J. Eng. Math. 65:25–52.（界面运动不稳定，切稀/屈服）｜本地 cankao/｜出自：下载清单 A2。
8. Yang, Z. & Yortsos, Y.C. (1997). Asymptotic solutions of miscible displacements in porous media.（黏度比 m 定义 m=μ_displaced/μ_displacing 的出处）｜出自：D2DGA三闭包公式速查。
9. Izadi & Frigaard (2026). GSE（本地已核；squeeze cementing/微环空黏塑性侵入流方向，a 后续应用）｜出自：下载清单 A1、独立文献调研 §2.1。

### B2. 偏心环空实验基准（居中度敏感性"应该长什么样"的证据）
10. Foolad, Bizhani & Frigaard (2021). 强偏心环空顶替实验（Carbopol）。Energies 14(6):1654. https://www.mdpi.com/1996-1073/14/6/1654 （OA 已全文提取）。**四条与本课题直接同构的结论**：①e=1.0 时残层占间隙 25–55% 而体积效率仍 70–91%——体积指标"具有欺骗性"（直接支撑 η_E 降级/η_N 并列）；②湍流（水）顶替下窄边几乎不动——"湍流扩散更有效"在强偏心不成立；③宽边突破后靠缓慢二次流清除，接触时间规则失效；④转引工业准则 e≤0.5–0.6（**尾管 0.4**）｜出自：环空失真修正调研 §2.4。
11. Bizhani, Foolad & Frigaard (2020). Turbulent displacement flow of viscoplastic fluids in eccentric annulus: Experiments. Phys. Fluids 32:045117.（强偏心湍流顶替实验，Foolad 2021 上游）｜出自：下载清单 B#5、归档[S20]。
12. Bizhani & Frigaard (2020). Buoyancy effects on turbulent displacement of viscoplastic fluids from strongly eccentric horizontal annuli. Phys. Fluids 32:125112.（浮力对窄边清除的定量作用——I3/横流修正依据）｜出自：下载清单 B#6。
13. Malekmohammadi, Carrasco-Teja, Storey, Frigaard & Martinez (2010). An experimental study of laminar displacement flows in narrow vertical eccentric annuli. JFM 649:371–398.（**垂直偏心环空层流顶替实验基准**，与本项目井型一致，Z22/Z23 上游）｜出自：下载清单 D#1。
14. Escudier, Oliveira, Pinho & Smith (2002). Fully developed laminar flow of non-Newtonian liquids through annuli. Exp. Fluids 32（免费直链 https://web.fe.up.pt/~fpinho/pdfs/ExpFluids2002.pdf ，已全文精读）.（e=0.8 宽边 LDA 实测与层流计算几乎一致、窄边仅定性——**层流侧不改形状、只改流态判定**的实验锚点；原文无 b³ 直接评价，只能间接背书）｜出自：环空失真修正调研 §3.1。
15. Walton & Bittleston (1991). The axial flow of a Bingham plastic in a narrow eccentric annulus. JFM 222:39–60.（b³ 律+屈服修正的解析基准——M2/M3 与单元回归的对照解）｜出自：下载清单 B#7、环空失真调研 §3.1。
16. Lund, Taghipour, Ytrehus & Saasen (2020). SINTEF 环空实验（ε=−0.42）。Energies 13(19):5201. https://www.mdpi.com/1996-1073/13/19/5201 （OA 全文）。**旋转使窄边与宽边同时或更早被顶替**；偏心→强倾斜界面、flooding wave 充填｜出自：环空失真调研 §3.3。
17. Ytrehus et al. (2017). ASME OMAE-62028.（SINTEF 系列效率-偏心定量数据最多的姊妹篇）｜出自：下载清单 D#7。
18. Zhang, Jung, Renteria & Frigaard (2022). Experimental study of displacement flows in a vertical eccentric annulus. ASME OMAE V010T11A051.（Z22/23 作者群 OMAE 实验篇）｜出自：下载清单 D#8。
19. Renteria & Frigaard (2020). Primary cementing of horizontal wells... Part 1. Experiments. JFM 905:A7；Sarmadi, Renteria & Frigaard (2021) Part 2. Computations. JFM 915.（水平偏心实验+数值对；D2DGA 家族横向对照）｜出自：下载清单 D#2/D#3。
20. Renteria et al. (2022). Effects of wellbore irregularity on primary cementing of horizontal wells, Part 1. JPSE 208:109581.（井眼不规则/扩径缩径效应——hu101 井径剖面相关）｜出自：下载清单 D#4（DOI 109581 曾被误挂给 Jung 论文，已勘误）。
21. Renteria, Maleki, Frigaard, Lund, Taghipour & Ytrehus (2019). Effects of irregularity on displacement flows... highly deviated wells. JPSE 172:662–680.｜出自：下载清单 D#5。
22. Skadsem, Kragset, Lund, Ytrehus & Taghipour (2019). Annular displacement in a highly inclined irregular wellbore: Experimental and three-dimensional numerical simulations. JPSE 172:998–1013.（**3D 数值+实验**，宽窄边分配的 3D 对照、M2 验证素材）｜出自：下载清单 D#6。
23. Tehrani, Bittleston & Long (1993). Flow instabilities during annular displacement of one non-Newtonian fluid by another. Exp. Fluids 14:246–256.（**层流顶替最小密度差 10–15% 规则的实验原文**；hu101 领浆/泥浆仅 7% 略低于该值——论文可讨论点；本地 cankao/，综述转引可升级为直接引用）｜出自：下载清单 A2、环空失真调研 §3.2。
24. Wu et al. (2023). Theoretical and Experimental Study on Cementing Displacement Interface for Highly Deviated Wells. Energies 16:733, doi:10.3390/en16020733.（83° 井斜界面理论与实验；e=0.4 降 9.4pp 量化锚点）｜本地 cankao/｜出自：下载清单 A2、归档[S4]。
25. Sun, X. et al. (2020). Numerical modeling of motion of displacement interface in eccentric annulus. Energy Sci. Eng. 8:1579–1591（⚠️注意与 B6 Sun J(孙金斐)2020 同名不同文，此条为 cankao/ 存量 Energy Science…Sun…pdf）｜出自：下载清单 A2。

### B3. 屈服门槛 / 静态残泥 / 窄边窜槽判据（"为什么窄边顶不动"的物理层）
26. McLean, Manry & Whitaker (1967). Displacement Mechanics in Primary Cementing. JPT 19(3):251–260, SPE-1488.（**McLean 规则原文**：水泥屈服强度>泥浆屈服强度×最大/最窄间隙比；水泥屈服须低到能进窄缝）｜出自：下载清单 C#1、环空失真调研 §3.2。
27. Lockyear & Hibbert (1989). Integrated Primary Cementing Study Defines Key Factors for Field Success. JPT 41(11):1320–1325.（standoff 效应+破胶判据 (b/2)Δp>τ_g Eq.1、含浮力项 Eq.2——M3 Lockyear 判据谱系）｜出自：下载清单 C#3。
28. Lockyear, Ryan & Gunningham (1990). Cement Channeling: How to Predict and Prevent. SPE Drill. Eng. 5(2):201–208.（窜槽预测/预防规则原文）｜出自：下载清单 C#4。
29. Couturier, Guillot, Hendriks & Callet (1990). Design rules and associated spacer properties for optimal mud removal in eccentric annuli. SPE-21594.（**四规则原文：①密度差≥+10%；②摩阻梯度≥+20%；③窄边剪应力须超过被顶替液屈服应力；④宽边顶替液速度≤窄边被顶替液速度**；"全周湍流才有益"判据）｜出自：下载清单 C#2、M-F19 精读确认。
30. Roustaei, Gosselin & Frigaard (2015). Residual drilling mud during conditioning of uneven boreholes... Part 1: Rheology and geometry effects in non-inertial flows. JNNFM 220:87–98.（**临界压差动员判据；无屈服应力则无静态残泥只有 Moffatt 涡；e=0.6±0.2 定量图版**；两篇均无浮力项——浮力增强须标注为模型扩展）｜本地开放直链已核｜出自：下载清单 B#8、环空失真调研 §3.2。
31. Roustaei, Gosselin & Frigaard (2015). Part 2（惯性警告）. JNNFM 226:1–15.（**固定 He 下中 Re 段加大排量反而使静区增大——"泵快=洗得净"不成立；"动了≠清走"**（回流区只剪不置换））｜已全文精读｜出自：环空失真调研 §3.2。
32. Wang et al. (2024). Processes 12(6):1176. https://www.mdpi.com/2227-9717/12/6/1176 （OA 全文提取）。泥浆滞留判据=泥浆屈服>水泥屈服 且 壁面剪应力<泥浆屈服；层厚随泥浆 τy/井斜增大、随水泥 τy/密度差减小；e=0–0.3 使宽窄滞留差增大｜出自：环空失真调研 §3.2。
33. Foroushan, Ozbayoglu, Gomes & Yu (2020). Mud/cement displacement in vertical eccentric annuli. SPE Drill. Complet. 35(2):297–316.（垂直偏心逐截面流量分配工程实现；**常数表观黏度使窄边速度高估 43%（e=0.5）/68.6%（e=0.75）**——局部流变修正必须与流态修正一起做；批评"逐切片+同心动量方程"解耦做法）｜已全文精读｜出自：下载清单 B#1、环空失真调研 §3.1。
34. Foroushan, Ozbayoglu & Gomes (2020). How Realistic is the Calculated Cementing Displacement Efficiency? IADC/SPE-199553-MS.（**"计算顶替效率有多真实"批判论文**，与 η_E 欺骗性主题直接同构；本地 cankao/）｜出自：下载清单 A2。
35. Haut & Crook (1979). Primary Cementing: The Mud Displacement Process. SPE ATCE.（泥浆不流动时密度差无用的实验证据）｜出自：下载清单 C#6。
36. Clark & Carter (1973). Mud Displacement with Cement Slurries. JPT 25:775–783, SPE-4090.（偏心/接触时间经典缩尺实验）｜出自：下载清单 C#5。
37. Brice & Holmes (1964). Engineered Casing Cementing Programs Using Turbulent Flow Techniques. JPT 16:503–508, SPE-742-PA, doi:10.2118/742-PA.（**接触时间 10 分钟规则源头**；Foolad 2021 证明强偏心下失效）｜出自：下载清单 C#9。
38. Jakobsen et al. (1991). Displacements in Eccentric Annuli During Primary Cementing in Deviated Wells. SPE Prod. Ops. Symp.（**密度差帮助清除低边窄间隙滞留泥浆的实验原文**——浮力救援论证）｜出自：下载清单 C#7。
39. Sauer (1987). Mud Displacement During Cementing: A State of the Art. JPT 39(9), SPE-14197（已勘误）.（经典判据汇总）｜出自：下载清单 C#8。
40. Daccord, Guillot & Nilsson (2006). Mud Removal. In: Nelson & Guillot (eds), Well Cementing, 2nd ed., ch.5, pp.143–189, Schlumberger.（行业教科书级权威综述）｜出自：下载清单 C#10。
41. Kragset & Skadsem（年份/出处不全，待主编补核）. 屈服应力偏心环空顶替、稳态解非唯一.（屈服正则化失真机理线索）｜出自：归档综合修复调研[S10]。

### B4. 流态分流 / 湍流-层流 / 数值方法（回答"是否数值计算问题"）
42. Maleki & Frigaard (2017). Primary cementing of oil and gas wells in turbulent and mixed regimes. J. Eng. Math. 107:201–230.（已全文精读。**M2 闭包来源：局部 Metzner-Reed Rep 式 58+Hedström 修正式 59；层流闭式 60–61；屈服静止写进流动度式 63；Dodge–Metzner 全湍流式 65；过渡区摩擦因子对数插值**；"宽边湍流、窄边层流甚至静止可同截面并存"；实现硬约束=固定点自洽迭代防"宽边自证湍流"）｜出自：环空失真调研 §3.1。
43. Maleki & Frigaard (2019). Comparing laminar and turbulent primary cementing flows. JPSE 177:808–821, doi:10.1016/j.petrol.2019.02.054.（已全文精读。Couturier 四规则原文+定量校准锚点 **e=0.3→η_N≥95%；e=0.6→η_N≈30–35%**；η_N 定义（最窄 1/4 方位 H 加权）与本项目 displacement_metrics 完全同构；湍流不必然更优——高偏心全湍流 η_N 最高仅 ~35%，等密度湍流 η_N=0）｜出自：环空失真调研 §3.1。
44. Maleki & Frigaard (2018). PoF 30:123101（下载清单 B 节列出）；Maleki & Frigaard (2020). Chem. Eng. Sci. 219（流态快速分类，归档[S12]）｜出自：下载清单 B 节、归档[S12]。
45. Boniou, Schmitt & Vié (2022). 界面捕捉格式对比. Int. J. Multiphase Flow（HAL 预印本 https://hal.science/hal-03241460 ，OA 全文）。**代数 VOF"受数值扩散困扰"；欠分辨结构产生"数值雾化"——直接印证"常数拉普拉斯平滑=过弥散"病灶定性**｜出自：环空失真调研 §3.4。
46. Sweby (1984). SIAM J. Numer. Anal. 21(5)（限幅器族 φ(r)）；Boris & Book (1973)、Zalesak (1979)（FCT 原始文献）。（TVD/FCT 实现参考；建议 MC 限幅器 θ≈1.4–1.5 起步）｜出自：环空失真调研 §3.4、改进方向 T2-3。
47. Ezeakacha/Salehi 等 (2021). gap-average 模型系统批评综述. Energies 14(3):573（⚠️出处冲突：主调研报告把同一卷期 14(3):573 归给 Foroushan et al. 2021 综述《Cement Placement: An Overview of Fluid Displacement Techniques and Modelling》——**两处必有一误，主编核对**）｜出自：归档[S7] 与 环空失真调研 §参考文献。

### B5. 居中度标准 / 现场准则（阈值型响应的规范锚点）
48. API Spec 10D.（弓簧扶正器规范；**67% SOR 是规范值非物理推导**，存在"Centralization Myth"反方——模型不得向 66.7% 校准）｜出自：环空失真调研 §3.3、下载清单 G。
49. SY/T 5334-1996《套管扶正器安装间距计算方法》.（**居中度>0.67 合格线**国内出处）｜出自：环空失真调研 §3.3、下载清单 G#1。
50. Lee, Smith & Tighe (1986). Optimal spacing for casing centralizers. SPE Drill. Eng. 1(2):122–130.（**67% standoff 规则的原始出处**，Jung&Frigaard 精读后锁定）｜出自：下载清单 H#1。
51. Jung & Frigaard (2022). Evaluation of common cementing practices affecting primary cementing quality. JPSE 208, doi:10.1016/j.petrol.2021.109622（⚠️已勘误，109581 是 Renteria）。（67% 谱系 Lee86→Bottiglieri14→Adesanya18；10 分钟接触时间规则仅基于现场观察）｜出自：下载清单 B#2/H 节。
52. Gorokhova, Parry & Flamant (2014). Soft vs stiff-string centralization. SPE Drill. Compl. 29(1):106–114；Guillot et al. (2008). Are current casing centralization calculations really conservative? SPE-112725；Sanchez, Brown & Adams (2012). SPE-150317.（居中度计算方法学三件套）｜出自：下载清单 H#2-4。
53. 丁保刚等 (2009). 套管居中设计与校核. 石油钻探技术 37(1).（SY/T 5334 66.7% 引证与计算方法）｜出自：下载清单 F#2。
54. 埕海油田大斜度井超短尾管固井 (2020). 石油钻采工艺 42(1).（**实测居中度：重叠段 75.4%、裸眼段仅 33.4%**——实测低尾普遍的直接现场证据；密度差≥0.24–0.25 kg/L 才有较好顶替）｜出自：环空失真调研 §3.3、下载清单 F#5。
55. 东海X7 长裸眼全封固井 (2022). 海洋石油 2022(2), doi:10.3969/j.issn.1008-2336.2022.02.077.（居中度>67% 设计实践、半刚性旋流扶正器）｜出自：下载清单 F#6。
56. 长江大学/川庆（刊名待查）. 套管居中度计算模型改进.（**理论居中度 vs 成像测井误差 ~11%——standoff ±0.1 区间标注的依据**）｜出自：环空失真调研 §3.3、下载清单 F#7。
57. 上海交大 (2014). JOES 大位移井.（建议偏心度控制在 0.5 以内；旋流扶正器 60° 螺旋角提效）｜出自：环空失真调研 §3.3。
58. Olkaria 地热井模拟 MSc (2025). skemman.is.（standoff 50% 显著低效、70% 接近完全顶替——单点案例，低-中证据级）｜出自：环空失真调研 §3.3、下载清单 E#9。

### B6. 中文/国内偏心顶替数值模拟谱系（SWPU 汇总收割，均"URL 佐证存在、未人工全文复核"）
59. 孙劲飞, 李早元, 罗平亚, 等 (2019). 水平井偏心环空低速顶替运移机制研究. 西南石油大学学报 41(1):111–118.（CFD+VOF：**e>0.5 时加隔离液无效；严重偏心下降速至 0.2 m/s 反而多置换 6.8% 钻井液**——"降速救窄边"反直觉锚点）｜URL: https://html.rhhz.net/XNSYDXXBZRB/HTML/2019-1-111.htm ｜出自：SWPU 汇总 ★1。
60. Sun J(孙金斐), Li Z(李早元), Luo P, et al. (2020). Numerical modeling of motion of displacement interface in eccentric annulus. Energy Sci. Eng. 8(5):1579–1591, doi:10.1002/ese3.615.（**径向+周向分解为两个 2D 平面建界面运动模型**——与本项目 D2DGA 同属降维谱系的最直接竞品；修正密度差 F*、Bingham 数、黏度比、偏心率对窜槽/残留的影响准则）｜出自：SWPU 汇总 ★2。
61. Wang J(王敬朋), Xiong Y(熊友明), Lu Z, et al. (2022). Research on Key Technologies to Improve Cementing Displacement Efficiency. ACS Omega 7(42):37039, doi:10.1021/acsomega.2c00419（全文 PMC9607674）。**五力平衡滞留层模型**（顶替压力/密度差浮力/隔离液壁面剪切/钻井液屈服应力/粘附力）→滞留层厚度→顶替效率与排量设计——与本项目屈服门槛修正直接公式化对照｜出自：SWPU 汇总 ★3。
62. 王敬朋, 张伟, 吴继伟, 等 (2022). 呼探1井 Φ139.7mm 尾管精细动态控压固井技术. 石油钻探技术 50(6):92–97, doi:10.11911/syztjs.2022021.（**同区块同井同场景，论文综述必须对标**；全过程动态控压+扶正器安放优化）｜出自：SWPU 汇总 ★11、下载清单 F#1。
63. 郑永刚 (1993). 非牛顿流体渗透性井壁偏心环空二维流动. 天然气工业 13(5).（SWPU 偏心环空顶替流体力学源头作；牛顿/宾汉/幂律解析解）｜出自：SWPU 汇总 #6。
64. 杨建波, 邓建民, 黎泽寒, 陈英 (2007). 低速注水泥过程中密度差对顶替效率的影响. 石油钻探技术 35(5)（起页 18）；杨建波等 (2008). 数值模拟研究. 同刊 36(5):62–65.（**η正>η无>η负；正密度差塞流最优；负密度差需提流速补偿**——Atwood 浮力修正 T1-11 中文支撑）｜出自：SWPU 汇总 #20/21。
65. 陈力力, 郭建华, 刘森, …, 杨谋(通讯) (2023). 提升偏心环空注水泥顶替效率的浆柱结构优化分析. 钻井液与完井液 40(1):103–110.（Fluent 四种浆柱次序；偏心 0.1 时充填 93.68%）｜出自：SWPU 汇总 #4。
66. Yang X(杨鑫), Sun J, Zheng G, Li Z, et al. (2025). Numerical simulation of liquid-solid two-phase flow in cementing displacement based on kinetic theory of granular flow. Geoenergy Sci. Eng. 252:213956.（**偏心 0.6 下套管转速 0→30 rpm 顶替效率 88.2%→94.21%**——旋转提效最新中文定量锚点）｜出自：SWPU 汇总 #5。
67. Yang M(杨谋), Che S(车双苗), et al. (2021). Optimizing cementing displacement efficiency by optimizing pad fluid injection sequence. JPSE 204:108691.（偏心 0.4 水平井四相序次 CFD，85.88%；=内部调研 Yang 2021 同篇）｜出自：SWPU 汇总 #15、改进方向 §1.6。
68. Yang M(杨谋) et al. (2025). Optimizing cementing displacement efficiency using a machine learning ensemble models. Geoenergy Sci. Eng. 214149.（RF+Extra Trees+Stacking：宽/窄间隙效率 +3.74%/+5.46%，偏差<3.5%）｜出自：SWPU 汇总 ★25。
69. 李明, 杨雨佳, 李早元, 程小伟, 郑友志, 郭小阳 (2014). 固井水泥浆与钻井液接触污染作用机理. 石油学报 35(6):1188–1196.（混浆增稠/絮凝机理——T1-10 混浆建模中文支撑）｜出自：SWPU 汇总 #22。
70. 高永海, 孙宝江等 (2005). 环空水泥浆顶替界面稳定性数值模拟研究. 石油学报 26(5):119–122.（国内顶替界面稳定性早期数模）｜出自：独立文献调研附录。
71. 方春飞, 周仕明 (2016). 井径不规则性对固井顶替效率影响规律研究. 石油钻采工艺.（**井径扩大率>1.1 时顶替效率显著下降**——模型用平均井径的缺陷依据）｜出自：独立文献调研 §1.3/附录。
72. 郭小阳, 刘崇建, 张明深, 等 (1998). 提高注水泥质量的综合因素. 西南石油学院学报 1998(3):57–59.（居中度/流变/顶替多因素系统论述，谱系源头）｜出自：SWPU 汇总 #29。
73. Sun J(孙金斐), Yang F(杨富杰), Qi B, Li Z, et al. (2024). A New Calculation Model for ECD Considering Interface Effect between Various Fluids during Cementing. SPE Journal 29(5):2242–2256.（**界面形态随时间变化修正 ECD**——"顶替界面-环空压力动态耦合"获国际期刊认可，控压固井差异化卖点印证）｜出自：SWPU 汇总 ★7。

### B7. 套管旋转 / 活动套管（居中度缺陷的工程补救通道，未决 P2 项的文献基础）
74. Tardy & Bittleston (2015). A model for annular displacements of wellbore completion fluids involving casing movement. J. Pet. Sci. Eng.（活动套管顶替模型；本地 cankao/）｜出自：下载清单 A2、归档[S17]。
75. Carrasco-Teja & Frigaard (2009). Phys. Fluids 21:073102；(2010) JFM 653.（移动内管 Hele-Shaw 框架）｜出自：归档[S18]。
76. Tardy (2018). JPSE 162.（含套管运动环空顶替模型续篇）｜出自：归档[S17]。
77. Bu et al. (2018). JNGSE 52:317–324.（高斜井套管旋转提效）｜出自：归档[S19]。
78. Dai & Liu (Pegasus Vertex) (2018). AADE-18-FTCE-094.（standoff×旋转/往复/排量矩阵；**70% standoff + 1 rpm 即大幅清通道，但存在最优转速**）https://www.aade.org/download_file/1335/394 ｜出自：独立文献调研 §2.3、归档[S1]。

### B8. 其余值得主编留意的条目
79. Frigaard 组 (2018). Nordic Rheology Society 26.（不规则井眼 τy/几何效应、静态通道形成条件）https://nrs.blob.core.windows.net/pdfs/nrspdf-2f7ee2dc-2158-41d5-9972-a5a9b0a602de.pdf （注：此 URL 在下载清单 E#8 亦被标为 Roustaei 2015 开放直链——⚠️同一 URL 两处归属，主编核对）｜出自：归档[S3]。
80. Foroushan, Lund, Ytrehus & Saasen (2021). Cement Placement: An Overview of Fluid Displacement Techniques and Modelling. Energies 14(3):573. https://www.mdpi.com/1996-1073/14/3/573 （OA 已全文提取；**确认 2DGA 是工业标准的综述**，M3 Lockyear 判据 Eq.2 的转引源）｜出自：环空失真调研 §参考文献、下载清单 E#2。
81. Malekzadeh & Heydari（出处待核，疑为 Foroushan et al. 2020 误写）｜出自：归档综合修复调研 Gaps#4。

**收割小结**：本次共收割 **81 条**去重条目（其中约 8 条带"⚠️出处冲突/待核"标记，其余为各报告已核实或带 URL/DOI/SPE 号佐证）；与"偏心环空顶替/居中度敏感性/窄边窜槽/屈服静态残留/Hele-Shaw/数值方法"六主题直接相关的核心 40 条集中在 B1–B4。另有 SWPU 汇总中前置液/冲洗液/控压固井设计类条目（★9/★10/#12–19/#27–28/#30–33 等）与两份事故调研报告文献未纳入本主题摘要，主编需要时可直接回源文件。

---

## C. 未决修复项汇总（与居中度/环空闭合/弥散/指标口径相关）

> 标注：【已修】=B1/B2/B3/M 系列已落地；【未修】=各报告立案至今无落地证据；【口径规避】=代码不动、用报告口径绕开；【待重验】=修复后需重新评估权重。

### C1. 环空闭合 / 横流物理（当前最大物理未决项）
1. 【未修·P2-1/M5】**真实 Hele-Shaw 方位横流**：现 `v` 场是连续性回推非真实横流（初诊报告 §5.2/§8 P2-1）；环空失真调研 M5 建议解 Zhang22 式(2.2–2.3) 椭圆流函数实现宽边→窄边浮力再充填（3–5 天，归档综合修复调研 P2-1 重申）。这是窄边"救回"通道，直接关系居中度响应形状。
2. 【未修·M5 子项】**I3 浮力救援通道弱**：eta2/delta_rho 已局部化（enable_local_i3，M 系列），但 I3 通量仍 ∝Δρ·H³（窄边 H 小通量小）+ clip ±0.5 限幅；R2 消融曾 +2e-10；I3 公式缺 g、Froude F=1 硬归一、浮力到达 factor=1.0 未反标定（09-02 独立调研 R4/P2-1，【未修】标定链不闭合）。
3. 【未修·P2-3】**窄缝穿透二值门槛**（窄缝宽度 vs 水泥颗粒/屈服应力"能进则进、不能进则完全不进"，初诊报告 P2-3）——屈服门槛已由 B2 物理屈服门部分承载，但颗粒级穿透判据无对应物。
4. 【未修·P2 扩展】**套管旋转/往复建模**（归档综合修复调研 P2-2：Tardy15/18、Carrasco-Teja09/10、Bu18；SINTEF/杨鑫25 提供定量锚点）——前提是现场有活动套管作业。
5. 【未修·P2-3'】**standoff 沿深非均匀分布**：扶正器间套管垂曲计算+7 井 standoff 实测/估算（归档 P2-3）；**7 井 standoff 仍全部 model_assumption**（08-20 报告层 2.1、归档 R6），仅 hu101 实测剖面已接入 loader。
6. 【未修·M6】**流变井下温压修正**：全部流变为井口单点常值（65/93℃ 化验）；YP 温度减半即把 Lockyear 门槛 standoff 从 0.72 拉到 0.12–0.36——M3 门槛的井下定量结论因此只能做区间敏感性（环空失真修正总结 §6.1、归档 P2-4）。

### C2. 弥散 / 数值方法
7. 【已修但留尾巴】**弥散物理化**：M1 dt 归一已实施，但 B2 证明弥散开关对结果几乎无影响（非承重项）；终局报告 §6.4 把"弥散物理化（D·dt/dz² 有量纲形式）"列为可选增强=纯审稿健壮性，不阻塞论文。历史教训：弥散系数 0.018/0.015 无文献出处且曾主导结果（scale 0→1 时 ht1_004 效率 1.7%→66.5%，归档 R4）。
8. 【未修·T2-3】**TVD/minmod 激波捕捉替换拉普拉斯平滑**（Tier 2 待做；Z22 用 FCT、B25 用 TVD+minmod；Boniou 2022 印证常数平滑=过弥散）。
9. 【未修·M2 标定】**re_crit(He)=2100(1+0.1He) 是临时标定钮**（代码自注未标定；环空失真修正总结 §6.2）；需 Walton&Bittleston 解析解+Z22 Table3+Foolad 定性三重回归；当前 8 井因 He 高（24–76）全层流、M2 中性，高排量/低黏井启用前必须标定。
10. 【未修·T2-1/T2-2】**HB 闭包查表 + augmented Lagrangian 间隙求解器**（Bararpour25 路线，各 2 周/2–3 周）——真实 plug region 与 static wall layer 的严格解法，当前用 Papanastasiou M=100 正则化软化屈服。
11. 【未修·T2-4】**瞬态流函数演化**（脉动效应，Moyers-González07）——准静态本质缺陷，仅当需建模脉动固井时做。
12. 【待重验】**b³ 分流高偏放大的权重**：A4 证伪报告定位的 529.7 倍放大机制，在 B1/B2 修复后口径下"需重验"（原文状态声明明示）；M2 已实现但对当前 8 井中性。
13. 【未修·小项】e 下限 0.05 无注释依据（完全居中井被强加 5% 偏心，09-02 独立调研 R4#3）。

### C3. 指标口径 / 论文叙事
14. 【硬要求】**验证叙事换轴**：修复后 8 井 η_E 挤在 0.59–0.999 无区分度（ht1_004 η_E=0.9985 ↔ CBL 0.3%），必须执行两轴叙事=顶替充分性（η_E/η_N）+候凝风险（失稳/窜槽），并并列报告 η_E（库存）/入库完成度/η_N 三口径（终局报告 §6.1–6.2，"硬要求，不是可选项"）。
15. 【硬要求】**η_E≡水泥库存比的定义声明必须写进论文**，否则"效率 0.4"被审稿人读成"60% 没顶干净"；禁止写模型-现场相关系数（无落盘依据、CBL 口径逐井不统一）（09-02 两报告+终局报告一致）。
16. 【未修】**评价窗效率 vs CBL 窗口对照的重测与图注**：hu103 139.7mm 段 CBL 12.06% 与模型域对应关系需图注写清（终局报告遗留⑤）；CBL 数字化窗口/官方判据双口径矛盾仍未裁定。
17. 【未修·M0 后续】**standoff ±0.1 区间标注**（长江大学 11% 误差依据）：08-20 报告 P1-1"每井给 standoff ±0.1 → η_N 区间"至今未落地（审稿必问）。
18. 【未修·低优先】失稳指数仍为启发式（1−exp(−...)），未替换为 Zhang23 三分类判据/Muskat regime（Tier 0 判据模块已建 flow_classification/muskat_regime 但作为诊断存在，summary 失稳指数本体未换）；失稳指数去饱和（线性/对数列）已由 M0 做。

### C4. 质量守恒 / 域口径（B1 修复后的残留）
19. 【口径规避】**2D 场体积记账结构性近似**：半拉格朗日在索引空间重采样不带体积权重 b、顶部吸收式饱和、五相过填归一销毁——"域内水泥 m³"与"入环空 m³"不可直接做质量平衡论述；论文改用"入库完成度+域满饱和度"两口径（终局报告 §2 守恒记账说明、遗留③）。原 P1-1（过填归一化守恒化/销毁记账回补）与 P1-2（顶部真出流 BC 替代反射钳位）未以代码落地，降级为可选增强（终局 §6.4）。
20. 【未修·小项】重合段井径口径偏差 13–28%（hu103 技套内径为 proxy 245.42；hu103/ht1_001/hu101 改技套 ID×尾管物理外径未落地）（09-01 报告 §3.4/遗留）。
21. 【未修·资料级】hu103 管容 88.55（实测）vs 20313 分段理论 90.88 差 2.3 m³ 无逐字构成记录，loader notes 双口径并列；<1 m³ 级精度需调 0708 原件逐段复核（终局遗留②）。
22. 【配置敏感】hu101 窄边几乎未推进（η_N=0.03）的量级对 `yield_gate_f_safety=1.15` 敏感——单参数敏感性分析列为可选增强（终局遗留④/§6.4）。
23. 【资产】历史 GIF/PNG 全部是 bug 时代产物，论文配图一律用终跑轮重新生成（终跑轮已含 hu101/hu103/呼102 新 GIF）（终局遗留⑥）。

### C5. 明确"不做/已排除"清单（防止主编重复立项）
- ❌ 解除 e 截断当"居中度死区解药"（A4 证伪，反物理）；e 上限已护栏式参数化 0.90。
- ❌ b_narrow≥5mm 物理下限（hu101 实测窄边 8.4mm 本就大于 5mm，下限只会掩盖真实风险）。
- ❌ "湍流救窄边"作为普适结论（Foolad21/Couturier/M-F19 三方反证）。
- ❌ 全局降幂律分流指数（幂律层流指数 2+1/n≈3.39>3，两侧都层流时 b³ 反而偏缓）。
- ❌ 刮泥器定量建模（无文献支持）。
- ❌ 把任何 8-29 之前的效率表"修回来"或进论文正文（表4/08-16/08-24/08-29/09-02 各轮全部作废）。
- ❌ U 型管建模（已删）；ML/PINN 数据驱动作论文主卖点（SWPU 已抢先）。

---

## 给主编的三句话导读

1. **"居中度体现不出来"在当前口径下已不成立**：B1 截面守恒 bug（半速前缘 4 个月）才是历史弱响应主因；修复后 8 井 η_E 0.586→0.9985 随居中度单调，hu101 低值=真实窄边窜槽。历史"三根因（e 截断/指标稀释/弥散钉死）"叙事中仅指标稀释层仍有效，且已升级为"阈值型响应+三件套指标"方案（A3 §5）。
2. **1D+2D 模型本身在框架层被两轮独立调研判定为合理**（B1 组文献），遗留的是：横流/I3 救援通道（M5）、流变温压修正与 7 井 standoff 假设（M6/数据层）、HB 真屈服（Tier2）、I3 分母 6/m 待核——这些是"改进空间"而非"模型错误"。
3. **数值计算问题已收敛到三条**：弥散（B2 证明非承重，留审稿健壮性收尾）、过填归一化/顶部吸收（口径规避，未代码化）、re_crit/标定链临时钮（高排量井前必标定）。论文侧的硬约束是 C3 两条叙事要求。
