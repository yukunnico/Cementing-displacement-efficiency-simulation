# 模型求解流程与流变参数取证报告（agent2，2026-09-04）

> 范围：`cemdisp` 包（HEAD=4eb3d85"9_4修改1"——该提交仅落盘 .tmp_research/取证产物与 results 重跑，**代码与 396dc77 逐位一致**，下述行号对两 HEAD 均有效）。只读取证，未改任何代码。
> 方法：逐文件精读 annulus_d2dga.py(1404行)/casing_flow.py(1258行)/boundary_bridge/d2dga_flux/regime_closure/shoe_timeline/pipe_exit_state/interface_tracking/fluid_spec/pumping_schedule/well_spec 全文，8 井 loader 流变段精读，diagnostics 关键函数精读，并用 conda env 实例化 8 井 fluids 实证本构分布、表观黏度触限与 1D/2D 口径差。
> 前次锚点核对：任务清单所列锚点（:784/:778-779/:862-863/:526-560/:240/:612/:628-634/:257-258/:802-811/:1233-1236/:253/:878-906 等）**全部核实成立**，行号有 ±3 行内的漂移已在下文按当前代码标注。

---

## 1. 端到端求解流程（文字版阶段图）

```
【阶段 A · 数据装载】loader（cemdisp/data/loaders/<well>_loader.py）
  ├─ CSV 井径/井斜实测剖面（hu102: 62点现场实测+顶段等效点，hu102_loader.py:160-178,246-251）
  ├─ WellSpec：井段/鞋深/管容链 shoe_lag_volume_m3/居中度分段代理(:181-203)/评价窗(:268-278)
  ├─ fluids: 元组[FluidSpec]——每流体=角色+密度+本构(四选一)+参数（hu102_loader.py:297-380）
  └─ PumpingSchedule：泵注步骤序列（体积/排量/event_tag），末步可带 RESTART 标签(:443-450)
                 │
【阶段 B · 1D 套管内输运】CasingFlowSolver.run()（casing_flow.py:225-354）
  B1 管容：_pipe_cross_section_area(:788-817) 优先级 shoe_lag_m3 > pipe_id_profile > liner_id；
      到达鞋口迟到体积 _timeline_pipe_volume(:771-786)
  B2 胶塞截断 _displacement_sequence_cutoff(:843-889)：顶替序列=首个 RESTART 步前全部步骤；
      pumping_end=首个 RESTART 的 start_time（无 RESTART 井=全日程最晚 end）；
      RESTART 及其后体积不入累计 → 替浆前缘永不到达 → "替浆不进环空"由构造保证
  B3 界面推进（纯体积坐标法）：_front_arrival_time(:902-917)=该步累计起点+管容 在日程中
      插值出时刻；_rear_arrival_time(:919-934)=累计终点+管容。不可达→None
  B4 重力修正 _gravity_corrected_arrival_time(:944-1016)：Atwood 数 at=|Δρ|/(ρ1+ρ2)(:991)
      ×井斜投影 cos(β̄)(:998-1000)×屈服抑制(1−0.8·min(τy/τcrit,1))(:1004-1010)；
      重驱轻加速/轻驱重减速(:1012-1016)。旧经验乘子路径 :1018-1063（enable_buoyancy_physics=False）
  B5 stop ≡ cement_end_time_s(:287-318)：最后一段水泥的尾缘过鞋口时刻
      （尾缘≡后继流体前缘，重力修正按同一配对 :300-310）；只上抬不下压约束 :316-318；
      胶塞截断后 = min(尾浆尾缘过鞋口, 顶替序列终点)
  B6 鞋口时间线 _build_shoe_timeline(:672-731)：切换/FRONT_ARRIVAL/REAR_EXIT/END 事件
      → 轴向弥散 _apply_dispersion_to_timeline(:562-670)：Taylor-Aris 系列弥散系数
      _compute_dispersion_coefficient(:394-471) + 混浆增强因子 _interface_instability_factor
      (:508-560, At>0 且 Re>100 → 1+5·At·√(Re/100)，上限10) → 每前缘生成 5 个 erf 过渡子事件(:637-652)
      + frac=1.0 收尾事件(:659-668, F4 修复)
                 │
【阶段 C · 1D→2D 耦合】boundary_bridge.build_coupled_annulus_inlet_provider(:139-206)
  生产链用 _legacy_provider(:179-206)：
    泵注期(Q>1e-9)：shoe_timeline.at(t)=最近过去事件 → 多相过渡带相分数
    停泵期：回退 pipe_exit_state_at(:356-392)→停泵沉降增强 _settled_exit_fluid_name_enhanced
      (:1136-1202, 凝胶指数发展 :1184+屈服抑制 :1188-1193)
  相映射 _phase_fractions_for_fluid(:72-97)：split_cement_phases=True 时
    LEAD/INTERMEDIATE→"lead"、TAIL→"tail"、其余水泥→"cement"、WASH/SPACER→"spacer"、
    FLUSHER→"flusher"、其余（含替浆）→"mud"
  stop/碰压：runner/scripts 的 build_case 直接取 cement_end_time_s 作 total_t
    （scripts/_mass_balance_diag_20260902.py:60；runners/hu102_tailpipe.py:148-183 的
    annulus_stop_time_s 仅作 cement_end 缺失回退）。2D 求解窗=[0, stop]
                 │
【阶段 D · 2D 环空 D2DGA】AnnulusD2DGASolver.run()（annulus_d2dga.py:938-1404）
  D1 几何 _build_geom(:362-438)：s∈[0,域长] 140→250 格、y∈[0,π·r̄] 40 格；
      e=clip(1−standoff, 0.05, 0.55)(:404)；b=2·h̄(1+e·cos(πφ))(:412-416)；
      体积校正 scale→物理环空体积(:431-437)；域顶 open_outlet=True 默认(:246)
  D2 _pick_fluids(:451-467)：mud/lead/tail/spacer/flusher 五槽位（第一个 MUD 角色当选）
  D3 时间循环 while True(:998-1289)，enable_cfl_adaptive=True 默认：
     每步先问边界 provider(t)（:1006）→ pump_active = Q>1e-9(:1016)
     ├─ 泵注分支(:1018-1188)：
     │   a. _compute_velocity(:698-876)：γ̇=6|w_prev|/b(:612,755) → 四本构表观黏度
     │      _apparent_viscosity(:476-509) → 相体积加权 μ(:613-619)+ρ(:620-626)
     │      +τy 体积加权(:628-634)+auto-m 场(:637-648)+η1/η2 相黏度(:649-651)
     │      +n_mix/kappa_mix(:652-660) → Papanastasiou 正则化 M=100(:756-765, mu_reg)
     │      → 两层黏度闭包 η_mix=1/(c̄³/η₂+(1−c̄³)/η₁)(:777-779)
     │      → base=(b/b̄)²/η_mix×I1(c̄,m)(:780-785) → 浮力方位修正(式4.24,
     │      correction=clip(Δρ·I2/I1,±0.5)·f_φ,:800-811) → pref(:817)×(1−wall)(:819-820)
     │      → w=q_half·pref/Σ(pref·b·dy)(:862-863)（2026-09-02 守恒修正，全环空 2∫wb dy=Q）
     │      → 横向 v 由连续性(:865-873)
     │   b. CFL 自适应 dt(:1037-1040→_compute_cfl_dt_step:878-906, cfl=0.5, dt_min=0.1)
     │   c. D2DGA 通量放大 f(c̄,m)(:1043-1050)→半拉格朗日平流 w·f(:1052-1063)
     │   d. 五相过填修正(:1065-1070)
     │   e. 显式拉普拉斯弥散(:1079-1098)：轴向0.018/方位0.015（lead/tail），
     │      0.012/0.012（spacer/flusher），按 dt_step/dt_ref 归一(M1)
     │   f. I3 浮力弥散通量(:1101-1142)：q=(ΔρH³/6η₂)·I3(c̄,m)·(−f_ξ,f_φ)(d2dga_flux.py:120-142)
     │      默认全场均值口径(enable_local_i3=False, :1112-1114)，散度→水泥相分配(:1130-1142)
     │   g. 累计入环空体积上限(:1146-1153, open_outlet=True 时不限)
     │   h. 壁面静止层 wall(:1159-1188)：M3 屈服门 _yield_gate_wall(:526-560) 默认启用
     │      （τw=G·b/2 外推≤1.15·τy 冻结，参考元永不冻结；整列无流动冻结兜底）
     ├─ 泵停分支(:1190-1211)：浓度场冻结（凝胶强度足以抗浮力滑塌），只算 dt 推进
     └─ 指标记录(:1214-1271)：bulk_fill=∫b·cement/∫b(:1235)、eff≡cement(:1233-1234)、
        前缘/窜槽/混浆/失稳指数(:1238-1250, mobility=b³/μ)、mean_wall
  D4 汇总(:1291-1404)：metrics DataFrame、depth_profiles(:908-936)、summary
      （η_E=末步 effective_efficiency :1322,1348；η_N 窄四分位 :1349；评价窗效率 :1358→
      _evaluation_window_efficiencies:75-92 即 b 加权浓度均值；低尾指标 :1359；
      浮力数 b=(Δρ)g d²/(μ_mud·ŵ₀) :1334-1340→_compute_buoyancy_number:686-696）
  D5 Tier0 诊断聚合(:1395-1402, try/except 不拖垮主求解)
                 │
【阶段 E · 输出】runners/<well>_tailpipe.py + scripts/rerun8_*.py
  时间序列 CSV/深度剖面 CSV/摘要 json+md/PNG 图族/NPZ(2D场快照+网格)/GIF 动画
  （runners/hu102_tailpipe.py:76-141；终跑脚本 :199-206 单井产物子集）
```

**生产链一锤定音**：三连终跑（守恒修复 0902→管容链修复 0902→胶塞语义修复 0903，HEAD 已落盘最后一轮）全部用 `AnnulusD2DGASolver(total_t=stop, nz=250, ny=40)` 纯默认参数（scripts/rerun8_after_plug_semantics_fix_20260903.py:148、rerun8_after_pipe_capacity_fix_20260902.py:138、_rerun8_after_conservation_fix_20260902.py:28），build_case 在 scripts/_mass_balance_diag_20260902.py:52-64。即生产口径=：**e_clip 0.55、屈服门 B2 默认开(f=1.15)、regime_split=False、弥散 0.018/0.015、CFL 自适应开、open_outlet=True、enable_local_i3=False、M=100、c_min=0.05、wall_seed=0.005**。

---

## 2. 模块功能清单

### 2.1 transport1d/（套管内 1D）
| 文件 | 功能 | 关键点 |
|---|---|---|
| casing_flow.py | 1D 体积推进求解器（生产主链） | 见 §1 阶段 B；弥散系数 :394-471（牛顿 Taylor-Aris :435 / 幂律 Batot2016 式28 :437-441 / Bingham Fan&Wang1966 k(ξ₀) 多项式 :443-468 / HB→等效 Bingham :446-450）；对流尺度上限 d_cap=α·U·R 截断 :431-432（高 Pe 下 TA 渐近不可用，hu101 3% 崩塌教训 :405-411） |
| interface_tracking.py | 数据类 InterfaceFront(fluid,depth,time)（:18-28） | 纯结构 |
| shoe_timeline.py | ShoeEventKind 六事件(:27-38)+ShoeTimeline.at()="最近过去事件"查询(:111-157) | 纯结构+查询 |
| pipe_exit_state.py | 数据类 PipeExitState(time,Q,stage,phase_fractions)（:23-34） | 纯结构 |
| boundary_bridge.py（物理上属 1D→2D 桥，文件在 models2d/） | AnnulusInletState(:47-58)、相映射(:72-97)、provider 两形态(:139-206) | 生产=_legacy_provider(:179-206) |

### 2.2 models2d/（环空 2D）
| 文件 | 角色 | 生产/备用 |
|---|---|---|
| annulus_d2dga.py | **生产主求解器**（几何/速度/平流/弥散/I3/wall/指标/summary 全在此） | 生产 |
| d2dga_flux.py | D2DGA 纯函数库：f(c̄,m)= [m c̄²+1.5(1−c̄²)]/[m c̄³+(1−c̄³)](:45-83)、I3(:86-117)、q_buoy(:120-142)、I1=Bararpour 式2.24(:145-155)、I2=式2.25(:158-168) | **生产**（被 annulus_d2dga 直接调用 :35-40） |
| regime_closure.py | M2 流态修正纯函数：Metzner-Reed Re_p(:9-15)、Hedstrom He(:18-26)、层流 f=24/Re_p(:29-37)、Dodge-Metzner 湍流迭代(:40-53)、过渡插值(:56-64)、drag_weight R(:67-84) | **默认休眠**——仅 enable_regime_split=True 时被 annulus_d2dga.py:830 动态 import 消费；默认 False → 生产 8 井 Re/He/f/R 一个都没算 |
| boundary_bridge.py | 1D→2D 桥（见上） | 生产 |

**enable_regime_split=False 的含义**：M2（Maleki & Frigaard 2017 式58-66 的层流/过渡/湍流阻力权重 R 修正流量分配）整条链只存在于 `_compute_velocity` 的 if 分支(:825-860)；默认 False 走 else 分支(:861-863) 纯 b²/η_mix 流动度分配——即生产 8 井是**全层流假设**，re_crit=2100(1+0.1·He) 的 provisional 公式(:833-836) 也从未在生产中执行。

### 2.3 diagnostics/（纯后处理，不反馈求解）
| 文件 | 功能 |
|---|---|
| tier0_diagnostics.py | 聚合入口 compute_all_tier0_diagnostics(:130-214)：T0-1 flow_classification/T0-2 muskat_regime/T0-4-5-7 displacement_metrics/T0-3 buoyancy_regime/T0-6 shutdown_decay，各子项独立 try/except |
| flow_classification.py | 流动分类（非弥散稳态/弥散稳态/非稳态弥散，:326-330），ŝ0 由水泥场等值面估计 |
| muskat_regime.py | Muskat 三 regime 判别(:345+)；黏度输入 _effective_viscosity_pa_s(:470-479)=优先 PV、否则 K（**幂律流体的"牛顿等效黏度"直接用 K 数值，无 γ̇ 折算**——诊断口径粗化） |
| displacement_metrics.py | 泥浆滞留率/界面长度比/窄四分位效率 _narrow_quarter_efficiency(:211-226)/突破时间/ŵ₀ 估计 |
| regime_classifiers.py | 浮力数分类(:125-172)：b<0 forbidden、<20 highly_dispersive、<80 steady_capable、≥80 non_dispersive；停泵衰减 T0-6(Moyers-González 2007 式3.35/3.40, :205+)，τy 取各流体最小屈服应力 |

### 2.4 data 层
| 文件 | 功能 |
|---|---|
| fluid_spec.py | RheologyModel 四枚举(:39-45)+FluidRole 九枚举(:48-59)+FluidSpec 冻结数据类与按本构必填校验(:67-109) |
| well_spec.py | WellSpec(几何剖面/管容链/双径字段/评价窗/standoff_profile)（:97-137） |
| pumping_schedule.py | PumpingScheduleStep(step_name/fluid/volume/rate/start_time_s/event_tag)+Schedule（:45-136） |
| loaders/<well>_loader.py | 8 井现场值装载（流变常量→FluidSpec→schedule→ValidationData） |
| provenance.py / fluid_provenance.py | 数据源画像/注入流体符合性摘要（仅报告用） |
| validation_data.py | CBL 路径/合格率句柄（仅报告用，validation/ 目录已于 2026-07-29 删除） |

---

## 3. 流变参数专项（核心）

### 3a. 定义层：每流体本构与参数（8 井实证，conda 实例化输出）

各井流变常量定义位置：hu1_loader.py:127-156、hu101_loader.py:126-145（含主检/复检双套 :142-145,295-298）、hu102_loader.py:120-157、hu103_loader.py:111-140+、ht1_001_loader.py:178-193、ht1_003_loader.py:252-280、ht1_004_loader.py:272-315、hu2_loader.py（与 ht1_001 同体系）。

| 井 | 钻井液 | 隔离液/前置液 | 领浆 | 尾浆 | 替浆 | 备注 |
|---|---|---|---|---|---|---|
| hu1 | Bingham 0.076/8 | **幂律** n=0.649/K=0.585 | PL 0.732/0.933 | PL 0.666/0.906 | Bingham 0.066/10 | 204121/204131 实测 |
| hu2 | Bingham 0.058/8.5 | PL 0.545/1.338 | PL 0.811/0.876 | PL 0.886/0.453 | 0.058/8.5 | |
| hu101 | Bingham 0.058/9.2 | Bingham 0.030/5(YP代理) | PL 0.844/0.381 | PL 0.830/0.352 | =泥浆 | 主检89℃口径 |
| hu102 | Bingham 0.066/10 | Bingham 0.030/8(YP无实测) | PL 0.737/0.947 | PL 0.737/0.947 | =泥浆 | 20211/20234 实测 |
| hu103 | Bingham 0.054/10 | **PL 0.54/2.12**（三段同） | PL 0.82/0.67 | PL 0.76/1.14 | =泥浆 | 203111 表7 实测；平衡液 role=MUD(0.025/1.5) |
| ht1_001 | Bingham 0.051/6 | PL 0.545/1.338 | PL 0.811/0.876 | PL 0.886/0.453 | =泥浆 | |
| ht1_003 | **幂律** n=0.631/K=0.751 | PL 0.668/1.245 | PL 0.597/1.622 | PL 0.585/1.673 | Bingham 0.062/10 | 全 8 井唯一非 Bingham 泥浆→**τy=0** |
| ht1_004 | Bingham 0.053/8.5 | Bingham 0.058/9.8、0.065/10 | **Bingham** 0.17/13 | **Bingham** 0.18/14 | 0.050/9.5 | 化验幂律(0.853/0.746,0.869/0.669)仅在注释，[optimized]=旧 MATLAB 反演值 |

要点：
- **三种本构混用**（Bingham/幂律为主，HB 枚举存在但 8 井未用）；本构选择是逐流体逐井的，同一流体角色跨井可不同本构。
- 全部从 loader 常量进 FluidSpec——数据源=0708 现场提取包（各常量行内注 field_measured/化验报告/20211.doc 等），其中部分字段显式标注 proxy/无实测（hu101 隔离液 YP、ht1_004 全套为 optimized 反演）。温度口径：实测值注明测温（65℃/89℃/93℃/140↘93℃），但**模型本身无温变**（见 3d）。

### 3b. 1D 套管侧：流变参数进哪些计算
| 计算 | 位置 | 消费的流变 |
|---|---|---|
| 轴向弥散系数（Taylor-Aris 系列） | casing_flow.py:394-471 | 牛顿：仅 PV；幂律：n（Batot 式28 :440）；Bingham：τy+PV（Fan&Wang ξ₀=τy R/(4 μ_app U) :454）；HB：τy+K+n（:446-450）。**弥散宽度 σ_t 直接受本构影响** |
| 混浆增强因子（界面失稳） | :508-560 | 有效黏度 _effective_viscosity(:473-506)：通用式 τy/γ̇+K·γ̇^(n-1)，γ̇=8U/2R。**注意 Bingham 的 K=None→fallback 0.01(:503)，剪切项用 0.01 而非 PV**；At 用密度(:539-541)、Re 用几何平均黏度(:546-554) |
| 重力到达时间修正 | :944-1016 | 仅密度（At）+**τy**（屈服抑制 :1004-1010）；PV/n/K 不参与 |
| 停泵沉降增强 | :1136-1202 | 密度差+τy（:1188-1193）+凝胶经验时间常数（非流变参数） |
| 1D 无摩阻/Re/泵压计算 | — | 1D 是纯体积推进，无动量方程；流变只通过"弥散+重力修正"间接影响鞋口边界 |

### 3c. 2D 环空侧
- **η(γ̇) 曲线**：`_apparent_viscosity`（annulus_d2dga.py:476-509）支持全部四种本构，γ̇=6|w_prev|/b（:612,755，w 一步滞后）。**裁剪 μ∈[1e-5, 3.0] Pa·s(:509)**——实证：hu102 钻井液(Y=10Pa)在 γ̇≤2 s⁻¹ 即触 3.0 上限（窄边死区实际发生）；hu103 隔离液 γ̇<1.3 时触限。
- **η_mix 两层闭包**(:777-779)：η₁=泥浆相表观黏度场、η₂=水泥相表观黏度场（lead 优先，无 lead 取 tail，:637-651）。**隔离液/冲洗液黏度不进闭包**——spacer 占据区 c̄≈0 → 按 η₁（泥浆黏度）处理；spacer 的 PV/YP/K 只进入 (i) 线性混合 μ(:613-619)→该 μ 只用于 Re 诊断(:769)、I3 的 η₂ 场均值口径(:1113)、mobility 指标(:1247-1248)；(ii) τy 混合(:628-634)；(iii) 密度混合→浮力修正与 I3 的 Δρ(:620-626,807,1111-1114)。
- **屈服门的 τy**：`_yield_gate_wall`(:526-560) 用 **τy 体积加权混合场**（:628-634，含 mud/lead/tail/spacer 各自屈服；幂律/牛顿流体 τy≡0，:469-474），判据 τw_extrap=G·b/2≤1.15·τy。不是"取被顶液"也不是"取顶替液"，是当地混合物。同一 τy 场还参与 Papanastasiou 正则化(:756-765)。
- **浮力数 b**(:686-696, 调用 :1334-1340)：μ_mud 用 `mud_fluid.plastic_viscosity_pa_s or 0.05`（**用 PV 标量，不含 τy**；幂律泥浆(仅 ht1_003)时 PV=None→fallback 0.05）；Δρ 用密度。失稳指数 mobility=b³/μ（:1247-1248）用正则化混合黏度 μ_reg。
- **auto-m 黏度比场**(:637-648)：m=μ_mud(γ̇)/μ_cement(γ̇)，clip[0.1,10]，供 f(c̄,m)/I1/I2/I3。

### 3d. 死参数排查（读入但未进求解器/未进生产路径）
| 参数 | 状态 | 证据 |
|---|---|---|
| **n/K 的 M2 消费（Re_p/He/流态权重）** | **默认死** | n_mix/kappa_mix 在 :652-660 算出并随 `_compute_velocity` 返回，但 run() :1020 接收为 `_n_mix/_kappa_mix` 直接丢弃；唯一消费者是 regime_closure 链(:832,840)，被 enable_regime_split=False(:241) 关闭。**生产 8 井中 n/K 只活在 _apparent_viscosity 的表观黏度里（这已经足够"活"）** |
| enable_local_i3 的 η₂ 透传 | 默认死 | :1109-1114，False 时 η₂ 用全场 mean(μ)、Δρ 用全场均值 |
| e_clip_max=0.90 | **非生产** | docstring :312 称"生产跑道显式设 0.90"，但三连终跑脚本均纯默认 → 生产实际 **0.55**（standoff>0.45 一律压到 0.45 偏心）。仅 09-02 机制核查脚本(bisect_hu103/c_verify_convergence 等)用 0.90 |
| has_plug 混浆增强=1 语义 | 生产未启用 | 8 井无任何 runner/scripts 设 has_plug=True（前次取证复核仍成立） |
| f_amp 的 [0.5,2] clip | **注释与实现不符** | :647 注释称"d2dga_flux 内还有 [0.5,2] clip"，实际 d2dga_flux.py:29-30 默认 min/max=±inf，调用处 :1048 未传 → **无裁剪** |
| mu_turbulent | 占位零场 | :770, 只为兼容旧结果对象 |
| d2dga_viscosity_ratio=1.0 | 被绕开 | enable_d2dga_auto_m=True 默认走 m 场(:1044-1047) |
| 温度/压力修正 | **不存在** | grep 全 models2d/transport1d：仅 docstring 声明"泥饼、温度、凝胶强度、湍流修正不再影响求解"（annulus_d2dga.py:17,1387；d2dga_flux.py:58）。流变参数是**常温常数**，井温 6500m+ 无任何折算 |
| spacer/flusher 黏度进 mobility | **结构性死区** | 两层闭包(:777-779)无第三相槽位，隔离液黏度对流动度 base 无直接贡献（仅经 μ 的间接诊断/τy/I3-均值口径） |

### 3e. 1D 与 2D 流变口径一致性
**同一数据源、两套实现、两处口径差**：
- 数据源同一：同一 `fluids` 元组既进 CasingFlowSolver.run() 又进 AnnulusD2DGASolver.run()，无第二套参数表。
- 实现各写各的：2D `_apparent_viscosity`（四本构，clip[1e-5,3.0]，γ̇=6|w|/b）vs 1D `_effective_viscosity`（τy/γ̇+K·γ̇^(n-1)，γ̇=8U/2R，无 clip）。
- **实证口径差**（hu102 钻井液，γ̇=80 s⁻¹）：2D μ=PV+τy/γ̇=0.191 Pa·s；1D=τy/γ̇+0.01=0.135 Pa·s（**1D 对 Bingham 用 K fallback 0.01 顶替 PV，差 41%**）。该黏度只进 1D 混浆增强判据（影响鞋口过渡带宽度→环空入口相分数），不影响 2D 场。
- 1D Taylor 弥散的 Bingham 分支用 μ_app=PV（:452，正确口径），与 _effective_viscosity 不一致——同文件内两处 Bingham 表观黏度公式不同。

---

## 4. 三行裁定

**① 流变参数纳入程度总评：部分纳入，但"活"的部分是主链**。PV/τy（Bingham 8 井主力）与 n/K（6 井水泥浆+4 井隔离液）**全部经 `_apparent_viscosity` 进入 2D 主求解**（表观黏度→μ/τy 混合→η_mix 闭包→mobility→速度场→平流）与 1D 弥散/重力修正——**没有"读进来只进报告"的流变参数**。但有三层折扣：(a) n/K 的流态学消费（Metzner-Reed/He/阻力权重）被 enable_regime_split=False 整体关闭，生产是全层流 b² 假设；(b) 隔离液黏度被两层闭包结构性排除在流动度之外；(c) ht1_003 泥浆是幂律 → τy=0，该井屈服门/屈服混合对泥浆相失效（其 n/K 活）。温度/压力修正不存在，所有流变是地面实测常温常数外推到 6823-7735m 井段。数值口径上 μ 上限 clip 3.0 Pa·s 在低剪切死区实际生效（屈服流体 γ̇≲2 s⁻¹ 全部压平）。

**② 求解器物理闭合与 D2DGA 谱系对应表**：
| 要素 | 本代码 | 谱系出处（代码内注） |
|---|---|---|
| 间隙几何 | b=2h̄(1+e cosπφ)，e=clip(1−standoff,0.05,0.55)（:404,412-416） | Bittleston/Walton 偏心槽 |
| 流动度 | base=(b/b̄)²/η_mix×I1(c̄,m)（:780-785） | Hele-Shaw；I1=Bararpour 2025 式2.24 |
| 两层黏度闭包 | η_mix=1/(c̄³/η₂+(1−c̄³)/η₁)（:777-779） | Zhang & Frigaard 2022 式4.23 |
| 通量放大 | f(c̄,m)=[m c̄²+1.5(1−c̄²)]/[m c̄³+(1−c̄³)]（d2dga_flux.py:45-83） | Z22（m 场化=R1 auto-m :637-648） |
| 浮力方位修正 | pref×(1+clip(Δρ·I2/I1,±0.5)·f_φ)（:800-811） | 式2.5b/4.24；f_φ=r_a sin(πφ)sinβ(:663-684) |
| 浮力弥散通量 I3 | q=(ΔρH³/6η₂)·I3(c̄,m)·(−f_ξ,f_φ)，散度进水泥相（:1101-1142;d2dga_flux.py:120-142） | Z22 式4.25 第二项/I3=式4.26 |
| 屈服正则化 | Papanastasiou M=100 替换屈服项（:756-765） | Papanastasiou 1987 |
| 壁面静止层 | M3 可逆屈服门 τw≤1.15·τy（默认，:526-560）或 T1-5 c_min=0.05 浓度门（Bararpour 式2.35-2.41，OFF 兜底 :1172-1188） | Bararpour 2025 |
| 流量分配 | w=q_half·pref/Σ(pref·b·dy)，全环空守恒（:862-863，09-02 修复后 ∫wb dy=Q） | 幂律近似 b³/η 分流 |
| 流态修正 | M2 固定点迭代（regime_closure，Maleki 2017 式58-66）——**默认关**（:825-860） | Maleki & Frigaard 2017 |
| 人工弥散 | 显式拉普拉斯 0.018/0.015（水泥）、0.012/0.012（前置液），dt 归一（:1079-1098） | 经验系数，非谱系公式 |
| 输运格式 | 半拉格朗日双线性回溯+五相过填归一（:1054-1070） | 数值格式 |
| 停泵 | 浓度场冻结（:1190-1206），无停泵滑移 | 工程假设 |
| 时间步 | CFL 自适应 cfl=0.5，dt_min=0.1s（:878-906） | T1-7 |

**③ 论文描述模型时必须写的公式清单**（供对照 agent 核对）：
1. 偏心间隙 b(φ,z) 与 e=clip(1−standoff,0.05,0.55)、体积校正；
2. 剪切率 γ̇=6|w|/b（含 w 一步滞后约定）与四本构表观黏度 η_i(γ̇)（含 μ∈[1e-5,3.0] 裁剪）；
3. 相混合三式：μ=Σφ_iη_i、τy=Σφ_iτy,i、ρ=Σφ_iρ_i（:613-634）；
4. 两层黏度闭包 η_mix=1/(c̄³/η₂+(1−c̄³)/η₁)（式4.23）；
5. 流动度与速度：pref=(b/b̄)²η_mix⁻¹I1(c̄,m)·[1+clip(ΔρI2/I1,±0.5)f_φ]·(1−wall)，w=q_half·pref/Σ(pref·b·dy)；
6. D2DGA 通量放大 f(c̄,m) 及 auto-m=η₁/η₂（clip[0.1,10]）；
7. I3 浮力通量 q=(ΔρH³/6η₂)I3(c̄,m)(−f_ξ,f_φ)+散度形式 dc/dt；
8. Papanastasiou 正则化 μ_reg=μ−τy/γ̇+τy(1−e^{−Mγ̇})/γ̇，M=100；
9. 屈服门：γ̇=6|w|/b、τw=G·b/2、immobile ⇔ τw≤1.15τy（+参考元不变量、整列冻结兜底）；
10. 人工弥散项（系数 0.018/0.015/0.012，每物理秒归一）与半拉格朗日格式；
11. 指标定义：η_E≡∫b·c/∫b（**即水泥占据率，非"被驱泥浆分数"**）、η_N 窄四分位、评价窗 b 加权均值、浮力数 b=(Δρ)gd²/(μ_PV·ŵ₀)、失稳/窜槽/混浆指数；
12. 1D 侧：体积推进到达时间、胶塞截断语义（stop≡cement_end=min(尾浆尾缘,顶替序列终点)）、Taylor 弥散三分支（TA/Batot/Fan&Wang+对流上限 α·U·R）、混浆增强 1+5·At·√(Re/100)、Atwood 重力修正+屈服抑制；
13. 边界耦合：鞋口时间线→相映射(lead/tail/spacer/mud/flusher)→环空入口；泵停冻结判据 Q≤1e-9。

**论文不可写/须改口径的**：流态（层/紊流）判别与 Re 修正（生产未启用）；隔离液参与流动度（闭包无此相）；紊流黏度项（占位零）；温度相关流变（无此物理）；"mobility clip [0.5,2]"（实际无裁剪）；e 截断 0.90（生产 0.55）。
