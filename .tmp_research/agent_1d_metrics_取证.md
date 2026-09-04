# 1D 管内模型 + 胶塞语义 + 指标管线取证报告（2026-09-03）

> 取证范围：仓库 `D:\users\desktop\research\控压固井项目\cement model`（main @ 396dc77）。
> 只读取证，未改任何代码。所有行号以当前工作区代码为准。
> 结论速览见 §6；对"另一会话 09-03 结论"的逐条裁定见 §3.4。

---

## 1. 胶塞语义取证（git show 21d4490）

### 1.1 胶塞在模型里是什么 —— 不是实体，是"序列截断"

`git show 21d4490` 共改 4 文件：`casing_flow.py`(+112/-19)、`tests/test_plug_semantics_restart.py`（新增 312 行 5 例）、`tests/test_casing_flow.py`、`tests/test_pipe_capacity_chain_fix.py`。

**胶塞不是界面、不是边界条件、没有任何动量/几何实体。** 它唯一的存在形式是
`cemdisp/transport1d/casing_flow.py:843-889` 的 `_displacement_sequence_cutoff`（纯静态方法）：

- 顶替序列 = 首个 `PumpingStageEvent.RESTART` 步**之前**的全部注入步骤（casing_flow.py:870-877，按 `s.step.event_tag == PumpingStageEvent.RESTART` 找索引）；
- 截断点时间 = 首个 RESTART 步的 `start_time_s`（casing_flow.py:881）；无 RESTART 井 = 全日程最晚结束时刻（878-880，与旧口径一致）；
- 方法幂等（883-888 为"首步即 RESTART"兜底）。

即：**胶塞面 ≡ 尾浆尾缘**（尾浆与其上方注入流体的体积坐标界面）。模型里不画胶塞、不判胶塞位置，只把 RESTART 及其后的体积从"累计泵入体积"中剔除。

### 1.2 "替浆顶胶塞驱动尾浆"如何驱动 —— 体积坐标链式推进

1D 是不可压缩体积推进法：每个后注入步的累计体积坐标推动前方所有界面（casing_flow.py:903-934 `_front_arrival_time` / `_rear_arrival_time`，目标体积 = 该步累计体积 ± 管容）。尾浆尾缘 = 后继步（压塞液/后置液/替浆）前缘，同一物理界面、同一体积坐标。胶塞语义没有改变驱动机制——替浆仍通过体积坐标链推动尾浆尾缘，只是把 RESTART 段的体积剔除。

### 1.3 "替浆不进环空"如何保证 —— 累计封顶 → 替浆前缘永不到达

截断后累计泵入体积封顶于顶替序列总量。替浆前缘的到达体积坐标 = cum_start(替浆) + 管容 > 封顶值，故：

- `_front_arrival_time` 返回 None → `_build_shoe_timeline` **不生成替浆 FRONT_ARRIVAL 鞋口事件**（casing_flow.py:695-701，`if front_time_s is not None` 才 append）；
- run() 的 fronts 列表里替浆前缘取 `pumping_end_time_s`（=截断终点）作为占位上界（casing_flow.py:258-261）——该占位只被 runner 的 stop 回退扫描消费，而回退仅在 cement_end 缺失时生效（终跑 build_case 直接取 cement_end，`_mass_balance_diag_20260902.py:60`）。

因此鞋口时间轴里替浆相分数恒为 0，"替浆不进环空"由时间线构造自动保证，不依赖 has_plug。

### 1.4 RESTART 后处理步如何与主泵注解耦

三处一致截断（幂等防御）：
1. `run()`：casing_flow.py:245 `scheduled_steps, pumping_end_time_s = self._displacement_sequence_cutoff(...)` —— 界面推进、cement_end、时间轴全部用截断序列；pumping_end = 顶替序列终点（245-248）；
2. `pipe_exit_state_at()`：casing_flow.py:363-365 再次截断（防存量调用方传入全序列）；
3. `_build_shoe_timeline()`：casing_flow.py:688-690 再次截断。

RESTART 步不产生任何鞋口事件（已不在序列里）；时间轴 END 事件 = 最后一个截断步的 end_time（casing_flow.py:709-710）。实际发生截断时 result.notes 写入胶塞注记（casing_flow.py:320-328；hu102 复核确认存在）。

### 1.5 has_plug 与胶塞语义完全解耦

`has_plug`（casing_flow.py:144,178-187,536-537）只承载"混浆增强因子=1"（胶塞刮拭阻止界面混浆 → `_interface_instability_factor` 直接 return 1.0）。**8 井生产链没有任何 runner/脚本设 has_plug=True**（runners/ 目录 grep 无匹配），即 1D 过渡带生成时混浆增强**仍然生效**。胶塞隔离语义由截断无条件保证，与 has_plug 无关——这是有意设计（21d4490 提交信息明示），但注意：**"有胶塞 → 混浆增强=1"这条语义在 8 井生产配置里并未启用**，1D 过渡带宽度仍受混浆增强（≤10×）加宽。

### 1.6 hu102 数值复核（本会话实测，1D-only）

hu102 日程末 5 步：尾浆 7.0m³@0.7 → 压塞液 7.0@0.6 → 后置液 6.0@0.7 → 替浆 74.0@0.841 → **循环排混浆 41.0@1.3 (RESTART)**。
实测：`pumping_end_time_s = cement_end_time_s = 11737.14s` = 替浆步末 = RESTART 起点；管容 88.9m³；notes 含胶塞注记。
尾浆尾缘目标体积 = cum_end(尾浆)+88.9 > 截断序列总量（= cum_end(尾浆)+87）→ 尾缘到达为 None → stop 封顶于替浆步末，尾浆尾段 1.9m³ 滞留管内。与提交信息"11824.84→11737.14s（-87.69s）"逐位吻合。

---

## 2. 1D→2D 耦合点与 stop 时刻定义

### 2.1 四个文件的分工

| 文件 | 角色 |
|---|---|
| `pipe_exit_state.py` | 纯数据类 `PipeExitState(time_s, flow_rate_m3_s, stage_name, phase_fractions)`（pipe_exit_state.py:23-34） |
| `shoe_timeline.py` | 事件枚举 + `ShoeTimeline.at(t)` = "最近过去事件"查询（shoe_timeline.py:111-157） |
| `casing_flow.py` | 事件构造（截断序列 → 鞋口事件 + 弥散过渡带）+ stop 计算 |
| `boundary_bridge.py` | 8 井链实际消费的桥：`_legacy_provider`（boundary_bridge.py:179-206） |

`_legacy_provider` 逻辑（boundary_bridge.py:179-206）：泵注期（Q>1e-9）取 `shoe_timeline.at(t)`（含轴向弥散/混浆增强生成的多相过渡带）；停泵期回退 `casing_solver.pipe_exit_state_at()`（带停泵沉降增强，casing_flow.py:1136-1202）。相映射 `split_cement_phases=True`：LEAD/INTERMEDIATE→lead、TAIL→tail、其余水泥→cement（boundary_bridge.py:72-97），2D 侧 `inlet_tail = tail+cement`（annulus_d2dga.py:1009）。

### 2.2 stop（停泵/碰压）时刻定义

**stop ≡ `cement_end_time_s`**（casing_flow.py:287-318）：

- 对最后一段水泥（LEAD/INTERMEDIATE/TAIL 步）取尾缘过鞋口时刻：尾缘 = 后继流体前缘（同一物理界面、同体积坐标），重力修正按 (后继流体， 尾浆) 配对（casing_flow.py:300-310，F2 修复）；尾缘不可达时取 pumping_end 上界（296-299）并再封顶（311-314）；
- 安全约束"只上抬不下压"：cement_end ≥ max(水泥前缘到达时刻)（316-318）；
- 胶塞截断后 pumping_end = 顶替序列终点，故 **cement_end = min(尾浆尾缘过鞋口， 顶替序列终点)** ——碰压成功井不变，未碰压井停在替浆步末。

### 2.3 hu102 未碰压分支

替浆 74m³ < 管容 88.9m³（hu102_loader.py:82,95；hu102_loader.py:19-20 注记"到量未碰压（单流阀失效）"）→ 尾浆尾缘目标体积超出封顶 → stop = 替浆步末 = 11737.14s，尾浆尾段 1.9m³ 滞留管内。2D 入环空水泥 18.74→16.86m³，η_E 0.9846→0.9836（汇总.csv 对照列，-0.0010）。

### 2.4 2D 侧停泵冻结

2D 内部另有一层停泵判定 `pump_active = inlet_state.flow_rate_m3_s > 1e-9`（annulus_d2dga.py:1016），停泵即冻结浓度场（1190-1206，注释"水泥静凝胶强度在短时间（<2h）内足以抵抗浮力滑塌"）。在 stop=cement_end 口径下，2D 求解窗 [0, stop] 内的日程内停泵段（如 ht1 井）会被冻结，停泵后无进一步顶替演化。

---

## 3. η_E/η_N 指标管线

### 3.1 η_E 精确定义（annulus_d2dga.py:1233-1236，确认）

```python
cement = np.clip(lead + tail, 0.0, 1.0)   # 1233
eff = cement                                # 1234
bulk_fill = _trapez2d(geom["b"] * cement, geom) / half_volume         # 1235
effective_efficiency = _trapez2d(geom["b"] * eff, geom) / half_volume # 1236
```

- η_E ≡ 全域 b 加权水泥占据率（库存比）。`half_volume = _trapez2d(geom["b"], geom)`（annulus_d2dga.py:969，体积重标定后 ≡ 0.5×物理环空体积）。
- **确认**：公式位置 1233-1236 正确。
- **纠正 1（时刻）**：1233-1236 是**每个时间步**都执行的行内指标；"停泵时刻"并非公式自带，而是两个条件合成——"摘要取 metrics 末行"（annulus_d2dga.py:1314 `final = metrics.iloc[-1]`）× "8 井链 total_t = cement_end_time_s = stop"（`_mass_balance_diag_20260902.py:60` + `rerun8_after_plug_semantics_fix_20260903.py:148` 的 `AnnulusD2DGASolver(total_t=stop, ...)`）。求解器自身口径是"total_t 时刻"（默认 12000s）。
- **纠正 2（口径）**：`eff = cement`（1234）——η_E 与 `bulk_cement_fill`（最终水泥浆占据率）是**同一个数**：摘要里"全井段最终有效顶替效率"与"最终水泥浆占据率"逐位相等（annulus_d2dga.py:1348 vs 1350；hu1 摘要.json 实测两者均 0.9659）。η_E 是**体积占据率**口径，不是浓度质量口径：任何浓度 c>0 的水泥（含弥散光晕/稀释过渡带）都全额计入，只有泥浆才拉低它。
- 深度剖面/评价窗同理：`_evaluation_window_efficiencies`（75-92）与 `_depth_profiles`（908-936）都在 total_t 时刻的最终场上计算。

**占据率的 b 加权：窄边残留被降权。** η_E = ∫b·c/∫b，窄边 b 小 → 窄边残留泥浆对 η_E 的拉低按 b 加权后被稀释。这是 η_N 存在的理由；8 井报告口径以 η_E 为主指标时，居中度效应天然被稀释（A4 报告所称"η_E 指标稀释"的机制，在当前代码就是 b 加权 + 占据率口径本身）。

### 3.2 η_N（窄四分位效率）

`_narrow_quarter_efficiency`（displacement_metrics.py:211-225）：方位角最后 1/4 行（ny=40→10 行）的 ∫b·c/∫b；窄边=第 0 行的对面（h=H(1+e·cos(πφ))，φ=1 间隙最小，与 front_narrow=cement[-1] 约定一致）。summary 字段 `eta_narrow`（annulus_d2dga.py:1361）取**最终场**（total_t 时刻）。

### 3.3 其他指标位置

| 指标 | 位置 | 口径 |
|---|---|---|
| 窄四分位 η_N | annulus_d2dga.py:1349/1361 ← displacement_metrics.py:211-225 | 最终场窄 1/4 方位 ∫b·c/∫b |
| 窜槽（channeling） | annulus_d2dga.py:1242-1245 | \|front_wide−front_narrow\|/域长；前缘=c≥0.5 最深位置（1238-1240） |
| 混浆指数 | annulus_d2dga.py:1246 | ∫b·4c(1−c)/域体积 |
| 失稳指数 | annulus_d2dga.py:1247-1250 | channeling×流度比×(1+0.4·mixing) 的 1−exp(−·/5) |
| 前缘到位率 | rerun8 脚本:158,170-174 ← metrics front_wide/narrow_m | c≥0.5 前缘/域长 |
| 壁面冻结占比 | annulus_d2dga.py:1266 | mean(wall) 末行 |
| 低尾指标 | annulus_d2dga.py:95-105 | standoff<0.5 段占比 + 窄边 c<0.05 域占比 |
| 评价窗效率 | annulus_d2dga.py:75-92 | 最终场上每评价窗 eta_E/eta_N |
| tier0 诊断 | annulus_d2dga.py:1395-1402（try/except 记 error） | 纯后处理注入 summary |

### 3.4 对"另一会话 09-03 结论"的逐条裁定

| 断言 | 裁定 | 证据 |
|---|---|---|
| 胶塞修复对 7 井 η_E 逐位无影响 | **成立** | 汇总.csv Δη_E 列 7 井全为 0.0 |
| hu102 η_E −0.001pp、stop −87.7s | **成立**（-0.0010 绝对值 ≈ -0.10pp） | 汇总.csv；本会话 1D 复核 stop=11737.14s=替浆步末 |
| η_E = 停泵时刻全域 clip(lead+tail) 库存比 @1233-1236 | **部分成立**：公式对（1233-1236），但"停泵时刻"由"metrics 末行 × total_t=stop"两条件合成（1314 + `_mass_balance_diag_20260902.py:60` + rerun8:148）；且 η_E≡水泥占据率（双名同值），占据率口径正是居中度"体现不出来"的第一机制 |
| "runner 取 metrics 末行" | **措辞纠正**：取 metrics 末行的是求解器（annulus_d2dga.py:1314）与终跑脚本（rerun8:155,163）；runner 本身只转写 result.summary（hu102_tailpipe.py:83-90） |

---

## 4. 居中度输入链

### 4.1 8 井 standoff 原始值来源

| 井 | 来源 | 值 | 证据 |
|---|---|---|---|
| hu1 | model_assumption 均一 | 0.65 | hu1_loader.py:295-308 |
| hu2 | 设计值（斯伦贝谢模拟），鞋底略降 | 0.78→0.76 | hu2_loader.py:281-295 |
| hu101 | 名义假设（悬挂器失效座底叙事）；实测剖面存在但默认不用 | 0.38–0.48 | hu101_loader.py:259-262（assumed）、249-258（实测 0.22–0.88）、265-274（默认 None→assumed） |
| hu102 | 分段代理 × 间隙因子 clamp[0.30,0.85] | 0.65/0.70/0.68 | hu102_loader.py:181-203 |
| hu103 | 设计代理 77.8% − 深度修正 × 间隙因子 clamp[0.30,0.90] | ≈0.748–0.778 | hu103_loader.py:74,183-204 |
| ht1_001 | 设计模拟值 | 0.804 | ht1_001_loader.py:101-102 |
| ht1_003 | 设计模拟值 | 0.83 | ht1_003_loader.py:380-391,417 |
| ht1_004 | 设计文档 6.3 节模拟值，全井段固定 | 0.83 | ht1_004_loader.py:81,173-180 |

### 裁定：standoff 全部为假设/设计值（hu101 实测未用）

8 井中 7 井为"设计软件模拟值/邻井代理/分段代理"；唯一有实测的 hu101（从居中度检测图 Pipe Standoff 读数：扶正器间 0.22–0.78、扶正器处 0.60–0.88）**默认不用**（`measured_standoff=None` → `_ASSUMED_STANDOFF`，hu101_loader.py:265-274；终跑 build_case 调 `load_hu101_tailpipe()` 无参）。**09-03 终跑 8 井无一使用实测 standoff。**

### 4.2 从原始值到模型 e 的映射链（"六层栈"实为 7 层）

| 层 | 位置 | 变换 |
|---|---|---|
| L1 loader 构造 | 各 loader | 假设/设计值 → DepthValuePoint 剖面 |
| L2 well_spec 传递 | well_spec.py:121 | standoff_profile 元组 |
| L3 深度插值 | annulus_d2dga.py:393,397 | np.interp 到 nz=250 网格 |
| L4 e 截断 | annulus_d2dga.py:404 | e = clip(1−standoff, 0.05, e_clip_max=0.55) |
| L5 方位角形状 | annulus_d2dga.py:412-416 | h = half_gap·(1+e·cos(πφ))，b = 2h |
| L6 体积重标定 | annulus_d2dga.py:431-435 | scale→H,b 同乘（e 对平均间隙的贡献被消掉，只留对比度） |
| L7 动力学消费 | annulus_d2dga.py:784 + 810-811 | base=(b/b̄)²/η_mix；浮力方位修正 clip(±0.5) |

补充：hu101 实测剖面未进入 09-03 终跑；hu102/hu103 的间隙因子（hu102_loader.py:202、hu103_loader.py:203）在 L1 内部已混入；e 截断下限 0.05（同心下限）。

### 4.3 8 井换算 e 值（1−standoff，clip 前）

hu101 0.52–0.62（clip 到 0.55，唯一大 e 井）；hu1 0.35；hu102 ≈0.30–0.35；hu2 0.22；hu103 ≈0.22；ht1_001 0.196；ht1_003 0.17；ht1_004 0.17。
即 **7 井 e≈0.17–0.35，只有 hu101 e 饱和在 0.55**。e_clip_max=0.55 为默认值（annulus_d2dga.py:253）——终跑未设 0.90（"生产跑道 0.90"仅出现在 docstring 与消融脚本中）。

---

## 5. runner 口径差异 + 09-03 终跑配置

### 5.1 9 个 runner：8+1 结构

9 个 runner = 8 个 `*_tailpipe.py` + `ht1_004_ablation.py`（消融 runner）。8 个 tailpipe runner **同构**（唯一差异：hu101 有 `_schedule_total_time_s + 20min` 回退分支 hu101_tailpipe.py:48-51,107，但终跑路径 total_t_s 恒传入不触发）：

- `CasingFlowSolver(enable_gravity=True)`（其余默认：dt=2.0、弥散开 α=0.25、Atwood 浮力开、混浆增强开、has_plug=False）；
- `build_coupled_annulus_inlet_provider(..., split_cement_phases=True)`；
- `annulus_stop_time_s`：F2 口径优先 cement_end_time_s，前缘扫描回退（hu102_tailpipe.py:148-183 等；ht1_004_tailpipe.py:171-204 与其完全一致）；
- `AnnulusD2DGASolver(total_t=stop, nz=250)`（其余全默认，ny 默认 40）。

### 5.2 09-03 终跑不经过任何 runner，而是复刻 runner 配置

09-03 终跑 = `scripts/rerun8_after_plug_semantics_fix_20260903.py` + `scripts/_mass_balance_diag_20260902.py:52-64 build_case`：

- `CasingFlowSolver(enable_gravity=True)`（_mass_balance_diag_20260902.py:56）；
- `stop = float(casing_result.cement_end_time_s)`（:60）；
- provider split_cement_phases=True（:58-59）；
- `AnnulusD2DGASolver(total_t=stop, nz=250, ny=40)` 其余全默认（rerun8:148），即：
  - **弥散**：2D 弥散 0.018/0.015 每 4s（dispersion_dt_scale=1.0）；CFL 自适应开（enable_cfl_adaptive=True，dt 4s→自适应，dt_min=0.1）；
  - **屈服门**：enable_yield_gate=True（09-02 起默认，可逆 τw 重建壁面层），c_min=0.05、wall_seed_c_min=0.005；
  - **e_clip_max=0.55**（默认；未设 0.90）；
  - open_outlet=True（开放出口，水泥可流出域顶）；
  - enable_d2dga / i3_flux / auto_m / true_buoyancy 全 True；enable_regime_split=False；enable_local_i3=False。

### 5.3 runner 与终跑口径等价性

终跑配置 = runner 配置（除 hu101 未触发回退）。准确表述：**终跑复刻 runner 配置但不调用 runner 代码**（提交信息明言"复用 _mass_balance_diag build_case（stop=cement_end_time_s）… 不改 cemdisp 包内任何代码"）；8 个 tailpipe runner 的直接产物在 results/呼X尾管_1D2D耦合模型/ 目录，与终跑目录相互独立（旧 runner 产物目录 396dc77 另有一次 chore 落盘）。

---

## 6. 结论：可能掩盖居中度效应的候选机制（按证据强度排序）

09-03 终跑事实矩阵（汇总.csv + 前缘到位率.csv）：

| 井 | η_E | η_N | 窄边前缘/域长 | 壁面冻结占比 | 库存比(设计水泥/环空) | e≈(1−standoff) |
|---|---|---|---|---|---|---|
| hu101 | 0.5859 | 0.0456 | 49.6/2468 = 0.020 | 0.4427 | 1.068 | ~0.55（饱和） |
| hu1 | 0.9659 | 0.8729 | 3963.1/4077.7 = 0.972 | 0.0093 | 1.002 | 0.35 |
| hu102 | 0.9836 | 0.8853 | 1.0 | 0.0 | 1.374 | ~0.30–0.35 |
| hu103 | 0.9926 | 0.9667 | 1.0 | 0.0 | 1.107 | ~0.22 |
| hu2 | 0.9846 | 0.9910 | 1.0 | 0.0 | 1.134 | 0.22 |
| ht1_001 | 0.9898 | 0.9780 | 1.0 | 0.0 | 0.981 | 0.196 |
| ht1_003 | 0.9984 | 0.9943 | 1.0 | 0.0 | 1.161 | 0.17 |
| ht1_004 | 0.9985 | 0.9938 | 1.0 | 0.0 | 1.271 | 0.17 |

（6/8 井宽边前缘=窄边前缘=域长，即前缘指标饱和、channeling=0；hu1 窄边 0.972；hu101 窄边 0.020。）

**hu101 是活的对照井**：e 饱和 0.55 + 悬挂器失效叙事 → η_N=0.046、窄边前缘只到 49.6m/2468m、壁面冻结 0.44 → **居中度效应在 e 够大时是能出现的**。所以问题不是"模型里 e 没接上"，而是"7 井 e 太小 + 指标口径把剩余效应稀释掉"。

### Tier 1（证据充分，直接从代码+数据读出）

1. **[输入数据问题] 8 井 standoff 无一实测**（§4.1）：7 井为设计/代理值（0.65–0.83），hu101 实测（0.22–0.88）存在但 09-03 终跑未用。输入端 e 对比度先天不足，且"实测 vs 假设"从未在同一链路里对照进 8 井口径。
2. **[指标定义问题] η_E ≡ 水泥占据率（双名同值）**（annulus_d2dga.py:1233-1236,1348,1350）：任意浓度 c>0 全额计入，壁面残泥/弥散光晕都算"顶替成功"；b 加权让窄边（b 小）残留降权。居中度主要影响"窄边何时被扫过"，而停泵时刻域内已无大片未扫区 → η_E 不敏感是口径的必然结果。η_N/窄边前缘/壁面冻结才是 e 的敏感指标（hu101 η_N=0.046 vs 其他 ≥0.87）。
3. **[输入数据问题] 库存比 V_cem_design/V_ann ≥1（除 ht1_001 0.981）**：水泥设计量≥环空体积，open_outlet=True（annulus_d2dga.py:246,120-127 open_outlet 时 `_limit_phase_volume` 直接 clip 放行）→ 停泵时刻域内几乎全被扫过，η_E→1。评价域相对水泥量"太小"，是 η_E 接近 100% 的第一性原因（汇总.csv 库存比列 0.981–1.374）。
4. **[物理上本就如此] stop=尾缘过鞋口≡碰压**：碰压成功井（7 井）停泵时刻=全部水泥已入环空，宏观体扫效率物理上就应接近 1；CBL 反映的是候凝后评价段胶结质量（口径错位是项目既有结论）。hu102 未碰压如实降 0.001——链路对"少进 1.9m³ 水泥"响应方向正确，但幅度被占据率口径压扁。

### Tier 2（机制在代码里，量级未独立取证）

5. **[模型结构问题] e 只以"对比度"进入动力学**：L6 体积重标定（431-435）消掉 e 对平均间隙的影响，e 仅通过 (b/b̄)²（784）与浮力方位修正 clip(±0.5)（810-811）起作用；7 井 e=0.17–0.35 → 宽窄流度比 (1+e)²/(1−e)² ≈ 2.0（ht1_004）至 4.3（hu1），瞬时窜槽存在但 stop 时刻已被扫平。
6. **[模型结构问题] η_N"被吞"链路仍在**：双层黏度闭包 1/η_mix = c̄³/η₂+(1−c̄³)/η₁（778-779），稀释水泥具有泥浆流度 → 该机制放大宽窄差异（本应放大 e 敏感性），但被 Tier1-2 的占据率口径盖住。
7. **[数值计算问题] 弥散+过填归一+屈服门解冻把残留压到极小**：2D 弥散 0.018/0.015（257-258）+ 五相过填归一（1093-1098）+ wall_seed_c_min=0.005（252,304-309）+ M3 可逆屈服门（254,1164-1167）→ 6 井壁面冻结占比=0（汇总.csv），残留泥浆薄到只剩弥散光晕（而光晕计入水泥）。
8. **[数值计算问题] 前缘指标饱和**：`_front` 取 c≥0.5 的最深位置（1238-1240），6/8 井停泵时宽窄边前缘都精确=域长 → channeling=0，窜槽/失稳指数在 η_E≈1 井上不携带 e 信息（仅 hu1 0.0281 与 hu101 有分离）。
9. **[指标定义问题] 1D 过渡带以"最近过去事件"进入 2D**：ShoeTimeline.at() 最近过去语义（shoe_timeline.py:111-157）+ 1D 弥散过渡带（casing_flow.py:562-670）把尖锐前缘抹成 S 形，入库完成度 <1（入库完成度=V_inj/V_design，`_mass_balance_diag_20260902.py:97`）——但占据率口径对入库稀释不敏感（稀释水泥仍计入）。

### Tier 3（已知坑的现状裁定）

10. **[数值计算问题] 鞋口时间轴 REAR_EXIT 尾缘配对残坑**：run() 的 cement_end 尾缘配对已由 F2 修复（300-310），但 `_build_shoe_timeline` 的 REAR_EXIT 事件仍按"本步流体 vs 上一步流体"配对（casing_flow.py:702-708 用 `displaced_fluid`=上一步流体），与同界面 FRONT_ARRIVAL（697-701）和 cement_end 配对（300-310）不一致。影响面：**仅时间轴事件时刻**——REAR_EXIT 事件相位由体积查询覆盖（phase_override=None → `_pipe_exit_state_from_volume`），不改变 2D 入口相分数/排量/stop；全仓 grep 确认 REAR_EXIT 无任何专门消费者（仅 casing_flow.py:708 产生、shoe_timeline.py 定义），影响限于事件排序与诊断读取，**不改变任何 8 井结果数字**。项目记忆中"cement_end_time_s 尾缘配对 bug 未修"为 08-22 旧状态，F2（09-01）已修，现残坑降级为本条。
11. **[模型结构问题] 8 井链 has_plug=False**：胶塞刮拭（混浆增强=1）语义未启用，1D 过渡带仍被混浆增强加宽（casing_flow.py:536-560,611-618）。与胶塞隔离语义（21d4490）实现独立，后者已生效。
12. **[物理上本就如此] 日程内停泵冻结**：pump_active=Q>1e-9（1016）→ 停泵段浓度场冻结（1190-1206），停泵后窄边不再演化；对日程内有长停泵的井（ht1 系列）窄边差距被冻结，但 stop=全部水泥入环空后域已被扫过，8 井上量级小。

### 面向用户的直接回答

- **"管内段 1D 模型有没有问题？"**——1D 侧没有发现会系统性抬高 η_E 的计算错误：胶塞语义修复语义正确且 7 井为无操作（Δη_E=0 证实）；遗留坑只剩时间轴 REAR_EXIT 配对（不影响结果）。1D 的"局限"在于它是体积坐标活塞流+经验修正（重力乘子、弥散过渡带），不模化管内偏心、不模化塞面初始扰动，进 2D 的边界与居中度无关。
- **"η_E 指标口径有没有问题？"**——**有，且是主因**：η_E≡水泥占据率（任意浓度计入+b 加权降权窄边+开放出口+水泥量≥环空量）在"库存比≥1+碰压停泵"设定下必然趋向 1。居中度效应应改看 η_N/窄边前缘/壁面冻结/评价窗效率；窜槽指标需修前缘饱和（加"前缘到达时间差"或浓度剖面指标）。分离居中度效应的具体抓手：① hu101 用实测 standoff（between_centralizers 剖面 0.22–0.78）重算对照；② 对 e 做参数扫描（e_clip_max 放开或 standoff 缩放）在同链路验证 η_N 响应；③ 报告口径并列 η_E 与 η_N（P0 口径并列既有结论）。
- **"体现不出来"的另一重含义**：8 井 η_E 0.97–0.999（hu101 除外）与 CBL 评价段质量（通过率 0.63–0.79）严重错位，根源是"停泵时刻体扫占据率 vs 候凝后声幅评价段质量"的口径错位（项目既有结论：结构性必然，非简单相关）。

---

## 附 A：本次取证产物

- hu102 1D 复核输出（本会话实测）：pumping_end=cement_end=11737.14s，notes 含胶塞注记，管容 88.9m³。
- 全部行号基于 396dc77 工作区。

## 附 B：关键 file:line 索引

- 胶塞截断：cemdisp/transport1d/casing_flow.py:843-889（`_displacement_sequence_cutoff`）；run() 接线 245-248；pipe_exit_state_at 363-365；时间轴 688-690；notes 320-328
- 停泵/碰压 stop：casing_flow.py:279-318（cement_end_time_s，F2 配对 300-310）
- 鞋口时间轴构造+弥散：casing_flow.py:672-731（REAR_EXIT 配对残坑 702-708）；弥散 562-670；混浆增强 508-560
- 1D→2D 桥：cemdisp/models2d/boundary_bridge.py:179-206（_legacy_provider）；相映射 72-97
- stop→total_t：scripts/_mass_balance_diag_20260902.py:52-64（build_case）；scripts/rerun8_after_plug_semantics_fix_20260903.py:148
- η_E 公式：cemdisp/models2d/annulus_d2dga.py:1233-1236；摘要取末行 1314；summary 1342-1365；评价窗 75-92；η_N 1361 ← cemdisp/diagnostics/displacement_metrics.py:211-225
- 2D 停泵冻结：annulus_d2dga.py:1016,1190-1206
- 速度场/b² 流动度/B1 守恒权重：annulus_d2dga.py:772-876（784、822-863）
- 几何/e/standoff 七层栈：annulus_d2dga.py:362-438（393,397,404,412-416,431-435,437）
- 8 井 standoff 来源：hu1_loader.py:295-308；hu2_loader.py:281-295；hu101_loader.py:249-274；hu102_loader.py:181-203；hu103_loader.py:74,183-204；ht1_001_loader.py:101-102；ht1_003_loader.py:380-391,417；ht1_004_loader.py:81,173-180
- 8 井 loader 注册：scripts/_mass_balance_diag_20260902.py:39-48
- 终跑脚本：scripts/rerun8_after_plug_semantics_fix_20260903.py
- 测试：tests/test_plug_semantics_restart.py（312 行 5 例）
