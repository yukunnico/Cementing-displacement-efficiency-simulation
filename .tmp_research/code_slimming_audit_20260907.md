# 环空模型代码精简审计（2026-09-07，取证中间稿）

> 用户指令："仔细调研一下目前环空段模型代码有哪些是可以去掉的，模型要保持简练。"
> 审计对象：cemdisp/models2d/（2033 行）+ reporting/（1012 行）+ 关联 1D 侧。
> 审计维度：①生产口径是否消费 ②物理作用是否已被取代 ③包外引用面（删除破坏面）。

## 一、可直接删除（生产口径零作用、无包外实质依赖）

### 1.1 FLUSHER 整条链（约 60 行 + 测试面）
- 证据：hu102/ht1_004 入库 flusher=0.00 m³（`integrate_injection` 实测）；8 井 mean_flusher=0；
- ht1_004 loader 注释明说"合成 FLUSHER 仅用于验证序列可表达，现场无此流体"（ht1_004_loader.py:317-318）；
- 位置：annulus_d2dga.py（flusher 场初始化/入库累计/过填修正/平滑/快照五处 T1-6 标记）+ reporting/animation.py:58-65（spacer 回退零场）+ 摘要 mean_flusher 列；
- 风险：hu102/ht1_004 loader 返回 FLUSHER 流体对象（不进泵注序列）——删链需同时把 loader 的 FLUSHER 角色改 OTHER 或保留 _pick_fluids 的 flusher 槽位为 None。建议做法：2D 内不再建独立浓度场，FLUSHER 流体直接按 mud 相映射（与 DISPLACEMENT 同口径），loader 不动。
- 预期收益：run() 主循环减 1 个场、过填修正从五相回四相、npz/摘要各减 1 键。

### 1.2 R0 时代兼容参数/分支（构造器与 run() 内 3 处）
- `d2dga_viscosity_ratio` + `enable_d2dga_auto_m=False` 分支（annulus_d2dga.py:234,343,901,1192,1261）：生产恒 auto_m=True，标量 m 路径仅 3 个旧测试引用；
- `enable_d2dga=False` 主分支（:936,1188,1246 的三重 and 门）：只被消融脚本（closure_contribution_scan/_mass_balance stage2）使用——**保留但可合并为单一闭包开关**；
- `instability_decay_scale`（:238）：仅缩放失稳指数代理，非物理闭合，可用常数替换；
- 注意：删除会破坏 tests/test_improved_d2dga_annulus.py 若干断言与 2 个消融脚本——按"测试同步改+脚本标注 legacy"处理。

### 1.3 占位死量
- `mu_turbulent`（:770 零场, 850 注释"占位"）：12 元组返回中 1 席，仅兼容旧结果对象——可与 `Re` 一起降级为 summary 诊断（不随每步返回）；
- `mud_cake_thickness` 参数（:401 `del`）：纯兼容占位，删除需同步改 runner 调用（无）；
- 死 import：已清（Path/_dataclass_replace，3bd6bb0 后续清理）。

### 1.4 reporting 层可精简项（对 runner 产物瘦身，非求解器）
- `plot_final_fields_contour` 的"壁面泥饼"面板（contour_plots.py:271-272 有 >1e-12 才加面板的守卫，生产恒不触发）——保留守卫即可，无需动；
- `wall_field` 恒零时 reporting 仍做全形状校验（contour_plots.py:259-264）——微收益；
- `animation.py:58-65` spacer 快照回退零场逻辑——删 FLUSHER 后 spacer 快照可并入可选槽位。
- reference_figures（186 行）是"参考项目风格图件"，非论文必需——降为可选导出。

## 二、可合并（功能重复、口径统一后减参）

### 2.1 wall 判据三轨合一（最大简化点）
现状三轨：
1. `enable_yield_gate=True`（B2 物理屈服门，生产口径）→ `_yield_gate_wall`；
2. `enable_yield_gate=False` → c_min 浓度判据 + wall_seed_c_min + 整列解冻防护（:1312-1332，~30 行）；
3. `yield_gate_c_min_residual` 参数在 B2 分支其实不消费（判据里没有它，仅 _yield_gate_wall 签名收着）。
建议：生产口径只剩 B2；c_min 兜底轨是 08-23 时代遗产（B2 取代它已被 09-02 取证），改为 B2 关闭时直接报错或恒 wall=0，删 ~30 行 + 2 参数。
`yield_gate_c_min_residual` 直接删（签名与构造器同步）。

### 2.2 M2 流态迭代参数收敛
`regime_relax_alpha/max_iter/tol_rel/re_turb_ratio` 四参数只服务 enable_regime_split（生产**关闭**、层流元 R=1 中性）。若论文不写 M2，整个分支（annulus_d2dga.py:825-860 内 M2 段 + regime_closure.py 84 行）可标 legacy 移除；若写，则保留但把四参数收进一个 dataclass。
⚠️ 与论文口径联动：论文 1.7 表写 M2 开——若保持 M2 在论文叙事中，代码留着（生产不开、论文口径单独跑）。

### 2.3 指标重复计算
metrics 每步算 `bulk_fill` 与 `effective_efficiency` 两个同值量（:1379-1380，η_E≡cement_occ 恒等式）——终态保留一个即可，时序保留一列（另列由恒等式派生）。省一次全场均分（每步 O(ny·nz) trapezoid ×2→×1）。

### 2.4 snapshot 五元组
`lead/tail/spacer/flusher/wall` 五组快照各自 copy 全场——npz 落盘只消费 cement/spacer/wall（rerun8 系列 :203-212）。lead/tail 快照可关（论文用 cement 场已足）。省内存 ~40%（5→3 场）。

## 三、不能删（看似冗余实则有据）

| 项 | 表面冗余 | 保留依据 |
|---|---|---|
| `enable_cfl_adaptive`/固定 dt 双轨 | 双轨重复 | 论文 CFL 口径声明（三类数字引用规则）需要固定 dt 复现附表口径 |
| `wall_seed_c_min` | 似只服务兜底轨 | B2 轨虽不消费，但它是 C 根因修复的历史锚——随兜底轨一起删 |
| `alpha_cfl` I3 局部裁剪 | 与全局 CFL 并存 | 语义不同（I3 通量散度裁剪 vs 平流步长），T1-2 有单测 |
| `_low_tail_indicators` | 诊断小函数 | 论文窄边低尾指标来源（09-04 报告 §低尾） |
| `dispersion_dt_ref/dt_scale` | 双参数 | M1 归一开关（固定 dt 逐位复现基线的关键），论文弥散声明需要 |
| pump 停分支 `_compute_velocity` 调用 | "冻结了还算速度" | dt_step 需要 v 场做 CFL；且 T0-6 停泵诊断消费 w/mu |

## 四、执行建议（分两批）

**第一批（零风险，半天）**：
1. FLUSHER 链降级（2D 不建独立场，映射 mud 相）；
2. `yield_gate_c_min_residual` 删除 + c_min 兜底轨改恒 wall=0（保留开关行为文档）；
3. metrics 恒等列去重（bulk_fill 主留）；
4. 死 import 清理（已完成）。
验证：363 测试跑通 + hu103/hu101 两井重跑与 09-06 基线逐位对照（expected：η_E/η_N 不变——第一批全部无物理效应）。

**第二批（需裁定，1-2 天）**：
5. R0 兼容分支（auto_m=False 路径 + 标量 m）删除或移 legacy 模块（3 个旧测试同步改）；
6. M2 分支去留（与论文 1.7 口径联动裁定）；
7. lead/tail 快照关闭（npz 减两键，报告脚本同步）；
8. reporting 层 wall 守卫微调 + reference_figures 降可选。

**不建议动**：`enable_d2dga` 主开关（消融链依赖）、CFL 双轨、弥散双参数、I3 局部化开关（09-06 后生产默认虽关，但 I3 局部化是 corrected 口径成员，论文对照需）。

## 五、量化预期

- annulus_d2dga.py：1549 行 → 约 1380 行（−11%）；构造参数 36 → 约 30；
- 每步计算量：−1 个浓度场平流/过填/平滑（FLUSHER）+ 恒等列去重 ≈ 每步省 8~12%（单井 70s → ~62s，8 井重跑省 ~1.5 分钟——绝对值小，主要是代码清晰度收益）；
- reporting：不动求解，仅产物瘦身。
