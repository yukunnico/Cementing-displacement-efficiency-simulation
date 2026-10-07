# Phase 4d 设计规格：七井温度场 + 表档重跑（2026-10-07）

> 上游 = 《详细执行计划_4d主批起_续作_2026-10-07.md》§4「4d 主批」+ Phase 4 设计规格
> `docs/superpowers/specs/2026-10-07-phase4-depthwise-temperature-design.md` §2「4d（七井场）」。
> 起点 HEAD `b70593c`；分支 `feat/temperature-coupling`；测试基线 17F/1202P/241sub。
> 依据纪律（§0.2）：锚点 = 现场资料原文（`temperature_pressure_profile.csv`）+ 用户裁定；
> 无新物理发明。

---

## §0 前置（4d-1，已完成）

`TableTemperatureField` 补 F-4 聚合入口（计划 §3-3 / §4-4d-1）：

- 新增 `T_column(md_values, t_s=0.0) -> np.ndarray`（与标量 `T` **逐位同值**；每批至多 1 条代表
  `ClampEvent`），只读属性 `oob_column_clamped_total`，`reset_audit()` 同时清该计数。
- `T()` / `oob_count` / `oob_events` 行为**一字不动**；`Constant`/`Geothermal`/`AnchoredProfileField`
  与 `load_delivered_pair`/`load_extended_pair_4d`/`assert_time_table_coverage` 全部不改。
- 测试：`tests/contract/test_temperature_field.py::TestColumnQuery`（9 条）；原
  `test_anchored_profile_field.py::TestExistingUnchanged::test_table_field_scalar_clamp_behavior`
  末行「无 T_column」状态快照断言随该前置翻转为「存在且与标量同值」。

---

## §1 井位与系数（用户裁定 2026-10-07，本窗口闭环）

**八井全集** = `hu101 / hu102 / hu103 / hu1 / hu2 / ht1_001 / ht1_003 / ht1_004`；
其中 **hu2 ≡ 呼探1-002（HT1-002）**（`hu2_loader.py` docstring 实证，现场资料目录名 `ht1_002_呼探1-002`），
仓库内无独立「呼2」井。

**七井口径裁定**：计划 §4-4d-2 写「其余四井」却列 5 名（含已裁定「新批一律不跑」的呼103）。
用户裁定（2026-10-07）：「**重点为呼101，呼1-003，呼1-004 这几口井**」，且非重点井 k 取
「notes 明示优先」（该问仅在非重点井纳入批跑时成立）⇒ **七井 = 八井 − 呼103** =
`hu101 / ht1_003 / ht1_004`（三重点井，叙事重点）+ `hu102 / hu1 / hu2 / ht1_001`（其余四井）。
**呼103 批外**（历史资产 T2/09-16/09-27 批照常引用，不重跑）。

**温度系数 k 与静温锚点集**（族 = 场景 A 电测/作业史静温；逐行原文见
`results/_probe_4d七井_2026-10-07/七井锚点盘点表.md`）：

| 井 | k | 静温锚点集 (md_m @ °C) | 依据 |
|---|---|---|---|
| 呼101 hu101 | **0.90** | 5700@123、7868@150 | 无明示系数；观测区间 argmin（用户裁定 2026-10-07） |
| 呼1-003 ht1_003 | **0.85** | 5290@123、7618@150 | notes 明示「领浆温度系数0.85」（用户裁定） |
| 呼1-004 ht1_004 | **0.85** | 5241@124、7660@155 | notes 明示「领浆温度系数0.85」（用户裁定） |
| 呼102 hu102 | **0.90** | 7120@147.8、7735@149 | notes 明示「水泥浆试验温度(0.9x温度系数)」（用户裁定：明示优先） |
| 呼探1-002 hu2 | **0.90** | 5292.5@111、7554@148 | notes 明示「循环温度133℃(系数0.90)」井底主段（用户裁定：明示优先；回接段另有 0.85，见 §1.1） |
| 呼探1-001 ht1_001 | **0.85** | 5460.159@110、5900@118、7000@137、7746@150 | notes 明示「温度系数0.85」（用户裁定：明示优先） |
| 呼探1 hu1 | **—** | **无锚** | 仅 1 个唯一锚深度（7601 m 三行 153.8/159/167 °C）⇒ 不可构造 ⇒ **Geothermal 回退（口径声明，不硬造锚）** |

**一致性自证**：上述归族规则（`anchor_inventory_20261007.py`，按序判定
邻井 > 出口 > 实验 > 静止/BHST > 电测）在三重点井上**逐字复现**其已裁定的场景 A 锚点集
（hu101 {5700,7868}、ht1_003 {5290,7618}、ht1_004 {5241,7660}）——与 4d 前期
`_probe_4d锚点_2026-10-07/三井锚点对账表.md` §2 一致。
**k 全部落在 [0.80, 0.95] 观测带内**（计划 §8-13 未触发）。

### §1.1 登记项（不阻塞，入册备查）

1. **呼探1-002 段落系数冲突**：井底主段明示 0.90、回接段明示 0.85（经验比值 94/111=0.847）。
   用户裁定「明示优先」，本规格取**井底主段 0.90**（域底 = 评价主体深度）；回接段差异
   （域顶 5292.5 m 一点）不单独建模。**若后续要分段 k，须扩 `AnchoredProfileField` 接口
   （超出 Phase 4 §1 改动面），届时先报后动。**
2. **呼探1-002 / 呼探1-001 的「顶部温度」行**（无「静止」措辞）未纳入锚点集，如实列在盘点表
   §2；其地温式残差 −1.70/−1.00/−0.60 °C（呼探1-002）与 −2.09/−1.83/−2.19 °C（呼探1-001）
   显示与地温线接近，**是否纳入属待裁**（本批不纳入）。
3. **呼103 数据完整但批外**：5750@125、7770@156（明示系数 0.85/0.9 双值）。若日后重启该井，
   锚点与系数材料已备。

---

## §2 改动面（超清单 = 先报后动）

| # | 文件 | 动作 |
|---|---|---|
| 1 | `cemdisp/data/temperature_field.py` | `TableTemperatureField.T_column` + `oob_column_clamped_total`（**已完成，§0**） |
| 2 | `scripts/entrypoints/run_sensitivity_current_20260916.py` | ①`TEMPERATURE_MODES` 增 `"anchored"` 与 `"table_ext"`；②新增模块级 `ANCHORED_WELLS` 注册表（§1 表）+ `EXT4D_*` 路径常量；③`build_temperature_fields` 增两档分支（anchored：1D/2D **各造独立实例**保审计分离；无锚井 ⇒ Geothermal + 备注；table_ext：`load_extended_pair_4d` + 仅 呼1-004） |
| 3 | `results/_probe_4d七井_2026-10-07/probe_seven_wells_20261007.py` | 七井批驱动（范式 = `_probe_屈服门四角_2026-10-07/probe_four_corner_20261007.py`：job-JSON + 子进程池 + hard 锚失败 exit 2） |
| 4 | `tests/contract/test_run_opts_wiring.py` 等 | 新档接线契约（全键穿透 / 无锚回退 / 非呼1-004 拒 table_ext / anchored 逐位可复现） |
| 5 | `cemdisp/models2d/annulus_d2dga.py` | **`_temperature_array` 接线 `T_column`**（4d-1 的消费方，见下方「改动面增补」） |

### 改动面增补：`_temperature_array` 接线聚合入口（0.5 行级，实施前声明）

**事实**（实测）：`_temperature_array`（`annulus_d2dga.py:1307-1320`）以
`np.array([field.T(float(md_i), float(t_s)) for md_i in md])` **逐点调标量 `T`**；该函数由
`geom` 初值（`:963`）与逐步刷新（`:1324`）调用。**4d 之前 `T_column` 零消费方**——即
F-4 聚合入口若只加不接，等于死代码，且锚域外节点仍会逐步逐点 append `ClampEvent`。

**本批的真实暴露**：呼101 锚域 `[5700, 7868]` 而环空域 `[5400, 7868]`、呼102 锚域
`[7120, 7735]` 而域 `[6823.10, 7735]` ⇒ 每步各有 ~30 / ~81 个锚域外节点（`fallback_geothermal=True`
按地温线取值，但**仍逐点记账**）。nz=250 × 数千步 ⇒ 10⁵ 量级 `ClampEvent` 无界增长，
正是 F-4 要防的形态（计划 §3-3「4d 主批表格档列化前必须补」）。

**改动**：`_temperature_array` 优先走 `col = getattr(field, "T_column", None)`；命中则
`col(md, float(t_s))`（**逐批聚合：每批至多 1 条代表事件**），未命中回退原逐点实现
（`Constant`/`Geothermal` 无该方法 ⇒ 逐位不变）。

**关2 依据**：`AnchoredProfileField.T_column` 与 `TableTemperatureField.T_column` 均已被契约测试钉为
与标量路径**逐位同值**（`test_anchored_profile_field.py::test_column_values_match_scalar_path`
`atol=0, rtol=0`；`test_temperature_field.py::TestColumnQuery::test_column_matches_scalar_bitwise`
`np.array_equal`）⇒ 数组内容不变，**仅审计记账粒度变化**（诊断字段，不进 summary ⇒ 关1 不破）。


**明确不动**：`_SWITCH_DEFAULTS`（一字不动）、`runners/`、`AnnulusInletState.temperature_c`、
Phase 1/2/3 已交付接口签名、`T()` 标量路径。

**不在本规格**：`enable_depthwise_temperature` 的 run_opts 接线（计划 §3-2 明定留到 **5d**）。
本批以探针 ctor 直装范式使用该开关（见 §3）。

---

## §3 批跑矩阵（r1.0，nz=250，CFL 自适应开）

统一底座 = `CORRECTED_KW`（`enable_yield_gate=True`/`enable_regime_split=True`/`enable_local_i3=True`）
+ `(F,F)` 角（`include_yield_term=False`、`enable_stream_yield_gate=False`）+ `pressure_mode="off"`
+ `pressure_caliber="shoe"` —— 与「后4a 口径排量锚」同底座，使数字可与
`results/_probe_排量锚后4a_2026-10-07/` 直接并列。

| 组 | 井 | 温度档 | 2D 逐列 | 用途 |
|---|---|---|---|---|
| A | 七井 | `off`（T-off） | — | 对照（T-off 逐位不受温度改动影响） |
| B | 七井 | `static`（统一地温式） | **开** | 无锚场对照（前 4d 口径） |
| C | 七井 | `anchored`（§1 表；呼探1 ⇒ Geothermal 回退） | **开** | **主产物**（锚定场逐井不同 ⇒ L3） |
| D | 呼101 | `static` / `anchored` | **关** | 与后4a 基线逐位对桥（hard 锚） |

- 组 A 7 + 组 B 7 + 组 C 7 + 组 D 2 = **23 run**；单 run 量级 50–115 s ⇒ 约 30–40 min 墙钟。
- **hard 锚**（组 D，判据 = 逐位）：呼101 × `static` × 逐列关 × r1.0 × (F,F) × pressure off
  ⇒ **η_N = 0.7709394595515454**（读
  `results/_probe_排量锚后4a_2026-10-07/jobs/hu101_r1.0_off_FF_post4a.json`，**禁手抄**；
  锚读文件、失败 exit 2 停批）。该锚证明新驱动零漂移。
- 组 B/C 的 2D 逐列开关 = ctor `enable_depthwise_temperature=True`（探针直装，计划 §3-2）。
- 口径标签：组 B/C/D 全属**后4a 口径**；引用 2e 四角/10-02 排量探针时必须标「前4a 口径」。

### §3.1 表档重跑（4d-4）

| 组 | 井 | 温度档 | 排量 | 说明 |
|---|---|---|---|---|
| E | 呼1-004 | `table_ext`（扩展表 333×362，表末 21600 s） | r0.6 / r0.8 / r1.0 | r0.6/0.8 解禁后重跑 |

- **`assert_time_table_coverage(…)` 在 `stop_t` 求得后、环空求解前调用**（批次级前置门）；
  报错即停（计划 §8-14，禁静默 clamp）。
- 附加 hard 判据：表档 run 的 1D/2D 场 `oob_count == 0`（真覆盖，非 clamp）。
- **必附声明**：「解禁 ≠ 洗脱名义排量热史错配（扩展段按末档 0.7 m³/min 外推，
  r0.6 实际为 0.6× 排量持续循环）」。
- **r0.5 维持禁引**（stop_t≈22819 s > 21600 s）。
- 扩展表文件位于 `results/_probe_4d扩表_2026-10-07/sandbox/HT1-004压力计算/`
  （`*_ext4d.xlsx` 受 `.gitignore` 的 `**/*.xlsx` 管辖不入库；由已入库的 `run_ext4d.m`
  + `HT1_004_T.m` 可重生成）。

---

## §4 产物（`results/_probe_4d七井_2026-10-07/`，目录带日期）

1. `七井锚点盘点表.md/.csv` + `anchor_inventory_20261007.py`（盘点，已产）。
2. `jobs/*.json` / `*.log` / `*.npz`（逐 run 观测）。
3. `七井温度场对账表.md/.csv`：逐井「锚点集 / k / 锚点处 T_anchored vs 实测 vs 地温式残差」。
4. `汇总表_七井温度场.csv/.md`：七井 × 组 A/B/C 的 η_E/η_N + 饥饿份额 + wall 占比 + stop_t +
   温度审计 + 关5 健康度（memo 尺寸 / oob / cfl_clip / dt 中位 / 步数）。
5. `一页解读.md`：井间分化叙事（接 2e）——组 C−组 B（锚定场相对统一地温式的净效应）、
   组 B/A（温度档效应）、组 D 与后4a 基线逐位对桥结论。

---

## §5 验收（计划 §5 五关 + L3）

| 关 | 判据 |
|---|---|
| 关1 | T-off 路径逐位：新增不得使 `results/_baseline_T_off/基线摘要_T_off.json` 位移（**读 JSON 比对，禁手抄**）；全量 `pytest tests/ -q` 预存 **17F/1202P/241sub** 不增，FAILED 清单与 `pytest_phase4b_full.log` 逐条一致 |
| 关2 | 新档缺省/不传场 ⇒ 逐位：`temperature_mode` 缺省 `"off"` 行为不变；`enable_depthwise_temperature=False` 缺省下逐位=HEAD；组 D hard 锚逐位命中（0.7709394595515454） |
| 关3 | 方向性：**锚定场逐井分化**——组 C 的 η 在七井间分化；三重点井锚点集与 `_probe_4d锚点` 场景 A 逐字一致；呼探1 = Geothermal 回退与其余井可区分 |
| 关4 | 专项锚：§1 锚点对账 + 表档 `assert_time_table_coverage` 通过 + 表档 oob=0；**措辞 = 四角 + 交互残差（禁「贡献分解」）** |
| 关5 | 数值健康度：memo 尺寸不随步数增长、oob 聚合计数有界、`cfl_clip=0`、dt 中位漂移 ≤2×；**红灯×关3 绿灯 ⇒ 停 + 报** |
| 两锚 | 排量响应锚（后4a 新基线三点单调 0.8205255446 / 0.7709394596 / 0.7224794436）不经本批重跑，但**组 D 逐位复核 r1.0 点**；**L3 判定** = 不同井的流体经公式分派 + 锚定场逐井不同 ⇒ 产出不同 η |

---

## §6 回退

全部 opt-in：`anchored`/`table_ext` 是新温度档名，缺省 `off` 行为不变；探针
`enable_depthwise_temperature` 默认 `False`。`git revert` 逐提交可回；批驱动与产物在
`results/` 日期目录内，不覆盖任何既有权威目录。

## §7 风险登记

1. **呼探1-002 段落系数**（§1.1-1）——取井底主段 0.90，回接段差异不建模；已登记。
2. **`table_ext` 依赖 gitignored 的 xlsx**——缺失即响亮报错（不静默降级到交付表，
   否则 r0.6/0.8 会被 clamp 污染，正是本波要消除的缺陷）。
3. **无锚井 Geothermal 回退**使呼探1 的「锚定场」与 `static` 档逐位相同 ⇒ 组 B/C 在该井
   必然重合；报表须显式标注（否则读成「锚定无效」）。
4. **组 B/C 开 2D 逐列**而组 D 关 ⇒ 跨组相减不可解释；组内比较有效。
5. **并行度**：Windows spawn 子进程池，`--workers 4`；任一 hard 锚失败即停批。
