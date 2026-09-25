# 官方口径对照：`rerun_all_wells_corrected.py` vs 生产 runner

> 内部自洽加固计划 **Task 7 Step 4** 交付物。**只写文档，不改任何脚本、不改任何口径**。
> 术语按仓库根 `CONTEXT.md`：**管内段**（不写"套管段/套管内"）、**环空段**（不写"环空"）、
> **环空二维 = 轴向 × 方位角、无径向维**、**前缘/尾缘成对且定义在管内段**
> （故本文描述环空水泥尖端时**不**使用"前缘"）。
> 按附加裁定 5：**行号以本文写作时实际读到的代码为准**，每条差异给出两侧可核对引用。

- 生成日期：2026-09-26；写作时 HEAD `93e2ad7`（工作树另有他人未提交改动）
- **官方侧**：`scripts/entrypoints/rerun_all_wells_corrected.py`（自述"论文/正式 8 井数字的官方入口"，`:1-17`）
- **生产侧**：`cemdisp/runners/<井>_tailpipe.py` 共 8 个，入口函数 `run_<井>_tailpipe_initial()`

---

## 0. 结论（先说结论）

**官方脚本产出的八井数字与 8 个生产 runner 不是同一口径，两者不可横向比较。**
逐行核到 **4 处**口径差（与计划所述处数一致），但**其中两处的实际后果与 brief 的暗示不同**：

| # | 差异 | 本文核实结论 |
|---|---|---|
| ① | 1D 三开关（管内段混浆口径） | **真差异**；量级已知但**本次未重测**（引用外部记录：7 井 ≤0.07 pp / hu2 −0.66 pp） |
| ② | 未传 `schedule=` | **真差异**；且**不是**求解器 docstring 所称的"只影响诊断"——它会经隔离液等效物性进入环空二维求解；**8/8 井进入该分支、6/8 井物性真的改变**（hu101/hu1 逐位相同，见 §2.2）；η 量级**未量化** |
| ③ | `total_t` 公式不同 | **公式确实不同，但 8 口井实测取值相同**（`min(...)` 在 8/8 井退化为同一个值）⇒ 数值上**不是**差异 |
| ④ | 额外 `CORRECTED_KW` | **名义 3 个开关，默认路径上净效应只有 1 项**（`enable_local_i3=True`）；另两项为空转 |

**是否把官方脚本收口到生产口径：待裁定**（Spec §4 R-A，建议默认"**收口**"，代价是会再次移动
G3 八井数字）。本任务**未执行**任何收口动作。

---

## 1. 差异 ① 1D 三开关（管内段混浆口径）

**官方侧**（`rerun_all_wells_corrected.py:73`，逐字）：

```python
    cr=CasingFlowSolver(enable_gravity=True).run(well,fluids,schedule)
```

只传 `enable_gravity`，**不传** `mixing_contact_time` / `plug_face_zero_mixing` / `has_plug`。

**生产侧**（8 个 runner 一致，示例 `cemdisp/runners/hu101_tailpipe.py:180-190`）：

```python
    casing_solver = CasingFlowSolver(
        enable_gravity=True,
        # T1 生产口径（2026-09-09 用户裁定）：双开关全开 + has_plug=True。
        mixing_contact_time=True,
        plug_face_zero_mixing=True,
        has_plug=True,
    )
    casing_result = casing_solver.run(well_spec, fluids, schedule)
```

同款构造在其余 7 个 runner 的行号：`hu102_tailpipe.py:289-299`、`hu103_tailpipe.py:292-302`、
`hu1_tailpipe.py:331-341`、`hu2_tailpipe.py:337-347`、`ht1_001_tailpipe.py:348-358`、
`ht1_003_tailpipe.py:293-303`、`ht1_004_tailpipe.py:311-321`。

**两侧都核实过默认值**（`cemdisp/transport1d/casing_flow.py:128-148`）：

```python
        has_plug: bool = False,              # :144
        mixing_contact_time: bool = False,   # :145
        plug_face_zero_mixing: bool = False, # :146
```

⇒ 生产 runner 的三个开关**全部偏离默认**，官方脚本则**全取默认**。这是真差异。

**量级（未在本任务重测）**：1D 三开关的已知效应 = **7 井 ≤0.07 pp、呼探1-002（hu2）−0.66 pp**。
出处为路线 B 落地记录（仓外 obsidian）：`论文撰写/论文整篇可写性终判与完成路线_2026-09-09.md:236`
（原文："关（现生产）vs 开（k_mix=5）差 7 井 ≤0.07pp / 呼探1-002 −0.66pp"）与同文 `:32`
（路线 B 落地为 T1 生产口径，`has_plug=True` ⇒ 混浆增强因子=1，"呼探1-002 受影响 −0.66pp"）。
⚠️ 该对数字是**三开关整组**（等价于"混浆增强关/开"）的效应，**未拆分到单个开关**；
`模型模块现状与环空段模拟调研_2026-09-25.md:176` 复述了同一对数字。**本任务未重跑，故不作二次确认。**

---

## 2. 差异 ② 未传 `schedule=`

**官方侧**（`rerun_all_wells_corrected.py:76-79`，逐字）：

```python
    kw=dict(total_t=tt,nz=NZ,enable_cfl_adaptive=cfl_on)
    if corrected:
        kw.update(**CORRECTED_KW)
    res=AnnulusD2DGASolver(**kw).run(well,fluids,inlet)
```

`:79` **不传** `schedule=`。

**生产侧**（8 个 runner 一致，示例 `hu101_tailpipe.py:113`）：

```python
    result = solver.run(well_spec, fluids, inlet_provider, schedule=schedule)
```

其余 7 处：`hu102:74`、`hu103:73`、`hu1:88`、`hu2:87`、`ht1_001:92`、`ht1_003:74`、`ht1_004:88`。

### 2.1 求解器 docstring 的口径描述**不准确**（本次读码核实）

`cemdisp/models2d/annulus_d2dga.py:2143-2146` 的 `run()` docstring 称：

```
            schedule: 泵注程序（可选，默认 None）。仅用于末尾 Tier0 诊断聚合：
                提供后 T0-6 停泵有限时间衰减诊断（shutdown_decay）可用；
                为 None 时诊断层优雅降级（记 notes "未提供 schedule"），
                不影响求解结果与既有调用方（向后兼容）。
```

但 `run()` 体内还有**第二处**消费（`annulus_d2dga.py:2160-2166`，逐字）：

```python
        _wash_spacer_fluids = [f for f in fluids
                               if f.role in {FluidRole.WASH, FluidRole.SPACER}]
        if len(_wash_spacer_fluids) > 1:
            _ws_weights = self._wash_spacer_volume_weights(_wash_spacer_fluids, schedule)
            spacer_fluid = self._composite_spacer_fluid(_wash_spacer_fluids, _ws_weights)
        # 诊断暴露：最近一次 run 实际使用的等效隔离液（rerun/报告脚本读取）
        self._active_spacer_fluid = spacer_fluid
```

数据链（全部读过）：

1. `_wash_spacer_volume_weights`（`:827-848`）：`schedule is None` → 返回 `None`；否则按泵注步的正向体积
   给出各 WASH/SPACER 的体积权重。
2. `_composite_spacer_fluid`（`:851` 起）：按体积分数**加权合成**等效隔离液的密度/黏度/屈服应力等
   （`volume_fractions=None` 时退化为等权）。
3. 该 `spacer_fluid` 作为环空二维的隔离液相参与 `_compute_props`（定义 `:1107`，调用 `:1917-1927`），
   消费点逐条为 `:1132`(μ)、`:1139`(ρ)、`:1150`(τy)、`:1174`+`:1175-1178`(n/K)，返回 `:1179`。
   **其产物在默认（流函数）路径上的实际去向只有两处**：
   - **I₃ 浮力通量**（`:2291-2314`，默认路径执行——上游 `enable_d2dga_i3_flux`/`enable_d2dga` 默认 `True`；
     `q_phi, q_xi = d2dga_buoyancy_flux(` 在 `:2311` 起、闭括号 `)` 在 `:2314`）：
     `rho`（含隔离液相贡献）**恒进入** `delta_rho`（`:2306` 或 `:2309`）；
     `mu` **只在 `enable_local_i3=False` 时**进入 `eta2`（`:2308`），
     开启局部化时 `eta2` 取水泥相黏度场 `_eta2`（`:2305`）
     ⇒ **隔离液 μ 是否进入 I₃ 取决于轴 ④（`CORRECTED_KW`）**；
   - **屈服门 `wall`**（`:2358-2361`，把 `mu`/`_tau_y` 喂给 `_yield_gate_wall`）：
     默认路径上 `wall` 只在 `enable_stream_yield_gate`（默认 `False`）时才传给速度场（`:1958`），
     故**默认路径上是仅诊断量**（进入 `wall_field` / `wall_snapshots`）。

   `n_mix`/`kappa_mix` 虽由 `_compute_props` 一并返回，但其动力学消费点在旧代数路径的流动度构造
   （`_mobility_base`/`_mobility_profile`，调用点 `:1481`/`:1992`/`:1999`），位于 `:1960` 早退之后
   ⇒ **默认路径上不参与求解**。

   ⚠️ **不要**写成"进入每步速度场构造"：默认路径的速度场函数 `_velocity_stream_function`
   （调用点 `:1956-1960`）**不接收 spacer 形参**（只收 `mud/lead/tail` 三个流体），
   按该函数验证会找不到该参数。（`spacer_fluid` 确实被传给 `_compute_velocity`——
   `:2229`（传参 `:2239`）与 `:2365`（传参 `:2375`）——但该函数在默认路径经
   `_velocity_stream_function` 返回，spacer 由此不再下行。）

### 2.2 8 口井**全部进入**该分支，其中 **6 口物性真的改变**（实测）

只跑 loader（不跑求解器）统计各井 WASH/SPACER 流体数：

| 井 | #WASH+SPACER | 名称 |
|---|---|---|
| hu101 | 2 | 平衡液(wash), 驱油隔离液(spacer) |
| hu102 | 3 | 平衡液(wash), 隔离液(spacer), 冲洗液(wash) |
| hu103 | 3 | 隔离液1(spacer), 隔离液2(spacer), 隔离液(实际)(spacer) |
| hu1 | 2 | 平衡液(先导泥浆)(wash), 驱油隔离液(spacer) |
| hu2 | 2 | 平衡液(wash), 隔离液(spacer) |
| ht1_001 | 2 | 平衡液(wash), 隔离液(spacer) |
| ht1_003 | 3 | 平衡液(wash), 隔离液1(spacer), 隔离液2(spacer) |
| ht1_004 | 3 | 先导浆(wash), 隔离液1(spacer), 隔离液2(spacer) |

⇒ `len(_wash_spacer_fluids) > 1` 对 **8/8** 井成立 ⇒ 官方脚本丢 `schedule=` 时，
**分支必被进入**。

**但"进入分支"不等于"物性被改变"（修复轮 1 订正）。** 直接对两口合成口径求值
（`_wash_spacer_volume_weights` + `_composite_spacer_fluid`，不跑求解器）实测：

| 井 | 各组分泵注体积 m³ | 归一化权重（生产/体积口径） | 官方等权 vs 生产体积权 | 实测物性变化 |
|---|---|---|---|---|
| hu101 | 25 / 25 | **0.5 / 0.5** | **逐位相同** | 无 |
| hu1 | 15 / 15 | **0.5 / 0.5** | **逐位相同** | 无 |
| hu102 | 15 / 20 / **0** | 0.4286 / 0.5714 / 0 | 不同 | ρ 2000.0 → 1985.714；YP 11.333 → 12.286 |
| hu103 | 17.5 / 17.5 / **0** | **0.5 / 0.5** / 0 | 不同 | ρ 1956.667 → 1975.0 |
| hu2 | 20 / 15 | 0.5714 / 0.4286 | 不同 | ρ 1950.0 → 1935.714；YP 1.5 → 1.714 |
| ht1_001 | 40 / 20 | 0.6667 / 0.3333 | 不同 | **ρ 1865.0 → 1826.667**；**YP 1.5 → 2.0** |
| ht1_003 | 28 / 16 / 10 | 0.5185 / 0.2963 / 0.1852 | 不同 | ρ 1916.667 → 1875.926；**YP 3.067 → 4.770** |
| ht1_004 | 25 / 16 / 10 | 0.4902 / 0.3137 / 0.1961 | 不同 | ρ 1816.667 → 1812.745；YP 9.867 → 9.839 |

⇒ **8/8 井进入分支，其中 6/8 井的等效隔离液物性真的改变**；**hu101 与 hu1 不改变**
（两组分泵注体积相等 25/25、15/15 ⇒ 等权与体积权归一化后都是恰好 0.5/0.5 ⇒ 合成逐位相同）。
另注意：hu103 的 `隔离液(实际)` 与 hu102 的 `冲洗液` **泵注体积为 0**，官方等权口径仍给它们
1/3 权重把它们平均进来（体积口径权重为 0 ⇒ 排除）——这是等权口径的**独立缺陷面**。

⚠️ 上一版本文写作"会改变**每口井**的等效隔离液物性"是**错的**（据 WASH/SPACER **计数**表
推出"物性改变"属推论越界；计数只能证明分支被进入）。严重性方向不变：
**差异 ② 不是"只影响诊断"**——它改变 6/8 井隔离液相物性，其中 `rho` 恒进入 I₃ 通量
（`:2306`/`:2309`），故会进入环空二维求解。

**对其余量级仍不量化**：上表是**合成物性**的实测差，**不是** η（位移）的量化；
要量化须同窗 A/B 重跑，不在本任务范围。

---

## 3. 差异 ③ `total_t` 公式

**官方侧**（`:52-53` 与 `:75`，逐字）：

```python
def _total_t(schedule):
    return sum(0. if s.rate_m3_min<=0 else s.volume_m3/s.rate_m3_min*60. for s in schedule.steps)
...
    tt=min(_total_t(schedule)+1200.,_stop_t(cr,fluids))
```

`_stop_t`（`:55-68`）优先返回 `float(cr.cement_end_time_s)`，缺失时才回退前缘扫描。

**生产侧**（`hu101_tailpipe.py:107` 与 `:202`；`annulus_stop_time_s` 定义在 `:54-89`，`:69`
返回 `float(casing_result.cement_end_time_s)`）：

```python
    total_t = total_t_s if total_t_s is not None else _schedule_total_time_s(schedule) + 20.0 * 60.0
...
        total_t_s=annulus_stop_time_value_s,   # = annulus_stop_time_s(casing_result, fluids)
```

⇒ 生产 `total_t` = `cement_end_time_s`，**没有 `min` 包装、没有 +1200s**。
8 个 runner 的 `annulus_stop_time_s` 逐行比对为同一逻辑（优先 `cement_end_time_s`，否则前缘扫描）；
行号：`hu102:148`、`hu103:151`、`hu1:166`、`hu2:165`、`ht1_001:171`、`ht1_003:152`、`ht1_004:171`。

**公式确实不同**，但**可由代码证明方向**：官方 `= min(泵注总时长+1200s, 尾浆入库时刻) ≤ 尾浆入库时刻
= 生产`。**是否取到实值差异，本次实测（只跑管内段 1D，不跑环空二维）**：

| 井 | 泵注总时长 s | 尾浆过鞋时刻 s | 泵注+1200 s | 官方 `total_t` = min(...) | 生产 `total_t` | 相等 |
|---|---|---|---|---|---|---|
| hu101 | 13742.7 | 13742.7 | 14942.7 | 13742.7 | 13742.7 | 是 |
| hu102 | 13629.5 | 11737.1 | 14829.5 | 11737.1 | 11737.1 | 是 |
| hu103 | 8135.8 | 8012.9 | 9335.8 | 8012.9 | 8012.9 | 是 |
| hu1 | 20944.4 | 20944.4 | 22144.4 | 20944.4 | 20944.4 | 是 |
| hu2 | 13525.0 | 13525.0 | 14725.0 | 13525.0 | 13525.0 | 是 |
| ht1_001 | 13219.7 | 13219.7 | 14419.7 | 13219.7 | 13219.7 | 是 |
| ht1_003 | 10340.9 | 10255.2 | 11540.9 | 10255.2 | 10255.2 | 是 |
| ht1_004 | 11409.6 | 11409.6 | 12609.6 | 11409.6 | 11409.6 | 是 |

⇒ **8/8 井 `min` 包装不激活，两侧 `total_t` 逐位相同**（"泵注+1200s"恒大于"尾浆过鞋时刻"）。
这与 `results/源模型口径重跑_2026-09-14/口径说明.md` 的自述一致（原文："旧 runner 停止时刻=尾浆入库时刻，
与上同口径——5 井已验证 min 包装不激活"）；本任务把该结论**扩到 8/8 井**。

⚠️ 该结论只覆盖 `total_t` **取值**；差异 ②（`schedule=`）仍**未量化**。

---

## 4. 差异 ④ 额外 `CORRECTED_KW`

**官方侧**（`:39-43` 定义、`:77-78` 施加）：

```python
CORRECTED_KW = dict(
    enable_yield_gate=True,    # M3: 屈服门槛
    enable_regime_split=True,  # M2: 局部流态修正（层流元 R=1，中性）
    enable_local_i3=True,      # I3: 浮力弥散通量局部化
)
...
    if corrected:
        kw.update(**CORRECTED_KW)
```

**生产侧**：`AnnulusD2DGASolver(total_t=total_t_s, nz=250)`（如 `hu101_tailpipe.py:108-111`），
不传任何 CORRECTED_KW ⇒ 取签名默认。

**默认值**（`annulus_d2dga.py:165` 起 `_SWITCH_DEFAULTS`）：

| 开关 | CORRECTED_KW | 默认（`:166/:172/:173/:185`） | 默认路径上是否真生效 |
|---|---|---|---|
| `enable_yield_gate` | `True` | `True`（`:172`） | ❌ **空转**：与默认同值 |
| `enable_regime_split` | `True` | `False`（`:173`） | ❌ **空转**：旧代数路径专属（见下） |
| `enable_local_i3` | `True` | `False`（`:185`） | ✅ **真生效**（见下） |

- `enable_regime_split` 的消费点在 `_compute_velocity`（`:2014  if self.enable_regime_split:`），
  但该函数在默认路径**提前返回**：`enable_stream_function` 默认 `True`（`:166`），
  函数体 `:1955-1960` 在新路径 `return` ⇒ `:2014` **不可达**。守卫文档同此口径
  （`:205-208`："旧代数路径专属开关：仅在 `enable_stream_function=False` 时被消费
  （`_compute_velocity` 在新路径早退…）"）。
- `enable_local_i3` 的消费点在 `:2304`，位于 I3 通量块 `:2291 if self.enable_d2dga_i3_flux and self.enable_d2dga:`
  之内；两个上游开关默认均为 `True`（`:183` `enable_d2dga` / `:184` `enable_d2dga_i3_flux`）⇒ 该块在默认路径执行，
  **`enable_local_i3=True` 真生效**。

⇒ **CORRECTED_KW 名义 3 开关，默认（`enable_stream_function=True` 流函数）路径上净效应 = 仅
`enable_local_i3=1` 一项。** 这是**独立读码结论**，与仓外 `模型模块现状与环空段模拟调研_2026-09-25.md:93`
（"按源码，`enable_regime_split` 默认路径不可达 ⇒ 该差应全部归 `enable_local_i3`"，该文自述"未重跑验证"）
一致。**量级未量化**（未重跑）。

---

## 5. 读了代码后**排除**的候选差异

1. **`enable_cfl_adaptive`**：官方 `:76` 显式传 `cfl_on`（`--cfl-mode` 默认 `on` ⇒ `True`）；
   生产不传 ⇒ 取签名默认 `enable_cfl_adaptive: bool = True`（`annulus_d2dga.py:399`）。**同值，不是差异。**
   （官方另有 `--cfl-mode fixed_dt` 模式，生产 runner 无对应开关——属**官方独有能力**，非口径差。）
2. **`nz` / `ny`**：官方 `:32 NZ=250` + `:76 nz=NZ`；生产 `nz=250`。`ny` 两侧都不传 ⇒ 默认 40。**同值。**
3. **`split_cement_phases`**：官方 `:74` 传 `split_cement_phases=True`；生产同（如 `hu101:196`）。**同值。**
4. **官方 `:73` 与 `:74` 各建一个 `CasingFlowSolver` 实例**（生产复用同一实例）：
   两处都是默认构造，**当前不构成数值差异**；但若将来在 `:73` 加开关而漏改 `:74`，两处会分裂——属**维护风险**
   （记为观察，不作断言）。
5. **入口桥只在停泵回退分支咨询 solver**：`boundary_bridge.py:179-189` 的 legacy provider 仅在
   `pipe_exit_state.flow_rate_m3_s < 1e-9` 时调用 `casing_solver.pipe_exit_state_at(...)`。
   该分支是否受 1D 三开关影响，**本文未验证、不声称**。

---

## 6. 量化一览（哪些量化了、哪些没有）

| 差异 | 量化？ | 数值 | 依据 |
|---|---|---|---|
| ① 1D 三开关 | **外部已量化，本次未重测** | 7 井 ≤0.07 pp；hu2 −0.66 pp | 仓外 obsidian `论文整篇可写性终判与完成路线_2026-09-09.md:32/:236`；本任务**未重跑** |
| ② `schedule=` | ❌ 未量化（只到"合成物性"层） | 8/8 井进入分支、**6/8 井隔离液物性改变**（最大 ρ 差 38.3 kg/m³ @ht1_001、YP 差 1.70 Pa @ht1_003）；η 未量化 | §2（直接求值 `_composite_spacer_fluid` 两种权重口径） |
| ③ `total_t` 公式 | ✅ **本次实测为 0** | 8/8 井两侧取值逐位相同 | §3 表（只跑管内段 1D 的实测） |
| ④ `CORRECTED_KW` | ❌ 未量化（但**开关活性已定**） | 净效应 1 项（`enable_local_i3`） | §4（读码 + 默认值表） |

---

## 7. 官方口径快照的可追溯缺口（本次读码新增发现）

官方脚本用 `write_adopted_config()`（`:141-172`）导出口径快照。实读 G3 输出目录的快照
（`results/源模型口径重跑_2026-09-14/adopted_config.json`）：

- `solver` = `{enable_yield_gate: True, enable_regime_split: True, enable_local_i3: True,
  enable_cfl_adaptive_default: True}`、`grid = {nz: 250, ny: 40}`、`cfl = {mode: on}`、
  `boundary = {split_cement_phases: True, casing_enable_gravity: True}`、
  `git_commit = 98ab9245…`、`generated_at_utc = 2026-09-15T04:56:04Z`。

⚠️ **`boundary` 块只记录 `split_cement_phases` 与 `casing_enable_gravity` 两项**，而官方脚本与生产
runner 在 `boundary` 上分歧的正是 **1D 三开关**（`mixing_contact_time` / `plug_face_zero_mixing` /
`has_plug`）——它**不在**快照里；`schedule` 与 `total_t`（轴 ②③）也都不在快照里。
⇒ **轴 ①②③ 在快照中缺席**：仅凭 `adopted_config.json` 无法判定一套八井数字用的是哪一套管内段
1D 口径、是否传了 `schedule`、`total_t` 取哪个公式。**建议把这三轴补进快照。**

⚠️ 上一版本文写作"这三者恰是官方脚本与生产 runner 分歧的**全部轴** ⇒ 快照无法自证它自己不是
生产口径"——**两句都错**（修复轮 1 订正）：

- 分歧有 **4 轴**（见 §1-§4），不是 3 轴；
- **轴 ④（`CORRECTED_KW`）恰恰是被记录了的**：`:161-164` 的 `solver` 块逐字展开 `**CORRECTED_KW`
  （本段上方已引用该 JSON）。而且它**正判别**官方与生产——记录值含 `enable_local_i3: true` 与
  `enable_regime_split: true`，而**没有任何生产 runner 设置这两项**：
  8 个 `*_tailpipe.py` 只传 `total_t`/`nz`（`hu101_tailpipe.py:108-111` 逐字），
  全文搜索 `cemdisp/runners/` 只有 `ht1_004_ablation.py`（消融入口，非常规生产 runner）设置它们。
  ⇒ 快照**能**把官方运行与生产口径区分开，分的正是轴 ④。

✅ 可辩护的结论收窄为：**快照记录了轴 ④，缺 轴 ①②③** ⇒ 这不是"快照无用"，而是"**快照不全**"；
若要保留双口径，应把这三轴写进快照（这仍支持"保留双口径需补记口径"的建议）。

（附带一条已核实的更细缺陷）`write_adopted_config` 在 `main()` 末尾**无条件**调用（`:239`），
而其 `solver` 块是**硬编码的 `**CORRECTED_KW`**（`:161-164`）。若调用时**未**加 `--skip-baseline`，
`汇总.csv` 会同时含 `corrected=False` 的基线行，而这些行**没有**用 `CORRECTED_KW`；
同一份快照却声称 `CORRECTED_KW` ⇒ **会高估其自身基线行的口径**。
G3 那次的 `汇总.csv` 里 `corrected` 列全为 `True`（确实用了 `--skip-baseline`）⇒ 该缺陷在 G3 上
**未触发**，但下次不带该开关跑就会触发。

另核实（可用于复核）：
- 官方脚本默认 `--out-dir`（`results/全井修正前后`，`:31`）**在工作树中不存在**；
  G3 实测输出在 `results/源模型口径重跑_2026-09-14/`（用 `--out-dir` 指定）。
- 该目录 `汇总.csv` 的表头 = `cement_occ,cfl_mode,channeling,corrected,elapsed_s,eta_E,eta_N,
  instability,mixing,wall_frac,well`，与官方 `run_one` 的返回键（`:81-87`）**逐字一致** ⇒ 出处可核；
  其 `corrected` 列全为 `True`、`cfl_mode` 全为 `on` ⇒ 与 `--skip-baseline` + 默认 CFL 一致。
- 该 `汇总.csv` 的 `eta_E/eta_N` 与 Task 1 基线 md（`docs/superpowers/plans/baseline-2026-09-25.md`）
  逐井一致（如 hu101 0.9710753804210159 / 0.82943572033394 ⇒ md 的 0.971075 / 0.829436）
  ⇒ 基线行的**数值来源 = 官方脚本口径**（其非生产性质正是本报告的结论）。

---

## 8. 复核方法（可复现）

```bash
# 差异 ①②④：读两侧代码即可
sed -n '39,43p;52,53p;70,87p' scripts/entrypoints/rerun_all_wells_corrected.py
sed -n '107,113p;180,203p'     cemdisp/runners/hu101_tailpipe.py
sed -n '165,187p;1955,1962p;2014p;2129,2146p;2160,2166p;2291,2314p;2358,2361p' cemdisp/models2d/annulus_d2dga.py

# 差异 ② 的"8/8 进入分支、6/8 物性改变"：只跑 loader，对 _composite_spacer_fluid 两种权重直接求值
#   eq  = AnnulusD2DGASolver._composite_spacer_fluid(ws, None)                     # schedule=None 口径
#   vol = AnnulusD2DGASolver._composite_spacer_fluid(ws, S._wash_spacer_volume_weights(ws, schedule))
#   比较 density_kg_m3 / yield_stress_pa / plastic_viscosity_pa_s / power_law_n / consistency_k（见 §2.2 表）
# 差异 ③ 的 total_t 实测：只跑管内段 1D 取 cement_end_time_s，与 schedule 总时长比（见 §3 表）
```

三处实测均为**只读**（不写任何文件、不改任何口径）；脚本片段见本任务的
`.superpowers/sdd/2026-09-25-internal-consistency-hardening/task-7-report.md`。
