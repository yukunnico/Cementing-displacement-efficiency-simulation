# 续跑交接书 —— 源模型口径重构 + 仓库清理（Task 5 起）

> **新窗口的第一条指令**：读本文件，然后按文中「执行流程」继续跑 Task 5 → Task 14。
> 生成时间：2026-09-14（上一窗口因上下文将满而交接）
> 本文件是**唯一入口**：状态、路径、裁定、命令、陷阱全在这里；其余文件按文中路径按需读。

---

# 0. 一句话现状

分支 `refactor/d2dga-source-fidelity`，**HEAD = `e762485`**，全量测试 **405 passed**。
**已完成 Task 0–4（仓库整理 + 浮力口径统一 + F² 定标），Task 5 刚派发即被主动停止、未留下任何改动。**
剩余 **Task 5–14**。

---

# 1. 关键路径（全部相对仓库根 `D:\users\desktop\research\控压固井项目\cement model`）

| 用途 | 路径 |
|---|---|
| **实施计划**（14 个 Task 的权威文本） | `docs/superpowers/plans/2026-09-14-d2dga-source-fidelity-and-repo-cleanup.md` |
| **基线**（重构前测试数 / 8 井 η_E / 基准算例 10 例 / CBL 靶值 / b 表） | `docs/superpowers/plans/baseline-2026-09-14.md` |
| **SDD 账本**（进度 + 全部裁定 + 预检扫描 + deferred minors） | `.superpowers/sdd/2026-09-14-d2dga-source-fidelity-and-repo-cleanup/progress.md` |
| **任务简报**（`task-N-brief.md`，派发用） | 同上目录 |
| **实现者报告**（`task-N-report.md`） | 同上目录 |
| **审查包**（`review-<base>..<head>.diff`） | 同上目录 |
| **两份源文献调研报告**（方法学依据） | obsidian `固井顶替效率改进/半环空设定与源文献适用域审查_2026-09-14.md`、`模型全面调研_早期版本对照与敏感性根因_2026-09-14.md` |

---

# 2. 已完成任务与提交（供你校验）

| Task | 内容 | 提交 |
|---|---|---|
| 0 | 分支 + 基线冻结 | `3e187bd` |
| 1 | scripts/ 归档（76 → entrypoints14/probes44/reruns7/plots7/lib2 + 2 个 `.DEPRECATED.py` + README） | `0d38927`→`78112ac`→`80f597b`（修复轮） |
| 2 | tests/ 分层（contract 22 / history 9）+ results/ 归档（53→9 保留） | `681dcd5`→`33b7750`（修复轮） |
| 3 | `buoyancy.py` 统一浮力口径 | `7bda971`→`be0447e`（修复轮） |
| 4 | `F²` 按 Z&F22 (2.6) 定标 | `f10a8b9`→`1ec6260`（修复轮） |
| — | controller 的文档/计划提交 | `a0519ab`、`3b20adb`、`9cf0c78`、`99c524d`、`e762485`、`6e6cd97` |

**Task 4 的实测效果**（供你判断后续位移是否合理）：`F²` 1.0→2.0e-3~5.5e-3；方位浮力 `f_φ(max)` 4.02e-3→**0.74~2.01**；`|1−buoyancy_shape|max` 8.3e-5→**4.2e-2**；**呼101 `η_E` 0.6805414→0.6938376（+1.33pp）**。回归锚：`f2` 强制回 1.0 ⇒ 与冻结基线**逐位相等**。

---

# 3. 执行流程（严格照做）

本计划的执行方式是 **superpowers:subagent-driven-development**。每个 Task 一轮循环：

```
1. BASE=$(git rev-parse --short HEAD)                    # 记录派发前的 HEAD
2. 读 plan 里该 Task 的文本（或直接用已生成好的 task-N-brief.md）
3. 派 implementer subagent（general-purpose，model 见下）
   —— 派发提示词必须包含：①该 Task 的一行定位 ②brief 路径 ③你为它消解的歧义
      ④全局约束 ⑤报告文件路径与回报契约 ⑥"你不得再派子 agent"
4. implementer 回报（DONE / DONE_WITH_CONCERNS / BLOCKED / NEEDS_CONTEXT）
5. 生成审查包并派 task reviewer（spec 合规 + 质量双裁）：
   bash "C:/Users/30525/.claude/plugins/cache/claude-plugins-official/superpowers/6.3.0/skills/subagent-driven-development/scripts/review-package" \
        "docs/superpowers/plans/2026-09-14-d2dga-source-fidelity-and-repo-cleanup.md" <BASE> HEAD
6. 有 Critical/Important → fix loop（≤5 轮）：resume 原 implementer（SendMessage）
   → 生成 scoped 复审包（FIX_BASE=上一轮审查看到的 HEAD）→ 派 re-reviewer
7. 全部 addressed → 记账 `Task N: complete (commits a..b, review clean)` → 下一个 Task
8. 全部完成后：整分支终审（最强模型）+ 按 superpowers:finishing-a-development-branch 收尾
```

**模型选择**：机械/照抄类 → sonnet；需设计判断（Task 5/8/9/11）→ **opus**；reviewer → sonnet；终审 → opus。

**⚠️ 必须串行**：`cemdisp/models2d/annulus_d2dga.py` 被 Task 3/4/5/7/9/10/11/14 触及 → **一个 Task 一次，不得并行派 implementer**。

**审查包太大时**（如 Task 2 的 787 文件）：自己用 `git diff --name-status` + `git diff -- <实质目录>` 拼一个聚焦包，把批量移动折叠成计数。

---

# 4. 剩余任务详细步骤

> 每个 Task 的**完整原文**在计划文件里；`task-N-brief.md` 已全部预生成。下面只补「计划没写、但我已经知道」的东西。

## Task 5：领先阶轴向浮力数 `b` 接进动力学 ★最关键

**目标**：Z&F22 (4.14) 的 `I₂·(Δρ/(H·r_a))` 浮力项进 `pref`。目前 `b` 只写进 summary，从不参与求解。

**我已预置的三条**：
- **R26**：`_mobility_profile` 签名**必须含 `f2`** 并实际使用（否则 Task 4 的 F² 定标二次静默消失）。**验收要有测试断言输出随 `f2` 变化。**
- **R2**：`K_AXIAL` 是 provisional —— 实现者须**自己完成量纲推导**并写进 docstring；取值可非 1，但必须在报告里写明推导与量级后果。
- **R27**：`buoyancy_shape` **必须始终 > 0**；须报告「八井 × 各排量段」的 min/max 实测。按八井 `b ∈ [−2.9, +20.7]`、`i2/i1 ~ O(0.1)`，`K_AXIAL=1` **可能越界**——若越界，重新推导，**不要硬压**。

**顺手修**：`tests/contract/test_buoyancy.py` 一处注释仍写 "物理 F² 取 O(10⁻²)"（应为 O(10⁻³)）；Task 4 的调用点守卫用严格集合相等，若新增第三个 `_buoyancy_force_vector` 调用点会误红 → 同步更新但**不得削弱其判别力**。

**预期**：会显著改变 `η_E`。**如实报告位移，不得为数字好看调 `K_AXIAL`。**

---

## Task 6：新建 `cemdisp/models2d/two_layer.py`

把 `d2dga_flux.py` 的 I₁/I₂/I₃ **迁移**过去（`d2dga_flux.py` 保留为 re-export 薄层，避免破坏既有 import），并新增 `isotropic_flux_q0`（Z&F22 (4.28)）与 `layer_thickness_fraction`。
**测试**：用 `m ∈ {0.2,0.5,1,2,5}` × `c̄ ∈ {0.25,0.5,0.75}` 全网格验证 (4.26) 与独立解析复算一致到 `1e-12`。
**⚠️ 硬约束**：`I₃` **必须**用 Z&F22 (4.26) —— B&F25 (2.27) 的印刷式（括号 `3(1−c̄²)`）已核实与 Z&F22/Z&F23/解析解三者皆不自洽，**不得替换**。

---

## Task 7：删自创弥散，改用分层通量闭合

删 `_smooth_dispersion`（`:658-687` 附近）及其两处调用；`dispersion_axial/azimuthal/dt_ref/dt_scale` 保留形参但置 `None` + `warnings.warn`（避免破坏 8 个 runner）。
**文献依据**：Z&F22 p.11 "we have no diffusive terms"；弥散由 `q₀` + `I₃` 承载。
**同时**：更正模块 docstring 里"自创弥散=论文口径"的错误表述。
**测试**：断言 `run()` 的源码里不再出现 `_smooth_dispersion`。

---

## Task 8：新建 `cemdisp/models2d/stream_function.py` ★设计判断

解 Z&F22 (4.22) `∇a·[(r_a/2I₁)∇aΨ + b] = 0`。
**关键简化（已验证）**：I₁/I₂ 只依赖 `c̄` 与 `m`（本步已知）⇒ 方程对 Ψ **线性**，每步是**变系数线性 Poisson**，`scipy.sparse.linalg.spsolve` 直接解（40×250=10000 未知数，实测 <0.15 s）。**论文的增广拉格朗日只在 HB 非线性 `S(|∇Ψ|)` 时才需要**——不要照抄那套。
BC：`Ψ(φ=0)=0`、`Ψ(φ=1)=1`（单位通量，调用方按 Q 缩放）；`v̄ = −∂ξΨ/(2H)`、`w̄ = ∂φΨ/(2 r_a H)`（Z&F22 (2.2)）。
**测试**：BC 与单调性；**每列通量与深度无关**；`b>0` 时窄边流速份额上升。

---

## Task 9：求解器接入 Ψ，替换代数流动度 ★风险最高

加开关 `enable_stream_function`（默认 True；False 走旧代数路径，保证可回退）。
**废止 `w_d2dga = w * f_amp`**（`annulus_d2dga.py:1168-1169` 附近，行号可能已漂移）——`f_amp` 只保留在通量层（由 Task 6 的 `q₀` 承担）。
**R7（预检裁定）**：补一条回归护栏测试，断言 `enable_stream_function=False` 时输出与重构前基线一致。
**验收**：基准算例 `mass_conservation_error` 从 0.33–0.48 降到 **< 0.05**；`tests/history/` 预期大面积变红（属预期，逐条标注）。

---

## Task 10：消 `e_clip` 截断

`annulus_d2dga.py:395-401` 的 `e = clip(1−standoff, 0.05, 0.55)` → `clip(1−standoff, 1e-6, 1−1e-6)`（文献 `e ∈ [0,1)`，Pelipenko04 (2.1)）。
`e_clip_max`/`e_clip_measured_max`/`enable_e_clip_ruling` 保留形参 + `warnings.warn`。
**基准算例里 `e_clip_max=1.0` 是显式传的**，改后仍应工作。
**验收**：`standoff=0.35` 与 `0.45` 给出**不同**输出（现状逐位相同 = 死区）。

### ⚠️ Task 10 必须裁定的遗留问题：hu101 的 `standoff = 0.80`

工作区当前（已随 Task 2 提交）`cemdisp/data/loaders/hu101_loader.py` 的 `_ASSUMED_STANDOFF` 是**全井常数 0.80**，来自**用户 2026-09-11 的既有指令**（"直接改他的居中度为0.8"）。但现场记录是"悬挂器**坐挂失败、座底固井**，居中度**下降**"——两者方向相反，且代码注释自认"**不得作为模型验证数字引用**"。
**用户在 2026-09-14 的指令是"修复要按照现场的实际情况来"**。
**→ 在 Task 10 里做一次显式裁定**：是保留 0.80（尊重用户既有指令）还是回退到 LEGACY 的 0.38–0.48（贴现场记录）。**两条都要记账**，并在最终交付时明确告诉用户这个选择及其后果。
（另：`tests/history/test_hu101_loader_standoff.py` 断言 0.80，`hu101_loader.py:14` 模块 docstring 仍写 0.38–0.48 —— 若回退，两处同步改。）

---

## Task 11：屈服门二值 → 连续

`annulus_d2dga.py:624-655` 附近 `wall_new = np.where(immobile, 1.0, 0.0)` → `wall = clip(1 − τw/(f_safety·τy), 0, 1)`（Pelipenko04 (2.6)-(2.8) 停流区判据的连续近似）。
**R3（预检裁定）**：**必须保留** `ref_mask` 与 `col_freeze` 两个不变量（参考元永不冻结 / 整列无流动则冻结），否则 `w=0` 的列会 `0/0`。

---

## Task 12：基准算例复现（数值正确性验收）

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m cemdisp.runners.zhang2022_benchmark \
  --out-dir results/基准算例对照_2026-09-14
```
⚠️ **必须显式传 `--out-dir`**（R18：`DEFAULT_OUTPUT_DIR` 指向权威目录 `results/基准算例对照_2026-09-10`，不传会覆写）。

**三条验收门槛**（对比 `baseline-2026-09-14.md` §4）：

| 门槛 | 重构前 | 目标 |
|---|---|---|
| `\|Δη_E\| ≤ 0.05` 的算例数 | **2/10** | **≥ 8/10** |
| `mass_conservation_error` | 0.33–0.48 | **< 0.05** |
| e 从 0.8→0.1 的 η_E 响应跨度 | **1.3pp** | **≥ 20pp**（论文 31pp） |

未达标 → 写归因报告、回 Task 8/9 迭代（最多 2 轮）；仍不达标如实报告并停在可回退点。

---

## Task 13：8 井重跑 + 现场对照

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/entrypoints/rerun_all_wells_corrected.py \
  --out-dir results/源模型口径重跑_2026-09-14
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/entrypoints/cbl_window_comparison.py \
  --results results/源模型口径重跑_2026-09-14
```
**默认不替换权威目录 `results/<井名>_1D2D耦合模型/`**——先出新目录，**经用户确认后再替换**。
**现场 CBL 靶值**（用 **CBL 评价窗口口径**，不是全井口径）：呼101 62.77% / 呼102 66.65% / 呼103 12.06% / 呼1-003 78.7% / 呼1-004 0.3%。目标 Spearman **ρ ≥ +0.7**（重构前 8 井全井口径 ρ≈−0.30，CBL 窗口口径 ρ≈+0.40）。
**出对照表**：重构前 η_E/η_N vs 重构后 vs CBL 靶值 vs 窗口 η_E vs 偏差 pp。

---

## Task 14：论文边界声明与配置表同步

写 `docs/源模型口径与适用域声明.md`，六条（每条对应文献）：
1. 半环空 = 人为 Ψ 对称条件（Z&F22 p.25/p.32），抹掉 3-D 方位非对称
2. 窄边=低边假设（Pelipenko04 p.4）
3. **源模型严格适用域 = 竖直井 + 牛顿**（Z&F23 p.34）；本项目 HB/斜井为**扩展应用，未获外部验证**
4. 8 井 `b ∈ [−2.9, +20.7]`，**无一口达 `b ≳ 80`**；呼探1/呼探1-002 密度倒置（文献警告区）
5. 弥散为 `q₀/I₃` 解析闭包，**无自由弥散参数**
6. `e ∈ [0,1)`，无截断

外加「式号 → 代码锚点」对照表。

---

# 5. 全部裁定（R1–R27）—— 约束后续所有工作

| # | 裁定 | 代价若错 |
|---|---|---|
| R1 | 用 in-place 分支（非 worktree），因工作区有用户大量未提交改动 | 失去隔离，由"同文件串行"覆盖 |
| R2 | Task 5 的 `K_AXIAL` 为 provisional，须写推导 | 轴向浮力幅值偏差，T12/13 需重标定 |
| R3 | Task 11 连续化须保留 `ref_mask` 与 `col_freeze` | 静止列 NaN 或错误解冻 |
| R4 | scripts 归档以派发清单为准 | 个别脚本留原地 |
| R5 | `tests/__init__.py` 与锚 json 留 `tests/` 根 | pytest 发现路径变化 |
| R6 | results 归档用**显式白名单**、禁通配符 | 误伤权威目录 |
| R7 | Task 9 补 `enable_stream_function=False` 回归护栏 | 无回退验证 |
| **R8** | **禁用 `scripts/entrypoints/smoke_all_wells.py`**（覆写权威目录）；冒烟改用只读 pytest | 覆写权威结果 |
| R9 | 接受 Task 1 的两处 1 行越界改动 | 无 |
| R10 | Task 1 的 265→284 行重写是**用户 2026-09-11 的既有改动**，不回退 | 无 |
| R11 | `!scripts/lib/` 修 `.gitignore` 的 `lib/` 误伤 | 新模块被静默忽略 |
| R12 | 修 R8 的错误路径 + Task 1 Step 6 冒烟命令 | 后续跑错路径 |
| R13/R14/R15 | 落地两个未跟踪/未提交的 `cemdisp/` 文件恢复分支自洽 | 分支不可用 |
| R16 | 31 条既有 `results/` 删除**披露不恢复** | archive 完整性略降 |
| R17 | **不把 `参考文档/`（716MB）入库**；声明"干净检出限制" | 分支无法仅凭 git 全绿 |
| R18 | `zhang2022_benchmark.py` 默认 out-dir 指向权威目录 → **必须显式传 `--out-dir`** | 覆写权威基准结果 |
| R19 | `_compute_buoyancy_number` 保留为薄委托（API 兼容） | 终审可能判死代码 |
| R20 | `half_gap_m` 用 **`float(np.mean(geom["H"]))`**；**禁止** `0.5*mean((hole−od)/1000)`（2× 错误） | b 差 4×，T4/5 全偏 |
| R21 | `froude_squared` 的 docstring 曾写成 (2.6) 的**倒数**，已订正（函数体本就对） | 反向数万倍 |
| R22 | controller 先前两份文档的 b 值同样偏大 4×，已更正 | 引用错值 |
| R23 | Task 4 断言改正比关系（非 `> 50×` 阈值） | 断言随 F² 失效 |
| R24 | 3 条 history 红灯属 API 变更，机械补 `f2=1.0`，**不得删改断言** | 掩盖真实回归 |
| R25 | δ₀ 引文改为论文的 `δ = d̂/r̂ₐ*`（**不要再乘除 π**） | 下一任务"修正"加 π，η_E 再动 |
| **R26** | **Task 5 的 `_mobility_profile` 必须显式透传 `f2`** | F² 定标二次静默消失 |
| **R27** | **Task 5 的 `buoyancy_shape` 必须始终 > 0；须报 min/max 实测** | 非物理负流动度 |

---

# 6. 已记录的 deferred minors（终审时统一 triage）

- `scripts/entrypoints/smoke_all_wells.py` 的 stdout-JSON 解析缺陷 → 汇总指标恒为 N/A（既有）
- `cemdisp/` 内 8 个 runner 的注释 + `PACKAGE_REFERENCE.md` §7 仍指向旧 `scripts/` 路径
- `cemdisp/transport1d/casing_flow.py:519/525` 有 `or 0.01` 静默回退（1D 路径）
- `annulus_d2dga.py` 里 `_compute_buoyancy_number` 死代码 + `gap_m`(全间隙)/`half_gap_m`(d̂) 双命名并存
- `annulus_d2dga.py:1432` 注释笔误（"H=2b 一半"）
- Task 4 M4 守卫用严格集合相等，新增调用点会误红
- `tests/contract/test_buoyancy.py` 注释仍写 O(10⁻²)
- Task 4 退化几何分支 ×2 在真实井不可达、无实测佐证

---

# 7. 已知陷阱（给后续所有 agent）

1. **`d̂` 口径**：`geom["H"]` 就是 Z&F22 的半间隙 `(hole−od)/4`；`geom["b"] = 2H`。**写成 `(hole−od)/2` 会差 2 倍、`b` 差 4 倍**——本计划已因它返工一次（Task 3），controller 自己也犯过一次。
2. **`F²` 不是它的倒数**：`F² = τ̂₀/(ρ̂₁ĝδ₀r̂ₐ*)`，`τ̂₀ = μ̂₁ŵ₀/d̂`。
3. **`δ` 不要再乘除 π**：论文 `δ = d̂/r̂ₐ*`（原文写作 `d̂/(πr̂ₐ*) = δ/π ≪ 1`）。
4. **`I₃` 用 Z&F22 (4.26)**，不是 B&F25 (2.27)。
5. **`f_φ/f_ξ` 的角标**：论文的 `f` 角标与它作用的动量方程**交叉**（(4.11)(4.12)），代码已换成物理命名，**是等价的**，不要去"修"。
6. **别用宽泛 `git add`**：Task 1/2 都因此把用户既有改动裹进提交。
7. **`smoke_all_wells.py` / `zhang2022_benchmark.py` 默认 out-dir 都会写权威目录**。

---

# 8. 全部跑完要交给用户什么

1. **基准算例复现表**（10 例：模型 η_E vs 论文 η_E、`mass_conservation_error`、e 响应跨度）与三条门槛的达标情况
2. **8 井重跑对照表**（重构前 vs 重构后 η_E/η_N vs CBL 靶值 vs 窗口 η_E vs 偏差 pp）+ Spearman ρ
3. **「我做的裁定」完整清单**（R1–R27 + 后续新增，每条附"代价若错"）—— 这是用户唯一能看到你替他做决定的窗口，**不得只列一部分**
4. 未修复的 deferred minors 清单与 triage 建议
5. 收尾按 **superpowers:finishing-a-development-branch**

---

# 9. 新窗口的开场指令（可直接粘贴）

> 继续跑 `docs/superpowers/plans/2026-09-14-CONTINUATION-HANDOFF.md` 里交接的计划：先读它，再读 `docs/superpowers/plans/2026-09-14-d2dga-source-fidelity-and-repo-cleanup.md` 与 `.superpowers/sdd/2026-09-14-d2dga-source-fidelity-and-repo-cleanup/progress.md`，确认 HEAD 是 `e762485`、分支是 `refactor/d2dga-source-fidelity`、全量测试 405 passed，然后按交接书的「执行流程」从 **Task 5** 开始逐 Task 执行（用 superpowers:subagent-driven-development），不要中途停下来问我，全部跑完把结果给我。
