# tests/ 分层说明

本目录的测试按**"期望在源模型口径重构（2026-09-14 起）中如何反应"**分为两层。
判据不是文件名、不是测试粒度，而是**这条断言在重构后应该保持绿还是应该变红**。

```
tests/
├── __init__.py                       # 包标记（测试以 tests.* 导入）
├── _baseline_shoe_timeline_hu103.json  # 共享基线锚：hu103 shoe_timeline 逐事件冻结值
├── contract/                         # 契约测试：守护当前正确性
├── history/                          # 历史锁定测试：锁定重构前行为
└── README.md
```

`__init__.py` 与 `_baseline_shoe_timeline_hu103.json` 不属任何一层：前者是包标记，
后者是被 `history/test_casing_mixing_contact_time.py` 消费的共享锚，留在根。

---

## `tests/contract/` —— 契约测试（守护当前正确性）

**判据**：断言描述的是"模型**应该**有的行为"——数据契约、物理闭包、
外部接口契约、已经裁决过的口径。它不记录某次实现的输出数值，
而是记录**规格**。

代表性判据：

- 数据结构/加载契约：`WellSpec`/`FluidSpec`/`PumpingSchedule` 的字段语义与校验、
  loader 的现场值口径、`provenance` 出处标注；
- 已裁决的物理闭包：`d2dga_flux` 的 I₁/I₂ 闭式、`regime_closure`/`regime_classifiers`
  的流态判据、`muskat_regime`、`axial_dispersion`、`yield_deadzone`；
- 外部接口契约：`export_depth_time_shares`（导出文件清单与列口径）、
  `sync_cards`（同步卡片）、`tier0_schedule_wiring`（诊断接线）、
  `boundary_bridge`/`outlet_boundary`（边界耦合）、
  `evaluation_windows`（评价窗定义）、`displacement_metrics`/`m0_metrics`（指标定义）、
  `shoe_timeline`/`casing_flow`（套管 1D 时序与流动）、
  `six_well_integration`（六井端到端只读集成）。

**期望**：在源模型口径重构全过程中**保持全绿**。
重构若使某个 contract 测试变红，只有两种可能——要么重构引入了回归，
要么这条"契约"其实应该被重新裁决并**显式改写**。不允许静默变红或直接删除。

## `tests/history/` —— 历史锁定测试（改模型预期变红）

**判据**：断言里钉着**某个具体实现路径的输出数值或开关行为**——
逐位（bitwise）对照锚文件、某次修复前后的差分、"默认关 = 逐位复现旧路径"、
某篇论文基准算例在当前实现下的具体取值。

代表性判据：

- 逐位锚：`casing_mixing_contact_time` 的"双开关默认关 → shoe_timeline 逐事件
  等于 `_baseline_shoe_timeline_hu103.json`"；
- 某次修复的差分固化：`three_fixes_20260906`、`pipe_capacity_chain_fix`、
  `plug_semantics_restart`、`m1_dispersion`、`casing_mixing_buoyancy`；
- 求解器级现状快照：`improved_d2dga_annulus`；
- 实测口径特例：`hu101_loader_standoff`；
- 论文基准算例的当前输出：`zhang2022_benchmark`（Table 1/2/3 对照）。

每个 `history/` 文件头都有 `# STATUS: history` 标注块，说明同一件事。

**期望**：重构（浮力口径统一、F² 定标、流函数速度场替换、废止 `f_amp` 速度乘子、
删自创弥散、消 `e_clip`、屈服门连续化）**会让其中一部分断言变红**——
这是设计意图，**红灯不等同回归失败**。

处理规则：

1. 变红时必须**逐条**对照
   `docs/superpowers/plans/2026-09-14-d2dga-source-fidelity-and-repo-cleanup.md`
   判断：是"预期变红"（旧口径被有意替换），还是回归；
2. 属"预期变红"的，**改写为新口径下的契约或新锚**，并把旧断言注释保留在原处
   说明它锁的是哪个已废弃口径；不得直接删除；
3. 属回归的，按回归处理，不得借"history 层可以红"蒙混过关；
4. 全量回归不得因为 `history/` 变红而整体判定失败——但必须逐条给出结论。

---

## 运行

```bash
# 全量（基线：396 passed）
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest tests/ -q

# 只跑守护层（重构期间应保持全绿）
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest tests/contract/ -q

# 只跑历史层（看哪些旧口径被动了）
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest tests/history/ -q
```

只读冒烟（**不要**用 `scripts/entrypoints/smoke_all_wells.py`——它会写入权威结果目录
`results/<井名>_1D2D耦合模型/`）：

```bash
PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python -m pytest tests/contract/test_six_well_integration.py -q
```

## 新增测试放哪一层？

问自己一句：**"如果源模型口径改了，这条断言应该跟着改吗？"**

- 应该跟着改 → `history/`（它锁的是实现，不是规格），并在文件头加 `STATUS` 标注；
- 不应该改、改了就是 bug → `contract/`。

新写契约测试时优先放 `contract/`；只有在固化"某次修复的前后差分"或
"某篇论文在当前实现下的具体数值"时才放 `history/`。
