# Phase 3 设计规格：压力全链 P-2→P-4（水力学内核 + 产品层）（2026-10-07）

> 上游 = 总纲《详细执行计划_真温压响应_2026-10-06》§3 Phase 3 + 续作计划 §7；
> 靶权威 = `results/_probe_matlab靶_2026-10-07/3.0_靶复现判定.md`（用户裁定 2026-10-07：
> **以本次沙箱运行为权威靶**）；结构权威 = 同目录 `P30_移植地图.md`（S1–S25 逐段公式/单位/常数）。
> **依据纪律（§0.2）**：本阶段一切数值行为以 MATLAB 靶为准（现场设计脚本 = 依据本体）；
> 方法学外部锚 = API RP 13D 7th ed. (2017)（水力学/摩阻方法学）；控压窗与 MPD 回压 = 呼1-004 现场设计口径。

---

## 0. 改动面总览（超出此清单 = 先报后动）

| # | 文件 | 动作 |
|---|---|---|
| 1 | `cemdisp/transport1d/hydraulics.py` | **新建**：P-2 Python 水力学内核（纯函数 + 一个装配入口，不接线生产 runner） |
| 2 | `cemdisp/data/pressure_field.py` | **追加** `TablePressureField`（二期；只加新类，既有 `PressureField/Constant/Hydrostatic` 零改动） |
| 3 | `tests/contract/test_hydraulics_ht1004.py` | **新建**：对靶逐点复现测试 |
| 4 | `tests/contract/test_table_pressure_field.py` | **新建**：TablePressureField 契约测试（含 oob 三件套同型） |
| 5 | `scripts/analysis/p3_product_report_20261007.py`（3.2 时建） | **新建**：P-4 产品报告驱动 |
| 6 | `results/P4产品_压力全链_20261007/`（3.2 时建） | **新目录**：三井产品报告 + 稠化窗校核 + 安全论 |

**明确不动**：`cemdisp/runners/*`、`annulus_d2dga.py`、`casing_flow.py`、敏感性 entrypoints、
`_SWITCH_DEFAULTS`、既有 pressure_field 三类。hydraulics 定位 = **探针/产品口径**（C-14 同型声明：
不进生产 runner；产品化另裁）。

## 1. 3.1 P-2 水力学内核设计

### 1.1 模块骨架（`cemdisp/transport1d/hydraulics.py`）

以「对靶逐点复现」倒推最小函数集（细节公式一律以 `P30_移植地图.md` + `p_jaifang1.m` 原文为准，
**逐字对齐、不改写**）：

```
load_structure(csv_path) -> StructureTable          # S1: 333 段，radius 列实为直径（已知陷阱）
build_segments(structure, thresholds7) -> Segments    # S2: 七段阈值 4025.73/5243.21/5578/7378.05/7521/7660
build_fluid_table() -> FluidTable                     # S4: 14 流体（L/min、L、g/cm³、mPa·s、Pa）
build_pump_schedule(...)                              # S7: 泵时表（节点非整数）
friction_pr(...)                                      # PR_Friction 移植（E=0 ⇒ 恒 1 空转，原样）
friction_annulus_bh(...)                              # Friction_annulus_bh 移植（30 步不动点，seed 100 Pa，tol 1e-3）
friction_casing_bh(...)                               # Friction_casing_bh 移植
track_annulus_interfaces(...)                         # S8-S9: 5 界面（自底累计容定位）+ VOF 线性混合
track_casing_interfaces(...)                          # S11-S13: 管内 14 界面（方向对偶、-1/10000 哨兵互换）
hydrostatic_and_friction_chain(...)                   # S10: 井口边界 bp·1e6 起算的静液/摩阻双链递推
pump_pressure_backtrack(...)                          # S13: 井底(=环空井底+P_bit_drop≡0)向上反推，负值钳 0
compute_ecd_esd(...)                                  # S15: ECD=压力/(9.81·TVD)、ESD 同式取静液压
control_windows(...)                                  # S16: 四点窗 [1.940,1.975]×3 + 井底 [2.010,2.050]；0.00981 因子
four_point_backpressure(...)                          # S21: 逐时走廊 BP_lower/upper、冲突取上界
bottom_window_backpressure(...)                       # S25: 持稳+冻结状态机（ECD_setpoint=2.025、≤5min 游程合并）
run_ht1004_target(...) -> dict[str, np.ndarray]       # 装配入口：复现靶全套 out_* 序列
```

### 1.2 原样复现红线（**修 bug = 失配 = 验收失败**）

移植地图登记的脚本自带缺陷与隐式状态**一律原样保留**，并在代码注释挂「MATLAB :行号 + 原样复现声明」：
- `:1153` 管内界面 2/3 的 τ 混合用 `tau1`（环空对偶处正确）——**照抄**；
- `:1464` `safety_margin=0.003` 只进 fprintf 不进走廊——**照抄**；
- `:733` 循环内反复赋 `volume_injected_casing_1_list(1)=1200` ⇒ 首界面 +1200 L 偏置——**照抄**；
- `n_time = floor(189.059) = 189`（尾丢 ≈41 L）——**照抄**；
- 摩擦不动点：30 步上限、seed 100 Pa、tol 1e-3——**照抄**；
- `E=0 ⇒ PR_Friction 恒 1`（空转）——**照抄**；
- 单位混用点（mm/m、MPa/Pa、g/cm³/kg/m³、0.00981 与 9.81 两族因子）按地图 §4 单位总表逐处对齐。

### 1.3 `TablePressureField`（pressure_field.py 二期）

- 协议同 `PressureField`：`P(md_m, t_s) -> MPa`；数据 = (md 节点, 时间列) 二维表 + 双线性插值；
- **oob 三件套同型**（越界计数/审计/钳位语义与 Hydrostatic 一致）；
- 用途：Phase 4 逐深温度/压力产品位；本阶段只落地 + 单测（用靶的 `out_annuli_bottom_pressure.csv`
  作冒烟数据源），**不注入生产路径**。

### 1.4 验收（3.1）

- **逐点复现**（`test_hydraulics_ht1004.py`，靶 = `sandbox/out_*.csv` 直读）：
  | 对象 | 判据 |
  |---|---|
  | 泵压时程（189） | max 相对差 ≤ 1e-9（超限逐点归因，禁调参） |
  | 井底 ECD 时程（189） | max 绝对差 ≤ 1e-9 g/cm³ |
  | ECD 全场（189×333） | max 绝对差 ≤ 1e-9 |
  | 四点走廊 BP_lower/upper/conflict（189×3） | conflict 逐位一致；BP ≤ 1e-9 MPa |
  | 井底窗持稳回压（189） | ≤ 1e-9 MPa + 状态机事件时刻逐位（交点 189min/冻结 5.117） |
  | 对比模块插值（reference 116 行） | unique-stable+钳位+linear 口径一致，设计列逐位 |
- **关1**：全量 pytest 17F 不增；`test_default_path_bitwise_anchor` 绿（新模块不接线 ⇒ 平凡成立，仍必跑）。
- **关2**：`pressure_field.py` 既有三类零改动（git diff 审查 + 既有 20 条测试绿）；hydraulics 无任何
  生产调用点（grep 审查）。
- **关3**：方向性——泵压随排量**正相关**（Ruling R-3.1-1，2026-10-07：判据 = 相关系数 > 0.5，
  阶跃界面期存在物理性局部回勾，严格逐点单调不成立；对抗评审 MINOR-1 登记）、ECD 随静液柱深度增、
  回压走廊下界≤上界处无冲突（靶 189/189 冲突如实）。
- **关5**：全程无 NaN/Inf；体积链闭合（累计注入 vs 靶 `Cumulative_Volume_m3` 逐位）。
- **排量响应锚**：本阶段零数值路径改动 ⇒ 引用 Phase 2e 实跑结果（两套单调 ✓）+ 全量回归位级锚，
  不重跑 2D（登记于台账）。

## 2. 3.2 P-4 产品层设计

- **三井泵压/ECD/控压窗余量报告**：呼1-004 = 靶全链（190 min 泵注程序）；呼101/呼1-003 = 内核跑其
  井身结构+泵注程序（`现场资料提取/<井>/` CSV），**控压窗数据缺口如实报告**（窗口值为呼1-004 设计口径，
  另两井无设计窗 ⇒ 只出泵压/ECD 时程 + 当量密度余量对地层压力当量 1.933（呼1-004 已有锚））。
- **稠化窗校核（D7-A）**：泵注程序 vs 井下 T,P 稠化余量；P 缩稠化至 1/3 口径：2.51→0.85 h @75°C、
  5→200 MPa；>110°C 水化产物警示区加不确定度标注。
- **作业史逐段文字对账**：⚠「9–24.6 MPa」在 `pumping_schedule.csv` **notes 文字列**（非数值列）。
- **安全论成文**：「忽略温压流变低估 ECD 0.06 s.g.（OBM 口径）≈ 1.7× 窗宽 0.035」+ 对照
  `bottom_ecd_reference.csv` + 189/189 冲突与控压后鞋 ECD 超窗 [2.012,2.073] vs [1.940,1.975] 如实呈现。
- 验收 = 三井报告齐 + 稠化窗表 + 安全论一段 + 引用注记全（靶运行注记/历史产物注记/窗口缺口注记）。

## 3. 3.3 P-3 ρ(T,P)：**本期不启动**（D6-A 二期另案；启动 = 基线作废级，§8-9）。

## 4. 回退

全部加法式（新文件 + pressure_field.py 只追加类）：`git revert` 单提交即回；测试文件独立可删。

## 5. 风险与开放问题

1. 不动点迭代/插值的浮点次序差异可能造成 1e-12~1e-9 级偏差——判据已留 1e-9；超限**先归因后放宽**，
   放宽须登记台账（不得静默）。
2. 呼101/呼1-003 井身结构 CSV 的「radius 列实为直径」陷阱（legacy 已知）——load_structure 必须核列语义。
3. MATLAB `unique(...,'stable')`/`interp1` 钳位语义与 numpy 的等价实现须显式单测钉住。
