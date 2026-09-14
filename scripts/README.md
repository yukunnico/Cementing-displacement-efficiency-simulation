# `scripts/` 索引（2026-09-14 归档）

本目录 2026-09-14 由平铺的 75 个脚本重组为 5 个子目录 + 本索引。**只做归档与索引，未改任何脚本的算法逻辑。**

## 目录职责

| 目录 | 职责 | 数量 |
|---|---|---|
| `entrypoints/` | **活跃入口**：会产生论文/正式口径数字、或位于生产链路上的脚本 | 14 |
| `lib/` | **活跃共享库**：被多个脚本 import 的运行骨架（不单独运行） | 2 |
| `reruns/` | 历史：各次修复后的 8 井重跑（一代一文件，按日期冻结） | 7 |
| `probes/` | 历史：诊断/取证/敏感性/消融探针（一次性，多为 spike） | 44 |
| `plots/` | 出图脚本：效率沿深度分布、CBL 对照、专利附图 | 7 |
| （根） | 2 个已废弃脚本，改名为 `*.DEPRECATED.py` 以断开 import | 2 |

## 运行约定

- **必须在仓库根目录（`cement model/`）下运行**：脚本用相对路径读写 `results/`，且自身把仓库根注入 `sys.path`。
  ```bash
  PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/entrypoints/<脚本名>.py
  ```
- 跨脚本引用统一用**完全限定名**（`scripts/` 是隐式命名空间包，无需 `__init__.py`）：
  ```python
  from scripts.lib.mass_balance_diag import WELLS, build_case   # 共享骨架
  from scripts.entrypoints.rerun_all_wells_corrected import CORRECTED_KW, WELLS
  ```
- 每个脚本内用 `Path(__file__).resolve().parents[2]` 定位仓库根（`scripts/<子目录>/x.py` 往上三级）；
  少数脚本（如 `entrypoints/run_grid_convergence.py`、`entrypoints/closure_contribution_scan.py`、
  `plots/plot_segment_quality.py`）不自算根目录，直接按相对路径读写 `results/`，**因此更必须在仓库根运行**。
- 历史脚本的输出目录名（如 `results/守恒修复后全量重跑_2026-09-02`）是冻结的历史口径，**不要与新基线数字混用**；权威口径见 `results/<井名>_1D2D耦合模型/`。
- `probes/ht1_004_sensitivity.py` 文件带 UTF-8 BOM（归档前既有），用 `encoding="utf-8"` 读取会报 `U+FEFF`，请用 `utf-8-sig`。

## 活跃面 —— `entrypoints/`

| 脚本 | 一句话用途 | 状态 | 输出目录 | 依赖 |
|---|---|---|---|---|
| `rerun_all_wells_corrected.py` | **8 井 baseline vs corrected 全量重跑（nz=250）**，论文/正式数字的官方口径入口，并导出 `adopted_config.json` 开关快照 | 活跃 | `results/全井修正前后/` | `cemdisp.{data,models2d,transport1d}`、numpy |
| `cbl_window_comparison.py` | 模型评价窗 × 现场 CBL 窗口逐井对照（论文素材，B1）；支持 `--all/--well/--tables-only/--force` | 活跃 | `results/最终基线_2026-08-29/` | 同上 + `scripts.entrypoints.rerun_all_wells_corrected` |
| `calibrate_dispersion_yield.py` | F1 弥散量级 + F2 屈服门槛安全系数的二维标定扫描（8 井，`ProcessPoolExecutor` 并行） | 活跃 | `results/弥散产额标定_2026-08-31/` | 同上 + numpy |
| `closure_contribution_scan.py` | 三闭包 R0→R3 逐级贡献分解（另含 NONE 纯平流基线，共 5 级） | 活跃 | `results/三闭包贡献分解/` | `cemdisp.{data,models2d,transport1d,runners}` |
| `verify_loader_against_field_v1.py` | `cemdisp` loader 载入值 vs 现场提取 CSV 的逐井校验（v1，2026-09-01） | 活跃 | `results/数据校验_2026-09-01/loader_dump/` | `cemdisp.data` |
| `export_depth_time_concentration.py` | 环空各深度×各流体浓度随时间变化的导出（长/宽格式 CSV + 说明 md）；`--fullwell` 模式写 `results/<井名>_1D2D耦合模型/`，也被 ht1-004 runner 作为库导入 | 活跃 | `results/<井名>_1D2D耦合模型/`（或调用方指定目录） | `cemdisp.{data,models2d,transport1d}`、numpy |
| `run_ht1_004_ablation.py` | HT1-004 R0→R3 消融运行入口（产出各级摘要 JSON + 汇总 CSV） | 活跃 | `results/ht1_004_ablation/` | `cemdisp.runners` |
| `run_grid_convergence.py` | 网格收敛（nz=140/280/500）+ 时间步敏感性（dt=2/8）验证，支撑论文图9/表6 | 活跃 | `results/ht1_004_ablation/` | `cemdisp.runners` |
| `run_theoretical_case.py` | 理论算例：偏心度放大到 35%/45%，检验 I3/真体力的改进能力 | 活跃 | `results/ht1_004_ablation/` | `cemdisp.{data,runners}` |
| `run_tier0_diagnostics.py` | Tier 0 诊断独立验证（默认注册 hu102/ht1_004，`--well` 选择） | 活跃 | `results/tier0_diagnostics/` | `cemdisp.{data,diagnostics,models2d,runners,transport1d}` |
| `run_ablation_variants.py` | runner 口径单变量变体实验（消融组，2026-09-09 裁定） | 活跃 | `results/消融变体_runner口径_2026-09-09/` | `cemdisp.{data,models2d,transport1d}`、numpy |
| `run_gap_fill_variants.py` | 论文正文补齐变体实验（runner 口径，2026-09-10 裁定） | 活跃 | `results/正文补齐变体_runner口径_2026-09-10/` | 同上 + numpy |
| `run_sensitivity_variants.py` | 现场杠杆单变量敏感性组实验（runner 口径，2026-09-10 裁定） | 活跃 | `results/敏感性变体_runner口径_2026-09-10/` | 同上 + numpy |
| `smoke_all_wells.py` | 冒烟：依次调 7 井 runner 的 `*_initial` 入口，汇总成功/失败与关键指标（原 `test_all_wells.py`） | 活跃 | `results/all_wells_test_result.json` | `cemdisp.runners` |

## 活跃面 —— `lib/`

| 脚本 | 一句话用途 | 状态 | 输出目录 | 依赖 |
|---|---|---|---|---|
| `mass_balance_diag.py` | 质量平衡（库存）取证共享骨架：`WELLS` 井表 + `build_case` / `run_variant` / `integrate_injection` + 结果落盘（原 `_mass_balance_diag_20260902.py`） | 活跃（库） | `results/_质量平衡取证_2026-09-02/` | `cemdisp.{data,models2d,transport1d}`、numpy |
| `sensitivity_common.py` | C4/弥散敏感性脚本共享运行骨架：单 case = 基线同款流水线 + 可覆盖开关/井规格（原 `_sensitivity_common.py`） | 活跃（库） | 由调用方指定 | 同上 + `scripts.entrypoints.rerun_all_wells_corrected` |

## 历史面 —— `reruns/`（各次修复后的 8 井重跑，一代一文件）

| 脚本 | 一句话用途 | 状态 | 输出目录 | 依赖 |
|---|---|---|---|---|
| `rerun8_clean_20260902.py` | B1/B2/B3 三项修复后的 8 井 nz=250 干净全量重跑（当时的权威表） | 历史 | `results/守恒修复后全量重跑_2026-09-02/` | `scripts.lib.mass_balance_diag`、numpy、pandas |
| `rerun8_after_pipe_capacity_fix_20260902.py` | 管容双链分裂修复后的 8 井终跑（权威表第二轮） | 历史 | `results/管容链修复后终跑_2026-09-02/` | 同上 |
| `rerun8_after_plug_semantics_fix_20260903.py` | 胶塞语义修复（替浆顶胶塞、不进环空）后的 8 井终跑 | 历史 | `results/胶塞语义修复后终跑_2026-09-03/` | 同上 |
| `rerun8_three_fixes_20260906.py` | 三项修复（WASH/SPACER 选相 + e_clip 裁定 + 幂律缝隙律）的 8 井量化重跑 | 历史 | `results/三项修复重跑_2026-09-06/` | 同上 |
| `rerun8_mixing_contact_time_20260908.py` | 路线 B（管内浓度剖面真解）Task 3：8 井双开关验证 + alpha 消融前置；支持 `smoke`/`full` | 历史 | `results/浓度剖面真解路线B_2026-09-08/` | 同上 |
| `rerun_stop_fix_20260901.py` | F2 停止标志修复的量化：去除 +600s 尾窗 与 环空域口径 两个 2D 效应 | 历史 | `results/停止标志修复与裸眼域对比_2026-09-01/` | `cemdisp.*`、numpy |
| `_rerun8_after_conservation_fix_20260902.py` | 阶段4：因子2守恒修正后按 runner 生产口径（nz=250）重算 8 井 | 历史 | `results/_质量平衡取证_2026-09-02/` | `scripts.lib.mass_balance_diag`、numpy、pandas |

## 历史面 —— `probes/`（诊断 / 取证 / 敏感性 / 消融）

### 质量平衡与守恒取证（2026-09-02 系列）

| 脚本 | 一句话用途 | 状态 | 输出目录 | 依赖 |
|---|---|---|---|---|
| `_audit_pipe_volume_split_20260902.py` | 独立审计：管容双链分裂逐井核查 | 历史 | 控制台 | `cemdisp.*` |
| `_pure_plug_conservation_20260902.py` | 决定性守恒实验：合成"纯水泥单相塞流"入口，剥离多相/过量/偏心干扰 | 历史 | `results/_质量平衡取证_2026-09-02/` | `cemdisp.*`、numpy |
| `_wallfreeze_grid_mechanism_20260902.py` | 定位"残余损耗 + 网格不收敛 + 关弥散崩塌"的共同病根（阶段5） | 历史 | 同上 | `scripts.lib.mass_balance_diag`、numpy |
| `_yieldgate_verify_20260902.py` | 验证物理 τw 屈服门（去 `residual_wall`）的网格收敛性与物理窜槽保留（阶段6） | 历史 | 同上 | 同上 |
| `_make_comparison_chart_20260902.py` | 生成八井修正前后 η_E 对比图与汇总 CSV | 历史 | 同上 | numpy、pandas、matplotlib |
| `bisect_hu103_20260902.py` | hu103 η≈0 二分定位：逐一叠加 `CORRECTED_KW` 成员，找"救活"开关 | 历史 | 控制台 | `cemdisp.*` |
| `debug_hu1_hu103_eta_zero_20260902.py` | hu1/hu103 η≈0 回归定位 | 历史 | 控制台 | `cemdisp.data` |
| `verify_fix_hu1_hu103_20260902.py` | 修复验证：hu1/hu103 baseline nz=60 复跑，确认死锁解除 | 历史 | 控制台 | `cemdisp.data` |
| `corrected_ref_hu1_hu103_20260902.py` | `corrected_full` 配置参照值（hu1/hu103 nz=60），确保修复不改坏已验收口径 | 历史 | 控制台 | `cemdisp.data` |
| `c_verify_convergence_20260902.py` | C 修复验收 2-4：corrected 逐位回归 + baseline nz 收敛（nz=60/150/250） | 历史 | 控制台 | `cemdisp.data` |
| `isolate_fix_mechanisms.py` | 隔离各修正机制，找出 corrected 配置 wall=100%/η_E=0.03 的元凶 | 历史 | 控制台 | numpy、`cemdisp.*` |
| `compare_old_new_model.py` | 对比新旧套管段模型的顶替效率差异 | 历史 | 控制台 | `cemdisp.*` |

### 停止时刻 / 尾窗 / 管容口径

| 脚本 | 一句话用途 | 状态 | 输出目录 | 依赖 |
|---|---|---|---|---|
| `dump_stop_times_v1.py` | 顶替结束标志取证：1D 时间轴实算（只读，不改仓库数据） | 历史 | 控制台 | `cemdisp.data` |
| `dump_stop_times_v2.py` | F2 修复验证：逐井对比停止标志修复前后的 1D 时间轴 | 历史 | 控制台 | `cemdisp.*` |
| `dump_tailwindow_2d_v1.py` | 验证 RR +600s 尾窗对 η_E 的实际影响（2D 小网格实跑对比） | 历史 | 控制台 | `cemdisp.data` |
| `_diag_hu103_stop_20260902.py` | 诊断 hu103（对照 hu102）1D 停止时刻：为何 `cement_end_time` 时只有部分水泥跨鞋口 | 历史 | 控制台 | numpy、`cemdisp.*` |
| `_hu103_stop_diag_20260902.py` | 同上（早期版本，重复实现，保留备查） | 历史 | 控制台 | numpy、`cemdisp.*` |
| `_diag_hu103_gravity_20260902.py` | 对比 hu103 重力修正/弥散开关下尾浆跨鞋口体积，定位尾浆被截原因 | 历史 | 控制台 | numpy、`cemdisp.*` |

### 敏感性扫描

| 脚本 | 一句话用途 | 状态 | 输出目录 | 依赖 |
|---|---|---|---|---|
| `dispersion_scale_sensitivity_scan.py` | M1 弥散系数 κ 缩放敏感性扫描（CFL on/off 双组存档） | 历史 | 控制台 | numpy、`cemdisp.*` |
| `dispersion_scale_sensitivity_nz250.py` | NZ=250 弥散敏感性附表（C1 方案 A）：8 井 × `dispersion_dt_scale ∈ {0.0, 0.5}` | 历史 | `results/最终基线_2026-08-29/` | `scripts.lib.sensitivity_common`、`scripts.entrypoints.rerun_all_wells_corrected` |
| `standoff_sensitivity_c4.py` | C4：standoff ±0.1 敏感性（修正后配置 nz=250 CFL on），论文素材 | 历史 | 同上 | `scripts.lib.sensitivity_common`、`cemdisp.data` |
| `standoff_sensitivity_scan.py` | 偏心度（居中度）敏感性参数扫描 | 历史 | 控制台 | numpy、`cemdisp.*` |
| `density_contrast_sensitivity_scan.py` | 密度差（水泥浆 vs 钻井液）敏感性参数扫描 | 历史 | 控制台 | numpy、`cemdisp.*` |
| `ht1_003_sensitivity.py` | 呼1-003 井尾浆参数敏感性分析（支持 `--resume`/`--plot-only`） | 历史 | `results/呼1-003_敏感性分析/` | matplotlib、pandas、numpy、tqdm |
| `ht1_004_sensitivity.py` | 呼1-004 井尾浆参数敏感性分析（支持 `--resume`/`--plot-only`；文件带 BOM） | 历史 | `results/呼1-004_敏感性分析/` | matplotlib、pandas、numpy、tqdm |
| `i3_ablation_rerun_20260829.py` | I3 修复后 R2 消融重测 | 历史 | `results/i3_ablation_rerun_2026-08-29/` | `cemdisp.*` |
| `analyze_calibrated_rerun_20260829.py` | 校准后重跑 vs 校准前基线 + CBL 对比分析 | 历史 | `results/校准后重跑_2026-08-29/` | `cemdisp.*`、numpy |
| `analyze_distortion_fix.py` | 环空失真修正前后对比分析（M0-M4 + I3） | 历史 | `results/失真修正前后对比/` | `cemdisp.*`、numpy |

### 剖面提取与可视化探针

| 脚本 | 一句话用途 | 状态 | 输出目录 | 依赖 |
|---|---|---|---|---|
| `extract_profiles.py` | 快速提取窄边水泥浓度剖面和鞋口出流时序（Step 3.3 & 3.4） | 历史 | `results/p3_p4_integration/` | `cemdisp.*` |
| `analyze_profiles.py` | 分析窄边水泥浓度剖面和鞋口出流时序并出对比图（承接 `extract_profiles` 产物） | 历史 | 同上 | `cemdisp.*`、matplotlib |

### 居中度 / 低分归因探针（hu101 系列，2026-09-10~11 spike）

| 脚本 | 一句话用途 | 状态 | 输出目录 | 依赖 |
|---|---|---|---|---|
| `hu101_standoff_finetune_probe.py` | 呼101 井居中度微调探针（找阈值位置） | 历史 | `results/呼101_居中度微调探针_2026-09-10/` | `cemdisp.*`、matplotlib |
| `hu101_standoff_measured_vs_assumed.py` | 呼101 井居中度"实测值 vs 模型假设值"敏感性对比（八井居中度出处取证用） | 历史 | `results/呼101_居中度实测对比/` | `cemdisp.*`、numpy |
| `hu101_vs_hu102_diff_attribution_20260911.py` | 呼101 vs 呼102 差异归因探针（第三轮） | 历史 | `results/呼101_vs_呼102_差异归因_2026-09-11/` | numpy、`cemdisp.*` |
| `hu101_gap_vs_inventory_decouple_20260911.py` | 间隙 vs 库存比 解耦探针（第四轮） | 历史 | `results/呼101_间隙vs库存比解耦_2026-09-11/` | numpy、`cemdisp.*` |
| `hu101_standoff_response_tuning_survey_20260911.py` | 居中度→顶替效率 响应强度调研：曲线是斜坡还是悬崖（第五轮） | 历史 | `results/居中度响应强度调研_2026-09-11/` | numpy、`cemdisp.*` |
| `hu101_low_score_attribution_probe_20260911.py` | 呼101 低分归因探针 | 历史 | `results/呼101_低分归因探针_2026-09-11/` | numpy、`cemdisp.*` |
| `hu101_low_score_probe2_20260911.py` | 呼101 第二次归因探针：阈值锐度 + 结构开关复检 | 历史 | `results/呼101_低分归因探针2_2026-09-11/` | numpy、`cemdisp.*` |
| `hu101_standoff_and_temperature_scenarios_20260911.py` | 呼101 居中度口径 × 温度口径 全场景探针（第六轮） | 历史 | `results/呼101_居中度与温度口径场景_2026-09-11/` | numpy、`cemdisp.*` |
| `hu101_temp_and_eclip_followup_20260911.py` | 呼101 温度口径 + `e_clip` 分叉补测（第六轮续） | 历史 | 同上 | numpy、`cemdisp.*` |
| `hu101_rheology_sensitivity_20260829.py` | hu101 领/尾浆"主检 vs 复检"双口径流变敏感性 | 历史 | `results/校准后重跑_2026-08-29/` | `scripts.entrypoints.rerun_all_wells_corrected`、`cemdisp.data` |

### 井间对照实验（名为 test_，实为对照探针，**不是** pytest 用例）

| 脚本 | 一句话用途 | 状态 | 输出目录 | 依赖 |
|---|---|---|---|---|
| `test_dual_factor_ht1_004.py` | 双因素对照：在呼1-004 上分离"居中度"与"流体体系流变"各自贡献 | 历史 | `results/_test_standoff/` | matplotlib、`cemdisp.*` |
| `test_geometry_rate_factor.py` | 呼1-004 几何/泵速影响因素分析 + 呼1-002 反向对照 | 历史 | 同上 | 同上 |
| `test_standoff_reverse_ht1_004.py` | 反对照实验：呼1-004 居中度 0.83→0.60，看效率是否暴跌 | 历史 | 同上 | 同上 |
| `test_standoff_sensitivity.py` | 把呼探1-001/002 居中度改为设计书 6.3 节设计模拟值的测试 | 历史 | 同上 | 同上 |

## 历史面 —— `plots/`

| 脚本 | 一句话用途 | 状态 | 输出目录 | 依赖 |
|---|---|---|---|---|
| `plot_ht1_001_efficiency_depth_distribution.py` | 呼探1-001 尾管段顶替效率沿深度分布图（学术风格） | 历史 | `results/呼探1-001尾管_1D2D耦合模型/` | matplotlib、numpy |
| `plot_ht1_004_efficiency_depth_distribution.py` | 呼1-004 尾管段顶替效率沿深度分布图（学术风格） | 历史 | `results/呼1-004_1D2D耦合模型/` | 同上 |
| `plot_wells_efficiency_depth_distribution.py` | 各井尾管段顶替效率沿深度分布图（学术风格，多井批量） | 历史 | `results/<井名>尾管_1D2D耦合模型/` | 同上 |
| `plot_hu101_cbl_quality_vs_efficiency.py` | 呼101 尾管井 CBL 胶结质量 vs 模型顶替效率分段对比图 | 历史 | `results/呼101尾管_1D2D耦合模型/` | matplotlib、pandas、numpy |
| `plot_hu101_full_depth_distribution.py` | 呼101 尾管井：模型计算顶替情况 vs 现场固井质量 逐段对照图 | 历史 | 同上 | matplotlib、pandas、numpy |
| `plot_patent_figures.py` | 专利技术交底书 12 张附图批量生成 | 历史 | `results/`（成组图集） | matplotlib、numpy |
| `plot_segment_quality.py` | 呼1-003 尾管段按 200m 分段顶替质量评价（图 + CSV，合格/良好线判定） | 历史 | `results/呼1-003_1D2D耦合模型/` | pandas、numpy、matplotlib |

注：`plot_segment_quality.py` 原在仓库根，2026-09-14 迁入本目录；它按**当前工作目录**解析 `results/` 相对路径，须在仓库根运行。

## `scripts/` 根：已废弃（仅改名断 import，未删除）

| 脚本 | 说明 | 状态 |
|---|---|---|
| `p3_p4_integration.DEPRECATED.py` | 标定/对比逻辑锚定已删除的 CBL 评价段指标列 `cbl_eval_interval_efficiency`，运行即 KeyError。替代方案：`entrypoints/cbl_window_comparison.py`、`entrypoints/rerun_all_wells_corrected.py` | 废弃 |
| `run_calibrated_all_wells.DEPRECATED.py` | 经 `from scripts.p3_p4_integration import ...` 依赖同一已删列；标定参数 M=500/alpha=0.05 来自已废弃旧扫描 | 废弃 |

## 备份

- `archive/scripts_untracked_2026-09-14/` —— 归档前 11 个**未被 git 跟踪**的脚本原件备份（8 个 `hu101_*` 探针 + 3 个 `run_*_variants.py`），以及 2 个当时已跟踪的 hu101 探针副本；移动操作不可逆，此备份为唯一还原来源。
- `archive/scripts_pycache_stale_2026-09-14/` —— 归档前的过期 `scripts/__pycache__`（字节码，已被 .gitignore 排除，未入库）。
