# D2DGA Tier 1 T1-7 CFL 自适应时间步设计

> **日期**：2026-07-20
> **范围**：Tier 1 剩余 T1-7（CFL 自适应时间步，L4 单文件）
> **关联**：`2026-07-18-d2dga-tier1-design.md`（T1-1~T1-5）、`2026-07-19-d2dga-tier1-t1-6-flusher-design.md`（T1-6）、调研报告
> **执行模式**：主 agent 统筹，子 agent 执行（haiku，spec §6 作 plan 用）
> **验证约束**：⚠️ 不对标 CBL，enable_cfl_adaptive 开关保基线对照

---

## 0. 背景与设计约束

### 0.1 前置状态

- **T1-1~T1-6 已完成**（commits `5cb03b0..c570a3e`，255 测试通过）：去通量限幅、I3 物理化（+局部 CFL 裁剪 `alpha_cfl`）、体力向量注入、两层闭包、static wall c_min、FLUSHER 独立浓度场。
- **T1-7 是 Tier 1 剩余第二项**：加全局 CFL 自适应时间步（每步算 `max|w|/Δs + max|v|/Δy` 动态调 dt）。

### 0.2 现状（审计核实）

- `dt` 固定 `4.0`（行 194），`total_t=12000`（行 197），`save_interval=60` 按步数（行 204）。
- `final_step_index = int(total_t/dt)`（行 730），`for step_index in range(...)`（行 731-732）。
- 快照 `step_index % save_interval`（行 887）。
- T1-2 `alpha_cfl=0.5`（行 207）是 **I3 局部通量裁剪**（`|div_q|·dt ≤ α·Δs`，行 840-842），**非全局时间步**。
- 半拉格朗日平流（行 775-780）；`_compute_velocity` 返回 `w/v`（行 649/652-659）。

### 0.3 设计原则

1. **加开关**：`enable_cfl_adaptive`（默认 `True` 自适应，`False`=固定 dt 复现旧 R0-R3 基线）。与 R0-R3 消融框架一致，可对照。
2. **T1-2 局部 CFL 保留**：`alpha_cfl`（I3 div_q 裁剪）与 T1-7 全局 `cfl_number` 独立，文档区分。
3. **save_interval 保留按步数**：`snapshot_times_s` 记录实际时间，下游不破坏。
4. **dt 上下限**：`dt ∈ [dt_min, dt_user]`，`dt_user=4s` 上限，`dt_min=0.1` 防爆炸。
5. **半拉格朗日 CFL 保守估计**：`cfl_number < 1`（默认 0.5）。
6. **不对标 CBL**：开关 False 复现基线 + True 不崩溃。

---

## 1. 改动总览

| 层 | 文件 | 改动 | 风险 |
|----|------|------|------|
| L4 | `cemdisp/models2d/annulus_d2dga.py` | `enable_cfl_adaptive`/`cfl_number`/`dt_min` 参数 + `_compute_velocity` 后 dt 自适应 + `run` 循环 while 分支 + 最后步裁剪 | 中 |

单文件，不跨层。T1-2 局部 `alpha_cfl` 保留不动。

---

## 2. 详细设计

### 2.1 构造参数（`__init__`，行 190-210 区）

新增：
- `enable_cfl_adaptive: bool = True`（默认自适应，False 固定 dt 复现基线）
- `cfl_number: float = 0.5`（全局 CFL 数，独立于 T1-2 `alpha_cfl`）
- `dt_min: float = 0.1`（dt 下限，防 CFL 过小步数爆炸）

存储 `self.enable_cfl_adaptive` / `self.cfl_number` / `self.dt_min`。`dt`/`alpha_cfl` 保留。

### 2.2 dt 自适应逻辑（`_compute_velocity` 返回 w/v 后，`run` 循环内）

```python
if self.enable_cfl_adaptive:
    ds = float(np.min(np.diff(geom["s"])))
    dy_arr = np.gradient(geom["y"])[:, None]
    dy_min = float(np.min(dy_arr)) if dy_arr.size else ds
    denom = max(float(np.max(np.abs(w))) / max(ds, 1e-12)
                + float(np.max(np.abs(v))) / max(dy_min, 1e-12), 1e-12)
    dt_cfl = self.cfl_number / denom
    dt_step = min(dt_cfl, self.dt, self.total_t - current_time_s)
    dt_step = max(dt_step, self.dt_min)  # 防爆炸
    # 最后步裁剪精确到 total_t
    if current_time_s + dt_step >= self.total_t:
        dt_step = self.total_t - current_time_s
else:
    dt_step = self.dt  # 固定 dt 复现基线
```

### 2.3 `run` 循环改造（行 730-732 区）

- `enable_cfl_adaptive=True`：改 `while current_time_s < self.total_t - 1e-9:`，每步 `current_time_s += dt_step`，`step_index += 1`。
- `False`：保留 `final_step_index = int(total_t/dt)` + `for step_index in range(...)` + `current_time_s = step_index * dt`（复现基线）。
- 快照触发（行 887）：`if step_index % self.save_interval == 0 or current_time_s >= self.total_t - 1e-9:`（`snapshot_times_s` 记录 `current_time_s` 实际值）。

### 2.4 T1-2 局部 CFL 关系

T1-2 `alpha_cfl`（I3 `div_q` 裁剪，行 840-842）保留。`dt` 变化时 `step_limit = alpha_cfl·Δs/dt` 随 `dt` 变（一致）。`alpha_cfl` 与 `cfl_number` 独立参数，文档区分（局部通量裁剪 vs 全局时间步）。

---

## 3. 集成与数据流

- **`snapshot_times_s`**：记录每快照实际 `current_time_s`（dt 变化致间隔不均，但下游用 `snapshot_times_s` 不破坏；Tier0 `flow_classification`/`c_bar_st` 用之）。
- **`save_interval` 按步数**：dt 自适应通常 dt 变小→步多→快照更密（好）；dt 大时快照疏（`snapshot_times_s` 仍记录实际时间）。
- **ablation 基线**：`enable_cfl_adaptive=False` 复现 T1-6 后基线（R0-R3 效率一致）；`True` 重跑（数值扩散/前锋略变）。`run_one_level`（`runners/ht1_004_ablation.py`）加 `enable_cfl_adaptive` 透传（默认 False 保 ablation 基线复现，或 True 重跑）。
- **Tier0 诊断复跑**：`snapshot_times_s` 实际时间，`compute_all_tier0_diagnostics` 不改（消费 result）。
- **`AnnulusSimulationResult` 不改结构**：`snapshot_times_s` 已存在。

---

## 4. 测试策略

| 测试 | 内容 | 文件 |
|------|------|------|
| dt 自适应 | `enable_cfl_adaptive=True` 时 dt≤dt_user、dt≥dt_min、CFL=cfl_number/denom<1 | `tests/test_improved_d2dga_annulus.py` |
| 最后步裁剪 | 循环结束 `current_time_s ≈ total_t`（误差<1e-6） | 同上 |
| 开关 False 复现 | `enable_cfl_adaptive=False` 与 T1-6 后基线一致（固定 dt for 循环） | 同上 |
| snapshot 时间 | `snapshot_times_s` 记录实际时间（非等间隔） | 同上 |
| 回归 | R0-R3 消融（False 复现）+ 六井集成（不对标 CBL） | `test_improved_d2dga_annulus.py` + `test_six_well_integration.py` |

---

## 5. 验证门（无 CBL）

| 验证类型 | 内容 |
|---------|------|
| 自适应稳定 | `enable_cfl_adaptive=True` 不崩溃 + CFL<1 + dt∈[dt_min, dt_user] + 最后步精确到 total_t |
| 基线复现 | `False` R0-R3 效率与 T1-6 后一致（固定 dt 复现） |
| 回归 | R0-R3 消融（False）+ 六井集成，不崩溃+物理合理 |
| 诊断复跑 | hu102/ht1_004 Tier0 诊断：`snapshot_times_s` 实际时间，muskat/浮力/η_N 不变 |

---

## 6. 执行计划（作 plan 用，haiku 子代理）

### 6.1 Task 分解（bite-sized TDD）

- **Task 1**（参数 + dt 自适应逻辑）：`annulus_d2dga.py` `__init__` 加 `enable_cfl_adaptive=True`/`cfl_number=0.5`/`dt_min=0.1` + 存储；`run` 内 `_compute_velocity` 返回 w/v 后加 dt 自适应逻辑（§2.2 代码）；测试 `enable_cfl_adaptive=True` 时 dt≤dt_user、dt≥dt_min、CFL<1。
- **Task 2**（run 循环 while 分支，依赖 Task1）：`run` 循环改——`True` 用 `while current_time_s < total_t` + `current += dt_step` + 最后步裁剪；`False` 保留固定 dt for 循环；快照触发改 `current_time_s >= total_t` 末步；测试最后步 `current_time_s≈total_t` + `enable_cfl_adaptive=False` 复现固定 dt。
- **Task 3**（集成验证）：`runners/ht1_004_ablation.py` `run_one_level` 加 `enable_cfl_adaptive` 透传（默认 False 保 ablation 基线）；全量 `pytest tests/` + R0-R3 消融（False 复现 + True 不崩溃）+ 六井 + Tier0 诊断复跑 + commit。

### 6.2 子 agent 约束

- haiku implementer + haiku reviewer（用户指定）。
- 中文 docstring + 类型注解 + 公式标注。
- 数值安全：dt 下限 `dt_min`、上限 `dt_user`、CFL<1、最后步裁剪、除零守卫。
- 环境：`conda activate shenjingwangluo`；`/d/apps/Anaconda/python` 加 `PYTHONIOENCODING=utf-8 PYTHONUTF8=1`。
- TDD：每 Task 写失败测试→跑 FAIL→实现→跑 PASS→commit。

---

## 7. 风险与缓解

| 风险 | 缓解 |
|------|------|
| snapshot 时间对齐（间隔不均） | `snapshot_times_s` 记录实际 `current_time_s`，下游用之不破坏 |
| 停止时间精度 | 最后步裁剪 `dt_step = total_t - current_time_s`，断言 `current≈total_t` |
| ablation 基线漂移 | `enable_cfl_adaptive=False` 复现固定 dt 基线；`run_one_level` 透传开关 |
| I3 局部裁剪与全局 CFL 耦合 | `step_limit=alpha_cfl·Δs/dt` 随 dt 变（一致）；`alpha_cfl`/`cfl_number` 独立 |
| 半拉格朗日 CFL 保守 | `cfl_number=0.5<1`；`dt_min=0.1` 防爆炸 |
| dt_min 过大欠解 | `dt_min=0.1` 可配；测试 dt 不被 dt_min 截断（除非 CFL 要求更小） |

---

## 8. 后续衔接

- T1-7 完成 → CFL 自适应就绪，去限幅后数值稳定性提升。
- Tier 1 剩余：T1-8 1D Taylor 弥散（L2 独立）/ T1-9 入口方位加权（L3+L4）。
- Tier 2：HB D2DGA 闭包查表 / augmented Lagrangian / TVD / 瞬态流函数 / T2-5 序列优化器（用 T1-6 flusher）。

