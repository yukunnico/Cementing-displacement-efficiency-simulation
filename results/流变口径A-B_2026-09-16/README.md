# 流变口径 A-B（阶段 B-4，2026-09-16）

本目录是阶段 B 开关效应的 **nz=60 粗网格相对对照**产物（`ablation_8wells.csv`：
8 井 × 5 变体 V0_baseline / V1_yield_gate / V2_power_law / V3_both / V4_legacy 的
η_E、η_N 等）。

- ⚠️ **不得作论文/验证绝对数字**：权威网格为 `nz=250`，本目录仅用于**横向比较开关
  效应大小**；其 η_E/η_N 数值不可直接引用。
- 生成脚本：`scripts/entrypoints/run_hb_ablation.py`（`nz=60`，相对对照口径）。
- 口径与适用域声明见 `docs/源模型口径与适用域声明.md` **声明 3b**（`enable_power_law_gap_correction`
  为一阶近似、非 B&F25 数值闭包；开关空转 ≠ 物理中性）。
