# D2DGA T1 收尾 + T1-8 重定向 + T2-3 设计

> **日期**：2026-07-20 | **状态**：等待确认后执行  
> **覆盖**：S1 Bug修复 + S2 清理文档 + S3 T1-8 重定向 + S5 T2-3 TVD

**导入者**：无（设计文档，非代码文件）。**API 影响**：S1 修改 `_compute_velocity` 内部逻辑不改变签名；S3 修改 `_compute_dispersion_coefficient` 内部逻辑不改变签名。**数据**：不涉及数据文件。  
**用户指令**：`T2 缩减到只剩 T2-3，写一下计划然后直接执行，最后再跑一下目前模型的结果分析一下目前模型的合理性`

---

## S1: Bug 修复（2 处）

### Bug 1: I3 局部 CFL 裁剪
- **位置**：`annulus_d2dga.py:904`
- `step_limit = self.alpha_cfl * ds / max(self.dt, 1e-9)` → `max(dt_step, 1e-9)`

### Bug 2: base 未乘 I₁
- **位置**：`annulus_d2dga.py:618` 后插入 `base = base * i1`

---

## S2: 清理
1. 删 `correlations/`、`utils/` 空占位包
2. 修过期注释、记忆
3. Tier 0 诊断接入 runner（`run_and_export` 末尾调用 `compute_all_tier0_diagnostics`）

---

## S3: T1-8
| 牛顿 | 保留 |
| 幂律 | Batot 2016 式(28): `D_m[1+Pe²/192·24n²/((3n+1)(5n+1))]` |
| Bingham | Fan & Wang 1966: `k(ξ₀)` 显式多项式 |
| HB | 等效 Bingham 近似 |

---

## S4: T1-9 → 废弃

## S5: T2-3 TVD（后续单独 spec，本 spec 不做详细设计）

## 执行顺序: Bug1 → Bug2 → 清理 → T1-8 → 全量测试 → 模型分析
