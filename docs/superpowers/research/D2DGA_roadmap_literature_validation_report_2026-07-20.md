# D2DGA 优化路线图 文献依据与方向合理性论证报告

> **日期**：2026-07-20 | **来源数**：5 篇本地 PDF + Tavily 联网检索 | **置信度**：中

---

## Executive Summary

- **T1-1 ~ T1-5**：全部有 **Zhang & Frigaard (2022) JFM 947 A32** 的直接公式支撑，方向正确，不构成反向优化。
- **T1-6**：**Yang 等 (2021) JPSE 204** 的 CFD 实验支撑（四相效率 89.67% vs 三相 85.88%），但当前被动平流实现无法复现论文中的流体力学交互效应。**无反向优化风险**。
- **T1-7**：半拉格朗日 CFL 自适应属于数值优化，物理合理。
- **T1-8**：权威来源 **Maleki & Frigaard (2016) JNNFM 235, 1–19** 已定位，但针对**弱湍流**。当前层流经验式 `α·U·R` 不可直接删除。建议先确定固井管流 Re 区间再做决策。
- **T1-9**：**未在任何已发表 D2DGA 论文中找到直接论述**，属于工程扩展。建议降级为实验性开关。

---

## 1. 文献支持总表

| 任务 | 方向 | 支持论文 | 支持强度 |
|------|------|----------|----------|
| T1-1 | 去 [0.5,2.0] 限幅 | Zhang & Frigaard 2022 式 4.28 | **强** |
| T1-2 | 物理 I3 系数 | Zhang & Frigaard 2022 式 4.25–4.26 | **强** |
| T1-3 | 体力向量 + I1/I2 | Zhang & Frigaard 2022 式 4.14–4.24 | **强** |
| T1-4 | 两层黏度闭包 | Zhang & Frigaard 2022 式 4.23 | **强** |
| T1-5 | static wall c_min | Bararpour & Frigaard 2025 式 2.35–2.41 | **中** |
| T1-6 | FLUSHER 独立相 | Yang 等 2021, JPSE 204 | **中** |
| T1-7 | CFL 自适应 | 数值 PDE 通识 | **中** |
| T1-8 | 屈服 Taylor 弥散 | **Maleki & Frigaard 2016 JNNFM 235** | **弱—中**（弱湍流） |
| T1-9 | 入口方位加权 | **无** | **无** |

---

## 2. T1-8 关键发现

通过 Tavily 联网检索定位到的 JNNFM 2016 论文：

> **Maleki, A. & Frigaard, I.A.** (2016). Axial dispersion in weakly turbulent flows of yield stress fluids. *J. Non-Newtonian Fluid Mech.*, 235, 1–19. doi: `10.1016/j.jnnfm.2016.07.002`.

**核心结论**：
- 针对**弱湍流管流**（5000 < Re < 10000），D_eff 可比高湍流大 O(10) 倍。
- 论文明确声明"for laminar flows primary cementing does not typically fall into the Taylor-regime"——**层流固井不在 Taylor 弥散区间**。
- 若当前固井管流多数在层流区间，用此论文公式反而可能**反向优化**。

**⚠️ 不可直接删除当前 `α·U·R` 经验式**。建议先确认固井管流典型 Re 区间后，再决定是否引入 Maleki 2016 作为弱湍流补充。

---

## 3. T1-9 判定

在 Zhang 2022、Bararpour 2025、Yang 2021 全文搜索均**未找到**入口按偏心间隙宽度方位加权的直接论述。D2DGA 的 gap-averaging 假设本身就不区分子方向的入口分布。方位加权可能与框架矛盾。建议**降级为实验性开关**，默认关闭。

---

## 4. 实现简化风险评估

| 差异 | 反向优化风险 | 建议 |
|------|-------------|------|
| T1-3b 标量 m 代替局部 m 场 | 低 | 保留，标注为工程简化 |
| T1-3b 未对 base 乘 I₁ | 低 | 保留 η_mix 闭包，Tier 2 统一 |
| I3 局部 CFL 用 self.dt 非 dt_step | 极低 | 保守，可后续修正 |
| FLUSHER 被动平流 | 无 | 合规，spec 已接受 |
| c_min=0.05 固定 | 低 | 标注为可调工程参数 |

---

## 5. 建议优先级

| 优先级 | 行动 |
|--------|------|
| P0 | T1-7 收尾 + 清理空占位 + Tier 0 runner 接线 |
| P1 | T1-8 前置研究：确定固井管流 Re 区间 + 获取 Maleki 2016 全文公式 |
| P2 | T1-8 实现（仅当前置研究支持时） |
| P3 | T1-9（实验性开关，默认关闭） |

---

## 6. 核心参考文献

1. **Zhang & Frigaard (2022)** — Primary cementing of vertical wells: displacement and dispersion effects in narrow eccentric annuli. *J. Fluid Mech.*, 947, A32.
2. **Zhang & Frigaard (2023)** — Part 2. Flow behaviour and classification. *J. Fluid Mech.*, 972, A38.
3. **Bararpour & Frigaard (2025)** — Capturing dispersion of Herschel–Bulkley fluids in miscible primary cementing displacement flows. *J. Fluid Mech.*, 1022, A15.
4. **Yang et al. (2021)** — Improving displacement efficiency by optimizing pad fluid injection sequence... *J. Pet. Sci. Eng.*, 204, 108691.
5. **Maleki & Frigaard (2016)** — Axial dispersion in weakly turbulent flows of yield stress fluids. *J. Non-Newtonian Fluid Mech.*, 235, 1–19.
6. **Pelipenko & Frigaard (2004)** — Visco-plastic fluid displacements in near-vertical narrow eccentric annuli. *J. Fluid Mech.*, 520.
7. **Moyers-González et al. (2007)** — Transient effects in oilfield cementing flows. *Euro. J. Appl. Math.*, 18.
