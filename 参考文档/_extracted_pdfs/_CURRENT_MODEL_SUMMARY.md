# 当前模型概要（供子代理对照用）

## 项目
模拟**尾管（liner）固井**的环空顶替效率。5 层 1D-2D 耦合 **D2DGA** 模型。

## 5 层架构
1. **数据层** `cemdisp/data/`：WellSpec（井径/井斜/居中度剖面、尾管尺寸、鞋深）、FluidSpec（密度+流变，支持 newtonian/bingham/power_law/herschel_bulkley 四种）、PumpingSchedule（泵注序列，支持停泵/续泵/换排量 event_tag）。
2. **套管内 1D 输运** `cemdisp/transport1d/casing_flow.py`：**解析体积推进法**（无显式时间积分循环），`distance = cumulative_volume / pipe_area`。叠加：重力沉降修正（密度差+井斜投影+屈服抑制）、轴向弥散（牛顿 Taylor-Aris `D_eff=U²R²/(192·D_mol)`；幂律 `k=48n+144`；Bingham/HB `D_eff=α·U·R`，α=0.25）、停泵凝胶发展。
3. **鞋口桥接** `boundary_bridge.py`：1D 鞋口时序 → 环空入口相分数+排量函数。
4. **环空 2D D2DGA 核心** `cemdisp/models2d/annulus_d2dga.py` + `d2dga_flux.py`。
5. **诊断/验证/报告**。

## 第 4 层 2D D2DGA 核心现状（最关键）
- **控制方程**：四相**纯对流** `∂c_k/∂t + ∇⋅(u·c_k)=0`（k=lead,tail,spacer；mud=1-Σ）。扩散由 D2DGA 通量放大 + 小系数拉普拉斯平滑（轴向0.018/方位0.015）替代，**无显式扩散项**。
- **几何**：偏心环空间隙展开 `H(φ,s)=half_gap_mean·(1+e·cos(πφ))`，e=clip(1-standoff,0.05,0.55)，体积校正。
- **速度场**（Hele-Shaw 风格，**不解 Navier-Stokes**）：基础流动度 `base=(b/b_mean)²/μ_reg`（宽边主导）；浮力修正形状 `1+density_contrast·e·(2φ-1)`；轴向速度 `w=q_half·pref/Σ(pref·b·dy)`（截面排量约束）；横向速度 `v` 由连续性数值积分。
- **时间推进**：**半拉格朗日反演追踪**（bilinear 插值）。
- **流变处理**：`_apparent_viscosity` 对 HB 用 **Papanastasiou 屈服正则化**（M=100），`μ_app = (γ̇^n + (τ_y/M)·(1-exp(-M·γ̇)))/γ̇` 之类。**没有用 augmented Lagrangian 求真实屈服速度场**。
- **三闭包**（生产默认全开 R3）：
  - R1 auto-m：局部黏度比场 `m_field=μ_mud(γ̇)/max(μ_cement(γ̇),1e-6)`，裁剪[0.1,10]。
  - R2 I3 浮力弥散通量：`q_buoy=(Δρ·H³/(6η2))·I3(c̄,m)·[-f_xi,f_phi]`，I3=`c̄²(1-c̄)³[4m·c̄+3(1-c̄)]/{2m[m·c̄³+1-c̄³]}`，flux_strength=0.05 限幅。
  - R3 真浮力：局部密度场→density_contrast→buoyancy_shape→速度。**注意：体力向量 _buoyancy_force_vector 是死代码，未注入流动度**。
  - 通量放大 `f=[m·c²+1.5(1-c²)]/[m·c³+(1-c³)]`，c∈[0.01,0.99]，f∈[0.5,2.0]。
- **诊断指标**：顶替效率=`∫b·cement/half_volume`；窜槽=`|front_wide-front_narrow|/s_max`；混浆=`∫b·4c(1-c)/half_volume`；失稳=`1-exp(-channeling·(mobility_ratio-1)·(1+0.4·mixing)/decay_scale)`（启发式）。
- **数值**：固定 `dt=4s`，**无 CFL 自适应**；网格 ny=40(方位)×nz=500(轴向,生产)；通量放大限幅≤2.0；open_outlet=True。

## 已知功能缺口（代码现状）
- **壁面泥饼 `wall` 场恒为零**，未参与求解 → 无泥饼清除模型。
- **温度耦合参数被 `del`**（占位未接）。
- **凝胶强度参数被 `del`**（占位未接，仅 1D 有凝胶，2D 无）。
- **无显式 CFL**，固定 dt 靠限幅保稳。
- **质量守恒容差 25%**（偏松）。
- **1D 输运是解析的**，无瞬态时间步进循环。
- **HB 流体仅用表观黏度正则化**，非 D2DGA 论文的 HB 闭包。
- **体力向量死代码**，真浮力仅通过 density_contrast 形状因子间接作用。
- **无流动脉动建模**（1D 有停泵/续泵事件，2D 速度场按当前排量准静态重算，无瞬态惯性）。
- **失稳指数是启发式**，非论文的 Muskat/线性稳定性分析。
- **无指进（fingering）类型不稳定性预测**，无"静态壁面层/static wall layer"判据。

## 论文与模型的关系
当前模型的 D2DGA 理论源自 **Zhang & Frigaard (2022) JFM 947 A32**。子代理需对照所读论文，指出当前模型**遗漏或简化了论文中的哪些物理/数学机制**，以及**如何改进**。
