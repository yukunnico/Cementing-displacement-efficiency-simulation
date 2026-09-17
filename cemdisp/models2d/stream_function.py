"""D2DGA 流函数椭圆方程求解（Zhang & Frigaard 2022, (2.3)/(4.22)）。

方程与文献锚点
--------------
解 (4.22)（与 §2 概览式 (2.3) 同型；Newton 两层下 S 线化为 (r_a/2I₁)∇aΨ）::

    ∇a·[ (r_a/(2I₁))·∇aΨ + b ] = 0，    ∇a = ( (1/r_a)∂φ, ∂ξ )    （Z&F22 (2.4)）

其中 I₁ 为两层牛顿流动度闭包（Z&F22 (4.21a)，经
:func:`cemdisp.models2d.two_layer.mobility_i1` 消费，本模块不重写闭式）。
关键性质：每个时间步 I₁、I₂ 只依赖 c̄ 与 m（均为已知场），算子矩阵与
浮力场无关 ⇒ 方程对 Ψ **线性**——每步只需解一次变系数线性 Poisson
（``scipy.sparse.linalg.spsolve`` 直接解）。论文 §2.1 的增广拉格朗日
迭代只用于 HB 流体 S(|∇aΨ|) 的非线性，本模块**不实现**该套算法。

H 口径核实结论（Task 6 遗留顾虑，2026-09-15 以源文献逐句核实）
--------------------------------------------------------------
**Z&F22 的 H = 半间隙（scaled half-gap-width），与本仓 geom["H"] =
(井径−外径)/4（Task 3/R20 钉定）同口径，不存在 8/16 倍速度换算风险。**
证据（Z&F22 JFM 947 A32, doi:10.1017/jfm.2022.626）：

1. §2.1（A32-5 页）："the mean annular gap (**2d̂**) is much narrower
   compared with both the mean circumference (2πr̂ₐ*)..."——平均间隙
   写作 2d̂，即 d̂ 为半间隙；
2. 同页末跨页句（A32-5→A32-6）："The scaled **half-gap-width** is
   H(φ, ξ) and the mean radius at depth ξ is r_a(ξ)."——H 的原文定义
   就是（无量纲化后的）**半隙宽度**；
3. 因子自洽：(2.2) 的 ∂φΨ = 2Hr_a·w̄ 中因子 2 恰把半隙 H 复原为全隙
   2H 的轴向通量（弧长 r_a dφ × 全隙 2H × w̄）；单流体极限
   I₁ = H³/(3η)（(4.21a) 取 y_i=0）对应半域 y∈[0,H]（y=0 对称面、
   y=H 壁面，(4.19) 分部积分要求 u(H)=0）的半通量——Task 6 实现者
   "H³/3η 对应 2H 全隙槽流的半域" 的观察正确，而结论是 H 本身就是
   半隙，故全部因子一致，无需换算。

坐标、网格与 r_a 约定
---------------------
- φ∈[0,1]：半环空方位坐标（图 2 的 Hele-Shaw 展开：πφ∈[0,π]，φ=0 低边
  → φ=1 高边；(2.5b) 方位重力投影 sin πφ 与之配套），均匀网格；
- ξ：轴向深度（取 ``geom["s"]``），均匀网格；
- ``r_a``：环空平均半径 = ``mean((hole_mm+od_mm)/4)/1000`` 米——与
  ``annulus_d2dga._buoyancy_force_vector``/``buoyancy.froude_squared``
  的 r̂ₐ* 口径一致（(2.4) 的 1/r_a 度量因子与 (2.2) 的速度换算均用它）。
  边界声明：r_a 取**沿深度平均的常数标量**（既有 r̂ₐ* 口径）；论文的
  r_a(ξ) 逐深度量本模块不支持——若未来需要，(2.4) 的 1/r_a 因子与
  (2.2) 速度换算须同步改造。

b_field 约定（controller 裁定 R29，2026-09-15；T9 按此派发）
--------------------------------------------------------------
``b_field`` 是 (4.22) 浮力向量场的**无量纲**入口，**唯一口径**为
**(2, ny, nz) 全向量、(4.22) 字面分组**：

- ``b_field[0] = b_φ``（φ-散度槽）、``b_field[1] = b_ξ``（ξ-散度槽），
  源项 = −(1/r_a)∂φ b_φ − ∂ξ b_ξ（面值取相邻均值，ξ 端一阶单侧）；
- 字面分组由本模块实现者从 (4.14)+(2.2) 交叉微分兼容条件独立复推核实：
  φ-槽配 cos β（轴向重力，竖直井机制——源项需 χ 的 φ-梯度，见 A32-22
  页"the elliptic equation for the stream function is driven by
  gradients in b..."）、ξ-槽配 sin πφ·sin β（方位重力，斜井前缘机制，
  本仓井斜 1.8–15° 为 O(sinβ) 二阶）；
- 竖直井装配示意：``b_φ = χ·r_a·cosβ/F²``、``b_ξ = χ·r_a·sin(πφ)·sinβ/F²``，
  ``χ = ρ − Δρ·I₂/(H·I₁)``（等黏度时 I₂/(H·I₁) = c̄(1−c̄²)/2——因子 2 由
  Z&F22 (4.24) 取 m=1 与 two_layer 闭式 I₂/H⁴÷I₁H/H⁴ 双证，T9 订正，
  2026-09-15）；密度均匀（c̄ 的 φ-梯度为零）时源项恒为零——这是论文
  自身的机制属性，不是缺陷。

⚠️ **均匀幅值 + 字面分组 ⇒ 源项为零**：均匀（或任一散度为零）的 b 场
不改变 Ψ——调用方（T9）必须让 b 携带 χ 的 φ-梯度（竖直井）或 ξ-梯度
（斜井前缘）才有浮力响应。

📜 R29 裁定记录：初版曾提供 (ny,nz) 标量入口并把 M·sin(πφ) 组装进
φ-槽（标为"默认/推荐"）。审查指出均匀 M 下该组装产生人为强迫
−Mπ·cos(πφ)/r_a——不是论文任何机制（竖直井机制要求 χ 的 φ-梯度，均匀
c̄ 下应为零；A32-20 页明言 f 的梯度仅来自缓变井向）。controller 裁定
取源文献忠实（R29，2026-09-15）：(4.22) 字面分组为唯一口径，该标量
路径**已删除**（无任何剩余消费方，YAGNI）。

数值格式（controller 裁定第 3 条）
----------------------------------
变系数线性 Poisson，5 点格式，行主序编号 ``idx = i*nz + j``：

- φ-槽系数 1/(2I₁)、ξ-槽系数 r_a/(2I₁)（(4.22) 展开后的面通量系数），
  面值取相邻格点的**调和平均**（面通量守恒）；
- φ 两端 Dirichlet：Ψ(φ=0)=0、Ψ(φ=1)=1——单位通量（Ψ∈[0,1] 归一），
  调用方按实际排量 Q 缩放；
- ξ 两端零梯度（一阶单侧，ghost=本格 ⇒ 边界面通量为零）；
- ``scipy.sparse.linalg.spsolve`` 直接求解。

速度换算（Z&F22 (2.2) 逐字）
----------------------------
``velocity_from_stream_function`` 返回 ``(w, v)``::

    w̄ = ∂φΨ/(2 r_a H)（轴向）、    v̄ = −∂ξΨ/(2H)（方位）

梯度内点中心差分、边界一阶单侧。该组合下
``trapezoid(2H·w̄, x=φ)`` 逐列严格等于 (Ψ(1)−Ψ(0))/r_a = 1/r_a
（差分-梯形恒等式，与 b_field 无关）——浮力只在列内重分配轴向通量，
不改变每列总通量。

**Q 缩放锚（不变量）**：单位 BC（Ψ∈[0,1]）下每列满足
``∫₀¹ 2·r_a·H·w̄ dφ = 1``（半环空轴向体积通量 = 1）。调用方把 Ψ（或
w、v）整体乘以**半环空物理通量 Q_half = Q/2**（Q 为全环空排量）即得
物理场——缩放后不变量变为 Q_half。禁止按其它倍数缩放（初版 docstring
的"按 Q/0.5 缩放"表述已被 R29 修复轮替换为本不变量锚，避免双倍计因子）。

⚠️ **T9 度量换算订正（2026-09-15，数值双验证）**：上段"∫2r_aH·w̄ dφ =
半环空轴向体积通量"是 **φ-度量**口径；模型输运/体积层用**弧长度量**
（geom["y"] = 半周长弧、geom["b"] 经体积标定），同一速度场的弧长列通量
``∫w·b·dy = π·∫₀¹2r_aH·w̄ dφ``（弧元 dy = π·r_a·dφ），数值验证逐列比恰
为 π。故消费方（``annulus_d2dga._velocity_stream_function``，T9）取

    w = w̄_unit·(Q/2)/π，  v = v̄_unit·(Q/2)，

其中 v 多乘 π 源于模块 (w,v) 对满足 ``π·∂y(bv) + ∂s(bw) = 0``（论文 ξ 以
πr̂ₐ* 标定、模块 ξ 用米）而非物理连续性 ``∂y(bv)+∂s(bw) = 0``——v×(Q/2)
后物理连续性残差 ~8e-15（.tmp_research/task9_probe/verify_pi_scaling.py
[1]-[5] 段，2026-09-15）。本模块代码不变，仅订正本契约表述。

非线性外迭代（A-3a，Task 5）
----------------------------
HB 流体下 ``I₁`` 依赖局部应力尺度 ⇒ (4.22) 不再线性，须外迭代：
:func:`solve_stream_function_nonlinear` 以"冻结 Ī₁ → 解线性 Ψ → 由 Ψ 取局部均速 →
逐格反求 G（闭包 (2.13) 的反演）→ 更新 Ī₁ → 欠松弛 → 重解"的循环逼近不动点。
**线性路径（:func:`solve_stream_function`）一字未动**，外迭代首轮即原样复用它
（R2）：牛顿极限下 Ī₁ 与 G 无关 ⇒ 首轮检测到 ``ΔĪ₁ ≡ 0`` 即刻返回该首轮 Ψ，
与线性调用**逐位一致**（硬验收）。

⚠️ **标度口径（单一标度；由接线方对齐，见 task-5-report.md §7.1 订正版——其中偏差
方向以 fix round 2 生产锚实测为准）**：外迭代的
反求与闭包求值都在本模块的**单位通量 Ψ 口径**内自洽（``ū`` 由 (2.2) 同口径给出、
``G`` 由 ``H·ū = I₁·G`` 反算）。闭包关系的**标度不变性**（速度缩放 ``s``、应力缩放
``σ``；本仓数值核验 9 组 rel ≤ 1.2e-14，见 ``tests/contract/``
``test_flow_curve_scale_covariance``）::

    κ' = κ·sⁿ/σ,   τ_Y' = τ_Y/σ,   ū' = ū/s  ⇒  G' = G/σ,  I₁' = (σ/s)·I₁

取 **``σ = s``**（同一物理问题换单位：长度不动 ⇒ 压力梯度与速度同步缩放）⇒
**``I₁' = I₁``：``I₁`` 对标度不变 ⇒ 算子与反求要的是同一个 ``I₁``，不存在冲突**。
⇒ **只有单一自由标度**（生产装配的 ``ŵ = q_half/π``，T9 推导）；接线只需让"喂进反求的
``ū``"与"闭包的 ``κ``/``τ_Y`` 口径"使用**同一**标度：

- **推荐**：``ū`` 乘 ``ŵ``（物理速度）+ 闭包用物理 ``(κ, τ_Y)``——该乘法由
  ``solve_stream_function_nonlinear`` 的 ``velocity_scale`` 形参承载（Task 6 生产
  接线传 ``ŵ = q_half/π``）；或
- 等价地：闭包改用模块口径参数 ``(κ·ŵ^{n−1}, τ_Y/ŵ)``、``ū`` 保持模块口径
  ——二者给出**同一个** ``I₁``（不变性）。模块口径速度比物理**大** ``1/ŵ`` 倍
  （生产锚 ``annulus_d2dga.py:1433``：``w = w_unit·(q_half/π)``），故模块口径的
  应力与 ``τ_Y`` 同比大 ``1/ŵ`` 倍、``κ`` 大 ``ŵ^{n−1}`` 倍——无量纲群不变。
  （⚠️ fix round 1 曾误写 ``(κ·ŵ^{1−n}, τ_Y·ŵ)``、并以反向口径探针"证伪"生产
  约定；fix round 2 以生产锚实测订正（正确式 rel ≈ 1e-14，旧式 rel ≈ 40），
  由 ``test_velocity_scale_alignment_equivalence`` 按生产口径钉住。）

**口径不一致的后果（方向以生产代码为锚实测，勿再写反）**：模块单位通量口径的
``ū`` 比物理**大** ``1/ŵ`` 倍（呼101 量级实测：模块 ū 中位 86.9208、物理速度
中位 0.230601 m/s、``ŵ = q_half/π = 2.653e-3``，1/ŵ = 376.9，探针
probe_direction_chain.py [1]-[2]）。若闭包用物理 ``(κ, τ_Y)`` 而 ``ū`` 未乘
``ŵ``（喂的是模块口径速度），闭包看到的应力比物理**大** ``ŵ^{-n}`` 倍（呼101
量级实测 ≈ 53×：错接线壁面应力中位 332.50 Pa vs 正确接线 6.29 Pa）⇒ 屈服门槛
``τ_Y`` 形同虚设 ⇒ 几乎无格点落 ``undefined`` ⇒ **屈服效应被削弱（weakened）**；
反之若某工况 ``ŵ > 1``（模块口径速度反而小于物理），方向反转（应力偏小 ⇒
静止区被夸大 ⇒ 屈服效应被增强）。
"""

from __future__ import annotations

from typing import Dict, Tuple

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numpy.typing import NDArray

from cemdisp.models2d.gap_solver import solve_g_from_mean_velocity_batch
from cemdisp.models2d.hb_closure import ClosureProvider, NewtonianClosure

Array = NDArray[np.float64]

__all__ = ["solve_stream_function", "solve_stream_function_nonlinear",
           "velocity_from_stream_function"]

# B-2：屈服门冻结区的电导地板——防止 wall≡1 时系数全零导致矩阵奇异。
_WALL_CONDUCTANCE_FLOOR = 1.0e-6

# A-3a 外迭代的量级常量（与 gap_solver 的求根容差/上限；见其 docstring）
_INVERSE_RTOL = 1.0e-12
_INVERSE_MAX_ITER = 200

# (2.2) 两分量的物理换算因子之比：w_phys ∝ w̄/π、v_phys ∝ v̄（T9 π 推导）⇒ v̄ 需乘 π
# 才与 w̄ 同口径（Important-4；本仓 |v|/|w| ~ 1e-6，数值影响极小，属正确性硬化）。
_VELOCITY_COMPONENT_RATIO = float(np.pi)


def _mean_radius_m(geom: Dict) -> float:
    """环空平均半径 r_a（米）= mean((hole_mm+od_mm)/4)/1000。

    与 ``annulus_d2dga._buoyancy_force_vector`` 的 r̂ₐ* 口径一致：
    (hole+od)/4 = (r_o+r_i)/2 为环空平均半径，沿深度取均值。
    """
    if "hole_mm" not in geom or "od_mm" not in geom:
        raise ValueError("geom 缺少 hole_mm/od_mm，无法确定环空平均半径 r_a（Z&F22 (2.4)）")
    hole = np.asarray(geom["hole_mm"], dtype=float)
    od = np.asarray(geom["od_mm"], dtype=float)
    return float(np.mean((hole + od) / 4.0)) / 1000.0


def _uniform_spacing(x, name: str) -> float:
    """均匀网格步长；非均匀网格直接报错（5 点格式按均匀网格推导，本仓 2D 网格为均匀 linspace）。"""
    arr = np.asarray(x, dtype=float).reshape(-1)
    if arr.size < 2:
        raise ValueError(f"geom[{name!r}] 至少需要 2 个网格点")
    d = (arr[-1] - arr[0]) / (arr.size - 1)
    if not np.allclose(np.diff(arr), d, rtol=1e-6, atol=1e-14):
        raise ValueError(f"geom[{name!r}] 网格必须均匀（本模块 5 点格式按均匀网格推导）")
    return float(d)


def _scalar_viscosity(eta, name: str) -> float:
    """两层闭包前提：每流体单一黏度（Z&F22 (4.21)）。数组输入必须全场一致，禁止静默平均。"""
    arr = np.asarray(eta, dtype=float)
    if arr.ndim == 0:
        return float(arr)
    first = float(arr.reshape(-1)[0])
    if not np.allclose(arr, first, rtol=1e-12, atol=0.0):
        raise ValueError(f"{name} 必须为标量或全场一致（两层闭包 (4.21) 以每流体单一黏度为前提）")
    return first


def solve_stream_function(geom: Dict, c_bar, eta1, eta2, m: float, b_field,
                          *, closure: ClosureProvider | None = None,
                          wall: Array | None = None,
                          ny=None, nz=None) -> Array:
    """解 (4.22) 流函数椭圆方程，返回 Ψ 场 ``(ny, nz)``（单位通量口径）。

    Args:
        geom: 几何字典，需含 ``phi``(ny,)、``H``(ny,nz)（Z&F22 半隙口径，
            见模块 docstring 的核实结论）、``s``(nz,)、``hole_mm``/``od_mm``。
        c_bar: (ny,nz) 间隙平均水泥体积分数 c̄ = y_i/H ∈ [0,1]。
        eta1: 流体 1（被顶替液，壁面带）黏度，标量或全场一致数组。
        eta2: 流体 2（顶替液，中线带）黏度，标量或全场一致数组。
        m: 黏度比 η₁/η₂ > 0。
        b_field: (2,ny,nz) 浮力全向量（(4.22) 字面分组，唯一口径）：
            b_field[0]=b_φ（φ-散度槽）、b_field[1]=b_ξ（ξ-散度槽），无量纲；
            竖直井装配 b_φ = χ·r_a·cosβ/F²、b_ξ = χ·r_a·sin(πφ)·sinβ/F²。
            均匀幅值输入源项为零——浮力响应须由 χ 的 φ-/ξ-梯度携带，
            约定与 R29 裁定记录见模块 docstring。
        closure: 闭包提供者（B-1，``ClosureProvider``）。``None`` ⇒
            ``NewtonianClosure()``，逐位等于 HEAD 的 ``two_layer.mobility_i1``。
            Phase A 的 HB 查表闭包实现同一协议后从此注入。
        wall: (ny,nz) 屈服门冻结度 ∈ [0,1]（B-2，Pelipenko04 (2.6)-(2.8) 口径）。
            ``None`` ⇒ 不施加（逐位 = HEAD）。非 ``None`` 时把有效流动度取为
            ``I₁_eff = I₁·max(1−wall, _WALL_CONDUCTANCE_FLOOR)``（算子 cell 系数
            1/(2I₁_eff) 随之放大，地板防止 wall≡1 时矩阵奇异）：冻结区流动度→0，
            Ψ 局部趋于常数 ⇒ 该处速度→0（static wall layer），且保持椭圆结构。
            注意方向：I₁ 是流动度、1/(2I₁) 是阻力——(1−wall) 必须作用在 I₁ 上。
            冻结不改变总通量（BC 单位通量口径仍是 Ψ(1)−Ψ(0)=1），只把通量从冻结
            区重新分配到活跃区（窄边冻结 ⇒ 宽边流速上升）。
        ny: 可选形状自校验（与 geom["H"].shape[0] 不符即抛错）。
        nz: 可选形状自校验（与 geom["H"].shape[1] 不符即抛错）。

    Returns:
        Ψ 场 (ny, nz)，float64；Ψ(φ=0)=0、Ψ(φ=1)=1（单位通量，调用方按 Q 缩放）。
    """
    H = np.asarray(geom["H"], dtype=float)
    if H.ndim != 2:
        raise ValueError("geom['H'] 必须是 (ny,nz) 二维数组")
    ny_, nz_ = H.shape
    if ny is not None and int(ny) != ny_:
        raise ValueError(f"ny={ny} 与 geom['H'].shape[0]={ny_} 不符")
    if nz is not None and int(nz) != nz_:
        raise ValueError(f"nz={nz} 与 geom['H'].shape[1]={nz_} 不符")
    if ny_ < 3 or nz_ < 2:
        raise ValueError("网格过细：需 ny>=3（φ 两端 Dirichlet + 至少 1 个内点）、nz>=2")
    if not np.all(H > 0.0):
        raise ValueError("geom['H'] 必须处处为正（Z&F22 半隙宽度）")

    if "phi" not in geom:
        raise ValueError("geom 缺少 phi（半环空方位坐标 (ny,)，Z&F22 图 2 口径）")
    phi = np.asarray(geom["phi"], dtype=float).reshape(-1)
    if phi.size != ny_:
        raise ValueError("geom['phi'] 长度须与 geom['H'].shape[0] 一致")
    dphi = _uniform_spacing(phi, "phi")
    if "s" not in geom:
        raise ValueError("geom 缺少 s（轴向深度坐标 (nz,)）")
    dxi = _uniform_spacing(geom["s"], "s")
    r_a = _mean_radius_m(geom)

    c = np.asarray(c_bar, dtype=float)
    if c.shape != H.shape:
        raise ValueError("c_bar 形状须与 geom['H'] 相同 (ny,nz)")
    e1 = _scalar_viscosity(eta1, "eta1")
    e2 = _scalar_viscosity(eta2, "eta2")
    m = float(m)
    if not m > 0.0:
        raise ValueError("m = η₁/η₂ 必须为正")

    # I₁ 闭包（Z&F22 (4.21a)；c̄/H 传场、η 传标量）。B-1：经 ClosureProvider 注入；
    # 缺省 NewtonianClosure 逐位等于 two_layer.mobility_i1（L1 硬约束）。
    if closure is None:
        closure = NewtonianClosure()
    I1 = np.asarray(closure.mobility(c, m, e1, e2, H), dtype=float)

    # ---- b 场入口（(4.22) 字面分组全向量，唯一口径，R29）----------------
    b = np.asarray(b_field, dtype=float)
    if b.shape != (2, ny_, nz_):
        raise ValueError(
            "b_field 形状须为 (2,ny,nz)——(4.22) 字面分组全向量（b_field[0]=b_φ、"
            "b_field[1]=b_ξ，唯一忠实口径，controller 裁定 R29）；"
            "(ny,nz) 标量路径已删除（均匀幅值+字面分组下源项本应为零，"
            "浮力响应须由 χ 的 φ-/ξ-梯度经 b 场携带，见模块 docstring）"
        )
    b_phi, b_xi = b[0], b[1]
    bphi_face = 0.5 * (b_phi[:-1] + b_phi[1:])

    # ---- 面迁移率（调和平均，面通量守恒）--------------------------------
    # B-2：屈服门冻结度 → 有效流动度 I₁_eff = I₁·max(1−wall, 地板)。
    # 物理：冻结层不可流动 ⇔ 该处流动度 I₁→0。Z&F22 (4.22) 的 S_φ = ∂φΨ/(2I₁)
    # = −r_a·∂ξp 为压力梯度、轴向流密度 ∂φΨ/r_a = 2·I₁·S_φ ∝ I₁ ⇒ I₁→0 处流密度
    # →0、Ψ 趋于常数 ⇒ w = ∂φΨ/(2r_aH) → 0（static wall layer），且椭圆结构不变。
    # ⚠️ 方向警示：a_cell = 1/(2I₁) 是 **阻力的倒数**——若把 (1−wall) 乘到 a_cell
    # 上（brief 字面式），冻结区反而变成 Ψ-方程的高阻区，Ψ 落差被挤进冻结带、
    # 该处速度不降反升（与窄边冻结物理及本设计验收判据「窄边速度下降」相反）。
    # 实测对照见 .tmp_research/task2_probe/probe_wall_direction.py。
    if wall is None:
        I1_eff = I1                       # 逐位无扰（不参与任何算术）
    else:
        w_arr = np.asarray(wall, dtype=float)
        if w_arr.shape != H.shape:
            raise ValueError("wall 形状须与 geom['H'] 相同 (ny,nz)")
        conductance = np.maximum(1.0 - w_arr, _WALL_CONDUCTANCE_FLOOR)
        I1_eff = I1 * conductance
    a_cell = 1.0 / (2.0 * I1_eff)        # φ-槽系数 1/(2I₁_eff)
    c_cell = r_a / (2.0 * I1_eff)        # ξ-槽系数 r_a/(2I₁_eff)
    a_face = 2.0 * a_cell[:-1] * a_cell[1:] / (a_cell[:-1] + a_cell[1:])       # (ny-1, nz)
    c_face = 2.0 * c_cell[:, :-1] * c_cell[:, 1:] / (c_cell[:, :-1] + c_cell[:, 1:])  # (ny, nz-1)

    # ---- 右端 −∇a·b = −[(1/r_a)∂φ b_φ + ∂ξ b_ξ]（面散度）----------------
    src = np.zeros((ny_, nz_))
    src[1:-1] += (bphi_face[1:ny_ - 1] - bphi_face[0:ny_ - 2]) / (r_a * dphi)
    bxi_face = 0.5 * (b_xi[:, :-1] + b_xi[:, 1:])    # 内部 ξ 面值 (ny, nz-1)
    flux_plus = np.zeros((ny_, nz_))                 # j+½ 面浮力通量
    flux_plus[:, :nz_ - 1] = bxi_face
    flux_plus[:, -1] = b_xi[:, -1]                   # ξ 上端一阶单侧（ghost=本格）
    flux_minus = np.zeros((ny_, nz_))                # j-½ 面浮力通量
    flux_minus[:, 1:] = bxi_face
    flux_minus[:, 0] = b_xi[:, 0]                    # ξ 下端一阶单侧（ghost=本格）
    src += (flux_plus - flux_minus) / dxi

    # ---- 5 点矩阵装配（行主序 idx = i*nz + j，全网格含 Dirichlet 行）----
    n = ny_ * nz_
    i_int = np.arange(1, ny_ - 1)
    ni = i_int.size
    II, JJ = np.meshgrid(i_int, np.arange(nz_), indexing="ij")   # (ni, nz)
    row_int = (II * nz_ + JJ).ravel()

    east_val = a_face[1:ny_ - 1] / (r_a * dphi * dphi)   # 面 i（i↔i+1）
    west_val = a_face[0:ny_ - 2] / (r_a * dphi * dphi)   # 面 i-1
    east_col = ((II + 1) * nz_ + JJ).ravel()
    west_col = ((II - 1) * nz_ + JJ).ravel()

    north_val = np.zeros((ni, nz_))
    north_val[:, :nz_ - 1] = c_face[1:ny_ - 1] / (dxi * dxi)
    south_val = np.zeros((ni, nz_))
    south_val[:, 1:] = c_face[1:ny_ - 1] / (dxi * dxi)
    # 越界邻居（j=0 无南、j=nz-1 无北）值已为 0；列号钳到自身，COO 累加 0 无害。
    north_col = (II * nz_ + np.minimum(JJ + 1, nz_ - 1)).ravel()
    south_col = (II * nz_ + np.maximum(JJ - 1, 0)).ravel()
    center_val = -(east_val + west_val) - (north_val + south_val)

    # Dirichlet 行（i=0 与 i=ny-1 的全部 j）：对角 1.0，BC 值进右端
    i_bc = np.concatenate([np.zeros(nz_, dtype=int), np.full(nz_, ny_ - 1, dtype=int)])
    j_bc = np.tile(np.arange(nz_), 2)
    row_bc = i_bc * nz_ + j_bc
    val_bc = np.where(i_bc == 0, 0.0, 1.0)   # BC 值：Ψ(φ=0)=0、Ψ(φ=1)=1（进右端）

    rows = np.concatenate([row_int, row_int, row_int, row_int, row_int, row_bc])
    cols = np.concatenate([east_col, west_col, row_int, north_col, south_col, row_bc])
    vals = np.concatenate([east_val.ravel(), west_val.ravel(), center_val.ravel(),
                           north_val.ravel(), south_val.ravel(), np.ones(2 * nz_)])
    matrix = sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()

    rhs = np.zeros(n)
    rhs[row_bc] = val_bc
    rhs[row_int] = -src[1:-1].ravel()

    psi = spla.spsolve(matrix, rhs)
    psi = np.asarray(psi, dtype=float).reshape(ny_, nz_)
    if not np.all(np.isfinite(psi)):
        raise RuntimeError("流函数方程求解失败（出现非有限值），请检查输入场的量级与正则性")
    return psi


def velocity_from_stream_function(psi, geom: Dict) -> Tuple[Array, Array]:
    """由 Ψ 场换算间隙平均速度（Z&F22 (2.2) 逐字）：返回 ``(w, v)``。

    ::

        w̄ = ∂φΨ/(2 r_a H)   （轴向）
        v̄ = −∂ξΨ/(2H)       （方位）

    梯度内点中心差分、边界一阶单侧。该组合下 trapezoid(2H·w̄, x=φ) 逐列
    严格等于 (Ψ(1)−Ψ(0))/r_a（差分-梯形恒等式），与浮力场无关。

    Args:
        psi: (ny, nz) 流函数场（``solve_stream_function`` 的返回值）。
        geom: 与求解时同一几何字典。

    Returns:
        ``(w, v)``：轴向与方位间隙平均速度场，各 (ny, nz)，
        "单位通量"口径（调用方按实际排量 Q 缩放）。
    """
    psi = np.asarray(psi, dtype=float)
    H = np.asarray(geom["H"], dtype=float)
    if psi.shape != H.shape:
        raise ValueError("psi 形状须与 geom['H'] 相同 (ny,nz)")
    phi = np.asarray(geom["phi"], dtype=float).reshape(-1)
    if phi.size != psi.shape[0]:
        raise ValueError("geom['phi'] 长度须与 psi.shape[0] 一致")
    dphi = _uniform_spacing(phi, "phi")
    dxi = _uniform_spacing(geom["s"], "s")
    r_a = _mean_radius_m(geom)

    # ∂φΨ：内点中心差分，φ 两端一阶单侧
    dpsi_dphi = np.empty_like(psi)
    dpsi_dphi[1:-1] = (psi[2:] - psi[:-2]) / (2.0 * dphi)
    dpsi_dphi[0] = (psi[1] - psi[0]) / dphi
    dpsi_dphi[-1] = (psi[-1] - psi[-2]) / dphi
    w = dpsi_dphi / (2.0 * r_a * H)          # (2.2) 第一式：轴向

    # ∂ξΨ：内点中心差分，ξ 两端一阶单侧
    dpsi_dxi = np.empty_like(psi)
    dpsi_dxi[:, 1:-1] = (psi[:, 2:] - psi[:, :-2]) / (2.0 * dxi)
    dpsi_dxi[:, 0] = (psi[:, 1] - psi[:, 0]) / dxi
    dpsi_dxi[:, -1] = (psi[:, -1] - psi[:, -2]) / dxi
    v = -dpsi_dxi / (2.0 * H)                # (2.2) 第二式：方位
    return w, v


# --------------------------------------------------------------------------- #
# A-3a（Task 5）：HB 闭包的非线性外迭代
# --------------------------------------------------------------------------- #


class _FrozenMobilityClosure:
    """把"冻结的 I₁ 场"适配到 :class:`ClosureProvider`（外迭代的线性化步）。

    外迭代把 Ī₁ 冻结（上一轮欠松弛后的场）后解一次变系数线性 Poisson；本适配器
    就是那个冻结场：``mobility`` 逐字返回缓存场（**量纲** I₁，与求解器消费的口径
    相同），``buoyant_mobility`` 委托底层闭包（(4.22) 的算子只消费 I₁，I₂ 走
    闭包自身的 χ/b 通路，不经本适配器）。
    """

    def __init__(self, I1_field: Array, base: ClosureProvider) -> None:
        self._I1 = np.asarray(I1_field, dtype=float)
        self._base = base

    def mobility(self, c_bar, m: float, eta1: float, eta2: float, H) -> Array:
        out_shape = np.broadcast_shapes(np.shape(c_bar), np.shape(H))
        if self._I1.shape != out_shape:
            raise ValueError(
                f"冻结流动度场形状 {self._I1.shape} 与闭包调用形状 {out_shape} 不符"
            )
        return self._I1

    def buoyant_mobility(self, c_bar, m: float, eta1: float, eta2: float, H) -> Array:
        return self._base.buoyant_mobility(c_bar, m, eta1, eta2, H)


def _rheology_from_closure(closure) -> Tuple[float, float, float]:
    """从 HB 闭包对象取线性求解器需要的 ``(η₁, η₂, m)``（供 A-3a 复用线性求解器）。

    ``solve_stream_function`` 的每个闭包调用都要 ``(m, η₁, η₂)``。对 HB 闭包而言：

    - **非牛顿路径**：这三个量只进 R-T1-6 的**地板口径参照**（``1e-6·I₁_牛顿``），
      HB 的真实黏度由 ``κ``/``n`` 携带 ⇒ 取 ``η₁ = κ₁``、``η₂ = κ₂``（牛顿极限下
      ``κ = η`` 精确，非牛顿下只作地板量级参照——地板是算子适定性正则化，
      见 ``hb_closure`` 模块 docstring）。
    - **牛顿极限路径**（``n ≡ 1、τ_y ≡ 0``）：闭包结构性短路到 ``NewtonianClosure``，
      此时 ``(η₁, η₂, m)`` **完整决定闭包值**，故必须取真实黏度口径——而
      ``n=1 ⇔ κ = η`` 恰使其可从闭包的 ``kappa`` 精确还原（``m`` 取闭包的 ``m``
      属性，缺省用 ``κ₁/κ₂``，二者在牛顿极限下相等）。
    """
    kappa = getattr(closure, "kappa", None)
    m_attr = getattr(closure, "m", None)
    if kappa is not None:
        try:
            k0, k1 = float(kappa[0]), float(kappa[1])
            if k0 > 0.0 and k1 > 0.0:
                return k0, k1, (float(m_attr) if m_attr is not None else k0 / k1)
        except (TypeError, IndexError, ValueError):
            pass
    if m_attr is not None and float(m_attr) > 0.0:
        return 1.0, 1.0, float(m_attr)
    return 1.0, 1.0, 1.0


def solve_stream_function_nonlinear(geom: Dict, c_bar, hb_closure, b_field,
                                    Gb=(0.0, 0.0), *,
                                    velocity_scale: float = 1.0,
                                    omega: float = 0.5, tol: float = 1e-6,
                                    max_outer: int = 50) -> Array:
    """解 HB 流体 (4.22) 流函数椭圆方程——**非线性外迭代**（A-3a）。

    算法（brief 的五步，controller 裁定的路径）
    -------------------------------------------
    1. **首轮原样复用线性求解器**（R2）：``solve_stream_function(..., closure=hb_closure)``
       ——闭包取调用方当前注入的 G 状态（``HBClosure`` 未注入 G 时牛顿极限短路仍
       可用；非牛顿闭包未注入会按 ``hb_closure`` 的规定抛 ``RuntimeError``）。
    2. 由 ``∇aΨ`` 经 (2.2) 取局部**间隙平均速度模** ``ū = |(w̄, v̄)|``。
    3. **逐格独立反求** ``G``：解 1 维闭包方程 ``F(g̃) ≡ Ī₁(g̃)·g̃ − ū = 0``
       （Gb=0 的生产口径；``Ī₁`` 走 :func:`gap_solver.solve_g_from_mean_velocity_batch`
       的批量闭包 + 标量求根，**不**逐格跑定均速 Uzawa）。无正根的格（``ū=0`` 或
       全场未屈服）注入 ``G=0`` ⇒ 交由闭包自身的 **R-T1-6 地板**（单次告警），
       整场不失败。
    4. 在反求的 G 处重取闭包 ``I₁``（**同一函数/同一口径**：闭包内部即
       ``closure_integrals_batch``；R-T4-4 要求反求与更新同口径）。
    5. **收敛判定 → 欠松弛 → 重解**：``ΔĪ₁ ≡ 0``（逐位，牛顿极限）⇒ 立即返回首轮 Ψ；
       否则以 ``‖ΔI₁‖_∞/‖I₁‖_∞ < tol`` 判定，欠松弛
       ``I₁ ← (1−ω)·I₁ + ω·I₁_new`` 后用**冻结流动度**重解线性 Poisson。
       收敛（非逐位分支）后再做一次终解，使返回的 Ψ 与闭包当前状态严格对应。
    6. ``max_outer`` 轮内不收敛 ⇒ **``RuntimeError``**（不静默回退牛顿闭包——回退会
       静默改变物理）。⚠️ **"回退牛顿闭包"由调用方负责**——只有调用方知道回退后如何
       继续该时间步：Task 6 在 ``enable_hb_closure=True`` 时对该异常 ``try/except`` →
       告警 + 回退牛顿路径（本函数只负责"抛错不静默"，见 Raises 段）。

    ⚠️ **适用域**：HB 闭包属**扩展应用，未获外部验证**（B&F25 自述尚无外部验证）；
    Z&F22/23 的判据与结论仅严格适用于**竖直井 + 牛顿流体**，不得引其作 HB/斜井依据。

    ⚠️ **G 反求的量/方向口径**：``Gb=0`` 时闭包各向同性 ⇒ ``I₁`` 只依赖 ``|G|``，
    故本函数注入的是**量值场**（第 1 分量为正的标量场，``gap_solver._as_2vec`` 口径）；
    方向/符号不进闭包、也不进算子（算子只消费标量 I₁），故不作符号判决。若下游
    需要带方向/符号的 G，须另立口径（本函数不支持）。

    ⚠️ **标度口径（单一标度，见模块 docstring 的"非线性外迭代"段）**：反求与闭包
    求值都在 ``solve_stream_function`` 的**单位通量 Ψ 口径**内自洽；接线须让喂进反求的
    ``ū`` 与闭包的 ``κ``/``τ_Y`` 口径用**同一**标度 ``ŵ = q_half/π``（T9）——该标度由
    **``velocity_scale`` 形参显式承载**（Task 6 生产接线传 ``ŵ = q_half/π``，推荐路径 =
    物理 (κ, τ_Y) + ū×ŵ；``velocity_scale=1.0`` 默认 = 既有行为）。口径不一致（闭包用
    物理参数而 ``ū`` 未乘 ``ŵ``；模块口径速度比物理**大** ``1/ŵ`` 倍，生产井
    ``ŵ = q_half/π ≪ 1``）会让闭包看到的应力比物理**大** ``ŵ^{-n}`` 倍（呼101 量级
    实测 ≈ 53×，探针 probe_direction_chain.py [3]）⇒ 屈服门槛形同虚设 ⇒
    **屈服效应被削弱（weakened）**；若某工况 ``ŵ > 1`` 则方向反转。

    Args:
        geom: 几何字典（同 :func:`solve_stream_function`）。
        c_bar: (ny,nz) 间隙平均水泥体积分数。
        hb_closure: HB 闭包提供者——须带 ``n``/``kappa``/``tau_y`` 属性
            （如 :class:`hb_closure.HBClosure`）且实现 ``set_pressure_gradient``。
        b_field: (2,ny,nz) 浮力全向量（(4.22) 字面分组，唯一口径）。
        Gb: 反演口径下的浮力向量（(2.12)）；**只支持零向量**（生产口径：浮力经
            ``b_field`` 承担）。非零 ⇒ ``NotImplementedError``（见 Raises）。
        velocity_scale: 反求所喂速度的显式标度（A-3b，R-T5-1 REVISED）：反求循环把
            (2.2) 单位通量口径的 ``ū`` 乘以该值后再解 ``H·ū = I₁·G``。语义
            ``ū_fed = velocity_scale·ū_module``；生产接线（Task 6）传生产换算锚
            ``ŵ = q_half/π``（``annulus_d2dga``：``w = w_unit·(q_half/π)``）+ 物理
            (κ, τ_Y) 参数（推荐路径）；``1.0``（默认）= 既有行为（模块口径自洽，
            牛顿极限等既有测试逐位不变）。
        omega: 欠松弛因子 ∈ (0,1]（``1`` ⇒ 不松弛）。
        tol: 外迭代收敛容差（``Ī₁`` 场的 ∞-范数相对变化；``0`` ⇒ 只认逐位相等）。
        max_outer: 最大外迭代轮数（含首轮）≥1。

    Returns:
        Ψ 场 (ny,nz)（单位通量口径，同 :func:`solve_stream_function`）；收敛后与
        闭包**当前**状态的线性解逐位一致。

    Raises:
        ValueError: 形状/取值非法；闭包缺少 ``n``/``kappa``/``tau_y`` 属性。
        NotImplementedError: ``Gb ≠ 0``（本函数只实现 Gb=0 反演口径，见下）。
        RuntimeError: ``max_outer`` 轮内未收敛（消息含残差轨迹与补救方向）。

    ⚠️ **Gb 口径（显式，不静默忽略）**：反演只实现 ``Gb = 0``——即 **Phase A 的生产
    口径**（浮力由 (4.22) 的 ``b`` 向量承担、闭包 ``Gb`` 恒 0，见 ``hb_closure`` 的
    "Gb 分工"段；本仓 ``b_field`` ↔ ``Gb`` 的换算未定义）。故：形参 ``Gb`` 非零即抛
    ``NotImplementedError``（忽略它会把 ``−Ī₂Gb̃`` 混进 ``Ī₁G̃``、系统性高估局部应力
    ⇒ 低估屈服效应）；同时以 duck-typing 校验**闭包实际注入的 ``Gb``**
    （``HBClosure._Gb``），与形参不一致同样抛错（防口径分裂）。
    """
    H = np.asarray(geom["H"], dtype=float)
    if H.ndim != 2:
        raise ValueError("geom['H'] 必须是 (ny,nz) 二维数组")
    ny, nz = H.shape
    c = np.asarray(c_bar, dtype=float)
    if c.shape != H.shape:
        raise ValueError("c_bar 形状须与 geom['H'] 相同 (ny,nz)")
    omega = float(omega)
    if not (0.0 < omega <= 1.0):
        raise ValueError(f"omega 须在 (0,1] 内（欠松弛因子），得到 omega={omega!r}")
    tol = float(tol)
    if not (tol >= 0.0):
        raise ValueError(f"tol 须为非负（相对收敛容差），得到 tol={tol!r}")
    if not isinstance(max_outer, int) or isinstance(max_outer, bool) or max_outer < 1:
        raise ValueError(f"max_outer 须为 ≥1 的整数，得到 max_outer={max_outer!r}")

    Gb_arr = np.asarray(Gb, dtype=float)
    if not np.all(np.isfinite(Gb_arr)):
        raise ValueError(f"Gb 须为有限值，得到 {Gb!r}")
    # velocity_scale（A-3b，Task 6 接线）：外迭代反求所喂速度的显式标度——
    # 语义 ū_fed = velocity_scale·ū_module（R-T5-1 REVISED：生产接线传 ŵ = q_half/π，
    # 使闭包在物理应力标度上求值；默认 1.0 = 既有行为，×1.0 对 IEEE 浮点逐位无扰）。
    velocity_scale = float(velocity_scale)
    if not (np.isfinite(velocity_scale) and velocity_scale > 0.0):
        raise ValueError(f"velocity_scale 须为正的有限值，得到 {velocity_scale!r}")
    if np.any(Gb_arr != 0.0):
        raise NotImplementedError(
            "solve_stream_function_nonlinear 只实现 **Gb=0** 的反演口径（Phase A 生产口径："
            "浮力由 (4.22) 的 b 向量承担）。Gb≠0 时 (2.13) 的 −Ī₂·Gb̃ 项使应力场非共线、"
            "需 2 维反演——本函数不支持且**不静默忽略**（忽略会系统性高估局部应力）。"
            "量级实测与建议见 task-5-report.md §Important-3。"
        )
    gb_inj = getattr(hb_closure, "_Gb", None)      # duck-typing（非协议成员；HBClosure 有）
    if gb_inj is not None and np.any(np.asarray(gb_inj, dtype=float) != 0.0):
        raise NotImplementedError(
            f"闭包已注入非零 Gb（{gb_inj!r}）：反演口径必须与闭包一致，本函数只实现 Gb=0。"
            "请把浮力交由 (4.22) 的 b_field 承担，或先 set_pressure_gradient(G, Gb=(0,0))。"
        )

    n_hb = getattr(hb_closure, "n", None)
    kappa_hb = getattr(hb_closure, "kappa", None)
    tau_y_hb = getattr(hb_closure, "tau_y", None)
    if n_hb is None or kappa_hb is None or tau_y_hb is None:
        raise ValueError(
            "solve_stream_function_nonlinear 需要 HB 闭包参数（对象属性 n/kappa/tau_y）"
            "做逐格反求 G，当前闭包缺少这些属性（非 HBClosure？）。若闭包的 I₁ 与 G "
            "无关（NewtonianClosure/PowerLawGapClosure 等），线性解即精确解——"
            "请直接调用 solve_stream_function。"
        )
    eta1, eta2, m = _rheology_from_closure(hb_closure)

    def _solve_linear(closure) -> Array:
        return solve_stream_function(geom, c, eta1, eta2, m, b_field,
                                     closure=closure, ny=ny, nz=nz)

    def _mobility_now() -> Array:
        I1 = np.asarray(hb_closure.mobility(c, m, eta1, eta2, H), dtype=float)
        if I1.shape != H.shape:
            raise ValueError(
                f"闭包返回的 I₁ 形状 {I1.shape} 与 c̄/H 形状 {H.shape} 不符"
            )
        return I1

    # ---- 第 1 轮：原样复用线性求解器（R2）--------------------------------
    psi = _solve_linear(hb_closure)
    I1_cur = _mobility_now()                 # 与本轮 Ψ 同口径的 Ī₁ 场
    rel = np.inf
    for _round in range(1, int(max_outer) + 1):
        # ① (2.2) 局部间隙平均速度模——**按各分量的物理换算因子加权**（Important-4）：
        # 生产换算 w_phys = (q/2)·w̄/π、v_phys = (q/2)·v̄（T9 π 推导）⇒ 模块 (2.2) 的
        # w̄ 与 v̄ 相差 π 倍口径，直接 hypot 会把两个口径混在一起；把 v̄ 乘 π 归一到 w̄
        # 口径后取模（整体标度 ŵ 由接线方对齐，见模块 docstring 的"标度口径"段）。
        w, v = velocity_from_stream_function(psi, geom)
        u_mag = np.hypot(w, _VELOCITY_COMPONENT_RATIO * v)
        # velocity_scale（R-T5-1 REVISED）：ū_fed = ŵ·ū_module——把模块单位通量口径的
        # ū 缩放到物理速度口径再反求（生产接线 ŵ = q_half/π）。默认 1.0 ⇒ ×1.0 对
        # IEEE 浮点逐位无扰（含 ±0/NaN），既有调用方行为不变。
        u_mag = u_mag * velocity_scale
        # ② 逐格反求 G（批量闭包 + 标量求根；Gb=0 生产口径）
        inv = solve_g_from_mean_velocity_batch(
            c.reshape(-1), n_hb, kappa_hb, tau_y_hb, u_mag.reshape(-1),
            H=H.reshape(-1), rtol=_INVERSE_RTOL, max_iter=_INVERSE_MAX_ITER,
        )
        # ③ 注入（无正根格 G=0 ⇒ 闭包侧 R-T1-6 地板 + 单次告警）
        hb_closure.set_pressure_gradient(inv.G.reshape(H.shape))
        I1_new = _mobility_now()
        # ④ 收敛判定：ΔĪ₁ ≡ 0（逐位短路，牛顿极限 R2）→ 直接返回首轮 Ψ
        if np.array_equal(I1_new, I1_cur):
            return psi
        den = max(float(np.max(np.abs(I1_new))), np.finfo(float).tiny)
        rel = float(np.max(np.abs(I1_new - I1_cur))) / den
        if rel < tol:
            # 终解：返回的 Ψ 与闭包**当前**状态严格对应（多一次线性求解）
            I1_cur = I1_new
            return _solve_linear(_FrozenMobilityClosure(I1_cur, hb_closure))
        # ⑤ 欠松弛 + 冻结流动度重解
        I1_cur = (1.0 - omega) * I1_cur + omega * I1_new
        psi = _solve_linear(_FrozenMobilityClosure(I1_cur, hb_closure))
    raise RuntimeError(
        f"非线性外迭代在 max_outer={max_outer} 轮内未收敛：Ī₁ 场 ∞-范数相对变化 = "
        f"{rel:.3e} > tol={tol:.1e}（ω={omega}）。补救方向：增大 max_outer 或降低 ω"
        "（ω→1 时非线性外迭代可能振荡）；本函数**不静默回退牛顿闭包**——需要回退"
        "请在调用方显式捕获本异常并自行告警。"
    )
