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
"""

from __future__ import annotations

from typing import Dict, Tuple

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numpy.typing import NDArray

from cemdisp.models2d.hb_closure import ClosureProvider, NewtonianClosure
from cemdisp.models2d.two_layer import mobility_i1

Array = NDArray[np.float64]

__all__ = ["solve_stream_function", "velocity_from_stream_function"]


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
            ``None`` ⇒ 不施加（逐位 = HEAD）。非 ``None`` 时在算子系数上乘
            ``max(1−wall, _WALL_CONDUCTANCE_FLOOR)``：冻结区流动度→0，Ψ 局部
            趋于常数 ⇒ 该处速度→0（static wall layer），且保持椭圆结构。
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
    a_cell = 1.0 / (2.0 * I1)        # φ-槽系数 1/(2I₁)
    c_cell = r_a / (2.0 * I1)        # ξ-槽系数 r_a/(2I₁)
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
