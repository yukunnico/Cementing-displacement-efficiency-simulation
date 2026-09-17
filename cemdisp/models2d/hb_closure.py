"""流变闭包提供者协议（Phase B-1）与 HB 闭包实现（Phase A-2）。

把 ``stream_function`` 求解器与具体闭包解耦：求解器只依赖 :class:`ClosureProvider`
协议。:class:`NewtonianClosure` 逐位等于 ``two_layer.mobility_i1/i2``
（Z&F22 (4.21a)/(4.21b)）。Phase A 的 :class:`HBClosure` 实现同一协议，额外依赖
局部压力梯度 G（由外层迭代注入）。

适用域提醒（务必随论文口径一致）
--------------------------------
源模型 Z&F22/23 的判据与结论**仅严格适用于竖直井 + 牛顿流体**；本协议的任何
幂律/HB 实现均属**扩展应用，未获外部验证**（B&F25 自述无外部验证）。见
``docs/源模型口径与适用域声明.md`` §1 声明 3。

References
----------
Zhang & Frigaard (2022), *JFM* **947**, A32：(4.21a)/(4.21b)。
Bararpour & Frigaard (2025), *JFM* **1022**, A15：(2.13)—(2.15) HB 闭包定义、
(2.29)—(2.34) HB 无量纲群、附录 A.2.1/A.2.2（一维间隙弱解，``gap_solver``）。
"""

from __future__ import annotations

import warnings
from typing import Protocol, runtime_checkable

import numpy as np
from numpy.typing import NDArray

from cemdisp.models2d.gap_solver import closure_integrals_batch
from cemdisp.models2d.two_layer import mobility_i1, mobility_i2

Array = NDArray[np.float64]

# R-T1-6：全场未屈服（static wall layer）点的流动度地板。
# 口径沿用 Phase B 先例 ``stream_function._WALL_CONDUCTANCE_FLOOR``（那里是
# ``I₁_eff = I₁·max(1−wall, 1e-6)``，防止 wall≡1 时算子系数全零 ⇒ 矩阵奇异）。
# 此处取**相对**地板 ``I₁_eff = max(I₁_HB, 1e-6·I₁_牛顿)``——I₁ 量纲为
# H³/√(η₁η₂)，绝对地板会随井/流体量级失配；相对牛顿流动度则与量纲无关，物理上
# 等价于"冻结格的电导最多降到开放值的 1e-6"。⚠️ 这是**算子适定性正则化**，不是
# 物理值：真实物理值是 I₁=0（间隙无流动）。
_STATIC_WALL_MOBILITY_FLOOR = 1.0e-6

_MISSING = object()  # 缓存缺省哨兵（与"缓存了未屈服点"区分）


@runtime_checkable
class ClosureProvider(Protocol):
    """两层闭包协议：返回间隙平均流动度 I₁ 与浮力流动度 I₂。

    约定（Z&F22 (4.21a)/(4.21b) 口径）：``c_bar`` 为 (ny,nz) 场或标量；
    ``m``/``eta1``/``eta2`` 为标量；``H`` 为 (ny,nz) 场或标量；返回与
    ``c_bar``/``H`` 广播后的 float 数组。
    """

    def mobility(self, c_bar, m: float, eta1: float, eta2: float, H) -> Array:
        """间隙平均流动度 I₁(c̄, m; η₁, η₂, H)（式 (4.21a)）。"""
        ...

    def buoyant_mobility(self, c_bar, m: float, eta1: float, eta2: float, H) -> Array:
        """浮力流动度 I₂(c̄, m; η₁, η₂, H)（式 (4.21b)）。"""
        ...


class NewtonianClosure:
    """牛顿两层闭包（Z&F22 (4.21a)/(4.21b)）——逐位等于 ``two_layer`` 原函数。

    本类不做任何数值加工，仅把 ``mobility_i1``/``mobility_i2`` 适配到
    :class:`ClosureProvider` 协议。它的存在保证「closure 缺省注入」与 HEAD
    逐位一致（B-1 的 L1 硬约束）。
    """

    def mobility(self, c_bar, m: float, eta1: float, eta2: float, H) -> Array:
        return np.asarray(mobility_i1(c_bar, m, eta1=eta1, eta2=eta2, H=H), dtype=float)

    def buoyant_mobility(self, c_bar, m: float, eta1: float, eta2: float, H) -> Array:
        return np.asarray(mobility_i2(c_bar, m, eta1=eta1, eta2=eta2, H=H), dtype=float)


class PowerLawGapClosure:
    """牛顿闭包 + 幂律间隙一阶修正（B-3，**近似**，非 B&F25 闭包口径）。

    I₁ → I₁ · (H/H̄)^{1/n − 1}（H̄ = 全场均值）。依据：单流体幂律槽流
    ū ∝ H^{1+1/n}G^{1/n} ⇒ I₁ = Hū/G ∝ H^{2+1/n}G^{1/n−1}；相对牛顿 H³ 的
    场修正即 (H/H̄)^{1/n−1}。G 的纯标量因子在椭圆解 flux 归一化中消去，
    故只有 H 依赖进入算子形状。n=1 时因子恒为 1（逐位退化为牛顿）。

    ⚠️ 适用域：一阶近似，不得引 B&F25 作为方法学依据（见设计规格 §1.5）。
    """

    def __init__(self, n: float, base: "ClosureProvider | None" = None) -> None:
        self._n = float(n)
        if not (self._n > 0.0):
            # n ≤ 0（含 0.0/NaN）会使 1/n 无定义或非物理——构造即拒绝。
            raise ValueError(f"幂律指数 n 必须为正（1/n 需有定义），得到 n={n!r}")
        self._base = base if base is not None else NewtonianClosure()

    def mobility(self, c_bar, m: float, eta1: float, eta2: float, H) -> Array:
        I1 = np.asarray(self._base.mobility(c_bar, m, eta1, eta2, H), dtype=float)
        if abs(self._n - 1.0) < 1e-12:
            return I1
        Harr = np.asarray(H, dtype=float)
        Hbar = float(np.mean(Harr))
        return I1 * (Harr / Hbar) ** (1.0 / self._n - 1.0)

    def buoyant_mobility(self, c_bar, m: float, eta1: float, eta2: float, H) -> Array:
        """浮力流动度 I₂ —— **故意不修正**。

        本类只覆盖 I₁ 的 H 依赖（I₁ ∝ H^{2+1/n}）；I₂ = H⁴ 的幂律标度未推导，
        故直接委托 base（牛顿闭包）。当前 solver 只调 :meth:`mobility`（I₂ 经
        ``χ`` 的 ``i2_field`` 另有牛顿闭式路径消费），该口子暂无影响；若未来
        I₂ 也进算子，须先补幂律标度推导。
        """
        return self._base.buoyant_mobility(c_bar, m, eta1, eta2, H)


class HBClosure:
    """B&F25 HB 两层闭包（附录 A 一维间隙弱解），实现冻结的 :class:`ClosureProvider`。

    与 :class:`NewtonianClosure` 的关系
    ----------------------------------
    ``n=(1,1)`` 且 ``τ_y≡0`` ⇒ :attr:`is_newtonian_limit` 为 True，两个方法
    **结构性短路**、逐位委托 :class:`NewtonianClosure`（不查表、不进
    ``gap_solver``）。理由：Task 5 与全局 L1 硬约束要求"牛顿极限逐位 = HEAD"，
    任何数值路径（Uzawa 解 vs 解析闭式，差 ~1e-9）都做不到逐位。

    返回口径（**量纲**，= ``two_layer.mobility_i1/mobility_i2``）
    -----------------------------------------------------------
    ``mobility``/``buoyant_mobility`` 返回量纲 ``I₁ = H²Ĩ₁``/``I₂ = H³Ĩ₂``
    （即 ``gap_solver`` 的 ``GapSolution.I1/I2``），满足 (2.13) 的唯一自洽对::

        H·ū = I₁·G − ₂·(Gb/H)    （等密度时即 H·ū = I₁·G）

    ⚠️ **注入的 ``G`` 必须与 ``solve_fixed_mean_velocity`` 返回的
    ``GapSolution.G`` 同口径**（量纲 Pa/m，已含 ``/H`` 换算）——Task 5 外迭代由
    ``∇aΨ`` 经定均速模式反求 ``G``，本类消费的就是它。数值核验（``H=0.008``、
    ``ny_gap=201``）：``ū/(I₁G/H)`` 偏离 1 ≤3e-5，全部来自 ``u_bar`` 的
    trapezoid 离散误差（O(h²)），非口径差异。

    计算路径：为什么是"真实 G 求值"而非"查表 (c̄, B)"
    ------------------------------------------------
    每个格点的闭包值一律由 :func:`gap_solver.closure_integrals_batch` 在**注入的 G/Gb** 与
    真实应力场 ``τ̃(ỹ)=c·ỹ+d``（含屈服面切段 ⇒ 真实塞流 / static wall layer）上
    现算。``(c̄, B)`` 表**只在固定 G 下自洽**：``B = max{τ̂_Y}/(μ̂eγ̇₀)`` 里的
    ``γ₀`` 与应力尺度绑定，而应力尺度随局部 ``G`` 变（``B`` 与 ``G`` 一一对应，
    故用固定 ``B`` 的二维表会系统性偏离——实测见 ``task-4-report.md`` §「表 vs
    真实计算」）。构造期传入的 ``m``/``B`` 只是 (2.32)/(2.33) 的**无量纲参数化**
    （论文报告 / 表索引），**不参与**运行时求值。

    ⚠️ **本类不求解 ``u``**（Task 4.5）：只需 ``I₁/I₂``，而闭包量由**解析应力场**
    求积得到（见 ``gap_solver`` 模块 docstring：与 Uzawa 解无关）⇒ 4.5 起改走
    ``closure_integrals_batch``，**不再**跑 Uzawa。语义影响只有一处：Task 4 的逐点版
    会把 ``gap_solver`` 的 ``RuntimeError``（Uzawa 不收敛/发散）传播出来，现已不会
    发生（不再迭代）；``I₁/I₂`` 的取值、缓存计数、告警与地板行为**逐位/逐语义不变**
    （同场实测 10000/10000 点逐位相同，见 task-4-report.md §10）。

    ️ 为何用**定 G** 而非定均速模式求值（controller 接口澄清的偏差，已上报）：
    注入量就是 ``G``，定 G 模式无需任何**任意的**参考均速 ``ū*`` 来定尺度；而
    两模式的不动点逐式相同（见 ``solve_fixed_mean_velocity`` docstring）⇒ 在注入
    的 ``G`` 处二者给出同一 ``I₁``（交叉核对见
    ``tests/contract/test_hb_closure_tabulated.py::test_hb_closure_matches_gap_solver_fixed_G_and_mean_velocity``，
    实测 rel ≤1e-8）。实测代价差 ~100×（同参数同 ``ny=201``：定 G **7 步/17 ms**
    @ ``tol=1e-10``；定均速 **tol 放宽到 1e-6 仍需 1506 步/2.0 s**，``tol=1e-8``
    在 ``max_iter=2000`` 内不收敛）——Task 5 若用定均速模式反求 G，须预期该量级。

    ``Gb`` 的分工（R-T1-3 余波）
    ---------------------------
    ``set_pressure_gradient(G, Gb=...)`` 的 ``Gb`` 是 (2.18) 的浮力**向量**，
    量级由 ``solve_stream_function`` 的 ``b_field`` 装配、**不由本类决定**。默认
    ``Gb=(0,0)``：Phase A 的生产接线把浮力留给 (4.22) 的 ``b`` 向量承担，本类只
    消费 ``G``（分工须由 Task 5 复核）；调用方若显式注入非零 ``Gb``，则 (2.13) 的
    ``−Ī₂Gb/H`` 项与本类求值同时计入——**不得**两处同时算同一份浮力。非共线
    ``G``/``Gb`` 由向量版求解器支持（R-T1-3），本类不做共线假设，``Gb`` 也不是
    报告量（注入即进真实应力场）。

    R-T1-6 地板（全场未屈服 / ``G=Gb=0`` 的退化格点）
    ------------------------------------------------
    这类格点闭包**合法地**无定义（``gap_solver`` 抛 ``ValueError``：无流动 ⇒ ``I₁=0``、
    ``q₀`` 无定义）。本类**捕获该异常**（不让它逃逸到求解器），返回::

        I₁_eff = max(I₁_HB, 1e-6·I₁_牛顿),  I₂ = 0

    ``I₂=0`` 是精确值（``1/η̃≡0`` 处处 ⇒ ``Ĩ₂`` 的积分恒为 0）；``I₁`` 的地板是
    **算子适定性正则化**（复用 Phase B 先例 ``stream_function._WALL_CONDUCTANCE_FLOOR``
    的同式思路：相对地板，与量纲无关），不是物理值。地板化点数记入
    :attr:`n_static_wall_points`（**最近一次调用**的口径），并**每实例告警一次**
    （``RuntimeWarning``，不每格刷屏）。⚠️ 输入校验（``c̄``/``H``/``G`` 形状与有限性）在本类内先做 ⇒ 只有
    "退化格点"会走到该 ``except``，真输入错误仍是显式 ``ValueError``。
    ⚠️ ``gap_solver`` 的 ``RuntimeError``（Uzawa 不收敛/发散）**不捕获**——按全局
    约定"不收敛不得静默回退"。

    ``q₀`` 不在本类（协议边界）
    --------------------------
    冻结协议只有 ``mobility``/``buoyant_mobility``（I₁/I₂）；(4.25) 的各向同性通量
    仍由 ``two_layer.isotropic_flux_q0`` 的**牛顿闭式**承担，本类不算 HB 口径的
    ``q``（``gap_solver`` 其实已给出该量，但协议无对应出口）。若 Task 6 需要 HB 的
    ``q``，须先扩协议（Phase A 未授权改协议签名）。

    缓存
    ----
    ``solve_fixed_G`` 结果按**精确键** ``(c̄, Gx, Gy, Gbx, Gby, H)`` 缓存（一个
    解同时给出 ``I₁``/``I₂`` ⇒ :meth:`buoyant_mobility` 紧跟 :meth:`mobility`
    时全命中）。键含 ``G`` 是必需的（同 ``c̄`` 不同 ``G`` 必须重解）；``G`` 场均匀
    时该缓存自动退化为"按唯一 ``c̄`` 值缓存"。缓存生命周期 = 一次
    ``set_pressure_gradient`` 注入（那里清空），因此规模被钉在单个场（≈``ny·nz``
    条），长跑不会无界增长。计数见 :attr:`cache_stats`（``misses`` = 未被缓存覆盖的
    **格点数**；这些键先去重再整批送 ``closure_integrals_batch``，故实际求解的唯一键数
    ≤ ``misses``；同一次调用内重复键算命中，与逐点版时序等价）。
    ⚠️ 缓存只存**原始闭包值**：地板在每次调用时按**当次**传入的 ``m``/``η₁``/``η₂``
    现算（改 ``m`` 不会命中旧地板）。
    ⚠️ 代价量级（实测，呼101 真实场 ``cement_final`` ``(40,250)``=10000 点、
    ``ny_gap=201``、唯一 ``c̄``=6647）：Task 4.5 起闭包值走
    :func:`gap_solver.closure_integrals_batch`（**批量**；闭包量由解析应力场求积、
    与 Uzawa 解无关 ⇒ 本类不需要 ``u``）⇒ 一次 ``mobility()`` **0.339 s**、
    ``buoyant_mobility()`` 0.013 s（Task 4 逐点版同场 **72.9 s** ⇒ **215×**），
    且输出与逐点版**逐位相同**（10000/10000 点，见 task-4-report.md §10）。
    （Task 4 的历史实测：均匀 ``G`` 82 s / 命中 33.5%；非均匀 ``G`` 116 s / 26.1%。）

    ⚠️ 适用域
    ----------
    HB 闭包属**扩展应用，未获外部验证**（B&F25 自述尚无外部验证；Z&F22/23 的判据
    与结论仅严格适用于**竖直井 + 牛顿流体**，不得引其作 HB/斜井依据）。默认开关
    ``enable_hb_closure=False``（Task 6 接线前不影响生产路径）。

    References
    ----------
    Bararpour & Frigaard (2025), *JFM* **1022**, A15：(2.13)—(2.15) 闭包定义、
    (2.29)—(2.34) 无量纲群（见 ``two_layer.hb_groups``）、附录 A.2.1（定 G 弱解）。
    Zhang & Frigaard (2022), *JFM* **947**, A32：(4.21a)/(4.21b) 牛顿极限
    （``NewtonianClosure`` 逐位锚）。
    """

    def __init__(self, n, kappa, tau_y, m: float, B: float, *, ny_gap: int = 201) -> None:
        """构造 HB 闭包提供者。

        Args:
            n: 幂律指数 ``(n₁, n₂)``——流体 1（被顶替液，壁面带）/ 流体 2（顶替液，
                水泥，中线带），与 ``gap_solver``/``two_layer`` 同下标约定。
            kappa: 稠度系数 ``(κ̂₁, κ̂₂)``（量纲，Pa·sⁿ）。
            tau_y: 屈服应力 ``(τ̂_{Y,1}, τ̂_{Y,2})``（Pa，≥ 0）。
            m: 黏度比——(2.32) 口径（``γ̇₀`` 处表观黏度比，见 ``hb_groups``）。
                **仅供报告/索引**：运行时求值用 ``kappa``，运行时 ``mobility(m=...)``
                的 ``m`` 只作告警地板的口径参照。
            B: Bingham 数（(2.33)，由 ``hb_groups`` 在某个参考 ``γ̇₀`` 下算出）。
                **仅供报告/索引**，不参与运行时求值（局部应力尺度随 ``G`` 变）。
            ny_gap: 一维间隙求解器的格边点数（``gap_solver`` 的 ``ny``；越大越准、
                越慢，实测 ``ny=41/201`` 单解 ~6.7/9.8 ms）。

        Raises:
            ValueError: ``n ≤ 0``、``κ̂ ≤ 0``、``τ̂_Y < 0``、``m ≤ 0``、``B < 0``
                或 ``ny_gap`` 非 ≥3 的整数。
        """
        self.n = _as_pair(n, "n", positive=True)
        self.kappa = _as_pair(kappa, "kappa", positive=True)
        self.tau_y = _as_pair(tau_y, "tau_y", positive=False)
        self.m = float(m)
        self.B = float(B)
        if not (np.isfinite(self.m) and self.m > 0.0):
            raise ValueError(f"m（黏度比）须为正的有限值，得到 m={m!r}")
        if not (np.isfinite(self.B) and self.B >= 0.0):
            raise ValueError(f"B（Bingham 数）须为非负的有限值，得到 B={B!r}")
        if not isinstance(ny_gap, int) or isinstance(ny_gap, bool) or ny_gap < 3:
            raise ValueError(f"ny_gap 须为 ≥3 的整数（格心数 = ny−1），得到 {ny_gap!r}")
        self.ny_gap = ny_gap

        # R2：牛顿极限结构性短路（不查表、不进 gap_solver ⇒ 逐位 = NewtonianClosure）
        self._newtonian: NewtonianClosure | None = (
            NewtonianClosure()
            if (self.n == (1.0, 1.0) and self.tau_y == (0.0, 0.0))
            else None
        )
        self._G: Array | None = None
        self._Gb: Array = np.zeros(2)
        self._cache: dict = {}
        self._hits = 0
        self._misses = 0
        self._static_wall_warned = False
        self.n_static_wall_points = 0

    # ------------------------------------------------------------------ #
    # 接口
    # ------------------------------------------------------------------ #
    @property
    def is_newtonian_limit(self) -> bool:
        """``n=(1,1)`` 且 ``τ_y≡0`` ⇒ True（此时两方法逐位委托 ``NewtonianClosure``）。"""
        return self._newtonian is not None

    @property
    def cache_stats(self) -> dict:
        """缓存查找计数（累计）：``{"hits": int, "misses": int}``（misses = 实际解数）。"""
        return {"hits": self._hits, "misses": self._misses}

    @property
    def current_G(self):
        """最近一次注入/反求收敛的 G（只读**副本**；``None`` = 尚未注入）。

        Task 6 修复轮 1（R-T6-2 阶梯 (i)）：供外层接线做时间步间 warm-start
        （``solve_stream_function_nonlinear(initial_G=...)``）。返回 ``.copy()``
        （评审 Minor⑤）：调用方持有快照，后续 ``set_pressure_gradient`` 原地
        清空/改写 ``_G`` 不会穿透到已缓存的上一步 G 场。
        """
        return None if self._G is None else self._G.copy()

    def set_pressure_gradient(self, G, Gb=(0.0, 0.0)) -> None:
        """注入当前局部压力梯度场 ``G``（与浮力向量 ``Gb``），并清空解缓存。

        Args:
            G: 修改压力梯度（(2.11)/(2.13) 口径，量纲，**须与
                ``solve_fixed_mean_velocity`` 返回的 ``GapSolution.G`` 同口径**）。
                接受标量（按第 1 分量 = φ 解释，``gap_solver._as_2vec`` 口径）、
                ``(ny,nz)`` 标量场或 ``(2,ny,nz)`` 逐点向量场。
            Gb: 浮力向量（(2.12)）；形状规则同 ``G``。默认 ``(0,0)``（Phase A 生产
                接线把浮力交给 (4.22) 的 ``b`` 向量，见类 docstring 的"Gb 分工"）。

        ⚠️ 形状/有限性在 :meth:`mobility` 调用时按 ``c̄``/``H`` 的广播形状校验。
        """
        self._G = np.asarray(G, dtype=float)
        self._Gb = np.asarray(Gb, dtype=float)
        # 缓存生命周期 = 一次注入：既保证"同 G 内复用"（mobility→buoyant_mobility、
        # 同 G 重复调用），又把缓存规模钉在单场（≈ny·nz 条），防长跑内存无界增长。
        self._cache.clear()

    def mobility(self, c_bar, m: float, eta1: float, eta2: float, H) -> Array:
        """间隙平均流动度 ``I₁``（量纲；(2.14) 的一维间隙弱解求积）。

        与 :class:`NewtonianClosure` 同名同义；``η₁/η₂`` 只用于地板口径参照
        （HB 的真实黏度由 ``kappa``/``n`` 携带），``m`` 同上（见类 docstring）。
        """
        if self._newtonian is not None:
            return self._newtonian.mobility(c_bar, m, eta1, eta2, H)
        return self._evaluate(c_bar, m, eta1, eta2, H)[0]

    def buoyant_mobility(self, c_bar, m: float, eta1: float, eta2: float, H) -> Array:
        """浮力流动度 ``I₂``（量纲；(2.15) 的一维间隙弱解求积）。

        与 :meth:`mobility` 共用同一批闭包求积（缓存）⇒ 紧随其后调用时全命中。
        地板格点返回精确值 ``0``（``1/η̃≡0`` ⇒ ``Ĩ₂`` 恒 0，见类 docstring）。
        """
        if self._newtonian is not None:
            return self._newtonian.buoyant_mobility(c_bar, m, eta1, eta2, H)
        return self._evaluate(c_bar, m, eta1, eta2, H)[1]

    # ------------------------------------------------------------------ #
    # 内部：真实 G 求值
    # ------------------------------------------------------------------ #
    def _evaluate(self, c_bar, m: float, eta1, eta2, H):
        """逐格点求 ``(I₁, I₂)``（真实 G + 真实应力场；带精确键缓存与地板）。"""
        if self._G is None:
            raise RuntimeError(
                "HBClosure 需要局部压力梯度 G 才能求值：请先调用"
                " set_pressure_gradient(G, Gb=...)（Task 5 外迭代由 ∇aΨ 经定均速模式"
                "反求 G 后注入）。不静默回退牛顿闭包——回退会静默改变结果。"
            )
        out_shape = np.broadcast_shapes(np.shape(c_bar), np.shape(H))
        c = np.broadcast_to(np.asarray(c_bar, dtype=float), out_shape)
        Hf = np.broadcast_to(np.asarray(H, dtype=float), out_shape)
        if np.any(~np.isfinite(c)) or np.any(c < 0.0) or np.any(c > 1.0):
            raise ValueError(
                f"c_bar 须为 [0,1] 内的有限值（gap_solver 定义域），得到 {c_bar!r}"
            )
        if np.any(~np.isfinite(Hf)) or np.any(Hf <= 0.0):
            raise ValueError(f"H 须为正的有限值（物理半隙），得到 {H!r}")
        gx, gy = _components(self._G, out_shape, "G")
        gbx, gby = _components(self._Gb, out_shape, "Gb")
        if not all(np.all(np.isfinite(v)) for v in (gx, gy, gbx, gby)):
            raise ValueError(f"注入的 G/Gb 须为有限值，得到 G={self._G!r}、Gb={self._Gb!r}")

        n_pts = c.size
        cf, Hflat = c.reshape(-1), Hf.reshape(-1)
        gxf, gyf = gx.reshape(-1), gy.reshape(-1)
        gbxf, gbyf = gbx.reshape(-1), gby.reshape(-1)
        i1 = np.empty(n_pts, dtype=float)
        i2 = np.empty(n_pts, dtype=float)
        floored = np.zeros(n_pts, dtype=bool)

        # Task 4.5：先按精确键查缓存，未命中的键**去重**后整批送
        # ``closure_integrals_batch``（闭包量由解析应力场求积，与 Uzawa 解无关 ⇒
        # HBClosure 不需要解流速剖面，省掉迭代）。计数口径与逐点版逐字相同：
        # 同一次调用内重复出现的键算**命中**（等价于逐点缓存时的时序）。
        results: list = [_MISSING] * n_pts
        pending: dict = {}
        for i in range(n_pts):
            key = (cf[i], gxf[i], gyf[i], gbxf[i], gbyf[i], Hflat[i])
            got = self._cache.get(key, _MISSING)
            if got is not _MISSING:
                self._hits += 1
                results[i] = got
            elif key in pending:
                self._hits += 1
                pending[key].append(i)
            else:
                self._misses += 1
                pending[key] = [i]
        if pending:
            uniq = np.array(list(pending), dtype=float)          # (n_unique, 6)
            cb = closure_integrals_batch(
                uniq[:, 0], self.n, self.kappa, self.tau_y,
                uniq[:, 1:3], Gb=uniq[:, 3:5], H=uniq[:, 5],
            )
            for j in range(uniq.shape[0]):
                # 无定义（全场未屈服 / G=Gb=0）⇒ 缓存 None，由下方地板处置
                got = None if cb.undefined[j] else (float(cb.I1[j]), float(cb.I2[j]))
                self._cache[tuple(uniq[j])] = got
                for i in pending[tuple(uniq[j])]:
                    results[i] = got
        for i in range(n_pts):
            got = results[i]
            if got is None:
                floored[i] = True
            else:
                i1[i], i2[i] = got

        n_static = int(np.count_nonzero(floored))
        self.n_static_wall_points = n_static
        # 成因拆分（诊断用）：复算标量路径的两条 ValueError 判据
        zero_drive = (np.hypot(gxf, gyf) == 0.0) & (np.hypot(gbxf, gbyf) == 0.0) & floored
        n_zero = int(np.count_nonzero(zero_drive))
        n_unyield = n_static - n_zero
        if n_static:
            i1_ref = np.broadcast_to(
                np.asarray(mobility_i1(c, float(m), eta1=eta1, eta2=eta2, H=Hf), dtype=float),
                out_shape,
            ).reshape(-1)
            i1 = np.where(floored, _STATIC_WALL_MOBILITY_FLOOR * i1_ref, i1)
            i2 = np.where(floored, 0.0, i2)
            if not self._static_wall_warned:
                self._static_wall_warned = True
                warnings.warn(
                    f"HBClosure：{n_static} 个格点闭包无定义（全场未屈服 ⇒ static wall "
                    f"layer，或 G=Gb=0），已取 I₁ 地板 = {_STATIC_WALL_MOBILITY_FLOOR:g}"
                    f"×I₁_牛顿、I₂=0。成因拆分：全场未屈服 {n_unyield} 个、"
                    f"G 与 Gb 同时为零（零驱动）{n_zero} 个。地板是算子适定性正则化"
                    "（先例 stream_function._WALL_CONDUCTANCE_FLOOR），不是物理值"
                    "（真实物理值为 I₁=0、I₂=0）。每实例只告警一次；"
                    "计数见 n_static_wall_points。",
                    RuntimeWarning,
                    stacklevel=2,
                )
        return i1.reshape(out_shape), i2.reshape(out_shape)


def _as_pair(values, name: str, *, positive: bool):
    """``(值, 值)`` 校验并转为 2 元组 float（``positive=True`` 时要求 > 0，否则 ≥ 0）。"""
    arr = np.asarray(values, dtype=float).reshape(-1)
    if arr.size != 2:
        raise ValueError(f"{name} 须为 2 元组（逐相一个），得到 {values!r}")
    lo_ok = (arr > 0.0) if positive else (arr >= 0.0)
    if not (np.all(np.isfinite(arr)) and np.all(lo_ok)):
        bound = "正" if positive else "非负"
        raise ValueError(f"{name} 须为{bound}的有限值，得到 {values!r}")
    return (float(arr[0]), float(arr[1]))


def _components(arr: Array, out_shape, name: str):
    """把注入的 ``G``/``Gb`` 解析为逐点分量数组（返回两个 ``out_shape`` 形状的数组）。

    接受（与 ``gap_solver._as_2vec`` 的二维向量口径一致）：

    - ``(2,) + out_shape`` ⇒ 逐点向量场，按分量拆开；
    - ``out_shape`` 或标量 ⇒ 沿第 1 分量（φ 方向）的标量场，第 2 分量取 0；
    - ``(2,)`` 且 ``out_shape == ()`` ⇒ 常向量。

    ⚠️ 歧义提示：``out_shape == (2,)`` 时 ``(2,)`` 形状的输入按**标量场**解释
    （此为项目既有标量口径；逐点向量场请传 ``(2,)+out_shape``）。
    """
    shape = tuple(out_shape)
    if arr.shape == (2,) + shape:
        return arr[0], arr[1]
    if arr.shape == shape or arr.ndim == 0:
        return np.broadcast_to(arr, shape), np.zeros(shape)
    if arr.shape == (2,):
        return np.full(shape, float(arr[0])), np.full(shape, float(arr[1]))
    raise ValueError(
        f"{name} 形状非法：须为标量、``c̄``/``H`` 广播形状 {shape} 的场，或 (2,)+{shape} "
        f"的逐点向量场，得到 shape={arr.shape}"
    )
