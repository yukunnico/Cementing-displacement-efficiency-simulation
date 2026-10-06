"""
温变流变公式全集（T0-2：温压耦合的公式层 + 流体分派器）

把用户裁定的温变流变公式实现为 ``fluid_at(fluid, T_c, P_mpa) -> FluidSpec``
派生器，供求解器 T-on 时整体替换流体参数（**绝对替换**口径：有公式的相，
τy→``yield_stress_pa``、μp→宾汉塑性粘度，模型族切 ``RheologyModel.BINGHAM``，
幂律参数清空；密度/名字/角色不变）。无公式/不替换的相原样返回同一对象。

公式来源（逐字录入，系数与计划 §3 公式全集一致）：
温压耦合改进计划_2026-09-30.md §3 + task-3-brief.md（含 2026-09-30 四项裁定）：

- 组 A 水泥（2.1 g/cm³，20–170°C）/ 组 B 水泥（1.9 g/cm³，20–200°C）
- 中间密度水泥：**就近取**组 B/组 A **整式**（Q16，2026-10-06 用户裁定）——
  捕获区 [1.88,2.12]，以 2.000 为界（平局取高密度端=组 A）；**不再做密度插值**
- 隔离液 1.95/2.05 含压二次曲面；**就近取**端点**整式**（Q16）——捕获区 [1.90,2.10]
  （= [1.95,2.05] + 半档 0.05 宽容带），以 2.000 为界（平局取 2.05 式）；
  **不再做六系数密度插值**；捕获区外维持"就近借用端点整式 + borrow 审计"（系数不外推）
- 钻井液（呼101 拟合，域 [40,80]°C clamp；``mud_extrapolate`` 可外推对比）

单位：T=°C，P=MPa，τ/τ₀/τy=Pa，μp=**Pa·s**（锚点表 mPa·s 由测试侧 ÷1000）。
断点（100°C）忠实两段式：T≤100 低温式、T>100 高温式（``smooth_break`` 占位，
默认 False=忠实；True 尚未实现抛 NotImplementedError）。
术语：τ₀≡τy≡YP≡``yield_stress_pa``。

越界一律 clamp 到域端值+审计；``get_audit()`` 查 clamp/借用/不替换等事件，
``reset_audit()`` 清空。审计列表定长 10000（超限丢最旧），接线后按步清空。
"""

from __future__ import annotations

import math
from collections import deque
from dataclasses import dataclass, replace
from typing import Optional

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel

__all__ = ["RheologyFormulaParams", "fluid_at", "get_audit", "reset_audit"]


# ---------------------------------------------------------------------------
# R1（2026-10-06）：公式系数参数化
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class RheologyFormulaParams:
    """温变流变公式的**系数与域界**（R1 参数化；`fluid_at(..., params=...)`）。

    默认值 = 2026-10-01 逐字录入的公式系数，**逐位不变**：`params=None` ⇒ 用模块级
    `_DEFAULT_PARAMS` 单例 ⇒ 与 HEAD 数值路径逐位一致（关2 红线）。

    本类只承载**数值**，不改任何运算次序。改系数一律经本对象显式传入
    （`dataclasses.replace` 派生或直接构造），默认路径永不吃扰动。

    ⚠️ 表达式次序约束（逐位红线的实现前提，改动前必读）
    ---------------------------------------------------
    各 `_*_values` 求值函数的表达式次序**不得重写**：多项式必须保持
    ``q0*T*T + q1*T + q2``（原式 ``a*T*T - b*T + c`` 的 ``-b`` 存为 ``q1``；
    IEEE754 下 ``x - y ≡ x + (-y)``，故逐位相同）。**禁止 Horner 化**——
    实测（210 温点网格）：加法形式 0 处变位，Horner 形式 109/210 处变位。
    """

    # --- 断点（忠实两段式：T≤break_T 低温式、T>break_T 高温式；100°C 不平滑）---
    break_T: float = 100.0
    # --- 组 A 水泥（2.1 g/cm³，20–170°C，纯温度）---
    ca_tau0_q: tuple = (0.002793, -0.391702, 21.6378)   # T≤100 二次
    ca_tau0_l: tuple = (0.624883, -51.7949)             # T>100 一次
    ca_mup: tuple = (0.53129, -0.012495, -0.00204)      # amp·e^(rate·T) + off
    ca_T_lo: float = 20.0
    ca_T_hi: float = 170.0
    # --- 组 B 水泥（1.9 g/cm³，20–200°C，纯温度）---
    cb_tau0_q: tuple = (4.433e-4, -0.07843, 5.962)
    cb_tau0_l: tuple = (0.03206, 1.979)
    cb_mup: tuple = (1.355e-5, -6.570e-4, 0.08215)
    cb_T_lo: float = 20.0
    cb_T_hi: float = 200.0
    # --- 隔离液二次曲面（含压力）：系数序 (常数, T, P, T², T·P, P²) ---
    sp_1p95_ty: tuple = (7.379374, -0.001857, 0.0225351,
                         -0.000157779, -0.00000615055, -0.0000993563)
    sp_1p95_mp: tuple = (0.0912039, -0.000719380, 0.000167466,
                         0.00000277807, -0.00000322380, 0.00000119469)
    sp_2p05_ty: tuple = (10.787391, -0.0206255, 0.0320167,
                         -0.000247286, 0.000219522, -0.000253669)
    sp_2p05_mp: tuple = (0.125694, -0.00108616, 0.000323974,
                         0.00000341573, -0.00000316667, 0.000000714405)
    sp_T_lo: float = 20.0
    sp_T_hi: float = 200.0
    sp_P_lo: float = 0.1
    sp_P_hi: float = 200.0
    sp_p_default: float = 0.1          # P 未传入时的常压缺省（记 p_default 审计）
    # Q16（2026-10-06 用户裁定）「密度就近取」捕获区与平局界：
    #   捕获区 = 锚点区间 + 半档宽容带（带内仍算"近"，整式就近取，不插值）；
    #   捕获区以 _mid 为界左闭右开 ⇒ 平局（恰在中点）取**高密度端**。
    #   捕获区外 = 「远」：维持现状（隔离液借用端点整式+borrow 审计；水泥不替换）。
    sp_lo: float = 1.90
    sp_hi: float = 2.10
    sp_mid: float = 2.000
    sp_anchor_lo: float = 1.95
    sp_anchor_hi: float = 2.05
    cm_lo: float = 1.88
    cm_hi: float = 2.12
    cm_mid: float = 2.000
    cm_anchor_a: float = 2.10           # 组 A 锚点（平局取高端）
    cm_anchor_b: float = 1.90           # 组 B 锚点
    # --- 钻井液（呼101 实测拟合，纯温度，域 [40,80]°C clamp）---
    mud_mup: tuple = (0.297154, -0.02310525)
    mud_tauy: tuple = (21.346710, -0.00830631)
    mud_T_lo: float = 40.0
    mud_T_hi: float = 80.0


_DEFAULT_PARAMS = RheologyFormulaParams()
"""默认参单例（`params=None` ⇒ 本对象；与 HEAD 逐位一致）。"""

# ---------------------------------------------------------------------------
# 审计（模块级；fluid_at 为自由函数，同型于 temperature_field 的实例审计）
# ---------------------------------------------------------------------------
_AUDIT_MAX = 10000  # 上限：deque 定长，超限丢最旧（review minor④；接线后按步 reset）
_AUDIT: deque = deque(maxlen=_AUDIT_MAX)


def get_audit() -> list[dict]:
    """返回审计事件列表副本（kind ∈ clamp / borrow / nearest / no_replace /
    model_assumption / extrapolate / p_default）。

    事件数超上限（10000）时丢弃最旧、保留最新（deque 定长，内存有界）。
    """
    return list(_AUDIT)


def reset_audit() -> None:
    """清空审计事件。"""
    _AUDIT.clear()


def _record(kind: str, fluid: str, detail: str, **extra) -> None:
    event = {"kind": kind, "fluid": fluid, "detail": detail}
    event.update(extra)
    _AUDIT.append(event)


def _clamp_record(value: float, lo: float, hi: float, name: str, param: str) -> float:
    """clamp 到 [lo, hi]；越界时记一条 clamp 审计。"""
    clamped = min(max(value, lo), hi)
    if clamped != value:
        _record(
            "clamp", name,
            f"{param} 越界 {value} → clamp {clamped}（域 [{lo}, {hi}]）",
            param=param, requested=value, clamped=clamped,
        )
    return clamped


# ---------------------------------------------------------------------------
# 组 A 水泥（2.1 g/cm³，20–170°C，纯温度）
# ---------------------------------------------------------------------------
def _cement_a_tau0(T: float, p: RheologyFormulaParams) -> float:
    """τ₀(T)：20≤T≤100 二次式；100<T≤170 一次式（T 已 clamp 进域）。"""
    if T <= p.break_T:
        q0, q1, q2 = p.ca_tau0_q
        return q0 * T * T + q1 * T + q2
    l1, l0 = p.ca_tau0_l
    return l1 * T + l0


def _cement_a_mup(T: float, p: RheologyFormulaParams) -> float:
    """μp(T) = 0.53129·e^(−0.012495·T) − 0.00204 [Pa·s]。"""
    amp, rate, off = p.ca_mup
    return amp * math.exp(rate * T) + off


# ---------------------------------------------------------------------------
# 组 B 水泥（1.9 g/cm³，20–200°C，纯温度）
# ---------------------------------------------------------------------------
def _cement_b_tau0(T: float, p: RheologyFormulaParams) -> float:
    """τ₀(T)：20≤T≤100 二次式；100<T≤200 一次式（T 已 clamp 进域）。"""
    if T <= p.break_T:
        q0, q1, q2 = p.cb_tau0_q
        return q0 * T * T + q1 * T + q2
    l1, l0 = p.cb_tau0_l
    return l1 * T + l0


def _cement_b_mup(T: float, p: RheologyFormulaParams) -> float:
    """μp(T) = 1.355×10⁻⁵·T² − 6.570×10⁻⁴·T + 0.08215 [Pa·s]。"""
    q0, q1, q2 = p.cb_mup
    return q0 * T * T + q1 * T + q2


# ---------------------------------------------------------------------------
# 隔离液二次曲面（含压力）：系数序 (常数, T, P, T², T·P, P²)
# ---------------------------------------------------------------------------
# R1（2026-10-06）：四组六系数 + Q16 密度捕获区/平局界 + 两类温区，全部收进
# `RheologyFormulaParams`（默认值逐字不变）；此处不再保留模块级副本，
# 以免"改参只改一处、另一处漂移"（同 `_SWITCH_DEFAULTS` 单一真源教训）。
# 原 Q16 注释随字段迁入 dataclass docstring 与字段注释。


def _quad6(c: tuple, T: float, P: float) -> float:
    """六系数二次曲面：c0 + c1·T + c2·P + c3·T² + c4·T·P + c5·P²。"""
    return c[0] + c[1] * T + c[2] * P + c[3] * T * T + c[4] * T * P + c[5] * P * P


# ---------------------------------------------------------------------------
# 钻井液（呼101 实测拟合，纯温度）
# ---------------------------------------------------------------------------
def _mud_mup(T: float, p: RheologyFormulaParams) -> float:
    """μp(T) = 0.297154·e^(−0.02310525·T) [Pa·s]。"""
    amp, rate = p.mud_mup
    return amp * math.exp(rate * T)


def _mud_tauy(T: float, p: RheologyFormulaParams) -> float:
    """τy(T) = 21.346710·e^(−0.00830631·T) [Pa]。"""
    amp, rate = p.mud_tauy
    return amp * math.exp(rate * T)


# ---------------------------------------------------------------------------
# 分派表（裁定定稿，含 2026-09-30 补裁）
# ---------------------------------------------------------------------------
# 行2「先导浆、平衡液」按**名字子串**匹配（不全等）：覆盖真实井况别名
# hu1「平衡液(先导泥浆)」（role=WASH，2026-08-29 由旧名"冲洗液"更名）。
# 子串刻意取窄（「先导浆」「平衡液」），不误伤替浆链/隔离液等其他相。
# ⚠️ 子串判定只对**非 SPACER** role 生效——role==SPACER（含 2D 合成相
# "先导浆+隔离液1+隔离液2"）优先走隔离液族（C1，见 _route docstring）。
_ASSUMPTION_NAME_PARTS = ("先导浆", "平衡液")
_CHAIN_NAMES = frozenset({"压塞液"})                  # 替浆链：其余由 role=DISPLACEMENT 覆盖

_NO_REPLACE_DETAIL = {
    "flusher": "冲洗液不替换（返回原 FluidSpec，审计标记）",
    "displacement_chain": "替浆链常数现状（不替换）",
    "unmatched": "不匹配任何档（不替换）",
}


def _route(fluid: FluidSpec) -> tuple[str, Optional[str]]:
    """按分派表返回 (族, 附注)。族 ∈ no_replace / mud / spacer / cement。

    路由优先级裁定（2026-10-01 终审 C1 修复）：
    - **行4/5 `role==SPACER` 优先于行2 名字子串**——2D 合成等效隔离液
      （`annulus_d2dga._composite_spacer_fluid`，名="先导浆+隔离液1+隔离液2"、
      role=SPACER、ρ≈1.82 域外）名字含「先导浆」子串，若子串先判会误落泥浆式
      （τy 差 ~2.2×、T 被误 clamp 到 [40,80]，Q2b 隔离液式在 2D 生产路径失效）；
      凡 role=SPACER 一律走隔离液族（密度档→1.95/2.05/插值/域外就近借+审计）；
    - 行2 名字子串「先导浆/平衡液」只对**非 SPACER** role 生效：单独先导浆
      （role=WASH）、平衡液（role=WASH/MUD）照旧泥浆式+model_assumption；
    - hu103「平衡液」role=MUD（轻泥浆，20313.doc）：仍按行2 带 model_assumption
      标注（公式同为钻井液式，差异仅在审计标注；role 不吞掉名字档）；
    - 行1「钻井液」（role=MUD，八井通用）名字不含行2 子串 → 无标注；
    - 其余 role=MUD 相同理：名字无「先导浆/平衡液」即按行1 无标注。
    """
    name = fluid.name
    role = fluid.role
    if role == FluidRole.FLUSHER or "冲洗" in name:
        return "no_replace", "flusher"
    if role == FluidRole.SPACER:
        return "spacer", None                  # 行4/5：role 优先于名字子串（C1）
    if any(part in name for part in _ASSUMPTION_NAME_PARTS):
        return "mud", "model_assumption"        # 行2：先导浆、平衡液（含别名，非 SPACER）
    if role == FluidRole.MUD or name == "钻井液":
        return "mud", None                      # 行1：八井通用，无标注
    if role == FluidRole.DISPLACEMENT or name in _CHAIN_NAMES:
        return "no_replace", "displacement_chain"  # 压塞液/替钻井液/井浆/基液/保护液…
    if role in (FluidRole.LEAD, FluidRole.TAIL, FluidRole.INTERMEDIATE):
        return "cement", None                   # 中间浆并入水泥档（补裁③）
    return "no_replace", "unmatched"


# ---------------------------------------------------------------------------
# 各族求值（返回 (τy [Pa], μp [Pa·s])；域越界已 clamp+审计）
# ---------------------------------------------------------------------------
def _mud_values(T: float, extrapolate: bool, name: str,
                p: RheologyFormulaParams) -> tuple[float, float]:
    if T < p.mud_T_lo or T > p.mud_T_hi:
        if extrapolate:
            _record("extrapolate", name,
                    f"T={T} 超出实测域 [{p.mud_T_lo:g}, {p.mud_T_hi:g}]，"
                    "按公式外推（mud_extrapolate=True）")
            return _mud_tauy(T, p), _mud_mup(T, p)
        T = _clamp_record(T, p.mud_T_lo, p.mud_T_hi, name, "T")
    return _mud_tauy(T, p), _mud_mup(T, p)


def _spacer_values(fluid: FluidSpec, T: float, P_mpa: Optional[float],
                   p: RheologyFormulaParams) -> tuple[float, float]:
    name = fluid.name
    d = fluid.density_kg_m3 / 1000.0  # g/cm³

    # P：缺省=常压 0.1 MPa（温度阶段；压力模块接入后由调用方传 P(z,t)）
    if P_mpa is None:
        P = p.sp_p_default
        _record("p_default", name, f"P 未传入，按常压 {p.sp_p_default} MPa 求值")
    else:
        P = float(P_mpa)

    T = _clamp_record(T, p.sp_T_lo, p.sp_T_hi, name, "T")
    P = _clamp_record(P, p.sp_P_lo, p.sp_P_hi, name, "P")

    # 密度分派（Q16 就近取，2026-10-06）：整式取最近锚点，不做系数插值。
    # 捕获区左闭右开（平局 d==sp_mid 取高密度端 2.05 式）。
    if d < p.sp_mid:
        ty_c, mp_c = p.sp_1p95_ty, p.sp_1p95_mp
        anchor = p.sp_anchor_lo
    else:
        ty_c, mp_c = p.sp_2p05_ty, p.sp_2p05_mp
        anchor = p.sp_anchor_hi
    if not (p.sp_lo <= d <= p.sp_hi):
        # 远（捕获区外）：维持现状——就近借用端点整式 + borrow 审计（处置待裁）
        _record("borrow", name,
                f"隔离液密度 {d:.3f} g/cm³ 超出就近取捕获区"
                f" [{p.sp_lo},{p.sp_hi}]，就近借用端点公式"
                f"（model_assumption，系数不外推）",
                note="model_assumption")
    elif abs(d - anchor) > 1e-12:
        # 近但不精确命中锚点：记 nearest（供口径对照计数；锚点命中不记，避免逐步刷屏）
        _record("nearest", name,
                f"隔离液密度 {d:.3f} g/cm³ 就近取 {anchor} 式（不插值）",
                requested=d, anchor=anchor)

    return _quad6(ty_c, T, P), _quad6(mp_c, T, P)


def _cement_values(fluid: FluidSpec, T: float,
                   p: RheologyFormulaParams) -> Optional[tuple[float, float]]:
    """按密度档求值（Q16 就近取）；密度在捕获区外返回 None（调用方记不替换审计）。

    难点口径（2026-10-06 Q16）：**整式就近取**锚点公式，不再做组 B↔组 A 密度插值；
    每个锚点用**自己的**温区——组 B [20,200]°C、组 A [20,170]°C（原插值支统一取
    交集 [20,170] 的口径随之取消）。
    """
    name = fluid.name
    d = fluid.density_kg_m3 / 1000.0
    if not (p.cm_lo <= d <= p.cm_hi):
        return None
    if d < p.cm_mid:                            # 组 B（1.90 g/cm³，温区 [20,200]）
        T = _clamp_record(T, p.cb_T_lo, p.cb_T_hi, name, "T")
        if abs(d - p.cm_anchor_b) > 1e-12:
            _record("nearest", name,
                    f"水泥密度 {d:.3f} g/cm³ 就近取 {p.cm_anchor_b:.2f}（组 B）式（不插值）",
                    requested=d, anchor=p.cm_anchor_b)
        return _cement_b_tau0(T, p), _cement_b_mup(T, p)
    # 组 A（2.10 g/cm³，温区 [20,170]）
    T = _clamp_record(T, p.ca_T_lo, p.ca_T_hi, name, "T")
    if abs(d - p.cm_anchor_a) > 1e-12:
        _record("nearest", name,
                f"水泥密度 {d:.3f} g/cm³ 就近取 {p.cm_anchor_a:.2f}（组 A）式（不插值）",
                requested=d, anchor=p.cm_anchor_a)
    return _cement_a_tau0(T, p), _cement_a_mup(T, p)


# ---------------------------------------------------------------------------
# 绝对替换派生
# ---------------------------------------------------------------------------
def _derive(fluid: FluidSpec, tauy: float, mup: float) -> FluidSpec:
    """dataclasses.replace 派生 Bingham 版本：τy→yield_stress_pa、μp→PV、
    幂律参数清空、密度/名字/角色不变（frozen，原对象不动）。"""
    return replace(
        fluid,
        rheology_model=RheologyModel.BINGHAM,
        plastic_viscosity_pa_s=float(mup),
        yield_stress_pa=float(tauy),
        power_law_n=None,
        consistency_k=None,
    )


# ---------------------------------------------------------------------------
# 公共接口
# ---------------------------------------------------------------------------
def fluid_at(
    fluid: FluidSpec,
    T_c: float,
    P_mpa: Optional[float] = None,
    *,
    mud_extrapolate: bool = False,
    smooth_break: bool = False,
    params: Optional[RheologyFormulaParams] = None,
) -> FluidSpec:
    """按温度（及压力）派生该相的流变参数（绝对替换口径）。

    参数
    ----
    fluid : 目标相（按分派表路由；不替换的相原样返回同一对象）
    T_c : 温度 [°C]
    P_mpa : 压力 [MPa]，仅隔离液消费；None=常压 0.1 MPa（记 p_default 审计）
    mud_extrapolate : True 时钻井液越出 [40,80]°C 按公式外推（对比用）
    smooth_break : 断点平滑开关占位——默认 False=忠实两段式；
                   True 尚未实现，抛 NotImplementedError
    params : 公式系数/域界（R1，2026-10-06）。``None`` ⇒ 模块级默认参单例
             ``_DEFAULT_PARAMS``，与 2026-10-01 逐字录入的公式**逐位一致**；
             扰动一律显式传入（`dataclasses.replace` 派生），默认路径永不吃扰动。

    返回
    ----
    派生的新 FluidSpec（Bingham：yield_stress_pa=τy、plastic_viscosity_pa_s=μp，
    power_law_n/consistency_k 清空），或不替换相的原对象。
    """
    if smooth_break:
        raise NotImplementedError(
            "smooth_break 为占位开关（默认 False=忠实两段式；True 尚未实现）"
        )
    p = _DEFAULT_PARAMS if params is None else params
    T = float(T_c)
    family, note = _route(fluid)

    if family == "no_replace":
        _record("no_replace", fluid.name, _NO_REPLACE_DETAIL[note])
        return fluid

    if family == "mud":
        tauy, mup = _mud_values(T, mud_extrapolate, fluid.name, p)
        if note == "model_assumption":
            _record("model_assumption", fluid.name, "暂按钻井液公式（model_assumption）")
        return _derive(fluid, tauy, mup)

    if family == "spacer":
        tauy, mup = _spacer_values(fluid, T, P_mpa, p)
        return _derive(fluid, tauy, mup)

    # cement（LEAD/TAIL/INTERMEDIATE）
    values = _cement_values(fluid, T, p)
    if values is None:
        _record("no_replace", fluid.name,
                f"水泥相密度 {fluid.density_kg_m3 / 1000.0:.3f} g/cm³ 不匹配任何档（不替换）")
        return fluid
    tauy, mup = values
    return _derive(fluid, tauy, mup)
