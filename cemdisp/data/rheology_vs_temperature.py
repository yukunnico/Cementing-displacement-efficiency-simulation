"""
温变流变公式全集（T0-2：温压耦合的公式层 + 流体分派器）

把用户裁定的温变流变公式实现为 ``fluid_at(fluid, T_c, P_mpa) -> FluidSpec``
派生器，供求解器 T-on 时整体替换流体参数（**绝对替换**口径：有公式的相，
τy→``yield_stress_pa``、μp→宾汉塑性粘度，模型族切 ``RheologyModel.BINGHAM``，
幂律参数清空；密度/名字/角色不变）。无公式/不替换的相原样返回同一对象。

公式来源（逐字录入，系数与计划 §3 公式全集一致）：
温压耦合改进计划_2026-09-30.md §3 + task-3-brief.md（含 2026-09-30 四项裁定）：

- 组 A 水泥（2.1 g/cm³，20–170°C）/ 组 B 水泥（1.9 g/cm³，20–200°C）
- 中间密度水泥：组 B↔组 A 按密度线性插值（温区交集 [20,170]°C）
- 隔离液 1.95/2.05 含压二次曲面；ρ∈[1.95,2.05] 六系数密度插值（Q2b），
  域外 clamp 到最近端公式+借用审计（系数不外推）
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
from dataclasses import replace
from typing import Optional

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel

__all__ = ["fluid_at", "get_audit", "reset_audit"]

# ---------------------------------------------------------------------------
# 审计（模块级；fluid_at 为自由函数，同型于 temperature_field 的实例审计）
# ---------------------------------------------------------------------------
_AUDIT_MAX = 10000  # 上限：deque 定长，超限丢最旧（review minor④；接线后按步 reset）
_AUDIT: deque = deque(maxlen=_AUDIT_MAX)


def get_audit() -> list[dict]:
    """返回审计事件列表副本（kind ∈ clamp / borrow / no_replace /
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
def _cement_a_tau0(T: float) -> float:
    """τ₀(T)：20≤T≤100 二次式；100<T≤170 一次式（T 已 clamp 进域）。"""
    if T <= 100.0:
        return 0.002793 * T * T - 0.391702 * T + 21.6378
    return 0.624883 * T - 51.7949


def _cement_a_mup(T: float) -> float:
    """μp(T) = 0.53129·e^(−0.012495·T) − 0.00204 [Pa·s]。"""
    return 0.53129 * math.exp(-0.012495 * T) - 0.00204


# ---------------------------------------------------------------------------
# 组 B 水泥（1.9 g/cm³，20–200°C，纯温度）
# ---------------------------------------------------------------------------
def _cement_b_tau0(T: float) -> float:
    """τ₀(T)：20≤T≤100 二次式；100<T≤200 一次式（T 已 clamp 进域）。"""
    if T <= 100.0:
        return 4.433e-4 * T * T - 0.07843 * T + 5.962
    return 0.03206 * T + 1.979


def _cement_b_mup(T: float) -> float:
    """μp(T) = 1.355×10⁻⁵·T² − 6.570×10⁻⁴·T + 0.08215 [Pa·s]。"""
    return 1.355e-5 * T * T - 6.570e-4 * T + 0.08215


# ---------------------------------------------------------------------------
# 隔离液二次曲面（含压力）：系数序 (常数, T, P, T², T·P, P²)
# ---------------------------------------------------------------------------
_SPACER_1P95_TY = (7.379374, -0.001857, 0.0225351,
                   -0.000157779, -0.00000615055, -0.0000993563)
_SPACER_1P95_MP = (0.0912039, -0.000719380, 0.000167466,
                   0.00000277807, -0.00000322380, 0.00000119469)
_SPACER_2P05_TY = (10.787391, -0.0206255, 0.0320167,
                   -0.000247286, 0.000219522, -0.000253669)
_SPACER_2P05_MP = (0.125694, -0.00108616, 0.000323974,
                   0.00000341573, -0.00000316667, 0.000000714405)


def _quad6(c: tuple, T: float, P: float) -> float:
    """六系数二次曲面：c0 + c1·T + c2·P + c3·T² + c4·T·P + c5·P²。"""
    return c[0] + c[1] * T + c[2] * P + c[3] * T * T + c[4] * T * P + c[5] * P * P


# ---------------------------------------------------------------------------
# 钻井液（呼101 实测拟合，纯温度）
# ---------------------------------------------------------------------------
def _mud_mup(T: float) -> float:
    """μp(T) = 0.297154·e^(−0.02310525·T) [Pa·s]。"""
    return 0.297154 * math.exp(-0.02310525 * T)


def _mud_tauy(T: float) -> float:
    """τy(T) = 21.346710·e^(−0.00830631·T) [Pa]。"""
    return 21.346710 * math.exp(-0.00830631 * T)


# ---------------------------------------------------------------------------
# 分派表（裁定定稿，含 2026-09-30 补裁）
# ---------------------------------------------------------------------------
# 行2「先导浆、平衡液」按**名字子串**匹配（不全等）：覆盖真实井况别名
# hu1「平衡液(先导泥浆)」（role=WASH，2026-08-29 由旧名"冲洗液"更名）。
# 子串刻意取窄（「先导浆」「平衡液」），不误伤替浆链/隔离液等其他相。
_ASSUMPTION_NAME_PARTS = ("先导浆", "平衡液")
_CHAIN_NAMES = frozenset({"压塞液"})                  # 替浆链：其余由 role=DISPLACEMENT 覆盖

_NO_REPLACE_DETAIL = {
    "flusher": "冲洗液不替换（返回原 FluidSpec，审计标记）",
    "displacement_chain": "替浆链常数现状（不替换）",
    "unmatched": "不匹配任何档（不替换）",
}


def _route(fluid: FluidSpec) -> tuple[str, Optional[str]]:
    """按分派表返回 (族, 附注)。族 ∈ no_replace / mud / spacer / cement。

    路由优先级裁定（review 修复）：**行2 名字命中优先于 role==MUD**——
    - hu103「平衡液」role=MUD（轻泥浆，20313.doc）：仍按行2 带 model_assumption
      标注（公式同为钻井液式，差异仅在审计标注；role 不吞掉名字档）；
    - 行1「钻井液」（role=MUD，八井通用）名字不含行2 子串 → 无标注；
    - 其余 role=MUD 相同理：名字无「先导浆/平衡液」即按行1 无标注。
    """
    name = fluid.name
    role = fluid.role
    if role == FluidRole.FLUSHER or "冲洗" in name:
        return "no_replace", "flusher"
    if any(part in name for part in _ASSUMPTION_NAME_PARTS):
        return "mud", "model_assumption"        # 行2：先导浆、平衡液（含别名）
    if role == FluidRole.MUD or name == "钻井液":
        return "mud", None                      # 行1：八井通用，无标注
    if role == FluidRole.DISPLACEMENT or name in _CHAIN_NAMES:
        return "no_replace", "displacement_chain"  # 压塞液/替钻井液/井浆/基液/保护液…
    if role == FluidRole.SPACER:
        return "spacer", None
    if role in (FluidRole.LEAD, FluidRole.TAIL, FluidRole.INTERMEDIATE):
        return "cement", None                   # 中间浆并入水泥档（补裁③）
    return "no_replace", "unmatched"


# ---------------------------------------------------------------------------
# 各族求值（返回 (τy [Pa], μp [Pa·s])；域越界已 clamp+审计）
# ---------------------------------------------------------------------------
def _mud_values(T: float, extrapolate: bool, name: str) -> tuple[float, float]:
    if T < 40.0 or T > 80.0:
        if extrapolate:
            _record("extrapolate", name,
                    f"T={T} 超出实测域 [40, 80]，按公式外推（mud_extrapolate=True）")
            return _mud_tauy(T), _mud_mup(T)
        T = _clamp_record(T, 40.0, 80.0, name, "T")
    return _mud_tauy(T), _mud_mup(T)


def _spacer_values(fluid: FluidSpec, T: float, P_mpa: Optional[float]) -> tuple[float, float]:
    name = fluid.name
    d = fluid.density_kg_m3 / 1000.0  # g/cm³

    # P：缺省=常压 0.1 MPa（温度阶段；压力模块接入后由调用方传 P(z,t)）
    if P_mpa is None:
        P = 0.1
        _record("p_default", name, "P 未传入，按常压 0.1 MPa 求值")
    else:
        P = float(P_mpa)

    T = _clamp_record(T, 20.0, 200.0, name, "T")
    P = _clamp_record(P, 0.1, 200.0, name, "P")

    # 密度分派（Q2b）：[1.95,2.05] 六系数插值；端点=专属式；域外就近端借用
    if 1.95 <= d <= 2.05:
        if d == 1.95:
            ty_c, mp_c = _SPACER_1P95_TY, _SPACER_1P95_MP
        elif d == 2.05:
            ty_c, mp_c = _SPACER_2P05_TY, _SPACER_2P05_MP
        else:
            w = (d - 1.95) / 0.10
            ty_c = tuple(a + w * (b - a) for a, b in zip(_SPACER_1P95_TY, _SPACER_2P05_TY))
            mp_c = tuple(a + w * (b - a) for a, b in zip(_SPACER_1P95_MP, _SPACER_2P05_MP))
    else:
        # 最近端：d<1.95 → 1.95 式；d>2.05 → 2.05 式（系数不外推）
        if d < 1.95:
            ty_c, mp_c = _SPACER_1P95_TY, _SPACER_1P95_MP
        else:
            ty_c, mp_c = _SPACER_2P05_TY, _SPACER_2P05_MP
        _record("borrow", name,
                f"隔离液密度 {d:.3f} g/cm³ 不在 [1.95,2.05]，就近借用端点公式"
                f"（model_assumption，系数不外推）",
                note="model_assumption")

    return _quad6(ty_c, T, P), _quad6(mp_c, T, P)


def _cement_values(fluid: FluidSpec, T: float) -> Optional[tuple[float, float]]:
    """按密度档求值；密度不匹配任何档返回 None（调用方记不替换审计）。"""
    name = fluid.name
    d = fluid.density_kg_m3 / 1000.0

    if 1.88 <= d <= 1.92:                       # 尾浆 1.90 档 ±0.02
        T = _clamp_record(T, 20.0, 200.0, name, "T")
        return _cement_b_tau0(T), _cement_b_mup(T)
    if 2.08 <= d <= 2.12:                       # 领/尾浆 2.10 档 ±0.02
        T = _clamp_record(T, 20.0, 170.0, name, "T")
        return _cement_a_tau0(T), _cement_a_mup(T)
    if 1.92 < d < 2.08:                         # 领浆 1.93/1.95/2.05 等 → 密度插值
        T = _clamp_record(T, 20.0, 170.0, name, "T")   # 组 A/B 温区交集
        w = (d - 1.9) / 0.2
        tb, ab = _cement_b_tau0(T), _cement_b_mup(T)
        ta, am = _cement_a_tau0(T), _cement_a_mup(T)
        return tb + w * (ta - tb), ab + w * (am - ab)
    return None


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

    返回
    ----
    派生的新 FluidSpec（Bingham：yield_stress_pa=τy、plastic_viscosity_pa_s=μp，
    power_law_n/consistency_k 清空），或不替换相的原对象。
    """
    if smooth_break:
        raise NotImplementedError(
            "smooth_break 为占位开关（默认 False=忠实两段式；True 尚未实现）"
        )
    T = float(T_c)
    family, note = _route(fluid)

    if family == "no_replace":
        _record("no_replace", fluid.name, _NO_REPLACE_DETAIL[note])
        return fluid

    if family == "mud":
        tauy, mup = _mud_values(T, mud_extrapolate, fluid.name)
        if note == "model_assumption":
            _record("model_assumption", fluid.name, "暂按钻井液公式（model_assumption）")
        return _derive(fluid, tauy, mup)

    if family == "spacer":
        tauy, mup = _spacer_values(fluid, T, P_mpa)
        return _derive(fluid, tauy, mup)

    # cement（LEAD/TAIL/INTERMEDIATE）
    values = _cement_values(fluid, T)
    if values is None:
        _record("no_replace", fluid.name,
                f"水泥相密度 {fluid.density_kg_m3 / 1000.0:.3f} g/cm³ 不匹配任何档（不替换）")
        return fluid
    tauy, mup = values
    return _derive(fluid, tauy, mup)
