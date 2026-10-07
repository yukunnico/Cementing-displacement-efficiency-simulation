# -*- coding: utf-8 -*-
"""静液柱压力场数据层（P-1：温压耦合压力侧第一步）

把井筒**静液柱压力**封装成 ``P(md_m, t_s) -> MPa`` 查询对象，供求解器每步取值喂
`cemdisp.data.rheology_vs_temperature.fluid_at`（隔离液族的含压二次曲面；水泥/钻井液
公式纯温度，不吃 P —— 见 `_spacer_values` vs `_cement_values`/`_mud_values`）。

三类场
------
- ``PressureField``：协议（``P(md_m, t_s) -> MPa``）
- ``ConstantPressureField``：常数场（关2 常数场自洽 / 对照档）
- ``HydrostaticPressureField``：静液柱压力场（**TVD 口径** + **在场相密度摊派**）

口径（用户裁定 Q9/D5，2026-10-05；素材C §2.3）
----------------------------------------------
① **静压先行**：只到静液柱压力，**不含**循环压耗 / MPD 回压 / ECD（属 Phase P-2..P-4）；
② 取 P 的深度：由**调用方**决定（求解器侧 ``pressure_caliber ∈ {"shoe","mean"}``）；
③ **在场相密度**：:func:`insitu_column_density` = 各在场相按**设计泵注体积加权**的均值
   （近似口径——见该函数 docstring 的口径声明与已知偏差）；
④ **TVD**：``TVD(md) = ∫₀^md cosθ dmd'``，θ 取自 ``well.inclination_profile``
   （与 loader 读的 ``inclination_profile.csv`` **同源**）；剖面首点以上按 **0°（直井）**
   外推（三重点井该段实测误差 ≤0.05%）。

参照实现（对齐目标）
--------------------
- MATLAB ``参考文档/温压耦合数据、/HT1-004压力计算/p_jaifang1.m:701-715``：
  ``P(i) = P(i-1) + 9.81·ρ(t,i)·Δz_垂直(i)``（**变密度逐段累加**，g=9.81）；
- 本仓停用绘图段 ``scripts/plots/plot_jieti_timing_figs_20260928.py:353-466``：
  ``p(z) = ∫ρ·g dz``（**g=9.80665**）。

本模块默认 ``g = G_STANDARD = 9.80665``（SI 标准 + 本仓既有先例），**可注入**；
两种 g 的差 0.03%（呼1-003 鞋深：145.63 vs 145.68 MPa）。

审计
----
数据层同型三件套 ``oob_count`` / ``oob_events`` / ``reset_audit``。静压无插值 ⇒
越界只记"查询深度超出井底/负深度"，按端点延拓返回，同时记一条事件
（与 ``ConstantTemperatureField`` 的"无越界概念 ⇒ 恒 0"先例同型）。
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Protocol, Sequence, Tuple, runtime_checkable

import numpy as np

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.pumping_schedule import PumpingSchedule
from cemdisp.data.well_spec import WellSpec

__all__ = [
    "PressureField",
    "ConstantPressureField",
    "HydrostaticPressureField",
    "TablePressureField",
    "PressureOutOfRangeEvent",
    "insitu_column_density",
    "ANNULUS_ROLES",
    "G_STANDARD",
]

G_STANDARD: float = 9.80665
"""默认重力加速度 [m/s²]（SI 标准值；与 ``plot_jieti`` 同源，与 MATLAB 的 9.81 差 0.03%）。"""

# 环空柱的流体角色（C-11 口径声明用；`insitu_column_density` 的 `roles=` 参数）
ANNULUS_ROLES: Tuple[FluidRole, ...] = (
    FluidRole.MUD, FluidRole.WASH, FluidRole.SPACER,
    FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL,
)


@runtime_checkable
class PressureField(Protocol):
    """压力场协议：``P(md_m, t_s) -> MPa``（相对井口，静液柱压力）。"""

    def P(self, md_m: float, t_s: float) -> float:  # pragma: no cover - 协议
        ...


@dataclass(frozen=True)
class PressureOutOfRangeEvent:
    """越界查询事件（md 超出 [0, bottom_md_m]）。"""

    md_m: float
    t_s: float
    bottom_md_m: float


class ConstantPressureField:
    """常数压力场（P ≡ 常数，与 md/t 无关）——关2 常数场自洽 / 对照档用。"""

    def __init__(self, p_mpa: float, *, bottom_md_m: float = 0.0) -> None:
        self.p_mpa: float = float(p_mpa)
        self.bottom_md_m: float = float(bottom_md_m)
        self._oob: list = []

    def P(self, md_m: float, t_s: float) -> float:
        md = float(md_m)
        if md < 0.0 or (self.bottom_md_m > 0.0 and md > self.bottom_md_m):
            self._oob.append(PressureOutOfRangeEvent(md, float(t_s), self.bottom_md_m))
        return self.p_mpa

    @property
    def oob_count(self) -> int:
        return len(self._oob)

    @property
    def oob_events(self) -> Tuple[PressureOutOfRangeEvent, ...]:
        return tuple(self._oob)

    def reset_audit(self) -> None:
        self._oob.clear()

    def __repr__(self) -> str:  # pragma: no cover - 仅调试可读性
        return f"ConstantPressureField(p_mpa={self.p_mpa!r})"


class HydrostaticPressureField:
    """静液柱压力场：``P(md,t) = g·∫₀^{TVD(md)} ρ dTVD'``（静压与 t 无关 ⇒ t 被忽略）。

    密度剖面 = **分层**（``((top_md_m, rho_kg_m3), ...)``，首层必须自 0.0 起）：
    第 k 层密度在其 ``top_md`` 之下、下一层 ``top_md`` 之上生效。层数=1 即均匀柱。

    TVD 由 ``well.inclination_profile`` 分段线性积分 ``∫cosθ dmd`` 得到；剖面首点以上
    按 0°（直井）外推。
    """

    def __init__(
        self,
        well_name: str,
        md_grid: Sequence[float],
        tvd_grid: Sequence[float],
        layers: Sequence[Tuple[float, float]],
        *,
        bottom_md_m: float,
        g: float = G_STANDARD,
    ) -> None:
        if not layers:
            raise ValueError("layers 不得为空")
        if float(layers[0][0]) != 0.0:
            raise ValueError(f"layers 首层必须自 0.0 m 起，实际 {layers[0][0]!r}")
        self.well_name = str(well_name)
        self._md_grid = np.asarray(md_grid, dtype=float)
        self._tvd_grid = np.asarray(tvd_grid, dtype=float)
        if self._md_grid.ndim != 1 or self._md_grid.size < 2:
            raise ValueError("md/tvd 网格须为长度 ≥2 的一维序列")
        if np.any(np.diff(self._md_grid) <= 0.0) or np.any(np.diff(self._tvd_grid) < 0.0):
            raise ValueError("md 网格须严格升序；tvd 网格须非降")
        self.layers: Tuple[Tuple[float, float], ...] = tuple(
            (float(top), float(rho)) for top, rho in layers
        )
        self.bottom_md_m: float = float(bottom_md_m)
        self.g: float = float(g)
        self._oob: list = []

    # ---------------------------------------------------------------- 构造器
    @staticmethod
    def _tvd_axis(well: WellSpec) -> Tuple[np.ndarray, np.ndarray]:
        """由 ``well.inclination_profile`` 构 ``(md 轴, 累计 TVD)``；首点以上按 0°。"""
        pts = well.inclination_profile
        if not pts:
            raise ValueError(f"{well.well_name}: 无 inclination_profile，无法做 TVD 口径")
        md = np.array([float(p.depth_md_m) for p in pts], dtype=float)
        inc = np.deg2rad(np.array([float(p.value) for p in pts], dtype=float))
        grid = np.concatenate(([0.0], md))
        cosv = np.concatenate(([1.0], np.cos(inc)))
        seg = 0.5 * (cosv[1:] + cosv[:-1]) * np.diff(grid)
        cum = np.concatenate(([0.0], np.cumsum(seg)))
        return grid, cum

    @classmethod
    def from_well(cls, well: WellSpec, rho_kg_m3: float, *,
                  g: float = G_STANDARD) -> "HydrostaticPressureField":
        """均匀密度柱（``insitu_column_density`` 的典型输出）。"""
        grid, cum = cls._tvd_axis(well)
        return cls(well.well_name, grid, cum, ((0.0, float(rho_kg_m3)),),
                   bottom_md_m=well.bottom_md_m, g=g)

    @classmethod
    def from_layers(cls, well: WellSpec, layers: Sequence[Tuple[float, float]], *,
                    g: float = G_STANDARD) -> "HydrostaticPressureField":
        """分层密度柱（P-2 细化口径的预留入口）。"""
        grid, cum = cls._tvd_axis(well)
        return cls(well.well_name, grid, cum, layers,
                   bottom_md_m=well.bottom_md_m, g=g)

    # ---------------------------------------------------------------- 查询
    def _tvd(self, md_m: float) -> float:
        return float(np.interp(md_m, self._md_grid, self._tvd_grid))

    def P(self, md_m: float, t_s: float) -> float:
        """静液柱压力 [MPa]（t 被忽略——静压口径与时间无关）。"""
        md = float(md_m)
        if md < 0.0 or md > self.bottom_md_m:
            self._oob.append(PressureOutOfRangeEvent(md, float(t_s), self.bottom_md_m))
        if md <= 0.0:
            return 0.0
        total = 0.0
        for k, (top, rho) in enumerate(self.layers):
            lo = top
            hi = self.layers[k + 1][0] if k + 1 < len(self.layers) else md
            hi = min(hi, md)
            if hi <= lo:
                break
            total += rho * (self._tvd(hi) - self._tvd(lo))
        return self.g * total / 1.0e6

    @property
    def oob_count(self) -> int:
        return len(self._oob)

    @property
    def oob_events(self) -> Tuple[PressureOutOfRangeEvent, ...]:
        return tuple(self._oob)

    def reset_audit(self) -> None:
        self._oob.clear()

    def __repr__(self) -> str:  # pragma: no cover - 仅调试可读性
        return (f"HydrostaticPressureField({self.well_name!r}, "
                f"{len(self.layers)} 层, g={self.g})")


class TablePressureField:
    """(md, t) 二维表压力场（Phase 3.1 P-2 二期，spec §1.3）：``P(md_m, t_s) -> MPa`` 双线性插值。

    数据 = ``md 节点列``（m，严格升序）× ``时间列``（s，严格升序）的二维表
    ``table_mpa``（形状 ``(n_md, n_t)``，单位 MPa）。用途 = Phase 4 逐深温度/压力产品位；
    **本阶段只落地 + 单测，不注入任何生产路径**（C-14 同型声明）。

    插值口径
    --------
    先对相邻两时间列在 md 方向各做线性插值（``np.interp``），再在 t 方向线性加权
    （= 矩形网格双线性）。查询点落在网格/列上时返回精确节点值。

    oob 三件套**同型 HydrostaticPressureField**（spec §1.3）
    --------------------------------------------------------
    - 计数条件与静液柱类逐字一致：**只统计 md 轴**（``md<0`` 或 ``md>bottom_md_m``），
      记 :class:`PressureOutOfRangeEvent`；``t`` 轴不属"井底/负深度"语义 ⇒ 端点延拓
      **钳位不计数**（对位 Hydrostatic"静压与 t 无关 ⇒ t 被忽略"——本类钳位而非外推）。
    - 返回语义 = 钳位端点延拓（两轴都不外推）。
    """

    def __init__(
        self,
        well_name: str,
        md_grid_m: Sequence[float],
        t_grid_s: Sequence[float],
        table_mpa: Sequence[Sequence[float]],
        *,
        bottom_md_m: Optional[float] = None,
    ) -> None:
        self._md = np.asarray(md_grid_m, dtype=float)
        self._t = np.asarray(t_grid_s, dtype=float)
        self._tab = np.asarray(table_mpa, dtype=float)
        if self._md.ndim != 1 or self._t.ndim != 1 or self._tab.ndim != 2:
            raise ValueError("md/t 网格须为一维，table 须为二维 (n_md, n_t)")
        if self._md.size < 2 or self._t.size < 2:
            raise ValueError("md/t 网格须为长度 ≥2 的一维序列")
        if self._tab.shape != (self._md.size, self._t.size):
            raise ValueError(
                f"table 形状 {self._tab.shape} 与网格 {(self._md.size, self._t.size)} 不符")
        if np.any(np.diff(self._md) <= 0.0):
            raise ValueError("md 网格须严格升序")
        if np.any(np.diff(self._t) <= 0.0):
            raise ValueError("t 网格须严格升序")
        self.well_name = str(well_name)
        self.bottom_md_m = float(bottom_md_m) if bottom_md_m is not None else float(self._md[-1])
        self._oob: list = []

    def P(self, md_m: float, t_s: float) -> float:
        """双线性插值压力 [MPa]；md 越界计数并钳位，t 钳位不计数（口径见类 docstring）。"""
        md = float(md_m)
        t = float(t_s)
        if md < 0.0 or md > self.bottom_md_m:
            self._oob.append(PressureOutOfRangeEvent(md, t, self.bottom_md_m))
        md_c = min(max(md, float(self._md[0])), float(self._md[-1]))
        t_c = min(max(t, float(self._t[0])), float(self._t[-1]))
        j = int(np.searchsorted(self._t, t_c, side="right")) - 1
        j = min(max(j, 0), self._t.size - 2)
        p0 = float(np.interp(md_c, self._md, self._tab[:, j]))
        p1 = float(np.interp(md_c, self._md, self._tab[:, j + 1]))
        w = (t_c - float(self._t[j])) / (float(self._t[j + 1]) - float(self._t[j]))
        return p0 + w * (p1 - p0)

    @property
    def oob_count(self) -> int:
        return len(self._oob)

    @property
    def oob_events(self) -> Tuple[PressureOutOfRangeEvent, ...]:
        return tuple(self._oob)

    def reset_audit(self) -> None:
        self._oob.clear()

    def __repr__(self) -> str:  # pragma: no cover - 仅调试可读性
        return (f"TablePressureField({self.well_name!r}, "
                f"{self._md.size}×{self._t.size}, bottom_md_m={self.bottom_md_m})")


# ---------------------------------------------------------------------------
# 在场相密度（口径③）
# ---------------------------------------------------------------------------
def insitu_column_density(
    fluids: Sequence[FluidSpec],
    schedule: Optional[PumpingSchedule] = None,
    *,
    roles: Optional[Sequence[FluidRole]] = None,
    fallback_kg_m3: Optional[float] = None,
) -> float:
    """**在场相密度按设计泵注体积加权**的柱平均密度 [kg/m³]（P-1 口径③）。

    参数
    ----
    fluids : 该井流体清单（取各相 ``density_kg_m3`` —— 在场常数，不做 ρ(T,P)）
    schedule : 泵注程序；``None`` 或无可匹配步骤 ⇒ 退化为 ``fluids`` 的**等权**均值
    roles : 参与摊派的角色白名单。``None``（默认）= **全部步骤**（口径见下）；
            传 :data:`ANNULUS_ROLES` 则只算环空柱（排除管内替浆/压塞液/冲洗液）
    fallback_kg_m3 : 全无可用步骤且 roles 过滤后为空时的兜底密度
        （``None`` ⇒ 抛 ``ValueError``）

    口径声明（**须追认**，见 spec §3 Δ9 / §5 C-11）
    ----------------------------------------------
    默认口径 = **全泵注体积加权**：把该井设计泵注的**所有**相按体积权重平均，作为
    "井筒内在场相"的代表密度。它是**近似**——泵注程序里含**管内**流体（替钻井液、
    压塞液、基液；按本仓工艺裁定"替浆不进环空"），严格说是两种物理量的混合。

    选择该口径的理由（实测，2026-10-06）：

    - 是三重点井**唯一**可无损实现的口径——"环空柱"口径需要环空几何才能定出各相在
      环空中的长度，本阶段（P-1）不做；且多数井的泵注程序**没有**独立的钻井液步骤
      （ht1_003/ht1_004 均无），环空口径反而会漏掉在环空占大头的钻井液。
    - 现场锚点吻合：呼101 P(shoe)=**152.34 MPa** vs 现场记录环空静压 **152.16 MPa**
      ⇒ **+0.12%**；呼1-003 P(shoe)=**145.63 MPa** vs 设计锚 **145.96 MPa** ⇒ **−0.23%**。
    - 两口径差异极小：呼1-003 全泵注 1949.53 vs 环空角色 1949.17（**0.02%**）。

    ⚠️ 已知偏差：呼1-004 含 1 m³ 基液（ρ=1.02 g/cm³）+ 轻先导浆 ⇒ ρ̄ 被下拉
    （1881.05；剔除基液为 1884.87）；该井本阶段无鞋深压力锚，偏差**如实记录不标定**。
    """
    by_name = {f.name.strip(): f for f in fluids}
    volumes: dict = {}
    if schedule is not None:
        for step in schedule.steps:
            f = by_name.get(str(step.fluid_name).strip())
            if f is None or step.volume_m3 is None:
                continue
            if roles is not None and f.role not in roles:
                continue
            if float(step.volume_m3) <= 0.0:
                continue
            key = f.name.strip()
            volumes[key] = volumes.get(key, 0.0) + float(step.volume_m3)

    if volumes:
        tot = sum(volumes.values())
        acc = sum(v * float(by_name[n].density_kg_m3) for n, v in volumes.items())
        if tot > 0.0:
            return acc / tot

    pool = [f for f in fluids if roles is None or f.role in roles]
    if not pool:
        if fallback_kg_m3 is None:
            raise ValueError("无可匹配的泵注步骤且 roles 过滤后无流体；请给 fallback_kg_m3")
        return float(fallback_kg_m3)
    return float(np.mean([float(f.density_kg_m3) for f in pool]))
