"""
温度场表数据层（T0-1：温压耦合的温度查询模块）

把交付的 333×200 温度场表封装成 ``T(md_m, t_s)`` 查询对象，供后续求解器每步取温：

- ``ConstantTemperatureField``: 恒温场（关温耦合时的回退默认）
- ``GeothermalTemperatureField``: 地温静温剖面场 T(z)=T0+grad·z（无瞬态表井的静温档）
- ``AnchoredProfileField``: 静温点锚分段线性场（4d 新增）——循环准稳态温度系数
  k 必给无默认 + Ramey 型瞬态接口位 + 锚域外回退地温式（F-4 聚合越界计数）
- ``TableTemperatureField``: 表格温度场，双线性插值（深度×时间），越界 clamp + 审计
- ``load_delivered_pair``: 加载交付的管内/环空两表，并校验井底（最深行）共享同值
- ``load_extended_pair_4d``: 呼1-004 时程扩展表（333x362，表末 21600 s）的版本感知
  加载（Phase 4d，F-4 登记）；``assert_time_table_coverage``: 批次级 stop_t>表末 前置断言

数据源（交付件）：
- 管内表: 参考文档/温压耦合数据、/T_in.xlsx（333 行 × 200 列纯数值矩阵，无表头）
- 环空表: 参考文档/温压耦合数据、/T_out.xlsx（同上；井底行与 T_in 共享同值）
- 深度轴: 参考文档/温压耦合数据、/HT1-004压力计算/呼1-004井身结构.csv
  第 2 列 ``depth_well_logging_m_``（333 点，30→7660 m）
- 时间轴: 200 列 = 200 min，col0 = 初始时刻，t_j = j min × 60 s（裁定总表）

xlsx 读取不依赖 openpyxl（cementT 环境无该包）：直接用标准库 zipfile + xml 解析
首个工作表的数值单元格。读表后缓存为 npz（默认落在 xlsx 同目录，可用
``cache_path`` 指定），二次加载免解析 xlsx。
"""

from __future__ import annotations

import re
import zipfile
from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd

__all__ = [
    "ConstantTemperatureField",
    "GeothermalTemperatureField",
    "AnchoredProfileField",
    "TableTemperatureField",
    "ClampEvent",
    "load_delivered_pair",
    "load_extended_pair_4d",
    "assert_time_table_coverage",
    "TemperatureTableCoverageError",
    "EXTENDED_SHAPE_4D",
    "EXTENDED_TABLE_MARKER",
    "TABLE_VERSIONS",
    "GEO_T0_C",
    "GEO_GRAD_C_PER_M",
    "EXPECTED_SHAPE",
    "TIME_STEP_S",
    "DEFAULT_T_IN_XLSX",
    "DEFAULT_T_OUT_XLSX",
    "DEFAULT_DEPTH_CSV",
]

# 仓库根：cemdisp/data/temperature_field.py → parents[2]
PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA_DIR = PROJECT_ROOT / "参考文档" / "温压耦合数据、"
DEFAULT_T_IN_XLSX = DEFAULT_DATA_DIR / "T_in.xlsx"          # 管内
DEFAULT_T_OUT_XLSX = DEFAULT_DATA_DIR / "T_out.xlsx"         # 环空
DEFAULT_DEPTH_CSV = (
    DEFAULT_DATA_DIR / "HT1-004压力计算" / "呼1-004井身结构.csv"
)

# 越界审计事件列表上限（Phase 5b-①：deque 定长同型化，与
# `rheology_vs_temperature._AUDIT_MAX` 同一策略——超限丢最旧，内存有界）。
# 语义说明：诊断事件列表是**观测窗**而非完整日志；批量列路径本就"每批至多 1 条"，
# 定长只影响极端越界场景下的保留窗口，不影响 `oob_column_clamped_total`（累加计数不丢）。
_OOB_AUDIT_MAX = 10000
_CACHE_VERSION = 2   # npz 缓存版本（Phase 5b-②：加指纹后升版）
EXPECTED_SHAPE: Tuple[int, int] = (333, 200)  # 深度 333 点 × 时间 200 列
# ---- Phase 4d 版本感知登记表（交付路径逐位不变；扩展表须经 load_extended_pair_4d）----
# 扩展表 = HT1_004_T.m 沙箱扩时程产物：0..198 min 逐分钟 + 施工终点 198.792801 min
# + 199..360 min 逐分钟，共 362 列，表末 21600 s ≥ stop_t(r0.6)=19016.0 s。
# 前 200 列与交付表逐元素位级一致（对账证据：results/_probe_4d扩表_2026-10-07/compare_result.json）。
EXTENDED_SHAPE_4D: Tuple[int, int] = (333, 362)
EXTENDED_TABLE_MARKER = "_ext4d"  # 扩展表文件名必备溯源标识
TABLE_VERSIONS = {"delivered": EXPECTED_SHAPE, "extended_4d": EXTENDED_SHAPE_4D}
TIME_STEP_S: float = 60.0                     # 每列 1 min，col0=初始时刻
DEPTH_CSV_COLUMN = "depth_well_logging_m_"
_SHEET_XML = "xl/worksheets/sheet1.xml"
_CELL_REF_RE = re.compile(r"([A-Z]+)(\d+)")


# ---------------------------------------------------------------------------
# 审计记录
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ClampEvent:
    """一次越界查询的审计记录（记录请求位置与 clamp 后的落点）。"""

    md_m: float          # 请求深度
    t_s: float           # 请求时间
    md_clamped_m: float  # clamp 后深度
    t_clamped_s: float   # clamp 后时间


# ---------------------------------------------------------------------------
# 恒温场
# ---------------------------------------------------------------------------
class ConstantTemperatureField:
    """恒温场：任意 (md, t) 返回常数 T_c（°C）。关温耦合时的回退默认。"""

    def __init__(self, T_c: float) -> None:
        self.T_c = float(T_c)

    def T(self, md_m: float, t_s: float) -> float:
        return self.T_c

    @property
    def oob_count(self) -> int:
        """越界次数（恒温场永不越界，恒为 0，与 TableTemperatureField 同型）。"""
        return 0

    @property
    def oob_events(self) -> Tuple[ClampEvent, ...]:
        return ()

    def reset_audit(self) -> None:
        """清空审计（恒温场无审计，保留同型接口）。"""


# ---------------------------------------------------------------------------
# 地温静温剖面（统一式，单一真源）
# ---------------------------------------------------------------------------
# 温压耦合改进计划_2026-09-30 §1 裁定：静温剖面统一取线性地温线
# T(z) = 16.006 + 1.7598e-2·z（°C，z 单位 m）。对呼1-004 交付温度表首列
# （col0 = 初始时刻）拟合的最大残差 ~0.04–0.05 °C。
# 呼101/呼103 无瞬态温度表 ⇒ T-on 静温档用本剖面（变体备注打「无瞬态表」）。
# ⚠️ 本两常量是全仓唯一真源：其它脚本一律 import，不得就地复制字面量。
GEO_T0_C = 16.006             # °C，地表（z=0）截距
GEO_GRAD_C_PER_M = 1.7598e-2  # °C/m，地温梯度


class GeothermalTemperatureField:
    """地温静温剖面场：``T(md_m, t_s) = GEO_T0_C + GEO_GRAD_C_PER_M·md_m``（°C）。

    与 :class:`ConstantTemperatureField` **完全同型**的零变查询接口——
    只消费 ``T(md, t) -> float``（时间维不参与计算，静温与时刻无关），
    永不越界（``oob_count`` 恒 0、``oob_events`` 恒 ``()``）、``reset_audit()`` 空操作，
    故可与 Constant / Table 互换注入 ``CasingFlowSolver.run`` / ``AnnulusD2DGASolver.run``。

    默认参数即统一式 ``GEO_T0_C``/``GEO_GRAD_C_PER_M``（模块常量，单一真源）；
    显式传参仅用于测试与敏感性档位构造。
    """

    def __init__(
        self,
        T0_c: float = GEO_T0_C,
        grad_c_per_m: float = GEO_GRAD_C_PER_M,
    ) -> None:
        self.T0_c = float(T0_c)
        self.grad_c_per_m = float(grad_c_per_m)

    def T(self, md_m: float, t_s: float) -> float:
        # t_s 不消费：静温剖面无时间维（与 Constant 同，查询签名保持一致）。
        return self.T0_c + self.grad_c_per_m * float(md_m)

    @property
    def oob_count(self) -> int:
        """越界次数（线性式对任意实数深度有定义，永不越界，恒为 0）。"""
        return 0

    @property
    def oob_events(self) -> Tuple[ClampEvent, ...]:
        return ()

    def reset_audit(self) -> None:
        """清空审计（静温场无审计，保留同型接口）。"""

    def __repr__(self) -> str:  # pragma: no cover —— 仅调试可读性
        return (f"GeothermalTemperatureField(T0_c={self.T0_c!r}, "
                f"grad_c_per_m={self.grad_c_per_m!r})")


# ---------------------------------------------------------------------------
# 静温锚点剖面场（4d 新增，纯追加；Phase 4 设计规格 §1 项 5 / 测绘项 14）
# ---------------------------------------------------------------------------
class AnchoredProfileField:
    """静温点锚分段线性场：``T(md_m, t_s) -> °C``，与既有三场同型可互换注入。

    物理口径（如实声明，无新物理发明）
    ----------------------------------
    - **基准线**：静温族点锚 ``(md_m, T_c)`` 分段线性（``np.interp``），代表
      静止地温剖面。锚点的族归属（静温/循环/出口/邻井四族混装是已知数据现状，
      测绘项 16）由调用方**构造前裁定**——同一 md 出现不同 T 视为族混装，
      构造期直接 ``ValueError``（同 md 同值的跨文档重复行自动合并）。
    - **循环准稳态修正**：``regime="circulating"`` 时返回 ``T = k · T_static(md)``。
      k 即施工设计「领浆/尾浆温度系数」口径——对摄氏绝对值的乘性系数
      （锚点 notes 实证：152 °C × 0.85 = 129.2 °C；155 °C × 0.85 ≈ 131.75 ≈ 132 °C）。
      ``k ∈ (0, 1]`` **必给、无默认**：系数取值 = 用户硬停点（先报后动），
      本类不携带任何缺省系数；k=1 ⇒ 纯静温锚。
    - **锚域外回退**：``fallback_geothermal=True``（默认，已声明）时锚点 md 域外
      回退 :class:`GeothermalTemperatureField` 统一地温式（复用模块常量单一真源，
      不复制字面量）；``False`` 时钳位到端点锚值。域外查询一律记越界审计。
    - **瞬态项（接口位，实现从简）**：``transient_tau_s`` 给定时按一阶集总
      Ramey/Hasan-Kabir 型混合
      ``T(md,t) = T_qs(md) + [T_static(md) − T_qs(md)] · exp(−t/τ)``，
      t=0 出发于静温、t→∞ 收敛于循环准稳态。这是对非稳态井筒热交换的
      **一阶集总近似**（无径向导热/传输线解析细节，τ 为 lumped 时间常数），
      仅作接口位；缺省 ``None`` = 不启用（时间维不消费，与 Geothermal 同型）。

    审计（Phase 4 设计规格 §0 F-4）
    ------------------------------
    三件套 ``oob_count`` / ``oob_events`` / ``reset_audit()`` 与既有场同型：
    标量路径 ``T()`` **逐查询**记 ``ClampEvent``（与 Table 一致）；批量列查询
    :meth:`T_column` **逐批聚合计数**——每批至多追加 1 条代表事件（偏移锚域
    边界最远的域外点），域外点数另累计于 :attr:`oob_column_clamped_total`，
    防逐查询 append 爆表。
    """

    def __init__(
        self,
        anchors: Sequence[Tuple[float, float]],
        temperature_factor_k: float,
        *,
        regime: str = "circulating",
        fallback_geothermal: bool = True,
        transient_tau_s: Optional[float] = None,
        source: Optional[str] = None,
    ) -> None:
        pts = sorted((float(md), float(t)) for md, t in anchors)
        if len(pts) < 2:
            raise ValueError(f"静温锚点至少需 2 点（分段线性），实际 {len(pts)}")
        if not np.isfinite(np.asarray(pts, dtype=float)).all():
            raise ValueError("静温锚点含 NaN/Inf，须先剔除非有限值")
        # 同 md：同值合并（跨文档重复行），异值抛错（族混装须构造前裁定）
        mds: List[float] = []
        ts: List[float] = []
        for md, t in pts:
            if mds and abs(md - mds[-1]) < 1e-9:
                if abs(t - ts[-1]) > 1e-9:
                    raise ValueError(
                        f"锚点冲突：md={md:.1f} m 出现不同静温值 {ts[-1]:.2f} / "
                        f"{t:.2f} °C（温度族混装须在构造前裁定，测绘项16/风险9）"
                    )
                continue
            mds.append(md)
            ts.append(t)
        if len(mds) < 2:
            raise ValueError(f"去重后锚点不足 2 点，无法分段线性：{pts}")

        k = float(temperature_factor_k)
        if not np.isfinite(k) or not (0.0 < k <= 1.0):
            raise ValueError(
                "温度系数 k 须属于 (0, 1]（必给、无默认——取值属用户硬停点"
                f"先报后动），实际 {temperature_factor_k!r}"
            )
        if regime not in ("static", "circulating"):
            raise ValueError(
                "regime 须为 static（纯静温锚）或 circulating（k 乘性修正），"
                f"实际 {regime!r}"
            )
        tau: Optional[float] = None
        if transient_tau_s is not None:
            tau = float(transient_tau_s)
            if not np.isfinite(tau) or tau <= 0.0:
                raise ValueError(
                    f"transient_tau_s 须为正有限秒数，实际 {transient_tau_s!r}"
                )

        self.md_anchor_m = np.asarray(mds, dtype=float)
        self.T_anchor_c = np.asarray(ts, dtype=float)
        self.temperature_factor_k = k
        self.regime = regime
        self.fallback_geothermal = bool(fallback_geothermal)
        self.transient_tau_s = tau
        self.source = source
        # 域外回退线复用统一地温场（单一真源常量，不复制字面量）
        self._geo = GeothermalTemperatureField()
        self._oob_events: deque = deque(maxlen=_OOB_AUDIT_MAX)
        self._oob_column_clamped_total = 0

    # -- 内部：静温基准线（域外按回退开关） ----------------------------------
    def _static_line(self, md: np.ndarray) -> np.ndarray:
        """静温基准线：锚域内分段线性；域外按 fallback 开关取地温式/端点值。"""
        mds = self.md_anchor_m
        vals = np.interp(md, mds, self.T_anchor_c)  # np.interp 域外默认钳端点值
        if self.fallback_geothermal:
            outside = (md < float(mds[0])) | (md > float(mds[-1]))
            if np.any(outside):
                vals = np.where(
                    outside,
                    self._geo.T0_c + self._geo.grad_c_per_m * md,
                    vals,
                )
        return vals

    def _apply_regime(self, t_static: np.ndarray, t_s: float) -> np.ndarray:
        """按 regime/瞬态开关把静温线映射为返回温度（static 原样；circulating 乘 k，
        启用瞬态时再按 exp(-t/tau) 从静温向准稳态一阶混合）。"""
        if self.regime == "static":
            return t_static
        t_qs = self.temperature_factor_k * t_static
        if self.transient_tau_s is None:
            return t_qs
        blend = float(np.exp(-max(t_s, 0.0) / self.transient_tau_s))
        return t_qs + (t_static - t_qs) * blend

    # -- 越界审计三件套（同型）+ F-4 聚合计数 --------------------------------
    @property
    def oob_count(self) -> int:
        """越界审计计数（F-4 口径：标量路径逐查询计数；列批量路径逐批计数，
        每批至多 1 条代表事件，域外点数另见 oob_column_clamped_total）。"""
        return len(self._oob_events)

    @property
    def oob_events(self) -> Tuple[ClampEvent, ...]:
        """越界查询记录（标量路径逐条；批量路径为逐批代表事件）。"""
        return tuple(self._oob_events)

    @property
    def oob_column_clamped_total(self) -> int:
        """批量列查询累计域外（被回退/钳位）查询点数——F-4 聚合计数。"""
        return self._oob_column_clamped_total

    def reset_audit(self) -> None:
        """清空越界审计（事件列表与批量聚合计数）。"""
        self._oob_events.clear()
        self._oob_column_clamped_total = 0

    # -- 查询接口（与既有场同型） --------------------------------------------
    def T(self, md_m: float, t_s: float) -> float:
        """取 (md_m [m], t_s [s]) 处温度（°C）。

        越界语义：md 出锚域 ⇒ 回退线取值 + 逐查询 ClampEvent；
        t<0 仅在瞬态启用时钳 0 并审计（瞬态关闭时时间维不消费，同 Geothermal）。
        """
        md = float(md_m)
        tt = float(t_s)
        mds = self.md_anchor_m
        lo, hi = float(mds[0]), float(mds[-1])
        md_out = md < lo or md > hi
        t_out = self.transient_tau_s is not None and tt < 0.0
        if md_out or t_out:
            self._oob_events.append(
                ClampEvent(
                    md_m=md,
                    t_s=tt,
                    md_clamped_m=min(max(md, lo), hi),
                    t_clamped_s=(max(tt, 0.0) if self.transient_tau_s is not None else tt),
                )
            )
        return float(self._apply_regime(self._static_line(np.array([md])), tt)[0])

    def T_column(self, md_values, t_s: float = 0.0) -> np.ndarray:
        """批量列查询（F-4 聚合计数入口）：返回与输入同形状的 °C float 数组。

        与逐点调用 T **同值**；越界一次计数——每批至多追加 1 条代表
        ClampEvent（取偏移锚域边界最远的域外点），域外点数累加进
        oob_column_clamped_total。既有标量路径不受影响。
        """
        arr = np.asarray(md_values, dtype=float)
        flat = arr.ravel()
        mds = self.md_anchor_m
        lo, hi = float(mds[0]), float(mds[-1])
        tt = float(t_s)
        below = flat < lo
        above = flat > hi
        n_oob_md = int(np.count_nonzero(below)) + int(np.count_nonzero(above))
        t_out = self.transient_tau_s is not None and tt < 0.0
        if n_oob_md or t_out:
            if n_oob_md:
                dist = np.where(below, lo - flat, 0.0) + np.where(above, flat - hi, 0.0)
                rep_md = float(flat[int(np.argmax(dist))])
            else:
                rep_md = float(flat[0]) if flat.size else 0.0
            self._oob_events.append(
                ClampEvent(
                    md_m=rep_md,
                    t_s=tt,
                    md_clamped_m=min(max(rep_md, lo), hi),
                    t_clamped_s=(
                        max(tt, 0.0) if self.transient_tau_s is not None else tt
                    ),
                )
            )
            self._oob_column_clamped_total += n_oob_md
        return self._apply_regime(self._static_line(flat), tt).reshape(arr.shape)

    def __repr__(self) -> str:  # pragma: no cover —— 仅调试可读性
        return (
            f"AnchoredProfileField(n_anchors={len(self.md_anchor_m)}, "
            f"temperature_factor_k={self.temperature_factor_k!r}, "
            f"regime={self.regime!r}, fallback_geothermal={self.fallback_geothermal!r}, "
            f"transient_tau_s={self.transient_tau_s!r}, source={self.source!r})"
        )


# ---------------------------------------------------------------------------
# 表格温度场
# ---------------------------------------------------------------------------
class TableTemperatureField:
    """表格温度场：深度×时间双线性插值，越界 clamp 到表域边界并记审计。

    参数
    ----------
    table : (n_depth, n_time) 纯数值矩阵
    depth_m : 深度轴（m），长度须等于 n_depth，严格单调增
    time_s : 时间轴（s），长度须等于 n_time；缺省按 ``j * TIME_STEP_S`` 生成
    expected_shape : 期望形状，默认交付口径 (333, 200)；仅测试可用小形状
    source : 来源文件路径（可选，仅用于审计/溯源）

    审计（Phase 4 设计规格 §0 F-4）
    ------------------------------
    标量路径 :meth:`T` **逐查询**记 ``ClampEvent``（既有行为，未改动）；
    批量列查询 :meth:`T_column` **逐批聚合计数**——每批至多追加 1 条代表事件
    （偏移表域边界最远的域外深度点），域外点数另累计于
    :attr:`oob_column_clamped_total`，防逐列×逐步 append 爆表。
    """

    def __init__(
        self,
        table: np.ndarray,
        depth_m: np.ndarray,
        time_s: Optional[np.ndarray] = None,
        *,
        expected_shape: Tuple[int, int] = EXPECTED_SHAPE,
        source: Optional[Path] = None,
    ) -> None:
        arr = np.asarray(table, dtype=float)
        if arr.ndim != 2:
            raise ValueError(f"温度表须为二维矩阵，实际 ndim={arr.ndim}")
        if arr.shape != tuple(expected_shape):
            raise ValueError(
                f"温度表形状应为 {tuple(expected_shape)}（深度×时间），"
                f"实际 {arr.shape}"
            )

        z = np.asarray(depth_m, dtype=float)
        if z.ndim != 1:
            raise ValueError(f"深度轴须为一维，实际 ndim={z.ndim}")
        if len(z) != arr.shape[0]:
            raise ValueError(
                f"深度轴长度应为 {arr.shape[0]}，实际 {len(z)}"
            )
        if len(z) < 2:
            raise ValueError(f"深度轴至少需 2 个节点，实际 {len(z)}")
        if np.isnan(z).any():
            raise ValueError(f"深度轴含 NaN（{int(np.isnan(z).sum())} 个）")
        if not np.all(np.diff(z) > 0):
            raise ValueError("深度轴非严格单调递增（须随深度增加）")

        if time_s is None:
            t = np.arange(arr.shape[1], dtype=float) * TIME_STEP_S
        else:
            t = np.asarray(time_s, dtype=float)
            if t.ndim != 1:
                raise ValueError(f"时间轴须为一维，实际 ndim={t.ndim}")
            if len(t) != arr.shape[1]:
                raise ValueError(
                    f"时间轴长度应为 {arr.shape[1]}（=温度表列数），实际 {len(t)}"
                )
            if len(t) < 2:
                raise ValueError(f"时间轴至少需 2 个节点，实际 {len(t)}")
            if np.isnan(t).any():
                raise ValueError(f"时间轴含 NaN（{int(np.isnan(t).sum())} 个）")
            if not np.all(np.diff(t) > 0):
                raise ValueError("时间轴非严格单调递增（须随时间增加）")

        if np.isnan(arr).any():
            raise ValueError(f"温度表含 NaN（{int(np.isnan(arr).sum())} 个）")

        self.table = arr
        self.depth_m = z
        self.time_s = t
        self.source = Path(source) if source is not None else None
        self._oob_events: deque = deque(maxlen=_OOB_AUDIT_MAX)
        self._oob_column_clamped_total: int = 0

    # -- 审计 ---------------------------------------------------------------
    @property
    def oob_count(self) -> int:
        """越界（被 clamp）查询次数（F-4 口径：标量路径逐查询计数；列批量路径
        逐批计数，每批至多 1 条代表事件，域外点数另见 oob_column_clamped_total）。"""
        return len(self._oob_events)

    @property
    def oob_events(self) -> Tuple[ClampEvent, ...]:
        """越界查询记录（标量路径逐条；批量路径为逐批代表事件）。"""
        return tuple(self._oob_events)

    @property
    def oob_column_clamped_total(self) -> int:
        """批量列查询累计域外（被 clamp）查询点数——F-4 聚合计数。"""
        return self._oob_column_clamped_total

    def reset_audit(self) -> None:
        """清空越界审计（事件列表与批量聚合计数）。"""
        self._oob_events.clear()
        self._oob_column_clamped_total = 0

    # -- 查询 ---------------------------------------------------------------
    def T(self, md_m: float, t_s: float) -> float:
        """取 (md_m [m], t_s [s]) 处温度（°C）：双线性插值，越界 clamp+审计。"""
        md = float(md_m)
        tt = float(t_s)
        z = self.depth_m
        ts = self.time_s

        md_c = min(max(md, float(z[0])), float(z[-1]))
        tt_c = min(max(tt, float(ts[0])), float(ts[-1]))
        if md_c != md or tt_c != tt:
            self._oob_events.append(
                ClampEvent(
                    md_m=md,
                    t_s=tt,
                    md_clamped_m=md_c,
                    t_clamped_s=tt_c,
                )
            )

        i = int(np.searchsorted(z, md_c, side="right")) - 1
        i = min(max(i, 0), len(z) - 2)
        j = int(np.searchsorted(ts, tt_c, side="right")) - 1
        j = min(max(j, 0), len(ts) - 2)

        fz = (md_c - z[i]) / (z[i + 1] - z[i])
        ft = (tt_c - ts[j]) / (ts[j + 1] - ts[j])

        t00 = self.table[i, j]
        t01 = self.table[i, j + 1]
        t10 = self.table[i + 1, j]
        t11 = self.table[i + 1, j + 1]
        top = t00 + ft * (t01 - t00)
        bot = t10 + ft * (t11 - t10)
        return float(top + fz * (bot - top))

    def T_column(self, md_values, t_s: float = 0.0) -> np.ndarray:
        """批量列查询（F-4 聚合计数入口）：返回与输入同形状的 °C float 数组。

        与逐点调用 :meth:`T` **逐位同值**（同一插值算式与运算次序，向量化后
        逐元素 IEEE-754 相同）；越界一次计数——每批至多追加 1 条代表
        ``ClampEvent``（取偏离表域边界最远的域外深度点），深度域外点数累加进
        :attr:`oob_column_clamped_total`。既有标量路径不受影响。

        与 :meth:`AnchoredProfileField.T_column` 同型（同一 F-4 聚合口径）；
        供二维逐列温度（`enable_depthwise_temperature`）批量取温使用，避免
        逐步×逐列调用 :meth:`T` 使 ``_oob_events`` 无界增长。
        """
        arr = np.asarray(md_values, dtype=float)
        flat = arr.ravel()
        z = self.depth_m
        ts = self.time_s
        lo_md, hi_md = float(z[0]), float(z[-1])
        lo_t, hi_t = float(ts[0]), float(ts[-1])
        tt = float(t_s)

        md_c = np.clip(flat, lo_md, hi_md)
        tt_c = min(max(tt, lo_t), hi_t)
        below = flat < lo_md
        above = flat > hi_md
        n_oob_md = int(np.count_nonzero(below)) + int(np.count_nonzero(above))
        t_out = tt_c != tt
        if n_oob_md or t_out:
            if n_oob_md:
                dist = np.where(below, lo_md - flat, 0.0) + np.where(above, flat - hi_md, 0.0)
                rep_md = float(flat[int(np.argmax(dist))])
            else:
                rep_md = float(flat[0]) if flat.size else 0.0
            self._oob_events.append(
                ClampEvent(
                    md_m=rep_md,
                    t_s=tt,
                    md_clamped_m=min(max(rep_md, lo_md), hi_md),
                    t_clamped_s=tt_c,
                )
            )
            self._oob_column_clamped_total += n_oob_md

        i = np.clip(np.searchsorted(z, md_c, side="right") - 1, 0, len(z) - 2)
        j = min(max(int(np.searchsorted(ts, tt_c, side="right")) - 1, 0), len(ts) - 2)

        fz = (md_c - z[i]) / (z[i + 1] - z[i])
        ft = (tt_c - ts[j]) / (ts[j + 1] - ts[j])

        t00 = self.table[i, j]
        t01 = self.table[i, j + 1]
        t10 = self.table[i + 1, j]
        t11 = self.table[i + 1, j + 1]
        top = t00 + ft * (t01 - t00)
        bot = t10 + ft * (t11 - t10)
        return (top + fz * (bot - top)).reshape(arr.shape)

    # -- 构造（文件 + npz 缓存） --------------------------------------------
    @classmethod
    def from_files(
        cls,
        xlsx_path: Path = DEFAULT_T_IN_XLSX,
        depth_csv_path: Path = DEFAULT_DEPTH_CSV,
        *,
        cache_path: Optional[Path] = None,
        use_cache: bool = True,
        expected_shape: Tuple[int, int] = EXPECTED_SHAPE,
    ) -> "TableTemperatureField":
        """从交付 xlsx + 井身结构 CSV 加载（默认管内 T_in）。

        缓存：``cache_path`` 缺省时落在 xlsx 同目录（``*.npz``，已被 .gitignore
        的 ``*.npz`` 规则覆盖）。命中缓存则免解析 xlsx/CSV；缓存内记录源路径，
        与请求不符时视为失效、重新解析并覆盖。
        """
        xlsx_path = Path(xlsx_path)
        depth_csv_path = Path(depth_csv_path)
        if cache_path is None:
            cache_path = xlsx_path.with_suffix(".npz")
        cache_path = Path(cache_path)

        if use_cache and cache_path.exists():
            cached = cls._load_cache(cache_path)
            if cached is not None and cached[2] == str(xlsx_path):
                table, depth_m, _src, time_s = cached
                return cls(
                    table,
                    depth_m,
                    time_s,
                    expected_shape=expected_shape,
                    source=xlsx_path,
                )

        table = _read_xlsx_numeric(xlsx_path)
        depth_m = _read_depth_axis(depth_csv_path)
        field = cls(
            table,
            depth_m,
            expected_shape=expected_shape,
            source=xlsx_path,
        )
        if use_cache:
            cls._save_cache(cache_path, field, xlsx_path)
        return field

    @staticmethod
    def _cache_fingerprint(table: np.ndarray, depth_m: np.ndarray,
                           time_s: np.ndarray) -> str:
        """内容指纹：三个数组的字节流 sha256（截断 16 位十六进制）。"""
        import hashlib
        h = hashlib.sha256()
        for arr in (table, depth_m, time_s):
            a = np.ascontiguousarray(np.asarray(arr, dtype=np.float64))
            h.update(str(a.shape).encode("ascii"))
            h.update(a.tobytes())
        return h.hexdigest()[:16]

    @staticmethod
    def _load_cache(
        cache_path: Path,
    ) -> Optional[Tuple[np.ndarray, np.ndarray, str, np.ndarray]]:
        try:
            with np.load(cache_path, allow_pickle=False) as data:
                keys = set(data.files)
                if not {"table", "depth_m", "time_s", "source",
                        "cache_version", "fingerprint"} <= keys:
                    return None        # 旧版缓存（无指纹）⇒ 视为未命中
                if int(np.asarray(data["cache_version"]).ravel()[0]) != _CACHE_VERSION:
                    return None
                t = np.asarray(data["table"], dtype=float)
                z = np.asarray(data["depth_m"], dtype=float)
                ts = np.asarray(data["time_s"], dtype=float)
                if TableTemperatureField._cache_fingerprint(t, z, ts) != str(
                        np.asarray(data["fingerprint"]).ravel()[0]):
                    return None        # 指纹不符（损坏/被改写）⇒ 回退解析源文件
                return t, z, str(data["source"]), ts
        except Exception:
            # 缓存损坏/不可读 → 视为未命中，回退解析源文件
            return None

    @staticmethod
    def _save_cache(
        cache_path: Path, field: "TableTemperatureField", xlsx_path: Path
    ) -> None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            cache_path,
            table=field.table,
            depth_m=field.depth_m,
            time_s=field.time_s,
            source=np.str_(str(xlsx_path)),
            cache_version=np.int64(_CACHE_VERSION),
            fingerprint=np.str_(TableTemperatureField._cache_fingerprint(
                field.table, field.depth_m, field.time_s)),
        )


# ---------------------------------------------------------------------------
# 交付件读取辅助
# ---------------------------------------------------------------------------
def _col_letters_to_index(letters: str) -> int:
    """Excel 列名（A/B/.../GV）→ 0 基列号。"""
    col = 0
    for ch in letters:
        col = col * 26 + (ord(ch) - ord("A") + 1)
    return col - 1


def _read_xlsx_numeric(xlsx_path: Path) -> np.ndarray:
    """读取 xlsx 首个工作表为二维 float 矩阵（纯数值、无表头；缺失单元格=NaN）。

    不依赖 openpyxl：直接解析 sheet XML 的 ``<c><v>`` 数值单元格。
    """
    xlsx_path = Path(xlsx_path)
    if not xlsx_path.exists():
        raise FileNotFoundError(f"温度场表文件不存在: {xlsx_path}")
    try:
        with zipfile.ZipFile(xlsx_path) as zf:
            xml_bytes = zf.read(_SHEET_XML)
    except KeyError as exc:
        raise ValueError(
            f"温度场表缺少工作表 {_SHEET_XML}: {xlsx_path}"
        ) from exc

    root = ET.fromstring(xml_bytes)
    ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    cells = {}
    max_row = max_col = -1
    for cell in root.iter(ns + "c"):
        ref = cell.get("r")
        if not ref:
            continue
        m = _CELL_REF_RE.match(ref)
        if m is None:
            continue
        col = _col_letters_to_index(m.group(1))
        row = int(m.group(2)) - 1
        cell_type = cell.get("t")
        if cell_type not in (None, "n"):
            # 交付表为纯数值；出现共享字符串/布尔等非数值单元格须显式报错，
            # 避免把共享字符串索引静默当温度读入
            raise ValueError(
                f"温度场表存在非数值单元格 {ref}（t={cell_type}）: {xlsx_path}"
            )
        v_node = cell.find(ns + "v")
        if v_node is None or v_node.text is None:
            continue
        cells[(row, col)] = float(v_node.text)
        max_row = max(max_row, row)
        max_col = max(max_col, col)

    if max_row < 0 or max_col < 0:
        raise ValueError(f"温度场表无有效数值单元格: {xlsx_path}")

    table = np.full((max_row + 1, max_col + 1), np.nan)
    for (r, c), val in cells.items():
        table[r, c] = val
    return table


def _read_depth_axis(depth_csv_path: Path) -> np.ndarray:
    """读井身结构 CSV 第 2 列 ``depth_well_logging_m_`` 作深度轴。"""
    depth_csv_path = Path(depth_csv_path)
    if not depth_csv_path.exists():
        raise FileNotFoundError(f"井身结构 CSV 不存在: {depth_csv_path}")

    try:
        series = pd.read_csv(depth_csv_path, usecols=[DEPTH_CSV_COLUMN])[
            DEPTH_CSV_COLUMN
        ]
    except ValueError as exc:
        raise ValueError(
            f"井身结构 CSV 缺少列 {DEPTH_CSV_COLUMN}: {depth_csv_path}"
        ) from exc
    return series.to_numpy(dtype=float)


# ---------------------------------------------------------------------------
# 交付两表加载（管内 + 环空）
# ---------------------------------------------------------------------------
def load_delivered_pair(
    *,
    t_in_xlsx: Path = DEFAULT_T_IN_XLSX,
    t_out_xlsx: Path = DEFAULT_T_OUT_XLSX,
    depth_csv_path: Path = DEFAULT_DEPTH_CSV,
    cache_dir: Optional[Path] = None,
    use_cache: bool = True,
    expected_shape: Tuple[int, int] = EXPECTED_SHAPE,
) -> Tuple[TableTemperatureField, TableTemperatureField]:
    """加载交付的管内（T_in）与环空（T_out）两表。

    校验井底（最深行）两表共享同值（浮点容差 1e-9）。
    ``cache_dir`` 给定时 npz 缓存放该目录，否则落各自 xlsx 同目录。
    """

    def _cache_for(xlsx: Path) -> Optional[Path]:
        if cache_dir is None:
            return None
        return Path(cache_dir) / (xlsx.stem + ".npz")

    t_in = TableTemperatureField.from_files(
        t_in_xlsx,
        depth_csv_path,
        cache_path=_cache_for(t_in_xlsx),
        use_cache=use_cache,
        expected_shape=expected_shape,
    )
    t_out = TableTemperatureField.from_files(
        t_out_xlsx,
        depth_csv_path,
        cache_path=_cache_for(t_out_xlsx),
        use_cache=use_cache,
        expected_shape=expected_shape,
    )
    if t_in.table.shape[0] != t_out.table.shape[0]:
        raise ValueError(
            f"管内/环空两表深度节点数不一致: {t_in.table.shape} vs {t_out.table.shape}"
        )
    bottom_in = t_in.table[-1]
    bottom_out = t_out.table[-1]
    if not np.allclose(bottom_in, bottom_out, atol=1e-9, rtol=0):
        max_diff = float(np.max(np.abs(bottom_in - bottom_out)))
        raise ValueError(
            f"井底（最深行）两表应共享同值，实际最大偏差 {max_diff:.3e} °C"
        )
    return t_in, t_out


# ---------------------------------------------------------------------------
# Phase 4d：时程扩展表加载 + 批次级 stop_t 覆盖断言
# ---------------------------------------------------------------------------
class TemperatureTableCoverageError(ValueError):
    """批次级前置断言失败：要求覆盖的时刻超过温度表时间轴末端。

    语义 = 响亮报错、不静默 clamp。``TableTemperatureField.T`` 自身的
    越界 clamp + 逐条审计行为不变——本异常只在批驱动的前置门抛出。
    """


def assert_time_table_coverage(
    fields,
    stop_t_s: float,
    *,
    label: str = "",
) -> None:
    """批次级前置断言：温度表时间轴末端须覆盖 ``stop_t_s``，否则抛错。

    供批驱动在起跑前调用（如敏感性补跑的低排量档 stop_t 校核），把历史上
    "表格档 r0.6/0.8 被静默 clamp 污染"的教训固化为硬门。

    参数
    ----------
    fields : TableTemperatureField 或其可迭代（管内/环空两表一起传）
    stop_t_s : 批次需要表覆盖到的时刻（s），通常为水泥顶替停泵时刻
    label : 批次标识，仅用于报错信息
    """
    seq = [fields] if hasattr(fields, "time_s") else list(fields)
    for f in seq:
        ts = getattr(f, "time_s", None)
        if ts is None:
            raise TypeError(f"字段 {type(f).__name__} 无 time_s 轴，不适用表格覆盖断言")
        ts_end = float(np.asarray(ts)[-1])
        if stop_t_s > ts_end:
            src = getattr(f, "source", None)
            raise TemperatureTableCoverageError(
                f"[{label or 'batch'}] stop_t={stop_t_s:.1f} s 超出温度表时间轴末端 "
                f"{ts_end:.1f} s（超出 {stop_t_s - ts_end:.1f} s；表源={src}）。"
                "禁止静默 clamp：请改用扩展表（load_extended_pair_4d）或裁短批次。"
            )


def load_extended_pair_4d(
    t_in_xlsx: Path,
    t_out_xlsx: Path,
    time_axis_csv_min: Path,
    *,
    depth_csv_path: Path = DEFAULT_DEPTH_CSV,
    marker: str = EXTENDED_TABLE_MARKER,
    expected_shape: Tuple[int, int] = EXTENDED_SHAPE_4D,
) -> Tuple[TableTemperatureField, TableTemperatureField]:
    """加载呼1-004 时程扩展温度表（333x362，表末 21600 s）。

    溯源要求（F-4 登记的硬条件，缺任一即报错）：
    - 两表文件名须含 ``_ext4d`` 标识（防把交付表/中间产物误当扩展表）；
    - 须附带时间轴侧车 CSV（单位 min，非均匀网格：0..198, 198.792801, 199..360）；
    - 形状硬校验 ``expected_shape``（默认 EXTENDED_SHAPE_4D）；
    - 井底（最深行）两表同值，容差与 ``load_delivered_pair`` 一致（1e-9）。

    实现只走 ``TableTemperatureField`` 既有的显式参数路径（time_s、
    expected_shape 均为原有入参），交付表路径 ``from_files`` /
    ``load_delivered_pair`` 零改动。扩展表不做 npz 缓存（批驱动一次性加载，
    zipfile+XML 解析成本可接受）。
    """
    for p in (t_in_xlsx, t_out_xlsx):
        if marker not in Path(p).name:
            raise ValueError(
                f"扩展表文件名须含溯源标识 {marker!r}：{p}"
            )

    axis_min = np.asarray(
        pd.read_csv(Path(time_axis_csv_min), header=None).to_numpy(dtype=float)
    ).ravel()
    if axis_min.size != expected_shape[1]:
        raise ValueError(
            f"时间轴侧车长度应为 {expected_shape[1]}（=扩展表列数），实际 {axis_min.size}"
        )
    if len(axis_min) < 2:
        raise ValueError("时间轴侧车至少需 2 个节点")
    if np.isnan(axis_min).any():
        raise ValueError("时间轴侧车含 NaN")
    if not np.all(np.diff(axis_min) > 0):
        raise ValueError("时间轴侧车非严格单调递增")
    time_s = axis_min * 60.0

    depth_m = _read_depth_axis(depth_csv_path)
    fields = []
    for p in (Path(t_in_xlsx), Path(t_out_xlsx)):
        table = _read_xlsx_numeric(p)
        fields.append(
            TableTemperatureField(
                table, depth_m, time_s, expected_shape=expected_shape, source=p
            )
        )
    t_in, t_out = fields
    if not np.allclose(t_in.table[-1], t_out.table[-1], atol=1e-9, rtol=0):
        max_diff = float(np.max(np.abs(t_in.table[-1] - t_out.table[-1])))
        raise ValueError(
            f"扩展表井底（最深行）两表应共享同值，实际最大偏差 {max_diff:.3e} °C"
        )
    return t_in, t_out
