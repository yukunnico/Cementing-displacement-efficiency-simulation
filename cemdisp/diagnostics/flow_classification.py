"""
流动分类判据 (flow_classification)
=========================================

本模块实现 Zhang & Frigaard (2023, JFM 972, A38) §3.1 的流动分类诊断，
将环空顶替流分类为三种 regime：

1. **unsteady_dispersive**（非稳态弥散）：前缘不满足相似坍缩，对流不占主导；
2. **dispersive_steady**（稳态弥散）：主前沿稳态推进，但伴随显著间隙尺度弥散
   （前缘尖峰 leading front / 壁面残留 residual wall layers）；
3. **non_dispersive_steady**（稳态非弥散）：接近活塞流的理想顶替。

物理方法
--------
用相似变换判断对流主导性：若流动由对流主导，截面平均浓度 c̄(ẑ, t̂) 在不同时刻
绘于相似变量 ẑ/t̂ 上应坍缩为一条主曲线 c̄(w_f)，其中归一化前缘速度

    w_f(c̄) = (ẑ/t̂) / ŵ₀

判据（论文式号）：

- 式 3.1：Δw_f = w_f(c̄=0.3) − w_f(c̄=0.7) ≤ 0.1 → steady front，否则 unsteady；
- 相对速度 w_r(c̄) = w_f(c̄) − 1，分解为正部 w_r+（下游弥散，快于平均泵速）
  与负部 w_r−（上游残留，如壁面层缓慢清除）；
- 式 3.2：σ_{w_r+} = sqrt(mean(w_r+²))（对 0<c̄<1 内全部采样点取均方根）；
- 式 3.4：|w̄_{r+}| = ∫₀¹ w_r+ dc̄（前缘下游总弥散面积）；
- 式 3.5：|w̄_{r−}| = −∫₀¹ w_r− dc̄（上游残留流体面积，随偏心度增大）；
- 式 3.6a,b：σ_{w_r+} > 0.08 且 |w̄_{r+}| > 0.05 同时成立 → dispersive。

输入数据约定
------------
截面平均浓度场从二维水泥快照重建（间隙宽度 b 加权方位平均）：

    c̄(s, t) = Σ_y b(y, s)·c(y, s, t) / Σ_y b(y, s)

**ŵ₀ 约定**：``AnnulusSimulationResult`` 未直接保存平均泵速，本模块默认用
c̄=0.5 等值线首末位置的前缘平均速度近似 ŵ₀（对活塞流与自相似前缘，该估计
等于真实平均泵速）；亦可通过 ``w0_m_s`` 显式传入（如由排量/环空截面积
估算的真实泵速）。

与论文的两点实现偏差（披露）：

1. 论文在**晚期时刻**评估 w_f(c̄) 收敛主曲线；本模块对所有有效快照的
   ξ=ẑ/t̂ 取算术平均作为主曲线估计，以抑制单时刻噪声。对相似坍缩良好的
   流动两者等价；对强非稳态流动，本约定给出时间平均意义下的分类。
2. 式 3.4/3.5 的积分在采样浓度网格 [0.02, 0.98] 上以梯形法则计算，
   端点值外延至 c̄=0 与 c̄=1，近似 ∫₀¹。

数值安全
--------
- 快照数 < 3（去除 t=0 初始快照后不足 2 个有效时刻）时返回 NaN 守卫结果，
  ``flow_class="insufficient_data"``；
- 浓度等值线在某时刻不存在（剖面整体位于阈值同侧）时该时刻跳过；
- 全部除法均有零除保护；间隙宽度权重 b 截断至 ≥1e-12。

使用示例
--------
>>> from cemdisp.diagnostics.flow_classification import compute_flow_classification
>>> fc = compute_flow_classification(result)
>>> print(fc.flow_class)   # "non_dispersive_steady"
>>> print(fc.delta_w_f)    # 0.03
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Dict, Optional

import numpy as np

if TYPE_CHECKING:  # 仅类型注解，避免运行时依赖 models2d 求解逻辑
    from cemdisp.models2d.annulus_d2dga import AnnulusSimulationResult


# --- 判据阈值常数（Zhang & Frigaard 2023） -----------------------------------
STEADY_FRONT_THRESHOLD = 0.1  # 式 3.1：Δw_f ≤ 0.1 → steady front
SIGMA_WR_PLUS_THRESHOLD = 0.08  # 式 3.6a：σ_{w_r+} > 0.08
ABS_WR_PLUS_THRESHOLD = 0.05  # 式 3.6b：|w̄_{r+}| > 0.05
MIN_SNAPSHOTS = 3  # 去除 t=0 初始快照后至少保留 2 个有效时刻

FLOW_CLASS_UNSTEADY_DISPERSIVE = "unsteady_dispersive"
FLOW_CLASS_DISPERSIVE_STEADY = "dispersive_steady"
FLOW_CLASS_NON_DISPERSIVE_STEADY = "non_dispersive_steady"
FLOW_CLASS_INSUFFICIENT_DATA = "insufficient_data"


@dataclass(frozen=True)
class FlowClassificationResult:
    """流动分类判据结果（Zhang & Frigaard 2023 式 3.1-3.6）。

    Attributes:
        delta_w_f: 前缘速度展布 Δw_f = w_f(c̄=0.3) − w_f(c̄=0.7)（式 3.1）
        sigma_wr_plus: 正相对速度均方根 σ_{w_r+}（式 3.2）
        abs_wr_plus: 前缘下游弥散面积 |w̄_{r+}| = ∫₀¹ w_r+ dc̄（式 3.4）
        abs_wr_minus: 上游残留面积 |w̄_{r−}| = −∫₀¹ w_r− dc̄（式 3.5）
        is_steady: 是否稳态前缘（式 3.1，Δw_f ≤ 0.1）
        is_dispersive: 是否弥散（式 3.6a,b 同时满足）
        flow_class: 三态分类 "unsteady_dispersive" / "dispersive_steady" /
            "non_dispersive_steady"；数据不足时为 "insufficient_data"
        w0_m_s: 所用归一化参考速度 ŵ₀ [m/s]（约定见模块 docstring）
    """

    delta_w_f: float
    sigma_wr_plus: float
    abs_wr_plus: float
    abs_wr_minus: float
    is_steady: bool
    is_dispersive: bool
    flow_class: str
    w0_m_s: float

    def to_dict(self) -> Dict[str, object]:
        """导出为可 JSON 序列化的字典。"""
        return {
            "delta_w_f": self.delta_w_f,
            "sigma_wr_plus": self.sigma_wr_plus,
            "abs_wr_plus": self.abs_wr_plus,
            "abs_wr_minus": self.abs_wr_minus,
            "is_steady": self.is_steady,
            "is_dispersive": self.is_dispersive,
            "flow_class": self.flow_class,
            "w0_m_s": self.w0_m_s,
        }


def _nan_result(w0_m_s: float = float("nan")) -> FlowClassificationResult:
    """数据不足时的 NaN 守卫结果（flow_class="insufficient_data"）。"""
    return FlowClassificationResult(
        delta_w_f=float("nan"),
        sigma_wr_plus=float("nan"),
        abs_wr_plus=float("nan"),
        abs_wr_minus=float("nan"),
        is_steady=False,
        is_dispersive=False,
        flow_class=FLOW_CLASS_INSUFFICIENT_DATA,
        w0_m_s=w0_m_s,
    )


def _contour_position(profile: np.ndarray, s: np.ndarray, level: float) -> float:
    """在剖面 profile(s) 上求浓度等值线 profile = level 的位置（线性插值）。

    水泥浆自入口（s 小端）向下游顶替，profile 总体随 s 递减（1→0）。
    对存在多个交叉的非单调剖面（如弥散尖峰），取**最下游**交叉点，
    对应前缘包络（leading front），与论文"所有波速 w_f 的并集"定义一致。

    Args:
        profile: 截面平均浓度剖面 c̄(s)，形状 (nz,)
        s: 轴向坐标 [m]，形状 (nz,)，单调递增
        level: 浓度阈值，取值 (0, 1)

    Returns:
        等值线位置 [m]；若剖面整体位于 level 同侧（无交叉）则返回 np.nan
    """
    above = profile >= level
    transitions = np.flatnonzero(above[:-1] != above[1:])
    if transitions.size == 0:
        return float("nan")
    j = int(transitions[-1])  # 最下游交叉点
    p0, p1 = float(profile[j]), float(profile[j + 1])
    s0, s1 = float(s[j]), float(s[j + 1])
    dp = p1 - p0
    if abs(dp) < 1.0e-12:  # 零除保护：平台期取左端点
        return s0
    return s0 + (level - p0) * (s1 - s0) / dp


def _build_c_bar_st(result: "AnnulusSimulationResult") -> np.ndarray:
    """从水泥快照重建截面平均浓度场 c̄(s,t)，间隙宽度 b 加权（spec §2.0）。

    Args:
        result: 环空二维求解结果，消费 cement_snapshots 与 geom["b"]

    Returns:
        (nz, n_snapshots) 数组，第 i 列为第 i 个快照时刻的 c̄(s)
    """
    b = np.asarray(result.geom["b"], dtype=float)
    columns = []
    for snap in result.cement_snapshots:
        snap_arr = np.asarray(snap, dtype=float)
        b_w = b
        if b_w.shape != snap_arr.shape:
            # 兼容 b 以 (ny,) 一维形式给出的情况
            b_w = np.broadcast_to(b_w.reshape(-1, 1), snap_arr.shape)
        b_w = np.clip(b_w, 1.0e-12, None)  # 零除保护
        columns.append(np.average(snap_arr, axis=0, weights=b_w))
    return np.stack(columns, axis=1)


def _estimate_w0(times: np.ndarray, z_half: np.ndarray) -> float:
    """由 c̄=0.5 等值线首末位置估计前缘平均速度 ŵ₀ [m/s]。

    约定：ŵ₀ = (z_0.5(t_last) − z_0.5(t_first)) / (t_last − t_first)。
    对活塞流与自相似前缘，该估计等于真实平均泵速。

    Args:
        times: 有效快照时刻 [s]（已过滤 t≤0），形状 (nt,)
        z_half: c̄=0.5 等值线位置 [m]，形状 (nt,)，可含 NaN

    Returns:
        ŵ₀ [m/s]；有效点不足或时间跨度为零时返回 np.nan
    """
    valid = ~np.isnan(z_half)
    if int(valid.sum()) < 2:
        return float("nan")
    t_v = times[valid]
    z_v = z_half[valid]
    dt = float(t_v[-1] - t_v[0])
    if dt <= 0.0:
        return float("nan")
    return float((z_v[-1] - z_v[0]) / dt)


def _trapz(y: np.ndarray, x: np.ndarray) -> float:
    """梯形法则积分（兼容 numpy<2.0，无 np.trapezoid 依赖）。"""
    return float(np.sum(0.5 * (y[:-1] + y[1:]) * np.diff(x)))


def compute_flow_classification(
    result: "AnnulusSimulationResult",
    *,
    w0_m_s: Optional[float] = None,
    n_levels: int = 49,
) -> FlowClassificationResult:
    """计算流动分类判据（主入口，Zhang & Frigaard 2023 式 3.1-3.6）。

    Args:
        result: 环空二维求解结果（消费 cement_snapshots / snapshot_times_s /
            geom["s"] / geom["b"]，不调用任何求解逻辑）
        w0_m_s: 可选，显式指定归一化参考速度 ŵ₀ [m/s]；
            None 时用 c̄=0.5 等值线前缘平均速度近似（约定见模块 docstring）
        n_levels: 浓度采样点数，默认 49（c̄ ∈ [0.02, 0.98]，避开 0/1 端点噪声）

    Returns:
        FlowClassificationResult；快照不足、等值线无法重建或 ŵ₀ 非法时
        返回 NaN 守卫结果（flow_class="insufficient_data"）

    Raises:
        ValueError: 快照数与时间戳数不一致
    """
    snapshots = tuple(result.cement_snapshots)
    times = tuple(result.snapshot_times_s)
    if len(snapshots) != len(times):
        raise ValueError(
            f"cement_snapshots 数量 ({len(snapshots)}) 与 "
            f"snapshot_times_s 数量 ({len(times)}) 不一致"
        )
    # NaN 守卫：快照数过少时相似变量统计无意义
    # （求解器通常含 t=0 初始快照，过滤后需至少 2 个有效时刻）
    if len(snapshots) < MIN_SNAPSHOTS:
        return _nan_result()

    s = np.asarray(result.geom["s"], dtype=float)
    t_all = np.asarray(times, dtype=float)
    c_bar_st = _build_c_bar_st(result)  # (nz, nt)

    # 相似变量 ξ = s/t 在 t=0 无定义，过滤初始快照
    valid_t = t_all > 0.0
    if int(valid_t.sum()) < 2:
        return _nan_result()
    t = t_all[valid_t]
    c_bar = c_bar_st[:, valid_t]

    levels = np.linspace(0.02, 0.98, int(n_levels))

    # 各时刻、各浓度阈值的等值线位置 z(c̄, t)（找不到则为 NaN，统计时跳过）
    z_ct = np.full((t.size, levels.size), np.nan)
    for i in range(t.size):
        profile = c_bar[:, i]
        for k, lev in enumerate(levels):
            z_ct[i, k] = _contour_position(profile, s, float(lev))

    # 归一化参考速度 ŵ₀
    if w0_m_s is not None:
        w0 = float(w0_m_s)
    else:
        idx_half = int(np.argmin(np.abs(levels - 0.5)))
        w0 = _estimate_w0(t, z_ct[:, idx_half])
    if not np.isfinite(w0) or w0 <= 1.0e-12:
        return _nan_result()

    # 主曲线 w_f(c̄)：各时刻 ξ = z/t 的算术平均，再除以 ŵ₀
    xi = z_ct / t[:, None]
    valid_xi = ~np.isnan(xi)
    counts = valid_xi.sum(axis=0)
    sums = np.where(valid_xi, xi, 0.0).sum(axis=0)
    w_f = np.where(counts > 0, sums / np.maximum(counts, 1), np.nan) / w0

    valid_lev = ~np.isnan(w_f)
    # 需要足够采样点且覆盖主前沿区间 [0.3, 0.7]（式 3.1 所需）
    if (
        int(valid_lev.sum()) < 5
        or float(levels[valid_lev].min()) > 0.3
        or float(levels[valid_lev].max()) < 0.7
    ):
        return _nan_result(w0_m_s=w0)

    lev_v = levels[valid_lev]
    w_f_v = w_f[valid_lev]

    # 式 3.1：Δw_f = w_f(c̄=0.3) − w_f(c̄=0.7)
    delta_w_f = float(np.interp(0.3, lev_v, w_f_v) - np.interp(0.7, lev_v, w_f_v))

    # 相对速度 w_r = w_f − 1 及其正/负部
    w_r = w_f_v - 1.0
    w_r_plus = np.clip(w_r, 0.0, None)
    w_r_minus = np.clip(w_r, None, 0.0)

    # 式 3.2：σ_{w_r+} = sqrt(mean(w_r+²))（对全部采样浓度点取平均，含零点）
    sigma_wr_plus = float(np.sqrt(np.mean(w_r_plus**2)))

    # 式 3.4/3.5：∫₀¹ w_r± dc̄，梯形法则 + 端点值外延至 c̄=0 与 c̄=1
    lev_ext = np.concatenate(([0.0], lev_v, [1.0]))
    w_rp_ext = np.concatenate(([w_r_plus[0]], w_r_plus, [w_r_plus[-1]]))
    w_rm_ext = np.concatenate(([w_r_minus[0]], w_r_minus, [w_r_minus[-1]]))
    abs_wr_plus = _trapz(w_rp_ext, lev_ext)
    abs_wr_minus = -_trapz(w_rm_ext, lev_ext)

    is_steady = bool(delta_w_f <= STEADY_FRONT_THRESHOLD)  # 式 3.1
    is_dispersive = bool(
        sigma_wr_plus > SIGMA_WR_PLUS_THRESHOLD  # 式 3.6a
        and abs_wr_plus > ABS_WR_PLUS_THRESHOLD  # 式 3.6b
    )
    # 三态分类：非稳态必伴随显著弥散（论文 §3.1.2），故 unsteady 直接归一类
    if not is_steady:
        flow_class = FLOW_CLASS_UNSTEADY_DISPERSIVE
    elif is_dispersive:
        flow_class = FLOW_CLASS_DISPERSIVE_STEADY
    else:
        flow_class = FLOW_CLASS_NON_DISPERSIVE_STEADY

    return FlowClassificationResult(
        delta_w_f=delta_w_f,
        sigma_wr_plus=sigma_wr_plus,
        abs_wr_plus=abs_wr_plus,
        abs_wr_minus=abs_wr_minus,
        is_steady=is_steady,
        is_dispersive=is_dispersive,
        flow_class=flow_class,
        w0_m_s=w0,
    )
