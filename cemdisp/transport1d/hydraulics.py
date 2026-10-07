# -*- coding: utf-8 -*-
"""Phase 3.1 P-2 水力学内核 —— MATLAB 甲方脚本 ``p_jaifang1.m`` 的原样复现器（HT1-004 口径）。

定位（spec §0 / §1.1；C-14 同型声明）
--------------------------------------
- **探针 / 甲方口径复现器**：不接线任何生产 runner；与 D2DGA 生产链（models2d 闭式 Ĩ₁、
  transport1d/casing_flow 幂律口径）**并列、不共享常数**。本文件是甲方 Bingham + mPa·s +
  经验常数（6895/216/5.33355/0.158278/0.81619/16800）口径，混用即污染对靶。
- 验收 = 对靶逐点复现 ``results/_probe_matlab靶_2026-10-07/sandbox/out_*.csv``
  （2026-10-07 沙箱 -batch 运行，wrapper_p30.m 捕获，原脚本零改动）。

原样复现红线（spec §1.2 —— 修 bug = 逐点失配 = 验收失败）
----------------------------------------------------------
- MATLAB :1153 管内界面 2/3 屈服混合误用 ``tau1``（应为 tau2）——**照抄**，见
  :func:`assign_casing_properties` 的 ``[原样复现 :1153]`` 注释（数值上 tau1==tau2==9.8，
  字面保留以便审查）。
- MATLAB :1464 ``safety_margin=0.003`` 只进 fprintf 不进走廊——**照抄**，见
  :func:`control_windows`（声明即终点，全链无消费方）。
- MATLAB :733 循环内反复赋 ``volume_injected_casing_1_list(1)=1200`` ⇒ 管内首界面 +1200 L
  偏置——**照抄**，见 :func:`track_casing_interfaces` 的 ``[原样复现 :733]`` 注释。
- MATLAB :132 ``n_time = floor(189.059) = 189``（尾丢 ≈41 L）——**照抄**，见
  :func:`build_pump_schedule`。
- 摩擦不动点迭代：seed 100 Pa、相对容差 1e-3、30 步强制退出后仍取末值——**照抄**，见
  :func:`_friction_bh`。
- MATLAB ``PR_Friction``：``E=0 ⇒ PR≡1`` 数值空转链——**照抄**（保留函数本体与接口位），
  :func:`friction_pr`；``Ff*PR`` 乘法保留。
- 单位混用点（mm/m、MPa/Pa、g/cm³↔kg/m³、0.00981 与 9.81 两族因子、L/min/dm²/m·s、
  dt=分钟）按移植地图 §4 单位总表逐处对齐，注释挂 MATLAB 行号。

索引约定（MATLAB 1 基镜像）
---------------------------
几何/时序内部数组统一做成 **1 基 padded**（index 0 空置，见 :func:`_pad1`），使核心循环
可逐行对照 ``p_jaifang1.m`` 审查；对外的 (189×333)/(189×5) 矩阵在返回前显式转 0 基。

MATLAB 来源（只读）：p_jaifang1.m(2057 行) / Friction_annulus_bh.m / Friction_casing_bh.m /
PR_Friction.m / 呼1-004井身结构.csv（UTF-8 BOM；**``annulus_radius_array_cm_`` 名义 radius
列实际按井眼直径使用**，移植地图 S1/R7）。

浮点次序纪律：MATLAB ``sum/cumsum/mean`` 按**逐项左结合**复刻（:func:`_sum_seq` 等），求和
顺序即双舍入顺序，禁换成 numpy pairwise 归约；表达式求值顺序逐字对齐（左到右）。
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

__all__ = [
    "MATLAB_EPS",
    "G_MPS2",
    "FACTOR_00981",
    "SENTINEL_UNFILLED",
    "SENTINEL_NOT_ENTERED",
    "SENTINEL_PUSHED_OUT",
    "SENTINEL_EXITED_SHOE",
    "DEPTH_THRESHOLDS_M",
    "StructureTable",
    "load_structure",
    "Segments",
    "build_segments",
    "FluidTable",
    "build_fluid_table",
    "PumpSchedule",
    "build_pump_schedule",
    "friction_pr",
    "friction_annulus_bh",
    "friction_casing_bh",
    "annulus_entry_volumes",
    "track_annulus_interfaces",
    "assign_annulus_properties",
    "track_casing_interfaces",
    "assign_casing_properties",
    "hydrostatic_and_friction_chain",
    "casing_static_chain",
    "pump_pressure_backtrack",
    "compute_ecd_esd",
    "control_windows",
    "four_point_backpressure",
    "bottom_window_backpressure",
    "compare_with_reference",
    "run_ht1004_target",
    # P-4 泛化附加（Phase 3.2，纯加法入口）
    "fluid_table_from",
    "build_pump_schedule_custom",
    "build_segments_custom",
    "track_casing_interfaces_custom",
    "assign_casing_properties_custom",
    "run_well_forward",
]

# --------------------------------------------------------------------------
# 常数（全部有 MATLAB 行号出处；禁止与生产链常量互换）
# --------------------------------------------------------------------------
G_MPS2: float = 9.81  # MATLAB :701/:1297/:1431 —— 本脚本重力加速度（≠生产链 9.80665）
FACTOR_00981: float = 0.00981  # MATLAB :1482 起 —— (g/cm³)×m→MPa 一步折算 = 9.81/1000
MATLAB_EPS: float = 2.220446049250313e-16  # MATLAB :1423 `eps`（=eps(1)）

# 界面哨兵（MATLAB :381-387 环空 / :794-801 管内；方向对偶，移植地图 R5）
SENTINEL_UNFILLED: float = 10000.0  # 环空：未进环空；管内：已出管鞋
SENTINEL_PUSHED_OUT: float = -1.0  # 环空：已顶出；管内：未注入
SENTINEL_NOT_ENTERED: float = 10000.0
SENTINEL_EXITED_SHOE: float = 10000.0

# S2 七段阈值（MATLAB :12-18 硬编码）
DEPTH_THRESHOLDS_M: Tuple[float, float, float, float, float, float] = (
    4025.73, 5243.21, 5578.00, 7378.05, 7521.00, 7660.00,
)
BORE_DIAMETER_MM = 245.37  # MATLAB :37/:41/:45 段1-3 电测井径常数
PIPE_OD_MM: Tuple[float, ...] = (149.2, 127.0, 168.3, 168.3, 139.7, 139.7, 0.0)
PIPE_WALL_MM: Tuple[float, ...] = (9.65, 9.65, 15.88, 15.88, 15.88, 15.88, 0.0)  # :39-59
PIANXIN: bool = True  # MATLAB :67 `pianxin = true`

# S16 关键点目标深度（MATLAB :1467-1469 最近网格点）
KEY_DEPTH_M: Tuple[float, float, float] = (5578.0, 6600.0, 7498.0)

_STRUCTURE_COLUMNS: Tuple[str, ...] = (
    "VarName1",
    "depth_well_logging_m_",
    "vertical_depth_for_logging_m_",
    "length_segment_array_m_",
    "annulus_radius_array_cm_",
    "volume_annulus_L_",
    "deg_for_logging_degree_",
    "square_annulus_dm2_",
)

# --------------------------------------------------------------------------
# 浮点次序工具（MATLAB sum/cumsum/mean 逐项左结合镜像）
# --------------------------------------------------------------------------


def _sum_seq(vals: Sequence[float]) -> float:
    """MATLAB ``sum(vector)``：左到右逐项累加（求和顺序=舍入顺序）。"""
    acc = 0.0
    for v in vals:
        acc += float(v)
    return acc


def _cumsum_seq(vals: Sequence[float]) -> List[float]:
    """MATLAB ``cumsum``：前缀逐项累加。"""
    out: List[float] = []
    acc = 0.0
    for v in vals:
        acc += float(v)
        out.append(acc)
    return out


def _mean_seq(vals: Sequence[float]) -> float:
    """MATLAB ``mean``：sum/n（n=len）。"""
    return _sum_seq(vals) / float(len(vals))


def _pad1(vals: Sequence[float]) -> np.ndarray:
    """1 基 padded 数组：a[0]=0.0 空置，a[i]=第 i 个元素（MATLAB 索引语义）。"""
    a = np.zeros(len(vals) + 1, dtype=float)
    a[1:] = np.asarray(vals, dtype=float)
    return a


def _div(a: float, b: float) -> float:
    """IEEE 语义除法（MATLAB x/0 → Inf/NaN；Python float 会抛 ZeroDivisionError）。"""
    if b == 0.0:
        if a == 0.0:
            return float("nan")
        return math.copysign(float("inf"), a) * (1.0 if math.copysign(1.0, b) > 0 else -1.0)
    return a / b


# --------------------------------------------------------------------------
# S1：井身结构加载
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class StructureTable:
    """``呼1-004井身结构.csv`` 原表（S1 :4-7）。列语义陷阱：hole_dia_cm 名为 radius 实为直径。"""

    n_segment: int
    c_depth_m: np.ndarray  # depth_well_logging_m_（测深 MD）
    tvd_csv_m: np.ndarray  # vertical_depth_for_logging_m_（CSV 电测垂深，≠cumsum 几何 TVD）
    seg_len_m: np.ndarray  # length_segment_array_m_
    hole_dia_cm: np.ndarray  # annulus_radius_array_cm_（**直径** [cm]）
    vol_annulus_L: np.ndarray  # volume_annulus_L_（甲方表值环空段容 [L]）
    deg: np.ndarray  # deg_for_logging_degree_ [°]
    sq_annulus_dm2: np.ndarray  # square_annulus_dm2_ [dm²]


def load_structure(csv_path: Path) -> StructureTable:
    """读井身结构 CSV（MATLAB :4 `readtable`）。UTF-8 BOM；列名白名单逐字符核对（地图 §5 R7）。"""
    path = Path(csv_path)
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f)
        header = tuple(next(reader))
        if header != _STRUCTURE_COLUMNS:
            raise ValueError(f"井身结构 CSV 列名与 MATLAB readtable 口径不符：{header!r}")
        rows = [[float(x) for x in r] for r in reader if r and any(x.strip() for x in r)]
    if not rows:
        raise ValueError("井身结构 CSV 无数据行")
    arr = np.asarray(rows, dtype=float)  # float() = 正确舍入 strtod，与 readtable 同值
    return StructureTable(
        n_segment=arr.shape[0],
        c_depth_m=arr[:, 1],
        tvd_csv_m=arr[:, 2],
        seg_len_m=arr[:, 3],
        hole_dia_cm=arr[:, 4],
        vol_annulus_L=arr[:, 5],
        deg=arr[:, 6],
        sq_annulus_dm2=arr[:, 7],
    )


# --------------------------------------------------------------------------
# S2/S3/S6/S10-PR：几何、截面、两条独立容积链、偏心 PR
# --------------------------------------------------------------------------


@dataclass
class Segments:
    """几何派生量。容积链纪律（地图 R13）：vab（CSV 链，界面定位）≠ area×len（重算链）≠ vtop（管内链）。"""

    n: int
    c_depth_m1b: np.ndarray
    seg_len_m1b: np.ndarray  # MATLAB structure_data.length_segment_array_m_
    sq_annulus_dm2_1b: np.ndarray
    d_bit_m1b: np.ndarray  # diameter_bit_out [m]（:37-63）
    d_cas_out_m1b: np.ndarray  # diameter_casing_out [m]
    d_cas_in_m1b: np.ndarray  # diameter_casing_in [m]
    out_diam_bole_mm1b: np.ndarray  # 井眼直径 mm（:66 `*10`，cm→mm）
    area_cout_m2_1b: np.ndarray  # :64
    area_cin_m2_1b: np.ndarray  # :65
    d_i_m1b: np.ndarray  # :171 = casing_out
    d_o_m1b: np.ndarray  # :172 = 名义 radius 列 /100（**直径** m）
    A_annulus_m2_1b: np.ndarray  # :664 sq/100（摩擦函数传入但未消费，原样）
    cos_deg_1b: np.ndarray  # :158
    vert_len_m1b: np.ndarray  # :159 abs(len*cos)
    tvd_cum_m1b: np.ndarray  # :161 cumsum（S15 :1423 会就地 eps 替换——由 compute_ecd_esd 处理）
    vab_L1b: np.ndarray  # :264 volume_all_from_bottom（CSV 链）
    vtop_L1b: np.ndarray  # :728 volume_all_from_top_casing（管内链）
    capacity_pipe_L: float  # :70/:75-77 volume_in_drilling_casing_L
    seg_masks: List[np.ndarray]  # S2 idx1..idx7（0 基布尔）
    pr_seg: List[float]  # PR_segment1..7（:582-643；E=0 ⇒ 非空段恒 1.0，空段 0）
    pr_1b: np.ndarray  # :645-662 按阈值填的 PR(i)


def build_segments(structure: StructureTable) -> Segments:
    """S2 七段阈值 + S3 截面/容积 + S6 几何 + S10 的 PR 七段。逐字对齐 :9-81/:155-173/:263-266/:581-664。"""
    n = structure.n_segment
    cd = structure.c_depth_m
    # ---- S2 七段掩码（:27-33；idx7 恒空——c_depth 最大 7660，地图 R16，分支保留）
    t1, t2, t3, t4, t5, t6 = DEPTH_THRESHOLDS_M
    m1 = cd <= t1
    m2 = (cd > t1) & (cd <= t2)
    m3 = (cd > t2) & (cd <= t3)
    m4 = (cd > t3) & (cd <= t4)
    m5 = (cd > t4) & (cd <= t5)
    m6 = (cd > t5) & (cd <= t6)
    m7 = cd > t6
    masks = [m1, m2, m3, m4, m5, m6, m7]

    d_bit = np.zeros(n)
    d_out = np.zeros(n)
    d_in = np.zeros(n)
    for k, m in enumerate(masks[:3]):
        d_bit[m] = BORE_DIAMETER_MM * 0.001  # :37/:41/:45 mm→m
    for k, m in enumerate(masks[3:], start=3):
        d_bit[m] = structure.hole_dia_cm[m] * 0.01  # :49/:53/:57/:61 "cm"→m（直径口径）
    for k, m in enumerate(masks):
        od = PIPE_OD_MM[k]
        wt = PIPE_WALL_MM[k]
        d_out[m] = od * 0.001  # :38…
        d_in[m] = (od - wt * 2) * 0.001  # :39… (OD-2*wall)
        # 段7（:62-63）od=wt=0 ⇒ d_out=d_in=0（无内管）——恒等空段
    # ---- S3 截面/分段容积（:64-71；乘序逐字：pi*(a²-b²)/4）
    area_cout = math.pi * (d_bit**2 - d_out**2) / 4.0
    area_cin = math.pi * (d_in**2) / 4.0
    out_diam_bole = structure.hole_dia_cm * 10.0  # :66 cm→mm
    vol_cas_seg = area_cin * structure.seg_len_m  # :70
    vol_ann_seg = area_cout * structure.seg_len_m  # :71
    # ---- :74 volume_in_fenggu / :80 volume_in_allfluid / :265 length_all_from_bottom：
    #      死变量（计算后全程未消费，地图 S3/S6），不复刻数值路径、不影响任何输出。
    capacity_m3 = _sum_seq(vol_cas_seg)  # :75 sum()
    capacity_L = capacity_m3 * 1000.0  # :77 m³→L
    # ---- S6 几何（:157-173）
    cos_deg = np.array([math.cos(math.radians(float(dg))) for dg in structure.deg])
    vert_len = np.abs(structure.seg_len_m * cos_deg)  # :159 abs(len*cos)（先乘后 abs）
    tvd_cum = np.asarray(_cumsum_seq(vert_len))  # :161 cumsum
    d_i = d_out.copy()  # :171
    d_o = structure.hole_dia_cm / 100.0  # :172（名义 radius 列 → m，直径语义）
    # ---- 两条独立容积链（地图 R13 禁互换）
    vab = np.zeros(n)
    for i in range(n):  # :264 每 i 重做 sum(CSV_vol(i:end))——左结合镜像
        vab[i] = _sum_seq(structure.vol_annulus_L[i:])
    vtop = np.zeros(n)
    prod_cin = area_cin * structure.seg_len_m  # :728 area_cin(1:i).*len(1:i)
    for i in range(n):
        vtop[i] = _sum_seq(prod_cin[: i + 1]) * 1000.0
    # ---- S10 PR 七段（:581-662；E=0 ⇒ 非空段 PR≡1.0 —— 原样复现红线，空段 PR=0 分支保留）
    pr_seg: List[float] = []
    for m in masks:
        if np.any(m):
            mean_dole = _mean_seq(out_diam_bole[m])  # mm
            mean_cas = _mean_seq(d_out[m]) * 1000.0  # m→mm
            pr_seg.append(friction_pr(mean_dole, mean_cas))
        else:
            pr_seg.append(0.0)  # :587-588 空段 ⇒ PR_segmentK = 0
    pr = np.zeros(n)
    for i in range(n):  # :646-661 阈值 if/elseif 链
        c = float(cd[i])
        if c <= t1:
            pr[i] = pr_seg[0]
        elif c <= t2:
            pr[i] = pr_seg[1]
        elif c <= t3:
            pr[i] = pr_seg[2]
        elif c <= t4:
            pr[i] = pr_seg[3]
        elif c <= t5:
            pr[i] = pr_seg[4]
        elif c <= t6:
            pr[i] = pr_seg[5]
        else:
            pr[i] = pr_seg[6]
    a_annulus = structure.sq_annulus_dm2 / 100.0  # :664 dm²→m²（摩擦函数内未消费，原样传入）

    return Segments(
        n=n,
        c_depth_m1b=_pad1(cd),
        seg_len_m1b=_pad1(structure.seg_len_m),
        sq_annulus_dm2_1b=_pad1(structure.sq_annulus_dm2),
        d_bit_m1b=_pad1(d_bit),
        d_cas_out_m1b=_pad1(d_out),
        d_cas_in_m1b=_pad1(d_in),
        out_diam_bole_mm1b=_pad1(out_diam_bole),
        area_cout_m2_1b=_pad1(area_cout),
        area_cin_m2_1b=_pad1(area_cin),
        d_i_m1b=_pad1(d_i),
        d_o_m1b=_pad1(d_o),
        A_annulus_m2_1b=_pad1(a_annulus),
        cos_deg_1b=_pad1(cos_deg),
        vert_len_m1b=_pad1(vert_len),
        tvd_cum_m1b=_pad1(tvd_cum),
        vab_L1b=_pad1(vab),
        vtop_L1b=_pad1(vtop),
        capacity_pipe_L=capacity_L,
        seg_masks=masks,
        pr_seg=pr_seg,
        pr_1b=_pad1(pr),
    )


# --------------------------------------------------------------------------
# S4/S5：14 流体表 + MPD 回压数组
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class FluidTable:
    """S4 :87-133 / S5 :138-152。顺序 = [1,2,3,4,5,6,7,8,9,91,92,93,94,95]（14 路）+ 背景 0。"""

    dt_min: float  # :131 dt=1（**分钟**）
    rates_Lmin: Tuple[float, ...]  # :88-102 pump_rateN = x*1000
    vols_L: Tuple[float, ...]  # :115-128
    rou0: float
    miu0: float
    tau0: float
    rou_gcc: Tuple[float, ...]
    miu_mPas: Tuple[float, ...]
    tau_Pa: Tuple[float, ...]
    bp_MPa: Tuple[float, ...]  # S5（当前全 0 = 正演基准）


def build_fluid_table() -> FluidTable:
    """逐字常量表（地图 S4）。浮点面值按 IEEE 求值：如 1.4*1000、0.9*1000 与 MATLAB 同值。"""
    rates = (
        1.4 * 1000, 1.2 * 1000, 1.2 * 1000, 1.2 * 1000, 1.25 * 1000, 1.2 * 1000,
        1.5 * 1000, 1.4 * 1000, 1.4 * 1000, 1.2 * 1000, 1.0 * 1000, 0.9 * 1000,
        0.8 * 1000, 0.7 * 1000,
    )
    vols = (
        25 * 1000, 16 * 1000, 10 * 1000, 48 * 1000, 28 * 1000, 2 * 1000,
        29 * 1000, 14 * 1000, 1 * 1000, 14 * 1000, 10 * 1000, 10 * 1000,
        10 * 1000, 7.1 * 1000,
    )
    rou = (1.75, 1.95, 1.75, 1.93, 1.90, 1.70, 1.90, 1.90, 1.02, 1.90, 1.90, 1.90, 1.90, 1.90)
    miu = (58.0, 58.0, 65.0, 200.0, 180.0, 50.0, 50.0, 50.0, 50.0, 55.0, 55.0, 55.0, 55.0, 55.0)
    tau = (9.8, 9.8, 10.0, 14.0, 14.0, 9.0, 9.5, 9.2, 9.0, 9.5, 9.5, 9.5, 9.5, 9.5)
    return FluidTable(
        dt_min=1.0,
        rates_Lmin=rates,
        vols_L=vols,
        rou0=1.90,
        miu0=53.0,
        tau0=8.5,
        rou_gcc=rou,
        miu_mPas=miu,
        tau_Pa=tau,
        bp_MPa=(0.0,) * 14,  # S5 :138-151 全 0
    )


# --------------------------------------------------------------------------
# S7：泵时表
# --------------------------------------------------------------------------


@dataclass
class PumpSchedule:
    """S7 :268-341。n_time=floor(189.059)=189（**原样复现 :132**：尾丢 ≈41 L，停泵分支不触发但保留）。"""

    n_time: int
    dt_min: float
    nodes_min: List[float]  # pump_time_node(1..14)，非整数
    pump_Lmin_1b: np.ndarray  # Pump_values_time_list（1 基，index0 空置）
    bp_MPa_1b: np.ndarray
    vol_all_L_1b: np.ndarray  # volume_injected_all_list
    q_m3s_1b: np.ndarray  # :341 (Pump/60)/1000


def build_pump_schedule(fluids: FluidTable) -> PumpSchedule:
    """泵率/回压按**时间节点**展开（非体积节点，地图 R10）。边界比较符链序逐字镜像 :291-332。"""
    rates = fluids.rates_Lmin
    vols = fluids.vols_L
    # ---- :132 n_time = floor(Σ v_k/PV(k) / dt)——左结合
    total_min = 0.0
    for v, r in zip(vols, rates):
        total_min += v / r
    n_time = int(math.floor(total_min / fluids.dt_min))  # = 189（尾丢 ≈0.059 min ≈ 41 L）
    # ---- :274-288 pump_time_node(1..14)：每行 = 前缀左结合和
    nodes: List[float] = []
    acc = 0.0
    for v, r in zip(vols, rates):
        acc += v / r
        nodes.append(acc / fluids.dt_min)
    # ---- 1 基展开
    pump = np.zeros(n_time + 1)
    bp = np.zeros(n_time + 1)
    vol_all = np.zeros(n_time + 1)
    q = np.zeros(n_time + 1)
    # :269-271 初始时刻（t=1）
    pump[1] = rates[0]
    bp[1] = fluids.bp_MPa[0]
    vol_all[1] = pump[1] * fluids.dt_min
    for t in range(2, n_time + 1):
        tv = float(t)
        # :291-335 分支链（第1支 `<`；第2/3支 `>=/<=`；第4支起 `>/<=`；先命中先赢——禁合并）
        if tv < nodes[0]:
            pump[t] = rates[0]; bp[t] = fluids.bp_MPa[0]
        elif tv >= nodes[0] and tv <= nodes[1]:
            pump[t] = rates[1]; bp[t] = fluids.bp_MPa[1]
        elif tv >= nodes[1] and tv <= nodes[2]:
            pump[t] = rates[2]; bp[t] = fluids.bp_MPa[2]
        elif tv > nodes[2] and tv <= nodes[3]:
            pump[t] = rates[3]; bp[t] = fluids.bp_MPa[3]
        elif tv > nodes[3] and tv <= nodes[4]:
            pump[t] = rates[4]; bp[t] = fluids.bp_MPa[4]
        elif tv > nodes[4] and tv <= nodes[5]:
            pump[t] = rates[5]; bp[t] = fluids.bp_MPa[5]
        elif tv > nodes[5] and tv <= nodes[6]:
            pump[t] = rates[6]; bp[t] = fluids.bp_MPa[6]
        elif tv > nodes[6] and tv <= nodes[7]:
            pump[t] = rates[7]; bp[t] = fluids.bp_MPa[7]
        elif tv > nodes[7] and tv <= nodes[8]:
            pump[t] = rates[8]; bp[t] = fluids.bp_MPa[8]
        elif tv > nodes[8] and tv <= nodes[9]:
            pump[t] = rates[9]; bp[t] = fluids.bp_MPa[9]
        elif tv > nodes[9] and tv <= nodes[10]:
            pump[t] = rates[10]; bp[t] = fluids.bp_MPa[10]
        elif tv > nodes[10] and tv <= nodes[11]:
            pump[t] = rates[11]; bp[t] = fluids.bp_MPa[11]
        elif tv > nodes[11] and tv <= nodes[12]:
            pump[t] = rates[12]; bp[t] = fluids.bp_MPa[12]
        elif tv > nodes[12] and tv <= nodes[13]:
            pump[t] = rates[13]; bp[t] = fluids.bp_MPa[13]
        elif tv > nodes[13]:  # :333-335 停泵分支（189 步内不触发，代码路径保留）
            pump[t] = 0.0; bp[t] = 0.0
        vol_all[t] = vol_all[t - 1] + pump[t] * fluids.dt_min  # :338
    # :341 L/min → m³/s：((Q/60)/1000) 逐元素
    q[1:] = np.asarray([((float(p) / 60.0) / 1000.0) for p in pump[1:]])
    return PumpSchedule(
        n_time=n_time, dt_min=fluids.dt_min, nodes_min=nodes,
        pump_Lmin_1b=pump, bp_MPa_1b=bp, vol_all_L_1b=vol_all, q_m3s_1b=q,
    )


# --------------------------------------------------------------------------
# 摩擦函数（移植地图 §3；逐字式，禁替换为闭式 Bingham/标准式）
# --------------------------------------------------------------------------


def friction_pr(Dw_mm: float, Dc_mm: float, E: float = 0.0, fai: float = 12.0) -> float:
    """PR_Friction.m 移植（汪海阁偏心环空压耗比）。**原样复现**：``E=0`` 硬编码 ⇒ PR≡1.0 空转
    （地图 R9）；保留函数本体与 E/fai 接口位供后续放开。输入 mm，比值无量纲。"""
    sigma = _div(Dc_mm, Dw_mm)
    term1 = ((0.1045 * fai + 0.073) * E) * (sigma**0.8454)
    term2 = ((0.8841 * (fai ** (-0.28879))) * (E**2)) * (sigma**0.1852)
    term3 = ((0.5658 * (fai ** (-0.28879))) * (E**3)) * (sigma**0.2527)
    return 1.0 - term1 - term2 + term3


def _friction_bh(rho_a: float, V_a: float, mu_p_a: float, tau_y: float,
                 D_w: float, D_do: float, QL: float, *, laminar_mult: float) -> Tuple[float, int]:
    """Friction_annulus_bh.m / Friction_casing_bh.m 共有体（两函数除层流乘子外逐字相同）。

    单位：rho_a [kg/m³]、V_a [m/s]、mu_p_a [**mPa·s**]（函数首行 /1000 转 Pa·s——头注释
    误标 Pa·s，地图 R7 以代码为准）、tau_y [Pa]、D_w/D_do [m]、QL [m³/s]；返回 (Pa/m, 流型)。

    **原样复现**：井壁切应力不动点迭代 seed=100 Pa、`abs(err)>1e-3 && iter<30` 强制退出、
    未收敛仍取末值 tau_w_new（地图 R12；禁闭式替换）。
    """
    mu = mu_p_a / 1000.0  # mPa·s → Pa·s（Friction_*.m:23）
    gap = D_w - D_do
    half_gap = gap / 2.0
    # tau_w_new = 8*mu*QL / (pi*r_h^2*(D_w+D_do)/2) + 3/2*tau_y - 1/2*tau_y^3/tau_w_old^2
    num = (8.0 * mu) * QL
    den = ((math.pi * (half_gap**2)) * (D_w + D_do)) / 2.0
    c1 = _div(num, den)
    c2 = (3.0 / 2.0) * tau_y
    c3 = (1.0 / 2.0) * (tau_y**3)
    tau_w_old = 100.0  # 种子（Friction_*.m:26）
    err_tau = 1.0
    it = 0
    tau_w_new = float("nan")
    while abs(err_tau) > 1e-3 and it < 30:  # 容差 1e-3 / maxIter 30（原样）
        tau_w_new = (c1 + c2) - _div(c3, tau_w_old**2)
        err_tau = _div(abs(tau_w_new - tau_w_old), tau_w_old)  # 除旧值（迭代路径敏感，原样）
        tau_w_old = tau_w_new
        it += 1
    tau_w = tau_w_new  # 未收敛也用最后一次值（原样）
    int1 = _div(tau_y, tau_w)
    He = _div(16800.0 * int1, (1.0 - int1) ** 3)
    re_c = _div((1.0 - (4.0 / 3.0) * int1) + (1.0 / 3.0) * (int1**4), 8.0 * int1) * He
    re_a = _div((((0.81619 * rho_a) * V_a) * gap), mu)
    if re_a <= re_c:  # 层流（注释"×0.75"名实不符，实码 0.8/1.25——地图 R8，以代码常数为准）
        t1 = _div((6895.0 * mu) * V_a, 216.0 * (gap**2))
        t2 = _div(5.33355 * tau_y, gap)
        return laminar_mult * (t1 + t2), 1
    # 紊流（两函数相同）：((((0.158278*rho^0.75)*V^1.75)*mu^0.25)/gap^1.25) 左到右乘序原样
    ff = 0.158278 * (rho_a**0.75)
    ff = ff * (V_a**1.75)
    ff = ff * (mu**0.25)
    ff = _div(ff, gap**1.25)
    return ff, 2


def friction_annulus_bh(rho_a: float, V_a: float, A_a: float, mu_p_a: float, tau_y: float,
                        D_w: float, D_do: float, QL: float) -> Tuple[float, int]:
    """Friction_annulus_bh.m：层流乘子 0.8（A_a 传入但函数内**从未消费**，原样保留签名）。"""
    del A_a  # MATLAB 形参占位（Friction_annulus_bh.m 内无引用——地图 §3.1）
    return _friction_bh(rho_a, V_a, mu_p_a, tau_y, D_w, D_do, QL, laminar_mult=0.8)


def friction_casing_bh(rho_a: float, V_a: float, A_a: float, mu_p_a: float, tau_y: float,
                       D_w: float, D_do: float, QL: float) -> Tuple[float, int]:
    """Friction_casing_bh.m：层流乘子 1.25；主脚本传 D_w=管内径、D_do=0 ⇒ 间隙=全径（:1308/:229 口径）。"""
    del A_a
    return _friction_bh(rho_a, V_a, mu_p_a, tau_y, D_w, D_do, QL, laminar_mult=1.25)


# --------------------------------------------------------------------------
# S8：环空进入体积（阈值链 = 管内总容 + 前序流体体积）
# --------------------------------------------------------------------------


def annulus_entry_volumes(sched: PumpSchedule, segs: Segments, fluids: FluidTable) -> np.ndarray:
    """S8 :343-360：返回 5×(nt+1) 1 基数组（流体 1..5 进环空后的逐时累计注入 [L]）。

    阈值链左结合求和：cap; cap+v1; cap+v1+v2; cap+v1+v2+v3; cap+v1+v2+v3+v4。
    t=1 时条件必假（管内远未充满），MATLAB `list(t-1)` 的 0 索引路径不会触发。
    """
    nt = sched.n_time
    out = np.zeros((5, nt + 1))
    cap = segs.capacity_pipe_L
    thr = [
        cap,
        cap + fluids.vols_L[0],
        cap + fluids.vols_L[0] + fluids.vols_L[1],
        cap + fluids.vols_L[0] + fluids.vols_L[1] + fluids.vols_L[2],
        cap + fluids.vols_L[0] + fluids.vols_L[1] + fluids.vols_L[2] + fluids.vols_L[3],
    ]
    dt = sched.dt_min
    for k in range(5):
        for t in range(1, nt + 1):
            if sched.vol_all_L_1b[t] > thr[k]:
                out[k, t] = out[k, t - 1] + sched.pump_Lmin_1b[t] * dt
    return out


# --------------------------------------------------------------------------
# S9：环空 5 界面追踪 + 属性分配 + 流速
# --------------------------------------------------------------------------


def _annulus_interface_scan(v: float, vab: np.ndarray, sq: np.ndarray, deg: np.ndarray,
                            n: int) -> Tuple[float, float, float, float]:
    """S9 :366-389 模板的**忠实仿真**：`for i=2:n` 重复扫描，未命中分支的 i 保留先前值
    （等价"查找唯一 i"，但保留 v==vab(1) 等病态值时列表停留 0 的原语义——地图 R5）。
    返回 (depth_tag, residual_volume, residual_height, vertical_residual_height)。
    """
    tag = 0.0
    rv = 0.0
    rh = 0.0
    vrh = 0.0
    for i in range(2, n + 1):
        if v >= vab[i] and v < vab[i - 1]:
            tag = float(i - 1)  # :367
            rv = v - vab[i]
            rh = (rv / sq[i - 1]) / 10.0  # L/dm²→dm→m（/10）
            vrh = rh * math.cos(math.radians(float(deg[i - 1])))  # cosd（:372）
        elif v > 0.0 and v < vab[n]:
            tag = float(n)
            rv = v
            rh = (rv / sq[n]) / 10.0
            vrh = rh * math.cos(math.radians(float(deg[n])))
        elif v > vab[1]:  # 充满整个环空（顶出）
            tag = SENTINEL_PUSHED_OUT
            rv = SENTINEL_PUSHED_OUT
            rh = SENTINEL_PUSHED_OUT
            vrh = SENTINEL_PUSHED_OUT
        elif v <= 0.0:  # 尚未进入环空
            tag = SENTINEL_UNFILLED
            rv = SENTINEL_UNFILLED
            rh = SENTINEL_UNFILLED
            vrh = SENTINEL_UNFILLED
    return tag, rv, rh, vrh


def track_annulus_interfaces(vol_inj: np.ndarray, segs: Segments,
                             deg_1b: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """对 5 界面逐时定位（S9 :363-494）。返回 tags (nt,5)、residual_height (nt,5)（0 基输出）。

    ``deg_1b`` = deg_for_logging_degree_ 的 1 基表（cosd 用，:372/:379）。
    """
    nt = vol_inj.shape[1] - 1
    tags = np.zeros((nt, 5))
    res_h = np.zeros((nt, 5))
    for t in range(1, nt + 1):
        for k in range(5):
            tag, _rv, rh, _vrh = _annulus_interface_scan(
                float(vol_inj[k, t]), segs.vab_L1b, segs.sq_annulus_dm2_1b, deg_1b, segs.n)
            tags[t - 1, k] = tag
            res_h[t - 1, k] = rh
    return tags, res_h


def assign_annulus_properties(tags: np.ndarray, res_h: np.ndarray, sched: PumpSchedule,
                              segs: Segments, fluids: FluidTable) -> Dict[str, np.ndarray]:
    """S9 :497-577：逐时逐段 elseif 链属性分配（VOF 线性混合 + clamp [0,1]）+ 流速。

    链序逐字（:499-558）：`i<depth_tag_0_1` 起、`i==tag` 混合、`tag_k<i<tag_{k+1}` 分区、
    `i>tag_5` 收尾。10000/-1 哨兵参与比较当 ±∞ 用（地图 R5）。
    流速 :576：`Pump/sq(i)/10/60`（L/min→(L/min)/dm²=dm/min→/10 m/min→/60 m/s）。
    """
    nt, n = sched.n_time, segs.n
    rou0, miu0, tau0 = fluids.rou0, fluids.miu0, fluids.tau0
    r1, r2, r3, r4, r5 = fluids.rou_gcc[:5]
    m1, m2, m3, m4, m5 = fluids.miu_mPas[:5]
    t1, t2, t3, t4, t5 = fluids.tau_Pa[:5]
    sq = segs.sq_annulus_dm2_1b
    plen = segs.seg_len_m1b
    pump = sched.pump_Lmin_1b
    rou = np.zeros((nt + 1, n + 1))
    miu = np.zeros((nt + 1, n + 1))
    tau = np.zeros((nt + 1, n + 1))
    velo = np.zeros((nt + 1, n + 1))
    for t in range(1, nt + 1):
        tg = tags[t - 1]  # [tag0_1, tag1_2, tag2_3, tag3_4, tag4_5]
        rh = res_h[t - 1]
        for i in range(1, n + 1):
            if i < tg[0]:
                rou[t, i], miu[t, i], tau[t, i] = rou0, miu0, tau0
            elif i == tg[0]:
                p = rh[0] / plen[i]
                p = min(max(p, 0.0), 1.0)  # :505 clamp
                q = 1.0 - p
                rou[t, i] = p * r1 + q * rou0
                miu[t, i] = p * m1 + q * miu0
                tau[t, i] = p * t1 + q * tau0
            elif tg[0] < i < tg[1]:
                rou[t, i], miu[t, i], tau[t, i] = r1, m1, t1
            elif i == tg[1]:
                p = min(max(rh[1] / plen[i], 0.0), 1.0)
                q = 1.0 - p
                rou[t, i] = p * r2 + q * r1
                miu[t, i] = p * m2 + q * m1
                tau[t, i] = p * t2 + q * t1
            elif tg[1] < i < tg[2]:
                rou[t, i], miu[t, i], tau[t, i] = r2, m2, t2
            elif i == tg[2]:
                p = min(max(rh[2] / plen[i], 0.0), 1.0)
                q = 1.0 - p
                rou[t, i] = p * r3 + q * r2
                miu[t, i] = p * m3 + q * m2
                tau[t, i] = p * t3 + q * t2
            elif tg[2] < i < tg[3]:
                rou[t, i], miu[t, i], tau[t, i] = r3, m3, t3
            elif i == tg[3]:
                p = min(max(rh[3] / plen[i], 0.0), 1.0)
                q = 1.0 - p
                rou[t, i] = p * r4 + q * r3
                miu[t, i] = p * m4 + q * m3
                tau[t, i] = p * t4 + q * t3
            elif tg[3] < i < tg[4]:
                rou[t, i], miu[t, i], tau[t, i] = r4, m4, t4
            elif i == tg[4]:
                p = min(max(rh[4] / plen[i], 0.0), 1.0)
                q = 1.0 - p
                rou[t, i] = p * r5 + q * r4
                miu[t, i] = p * m5 + q * m4
                tau[t, i] = p * t5 + q * t4
            elif tg[4] < i:  # :554 i>tag(t,5)
                rou[t, i], miu[t, i], tau[t, i] = r5, m5, t5
            # :576 流速（无条件，链外）
            velo[t, i] = ((pump[t] / sq[i]) / 10.0) / 60.0
    return {"rou_gcc": rou, "miu_mPas": miu, "tau_Pa": tau, "velo_m_s": velo}


# --------------------------------------------------------------------------
# S10：环空静液柱 + 循环压耗 + 井口边界递推
# --------------------------------------------------------------------------


def hydrostatic_and_friction_chain(segs: Segments, props: Dict[str, np.ndarray],
                                   sched: PumpSchedule) -> Dict[str, np.ndarray]:
    """S10 :666-716。井口边界 i=1 从第 1 段（30 m 网格）起算：`bp*1e6 + Ff*len + 9.81*ρ*Δz`；
    i≥2 递推。乘序逐字：`(9.81*ρ_kg)*vert`、`Ff*len`、`((bp*1e6 + fl) + hydro)`。

    PR 乘法（:697 `Ff*current_PR`）原样保留——E=0 时 PR≡1（空转，红线）。
    返回 Pa 三矩阵 + 流型（MPa 视图由调用方 /1e6，等价 :714-716）。
    """
    nt, n = sched.n_time, segs.n
    rou_kg = props["rou_gcc"] * 1000.0  # :665 g/cm³→kg/m³
    p_total = np.zeros((nt + 1, n + 1))
    p_static = np.zeros((nt + 1, n + 1))
    p_friction = np.zeros((nt + 1, n + 1))
    ff = np.zeros((nt + 1, n + 1))
    flow = np.zeros((nt + 1, n + 1))
    for t in range(1, nt + 1):
        for i in range(1, n + 1):
            ffv, flw = friction_annulus_bh(
                float(rou_kg[t, i]),
                float(props["velo_m_s"][t, i]),
                float(segs.A_annulus_m2_1b[i]),
                float(props["miu_mPas"][t, i]),
                float(props["tau_Pa"][t, i]),
                float(segs.d_o_m1b[i]),
                float(segs.d_i_m1b[i]),
                float(sched.q_m3s_1b[t]),
            )
            if PIANXIN:  # :681-698
                ffv = ffv * float(segs.pr_1b[i])  # current_PR 按 c_depth 阈值链 == pr_1b[i]
            ff[t, i], flow[t, i] = ffv, float(flw)
            fl = ffv * float(segs.seg_len_m1b[i])
            hydro = (9.81 * float(rou_kg[t, i])) * float(segs.vert_len_m1b[i])
            if i == 1:
                p_total[t, i] = (sched.bp_MPa_1b[t] * 1e6 + fl) + hydro
                p_static[t, i] = hydro
                p_friction[t, i] = fl
            else:
                p_total[t, i] = (p_total[t, i - 1] + fl) + hydro
                p_static[t, i] = p_static[t, i - 1] + hydro
                p_friction[t, i] = p_friction[t, i - 1] + fl
    return {"Ff_Pa_m": ff, "flow_pattern": flow, "P_total_Pa": p_total,
            "P_static_Pa": p_static, "P_friction_Pa": p_friction,
            "P_total_MPa": p_total / 1e6, "P_static_MPa": p_static / 1e6,
            "P_friction_MPa": p_friction / 1e6}


# --------------------------------------------------------------------------
# S11：管内 14 界面追踪（含 :733 +1200 L 偏置）+ 属性分配 + 流速
# --------------------------------------------------------------------------


def _casing_interface_scan(v: float, vtop: np.ndarray, area_cin: np.ndarray,
                           n: int) -> Tuple[float, float, float]:
    """S11 :782-803 模板（方向对偶：自顶累计容、tag=界面所在段 i）。

    哨兵互换：`10000`=已出管鞋、`-1`=未注入（地图 R5/R6）。分支链序与 0 初值保留语义同
    :func:`_annulus_interface_scan`。
    """
    tag = 0.0
    rv = 0.0
    rh = 0.0
    for i in range(2, n + 1):
        if v >= vtop[i - 1] and v < vtop[i]:
            tag = float(i)
            rv = vtop[i] - v  # 残余量 = 段容 − 已注入（自顶）
            rh = rv / (area_cin[i] * 1000.0)  # L→m³（*1000 面积化）后折 m
        elif v > 0.0 and v < vtop[1]:
            tag = 1.0
            rv = v
            rh = rv / (area_cin[1] * 1000.0)
        elif v >= vtop[n]:
            tag = SENTINEL_EXITED_SHOE
            rv = SENTINEL_EXITED_SHOE
            rh = SENTINEL_EXITED_SHOE
        elif v <= 0.0:
            tag = SENTINEL_PUSHED_OUT
            rv = SENTINEL_PUSHED_OUT
            rh = SENTINEL_PUSHED_OUT
    return tag, rv, rh


def track_casing_interfaces(sched: PumpSchedule, segs: Segments,
                            fluids: FluidTable) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """S11 :732-777 注入体积 + :779-1116 14 界面追踪。返回 (vol_cas 14×(nt+1), tags (nt,14), res_h (nt,14))。

    **[原样复现 :733]** `volume_injected_casing_1_list(1) = 1200;` 位于 `for t=2:n_time` **体内**，
    每轮覆盖下标 1 ⇒ 管内第一界面全序列带 +1200 L 隐式偏置（其余 13 路没有；地图 R4）。
    此偏置前移井浆→先导浆界面时刻/深度、影响泵压早期形状——禁止当笔误删除。
    """
    nt = sched.n_time
    dt = sched.dt_min
    vol_c = np.zeros((14, nt + 1))
    # 阈值链（:738-774）：casing_k(k≥2) 起算当 Σvol[0..k-2] < 总注入；casing_1 当总注入>0
    thr = [0.0]
    acc = 0.0
    for k in range(1, 14):
        acc += fluids.vols_L[k - 1]  # v1; v1+v2; …（左结合）
        thr.append(acc)
    for t in range(2, nt + 1):
        vol_c[0, 1] = 1200.0  # ← [原样复现 :733] 循环内反复赋 1200 L 偏置
        if sched.vol_all_L_1b[t] > thr[0]:
            vol_c[0, t] = vol_c[0, t - 1] + sched.pump_Lmin_1b[t] * dt
        for k in range(1, 14):
            if sched.vol_all_L_1b[t] > thr[k]:
                vol_c[k, t] = vol_c[k, t - 1] + sched.pump_Lmin_1b[t] * dt
    tags = np.zeros((nt, 14))
    res_h = np.zeros((nt, 14))
    for t in range(1, nt + 1):
        for k in range(14):
            tag, _rv, rh = _casing_interface_scan(
                float(vol_c[k, t]), segs.vtop_L1b, segs.area_cin_m2_1b, segs.n)
            tags[t - 1, k] = tag
            res_h[t - 1, k] = rh
    return vol_c, tags, res_h


def assign_casing_properties(tags: np.ndarray, res_h: np.ndarray, sched: PumpSchedule,
                             segs: Segments, fluids: FluidTable, *,
                             legacy_tau_bug: bool = True) -> Dict[str, np.ndarray]:
    """S11 :1119-1286：管内 14 界面 elseif 链（方向对偶：`i>tag_0_1→井浆` 起 … `i<tag_14→rou95` 终）
    + 流速（`((Pump/1000)/60)/area_cin`，Pump=0 置 0）。

    **[原样复现 :1153]** 界面 2/3 的 τ 混合原文为
    `tau = proportion*tau1 + (1-proportion)*tau3`（正确应为 tau2/tau3；同行 rou/miu 与环空对偶
    :531 均正确）。`legacy_tau_bug=True`（默认=复现态）照抄；tau1==tau2==9.8 故数值等价，但
    字面保留、禁"顺手修正"（spec §1.2/地图 R2；是否修由用户另裁）。
    """
    nt, n = sched.n_time, segs.n
    f = list(range(14))  # 位置 0..13 ↔ 流体 1,2,3,4,5,6,7,8,9,91,92,93,94,95
    R = list(fluids.rou_gcc)
    M = list(fluids.miu_mPas)
    T = list(fluids.tau_Pa)
    R0, M0, T0 = fluids.rou0, fluids.miu0, fluids.tau0
    plen = segs.seg_len_m1b
    acin = segs.area_cin_m2_1b
    pump = sched.pump_Lmin_1b
    rou = np.zeros((nt + 1, n + 1))
    miu = np.zeros((nt + 1, n + 1))
    tau = np.zeros((nt + 1, n + 1))
    velo = np.zeros((nt + 1, n + 1))

    def mix(t_: int, i: int, k: int, lo: int, hi: int):
        """界面 k 混合：p=residual/len clamp [0,1]；rou/miu/tau = p*lo + q*hi（管内 p=下方流体份额）。"""
        p = min(max(res_h[t_ - 1, k] / plen[i], 0.0), 1.0)
        q = 1.0 - p
        rou[t_, i] = p * R[lo] + q * R[hi]
        miu[t_, i] = p * M[lo] + q * M[hi]
        tau[t_, i] = p * T[lo] + q * T[hi]

    # 位置索引：流体 k 在列表位 pos：f1=0…f5=4, f6=5, f7=6, f8=7, f9=8, f91=9, f92=10, f93=11, f94=12, f95=13
    pos1, pos2, pos3, pos4, pos5, pos6, pos7, pos8, pos9, pos91, pos92, pos93, pos94, pos95 = range(14)
    for t in range(1, nt + 1):
        tc = tags[t - 1]
        for i in range(1, n + 1):
            if i > tc[0]:
                rou[t, i], miu[t, i], tau[t, i] = R0, M0, T0  # :1121
            elif i == tc[0]:
                p = min(max(res_h[t - 1, 0] / plen[i], 0.0), 1.0)
                q = 1.0 - p
                rou[t, i] = p * R0 + q * R[pos1]  # :1129（p=井浆份额）
                miu[t, i] = p * M0 + q * M[pos1]
                tau[t, i] = p * T0 + q * T[pos1]
            elif tc[1] < i < tc[0]:
                rou[t, i], miu[t, i], tau[t, i] = R[pos1], M[pos1], T[pos1]  # :1132-1135
            elif i == tc[1]:
                mix(t, i, 1, pos1, pos2)  # :1136-1142
            elif tc[2] < i < tc[1]:
                rou[t, i], miu[t, i], tau[t, i] = R[pos2], M[pos2], T[pos2]
            elif i == tc[2]:
                # :1147-1153：rou/miu 正确（2/3），**tau 原文混用 tau1/tau3**
                p = min(max(res_h[t - 1, 2] / plen[i], 0.0), 1.0)
                q = 1.0 - p
                rou[t, i] = p * R[pos2] + q * R[pos3]
                miu[t, i] = p * M[pos2] + q * M[pos3]
                if legacy_tau_bug:
                    # [原样复现 :1153] `proportion*tau1 + proportion_0*tau3`
                    tau[t, i] = p * T[pos1] + q * T[pos3]
                else:  # 修正档（未裁定启用；默认走上面复现态）
                    tau[t, i] = p * T[pos2] + q * T[pos3]
            elif tc[3] < i < tc[2]:
                rou[t, i], miu[t, i], tau[t, i] = R[pos3], M[pos3], T[pos3]
            elif i == tc[3]:
                mix(t, i, 3, pos3, pos4)
            elif tc[4] < i < tc[3]:
                rou[t, i], miu[t, i], tau[t, i] = R[pos4], M[pos4], T[pos4]
            elif i == tc[4]:
                mix(t, i, 4, pos4, pos5)
            elif tc[5] < i < tc[4]:
                rou[t, i], miu[t, i], tau[t, i] = R[pos5], M[pos5], T[pos5]
            elif i == tc[5]:
                mix(t, i, 5, pos5, pos6)
            elif tc[6] < i < tc[5]:
                rou[t, i], miu[t, i], tau[t, i] = R[pos6], M[pos6], T[pos6]
            elif i == tc[6]:
                mix(t, i, 6, pos6, pos7)
            elif tc[7] < i < tc[6]:
                rou[t, i], miu[t, i], tau[t, i] = R[pos7], M[pos7], T[pos7]
            elif i == tc[7]:
                mix(t, i, 7, pos7, pos8)
            elif tc[8] < i < tc[7]:
                rou[t, i], miu[t, i], tau[t, i] = R[pos8], M[pos8], T[pos8]
            elif i == tc[8]:
                mix(t, i, 8, pos8, pos9)
            elif tc[9] < i < tc[8]:
                rou[t, i], miu[t, i], tau[t, i] = R[pos9], M[pos9], T[pos9]
            elif i == tc[9]:
                mix(t, i, 9, pos9, pos91)
            elif tc[10] < i < tc[9]:
                rou[t, i], miu[t, i], tau[t, i] = R[pos91], M[pos91], T[pos91]
            elif i == tc[10]:
                mix(t, i, 10, pos91, pos92)
            elif tc[11] < i < tc[10]:
                rou[t, i], miu[t, i], tau[t, i] = R[pos92], M[pos92], T[pos92]
            elif i == tc[11]:
                mix(t, i, 11, pos92, pos93)
            elif tc[12] < i < tc[11]:
                rou[t, i], miu[t, i], tau[t, i] = R[pos93], M[pos93], T[pos93]
            elif i == tc[12]:
                mix(t, i, 12, pos93, pos94)
            elif tc[13] < i < tc[12]:
                rou[t, i], miu[t, i], tau[t, i] = R[pos94], M[pos94], T[pos94]
            elif i == tc[13]:
                mix(t, i, 13, pos94, pos95)
            elif i < tc[13]:  # :1275
                rou[t, i], miu[t, i], tau[t, i] = R[pos95], M[pos95], T[pos95]
            # :1282-1286 流速
            if pump[t] > 0:
                velo[t, i] = ((pump[t] / 1000.0) / 60.0) / acin[i]
            else:
                velo[t, i] = 0.0
    return {"rou_gcc": rou, "miu_mPas": miu, "tau_Pa": tau, "velo_m_s": velo}


# --------------------------------------------------------------------------
# S12/S13：管内静液柱 + 井底向上反推 + 泵压
# --------------------------------------------------------------------------


def casing_static_chain(props: Dict[str, np.ndarray], segs: Segments,
                        nt: int) -> Dict[str, np.ndarray]:
    """S12 :1294-1303：`i==1: 9.81*ρ_kg*vert`；`else: 前值 + 9.81*ρ*vert`；MPa 视图 /1e6。"""
    n = segs.n
    rou_kg = props["rou_gcc"] * 1000.0  # :1290
    st = np.zeros((nt + 1, n + 1))
    for t in range(1, nt + 1):
        for i in range(1, n + 1):
            hydro = (9.81 * float(rou_kg[t, i])) * float(segs.vert_len_m1b[i])
            st[t, i] = hydro if i == 1 else st[t, i - 1] + hydro
    return {"P_static_Pa": st, "P_static_MPa": st / 1e6}


def pump_pressure_backtrack(segs: Segments, props: Dict[str, np.ndarray], sched: PumpSchedule,
                            p_ann_total_Pa_1b: np.ndarray, static_Pa_1b: np.ndarray,
                            *, p_bit_drop: float = 0.0) -> Dict[str, np.ndarray]:
    """S13 :1305-1349：管内摩阻 → 井底(=环空井底列 + ``P_bit_drop≡0`` :251) → 自底向上反推 → 泵压。

    原样口径（地图 R14）：
    - `pressure_casing_friction` 累计链从 i=2 起（i=1 摩阻不入链）；
    - 井口泵压式 :1333 单独补加 i=1 摩阻 —— 两处口径不同是**原样语义**，不得统一；
    - 负值钳 0（:1336-1338 自由下落修正）；
    - 停泵分支（:1339-1345）保留代码路径（189 步内不触发）。
    """
    nt, n = sched.n_time, segs.n
    rou_kg = props["rou_gcc"] * 1000.0
    ff = np.zeros((nt + 1, n + 1))
    flow = np.zeros((nt + 1, n + 1))
    p_c = np.zeros((nt + 1, n + 1))  # pressure_casing [Pa]
    p_c_friction = np.zeros((nt + 1, n + 1))
    pump_surf = np.zeros(nt + 1)  # MPa（1 基）
    for t in range(1, nt + 1):
        for i in range(1, n + 1):
            ffv, flw = friction_casing_bh(
                float(rou_kg[t, i]),
                float(props["velo_m_s"][t, i]),
                float(segs.area_cin_m2_1b[i]),
                float(props["miu_mPas"][t, i]),
                float(props["tau_Pa"][t, i]),
                float(segs.d_cas_in_m1b[i]),
                0.0,  # :1315 D_do=0 ⇒ 间隙=全径
                float(sched.q_m3s_1b[t]),
            )
            ff[t, i], flow[t, i] = ffv, float(flw)
        if sched.pump_Lmin_1b[t] > 0:  # :1321
            bottom = float(p_ann_total_Pa_1b[t, n]) + p_bit_drop  # :1322（P_bit_drop≡0，R15）
            p_c[t, n] = bottom
            for i in range(n - 1, 0, -1):  # :1325-1327
                p_c[t, i] = (p_c[t, i + 1] + ff[t, i + 1] * float(segs.seg_len_m1b[i + 1])) \
                    - (static_Pa_1b[t, i + 1] - static_Pa_1b[t, i])
            for i in range(2, n + 1):  # :1329-1331（i=1 不入链——原样）
                p_c_friction[t, i] = p_c_friction[t, i - 1] + ff[t, i] * float(segs.seg_len_m1b[i])
            true_surf = (p_c[t, 1] + ff[t, 1] * float(segs.seg_len_m1b[1])) - static_Pa_1b[t, 1]
            pv = true_surf / 1e6  # :1334
            if pv < 0.0:  # :1336-1338 负值钳 0（原样）
                pv = 0.0
            pump_surf[t] = pv
        else:  # :1339-1345 停泵分支（不触发；保留）
            p_c[t, n] = float(p_ann_total_Pa_1b[t, n]) + p_bit_drop
            for i in range(n - 1, 0, -1):
                p_c[t, i] = p_c[t, i + 1] - (9.81 * float(rou_kg[t, i + 1])) * float(segs.vert_len_m1b[i + 1])
            pump_surf[t] = p_c[t, 1] / 1e6
    return {"Ff_Pa_m": ff, "flow_pattern": flow, "P_casing_Pa": p_c,
            "P_casing_MPa": p_c / 1e6,  # :1349
            "P_casing_friction_Pa": p_c_friction,
            "P_casing_friction_MPa": p_c_friction / 1e6,
            "pump_surface_MPa": pump_surf}


# --------------------------------------------------------------------------
# S15：ECD/ESD
# --------------------------------------------------------------------------


def compute_ecd_esd(p_ann_Pa: np.ndarray, p_ann_static_Pa: np.ndarray,
                    p_cas_Pa: np.ndarray, tvd_cum_m1b: np.ndarray, nt: int, n: int) -> Dict[str, np.ndarray]:
    """S15 :1422-1457。**两步除法逐字**：`P/(9.81*TVD)` 再 `/1000`（与一步折算差 ~1 ulp，
    对靶要求保序）。`TVD_cum(TVD_cum==0)=eps`（:1423，eps=eps(1)）就地镜像；`TVD>0` 分支守卫。
    """
    tvd = tvd_cum_m1b.copy()
    for i in range(1, n + 1):  # :1423
        if tvd[i] == 0.0:
            tvd[i] = MATLAB_EPS
    ecd_ann_kg = np.zeros((nt + 1, n + 1))
    ecd_ann_g = np.zeros((nt + 1, n + 1))
    esd_g = np.zeros((nt + 1, n + 1))
    ecd_c_kg = np.zeros((nt + 1, n + 1))
    ecd_c_g = np.zeros((nt + 1, n + 1))
    for t in range(1, nt + 1):
        for i in range(1, n + 1):
            if tvd[i] > 0:  # :1430/:1448
                ecd_ann_kg[t, i] = _div(float(p_ann_Pa[t, i]), 9.81 * float(tvd[i]))
                ecd_ann_g[t, i] = ecd_ann_kg[t, i] / 1000.0  # 两步（:1432）
                esd_g[t, i] = _div(float(p_ann_static_Pa[t, i]), 9.81 * float(tvd[i])) / 1000.0  # :1433 两步
                ecd_c_kg[t, i] = _div(float(p_cas_Pa[t, i]), 9.81 * float(tvd[i]))
                ecd_c_g[t, i] = ecd_c_kg[t, i] / 1000.0  # :1450 两步
            else:  # 0 分支保留（tvd 已 eps 化，不达）
                ecd_ann_kg[t, i] = 0.0
                ecd_ann_g[t, i] = 0.0
                ecd_c_kg[t, i] = 0.0
                ecd_c_g[t, i] = 0.0
    return {"ECD_annulus_kg_m3": ecd_ann_kg, "ECD_annulus_g_cm3": ecd_ann_g,
            "ESD_annulus_g_cm3": esd_g, "ECD_casing_kg_m3": ecd_c_kg,
            "ECD_casing_g_cm3": ecd_c_g, "TVD_cum_guarded_m1b": tvd}


# --------------------------------------------------------------------------
# S16：控压窗与关键点索引
# --------------------------------------------------------------------------


def control_windows(structure: StructureTable, segs: Segments, tvd_cum_m1b: np.ndarray,
                    nt: int) -> Dict[str, object]:
    """S16 :1459-1498。窗值 = 呼1-004 现场设计口径（代码为准：:1460 上限 1.975，
    注释块 :1682 的 1.970 是过期注释——地图 S21 注记）。

    **[原样复现 :1464]** `safety_margin = 0.003`：MATLAB 中只进 fprintf（:1753）不进走廊公式
    （:1716-1727）⇒ 本函数**声明即终点、无任何消费方**；照文补码会内缩走廊、S21/S22 全变。

    关键点 = 最近网格点（`min(abs(c_depth - d))`，并列取小索引，:1467-1469）；
    CSV 缺关键点 ⇒ error()（:1471-1473）——空网格时抛错镜像。
    """
    if nt == 0 or segs.n == 0:
        raise RuntimeError("井身结构 CSV 缺少 5578m、6600m 或 7498m 关注点（MATLAB :1471-1473）")
    idx: List[int] = []
    for d in KEY_DEPTH_M:
        diffs = np.abs(structure.c_depth_m - d)
        idx.append(int(np.argmin(diffs)) + 1)  # 1 基；argmin 并列取小索引 == MATLAB min 首指标
    tvd_bottom = float(structure.tvd_csv_m[-1])  # :1479 **CSV 电测垂深**（≠TVD_cum(end)，原样）
    out: Dict[str, object] = {
        "safe_ECD_5568_lower": 1.940, "safe_ECD_5568_upper": 1.975,  # :1460
        "safe_ECD_6880_lower": 1.940, "safe_ECD_6880_upper": 1.975,  # :1461
        "safe_ECD_7463_lower": 1.940, "safe_ECD_7463_upper": 1.975,  # :1462
        "safe_ECD_bottom_lower": 2.010, "safe_ECD_bottom_upper": 2.050,  # :1463
        "safety_margin": 0.003,  # ← [原样复现 :1464] 装饰品：声明即终点，全链无消费方
        "idx_5568": idx[0], "idx_6880": idx[1], "idx_7463": idx[2],  # 1 基
        "idx_bottom_critical": segs.n,  # :1470
        "TVD_bottom_m": tvd_bottom,  # :1479
        "time_minutes_real_1b": _pad1([t * segs_dt_min(nt) for t in range(1, nt + 1)]),
    }
    return out


def segs_dt_min(nt: int) -> float:
    """time_minutes_real = (1:n_time)'*dt，dt=1（:1480）。独立小函数避免闭包依赖。"""
    return 1.0


# --------------------------------------------------------------------------
# S21/S22：四点走廊 + 控压后 ECD
# --------------------------------------------------------------------------


def four_point_backpressure(p_static_MPa: np.ndarray, p_friction_MPa: np.ndarray,
                            win: Dict[str, object], tvd_cum_m1b: np.ndarray,
                            segs: Segments, nt: int, tol_MPa: float = 1e-6) -> Dict[str, np.ndarray]:
    """S21 :1692-1741 + S22 :1764-1774（捕获列）。

    原样要点：下界 max 列表顺序 [6880, 7463, 5568, bottom, 0]（:1720）；上界 min 列表顺序
    [7463, 6880, 5568, bottom]（:1727）；冲突判 `BP_lower > BP_upper + 1e-6` ⇒ conflict=1、
    取**上界**（防漏优先 :1734）；否则取中点。井底 TVD 用 `TVD_cum(end)`（几何链，非 CSV 列）。
    `safety_margin` 不参与（红线 :1464）。
    """
    i5 = int(win["idx_5568"]); i6 = int(win["idx_6880"]); i7 = int(win["idx_7463"])
    n = segs.n
    tvd5, tvd6, tvd7, tvdB = (float(tvd_cum_m1b[i5]), float(tvd_cum_m1b[i6]),
                              float(tvd_cum_m1b[i7]), float(tvd_cum_m1b[n]))
    P_safe = {
        "5568_l": win["safe_ECD_5568_lower"] * FACTOR_00981 * tvd5,  # :1699-1706
        "5568_u": win["safe_ECD_5568_upper"] * FACTOR_00981 * tvd5,
        "6880_l": win["safe_ECD_6880_lower"] * FACTOR_00981 * tvd6,
        "6880_u": win["safe_ECD_6880_upper"] * FACTOR_00981 * tvd6,
        "7463_l": win["safe_ECD_7463_lower"] * FACTOR_00981 * tvd7,
        "7463_u": win["safe_ECD_7463_upper"] * FACTOR_00981 * tvd7,
        "bot_l": win["safe_ECD_bottom_lower"] * FACTOR_00981 * tvdB,
        "bot_u": win["safe_ECD_bottom_upper"] * FACTOR_00981 * tvdB,
    }
    required = np.zeros(nt + 1)
    new_bottom = np.zeros(nt + 1)
    bp_lower = np.zeros(nt + 1)
    bp_upper = np.zeros(nt + 1)
    conflict = np.zeros(nt + 1)
    for t in range(1, nt + 1):
        pb5 = float(p_static_MPa[t, i5]) + float(p_friction_MPa[t, i5])  # :1710（不含回压）
        pb6 = float(p_static_MPa[t, i6]) + float(p_friction_MPa[t, i6])
        pb7 = float(p_static_MPa[t, i7]) + float(p_friction_MPa[t, i7])
        pbb = float(p_static_MPa[t, n]) + float(p_friction_MPa[t, n])
        bp_l = max([P_safe["6880_l"] - pb6, P_safe["7463_l"] - pb7,
                    P_safe["5568_l"] - pb5, P_safe["bot_l"] - pbb, 0.0])  # :1720 顺序原样
        bp_u = min([P_safe["7463_u"] - pb7, P_safe["6880_u"] - pb6,
                    P_safe["5568_u"] - pb5, P_safe["bot_u"] - pbb])  # :1727 顺序原样
        bp_lower[t], bp_upper[t] = bp_l, bp_u
        if bp_l > bp_u + tol_MPa:  # :1732 冲突取上界（防漏优先）
            conflict[t] = 1.0
            required[t] = bp_u
        else:
            required[t] = (bp_l + bp_u) / 2.0
        new_bottom[t] = pbb + required[t]  # :1740
    # S22 捕获列（:1764-1774）：回压沿全环空平移
    ecd5_new = np.zeros(nt + 1)
    ecd6_new = np.zeros(nt + 1)
    ecd7_new = np.zeros(nt + 1)
    ecd_b_new = np.zeros(nt + 1)
    d5, d6, d7, dB = (FACTOR_00981 * tvd5, FACTOR_00981 * tvd6,
                      FACTOR_00981 * tvd7, FACTOR_00981 * tvdB)
    for t in range(1, nt + 1):
        ecd5_new[t] = (float(p_static_MPa[t, i5]) + float(p_friction_MPa[t, i5]) + required[t]) / d5
        ecd6_new[t] = (float(p_static_MPa[t, i6]) + float(p_friction_MPa[t, i6]) + required[t]) / d6
        ecd7_new[t] = (float(p_static_MPa[t, i7]) + float(p_friction_MPa[t, i7]) + required[t]) / d7
        ecd_b_new[t] = new_bottom[t] / dB
    return {"required_backpressure_MPa": required, "BP_lower_MPa": bp_lower,
            "BP_upper_MPa": bp_upper, "conflict_flag": conflict,
            "new_bottom_pressure_MPa": new_bottom, "ECD_5568_new": ecd5_new,
            "ECD_6880_new": ecd6_new, "ECD_7463_new": ecd7_new, "ECD_bottom_new": ecd_b_new}


# --------------------------------------------------------------------------
# S25：井底窗"持稳+冻结"状态机
# --------------------------------------------------------------------------


@dataclass
class FreezeResult:
    """S25 产出。`bp_apply_MPa` 等均为 1 基数组（index0 空置）。"""

    ecd_base_bottom: np.ndarray
    bp_ecd_win: np.ndarray
    bp_apply_MPa: np.ndarray
    p_bottom_ctrl: np.ndarray
    ecd_bottom_ctrl: np.ndarray
    p_shoe_ctrl: np.ndarray
    ecd_shoe_ctrl: np.ndarray
    p_wh_ctrl: np.ndarray
    seg_s: List[int]
    seg_e: List[int]
    overpress: np.ndarray
    bp_min_win: np.ndarray
    bp_max_win: np.ndarray


def bottom_window_backpressure(p_static_MPa: np.ndarray, p_friction_MPa: np.ndarray,
                               win: Dict[str, object], tvd_cum_m1b: np.ndarray,
                               segs: Segments, nt: int, setpoint_offset: float = 0.015,
                               merge_gap_min: float = 5.0) -> FreezeResult:
    """S25 :1952-2005。**前向状态机不可一拍向量化**（地图 R17）：跌破段逐时补量持稳至
    `ECD_setpoint_win = 2.025`；段外沿用最近一次跌破段的滚动 `freeze_ecd_win`；
    ≤5 min 游程合并（:1972-1978，按 time_minutes_real 差值比较）。

    乘序原样：`BP_apply = ((BP_ecd*0.00981)*TVD_cum(end))`（:1998）。
    分支优先级原样：overpress（自然超上限 ⇒ 0）→ below_win（补量+更新冻结）→ else（冻结沿用）。
    """
    n = segs.n
    i5 = int(win["idx_5568"])
    tvdB = float(tvd_cum_m1b[n])
    d_b = FACTOR_00981 * tvdB
    lo = float(win["safe_ECD_bottom_lower"]); hi = float(win["safe_ECD_bottom_upper"])
    p_base = np.zeros(nt + 1)
    ecd_base = np.zeros(nt + 1)
    for t in range(1, nt + 1):  # :1953-1954
        p_base[t] = float(p_static_MPa[t, n]) + float(p_friction_MPa[t, n])
        ecd_base[t] = p_base[t] / d_b
    p_safe_l = lo * FACTOR_00981 * tvdB  # :1957（乘序 ((lo*0.00981)*tvd)）
    p_safe_u = hi * FACTOR_00981 * tvdB
    setpoint = lo + setpoint_offset  # :1959 = 2.025
    bp_min = np.zeros(nt + 1)
    bp_max = np.zeros(nt + 1)
    over = np.zeros(nt + 1, dtype=bool)
    for t in range(1, nt + 1):  # :1962-1964
        bp_min[t] = max(p_safe_l - p_base[t], 0.0)
        bp_max[t] = p_safe_u - p_base[t]
        over[t] = bp_max[t] < 0.0
    # ---- 跌破段识别与合并（:1967-1983；0/1 游程 find 语义镜像）
    below_raw = [ecd_base[t] < lo for t in range(1, nt + 1)]
    seq = [0] + [1 if b else 0 for b in below_raw] + [0]
    d0 = [seq[k + 1] - seq[k] for k in range(len(seq) - 1)]  # diff，1 基游程头尾
    st0 = [k + 1 for k, x in enumerate(d0) if x == 1]  # find(d0==1)（1 基时间号）
    en0 = [k for k, x in enumerate(d0) if x == -1]  # find(d0==-1)-1
    seg_s: List[int] = []
    seg_e: List[int] = []
    for k in range(len(st0)):
        if seg_e and (float(st0[k]) - float(seg_e[-1])) <= merge_gap_min:  # :1973（time*1）
            seg_e[-1] = en0[k]  # 并入前一段（仅改末点）
        else:
            seg_s.append(st0[k])
            seg_e.append(en0[k])
    below_win = [False] * (nt + 1)
    for k in range(len(seg_s)):
        for t in range(seg_s[k], seg_e[k] + 1):
            below_win[t] = True
    # ---- 回压施加状态机（:1988-1997）
    bp_ecd = np.zeros(nt + 1)
    freeze = 0.0
    for t in range(1, nt + 1):
        if over[t]:
            bp_ecd[t] = 0.0
        elif below_win[t]:
            bp_ecd[t] = setpoint - ecd_base[t]
            freeze = bp_ecd[t]  # 段内滚动更新（交点=段末值即冻结值）
        else:
            bp_ecd[t] = freeze  # 交点后冻结沿用
    bp_apply = np.zeros(nt + 1)
    p_bot = np.zeros(nt + 1)
    ecd_bot = np.zeros(nt + 1)
    p_shoe = np.zeros(nt + 1)
    ecd_shoe = np.zeros(nt + 1)
    p_wh = np.zeros(nt + 1)
    d5 = FACTOR_00981 * float(tvd_cum_m1b[i5])
    for t in range(1, nt + 1):
        bp_apply[t] = (bp_ecd[t] * FACTOR_00981) * tvdB  # :1998 乘序原样
        p_bot[t] = p_base[t] + bp_apply[t]  # :2001-2002
        ecd_bot[t] = p_bot[t] / d_b
        p_shoe[t] = (float(p_static_MPa[t, i5]) + float(p_friction_MPa[t, i5])) + bp_apply[t]  # :2003
        ecd_shoe[t] = p_shoe[t] / d5
        p_wh[t] = (float(p_static_MPa[t, 1]) + float(p_friction_MPa[t, 1])) + bp_apply[t]  # :2005
    return FreezeResult(
        ecd_base_bottom=ecd_base, bp_ecd_win=bp_ecd, bp_apply_MPa=bp_apply,
        p_bottom_ctrl=p_bot, ecd_bottom_ctrl=ecd_bot, p_shoe_ctrl=p_shoe,
        ecd_shoe_ctrl=ecd_shoe, p_wh_ctrl=p_wh, seg_s=seg_s, seg_e=seg_e,
        overpress=over, bp_min_win=bp_min, bp_max_win=bp_max,
    )


# --------------------------------------------------------------------------
# S19/S20：对比模块（reference 插值）
# --------------------------------------------------------------------------


def _read_table_csv(path: Path) -> Dict[str, List[float]]:
    """readtable 镜像：utf-8-sig，列名→浮点值列表。"""
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        cols = {h: [] for h in header}
        for row in reader:
            if not row or not any(x.strip() for x in row):
                continue
            for h, x in zip(header, row):
                cols[h].append(float(x))
    return cols


def _unique_stable(vals: Sequence[float]) -> Tuple[List[float], List[int]]:
    """MATLAB `unique(x,'stable')`：按首次出现序去重 + 首次索引（0 基）。"""
    seen: Dict[float, int] = {}
    uniq: List[float] = []
    first: List[int] = []
    for j, v in enumerate(vals):
        if v not in seen:
            seen[v] = j
            uniq.append(v)
            first.append(j)
    return uniq, first


def _interp1_linear(x: List[float], y: List[float], q: Sequence[float]) -> List[float]:
    """MATLAB `interp1(x,y,q,'linear')` 镜像：x 严格升序（unique 后天然成立）；
    q 已在界内（钳位后）。区间内 `y0 + ((y1-y0)/(x1-x0))*(q-x0)`。"""
    out: List[float] = []
    m = len(x)
    if m == 1:
        return [float(y[0]) for _ in q]
    import bisect
    for v in q:
        j = bisect.bisect_right(x, v) - 1
        if j < 0:
            j = 0
        if j >= m - 1:
            j = m - 2
        if v == x[m - 1]:
            out.append(float(y[m - 1]))
            continue
        slope = (float(y[j + 1]) - float(y[j])) / (float(x[j + 1]) - float(x[j]))
        out.append(float(y[j]) + slope * (float(v) - float(x[j])))
    return out


def compare_with_reference(calculated_1b: np.ndarray, cum_vol_m3_1b: np.ndarray,
                           ref_csv: Path, value_col: str, nt: int) -> np.ndarray:
    """S19/S20 :1611-1678 通用体：unique-stable（取首现压力/ECD）→ 累计量钳位（禁外推 :1621）→
    interp1 linear → 误差列。返回 (nt,5)：[time, cumvol, calculated, design_interp, error]。
    缺文件 ⇒ error() 镜像（:1613-1615/:1647-1649）。
    """
    if not Path(ref_csv).is_file():
        raise RuntimeError(f"缺少 8.11 参考表：{ref_csv}")
    tbl = _read_table_csv(Path(ref_csv))
    vol_raw = tbl["volume_m3"]
    val_raw = tbl[value_col]
    uniq, first = _unique_stable(vol_raw)
    vals = [val_raw[j] for j in first]
    x_lo, x_hi = uniq[0], uniq[-1]
    clamped = [min(max(float(cum_vol_m3_1b[t]), x_lo), x_hi) for t in range(1, nt + 1)]
    interp = _interp1_linear(uniq, vals, clamped)
    out = np.zeros((nt, 5))
    for t in range(1, nt + 1):
        out[t - 1, 0] = float(t)  # time_minutes_real（dt=1）
        out[t - 1, 1] = float(cum_vol_m3_1b[t])
        out[t - 1, 2] = float(calculated_1b[t])
        out[t - 1, 3] = interp[t - 1]
        out[t - 1, 4] = float(calculated_1b[t]) - interp[t - 1]
    return out


# --------------------------------------------------------------------------
# 装配入口：复现靶全套 out_* 序列
# --------------------------------------------------------------------------


def run_ht1004_target(sandbox_dir: Path) -> Dict[str, object]:
    """顺序 = 移植地图 §2 数据流图。返回 dict（键与 out_*.csv 同名，0 基矩阵；附标量）。

    - ``out_result_matrix_volume_pressure`` (189,5)：[time, CumVol_m3, ECD_casing_g(:,end),
      pump_MPa, annuli_bottom_MPa]（S18 :1608-1609）
    - ``out_pump_pressure_comparison`` / ``out_bottom_ecd_comparison`` (189,5)（S19/S20）
    - ``out_four_point_backpressure`` (189,5)：[time, required, BP_lower, BP_upper, conflict]
    - ``out_four_point_ecd_new`` (189,5)：[time, ECD_5568/6880/7463/bottom_new]
    - ``out_pump_pressure_surface`` (189,2)、``out_annuli_bottom_pressure`` (189,2)
    - ``out_ecd_casing_full`` (189,333)
    - ``out_window_backpressure_ctrl`` (189,5)：[time, ECD_base, ECD_bottom_ctrl, ECD_shoe_ctrl, P_wh]
    - ``bp_apply_win_MPa`` (189,)、``cumulative_volume_m3`` (189,)、``summary`` 标量组
    """
    sandbox = Path(sandbox_dir)
    st = load_structure(sandbox / "呼1-004井身结构.csv")
    segs = build_segments(st)
    fluids = build_fluid_table()
    sched = build_pump_schedule(fluids)
    nt = sched.n_time

    vol_a = annulus_entry_volumes(sched, segs, fluids)
    tags_a, res_h_a = track_annulus_interfaces(vol_a, segs, _pad1(st.deg))
    props_a = assign_annulus_properties(tags_a, res_h_a, sched, segs, fluids)
    ann = hydrostatic_and_friction_chain(segs, props_a, sched)

    vol_c, tags_c, res_h_c = track_casing_interfaces(sched, segs, fluids)
    props_c = assign_casing_properties(tags_c, res_h_c, sched, segs, fluids)
    st_c = casing_static_chain(props_c, segs, nt)
    cas = pump_pressure_backtrack(segs, props_c, sched, ann["P_total_Pa"], st_c["P_static_Pa"])

    ecd = compute_ecd_esd(ann["P_total_Pa"], ann["P_static_Pa"], cas["P_casing_Pa"],
                          segs.tvd_cum_m1b, nt, segs.n)
    win = control_windows(st, segs, ecd["TVD_cum_guarded_m1b"], nt)
    corridor = four_point_backpressure(ann["P_static_MPa"], ann["P_friction_MPa"], win,
                                       ecd["TVD_cum_guarded_m1b"], segs, nt)
    fz = bottom_window_backpressure(ann["P_static_MPa"], ann["P_friction_MPa"], win,
                                    ecd["TVD_cum_guarded_m1b"], segs, nt)

    # S18（:1608-1609）；Cumulative_Volume_m3 = (Σ 逐时累加)/1000
    cumvol1 = np.asarray([sched.vol_all_L_1b[t] / 1000.0 for t in range(1, nt + 1)])
    rm = np.zeros((nt, 5))
    rm[:, 0] = np.arange(1, nt + 1, dtype=float)
    rm[:, 1] = cumvol1
    rm[:, 2] = ecd["ECD_casing_g_cm3"][1:, segs.n]
    rm[:, 3] = cas["pump_surface_MPa"][1:]
    rm[:, 4] = ann["P_total_MPa"][1:, segs.n]

    pump_cmp = compare_with_reference(cas["pump_surface_MPa"], _pad1(cumvol1),
                                      sandbox / "HT1-004_8_11_pump_pressure_reference.csv",
                                      "design_pressure_MPa", nt)
    ecd_bot1 = _pad1(ecd["ECD_annulus_g_cm3"][1:, segs.n])
    ecd_cmp = compare_with_reference(ecd_bot1, _pad1(cumvol1),
                                     sandbox / "HT1-004_8_11_bottom_ecd_reference.csv",
                                     "design_bottom_ECD_g_cm3", nt)

    fp = np.zeros((nt, 5))
    fp[:, 0] = np.arange(1, nt + 1, dtype=float)
    fp[:, 1] = corridor["required_backpressure_MPa"][1:]
    fp[:, 2] = corridor["BP_lower_MPa"][1:]
    fp[:, 3] = corridor["BP_upper_MPa"][1:]
    fp[:, 4] = corridor["conflict_flag"][1:]
    fe = np.zeros((nt, 5))
    fe[:, 0] = np.arange(1, nt + 1, dtype=float)
    fe[:, 1] = corridor["ECD_5568_new"][1:]
    fe[:, 2] = corridor["ECD_6880_new"][1:]
    fe[:, 3] = corridor["ECD_7463_new"][1:]
    fe[:, 4] = corridor["ECD_bottom_new"][1:]
    pps = np.column_stack([np.arange(1, nt + 1, dtype=float), cas["pump_surface_MPa"][1:]])
    abp = np.column_stack([np.arange(1, nt + 1, dtype=float), ann["P_total_MPa"][1:, segs.n]])
    wbc = np.zeros((nt, 5))
    wbc[:, 0] = np.arange(1, nt + 1, dtype=float)
    wbc[:, 1] = fz.ecd_base_bottom[1:]
    wbc[:, 2] = fz.ecd_bottom_ctrl[1:]
    wbc[:, 3] = fz.ecd_shoe_ctrl[1:]
    wbc[:, 4] = fz.p_wh_ctrl[1:]

    pump1 = cas["pump_surface_MPa"][1:]
    summary = {
        "n_time": nt,
        "n_segment": segs.n,
        "TVD_bottom_m": float(win["TVD_bottom_m"]),
        "max_pump_MPa": float(np.max(pump1)),
        "min_pump_MPa": float(np.min(pump1)),
        "max_annuli_bottom_MPa": float(np.max(ann["P_total_MPa"][1:, segs.n])),
        "max_bottom_ECD_casing_g_cm3": float(np.max(rm[:, 2])),
        "min_bottom_ECD_casing_g_cm3": float(np.min(rm[:, 2])),
        "conflict_steps": int(np.sum(corridor["conflict_flag"][1:] != 0)),
        "freeze_segments_1b": list(zip(fz.seg_s, fz.seg_e)),
        "bp_apply_at_last_min_MPa": float(fz.bp_apply_MPa[nt]),
        "ECD_setpoint_win": 2.010 + 0.015,
    }
    return {
        "out_result_matrix_volume_pressure": rm,
        "out_pump_pressure_comparison": pump_cmp,
        "out_bottom_ecd_comparison": ecd_cmp,
        "out_four_point_backpressure": fp,
        "out_four_point_ecd_new": fe,
        "out_pump_pressure_surface": pps,
        "out_annuli_bottom_pressure": abp,
        "out_ecd_casing_full": ecd["ECD_casing_g_cm3"][1:, 1:],
        "out_window_backpressure_ctrl": wbc,
        "cumulative_volume_m3": cumvol1,
        "bp_apply_win_MPa": fz.bp_apply_MPa[1:],
        "ecd_annulus_g_cm3_full": ecd["ECD_annulus_g_cm3"][1:, 1:],
        "summary": summary,
        # 方向性/审计辅助
        "pump_Lmin_1b": sched.pump_Lmin_1b,
        "annulus_static_MPa": ann["P_static_MPa"],
        "tvd_cum_m1b": ecd["TVD_cum_guarded_m1b"],
        "pr_seg": segs.pr_seg,
    }


# ==========================================================================
# P-4 泛化附加块（Phase 3.2，2026-10-07）——纯加法：
# 不改动上方任何函数；run_ht1004_target 数值路径零改动（契约测试钉住）。
# 用途 = 呼101/呼1-003 正演产品口径（无 MATLAB 靶，不背原样复现红线，
# 但公式/乘序/单位镜像内核，与靶链同口径可比）。
# ==========================================================================


def fluid_table_from(
    rates_Lmin, vols_L, rou_gcc, miu_mPas, tau_Pa, rou0, miu0, tau0,
    dt_min=1.0, bp_MPa=None,
):
    """通用流体表构造（N 路；前 5 路必须 = 进环空流体，S8/S9 五界面模型硬约束）。"""
    n = len(rates_Lmin)
    if n < 5:
        raise ValueError("环空侧按 5 界面建模（S8/S9），流体路数必须 >= 5")
    for name, seq in (("vols_L", vols_L), ("rou_gcc", rou_gcc),
                      ("miu_mPas", miu_mPas), ("tau_Pa", tau_Pa)):
        if len(seq) != n:
            raise ValueError(name + " 长度与 rates 不符")
    bp = tuple(0.0 for _ in range(n)) if bp_MPa is None else tuple(float(x) for x in bp_MPa)
    if len(bp) != n:
        raise ValueError("bp_MPa 长度与路数不符")
    return FluidTable(
        dt_min=dt_min,
        rates_Lmin=tuple(float(x) for x in rates_Lmin),
        vols_L=tuple(float(x) for x in vols_L),
        rou0=float(rou0), miu0=float(miu0), tau0=float(tau0),
        rou_gcc=tuple(float(x) for x in rou_gcc),
        miu_mPas=tuple(float(x) for x in miu_mPas),
        tau_Pa=tuple(float(x) for x in tau_Pa),
        bp_MPa=bp,
    )


def build_pump_schedule_custom(fluids):
    """S7 通用版（N 路）。分支序镜像 :291-335 边界归属习惯：
    第 1 路 t<node1；第 2/3 路以 >= 起；第 4 路起以 > 起（先命中先赢）；
    越过末节点 = 停泵（pump=0，与原样一致）。无 :733 类靶专属偏置。"""
    rates = fluids.rates_Lmin
    vols = fluids.vols_L
    n_f = len(rates)
    total_min = 0.0
    for v, r in zip(vols, rates):
        total_min += v / r
    n_time = int(math.floor(total_min / fluids.dt_min))
    nodes = []
    acc = 0.0
    for v, r in zip(vols, rates):
        acc += v / r
        nodes.append(acc / fluids.dt_min)
    pump = np.zeros(n_time + 1)
    bp = np.zeros(n_time + 1)
    vol_all = np.zeros(n_time + 1)
    q = np.zeros(n_time + 1)
    pump[1] = rates[0]
    bp[1] = fluids.bp_MPa[0]
    vol_all[1] = pump[1] * fluids.dt_min
    for t in range(2, n_time + 1):
        tv = float(t)
        k = None
        if tv < nodes[0]:
            k = 0
        elif n_f >= 2 and tv >= nodes[0] and tv <= nodes[1]:
            k = 1
        elif n_f >= 3 and tv >= nodes[1] and tv <= nodes[2]:
            k = 2
        else:
            for j in range(3, n_f):
                if tv > nodes[j - 2] and tv <= nodes[j - 1]:
                    k = j
                    break
        if k is None:
            pump[t] = 0.0
            bp[t] = 0.0
        else:
            pump[t] = rates[k]
            bp[t] = fluids.bp_MPa[k]
        vol_all[t] = vol_all[t - 1] + pump[t] * fluids.dt_min
    q[1:] = np.asarray([((float(p) / 60.0) / 1000.0) for p in pump[1:]])
    return PumpSchedule(
        n_time=n_time, dt_min=fluids.dt_min, nodes_min=nodes,
        pump_Lmin_1b=pump, bp_MPa_1b=bp, vol_all_L_1b=vol_all, q_m3s_1b=q,
    )


def build_segments_custom(structure, thresholds_m, bore_const_mm, pipe_od_mm,
                          pipe_wall_mm, pianxin=True):
    """S2/S3/S6 通用版。K = len(thresholds)+1 段（靶 = 7 段硬编码；此处参数化）。

    - thresholds_m：升序变径/变段分界 MD [m]；
    - bore_const_mm：逐段井眼/套管内径常数 [mm]，None = 取 CSV 名义直径列（靶段 4-7 语义）。
      靶语义注记：段 1-3 常数 245.37 = 273.1 套管 ID —— 传值一律按**直径**口径
      （legacy "radius 列实为直径"陷阱不在此参数上）；
    - pipe_od_mm / pipe_wall_mm：逐段管串外径/壁厚 [mm]（ID = OD - 2*wall）。
    其余几何/容积链/PR 运算与 build_segments 同式同序。
    """
    thr = [float(x) for x in thresholds_m]
    K = len(thr) + 1
    if not (len(bore_const_mm) == len(pipe_od_mm) == len(pipe_wall_mm) == K):
        raise ValueError("thresholds/常数表段数不符（K = len(thresholds)+1）")
    n = structure.n_segment
    cd = structure.c_depth_m
    masks = []
    for j, t in enumerate(thr):
        lo = thr[j - 1] if j > 0 else float("-inf")
        masks.append((cd > lo) & (cd <= t))
    masks.append(cd > thr[-1])
    d_bit = np.zeros(n)
    d_out = np.zeros(n)
    d_in = np.zeros(n)
    for k, m in enumerate(masks):
        bc = bore_const_mm[k]
        if bc is not None:
            d_bit[m] = float(bc) * 0.001
        else:
            d_bit[m] = structure.hole_dia_cm[m] * 0.01
        od = float(pipe_od_mm[k])
        wt = float(pipe_wall_mm[k])
        d_out[m] = od * 0.001
        d_in[m] = (od - wt * 2) * 0.001
    area_cout = math.pi * (d_bit**2 - d_out**2) / 4.0
    area_cin = math.pi * (d_in**2) / 4.0
    out_diam_bole = structure.hole_dia_cm * 10.0
    vol_cas_seg = area_cin * structure.seg_len_m
    capacity_m3 = _sum_seq(vol_cas_seg)
    capacity_L = capacity_m3 * 1000.0
    cos_deg = np.array([math.cos(math.radians(float(dg))) for dg in structure.deg])
    vert_len = np.abs(structure.seg_len_m * cos_deg)
    tvd_cum = np.asarray(_cumsum_seq(vert_len))
    d_i = d_out.copy()
    d_o = structure.hole_dia_cm / 100.0
    vab = np.zeros(n)
    for i in range(n):
        vab[i] = _sum_seq(structure.vol_annulus_L[i:])
    vtop = np.zeros(n)
    prod_cin = area_cin * structure.seg_len_m
    for i in range(n):
        vtop[i] = _sum_seq(prod_cin[: i + 1]) * 1000.0
    pr_seg = []
    for m in masks:
        if np.any(m):
            mean_dole = _mean_seq(out_diam_bole[m])
            mean_cas = _mean_seq(d_out[m]) * 1000.0
            pr_seg.append(friction_pr(mean_dole, mean_cas))
        else:
            pr_seg.append(0.0)
    pr = np.zeros(n)
    for i in range(n):
        c = float(cd[i])
        pr[i] = pr_seg[K - 1]
        for j, t in enumerate(thr):
            if c <= t:
                pr[i] = pr_seg[j]
                break
    a_annulus = structure.sq_annulus_dm2 / 100.0
    return Segments(
        n=n,
        c_depth_m1b=_pad1(cd),
        seg_len_m1b=_pad1(structure.seg_len_m),
        sq_annulus_dm2_1b=_pad1(structure.sq_annulus_dm2),
        d_bit_m1b=_pad1(d_bit),
        d_cas_out_m1b=_pad1(d_out),
        d_cas_in_m1b=_pad1(d_in),
        out_diam_bole_mm1b=_pad1(out_diam_bole),
        area_cout_m2_1b=_pad1(area_cout),
        area_cin_m2_1b=_pad1(area_cin),
        d_i_m1b=_pad1(d_i),
        d_o_m1b=_pad1(d_o),
        A_annulus_m2_1b=_pad1(a_annulus),
        cos_deg_1b=_pad1(cos_deg),
        vert_len_m1b=_pad1(vert_len),
        tvd_cum_m1b=_pad1(tvd_cum),
        vab_L1b=_pad1(vab),
        vtop_L1b=_pad1(vtop),
        capacity_pipe_L=capacity_L,
        seg_masks=masks,
        pr_seg=pr_seg,
        pr_1b=_pad1(pr),
    )


def track_casing_interfaces_custom(sched, segs, fluids):
    """S11 通用版（N 界面）。阈值链同 :738-774；不携带 :733 +1200 L 靶专属偏置
    （原样复现红线仅约束 run_ht1004_target；正演产品口径按干净语义）。"""
    nt = sched.n_time
    n_f = len(fluids.vols_L)
    dt = sched.dt_min
    vol_c = np.zeros((n_f, nt + 1))
    thr = [0.0]
    acc = 0.0
    for k in range(1, n_f):
        acc += fluids.vols_L[k - 1]
        thr.append(acc)
    for t in range(2, nt + 1):
        for k in range(n_f):
            if sched.vol_all_L_1b[t] > thr[k]:
                vol_c[k, t] = vol_c[k, t - 1] + sched.pump_Lmin_1b[t] * dt
    tags = np.zeros((nt, n_f))
    res_h = np.zeros((nt, n_f))
    for t in range(1, nt + 1):
        for k in range(n_f):
            tag, _rv, rh = _casing_interface_scan(
                float(vol_c[k, t]), segs.vtop_L1b, segs.area_cin_m2_1b, segs.n)
            tags[t - 1, k] = tag
            res_h[t - 1, k] = rh
    return vol_c, tags, res_h


def assign_casing_properties_custom(tags, res_h, sched, segs, fluids):
    """S11 属性分配通用版（N 流体链，方向对偶同 :1119-1279；无 :1153 legacy 分支）。"""
    nt, n = sched.n_time, segs.n
    N = len(fluids.rou_gcc)
    R = list(fluids.rou_gcc)
    M = list(fluids.miu_mPas)
    T = list(fluids.tau_Pa)
    R0, M0, T0 = fluids.rou0, fluids.miu0, fluids.tau0
    plen = segs.seg_len_m1b
    acin = segs.area_cin_m2_1b
    pump = sched.pump_Lmin_1b
    rou = np.zeros((nt + 1, n + 1))
    miu = np.zeros((nt + 1, n + 1))
    tau = np.zeros((nt + 1, n + 1))
    velo = np.zeros((nt + 1, n + 1))
    for t in range(1, nt + 1):
        tc = tags[t - 1]
        for i in range(1, n + 1):
            if i > tc[0]:
                rou[t, i], miu[t, i], tau[t, i] = R0, M0, T0
            elif i == tc[0]:
                p = min(max(res_h[t - 1, 0] / plen[i], 0.0), 1.0)
                q = 1.0 - p
                rou[t, i] = p * R0 + q * R[0]
                miu[t, i] = p * M0 + q * M[0]
                tau[t, i] = p * T0 + q * T[0]
            else:
                placed = False
                for k in range(1, N):
                    if tc[k] < i < tc[k - 1]:
                        rou[t, i], miu[t, i], tau[t, i] = R[k], M[k], T[k]
                        placed = True
                        break
                    if i == tc[k]:
                        p = min(max(res_h[t - 1, k] / plen[i], 0.0), 1.0)
                        q = 1.0 - p
                        rou[t, i] = p * R[k] + q * R[k + 1]
                        miu[t, i] = p * M[k] + q * M[k + 1]
                        tau[t, i] = p * T[k] + q * T[k + 1]
                        placed = True
                        break
                if not placed and i < tc[N - 1]:
                    rou[t, i], miu[t, i], tau[t, i] = R[N - 1], M[N - 1], T[N - 1]
            if pump[t] > 0:
                velo[t, i] = ((pump[t] / 1000.0) / 60.0) / acin[i]
            else:
                velo[t, i] = 0.0
    return {"rou_gcc": rou, "miu_mPas": miu, "tau_Pa": tau, "velo_m_s": velo}


def run_well_forward(structure, segs, fluids, sched):
    """P-4 正演装配（复用靶内核全链；无 S16/S19-S22/S25——那些是呼1-004 靶专属窗口/对比）。

    返回：泵压时程 [MPa]、环空井底压力 [MPa]、井底环空/管内 ECD [g/cm3]、
    cumvol [m3]、summary 标量组。
    """
    nt, n = sched.n_time, segs.n
    vol_a = annulus_entry_volumes(sched, segs, fluids)
    tags_a, res_h_a = track_annulus_interfaces(vol_a, segs, _pad1(structure.deg))
    props_a = assign_annulus_properties(tags_a, res_h_a, sched, segs, fluids)
    ann = hydrostatic_and_friction_chain(segs, props_a, sched)
    vol_c, tags_c, res_h_c = track_casing_interfaces_custom(sched, segs, fluids)
    props_c = assign_casing_properties_custom(tags_c, res_h_c, sched, segs, fluids)
    st_c = casing_static_chain(props_c, segs, nt)
    cas = pump_pressure_backtrack(segs, props_c, sched, ann["P_total_Pa"], st_c["P_static_Pa"])
    ecd = compute_ecd_esd(ann["P_total_Pa"], ann["P_static_Pa"], cas["P_casing_Pa"],
                          segs.tvd_cum_m1b, nt, n)
    pump1 = cas["pump_surface_MPa"][1:]
    ecd_bot = ecd["ECD_annulus_g_cm3"][1:, n]
    ecd_cas_bot = ecd["ECD_casing_g_cm3"][1:, n]
    p_bot = ann["P_total_MPa"][1:, n]
    cumvol1 = np.asarray([sched.vol_all_L_1b[t] / 1000.0 for t in range(1, nt + 1)])
    summary = {
        "n_time": nt,
        "n_segment": n,
        "n_fluid_paths": len(fluids.vols_L),
        "min_pump_MPa": float(np.min(pump1)),
        "max_pump_MPa": float(np.max(pump1)),
        "max_annuli_bottom_MPa": float(np.max(p_bot)),
        "min_bottom_ECD_g_cm3": float(np.min(ecd_bot)),
        "max_bottom_ECD_g_cm3": float(np.max(ecd_bot)),
        "min_bottom_ECD_casing_g_cm3": float(np.min(ecd_cas_bot)),
        "max_bottom_ECD_casing_g_cm3": float(np.max(ecd_cas_bot)),
        "capacity_pipe_L": float(segs.capacity_pipe_L),
        "TVD_bottom_m": float(segs.tvd_cum_m1b[n]),
        "nan_inf_count": int(np.sum(~np.isfinite(pump1)) + np.sum(~np.isfinite(ecd_bot))),
    }
    return {
        "pump_MPa": pump1, "ann_bottom_MPa": p_bot,
        "ecd_bottom_annulus_g_cm3": ecd_bot, "ecd_bottom_casing_g_cm3": ecd_cas_bot,
        "cumvol_m3": cumvol1, "ecd_ann_full": ecd["ECD_annulus_g_cm3"][1:, 1:],
        "ecd_cas_full": ecd["ECD_casing_g_cm3"][1:, 1:],
        "flow_pattern_ann": ann["flow_pattern"][1:, 1:],
        "sched_nodes_min": sched.nodes_min, "summary": summary,
    }
