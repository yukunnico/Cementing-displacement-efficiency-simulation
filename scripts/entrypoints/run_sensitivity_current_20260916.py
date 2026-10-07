"""当前口径敏感性重跑（2026-09-16）——为论文第 4 章取当前模型数字。

背景
----
论文第 3/4 章现有数字全部来自 2026-09-09/09-10 口径（求解器旧路径，HEAD 早于
2026-09-14 源模型口径重构 Task 0–14）。重构把环空二维求解切到流函数新路径
（`enable_stream_function=True` 默认），权威 8 井数字随之改变
（`results/源模型口径重跑_2026-09-14/汇总.csv`）。本脚本在**当前口径**下重跑
同一批单变量变体，供论文数字替换。

口径（与 `results/源模型口径重跑_2026-09-14` 严格同源）
------------------------------------------------------
- 求解器：AnnulusD2DGASolver(total_t, nz=250, enable_cfl_adaptive=True, **CORRECTED_KW)
  —— CORRECTED_KW / NZ / _stop_t / _total_t 直接 import 自
  `scripts/entrypoints/rerun_all_wells_corrected.py`，避免口径漂移。
- 停止时刻：tt = min(总泵注时长 + 1200 s, 尾浆入库时刻)（F2 口径，无 +600 s 尾窗）。
- 套管内 1D：CasingFlowSolver(enable_gravity=True, mixing_contact_time=True,
  plug_face_zero_mixing=True, has_plug=True)（T1 生产口径）。
- 基线：`results/源模型口径重跑_2026-09-14/单井结果/<well>_corrected_on.json`
  （当前口径权威基线，**不**用 results/<井名>_1D2D耦合模型/ 的旧路径数字）。

变体
----
12 个现场杠杆（与 2026-09-10 敏感性组逐个同名，便于新旧逐格对照）：
    rate_x0.6 / rate_x0.8 / rate_x1.2 / rate_x1.4      排量缩放（泵序逐段）
    cement_n_p0.1 / cement_n_m0.1                       水泥浆 n ±0.1
    spacer_dens_p100 / spacer_dens_m100                 隔离液密度 ±100 kg/m³
    mud_pv_x0.7 / mud_pv_x1.3                           钻井液 PV ±30%
    standoff_p0.1 / standoff_m0.1                       standoff 剖面整体 ±0.1
另加 standoff 连续扫描 7 档（±0.05 / ±0.15 / ±0.20 / ±0.30），与 ±0.1 两档合成
9 档连续曲线，用于判定「居中度阈值悬崖」在当前口径下是否仍然存在。

输出
----
results/敏感性变体_当前口径_2026-09-16/（逐变体摘要 JSON + 汇总 CSV/MD）。
⚠️ **该目录已冻结**（2026-10-01 起）：main() 只保留给历史复现，实际执行会被
目录内容守卫拒绝（见 main() 内 RuntimeError）。本文件的**装配层**（变体定义、
run_variant、温度场构造、判别量）是 09-27 两下游批与 T2 温压敏感性批的共享源。

T2 温压敏感性扩展（2026-10-01，T2-1a 装配面）
-----------------------------------------------
- ``build_variants()`` 元组由 4 元扩为 **5 元**，末位 ``run_opts: dict | None``
  （legacy 19 杠杆全为 ``None`` = 纯 T-off）。**故意不给默认值**：旧的 4 元解包
  会响亮 ValueError，防下游漏适配。
- ``run_variant(loader, well_fn, fluid_fn, sched_fn, run_opts, well_key="")``
  按 run_opts 决定温度开关/温度档/屈服门，并返回 **(summary, extra_metrics)**
  ——summary 原样（一个键都不加，保 T-off 字节逐位），判别量走敏感性层新表。
- 2026-10-01（T2 时程观察层）：内核抽出为
  ``run_variant_res(...) -> (res, cr, schedule2, extra)``，``run_variant`` 改为其
  薄壳（数值路径逐行不变，供过程观察层取时程快照与鞋口事件时刻）。
- 新增 ``build_temperature_fields(well_key, mode)`` / ``extra_metrics(...)``。
  详见 T2 批 scripts/entrypoints/run_sensitivity_temperature_t2_20261001.py。

用法
----
    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/entrypoints/run_sensitivity_current_20260916.py
"""
from __future__ import annotations

import csv
import importlib
import json
import sys
import time
from dataclasses import fields, replace
from pathlib import Path

import numpy as np

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.pressure_field import HydrostaticPressureField, insitu_column_density
from cemdisp.data.rheology_vs_temperature import RheologyFormulaParams
from cemdisp.data.pumping_schedule import PumpingSchedule
from cemdisp.data.temperature_field import (
    AnchoredProfileField,
    ConstantTemperatureField,
    GeothermalTemperatureField,
    load_delivered_pair,
    load_extended_pair_4d,
)
from cemdisp.data.well_spec import WellSpec
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver

# 口径唯一来源：与权威 8 井重跑共用同一组常量与函数，防止漂移。
_SCRIPTS_DIR = Path(__file__).resolve().parents[1]
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))
from entrypoints.rerun_all_wells_corrected import (  # noqa: E402
    CORRECTED_KW,
    NZ,
    _stop_t,
    _total_t,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = PROJECT_ROOT / "results" / "敏感性变体_当前口径_2026-09-16"
BASELINE_DIR = PROJECT_ROOT / "results" / "源模型口径重跑_2026-09-14" / "单井结果"

# 井名 → (loader 模块, loader 函数, 当前口径基线 JSON 名)
WELLS: dict[str, tuple[str, str, str]] = {
    "呼103": ("cemdisp.data.loaders.hu103_loader", "load_hu103_tailpipe", "hu103_corrected_on.json"),
    "呼101": ("cemdisp.data.loaders.hu101_loader", "load_hu101_tailpipe", "hu101_corrected_on.json"),
}

CEMENT_ROLES = {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}


# ---------------------------------------------------------------- 输入变换
def scale_schedule(schedule: PumpingSchedule, factor: float) -> PumpingSchedule:
    """泵序逐段排量缩放（体积不变，时长随 1/factor 拉长）。"""
    return replace(schedule, steps=tuple(
        replace(s, rate_m3_min=s.rate_m3_min * factor) if s.rate_m3_min > 0 else s
        for s in schedule.steps
    ))


def shift_cement_n(fluids: tuple[FluidSpec, ...], delta: float) -> tuple[FluidSpec, ...]:
    """水泥浆幂律指数 n ± delta（仅水泥三角色且 power_law_n 非空时生效）。"""
    return tuple(
        replace(f, power_law_n=f.power_law_n + delta)
        if f.role in CEMENT_ROLES and f.power_law_n is not None else f
        for f in fluids
    )


def shift_spacer_density(fluids: tuple[FluidSpec, ...], d_kg_m3: float) -> tuple[FluidSpec, ...]:
    """隔离液密度 ± d_kg_m3（仅 SPACER 角色）。"""
    return tuple(
        replace(f, density_kg_m3=f.density_kg_m3 + d_kg_m3) if f.role == FluidRole.SPACER else f
        for f in fluids
    )


def scale_mud_pv(fluids: tuple[FluidSpec, ...], factor: float) -> tuple[FluidSpec, ...]:
    """钻井液塑性黏度 × factor（仅 MUD 角色）。"""
    return tuple(
        replace(f, plastic_viscosity_pa_s=f.plastic_viscosity_pa_s * factor)
        if f.role == FluidRole.MUD and f.plastic_viscosity_pa_s is not None else f
        for f in fluids
    )


def shift_standoff(well: WellSpec, delta: float) -> WellSpec:
    """standoff 剖面整体平移 delta（clip [0,1]），形状与深度点不变。"""
    from cemdisp.data.well_spec import DepthValuePoint

    if not well.standoff_profile:
        raise ValueError(f"{well.well_name}: 无 standoff_profile")
    return replace(well, standoff_profile=tuple(
        DepthValuePoint(p.depth_md_m, float(np.clip(p.value + delta, 0.0, 1.0)))
        for p in well.standoff_profile
    ))


def _identity(x):
    return x


# ---------------------------------------------------------------- 变体表
def build_variants() -> list[tuple[str, object, object, object, dict | None]]:
    """(变体名, well 变换, fluids 变换, schedule 变换, run_opts)。

    5 元组：末位 ``run_opts``（legacy 19 杠杆全为 ``None`` ⇒ 纯 T-off，
    与本文件 2026-09-16 首版行为逐位一致）。刻意不给默认值，令旧 4 元解包
    在下游响亮 ValueError（防漏适配）。
    """
    return [
        # 12 杠杆（与 2026-09-10 组同名，逐格可比）
        ("rate_x1.2", _identity, _identity, lambda s: scale_schedule(s, 1.2), None),
        ("rate_x0.8", _identity, _identity, lambda s: scale_schedule(s, 0.8), None),
        ("rate_x1.4", _identity, _identity, lambda s: scale_schedule(s, 1.4), None),
        ("rate_x0.6", _identity, _identity, lambda s: scale_schedule(s, 0.6), None),
        ("cement_n_p0.1", _identity, lambda fs: shift_cement_n(fs, +0.1), _identity, None),
        ("cement_n_m0.1", _identity, lambda fs: shift_cement_n(fs, -0.1), _identity, None),
        ("spacer_dens_p100", _identity, lambda fs: shift_spacer_density(fs, +100.0), _identity, None),
        ("spacer_dens_m100", _identity, lambda fs: shift_spacer_density(fs, -100.0), _identity, None),
        ("mud_pv_x1.3", _identity, lambda fs: scale_mud_pv(fs, 1.3), _identity, None),
        ("mud_pv_x0.7", _identity, lambda fs: scale_mud_pv(fs, 0.7), _identity, None),
        ("standoff_p0.1", lambda w: shift_standoff(w, +0.1), _identity, _identity, None),
        ("standoff_m0.1", lambda w: shift_standoff(w, -0.1), _identity, _identity, None),
        # standoff 连续扫描补充档
        ("standoff_m0.05", lambda w: shift_standoff(w, -0.05), _identity, _identity, None),
        ("standoff_p0.05", lambda w: shift_standoff(w, +0.05), _identity, _identity, None),
        ("standoff_m0.15", lambda w: shift_standoff(w, -0.15), _identity, _identity, None),
        ("standoff_m0.20", lambda w: shift_standoff(w, -0.20), _identity, _identity, None),
        ("standoff_p0.20", lambda w: shift_standoff(w, +0.20), _identity, _identity, None),
        ("standoff_m0.30", lambda w: shift_standoff(w, -0.30), _identity, _identity, None),
        ("standoff_p0.30", lambda w: shift_standoff(w, +0.30), _identity, _identity, None),
    ]


# ---------------------------------------------------------------- T2 装配面
# run_opts 键（R1 2026-10-06 由三键扩为五键；**默认值 = 旧行为**）：
#   enable_temperature_rheology : 总开关（False ⇒ 数值路径与 HEAD 逐位）
#   temperature_mode            : "off" | "static" | "table" | "const60"
#   enable_yield_gate           : None=沿用 CORRECTED_KW（True）/ False=Q10乙通道代理
#   rheology_formula            : None=默认公式系数；dict 见 `params_from_spec`
#                                 （R1 公式系数扰动档；L2 流变响应验收的操作面）
#   mud_extrapolate             : False=钻井液 T>80°C clamp 到域端（现状）；
#                                 True=按公式外推（深段泥浆 clamp 对照档）
#   pressure_mode               : "off"=不注入压力场（P=None ⇒ p_default 审计，逐位=HEAD）；
#                                 "hydrostatic"=注入 `HydrostaticPressureField`（P-1）
#   pressure_caliber            : "shoe"（默认，域底单点）/ "mean"（域内均值对照档）
#   enable_stream_yield_gate    : None=**不给键**（沿用 annulus 构造默认 False，逐位=现状）；
#                                 True/False=显式出键覆盖。
#     ⚠ 语义区分（勿与 enable_yield_gate 混）：
#       `enable_yield_gate`       = **生产端**：wall 场算不算（默认 True；T2 的 gateoff 对照=本键）
#       `enable_stream_yield_gate`= **消费端**：wall 进不进流函数算子（默认 False；
#                                   唯一扯断点 annulus 三目，开启即通）
#     二者组合语义见 annulus `_dead_switches`（置真但 wall 到不了算子会判"死开关"告警）。
#     ⚠ 仅 annulus 有本开关（casing 无）⇒ 只准 `annulus_kwargs_from_opts` 出此键。
#   include_yield_term          : **R2（Phase 2）μp/τy 真拆分**（默认 False ⇒ 逐位=HEAD）；
#                                 True ⇒ 标量黏度口径补 τy/γ̇ 并与场口径同构（含 3.0 Pa·s 上限）。
#                                 仅 annulus 有该形参 ⇒ 只准 `annulus_kwargs_from_opts` 出此键。
#                                 依据（文献+资料）见 Phase 2 spec §1。
RUN_OPTS_KEYS = ("enable_temperature_rheology", "temperature_mode",
                 "enable_yield_gate", "rheology_formula", "mud_extrapolate",
                 "pressure_mode", "pressure_caliber", "enable_stream_yield_gate",
                 "include_yield_term")
TEMPERATURE_MODES = ("off", "static", "table", "const60", "anchored", "table_ext")
PRESSURE_MODES = ("off", "hydrostatic")
PRESSURE_CALIBERS = ("shoe", "mean")
# 有交付瞬态温度表的井（table 档专用；呼101/呼103 无表 ⇒ 只能 static/const60/off）
TABLE_WELLS = frozenset({"呼1-004"})

# ---------------------------------------------------------------- 4d 锚定静温场
# 「anchored」档（Phase 4d 主批，spec `2026-10-07-phase4d-seven-well-batch-design.md` §1）
# 锚点族 = **场景 A 电测/作业史静温**；逐行出处见
# `results/_probe_4d七井_2026-10-07/七井锚点盘点表.md`（机读盘点脚本
# `anchor_inventory_20261007.py`，在三重点井上逐字复现其在 `_probe_4d锚点_2026-10-07`
# 的三井对账表场景 A 锚点集）。锚点**内联为常量**：不在运行时读 `参考文档/`（版本控制外）。
#
# k 取值：三重点井 = 用户 2026-10-07 裁定；非重点井 = 用户 2026-10-07 裁定
# 「notes 明示优先」（呼探1-002 取井底主段明示值 0.90，见 spec §1.1-1 登记）。
# 值一律落在 [0.80, 0.95] 观测带内（计划 §8-13 未触发）。
ANCHORED_WELLS: dict[str, tuple[float, tuple[tuple[float, float], ...], str]] = {
    "呼101": (0.90, ((5700.0, 123.0), (7868.0, 150.0)),
              "无 notes 明示系数；观测区间 argmin（用户裁定）"),
    "呼1-003": (0.85, ((5290.0, 123.0), (7618.0, 150.0)),
                "notes 明示「领浆温度系数0.85」（用户裁定）"),
    "呼1-004": (0.85, ((5241.0, 124.0), (7660.0, 155.0)),
                "notes 明示「领浆温度系数0.85」（用户裁定）"),
    "呼102": (0.90, ((7120.0, 147.8), (7735.0, 149.0)),
              "notes 明示「水泥浆试验温度(0.9x温度系数)」（用户裁定：明示优先）"),
    "呼探1-002": (0.90, ((5292.5, 111.0), (7554.0, 148.0)),
                  "notes 明示「循环温度133℃(系数0.90)」井底主段（用户裁定：明示优先）"),
    "呼探1-001": (0.85, ((5460.159, 110.0), (5900.0, 118.0), (7000.0, 137.0), (7746.0, 150.0)),
                  "notes 明示「温度系数0.85」（用户裁定：明示优先）"),
}
# 无锚井（唯一锚深度 < 2，AnchoredProfileField 不可构造）⇒ Geothermal 回退 + 口径声明。
# 呼探1 的 `temperature_pressure_profile.csv` 只有 3 个带温度行且同落在 md=7601 m
# （153.8/159/167 °C），无 nd 维可分段（盘点表 §2「呼探1」节）。
ANCHORED_NO_ANCHOR_WELLS = frozenset({"呼探1"})

# 呼1-004 时程扩展温度表（333×362，表末 21600 s；4d 扩表产物）。
# ⚠ `*_ext4d.xlsx` 受 `.gitignore` 的 `**/*.xlsx` 管辖**不入库**——由已入库的
# `run_ext4d.m` + `HT1_004_T.m` 可重生成（MATLAB R2025b，见 `_probe_4d扩表_2026-10-07/对账报告_4d扩表.md`）。
# 缺失即响亮报错：**不静默降级到交付表**（交付表末 11940 s，会把 r0.6/0.8 的 stop_t
# clamp 掉——正是 4d 要消除的污染）。
EXT4D_DIR = (Path(__file__).resolve().parents[2] / "results"
             / "_probe_4d扩表_2026-10-07" / "sandbox" / "HT1-004压力计算")
EXT4D_T_IN_XLSX = EXT4D_DIR / "T_in_ext4d.xlsx"
EXT4D_T_OUT_XLSX = EXT4D_DIR / "T_out_ext4d.xlsx"
EXT4D_TIME_AXIS_MIN = EXT4D_DIR / "T_ext4d_time_axis_min.csv"
EXT4D_WELLS = frozenset({"呼1-004"})   # 扩展表只覆盖呼1-004（同交付表）
TABLE_EXT_WELLS = EXT4D_WELLS


_PARAMS_FIELDS = frozenset(f.name for f in fields(RheologyFormulaParams))
"""`RheologyFormulaParams` 的合法字段名（`params_from_spec` 校验用）。"""


def params_from_spec(spec: dict | None) -> RheologyFormulaParams | None:
    """``run_opts["rheology_formula"]`` 载荷 → `RheologyFormulaParams`（R1）。

    载荷（均可省略，``None`` ⇒ 返回 ``None`` = 用 `fluid_at` 默认参）::

        {"scale": {"<字段名>": 乘子, ...},   # 逐字段缩放：tuple 字段逐元素乘、float 字段直接乘
         "set":   {"<字段名>": 值, ...}}      # 逐字段直接替换（在 scale 之后施加，优先级更高）

    未知字段名 / 空载荷 / 非数值一律**响亮报错**（防"配错了但静默按默认跑"）。
    例：组 A τ₀ 系数 +10% ⇒ ``{"scale": {"ca_tau0_q": 1.1}}``。
    """
    if spec is None:
        return None
    if not isinstance(spec, dict):
        raise TypeError(f"rheology_formula 须为 dict|None，实际 {type(spec).__name__}")
    unknown_keys = sorted(set(spec) - {"scale", "set"})
    if unknown_keys:
        raise KeyError(f"rheology_formula 含未知子键 {unknown_keys}，允许：['scale', 'set']")
    changes: dict = {}
    for sub in ("scale", "set"):
        table = spec.get(sub)
        if table is None:
            continue
        if not isinstance(table, dict):
            raise TypeError(f"rheology_formula[{sub!r}] 须为 dict，实际 {type(table).__name__}")
        for field_name, value in table.items():
            if field_name not in _PARAMS_FIELDS:
                raise KeyError(
                    f"rheology_formula[{sub!r}] 含未知字段 {field_name!r}；"
                    f"合法字段：{sorted(_PARAMS_FIELDS)}"
                )
            current = getattr(RheologyFormulaParams(), field_name)
            if sub == "scale":
                factor = float(value)
                if isinstance(current, tuple):
                    changes[field_name] = tuple(float(c) * factor for c in current)
                else:
                    changes[field_name] = float(current) * factor
            else:
                if isinstance(current, tuple):
                    changes[field_name] = tuple(float(c) for c in value)
                else:
                    changes[field_name] = float(value)
    if not changes:
        raise ValueError("rheology_formula 载荷为空（既无 scale 也无 set）")
    return replace(RheologyFormulaParams(), **changes)


def normalize_run_opts(run_opts: dict | None) -> dict:
    """补默认 + 校验 run_opts（未知键/非法温度档/开关与档位双向矛盾响亮报错）。"""
    opts: dict = {
        "enable_temperature_rheology": False,
        "temperature_mode": "off",
        "enable_yield_gate": None,
        "rheology_formula": None,
        "mud_extrapolate": False,
        "pressure_mode": "off",
        "pressure_caliber": "shoe",
        "enable_stream_yield_gate": None,
        "include_yield_term": False,
    }
    if run_opts:
        unknown = sorted(set(run_opts) - set(opts))
        if unknown:
            raise KeyError(f"run_opts 含未知键 {unknown}，允许键：{list(opts)}")
        opts.update(run_opts)
    mode = opts["temperature_mode"]
    if mode not in TEMPERATURE_MODES:
        raise ValueError(f"temperature_mode 非法：{mode!r}，允许：{list(TEMPERATURE_MODES)}")
    if mode != "off" and not opts["enable_temperature_rheology"]:
        # 温度档被静默忽略是最危险的配错（跑出来是 T-off 却当 T-on 解读）
        raise ValueError(
            f"temperature_mode={mode!r} 但 enable_temperature_rheology=False："
            "温度档不会被消费，请改为 True 或把 mode 设为 'off'"
        )
    if mode == "off" and opts["enable_temperature_rheology"]:
        # 反向陷阱：标注 off 却打开 T 开关 ⇒ 场不注入，solver 构造层
        # temperature_rheology_t_c=60 回退恒温场（实跑 60 °C 温变流变），
        # 而汇总表会写 温度档=off / T开关=on 的矛盾标签（标签失真）
        raise ValueError(
            f"temperature_mode='off' 但 enable_temperature_rheology="
            f"{opts['enable_temperature_rheology']!r}：标注 off 实跑会回退 "
            "60 °C 恒温场（标签失真），请改 mode 或把开关设为 False"
        )
    yg = opts["enable_yield_gate"]
    if yg is not None and not isinstance(yg, bool):
        raise TypeError(f"enable_yield_gate 须为 bool|None，实际 {type(yg).__name__}")
    if not isinstance(opts["mud_extrapolate"], bool):
        raise TypeError(
            f"mud_extrapolate 须为 bool，实际 {type(opts['mud_extrapolate']).__name__}"
        )
    # R1 交叉守卫：温度关时公式参**不会被消费**（fluid_at 根本不被调用）——
    # 与 temperature_mode 的"标签失真"同型，一律响亮报错而非静默空转。
    if opts["rheology_formula"] is not None and not opts["enable_temperature_rheology"]:
        raise ValueError(
            "rheology_formula 已给但 enable_temperature_rheology=False："
            "T-off 路径不调用 fluid_at ⇒ 公式参不会被消费，请开温度开关"
        )
    params_from_spec(opts["rheology_formula"])   # 载荷合法性前置校验（错字段/空载荷响亮报错）
    syg = opts["enable_stream_yield_gate"]
    if syg is not None and not isinstance(syg, bool):
        raise TypeError(
            f"enable_stream_yield_gate 须为 bool|None，实际 {type(syg).__name__}"
        )
    # R2（Phase 2）μp/τy 真拆分：纯 bool（无 None 态——默认即 False，语义无歧义）
    if not isinstance(opts["include_yield_term"], bool):
        raise TypeError(
            f"include_yield_term 须为 bool，实际 {type(opts['include_yield_term']).__name__}"
        )
    pmode = opts["pressure_mode"]
    if pmode not in PRESSURE_MODES:
        raise ValueError(f"pressure_mode 非法：{pmode!r}，允许：{list(PRESSURE_MODES)}")
    pcal = opts["pressure_caliber"]
    if pcal not in PRESSURE_CALIBERS:
        raise ValueError(f"pressure_caliber 非法：{pcal!r}，允许：{list(PRESSURE_CALIBERS)}")
    # P-1 交叉守卫（同 temperature_mode 的"标签失真"陷阱）：压力场只在 T-on 路径被
    # 消费（`_phase_props` 的派生分支内）⇒ 温度关时给 pressure_mode 会静默空转。
    if pmode != "off" and not opts["enable_temperature_rheology"]:
        raise ValueError(
            f"pressure_mode={pmode!r} 但 enable_temperature_rheology=False："
            "压力场不会被消费（T-off 不调用 fluid_at），请开温度开关"
        )
    return opts


def casing_kwargs_from_opts(opts: dict) -> dict:
    """1D CasingFlowSolver 的 run_opts 消费（温度总开关 + R1 公式参）。

    ⚠️ **键名映射（必须在此显式改名）**：run_opts 键 ``rheology_formula`` →
    求解器形参 ``rheology_formula_params``。本函数返回的 dict 由 :func:`run_variant_res`
    ``**`` 直接展开进构造函数，键名不等即运行期 ``TypeError: unexpected keyword``。
    仓库签名闸门（``scripts/entrypoints/check_call_signatures.py``）对"函数返回的 dict"
    是盲区（实测 ``find_bad_kwargs`` 对 ``**f(opts)`` 形态返回空）⇒ 由
    ``tests/contract/test_run_opts_wiring.py`` 补闸（返回键集 ⊆ ``__init__`` 形参集）。

    ⚠️ **不得**出 ``enable_stream_yield_gate`` 键：该消费端开关只在 `AnnulusD2DGASolver`
    上存在（`CasingFlowSolver` 无）——本函数返回的 dict 同样被 ``**`` 展开进 casing 构造函数。
    ⚠️ **同样不得**出 ``include_yield_term``（R2 拆分）：`CasingFlowSolver` 无该形参，
    且 1D 套管路径根本不存在 ``fluid_apparent_viscosity`` 调用点（全仓 4 处站点全在 annulus）。
    """
    return {
        "enable_temperature_rheology": bool(opts["enable_temperature_rheology"]),
        "rheology_formula_params": params_from_spec(opts["rheology_formula"]),
        "mud_extrapolate": bool(opts["mud_extrapolate"]),
        # P-1：只映射**标量口径**；`pressure_field` 是需 (well, fluids, schedule) 的
        # **对象**，无法在只拿到 opts 的本函数里构造 ⇒ 由 `run_variant_res` 用
        # `build_pressure_field(...)` 另路合并（见该函数）。
        "pressure_caliber": str(opts["pressure_caliber"]),
    }


def annulus_kwargs_from_opts(opts: dict) -> dict:
    """2D AnnulusD2DGASolver 的 run_opts 消费（温度总开关 + R1 公式参 + 屈服门覆盖）。

    屈服门：``enable_yield_gate=None`` ⇒ **不给键**（沿用调用方既有口径——
    09-16 链的 CORRECTED_KW=True、runner 链的构造默认 True）；给 False ⇒ 出键
    覆盖为 False（Q10 乙通道代理）。调用方须先把本返回值与自己的基线 kw 合并
    成一个 dict 再 ``**`` 展开（两处都含同名键会触发 multiple-values TypeError）。

    ⚠️ 键名映射同 :func:`casing_kwargs_from_opts`（``rheology_formula`` →
    ``rheology_formula_params``）。
    """
    kw: dict = {
        "enable_temperature_rheology": bool(opts["enable_temperature_rheology"]),
        "rheology_formula_params": params_from_spec(opts["rheology_formula"]),
        "mud_extrapolate": bool(opts["mud_extrapolate"]),
        "pressure_caliber": str(opts["pressure_caliber"]),   # pressure_field 由调用方合并
    }
    if opts["enable_yield_gate"] is not None:
        kw["enable_yield_gate"] = bool(opts["enable_yield_gate"])
    # Phase 1.5（2026-10-06）：消费端开关（wall 进不进流函数算子）。
    # None ⇒ **不给键**（沿用 annulus 构造默认 False ⇒ 逐位=现状）；给值 ⇒ 显式出键。
    # ⚠️ 只有 **annulus** 有本开关；casing 侧**不得**出此键（出键=运行期 TypeError，
    # 因为 casing kwargs 也被 `**` 展开进 `CasingFlowSolver(...)`）。
    if opts["enable_stream_yield_gate"] is not None:
        kw["enable_stream_yield_gate"] = bool(opts["enable_stream_yield_gate"])
    # R2（Phase 2，2026-10-06）：μp/τy 真拆分。**仅 True 时出键**——沿用上面
    # `enable_stream_yield_gate` 的"默认值不产生额外 kwarg"模式（既有契约测试
    # `test_default_and_t_off_opts_map_to_legacy_four_key_caliber` 的 docstring 不变量的字面口径）。
    # False/缺键 ⇒ 不出键 ⇒ solver 用构造默认 False ⇒ 与 HEAD 逐位（关 2 红线）。
    if opts["include_yield_term"]:
        kw["include_yield_term"] = True
    return kw


# 交付两表进程内缓存：免 15 个 table 变体反复解析 xlsx；
# ⚠️ 共享实例 ⇒ 每次取出必须 reset_audit，否则 oob_count 跨变体累积。
_TABLE_PAIR_CACHE: dict = {}


def build_temperature_fields(
    well_key: str, mode: str
) -> tuple[ConstantTemperatureField | GeothermalTemperatureField
           | AnchoredProfileField | None,
           ConstantTemperatureField | GeothermalTemperatureField
           | AnchoredProfileField | None,
           str]:
    """按 (井, 温度档) 造 (1D 场, 2D 场, 备注)。

    - ``off``   → (None, None, "")：不注入，solver 走 T-off 逐位路径
    - ``const60`` → (Constant(60), Constant(60), "")：显式注入恒温场（与
      ``temperature_rheology_t_c`` 默认 60 的回退同值；**显式注入**意图更清楚，
      避免依赖"忘了传就是 60"的隐式默认）
    - ``static`` → (Geothermal, Geothermal, 备注)：统一地温式 T(z)=16.006+1.7598e-2·z；
      呼101/呼103 无瞬态表 ⇒ 备注打「无瞬态表」
    - ``table`` → (T_in, T_out, "")：仅呼1-004（load_delivered_pair）
    - ``anchored`` → (Anchored, Anchored, 备注)：**4d 七井静温锚定场**（注册表
      ``ANCHORED_WELLS``，k 与锚点集见 Phase 4d spec §1）；无锚井
      （``ANCHORED_NO_ANCHOR_WELLS``）⇒ Geothermal 回退 + 备注（口径声明，不硬造锚）
    - ``table_ext`` → (T_in, T_out, 备注)：仅呼1-004 的**时程扩展表**
      （333×362，表末 21600 s；`load_extended_pair_4d`）——r0.6/0.8 解禁档专用
    """
    if mode == "off":
        return None, None, ""
    if mode == "const60":
        return ConstantTemperatureField(60.0), ConstantTemperatureField(60.0), ""
    if mode == "static":
        note = "无瞬态表" if well_key not in TABLE_WELLS else ""
        return GeothermalTemperatureField(), GeothermalTemperatureField(), note
    if mode == "table":
        if well_key not in TABLE_WELLS:
            raise ValueError(
                f"{well_key}: 无交付瞬态温度表，table 档仅支持 {sorted(TABLE_WELLS)}"
            )
        if "pair" not in _TABLE_PAIR_CACHE:
            _TABLE_PAIR_CACHE["pair"] = load_delivered_pair()
        t_in, t_out = _TABLE_PAIR_CACHE["pair"]
        # 共享实例：清审计 ⇒ oob_count 只统计本次 run
        t_in.reset_audit()
        t_out.reset_audit()
        return t_in, t_out, ""
    if mode == "anchored":
        return _build_anchored_pair(well_key)
    if mode == "table_ext":
        if well_key not in TABLE_EXT_WELLS:
            raise ValueError(
                f"{well_key}: 无时程扩展温度表，table_ext 档仅支持 "
                f"{sorted(TABLE_EXT_WELLS)}"
            )
        for p in (EXT4D_T_IN_XLSX, EXT4D_T_OUT_XLSX, EXT4D_TIME_AXIS_MIN):
            if not p.exists():
                raise FileNotFoundError(
                    f"扩展温度表缺件：{p}（不入库，须由 run_ext4d.m 重生成；"
                    "禁止静默降级到交付表——交付表末 11940 s 会 clamp r0.6/0.8）"
                )
        t_in, t_out = load_extended_pair_4d(
            EXT4D_T_IN_XLSX, EXT4D_T_OUT_XLSX, EXT4D_TIME_AXIS_MIN
        )
        return t_in, t_out, "扩展表(表末21600s)"
    raise ValueError(f"未知温度档 {mode!r}，允许：{list(TEMPERATURE_MODES)}")


def _build_anchored_pair(
    well_key: str,
) -> tuple[AnchoredProfileField | GeothermalTemperatureField,
           AnchoredProfileField | GeothermalTemperatureField,
           str]:
    """「anchored」档：按井造 (1D 场, 2D 场, 备注)。

    - 注册表内有该井 ⇒ 造 `AnchoredProfileField(anchors, k, regime="circulating")`，
      1D/2D **各造独立实例**（审计计数互不污染；与 table 档的共享+reset 口径不同，
      因 anchored 无跨 run 共享对象）。
    - 无锚井（`ANCHORED_NO_ANCHOR_WELLS`）⇒ **Geothermal 回退 + 备注**（口径声明，
      不硬造锚；计划 §4-4d-2）。
    - 注册表外且不在无锚名单 ⇒ 抛 ValueError（防"配了井名但静默回退"，同 table 档口径）。
    """
    if well_key in ANCHORED_NO_ANCHOR_WELLS:
        note = "无静温锚(<2 唯一 md)⇒Geothermal 回退"
        return GeothermalTemperatureField(), GeothermalTemperatureField(), note
    spec = ANCHORED_WELLS.get(well_key)
    if spec is None:
        raise ValueError(
            f"{well_key}: 无 4d 锚定场注册项（ANCHORED_WELLS 未含，且不在 "
            f"ANCHORED_NO_ANCHOR_WELLS={sorted(ANCHORED_NO_ANCHOR_WELLS)}）"
        )
    k, anchors, _basis = spec
    mk = lambda: AnchoredProfileField(  # noqa: E731
        anchors, k, regime="circulating", source=f"4d七井锚点:{well_key}"
    )
    return mk(), mk(), ""


def build_pressure_field(well, fluids, schedule, mode: str):
    """按压力档造压力场对象（P-1，2026-10-06）。

    - ``off`` → ``None``：不注入 ⇒ solver 侧 ``pressure_field is None`` ⇒ ``fluid_at``
      收到 ``P=None`` ⇒ 记 ``p_default`` 审计 ⇒ **逐位 = HEAD**（关2 红线）。
    - ``hydrostatic`` → ``HydrostaticPressureField.from_well(well, ρ̄)``，其中
      ``ρ̄ = insitu_column_density(fluids, schedule)`` = **在场相密度按设计泵注体积加权**
      （Q9/D5 口径③；口径声明与已知偏差见 `cemdisp.data.pressure_field` docstring）。
    """
    if mode == "off":
        return None
    if mode == "hydrostatic":
        rho_bar = insitu_column_density(fluids, schedule)
        return HydrostaticPressureField.from_well(well, rho_bar)
    raise ValueError(f"未知压力档 {mode!r}，允许：{list(PRESSURE_MODES)}")


def extra_metrics(
    res,
    summary: dict,
    field_1d=None,
    field_2d=None,
) -> dict:
    """敏感性层判别量（**绝不写进 solver summary**，保 T-off 字节逐位）。

    返回键（与 T2 汇总表列同名）：
    饥饿份额 / front_narrow_m / front_wide_m / interface_length_ratio /
    屈服门活化率_b加权 / 屈服门_wall占比 / stop_t_s / 温度审计_oob。
    另附非表列信息：温度流变审计（fluid_at 计数，T-on 才有）。
    """
    from cemdisp.diagnostics.internal_consistency import starved_volume_fraction
    from cemdisp.models2d.annulus_d2dga import _trapez2d  # 同 starved 的口径来源

    geom = res.geom
    wall = res.wall_field
    metrics = res.metrics
    last = metrics.iloc[-1]

    # ① 饥饿份额（c<0.5 体积占比，b 加权；口径同 2026-09-11 探针）
    hungry = starved_volume_fraction(res.cement_field, geom)

    # ② 宽/窄边前缘（metrics 末行，m，自鞋口沿环空）
    front_wide = float(last["front_wide_m"])
    front_narrow = float(last["front_narrow_m"])

    # ③ 界面长度比（Tier0 displacement_metrics；tier0 失效 ⇒ None，见 schema 空值）
    interface_ratio = None
    tier0 = summary.get("tier0_diagnostics")
    if isinstance(tier0, dict):
        dm = tier0.get("displacement_metrics")
        if isinstance(dm, dict):
            interface_ratio = dm.get("interface_length_ratio")

    # ④ 屈服门活化率双口径（主=b 加权 mean(wall)，辅=wall>0 单元占比）
    if wall is None:
        gate_b_weighted = None
        gate_wall_frac = None
    else:
        w = np.asarray(wall, dtype=float)
        b = np.asarray(geom["b"], dtype=float)
        # b 加权 mean(wall)：与 starved_volume_fraction 同一 2D 梯形口径
        gate_b_weighted = _trapez2d(b * w, geom) / max(_trapez2d(b, geom), 1e-12)
        gate_wall_frac = float(np.mean(w > 0.0))

    # ⑤ 停算时刻（metrics 末行 time_s = total_t = min(泵总+1200, stop_t)）
    stop_t = float(last["time_s"])

    # ⑥ 温度场越界 clamp 计数（1D + 2D；T-off 无场 ⇒ 0）
    oob = (int(field_1d.oob_count) if field_1d is not None else 0) + \
          (int(field_2d.oob_count) if field_2d is not None else 0)

    out = {
        "饥饿份额": hungry,
        "front_narrow_m": front_narrow,
        "front_wide_m": front_wide,
        "interface_length_ratio": interface_ratio,
        "屈服门活化率_b加权": gate_b_weighted,
        "屈服门_wall占比": gate_wall_frac,
        "stop_t_s": stop_t,
        "温度审计_oob": oob,
        # 非表列（sidecar JSON 用）：fluid_at 绝对替换审计计数，T-on 才有
        "温度流变审计": summary.get("temperature_rheology_audit"),
    }
    return out


def run_variant_res(
    loader,
    well_fn,
    fluid_fn,
    sched_fn,
    run_opts: dict | None = None,
    well_key: str = "",
):
    """单变体一次完整流水线（1D 重跑 + 环空二维），返回 ``(res, cr, schedule2, extra)``。

    2026-10-01 从 :func:`run_variant` 抽出的内核（**逐行搬移、数值路径零改动**）：
    过程观察层（T2 时程批）需要 ``res``（snapshots/metrics）与 ``cr``
    （shoe_timeline 事件时刻），而 :func:`run_variant` 只回 summary/extra。

    返回
    ----
    res : AnnulusSimulationResult   二维求解全量结果（含时程快照）
    cr  : CasingFlowResult          1D 套管解（含 shoe_timeline / cement_end_time_s）
    schedule2 : PumpingSchedule     本变体变换后的泵序（算泵注 50% 用）
    extra : dict                    敏感性层判别量（口径同 :func:`extra_metrics`）
    """
    opts = normalize_run_opts(run_opts)
    well, fluids, schedule, _ = loader()
    well2, fluids2, schedule2 = well_fn(well), fluid_fn(fluids), sched_fn(schedule)

    field_1d, field_2d, note = build_temperature_fields(
        well_key, opts["temperature_mode"]
    )
    # P-1：压力场是**对象**（需 well/fluids/schedule）⇒ 不能进 kwargs 映射函数，
    # 在此另路构造并传给两个 solver（评审 C-08：否则 pressure_mode 静默丢弃）。
    p_field = build_pressure_field(well2, fluids2, schedule2, opts["pressure_mode"])

    casing = CasingFlowSolver(
        enable_gravity=True,
        mixing_contact_time=True,
        plug_face_zero_mixing=True,
        has_plug=True,
        pressure_field=p_field,
        **casing_kwargs_from_opts(opts),
    )
    if opts["enable_temperature_rheology"] and field_1d is not None:
        # 1D 用 T_in 表 / 静温剖面（管内语义由场对象承载）
        cr = casing.run(well2, fluids2, schedule2, temperature_field=field_1d)
    else:
        cr = casing.run(well2, fluids2, schedule2)   # T-off：不读该参，零痕迹
    inlet = build_coupled_annulus_inlet_provider(
        cr, casing, fluids2, split_cement_phases=True,
    )
    tt = min(_total_t(schedule2) + 1200.0, _stop_t(cr, fluids2))
    # 先合并成单个 dict 再展开：CORRECTED_KW 与 run_opts 都可能含 enable_yield_gate
    solver_kw = {**CORRECTED_KW, **annulus_kwargs_from_opts(opts),
                 "pressure_field": p_field}
    solver = AnnulusD2DGASolver(
        total_t=tt, nz=NZ, enable_cfl_adaptive=True, **solver_kw,
    )
    if opts["enable_temperature_rheology"] and field_2d is not None:
        # 2D 用 T_out 表 / 静温剖面（环空语义由场对象承载）
        res = solver.run(
            well2, fluids2, inlet, schedule=schedule2, temperature_field=field_2d,
        )
    else:
        res = solver.run(well2, fluids2, inlet, schedule=schedule2)

    extra = extra_metrics(res, res.summary, field_1d, field_2d)
    extra["温度场备注"] = note
    return res, cr, schedule2, extra


def run_variant(
    loader,
    well_fn,
    fluid_fn,
    sched_fn,
    run_opts: dict | None = None,
    well_key: str = "",
) -> tuple[dict, dict]:
    """单变体一次完整流水线，返回 ``(summary, extra)``——:func:`run_variant_res` 的薄壳。

    - ``summary``：``res.summary`` **原样**（一个键都不加 ⇒ T-off 字节逐位红线）
    - ``extra``：敏感性层判别量（见 :func:`extra_metrics`），另附
      ``温度场备注``（如「无瞬态表」）——供下游写新 schema 汇总表。
    - ``run_opts=None`` ⇒ 纯 T-off（温度开关 False、不注入场、屈服门沿用
      CORRECTED_KW），与本文件 2026-09-16 首版行为逐位一致。
    - ``well_key`` 供 :func:`build_temperature_fields` 判定 table 档可用性
      与「无瞬态表」备注；run_opts=None（T-off）时恒不消费。
    - 需要 res/cr/schedule2 的下游（T2 时程批）直接用 :func:`run_variant_res`。
    """
    res, _cr, _schedule2, extra = run_variant_res(
        loader, well_fn, fluid_fn, sched_fn, run_opts, well_key=well_key,
    )
    return res.summary, extra


def main() -> None:
    # 冻结目录保护（红线：09-16 产物冻结零写入）。装配层/变体定义照常供
    # 09-27 两批与 T2 批 import；本入口仅作历史复现且被内容守卫拒绝。
    if any(OUT_DIR.glob("*_结果摘要.json")):
        raise RuntimeError(
            f"输出目录已冻结（存在既有产物）：{OUT_DIR}\n"
            "09-16 批为冻结基线，禁止重跑覆盖；新批次请用各自的新目录"
            "（09-27 / 温压T2）。"
        )
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    variants = build_variants()
    rows: list[dict] = []

    for well_key, (mod_name, fn_name, baseline_file) in WELLS.items():
        loader = getattr(importlib.import_module(mod_name), fn_name)
        base = json.loads((BASELINE_DIR / baseline_file).read_text(encoding="utf-8"))
        base_eta_e, base_eta_n = float(base["eta_E"]), float(base["eta_N"])
        print(f"\n=== {well_key}  基线（当前口径） η_E={base_eta_e:.4f} η_N={base_eta_n:.4f} ===",
              flush=True)

        for var_name, well_fn, fluid_fn, sched_fn, run_opts in variants:
            case_json = OUT_DIR / f"{well_key}_{var_name}_结果摘要.json"
            if case_json.exists():                       # 断点续跑
                summary = json.loads(case_json.read_text(encoding="utf-8"))
                final, elapsed = summary["最终结果"], None
                print(f"  [复用] {well_key} × {var_name}", flush=True)
            else:
                t0 = time.perf_counter()
                summary, _extra = run_variant(
                    loader, well_fn, fluid_fn, sched_fn, run_opts, well_key=well_key,
                )
                final = summary["最终结果"]
                elapsed = round(time.perf_counter() - t0, 1)
                case_json.write_text(
                    json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8",
                )
                print(f"  [计算] {well_key} × {var_name}: "
                      f"η_E={float(final['全井段最终有效顶替效率']):.4f} "
                      f"η_N={float(final['窄四分位效率']):.4f} ({elapsed}s)", flush=True)

            eta_e = float(final["全井段最终有效顶替效率"])
            eta_n = float(final["窄四分位效率"])
            rows.append({
                "井名": well_key,
                "变体": var_name,
                "η_E": eta_e,
                "η_N": eta_n,
                "窜槽": float(final["最终窜槽指数"]),
                "混浆": float(final["最终混浆指数"]),
                "失稳": float(final["最终失稳指数"]),
                "Δη_E_pp": (eta_e - base_eta_e) * 100.0,
                "Δη_N_pp": (eta_n - base_eta_n) * 100.0,
                "基线η_E": base_eta_e,
                "基线η_N": base_eta_n,
                "耗时_s": elapsed if elapsed is not None else "",
            })

    csv_path = OUT_DIR / "汇总表_敏感性变体_当前口径.csv"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    md = [
        "# 敏感性变体汇总（当前口径 2026-09-16，流函数新路径）",
        "",
        f"口径：nz={NZ}、CFL 自适应 on、CORRECTED_KW、T1 套管生产开关；"
        "与 `results/源模型口径重跑_2026-09-14` 严格同源。",
        "基线：`results/源模型口径重跑_2026-09-14/单井结果/<well>_corrected_on.json`。",
        "",
        "| 井名 | 变体 | η_E | η_N | Δη_E/pp | Δη_N/pp | 耗时/s |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        md.append(f"| {r['井名']} | {r['变体']} | {r['η_E']:.4f} | {r['η_N']:.4f} "
                  f"| {r['Δη_E_pp']:+.2f} | {r['Δη_N_pp']:+.2f} | {r['耗时_s']} |")
    (OUT_DIR / "汇总表_敏感性变体_当前口径.md").write_text("\n".join(md), encoding="utf-8")

    print(f"\n完成：{len(rows)} 个变体 → {OUT_DIR}", flush=True)


if __name__ == "__main__":
    main()
