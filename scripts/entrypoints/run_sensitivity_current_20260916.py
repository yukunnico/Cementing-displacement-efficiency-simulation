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
from dataclasses import replace
from pathlib import Path

import numpy as np

from cemdisp.data.fluid_spec import FluidRole, FluidSpec
from cemdisp.data.pumping_schedule import PumpingSchedule
from cemdisp.data.temperature_field import (
    ConstantTemperatureField,
    GeothermalTemperatureField,
    load_delivered_pair,
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
# run_opts 三键（T2 施工口径 2026-10-01）：
#   enable_temperature_rheology : 总开关（False ⇒ 数值路径与 HEAD 逐位）
#   temperature_mode            : "off" | "static" | "table" | "const60"
#   enable_yield_gate           : None=沿用 CORRECTED_KW（True）/ False=Q10乙通道代理
RUN_OPTS_KEYS = ("enable_temperature_rheology", "temperature_mode", "enable_yield_gate")
TEMPERATURE_MODES = ("off", "static", "table", "const60")
# 有交付瞬态温度表的井（table 档专用；呼101/呼103 无表 ⇒ 只能 static/const60/off）
TABLE_WELLS = frozenset({"呼1-004"})


def normalize_run_opts(run_opts: dict | None) -> dict:
    """补默认 + 校验 run_opts（未知键/非法温度档/开关与档位双向矛盾响亮报错）。"""
    opts: dict = {
        "enable_temperature_rheology": False,
        "temperature_mode": "off",
        "enable_yield_gate": None,
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
    return opts


def casing_kwargs_from_opts(opts: dict) -> dict:
    """1D CasingFlowSolver 的 run_opts 消费（仅温度总开关）。"""
    return {"enable_temperature_rheology": bool(opts["enable_temperature_rheology"])}


def annulus_kwargs_from_opts(opts: dict) -> dict:
    """2D AnnulusD2DGASolver 的 run_opts 消费（温度总开关 + 屈服门覆盖）。

    屈服门：``enable_yield_gate=None`` ⇒ **不给键**（沿用调用方既有口径——
    09-16 链的 CORRECTED_KW=True、runner 链的构造默认 True）；给 False ⇒ 出键
    覆盖为 False（Q10 乙通道代理）。调用方须先把本返回值与自己的基线 kw 合并
    成一个 dict 再 ``**`` 展开（两处都含同名键会触发 multiple-values TypeError）。
    """
    kw: dict = {
        "enable_temperature_rheology": bool(opts["enable_temperature_rheology"]),
    }
    if opts["enable_yield_gate"] is not None:
        kw["enable_yield_gate"] = bool(opts["enable_yield_gate"])
    return kw


# 交付两表进程内缓存：免 15 个 table 变体反复解析 xlsx；
# ⚠️ 共享实例 ⇒ 每次取出必须 reset_audit，否则 oob_count 跨变体累积。
_TABLE_PAIR_CACHE: dict = {}


def build_temperature_fields(
    well_key: str, mode: str
) -> tuple[ConstantTemperatureField | GeothermalTemperatureField | None,
           ConstantTemperatureField | GeothermalTemperatureField | None,
           str]:
    """按 (井, 温度档) 造 (1D 场, 2D 场, 备注)。

    - ``off``   → (None, None, "")：不注入，solver 走 T-off 逐位路径
    - ``const60`` → (Constant(60), Constant(60), "")：显式注入恒温场（与
      ``temperature_rheology_t_c`` 默认 60 的回退同值；**显式注入**意图更清楚，
      避免依赖"忘了传就是 60"的隐式默认）
    - ``static`` → (Geothermal, Geothermal, 备注)：统一地温式 T(z)=16.006+1.7598e-2·z；
      呼101/呼103 无瞬态表 ⇒ 备注打「无瞬态表」
    - ``table`` → (T_in, T_out, "")：仅呼1-004（load_delivered_pair）
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
    raise ValueError(f"未知温度档 {mode!r}，允许：{list(TEMPERATURE_MODES)}")


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


def run_variant(
    loader,
    well_fn,
    fluid_fn,
    sched_fn,
    run_opts: dict | None = None,
    well_key: str = "",
) -> tuple[dict, dict]:
    """单变体一次完整流水线（1D 重跑 + 环空二维），返回 ``(summary, extra)``。

    - ``summary``：``res.summary`` **原样**（一个键都不加 ⇒ T-off 字节逐位红线）
    - ``extra``：敏感性层判别量（见 :func:`extra_metrics`），另附
      ``温度场备注``（如「无瞬态表」）——供下游写新 schema 汇总表。
    - ``run_opts=None`` ⇒ 纯 T-off（温度开关 False、不注入场、屈服门沿用
      CORRECTED_KW），与本文件 2026-09-16 首版行为逐位一致。
    - ``well_key`` 供 :func:`build_temperature_fields` 判定 table 档可用性
      与「无瞬态表」备注；run_opts=None（T-off）时恒不消费。
    """
    opts = normalize_run_opts(run_opts)
    well, fluids, schedule, _ = loader()
    well2, fluids2, schedule2 = well_fn(well), fluid_fn(fluids), sched_fn(schedule)

    field_1d, field_2d, note = build_temperature_fields(
        well_key, opts["temperature_mode"]
    )

    casing = CasingFlowSolver(
        enable_gravity=True,
        mixing_contact_time=True,
        plug_face_zero_mixing=True,
        has_plug=True,
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
    solver_kw = {**CORRECTED_KW, **annulus_kwargs_from_opts(opts)}
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
