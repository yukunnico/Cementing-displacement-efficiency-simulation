"""T2 温压敏感性批（2026-10-01）——温度开关/温度档/排量/屈服门联合敏感性。

用户裁定（2026-10-01，全部生效，本脚本只照单执行）
--------------------------------------------------
1. **链口径 = 09-16 脚本链**：CORRECTED_KW + CFL 自适应、
   ``tt = min(_total_t + 1200, _stop_t)``；**不跑 runner 链**。
2. **输出新目录** ``results/敏感性变体_温压T2_2026-10-01/``（三旧敏感性目录冻结
   零写入，本脚本有目录守卫）；CSV 走新 schema；**基线 = 本批 ``Toff_zero``**；
   断点续跑（同名 JSON 复用）。
3. **变体全新命名**（``Toff_*`` / ``Ton_*`` 前缀），绝不复用 ``rate_x0.6`` 等旧名。
4. **Q10 乙 = 通道代理**：μp/τy 不拆公式层；τy 通道用
   ``enable_yield_gate=False`` 对照行（run_opts 传入）。
5. **静温剖面 = 统一式** ``T(z)=16.006+1.7598e-2·z``（°C，z 单位 m），
   呼101/呼103 用它并打「无瞬态表」标注（``GeothermalTemperatureField`` 单一真源）。
6. **判别量写敏感性层新表，绝不加进 solver summary**（保 T-off 字节逐位）。
   屈服门活化率主口径 = 末态 b 加权 mean(wall)，辅口径 = wall>0 单元占比，双列。
7. **绝对替换不动**（红线：不回缩放、不拆 μp/τy、公式层零改动）；
   ``enable_hb_closure`` / ``hb_fix_cement_tau_y`` / ``enable_stream_yield_gate``
   **全程不开**（本脚本 run_opts 只有三键，压根传不进这三个开关）。

变体矩阵（45 个）
-----------------
- 呼1-004（23）：``Toff_zero``；``Toff_rate_x{0.6,0.8,1.2,1.4}``；
  ``Ton_{static,const60,table}_rate_x{0.6,0.8,1.0,1.2,1.4}``（15）；
  ``Ton_table_gateoff_rate_x1.0``（Q10 乙）；
  ``Ton_static_cement_n_p0.1_rate_x1.0`` / ``Ton_static_mud_pv_x1.3_rate_x1.0``
  （Q13 恒等检验，预期 Δ≈0）。
- 呼101、呼103 各 11：``Toff_zero``；``Toff_rate_x{0.6,1.4}``；
  ``Ton_{static,const60}_rate_x{0.6,1.0,1.4}``（6，static 备注「无瞬态表」）；
  ``Ton_static_gateoff_rate_x1.0``；``Ton_static_cement_n_p0.1_rate_x1.0``。

双批次（2026-10-06 执行窗口裁定，用户选 A 案）
-----------------------------------------------
--batch t2（默认）：2026-10-01 原 45 变体（呼1-004 23 + 呼101 11 + 呼103 11），
**逐位不变**——不加 --out-dir 时行为与本脚本 2026-10-01 版完全相同。
--batch phase0：2026-10-06 三重点井补跑矩阵 54 变体（呼101 14 + 呼1-003 20 +
呼1-004 20），**必须**配 --out-dir 指向新日期目录（红线：补跑新目录带日期后缀）。
规格 = docs/superpowers/specs/2026-10-06-phase0-supplementary-runs-design.md。
呼103 **不入 phase0 批**（2026-10-06 裁定：退出重点井，历史资产照常引用）。
前置 = Phase 0.0 密度「就近取」口径（Q16，2026-10-06 用户裁定）。

输出
----
results/敏感性变体_温压T2_2026-10-01/        （batch=t2 默认）
    {井}_{变体}_结果摘要.json   res.summary 原样（schema 同 09-16 批，一个键不加）
    {井}_{变体}_判别量.json     敏感性层判别量 + 耗时_s（断点续跑据此复用）
    汇总表_温压T2.csv/.md       新 schema 主表
    分解表_温压T2.csv/.md       同排量两分量：口径差=const60−Toff、温度场效应=mode−const60

用法
----
    PYTHONIOENCODING=utf-8 PYTHONUTF8=1 conda run -n cementT \
        python scripts/entrypoints/run_sensitivity_temperature_t2_20261001.py
    # 子集（冒烟/补跑）：
    #   --well 呼1-004 --variant Ton_const60_rate_x1.0
    # Phase 0 三重点井补跑批：
    #   --batch phase0 --out-dir results/敏感性补跑_三重点井_20261006

断点续跑与强制重算
------------------
同名 ``*_判别量.json`` 存在且 ``schema_version`` 相符 ⇒ 复用不重跑。
**改码后须 ``--force`` 或删除 ``*_判别量.json`` 才会重算**（判别量由代码
公式算出，改码不改档名不会触发重跑，必须显式强制）。``--force`` 时忽略
已存在 JSON、全部重算并覆盖。
"""
from __future__ import annotations

import argparse
import csv
import importlib
import json
import re
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
_SCRIPTS_DIR = PROJECT_ROOT / "scripts"
for p in (str(PROJECT_ROOT), str(_SCRIPTS_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

# 变体变换与流水线与 09-16 批同源（scale_schedule/shift_cement_n/scale_mud_pv 复用）
from entrypoints.run_sensitivity_current_20260916 import (  # noqa: E402
    _identity,
    run_variant,
    scale_mud_pv,
    scale_schedule,
    shift_cement_n,
    shift_spacer_density,
    shift_standoff,
)

OUT_DIR = PROJECT_ROOT / "results" / "敏感性变体_温压T2_2026-10-01"

# Phase 0 补跑批输出目录（2026-10-06 执行窗口；红线：补跑新目录带日期后缀）
PHASE0_DIR = PROJECT_ROOT / "results" / "敏感性补跑_三重点井_20261006"

# 批次定义：out_dir=默认输出目录，table_stem=汇总/分解表文件名主干，title=表头标题。
# **默认批次 t2 的一切取值与 2026-10-01 原版逐字相同**（关2：默认路径行为不变）。
BATCHES: dict[str, dict] = {
    "t2": {
        "out_dir": OUT_DIR,
        "table_stem": "温压T2",
        "title": "T2 温压敏感性",
        "date": "2026-10-01",
    },
    "phase0": {
        "out_dir": PHASE0_DIR,
        "table_stem": "三重点井补跑",
        "title": "三重点井温压补跑",
        "date": "2026-10-06",
    },
}

# 批次变体数期望值（响亮失败，防手滑改坏矩阵；键即井名）
EXPECTED_VARIANTS: dict[str, dict[str, int]] = {
    "t2": {"呼1-004": 23, "呼101": 11, "呼103": 11},
    "phase0": {"呼101": 14, "呼1-003": 20, "呼1-004": 20},
}

# 三旧敏感性目录（冻结产物，红线：零写入）
FROZEN_DIRS = (
    PROJECT_ROOT / "results" / "敏感性变体_当前口径_2026-09-16",
    PROJECT_ROOT / "results" / "敏感性变体_呼1-004_2026-09-27",
    PROJECT_ROOT / "results" / "敏感性变体_runner链_2026-09-27",
)
for _d in FROZEN_DIRS:
    if OUT_DIR == _d or _d in OUT_DIR.parents:
        raise RuntimeError(f"输出目录防护触发：{OUT_DIR} 落在冻结目录 {_d} 内")

# 井名 → (loader 模块, loader 函数)。顺序即跑批顺序（呼1-004 全量在前）。
WELLS: dict[str, tuple[str, str]] = {
    "呼1-004": ("cemdisp.data.loaders.ht1_004_loader", "load_ht1_004_tailpipe"),
    "呼101": ("cemdisp.data.loaders.hu101_loader", "load_hu101_tailpipe"),
    "呼103": ("cemdisp.data.loaders.hu103_loader", "load_hu103_tailpipe"),
    # 2026-10-06：新增重点井呼1-003（static 档 only；无瞬态表 ⇒ 备注「无瞬态表」）
    "呼1-003": ("cemdisp.data.loaders.ht1_003_loader", "load_ht1_003_tailpipe"),
}

BASELINE_VARIANT = "Toff_zero"

# 判别量 JSON 的 schema 版本（写入每个 *_判别量.json；_load_case 要求版本
# 相符才复用 ⇒ 缺版本号或不符一律重算，防改码后吃到旧公式产物）
EXTRA_SCHEMA_VERSION = 2

# 汇总表 schema（列序即写盘序）
CSV_COLUMNS = [
    "井名", "变体", "温度档", "T开关", "yield_gate", "备注",
    "η_E", "η_N", "饥饿份额",
    "front_narrow_m", "front_wide_m", "interface_length_ratio",
    "屈服门活化率_b加权", "屈服门_wall占比",
    "stop_t_s", "温度审计_oob",
    "Δη_E_pp", "Δη_N_pp", "基线η_E", "基线η_N", "耗时_s",
]

# 分解表 schema（脚本层算术，从行数据聚合；差值列以 Δ_ 前缀）
DECOMP_METRICS = ["η_E", "η_N", "饥饿份额", "front_narrow_m", "front_wide_m",
                  "interface_length_ratio", "屈服门活化率_b加权", "屈服门_wall占比"]
DECOMP_COLUMNS = ["井名", "排量档", "分量", "对照"] + [f"Δ_{m}" for m in DECOMP_METRICS]

# 变体级备注（与建表时的温度场备注「无瞬态表」合并）
VARIANT_REMARKS = {
    "Ton_table_gateoff_rate_x1.0": "Q10乙通道代理（屈服门关）",
    "Ton_static_gateoff_rate_x1.0": "Q10乙通道代理（屈服门关）",
    "Ton_static_cement_n_p0.1_rate_x1.0": "Q13恒等检验（预期Δ≈0）",
    "Ton_static_mud_pv_x1.3_rate_x1.0": "Q13恒等检验（预期Δ≈0）",
}

# 变体名 → 温度档/排量 的解析（分解表用；只认纯温度×排量变体）
_RE_TOFF_RATE = re.compile(r"^Toff_rate_x([0-9.]+)$")
_RE_TON_RATE = re.compile(r"^Ton_(static|const60|table)_rate_x([0-9.]+)$")


# ---------------------------------------------------------------- 变体注册表
def _rate(factor: float):
    """排量缩放变换（`scale_schedule` 同一实现；1.0 恒等）。"""
    return lambda s: scale_schedule(s, factor)


def _opts(mode: str, *, on: bool = True, yield_gate: bool | None = None) -> dict:
    """run_opts（09-16 共享装配层三键）。"""
    return {
        "enable_temperature_rheology": on,
        "temperature_mode": mode,
        "enable_yield_gate": yield_gate,
    }


# standoff 9 档（与 09-16 矩阵同名档位，可跨批对照）：±0.05 / ±0.10 / −0.15 / ±0.20 / ±0.30
#
# ⚠️ 退化声明（Phase 0 预检发现）：呼1-003 / 呼1-004 的 standoff 剖面恒 0.83，
# 而 shift_standoff 用 np.clip(v+δ, 0, 1) ⇒ **+0.20 与 +0.30 都被 clip 成 1.0，
# 两变体剖面指纹完全相同**。本注册表按计划保留完整 9 档以维持三井可比，
# 但产物/分析必须声明「27 名义格 = 25 个互异剖面」，不得把这两档当独立点。
STANDOFF_DELTAS: tuple[tuple[str, float], ...] = (
    ("m0.30", -0.30), ("m0.20", -0.20), ("m0.15", -0.15), ("m0.10", -0.10),
    ("m0.05", -0.05), ("p0.05", +0.05), ("p0.10", +0.10), ("p0.20", +0.20),
    ("p0.30", +0.30),
)
SPACER_DENS_DELTAS: tuple[tuple[str, float], ...] = (("p100", +100.0), ("m100", -100.0))


def _standoff_variants() -> list[tuple]:
    """standoff 9 档（T-on static × rate1.0），变体名与 09-16 档位一一对应。"""
    return [
        (f"Ton_static_standoff_{tag}_rate_x1.0",
         (lambda w, _d=d: shift_standoff(w, _d)), _identity, _rate(1.0), _opts("static"))
        for tag, d in STANDOFF_DELTAS
    ]


def _spacer_dens_variants() -> list[tuple]:
    """隔离液密度 ±100 kg/m³（T-on static × rate1.0）。"""
    return [
        (f"Ton_static_spacer_dens_{tag}_rate_x1.0",
         _identity, (lambda fs, _d=d: shift_spacer_density(fs, _d)), _rate(1.0),
         _opts("static"))
        for tag, d in SPACER_DENS_DELTAS
    ]


def _phase0_variants() -> dict[str, list[tuple]]:
    """Phase 0 三重点井补跑矩阵（2026-10-06 执行窗口裁定）。"
    """
    reg: dict[str, list[tuple]] = {}

    # ---- 呼101（14）：判别井；T-on static 主干 + standoff 9 档 + 密度 ±100 ----
    reg["呼101"] = [
        ("Toff_zero", _identity, _identity, _identity, _opts("off", on=False)),
        ("Ton_const60_rate_x1.0", _identity, _identity, _rate(1.0), _opts("const60")),
        ("Ton_static_rate_x1.0", _identity, _identity, _rate(1.0), _opts("static")),
    ] + _standoff_variants() + _spacer_dens_variants()

    # ---- 呼1-003（20）：首建；static only（无瞬态表）+ 排量三点 + standoff + 密度 ----
    v: list[tuple] = [
        ("Toff_zero", _identity, _identity, _identity, _opts("off", on=False)),
    ]
    # 排量非 1.0 档补同排量 T-off 对照（红线：T-on 数字必附口径差分解；
    # 分解表的「口径差」需 Ton_const60_rate_x<r> − Toff_rate_x<r> 配对）
    for r in (0.6, 1.4):
        v.append((f"Toff_rate_x{r}", _identity, _identity, _rate(r), _opts("off", on=False)))
    for r in (0.6, 1.0, 1.4):
        v.append((f"Ton_const60_rate_x{r}", _identity, _identity, _rate(r), _opts("const60")))
    for r in (0.6, 1.0, 1.4):
        v.append((f"Ton_static_rate_x{r}", _identity, _identity, _rate(r), _opts("static")))
    reg["呼1-003"] = v + _standoff_variants() + _spacer_dens_variants()

    # ---- 呼1-004（20）：饱和井；主干 + standoff + 密度 + 排量加密 x0.8/x1.2 ----
    v = [
        ("Toff_zero", _identity, _identity, _identity, _opts("off", on=False)),
        ("Ton_const60_rate_x1.0", _identity, _identity, _rate(1.0), _opts("const60")),
        ("Ton_static_rate_x1.0", _identity, _identity, _rate(1.0), _opts("static")),
    ]
    for r in (0.8, 1.2):
        v.append((f"Toff_rate_x{r}", _identity, _identity, _rate(r), _opts("off", on=False)))
    for r in (0.8, 1.2):
        v.append((f"Ton_const60_rate_x{r}", _identity, _identity, _rate(r), _opts("const60")))
    for r in (0.8, 1.2):
        v.append((f"Ton_static_rate_x{r}", _identity, _identity, _rate(r), _opts("static")))
    reg["呼1-004"] = v + _standoff_variants() + _spacer_dens_variants()

    return reg


def build_temperature_variants(batch: str = "t2") -> dict[str, list[tuple]]:
    """按批次返回注册表：井名 → [(变体名, well_fn, fluid_fn, sched_fn, run_opts), ...]。

    ``batch="t2"``（默认）= 2026-10-01 原 45 变体，**逐位不变**；
    ``batch="phase0"`` = 2026-10-06 三重点井补跑矩阵（52 变体）。
    与 09-16 legacy `build_variants()` 完全分离（命名 Toff_/Ton_ 前缀，绝不复用）。
    """
    if batch == "phase0":
        return _phase0_variants()
    if batch != "t2":
        raise ValueError(f"未知批次 {batch!r}，允许：{sorted(BATCHES)}")
    reg: dict[str, list[tuple]] = {}

    # ---- 呼1-004（23）：有交付瞬态表 ⇒ 含 table 档全排量扫描 ----
    v: list[tuple] = [
        ("Toff_zero", _identity, _identity, _identity, _opts("off", on=False)),
    ]
    for r in (0.6, 0.8, 1.2, 1.4):
        v.append((f"Toff_rate_x{r}", _identity, _identity, _rate(r),
                  _opts("off", on=False)))
    for mode in ("static", "const60", "table"):
        for r in (0.6, 0.8, 1.0, 1.2, 1.4):
            v.append((f"Ton_{mode}_rate_x{r}", _identity, _identity, _rate(r),
                      _opts(mode)))
    v.append(("Ton_table_gateoff_rate_x1.0", _identity, _identity, _rate(1.0),
              _opts("table", yield_gate=False)))
    v.append(("Ton_static_cement_n_p0.1_rate_x1.0",
              _identity, lambda fs: shift_cement_n(fs, +0.1), _rate(1.0),
              _opts("static")))
    v.append(("Ton_static_mud_pv_x1.3_rate_x1.0",
              _identity, lambda fs: scale_mud_pv(fs, 1.3), _rate(1.0),
              _opts("static")))
    reg["呼1-004"] = v

    # ---- 呼101 / 呼103 各 11：无瞬态表 ⇒ 只有 static/const60 ----
    for well in ("呼101", "呼103"):
        v = [
            ("Toff_zero", _identity, _identity, _identity, _opts("off", on=False)),
        ]
        for r in (0.6, 1.4):
            v.append((f"Toff_rate_x{r}", _identity, _identity, _rate(r),
                      _opts("off", on=False)))
        for mode in ("static", "const60"):
            for r in (0.6, 1.0, 1.4):
                v.append((f"Ton_{mode}_rate_x{r}", _identity, _identity, _rate(r),
                          _opts(mode)))
        v.append(("Ton_static_gateoff_rate_x1.0", _identity, _identity, _rate(1.0),
                  _opts("static", yield_gate=False)))
        v.append(("Ton_static_cement_n_p0.1_rate_x1.0",
                  _identity, lambda fs: shift_cement_n(fs, +0.1), _rate(1.0),
                  _opts("static")))
        reg[well] = v

    return reg


# ---------------------------------------------------------------- 落盘/复用
def _summary_path(well: str, variant: str, out_dir: Path = OUT_DIR) -> Path:
    return out_dir / f"{well}_{variant}_结果摘要.json"


def _extra_path(well: str, variant: str, out_dir: Path = OUT_DIR) -> Path:
    return out_dir / f"{well}_{variant}_判别量.json"


def _jsonable(obj):
    """numpy 标量 → Python 原生（json 无法序列化 np.int64）。"""
    if isinstance(obj, dict):
        return {k: _jsonable(x) for k, x in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(x) for x in obj]
    if hasattr(obj, "item") and not isinstance(obj, (str, bytes)):
        try:
            return obj.item()
        except Exception:  # noqa: BLE001 —— 非标量数组等原样交给 json 报错
            return obj
    return obj


def _load_case(well: str, variant: str, out_dir: Path = OUT_DIR) -> tuple[dict, dict] | None:
    """summary + 判别量 JSON 齐在**且 schema_version 相符**才视为可复用。

    缺任一文件、缺 ``schema_version`` 键或版本号不符 ⇒ 返回 None（须重跑）。
    改码后判别量公式变了但档名没变时，靠 bump :data:`EXTRA_SCHEMA_VERSION`
    或 ``--force`` / 删除 ``*_判别量.json`` 触发重算。
    """
    sj, ej = _summary_path(well, variant, out_dir), _extra_path(well, variant, out_dir)
    if not (sj.exists() and ej.exists()):
        return None
    summary = json.loads(sj.read_text(encoding="utf-8"))
    extra = json.loads(ej.read_text(encoding="utf-8"))
    if extra.get("schema_version") != EXTRA_SCHEMA_VERSION:
        return None
    return summary, extra


def _baseline_of(well: str, out_dir: Path = OUT_DIR) -> dict | None:
    """本批 Toff_zero 基线（未跑则 None ⇒ Δ 列留空）。"""
    case = _load_case(well, BASELINE_VARIANT, out_dir)
    if case is None:
        return None
    final = case[0]["最终结果"]
    return {
        "eta_E": float(final["全井段最终有效顶替效率"]),
        "eta_N": float(final["窄四分位效率"]),
    }


def _run_or_reuse(well: str, loader, variant: str, transforms, run_opts,
                  force: bool = False, out_dir: Path = OUT_DIR) -> tuple[dict, dict, object]:
    """跑（或复用）一个变体，返回 (summary, extra, elapsed_s)。

    ``force=True``（CLI ``--force``）⇒ 忽略已存在 JSON，全部重算并覆盖。
    """
    well_fn, fluid_fn, sched_fn = transforms
    case = None if force else _load_case(well, variant, out_dir)
    if case is not None:
        summary, extra = case
        print(f"  [复用] {well} × {variant}", flush=True)
        return summary, extra, extra.get("耗时_s", "")

    t0 = time.perf_counter()
    summary, extra = run_variant(
        loader, well_fn, fluid_fn, sched_fn, run_opts, well_key=well,
    )
    elapsed = round(time.perf_counter() - t0, 1)
    extra = dict(extra)
    extra["耗时_s"] = elapsed
    extra["schema_version"] = EXTRA_SCHEMA_VERSION
    _summary_path(well, variant, out_dir).write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8",
    )
    _extra_path(well, variant, out_dir).write_text(
        json.dumps(_jsonable(extra), ensure_ascii=False, indent=2), encoding="utf-8",
    )
    final = summary["最终结果"]
    print(f"  [计算] {well} × {variant}: "
          f"η_E={float(final['全井段最终有效顶替效率']):.4f} "
          f"η_N={float(final['窄四分位效率']):.4f} ({elapsed}s)", flush=True)
    return summary, extra, elapsed


# ---------------------------------------------------------------- 汇总表
def _row_of(well: str, variant: str, run_opts: dict, summary: dict,
            extra: dict, baseline: dict | None) -> dict:
    final = summary["最终结果"]
    eta_e = float(final["全井段最终有效顶替效率"])
    eta_n = float(final["窄四分位效率"])
    mode = run_opts["temperature_mode"]
    yg = run_opts["enable_yield_gate"]
    notes = [x for x in (
        VARIANT_REMARKS.get(variant, ""),
        str(extra.get("温度场备注") or ""),
    ) if x]
    return {
        "井名": well,
        "变体": variant,
        "温度档": mode,
        "T开关": "on" if run_opts["enable_temperature_rheology"] else "off",
        "yield_gate": "off" if yg is False else "on",
        "备注": "；".join(notes),
        "η_E": eta_e,
        "η_N": eta_n,
        "饥饿份额": extra.get("饥饿份额"),
        "front_narrow_m": extra.get("front_narrow_m"),
        "front_wide_m": extra.get("front_wide_m"),
        "interface_length_ratio": extra.get("interface_length_ratio"),
        "屈服门活化率_b加权": extra.get("屈服门活化率_b加权"),
        "屈服门_wall占比": extra.get("屈服门_wall占比"),
        "stop_t_s": extra.get("stop_t_s"),
        "温度审计_oob": extra.get("温度审计_oob"),
        "Δη_E_pp": "" if baseline is None else (eta_e - baseline["eta_E"]) * 100.0,
        "Δη_N_pp": "" if baseline is None else (eta_n - baseline["eta_N"]) * 100.0,
        "基线η_E": "" if baseline is None else baseline["eta_E"],
        "基线η_N": "" if baseline is None else baseline["eta_N"],
        "耗时_s": extra.get("耗时_s", ""),
    }


def rebuild_tables(registry: dict[str, list[tuple]], out_dir: Path = OUT_DIR,
                   table_stem: str = "温压T2", title: str = "T2 温压敏感性",
                   date: str = "2026-10-01") -> None:
    """从落盘 JSON 幂等重建汇总表 + 分解表（缺档跳过 ⇒ 支持部分批）。"""
    rows: list[dict] = []
    metrics_by_key: dict[tuple[str, str], dict] = {}
    for well, variants in registry.items():
        baseline = _baseline_of(well, out_dir)
        if baseline is None:
            print(f"  [提示] {well}: {BASELINE_VARIANT} 未跑，Δ 列留空", flush=True)
        for name, *_rest, run_opts in variants:
            case = _load_case(well, name, out_dir)
            if case is None:
                continue
            summary, extra = case
            row = _row_of(well, name, run_opts, summary, extra, baseline)
            rows.append(row)
            metrics_by_key[(well, name)] = row
    if rows:
        with (out_dir / f"汇总表_{table_stem}.csv").open(
            "w", encoding="utf-8-sig", newline=""
        ) as fh:
            writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
            writer.writeheader()
            writer.writerows(rows)
        _write_summary_md(rows, out_dir, table_stem, title, date)

    decomp = _decomposition_rows(registry, metrics_by_key)
    if decomp:
        with (out_dir / f"分解表_{table_stem}.csv").open(
            "w", encoding="utf-8-sig", newline=""
        ) as fh:
            writer = csv.DictWriter(fh, fieldnames=DECOMP_COLUMNS)
            writer.writeheader()
            writer.writerows(decomp)
        _write_decomp_md(decomp, out_dir, table_stem, title, date)
    print(f"重建汇总：{len(rows)} 行、分解：{len(decomp)} 行 → {out_dir}", flush=True)


def _fmt(v, nd=4) -> str:
    if v is None or v == "":
        return "—"
    if isinstance(v, str):
        return v
    return f"{float(v):.{nd}f}"


def _write_summary_md(rows: list[dict], out_dir: Path = OUT_DIR,
                      table_stem: str = "温压T2", title: str = "T2 温压敏感性",
                      date: str = "2026-10-01") -> None:
    md = [
        f"# {title}汇总（{date}）",
        "",
        "口径：09-16 脚本链（CORRECTED_KW + CFL、tt=min(泵总+1200, stop_t)）；"
        "**基线 = 本批 `Toff_zero`**；判别量只在本表（solver summary 零新增键）。",
        "静温统一式 T(z)=16.006+1.7598e-2·z（呼101/呼103 标「无瞬态表」）；"
        "Q10乙=enable_yield_gate=False 通道代理；绝对替换不动、HB/水泥τy/屈服门进方程三开关未开。",
        "",
        "| 井名 | 变体 | 温度档 | T | gate | η_E | η_N | Δη_N/pp | 饥饿份额 | "
        "wall_b加权 | 备注 | 耗时/s |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        dn = f"{r['Δη_N_pp']:+.2f}" if r["Δη_N_pp"] != "" else "—"
        md.append(
            f"| {r['井名']} | {r['变体']} | {r['温度档']} | {r['T开关']} | "
            f"{r['yield_gate']} | {r['η_E']:.4f} | {r['η_N']:.4f} | {dn} | "
            f"{_fmt(r['饥饿份额'])} | {_fmt(r['屈服门活化率_b加权'])} | "
            f"{r['备注'] or '—'} | {r['耗时_s']} |"
        )
    (out_dir / f"汇总表_{table_stem}.md").write_text("\n".join(md), encoding="utf-8")


# ---------------------------------------------------------------- 分解表
def _decomposition_rows(registry, metrics_by_key) -> list[dict]:
    """同排量两分量：口径差 = Ton_const60 − Toff；温度场效应 = Ton_mode − Ton_const60。

    排量 1.0 的 Toff 对照取 ``Toff_zero``（rate_x1.0 与恒等变换逐位同值）。
    差值列口径：η_E / η_N 用 pp（×100），判别量用原单位。
    const60 分母锚（``a``）缺档时不追加空 Δ 行（``a``/``c`` 任一为 None 即跳过）。
    """
    out: list[dict] = []

    def _vals(well: str, name: str) -> dict | None:
        return metrics_by_key.get((well, name))

    def _diff(a: dict | None, b: dict | None, metric: str):
        if a is None or b is None:
            return ""
        va, vb = a.get(metric), b.get(metric)
        if va is None or vb is None or va == "" or vb == "":
            return ""
        d = float(va) - float(vb)
        return d * 100.0 if metric in ("η_E", "η_N") else d

    for well, variants in registry.items():
        names = {v[0] for v in variants}
        # 该井存在的纯温度×排量档（const60 档为分母锚）
        rates = sorted(
            {m.group(2) for n in names if (m := _RE_TON_RATE.match(n))
             and m.group(1) == "const60"},
            key=float,
        )
        for rate in rates:
            const60 = f"Ton_const60_rate_x{rate}"
            toff = BASELINE_VARIANT if float(rate) == 1.0 else f"Toff_rate_x{rate}"
            a, b = _vals(well, const60), _vals(well, toff)
            if a is not None and b is not None:
                row = {"井名": well, "排量档": rate, "分量": "口径差",
                       "对照": f"{const60} − {toff}"}
                for m in DECOMP_METRICS:
                    row[f"Δ_{m}"] = _diff(a, b, m)
                out.append(row)
            for mode in ("static", "table"):
                name = f"Ton_{mode}_rate_x{rate}"
                c = _vals(well, name)
                if c is None or a is None:   # const60 缺档 ⇒ 无分母锚，不追加空 Δ 行
                    continue
                row = {"井名": well, "排量档": rate, "分量": "温度场效应",
                       "对照": f"{name} − {const60}"}
                for m in DECOMP_METRICS:
                    row[f"Δ_{m}"] = _diff(c, a, m)
                out.append(row)
    return out


def _write_decomp_md(rows: list[dict], out_dir: Path = OUT_DIR,
                     table_stem: str = "温压T2", title: str = "T2 温压敏感性",
                     date: str = "2026-10-01") -> None:
    metric_cols = [f"Δ_{m}" for m in DECOMP_METRICS]
    md = [
        f"# {title}分解表（{date}）",
        "",
        "同排量两分量（脚本层算术，从行数据聚合）：",
        "- **口径差** = `Ton_const60 − Toff`（温度开关打开但取常数 60 °C 的净效应）",
        "- **温度场效应** = `Ton_{static,table} − Ton_const60`（真正来自温度场形状的部分）",
        "",
        "Δ_η_E、Δ_η_N 单位 pp；其余 Δ 为原单位。",
        "",
        "| " + " | ".join(["井名", "排量档", "分量", "对照"] + metric_cols) + " |",
        "|" + "---|" * (4 + len(metric_cols)),
    ]
    for r in rows:
        cells = [str(r["井名"]), str(r["排量档"]), r["分量"], r["对照"]]
        for m in metric_cols:
            v = r.get(m, "")
            cells.append("—" if v == "" or v is None else f"{float(v):+.4f}")
        md.append("| " + " | ".join(cells) + " |")
    (out_dir / f"分解表_{table_stem}.md").write_text("\n".join(md), encoding="utf-8")


# ---------------------------------------------------------------- 主流程
def _assert_out_dir_allowed(out_dir: Path) -> None:
    """冻结目录防护：输出目录不得落在三个历史敏感性目录内（红线：零写入）。"""
    for _d in FROZEN_DIRS:
        if out_dir == _d or _d in out_dir.parents:
            raise RuntimeError(f"输出目录防护触发：{out_dir} 落在冻结目录 {_d} 内")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="T2 温压敏感性批（双批次，断点续跑）"
    )
    ap.add_argument("--batch", choices=sorted(BATCHES), default="t2",
                    help="批次：t2=2026-10-01 原 45 变体（默认，逐位不变）；"
                         "phase0=2026-10-06 三重点井补跑矩阵（52 变体）")
    ap.add_argument("--out-dir", default=None,
                    help="覆盖输出目录（默认取批次定义）。phase0 应指向"
                         "新日期目录（红线：补跑新目录带日期后缀）；"
                         "不带本参数时 t2 批行为与 2026-10-01 原版完全相同")
    ap.add_argument("--well", choices=[*WELLS.keys(), "all"], default="all")
    ap.add_argument("--variant", default=None,
                    help="只跑/只列该变体名（冒烟与补跑用）")
    ap.add_argument("--force", action="store_true",
                    help="忽略已存在的 *_结果摘要.json/*_判别量.json，全部重算并"
                         "覆盖。改码后须 --force 或删除 *_判别量.json 才会重算"
                         "（判别量按代码公式算出，改码不改档名不会自动触发）")
    args = ap.parse_args()

    batch = BATCHES[args.batch]
    out_dir = Path(args.out_dir) if args.out_dir else batch["out_dir"]
    _assert_out_dir_allowed(out_dir)
    registry = build_temperature_variants(args.batch)
    well_keys = list(WELLS) if args.well == "all" else [args.well]

    # 注册表自检（响亮失败，防手滑改坏矩阵）
    actual = {w: len(v) for w, v in registry.items()}
    if actual != EXPECTED_VARIANTS[args.batch]:
        raise RuntimeError(
            f"{args.batch} 变体矩阵应为 {EXPECTED_VARIANTS[args.batch]}，"
            f"实际 {actual}"
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    for well in well_keys:
        if well not in registry:
            raise SystemExit(
                f"{well} 不在批次 {args.batch} 的注册表内"
                f"（该批井：{sorted(registry)}）"
            )
        mod, fn = WELLS[well]
        loader = getattr(importlib.import_module(mod), fn)
        variants = registry[well]
        if args.variant is not None:
            variants = [v for v in variants if v[0] == args.variant]
            if not variants:
                raise SystemExit(f"{well} 无变体 {args.variant!r}（见注册表）")
        print(f"\n=== [{args.batch}] {well}  {len(variants)} 变体 ===", flush=True)
        for name, well_fn, fluid_fn, sched_fn, run_opts in variants:
            _run_or_reuse(
                well, loader, name, (well_fn, fluid_fn, sched_fn), run_opts,
                force=args.force, out_dir=out_dir,
            )

    # 幂等重建（从本批全量落盘 JSON 收；子集/断点跑只出已有行）
    rebuild_tables(
        registry, out_dir, batch["table_stem"], batch["title"], batch["date"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
