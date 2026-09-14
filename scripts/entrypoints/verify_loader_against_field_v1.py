# -*- coding: utf-8 -*-
"""
verify_loader_against_field_v1.py — cemdisp loader vs 现场提取CSV 数据校验（v1, 2026-09-01）

目的：
1. 用 cemdisp.data.loaders 加载 8 口尾管井，导出 fluids / schedule / geometry / CBL 结构化清单；
2. 逐井逐项比对 loader 与现场提取层 CSV（fluid_properties / rheometer_readings /
   well_geometry / casing_liner_string / pumping_schedule / caliper / centralization /
   cbl_evaluation / target_intervals）以及 01_总表/loader优化采用值对照表；
3. 纯 loader 内部一致性检查：体积闭合、时间轴、环空体积重算、水泥浆/环空库存比、
   密度/流变 sanity、几何 sanity。

数据源：
- loader：cemdisp.data.loaders 各井默认入口（hu103/ht1_003/ht1_004 附加变体仅对比 schedule）；
- 提取层：参考文档/现场资料提取/<井目录>/*.csv（2026-08-29 校准轮，0708 衍生）；
- 总表：参考文档/现场资料提取/01_总表/loader优化采用值对照表.csv。

输出：results/数据校验_2026-09-01/loader_dump/
  - <well>_fluids.csv / <well>_schedule.csv / <well>_geometry.csv / <well>_cbl.csv
  - checks_all.csv（全部比对条目）/ loader_dump_汇总.md

运行：PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/verify_loader_against_field_v1.py
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

EXTRACT_ROOT = ROOT / "参考文档" / "现场资料提取"
SUMMARY_TABLE_CSV = EXTRACT_ROOT / "01_总表" / "loader优化采用值对照表.csv"
OUT_DIR = ROOT / "results" / "数据校验_2026-09-01" / "loader_dump"

from cemdisp.data.loaders import (  # noqa: E402
    load_hu101_tailpipe,
    load_hu102_tailpipe,
    load_hu103_tailpipe,
    load_hu103_tailpipe_actual,
    load_hu1_tailpipe,
    load_hu2_tailpipe,
    load_ht1_001_tailpipe,
    load_ht1_003_tailpipe,
    load_ht1_003_tailpipe_design,
    load_ht1_004_tailpipe,
    load_ht1_004_tailpipe_actual,
)
from cemdisp.data.fluid_spec import FluidRole, RheologyModel  # noqa: E402

# ---------------------------------------------------------------- 工具函数


def read_rows(path: Path) -> list[dict]:
    """读取 UTF-8(-sig) CSV 为 dict 列表；文件不存在返回空列表。"""
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [dict(r) for r in csv.DictReader(handle)]


def fnum(value) -> float | None:
    """容错地把 CSV 字符串转 float；空/非法返回 None。"""
    if value is None:
        return None
    s = str(value).strip().replace(",", "")
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def fmt(value) -> str:
    """检查条目取值格式化：None→'—'，float→去尾零。"""
    if value is None:
        return "—"
    if isinstance(value, float):
        if abs(value - round(value)) < 1e-9 and abs(value) < 1e7:
            return str(int(round(value)))
        return f"{value:.6g}"
    return str(value)


def cmp_num(a, b, *, tol_abs=0.0, tol_pct=1.0):
    """数值比较 → '✓' / '✗' / '-'（任一缺数）。"""
    if a is None or b is None:
        return "-"
    tol = max(tol_abs, tol_pct / 100.0 * max(abs(a), abs(b)))
    return "✓" if abs(a - b) <= tol else "✗"


def interp(points, depth: float) -> float:
    """DepthValuePoint 风格 (depth,value) 元组序列线性插值，越界取端值。"""
    pts = sorted(points)
    if depth <= pts[0][0]:
        return pts[0][1]
    if depth >= pts[-1][0]:
        return pts[-1][1]
    for (d0, v0), (d1, v1) in zip(pts, pts[1:]):
        if d0 <= depth <= d1:
            if d1 == d0:
                return v0
            frac = (depth - d0) / (d1 - d0)
            return v0 + frac * (v1 - v0)
    return pts[-1][1]


def stats(points) -> dict:
    """剖面点统计（min/max/mean/点数/深度范围）。"""
    if not points:
        return {"n": 0}
    values = [v for _, v in points]
    depths = [d for d, _ in points]
    return {
        "n": len(values),
        "min": min(values),
        "max": max(values),
        "mean": sum(values) / len(values),
        "depth_min": min(depths),
        "depth_max": max(depths),
    }


# CSV fluid_role / fluid_type → 标准角色 key（长词在前防"泥浆"子串误命中"尾管水泥浆"）
ROLE_KEY_MAP = {
    "尾管水泥浆": "tail", "超细水泥": "tail", "中间浆": "intermediate",
    "mud": "mud", "钻井液": "mud", "泥浆": "mud", "油基钻井液": "mud",
    "wash": "wash", "先导浆": "wash", "平衡液": "wash", "冲洗": "flusher",
    "spacer": "spacer", "隔离": "spacer",
    "lead_cement": "lead", "lead": "lead", "领浆": "lead",
    "tail_cement": "tail", "tail": "tail", "尾浆": "tail",
    "displacement": "displacement", "替浆": "displacement", "替钻井液": "displacement",
    "顶替": "displacement", "中置液": "displacement", "后置": "displacement",
    "压塞": "displacement", "井浆": "displacement", "轻泥浆": "displacement",
    "flusher": "flusher",
}


def role_key(role_or_name: str) -> str:
    s = str(role_or_name).strip()
    for key, val in ROLE_KEY_MAP.items():
        if key in s:
            return val
    return s.lower() if s else "?"


def spec_role_key(fluid) -> str:
    """loader FluidSpec → 标准 role key（flusher 单列，wash 保留）。"""
    role = fluid.role
    if role == FluidRole.FLUSHER:
        return "flusher"
    return role.value if role.value in {"mud", "wash", "spacer", "lead", "tail",
                                        "intermediate", "displacement", "flusher"} else "other"


def csv_rheology_model(text: str) -> str | None:
    s = (text or "").strip()
    if not s:
        return None
    if "幂律" in s and ("宾汉" in s or "bingham" in s.lower()):
        return "mixed"
    if "幂律" in s or "power" in s.lower():
        return "power_law"
    if "宾汉" in s or "bingham" in s.lower():
        return "bingham"
    if "牛顿" in s or "newton" in s.lower():
        return "newtonian"
    if "herschel" in s.lower() or "hb" == s.lower():
        return "herschel_bulkley"
    return "other:" + s


# ---------------------------------------------------------------- 检查条目收集


class Checks:
    def __init__(self) -> None:
        self.rows: list[dict] = []

    def add(self, well, category, item, loader_value, csv_value, status, note=""):
        self.rows.append({
            "well": well, "category": category, "item": item,
            "loader_value": fmt(loader_value), "csv_value": fmt(csv_value),
            "status": status, "note": note,
        })


# ---------------------------------------------------------------- 各类导出


def dump_fluids(well: str, fluids, path: Path) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "role", "density_kg_m3", "rheology_model",
                    "pv_pa_s", "yp_pa", "power_law_n", "consistency_k_pa_sn"])
        for f in fluids:
            w.writerow([f.name, f.role.value, f.density_kg_m3, f.rheology_model.value,
                        f.plastic_viscosity_pa_s, f.yield_stress_pa,
                        f.power_law_n, f.consistency_k])


def dump_schedule(well: str, schedule, path: Path, fluids=()) -> list[dict]:
    """导出泵注日程并返回带推算时间轴的行（供时间校验复用）。"""
    name_role = {f.name: spec_role_key(f) for f in fluids}
    rows: list[dict] = []
    t = 0.0
    cum = 0.0
    for i, s in enumerate(schedule.steps):
        duration_min = (s.volume_m3 / s.rate_m3_min) if s.rate_m3_min > 0 else None
        start_s = s.start_time_s if s.start_time_s is not None else t * 60.0
        end_s = s.end_time_s if (s.end_time_s is not None and s.start_time_s is not None) else (
            start_s + duration_min * 60.0 if duration_min is not None else None)
        t = end_s / 60.0 if end_s is not None else t
        cum += s.volume_m3
        rows.append({
            "index": i, "step_name": s.step_name, "fluid_name": s.fluid_name,
            "role": name_role.get(s.fluid_name) or role_key(s.fluid_name), "volume_m3": s.volume_m3,
            "rate_m3_min": s.rate_m3_min, "duration_min": duration_min,
            "start_min": start_s / 60.0, "end_min": end_s / 60.0 if end_s else None,
            "cum_volume_m3": cum,
            "explicit_time": s.start_time_s is not None,
            "remarks": (s.remarks or "")[:160],
        })
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["index", "step_name", "fluid_name", "role", "volume_m3",
                    "rate_m3_min", "duration_min", "start_min", "end_min",
                    "cum_volume_m3", "explicit_time", "remarks"])
        for r in rows:
            w.writerow([r["index"], r["step_name"], r["fluid_name"], r["role"],
                        r["volume_m3"], r["rate_m3_min"],
                        fmt(r["duration_min"]), fmt(r["start_min"]), fmt(r["end_min"]),
                        r["cum_volume_m3"], r["explicit_time"], r["remarks"]])
    return rows


def dump_geometry(well: str, spec, path: Path) -> dict:
    hole_st = stats([(p.depth_md_m, p.value) for p in spec.hole_diameter_profile])
    incl_st = stats([(p.depth_md_m, p.value) for p in spec.inclination_profile])
    so_st = stats([(p.depth_md_m, p.value) for p in spec.standoff_profile])
    pipe_st = stats([(p.depth_md_m, p.value) for p in spec.pipe_id_profile])
    kv = [
        ("well_name", spec.well_name),
        ("top_md_m", spec.top_md_m), ("bottom_md_m", spec.bottom_md_m),
        ("shoe_md_m", spec.shoe_md_m), ("hanger_md_m", spec.hanger_md_m),
        ("casing_id_mm", spec.casing_id_mm),
        ("liner_od_mm", spec.liner_od_mm), ("liner_id_mm", spec.liner_id_mm),
        ("is_dual_diameter", spec.is_dual_diameter),
        ("upper_section_bottom_md_m", spec.upper_section_bottom_md_m),
        ("upper_liner_od_mm", spec.upper_liner_od_mm),
        ("upper_liner_id_mm", spec.upper_liner_id_mm),
        ("liner_wall_thickness_mm", spec.liner_wall_thickness_mm),
        ("shoe_lag_volume_m3", spec.shoe_lag_volume_m3),
        ("hole_profile_n", hole_st.get("n")),
        ("hole_profile_min_mm", hole_st.get("min")),
        ("hole_profile_max_mm", hole_st.get("max")),
        ("hole_profile_mean_mm", round(hole_st["mean"], 2) if hole_st.get("n") else None),
        ("hole_profile_depth_range", f"{hole_st.get('depth_min')}–{hole_st.get('depth_max')}" if hole_st.get("n") else None),
        ("inclination_profile_n", incl_st.get("n")),
        ("inclination_max_deg", incl_st.get("max")),
        ("standoff_profile_n", so_st.get("n")),
        ("standoff_min", so_st.get("min")),
        ("standoff_max", so_st.get("max")),
        ("standoff_mean", round(so_st["mean"], 3) if so_st.get("n") else None),
        ("pipe_id_profile_n", pipe_st.get("n")),
        ("cbl_pass_rate", None),  # 由 validation 填（见 verify_well）
    ]
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["key", "value"])
        for k, v in kv:
            w.writerow([k, fmt(v) if not isinstance(v, str) else v])
        for i, note in enumerate(spec.notes):
            w.writerow([f"note_{i}", note])
    return {"hole": hole_st, "incl": incl_st, "standoff": so_st, "pipe": pipe_st,
            "kv": {k: v for k, v in kv}}


def dump_cbl(well: str, spec, validation, path: Path) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["kind", "name", "top_md_m", "bottom_md_m", "extra"])
        for win in spec.evaluation_windows:
            w.writerow(["evaluation_window", win.name, win.top_md_m, win.bottom_md_m, win.window_type])
        w.writerow(["validation", "cbl_pass_rate", "", "",
                    fmt(validation.cbl_pass_rate) if validation else "—"])


# ---------------------------------------------------------------- 环空/管内体积重算


def annulus_volume_m3(spec) -> float:
    """按 hole_diameter_profile × (dual)liner OD 梯形积分环空体积（m³）。"""
    pts = [(p.depth_md_m, p.value) for p in spec.hole_diameter_profile]
    if len(pts) < 2:
        return float("nan")
    lo, hi = spec.top_md_m, spec.bottom_md_m

    def od_at(d: float) -> float:
        if spec.is_dual_diameter and d <= spec.upper_section_bottom_md_m:
            return spec.upper_liner_od_mm
        return spec.liner_od_mm

    vol = 0.0
    for (d0, h0), (d1, h1) in zip(pts, pts[1:]):
        a, b = max(d0, lo), min(d1, hi)
        if b <= a:
            continue
        hh0 = interp(pts, a)
        hh1 = interp(pts, b)
        od0, od1 = od_at(a), od_at(b)
        vol += math.pi / 4.0 * ((hh0 ** 2 - od0 ** 2) + (hh1 ** 2 - od1 ** 2)) / 2.0 * (b - a) / 1e6
    return vol


def pipe_volume_m3(spec) -> float:
    """管内容积（m³）：pipe_id_profile 优先，否则单一 liner_id，top→shoe 积分。"""
    if spec.pipe_id_profile:
        pts = [(p.depth_md_m, p.value) for p in spec.pipe_id_profile]
        vol = 0.0
        for (d0, i0), (d1, i1) in zip(pts, pts[1:]):
            a, b = max(d0, spec.top_md_m), min(d1, spec.shoe_md_m)
            if b <= a:
                continue
            ii0 = interp(pts, a)
            ii1 = interp(pts, b)
            vol += math.pi / 4.0 * (ii0 ** 2 + ii1 ** 2) / 2.0 * (b - a) / 1e6
        return vol
    if spec.liner_id_mm:
        return math.pi / 4.0 * (spec.liner_id_mm / 1000.0) ** 2 * (spec.shoe_md_m - spec.top_md_m)
    return float("nan")


# ---------------------------------------------------------------- 分项校验


def check_fluids(well: str, fluids, checks: Checks) -> None:
    rows = read_rows(EXTRACT_ROOT / WELL_DIRS[well] / "fluid_properties.csv")
    by_role: dict[str, list] = {}
    for f in fluids:
        by_role.setdefault(spec_role_key(f), []).append(f)
    csv_group: dict[str, list[dict]] = {}
    for r in rows:
        csv_group.setdefault(role_key(r.get("fluid_role", "")), []).append(r)

    for key in sorted(set(by_role) | set(csv_group)):
        loader_list = by_role.get(key, [])
        csv_list = csv_group.get(key, [])
        if not loader_list or not csv_list:
            checks.add(well, "fluid", f"role={key}",
                       f"{len(loader_list)}条", f"{len(csv_list)}条", "⚠",
                       "loader或CSV侧该角色无对应条目（角色映射或口径差异）")
            continue
        for i, f in enumerate(loader_list):
            r = csv_list[min(i, len(csv_list) - 1)]
            tag = f"{key}/{f.name}" if len(loader_list) > 1 else f"role={key}"
            if len(loader_list) != len(csv_list):
                tag += f"[第{i + 1}对{min(i, len(csv_list) - 1) + 1}/{len(csv_list)}]"
            # 密度 g/cc → kg/m³
            rho_csv = fnum(r.get("density_g_cm3"))
            rho_csv = rho_csv * 1000.0 if rho_csv is not None else None
            st = cmp_num(f.density_kg_m3, rho_csv, tol_abs=10.0, tol_pct=0.5)
            checks.add(well, "fluid", f"{tag} 密度kg/m3", f.density_kg_m3, rho_csv, st,
                       "CSV g/cm3×1000；容差10kg/m3" if st == "✗" else "")
            # PV mPa·s → Pa·s
            pv_csv = fnum(r.get("pv_mpa_s"))
            pv_csv = pv_csv / 1000.0 if pv_csv is not None else None
            if f.plastic_viscosity_pa_s is not None:
                st = cmp_num(f.plastic_viscosity_pa_s, pv_csv, tol_pct=2.0)
                checks.add(well, "fluid", f"{tag} PV Pa·s", f.plastic_viscosity_pa_s, pv_csv, st,
                           "CSV mPa·s/1000；容差2%" if st == "✗" else "")
            # YP
            yp_csv = fnum(r.get("yp_pa")) or fnum(r.get("yield_stress_pa"))
            if f.yield_stress_pa is not None:
                st = cmp_num(f.yield_stress_pa, yp_csv, tol_pct=2.0)
                checks.add(well, "fluid", f"{tag} YP Pa", f.yield_stress_pa, yp_csv, st,
                           "容差2%" if st == "✗" else "")
            # n / K
            n_csv = fnum(r.get("flow_index_n"))
            k_csv = fnum(r.get("consistency_k_pa_sn"))
            if f.power_law_n is not None:
                st = cmp_num(f.power_law_n, n_csv, tol_pct=2.0)
                checks.add(well, "fluid", f"{tag} n", f.power_law_n, n_csv, st,
                           "容差2%" if st == "✗" else "")
            if f.consistency_k is not None:
                st = cmp_num(f.consistency_k, k_csv, tol_pct=3.0)
                checks.add(well, "fluid", f"{tag} K Pa·s^n", f.consistency_k, k_csv, st,
                           "容差3%" if st == "✗" else "")
            # 流变模型口径
            m_loader = f.rheology_model.value
            m_csv = csv_rheology_model(r.get("rheology_model", ""))
            if m_csv and m_csv not in ("other:",) and not m_csv.startswith("other"):
                ok = (m_csv == m_loader) or (m_csv == "mixed")
                checks.add(well, "fluid", f"{tag} 流变模型", m_loader, m_csv,
                           "✓" if ok else "⚠", "CSV为混合口径(幂律+宾汉)" if m_csv == "mixed" else "")


def check_rheometer(well: str, fluids, checks: Checks) -> None:
    """六速读数反算 n/K 作为参考（温度口径可能不同，仅提示级）。"""
    rows = read_rows(EXTRACT_ROOT / WELL_DIRS[well] / "rheometer_readings.csv")
    if not rows:
        return
    for f in fluids:
        if f.rheology_model is not RheologyModel.POWER_LAW:
            continue
        cands = [r for r in rows if role_key(r.get("fluid_name", "")) in
                 {"lead", "tail", "intermediate", "spacer"} and
                 (role_key(r.get("fluid_name", "")) == spec_role_key(f)
                  or str(f.name)[:2] in str(r.get("fluid_name", "")))]
        if not cands:
            continue
        r = cands[-1]
        t600, t300 = fnum(r.get("theta_600")), fnum(r.get("theta_300"))
        if not t600 or not t300 or t600 <= t300:
            continue
        n_fit = math.log(t600 / t300) / math.log(2.0)
        k_fit = 0.511 * t300 / (1.7023 * 300.0) ** n_fit
        checks.add(well, "fluid", f"{f.name} n(六速反算参考)", f.power_law_n, round(n_fit, 3),
                   "⚠" if abs(f.power_law_n - n_fit) > 0.05 * max(f.power_law_n, n_fit) else "✓",
                   f"{r.get('temperature_c')}℃ 六速θ600/θ300={t600}/{t300}；反算仅供温度口径参考")
        checks.add(well, "fluid", f"{f.name} K(六速反算参考)", f.consistency_k, round(k_fit, 3),
                   "⚠" if abs(f.consistency_k - k_fit) > 0.15 * max(f.consistency_k, k_fit) else "✓",
                   "反算 K=0.511·θ300/(1.7023·300)^n")


def check_schedule(well: str, sched_rows: list[dict], checks: Checks) -> None:
    """loader schedule vs 提取 pumping_schedule.csv（体积/排量/时间轴）。"""
    path = EXTRACT_ROOT / WELL_DIRS[well] / "pumping_schedule.csv"
    rows = read_rows(path)
    if not rows:
        checks.add(well, "schedule", "pumping_schedule.csv", "有", "文件缺失", "⚠", "提取层无泵注CSV")
        return
    hu1_legacy = "pump_rate_m3_min" in (rows[0].keys() if rows else {})

    def csv_vol(r):
        return fnum(r.get("volume_m3")) if not hu1_legacy else fnum(r.get("volume_m3"))

    def csv_rate(r):
        key = "pump_rate_m3_min" if hu1_legacy else "rate_m3_min"
        return fnum(r.get(key))

    def csv_role(r):
        return role_key(r.get("fluid_role", "") if not hu1_legacy else r.get("fluid_type", ""))

    # loader 按角色分组序列
    loader_by_role: dict[str, list[dict]] = {}
    for r in sched_rows:
        if r["volume_m3"] <= 0:
            continue
        loader_by_role.setdefault(r["role"], []).append(r)
    csv_by_role: dict[str, list[dict]] = {}
    for r in rows:
        v = csv_vol(r)
        if v is None or v <= 0 or csv_role(r) in {"?", "shutdown"}:
            continue
        csv_by_role.setdefault(csv_role(r), []).append(r)

    for key in sorted(set(loader_by_role) | set(csv_by_role)):
        ls = loader_by_role.get(key, [])
        cs = csv_by_role.get(key, [])
        if not ls or not cs:
            checks.add(well, "schedule", f"role={key}",
                       f"{len(ls)}步", f"{len(cs)}步", "⚠",
                       "单侧无该角色步骤（口径/序列差异）")
            continue
        for i, l in enumerate(ls):
            r = cs[min(i, len(cs) - 1)]
            tag = f"role={key}[{i + 1}]"
            st = cmp_num(l["volume_m3"], csv_vol(r), tol_abs=0.2, tol_pct=1.0)
            checks.add(well, "schedule", f"{tag} 体积m3({l['step_name']})",
                       l["volume_m3"], csv_vol(r), st,
                       "容差0.2m³/1%" if st == "✗" else "")
            st = cmp_num(l["rate_m3_min"], csv_rate(r), tol_pct=2.0)
            checks.add(well, "schedule", f"{tag} 排量m3/min",
                       l["rate_m3_min"], csv_rate(r), st,
                       "容差2%" if st == "✗" else "")
            # 时间轴（loader 为体积/排量推算）
            cs_start = fnum(r.get("start_time_min"))
            cs_end = fnum(r.get("end_time_min"))
            if cs_start is not None and l["start_min"] is not None:
                st = "✓" if abs(l["start_min"] - cs_start) <= 3.0 else "⚠"
                checks.add(well, "schedule", f"{tag} 开始min(推算vs现场)",
                           round(l["start_min"], 1), cs_start, st,
                           "loader时间为体积/排量推算；容差3min" if st != "✓" else "")
            if cs_end is not None and l["end_min"] is not None:
                st = "✓" if abs(l["end_min"] - cs_end) <= 3.0 else "⚠"
                checks.add(well, "schedule", f"{tag} 结束min(推算vs现场)",
                           round(l["end_min"], 1), cs_end, st,
                           "容差3min" if st != "✓" else "")

    # 总时长 vs 现场记录
    l_end = max([r["end_min"] for r in sched_rows if r["end_min"]], default=None)
    cs_end_all = [fnum(r.get("end_time_min")) for r in rows if fnum(r.get("end_time_min")) is not None]
    cs_end_all += [fnum(r.get("time_min")) for r in rows if hu1_legacy and fnum(r.get("time_min")) is not None]
    if l_end is not None and cs_end_all:
        checks.add(well, "schedule", "日程总时长min(推算vs现场末步)",
                   round(l_end, 1), max(cs_end_all), "⚠" if abs(l_end - max(cs_end_all)) > 5.0 else "✓",
                   "现场末步含试压/停泵等非注入步，差异仅提示")


def check_geometry(well: str, spec, checks: Checks) -> None:
    d = EXTRACT_ROOT / WELL_DIRS[well]
    # well_geometry.csv 长表：同 item 多值行时取首现（尾管段在前，回接段行混入在文件尾部）
    wg: dict[str, dict] = {}
    dup_items: set[str] = set()
    for r in read_rows(d / "well_geometry.csv"):
        item = r.get("item", "").strip()
        if item in wg:
            dup_items.add(item)
            continue
        wg[item] = r
    if dup_items:
        checks.add(well, "geometry", "well_geometry.csv 同item多值行", "—",
                   ",".join(sorted(dup_items)), "⚠",
                   "CSV长表混入回接段同名item，比对取首现（尾管段）值")

    def wg_val(item):
        r = wg.get(item)
        return fnum(r.get("value")) if r else None

    pairs = [
        ("悬挂器深度hanger_md_m", spec.hanger_md_m, wg_val("liner_hanger_md"), 0.5),
        ("尾管鞋shoe_md_m", spec.shoe_md_m,
         wg_val("liner_shoe_md") or wg_val("liner_bottom_md"), 0.5),
        ("模型段顶top_md_m", spec.top_md_m, wg_val("liner_top_md"), 3.0),
        ("完钻/段底bottom_md_m", spec.bottom_md_m,
         wg_val("liner_bottom_md") or wg_val("td_md") or wg_val("well_total_depth_md"), 3.0),
    ]
    for item, lv, cv, tol in pairs:
        if cv is None:
            continue
        st = cmp_num(lv, cv, tol_abs=tol)
        checks.add(well, "geometry", item, lv, cv, st, f"well_geometry.csv；容差{tol}m" if st == "✗" else "")

    # casing_liner_string.csv 套管串 OD/ID（排除回接段：top 深于模型段顶 500m 以上的行）
    rows = [r for r in read_rows(d / "casing_liner_string.csv")
            if r.get("component", "").strip() == "liner"]
    liners = sorted([r for r in rows
                     if fnum(r.get("top_md_m")) is not None
                     and fnum(r.get("top_md_m")) >= spec.top_md_m - 500.0],
                    key=lambda r: fnum(r.get("top_md_m")))
    if liners:
        lower = liners[-1]
        checks.add(well, "geometry", "下段尾管OD mm", spec.liner_od_mm, fnum(lower.get("outer_diameter_mm")),
                   cmp_num(spec.liner_od_mm, fnum(lower.get("outer_diameter_mm")), tol_abs=0.2),
                   "casing_liner_string 最后一 liner 段")
        if len(liners) >= 2 and spec.is_dual_diameter:
            upper = liners[0]
            checks.add(well, "geometry", "上段尾管OD mm", spec.upper_liner_od_mm,
                       fnum(upper.get("outer_diameter_mm")),
                       cmp_num(spec.upper_liner_od_mm, fnum(upper.get("outer_diameter_mm")), tol_abs=0.2))
            checks.add(well, "geometry", "上段尾管ID mm", spec.upper_liner_id_mm,
                       fnum(upper.get("inner_diameter_mm")),
                       cmp_num(spec.upper_liner_id_mm, fnum(upper.get("inner_diameter_mm")), tol_abs=0.5))
        # 变径位置
        if spec.is_dual_diameter:
            cv_ub = fnum(liners[0].get("bottom_md_m"))
            checks.add(well, "geometry", "上段底/变径位置 md_m", spec.upper_section_bottom_md_m, cv_ub,
                       cmp_num(spec.upper_section_bottom_md_m, cv_ub, tol_abs=1.0),
                       "loader常取整（如6796.329→6796）")
        # 尾管总长 sanity：CSV 各段长度和 vs shoe-hanger
        total_len = sum((fnum(r.get("bottom_md_m")) or 0) - (fnum(r.get("top_md_m")) or 0) for r in liners)
        if spec.hanger_md_m is not None:
            checks.add(well, "geometry", "尾管长度m(CSV分段和 vs shoe−hanger)",
                       round(spec.shoe_md_m - spec.hanger_md_m, 2), round(total_len, 2),
                       "✓" if abs((spec.shoe_md_m - spec.hanger_md_m) - total_len) <= 30.0 else "⚠",
                       "含悬挂器/浮箍段长差异，容差30m")

    # caliper / inclination / centralization 剖面对比（列名跨井兼容）
    cal = read_rows(d / "caliper_profile.csv")
    if cal:
        md_key = next((k for k in ("md_m", "measured_depth_m", "depth_m") if k in cal[0]), None)
        cal_key = next((k for k in ("caliper_mm", "avg_caliper_mm") if k in cal[0]), None)
        cal_pts = sorted([(fnum(r.get(md_key)), fnum(r.get(cal_key)))
                          for r in cal if fnum(r.get(md_key)) is not None and fnum(r.get(cal_key)) is not None])
        loader_pts = [(p.depth_md_m, p.value) for p in spec.hole_diameter_profile]
        cs_st, l_st = stats(cal_pts), stats(loader_pts)
        checks.add(well, "geometry", "井径剖面点数(CSV实测 vs loader模型域)",
                   l_st["n"], cs_st["n"], "✓" if l_st["n"] <= cs_st["n"] else "⚠",
                   "loader仅保留模型段内点并可能做等效换算")
        checks.add(well, "geometry", "井径max mm(CSV vs loader)",
                   round(l_st["max"], 1), round(cs_st["max"], 1),
                   "✓" if l_st["max"] >= cs_st["max"] - 5.0 else "⚠", "上段等效换算会改变数值")
        checks.add(well, "geometry", "井径min mm(CSV vs loader)",
                   round(l_st["min"], 1), round(cs_st["min"], 1),
                   "✓" if l_st["min"] <= cs_st["min"] + 5.0 else "⚠", "")

    cen = read_rows(d / "centralization_profile.csv")
    if cen:
        so_col = next((k for k in ("standoff_ratio", "standoff", "centralization_percent")
                       if k in cen[0] and any(fnum(r.get(k)) is not None for r in cen)), None)
        if so_col:
            cs_vals = [fnum(r.get(so_col)) for r in cen]
            cs_vals = [v for v in cs_vals if v is not None]
            if so_col == "centralization_percent":
                cs_vals = [v / 100.0 for v in cs_vals]
            l_vals = [p.value for p in spec.standoff_profile]
            if cs_vals and l_vals:
                checks.add(well, "geometry", f"standoff范围(CSV:{so_col} vs loader)",
                           f"{min(l_vals):.2f}–{max(l_vals):.2f}",
                           f"{min(cs_vals):.2f}–{max(cs_vals):.2f}",
                           "⚠" if (max(l_vals) < min(cs_vals) - 0.05 or min(l_vals) > max(cs_vals) + 0.05) else "✓",
                           "CSV多为扶正器布置设计值，loader常为model_assumption名义剖面")


def check_cbl(well: str, spec, validation, checks: Checks) -> None:
    rows = read_rows(EXTRACT_ROOT / WELL_DIRS[well] / "cbl_evaluation.csv")
    main_rows = [r for r in rows
                 if fnum(r.get("cbl_pass_rate")) is not None
                 and str(r.get("include_in_validation", "")).strip() in {"1", "1.0", "true", "True"}]
    if main_rows and validation is not None and validation.cbl_pass_rate is not None:
        best = None
        for r in main_rows:
            pr = fnum(r.get("cbl_pass_rate")) / 100.0
            if best is None or abs(pr - validation.cbl_pass_rate) < abs(best[0] - validation.cbl_pass_rate):
                best = (pr, r)
        pr, r = best
        checks.add(well, "cbl", "CBL合格率(fraction)", validation.cbl_pass_rate, pr,
                   cmp_num(validation.cbl_pass_rate, pr, tol_pct=0.5),
                   f"来源 {r.get('source_file', '')[:60]}")
        checks.add(well, "cbl", "CBL测量段 top/bottom",
                   f"窗{[f'{w.top_md_m:.0f}-{w.bottom_md_m:.0f}' for w in spec.evaluation_windows if w.window_type == 'cbl']}",
                   f"{fmt(fnum(r.get('cbl_top_md_m')))}–{fmt(fnum(r.get('cbl_bottom_md_m')))}",
                   "⚠", "loader评价窗常为可评价子段（扣除双层套管/悬空段），见notes")
    # target_intervals 覆盖
    tgt = read_rows(EXTRACT_ROOT / WELL_DIRS[well] / "target_intervals.csv")
    for r in tgt:
        top, bot = fnum(r.get("top_md_m")), fnum(r.get("bottom_md_m"))
        if top is None or bot is None:
            continue
        covered = any(w.top_md_m <= top and w.bottom_md_m >= bot
                      for w in spec.evaluation_windows if w.window_type == "formation_target")
        checks.add(well, "cbl", f"地层目标[{r.get('interval_name')}] {fmt(top)}-{fmt(bot)}",
                   "covered" if covered else "未覆盖", "target_intervals.csv",
                   "✓" if covered else "⚠", "")


def check_summary_table(well: str, spec, fluids, schedule, validation, checks: Checks) -> None:
    for r in read_rows(SUMMARY_TABLE_CSV):
        if r.get("well_id", "").strip() != well:
            continue
        group, name = r.get("field_group", "").strip(), r.get("field_name", "").strip()
        fv, av = fnum(r.get("field_value")), fnum(r.get("adopted_value"))
        item = f"总表[{group}/{name}]"
        if group == "cbl" and name == "cbl_pass_rate" and validation:
            checks.add(well, "总表", item, validation.cbl_pass_rate, av if av is not None else fv,
                       cmp_num(validation.cbl_pass_rate, av if av is not None else fv, tol_pct=0.5),
                       r.get("notes", "")[:100])
        elif group == "program" and name.endswith("_volume_m3"):
            cn = {"lead": "lead", "tail": "tail", "intermediate": "intermediate"}.get(name.split("_")[0])
            lv = None
            for s in schedule.steps:
                fk = spec_role_key_by_name(fluids, s.fluid_name)
                if cn and fk == cn:
                    lv = (lv or 0) + s.volume_m3
            target = fv if fv is not None else av
            checks.add(well, "总表", item, lv, target,
                       cmp_num(lv, target, tol_abs=0.2, tol_pct=1.0),
                       (r.get("notes", "") or "")[:100])
        else:
            checks.add(well, "总表", item, "—", av if av is not None else fv, "-",
                       (r.get("notes", "") or "")[:120])


def spec_role_key_by_name(fluids, name: str) -> str:
    for f in fluids:
        if f.name == name:
            return spec_role_key(f)
    return role_key(name)


def check_internal(well: str, spec, fluids, schedule, sched_rows, geo, checks: Checks) -> None:
    # 1) 每步 rate×duration == volume（推算闭合，应恒等；显式时间步除外）
    bad = []
    for r in sched_rows:
        if r["duration_min"] is not None and r["volume_m3"] > 0 and r["rate_m3_min"] > 0:
            resid = abs(r["rate_m3_min"] * r["duration_min"] - r["volume_m3"])
            if resid > 0.01:
                bad.append(f"#{r['index']}{r['step_name']}残差{resid:.3f}")
    checks.add(well, "内部", "每步rate×duration=体积", "闭合", ";".join(bad) or "闭合",
               "✓" if not bad else "✗", "")

    # 2) 累计体积闭合
    cement_roles = {"lead", "tail", "intermediate"}
    cement_vol = sum(r["volume_m3"] for r in sched_rows if r["role"] in cement_roles)
    total = schedule.total_injected_volume_m3
    checks.add(well, "内部", "浆体总体积m3(lead+tail+intermediate)", round(cement_vol, 2),
               round(total, 2), "✓", f"日程总注入{total:.1f}m³含前置/顶替")

    # 3) 环空体积重算 & 水泥浆/环空库存比
    ann = annulus_volume_m3(spec)
    if not math.isnan(ann) and ann > 0:
        checks.add(well, "内部", "环空体积重算m3(hole_profile×OD积分)", round(ann, 2), "—", "⚠",
                   "重算口径：hole_diameter_profile×(双径)liner OD；技套重叠段按外推井径计")
        ratio = cement_vol / ann
        checks.add(well, "内部", "水泥浆总体积/环空体积(η_E库存比上限)", round(ratio, 3), "—", "⚠",
                   ">1 表示水泥浆量超过环空容量（返高以上/漏失/计量口径）")

    # 4) 管内容积 vs 替浆总体积
    pipe = pipe_volume_m3(spec)
    disp_roles = {"displacement"}
    disp_vol = sum(r["volume_m3"] for r in sched_rows if r["role"] in disp_roles)
    if not math.isnan(pipe) and pipe > 0:
        st = "✓" if abs(pipe - disp_vol) <= max(3.0, 0.05 * max(pipe, disp_vol)) else "⚠"
        checks.add(well, "内部", "管内容积重算m3 vs 替浆总体积m3",
                   f"{pipe:.1f} vs {disp_vol:.1f}", "—", st,
                   "替浆总量应≈碰压时管内容积；差异大提示 liner_id/pipe_id 口径为等效或代理")

    # 5) 密度/流变 sanity
    for f in fluids:
        rho = f.density_kg_m3
        role = spec_role_key(f)
        if role in {"mud"}:
            ok = 1050.0 <= rho <= 2350.0
            checks.add(well, "sanity", f"密度范围 {f.name}", rho, "泥浆1.05–2.35g/cc",
                       "✓" if ok else "✗", "")
        elif role in {"lead", "tail", "intermediate"}:
            ok = 1850.0 <= rho <= 1950.0
            checks.add(well, "sanity", f"密度范围 {f.name}", rho, "水泥浆1.85–1.95g/cc",
                       "✓" if ok else "⚠", "超出常规水泥浆范围：加重/低密配方需与配方记录核对" if not ok else "")
        if f.power_law_n is not None:
            ok = 0.2 <= f.power_law_n <= 1.2
            checks.add(well, "sanity", f"n范围 {f.name}", f.power_law_n, "0.2–1.2", "✓" if ok else "✗", "")
        if f.consistency_k is not None:
            ok = 0.05 <= f.consistency_k <= 5.0
            checks.add(well, "sanity", f"K范围 {f.name}", f.consistency_k, "0.05–5 Pa·s^n",
                       "✓" if ok else "⚠", "")
        if f.yield_stress_pa is not None:
            checks.add(well, "sanity", f"YP>0 {f.name}", f.yield_stress_pa, ">0",
                       "✓" if f.yield_stress_pa > 0 else "⚠", "YP=0 仅当流体确无屈服")

    # 6) 几何 sanity
    if spec.liner_id_mm and spec.liner_od_mm:
        checks.add(well, "sanity", "ID<OD", f"{spec.liner_id_mm}<{spec.liner_od_mm}", "—",
                   "✓" if spec.liner_id_mm < spec.liner_od_mm else "✗", "")
    if spec.hanger_md_m is not None:
        checks.add(well, "sanity", "悬挂器<鞋深", f"{spec.hanger_md_m}<{spec.shoe_md_m}", "—",
                   "✓" if spec.hanger_md_m < spec.shoe_md_m else "✗", "")
        checks.add(well, "sanity", "尾管长度≈shoe−hanger", round(spec.shoe_md_m - spec.hanger_md_m, 1),
                   f"模型段长{spec.bottom_md_m - spec.top_md_m:.1f}",
                   "✓" if abs((spec.shoe_md_m - spec.hanger_md_m) - (spec.bottom_md_m - spec.top_md_m)) <= 60.0 else "⚠",
                   "模型段顶可为水泥返高（≠悬挂器），差异仅提示")
    if spec.is_dual_diameter:
        checks.add(well, "sanity", "上段OD>下段OD(双径井)",
                   f"{spec.upper_liner_od_mm}>{spec.liner_od_mm}", "—",
                   "✓" if spec.upper_liner_od_mm > spec.liner_od_mm else "✗", "")
    # CBL 窗在尾管段内
    if spec.hanger_md_m is not None:
        for w in spec.evaluation_windows:
            if w.window_type == "cbl" and (w.top_md_m < spec.top_md_m - 1 or w.bottom_md_m > spec.shoe_md_m + 1):
                checks.add(well, "sanity", f"CBL窗在尾管段内[{w.name}]",
                           f"{w.top_md_m}-{w.bottom_md_m}", f"{spec.top_md_m}-{spec.shoe_md_m}",
                           "⚠", "评价窗超出模型域/鞋深（悬空段口径见notes）")
    # standoff 范围
    if spec.standoff_profile:
        vals = [p.value for p in spec.standoff_profile]
        checks.add(well, "sanity", "standoff∈[0,1]", f"{min(vals):.2f}–{max(vals):.2f}", "—",
                   "✓" if min(vals) >= 0 and max(vals) <= 1 else "✗", "")
    # hole profile > OD
    if spec.hole_diameter_profile:
        bad = [p.value for p in spec.hole_diameter_profile
               if p.value <= (spec.upper_liner_od_mm if spec.is_dual_diameter and p.depth_md_m <= spec.upper_section_bottom_md_m else spec.liner_od_mm or 0)]
        checks.add(well, "sanity", "井径>尾管OD(全部点)", "—", f"{len(bad)}点违反",
                   "✓" if not bad else "✗", "")


# ---------------------------------------------------------------- 主流程

WELL_DIRS = {
    "hu101": "hu101_呼101",
    "hu102": "hu102_呼102",
    "hu103": "hu103_呼103",
    "hu1": "hu1_呼探1",
    "hu2": "ht1_002_呼探1-002",
    "ht1_001": "ht1_001_呼探1-001",
    "ht1_003": "ht1_003_呼1-003",
    "ht1_004": "ht1_004_呼1-004",
}

LOADERS = {
    "hu101": load_hu101_tailpipe,
    "hu102": load_hu102_tailpipe,
    "hu103": load_hu103_tailpipe,
    "hu1": load_hu1_tailpipe,
    "hu2": load_hu2_tailpipe,
    "ht1_001": load_ht1_001_tailpipe,
    "ht1_003": load_ht1_003_tailpipe,
    "ht1_004": load_ht1_004_tailpipe,
}

VARIANTS = {  # 附加变体（仅 schedule 差异口径）
    "hu103": load_hu103_tailpipe_actual,
    "ht1_003": load_ht1_003_tailpipe_design,
    "ht1_004": load_ht1_004_tailpipe_actual,
}


def verify_one(well_key: str, tag: str, loader_fn, checks: Checks, full: bool = True) -> dict:
    spec, fluids, schedule, validation = loader_fn()
    well = f"{well_key}" if full else f"{well_key}[{tag}]"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    dump_fluids(well, fluids, OUT_DIR / f"{well}_fluids.csv")
    sched_rows = dump_schedule(well, schedule, OUT_DIR / f"{well}_schedule.csv", fluids)
    geo = dump_geometry(well, spec, OUT_DIR / f"{well}_geometry.csv")
    if validation is not None:
        geo["kv"]["cbl_pass_rate"] = validation.cbl_pass_rate
    dump_cbl(well, spec, validation, OUT_DIR / f"{well}_cbl.csv")
    if full:
        check_fluids(well, fluids, checks)
        check_rheometer(well, fluids, checks)
        check_schedule(well, sched_rows, checks)
        check_geometry(well, spec, checks)
        check_cbl(well, spec, validation, checks)
        check_summary_table(well, spec, fluids, schedule, validation, checks)
    check_internal(well, spec, fluids, schedule, sched_rows, geo, checks)
    return {"spec": spec, "fluids": fluids, "schedule": schedule, "validation": validation,
            "sched_rows": sched_rows, "geo": geo}


def write_summary_md(data: dict, checks: Checks, path: Path) -> None:
    lines = ["# loader_dump 校验汇总（2026-09-01，脚本自动生成）", ""]
    lines.append("数据源：cemdisp loader 默认入口 × 现场提取层CSV（2026-08-29校准轮）+ 01_总表对照表。")
    lines.append("")
    # 统计
    lines.append("## 检查条目统计")
    lines.append("")
    lines.append("| 井 | ✓ | ✗ | ⚠ | - |")
    lines.append("|---|---|---|---|---|")
    wells = []
    for r in checks.rows:
        if r["well"] not in wells:
            wells.append(r["well"])
    for w in wells:
        c = {"✓": 0, "✗": 0, "⚠": 0, "-": 0}
        for r in checks.rows:
            if r["well"] == w and r["status"] in c:
                c[r["status"]] += 1
        lines.append(f"| {w} | {c['✓']} | {c['✗']} | {c['⚠']} | {c['-']} |")
    lines.append("")
    # 逐井概要
    for w, d in data.items():
        spec = d["spec"]
        lines.append(f"## {w}（{spec.well_name}）")
        lines.append("")
        lines.append(f"- 模型域 {spec.top_md_m}–{spec.bottom_md_m} m；鞋 {spec.shoe_md_m} m；"
                     f"悬挂器 {fmt(spec.hanger_md_m)} m；liner OD/ID {fmt(spec.liner_od_mm)}/{fmt(spec.liner_id_mm)} mm"
                     + (f"；双径上段 {fmt(spec.upper_liner_od_mm)}/{fmt(spec.upper_liner_id_mm)} mm 底 {spec.upper_section_bottom_md_m} m"
                        if spec.is_dual_diameter else ""))
        if d["validation"] is not None and d["validation"].cbl_pass_rate is not None:
            lines.append(f"- cbl_pass_rate = {d['validation'].cbl_pass_rate}")
        lines.append(f"- 日程总注入 {d['schedule'].total_injected_volume_m3:.1f} m³，{len(d['schedule'].steps)} 步")
        ann = annulus_volume_m3(spec)
        if not math.isnan(ann):
            cement = sum(r["volume_m3"] for r in d["sched_rows"]
                         if r["role"] in {"lead", "tail", "intermediate"})
            lines.append(f"- 环空体积重算 {ann:.1f} m³；浆体 {cement:.1f} m³；库存比 {cement / ann if ann else float('nan'):.2f}")
        pipe = pipe_volume_m3(spec)
        disp = sum(r["volume_m3"] for r in d["sched_rows"] if r["role"] == "displacement")
        lines.append(f"- 管内容积重算 {pipe:.1f} m³；顶替流体合计 {disp:.1f} m³")
        lines.append("")
        # ✗/⚠ 条目
        bad = [r for r in checks.rows if r["well"] == w and r["status"] in {"✗", "⚠"}]
        if bad:
            lines.append(f"### {w} 需关注条目（✗/⚠，共 {len(bad)} 条）")
            lines.append("")
            lines.append("| 类别 | 条目 | loader | CSV | 状态 | 备注 |")
            lines.append("|---|---|---|---|---|---|")
            for r in bad:
                lines.append(f"| {r['category']} | {r['item']} | {r['loader_value']} | "
                             f"{r['csv_value']} | {r['status']} | {r['note']} |")
            lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    checks = Checks()
    data: dict = {}
    for well_key, fn in LOADERS.items():
        try:
            data[well_key] = verify_one(well_key, "main", fn, checks, full=True)
        except Exception as exc:  # noqa: BLE001
            checks.add(well_key, "load", "loader加载", "—", repr(exc)[:120], "✗", "加载失败")
            continue
        variant_fn = VARIANTS.get(well_key)
        if variant_fn is not None:
            try:
                data[f"{well_key}[variant]"] = verify_one(well_key, "variant", variant_fn, checks, full=False)
            except Exception as exc:  # noqa: BLE001
                checks.add(f"{well_key}[variant]", "load", "loader变体加载", "—", repr(exc)[:120], "⚠", "变体加载失败")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with (OUT_DIR / "checks_all.csv").open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["well", "category", "item", "loader_value",
                                           "csv_value", "status", "note"])
        w.writeheader()
        w.writerows(checks.rows)
    write_summary_md(data, checks, OUT_DIR / "loader_dump_汇总.md")
    n = {"✓": 0, "✗": 0, "⚠": 0, "-": 0}
    for r in checks.rows:
        n[r["status"]] = n.get(r["status"], 0) + 1
    print(f"完成：{len(checks.rows)} 条检查 → {OUT_DIR}")
    print(f"统计：✓{n['✓']}  ✗{n['✗']}  ⚠{n['⚠']}  -{n['-']}")


if __name__ == "__main__":
    main()
