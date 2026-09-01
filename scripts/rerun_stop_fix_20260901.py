"""F2 停止标志修复（2026-09-01）后续量化：+600s 尾窗去除 与 环空域口径 两个 2D 效应。

背景：
- casing_flow.py 尾缘重力修正已改为按"后继流体(pusher) vs 尾浆"配对 + max 约束 +
  "不超过泵注结束"封顶；8 井 1D stop 数值与修复前一致（Δ=0），334 测试全过。1D 侧已闭环。
- 待量化的是两个 2D 侧效应：
  (1) +600s 尾窗去除：旧基线（results/最终基线_2026-08-29，tt=stop+600）vs 修复后 tt=stop；
  (2) 环空域口径：A"当前域"（2D 域含技套重叠段）vs B"仅裸眼域"（2D 域顶截到技套鞋）。

变体定义（每井 1D cr / inlet / tt 公共，只换 2D 求解器收到的 WellSpec）：
  A 当前域  : well 原样，tt = 修复后 stop（cr.cement_end_time_s 优先，None 回退 pumping_end_time_s）
  B 仅裸眼域: dataclasses.replace 构造 B WellSpec——top_md_m=技套鞋；hole_diameter_profile 截断到
              md≥技套鞋并在技套鞋处补线性插值端点；evaluation_windows 逐窗裁剪到 [技套鞋, 鞋深]，
              完全在域外的窗丢弃；悬挂器位于技套鞋之上（域外），按 WellSpec 校验约束（hanger∈[top,bottom]）
              夹取 hanger_md_m=域顶；standoff/inclination/liner_od 剖面不动（求解器 np.interp 自动夹取）。
  口径注意：B 复用与 A 同源的入口时线（重叠段环空输运时滞未重算），是域截断效应的近似量化。

输出 results/停止标志修复与裸眼域对比_2026-09-01/：
  对比.csv / 域几何量化.csv / 摘要.md / 单井结果/*.json（逐井落盘，断点续跑）

运行（nz=250 生产网格，16 次 2D 重跑约 1-2 小时）：
  PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/rerun_stop_fix_20260901.py
"""
from __future__ import annotations

import csv
import dataclasses
import json
import time
import traceback
from pathlib import Path
from typing import cast

import numpy as np

from cemdisp.data.fluid_spec import FluidRole
from cemdisp.data.well_spec import DepthValuePoint, EvaluationWindow, WellSpec
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver
import cemdisp.data.loaders as L

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT = PROJECT_ROOT / "results" / "停止标志修复与裸眼域对比_2026-09-01"
RESULTS_DIR = OUT / "单井结果"
NZ = 250
BASELINE_CSV = PROJECT_ROOT / "results" / "最终基线_2026-08-29" / "cfl_on" / "汇总.csv"

# corrected（adopted）口径，照抄 scripts/rerun_all_wells_corrected.py:33-39（RR）
CORRECTED_KW = dict(
    dispersion_dt_scale=1.0,   # M1: 弥散 dt 归一（不再随 dt 缩放）
    enable_yield_gate=True,    # M3: 屈服门槛
    enable_regime_split=True,  # M2: 局部流态修正（层流元 R=1，中性）
    enable_local_i3=True,      # I3: 浮力弥散通量局部化
    e_clip_max=0.90,           # M4: 效率截断上限 0.55 -> 0.90
)

WELLS = [
    ("hu101", L.load_hu101_tailpipe), ("hu102", L.load_hu102_tailpipe),
    ("hu103", L.load_hu103_tailpipe), ("hu1", L.load_hu1_tailpipe),
    ("hu2", L.load_hu2_tailpipe), ("ht1_001", L.load_ht1_001_tailpipe),
    ("ht1_003", L.load_ht1_003_tailpipe), ("ht1_004", L.load_ht1_004_tailpipe),
]

# ---------------------------------------------------------------------------
# 技套鞋 / 上层套管鞋深度（B 变体的域顶）。逐井来源：
# - hu101  : 273.1 技套鞋 5699.8m（loader 注释 L70/L178；CBL"单层套管可评价段"顶 5699.8）
# - hu102  : 219.1 技套鞋 7119.80m（6823.10–7119.80m 双层套管段底=技套鞋，loader L74；裸眼实测自 7120）
# - hu103  : 273.1 技套鞋 5750m（现场提取包 well_geometry.csv openhole_top_md=5750"上层273.1mm套管鞋"；
#            10053.txt 四开"技术套管273.1mm，下深5750m"；extraction_notes 重叠段~200m 在 273.1 技套内）
# - hu1    : 219.1 技术尾管鞋 5693.89m，裸眼段顶取 HU1_OPEN_HOLE_TOP_MD_M=5694.1（井径首测点）
# - hu2    : 219.1 技套鞋 5630m（loader 注释 L10/L57/L207）
# - ht1_001: 273.1 技套鞋 5670m（重叠段 5469.711–5670m 在技套内，loader L96；裸眼实测自 5670）
# - ht1_003: HT1_003_CASING_SHOE_MD_M = 5568.0
# - ht1_004: HT1_004_CASING_SHOE_MD_M = 5578.0
TECH_SHOE_MD_M = {
    "hu101": 5699.8,
    "hu102": 7119.80,
    "hu103": 5750.0,
    "hu1": 5694.1,
    "hu2": 5630.0,
    "ht1_001": 5670.0,
    "ht1_003": 5568.0,
    "ht1_004": 5578.0,
}

# 重叠段技套内径（"技套ID口径"体积用）。来源同各 loader 常量；
# hu103 loader 未存 273.1 技套 ID，取计算值 245.42（273.1−2×13.84，同 ht1_001 口径，proxy）。
TECH_CASING_ID_MM = {
    "hu101": 245.37,   # HU101_TECH_CASING_EQUIV_ID_MM
    "hu102": 193.70,   # HU102_CASING_ID_MM
    "hu103": 245.42,   # 计算值 proxy（273.1−2×13.84）
    "hu1": None,       # 两段式，见下方 OVERLAP_SEGMENTS_SPECIAL
    "hu2": 193.7,      # 219.1−2×12.7
    "ht1_001": 245.42, # HT1_001_CASING_ID_MM
    "ht1_003": 245.37, # HT1_003_CASING_ID_MM
    "ht1_004": 245.37, # HT1_004_CASING_INNER_DIAMETER_MM
}

# 重叠段内尾管外径（技套 ID 口径体积用）。注意：hu101/ht1_001 模型为单一 139.7 外径表示
# （无 upper_section 字段），ht1_003 有 liner_od_profile 但重叠段同样属 168.3 上段——
# 物理口径按管串表取 168.3，不能回落到模型的 139.7 表示外径。
OVERLAP_LINER_OD_MM = {
    "hu101": 168.3,   # 上段 168.3mm 尾管（5402.9–6796.3m 段覆盖重叠段）
    "hu102": 139.7,   # 139.7mm 尾管在 219.1 技套内
    "hu103": 168.3,   # 上段 168.3mm 尾管（5546–7330.7m 段覆盖重叠段）
    "hu1": 139.7,     # 两段式特例见 OVERLAP_SEGMENTS_SPECIAL
    "hu2": 139.7,
    "ht1_001": 168.3, # 上段 168.3mm 尾管（5469.7–7174.9m 段覆盖重叠段）
    "ht1_003": 168.3, # 上段 168.3mm 尾管（5307.5–7089.6m 段覆盖重叠段）
    "ht1_004": 168.3, # 上段 168.3mm 尾管（5243.2–7376.7m 段覆盖重叠段）
}

# hu1 重叠段两段式：139.7 尾管在 273.05 技套内（ID 245.37，3523.27–3623.27），
# 再在 219.1 技术尾管内（ID 193.7，3623.27–技套鞋）。
OVERLAP_SEGMENTS_SPECIAL = {
    "hu1": lambda top, shoe: [
        (top, 3623.27, 245.37, 139.7),
        (3623.27, shoe, 193.7, 139.7),
    ],
}


# ---------------------------------------------------------------------------
def load_old_baseline_eta_e() -> dict[str, float]:
    """读 results/最终基线_2026-08-29/cfl_on/汇总.csv 的 eta_E（旧口径 tt=stop+600）。"""
    eta: dict[str, float] = {}
    with BASELINE_CSV.open("r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            eta[row["well"].strip()] = float(row["eta_E"])
    return eta


def _stop_t(cr) -> float:
    """修复后停止时刻：优先 cement_end_time_s（尾缘≡后继流体前缘，同刻修正），
    None 时回退泵注结束时刻。"""
    if cr.cement_end_time_s is not None:
        return float(cr.cement_end_time_s)
    return float(cr.pumping_end_time_s)


def build_openhole_well(well: WellSpec, tech_shoe: float) -> WellSpec:
    """构造 B 变体"仅裸眼域" WellSpec（frozen dataclass → dataclasses.replace）。"""
    if not (well.top_md_m < tech_shoe < well.bottom_md_m):
        raise ValueError(f"技套鞋 {tech_shoe} 不在井段 ({well.top_md_m}, {well.bottom_md_m}) 内")
    return dataclasses.replace(
        well,
        top_md_m=tech_shoe,
        # 悬挂器在技套鞋之上（域外）；WellSpec 校验要求 hanger∈[top,bottom]，夹取到域顶（元数据口径）
        hanger_md_m=max(well.hanger_md_m, tech_shoe) if well.hanger_md_m is not None else None,
        hole_diameter_profile=_truncate_hole_profile(well.hole_diameter_profile, tech_shoe),
        evaluation_windows=_clip_windows(well.evaluation_windows, tech_shoe, well.bottom_md_m),
    )


def _hole_profile_value(profile: tuple[DepthValuePoint, ...], md: float) -> float:
    """模型井径剖面在 md 处的线性插值（无下方点时 np.interp 夹取到首点值，同求解器边界行为）。"""
    mds = np.array([p.depth_md_m for p in profile], dtype=float)
    vals = np.array([p.value for p in profile], dtype=float)
    return float(np.interp(md, mds, vals))


def _truncate_hole_profile(profile, tech_shoe: float):
    """截断为 md≥技套鞋，并在技套鞋处补一个线性插值端点（域顶点必须存在）。"""
    kept = [p for p in profile if p.depth_md_m >= tech_shoe]
    if not kept:
        raise ValueError("井径剖面截断后为空")
    if kept[0].depth_md_m > tech_shoe:
        kept.insert(0, DepthValuePoint(depth_md_m=tech_shoe, value=_hole_profile_value(profile, tech_shoe)))
    return tuple(kept)


def _clip_windows(windows, tech_shoe: float, bottom: float):
    """评价窗逐个裁剪到 [技套鞋, 鞋深]，完全在域外的窗丢弃。"""
    out = []
    for w in windows:
        top2 = max(w.top_md_m, tech_shoe)
        bot2 = min(w.bottom_md_m, bottom)
        if top2 >= bot2:
            continue
        if top2 == w.top_md_m and bot2 == w.bottom_md_m:
            out.append(w)
        else:
            out.append(EvaluationWindow(name=w.name, top_md_m=top2, bottom_md_m=bot2, window_type=w.window_type))
    return tuple(out)


# ---------------------------------------------------------------------------
def _model_od_arrays(well: WellSpec):
    """模型外径剖面：liner_od_profile 存在则用之，否则全域单值 liner_od_mm（同 2D 求解器口径）。"""
    if well.liner_od_profile:
        md = np.array([p.depth_md_m for p in well.liner_od_profile], dtype=float)
        v = np.array([p.value for p in well.liner_od_profile], dtype=float)
        return md, v
    return None, float(well.liner_od_mm or 0.0)


def _annulus_volume_m3(well: WellSpec, md_from: float, md_to: float) -> float:
    """环空体积（模型口径）：∫ π/4·(D_hole²−D_od²) dz，井径/外径取模型剖面插值。

    与 2D 求解器几何一致：井径用 hole_diameter_profile（168.3 段已保面积等效到
    139.7 参考外径的井，等效值直接进积分），外径用 liner_od_profile 或单值 liner_od_mm。
    """
    hole_md = np.array([p.depth_md_m for p in well.hole_diameter_profile], dtype=float)
    hole_v = np.array([p.value for p in well.hole_diameter_profile], dtype=float)
    od_md, od_v = _model_od_arrays(well)
    n = 4000
    zs = np.linspace(md_from, md_to, n + 1)
    dz = (md_to - md_from) / n
    hole = np.interp(zs, hole_md, hole_v)
    od = od_v if isinstance(od_v, float) else np.interp(zs, od_md, od_v)
    area_m2 = np.pi / 4.0 * ((hole / 1000.0) ** 2 - (od / 1000.0) ** 2)
    return float(np.sum(area_m2) * dz)


def _overlap_tech_id_volume_m3(key: str, well: WellSpec, tech_shoe: float) -> float | None:
    """技套 ID 口径的重合段体积：Σ π/4·(D_tech_id²−D_liner_od²)·dz。

    返回 None 表示该井无法计算（缺技套内径）。
    """
    top_md = well.top_md_m
    if key in OVERLAP_SEGMENTS_SPECIAL:
        segs = OVERLAP_SEGMENTS_SPECIAL[key](top_md, tech_shoe)
    else:
        tech_id = TECH_CASING_ID_MM.get(key)
        if tech_id is None:
            return None
        od = OVERLAP_LINER_OD_MM.get(key, float(well.liner_od_mm or 0.0))
        segs = [(top_md, tech_shoe, tech_id, od)]
    vol = 0.0
    for a, b, tid, lod in segs:
        if b > a:
            vol += np.pi / 4.0 * ((tid / 1000.0) ** 2 - (lod / 1000.0) ** 2) * (b - a)
    return float(vol)


def _cement_volume_m3(fluids, schedule) -> float:
    """水泥浆总体积：泵注程序中 LEAD/INTERMEDIATE/TAIL 角色步骤的体积和（排量>0）。"""
    cement_roles = {FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL}
    role_by_name = {f.name: f.role for f in fluids}
    return float(sum(
        s.volume_m3 for s in schedule.steps
        if s.rate_m3_min > 0 and role_by_name.get(s.fluid_name) in cement_roles
    ))


# ---------------------------------------------------------------------------
def run_variant(loader, well_b: WellSpec | None):
    """复刻 RR run_one：loader → 1D cr → 修复后 stop → inlet → 2D（A 原 well / B 截域 well）。"""
    t0 = time.perf_counter()
    well, fluids, schedule, _ = loader()
    cr = CasingFlowSolver(enable_gravity=True).run(well, fluids, schedule)
    tt = _stop_t(cr)
    inlet = build_coupled_annulus_inlet_provider(
        cr, CasingFlowSolver(enable_gravity=True), fluids, split_cement_phases=True
    )
    solver_well = well if well_b is None else well_b
    res = AnnulusD2DGASolver(total_t=tt, nz=NZ, enable_cfl_adaptive=True, **CORRECTED_KW).run(
        solver_well, fluids, inlet
    )
    fr = cast(dict, res.summary["最终结果"])
    return {
        "stop_s": tt,
        "eta_E": float(fr["全井段最终有效顶替效率"]),
        "eta_N": float(fr.get("窄四分位效率", float("nan"))),
        "mixing": float(fr["最终混浆指数"]),
        "channeling": float(fr["最终窜槽指数"]),
        "instability": float(fr["最终失稳指数"]),
        "cement_occ": float(fr["最终水泥浆占据率"]),
        "wall_frac": float(np.mean(res.wall_field)) if res.wall_field is not None else float("nan"),
        "elapsed_s": round(time.perf_counter() - t0, 1),
    }


def _save_json(label: str, domain: str, row: dict):
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    path = RESULTS_DIR / f"{label}_{domain}.json"
    path.write_text(json.dumps(row, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"  落盘: {path.name}", flush=True)


def _load_json(label: str, domain: str) -> dict | None:
    path = RESULTS_DIR / f"{label}_{domain}.json"
    if path.exists():
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            return None
    return None


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    baseline = load_old_baseline_eta_e()
    print(f"旧基线对照: {BASELINE_CSV.name}（{len(baseline)} 井）", flush=True)

    rows: list[dict] = []
    geom_rows: list[dict] = []
    errors: list[str] = []
    started = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"=== 开始：{started} | 共 {len(WELLS)} 井 × 2 域变体（nz={NZ}）===", flush=True)

    for idx, (label, loader) in enumerate(WELLS, 1):
        print(f"\n[{idx}/{len(WELLS)}] === {label} ===", flush=True)
        well = fluids = schedule = None
        try:
            well, fluids, schedule, _ = loader()
        except Exception as e:
            errors.append(f"{label}: loader 失败 {e}")
            traceback.print_exc()
            continue

        tech_shoe = TECH_SHOE_MD_M.get(label)
        try:
            well_b = build_openhole_well(well, tech_shoe) if tech_shoe is not None else None
        except Exception as e:
            errors.append(f"{label}: B 变体 WellSpec 构造失败 {e}")
            traceback.print_exc()
            well_b = None

        # ---- 几何量化（8 井全出，无论 B 是否跑）----
        try:
            cement_vol = _cement_volume_m3(fluids, schedule)
            top, bottom = well.top_md_m, well.bottom_md_m
            vol_overlap_model = _annulus_volume_m3(well, top, tech_shoe) if tech_shoe else 0.0
            vol_openhole = _annulus_volume_m3(well, tech_shoe, bottom) if tech_shoe else float("nan")
            vol_overlap_techid = (
                _overlap_tech_id_volume_m3(label, well, tech_shoe) if tech_shoe else None
            )
            vol_cur = (vol_overlap_model or 0.0) + (vol_openhole if vol_openhole == vol_openhole else 0.0)
            geom_rows.append({
                "well": label,
                "top_md": top,
                "hanger_md": well.hanger_md_m,
                "tech_shoe_md": tech_shoe if tech_shoe is not None else "",
                "shoe_md": well.shoe_md_m,
                "域总长_m": round(bottom - top, 2),
                "裸眼段体积_m3": round(vol_openhole, 3) if tech_shoe else "",
                "重合段体积_m3_模型口径": round(vol_overlap_model, 3) if tech_shoe else "",
                "重合段体积_m3_技套ID口径": (
                    round(vol_overlap_techid, 3) if vol_overlap_techid is not None else ""
                ),
                "水泥浆总体积_m3": round(cement_vol, 2),
                "库存比_当前域": round(cement_vol / vol_cur, 4) if vol_cur > 0 else "",
                "库存比_裸眼域": round(cement_vol / vol_openhole, 4) if vol_openhole and vol_openhole > 0 else "",
            })
            print(f"  几何: 技套鞋={tech_shoe} 裸眼段={geom_rows[-1]['裸眼段体积_m3']}m³ "
                  f"重合段(模型)={geom_rows[-1]['重合段体积_m3_模型口径']}m³ 水泥浆={geom_rows[-1]['水泥浆总体积_m3']}m³",
                  flush=True)
        except Exception as e:
            errors.append(f"{label}: 几何量化失败 {e}")
            traceback.print_exc()

        # ---- A 变体（当前域，tt=修复后 stop）----
        row_a = _load_json(label, "A")
        if row_a is None:
            print(f"  A(当前域) 开始…", flush=True)
            try:
                row_a = run_variant(loader, None)
                row_a["domain"] = "A_当前域"
                _save_json(label, "A", row_a)
            except Exception as e:
                errors.append(f"{label}: A 变体失败 {e}")
                traceback.print_exc()
                row_a = None
        if row_a is not None:
            old = baseline.get(label, float("nan"))
            rec = {
                "well": label, "domain": "A_当前域", "stop_s": row_a["stop_s"],
                "eta_E": row_a["eta_E"], "eta_N": row_a["eta_N"], "mixing": row_a["mixing"],
                "channeling": row_a["channeling"], "instability": row_a["instability"],
                "cement_occ": row_a["cement_occ"], "wall_frac": row_a["wall_frac"],
                "eta_E_旧基线": old, "delta_vs_旧基线": row_a["eta_E"] - old,
                "elapsed_s": row_a["elapsed_s"],
            }
            rows.append(rec)
            print(f"  A(当前域): tt={row_a['stop_s']:.1f}s eta_E={row_a['eta_E']:.4f} "
                  f"(旧基线 {old:.4f}, Δ{rec['delta_vs_旧基线']:+.4f}) "
                  f"mix={row_a['mixing']:.4f} occ={row_a['cement_occ']:.4f} ({row_a['elapsed_s']}s)", flush=True)

        # ---- B 变体（仅裸眼域）----
        if well_b is None:
            print(f"  B(仅裸眼域): 跳过（技套鞋不可得或构造失败）", flush=True)
            continue
        row_b = _load_json(label, "B")
        if row_b is None:
            print(f"  B(仅裸眼域) 开始…", flush=True)
            try:
                row_b = run_variant(loader, well_b)
                row_b["domain"] = "B_仅裸眼域"
                _save_json(label, "B", row_b)
            except Exception as e:
                errors.append(f"{label}: B 变体失败 {e}")
                traceback.print_exc()
                row_b = None
        if row_b is not None:
            old = baseline.get(label, float("nan"))
            rec = {
                "well": label, "domain": "B_仅裸眼域", "stop_s": row_b["stop_s"],
                "eta_E": row_b["eta_E"], "eta_N": row_b["eta_N"], "mixing": row_b["mixing"],
                "channeling": row_b["channeling"], "instability": row_b["instability"],
                "cement_occ": row_b["cement_occ"], "wall_frac": row_b["wall_frac"],
                "eta_E_旧基线": old, "delta_vs_旧基线": row_b["eta_E"] - old,
                "elapsed_s": row_b["elapsed_s"],
            }
            rows.append(rec)
            delta_ba = (row_b["eta_E"] - row_a["eta_E"]) if row_a is not None else float("nan")
            print(f"  B(仅裸眼域): tt={row_b['stop_s']:.1f}s eta_E={row_b['eta_E']:.4f} "
                  f"(vs A {delta_ba:+.4f}) mix={row_b['mixing']:.4f} occ={row_b['cement_occ']:.4f} "
                  f"({row_b['elapsed_s']}s)", flush=True)

    # ---- 写 对比.csv ----
    csv_path = OUT / "对比.csv"
    keys = ["well", "domain", "stop_s", "eta_E", "eta_N", "mixing", "channeling",
            "instability", "cement_occ", "wall_frac", "eta_E_旧基线", "delta_vs_旧基线", "elapsed_s"]
    with csv_path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in keys})
    print(f"\n对比表: {csv_path}（{len(rows)} 行数据）", flush=True)

    # ---- 写 域几何量化.csv（8 井全出）----
    geom_path = OUT / "域几何量化.csv"
    gkeys = ["well", "top_md", "hanger_md", "tech_shoe_md", "shoe_md", "域总长_m",
             "裸眼段体积_m3", "重合段体积_m3_模型口径", "重合段体积_m3_技套ID口径",
             "水泥浆总体积_m3", "库存比_当前域", "库存比_裸眼域"]
    with geom_path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=gkeys)
        w.writeheader()
        for r in geom_rows:
            w.writerow({k: r.get(k, "") for k in gkeys})
    print(f"几何量化表: {geom_path}（{len(geom_rows)} 行数据）", flush=True)

    # ---- 写 摘要.md ----
    write_summary(rows, geom_rows, errors, started)

    # ---- 校验 ----
    print("\n" + "=" * 78)
    n_a = sum(1 for r in rows if r["domain"] == "A_当前域")
    n_b = sum(1 for r in rows if r["domain"] == "B_仅裸眼域")
    print(f"完成: 对比 {len(rows)} 行（A={n_a}, B={n_b}）| 几何 {len(geom_rows)} 行 | 异常 {len(errors)} 项")
    for e in errors:
        print(f"  [异常] {e}")
    print(f"输出目录: {OUT}")


def write_summary(rows, geom_rows, errors, started):
    """摘要.md：两张表 + 结论要点。"""
    by = {(r["well"], r["domain"]): r for r in rows}
    wells = [w for w, _ in WELLS]

    def fmt(v, nd=4):
        return f"{v:.{nd}f}" if isinstance(v, (int, float)) and v == v else "—"

    lines = [
        "# 停止标志修复与裸眼域对比摘要（2026-09-01）",
        "",
        f"- 生成：{started}，脚本 `scripts/rerun_stop_fix_20260901.py`，生产网格 nz={NZ}，CFL 自适应，CORRECTED 口径",
        "  （M1 弥散 dt 归一 + M2 流态修正 + M3 屈服门槛 + I3 局部化 + M4 e_clip_max=0.90）。",
        "- 停止口径：修复后 `tt = cr.cement_end_time_s`（尾浆全部入库时刻，None 回退泵注结束），无 +600s 尾窗。",
        "- 旧基线对照：`results/最终基线_2026-08-29/cfl_on/汇总.csv` 的 eta_E（旧口径 tt=stop+600）。",
        "- A=当前域（2D 域含技套重叠段）；B=仅裸眼域（2D 域顶截到技套鞋，重叠段剔除）。",
        "- 口径注意：B 复用与 A 同源的入口时线（重叠段环空输运时滞未重算），是域截断效应的近似量化。",
        "- 重合段体积'技套ID口径'：按技套内径×重叠段内尾管物理外径（hu101/ht1_001/ht1_003 等双径井取 168.3，",
        "  非模型的 139.7 表示外径）；hu103 技套内径为计算值 proxy（273.1−2×13.84=245.42，loader 未存）；",
        "  hu1 重叠段两段式（273.05 技套 ID 245.37 → 219.1 技术尾管 ID 193.7）。裸眼段/重合段(模型口径)体积",
        "  按模型井径剖面+模型外径积分（168.3 段保面积等效到 139.7 参考外径），与 2D 求解器几何一致。",
        "  库存比 = 水泥浆总体积 ÷ 对应域环空体积（当前域=重合段+裸眼段；裸眼域=仅裸眼段）。",
        "",
        "## 表1 尾窗去除效应：A 变体（tt=stop）vs 旧基线（tt=stop+600）",
        "",
        "| 井 | stop_s | eta_E(修复后) | eta_E(旧基线) | Δ |",
        "|---|---|---|---|---|",
    ]
    dA = []
    for w in wells:
        r = by.get((w, "A_当前域"))
        if not r:
            lines.append(f"| {w} | — | — | — | — |")
            continue
        lines.append(f"| {w} | {r['stop_s']:.1f} | {fmt(r['eta_E'])} | {fmt(r['eta_E_旧基线'])} | {r['delta_vs_旧基线']:+.4f} |")
        dA.append(r["delta_vs_旧基线"])
    if dA:
        lines += ["", f"A vs 旧基线：均值 {np.mean(dA):+.4f}，范围 [{np.min(dA):+.4f}, {np.max(dA):+.4f}] pp。"]

    lines += [
        "",
        "## 表2 域口径效应：B 变体（仅裸眼域）vs A（当前域）",
        "",
        "| 井 | eta_E(A 当前域) | eta_E(B 仅裸眼域) | Δ(B−A) |",
        "|---|---|---|---|",
    ]
    dB = []
    for w in wells:
        ra, rb = by.get((w, "A_当前域")), by.get((w, "B_仅裸眼域"))
        if not ra or not rb:
            lines.append(f"| {w} | {fmt(ra['eta_E']) if ra else '—'} | {fmt(rb['eta_E']) if rb else '—'} | 跳过 |")
            continue
        d = rb["eta_E"] - ra["eta_E"]
        dB.append(d)
        lines.append(f"| {w} | {fmt(ra['eta_E'])} | {fmt(rb['eta_E'])} | {d:+.4f} |")
    if dB:
        lines += ["", f"B vs A：均值 {np.mean(dB):+.4f}，范围 [{np.min(dB):+.4f}, {np.max(dB):+.4f}] pp。"]

    lines += [
        "",
        "## 表3 域几何量化（8 井全出）",
        "",
        "| 井 | 域顶 | 悬挂器 | 技套鞋 | 鞋深 | 域总长_m | 裸眼段体积_m3 | 重合段_m3(模型) | 重合段_m3(技套ID) | 水泥浆_m3 | 库存比(当前域) | 库存比(裸眼域) |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for g in geom_rows:
        lines.append(
            f"| {g['well']} | {g['top_md']:.1f} | {g['hanger_md']:.1f} | {g['tech_shoe_md']} | {g['shoe_md']:.1f} "
            f"| {g['域总长_m']} | {g['裸眼段体积_m3']} | {g['重合段体积_m3_模型口径']} "
            f"| {g['重合段体积_m3_技套ID口径']} | {g['水泥浆总体积_m3']} | {g['库存比_当前域']} | {g['库存比_裸眼域']} |"
        )

    lines += ["", "## 结论要点", ""]
    if dA:
        lines.append(f"1. 尾窗去除（tt=stop+600→stop）对 eta_E 的生产网格实测影响：均值 {np.mean(dA):+.4f} pp"
                     f"（范围 {np.min(dA):+.4f}~{np.max(dA):+.4f}），逐井见表1。")
    else:
        lines.append("1. 尾窗去除效应：无可用数据。")
    if dB:
        lines.append(f"2. 域口径（剔除技套重叠段）对 eta_E 的影响：均值 {np.mean(dB):+.4f} pp"
                     f"（范围 {np.min(dB):+.4f}~{np.max(dB):+.4f}），逐井见表2；几何成因见表3（重合段占比与库存比变化）。")
    else:
        lines.append("2. 域口径效应：无可用数据。")
    if errors:
        lines += ["", "## 异常 / 跳过", ""]
        lines += [f"- {e}" for e in errors]

    path = OUT / "摘要.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"摘要: {path}", flush=True)


if __name__ == "__main__":
    main()
