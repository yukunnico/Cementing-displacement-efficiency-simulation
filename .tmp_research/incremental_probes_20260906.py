"""增量取证（2026-09-06）：η_E 饱和结构分解 + 弥散/居中度消融。

实验1：8 井终跑 2D 场分解——域内水泥 / 隔离液 / 残泥体积各占多少，
       检验 η_E≈1 时"1−η_E"由什么构成（残泥=效率真损失；隔离液=口径缺口），
       并给出窄四分位的同口径分解与残泥轴向三段分布。
实验2：hu103/ht1_003 当前生产口径（纯默认 nz250，与胶塞语义终跑一致）下
       消融方位角弥散 + standoff ±0.1，检验
       （a）方位角弥散常数 0.015 是否是窄边被人工填满的来源；
       （b）当前口径下居中度敏感性是否存在。

只读仓库、不改动 cemdisp 包；结果落 results/_增量取证_2026-09-06/。
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from scripts._mass_balance_diag_20260902 import WELLS, build_case  # noqa: E402
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver, _trapez2d  # noqa: E402
from cemdisp.data.well_spec import DepthValuePoint  # noqa: E402

OUT = PROJECT_ROOT / "results" / "_增量取证_2026-09-06"
OUT.mkdir(parents=True, exist_ok=True)
FINAL_DIR = PROJECT_ROOT / "results" / "胶塞语义修复后终跑_2026-09-03"
ABLATION_WELLS = ("hu103", "ht1_003")


def _rebuild_geom(name: str):
    """按终跑口径重建几何（nz=250, ny=40, 纯默认 e_clip=0.55）。"""
    well_spec, _f, _s, _p, _stop, _v = build_case(WELLS[name])
    geom = AnnulusD2DGASolver(nz=250, ny=40)._build_geom(well_spec)
    return geom


def exp1_field_decomposition() -> list[dict]:
    """实验1：读 8 井终跑 npz 场数据，做水泥/隔离液/残泥占据分解。"""
    rows = []
    for name in WELLS:
        npz_path = FINAL_DIR / f"{name}_2D场数据.npz"
        if not npz_path.exists():
            rows.append({"井": name, "状态": "npz缺失"})
            continue
        data = np.load(npz_path)
        cement = np.clip(data["cement_final"], 0.0, 1.0)   # (ny, nz)
        spacer = np.clip(data["spacer_final"], 0.0, 1.0)
        geom = _rebuild_geom(name)
        b = geom["b"]
        y = geom["y"]
        half_v = _trapez2d(b, geom)
        mud = np.clip(1.0 - cement - spacer, 0.0, 1.0)

        eta_E = _trapez2d(b * cement, geom) / half_v
        spacer_occ = _trapez2d(b * spacer, geom) / half_v
        mud_occ = _trapez2d(b * mud, geom) / half_v

        # 窄四分位同口径分解（与 _narrow_quarter_efficiency 一致：最后 ny//4 行）
        nq = max(1, cement.shape[0] // 4)
        b_q = b[-nq:, :]
        geom_q = {**geom, "b": b_q, "y": y[-nq:]}
        denom_q = _trapez2d(b_q, geom_q)
        eta_N = _trapez2d(b_q * cement[-nq:, :], geom_q) / denom_q
        mud_occ_N = _trapez2d(b_q * mud[-nq:, :], geom_q) / denom_q
        spacer_occ_N = _trapez2d(b_q * spacer[-nq:, :], geom_q) / denom_q

        # 残泥轴向三段分布（s 从鞋口=0 到域顶=s_max；域顶段是最难顶替的上部）
        nz_len = cement.shape[1]
        segs = {
            "域顶段": slice(0, nz_len // 3),
            "域中段": slice(nz_len // 3, 2 * nz_len // 3),
            "域底段": slice(2 * nz_len // 3, None),
        }
        mud_by_seg = {}
        for sname, seg in segs.items():
            b_seg = b[:, seg]
            g_seg = {**geom, "b": b_seg, "s": geom["s"][seg]}
            mud_seg = np.clip(1.0 - cement[:, seg] - spacer[:, seg], 0.0, 1.0)
            mud_by_seg[sname] = round(
                float(_trapez2d(b_seg * mud_seg, g_seg) / max(_trapez2d(b_seg, g_seg), 1e-12)), 4)

        rows.append({
            "井": name,
            "η_E": round(float(eta_E), 4),
            "隔离液占据率": round(float(spacer_occ), 4),
            "残泥占据率": round(float(mud_occ), 4),
            "η_N": round(float(eta_N), 4),
            "窄边残泥占据": round(float(mud_occ_N), 4),
            "窄边隔离液占据": round(float(spacer_occ_N), 4),
            **mud_by_seg,
        })
        print(f"[exp1] {name}: " + json.dumps(rows[-1], ensure_ascii=False), flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "场分解_8井.csv", index=False, encoding="utf-8-sig")
    return rows


def exp2_ablation() -> list[dict]:
    """实验2：方位角弥散消融 + standoff ±0.1 敏感性（当前生产口径）。"""
    rows = []
    for name in ABLATION_WELLS:
        well_spec, fluids, schedule, provider, stop, _v = build_case(WELLS[name])
        variants: dict[str, dict] = {
            "生产口径": dict(),
            "无方位角弥散": dict(dispersion_azimuthal=0.0),
            "无弥散": dict(dispersion_axial=0.0, dispersion_azimuthal=0.0),
            "双倍方位角弥散": dict(dispersion_azimuthal=0.030),
        }
        # standoff ±0.1：改设计剖面（clip 到 [0.2,1.0]）
        for delta in (-0.1, 0.1):
            new_pts = tuple(
                DepthValuePoint(depth_md_m=p.depth_md_m,
                                value=float(np.clip(p.value + delta, 0.2, 1.0)))
                for p in well_spec.standoff_profile
            )
            ws2 = dataclasses_replace_standoff(well_spec, new_pts)
            key = f"standoff{delta:+.1f}"
            variants[key] = {"_well_override": ws2}

        for vname, kw in variants.items():
            ws = kw.pop("_well_override", None)
            spec = ws if ws is not None else well_spec
            t0 = time.time()
            solver = AnnulusD2DGASolver(total_t=stop, nz=250, ny=40, **kw)
            res = solver.run(spec, fluids, provider, schedule=schedule)
            fin = res.metrics.iloc[-1]
            g = res.geom
            s_max = float(g["s"][-1])
            row = {
                "井": name,
                "配置": vname,
                "eta_E": round(float(fin["effective_efficiency"]), 4),
                "eta_N": round(float(res.summary["eta_narrow"]), 4),
                "窄边前缘m": round(float(fin["front_narrow_m"]), 1),
                "宽边前缘m": round(float(fin["front_wide_m"]), 1),
                "混浆指数": round(float(fin["mixing_index"]), 4),
                "壁面冻结占比": round(float(fin["mean_wall_mud"]), 4),
                "耗时s": round(time.time() - t0, 1),
            }
            rows.append(row)
            print(f"[exp2] {vname}: " + json.dumps(row, ensure_ascii=False), flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "消融_弥散与standoff.csv", index=False, encoding="utf-8-sig")
    return rows


def dataclasses_replace_standoff(well_spec, new_pts):
    """返回 standoff_profile 被替换的 WellSpec 副本（dataclasses.replace）。"""
    import dataclasses
    return dataclasses.replace(well_spec, standoff_profile=new_pts)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode in ("all", "exp1"):
        r1 = exp1_field_decomposition()
    if mode in ("all", "exp2"):
        r2 = exp2_ablation()
    print("[done] 结果已落盘", OUT)
