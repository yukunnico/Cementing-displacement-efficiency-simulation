"""2026-09-02 管容链修复后 8 井终跑（权威表第二轮）。

修复内容（commit 3d88385）：
- hu103/hu101/hu102 管容链统一由 shoe_lag_volume_m3 现场核实值驱动（88.55/102.6/88.9 m³）；
- F4 erf 收尾：弥散过渡带追加 frac=1.0 收尾事件，正体期入库浓度不再封顶 92.14%（全井生效）。

复用 rerun8_clean_20260902.py 的模式（每井跑完立即增量落盘，单井失败不中断），
输出目录改为 results/管容链修复后终跑_2026-09-02/，汇总含与上一轮
results/守恒修复后全量重跑_2026-09-02/汇总.csv 的逐井对照列（上轮η_E/Δη_E/
上轮stop/上轮入库水泥等）+ 前缘到位率.csv。不改 cemdisp 包内任何代码。

用法：
  PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/reruns/rerun8_after_pipe_capacity_fix_20260902.py [井名逗号列表]
"""
from __future__ import annotations

import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib  # noqa: E402

matplotlib.use("Agg")

from scripts.lib.mass_balance_diag import (  # noqa: E402
    WELLS,
    build_case,
    integrate_injection,
)
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver, _trapez2d  # noqa: E402
from cemdisp.reporting.animation import animate_cement_field  # noqa: E402

OUT = PROJECT_ROOT / "results" / "管容链修复后终跑_2026-09-02"
OUT.mkdir(parents=True, exist_ok=True)
PREV_DIR = PROJECT_ROOT / "results" / "守恒修复后全量重跑_2026-09-02"
NZ = 250
GIF_WELLS = ("hu101", "hu103")
# 本轮修复注记（逐井）
FIX_NOTE = {
    "hu103": "管容链 88.55m³(20315灌水+流量计双证;上轮 dual-diameter 仅尾管段≈31m³)",
    "hu101": "管容链 102.6m³(0708真实内径链;2011114 理论碰压102.2/实泵101.8;上轮52m³无出处)",
    "hu102": "管容链 88.9m³(149.2+114.3钻杆+尾管;上轮≈71m³漏计送入钻杆段)",
    "hu1": "仅 F4 erf 收尾(管容链未动)",
    "hu2": "仅 F4 erf 收尾(管容链未动)",
    "ht1_001": "仅 F4 erf 收尾(管容链未动)",
    "ht1_003": "仅 F4 erf 收尾(管容链未动)",
    "ht1_004": "仅 F4 erf 收尾(管容链未动)",
}
CEMENT_KEYS = ("lead", "tail", "cement")  # 与 integrate_injection / 2D inlet_tail 同一口径


def _load_prev_rows() -> dict:
    """读取上一轮（守恒修复后全量重跑）汇总，供逐井对照列。"""
    path = PREV_DIR / "汇总.csv"
    if not path.exists():
        return {}
    df = pd.read_csv(path)
    return {str(r["井"]): r for _, r in df.iterrows()}


def _load_stage1_stops() -> dict:
    """B3 前停止时刻（阶段1_库存核算.json），延续上轮对照口径。"""
    path = PROJECT_ROOT / "results" / "_质量平衡取证_2026-09-02" / "阶段1_库存核算.json"
    if path.exists():
        data = json.loads(path.read_text(encoding="utf-8"))
        return {k: float(v.get("停止时刻_s", np.nan)) for k, v in data.items()}
    return {}


class _RecordingProvider:
    """包装入口 provider，按求解器实际调用序列记录 (t, Q, 各相分数)。"""

    def __init__(self, provider):
        self._provider = provider
        self.calls: list[tuple[float, float, dict]] = []

    def __call__(self, t):
        st = self._provider(t)
        try:
            self.calls.append(
                (float(t), float(st.flow_rate_m3_s), dict(st.phase_fractions)))
        except Exception:
            pass
        return st


def replay_intake(calls, total_t: float) -> tuple[float, float]:
    """按主链（严格递增时间、左端点黎曼和）复算 2D 侧水泥/隔离液入库体积。"""
    cement = 0.0
    spacer = 0.0
    prev = None
    for t, q, frac in calls:
        if prev is not None and t > prev[0] + 1e-9:
            dt = min(t, total_t) - prev[0]
            if dt > 0.0:
                cement += prev[1] * sum(prev[2].get(k, 0.0) for k in CEMENT_KEYS) * dt
                spacer += prev[1] * prev[2].get("spacer", 0.0) * dt
        prev = (t, q, frac)
    if prev is not None and total_t > prev[0] + 1e-9:
        dt = total_t - prev[0]
        cement += prev[1] * sum(prev[2].get(k, 0.0) for k in CEMENT_KEYS) * dt
        spacer += prev[1] * prev[2].get("spacer", 0.0) * dt
    return cement, spacer


class _Tee:
    """同时写终端与每井日志文件。"""

    def __init__(self, *streams):
        self._streams = streams

    def write(self, s):
        for st in self._streams:
            st.write(s)
        return len(s)

    def flush(self):
        for st in self._streams:
            st.flush()


def run_one(name: str, nz: int = NZ, gen_gif: bool = True) -> dict:
    """跑单井并落盘全部单井产物，返回汇总行（含与上一轮对照列）。"""
    t0 = time.time()
    well_spec, fluids, schedule, provider, stop, v_cem_design = build_case(WELLS[name])
    v_ann = AnnulusD2DGASolver()._physical_annular_volume(well_spec)
    inj = integrate_injection(provider, stop)

    rec = _RecordingProvider(provider)
    solver = AnnulusD2DGASolver(total_t=stop, nz=nz, ny=40)  # runner 生产口径（B2 默认屈服门开）
    res = solver.run(well_spec, fluids, rec, schedule=schedule)
    elapsed = time.time() - t0

    g = res.geom
    v_dom = 2.0 * _trapez2d(g["b"] * res.cement_field, g)     # 域内最终水泥体积（2D 口径）
    v_domain_full = 2.0 * _trapez2d(g["b"], g)                # 域满体积（与物理环空体积按构造一致）
    fin = res.metrics.iloc[-1]
    sm = res.summary
    s_max = float(g["s"][-1])
    fw, fn = float(fin["front_wide_m"]), float(fin["front_narrow_m"])
    replay_cem, _ = replay_intake(rec.calls, stop)

    row = {
        "井": name,
        "本轮η_E": round(float(fin["effective_efficiency"]), 4),
        "η_N": round(float(sm.get("eta_narrow", np.nan)), 4),
        "stop_s": round(stop, 1),
        "入环空水泥_m3": round(inj["cement"], 2),
        "域内水泥_m3": round(v_dom, 2),
        "守恒率": round(v_dom / max(inj["cement"], 1e-9), 4),
        "库存比": round(v_cem_design / v_ann, 3),
        "宽边前缘m": round(fw, 1),
        "窄边前缘m": round(fn, 1),
        "域长m": round(s_max, 1),
        "宽边到位率": round(fw / s_max, 4),
        "窄边到位率": round(fn / s_max, 4),
        "壁面冻结占比": round(float(fin["mean_wall_mud"]), 4),
        "混浆指数": round(float(fin.get("mixing_index", np.nan)), 4),
        "耗时s": round(elapsed, 1),
        "shoe_lag_m3": well_spec.shoe_lag_volume_m3,
        # —— 与上一轮（守恒修复后全量重跑）对照列 ——
        "上轮η_E": _prev_val(name, "修正后η_E"),
        "Δη_E": _delta(float(fin["effective_efficiency"]), _prev_val(name, "修正后η_E")),
        "上轮stop_s": _prev_val(name, "stop_s"),
        "stop增量_s": _delta(stop, _prev_val(name, "stop_s")),
        "上轮入库水泥_m3": _prev_val(name, "入环空水泥_m3"),
        "Δ入库水泥_m3": _delta(inj["cement"], _prev_val(name, "入环空水泥_m3")),
        "上轮域内水泥_m3": _prev_val(name, "域内水泥_m3"),
        "Δ域内水泥_m3": _delta(v_dom, _prev_val(name, "域内水泥_m3")),
        "本轮修复注记": FIX_NOTE.get(name, ""),
        # —— 守恒记账口径取证列 ——
        "2D入库积分_m3": round(replay_cem, 2),
        "2D入库减1s积分_m3": round(replay_cem - inj["cement"], 2),
        "场体积漂移_m3": round(v_dom - replay_cem, 2),
        "域满体积_m3": round(v_domain_full, 2),
        "物理环空体积_m3": round(v_ann, 3),
    }

    # 单井产物：时间序列 / 深度剖面 / 摘要 / 2D场NPZ
    res.metrics.to_csv(OUT / f"{name}_时间序列.csv", index=False, encoding="utf-8-sig")
    res.depth_profiles.to_csv(OUT / f"{name}_深度剖面.csv", index=False, encoding="utf-8-sig")
    (OUT / f"{name}_摘要.json").write_text(
        json.dumps(res.summary, ensure_ascii=False, indent=2, default=float),
        encoding="utf-8")
    np.savez(
        OUT / f"{name}_2D场数据.npz",
        cement_snapshots=np.array(res.cement_snapshots),
        spacer_snapshots=np.array(res.spacer_snapshots),
        wall_snapshots=np.array(res.wall_snapshots),
        snapshot_times_s=np.array(res.snapshot_times_s),
        md=g["md"], y=g["y"],
        cement_final=res.cement_field, spacer_final=res.spacer_field,
        wall_final=res.wall_field,
    )

    # GIF（仅 GIF_WELLS；不查看内容，只记录文件/大小/帧数）
    if gen_gif and name in GIF_WELLS:
        try:
            animate_cement_field(res, output_dir=OUT, save_format="gif")
            candidates = sorted(OUT.glob("*_顶替过程动画.gif"),
                                key=lambda p: p.stat().st_mtime, reverse=True)
            gif_path = candidates[0] if candidates else OUT / "缺失.gif"
            manifest = {}
            mpath = OUT / "动画清单.json"
            if mpath.exists():
                manifest = json.loads(mpath.read_text(encoding="utf-8"))
            manifest[name] = {
                "文件": gif_path.name,
                "大小MB": round(gif_path.stat().st_size / 1e6, 2),
                "帧数": len(res.cement_snapshots),
                "快照时间点数": len(res.snapshot_times_s),
            }
            mpath.write_text(json.dumps(manifest, ensure_ascii=False, indent=2),
                             encoding="utf-8")
            print(f"[GIF] {name}: {manifest[name]}", flush=True)
        except Exception:
            print(f"[GIF] {name} 生成失败：", flush=True)
            traceback.print_exc()
    return row


def _prev_val(name: str, col: str):
    """上一轮汇总值（缺井/缺列返回 None）。"""
    row = _PREV_ROWS.get(name)
    if row is None or col not in row:
        return None
    v = row[col]
    return None if pd.isna(v) else round(float(v), 4)


def _delta(new, old):
    """差值（任一侧缺失返回 None）。"""
    if new is None or old is None:
        return None
    return round(float(new) - float(old), 4)


def _flush_tables(rows: list) -> None:
    """每井完成后立即增量重写汇总表（中断也保留已完成井；单井/全量调用同一文件）。"""
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "汇总.csv", index=False, encoding="utf-8-sig")
    (OUT / "汇总.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    ok = df[df.get("状态", "ok") == "ok"] if "状态" in df.columns else df
    if len(ok) and "宽边到位率" in ok.columns:
        front = ok[["井", "宽边前缘m", "窄边前缘m", "域长m",
                    "宽边到位率", "窄边到位率"]].copy()
        front["宽边到顶(≥域长-1m)"] = (ok["宽边前缘m"] >= ok["域长m"] - 1.0).map({True: "是", False: "否"})
        front["窄边到顶(≥域长-1m)"] = (ok["窄边前缘m"] >= ok["域长m"] - 1.0).map({True: "是", False: "否"})
        front.to_csv(OUT / "前缘到位率.csv", index=False, encoding="utf-8-sig")


def main(names: list[str], nz: int = NZ, gen_gif: bool = True) -> None:
    global _PREV_ROWS
    _PREV_ROWS = _load_prev_rows()
    stage1_stops = _load_stage1_stops()
    print(f"[rerun8_pipefix] 输出目录: {OUT}")
    print(f"[rerun8_pipefix] 上轮对照: {PREV_DIR}")
    print(f"[rerun8_pipefix] 井序: {names}  nz={nz}  GIF={GIF_WELLS if gen_gif else '无'}")
    rows: list = []
    for name in names:
        log_path = OUT / f"{name}_console.log"
        print(f"\n===== {name} =====", flush=True)
        with log_path.open("w", encoding="utf-8") as fh:
            old_stdout = sys.stdout
            sys.stdout = _Tee(old_stdout, fh)
            try:
                row = run_one(name, nz=nz, gen_gif=gen_gif)
                row["状态"] = "ok"
            except Exception as exc:
                row = {"井": name, "状态": "error",
                       "错误": f"{type(exc).__name__}: {exc}"}
                traceback.print_exc()
            finally:
                sys.stdout = old_stdout
        if "stop_s" in row and name in stage1_stops:
            row["stop_B3前_s"] = stage1_stops[name]
        rows.append(row)
        print(json.dumps(row, ensure_ascii=False, default=str), flush=True)
        _flush_tables(rows)  # 关键：每井后立即落盘
    print("\n[rerun8_pipefix] 完成。汇总：")
    print(pd.DataFrame(rows).to_string(index=False))


_PREV_ROWS: dict = {}

if __name__ == "__main__":
    names_arg = sys.argv[1].split(",") if len(sys.argv) > 1 else list(WELLS)
    main(names_arg)
