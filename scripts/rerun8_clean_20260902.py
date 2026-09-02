"""2026-09-02 B1/B2/B3 修复后 8 井 nz=250 干净全量重跑（权威表）。

修复 _rerun8_after_conservation_fix_20260902.py（阶段4）的两个落盘缺陷：
  1) CSV/JSON 只在全部井跑完后一次性写出——中途中断即全丢；
  2) 单井调用（如 `... 250 hu103`）会整体覆盖 CSV——阶段4 nz250 CSV 只剩
     hu103 一行（B3 后口径）而控制台日志是 B3 前口径的直接原因。
本脚本：每井跑完立即增量重写 汇总.csv / 汇总.json（保留已完成井），单井失败
不中断循环；另落盘每井 console 日志 / 时间序列 / 深度剖面 / 摘要 / 2D场NPZ，
hu101+hu103 生成顶替动画 GIF，并按 provider 调用序列复算 2D 实际入库积分
（守恒记账口径取证，不改 cemdisp 包内任何代码）。

用法：
  PYTHONIOENCODING=utf-8 PYTHONUTF8=1 python scripts/rerun8_clean_20260902.py [井名逗号列表]
"""
from __future__ import annotations

import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib  # noqa: E402

matplotlib.use("Agg")

from scripts._mass_balance_diag_20260902 import (  # noqa: E402
    WELLS,
    build_case,
    integrate_injection,
)
from cemdisp.models2d.annulus_d2dga import AnnulusD2DGASolver, _trapez2d  # noqa: E402
from cemdisp.reporting.animation import animate_cement_field  # noqa: E402

OUT = PROJECT_ROOT / "results" / "守恒修复后全量重跑_2026-09-02"
OUT.mkdir(parents=True, exist_ok=True)
NZ = 250
GIF_WELLS = ("hu101", "hu103")
# 09-02 修正前 runner 基线（nz=250，results/全井runner重跑_停止修复_2026-09-02/汇总.md）
OLD_ETA = {
    "hu101": 0.3992, "hu102": 0.3814, "hu103": 0.1940, "hu1": 0.5014,
    "hu2": 0.5690, "ht1_001": 0.5195, "ht1_003": 0.6132, "ht1_004": 0.6284,
}
# B3 前停止时刻（阶段1_库存核算.json，B3 提交前口径，用于对照）
STAGE1_JSON = PROJECT_ROOT / "results" / "_质量平衡取证_2026-09-02" / "阶段1_库存核算.json"
CEMENT_KEYS = ("lead", "tail", "cement")  # 与 integrate_injection / 2D inlet_tail 同一口径


def _load_stage1_stops() -> dict:
    if STAGE1_JSON.exists():
        data = json.loads(STAGE1_JSON.read_text(encoding="utf-8"))
        return {k: float(v.get("停止时刻_s", np.nan)) for k, v in data.items()}
    return {}


class _RecordingProvider:
    """包装入口 provider，按求解器实际调用序列记录 (t, Q, 各相分数)。

    用于按 2D 自己的时间步序列复算"实际入库体积"，与 integrate_injection 的
    1s 采样积分互为对照（守恒记账口径取证）。不改求解器行为。
    """

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
    # 主链最后一步 [t_K, total_t]（solver 末步裁剪到 total_t，但不再调用 provider）
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
    """跑单井并落盘全部单井产物，返回汇总行（不含状态列）。"""
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
        "修正前η_E": OLD_ETA.get(name),
        "修正后η_E": round(float(fin["effective_efficiency"]), 4),
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
            # 文件名以 plots._safe_filename_component(井名) 为准，按模式匹配实际产物
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
    stage1_stops = _load_stage1_stops()
    print(f"[rerun8_clean] 输出目录: {OUT}")
    print(f"[rerun8_clean] 井序: {names}  nz={nz}  GIF={GIF_WELLS if gen_gif else '无'}")
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
            row["stop增量_s"] = round(row["stop_s"] - stage1_stops[name], 1)
        rows.append(row)
        print(json.dumps(row, ensure_ascii=False), flush=True)
        _flush_tables(rows)  # 关键：每井后立即落盘
    print("\n[rerun8_clean] 完成。汇总：")
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    names_arg = sys.argv[1].split(",") if len(sys.argv) > 1 else list(WELLS)
    main(names_arg)
