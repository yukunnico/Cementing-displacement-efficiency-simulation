"""Task 13: 全量 8 井 baseline vs corrected 重跑（生产网格 nz=250）。

本脚本是论文/正式 8 井数字的官方入口（A2 裁定口径 (b)）：各 runner 保持默认基线口径，
正式数字一律以本脚本 corrected 配置为准。

corrected 配置：M1弥散dt归一 + M3屈服门槛 + I3局部化 + M4 e=0.90；M2流态修正统一开启
（层流元 R=1 不改变结果，对高Re/低黏井才生效）。输出 results/全井修正前后/汇总.csv，
并同步导出 adopted_config.json（solver 开关快照 + git_commit + data_state 数据状态 + 生成时间）。
nz=250+M2迭代，单井约3-10分钟。

Task A1/A3（2026-08-29）扩展：新增命令行参数，默认行为不变——
  --out-dir PATH   输出目录（默认 results/全井修正前后）
  --cfl-mode on|fixed_dt   CFL 口径：on=自适应（默认）；fixed_dt=固定 dt=4s（enable_cfl_adaptive=False）
  --skip-baseline  只跑 corrected 组（A1/A3 口径裁定与最终基线重跑用）
"""
from __future__ import annotations
import argparse, csv, inspect, json, subprocess, time
from pathlib import Path
from typing import cast
import numpy as np

from cemdisp.data.fluid_spec import FluidRole
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver
import cemdisp.data.loaders as L

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT = PROJECT_ROOT / "results" / "全井修正前后"
NZ = 250

# corrected（adopted）口径唯一来源：run_one 与 adopted_config.json 快照共用，防止漂移。
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

def _total_t(schedule):
    return sum(0. if s.rate_m3_min<=0 else s.volume_m3/s.rate_m3_min*60. for s in schedule.steps)

def _stop_t(cr, fluids):
    roles={FluidRole.LEAD,FluidRole.INTERMEDIATE,FluidRole.TAIL}; by={f.name:f for f in fluids}
    fs=sorted(cr.fronts,key=lambda f:f.time_s)
    last=next((f.time_s for f in fs if by.get(f.fluid_name) and by[f.fluid_name].role in roles),None)
    if last is not None:
        for f in fs:
            fl=by.get(f.fluid_name)
            if fl is None or fl.role in roles: continue
            if f.time_s>=last-1e-9: return float(f.time_s)
    return float(cr.cement_end_time_s)

def run_one(label, loader, *, corrected, cfl_on=True):
    t0=time.perf_counter()
    well,fluids,schedule,_=loader()
    cr=CasingFlowSolver(enable_gravity=True).run(well,fluids,schedule)
    inlet=build_coupled_annulus_inlet_provider(cr,CasingFlowSolver(enable_gravity=True),fluids,split_cement_phases=True)
    tt=min(_total_t(schedule)+1200.,_stop_t(cr,fluids)+600.)
    kw=dict(total_t=tt,nz=NZ,enable_cfl_adaptive=cfl_on)
    if corrected:
        kw.update(**CORRECTED_KW)
    res=AnnulusD2DGASolver(**kw).run(well,fluids,inlet)
    fr=cast(dict,res.summary["最终结果"])
    return {"well":label,"corrected":corrected,"cfl_mode":"on" if cfl_on else "fixed_dt",
        "eta_E":float(fr["全井段最终有效顶替效率"]),
        "eta_N":float(fr.get("窄四分位效率",float("nan"))),
        "mixing":float(fr["最终混浆指数"]),"channeling":float(fr["最终窜槽指数"]),
        "instability":float(fr["最终失稳指数"]),"cement_occ":float(fr["最终水泥浆占据率"]),
        "wall_frac":float(np.mean(res.wall_field)) if res.wall_field is not None else float("nan"),
        "elapsed_s":round(time.perf_counter()-t0,1)}

def _git_commit():
    """返回当前 HEAD 提交号；git 不可用时返回 'unknown'。"""
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT,
            capture_output=True, text=True, check=True,
        ).stdout.strip()
    except Exception:
        return "unknown"


def _data_state():
    """数据状态摘要（Task A1/A3 数字可追溯要求）：git HEAD + 未提交的 loader/现场资料改动行 + 内容 hash。

    本次重跑消费工作区中另一会话未提交的 loader（cemdisp/data/loaders/*）与
    参考文档/现场资料提取/* 数据——这些改动不在任何 commit 里，必须显式记录
    （git status 行 + git hash-object 内容 hash），否则结果数字无法用 git_commit 单独复现。
    core.quotepath=false 让中文路径保持原文而非 \\3xx 八进制转义。
    """
    try:
        st = subprocess.run(
            ["git", "-c", "core.quotepath=false", "status", "--short"], cwd=PROJECT_ROOT,
            capture_output=True, text=True, check=True,
        ).stdout.splitlines()
    except Exception:
        st = []
    keys = ("cemdisp/data/", "参考文档/现场资料提取")
    relevant = [l.strip() for l in st if any(k in l for k in keys)]
    # 对被修改（M，不含删除/未跟踪）的数据文件算内容 hash（git hash-object，与 git 生态可互核）
    hashes = {}
    for line in relevant:
        parts = line.split(maxsplit=1)
        if len(parts) == 2 and parts[0] in ("M", "MM") and parts[1].strip('"'):
            p = PROJECT_ROOT / parts[1].strip().strip('"')
            if p.is_file():
                try:
                    h = subprocess.run(
                        ["git", "hash-object", str(p)], cwd=PROJECT_ROOT,
                        capture_output=True, text=True, check=True,
                    ).stdout.strip()
                    hashes[parts[1].strip('"')] = h
                except Exception:
                    hashes[parts[1].strip('"')] = "hash_failed"
    return {
        "git_commit": _git_commit(),
        "说明": "重跑消费的工作区数据状态：以下 loader/现场资料改动存在于工作区但未必已提交，"
               "复现数字需同时核对 commit 与这些未提交改动（hash 为 git hash-object 内容指纹）",
        "uncommitted_or_modified": relevant,
        "modified_file_hashes": hashes,
    }


def write_adopted_config(out_dir: Path, *, cfl_mode: str):
    """导出 adopted 口径（corrected 配置）快照到 <out_dir>/adopted_config.json。

    内容：CORRECTED_KW solver 开关 + 边界口径 + 网格 + CFL 口径（实际运行模式）+
    git_commit + data_state 数据状态 + 生成时间，供论文数据溯源与复现。
    """
    params = inspect.signature(AnnulusD2DGASolver.__init__).parameters
    default = lambda name: params[name].default  # noqa: E731
    snapshot = {
        "说明": "论文/正式 8 井数字采用口径（corrected 配置）快照，由 scripts/rerun_all_wells_corrected.py 自动生成",
        "generated_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "git_commit": _git_commit(),
        "data_state": _data_state(),
        "grid": {"nz": NZ, "ny": default("ny")},
        "cfl": {
            "mode": cfl_mode,  # 本次运行实际口径：on=CFL 自适应 / fixed_dt=固定 dt
            "enable_cfl_adaptive": cfl_mode == "on",
            "cfl_number": default("cfl_number"),
            "dt_fixed_s": default("dt") if cfl_mode != "on" else None,
        },
        "solver": {
            **CORRECTED_KW,
            "enable_cfl_adaptive_default": default("enable_cfl_adaptive"),
        },
        "boundary": {
            "split_cement_phases": True,       # 环空入口桥接：领浆/中间浆并入前导水泥相，尾浆单独成相
            "casing_enable_gravity": True,     # 套管 1D 启用重力项
        },
    }
    path = out_dir / "adopted_config.json"
    path.write_text(json.dumps(snapshot, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"配置快照: {path}", flush=True)


def _save_well_json(results_dir: Path, row: dict):
    """单井结果立即落盘（断点续跑依据）：每井跑完即写 JSON，不等全量汇总。

    背景：首次后台重跑在 hu102 期间被会话中断杀掉，hu101 结果只存在于日志；
    改为逐井落盘后，中断恢复时已完成井自动跳过。
    """
    results_dir.mkdir(parents=True, exist_ok=True)
    tag = "base" if row.get("corrected") is False else ("corrected" if row.get("corrected") is True else "error")
    path = results_dir / f"{row['well']}_{tag}_{row.get('cfl_mode','?')}.json"
    path.write_text(json.dumps(row, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"  单井落盘: {path.name}", flush=True)


def _load_saved_well(results_dir: Path, label: str, tag: str, cfl_mode: str) -> dict | None:
    """读取已落盘单井结果；不存在返回 None。"""
    path = results_dir / f"{label}_{tag}_{cfl_mode}.json"
    if path.exists():
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            return None
    return None


def main():
    ap = argparse.ArgumentParser(description="8 井 corrected 口径重跑（官方论文数字入口）")
    ap.add_argument("--out-dir", type=Path, default=OUT, help="输出目录（默认 results/全井修正前后）")
    ap.add_argument("--cfl-mode", choices=["on", "fixed_dt"], default="on",
                    help="CFL 口径：on=自适应（默认）；fixed_dt=固定 dt=4s（Task A1 附表）")
    ap.add_argument("--skip-baseline", action="store_true",
                    help="只跑 corrected 组（Task A1/A3：口径对比与最终基线不需要 base 组）")
    args = ap.parse_args()
    out_dir = args.out_dir if args.out_dir.is_absolute() else PROJECT_ROOT / args.out_dir
    cfl_on = args.cfl_mode == "on"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"输出目录: {out_dir} | CFL 模式: {args.cfl_mode} | baseline 组: {'跳过' if args.skip_baseline else '运行'}", flush=True)
    rows=[]
    results_dir = out_dir / "单井结果"
    results_dir.mkdir(parents=True, exist_ok=True)
    for label,loader in WELLS:
        print(f"\n=== {label} ===",flush=True)
        try:
            b = c = None
            if not args.skip_baseline:
                # 断点续跑：已落盘的单井结果直接复用（日志曾因会话中断丢失 hu101 结果而设）
                b = _load_saved_well(results_dir, label, "base", args.cfl_mode)
                if b is None:
                    b = run_one(label,loader,corrected=False,cfl_on=cfl_on)
                    _save_well_json(results_dir, b)
                rows.append(b)
                print(f"  base: eta_E={b['eta_E']:.4f} eta_N={b['eta_N']:.4f} mix={b['mixing']:.4f} ({b['elapsed_s']}s)",flush=True)
            c = _load_saved_well(results_dir, label, "corrected", args.cfl_mode)
            if c is None:
                c = run_one(label,loader,corrected=True,cfl_on=cfl_on)
                _save_well_json(results_dir, c)
            rows.append(c)
            print(f"  corr: eta_E={c['eta_E']:.4f} eta_N={c['eta_N']:.4f} mix={c['mixing']:.4f} ({c['elapsed_s']}s)",flush=True)
        except Exception as e:
            import traceback; traceback.print_exc()
            rows.append({"well":label,"corrected":"ERROR","cfl_mode":args.cfl_mode,"error":str(e)})
    csvp=out_dir/"汇总.csv"
    keys=sorted(set().union(*[r.keys() for r in rows]))
    with open(csvp,"w",newline="",encoding="utf-8-sig") as f:
        w=csv.DictWriter(f,fieldnames=keys); w.writeheader(); w.writerows(rows)
    write_adopted_config(out_dir, cfl_mode=args.cfl_mode)
    print(f"\nCSV: {csvp}",flush=True)
    print("\n"+"="*78)
    has_base = any(r.get('corrected') is False for r in rows)
    print(f"  {'井':<10}{'eta_E base->corr':>18}{'eta_N base->corr':>18}{'mix base->corr':>16}" if has_base
          else f"  {'井':<10}{'eta_E corrected':>16}{'eta_N':>10}{'mixing':>10}{'wall%':>8}")
    for label,_ in WELLS:
        rb=next((r for r in rows if r['well']==label and r.get('corrected') is False),None)
        rc=next((r for r in rows if r['well']==label and r.get('corrected') is True),None)
        if rb and rc:
            print(f"  {label:<10}{rb['eta_E']:.3f}->{rc['eta_E']:.3f}        {rb['eta_N']:.3f}->{rc['eta_N']:.3f}        {rb['mixing']:.3f}->{rc['mixing']:.3f}")
        elif rc:
            print(f"  {label:<10}{rc['eta_E']:.4f}        {rc['eta_N']:.4f}    {rc['mixing']:.4f}  {rc['wall_frac']:.4f}")

if __name__=="__main__":
    main()
