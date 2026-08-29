"""Task 13: 全量 8 井 baseline vs corrected 重跑（生产网格 nz=250）。

本脚本是论文/正式 8 井数字的官方入口（A2 裁定口径 (b)）：各 runner 保持默认基线口径，
正式数字一律以本脚本 corrected 配置为准。

corrected 配置：M1弥散dt归一 + M3屈服门槛 + I3局部化 + M4 e=0.90；M2流态修正统一开启
（层流元 R=1 不改变结果，对高Re/低黏井才生效）。输出 results/全井修正前后/汇总.csv，
并同步导出 adopted_config.json（solver 开关快照 + git_commit + 生成时间）。
nz=250+M2迭代，单井约3-10分钟。
"""
from __future__ import annotations
import csv, inspect, json, subprocess, time
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

def run_one(label, loader, *, corrected):
    t0=time.perf_counter()
    well,fluids,schedule,_=loader()
    cr=CasingFlowSolver(enable_gravity=True).run(well,fluids,schedule)
    inlet=build_coupled_annulus_inlet_provider(cr,CasingFlowSolver(enable_gravity=True),fluids,split_cement_phases=True)
    tt=min(_total_t(schedule)+1200.,_stop_t(cr,fluids)+600.)
    kw=dict(total_t=tt,nz=NZ)
    if corrected:
        kw.update(**CORRECTED_KW)
    res=AnnulusD2DGASolver(**kw).run(well,fluids,inlet)
    fr=cast(dict,res.summary["最终结果"])
    return {"well":label,"corrected":corrected,
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


def write_adopted_config():
    """导出 adopted 口径（corrected 配置）快照到 results/全井修正前后/adopted_config.json。

    内容：CORRECTED_KW solver 开关 + 边界口径 + 网格 + CFL 相关实际默认值 +
    git_commit + 生成时间，供论文数据溯源与复现。
    """
    params = inspect.signature(AnnulusD2DGASolver.__init__).parameters
    default = lambda name: params[name].default  # noqa: E731
    snapshot = {
        "说明": "论文/正式 8 井数字采用口径（corrected 配置）快照，由 scripts/rerun_all_wells_corrected.py 自动生成",
        "generated_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "git_commit": _git_commit(),
        "grid": {"nz": NZ, "ny": default("ny")},
        "solver": {
            **CORRECTED_KW,
            "enable_cfl_adaptive": default("enable_cfl_adaptive"),
            "cfl_number": default("cfl_number"),
        },
        "boundary": {
            "split_cement_phases": True,       # 环空入口桥接：领浆/中间浆并入前导水泥相，尾浆单独成相
            "casing_enable_gravity": True,     # 套管 1D 启用重力项
        },
    }
    path = OUT / "adopted_config.json"
    path.write_text(json.dumps(snapshot, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"配置快照: {path}", flush=True)


def main():
    OUT.mkdir(parents=True,exist_ok=True); rows=[]
    for label,loader in WELLS:
        print(f"\n=== {label} ===",flush=True)
        try:
            b=run_one(label,loader,corrected=False); rows.append(b)
            print(f"  base: eta_E={b['eta_E']:.4f} eta_N={b['eta_N']:.4f} mix={b['mixing']:.4f} ({b['elapsed_s']}s)",flush=True)
            c=run_one(label,loader,corrected=True); rows.append(c)
            print(f"  corr: eta_E={c['eta_E']:.4f} eta_N={c['eta_N']:.4f} mix={c['mixing']:.4f} ({c['elapsed_s']}s)",flush=True)
        except Exception as e:
            import traceback; traceback.print_exc()
            rows.append({"well":label,"corrected":"ERROR","error":str(e)})
    csvp=OUT/"汇总.csv"
    keys=sorted(set().union(*[r.keys() for r in rows]))
    with open(csvp,"w",newline="",encoding="utf-8-sig") as f:
        w=csv.DictWriter(f,fieldnames=keys); w.writeheader(); w.writerows(rows)
    write_adopted_config()
    print(f"\nCSV: {csvp}",flush=True)
    print("\n"+"="*78)
    print(f"  {'井':<10}{'eta_E base->corr':>18}{'eta_N base->corr':>18}{'mix base->corr':>16}")
    for label,_ in WELLS:
        rb=next((r for r in rows if r['well']==label and r.get('corrected') is False),None)
        rc=next((r for r in rows if r['well']==label and r.get('corrected') is True),None)
        if rb and rc:
            print(f"  {label:<10}{rb['eta_E']:.3f}->{rc['eta_E']:.3f}        {rb['eta_N']:.3f}->{rc['eta_N']:.3f}        {rb['mixing']:.3f}->{rc['mixing']:.3f}")

if __name__=="__main__":
    main()
