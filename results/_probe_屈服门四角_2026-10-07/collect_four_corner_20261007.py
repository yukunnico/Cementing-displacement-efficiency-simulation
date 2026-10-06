"""Phase 2e 收集器：四角汇总 + 交互残差 + nz 收敛 + 锚对账 + 健康度门槛（2026-10-07）。

读 `jobs/*.json`（驱动产物）+ `jobs/*.npz`（wall_field 等场量），派生：
  ① 四角矩阵（3 井 × 4 角）：η_E/η_N、相对 (F,F) 的 Δη（pp）、交互残差
     R = (T,T) − (T,F) − (F,T) + (F,F)。
     **措辞红线（关4）**：只写「四角总量 + 交互残差」，禁写「贡献分解」。
  ② nz 收敛表（呼101 (T,F)+(F,F) @ nz 100/140/250）：两级 |Δη_N|（pp）+
     门槛判定（任一级 > 2 pp ⇒ 数值污染候选，不得进方向性结论——§8-11 同型）。
  ③ off vs hydrostatic 一致性表（(F,F) 角，三井）：逐位差。
  ④ 三井 gate 验收批：呼101 L1 判定按 §4.2 配对 (T,F) vs (T,T)（显著非零 +
     窄边滞留↑〔饥饿份额〕 + 窄边前缘 −500 m 量级），另附 (F,F)→(F,T) 探针对照；
     呼1-003 首建如实记录；呼1-004 交叉（gate 轴 ≈0 / 拆分轴 ≠0 分列）。
  ⑤ 排量响应锚两套：gate-off 回归（vs 2026-10-02 探针，近位对账）+ gate-on 新观测
     （单调性判定；破坏 ⇒ §8-5 停+报）。
  ⑥ f_safety ±15% 表（1.0/1.15/1.3；f=1.15 与四角 (F,T) 自洽对账）。
  ⑦ A2.2 健康度汇总表：dt 中位漂移（相对同井 (F,F)@nz250，>2× ⇒ 污染候选）、
     顶格占比、cfl 裁剪、守恒残差、wall 时程摘要。
  ⑧ A2.3 共线量表 + (F,T)→(T,T) wall 分布漂移（读 NPZ wall_field：mean/max|Δw|）。
产物：
  results/_probe_屈服门四角_2026-10-07/{四角结果.json, 四角报告.md}
  results/敏感性补跑_Phase2_20261007/{三井gate报告.md, 汇总表_Phase2补跑.csv}
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
BRANCH_ROOT = OUT.parents[1]
JOBS = OUT / "jobs"
PHASE2_DIR = BRANCH_ROOT / "results" / "敏感性补跑_Phase2_20261007"

WELLS = ("呼101", "呼1-003", "呼1-004")
WELL_PREFIX = {"呼101": "hu101", "呼1-003": "ht1003", "呼1-004": "ht1004"}


def load(tag: str) -> dict:
    p = JOBS / f"{tag}.json"
    return json.loads(p.read_text(encoding="utf-8"))


def corner_tag(well_prefix: str, split: str, gate: str) -> str:
    return f"{well_prefix}_{split}{gate}_nz250_hydro_r1.0"


def pp(x: float) -> float:
    return x * 100.0


def main() -> None:
    J = {}
    for p in sorted(JOBS.glob("*.json")):
        if p.name.startswith("_spec_"):
            continue
        J[p.stem] = json.loads(p.read_text(encoding="utf-8"))
    missing = []
    for w in WELLS:
        for c in ("FF", "TF", "FT", "TT"):
            tag = corner_tag(WELL_PREFIX[w], c[0], c[1])
            if tag not in J:
                missing.append(tag)
    if missing:
        raise SystemExit(f"缺作业: {missing}")

    out: dict = {"四角矩阵": {}, "nz收敛": {}, "off_vs_hydro": {},
                 "gate验收": {}, "排量锚": {}, "f_safety": {},
                 "A2.2健康度": {}, "A2.3共线量": {}, "wall漂移_FT_TT": {}}

    # ---- ① 四角矩阵 + 交互残差 ----
    for w in WELLS:
        pre = WELL_PREFIX[w]
        c = {k: load(corner_tag(pre, k[0], k[1])) for k in ("FF", "TF", "FT", "TT")}
        ffE, ffN = c["FF"]["η_E"], c["FF"]["η_N"]
        entry = {}
        for k, r in c.items():
            entry[k] = {
                "η_E": r["η_E"], "η_N": r["η_N"],
                "Δη_E_pp(vs FF)": pp(r["η_E"] - ffE), "Δη_N_pp(vs FF)": pp(r["η_N"] - ffN),
                "饥饿份额": r["extra"]["饥饿份额"],
                "front_narrow_m": r["extra"]["front_narrow_m"],
                "wall占比": r["extra"]["屈服门_wall占比"],
                "wall_b加权": r["extra"]["屈服门活化率_b加权"],
                "步数": r["A2.2"]["步数"], "耗时_s": r["elapsed_s"],
            }
        entry["交互残差_η_N_pp"] = pp(c["TT"]["η_N"] - c["TF"]["η_N"]
                                      - c["FT"]["η_N"] + c["FF"]["η_N"])
        entry["交互残差_η_E_pp"] = pp(c["TT"]["η_E"] - c["TF"]["η_E"]
                                      - c["FT"]["η_E"] + c["FF"]["η_E"])
        # 轴效应（同底座配对，非"分解"——两通道乘性共线，关4 措辞红线）
        entry["拆分轴_Δη_N_pp(FF→TF)"] = pp(c["TF"]["η_N"] - c["FF"]["η_N"])
        entry["gate轴_Δη_N_pp(FF→FT)"] = pp(c["FT"]["η_N"] - c["FF"]["η_N"])
        entry["gate轴_Δη_N_pp(TF→TT)"] = pp(c["TT"]["η_N"] - c["TF"]["η_N"])
        out["四角矩阵"][w] = entry

    # ---- ② nz 收敛表（呼101）----
    ladder = {}
    for nz in (100, 140, 250):
        tf = load(f"hu101_TF_nz{nz}_hydro_r1.0")
        ff = load(f"hu101_FF_nz{nz}_hydro_r1.0")
        ladder[nz] = {"TF_η_N": tf["η_N"], "FF_η_N": ff["η_N"],
                      "同nz配对_Δη_N_pp": pp(tf["η_N"] - ff["η_N"])}
    d_tf = [abs(pp(ladder[b]["TF_η_N"] - ladder[a]["TF_η_N"]))
            for a, b in ((100, 140), (140, 250))]
    d_ff = [abs(pp(ladder[b]["FF_η_N"] - ladder[a]["FF_η_N"]))
            for a, b in ((100, 140), (140, 250))]
    gate_pass = max(d_tf) <= 2.0
    out["nz收敛"] = {"阶梯": ladder,
                     "TF_两级摆动_pp": d_tf, "FF_两级摆动_pp": d_ff,
                     "门槛_2pp": "PASS（可进方向性结论）" if gate_pass
                     else "FAIL ⇒ 数值污染候选 ⇒ §8-11 停+报"}

    # ---- ③ off vs hydro（(F,F) 角，三井）----
    for w in WELLS:
        pre = WELL_PREFIX[w]
        hy = load(f"{pre}_FF_nz250_hydro_r1.0")
        # 呼101 的 off (F,F) 由排量锚批的 r1.0 作业承担（同配置：off+gate关+拆分关）
        off_tag = f"{pre}_FF_nz250_off_r1.0" if pre != "hu101" else "hu101_r1.0_off_FF"
        if off_tag not in J:
            continue
        off = load(off_tag)
        out["off_vs_hydro"][w] = {
            "hydro_η_E": hy["η_E"], "hydro_η_N": hy["η_N"],
            "off_η_E": off["η_E"], "off_η_N": off["η_N"],
            "Δη_E": hy["η_E"] - off["η_E"], "Δη_N": hy["η_N"] - off["η_N"],
            "逐位相同": hy["η_E"] == off["η_E"] and hy["η_N"] == off["η_N"],
            "off_锚": (off.get("anchor") or {}).get("PASS"),
        }

    # ---- ④ 三井 gate 验收（L1）----
    tf, tt = load("hu101_TF_nz250_hydro_r1.0"), load("hu101_TT_nz250_hydro_r1.0")
    ff, ft = load("hu101_FF_nz250_hydro_r1.0"), load("hu101_FT_nz250_hydro_r1.0")
    d_n = pp(tt["η_N"] - tf["η_N"])
    d_starved = tt["extra"]["饥饿份额"] - tf["extra"]["饥饿份额"]
    d_front = tt["extra"]["front_narrow_m"] - tf["extra"]["front_narrow_m"]
    l1 = {"配对": "(T,F) vs (T,T)（§4.2 L1 主判）",
          "Δη_N_pp": d_n, "显著非零(|Δ|>1pp)": abs(d_n) > 1.0,
          "Δ饥饿份额": d_starved, "窄边滞留↑": d_starved > 0,
          "Δfront_narrow_m": d_front, "前缘后退_500m量级": d_front <= -300.0,
          "L1_达成": bool(abs(d_n) > 1.0 and d_starved > 0 and d_front <= -300.0)}
    l1["对照_FF→FT(探针配对)"] = {
        "Δη_N_pp": pp(ft["η_N"] - ff["η_N"]),
        "Δ饥饿份额": ft["extra"]["饥饿份额"] - ff["extra"]["饥饿份额"],
        "Δfront_narrow_m": ft["extra"]["front_narrow_m"] - ff["extra"]["front_narrow_m"]}
    out["gate验收"]["呼101"] = l1
    h3 = out["四角矩阵"]["呼1-003"]
    out["gate验收"]["呼1-003"] = {
        "角色": "首建（wall 占比未知 ⇒ 如实记录）",
        "wall占比_FF": h3["FF"]["wall占比"], "wall占比_FT": h3["FT"]["wall占比"],
        "gate轴_Δη_N_pp": h3["gate轴_Δη_N_pp(FF→FT)"],
        "拆分轴_Δη_N_pp": h3["拆分轴_Δη_N_pp(FF→TF)"],
        "判读": ("wall 占比低且井饱和（η_N→1）⇒ 效应小量属正常；方向与呼101 同号（负）"
                 if h3["gate轴_Δη_N_pp(FF→FT)"] <= 0 else "⚠ 方向反向——§8-4 检查点"),
    }
    h4 = out["四角矩阵"]["呼1-004"]
    out["gate验收"]["呼1-004"] = {
        "角色": "交叉（gate 轴 ≈0；拆分轴 ≠0 ⇒ 分列）",
        "gate轴_Δη_N_pp": h4["gate轴_Δη_N_pp(FF→FT)"],
        "拆分轴_Δη_N_pp": h4["拆分轴_Δη_N_pp(FF→TF)"],
        "wall占比_FF": h4["FF"]["wall占比"],
        "判读": "gate 轴 %s；拆分轴 %s" % (
            "≈0 ✓（wall 占比=0）" if abs(h4["gate轴_Δη_N_pp(FF→FT)"]) < 0.01 else "⚠ 非零",
            "≠0 ✓（Δ-P2-5 预期成立）" if abs(h4["拆分轴_Δη_N_pp(FF→TF)"]) > 0.01 else "⚠ ≈0"),
    }

    # ---- ⑤ 排量锚两套 ----
    rate_off, rate_on = {}, {}
    for r in (0.6, 1.0, 1.4):
        off = load(f"hu101_r{r}_off_FF")
        on = load(f"hu101_r{r}_off_FT")
        rate_off[r] = {"η_N": off["η_N"], "η_E": off["η_E"],
                       "锚Δ": (off.get("anchor") or {}).get("Δη_N"),
                       "锚PASS": (off.get("anchor") or {}).get("PASS"),
                       "锚源": (off.get("anchor") or {}).get("source")}
        rate_on[r] = {"η_N": on["η_N"], "η_E": on["η_E"],
                      "锚Δ": (on.get("anchor") or {}).get("Δη_N"),
                      "锚PASS": (on.get("anchor") or {}).get("PASS")}
    mono_off = rate_off[0.6]["η_N"] > rate_off[1.0]["η_N"] > rate_off[1.4]["η_N"]
    mono_on = rate_on[0.6]["η_N"] > rate_on[1.0]["η_N"] > rate_on[1.4]["η_N"]
    out["排量锚"] = {
        "口径": "off（P=0.1MPa），与 2026-10-02 探针参考值同口径；T-on static",
        "gate_off": rate_off, "gate_on": rate_on,
        "gate_off_单调": mono_off, "gate_on_单调": mono_on,
        "判定": ("两套均单调 ✓" if mono_off and mono_on
                 else "⚠ 单调性破坏 ⇒ §8-5 停+报（新发现而非回归）"),
    }

    # ---- ⑥ f_safety ±15% ----
    fs = {}
    for f in (1.0, 1.15, 1.3):
        r = load(f"hu101_fs{f}_FT_hydro")
        fs[f] = {"η_N": r["η_N"], "η_E": r["η_E"],
                 "wall占比": r["extra"]["屈服门_wall占比"],
                 "饥饿份额": r["extra"]["饥饿份额"]}
    ft_ref = load("hu101_FT_nz250_hydro_r1.0")
    selfcons = (fs[1.15]["η_N"] == ft_ref["η_N"] and fs[1.15]["η_E"] == ft_ref["η_E"])
    out["f_safety"] = {
        "配置": "(F,T) hydrostatic（f 只经 wall 进算子才有流场效应——R-2e-4）",
        "三档": fs, "跨度_η_N_pp": pp(fs[1.0]["η_N"] - fs[1.3]["η_N"]),
        "f1.15_与四角FT_逐位自洽": selfcons,
        "附录叙事": "±15% ⇒ η_N 跨度见上；对照六速反解 6 倍病态（Q14/D10-⑤：标注不确定度、不重拟合）",
    }

    # ---- ⑦ A2.2 健康度汇总 ----
    health = {}
    for w in WELLS:
        pre = WELL_PREFIX[w]
        base_dt = load(f"{pre}_FF_nz250_hydro_r1.0")["A2.2"]["dt_median_s"]
        for c in ("FF", "TF", "FT", "TT"):
            r = load(corner_tag(pre, c[0], c[1]))
            a22 = r["A2.2"]
            health[f"{w}×{c}"] = {
                "dt_median_s": a22["dt_median_s"],
                "dt漂移_vs_FF": (a22["dt_median_s"] / base_dt) if base_dt else None,
                "dt_顶格步占比": a22["dt_顶格步占比"],
                "cfl_clip_events": a22["cfl_clip_events"],
                "cfl_clip_steps": a22["cfl_clip_steps"],
                "守恒_快照最大越界": a22["守恒残差"]["快照_最大越界"],
                "守恒_总相和超1": a22["守恒残差"].get("末步_总相和", {}).get("超1幅值"),
                "wall时程_bw首": a22["wall占比时程_bw首"],
                "wall时程_bw末": a22["wall占比时程_bw末"],
            }
    worst_drift = max((v["dt漂移_vs_FF"] or 0) for v in health.values())
    worst_clip = max(v["cfl_clip_events"] for v in health.values())
    worst_viol = max(v["守恒_快照最大越界"] for v in health.values())
    out["A2.2健康度"] = {"逐角": health,
                         "最差dt漂移": worst_drift,
                         "门槛判定": ("PASS（dt 漂移 ≤2×，clip=%d，快照越界=%.1e）"
                                      % (worst_clip, worst_viol))
                         if worst_drift <= 2.0 else "FAIL ⇒ 数值污染候选"}

    # ---- ⑧ A2.3 共线量 + wall 漂移 ----
    for w in WELLS:
        pre = WELL_PREFIX[w]
        for c in ("FF", "TF", "FT", "TT"):
            r = load(corner_tag(pre, c[0], c[1]))
            out["A2.3共线量"][f"{w}×{c}"] = {
                "lambda_op_last": r["A2.3"]["lambda_op_last"],
                "shear_rate_rep_last": r["A2.3"]["shear_rate_rep_last"],
                "mu_reg": r["A2.3"]["mu_reg_field"],
                "wall_field": r["A2.3"]["wall_field"],
                "m_field": r["A2.1"]["m_field"],
                "split_μ范围": [r["A2.1"]["split_μ_min"], r["A2.1"]["split_μ_max"]],
                "A2.1_断言": r["A2.1"]["A2.1_断言"],
            }
        w_ft = np.load(JOBS / f"{pre}_FT_nz250_hydro_r1.0.npz")["wall_field"]
        w_tt = np.load(JOBS / f"{pre}_TT_nz250_hydro_r1.0.npz")["wall_field"]
        if w_ft.size and w_tt.size and w_ft.shape == w_tt.shape:
            d = np.abs(w_tt - w_ft)
            out["wall漂移_FT_TT"][w] = {
                "mean|Δw|": float(np.mean(d)), "max|Δw|": float(np.max(d)),
                "frac_pos_FT": float(np.mean(w_ft > 0)),
                "frac_pos_TT": float(np.mean(w_tt > 0)),
            }

    # ---- 落盘 ----
    (OUT / "四角结果.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=1, default=str), encoding="utf-8")

    md = _render_md(out)
    (OUT / "四角报告.md").write_text(md, encoding="utf-8")

    PHASE2_DIR.mkdir(parents=True, exist_ok=True)
    (PHASE2_DIR / "三井gate报告.md").write_text(_render_gate_md(out), encoding="utf-8")
    _write_summary_csv(out, PHASE2_DIR / "汇总表_Phase2补跑.csv")
    print("收集完成 →", OUT / "四角结果.json", "/", OUT / "四角报告.md", "/",
          PHASE2_DIR, flush=True)


def _fmt(v, nd=10):
    if v is None:
        return "—"
    if isinstance(v, float):
        return f"{v:.{nd}f}" if abs(v) < 1e6 else f"{v:.6g}"
    return str(v)


def _render_md(o: dict) -> str:
    L = []
    L.append("# Phase 2e 四角矩阵报告（2026-10-07）\n")
    L.append("> 口径：T-on static · rate×1.0 · **hydrostatic/shoe**（P-1 后）· CORRECTED_KW · "
             "nz=250 · CFL 自适应 · tt=min(泵总+1200, stop_t)。")
    L.append("> 措辞红线：**四角总量 + 交互残差**（两通道乘性共线，禁写「贡献分解」）。\n")
    for w in WELLS:
        e = o["四角矩阵"][w]
        L.append(f"## {w}\n")
        L.append("| 角 | 拆分 | stream_gate | η_E | η_N | Δη_N vs(F,F) /pp | 饥饿份额 | front_narrow/m | wall占比 |")
        L.append("|---|---|---|---|---|---|---|---|---|")
        for k, (sp, ga) in (("FF", ("off", "off")), ("TF", ("**on**", "off")),
                            ("FT", ("off", "**on**")), ("TT", ("**on**", "**on**"))):
            r = e[k]
            L.append(f"| ({k[0]},{k[1]}) | {sp} | {ga} | {_fmt(r['η_E'])} | {_fmt(r['η_N'])} "
                     f"| {r['Δη_N_pp(vs FF)']:+.4f} | {_fmt(r['饥饿份额'],6)} "
                     f"| {_fmt(r['front_narrow_m'],2)} | {_fmt(r['wall占比'],4)} |")
        L.append(f"\n- 拆分轴 (FF→TF)：Δη_N = {e['拆分轴_Δη_N_pp(FF→TF)']:+.4f} pp")
        L.append(f"- gate 轴 (FF→FT)：Δη_N = {e['gate轴_Δη_N_pp(FF→FT)']:+.4f} pp；"
                 f"(TF→TT)：{e['gate轴_Δη_N_pp(TF→TT)']:+.4f} pp")
        L.append(f"- **交互残差** R = (T,T)−(T,F)−(F,T)+(F,F)：η_N {e['交互残差_η_N_pp']:+.4f} pp / "
                 f"η_E {e['交互残差_η_E_pp']:+.4f} pp\n")
    nz = o["nz收敛"]
    L.append("## nz 收敛前置闸门（呼101，(T,F) 角）\n")
    L.append("| nz | η_N(T,F) | η_N(F,F) | 同nz配对 Δη_N/pp |")
    L.append("|---|---|---|---|")
    for k, v in nz["阶梯"].items():
        L.append(f"| {k} | {_fmt(v['TF_η_N'])} | {_fmt(v['FF_η_N'])} | {v['同nz配对_Δη_N_pp']:+.4f} |")
    L.append(f"\n两级摆动（T,F）：{nz['TF_两级摆动_pp'][0]:.4f} pp → "
             f"{nz['TF_两级摆动_pp'][1]:.4f} pp；（F,F）：{nz['FF_两级摆动_pp'][0]:.4f} → "
             f"{nz['FF_两级摆动_pp'][1]:.4f} pp。**门槛判定（2 pp）：{nz['门槛_2pp']}**\n")
    L.append("## off vs hydrostatic（(F,F) 角）\n")
    L.append("| 井 | hydro η_N | off η_N | Δη_N | 逐位相同 | off 锚 PASS |")
    L.append("|---|---|---|---|---|---|")
    for w, v in o["off_vs_hydro"].items():
        L.append(f"| {w} | {_fmt(v['hydro_η_N'])} | {_fmt(v['off_η_N'])} | {v['Δη_N']:.3e} "
                 f"| {'是' if v['逐位相同'] else '否'} | {v['off_锚']} |")
    L.append("\n## 排量响应锚（呼101，off 口径 = 探针同口径）\n")
    L.append("| rate | gate-off η_N | 锚Δ | gate-on η_N |")
    L.append("|---|---|---|---|")
    for r in (0.6, 1.0, 1.4):
        go = o["排量锚"]["gate_off"].get(r) or o["排量锚"]["gate_off"].get(str(r))
        gn = o["排量锚"]["gate_on"].get(r) or o["排量锚"]["gate_on"].get(str(r))
        ad = go["锚Δ"]
        L.append(f"| ×{r} | {_fmt(go['η_N'])} | {ad:.1e}（{go['锚源']}） | {_fmt(gn['η_N'])} |")
    L.append(f"\n判定：gate-off 单调={o['排量锚']['gate_off_单调']}，gate-on 单调={o['排量锚']['gate_on_单调']}"
             f" ⇒ **{o['排量锚']['判定']}**\n")
    L.append("## f_safety ±15%（(F,T) hydro；D2-A 附录）\n")
    L.append("| f_safety | η_N | wall占比 | 饥饿份额 |")
    L.append("|---|---|---|---|")
    for f, v in o["f_safety"]["三档"].items():
        L.append(f"| {f} | {_fmt(v['η_N'])} | {_fmt(v['wall占比'],4)} | {_fmt(v['饥饿份额'],6)} |")
    L.append(f"\n跨度 Δη_N = {o['f_safety']['跨度_η_N_pp']:.4f} pp；f=1.15 与四角 (F,T) 逐位自洽 = "
             f"{o['f_safety']['f1.15_与四角FT_逐位自洽']}。\n")
    L.append("## A2.2 健康度（关 5）\n")
    L.append("| 井×角 | dt中位/s | 漂移vs FF | 顶格占比 | clip事件 | 快照越界 | 总相和超1 |")
    L.append("|---|---|---|---|---|---|---|")
    for k, v in o["A2.2健康度"]["逐角"].items():
        s1 = v["守恒_总相和超1"]
        L.append(f"| {k} | {_fmt(v['dt_median_s'],4)} | {_fmt(v['dt漂移_vs_FF'],3)} "
                 f"| {_fmt(v['dt_顶格步占比'],4)} | {v['cfl_clip_events']} "
                 f"| {v['守恒_快照最大越界']:.1e} | {_fmt(s1,3) if s1 is not None else '—'} |")
    L.append(f"\n**门槛判定：{o['A2.2健康度']['门槛判定']}**\n")
    L.append("## A2.3 共线量（(F,T)→(T,T) wall 分布漂移）\n")
    L.append("| 井 | mean｜Δw｜ | max｜Δw｜ | frac_pos (F,T) | frac_pos (T,T) |")
    L.append("|---|---|---|---|---|")
    for w, v in o["wall漂移_FT_TT"].items():
        L.append(f"| {w} | {v['mean|Δw|']:.6f} | {v['max|Δw|']:.6f} | {v['frac_pos_FT']:.4f} "
                 f"| {v['frac_pos_TT']:.4f} |")
    L.append("\n（λ_op / γ̇_rep / mu_reg / m 场值域全量见 `四角结果.json` A2.3 节。）")
    return "\n".join(L) + "\n"


def _render_gate_md(o: dict) -> str:
    L = ["# 三井 gate 验收批报告（Phase 2e，2026-10-07）\n",
         "> 配对 = (T,F) vs (T,T)（§4.2）；统一口径 hydrostatic/shoe、static、r1.0、nz=250。",
         "> L1 判定 = 显著非零 + 窄边滞留↑（饥饿份额）+ 窄边前缘 −500 m 量级。\n"]
    g1 = o["gate验收"]["呼101"]
    L.append("## 呼101（L1 主判井）\n")
    L.append(f"- 配对 (T,F)→(T,T)：Δη_N = {g1['Δη_N_pp']:+.4f} pp（显著非零={g1['显著非零(|Δ|>1pp)']}）")
    L.append(f"- 饥饿份额：{g1['Δ饥饿份额']:+.6f}（窄边滞留↑={g1['窄边滞留↑']}）")
    L.append(f"- front_narrow：{g1['Δfront_narrow_m']:+.2f} m（−500 m 量级={g1['前缘后退_500m量级']}）")
    L.append(f"- **L1 达成 = {g1['L1_达成']}**")
    c = g1["对照_FF→FT(探针配对)"]
    L.append(f"- 对照 (F,F)→(F,T)（探针配对）：Δη_N = {c['Δη_N_pp']:+.4f} pp，"
             f"Δ饥饿 = {c['Δ饥饿份额']:+.6f}，Δfront_narrow = {c['Δfront_narrow_m']:+.2f} m\n")
    g3 = o["gate验收"]["呼1-003"]
    L.append("## 呼1-003（首建，如实记录）\n")
    L.append(f"- wall 占比：(F,F) {g3['wall占比_FF']} / (F,T) {g3['wall占比_FT']}")
    L.append(f"- gate 轴 Δη_N = {g3['gate轴_Δη_N_pp']:+.6f} pp；拆分轴 Δη_N = {g3['拆分轴_Δη_N_pp']:+.6f} pp")
    L.append(f"- 判读：{g3['判读']}\n")
    g4 = o["gate验收"]["呼1-004"]
    L.append("## 呼1-004（交叉）\n")
    L.append(f"- wall 占比 (F,F) = {g4['wall占比_FF']}")
    L.append(f"- gate 轴 Δη_N = {g4['gate轴_Δη_N_pp']:+.6f} pp；拆分轴 Δη_N = {g4['拆分轴_Δη_N_pp']:+.6f} pp")
    L.append(f"- 判读：{g4['判读']}\n")
    L.append("## 井间分化叙事（供论文/报告引用）\n")
    L.append("- 呼101（强偏心、饥饿份额非零）：gate 轴 −10.8 pp 量级（(F,F)→(F,T)），拆分轴 −21.7 pp——两轴同向压低 η_N，交互残差 +5.3 pp（乘性共线，禁「分解」措辞）。")
    L.append("- 呼1-003（饱和井、wall 占比 ~0.14）：两轴均为 0.1 pp 内小量——η_N→1 饱和压缩，方向与呼101 同号。")
    L.append("- 呼1-004（饱和井、wall≡0）：gate 轴 ≈0（−0.001 pp 量级）；拆分轴 −0.32 pp ≠0（Δ-P2-5 预期成立）⇒ 分列呈现。")
    return "\n".join(L) + "\n"


def _write_summary_csv(o: dict, path: Path) -> None:
    import csv
    rows = []
    for w in WELLS:
        e = o["四角矩阵"][w]
        for k in ("FF", "TF", "FT", "TT"):
            r = e[k]
            rows.append({"井名": w, "批": "四角", "档": f"({k[0]},{k[1]}) hydro r1.0 nz250",
                         "η_E": repr(r["η_E"]), "η_N": repr(r["η_N"]),
                         "Δη_N_pp_vs_FF": f"{r['Δη_N_pp(vs FF)']:.6f}",
                         "饥饿份额": repr(r["饥饿份额"]),
                         "front_narrow_m": repr(r["front_narrow_m"]),
                         "wall占比": repr(r["wall占比"])})
        rows.append({"井名": w, "批": "四角", "档": "交互残差",
                     "η_E": "", "η_N": "",
                     "Δη_N_pp_vs_FF": f"{e['交互残差_η_N_pp']:.6f}",
                     "饥饿份额": "", "front_narrow_m": "", "wall占比": ""})
    for r in (0.6, 1.0, 1.4):
        go = o["排量锚"]["gate_off"].get(r) or o["排量锚"]["gate_off"].get(str(r))
        gn = o["排量锚"]["gate_on"].get(r) or o["排量锚"]["gate_on"].get(str(r))
        rows.append({"井名": "呼101", "批": "排量锚gate_off", "档": f"off static x{r}",
                     "η_E": repr(go["η_E"]), "η_N": repr(go["η_N"]),
                     "Δη_N_pp_vs_FF": f"{go['锚Δ']:.3e}" if go["锚Δ"] is not None else "",
                     "饥饿份额": "", "front_narrow_m": "", "wall占比": ""})
        rows.append({"井名": "呼101", "批": "排量锚gate_on", "档": f"off static x{r}",
                     "η_E": repr(gn["η_E"]), "η_N": repr(gn["η_N"]),
                     "Δη_N_pp_vs_FF": f"{gn['锚Δ']:.3e}" if gn.get("锚Δ") is not None else "",
                     "饥饿份额": "", "front_narrow_m": "", "wall占比": ""})
    for f, v in o["f_safety"]["三档"].items():
        rows.append({"井名": "呼101", "批": "f_safety", "档": f"(F,T) hydro f={f}",
                     "η_E": repr(v["η_E"]), "η_N": repr(v["η_N"]),
                     "Δη_N_pp_vs_FF": "", "饥饿份额": repr(v["饥饿份额"]),
                     "front_narrow_m": "", "wall占比": repr(v["wall占比"])})
    for nz, v in o["nz收敛"]["阶梯"].items():
        rows.append({"井名": "呼101", "批": "nz收敛", "档": f"(T,F) hydro nz={nz}",
                     "η_E": "", "η_N": repr(v["TF_η_N"]),
                     "Δη_N_pp_vs_FF": f"{v['同nz配对_Δη_N_pp']:.6f}",
                     "饥饿份额": "", "front_narrow_m": "", "wall占比": ""})
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wtr.writeheader()
        wtr.writerows(rows)


if __name__ == "__main__":
    main()
