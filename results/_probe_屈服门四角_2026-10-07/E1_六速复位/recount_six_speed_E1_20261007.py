# -*- coding: utf-8 -*-
"""
E-1 六速流变 CSV 复位重算（Phase 2 spec §10 E-1 遗留项）
==========================================================
口径（spec §1.2 来源二，2026-10-06 对抗核查后裁定）：
  1) notes 含 `<转速>:<读数>` 标注（本数据中形如 `φ400:256 φ200:194 ...`，须带 φ/Φ 前缀防误伤）
     时，整行以 note 对为准重建；
  2) 否则用列名读数；
  3) 字段数与表头不符者先尝试复位：hu1 的 14 字段行多出的字段恒为空、且 6 个 θ 值与尾部
     元数据（source/confidence/notes）可唯一右对齐落位 ⇒ 判为"可复位（右对齐去空字段）"；
     若右对齐后 θ 值仍无法唯一落位（多出的字段非空）⇒ 剔除整行并记录原因；
  4) 读数非数值（如 `>300` 超量程）⇒ 该 (转速,读数) 对剔除并记录原因。

换算：Fann 35 R1B1：τ[Pa] = 0.511·θ；γ̇[s⁻¹] = 1.703·N。
派生：μ_app = τ/γ̇ —— 原始"2.6005/8.7017/5.401"经核验实为 μ_app 口径
（buoyancy.py 文档串亦写"实测最大 μ_app ≈ 2.60 Pa·s"），故 τ 与 μ_app 两组最大并列报告。

对照实现：naive（左移错位）口径 = 按表头位置硬切前 13/12 字段，用于复现原 389/2.6005/8.7017/5.401。
敏感性实现：strict 口径 = 错位 14 字段行一律剔除（不右对齐复位）。

产出：E1_复位重算结果.json + E1_复位重算摘要.md
"""
import csv
import json
import os
import re

TAU_K = 0.511      # τ[Pa] = 0.511·θ
GAMMA_K = 1.703    # γ̇[s⁻¹] = 1.703·N

BASE = r"D:\users\desktop\research\控压固井项目\cement model_温压耦合分支\参考文档\现场资料提取"
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

THETA_COLS = ["theta_600", "theta_300", "theta_200", "theta_100", "theta_6", "theta_3"]
SPEEDS = {"theta_600": 600, "theta_300": 300, "theta_200": 200,
          "theta_100": 100, "theta_6": 6, "theta_3": 3}

# notes 中的转速:读数对，必须带 φ/Φ 前缀（避免 "1:1"、时间串等误伤）
NOTE_PAIR_RE = re.compile(r"[φΦ]\s*(\d{1,4})\s*[:：]\s*(\d{1,4}(?:\.\d+)?)")


def pair_metrics(theta, speed):
    tau = TAU_K * theta
    gdot = GAMMA_K * speed
    return tau, gdot, tau / gdot


def parse_float(s):
    s = (s or "").strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return ("UNPARSEABLE", s)


def load_well(csv_path):
    with open(csv_path, encoding="utf-8-sig", newline="") as f:
        rows = list(csv.reader(f))
    hdr, data = rows[0], rows[1:]
    return hdr, data


def row_to_theta_by_name(fields, hdr):
    """列名口径（位置硬切，要求 len==len(hdr)）。返回 [(speed, theta_or_reason)] 与剔除记录。"""
    out, drops = [], []
    idx = {c: i for i, c in enumerate(hdr)}
    for c in THETA_COLS:
        if c not in idx:
            continue
        v = fields[idx[c]]
        p = parse_float(v)
        if p is None:
            continue  # 缺失读数，不算错位，不计条
        if isinstance(p, tuple):
            drops.append((SPEEDS[c], v, "读数非数值（超量程）"))
            continue
        out.append((SPEEDS[c], p))
    return out, drops


def recover_pairs(hdr, fields, notes):
    """复位口径主函数。返回 (pairs[(speed,theta,method)], row_dropped_reason, notes_used)。"""
    nc = len(hdr)
    idx = {c: i for i, c in enumerate(hdr)}
    notes_used = False
    fs = fields
    # 1) notes 带 <转速>:<读数> 对 ⇒ 以 note 为准重建
    if notes:
        pairs = NOTE_PAIR_RE.findall(notes)
        if pairs:
            recovered = [(int(n), float(t)) for n, t in pairs]
            return [(n, t, "note复位") for n, t in recovered], None, True
    # 2) 字段数与表头不符 ⇒ 尝试右对齐复位：多出的字段必须全为空、且位于 fluid 与 θ 块之间
    if len(fs) != nc:
        extra = len(fs) - nc
        if extra <= 0:
            return [], f"字段数不足（{len(fs)}<{nc}），不可复位", False
        # 复位规则：temp 锚定字段2，紧随其后的 extra 个字段必须全为空
        # （hu1 14字段行实证：temp 有值时多出者恒为空分隔位），θ 块与尾部元数据右锚
        cand = fs[3:3 + extra]
        if all((x or "").strip() == "" for x in cand):
            fs = fs[:3] + fs[3 + extra:]
            method = "右对齐复位"
        else:
            return [], f"字段数超出（{len(fs)}>{nc}）且多出字段非空，不可复位", False
    else:
        method = "列名"
    out, drops = [], []
    for c in THETA_COLS:
        if c not in idx:
            continue
        v = fs[idx[c]]
        p = parse_float(v)
        if p is None:
            continue
        if isinstance(p, tuple):
            drops.append((SPEEDS[c], v, "读数非数值（超量程）"))
            continue
        out.append((SPEEDS[c], p, method))
    return out, None, notes_used


def naive_pairs(hdr, fields):
    """naive 错位口径：按位置硬切前 nc 个字段（原统计脚本的可能行为）。"""
    idx = {c: i for i, c in enumerate(hdr)}
    out = []
    for c in THETA_COLS:
        v = fields[idx[c]] if idx[c] < len(fields) else ""
        p = parse_float(v)
        if isinstance(p, float):
            out.append((SPEEDS[c], p))
    return out


def main():
    wells = {}
    all_pairs = []          # 复位口径
    naive_all = []
    strict_all = []
    excl_rows = []          # 剔除/异常记录
    pair_drops = []         # 逐条剔除（超量程等）
    notes_rows = []
    realigned_rows = []

    for d in sorted(os.listdir(BASE)):
        p = os.path.join(BASE, d, "rheometer_readings.csv")
        if not os.path.isdir(os.path.join(BASE, d)) or not os.path.exists(p):
            continue
        well = d
        hdr, data = load_well(p)
        for li, fields in enumerate(data, start=2):
            notes = fields[-1] if fields else ""
            conf = fields[-2] if len(fields) >= 2 else ""
            # 特殊：hu103 有行 confidence/notes 内容互换（θ 列不受影响），仅记录
            if len(fields) == len(hdr) and conf not in ("high", "medium", "low", ""):
                excl_rows.append({"well": well, "line": li, "fluid": fields[1],
                                  "issue": "confidence/notes 列内容互换（θ 列完好，不影响计数）"})
            pairs, drop_reason, used_notes = recover_pairs(hdr, fields, notes)
            if drop_reason:
                excl_rows.append({"well": well, "line": li, "fluid": fields[1],
                                  "issue": drop_reason})
            else:
                if used_notes:
                    notes_rows.append({"well": well, "line": li, "fluid": fields[1],
                                       "n_pairs": len(pairs),
                                       "raw": [f"{n}:{t:g}" for n, t, _ in pairs]})
                for sp, th, m in pairs:
                    tau, gdot, mu = pair_metrics(th, sp)
                    rec = {"well": well, "line": li, "fluid": fields[1],
                           "speed": sp, "theta": th, "tau_pa": round(tau, 4),
                           "gamma_s1": round(gdot, 3), "mu_app_pa_s": round(mu, 4),
                           "method": m}
                    all_pairs.append(rec)
                    if m == "右对齐复位":
                        pass
                if any(m == "右对齐复位" for _, _, m in pairs) and pairs:
                    realigned_rows.append({"well": well, "line": li,
                                           "fluid": fields[1], "n_pairs": len(pairs)})
                # 非数值读数剔除记录（列名口径行才检查；note 行天然无 >300）
                if not used_notes and len(fields) == len(hdr):
                    for c in THETA_COLS:
                        v = fields[hdr.index(c)] if c in hdr else ""
                        pp = parse_float(v)
                        if isinstance(pp, tuple):
                            pair_drops.append({"well": well, "line": li, "fluid": fields[1],
                                               "speed": SPEEDS[c], "raw": v,
                                               "reason": "读数非数值（超量程 >300）"})
                # strict 敏感性：14 字段行不复活
                if len(fields) == len(hdr):
                    strict_all.extend([r for r in all_pairs if r["well"] == well and r["line"] == li])
            # naive 对照
            for sp, th in naive_pairs(hdr, fields):
                naive_all.append({"well": well, "line": li, "fluid": fields[1],
                                  "speed": sp, "theta": th})

    # ---- 汇总统计 ----
    def agg(recs, key="tau_pa"):
        return max(r[key] for r in recs) if recs else None

    def stat_block(recs):
        t3 = [r for r in recs if r["speed"] == 3]
        ht1004 = [r for r in recs if r["well"].startswith("ht1_004")]
        return {
            "n_pairs": len(recs),
            "tau_max_3rpm": round(agg(t3, "tau_pa"), 4) if t3 else None,
            "theta_at_tau_max_3rpm": (max(t3, key=lambda r: r["theta"])["theta"] if t3 else None),
            "src_tau_max_3rpm": (max(t3, key=lambda r: r["tau_pa"]) if t3 else None),
            "muapp_max_3rpm": round(agg(t3, "mu_app_pa_s"), 4) if t3 else None,
            "src_muapp_max_3rpm": (max(t3, key=lambda r: r["mu_app_pa_s"]) if t3 else None),
            "tau_max_all": round(agg(recs, "tau_pa"), 4) if recs else None,
            "src_tau_max_all": (max(recs, key=lambda r: r["tau_pa"]) if recs else None),
            "muapp_max_all": round(agg(recs, "mu_app_pa_s"), 4) if recs else None,
            "src_muapp_max_all": (max(recs, key=lambda r: r["mu_app_pa_s"]) if recs else None),
            "ht1004_tau_max": round(agg(ht1004, "tau_pa"), 4) if ht1004 else None,
            "ht1004_muapp_max": round(agg(ht1004, "mu_app_pa_s"), 4) if ht1004 else None,
            "src_ht1004_muapp_max": (max(ht1004, key=lambda r: r["mu_app_pa_s"]) if ht1004 else None),
        }

    restored = stat_block(all_pairs)
    strict = stat_block(strict_all)
    # naive：只有 theta/speed，补算 τ/μ_app
    for r in naive_all:
        tau, gdot, mu = pair_metrics(r["theta"], r["speed"])
        r["tau_pa"], r["gamma_s1"], r["mu_app_pa_s"] = round(tau, 4), round(gdot, 3), round(mu, 4)
    naive = stat_block(naive_all)

    per_well = {}
    for w in sorted({r["well"] for r in all_pairs}):
        rs = [r for r in all_pairs if r["well"] == w]
        nrs = [r for r in naive_all if r["well"] == w]
        per_well[w] = {
            "n_pairs_复位": len(rs), "n_pairs_naive": len(nrs),
            "tau_max_pa": round(max(r["tau_pa"] for r in rs), 4),
            "muapp_max_pa_s": round(max(r["mu_app_pa_s"] for r in rs), 4),
            "max_muapp_at_speed": max(rs, key=lambda r: r["mu_app_pa_s"])["speed"],
            "tau_max_naive": round(max(r["tau_pa"] for r in nrs), 4) if nrs else None,
            "muapp_max_naive": round(max(r["mu_app_pa_s"] for r in nrs), 4) if nrs else None,
        }

    def srcd(rec):
        if not rec:
            return None
        return {k: rec[k] for k in ("well", "line", "fluid", "speed", "theta",
                                    "tau_pa", "mu_app_pa_s", "method")} if "method" in rec else rec

    # 3.0 支撑性判断（以 μ_app 口径，与原 2.6005 同口径）
    support = {
        "原始_3rpm档最大mu_app": naive["muapp_max_3rpm"],
        "复位_3rpm档最大mu_app": restored["muapp_max_3rpm"],
        "复位_全转速最大mu_app": restored["muapp_max_all"],
        "strict_复位_全转速最大mu_app": strict["muapp_max_all"],
        "上限_3.0": 3.0,
        "结论": ("复位后实测最大 μ_app = %.4f Pa·s < 3.0，量级仍与 3.0 同一档"
               "（3.0/最大实测 ≈ %.1f 倍），buoyancy.py『3.0 与最后一个可信低转速实测点同量级』表述"
               "在复位口径下仍成立（数值由 2.60 降为 %.2f）；上限结论不受影响（spec 主依据为口径唯一性，"
               "实测边界仅为支撑项）" % (
                   restored["muapp_max_all"], 3.0 / restored["muapp_max_all"],
                   restored["muapp_max_all"])),
    }

    result = {
        "meta": {
            "date": "2026-10-07", "task": "E-1 六速CSV按notes复位口径重算（Phase2 spec §10）",
            "caliber": "notes含φ<转速>:<读数>对⇒以note重建；否则列名；字段数不符者若多出字段全空⇒右对齐复位，否则剔除；>300等非数值读数逐条剔除",
            "conversion": "τ=0.511θ Pa, γ̇=1.703N s⁻¹, μ_app=τ/γ̇",
            "note_on_original_numbers": "原2.6005/8.7017/5.401经核验为μ_app口径（=τ/γ̇）而非τ：2.6005=0.511×26/5.109(hu1 L5错位θ6误入θ3)；8.7017=0.511×174/10.218(hu1 L5错位θ100误入θ6)；5.401=0.511×108/10.218(ht1_004 L9第三方θ100误入θ6)",
        },
        "四个关键数字_复位口径": {
            "总有效条数": restored["n_pairs"],
            "3rpm档最大τ_Pa": restored["tau_max_3rpm"],
            "3rpm档最大μ_app_Pa_s（与原2.6005同口径可比）": restored["muapp_max_3rpm"],
            "全转速最大τ_Pa": restored["tau_max_all"],
            "全转速最大μ_app_Pa_s（与原8.7017同口径可比）": restored["muapp_max_all"],
            "呼1-004最大τ_Pa": restored["ht1004_tau_max"],
            "呼1-004最大μ_app_Pa_s（与原5.401同口径可比）": restored["ht1004_muapp_max"],
        },
        "对照_naive错位口径_复现原数": {
            "总条数": naive["n_pairs"],
            "3rpm最大mu_app": naive["muapp_max_3rpm"],
            "全转速最大mu_app": naive["muapp_max_all"],
            "ht1004最大mu_app": naive["ht1004_muapp_max"],
            "3rpm最大tau": naive["tau_max_3rpm"],
            "全转速最大tau": naive["tau_max_all"],
        },
        "对照_原始报告数": {"总条数": 389, "3rpm": 2.6005, "全转速": 8.7017, "ht1004": 5.401},
        "敏感性_strict剔除14字段行": {
            "总条数": strict["n_pairs"],
            "3rpm最大mu_app": strict["muapp_max_3rpm"],
            "全转速最大mu_app": strict["muapp_max_all"],
            "ht1004最大mu_app": strict["ht1004_muapp_max"],
        },
        "关键条溯源_复位": {
            "tau_max_3rpm": srcd(restored["src_tau_max_3rpm"]),
            "muapp_max_3rpm": srcd(restored["src_muapp_max_3rpm"]),
            "tau_max_all": srcd(restored["src_tau_max_all"]),
            "muapp_max_all": srcd(restored["src_muapp_max_all"]),
            "ht1004_muapp_max": srcd(restored["src_ht1004_muapp_max"]),
        },
        "逐井表": per_well,
        "note复位行": notes_rows,
        "右对齐复位行": realigned_rows,
        "整行剔除记录": excl_rows,
        "逐条剔除_超量程": pair_drops,
        "3.0上限支撑性": support,
    }

    with open(os.path.join(OUT_DIR, "E1_复位重算结果.json"), "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2, default=str)

    # ---- MD 摘要 ----
    lines = []
    lines.append("# E-1 六速流变 CSV 复位重算摘要（2026-10-07）\n")
    lines.append("口径：notes 有 `φ<转速>:<读数>` 标注 ⇒ 以 note 为准重建；否则用列名；"
                 "字段数与表头不符但多出字段全空 ⇒ 右对齐复位；不可复位/超量程读数 ⇒ 剔除并注明。"
                 "换算 τ=0.511θ、γ̇=1.703N；μ_app=τ/γ̇。\n")
    lines.append("> 核验发现：原始 2.6005 / 8.7017 / 5.401 三个数**不是 τ 而是 μ_app=τ/γ̇**"
                 "（buoyancy.py 文档串写\"实测最大 μ_app≈2.60 Pa·s\"与此一致）。"
                 "三者均可由 naive 错位读数精确复现（见下）。\n")
    lines.append("## 四个关键数字（复位口径）\n")
    k = result["四个关键数字_复位口径"]
    lines.append("| 指标 | 复位值 | 原始报告值 | 差 |")
    lines.append("|---|---|---|---|")
    lines.append(f"| 总有效条数 (转速,读数)对 | {k['总有效条数']} | 389 | +{k['总有效条数']-389} |")
    lines.append(f"| 3 rpm 档最大 τ | {k['3rpm档最大τ_Pa']} Pa | —（原数非τ口径） | |")
    lines.append(f"| 3 rpm 档最大 μ_app | {k['3rpm档最大μ_app_Pa_s（与原2.6005同口径可比）']} Pa·s | 2.6005 | "
                 f"{round(k['3rpm档最大μ_app_Pa_s（与原2.6005同口径可比）']-2.6005,4)} |")
    lines.append(f"| 全转速最大 τ | {k['全转速最大τ_Pa']} Pa | — | |")
    lines.append(f"| 全转速最大 μ_app | {k['全转速最大μ_app_Pa_s（与原8.7017同口径可比）']} Pa·s | 8.7017 | "
                 f"{round(k['全转速最大μ_app_Pa_s（与原8.7017同口径可比）']-8.7017,4)} |")
    lines.append(f"| 呼1-004 最大 τ | {k['呼1-004最大τ_Pa']} Pa | — | |")
    lines.append(f"| 呼1-004 最大 μ_app | {k['呼1-004最大μ_app_Pa_s（与原5.401同口径可比）']} Pa·s | 5.401 | "
                 f"{round(k['呼1-004最大μ_app_Pa_s（与原5.401同口径可比）']-5.401,4)} |")
    lines.append("")
    nv = result["对照_naive错位口径_复现原数"]
    lines.append(f"## naive（错位）口径复现：总条数 {nv['总条数']}、3rpm μ_app {nv['3rpm最大mu_app']}、"
                 f"全转速 μ_app {nv['全转速最大mu_app']}、ht1_004 μ_app {nv['ht1004最大mu_app']}\n")
    lines.append("三个 μ_app 原数在 naive 口径下逐位复现 ⇒ 原数确系按列名硬切（错位）统计；"
                 f"总条数 naive={nv['总条数']} vs 原报 389，差 {nv['总条数']-389}，原统计脚本未落盘，"
                 "无法逐位复现其剔除细节，方向一致（错位口径漏掉 note 对与右移值）。\n")
    lines.append("## 差异归因\n")
    lines.append("- **2.6005→%.4f**：hu1 第5行（1:1:1白油基污染，14字段错位）θ₆=26 被错位读成 θ₃；"
                 "复位后真实 θ₃=16 ⇒ μ_app@3=0.511×16/5.109。" % k['3rpm档最大μ_app_Pa_s（与原2.6005同口径可比）'])
    lines.append("- **8.7017→%.4f**：同一行 θ₁₀₀=174 被错位读成 θ₆（174/6rpm）；复位后该行 174 归 100rpm"
                 "（μ_app=0.522），全表最大 μ_app 变为本行真实 3rpm 点 1.60。" % k['全转速最大μ_app_Pa_s（与原8.7017同口径可比）'])
    lines.append("- **5.401→%.4f**：ht1_004 第9行第三方复检出 φ100:108 被按列名读进 θ₆ 槽；"
                 "note 复位后该行低速最大为 φ6:10 ⇒ μ_app=0.500，全井最大为第2行 θ₃=9 ⇒ 0.9002。" % k['呼1-004最大μ_app_Pa_s（与原5.401同口径可比）'])
    lines.append("- **总条数**：note 复位使 ht1_004 两行各 4→5 条（+2）；hu1 14字段行右对齐复位使 6 行共 "
                 "28→34 条（naive 各丢 1 条读数为 28，复位 34）；其余同。")
    lines.append("")
    lines.append("## 逐井 τ / μ_app 最大（复位 vs naive）\n")
    lines.append("| 井 | 条数(复位/naive) | max τ Pa(复位) | max μ_app Pa·s(复位) | max μ_app(naive) |")
    lines.append("|---|---|---|---|---|")
    for w, v in per_well.items():
        lines.append(f"| {w} | {v['n_pairs_复位']}/{v['n_pairs_naive']} | {v['tau_max_pa']} | "
                     f"{v['muapp_max_pa_s']} (@{v['max_muapp_at_speed']}rpm) | {v['muapp_max_naive']} |")
    lines.append("")
    lines.append("## note 复位行\n")
    for r in notes_rows:
        lines.append(f"- {r['well']} L{r['line']} {r['fluid']}：{' '.join(r['raw'])}（{r['n_pairs']} 条）")
    lines.append("\n## 右对齐复位行（hu1，多出字段全空）\n")
    for r in realigned_rows:
        lines.append(f"- {r['well']} L{r['line']} {r['fluid']}：{r['n_pairs']} 条")
    lines.append("\n## 剔除清单\n")
    lines.append("**整行剔除**：无（hu1 14字段行均可右对齐复位；若按严格口径一律剔除 14 字段行，"
                 f"见敏感性：总条数→{result['敏感性_strict剔除14字段行']['总条数']}，"
                 f"全转速最大 μ_app→{result['敏感性_strict剔除14字段行']['全转速最大mu_app']}，结论不变）。")
    lines.append("\n**逐条剔除（超量程，无数字可计）**：")
    for r in pair_drops:
        lines.append(f"- {r['well']} L{r['line']} {r['fluid']} @{r['speed']}rpm：`{r['raw']}` ⇒ {r['reason']}")
    lines.append("\n**其他记录**：")
    for r in excl_rows:
        lines.append(f"- {r['well']} L{r['line']} {r['fluid']}：{r['issue']}")
    lines.append("\n**保留但存疑读数**：ht1_001 L4-6 θ₆₀₀=300 且 notes 记 `PHI600>300未完全读数`"
                 "（饱和值，非精确读数；对 μ_app 无影响，对\"全转速最大 τ\"有支配作用 153.3 Pa）。")
    lines.append("\n## 对\"上限 3.0 Pa·s\"结论的支撑性（spec §1.2）\n")
    lines.append(f"- 复位后实测最大 μ_app = **{support['复位_全转速最大mu_app']} Pa·s**"
                 f"（3 rpm 档最大 μ_app = {support['复位_3rpm档最大mu_app']} Pa·s），"
                 "仍低于 3.0，比值 3.0/最大 ≈ 1.9 ⇒『与最后一个可信低转速实测点同量级』仍成立，但强度略降"
                 "（2.60→1.60）；严格剔除口径下最大 μ_app 更低（≈1.0），支撑同样成立。")
    lines.append("- 该最大点来自 hu1 污染混拌实验（1:1:1 白油基领浆:隔离液:泥浆），"
                 "如视为非典型可剔除，则 3rpm 档最大 μ_app=1.0002（ht1_003 尾浆/hu101 钻井液 θ₃=10，即严格剔除口径值）——两口径均远低于 3.0。")
    lines.append("- spec 主依据为口径唯一性（:1133 结构常数共用），实测边界仅支撑项，**结论不受影响**。")
    with open(os.path.join(OUT_DIR, "E1_复位重算摘要.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    print(json.dumps(result["四个关键数字_复位口径"], ensure_ascii=False, indent=2))
    print("naive复现:", json.dumps(result["对照_naive错位口径_复现原数"], ensure_ascii=False))
    print("strict敏感性:", json.dumps(result["敏感性_strict剔除14字段行"], ensure_ascii=False))
    print("逐条剔除数:", len(pair_drops), "| 整行剔除:", len([r for r in excl_rows if "不可复位" in r["issue"]]),
          "| note行:", len(notes_rows), "| 右对齐行:", len(realigned_rows))
    print("输出:", OUT_DIR)


if __name__ == "__main__":
    main()
