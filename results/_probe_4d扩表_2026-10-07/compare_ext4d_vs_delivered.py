# -*- coding: utf-8 -*-
"""Phase 4d 扩表对账：沙箱扩展表(333x362) vs 交付表(333x200)。

硬判据：新表前 200 列与交付表逐元素一致（同深度行、同时间列，网格完全同款：
0:1:198 min + 施工终点 198.792801 min，其后 199:1:360 min 为扩展段）。
>11940 s 段行为如实报告（由 HT1_004_T.m 物理决定：界面冻结+末阶段排量续循环）。
"""
import json
import sys
from pathlib import Path

import numpy as np

BRANCH = Path(r"D:\users\desktop\research\控压固井项目\cement model_温压耦合分支")
sys.path.insert(0, str(BRANCH))  # editable 安装指向冻结源目录，必须显式前插

from cemdisp.data.temperature_field import _read_xlsx_numeric, _read_depth_axis  # noqa: E402

SB = BRANCH / "results" / "_probe_4d扩表_2026-10-07" / "sandbox" / "HT1-004压力计算"
DELIV = BRANCH / "参考文档" / "温压耦合数据、"

out = {}

# ---- 读入 ----
t_in_ref = _read_xlsx_numeric(DELIV / "T_in.xlsx")
t_out_ref = _read_xlsx_numeric(DELIV / "T_out.xlsx")
t_in_ext = _read_xlsx_numeric(SB / "T_in_ext4d.xlsx")
t_out_ext = _read_xlsx_numeric(SB / "T_out_ext4d.xlsx")
axis_min = np.loadtxt(SB / "T_ext4d_time_axis_min.csv", dtype=float)
depth = _read_depth_axis(SB / "呼1-004井身结构.csv")

out["shapes"] = {
    "delivered_T_in": list(t_in_ref.shape), "delivered_T_out": list(t_out_ref.shape),
    "ext_T_in": list(t_in_ext.shape), "ext_T_out": list(t_out_ext.shape),
    "axis_len": int(axis_min.size), "depth_len": int(depth.size),
}
assert t_in_ext.shape[1] == axis_min.size == 362

# ---- 时间轴核验 ----
total_time_min = 198.792801055535  # ext_scalars.txt（= sum(stage_volume/Q)）
ref_axis = np.append(np.arange(199, dtype=float), total_time_min)  # 交付表真实轴 0..198,198.7928
ext_first200 = axis_min[:200]
out["axis"] = {
    "first200_equals_delivered_true_axis": bool(np.array_equal(ext_first200, ref_axis)),
    "axis_max_s": float(axis_min[-1] * 60.0),
    "delivered_true_end_s": float(total_time_min * 60.0),
    "monotonic": bool(np.all(np.diff(axis_min) > 0)),
    "steps_after_col199_min": [float(x) for x in np.diff(axis_min[199:203])],
}

# ---- 硬判据：前 200 列逐元素 ----
for name, ref, ext in (("T_in", t_in_ref, t_in_ext), ("T_out", t_out_ref, t_out_ext)):
    d = np.abs(ext[:, :200] - ref)
    out[name] = {
        "n_exact": int(np.sum(d == 0.0)),
        "n_gt_1e-12": int(np.sum(d > 1e-12)),
        "n_gt_1e-9": int(np.sum(d > 1e-9)),
        "max_abs_diff": float(d.max()),
        "mean_abs_diff": float(d.mean()),
    }
    # 扩展段（>11940 s 名义 / >198.79 min 真实）物理行为
    ext_seg = ext[:, 200:]
    out[name]["extension_segment"] = {
        "col200_end_s_note": "col199=施工终点198.7928min；col200 起为 199..360 min 扩展段",
        "depth_last_row_over_ext": {  # 井底温度演变（持稳/降温/升温）
            "at_199min": float(ext[-1, 200]),
            "at_250min": float(ext[-1, 251]),
            "at_300min": float(ext[-1, 301]),
            "at_360min": float(ext[-1, 361]),
        },
        "depth_30m_row_over_ext": {
            "at_199min": float(ext[0, 200]),
            "at_360min": float(ext[0, 361]),
        },
        "colmonotone_last_row": bool(np.all(np.diff(ext[-1, 200:]) < 0)),
    }

# ---- 井底（最深行）两表同值（交付约定在扩展表延续；注意只比最后一行） ----
d_bottom = np.abs(t_in_ext[-1, :] - t_out_ext[-1, :])
out["bottom_row_shared"] = bool(d_bottom.max() == 0.0)
out["bottom_row_max_abs_diff"] = float(d_bottom.max())

# ---- stop_t 覆盖评估（解禁条件） ----
axis_end_s = float(axis_min[-1] * 60.0)
stop_t = {  # 来源：results/敏感性变体_温压T2_2026-10-01/呼1-004*_判别量.json
    "r0.6_Toff": 19016.017281054297,
    "r0.6_Ton_const60": 19016.017281054297,
    "r0.6_Ton_static": 19016.017281054297,
    "r0.6_Ton_table": 19016.017281054297,
    "r0.8_Toff": 14262.0,
    "r0.8_Ton_static": 14262.0,
    "r1.0_baseline": 11409.6,
}
out["stop_t_coverage"] = {
    k: {"stop_t_s": v, "covered_by_ext": bool(v <= axis_end_s),
        "covered_by_delivered": bool(v <= 11940.0),
        "margin_s_vs_ext_end": axis_end_s - v}
    for k, v in stop_t.items()
}
out["verdict"] = {
    "history_unchanged_T_in": out["T_in"]["n_gt_1e-9"] == 0,
    "history_unchanged_T_out": out["T_out"]["n_gt_1e-9"] == 0,
    "r06_covered": axis_end_s >= max(stop_t["r0.6_Toff"], 0),
    "r08_covered": axis_end_s >= stop_t["r0.8_Toff"],
}

with open(Path(__file__).with_name("compare_result.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print(json.dumps(out["verdict"], ensure_ascii=False))
print("T_in : max|Δ|前200列 =", out["T_in"]["max_abs_diff"], " 精确等=", out["T_in"]["n_exact"], "/ 66600")
print("T_out: max|Δ|前200列 =", out["T_out"]["max_abs_diff"], " 精确等=", out["T_out"]["n_exact"], "/ 66600")
print("表末 =", axis_end_s, "s; r0.6 stop_t 覆盖 =", out["stop_t_coverage"]["r0.6_Toff"]["covered_by_ext"])
