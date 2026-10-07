# -*- coding: utf-8 -*-
"""Phase 4d 库侧验证：真实扩展表经 load_extended_pair_4d 加载 + stop_t 覆盖断言。

对照：交付表上 assert_time_table_coverage(r0.6/r0.8 stop_t) 必须响亮报错（旧污染场景）；
扩展表上必须放行且查询零 clamp（oob_count==0）。
"""
import json, sys
from pathlib import Path
BRANCH = Path(r"D:\users\desktop\research\控压固井项目\cement model_温压耦合分支")
sys.path.insert(0, str(BRANCH))  # editable 安装指向冻结源目录，必须前插
import numpy as np
from cemdisp.data.temperature_field import (
    TemperatureTableCoverageError, TableTemperatureField,
    assert_time_table_coverage, load_delivered_pair, load_extended_pair_4d,
)

SB = BRANCH / "results" / "_probe_4d扩表_2026-10-07" / "sandbox" / "HT1-004压力计算"
out = {}

# 1) 真实扩展表加载
f_in, f_out = load_extended_pair_4d(
    SB / "T_in_ext4d.xlsx", SB / "T_out_ext4d.xlsx", SB / "T_ext4d_time_axis_min.csv",
)
out["extended_load"] = {
    "shape": list(f_in.table.shape),
    "axis_end_s": float(f_in.time_s[-1]),
    "axis_col199_s": float(f_in.time_s[199]),
    "first200_equals_delivered_bitwise": None,  # 由 compare_result.json 判定位级；此处校验 npz 一致
}

# 2) 前 200 列与交付表位级一致（经库两条路径各取一次）
del_in, del_out = load_delivered_pair(use_cache=False)
out["extended_load"]["first200_equals_delivered_bitwise"] = bool(
    np.array_equal(f_in.table[:, :200], del_in.table)
    and np.array_equal(f_out.table[:, :200], del_out.table)
)

# 3) 交付表上的旧污染场景：断言必须响亮报错
stop = {"r0.6": 19016.017281054297, "r0.8": 14262.0, "r1.0": 11409.6, "r1.2": 9508.0}
raises, passes = {}, {}
for k, v in stop.items():
    try:
        assert_time_table_coverage([del_in, del_out], v, label=f"delivered_{k}")
        raises[k] = False
    except TemperatureTableCoverageError:
        raises[k] = True
    assert_time_table_coverage([f_in, f_out], v, label=f"ext4d_{k}")
    passes[k] = True
out["assert_on_delivered_raises"] = raises
out["assert_on_extended_passes"] = passes

# 4) r0.6 stop_t 处采样查询：扩展表零 clamp，交付表仍 clamp+计数（库行为不变）
t_q = stop["r0.6"]
v_ext = f_in.T(7660.0 * 0.99, t_q)  # 井底附近深度（域内）
out["query_at_r06_stop"] = {
    "extended_value_C": v_ext,
    "extended_oob": f_in.oob_count,
    "delivered_value_C": del_in.T(7660.0 * 0.99, t_q),
    "delivered_oob": del_in.oob_count,
    "note": "扩展表 oob=0（真覆盖）；交付表 oob=1（clamp 计数行为保持原样）",
}

with open(Path(__file__).with_name("validate_result.json"), "w", encoding="utf-8") as fp:
    json.dump(out, fp, ensure_ascii=False, indent=1)
print(json.dumps(out, ensure_ascii=False, indent=1))
