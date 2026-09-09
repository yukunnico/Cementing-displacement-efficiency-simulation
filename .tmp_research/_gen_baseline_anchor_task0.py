# -*- coding: utf-8 -*-
"""Task 0 基线锚生成脚本（一次性）：hu103 shoe_timeline.events 冻结为 JSON。

产出 tests/_baseline_shoe_timeline_hu103.json，供 Task 1 逐位回归对比。
格式：
{
  "_meta": {...},
  "events": [ {"time_s":..., "kind":..., "flow_rate_m3_s":..., "stage_name":..., "phase_fractions":[[name,frac],...]}, ... ]
}
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from cemdisp.data.loaders import load_hu103_tailpipe
from cemdisp.transport1d import CasingFlowSolver

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "tests" / "_baseline_shoe_timeline_hu103.json"

well_spec, fluids, schedule, _validation = load_hu103_tailpipe()

solver = CasingFlowSolver(enable_gravity=True)
res = solver.run(well_spec, fluids, schedule)
events = res.shoe_timeline.events

payload = {
    "_meta": {
        "purpose": "Task0 基线冻结锚：hu103 shoe_timeline.events 逐位回归基准（路线B Task1 用）",
        "well_name": well_spec.well_name,
        "loader": "cemdisp.data.loaders.load_hu103_tailpipe",
        "loader_kwargs": {},
        "solver_class": "cemdisp.transport1d.CasingFlowSolver",
        "solver_params": {
            "enable_gravity": True,
            "其余参数全部默认": True,
        },
        "event_count": len(events),
        "cement_end_time_s": res.cement_end_time_s,
        "baseline_commit": "2190bb3",
        "branch": "feat/casing-mixing-contact-time-routeB",
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "format_notes": (
            "events 为按时间升序的 ShoeEvent 序列；kind 取 ShoeEventKind.value 字符串；"
            "phase_fractions 序列化为 [[name, frac], ...]（原为 tuple of (name, frac)）；"
            "float 全部 JSON 原样精度（repr 双射），Task1 对比用 == 或 pytest.approx(rel=1e-12)"
        ),
    },
    "events": [
        {
            "time_s": ev.time_s,
            "kind": ev.kind.value,
            "flow_rate_m3_s": ev.flow_rate_m3_s,
            "stage_name": ev.stage_name,
            "phase_fractions": [[name, frac] for name, frac in ev.phase_fractions],
        }
        for ev in events
    ],
}

OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
print(f"[done] {len(events)} events -> {OUT}")
print("first:", json.dumps(payload["events"][0], ensure_ascii=False))
print("last :", json.dumps(payload["events"][-1], ensure_ascii=False))
