"""开关关逐位验证：hu103 全事件 vs 基线锚 JSON"""
import json
from cemdisp.data.loaders import load_hu103_tailpipe
from cemdisp.transport1d.casing_flow import CasingFlowSolver
from cemdisp.transport1d.shoe_timeline import ShoeEventKind

well, fluids, schedule, _ = load_hu103_tailpipe()
solver = CasingFlowSolver(enable_gravity=True)
result = solver.run(well, fluids, schedule)
with open("tests/_baseline_shoe_timeline_hu103.json", encoding="utf-8") as f:
    anchor = json.load(f)
ev = result.shoe_timeline.events
base = anchor["events"]
assert len(ev) == len(base) == anchor["_meta"]["event_count"], (len(ev), len(base))
mism = 0
for i, (a, b) in enumerate(zip(ev, base)):
    ok = (a.time_s == b["time_s"] and a.kind.value == b["kind"]
          and a.flow_rate_m3_s == b["flow_rate_m3_s"] and a.stage_name == b["stage_name"]
          and list(map(list, a.phase_fractions)) == b["phase_fractions"])
    if not ok:
        mism += 1
        print(f"MISMATCH event {i}: {a} vs {b}")
print(f"events={len(ev)} mismatch={mism} cement_end={result.cement_end_time_s} anchor={anchor['_meta']['cement_end_time_s']}")
print("BITWISE OFF-PATH:", "PASS" if mism == 0 and result.cement_end_time_s == anchor["_meta"]["cement_end_time_s"] else "FAIL")
