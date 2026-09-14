"""阶段6：验证物理τw屈服门(去residual_wall)的网格收敛性与物理窜槽保留。"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.lib.mass_balance_diag import WELLS, run_variant

V = {
  "物理屈服门": dict(enable_yield_gate=True),
  "对照_浓度冻结基线": dict(),
  "对照_完全无壁面层": dict(c_min=0.0),
}
res = {}
for w in (sys.argv[1].split(",") if len(sys.argv)>1 else ["hu102","hu101"]):
    res[w] = {}
    for nz in ([int(x) for x in sys.argv[2].split(",")] if len(sys.argv)>2 else [120,250]):
        res[w][nz] = {vn: run_variant(vn, WELLS[w], nz, **kw) for vn,kw in V.items()}
out = Path(__file__).resolve().parents[2]/"results"/"_质量平衡取证_2026-09-02"/"阶段6_物理屈服门验证.json"
out.write_text(json.dumps(res, ensure_ascii=False, indent=2), encoding="utf-8")
