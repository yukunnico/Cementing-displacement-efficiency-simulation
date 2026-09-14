"""2026-09-02 阶段5：定位“残余损耗 + 网格不收敛 + 关弥散崩塌”的共同病根。
对 hu102（及 hu101）在 nz=120/250 下跑若干配置，比较 η_E、守恒率、前缘、壁面冻结。
假设：不可逆残余壁面冻结(cement_ever>0 且 c<c_min 永久冻结)在尖锐前缘自我扩张卡住前缘，
人工无量纲拉普拉斯弥散只是把它糊住，且其等效扩散系数 D∝dz，加密反而趋近崩塌。"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.lib.mass_balance_diag import WELLS, run_variant

VARIANTS = {
    "BASE": dict(),
    "关壁面冻结": dict(c_min=0.0),
    "关壁面冻结+关弥散": dict(c_min=0.0, dispersion_axial=0.0, dispersion_azimuthal=0.0),
    "关壁面冻结+关弥散+关D2DGA": dict(c_min=0.0, dispersion_axial=0.0, dispersion_azimuthal=0.0,
                                  enable_d2dga=False, enable_true_buoyancy=False, enable_d2dga_i3_flux=False),
    "关壁面冻结+关弥散+近同心": dict(c_min=0.0, dispersion_axial=0.0, dispersion_azimuthal=0.0, e_clip_max=0.05),
}

def main(wells, nzs):
    allr = {}
    for w in wells:
        allr[w] = {}
        for nz in nzs:
            allr[w][nz] = {}
            for vn, kw in VARIANTS.items():
                allr[w][nz][vn] = run_variant(f"nz{nz}_{vn}", WELLS[w], nz, **kw)
    out = Path(__file__).resolve().parents[2]/"results"/"_质量平衡取证_2026-09-02"/"阶段5_壁面冻结网格机制.json"
    out.write_text(json.dumps(allr, ensure_ascii=False, indent=2), encoding="utf-8")

if __name__ == "__main__":
    wells = sys.argv[1].split(",") if len(sys.argv)>1 else ["hu102"]
    nzs = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv)>2 else [120,250]
    main(wells, nzs)
