"""
分析窄边水泥浓度剖面和鞋口出流时序（Step 3.3 & 3.4）
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

_PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

# 设置中文字体
plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False


def analyze_well(well_key: str):
    path = _PROJECT_ROOT / "results" / "p3_p4_integration" / f"{well_key}_profiles.json"
    data = json.loads(path.read_text(encoding="utf-8"))

    baseline = data["baseline"]
    improved = data["improved"]

    depth = np.array(baseline["profile"]["depth_m"])

    # Step 3.3: 窄边水泥浓度剖面对比图
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    for ax, label, result in [
        (axes[0], "基准配置 (M=0, 无弥散)", baseline),
        (axes[1], "改进配置 (M=100, 有弥散)", improved),
    ]:
        profile = result["profile"]
        ax.plot(profile["wide_cement"], depth, label="宽边水泥浓度", linewidth=2)
        ax.plot(profile["narrow_cement"], depth, label="窄边水泥浓度", linewidth=2)
        ax.plot(profile["wide_eff"], depth, label="宽边有效效率", linewidth=2, linestyle="--")
        ax.plot(profile["narrow_eff"], depth, label="窄边有效效率", linewidth=2, linestyle="--")
        ax.set_xlabel("浓度 / 效率")
        ax.set_ylabel("井深 (m)")
        ax.set_title(f"{well_key.upper()} {label}\nCBL效率={result['cbl_eval_eff']:.4f}")
        ax.invert_yaxis()
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    fig.tight_layout()
    out_png = _PROJECT_ROOT / "results" / "p3_p4_integration" / f"{well_key}_窄边水泥浓度剖面对比.png"
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"Saved {out_png}")

    # Step 3.4: 鞋口出流时序对比
    fig, ax = plt.subplots(figsize=(10, 5))

    # 基准时序
    base_events = baseline["shoe_events"]
    times_b = [ev["time_s"] / 60.0 for ev in base_events]
    phases_b = []
    for ev in base_events:
        pf = ev["phase_fractions"]
        if pf:
            phases_b.append(pf[0][0])
        else:
            phases_b.append("未知")

    # 改进时序
    imp_events = improved["shoe_events"]
    times_i = [ev["time_s"] / 60.0 for ev in imp_events]
    # 计算水泥相分数
    cement_frac_i = []
    for ev in imp_events:
        pf = ev["phase_fractions"]
        frac = 0.0
        for name, f in pf:
            if "水泥" in name or name == "cement":
                frac += f
        cement_frac_i.append(frac)

    ax.step(times_b, range(len(times_b)), where="post", label="基准 (离散阶跃)", linewidth=2)
    ax.plot(times_i, cement_frac_i, "o-", label="改进 (水泥相分数)", linewidth=2, markersize=4)
    ax.set_xlabel("地面累计时间 (min)")
    ax.set_ylabel("事件序号 / 水泥相分数")
    ax.set_title(f"{well_key.upper()} 鞋口出流时序对比")
    ax.legend()
    ax.grid(True, alpha=0.3)

    out_png2 = _PROJECT_ROOT / "results" / "p3_p4_integration" / f"{well_key}_鞋口出流时序对比.png"
    fig.savefig(out_png2, dpi=200)
    plt.close(fig)
    print(f"Saved {out_png2}")

    # 打印关键统计
    print(f"\n=== {well_key.upper()} 窄边水泥浓度统计 ===")
    for label, result in [("基准", baseline), ("改进", improved)]:
        nc = np.array(result["profile"]["narrow_cement"])
        wc = np.array(result["profile"]["wide_cement"])
        print(f"  {label}: 窄边均值={nc.mean():.4f}, 宽边均值={wc.mean():.4f}, 宽窄比={nc.mean()/wc.mean():.4f}")

    print(f"\n=== {well_key.upper()} 鞋口出流事件数 ===")
    print(f"  基准: {len(base_events)} 个事件")
    print(f"  改进: {len(imp_events)} 个事件")
    if len(imp_events) > len(base_events):
        print(f"  弥散增加了 {len(imp_events) - len(base_events)} 个过渡事件")


if __name__ == "__main__":
    for well in ["hu101", "hu102"]:
        analyze_well(well)
