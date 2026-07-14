"""HT1-004 R0->R3 消融运行入口。产出 results/ht1_004_ablation/ 下的各级摘要。"""
from __future__ import annotations

from cemdisp.runners.ht1_004_ablation import (
    run_full_ablation,
    extract_ablation_metrics,
)

if __name__ == "__main__":
    results = run_full_ablation(
        nz=500,
        dt=4.0,
        total_t=None,
        output_dir="results/ht1_004_ablation",
    )
    for name, res in results.items():
        m = extract_ablation_metrics(res)
        eff = m.get("effective_efficiency")
        print(f"{name}: effective_efficiency={eff}")