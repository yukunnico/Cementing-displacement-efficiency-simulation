"""呼探1-001井尾管段模型运行入口 — 薄wrapper，核心逻辑在 cemdisp.runners.ht1_001_tailpipe"""

from cemdisp.runners.ht1_001_tailpipe import run_ht1_001_tailpipe_initial

if __name__ == "__main__":
    run_ht1_001_tailpipe_initial()