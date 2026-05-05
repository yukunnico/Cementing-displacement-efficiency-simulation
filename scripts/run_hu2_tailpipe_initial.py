"""呼探1-002井尾管段模型运行入口 — 薄wrapper，核心逻辑在 cemdisp.runners.hu2_tailpipe"""

from cemdisp.runners.hu2_tailpipe import run_hu2_tailpipe_initial

if __name__ == "__main__":
    run_hu2_tailpipe_initial()