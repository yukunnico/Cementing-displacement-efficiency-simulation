"""呼探1井尾管段模型运行入口 — 薄wrapper，核心逻辑在 cemdisp.runners.hu1_tailpipe"""

from cemdisp.runners.hu1_tailpipe import run_hu1_tailpipe_initial

if __name__ == "__main__":
    run_hu1_tailpipe_initial()