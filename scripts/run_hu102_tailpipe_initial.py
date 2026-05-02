"""呼102尾管段模型运行入口 — 薄wrapper，核心逻辑在 cemdisp.runners.hu102_tailpipe"""

from cemdisp.runners.hu102_tailpipe import run_hu102_tailpipe_initial

if __name__ == "__main__":
    run_hu102_tailpipe_initial()