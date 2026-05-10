"""呼101尾管段模型运行入口 — 薄wrapper，核心逻辑在 cemdisp.runners.hu101_tailpipe"""

from cemdisp.runners.hu101_tailpipe import run_hu101_tailpipe_initial

if __name__ == "__main__":
    run_hu101_tailpipe_initial()
