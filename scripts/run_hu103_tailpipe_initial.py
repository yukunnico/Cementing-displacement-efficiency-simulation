"""呼103尾管段模型运行入口 — 薄wrapper，核心逻辑在 cemdisp.runners.hu103_tailpipe"""

from cemdisp.runners.hu103_tailpipe import run_hu103_tailpipe_initial

if __name__ == "__main__":
    run_hu103_tailpipe_initial()