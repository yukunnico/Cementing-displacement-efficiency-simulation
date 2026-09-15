# STATUS: contract
"""zhang2022_benchmark main() 权威目录防覆写守卫。

终审 I-2（2026-09-15）：``main()`` 未传 ``--out-dir`` 时必须报错退出（SystemExit），
绝不回落 ``DEFAULT_OUTPUT_DIR``（``results/基准算例对照_2026-09-10``，R18 权威目录）。
本测试只断言报错与权威目录无写入，不真跑算例、不真覆写。
"""

from __future__ import annotations

import pytest

from cemdisp.runners.zhang2022_benchmark import DEFAULT_OUTPUT_DIR, main


def _snapshot():
    """权威目录文件清单+mtime 快照；目录不存在时返回 None。"""
    if not DEFAULT_OUTPUT_DIR.exists():
        return None
    return sorted(
        (str(p.relative_to(DEFAULT_OUTPUT_DIR)), p.stat().st_mtime)
        for p in DEFAULT_OUTPUT_DIR.rglob("*")
        if p.is_file()
    )


def test_main_without_out_dir_exits_and_never_touches_authoritative_dir(capsys):
    """无 --out-dir 调 main 必须抛 SystemExit(2)，且权威目录逐位无写入。"""
    before = _snapshot()
    with pytest.raises(SystemExit) as exc_info:
        main([])
    # argparse parser.error 的约定退出码
    assert exc_info.value.code == 2
    # 报错信息须指引用户显式传 --out-dir
    captured = capsys.readouterr()
    assert "--out-dir" in captured.err
    # 权威目录零改动（文件集合与 mtime 均不变；目录不存在也不得被创建）
    assert _snapshot() == before
