# -*- coding: utf-8 -*-
"""默认路径逐位锚的**重锚脚本**（R25）。

存在理由
--------
主测试 `test_default_path_bitwise_anchor.py` 在数值位移时的失败提示会指向本脚本。
在本次修复之前该提示指向一个**不存在的文件**——维护者照做会得到 `can't open file`，
随后最自然的动作就是"手工改 JSON 里的数值"，而那正是锚加固要防的路径：
防线在最需要它的那一刻给出错误指引。故本脚本必须真实存在且可跑。

用法
----
    python tests/contract/_regenerate_default_path_anchor.py            # 干跑：只打印差异
    python tests/contract/_regenerate_default_path_anchor.py --confirm  # 确认后写盘

- 无 `--confirm` 时**只打印每个键的 before/after，绝不写任何文件**；
- 有 `--confirm` 时写盘，且**保留** `_note` / `_generated_from` / `_env` 三键、
  把 `_generated_from` 里引用的生成批次更新为当前 HEAD（git 不可用则写 `unknown`）；
- 重锚后**必须**把本脚本输出贴进
  `results/内部自洽加固_2026-09-25/位移台账.csv` 留迹（本脚本只提示，不写 results/）。

退出码：0=无差异或写盘成功；1=干跑发现差异（便于 CI/脚本判"该重锚了"）。
"""
import argparse
import json
import sys
from pathlib import Path

# 允许 `python tests/contract/_regenerate_default_path_anchor.py` 直接运行：
# 运行脚本时 sys.path[0] 是脚本所在目录（tests/contract/），仓库根不在路径上，
# 故把仓库根插到最前，使 `tests.contract.*` 可导入（tests/contract/__init__.py 存在）。
_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

# 优先走包导入（可享受 pytest 同一套导入语义）；失败再按文件路径加载（R25 要求的回退）。
try:
    from tests.contract.test_default_path_bitwise_anchor import (  # noqa: E402
        ANCHOR, FINGERPRINT_KEYS, POINTER_KEYS, _run_default_case,
        build_fingerprint, current_head_short, expected_keys)
except ImportError:  # pragma: no cover —— 仅在包导入不可用时走到
    import importlib.util

    _MAIN = Path(__file__).with_name("test_default_path_bitwise_anchor.py")
    _spec = importlib.util.spec_from_file_location("_anchor_main_module", _MAIN)
    _mod = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(_mod)
    ANCHOR = _mod.ANCHOR
    FINGERPRINT_KEYS = _mod.FINGERPRINT_KEYS
    POINTER_KEYS = _mod.POINTER_KEYS
    _run_default_case = _mod._run_default_case
    build_fingerprint = _mod.build_fingerprint
    current_head_short = _mod.current_head_short
    expected_keys = _mod.expected_keys

LEDGER_HINT = "results/内部自洽加固_2026-09-25/位移台账.csv"


def _brief(key: str, value: object) -> str:
    """短展示：sha 只给前 12 字符，其余原样。"""
    if key.startswith("sha_") and isinstance(value, str) and len(value) > 12:
        return f"{value[:12]}…（len={len(value)}）"
    return repr(value)


def _new_generated_from(old: str, head: str) -> str:
    """把 `_generated_from` 里引用的生成批次更新为当前 HEAD，其余说明文字保留。

    Task 10 写的值形如 `HEAD bc155b0 + 4 个未提交 HB 改动（工作树脏态）；生成命令见 …`；
    这里只替换 sha 段（正则找到第一个"HEAD <sha7位左右>"）。
    """
    import re
    if re.search(r"HEAD [0-9a-f]{7,40}", old):
        return re.sub(r"HEAD [0-9a-f]{7,40}", f"HEAD {head}", old, count=1)
    return f"HEAD {head}；原指纹：{old}"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="重锚默认路径逐位锚（无 --confirm 时只打印差异，不写盘）")
    parser.add_argument("--confirm", action="store_true",
                        help="确认写盘（须已人工核对差异，并在位移台账留迹）")
    args = parser.parse_args(argv)

    print(f"[锚文件] {ANCHOR}")
    if not ANCHOR.is_file():
        print(f"[错误] 锚文件不存在：{ANCHOR}")
        print(f"       首次写锚可加 --confirm 由本脚本创建。")
    else:
        print(f"[旧锚] 已读取（{ANCHOR.stat().st_size} 字节）")

    print("[实测] 正在跑默认路径端到端算例（微网格，约 10s）…")
    measured = _run_default_case()

    old: dict = {}
    if ANCHOR.is_file():
        old = json.loads(ANCHOR.read_text(encoding="utf-8"))

    # 位移只算**实测键**（`expected_keys()`）；指纹键与可解析指针键由本脚本
    # 从旧锚继承后更新，不参与实测比对（否则会误报成"位移"）。
    keys = sorted(expected_keys() | set(POINTER_KEYS) | set(FINGERPRINT_KEYS))
    diffs = [k for k in sorted(expected_keys())
             if old.get(k) != measured.get(k)]

    print("\n=== 逐键 before/after ===")
    for k in keys:
        if k in FINGERPRINT_KEYS:
            continue
        before, after = old.get(k, "<缺>"), measured.get(k)
        if k in POINTER_KEYS:
            print(f"·  {k:22s} = {before!r}（指针键，非实测值）")
        elif before == after:
            print(f"   {k:22s} = {_brief(k, before)}")
        else:
            print(f"≠  {k:22s} before = {_brief(k, before)}")
            print(f"   {' ' * 22} after  = {_brief(k, after)}")

    if not diffs:
        print("\n[结论] 与实测逐位一致：无需重锚（未写任何文件）。")
        return 0

    head = current_head_short()
    print(f"\n[结论] 发现 {len(diffs)} 个键位移：{diffs}")
    print(f"[指纹] 当前 HEAD：{head}")

    if not args.confirm:
        print("\n[干跑] 未加 --confirm ⇒ **不写任何文件**。")
        print(f"       如确认重锚，请加 --confirm，并在 {LEDGER_HINT} 留迹"
              "（重锚须与改动同批次提交）。")
        return 1

    fingerprint = build_fingerprint(
        env=old.get("_env", "unknown"),
        generated_from=_new_generated_from(old.get("_generated_from", ""), head),
        reconcile_source=old.get("_env_reconcile_source", ""))
    new_anchor = {**fingerprint, **measured}
    ANCHOR.write_text(
        json.dumps(new_anchor, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8")
    print(f"\n[写盘] 已写入 {ANCHOR}（保留并更新指纹键：{sorted(fingerprint)}）")
    print(f"[留迹] 请把以上差异贴进 {LEDGER_HINT}，并与改动同批次提交。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
