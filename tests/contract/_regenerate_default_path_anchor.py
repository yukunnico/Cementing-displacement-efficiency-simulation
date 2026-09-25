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
- 写盘序列化**固定 `sort_keys=True`**（见 `serialize_anchor`）：合法重锚不得重排键序，
  否则 diff 会把每个数值行写成"删+加"，抹掉"diff 里没有数值行被改动"这一
  可核对形式（R33）；
- 锚文件**不存在**（首建分支）时，**必须**显式给出 `--env` 与 `--reconcile-source`，
  否则本脚本**拒绝写盘**并以用法错误退出（码 2）——否则会写出
  `_env="unknown"`、指针键为空的锚，当场被 `test_anchor_integrity.py` 的指纹
  守卫判红，逼维护者手工改 JSON，正是锚加固要消灭的路径（R34）；
- 重锚后**必须**把本脚本输出贴进
  `results/内部自洽加固_2026-09-25/位移台账.csv` 留迹（本脚本只提示，不写 results/）。

退出码：0=无差异或写盘成功；1=干跑发现差异（便于 CI/脚本判"该重锚了"）；
2=首建缺指纹参数等用法错误（未写任何文件）。
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
    首建锚时 `old` 为空 ⇒ 直接写 `HEAD <sha>`（不留半截"原指纹："）。
    """
    import re
    if not (old or "").strip():
        return f"HEAD {head}"
    if re.search(r"HEAD [0-9a-f]{7,40}", old):
        return re.sub(r"HEAD [0-9a-f]{7,40}", f"HEAD {head}", old, count=1)
    return f"HEAD {head}；原指纹：{old}"


def serialize_anchor(anchor: dict) -> str:
    """锚的**唯一**序列化口径（R33）。

    `sort_keys=True` 是硬要求：`{**fingerprint, **measured}` 的插入序是
    "指纹键在前、实测键按 `_run_default_case()` 的产出序在后"，与锚文件里
    `sorted()` 的键序不同。少了 `sort_keys`，一次**合法**重锚就会在 diff 里
    把 10 个实测键全部写成"删+加"，正好抹掉"diff 里没有数值行被改动"这一
    本轮赖以复核的硬约束。`test_anchor_integrity.py` 有测试锁住此性质。
    """
    return json.dumps(anchor, ensure_ascii=False, indent=2, sort_keys=True) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="重锚默认路径逐位锚（无 --confirm 时只打印差异，不写盘）")
    parser.add_argument("--confirm", action="store_true",
                        help="确认写盘（须已人工核对差异，并在位移台账留迹）")
    parser.add_argument("--env", default=None,
                        help="首建锚时必须显式给出：环境指纹，如 "
                             "'Python 3.13.7 / numpy 2.3.3 / scipy 1.16.2'")
    parser.add_argument("--reconcile-source", default=None,
                        help="首建锚时必须显式给出：对账源（仓内 POSIX 相对路径，"
                             "如 docs/superpowers/plans/baseline-2026-09-25.md）")
    args = parser.parse_args(argv)

    # 首建分支（锚文件不存在 ⇒ `old` 为空）：必须显式给出两项指纹，否则会写出
    # `_env="unknown"` + 指针键为空的锚，当场被 test_anchor_integrity.py 的指纹
    # 守卫判红，逼维护者手工改 JSON——正是锚加固要消灭的路径（R34）。
    # `parser.error` = 打印用法 + 非零退出（2），且在任何写盘/算例之前，
    # 故"缺参数 ⇒ 不写盘"是结构性保证，不依赖后续分支的自觉。
    first_build = not ANCHOR.is_file()
    if first_build:
        missing = [opt for opt, value in (("--env", args.env),
                                          ("--reconcile-source", args.reconcile_source))
                   if not (value or "").strip()]
        if missing:
            parser.error(
                f"锚文件不存在（首建：{ANCHOR}）：必须显式提供 {' 与 '.join(missing)}，"
                "否则写出的锚会缺来源指纹（_env='unknown'、_env_reconcile_source 为空），"
                "被 tests/contract/test_anchor_integrity.py 判红；本脚本拒绝写盘。")

    print(f"[锚文件] {ANCHOR}")
    if first_build:
        print("[首建] 锚文件不存在：--env/--reconcile-source 已给出并校验非空，"
              "加 --confirm 才会写盘。")
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

    if first_build:
        env, reconcile_source = args.env, args.reconcile_source
    else:
        env = old.get("_env", "unknown")
        reconcile_source = old.get("_env_reconcile_source", "")
    fingerprint = build_fingerprint(
        env=env,
        generated_from=_new_generated_from(old.get("_generated_from", ""), head),
        reconcile_source=reconcile_source)
    new_anchor = {**fingerprint, **measured}
    ANCHOR.write_text(serialize_anchor(new_anchor), encoding="utf-8")
    print(f"\n[写盘] 已写入 {ANCHOR}（保留并更新指纹键：{sorted(fingerprint)}）")
    print(f"[留迹] 请把以上差异贴进 {LEDGER_HINT}，并与改动同批次提交。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
