# -*- coding: utf-8 -*-
"""锚完整性守卫（R17/R18/R19）：缺锚必须失败而非跳过；键集必须齐备；必须有来源指纹。

本文件是"默认路径零数值影响"证明的元守卫：主锚测试
（`test_default_path_bitwise_anchor.py`）只回答"数字变了没有"，
本文件回答"它是否真的在检查、检查得够不够全、基准本身有没有被换掉"。
"""
import json
import re
import subprocess
from pathlib import Path

from tests.contract.test_default_path_bitwise_anchor import (
    ANCHOR, FINGERPRINT_KEYS, POINTER_KEYS, _FIELD_ATTRS, _SCALAR_FIELDS,
    expected_keys)

# R41（单一真源）：``FINGERPRINT_KEYS`` **不再在本文件重复定义**，改为从主锚测试
# import（那里是唯一定义处）。此前本文件写集合、主锚测试写元组 ⇒ 两处各一份，
# 改一处不会让另一处变红，正是本计划所防的"一事实两处写"漂移。
ANCHOR_TEST_MODULE = "tests.contract.test_default_path_bitwise_anchor"

# 锚里 sha_* 的摘要口径：sha256 = 64 字符小写 hex（R27-1）。
# 守卫存在理由：若将来有人把摘要算法换成 ≤16 字节（如 md5/截断 sha），
# 主测试的"逐位"比对将不再对场数组的微小扰动敏感（抗碰撞/截断），
# 锚会在毫无察觉的情况下失去判别力。
_SHA_LEN = 64
_SHA_RE = re.compile(r"^[0-9a-f]{64}$")


def test_anchor_file_must_exist():
    assert ANCHOR.is_file(), (
        f"逐位锚文件缺失：{ANCHOR} —— 缺锚会让默认路径的零数值影响证明静默失效，"
        "严禁用 skip 掩盖")


def test_anchor_covers_exactly_expected_keys():
    want = json.loads(ANCHOR.read_text(encoding="utf-8"))
    keys = set(want) - set(FINGERPRINT_KEYS) - set(POINTER_KEYS)
    assert keys == expected_keys(), (
        f"锚键集不匹配：缺 {expected_keys() - keys}；多 {keys - expected_keys()}")


def test_expected_keys_derived_not_hardcoded():
    """键集必须由字段表推导，避免两处各写一份而漂移。"""
    assert expected_keys() == {f"sha_{n}" for n in _FIELD_ATTRS} | set(_SCALAR_FIELDS)


def test_scalar_fields_are_measured_not_declared_silently():
    """`_SCALAR_FIELDS` 必须是非空映射，且键与锚里的标量键一一对应。

    Task 10 实测口径为 `{锚键: summary["最终结果"] 下的中文列名}` 字典；
    本断言只在"被改成裸元组/空表"时报警，不绑定具体列名（列名由主测试覆盖）。
    """
    assert isinstance(_SCALAR_FIELDS, dict) and _SCALAR_FIELDS, (
        "_SCALAR_FIELDS 必须是非空映射（锚键 -> summary 列名）")
    want = json.loads(ANCHOR.read_text(encoding="utf-8"))
    assert set(_SCALAR_FIELDS) <= set(want), (
        f"字段表声明了锚里没有的标量键：{set(_SCALAR_FIELDS) - set(want)}")


def test_anchor_carries_provenance_note():
    want = json.loads(ANCHOR.read_text(encoding="utf-8"))
    assert set(FINGERPRINT_KEYS) <= set(want), (
        f"锚必须带来源指纹 {sorted(FINGERPRINT_KEYS)}："
        "微网格非生产数字的警告 + 生成条件 + 环境版本，"
        f"当前缺 {sorted(set(FINGERPRINT_KEYS) - set(want))}")
    for key in sorted(FINGERPRINT_KEYS):
        assert isinstance(want[key], str) and want[key].strip(), (
            f"锚指纹 {key} 必须是非空字符串")


def test_sha_values_are_full_sha256_hex():
    """R27-1：每个 `sha_*` 必须是 64 字符全小写 hex（sha256）。

    防的是"摘要算法被换成 ≤16 字节 / 被截断"这类**静默失牙**：
    它不会让任何测试变红（改算法的同批次自然会把锚一起重写），
    却会让主测试对场数组的微扰不再敏感——锚看起来还在检查，实际已不咬人。
    """
    want = json.loads(ANCHOR.read_text(encoding="utf-8"))
    sha_keys = sorted(k for k in want if k.startswith("sha_"))
    assert sha_keys, "锚里没有任何 sha_* 键（字段表与锚脱节）"
    for key in sha_keys:
        value = want[key]
        assert isinstance(value, str), f"{key} 必须是字符串，实为 {type(value).__name__}"
        assert len(value) == _SHA_LEN, (
            f"{key} 长度应为 {_SHA_LEN}（sha256），实为 {len(value)}——"
            "摘要算法疑似被截断/换短，锚会静默失去微扰判别力")
        assert _SHA_RE.match(value), (
            f"{key} 含非小写 hex 字符（大小写混用或非 hex）：{value[:16]}…——"
            "口径必须是 sha256 的小写十六进制字符串")


def test_anchor_run_does_not_skip_when_anchor_missing(monkeypatch):
    """R17 行为守卫：锚文件缺失时，主测试必须 **FAIL**，绝不能 SKIP。

    不用"读源码找 pytest.skip 字样"（脆：注释里出现该词即误判），
    改为把 `ANCHOR` 指向一个不存在路径后真跑一次主测试函数：
    若它抛 `Skipped` ⇒ 假绿路径仍在；若它抛 `AssertionError`/`FileNotFoundError`
    ⇒ 缺锚即硬失败，R17 意图成立。

    **失效前提（R27-2，务必知晓）**：本守卫依赖主测试以**模块全局 `ANCHOR`**
    读取路径（monkeypatch 改的就是这个全局）。若主测试将来被重构为函数内
    `Path(__file__).with_name("_default_path_anchor_hu101.json")` 之类内联写法，
    monkeypatch 将不再影响被测函数——**本守卫会静默失去判别力**（仍 PASS，
    但什么都没验证）。彼时须改回源码字符串检查，或同步调整主测试保留全局常量。
    """
    missing = ANCHOR.with_name("__不存在的锚文件__task11__.json")
    assert not missing.exists(), f"占位路径意外存在：{missing}"
    mod = __import__(ANCHOR_TEST_MODULE, fromlist=["*"])
    monkeypatch.setattr(mod, "ANCHOR", missing)
    try:
        mod.test_default_path_matches_bitwise_anchor()
    except Exception as exc:  # noqa: BLE001 —— 需要按类型分流
        assert not isinstance(exc, __import__("pytest").skip.Exception), (
            "主测试在缺锚时走了 pytest.skip（R17 假绿路径仍在）：缺锚必须直接 FAIL")
    else:
        raise AssertionError(
            "主测试在锚文件缺失时未报错——说明它既没检查锚也没比对，属假绿")


def test_fingerprint_helper_records_env_source():
    """R27-3：对账源必须写在**可独立解析**的指针键里，且该文件真实存在。

    原先对账源被全角括号包在 `_env` 那句自然语言中间（"…（对账源：path）"），
    只能靠子串包含判断，且无法被脚本/工具当作路径直接解析。现改为独立的
    `_env_reconcile_source` 键：值本身就是一个仓内相对路径（纯 POSIX 相对路径）。
    """
    want = json.loads(ANCHOR.read_text(encoding="utf-8"))
    for key in POINTER_KEYS:
        assert key in want, f"锚缺少可解析指针键 {key}（R27-3）"
        pointer = want[key]
        assert isinstance(pointer, str) and pointer.strip(), (
            f"{key} 必须是非空字符串，实为 {pointer!r}")
        # 必须是**可独立解析**的相对路径：不夹在句子里，可被纯逻辑还原成 Path。
        assert pointer == pointer.strip(), f"{key} 首尾有空白，无法直接解析：{pointer!r}"
        assert not pointer.startswith("/") and "\\" not in pointer, (
            f"{key} 必须是 POSIX 相对路径（正斜杠、非绝对路径）：{pointer!r}")
        assert re.match(r"^[\w./\-]+\.md$", pointer), (
            f"{key} 不是可解析的 markdown 相对路径：{pointer!r}")
        root = Path(__file__).resolve().parents[2]
        assert (root / pointer).is_file(), f"对账源文件不存在：{root / pointer}"
    # `_env` 仍需保留环境版本信息（指纹可读性），但不再承载路径。
    assert "Python" in want["_env"], f"_env 应保留环境版本信息：{want['_env']!r}"


def test_anchor_values_are_committed_not_just_working_tree():
    """R19 残余：锚若只存在于工作树（未入库），CI/他人 clone 仍会缺锚。

    用 `git ls-files` 确认锚文件已被版本控制跟踪。

    无 git 的环境（镜像/沙箱/只 checkout 子路径的 CI 容器）必须**显式失败**，
    不能用 `pytest.skip`：skip 会重开 R17 的假绿路径（"验不了就变绿"）。
    报错文案必须能区分"环境无 git"与"锚未入库"两种原因，否则读者会把
    `FileNotFoundError` 误判为环境坏了、而真实原因（锚未入库）被掩盖。
    """
    root = Path(__file__).resolve().parents[2]
    rel = ANCHOR.relative_to(root).as_posix()
    try:
        proc = subprocess.run(
            ["git", "ls-files", "--error-unmatch", rel],
            cwd=root, capture_output=True, text=True)
    except FileNotFoundError:
        assert False, (
            f"环境无 git，无法验证锚已入库：{rel}——"
            "请在有 git 的环境运行，或人工确认锚已在版本控制内。"
            "（本守卫不 skip：skip 会让'验不了'伪装成'验过了'，即 R17 假绿路径）"
        )
    assert proc.returncode == 0, (
        f"锚文件未被 git 跟踪（尚未 add/commit）：{rel}\n{proc.stderr.strip()}")


def test_committed_guard_fails_loudly_when_git_absent(monkeypatch):
    """R26 行为守卫：环境无 git 时，入库守卫必须**显式 FAIL**，不得 skip/xfail。

    真实场景：镜像/沙箱/只 checkout 子路径的 CI 容器 PATH 无 `git`。若不接住
    `FileNotFoundError`，pytest 会把该用例报成 **FAILED**（异常点就在用例体内），
    只是 traceback 只剩 `FileNotFoundError`，与真实原因（锚未入库）毫无关系，
    读者会误判为环境坏了；若改用 skip，则是重开 R17 的假绿路径。
    这里把 `subprocess.run` 替成必抛 `FileNotFoundError` 的桩，断言结果是一条
    **区分性 AssertionError**。

    **禁止用 skip/xfail 兜底（I-2 子性质）**：`Skipped` / `XFailed` 是
    `BaseException` 子类，只 `except AssertionError` 会让它们**穿透本用例**，
    把**守卫自己**报成 SKIPPED/XFAIL（可见但非失败）——那正好重开假绿路径。
    故下面对 `AssertionError` 之外的一切异常一律判 FAIL，并把文案写死。
    """
    import subprocess as _sp

    def _no_git(*args, **kwargs):
        raise FileNotFoundError("[WinError 2] 系统找不到指定的文件")

    monkeypatch.setattr(_sp, "run", _no_git)
    import tests.contract.test_anchor_integrity as _self
    monkeypatch.setattr(_self.subprocess, "run", _no_git, raising=False)
    try:
        _self.test_anchor_values_are_committed_not_just_working_tree()
    except AssertionError as exc:
        msg = str(exc)
        assert "环境无 git" in msg and "无法验证锚已入库" in msg, (
            f"缺 git 的失败文案必须能区分两种原因（环境 vs 未入库），实为：{msg}")
    except BaseException as exc:  # noqa: BLE001 —— skip/xfail 也要判为失败
        raise AssertionError(
            f"环境无 git 时入库守卫抛出 {type(exc).__name__}：{exc}——"
            "禁止用 skip/xfail 兜底：那会把'验不了'伪装成'验过了'（R17 假绿路径），"
            "并会把本守卫自己报成 SKIPPED/XFAIL 而非 FAILED") from exc
    else:
        raise AssertionError(
            "环境无 git 时入库守卫未失败——它要么静默通过（假绿），"
            "要么异常没被接住而只留下 FileNotFoundError 的 traceback"
            "（读者会误判成环境坏了，真实原因'锚未入库'被掩盖）")


def _measured_key_lines(text: str) -> list[str]:
    """抽出锚文本里 10 个**实测键**（`expected_keys()`）的整行，保持出现顺序。"""
    keys = expected_keys()
    return [line for line in text.splitlines()
            if line.split(":", 1)[0].strip().strip('"') in keys]


def test_regenerate_script_preserves_measured_key_lines(tmp_path):
    """R33：一次**合法**重锚不得改动/重排 10 个实测键的文本行。

    重锚脚本写盘唯一走 `serialize_anchor()`。`{**fingerprint, **measured}` 的插入序
    是"指纹键在前、实测键按 `_run_default_case()` 产出序在后"，与锚文件里 `sorted()`
    的键序**不同**；若序列化漏了 `sort_keys=True`，Task 8 真重锚时 diff 会把 10 个
    实测键全部写成"删+加"——"diff 里没有数值行被改动"这一本轮赖以复核的硬约束
    当场崩塌。故此处用同一份数据复现写盘路径，逐行比对重写前后。
    """
    import tests.contract._regenerate_default_path_anchor as regen

    original = ANCHOR.read_text(encoding="utf-8")
    data = json.loads(original)

    # 按脚本写盘时的真实构造顺序复现：fingerprint 先入、measured 后入；
    # measured 的顺序与 `_run_default_case()` 一致（先五场 sha_*，后五标量）。
    fingerprint = regen.build_fingerprint(
        env=data["_env"],
        generated_from=data["_generated_from"],
        reconcile_source=data["_env_reconcile_source"])
    measured = {f"sha_{name}": data[f"sha_{name}"] for name in _FIELD_ATTRS}
    measured.update({key: data[key] for key in _SCALAR_FIELDS})

    written = tmp_path / ANCHOR.name
    written.write_text(regen.serialize_anchor({**fingerprint, **measured}),
                       encoding="utf-8")
    rewritten = written.read_text(encoding="utf-8")

    before, after = _measured_key_lines(original), _measured_key_lines(rewritten)
    assert len(before) == len(expected_keys()), (
        f"锚里实测键行数应为 {len(expected_keys())}，实为 {len(before)}（字段表与锚脱节）")
    assert after == before, (
        "重锚写盘改变了 10 个实测键的文本行（内容或顺序，R33）：\n"
        f"  重写前：{before}\n  重写后：{after}\n"
        "serialize_anchor() 必须以 sort_keys=True 序列化，否则合法重锚会在 diff 里"
        "把每个数值行写成『删+加』，抹掉『数值行未被改动』这一可核对形式")


def test_regenerate_refuses_first_build_without_fingerprint(monkeypatch, tmp_path):
    """R34：首建锚缺 `--env`/`--reconcile-source` 时应**非零退出且不写盘**。

    锚不存在 ⇒ `old = {}` ⇒ 旧代码会写 `_env="unknown"`、指针键为空的锚，
    当场被本文件的指纹守卫判红，逼维护者手工改 JSON——正是锚加固要消灭的路径。
    正向（给全两参数）也必须仍然可写，否则"拒绝一切"式修复会假绿通过。
    """
    import pytest as _pytest

    import tests.contract._regenerate_default_path_anchor as regen

    anchor = tmp_path / "_default_path_anchor_hu101.json"
    monkeypatch.setattr(regen, "ANCHOR", anchor)
    # 用桩替换 10s 实测端到端算例：本条只测"首建参数校验"，不该跑模型；
    # 也让"守卫若失效"的失败保持廉价（否则负向路径会白跑 4 次 10s 算例）。
    monkeypatch.setattr(
        regen, "_run_default_case",
        lambda: {**{f"sha_{name}": "0" * 64 for name in _FIELD_ATTRS},
                 **{key: 1.0 for key in _SCALAR_FIELDS}})

    def _run(argv: list[str]) -> int:
        """跑一次脚本入口并返回退出码（正常返回 0；`parser.error` ⇒ SystemExit(2)）。"""
        with _pytest.raises(SystemExit) as excinfo:
            regen.main(argv)
        code = excinfo.value.code
        assert isinstance(code, int), f"退出码应为 int，实为 {code!r}"
        return code

    # 负向：三种缺法都必须在写盘前被拦下（usage 错误 ⇒ 码 2）。
    assert _run([]) == 2
    assert _run(["--confirm"]) == 2
    assert _run(["--confirm", "--env", "Python 3.13.7"]) == 2
    assert _run(["--confirm", "--reconcile-source",
                 "docs/superpowers/plans/baseline-2026-09-25.md"]) == 2
    assert not anchor.exists(), "首建缺指纹时脚本仍写了盘"

    # 正向：给全两参数后必须能写出结构合法的锚。
    # 负向靠 `parser.error` ⇒ `SystemExit(2)`；正向是 `main()` 正常 **return 0**，
    # 两者不是同一种结束方式，故正向直接调 `main()`。
    assert regen.main(["--confirm", "--env", "Python 3.13.7 / numpy 2.3.3",
                       "--reconcile-source",
                       "docs/superpowers/plans/baseline-2026-09-25.md"]) == 0
    written = json.loads(anchor.read_text(encoding="utf-8"))
    assert set(written) == expected_keys() | set(POINTER_KEYS) | set(FINGERPRINT_KEYS), (
        f"首建锚键集不对：{sorted(set(written))}")
    assert written["_env"] == "Python 3.13.7 / numpy 2.3.3", (
        f"首建锚必须用 --env 的值，而非 'unknown'：{written['_env']!r}")
    assert written["_env_reconcile_source"] == (
        "docs/superpowers/plans/baseline-2026-09-25.md"), (
        f"首建锚必须用 --reconcile-source 的值，而非空串："
        f"{written['_env_reconcile_source']!r}")
    assert written["_generated_from"].startswith("HEAD "), written["_generated_from"]
