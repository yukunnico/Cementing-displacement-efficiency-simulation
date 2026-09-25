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
    ANCHOR, POINTER_KEYS, _FIELD_ATTRS, _SCALAR_FIELDS, expected_keys)

FINGERPRINT_KEYS = {"_note", "_generated_from", "_env"}
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
    keys = set(want) - FINGERPRINT_KEYS - set(POINTER_KEYS)
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
    assert FINGERPRINT_KEYS <= set(want), (
        f"锚必须带来源指纹 {sorted(FINGERPRINT_KEYS)}："
        "微网格非生产数字的警告 + 生成条件 + 环境版本，"
        f"当前缺 {sorted(FINGERPRINT_KEYS - set(want))}")
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
    """R26 行为守卫：环境无 git 时，入库守卫必须**显式 FAIL**，不得 skip/ERROR。

    真实场景：镜像/沙箱/只 checkout 子路径的 CI 容器 PATH 无 `git`。若不接住
    `FileNotFoundError`，pytest 会报 **ERROR**，文本与真实原因（锚未入库）毫无
    关系，读者会误判为环境坏了；若改用 skip，则是重开 R17 的假绿路径。
    这里把 `subprocess.run` 替成必抛 `FileNotFoundError` 的桩，断言结果是一条
    **区分性 AssertionError**。
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
    else:
        raise AssertionError(
            "环境无 git 时入库守卫未失败——它要么静默通过（假绿），"
            "要么抛出非 AssertionError（会被报成 ERROR 而掩盖真实原因）")
