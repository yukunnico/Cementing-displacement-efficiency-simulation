# -*- coding: utf-8 -*-
"""锚完整性守卫（R17/R18/R19）：缺锚必须失败而非跳过；键集必须齐备；必须有来源指纹。

本文件是"默认路径零数值影响"证明的元守卫：主锚测试
（`test_default_path_bitwise_anchor.py`）只回答"数字变了没有"，
本文件回答"它是否真的在检查、检查得够不够全、基准本身有没有被换掉"。
"""
import json
import subprocess
from pathlib import Path

from tests.contract.test_default_path_bitwise_anchor import (
    ANCHOR, _FIELD_ATTRS, _SCALAR_FIELDS, expected_keys)

FINGERPRINT_KEYS = {"_note", "_generated_from", "_env"}
ANCHOR_TEST_MODULE = "tests.contract.test_default_path_bitwise_anchor"


def test_anchor_file_must_exist():
    assert ANCHOR.is_file(), (
        f"逐位锚文件缺失：{ANCHOR} —— 缺锚会让默认路径的零数值影响证明静默失效，"
        "严禁用 skip 掩盖")


def test_anchor_covers_exactly_expected_keys():
    want = json.loads(ANCHOR.read_text(encoding="utf-8"))
    keys = set(want) - FINGERPRINT_KEYS
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


def test_anchor_run_does_not_skip_when_anchor_missing(monkeypatch):
    """R17 行为守卫：锚文件缺失时，主测试必须 **FAIL**，绝不能 SKIP。

    不用"读源码找 pytest.skip 字样"（脆：注释里出现该词即误判），
    改为把 `ANCHOR` 指向一个不存在路径后真跑一次主测试函数：
    若它抛 `Skipped` ⇒ 假绿路径仍在；若它抛 `AssertionError`/`FileNotFoundError`
    ⇒ 缺锚即硬失败，R17 意图成立。
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
    """`_env` 必须指向对账源文件，且该文件真实存在（防指纹变成一句空话）。"""
    want = json.loads(ANCHOR.read_text(encoding="utf-8"))
    note = want["_env"]
    source = "docs/superpowers/plans/baseline-2026-09-25.md"
    assert source.replace("/", "\\") in note or source in note, (
        f"_env 未指向对账源 {source}：{note!r}")
    root = Path(__file__).resolve().parents[2]
    assert (root / source).is_file(), f"对账源文件不存在：{root / source}"


def test_anchor_values_are_committed_not_just_working_tree():
    """R19 残余：锚若只存在于工作树（未入库），CI/他人 clone 仍会缺锚。

    用 `git ls-files` 确认锚文件已被版本控制跟踪。
    """
    root = Path(__file__).resolve().parents[2]
    rel = ANCHOR.relative_to(root).as_posix()
    proc = subprocess.run(
        ["git", "ls-files", "--error-unmatch", rel],
        cwd=root, capture_output=True, text=True)
    assert proc.returncode == 0, (
        f"锚文件未被 git 跟踪（尚未 add/commit）：{rel}\n{proc.stderr.strip()}")
