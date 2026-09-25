#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""扫描全仓对 `AnnulusD2DGASolver` / `CasingFlowSolver` 的调用关键字，比对真实签名。

背景（Task 5）：形参一旦删除（如 `enable_d2dga_auto_m`、`dispersion_dt_scale`），
仍照旧传参的调用点只在**运行期**抛 `TypeError` —— 静态无人察觉，直到有人真的跑
那个脚本。本扫描器把这条腐烂路径提前到测试期。

用法::

    python scripts/entrypoints/check_call_signatures.py    # 有坏调用 ⇒ 退出码 1

规则：
  1. 只认 AST 里 `func` 名为两个目标类（`Name` 或 `Attribute`，故 `m.AnnulusD2DGASolver`
     与 `AnnulusD2DGASolver` 都算）的 `Call` 节点；
  2. 关键字的 `arg` 不在该类 `__init__` 的真实形参集合内 ⇒ 报错；
  3. **唯一豁免**：调用位于 `with pytest.raises(...)` 块内 —— 那是"传值必须显式失败"
     的故意反例（如 `tests/contract/test_no_invented_dispersion.py`），
     不是坏调用点。豁免数会在输出中显式打印，不做静默吞掉。
  4. 已知盲区：形参经 `**extra_kw` / dict 展开传参的调用点无法静态判定。
     该盲区**是实打实有内容的**：至少 9 个脚本仍以 dict 字面量 + `**kw` 展开形式
     传已删形参 `dispersion_dt_scale`
     （bisect_hu103_20260902 / corrected_ref_hu1_hu103_20260902 /
     c_verify_convergence_20260902 / debug_hu1_hu103_eta_zero_20260902 /
     dump_tailwindow_2d_v1 / rerun_stop_fix_20260901 / analyze_distortion_fix /
     dispersion_scale_sensitivity_nz250 / isolate_fix_mechanisms），运行期仍会 TypeError。
     本闸门管不到它们；如需覆盖须另做"dict 字面量 → **splat"数据流分析。
"""
from __future__ import annotations

import ast
import inspect
import sys
import warnings
from pathlib import Path

TARGETS = {"AnnulusD2DGASolver", "CasingFlowSolver"}

# 不扫：版本库内部目录、legacy 原型、本计划工作区（含计划自身的示例代码会误报）
SKIP_DIRS = {
    ".git",
    "__pycache__",
    "archive",
    "hu101model",
    "hu102model",
    ".tmp_research",
    ".superpowers",
}


def _signatures() -> dict[str, set[str]]:
    """取两个目标类 `__init__` 的真实形参名（单一真源 = 代码本身，不写死清单）。"""
    from cemdisp.models2d import AnnulusD2DGASolver
    from cemdisp.transport1d import CasingFlowSolver

    out: dict[str, set[str]] = {}
    for cls in (AnnulusD2DGASolver, CasingFlowSolver):
        params = inspect.signature(cls.__init__).parameters
        out[cls.__name__] = {n for n in params if n != "self"}
    return out


def _call_name(node: ast.AST) -> str | None:
    """取调用目标名：`f(...)` 与 `obj.f(...)` 都还原成 `f`。"""
    if isinstance(node, ast.Name):
        return node.id
    return getattr(node, "attr", None)


def _in_raises_block(node: ast.Call, parents: dict[int, ast.AST]) -> bool:
    """调用是否位于 `with pytest.raises(...)` 块内（故意反例的唯一豁免形态）。"""
    cur = parents.get(id(node))
    while cur is not None:
        if isinstance(cur, ast.With):
            for item in cur.items:
                ctx = item.context_expr
                if isinstance(ctx, ast.Call) and _call_name(ctx.func) == "raises":
                    return True
        cur = parents.get(id(cur))
    return False


def _scan(root: Path) -> tuple[list[tuple[str, int, str, str]],
                               list[tuple[str, int, str, str]]]:
    """返回 (待修的坏调用, 已豁免的 raises 反例)，均为 (文件, 行号, 类名, 未知形参)。"""
    sigs = _signatures()
    bad: list[tuple[str, int, str, str]] = []
    exempt: list[tuple[str, int, str, str]] = []
    for path in sorted(Path(root).rglob("*.py")):
        if any(part in SKIP_DIRS for part in path.parts):
            continue
        try:
            with warnings.catch_warnings():
                # 仓库里若干 legacy 文件含非法转义序列，ast.parse 会刷 SyntaxWarning
                # （与本闸门职责无关，且会污染 pytest 输出），此处静音。
                warnings.simplefilter("ignore", SyntaxWarning)
                tree = ast.parse(path.read_text(encoding="utf-8"))
        except (SyntaxError, UnicodeDecodeError):
            continue  # 语法不合法的文件不是本闸门的职责（另有编译检查）
        parents = {id(child): node for node in ast.walk(tree)
                   for child in ast.iter_child_nodes(node)}
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            name = _call_name(node.func)
            if name not in TARGETS:
                continue
            for kw in node.keywords:
                if kw.arg and kw.arg not in sigs[name]:
                    hit = (str(path), node.lineno, name, kw.arg)
                    (exempt if _in_raises_block(node, parents) else bad).append(hit)
    return bad, exempt


def find_bad_kwargs(root: Path) -> list[tuple[str, int, str, str]]:
    """全仓扫描，返回无法解释的未知形参调用点：(文件, 行号, 类名, 未知形参)。"""
    return _scan(Path(root))[0]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    bad, exempt = _scan(root)
    for f, n, c, k in bad:
        try:
            shown = Path(f).relative_to(root)
        except ValueError:
            shown = Path(f)
        print(f"{shown}:{n}: {c}({k}=...) 形参不存在")
    if exempt:
        print(f"[signatures] 已豁免 {len(exempt)} 处 pytest.raises 故意反例"
              f"（{', '.join(sorted({f'{Path(f).name}:{n}' for f, n, _, _ in exempt}))}）")
    print(f"[signatures] {len(bad)} 处问题")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
