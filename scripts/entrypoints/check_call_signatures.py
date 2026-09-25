#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""扫描全仓对 `AnnulusD2DGASolver` / `CasingFlowSolver` 的调用关键字，比对真实签名。

背景（Task 5）：形参一旦删除（如 `enable_d2dga_auto_m`、`dispersion_dt_scale`），
仍照旧传参的调用点只在**运行期**抛 `TypeError` —— 静态无人察觉，直到有人真的跑
那个脚本。本扫描器把这条腐烂路径提前到测试期。

Task 14 扩展（闸门"绿而无用"的修复）
------------------------------------
原闸门只认**直接关键字调用**，于是 `dict(...)` + `**kw` 展开形态（含 dict 字面量在
**另一个文件**里构造、经函数形参中转的情形）全成盲区 ⇒ 闸门报绿而脚本运行即崩。
现增加**仓储级 `dict 键 → ** splat → 目标构造函数`** 追踪，覆盖：

  * 同文件同作用域的绑定：`name = {...}` / `name = dict(...)` /
    `name.update(...)` / `name = a | {...}` / `name[key] = ...`（后者视为"任意键"）；
  * 跨文件形参中转：`f(..., extra_kw={...})` 与 `f(..., **compute)` 递归解析；
  * `for k, v in MODULE_LEVEL_DICT.items():` —— `v` 取该字典**全部值**的并集。

解析是**并集**近似（同名形参的多调用点合并、字典值并集、`_MAX_DEPTH` 截断），
因此**只可能多报**（且多报也只报"在某个真实调用路径上确定会出现"的坏键），
绝不会因为"解析不到"而漏掉一个确定会出现的坏键。

用法::

    python scripts/entrypoints/check_call_signatures.py    # 有坏调用 ⇒ 退出码 1

规则：
  1. 只认 AST 里 `func` 名为两个目标类（`Name` 或 `Attribute`，故 `m.AnnulusD2DGASolver`
     与 `AnnulusD2DGASolver` 都算）的 `Call` 节点；
  2. 关键字的 `arg` / 展开 dict 的键不在该类 `__init__` 的真实形参集合内 ⇒ 报错；
  3. **唯一豁免**须三条**同时**成立 —— 调用位于 `pytest.raises(...)` 上下文内 ∧ 该
     `raises` **解析到 pytest**（`import pytest` / `import pytest as pt` / `from pytest
     import raises`，即接收者有 import 证据）∧ 其**首个位置参数是 `TypeError`** ∧
     文件位于 `tests/` 目录下。那是"传值必须显式失败"的故意反例
     （`tests/contract/test_no_invented_dispersion.py`）。三条缺一即算坏调用：
     自定义的 `def raises(e)` 上下文管理器（无论文件在不在 `tests/` 下）、
     `helpers.raises(...)`、`raises(ValueError)` 都**不能**蹭豁免。
     豁免数会在输出中显式打印，不做静默吞掉。
  4. 盲区显式化：`SyntaxError` / `UnicodeDecodeError` 的文件**打印跳过计数**（不静默）；
     `with` 与 `async with` 一视同仁；目标类日后若新增 `**kwargs`（`VAR_KEYWORD`），
     该类关键字校验**显式停用并告警**，而不是把全部关键字误报成"形参不存在"。

开销（本仓 279 个可解析文件、0 跳过）：
  全仓扫描是**秒级**操作，支配项是 `Repo.load`（`ast.parse` + 逐节点绑定事实收集），
  不是关键字比对本身；纯 `ast.parse` 全仓约 0.8 s，`Repo.load` 约 2.6–2.9 s，
  `_scan` 合计约 **3.6–4.7 s**（两次独立测量：3.62–4.69 s / 3.77–4.10 s）。
  ⚠️ **墙钟随机器负载浮动**（同一提交态在不同负载下可差近 2 倍），因此请预期一个**区间**
  而不是某个固定常数；`visible_calls` / `_signatures` / `_var_keyword_targets` 已做进程内
  缓存（缓存只影响耗时，实测不改变判定结果）。本闸门挂在合并前路径上
  （`tests/contract/test_entrypoint_signatures.py::test_no_bad_keyword_arguments_anywhere`），
  改动前请把这一量级计入预算。
"""
from __future__ import annotations

import ast
import inspect
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import NamedTuple

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

# 跨文件形参中转的解析深度上限（防环、防病态链式调用把扫描拖死）
_MAX_DEPTH = 8

# 目标类签名 / `**kwargs` 检出结果缓存（每次全仓扫描都会调；re-inspect 无意义开销）
_SIGNATURE_CACHE: dict[str, set[str]] | None = None
_VAR_KEYWORD_CACHE: set[str] | None = None


# --------------------------------------------------------------------------- #
# 仓储级 AST 事实
# --------------------------------------------------------------------------- #

@dataclass
class _Scope:
    """一个词法作用域（模块体 / 类体 / 函数体）内与字典相关的绑定事实。"""

    key: tuple[str, ...] = ()
    kind: str = "module"                      # module | class | function
    func_name: str = ""                       # 函数作用域的简单名（用于调用点索引）
    params: tuple[str, ...] = ()
    n_positional: int = 0                     # 位置形参个数（含 posonly）
    kwarg: str = ""                           # `**name` 形参名（空 = 无）
    # name -> 赋值右值表达式（多次赋值全部保留，取并集）
    assigns: dict[str, list[ast.expr]] = field(default_factory=dict)
    # name -> [(显式新增键, 展开表达式列表)]，来自 name.update(...)
    merges: dict[str, list[tuple[dict[str, ast.expr], list[ast.expr]]]] = field(
        default_factory=dict)
    # name -> 该名字可能含任意键（下标赋值 / 增强赋值）
    open_defs: set[str] = field(default_factory=set)
    # name -> (槽位, 迭代表达式)：for 循环目标名；槽位 -1=整体, i=元组第 i 项
    foreach: dict[str, tuple[int, ast.expr]] = field(default_factory=dict)
    # name -> 列表字面量元素（含 name.append(...) / name.extend([...]) 累加）
    list_elems: dict[str, list[ast.expr]] = field(default_factory=dict)


@dataclass
class _File:
    """单文件解析结果。"""

    path: Path
    tree: ast.AST
    scopes: dict[tuple[str, ...], _Scope] = field(default_factory=dict)
    owners: dict[tuple[str, ...], ast.AST] = field(default_factory=dict)
    scope_of: dict[int, tuple[str, ...]] = field(default_factory=dict)


class _Hit(NamedTuple):
    """一处坏调用点。`origin` 为 `direct` 或 `splat`；`via_line` 为 splat 所在行。"""

    path: str
    lineno: int
    cls: str
    key: str
    origin: str
    via_line: int


@dataclass
class _ScanResult:
    bad: list[_Hit] = field(default_factory=list)
    exempt: list[_Hit] = field(default_factory=list)
    skipped: list[tuple[str, str]] = field(default_factory=list)
    disabled_targets: list[str] = field(default_factory=list)


def _build_file(path: Path, tree: ast.AST) -> _File:
    """把一棵 AST 拆成作用域 + 绑定事实。"""
    facts = _File(path=path, tree=tree)
    facts.scopes[()] = _Scope(key=(), kind="module")

    def descend(node: ast.AST, scope_key: tuple[str, ...]) -> None:
        facts.scope_of[id(node)] = scope_key
        for child in ast.iter_child_nodes(node):
            if not isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef,
                                      ast.ClassDef)):
                descend(child, scope_key)
                continue
            key = scope_key + (child.name,)
            if key in facts.scopes:                      # 同名嵌套定义：加后缀避免覆盖
                suffix = 2
                while key in facts.scopes:
                    key = scope_key + (f"{child.name}#{suffix}",)
                    suffix += 1
            if isinstance(child, ast.ClassDef):
                facts.scopes[key] = _Scope(key=key, kind="class", func_name=child.name)
            else:
                args = child.args
                names = [a.arg for a in (*args.posonlyargs, *args.args, *args.kwonlyargs)]
                if args.vararg:
                    names.append(args.vararg.arg)
                kwarg = args.kwarg.arg if args.kwarg else ""
                if kwarg:
                    names.append(kwarg)
                facts.scopes[key] = _Scope(key=key, kind="function",
                                           func_name=child.name, params=tuple(names),
                                           n_positional=len(args.posonlyargs) + len(args.args),
                                           kwarg=kwarg)
            facts.owners[key] = child
            descend(child, key)

    descend(tree, ())
    _collect_bindings(facts)
    return facts


def _collect_bindings(facts: _File) -> None:
    """逐节点登记字典绑定 / 展开 / for-别名。"""
    for node in ast.walk(facts.tree):
        scope = facts.scopes[facts.scope_of[id(node)]]
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            value = node.value
            for target in targets:
                if isinstance(target, ast.Name) and value is not None:
                    scope.assigns.setdefault(target.id, []).append(value)
                    if isinstance(value, (ast.List, ast.Tuple)):
                        scope.list_elems.setdefault(target.id, []).extend(value.elts)
                elif isinstance(target, ast.Subscript) and isinstance(target.value, ast.Name):
                    # `name["k"] = v` 是**已知键**的新增；键名不是字符串字面量才记盲区
                    if (isinstance(target.slice, ast.Constant)
                            and isinstance(target.slice.value, str) and value is not None):
                        scope.merges.setdefault(target.value.id, []).append(
                            ({target.slice.value: value}, []))
                    else:
                        scope.open_defs.add(target.value.id)
        elif isinstance(node, ast.AugAssign):
            if isinstance(node.target, ast.Subscript) and isinstance(node.target.value, ast.Name):
                scope.open_defs.add(node.target.value.id)
        elif isinstance(node, ast.Call):
            func = node.func
            if isinstance(func, ast.Attribute) and isinstance(func.value, ast.Name):
                if func.attr == "update":
                    added: dict[str, ast.expr] = {}
                    splats: list[ast.expr] = list(node.args)
                    for kw in node.keywords:
                        if kw.arg is None:
                            splats.append(kw.value)
                        else:
                            added[kw.arg] = kw.value
                    scope.merges.setdefault(func.value.id, []).append((added, splats))
                elif func.attr == "append" and node.args:
                    # `.append(x)` 把 x 作为**一个**元素加进去（不要摊平元组）
                    for arg in node.args:
                        scope.list_elems.setdefault(func.value.id, []).append(arg)
                elif func.attr == "extend" and node.args:
                    for arg in node.args:
                        if isinstance(arg, (ast.List, ast.Tuple)):
                            scope.list_elems.setdefault(func.value.id, []).extend(arg.elts)
                        else:
                            scope.open_defs.add(func.value.id)
        elif isinstance(node, (ast.For, ast.AsyncFor)):
            target = node.target
            if isinstance(target, ast.Name):
                scope.foreach[target.id] = (-1, node.iter)
            elif isinstance(target, (ast.Tuple, ast.List)) and target.elts:
                if len(target.elts) == 1:
                    if isinstance(target.elts[0], ast.Name):
                        scope.foreach[target.elts[0].id] = (-1, node.iter)
                else:
                    for index, element in enumerate(target.elts):
                        if isinstance(element, ast.Name):
                            scope.foreach[element.id] = (index, node.iter)


# --------------------------------------------------------------------------- #
# 仓储级解析
# --------------------------------------------------------------------------- #

def _module_name(path: Path, root: Path) -> str:
    """文件 → 点分模块名（`__init__.py` 退成包名）。"""
    try:
        rel = Path(path).relative_to(root)
    except ValueError:
        rel = Path(path)
    parts = list(rel.with_suffix("").parts)
    if parts and parts[-1] == "__init__":
        parts.pop()
    return ".".join(parts)


def _mod_match(callee: str, imported: str) -> bool:
    """模块名比对：允许 `scripts.foo` 与 `foo` 这类 sys.path 前缀差异。"""
    if not callee or not imported:
        return False
    return callee == imported or callee.endswith("." + imported) or imported.endswith("." + callee)


def _collect_imports(tree: ast.AST) -> dict[str, tuple[str, str | None]]:
    """局部名 → (模块名, 原名)；`import M` 的原名为 None（表示模块对象）。"""
    out: dict[str, tuple[str, str | None]] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            if node.module is None or node.level:
                continue                                # 相对导入不在本仓形态内
            for alias in node.names:
                out[alias.asname or alias.name] = (node.module, alias.name)
        elif isinstance(node, ast.Import):
            for alias in node.names:
                out[alias.asname or alias.name.split(".")[0]] = (alias.name, None)
    return out


class _Repo:
    """全仓事实 + 惰性解析（含跨文件形参中转）。

    跨文件跳转**只在有 import 证据时**发生：同名 `run_case` 散布在十几个探针里，
    若不分模块合并调用点，A 文件的坏键会被算到 B 文件头上（假报）。因此调用点按
    "定义处可见性"索引：同文件局部名，或 `from M import f` / `import M` + `M.f`。
    """

    def __init__(self) -> None:
        self.files: list[_File] = []
        self.skipped: list[tuple[str, str]] = []
        self.module_of: dict[str, str] = {}
        self.file_by_module: dict[str, _File] = {}
        self.imports: dict[str, dict[str, tuple[str, str | None]]] = {}
        self.calls: dict[tuple[str, ...], list[tuple[_File, tuple[str, ...], ast.Call]]] = {}
        self.calls_by_name: dict[str, list[tuple[str, ...]]] = {}
        self.visible_cache: dict[tuple[str, str], list[tuple[_File, tuple[str, ...], ast.Call]]] = {}
        self.memo: dict[tuple, tuple[dict[str, ast.expr | None], bool]] = {}
        self.busy: set[tuple] = set()

    # ---------------------------------------------------------------- 装载
    def load(self, root: Path) -> None:
        root = Path(root)
        for path in sorted(root.rglob("*.py")):
            if any(part in SKIP_DIRS for part in path.parts):
                continue
            try:
                with warnings.catch_warnings():
                    # 仓库里若干 legacy 文件含非法转义序列，ast.parse 会刷
                    # SyntaxWarning（与本闸门职责无关，且会污染 pytest 输出），静音。
                    warnings.simplefilter("ignore", SyntaxWarning)
                    # utf-8-sig：仓库里有带 BOM 的脚本（如 scripts/probes/
                    # ht1_004_sensitivity.py），用 utf-8 读会把 BOM 留成 ﻿ ⇒ 假
                    # SyntaxError ⇒ 整份文件被静默跳过。BOM 兼容后仍留跳过计数兜底。
                    tree = ast.parse(path.read_text(encoding="utf-8-sig"))
            except (SyntaxError, UnicodeDecodeError) as exc:
                self.skipped.append((str(path), type(exc).__name__))
                continue
            self.files.append(_build_file(path, tree))
        for facts in self.files:
            self.module_of[str(facts.path)] = _module_name(facts.path, root)
            self.file_by_module.setdefault(self.module_of[str(facts.path)], facts)
            self.imports[str(facts.path)] = _collect_imports(facts.tree)
        for facts in self.files:
            for node in ast.walk(facts.tree):
                if not isinstance(node, ast.Call):
                    continue
                key = self._call_key(facts, node)
                if key is not None:
                    self.calls.setdefault(key, []).append(
                        (facts, facts.scope_of[id(node)], node))
        # 形参解析要按"函数名"取调用点；预建 名字 → 键 的倒排索引，
        # 否则每次解析都要线性扫 self.calls.items()（全仓 279 文件量级下这是热点）。
        for key in self.calls:
            self.calls_by_name.setdefault(key[-1], []).append(key)

    def _call_key(self, findings: _File, node: ast.Call) -> tuple[str, ...] | None:
        """调用点在"定义处可见性"坐标系里的键。"""
        func = node.func
        imports = self.imports.get(str(findings.path), {})
        if isinstance(func, ast.Name):
            imported = imports.get(func.id)
            if imported and imported[1] is not None:
                return ("import", imported[0], imported[1])
            return ("local", str(findings.path), func.id)
        if isinstance(func, ast.Attribute) and isinstance(func.value, ast.Name):
            module = imports.get(func.value.id)
            if module is not None and module[1] is None:
                return ("import", module[0], func.attr)
        return None                                     # 接收者不可判定 ⇒ 不索引

    def visible_calls(self, facts: _File, func_name: str
                      ) -> list[tuple[_File, tuple[str, ...], ast.Call]]:
        """能真正解析到 `facts` 里这个函数的调用点（同文件 + 显式 import 进来的）。"""
        cache_key = (str(facts.path), func_name)
        cached = self.visible_cache.get(cache_key)
        if cached is not None:
            return cached
        out = list(self.calls.get(("local", str(facts.path), func_name), ()))
        module = self.module_of.get(str(facts.path), "")
        for key in self.calls_by_name.get(func_name, ()):
            if key[0] == "import" and _mod_match(module, key[1]):
                out.extend(self.calls.get(key, ()))
        self.visible_cache[cache_key] = out
        return out

    # ---------------------------------------------------------------- 解析
    def resolve_expr(self, facts: _File, scope: tuple[str, ...],
                     expr: ast.expr | None, depth: int = 0
                     ) -> tuple[dict[str, ast.expr | None], bool]:
        """返回 (键 -> 取值表达式或 None, 是否可能还有看不见的键)。"""
        if expr is None or depth > _MAX_DEPTH:
            return {}, True
        if isinstance(expr, ast.Name):
            return self.resolve_name(facts, scope, expr.id, depth)
        if isinstance(expr, ast.Dict):
            out: dict[str, ast.expr | None] = {}
            open_ended = False
            for key, value in zip(expr.keys, expr.values):
                if key is None:                             # **展开
                    sub, sub_open = self.resolve_expr(facts, scope, value, depth + 1)
                    out.update(sub)
                    open_ended |= sub_open
                    continue
                names = self.resolve_key_names(facts, scope, key, depth)
                if names is not None:
                    for name in names:
                        out[name] = value
                elif isinstance(key, ast.Constant) and isinstance(key.value, str):
                    out[key.value] = value
                else:
                    open_ended = True
            return out, open_ended
        if isinstance(expr, ast.BinOp) and isinstance(expr.op, ast.BitOr):
            left, left_open = self.resolve_expr(facts, scope, expr.left, depth + 1)
            right, right_open = self.resolve_expr(facts, scope, expr.right, depth + 1)
            return {**left, **right}, left_open or right_open
        if isinstance(expr, ast.IfExp) or isinstance(expr, ast.BoolOp):
            # `a if c else b` / `a or b` / `a and b`：取值是其中之一 ⇒ 并集（近似）
            branches = ([expr.body, expr.orelse] if isinstance(expr, ast.IfExp)
                        else list(expr.values))
            merged: dict[str, ast.expr | None] = {}
            merged_open = False
            for branch in branches:
                sub, sub_open = self.resolve_expr(facts, scope, branch, depth + 1)
                merged.update(sub)
                merged_open |= sub_open
            return merged, merged_open
        if isinstance(expr, ast.Call):
            func = expr.func
            name = _call_name(func)
            if name == "dict":
                out = {}
                open_ended = False
                for arg in expr.args:                       # dict(other_mapping)
                    sub, sub_open = self.resolve_expr(facts, scope, arg, depth + 1)
                    out.update(sub)
                    open_ended |= sub_open
                for kw in expr.keywords:
                    if kw.arg is None:
                        sub, sub_open = self.resolve_expr(facts, scope, kw.value, depth + 1)
                        out.update(sub)
                        open_ended |= sub_open
                    else:
                        out[kw.arg] = kw.value
                return out, open_ended
            if name == "copy" and isinstance(func, ast.Attribute):
                return self.resolve_expr(facts, scope, func.value, depth + 1)
        return {}, True

    def resolve_key_names(self, facts: _File, scope: tuple[str, ...], expr: ast.expr,
                          depth: int) -> set[str] | None:
        """表达式作为**字典键**时的可能字符串名集合。

        仅支持一种形态：`for k, v in SOME_DICT.items(): {k: v, ...}` —— 此时键名
        就是 `SOME_DICT` 的键集。其余返回 `None`（不假定）。
        """
        if isinstance(expr, ast.Constant) and isinstance(expr.value, str):
            return {expr.value}
        if not isinstance(expr, ast.Name):
            return None
        current: tuple[str, ...] | None = scope
        while current is not None:
            facts_scope = facts.scopes[current]
            if expr.id in facts_scope.foreach:
                slot, iter_expr = facts_scope.foreach[expr.id]
                base = None
                if (slot == 0 and isinstance(iter_expr, ast.Call)
                        and isinstance(iter_expr.func, ast.Attribute)
                        and iter_expr.func.attr == "items"
                        and isinstance(iter_expr.func.value, ast.Name)):
                    base = iter_expr.func.value.id
                if base is None:
                    return None
                mapping, _ = self.resolve_name(facts, current, base, depth + 1)
                return set(mapping)
            if (expr.id in facts_scope.assigns or expr.id in facts_scope.merges
                    or expr.id in facts_scope.open_defs or expr.id in facts_scope.params):
                return None
            current = current[:-1] if current else None
        return None

    def iter_elements(self, facts: _File, scope: tuple[str, ...], iter_expr: ast.expr,
                      depth: int) -> tuple[list[ast.expr | None], bool]:
        """`for ... in iter_expr` 每次迭代拿到的元素表达式列表（None = 取值未知）。"""
        if depth > _MAX_DEPTH:
            return [], True
        if isinstance(iter_expr, (ast.Tuple, ast.List)):
            return list(iter_expr.elts), False
        if isinstance(iter_expr, ast.Call) and isinstance(iter_expr.func, ast.Attribute):
            attr = iter_expr.func.attr
            base = iter_expr.func.value
            if (attr in ("items", "values") and isinstance(base, ast.Name)
                    and not iter_expr.args):
                mapping, open_ended = self.resolve_name(facts, scope, base.id, depth + 1)
                if attr == "values":
                    return list(mapping.values()), open_ended
                pairs: list[ast.expr | None] = [
                    ast.Tuple(elts=[ast.Constant(value=key, kind=None),
                                    (value if value is not None
                                     else ast.Constant(value=None, kind=None))],
                              ctx=ast.Load())
                    for key, value in mapping.items()]
                return pairs, open_ended
        if isinstance(iter_expr, ast.Name):
            current: tuple[str, ...] | None = scope
            while current is not None:
                facts_scope = facts.scopes[current]
                if iter_expr.id in facts_scope.list_elems:
                    return list(facts_scope.list_elems[iter_expr.id]), False
                if iter_expr.id in facts_scope.assigns:
                    elements: list[ast.expr] = []
                    for value in facts_scope.assigns[iter_expr.id]:
                        if isinstance(value, (ast.List, ast.Tuple)):
                            elements.extend(value.elts)
                    return elements, not elements
                if (iter_expr.id in facts_scope.merges or iter_expr.id in facts_scope.open_defs
                        or iter_expr.id in facts_scope.params):
                    return [], True
                current = current[:-1] if current else None
        return [], True

    def resolve_name(self, facts: _File, scope: tuple[str, ...], name: str,
                     depth: int = 0) -> tuple[dict[str, ast.expr | None], bool]:
        if depth > _MAX_DEPTH:
            return {}, True
        cache_key = (id(facts), scope, name)
        if cache_key in self.memo:
            return self.memo[cache_key]
        if cache_key in self.busy:                          # 环 ⇒ 不假定
            return {}, True
        self.busy.add(cache_key)
        try:
            result = self._resolve_name_inner(facts, scope, name, depth)
        finally:
            self.busy.discard(cache_key)
        self.memo[cache_key] = result
        return result

    def _resolve_name_inner(self, facts: _File, scope: tuple[str, ...], name: str,
                            depth: int) -> tuple[dict[str, ast.expr | None], bool]:
        out: dict[str, ast.expr | None] = {}
        open_ended = False
        current: tuple[str, ...] | None = scope
        while current is not None:
            facts_scope = facts.scopes[current]
            handled = False
            if name in facts_scope.foreach:
                slot, iter_expr = facts_scope.foreach[name]
                handled = True
                if slot == 0:                           # 键名：交给 resolve_key_names
                    open_ended = True
                else:
                    elements, sub_open = self.iter_elements(
                        facts, current, iter_expr, depth + 1)
                    open_ended |= sub_open
                    for element in elements:
                        if slot > 0:
                            if (isinstance(element, (ast.Tuple, ast.List))
                                    and len(element.elts) > slot):
                                element = element.elts[slot]
                            else:
                                open_ended = True
                                continue
                        if element is None:
                            open_ended = True
                            continue
                        sub, sub_open = self.resolve_expr(
                            facts, current, element, depth + 1)
                        out.update(sub)
                        open_ended |= sub_open
            if (name in facts_scope.assigns or name in facts_scope.merges
                    or name in facts_scope.open_defs):
                handled = True
                open_ended |= name in facts_scope.open_defs
                for value in facts_scope.assigns.get(name, ()):
                    sub, sub_open = self.resolve_expr(facts, current, value, depth + 1)
                    out.update(sub)
                    open_ended |= sub_open
                for added, splats in facts_scope.merges.get(name, ()):
                    out.update(added)
                    for splat in splats:
                        sub, sub_open = self.resolve_expr(facts, current, splat, depth + 1)
                        out.update(sub)
                        open_ended |= sub_open
            if name in facts_scope.params and facts_scope.kind == "function":
                # 形参 + 局部重绑定（如 `kw = dict(kw or {})` 的防御性拷贝）取并集：
                # 只看局部赋值会把调用点传进来的键整段丢掉。
                handled = True
                sub, sub_open = self.resolve_param(facts, current, facts_scope, name, depth)
                out.update(sub)
                open_ended |= sub_open
            if handled:
                return out, open_ended
            current = current[:-1] if current else None
        return self.resolve_imported(facts, name, depth)

    def resolve_imported(self, facts: _File, name: str,
                         depth: int) -> tuple[dict[str, ast.expr | None], bool]:
        """本文件没绑定的名字：若是 `from M import X`，去 M 的模块作用域解析 X。

        `CORRECTED_KW` 这类"全仓唯一口径袋"就是靠这条找到的 —— 否则每个用它的
        脚本都成了看不见的展开点。
        """
        imported = self.imports.get(str(facts.path), {}).get(name)
        if not imported or imported[1] is None:
            return {}, True
        module, original = imported
        target = self.file_by_module.get(module)
        if target is None:
            for mod, candidate in self.file_by_module.items():
                if _mod_match(mod, module):
                    target = candidate
                    break
        if target is None or target is facts:
            return {}, True
        return self.resolve_name(target, (), original, depth + 1)

    def resolve_param(self, facts: _File, scope: tuple[str, ...], facts_scope: _Scope,
                      name: str, depth: int) -> tuple[dict[str, ast.expr | None], bool]:
        """形参取值 = 能解析到本函数的调用点上该实参的并集（含 `**splat` 命中）。

        `def f(**kw)` 的 `kw` 形参特殊：调用点上**所有**展开字典的键都落进它，
        外加"未在形参表里显式声明"的关键字实参（它们成为以实参名为名的普通键）。
        """
        is_kwarg = bool(facts_scope.kwarg) and name == facts_scope.kwarg
        out: dict[str, ast.expr | None] = {}
        open_ended = False
        found = False
        for caller, caller_scope, call in self.visible_calls(facts, facts_scope.func_name):
            for index, arg in enumerate(call.args):          # 位置实参 → 形参按位对应
                if index >= facts_scope.n_positional:
                    break
                if isinstance(arg, ast.Starred):
                    open_ended = True
                    continue
                if facts_scope.params[index] == name:
                    found = True
                    sub, sub_open = self.resolve_expr(caller, caller_scope, arg, depth + 1)
                    out.update(sub)
                    open_ended |= sub_open
            for kw in call.keywords:
                if kw.arg is None:
                    mapping, splat_open = self.resolve_expr(
                        caller, caller_scope, kw.value, depth + 1)
                    if is_kwarg:
                        found = True
                        out.update(mapping)
                        open_ended |= splat_open
                    elif name in mapping:
                        found = True
                        value = mapping[name]
                        if value is None:
                            open_ended = True
                        else:
                            sub, sub_open = self.resolve_expr(
                                caller, caller_scope, value, depth + 1)
                            out.update(sub)
                            open_ended |= sub_open
                    elif splat_open:
                        open_ended = True
                elif is_kwarg and kw.arg not in facts_scope.params:
                    found = True
                    out[kw.arg] = kw.value
                elif kw.arg == name:
                    found = True
                    sub, sub_open = self.resolve_expr(
                        caller, caller_scope, kw.value, depth + 1)
                    out.update(sub)
                    open_ended |= sub_open
        if not found:
            return {}, True                                 # 看不见调用点 ⇒ 不假定
        return out, open_ended


# --------------------------------------------------------------------------- #
# 目标类签名 + 调用点扫描
# --------------------------------------------------------------------------- #

def _signatures() -> dict[str, set[str]]:
    """取两个目标类 `__init__` 的真实形参名（单一真源 = 代码本身，不写死清单）。"""
    global _SIGNATURE_CACHE
    if _SIGNATURE_CACHE is None:
        from cemdisp.models2d import AnnulusD2DGASolver
        from cemdisp.transport1d import CasingFlowSolver

        out: dict[str, set[str]] = {}
        for cls in (AnnulusD2DGASolver, CasingFlowSolver):
            params = inspect.signature(cls.__init__).parameters
            out[cls.__name__] = {n for n in params if n != "self"}
        _SIGNATURE_CACHE = out
    return _SIGNATURE_CACHE


def _var_keyword_targets() -> set[str]:
    """接受 `**kwargs` 的目标类：其关键字校验必须显式停用（否则全量误报）。"""
    global _VAR_KEYWORD_CACHE
    if _VAR_KEYWORD_CACHE is None:
        from cemdisp.models2d import AnnulusD2DGASolver
        from cemdisp.transport1d import CasingFlowSolver

        out: set[str] = set()
        for cls in (AnnulusD2DGASolver, CasingFlowSolver):
            params = inspect.signature(cls.__init__).parameters
            if any(p.kind is inspect.Parameter.VAR_KEYWORD for p in params.values()):
                out.add(cls.__name__)
        _VAR_KEYWORD_CACHE = out
    return _VAR_KEYWORD_CACHE


def _call_name(node: ast.AST) -> str | None:
    """取调用目标名：`f(...)` 与 `obj.f(...)` 都还原成 `f`。"""
    if isinstance(node, ast.Name):
        return node.id
    return getattr(node, "attr", None)


def _is_type_error_name(node: ast.AST) -> bool:
    """`TypeError` / `builtins.TypeError` 都认。"""
    if isinstance(node, ast.Name):
        return node.id == "TypeError"
    return isinstance(node, ast.Attribute) and node.attr == "TypeError"


def _under_tests(path: Path, root: Path) -> bool:
    """文件是否位于被扫描根下的 `tests/` 目录内（豁免的第三条）。"""
    try:
        rel = Path(path).resolve().relative_to(Path(root).resolve())
    except ValueError:
        rel = Path(path)
    return "tests" in rel.parts[:-1]


def _is_pytest_raises(ctx: ast.AST, imports: dict[str, tuple[str, str | None]]) -> bool:
    """该上下文管理器调用是否**解析到** `pytest.raises`（豁免的第一要件）。

    认两种形态，且都要有 import 证据：
      * `pytest.raises(...)` / `pt.raises(...)`（`import pytest` / `import pytest as pt`）
        —— 接收者必须解析到 **pytest 模块本身**；
      * `raises(...)`（`from pytest import raises`）—— 局部名必须来自 pytest。
    其余一律不认：`with raises(...)` 里那个 `raises` 若是本文件自定义的
    `def raises(e)`（或 `helpers.raises` / `self.raises` 等），即使文件在 `tests/` 下
    也**不豁免**。
    """
    if not isinstance(ctx, ast.Call):
        return False
    func = ctx.func
    if isinstance(func, ast.Attribute):
        if func.attr != "raises" or not isinstance(func.value, ast.Name):
            return False
        return imports.get(func.value.id) == ("pytest", None)
    if isinstance(func, ast.Name):
        return imports.get(func.id) == ("pytest", "raises")
    return False


def _find_exempt_raises(node: ast.Call, parents: dict[int, ast.AST], path: Path,
                        root: Path, imports: dict[str, tuple[str, str | None]]) -> bool:
    """三条**同时**成立才算豁免：

    ① 调用位于 `pytest.raises(...)` 上下文内（`with` / `async with`），且该 `raises`
       **解析到 pytest**（含 `import pytest as pt`、`from pytest import raises` 两种别名）；
    ② 该 `raises` 的**首个位置参数**是 `TypeError`；
    ③ 文件位于 `tests/` 目录下。

    缺一即算坏调用：非测试文件里自定义 `def raises(e)` 的上下文管理器只满足"在 with 内"；
    `tests/` 目录下自定义的 `def raises(e)` 同样**不豁免**（接收者解析不到 pytest）。
    """
    if not _under_tests(path, root):
        return False
    current = parents.get(id(node))
    while current is not None:
        if isinstance(current, (ast.With, ast.AsyncWith)):
            for item in current.items:
                ctx = item.context_expr
                if (_is_pytest_raises(ctx, imports)
                        and isinstance(ctx, ast.Call) and ctx.args
                        and _is_type_error_name(ctx.args[0])):
                    return True
        current = parents.get(id(current))
    return False


def _scan(root: Path) -> _ScanResult:
    """仓储级扫描：直接关键字 + dict→splat 双重数据流。"""
    root = Path(root)
    sigs = _signatures()
    # ⚠️ 仅测试 monkeypatch 下可达：两个目标类当前都无 `**kwargs`（`_var_keyword_targets()`
    # 实测恒为空，见 tests/contract/test_entrypoint_signatures.py::
    # test_target_classes_accept_no_var_keyword）⇒ 正常运行下 `disabled` 恒为空集。
    # R117③ 要求"目标类日后新增 **kwargs 时显式停用 + 告警"，故这条分支必须保留。
    disabled = _var_keyword_targets() & TARGETS
    repo = _Repo()
    repo.load(root)

    result = _ScanResult(disabled_targets=sorted(disabled), skipped=list(repo.skipped))
    seen: set[tuple] = set()

    def record(facts: _File, node: ast.Call, cls: str, key: str,
               origin: str, via_line: int, parents: dict[int, ast.AST],
               imports: dict[str, tuple[str, str | None]]) -> None:
        marker = (str(facts.path), node.lineno, cls, key)
        if marker in seen:
            return
        seen.add(marker)
        hit = _Hit(str(facts.path), node.lineno, cls, key, origin, via_line)
        (result.exempt if _find_exempt_raises(node, parents, facts.path, root, imports)
         else result.bad).append(hit)

    for facts in repo.files:
        parents = {id(child): node for node in ast.walk(facts.tree)
                   for child in ast.iter_child_nodes(node)}
        file_imports = repo.imports.get(str(facts.path), {})
        for node in ast.walk(facts.tree):
            if not isinstance(node, ast.Call):
                continue
            cls = _call_name(node.func)
            if cls not in TARGETS or cls in disabled:
                # `cls in disabled` 仅测试 monkeypatch 下可达（见 `_scan` 顶部注释）
                continue
            scope = facts.scope_of[id(node)]
            for kw in node.keywords:
                if kw.arg is None:
                    mapping, _ = repo.resolve_expr(facts, scope, kw.value)
                    for key in sorted(mapping):
                        if key not in sigs[cls]:
                            record(facts, node, cls, key, "splat", kw.value.lineno,
                                   parents, file_imports)
                elif kw.arg not in sigs[cls]:
                    record(facts, node, cls, kw.arg, "direct", node.lineno,
                           parents, file_imports)
    return result


def find_bad_kwargs(root: Path) -> list[tuple[str, int, str, str]]:
    """全仓扫描，返回无法解释的未知形参调用点：(文件, 行号, 类名, 未知形参)。"""
    return [(h.path, h.lineno, h.cls, h.key) for h in _scan(Path(root)).bad]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    result = _scan(root)
    for hit in result.bad:
        try:
            shown = Path(hit.path).relative_to(root)
        except ValueError:
            shown = Path(hit.path)
        how = ("" if hit.origin == "direct"
               else f"（经第 {hit.via_line} 行的 dict/** 展开传入）")
        print(f"{shown}:{hit.lineno}: {hit.cls}({hit.key}=...) 形参不存在{how}")
    if result.exempt:
        print(f"[signatures] 已豁免 {len(result.exempt)} 处 pytest.raises(TypeError) "
              f"故意反例（{', '.join(sorted({f'{Path(h.path).name}:{h.lineno}' for h in result.exempt}))}）")
    if result.skipped:
        kinds = ", ".join(sorted({kind for _, kind in result.skipped}))
        print(f"[signatures] 跳过 {len(result.skipped)} 个无法解析的文件（{kinds}）"
              f"：{', '.join(sorted({Path(p).name for p, _ in result.skipped}))}")
    for name in result.disabled_targets:   # 仅测试 monkeypatch 下非空（见 `_scan` 顶部注释）
        print(f"[signatures] ⚠️ 目标类 {name} 接受 **kwargs ⇒ 该类关键字校验已停用"
              f"（显式处理，避免全量误报）；若系误增请复查构造签名")
    print(f"[signatures] {len(result.bad)} 处问题")
    return 1 if result.bad else 0


if __name__ == "__main__":
    sys.exit(main())
