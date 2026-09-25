# -*- coding: utf-8 -*-
"""内部自洽校验量：域内效率、饥饿体积份额、前缘位置（Task 6）。

三个量都为**纯后处理**：只消费 (cement 浓度场, geom 几何字典)，不改任何数值路径。

口径与声明
----------
1. ``domain_eta_e`` 与求解器 ``annulus_d2dga.py`` 内的 ``effective_efficiency``
   （= ``bulk_cement_fill`` = ``_trapez2d(b·cement, geom)/_trapez2d(b, geom)``）
   **同一口径**，故可用它反向复核 metrics 末行——二者是同一条表达式的两次求值，
   实测浮点一致到 ~1e-16。
2. ``η_E ≡ 1 − 饥饿体积份额`` 是**代数恒等**，对任意连续浓度场可严格展开为

       (1 − η_E) − 饥饿份额 ≡ Σ_{c≥0.5} b(1−c)/∬b − Σ_{c<0.5} b·c/∬b

   （即过渡带的代数和，可正可负）；二值场（c̄ ∈ {0,1}）时该代数和恒为 0，
   故二值场下恒等式严格成立。独立复算的浮点差约 1e-17~1e-16
   （hu101 现场观测到的 0.0 是该场的巧合）——**不是**可对连续场写硬相等断言的等式。
   历史口径"偏差 ≤0.007"**出自基准算例**（纯水泥恒定入口 ⇒ 场近乎二值），
   **不是**现场井；现场井实测 0.006~0.043。见
   ``tests/contract/test_internal_consistency.py`` 的"二值场严格 / 非二值场有空隙"两条断言。
3. ``front_position_m`` 返回**模型 s 口径**的前缘位置（自鞋口沿环空的距离，m）：
   b 加权方位列均值首达 ``level`` 的最大 ``s``；未达时返回 ``0.0``
   （与求解器 ``_front``（``annulus_d2dga.py``）的 ``else 0.0`` 同约定，
   只是把单条方位线换成全环空的 b 加权列均值）。**方向约定**：环空求解域自
   ``s=0``（鞋口、入口）向 ``s=L``（悬挂器侧、出口）展开，故"前缘推进"表现为 ``s`` 增大；
   测深 ``md = bottom_md − s`` 只是**派生标签**（``md`` 随 ``s`` 递减），
   本函数**不再**返回 md 口径的量（旧 md 口径会恒等于域底，见 task-6-report.md 修复轮 1）。
4. 本模块**不**在顶层导入 ``annulus_d2dga``：``annulus_d2dga`` 会反向导入
   ``cemdisp.diagnostics.displacement_metrics``，顶层互导会在把本模块加入
   ``cemdisp/diagnostics/__init__.py`` 时触发循环 ImportError。故 ``_trapez2d``
   一律**函数内延迟导入**。
"""
from __future__ import annotations

from typing import Any, Mapping

import numpy as np
from numpy.typing import NDArray


def domain_eta_e(cement: NDArray, geom: Mapping[str, Any]) -> float:
    """域内体积加权效率 η_E = ∬b·c / ∬b（与 ``_evaluation_window_efficiencies`` 同口径）。"""
    from cemdisp.models2d.annulus_d2dga import _trapez2d  # 延迟导入，见模块 docstring 第 4 条

    b = geom["b"]
    return float(_trapez2d(b * np.asarray(cement, float), geom)
                 / max(_trapez2d(b, geom), 1e-12))


def starved_volume_fraction(cement: NDArray, geom: Mapping[str, Any],
                            threshold: float = 0.5) -> float:
    """饥饿体积份额 ∬b·1[c<threshold] / ∬b（口径同 2026-09-11 探针脚本）。"""
    from cemdisp.models2d.annulus_d2dga import _trapez2d  # 延迟导入，见模块 docstring 第 4 条

    b = geom["b"]
    ind = (np.asarray(cement, float) < threshold).astype(float)
    return float(_trapez2d(b * ind, geom) / max(_trapez2d(b, geom), 1e-12))


def front_position_m(cement: NDArray, geom: Mapping[str, Any], level: float = 0.5) -> float:
    """模型 s 口径前缘位置（m，**自鞋口沿环空**）：b 加权方位列均值首达 ``level`` 的最大 ``s``。

    未达阈值时返回 ``0.0``（"前缘尚未离开鞋口"，与求解器 ``_front`` 的兜底同约定）。
    方向约定见模块 docstring 第 3 条：``md = bottom_md − s`` 只是派生标签；
    环空入口在 ``s=0``（鞋口）侧，故水泥入环空后前缘自 0 向上增长。

    注：列均值的权重用 ``b`` 的方位向**和**（生产/诊断网格的 ``y`` 均为均匀 ``linspace``，
    与梯形权重只差同一常数因子，比值不变）。
    """
    b = np.asarray(geom["b"], float)
    c = np.asarray(cement, float)
    col = (b * c).sum(axis=0) / np.maximum(b.sum(axis=0), 1e-30)
    s = np.asarray(geom["s"], float)
    reached = np.where(col >= level)[0]
    return float(s[reached.max()]) if reached.size else 0.0
