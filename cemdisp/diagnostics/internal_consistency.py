# -*- coding: utf-8 -*-
"""内部自洽校验量：域内效率、饥饿体积份额、前缘位置（Task 6）。

三个量都为**纯后处理**：只消费 (cement 浓度场, geom 几何字典)，不改任何数值路径。

口径与声明
----------
1. ``domain_eta_e`` 与求解器 ``annulus_d2dga.py`` 内的 ``effective_efficiency``
   （= ``bulk_cement_fill`` = ``_trapez2d(b·cement, geom)/_trapez2d(b, geom)``）
   **同一口径**，故可用它反向复核 metrics 末行，二者应逐位一致。
2. ⚠️ 恒等式 ``η_E ≡ 1 − 饥饿体积份额`` **只在 c̄ ∈ {0,1} 时严格成立**；
   生产井浓度场连续，故生产口径下该式是**带阈值的记录**（历史实测偏差 ≤0.007），
   **不是**可硬断言的等式。见 ``tests/contract/test_internal_consistency.py``
   的"二值场严格 / 非二值场有空隙"两条断言。
3. ``front_position_m`` 返回的是**钻井液度量（md）下已到水泥的最深点**，
   见其 docstring 的方向说明——环空求解域自 ``s=0``（鞋口、入口）向 ``s=L``
   （悬挂器侧、出口）展开，故该量在"水泥自鞋口入环空"的物理下**恒等于域底 md**
   （水泥一旦入环空，``z=0`` 列即被点亮）。若需"前缘推进距离"，
   应取求解器 metrics 的 ``front_wide_m``/``front_mid_m``/``front_narrow_m``
   （口径为 ``s = max{s : c̄ ≥ 0.5}``，单位 m，自鞋口起算）。
"""
from __future__ import annotations

from typing import Any, Mapping

import numpy as np
from numpy.typing import NDArray

from cemdisp.models2d.annulus_d2dga import _trapez2d


def domain_eta_e(cement: NDArray, geom: Mapping[str, Any]) -> float:
    """域内体积加权效率 η_E = ∬b·c / ∬b（与 ``_evaluation_window_efficiencies`` 同口径）。"""
    b = geom["b"]
    return float(_trapez2d(b * np.asarray(cement, float), geom)
                 / max(_trapez2d(b, geom), 1e-12))


def starved_volume_fraction(cement: NDArray, geom: Mapping[str, Any],
                            threshold: float = 0.5) -> float:
    """饥饿体积份额 ∬b·1[c<threshold] / ∬b（口径同 2026-09-11 探针脚本）。"""
    b = geom["b"]
    ind = (np.asarray(cement, float) < threshold).astype(float)
    return float(_trapez2d(b * ind, geom) / max(_trapez2d(b, geom), 1e-12))


def front_position_m(cement: NDArray, geom: Mapping[str, Any], level: float = 0.5) -> float:
    """前缘位置（m）：b 加权列均值首达 ``level`` 的最深 md；未达则返回域底 md。

    方向约定（勿误读）：``md = bottom_md − s``，故 ``md`` 随 ``z`` 索引**递减**，
    ``md[0] = bottom_md`` 为鞋口（入口）侧、``md[-1] = top_md`` 为悬挂器（出口）侧。
    本函数取 ``md[reached.min()]``＝已点亮水泥的**最深**深度坐标。
    """
    b = np.asarray(geom["b"], float)
    col = (b * np.asarray(cement, float)).sum(axis=0) / np.maximum(b.sum(axis=0), 1e-30)
    reached = np.where(col >= level)[0]
    md = np.asarray(geom["md"], float)
    return float(md[reached.min()]) if reached.size else float(md.max())
