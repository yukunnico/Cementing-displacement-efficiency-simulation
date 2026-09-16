"""B-4：L1 硬约束——牛顿极限（τ_y=0、n=1）下新开关不改变结果。

本文件是阶段 B（HB 闭包改造）的**回归护栏**，守护设计规格 §6 的 L1 锚：

1. ``wall`` 恒零（等价于 τ_y=0 ⇒ 无冻结）时，``solve_stream_function`` 的
   输出与不传 ``wall`` 逐位一致——Task 2 的屈服门是 opt-in、默认无扰；
2. ``closure=None``（⇒ ``NewtonianClosure``）时，默认参数路径的输出与一份
   **测试内固化的期望数组**逐位一致。

第 2 条是 Task 1 评审新增要求：既有测试只比"新代码内部两条路径互等"
（如 ``test_stream_yield_gate.py`` 里 base vs wall=0），若有人改了 I₁ 闭包
**之外**的东西（源项装配、面调和平均、BC 行、速度换算前的 Ψ 归一…），
两条路径会**同步**漂移而不报警。固化快照把"当前 solver 的绝对输出"钉死，
任何非预期改动都会立刻变红。

⚠️ 快照锚的层次说明：按 ``tests/README.md``，'记录某次实现的输出数值' 属
``history/`` 层判据。本文件整体放 ``contract/``（L1 是规格级硬约束），但
``test_frozen_snapshot_*`` 这一条本质是 history 风格的漂移探测器——**预期
在口径有意变更时改写**（改写时须更新下方 ``_FROZEN_GENERATED_AT`` 注记）。

References
----------
Zhang & Frigaard (2022), *JFM* **947**, A32：(4.21a)/(4.22)/(2.2)。
"""
from __future__ import annotations

import numpy as np

from cemdisp.models2d.stream_function import solve_stream_function


def _geom_arrays(ny=9, nz=4):
    return {"phi": np.linspace(0.0, 1.0, ny), "s": np.linspace(0.0, 30.0, nz),
            "H": np.full((ny, nz), 0.01), "hole_mm": np.full(nz, 215.9),
            "od_mm": np.full(nz, 168.3), "y": np.linspace(0.0, np.pi * 0.1, ny),
            "b": np.full((ny, nz), 0.02), "inc_deg": np.zeros(nz)}


def test_zero_yield_wall_is_bitwise_unchanged():
    """τ_y=0 ⇒ wall 恒 0 ⇒ 与不传 wall 逐位一致。"""
    g = _geom_arrays(); c = np.full((9, 4), 0.5); b = np.zeros((2, 9, 4))
    base = solve_stream_function(g, c, 0.058, 0.171, 0.34, b)
    gated = solve_stream_function(g, c, 0.058, 0.171, 0.34, b, wall=np.zeros((9, 4)))
    assert np.array_equal(base, gated)


def _synthetic_case():
    """固定合成输入：非均匀半隙场 + 携带 φ-/ξ-梯度的 b 场（确定性，无 RNG）。"""
    ny, nz = 9, 4
    phi = np.linspace(0.0, 1.0, ny)
    s = np.linspace(0.0, 30.0, nz)
    H = 0.008 + 0.004 * phi[:, None] + 0.001 * np.arange(nz)[None, :] / float(nz)
    geom = {"phi": phi, "s": s, "H": H,
            "hole_mm": np.linspace(215.0, 220.0, nz),
            "od_mm": np.full(nz, 168.3),
            "y": phi * np.pi * 0.1,
            "b": 2.0 * H,
            "inc_deg": np.linspace(0.0, 5.0, nz)}
    # c̄ 沿 φ 线性 0.1→0.9 ⇒ I₁ 场非平凡；两层黏度经标量传入
    c = np.linspace(0.1, 0.9, ny)[:, None] * np.ones((1, nz))
    b = np.zeros((2, ny, nz))
    b[0] = 0.02 * (1.0 + phi[:, None]) * (1.0 + 0.1 * np.arange(nz)[None, :])   # b_φ（竖直井机制）
    b[1] = 0.03 * np.sin(np.pi * phi)[:, None] * np.ones((1, nz))               # b_ξ（方位机制）
    return geom, c, b


# 生成条件（Task 4 评审要求注明）：
#   commit 611dcce（B-3 修复轮 1 / 阶段 B HEAD），closure=None、wall=None 的默认路径，
#   输入 = _synthetic_case()，调用 = solve_stream_function(geom, c, 0.058, 0.171, 0.34, b)。
#   数值为 float64 的 repr 全精度往返（Python float repr 保证 round-trip）。
_FROZEN_GENERATED_AT = "611dcce (feat/hb-closure-phase-b, B-3 修复轮 1)"
_FROZEN_PSI = (
    (0.0, 0.0, 0.0, 0.0),
    (0.08255340426295936, 0.08387390830467853, 0.0851417976512296, 0.0863599486888221),
    (0.1801347625608875, 0.18252049798037448, 0.184806594187366, 0.1869988456349786),
    (0.29318016639204963, 0.29628958701424085, 0.2992635252349338, 0.30211027826673426),
    (0.4210434589368052, 0.4244550985240726, 0.4277122967677184, 0.43082493671895844),
    (0.5615142604723655, 0.5647477753041077, 0.567829773771516, 0.5707703241456674),
    (0.7102355217233317, 0.7127925162396236, 0.7152259118819857, 0.7175441933760892),
    (0.8600118683828302, 0.8614422173217785, 0.8628014772631848, 0.864094667249213),
    (1.0, 1.0, 1.0, 1.0),
)


def test_frozen_snapshot_stream_function_defaults():
    """当前 solver（默认 closure=None/wall=None）逐位等于冻结锚。

    冻结条件见 ``_FROZEN_GENERATED_AT``。本断言守护的不只是"新开关无扰"，
    而是**整条默认路径的绝对输出**——I₁ 闭包、源项装配、面调和平均、BC 行、
    单位通量归一任一被改动都会逐位变红。
    """
    geom, c, b = _synthetic_case()
    psi = solve_stream_function(geom, c, 0.058, 0.171, 0.34, b)
    frozen = np.array(_FROZEN_PSI, dtype=np.float64)
    assert psi.shape == frozen.shape
    # 逐位：不接受任何浮点容差
    assert np.array_equal(psi, frozen), (
        "默认路径输出偏离冻结锚（生成于 " + _FROZEN_GENERATED_AT + "）："
        "若本次变更**有意**改动了 (4.22) 装配/闭包/归一化，请同步更新锚并注明新条件；"
        "否则即为回归。"
    )
