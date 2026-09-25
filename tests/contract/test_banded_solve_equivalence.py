# -*- coding: utf-8 -*-
"""带状 Cholesky 与 spsolve 解同一离散系统 ⇒ Ψ 相对差应为舍入级。

口径（Task 8 简报 + 协调者裁定 R5）
----------------------------------
- 本文件是**验收权威**：``banded=True``（新默认）与 ``banded=False``（历史
  spsolve 口径）必须解**同一个**离散系统，Ψ 相对差 < 1e-12、且端点 Dirichlet
  口径（φ=0 ⇒ Ψ=0、φ=1 ⇒ Ψ=1）逐位不变。达不到该容差即停下上报，不放宽阈值。
- 简报给出的 ``_solve_banded_interior`` 是**草图**（R5）：其 a_face/c_face 少了
  ``1/r_a`` 因子、Dirichlet 反代项亦同 ⇒ 本条按 R5 允许"改用更稳的带状装配"，
  并以 ``test_banded_storage_matches_sparse_matrix_bands`` 做**逐带核对**——
  把带状装配与现行稀疏装配的内部子块**逐元素**比对（不重写一份系数公式，
  而是捕获 ``solve_stream_function`` 真实装配出的稀疏矩阵作为参考）。
- 计时（R5 第 5 条）自带于 ``test_timing_spsolve_vs_banded_single_solve``：
  打印同机同参数单次耗时对照（生产网格 40×250），并断言带状严格更快
  （≥3× 的目标由报告里的人工核对给出；此处只设一个抗噪声的单调性下限，
  避免把性能抖动做成 CI 的假红）。
"""
import time

import numpy as np
import pytest

import cemdisp.models2d.stream_function as sf
from cemdisp.models2d.stream_function import solve_stream_function


def _geom(ny=40, nz=250, seed=0):
    """均匀 ξ/φ 网格、随机（非均匀）半隙 H —— 与生产口径同构的最小几何。"""
    rng = np.random.default_rng(seed)
    phi = np.linspace(0.0, 1.0, ny)
    H = 0.01 + 0.002 * rng.random((ny, nz))
    return {"phi": phi, "s": np.linspace(0.0, 2468.0, nz), "H": H,
            "hole_mm": np.full(nz, 215.9), "od_mm": np.full(nz, 168.3),
            "y": phi * np.pi * 0.1, "b": 2.0 * H}


def _c_field(g):
    ny, nz = g["H"].shape
    return np.clip(0.5 + 0.4 * np.sin(np.linspace(0, 3, nz))[None, :]
                   * np.ones((ny, 1)), 0, 1)


def _b_phi_only(g):
    """仅 φ-槽（简报版）：源项只有 φ-梯度 ⇒ 主要激励 φ-向耦合。"""
    ny, nz = g["H"].shape
    b = np.zeros((2, ny, nz))
    b[0] = 0.1 * np.cos(np.pi * np.asarray(g["phi"]))[:, None]
    return b


def _b_full(g):
    """两槽都随 (φ, ξ) 变化 ⇒ 同时激励 φ-向与 ξ-向（最远一条带）耦合。"""
    ny, nz = g["H"].shape
    phi = np.asarray(g["phi"])[:, None]
    xi = np.linspace(0.0, 1.0, nz)[None, :]
    b = np.zeros((2, ny, nz))
    b[0] = 0.1 * np.cos(np.pi * phi) * (1.0 + 0.7 * xi)
    b[1] = 0.08 * np.sin(np.pi * phi) * np.cos(2.0 * np.pi * xi)
    return b


@pytest.mark.parametrize("ny,nz", [(9, 8), (40, 60)])
def test_banded_matches_sparse(ny, nz):
    """简报 Step 1 的逐字验收（两组网格规模）。"""
    g = _geom(ny, nz)
    c = _c_field(g)
    psi_sp = solve_stream_function(g, c, 0.058, 0.171, 0.34, _b_phi_only(g),
                                   banded=False)
    psi_bd = solve_stream_function(g, c, 0.058, 0.171, 0.34, _b_phi_only(g),
                                   banded=True)
    assert np.all(np.isfinite(psi_bd))
    rel = np.max(np.abs(psi_bd - psi_sp)) / max(np.max(np.abs(psi_sp)), 1e-30)
    assert rel < 1e-12
    # 端点 Dirichlet 口径不变
    assert psi_bd[0, :] == pytest.approx(0.0, abs=1e-14)
    assert psi_bd[-1, :] == pytest.approx(1.0, abs=1e-14)


@pytest.mark.parametrize("ny,nz", [(9, 8), (40, 60), (40, 250)])
def test_banded_matches_sparse_general_source(ny, nz):
    """两槽都随 (φ, ξ) 变化（最远带耦合被真实激励）下的等价性。"""
    g = _geom(ny, nz)
    c = _c_field(g)
    b = _b_full(g)
    psi_sp = solve_stream_function(g, c, 0.058, 0.171, 0.34, b, banded=False)
    psi_bd = solve_stream_function(g, c, 0.058, 0.171, 0.34, b, banded=True)
    rel = np.max(np.abs(psi_bd - psi_sp)) / max(np.max(np.abs(psi_sp)), 1e-30)
    assert rel < 1e-12


@pytest.mark.parametrize("ny,nz", [(9, 8), (40, 60)])
def test_banded_matches_sparse_with_wall(ny, nz):
    """屈服门冻结（wall≠None）⇒ 算子各向异性被拉到 1/地板 = 1e6，
    正是生产强偏心情形；此路径下等价性同样必须成立。"""
    g = _geom(ny, nz)
    ny_, nz_ = g["H"].shape
    c = _c_field(g)
    b = _b_full(g)
    # φ 靠近窄边（φ→1）的一整块被冻结到地板 + 随机部分冻结
    rng = np.random.default_rng(1)
    wall = np.zeros((ny_, nz_))
    wall[3 * ny_ // 4:, :] = 1.0
    wall += 0.3 * rng.random((ny_, nz_))
    wall = np.clip(wall, 0.0, 1.0)
    psi_sp = solve_stream_function(g, c, 0.058, 0.171, 0.34, b, wall=wall,
                                   banded=False)
    psi_bd = solve_stream_function(g, c, 0.058, 0.171, 0.34, b, wall=wall,
                                   banded=True)
    rel = np.max(np.abs(psi_bd - psi_sp)) / max(np.max(np.abs(psi_sp)), 1e-30)
    assert rel < 1e-12


def test_banded_storage_matches_sparse_matrix_bands(monkeypatch):
    """逐带核对（R5）：带状装配的 (ab, rhs) 必须与现行稀疏装配的内部子块逐元素一致。

    做法：先让 ``banded=False`` 跑一遍并用 spy 捕获它真实装配出的稀疏矩阵与右端，
    再让 ``banded=True`` 跑一遍并用 spy 捕获带状装配产物；两者按
    **同一套 COO 行序** ``row = i*nz + j`` ↔ ``m = j*(ny-2) + (i-1)`` 对照。

    参考量（全部来自现行稀疏装配，不重写任何系数公式）：
      * 内部子块 ``A_int = A_full[row_int][:, row_int]``；带状系统解的是**取负**形式
        （现行全矩阵对角为负的 M-矩阵 ⇒ 取负后为对称正定 Laplacian）；
      * Dirichlet 反代：``rhs_banded = -rhs_sparse[row_int] + Σ_dir A_full[row_int, col_dir]·bc``。
    """
    captured = {}
    real_spsolve = sf.spla.spsolve

    def _spy_spsolve(A, b, *args, **kwargs):
        captured["A"] = A.tocsr()
        captured["b"] = np.asarray(b, dtype=float).copy()
        return real_spsolve(A, b, *args, **kwargs)

    monkeypatch.setattr(sf.spla, "spsolve", _spy_spsolve)

    real_assemble = sf._assemble_banded_interior

    def _spy_assemble(*args, **kwargs):
        ab, rhs, n_int = real_assemble(*args, **kwargs)
        captured["ab"] = ab.copy()          # 必须在 cholesky 就地改写前拷出
        captured["rhs"] = rhs.copy()
        captured["n_int"] = n_int
        return ab, rhs, n_int

    monkeypatch.setattr(sf, "_assemble_banded_interior", _spy_assemble)

    ny, nz = 7, 5
    g = _geom(ny, nz)
    c = _c_field(g)
    b = _b_full(g)
    solve_stream_function(g, c, 0.058, 0.171, 0.34, b, banded=False)
    A_full = captured["A"].toarray()
    rhs_sparse = captured["b"]
    solve_stream_function(g, c, 0.058, 0.171, 0.34, b, banded=True)
    ab, rhs_bd, n_int = captured["ab"], captured["rhs"], captured["n_int"]

    assert n_int == ny - 2
    row_int = np.ravel([i * nz + j for j in range(nz) for i in range(1, ny - 1)])
    A_int = A_full[np.ix_(row_int, row_int)]

    # ① 稀疏装配的内部子块必须精确对称（带状 Cholesky 的前提）
    assert np.array_equal(A_int, A_int.T), "内部子块不精确对称，带状 Cholesky 前提不成立"

    # ② 逐带核验：下三角（含对角）逐元素等于 −A_int
    band = n_int
    n_unknown = n_int * nz
    assert ab.shape == (band + 1, n_unknown)
    checked = 0
    for m in range(n_unknown):
        for k in range(0, min(band, m) + 1):
            want = -A_int[m, m - k]
            got = ab[k, m - k]
            assert got == pytest.approx(want, rel=1e-14, abs=1e-18), (
                f"带 k={k} 列 {m - k} 处不一致：got={got!r} want={want!r}")
            checked += 1
    assert checked == n_unknown * (band + 1) - band * (band + 1) // 2

    # ③ 上三角槽位（``m+k < n_unknown``）必须保持初值 0。
    #
    # ⚠️ 按 LAPACK 下三角存法 ``ab[i-j, j] == a[i,j] (i>=j)``，只有 ``m < k`` 的
    # (k, m) 槽位才**该**未写；但**构造时不能断言它们为 0**——现行装配用的是
    # "写入切片 ``ab[1, :n-1] = band1.ravel('F')[1:]``"这一形式，它把
    # ``band1[0]``（i=1 行的 φ⁻ 槽，结构上应为 0）挤到了 ``ab[1, 0]``…
    # 即 **m < k 的外侧填充区里会残留合法但无意义的邻带值**（Cholesky 只读下三角、
    # 不读外侧填充，故不影响解）。真正有判别力的是 ②——**逐带逐元素**核对
    # ``ab[k, m-k] == −A_int[m, m-k]``（已覆盖全部下三角元素，含必须为 0 的
    # 结构性零槽位：如偏移 1 在 i=1 行、偏移 2..n_int-1 的全部格子）。
    # 故此处只做一条弱断言：至少有一条真带被写过（防"全 0 也能过 ②"的退化解），
    # 且对角槽位全非零。
    assert np.all(ab[0, :] != 0.0), "对角带存在 0（内部子块不可能有零对角）"
    assert np.all(ab[n_int, :n_unknown - n_int] != 0.0), (
        "最远带（ξ⁻ 邻格）全 0：排序/偏移装配可能整条丢失")

    # ④ 右端：内部源项 + φ=1 端 Dirichlet（值 1）反代
    src_int = -rhs_sparse[row_int]
    dir_cols = np.concatenate([np.arange(nz), (ny - 1) * nz + np.arange(nz)])
    bc_vals = np.concatenate([np.zeros(nz), np.ones(nz)])
    want_rhs = src_int + A_full[np.ix_(row_int, dir_cols)] @ bc_vals
    assert np.allclose(rhs_bd, want_rhs, rtol=1e-14, atol=1e-18)
    # φ=0 端（bc=0）不贡献右端项——显式确认该端确实进了对角而非右端
    assert np.allclose(A_full[np.ix_(row_int, np.arange(nz))] @ np.zeros(nz), 0.0)


def test_banded_endpoint_dirichlet_is_exact():
    """端点口径逐位：Ψ(φ=0) 恰为 0.0、Ψ(φ=1) 恰为 1.0（不是 "接近"）。"""
    g = _geom(13, 9)
    psi = solve_stream_function(g, _c_field(g), 0.058, 0.171, 0.34,
                                _b_full(g), banded=True)
    assert np.array_equal(psi[0, :], np.zeros(psi.shape[1]))
    assert np.array_equal(psi[-1, :], np.ones(psi.shape[1]))


def _time_once(fn, reps):
    best = float("inf")
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best


def test_timing_spsolve_vs_banded_single_solve(capsys=None):
    """计时对照（R5 第 5 条）：生产网格 40×250 单次求解 spsolve vs 带状 Cholesky。

    取 5 次中的**最小值**（抗调度噪声），打印原始数字；断言带状严格更快
    （单调性下限，抗抖动）。≥3× 的目标由人工核对报告里的原始输出。
    """
    ny, nz = 40, 250
    g = _geom(ny, nz)
    c = _c_field(g)
    b = _b_full(g)
    args = (g, c, 0.058, 0.171, 0.34, b)

    # 预热（首次调用含 scipy 惰性导入/缓存，不得计入）
    solve_stream_function(*args, banded=False)
    solve_stream_function(*args, banded=True)

    t_sp = _time_once(lambda: solve_stream_function(*args, banded=False), 5)
    t_bd = _time_once(lambda: solve_stream_function(*args, banded=True), 5)
    ratio = t_sp / t_bd
    print(f"\n[计时] ny={ny} nz={nz}（各取 5 次最小值）")
    print(f"       spsolve        : {t_sp * 1e3:8.3f} ms")
    print(f"       带状 Cholesky  : {t_bd * 1e3:8.3f} ms")
    print(f"       加速比         : {ratio:8.2f}×")
    assert t_bd < t_sp, (
        f"带状 Cholesky 未快于 spsolve（{t_bd * 1e3:.3f} ms vs {t_sp * 1e3:.3f} ms）"
        "——接线可能回退到了稠密/稀疏路径")
