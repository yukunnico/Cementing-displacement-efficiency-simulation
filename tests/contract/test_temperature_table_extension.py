# -*- coding: utf-8 -*-
"""Phase 4d 合同测试：温度表版本感知加载 + stop_t>表末 批次级前置断言。

原则：
- 交付表加载路径（from_files / load_delivered_pair）逐位不变 => 用真实交付表钉桩；
- 扩展表加载走新增 load_extended_pair_4d（_ext4d 溯源标识 + 时间轴侧车 + 形状硬校验）；
- assert_time_table_coverage 只在批驱动前置门响亮报错，不改 T() 的 clamp+审计行为。
"""
from __future__ import annotations

import zipfile
from pathlib import Path

import numpy as np
import pytest

from cemdisp.data.temperature_field import (
    EXPECTED_SHAPE,
    EXTENDED_SHAPE_4D,
    EXTENDED_TABLE_MARKER,
    TABLE_VERSIONS,
    TemperatureTableCoverageError,
    TableTemperatureField,
    assert_time_table_coverage,
    load_delivered_pair,
    load_extended_pair_4d,
)

_STOP_T_R06_S = 19016.017281054297  # 呼1-004 Toff_rate_x0.6 判别量 stop_t（T2 批）
_N_DEPTH, _N_EXT = 333, 362


# ---------------------------------------------------------------------------
# 合成件工具（无 openpyxl 环境：手写最小 OOXML）
# ---------------------------------------------------------------------------
def _col_letter(idx0: int) -> str:
    s = ""
    n = idx0 + 1
    while n:
        n, rem = divmod(n - 1, 26)
        s = chr(65 + rem) + s
    return s


def _write_minimal_xlsx(path: Path, matrix: np.ndarray) -> None:
    rows = []
    for r in range(matrix.shape[0]):
        cells = "".join(
            '<c r="' + _col_letter(c) + str(r + 1) + '"><v>' + repr(float(matrix[r, c])) + "</v></c>"
            for c in range(matrix.shape[1])
        )
        rows.append('<row r="' + str(r + 1) + '">' + cells + "</row>")
    ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    sheet = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<worksheet xmlns="' + ns + '"><sheetData>' + "".join(rows) + "</sheetData></worksheet>"
    )
    ct = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        '<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        "</Types>"
    )
    rels = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
        "</Relationships>"
    )
    wb = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<workbook xmlns="' + ns + '"><sheets><sheet name="Sheet1" sheetId="1"/></sheets></workbook>'
    )
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("[Content_Types].xml", ct)
        z.writestr("_rels/.rels", rels)
        z.writestr("xl/workbook.xml", wb)
        z.writestr("xl/worksheets/sheet1.xml", sheet)


def _synthetic_ext_bundle(tmp_path: Path, marker_ok: bool = True):
    """构造 333x362 扩展表合成件：两表井底行严格同值（=最深行本身相同）。"""
    depth = np.linspace(30.0, 7660.0, _N_DEPTH)
    axis_min = np.concatenate(
        [np.arange(199.0), [198.792801055535], np.arange(199.0, 361.0)]
    )
    assert axis_min.size == _N_EXT
    t_in = depth[:, None] * 0.01 + axis_min[None, :] * 0.05 + 16.0
    t_out = t_in.copy()  # 井底行同值天然成立；上部行两表本可异值，同值也不违例
    tag = "ext4d" if marker_ok else "wrong"
    p_in = tmp_path / ("T_in_" + tag + ".xlsx")
    p_out = tmp_path / ("T_out_" + tag + ".xlsx")
    _write_minimal_xlsx(p_in, t_in)
    _write_minimal_xlsx(p_out, t_out)
    ax = tmp_path / "T_ext4d_time_axis_min.csv"
    ax.write_text("\n".join(repr(float(v)) for v in axis_min), encoding="utf-8")
    dp = tmp_path / "depth.csv"
    dp.write_text(
        "depth_well_logging_m_\n" + "\n".join(repr(float(v)) for v in depth),
        encoding="utf-8",
    )
    return p_in, p_out, ax, dp, t_in, t_out, depth, axis_min


# ---------------------------------------------------------------------------
# 登记表
# ---------------------------------------------------------------------------
def test_version_registry_pinned():
    assert EXPECTED_SHAPE == (333, 200)          # 交付口径原样保留
    assert EXTENDED_SHAPE_4D == (333, 362)
    assert TABLE_VERSIONS == {"delivered": (333, 200), "extended_4d": (333, 362)}
    assert EXTENDED_TABLE_MARKER == "_ext4d"


# ---------------------------------------------------------------------------
# 交付表路径逐位不变（真实交付件钉桩）
# ---------------------------------------------------------------------------
def test_delivered_load_path_unchanged():
    t_in, t_out = load_delivered_pair(use_cache=False)
    assert t_in.table.shape == EXPECTED_SHAPE == t_out.table.shape
    np.testing.assert_allclose(t_in.time_s, np.arange(200) * 60.0)  # 均匀理想化轴不变
    assert float(t_in.time_s[-1]) == 11940.0
    assert t_in.depth_m[0] == pytest.approx(30.0)
    assert t_in.depth_m[-1] == pytest.approx(7660.0)


def test_delivered_clamp_audit_unchanged():
    """stop_t 超表末时：T() 仍静默 clamp+计数（行为不变），响亮报错只在批次门。"""
    t_in, _ = load_delivered_pair(use_cache=False)
    v_clamped = t_in.T(7000.0, _STOP_T_R06_S)
    v_at_end = t_in.T(7000.0, 11940.0)
    assert v_clamped == v_at_end          # clamp 到末列
    assert t_in.oob_count == 1            # 计数行为不变
    t_in.reset_audit()


# ---------------------------------------------------------------------------
# 扩展表加载
# ---------------------------------------------------------------------------
def test_extended_load_roundtrip(tmp_path):
    p_in, p_out, ax, dp, t_in_m, t_out_m, depth, axis_min = _synthetic_ext_bundle(tmp_path)
    f_in, f_out = load_extended_pair_4d(p_in, p_out, ax, depth_csv_path=dp)
    assert f_in.table.shape == EXTENDED_SHAPE_4D
    np.testing.assert_allclose(f_in.table, t_in_m, rtol=0, atol=0)   # 逐位回读
    np.testing.assert_allclose(f_out.table, t_out_m, rtol=0, atol=0)
    np.testing.assert_allclose(f_in.time_s, axis_min * 60.0, rtol=0, atol=0)
    assert float(f_in.time_s[-1]) == 21600.0
    assert float(f_in.time_s[199]) == pytest.approx(198.792801055535 * 60.0)
    # 跨施工终点节点插值可用且不越界（非均匀轴被正确消费）
    v = f_in.T(7000.0, 198.8 * 60.0)
    assert f_in.oob_count == 0
    assert 16.0 < v < 200.0


def test_extended_requires_marker(tmp_path):
    p_in, p_out, ax, dp, *_ = _synthetic_ext_bundle(tmp_path, marker_ok=False)
    with pytest.raises(ValueError, match="溯源标识"):
        load_extended_pair_4d(p_in, p_out, ax, depth_csv_path=dp)


def test_extended_rejects_illegal_shape(tmp_path):
    depth = np.linspace(30.0, 7660.0, _N_DEPTH)
    bad = np.zeros((_N_DEPTH, 361))          # 既非交付也非扩展口径
    p_in = tmp_path / "T_in_ext4d.xlsx"
    p_out = tmp_path / "T_out_ext4d.xlsx"
    _write_minimal_xlsx(p_in, bad)
    _write_minimal_xlsx(p_out, bad)
    ax = tmp_path / "T_ext4d_time_axis_min.csv"
    # 侧车给足 362（合法），让"形状"成为被检门（侧车长度门在其前，另行覆盖）
    ax.write_text("\n".join(repr(float(v)) for v in np.arange(362.0)), encoding="utf-8")
    dp = tmp_path / "depth.csv"
    dp.write_text(
        "depth_well_logging_m_\n" + "\n".join(repr(float(v)) for v in depth),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="形状"):
        load_extended_pair_4d(p_in, p_out, ax, depth_csv_path=dp)
    # 交付口径硬校验同样拒绝扩展形状（EXPECTED_SHAPE 单点不放宽，F-4）
    with pytest.raises(ValueError, match="形状"):
        TableTemperatureField(np.zeros((_N_DEPTH, _N_EXT)), depth)


def test_extended_rejects_bad_axis(tmp_path):
    p_in, p_out, ax, dp, *_ = _synthetic_ext_bundle(tmp_path)
    # 非单调
    ax.write_text("\n".join(repr(float(v)) for v in np.arange(362.0)[::-1]), encoding="utf-8")
    with pytest.raises(ValueError, match="单调"):
        load_extended_pair_4d(p_in, p_out, ax, depth_csv_path=dp)
    # 长度不匹配
    ax.write_text("0\n1\n2", encoding="utf-8")
    with pytest.raises(ValueError, match="侧车长度"):
        load_extended_pair_4d(p_in, p_out, ax, depth_csv_path=dp)


# ---------------------------------------------------------------------------
# 批次级前置断言
# ---------------------------------------------------------------------------
def test_coverage_assert_passes_for_extended(tmp_path):
    p_in, p_out, ax, dp, *_ = _synthetic_ext_bundle(tmp_path)
    f_in, f_out = load_extended_pair_4d(p_in, p_out, ax, depth_csv_path=dp)
    assert_time_table_coverage([f_in, f_out], _STOP_T_R06_S, label="呼1-004_r0.6")  # 不抛
    assert_time_table_coverage(f_in, _STOP_T_R06_S)  # 单场对象也可


def test_coverage_assert_raises_for_delivered(tmp_path):
    # 造 200 列表末 11940 s 的场：模拟交付表覆盖不住 r0.6
    f = TableTemperatureField(
        np.zeros((_N_DEPTH, 200)),
        np.linspace(30.0, 7660.0, _N_DEPTH),
        np.arange(200) * 60.0,
    )
    with pytest.raises(TemperatureTableCoverageError) as ei:
        assert_time_table_coverage([f], _STOP_T_R06_S, label="呼1-004_r0.6")
    msg = str(ei.value)
    assert "19016.0" in msg and "11940.0" in msg and "呼1-004_r0.6" in msg
    # 报错后 T() clamp 行为不受影响
    assert f.T(7000.0, _STOP_T_R06_S) == f.T(7000.0, 11940.0)


def test_coverage_assert_rejects_non_table_field():
    class _NoAxis:
        pass

    with pytest.raises(TypeError):
        assert_time_table_coverage([_NoAxis()], 100.0)
