"""
温度场表数据层（T0-1：温压耦合的温度查询模块）

把交付的 333×200 温度场表封装成 ``T(md_m, t_s)`` 查询对象，供后续求解器每步取温：

- ``ConstantTemperatureField``: 恒温场（关温耦合时的回退默认）
- ``TableTemperatureField``: 表格温度场，双线性插值（深度×时间），越界 clamp + 审计
- ``load_delivered_pair``: 加载交付的管内/环空两表，并校验井底（最深行）共享同值

数据源（交付件）：
- 管内表: 参考文档/温压耦合数据、/T_in.xlsx（333 行 × 200 列纯数值矩阵，无表头）
- 环空表: 参考文档/温压耦合数据、/T_out.xlsx（同上；井底行与 T_in 共享同值）
- 深度轴: 参考文档/温压耦合数据、/HT1-004压力计算/呼1-004井身结构.csv
  第 2 列 ``depth_well_logging_m_``（333 点，30→7660 m）
- 时间轴: 200 列 = 200 min，col0 = 初始时刻，t_j = j min × 60 s（裁定总表）

xlsx 读取不依赖 openpyxl（cementT 环境无该包）：直接用标准库 zipfile + xml 解析
首个工作表的数值单元格。读表后缓存为 npz（默认落在 xlsx 同目录，可用
``cache_path`` 指定），二次加载免解析 xlsx。
"""

from __future__ import annotations

import re
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd

__all__ = [
    "ConstantTemperatureField",
    "TableTemperatureField",
    "ClampEvent",
    "load_delivered_pair",
    "EXPECTED_SHAPE",
    "TIME_STEP_S",
    "DEFAULT_T_IN_XLSX",
    "DEFAULT_T_OUT_XLSX",
    "DEFAULT_DEPTH_CSV",
]

# 仓库根：cemdisp/data/temperature_field.py → parents[2]
PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA_DIR = PROJECT_ROOT / "参考文档" / "温压耦合数据、"
DEFAULT_T_IN_XLSX = DEFAULT_DATA_DIR / "T_in.xlsx"          # 管内
DEFAULT_T_OUT_XLSX = DEFAULT_DATA_DIR / "T_out.xlsx"         # 环空
DEFAULT_DEPTH_CSV = (
    DEFAULT_DATA_DIR / "HT1-004压力计算" / "呼1-004井身结构.csv"
)

EXPECTED_SHAPE: Tuple[int, int] = (333, 200)  # 深度 333 点 × 时间 200 列
TIME_STEP_S: float = 60.0                     # 每列 1 min，col0=初始时刻
DEPTH_CSV_COLUMN = "depth_well_logging_m_"
_SHEET_XML = "xl/worksheets/sheet1.xml"
_CELL_REF_RE = re.compile(r"([A-Z]+)(\d+)")


# ---------------------------------------------------------------------------
# 审计记录
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ClampEvent:
    """一次越界查询的审计记录（记录请求位置与 clamp 后的落点）。"""

    md_m: float          # 请求深度
    t_s: float           # 请求时间
    md_clamped_m: float  # clamp 后深度
    t_clamped_s: float   # clamp 后时间


# ---------------------------------------------------------------------------
# 恒温场
# ---------------------------------------------------------------------------
class ConstantTemperatureField:
    """恒温场：任意 (md, t) 返回常数 T_c（°C）。关温耦合时的回退默认。"""

    def __init__(self, T_c: float) -> None:
        self.T_c = float(T_c)

    def T(self, md_m: float, t_s: float) -> float:
        return self.T_c

    @property
    def oob_count(self) -> int:
        """越界次数（恒温场永不越界，恒为 0，与 TableTemperatureField 同型）。"""
        return 0

    @property
    def oob_events(self) -> Tuple[ClampEvent, ...]:
        return ()

    def reset_audit(self) -> None:
        """清空审计（恒温场无审计，保留同型接口）。"""


# ---------------------------------------------------------------------------
# 表格温度场
# ---------------------------------------------------------------------------
class TableTemperatureField:
    """表格温度场：深度×时间双线性插值，越界 clamp 到表域边界并记审计。

    参数
    ----------
    table : (n_depth, n_time) 纯数值矩阵
    depth_m : 深度轴（m），长度须等于 n_depth，严格单调增
    time_s : 时间轴（s），长度须等于 n_time；缺省按 ``j * TIME_STEP_S`` 生成
    expected_shape : 期望形状，默认交付口径 (333, 200)；仅测试可用小形状
    source : 来源文件路径（可选，仅用于审计/溯源）
    """

    def __init__(
        self,
        table: np.ndarray,
        depth_m: np.ndarray,
        time_s: Optional[np.ndarray] = None,
        *,
        expected_shape: Tuple[int, int] = EXPECTED_SHAPE,
        source: Optional[Path] = None,
    ) -> None:
        arr = np.asarray(table, dtype=float)
        if arr.ndim != 2:
            raise ValueError(f"温度表须为二维矩阵，实际 ndim={arr.ndim}")
        if arr.shape != tuple(expected_shape):
            raise ValueError(
                f"温度表形状应为 {tuple(expected_shape)}（深度×时间），"
                f"实际 {arr.shape}"
            )

        z = np.asarray(depth_m, dtype=float)
        if z.ndim != 1:
            raise ValueError(f"深度轴须为一维，实际 ndim={z.ndim}")
        if len(z) != arr.shape[0]:
            raise ValueError(
                f"深度轴长度应为 {arr.shape[0]}，实际 {len(z)}"
            )
        if len(z) < 2:
            raise ValueError(f"深度轴至少需 2 个节点，实际 {len(z)}")
        if np.isnan(z).any():
            raise ValueError(f"深度轴含 NaN（{int(np.isnan(z).sum())} 个）")
        if not np.all(np.diff(z) > 0):
            raise ValueError("深度轴非严格单调递增（须随深度增加）")

        if time_s is None:
            t = np.arange(arr.shape[1], dtype=float) * TIME_STEP_S
        else:
            t = np.asarray(time_s, dtype=float)
            if t.ndim != 1:
                raise ValueError(f"时间轴须为一维，实际 ndim={t.ndim}")
            if len(t) != arr.shape[1]:
                raise ValueError(
                    f"时间轴长度应为 {arr.shape[1]}（=温度表列数），实际 {len(t)}"
                )
            if len(t) < 2:
                raise ValueError(f"时间轴至少需 2 个节点，实际 {len(t)}")
            if np.isnan(t).any():
                raise ValueError(f"时间轴含 NaN（{int(np.isnan(t).sum())} 个）")
            if not np.all(np.diff(t) > 0):
                raise ValueError("时间轴非严格单调递增（须随时间增加）")

        if np.isnan(arr).any():
            raise ValueError(f"温度表含 NaN（{int(np.isnan(arr).sum())} 个）")

        self.table = arr
        self.depth_m = z
        self.time_s = t
        self.source = Path(source) if source is not None else None
        self._oob_events: List[ClampEvent] = []

    # -- 审计 ---------------------------------------------------------------
    @property
    def oob_count(self) -> int:
        """越界（被 clamp）查询次数。"""
        return len(self._oob_events)

    @property
    def oob_events(self) -> Tuple[ClampEvent, ...]:
        """越界查询记录（按发生顺序）。"""
        return tuple(self._oob_events)

    def reset_audit(self) -> None:
        """清空越界审计。"""
        self._oob_events.clear()

    # -- 查询 ---------------------------------------------------------------
    def T(self, md_m: float, t_s: float) -> float:
        """取 (md_m [m], t_s [s]) 处温度（°C）：双线性插值，越界 clamp+审计。"""
        md = float(md_m)
        tt = float(t_s)
        z = self.depth_m
        ts = self.time_s

        md_c = min(max(md, float(z[0])), float(z[-1]))
        tt_c = min(max(tt, float(ts[0])), float(ts[-1]))
        if md_c != md or tt_c != tt:
            self._oob_events.append(
                ClampEvent(
                    md_m=md,
                    t_s=tt,
                    md_clamped_m=md_c,
                    t_clamped_s=tt_c,
                )
            )

        i = int(np.searchsorted(z, md_c, side="right")) - 1
        i = min(max(i, 0), len(z) - 2)
        j = int(np.searchsorted(ts, tt_c, side="right")) - 1
        j = min(max(j, 0), len(ts) - 2)

        fz = (md_c - z[i]) / (z[i + 1] - z[i])
        ft = (tt_c - ts[j]) / (ts[j + 1] - ts[j])

        t00 = self.table[i, j]
        t01 = self.table[i, j + 1]
        t10 = self.table[i + 1, j]
        t11 = self.table[i + 1, j + 1]
        top = t00 + ft * (t01 - t00)
        bot = t10 + ft * (t11 - t10)
        return float(top + fz * (bot - top))

    # -- 构造（文件 + npz 缓存） --------------------------------------------
    @classmethod
    def from_files(
        cls,
        xlsx_path: Path = DEFAULT_T_IN_XLSX,
        depth_csv_path: Path = DEFAULT_DEPTH_CSV,
        *,
        cache_path: Optional[Path] = None,
        use_cache: bool = True,
        expected_shape: Tuple[int, int] = EXPECTED_SHAPE,
    ) -> "TableTemperatureField":
        """从交付 xlsx + 井身结构 CSV 加载（默认管内 T_in）。

        缓存：``cache_path`` 缺省时落在 xlsx 同目录（``*.npz``，已被 .gitignore
        的 ``*.npz`` 规则覆盖）。命中缓存则免解析 xlsx/CSV；缓存内记录源路径，
        与请求不符时视为失效、重新解析并覆盖。
        """
        xlsx_path = Path(xlsx_path)
        depth_csv_path = Path(depth_csv_path)
        if cache_path is None:
            cache_path = xlsx_path.with_suffix(".npz")
        cache_path = Path(cache_path)

        if use_cache and cache_path.exists():
            cached = cls._load_cache(cache_path)
            if cached is not None and cached[2] == str(xlsx_path):
                table, depth_m, _src, time_s = cached
                return cls(
                    table,
                    depth_m,
                    time_s,
                    expected_shape=expected_shape,
                    source=xlsx_path,
                )

        table = _read_xlsx_numeric(xlsx_path)
        depth_m = _read_depth_axis(depth_csv_path)
        field = cls(
            table,
            depth_m,
            expected_shape=expected_shape,
            source=xlsx_path,
        )
        if use_cache:
            cls._save_cache(cache_path, field, xlsx_path)
        return field

    @staticmethod
    def _load_cache(
        cache_path: Path,
    ) -> Optional[Tuple[np.ndarray, np.ndarray, str, np.ndarray]]:
        try:
            with np.load(cache_path, allow_pickle=False) as data:
                keys = set(data.files)
                if not {"table", "depth_m", "time_s", "source"} <= keys:
                    return None
                return (
                    np.asarray(data["table"], dtype=float),
                    np.asarray(data["depth_m"], dtype=float),
                    str(data["source"]),
                    np.asarray(data["time_s"], dtype=float),
                )
        except Exception:
            # 缓存损坏/不可读 → 视为未命中，回退解析源文件
            return None

    @staticmethod
    def _save_cache(
        cache_path: Path, field: "TableTemperatureField", xlsx_path: Path
    ) -> None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            cache_path,
            table=field.table,
            depth_m=field.depth_m,
            time_s=field.time_s,
            source=np.str_(str(xlsx_path)),
        )


# ---------------------------------------------------------------------------
# 交付件读取辅助
# ---------------------------------------------------------------------------
def _col_letters_to_index(letters: str) -> int:
    """Excel 列名（A/B/.../GV）→ 0 基列号。"""
    col = 0
    for ch in letters:
        col = col * 26 + (ord(ch) - ord("A") + 1)
    return col - 1


def _read_xlsx_numeric(xlsx_path: Path) -> np.ndarray:
    """读取 xlsx 首个工作表为二维 float 矩阵（纯数值、无表头；缺失单元格=NaN）。

    不依赖 openpyxl：直接解析 sheet XML 的 ``<c><v>`` 数值单元格。
    """
    xlsx_path = Path(xlsx_path)
    if not xlsx_path.exists():
        raise FileNotFoundError(f"温度场表文件不存在: {xlsx_path}")
    try:
        with zipfile.ZipFile(xlsx_path) as zf:
            xml_bytes = zf.read(_SHEET_XML)
    except KeyError as exc:
        raise ValueError(
            f"温度场表缺少工作表 {_SHEET_XML}: {xlsx_path}"
        ) from exc

    root = ET.fromstring(xml_bytes)
    ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    cells = {}
    max_row = max_col = -1
    for cell in root.iter(ns + "c"):
        ref = cell.get("r")
        if not ref:
            continue
        m = _CELL_REF_RE.match(ref)
        if m is None:
            continue
        col = _col_letters_to_index(m.group(1))
        row = int(m.group(2)) - 1
        cell_type = cell.get("t")
        if cell_type not in (None, "n"):
            # 交付表为纯数值；出现共享字符串/布尔等非数值单元格须显式报错，
            # 避免把共享字符串索引静默当温度读入
            raise ValueError(
                f"温度场表存在非数值单元格 {ref}（t={cell_type}）: {xlsx_path}"
            )
        v_node = cell.find(ns + "v")
        if v_node is None or v_node.text is None:
            continue
        cells[(row, col)] = float(v_node.text)
        max_row = max(max_row, row)
        max_col = max(max_col, col)

    if max_row < 0 or max_col < 0:
        raise ValueError(f"温度场表无有效数值单元格: {xlsx_path}")

    table = np.full((max_row + 1, max_col + 1), np.nan)
    for (r, c), val in cells.items():
        table[r, c] = val
    return table


def _read_depth_axis(depth_csv_path: Path) -> np.ndarray:
    """读井身结构 CSV 第 2 列 ``depth_well_logging_m_`` 作深度轴。"""
    depth_csv_path = Path(depth_csv_path)
    if not depth_csv_path.exists():
        raise FileNotFoundError(f"井身结构 CSV 不存在: {depth_csv_path}")

    try:
        series = pd.read_csv(depth_csv_path, usecols=[DEPTH_CSV_COLUMN])[
            DEPTH_CSV_COLUMN
        ]
    except ValueError as exc:
        raise ValueError(
            f"井身结构 CSV 缺少列 {DEPTH_CSV_COLUMN}: {depth_csv_path}"
        ) from exc
    return series.to_numpy(dtype=float)


# ---------------------------------------------------------------------------
# 交付两表加载（管内 + 环空）
# ---------------------------------------------------------------------------
def load_delivered_pair(
    *,
    t_in_xlsx: Path = DEFAULT_T_IN_XLSX,
    t_out_xlsx: Path = DEFAULT_T_OUT_XLSX,
    depth_csv_path: Path = DEFAULT_DEPTH_CSV,
    cache_dir: Optional[Path] = None,
    use_cache: bool = True,
    expected_shape: Tuple[int, int] = EXPECTED_SHAPE,
) -> Tuple[TableTemperatureField, TableTemperatureField]:
    """加载交付的管内（T_in）与环空（T_out）两表。

    校验井底（最深行）两表共享同值（浮点容差 1e-9）。
    ``cache_dir`` 给定时 npz 缓存放该目录，否则落各自 xlsx 同目录。
    """

    def _cache_for(xlsx: Path) -> Optional[Path]:
        if cache_dir is None:
            return None
        return Path(cache_dir) / (xlsx.stem + ".npz")

    t_in = TableTemperatureField.from_files(
        t_in_xlsx,
        depth_csv_path,
        cache_path=_cache_for(t_in_xlsx),
        use_cache=use_cache,
        expected_shape=expected_shape,
    )
    t_out = TableTemperatureField.from_files(
        t_out_xlsx,
        depth_csv_path,
        cache_path=_cache_for(t_out_xlsx),
        use_cache=use_cache,
        expected_shape=expected_shape,
    )
    if t_in.table.shape[0] != t_out.table.shape[0]:
        raise ValueError(
            f"管内/环空两表深度节点数不一致: {t_in.table.shape} vs {t_out.table.shape}"
        )
    bottom_in = t_in.table[-1]
    bottom_out = t_out.table[-1]
    if not np.allclose(bottom_in, bottom_out, atol=1e-9, rtol=0):
        max_diff = float(np.max(np.abs(bottom_in - bottom_out)))
        raise ValueError(
            f"井底（最深行）两表应共享同值，实际最大偏差 {max_diff:.3e} °C"
        )
    return t_in, t_out
