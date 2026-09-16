# -*- coding: utf-8 -*-
"""水泥屈服应力 τ_y 有源数据契约（Phase A Task 0，2026-09-16）。

本文件守护四件事，任一条被后人改动都会立刻变红：

1. **两路线登记完整性**：`cemdisp/data/cement_yield_stress.py` 的 21 条记录与
   8 个 loader 的**两套**模块级常数（`..._REPORT_GIVEN_PA` / `..._FITTED_PA`）
   一一对应；覆盖范围恰等于各 loader 实际消费的水泥相（不多取、不少取）；
   A 路线全井登记 `MISSING`、B 路线为有限非负值。
2. **出处可追溯**：每条记录的 `source_file` / `source_location` / 温度被断言成
   **字面字符串**，并区分"读数能否在 `rheometer_readings.csv` 直接查到"
   （`in_rheometer_csv`）——防后人改数据源或悄悄换出处。
3. **可复算**：用模块自带的 `fann_herschel_bulkley` 从**冻结读数**重算 B 路线 τ_y、
   `fann_bingham_two_point` 重算 API 口径旁枝、`fann_bingham_yield_stress` 重算
   口径旁证；与 loader 常数逐条比对。另钉死该口径的**现场锚**：
   `2061144.xls / 流变曲线1` 三组读数逐位复现。
4. **R7 护栏（本任务最容易被破坏处）**：水泥相 `FluidSpec` 的流变模型与屈服
   字段**必须与 HEAD 逐项一致**——7 井保持 `POWER_LAW` 且 `yield_stress_pa is None`，
   ht1_004 保持既有 `BINGHAM`（YP=13/14）。任何"顺手把新 τ_y 填进 FluidSpec"
   的改动都会在这里变红，因为它会让默认参数下的既有屈服门
   （`annulus_d2dga.py:857-864` 相体积加权 τ_y 场）立刻改变行为。

口径（与模块 docstring 一致）
-----------------------------
- Fann 35：γ [s⁻¹] = 1.703 × RPM、τ [Pa] = 0.511 × θ。
- 路线 A `report_given` = 化验报告原文直接给出的屈服应力（**本语料全部 MISSING**）。
- 路线 B `fitted_new` = HB 三参数 τ_y（τy≥0 约束）；另记 API 口径
  `0.511·(2θ₃₀ − θ₆₀)`（仅 ht1_003 可算）与 Bingham 线性最小二乘截距（口径旁证）。
- **两条路线并列、不合并、不取平均、不互相填补**；生产消费哪条由 controller 裁定。

Step 1「n/K 可复现性核实」的结论也钉在本文件里（见
`test_lab_printed_power_law_is_not_reproducible_by_a_single_convention`）：
**不存在能复现全部化验给定 n/K 的单一拟合口径**——语料内部本身就并存两套
（多数井 5 点对数-对数最小二乘；hu101 主检 2011122.pdf 用 300/100 两点式），
ht1_003 的给定 n/K 用两套都复现不出来。故 loader 只记出处与温度、不记拟合
方法，是**正确的历史取舍**；τ_y 只能作为 `fitted_new`（新推导值）登记，
不得冒充化验给定值。
"""
from __future__ import annotations

import importlib
import math
from pathlib import Path

import pytest

from cemdisp.data.cement_yield_stress import (
    CEMENT_YIELD_STRESS,
    FANN_SHEAR_RATE_PER_RPM,
    FANN_STRESS_PA_PER_DEGREE,
    MISSING,
    ROUTE_FITTED,
    ROUTE_REPORT_GIVEN,
    ROUTES,
    fann_bingham_two_point,
    fann_bingham_yield_stress,
    fann_herschel_bulkley,
    fann_power_law,
    get_cement_yield_stress,
    yield_stress_by_role,
    yield_stress_by_route,
)
from cemdisp.data.fluid_spec import FluidRole, RheologyModel

_CEMENT_ROLES = (FluidRole.LEAD, FluidRole.INTERMEDIATE, FluidRole.TAIL)

# --------------------------------------------------------------------------
# 冻结的 loader 入口表：(井键, 模块, 加载函数, 该井实际消费的水泥相)
# --------------------------------------------------------------------------
_LOADERS: tuple[tuple[str, str, str, tuple[str, ...]], ...] = (
    ("hu101", "cemdisp.data.loaders.hu101_loader", "load_hu101_tailpipe", ("lead", "tail")),
    ("hu102", "cemdisp.data.loaders.hu102_loader", "load_hu102_tailpipe", ("lead", "tail")),
    ("hu103", "cemdisp.data.loaders.hu103_loader", "load_hu103_tailpipe",
     ("lead", "intermediate", "tail")),
    ("hu1", "cemdisp.data.loaders.hu1_loader", "load_hu1_tailpipe", ("lead", "tail")),
    ("hu2", "cemdisp.data.loaders.hu2_loader", "load_hu2_tailpipe",
     ("lead", "intermediate", "tail")),
    ("ht1_001", "cemdisp.data.loaders.ht1_001_loader", "load_ht1_001_tailpipe",
     ("lead", "intermediate", "tail")),
    ("ht1_003", "cemdisp.data.loaders.ht1_003_loader", "load_ht1_003_tailpipe",
     ("lead", "tail")),
    ("ht1_004", "cemdisp.data.loaders.ht1_004_loader", "load_ht1_004_tailpipe",
     ("lead", "tail")),
)

# 常数名主干（命名契约 `<主干>_YIELD_STRESS_REPORT_GIVEN_PA` / `<主干>_YIELD_STRESS_FITTED_PA`）
_CONSTANT_STEMS: dict[str, dict[str, str]] = {
    "hu101": {
        "lead": "HU101_LEAD", "tail": "HU101_TAIL",
        "lead_recheck": "HU101_LEAD_RECHECK", "tail_recheck": "HU101_TAIL_RECHECK",
    },
    "hu102": {"lead": "HU102_LEAD", "tail": "HU102_TAIL"},
    "hu103": {"lead": "HU103_LEAD", "intermediate": "HU103_INTERMEDIATE",
              "tail": "HU103_TAIL"},
    "hu1": {"lead": "HU1_LEAD", "tail": "HU1_TAIL"},
    "hu2": {"lead": "HU2_LEAD", "intermediate": "HU2_INTERMEDIATE", "tail": "HU2_TAIL"},
    "ht1_001": {"lead": "HT1_001_LEAD", "intermediate": "HT1_001_INTERMEDIATE",
                "tail": "HT1_001_TAIL"},
    "ht1_003": {"lead": "HT1_003_LEAD", "tail": "HT1_003_TAIL"},
    "ht1_004": {"lead": "HT1_004_LEAD", "tail": "HT1_004_TAIL"},
}

# 冻结出处（source_file, source_location 字面值 + 是否可在 CSV 中直接查到）
_FROZEN_PROVENANCE: dict[tuple[str, str], tuple[str, str, bool]] = {
    ("hu101", "lead"): (
        "2/201/2011/20111/201112/2011122.pdf",
        "检测分析结果·领浆·流变性能（委托书 W301-22094）", False),
    ("hu101", "tail"): (
        "2/201/2011/20111/201112/2011122.pdf",
        "检测分析结果·尾浆·流变性能（委托书 W301-22094）", False),
    ("hu101", "lead_recheck"): (
        "2/201/2011/20111/201112/2011121.doc", "领浆流变", True),
    ("hu101", "tail_recheck"): (
        "2/201/2011/20111/201112/2011121.doc", "尾浆流变", True),
    ("hu102", "lead"): ("2021/20234.doc", "slurry rheology", True),
    ("hu102", "tail"): ("2021/20234.doc", "slurry rheology", True),
    ("hu103", "lead"): ("2\\203\\2031\\20311\\203111.docx", "Table7", True),
    ("hu103", "intermediate"): ("2\\203\\2031\\20311\\203111.docx", "Table7", True),
    ("hu103", "tail"): ("2\\203\\2031\\20311\\203111.docx", "Table7", True),
    ("hu1", "lead"): (
        "2/204/2041/20413/204131.doc",
        "油井水泥浆物理性能试验结果·第1个浆体块·流变性能", False),
    ("hu1", "tail"): (
        "2/204/2041/20413/204131.doc",
        "油井水泥浆物理性能试验结果·第2个浆体块·流变性能", False),
    ("hu2", "lead"): ("化验报告 Table7", "化验报告 Table7", True),
    ("hu2", "intermediate"): ("化验报告 Table7", "化验报告 Table7", True),
    ("hu2", "tail"): ("化验报告 Table7", "化验报告 Table7", True),
    ("ht1_001", "lead"): ("化验报告/Table 8", "128->93C降温测试", True),
    ("ht1_001", "intermediate"): ("化验报告/Table 8", "128->93C降温测试", True),
    ("ht1_001", "tail"): ("化验报告/Table 8", "128->93C降温测试", True),
    ("ht1_003", "lead"): (
        "化验报告/HT1-003 油层尾管 化验报告.docx 表7", "表7 流变性能 · 领浆", True),
    ("ht1_003", "tail"): (
        "化验报告/HT1-003 油层尾管 化验报告.docx 表7", "表7 流变性能 · 尾浆", True),
    ("ht1_004", "lead"): (
        "化验报告/HT1-004 油层尾管 化验报告.docx 表7", "表7 流变性能 · 领浆", True),
    ("ht1_004", "tail"): (
        "化验报告/HT1-004 油层尾管 化验报告.docx 表7", "表7 流变性能 · 尾浆", True),
}

# 冻结的 B 路线值（4 位小数；HB 三参数 τy，τy≥0 约束；0.0 = 落在 τy=0 边界）
_FROZEN_FITTED_TAU_Y_PA: dict[tuple[str, str], float] = {
    ("hu101", "lead"): 0.9697, ("hu101", "tail"): 1.5111,
    ("hu101", "lead_recheck"): 1.5182, ("hu101", "tail_recheck"): 1.2804,
    ("hu102", "lead"): 0.0000, ("hu102", "tail"): 0.0000,
    ("hu103", "lead"): 0.0000, ("hu103", "intermediate"): 0.4961,
    ("hu103", "tail"): 0.8967,
    ("hu1", "lead"): 0.0000, ("hu1", "tail"): 1.2067,
    ("hu2", "lead"): 0.0000, ("hu2", "intermediate"): 0.0000, ("hu2", "tail"): 0.0000,
    ("ht1_001", "lead"): 0.0000, ("ht1_001", "intermediate"): 0.0000,
    ("ht1_001", "tail"): 0.0000,
    ("ht1_003", "lead"): 4.0132, ("ht1_003", "tail"): 3.8722,
    ("ht1_004", "lead"): 0.0000, ("ht1_004", "tail"): 0.0000,
}

_FROZEN_REPORT_GIVEN_TAU_Y_PA: dict[tuple[str, str], object] = {
    key: MISSING for key in _FROZEN_FITTED_TAU_Y_PA
}
"""路线 A（报告原文直接给出）：**8 井全部 MISSING**——化验报告只给 n/K。

不是"没查"：逐井取证说明见 `CementYieldStress.report_given_note`
（含 hu101 2011122.pdf 扫描件无文字层的取证受限说明）。**不许跨路线借值填补**。
"""

# 路线 B 旁枝：API 口径 0.511·(2θ₃₀−θ₆₀)，仅 ht1_003 两相可算
_FROZEN_BINGHAM_TWO_POINT_PA: dict[tuple[str, str], object] = {
    key: MISSING for key in _FROZEN_FITTED_TAU_Y_PA
}
_FROZEN_BINGHAM_TWO_POINT_PA[("ht1_003", "lead")] = 17.885
_FROZEN_BINGHAM_TWO_POINT_PA[("ht1_003", "tail")] = 20.440

# 冻结的 (RPM, θ) 读数——直读，改动即红
_FROZEN_READINGS: dict[tuple[str, str], tuple[tuple[int, float], ...]] = {
    ("hu101", "lead"): ((300, 144), (200, 104), (100, 57), (6, 8), (3, 4)),
    ("hu101", "tail"): ((300, 122), (200, 82), (100, 49), (6, 7), (3, 4)),
    ("hu101", "lead_recheck"): ((300, 149), (200, 106), (100, 60), (6, 9), (3, 5)),
    ("hu101", "tail_recheck"): ((300, 126), (200, 87), (100, 53), (6, 8), (3, 4)),
    ("hu102", "lead"): ((300, 163), (200, 128), (100, 106), (6, 10), (3, 6)),
    ("hu102", "tail"): ((300, 163), (200, 128), (100, 106), (6, 10), (3, 6)),
    ("hu103", "lead"): ((300, 224), (200, 158), (100, 93), (6, 9), (3, 5)),
    ("hu103", "intermediate"): ((300, 254), (200, 179), (100, 112), (6, 14), (3, 7)),
    ("hu103", "tail"): ((300, 265), (200, 186), (100, 115), (6, 15), (3, 7)),
    ("hu1", "lead"): ((300, 171), (200, 132), (100, 81), (6, 10), (3, 6)),
    ("hu1", "tail"): ((300, 117), (200, 86), (100, 52), (6, 9), (3, 5)),
    ("hu2", "lead"): ((300, 268), (200, 194), (100, 113), (6, 10), (3, 7)),
    ("hu2", "intermediate"): ((300, 221), (200, 165), (100, 92), (6, 7), (3, 4)),
    ("hu2", "tail"): ((300, 218), (200, 163), (100, 88), (6, 6), (3, 4)),
    ("ht1_001", "lead"): ((300, 268), (200, 194), (100, 113), (6, 10), (3, 7)),
    ("ht1_001", "intermediate"): ((300, 221), (200, 165), (100, 92), (6, 7), (3, 4)),
    ("ht1_001", "tail"): ((300, 218), (200, 163), (100, 88), (6, 6), (3, 4)),
    ("ht1_003", "lead"): ((600, 207), (300, 121), (200, 93), (100, 60), (6, 15), (3, 9)),
    ("ht1_003", "tail"): ((600, 196), (300, 118), (200, 89), (100, 56), (6, 14), (3, 10)),
    ("ht1_004", "lead"): ((300, 296), (200, 213), (100, 120), (6, 10), (3, 6)),
    ("ht1_004", "tail"): ((300, 293), (200, 211), (100, 120), (6, 10), (3, 5)),
}

_ALL_KEYS = sorted(_FROZEN_FITTED_TAU_Y_PA)


def _loader_module(well_key: str):
    return importlib.import_module(next(m for k, m, _, _ in _LOADERS if k == well_key))


def _load_cement_fluids(well_key: str):
    """返回该井 loader 实际构造的水泥相 FluidSpec 列表。"""
    for key, module_name, function_name, _ in _LOADERS:
        if key != well_key:
            continue
        module = importlib.import_module(module_name)
        _, fluids, _, _ = getattr(module, function_name)()
        return [f for f in fluids if f.role in _CEMENT_ROLES]
    raise AssertionError(f"未登记的井键：{well_key!r}")


# --------------------------------------------------------------------------
# 1. 两路线登记完整性
# --------------------------------------------------------------------------
def test_every_loader_constant_exists_for_both_routes():
    """两套常数均存在：A 路线必须登记 MISSING 哨兵，B 路线为有限非负 float。"""
    for well_key, phases in _CONSTANT_STEMS.items():
        module = _loader_module(well_key)
        for phase, stem in phases.items():
            given = getattr(module, f"{stem}_YIELD_STRESS_REPORT_GIVEN_PA", None)
            fitted = getattr(module, f"{stem}_YIELD_STRESS_FITTED_PA", None)
            assert given == MISSING, (
                f"{well_key}/{phase} A 路线应为显式 MISSING，得到 {given!r}")
            assert isinstance(fitted, float), f"{well_key}/{phase} B 路线必须是 float"
            assert math.isfinite(fitted), f"{well_key}/{phase} B 路线必须有限"
            assert fitted >= 0.0, f"{well_key}/{phase} B 路线 τy 不可为负"


def test_loader_missing_sentinel_matches_registry():
    """各 loader 的缺项哨兵字面值与汇总模块一致（防两处哨兵漂移）。"""
    for well_key in _CONSTANT_STEMS:
        module = _loader_module(well_key)
        assert getattr(module, "MISSING") == MISSING


def test_registry_covers_exactly_the_phases_each_loader_consumes():
    """登记范围恰等于各 loader 实际消费的水泥相——不多取、不少取。"""
    assert sorted(CEMENT_YIELD_STRESS) == sorted(t for t, _, _, _ in _LOADERS)
    for well_key, _, _, phases in _LOADERS:
        registered = set(CEMENT_YIELD_STRESS[well_key])
        expected = set(phases)
        if well_key == "hu101":
            expected = expected | {"lead_recheck", "tail_recheck"}
        assert registered == expected, f"{well_key} 登记相不匹配：{registered} vs {expected}"


def test_fluid_name_matches_the_loader_fluid_spec():
    """记录里的 fluid_name 必须就是该井 loader 构造的流体名（防井相错配）。"""
    for well_key, _, _, phases in _LOADERS:
        fluids = _load_cement_fluids(well_key)
        by_role = {f.role.value: f.name for f in fluids}
        for phase in phases:
            record = CEMENT_YIELD_STRESS[well_key][phase]
            assert record.fluid_name == by_role[phase], (
                f"{well_key}/{phase} 流体名不符：记录={record.fluid_name!r} loader={by_role[phase]!r}")


def test_no_silent_none_registration():
    """不许静默 None：按路线取值必须返回值 / MISSING / 抛错。"""
    assert MISSING == "MISSING"
    assert ROUTES == (ROUTE_REPORT_GIVEN, ROUTE_FITTED)
    for well_key, phase in _ALL_KEYS:
        record = get_cement_yield_stress(well_key, phase)
        assert record.report_given_pa == MISSING
        assert record.fitted_hb_tau_y_pa is not None
        assert yield_stress_by_route(well_key, phase, ROUTE_REPORT_GIVEN) == MISSING
        assert (yield_stress_by_route(well_key, phase, ROUTE_FITTED)
                == record.fitted_hb_tau_y_pa)
    with pytest.raises(KeyError):
        get_cement_yield_stress("hu101", "不存在相")
    with pytest.raises(KeyError):
        get_cement_yield_stress("不存在的井", "lead")
    with pytest.raises(ValueError):
        yield_stress_by_route("hu101", "lead", "不存在的路线")


def test_by_role_accessor_matches_records():
    """Task 6 的按角色取用接口按路线逐位一致，且不含复检口径。"""
    for well_key in CEMENT_YIELD_STRESS:
        for route in ROUTES:
            mapping = yield_stress_by_role(well_key, route)
            assert set(mapping) <= {"LEAD", "INTERMEDIATE", "TAIL"}
            for role_name, value in mapping.items():
                record = CEMENT_YIELD_STRESS[well_key][role_name.lower()]
                expected = yield_stress_by_route(well_key, role_name.lower(), route)
                assert value == expected


# --------------------------------------------------------------------------
# 2. 出处可追溯（字面字符串 + CSV 直查标记）
# --------------------------------------------------------------------------
def test_provenance_strings_are_frozen():
    """每条记录的出处字符串本身被断言——防后人改动数据源或悄悄换出处。"""
    for well_key, phase in _ALL_KEYS:
        record = CEMENT_YIELD_STRESS[well_key][phase]
        expected_file, expected_location, expected_in_csv = _FROZEN_PROVENANCE[(well_key, phase)]
        assert record.source_file == expected_file, f"{well_key}/{phase} source_file 被改动"
        assert record.source_location == expected_location, (
            f"{well_key}/{phase} source_location 被改动")
        assert record.in_rheometer_csv is expected_in_csv, (
            f"{well_key}/{phase} in_rheometer_csv 被改动")
        assert record.temperature_c, f"{well_key}/{phase} 缺温度登记"
        assert record.report_given_note, f"{well_key}/{phase} 缺 A 路线取证说明"


def test_temperature_is_registered_for_every_record():
    """温度必须逐条登记（本任务硬约束：出处 = 文件 + 小节 + 温度）。"""
    frozen_temperatures = {
        ("hu101", "lead"): "89", ("hu101", "tail"): "89",
        ("hu101", "lead_recheck"): "93", ("hu101", "tail_recheck"): "93",
        ("hu102", "lead"): "93", ("hu102", "tail"): "93",
        ("hu103", "lead"): "140→93", ("hu103", "intermediate"): "140→93",
        ("hu103", "tail"): "133→93",
        ("hu1", "lead"): "93", ("hu1", "tail"): "93",
        ("hu2", "lead"): "133→93", ("hu2", "intermediate"): "133→93",
        ("hu2", "tail"): "133→93",
        ("ht1_001", "lead"): "128→93", ("ht1_001", "intermediate"): "128→93",
        ("ht1_001", "tail"): "128→93",
        ("ht1_003", "lead"): "129→93", ("ht1_003", "tail"): "129→93",
        ("ht1_004", "lead"): "132→93", ("ht1_004", "tail"): "132→93",
    }
    assert sorted(frozen_temperatures) == _ALL_KEYS
    for key, expected in frozen_temperatures.items():
        assert CEMENT_YIELD_STRESS[key[0]][key[1]].temperature_c == expected


def test_out_of_range_theta600_is_excluded_not_read_as_300():
    """两处数据陷阱的护栏：θ₆₀₀ 超量程 / 未记录一律剔除，不得当 300 用。

    - ht1_001 化验报告原文记 `>300`，提取包 CSV 却写作 300；
    - hu2 记 `>300`；hu103 记 `/`；hu101/ht1_004 记 `>300` 或未测。
    因此除 ht1_003（真实读数 207/196）外，**任何记录都不得含 600 rpm 点**。
    """
    for well_key, phase in _ALL_KEYS:
        readings = CEMENT_YIELD_STRESS[well_key][phase].readings
        speeds = {rpm for rpm, _ in readings}
        if well_key == "ht1_003":
            assert 600 in speeds, "ht1_003 的 θ600 是真实读数，应参与拟合"
        else:
            assert 600 not in speeds, (
                f"{well_key}/{phase} 含 600 rpm 点：超量程/未记录的 θ600 不得当读数用")
        assert 300 in speeds and len(readings) >= 5


# --------------------------------------------------------------------------
# 3. 可复算（从冻结读数重算）
# --------------------------------------------------------------------------
def test_frozen_readings_match_registry():
    """读数本身钉死：登记读数与冻结表逐位一致。"""
    for key, expected in _FROZEN_READINGS.items():
        assert CEMENT_YIELD_STRESS[key[0]][key[1]].readings == expected, (
            f"{key} 冻结读数被改动")


def test_fitted_tau_y_is_recomputed_from_frozen_readings():
    """B 路线值必须能由冻结读数用模块自带函数复算（容差 = 4 位取整半量）。"""
    for well_key, phase in _ALL_KEYS:
        record = CEMENT_YIELD_STRESS[well_key][phase]
        hb_tau_y, hb_k, hb_n = fann_herschel_bulkley(record.readings)
        assert hb_tau_y == pytest.approx(record.fitted_hb_tau_y_pa, abs=5e-5)
        assert hb_k == pytest.approx(record.fitted_hb_consistency_k_pa_s_n, abs=5e-5)
        assert hb_n == pytest.approx(record.fitted_hb_power_law_n, abs=5e-5)
        assert record.fitted_bingham_ls_r_squared > 0.90, (
            f"{well_key}/{phase} Bingham 口径旁证的 R² 过低")


def test_loader_constants_equal_registry_for_both_routes():
    """loader 两套常数与汇总登记逐位一致（单一真值来源）。"""
    for well_key, phase in _ALL_KEYS:
        module = _loader_module(well_key)
        stem = _CONSTANT_STEMS[well_key][phase]
        assert getattr(module, f"{stem}_YIELD_STRESS_REPORT_GIVEN_PA") == MISSING
        fitted = getattr(module, f"{stem}_YIELD_STRESS_FITTED_PA")
        assert fitted == pytest.approx(_FROZEN_FITTED_TAU_Y_PA[(well_key, phase)], abs=1e-12)
        # 记录里存的是全精度复算值，loader 常数按 4 位小数取整 → 容差取取整半量
        assert fitted == pytest.approx(
            CEMENT_YIELD_STRESS[well_key][phase].fitted_hb_tau_y_pa, abs=5e-5)
    assert sorted(_FROZEN_REPORT_GIVEN_TAU_Y_PA) == _ALL_KEYS
    assert all(v == MISSING for v in _FROZEN_REPORT_GIVEN_TAU_Y_PA.values())


def test_herschel_bulkley_boundary_hits_are_explicitly_registered():
    """HB 不可辨识性的显式登记：21 条里 12 条落在 τy = 0 边界。

    边界命中是**取值**不是缺失；数量变化（新井、改读数、改方法）一律变红，
    逼迫复核。这也是 HB 不宜直接作为唯一生产值的方法学代价。
    """
    at_boundary = []
    for well_key, phase in _ALL_KEYS:
        record = CEMENT_YIELD_STRESS[well_key][phase]
        upper = FANN_STRESS_PA_PER_DEGREE * min(theta for _, theta in record.readings) + 1e-9
        assert 0.0 <= record.fitted_hb_tau_y_pa <= upper, (
            f"{well_key}/{phase} HB τy 超出 [0, min τ] 物理可行域")
        if record.fitted_hb_tau_y_pa <= 1e-6:
            at_boundary.append((well_key, phase))
    assert len(at_boundary) == 12, f"HB 边界命中数变化：{at_boundary}"
    assert ("hu102", "lead") in at_boundary and ("ht1_004", "tail") in at_boundary
    assert ("ht1_003", "lead") not in at_boundary


def test_bingham_two_point_is_computed_only_where_theta600_is_a_real_reading():
    """API 口径 `0.511·(2θ₃₀−θ₆)` 只在 θ₀₀ 是真实读数时给出；其余登记 MISSING。

    本语料只有 ht1_003 两相有真实 θ₆₀（207/196）；其余井 θ₆₀ 是 `>300` 超量程
    或 `/` 未记录，**不许拿 300 顶上**（那会把 YP 系统性算小）。
    """
    for well_key, phase in _ALL_KEYS:
        record = CEMENT_YIELD_STRESS[well_key][phase]
        expected = _FROZEN_BINGHAM_TWO_POINT_PA[(well_key, phase)]
        if expected == MISSING:
            assert record.fitted_bingham_two_point_pa == MISSING
        else:
            assert record.fitted_bingham_two_point_pa == pytest.approx(expected, abs=5e-4)
            assert fann_bingham_two_point(record.readings) == pytest.approx(expected, abs=5e-4)


def test_two_routes_are_side_by_side_and_never_merged():
    """两路线并列呈现：A 全 MISSING、B 有值；三条 B 口径各自独立保留、未合并。

    本测试只固定"并列不合并"这一契约，不对"哪条更权威"表态——那属 controller 裁定。
    """
    for well_key, phase in _ALL_KEYS:
        record = CEMENT_YIELD_STRESS[well_key][phase]
        assert record.report_given_pa == MISSING, "A 路线被填了值：必须先有报告出处"
        # B 路线三条口径（HB τy / Bingham 截距 / API 两点式）都保留在原位
        assert hasattr(record, "fitted_hb_tau_y_pa")
        assert math.isfinite(record.fitted_bingham_ls_intercept_pa)
        assert record.fitted_bingham_two_point_pa is not None


# --------------------------------------------------------------------------
# 3b. 口径锚：现场固井工程计算表逐位复现（本口径不是本项目发明的证据）
# --------------------------------------------------------------------------
_TOOLKIT_CASES = (
    # (读数, 表内"屈服值，Pa", 表内"塑性粘度，mPa.s", 表内 R², 表内 n, 表内 k)
    # 领浆：表内 600 rpm 行为空，回归输入集为 300/200/100/6/3。
    (((300, 289), (200, 232), (100, 135), (6, 11), (3, 8)),
     8.062209467911948, 293.70940363180233, 0.9764922350398684,
     0.8160878346208227, 0.9744824295628708),
    # 尾浆：表内 600 rpm 行填了 48，但其回归区间**不含**该行（用 300/200/100/6/3），
    # 故此处按工作簿实际输入集给出——含 600 会算成 τy=1.520、n=0.572，与表内不符。
    (((300, 28), (200, 19), (100, 11), (6, 3), (3, 2)),
     1.0993162438058675, 25.74074224368921, 0.9990531341085161,
     0.5427129654442487, 0.4169271535084223),
    (((600, 70), (300, 38), (200, 27), (100, 15), (6, 4), (3, 2)),
     1.6832133532736824, 33.81213606910683, 0.9982058892000196,
     0.6221678228612535, 0.3975288752750177),
)
"""`0708/2/206/2061/20611/206114/2061144.xls` 的 `流变曲线1` 表两组读数。

该工作簿是现场固井工程计算工具箱，对领浆/尾浆/钻井液分别给出"屈服值，Pa /
塑性粘度，相关系数"，以及 PowerLaw 的 n / k / 相关系数。上列为表内原值
（xlrd 直读，未经手工誊抄）。**注意**：该工作簿不含井号/井段标识，其读数与
8 口井任一浆体都不对应，故其数值**不可归属、不作为路线 A 的 report_given 值**；
此处只钉死**计算口径**。
"""



def test_toolkit_yield_stress_convention_is_bitwise_reproduced():
    """现场工程计算表的"屈服值" = Bingham 线性最小二乘截距；读数逐位复现。

    这是 B 路线旁证口径的锚：证明它不是本项目发明，而是现场既有工具的口径。
    若后人把该旁证口径改成别的估计量，此测试立刻变红。
    """
    for readings, tau_y, mu_p_mpa_s, r_squared, n, toolkit_k in _TOOLKIT_CASES:
        got_tau_y, got_mu_p, got_r2 = fann_bingham_yield_stress(readings)
        assert got_tau_y == pytest.approx(tau_y, rel=1e-9)
        assert got_mu_p * 1000.0 == pytest.approx(mu_p_mpa_s, rel=1e-9)
        assert got_r2 == pytest.approx(r_squared, rel=1e-9)
        got_n, got_k, got_pl_r2 = fann_power_law(readings)
        assert got_n == pytest.approx(n, rel=1e-9)
        assert got_pl_r2 > 0.98
        assert got_k == pytest.approx(toolkit_k, rel=1e-9)


def test_fann_conversion_constants_are_frozen_and_no_extra_scaling():
    """换算常量钉死，且 power-law K 与现场工具箱三列**逐位一致**（无额外 0.511）。

    `fann_power_law` 的 K = 0.511·exp(b)（SI，Pa·sⁿ）。曾一度误判"工具箱领浆列
    漏乘 0.511"，实为手工换算口误：该表三列（领浆/尾浆/钻井液）与本函数逐位一致。
    此测试把"三列全一致"钉死，防后人按错误结论改口径。
    """
    for readings, _tau_y, _mu, _r2, _n, toolkit_k in _TOOLKIT_CASES:
        got_n, got_k, _ = fann_power_law(readings)
        assert got_k == pytest.approx(toolkit_k, rel=1e-9)
        assert got_k == pytest.approx(
            FANN_STRESS_PA_PER_DEGREE * math.exp(
                np_polyfit_intercept(readings)), rel=1e-9)


def np_polyfit_intercept(readings):
    """返回 (ln γ̇, ln θ) 回归的截距（仅供上面的口径断言用）。"""
    import numpy as np
    gamma = np.array([FANN_SHEAR_RATE_PER_RPM * rpm for rpm, _ in readings], dtype=float)
    theta = np.array([value for _, value in readings], dtype=float)
    return float(np.polyfit(np.log(gamma), np.log(theta), 1)[1])


def test_fann_conversion_constants_are_frozen():
    """换算常量钉死：0.511 与 1.703 不得被换成另一套（0.4788）或混用。"""
    assert FANN_STRESS_PA_PER_DEGREE == 0.511
    assert FANN_SHEAR_RATE_PER_RPM == 1.703


# --------------------------------------------------------------------------
# 3c. Step 1「n/K 可复现性核实」结论
# --------------------------------------------------------------------------
_LAB_PRINTED_NK = (
    # (井, 相, 化验给定 n, K, 报告口径)
    ("hu101", "lead", 0.844, 0.381, "2011122.pdf 89℃（主检）"),
    ("hu101", "tail", 0.830, 0.352, "2011122.pdf 89℃（主检）"),
    ("hu101", "lead_recheck", 0.719, 0.815, "2011121.doc 93℃（复检）"),
    ("hu101", "tail_recheck", 0.722, 0.684, "2011121.doc 93℃（复检）"),
    ("hu102", "lead", 0.737, 0.947, "20234.doc 93℃"),
    ("hu103", "lead", 0.82, 0.67, "203111.docx 表7"),
    ("hu103", "intermediate", 0.76, 1.11, "203111.docx 表7"),
    ("hu103", "tail", 0.76, 1.14, "203111.docx 表7"),
    ("hu1", "lead", 0.732, 0.933, "204131.doc 93℃"),
    ("hu1", "tail", 0.666, 0.906, "204131.doc 93℃"),
    ("hu2", "lead", 0.811, 0.876, "化验报告 表8"),
    ("hu2", "intermediate", 0.871, 0.504, "化验报告 表8"),
    ("hu2", "tail", 0.886, 0.453, "化验报告 表8"),
    ("ht1_001", "lead", 0.811, 0.876, "化验报告 表7"),
    ("ht1_001", "intermediate", 0.871, 0.504, "化验报告 表7"),
    ("ht1_001", "tail", 0.886, 0.453, "化验报告 表7"),
    ("ht1_003", "lead", 0.597, 1.622, "化验报告 表7"),
    ("ht1_003", "tail", 0.585, 1.673, "化验报告 表7"),
    ("ht1_004", "lead", 0.853, 0.746, "化验报告 表7"),
    ("ht1_004", "tail", 0.869, 0.669, "化验报告 表7"),
)

_NK_LS_TOLERANCE = 0.02
_NK_PAIR_300_100_TOLERANCE = 0.002
_NK_LS_NON_REPRODUCIBLE = {
    # 5 点对数-对数最小二乘复现不出的：hu101 主检两相（它们用 300/100 两点式）
    # 与 ht1_003 两相（两套口径都复现不出，属源报告内部不自洽）
    ("hu101", "lead"), ("hu101", "tail"),
    ("ht1_003", "lead"), ("ht1_003", "tail"),
}
_NK_PAIR_300_100_NON_REPRODUCIBLE = {("ht1_003", "lead"), ("ht1_003", "tail")}
_NK_PAIR_300_100_REPRODUCIBLE = {("hu101", "lead"), ("hu101", "tail")}


def test_lab_printed_power_law_is_not_reproducible_by_a_single_convention():
    """Step 1 结论：**不存在**能复现全部化验给定 n/K 的单一拟合口径。

    钉死三点，供后人复核（这也是 Task 0 把 τ_y 登记为 B 路线新推导值的依据）：
    1. 5 点对数-对数最小二乘复现 20 条里的 **16 条**（|Δn| ≤ 0.02、|ΔK| ≤ 0.05）；
    2. hu101 主检 2011122.pdf 的两相只能由 **300/100 两点式** 复现（附逐位验算），
       该两行同时被 5 点 LS 判为不可复现——即**同一语料内部并存两套口径**；
    3. ht1_003 两相的给定 n/K 用**上述两套口径都**复现不出来 —— 属源报告内部
       不自洽，不是本项目的方法问题。
    """
    ls_reproducible, ls_failed = set(), set()
    for well_key, phase, printed_n, printed_k, _label in _LAB_PRINTED_NK:
        readings = CEMENT_YIELD_STRESS[well_key][phase].readings
        n, k, _ = fann_power_law(readings)
        if abs(n - printed_n) <= _NK_LS_TOLERANCE:
            ls_reproducible.add((well_key, phase))
            assert k == pytest.approx(printed_k, abs=0.05), f"{well_key}/{phase} K 复现超差"
        else:
            ls_failed.add((well_key, phase))
    assert ls_failed == _NK_LS_NON_REPRODUCIBLE, f"不可复现集合变化：{ls_failed}"
    assert len(ls_reproducible) == 16, f"可复现条数变化：{sorted(ls_reproducible)}"
    assert len(ls_reproducible) + len(ls_failed) == len(_LAB_PRINTED_NK)

    # hu101 主检：300/100 两点式复现（θ₃₀=144、θ₁₀=57 → n=0.844；尾浆 122/49 → 0.830）
    for well_key, phase, printed_n, printed_k, _ in (
            c for c in _LAB_PRINTED_NK if (c[0], c[1]) in _NK_PAIR_300_100_REPRODUCIBLE):
        readings = dict(CEMENT_YIELD_STRESS[well_key][phase].readings)
        n_pair = math.log(readings[300] / readings[100]) / math.log(300.0 / 100.0)
        assert n_pair == pytest.approx(printed_n, abs=0.001), (
            f"{well_key}/{phase} 300/100 两点式未复现给定 n")
        k_pair = (FANN_STRESS_PA_PER_DEGREE * readings[100]
                  / (FANN_SHEAR_RATE_PER_RPM * 100.0) ** n_pair)
        assert k_pair == pytest.approx(printed_k, abs=0.002), (
            f"{well_key}/{phase} 300/100 两点式未复现给定 K")

    # ht1_003：两套口径都复现不了（源报告内部不自洽）
    for phase in ("lead", "tail"):
        readings = CEMENT_YIELD_STRESS["ht1_003"][phase].readings
        by_rpm = dict(readings)
        printed_n = next(n for w, p, n, _, _ in _LAB_PRINTED_NK
                         if (w, p) == ("ht1_003", phase))
        n_pair = math.log(by_rpm[300] / by_rpm[100]) / math.log(300.0 / 100.0)
        n_ls, _, _ = fann_power_law(readings)
        assert abs(n_pair - printed_n) > 0.02
        assert abs(n_ls - printed_n) > 0.02


# --------------------------------------------------------------------------
# 4. R7 护栏：水泥 FluidSpec 与 HEAD 逐项一致
# --------------------------------------------------------------------------
_HEAD_CEMENT_RHEOLOGY: dict[str, dict[str, tuple[RheologyModel, float, float]]] = {
    # 井键 → {相: (模型, 幂律 n 或 Bingham PV, K 或 Bingham YP)}
    "hu101": {"lead": (RheologyModel.POWER_LAW, 0.844, 0.381),
              "tail": (RheologyModel.POWER_LAW, 0.830, 0.352)},
    "hu102": {"lead": (RheologyModel.POWER_LAW, 0.737, 0.947),
              "tail": (RheologyModel.POWER_LAW, 0.737, 0.947)},
    "hu103": {"lead": (RheologyModel.POWER_LAW, 0.82, 0.67),
              "intermediate": (RheologyModel.POWER_LAW, 0.76, 1.11),
              "tail": (RheologyModel.POWER_LAW, 0.76, 1.14)},
    "hu1": {"lead": (RheologyModel.POWER_LAW, 0.732, 0.933),
            "tail": (RheologyModel.POWER_LAW, 0.666, 0.906)},
    "hu2": {"lead": (RheologyModel.POWER_LAW, 0.811, 0.876),
            "intermediate": (RheologyModel.POWER_LAW, 0.871, 0.504),
            "tail": (RheologyModel.POWER_LAW, 0.886, 0.453)},
    "ht1_001": {"lead": (RheologyModel.POWER_LAW, 0.811, 0.876),
                "intermediate": (RheologyModel.POWER_LAW, 0.871, 0.504),
                "tail": (RheologyModel.POWER_LAW, 0.886, 0.453)},
    "ht1_003": {"lead": (RheologyModel.POWER_LAW, 0.597, 1.622),
                "tail": (RheologyModel.POWER_LAW, 0.585, 1.673)},
    # ht1_004 生产口径的水泥相本来就是 BINGHAM（优化参数化），不是 POWER_LAW。
    "ht1_004": {"lead": (RheologyModel.BINGHAM, 0.17, 13.0),
                "tail": (RheologyModel.BINGHAM, 0.18, 14.0)},
}


def test_R7_cement_fluid_spec_rheology_is_unchanged_from_head():
    """R7 最高优先护栏：新增 τ_y 常数**不得**改动任何水泥相的 FluidSpec。

    本任务最容易被后人破坏的方式，就是"顺手"把新 τ_y 填进 `FluidSpec.yield_stress_pa`
    或把 `rheology_model` 改成 HERSCHEL_BULKLEY——那会让 `_fluid_yield_stress()`
    (`annulus_d2dga.py:706-711`) 立刻返回非零 τ_y，经 `annulus_d2dga.py:857-864`
    的相体积加权污染**既有屈服门**，破坏 Phase A「默认关逐位 = HEAD」的硬约束。
    """
    assert sorted(_HEAD_CEMENT_RHEOLOGY) == sorted(t for t, _, _, _ in _LOADERS)
    for well_key, expected in _HEAD_CEMENT_RHEOLOGY.items():
        fluids = _load_cement_fluids(well_key)
        by_role = {f.role.value: f for f in fluids}
        assert sorted(by_role) == sorted(expected), f"{well_key} 水泥相集合变化"
        for phase, (model, first, second) in expected.items():
            fluid = by_role[phase]
            assert fluid.rheology_model == model, (
                f"{well_key}/{phase} 流变模型被改动：{fluid.rheology_model} != {model}")
            if model is RheologyModel.POWER_LAW:
                assert fluid.power_law_n == pytest.approx(first, abs=1e-9)
                assert fluid.consistency_k == pytest.approx(second, abs=1e-9)
                assert fluid.yield_stress_pa is None, (
                    f"{well_key}/{phase} 幂律水泥相被填了 yield_stress_pa，违反 R7")
                assert fluid.plastic_viscosity_pa_s is None
            else:
                assert fluid.plastic_viscosity_pa_s == pytest.approx(first, abs=1e-9)
                assert fluid.yield_stress_pa == pytest.approx(second, abs=1e-9)


def test_R7_new_constants_are_not_consumed_by_any_fluid_spec():
    """新增 τ_y 与 FluidSpec 的屈服字段**必须不同源**：值不得等于任何水泥相 τ_y。"""
    for well_key in _HEAD_CEMENT_RHEOLOGY:
        for fluid in _load_cement_fluids(well_key):
            if fluid.yield_stress_pa is None:
                continue
            for phase in CEMENT_YIELD_STRESS[well_key]:
                for route in ROUTES:
                    value = yield_stress_by_route(well_key, phase, route)
                    if value == MISSING:
                        continue
                    assert fluid.yield_stress_pa != pytest.approx(value, abs=1e-9), (
                        f"{well_key}/{fluid.name} 的 yield_stress_pa 被新常数顶替，违反 R7")


def test_ht1_004_cement_stays_bingham():
    """ht1_004 是本任务唯一"水泥相本来就有非零 τ_y"的井，必须原样保留。

    该井生产口径（`load_ht1_004_tailpipe` → `_build_fluids`）的水泥相是 Bingham
    优化参数化 PV0.17/YP13.0、PV0.18/YP14.0（来源『优化参数.docx 2026-06-11』），
    与本次新增的 B 路线值（4.0132 量级来自化验表7 读数新推导）**不同源**，
    两者不得互相替代。本测试同时锁住 loader 的元数据声明，防止后人误当化验
    实测引用。
    """
    module = importlib.import_module("cemdisp.data.loaders.ht1_004_loader")
    _, fluids, _, validation = module.load_ht1_004_tailpipe()
    by_role = {f.role.value: f for f in fluids if f.role in _CEMENT_ROLES}
    assert by_role["lead"].rheology_model is RheologyModel.BINGHAM
    assert by_role["lead"].yield_stress_pa == pytest.approx(13.0)
    assert by_role["tail"].yield_stress_pa == pytest.approx(14.0)
    # 新推导值只作新常数存在，不进入 FluidSpec
    assert module.HT1_004_LEAD_YIELD_STRESS_FITTED_PA == pytest.approx(0.0000)
    assert module.HT1_004_TAIL_YIELD_STRESS_FITTED_PA == pytest.approx(0.0000)
    assert module.HT1_004_LEAD_YIELD_STRESS_REPORT_GIVEN_PA == MISSING
    assert any("优化参数化" in note for note in validation.notes)


# --------------------------------------------------------------------------
# 5. 数据源在场性（防"引用一个不存在的文件"）
# --------------------------------------------------------------------------
_PROJECT_ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize("well_dir", [
    "hu101_呼101", "hu102_呼102", "hu103_呼103", "hu1_呼探1",
    "ht1_002_呼探1-002", "ht1_001_呼探1-001", "ht1_003_呼1-003", "ht1_004_呼1-004",
])
def test_rheometer_csv_is_present_and_has_the_promised_cement_rows(well_dir):
    """`rheometer_readings.csv` 在场且含承诺的水泥行（该目录被 gitignore，仅本地可读）。

    文件缺失时**跳过**（本目录不入库，CI/克隆环境无此文件）；文件在场时要求
    含水泥行——防"引用一个本地不存在的数据源"。
    """
    csv_path = _PROJECT_ROOT / "参考文档" / "现场资料提取" / well_dir / "rheometer_readings.csv"
    if not csv_path.exists():
        pytest.skip(f"现场提取包不在本地：{csv_path}")
    text = csv_path.read_text(encoding="utf-8-sig")
    assert "领浆" in text or "尾浆" in text, f"{well_dir} CSV 里没有水泥行"