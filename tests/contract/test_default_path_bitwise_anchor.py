# -*- coding: utf-8 -*-
"""默认路径端到端逐位锚（R10，由 Task 1 评审核实后新建）。

口径：hu101 生产 loader + T1 生产 1D 三开关 + 默认环空开关（不含 CORRECTED_KW），
微网格 nz=30/ny=12（只求逐位稳定，不求数值精度）。
total_t=14000s：brief 原定 200s 时入口尚为钻井液、五场恒为全零、五标量恒 0
（零判别力），故上调至覆盖整个顶替序列（泵注终点 13743s），使锚含活前缘。
锁：cement/lead/tail/spacer/wall 五场 sha256 + η_E/η_N/窜槽/混浆/失稳五个标量。
用途：任何声称"零数值影响"的改动都必须让本测试保持全绿。

缺锚语义（R17）：**锚文件缺失 ⇒ 直接 FAIL**，不走 `pytest.skip`。
skip 会让 CI 变绿而零验证（`pyproject.toml` 无 `addopts`、全仓无 `conftest.py` 兜底），
而本锚是后续所有"零数值影响"结论的唯一证据载体。缺锚时的正确动作是
`python tests/contract/_regenerate_default_path_anchor.py` 重锚并在提交中一并留痕。
键集/指纹/缺锚语义的元守卫见 `tests/contract/test_anchor_integrity.py`。

字段口径说明（2026-09-25 实测）：
  `AnnulusSimulationResult` 上不存在 `*_final` 字段，末帧位于非快照的场属性：
  cement→`cement_field`、lead→`lead_field`、tail→`tail_field`、
  spacer→`spacer_field`、wall→`wall_field`。解包时 `cement = cement_snapshots[-1]`
  等快照序列只承载历史帧，`*_field` 即求解结束时的末帧（见 annulus_d2dga.py
  `AnnulusSimulationResult` 构造处），故锚取 `*_field`。
"""
import hashlib
import json
from pathlib import Path

import numpy as np

from cemdisp.data.loaders import hu101_loader
from cemdisp.models2d import AnnulusD2DGASolver
from cemdisp.models2d.boundary_bridge import build_coupled_annulus_inlet_provider
from cemdisp.transport1d import CasingFlowSolver

ANCHOR = Path(__file__).with_name("_default_path_anchor_hu101.json")

# 五场在 `AnnulusSimulationResult` 上的末帧属性名（实测，无 `*_final` 字段）
_FIELD_ATTRS = {
    "cement": "cement_field",
    "lead": "lead_field",
    "tail": "tail_field",
    "spacer": "spacer_field",
    "wall": "wall_field",
}

# summary["最终结果"] 中五个标量的字段名（锚键 -> 中文列名）
_SCALAR_FIELDS = {
    "eta_E": "全井段最终有效顶替效率",
    "eta_N": "窄四分位效率",
    "channeling": "最终窜槽指数",
    "mixing": "最终混浆指数",
    "instability": "最终失稳指数",
}

# 锚 JSON 里的来源指纹键（不参与数值比对，见 test_anchor_integrity.py）
FINGERPRINT_KEYS = ("_note", "_generated_from", "_env")


def expected_keys() -> set[str]:
    """锚必须覆盖的键集（由字段表推导，禁止与断言行各写一份）。"""
    return {f"sha_{n}" for n in _FIELD_ATTRS} | set(_SCALAR_FIELDS)


def _run_default_case() -> dict:
    """跑一次默认路径端到端算例，返回五场 sha256 + 五个标量的指纹字典。"""
    well, fluids, schedule, _ = hu101_loader.load_hu101_tailpipe()
    cr = CasingFlowSolver(enable_gravity=True, mixing_contact_time=True,
                          plug_face_zero_mixing=True, has_plug=True).run(
        well, fluids, schedule)
    inlet = build_coupled_annulus_inlet_provider(
        cr, CasingFlowSolver(enable_gravity=True), fluids, split_cement_phases=True)
    solver = AnnulusD2DGASolver(total_t=14000.0, nz=30, ny=12)
    res = solver.run(well, fluids, inlet, schedule=schedule)
    out = {}
    for name, attr in _FIELD_ATTRS.items():
        arr = np.asarray(getattr(res, attr), dtype=float)
        out[f"sha_{name}"] = hashlib.sha256(arr.tobytes()).hexdigest()
    fr = res.summary["最终结果"]
    for key, field in _SCALAR_FIELDS.items():
        out[key] = float(fr[field])
    return out


def test_default_path_matches_bitwise_anchor():
    """默认路径数值必须与锚逐位一致；缺锚 ⇒ FAIL（R17，禁止 skip）。"""
    assert ANCHOR.is_file(), (
        f"逐位锚文件缺失：{ANCHOR}（R17：缺锚必须硬失败，严禁用 skip 掩盖；"
        "重锚请运行 tests/contract/_regenerate_default_path_anchor.py）")
    want = json.loads(ANCHOR.read_text(encoding="utf-8"))
    assert set(want) >= expected_keys(), (
        f"锚键集不齐：缺 {sorted(expected_keys() - set(want))}（R18）")
    got = _run_default_case()
    mismatched = {k for k in expected_keys() if got.get(k) != want.get(k)}
    assert not mismatched, (
        f"默认路径数值位移：{sorted(mismatched)}（若为有意改动，须走重锚："
        "与改动同批次、并在 results/内部自洽加固_2026-09-25/位移台账.csv 留迹）")
