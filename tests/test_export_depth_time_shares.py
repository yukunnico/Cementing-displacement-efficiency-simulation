# -*- coding: utf-8 -*-
"""全井深度×时间×流体占比导出（环空2D插值 + 套管1D sharp重建）单元测试。

合成数据单测，不依赖真实 NPZ / 现场资料：
- V(t) 分段积分与停泵段、碰压/序列终点冻结
- sharp 界面定位解析解（管内容积坐标查表，半开区间口径）
- 井身结构表深度域判定与 2D 格点逐位插值
- 泵序流体名 → 四通道映射
- utf-8-sig 写读 roundtrip、防覆盖守卫、到达时刻校验器
- 全管线集成：合成 NPZ + 表格 CSV + 鞋口时序 CSV → 6 文件导出
"""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from cemdisp.data.fluid_spec import FluidRole, FluidSpec, RheologyModel
from cemdisp.data.pumping_schedule import PumpingSchedule, PumpingScheduleStep

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SCRIPT_PATH = _REPO_ROOT / "scripts" / "entrypoints" / "export_depth_time_concentration.py"


@pytest.fixture(scope="module")
def mod():
    """按文件路径加载被测脚本模块（scripts/ 非 package）。"""
    spec = importlib.util.spec_from_file_location("_export_depth_time_concentration_test", _SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ── 合成数据工具 ──────────────────────────────────────────────────

def _fluids():
    return (
        FluidSpec("井浆", FluidRole.MUD, 1000.0, RheologyModel.NEWTONIAN, plastic_viscosity_pa_s=0.01),
        FluidSpec("平衡液", FluidRole.SPACER, 1200.0, RheologyModel.NEWTONIAN, plastic_viscosity_pa_s=0.02),
        FluidSpec("尾浆", FluidRole.TAIL, 1900.0, RheologyModel.NEWTONIAN, plastic_viscosity_pa_s=0.05),
        FluidSpec("压塞液", FluidRole.OTHER, 1200.0, RheologyModel.NEWTONIAN, plastic_viscosity_pa_s=0.01),
        FluidSpec("替钻井液", FluidRole.DISPLACEMENT, 1100.0, RheologyModel.NEWTONIAN, plastic_viscosity_pa_s=0.015),
    )


def _schedule_4step():
    """四步合成泵序：平衡液50 → 尾浆50 → 压塞液10 → 替钻井液120（总230 m³，排量全 1.0 m³/min）。"""
    return PumpingSchedule(
        steps=(
            PumpingScheduleStep("注入平衡液", "平衡液", 50.0, 1.0),
            PumpingScheduleStep("注入尾浆", "尾浆", 50.0, 1.0),
            PumpingScheduleStep("注入压塞液", "压塞液", 10.0, 1.0),
            PumpingScheduleStep("替钻井液", "替钻井液", 120.0, 1.0),
        )
    )


def _prepared_steps():
    """泵序 → 截断序列（复用 1D 求解器静态方法）。"""
    from cemdisp.transport1d.casing_flow import CasingFlowSolver

    steps_full = CasingFlowSolver._build_scheduled_steps(_schedule_4step())
    steps, _ = CasingFlowSolver._displacement_sequence_cutoff(steps_full)
    return steps


# ── V(t) 分段积分 / 停泵 / 冻结 ───────────────────────────────────

def test_volume_of_time_with_shutdown_and_freeze(mod):
    """停泵段（显式时间、体积0）V 不增；超过冻结体积后 V 恒定。"""
    schedule = PumpingSchedule(
        steps=(
            PumpingScheduleStep("注A", "平衡液", 60.0, 1.0,
                                start_time_s=0.0, end_time_s=3600.0),              # 0–60min
            PumpingScheduleStep("停泵", "平衡液", 0.0, 0.0,
                                start_time_s=3600.0, end_time_s=7200.0),           # 停泵60min
            PumpingScheduleStep("注B", "替钻井液", 30.0, 1.0,
                                start_time_s=7200.0, end_time_s=9000.0),           # 120–150min
        )
    )
    from cemdisp.transport1d.casing_flow import CasingFlowSolver

    steps_full = CasingFlowSolver._build_scheduled_steps(schedule)
    steps, _ = CasingFlowSolver._displacement_sequence_cutoff(steps_full)
    v_fn = mod.make_volume_of_time_fn(steps, v_frozen=90.0)
    assert v_fn(0.0) == pytest.approx(0.0)
    assert v_fn(1800.0) == pytest.approx(30.0)     # 注A中段
    assert v_fn(5400.0) == pytest.approx(60.0)     # 停泵段中点：V 不增
    assert v_fn(7200.0) == pytest.approx(60.0)
    assert v_fn(8400.0) == pytest.approx(80.0)     # 恢复泵注20min
    assert v_fn(1.0e9) == pytest.approx(90.0)      # 泵注结束后 = 总量
    v_fn_frozen = mod.make_volume_of_time_fn(steps, v_frozen=75.0)
    assert v_fn_frozen(8040.0) == pytest.approx(74.0)   # 冻结前正常（raw=74）
    assert v_fn_frozen(9000.0) == pytest.approx(75.0)   # 碰压冻结：V 钳在 v_frozen
    assert v_fn_frozen(1.0e9) == pytest.approx(75.0)


def test_volume_freeze_at_bump_matches_tail_plus_pipe_volume(mod):
    """冻结体积 = S_尾浆 + Vp_shoe 时，冻结后鞋口体积坐标 u_shoe = S_尾浆（胶塞面停在鞋口）。"""
    steps = _prepared_steps()
    vp_shoe = 100.0
    s_tail = 100.0  # 平衡液+尾浆累计
    v_fn = mod.make_volume_of_time_fn(steps, v_frozen=s_tail + vp_shoe)
    t_bump = mod.invert_volume_to_time(steps, s_tail + vp_shoe)
    assert t_bump == pytest.approx(12000.0)  # V=200 落在替浆步 [110,230)
    assert v_fn(t_bump) == pytest.approx(200.0)
    assert v_fn(t_bump + 3600.0) == pytest.approx(200.0)  # 碰压冻结


# ── sharp 界面定位解析解 ──────────────────────────────────────────

def test_sharp_interface_lookup_analytic(mod):
    """u = V(t) − v(z) 定位：解析解逐点核对（半开区间 [S_{j-1}, S_j)，u<0 → 初始钻井液）。"""
    steps = _prepared_steps()
    initial = "井浆"
    v = 60.0  # V(t)=60：平衡液已注完(50)，尾浆注入中段
    # 平衡液/尾浆界面在 v = V − S_1 = 10；井浆/平衡液界面在 v = V = 60
    assert mod.fluid_name_at_volume_coordinate(v - 5.0, steps, initial) == "尾浆"       # u=55
    assert mod.fluid_name_at_volume_coordinate(v - 10.0, steps, initial) == "尾浆"      # u=50 边界半开
    assert mod.fluid_name_at_volume_coordinate(v - 30.0, steps, initial) == "平衡液"    # u=30
    assert mod.fluid_name_at_volume_coordinate(v - 50.0, steps, initial) == "平衡液"    # u=10
    assert mod.fluid_name_at_volume_coordinate(v - 60.0, steps, initial) == "平衡液"    # u=0 边界半开
    assert mod.fluid_name_at_volume_coordinate(v - 60.001, steps, initial) == "井浆"    # u<0


def test_sharp_state_at_bump_freeze(mod):
    """碰压冻结后：胶塞面（尾浆/压塞液界面）停在鞋口，压塞液/替浆保持在胶塞上方。"""
    steps = _prepared_steps()
    vp_shoe = 100.0
    v_fn = mod.make_volume_of_time_fn(steps, v_frozen=100.0 + vp_shoe)  # V_bump = 200
    t_bump = mod.invert_volume_to_time(steps, 100.0 + vp_shoe)
    v_after = v_fn(t_bump + 600.0)
    assert v_after == pytest.approx(200.0)
    # u = V − v：鞋口处 u = 100（胶塞面），管内自上而下：替浆 → 压塞液 → 胶塞面
    assert mod.fluid_name_at_volume_coordinate(v_after - 5.0, steps, "井浆") == "替钻井液"   # u=195
    assert mod.fluid_name_at_volume_coordinate(v_after - 95.0, steps, "井浆") == "压塞液"    # u=105
    assert mod.fluid_name_at_volume_coordinate(v_after - 99.9999, steps, "井浆") == "压塞液"
    assert mod.fluid_name_at_volume_coordinate(v_after - 100.0, steps, "井浆") == "压塞液"   # u=100 半开
    # 替浆永不越过鞋口：鞋口处 u 最大 = V_frozen − Vp_shoe = S_tail
    assert mod.fluid_name_at_volume_coordinate(v_after - vp_shoe, steps, "井浆") != "替钻井液"


def test_no_bump_well_tail_stalls_in_casing(mod):
    """未碰压井（V_bump > 泵序总量）：V 冻结在序列终点，尾浆尾段滞留管内。"""
    steps = _prepared_steps()
    vp_shoe = 150.0  # V_bump = 250 > 230 → 无碰压
    v_fn = mod.make_volume_of_time_fn(steps, v_frozen=steps[-1].cumulative_volume_end_m3)
    assert v_fn(1.0e12) == pytest.approx(230.0)
    # 鞋口处 u_max = 230 − 150 = 80 ∈ 尾浆步 [50, 100) → 尾浆滞留
    assert mod.fluid_name_at_volume_coordinate(80.0, steps, "井浆") == "尾浆"


def test_arrival_times_inversion(mod):
    """前缘到达鞋口时刻 = V 反解（与 CasingFlowSolver._front_arrival_time 同口径）。"""
    from cemdisp.transport1d.casing_flow import CasingFlowSolver

    steps = _prepared_steps()
    vp_shoe = 100.0
    # 排量全 1.0 m³/min：平衡液前缘 target = 0+100 → 尾浆步内 50 m³ → 3000+3000 = 6000s
    assert CasingFlowSolver._front_arrival_time(steps[0], steps, vp_shoe) == pytest.approx(6000.0)
    # 尾浆前缘 target = 50+100 = 150 → 替浆步内 40 m³ → 6600+2400 = 9000s
    assert CasingFlowSolver._front_arrival_time(steps[1], steps, vp_shoe) == pytest.approx(9000.0)
    # 压塞液前缘 target = 100+100 = 200 → 替浆步内 90 m³ → 6600+5400 = 12000s
    assert CasingFlowSolver._front_arrival_time(steps[2], steps, vp_shoe) == pytest.approx(12000.0)
    # 超过泵序总量 → None
    assert CasingFlowSolver._front_arrival_time(steps[3], steps, 250.0) is None


# ── 管容剖面 ──────────────────────────────────────────────────────

def test_pipe_volume_profile_scaled_to_target(mod):
    """Vp(z) = pipe_id_profile 梯形积分 × k，鞋口处 Vp(shoe) = pipe_volume_target。"""
    depths = [0.0, 50.0, 100.0]
    ids = [100.0, 100.0, 100.0]  # 恒定 100mm 内径
    target = 0.6  # 剖面梯形积分 = π×0.05²×100 ≈ 0.7854 → k ≈ 0.7639
    d_out, vp_out, k = mod.build_pipe_volume_profile(depths, ids, 100.0, target)
    assert d_out[0] == pytest.approx(0.0) and d_out[-1] == pytest.approx(100.0)
    assert vp_out[0] == pytest.approx(0.0)
    assert vp_out[-1] == pytest.approx(target)
    assert k == pytest.approx(target / (np.pi * 0.05 ** 2 * 100.0))
    # 段内线性
    assert float(np.interp(50.0, d_out, vp_out)) == pytest.approx(target / 2.0)


def test_pipe_volume_profile_two_diameter(mod):
    """双内径剖面：段内体积按各自内径梯形积分（变径用 0.001m 过渡点，与真实 loader 同构）。"""
    depths = [0.0, 100.0, 100.001, 200.0]
    ids = [200.0, 200.0, 100.0, 100.0]  # 0–100m: 200mm; 100.001–200m: 100mm
    seg1 = np.pi * 0.1 ** 2 * 100.0
    taper = 0.001 * (np.pi * 0.1 ** 2 + np.pi * 0.05 ** 2) / 2
    seg2 = np.pi * 0.05 ** 2 * 99.999
    target = seg1 + taper + seg2  # 与剖面积分一致 → k=1
    d_out, vp_out, k = mod.build_pipe_volume_profile(depths, ids, 200.0, target)
    assert k == pytest.approx(1.0, rel=1e-9)
    assert float(np.interp(100.0, d_out, vp_out)) == pytest.approx(seg1, rel=1e-9)
    assert vp_out[-1] == pytest.approx(target, rel=1e-9)


# ── 深度域判定与 2D 插值 ──────────────────────────────────────────

def _synthetic_means(n_t=2, nz=4):
    """合成周向平均场：lead/tail/spacer 各为深度的线性函数。"""
    md_desc = np.array([130.0, 120.0, 110.0, 100.0])  # 降序（井底在前）
    fields = {}
    for i, ch in enumerate(("lead", "tail", "spacer")):
        base = np.linspace(0.1 * (i + 1), 0.5 * (i + 1), nz)
        fields[ch] = np.tile(base, (n_t, 1))
    return md_desc, fields


def test_interp_annulus_exact_at_nodes(mod):
    """表格深度 = 模型网格节点 → 逐位一致；节点间深度 → 线性插值解析解。

    合成约定与真实 NPZ 一致：md_desc 降序（index 0=最深），fields 逐位对应 md_desc，
    即 fields[:, 0] 是最深点的值、fields[:, -1] 是最浅点的值。
    （旧版此测试误把"浅深度配 fields[:,0]"的镜像期望当基准，2026-09-09 修正。）
    """
    md_desc, fields = _synthetic_means()
    table = np.array([100.0, 110.0, 120.0, 130.0])
    in_mask, in_depths, in_vals = mod.interp_annulus_profiles(fields, md_desc, table)
    assert in_mask.tolist() == [True] * 4
    for ch in ("lead", "tail", "spacer"):
        np.testing.assert_array_equal(in_vals[ch], fields[ch][:, ::-1])  # 节点逐位
    # 中点深度 105：位于 md_asc [100,110] 之间 → (fields[:, -1] + fields[:, -2]) / 2
    table2 = np.array([105.0])
    _, _, vals2 = mod.interp_annulus_profiles(fields, md_desc, table2)
    for ch in ("lead", "tail", "spacer"):
        expected = 0.5 * (fields[ch][:, -1] + fields[ch][:, -2])
        np.testing.assert_allclose(vals2[ch][:, 0], expected)


def test_interp_annulus_mirror_regression(mod):
    """镜像回归（2026-09-09 实锤缺陷）：md 降序喂入、场取非对称单调剖面，
    浅深度值必须出现在浅表格深度（而非镜像到深端）。
    缺陷实现（md 翻转而场未翻转）在本断言下必然失败。"""
    # 场沿"升序 md"递增：浅(100)=0.1 → 深(130)=0.9（非对称，镜像必错）
    md_desc = np.array([130.0, 120.0, 110.0, 100.0])
    field = np.tile(np.array([0.9, 0.7, 0.3, 0.1]), (2, 1))  # (n_t=2, nz=4)，与 md_desc 逐位对应
    table = np.array([100.0, 115.0, 130.0])
    _, in_depths, vals = mod.interp_annulus_profiles({"tail": field}, md_desc, table)
    np.testing.assert_allclose(vals["tail"][0], [0.1, 0.5, 0.9], atol=1e-12)


def test_domain_mask_tolerance(mod):
    """域边界 ±1e-6 容差；域外深度不外推。"""
    md_desc, fields = _synthetic_means()
    table = np.array([99.9, 99.9999995, 100.0, 130.0000005, 130.001])
    in_mask, in_depths, _ = mod.interp_annulus_profiles(fields, md_desc, table)
    assert in_mask.tolist() == [False, True, True, True, False]
    assert in_depths.tolist() == [99.9999995, 100.0, 130.0000005]


# ── 四通道映射 ────────────────────────────────────────────────────

def test_channel_mapping(mod):
    cases = {
        "平衡液": "spacer", "先导浆": "spacer", "先导液": "spacer",
        "隔离液1": "spacer", "隔离液2": "spacer", "冲洗液（FLUSHER）": "spacer",
        "领浆": "lead", "中间浆": "lead",
        "尾浆": "tail", "尾管水泥浆": "tail",
        "井浆": "mud", "钻井液": "mud", "替钻井液": "mud",
        "压塞液": "mud", "保护液": "mud", "基液": "mud", "泥浆": "mud",
    }
    for name, ch in cases.items():
        assert mod.map_fluid_to_channel(name) == ch, name
    with pytest.raises(ValueError):
        mod.map_fluid_to_channel("未知流体X")


# ── 防覆盖 / 校验器 / utf-8-sig ───────────────────────────────────

def test_check_no_overwrite(mod, tmp_path):
    existing = tmp_path / "x.csv"
    existing.write_text("", encoding="utf-8")
    with pytest.raises(RuntimeError):
        mod.check_no_overwrite([tmp_path / "new.csv", existing])
    mod.check_no_overwrite([tmp_path / "new.csv"])  # 不存在 → 通过


def test_validate_shoe_arrivals(mod):
    mine = [("平衡液", 100.0), ("尾浆", 200.0), ("替钻井液", None)]
    rows = [(0.0, "井浆"), (100.001, "平衡液"), (200.5, "尾浆"), (300.0, "替钻井液")]
    assert mod.validate_shoe_arrivals(mine, rows, tol_s=1.0) == []
    rows_bad = [(0.0, "井浆"), (102.0, "平衡液"), (200.0, "尾浆")]
    mismatches = mod.validate_shoe_arrivals(mine, rows_bad, tol_s=1.0)
    assert len(mismatches) == 1
    assert mismatches[0]["fluid"] == "平衡液"


def test_utf8sig_roundtrip(mod, tmp_path):
    rows = [("平衡液", 1.0), ("尾浆", 1.0), ("压塞液", 1.0)]
    path = tmp_path / "rt.csv"
    mod.write_csv_utf8sig(
        path,
        ("time_s", "depth_m", "fluid", "share", "source"),
        ((1.0, 30.0, name, share, "1D_casing") for name, share in rows),
    )
    assert path.read_bytes().startswith(b"\xef\xbb\xbf")  # BOM
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        data = list(csv.DictReader(fh))
    assert [r["fluid"] for r in data] == ["平衡液", "尾浆", "压塞液"]


# ── 全管线集成（合成 NPZ + 表格 CSV + 鞋口 CSV）────────────────────

def _write_synthetic_inputs(tmp_path):
    from cemdisp.transport1d.casing_flow import CasingFlowSolver

    npz_path = tmp_path / "synthetic_2D场数据.npz"
    table_csv = tmp_path / "synthetic_井身结构.csv"
    shoe_csv = tmp_path / "synthetic_鞋口出流时序.csv"
    out_dir = tmp_path / "out"

    md_desc = np.array([130.0, 120.0, 110.0, 100.0])
    times = np.array([5.0, 4000.0, 14000.0])
    rng = np.random.default_rng(7)
    snaps = {ch: rng.random((3, 5, 4)) * 0.3 for ch in ("lead", "tail", "spacer")}
    snaps["wall"] = np.zeros_like(snaps["lead"])
    np.savez(
        npz_path,
        md=md_desc, y=np.linspace(0, 0.3, 5), snapshot_times_s=times,
        lead_snapshots=snaps["lead"], tail_snapshots=snaps["tail"],
        spacer_snapshots=snaps["spacer"], wall_snapshots=snaps["wall"],
        cement_snapshots=snaps["lead"] + snaps["tail"],
        lead_final=snaps["lead"][-1], tail_final=snaps["tail"][-1],
        spacer_final=snaps["spacer"][-1], cement_final=snaps["lead"][-1] + snaps["tail"][-1],
        wall_final=snaps["wall"][-1],
    )

    # 井身结构表：20/50 域外（1D 套管段），100/105/110/120/130 域内（2D）
    depths = [20.0, 50.0, 100.0, 105.0, 110.0, 120.0, 130.0]
    with table_csv.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["VarName1", "depth_well_logging_m_", "other"])
        for i, d in enumerate(depths, start=1):
            writer.writerow([i, d, 0.0])

    # 鞋口出流时序 CSV：按同口径体积反解的到达时刻（排量全 1.0 m³/min）
    steps_full = CasingFlowSolver._build_scheduled_steps(_schedule_4step())
    steps, _ = CasingFlowSolver._displacement_sequence_cutoff(steps_full)
    vp_shoe = 0.7853981633974483  # 100mm 内径 × 100m
    shoe_rows = [(0.0, "井浆")]
    for step in steps:
        t = CasingFlowSolver._front_arrival_time(step, steps, vp_shoe)
        if t is not None:
            shoe_rows.append((round(t, 3), step.step.fluid_name))
    with shoe_csv.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["地面累计时间_s", "鞋口出流流体"])
        for t, name in shoe_rows:
            writer.writerow([f"{t:.3f}", name])
    return npz_path, table_csv, shoe_csv, out_dir, steps, vp_shoe


def test_full_pipeline_synthetic(mod, tmp_path):
    npz_path, table_csv, shoe_csv, out_dir, steps, vp_shoe = _write_synthetic_inputs(tmp_path)
    summary = mod.export_fullwell_depth_time_shares(
        npz_path, table_csv, out_dir, shoe_csv,
        well_label="合成井", schedule=_schedule_4step(), fluids=_fluids(),
        pipe_profile_depths=[0.0, 100.0],
        pipe_profile_ids_mm=[100.0, 100.0],
        shoe_md_m=100.0,
        pipe_volume_target=vp_shoe,
    )

    # 6 个产出文件
    long_csv = out_dir / "合成井_1D2D耦合模型_全井深度时间占比_长格式.csv"
    doc_md = out_dir / "合成井_1D2D耦合模型_全井深度时间占比_说明.md"
    wide = {ch: out_dir / f"合成井_1D2D耦合模型_宽格式_{ch}.csv" for ch in ("lead", "tail", "spacer", "mud")}
    for p in (long_csv, doc_md, *wide.values()):
        assert p.exists(), p
    assert long_csv.read_bytes().startswith(b"\xef\xbb\xbf")

    # 长格式行数 = n_t × (n_in×4 + n_out) = 3 × (5×4 + 2)
    with long_csv.open("r", encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 3 * (5 * 4 + 2)
    assert set(r["source"] for r in rows) == {"2D_annulus", "1D_casing"}
    ann_fluids = set(r["fluid"] for r in rows if r["source"] == "2D_annulus")
    assert ann_fluids == {"lead", "tail", "spacer", "mud"}

    # t=5s：V ≈ 0.083 < Vp(20)=0.157 → 全部套管深度为初始井浆
    casing_t0 = [r for r in rows if r["source"] == "1D_casing" and float(r["time_s"]) == 5.0]
    assert len(casing_t0) == 2  # 深度 20/50 各 1 行（sharp 单相）
    assert all(r["fluid"] == "井浆" and float(r["share"]) == 1.0 for r in casing_t0)

    # t=4000s：V ≈ 66.67 → 深度20/50（Vp≈0.157/0.393）u≈66.5/66.3 ∈ 尾浆步 [50,100)
    casing_t1 = {r["depth_m"]: r["fluid"] for r in rows
                 if r["source"] == "1D_casing" and float(r["time_s"]) == 4000.0}
    assert set(casing_t1.values()) == {"尾浆"}

    # 冻结末态（t=14000s）：raw V=233 冻结在 V_bump = 100 + Vp_shoe → 套管段全为压塞液（胶塞面上方）
    casing_t2 = {r["depth_m"]: r["fluid"] for r in rows
                 if r["source"] == "1D_casing" and float(r["time_s"]) == 14000.0}
    assert set(casing_t2.values()) == {"压塞液"}

    # 宽格式：行 = 全部表格深度，列 = 快照时刻，无空值
    for ch, p in wide.items():
        with p.open("r", encoding="utf-8-sig", newline="") as fh:
            wrows = list(csv.DictReader(fh))
        assert len(wrows) == 7
        assert list(wrows[0].keys()) == ["深度_m", "5.000_s", "4000.000_s", "14000.000_s"]
        assert all(all(v != "" for v in r.values()) for r in wrows)
    # 深度50、t=4000 → 尾浆通道 = 1.0（套管段 sharp 映射）
    with wide["tail"].open("r", encoding="utf-8-sig", newline="") as fh:
        tail_row50 = [r for r in csv.DictReader(fh) if r["深度_m"] == "50.0"][0]
    assert float(tail_row50["4000.000_s"]) == 1.0

    # 校验摘要回传
    assert summary["arrival_mismatches"] == []
    assert summary["n_in_domain"] == 5 and summary["n_out_domain"] == 2
    assert summary["t_bump_s"] == pytest.approx(6000.0 + vp_shoe * 60.0)
